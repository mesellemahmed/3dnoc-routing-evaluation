"""ERA structural signature of a typed NoC trace.

Implements the Entropic Structural Observability framework of
A. Mesellem and M. Mana, "A Trace-Based Structural Observability Framework
for Network-on-Chip Routing Evaluation", Algorithms 19(9):765, 2026.
Equation numbers in comments refer to that paper.

Default policy (paper, Table 2): warm-up excluded; router activity from
inject/forward/receive; directed-link activity from forward only; buffer-based
critical resources (Eq. 66-69); plug-in entropy with Miller-Madow sensitivity;
kappa_ST effect size with permutation p-value (B = 200); centre-of-mass
geometric drift. In addition, the TSV critical-resource mapping of Eq. 71 is
reported for 3D topologies.
"""
from __future__ import annotations

import json
import math
from collections import deque
from dataclasses import dataclass
from pathlib import Path

import numpy as np

TRACE_DTYPE = np.dtype([
    ("t", "<i4"), ("p", "<i4"), ("r", "<i2"), ("next", "<i2"),
    ("xi", "u1"), ("c", "u1"), ("port", "<u2"),
])
INJECT, FORWARD, RECEIVE, DROP = 0, 1, 2, 3
N_PORTS = 7
LOCAL, UP, DOWN = 0, 5, 6


@dataclass
class Topology:
    kx: int
    ky: int
    kz: int
    n_windows: int
    window: int
    links: np.ndarray       # (L, 2) directed links
    vertical: np.ndarray    # (L,) bool

    @property
    def n(self) -> int:
        return self.kx * self.ky * self.kz

    def coords(self) -> np.ndarray:
        r = np.arange(self.n)
        return np.stack([r % self.kx, (r // self.kx) % self.ky, r // (self.kx * self.ky)], axis=1)


def load_topology(path: Path) -> Topology:
    lines = path.read_text().split("\n")
    kx, ky, kz, nw, win = map(int, lines[0].split())
    rows = np.array([list(map(int, ln.split())) for ln in lines[1:] if ln.strip()], dtype=np.int64)
    return Topology(kx, ky, kz, nw, win, rows[:, :2], rows[:, 2].astype(bool))


# ---------------------------------------------------------------------------
# Information-theoretic primitives (Eqs. 31-37, 40-41, 44-45, 49)
# ---------------------------------------------------------------------------

def entropy_bits(counts: np.ndarray) -> float:
    total = counts.sum()
    if total <= 0:
        return float("nan")
    p = counts[counts > 0] / total
    return float(-(p * np.log2(p)).sum())


def entropy_report(counts: np.ndarray, m_full: int, n_events: float | None = None) -> dict:
    counts = np.asarray(counts, dtype=float)
    h = entropy_bits(counts)
    m_act = int((counts > 0).sum())
    out = {
        "H": h,
        "h_act": h / math.log2(m_act) if m_act > 1 else float("nan"),
        "h_full": h / math.log2(m_full) if m_full > 1 else float("nan"),
        "S_eff": 2.0 ** h if not math.isnan(h) else float("nan"),
        "m_act": m_act,
        "m_full": m_full,
    }
    out["C"] = 1.0 - out["h_act"]
    if n_events is not None and n_events > 0:
        out["MM_delta"] = (m_act - 1) / (2.0 * n_events * math.log(2.0))
    return out


def mutual_information(joint: np.ndarray) -> tuple[float, float, float, float]:
    h_r = entropy_bits(joint.sum(axis=1))
    h_w = entropy_bits(joint.sum(axis=0))
    h_rw = entropy_bits(joint.ravel())
    mi = h_r + h_w - h_rw
    den = min(h_r, h_w)
    kappa = mi / den if den > 0 else float("nan")
    return mi, kappa, h_r, h_w


def permutation_pvalue(joint: np.ndarray, mi_obs: float, n_perm: int, rng: np.random.Generator) -> float:
    """Eq. 47. Shuffling window labels across events with router labels fixed
    is equivalent to drawing a contingency table with the observed router and
    window margins: each router row is a multivariate hypergeometric draw from
    the remaining window totals. The draw is decomposed into sequential
    univariate hypergeometric draws over windows, vectorised over the B
    permutations. Since both margins are fixed, only the joint entropy varies
    across permutations."""
    rows = joint.sum(axis=1).astype(np.int64)
    cols = joint.sum(axis=0).astype(np.int64)
    n_rows, n_cols = joint.shape
    remaining = np.tile(cols, (n_perm, 1))
    h_joint = np.zeros(n_perm)
    total = float(rows.sum())
    for r in np.argsort(-rows):
        k = int(rows[r])
        if k == 0:
            continue
        need = np.full(n_perm, k, dtype=np.int64)
        rest = remaining.sum(axis=1)
        row = np.zeros((n_perm, n_cols), dtype=np.int64)
        for j in range(n_cols - 1):
            good = remaining[:, j]
            rest = rest - good
            x = np.zeros(n_perm, dtype=np.int64)
            live = need > 0
            if live.any():
                x[live] = rng.hypergeometric(good[live], rest[live], need[live])
            row[:, j] = x
            need = need - x
        row[:, n_cols - 1] = need
        remaining -= row
        with np.errstate(divide="ignore", invalid="ignore"):
            pr = row / total
            h_joint -= np.where(row > 0, pr * np.log2(pr), 0.0).sum(axis=1)
    h_r = entropy_bits(rows.astype(float))
    h_w = entropy_bits(cols.astype(float))
    mi_perm = h_r + h_w - h_joint
    exceed = int((mi_perm >= mi_obs - 1e-12).sum())
    return (1.0 + exceed) / (n_perm + 1.0)


def jsd_bits(p: np.ndarray, q: np.ndarray) -> float:
    m = 0.5 * (p + q)

    def kl(a: np.ndarray, b: np.ndarray) -> float:
        mask = a > 0
        return float((a[mask] * np.log2(a[mask] / b[mask])).sum())

    return 0.5 * kl(p, m) + 0.5 * kl(q, m)


# ---------------------------------------------------------------------------
# Classical imbalance baselines (Eqs. 56-60)
# ---------------------------------------------------------------------------

def imbalance(u: np.ndarray) -> dict:
    u = np.asarray(u, dtype=float)
    m, s = len(u), u.sum()
    mu = u.mean()
    srt = np.sort(u)
    idx = np.arange(1, m + 1)
    gini = float((2.0 * (idx * srt).sum() / (m * s)) - (m + 1.0) / m) if s > 0 else float("nan")
    ratio = u / mu if mu > 0 else np.zeros_like(u)
    with np.errstate(divide="ignore", invalid="ignore"):
        theil_terms = np.where(ratio > 0, ratio * np.log(ratio), 0.0)
    return {
        "jain": float(s * s / (m * (u * u).sum())) if s > 0 else float("nan"),
        "cv": float(u.std() / mu) if mu > 0 else float("nan"),
        "hotspot": float(u.max() / s) if s > 0 else float("nan"),
        "gini": gini,
        "theil": float(theil_terms.mean()),
    }


# ---------------------------------------------------------------------------
# Spatial baselines (Eqs. 61-64)
# ---------------------------------------------------------------------------

def adjacency(topo: Topology) -> np.ndarray:
    a = np.zeros((topo.n, topo.n))
    a[topo.links[:, 0], topo.links[:, 1]] = 1.0
    return np.maximum(a, a.T)


_DIAMETER_CACHE: dict = {}


def graph_diameter(adj: np.ndarray) -> int:
    key = hash(adj.tobytes())
    if key in _DIAMETER_CACHE:
        return _DIAMETER_CACHE[key]
    _DIAMETER_CACHE[key] = _graph_diameter(adj)
    return _DIAMETER_CACHE[key]


def _graph_diameter(adj: np.ndarray) -> int:
    n = adj.shape[0]
    nbrs = [np.nonzero(adj[i])[0] for i in range(n)]
    best = 0
    for s in range(n):
        dist = np.full(n, -1)
        dist[s] = 0
        q = deque([s])
        while q:
            u = q.popleft()
            for v in nbrs[u]:
                if dist[v] < 0:
                    dist[v] = dist[u] + 1
                    q.append(v)
        best = max(best, int(dist.max()))
    return best


def morans_i(x: np.ndarray, adj: np.ndarray) -> float:
    dev = x - x.mean()
    den = (dev * dev).sum()
    if den <= 0:
        return float("nan")
    return float(len(x) / adj.sum() * (dev @ adj @ dev) / den)


def smoothness(x: np.ndarray, lap: np.ndarray) -> float:
    den = x @ x
    return float(x @ lap @ x / den) if den > 0 else float("nan")


# ---------------------------------------------------------------------------
# Signature
# ---------------------------------------------------------------------------

def era_signature(prefix: Path, warmup: int, cycles: int, n_perm: int, seed: int) -> dict:
    topo = load_topology(prefix.with_suffix(".topo"))
    trace = np.fromfile(prefix.with_suffix(".trace"), dtype=TRACE_DTYPE)
    trace = trace[(trace["t"] >= warmup) & (trace["t"] < cycles)]
    n, w_count = topo.n, topo.n_windows
    win = ((trace["t"] - warmup) // topo.window).astype(np.int64)

    # Router activity, Eq. 12-13 (omega_R = 1 for inject, forward, receive).
    sel_r = np.isin(trace["xi"], (INJECT, FORWARD, RECEIVE))
    r_idx = trace["r"][sel_r].astype(np.int64)
    X = np.bincount(r_idx * w_count + win[sel_r], minlength=n * w_count).reshape(n, w_count).astype(float)

    # Directed-link activity, Eq. 15-16 and reconstruction policy, Eq. 19-22.
    fwd = trace[trace["xi"] == FORWARD]
    link_id = {(int(a), int(b)): i for i, (a, b) in enumerate(topo.links)}
    keys = fwd["r"].astype(np.int64) * n + fwd["next"].astype(np.int64)
    lut = np.full(n * n, -1, dtype=np.int64)
    for (a, b), i in link_id.items():
        lut[a * n + b] = i
    eid = np.where(fwd["next"] >= 0, lut[np.clip(keys, 0, n * n - 1)], -1)
    resolved = eid >= 0
    Y = np.bincount(eid[resolved], minlength=len(topo.links)).astype(float)
    rho_unres = float(1.0 - resolved.mean()) if len(fwd) else float("nan")

    # Temporal distribution, Eq. 27-28.
    A = X.sum(axis=0)

    # Buffer-based critical resources, Eq. 66-69.
    buf = np.fromfile(prefix.with_suffix(".buf"), dtype="<f4").reshape(w_count, n, N_PORTS)
    exists = np.zeros((n, N_PORTS), dtype=bool)
    exists[:, LOCAL] = True
    port_of = {(1, 0, 0): 2, (-1, 0, 0): 1, (0, 1, 0): 4, (0, -1, 0): 3, (0, 0, 1): 6, (0, 0, -1): 5}
    xyz = topo.coords()
    for a, b in topo.links:
        d = tuple(int(v) for v in (xyz[b] - xyz[a]))
        exists[b, port_of[d]] = True  # input port of b receiving from a
    B = buf[:, exists]
    Zbar = B.sum(axis=0)

    # TSV critical resources, Eq. 71: c = (x, y, TSV) column groups.
    vert = fwd[fwd["c"] == 1]
    col = (vert["r"].astype(np.int64) % (topo.kx * topo.ky))
    columns = np.unique(xyz[topo.links[topo.vertical, 0], 0] + topo.kx * xyz[topo.links[topo.vertical, 0], 1])
    tsv_counts = np.array([(col == c).sum() for c in columns], dtype=float)

    sig: dict = {}
    sig["n_events_router"] = int(sel_r.sum())
    sig["n_forward"] = int(len(fwd))
    sig["rho_unres"] = rho_unres
    sig["R"] = entropy_report(X.sum(axis=1), n, n_events=X.sum())
    sig["E"] = entropy_report(Y, len(topo.links), n_events=Y.sum())
    sig["T"] = entropy_report(A, w_count)
    sig["C_buf"] = entropy_report(Zbar, int(exists.sum()))
    sig["C_tsv"] = entropy_report(tsv_counts, len(columns), n_events=tsv_counts.sum()) if len(columns) else None
    sig["vertical_share"] = float(len(vert) / len(fwd)) if len(fwd) else float("nan")
    sig["tsv_jain"] = imbalance(tsv_counts)["jain"] if tsv_counts.sum() > 0 else float("nan")

    # Spatiotemporal mutual information, Eq. 42-47.
    mi, kappa, _, _ = mutual_information(X)
    rng = np.random.default_rng(seed)
    sig["MI"] = mi
    sig["kappa_ST"] = kappa
    sig["MI_p"] = permutation_pvalue(X, mi, n_perm, rng) if n_perm > 0 else float("nan")
    sig["MI_B"] = n_perm

    # Jensen-Shannon drift, Eq. 48-51; geometric drift, Eq. 53-55.
    adj = adjacency(topo)
    d_max = graph_diameter(adj)
    nonempty = [w for w in range(w_count) if X[:, w].sum() > 0]
    maps = {w: X[:, w] / X[:, w].sum() for w in nonempty}
    jsd = [jsd_bits(maps[w], maps[w + 1]) for w in range(w_count - 1) if w in maps and w + 1 in maps]
    sig["JSD_mean"] = float(np.mean(jsd)) if jsd else float("nan")
    sig["JSD_p95"] = float(np.quantile(jsd, 0.95)) if jsd else float("nan")
    sig["JSD_max"] = float(np.max(jsd)) if jsd else float("nan")
    centres = [maps[w] @ xyz for w in nonempty]
    dcom = [np.abs(centres[i + 1] - centres[i]).sum() / d_max for i in range(len(centres) - 1)]
    sig["geo_drift"] = float(np.mean(dcom)) if dcom else float("nan")
    sig["d_max"] = d_max

    # Classical imbalance on router totals (Eqs. 56-60) and spatial baselines.
    sig["imbalance_R"] = imbalance(X.sum(axis=1))
    lap = np.diag(adj.sum(axis=1)) - adj
    sig["moran"] = float(np.nanmean([morans_i(X[:, w], adj) for w in nonempty])) if nonempty else float("nan")
    sig["smoothness"] = float(np.nanmean([smoothness(X[:, w], lap) for w in nonempty])) if nonempty else float("nan")

    # Buffer pressure, Eq. 67-68.
    sig["buffer"] = {
        "mean": float(B.mean()), "p95": float(np.quantile(B, 0.95)),
        "p99": float(np.quantile(B, 0.99)), "max": float(B.max()),
    }
    return sig


PKT_COLS = 11  # src, dst, t_create, t_head_out, t_eject, hops, vhops, dropped, elevator, len, attacker


def load_packets(prefix: Path) -> np.ndarray:
    return np.fromfile(prefix.with_suffix(".pkt"), dtype=np.int32).reshape(-1, PKT_COLS)


def _jain(u: np.ndarray) -> float:
    u = np.asarray(u, dtype=float)
    return float(u.sum() ** 2 / (len(u) * (u * u).sum())) if len(u) and (u * u).sum() > 0 else float("nan")


def classical_metrics(prefix: Path, warmup: int, cycles: int, tail_guard: int, n_nodes: int) -> dict:
    """Packet-level metrics over packets created in [warmup, cycles - tail_guard).

    Packets ejected during the drain phase count as delivered; packets still in
    source queues when the run ends count as undelivered. Attack packets (DoS
    experiments) are excluded from the benign statistics and reported apart.
    """
    pk = load_packets(prefix)
    in_win = (pk[:, 2] >= warmup) & (pk[:, 2] < cycles - tail_guard)
    benign = pk[:, 10] == 0
    measured = pk[in_win & benign]
    dropped = measured[:, 7] == 1
    live = measured[~dropped]
    delivered = live[live[:, 4] >= 0]
    lat = (delivered[:, 4] - delivered[:, 2]).astype(float)
    net = (delivered[:, 4] - delivered[:, 3]).astype(float)
    span = cycles - warmup
    ej_win = (pk[:, 4] >= warmup) & (pk[:, 4] < cycles)
    acc_flits = float(pk[ej_win & benign, 9].sum())
    created_flits = float(pk[(pk[:, 2] >= warmup) & (pk[:, 2] < cycles) & benign, 9].sum())
    out = {
        "n_measured": int(len(measured)),
        "n_delivered": int(len(delivered)),
        "drop_ratio": float(dropped.mean()) if len(measured) else float("nan"),
        "undelivered": float(1.0 - len(delivered) / len(live)) if len(live) else float("nan"),
        "throughput_flits_per_cycle": acc_flits / span,
        "accepted_flits_per_node_cycle": acc_flits / span / n_nodes,
        "offered_flits_per_node_cycle": created_flits / span / n_nodes,
    }
    if len(lat):
        q = np.quantile(lat, [0.5, 0.95, 0.99, 0.999])
        out.update({
            "delay_mean": float(lat.mean()), "delay_p50": float(q[0]), "delay_p95": float(q[1]),
            "delay_p99": float(q[2]), "delay_p999": float(q[3]), "delay_max": float(lat.max()),
            "net_delay_mean": float(net.mean()), "queue_delay_mean": float((lat - net).mean()),
            "hops_mean": float(delivered[:, 5].mean()), "vhops_mean": float(delivered[:, 6].mean()),
        })
        src = delivered[:, 0]
        cnt = np.bincount(src, minlength=n_nodes)
        tot = np.bincount(src, weights=lat, minlength=n_nodes)
        ok = cnt >= 5
        per_src = tot[ok] / cnt[ok]
        out["src_lat_max_over_mean"] = float(per_src.max() / per_src.mean()) if ok.any() else float("nan")
        out["src_lat_jain"] = _jain(per_src) if ok.any() else float("nan")
        sources = np.unique(live[:, 0])
        thr = np.bincount(delivered[:, 0], weights=delivered[:, 9], minlength=n_nodes)[sources]
        out["src_thr_jain"] = _jain(thr)
    else:
        for k in ("delay_mean", "delay_p50", "delay_p95", "delay_p99", "delay_p999", "delay_max",
                  "net_delay_mean", "queue_delay_mean", "hops_mean", "vhops_mean",
                  "src_lat_max_over_mean", "src_lat_jain", "src_thr_jain"):
            out[k] = None
    att = pk[in_win & ~benign]
    if len(att):
        att_del = att[att[:, 4] >= 0]
        out["attack_packets"] = int(len(att))
        out["attack_delivered_frac"] = float(len(att_del) / len(att))
    return out


def activity(prefix: Path, n_nodes: int, kx: int, ky: int, keep_vector: bool) -> dict:
    act = np.fromfile(prefix.with_suffix(".act"), dtype=np.int64).reshape(n_nodes, 4)
    plane = kx * ky
    layers = n_nodes // plane
    out = {
        "xbar": int(act[:, 0].sum()), "planar": int(act[:, 1].sum()), "vertical": int(act[:, 2].sum()),
        "bufw": int(act[:, 3].sum()),
        "xbar_layer": [int(act[z * plane:(z + 1) * plane, 0].sum()) for z in range(layers)],
        "xbar_max": int(act[:, 0].max()),
    }
    if keep_vector:
        out["xbar_vec"] = act[:, 0].astype(int).tolist()
        out["vert_vec"] = act[:, 2].astype(int).tolist()
    return out


def analyse(prefix: Path, warmup: int = 3000, cycles: int = 18000, n_perm: int = 200,
            seed: int = 0, tail_guard: int = 1000, traced: bool = True, keep_vector: bool = True) -> dict:
    topo_line = prefix.with_suffix(".topo").read_text().split("\n", 1)[0].split()
    kx, ky, kz = int(topo_line[0]), int(topo_line[1]), int(topo_line[2])
    n = kx * ky * kz
    res = {"classical": classical_metrics(prefix, warmup, cycles, tail_guard, n),
           "activity": activity(prefix, n, kx, ky, keep_vector)}
    if traced:
        res["era"] = era_signature(prefix, warmup, cycles, n_perm, seed)
    return res


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("prefix", type=Path)
    ap.add_argument("--perm", type=int, default=200)
    args = ap.parse_args()
    print(json.dumps(analyse(args.prefix, n_perm=args.perm, keep_vector=False), indent=2, default=float))
