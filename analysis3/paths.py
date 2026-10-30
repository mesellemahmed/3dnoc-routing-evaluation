"""Exact path diversity of the evaluated routing relations on the 8x8x4 mesh.

For every ordered pair (s, d), s != d, we count the distinct router paths the
routing relation admits and the fraction alpha of all minimal paths it admits.
"""
from __future__ import annotations

import itertools
import math
from functools import lru_cache

import numpy as np

KX, KY, KZ = 8, 8, 4
N = KX * KY * KZ


def xyz(r):
    return r % KX, (r // KX) % KY, r // (KX * KY)


def rid(x, y, z):
    return x + KX * (y + KY * z)


def multinomial(*ks):
    n = sum(ks)
    out = math.factorial(n)
    for k in ks:
        out //= math.factorial(k)
    return out


def n_minimal(s, d):
    (sx, sy, sz), (dx, dy, dz) = xyz(s), xyz(d)
    return multinomial(abs(dx - sx), abs(dy - sy), abs(dz - sz))


def nf_paths(s, d):
    (sx, sy, sz), (dx, dy, dz) = xyz(s), xyz(d)
    neg = [max(0, sx - dx), max(0, sy - dy), max(0, sz - dz)]
    pos = [max(0, dx - sx), max(0, dy - sy), max(0, dz - sz)]
    return multinomial(*neg) * multinomial(*pos)


def ham_label(r):
    x, y, z = xyz(r)
    idx = y * KX + x if y % 2 == 0 else y * KX + (KX - 1 - x)
    if z % 2 == 1:
        idx = KX * KY - 1 - idx
    return z * KX * KY + idx


LAB = [ham_label(r) for r in range(N)]


def neighbours(r):
    x, y, z = xyz(r)
    for dx, dy, dz in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)):
        nx, ny, nz = x + dx, y + dy, z + dz
        if 0 <= nx < KX and 0 <= ny < KY and 0 <= nz < KZ:
            yield rid(nx, ny, nz)


def manhattan(a, b):
    (ax, ay, az), (bx, by, bz) = xyz(a), xyz(b)
    return abs(ax - bx) + abs(ay - by) + abs(az - bz)


def ham_counts(d):
    """Number of label-monotone paths and of minimal label-monotone paths from
    every node to d (dynamic programming in label order)."""
    ld = LAB[d]
    order = sorted(range(N), key=lambda u: abs(LAB[u] - ld))
    allp = np.zeros(N, dtype=object)
    minp = np.zeros(N, dtype=object)
    for u in order:
        if u == d:
            allp[u], minp[u] = 1, 1
            continue
        lu = LAB[u]
        a = m = 0
        for n in neighbours(u):
            ln = LAB[n]
            if (ld > lu and lu < ln <= ld) or (ld < lu and ld <= ln < lu):
                a += allp[n]
                if manhattan(n, d) < manhattan(u, d):
                    m += minp[n]
        allp[u], minp[u] = a, m
    return allp, minp


def dor_path(s, d, order):
    path = [s]
    cur = list(xyz(s))
    tgt = xyz(d)
    for dim in order:
        while cur[dim] != tgt[dim]:
            cur[dim] += 1 if tgt[dim] > cur[dim] else -1
            path.append(rid(*cur))
    return tuple(path)


def compute() -> dict:
    pairs = [(s, d) for s in range(N) for d in range(N) if s != d]
    res = {}
    nmin = np.array([n_minimal(s, d) for s, d in pairs], dtype=float)
    one = np.ones(len(pairs))
    # deterministic
    res["xyz"] = res["zxy"] = {"paths": one, "alpha": 1.0 / nmin}
    # negative-first (and DyAD, which shares its relation)
    nf = np.array([nf_paths(s, d) for s, d in pairs], dtype=float)
    res["nf"] = res["dyad"] = {"paths": nf, "alpha": nf / nmin}
    # minimal fully adaptive relations
    for c in ("fa", "rca", "qr", "ta"):
        res[c] = {"paths": nmin, "alpha": np.ones(len(pairs))}
    # Hamiltonian
    ham_all, ham_min = {}, {}
    for d in range(N):
        a, m = ham_counts(d)
        for s in range(N):
            ham_all[(s, d)] = a[s]
            ham_min[(s, d)] = m[s]
    res["ham"] = {"paths": np.array([float(ham_all[p]) for p in pairs]),
                  "alpha": np.array([float(ham_min[p]) for p in pairs]) / nmin}
    # O1TURN-3D: XYZ or ZYX
    o1 = np.array([len({dor_path(s, d, (0, 1, 2)), dor_path(s, d, (2, 1, 0))}) for s, d in pairs], dtype=float)
    res["o1t"] = {"paths": o1, "alpha": o1 / nmin}
    # ROMM: two-phase XYZ through any node of the minimal box
    romm = []
    for s, d in pairs:
        (sx, sy, sz), (dx, dy, dz) = xyz(s), xyz(d)
        seen = set()
        for mx in range(min(sx, dx), max(sx, dx) + 1):
            for my in range(min(sy, dy), max(sy, dy) + 1):
                for mz in range(min(sz, dz), max(sz, dz) + 1):
                    m = rid(mx, my, mz)
                    seen.add(dor_path(s, m, (0, 1, 2)) + dor_path(m, d, (0, 1, 2))[1:])
        romm.append(len(seen))
    romm = np.array(romm, dtype=float)
    res["romm"] = {"paths": romm, "alpha": romm / nmin}
    # Valiant: any intermediate node (non-minimal); distinct paths counted
    val = []
    for s, d in pairs:
        seen = set()
        for m in range(N):
            seen.add(dor_path(s, m, (0, 1, 2)) + dor_path(m, d, (0, 1, 2))[1:])
        val.append(len(seen))
    val = np.array(val, dtype=float)
    minimal_val = []
    res["val"] = {"paths": val, "alpha": np.full(len(pairs), np.nan)}
    out = {}
    for c, r in res.items():
        out[c] = {"multi_frac": float((r["paths"] > 1).mean()),
                  "paths_median": float(np.median(r["paths"])),
                  "alpha_mean": float(np.nanmean(r["alpha"])) if not np.all(np.isnan(r["alpha"])) else float("nan")}
    return out


if __name__ == "__main__":
    import json
    print(json.dumps(compute(), indent=1))
