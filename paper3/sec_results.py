"""Results by dimension (D1-D12), paired statistics and verdict."""
from __future__ import annotations

import json
from collections import defaultdict

import numpy as np
import pandas as pd

from common import ALL, ESCAPE, FC, GEN, M, MIN_VCS, NAME, ROOT, TR_NAME, TRAFFIC, VPC, cls, fmt_p, rank_marks, write
from energy_thermal import THERMAL_VARIANTS, energy_grid, energy_per_flit, temperatures
from stats import bootstrap_ci, holm, mean_ci, paired_delta, wilcoxon_p

REF = {c: ("xyz" if c in FC else "ef") for c in ALL}
VERDICT: dict = defaultdict(dict)   # VERDICT[class][dim] = "pass" | "fail" | "partial" | "na" | "ref"
RAW: dict = defaultdict(dict)       # RAW[class][dim] = statistic behind the verdict (threshold sensitivity)
BASE_GROUPS = ("era", "era_vpc", "sweep", "sweep_vpc", "stress", "dos", "burst", "vpcfault", "vpc4", "scal",
               "window", "period")


def mark(v: str) -> str:
    return {"pass": r"\cmark", "fail": r"\xmark", "partial": r"\pmark", "na": "--", "ref": "ref",
            "dl": r"\xmark$^{\ddagger}$", "cond": r"\xmark$^{\S}$"}[v]


FAIL_LIKE = ("fail", "partial", "dl", "cond")
NUM_WORD = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine"}


def counted(n: int, noun: str) -> str:
    """'one verdict', 'three verdicts' (number in words up to nine)."""
    return f"{NUM_WORD.get(n, n)} {noun}{'' if n == 1 else 's'}"


def dimension_runs(df: pd.DataFrame) -> dict:
    """Run set behind each simulated dimension (registered criteria, Table protocol)."""
    era = df.group.isin(["era", "era_vpc"])
    return {
        "D3": df.group.isin(["sweep", "sweep_vpc", "sweepref", "sweepref_vpc"]),
        "D4": era,
        "D5": df.group == "vpcfault",
        "D6": (era & (df.traffic == "random")) | ((df.group == "stress") & (df.traffic == "random")),
        "D9": era & (df.traffic == "random") & (df.pir == 0.006),
        "D10": era & df.traffic.isin(["random", "hotspot"]) & (df.pir == 0.006),
        "D11": df.group == "dos",
        "D12": (df.group == "sens") & (df.load == "over"),
    }


DLK: dict = {}   # DLK[dim][class] = deadlocked runs of the class in the run set of the dimension
DL_CHANGED: list = []


def prepare_dl(df: pd.DataFrame) -> None:
    for d, mask in dimension_runs(df).items():
        DLK[d] = {c: int(k) for c, k in df[mask].groupby("routing").deadlocked.sum().items()}


def dlmark(c: str, d: str) -> None:
    """Deadlock rule: a class that deadlocks in any run of a dimension's run set fails that
    dimension (double dagger), because its statistics would be computed on surviving runs only."""
    if DLK.get(d, {}).get(c, 0) > 0 and VERDICT[c].get(d) not in ("ref", "na", "dl"):
        DL_CHANGED.append((c, d, VERDICT[c].get(d)))
        VERDICT[c][d] = "dl"


# ---------------------------------------------------------------------------
# D1 correctness
# ---------------------------------------------------------------------------
def d1(df: pd.DataFrame, cdg: dict) -> None:
    rows = []
    sens = df[df.group == "sens"]
    for c in ALL:
        d = cdg[c]
        proven = (not d["cyclic"]) or (d["escape"] and not d["escape_cyclic"] and d["escape_to_adaptive"] == 0
                                        and d["adaptive_states_without_escape"] == 0)
        sub = df[df.routing == c]
        atomic = sub[sub.atomic == 1]
        nonatomic = sub[sub.atomic == 0]
        k_at, n_at = int(atomic.deadlocked.sum()), len(atomic)
        k_na, n_na = int(nonatomic.deadlocked.sum()), len(nonatomic)
        rows.append((c, proven, k_at, n_at, k_na, n_na))
        M.set(f"dl.{c}.k", k_at)
        M.set(f"dl.{c}.n", n_at)
        M.set(f"dl.{c}.na.k", k_na)
        M.set(f"dl.{c}.na.n", n_na)
        # Registered criterion: proven free AND no run of the class, in any group, deadlocked.
        # "cond": proven free and deadlock-free under atomic reallocation, but deadlocked under
        # non-atomic reallocation (Proposition 1) -- a fail under the registered criterion.
        if proven and k_at == 0 and k_na == 0:
            VERDICT[c]["D1"] = "pass"
        elif proven and k_at == 0:
            VERDICT[c]["D1"] = "cond"
        else:
            VERDICT[c]["D1"] = "fail"
    cond = [c for c in ALL if VERDICT[c]["D1"] == "cond"]
    M.set("d1.cond.classes", ", ".join(NAME[c] for c in cond))
    M.set("d1.cond.n", len(cond))
    M.set("d1.pass.n", sum(1 for c in ALL if VERDICT[c]["D1"] == "pass"))
    M.set("d1.pass.classes", ", ".join(NAME[c] for c in ALL if VERDICT[c]["D1"] == "pass"))
    lines = [r"\begin{tabular}{@{}lccrr@{}}", r"\toprule",
             r"Class & CDG verdict & Min.\ VCs & \multicolumn{2}{c}{Deadlocked runs} \\",
             r"\cmidrule(lr){4-5}", r" & & & atomic VC realloc. & non-atomic \\", r"\midrule"]
    for c, proven, k_at, n_at, k_na, n_na in rows:
        v = r"proven free" if proven else r"\textbf{cyclic}"
        na = f"{k_na}/{n_na}" if n_na else "--"
        lines.append(f"{cls(c)} & {v} & {MIN_VCS[c]} & {k_at}/{n_at} & {na} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    write("tab_d1.tex", "\n".join(lines) + "\n")
    # Escape-based classes under non-atomic reallocation.
    esc = sens[(sens.routing.isin(ESCAPE)) & (sens.atomic == 0)]
    M.set("dl.esc.na.k", int(esc.deadlocked.sum()))
    M.set("dl.esc.na.n", len(esc))
    esc_over = esc[esc.load == "over"]
    M.set("dl.esc.na.over.k", int(esc_over.deadlocked.sum()))
    M.set("dl.esc.na.over.n", len(esc_over))
    esc_low = esc[esc.load == "low"]
    M.set("dl.esc.na.low.k", int(esc_low.deadlocked.sum()))
    stuck = esc[esc.deadlocked]["st.stuck_flits"]
    M.set("dl.esc.na.stuck.min", int(stuck.min()) if len(stuck) else 0)
    M.set("dl.esc.na.stuck.max", int(stuck.max()) if len(stuck) else 0)
    esc_at = df[(df.routing.isin(ESCAPE)) & ~((df.group == "sens") & (df.atomic == 0))]
    M.set("dl.esc.atomic.k", int(esc_at.deadlocked.sum()))
    M.set("dl.esc.atomic.n", len(esc_at))
    partial_all = df[df.deadlocked & (df["st.deadlock_cycle"] < 0)]
    M.set("dl.partial.frac", f"{100 * len(partial_all) / max(1, int(df.deadlocked.sum())):.0f}")
    M.set("dl.total", int(df.deadlocked.sum()))
    M.set("dl.partial.n", len(partial_all))


# ---------------------------------------------------------------------------
# D2 path diversity
# ---------------------------------------------------------------------------
def d2() -> None:
    p = json.loads((ROOT / "analysis3" / "paths_result.json").read_text())
    vpc = json.loads((ROOT / "analysis3" / "paths_vpc.json").read_text())
    p.update(vpc)
    lines = [r"\begin{tabular}{@{}lrrr@{}}", r"\toprule",
             r"Class & $\alpha$ (minimal paths admitted) & Pairs with $>1$ path & Median paths \\", r"\midrule"]
    for c in ALL:
        r = p[c]
        a = "--" if not np.isfinite(r["alpha_mean"]) else f"{r['alpha_mean']:.3f}"
        med = r["paths_median"]
        meds = f"{med:,.0f}".replace(",", "{,}") if med < 1e7 else f"${med:.1e}$"
        lines.append(f"{cls(c)} & {a} & {100 * r['multi_frac']:.1f}\\% & {meds} \\\\")
        M.set(f"d2.{c}.alpha", r["alpha_mean"], ".3f")
        M.set(f"d2.{c}.multi", 100 * r["multi_frac"], ".1f")
        VERDICT[c]["D2"] = "pass" if r["multi_frac"] >= 0.5 else "fail"
        RAW[c]["D2"] = r["multi_frac"]
    lines += [r"\bottomrule", r"\end{tabular}"]
    write("tab_d2.tex", "\n".join(lines) + "\n")
    _m = p['ham']['paths_median']
    _e = int(np.floor(np.log10(_m)))
    M.set("d2.ham.median", rf"${_m / 10 ** _e:.1f}{{\times}}10^{{{_e}}}$")


# ---------------------------------------------------------------------------
# D3 performance to saturation
# ---------------------------------------------------------------------------
SAT_GROUPS = {"sweep": ("sweep", "sweepref"), "sweep_vpc": ("sweep_vpc", "sweepref_vpc")}
SAT_BRACKET: dict = {}   # (class, pattern) -> (saturation throughput, offered load of next unstable point)


def saturation_table(df: pd.DataFrame, group: str) -> pd.DataFrame:
    """Registered definition: highest accepted throughput among sweep points (coarse grid and
    bisection refinement) whose mean undelivered fraction is <= 1 % and with no deadlocked seed."""
    sw = df[df.group.isin(SAT_GROUPS[group])]
    m = sw.groupby(["routing", "traffic", "pir"]).agg(acc=("c.accepted_flits_per_node_cycle", "mean"),
                                                      off=("c.offered_flits_per_node_cycle", "mean"),
                                                      und=("c.undelivered", "mean"),
                                                      dl=("deadlocked", "sum"),
                                                      lat=("c.delay_mean", "mean")).reset_index()
    m["stable"] = (m.und <= 0.01) & (m.dl == 0)
    for (c, tr), g in m.groupby(["routing", "traffic"]):
        g = g.sort_values("pir")
        st = g[g.stable]
        if st.empty:
            continue
        top = st.loc[st.acc.idxmax()]
        above = g[(g.pir > top.pir) & ~g.stable]
        SAT_BRACKET[(c, tr)] = (float(top.acc), float(above.off.min()) if len(above) else np.nan)
    ok = m[m.stable]
    return ok.groupby(["routing", "traffic"])["acc"].max().unstack()


def d3(df: pd.DataFrame) -> dict:
    sat_fc = saturation_table(df, "sweep")
    sat_vpc = saturation_table(df, "sweep_vpc")
    sat = pd.concat([sat_fc, sat_vpc])
    zl = pd.concat([df[(df.group == "sweep") & (df.pir == 0.002)],
                    df[(df.group == "sweep_vpc") & (df.pir == 0.002)]]).groupby(["routing", "traffic"])[
        "c.delay_mean"].mean().unstack()
    sens = df[(df.group == "sens") & (df.variant == "base") & (df.load == "over")]
    mst = sens.groupby(["routing", "traffic"])["c.accepted_flits_per_node_cycle"].mean().unstack()
    out = {}
    for fam, classes in (("fc", FC), ("vpc", VPC)):
        cols = ["random", "transpose1", "bitreversal"]
        vals_sat = {tr: {c: sat.loc[c, tr] if c in sat.index else np.nan for c in classes if c != "efr"} for tr in cols}
        vals_zl = {tr: {c: zl.loc[c, tr] if c in zl.index else np.nan for c in classes if c != "efr"} for tr in cols}
        marks_sat = {tr: rank_marks(vals_sat[tr], lower_better=False, nd=3) for tr in cols}
        marks_zl = {tr: rank_marks(vals_zl[tr], lower_better=True, nd=1) for tr in cols}
        lines = [r"\begin{tabular}{@{}l" + "r" * 8 + r"@{}}", r"\toprule",
                 r"Class & \multicolumn{3}{c}{Zero-load latency (cycles)} & \multicolumn{3}{c}{Saturation throughput} & \multicolumn{2}{c}{Overload throughput} \\",
                 r"\cmidrule(lr){2-4}\cmidrule(lr){5-7}\cmidrule(lr){8-9}",
                 r" & Unif. & Transp. & Bit-rev. & Unif. & Transp. & Bit-rev. & Unif. & Transp. \\", r"\midrule"]
        for c in classes:
            if c == "efr":
                continue
            mrow = []
            for tr in ("random", "transpose1"):
                v = mst.loc[c, tr] if c in mst.index else np.nan
                mrow.append("--" if not np.isfinite(v) else f"{v:.3f}")
            lines.append(f"{cls(c)} & " + " & ".join(marks_zl[tr][c] for tr in cols) + " & " +
                         " & ".join(marks_sat[tr][c] for tr in cols) + " & " + " & ".join(mrow) + r" \\")
            for tr in cols:
                M.set(f"d3.sat.{c}.{tr}", vals_sat[tr][c], ".3f")
                M.set(f"d3.zl.{c}.{tr}", vals_zl[tr][c], ".1f")
            for tr in ("random", "transpose1"):
                if c in mst.index:
                    M.set(f"d3.mst.{c}.{tr}", mst.loc[c, tr], ".3f")
        lines += [r"\bottomrule", r"\end{tabular}"]
        write(f"tab_d3_{fam}.tex", "\n".join(lines) + "\n")
        for c in classes:
            ref = REF[c]
            if c == ref:
                VERDICT[c]["D3"] = "ref"
                continue
            if c == "efr":
                VERDICT[c]["D3"] = "na"
                continue
            s = {tr: (sat.loc[c, tr] if c in sat.index else np.nan) for tr in cols}
            r = {tr: sat.loc[ref, tr] for tr in cols}
            ok = s["random"] >= 0.9 * r["random"] and (s["transpose1"] > r["transpose1"] or s["bitreversal"] > r["bitreversal"])
            VERDICT[c]["D3"] = "pass" if ok else "fail"
            RAW[c]["D3"] = (s["random"] / r["random"] if np.isfinite(s["random"]) else 0.0,
                            bool(s["transpose1"] > r["transpose1"] or s["bitreversal"] > r["bitreversal"]))
            dlmark(c, "D3")
            M.set(f"d3.unif.ratio.{c}", s["random"] / r["random"] if np.isfinite(s["random"]) else np.nan, ".3f")
    VERDICT["efr"]["D3"] = "na"
    # Ratios used in the text, with the range allowed by the remaining load bracket.
    def ratio_range(a: str, b: str, tr: str) -> tuple:
        sa, ha = SAT_BRACKET.get((a, tr), (np.nan, np.nan))
        sb, hb = SAT_BRACKET.get((b, tr), (np.nan, np.nan))
        return sa / sb, sa / hb, ha / sb

    for key, a, b, tr, fmt in (("d3.br.fa.over.xyz", "fa", "xyz", "bitreversal", ".1f"),
                               ("d3.br.zxy.over.xyz", "zxy", "xyz", "bitreversal", ".1f"),
                               ("d3.tr.fa.over.xyz", "fa", "xyz", "transpose1", ".2f")):
        pt, lo, hi = ratio_range(a, b, tr)
        M.set(key, pt, fmt)
        M.set(key + ".lo", lo, fmt)
        M.set(key + ".hi", hi, fmt)
    # The same ratio on the coarse grid alone, before the bisection refinement (quoted in the protocol deviations).
    co = df[df.group == "sweep"].groupby(["routing", "traffic", "pir"]).agg(
        acc=("c.accepted_flits_per_node_cycle", "mean"), und=("c.undelivered", "mean"),
        dl=("deadlocked", "sum")).reset_index()
    cs = co[(co.und <= 0.01) & (co.dl == 0)].groupby(["routing", "traffic"])["acc"].max()
    M.set("d3.br.fa.over.xyz.coarse", cs[("fa", "bitreversal")] / cs[("xyz", "bitreversal")], ".1f")
    widths = [(h - s) / s for (s, h) in SAT_BRACKET.values() if np.isfinite(h) and s > 0]
    M.set("d3.bracket.max", 100 * max(widths), ".1f")
    M.set("d3.bracket.median", 100 * float(np.median(widths)), ".1f")
    M.set("d3.thr.unif", 0.9 * sat.loc["xyz", "random"], ".3f")
    for c in FC + VPC:
        for tr in cols:
            if (c, tr) in SAT_BRACKET:
                M.set(f"d3.sat.{c}.{tr}.next", SAT_BRACKET[(c, tr)][1], ".3f")
    out["sat"] = sat
    out["zl"] = zl
    # Latency-load curves (pgfplots coordinates), stable points only.
    sw = df[df.group == "sweep"]
    m = sw.groupby(["routing", "traffic", "pir"]).agg(lat=("c.delay_mean", "mean"), und=("c.undelivered", "mean"),
                                                      dl=("deadlocked", "sum"),
                                                      off=("c.offered_flits_per_node_cycle", "mean")).reset_index()
    curves = []
    for tr in ("random", "transpose1", "bitreversal"):
        for c in ("xyz", "zxy", "nf", "ham", "o1t", "romm", "val", "fa", "qr"):
            pts = m[(m.routing == c) & (m.traffic == tr) & (m.und <= 0.01) & (m.dl == 0) & (m.lat <= 150)]
            coords = " ".join(f"({r.off:.3f},{r.lat:.1f})" for r in pts.itertuples())
            curves.append(f"\\expandafter\\def\\csname curve@{c}@{tr}\\endcsname{{{coords}}}")
    write("curves.tex", "\n".join(curves) + "\n")
    return out


# ---------------------------------------------------------------------------
# D4 traffic robustness (ERA grid) and bursty traffic
# ---------------------------------------------------------------------------
def d4(df: pd.DataFrame) -> None:
    for fam, grp, classes in (("fc", "era", FC), ("vpc", "era_vpc", [c for c in VPC if c != "efr"])):
        era = df[df.group == grp]
        # "Best class": lowest mean latency among the classes that never deadlocked in these runs
        # (a deadlocked class's mean over its surviving runs is biased towards light loads).
        kdl = era.groupby("routing").deadlocked.sum()
        safe = [c for c in classes if int(kdl.get(c, 0)) == 0]
        lat = era[~era.deadlocked].groupby(["routing", "traffic"])["c.delay_mean"].mean().unstack()
        best = lat.loc[safe].min()
        lines = [r"\begin{tabular}{@{}l" + "r" * 7 + r"c@{}}", r"\toprule",
                 "Class & " + " & ".join(TR_NAME[t] for t in TRAFFIC) + r" & Worst ratio & D4 \\", r"\midrule"]
        marks = {tr: rank_marks({c: lat.loc[c, tr] for c in safe}, lower_better=True, nd=1) for tr in TRAFFIC}
        for c in classes:
            cells = []
            for tr in TRAFFIC:
                v = lat.loc[c, tr]
                if c in safe:
                    cells.append(marks[tr][c])
                else:
                    k = int(era[(era.routing == c) & (era.traffic == tr)].deadlocked.sum())
                    cells.append(f"{v:.1f}" + (r"$^{\dagger}$" if k else ""))
            if c in safe:
                ratio = max(lat.loc[c, tr] / best[tr] for tr in TRAFFIC)
                worst_tr = max(TRAFFIC, key=lambda tr: lat.loc[c, tr] / best[tr])
                VERDICT[c]["D4"] = "pass" if ratio <= 2.0 else "fail"
                RAW[c]["D4"] = ratio
                rs = f"{ratio:.2f} ({TR_NAME[worst_tr]})"
                M.set(f"d4.{c}.ratio", ratio, ".2f")
                M.set(f"d4.{c}.worst", TR_NAME[worst_tr])
            else:
                VERDICT[c]["D4"] = "dl"
                rs = "--"
            lines.append(f"{cls(c)} & " + " & ".join(cells) + f" & {rs} & {mark(VERDICT[c]['D4'])} \\\\")
            for tr in TRAFFIC:
                M.set(f"d4.{c}.{tr}", lat.loc[c, tr], ".1f")
        lines += [r"\bottomrule", r"\end{tabular}"]
        write(f"tab_d4_{fam}.tex", "\n".join(lines) + "\n")
    VERDICT["efr"]["D4"] = "na"
    # Censoring diagnostics: hotspot runs leave packets undelivered from the lowest loads on.
    for fam, grp, classes in (("fc", "era", FC), ("vpc", "era_vpc", [c for c in VPC if c != "efr"])):
        e_ = df[df.group == grp]
        h_ = e_[e_.traffic == "hotspot"]
        for pir, u in h_.groupby("pir")["c.undelivered"].mean().items():
            M.set(f"d4.{fam}.hot.und.{pir:g}", 100 * u, ".0f")
        M.set(f"d4.{fam}.hot.sat.frac", 100 * float((h_["c.undelivered"] > 0.01).mean()), ".0f")
        M.set(f"d4.{fam}.sat.frac", 100 * float((e_["c.undelivered"] > 0.01).mean()), ".0f")
        kd = e_.groupby("routing").deadlocked.sum()
        sf = [c for c in classes if int(kd.get(c, 0)) == 0]
        lt = e_[~e_.deadlocked].groupby(["routing", "traffic"])["c.delay_mean"].mean().unstack()
        bst = lt.loc[sf].min()
        flips = []
        for c in sf:
            rn = max(lt.loc[c, t] / bst[t] for t in TRAFFIC if t != "hotspot")
            M.set(f"d4.{c}.ratio.nohot", rn, ".2f")
            if VERDICT[c].get("D4") == "fail" and rn <= 2.0:
                flips.append(NAME[c])
        M.set(f"d4.{fam}.nohot.flips", ", ".join(flips) if flips else "none")
        M.set(f"d4.{fam}.nohot.nflips", len(flips))
        M.set(f"d4.{fam}.nfail", sum(1 for c in classes if VERDICT[c].get("D4") in ("fail", "dl")))
    era = df[df.group == "era"]
    lat = era[~era.deadlocked].groupby(["routing", "traffic"])["c.delay_mean"].mean().unstack()
    M.set("d4.hot.fa.over.xyz", lat.loc["fa", "hotspot"] / lat.loc["xyz", "hotspot"], ".1f")
    M.set("d4.hot.nf.over.xyz", lat.loc["nf", "hotspot"] / lat.loc["xyz", "hotspot"], ".1f")
    M.set("d4.hot.o1t.over.xyz", lat.loc["o1t", "hotspot"] / lat.loc["xyz", "hotspot"], ".1f")
    M.set("d4.hot.zxy.over.xyz", lat.loc["zxy", "hotspot"] / lat.loc["xyz", "hotspot"], ".2f")

    M.set("d4.br.xyz.over.zxy", lat.loc["xyz", "bitreversal"] / lat.loc["zxy", "bitreversal"], ".1f")
    M.set("d4.br.xyz.over.fa", lat.loc["xyz", "bitreversal"] / lat.loc["fa", "bitreversal"], ".1f")
    M.set("d4.br.o1t.over.fa", lat.loc["o1t", "bitreversal"] / lat.loc["fa", "bitreversal"], ".1f")
    # Hotspot structural signature (ERA) relative to XYZ, all hotspot runs that did not deadlock.
    hot = era[(era.traffic == "hotspot") & ~era.deadlocked]
    sig = hot.groupby("routing").agg(hc=("e.C_buf.h_full", "mean"), mi=("e.MI", "mean"), jsd=("e.JSD_mean", "mean"),
                                     jain=("e.imbalance_R.jain", "mean"))
    many = [c for c in ("nf", "dyad", "romm", "fa", "rca", "qr", "ta") if c in sig.index]
    mir = {c: sig.loc[c, "mi"] / sig.loc["xyz", "mi"] for c in many}
    M.set("d7.mi.adapt.min", min(mir.values()), ".1f")
    M.set("d7.mi.adapt.max", max(mir.values()), ".1f")
    M.set("d7.mi.adapt.min.class", ", ".join(NAME[c] for c in many if round(mir[c], 1) == round(min(mir.values()), 1)))
    M.set("d7.mi.adapt.max.class", ", ".join(NAME[c] for c in many if round(mir[c], 1) == round(max(mir.values()), 1)))
    for c in FC:
        if c in sig.index:
            M.set(f"d7.hot.{c}.hc", sig.loc[c, "hc"], ".3f")
            M.set(f"d7.hot.{c}.mi.ratio", sig.loc[c, "mi"] / sig.loc["xyz", "mi"], ".1f")
            M.set(f"d7.hot.{c}.jsd.ratio", sig.loc[c, "jsd"] / sig.loc["xyz", "jsd"], ".1f")
    # Bursty (self-similar) traffic: p99 latency, PIR 0.009 and 0.012.
    bur = df[df.group == "burst"]
    lines = [r"\begin{tabular}{@{}lrrrr@{}}", r"\toprule",
             r"Class & \multicolumn{2}{c}{Mean latency} & \multicolumn{2}{c}{99th-percentile latency} \\",
             r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}", r" & PIR 0.009 & PIR 0.012 & PIR 0.009 & PIR 0.012 \\", r"\midrule"]
    tab = bur.groupby(["routing", "pir"]).agg(lat=("c.delay_mean", "mean"), p99=("c.delay_p99", "mean"),
                                              und=("c.undelivered", "mean")).unstack()
    nsat = 0
    for fam_classes in (FC, [c for c in VPC if c != "efr"]):
        mk = {}
        for col in (("lat", 0.009), ("lat", 0.012), ("p99", 0.009), ("p99", 0.012)):
            # saturated points (> 1 % undelivered) are shown but not ranked, as in the stride-4 table
            stable = {c: tab.loc[c, col] for c in fam_classes if tab.loc[c, ("und", col[1])] <= 0.01}
            mk[col] = rank_marks(stable, lower_better=True, nd=1) if len(stable) > 1 else {}
            for c in fam_classes:
                if c not in mk[col]:
                    mk[col][c] = f"{tab.loc[c, col]:.1f}" + ("$^{\\ast}$" if tab.loc[c, ("und", col[1])] > 0.01 else "")
        for c in fam_classes:
            nsat += sum(1 for p_ in (0.009, 0.012) if tab.loc[c, ("und", p_)] > 0.01)
            lines.append(f"{cls(c)} & " + " & ".join(mk[col][c] for col in mk) + r" \\")
            M.set(f"d4.burst.{c}.p99", tab.loc[c, ("p99", 0.012)], ".1f")
            M.set(f"d4.burst.{c}.lat", tab.loc[c, ("lat", 0.012)], ".1f")
            M.set(f"d4.burst.{c}.lat9", tab.loc[c, ("lat", 0.009)], ".1f")
            M.set(f"d4.burst.{c}.und9", 100 * tab.loc[c, ("und", 0.009)], ".0f")
            M.set(f"d4.burst.{c}.und12", 100 * tab.loc[c, ("und", 0.012)], ".0f")
        if fam_classes is FC:
            lines.append(r"\midrule")
    lines += [r"\bottomrule", r"\end{tabular}"]
    write("tab_burst.tex", "\n".join(lines) + "\n")
    M.set("d4.burst.nsat", nsat)


# ---------------------------------------------------------------------------
# D5 fault tolerance
# ---------------------------------------------------------------------------
POLICY = {"xyz": "xyz", "zxy": "zxy", "o1t": "o1t", "romm": "romm", "val": "val", "nf": "nf_greedy",
          "dyad": "nf_greedy", "ham": "ham_greedy", "fa": "min_greedy", "rca": "min_greedy", "qr": "min_greedy",
          "ta": "min_greedy"}


def d5(df: pd.DataFrame, fr: pd.DataFrame) -> None:
    pv = fr.pivot_table(index="policy", columns=["p_h", "p_v", "xi"], values="reach_mean")
    pany = fr.pivot_table(index="policy", columns=["p_h", "p_v", "xi"], values="p_any_unreachable")
    confs = [(0.02, 0.0, 0.0), (0.05, 0.0, 0.0), (0.10, 0.0, 0.0), (0.0, 0.10, 0.0), (0.0, 0.10, 4.0),
             (0.0, 0.20, 0.0), (0.0, 0.20, 4.0), (0.0, 0.30, 0.0), (0.0, 0.30, 4.0)]
    order = ["xyz", "zxy", "o1t", "romm", "val", "nf_greedy", "nf_oracle", "ham_greedy", "ham_oracle", "min_order",
             "min_greedy", "min_oracle", "connected"]
    pname = {"xyz": r"\cls{XYZ}", "zxy": r"\cls{ZXY}", "o1t": r"\cls{O1T}", "romm": r"\cls{ROMM}",
             "val": r"\cls{VAL}", "nf_greedy": r"\cls{NF}, \cls{DyAD} (local)", "nf_oracle": r"NF (oracle)",
             "ham_greedy": r"\cls{HAM} (local)", "ham_oracle": r"HAM (oracle)",
             "min_order": r"min.\ adaptive, X--Y--Z port order (local)",
             "min_greedy": r"min.\ adaptive, local (\cls{FA}, \cls{RCA}, \cls{QR}, \cls{TA})",
             "min_oracle": r"min.\ adaptive (oracle)", "connected": r"connectivity bound (any routing)"}
    sub = " & ".join((f"{ph:g}" if pvv == 0 else f"{pvv:g} ({xi:g})") for ph, pvv, xi in confs)
    lines = [r"\begin{tabular}{@{}l" + "r" * len(confs) + r"@{}}", r"\toprule",
             r"Policy & \multicolumn{3}{c}{Planar faults, $p_h$} & \multicolumn{6}{c}{Vertical faults, $p_v$ ($\xi$)} \\",
             r"\cmidrule(lr){2-4}\cmidrule(lr){5-10}", " & " + sub + r" \\", r"\midrule"]
    for pol in order:
        vals = {cf: pv.loc[pol, cf] for cf in confs}
        lines.append(pname[pol] + " & " + " & ".join(f"{vals[cf]:.3f}" for cf in confs) + r" \\")
        if pol in ("min_oracle",):
            lines.append(r"\midrule")
    lines += [r"\bottomrule", r"\end{tabular}"]
    write("tab_d5_fc.tex", "\n".join(lines) + "\n")
    for pol in order:
        for (ph, pvv, xi) in confs:
            M.set(f"d5.{pol}.{ph:g}.{pvv:g}.{xi:g}", pv.loc[pol, (ph, pvv, xi)], ".3f")
    M.set("d5.pany.ham.0.05", pany.loc["ham_greedy", (0.0, 0.05, 0.0)], ".3f")
    # Lower bound (rounded down to 0.001) of the connectivity bound over the configurations of the table.
    M.set("d5.conn.min", np.floor(1000 * min(pv.loc["connected", cf] for cf in confs)) / 1000, ".3f")
    M.set("d5.pany.min.0.05", pany.loc["min_greedy", (0.0, 0.05, 0.0)], ".3f")
    # verdicts (FC)
    for c, pol in POLICY.items():
        if c == "xyz":
            VERDICT[c]["D5"] = "ref"
            continue
        checks = [(0.05, 0.0, 0.0), (0.0, 0.2, 0.0), (0.0, 0.2, 4.0)]
        gains = [pv.loc[pol, cf] - pv.loc["xyz", cf] for cf in checks]
        RAW[c]["D5"] = ("fc", gains)
        if all(g >= 0.05 for g in gains):
            VERDICT[c]["D5"] = "pass"
        elif any(g > 1e-3 for g in gains):
            VERDICT[c]["D5"] = "partial"
        else:
            VERDICT[c]["D5"] = "fail"
    # correlation penalty
    for pol in ("min_greedy", "ham_greedy", "nf_greedy", "xyz"):
        for p in (0.2, 0.3):
            M.set(f"d5.corr.{pol}.{p:g}", pv.loc[pol, (0.0, p, 4.0)] - pv.loc[pol, (0.0, p, 0.0)], "+.3f")
    # VPC simulation under correlated TSV faults.
    vf = df[df.group == "vpcfault"]
    t = vf.groupby(["routing", "pv", "xi"]).agg(drop=("c.drop_ratio", "mean"), lat=("c.delay_mean", "mean"),
                                                cols=("e.C_tsv.S_eff", "mean")).reset_index()
    lines = [r"\begin{tabular}{@{}l" + "rrr" * 2 + r"c@{}}", r"\toprule",
             r"Class & \multicolumn{3}{c}{$p_v=0.2$, $\xi=0$} & \multicolumn{3}{c}{$p_v=0.2$, $\xi=4$} & D5 \\",
             r"\cmidrule(lr){2-4}\cmidrule(lr){5-7}",
             r" & loss & latency & $S_{\mathrm{eff}}^{\mathrm{TSV}}$ & loss & latency & $S_{\mathrm{eff}}^{\mathrm{TSV}}$ & \\",
             r"\midrule"]
    sdt = vf.groupby(["routing", "pv", "xi"])["c.drop_ratio"].std()
    for c in VPC:
        for pvv in (0.1, 0.2):
            for xi in (0.0, 4.0):
                M.set(f"d5.vpc.{c}.{pvv:g}.{xi:g}.sd", sdt.get((c, pvv, xi), np.nan), ".3f")
    for c in VPC:
        cells = []
        for xi in (0.0, 4.0):
            r = t[(t.routing == c) & (t.pv == 0.2) & (t.xi == xi)].iloc[0]
            cells += [f"{r['drop']:.3f}", f"{r['lat']:.1f}", f"{r['cols']:.1f}"]
            M.set(f"d5.vpc.{c}.0.2.{xi:g}.drop", r["drop"], ".3f")
            M.set(f"d5.vpc.{c}.0.2.{xi:g}.lat", r["lat"], ".1f")
            M.set(f"d5.vpc.{c}.0.2.{xi:g}.cols", r["cols"], ".1f")
        maxdrop = t[t.routing == c]["drop"].max()
        VERDICT[c]["D5"] = "pass" if maxdrop == 0 else "fail"
        RAW[c]["D5"] = ("vpc", maxdrop)
        dlmark(c, "D5")
        lines.append(f"{cls(c)} & " + " & ".join(cells) + f" & {mark(VERDICT[c]['D5'])} \\\\")
        M.set(f"d5.vpc.{c}.maxdrop", maxdrop, ".3f")
        M.set(f"d5.vpc.{c}.runmax", float(vf[vf.routing == c]["c.drop_ratio"].max()), ".2f")
    lines += [r"\bottomrule", r"\end{tabular}"]
    write("tab_d5_vpc.tex", "\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# D6 load balance and fairness
# ---------------------------------------------------------------------------
def d6(df: pd.DataFrame) -> None:
    lines = [r"\begin{tabular}{@{}lrrrrc@{}}", r"\toprule",
             r"Class & Router-load Jain & Source-throughput Jain & Max/mean source latency & Max VA wait & D6 \\",
             r" & (uniform, grid) & (stress, min.) & (stress, mean) & (cycles) & \\", r"\midrule"]
    era = pd.concat([df[df.group == "era"], df[df.group == "era_vpc"]])
    jr = era[(era.traffic == "random") & ~era.deadlocked].groupby("routing")["e.imbalance_R.jain"].mean()
    st_all = df[(df.group == "stress") & (df.traffic == "random")]
    ndl = st_all.groupby("routing")["deadlocked"].sum()
    st = st_all[~st_all.deadlocked]
    sj = st.groupby("routing")["c.src_thr_jain"].min()
    sl = st.groupby("routing")["c.src_lat_max_over_mean"].mean()
    va = st.groupby("routing")["st.max_va_wait"].max()
    for fam in (FC, VPC):
        for c in fam:
            if c == "efr":
                VERDICT[c]["D6"] = "na"
                continue
            ref = REF[c]
            cond1 = jr.get(c, np.nan) >= jr.get(ref, np.nan) - 1e-12
            cond2 = sj.get(c, np.nan) >= 0.95
            RAW[c]["D6"] = (float(jr.get(c, np.nan) - jr.get(ref, np.nan)), float(sj.get(c, np.nan)))
            if c == ref:
                VERDICT[c]["D6"] = "ref" if cond2 else "fail"
            else:
                VERDICT[c]["D6"] = "pass" if (cond1 and cond2) else "fail"
            dlmark(c, "D6")
            k = int(ndl.get(c, 0))
            dag = r"$^{\ddagger}$" if k else ""
            lines.append(f"{cls(c)} & {jr.get(c, np.nan):.3f} & {sj.get(c, np.nan):.3f}{dag} & "
                         f"{sl.get(c, np.nan):.2f}{dag} & {int(va.get(c, 0)):,}{dag} & "
                         f"{mark(VERDICT[c]['D6'])} \\\\".replace(",", "{,}"))
            M.set(f"d6.{c}.jain", jr.get(c, np.nan), ".3f")
            M.set(f"d6.{c}.srcjain", sj.get(c, np.nan), ".3f")
            M.set(f"d6.{c}.srclat", sl.get(c, np.nan), ".2f")
            M.set(f"d6.{c}.va", int(va.get(c, 0)))
            M.set(f"d6.{c}.stressdl", k)
        if fam is FC:
            lines.append(r"\midrule")
    lines += [r"\bottomrule", r"\end{tabular}"]
    write("tab_d6.tex", "\n".join(lines) + "\n")
    best_src = max((c for c in ALL if c != "efr"), key=lambda c: sj.get(c, -1))
    M.set("d6.best.src", NAME[best_src])
    M.set("d6.best.src.v", sj.get(best_src), ".3f")
    # Range of the per-source Jain index over individual stress runs (all loads).
    sjr = st.groupby("routing")["c.src_thr_jain"].max()
    M.set("d6.srcjain.maxrun", float(sjr.max()), ".3f")
    M.set("d6.srcjain.min.lo", float(sj.min()), ".3f")
    lo_c = min((c for c in ALL if c != "efr" and c in sj.index), key=lambda c: sj[c])
    M.set("d6.srcjain.min.lo.class", NAME[lo_c])
    M.set("d6.xyz.srcjain.min", sj.get("xyz"), ".3f")
    M.set("d6.va.over10k", sum(1 for c in ALL if c != "efr" and va.get(c, 0) > 10000))


# ---------------------------------------------------------------------------
# D7 structure (ERA): window sensitivity, permutation tests, VPC TSV usage
# ---------------------------------------------------------------------------
def d7(df: pd.DataFrame) -> None:
    traced = df[df["e.MI"].notna()]
    p = traced["e.MI_p"].dropna().values
    # Benjamini-Hochberg at q = 0.05 over all traced runs.
    order = np.sort(p)
    m = len(order)
    k = 0
    for i, v in enumerate(order, 1):
        if v <= 0.05 * i / m:
            k = i
    M.set("d7.perm.n", m)
    M.set("d7.perm.bh.frac", k / m, ".3f")
    M.set("d7.perm.bh.pct", 100 * k / m, ".1f")
    M.set("d7.perm.min.p", float(order.min()), ".4f")
    M.set("d7.rho.max", float(traced["e.rho_unres"].max()), ".3f")
    # Window sensitivity: same runs re-analysed with 250/500/1000/2000-cycle windows.
    win = df[df.group == "window"]
    base = df[(df.group == "era") & (df.pir == 0.012) & (df.routing.isin(["xyz", "fa"])) &
              (df.traffic.isin(["random", "hotspot"])) & (df.seed <= 3)]
    allw = pd.concat([win, base])
    g = allw[~allw.deadlocked].groupby(["routing", "traffic", "window"]).agg(hR=("e.R.h_act", "mean"),
                                                                             mi=("e.MI", "mean"),
                                                                             kap=("e.kappa_ST", "mean"),
                                                                             jsd=("e.JSD_mean", "mean")).reset_index()
    rng = g.groupby(["routing", "traffic"]).agg(hR_lo=("hR", "min"), hR_hi=("hR", "max"), mi_lo=("mi", "min"),
                                                mi_hi=("mi", "max"), jsd_lo=("jsd", "min"), jsd_hi=("jsd", "max"))
    M.set("d7.win.hR.range", float((rng.hR_hi - rng.hR_lo).max()), ".3f")
    M.set("d7.win.mi.ratio", float((rng.mi_hi / rng.mi_lo).max()), ".1f")
    M.set("d7.win.jsd.ratio", float((rng.jsd_hi / rng.jsd_lo).max()), ".1f")
    # VPC: effective number of TSV columns used (fault-free grid, uniform traffic).
    ev = df[(df.group == "era_vpc") & (df.traffic == "random")]
    cols = ev.groupby("routing")["e.C_tsv.S_eff"].mean()
    for c in cols.index:
        M.set(f"d7.vpc.{c}.cols", cols[c], ".1f")
    hot = df[(df.group == "era_vpc") & (df.traffic == "hotspot")].groupby("routing")["e.C_tsv.S_eff"].mean()
    for c in hot.index:
        M.set(f"d7.vpchot.{c}.cols", hot[c], ".1f")


# ---------------------------------------------------------------------------
# D8 hardware cost (Yosys / NanGate 45 nm)
# ---------------------------------------------------------------------------
def read_area(module: str) -> float:
    txt = (ROOT / "hw" / f"rpt_{module}.txt").read_text()
    import re
    return float(re.search(r"Chip area for module '\\" + module + r"': ([0-9.]+)", txt).group(1))


D8_UNITS = ["rc_xyz", "rc_zxy", "rc_nf", "rc_dyad", "rc_ham", "rc_o1t", "rc_twophase", "rc_fa", "rca_select",
            "rca_agg", "qr_select", "qr_update", "rc_ta", "rc_ef", "rc_fl", "adele_sel", "vc_fifo", "xbar7"]
UNIT_DESC = {"rc_xyz": "the XYZ route unit", "rc_zxy": "the ZXY route unit", "rc_nf": "the negative-first route unit",
             "rc_dyad": "the DyAD route unit", "rc_ham": "the Hamiltonian route unit", "rc_o1t": "the O1TURN route unit",
             "rc_twophase": "the two-phase route unit", "rc_fa": "the minimal adaptive route unit",
             "rca_select": "the RCA selection unit", "rca_agg": "the RCA aggregation unit",
             "qr_select": "the Q-routing selection unit", "qr_update": "the Q-update unit",
             "rc_ta": "the thermal-aware route unit", "rc_ef": "the Elevator-First route unit",
             "rc_fl": "the First-Last route unit", "adele_sel": "AdEle's selection unit", "vc_fifo": "the VC buffer",
             "xbar7": "the crossbar"}
DFF_AREA = 4.522          # um^2 per stored bit (DFF_X1)
Q_BITS = 256 * 6 * 16     # tabular Q-router: destinations x output ports x 16-bit entries


def d8_cost(a: dict) -> tuple:
    """Per-class route/selection logic, storage, header and side-band bits, overhead and D8 verdict for unit areas a."""
    comp = {
        "xyz": (7 * a["rc_xyz"], 0, 0, 0), "zxy": (7 * a["rc_zxy"], 0, 0, 0), "nf": (7 * a["rc_nf"], 0, 0, 0),
        "dyad": (7 * a["rc_dyad"], 0, 0, 0), "ham": (7 * a["rc_ham"], 0, 0, 0), "o1t": (7 * a["rc_o1t"], 0, 1, 0),
        "romm": (7 * a["rc_twophase"], 0, 9, 0), "val": (7 * a["rc_twophase"], 0, 9, 0),
        "fa": (7 * a["rc_fa"], 0, 0, 0), "rca": (7 * a["rca_select"] + a["rca_agg"], 6 * 7, 0, 7),
        "qr": (7 * a["qr_select"] + 6 * a["qr_update"], Q_BITS, 0, 16), "ta": (7 * a["rc_ta"], 6 * 8, 0, 8),
        "ef": (7 * a["rc_ef"], 4, 4, 0), "efr": (7 * a["rc_ef"], 4 + 16, 4, 0),
        "aesg": (7 * a["rc_ef"], 16 * 8, 4, None), "adele": (7 * a["rc_ef"] + a["adele_sel"], 82, 4, 0),
        "fl": (7 * a["rc_fl"], 12, 4, 0),
    }
    base = 14 * a["vc_fifo"] + a["xbar7"] + 7 * a["rc_xyz"]
    out = {}
    for c in ALL:
        logic, bits, hdr, side = comp[c]
        over = (logic - 7 * a["rc_xyz"] + bits * DFF_AREA) / base
        glob = side is None
        out[c] = (logic, bits, hdr, side, over, glob, "pass" if (over <= 0.05 and not glob) else "fail")
    return out, base


def d8() -> None:
    a = {m: read_area(m) for m in D8_UNITS}
    cost, base = d8_cost(a)
    M.set("hw.fifo", a["vc_fifo"], ".0f")
    M.set("hw.xbar", a["xbar7"], ".0f")
    M.set("hw.base", base, ",.0f")
    M.set("hw.qtable.bits", Q_BITS)
    M.set("hw.qtable.area", Q_BITS * DFF_AREA, ",.0f")
    M.set("hw.qtable.ratio", Q_BITS * DFF_AREA / base, ".1f")
    M.set("hw.extra.vc", 7 * a["vc_fifo"] / (7 * a["vc_fifo"] + a["xbar7"] + 7 * a["rc_xyz"]), ".2f")
    lines = [r"\begin{tabular}{@{}lrrrrrc@{}}", r"\toprule",
             r"Class & Route/selection logic & Storage & Header & Side-band & Overhead & D8 \\",
             r" & ($\mu$m$^2$, per router) & (bits) & (bits) & (bits/link) & (\% of router) & \\", r"\midrule"]
    for c in ALL:
        logic, bits, hdr, side, over, glob, v = cost[c]
        sb = "global" if glob else (str(side) if side else "--")
        VERDICT[c]["D8"] = v
        RAW[c]["D8"] = (over, glob)
        ov = f"{100 * over:.1f}".replace("-0.0", "0.0")
        lines.append(f"{cls(c)} & {logic:,.0f} & {bits:,} & {hdr} & {sb} & {ov} & "
                     f"{mark(VERDICT[c]['D8'])} \\\\".replace(",", "{,}"))
        M.set(f"hw.{c}.logic", logic, ",.0f")
        M.set(f"hw.{c}.over", ov)
    lines += [r"\bottomrule", r"\end{tabular}"]
    write("tab_d8.tex", "\n".join(lines) + "\n")
    for m_ in ("rc_xyz", "rc_fa", "qr_select", "adele_sel", "rc_ham", "rc_nf"):
        M.set(f"hw.mod.{m_}", a[m_], ",.0f")
    # Sensitivity to the synthesis context (hw/synth_context.py): each unit synthesised alone with its submodules,
    # and from the full file with the modules in reverse order.
    ctx = json.loads((ROOT / "hw" / "synth_context.json").read_text())
    for u in D8_UNITS:
        if abs(ctx[u]["committed"] - a[u]) > 1e-6:
            raise SystemExit(f"d8: hw/synth_context.json is stale for {u}; re-run hw/synth_context.py")
    dev = {(u, var): ctx[u][var] / a[u] - 1 for u in D8_UNITS for var in ("alone", "reversed")}
    (u_max, _), d_max = max(dev.items(), key=lambda kv: abs(kv[1]))
    if max(abs(dev[(u, "reversed")]) for u in D8_UNITS) > 0:
        raise SystemExit("d8: reversing the module order changes an area; update the text of Section 7.6")
    alt, _ = d8_cost({u: ctx[u]["alone"] for u in D8_UNITS})
    changed = [NAME[c] for c in ALL if alt[c][6] != cost[c][6]]
    if changed:
        raise SystemExit(f"d8: the synthesis context changes the D8 verdict of {changed}; update the text")
    M.set("hw.ctx.max", 100 * abs(d_max), ".0f")
    M.set("hw.ctx.unit", UNIT_DESC[u_max])
    M.set("hw.ctx.pp", 100 * max(abs(alt[c][4] - cost[c][4]) for c in ALL), ".1f")


# ---------------------------------------------------------------------------
# D9 energy, D10 thermal
# ---------------------------------------------------------------------------
def d9_d10(df: pd.DataFrame) -> None:
    grid = energy_grid()
    lines9 = [r"\begin{tabular}{@{}lrrrc@{}}", r"\toprule",
              r"Class & min & median & max & D9 \\", r"\midrule"]
    res = {}
    for fam, grp in (("fc", "era"), ("vpc", "era_vpc")):
        sub = df[(df.group == grp) & (df.traffic == "random") & (df.pir == 0.006)]
        per = defaultdict(lambda: defaultdict(list))
        for r in sub.to_dict("records"):
            for prm in grid:
                per[r["routing"]][prm].append(energy_per_flit(r, *prm))
        ref = "xyz" if fam == "fc" else "ef"
        for c in per:
            ratios = [np.mean(per[c][prm]) / np.mean(per[ref][prm]) for prm in grid]
            res[c] = (min(ratios), float(np.median(ratios)), max(ratios))
    for c in ALL:
        if c not in res:
            continue
        lo, med, hi = res[c]
        ref = REF[c]
        if c == ref:
            VERDICT[c]["D9"] = "ref"
        else:
            VERDICT[c]["D9"] = "pass" if hi <= 1.05 else "fail"
            RAW[c]["D9"] = hi
            dlmark(c, "D9")
        lines9.append(f"{cls(c)} & {lo:.3f} & {med:.3f} & {hi:.3f} & {mark(VERDICT[c]['D9'])} \\\\")
        M.set(f"d9.{c}.max", hi, ".3f")
        M.set(f"d9.{c}.med", med, ".3f")
        M.set(f"d9.{c}.pct", 100 * (med - 1), ".1f")
    VERDICT["efr"]["D9"] = "na"
    lines9 += [r"\bottomrule", r"\end{tabular}"]
    write("tab_d9.tex", "\n".join(lines9) + "\n")
    # Thermal.
    lines10 = [r"\begin{tabular}{@{}l" + "r" * 8 + r"c@{}}", r"\toprule",
               r"Class & \multicolumn{4}{c}{Uniform traffic} & \multicolumn{4}{c}{Hotspot traffic} & D10 \\",
               r"\cmidrule(lr){2-5}\cmidrule(lr){6-9}",
               r" & base & NoC-heavy & weak vert. & strong lat. & base & NoC-heavy & weak vert. & strong lat. & \\",
               r"\midrule"]
    peaks = {}
    for fam, grp in (("fc", "era"), ("vpc", "era_vpc")):
        for tr in ("random", "hotspot"):
            sub = df[(df.group == grp) & (df.traffic == tr) & (df.pir == 0.006) & ~df.deadlocked]
            for var in THERMAL_VARIANTS:
                for c, g in sub.groupby("routing"):
                    vals, sds = [], []
                    for r in g.to_dict("records"):
                        T = temperatures(r["a.xbar_vec"], r["cycles"] - r["warmup"], var)
                        vals.append(T.max())
                        sds.append(T.std())
                    peaks[(c, tr, var)] = (float(np.mean(vals)), float(np.mean(sds)))
    for c in ALL:
        if c == "efr":
            VERDICT[c]["D10"] = "na"
            continue
        ref = REF[c]
        cells, ok, worst = [], True, -np.inf
        for tr in ("random", "hotspot"):
            for var in THERMAL_VARIANTS:
                if (c, tr, var) not in peaks:
                    cells.append("--")
                    continue
                dlt = peaks[(c, tr, var)][0] - peaks[(ref, tr, var)][0]
                cells.append(f"{dlt:+.3f}")
                ok = ok and dlt <= 0.05
                worst = max(worst, dlt)
                M.set(f"d10.{c}.{tr}.{var}", dlt, "+.2f")
        if c == ref:
            VERDICT[c]["D10"] = "ref"
        elif any(x == "--" for x in cells):
            VERDICT[c]["D10"] = "dl"
        else:
            VERDICT[c]["D10"] = "pass" if ok else "fail"
            RAW[c]["D10"] = worst
            dlmark(c, "D10")
        M.set(f"d10.{c}.worst", worst, "+.2f")
        M.set(f"d10.{c}.worstabs", abs(worst), ".2f")
        bestv = min(float(v_) for v_ in cells if v_ != "--") if any(v_ != "--" for v_ in cells) else np.nan
        M.set(f"d10.{c}.best", bestv, "+.2f")
        M.set(f"d10.{c}.bestabs", abs(bestv), ".2f")
        lines10.append(f"{cls(c)} & " + " & ".join("0.000" if c_ in ("-0.000", "+0.000") else c_ for c_ in cells) +
                       f" & {mark(VERDICT[c]['D10'])} \\\\")
    # The caption of tab:d10 has no legend for missing (all-deadlocked) cells; fail loudly if one appears.
    if any("--" in l_ for l_ in lines10[4:]):
        raise SystemExit("d10: a class has no non-deadlocked run; add a legend for '--' to the caption of tab:d10")
    lines10 += [r"\bottomrule", r"\end{tabular}"]
    write("tab_d10.tex", "\n".join(lines10) + "\n")
    M.set("d10.xyz.abs", peaks[("xyz", "random", "base")][0], ".1f")
    M.set("d10.ta.sd.ratio", peaks[("ta", "random", "base")][1] / peaks[("xyz", "random", "base")][1], ".2f")
    M.set("d10.val.sd.ratio", peaks[("val", "random", "base")][1] / peaks[("xyz", "random", "base")][1], ".2f")


# ---------------------------------------------------------------------------
# D11 security (flooding DoS)
# ---------------------------------------------------------------------------
def d11(df: pd.DataFrame) -> None:
    dos = df[df.group == "dos"]
    ok_runs = dos[~dos.deadlocked]
    p = ok_runs.groupby(["routing", "attackers", "pir"]).agg(lat=("c.delay_mean", "mean")).reset_index()
    dl = dos.groupby(["routing", "attackers", "pir"]).agg(dl=("deadlocked", "sum"),
                                                         n=("deadlocked", "size")).reset_index()
    # The caption of tab:d11 states that no flooding run deadlocked; fail loudly if new data contradict it.
    if int(dos.deadlocked.sum()):
        raise SystemExit("d11: deadlocked flooding runs; restore a deadlock column in tab_d11 and its caption")
    lines = [r"\begin{tabular}{@{}lrrrrc@{}}", r"\toprule",
             r"Class & \multicolumn{2}{c}{Benign latency, no attack} & \multicolumn{2}{c}{Inflation factor under flooding} & D11 \\",
             r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}", r" & PIR 0.003 & PIR 0.006 & PIR 0.003 & PIR 0.006 & \\",
             r"\midrule"]
    infl = {}

    def get(frame, c, att, pir, col):
        r = frame[(frame.routing == c) & (frame.attackers == att) & (frame.pir == pir)]
        return float(r[col].iloc[0]) if len(r) else np.nan

    for c in ALL:
        if c == "efr":
            continue
        row = {}
        for pir in (0.003, 0.006):
            b = get(p, c, 0, pir, "lat")
            a = get(p, c, 8, pir, "lat")
            k = int(get(dl, c, 8, pir, "dl")) + int(get(dl, c, 0, pir, "dl"))
            n = int(get(dl, c, 8, pir, "n")) + int(get(dl, c, 0, pir, "n"))
            row[pir] = (b, a / b if np.isfinite(a) and np.isfinite(b) else np.nan, k, n)
        infl[c] = row
    for fam in (FC, [c for c in VPC if c != "efr"]):
        mk = {pir: rank_marks({c: infl[c][pir][1] for c in fam if infl[c][pir][2] == 0}, lower_better=True, nd=1)
              for pir in (0.003, 0.006)}
        for c in fam:
            ref = REF[c]
            dlr = infl[c][0.003][2] + infl[c][0.006][2]
            nr = infl[c][0.003][3] + infl[c][0.006][3]
            if c == ref:
                VERDICT[c]["D11"] = "ref"
            elif dlr:
                VERDICT[c]["D11"] = "dl"
            else:
                ok = all(infl[c][pir][1] <= infl[ref][pir][1] for pir in (0.003, 0.006))
                VERDICT[c]["D11"] = "pass" if ok else "fail"
                RAW[c]["D11"] = max(infl[c][pir][1] / infl[ref][pir][1] for pir in (0.003, 0.006))
                dlmark(c, "D11")
            cells = []
            for pir in (0.003, 0.006):
                v = infl[c][pir]
                if v[2]:
                    cells.append(("--" if not np.isfinite(v[1]) else f"{v[1]:.1f}") + r"$^{\ddagger}$")
                else:
                    cells.append(mk[pir][c])
            lines.append(f"{cls(c)} & {infl[c][0.003][0]:.1f} & {infl[c][0.006][0]:.1f} & {cells[0]} & "
                         f"{cells[1]} & {mark(VERDICT[c]['D11'])} \\\\")
            for pir in (0.003, 0.006):
                M.set(f"d11.{c}.{pir}", infl[c][pir][1], ".1f")
            M.set(f"d11.{c}.dl", f"{dlr}/{nr}")
        if fam is FC:
            lines.append(r"\midrule")
    VERDICT["efr"]["D11"] = "na"
    lines += [r"\bottomrule", r"\end{tabular}"]
    write("tab_d11.tex", "\n".join(lines) + "\n")
    und = dos[dos.attackers == 8].groupby(["routing", "pir"])["c.undelivered"].mean()
    alt_diff = []
    for c in ALL:
        if c == "efr":
            continue
        for pir in (0.003, 0.006):
            M.set(f"d11.und.{c}.{pir}", und.get((c, pir), np.nan), ".2f")
        ref = REF[c]
        if c != ref and VERDICT[c]["D11"] in ("pass", "fail"):
            alt = "pass" if all(und[(c, p_)] <= und[(ref, p_)] for p_ in (0.003, 0.006)) else "fail"
            if alt != VERDICT[c]["D11"]:
                alt_diff.append(NAME[c])
                # The text says that these classes "would pass" under the alternative metric.
                if alt != "pass":
                    raise SystemExit(f"d11: {NAME[c]} would fail under the undelivered fraction; update the text")
    if not alt_diff:
        raise SystemExit("d11: no verdict depends on the metric; update the text of RQ5 and of the threats")
    M.set("d11.alt.diff", ", ".join(alt_diff[:-1]) + (" and " if len(alt_diff) > 1 else "") + alt_diff[-1])
    M.set("d11.alt.ndiff", len(alt_diff))
    M.set("d11.alt.nverdicts", counted(len(alt_diff), "verdict"))
    M.set("d11.alt.except", ("that of " if len(alt_diff) == 1 else "those of ") + M.v["d11.alt.diff"])
    M.set("d11.alt.which", "this class" if len(alt_diff) == 1 else f"these {NUM_WORD.get(len(alt_diff), len(alt_diff))} classes")
    rel = {}
    for c in FC:
        if c == "xyz" or infl[c][0.003][2] + infl[c][0.006][2]:
            continue
        for pir in (0.003, 0.006):
            rel[(c, pir)] = infl[c][pir][1] / infl["xyz"][pir][1]
            M.set(f"d11.{c}.rel.{pir}", rel[(c, pir)], ".1f")
    (cm, pm) = max(rel, key=rel.get)
    M.set("d11.max.rel", rel[(cm, pm)], ".1f")
    M.set("d11.max.rel.class", NAME[cm])
    M.set("d11.max.rel.pir", f"{pm:.3f}")
    fa_rel = [rel[("fa", p_)] for p_ in (0.003, 0.006)]
    M.set("d11.fa.rel.min", min(fa_rel), ".1f")
    M.set("d11.fa.rel.max", max(fa_rel), ".1f")
    M.set("d11.fa.over.xyz", infl["fa"][0.003][1] / infl["xyz"][0.003][1], ".1f")
    M.set("d11.fa.over.xyz.006", infl["fa"][0.006][1] / infl["xyz"][0.006][1], ".1f")
    M.set("d11.nf.over.xyz", infl["nf"][0.003][1] / infl["xyz"][0.003][1], ".1f")
    M.set("d11.o1t.over.xyz", infl["o1t"][0.003][1] / infl["xyz"][0.003][1], ".2f")


# ---------------------------------------------------------------------------
# D12 robustness of conclusions: router micro-architecture and network size
# ---------------------------------------------------------------------------
VARIANTS = [("base", "baseline"), ("v4", "4 VCs"), ("na", "non-atomic"), ("d4", "depth 4"), ("d16", "depth 16"),
            ("p4", "4-flit pkts"), ("p16", "16-flit pkts"), ("h3", "3-cycle hop")]


def d12(df: pd.DataFrame) -> None:
    sens = df[(df.group == "sens") & (df.load == "over")]
    mst = sens.groupby(["routing", "variant", "traffic"])["c.accepted_flits_per_node_cycle"].mean()
    dlk = sens.groupby(["routing", "variant"])["deadlocked"].sum()
    lines = [r"\begin{tabular}{@{}l" + "r" * len(VARIANTS) + r"ccc@{}}", r"\toprule",
             r"Class & \multicolumn{" + str(len(VARIANTS)) + r"}{c}{Accepted throughput at overload, uniform (flits/node/cycle)} & \multicolumn{2}{c}{Sign vs.\ ref.\ ($+/-$)} & D12 \\",
             r"\cmidrule(lr){2-" + str(1 + len(VARIANTS)) + r"}\cmidrule(lr){" + str(2 + len(VARIANTS)) + "-" +
             str(3 + len(VARIANTS)) + "}",
             " & " + " & ".join(v[1] for v in VARIANTS) + r" & Unif. & Transp. & \\", r"\midrule"]
    flips = []
    for fam in (FC, [c for c in VPC if c != "efr"]):
        for c in fam:
            ref = REF[c]
            cells, signs = [], {"random": [], "transpose1": []}
            for v, _ in VARIANTS:
                x = mst.get((c, v, "random"), np.nan)
                d = int(dlk.get((c, v), 0))
                cells.append(f"{x:.3f}" + (r"$^{\ddagger}$" if d else ""))
                M.set(f"d12.{c}.{v}", x, ".3f")
                for tr in ("random", "transpose1"):
                    a_ = mst.get((c, v, tr), np.nan)
                    b_ = mst.get((ref, v, tr), np.nan)
                    if np.isfinite(a_) and np.isfinite(b_) and c != ref:
                        signs[tr].append(float(a_ - b_))   # exact difference; an exact tie has no sign
            if c == ref:
                VERDICT[c]["D12"] = "ref"
                sg = ["--", "--"]
            else:
                RAW[c]["D12"] = {tr: list(signs[tr]) for tr in signs}
                consistent = all(len({int(np.sign(x)) for x in signs[tr] if x != 0}) <= 1
                                 for tr in ("random", "transpose1"))
                VERDICT[c]["D12"] = "pass" if consistent else "fail"
                sg = []
                for tr in ("random", "transpose1"):
                    npos = sum(1 for s in signs[tr] if s > 0)
                    nneg = sum(1 for s in signs[tr] if s < 0)
                    sg.append(f"{npos}/{nneg}")
                if not consistent:
                    flips.append(NAME[c])
                dlmark(c, "D12")
            lines.append(f"{cls(c)} & " + " & ".join(cells) + f" & {sg[0]} & {sg[1]} & {mark(VERDICT[c]['D12'])} \\\\")
        if fam is FC:
            lines.append(r"\midrule")
    VERDICT["efr"]["D12"] = "na"
    lines += [r"\bottomrule", r"\end{tabular}"]
    write("tab_d12.tex", "\n".join(lines) + "\n")
    M.set("d12.flips", ", ".join(flips))
    M.set("d12.nflips", len(flips))
    # Degenerate variant: depth 16 versus baseline (atomic reallocation, 8-flit packets).
    sall = df[df.group == "sens"]
    b = sall[sall.variant == "base"].set_index(["routing", "traffic", "load", "seed"])["c.accepted_flits_per_node_cycle"]
    d16 = sall[sall.variant == "d16"].set_index(["routing", "traffic", "load", "seed"])["c.accepted_flits_per_node_cycle"]
    j = pd.concat([b.rename("b"), d16.rename("d")], axis=1).dropna()
    M.set("d12.d16.same", int((j.b == j.d).sum()))
    M.set("d12.d16.n", len(j))
    M.set("d12.d16.diff.classes", ", ".join(sorted({NAME[i[0]] for i in j.index[j.b != j.d]})))
    # Noise level: seed-to-seed standard deviation of overload throughput.
    sd = sens.groupby(["routing", "variant", "traffic"])["c.accepted_flits_per_node_cycle"].std()
    M.set("d12.seed.sd.median", float(sd.median()), ".4f")
    # Scalability: overload throughput and zero-load latency per mesh size.
    sc = df[df.group == "scal"]
    sizes = [(4, 4, 4), (8, 8, 4), (8, 8, 8), (16, 16, 4)]
    over = sc[sc.load == "over"].groupby(["routing", "kx", "ky", "kz", "traffic"])[
        "c.accepted_flits_per_node_cycle"].mean()
    low = sc[sc.load == "low"].groupby(["routing", "kx", "ky", "kz"])["c.delay_mean"].mean()
    short = {(4, 4, 4): r"$4^3$", (8, 8, 4): r"$8^2{\times}4$", (8, 8, 8): r"$8^3$", (16, 16, 4): r"$16^2{\times}4$"}
    sz = " & ".join(short[x] for x in sizes)
    trs = ("random", "transpose1", "bitreversal")
    lines = [r"\begin{tabular}{@{}l" + "r" * (3 * len(sizes)) + r"@{}}", r"\toprule",
             r"Class & \multicolumn{4}{c}{Uniform} & \multicolumn{4}{c}{Transpose} & \multicolumn{4}{c}{Bit-reversal} \\",
             r"\cmidrule(lr){2-5}\cmidrule(lr){6-9}\cmidrule(lr){10-13}",
             " & " + " & ".join([sz] * 3) + r" \\", r"\midrule"]
    vals = {tr: {} for tr in trs}
    for c in FC:
        for tr in trs:
            for (a_, b_, c_) in sizes:
                vals[tr][(c, a_, b_, c_)] = over.get((c, a_, b_, c_, tr), np.nan)
    marks = {}
    for tr in trs:
        for (a_, b_, c_) in sizes:
            mk = rank_marks({c: vals[tr][(c, a_, b_, c_)] for c in FC}, lower_better=False, nd=3)
            for c in FC:
                marks[(c, tr, a_, b_, c_)] = mk[c]
    for c in FC:
        cells = []
        for tr in trs:
            for (a_, b_, c_) in sizes:
                v = vals[tr][(c, a_, b_, c_)]
                cells.append("--" if not np.isfinite(v) else marks[(c, tr, a_, b_, c_)])
                M.set(f"d12.scal.{c}.{tr}.{a_}{b_}{c_}", v, ".3f")
        lines.append(f"{cls(c)} & " + " & ".join(cells) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    write("tab_scal.tex", "\n".join(lines) + "\n")
    na = sens[(sens.variant == "na")]
    M.set("d12.na.dl.classes", ", ".join(NAME[c] for c in sorted(set(na[na.deadlocked].routing))))


# ---------------------------------------------------------------------------
# Paired statistics
# ---------------------------------------------------------------------------
def paired(df: pd.DataFrame) -> None:
    era = df[df.group == "era"]
    ev = df[df.group == "era_vpc"]
    dos = df[(df.group == "dos") & (df.attackers == 8)]
    keys = ["traffic", "pir", "seed"]

    def pairs(frame, a, b, metric, where=None):
        x = frame[frame.routing == a].set_index(keys)
        y = frame[frame.routing == b].set_index(keys)
        idx = x.index.intersection(y.index)
        if where is not None:
            idx = [i for i in idx if where(i)]
        ok = [i for i in idx if not x.loc[i, "deadlocked"] and not y.loc[i, "deadlocked"]]
        return x.loc[ok, metric].values.astype(float), y.loc[ok, metric].values.astype(float)

    nonhot = lambda i: i[0] != "hotspot"
    comps = [
        ("FA vs XYZ", "fa", "xyz", era, "c.delay_mean", lambda i: i[0] == "random"),
        ("FA vs XYZ", "fa", "xyz", era, "c.delay_mean", lambda i: i[0] == "hotspot"),
        ("ZXY vs XYZ", "zxy", "xyz", era, "c.delay_mean", lambda i: i[0] == "bitreversal"),
        ("QR vs FA", "qr", "fa", era, "c.delay_mean", nonhot),
        ("RCA vs FA", "rca", "fa", era, "c.delay_mean", nonhot),
        ("TA vs FA", "ta", "fa", era, "c.delay_mean", nonhot),
        ("HAM vs NF", "ham", "nf", era, "c.delay_mean", nonhot),
        ("O1T vs XYZ", "o1t", "xyz", dos, "c.delay_mean", None),
        ("FA vs XYZ", "fa", "xyz", dos, "c.delay_mean", None),
        ("AES-G vs EF", "aesg", "ef", ev, "c.delay_mean", nonhot),
        ("AdEle vs EF", "adele", "ef", ev, "c.delay_mean", nonhot),
        ("FL vs EF", "fl", "ef", ev, "c.delay_mean", nonhot),
    ]
    labels = ["uniform", "hotspot", "bit-reversal", "non-hotspot", "non-hotspot", "non-hotspot",
              "non-hotspot", "flooding DoS", "flooding DoS", "non-hotspot", "non-hotspot", "non-hotspot"]
    rows, ps = [], {}
    for (name, a, b, frame, metric, where), lab in zip(comps, labels):
        x, y = pairs(frame, a, b, metric, where)
        if len(x) < 6:
            continue
        p = wilcoxon_p(x, y)
        dlt = paired_delta(x, y)
        lo, hi = bootstrap_ci(paired_delta, x, y)
        med = float(np.median((x - y) / y))
        key = f"{name} ({lab})"
        ps[key] = p
        rows.append((key, len(x), med, dlt, lo, hi))
    adj = holm(ps)
    lines = [r"\begin{tabular}{@{}lrrrr@{}}", r"\toprule",
             r"Comparison (latency) & Pairs & Median rel.\ diff. & Paired $\delta$ [95\% CI] & Holm $p$ \\",
             r"\midrule"]
    for key, n, med, dlt, lo, hi in rows:
        lines.append(f"{key} & {n} & {100 * med:+.1f}\\% & {dlt:+.2f} [{lo:+.2f}, {hi:+.2f}] & {fmt_p(adj[key])} \\\\")
        tag = key.split(" (")[0].replace(" vs ", ".").replace("-", "").lower() + "." + key.split("(")[1].rstrip(")").replace(" ", "").replace("-", "")
        M.set(f"st.{tag}.med", 100 * med, "+.1f")
        M.set(f"st.{tag}.medabs", abs(100 * med), ".1f")
        M.set(f"st.{tag}.delta", dlt, "+.2f")
        M.set(f"st.{tag}.p", fmt_p(adj[key]))
        M.set(f"st.{tag}.n", n)
    lines += [r"\bottomrule", r"\end{tabular}"]
    write("tab_stats.tex", "\n".join(lines) + "\n")
    M.set("st.n.tests", len(rows))


# ---------------------------------------------------------------------------
# Verdict
# ---------------------------------------------------------------------------
DIMS = ["D1", "D2", "D3", "D4", "D5", "D6", "D8", "D9", "D10", "D11", "D12"]


def verdict() -> None:
    lines = [r"\begin{tabular}{@{}l" + "c" * len(DIMS) + r"cr@{}}", r"\toprule",
             "Class & " + " & ".join(r"\D{" + d[1:] + "}" for d in DIMS) + r" & Passed & VCs \\", r"\midrule"]
    for fam in (FC, VPC):
        for c in fam:
            cells = [mark(VERDICT[c].get(d, "na")) for d in DIMS]
            npass = sum(1 for d in DIMS if VERDICT[c].get(d) == "pass")
            nall = sum(1 for d in DIMS if VERDICT[c].get(d) in ("pass",) + FAIL_LIKE)
            lines.append(f"{cls(c)} & " + " & ".join(cells) + f" & {npass}/{nall} & {MIN_VCS[c]} \\\\")
            M.set(f"vd.{c}.pass", npass)
            M.set(f"vd.{c}.n", nall)
        if fam is FC:
            lines.append(r"\midrule")
    lines += [r"\bottomrule", r"\end{tabular}"]
    write("tab_verdict.tex", "\n".join(lines) + "\n")
    allpass = [c for c in ALL if all(VERDICT[c].get(d) in ("pass", "ref", "na") for d in DIMS)]
    M.set("vd.allpass", ", ".join(NAME[c] for c in allpass) if allpass else "none")
    best = max(ALL, key=lambda c: sum(1 for d in DIMS if VERDICT[c].get(d) == "pass"))
    M.set("vd.best", NAME[best])
    mx = max(sum(1 for d in DIMS if VERDICT[c].get(d) == "pass") for c in ALL)
    M.set("vd.top", ", ".join(NAME[c] for c in ALL if sum(1 for d in DIMS if VERDICT[c].get(d) == "pass") == mx))
    M.set("vd.topn", mx)
    M.set("dlrule.n", len(DL_CHANGED))
    M.set("dlrule.list", "; ".join(f"{NAME[c]} {d}" for c, d, _ in DL_CHANGED))
    (GEN / "verdict.json").write_text(json.dumps(VERDICT, indent=1))
    (GEN / "deadlock_rule.json").write_text(json.dumps(DL_CHANGED))


# ---------------------------------------------------------------------------
# Sensitivity of the verdict to the pre-registered thresholds
# ---------------------------------------------------------------------------
# (strict, registered, lenient): every tolerance halved / doubled; criteria
# whose registered tolerance is zero receive a fixed lenient tolerance.
THRESH = {
    "D2": ((0.75,), (0.50,), (0.25,)),
    "D3": ((0.95,), (0.90,), (0.80,)),
    "D4": ((1.5,), (2.0,), (3.0,)),
    "D5": ((0.10, 0.0), (0.05, 0.0), (0.025, 0.01)),
    "D6": ((0.0, 0.975), (0.0, 0.95), (0.01, 0.90)),
    "D8": ((0.025,), (0.05,), (0.10,)),
    "D9": ((1.025,), (1.05,), (1.10,)),
    "D10": ((0.025,), (0.05,), (0.10,)),
    "D11": ((1.0,), (1.0,), (1.25,)),
    "D12": ((0.0,), (0.0,), (0.01,)),
}
THRESH_TXT = {
    "D2": ("pairs with $>1$ path", r"$\ge 75\%$", r"$\ge 50\%$", r"$\ge 25\%$"),
    "D3": ("uniform sat.\\ / ref.", r"$\ge 0.95$", r"$\ge 0.90$", r"$\ge 0.80$"),
    "D4": ("worst latency ratio", r"$\le 1.5$", r"$\le 2$", r"$\le 3$"),
    "D5": ("reach.\\ gain; VPC loss", r"$\ge 0.10$; $0$", r"$\ge 0.05$; $0$", r"$\ge 0.025$; $\le 0.01$"),
    "D6": ("source Jain", r"$\ge 0.975$", r"$\ge 0.95$", r"$\ge 0.90^{a}$"),
    "D8": ("area overhead", r"$\le 2.5\%$", r"$\le 5\%$", r"$\le 10\%$"),
    "D9": ("energy ratio (max)", r"$\le 1.025$", r"$\le 1.05$", r"$\le 1.10$"),
    "D10": (r"$\Delta T_{\max}$ vs.\ ref.", r"$\le 0.025$\,K", r"$\le 0.05$\,K", r"$\le 0.10$\,K"),
    "D11": ("inflation / ref.", r"$\le 1$", r"$\le 1$", r"$\le 1.25$"),
    "D12": ("sign of throughput diff.", "exact", "exact", r"ties $<0.01^{b}$"),
}


def judge(c: str, d: str, th: tuple) -> str:
    """Re-evaluate one verdict from its raw statistic under thresholds th."""
    v = VERDICT[c].get(d)
    if v in ("ref", "na", "dl") or d not in RAW[c]:
        return v or "na"
    r = RAW[c][d]
    if d == "D2":
        return "pass" if r >= th[0] else "fail"
    if d == "D3":
        return "pass" if (r[0] >= th[0] and r[1]) else "fail"
    if d == "D4":
        return "pass" if r <= th[0] else "fail"
    if d == "D5":
        kind, x = r
        if kind == "fc":
            if all(g >= th[0] for g in x):
                return "pass"
            return "partial" if any(g > 1e-3 for g in x) else "fail"
        return "pass" if x <= th[1] + 1e-12 else "fail"
    if d == "D6":
        dj, sj = r
        if c == REF[c]:
            return v
        return "pass" if (dj >= -th[0] - 1e-12 and sj >= th[1]) else "fail"
    if d == "D8":
        over, glob = r
        return "pass" if (over <= th[0] and not glob) else "fail"
    if d == "D9":
        return "pass" if r <= th[0] else "fail"
    if d == "D10":
        return "pass" if r <= th[0] + 1e-12 else "fail"
    if d == "D11":
        return "pass" if r <= th[0] + 1e-12 else "fail"
    if d == "D12":
        ok = all(len({int(np.sign(x)) for x in r[tr] if abs(x) > th[0]}) <= 1 for tr in r)
        return "pass" if ok else "fail"
    return v


def threshold_sensitivity() -> None:
    regimes = ("strict", "registered", "lenient")
    # consistency check: the registered regime must reproduce the verdict
    for c in ALL:
        for d in THRESH:
            got = judge(c, d, THRESH[d][1])
            if VERDICT[c].get(d) not in ("ref", "na", "dl") and d in RAW[c]:
                assert got == VERDICT[c][d], (c, d, got, VERDICT[c][d])
    lines = [r"\begin{tabular}{@{}llccccl@{}}", r"\toprule",
             r"Criterion & Statistic & Strict & Registered & Lenient & Passes (S/R/L) & Changes from registered \\",
             r"\midrule"]
    totals = {rg: defaultdict(int) for rg in regimes}
    nflip = 0
    for d in THRESH:
        counts, gained, lost = [], [], []
        for i, rg in enumerate(regimes):
            n = 0
            for c in ALL:
                v = judge(c, d, THRESH[d][i])
                if v == "pass":
                    n += 1
                    totals[rg][c] += 1
            counts.append(n)
        for c in ALL:
            reg = judge(c, d, THRESH[d][1])
            if judge(c, d, THRESH[d][2]) == "pass" and reg != "pass":
                gained.append(NAME[c])
            if judge(c, d, THRESH[d][0]) != "pass" and reg == "pass":
                lost.append(NAME[c])
        nflip += len(gained) + len(lost)
        ch = []
        if lost:
            ch.append(r"S: $-$" + ", ".join(lost))
        if gained:
            ch.append(r"L: $+$" + ", ".join(gained))
        t = THRESH_TXT[d]
        lines.append(rf"\D{{{d[1:]}}} & {t[0]} & {t[1]} & {t[2]} & {t[3]} & {counts[0]}/{counts[1]}/{counts[2]} & "
                     + ("; ".join(ch) if ch else "none") + r" \\")
        M.set(f"th.{d}.counts", f"{counts[0]}/{counts[1]}/{counts[2]}")
    lines += [r"\bottomrule", r"\end{tabular}"]
    write("tab_thresh.tex", "\n".join(lines) + "\n")
    M.set("th.nflip", nflip)
    for rg in regimes:
        for c in ALL:
            base = sum(1 for d in DIMS if d not in THRESH and VERDICT[c].get(d) == "pass")
            M.set(f"vd.{c}.pass.{rg}", base + totals[rg][c])
        top = max(ALL, key=lambda c: (sum(1 for d in DIMS if d not in THRESH and VERDICT[c].get(d) == "pass")
                                      + totals[rg][c]))
        mx = sum(1 for d in DIMS if d not in THRESH and VERDICT[top].get(d) == "pass") + totals[rg][top]
        tops = [NAME[c] for c in ALL if (sum(1 for d in DIMS if d not in THRESH and VERDICT[c].get(d) == "pass")
                                         + totals[rg][c]) == mx]
        M.set(f"th.top.{rg}", ", ".join(tops))
        M.set(f"th.topn.{rg}", mx)


# ---------------------------------------------------------------------------
# Partition risk under correlated TSV faults (analytical / Monte Carlo)
# ---------------------------------------------------------------------------
def _sci(p: float, ns: int = 2_000_000, exact: bool = False) -> str:
    """Scientific notation. A Monte Carlo zero is reported as the 95 % upper bound 3/ns (rule of
    three); a Monte Carlo value resting on fewer than 20 events gets one significant digit."""
    if p <= 0:
        b = 3 / ns
        e = int(np.floor(np.log10(b)))
        return rf"$<{b / 10 ** e:.1f}{{\times}}10^{{{e}}}$"
    e = int(np.floor(np.log10(p)))
    m = p / 10 ** e
    few = (not exact) and p * ns < 20
    if e >= -2:
        return f"{p:.3f}" if not few else f"{p:.0e}"
    if few:
        m = round(m)
        if m == 10:
            m, e = 1, e + 1
        return rf"${m:d}{{\times}}10^{{{e}}}$"
    return rf"${m:.1f}{{\times}}10^{{{e}}}$"


def partition() -> None:
    res = json.loads((ROOT / "analysis3" / "partition.json").read_text())
    lines = [r"\begin{tabular}{@{}lrrrrrrr@{}}", r"\toprule",
             r"Elevators per boundary & $p_v$ & \multicolumn{4}{c}{$\Prob(\text{boundary severed})$} & "
             r"\multicolumn{2}{c}{s.d.\ of failed fraction} \\",
             r"\cmidrule(lr){3-6}\cmidrule(lr){7-8}",
             r" & & $\xi=0$ & $\xi=2$ & $\xi=4$ & $\xi=8$ & $\xi=0$ & $\xi=4$ \\", r"\midrule"]
    names = {"full": "64 (fully connected)", "stride2": "16 (stride 2)", "stride4": "4 (stride 4)"}
    for pl in ("full", "stride2", "stride4"):
        for pv in (0.1, 0.2, 0.3):
            r = {xi: res[f"{pl}|{pv:g}|{xi:g}"] for xi in (0, 2, 4, 8)}
            cells = [_sci(r[xi]["p_cut"], exact=(xi == 0)) for xi in (0, 2, 4, 8)]
            sds = [f"{r[xi]['sd_frac']:.3f}" for xi in (0, 4)]
            lab = names[pl] if pv == 0.1 else ""
            lines.append(f"{lab} & {pv:.1f} & " + " & ".join(cells) + " & " + " & ".join(sds) + r" \\")
            for xi in (0, 2, 4, 8):
                M.set(f"part.{pl}.{pv:g}.{xi}", _sci(r[xi]["p_cut"], exact=(xi == 0)))
                M.set(f"part3.{pl}.{pv:g}.{xi}", _sci(r[xi]["p_cut3"], exact=(xi == 0)))
            M.set(f"part.sd.{pl}.{pv:g}.0", r[0]["sd_frac"], ".3f")
            M.set(f"part.sd.{pl}.{pv:g}.4", r[4]["sd_frac"], ".3f")
        if pl != "stride4":
            lines.append(r"\addlinespace")
    lines += [r"\bottomrule", r"\end{tabular}"]
    write("tab_partition.tex", "\n".join(lines) + "\n")
    ratio = res["stride4|0.2|4"]["p_cut"] / res["stride4|0.2|0"]["p_cut"]
    M.set("part.stride4.ratio", ratio, ".0f")
    r4 = [res[f"stride4|{pv:g}|4"]["p_cut"] / res[f"stride4|{pv:g}|0"]["p_cut"] for pv in (0.1, 0.2, 0.3)]
    M.set("part.ratio.stride4.min", min(r4), ".1f")
    M.set("part.ratio.stride4.max", max(r4), ".0f")
    for pl in ("full", "stride2", "stride4"):
        sdr = [res[f"{pl}|{pv:g}|4"]["sd_frac"] / res[f"{pl}|{pv:g}|0"]["sd_frac"] for pv in (0.1, 0.2, 0.3)]
        M.set(f"part.sdr.{pl}.min", min(sdr), ".1f")
        M.set(f"part.sdr.{pl}.max", max(sdr), ".1f")


# ---------------------------------------------------------------------------
# VPC with four elevators (stride 4)
# ---------------------------------------------------------------------------
def vpc4(df: pd.DataFrame) -> None:
    v4 = df[df.group == "vpc4"]
    g = v4.groupby(["routing", "traffic", "pir"]).agg(lat=("c.delay_mean", "mean"),
                                                     und=("c.undelivered", "mean")).reset_index()
    pirs = sorted(g.pir.unique())
    classes = [c for c in VPC if c in set(g.routing)]
    lines = [r"\begin{tabular}{@{}l" + "r" * (2 * len(pirs)) + r"@{}}", r"\toprule",
             r"Class & \multicolumn{" + str(len(pirs)) + r"}{c}{Uniform} & \multicolumn{" + str(len(pirs)) +
             r"}{c}{Transpose} \\",
             r"\cmidrule(lr){2-" + str(1 + len(pirs)) + r"}\cmidrule(lr){" + str(2 + len(pirs)) + "-" +
             str(1 + 2 * len(pirs)) + "}",
             r"PIR & " + " & ".join(f"{p:.3f}" for p in pirs) * 1 + " & " + " & ".join(f"{p:.3f}" for p in pirs) +
             r" \\", r"\midrule"]
    marks = {}
    for tr in ("random", "transpose1"):
        for p in pirs:
            sub = g[(g.traffic == tr) & (g.pir == p)]
            vals = {c: float(sub[sub.routing == c].lat.iloc[0]) for c in classes
                    if float(sub[sub.routing == c].und.iloc[0]) <= 0.01}
            mk = rank_marks(vals, lower_better=True, nd=1) if len(vals) > 1 else {}
            for c in classes:
                if c not in mk:
                    mk[c] = f"{float(sub[sub.routing == c].lat.iloc[0]):.1f}"
            marks[(tr, p)] = mk
    for c in classes:
        cells = []
        for tr in ("random", "transpose1"):
            for p in pirs:
                r = g[(g.routing == c) & (g.traffic == tr) & (g.pir == p)].iloc[0]
                s = marks[(tr, p)][c]
                if r.und > 0.01:
                    s += rf"$^{{\ast}}$"
                cells.append(s)
                M.set(f"vpc4.{c}.{tr}.{p:g}", r.lat, ".1f")
                M.set(f"vpc4.{c}.{tr}.{p:g}.und", 100 * r.und, ".0f")
        lines.append(f"{cls(c)} & " + " & ".join(cells) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    write("tab_vpc4.tex", "\n".join(lines) + "\n")
    und6 = g[g.pir == max(pirs)].und
    M.set("vpc4.und.min", 100 * und6.min(), ".0f")
    M.set("vpc4.und.max", 100 * und6.max(), ".0f")


# ---------------------------------------------------------------------------
# Exact check of the VPC fault path: Elevator-First loss given each fault map
# ---------------------------------------------------------------------------
def faultcheck() -> None:
    s = json.loads((ROOT / "sim3" / "experiments" / "fault_check" / "summary.json").read_text())
    for pv in ("0.1", "0.2"):
        for xi in ("0.0", "4.0"):
            r = s[f"{pv}|{xi}"]
            M.set(f"fcheck.{pv}.{xi[0]}.mean", r["loss_mean"], ".3f")
            M.set(f"fcheck.{pv}.{xi[0]}.sd", r["loss_sd"], ".3f")
            M.set("fcheck.maps", r["maps"])
    M.set("fcheck.run.max", s["per_run_max_absdiff"], ".4f")


# ---------------------------------------------------------------------------
# Out-of-order delivery per flow (side effect of path diversity)
# ---------------------------------------------------------------------------
def reorder() -> None:
    f = ROOT / "sim3" / "experiments" / "reorder" / "reorder.json"
    if not f.exists():
        return
    d = pd.DataFrame(json.loads(f.read_text()))
    g = d.groupby(["routing", "traffic", "pir"])["ooo"].mean()
    classes = [c for c in ALL if c in set(d.routing)]
    cols = [("transpose1", 0.006), ("transpose1", 0.012), ("bitreversal", 0.006), ("bitreversal", 0.012)]
    lines = [r"\begin{tabular}{@{}lrrrr@{}}", r"\toprule",
             r"Class & \multicolumn{2}{c}{Transpose} & \multicolumn{2}{c}{Bit-reversal} \\",
             r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}", r"PIR & 0.006 & 0.012 & 0.006 & 0.012 \\", r"\midrule"]
    for fam in (FC, VPC):
        rows = [c for c in fam if c in classes]
        for c in rows:
            vals = [100 * g.get((c, tr, p), np.nan) for tr, p in cols]
            lines.append(f"{cls(c)} & " + " & ".join("--" if not np.isfinite(v) else f"{v:.2f}" for v in vals) + r" \\")
            for (tr, p), v in zip(cols, vals):
                M.set(f"ooo.{c}.{tr}.{p:g}", v, ".2f")
        if fam is FC and any(c in classes for c in VPC):
            lines.append(r"\midrule")
    lines += [r"\bottomrule", r"\end{tabular}"]
    write("tab_reorder.tex", "\n".join(lines) + "\n")
    mx = d.groupby("routing")["ooo"].max()
    for c in classes:
        M.set(f"ooo.{c}.max", 100 * mx[c], ".2f")
    top = max(classes, key=lambda c: mx[c])
    M.set("ooo.max.class", NAME[top])
    M.set("ooo.max.v", 100 * mx[top], ".1f")
    inorder = [NAME[c] for c in classes if mx[c] == 0]
    M.set("ooo.inorder", ", ".join(inorder))
    M.set("ooo.inorder.n", len(inorder))
    adaptive_fc = [c for c in ("nf", "dyad", "ham", "fa", "rca", "qr", "ta") if c in classes]
    M.set("ooo.adapt.maxv", 100 * max(mx[c] for c in adaptive_fc), ".2f")
