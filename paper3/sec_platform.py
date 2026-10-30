"""Counts, simulator cross-validation and channel-dependency verdicts."""
from __future__ import annotations

import json
from collections import defaultdict

import numpy as np
import pandas as pd

from common import ALL, FC, GEN, M, NAME, ROOT, VPC, cls, write
from load import load

# ---------------------------------------------------------------------------
def load_all() -> pd.DataFrame:
    df = load(ROOT / "sim3" / "results" / "campaign3.jsonl")
    ref = ROOT / "sim3" / "results" / "sweep_refine.jsonl"
    if ref.exists() and ref.stat().st_size:
        df = pd.concat([df, load(ref)], ignore_index=True)
    return df


def section_counts(df: pd.DataFrame, val: list, fr: pd.DataFrame) -> None:
    M.set("n.sim", len(df))
    M.set("n.val", len(val))
    M.set("n.grand", len(df) + len(val))
    M.set("n.traced", int(df["e.MI"].notna().sum()))
    confs = fr[["p_h", "p_v", "xi"]].drop_duplicates()
    faulty = confs[(confs.p_h > 0) | (confs.p_v > 0)]
    M.set("n.faultmaps", 2000 * len(faulty))          # maps with at least one fault probability > 0
    M.set("fr.npol", fr["policy"].nunique() - 1)
    M.set("fr.nconf", len(confs))
    M.set("fr.nconf.faulty", len(faulty))
    M.set("cls.n", len(ALL))
    M.set("cls.nfc", len(FC))
    M.set("cls.nvpc", len(VPC))
    for g, n in df.groupby("group").size().items():
        M.set(f"n.group.{g}", int(n))
    M.set("n.group.sweepref.all", int(df.group.isin(["sweepref", "sweepref_vpc"]).sum()))
    M.set("n.groups", int(df.group.replace({"sweepref_vpc": "sweepref"}).nunique()))
    M.set("n.deadlockruns", int(df["deadlocked"].sum()))
    # Supplementary experiments outside the campaign groups.
    ex = ROOT / "sim3" / "experiments"
    nlines = lambda f: sum(1 for line in f.read_text().splitlines() if line.strip()) if f.exists() else 0
    n_re = len(json.loads((ex / "reorder" / "reorder.json").read_text()))
    n_in = nlines(ex / "inject_policy" / "res.txt")
    n_sz = nlines(ex / "inject_policy" / "res_sizes.txt")
    n_pi = nlines(ex / "pipeline" / "res.txt") + nlines(ex / "pipeline" / "res2.txt")
    n_rs = nlines(ex / "rand_selection" / "res.txt")
    n_co = len(list((ex / "commitment").glob("e_*.txt")))
    bs = json.loads((ROOT / "validation" / "nonatomic" / "summary.json").read_text())
    n_bs = len(bs["booksim_batch_nonatomic"]) + len(list((ex / "commitment").glob("instr*.err")))
    for k, v in (("reorder", n_re), ("inject", n_in), ("sizes", n_sz), ("pipe", n_pi), ("rsel", n_rs),
                 ("commit.ours", n_co), ("bs.n", n_bs)):
        M.set(f"supp.{k}", v)
    M.set("supp.total.ours", n_re + n_in + n_sz + n_pi + n_rs + n_co)


# ---------------------------------------------------------------------------
def section_validation(val: list) -> None:
    """BookSim 2 cross-validation. Saturation load: highest offered load whose
    mean latency stays below three times the zero-load latency, applied
    identically to both simulators."""
    by = defaultdict(lambda: defaultdict(list))
    for r in val:
        by[(r["k"], r["routing"], r["rate"])][r["sim"]].append(r)
    offsets, tables = [], []
    for k in (4, 8):
        for alg in ("xyz", "romm", "val", "fa"):
            rates = sorted({key[2] for key in by if key[0] == k and key[1] == alg})

            def lat(sim, rt):
                xs = [x["latency"] for x in by[(k, alg, rt)][sim]]
                if any(x is None for x in xs):
                    return float("inf")
                return float(np.mean(xs))

            z_b, z_o = lat("booksim", rates[0]), lat("ours", rates[0])

            def sat(sim, z):
                ok = [rt for rt in rates if lat(sim, rt) <= 3 * z]
                return max(ok)

            s_b, s_o = sat("booksim", z_b), sat("ours", z_o)
            offsets.append(z_b - z_o)
            lim = 0.8 * min(s_b, s_o)
            dev = [abs((lat("booksim", rt) - z_b) - (lat("ours", rt) - z_o)) / lat("ours", rt)
                   for rt in rates if rt <= lim]
            tables.append(dict(k=k, alg=alg, zb=z_b, zo=z_o, off=z_b - z_o, sb=s_b, so=s_o,
                               dev=max(dev) if len(dev) > 1 else float("nan")))
    lines = [r"\begin{tabular}{@{}llrrrrrr@{}}", r"\toprule",
             r"Mesh & Routing & \multicolumn{2}{c}{Zero-load latency} & Offset & \multicolumn{2}{c}{Saturation load} & Max.\ dev. \\",
             r"\cmidrule(lr){3-4}\cmidrule(lr){6-7}",
             r" & & BookSim & ours & (cycles) & BookSim & ours & (\%) \\", r"\midrule"]
    names = {"xyz": "DOR (XYZ)", "romm": "ROMM", "val": "Valiant", "fa": "Min.\\ adaptive + escape"}
    for t in tables:
        dev = "--" if not np.isfinite(t["dev"]) else f"{100 * t['dev']:.1f}"
        lines.append(f"${t['k']}\\times{t['k']}\\times{t['k']}$ & {names[t['alg']]} & {t['zb']:.1f} & {t['zo']:.1f} & "
                     f"{t['off']:.1f} & {t['sb']:.2f} & {t['so']:.2f} & {dev} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    write("tab_validation.tex", "\n".join(lines) + "\n")
    M.set("val.off.min", min(offsets), ".1f")
    M.set("val.off.max", max(offsets), ".1f")
    agree = sum(1 for t in tables if abs(t["sb"] - t["so"]) < 1e-9)
    M.set("val.sat.agree", agree)
    M.set("val.sat.n", len(tables))
    nonfa = [t["dev"] for t in tables if t["alg"] != "fa" and np.isfinite(t["dev"])]
    fa = [t["dev"] for t in tables if t["alg"] == "fa" and np.isfinite(t["dev"])]
    M.set("val.nonfa.dev", 100 * max(nonfa), ".1f")
    M.set("val.fa.dev", 100 * max(fa), ".1f")
    mism = [t for t in tables if abs(t["sb"] - t["so"]) >= 1e-9]
    M.set("val.mismatch", ", ".join(f"{t['alg']}@{t['k']}" for t in mism) if mism else "none")


# ---------------------------------------------------------------------------
def section_cdg() -> dict:
    res = {}
    for c in FC:
        res[c] = json.loads((ROOT / "sim3" / "results" / f"cdg_{c}.json").read_text())
    for c in VPC:
        res[c] = json.loads((ROOT / "sim3" / "results" / f"cdg_vpc_stride2_{c}.json").read_text())
        res[c + "@4"] = json.loads((ROOT / "sim3" / "results" / f"cdg_vpc_stride4_{c}.json").read_text())
    lines = [r"\begin{tabular}{@{}lrrcccl@{}}", r"\toprule",
             r"Class & Channels & Dependencies & Cyclic & Escape acyclic & Shortest cycle & Verdict \\",
             r"\midrule"]
    for c in ALL:
        d = res[c]
        if not d["cyclic"]:
            verdict = r"deadlock-free (acyclic CDG)"
        elif d["escape"] and not d["escape_cyclic"] and d["escape_to_adaptive"] == 0 \
                and d["adaptive_states_without_escape"] == 0:
            verdict = r"deadlock-free (Duato)"
        else:
            verdict = r"\textbf{deadlock possible}"
        esc = (r"\cmark" if not d["escape_cyclic"] else r"\xmark") if d["escape"] else "--"
        lines.append(f"{cls(c)} & {d['channels_with_deps']:,} & {d['edges']:,} & "
                     f"{r'\cmark' if d['cyclic'] else r'\xmark'} & {esc} & "
                     f"{d['cycle_length'] if d['cyclic'] else '--'} & {verdict} \\\\".replace(",", "{,}"))
        M.set(f"cdg.{c}.edges", d["edges"])
        M.set(f"cdg.{c}.chan", d["channels_with_deps"])
    lines += [r"\bottomrule", r"\end{tabular}"]
    write("tab_cdg.tex", "\n".join(lines) + "\n")
    neg = json.loads((ROOT / "sim3" / "results" / "cdg_negative_ef_mixed.json").read_text())
    M.set("cdg.neg.len", neg["cycle_length"])
    M.set("cdg.fa.len", res["fa"]["cycle_length"])
    M.set("cdg.ham.nonmin", res["ham"]["nonminimal_moves"])
    proven = sum(1 for c in ALL if not res[c]["cyclic"] or (res[c]["escape"] and not res[c]["escape_cyclic"]))
    M.set("cdg.proven", proven)
    M.set("cdg.states.max", max(res[c]["states"] for c in ALL))
    return res




# ---------------------------------------------------------------------------
def section_nonatomic() -> None:
    """Cross-simulator check of the non-atomic deadlock: BookSim batch-mode drain
    tests and the injection-policy experiment of the in-house simulator."""
    s = json.loads((ROOT / "validation" / "nonatomic" / "summary.json").read_text())
    bs = s["booksim_batch_nonatomic"]
    M.set("na.bs.n", len(bs))
    M.set("na.bs.dl", sum(1 for r in bs if r["warn"] > 0 or not r["completed"]))
    cyc = [r.get("cycles") for r in bs if r.get("cycles")]
    if cyc:
        M.set("na.bs.maxcycles", max(cyc))
    inj = s["inject_policy"]
    for ie in (0, 1):
        x = [r for r in inj if r["inject_escape"] == ie]
        M.set(f"na.inj{ie}.k", sum(1 for r in x if not r["drained"]))
        M.set(f"na.inj{ie}.n", len(x))
        st = [r["stuck"] for r in x if not r["drained"]]
        if st:
            M.set(f"na.inj{ie}.stuck.max", max(st))
    for h in (3, 4):
        for size in ("8x8x4", "8x8x8"):
            x = [r for r in s["pipeline_default_injection"] if r["hop_delay"] == h and r.get("size", "8x8x4") == size]
            if x:
                M.set(f"na.pipe.h{h}.{size}", f"{sum(1 for r in x if not r['drained'])}/{len(x)}")
    for size in ("8x8x4", "8x8x8"):
        x = [r for r in s.get("equal_priority_selection", []) if r["size"] == size]
        if x:
            M.set(f"na.rsel.{size}", f"{sum(1 for r in x if not r['drained'])}/{len(x)}")
    cr = s.get("commitment_rates", {})
    for k, v in cr.items():
        M.set(f"na.commit.{k.replace('_', '.')}", v["pct_full_adaptive"], ".1f")
    sz = s["sizes_default_injection"]
    for size in ("4x4x4", "8x8x8", "8x8x4"):
        x = [r for r in sz if r["size"] == size]
        M.set(f"na.size.{size}", f"{sum(1 for r in x if not r['drained'])}/{len(x)}")
