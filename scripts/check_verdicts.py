"""Independent check of the verdict table (Table "Verdict" of the article).

Recomputes, directly from the raw results and without the article's generator (paper3/sec_results.py), the
verdicts that depend on simulation runs: D1 (checker verdict and deadlocks in every run), D3 (saturation over the
coarse sweep and its bisection refinement), the deadlock rule on the run set of every simulated dimension
(D4, D6, D9, D10, D11, D12), the exact-sign test of D12, and the pass counts of the table.
Usage (from the repository root): python3 scripts/check_verdicts.py
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis3"))
from load import load
df = pd.concat([load(ROOT/"sim3/results/campaign3.jsonl"), load(ROOT/"sim3/results/sweep_refine.jsonl")], ignore_index=True)
V = json.loads((ROOT/"paper3/gen/verdict.json").read_text())
FC = ["xyz","zxy","nf","dyad","ham","o1t","romm","val","fa","rca","qr","ta"]
VPC = ["ef","efr","aesg","adele","fl"]
errs = []
# D1
for c in FC + VPC:
    sub = "vpc_stride2_" + c if c in VPC else c
    d = json.loads((ROOT/f"sim3/results/cdg_{sub}.json").read_text())
    proven = (not d["cyclic"]) or (d["escape"] and not d["escape_cyclic"] and d["escape_to_adaptive"] == 0 and d["adaptive_states_without_escape"] == 0)
    k = int(df[df.routing == c].deadlocked.sum())
    na = df[(df.routing == c) & (df.atomic == 0)]
    exp = "pass" if proven and k == 0 else ("cond" if proven and int(df[(df.routing==c)&(df.atomic==1)].deadlocked.sum()) == 0 else "fail")
    if V[c]["D1"] != exp: errs.append(("D1", c, V[c]["D1"], exp))
# D3 saturation
for fam, grps, classes, ref in (("fc", ["sweep","sweepref"], FC, "xyz"), ("vpc", ["sweep_vpc","sweepref_vpc"], [c for c in VPC if c != "efr"], "ef")):
    sw = df[df.group.isin(grps)]
    m = sw.groupby(["routing","traffic","pir"]).agg(acc=("c.accepted_flits_per_node_cycle","mean"), und=("c.undelivered","mean"), dl=("deadlocked","sum")).reset_index()
    ok = m[(m.und <= 0.01) & (m.dl == 0)]
    sat = ok.groupby(["routing","traffic"]).acc.max().unstack()
    for c in classes:
        if c == ref: continue
        kdl = int(sw[sw.routing == c].deadlocked.sum())
        p = sat.loc[c,"random"] >= 0.9*sat.loc[ref,"random"] and (sat.loc[c,"transpose1"] > sat.loc[ref,"transpose1"] or sat.loc[c,"bitreversal"] > sat.loc[ref,"bitreversal"])
        exp = "dl" if kdl else ("pass" if p else "fail")
        if V[c]["D3"] != exp: errs.append(("D3", c, V[c]["D3"], exp))
# deadlock rule per dimension
era = df.group.isin(["era","era_vpc"])
sets = {"D4": era, "D6": (era & (df.traffic=="random")) | ((df.group=="stress") & (df.traffic=="random")),
        "D10": era & df.traffic.isin(["random","hotspot"]) & (df.pir==0.006), "D11": df.group=="dos",
        "D12": (df.group=="sens") & (df.load=="over"), "D9": era & (df.traffic=="random") & (df.pir==0.006)}
for d, mask in sets.items():
    k = df[mask].groupby("routing").deadlocked.sum()
    for c in FC + VPC:
        if V[c].get(d) in ("ref","na"): continue
        if k.get(c, 0) > 0 and V[c].get(d) != "dl": errs.append((d, c, V[c].get(d), "dl"))
        if k.get(c, 0) == 0 and V[c].get(d) == "dl": errs.append((d, c, "dl", "not dl"))
# D12 exact signs (non-dl classes)
sens = df[(df.group=="sens") & (df.load=="over")]
mst = sens.groupby(["routing","variant","traffic"])["c.accepted_flits_per_node_cycle"].mean()
for c in FC + [x for x in VPC if x != "efr"]:
    ref = "ef" if c in VPC else "xyz"
    if c == ref or V[c]["D12"] == "dl": continue
    ok = True
    for tr in ("random","transpose1"):
        s = {int(np.sign(mst[(c,v,tr)] - mst[(ref,v,tr)])) for v in ["base","v4","na","d4","d16","p4","p16","h3"]} - {0}
        ok &= len(s) <= 1
    if V[c]["D12"] != ("pass" if ok else "fail"): errs.append(("D12", c, V[c]["D12"], ok))
# pass counts in the verdict table
tab = (ROOT/"paper3/gen/tab_verdict.tex").read_text()
D = ["D1","D2","D3","D4","D5","D6","D8","D9","D10","D11","D12"]
NAME = {"xyz":"XYZ","zxy":"ZXY","nf":"NF","dyad":"DyAD","ham":"HAM","o1t":"O1T","romm":"ROMM","val":"VAL","fa":"FA","rca":"RCA","qr":"QR","ta":"TA","ef":"EF","efr":"EF-R","aesg":"AES-G","adele":"AdEle","fl":"FL"}
for c in FC + VPC:
    n = sum(V[c].get(d) == "pass" for d in D)
    line = [l for l in tab.splitlines() if l.startswith(f"\\cls{{{NAME[c]}}} &")][0]
    if f" & {n}/" not in line: errs.append(("count", c, line[-20:], n))
print("verdict table consistent with the raw results" if not errs else f"MISMATCHES: {errs}")
sys.exit(1 if errs else 0)
