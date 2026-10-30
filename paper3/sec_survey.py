"""Corpus tables and the audit of evaluation practice."""
from __future__ import annotations

from collections import Counter

from common import M, write
from corpus import CORPUS, HISTORICAL

DL_TXT = {"F": "formal", "R": "inherited", "C": "claimed", "I": "mitigated", "N": "none", "n/v": "n/v"}
NOXIM = ("Access Noxim", "Noxim (ext.)", "Noxim")


def esc(s: str) -> str:
    return s.replace("&", r"\&").replace("%", r"\%").replace("ChipletSimulator", r"Chiplet\-Simulator")


def corpus_tables() -> None:
    rows = sorted(CORPUS, key=lambda d: (d["year"], d["name"].lower()))
    lines = [r"\begin{tabularx}{\linewidth}{@{}>{\raggedright\arraybackslash}p{2.35cm}c>{\raggedright\arraybackslash}p{1.9cm}"
             r">{\raggedright\arraybackslash}p{1.35cm}>{\raggedright\arraybackslash}X>{\raggedright\arraybackslash}p{2.0cm}"
             r">{\raggedright\arraybackslash}p{1.3cm}c>{\raggedright\arraybackslash}p{1.45cm}c@{}}",
             r"\toprule",
             r"Algorithm & Year & Venue & Topology & Mechanism (axis B) & Objective (axis C) & Info.\\ (axis D) & DL & Simulator & Base. \\",
             r"\midrule"]
    for d in rows:
        base = {"det": "det.", "adapt": "adapt.", "n/v": "n/v"}[d["base"]]
        lines.append(f"{d['name']}~\\cite{{{d['key']}}} & {d['year']} & {esc(d['venue'])} & {esc(d['topo'])} & "
                     f"{esc(d['mech'])} & {esc(d['obj'])} & {esc(d['info'])} & {d['dl']} & {esc(d['sim'])} & {base} \\\\")
    lines += [r"\bottomrule", r"\end{tabularx}"]
    write("tab_corpus.tex", "\n".join(lines) + "\n")
    hist = sorted(HISTORICAL, key=lambda d: (d["year"], d["name"].lower()))
    lines = [r"\begin{tabularx}{\linewidth}{@{}>{\raggedright\arraybackslash}p{3.6cm}c>{\raggedright\arraybackslash}p{2.0cm}"
             r">{\raggedright\arraybackslash}X c@{}}", r"\toprule",
             r"Design & Year & Topology & Mechanism & Deadlock argument \\", r"\midrule"]
    for d in hist:
        lines.append(f"{d['name']}~\\cite{{{d['key']}}} & {d['year']} & {d['topo']} & {esc(d['mech'])} & {d['dl']} \\\\")
    lines += [r"\bottomrule", r"\end{tabularx}"]
    write("tab_hist.tex", "\n".join(lines) + "\n")


def audit() -> None:
    n = len(CORPUS)
    dl = Counter(d["dl"] for d in CORPUS)
    sim = Counter(d["sim"] for d in CORPUS)
    base = Counter(d["base"] for d in CORPUS)
    topo = Counter(d["topo"] for d in CORPUS)
    M.set("corpus.n", n)
    M.set("corpus.hist", len(HISTORICAL))
    M.set("corpus.dl.verif", n - dl["n/v"])
    M.set("corpus.formal", dl["F"])
    M.set("corpus.inherit", dl["R"])
    M.set("corpus.claim", dl["C"])
    M.set("corpus.mitig", dl["I"])
    M.set("corpus.nodl", dl["N"])
    M.set("corpus.dl.nv", dl["n/v"])
    M.set("corpus.weak", dl["C"] + dl["I"] + dl["N"])
    known_sim = {k: v for k, v in sim.items() if k != "n/v"}
    M.set("corpus.sim.verif", sum(known_sim.values()))
    M.set("corpus.simdistinct", len(known_sim))
    M.set("corpus.noxim", sum(v for k, v in known_sim.items() if k in NOXIM))
    M.set("corpus.base.verif", n - base["n/v"])
    M.set("corpus.detonly", base["det"])
    M.set("corpus.adapt", base["adapt"])
    M.set("corpus.vpc", topo["VPC"])
    M.set("corpus.25d", topo["2.5D"])
    M.set("corpus.fc", topo["3D-FC"] + topo["3D-FC + wireless"] + topo["2D/3D"] + topo["3D + wireless"])
    yrs = Counter(d["year"] for d in CORPUS)
    for y in range(2021, 2027):
        M.set(f"corpus.y{y}", yrs.get(y, 0))
    # Each work counts once per objective whose stem appears anywhere in its objective field
    # ("VL faults", "TSV faults" and "HL and VL faults" are fault objectives).
    stems = {"faults": "fault", "congestion": "congestion", "thermal": "thermal", "security": "security",
             "energy": "energy"}
    obj = Counter()
    for d in CORPUS:
        for k, stem in stems.items():
            if stem in d["obj"].lower():
                obj[k] += 1
    for k in stems:
        M.set(f"corpus.obj.{k}", obj.get(k, 0))
    # Simulator families: the Noxim variants count as one family.
    fam = {("Noxim" if k in NOXIM else k) for k in known_sim}
    M.set("corpus.simfamilies", len(fam))
    info = Counter(d["info"] for d in CORPUS)
    M.set("corpus.info.global", sum(v for k, v in info.items() if "global" in k))
    M.set("corpus.info.learned", info.get("learned", 0))
    hdl = Counter(d["dl"] for d in HISTORICAL)
    M.set("hist.formal", hdl["F"])


def build() -> None:
    corpus_tables()
    audit()
