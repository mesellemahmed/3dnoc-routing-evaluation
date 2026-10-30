#!/usr/bin/env python3
"""Assemble the standalone LaTeX source of the v3 article.

Inlines the generated macros, plot data and tables (gen/*.tex), concatenates the
hand-written parts p01..p07, builds the bibliography from the keys actually
cited (order of first citation) and reports missing keys or undefined macros.
Output: routing_3dnoc_survey_v3.tex (self-contained, compiles with pdflatex).
"""
from __future__ import annotations

import json
import re
import sys

from common import GEN, HERE
from refs import REFS

PARTS = ["p01_front.tex", "p02_intro.tex", "p03_background.tex", "p04_survey.tex", "p05_protocol.tex",
         "p06_results.tex", "p07_verdict.tex"]
OUT = HERE / "routing_3dnoc_survey_v3.tex"


def inline_inputs(text: str) -> str:
    def rep(m: re.Match) -> str:
        path = HERE / m.group(1)
        if not path.exists():
            sys.exit(f"missing generated file {path}")
        return path.read_text().rstrip("\n")
    return re.sub(r"\\input\{(gen/[A-Za-z0-9_\-]+\.tex)\}", rep, text)


def main() -> None:
    body = "\n".join((HERE / p).read_text() for p in PARTS)
    body = body.replace("%%GEN:macros%%", (GEN / "macros.tex").read_text())
    body = body.replace("%%GEN:plotdata%%", (GEN / "curves.tex").read_text())
    body = body.replace("%%GEN:plotdata2%%", "")
    body = inline_inputs(body)
    # Citations in order of first appearance.
    keys: list[str] = []
    for grp in re.findall(r"\\cite\{([^}]*)\}", body):
        for k in grp.split(","):
            k = k.strip()
            if k and k not in keys:
                keys.append(k)
    missing = [k for k in keys if k not in REFS]
    if missing:
        sys.exit(f"missing bibliography entries: {missing}")
    bib = ["", r"\begin{thebibliography}{" + "9" * len(str(len(keys))) + "}", r"\small",
           r"\setlength{\itemsep}{1pt}"]
    for k in keys:
        bib.append(r"\bibitem{" + k + "} " + REFS[k])
    bib.append(r"\end{thebibliography}")
    body = body.rstrip() + "\n" + "\n".join(bib) + "\n\n\\end{document}\n"
    # Macro check: every \val{...} must be defined.
    values = json.loads((HERE / "generated_values.json").read_text())
    used = set(re.findall(r"\\val\{([^}]*)\}", body))
    undefined = sorted(u for u in used if u not in values)
    if undefined:
        sys.exit(f"undefined macros: {undefined}")
    OUT.write_text(body)
    unused = len(REFS) - len(keys)
    print(f"{OUT.name}: {len(body.splitlines())} lines, {len(keys)} references cited "
          f"({unused} entries of refs.py unused), {len(used)} distinct macros used")


if __name__ == "__main__":
    main()
