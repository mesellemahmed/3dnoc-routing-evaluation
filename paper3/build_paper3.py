#!/usr/bin/env python3
"""Generate every number, table and figure of the v3 article from the raw
results, then assemble the standalone LaTeX source.

Inputs: sim3/results (campaign v3, CDG checker), bench3/fault_reach.csv,
validation/results/validation.json, hw/rpt_*.txt, analysis3 (corpus, refs,
path counts). Output: gen/*.tex, generated_values.json and
routing_3dnoc_survey_v3.tex.
"""
from __future__ import annotations

import json

import pandas as pd

from common import GEN, HERE, M, ROOT
import sec_platform
import sec_results
import sec_survey


def main() -> None:
    df = sec_platform.load_all()
    val = json.loads((ROOT / "validation" / "results" / "validation.json").read_text())
    fr = pd.read_csv(ROOT / "bench3" / "fault_reach.csv")
    sec_platform.section_counts(df, val, fr)
    sec_survey.build()
    sec_platform.section_validation(val)
    sec_platform.section_nonatomic()
    cdg = sec_platform.section_cdg()
    sec_results.prepare_dl(df)
    sec_results.d1(df, cdg)
    sec_results.d2()
    sec_results.d3(df)
    sec_results.d4(df)
    sec_results.d5(df, fr)
    sec_results.d6(df)
    sec_results.d7(df)
    sec_results.d8()
    sec_results.d9_d10(df)
    sec_results.d11(df)
    sec_results.d12(df)
    sec_results.paired(df)
    sec_results.partition()
    sec_results.vpc4(df)
    sec_results.faultcheck()
    sec_results.reorder()
    sec_results.verdict()
    sec_results.threshold_sensitivity()
    (GEN / "macros.tex").write_text(M.tex())
    (HERE / "generated_values.json").write_text(json.dumps(M.v, indent=1, sort_keys=True))
    print(f"{len(M.v)} macros written")


if __name__ == "__main__":
    main()
