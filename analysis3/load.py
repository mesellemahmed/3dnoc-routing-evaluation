"""Load campaign v3 results into a flat pandas DataFrame."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

RESULTS = Path(__file__).resolve().parent.parent / "sim3" / "results" / "campaign3.jsonl"


def flatten(d: dict, prefix: str = "") -> dict:
    out = {}
    for k, v in d.items():
        key = f"{prefix}{k}"
        if isinstance(v, dict):
            out.update(flatten(v, key + "."))
        elif isinstance(v, list):
            out[key] = v
        else:
            out[key] = v
    return out


def load(path: Path = RESULTS) -> pd.DataFrame:
    rows = []
    with open(path) as fh:
        for line in fh:
            r = json.loads(line)
            flat = flatten(r["config"], "")
            flat.update(flatten(r["stats"], "st."))
            flat.update(flatten(r["classical"], "c."))
            flat.update(flatten(r.get("activity", {}), "a."))
            if "era" in r:
                flat.update(flatten(r["era"], "e."))
            flat["t_total"] = r.get("t_total")
            rows.append(flat)
    df = pd.DataFrame(rows)
    df["deadlocked"] = (df["st.deadlock_cycle"] >= 0) | (~df["st.drained"].astype(bool))
    return df


if __name__ == "__main__":
    df = load()
    print(df.shape)
    print(df.groupby("group").size())
