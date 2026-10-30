#!/usr/bin/env python3
"""Rebuild summary.json of the exact VPC fault-path check.

Inputs: replica_20000maps.txt and replica_5maps.txt (output of `./replica 20000` and `./replica 5`; columns
p_v, xi, seed, realised fraction of failed elevator segments, exact expected Elevator-First loss under uniform
traffic) and the vpcfault runs of results/campaign3.jsonl (EF, uniform traffic, seeds 1-5, whose fault maps are
the maps of seeds 1-5 of the replica).

Usage: g++ -O2 -std=c++17 -o replica replica.cpp && ./replica 20000 > replica_20000maps.txt &&
       ./replica 5 > replica_5maps.txt && python3 summarize.py
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "analysis3"))
from load import load  # noqa: E402


def read(name: str) -> dict:
    out = {}
    for line in (HERE / name).read_text().split("\n"):
        if line.strip():
            pv, xi, seed, seg, loss = line.split()
            out.setdefault(f"{float(pv):.1f}|{float(xi):.1f}", []).append((int(seed), float(seg), float(loss)))
    return out


def main() -> None:
    big, five = read("replica_20000maps.txt"), read("replica_5maps.txt")
    res = {}
    for key, rows in big.items():
        loss = np.array([r[2] for r in rows])
        res[key] = {"maps": len(rows), "loss_mean": float(loss.mean()), "loss_sd": float(loss.std(ddof=1)),
                    "seg_fail_mean": float(np.mean([r[1] for r in rows]))}
    df = load(HERE.parents[1] / "results" / "campaign3.jsonl")
    sim = df[(df.group == "vpcfault") & (df.routing == "ef") & (df.traffic == "random")]
    res["five_maps_mean"] = {k: float(np.mean([r[2] for r in v])) for k, v in five.items()}
    run_diff, conf_diff = [], []
    for key, rows in five.items():
        pv, xi = (float(x) for x in key.split("|"))
        exact = {r[0]: r[2] for r in rows}
        s = sim[(sim.pv == pv) & (sim.xi == xi)]
        run_diff += [abs(r["c.drop_ratio"] - exact[int(r["seed"])]) for _, r in s.iterrows()]
        conf_diff.append(abs(s["c.drop_ratio"].mean() - np.mean(list(exact.values()))))
    res["per_run_max_absdiff"] = float(max(run_diff))
    res["per_config_max_absdiff"] = float(max(conf_diff))
    print(f"{len(run_diff)} simulated runs compared with their exact expected loss", file=sys.stderr)
    print(json.dumps(res, indent=1))
    if "--write" in sys.argv:
        (HERE / "summary.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
