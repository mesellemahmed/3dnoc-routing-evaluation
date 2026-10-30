"""Saturation refinement by bisection (added after the internal audit).

The coarse load sweep (groups 'sweep' and 'sweep_vpc') brackets the saturation
load of every (class, pattern) pair between the highest stable grid point
(mean undelivered fraction <= 1 % over the three seeds, no deadlocked seed) and
the next grid point. This script halves that bracket ROUNDS times. Every probe
uses the three seeds, the run length and the router model of the coarse sweep;
its runs are appended to results/sweep_refine.jsonl with group 'sweepref'
(3D-FC) or 'sweepref_vpc' (3D-VPC). Saturation throughput keeps its registered
definition: the highest accepted throughput among all stable sweep points,
coarse or refined.

Usage: python3 run_refine.py [--rounds 5] [--workers 2] [--coarse results/campaign3.jsonl] [--out results/sweep_refine.jsonl]
(paths relative to sim3/). The refinement is resumable: probes already present in --out are not re-run.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

SIM3 = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(SIM3))
sys.path.insert(0, str(SIM3.parent / "analysis3"))
from campaign3 import BASE, FC, SWEEP_PIR, SWEEP_VPC_PIR, VPC, job_key, run_one  # noqa: E402
from load import load  # noqa: E402

COARSE = SIM3 / "results" / "campaign3.jsonl"
OUT = SIM3 / "results" / "sweep_refine.jsonl"
PATTERNS = ["random", "transpose1", "bitreversal"]
SEEDS = (1, 2, 3)


def cfg(fam: str, alg: str, tr: str, pir: float, sd: int) -> dict:
    c = {**BASE, "group": "sweepref" if fam == "fc" else "sweepref_vpc", "routing": alg, "traffic": tr,
         "pir": round(pir, 7), "seed": sd, "trace": 0, "cycles": 12000, "warmup": 2000}
    if fam == "vpc":
        c["topo"] = "vpc"
    return c


def stable_table(frames) -> dict:
    """(fam, alg, tr, pir) -> stable?  over the three seeds."""
    import pandas as pd
    df = pd.concat(frames, ignore_index=True)
    out = {}
    for (grp, alg, tr, pir), g in df.groupby(["group", "routing", "traffic", "pir"]):
        fam = "vpc" if grp.endswith("vpc") else "fc"
        out[(fam, alg, tr, round(pir, 7))] = bool(g["c.undelivered"].mean() <= 0.01 and not g.deadlocked.any())
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rounds", type=int, default=5)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--coarse", type=Path, default=COARSE, help="coarse campaign results (jsonl)")
    ap.add_argument("--out", type=Path, default=OUT, help="refinement results (jsonl, appended)")
    ap.add_argument("--dry-run", action="store_true", help="print the probes of the next round and stop")
    args = ap.parse_args()
    out_path = args.out if args.out.is_absolute() else SIM3 / args.out
    coarse_path = args.coarse if args.coarse.is_absolute() else SIM3 / args.coarse
    coarse = load(coarse_path)
    coarse = coarse[coarse.group.isin(["sweep", "sweep_vpc"])]
    pairs = [("fc", a, t) for a in FC for t in PATTERNS] + [("vpc", a, t) for a in VPC for t in PATTERNS]
    grids = {"fc": SWEEP_PIR, "vpc": SWEEP_VPC_PIR}
    done = set()
    if out_path.exists():
        done = {json.loads(line)["key"] for line in out_path.open()}
    def bracket(st: dict, fam: str, alg: str, tr: str):
        pts = sorted(p for (f, a, t, p) in st if (f, a, t) == (fam, alg, tr))
        stab = [p for p in pts if st[(fam, alg, tr, p)]]
        if not stab:
            return None
        lo = max(stab)
        above = [p for p in pts if p > lo]
        return (lo, min(above)) if above else None   # None: stable at the highest grid point

    # Width of the coarse bracket of every pair; the refinement stops once a bracket has been halved
    # args.rounds times, so that re-running the script on complete results launches no new run.
    st0 = stable_table([coarse])
    width0 = {pr: (b[1] - b[0]) for pr in pairs if (b := bracket(st0, *pr))}
    for rnd in range(1, args.rounds + 1):
        frames = [coarse]
        if out_path.exists() and out_path.stat().st_size:
            frames.append(load(out_path))
        st = stable_table(frames)
        jobs, brackets = [], {}
        for fam, alg, tr in pairs:
            b = bracket(st, fam, alg, tr)
            if b is None or (fam, alg, tr) not in width0:
                continue
            lo, hi = b
            if hi - lo <= width0[(fam, alg, tr)] / 2 ** args.rounds + 1e-6:
                continue  # already refined to the requested resolution
            mid = round((lo + hi) / 2, 7)
            brackets[(fam, alg, tr)] = (lo, hi)
            for sd in SEEDS:
                c = cfg(fam, alg, tr, mid, sd)
                if job_key(c) not in done:
                    jobs.append(c)
        print(f"round {rnd}: {len(brackets)} brackets, {len(jobs)} runs {time.strftime('%H:%M:%S')}", flush=True)
        if args.dry_run:
            for j in sorted(jobs, key=lambda c: (c["group"], c["routing"], c["traffic"], c["seed"])):
                print(j["group"], j["routing"], j["traffic"], j["pir"], j["seed"])
            return
        if not jobs:
            continue
        with Pool(args.workers) as pool, out_path.open("a") as fh:
            for res in pool.imap_unordered(run_one, jobs, chunksize=1):
                fh.write(json.dumps(res, default=float) + "\n")
                fh.flush()
                done.add(res["key"])
    print("done", time.strftime('%H:%M:%S'), flush=True)


if __name__ == "__main__":
    main()
