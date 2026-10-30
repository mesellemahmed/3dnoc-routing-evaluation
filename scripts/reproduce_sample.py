#!/usr/bin/env python3
"""Re-execute a random sample of campaign runs and compare the simulator statistics bit for bit.

Every run of sim3/results/*.jsonl stores its full configuration and the statistics written by the simulator.
This script rebuilds the simulator command line of N randomly chosen runs, re-runs them, and checks that the
new statistics are identical to the stored ones (traffic and routing use separate, seeded random streams).
Usage (from the repository root, after building sim3/noc3d_sim):
    python3 scripts/reproduce_sample.py [N] [--seed S] [--file sim3/results/campaign3.jsonl]
"""
import argparse
import json
import random
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SIM = ROOT / "sim3" / "noc3d_sim"


def command(c: dict, out: str) -> list:
    a = {"--topo": c["topo"], "--elevators": c["elevators"], "--routing": c["routing"], "--traffic": c["traffic"],
         "--pir": c["pir"], "--seed": c["seed"], "--window": c["window"], "--kx": c["kx"], "--ky": c["ky"],
         "--kz": c["kz"], "--flits": c["flits"], "--vcs": c["vcs"], "--depth": c["depth"], "--atomic": c["atomic"],
         "--hop-delay": c["hop_delay"], "--credit-delay": c["credit_delay"], "--cycles": c["cycles"],
         "--warmup": c["warmup"], "--trace": 0, "--pv": c["pv"], "--xi": c["xi"],
         "--attackers": c["attackers"], "--attack-pir": c["attack_pir"], "--out": out}
    return [str(SIM)] + [str(x) for kv in a.items() for x in kv]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("n", nargs="?", type=int, default=20)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--file", default=str(ROOT / "sim3" / "results" / "campaign3.jsonl"))
    args = ap.parse_args()
    lines = Path(args.file).read_text().splitlines()
    picks = random.Random(args.seed).sample(range(len(lines)), args.n)
    bad = 0
    for i in picks:
        r = json.loads(lines[i])
        c = r["config"]
        with tempfile.TemporaryDirectory() as t:
            subprocess.run(command(c, f"{t}/o"), check=True, capture_output=True)
            st = json.loads(Path(f"{t}/o.stats").read_text())
        same = st == r["stats"]
        bad += not same
        print(f"line {i:5d}  {c['group']:12s} {c['routing']:6s} {c['traffic']:11s} PIR {c['pir']:<8g} "
              f"{'identical' if same else 'DIFFERENT'}", flush=True)
    print(f"{args.n - bad}/{args.n} runs reproduced bit for bit")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
