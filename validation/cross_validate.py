"""Cross-validation of the in-house simulator against BookSim 2.0.

Both simulators model a k-ary 3-mesh with input-queued virtual-channel routers,
2 VCs of 8 flits, atomic VC reallocation (BookSim: wait_for_tail_credit = 1),
8-flit packets, Bernoulli injection and uniform random traffic. BookSim uses
its default iq_router pipeline (VC and switch allocation, switch and link
traversal); the in-house simulator is run with a matching per-hop delay of 4
cycles and a credit delay of 1 cycle. Four routing functions exist in both:
dimension order, minimal adaptive with a dimension-order escape VC, ROMM and
Valiant.

Outputs results/validation.json with, per (k, routing, load), the mean packet
latency and the accepted throughput of both simulators.
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import tempfile
from multiprocessing import Pool
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
BOOKSIM = Path(os.environ.get("BOOKSIM", HERE / "booksim"))  # BookSim 2 binary, commit 28f4329
SIM = HERE.parent / "sim3" / "noc3d_sim"
ROUTING = {"xyz": "dim_order", "fa": "min_adapt", "romm": "romm", "val": "valiant"}
BASE_CFG = """topology = mesh;
k = {k};
n = 3;
routing_function = {rf};
num_vcs = 2;
vc_buf_size = 8;
wait_for_tail_credit = 1;
vc_allocator = separable_input_first;
sw_allocator = separable_input_first;
alloc_iters = 1;
routing_delay = 0;
vc_alloc_delay = 1;
sw_alloc_delay = 1;
st_final_delay = 1;
credit_delay = 0;
traffic = uniform;
packet_size = 8;
injection_rate_uses_flits = 1;
injection_rate = {rate};
sim_type = latency;
warmup_periods = 3;
sample_period = 3000;
max_samples = 5;
latency_thres = 2000.0;
seed = {seed};
"""


def run_booksim(k: int, alg: str, rate: float, seed: int) -> dict:
    with tempfile.TemporaryDirectory() as tmp:
        cfg = Path(tmp) / "c.cfg"
        cfg.write_text(BASE_CFG.format(k=k, rf=ROUTING[alg], rate=rate, seed=seed))
        out = subprocess.run([str(BOOKSIM), str(cfg)], capture_output=True, text=True, timeout=3600).stdout
    unstable = "unstable" in out.lower()
    lat = re.findall(r"Packet latency average = ([0-9.]+) \(", out)
    acc = re.findall(r"Accepted flit rate average = ([0-9.]+) \(", out)
    hops = re.findall(r"Hops average = ([0-9.]+) \(", out)
    return {"latency": float(lat[-1]) if lat and not unstable else None,
            "accepted": float(acc[-1]) if acc and not unstable else None,
            "hops": float(hops[-1]) if hops and not unstable else None,
            "unstable": unstable}


def run_ours(k: int, alg: str, rate: float, seed: int) -> dict:
    with tempfile.TemporaryDirectory() as tmp:
        prefix = Path(tmp) / "r"
        cycles, warmup = 25000, 5000
        subprocess.run([str(SIM), "--kx", str(k), "--ky", str(k), "--kz", str(k), "--routing", alg,
                        "--traffic", "random", "--pir", str(rate / 8.0), "--seed", str(seed),
                        "--cycles", str(cycles), "--warmup", str(warmup), "--hop-delay", "4",
                        "--credit-delay", "1", "--trace", "0", "--out", str(prefix)], check=True)
        pk = np.fromfile(prefix.with_suffix(".pkt"), dtype=np.int32).reshape(-1, 11)
    n = k ** 3
    meas = pk[(pk[:, 2] >= warmup) & (pk[:, 2] < cycles - 2000)]
    dl = meas[meas[:, 4] >= 0]
    lat = (dl[:, 4] - dl[:, 2]).astype(float)
    acc = ((pk[:, 4] >= warmup) & (pk[:, 4] < cycles)).sum() * 8.0 / (cycles - warmup) / n
    undelivered = 1.0 - len(dl) / max(1, len(meas))
    saturated = undelivered > 0.01
    return {"latency": float(lat.mean()) if len(lat) and not saturated else None,
            "accepted": float(acc), "hops": float(dl[:, 5].mean()) if len(dl) else None,
            "unstable": saturated}


def job(args: tuple) -> dict:
    sim, k, alg, rate, seed = args
    res = run_booksim(k, alg, rate, seed) if sim == "booksim" else run_ours(k, alg, rate, seed)
    return {"sim": sim, "k": k, "routing": alg, "rate": rate, "seed": seed, **res}


def main() -> None:
    rates = {4: [0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60],
             8: [0.01, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35]}
    seeds = {4: (1, 2), 8: (1,)}
    jobs = [(sim, k, alg, r, seed) for k in (4, 8) for alg in ROUTING for r in rates[k]
            for seed in seeds[k] for sim in ("booksim", "ours")]
    jobs.sort(key=lambda j: -j[1])
    with Pool(2) as pool:
        out = pool.map(job, jobs, chunksize=1)
    (HERE / "results").mkdir(exist_ok=True)
    (HERE / "results" / "validation.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
