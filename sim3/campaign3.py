"""Evaluation campaign v3 for the survey (simulator v3).

Every run is identified by a key built from all of its parameters, so the
campaign is resumable and each group can be launched separately.

Groups
  era       ERA-traced grid on 3D-FC: 12 classes x 6 patterns x 5 PIR x 5 seeds
  era_vpc   ERA-traced grid on 3D-VPC (stride-2 elevators): 4 classes
  sweep     latency/throughput versus load up to saturation (untraced)
  sweep_vpc the same for 3D-VPC classes
  stress    overload with deadlock watchdog and drain test
  sens      router micro-architecture sensitivity (VCs, atomicity, buffer
            depth, packet length, pipeline depth) at low load and overload
  scal      network-size scalability (4x4x4, 8x8x4, 8x8x8, 16x16x4)
  dos       flooding denial-of-service: 8 attackers flood one router
  burst     self-similar (Pareto ON/OFF) traffic, ERA-traced
  vpcfault  correlated TSV faults on 3D-VPC, ERA-traced
  vpc4      sparse (stride-4) elevator layout
  window    sensitivity of ERA coupling/drift to the analysis window

The run key (resume, and seed of the ERA permutation test) is built from KEY_FIELDS. It keeps a constant
field 'bf_period' (always 100) inherited from the original driver: it is not passed to the simulator and
has no effect on any run, but removing it would change the keys, and with them the permutation seeds, of
the released results.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import shutil
import subprocess
import tempfile
import time
from multiprocessing import Pool
from pathlib import Path

from era_analyzer import analyse

SIM = Path(__file__).resolve().parent / "noc3d_sim"
FC = ["xyz", "zxy", "nf", "dyad", "ham", "o1t", "romm", "val", "fa", "rca", "qr", "ta"]
VPC = ["ef", "aesg", "adele", "fl"]
TRAFFIC = ["random", "transpose1", "shuffle", "bitreversal", "butterfly", "hotspot"]
PIR = [0.001, 0.003, 0.006, 0.009, 0.012]
SWEEP_PIR = [0.002, 0.006, 0.010, 0.014, 0.018, 0.022, 0.026, 0.030, 0.035, 0.040, 0.050]
SWEEP_VPC_PIR = [0.002, 0.004, 0.006, 0.008, 0.010, 0.012, 0.015, 0.020]

BASE = dict(topo="fc", elevators="stride2", pv=0.0, xi=0.0, window=1000, kx=8, ky=8, kz=4, flits=8,
            vcs=2, depth=8, atomic=1, hop_delay=1, credit_delay=0, cycles=18000, warmup=3000,
            trace=1, bf_period=100, attackers=0, attack_pir=0.0, traffic="random")
KEY_FIELDS = ("group", "topo", "elevators", "routing", "traffic", "pir", "seed", "kx", "ky", "kz", "flits",
              "vcs", "depth", "atomic", "hop_delay", "credit_delay", "pv", "xi", "attackers", "attack_pir",
              "bf_period", "window", "cycles")

SENS_VARIANTS = {
    "base": {},
    "v4": {"vcs": 4},
    "na": {"atomic": 0},
    "d4": {"depth": 4},
    "d16": {"depth": 16},
    "p4": {"flits": 4},
    "p16": {"flits": 16},
    "h3": {"hop_delay": 3, "credit_delay": 1},
}


def configs() -> list[dict]:
    jobs: list[dict] = []

    def add(group: str, **kw) -> None:
        jobs.append({**BASE, "group": group, **kw})

    for alg, tr, pir, sd in itertools.product(FC, TRAFFIC, PIR, range(1, 6)):
        add("era", routing=alg, traffic=tr, pir=pir, seed=sd)
    for alg, tr, pir, sd in itertools.product(VPC, TRAFFIC, PIR, range(1, 6)):
        add("era_vpc", topo="vpc", routing=alg, traffic=tr, pir=pir, seed=sd)
    for alg, tr, pir, sd in itertools.product(FC, ["random", "transpose1", "bitreversal"], SWEEP_PIR, range(1, 4)):
        add("sweep", routing=alg, traffic=tr, pir=pir, seed=sd, trace=0, cycles=12000, warmup=2000)
    for alg, tr, pir, sd in itertools.product(VPC, ["random", "transpose1", "bitreversal"], SWEEP_VPC_PIR,
                                              range(1, 4)):
        add("sweep_vpc", topo="vpc", routing=alg, traffic=tr, pir=pir, seed=sd, trace=0, cycles=12000,
            warmup=2000)
    for alg, tr, pir, sd in itertools.product(FC, ["random", "transpose1", "hotspot"], [0.03, 0.06, 0.10],
                                              range(1, 6)):
        add("stress", routing=alg, traffic=tr, pir=pir, seed=sd, cycles=15000, trace=0)
    for alg, tr, pir, sd in itertools.product(VPC, ["random", "transpose1", "hotspot"], [0.03, 0.06, 0.10],
                                              range(1, 6)):
        add("stress", topo="vpc", routing=alg, traffic=tr, pir=pir, seed=sd, cycles=15000, trace=0)
    for (var, over), alg, tr, load, sd in itertools.product(SENS_VARIANTS.items(), FC + VPC,
                                                           ["random", "transpose1"], ["low", "over"],
                                                           range(1, 3)):
        cfg = {**BASE, **over}
        offered = 0.016 if load == "low" else 0.6
        add("sens", routing=alg, topo="vpc" if alg in VPC else "fc", traffic=tr,
            pir=round(offered / cfg["flits"], 6), seed=sd, trace=0, cycles=12000, warmup=2000,
            variant=var, load=load, **over)
    for (kx, ky, kz), alg, tr, load, sd in itertools.product([(4, 4, 4), (8, 8, 4), (8, 8, 8), (16, 16, 4)], FC,
                                                             ["random", "transpose1", "bitreversal"],
                                                             ["low", "over"], range(1, 3)):
        add("scal", routing=alg, traffic=tr, pir=0.002 if load == "low" else 0.075, seed=sd, trace=0,
            cycles=12000, warmup=2000, kx=kx, ky=ky, kz=kz, load=load)
    for alg, att, pir, sd in itertools.product(FC + VPC, [0, 8], [0.003, 0.006], range(1, 6)):
        add("dos", routing=alg, topo="vpc" if alg in VPC else "fc", traffic="dos", attackers=att,
            attack_pir=0.05 if att else 0.0, pir=pir, seed=sd, trace=0)
    for alg, pir, sd in itertools.product(FC + VPC, [0.003, 0.006, 0.009, 0.012], range(1, 6)):
        add("burst", routing=alg, topo="vpc" if alg in VPC else "fc", traffic="selfsim", pir=pir, seed=sd)
    for alg, pv, xi, tr, sd in itertools.product(["ef", "efr", "aesg", "adele", "fl"], [0.1, 0.2], [0.0, 4.0],
                                                 ["random", "transpose1"], range(1, 6)):
        add("vpcfault", topo="vpc", routing=alg, pv=pv, xi=xi, traffic=tr, pir=0.006, seed=sd)
    for alg, tr, pir, sd in itertools.product(VPC, ["random", "transpose1"], [0.001, 0.003, 0.006], range(1, 4)):
        add("vpc4", topo="vpc", elevators="stride4", routing=alg, traffic=tr, pir=pir, seed=sd, trace=0)
    for alg, win, tr, sd in itertools.product(["xyz", "fa"], [250, 500, 2000], ["random", "hotspot"],
                                              range(1, 4)):
        add("window", routing=alg, window=win, traffic=tr, pir=0.012, seed=sd)
    return jobs


def job_key(c: dict) -> str:
    return "|".join(f"{k}={c[k]}" for k in KEY_FIELDS)


def run_one(cfg: dict) -> dict:
    key = job_key(cfg)
    seed = int(hashlib.sha256(key.encode()).hexdigest()[:8], 16)
    work = Path(tempfile.mkdtemp(prefix="c3_"))
    prefix = work / "run"
    args = {"--topo": cfg["topo"], "--elevators": cfg["elevators"], "--routing": cfg["routing"],
            "--traffic": cfg["traffic"], "--pir": cfg["pir"], "--seed": cfg["seed"], "--window": cfg["window"],
            "--kx": cfg["kx"], "--ky": cfg["ky"], "--kz": cfg["kz"], "--flits": cfg["flits"], "--vcs": cfg["vcs"],
            "--depth": cfg["depth"], "--atomic": cfg["atomic"], "--hop-delay": cfg["hop_delay"],
            "--credit-delay": cfg["credit_delay"], "--cycles": cfg["cycles"], "--warmup": cfg["warmup"],
            "--trace": cfg["trace"], "--pv": cfg["pv"], "--xi": cfg["xi"],
            "--attackers": cfg["attackers"], "--attack-pir": cfg["attack_pir"], "--out": prefix}
    cmd = [str(SIM)] + [str(x) for kv in args.items() for x in kv]
    t0 = time.time()
    subprocess.run(cmd, check=True)
    t_sim = time.time() - t0
    stats = json.loads(prefix.with_suffix(".stats").read_text())
    keep_vec = cfg["kx"] == 8 and cfg["ky"] == 8 and cfg["kz"] == 4 and cfg["group"] in (
        "era", "era_vpc", "burst", "dos", "vpcfault")
    res = analyse(prefix, warmup=cfg["warmup"], cycles=cfg["cycles"], n_perm=200, seed=seed,
                  traced=bool(cfg["trace"]), keep_vector=keep_vec,
                  tail_guard=max(1000, (cfg["cycles"] - cfg["warmup"]) // 15))
    shutil.rmtree(work)
    return {"config": cfg, "key": key, "stats": stats, "t_sim": round(t_sim, 2),
            "t_total": round(time.time() - t0, 2), **res}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=Path("results/campaign3.jsonl"))
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--groups", default="era,era_vpc,sweep,sweep_vpc,stress,sens,scal,dos,burst,vpcfault,vpc4,"
                                        "window")
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()
    groups = args.groups.split(",")
    all_jobs = configs()
    keys = [job_key(c) for c in all_jobs]
    assert len(keys) == len(set(keys)), "duplicate job keys"
    if args.list:
        from collections import Counter
        print(Counter(c["group"] for c in all_jobs), len(all_jobs))
        return
    done = set()
    if args.out.exists():
        done = {json.loads(line)["key"] for line in args.out.open()}
    order = {g: i for i, g in enumerate(groups)}
    todo = [c for c in all_jobs if c["group"] in order and job_key(c) not in done]
    todo.sort(key=lambda c: (order[c["group"]], -c["kx"] * c["ky"] * c["kz"] * c["cycles"]))
    print(f"{len(todo)} runs to do", flush=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with Pool(args.workers) as pool, args.out.open("a") as fh:
        for i, res in enumerate(pool.imap_unordered(run_one, todo, chunksize=1), 1):
            fh.write(json.dumps(res, default=float) + "\n")
            fh.flush()
            if i % 100 == 0:
                print(f"{i}/{len(todo)} {time.strftime('%H:%M:%S')}", flush=True)


if __name__ == "__main__":
    main()
