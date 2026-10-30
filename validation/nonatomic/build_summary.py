#!/usr/bin/env python3
"""Rebuild summary.json of the non-atomic reallocation experiments from their raw outputs.

Sections and sources:
  booksim_batch_nonatomic     batch_results.txt, batch_results2.txt (run_batch.sh, run_batch2.sh; BookSim 2)
  inject_policy               ../../sim3/experiments/inject_policy/res.txt        (run.sh)
  sizes_default_injection     ../../sim3/experiments/inject_policy/res_sizes.txt  (run_sizes.sh)
  pipeline_default_injection  ../../sim3/experiments/pipeline/res.txt, res2.txt   (run.sh, run2.sh)
  equal_priority_selection    ../../sim3/experiments/rand_selection/res.txt       (run.sh)
  commitment_rates            ../../sim3/experiments/commitment/e_*.txt, instr*.err (run.sh)
Each in-house line reads: <config> <seed> <deadlock_cycle> <drained> <stuck flits>.
Usage: python3 build_summary.py [--write]
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXP = HERE.parents[1] / "sim3" / "experiments"


def lines(path: Path):
    return [l.split() for l in path.read_text().splitlines() if l.strip()]


def booksim() -> list:
    runs = []
    for name in ("batch_results.txt", "batch_results2.txt"):
        for line in (HERE / name).read_text().splitlines():
            if not line.startswith("k="):
                continue
            d = dict(kv.split("=", 1) for kv in line.split() if "=" in kv)
            runs.append(dict(k=int(d["k"]), traffic=d["tr"], wtc=int(d["wtc"]), batch=int(d["batch"]),
                             seed=int(d["seed"]), warn=int(d["warn"]), completed=d.get("dur", "") != "",
                             cycles=int(d["dur"]) if d.get("dur") else None))
    return runs


def counters(path: Path) -> dict:
    t = path.read_text().split()
    g, ne, fu, fa = (int(t[i]) for i in (1, 3, 5, 7))
    return {"grants": g, "nonempty": ne, "full": fu, "full_adaptive": fa,
            "pct_full_adaptive": 100 * fa / g, "pct_full": 100 * fu / g}


def main() -> None:
    s = {"booksim_batch_nonatomic": booksim()}
    s["inject_policy"] = [dict(size=a, inject_escape=int(b.split("=")[1]), traffic=c, seed=int(d),
                               drained=f == "True", stuck=int(g))
                          for a, b, c, d, e, f, g in lines(EXP / "inject_policy" / "res.txt")]
    s["sizes_default_injection"] = [dict(size=a, traffic=b, seed=int(c), drained=e == "True", stuck=int(f))
                                    for a, b, c, d, e, f in lines(EXP / "inject_policy" / "res_sizes.txt")]
    pipe = [dict(hop_delay=int(a.split("=")[1]), traffic=b, seed=int(c), drained=e == "True", stuck=int(f),
                 size="8x8x4") for a, b, c, d, e, f in lines(EXP / "pipeline" / "res.txt")]
    pipe += [dict(hop_delay=int(a.split("=")[1]), size=b, traffic=c, seed=int(d), drained=f == "True", stuck=int(g))
             for a, b, c, d, e, f, g in lines(EXP / "pipeline" / "res2.txt")]
    s["pipeline_default_injection"] = pipe
    s["equal_priority_selection"] = [dict(size=a, hop_delay=int(b.split("=")[1]), traffic=c, seed=int(d),
                                          drained=f == "True", stuck=int(g))
                                     for a, b, c, d, e, f, g in lines(EXP / "rand_selection" / "res.txt")]
    com = EXP / "commitment"
    s["commitment_rates"] = {"booksim_444": counters(com / "instr4.err"), "booksim_888": counters(com / "instr8.err")}
    for k in ("444", "888"):
        for h in ("1", "4"):
            s["commitment_rates"][f"ours_{k}_h{h}"] = counters(com / f"e_{k}_h{h}.txt")
    out = json.dumps(s, indent=1)
    if "--write" in sys.argv:
        (HERE / "summary.json").write_text(out)
    print(out)


if __name__ == "__main__":
    main()
