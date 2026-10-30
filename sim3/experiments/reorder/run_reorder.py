#!/usr/bin/env python3
"""Out-of-order delivery per (source, destination) flow: fraction of delivered packets that are ejected before an
earlier-created packet of the same flow. Permutation traffic (transpose, bit-reversal: one persistent flow per
source), 8x8x4 (FC) or 16-elevator VPC, PIR 0.006 and 0.012."""
import json
import subprocess
import tempfile
from pathlib import Path

import numpy as np

SIM = Path(__file__).resolve().parents[2] / "noc3d_sim"
FC = ["xyz", "zxy", "nf", "dyad", "ham", "o1t", "romm", "val", "fa", "rca", "qr", "ta"]
VPC = ["ef", "aesg", "adele", "fl"]


def reorder_fraction(pkt: np.ndarray, warmup: int) -> float:
    pk = pkt[(pkt[:, 2] >= warmup) & (pkt[:, 4] >= 0) & (pkt[:, 7] == 0)]
    order = np.lexsort((pk[:, 2], pk[:, 1], pk[:, 0]))
    pk = pk[order]
    ooo, n = 0, 0
    flow = pk[:, 0] * 100000 + pk[:, 1]
    start = 0
    for i in range(1, len(pk) + 1):
        if i == len(pk) or flow[i] != flow[start]:
            ej = pk[start:i, 4]
            run_max = np.maximum.accumulate(ej)
            ooo += int((ej[1:] < run_max[:-1]).sum())
            n += i - start
            start = i
    return ooo / n if n else float("nan")


def main() -> None:
    out = []
    for alg in FC + VPC:
        for tr in ("transpose1", "bitreversal"):
          for pir in (0.006, 0.012):
            for seed in (1, 2, 3):
                with tempfile.TemporaryDirectory() as tmp:
                    pre = Path(tmp) / "r"
                    cmd = [str(SIM), "--topo", "vpc" if alg in VPC else "fc", "--routing", alg, "--traffic", tr,
                           "--pir", str(pir), "--seed", str(seed), "--cycles", "18000", "--warmup", "3000",
                           "--trace", "0", "--out", str(pre)]
                    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                    pkt = np.fromfile(pre.with_suffix(".pkt"), dtype="<i4").reshape(-1, 11)
                    out.append(dict(routing=alg, traffic=tr, pir=pir, seed=seed, ooo=reorder_fraction(pkt, 3000)))
                    print(alg, tr, pir, seed, round(out[-1]["ooo"], 4), flush=True)
    Path(__file__).with_name("reorder.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
