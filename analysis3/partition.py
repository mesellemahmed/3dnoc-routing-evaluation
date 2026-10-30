"""Partition risk and loss variance of VPC stacks under the Gaussian-copula
TSV fault model (exponential covariance exp(-d/xi) over the planar distance
between columns, tier boundaries independent, exact marginal p_v).

For each elevator placement (stride2: 16 columns, stride4: 4 columns, full: 64
columns of an 8x8 tier) and each (p_v, xi):
  p_cut     P(no healthy elevator across one boundary)   -- Slepian-monotone
  p_cut3    P(at least one of the 3 boundaries severed)
  sd_frac   standard deviation of the failed fraction of elevators per boundary
  mean_frac mean failed fraction (= p_v, invariance check)
Monte Carlo with 2,000,000 boundary samples per configuration plus the exact
independent value for xi = 0. Output: partition.json.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.stats import norm

KX = KY = 8
PLACEMENTS = {
    "stride2": [(x, y) for y in range(KY) for x in range(KX) if x % 2 == 1 and y % 2 == 1],
    "stride4": [(x, y) for y in range(KY) for x in range(KX) if x % 4 == 1 and y % 4 == 1],
    "full": [(x, y) for y in range(KY) for x in range(KX)],
}
PV = (0.1, 0.2, 0.3)
XI = (0.0, 1.0, 2.0, 4.0, 8.0)
NS = 2_000_000


def run() -> dict:
    rng = np.random.default_rng(20261025)
    out = {}
    for name, cols in PLACEMENTS.items():
        c = np.array(cols, dtype=float)
        d = np.sqrt(((c[:, None, :] - c[None, :, :]) ** 2).sum(-1))
        for xi in XI:
            S = np.exp(-d / xi) if xi > 0 else np.eye(len(cols))
            L = np.linalg.cholesky(S + 1e-12 * np.eye(len(cols)))
            fails_all = {pv: 0 for pv in PV}
            s1 = {pv: 0.0 for pv in PV}
            s2 = {pv: 0.0 for pv in PV}
            done = 0
            while done < NS:
                b = min(250_000, NS - done)
                z = rng.standard_normal((b, len(cols))) @ L.T
                for pv in PV:
                    f = z < norm.ppf(pv)
                    k = f.sum(1)
                    fails_all[pv] += int((k == len(cols)).sum())
                    frac = k / len(cols)
                    s1[pv] += float(frac.sum())
                    s2[pv] += float((frac ** 2).sum())
                done += b
            for pv in PV:
                p_cut = fails_all[pv] / NS
                if xi == 0:
                    p_cut = pv ** len(cols)          # exact
                mean = s1[pv] / NS
                sd = float(np.sqrt(max(0.0, s2[pv] / NS - mean ** 2)))
                out[f"{name}|{pv:g}|{xi:g}"] = dict(n_cols=len(cols), p_cut=p_cut,
                                                     p_cut3=1 - (1 - p_cut) ** 3, sd_frac=sd, mean_frac=mean)
    return out


if __name__ == "__main__":
    res = run()
    Path(__file__).with_name("partition.json").write_text(json.dumps(res, indent=1))
    for k, v in res.items():
        print(k, f"cut={v['p_cut']:.2e} cut3={v['p_cut3']:.2e} sd={v['sd_frac']:.3f} mean={v['mean_frac']:.3f}")
