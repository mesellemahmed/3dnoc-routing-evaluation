"""Statistical helpers: Student-t intervals, paired Wilcoxon tests with Holm
correction, Cliff's delta with bootstrap intervals, Wilson intervals."""
from __future__ import annotations

import math

import numpy as np
from scipy import stats as st


def mean_ci(x, conf: float = 0.95) -> tuple[float, float, float]:
    """Mean and half-width of the Student-t confidence interval."""
    x = np.asarray([v for v in x if v is not None and np.isfinite(v)], dtype=float)
    if len(x) == 0:
        return float("nan"), float("nan"), 0
    m = float(x.mean())
    if len(x) < 2:
        return m, float("nan"), len(x)
    h = float(st.t.ppf(0.5 + conf / 2, len(x) - 1) * x.std(ddof=1) / math.sqrt(len(x)))
    return m, h, len(x)


def cliffs_delta(a, b) -> float:
    """Cliff's delta: P(a > b) - P(a < b) over all cross pairs."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) == 0 or len(b) == 0:
        return float("nan")
    gt = (a[:, None] > b[None, :]).sum()
    lt = (a[:, None] < b[None, :]).sum()
    return float((gt - lt) / (len(a) * len(b)))


def paired_delta(a, b) -> float:
    """Matched-pairs effect size: fraction of pairs with a > b minus fraction
    with a < b (the paired analogue of Cliff's delta)."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    d = a - b
    return float(((d > 0).sum() - (d < 0).sum()) / len(d)) if len(d) else float("nan")


def bootstrap_ci(fn, a, b, n_boot: int = 2000, seed: int = 20261025, conf: float = 0.95):
    """Percentile bootstrap interval of a paired statistic fn(a, b)."""
    rng = np.random.default_rng(seed)
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    n = len(a)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        vals.append(fn(a[idx], b[idx]))
    lo, hi = np.quantile(vals, [(1 - conf) / 2, 1 - (1 - conf) / 2])
    return float(lo), float(hi)


def wilcoxon_p(a, b) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    d = a - b
    if np.allclose(d, 0):
        return 1.0
    return float(st.wilcoxon(a, b, zero_method="wilcox", alternative="two-sided").pvalue)


def holm(pvals: dict) -> dict:
    """Holm step-down adjusted p-values for a dict name -> p."""
    items = sorted(pvals.items(), key=lambda kv: kv[1])
    m = len(items)
    out, running = {}, 0.0
    for i, (k, p) in enumerate(items):
        adj = min(1.0, (m - i) * p)
        running = max(running, adj)
        out[k] = running
    return out


def wilson(k: int, n: int, conf: float = 0.95) -> tuple[float, float]:
    if n == 0:
        return float("nan"), float("nan")
    z = st.norm.ppf(0.5 + conf / 2)
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return c - h, c + h


def jain(u) -> float:
    u = np.asarray(u, dtype=float)
    s2 = (u * u).sum()
    return float(u.sum() ** 2 / (len(u) * s2)) if s2 > 0 else float("nan")
