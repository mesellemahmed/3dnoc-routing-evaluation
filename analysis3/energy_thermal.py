"""Parametric energy model and compact steady-state thermal model.

Energy (normalised to the energy of one router crossbar traversal, e_r = 1):
    E = sum_r xbar_r * e_r * w + planar * e_l * w + vertical * e_v * w
        + aux + lambda * slots * N * cycles
per delivered flit, where w = (64 + extra bits) / 64 accounts for fields that
a class would piggyback on every flit (w = 1 for every evaluated class), aux
counts table updates and side-band messages, slots = 7 * V * depth buffer
slots per router.
The 27 parameter combinations are e_l in {0.25, 0.5, 1}, e_v in {0.05, 0.2,
0.5} and lambda in {0, 5e-4, 2e-3}.

Thermal: G T = P on the tile grid (lateral g_lat, inter-tier g_vert, sink g_sink
under tier 0); P = p_core + p_leak + e_flit * f * xbar_rate. Four variants.
"""
from __future__ import annotations

import itertools

import numpy as np

E_L = (0.25, 0.5, 1.0)
E_V = (0.05, 0.2, 0.5)
LAMBDA = (0.0, 5e-4, 2e-3)
E_AUX = 0.1          # one table write or side-band message, relative to e_r
RCA_TOGGLE = 0.1     # activity factor of the 7-bit RCA side-band wires


def energy_per_flit(row: dict, e_l: float, e_v: float, lam: float) -> float:
    """row: one run (activity counts in the measurement window, stats)."""
    routing = row["routing"]
    span = row["cycles"] - row["warmup"]
    n = row["kx"] * row["ky"] * row["kz"]
    width = 1.0
    dyn = width * (row["a.xbar"] + e_l * row["a.planar"] + e_v * row["a.vertical"])
    aux = 0.0
    if routing == "qr":
        aux += E_AUX * row["st.q_updates"] + e_l * (16 / 64) * row["st.q_updates"]
    if routing == "rca":
        n_links = 2 * ((row["kx"] - 1) * row["ky"] * row["kz"] + row["kx"] * (row["ky"] - 1) * row["kz"]
                       + row["kx"] * row["ky"] * (row["kz"] - 1))
        aux += RCA_TOGGLE * e_l * (7 / 64) * n_links * span
    slots = 7 * row["vcs"] * row["depth"]
    leak = lam * slots * n * span
    delivered = row["c.throughput_flits_per_cycle"] * span
    return (dyn + aux + leak) / delivered if delivered > 0 else float("nan")


def energy_grid():
    return list(itertools.product(E_L, E_V, LAMBDA))


# ---------------------------------------------------------------------------
THERMAL_VARIANTS = {
    "base": dict(g_lat=6.0e-3, g_vert=0.14, g_sink=0.0625, p_core=0.25, p_leak=0.01, e_flit=30e-12),
    "noc_heavy": dict(g_lat=6.0e-3, g_vert=0.14, g_sink=0.0625, p_core=0.10, p_leak=0.01, e_flit=60e-12),
    "weak_vertical": dict(g_lat=6.0e-3, g_vert=0.14 / 3, g_sink=0.0625, p_core=0.25, p_leak=0.01, e_flit=30e-12),
    "strong_lateral": dict(g_lat=1.8e-2, g_vert=0.14, g_sink=0.0625, p_core=0.25, p_leak=0.01, e_flit=30e-12),
}
FREQ = 1.0e9


def conductance(kx: int, ky: int, kz: int, g_lat: float, g_vert: float, g_sink: float) -> np.ndarray:
    n = kx * ky * kz
    G = np.zeros((n, n))

    def rid(x, y, z):
        return x + kx * (y + ky * z)

    for z in range(kz):
        for y in range(ky):
            for x in range(kx):
                r = rid(x, y, z)
                for (dx, dy, dz, g) in ((1, 0, 0, g_lat), (0, 1, 0, g_lat), (0, 0, 1, g_vert)):
                    nx, ny, nz = x + dx, y + dy, z + dz
                    if nx < kx and ny < ky and nz < kz:
                        q = rid(nx, ny, nz)
                        G[r, r] += g
                        G[q, q] += g
                        G[r, q] -= g
                        G[q, r] -= g
                if z == 0:
                    G[r, r] += g_sink
    return G


_RTH_CACHE: dict = {}


def temperatures(xbar_vec, span: int, variant: str, kx=8, ky=8, kz=4) -> np.ndarray:
    p = THERMAL_VARIANTS[variant]
    key = (variant, kx, ky, kz)
    if key not in _RTH_CACHE:
        _RTH_CACHE[key] = np.linalg.inv(conductance(kx, ky, kz, p["g_lat"], p["g_vert"], p["g_sink"]))
    rate = np.asarray(xbar_vec, dtype=float) / span
    power = p["p_core"] + p["p_leak"] + p["e_flit"] * FREQ * rate
    return _RTH_CACHE[key] @ power
