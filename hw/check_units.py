"""Functional spot check of the synthesised route units against reference
models that mirror the simulator's routing relation (yosys `eval`)."""
from __future__ import annotations

import random
import re
import subprocess
import sys

KX, KY, KZ = 8, 8, 4
LOCAL, E, W, N, S, U, D = range(7)


def onehot(p: int) -> int:
    return 1 << p


def dor(c, d, order):
    for dim in order:
        delta = d[dim] - c[dim]
        if delta:
            return [(E, W), (N, S), (U, D)][dim][0 if delta > 0 else 1]
    return LOCAL


def ham_label(x, y, z):
    idx = y * KX + x if y % 2 == 0 else y * KX + (KX - 1 - x)
    if z % 2 == 1:
        idx = KX * KY - 1 - idx
    return z * KX * KY + idx


def ham_ref(c, d):
    if c == d:
        return LOCAL
    lc, ld = ham_label(*c), ham_label(*d)
    nbrs = {E: (c[0] + 1, c[1], c[2]), W: (c[0] - 1, c[1], c[2]), N: (c[0], c[1] + 1, c[2]),
            S: (c[0], c[1] - 1, c[2]), U: (c[0], c[1], c[2] + 1), D: (c[0], c[1], c[2] - 1)}
    adm = []
    for p in (E, W, N, S, U, D):
        n = nbrs[p]
        if not (0 <= n[0] < KX and 0 <= n[1] < KY and 0 <= n[2] < KZ):
            continue
        ln = ham_label(*n)
        if (ld > lc and lc < ln <= ld) or (ld < lc and ld <= ln < lc):
            adm.append(p)
    dist = lambda a: sum(abs(a[i] - d[i]) for i in range(3))
    red = [p for p in adm if dist(nbrs[p]) < dist(c)]
    if red:
        return red[0]
    return min(adm, key=lambda p: (abs(ham_label(*nbrs[p]) - ld), p))


def run_eval(module: str, vectors: list[dict], outputs: list[str]) -> list[dict]:
    cmds = [f"read_verilog route_units.v", f"prep -top {module}", "flatten"]
    for v in vectors:
        sets = " ".join(f"-set {k} {val}" for k, val in v.items())
        shows = " ".join(f"-show {o}" for o in outputs)
        cmds.append(f"eval {sets} {shows}")
    out = subprocess.run(["yowasp-yosys", "-p", "; ".join(cmds)], capture_output=True, text=True).stdout
    vals = re.findall(r"Eval result: \\(\w+) = (\d+)'([01xz]+)\.", out)
    res, k = [], 0
    for _ in vectors:
        row = {}
        for o in outputs:
            name, width, bits = vals[k]
            row[name] = int(bits, 2)
            k += 1
        res.append(row)
    return res


def main() -> None:
    rng = random.Random(7)
    pts = [((rng.randrange(KX), rng.randrange(KY), rng.randrange(KZ)),
            (rng.randrange(KX), rng.randrange(KY), rng.randrange(KZ))) for _ in range(60)]
    base = lambda c, d: {"cx": c[0], "cy": c[1], "cz": c[2], "dx": d[0], "dy": d[1], "dz": d[2]}
    failures = 0
    checks = {
        "rc_xyz": lambda c, d: dor(c, d, (0, 1, 2)),
        "rc_zxy": lambda c, d: dor(c, d, (2, 0, 1)),
    }
    for mod, ref in checks.items():
        out = run_eval(mod, [base(c, d) for c, d in pts], ["port"])
        bad = sum(o["port"] != onehot(ref(c, d)) for o, (c, d) in zip(out, pts))
        print(f"{mod}: {len(pts) - bad}/{len(pts)} match")
        failures += bad
    vec = [{**base(c, d), "order": o} for (c, d), o in zip(pts, [i % 2 for i in range(len(pts))])]
    out = run_eval("rc_o1t", vec, ["port"])
    bad = sum(o["port"] != onehot(dor(c, d, (0, 1, 2) if v["order"] == 0 else (2, 1, 0)))
              for o, (c, d), v in zip(out, pts, vec))
    print(f"rc_o1t: {len(pts) - bad}/{len(pts)} match")
    failures += bad
    vec = [{**base(c, d), "credits": 0} for c, d in pts]
    out = run_eval("rc_ham", vec, ["port"])
    bad = 0
    for o, (c, d) in zip(out, pts):
        if o["port"] != onehot(ham_ref(c, d)):
            bad += 1
            print("  ham mismatch", c, d, bin(o["port"]), ham_ref(c, d))
    print(f"rc_ham: {len(pts) - bad}/{len(pts)} match")
    failures += bad
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
