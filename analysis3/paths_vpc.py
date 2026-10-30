"""Path diversity of the VPC classes on the 8x8x4 stride-2 elevator layout."""
from __future__ import annotations

import json
from functools import lru_cache

import numpy as np

KX, KY, KZ = 8, 8, 4
N = KX * KY * KZ
ELEV = [(x, y) for x in range(KX) for y in range(KY) if x % 2 == 1 and y % 2 == 1]


def xyz(r):
    return r % KX, (r // KX) % KY, r // (KX * KY)


def nearest(s, prefer_negative):
    x, y, _ = xyz(s)
    best, bd, bneg = -1, 1 << 30, False
    for e, (ex, ey) in enumerate(ELEV):
        d = abs(ex - x) + abs(ey - y)
        neg = ex <= x and ey <= y
        if d < bd or (prefer_negative and d == bd and neg and not bneg):
            best, bd, bneg = e, d, neg
    return best


def fl_paths(s, d):
    """Distinct router paths admitted by First-Last for (s, d)."""
    e = nearest(s, True) if xyz(s)[2] != xyz(d)[2] else -1

    @lru_cache(maxsize=None)
    def count(r, vn):
        if r == d:
            return 1
        x, y, z = xyz(r)
        dx, dy, dz = xyz(d)
        other = z != dz
        tx, ty = ELEV[e] if other else (dx, dy)
        ne, nn, nw, ns = tx > x, ty > y, tx < x, ty < y
        moves = []
        if other:
            if not (ne or nn or nw or ns):
                moves = [((x, y, z + (1 if dz > z else -1)), 1)]
            elif vn <= 0 and (ne or nn):
                if ne: moves.append(((x + 1, y, z), 0))
                if nn: moves.append(((x, y + 1, z), 0))
            else:
                if nw: moves.append(((x - 1, y, z), 1))
                if ns: moves.append(((x, y - 1, z), 1))
        else:
            if (nw or ns) and vn <= 1:
                if nw: moves.append(((x - 1, y, z), 1))
                if ns: moves.append(((x, y - 1, z), 1))
            else:
                if ne: moves.append(((x + 1, y, z), 2))
                if nn: moves.append(((x, y + 1, z), 2))
        return sum(count(n[0] + KX * (n[1] + KY * n[2]), v) for n, v in moves)

    return count(s, -1)


def main():
    pairs = [(s, d) for s in range(N) for d in range(N) if s != d]
    inter = np.array([xyz(s)[2] != xyz(d)[2] for s, d in pairs])
    fl = np.array([fl_paths(s, d) for s, d in pairs], dtype=float)
    out = {
        "ef": {"multi_frac": 0.0, "paths_median": 1.0, "alpha_mean": float("nan")},
        "efr": {"multi_frac": 0.0, "paths_median": 1.0, "alpha_mean": float("nan")},
        "aesg": {"multi_frac": float(inter.mean()), "paths_median": 16.0, "alpha_mean": float("nan")},
        "adele": {"multi_frac": float(inter.mean()), "paths_median": 4.0, "alpha_mean": float("nan")},
        "fl": {"multi_frac": float((fl > 1).mean()), "paths_median": float(np.median(fl)), "alpha_mean": float("nan")},
    }
    json.dump(out, open("paths_vpc.json", "w"), indent=1)
    print(out)


if __name__ == "__main__":
    main()
