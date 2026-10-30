#!/bin/bash
# Re-run the channel-dependency checker for every routing class and compare with the committed verdicts.
# Usage: bash scripts/regenerate_cdg.sh [outdir]      (run from the repository root; default outdir: cdg_check)
set -e
ROOT=$(cd "$(dirname "$0")/.." && pwd)
OUT=${1:-cdg_check}
mkdir -p "$OUT"
g++ -O3 -std=c++17 -o "$OUT/noc3d_sim" "$ROOT/sim3/noc3d_sim.cpp"
g++ -O3 -std=c++17 -DNEGATIVE_TEST_EF_MIXED -o "$OUT/noc3d_sim_neg" "$ROOT/sim3/noc3d_sim.cpp"
for c in xyz zxy nf dyad ham o1t romm val fa rca qr ta; do
  "$OUT/noc3d_sim" --routing $c --check 1 > "$OUT/cdg_$c.json"
done
for el in stride2 stride4; do
  for c in ef efr aesg adele fl; do
    "$OUT/noc3d_sim" --topo vpc --elevators $el --routing $c --check 1 > "$OUT/cdg_vpc_${el}_$c.json"
  done
done
"$OUT/noc3d_sim_neg" --topo vpc --elevators stride2 --routing ef --check 1 > "$OUT/cdg_negative_ef_mixed.json"
python3 - "$ROOT/sim3/results" "$OUT" <<'PY'
import json, sys
from pathlib import Path
ref, new = Path(sys.argv[1]), Path(sys.argv[2])
bad = [f.name for f in sorted(ref.glob("cdg_*.json")) if json.loads(f.read_text()) != json.loads((new / f.name).read_text())]
print("checker verdicts identical to the committed ones" if not bad else f"DIFFERENT: {bad}")
PY
