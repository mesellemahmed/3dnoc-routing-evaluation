#!/bin/bash
# Simulator variant with the --rand-sel option (equal-priority random choice among adaptive candidates, as in BookSim).
# Build first: g++ -O3 -std=c++17 -o noc3d_sim_randsel noc3d_sim_randsel.cpp
SIM=${SIM_RANDSEL:-./noc3d_sim_randsel}
for cfg in "8 8 4 1 0" "8 8 8 4 1"; do
 set -- $cfg
 for seed in 1 2 3 4; do
  for tr in random transpose1; do
   p=o_${1}${2}${3}_h$4_${tr}_${seed}
   $SIM --topo fc --routing fa --traffic $tr --pir 0.075 --seed $seed --cycles 12000 --warmup 2000 --kx $1 --ky $2 --kz $3 --vcs 2 --depth 8 --atomic 0 --hop-delay $4 --credit-delay $5 --rand-sel 1 --trace 0 --out $p > /dev/null 2>&1
   echo "$1x$2x$3 h=$4 $tr $seed $(python3 -c "import json;d=json.load(open('$p.stats'));print(d['deadlock_cycle'],d['drained'],d['stuck_flits'])")"
   rm -f $p.trace $p.buf $p.pkt $p.act $p.topo
  done
 done
done
