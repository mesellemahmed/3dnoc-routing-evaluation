#!/bin/bash
# Uses the main simulator (build it in sim3/ first).
SIM=${SIM:-../../noc3d_sim}
for size in "4 4 4" "8 8 8" "8 8 4"; do
 set -- $size
 for tr in random transpose1; do
  for seed in 1 2 3; do
   p=o_${1}${2}${3}_${tr}_${seed}
   $SIM --topo fc --routing fa --traffic $tr --pir 0.075 --seed $seed --cycles 12000 --warmup 2000 --kx $1 --ky $2 --kz $3 --vcs 2 --depth 8 --atomic 0 --trace 0 --out $p > /dev/null 2>&1
   echo "$1x$2x$3 $tr $seed $(python3 -c "import json;d=json.load(open('$p.stats'));print(d['deadlock_cycle'],d['drained'],d['stuck_flits'])")"
   rm -f $p.trace $p.buf $p.pkt $p.act
  done
 done
done
