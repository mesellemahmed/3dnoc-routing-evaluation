#!/bin/bash
# More runs with a four-cycle hop (BookSim-like pipeline depth): FA, non-atomic, offered 0.6 flits/node/cycle.
SIM=${SIM:-../../noc3d_sim}
for size in "8 8 4" "8 8 8"; do
 set -- $size
 for tr in random transpose1; do
  for seed in 3 4 5 6; do
   p=o2_${1}${2}${3}_${tr}_${seed}
   $SIM --topo fc --routing fa --traffic $tr --pir 0.075 --seed $seed --cycles 12000 --warmup 2000 --kx $1 --ky $2 --kz $3 --vcs 2 --depth 8 --atomic 0 --hop-delay 4 --credit-delay 1 --trace 0 --out $p > /dev/null 2>&1
   echo "h=4 $1x$2x$3 $tr $seed $(python3 -c "import json;d=json.load(open('$p.stats'));print(d['deadlock_cycle'],d['drained'],d['stuck_flits'])")"
   rm -f $p.trace $p.buf $p.pkt $p.act $p.topo
  done
 done
done
