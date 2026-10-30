#!/bin/bash
# Non-atomic reallocation with deeper router pipelines: FA, 8x8x4, offered 0.6 flits/node/cycle.
SIM=${SIM:-../../noc3d_sim}
for h in 3 4; do
 for tr in random transpose1; do
  for seed in 1 2; do
   p=o_h${h}_${tr}_${seed}
   $SIM --topo fc --routing fa --traffic $tr --pir 0.075 --seed $seed --cycles 12000 --warmup 2000 --kx 8 --ky 8 --kz 4 --vcs 2 --depth 8 --atomic 0 --hop-delay $h --credit-delay 1 --trace 0 --out $p > /dev/null 2>&1
   echo "h=$h $tr $seed $(python3 -c "import json;d=json.load(open('$p.stats'));print(d['deadlock_cycle'],d['drained'],d['stuck_flits'])")"
   rm -f $p.trace $p.buf $p.pkt $p.act $p.topo
  done
 done
done
