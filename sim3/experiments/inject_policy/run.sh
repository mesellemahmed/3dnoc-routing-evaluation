#!/bin/bash
# Simulator variant with the --inject-escape option (confines half of the packets to the escape VC at injection).
# Build first: g++ -O3 -std=c++17 -o noc3d_sim_injesc noc3d_sim_injesc.cpp
SIM=${SIM_INJESC:-./noc3d_sim_injesc}
for size in "8 8 4" "4 4 4"; do
 set -- $size
 for ie in 0 1; do
  for tr in random transpose1; do
   for seed in 1 2 3 4; do
    p=o_${1}${2}${3}_ie${ie}_${tr}_${seed}
    $SIM --topo fc --routing fa --traffic $tr --pir 0.075 --seed $seed --cycles 12000 --warmup 2000 --kx $1 --ky $2 --kz $3 --vcs 2 --depth 8 --atomic 0 --inject-escape $ie --trace 0 --out $p > /dev/null 2>&1
    echo "$1x$2x$3 ie=$ie $tr $seed $(python3 -c "import json;d=json.load(open('$p.stats'));print(d['deadlock_cycle'],d['drained'],d['stuck_flits'])")"
    rm -f $p.trace $p.buf $p.pkt $p.act $p.topo
   done
  done
 done
done
