#!/bin/bash
# Commitment counters under non-atomic VC reallocation (minimal adaptive routing with an escape VC, 2 VCs x 8 flits,
# 8-flit packets, uniform traffic, offered 1.0 flits/node/cycle, 6,000 cycles, seed 1).
# The in-house simulator prints on stderr: VA_GRANTS <g> NONEMPTY <n> FULL <f> FULL_ADAPTIVE <fa>.
# This script reproduces e_444_h1.txt, e_444_h4.txt, e_888_h1.txt and e_888_h4.txt exactly.
set -e
g++ -O3 -std=c++17 -o noc3d_sim_instr noc3d_sim_instr.cpp
for k in 4 8; do
  for h in 1 4; do
    cdelay=0; [ "$h" = 4 ] && cdelay=1
    ./noc3d_sim_instr --topo fc --routing fa --traffic random --pir 0.125 --seed 1 --cycles 6000 --warmup 1000 \
      --kx $k --ky $k --kz $k --vcs 2 --depth 8 --atomic 0 --hop-delay $h --credit-delay $cdelay --trace 0 \
      --out tmp_commit 2>&1 >/dev/null | grep VA_GRANTS > e_${k}${k}${k}_h${h}.txt
  done
done
rm -f tmp_commit.*
# BookSim 2 (public repository, commit 28f4329): replace src/routers/iq_router.cpp by
# booksim_iq_router_instrumented.cpp, rebuild, then
#   $BOOKSIM instr4.cfg 2>&1 >/dev/null | grep VA_GRANTS > instr4.err
#   $BOOKSIM instr8.cfg 2>&1 >/dev/null | grep VA_GRANTS > instr8.err
