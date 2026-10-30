#!/bin/bash
# Deadlock witness of Figure 3: minimal adaptive routing with an escape VC (fa), 8x8x4 mesh, uniform traffic at
# PIR 0.075 (offered 0.6 flits/node/cycle), non-atomic VC reallocation, seed 1: the run of group sens, variant na.
# --dump-stuck 1 writes one JSON line per non-empty input VC of the network that failed to drain (fa_na.stuck):
# router r, input port ip, VC v, allocated output port and VC, the packets in the buffer as [id, head flit
# present, destination], the candidates of an unallocated head, and the credits of the allocated output VC.
# Ports: 0 local, 1 east, 2 west, 3 north, 4 south, 5 up, 6 down; router r = x + 8 y + 64 z.
SIM=${SIM:-../../noc3d_sim}
$SIM --topo fc --routing fa --traffic random --pir 0.075 --seed 1 --cycles 12000 --warmup 2000 --kx 8 --ky 8 --kz 4 \
  --vcs 2 --depth 8 --atomic 0 --trace 0 --dump-stuck 1 --out fa_na > /dev/null 2>&1
rm -f fa_na.trace fa_na.buf fa_na.pkt fa_na.act fa_na.topo
