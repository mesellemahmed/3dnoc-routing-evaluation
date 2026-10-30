Commitment rates under non-atomic VC reallocation (minimal adaptive routing with escape VC, 2 VCs x 8 flits,
8-flit packets, uniform traffic, full injection). A "commitment to a full adaptive VC" is a VC grant on a network
output whose downstream adaptive VC has no free slot (the state of Proposition 1).

* BookSim 2 (commit 28f4329) instrumented in routers/iq_router.cpp (_VCAllocUpdate), batch mode:
  instr4.cfg / instr8.cfg -> instr4.err / instr8.err.
* In-house simulator instrumented in allocate() (noc3d_sim_instr.cpp), 6,000 cycles at 1.0 flits/node/cycle offered:
  e_444_h1.txt, e_444_h4.txt, e_888_h1.txt, e_888_h4.txt.

Reproduce: `bash run.sh` (in-house counters, exact) and the BookSim commands at the end of run.sh.
