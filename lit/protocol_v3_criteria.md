# Evaluation protocol v3 — pass criteria fixed before the campaign results
(written 2026-10-25, before any v3 result beyond the first 153 sweep runs was inspected)

Notation: "class" = routing class; XYZ is the reference for relative criteria.
All simulation statistics: means over seeds; 95% Student-t intervals; paired
Wilcoxon signed-rank tests with Holm correction; effect sizes as Cliff's delta.

D1 Correctness. PASS iff (i) the channel-dependency checker proves deadlock
freedom (acyclic CDG, or Duato: acyclic escape subnetwork without return and
connected escape) AND (ii) no run of the class, in any group, fails the
watchdog or the drain test. FAIL otherwise.

D2 Path diversity. PASS iff the routing relation admits more than one path for
at least 50% of ordered source-destination pairs (fault-free mesh).

D3 Performance to saturation. Saturation throughput = highest accepted
throughput among sweep points with undelivered fraction <= 1%. PASS iff
saturation throughput >= 0.9 x XYZ on uniform random AND > XYZ on at least one
of transpose1, bit-reversal.

D4 Traffic robustness. Over the ERA grid (5 PIR x 5 seeds), for each of the six
patterns the class's mean latency must be <= 2 x the best class's mean latency
on that pattern (within the same topology family). PASS iff true for all six.

D5 Fault tolerance. FC classes: PASS iff exact pair reachability exceeds XYZ's
by >= 0.05 at p_h = 0.05 and at p_v = 0.2 (xi = 0 and xi = 4); PARTIAL if it
exceeds XYZ in some but not all of these; FAIL otherwise. VPC classes: PASS
iff no packet is lost (drop ratio 0) in every vpcfault configuration.

D6 Load balance and fairness. PASS iff mean router-load Jain index >= XYZ's on
uniform random (ERA grid) AND per-source accepted-throughput Jain >= 0.95 in
every stress run on uniform random that did not deadlock.

D7 Structure (ERA). Descriptive, no pass/fail.

D8 Hardware cost. Synthesised route/selection logic (NanGate 45 nm) plus
storage, relative to a 2-VC baseline router. PASS iff overhead <= 5% of the
baseline router area and no global information network is required.

D9 Energy. Parametric model (27 parameter combinations). PASS iff energy per
delivered flit <= 1.05 x XYZ's in every combination (uniform random,
PIR 0.006).

D10 Thermal. Compact steady-state model. PASS iff peak temperature rise <= XYZ's
+ 0.05 K under uniform random and hotspot traffic (PIR 0.006) in every thermal
parameter variant.

D11 Security (flooding DoS). PASS iff the benign-traffic latency inflation
factor under flooding is <= XYZ's at both benign loads.

D12 Robustness of conclusions. PASS iff the sign of the class's saturation-
throughput difference with XYZ (uniform random and transpose1) is the same in
all eight router variants of the sensitivity study.
