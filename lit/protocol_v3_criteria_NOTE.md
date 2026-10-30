# History and status of the pass criteria (`protocol_v3_criteria.md`)

* The criteria file was written by the authors during the v3 campaign. The campaign started on 24 October 2026 at
  23:27, after two aborted launches whose runs are kept in `sim3/results/campaign3_aborted_v*.jsonl` (153 load-sweep
  runs of XYZ and ZXY, and 138 ERA-grid runs of XYZ; none of them is used in the article), and most runs had
  completed by 01:37 on 25 October; the file is dated 25 October 2026, 06:07. Its header states that no result
  beyond a 153-run pilot of the load sweep had been inspected when it was written. These times and this statement
  cannot be verified from the artefact, and the criteria were not deposited with any external registry. The
  article therefore does **not** claim pre-registration: it presents the criteria as fixed analysis rules, applied
  identically to every class, and reports the verdict under halved and doubled thresholds.
* The file is kept unmodified. One rule that it does not contain was added when the criteria were applied: a class
  that deadlocks in any run of a dimension fails that dimension. (The analysis code also takes the "best class" of
  D4 among the classes that never deadlocked in the D4 runs; no D4 run deadlocked, so this has no effect.)
* Other departures, all listed in Section "Protocol deviations" of the article: D12 uses overload throughput as a
  proxy for saturation throughput and exact differences; the effect size is the paired dominance statistic instead
  of Cliff's delta; EF-R is evaluated on the fault experiments only; the threshold-sensitivity analysis is post hoc;
  the bisection refinement of the load sweep was added after an internal audit; and D1 is applied as written
  (non-atomic runs counted), which an earlier version of the analysis did not do.
