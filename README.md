# Routing in 3D Networks-on-Chip — code, data and scripts

Artefact of the article

> **Beyond Performance: A Critical Survey and Verified Twelve-Dimensional Evaluation of 3D Network-on-Chip Routing**
> Ahmed Mesellem & Mohammed Mana — manuscript under review (2026).

This repository contains everything behind the article's numbers. The pieces are:

- the cycle-level simulator and its channel-dependency (CDG) deadlock checker;
- the analytical fault-reachability and partition models;
- the Verilog route units and the synthesis scripts;
- the campaign drivers, every raw result (7,864 simulations, 28,000 fault maps) and every supplementary experiment;
- the scripts that regenerate every number, table and figure input of the article.

Every number in the article's tables, and every measured value that its text quotes as a number, is a LaTeX macro
`\val{key}` written by `paper3/build_paper3.py` from the raw data. Counts and rounded comparisons written as words
(for example "a six-channel cycle", "about half the rate") and the router coordinates of the deadlock witness in
Figure 3 are written by hand; the witness comes from `sim3/experiments/witness/fa_na.stuck`.

---

## Checks a reviewer can run

All commands run from the repository root (requirements below). Times were measured on a 2-core machine.

| What is checked | Command | Time |
|---|---|---|
| Rebuild every number and table of the article from the raw results; the output should be identical to the committed one | `make numbers && git status --short paper3` (no change expected) | ~10 s |
| Recompute the verdict table (Table 32) directly from the raw data, independently of the article's generator | `make check-verdicts` | ~5 s |
| Re-run a random sample of the 7,864 simulations and compare the simulator statistics bit for bit | `make reproduce-sample N=20` | ~15 s |
| Re-run the deadlock checker for all 17 classes and the negative test, and compare with the stored verdicts | `make cdg` | ~10 s |
| Rebuild the PDF of the article | `make paper` (needs pdflatex) | ~15 s |
| Re-run the analytical models (fault reachability, partition risk, path counts); outputs are regenerated in place and should be identical | `make analytical && git status --short` | ~12 min |
| Re-run the whole campaign from scratch | `make campaign WORKERS=8` then `make refine WORKERS=8` | ~4 CPU-hours |

Before release, every check of this table except the full campaign was run on a fresh copy of the repository:
the regenerated numbers, tables, analytical outputs, checker verdicts and PDF were identical to the committed ones,
and 92 randomly chosen simulations of the campaign and of the refinement were reproduced bit for bit.

---

## Repository layout

| Path | Content |
|---|---|
| `sim3/noc3d_sim.cpp` | Cycle-level wormhole/VC simulator for fully connected 3D meshes and vertically partially connected (VPC) stacks. It covers 17 routing classes and includes the CDG checker (`--check 1`). One C++17 file, no dependencies. |
| `sim3/era_analyzer.py` | Classical metrics and ERA structural metrics computed from the simulator outputs. ERA is the method of Mesellem & Mana, *Algorithms* 19(9):765, 2026. |
| `sim3/campaign3.py` | Resumable campaign driver: 12 run groups, 7,144 runs (`--list` prints the groups). |
| `sim3/experiments/sat_refine/` | Bisection refinement of the saturation loads (groups `sweepref`, `sweepref_vpc`, 720 runs). It was added after an internal audit and is idempotent: re-running it on complete results launches no run. |
| `sim3/results/campaign3.jsonl` | Raw results of the campaign, one JSON object per run. |
| `sim3/results/sweep_refine.jsonl` | Raw results of the refinement. The article uses this file and `campaign3.jsonl` together. |
| `sim3/results/campaign3_aborted_v*.jsonl` | Runs of two aborted launches: 153 load-sweep runs (the pilot mentioned in `lit/protocol_v3_criteria_NOTE.md`) and 138 ERA-grid runs. They are kept to document the timeline of the pass criteria; the article does not use them. |
| `sim3/results/cdg_*.json` | Checker verdicts and shortest dependency cycles for every class, plus the negative test (`cdg_negative_ef_mixed.json`). |
| `sim3/experiments/` | Supplementary experiments, each with its sources, scripts and raw outputs (see below). |
| `bench3/` | Exact pair reachability of 12 routing policies, and of the connectivity bound, under planar faults and Gaussian-copula vertical faults: `fault_reach.cpp` and `fault_reach.csv`, 15 configurations (14 with faults) × 2,000 maps. |
| `analysis3/` | Shared analysis code: results loader, statistics, energy and thermal models, path enumeration (D2), partition Monte Carlo (Theorem 3), and the coded corpus and bibliography of the survey. |
| `hw/` | Verilog route and selection units (`route_units.v`), Yosys synthesis script for NanGate 45 nm, area reports and logs. Also the functional check of four route units against reference models (`check_units.py`, `check_units.log`) and the sensitivity of the areas to the synthesis context (`synth_context.py`, `synth_context.log`). |
| `validation/` | BookSim 2 cross-validation (`cross_validate.py`, `results/validation.json`) and the BookSim batch-mode drain tests under non-atomic reallocation (`nonatomic/`). |
| `paper3/` | The article itself (described in the next three rows). |
| `paper3/p01`–`p07*.tex` | Hand-written LaTeX parts of the article. |
| `paper3/build_paper3.py`, `sec_*.py`, `common.py` | Generators that write `gen/` and `generated_values.json` from the raw data. |
| `paper3/assemble.py` | Assembles the standalone `routing_3dnoc_survey_v3.tex`. It refuses to build if a macro or a citation is undefined. |
| `lit/` | The pass criteria as written during the campaign (`protocol_v3_criteria.md`, unmodified). The note `protocol_v3_criteria_NOTE.md` gives their history, the rule added when applying them and every deviation. |
| `scripts/` | Reviewer tools: `check_verdicts.py`, `reproduce_sample.py` and `regenerate_cdg.sh`. |

### Where each table and figure comes from

Table and figure numbers refer to the submitted PDF. All generator functions are in `paper3/`. Tables 1, 5, 6, 7
and 34 and the drawings of Figures 1–3 are written by hand in `paper3/p0*.tex`; the counts of Fig. 2 and the
results quoted in Table 34 are macros.

| Article | Generator | Input data |
|---|---|---|
| Tables 2, 3 (corpus, historical designs), Fig. 2 counts, Table 4 (audit) | `sec_survey.py` | `analysis3/corpus.py` (coded attributes) |
| Table 8 (campaigns) | `sec_platform.section_counts` | `campaign3.jsonl`, `sweep_refine.jsonl` |
| Table 9 (BookSim validation) | `sec_platform.section_validation` | `validation/results/validation.json` |
| Table 10 (CDG verdicts) | `sec_platform.section_cdg` | `sim3/results/cdg_*.json` |
| Table 11 (partition risk) | `sec_results.partition` | `analysis3/partition.json` |
| Table 12 (D1) | `sec_results.d1` | `campaign3.jsonl`, `cdg_*.json` |
| Fig. 3 (deadlock witness) | drawn by hand | `sim3/experiments/witness/fa_na.stuck` (re-run of the `sens` run FA, uniform, non-atomic, overload, seed 1) |
| Non-atomic and BookSim paragraphs (Section 8.1) | `sec_platform.section_nonatomic` | `validation/nonatomic/summary.json` (rebuilt by `build_summary.py` from the raw outputs) |
| Table 13 (D2) | `sec_results.d2` | `analysis3/paths_result.json`, `paths_vpc.json` |
| Tables 14, 15 (D3), Fig. 4 (latency–load curves) | `sec_results.d3` | `campaign3.jsonl` (sweep groups), `sweep_refine.jsonl` |
| Table 16 (four elevators) | `sec_results.vpc4` | `campaign3.jsonl` (group `vpc4`) |
| Tables 17–19 (D4, bursty traffic) | `sec_results.d4` | `campaign3.jsonl` (groups `era`, `era_vpc`, `burst`) |
| Table 20 (D6) | `sec_results.d6` | groups `era`, `stress` |
| Tables 21, 22 (D12, scalability) | `sec_results.d12` | groups `sens`, `scal` |
| Tables 23, 24 (D5) | `sec_results.d5` | `bench3/fault_reach.csv`, group `vpcfault` |
| VPC fault-path check (Section 8.3) | `sec_results.faultcheck` | `sim3/experiments/fault_check/summary.json` (rebuilt by `summarize.py`) |
| Table 25 (ERA signatures) | `sec_results.d4`, `sec_results.d7` | traced runs (`era` fields) |
| Table 26 (D8) | `sec_results.d8` | `hw/rpt_*.txt` |
| Tables 27, 28 (D9, D10) | `sec_results.d9_d10` | groups `era`, `era_vpc` (activity counters); `analysis3/energy_thermal.py` |
| Table 29 (D11) | `sec_results.d11` | group `dos` |
| Table 30 (out-of-order delivery) | `sec_results.reorder` | `sim3/experiments/reorder/reorder.json` |
| Table 31 (paired statistics) | `sec_results.paired` | groups `era`, `era_vpc`, `dos`; `analysis3/stats.py` |
| Table 32 (verdict), Table 33 (threshold sensitivity) | `sec_results.verdict`, `threshold_sensitivity` | all of the above |

---

## Requirements

- g++ ≥ 11 (C++17). Tested with g++ 13.3.
- Python ≥ 3.10 with numpy, scipy and pandas (`pip install -r requirements.txt`). Tested with Python 3.13, numpy 2.5, scipy 1.18 and pandas 3.0.
- To build the PDF: pdflatex with TikZ/pgfplots, cleveref, algorithmicx, threeparttable and enumitem (TeX Live 2023 or later).
- Optional:
  - BookSim 2 for the validation and drain tests. Use the public repository at commit `28f4329` and point the scripts to the binary with `BOOKSIM=/path/to/booksim`. BookSim is not redistributed here.
  - YoWASP Yosys ≥ 0.40 and the NanGate 45 nm Open Cell Library (`NangateOpenCellLibrary_typical.lib`) for synthesis. The library is not redistributed here; place it in `hw/`.

## Reproducing everything from scratch

```bash
make sim                                   # build sim3/noc3d_sim
make cdg                                   # deadlock verification of every class (+ negative test)
make campaign WORKERS=8                    # 7,144 runs -> sim3/results/campaign3_rerun.jsonl
make refine WORKERS=8                      # bisection refinement of the saturation loads
make analytical                            # fault reachability (~10 min), partition Monte Carlo, path counts (~1 min)
make numbers                               # regenerate paper3/gen, generated_values.json and the .tex
make paper                                 # compile the PDF into paper3/build/
```

`make campaign` and `make refine` write to `campaign3_rerun.jsonl` and `sweep_refine_rerun.jsonl`, so that the committed results are not overwritten; to analyse a re-run, point `paper3/sec_platform.load_all()` to these files. `make analytical` regenerates its outputs in place (they are deterministic; `git diff` shows any difference). Runs are
deterministic: traffic and routing draw from separate pseudo-random streams seeded from the configuration, and
the seed of the ERA permutation test is derived from the configuration key. A re-run therefore gives the same
statistics with the same compiler and flags. The committed results were produced with the flags of the `Makefile`
(`-O3 -march=native`) on an x86-64 processor with FMA instructions; with other flags or processors the compiler may
contract floating-point operations differently, and classes whose decisions read floating-point state (for example
the thermal-aware class) can then diverge after many cycles and agree only statistically.

### Supplementary experiments (`sim3/experiments/`)

| Directory | Purpose (section of the article) | Reproduce |
|---|---|---|
| `inject_policy/` | Non-atomic deadlock with and without confining half of the packets to the escape VC at injection, and on three mesh sizes (8.1) | build `noc3d_sim_injesc.cpp`, then `bash run.sh`, `bash run_sizes.sh` |
| `pipeline/` | Non-atomic deadlock with three- and four-cycle hops (8.1) | `bash run.sh`, `bash run2.sh` |
| `rand_selection/` | Equal-priority choice among adaptive candidates, as in BookSim (8.1) | build `noc3d_sim_randsel.cpp`, then `bash run.sh` |
| `commitment/` | Commitment counters in our simulator and in an instrumented BookSim (8.1) | `bash run.sh` (reproduces `e_*.txt` exactly); BookSim commands at its end |
| `fault_check/` | Exact Elevator-First loss per fault map, validating the VPC fault path (8.3) | `g++ -O2 -std=c++17 -o replica replica.cpp`, `./replica 20000 > replica_20000maps.txt`, `./replica 5 > replica_5maps.txt`, `python3 summarize.py --write` |
| `reorder/` | Out-of-order delivery per flow (8.5) | `python3 run_reorder.py` |
| `sat_refine/` | Bisection refinement of saturation (8.2, D3) | `python3 run_refine.py --rounds 5` (from `sim3/`) |
| `witness/` | Stuck-VC dump of the non-atomic run whose deadlock cycles are drawn in Figure 3 (8.1) | `bash run.sh` (reproduces `fa_na.stuck` and `fa_na.stats`) |

The three simulator variants (`noc3d_sim_injesc.cpp`, `noc3d_sim_randsel.cpp`, `noc3d_sim_instr.cpp`) are
`sim3/noc3d_sim.cpp` plus the options they add; `diff` against the main file shows the changes. Without these
options they behave exactly like the main simulator. `validation/nonatomic/build_summary.py --write` rebuilds
`summary.json` from all of these raw outputs.

---

## Simulator

```text
sim3/noc3d_sim --topo fc|vpc --elevators stride2|stride4 --routing <class>
  --traffic random|transpose1|shuffle|bitreversal|butterfly|hotspot|selfsim|dos
  --pir <packets/node/cycle> --kx --ky --kz --vcs --depth --flits --atomic 0|1 --hop-delay --credit-delay
  --cycles --warmup --window --trace 0|1 --pv --xi --attackers --attack-pir --seed --out <prefix>
  --check 1          # channel-dependency verification instead of simulation (JSON on stdout)
```

- **Classes, 3D-FC:** `xyz zxy nf dyad ham o1t romm val fa rca qr ta`.
- **Classes, VPC:** `ef efr aesg adele fl`.
- **Negative test of the checker:** compile with `-DNEGATIVE_TEST_EF_MIXED`, then run `--topo vpc --routing ef --check 1`.
- **Shared code:** the checker calls the same routing-relation function as the simulator.
- **Deadlock detection:** a watchdog (no flit moves for 2,000 cycles) and a drain test (the network must empty within 200,000 cycles after injection stops).

### Raw result records (`*.jsonl`)

Each line describes one run:

| Field | Content |
|---|---|
| `config` | Every parameter, including `group` and `seed`. The field `bf_period` is a constant (100) that the simulator does not read; it is kept because the run key, and with it the seed of the ERA permutation test, contains it. |
| `key` | Unique key of the run. |
| `stats` | Simulator statistics: deadlock cycle, drained flag, stuck flits, planar and vertical flit hops, Q-table updates, VC-allocation waits. |
| `classical` | Latency distribution, accepted and offered throughput, undelivered and dropped fractions, per-source Jain indices. |
| `activity` | Crossbar, planar and vertical traversals, per-router vectors. |
| `era` | Entropies, effective supports, space–time mutual information and permutation p-value, Jensen–Shannon drift, imbalance indicators. Present for traced runs only. |

`analysis3/load.py` flattens the records into a pandas DataFrame.

---

## Pass criteria and deviations

`lit/protocol_v3_criteria.md` holds the twelve pass criteria exactly as they were written during the campaign;
they were **not** pre-registered. `lit/protocol_v3_criteria_NOTE.md` gives their history and the rule added
when they were applied. The article lists every departure from the written criteria in Section 10.2, "Protocol
deviations".

## License

- **Code** (all source files and scripts): MIT License, see `LICENSE`.
- **Data** (`sim3/results/`, the outputs in `sim3/experiments/`, `bench3/fault_reach.csv`, `analysis3/*.json`, `validation/**/*.json`, `hw/*.txt`, `hw/*.log`, `hw/synth_context.json`): Creative Commons Attribution 4.0 International, see `LICENSE-DATA.md`.
- `sim3/experiments/commitment/booksim_iq_router_instrumented.cpp` is a modified BookSim 2 file and remains under BookSim's BSD license (header kept).

## Citation

Please cite the article (see `CITATION.cff`).