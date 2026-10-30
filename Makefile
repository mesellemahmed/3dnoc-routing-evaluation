# Reproduction targets for "Routing in Three-Dimensional Networks-on-Chip" (run from the repository root).
PY       ?= python3
CXX      ?= g++
CXXFLAGS ?= -O3 -march=native -std=c++17
WORKERS  ?= 2
N        ?= 20

.PHONY: all sim numbers check-verdicts reproduce-sample cdg paper campaign refine analytical hw clean

all: sim numbers check-verdicts

sim: sim3/noc3d_sim

sim3/noc3d_sim: sim3/noc3d_sim.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

# Regenerate every number, table and figure input of the article from the raw results.
numbers:
	cd paper3 && $(PY) build_paper3.py && $(PY) assemble.py

# Recompute the verdict table independently of the article's generator.
check-verdicts:
	$(PY) scripts/check_verdicts.py

# Re-run N randomly chosen simulations and compare their statistics bit for bit.
reproduce-sample: sim
	$(PY) scripts/reproduce_sample.py $(N)

# Re-run the channel-dependency checker for every class and the negative test.
cdg:
	bash scripts/regenerate_cdg.sh cdg_check

paper: numbers
	mkdir -p paper3/build && cp paper3/routing_3dnoc_survey_v3.tex paper3/build/
	cd paper3/build && pdflatex -interaction=nonstopmode -halt-on-error routing_3dnoc_survey_v3.tex > /dev/null \
	  && pdflatex -interaction=nonstopmode -halt-on-error routing_3dnoc_survey_v3.tex > /dev/null \
	  && pdflatex -interaction=nonstopmode -halt-on-error routing_3dnoc_survey_v3.tex > /dev/null
	@echo "PDF: paper3/build/routing_3dnoc_survey_v3.pdf"

# Full campaign from scratch (does not touch the committed results).
campaign: sim
	cd sim3 && $(PY) campaign3.py --workers $(WORKERS) --out results/campaign3_rerun.jsonl

refine: sim
	cd sim3 && $(PY) experiments/sat_refine/run_refine.py --rounds 5 --workers $(WORKERS) \
	  --coarse results/campaign3_rerun.jsonl --out results/sweep_refine_rerun.jsonl

# Analytical models (deterministic; outputs are regenerated in place).
analytical:
	cd bench3 && $(CXX) $(CXXFLAGS) -o fault_reach fault_reach.cpp && ./fault_reach 2000 > fault_reach.csv
	cd analysis3 && $(PY) partition.py && $(PY) paths.py > paths_result.json && $(PY) paths_vpc.py

# Synthesis (needs YoWASP Yosys and hw/NangateOpenCellLibrary_typical.lib), functional check and sensitivity of the
# areas to the synthesis context.
hw:
	cd hw && for m in rc_xyz rc_zxy rc_nf rc_dyad rc_ham rc_o1t rc_twophase rc_fa rca_select rca_agg qr_select \
	  qr_update rc_ta rc_ef rc_fl adele_sel vc_fifo xbar7; do \
	  bash synth.sh $$m; done && $(PY) check_units.py | tee check_units.log && $(PY) synth_context.py | tee synth_context.log

clean:
	rm -rf sim3/noc3d_sim cdg_check paper3/build bench3/fault_reach
	find . -name '__pycache__' -prune -exec rm -rf {} +
