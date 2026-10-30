#!/bin/bash
# Longer BookSim batch-mode drain tests (exposure comparable to the in-house overload runs).
# BOOKSIM: path to a BookSim 2 binary (commit 28f4329); exit code 255 is BookSim's normal completion code.
BS=${BOOKSIM:-../booksim}
run() { # k traffic wtc batch seed
  f=f_$1_$2_w$3_b$4_$5
  sed -e "s/^k = .*/k = $1;/" -e "s/^traffic = .*/traffic = $2;/" -e "s/^wait_for_tail_credit = .*/wait_for_tail_credit = $3;/" \
      -e "s/^batch_size = .*/batch_size = $4;/" -e "s/^seed = .*/seed = $5;/" batch.cfg > $f.cfg
  start=$(date +%s)
  timeout 1500 $BS $f.cfg > $f.log 2>&1
  rc=$?
  echo "k=$1 tr=$2 wtc=$3 batch=$4 seed=$5 rc=$rc t=$(( $(date +%s) - start ))s warn=$(grep -c 'Possible network deadlock' $f.log) dur=$(grep -m1 'Overall min batch duration' $f.log | awk '{print $6}')"
}
for s in 1 2 3 4; do run 8 uniform 0 600 $s; done
for s in 1 2 3 4; do run 4 uniform 0 2000 $s; done
