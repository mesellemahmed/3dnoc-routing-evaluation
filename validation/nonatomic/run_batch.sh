#!/bin/bash
# BookSim 2 batch mode = drain test: every node injects batch_size packets at full rate; the run ends only when
# all packets are delivered. A partial deadlock leaves packets undelivered -> "Possible network deadlock" warnings
# and a timeout.
# BOOKSIM: path to a BookSim 2 binary (public repository, commit 28f4329). BookSim returns exit code 255 (-1) on
# normal completion; a run is complete when its log reports the batch duration.
BS=${BOOKSIM:-../booksim}
run() { # k traffic wtc batch seed
  f=e_$1_$2_w$3_$5
  sed -e "s/^k = .*/k = $1;/" -e "s/^traffic = .*/traffic = $2;/" -e "s/^wait_for_tail_credit = .*/wait_for_tail_credit = $3;/" \
      -e "s/^batch_size = .*/batch_size = $4;/" -e "s/^seed = .*/seed = $5;/" batch.cfg > $f.cfg
  start=$(date +%s)
  timeout 600 $BS $f.cfg > $f.log 2>&1
  rc=$?
  echo "k=$1 tr=$2 wtc=$3 batch=$4 seed=$5 rc=$rc t=$(( $(date +%s) - start ))s warn=$(grep -c 'Possible network deadlock' $f.log) dur=$(grep -m1 'Overall min batch duration' $f.log | awk '{print $6}')"
}
for s in 1 2 3 4 5; do
  run 4 uniform 0 400 $s
  run 4 transpose 0 400 $s
  run 4 bitrev 0 400 $s
  run 8 uniform 0 150 $s
  run 8 bitrev 0 150 $s
done
