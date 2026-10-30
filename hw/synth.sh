#!/bin/bash
# Synthesise one top module onto NanGate 45 nm and print the chip area.
top=$1
yowasp-yosys -q -p "read_verilog route_units.v; synth -top $top -flatten; dfflibmap -liberty NangateOpenCellLibrary_typical.lib; abc -liberty NangateOpenCellLibrary_typical.lib; opt_clean; tee -o rpt_$top.txt stat -liberty NangateOpenCellLibrary_typical.lib" > log_$top.txt 2>&1
grep -E "Chip area|Number of cells" rpt_$top.txt | head -3 | tr '\n' ' '
echo
