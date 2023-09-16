#!/bin/bash

BASE_DIR="/home/psimmerl/LLP/mdc_analysis/reports/weekly/2023-09-14"
# alias python="/home/psimmerl/mambaforge/envs/ROOT/bin/python"

if [ $# -ge 1 ]; then
    NEV="$1"
else
    NEV="-1"
fi

if [ $# -ge 2 ]; then
    LAB="_$2"
else
    LAB=""
fi

# ************ #

LAB="DTOOT"

pkill python
time /home/psimmerl/mambaforge/envs/ROOT/bin/python -u train_bdt.py "$NEV" "standard_$LAB" "$LAB" | tee "$BASE_DIR/standard_$LAB.out"

pkill python
time /home/psimmerl/mambaforge/envs/ROOT/bin/python -u train_bdt.py "$NEV" "halo_$LAB" "$LAB" | tee "$BASE_DIR/halo_$LAB.out"

pkill python
time /home/psimmerl/mambaforge/envs/ROOT/bin/python -u train_bdt.py "$NEV" "dphi_$LAB" "$LAB" | tee "$BASE_DIR/dphi_$LAB.out"

pkill python
time /home/psimmerl/mambaforge/envs/ROOT/bin/python -u train_bdt.py "$NEV" "neither_$LAB" "$LAB" | tee "$BASE_DIR/neither_$LAB.out"

##

LAB="BLINDSR"

pkill python
time /home/psimmerl/mambaforge/envs/ROOT/bin/python -u train_bdt.py "$NEV" "standard_$LAB" "$LAB" | tee "$BASE_DIR/standard_$LAB.out"

pkill python
time /home/psimmerl/mambaforge/envs/ROOT/bin/python -u train_bdt.py "$NEV" "halo_$LAB" "$LAB" | tee "$BASE_DIR/halo_$LAB.out"

pkill python
time /home/psimmerl/mambaforge/envs/ROOT/bin/python -u train_bdt.py "$NEV" "dphi_$LAB" "$LAB" | tee "$BASE_DIR/dphi_$LAB.out"

pkill python
time /home/psimmerl/mambaforge/envs/ROOT/bin/python -u train_bdt.py "$NEV" "neither_$LAB" "$LAB" | tee "$BASE_DIR/neither_$LAB.out"

# ************ #


