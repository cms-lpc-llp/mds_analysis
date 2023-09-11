#!/bin/bash

BASE_DIR="/home/psimmerl/LLP/mdc_analysis/reports/weekly/2023-09-07"

if [ $# -ge 1 ]; then
    NEV="$1"
else
    NEV="-1"
fi

if [ $# -ge 2 ]; then
    LAB="_$2"
else
    # LAB="_test"
    LAB=""
fi

# pkill python
# time /home/psimmerl/mambaforge/envs/ROOT/bin/python -u train_bdt.py "$NEV" "$LAB" | tee "$BASE_DIR/$LAB.out"


# All Events

pkill python
time /home/psimmerl/mambaforge/envs/ROOT/bin/python -u "train_bdt copy.py" "$NEV" "standard$LAB"  | tee "$BASE_DIR/standard$LAB.out"

pkill python
time /home/psimmerl/mambaforge/envs/ROOT/bin/python -u "train_bdt copy.py" "$NEV" "neither$LAB"  | tee "$BASE_DIR/neither$LAB.out"

# pkill python
# time /home/psimmerl/mambaforge/envs/ROOT/bin/python -u "train_bdt copy.py" "$NEV" "both$LAB"  | tee "$BASE_DIR/both$LAB.out"

pkill python
time /home/psimmerl/mambaforge/envs/ROOT/bin/python -u "train_bdt copy.py" "$NEV" "halo$LAB"  | tee "$BASE_DIR/halo$LAB.out"

# pkill python
# time /home/psimmerl/mambaforge/envs/ROOT/bin/python -u "train_bdt copy.py" "$NEV" "dphi$LAB"  | tee "$BASE_DIR/dphi$LAB.out"

# # pkill python
# # time /home/psimmerl/mambaforge/envs/ROOT/bin/python -u train_bdt.py "$NEV" "bdt$LAB" | tee "$BASE_DIR/bdt$LAB.out"





