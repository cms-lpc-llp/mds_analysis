# CMS Run 3 LLP/MDS Analysis

## Skimming datasets
- skim_csccsc.py
- skim_cscdt.py

## Notebooks
- split_met.ipynb
- --cscdt_cuts.ipynb--
- closure_test.ipynb
- datacard_scan.ipynb

<!-- # Running the Analysis
## Setting Up the Environment

## Cut Flow

## Program Flow

Use skim.py to reduce the large MuonSystem TTree to a flattened TTree with only 1 CSC cluster and only 1 DT cluster. This allows for easy analysis using RDataFrames

Applies the following cuts:
0. Matches gen LLP with reco LLP (MC only)
1. HLT Selections for 1 CSC + 1 DT
2. L1 Trigger on the CSC clusters
3. Requires the CSC clusters to be in-time
4. Requires the DT clusters to be in-time
5. Rejects DT clusters with MB1 hits 
6. Rejects clusters matched to jets
7. Rejects clusters matched to muons
9. Beam halo selection

If there are more than 1 CSC or 1 DT clusters, new rows are appended combinatorially.
If there are fewer than 1 CSC or 1 DT clusters, the event is rejected.

# Data



## MC Event Generation

## Run 3 (2022) Data


# References -->
