# CMS Run 3 LLP/MDS Analysis

## Skimming datasets
- skim_csccsc.py
- skim_cscdt.py

## Notebooks
- split_met.ipynb
- --cscdt_cuts.ipynb--
- closure_test.ipynb
- datacard_scan.ipynb

## Setting Up the Environment on lxplus
- log in a lxplus8 machine
- Do cmsenv in a CMSSW_14_X release
- Run the .py script: python3 <script.py>


## Program Flow

First, run rdf_hlt_builder.py to load and add the HLTDecision branch to the ntuple.

Use skim_cscdt.py to reduce the large MuonSystem TTree to a flattened TTree with only 1 CSC cluster and only 1 DT cluster. This allows for easy analysis using RDataFrames

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

If there are more than 1 CSC or 1 DT clusters, the largest cluster selected to as the LLP candidate.
