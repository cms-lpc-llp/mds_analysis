# CMS Run 3 LLP/MDS Analysis
1 CSC cluster and 1 DT cluster


# Running the Analysis
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


cp mc_hlt569.root mc_hlt569.rdf
cp mc_rdf.root mc_rdf.rdf
cp r3_hlt569.root r3_hlt569.rdf
cp r3_rdf.root r3_rdf.rdf



equipment tracker account
    ask justas


Run 2 0.0018064907 1.0
Inclusivee 0.00894374077896743 4.950892220196501
ME_{T} 0.008357128675141124 4.626167547055642
high 0.010753027789912285 5.952440141582436
omb 0.005734159116954212 3.17419796292162



# Data



## MC Event Generation

## Run 3 (2022) Data


# References
