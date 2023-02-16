# Long Lived Particle Analysis @ CMS

contact: Paul Simmerling - psimmerl@caltech.edu

---

## Notes:

- Run3:
  - Simulated signal: ggH_HToSSTobbbb_MH-125_MS-15_ctau-1000_TuneCP5_13TeV-powheg-pythia8_59740pb_weighted.root
  - CMS data: DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi.root

---
## Tasks
- ctau is empty
- Using gLLP_ctau
- dt match gllp decay z is empty
### Week of 1/30/23
- [X] Blind signal data (Extract all clusters with nRechits > 100)
- [X] Just 'raw' signal
- [X] With MET > 50 cut
- [X] With MET > 200 cut
- [X] With event level cuts
  - [X] 0 leptons
  - [X] 0 jets(?)
  - [ ] ~~METNoMu (?)~~
- [ ] With cluster level cuts
  - [X] CBID (CSC only)
  - [X] size > 50
  - [ ] ~~dR(cluster, jet) > 0.4~~
  - [ ] ~~dR(cluster, muon) > 0.4~~
  - [X] no rechits in ME11 & ME12
  - [X] no rechits in MB1 & RB1
  - [X] -5ns < clusterTime < 12.5 ns
  - [X] time spread < 20ns
- [ ] Figure out how to do ctau 95% plot
- [ ] If I have time try BDT again?

Plots:
- [X] MET, dPhi
- [X] Efficiency (sorta)
- [X] nClusters, nCscClusters, nDtClusters, 
- [X] nDtClusters_nCscClusters
- [X] NStation10, AvgStation10
- [X] Size, Eta, Phi
- [ ] ctau 95%
- [ ] Acceptance vs ctau (christina)
- [ ] dPhi_size (ABCD background estimation)

### Week of 2/6/23:

<!-- ### Goals

Signal Study

- Reproduce cut-flow and efficiencies of the standard MET>200 analysis
- Remove MET>200 cut, then check which kinematic effects are changed:
  - Impact on LLP pT, and LLP energy (gen-level)
  - Change in Cluster efficiency
  - Change in ClusterSize shape
  - Change in deltaPhi (MET , cluster) distribution
  - What’s the MET distribution after requiring at least one cluster with clusterSize > 200 (or whatever passes the trigger)
  - What’s the fraction of events with CSC-CSC, DT-CSC, a single CSC cluster only
  - Plot distribution of cluster ID cut variables ( NStation, MaxStation ) -->
