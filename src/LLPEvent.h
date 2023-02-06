//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Dec 14 18:38:55 2022 by ROOT version 6.26/02
// from TTree MuonSystem/MuonSystem
// found on file: ggH_HToSSTobbbb_MH-125_MS-15_ctau-1000_TuneCP5_13TeV-powheg-pythia8_59740pb_weighted.root
//////////////////////////////////////////////////////////

#ifndef LLPEvent_h
#define LLPEvent_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class LLPEvent {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          runNum;
   UInt_t          MC_condition;
   UInt_t          lumiSec;
   UInt_t          evtNum;
   Int_t           mH;
   Int_t           mX;
   Int_t           ctau;
   UInt_t          ZCategory;
   UInt_t          category;
   UInt_t          npv;
   UInt_t          npu;
   Float_t         weight;
   Float_t         higgsPtWeight;
   Float_t         higgsPtWeightSys[9];
   Float_t         scaleWeights[9];
   Float_t         lepOverallSF;
   Float_t         sf_facScaleUp;
   Float_t         sf_facScaleDown;
   Float_t         sf_renScaleUp;
   Float_t         sf_renScaleDown;
   Float_t         sf_facRenScaleUp;
   Float_t         sf_facRenScaleDown;
   Float_t         metSF;
   Float_t         pileupWeight;
   Float_t         pileupWeightUp;
   Float_t         pileupWeightDown;
   Bool_t          Flag_HBHENoiseFilter;
   Bool_t          Flag_BadPFMuonFilter;
   Bool_t          Flag_HBHEIsoNoiseFilter;
   Bool_t          Flag_CSCTightHaloFilter;
   Bool_t          Flag_globalSuperTightHalo2016Filter;
   Bool_t          Flag_goodVertices;
   Bool_t          Flag_ecalBadCalibFilter;
   Bool_t          Flag_BadChargedCandidateFilter;
   Bool_t          Flag_eeBadScFilter;
   Bool_t          Flag_all;
   Bool_t          Flag2_HBHENoiseFilter;
   Bool_t          Flag2_HBHEIsoNoiseFilter;
   Bool_t          Flag2_BadPFMuonFilter;
   Bool_t          Flag2_globalSuperTightHalo2016Filter;
   Bool_t          Flag2_globalTightHalo2016Filter;
   Bool_t          Flag2_BadChargedCandidateFilter;
   Bool_t          Flag2_EcalDeadCellTriggerPrimitiveFilter;
   Bool_t          Flag2_ecalBadCalibFilter;
   Bool_t          Flag2_eeBadScFilter;
   Bool_t          Flag2_all;
   Bool_t          EE_prefiring;
   Float_t         rho;
   Float_t         met;
   Float_t         metNoMu;
   Float_t         metPhi;
   Float_t         metXYCorr;
   Float_t         metPhiXYCorr;
   Float_t         HT;
   Float_t         jetMet_dPhi;
   Float_t         jetMet_dPhiMin;
   Float_t         jetMet_dPhiMin4;
   Float_t         metJESUp;
   Float_t         metJESDown;
   Float_t         metPhiJESUp;
   Float_t         metPhiJESDown;
   Float_t         metJESUpSF;
   Float_t         metJESDownSF;
   Float_t         metEENoise;
   Float_t         metPhiEENoise;
   Float_t         metHEM;
   Float_t         metPhiHEM;
   Float_t         metEENoiseXYCorr;
   Float_t         metPhiEENoiseXYCorr;
   Float_t         metHEMXYCorr;
   Float_t         metPhiHEMXYCorr;
   Float_t         genMetPtTrue;
   Float_t         genMetPhiTrue;
   Float_t         genMetPtCalo;
   Float_t         genMetPhiCalo;
   Int_t           nGenParticle;
   Int_t           gParticleId[1];   //[nGenParticle]
   Int_t           gParticleStatus[1];   //[nGenParticle]
   Int_t           gParticleMotherId[1];   //[nGenParticle]
   Float_t         gParticleE[1];   //[nGenParticle]
   Float_t         gParticlePt[1];   //[nGenParticle]
   Float_t         gParticleEta[1];   //[nGenParticle]
   Float_t         gParticlePhi[1];   //[nGenParticle]
   Int_t           nGenJets;
   Float_t         genJetE[1];   //[nGenJets]
   Float_t         genJetPt[1];   //[nGenJets]
   Float_t         genJetEta[1];   //[nGenJets]
   Float_t         genJetPhi[1];   //[nGenJets]
   Float_t         genJetMET[1];   //[nGenJets]
   Float_t         gWPt;
   Int_t           gLepId;
   Float_t         gLepPt;
   Float_t         gLepE;
   Float_t         gLepEta;
   Float_t         gLepPhi;
   Float_t         gHiggsPt;
   Float_t         gHiggsE;
   Float_t         gHiggsEta;
   Float_t         gHiggsPhi;
   Int_t           nCscRechits;
   Int_t           nEarlyCscRechits;
   Int_t           nLateCscRechits;
   Int_t           nEarly2CscRechits;
   Int_t           nLate2CscRechits;
   Int_t           nCscRings;
   Int_t           nCscPositiveYRechits;
   Int_t           nCscNegativeYRechits;
   Float_t         cscPosTpeak;
   Float_t         cscNegTpeak;
   Int_t           nCscRechitsChamberPlus11;
   Int_t           nCscRechitsChamberPlus12;
   Int_t           nCscRechitsChamberPlus13;
   Int_t           nCscRechitsChamberPlus21;
   Int_t           nCscRechitsChamberPlus22;
   Int_t           nCscRechitsChamberPlus31;
   Int_t           nCscRechitsChamberPlus32;
   Int_t           nCscRechitsChamberPlus41;
   Int_t           nCscRechitsChamberPlus42;
   Int_t           nCscRechitsChamberMinus11;
   Int_t           nCscRechitsChamberMinus12;
   Int_t           nCscRechitsChamberMinus13;
   Int_t           nCscRechitsChamberMinus21;
   Int_t           nCscRechitsChamberMinus22;
   Int_t           nCscRechitsChamberMinus31;
   Int_t           nCscRechitsChamberMinus32;
   Int_t           nCscRechitsChamberMinus41;
   Int_t           nCscRechitsChamberMinus42;
   Int_t           nRpc;
   Int_t           nDtSeg;
   Int_t           nDTRechits;
   Int_t           nDtRings;
   Int_t           nDtWheels25;
   Int_t           nDtStations25;
   Int_t           nDTPositiveYRechits;
   Int_t           nDTNegativeYRechits;
   Int_t           nDTRechitsWheelMinus2;
   Int_t           nDTRechitsWheelMinus1;
   Int_t           nDTRechitsWheel0;
   Int_t           nDTRechitsWheelPlus1;
   Int_t           nDTRechitsWheelPlus2;
   Int_t           nDTRechitsStation1;
   Int_t           nDTRechitsStation2;
   Int_t           nDTRechitsStation3;
   Int_t           nDTRechitsStation4;
   Int_t           nDTRechitsChamberMinus12;
   Int_t           nDTRechitsChamberMinus11;
   Int_t           nDTRechitsChamber10;
   Int_t           nDTRechitsChamberPlus11;
   Int_t           nDTRechitsChamberPlus12;
   Int_t           nDTRechitsChamberMinus22;
   Int_t           nDTRechitsChamberMinus21;
   Int_t           nDTRechitsChamber20;
   Int_t           nDTRechitsChamberPlus21;
   Int_t           nDTRechitsChamberPlus22;
   Int_t           nDTRechitsChamberMinus32;
   Int_t           nDTRechitsChamberMinus31;
   Int_t           nDTRechitsChamber30;
   Int_t           nDTRechitsChamberPlus31;
   Int_t           nDTRechitsChamberPlus32;
   Int_t           nDTRechitsChamberMinus42;
   Int_t           nDTRechitsChamberMinus41;
   Int_t           nDTRechitsChamber40;
   Int_t           nDTRechitsChamberPlus41;
   Int_t           nDTRechitsChamberPlus42;
   Int_t           nDTRechitsSector[4][5][12];
   Int_t           nDTSegSector[4][5][12];
   Int_t           cscRechitsStation[2923];   //[nCscRechits]
   Int_t           cscRechitsChamber[2923];   //[nCscRechits]
   Float_t         cscRechitsPhi[2923];   //[nCscRechits]
   Float_t         cscRechitsEta[2923];   //[nCscRechits]
   Float_t         dtRechitsX[1];   //[nDTRechits]
   Float_t         dtRechitsY[1];   //[nDTRechits]
   Float_t         dtRechitsZ[1];   //[nDTRechits]
   Float_t         dtRechitsEta[1];   //[nDTRechits]
   Float_t         dtRechitsPhi[1];   //[nDTRechits]
   Int_t           dtRechitsStation[1];   //[nDTRechits]
   Int_t           dtRechitsWheel[1];   //[nDTRechits]
   Int_t           dtRechitsClusterId[1];   //[nDTRechits]
   Float_t         rpcEta[1];   //[nRpc]
   Float_t         rpcPhi[1];   //[nRpc]
   Bool_t          rpc_RE12[1];   //[nRpc]
   Bool_t          rpc_RB1[1];   //[nRpc]
   Float_t         dtSegX[270];   //[nDtSeg]
   Float_t         dtSegY[270];   //[nDtSeg]
   Float_t         dtSegZ[270];   //[nDtSeg]
   Float_t         dtSegEta[270];   //[nDtSeg]
   Float_t         dtSegPhi[270];   //[nDtSeg]
   Int_t           dtSegWheel[270];   //[nDtSeg]
   Int_t           dtSegStation[270];   //[nDtSeg]
   Int_t           nCscRechitClusters;
   Bool_t          cscRechitCluster_match_gLLP[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_minDeltaR[5];   //[nCscRechitClusters]
   Int_t           cscRechitCluster_match_gLLP_index[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_eta[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_phi[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_decay_r[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_decay_x[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_decay_y[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_decay_z[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_ctau[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_beta[5];   //[nCscRechitClusters]
   Bool_t          cscRechitCluster_match_gLLP_csc[5];   //[nCscRechitClusters]
   Bool_t          cscRechitCluster_match_gLLP_dt[5];   //[nCscRechitClusters]
   Int_t           cscRechitCluster_match_gLLP_multiplicity[5];   //[nCscRechitClusters]
   Int_t           cscRechitCluster_match_gLLP_EM_multiplicity[5];   //[nCscRechitClusters]
   Bool_t          cscRechitCluster_match_gLLP_daughterKaon[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_e[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_pt[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_EMFracE[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_EMFracEz[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_EMFracP[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_EMFracPz[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_visE[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_visEz[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_visP[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_visPz[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_lepdPhi[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_daughter0_deltaR[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_daughter1_deltaR[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_daughter2_deltaR[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_daughter3_deltaR[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_other_daughter_deltaR[5];   //[nCscRechitClusters]
   Int_t           cscRechitCluster_match_gLLP_other_daughter_index[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_other_eta[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_other_phi[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_other_decay_r[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_other_decay_x[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_other_decay_y[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_other_decay_z[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_other_ctau[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_other_beta[5];   //[nCscRechitClusters]
   Bool_t          cscRechitCluster_match_gLLP_other_csc[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_other_e[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gLLP_other_pt[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterX[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterY[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterZ[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterTime[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterTimeWeighted[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterTimeTotal[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterTimeSpread[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterTimeSpreadWeighted[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterTimeSpreadWeightedAll[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterGenMuonDeltaR[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterXYSpread[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMajorAxis[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMinorAxis[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterEtaPhiSpread[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterPhiSpread[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterEtaSpread[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterDeltaRSpread[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterXSpread[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterRSpread[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterYSpread[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterZSpread[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterPhi[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterEta[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterJetVetoPt[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterJetVetoPtJESDown[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterJetVetoPtJESUp[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterJetVetoEta[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterJetVetoPhi[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterTightJetVetoPt[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterTightJetVetoPtJESDown[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterTightJetVetoPtJESUp[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterTightJetVetoEta[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterTightJetVetoPhi[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterJetVetoE[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterGenJetVetoPt[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterGenJetVetoE[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMuonVetoPt[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMuonVetoE[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMuonVetoPhi[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMuonVetoEta[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterJetVetoElectronEnergyFraction[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterJetVetoPhotonEnergyFraction[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterJetVetoChargedHadronEnergyFraction[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterJetVetoNeutralHadronEnergyFraction[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterJetVetoMuonEnergyFraction[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterJetVetoPt_0p6[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterJetVetoPt_0p8[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterJetVetoE_0p6[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterJetVetoE_0p8[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMuonVetoPt_0p6[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMuonVetoPt_0p8[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMuonVetoE_0p6[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMuonVetoE_0p8[5];   //[nCscRechitClusters]
   Bool_t          cscRechitClusterMuonVetoLooseIso[5];   //[nCscRechitClusters]
   Bool_t          cscRechitClusterMuonVetoTightIso[5];   //[nCscRechitClusters]
   Bool_t          cscRechitClusterMuonVetoVTightIso[5];   //[nCscRechitClusters]
   Bool_t          cscRechitClusterMuonVetoVVTightIso[5];   //[nCscRechitClusters]
   Bool_t          cscRechitClusterMuonVetoTightId[5];   //[nCscRechitClusters]
   Bool_t          cscRechitClusterMuonVetoLooseId[5];   //[nCscRechitClusters]
   Bool_t          cscRechitClusterMuonVetoGlobal[5];   //[nCscRechitClusters]
   Bool_t          cscRechitClusterMuonVetoIso[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterIsoMuonVetoPt[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterIsoMuonVetoE[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterIsoMuonVetoPhi[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterIsoMuonVetoEta[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterGenMuonVetoPt[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterMuonVetoType[5];   //[nCscRechitClusters]
   Bool_t          cscRechitClusterZLep1[5];   //[nCscRechitClusters]
   Bool_t          cscRechitClusterZLep2[5];   //[nCscRechitClusters]
   Bool_t          cscRechitClusterZLep1Tag[5];   //[nCscRechitClusters]
   Bool_t          cscRechitClusterZLep2Tag[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterZLep1Id[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterZLep2Id[5];   //[nCscRechitClusters]
   Bool_t          cscRechitCluster2ZLep1LooseIso[5];   //[nCscRechitClusters]
   Bool_t          cscRechitClusterZLep1TightIso[5];   //[nCscRechitClusters]
   Bool_t          cscRechitClusterZLep1VTightIso[5];   //[nCscRechitClusters]
   Bool_t          cscRechitClusterZLep1VVTightIso[5];   //[nCscRechitClusters]
   Bool_t          cscRechitClusterZLep1TightId[5];   //[nCscRechitClusters]
   Bool_t          cscRechitClusterZLep2LooseIso[5];   //[nCscRechitClusters]
   Bool_t          cscRechitClusterZLep2TightIso[5];   //[nCscRechitClusters]
   Bool_t          cscRechitClusterZLep2VTightIso[5];   //[nCscRechitClusters]
   Bool_t          cscRechitClusterZLep2VVTightIso[5];   //[nCscRechitClusters]
   Bool_t          cscRechitClusterZLep2TightId[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterSize[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMe11Ratio[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMe12Ratio[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNStation[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNStation5[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNStation10[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNStation10perc[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterAvgStation[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterAvgStation5[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterAvgStation10[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterAvgStation10perc[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterMaxStation[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMaxStationRatio[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNChamber[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterMaxChamber[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMaxChamberRatio[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus11[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus12[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus13[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus21[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus22[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus31[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus32[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus41[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus42[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus11[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus12[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus13[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus21[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus22[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus31[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus32[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus41[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus42[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMet_dPhi[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMetJESUp_dPhi[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMetJESDown_dPhi[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMetXYCorr_dPhi[5];   //[nCscRechitClusters]
   Int_t           cscRechitCluster_match_cscRechits_0p4[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNLayersChamberPlus11[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNLayersChamberPlus12[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNLayersChamberPlus13[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNLayersChamberPlus21[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNLayersChamberPlus22[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNLayersChamberPlus31[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNLayersChamberPlus32[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNLayersChamberPlus41[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNLayersChamberPlus42[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNLayersChamberMinus11[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNLayersChamberMinus12[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNLayersChamberMinus13[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNLayersChamberMinus21[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNLayersChamberMinus22[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNLayersChamberMinus31[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNLayersChamberMinus32[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNLayersChamberMinus41[5];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNLayersChamberMinus42[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMetHEM_dPhi[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMetHEMXYCorr_dPhi[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMetEENoise_dPhi[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMetEENoiseXYCorr_dPhi[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMetJesUp_dPhi[5];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMetJesDown_dPhi[5];   //[nCscRechitClusters]
   Int_t           cscRechitCluster_match_dtRechits_phi0p2[5];   //[nCscRechitClusters]
   Int_t           cscRechitCluster_match_dtRechits_0p4[5];   //[nCscRechitClusters]
   Int_t           cscRechitCluster_match_MB1_0p4[5];   //[nCscRechitClusters]
   Int_t           cscRechitCluster_match_dtSeg_0p4[5];   //[nCscRechitClusters]
   Int_t           cscRechitCluster_match_MB1Seg_0p4[5];   //[nCscRechitClusters]
   Int_t           cscRechitCluster_match_RB1_0p4[5];   //[nCscRechitClusters]
   Int_t           cscRechitCluster_match_RE12_0p4[5];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_cluster_dR[5];   //[nCscRechitClusters]
   Int_t           cscRechitCluster_match_cluster_index[5];   //[nCscRechitClusters]
   Int_t           nDtRechitClusters;
   Int_t           nDtRechitClusters2;
   Float_t         dtRechitClusterMaxDPhi[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterMaxDPhi_index[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNSegStation1[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNSegStation2[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNSegStation3[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNSegStation4[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNOppositeSegStation1[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNOppositeSegStation2[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNOppositeSegStation3[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNOppositeSegStation4[4];   //[nDtRechitClusters]
   Bool_t          dtRechitCluster_match_gLLP[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_minDeltaR[4];   //[nDtRechitClusters]
   Int_t           dtRechitCluster_match_gLLP_index[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_eta[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_phi[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_decay_r[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_decay_x[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_decay_y[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_decay_z[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_ctau[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_beta[4];   //[nDtRechitClusters]
   Bool_t          dtRechitCluster_match_gLLP_csc[4];   //[nDtRechitClusters]
   Bool_t          dtRechitCluster_match_gLLP_dt[4];   //[nDtRechitClusters]
   Int_t           dtRechitCluster_match_gLLP_multiplicity[4];   //[nDtRechitClusters]
   Int_t           dtRechitCluster_match_gLLP_EM_multiplicity[4];   //[nDtRechitClusters]
   Bool_t          dtRechitCluster_match_gLLP_daughterKaon[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_e[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_pt[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_EMFracE[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_EMFracEz[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_EMFracP[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_EMFracPz[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_visE[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_visEz[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_visP[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_visPz[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_lepdPhi[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_daughter0_deltaR[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_daughter1_deltaR[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_daughter2_deltaR[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_daughter3_deltaR[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_other_daughter_deltaR[4];   //[nDtRechitClusters]
   Int_t           dtRechitCluster_match_gLLP_other_daughter_index[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_other_eta[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_other_phi[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_other_decay_r[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_other_decay_x[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_other_decay_y[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_other_decay_z[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_other_ctau[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_other_beta[4];   //[nDtRechitClusters]
   Bool_t          dtRechitCluster_match_gLLP_other_csc[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_other_e[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gLLP_other_pt[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterX[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterY[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterZ[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterTime[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterTimeWire[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterTimeWirePruned[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterTimeTotal[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterWheel[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterTimeSpread[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterTimeTotalSpread[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterTimeTotalSpreadPruned[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterTimeWireSpread[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterGenMuonDeltaR[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterXYSpread[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMajorAxis[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMinorAxis[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterEtaPhiSpread[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterPhiSpread[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterEtaSpread[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterDeltaRSpread[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterXSpread[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterRSpread[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterYSpread[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterZSpread[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterPhi[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterEta[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterJetVetoPt[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterJetVetoPtJESUp[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterJetVetoPtJESDown[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterJetVetoEta[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterJetVetoPhi[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterTightJetVetoPt[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterTightJetVetoPtJESUp[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterTightJetVetoPtJESDown[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterTightJetVetoEta[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterTightJetVetoPhi[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterJetVetoE[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterGenJetVetoPt[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterGenJetVetoE[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMuonVetoPt[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMuonVetoE[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMuonVetoPhi[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMuonVetoEta[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterJetVetoElectronEnergyFraction[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterJetVetoPhotonEnergyFraction[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterJetVetoChargedHadronEnergyFraction[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterJetVetoNeutralHadronEnergyFraction[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterJetVetoMuonEnergyFraction[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterJetVetoPt_0p6[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterJetVetoPt_0p8[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterJetVetoE_0p6[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterJetVetoE_0p8[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMuonVetoPt_0p6[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMuonVetoPt_0p8[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMuonVetoE_0p6[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMuonVetoE_0p8[4];   //[nDtRechitClusters]
   Bool_t          dtRechitClusterMuonVetoLooseIso[4];   //[nDtRechitClusters]
   Bool_t          dtRechitClusterMuonVetoTightIso[4];   //[nDtRechitClusters]
   Bool_t          dtRechitClusterMuonVetoVTightIso[4];   //[nDtRechitClusters]
   Bool_t          dtRechitClusterMuonVetoVVTightIso[4];   //[nDtRechitClusters]
   Bool_t          dtRechitClusterMuonVetoTightId[4];   //[nDtRechitClusters]
   Bool_t          dtRechitClusterMuonVetoLooseId[4];   //[nDtRechitClusters]
   Bool_t          dtRechitClusterMuonVetoGlobal[4];   //[nDtRechitClusters]
   Bool_t          dtRechitClusterMuonVetoIso[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterIsoMuonVetoPt[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterIsoMuonVetoE[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterIsoMuonVetoPhi[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterIsoMuonVetoEta[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterGenMuonVetoPt[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterMuonVetoType[4];   //[nDtRechitClusters]
   Bool_t          dtRechitClusterZLep1[4];   //[nDtRechitClusters]
   Bool_t          dtRechitClusterZLep2[4];   //[nDtRechitClusters]
   Bool_t          dtRechitClusterZLep1Tag[4];   //[nDtRechitClusters]
   Bool_t          dtRechitClusterZLep2Tag[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterZLep1Id[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterZLep2Id[4];   //[nDtRechitClusters]
   Bool_t          dtRechitCluster2ZLep1LooseIso[4];   //[nDtRechitClusters]
   Bool_t          dtRechitClusterZLep1TightIso[4];   //[nDtRechitClusters]
   Bool_t          dtRechitClusterZLep1VTightIso[4];   //[nDtRechitClusters]
   Bool_t          dtRechitClusterZLep1VVTightIso[4];   //[nDtRechitClusters]
   Bool_t          dtRechitClusterZLep1TightId[4];   //[nDtRechitClusters]
   Bool_t          dtRechitClusterZLep2LooseIso[4];   //[nDtRechitClusters]
   Bool_t          dtRechitClusterZLep2TightIso[4];   //[nDtRechitClusters]
   Bool_t          dtRechitClusterZLep2VTightIso[4];   //[nDtRechitClusters]
   Bool_t          dtRechitClusterZLep2VVTightIso[4];   //[nDtRechitClusters]
   Bool_t          dtRechitClusterZLep2TightId[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterSize[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNoiseHit[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNoiseHitStation1[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNoiseHitStation2[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNoiseHitStation3[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNoiseHitStation4[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMe11Ratio[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMe12Ratio[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNStation[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNStation5[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNStation10[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNStation10perc[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterAvgStation[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterAvgStation5[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterAvgStation10[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterAvgStation10perc[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterMaxStation[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMaxStationRatio[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNChamber[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterMaxChamber[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMaxChamberRatio[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNSegmentStation1[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNSegmentStation2[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNSegmentStation3[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNSegmentStation4[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNRechitChamberPlus11[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNRechitChamberPlus12[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNRechitChamberPlus13[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNRechitChamberPlus21[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNRechitChamberPlus22[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNRechitChamberPlus31[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNRechitChamberPlus32[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNRechitChamberPlus41[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNRechitChamberPlus42[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNRechitChamberMinus11[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNRechitChamberMinus12[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNRechitChamberMinus13[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNRechitChamberMinus21[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNRechitChamberMinus22[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNRechitChamberMinus31[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNRechitChamberMinus32[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNRechitChamberMinus41[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNRechitChamberMinus42[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMet_dPhi[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMetJESUp_dPhi[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMetJESDown_dPhi[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMetXYCorr_dPhi[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNLayersChamberPlus11[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNLayersChamberPlus12[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNLayersChamberPlus13[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNLayersChamberPlus21[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNLayersChamberPlus22[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNLayersChamberPlus31[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNLayersChamberPlus32[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNLayersChamberPlus41[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNLayersChamberPlus42[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNLayersChamberMinus11[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNLayersChamberMinus12[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNLayersChamberMinus13[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNLayersChamberMinus21[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNLayersChamberMinus22[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNLayersChamberMinus31[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNLayersChamberMinus32[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNLayersChamberMinus41[4];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNLayersChamberMinus42[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMetHEM_dPhi[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMetHEMXYCorr_dPhi[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMetEENoise_dPhi[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMetEENoiseXYCorr_dPhi[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMetJesUp_dPhi[4];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMetJesDown_dPhi[4];   //[nDtRechitClusters]
   Int_t           dtRechitCluster_match_dtSeg_0p5[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_dtSegTime_0p5[4];   //[nDtRechitClusters]
   Int_t           dtRechitCluster_match_dtSeg_0p4[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_dtSegTime_0p4[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_dtSegTimeSpread_0p5[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_dtSegTimeSpread_0p4[4];   //[nDtRechitClusters]
   Int_t           dtRechitCluster_match_dtSeg_sameStation_0p5[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_dtSegTime_sameStation_0p5[4];   //[nDtRechitClusters]
   Int_t           dtRechitCluster_match_dtSeg_sameStation_0p4[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_dtSegTime_sameStation_0p4[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_dtSegTimeSpread_sameStation_0p5[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_dtSegTimeSpread_sameStation_0p4[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_RPCTime_dPhi0p5[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_RPCTimeSpread_dPhi0p5[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_RPCTime_dR0p4[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_RPCTimeSpread_dR0p4[4];   //[nDtRechitClusters]
   Int_t           dtRechitCluster_match_RPChits_dR0p4[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_RPCTime_sameStation_dR0p4[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_RPCTimeSpread_sameStation_dR0p4[4];   //[nDtRechitClusters]
   Int_t           dtRechitCluster_match_RPChits_sameStation_dR0p4[4];   //[nDtRechitClusters]
   Int_t           dtRechitCluster_match_gParticle_Id[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gParticle_Pt[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gParticle_Eta[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gParticle_Phi[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gParticle_E[4];   //[nDtRechitClusters]
   Int_t           dtRechitCluster_match_gParticle_Status[4];   //[nDtRechitClusters]
   Int_t           dtRechitCluster_match_gParticle_MotherId[4];   //[nDtRechitClusters]
   Float_t         dtRechitCluster_match_gParticle_deltaR[4];   //[nDtRechitClusters]
   Int_t           dtRechitCluster_match_RPChits_dPhi0p5[4];   //[nDtRechitClusters]
   Int_t           dtRechitCluster_match_RPCBx_dPhi0p5[4];   //[nDtRechitClusters]
   Int_t           dtRechitCluster_match_RB1_0p4[4];   //[nDtRechitClusters]
   Int_t           dtRechitCluster_match_RB1_dPhi0p5[4];   //[nDtRechitClusters]
   Int_t           dtRechitCluster_match_MB1Seg_0p4[4];   //[nDtRechitClusters]
   Int_t           dtRechitCluster_match_MB1Seg_0p5[4];   //[nDtRechitClusters]
   Int_t           dtRechitCluster_match_MB1hits_0p4[4];   //[nDtRechitClusters]
   Int_t           dtRechitCluster_match_MB1hits_0p5[4];   //[nDtRechitClusters]
   Int_t           dtRechitCluster_match_MB1hits_cosmics_plus[4];   //[nDtRechitClusters]
   Int_t           dtRechitCluster_match_MB1hits_cosmics_minus[4];   //[nDtRechitClusters]
   Bool_t          dtRechitCluster2_match_gLLP[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_minDeltaR[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2_match_gLLP_index[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_eta[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_phi[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_decay_r[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_decay_x[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_decay_y[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_decay_z[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_ctau[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_beta[1];   //[nDtRechitClusters2]
   Bool_t          dtRechitCluster2_match_gLLP_csc[1];   //[nDtRechitClusters2]
   Bool_t          dtRechitCluster2_match_gLLP_dt[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2_match_gLLP_multiplicity[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2_match_gLLP_EM_multiplicity[1];   //[nDtRechitClusters2]
   Bool_t          dtRechitCluster2_match_gLLP_daughterKaon[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_e[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_pt[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_EMFracE[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_EMFracEz[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_EMFracP[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_EMFracPz[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_visE[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_visEz[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_visP[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_visPz[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_lepdPhi[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_daughter0_deltaR[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_daughter1_deltaR[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_daughter2_deltaR[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_daughter3_deltaR[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_other_daughter_deltaR[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2_match_gLLP_other_daughter_index[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_other_eta[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_other_phi[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_other_decay_r[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_other_decay_x[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_other_decay_y[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_other_decay_z[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_other_ctau[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_other_beta[1];   //[nDtRechitClusters2]
   Bool_t          dtRechitCluster2_match_gLLP_other_csc[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_other_e[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gLLP_other_pt[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2X[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2Y[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2Z[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2Time[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2TimeWire[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2TimeWirePruned[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2TimeTotal[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2Wheel[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2TimeSpread[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2TimeTotalSpread[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2TimeTotalSpreadPruned[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2TimeWireSpread[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2GenMuonDeltaR[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2XYSpread[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2MajorAxis[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2MinorAxis[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2EtaPhiSpread[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2PhiSpread[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2EtaSpread[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2DeltaRSpread[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2XSpread[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2RSpread[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2YSpread[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2ZSpread[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2Phi[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2Eta[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2JetVetoPt[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2JetVetoEta[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2JetVetoPhi[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2JetVetoE[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2GenJetVetoPt[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2GenJetVetoE[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2MuonVetoPt[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2MuonVetoE[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2MuonVetoPhi[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2MuonVetoEta[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2JetVetoElectronEnergyFraction[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2JetVetoPhotonEnergyFraction[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2JetVetoChargedHadronEnergyFraction[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2JetVetoNeutralHadronEnergyFraction[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2JetVetoMuonEnergyFraction[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2JetVetoPt_0p6[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2JetVetoPt_0p8[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2JetVetoE_0p6[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2JetVetoE_0p8[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2MuonVetoPt_0p6[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2MuonVetoPt_0p8[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2MuonVetoE_0p6[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2MuonVetoE_0p8[1];   //[nDtRechitClusters2]
   Bool_t          dtRechitCluster2MuonVetoLooseIso[1];   //[nDtRechitClusters2]
   Bool_t          dtRechitCluster2MuonVetoTightIso[1];   //[nDtRechitClusters2]
   Bool_t          dtRechitCluster2MuonVetoVTightIso[1];   //[nDtRechitClusters2]
   Bool_t          dtRechitCluster2MuonVetoVVTightIso[1];   //[nDtRechitClusters2]
   Bool_t          dtRechitCluster2MuonVetoTightId[1];   //[nDtRechitClusters2]
   Bool_t          dtRechitCluster2MuonVetoLooseId[1];   //[nDtRechitClusters2]
   Bool_t          dtRechitCluster2MuonVetoGlobal[1];   //[nDtRechitClusters2]
   Bool_t          dtRechitCluster2MuonVetoIso[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2IsoMuonVetoPt[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2IsoMuonVetoE[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2IsoMuonVetoPhi[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2IsoMuonVetoEta[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2GenMuonVetoPt[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2MuonVetoType[1];   //[nDtRechitClusters2]
   Bool_t          dtRechitCluster2ZLep1[1];   //[nDtRechitClusters2]
   Bool_t          dtRechitCluster2ZLep2[1];   //[nDtRechitClusters2]
   Bool_t          dtRechitCluster2ZLep1Tag[1];   //[nDtRechitClusters2]
   Bool_t          dtRechitCluster2ZLep2Tag[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2ZLep1Id[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2ZLep2Id[1];   //[nDtRechitClusters2]
   Bool_t          dtRechitCluster22ZLep1LooseIso[1];   //[nDtRechitClusters2]
   Bool_t          dtRechitCluster2ZLep1TightIso[1];   //[nDtRechitClusters2]
   Bool_t          dtRechitCluster2ZLep1VTightIso[1];   //[nDtRechitClusters2]
   Bool_t          dtRechitCluster2ZLep1VVTightIso[1];   //[nDtRechitClusters2]
   Bool_t          dtRechitCluster2ZLep1TightId[1];   //[nDtRechitClusters2]
   Bool_t          dtRechitCluster2ZLep2LooseIso[1];   //[nDtRechitClusters2]
   Bool_t          dtRechitCluster2ZLep2TightIso[1];   //[nDtRechitClusters2]
   Bool_t          dtRechitCluster2ZLep2VTightIso[1];   //[nDtRechitClusters2]
   Bool_t          dtRechitCluster2ZLep2VVTightIso[1];   //[nDtRechitClusters2]
   Bool_t          dtRechitCluster2ZLep2TightId[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2Size[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2Me11Ratio[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2Me12Ratio[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NStation[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NStation5[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NStation10[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NStation10perc[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2AvgStation[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2AvgStation5[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2AvgStation10[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2AvgStation10perc[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2MaxStation[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2MaxStationRatio[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NChamber[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2MaxChamber[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2MaxChamberRatio[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NSegmentStation1[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NSegmentStation2[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NSegmentStation3[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NSegmentStation4[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NRechitChamberPlus11[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NRechitChamberPlus12[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NRechitChamberPlus13[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NRechitChamberPlus21[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NRechitChamberPlus22[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NRechitChamberPlus31[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NRechitChamberPlus32[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NRechitChamberPlus41[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NRechitChamberPlus42[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NRechitChamberMinus11[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NRechitChamberMinus12[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NRechitChamberMinus13[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NRechitChamberMinus21[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NRechitChamberMinus22[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NRechitChamberMinus31[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NRechitChamberMinus32[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NRechitChamberMinus41[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NRechitChamberMinus42[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2Met_dPhi[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2MetXYCorr_dPhi[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NLayersChamberPlus11[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NLayersChamberPlus12[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NLayersChamberPlus13[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NLayersChamberPlus21[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NLayersChamberPlus22[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NLayersChamberPlus31[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NLayersChamberPlus32[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NLayersChamberPlus41[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NLayersChamberPlus42[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NLayersChamberMinus11[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NLayersChamberMinus12[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NLayersChamberMinus13[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NLayersChamberMinus21[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NLayersChamberMinus22[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NLayersChamberMinus31[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NLayersChamberMinus32[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NLayersChamberMinus41[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2NLayersChamberMinus42[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2MetHEM_dPhi[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2MetHEMXYCorr_dPhi[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2MetEENoise_dPhi[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2MetEENoiseXYCorr_dPhi[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2MetJesUp_dPhi[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2MetJesDown_dPhi[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2_match_dtSeg_0p5[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_dtSegTime_0p5[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2_match_dtSeg_0p4[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_dtSegTime_0p4[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_dtSegTimeSpread_0p5[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_dtSegTimeSpread_0p4[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2_match_dtSeg_sameStation_0p5[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_dtSegTime_sameStation_0p5[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2_match_dtSeg_sameStation_0p4[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_dtSegTime_sameStation_0p4[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_dtSegTimeSpread_sameStation_0p5[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_dtSegTimeSpread_sameStation_0p4[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_RPCTime_dPhi0p5[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_RPCTimeSpread_dPhi0p5[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_RPCTime_dR0p4[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_RPCTimeSpread_dR0p4[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2_match_RPChits_dR0p4[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_RPCTime_sameStation_dR0p4[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_RPCTimeSpread_sameStation_dR0p4[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2_match_RPChits_sameStation_dR0p4[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2_match_gParticle_Id[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gParticle_Pt[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gParticle_Eta[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gParticle_Phi[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gParticle_E[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2_match_gParticle_Status[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2_match_gParticle_MotherId[1];   //[nDtRechitClusters2]
   Float_t         dtRechitCluster2_match_gParticle_deltaR[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2_match_RPChits_dPhi0p5[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2_match_RPCBx_dPhi0p5[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2_match_RB1_0p4[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2_match_RB1_dPhi0p5[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2_match_MB1Seg_0p4[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2_match_MB1Seg_0p5[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2_match_MB1hits_0p4[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2_match_MB1hits_0p5[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2_match_MB1hits_cosmics_plus[1];   //[nDtRechitClusters2]
   Int_t           dtRechitCluster2_match_MB1hits_cosmics_minus[1];   //[nDtRechitClusters2]
   Int_t           gLLP_multiplicity[2];
   Int_t           gLLP_multiplicity20[2];
   Int_t           gLLP_EM_multiplicity[2];
   Float_t         gLLP_eta[2];
   Float_t         gLLP_phi[2];
   Float_t         gLLP_csc[2];
   Float_t         gLLP_dt[2];
   Float_t         gLLP_beta[2];
   Float_t         gLLP_maxMatchedDis[2];
   Float_t         gLLP_e[2];
   Float_t         gLLP_pt[2];
   Float_t         gLLP_lepdPhi[2];
   Bool_t          gLLP_daughterKaon[2];
   Float_t         gLLP_ctau[2];
   Float_t         gLLP_EMFracE[2];
   Float_t         gLLP_EMFracEz[2];
   Float_t         gLLP_EMFracP[2];
   Float_t         gLLP_EMFracPz[2];
   Float_t         gLLP_visE[2];
   Float_t         gLLP_visE20[2];
   Float_t         gLLP_visEz[2];
   Float_t         gLLP_visP[2];
   Float_t         gLLP_visPz[2];
   Float_t         gLLP_match_jet_pt[2];
   Int_t           gLLP_match_jet_index[2];
   Float_t         gLLP_match_jet_minDeltaR[2];
   Float_t         gLLP_decay_vertex_r[2];
   Float_t         gLLP_decay_vertex_x[2];
   Float_t         gLLP_decay_vertex_y[2];
   Float_t         gLLP_decay_vertex_z[2];
   Float_t         gLLP_deltaR;
   Float_t         gLLP_daughter_deltaR[2];
   Float_t         gLLP_daughter_pt[4];
   Int_t           gLLP_daughter_id[4];
   Float_t         gLLP_daughter_eta[4];
   Float_t         gLLP_daughter_phi[4];
   Float_t         gLLP_daughter_e[4];
   Float_t         gLLP_daughter_mass[4];
   Int_t           nGlobalMuons;
   Float_t         GlobalMuonPt[1];   //[nGlobalMuons]
   Float_t         GlobalMuonEta[1];   //[nGlobalMuons]
   Float_t         GlobalMuonPhi[1];   //[nGlobalMuons]
   Bool_t          GlobalMuonLooseId[1];   //[nGlobalMuons]
   Int_t           nLeptons;
   Float_t         lepE[2];   //[nLeptons]
   Float_t         lepPt[2];   //[nLeptons]
   Float_t         lepEta[2];   //[nLeptons]
   Float_t         lepPhi[2];   //[nLeptons]
   Int_t           lepPdgId[2];   //[nLeptons]
   Float_t         lepDZ[2];   //[nLeptons]
   Bool_t          lepPassId[2];   //[nLeptons]
   Bool_t          lepFromZ[2];   //[nLeptons]
   Float_t         lepEff[2];   //[nLeptons]
   Float_t         lepSF[2];   //[nLeptons]
   Float_t         lepTriggerSF[2];   //[nLeptons]
   Float_t         lepTightIdSF[2];   //[nLeptons]
   Float_t         lepLooseIdSF[2];   //[nLeptons]
   Float_t         lepTightIsoSF[2];   //[nLeptons]
   Float_t         lepLooseIsoSF[2];   //[nLeptons]
   Bool_t          lepTag[2];   //[nLeptons]
   Bool_t          lepPassLooseIso[2];   //[nLeptons]
   Bool_t          lepPassTightIso[2];   //[nLeptons]
   Bool_t          lepPassVTightIso[2];   //[nLeptons]
   Bool_t          lepPassVVTightIso[2];   //[nLeptons]
   Float_t         MT;
   Float_t         ZMass1;
   Float_t         ZMass;
   Float_t         ZPt;
   Float_t         ZEta;
   Float_t         ZPhi;
   Int_t           ZleptonIndex1;
   Int_t           ZleptonIndex2;
   Int_t           nJets;
   Float_t         jetE[11];   //[nJets]
   Float_t         jetPt[11];   //[nJets]
   Float_t         jetEta[11];   //[nJets]
   Float_t         jetPhi[11];   //[nJets]
   Float_t         jetTime[11];   //[nJets]
   Bool_t          jetPassId[11];   //[nJets]
   Float_t         jetPtJESUp[11];   //[nJets]
   Float_t         jetPtJESDown[11];   //[nJets]
   Float_t         jetEJESUp[11];   //[nJets]
   Float_t         jetEJESDown[11];   //[nJets]
   Float_t         JecUnc[11];   //[nJets]
   Float_t         jet_match_llp_pt[11];   //[nJets]
   Int_t           jet_match_llp_index[11];   //[nJets]
   Float_t         jet_match_llp_minDeltaR[11];   //[nJets]
   Float_t         jet_match_genJet_pt[11];   //[nJets]
   Int_t           jet_match_genJet_index[11];   //[nJets]
   Float_t         jet_match_genJet_minDeltaR[11];   //[nJets]
   Bool_t          jetTightPassId[11];   //[nJets]
   Bool_t          HLTDecision[982];
   Bool_t          METTrigger;
   Bool_t          METNoMuTrigger;
   Float_t         jetChargedEMEnergyFraction[11];   //[nJets]
   Float_t         jetChargedHadronEnergyFraction[11];   //[nJets]
   Float_t         jetNeutralEMEnergyFraction[11];   //[nJets]
   Float_t         jetNeutralHadronEnergyFraction[11];   //[nJets]

   // List of branches
   TBranch        *b_runNum;   //!
   TBranch        *b_MC_condition;   //!
   TBranch        *b_lumiSec;   //!
   TBranch        *b_evtNum;   //!
   TBranch        *b_mH;   //!
   TBranch        *b_mX;   //!
   TBranch        *b_ctau;   //!
   TBranch        *b_ZCategory;   //!
   TBranch        *b_category;   //!
   TBranch        *b_npv;   //!
   TBranch        *b_npu;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_higgsPtWeight;   //!
   TBranch        *b_higgsPtWeightSys;   //!
   TBranch        *b_scaleWeights;   //!
   TBranch        *b_lepOverallSF;   //!
   TBranch        *b_sf_facScaleUp;   //!
   TBranch        *b_sf_facScaleDown;   //!
   TBranch        *b_sf_renScaleUp;   //!
   TBranch        *b_sf_renScaleDown;   //!
   TBranch        *b_sf_facRenScaleUp;   //!
   TBranch        *b_sf_facRenScaleDown;   //!
   TBranch        *b_metSF;   //!
   TBranch        *b_pileupWeight;   //!
   TBranch        *b_pileupWeightUp;   //!
   TBranch        *b_pileupWeightDown;   //!
   TBranch        *b_Flag_HBHENoiseFilter;   //!
   TBranch        *b_Flag_BadPFMuonFilter;   //!
   TBranch        *b_Flag_HBHEIsoNoiseFilter;   //!
   TBranch        *b_Flag_CSCTightHaloFilter;   //!
   TBranch        *b_Flag_globalSuperTightHalo2016Filter;   //!
   TBranch        *b_Flag_goodVertices;   //!
   TBranch        *b_Flag_ecalBadCalibFilter;   //!
   TBranch        *b_Flag_BadChargedCandidateFilter;   //!
   TBranch        *b_Flag_eeBadScFilter;   //!
   TBranch        *b_Flag_all;   //!
   TBranch        *b_Flag2_HBHENoiseFilter;   //!
   TBranch        *b_Flag2_HBHEIsoNoiseFilter;   //!
   TBranch        *b_Flag2_BadPFMuonFilter;   //!
   TBranch        *b_Flag2_globalSuperTightHalo2016Filter;   //!
   TBranch        *b_Flag2_globalTightHalo2016Filter;   //!
   TBranch        *b_Flag2_BadChargedCandidateFilter;   //!
   TBranch        *b_Flag2_EcalDeadCellTriggerPrimitiveFilter;   //!
   TBranch        *b_Flag2_ecalBadCalibFilter;   //!
   TBranch        *b_Flag2_eeBadScFilter;   //!
   TBranch        *b_Flag2_all;   //!
   TBranch        *b_EE_prefiring;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_met;   //!
   TBranch        *b_metNoMu;   //!
   TBranch        *b_metPhi;   //!
   TBranch        *b_metXYCorr;   //!
   TBranch        *b_metPhiXYCorr;   //!
   TBranch        *b_HT;   //!
   TBranch        *b_jetMet_dPhi;   //!
   TBranch        *b_jetMet_dPhiMin;   //!
   TBranch        *b_jetMet_dPhiMin4;   //!
   TBranch        *b_metJESUp;   //!
   TBranch        *b_metJESDown;   //!
   TBranch        *b_metPhiJESUp;   //!
   TBranch        *b_metPhiJESDown;   //!
   TBranch        *b_metJESUpSF;   //!
   TBranch        *b_metJESDownSF;   //!
   TBranch        *b_metEENoise;   //!
   TBranch        *b_metPhiEENoise;   //!
   TBranch        *b_metHEM;   //!
   TBranch        *b_metPhiHEM;   //!
   TBranch        *b_metEENoiseXYCorr;   //!
   TBranch        *b_metPhiEENoiseXYCorr;   //!
   TBranch        *b_metHEMXYCorr;   //!
   TBranch        *b_metPhiHEMXYCorr;   //!
   TBranch        *b_genMetPtTrue;   //!
   TBranch        *b_genMetPhiTrue;   //!
   TBranch        *b_genMetPtCalo;   //!
   TBranch        *b_genMetPhiCalo;   //!
   TBranch        *b_nGenParticle;   //!
   TBranch        *b_gParticleId;   //!
   TBranch        *b_gParticleStatus;   //!
   TBranch        *b_gParticleMotherId;   //!
   TBranch        *b_gParticleE;   //!
   TBranch        *b_gParticlePt;   //!
   TBranch        *b_gParticleEta;   //!
   TBranch        *b_gParticlePhi;   //!
   TBranch        *b_nGenJets;   //!
   TBranch        *b_genJetE;   //!
   TBranch        *b_genJetPt;   //!
   TBranch        *b_genJetEta;   //!
   TBranch        *b_genJetPhi;   //!
   TBranch        *b_genJetMET;   //!
   TBranch        *b_gWPt;   //!
   TBranch        *b_gLepId;   //!
   TBranch        *b_gLepPt;   //!
   TBranch        *b_gLepE;   //!
   TBranch        *b_gLepEta;   //!
   TBranch        *b_gLepPhi;   //!
   TBranch        *b_gHiggsPt;   //!
   TBranch        *b_gHiggsE;   //!
   TBranch        *b_gHiggsEta;   //!
   TBranch        *b_gHiggsPhi;   //!
   TBranch        *b_nCscRechits;   //!
   TBranch        *b_nEarlyCscRechits;   //!
   TBranch        *b_nLateCscRechits;   //!
   TBranch        *b_nEarly2CscRechits;   //!
   TBranch        *b_nLate2CscRechits;   //!
   TBranch        *b_nCscRings;   //!
   TBranch        *b_nCscPositiveYRechits;   //!
   TBranch        *b_nCscNegativeYRechits;   //!
   TBranch        *b_cscPosTpeak;   //!
   TBranch        *b_cscNegTpeak;   //!
   TBranch        *b_nCscRechitsChamberPlus11;   //!
   TBranch        *b_nCscRechitsChamberPlus12;   //!
   TBranch        *b_nCscRechitsChamberPlus13;   //!
   TBranch        *b_nCscRechitsChamberPlus21;   //!
   TBranch        *b_nCscRechitsChamberPlus22;   //!
   TBranch        *b_nCscRechitsChamberPlus31;   //!
   TBranch        *b_nCscRechitsChamberPlus32;   //!
   TBranch        *b_nCscRechitsChamberPlus41;   //!
   TBranch        *b_nCscRechitsChamberPlus42;   //!
   TBranch        *b_nCscRechitsChamberMinus11;   //!
   TBranch        *b_nCscRechitsChamberMinus12;   //!
   TBranch        *b_nCscRechitsChamberMinus13;   //!
   TBranch        *b_nCscRechitsChamberMinus21;   //!
   TBranch        *b_nCscRechitsChamberMinus22;   //!
   TBranch        *b_nCscRechitsChamberMinus31;   //!
   TBranch        *b_nCscRechitsChamberMinus32;   //!
   TBranch        *b_nCscRechitsChamberMinus41;   //!
   TBranch        *b_nCscRechitsChamberMinus42;   //!
   TBranch        *b_nRpc;   //!
   TBranch        *b_nDtSeg;   //!
   TBranch        *b_nDTRechits;   //!
   TBranch        *b_nDtRings;   //!
   TBranch        *b_nDtWheels25;   //!
   TBranch        *b_nDtStations25;   //!
   TBranch        *b_nDTPositiveYRechits;   //!
   TBranch        *b_nDTNegativeYRechits;   //!
   TBranch        *b_nDTRechitsWheelMinus2;   //!
   TBranch        *b_nDTRechitsWheelMinus1;   //!
   TBranch        *b_nDTRechitsWheel0;   //!
   TBranch        *b_nDTRechitsWheelPlus1;   //!
   TBranch        *b_nDTRechitsWheelPlus2;   //!
   TBranch        *b_nDTRechitsStation1;   //!
   TBranch        *b_nDTRechitsStation2;   //!
   TBranch        *b_nDTRechitsStation3;   //!
   TBranch        *b_nDTRechitsStation4;   //!
   TBranch        *b_nDTRechitsChamberMinus12;   //!
   TBranch        *b_nDTRechitsChamberMinus11;   //!
   TBranch        *b_nDTRechitsChamber10;   //!
   TBranch        *b_nDTRechitsChamberPlus11;   //!
   TBranch        *b_nDTRechitsChamberPlus12;   //!
   TBranch        *b_nDTRechitsChamberMinus22;   //!
   TBranch        *b_nDTRechitsChamberMinus21;   //!
   TBranch        *b_nDTRechitsChamber20;   //!
   TBranch        *b_nDTRechitsChamberPlus21;   //!
   TBranch        *b_nDTRechitsChamberPlus22;   //!
   TBranch        *b_nDTRechitsChamberMinus32;   //!
   TBranch        *b_nDTRechitsChamberMinus31;   //!
   TBranch        *b_nDTRechitsChamber30;   //!
   TBranch        *b_nDTRechitsChamberPlus31;   //!
   TBranch        *b_nDTRechitsChamberPlus32;   //!
   TBranch        *b_nDTRechitsChamberMinus42;   //!
   TBranch        *b_nDTRechitsChamberMinus41;   //!
   TBranch        *b_nDTRechitsChamber40;   //!
   TBranch        *b_nDTRechitsChamberPlus41;   //!
   TBranch        *b_nDTRechitsChamberPlus42;   //!
   TBranch        *b_nDTRechitsSector;   //!
   TBranch        *b_nDTSegSector;   //!
   TBranch        *b_cscRechitsStation;   //!
   TBranch        *b_cscRechitsChamber;   //!
   TBranch        *b_cscRechitsPhi;   //!
   TBranch        *b_cscRechitsEta;   //!
   TBranch        *b_dtRechitsX;   //!
   TBranch        *b_dtRechitsY;   //!
   TBranch        *b_dtRechitsZ;   //!
   TBranch        *b_dtRechitsEta;   //!
   TBranch        *b_dtRechitsPhi;   //!
   TBranch        *b_dtRechitsStation;   //!
   TBranch        *b_dtRechitsWheel;   //!
   TBranch        *b_dtRechitsClusterId;   //!
   TBranch        *b_rpcEta;   //!
   TBranch        *b_rpcPhi;   //!
   TBranch        *b_rpc_RE12;   //!
   TBranch        *b_rpc_RB1;   //!
   TBranch        *b_dtSegX;   //!
   TBranch        *b_dtSegY;   //!
   TBranch        *b_dtSegZ;   //!
   TBranch        *b_dtSegEta;   //!
   TBranch        *b_dtSegPhi;   //!
   TBranch        *b_dtSegWheel;   //!
   TBranch        *b_dtSegStation;   //!
   TBranch        *b_nCscRechitClusters;   //!
   TBranch        *b_cscRechitCluster_match_gLLP;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_minDeltaR;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_index;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_eta;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_phi;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_decay_r;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_decay_x;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_decay_y;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_decay_z;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_ctau;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_beta;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_csc;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_dt;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_multiplicity;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_EM_multiplicity;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_daughterKaon;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_e;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_pt;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_EMFracE;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_EMFracEz;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_EMFracP;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_EMFracPz;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_visE;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_visEz;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_visP;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_visPz;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_lepdPhi;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_daughter0_deltaR;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_daughter1_deltaR;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_daughter2_deltaR;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_daughter3_deltaR;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_other_daughter_deltaR;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_other_daughter_index;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_other_eta;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_other_phi;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_other_decay_r;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_other_decay_x;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_other_decay_y;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_other_decay_z;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_other_ctau;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_other_beta;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_other_csc;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_other_e;   //!
   TBranch        *b_cscRechitCluster_match_gLLP_other_pt;   //!
   TBranch        *b_cscRechitClusterX;   //!
   TBranch        *b_cscRechitClusterY;   //!
   TBranch        *b_cscRechitClusterZ;   //!
   TBranch        *b_cscRechitClusterTime;   //!
   TBranch        *b_cscRechitClusterTimeWeighted;   //!
   TBranch        *b_cscRechitClusterTimeTotal;   //!
   TBranch        *b_cscRechitClusterTimeSpread;   //!
   TBranch        *b_cscRechitClusterTimeSpreadWeighted;   //!
   TBranch        *b_cscRechitClusterTimeSpreadWeightedAll;   //!
   TBranch        *b_cscRechitClusterGenMuonDeltaR;   //!
   TBranch        *b_cscRechitClusterXYSpread;   //!
   TBranch        *b_cscRechitClusterMajorAxis;   //!
   TBranch        *b_cscRechitClusterMinorAxis;   //!
   TBranch        *b_cscRechitClusterEtaPhiSpread;   //!
   TBranch        *b_cscRechitClusterPhiSpread;   //!
   TBranch        *b_cscRechitClusterEtaSpread;   //!
   TBranch        *b_cscRechitClusterDeltaRSpread;   //!
   TBranch        *b_cscRechitClusterXSpread;   //!
   TBranch        *b_cscRechitClusterRSpread;   //!
   TBranch        *b_cscRechitClusterYSpread;   //!
   TBranch        *b_cscRechitClusterZSpread;   //!
   TBranch        *b_cscRechitClusterPhi;   //!
   TBranch        *b_cscRechitClusterEta;   //!
   TBranch        *b_cscRechitClusterJetVetoPt;   //!
   TBranch        *b_cscRechitClusterJetVetoPtJESDown;   //!
   TBranch        *b_cscRechitClusterJetVetoPtJESUp;   //!
   TBranch        *b_cscRechitClusterJetVetoEta;   //!
   TBranch        *b_cscRechitClusterJetVetoPhi;   //!
   TBranch        *b_cscRechitClusterTightJetVetoPt;   //!
   TBranch        *b_cscRechitClusterTightJetVetoPtJESDown;   //!
   TBranch        *b_cscRechitClusterTightJetVetoPtJESUp;   //!
   TBranch        *b_cscRechitClusterTightJetVetoEta;   //!
   TBranch        *b_cscRechitClusterTightJetVetoPhi;   //!
   TBranch        *b_cscRechitClusterJetVetoE;   //!
   TBranch        *b_cscRechitClusterGenJetVetoPt;   //!
   TBranch        *b_cscRechitClusterGenJetVetoE;   //!
   TBranch        *b_cscRechitClusterMuonVetoPt;   //!
   TBranch        *b_cscRechitClusterMuonVetoE;   //!
   TBranch        *b_cscRechitClusterMuonVetoPhi;   //!
   TBranch        *b_cscRechitClusterMuonVetoEta;   //!
   TBranch        *b_cscRechitClusterJetVetoElectronEnergyFraction;   //!
   TBranch        *b_cscRechitClusterJetVetoPhotonEnergyFraction;   //!
   TBranch        *b_cscRechitClusterJetVetoChargedHadronEnergyFraction;   //!
   TBranch        *b_cscRechitClusterJetVetoNeutralHadronEnergyFraction;   //!
   TBranch        *b_cscRechitClusterJetVetoMuonEnergyFraction;   //!
   TBranch        *b_cscRechitClusterJetVetoPt_0p6;   //!
   TBranch        *b_cscRechitClusterJetVetoPt_0p8;   //!
   TBranch        *b_cscRechitClusterJetVetoE_0p6;   //!
   TBranch        *b_cscRechitClusterJetVetoE_0p8;   //!
   TBranch        *b_cscRechitClusterMuonVetoPt_0p6;   //!
   TBranch        *b_cscRechitClusterMuonVetoPt_0p8;   //!
   TBranch        *b_cscRechitClusterMuonVetoE_0p6;   //!
   TBranch        *b_cscRechitClusterMuonVetoE_0p8;   //!
   TBranch        *b_cscRechitClusterMuonVetoLooseIso;   //!
   TBranch        *b_cscRechitClusterMuonVetoTightIso;   //!
   TBranch        *b_cscRechitClusterMuonVetoVTightIso;   //!
   TBranch        *b_cscRechitClusterMuonVetoVVTightIso;   //!
   TBranch        *b_cscRechitClusterMuonVetoTightId;   //!
   TBranch        *b_cscRechitClusterMuonVetoLooseId;   //!
   TBranch        *b_cscRechitClusterMuonVetoGlobal;   //!
   TBranch        *b_cscRechitClusterMuonVetoIso;   //!
   TBranch        *b_cscRechitClusterIsoMuonVetoPt;   //!
   TBranch        *b_cscRechitClusterIsoMuonVetoE;   //!
   TBranch        *b_cscRechitClusterIsoMuonVetoPhi;   //!
   TBranch        *b_cscRechitClusterIsoMuonVetoEta;   //!
   TBranch        *b_cscRechitClusterGenMuonVetoPt;   //!
   TBranch        *b_cscRechitClusterMuonVetoType;   //!
   TBranch        *b_cscRechitClusterZLep1;   //!
   TBranch        *b_cscRechitClusterZLep2;   //!
   TBranch        *b_cscRechitClusterZLep1Tag;   //!
   TBranch        *b_cscRechitClusterZLep2Tag;   //!
   TBranch        *b_cscRechitClusterZLep1Id;   //!
   TBranch        *b_cscRechitClusterZLep2Id;   //!
   TBranch        *b_cscRechitCluster2ZLep1LooseIso;   //!
   TBranch        *b_cscRechitClusterZLep1TightIso;   //!
   TBranch        *b_cscRechitClusterZLep1VTightIso;   //!
   TBranch        *b_cscRechitClusterZLep1VVTightIso;   //!
   TBranch        *b_cscRechitClusterZLep1TightId;   //!
   TBranch        *b_cscRechitClusterZLep2LooseIso;   //!
   TBranch        *b_cscRechitClusterZLep2TightIso;   //!
   TBranch        *b_cscRechitClusterZLep2VTightIso;   //!
   TBranch        *b_cscRechitClusterZLep2VVTightIso;   //!
   TBranch        *b_cscRechitClusterZLep2TightId;   //!
   TBranch        *b_cscRechitClusterSize;   //!
   TBranch        *b_cscRechitClusterMe11Ratio;   //!
   TBranch        *b_cscRechitClusterMe12Ratio;   //!
   TBranch        *b_cscRechitClusterNStation;   //!
   TBranch        *b_cscRechitClusterNStation5;   //!
   TBranch        *b_cscRechitClusterNStation10;   //!
   TBranch        *b_cscRechitClusterNStation10perc;   //!
   TBranch        *b_cscRechitClusterAvgStation;   //!
   TBranch        *b_cscRechitClusterAvgStation5;   //!
   TBranch        *b_cscRechitClusterAvgStation10;   //!
   TBranch        *b_cscRechitClusterAvgStation10perc;   //!
   TBranch        *b_cscRechitClusterMaxStation;   //!
   TBranch        *b_cscRechitClusterMaxStationRatio;   //!
   TBranch        *b_cscRechitClusterNChamber;   //!
   TBranch        *b_cscRechitClusterMaxChamber;   //!
   TBranch        *b_cscRechitClusterMaxChamberRatio;   //!
   TBranch        *b_cscRechitClusterNRechitChamberPlus11;   //!
   TBranch        *b_cscRechitClusterNRechitChamberPlus12;   //!
   TBranch        *b_cscRechitClusterNRechitChamberPlus13;   //!
   TBranch        *b_cscRechitClusterNRechitChamberPlus21;   //!
   TBranch        *b_cscRechitClusterNRechitChamberPlus22;   //!
   TBranch        *b_cscRechitClusterNRechitChamberPlus31;   //!
   TBranch        *b_cscRechitClusterNRechitChamberPlus32;   //!
   TBranch        *b_cscRechitClusterNRechitChamberPlus41;   //!
   TBranch        *b_cscRechitClusterNRechitChamberPlus42;   //!
   TBranch        *b_cscRechitClusterNRechitChamberMinus11;   //!
   TBranch        *b_cscRechitClusterNRechitChamberMinus12;   //!
   TBranch        *b_cscRechitClusterNRechitChamberMinus13;   //!
   TBranch        *b_cscRechitClusterNRechitChamberMinus21;   //!
   TBranch        *b_cscRechitClusterNRechitChamberMinus22;   //!
   TBranch        *b_cscRechitClusterNRechitChamberMinus31;   //!
   TBranch        *b_cscRechitClusterNRechitChamberMinus32;   //!
   TBranch        *b_cscRechitClusterNRechitChamberMinus41;   //!
   TBranch        *b_cscRechitClusterNRechitChamberMinus42;   //!
   TBranch        *b_cscRechitClusterMet_dPhi;   //!
   TBranch        *b_cscRechitClusterMetJESUp_dPhi;   //!
   TBranch        *b_cscRechitClusterMetJESDown_dPhi;   //!
   TBranch        *b_cscRechitClusterMetXYCorr_dPhi;   //!
   TBranch        *b_cscRechitCluster_match_cscRechits_0p4;   //!
   TBranch        *b_cscRechitClusterNLayersChamberPlus11;   //!
   TBranch        *b_cscRechitClusterNLayersChamberPlus12;   //!
   TBranch        *b_cscRechitClusterNLayersChamberPlus13;   //!
   TBranch        *b_cscRechitClusterNLayersChamberPlus21;   //!
   TBranch        *b_cscRechitClusterNLayersChamberPlus22;   //!
   TBranch        *b_cscRechitClusterNLayersChamberPlus31;   //!
   TBranch        *b_cscRechitClusterNLayersChamberPlus32;   //!
   TBranch        *b_cscRechitClusterNLayersChamberPlus41;   //!
   TBranch        *b_cscRechitClusterNLayersChamberPlus42;   //!
   TBranch        *b_cscRechitClusterNLayersChamberMinus11;   //!
   TBranch        *b_cscRechitClusterNLayersChamberMinus12;   //!
   TBranch        *b_cscRechitClusterNLayersChamberMinus13;   //!
   TBranch        *b_cscRechitClusterNLayersChamberMinus21;   //!
   TBranch        *b_cscRechitClusterNLayersChamberMinus22;   //!
   TBranch        *b_cscRechitClusterNLayersChamberMinus31;   //!
   TBranch        *b_cscRechitClusterNLayersChamberMinus32;   //!
   TBranch        *b_cscRechitClusterNLayersChamberMinus41;   //!
   TBranch        *b_cscRechitClusterNLayersChamberMinus42;   //!
   TBranch        *b_cscRechitClusterMetHEM_dPhi;   //!
   TBranch        *b_cscRechitClusterMetHEMXYCorr_dPhi;   //!
   TBranch        *b_cscRechitClusterMetEENoise_dPhi;   //!
   TBranch        *b_cscRechitClusterMetEENoiseXYCorr_dPhi;   //!
   TBranch        *b_cscRechitClusterMetJesUp_dPhi;   //!
   TBranch        *b_cscRechitClusterMetJesDown_dPhi;   //!
   TBranch        *b_cscRechitCluster_match_dtRechits_phi0p2;   //!
   TBranch        *b_cscRechitCluster_match_dtRechits_0p4;   //!
   TBranch        *b_cscRechitCluster_match_MB1_0p4;   //!
   TBranch        *b_cscRechitCluster_match_dtSeg_0p4;   //!
   TBranch        *b_cscRechitCluster_match_MB1Seg_0p4;   //!
   TBranch        *b_cscRechitCluster_match_RB1_0p4;   //!
   TBranch        *b_cscRechitCluster_match_RE12_0p4;   //!
   TBranch        *b_cscRechitCluster_match_cluster_dR;   //!
   TBranch        *b_cscRechitCluster_match_cluster_index;   //!
   TBranch        *b_nDtRechitClusters;   //!
   TBranch        *b_nDtRechitClusters2;   //!
   TBranch        *b_dtRechitClusterMaxDPhi;   //!
   TBranch        *b_dtRechitClusterMaxDPhi_index;   //!
   TBranch        *b_dtRechitClusterNSegStation1;   //!
   TBranch        *b_dtRechitClusterNSegStation2;   //!
   TBranch        *b_dtRechitClusterNSegStation3;   //!
   TBranch        *b_dtRechitClusterNSegStation4;   //!
   TBranch        *b_dtRechitClusterNOppositeSegStation1;   //!
   TBranch        *b_dtRechitClusterNOppositeSegStation2;   //!
   TBranch        *b_dtRechitClusterNOppositeSegStation3;   //!
   TBranch        *b_dtRechitClusterNOppositeSegStation4;   //!
   TBranch        *b_dtRechitCluster_match_gLLP;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_minDeltaR;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_index;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_eta;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_phi;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_decay_r;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_decay_x;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_decay_y;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_decay_z;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_ctau;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_beta;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_csc;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_dt;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_multiplicity;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_EM_multiplicity;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_daughterKaon;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_e;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_pt;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_EMFracE;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_EMFracEz;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_EMFracP;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_EMFracPz;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_visE;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_visEz;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_visP;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_visPz;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_lepdPhi;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_daughter0_deltaR;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_daughter1_deltaR;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_daughter2_deltaR;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_daughter3_deltaR;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_other_daughter_deltaR;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_other_daughter_index;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_other_eta;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_other_phi;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_other_decay_r;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_other_decay_x;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_other_decay_y;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_other_decay_z;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_other_ctau;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_other_beta;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_other_csc;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_other_e;   //!
   TBranch        *b_dtRechitCluster_match_gLLP_other_pt;   //!
   TBranch        *b_dtRechitClusterX;   //!
   TBranch        *b_dtRechitClusterY;   //!
   TBranch        *b_dtRechitClusterZ;   //!
   TBranch        *b_dtRechitClusterTime;   //!
   TBranch        *b_dtRechitClusterTimeWire;   //!
   TBranch        *b_dtRechitClusterTimeWirePruned;   //!
   TBranch        *b_dtRechitClusterTimeTotal;   //!
   TBranch        *b_dtRechitClusterWheel;   //!
   TBranch        *b_dtRechitClusterTimeSpread;   //!
   TBranch        *b_dtRechitClusterTimeTotalSpread;   //!
   TBranch        *b_dtRechitClusterTimeTotalSpreadPruned;   //!
   TBranch        *b_dtRechitClusterTimeWireSpread;   //!
   TBranch        *b_dtRechitClusterGenMuonDeltaR;   //!
   TBranch        *b_dtRechitClusterXYSpread;   //!
   TBranch        *b_dtRechitClusterMajorAxis;   //!
   TBranch        *b_dtRechitClusterMinorAxis;   //!
   TBranch        *b_dtRechitClusterEtaPhiSpread;   //!
   TBranch        *b_dtRechitClusterPhiSpread;   //!
   TBranch        *b_dtRechitClusterEtaSpread;   //!
   TBranch        *b_dtRechitClusterDeltaRSpread;   //!
   TBranch        *b_dtRechitClusterXSpread;   //!
   TBranch        *b_dtRechitClusterRSpread;   //!
   TBranch        *b_dtRechitClusterYSpread;   //!
   TBranch        *b_dtRechitClusterZSpread;   //!
   TBranch        *b_dtRechitClusterPhi;   //!
   TBranch        *b_dtRechitClusterEta;   //!
   TBranch        *b_dtRechitClusterJetVetoPt;   //!
   TBranch        *b_dtRechitClusterJetVetoPtJESUp;   //!
   TBranch        *b_dtRechitClusterJetVetoPtJESDown;   //!
   TBranch        *b_dtRechitClusterJetVetoEta;   //!
   TBranch        *b_dtRechitClusterJetVetoPhi;   //!
   TBranch        *b_dtRechitClusterTightJetVetoPt;   //!
   TBranch        *b_dtRechitClusterTightJetVetoPtJESUp;   //!
   TBranch        *b_dtRechitClusterTightJetVetoPtJESDown;   //!
   TBranch        *b_dtRechitClusterTightJetVetoEta;   //!
   TBranch        *b_dtRechitClusterTightJetVetoPhi;   //!
   TBranch        *b_dtRechitClusterJetVetoE;   //!
   TBranch        *b_dtRechitClusterGenJetVetoPt;   //!
   TBranch        *b_dtRechitClusterGenJetVetoE;   //!
   TBranch        *b_dtRechitClusterMuonVetoPt;   //!
   TBranch        *b_dtRechitClusterMuonVetoE;   //!
   TBranch        *b_dtRechitClusterMuonVetoPhi;   //!
   TBranch        *b_dtRechitClusterMuonVetoEta;   //!
   TBranch        *b_dtRechitClusterJetVetoElectronEnergyFraction;   //!
   TBranch        *b_dtRechitClusterJetVetoPhotonEnergyFraction;   //!
   TBranch        *b_dtRechitClusterJetVetoChargedHadronEnergyFraction;   //!
   TBranch        *b_dtRechitClusterJetVetoNeutralHadronEnergyFraction;   //!
   TBranch        *b_dtRechitClusterJetVetoMuonEnergyFraction;   //!
   TBranch        *b_dtRechitClusterJetVetoPt_0p6;   //!
   TBranch        *b_dtRechitClusterJetVetoPt_0p8;   //!
   TBranch        *b_dtRechitClusterJetVetoE_0p6;   //!
   TBranch        *b_dtRechitClusterJetVetoE_0p8;   //!
   TBranch        *b_dtRechitClusterMuonVetoPt_0p6;   //!
   TBranch        *b_dtRechitClusterMuonVetoPt_0p8;   //!
   TBranch        *b_dtRechitClusterMuonVetoE_0p6;   //!
   TBranch        *b_dtRechitClusterMuonVetoE_0p8;   //!
   TBranch        *b_dtRechitClusterMuonVetoLooseIso;   //!
   TBranch        *b_dtRechitClusterMuonVetoTightIso;   //!
   TBranch        *b_dtRechitClusterMuonVetoVTightIso;   //!
   TBranch        *b_dtRechitClusterMuonVetoVVTightIso;   //!
   TBranch        *b_dtRechitClusterMuonVetoTightId;   //!
   TBranch        *b_dtRechitClusterMuonVetoLooseId;   //!
   TBranch        *b_dtRechitClusterMuonVetoGlobal;   //!
   TBranch        *b_dtRechitClusterMuonVetoIso;   //!
   TBranch        *b_dtRechitClusterIsoMuonVetoPt;   //!
   TBranch        *b_dtRechitClusterIsoMuonVetoE;   //!
   TBranch        *b_dtRechitClusterIsoMuonVetoPhi;   //!
   TBranch        *b_dtRechitClusterIsoMuonVetoEta;   //!
   TBranch        *b_dtRechitClusterGenMuonVetoPt;   //!
   TBranch        *b_dtRechitClusterMuonVetoType;   //!
   TBranch        *b_dtRechitClusterZLep1;   //!
   TBranch        *b_dtRechitClusterZLep2;   //!
   TBranch        *b_dtRechitClusterZLep1Tag;   //!
   TBranch        *b_dtRechitClusterZLep2Tag;   //!
   TBranch        *b_dtRechitClusterZLep1Id;   //!
   TBranch        *b_dtRechitClusterZLep2Id;   //!
   TBranch        *b_dtRechitCluster2ZLep1LooseIso;   //!
   TBranch        *b_dtRechitClusterZLep1TightIso;   //!
   TBranch        *b_dtRechitClusterZLep1VTightIso;   //!
   TBranch        *b_dtRechitClusterZLep1VVTightIso;   //!
   TBranch        *b_dtRechitClusterZLep1TightId;   //!
   TBranch        *b_dtRechitClusterZLep2LooseIso;   //!
   TBranch        *b_dtRechitClusterZLep2TightIso;   //!
   TBranch        *b_dtRechitClusterZLep2VTightIso;   //!
   TBranch        *b_dtRechitClusterZLep2VVTightIso;   //!
   TBranch        *b_dtRechitClusterZLep2TightId;   //!
   TBranch        *b_dtRechitClusterSize;   //!
   TBranch        *b_dtRechitClusterNoiseHit;   //!
   TBranch        *b_dtRechitClusterNoiseHitStation1;   //!
   TBranch        *b_dtRechitClusterNoiseHitStation2;   //!
   TBranch        *b_dtRechitClusterNoiseHitStation3;   //!
   TBranch        *b_dtRechitClusterNoiseHitStation4;   //!
   TBranch        *b_dtRechitClusterMe11Ratio;   //!
   TBranch        *b_dtRechitClusterMe12Ratio;   //!
   TBranch        *b_dtRechitClusterNStation;   //!
   TBranch        *b_dtRechitClusterNStation5;   //!
   TBranch        *b_dtRechitClusterNStation10;   //!
   TBranch        *b_dtRechitClusterNStation10perc;   //!
   TBranch        *b_dtRechitClusterAvgStation;   //!
   TBranch        *b_dtRechitClusterAvgStation5;   //!
   TBranch        *b_dtRechitClusterAvgStation10;   //!
   TBranch        *b_dtRechitClusterAvgStation10perc;   //!
   TBranch        *b_dtRechitClusterMaxStation;   //!
   TBranch        *b_dtRechitClusterMaxStationRatio;   //!
   TBranch        *b_dtRechitClusterNChamber;   //!
   TBranch        *b_dtRechitClusterMaxChamber;   //!
   TBranch        *b_dtRechitClusterMaxChamberRatio;   //!
   TBranch        *b_dtRechitClusterNSegmentStation1;   //!
   TBranch        *b_dtRechitClusterNSegmentStation2;   //!
   TBranch        *b_dtRechitClusterNSegmentStation3;   //!
   TBranch        *b_dtRechitClusterNSegmentStation4;   //!
   TBranch        *b_dtRechitClusterNRechitChamberPlus11;   //!
   TBranch        *b_dtRechitClusterNRechitChamberPlus12;   //!
   TBranch        *b_dtRechitClusterNRechitChamberPlus13;   //!
   TBranch        *b_dtRechitClusterNRechitChamberPlus21;   //!
   TBranch        *b_dtRechitClusterNRechitChamberPlus22;   //!
   TBranch        *b_dtRechitClusterNRechitChamberPlus31;   //!
   TBranch        *b_dtRechitClusterNRechitChamberPlus32;   //!
   TBranch        *b_dtRechitClusterNRechitChamberPlus41;   //!
   TBranch        *b_dtRechitClusterNRechitChamberPlus42;   //!
   TBranch        *b_dtRechitClusterNRechitChamberMinus11;   //!
   TBranch        *b_dtRechitClusterNRechitChamberMinus12;   //!
   TBranch        *b_dtRechitClusterNRechitChamberMinus13;   //!
   TBranch        *b_dtRechitClusterNRechitChamberMinus21;   //!
   TBranch        *b_dtRechitClusterNRechitChamberMinus22;   //!
   TBranch        *b_dtRechitClusterNRechitChamberMinus31;   //!
   TBranch        *b_dtRechitClusterNRechitChamberMinus32;   //!
   TBranch        *b_dtRechitClusterNRechitChamberMinus41;   //!
   TBranch        *b_dtRechitClusterNRechitChamberMinus42;   //!
   TBranch        *b_dtRechitClusterMet_dPhi;   //!
   TBranch        *b_dtRechitClusterMetJESUp_dPhi;   //!
   TBranch        *b_dtRechitClusterMetJESDown_dPhi;   //!
   TBranch        *b_dtRechitClusterMetXYCorr_dPhi;   //!
   TBranch        *b_dtRechitClusterNLayersChamberPlus11;   //!
   TBranch        *b_dtRechitClusterNLayersChamberPlus12;   //!
   TBranch        *b_dtRechitClusterNLayersChamberPlus13;   //!
   TBranch        *b_dtRechitClusterNLayersChamberPlus21;   //!
   TBranch        *b_dtRechitClusterNLayersChamberPlus22;   //!
   TBranch        *b_dtRechitClusterNLayersChamberPlus31;   //!
   TBranch        *b_dtRechitClusterNLayersChamberPlus32;   //!
   TBranch        *b_dtRechitClusterNLayersChamberPlus41;   //!
   TBranch        *b_dtRechitClusterNLayersChamberPlus42;   //!
   TBranch        *b_dtRechitClusterNLayersChamberMinus11;   //!
   TBranch        *b_dtRechitClusterNLayersChamberMinus12;   //!
   TBranch        *b_dtRechitClusterNLayersChamberMinus13;   //!
   TBranch        *b_dtRechitClusterNLayersChamberMinus21;   //!
   TBranch        *b_dtRechitClusterNLayersChamberMinus22;   //!
   TBranch        *b_dtRechitClusterNLayersChamberMinus31;   //!
   TBranch        *b_dtRechitClusterNLayersChamberMinus32;   //!
   TBranch        *b_dtRechitClusterNLayersChamberMinus41;   //!
   TBranch        *b_dtRechitClusterNLayersChamberMinus42;   //!
   TBranch        *b_dtRechitClusterMetHEM_dPhi;   //!
   TBranch        *b_dtRechitClusterMetHEMXYCorr_dPhi;   //!
   TBranch        *b_dtRechitClusterMetEENoise_dPhi;   //!
   TBranch        *b_dtRechitClusterMetEENoiseXYCorr_dPhi;   //!
   TBranch        *b_dtRechitClusterMetJesUp_dPhi;   //!
   TBranch        *b_dtRechitClusterMetJesDown_dPhi;   //!
   TBranch        *b_dtRechitCluster_match_dtSeg_0p5;   //!
   TBranch        *b_dtRechitCluster_match_dtSegTime_0p5;   //!
   TBranch        *b_dtRechitCluster_match_dtSeg_0p4;   //!
   TBranch        *b_dtRechitCluster_match_dtSegTime_0p4;   //!
   TBranch        *b_dtRechitCluster_match_dtSegTimeSpread_0p5;   //!
   TBranch        *b_dtRechitCluster_match_dtSegTimeSpread_0p4;   //!
   TBranch        *b_dtRechitCluster_match_dtSeg_sameStation_0p5;   //!
   TBranch        *b_dtRechitCluster_match_dtSegTime_sameStation_0p5;   //!
   TBranch        *b_dtRechitCluster_match_dtSeg_sameStation_0p4;   //!
   TBranch        *b_dtRechitCluster_match_dtSegTime_sameStation_0p4;   //!
   TBranch        *b_dtRechitCluster_match_dtSegTimeSpread_sameStation_0p5;   //!
   TBranch        *b_dtRechitCluster_match_dtSegTimeSpread_sameStation_0p4;   //!
   TBranch        *b_dtRechitCluster_match_RPCTime_dPhi0p5;   //!
   TBranch        *b_dtRechitCluster_match_RPCTimeSpread_dPhi0p5;   //!
   TBranch        *b_dtRechitCluster_match_RPCTime_dR0p4;   //!
   TBranch        *b_dtRechitCluster_match_RPCTimeSpread_dR0p4;   //!
   TBranch        *b_dtRechitCluster_match_RPChits_dR0p4;   //!
   TBranch        *b_dtRechitCluster_match_RPCTime_sameStation_dR0p4;   //!
   TBranch        *b_dtRechitCluster_match_RPCTimeSpread_sameStation_dR0p4;   //!
   TBranch        *b_dtRechitCluster_match_RPChits_sameStation_dR0p4;   //!
   TBranch        *b_dtRechitCluster_match_gParticle_Id;   //!
   TBranch        *b_dtRechitCluster_match_gParticle_Pt;   //!
   TBranch        *b_dtRechitCluster_match_gParticle_Eta;   //!
   TBranch        *b_dtRechitCluster_match_gParticle_Phi;   //!
   TBranch        *b_dtRechitCluster_match_gParticle_E;   //!
   TBranch        *b_dtRechitCluster_match_gParticle_Status;   //!
   TBranch        *b_dtRechitCluster_match_gParticle_MotherId;   //!
   TBranch        *b_dtRechitCluster_match_gParticle_deltaR;   //!
   TBranch        *b_dtRechitCluster_match_RPChits_dPhi0p5;   //!
   TBranch        *b_dtRechitCluster_match_RPCBx_dPhi0p5;   //!
   TBranch        *b_dtRechitCluster_match_RB1_0p4;   //!
   TBranch        *b_dtRechitCluster_match_RB1_dPhi0p5;   //!
   TBranch        *b_dtRechitCluster_match_MB1Seg_0p4;   //!
   TBranch        *b_dtRechitCluster_match_MB1Seg_0p5;   //!
   TBranch        *b_dtRechitCluster_match_MB1hits_0p4;   //!
   TBranch        *b_dtRechitCluster_match_MB1hits_0p5;   //!
   TBranch        *b_dtRechitCluster_match_MB1hits_cosmics_plus;   //!
   TBranch        *b_dtRechitCluster_match_MB1hits_cosmics_minus;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_minDeltaR;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_index;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_eta;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_phi;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_decay_r;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_decay_x;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_decay_y;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_decay_z;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_ctau;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_beta;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_csc;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_dt;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_multiplicity;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_EM_multiplicity;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_daughterKaon;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_e;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_pt;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_EMFracE;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_EMFracEz;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_EMFracP;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_EMFracPz;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_visE;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_visEz;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_visP;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_visPz;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_lepdPhi;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_daughter0_deltaR;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_daughter1_deltaR;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_daughter2_deltaR;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_daughter3_deltaR;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_other_daughter_deltaR;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_other_daughter_index;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_other_eta;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_other_phi;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_other_decay_r;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_other_decay_x;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_other_decay_y;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_other_decay_z;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_other_ctau;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_other_beta;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_other_csc;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_other_e;   //!
   TBranch        *b_dtRechitCluster2_match_gLLP_other_pt;   //!
   TBranch        *b_dtRechitCluster2X;   //!
   TBranch        *b_dtRechitCluster2Y;   //!
   TBranch        *b_dtRechitCluster2Z;   //!
   TBranch        *b_dtRechitCluster2Time;   //!
   TBranch        *b_dtRechitCluster2TimeWire;   //!
   TBranch        *b_dtRechitCluster2TimeWirePruned;   //!
   TBranch        *b_dtRechitCluster2TimeTotal;   //!
   TBranch        *b_dtRechitCluster2Wheel;   //!
   TBranch        *b_dtRechitCluster2TimeSpread;   //!
   TBranch        *b_dtRechitCluster2TimeTotalSpread;   //!
   TBranch        *b_dtRechitCluster2TimeTotalSpreadPruned;   //!
   TBranch        *b_dtRechitCluster2TimeWireSpread;   //!
   TBranch        *b_dtRechitCluster2GenMuonDeltaR;   //!
   TBranch        *b_dtRechitCluster2XYSpread;   //!
   TBranch        *b_dtRechitCluster2MajorAxis;   //!
   TBranch        *b_dtRechitCluster2MinorAxis;   //!
   TBranch        *b_dtRechitCluster2EtaPhiSpread;   //!
   TBranch        *b_dtRechitCluster2PhiSpread;   //!
   TBranch        *b_dtRechitCluster2EtaSpread;   //!
   TBranch        *b_dtRechitCluster2DeltaRSpread;   //!
   TBranch        *b_dtRechitCluster2XSpread;   //!
   TBranch        *b_dtRechitCluster2RSpread;   //!
   TBranch        *b_dtRechitCluster2YSpread;   //!
   TBranch        *b_dtRechitCluster2ZSpread;   //!
   TBranch        *b_dtRechitCluster2Phi;   //!
   TBranch        *b_dtRechitCluster2Eta;   //!
   TBranch        *b_dtRechitCluster2JetVetoPt;   //!
   TBranch        *b_dtRechitCluster2JetVetoEta;   //!
   TBranch        *b_dtRechitCluster2JetVetoPhi;   //!
   TBranch        *b_dtRechitCluster2JetVetoE;   //!
   TBranch        *b_dtRechitCluster2GenJetVetoPt;   //!
   TBranch        *b_dtRechitCluster2GenJetVetoE;   //!
   TBranch        *b_dtRechitCluster2MuonVetoPt;   //!
   TBranch        *b_dtRechitCluster2MuonVetoE;   //!
   TBranch        *b_dtRechitCluster2MuonVetoPhi;   //!
   TBranch        *b_dtRechitCluster2MuonVetoEta;   //!
   TBranch        *b_dtRechitCluster2JetVetoElectronEnergyFraction;   //!
   TBranch        *b_dtRechitCluster2JetVetoPhotonEnergyFraction;   //!
   TBranch        *b_dtRechitCluster2JetVetoChargedHadronEnergyFraction;   //!
   TBranch        *b_dtRechitCluster2JetVetoNeutralHadronEnergyFraction;   //!
   TBranch        *b_dtRechitCluster2JetVetoMuonEnergyFraction;   //!
   TBranch        *b_dtRechitCluster2JetVetoPt_0p6;   //!
   TBranch        *b_dtRechitCluster2JetVetoPt_0p8;   //!
   TBranch        *b_dtRechitCluster2JetVetoE_0p6;   //!
   TBranch        *b_dtRechitCluster2JetVetoE_0p8;   //!
   TBranch        *b_dtRechitCluster2MuonVetoPt_0p6;   //!
   TBranch        *b_dtRechitCluster2MuonVetoPt_0p8;   //!
   TBranch        *b_dtRechitCluster2MuonVetoE_0p6;   //!
   TBranch        *b_dtRechitCluster2MuonVetoE_0p8;   //!
   TBranch        *b_dtRechitCluster2MuonVetoLooseIso;   //!
   TBranch        *b_dtRechitCluster2MuonVetoTightIso;   //!
   TBranch        *b_dtRechitCluster2MuonVetoVTightIso;   //!
   TBranch        *b_dtRechitCluster2MuonVetoVVTightIso;   //!
   TBranch        *b_dtRechitCluster2MuonVetoTightId;   //!
   TBranch        *b_dtRechitCluster2MuonVetoLooseId;   //!
   TBranch        *b_dtRechitCluster2MuonVetoGlobal;   //!
   TBranch        *b_dtRechitCluster2MuonVetoIso;   //!
   TBranch        *b_dtRechitCluster2IsoMuonVetoPt;   //!
   TBranch        *b_dtRechitCluster2IsoMuonVetoE;   //!
   TBranch        *b_dtRechitCluster2IsoMuonVetoPhi;   //!
   TBranch        *b_dtRechitCluster2IsoMuonVetoEta;   //!
   TBranch        *b_dtRechitCluster2GenMuonVetoPt;   //!
   TBranch        *b_dtRechitCluster2MuonVetoType;   //!
   TBranch        *b_dtRechitCluster2ZLep1;   //!
   TBranch        *b_dtRechitCluster2ZLep2;   //!
   TBranch        *b_dtRechitCluster2ZLep1Tag;   //!
   TBranch        *b_dtRechitCluster2ZLep2Tag;   //!
   TBranch        *b_dtRechitCluster2ZLep1Id;   //!
   TBranch        *b_dtRechitCluster2ZLep2Id;   //!
   TBranch        *b_dtRechitCluster22ZLep1LooseIso;   //!
   TBranch        *b_dtRechitCluster2ZLep1TightIso;   //!
   TBranch        *b_dtRechitCluster2ZLep1VTightIso;   //!
   TBranch        *b_dtRechitCluster2ZLep1VVTightIso;   //!
   TBranch        *b_dtRechitCluster2ZLep1TightId;   //!
   TBranch        *b_dtRechitCluster2ZLep2LooseIso;   //!
   TBranch        *b_dtRechitCluster2ZLep2TightIso;   //!
   TBranch        *b_dtRechitCluster2ZLep2VTightIso;   //!
   TBranch        *b_dtRechitCluster2ZLep2VVTightIso;   //!
   TBranch        *b_dtRechitCluster2ZLep2TightId;   //!
   TBranch        *b_dtRechitCluster2Size;   //!
   TBranch        *b_dtRechitCluster2Me11Ratio;   //!
   TBranch        *b_dtRechitCluster2Me12Ratio;   //!
   TBranch        *b_dtRechitCluster2NStation;   //!
   TBranch        *b_dtRechitCluster2NStation5;   //!
   TBranch        *b_dtRechitCluster2NStation10;   //!
   TBranch        *b_dtRechitCluster2NStation10perc;   //!
   TBranch        *b_dtRechitCluster2AvgStation;   //!
   TBranch        *b_dtRechitCluster2AvgStation5;   //!
   TBranch        *b_dtRechitCluster2AvgStation10;   //!
   TBranch        *b_dtRechitCluster2AvgStation10perc;   //!
   TBranch        *b_dtRechitCluster2MaxStation;   //!
   TBranch        *b_dtRechitCluster2MaxStationRatio;   //!
   TBranch        *b_dtRechitCluster2NChamber;   //!
   TBranch        *b_dtRechitCluster2MaxChamber;   //!
   TBranch        *b_dtRechitCluster2MaxChamberRatio;   //!
   TBranch        *b_dtRechitCluster2NSegmentStation1;   //!
   TBranch        *b_dtRechitCluster2NSegmentStation2;   //!
   TBranch        *b_dtRechitCluster2NSegmentStation3;   //!
   TBranch        *b_dtRechitCluster2NSegmentStation4;   //!
   TBranch        *b_dtRechitCluster2NRechitChamberPlus11;   //!
   TBranch        *b_dtRechitCluster2NRechitChamberPlus12;   //!
   TBranch        *b_dtRechitCluster2NRechitChamberPlus13;   //!
   TBranch        *b_dtRechitCluster2NRechitChamberPlus21;   //!
   TBranch        *b_dtRechitCluster2NRechitChamberPlus22;   //!
   TBranch        *b_dtRechitCluster2NRechitChamberPlus31;   //!
   TBranch        *b_dtRechitCluster2NRechitChamberPlus32;   //!
   TBranch        *b_dtRechitCluster2NRechitChamberPlus41;   //!
   TBranch        *b_dtRechitCluster2NRechitChamberPlus42;   //!
   TBranch        *b_dtRechitCluster2NRechitChamberMinus11;   //!
   TBranch        *b_dtRechitCluster2NRechitChamberMinus12;   //!
   TBranch        *b_dtRechitCluster2NRechitChamberMinus13;   //!
   TBranch        *b_dtRechitCluster2NRechitChamberMinus21;   //!
   TBranch        *b_dtRechitCluster2NRechitChamberMinus22;   //!
   TBranch        *b_dtRechitCluster2NRechitChamberMinus31;   //!
   TBranch        *b_dtRechitCluster2NRechitChamberMinus32;   //!
   TBranch        *b_dtRechitCluster2NRechitChamberMinus41;   //!
   TBranch        *b_dtRechitCluster2NRechitChamberMinus42;   //!
   TBranch        *b_dtRechitCluster2Met_dPhi;   //!
   TBranch        *b_dtRechitCluster2MetXYCorr_dPhi;   //!
   TBranch        *b_dtRechitCluster2NLayersChamberPlus11;   //!
   TBranch        *b_dtRechitCluster2NLayersChamberPlus12;   //!
   TBranch        *b_dtRechitCluster2NLayersChamberPlus13;   //!
   TBranch        *b_dtRechitCluster2NLayersChamberPlus21;   //!
   TBranch        *b_dtRechitCluster2NLayersChamberPlus22;   //!
   TBranch        *b_dtRechitCluster2NLayersChamberPlus31;   //!
   TBranch        *b_dtRechitCluster2NLayersChamberPlus32;   //!
   TBranch        *b_dtRechitCluster2NLayersChamberPlus41;   //!
   TBranch        *b_dtRechitCluster2NLayersChamberPlus42;   //!
   TBranch        *b_dtRechitCluster2NLayersChamberMinus11;   //!
   TBranch        *b_dtRechitCluster2NLayersChamberMinus12;   //!
   TBranch        *b_dtRechitCluster2NLayersChamberMinus13;   //!
   TBranch        *b_dtRechitCluster2NLayersChamberMinus21;   //!
   TBranch        *b_dtRechitCluster2NLayersChamberMinus22;   //!
   TBranch        *b_dtRechitCluster2NLayersChamberMinus31;   //!
   TBranch        *b_dtRechitCluster2NLayersChamberMinus32;   //!
   TBranch        *b_dtRechitCluster2NLayersChamberMinus41;   //!
   TBranch        *b_dtRechitCluster2NLayersChamberMinus42;   //!
   TBranch        *b_dtRechitCluster2MetHEM_dPhi;   //!
   TBranch        *b_dtRechitCluster2MetHEMXYCorr_dPhi;   //!
   TBranch        *b_dtRechitCluster2MetEENoise_dPhi;   //!
   TBranch        *b_dtRechitCluster2MetEENoiseXYCorr_dPhi;   //!
   TBranch        *b_dtRechitCluster2MetJesUp_dPhi;   //!
   TBranch        *b_dtRechitCluster2MetJesDown_dPhi;   //!
   TBranch        *b_dtRechitCluster2_match_dtSeg_0p5;   //!
   TBranch        *b_dtRechitCluster2_match_dtSegTime_0p5;   //!
   TBranch        *b_dtRechitCluster2_match_dtSeg_0p4;   //!
   TBranch        *b_dtRechitCluster2_match_dtSegTime_0p4;   //!
   TBranch        *b_dtRechitCluster2_match_dtSegTimeSpread_0p5;   //!
   TBranch        *b_dtRechitCluster2_match_dtSegTimeSpread_0p4;   //!
   TBranch        *b_dtRechitCluster2_match_dtSeg_sameStation_0p5;   //!
   TBranch        *b_dtRechitCluster2_match_dtSegTime_sameStation_0p5;   //!
   TBranch        *b_dtRechitCluster2_match_dtSeg_sameStation_0p4;   //!
   TBranch        *b_dtRechitCluster2_match_dtSegTime_sameStation_0p4;   //!
   TBranch        *b_dtRechitCluster2_match_dtSegTimeSpread_sameStation_0p5;   //!
   TBranch        *b_dtRechitCluster2_match_dtSegTimeSpread_sameStation_0p4;   //!
   TBranch        *b_dtRechitCluster2_match_RPCTime_dPhi0p5;   //!
   TBranch        *b_dtRechitCluster2_match_RPCTimeSpread_dPhi0p5;   //!
   TBranch        *b_dtRechitCluster2_match_RPCTime_dR0p4;   //!
   TBranch        *b_dtRechitCluster2_match_RPCTimeSpread_dR0p4;   //!
   TBranch        *b_dtRechitCluster2_match_RPChits_dR0p4;   //!
   TBranch        *b_dtRechitCluster2_match_RPCTime_sameStation_dR0p4;   //!
   TBranch        *b_dtRechitCluster2_match_RPCTimeSpread_sameStation_dR0p4;   //!
   TBranch        *b_dtRechitCluster2_match_RPChits_sameStation_dR0p4;   //!
   TBranch        *b_dtRechitCluster2_match_gParticle_Id;   //!
   TBranch        *b_dtRechitCluster2_match_gParticle_Pt;   //!
   TBranch        *b_dtRechitCluster2_match_gParticle_Eta;   //!
   TBranch        *b_dtRechitCluster2_match_gParticle_Phi;   //!
   TBranch        *b_dtRechitCluster2_match_gParticle_E;   //!
   TBranch        *b_dtRechitCluster2_match_gParticle_Status;   //!
   TBranch        *b_dtRechitCluster2_match_gParticle_MotherId;   //!
   TBranch        *b_dtRechitCluster2_match_gParticle_deltaR;   //!
   TBranch        *b_dtRechitCluster2_match_RPChits_dPhi0p5;   //!
   TBranch        *b_dtRechitCluster2_match_RPCBx_dPhi0p5;   //!
   TBranch        *b_dtRechitCluster2_match_RB1_0p4;   //!
   TBranch        *b_dtRechitCluster2_match_RB1_dPhi0p5;   //!
   TBranch        *b_dtRechitCluster2_match_MB1Seg_0p4;   //!
   TBranch        *b_dtRechitCluster2_match_MB1Seg_0p5;   //!
   TBranch        *b_dtRechitCluster2_match_MB1hits_0p4;   //!
   TBranch        *b_dtRechitCluster2_match_MB1hits_0p5;   //!
   TBranch        *b_dtRechitCluster2_match_MB1hits_cosmics_plus;   //!
   TBranch        *b_dtRechitCluster2_match_MB1hits_cosmics_minus;   //!
   TBranch        *b_gLLP_multiplicity;   //!
   TBranch        *b_gLLP_multiplicity20;   //!
   TBranch        *b_gLLP_EM_multiplicity;   //!
   TBranch        *b_gLLP_eta;   //!
   TBranch        *b_gLLP_phi;   //!
   TBranch        *b_gLLP_csc;   //!
   TBranch        *b_gLLP_dt;   //!
   TBranch        *b_gLLP_beta;   //!
   TBranch        *b_gLLP_maxMatchedDis;   //!
   TBranch        *b_gLLP_e;   //!
   TBranch        *b_gLLP_pt;   //!
   TBranch        *b_gLLP_lepdPhi;   //!
   TBranch        *b_gLLP_daughterKaon;   //!
   TBranch        *b_gLLP_ctau;   //!
   TBranch        *b_gLLP_EMFracE;   //!
   TBranch        *b_gLLP_EMFracEz;   //!
   TBranch        *b_gLLP_EMFracP;   //!
   TBranch        *b_gLLP_EMFracPz;   //!
   TBranch        *b_gLLP_visE;   //!
   TBranch        *b_gLLP_visE20;   //!
   TBranch        *b_gLLP_visEz;   //!
   TBranch        *b_gLLP_visP;   //!
   TBranch        *b_gLLP_visPz;   //!
   TBranch        *b_gLLP_match_jet_pt;   //!
   TBranch        *b_gLLP_match_jet_index;   //!
   TBranch        *b_gLLP_match_jet_minDeltaR;   //!
   TBranch        *b_gLLP_decay_vertex_r;   //!
   TBranch        *b_gLLP_decay_vertex_x;   //!
   TBranch        *b_gLLP_decay_vertex_y;   //!
   TBranch        *b_gLLP_decay_vertex_z;   //!
   TBranch        *b_gLLP_deltaR;   //!
   TBranch        *b_gLLP_daughter_deltaR;   //!
   TBranch        *b_gLLP_daughter_pt;   //!
   TBranch        *b_gLLP_daughter_id;   //!
   TBranch        *b_gLLP_daughter_eta;   //!
   TBranch        *b_gLLP_daughter_phi;   //!
   TBranch        *b_gLLP_daughter_e;   //!
   TBranch        *b_gLLP_daughter_mass;   //!
   TBranch        *b_nGlobalMuons;   //!
   TBranch        *b_GlobalMuonPt;   //!
   TBranch        *b_GlobalMuonEta;   //!
   TBranch        *b_GlobalMuonPhi;   //!
   TBranch        *b_GlobalMuonLooseId;   //!
   TBranch        *b_nLeptons;   //!
   TBranch        *b_lepE;   //!
   TBranch        *b_lepPt;   //!
   TBranch        *b_lepEta;   //!
   TBranch        *b_lepPhi;   //!
   TBranch        *b_lepPdgId;   //!
   TBranch        *b_lepDZ;   //!
   TBranch        *b_lepPassId;   //!
   TBranch        *b_lepFromZ;   //!
   TBranch        *b_lepEff;   //!
   TBranch        *b_lepSF;   //!
   TBranch        *b_lepTriggerSF;   //!
   TBranch        *b_lepTightIdSF;   //!
   TBranch        *b_lepLooseIdSF;   //!
   TBranch        *b_lepTightIsoSF;   //!
   TBranch        *b_lepLooseIsoSF;   //!
   TBranch        *b_lepTag;   //!
   TBranch        *b_lepPassLooseIso;   //!
   TBranch        *b_lepPassTightIso;   //!
   TBranch        *b_lepPassVTightIso;   //!
   TBranch        *b_lepPassVVTightIso;   //!
   TBranch        *b_MT;   //!
   TBranch        *b_ZMass1;   //!
   TBranch        *b_ZMass;   //!
   TBranch        *b_ZPt;   //!
   TBranch        *b_ZEta;   //!
   TBranch        *b_ZPhi;   //!
   TBranch        *b_ZleptonIndex1;   //!
   TBranch        *b_ZleptonIndex2;   //!
   TBranch        *b_nJets;   //!
   TBranch        *b_jetE;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetTime;   //!
   TBranch        *b_jetPassId;   //!
   TBranch        *b_jetPtJESUp;   //!
   TBranch        *b_jetPtJESDown;   //!
   TBranch        *b_jetEJESUp;   //!
   TBranch        *b_jetEJESDown;   //!
   TBranch        *b_JecUnc;   //!
   TBranch        *b_jet_match_llp_pt;   //!
   TBranch        *b_jet_match_llp_index;   //!
   TBranch        *b_jet_match_llp_minDeltaR;   //!
   TBranch        *b_jet_match_genJet_pt;   //!
   TBranch        *b_jet_match_genJet_index;   //!
   TBranch        *b_jet_match_genJet_minDeltaR;   //!
   TBranch        *b_jetTightPassId;   //!
   TBranch        *b_HLTDecision;   //!
   TBranch        *b_METTrigger;   //!
   TBranch        *b_METNoMuTrigger;   //!
   TBranch        *b_jetChargedEMEnergyFraction;   //!
   TBranch        *b_jetChargedHadronEnergyFraction;   //!
   TBranch        *b_jetNeutralEMEnergyFraction;   //!
   TBranch        *b_jetNeutralHadronEnergyFraction;   //!

   LLPEvent(TTree *tree=0);
   virtual ~LLPEvent();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef LLPEvent_cxx
LLPEvent::LLPEvent(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ggH_HToSSTobbbb_MH-125_MS-15_ctau-1000_TuneCP5_13TeV-powheg-pythia8_59740pb_weighted.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("ggH_HToSSTobbbb_MH-125_MS-15_ctau-1000_TuneCP5_13TeV-powheg-pythia8_59740pb_weighted.root");
      }
      f->GetObject("MuonSystem",tree);

   }
   Init(tree);
}

LLPEvent::~LLPEvent()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t LLPEvent::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t LLPEvent::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void LLPEvent::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runNum", &runNum, &b_runNum);
   fChain->SetBranchAddress("MC_condition", &MC_condition, &b_MC_condition);
   fChain->SetBranchAddress("lumiSec", &lumiSec, &b_lumiSec);
   fChain->SetBranchAddress("evtNum", &evtNum, &b_evtNum);
   fChain->SetBranchAddress("mH", &mH, &b_mH);
   fChain->SetBranchAddress("mX", &mX, &b_mX);
   fChain->SetBranchAddress("ctau", &ctau, &b_ctau);
   fChain->SetBranchAddress("ZCategory", &ZCategory, &b_ZCategory);
   fChain->SetBranchAddress("category", &category, &b_category);
   fChain->SetBranchAddress("npv", &npv, &b_npv);
   fChain->SetBranchAddress("npu", &npu, &b_npu);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("higgsPtWeight", &higgsPtWeight, &b_higgsPtWeight);
   fChain->SetBranchAddress("higgsPtWeightSys", higgsPtWeightSys, &b_higgsPtWeightSys);
   fChain->SetBranchAddress("scaleWeights", scaleWeights, &b_scaleWeights);
   fChain->SetBranchAddress("lepOverallSF", &lepOverallSF, &b_lepOverallSF);
   fChain->SetBranchAddress("sf_facScaleUp", &sf_facScaleUp, &b_sf_facScaleUp);
   fChain->SetBranchAddress("sf_facScaleDown", &sf_facScaleDown, &b_sf_facScaleDown);
   fChain->SetBranchAddress("sf_renScaleUp", &sf_renScaleUp, &b_sf_renScaleUp);
   fChain->SetBranchAddress("sf_renScaleDown", &sf_renScaleDown, &b_sf_renScaleDown);
   fChain->SetBranchAddress("sf_facRenScaleUp", &sf_facRenScaleUp, &b_sf_facRenScaleUp);
   fChain->SetBranchAddress("sf_facRenScaleDown", &sf_facRenScaleDown, &b_sf_facRenScaleDown);
   fChain->SetBranchAddress("metSF", &metSF, &b_metSF);
   fChain->SetBranchAddress("pileupWeight", &pileupWeight, &b_pileupWeight);
   fChain->SetBranchAddress("pileupWeightUp", &pileupWeightUp, &b_pileupWeightUp);
   fChain->SetBranchAddress("pileupWeightDown", &pileupWeightDown, &b_pileupWeightDown);
   fChain->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, &b_Flag_HBHENoiseFilter);
   fChain->SetBranchAddress("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter, &b_Flag_BadPFMuonFilter);
   fChain->SetBranchAddress("Flag_HBHEIsoNoiseFilter", &Flag_HBHEIsoNoiseFilter, &b_Flag_HBHEIsoNoiseFilter);
   fChain->SetBranchAddress("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter, &b_Flag_CSCTightHaloFilter);
   fChain->SetBranchAddress("Flag_globalSuperTightHalo2016Filter", &Flag_globalSuperTightHalo2016Filter, &b_Flag_globalSuperTightHalo2016Filter);
   fChain->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices, &b_Flag_goodVertices);
   fChain->SetBranchAddress("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter, &b_Flag_ecalBadCalibFilter);
   fChain->SetBranchAddress("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter, &b_Flag_BadChargedCandidateFilter);
   fChain->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
   fChain->SetBranchAddress("Flag_all", &Flag_all, &b_Flag_all);
   fChain->SetBranchAddress("Flag2_HBHENoiseFilter", &Flag2_HBHENoiseFilter, &b_Flag2_HBHENoiseFilter);
   fChain->SetBranchAddress("Flag2_HBHEIsoNoiseFilter", &Flag2_HBHEIsoNoiseFilter, &b_Flag2_HBHEIsoNoiseFilter);
   fChain->SetBranchAddress("Flag2_BadPFMuonFilter", &Flag2_BadPFMuonFilter, &b_Flag2_BadPFMuonFilter);
   fChain->SetBranchAddress("Flag2_globalSuperTightHalo2016Filter", &Flag2_globalSuperTightHalo2016Filter, &b_Flag2_globalSuperTightHalo2016Filter);
   fChain->SetBranchAddress("Flag2_globalTightHalo2016Filter", &Flag2_globalTightHalo2016Filter, &b_Flag2_globalTightHalo2016Filter);
   fChain->SetBranchAddress("Flag2_BadChargedCandidateFilter", &Flag2_BadChargedCandidateFilter, &b_Flag2_BadChargedCandidateFilter);
   fChain->SetBranchAddress("Flag2_EcalDeadCellTriggerPrimitiveFilter", &Flag2_EcalDeadCellTriggerPrimitiveFilter, &b_Flag2_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("Flag2_ecalBadCalibFilter", &Flag2_ecalBadCalibFilter, &b_Flag2_ecalBadCalibFilter);
   fChain->SetBranchAddress("Flag2_eeBadScFilter", &Flag2_eeBadScFilter, &b_Flag2_eeBadScFilter);
   fChain->SetBranchAddress("Flag2_all", &Flag2_all, &b_Flag2_all);
   fChain->SetBranchAddress("EE_prefiring", &EE_prefiring, &b_EE_prefiring);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("met", &met, &b_met);
   fChain->SetBranchAddress("metNoMu", &metNoMu, &b_metNoMu);
   fChain->SetBranchAddress("metPhi", &metPhi, &b_metPhi);
   fChain->SetBranchAddress("metXYCorr", &metXYCorr, &b_metXYCorr);
   fChain->SetBranchAddress("metPhiXYCorr", &metPhiXYCorr, &b_metPhiXYCorr);
   fChain->SetBranchAddress("HT", &HT, &b_HT);
   fChain->SetBranchAddress("jetMet_dPhi", &jetMet_dPhi, &b_jetMet_dPhi);
   fChain->SetBranchAddress("jetMet_dPhiMin", &jetMet_dPhiMin, &b_jetMet_dPhiMin);
   fChain->SetBranchAddress("jetMet_dPhiMin4", &jetMet_dPhiMin4, &b_jetMet_dPhiMin4);
   fChain->SetBranchAddress("metJESUp", &metJESUp, &b_metJESUp);
   fChain->SetBranchAddress("metJESDown", &metJESDown, &b_metJESDown);
   fChain->SetBranchAddress("metPhiJESUp", &metPhiJESUp, &b_metPhiJESUp);
   fChain->SetBranchAddress("metPhiJESDown", &metPhiJESDown, &b_metPhiJESDown);
   fChain->SetBranchAddress("metJESUpSF", &metJESUpSF, &b_metJESUpSF);
   fChain->SetBranchAddress("metJESDownSF", &metJESDownSF, &b_metJESDownSF);
   fChain->SetBranchAddress("metEENoise", &metEENoise, &b_metEENoise);
   fChain->SetBranchAddress("metPhiEENoise", &metPhiEENoise, &b_metPhiEENoise);
   fChain->SetBranchAddress("metHEM", &metHEM, &b_metHEM);
   fChain->SetBranchAddress("metPhiHEM", &metPhiHEM, &b_metPhiHEM);
   fChain->SetBranchAddress("metEENoiseXYCorr", &metEENoiseXYCorr, &b_metEENoiseXYCorr);
   fChain->SetBranchAddress("metPhiEENoiseXYCorr", &metPhiEENoiseXYCorr, &b_metPhiEENoiseXYCorr);
   fChain->SetBranchAddress("metHEMXYCorr", &metHEMXYCorr, &b_metHEMXYCorr);
   fChain->SetBranchAddress("metPhiHEMXYCorr", &metPhiHEMXYCorr, &b_metPhiHEMXYCorr);
   fChain->SetBranchAddress("genMetPtTrue", &genMetPtTrue, &b_genMetPtTrue);
   fChain->SetBranchAddress("genMetPhiTrue", &genMetPhiTrue, &b_genMetPhiTrue);
   fChain->SetBranchAddress("genMetPtCalo", &genMetPtCalo, &b_genMetPtCalo);
   fChain->SetBranchAddress("genMetPhiCalo", &genMetPhiCalo, &b_genMetPhiCalo);
   fChain->SetBranchAddress("nGenParticle", &nGenParticle, &b_nGenParticle);
   fChain->SetBranchAddress("gParticleId", &gParticleId, &b_gParticleId);
   fChain->SetBranchAddress("gParticleStatus", &gParticleStatus, &b_gParticleStatus);
   fChain->SetBranchAddress("gParticleMotherId", &gParticleMotherId, &b_gParticleMotherId);
   fChain->SetBranchAddress("gParticleE", &gParticleE, &b_gParticleE);
   fChain->SetBranchAddress("gParticlePt", &gParticlePt, &b_gParticlePt);
   fChain->SetBranchAddress("gParticleEta", &gParticleEta, &b_gParticleEta);
   fChain->SetBranchAddress("gParticlePhi", &gParticlePhi, &b_gParticlePhi);
   fChain->SetBranchAddress("nGenJets", &nGenJets, &b_nGenJets);
   fChain->SetBranchAddress("genJetE", &genJetE, &b_genJetE);
   fChain->SetBranchAddress("genJetPt", &genJetPt, &b_genJetPt);
   fChain->SetBranchAddress("genJetEta", &genJetEta, &b_genJetEta);
   fChain->SetBranchAddress("genJetPhi", &genJetPhi, &b_genJetPhi);
   fChain->SetBranchAddress("genJetMET", &genJetMET, &b_genJetMET);
   fChain->SetBranchAddress("gWPt", &gWPt, &b_gWPt);
   fChain->SetBranchAddress("gLepId", &gLepId, &b_gLepId);
   fChain->SetBranchAddress("gLepPt", &gLepPt, &b_gLepPt);
   fChain->SetBranchAddress("gLepE", &gLepE, &b_gLepE);
   fChain->SetBranchAddress("gLepEta", &gLepEta, &b_gLepEta);
   fChain->SetBranchAddress("gLepPhi", &gLepPhi, &b_gLepPhi);
   fChain->SetBranchAddress("gHiggsPt", &gHiggsPt, &b_gHiggsPt);
   fChain->SetBranchAddress("gHiggsE", &gHiggsE, &b_gHiggsE);
   fChain->SetBranchAddress("gHiggsEta", &gHiggsEta, &b_gHiggsEta);
   fChain->SetBranchAddress("gHiggsPhi", &gHiggsPhi, &b_gHiggsPhi);
   fChain->SetBranchAddress("nCscRechits", &nCscRechits, &b_nCscRechits);
   fChain->SetBranchAddress("nEarlyCscRechits", &nEarlyCscRechits, &b_nEarlyCscRechits);
   fChain->SetBranchAddress("nLateCscRechits", &nLateCscRechits, &b_nLateCscRechits);
   fChain->SetBranchAddress("nEarly2CscRechits", &nEarly2CscRechits, &b_nEarly2CscRechits);
   fChain->SetBranchAddress("nLate2CscRechits", &nLate2CscRechits, &b_nLate2CscRechits);
   fChain->SetBranchAddress("nCscRings", &nCscRings, &b_nCscRings);
   fChain->SetBranchAddress("nCscPositiveYRechits", &nCscPositiveYRechits, &b_nCscPositiveYRechits);
   fChain->SetBranchAddress("nCscNegativeYRechits", &nCscNegativeYRechits, &b_nCscNegativeYRechits);
   fChain->SetBranchAddress("cscPosTpeak", &cscPosTpeak, &b_cscPosTpeak);
   fChain->SetBranchAddress("cscNegTpeak", &cscNegTpeak, &b_cscNegTpeak);
   fChain->SetBranchAddress("nCscRechitsChamberPlus11", &nCscRechitsChamberPlus11, &b_nCscRechitsChamberPlus11);
   fChain->SetBranchAddress("nCscRechitsChamberPlus12", &nCscRechitsChamberPlus12, &b_nCscRechitsChamberPlus12);
   fChain->SetBranchAddress("nCscRechitsChamberPlus13", &nCscRechitsChamberPlus13, &b_nCscRechitsChamberPlus13);
   fChain->SetBranchAddress("nCscRechitsChamberPlus21", &nCscRechitsChamberPlus21, &b_nCscRechitsChamberPlus21);
   fChain->SetBranchAddress("nCscRechitsChamberPlus22", &nCscRechitsChamberPlus22, &b_nCscRechitsChamberPlus22);
   fChain->SetBranchAddress("nCscRechitsChamberPlus31", &nCscRechitsChamberPlus31, &b_nCscRechitsChamberPlus31);
   fChain->SetBranchAddress("nCscRechitsChamberPlus32", &nCscRechitsChamberPlus32, &b_nCscRechitsChamberPlus32);
   fChain->SetBranchAddress("nCscRechitsChamberPlus41", &nCscRechitsChamberPlus41, &b_nCscRechitsChamberPlus41);
   fChain->SetBranchAddress("nCscRechitsChamberPlus42", &nCscRechitsChamberPlus42, &b_nCscRechitsChamberPlus42);
   fChain->SetBranchAddress("nCscRechitsChamberMinus11", &nCscRechitsChamberMinus11, &b_nCscRechitsChamberMinus11);
   fChain->SetBranchAddress("nCscRechitsChamberMinus12", &nCscRechitsChamberMinus12, &b_nCscRechitsChamberMinus12);
   fChain->SetBranchAddress("nCscRechitsChamberMinus13", &nCscRechitsChamberMinus13, &b_nCscRechitsChamberMinus13);
   fChain->SetBranchAddress("nCscRechitsChamberMinus21", &nCscRechitsChamberMinus21, &b_nCscRechitsChamberMinus21);
   fChain->SetBranchAddress("nCscRechitsChamberMinus22", &nCscRechitsChamberMinus22, &b_nCscRechitsChamberMinus22);
   fChain->SetBranchAddress("nCscRechitsChamberMinus31", &nCscRechitsChamberMinus31, &b_nCscRechitsChamberMinus31);
   fChain->SetBranchAddress("nCscRechitsChamberMinus32", &nCscRechitsChamberMinus32, &b_nCscRechitsChamberMinus32);
   fChain->SetBranchAddress("nCscRechitsChamberMinus41", &nCscRechitsChamberMinus41, &b_nCscRechitsChamberMinus41);
   fChain->SetBranchAddress("nCscRechitsChamberMinus42", &nCscRechitsChamberMinus42, &b_nCscRechitsChamberMinus42);
   fChain->SetBranchAddress("nRpc", &nRpc, &b_nRpc);
   fChain->SetBranchAddress("nDtSeg", &nDtSeg, &b_nDtSeg);
   fChain->SetBranchAddress("nDTRechits", &nDTRechits, &b_nDTRechits);
   fChain->SetBranchAddress("nDtRings", &nDtRings, &b_nDtRings);
   fChain->SetBranchAddress("nDtWheels25", &nDtWheels25, &b_nDtWheels25);
   fChain->SetBranchAddress("nDtStations25", &nDtStations25, &b_nDtStations25);
   fChain->SetBranchAddress("nDTPositiveYRechits", &nDTPositiveYRechits, &b_nDTPositiveYRechits);
   fChain->SetBranchAddress("nDTNegativeYRechits", &nDTNegativeYRechits, &b_nDTNegativeYRechits);
   fChain->SetBranchAddress("nDTRechitsWheelMinus2", &nDTRechitsWheelMinus2, &b_nDTRechitsWheelMinus2);
   fChain->SetBranchAddress("nDTRechitsWheelMinus1", &nDTRechitsWheelMinus1, &b_nDTRechitsWheelMinus1);
   fChain->SetBranchAddress("nDTRechitsWheel0", &nDTRechitsWheel0, &b_nDTRechitsWheel0);
   fChain->SetBranchAddress("nDTRechitsWheelPlus1", &nDTRechitsWheelPlus1, &b_nDTRechitsWheelPlus1);
   fChain->SetBranchAddress("nDTRechitsWheelPlus2", &nDTRechitsWheelPlus2, &b_nDTRechitsWheelPlus2);
   fChain->SetBranchAddress("nDTRechitsStation1", &nDTRechitsStation1, &b_nDTRechitsStation1);
   fChain->SetBranchAddress("nDTRechitsStation2", &nDTRechitsStation2, &b_nDTRechitsStation2);
   fChain->SetBranchAddress("nDTRechitsStation3", &nDTRechitsStation3, &b_nDTRechitsStation3);
   fChain->SetBranchAddress("nDTRechitsStation4", &nDTRechitsStation4, &b_nDTRechitsStation4);
   fChain->SetBranchAddress("nDTRechitsChamberMinus12", &nDTRechitsChamberMinus12, &b_nDTRechitsChamberMinus12);
   fChain->SetBranchAddress("nDTRechitsChamberMinus11", &nDTRechitsChamberMinus11, &b_nDTRechitsChamberMinus11);
   fChain->SetBranchAddress("nDTRechitsChamber10", &nDTRechitsChamber10, &b_nDTRechitsChamber10);
   fChain->SetBranchAddress("nDTRechitsChamberPlus11", &nDTRechitsChamberPlus11, &b_nDTRechitsChamberPlus11);
   fChain->SetBranchAddress("nDTRechitsChamberPlus12", &nDTRechitsChamberPlus12, &b_nDTRechitsChamberPlus12);
   fChain->SetBranchAddress("nDTRechitsChamberMinus22", &nDTRechitsChamberMinus22, &b_nDTRechitsChamberMinus22);
   fChain->SetBranchAddress("nDTRechitsChamberMinus21", &nDTRechitsChamberMinus21, &b_nDTRechitsChamberMinus21);
   fChain->SetBranchAddress("nDTRechitsChamber20", &nDTRechitsChamber20, &b_nDTRechitsChamber20);
   fChain->SetBranchAddress("nDTRechitsChamberPlus21", &nDTRechitsChamberPlus21, &b_nDTRechitsChamberPlus21);
   fChain->SetBranchAddress("nDTRechitsChamberPlus22", &nDTRechitsChamberPlus22, &b_nDTRechitsChamberPlus22);
   fChain->SetBranchAddress("nDTRechitsChamberMinus32", &nDTRechitsChamberMinus32, &b_nDTRechitsChamberMinus32);
   fChain->SetBranchAddress("nDTRechitsChamberMinus31", &nDTRechitsChamberMinus31, &b_nDTRechitsChamberMinus31);
   fChain->SetBranchAddress("nDTRechitsChamber30", &nDTRechitsChamber30, &b_nDTRechitsChamber30);
   fChain->SetBranchAddress("nDTRechitsChamberPlus31", &nDTRechitsChamberPlus31, &b_nDTRechitsChamberPlus31);
   fChain->SetBranchAddress("nDTRechitsChamberPlus32", &nDTRechitsChamberPlus32, &b_nDTRechitsChamberPlus32);
   fChain->SetBranchAddress("nDTRechitsChamberMinus42", &nDTRechitsChamberMinus42, &b_nDTRechitsChamberMinus42);
   fChain->SetBranchAddress("nDTRechitsChamberMinus41", &nDTRechitsChamberMinus41, &b_nDTRechitsChamberMinus41);
   fChain->SetBranchAddress("nDTRechitsChamber40", &nDTRechitsChamber40, &b_nDTRechitsChamber40);
   fChain->SetBranchAddress("nDTRechitsChamberPlus41", &nDTRechitsChamberPlus41, &b_nDTRechitsChamberPlus41);
   fChain->SetBranchAddress("nDTRechitsChamberPlus42", &nDTRechitsChamberPlus42, &b_nDTRechitsChamberPlus42);
   fChain->SetBranchAddress("nDTRechitsSector", nDTRechitsSector, &b_nDTRechitsSector);
   fChain->SetBranchAddress("nDTSegSector", nDTSegSector, &b_nDTSegSector);
   fChain->SetBranchAddress("cscRechitsStation", cscRechitsStation, &b_cscRechitsStation);
   fChain->SetBranchAddress("cscRechitsChamber", cscRechitsChamber, &b_cscRechitsChamber);
   fChain->SetBranchAddress("cscRechitsPhi", cscRechitsPhi, &b_cscRechitsPhi);
   fChain->SetBranchAddress("cscRechitsEta", cscRechitsEta, &b_cscRechitsEta);
   fChain->SetBranchAddress("dtRechitsX", &dtRechitsX, &b_dtRechitsX);
   fChain->SetBranchAddress("dtRechitsY", &dtRechitsY, &b_dtRechitsY);
   fChain->SetBranchAddress("dtRechitsZ", &dtRechitsZ, &b_dtRechitsZ);
   fChain->SetBranchAddress("dtRechitsEta", &dtRechitsEta, &b_dtRechitsEta);
   fChain->SetBranchAddress("dtRechitsPhi", &dtRechitsPhi, &b_dtRechitsPhi);
   fChain->SetBranchAddress("dtRechitsStation", &dtRechitsStation, &b_dtRechitsStation);
   fChain->SetBranchAddress("dtRechitsWheel", &dtRechitsWheel, &b_dtRechitsWheel);
   fChain->SetBranchAddress("dtRechitsClusterId", &dtRechitsClusterId, &b_dtRechitsClusterId);
   fChain->SetBranchAddress("rpcEta", &rpcEta, &b_rpcEta);
   fChain->SetBranchAddress("rpcPhi", &rpcPhi, &b_rpcPhi);
   fChain->SetBranchAddress("rpc_RE12", &rpc_RE12, &b_rpc_RE12);
   fChain->SetBranchAddress("rpc_RB1", &rpc_RB1, &b_rpc_RB1);
   fChain->SetBranchAddress("dtSegX", dtSegX, &b_dtSegX);
   fChain->SetBranchAddress("dtSegY", dtSegY, &b_dtSegY);
   fChain->SetBranchAddress("dtSegZ", dtSegZ, &b_dtSegZ);
   fChain->SetBranchAddress("dtSegEta", dtSegEta, &b_dtSegEta);
   fChain->SetBranchAddress("dtSegPhi", dtSegPhi, &b_dtSegPhi);
   fChain->SetBranchAddress("dtSegWheel", dtSegWheel, &b_dtSegWheel);
   fChain->SetBranchAddress("dtSegStation", dtSegStation, &b_dtSegStation);
   fChain->SetBranchAddress("nCscRechitClusters", &nCscRechitClusters, &b_nCscRechitClusters);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP", cscRechitCluster_match_gLLP, &b_cscRechitCluster_match_gLLP);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_minDeltaR", cscRechitCluster_match_gLLP_minDeltaR, &b_cscRechitCluster_match_gLLP_minDeltaR);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_index", cscRechitCluster_match_gLLP_index, &b_cscRechitCluster_match_gLLP_index);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_eta", cscRechitCluster_match_gLLP_eta, &b_cscRechitCluster_match_gLLP_eta);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_phi", cscRechitCluster_match_gLLP_phi, &b_cscRechitCluster_match_gLLP_phi);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_decay_r", cscRechitCluster_match_gLLP_decay_r, &b_cscRechitCluster_match_gLLP_decay_r);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_decay_x", cscRechitCluster_match_gLLP_decay_x, &b_cscRechitCluster_match_gLLP_decay_x);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_decay_y", cscRechitCluster_match_gLLP_decay_y, &b_cscRechitCluster_match_gLLP_decay_y);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_decay_z", cscRechitCluster_match_gLLP_decay_z, &b_cscRechitCluster_match_gLLP_decay_z);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_ctau", cscRechitCluster_match_gLLP_ctau, &b_cscRechitCluster_match_gLLP_ctau);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_beta", cscRechitCluster_match_gLLP_beta, &b_cscRechitCluster_match_gLLP_beta);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_csc", cscRechitCluster_match_gLLP_csc, &b_cscRechitCluster_match_gLLP_csc);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_dt", cscRechitCluster_match_gLLP_dt, &b_cscRechitCluster_match_gLLP_dt);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_multiplicity", cscRechitCluster_match_gLLP_multiplicity, &b_cscRechitCluster_match_gLLP_multiplicity);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_EM_multiplicity", cscRechitCluster_match_gLLP_EM_multiplicity, &b_cscRechitCluster_match_gLLP_EM_multiplicity);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_daughterKaon", cscRechitCluster_match_gLLP_daughterKaon, &b_cscRechitCluster_match_gLLP_daughterKaon);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_e", cscRechitCluster_match_gLLP_e, &b_cscRechitCluster_match_gLLP_e);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_pt", cscRechitCluster_match_gLLP_pt, &b_cscRechitCluster_match_gLLP_pt);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_EMFracE", cscRechitCluster_match_gLLP_EMFracE, &b_cscRechitCluster_match_gLLP_EMFracE);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_EMFracEz", cscRechitCluster_match_gLLP_EMFracEz, &b_cscRechitCluster_match_gLLP_EMFracEz);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_EMFracP", cscRechitCluster_match_gLLP_EMFracP, &b_cscRechitCluster_match_gLLP_EMFracP);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_EMFracPz", cscRechitCluster_match_gLLP_EMFracPz, &b_cscRechitCluster_match_gLLP_EMFracPz);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_visE", cscRechitCluster_match_gLLP_visE, &b_cscRechitCluster_match_gLLP_visE);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_visEz", cscRechitCluster_match_gLLP_visEz, &b_cscRechitCluster_match_gLLP_visEz);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_visP", cscRechitCluster_match_gLLP_visP, &b_cscRechitCluster_match_gLLP_visP);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_visPz", cscRechitCluster_match_gLLP_visPz, &b_cscRechitCluster_match_gLLP_visPz);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_lepdPhi", cscRechitCluster_match_gLLP_lepdPhi, &b_cscRechitCluster_match_gLLP_lepdPhi);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_daughter0_deltaR", cscRechitCluster_match_gLLP_daughter0_deltaR, &b_cscRechitCluster_match_gLLP_daughter0_deltaR);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_daughter1_deltaR", cscRechitCluster_match_gLLP_daughter1_deltaR, &b_cscRechitCluster_match_gLLP_daughter1_deltaR);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_daughter2_deltaR", cscRechitCluster_match_gLLP_daughter2_deltaR, &b_cscRechitCluster_match_gLLP_daughter2_deltaR);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_daughter3_deltaR", cscRechitCluster_match_gLLP_daughter3_deltaR, &b_cscRechitCluster_match_gLLP_daughter3_deltaR);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_other_daughter_deltaR", cscRechitCluster_match_gLLP_other_daughter_deltaR, &b_cscRechitCluster_match_gLLP_other_daughter_deltaR);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_other_daughter_index", cscRechitCluster_match_gLLP_other_daughter_index, &b_cscRechitCluster_match_gLLP_other_daughter_index);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_other_eta", cscRechitCluster_match_gLLP_other_eta, &b_cscRechitCluster_match_gLLP_other_eta);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_other_phi", cscRechitCluster_match_gLLP_other_phi, &b_cscRechitCluster_match_gLLP_other_phi);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_other_decay_r", cscRechitCluster_match_gLLP_other_decay_r, &b_cscRechitCluster_match_gLLP_other_decay_r);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_other_decay_x", cscRechitCluster_match_gLLP_other_decay_x, &b_cscRechitCluster_match_gLLP_other_decay_x);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_other_decay_y", cscRechitCluster_match_gLLP_other_decay_y, &b_cscRechitCluster_match_gLLP_other_decay_y);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_other_decay_z", cscRechitCluster_match_gLLP_other_decay_z, &b_cscRechitCluster_match_gLLP_other_decay_z);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_other_ctau", cscRechitCluster_match_gLLP_other_ctau, &b_cscRechitCluster_match_gLLP_other_ctau);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_other_beta", cscRechitCluster_match_gLLP_other_beta, &b_cscRechitCluster_match_gLLP_other_beta);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_other_csc", cscRechitCluster_match_gLLP_other_csc, &b_cscRechitCluster_match_gLLP_other_csc);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_other_e", cscRechitCluster_match_gLLP_other_e, &b_cscRechitCluster_match_gLLP_other_e);
   fChain->SetBranchAddress("cscRechitCluster_match_gLLP_other_pt", cscRechitCluster_match_gLLP_other_pt, &b_cscRechitCluster_match_gLLP_other_pt);
   fChain->SetBranchAddress("cscRechitClusterX", cscRechitClusterX, &b_cscRechitClusterX);
   fChain->SetBranchAddress("cscRechitClusterY", cscRechitClusterY, &b_cscRechitClusterY);
   fChain->SetBranchAddress("cscRechitClusterZ", cscRechitClusterZ, &b_cscRechitClusterZ);
   fChain->SetBranchAddress("cscRechitClusterTime", cscRechitClusterTime, &b_cscRechitClusterTime);
   fChain->SetBranchAddress("cscRechitClusterTimeWeighted", cscRechitClusterTimeWeighted, &b_cscRechitClusterTimeWeighted);
   fChain->SetBranchAddress("cscRechitClusterTimeTotal", cscRechitClusterTimeTotal, &b_cscRechitClusterTimeTotal);
   fChain->SetBranchAddress("cscRechitClusterTimeSpread", cscRechitClusterTimeSpread, &b_cscRechitClusterTimeSpread);
   fChain->SetBranchAddress("cscRechitClusterTimeSpreadWeighted", cscRechitClusterTimeSpreadWeighted, &b_cscRechitClusterTimeSpreadWeighted);
   fChain->SetBranchAddress("cscRechitClusterTimeSpreadWeightedAll", cscRechitClusterTimeSpreadWeightedAll, &b_cscRechitClusterTimeSpreadWeightedAll);
   fChain->SetBranchAddress("cscRechitClusterGenMuonDeltaR", cscRechitClusterGenMuonDeltaR, &b_cscRechitClusterGenMuonDeltaR);
   fChain->SetBranchAddress("cscRechitClusterXYSpread", cscRechitClusterXYSpread, &b_cscRechitClusterXYSpread);
   fChain->SetBranchAddress("cscRechitClusterMajorAxis", cscRechitClusterMajorAxis, &b_cscRechitClusterMajorAxis);
   fChain->SetBranchAddress("cscRechitClusterMinorAxis", cscRechitClusterMinorAxis, &b_cscRechitClusterMinorAxis);
   fChain->SetBranchAddress("cscRechitClusterEtaPhiSpread", cscRechitClusterEtaPhiSpread, &b_cscRechitClusterEtaPhiSpread);
   fChain->SetBranchAddress("cscRechitClusterPhiSpread", cscRechitClusterPhiSpread, &b_cscRechitClusterPhiSpread);
   fChain->SetBranchAddress("cscRechitClusterEtaSpread", cscRechitClusterEtaSpread, &b_cscRechitClusterEtaSpread);
   fChain->SetBranchAddress("cscRechitClusterDeltaRSpread", cscRechitClusterDeltaRSpread, &b_cscRechitClusterDeltaRSpread);
   fChain->SetBranchAddress("cscRechitClusterXSpread", cscRechitClusterXSpread, &b_cscRechitClusterXSpread);
   fChain->SetBranchAddress("cscRechitClusterRSpread", cscRechitClusterRSpread, &b_cscRechitClusterRSpread);
   fChain->SetBranchAddress("cscRechitClusterYSpread", cscRechitClusterYSpread, &b_cscRechitClusterYSpread);
   fChain->SetBranchAddress("cscRechitClusterZSpread", cscRechitClusterZSpread, &b_cscRechitClusterZSpread);
   fChain->SetBranchAddress("cscRechitClusterPhi", cscRechitClusterPhi, &b_cscRechitClusterPhi);
   fChain->SetBranchAddress("cscRechitClusterEta", cscRechitClusterEta, &b_cscRechitClusterEta);
   fChain->SetBranchAddress("cscRechitClusterJetVetoPt", cscRechitClusterJetVetoPt, &b_cscRechitClusterJetVetoPt);
   fChain->SetBranchAddress("cscRechitClusterJetVetoPtJESDown", cscRechitClusterJetVetoPtJESDown, &b_cscRechitClusterJetVetoPtJESDown);
   fChain->SetBranchAddress("cscRechitClusterJetVetoPtJESUp", cscRechitClusterJetVetoPtJESUp, &b_cscRechitClusterJetVetoPtJESUp);
   fChain->SetBranchAddress("cscRechitClusterJetVetoEta", cscRechitClusterJetVetoEta, &b_cscRechitClusterJetVetoEta);
   fChain->SetBranchAddress("cscRechitClusterJetVetoPhi", cscRechitClusterJetVetoPhi, &b_cscRechitClusterJetVetoPhi);
   fChain->SetBranchAddress("cscRechitClusterTightJetVetoPt", cscRechitClusterTightJetVetoPt, &b_cscRechitClusterTightJetVetoPt);
   fChain->SetBranchAddress("cscRechitClusterTightJetVetoPtJESDown", cscRechitClusterTightJetVetoPtJESDown, &b_cscRechitClusterTightJetVetoPtJESDown);
   fChain->SetBranchAddress("cscRechitClusterTightJetVetoPtJESUp", cscRechitClusterTightJetVetoPtJESUp, &b_cscRechitClusterTightJetVetoPtJESUp);
   fChain->SetBranchAddress("cscRechitClusterTightJetVetoEta", cscRechitClusterTightJetVetoEta, &b_cscRechitClusterTightJetVetoEta);
   fChain->SetBranchAddress("cscRechitClusterTightJetVetoPhi", cscRechitClusterTightJetVetoPhi, &b_cscRechitClusterTightJetVetoPhi);
   fChain->SetBranchAddress("cscRechitClusterJetVetoE", cscRechitClusterJetVetoE, &b_cscRechitClusterJetVetoE);
   fChain->SetBranchAddress("cscRechitClusterGenJetVetoPt", cscRechitClusterGenJetVetoPt, &b_cscRechitClusterGenJetVetoPt);
   fChain->SetBranchAddress("cscRechitClusterGenJetVetoE", cscRechitClusterGenJetVetoE, &b_cscRechitClusterGenJetVetoE);
   fChain->SetBranchAddress("cscRechitClusterMuonVetoPt", cscRechitClusterMuonVetoPt, &b_cscRechitClusterMuonVetoPt);
   fChain->SetBranchAddress("cscRechitClusterMuonVetoE", cscRechitClusterMuonVetoE, &b_cscRechitClusterMuonVetoE);
   fChain->SetBranchAddress("cscRechitClusterMuonVetoPhi", cscRechitClusterMuonVetoPhi, &b_cscRechitClusterMuonVetoPhi);
   fChain->SetBranchAddress("cscRechitClusterMuonVetoEta", cscRechitClusterMuonVetoEta, &b_cscRechitClusterMuonVetoEta);
   fChain->SetBranchAddress("cscRechitClusterJetVetoElectronEnergyFraction", cscRechitClusterJetVetoElectronEnergyFraction, &b_cscRechitClusterJetVetoElectronEnergyFraction);
   fChain->SetBranchAddress("cscRechitClusterJetVetoPhotonEnergyFraction", cscRechitClusterJetVetoPhotonEnergyFraction, &b_cscRechitClusterJetVetoPhotonEnergyFraction);
   fChain->SetBranchAddress("cscRechitClusterJetVetoChargedHadronEnergyFraction", cscRechitClusterJetVetoChargedHadronEnergyFraction, &b_cscRechitClusterJetVetoChargedHadronEnergyFraction);
   fChain->SetBranchAddress("cscRechitClusterJetVetoNeutralHadronEnergyFraction", cscRechitClusterJetVetoNeutralHadronEnergyFraction, &b_cscRechitClusterJetVetoNeutralHadronEnergyFraction);
   fChain->SetBranchAddress("cscRechitClusterJetVetoMuonEnergyFraction", cscRechitClusterJetVetoMuonEnergyFraction, &b_cscRechitClusterJetVetoMuonEnergyFraction);
   fChain->SetBranchAddress("cscRechitClusterJetVetoPt_0p6", cscRechitClusterJetVetoPt_0p6, &b_cscRechitClusterJetVetoPt_0p6);
   fChain->SetBranchAddress("cscRechitClusterJetVetoPt_0p8", cscRechitClusterJetVetoPt_0p8, &b_cscRechitClusterJetVetoPt_0p8);
   fChain->SetBranchAddress("cscRechitClusterJetVetoE_0p6", cscRechitClusterJetVetoE_0p6, &b_cscRechitClusterJetVetoE_0p6);
   fChain->SetBranchAddress("cscRechitClusterJetVetoE_0p8", cscRechitClusterJetVetoE_0p8, &b_cscRechitClusterJetVetoE_0p8);
   fChain->SetBranchAddress("cscRechitClusterMuonVetoPt_0p6", cscRechitClusterMuonVetoPt_0p6, &b_cscRechitClusterMuonVetoPt_0p6);
   fChain->SetBranchAddress("cscRechitClusterMuonVetoPt_0p8", cscRechitClusterMuonVetoPt_0p8, &b_cscRechitClusterMuonVetoPt_0p8);
   fChain->SetBranchAddress("cscRechitClusterMuonVetoE_0p6", cscRechitClusterMuonVetoE_0p6, &b_cscRechitClusterMuonVetoE_0p6);
   fChain->SetBranchAddress("cscRechitClusterMuonVetoE_0p8", cscRechitClusterMuonVetoE_0p8, &b_cscRechitClusterMuonVetoE_0p8);
   fChain->SetBranchAddress("cscRechitClusterMuonVetoLooseIso", cscRechitClusterMuonVetoLooseIso, &b_cscRechitClusterMuonVetoLooseIso);
   fChain->SetBranchAddress("cscRechitClusterMuonVetoTightIso", cscRechitClusterMuonVetoTightIso, &b_cscRechitClusterMuonVetoTightIso);
   fChain->SetBranchAddress("cscRechitClusterMuonVetoVTightIso", cscRechitClusterMuonVetoVTightIso, &b_cscRechitClusterMuonVetoVTightIso);
   fChain->SetBranchAddress("cscRechitClusterMuonVetoVVTightIso", cscRechitClusterMuonVetoVVTightIso, &b_cscRechitClusterMuonVetoVVTightIso);
   fChain->SetBranchAddress("cscRechitClusterMuonVetoTightId", cscRechitClusterMuonVetoTightId, &b_cscRechitClusterMuonVetoTightId);
   fChain->SetBranchAddress("cscRechitClusterMuonVetoLooseId", cscRechitClusterMuonVetoLooseId, &b_cscRechitClusterMuonVetoLooseId);
   fChain->SetBranchAddress("cscRechitClusterMuonVetoGlobal", cscRechitClusterMuonVetoGlobal, &b_cscRechitClusterMuonVetoGlobal);
   fChain->SetBranchAddress("cscRechitClusterMuonVetoIso", cscRechitClusterMuonVetoIso, &b_cscRechitClusterMuonVetoIso);
   fChain->SetBranchAddress("cscRechitClusterIsoMuonVetoPt", cscRechitClusterIsoMuonVetoPt, &b_cscRechitClusterIsoMuonVetoPt);
   fChain->SetBranchAddress("cscRechitClusterIsoMuonVetoE", cscRechitClusterIsoMuonVetoE, &b_cscRechitClusterIsoMuonVetoE);
   fChain->SetBranchAddress("cscRechitClusterIsoMuonVetoPhi", cscRechitClusterIsoMuonVetoPhi, &b_cscRechitClusterIsoMuonVetoPhi);
   fChain->SetBranchAddress("cscRechitClusterIsoMuonVetoEta", cscRechitClusterIsoMuonVetoEta, &b_cscRechitClusterIsoMuonVetoEta);
   fChain->SetBranchAddress("cscRechitClusterGenMuonVetoPt", cscRechitClusterGenMuonVetoPt, &b_cscRechitClusterGenMuonVetoPt);
   fChain->SetBranchAddress("cscRechitClusterMuonVetoType", cscRechitClusterMuonVetoType, &b_cscRechitClusterMuonVetoType);
   fChain->SetBranchAddress("cscRechitClusterZLep1", cscRechitClusterZLep1, &b_cscRechitClusterZLep1);
   fChain->SetBranchAddress("cscRechitClusterZLep2", cscRechitClusterZLep2, &b_cscRechitClusterZLep2);
   fChain->SetBranchAddress("cscRechitClusterZLep1Tag", cscRechitClusterZLep1Tag, &b_cscRechitClusterZLep1Tag);
   fChain->SetBranchAddress("cscRechitClusterZLep2Tag", cscRechitClusterZLep2Tag, &b_cscRechitClusterZLep2Tag);
   fChain->SetBranchAddress("cscRechitClusterZLep1Id", cscRechitClusterZLep1Id, &b_cscRechitClusterZLep1Id);
   fChain->SetBranchAddress("cscRechitClusterZLep2Id", cscRechitClusterZLep2Id, &b_cscRechitClusterZLep2Id);
   fChain->SetBranchAddress("cscRechitCluster2ZLep1LooseIso", cscRechitCluster2ZLep1LooseIso, &b_cscRechitCluster2ZLep1LooseIso);
   fChain->SetBranchAddress("cscRechitClusterZLep1TightIso", cscRechitClusterZLep1TightIso, &b_cscRechitClusterZLep1TightIso);
   fChain->SetBranchAddress("cscRechitClusterZLep1VTightIso", cscRechitClusterZLep1VTightIso, &b_cscRechitClusterZLep1VTightIso);
   fChain->SetBranchAddress("cscRechitClusterZLep1VVTightIso", cscRechitClusterZLep1VVTightIso, &b_cscRechitClusterZLep1VVTightIso);
   fChain->SetBranchAddress("cscRechitClusterZLep1TightId", cscRechitClusterZLep1TightId, &b_cscRechitClusterZLep1TightId);
   fChain->SetBranchAddress("cscRechitClusterZLep2LooseIso", cscRechitClusterZLep2LooseIso, &b_cscRechitClusterZLep2LooseIso);
   fChain->SetBranchAddress("cscRechitClusterZLep2TightIso", cscRechitClusterZLep2TightIso, &b_cscRechitClusterZLep2TightIso);
   fChain->SetBranchAddress("cscRechitClusterZLep2VTightIso", cscRechitClusterZLep2VTightIso, &b_cscRechitClusterZLep2VTightIso);
   fChain->SetBranchAddress("cscRechitClusterZLep2VVTightIso", cscRechitClusterZLep2VVTightIso, &b_cscRechitClusterZLep2VVTightIso);
   fChain->SetBranchAddress("cscRechitClusterZLep2TightId", cscRechitClusterZLep2TightId, &b_cscRechitClusterZLep2TightId);
   fChain->SetBranchAddress("cscRechitClusterSize", cscRechitClusterSize, &b_cscRechitClusterSize);
   fChain->SetBranchAddress("cscRechitClusterMe11Ratio", cscRechitClusterMe11Ratio, &b_cscRechitClusterMe11Ratio);
   fChain->SetBranchAddress("cscRechitClusterMe12Ratio", cscRechitClusterMe12Ratio, &b_cscRechitClusterMe12Ratio);
   fChain->SetBranchAddress("cscRechitClusterNStation", cscRechitClusterNStation, &b_cscRechitClusterNStation);
   fChain->SetBranchAddress("cscRechitClusterNStation5", cscRechitClusterNStation5, &b_cscRechitClusterNStation5);
   fChain->SetBranchAddress("cscRechitClusterNStation10", cscRechitClusterNStation10, &b_cscRechitClusterNStation10);
   fChain->SetBranchAddress("cscRechitClusterNStation10perc", cscRechitClusterNStation10perc, &b_cscRechitClusterNStation10perc);
   fChain->SetBranchAddress("cscRechitClusterAvgStation", cscRechitClusterAvgStation, &b_cscRechitClusterAvgStation);
   fChain->SetBranchAddress("cscRechitClusterAvgStation5", cscRechitClusterAvgStation5, &b_cscRechitClusterAvgStation5);
   fChain->SetBranchAddress("cscRechitClusterAvgStation10", cscRechitClusterAvgStation10, &b_cscRechitClusterAvgStation10);
   fChain->SetBranchAddress("cscRechitClusterAvgStation10perc", cscRechitClusterAvgStation10perc, &b_cscRechitClusterAvgStation10perc);
   fChain->SetBranchAddress("cscRechitClusterMaxStation", cscRechitClusterMaxStation, &b_cscRechitClusterMaxStation);
   fChain->SetBranchAddress("cscRechitClusterMaxStationRatio", cscRechitClusterMaxStationRatio, &b_cscRechitClusterMaxStationRatio);
   fChain->SetBranchAddress("cscRechitClusterNChamber", cscRechitClusterNChamber, &b_cscRechitClusterNChamber);
   fChain->SetBranchAddress("cscRechitClusterMaxChamber", cscRechitClusterMaxChamber, &b_cscRechitClusterMaxChamber);
   fChain->SetBranchAddress("cscRechitClusterMaxChamberRatio", cscRechitClusterMaxChamberRatio, &b_cscRechitClusterMaxChamberRatio);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberPlus11", cscRechitClusterNRechitChamberPlus11, &b_cscRechitClusterNRechitChamberPlus11);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberPlus12", cscRechitClusterNRechitChamberPlus12, &b_cscRechitClusterNRechitChamberPlus12);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberPlus13", cscRechitClusterNRechitChamberPlus13, &b_cscRechitClusterNRechitChamberPlus13);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberPlus21", cscRechitClusterNRechitChamberPlus21, &b_cscRechitClusterNRechitChamberPlus21);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberPlus22", cscRechitClusterNRechitChamberPlus22, &b_cscRechitClusterNRechitChamberPlus22);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberPlus31", cscRechitClusterNRechitChamberPlus31, &b_cscRechitClusterNRechitChamberPlus31);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberPlus32", cscRechitClusterNRechitChamberPlus32, &b_cscRechitClusterNRechitChamberPlus32);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberPlus41", cscRechitClusterNRechitChamberPlus41, &b_cscRechitClusterNRechitChamberPlus41);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberPlus42", cscRechitClusterNRechitChamberPlus42, &b_cscRechitClusterNRechitChamberPlus42);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberMinus11", cscRechitClusterNRechitChamberMinus11, &b_cscRechitClusterNRechitChamberMinus11);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberMinus12", cscRechitClusterNRechitChamberMinus12, &b_cscRechitClusterNRechitChamberMinus12);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberMinus13", cscRechitClusterNRechitChamberMinus13, &b_cscRechitClusterNRechitChamberMinus13);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberMinus21", cscRechitClusterNRechitChamberMinus21, &b_cscRechitClusterNRechitChamberMinus21);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberMinus22", cscRechitClusterNRechitChamberMinus22, &b_cscRechitClusterNRechitChamberMinus22);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberMinus31", cscRechitClusterNRechitChamberMinus31, &b_cscRechitClusterNRechitChamberMinus31);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberMinus32", cscRechitClusterNRechitChamberMinus32, &b_cscRechitClusterNRechitChamberMinus32);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberMinus41", cscRechitClusterNRechitChamberMinus41, &b_cscRechitClusterNRechitChamberMinus41);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberMinus42", cscRechitClusterNRechitChamberMinus42, &b_cscRechitClusterNRechitChamberMinus42);
   fChain->SetBranchAddress("cscRechitClusterMet_dPhi", cscRechitClusterMet_dPhi, &b_cscRechitClusterMet_dPhi);
   fChain->SetBranchAddress("cscRechitClusterMetJESUp_dPhi", cscRechitClusterMetJESUp_dPhi, &b_cscRechitClusterMetJESUp_dPhi);
   fChain->SetBranchAddress("cscRechitClusterMetJESDown_dPhi", cscRechitClusterMetJESDown_dPhi, &b_cscRechitClusterMetJESDown_dPhi);
   fChain->SetBranchAddress("cscRechitClusterMetXYCorr_dPhi", cscRechitClusterMetXYCorr_dPhi, &b_cscRechitClusterMetXYCorr_dPhi);
   fChain->SetBranchAddress("cscRechitCluster_match_cscRechits_0p4", cscRechitCluster_match_cscRechits_0p4, &b_cscRechitCluster_match_cscRechits_0p4);
   fChain->SetBranchAddress("cscRechitClusterNLayersChamberPlus11", cscRechitClusterNLayersChamberPlus11, &b_cscRechitClusterNLayersChamberPlus11);
   fChain->SetBranchAddress("cscRechitClusterNLayersChamberPlus12", cscRechitClusterNLayersChamberPlus12, &b_cscRechitClusterNLayersChamberPlus12);
   fChain->SetBranchAddress("cscRechitClusterNLayersChamberPlus13", cscRechitClusterNLayersChamberPlus13, &b_cscRechitClusterNLayersChamberPlus13);
   fChain->SetBranchAddress("cscRechitClusterNLayersChamberPlus21", cscRechitClusterNLayersChamberPlus21, &b_cscRechitClusterNLayersChamberPlus21);
   fChain->SetBranchAddress("cscRechitClusterNLayersChamberPlus22", cscRechitClusterNLayersChamberPlus22, &b_cscRechitClusterNLayersChamberPlus22);
   fChain->SetBranchAddress("cscRechitClusterNLayersChamberPlus31", cscRechitClusterNLayersChamberPlus31, &b_cscRechitClusterNLayersChamberPlus31);
   fChain->SetBranchAddress("cscRechitClusterNLayersChamberPlus32", cscRechitClusterNLayersChamberPlus32, &b_cscRechitClusterNLayersChamberPlus32);
   fChain->SetBranchAddress("cscRechitClusterNLayersChamberPlus41", cscRechitClusterNLayersChamberPlus41, &b_cscRechitClusterNLayersChamberPlus41);
   fChain->SetBranchAddress("cscRechitClusterNLayersChamberPlus42", cscRechitClusterNLayersChamberPlus42, &b_cscRechitClusterNLayersChamberPlus42);
   fChain->SetBranchAddress("cscRechitClusterNLayersChamberMinus11", cscRechitClusterNLayersChamberMinus11, &b_cscRechitClusterNLayersChamberMinus11);
   fChain->SetBranchAddress("cscRechitClusterNLayersChamberMinus12", cscRechitClusterNLayersChamberMinus12, &b_cscRechitClusterNLayersChamberMinus12);
   fChain->SetBranchAddress("cscRechitClusterNLayersChamberMinus13", cscRechitClusterNLayersChamberMinus13, &b_cscRechitClusterNLayersChamberMinus13);
   fChain->SetBranchAddress("cscRechitClusterNLayersChamberMinus21", cscRechitClusterNLayersChamberMinus21, &b_cscRechitClusterNLayersChamberMinus21);
   fChain->SetBranchAddress("cscRechitClusterNLayersChamberMinus22", cscRechitClusterNLayersChamberMinus22, &b_cscRechitClusterNLayersChamberMinus22);
   fChain->SetBranchAddress("cscRechitClusterNLayersChamberMinus31", cscRechitClusterNLayersChamberMinus31, &b_cscRechitClusterNLayersChamberMinus31);
   fChain->SetBranchAddress("cscRechitClusterNLayersChamberMinus32", cscRechitClusterNLayersChamberMinus32, &b_cscRechitClusterNLayersChamberMinus32);
   fChain->SetBranchAddress("cscRechitClusterNLayersChamberMinus41", cscRechitClusterNLayersChamberMinus41, &b_cscRechitClusterNLayersChamberMinus41);
   fChain->SetBranchAddress("cscRechitClusterNLayersChamberMinus42", cscRechitClusterNLayersChamberMinus42, &b_cscRechitClusterNLayersChamberMinus42);
   fChain->SetBranchAddress("cscRechitClusterMetHEM_dPhi", cscRechitClusterMetHEM_dPhi, &b_cscRechitClusterMetHEM_dPhi);
   fChain->SetBranchAddress("cscRechitClusterMetHEMXYCorr_dPhi", cscRechitClusterMetHEMXYCorr_dPhi, &b_cscRechitClusterMetHEMXYCorr_dPhi);
   fChain->SetBranchAddress("cscRechitClusterMetEENoise_dPhi", cscRechitClusterMetEENoise_dPhi, &b_cscRechitClusterMetEENoise_dPhi);
   fChain->SetBranchAddress("cscRechitClusterMetEENoiseXYCorr_dPhi", cscRechitClusterMetEENoiseXYCorr_dPhi, &b_cscRechitClusterMetEENoiseXYCorr_dPhi);
   fChain->SetBranchAddress("cscRechitClusterMetJesUp_dPhi", cscRechitClusterMetJesUp_dPhi, &b_cscRechitClusterMetJesUp_dPhi);
   fChain->SetBranchAddress("cscRechitClusterMetJesDown_dPhi", cscRechitClusterMetJesDown_dPhi, &b_cscRechitClusterMetJesDown_dPhi);
   fChain->SetBranchAddress("cscRechitCluster_match_dtRechits_phi0p2", cscRechitCluster_match_dtRechits_phi0p2, &b_cscRechitCluster_match_dtRechits_phi0p2);
   fChain->SetBranchAddress("cscRechitCluster_match_dtRechits_0p4", cscRechitCluster_match_dtRechits_0p4, &b_cscRechitCluster_match_dtRechits_0p4);
   fChain->SetBranchAddress("cscRechitCluster_match_MB1_0p4", cscRechitCluster_match_MB1_0p4, &b_cscRechitCluster_match_MB1_0p4);
   fChain->SetBranchAddress("cscRechitCluster_match_dtSeg_0p4", cscRechitCluster_match_dtSeg_0p4, &b_cscRechitCluster_match_dtSeg_0p4);
   fChain->SetBranchAddress("cscRechitCluster_match_MB1Seg_0p4", cscRechitCluster_match_MB1Seg_0p4, &b_cscRechitCluster_match_MB1Seg_0p4);
   fChain->SetBranchAddress("cscRechitCluster_match_RB1_0p4", cscRechitCluster_match_RB1_0p4, &b_cscRechitCluster_match_RB1_0p4);
   fChain->SetBranchAddress("cscRechitCluster_match_RE12_0p4", cscRechitCluster_match_RE12_0p4, &b_cscRechitCluster_match_RE12_0p4);
   fChain->SetBranchAddress("cscRechitCluster_match_cluster_dR", cscRechitCluster_match_cluster_dR, &b_cscRechitCluster_match_cluster_dR);
   fChain->SetBranchAddress("cscRechitCluster_match_cluster_index", cscRechitCluster_match_cluster_index, &b_cscRechitCluster_match_cluster_index);
   fChain->SetBranchAddress("nDtRechitClusters", &nDtRechitClusters, &b_nDtRechitClusters);
   fChain->SetBranchAddress("nDtRechitClusters2", &nDtRechitClusters2, &b_nDtRechitClusters2);
   fChain->SetBranchAddress("dtRechitClusterMaxDPhi", dtRechitClusterMaxDPhi, &b_dtRechitClusterMaxDPhi);
   fChain->SetBranchAddress("dtRechitClusterMaxDPhi_index", dtRechitClusterMaxDPhi_index, &b_dtRechitClusterMaxDPhi_index);
   fChain->SetBranchAddress("dtRechitClusterNSegStation1", dtRechitClusterNSegStation1, &b_dtRechitClusterNSegStation1);
   fChain->SetBranchAddress("dtRechitClusterNSegStation2", dtRechitClusterNSegStation2, &b_dtRechitClusterNSegStation2);
   fChain->SetBranchAddress("dtRechitClusterNSegStation3", dtRechitClusterNSegStation3, &b_dtRechitClusterNSegStation3);
   fChain->SetBranchAddress("dtRechitClusterNSegStation4", dtRechitClusterNSegStation4, &b_dtRechitClusterNSegStation4);
   fChain->SetBranchAddress("dtRechitClusterNOppositeSegStation1", dtRechitClusterNOppositeSegStation1, &b_dtRechitClusterNOppositeSegStation1);
   fChain->SetBranchAddress("dtRechitClusterNOppositeSegStation2", dtRechitClusterNOppositeSegStation2, &b_dtRechitClusterNOppositeSegStation2);
   fChain->SetBranchAddress("dtRechitClusterNOppositeSegStation3", dtRechitClusterNOppositeSegStation3, &b_dtRechitClusterNOppositeSegStation3);
   fChain->SetBranchAddress("dtRechitClusterNOppositeSegStation4", dtRechitClusterNOppositeSegStation4, &b_dtRechitClusterNOppositeSegStation4);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP", dtRechitCluster_match_gLLP, &b_dtRechitCluster_match_gLLP);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_minDeltaR", dtRechitCluster_match_gLLP_minDeltaR, &b_dtRechitCluster_match_gLLP_minDeltaR);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_index", dtRechitCluster_match_gLLP_index, &b_dtRechitCluster_match_gLLP_index);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_eta", dtRechitCluster_match_gLLP_eta, &b_dtRechitCluster_match_gLLP_eta);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_phi", dtRechitCluster_match_gLLP_phi, &b_dtRechitCluster_match_gLLP_phi);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_decay_r", dtRechitCluster_match_gLLP_decay_r, &b_dtRechitCluster_match_gLLP_decay_r);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_decay_x", dtRechitCluster_match_gLLP_decay_x, &b_dtRechitCluster_match_gLLP_decay_x);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_decay_y", dtRechitCluster_match_gLLP_decay_y, &b_dtRechitCluster_match_gLLP_decay_y);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_decay_z", dtRechitCluster_match_gLLP_decay_z, &b_dtRechitCluster_match_gLLP_decay_z);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_ctau", dtRechitCluster_match_gLLP_ctau, &b_dtRechitCluster_match_gLLP_ctau);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_beta", dtRechitCluster_match_gLLP_beta, &b_dtRechitCluster_match_gLLP_beta);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_csc", dtRechitCluster_match_gLLP_csc, &b_dtRechitCluster_match_gLLP_csc);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_dt", dtRechitCluster_match_gLLP_dt, &b_dtRechitCluster_match_gLLP_dt);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_multiplicity", dtRechitCluster_match_gLLP_multiplicity, &b_dtRechitCluster_match_gLLP_multiplicity);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_EM_multiplicity", dtRechitCluster_match_gLLP_EM_multiplicity, &b_dtRechitCluster_match_gLLP_EM_multiplicity);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_daughterKaon", dtRechitCluster_match_gLLP_daughterKaon, &b_dtRechitCluster_match_gLLP_daughterKaon);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_e", dtRechitCluster_match_gLLP_e, &b_dtRechitCluster_match_gLLP_e);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_pt", dtRechitCluster_match_gLLP_pt, &b_dtRechitCluster_match_gLLP_pt);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_EMFracE", dtRechitCluster_match_gLLP_EMFracE, &b_dtRechitCluster_match_gLLP_EMFracE);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_EMFracEz", dtRechitCluster_match_gLLP_EMFracEz, &b_dtRechitCluster_match_gLLP_EMFracEz);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_EMFracP", dtRechitCluster_match_gLLP_EMFracP, &b_dtRechitCluster_match_gLLP_EMFracP);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_EMFracPz", dtRechitCluster_match_gLLP_EMFracPz, &b_dtRechitCluster_match_gLLP_EMFracPz);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_visE", dtRechitCluster_match_gLLP_visE, &b_dtRechitCluster_match_gLLP_visE);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_visEz", dtRechitCluster_match_gLLP_visEz, &b_dtRechitCluster_match_gLLP_visEz);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_visP", dtRechitCluster_match_gLLP_visP, &b_dtRechitCluster_match_gLLP_visP);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_visPz", dtRechitCluster_match_gLLP_visPz, &b_dtRechitCluster_match_gLLP_visPz);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_lepdPhi", dtRechitCluster_match_gLLP_lepdPhi, &b_dtRechitCluster_match_gLLP_lepdPhi);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_daughter0_deltaR", dtRechitCluster_match_gLLP_daughter0_deltaR, &b_dtRechitCluster_match_gLLP_daughter0_deltaR);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_daughter1_deltaR", dtRechitCluster_match_gLLP_daughter1_deltaR, &b_dtRechitCluster_match_gLLP_daughter1_deltaR);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_daughter2_deltaR", dtRechitCluster_match_gLLP_daughter2_deltaR, &b_dtRechitCluster_match_gLLP_daughter2_deltaR);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_daughter3_deltaR", dtRechitCluster_match_gLLP_daughter3_deltaR, &b_dtRechitCluster_match_gLLP_daughter3_deltaR);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_other_daughter_deltaR", dtRechitCluster_match_gLLP_other_daughter_deltaR, &b_dtRechitCluster_match_gLLP_other_daughter_deltaR);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_other_daughter_index", dtRechitCluster_match_gLLP_other_daughter_index, &b_dtRechitCluster_match_gLLP_other_daughter_index);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_other_eta", dtRechitCluster_match_gLLP_other_eta, &b_dtRechitCluster_match_gLLP_other_eta);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_other_phi", dtRechitCluster_match_gLLP_other_phi, &b_dtRechitCluster_match_gLLP_other_phi);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_other_decay_r", dtRechitCluster_match_gLLP_other_decay_r, &b_dtRechitCluster_match_gLLP_other_decay_r);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_other_decay_x", dtRechitCluster_match_gLLP_other_decay_x, &b_dtRechitCluster_match_gLLP_other_decay_x);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_other_decay_y", dtRechitCluster_match_gLLP_other_decay_y, &b_dtRechitCluster_match_gLLP_other_decay_y);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_other_decay_z", dtRechitCluster_match_gLLP_other_decay_z, &b_dtRechitCluster_match_gLLP_other_decay_z);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_other_ctau", dtRechitCluster_match_gLLP_other_ctau, &b_dtRechitCluster_match_gLLP_other_ctau);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_other_beta", dtRechitCluster_match_gLLP_other_beta, &b_dtRechitCluster_match_gLLP_other_beta);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_other_csc", dtRechitCluster_match_gLLP_other_csc, &b_dtRechitCluster_match_gLLP_other_csc);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_other_e", dtRechitCluster_match_gLLP_other_e, &b_dtRechitCluster_match_gLLP_other_e);
   fChain->SetBranchAddress("dtRechitCluster_match_gLLP_other_pt", dtRechitCluster_match_gLLP_other_pt, &b_dtRechitCluster_match_gLLP_other_pt);
   fChain->SetBranchAddress("dtRechitClusterX", dtRechitClusterX, &b_dtRechitClusterX);
   fChain->SetBranchAddress("dtRechitClusterY", dtRechitClusterY, &b_dtRechitClusterY);
   fChain->SetBranchAddress("dtRechitClusterZ", dtRechitClusterZ, &b_dtRechitClusterZ);
   fChain->SetBranchAddress("dtRechitClusterTime", dtRechitClusterTime, &b_dtRechitClusterTime);
   fChain->SetBranchAddress("dtRechitClusterTimeWire", dtRechitClusterTimeWire, &b_dtRechitClusterTimeWire);
   fChain->SetBranchAddress("dtRechitClusterTimeWirePruned", dtRechitClusterTimeWirePruned, &b_dtRechitClusterTimeWirePruned);
   fChain->SetBranchAddress("dtRechitClusterTimeTotal", dtRechitClusterTimeTotal, &b_dtRechitClusterTimeTotal);
   fChain->SetBranchAddress("dtRechitClusterWheel", dtRechitClusterWheel, &b_dtRechitClusterWheel);
   fChain->SetBranchAddress("dtRechitClusterTimeSpread", dtRechitClusterTimeSpread, &b_dtRechitClusterTimeSpread);
   fChain->SetBranchAddress("dtRechitClusterTimeTotalSpread", dtRechitClusterTimeTotalSpread, &b_dtRechitClusterTimeTotalSpread);
   fChain->SetBranchAddress("dtRechitClusterTimeTotalSpreadPruned", dtRechitClusterTimeTotalSpreadPruned, &b_dtRechitClusterTimeTotalSpreadPruned);
   fChain->SetBranchAddress("dtRechitClusterTimeWireSpread", dtRechitClusterTimeWireSpread, &b_dtRechitClusterTimeWireSpread);
   fChain->SetBranchAddress("dtRechitClusterGenMuonDeltaR", dtRechitClusterGenMuonDeltaR, &b_dtRechitClusterGenMuonDeltaR);
   fChain->SetBranchAddress("dtRechitClusterXYSpread", dtRechitClusterXYSpread, &b_dtRechitClusterXYSpread);
   fChain->SetBranchAddress("dtRechitClusterMajorAxis", dtRechitClusterMajorAxis, &b_dtRechitClusterMajorAxis);
   fChain->SetBranchAddress("dtRechitClusterMinorAxis", dtRechitClusterMinorAxis, &b_dtRechitClusterMinorAxis);
   fChain->SetBranchAddress("dtRechitClusterEtaPhiSpread", dtRechitClusterEtaPhiSpread, &b_dtRechitClusterEtaPhiSpread);
   fChain->SetBranchAddress("dtRechitClusterPhiSpread", dtRechitClusterPhiSpread, &b_dtRechitClusterPhiSpread);
   fChain->SetBranchAddress("dtRechitClusterEtaSpread", dtRechitClusterEtaSpread, &b_dtRechitClusterEtaSpread);
   fChain->SetBranchAddress("dtRechitClusterDeltaRSpread", dtRechitClusterDeltaRSpread, &b_dtRechitClusterDeltaRSpread);
   fChain->SetBranchAddress("dtRechitClusterXSpread", dtRechitClusterXSpread, &b_dtRechitClusterXSpread);
   fChain->SetBranchAddress("dtRechitClusterRSpread", dtRechitClusterRSpread, &b_dtRechitClusterRSpread);
   fChain->SetBranchAddress("dtRechitClusterYSpread", dtRechitClusterYSpread, &b_dtRechitClusterYSpread);
   fChain->SetBranchAddress("dtRechitClusterZSpread", dtRechitClusterZSpread, &b_dtRechitClusterZSpread);
   fChain->SetBranchAddress("dtRechitClusterPhi", dtRechitClusterPhi, &b_dtRechitClusterPhi);
   fChain->SetBranchAddress("dtRechitClusterEta", dtRechitClusterEta, &b_dtRechitClusterEta);
   fChain->SetBranchAddress("dtRechitClusterJetVetoPt", dtRechitClusterJetVetoPt, &b_dtRechitClusterJetVetoPt);
   fChain->SetBranchAddress("dtRechitClusterJetVetoPtJESUp", dtRechitClusterJetVetoPtJESUp, &b_dtRechitClusterJetVetoPtJESUp);
   fChain->SetBranchAddress("dtRechitClusterJetVetoPtJESDown", dtRechitClusterJetVetoPtJESDown, &b_dtRechitClusterJetVetoPtJESDown);
   fChain->SetBranchAddress("dtRechitClusterJetVetoEta", dtRechitClusterJetVetoEta, &b_dtRechitClusterJetVetoEta);
   fChain->SetBranchAddress("dtRechitClusterJetVetoPhi", dtRechitClusterJetVetoPhi, &b_dtRechitClusterJetVetoPhi);
   fChain->SetBranchAddress("dtRechitClusterTightJetVetoPt", dtRechitClusterTightJetVetoPt, &b_dtRechitClusterTightJetVetoPt);
   fChain->SetBranchAddress("dtRechitClusterTightJetVetoPtJESUp", dtRechitClusterTightJetVetoPtJESUp, &b_dtRechitClusterTightJetVetoPtJESUp);
   fChain->SetBranchAddress("dtRechitClusterTightJetVetoPtJESDown", dtRechitClusterTightJetVetoPtJESDown, &b_dtRechitClusterTightJetVetoPtJESDown);
   fChain->SetBranchAddress("dtRechitClusterTightJetVetoEta", dtRechitClusterTightJetVetoEta, &b_dtRechitClusterTightJetVetoEta);
   fChain->SetBranchAddress("dtRechitClusterTightJetVetoPhi", dtRechitClusterTightJetVetoPhi, &b_dtRechitClusterTightJetVetoPhi);
   fChain->SetBranchAddress("dtRechitClusterJetVetoE", dtRechitClusterJetVetoE, &b_dtRechitClusterJetVetoE);
   fChain->SetBranchAddress("dtRechitClusterGenJetVetoPt", dtRechitClusterGenJetVetoPt, &b_dtRechitClusterGenJetVetoPt);
   fChain->SetBranchAddress("dtRechitClusterGenJetVetoE", dtRechitClusterGenJetVetoE, &b_dtRechitClusterGenJetVetoE);
   fChain->SetBranchAddress("dtRechitClusterMuonVetoPt", dtRechitClusterMuonVetoPt, &b_dtRechitClusterMuonVetoPt);
   fChain->SetBranchAddress("dtRechitClusterMuonVetoE", dtRechitClusterMuonVetoE, &b_dtRechitClusterMuonVetoE);
   fChain->SetBranchAddress("dtRechitClusterMuonVetoPhi", dtRechitClusterMuonVetoPhi, &b_dtRechitClusterMuonVetoPhi);
   fChain->SetBranchAddress("dtRechitClusterMuonVetoEta", dtRechitClusterMuonVetoEta, &b_dtRechitClusterMuonVetoEta);
   fChain->SetBranchAddress("dtRechitClusterJetVetoElectronEnergyFraction", dtRechitClusterJetVetoElectronEnergyFraction, &b_dtRechitClusterJetVetoElectronEnergyFraction);
   fChain->SetBranchAddress("dtRechitClusterJetVetoPhotonEnergyFraction", dtRechitClusterJetVetoPhotonEnergyFraction, &b_dtRechitClusterJetVetoPhotonEnergyFraction);
   fChain->SetBranchAddress("dtRechitClusterJetVetoChargedHadronEnergyFraction", dtRechitClusterJetVetoChargedHadronEnergyFraction, &b_dtRechitClusterJetVetoChargedHadronEnergyFraction);
   fChain->SetBranchAddress("dtRechitClusterJetVetoNeutralHadronEnergyFraction", dtRechitClusterJetVetoNeutralHadronEnergyFraction, &b_dtRechitClusterJetVetoNeutralHadronEnergyFraction);
   fChain->SetBranchAddress("dtRechitClusterJetVetoMuonEnergyFraction", dtRechitClusterJetVetoMuonEnergyFraction, &b_dtRechitClusterJetVetoMuonEnergyFraction);
   fChain->SetBranchAddress("dtRechitClusterJetVetoPt_0p6", dtRechitClusterJetVetoPt_0p6, &b_dtRechitClusterJetVetoPt_0p6);
   fChain->SetBranchAddress("dtRechitClusterJetVetoPt_0p8", dtRechitClusterJetVetoPt_0p8, &b_dtRechitClusterJetVetoPt_0p8);
   fChain->SetBranchAddress("dtRechitClusterJetVetoE_0p6", dtRechitClusterJetVetoE_0p6, &b_dtRechitClusterJetVetoE_0p6);
   fChain->SetBranchAddress("dtRechitClusterJetVetoE_0p8", dtRechitClusterJetVetoE_0p8, &b_dtRechitClusterJetVetoE_0p8);
   fChain->SetBranchAddress("dtRechitClusterMuonVetoPt_0p6", dtRechitClusterMuonVetoPt_0p6, &b_dtRechitClusterMuonVetoPt_0p6);
   fChain->SetBranchAddress("dtRechitClusterMuonVetoPt_0p8", dtRechitClusterMuonVetoPt_0p8, &b_dtRechitClusterMuonVetoPt_0p8);
   fChain->SetBranchAddress("dtRechitClusterMuonVetoE_0p6", dtRechitClusterMuonVetoE_0p6, &b_dtRechitClusterMuonVetoE_0p6);
   fChain->SetBranchAddress("dtRechitClusterMuonVetoE_0p8", dtRechitClusterMuonVetoE_0p8, &b_dtRechitClusterMuonVetoE_0p8);
   fChain->SetBranchAddress("dtRechitClusterMuonVetoLooseIso", dtRechitClusterMuonVetoLooseIso, &b_dtRechitClusterMuonVetoLooseIso);
   fChain->SetBranchAddress("dtRechitClusterMuonVetoTightIso", dtRechitClusterMuonVetoTightIso, &b_dtRechitClusterMuonVetoTightIso);
   fChain->SetBranchAddress("dtRechitClusterMuonVetoVTightIso", dtRechitClusterMuonVetoVTightIso, &b_dtRechitClusterMuonVetoVTightIso);
   fChain->SetBranchAddress("dtRechitClusterMuonVetoVVTightIso", dtRechitClusterMuonVetoVVTightIso, &b_dtRechitClusterMuonVetoVVTightIso);
   fChain->SetBranchAddress("dtRechitClusterMuonVetoTightId", dtRechitClusterMuonVetoTightId, &b_dtRechitClusterMuonVetoTightId);
   fChain->SetBranchAddress("dtRechitClusterMuonVetoLooseId", dtRechitClusterMuonVetoLooseId, &b_dtRechitClusterMuonVetoLooseId);
   fChain->SetBranchAddress("dtRechitClusterMuonVetoGlobal", dtRechitClusterMuonVetoGlobal, &b_dtRechitClusterMuonVetoGlobal);
   fChain->SetBranchAddress("dtRechitClusterMuonVetoIso", dtRechitClusterMuonVetoIso, &b_dtRechitClusterMuonVetoIso);
   fChain->SetBranchAddress("dtRechitClusterIsoMuonVetoPt", dtRechitClusterIsoMuonVetoPt, &b_dtRechitClusterIsoMuonVetoPt);
   fChain->SetBranchAddress("dtRechitClusterIsoMuonVetoE", dtRechitClusterIsoMuonVetoE, &b_dtRechitClusterIsoMuonVetoE);
   fChain->SetBranchAddress("dtRechitClusterIsoMuonVetoPhi", dtRechitClusterIsoMuonVetoPhi, &b_dtRechitClusterIsoMuonVetoPhi);
   fChain->SetBranchAddress("dtRechitClusterIsoMuonVetoEta", dtRechitClusterIsoMuonVetoEta, &b_dtRechitClusterIsoMuonVetoEta);
   fChain->SetBranchAddress("dtRechitClusterGenMuonVetoPt", dtRechitClusterGenMuonVetoPt, &b_dtRechitClusterGenMuonVetoPt);
   fChain->SetBranchAddress("dtRechitClusterMuonVetoType", dtRechitClusterMuonVetoType, &b_dtRechitClusterMuonVetoType);
   fChain->SetBranchAddress("dtRechitClusterZLep1", dtRechitClusterZLep1, &b_dtRechitClusterZLep1);
   fChain->SetBranchAddress("dtRechitClusterZLep2", dtRechitClusterZLep2, &b_dtRechitClusterZLep2);
   fChain->SetBranchAddress("dtRechitClusterZLep1Tag", dtRechitClusterZLep1Tag, &b_dtRechitClusterZLep1Tag);
   fChain->SetBranchAddress("dtRechitClusterZLep2Tag", dtRechitClusterZLep2Tag, &b_dtRechitClusterZLep2Tag);
   fChain->SetBranchAddress("dtRechitClusterZLep1Id", dtRechitClusterZLep1Id, &b_dtRechitClusterZLep1Id);
   fChain->SetBranchAddress("dtRechitClusterZLep2Id", dtRechitClusterZLep2Id, &b_dtRechitClusterZLep2Id);
   fChain->SetBranchAddress("dtRechitCluster2ZLep1LooseIso", dtRechitCluster2ZLep1LooseIso, &b_dtRechitCluster2ZLep1LooseIso);
   fChain->SetBranchAddress("dtRechitClusterZLep1TightIso", dtRechitClusterZLep1TightIso, &b_dtRechitClusterZLep1TightIso);
   fChain->SetBranchAddress("dtRechitClusterZLep1VTightIso", dtRechitClusterZLep1VTightIso, &b_dtRechitClusterZLep1VTightIso);
   fChain->SetBranchAddress("dtRechitClusterZLep1VVTightIso", dtRechitClusterZLep1VVTightIso, &b_dtRechitClusterZLep1VVTightIso);
   fChain->SetBranchAddress("dtRechitClusterZLep1TightId", dtRechitClusterZLep1TightId, &b_dtRechitClusterZLep1TightId);
   fChain->SetBranchAddress("dtRechitClusterZLep2LooseIso", dtRechitClusterZLep2LooseIso, &b_dtRechitClusterZLep2LooseIso);
   fChain->SetBranchAddress("dtRechitClusterZLep2TightIso", dtRechitClusterZLep2TightIso, &b_dtRechitClusterZLep2TightIso);
   fChain->SetBranchAddress("dtRechitClusterZLep2VTightIso", dtRechitClusterZLep2VTightIso, &b_dtRechitClusterZLep2VTightIso);
   fChain->SetBranchAddress("dtRechitClusterZLep2VVTightIso", dtRechitClusterZLep2VVTightIso, &b_dtRechitClusterZLep2VVTightIso);
   fChain->SetBranchAddress("dtRechitClusterZLep2TightId", dtRechitClusterZLep2TightId, &b_dtRechitClusterZLep2TightId);
   fChain->SetBranchAddress("dtRechitClusterSize", dtRechitClusterSize, &b_dtRechitClusterSize);
   fChain->SetBranchAddress("dtRechitClusterNoiseHit", dtRechitClusterNoiseHit, &b_dtRechitClusterNoiseHit);
   fChain->SetBranchAddress("dtRechitClusterNoiseHitStation1", dtRechitClusterNoiseHitStation1, &b_dtRechitClusterNoiseHitStation1);
   fChain->SetBranchAddress("dtRechitClusterNoiseHitStation2", dtRechitClusterNoiseHitStation2, &b_dtRechitClusterNoiseHitStation2);
   fChain->SetBranchAddress("dtRechitClusterNoiseHitStation3", dtRechitClusterNoiseHitStation3, &b_dtRechitClusterNoiseHitStation3);
   fChain->SetBranchAddress("dtRechitClusterNoiseHitStation4", dtRechitClusterNoiseHitStation4, &b_dtRechitClusterNoiseHitStation4);
   fChain->SetBranchAddress("dtRechitClusterMe11Ratio", dtRechitClusterMe11Ratio, &b_dtRechitClusterMe11Ratio);
   fChain->SetBranchAddress("dtRechitClusterMe12Ratio", dtRechitClusterMe12Ratio, &b_dtRechitClusterMe12Ratio);
   fChain->SetBranchAddress("dtRechitClusterNStation", dtRechitClusterNStation, &b_dtRechitClusterNStation);
   fChain->SetBranchAddress("dtRechitClusterNStation5", dtRechitClusterNStation5, &b_dtRechitClusterNStation5);
   fChain->SetBranchAddress("dtRechitClusterNStation10", dtRechitClusterNStation10, &b_dtRechitClusterNStation10);
   fChain->SetBranchAddress("dtRechitClusterNStation10perc", dtRechitClusterNStation10perc, &b_dtRechitClusterNStation10perc);
   fChain->SetBranchAddress("dtRechitClusterAvgStation", dtRechitClusterAvgStation, &b_dtRechitClusterAvgStation);
   fChain->SetBranchAddress("dtRechitClusterAvgStation5", dtRechitClusterAvgStation5, &b_dtRechitClusterAvgStation5);
   fChain->SetBranchAddress("dtRechitClusterAvgStation10", dtRechitClusterAvgStation10, &b_dtRechitClusterAvgStation10);
   fChain->SetBranchAddress("dtRechitClusterAvgStation10perc", dtRechitClusterAvgStation10perc, &b_dtRechitClusterAvgStation10perc);
   fChain->SetBranchAddress("dtRechitClusterMaxStation", dtRechitClusterMaxStation, &b_dtRechitClusterMaxStation);
   fChain->SetBranchAddress("dtRechitClusterMaxStationRatio", dtRechitClusterMaxStationRatio, &b_dtRechitClusterMaxStationRatio);
   fChain->SetBranchAddress("dtRechitClusterNChamber", dtRechitClusterNChamber, &b_dtRechitClusterNChamber);
   fChain->SetBranchAddress("dtRechitClusterMaxChamber", dtRechitClusterMaxChamber, &b_dtRechitClusterMaxChamber);
   fChain->SetBranchAddress("dtRechitClusterMaxChamberRatio", dtRechitClusterMaxChamberRatio, &b_dtRechitClusterMaxChamberRatio);
   fChain->SetBranchAddress("dtRechitClusterNSegmentStation1", dtRechitClusterNSegmentStation1, &b_dtRechitClusterNSegmentStation1);
   fChain->SetBranchAddress("dtRechitClusterNSegmentStation2", dtRechitClusterNSegmentStation2, &b_dtRechitClusterNSegmentStation2);
   fChain->SetBranchAddress("dtRechitClusterNSegmentStation3", dtRechitClusterNSegmentStation3, &b_dtRechitClusterNSegmentStation3);
   fChain->SetBranchAddress("dtRechitClusterNSegmentStation4", dtRechitClusterNSegmentStation4, &b_dtRechitClusterNSegmentStation4);
   fChain->SetBranchAddress("dtRechitClusterNRechitChamberPlus11", dtRechitClusterNRechitChamberPlus11, &b_dtRechitClusterNRechitChamberPlus11);
   fChain->SetBranchAddress("dtRechitClusterNRechitChamberPlus12", dtRechitClusterNRechitChamberPlus12, &b_dtRechitClusterNRechitChamberPlus12);
   fChain->SetBranchAddress("dtRechitClusterNRechitChamberPlus13", dtRechitClusterNRechitChamberPlus13, &b_dtRechitClusterNRechitChamberPlus13);
   fChain->SetBranchAddress("dtRechitClusterNRechitChamberPlus21", dtRechitClusterNRechitChamberPlus21, &b_dtRechitClusterNRechitChamberPlus21);
   fChain->SetBranchAddress("dtRechitClusterNRechitChamberPlus22", dtRechitClusterNRechitChamberPlus22, &b_dtRechitClusterNRechitChamberPlus22);
   fChain->SetBranchAddress("dtRechitClusterNRechitChamberPlus31", dtRechitClusterNRechitChamberPlus31, &b_dtRechitClusterNRechitChamberPlus31);
   fChain->SetBranchAddress("dtRechitClusterNRechitChamberPlus32", dtRechitClusterNRechitChamberPlus32, &b_dtRechitClusterNRechitChamberPlus32);
   fChain->SetBranchAddress("dtRechitClusterNRechitChamberPlus41", dtRechitClusterNRechitChamberPlus41, &b_dtRechitClusterNRechitChamberPlus41);
   fChain->SetBranchAddress("dtRechitClusterNRechitChamberPlus42", dtRechitClusterNRechitChamberPlus42, &b_dtRechitClusterNRechitChamberPlus42);
   fChain->SetBranchAddress("dtRechitClusterNRechitChamberMinus11", dtRechitClusterNRechitChamberMinus11, &b_dtRechitClusterNRechitChamberMinus11);
   fChain->SetBranchAddress("dtRechitClusterNRechitChamberMinus12", dtRechitClusterNRechitChamberMinus12, &b_dtRechitClusterNRechitChamberMinus12);
   fChain->SetBranchAddress("dtRechitClusterNRechitChamberMinus13", dtRechitClusterNRechitChamberMinus13, &b_dtRechitClusterNRechitChamberMinus13);
   fChain->SetBranchAddress("dtRechitClusterNRechitChamberMinus21", dtRechitClusterNRechitChamberMinus21, &b_dtRechitClusterNRechitChamberMinus21);
   fChain->SetBranchAddress("dtRechitClusterNRechitChamberMinus22", dtRechitClusterNRechitChamberMinus22, &b_dtRechitClusterNRechitChamberMinus22);
   fChain->SetBranchAddress("dtRechitClusterNRechitChamberMinus31", dtRechitClusterNRechitChamberMinus31, &b_dtRechitClusterNRechitChamberMinus31);
   fChain->SetBranchAddress("dtRechitClusterNRechitChamberMinus32", dtRechitClusterNRechitChamberMinus32, &b_dtRechitClusterNRechitChamberMinus32);
   fChain->SetBranchAddress("dtRechitClusterNRechitChamberMinus41", dtRechitClusterNRechitChamberMinus41, &b_dtRechitClusterNRechitChamberMinus41);
   fChain->SetBranchAddress("dtRechitClusterNRechitChamberMinus42", dtRechitClusterNRechitChamberMinus42, &b_dtRechitClusterNRechitChamberMinus42);
   fChain->SetBranchAddress("dtRechitClusterMet_dPhi", dtRechitClusterMet_dPhi, &b_dtRechitClusterMet_dPhi);
   fChain->SetBranchAddress("dtRechitClusterMetJESUp_dPhi", dtRechitClusterMetJESUp_dPhi, &b_dtRechitClusterMetJESUp_dPhi);
   fChain->SetBranchAddress("dtRechitClusterMetJESDown_dPhi", dtRechitClusterMetJESDown_dPhi, &b_dtRechitClusterMetJESDown_dPhi);
   fChain->SetBranchAddress("dtRechitClusterMetXYCorr_dPhi", dtRechitClusterMetXYCorr_dPhi, &b_dtRechitClusterMetXYCorr_dPhi);
   fChain->SetBranchAddress("dtRechitClusterNLayersChamberPlus11", dtRechitClusterNLayersChamberPlus11, &b_dtRechitClusterNLayersChamberPlus11);
   fChain->SetBranchAddress("dtRechitClusterNLayersChamberPlus12", dtRechitClusterNLayersChamberPlus12, &b_dtRechitClusterNLayersChamberPlus12);
   fChain->SetBranchAddress("dtRechitClusterNLayersChamberPlus13", dtRechitClusterNLayersChamberPlus13, &b_dtRechitClusterNLayersChamberPlus13);
   fChain->SetBranchAddress("dtRechitClusterNLayersChamberPlus21", dtRechitClusterNLayersChamberPlus21, &b_dtRechitClusterNLayersChamberPlus21);
   fChain->SetBranchAddress("dtRechitClusterNLayersChamberPlus22", dtRechitClusterNLayersChamberPlus22, &b_dtRechitClusterNLayersChamberPlus22);
   fChain->SetBranchAddress("dtRechitClusterNLayersChamberPlus31", dtRechitClusterNLayersChamberPlus31, &b_dtRechitClusterNLayersChamberPlus31);
   fChain->SetBranchAddress("dtRechitClusterNLayersChamberPlus32", dtRechitClusterNLayersChamberPlus32, &b_dtRechitClusterNLayersChamberPlus32);
   fChain->SetBranchAddress("dtRechitClusterNLayersChamberPlus41", dtRechitClusterNLayersChamberPlus41, &b_dtRechitClusterNLayersChamberPlus41);
   fChain->SetBranchAddress("dtRechitClusterNLayersChamberPlus42", dtRechitClusterNLayersChamberPlus42, &b_dtRechitClusterNLayersChamberPlus42);
   fChain->SetBranchAddress("dtRechitClusterNLayersChamberMinus11", dtRechitClusterNLayersChamberMinus11, &b_dtRechitClusterNLayersChamberMinus11);
   fChain->SetBranchAddress("dtRechitClusterNLayersChamberMinus12", dtRechitClusterNLayersChamberMinus12, &b_dtRechitClusterNLayersChamberMinus12);
   fChain->SetBranchAddress("dtRechitClusterNLayersChamberMinus13", dtRechitClusterNLayersChamberMinus13, &b_dtRechitClusterNLayersChamberMinus13);
   fChain->SetBranchAddress("dtRechitClusterNLayersChamberMinus21", dtRechitClusterNLayersChamberMinus21, &b_dtRechitClusterNLayersChamberMinus21);
   fChain->SetBranchAddress("dtRechitClusterNLayersChamberMinus22", dtRechitClusterNLayersChamberMinus22, &b_dtRechitClusterNLayersChamberMinus22);
   fChain->SetBranchAddress("dtRechitClusterNLayersChamberMinus31", dtRechitClusterNLayersChamberMinus31, &b_dtRechitClusterNLayersChamberMinus31);
   fChain->SetBranchAddress("dtRechitClusterNLayersChamberMinus32", dtRechitClusterNLayersChamberMinus32, &b_dtRechitClusterNLayersChamberMinus32);
   fChain->SetBranchAddress("dtRechitClusterNLayersChamberMinus41", dtRechitClusterNLayersChamberMinus41, &b_dtRechitClusterNLayersChamberMinus41);
   fChain->SetBranchAddress("dtRechitClusterNLayersChamberMinus42", dtRechitClusterNLayersChamberMinus42, &b_dtRechitClusterNLayersChamberMinus42);
   fChain->SetBranchAddress("dtRechitClusterMetHEM_dPhi", dtRechitClusterMetHEM_dPhi, &b_dtRechitClusterMetHEM_dPhi);
   fChain->SetBranchAddress("dtRechitClusterMetHEMXYCorr_dPhi", dtRechitClusterMetHEMXYCorr_dPhi, &b_dtRechitClusterMetHEMXYCorr_dPhi);
   fChain->SetBranchAddress("dtRechitClusterMetEENoise_dPhi", dtRechitClusterMetEENoise_dPhi, &b_dtRechitClusterMetEENoise_dPhi);
   fChain->SetBranchAddress("dtRechitClusterMetEENoiseXYCorr_dPhi", dtRechitClusterMetEENoiseXYCorr_dPhi, &b_dtRechitClusterMetEENoiseXYCorr_dPhi);
   fChain->SetBranchAddress("dtRechitClusterMetJesUp_dPhi", dtRechitClusterMetJesUp_dPhi, &b_dtRechitClusterMetJesUp_dPhi);
   fChain->SetBranchAddress("dtRechitClusterMetJesDown_dPhi", dtRechitClusterMetJesDown_dPhi, &b_dtRechitClusterMetJesDown_dPhi);
   fChain->SetBranchAddress("dtRechitCluster_match_dtSeg_0p5", dtRechitCluster_match_dtSeg_0p5, &b_dtRechitCluster_match_dtSeg_0p5);
   fChain->SetBranchAddress("dtRechitCluster_match_dtSegTime_0p5", dtRechitCluster_match_dtSegTime_0p5, &b_dtRechitCluster_match_dtSegTime_0p5);
   fChain->SetBranchAddress("dtRechitCluster_match_dtSeg_0p4", dtRechitCluster_match_dtSeg_0p4, &b_dtRechitCluster_match_dtSeg_0p4);
   fChain->SetBranchAddress("dtRechitCluster_match_dtSegTime_0p4", dtRechitCluster_match_dtSegTime_0p4, &b_dtRechitCluster_match_dtSegTime_0p4);
   fChain->SetBranchAddress("dtRechitCluster_match_dtSegTimeSpread_0p5", dtRechitCluster_match_dtSegTimeSpread_0p5, &b_dtRechitCluster_match_dtSegTimeSpread_0p5);
   fChain->SetBranchAddress("dtRechitCluster_match_dtSegTimeSpread_0p4", dtRechitCluster_match_dtSegTimeSpread_0p4, &b_dtRechitCluster_match_dtSegTimeSpread_0p4);
   fChain->SetBranchAddress("dtRechitCluster_match_dtSeg_sameStation_0p5", dtRechitCluster_match_dtSeg_sameStation_0p5, &b_dtRechitCluster_match_dtSeg_sameStation_0p5);
   fChain->SetBranchAddress("dtRechitCluster_match_dtSegTime_sameStation_0p5", dtRechitCluster_match_dtSegTime_sameStation_0p5, &b_dtRechitCluster_match_dtSegTime_sameStation_0p5);
   fChain->SetBranchAddress("dtRechitCluster_match_dtSeg_sameStation_0p4", dtRechitCluster_match_dtSeg_sameStation_0p4, &b_dtRechitCluster_match_dtSeg_sameStation_0p4);
   fChain->SetBranchAddress("dtRechitCluster_match_dtSegTime_sameStation_0p4", dtRechitCluster_match_dtSegTime_sameStation_0p4, &b_dtRechitCluster_match_dtSegTime_sameStation_0p4);
   fChain->SetBranchAddress("dtRechitCluster_match_dtSegTimeSpread_sameStation_0p5", dtRechitCluster_match_dtSegTimeSpread_sameStation_0p5, &b_dtRechitCluster_match_dtSegTimeSpread_sameStation_0p5);
   fChain->SetBranchAddress("dtRechitCluster_match_dtSegTimeSpread_sameStation_0p4", dtRechitCluster_match_dtSegTimeSpread_sameStation_0p4, &b_dtRechitCluster_match_dtSegTimeSpread_sameStation_0p4);
   fChain->SetBranchAddress("dtRechitCluster_match_RPCTime_dPhi0p5", dtRechitCluster_match_RPCTime_dPhi0p5, &b_dtRechitCluster_match_RPCTime_dPhi0p5);
   fChain->SetBranchAddress("dtRechitCluster_match_RPCTimeSpread_dPhi0p5", dtRechitCluster_match_RPCTimeSpread_dPhi0p5, &b_dtRechitCluster_match_RPCTimeSpread_dPhi0p5);
   fChain->SetBranchAddress("dtRechitCluster_match_RPCTime_dR0p4", dtRechitCluster_match_RPCTime_dR0p4, &b_dtRechitCluster_match_RPCTime_dR0p4);
   fChain->SetBranchAddress("dtRechitCluster_match_RPCTimeSpread_dR0p4", dtRechitCluster_match_RPCTimeSpread_dR0p4, &b_dtRechitCluster_match_RPCTimeSpread_dR0p4);
   fChain->SetBranchAddress("dtRechitCluster_match_RPChits_dR0p4", dtRechitCluster_match_RPChits_dR0p4, &b_dtRechitCluster_match_RPChits_dR0p4);
   fChain->SetBranchAddress("dtRechitCluster_match_RPCTime_sameStation_dR0p4", dtRechitCluster_match_RPCTime_sameStation_dR0p4, &b_dtRechitCluster_match_RPCTime_sameStation_dR0p4);
   fChain->SetBranchAddress("dtRechitCluster_match_RPCTimeSpread_sameStation_dR0p4", dtRechitCluster_match_RPCTimeSpread_sameStation_dR0p4, &b_dtRechitCluster_match_RPCTimeSpread_sameStation_dR0p4);
   fChain->SetBranchAddress("dtRechitCluster_match_RPChits_sameStation_dR0p4", dtRechitCluster_match_RPChits_sameStation_dR0p4, &b_dtRechitCluster_match_RPChits_sameStation_dR0p4);
   fChain->SetBranchAddress("dtRechitCluster_match_gParticle_Id", dtRechitCluster_match_gParticle_Id, &b_dtRechitCluster_match_gParticle_Id);
   fChain->SetBranchAddress("dtRechitCluster_match_gParticle_Pt", dtRechitCluster_match_gParticle_Pt, &b_dtRechitCluster_match_gParticle_Pt);
   fChain->SetBranchAddress("dtRechitCluster_match_gParticle_Eta", dtRechitCluster_match_gParticle_Eta, &b_dtRechitCluster_match_gParticle_Eta);
   fChain->SetBranchAddress("dtRechitCluster_match_gParticle_Phi", dtRechitCluster_match_gParticle_Phi, &b_dtRechitCluster_match_gParticle_Phi);
   fChain->SetBranchAddress("dtRechitCluster_match_gParticle_E", dtRechitCluster_match_gParticle_E, &b_dtRechitCluster_match_gParticle_E);
   fChain->SetBranchAddress("dtRechitCluster_match_gParticle_Status", dtRechitCluster_match_gParticle_Status, &b_dtRechitCluster_match_gParticle_Status);
   fChain->SetBranchAddress("dtRechitCluster_match_gParticle_MotherId", dtRechitCluster_match_gParticle_MotherId, &b_dtRechitCluster_match_gParticle_MotherId);
   fChain->SetBranchAddress("dtRechitCluster_match_gParticle_deltaR", dtRechitCluster_match_gParticle_deltaR, &b_dtRechitCluster_match_gParticle_deltaR);
   fChain->SetBranchAddress("dtRechitCluster_match_RPChits_dPhi0p5", dtRechitCluster_match_RPChits_dPhi0p5, &b_dtRechitCluster_match_RPChits_dPhi0p5);
   fChain->SetBranchAddress("dtRechitCluster_match_RPCBx_dPhi0p5", dtRechitCluster_match_RPCBx_dPhi0p5, &b_dtRechitCluster_match_RPCBx_dPhi0p5);
   fChain->SetBranchAddress("dtRechitCluster_match_RB1_0p4", dtRechitCluster_match_RB1_0p4, &b_dtRechitCluster_match_RB1_0p4);
   fChain->SetBranchAddress("dtRechitCluster_match_RB1_dPhi0p5", dtRechitCluster_match_RB1_dPhi0p5, &b_dtRechitCluster_match_RB1_dPhi0p5);
   fChain->SetBranchAddress("dtRechitCluster_match_MB1Seg_0p4", dtRechitCluster_match_MB1Seg_0p4, &b_dtRechitCluster_match_MB1Seg_0p4);
   fChain->SetBranchAddress("dtRechitCluster_match_MB1Seg_0p5", dtRechitCluster_match_MB1Seg_0p5, &b_dtRechitCluster_match_MB1Seg_0p5);
   fChain->SetBranchAddress("dtRechitCluster_match_MB1hits_0p4", dtRechitCluster_match_MB1hits_0p4, &b_dtRechitCluster_match_MB1hits_0p4);
   fChain->SetBranchAddress("dtRechitCluster_match_MB1hits_0p5", dtRechitCluster_match_MB1hits_0p5, &b_dtRechitCluster_match_MB1hits_0p5);
   fChain->SetBranchAddress("dtRechitCluster_match_MB1hits_cosmics_plus", dtRechitCluster_match_MB1hits_cosmics_plus, &b_dtRechitCluster_match_MB1hits_cosmics_plus);
   fChain->SetBranchAddress("dtRechitCluster_match_MB1hits_cosmics_minus", dtRechitCluster_match_MB1hits_cosmics_minus, &b_dtRechitCluster_match_MB1hits_cosmics_minus);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP", &dtRechitCluster2_match_gLLP, &b_dtRechitCluster2_match_gLLP);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_minDeltaR", &dtRechitCluster2_match_gLLP_minDeltaR, &b_dtRechitCluster2_match_gLLP_minDeltaR);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_index", &dtRechitCluster2_match_gLLP_index, &b_dtRechitCluster2_match_gLLP_index);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_eta", &dtRechitCluster2_match_gLLP_eta, &b_dtRechitCluster2_match_gLLP_eta);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_phi", &dtRechitCluster2_match_gLLP_phi, &b_dtRechitCluster2_match_gLLP_phi);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_decay_r", &dtRechitCluster2_match_gLLP_decay_r, &b_dtRechitCluster2_match_gLLP_decay_r);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_decay_x", &dtRechitCluster2_match_gLLP_decay_x, &b_dtRechitCluster2_match_gLLP_decay_x);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_decay_y", &dtRechitCluster2_match_gLLP_decay_y, &b_dtRechitCluster2_match_gLLP_decay_y);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_decay_z", &dtRechitCluster2_match_gLLP_decay_z, &b_dtRechitCluster2_match_gLLP_decay_z);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_ctau", &dtRechitCluster2_match_gLLP_ctau, &b_dtRechitCluster2_match_gLLP_ctau);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_beta", &dtRechitCluster2_match_gLLP_beta, &b_dtRechitCluster2_match_gLLP_beta);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_csc", &dtRechitCluster2_match_gLLP_csc, &b_dtRechitCluster2_match_gLLP_csc);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_dt", &dtRechitCluster2_match_gLLP_dt, &b_dtRechitCluster2_match_gLLP_dt);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_multiplicity", &dtRechitCluster2_match_gLLP_multiplicity, &b_dtRechitCluster2_match_gLLP_multiplicity);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_EM_multiplicity", &dtRechitCluster2_match_gLLP_EM_multiplicity, &b_dtRechitCluster2_match_gLLP_EM_multiplicity);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_daughterKaon", &dtRechitCluster2_match_gLLP_daughterKaon, &b_dtRechitCluster2_match_gLLP_daughterKaon);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_e", &dtRechitCluster2_match_gLLP_e, &b_dtRechitCluster2_match_gLLP_e);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_pt", &dtRechitCluster2_match_gLLP_pt, &b_dtRechitCluster2_match_gLLP_pt);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_EMFracE", &dtRechitCluster2_match_gLLP_EMFracE, &b_dtRechitCluster2_match_gLLP_EMFracE);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_EMFracEz", &dtRechitCluster2_match_gLLP_EMFracEz, &b_dtRechitCluster2_match_gLLP_EMFracEz);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_EMFracP", &dtRechitCluster2_match_gLLP_EMFracP, &b_dtRechitCluster2_match_gLLP_EMFracP);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_EMFracPz", &dtRechitCluster2_match_gLLP_EMFracPz, &b_dtRechitCluster2_match_gLLP_EMFracPz);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_visE", &dtRechitCluster2_match_gLLP_visE, &b_dtRechitCluster2_match_gLLP_visE);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_visEz", &dtRechitCluster2_match_gLLP_visEz, &b_dtRechitCluster2_match_gLLP_visEz);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_visP", &dtRechitCluster2_match_gLLP_visP, &b_dtRechitCluster2_match_gLLP_visP);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_visPz", &dtRechitCluster2_match_gLLP_visPz, &b_dtRechitCluster2_match_gLLP_visPz);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_lepdPhi", &dtRechitCluster2_match_gLLP_lepdPhi, &b_dtRechitCluster2_match_gLLP_lepdPhi);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_daughter0_deltaR", &dtRechitCluster2_match_gLLP_daughter0_deltaR, &b_dtRechitCluster2_match_gLLP_daughter0_deltaR);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_daughter1_deltaR", &dtRechitCluster2_match_gLLP_daughter1_deltaR, &b_dtRechitCluster2_match_gLLP_daughter1_deltaR);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_daughter2_deltaR", &dtRechitCluster2_match_gLLP_daughter2_deltaR, &b_dtRechitCluster2_match_gLLP_daughter2_deltaR);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_daughter3_deltaR", &dtRechitCluster2_match_gLLP_daughter3_deltaR, &b_dtRechitCluster2_match_gLLP_daughter3_deltaR);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_other_daughter_deltaR", &dtRechitCluster2_match_gLLP_other_daughter_deltaR, &b_dtRechitCluster2_match_gLLP_other_daughter_deltaR);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_other_daughter_index", &dtRechitCluster2_match_gLLP_other_daughter_index, &b_dtRechitCluster2_match_gLLP_other_daughter_index);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_other_eta", &dtRechitCluster2_match_gLLP_other_eta, &b_dtRechitCluster2_match_gLLP_other_eta);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_other_phi", &dtRechitCluster2_match_gLLP_other_phi, &b_dtRechitCluster2_match_gLLP_other_phi);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_other_decay_r", &dtRechitCluster2_match_gLLP_other_decay_r, &b_dtRechitCluster2_match_gLLP_other_decay_r);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_other_decay_x", &dtRechitCluster2_match_gLLP_other_decay_x, &b_dtRechitCluster2_match_gLLP_other_decay_x);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_other_decay_y", &dtRechitCluster2_match_gLLP_other_decay_y, &b_dtRechitCluster2_match_gLLP_other_decay_y);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_other_decay_z", &dtRechitCluster2_match_gLLP_other_decay_z, &b_dtRechitCluster2_match_gLLP_other_decay_z);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_other_ctau", &dtRechitCluster2_match_gLLP_other_ctau, &b_dtRechitCluster2_match_gLLP_other_ctau);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_other_beta", &dtRechitCluster2_match_gLLP_other_beta, &b_dtRechitCluster2_match_gLLP_other_beta);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_other_csc", &dtRechitCluster2_match_gLLP_other_csc, &b_dtRechitCluster2_match_gLLP_other_csc);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_other_e", &dtRechitCluster2_match_gLLP_other_e, &b_dtRechitCluster2_match_gLLP_other_e);
   fChain->SetBranchAddress("dtRechitCluster2_match_gLLP_other_pt", &dtRechitCluster2_match_gLLP_other_pt, &b_dtRechitCluster2_match_gLLP_other_pt);
   fChain->SetBranchAddress("dtRechitCluster2X", &dtRechitCluster2X, &b_dtRechitCluster2X);
   fChain->SetBranchAddress("dtRechitCluster2Y", &dtRechitCluster2Y, &b_dtRechitCluster2Y);
   fChain->SetBranchAddress("dtRechitCluster2Z", &dtRechitCluster2Z, &b_dtRechitCluster2Z);
   fChain->SetBranchAddress("dtRechitCluster2Time", &dtRechitCluster2Time, &b_dtRechitCluster2Time);
   fChain->SetBranchAddress("dtRechitCluster2TimeWire", &dtRechitCluster2TimeWire, &b_dtRechitCluster2TimeWire);
   fChain->SetBranchAddress("dtRechitCluster2TimeWirePruned", &dtRechitCluster2TimeWirePruned, &b_dtRechitCluster2TimeWirePruned);
   fChain->SetBranchAddress("dtRechitCluster2TimeTotal", &dtRechitCluster2TimeTotal, &b_dtRechitCluster2TimeTotal);
   fChain->SetBranchAddress("dtRechitCluster2Wheel", &dtRechitCluster2Wheel, &b_dtRechitCluster2Wheel);
   fChain->SetBranchAddress("dtRechitCluster2TimeSpread", &dtRechitCluster2TimeSpread, &b_dtRechitCluster2TimeSpread);
   fChain->SetBranchAddress("dtRechitCluster2TimeTotalSpread", &dtRechitCluster2TimeTotalSpread, &b_dtRechitCluster2TimeTotalSpread);
   fChain->SetBranchAddress("dtRechitCluster2TimeTotalSpreadPruned", &dtRechitCluster2TimeTotalSpreadPruned, &b_dtRechitCluster2TimeTotalSpreadPruned);
   fChain->SetBranchAddress("dtRechitCluster2TimeWireSpread", &dtRechitCluster2TimeWireSpread, &b_dtRechitCluster2TimeWireSpread);
   fChain->SetBranchAddress("dtRechitCluster2GenMuonDeltaR", &dtRechitCluster2GenMuonDeltaR, &b_dtRechitCluster2GenMuonDeltaR);
   fChain->SetBranchAddress("dtRechitCluster2XYSpread", &dtRechitCluster2XYSpread, &b_dtRechitCluster2XYSpread);
   fChain->SetBranchAddress("dtRechitCluster2MajorAxis", &dtRechitCluster2MajorAxis, &b_dtRechitCluster2MajorAxis);
   fChain->SetBranchAddress("dtRechitCluster2MinorAxis", &dtRechitCluster2MinorAxis, &b_dtRechitCluster2MinorAxis);
   fChain->SetBranchAddress("dtRechitCluster2EtaPhiSpread", &dtRechitCluster2EtaPhiSpread, &b_dtRechitCluster2EtaPhiSpread);
   fChain->SetBranchAddress("dtRechitCluster2PhiSpread", &dtRechitCluster2PhiSpread, &b_dtRechitCluster2PhiSpread);
   fChain->SetBranchAddress("dtRechitCluster2EtaSpread", &dtRechitCluster2EtaSpread, &b_dtRechitCluster2EtaSpread);
   fChain->SetBranchAddress("dtRechitCluster2DeltaRSpread", &dtRechitCluster2DeltaRSpread, &b_dtRechitCluster2DeltaRSpread);
   fChain->SetBranchAddress("dtRechitCluster2XSpread", &dtRechitCluster2XSpread, &b_dtRechitCluster2XSpread);
   fChain->SetBranchAddress("dtRechitCluster2RSpread", &dtRechitCluster2RSpread, &b_dtRechitCluster2RSpread);
   fChain->SetBranchAddress("dtRechitCluster2YSpread", &dtRechitCluster2YSpread, &b_dtRechitCluster2YSpread);
   fChain->SetBranchAddress("dtRechitCluster2ZSpread", &dtRechitCluster2ZSpread, &b_dtRechitCluster2ZSpread);
   fChain->SetBranchAddress("dtRechitCluster2Phi", &dtRechitCluster2Phi, &b_dtRechitCluster2Phi);
   fChain->SetBranchAddress("dtRechitCluster2Eta", &dtRechitCluster2Eta, &b_dtRechitCluster2Eta);
   fChain->SetBranchAddress("dtRechitCluster2JetVetoPt", &dtRechitCluster2JetVetoPt, &b_dtRechitCluster2JetVetoPt);
   fChain->SetBranchAddress("dtRechitCluster2JetVetoEta", &dtRechitCluster2JetVetoEta, &b_dtRechitCluster2JetVetoEta);
   fChain->SetBranchAddress("dtRechitCluster2JetVetoPhi", &dtRechitCluster2JetVetoPhi, &b_dtRechitCluster2JetVetoPhi);
   fChain->SetBranchAddress("dtRechitCluster2JetVetoE", &dtRechitCluster2JetVetoE, &b_dtRechitCluster2JetVetoE);
   fChain->SetBranchAddress("dtRechitCluster2GenJetVetoPt", &dtRechitCluster2GenJetVetoPt, &b_dtRechitCluster2GenJetVetoPt);
   fChain->SetBranchAddress("dtRechitCluster2GenJetVetoE", &dtRechitCluster2GenJetVetoE, &b_dtRechitCluster2GenJetVetoE);
   fChain->SetBranchAddress("dtRechitCluster2MuonVetoPt", &dtRechitCluster2MuonVetoPt, &b_dtRechitCluster2MuonVetoPt);
   fChain->SetBranchAddress("dtRechitCluster2MuonVetoE", &dtRechitCluster2MuonVetoE, &b_dtRechitCluster2MuonVetoE);
   fChain->SetBranchAddress("dtRechitCluster2MuonVetoPhi", &dtRechitCluster2MuonVetoPhi, &b_dtRechitCluster2MuonVetoPhi);
   fChain->SetBranchAddress("dtRechitCluster2MuonVetoEta", &dtRechitCluster2MuonVetoEta, &b_dtRechitCluster2MuonVetoEta);
   fChain->SetBranchAddress("dtRechitCluster2JetVetoElectronEnergyFraction", &dtRechitCluster2JetVetoElectronEnergyFraction, &b_dtRechitCluster2JetVetoElectronEnergyFraction);
   fChain->SetBranchAddress("dtRechitCluster2JetVetoPhotonEnergyFraction", &dtRechitCluster2JetVetoPhotonEnergyFraction, &b_dtRechitCluster2JetVetoPhotonEnergyFraction);
   fChain->SetBranchAddress("dtRechitCluster2JetVetoChargedHadronEnergyFraction", &dtRechitCluster2JetVetoChargedHadronEnergyFraction, &b_dtRechitCluster2JetVetoChargedHadronEnergyFraction);
   fChain->SetBranchAddress("dtRechitCluster2JetVetoNeutralHadronEnergyFraction", &dtRechitCluster2JetVetoNeutralHadronEnergyFraction, &b_dtRechitCluster2JetVetoNeutralHadronEnergyFraction);
   fChain->SetBranchAddress("dtRechitCluster2JetVetoMuonEnergyFraction", &dtRechitCluster2JetVetoMuonEnergyFraction, &b_dtRechitCluster2JetVetoMuonEnergyFraction);
   fChain->SetBranchAddress("dtRechitCluster2JetVetoPt_0p6", &dtRechitCluster2JetVetoPt_0p6, &b_dtRechitCluster2JetVetoPt_0p6);
   fChain->SetBranchAddress("dtRechitCluster2JetVetoPt_0p8", &dtRechitCluster2JetVetoPt_0p8, &b_dtRechitCluster2JetVetoPt_0p8);
   fChain->SetBranchAddress("dtRechitCluster2JetVetoE_0p6", &dtRechitCluster2JetVetoE_0p6, &b_dtRechitCluster2JetVetoE_0p6);
   fChain->SetBranchAddress("dtRechitCluster2JetVetoE_0p8", &dtRechitCluster2JetVetoE_0p8, &b_dtRechitCluster2JetVetoE_0p8);
   fChain->SetBranchAddress("dtRechitCluster2MuonVetoPt_0p6", &dtRechitCluster2MuonVetoPt_0p6, &b_dtRechitCluster2MuonVetoPt_0p6);
   fChain->SetBranchAddress("dtRechitCluster2MuonVetoPt_0p8", &dtRechitCluster2MuonVetoPt_0p8, &b_dtRechitCluster2MuonVetoPt_0p8);
   fChain->SetBranchAddress("dtRechitCluster2MuonVetoE_0p6", &dtRechitCluster2MuonVetoE_0p6, &b_dtRechitCluster2MuonVetoE_0p6);
   fChain->SetBranchAddress("dtRechitCluster2MuonVetoE_0p8", &dtRechitCluster2MuonVetoE_0p8, &b_dtRechitCluster2MuonVetoE_0p8);
   fChain->SetBranchAddress("dtRechitCluster2MuonVetoLooseIso", &dtRechitCluster2MuonVetoLooseIso, &b_dtRechitCluster2MuonVetoLooseIso);
   fChain->SetBranchAddress("dtRechitCluster2MuonVetoTightIso", &dtRechitCluster2MuonVetoTightIso, &b_dtRechitCluster2MuonVetoTightIso);
   fChain->SetBranchAddress("dtRechitCluster2MuonVetoVTightIso", &dtRechitCluster2MuonVetoVTightIso, &b_dtRechitCluster2MuonVetoVTightIso);
   fChain->SetBranchAddress("dtRechitCluster2MuonVetoVVTightIso", &dtRechitCluster2MuonVetoVVTightIso, &b_dtRechitCluster2MuonVetoVVTightIso);
   fChain->SetBranchAddress("dtRechitCluster2MuonVetoTightId", &dtRechitCluster2MuonVetoTightId, &b_dtRechitCluster2MuonVetoTightId);
   fChain->SetBranchAddress("dtRechitCluster2MuonVetoLooseId", &dtRechitCluster2MuonVetoLooseId, &b_dtRechitCluster2MuonVetoLooseId);
   fChain->SetBranchAddress("dtRechitCluster2MuonVetoGlobal", &dtRechitCluster2MuonVetoGlobal, &b_dtRechitCluster2MuonVetoGlobal);
   fChain->SetBranchAddress("dtRechitCluster2MuonVetoIso", &dtRechitCluster2MuonVetoIso, &b_dtRechitCluster2MuonVetoIso);
   fChain->SetBranchAddress("dtRechitCluster2IsoMuonVetoPt", &dtRechitCluster2IsoMuonVetoPt, &b_dtRechitCluster2IsoMuonVetoPt);
   fChain->SetBranchAddress("dtRechitCluster2IsoMuonVetoE", &dtRechitCluster2IsoMuonVetoE, &b_dtRechitCluster2IsoMuonVetoE);
   fChain->SetBranchAddress("dtRechitCluster2IsoMuonVetoPhi", &dtRechitCluster2IsoMuonVetoPhi, &b_dtRechitCluster2IsoMuonVetoPhi);
   fChain->SetBranchAddress("dtRechitCluster2IsoMuonVetoEta", &dtRechitCluster2IsoMuonVetoEta, &b_dtRechitCluster2IsoMuonVetoEta);
   fChain->SetBranchAddress("dtRechitCluster2GenMuonVetoPt", &dtRechitCluster2GenMuonVetoPt, &b_dtRechitCluster2GenMuonVetoPt);
   fChain->SetBranchAddress("dtRechitCluster2MuonVetoType", &dtRechitCluster2MuonVetoType, &b_dtRechitCluster2MuonVetoType);
   fChain->SetBranchAddress("dtRechitCluster2ZLep1", &dtRechitCluster2ZLep1, &b_dtRechitCluster2ZLep1);
   fChain->SetBranchAddress("dtRechitCluster2ZLep2", &dtRechitCluster2ZLep2, &b_dtRechitCluster2ZLep2);
   fChain->SetBranchAddress("dtRechitCluster2ZLep1Tag", &dtRechitCluster2ZLep1Tag, &b_dtRechitCluster2ZLep1Tag);
   fChain->SetBranchAddress("dtRechitCluster2ZLep2Tag", &dtRechitCluster2ZLep2Tag, &b_dtRechitCluster2ZLep2Tag);
   fChain->SetBranchAddress("dtRechitCluster2ZLep1Id", &dtRechitCluster2ZLep1Id, &b_dtRechitCluster2ZLep1Id);
   fChain->SetBranchAddress("dtRechitCluster2ZLep2Id", &dtRechitCluster2ZLep2Id, &b_dtRechitCluster2ZLep2Id);
   fChain->SetBranchAddress("dtRechitCluster22ZLep1LooseIso", &dtRechitCluster22ZLep1LooseIso, &b_dtRechitCluster22ZLep1LooseIso);
   fChain->SetBranchAddress("dtRechitCluster2ZLep1TightIso", &dtRechitCluster2ZLep1TightIso, &b_dtRechitCluster2ZLep1TightIso);
   fChain->SetBranchAddress("dtRechitCluster2ZLep1VTightIso", &dtRechitCluster2ZLep1VTightIso, &b_dtRechitCluster2ZLep1VTightIso);
   fChain->SetBranchAddress("dtRechitCluster2ZLep1VVTightIso", &dtRechitCluster2ZLep1VVTightIso, &b_dtRechitCluster2ZLep1VVTightIso);
   fChain->SetBranchAddress("dtRechitCluster2ZLep1TightId", &dtRechitCluster2ZLep1TightId, &b_dtRechitCluster2ZLep1TightId);
   fChain->SetBranchAddress("dtRechitCluster2ZLep2LooseIso", &dtRechitCluster2ZLep2LooseIso, &b_dtRechitCluster2ZLep2LooseIso);
   fChain->SetBranchAddress("dtRechitCluster2ZLep2TightIso", &dtRechitCluster2ZLep2TightIso, &b_dtRechitCluster2ZLep2TightIso);
   fChain->SetBranchAddress("dtRechitCluster2ZLep2VTightIso", &dtRechitCluster2ZLep2VTightIso, &b_dtRechitCluster2ZLep2VTightIso);
   fChain->SetBranchAddress("dtRechitCluster2ZLep2VVTightIso", &dtRechitCluster2ZLep2VVTightIso, &b_dtRechitCluster2ZLep2VVTightIso);
   fChain->SetBranchAddress("dtRechitCluster2ZLep2TightId", &dtRechitCluster2ZLep2TightId, &b_dtRechitCluster2ZLep2TightId);
   fChain->SetBranchAddress("dtRechitCluster2Size", &dtRechitCluster2Size, &b_dtRechitCluster2Size);
   fChain->SetBranchAddress("dtRechitCluster2Me11Ratio", &dtRechitCluster2Me11Ratio, &b_dtRechitCluster2Me11Ratio);
   fChain->SetBranchAddress("dtRechitCluster2Me12Ratio", &dtRechitCluster2Me12Ratio, &b_dtRechitCluster2Me12Ratio);
   fChain->SetBranchAddress("dtRechitCluster2NStation", &dtRechitCluster2NStation, &b_dtRechitCluster2NStation);
   fChain->SetBranchAddress("dtRechitCluster2NStation5", &dtRechitCluster2NStation5, &b_dtRechitCluster2NStation5);
   fChain->SetBranchAddress("dtRechitCluster2NStation10", &dtRechitCluster2NStation10, &b_dtRechitCluster2NStation10);
   fChain->SetBranchAddress("dtRechitCluster2NStation10perc", &dtRechitCluster2NStation10perc, &b_dtRechitCluster2NStation10perc);
   fChain->SetBranchAddress("dtRechitCluster2AvgStation", &dtRechitCluster2AvgStation, &b_dtRechitCluster2AvgStation);
   fChain->SetBranchAddress("dtRechitCluster2AvgStation5", &dtRechitCluster2AvgStation5, &b_dtRechitCluster2AvgStation5);
   fChain->SetBranchAddress("dtRechitCluster2AvgStation10", &dtRechitCluster2AvgStation10, &b_dtRechitCluster2AvgStation10);
   fChain->SetBranchAddress("dtRechitCluster2AvgStation10perc", &dtRechitCluster2AvgStation10perc, &b_dtRechitCluster2AvgStation10perc);
   fChain->SetBranchAddress("dtRechitCluster2MaxStation", &dtRechitCluster2MaxStation, &b_dtRechitCluster2MaxStation);
   fChain->SetBranchAddress("dtRechitCluster2MaxStationRatio", &dtRechitCluster2MaxStationRatio, &b_dtRechitCluster2MaxStationRatio);
   fChain->SetBranchAddress("dtRechitCluster2NChamber", &dtRechitCluster2NChamber, &b_dtRechitCluster2NChamber);
   fChain->SetBranchAddress("dtRechitCluster2MaxChamber", &dtRechitCluster2MaxChamber, &b_dtRechitCluster2MaxChamber);
   fChain->SetBranchAddress("dtRechitCluster2MaxChamberRatio", &dtRechitCluster2MaxChamberRatio, &b_dtRechitCluster2MaxChamberRatio);
   fChain->SetBranchAddress("dtRechitCluster2NSegmentStation1", &dtRechitCluster2NSegmentStation1, &b_dtRechitCluster2NSegmentStation1);
   fChain->SetBranchAddress("dtRechitCluster2NSegmentStation2", &dtRechitCluster2NSegmentStation2, &b_dtRechitCluster2NSegmentStation2);
   fChain->SetBranchAddress("dtRechitCluster2NSegmentStation3", &dtRechitCluster2NSegmentStation3, &b_dtRechitCluster2NSegmentStation3);
   fChain->SetBranchAddress("dtRechitCluster2NSegmentStation4", &dtRechitCluster2NSegmentStation4, &b_dtRechitCluster2NSegmentStation4);
   fChain->SetBranchAddress("dtRechitCluster2NRechitChamberPlus11", &dtRechitCluster2NRechitChamberPlus11, &b_dtRechitCluster2NRechitChamberPlus11);
   fChain->SetBranchAddress("dtRechitCluster2NRechitChamberPlus12", &dtRechitCluster2NRechitChamberPlus12, &b_dtRechitCluster2NRechitChamberPlus12);
   fChain->SetBranchAddress("dtRechitCluster2NRechitChamberPlus13", &dtRechitCluster2NRechitChamberPlus13, &b_dtRechitCluster2NRechitChamberPlus13);
   fChain->SetBranchAddress("dtRechitCluster2NRechitChamberPlus21", &dtRechitCluster2NRechitChamberPlus21, &b_dtRechitCluster2NRechitChamberPlus21);
   fChain->SetBranchAddress("dtRechitCluster2NRechitChamberPlus22", &dtRechitCluster2NRechitChamberPlus22, &b_dtRechitCluster2NRechitChamberPlus22);
   fChain->SetBranchAddress("dtRechitCluster2NRechitChamberPlus31", &dtRechitCluster2NRechitChamberPlus31, &b_dtRechitCluster2NRechitChamberPlus31);
   fChain->SetBranchAddress("dtRechitCluster2NRechitChamberPlus32", &dtRechitCluster2NRechitChamberPlus32, &b_dtRechitCluster2NRechitChamberPlus32);
   fChain->SetBranchAddress("dtRechitCluster2NRechitChamberPlus41", &dtRechitCluster2NRechitChamberPlus41, &b_dtRechitCluster2NRechitChamberPlus41);
   fChain->SetBranchAddress("dtRechitCluster2NRechitChamberPlus42", &dtRechitCluster2NRechitChamberPlus42, &b_dtRechitCluster2NRechitChamberPlus42);
   fChain->SetBranchAddress("dtRechitCluster2NRechitChamberMinus11", &dtRechitCluster2NRechitChamberMinus11, &b_dtRechitCluster2NRechitChamberMinus11);
   fChain->SetBranchAddress("dtRechitCluster2NRechitChamberMinus12", &dtRechitCluster2NRechitChamberMinus12, &b_dtRechitCluster2NRechitChamberMinus12);
   fChain->SetBranchAddress("dtRechitCluster2NRechitChamberMinus13", &dtRechitCluster2NRechitChamberMinus13, &b_dtRechitCluster2NRechitChamberMinus13);
   fChain->SetBranchAddress("dtRechitCluster2NRechitChamberMinus21", &dtRechitCluster2NRechitChamberMinus21, &b_dtRechitCluster2NRechitChamberMinus21);
   fChain->SetBranchAddress("dtRechitCluster2NRechitChamberMinus22", &dtRechitCluster2NRechitChamberMinus22, &b_dtRechitCluster2NRechitChamberMinus22);
   fChain->SetBranchAddress("dtRechitCluster2NRechitChamberMinus31", &dtRechitCluster2NRechitChamberMinus31, &b_dtRechitCluster2NRechitChamberMinus31);
   fChain->SetBranchAddress("dtRechitCluster2NRechitChamberMinus32", &dtRechitCluster2NRechitChamberMinus32, &b_dtRechitCluster2NRechitChamberMinus32);
   fChain->SetBranchAddress("dtRechitCluster2NRechitChamberMinus41", &dtRechitCluster2NRechitChamberMinus41, &b_dtRechitCluster2NRechitChamberMinus41);
   fChain->SetBranchAddress("dtRechitCluster2NRechitChamberMinus42", &dtRechitCluster2NRechitChamberMinus42, &b_dtRechitCluster2NRechitChamberMinus42);
   fChain->SetBranchAddress("dtRechitCluster2Met_dPhi", &dtRechitCluster2Met_dPhi, &b_dtRechitCluster2Met_dPhi);
   fChain->SetBranchAddress("dtRechitCluster2MetXYCorr_dPhi", &dtRechitCluster2MetXYCorr_dPhi, &b_dtRechitCluster2MetXYCorr_dPhi);
   fChain->SetBranchAddress("dtRechitCluster2NLayersChamberPlus11", &dtRechitCluster2NLayersChamberPlus11, &b_dtRechitCluster2NLayersChamberPlus11);
   fChain->SetBranchAddress("dtRechitCluster2NLayersChamberPlus12", &dtRechitCluster2NLayersChamberPlus12, &b_dtRechitCluster2NLayersChamberPlus12);
   fChain->SetBranchAddress("dtRechitCluster2NLayersChamberPlus13", &dtRechitCluster2NLayersChamberPlus13, &b_dtRechitCluster2NLayersChamberPlus13);
   fChain->SetBranchAddress("dtRechitCluster2NLayersChamberPlus21", &dtRechitCluster2NLayersChamberPlus21, &b_dtRechitCluster2NLayersChamberPlus21);
   fChain->SetBranchAddress("dtRechitCluster2NLayersChamberPlus22", &dtRechitCluster2NLayersChamberPlus22, &b_dtRechitCluster2NLayersChamberPlus22);
   fChain->SetBranchAddress("dtRechitCluster2NLayersChamberPlus31", &dtRechitCluster2NLayersChamberPlus31, &b_dtRechitCluster2NLayersChamberPlus31);
   fChain->SetBranchAddress("dtRechitCluster2NLayersChamberPlus32", &dtRechitCluster2NLayersChamberPlus32, &b_dtRechitCluster2NLayersChamberPlus32);
   fChain->SetBranchAddress("dtRechitCluster2NLayersChamberPlus41", &dtRechitCluster2NLayersChamberPlus41, &b_dtRechitCluster2NLayersChamberPlus41);
   fChain->SetBranchAddress("dtRechitCluster2NLayersChamberPlus42", &dtRechitCluster2NLayersChamberPlus42, &b_dtRechitCluster2NLayersChamberPlus42);
   fChain->SetBranchAddress("dtRechitCluster2NLayersChamberMinus11", &dtRechitCluster2NLayersChamberMinus11, &b_dtRechitCluster2NLayersChamberMinus11);
   fChain->SetBranchAddress("dtRechitCluster2NLayersChamberMinus12", &dtRechitCluster2NLayersChamberMinus12, &b_dtRechitCluster2NLayersChamberMinus12);
   fChain->SetBranchAddress("dtRechitCluster2NLayersChamberMinus13", &dtRechitCluster2NLayersChamberMinus13, &b_dtRechitCluster2NLayersChamberMinus13);
   fChain->SetBranchAddress("dtRechitCluster2NLayersChamberMinus21", &dtRechitCluster2NLayersChamberMinus21, &b_dtRechitCluster2NLayersChamberMinus21);
   fChain->SetBranchAddress("dtRechitCluster2NLayersChamberMinus22", &dtRechitCluster2NLayersChamberMinus22, &b_dtRechitCluster2NLayersChamberMinus22);
   fChain->SetBranchAddress("dtRechitCluster2NLayersChamberMinus31", &dtRechitCluster2NLayersChamberMinus31, &b_dtRechitCluster2NLayersChamberMinus31);
   fChain->SetBranchAddress("dtRechitCluster2NLayersChamberMinus32", &dtRechitCluster2NLayersChamberMinus32, &b_dtRechitCluster2NLayersChamberMinus32);
   fChain->SetBranchAddress("dtRechitCluster2NLayersChamberMinus41", &dtRechitCluster2NLayersChamberMinus41, &b_dtRechitCluster2NLayersChamberMinus41);
   fChain->SetBranchAddress("dtRechitCluster2NLayersChamberMinus42", &dtRechitCluster2NLayersChamberMinus42, &b_dtRechitCluster2NLayersChamberMinus42);
   fChain->SetBranchAddress("dtRechitCluster2MetHEM_dPhi", &dtRechitCluster2MetHEM_dPhi, &b_dtRechitCluster2MetHEM_dPhi);
   fChain->SetBranchAddress("dtRechitCluster2MetHEMXYCorr_dPhi", &dtRechitCluster2MetHEMXYCorr_dPhi, &b_dtRechitCluster2MetHEMXYCorr_dPhi);
   fChain->SetBranchAddress("dtRechitCluster2MetEENoise_dPhi", &dtRechitCluster2MetEENoise_dPhi, &b_dtRechitCluster2MetEENoise_dPhi);
   fChain->SetBranchAddress("dtRechitCluster2MetEENoiseXYCorr_dPhi", &dtRechitCluster2MetEENoiseXYCorr_dPhi, &b_dtRechitCluster2MetEENoiseXYCorr_dPhi);
   fChain->SetBranchAddress("dtRechitCluster2MetJesUp_dPhi", &dtRechitCluster2MetJesUp_dPhi, &b_dtRechitCluster2MetJesUp_dPhi);
   fChain->SetBranchAddress("dtRechitCluster2MetJesDown_dPhi", &dtRechitCluster2MetJesDown_dPhi, &b_dtRechitCluster2MetJesDown_dPhi);
   fChain->SetBranchAddress("dtRechitCluster2_match_dtSeg_0p5", &dtRechitCluster2_match_dtSeg_0p5, &b_dtRechitCluster2_match_dtSeg_0p5);
   fChain->SetBranchAddress("dtRechitCluster2_match_dtSegTime_0p5", &dtRechitCluster2_match_dtSegTime_0p5, &b_dtRechitCluster2_match_dtSegTime_0p5);
   fChain->SetBranchAddress("dtRechitCluster2_match_dtSeg_0p4", &dtRechitCluster2_match_dtSeg_0p4, &b_dtRechitCluster2_match_dtSeg_0p4);
   fChain->SetBranchAddress("dtRechitCluster2_match_dtSegTime_0p4", &dtRechitCluster2_match_dtSegTime_0p4, &b_dtRechitCluster2_match_dtSegTime_0p4);
   fChain->SetBranchAddress("dtRechitCluster2_match_dtSegTimeSpread_0p5", &dtRechitCluster2_match_dtSegTimeSpread_0p5, &b_dtRechitCluster2_match_dtSegTimeSpread_0p5);
   fChain->SetBranchAddress("dtRechitCluster2_match_dtSegTimeSpread_0p4", &dtRechitCluster2_match_dtSegTimeSpread_0p4, &b_dtRechitCluster2_match_dtSegTimeSpread_0p4);
   fChain->SetBranchAddress("dtRechitCluster2_match_dtSeg_sameStation_0p5", &dtRechitCluster2_match_dtSeg_sameStation_0p5, &b_dtRechitCluster2_match_dtSeg_sameStation_0p5);
   fChain->SetBranchAddress("dtRechitCluster2_match_dtSegTime_sameStation_0p5", &dtRechitCluster2_match_dtSegTime_sameStation_0p5, &b_dtRechitCluster2_match_dtSegTime_sameStation_0p5);
   fChain->SetBranchAddress("dtRechitCluster2_match_dtSeg_sameStation_0p4", &dtRechitCluster2_match_dtSeg_sameStation_0p4, &b_dtRechitCluster2_match_dtSeg_sameStation_0p4);
   fChain->SetBranchAddress("dtRechitCluster2_match_dtSegTime_sameStation_0p4", &dtRechitCluster2_match_dtSegTime_sameStation_0p4, &b_dtRechitCluster2_match_dtSegTime_sameStation_0p4);
   fChain->SetBranchAddress("dtRechitCluster2_match_dtSegTimeSpread_sameStation_0p5", &dtRechitCluster2_match_dtSegTimeSpread_sameStation_0p5, &b_dtRechitCluster2_match_dtSegTimeSpread_sameStation_0p5);
   fChain->SetBranchAddress("dtRechitCluster2_match_dtSegTimeSpread_sameStation_0p4", &dtRechitCluster2_match_dtSegTimeSpread_sameStation_0p4, &b_dtRechitCluster2_match_dtSegTimeSpread_sameStation_0p4);
   fChain->SetBranchAddress("dtRechitCluster2_match_RPCTime_dPhi0p5", &dtRechitCluster2_match_RPCTime_dPhi0p5, &b_dtRechitCluster2_match_RPCTime_dPhi0p5);
   fChain->SetBranchAddress("dtRechitCluster2_match_RPCTimeSpread_dPhi0p5", &dtRechitCluster2_match_RPCTimeSpread_dPhi0p5, &b_dtRechitCluster2_match_RPCTimeSpread_dPhi0p5);
   fChain->SetBranchAddress("dtRechitCluster2_match_RPCTime_dR0p4", &dtRechitCluster2_match_RPCTime_dR0p4, &b_dtRechitCluster2_match_RPCTime_dR0p4);
   fChain->SetBranchAddress("dtRechitCluster2_match_RPCTimeSpread_dR0p4", &dtRechitCluster2_match_RPCTimeSpread_dR0p4, &b_dtRechitCluster2_match_RPCTimeSpread_dR0p4);
   fChain->SetBranchAddress("dtRechitCluster2_match_RPChits_dR0p4", &dtRechitCluster2_match_RPChits_dR0p4, &b_dtRechitCluster2_match_RPChits_dR0p4);
   fChain->SetBranchAddress("dtRechitCluster2_match_RPCTime_sameStation_dR0p4", &dtRechitCluster2_match_RPCTime_sameStation_dR0p4, &b_dtRechitCluster2_match_RPCTime_sameStation_dR0p4);
   fChain->SetBranchAddress("dtRechitCluster2_match_RPCTimeSpread_sameStation_dR0p4", &dtRechitCluster2_match_RPCTimeSpread_sameStation_dR0p4, &b_dtRechitCluster2_match_RPCTimeSpread_sameStation_dR0p4);
   fChain->SetBranchAddress("dtRechitCluster2_match_RPChits_sameStation_dR0p4", &dtRechitCluster2_match_RPChits_sameStation_dR0p4, &b_dtRechitCluster2_match_RPChits_sameStation_dR0p4);
   fChain->SetBranchAddress("dtRechitCluster2_match_gParticle_Id", &dtRechitCluster2_match_gParticle_Id, &b_dtRechitCluster2_match_gParticle_Id);
   fChain->SetBranchAddress("dtRechitCluster2_match_gParticle_Pt", &dtRechitCluster2_match_gParticle_Pt, &b_dtRechitCluster2_match_gParticle_Pt);
   fChain->SetBranchAddress("dtRechitCluster2_match_gParticle_Eta", &dtRechitCluster2_match_gParticle_Eta, &b_dtRechitCluster2_match_gParticle_Eta);
   fChain->SetBranchAddress("dtRechitCluster2_match_gParticle_Phi", &dtRechitCluster2_match_gParticle_Phi, &b_dtRechitCluster2_match_gParticle_Phi);
   fChain->SetBranchAddress("dtRechitCluster2_match_gParticle_E", &dtRechitCluster2_match_gParticle_E, &b_dtRechitCluster2_match_gParticle_E);
   fChain->SetBranchAddress("dtRechitCluster2_match_gParticle_Status", &dtRechitCluster2_match_gParticle_Status, &b_dtRechitCluster2_match_gParticle_Status);
   fChain->SetBranchAddress("dtRechitCluster2_match_gParticle_MotherId", &dtRechitCluster2_match_gParticle_MotherId, &b_dtRechitCluster2_match_gParticle_MotherId);
   fChain->SetBranchAddress("dtRechitCluster2_match_gParticle_deltaR", &dtRechitCluster2_match_gParticle_deltaR, &b_dtRechitCluster2_match_gParticle_deltaR);
   fChain->SetBranchAddress("dtRechitCluster2_match_RPChits_dPhi0p5", &dtRechitCluster2_match_RPChits_dPhi0p5, &b_dtRechitCluster2_match_RPChits_dPhi0p5);
   fChain->SetBranchAddress("dtRechitCluster2_match_RPCBx_dPhi0p5", &dtRechitCluster2_match_RPCBx_dPhi0p5, &b_dtRechitCluster2_match_RPCBx_dPhi0p5);
   fChain->SetBranchAddress("dtRechitCluster2_match_RB1_0p4", &dtRechitCluster2_match_RB1_0p4, &b_dtRechitCluster2_match_RB1_0p4);
   fChain->SetBranchAddress("dtRechitCluster2_match_RB1_dPhi0p5", &dtRechitCluster2_match_RB1_dPhi0p5, &b_dtRechitCluster2_match_RB1_dPhi0p5);
   fChain->SetBranchAddress("dtRechitCluster2_match_MB1Seg_0p4", &dtRechitCluster2_match_MB1Seg_0p4, &b_dtRechitCluster2_match_MB1Seg_0p4);
   fChain->SetBranchAddress("dtRechitCluster2_match_MB1Seg_0p5", &dtRechitCluster2_match_MB1Seg_0p5, &b_dtRechitCluster2_match_MB1Seg_0p5);
   fChain->SetBranchAddress("dtRechitCluster2_match_MB1hits_0p4", &dtRechitCluster2_match_MB1hits_0p4, &b_dtRechitCluster2_match_MB1hits_0p4);
   fChain->SetBranchAddress("dtRechitCluster2_match_MB1hits_0p5", &dtRechitCluster2_match_MB1hits_0p5, &b_dtRechitCluster2_match_MB1hits_0p5);
   fChain->SetBranchAddress("dtRechitCluster2_match_MB1hits_cosmics_plus", &dtRechitCluster2_match_MB1hits_cosmics_plus, &b_dtRechitCluster2_match_MB1hits_cosmics_plus);
   fChain->SetBranchAddress("dtRechitCluster2_match_MB1hits_cosmics_minus", &dtRechitCluster2_match_MB1hits_cosmics_minus, &b_dtRechitCluster2_match_MB1hits_cosmics_minus);
   fChain->SetBranchAddress("gLLP_multiplicity", gLLP_multiplicity, &b_gLLP_multiplicity);
   fChain->SetBranchAddress("gLLP_multiplicity20", gLLP_multiplicity20, &b_gLLP_multiplicity20);
   fChain->SetBranchAddress("gLLP_EM_multiplicity", gLLP_EM_multiplicity, &b_gLLP_EM_multiplicity);
   fChain->SetBranchAddress("gLLP_eta", gLLP_eta, &b_gLLP_eta);
   fChain->SetBranchAddress("gLLP_phi", gLLP_phi, &b_gLLP_phi);
   fChain->SetBranchAddress("gLLP_csc", gLLP_csc, &b_gLLP_csc);
   fChain->SetBranchAddress("gLLP_dt", gLLP_dt, &b_gLLP_dt);
   fChain->SetBranchAddress("gLLP_beta", gLLP_beta, &b_gLLP_beta);
   fChain->SetBranchAddress("gLLP_maxMatchedDis", gLLP_maxMatchedDis, &b_gLLP_maxMatchedDis);
   fChain->SetBranchAddress("gLLP_e", gLLP_e, &b_gLLP_e);
   fChain->SetBranchAddress("gLLP_pt", gLLP_pt, &b_gLLP_pt);
   fChain->SetBranchAddress("gLLP_lepdPhi", gLLP_lepdPhi, &b_gLLP_lepdPhi);
   fChain->SetBranchAddress("gLLP_daughterKaon", gLLP_daughterKaon, &b_gLLP_daughterKaon);
   fChain->SetBranchAddress("gLLP_ctau", gLLP_ctau, &b_gLLP_ctau);
   fChain->SetBranchAddress("gLLP_EMFracE", gLLP_EMFracE, &b_gLLP_EMFracE);
   fChain->SetBranchAddress("gLLP_EMFracEz", gLLP_EMFracEz, &b_gLLP_EMFracEz);
   fChain->SetBranchAddress("gLLP_EMFracP", gLLP_EMFracP, &b_gLLP_EMFracP);
   fChain->SetBranchAddress("gLLP_EMFracPz", gLLP_EMFracPz, &b_gLLP_EMFracPz);
   fChain->SetBranchAddress("gLLP_visE", gLLP_visE, &b_gLLP_visE);
   fChain->SetBranchAddress("gLLP_visE20", gLLP_visE20, &b_gLLP_visE20);
   fChain->SetBranchAddress("gLLP_visEz", gLLP_visEz, &b_gLLP_visEz);
   fChain->SetBranchAddress("gLLP_visP", gLLP_visP, &b_gLLP_visP);
   fChain->SetBranchAddress("gLLP_visPz", gLLP_visPz, &b_gLLP_visPz);
   fChain->SetBranchAddress("gLLP_match_jet_pt", gLLP_match_jet_pt, &b_gLLP_match_jet_pt);
   fChain->SetBranchAddress("gLLP_match_jet_index", gLLP_match_jet_index, &b_gLLP_match_jet_index);
   fChain->SetBranchAddress("gLLP_match_jet_minDeltaR", gLLP_match_jet_minDeltaR, &b_gLLP_match_jet_minDeltaR);
   fChain->SetBranchAddress("gLLP_decay_vertex_r", gLLP_decay_vertex_r, &b_gLLP_decay_vertex_r);
   fChain->SetBranchAddress("gLLP_decay_vertex_x", gLLP_decay_vertex_x, &b_gLLP_decay_vertex_x);
   fChain->SetBranchAddress("gLLP_decay_vertex_y", gLLP_decay_vertex_y, &b_gLLP_decay_vertex_y);
   fChain->SetBranchAddress("gLLP_decay_vertex_z", gLLP_decay_vertex_z, &b_gLLP_decay_vertex_z);
   fChain->SetBranchAddress("gLLP_deltaR", &gLLP_deltaR, &b_gLLP_deltaR);
   fChain->SetBranchAddress("gLLP_daughter_deltaR", gLLP_daughter_deltaR, &b_gLLP_daughter_deltaR);
   fChain->SetBranchAddress("gLLP_daughter_pt", gLLP_daughter_pt, &b_gLLP_daughter_pt);
   fChain->SetBranchAddress("gLLP_daughter_id", gLLP_daughter_id, &b_gLLP_daughter_id);
   fChain->SetBranchAddress("gLLP_daughter_eta", gLLP_daughter_eta, &b_gLLP_daughter_eta);
   fChain->SetBranchAddress("gLLP_daughter_phi", gLLP_daughter_phi, &b_gLLP_daughter_phi);
   fChain->SetBranchAddress("gLLP_daughter_e", gLLP_daughter_e, &b_gLLP_daughter_e);
   fChain->SetBranchAddress("gLLP_daughter_mass", gLLP_daughter_mass, &b_gLLP_daughter_mass);
   fChain->SetBranchAddress("nGlobalMuons", &nGlobalMuons, &b_nGlobalMuons);
   fChain->SetBranchAddress("GlobalMuonPt", &GlobalMuonPt, &b_GlobalMuonPt);
   fChain->SetBranchAddress("GlobalMuonEta", &GlobalMuonEta, &b_GlobalMuonEta);
   fChain->SetBranchAddress("GlobalMuonPhi", &GlobalMuonPhi, &b_GlobalMuonPhi);
   fChain->SetBranchAddress("GlobalMuonLooseId", &GlobalMuonLooseId, &b_GlobalMuonLooseId);
   fChain->SetBranchAddress("nLeptons", &nLeptons, &b_nLeptons);
   fChain->SetBranchAddress("lepE", lepE, &b_lepE);
   fChain->SetBranchAddress("lepPt", lepPt, &b_lepPt);
   fChain->SetBranchAddress("lepEta", lepEta, &b_lepEta);
   fChain->SetBranchAddress("lepPhi", lepPhi, &b_lepPhi);
   fChain->SetBranchAddress("lepPdgId", lepPdgId, &b_lepPdgId);
   fChain->SetBranchAddress("lepDZ", lepDZ, &b_lepDZ);
   fChain->SetBranchAddress("lepPassId", lepPassId, &b_lepPassId);
   fChain->SetBranchAddress("lepFromZ", lepFromZ, &b_lepFromZ);
   fChain->SetBranchAddress("lepEff", lepEff, &b_lepEff);
   fChain->SetBranchAddress("lepSF", lepSF, &b_lepSF);
   fChain->SetBranchAddress("lepTriggerSF", lepTriggerSF, &b_lepTriggerSF);
   fChain->SetBranchAddress("lepTightIdSF", lepTightIdSF, &b_lepTightIdSF);
   fChain->SetBranchAddress("lepLooseIdSF", lepLooseIdSF, &b_lepLooseIdSF);
   fChain->SetBranchAddress("lepTightIsoSF", lepTightIsoSF, &b_lepTightIsoSF);
   fChain->SetBranchAddress("lepLooseIsoSF", lepLooseIsoSF, &b_lepLooseIsoSF);
   fChain->SetBranchAddress("lepTag", lepTag, &b_lepTag);
   fChain->SetBranchAddress("lepPassLooseIso", lepPassLooseIso, &b_lepPassLooseIso);
   fChain->SetBranchAddress("lepPassTightIso", lepPassTightIso, &b_lepPassTightIso);
   fChain->SetBranchAddress("lepPassVTightIso", lepPassVTightIso, &b_lepPassVTightIso);
   fChain->SetBranchAddress("lepPassVVTightIso", lepPassVVTightIso, &b_lepPassVVTightIso);
   fChain->SetBranchAddress("MT", &MT, &b_MT);
   fChain->SetBranchAddress("ZMass1", &ZMass1, &b_ZMass1);
   fChain->SetBranchAddress("ZMass", &ZMass, &b_ZMass);
   fChain->SetBranchAddress("ZPt", &ZPt, &b_ZPt);
   fChain->SetBranchAddress("ZEta", &ZEta, &b_ZEta);
   fChain->SetBranchAddress("ZPhi", &ZPhi, &b_ZPhi);
   fChain->SetBranchAddress("ZleptonIndex1", &ZleptonIndex1, &b_ZleptonIndex1);
   fChain->SetBranchAddress("ZleptonIndex2", &ZleptonIndex2, &b_ZleptonIndex2);
   fChain->SetBranchAddress("nJets", &nJets, &b_nJets);
   fChain->SetBranchAddress("jetE", jetE, &b_jetE);
   fChain->SetBranchAddress("jetPt", jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetEta", jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetTime", jetTime, &b_jetTime);
   fChain->SetBranchAddress("jetPassId", jetPassId, &b_jetPassId);
   fChain->SetBranchAddress("jetPtJESUp", jetPtJESUp, &b_jetPtJESUp);
   fChain->SetBranchAddress("jetPtJESDown", jetPtJESDown, &b_jetPtJESDown);
   fChain->SetBranchAddress("jetEJESUp", jetEJESUp, &b_jetEJESUp);
   fChain->SetBranchAddress("jetEJESDown", jetEJESDown, &b_jetEJESDown);
   fChain->SetBranchAddress("JecUnc", JecUnc, &b_JecUnc);
   fChain->SetBranchAddress("jet_match_llp_pt", jet_match_llp_pt, &b_jet_match_llp_pt);
   fChain->SetBranchAddress("jet_match_llp_index", jet_match_llp_index, &b_jet_match_llp_index);
   fChain->SetBranchAddress("jet_match_llp_minDeltaR", jet_match_llp_minDeltaR, &b_jet_match_llp_minDeltaR);
   fChain->SetBranchAddress("jet_match_genJet_pt", jet_match_genJet_pt, &b_jet_match_genJet_pt);
   fChain->SetBranchAddress("jet_match_genJet_index", jet_match_genJet_index, &b_jet_match_genJet_index);
   fChain->SetBranchAddress("jet_match_genJet_minDeltaR", jet_match_genJet_minDeltaR, &b_jet_match_genJet_minDeltaR);
   fChain->SetBranchAddress("jetTightPassId", jetTightPassId, &b_jetTightPassId);
   fChain->SetBranchAddress("HLTDecision", HLTDecision, &b_HLTDecision);
   fChain->SetBranchAddress("METTrigger", &METTrigger, &b_METTrigger);
   fChain->SetBranchAddress("METNoMuTrigger", &METNoMuTrigger, &b_METNoMuTrigger);
   fChain->SetBranchAddress("jetChargedEMEnergyFraction", jetChargedEMEnergyFraction, &b_jetChargedEMEnergyFraction);
   fChain->SetBranchAddress("jetChargedHadronEnergyFraction", jetChargedHadronEnergyFraction, &b_jetChargedHadronEnergyFraction);
   fChain->SetBranchAddress("jetNeutralEMEnergyFraction", jetNeutralEMEnergyFraction, &b_jetNeutralEMEnergyFraction);
   fChain->SetBranchAddress("jetNeutralHadronEnergyFraction", jetNeutralHadronEnergyFraction, &b_jetNeutralHadronEnergyFraction);
   Notify();
}

Bool_t LLPEvent::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void LLPEvent::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t LLPEvent::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef LLPEvent_cxx
