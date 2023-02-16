"""main.py

Main script for the CMS Run3 analysis
"""

import os
import sys
import pathlib

import awkward as awk
import numba as nb
import awkward.numba
import uproot

import numpy as np
from tqdm.auto import tqdm
from functools import partial
from argparse import ArgumentParser

import ROOT as rt
from ROOT import TChain, TCanvas, TH1F, TH2F, TLatex, TGraph, RDataFrame
# from src.candidates import LLPCandidate, DTCandidate, CSCCandidate

from multiprocessing import Pool

####################################################################################################

branches = [
    'runNum',
    'MC_condition',
    'lumiSec',
    'evtNum',
    # 'mH',
    # 'mX',
    'ctau',
    # 'npv',
    # 'npu',
    'weight',
    'pileupWeight',
    # 'Flag_HBHENoiseFilter',
    # 'Flag_BadPFMuonFilter',
    # 'Flag_HBHEIsoNoiseFilter',
    # 'Flag_CSCTightHaloFilter',
    # 'Flag_globalSuperTightHalo2016Filter',
    # 'Flag_goodVertices',
    # 'Flag_ecalBadCalibFilter',
    # 'Flag_BadChargedCandidateFilter',
    # 'Flag_eeBadScFilter',
    # 'Flag_all',
    # 'Flag2_HBHENoiseFilter',
    # 'Flag2_HBHEIsoNoiseFilter',
    # 'Flag2_BadPFMuonFilter',
    # 'Flag2_globalSuperTightHalo2016Filter',
    # 'Flag2_globalTightHalo2016Filter',
    # 'Flag2_BadChargedCandidateFilter',
    # 'Flag2_EcalDeadCellTriggerPrimitiveFilter',
    # 'Flag2_ecalBadCalibFilter',
    # 'Flag2_eeBadScFilter',
    # 'Flag2_all',
    # 'rho',
    'met',
    'metPhi',
    'gHiggsPt',
    'gHiggsE',
    'gHiggsEta',
    'gHiggsPhi',
    'nCscRings',
    'nDtRings',
    'nCscRechitClusters',
    'cscRechitCluster_match_gLLP',
    'cscRechitCluster_match_gLLP_minDeltaR',
    'cscRechitCluster_match_gLLP_index',
    'cscRechitCluster_match_gLLP_eta',
    'cscRechitCluster_match_gLLP_phi',
    'cscRechitCluster_match_gLLP_decay_r',
    'cscRechitCluster_match_gLLP_decay_z',
    'cscRechitCluster_match_gLLP_csc',
    'cscRechitCluster_match_gLLP_dt',
    'cscRechitCluster_match_gLLP_e',
    'cscRechitClusterX',
    'cscRechitClusterY',
    'cscRechitClusterZ',
    'cscRechitClusterTimeWeighted',
    'cscRechitClusterTimeSpreadWeightedAll',
    'cscRechitClusterPhi',
    'cscRechitClusterEta',
    'cscRechitClusterJetVetoPt',
    'cscRechitClusterJetVetoLooseId',
    'cscRechitClusterJetVetoTightId',
    'cscRechitClusterJetVetoE',
    'cscRechitClusterMuonVetoPt',
    'cscRechitClusterMuonVetoE',
    'cscRechitClusterMuonVetoLooseId',
    'cscRechitClusterMuonVetoGlobal',
    'cscRechitCluster_match_dtSeg_0p4',
    'cscRechitCluster_match_MB1Seg_0p4',
    'cscRechitCluster_match_RE12_0p4',
    'cscRechitCluster_match_RB1_0p4',
    'cscRechitClusterSize',
    'cscRechitClusterNStation10',
    'cscRechitClusterAvgStation10',
    'cscRechitClusterMaxStation',
    'cscRechitClusterMaxStationRatio',
    'cscRechitClusterNChamber',
    'cscRechitClusterMaxChamber',
    'cscRechitClusterMaxChamberRatio',
    'cscRechitClusterNRechitChamberPlus11',
    'cscRechitClusterNRechitChamberPlus12',
    'cscRechitClusterNRechitChamberPlus13',
    'cscRechitClusterNRechitChamberPlus21',
    'cscRechitClusterNRechitChamberPlus22',
    'cscRechitClusterNRechitChamberPlus31',
    'cscRechitClusterNRechitChamberPlus32',
    'cscRechitClusterNRechitChamberPlus41',
    'cscRechitClusterNRechitChamberPlus42',
    'cscRechitClusterNRechitChamberMinus11',
    'cscRechitClusterNRechitChamberMinus12',
    'cscRechitClusterNRechitChamberMinus13',
    'cscRechitClusterNRechitChamberMinus21',
    'cscRechitClusterNRechitChamberMinus22',
    'cscRechitClusterNRechitChamberMinus31',
    'cscRechitClusterNRechitChamberMinus32',
    'cscRechitClusterNRechitChamberMinus41',
    'cscRechitClusterNRechitChamberMinus42',
    'cscRechitClusterMet_dPhi',
    'nDtRechitClusters',
    'dtRechitClusterNSegStation1',
    'dtRechitClusterNSegStation2',
    'dtRechitClusterNSegStation3',
    'dtRechitClusterNSegStation4',
    'dtRechitClusterNOppositeSegStation1',
    'dtRechitClusterNOppositeSegStation2',
    'dtRechitClusterNOppositeSegStation3',
    'dtRechitClusterNOppositeSegStation4',
    'dtRechitCluster_match_gLLP',
    'dtRechitCluster_match_gLLP_minDeltaR',
    'dtRechitCluster_match_gLLP_index',
    'dtRechitCluster_match_gLLP_eta',
    'dtRechitCluster_match_gLLP_phi',
    'dtRechitCluster_match_gLLP_decay_r',
    'dtRechitCluster_match_gLLP_decay_z',
    'dtRechitCluster_match_gLLP_csc',
    'dtRechitCluster_match_gLLP_dt',
    'dtRechitCluster_match_gLLP_e',
    'dtRechitClusterX',
    'dtRechitClusterY',
    'dtRechitClusterZ',
    'dtRechitClusterWheel',
    'dtRechitClusterPhi',
    'dtRechitClusterEta',
    'dtRechitClusterJetVetoPt',
    'dtRechitClusterJetVetoLooseId',
    'dtRechitClusterJetVetoTightId',
    'dtRechitClusterJetVetoE',
    'dtRechitClusterMuonVetoPt',
    'dtRechitClusterMuonVetoE',
    'dtRechitClusterMuonVetoTightId',
    'dtRechitClusterMuonVetoLooseId',
    'dtRechitClusterMuonVetoGlobal',
    'dtRechitClusterOverlap',
    'dtRechitClusterSize',
    'dtRechitClusterNoiseHit',
    'dtRechitClusterNoiseHitStation1',
    'dtRechitClusterNoiseHitStation2',
    'dtRechitClusterNoiseHitStation3',
    'dtRechitClusterNoiseHitStation4',
    'dtRechitClusterNStation10',
    'dtRechitClusterAvgStation10',
    'dtRechitClusterMaxStation',
    'dtRechitClusterMaxStationRatio',
    'dtRechitClusterNChamber',
    'dtRechitClusterMaxChamber',
    'dtRechitClusterMaxChamberRatio',
    'dtRechitClusterNHitStation1',
    'dtRechitClusterNHitStation2',
    'dtRechitClusterNHitStation3',
    'dtRechitClusterNHitStation4',
    'dtRechitClusterMet_dPhi',
    'dtRechitCluster_match_RPChits_dPhi0p5',
    'dtRechitCluster_match_RPCBx_dPhi0p5',
    'dtRechitCluster_match_RB1_0p4',
    'dtRechitCluster_match_RB1_dPhi0p5',
    'dtRechitCluster_match_MB1Seg_0p4',
    'dtRechitCluster_match_MB1Seg_0p5',
    'dtRechitCluster_match_MB1hits_0p4',
    'dtRechitCluster_match_MB1hits_0p5',
    'dtRechitCluster_match_MB1hits_cosmics_plus',
    'dtRechitCluster_match_MB1hits_cosmics_minus',
    'nGLLP',
    'gLLP_eta',
    'gLLP_phi',
    'gLLP_csc',
    'gLLP_dt',
    'gLLP_beta',
    'gLLP_e',
    'gLLP_pt',
    'gLLP_ctau',
    'gLLP_decay_vertex_r',
    'gLLP_decay_vertex_x',
    'gLLP_decay_vertex_y',
    'gLLP_decay_vertex_z',
    'nLeptons',
    'lepE',
    'lepPt',
    'lepEta',
    'lepPhi',
    'lepPdgId',
    'lepDZ',
    'lepTightId',
    'lepPassLooseIso',
    'lepPassTightIso',
    'lepPassVTightIso',
    'lepPassVVTightIso',
    'nJets',
    'jetE',
    'jetPt',
    'jetEta',
    'jetPhi',
    'jetTightPassId',
    'HLTDecision',
]


class Lepton():

  def __init__(self, event=None, index=0) -> None:
    if event is not None:
      self.process_event(event, index=index)

  def process_event(self, event, index) -> None:
    self.e = event['lepE'][index]
    self.pt = event['lepPt'][index]
    self.eta = event['lepEta'][index]
    self.phi = event['lepPhi'][index]

    self.PdgId = event['lepPdgId'][index]
    # self.DZ = event['lepDZ'][index]
    # self.TightId = event['lepTightId'][index]
    # self.PassLooseIso = event['lepPassLooseIso'][index]
    # self.PassTightIso = event['lepPassTightIso'][index]
    # self.PassVTightIso = event['lepPassVTightIso'][index]
    # self.PassVVTightIso = event['lepPassVVTightIso'][index]


class Jet():

  def __init__(self, event=None, index=0) -> None:
    if event is not None:
      self.process_event(event, index=index)

  def process_event(self, event, index) -> None:
    self.e = event['jetE'][index]
    self.pt = event['jetPt'][index]
    self.eta = event['jetEta'][index]
    self.phi = event['jetPhi'][index]

    # self.TightPassId = event['jetTightPassId'][index]


class gLLP():

  def __init__(self, event=None, index=0) -> None:
    if event is not None:
      self.process_event(event, index=index)

  def process_event(self, event, index) -> None:
    self.e = event['gLLP_e'][index]
    self.pt = event['gLLP_pt'][index]
    self.eta = event['gLLP_eta'][index]
    self.phi = event['gLLP_phi'][index]

    self.csc = event['gLLP_csc'][index]
    self.dt = event['gLLP_dt'][index]

    self.beta = event['gLLP_beta'][index]
    self.ctau = event['gLLP_ctau'][index]

    self.decay_vertex_r = event['gLLP_decay_vertex_r'][index]
    self.decay_vertex_x = event['gLLP_decay_vertex_x'][index]
    self.decay_vertex_y = event['gLLP_decay_vertex_y'][index]
    self.decay_vertex_z = event['gLLP_decay_vertex_z'][index]


class DTCluster():

  def __init__(self, event=None, index=0) -> None:
    if event is not None:
      self.process_event(event, index=index)

  def process_event(self, event, index) -> None:
    self.size = event['dtRechitClusterSize'][index]
    self.phi = event['dtRechitClusterPhi'][index]
    self.eta = event['dtRechitClusterEta'][index]

    self.met_dPhi = event['dtRechitClusterMet_dPhi'][index]

    self.x = event['dtRechitClusterX'][index]
    self.y = event['dtRechitClusterY'][index]
    self.z = event['dtRechitClusterZ'][index]

    self.NStation10 = event['dtRechitClusterNStation10'][index]
    self.AvgStation10 = event['dtRechitClusterAvgStation10'][index]
    self.MaxStation = event['dtRechitClusterMaxStation'][index]
    # self.MaxStationRatio = event['dtRechitClusterMaxStationRatio'][index]
    # self.NChamber = event['dtRechitClusterNChamber'][index]
    # self.MaxChamber = event['dtRechitClusterMaxChamber'][index]
    # self.MaxChamberRatio = event['dtRechitClusterMaxChamberRatio'][index]

    self.match_gLLP = event['dtRechitCluster_match_gLLP'][index]
    self.match_gLLP_e = event['dtRechitCluster_match_gLLP_e'][index]
    self.match_gLLP_eta = event['dtRechitCluster_match_gLLP_eta'][index]
    self.match_gLLP_phi = event['dtRechitCluster_match_gLLP_phi'][index]
    self.match_gLLP_decay_r = event['dtRechitCluster_match_gLLP_decay_r'][index]
    self.match_gLLP_decay_z = event['dtRechitCluster_match_gLLP_decay_z'][index]
    self.match_gLLP_minDeltaR = event['dtRechitCluster_match_gLLP_minDeltaR'][index]
    # self.match_gLLP_index = event['dtRechitCluster_match_gLLP_index'][index]
    self.match_gLLP_csc = event['dtRechitCluster_match_gLLP_csc'][index]
    self.match_gLLP_dt = event['dtRechitCluster_match_gLLP_dt'][index]

    # self.match_RPChits_dPhi0p5 = event['dtRechitCluster_match_RPChits_dPhi0p5'][index]
    # self.match_RPCBx_dPhi0p5 = event['dtRechitCluster_match_RPCBx_dPhi0p5'][index]
    # self.match_RB1_0p4 = event['dtRechitCluster_match_RB1_0p4'][index]
    # self.match_RB1_dPhi0p5 = event['dtRechitCluster_match_RB1_dPhi0p5'][index]
    # self.match_MB1Seg_0p4 = event['dtRechitCluster_match_MB1Seg_0p4'][index]
    # self.match_MB1Seg_0p5 = event['dtRechitCluster_match_MB1Seg_0p5'][index]
    # self.match_MB1hits_0p4 = event['dtRechitCluster_match_MB1hits_0p4'][index]
    # self.match_MB1hits_0p5 = event['dtRechitCluster_match_MB1hits_0p5'][index]
    # self.match_MB1hits_cosmics_plus = event['dtRechitCluster_match_MB1hits_cosmics_plus'][index]
    # self.match_MB1hits_cosmics_minus = event['dtRechitCluster_match_MB1hits_cosmics_minus'][index]

    self.NSegStation1 = event['dtRechitClusterNSegStation1'][index]
    self.NSegStation2 = event['dtRechitClusterNSegStation2'][index]
    self.NSegStation3 = event['dtRechitClusterNSegStation3'][index]
    self.NSegStation4 = event['dtRechitClusterNSegStation4'][index]
    self.NOppositeSegStation1 = event['dtRechitClusterNOppositeSegStation1'][index]
    self.NOppositeSegStation2 = event['dtRechitClusterNOppositeSegStation2'][index]
    self.NOppositeSegStation3 = event['dtRechitClusterNOppositeSegStation3'][index]
    self.NOppositeSegStation4 = event['dtRechitClusterNOppositeSegStation4'][index]
    # self.Wheel = event['dtRechitClusterWheel'][index]

    # self.JetVetoPt = event['dtRechitClusterJetVetoPt'][index]
    # self.JetVetoLooseId = event['dtRechitClusterJetVetoLooseId'][index]
    # self.JetVetoTightId = event['dtRechitClusterJetVetoTightId'][index]
    # self.JetVetoE = event['dtRechitClusterJetVetoE'][index]
    # self.MuonVetoPt = event['dtRechitClusterMuonVetoPt'][index]
    # self.MuonVetoE = event['dtRechitClusterMuonVetoE'][index]
    # self.MuonVetoTightId = event['dtRechitClusterMuonVetoTightId'][index]
    # self.MuonVetoLooseId = event['dtRechitClusterMuonVetoLooseId'][index]
    # self.MuonVetoGlobal = event['dtRechitClusterMuonVetoGlobal'][index]

    # self.Overlap = event['dtRechitClusterOverlap'][index]
    # self.NoiseHit = event['dtRechitClusterNoiseHit'][index]
    # self.NoiseHitStation1 = event['dtRechitClusterNoiseHitStation1'][index]
    # self.NoiseHitStation2 = event['dtRechitClusterNoiseHitStation2'][index]
    # self.NoiseHitStation3 = event['dtRechitClusterNoiseHitStation3'][index]
    # self.NoiseHitStation4 = event['dtRechitClusterNoiseHitStation4'][index]
    # self.NHitStation1 = event['dtRechitClusterNHitStation1'][index]
    # self.NHitStation2 = event['dtRechitClusterNHitStation2'][index]
    # self.NHitStation3 = event['dtRechitClusterNHitStation3'][index]
    # self.NHitStation4 = event['dtRechitClusterNHitStation4'][index]


class CSCCluster():

  def __init__(self, event=None, index=0) -> None:
    if event is not None:
      self.process_event(event, index=index)

  def process_event(self, event, index) -> None:
    self.size = event['cscRechitClusterSize'][index]
    self.phi = event['cscRechitClusterPhi'][index]
    self.eta = event['cscRechitClusterEta'][index]

    self.met_dPhi = event['cscRechitClusterMet_dPhi'][index]

    self.x = event['cscRechitClusterX'][index]
    self.y = event['cscRechitClusterY'][index]
    self.z = event['cscRechitClusterZ'][index]

    self.TimeWeighted = event['cscRechitClusterTimeWeighted'][index]
    self.TimeSpreadWeightedAll = event['cscRechitClusterTimeSpreadWeightedAll'][index]

    self.NStation10 = event['cscRechitClusterNStation10'][index]
    self.AvgStation10 = event['cscRechitClusterAvgStation10'][index]
    self.MaxStation = event['cscRechitClusterMaxStation'][index]
    # self.MaxStationRatio = event['cscRechitClusterMaxStationRatio'][index]
    # self.NChamber = event['cscRechitClusterNChamber'][index]
    # self.MaxChamber = event['cscRechitClusterMaxChamber'][index]
    # self.MaxChamberRatio = event['cscRechitClusterMaxChamberRatio'][index]

    self.match_gLLP = event['cscRechitCluster_match_gLLP'][index]
    self.match_gLLP_minDeltaR = event['cscRechitCluster_match_gLLP_minDeltaR'][index]
    # self.match_gLLP_index = event['cscRechitCluster_match_gLLP_index'][index]
    self.match_gLLP_eta = event['cscRechitCluster_match_gLLP_eta'][index]
    self.match_gLLP_phi = event['cscRechitCluster_match_gLLP_phi'][index]
    self.match_gLLP_decay_r = event['cscRechitCluster_match_gLLP_decay_r'][index]
    self.match_gLLP_decay_z = event['cscRechitCluster_match_gLLP_decay_z'][index]
    self.match_gLLP_csc = event['cscRechitCluster_match_gLLP_csc'][index]
    self.match_gLLP_dt = event['cscRechitCluster_match_gLLP_dt'][index]
    self.match_gLLP_e = event['cscRechitCluster_match_gLLP_e'][index]
    # self.match_dtSeg_0p4 = event['cscRechitCluster_match_dtSeg_0p4'][index]
    # self.match_MB1Seg_0p4 = event['cscRechitCluster_match_MB1Seg_0p4'][index]
    # self.match_RE12_0p4 = event['cscRechitCluster_match_RE12_0p4'][index]
    # self.match_RB1_0p4 = event['cscRechitCluster_match_RB1_0p4'][index]

    # self.JetVetoPt = event['cscRechitClusterJetVetoPt'][index]
    # self.JetVetoLooseId = event['cscRechitClusterJetVetoLooseId'][index]
    # self.JetVetoTightId = event['cscRechitClusterJetVetoTightId'][index]
    # self.JetVetoE = event['cscRechitClusterJetVetoE'][index]
    # self.MuonVetoPt = event['cscRechitClusterMuonVetoPt'][index]
    # self.MuonVetoE = event['cscRechitClusterMuonVetoE'][index]
    # self.MuonVetoLooseId = event['cscRechitClusterMuonVetoLooseId'][index]
    # self.MuonVetoGlobal = event['cscRechitClusterMuonVetoGlobal'][index]

    self.NRechitChamberPlus11 = event['cscRechitClusterNRechitChamberPlus11'][index]
    self.NRechitChamberPlus12 = event['cscRechitClusterNRechitChamberPlus12'][index]
    # self.NRechitChamberPlus13 = event['cscRechitClusterNRechitChamberPlus13'][index]
    # self.NRechitChamberPlus21 = event['cscRechitClusterNRechitChamberPlus21'][index]
    # self.NRechitChamberPlus22 = event['cscRechitClusterNRechitChamberPlus22'][index]
    # self.NRechitChamberPlus31 = event['cscRechitClusterNRechitChamberPlus31'][index]
    # self.NRechitChamberPlus32 = event['cscRechitClusterNRechitChamberPlus32'][index]
    # self.NRechitChamberPlus41 = event['cscRechitClusterNRechitChamberPlus41'][index]
    # self.NRechitChamberPlus42 = event['cscRechitClusterNRechitChamberPlus42'][index]
    self.NRechitChamberMinus11 = event['cscRechitClusterNRechitChamberMinus11'][index]
    self.NRechitChamberMinus12 = event['cscRechitClusterNRechitChamberMinus12'][index]
    # self.NRechitChamberMinus13 = event['cscRechitClusterNRechitChamberMinus13'][index]
    # self.NRechitChamberMinus21 = event['cscRechitClusterNRechitChamberMinus21'][index]
    # self.NRechitChamberMinus22 = event['cscRechitClusterNRechitChamberMinus22'][index]
    # self.NRechitChamberMinus31 = event['cscRechitClusterNRechitChamberMinus31'][index]
    # self.NRechitChamberMinus32 = event['cscRechitClusterNRechitChamberMinus32'][index]
    # self.NRechitChamberMinus41 = event['cscRechitClusterNRechitChamberMinus41'][index]
    # self.NRechitChamberMinus42 = event['cscRechitClusterNRechitChamberMinus42'][index]


class LLPEventCandidate():

  def __init__(self, event=None) -> None:
    if event is not None:
      self.process_event(event)

  def process_event(self, event):

    self.runNum = event['runNum']
    self.MC_condition = event['MC_condition']
    self.lumiSec = event['lumiSec']
    self.evtNum = event['evtNum']

    self.HLTDecision = event['HLTDecision']

    self.weight = event['weight']
    self.pileupWeight = event['pileupWeight']

    self.met = event['met']
    self.metPhi = event['metPhi']

    # self.mH = event['mH']
    # self.mX = event['mX']
    self.ctau = event['ctau']

    # self.npv = event['npv']
    # self.npu = event['npu']

    # self.Flag_HBHENoiseFilter = event['Flag_HBHENoiseFilter']
    # self.Flag_BadPFMuonFilter = event['Flag_BadPFMuonFilter']
    # self.Flag_HBHEIsoNoiseFilter = event['Flag_HBHEIsoNoiseFilter']
    # self.Flag_CSCTightHaloFilter = event['Flag_CSCTightHaloFilter']
    # self.Flag_globalSuperTightHalo2016Filter = event['Flag_globalSuperTightHalo2016Filter']
    # self.Flag_goodVertices = event['Flag_goodVertices']
    # self.Flag_ecalBadCalibFilter = event['Flag_ecalBadCalibFilter']
    # self.Flag_BadChargedCandidateFilter = event['Flag_BadChargedCandidateFilter']
    # self.Flag_eeBadScFilter = event['Flag_eeBadScFilter']
    # self.Flag_all = event['Flag_all']
    # self.Flag2_HBHENoiseFilter = event['Flag2_HBHENoiseFilter']
    # self.Flag2_HBHEIsoNoiseFilter = event['Flag2_HBHEIsoNoiseFilter']
    # self.Flag2_BadPFMuonFilter = event['Flag2_BadPFMuonFilter']
    # self.Flag2_globalSuperTightHalo2016Filter = event['Flag2_globalSuperTightHalo2016Filter']
    # self.Flag2_globalTightHalo2016Filter = event['Flag2_globalTightHalo2016Filter']
    # self.Flag2_BadChargedCandidateFilter = event['Flag2_BadChargedCandidateFilter']
    # self.Flag2_EcalDeadCellTriggerPrimitiveFilter = event['Flag2_EcalDeadCellTriggerPrimitiveFilter']
    # self.Flag2_ecalBadCalibFilter = event['Flag2_ecalBadCalibFilter']
    # self.Flag2_eeBadScFilter = event['Flag2_eeBadScFilter']
    # self.Flag2_all = event['Flag2_all']

    # self.rho = event['rho']

    self.gHiggsPt = event['gHiggsPt']
    self.gHiggsE = event['gHiggsE']
    self.gHiggsEta = event['gHiggsEta']
    self.gHiggsPhi = event['gHiggsPhi']

    self.nCscRings = event['nCscRings']
    self.nDtRings = event['nDtRings']

    self.csc_clusters = [CSCCluster(event, idx) for idx in range(event['nCscRechitClusters'])]
    self.dt_clusters = [DTCluster(event, idx) for idx in range(event['nDtRechitClusters'])]
    self.gLLPs = [gLLP(event, idx) for idx in range(event['nGLLP'])]
    self.leptons = [Lepton(event, idx) for idx in range(event['nLeptons'])]
    self.jets = [Jet(event, idx) for idx in range(event['nJets'])]

  def met_cut(self, min_met=50, max_met=np.inf) -> bool:
    pass

  def csc_match_cut(self, in_place=True, inverse=False):
    pass


####################################################################################################


def fill_event_hists(candidate, hists, comment=None):
  if comment is not None:
    comment = '_' + comment


# @numba.njit
# def process_batch(batch, hists=None):
#   for event in batch:
#     cand = LLPEventCandidate(event)
#     # print(cand.met)


def process_event(event, hists=None) -> LLPEventCandidate:
  cand = LLPEventCandidate(event)
  # print(cand.met)
  # print(hists)
  # hists['ncsc_ndt'].Fill(len(cand.csc_clusters), len(cand.dt_clusters))
  return cand


if __name__ == '__main__':
  rt.EnableImplicitMT()

  out_dir = '/home/psimmerl/Documents/CMS/LLP/reports/weekly/feb9/'
  data_dir = '/home/psimmerl/Documents/CMS/LLP/data/raw/'
  MC_file = data_dir + 'ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted.root'
  run3_file = data_dir + 'DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi.root'

  parser = ArgumentParser(prog='LLPAnalysis', description='LLP skimming for Run3')
  parser.add_argument('-f', '--file', default=MC_file, help='ROOT file to run the analysis on.')
  parser.add_argument('-o', '--out', default=out_dir, help='Output directory to write to.')
  parser.add_argument('-t', '--tree', default='MuonSystem', help='TTree to perform the analysis on.')
  parser.add_argument('-d', '--data', default=True, help='Flag to treat files as data or simulation.')
  parser.add_argument('-j', '--jobs', default=16, help='Number of jobs to run.')
  args = parser.parse_args()

  file_in = args.file
  out_dir = args.out
  tree_name = args.tree
  is_MC = not args.data
  jobs = args.jobs

  pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)  # make out directory if it doesn't exist

  hists = {
      '': TH1F('', '', 10, 0, 10),
      'ncsc_ndt': TH2F('ncsc_ndt', ';N CSC Clusters; N DT Clusters', 7, -0.5, 6.5, 7, -0.5, 6.5)
  }

  print(f'Processing tree \'{tree_name}\' from \'{file_in}\' :')
  muon_system = uproot.open(file_in + ':' + tree_name)
  data = muon_system.arrays(branches)

  hists['ncsc_ndt'].Fill(data['nCscRechitClusters'], data['nDtRechitClusters'])

  # prog_bar = tqdm(total=muon_system.num_entries)
  # for event in data:
  #   cand = process_event(event)

  # rdf = RDataFrame(tree_name, file_in)
  # values = []
  # muon_system = rdf.AsNumpy(values)

  # branches = [b.name for b in muon_system.branches]

  # data = muon_system.arrays(library='np')
  # print(data)

  # # cands = np.vectorize(partial(process_event, hists=hists), data)
  # for batch in muon_system.iterate(step_size=10_000, library='np'):  #'5MB'):  # batch is an Awkward array
  # print(type(batch))
  # # data = map(partial(process_event, hists=hists), [b for b in batch])

  # # data = awk.vectorize(process_event)(batch)
  # # cands = p.map(partial(process_event, hists=hists), batch.to_list())

  # # pool = uproot.ThreadPoolExecutor(jobs)
  # # data = pool.submit(partial(process_event, hists=hists), batch).result()
  # # print(cands)

  # prog_bar.update(1)  #len(batch))

  # prog_bar.close()

  c = TCanvas('c', 'c', 800, 600)
  hists['ncsc_ndt'].Draw('coltext')
  c.Draw()
  c.Print(out_dir + 'ncsc_ndt.png')
  # # ctau vs acceptance TGraph

  # # ctau 95%

  # # background estimation
