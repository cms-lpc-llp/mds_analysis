"""main.py

Main script for the CMS Run3 analysis
"""

import os
import sys
import pathlib

import numpy as np
from argparse import ArgumentParser

import ROOT as rt
from ROOT import TChain, TCanvas, TH1F, TH2F, TLatex, TGraph, RDataFrame
# from src.candidates import LLPCandidate, DTCandidate, CSCCandidate
from src.histo_utilities import create_prof1D, create_TGraph, create_TH1D, create_TH2D, std_color_list
from src import CMS_lumi, tdrstyle
from src.helper_functions import getRecoTime, find_nearest, deltaPhi, deltaR
import uproot
import awkward as awk

from build_rdf import *

####################################################################################################
# branches = [
#     ('met', 'float'),
#     ('metPhi', 'float'),
#     ('MC_condition', 'float'),
#     ('HLTDecision', 'float'),
#     #
#     ('lepE', 'RVec<float>'),
#     ('lepPt', 'RVec<float>'),
#     ('lepEta', 'RVec<float>'),
#     ('lepPhi', 'RVec<float>'),
#     ('lepPdgId', 'RVec<float>'),
#     #
#     ('jetE', 'RVec<float>'),
#     ('jetPt', 'RVec<float>'),
#     ('jetEta', 'RVec<float>'),
#     ('jetPhi', 'RVec<float>'),
#     #
#     ('gLLP_e', 'RVec<float>'),
#     ('gLLP_pt', 'RVec<float>'),
#     ('gLLP_eta', 'RVec<float>'),
#     ('gLLP_phi', 'RVec<float>'),
#     ('gLLP_ctau', 'RVec<float>'),
#     ('gLLP_decay_vertex_r', 'RVec<float>'),
#     ('gLLP_decay_vertex_x', 'RVec<float>'),
#     ('gLLP_decay_vertex_y', 'RVec<float>'),
#     ('gLLP_decay_vertex_z', 'RVec<float>'),
#     #
#     ('cscRechitClusterSize', 'RVec<float>'),
#     ('cscRechitClusterEta', 'RVec<float>'),
#     ('cscRechitClusterPhi', 'RVec<float>'),
#     ('cscRechitClusterMet_dPhi', 'RVec<float>'),
#     ('cscRechitClusterX', 'RVec<float>'),
#     ('cscRechitClusterY', 'RVec<float>'),
#     ('cscRechitClusterZ', 'RVec<float>'),
#     ('cscRechitClusterTimeWeighted', 'RVec<float>'),
#     ('cscRechitClusterTimeSpreadWeightedAll', 'RVec<float>'),
#     ('cscRechitClusterNStation10', 'RVec<float>'),
#     ('cscRechitClusterAvgStation10', 'RVec<float>'),
#     ('cscRechitClusterMaxStation', 'RVec<float>'),
#     ('cscRechitClusterNRechitChamberPlus11', 'RVec<float>'),
#     ('cscRechitClusterNRechitChamberPlus12', 'RVec<float>'),
#     ('cscRechitClusterNRechitChamberMinus11', 'RVec<float>'),
#     ('cscRechitClusterNRechitChamberMinus12', 'RVec<float>'),
#     ('cscRechitCluster_match_gLLP', 'RVec<float>'),
#     ('cscRechitCluster_match_gLLP_eta', 'RVec<float>'),
#     ('cscRechitCluster_match_gLLP_phi', 'RVec<float>'),
#     ('cscRechitCluster_match_gLLP_decay_r', 'RVec<float>'),
#     ('cscRechitCluster_match_gLLP_decay_z', 'RVec<float>'),
#     #
#     ('dtRechitClusterSize', 'RVec<float>'),
#     ('dtRechitClusterEta', 'RVec<float>'),
#     ('dtRechitClusterPhi', 'RVec<float>'),
#     ('dtRechitClusterMet_dPhi', 'RVec<float>'),
#     ('dtRechitClusterX', 'RVec<float>'),
#     ('dtRechitClusterY', 'RVec<float>'),
#     ('dtRechitClusterZ', 'RVec<float>'),
#     ('dtRechitClusterNStation10', 'RVec<float>'),
#     ('dtRechitClusterAvgStation10', 'RVec<float>'),
#     ('dtRechitClusterMaxStation', 'RVec<float>'),
#     ('dtRechitClusterNSegStation1', 'RVec<float>'),
#     ('dtRechitClusterNOppositeSegStation1', 'RVec<float>'),
#     ('dtRechitCluster_match_gLLP', 'RVec<float>'),
#     ('dtRechitCluster_match_gLLP_eta', 'RVec<float>'),
#     ('dtRechitCluster_match_gLLP_phi', 'RVec<float>'),
#     ('dtRechitCluster_match_gLLP_decay_r', 'RVec<float>'),
#     ('dtRechitCluster_match_gLLP_decay_z', 'RVec<float>'),
# ]

# scalar_values = ['pass_met50', 'pass_met200']
# vector_values = ['cscSize50', 'cscSig', 'cscBkg', 'dtSize50', 'dtSig', 'dtBkg']

# @rt.Numba.Declare([b[1] for b in branches], 'RVec<float>')
# def make_scalar_cols(met, metPhi, MC_condition, HLTDecision, lepE, lepPt, lepEta, lepPhi, lepPdgId, jetE, jetPt, jetEta,
#                      jetPhi, gLLP_e, gLLP_pt, gLLP_eta, gLLP_phi, gLLP_ctau, gLLP_decay_vertex_r, gLLP_decay_vertex_x,
#                      gLLP_decay_vertex_y, gLLP_decay_vertex_z, cscSize, cscEta, cscPhi, cscMet_dPhi, cscX, cscY, cscZ,
#                      cscTimeWeighted, cscTimeSpreadWeightedAll, cscNStation10, cscAvgStation10, cscMaxStation,
#                      cscNRechitChamberPlus11, cscNRechitChamberPlus12, cscNRechitChamberMinus11,
#                      cscNRechitChamberMinus12, csc_match_gLLP, csc_match_gLLP_eta, csc_match_gLLP_phi,
#                      csc_match_gLLP_decay_r, csc_match_gLLP_decay_z, dtSize, dtEta, dtPhi, dtMet_dPhi, dtX, dtY, dtZ,
#                      dtNStation10, dtAvgStation10, dtMaxStation, dtNSegStation1, dtNOppositeSegStation1, dt_match_gLLP,
#                      dt_match_gLLP_eta, dt_match_gLLP_phi, dt_match_gLLP_decay_r, dt_match_gLLP_decay_z):

#   pass_met50 = met > 50
#   pass_met200 = met > 200

#   return [pass_met50, pass_met200]

# @rt.Numba.Declare([b[1] for b in branches], 'RVec<RVec<float>>')
# def make_vector_cols(met, metPhi, MC_condition, HLTDecision, lepE, lepPt, lepEta, lepPhi, lepPdgId, jetE, jetPt, jetEta,
#                      jetPhi, gLLP_e, gLLP_pt, gLLP_eta, gLLP_phi, gLLP_ctau, gLLP_decay_vertex_r, gLLP_decay_vertex_x,
#                      gLLP_decay_vertex_y, gLLP_decay_vertex_z, cscSize, cscEta, cscPhi, cscMet_dPhi, cscX, cscY, cscZ,
#                      cscTimeWeighted, cscTimeSpreadWeightedAll, cscNStation10, cscAvgStation10, cscMaxStation,
#                      cscNRechitChamberPlus11, cscNRechitChamberPlus12, cscNRechitChamberMinus11,
#                      cscNRechitChamberMinus12, csc_match_gLLP, csc_match_gLLP_eta, csc_match_gLLP_phi,
#                      csc_match_gLLP_decay_r, csc_match_gLLP_decay_z, dtSize, dtEta, dtPhi, dtMet_dPhi, dtX, dtY, dtZ,
#                      dtNStation10, dtAvgStation10, dtMaxStation, dtNSegStation1, dtNOppositeSegStation1, dt_match_gLLP,
#                      dt_match_gLLP_eta, dt_match_gLLP_phi, dt_match_gLLP_decay_r, dt_match_gLLP_decay_z):

#   cscSize50 = cscSize > 50
#   cscSig = cscSize > 100
#   cscBkg = cscSize < 100

#   dtSize50 = dtSize > 50
#   dtSig = dtSize > 100
#   dtBkg = dtSize < 100

#   return [cscSize50, cscSig, cscBkg, dtSize50, dtSig, dtBkg]

#####

# @rt.Numba.Declare(['RVec<bool>'], 'int')
# def _fix_n_col(cut_col):
#   return np.sum(cut_col)

# @rt.Numba.Declare(['RVec<float>', 'RVec<bool>'], 'RVec<float>')
# def _fix_vector_col_f(col, cut_col):
#   return col[cut_col]

# @rt.Numba.Declare(['RVec<int>', 'RVec<bool>'], 'RVec<int>')
# def _fix_vector_col_i(col, cut_col):
#   return col[cut_col]

# @rt.Numba.Declare(['RVec<bool>', 'RVec<bool>'], 'RVec<bool>')
# def _fix_vector_col_b(col, cut_col):
#   return col[cut_col]

# def update_cols(targ_rdf, cut: str, cut_col: str, col_type: str, n_cluster: int = 0):
#   if cut_col in targ_rdf.GetColumnNames():
#     targ_rdf = targ_rdf.Redefine(cut_col, cut)
#   else:
#     targ_rdf = targ_rdf.Define(cut_col, cut)

#   if col_type == 'jet':
#     n_col_name, cut_prefix = 'nJets', 'jet'
#   elif col_type == 'lep':
#     n_col_name, cut_prefix = 'nLeptons', 'lep'
#   elif col_type == 'csc':
#     n_col_name, cut_prefix = 'nCscRechitClusters', 'csc'  #RechitCluster'
#   elif col_type == 'dt':
#     n_col_name, cut_prefix = 'nDtRechitClusters', 'dt'  #RechitCluster'
#   else:
#     raise ValueError(col_type)

#   bad_cols = ('cscRechitCluster2ZLep1LooseIso', 'dtRechitCluster2ZLep1LooseIso', 'dtRechitCluster22ZLep1LooseIso',
#               'dtRechitClusterJetVetoPtJESDown', 'dtRechitClusterTightJetVetoPtJESDown')

#   for col in targ_rdf.GetColumnNames():
#     isvec = 'RVec' in targ_rdf.GetColumnType(col)
#     ctype = targ_rdf.GetColumnType(col).split('<')[-1].split('>')[0]
#     if str(col)[:len(cut_prefix)] == cut_prefix and isvec and ctype in ('Float_t', 'Int_t',
#                                                                         'Bool_t') and col not in bad_cols:
#       targ_rdf = targ_rdf.Redefine(col, f'Numba::_fix_vector_col_{ctype[0].lower()}({col}, {cut_col})')

#   # print(cut, cut_col, col_type, cut_prefix, n_col_name, targ_rdf.Sum(n_col_name).GetValue())
#   targ_rdf = targ_rdf.Redefine(n_col_name, f'Numba::_fix_n_col({cut_col})')
#   # print('\t', targ_rdf.Sum(n_col_name).GetValue())
#   # print('')

#   if n_cluster:
#     targ_rdf = targ_rdf.Filter(f'(nCscRechitClusters + nDtRechitClusters) == {n_cluster}')

#   return targ_rdf

####################################################################################################


def my_define(rdf, out_col, in_col):
  if out_col in rdf.GetColumnNames():
    rdf = rdf.Redefine(out_col, in_col)
  else:
    rdf = rdf.Define(out_col, in_col)
  return rdf


def add_data_columns(rdf):
  rdf = my_define(rdf, 'nRechitClusters', 'nCscRechitClusters + nDtRechitClusters')
  for v in ('z', 'r'):
    for det in ('cscRechitCluster_match_gLLP_decay', 'dtRechitCluster_match_gLLP_decay', 'gLLP_decay_vertex'):
      rdf = my_define(rdf, f'{det}_abs_{v}', f'abs({det}_{v})')
  return rdf


def apply_event_cuts(rdf):
  rdf = rdf.Filter('nLeptons == 0 && nJets == 0')
  rdf = add_data_columns(rdf)
  return rdf


def apply_cluster_cuts(rdf, n_clusters=0):
  rdf = apply_size_cut(rdf, col_type='csc')
  rdf = apply_time_cut(rdf, col_type='csc')
  rdf = apply_time_spread_cut(rdf, col_type='csc')

  rdf = apply_size_cut(rdf, col_type='dt')
  # rdf = apply_time_cut(rdf, col_type='dt') #! No time column in ntuple
  # rdf = apply_time_spread_cut(rdf, col_type='dt') #! No time column in ntuple

  rdf = apply_cbid_cut(rdf)
  rdf = apply_endcap_nrechit_cut(rdf)
  rdf = apply_barrel_nrechit_cut(rdf)

  rdf = add_data_columns(rdf).Filter(f'nRechitClusters > {n_clusters}')
  return rdf


def canvas(nrows=1, ncols=1, width=None, height=None, name='c', grid=True):
  if width is None:
    width = ncols * 1000
  if height is None:
    height = nrows * 800
  c = rt.TCanvas(name, name, width, height)
  c.Clear()
  if ncols * nrows > 1:
    c.Divide(ncols, nrows)
    if grid:
      for i in range(nrows * ncols):
        c.cd(i + 1).SetGrid()
    c.cd(1)
  else:
    if grid:
      c.SetGrid()

  c.Draw()
  return c


def weight_calc(llp_ct, new_ctau, old_ctau):
  if llp_ct.ndim > 1:
    llp_ct = np.array(np.sum(llp_ct, axis=1))
  source = np.exp(-1.0 * llp_ct / old_ctau) / old_ctau**2
  weight = 1.0 / new_ctau**2 * np.exp(-1.0 * llp_ct / new_ctau) / source
  return weight


def land(*args):
  out = args[0]
  for v in args[1:]:
    out = np.logical_and(out, v)
  return out


def lor(*args):
  out = args[0]
  for v in args[1:]:
    out = np.logical_or(out, v)
  return out


def lxor(*args):
  out = args[0]
  for v in args[1:]:
    out = np.logical_xor(out, v)
  return out


def lnot(arg):
  return np.logical_not(arg)


if __name__ == '__main__':
  rt.gROOT.SetBatch()
  out_dir = '/home/psimmerl/Documents/CMS/LLP/reports/weekly/feb9/'
  data_dir = '/home/psimmerl/Documents/CMS/LLP/data/raw/'
  MC_file = data_dir + 'ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted.root'
  run3_file = data_dir + 'DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi.root'

parser = ArgumentParser(prog='LLPAnalysis', description='LLP skimming for Run3')
parser.add_argument('-f', '--file', default=MC_file, help='ROOT file to run the analysis on.')
parser.add_argument('-o', '--out', default=out_dir, help='Output directory to write to.')
parser.add_argument('-t', '--tree', default='MuonSystem', help='TTree to perform the analysis on.')
parser.add_argument('-d', '--data', default=True, help='Flag to treat files as data or simulation.')
parser.add_argument('-j', '--jobs', default=0, help='Number of jobs to run.')
args = parser.parse_args()

file_in = args.file
out_dir = args.out
tree_name = args.tree
is_MC = not args.data
jobs = args.jobs

  pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)  # make out directory if it doesn't exist

  rt.EnableImplicitMT(jobs)

  print(f'Processing tree \'{tree_name}\' from \'{file_in}\' :')
  rdf = RDataFrame(tree_name, file_in)
  # rdf.AsNumpy([])

  # rdf = rdf.Declare('scalar_values', 'Numba::make_scalar_cols(' + ','.join([b[0] for b in branches] + ')'))
  # rdf = rdf.Declare('vector_values', 'Numba::make_vector_cols(' + ','.join([b[0] for b in branches] + ')'))

  # for i, v in enumerate(scalar_values):
  #   rdf = rdf.Define(v, f'scalar_values[{i}]')
  # for i, v in enumerate(vector_values):
  #   rdf = rdf.Define(v, f'vector_values[{i}]')

  # ctau vs acceptance TGraph, https://github.com/cms-lpc-llp/delayed_jet_analyzer/blob/master/plot_scripts/MuonSystem_Analysis/combination/2tag/MuonSystem_plots-accp_vs_ctau.ipynb

  a = tdrstyle.setTDRStyle()
  CMS_lumi.writeExtraText = 0
  met = True
  mass, ctaus = ['15'], ['1000']
  plot_ctaus = ['10', '100', '200', '500', '700', '1000', '2000', '5000', '10000', '50000', '100000', '1000000']

  # rdf = rdf.Define('gLLP_decay_vertex_abs_r', 'abs(gLLP_decay_vertex_r)')
  # rdf = rdf.Define('gLLP_decay_vertex_abs_z', 'abs(gLLP_decay_vertex_z)')
  # rdf = rdf.Define('gLLP_abs_eta', 'abs(gLLP_eta)')
  # data = rdf.AsNumpy([
  #     'gLLP_decay_vertex_abs_r',
  #     'gLLP_decay_vertex_abs_z',
  #     'gLLP_abs_eta',
  #     'weight',
  #     'gHiggsPt',
  # ])

  # for k, v in data.items():
  #   if not isinstance(v[0], (int, float, bool)):
  #     print('fixing k')
  #     data = np.array([np.asarray(vv) for vv in v])

  muon_system = uproot.open(file_in + ':' + tree_name)
  # data = muon_system.arrays(
  #     ['gLLP_decay_vertex_r', 'gLLP_decay_vertex_z', 'gLLP_eta', 'weight', 'gHiggsPt', 'gLLP_ctau', 'HLTDecision'],
  #     library='ak')

  data = muon_system.arrays(
      ['gLLP_decay_vertex_r', 'gLLP_decay_vertex_z', 'gLLP_eta', 'weight', 'gHiggsPt', 'gLLP_ctau'], library='ak')

  data = muon_system.arrays([
      'runNum', 'MC_condition', 'lumiSec', 'evtNum', 'mH', 'mX', 'ctau', 'npv', 'npu', 'weight', 'pileupWeight',
      'Flag_HBHENoiseFilter', 'Flag_BadPFMuonFilter', 'Flag_HBHEIsoNoiseFilter', 'Flag_CSCTightHaloFilter',
      'Flag_globalSuperTightHalo2016Filter', 'Flag_goodVertices', 'Flag_ecalBadCalibFilter',
      'Flag_BadChargedCandidateFilter', 'Flag_eeBadScFilter', 'Flag_all', 'Flag2_HBHENoiseFilter',
      'Flag2_HBHEIsoNoiseFilter', 'Flag2_BadPFMuonFilter', 'Flag2_globalSuperTightHalo2016Filter',
      'Flag2_globalTightHalo2016Filter', 'Flag2_BadChargedCandidateFilter', 'Flag2_EcalDeadCellTriggerPrimitiveFilter',
      'Flag2_ecalBadCalibFilter', 'Flag2_eeBadScFilter', 'Flag2_all', 'rho', 'met', 'metPhi', 'gHiggsPt', 'gHiggsE',
      'gHiggsEta', 'gHiggsPhi', 'nCscRings', 'nDtRings', 'nCscRechitClusters', 'cscRechitCluster_match_gLLP',
      'cscRechitCluster_match_gLLP_minDeltaR', 'cscRechitCluster_match_gLLP_index', 'cscRechitCluster_match_gLLP_eta',
      'cscRechitCluster_match_gLLP_phi', 'cscRechitCluster_match_gLLP_decay_r', 'cscRechitCluster_match_gLLP_decay_z',
      'cscRechitCluster_match_gLLP_csc', 'cscRechitCluster_match_gLLP_dt', 'cscRechitCluster_match_gLLP_e',
      'cscRechitClusterX', 'cscRechitClusterY', 'cscRechitClusterZ', 'cscRechitClusterTimeWeighted',
      'cscRechitClusterTimeSpreadWeightedAll', 'cscRechitClusterPhi', 'cscRechitClusterEta',
      'cscRechitClusterJetVetoPt', 'cscRechitClusterJetVetoLooseId', 'cscRechitClusterJetVetoTightId',
      'cscRechitClusterJetVetoE', 'cscRechitClusterMuonVetoPt', 'cscRechitClusterMuonVetoE',
      'cscRechitClusterMuonVetoLooseId', 'cscRechitClusterMuonVetoGlobal', 'cscRechitCluster_match_dtSeg_0p4',
      'cscRechitCluster_match_MB1Seg_0p4', 'cscRechitCluster_match_RE12_0p4', 'cscRechitCluster_match_RB1_0p4',
      'cscRechitClusterSize', 'cscRechitClusterNStation10', 'cscRechitClusterAvgStation10',
      'cscRechitClusterMaxStation', 'cscRechitClusterMaxStationRatio', 'cscRechitClusterNChamber',
      'cscRechitClusterMaxChamber', 'cscRechitClusterMaxChamberRatio', 'cscRechitClusterNRechitChamberPlus11',
      'cscRechitClusterNRechitChamberPlus12', 'cscRechitClusterNRechitChamberPlus13',
      'cscRechitClusterNRechitChamberPlus21', 'cscRechitClusterNRechitChamberPlus22',
      'cscRechitClusterNRechitChamberPlus31', 'cscRechitClusterNRechitChamberPlus32',
      'cscRechitClusterNRechitChamberPlus41', 'cscRechitClusterNRechitChamberPlus42',
      'cscRechitClusterNRechitChamberMinus11', 'cscRechitClusterNRechitChamberMinus12',
      'cscRechitClusterNRechitChamberMinus13', 'cscRechitClusterNRechitChamberMinus21',
      'cscRechitClusterNRechitChamberMinus22', 'cscRechitClusterNRechitChamberMinus31',
      'cscRechitClusterNRechitChamberMinus32', 'cscRechitClusterNRechitChamberMinus41',
      'cscRechitClusterNRechitChamberMinus42', 'cscRechitClusterMet_dPhi', 'nDtRechitClusters',
      'dtRechitClusterNSegStation1', 'dtRechitClusterNSegStation2', 'dtRechitClusterNSegStation3',
      'dtRechitClusterNSegStation4', 'dtRechitClusterNOppositeSegStation1', 'dtRechitClusterNOppositeSegStation2',
      'dtRechitClusterNOppositeSegStation3', 'dtRechitClusterNOppositeSegStation4', 'dtRechitCluster_match_gLLP',
      'dtRechitCluster_match_gLLP_minDeltaR', 'dtRechitCluster_match_gLLP_index', 'dtRechitCluster_match_gLLP_eta',
      'dtRechitCluster_match_gLLP_phi', 'dtRechitCluster_match_gLLP_decay_r', 'dtRechitCluster_match_gLLP_decay_z',
      'dtRechitCluster_match_gLLP_csc', 'dtRechitCluster_match_gLLP_dt', 'dtRechitCluster_match_gLLP_e',
      'dtRechitClusterX', 'dtRechitClusterY', 'dtRechitClusterZ', 'dtRechitClusterWheel', 'dtRechitClusterPhi',
      'dtRechitClusterEta', 'dtRechitClusterJetVetoPt', 'dtRechitClusterJetVetoLooseId',
      'dtRechitClusterJetVetoTightId', 'dtRechitClusterJetVetoE', 'dtRechitClusterMuonVetoPt',
      'dtRechitClusterMuonVetoE', 'dtRechitClusterMuonVetoTightId', 'dtRechitClusterMuonVetoLooseId',
      'dtRechitClusterMuonVetoGlobal', 'dtRechitClusterOverlap', 'dtRechitClusterSize', 'dtRechitClusterNoiseHit',
      'dtRechitClusterNoiseHitStation1', 'dtRechitClusterNoiseHitStation2', 'dtRechitClusterNoiseHitStation3',
      'dtRechitClusterNoiseHitStation4', 'dtRechitClusterNStation10', 'dtRechitClusterAvgStation10',
      'dtRechitClusterMaxStation', 'dtRechitClusterMaxStationRatio', 'dtRechitClusterNChamber',
      'dtRechitClusterMaxChamber', 'dtRechitClusterMaxChamberRatio', 'dtRechitClusterNHitStation1',
      'dtRechitClusterNHitStation2', 'dtRechitClusterNHitStation3', 'dtRechitClusterNHitStation4',
      'dtRechitClusterMet_dPhi', 'dtRechitCluster_match_RPChits_dPhi0p5', 'dtRechitCluster_match_RPCBx_dPhi0p5',
      'dtRechitCluster_match_RB1_0p4', 'dtRechitCluster_match_RB1_dPhi0p5', 'dtRechitCluster_match_MB1Seg_0p4',
      'dtRechitCluster_match_MB1Seg_0p5', 'dtRechitCluster_match_MB1hits_0p4', 'dtRechitCluster_match_MB1hits_0p5',
      'dtRechitCluster_match_MB1hits_cosmics_plus', 'dtRechitCluster_match_MB1hits_cosmics_minus', 'nGLLP', 'gLLP_eta',
      'gLLP_phi', 'gLLP_csc', 'gLLP_dt', 'gLLP_beta', 'gLLP_e', 'gLLP_pt', 'gLLP_ctau', 'gLLP_decay_vertex_r',
      'gLLP_decay_vertex_x', 'gLLP_decay_vertex_y', 'gLLP_decay_vertex_z', 'nLeptons', 'lepE', 'lepPt', 'lepEta',
      'lepPhi', 'lepPdgId', 'lepDZ', 'lepTightId', 'lepPassLooseIso', 'lepPassTightIso', 'lepPassVTightIso',
      'lepPassVVTightIso', 'nJets', 'jetE', 'jetPt', 'jetEta', 'jetPhi', 'jetTightPassId', 'HLTDecision'
  ])
  data['gLLP_decay_vertex_abs_r'] = np.abs(data['gLLP_decay_vertex_r'])
  data['gLLP_decay_vertex_abs_z'] = np.abs(data['gLLP_decay_vertex_z'])
  data['gLLP_abs_eta'] = np.abs(data['gLLP_eta'])

  # for k, v in data.items():
  #   print(k, type(v), v.shape)

  sel_dt = land(data['gLLP_decay_vertex_abs_z'] < 661, data['gLLP_decay_vertex_abs_r'] > 380,
                data['gLLP_decay_vertex_r'] < 738)

  sel_csc = land(data['gLLP_abs_eta'] < 2.4, data['gLLP_decay_vertex_abs_r'] < 695.5,
                 data['gLLP_decay_vertex_z'] < 1100, data['gLLP_decay_vertex_abs_z'] > 661)

  me11 = land(data['gLLP_abs_eta'] < 2.4, data['gLLP_decay_vertex_abs_r'] < 270, data['gLLP_decay_vertex_abs_z'] < 661,
              data['gLLP_decay_vertex_abs_z'] > 500)

  sel_csc = lor(sel_csc, me11)  # definition to make sure no overlap, BUG?
  sel_ev_cross = land(np.sum(sel_csc, axis=1) == 1, np.sum(sel_dt, axis=1) == 1)
  sel_ev_csc2 = np.sum(sel_csc, axis=1) == 2
  sel_ev_dt2 = np.sum(sel_dt, axis=1) == 2
  sel_ev_csc1 = np.sum(sel_csc, axis=1) >= 1
  sel_ev_dt1 = np.sum(sel_dt, axis=1) >= 1
  sel_ev_2tag = lor(sel_ev_cross, sel_ev_csc2, sel_ev_dt2)

  accep_dt1 = np.sum(data['weight'][sel_ev_dt1]) / np.sum(data['weight'])
  accep_csc1 = np.sum(data['weight'][sel_ev_csc1]) / np.sum(data['weight'])
  accep_dt2 = np.sum(data['weight'][sel_ev_dt2]) / np.sum(data['weight'])
  accep_csc2 = np.sum(data['weight'][sel_ev_csc2]) / np.sum(data['weight'])
  accep_cross = np.sum(data['weight'][sel_ev_cross]) / np.sum(data['weight'])

  accep_2tag = accep_dt2 + accep_csc2 + accep_cross

  sel_ev = data['met'] > 200  #np.sum(data['HLTDecision'][:, 562:771],
  # axis=1)  #  #(data['metEENoise'] >= 200), (data['METNoMuTrigger'])
  # sel_ev = data['HLTDecision']  #(data['metEENoise'] >= 200), (data['METNoMuTrigger'])

  print('weight', np.sum(sel_ev_dt1), '\n')
  accep_met_dt1 = np.sum(data['weight'][land(sel_ev, sel_ev_dt1)]) / np.sum(data['weight'])
  accep_met_csc1 = np.sum(data['weight'][land(sel_ev, sel_ev_csc1)]) / np.sum(data['weight'])
  accep_met_dt2 = np.sum(data['weight'][land(sel_ev, sel_ev_dt2)]) / np.sum(data['weight'])
  accep_met_csc2 = np.sum(data['weight'][land(sel_ev, sel_ev_csc2)]) / np.sum(data['weight'])
  accep_met_cross = np.sum(data['weight'][land(sel_ev, sel_ev_cross)]) / np.sum(data['weight'])

  accep_met_2tag = accep_met_dt2 + accep_met_csc2 + accep_met_cross

  variables = [accep_dt1, accep_csc1, accep_dt2, accep_csc2, accep_cross]
  if met:
    variables = [accep_met_dt1, accep_met_csc1, accep_met_dt2, accep_met_csc2, accep_met_cross]
  selections = [sel_ev_dt1, sel_ev_csc1, sel_ev_dt2, sel_ev_csc2, sel_ev_cross]
  names = ['1 tag DT', '1 tag CSC', '2 tag DT', '2 tag CSC', '1 tag DT + 1 tag CSC']
  names = ['1_tagDT', '1_tag_CSC', '2_tag_DT', '2_tag_CSC', '1_DT_1_CSC']
  names = ['dt1', 'csc1', 'dt2', 'csc2', 'dt1csc1']

  print(accep_2tag, accep_met_2tag)

  c = canvas()
  for j, var in enumerate(variables):
    # if not j == 2:
    #   continue
    leg = rt.TLegend(0.6, 0.2, 0.9, 0.35)
    leg.SetTextSize(0.025)
    leg.SetBorderSize(0)
    leg.SetEntrySeparation(0.01)
    c.SetRightMargin(0.04)

    rt.gStyle.SetOptFit(1011)
    h = {}

    for i, m in enumerate(mass):
      x = []
      y = []
      for ct in plot_ctaus:
        if ct in ctaus:
          x.append(int(ct))
          y.append(var)  #['m' + m + 'ctau' + ct])
        else:

          ctf = int(ct)
          ct_list = 10**int(np.log10(ctf))

          if ctf < int(ctaus[0]):
            ct_list = [int(ctaus[0])]
          elif ctf > int(ctaus[-1]):
            ct_list = [int(ctaus[-1])]
          elif ct in ctaus:
            ct_list = [ctf]
          else:
            ct_list = [ct_list, ct_list * 10]
          sig = 0
          for l, ct0 in enumerate(ct_list):

            k = 'm' + m + 'ctau' + str(ct0)
            if met:
              cond = land(sel_ev, selections[j])
            else:
              cond = selections[j]
            gLLP_ctau = data['gLLP_ctau'][cond]

            weight_ctau = weight_calc(gLLP_ctau, ctf / 10, int(ct0) / 10)  # convert everything to cm
            gLLP_ctau = np.sum(data['gLLP_ctau'], axis=1)[cond]
            if len(ct_list) == 1:
              weight_cond = gLLP_ctau >= 0
            else:
              if l == 0:
                weight_cond = gLLP_ctau < int(ct_list[0] / 2)
              else:
                weight_cond = gLLP_ctau >= int(ct_list[0] / 2)
            sig += np.sum((data['weight'])[cond][weight_cond] * weight_ctau[weight_cond])
            print(k, np.sum(data['weight']))
          x.append(int(ct))
          y.append(sig / np.sum(data['weight']))

      x = np.array(x)
      y = np.array(y)
      print(x, y)
      cond = y > 0.0
      h[m + str(j)] = create_TGraph(x[cond] / 1000.0, y[cond], axis_title=['c#tau [m]', 'acceptance'])

      h[m + str(j)].SetLineColor(std_color_list[i])
      h[m + str(j)].SetLineWidth(3)
      h[m + str(j)].GetXaxis().SetTitleOffset(1)
      h[m + str(j)].GetYaxis().SetTitleSize(0.05)
      h[m + str(j)].GetYaxis().SetTitleOffset(1.5)
      h[m + str(j)].GetXaxis().SetLimits(0.01, 1000.0)

      h[m + str(j)].GetYaxis().SetRangeUser(1e-7, 2)
      leg.AddEntry(h[m + str(j)], "m_{x} = " + str(m) + " GeV", "L")

    for i, m in enumerate(h.keys()):
      h[m].Draw('CA' if i == 0 else 'Csame')

    # with rt.TFile(out_dir + 'out.root', 'CREATE') as fout:
    #   for k, v in h.items():
    #     k = k + ('_MET' if met else '')
    #     fout.write(k, v)

    leg.Draw()
    c.SetLogy(True)
    c.SetLogx(True)
    c.SetTicky(1)
    c.SetTickx(1)

  c.Print(out_dir + 'ctau_accp' + ('_MET' if met else '') + '.png')

  #########################################################3
  c = canvas()
  names = ['1 tag DT', '1 tag CSC', '2 tag total', '2 tag DT', '2 tag CSC', '1 tag DT + 1 tag CSC']
  variables = [accep_dt1, accep_csc1, accep_2tag, accep_dt2, accep_csc2, accep_cross]
  if met:
    variables = [accep_met_dt1, accep_met_csc1, accep_met_2tag, accep_met_dt2, accep_met_csc2, accep_met_cross]
  selections = [sel_ev_dt1, sel_ev_csc1, sel_ev_2tag, sel_ev_dt2, sel_ev_csc2, sel_ev_cross]

  for j, m in enumerate(mass):

    if not m == '15':
      continue
    c.SetRightMargin(0.04)
    rt.gStyle.SetOptFit(1011)
    leg = rt.TLegend(0.6, 0.2, 0.9, 0.35)
    leg.SetTextSize(0.025)
    leg.SetBorderSize(0)
    leg.SetEntrySeparation(0.01)
    h = {}
    for i, var in enumerate(variables):
      x = []
      y = []
      for ct in plot_ctaus:
        if ct in ctaus:
          x.append(int(ct))
          y.append(var)  #['m' + m + 'ctau' + ct])
        else:

          ctf = int(ct)
          ct_list = 10**int(np.log10(ctf))

          if ctf < int(ctaus[0]):
            ct_list = [int(ctaus[0])]
          elif ctf > int(ctaus[-1]):
            ct_list = [int(ctaus[-1])]
          elif ct in ctaus:
            ct_list = [ctf]
          else:
            ct_list = [ct_list, ct_list * 10]
          sig = 0
          for l, ct0 in enumerate(ct_list):

            k = 'm' + m + 'ctau' + str(ct0)
            if met:
              cond = land(sel_ev, selections[i])
            else:
              cond = selections[i]
            gLLP_ctau = data['gLLP_ctau'][cond]

            weight_ctau = weight_calc(gLLP_ctau, ctf / 10, int(ct0) / 10)  # convert everything to cm
            gLLP_ctau = np.sum(data['gLLP_ctau'], axis=1)[cond]
            if len(ct_list) == 1:
              weight_cond = gLLP_ctau >= 0
            else:
              if l == 0:
                weight_cond = gLLP_ctau < int(ct_list[0] / 2)
              else:
                weight_cond = gLLP_ctau >= int(ct_list[0] / 2)
            sig += np.sum((data['weight'])[cond][weight_cond] * weight_ctau[weight_cond])
            print(k, np.sum(data['weight']))
          x.append(int(ct))
          y.append(sig / np.sum(data['weight']))

      x = np.array(x)
      y = np.array(y)
      cond = y > 0.0
      print(x, y)
      h[m + str(i)] = create_TGraph(x[cond] / 1000.0, y[cond], axis_title=['c#tau [m]', 'acceptance'])
      h[m + str(i)].SetLineColor(std_color_list[i])
      h[m + str(i)].SetLineWidth(3)
      h[m + str(i)].GetXaxis().SetTitleOffset(1)
      h[m + str(i)].GetYaxis().SetTitleSize(0.05)
      h[m + str(i)].GetYaxis().SetTitleOffset(1.5)
      h[m + str(i)].GetXaxis().SetLimits(0.01, 1000.0)
      h[m + str(j)].GetYaxis().SetRangeUser(1e-8, 1)
      leg.AddEntry(h[m + str(i)], names[i], "L")

    for i, k in enumerate(h.keys()):
      h[k].Draw('CA' if i == 0 else 'Csame')

    # with rt.TFile(out_dir + 'out.root', 'CREATE') as fout:
    #   for k, v in h.items():
    #     k = k + '_multitags' + ('_MET' if met else '')
    #     fout.write(k, v)

    leg.Draw()
    c.SetLogy()
    c.SetLogx()
    c.SetTicky(1)
    c.SetTickx(1)

  c.Print(out_dir + 'ctau_accp_multitags' + ('_MET' if met else '') + '.png')

  ###############################################################
  # backgro
