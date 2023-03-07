"""main.py

Main script for the CMS Run3 analysis
"""

import os
import sys
import pathlib

import numpy as np

import uproot
import awkward as ak
import awkward.numba

import ROOT as rt
from src import CMS_lumi, tdrstyle
from ROOT import TChain, TCanvas, TH1F, TH2F, TF1, TLatex, TGraph, RDataFrame, TLine, TBox
from src.histo_utilities import create_prof1D, create_TGraph, create_TH1D, create_TH2D, std_color_list
from src.helper_functions import (getRecoTime, find_nearest, deltaPhi, deltaR, lor, land, lxor, lnot, asum, canvas,
                                  weight_calc)

from main2 import MuonSystem, get_lat_leg

from main_030223 import H1D, multi_plot, fix_nbranch, match_clusters, pass_NCSC_NDT, pass_muon_veto, pass_jet_veto, pass_in_time, pass_L1

########################################################

if __name__ == '__main__':
  out_dir = '/home/psimmerl/Documents/CMS/LLP/reports/weekly/mar2/'
  data_dir = '/home/psimmerl/Documents/CMS/LLP/data/raw/'

  file_db_0p4 = data_dir + 'ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v4.root'
  file_ca_0p4 = data_dir + 'ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v5.root'
  file_ca_0p6 = data_dir + 'ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v6.root'
  # run3_file = data_dir + 'DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi.root'

  files = [file_db_0p4, file_ca_0p4, file_ca_0p6]
  tree_name = 'MuonSystem'
  ending = '.png'

  dets = ['CSC', 'DT']
  pi = rt.TMath.Pi()
  met = False
  gc = []

  is_mc = 'DisplacedJet' not in files[0]

  rt.gROOT.SetBatch()

  a = tdrstyle.setTDRStyle()
  # CMS_lumi.writeExtraText = 0
  rt.gStyle.SetOptFit(0)  #1011)

  pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)  # make out directory if it doesn't exist

  print(f'Using output directory \'{out_dir}\'')
  print(f'Processing tree \'{tree_name}\':')  # from \'{file_in}\':')
  print(f'\tTreating file as {"MC" if is_mc else "DATA"}')

  fuproots = [uproot.open(ff + ':' + tree_name) for ff in files]
  mss = [fup.arrays() for fup in fuproots]
  mss = [pass_NCSC_NDT(ms, ncsc=1, ndt=0) for ms in mss]
  mss = [fix_nbranch(ms) for ms in mss]
  mss_matched = [match_clusters(ms) for ms in mss]

  # if not is_mc:
  # print('REMOVING ALL EVENTS WITH N CLUSTERS > 2')
  # mss = [ms[(ms['nCscRechitClusters'] + ms['nDtRechitClusters']) <= 2] for ms in mss]

  labs = [
      'DBSCAN',
      'CA 0.4',
      'CA 0.6',
  ]
  c = canvas()
  title = ';gLLP Z Decay Vertex [cm];efficiency'
  bins = (100, 400, 1100)
  nums = [H1D(np.abs(ms['cscRechitCluster_match_gLLP_decay_z']), title, bins=bins) for ms in mss_matched]
  dens = [H1D(np.abs(ms['gLLP_decay_vertex_z']), title, bins=bins) for ms in mss]
  for num, den in zip(nums, dens):
    num.Divide(den)

  _ = multi_plot(nums, labs, legxy=(0.4, 0.2, 0.6, 0.3))
  # CMS_lumi.CMS_lumi(c.cd(), 9999, 10)
  c.Print(out_dir + 'eff_compare_DBSCAN-CA' + ending)

  # if is_mc:
  #   ms = match_clusters(ms)

#   ##########
#   # yapf: disable
#   ms_limited = ms[(ms['nCscRechitClusters'] == 1) & (ms['nDtRechitClusters'] == 1)]
#   ms_limited = apply_cluster_cut(ms_limited, ms_limited['cscRechitClusterSize'] >= 100, 'csc')
#   ms_limited = apply_cluster_cut(ms_limited, (50 <= ms_limited['dtRechitClusterSize']) & (ms_limited['dtRechitClusterSize'] < 60), 'dt')
#   ms_limited = pass_NCSC_NDT(ms_limited)

#   CSC___DT____ = pass_NCSC_NDT(ms_limited)
#   # CSCJ__DT____ = pass_NCSC_NDT(pass_jet_veto(ms_limited,csc=True,dt=False))
#   # CSC___DTJ___ = pass_NCSC_NDT(pass_jet_veto(ms_limited,csc=False,dt=True))
#   CSCJ__DTJ___ = pass_NCSC_NDT(pass_jet_veto(ms_limited,csc=True,dt=True))
#   # CSC_M_DT____ = pass_NCSC_NDT(pass_muon_veto(ms_limited,csc=True,dt=False))
#   # CSC___DT_M__ = pass_NCSC_NDT(pass_muon_veto(ms_limited,csc=False,dt=True))
#   CSC_M_DT_M__ = pass_NCSC_NDT(pass_muon_veto(ms_limited,csc=True,dt=True))
#   # CSCJM_DT____ = pass_NCSC_NDT(pass_muon_veto(CSCJ__DT____,csc=True,dt=False))
#   # CSC___DTJM__ = pass_NCSC_NDT(pass_muon_veto(CSC___DTJ___,csc=False,dt=True))
#   CSCJM_DTJM__ = pass_NCSC_NDT(pass_muon_veto(CSCJ__DTJ___,csc=True,dt=True))
#   # CSC__TDT____ = pass_NCSC_NDT(pass_in_time(ms_limited,csc=True,dt=False))
#   # CSC___DT__T_ = pass_NCSC_NDT(pass_in_time(ms_limited,csc=False,dt=True))
#   CSC__TDT__T_ = pass_NCSC_NDT(pass_in_time(ms_limited,csc=True,dt=True))
#   CSCJMTDTJM__ = pass_NCSC_NDT(pass_in_time(CSCJM_DTJM__,csc=True,dt=False))
#   CSCJM_DTJMT_ = pass_NCSC_NDT(pass_in_time(CSCJM_DTJM__,csc=False,dt=True))
#   CSCJMTDTJMT_ = pass_NCSC_NDT(pass_in_time(CSCJM_DTJM__,csc=True,dt=True))
#   # CSC___DT___1 = pass_NCSC_NDT(pass_L1(ms_limited))
#   # CSCJ__DTJ__1 = pass_NCSC_NDT(pass_L1(CSCJ__DTJ___))
#   # CSCJM_DTJM_1 = pass_NCSC_NDT(pass_L1(CSCJM_DTJM__))
#   # CSC__TDT__T1 = pass_NCSC_NDT(pass_L1(CSC__TDT__T_))
#   # CSCJMTDTJMT1 = pass_NCSC_NDT(pass_L1(CSCJMTDTJMT_))

#   CSC_MTDTJMT_ = pass_NCSC_NDT(pass_jet_veto(pass_muon_veto(CSC__TDT__T_,csc=True,dt=True),csc=False,dt=True))
#   CSCJ_TDTJMT_ = pass_NCSC_NDT(pass_jet_veto(pass_muon_veto(CSC__TDT__T_,csc=False,dt=True),csc=True,dt=True))
#   CSCJMTDT_MT_ = pass_NCSC_NDT(pass_jet_veto(pass_muon_veto(CSC__TDT__T_,csc=True,dt=True),csc=True,dt=False))
#   CSCJMTDTJ_T_ = pass_NCSC_NDT(pass_jet_veto(pass_muon_veto(CSC__TDT__T_,csc=True,dt=False),csc=True,dt=True))

#   print(f'CSC: ___, DT: ___ = {ak.sum(land(50 <= CSC___DT____["dtRechitClusterSize"], CSC___DT____["dtRechitClusterSize"] < 60)):>7,}')
#   # print(f'CSC: J__, DT: ___ = {ak.sum(land(50 <= CSCJ__DT____["dtRechitClusterSize"], CSCJ__DT____["dtRechitClusterSize"] < 60)):>7,}')
#   # print(f'CSC: ___, DT: J__ = {ak.sum(land(50 <= CSC___DTJ___["dtRechitClusterSize"], CSC___DTJ___["dtRechitClusterSize"] < 60)):>7,}')
#   print(f'CSC: J__, DT: J__ = {ak.sum(land(50 <= CSCJ__DTJ___["dtRechitClusterSize"], CSCJ__DTJ___["dtRechitClusterSize"] < 60)):>7,}')
#   # print(f'CSC: _M_, DT: ___ = {ak.sum(land(50 <= CSC_M_DT____["dtRechitClusterSize"], CSC_M_DT____["dtRechitClusterSize"] < 60)):>7,}')
#   # print(f'CSC: ___, DT: _M_ = {ak.sum(land(50 <= CSC___DT_M__["dtRechitClusterSize"], CSC___DT_M__["dtRechitClusterSize"] < 60)):>7,}')
#   print(f'CSC: _M_, DT: _M_ = {ak.sum(land(50 <= CSC_M_DT_M__["dtRechitClusterSize"], CSC_M_DT_M__["dtRechitClusterSize"] < 60)):>7,}')
#   # print(f'CSC: JM_, DT: ___ = {ak.sum(land(50 <= CSCJM_DT____["dtRechitClusterSize"], CSCJM_DT____["dtRechitClusterSize"] < 60)):>7,}')
#   # print(f'CSC: ___, DT: JM_ = {ak.sum(land(50 <= CSC___DTJM__["dtRechitClusterSize"], CSC___DTJM__["dtRechitClusterSize"] < 60)):>7,}')
#   print(f'CSC: JM_, DT: JM_ = {ak.sum(land(50 <= CSCJM_DTJM__["dtRechitClusterSize"], CSCJM_DTJM__["dtRechitClusterSize"] < 60)):>7,}')
#   # print(f'CSC: __T, DT: ___ = {ak.sum(land(50 <= CSC__TDT____["dtRechitClusterSize"], CSC__TDT____["dtRechitClusterSize"] < 60)):>7,}')
#   # print(f'CSC: ___, DT: __T = {ak.sum(land(50 <= CSC___DT__T_["dtRechitClusterSize"], CSC___DT__T_["dtRechitClusterSize"] < 60)):>7,}')
#   print(f'CSC: __T, DT: __T = {ak.sum(land(50 <= CSC__TDT__T_["dtRechitClusterSize"], CSC__TDT__T_["dtRechitClusterSize"] < 60)):>7,}')
#   # print(f'CSC: JMT, DT: JM_ = {ak.sum(land(50 <= CSCJMTDTJM__["dtRechitClusterSize"], CSCJMTDTJM__["dtRechitClusterSize"] < 60)):>7,}')
#   # print(f'CSC: JM_, DT: JMT = {ak.sum(land(50 <= CSCJM_DTJMT_["dtRechitClusterSize"], CSCJM_DTJMT_["dtRechitClusterSize"] < 60)):>7,}')
#   print(f'CSC: JMT, DT: JMT = {ak.sum(land(50 <= CSCJMTDTJMT_["dtRechitClusterSize"], CSCJMTDTJMT_["dtRechitClusterSize"] < 60)):>7,}')
#   # print('With level 1 trigger (569):')
#   # print(f'CSC: ___, DT: ___ = {ak.sum(land(50 <= CSC___DT___1["dtRechitClusterSize"], CSC___DT___1["dtRechitClusterSize"] < 60)):>7,}')
#   # print(f'CSC: JM_, DT: JM_ = {ak.sum(land(50 <= CSCJM_DTJM_1["dtRechitClusterSize"], CSCJM_DTJM_1["dtRechitClusterSize"] < 60)):>7,}')
#   # print(f'CSC: __T, DT: __T = {ak.sum(land(50 <= CSC__TDT__T1["dtRechitClusterSize"], CSC__TDT__T1["dtRechitClusterSize"] < 60)):>7,}')
#   # print(f'CSC: JMT, DT: JMT = {ak.sum(land(50 <= CSCJMTDTJMT1["dtRechitClusterSize"], CSCJMTDTJMT1["dtRechitClusterSize"] < 60)):>7,}')
#   # yapf: enable
#   ##########
#   mss = [
#       CSC___DT____,
#       # CSC___DT___1,
#       CSCJ__DTJ___,
#       # CSCJ__DTJ__1,
#       # CSC_M_DT_M__,
#       CSCJM_DTJM__,
#       # CSCJM_DTJM_1,
#       # CSC__TDT____,
#       # CSC___DT__T_,
#       # CSC__TDT__T_,
#       CSCJMTDTJMT_,
#       # CSCJMTDTJMT1,
#   ]
#   labels = [
#       '1 CSC + 1 DT',
#       # '1 CSC + 1 DT & L1',
#       'After Jet Vetoes',
#       # 'After Jet Vetoes & L1',
#       # 'After Muon Vetoes',
#       'After Jet+Muon Vetoes',
#       # 'After Jet+Muon Vetoes & L1',
#       # 'CSC In-time',
#       # 'DT In-time',
#       # 'CSC/DT In-time',
#       'After Jet+Muon Vetoes & In-time',
#       # 'After Jet+Muon Vetoes & In-time & L1',
#   ]

#   c = canvas()
#   c.SetLogy()
#   hhs = [H1D(s['dtRechitClusterSize'], ';DT Cluster Size;count', (10, 50, 60), ymin=99) for s in mss]
#   _ = multi_plot(hhs, labels, ymax_mult=1.5, legxy=(0.6, 0.4, 0.9, 0.6))
#   c.Print(out_dir + 'expfit' + ending)

#   c = canvas()
#   c.SetLogy()
#   lat, leg = get_lat_leg((0.6, 0.85, 0.9, 0.9))
#   ymax = max([h.GetMaximum() for h in hhs[-2:]]) * 1.05
#   hhs[-1].SetMinimum(100)  #ymax)
#   hhs[-2].SetMinimum(100)  #ymax)
#   hhs[-1].SetMaximum(10000)  #ymax)
#   hhs[-2].SetMaximum(10000)  #ymax)
#   hhs[-1].Draw()
#   hhs[-2].Draw('same')
#   leg.AddEntry(hhs[-1], labels[-1], 'L')
#   leg.AddEntry(hhs[-2], labels[-2], 'L')

#   c.Print(out_dir + 'expfit_jmt' + ending)

#   ##########

#   _ms = ms_limited
#   vals1 = [
#       _ms['cscRechitClusterJetVetoLooseId'],
#       _ms['cscRechitClusterJetVetoPt'],
#       _ms['cscRechitClusterMuonVetoLooseId'],
#       np.abs(_ms['cscRechitClusterEta']),
#       _ms['cscRechitClusterNRechitChamberPlus11'] + _ms['cscRechitClusterNRechitChamberMinus11'],
#       _ms['cscRechitClusterNRechitChamberMinus12'] + _ms['cscRechitClusterNRechitChamberPlus12'],
#       _ms['cscRechitClusterTimeWeighted'],
#       _ms['dtRechitClusterJetVetoLooseId'],
#       _ms['dtRechitClusterJetVetoPt'],
#       _ms['dtRechitClusterMuonVetoLooseId'],
#       _ms['dtRechitClusterNSegStation1'],
#       _ms['dtRechitCluster_match_RPCBx_dPhi0p5'],
#   ]
#   vals2 = [
#       CSC_MTDTJMT_['cscRechitClusterJetVetoLooseId'],
#       CSC_MTDTJMT_['cscRechitClusterJetVetoPt'],
#       CSCJ_TDTJMT_['cscRechitClusterMuonVetoLooseId'],
#       np.abs(CSCJ_TDTJMT_['cscRechitClusterEta']),
#       CSCJ_TDTJMT_['cscRechitClusterNRechitChamberPlus11'] + CSCJ_TDTJMT_['cscRechitClusterNRechitChamberMinus11'],
#       CSCJ_TDTJMT_['cscRechitClusterNRechitChamberMinus12'] + CSCJ_TDTJMT_['cscRechitClusterNRechitChamberPlus12'],
#       CSCJM_DTJMT_['cscRechitClusterTimeWeighted'],
#       CSCJMTDT_MT_['dtRechitClusterJetVetoLooseId'],
#       CSCJMTDT_MT_['dtRechitClusterJetVetoPt'],
#       CSCJMTDTJ_T_['dtRechitClusterMuonVetoLooseId'],
#       CSCJMTDTJ_T_['dtRechitClusterNSegStation1'],
#       CSCJMTDTJM__['dtRechitCluster_match_RPCBx_dPhi0p5'],
#   ]
#   labs = [
#       'CSC JetVetoLooseId',
#       'CSC JetVetoPt [GeV]',
#       'CSC MuonVetoLooseId',
#       'CSC |#eta|',
#       'CSC NRechit M11/P11',
#       'CSC NRechit M12/P12',
#       'CSC Time [ns]',
#       'DT JetVetoLooseId',
#       'DT JetVetoPt [GeV]',
#       'DT MuonVetoLooseId',
#       'DT NRechit MB1',
#       'DT RPCBx',
#   ]
#   bins = [
#       (2, -0.5, 1.5),
#       (50, -0.5, 100.5),
#       (2, -0.5, 1.5),
#       (20, 0, 3),
#       (21, -0.5, 20.5),
#       (21, -0.5, 20.5),
#       (100, -60, 60),
#       (2, -0.5, 1.5),
#       (50, -0.5, 100.5),
#       (2, -0.5, 1.5),
#       (41, -0.5, 40.5),
#       (7, -3.5, 3.5),
#   ]
#   threshs = [  #x1,y1,x2,_
#       [0.5],
#       [50],
#       [0.5],
#       [2],
#       [0.5],
#       [0.5],
#       [-5, 12.5],
#       [0.5],
#       [50],
#       [0.5],
#       [0.5],
#       [-0.5, 0.5],
#   ]

#   c = canvas(3, 4)
#   for i, (x1, x2, l, b, ths) in enumerate(zip(vals1, vals2, labs, bins, threshs)):
#     c.cd(i + 1)
#     hh1 = H1D(x1, f';{l};fraction of events', b, norm=True)
#     hh2 = H1D(x2, f';{l};fraction of events', b, norm=True)
#     hh2.SetLineColor(28)
#     ymax = max([hh1.GetMaximum(), hh2.GetMaximum()]) * 1.05
#     ymin = min([hh1.GetMinimum(0),
#                 hh2.GetMinimum(0)])  #min([1 / len(x1), 1 / len(x2)])  #min([hh1.GetMinimum(), hh2.GetMinimum()]) / 10
#     for hh in (hh1, hh2):
#       hh.SetMinimum(0)
#       hh.SetMaximum(ymax)
#       if b[0] > 10 and 'eta' not in l:
#         c.cd(i + 1).SetLogy()
#         hh.SetMinimum(ymin)
#       hh.Draw('histsame')
#       gc.append(hh)

#     for th in ths:
#       line = TLine(th, hh.GetMinimum(), th, ymax)
#       line.SetLineWidth(4)
#       line.SetLineColor(rt.kRed)
#       line.SetLineStyle(rt.kDashed)
#       line.Draw()
#       gc.append(line)

#   c.Print(out_dir + 'cut_vars' + ending)
