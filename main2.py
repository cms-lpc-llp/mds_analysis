"""main.py

Main script for the CMS Run3 analysis
"""

import os
import sys
import pathlib

import numpy as np
from argparse import ArgumentParser

import uproot
import awkward as ak
import awkward.numba

import ROOT as rt
from src import CMS_lumi, tdrstyle
from ROOT import TChain, TCanvas, TH1F, TH2F, TF1, TLatex, TGraph, RDataFrame
from src.histo_utilities import create_prof1D, create_TGraph, create_TH1D, create_TH2D, std_color_list
from src.helper_functions import (getRecoTime, find_nearest, deltaPhi, deltaR, lor, land, lxor, lnot, canvas,
                                  weight_calc)

###########################################


class MuonSystem():

  def __init__(self, file_name, tree_name='MuonSystem') -> None:
    self.gc = []  # Garbage collector for ROOT/C++
    self.fuproot = uproot.open(file_in + ':' + tree_name)
    self.arrs = self.fuproot.arrays()

    self.mh, self.ms, self.ct = 0, 0, 0
    for stat in file_in.split('_'):
      if 'MH' in stat[:2]:
        self.mh = int(stat[3:])
      elif 'MS' in stat[:2]:
        self.ms = int(stat[3:])
      elif 'ctau' in stat[:4]:
        self.ct = int(stat[5:])

  def __getitem__(self, key):
    return self.arrs[key]

  def __setitem__(self, key, value):
    self.arrs[key] = value

  def create_selection(self, vars, func):
    pass

  def parse_selection(self, vars):
    # stack = []

    # sel = np.ones_like()
    # while stack:
    #   v = stack.pop(0)

    #   if '(' in s:
    #     s = s[1:]
    #     pass
    #   elif ')' in s:

    #   if s not in self.arrs:
    #     raise ValueError(f'{s} not found in possible selections')
    pass

  def fix_nclusters(self):
    self.arrs['nCscRechitClusters'] = np.sum(self.arrs['cscRechitClusterSize'] > 0, axis=1)
    self.arrs['nDtRechitClusters'] = np.sum(self.arrs['dtRechitClusterSize'] > 0, axis=1)
    self.arrs['nRechitClusters'] = self.arrs['nCscRechitClusters'] + self.arrs['nDtRechitClusters']

  def Histo1D(self, xname, title=None, bins=None, **kwargs):
    if title is None:
      title = f'{xname};{xname};count'

    if 'x' in kwargs:
      xv = kwargs['x']
    else:
      if '[' in xname:
        xname, sel = xname.split('[')
        sel = sel[:-1]  #remove ']'

        if 'evt' in sel:
          xv = self.arrs[self.arrs[sel]][xname]
        else:
          xv = self.arrs[xname][self.arrs[sel]]
      else:
        xv = self.arrs[xname]

    x = np.asarray(ak.flatten(xv, axis=None))

    # if xname not in self.arrs:
    #   raise ValueError(f'{xname} not found in tree.')

    if isinstance(bins, int) or bins is None:
      xmin, xmax = np.quantile(x, (0.005, 0.995))
      bins = (bins if isinstance(bins, int) else 100, xmin, xmax)

    name = xname + '_' + str(np.random.randint(0, 10000))
    h = rt.TH1F(name, title, *bins)

    if 'c' in kwargs:
      h.SetLineColor(kwargs['c'])  #TODO: implement cases for e.g. 'r' or 'b'
    if 'lw' in kwargs:
      h.SetLineWidth(kwargs['lw'])
    if 'ls' in kwargs:
      h.SetLineStyle(kwargs['ls'])  #TODO: implement cases for '-', '--', and ':'

    for xx in x:
      h.Fill(xx)

    self.gc.append(h)
    return h

  def Histo2D(self, xname, yname, title=None, bins=None, fill='seq', **kwargs):
    """Defaults to sequential filling of the histogram (unless x and y have different shapes).
    Use fill=('seq' | 0) for sequential filling
    Use fill=('comb' | 1) for combinatorial filling"""

    if title is None:
      title = f'{xname}_{yname};{xname};{yname};count'

    if '[' in xname:
      xname, sel = xname.split('[')
      sel = sel[:-1]  #remove ']'

      if 'evt' in sel:
        xv = self.arrs[self.arrs[sel]][xname]
      else:
        xv = self.arrs[xname][self.arrs[sel]]
    else:
      xv = self.arrs[xname]

    if '[' in yname:
      yname, sel = yname.split('[')
      sel = sel[:-1]  #remove ']'

      if 'evt' in sel:
        yv = self.arrs[self.arrs[sel]][yname]
      else:
        yv = self.arrs[yname][self.arrs[sel]]
    else:
      yv = self.arrs[yname]

    x = np.asarray(ak.flatten(xv, axis=None))
    y = np.asarray(ak.flatten(yv, axis=None))

    if x.shape != y.shape:
      raise ValueError(f'x and y don\'t have the same shape. {x.shape} != {y.shape}')

    # if bins is None or len(bins) == 2:
    #   xmin, xmax = np.quantile(x, (0.005, 0.995))
    #   ymin, ymax = np.quantile(y, (0.005, 0.995))
    #   bins = (bins if isinstance(bins, int) else 100, xmin, xmax)

    name = xname + '_' + yname + '_' + str(np.random.randint(0, 10000))
    h = rt.TH2F(name, title, *bins)

    if 'c' in kwargs:
      h.SetLineColor(kwargs['c'])  #TODO: implement cases for e.g. 'r' or 'b'
    if 'lw' in kwargs:
      h.SetLineWidth(kwargs['lw'])
    if 'ls' in kwargs:
      h.SetLineStyle(kwargs['ls'])  #TODO: implement cases for '-', '--', and ':'

    for xx, yy in zip(x, y):
      h.Fill(xx, yy)

    self.gc.append(h)
    return h

  def ctau_scan(self):
    pass


###########################################


def get_lat_leg(leg_coords=(0.6, 0.7, 0.9, 0.85)):
  lat = rt.TLatex()
  lat.SetTextColor(rt.kRed)
  lat.SetTextSize(0.03)
  lat.SetTextAlign(33)

  leg = rt.TLegend(*leg_coords)
  leg.SetTextSize(0.025)
  leg.SetBorderSize(0)
  leg.SetEntrySeparation(0.01)

  return lat, leg


if __name__ == '__main__':

  out_dir = '/home/psimmerl/Documents/CMS/LLP/reports/weekly/feb16/'
  data_dir = '/home/psimmerl/Documents/CMS/LLP/data/raw/'
  MC_file = data_dir + 'ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted.root'
  run3_file = data_dir + 'DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi.root'

  tree_name = 'MuonSystem'
  file_in = run3_file

  dets = ['CSC', 'DT']
  pi = rt.TMath.Pi()

  parser = ArgumentParser(prog='LLPAnalysis', description='LLP skimming for Run3')
  parser.add_argument('-o', '--out', default=out_dir, help='Output directory to write to.')
  parser.add_argument('-f', '--file', default=file_in, help='ROOT file to run the analysis on.')
  parser.add_argument('-t', '--tree', default='MuonSystem', help='TTree to perform the analysis on.')
  parser.add_argument('-s', '--MC', default=True, help='Flag to treat files as data or simulation.')
  parser.add_argument('-j', '--jobs', default=0, help='Number of jobs to run.')
  parser.add_argument('-b', '--batch', default=True, help='Set batch mode (turn drawing off).')
  args = parser.parse_args()

  file_in, tree_name, is_mc = args.file, args.tree, args.MC
  out_dir, jobs = args.out, args.jobs
  met = False

  is_mc = 'DisplacedJet' not in file_in

  if args.batch:
    rt.gROOT.SetBatch()

  a = tdrstyle.setTDRStyle()
  CMS_lumi.writeExtraText = 0
  rt.gStyle.SetOptFit(0)  #1011)

  pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)  # make out directory if it doesn't exist

  print(f'Using output directory \'{out_dir}\'')
  print(f'Processing tree \'{tree_name}\' from \'{file_in}\':')
  print(f'\tTreating file as {"MC" if is_mc else "DATA"}')

  ms = MuonSystem(file_name=file_in, tree_name=tree_name)
  ms.fix_nclusters()

  # if not is_mc:
  print('REMOVING ALL EVENTS WITH N CLUSTERS > 2')
  ms.arrs = ms[ms['nRechitClusters'] <= 2]

  ##########

  ms['gLLP_abs_eta'] = np.abs(ms['gLLP_eta'])
  ms['gLLP_decay_vertex_abs_r'] = np.abs(ms['gLLP_decay_vertex_r'])
  ms['gLLP_decay_vertex_abs_z'] = np.abs(ms['gLLP_decay_vertex_z'])

  ##########

  ms['evt_pass_met200'] = ms['met'] > 200

  ms['csc_gLLP_decay'] = land(ms['gLLP_abs_eta'] < 2.4, ms['gLLP_decay_vertex_abs_r'] < 695.5,
                              ms['gLLP_decay_vertex_z'] < 1100, ms['gLLP_decay_vertex_abs_z'] > 661)
  ms['dt_gLLP_decay'] = land(ms['gLLP_decay_vertex_abs_z'] < 661, ms['gLLP_decay_vertex_abs_r'] > 380,
                             ms['gLLP_decay_vertex_r'] < 738)

  ms['ngLLP_decay_csc'] = np.sum(ms['csc_gLLP_decay'], axis=1)
  ms['ngLLP_decay_dt'] = np.sum(ms['dt_gLLP_decay'], axis=1)

  ms['evt_gLLP_2csc0dt'] = land(ms['ngLLP_decay_csc'] == 2, ms['ngLLP_decay_dt'] == 0)
  ms['evt_gLLP_1csc1dt'] = land(ms['ngLLP_decay_csc'] == 1, ms['ngLLP_decay_dt'] == 1)
  ms['evt_gLLP_0csc2dt'] = land(ms['ngLLP_decay_csc'] == 0, ms['ngLLP_decay_dt'] == 2)

  ms['csc_match'] = ms['cscRechitCluster_match_gLLP']  # just shorthand
  ms['dt_match'] = ms['dtRechitCluster_match_gLLP']  # just shorthand

  # ms['cscRechitClusterTimeWeighted1'] = ms['cscRechitClusterTimeWeighted'][:, 0]
  # ms['cscRechitClusterTimeWeighted2'] = ms['cscRechitClusterTimeWeighted'][:, 1]
  #

  ms['csc_pass_jet_veto'] = lnot(land(ms['cscRechitClusterJetVetoLooseId'], ms['cscRechitClusterJetVetoPt'] > 50))
  ms['dt_pass_jet_veto'] = lnot(land(ms['dtRechitClusterJetVetoLooseId'], ms['dtRechitClusterJetVetoPt'] > 50))

  ms['csc_pass_muon_veto'] = lnot(ms['cscRechitClusterMuonVetoLooseId'])
  ms['dt_pass_muon_veto'] = lnot(ms['dtRechitClusterMuonVetoLooseId'])

  ms['csc_pass_jet_muon_veto'] = land(ms['csc_pass_jet_veto'], ms['csc_pass_muon_veto'])
  ms['dt_pass_jet_muon_veto'] = land(ms['dt_pass_jet_veto'], ms['dt_pass_muon_veto'])

  ms['csc_fail_jet_muon_veto'] = lnot(ms['csc_pass_jet_muon_veto'])
  ms['dt_fail_jet_muon_veto'] = lnot(ms['dt_pass_jet_muon_veto'])

  ms['evt_2csc_cscL1'] = land(ms['nCscRechitClusters'] == 2, ms['HLTDecision'][:, 569])
  ms['csc_pass_jet_muon_veto_2csc_cscL1'] = land(ms['csc_pass_jet_muon_veto'], ms['evt_2csc_cscL1'])
  ms['csc_fail_jet_muon_veto_2csc_cscL1'] = land(ms['csc_fail_jet_muon_veto'], ms['evt_2csc_cscL1'])

  ms['evt_cscdt_cscL1'] = land(ms['nCscRechitClusters'] == 1, ms['nDtRechitClusters'] == 1, ms['HLTDecision'][:, 569])
  ms['csc_pass_jet_muon_veto_cscdt_cscL1'] = land(ms['csc_pass_jet_muon_veto'], ms['evt_cscdt_cscL1'])
  ms['dt_pass_jet_muon_veto_cscdt_cscL1'] = land(ms['dt_pass_jet_muon_veto'], ms['evt_cscdt_cscL1'])

  ms['csc_fail_jet_muon_veto_cscdt_cscL1'] = land(ms['csc_fail_jet_muon_veto'], ms['evt_cscdt_cscL1'])
  ms['dt_fail_jet_muon_veto_cscdt_cscL1'] = land(ms['dt_fail_jet_muon_veto'], ms['evt_cscdt_cscL1'])

  ms['csc_oot'] = ms['cscRechitClusterTimeWeighted'] < -12.5

  ms['csc_sig'] = land(ms['csc_pass_jet_muon_veto'], lnot(ms['csc_oot']))
  ms['dt_sig'] = ms['dt_pass_jet_muon_veto']  #lor(ms['dt_pass_jet_muon_veto'], lnot(ms['dt_oot']))

  ms['csc_bkg'] = lor(ms['csc_fail_jet_muon_veto'], ms['csc_oot'])
  ms['dt_bkg'] = ms['dt_fail_jet_muon_veto']  #lor(ms['dt_fail_jet_muon_veto'], ms['dt_oot'])

  ms['nCscRechitClusters_sig'] = ak.sum(ms['csc_sig'], axis=1)
  ms['nDtRechitClusters_sig'] = ak.sum(ms['dt_sig'], axis=1)
  ms['nRechitClusters_sig'] = ms['nCscRechitClusters_sig'] + ms['nDtRechitClusters_sig']

  ms['nCscRechitClusters_bkg'] = ak.sum(ms['csc_bkg'], axis=1)
  ms['nDtRechitClusters_bkg'] = ak.sum(ms['dt_bkg'], axis=1)
  ms['nRechitClusters_bkg'] = ms['nCscRechitClusters_bkg'] + ms['nDtRechitClusters_bkg']

  ms['csc_limited'] = land(50 <= ms['cscRechitClusterSize'], ms['cscRechitClusterSize'] < 100)
  ms['dt_limited'] = land(50 <= ms['dtRechitClusterSize'], ms['dtRechitClusterSize'] < 60)

  print('\nGenerating Histograms')
  ####################

  c = canvas()

  # cut_comment = 'nClusters #leq 2, #splitline{CSC: (Jet & P_{T,Jet}>50GeV) | Muon | (t<-12.5ns)}{DT: (Jet & P_{T,Jet}>50GeV) | Muon}'
  cut_comment = 'CSC+DT#leq2, CSC: (Jet & P_{T,Jet}>50GeV) | Muon | (t<-12.5ns), DT: (Jet & P_{T,Jet}>50GeV) | Muon'
  file_comment = '_bkg' + ('_MC' if is_mc else '') + '.png'

  lat, leg = get_lat_leg(leg_coords=(0.4, 0.7, 0.5, 0.85))

  hcsc_all = ms.Histo1D('cscRechitClusterEta', c=std_color_list[0])
  hcsc_sig = ms.Histo1D('cscRechitClusterEta[csc_sig]', c=std_color_list[1])
  hcsc_bkg = ms.Histo1D('cscRechitClusterEta[csc_bkg]', c=std_color_list[2])
  hdt_all = ms.Histo1D('dtRechitClusterEta', c=std_color_list[3])
  hdt_sig = ms.Histo1D('dtRechitClusterEta[dt_sig]', c=std_color_list[4])
  hdt_bkg = ms.Histo1D('dtRechitClusterEta[dt_bkg]', c=std_color_list[5])

  if is_mc:
    hhs = [hcsc_all, hcsc_sig, hcsc_bkg, hdt_all, hdt_sig, hdt_bkg]
    lls = [det + ' ' + reg for det in dets for reg in ['all', 'sig', 'bkg']]
  else:
    hhs = [hcsc_bkg, hdt_bkg]
    lls = [det + reg for det in dets for reg in ['bkg']]

  ymax = max([hh.GetMaximum() for hh in hhs])
  for hh, ll in zip(hhs, lls):
    hh.SetMaximum(ymax * 1.05)
    leg.AddEntry(hh, ll, 'L')
    hh.Draw('same')

  leg.Draw()
  lat.DrawLatexNDC(1, 1, cut_comment)
  c.Print(out_dir + 'eta' + file_comment)

  ####################

  rt.gStyle.SetPaintTextFormat('.0f')  #.2G')
  lat, leg = get_lat_leg(leg_coords=(0.4, 0.7, 0.5, 0.85))

  bins = (7, -0.5, 6.5, 7, -0.5, 6.5)
  hcsc_all = ms.Histo2D('nCscRechitClusters', 'nDtRechitClusters', bins=bins, c=std_color_list[0])
  hcsc_sig = ms.Histo2D('nCscRechitClusters_sig', 'nDtRechitClusters_sig', bins=bins, c=std_color_list[1])
  hcsc_bkg = ms.Histo2D('nCscRechitClusters_bkg', 'nDtRechitClusters_bkg', bins=bins, c=std_color_list[2])
  hdt_all = ms.Histo2D('nCscRechitClusters', 'nDtRechitClusters', bins=bins, c=std_color_list[3])
  hdt_sig = ms.Histo2D('nCscRechitClusters_sig', 'nDtRechitClusters_sig', bins=bins, c=std_color_list[4])
  hdt_bkg = ms.Histo2D('nCscRechitClusters_bkg', 'nDtRechitClusters_bkg', bins=bins, c=std_color_list[5])

  if is_mc:
    c = canvas(2, 3)
    hhs = [hcsc_all, hcsc_sig, hcsc_bkg, hdt_all, hdt_sig, hdt_bkg]
    lls = [det + ' ' + reg for det in dets for reg in ['all', 'sig', 'bkg']]
  else:
    c = canvas(1, 2)
    hhs = [hcsc_bkg, hdt_bkg]
    lls = [det + reg for det in dets for reg in ['bkg']]

  ymax = max([hh.GetMaximum() for hh in hhs])
  for i, (hh, ll) in enumerate(zip(hhs, lls)):
    c.cd(i + 1).SetLogz()
    # hh.SetMaximum(ymax)
    # hh.SetMaximum(1)
    hh.SetMarkerSize(1.3)
    hh.Draw('coltext')

    lat.DrawLatexNDC(1, 1, cut_comment)
  c.Print(out_dir + 'ncscndt' + file_comment)

  ####################

  # timing
  # All clusters

  cut_comment = 'CSC+DT#leq2, CSC: (Jet & P_{T,Jet}>50GeV) | Muon, DT: (Jet & P_{T,Jet}>50GeV) | Muon'
  file_comment = '_fail_jet_muon' + ('_MC' if is_mc else '') + '.png'
  for bins, stat in zip([(201, -80, 80), (21, -10, 10)], ['full', 'zoom']):
    c = canvas()
    c.SetLogy()

    lat, leg = get_lat_leg(leg_coords=(0.4, 0.7, 0.5, 0.85))

    hcsc_all = ms.Histo1D('cscRechitClusterTimeWeighted', bins=bins, c=std_color_list[0])
    hcsc_sig = ms.Histo1D('cscRechitClusterTimeWeighted[csc_pass_jet_muon_veto]', bins=bins, c=std_color_list[1])
    hcsc_bkg = ms.Histo1D('cscRechitClusterTimeWeighted[csc_fail_jet_muon_veto]', bins=bins, c=std_color_list[2])
    hdt_all = ms.Histo1D('dtRechitCluster_match_RPCBx_dPhi0p5', bins=bins, c=std_color_list[3])
    hdt_sig = ms.Histo1D('dtRechitCluster_match_RPCBx_dPhi0p5[dt_pass_jet_muon_veto]', bins=bins, c=std_color_list[4])
    hdt_bkg = ms.Histo1D('dtRechitCluster_match_RPCBx_dPhi0p5[dt_fail_jet_muon_veto]', bins=bins, c=std_color_list[5])

    if is_mc:
      hhs = [hcsc_all, hcsc_sig, hcsc_bkg, hdt_all, hdt_sig, hdt_bkg]
      lls = [det + ' ' + reg for det in dets for reg in ['all', 'sig', 'bkg']]
    else:
      hhs = [hcsc_bkg, hdt_bkg]
      lls = [det + reg for det in dets for reg in ['bkg']]

    ymax = max([hh.GetMaximum() for hh in hhs])
    for hh, ll in zip(hhs, lls):
      hh.SetMaximum(ymax * 1.05)
      hh.SetMinimum(0.1)
      leg.AddEntry(hh, ll, 'L')
      hh.Draw('same')

    leg.Draw()
    lat.DrawLatexNDC(1, 1, cut_comment)
    c.Print(out_dir + 'time_' + stat + file_comment)

  # 2 CSC (CSC passes HLT)

  # cut_comment = 'CSC+CSC, CSC: (Jet & P_{T,Jet}>50GeV) | Muon & HLT'
  # file_comment = '_fail_jet_muon_2csc_cscL1' + ('_MC' if is_mc else '') + '.png'
  # for bins, stat in zip([(201, -80, 80), (21, -10, 10)], ['full', 'zoom']):
  #   c = canvas()
  #   c.SetLogy()

  #   lat, leg = get_lat_leg(leg_coords=(0.4, 0.7, 0.5, 0.85))

  #   hcsc1_all = ms.Histo1D('cscRechitClusterTimeWeighted1[evt_2csc_cscL1]', bins=bins, c=std_color_list[0])
  #   hcsc1_sig = ms.Histo1D('cscRechitClusterTimeWeighted1[csc_pass_jet_muon_veto_2csc_cscL1]',
  #                          bins=bins,
  #                          c=std_color_list[1])
  #   hcsc1_bkg = ms.Histo1D('cscRechitClusterTimeWeighted1[csc_fail_jet_muon_veto_2csc_cscL1]',
  #                          bins=bins,
  #                          c=std_color_list[2])

  #   hcsc2_all = ms.Histo1D('cscRechitClusterTimeWeighted2[evt_2csc_cscL1]', bins=bins, c=std_color_list[0])
  #   hcsc2_sig = ms.Histo1D('cscRechitClusterTimeWeighted2[csc_pass_jet_muon_veto_2csc_cscL1]',
  #                          bins=bins,
  #                          c=std_color_list[1])
  #   hcsc2_bkg = ms.Histo1D('cscRechitClusterTimeWeighted2[csc_fail_jet_muon_veto_2csc_cscL1]',
  #                          bins=bins,
  #                          c=std_color_list[2])

  #   if is_mc:
  #     hhs = [hcsc_all, hcsc_sig, hcsc_bkg, hdt_all, hdt_sig, hdt_bkg]
  #     lls = [det + ' ' + reg for det in ['Big CSC', 'Small CSC'] for reg in ['all', 'sig', 'bkg']]
  #   else:
  #     hhs = [hcsc_bkg, hdt_bkg]
  #     lls = [det + reg for det in ['Big CSC', 'Small CSC'] for reg in ['bkg']]

  #   ymax = max([hh.GetMaximum() for hh in hhs])
  #   for hh, ll in zip(hhs, lls):
  #     hh.SetMaximum(ymax * 1.05)
  #     hh.SetMinimum(0.1)
  #     leg.AddEntry(hh, ll, 'L')
  #     hh.Draw('same')

  #   leg.Draw()
  #   lat.DrawLatexNDC(1, 1, cut_comment)
  #   c.Print(out_dir + 'time_' + stat + file_comment)

  # 1 CSC - 1 DT (CSC passes HLT)

  cut_comment = 'CSC+DT, CSC: (Jet & P_{T,Jet}>50GeV) | Muon & HLT, DT: (Jet & P_{T,Jet}>50GeV) | Muon'
  file_comment = '_fail_jet_muon_cscdt_cscL1' + ('_MC' if is_mc else '') + '.png'
  for bins, stat in zip([(201, -80, 80), (21, -10, 10)], ['full', 'zoom']):
    c = canvas()
    c.SetLogy()

    lat, leg = get_lat_leg(leg_coords=(0.4, 0.7, 0.5, 0.85))

    hcsc_all = ms.Histo1D('cscRechitClusterTimeWeighted[evt_cscdt_cscL1]', bins=bins, c=std_color_list[0])
    hcsc_sig = ms.Histo1D('cscRechitClusterTimeWeighted[csc_pass_jet_muon_veto_cscdt_cscL1]',
                          bins=bins,
                          c=std_color_list[1])
    hcsc_bkg = ms.Histo1D('cscRechitClusterTimeWeighted[csc_fail_jet_muon_veto_cscdt_cscL1]',
                          bins=bins,
                          c=std_color_list[2])
    hdt_all = ms.Histo1D('dtRechitCluster_match_RPCBx_dPhi0p5[evt_cscdt_cscL1]', bins=bins, c=std_color_list[3])
    hdt_sig = ms.Histo1D('dtRechitCluster_match_RPCBx_dPhi0p5[dt_pass_jet_muon_veto_cscdt_cscL1]',
                         bins=bins,
                         c=std_color_list[4])
    hdt_bkg = ms.Histo1D('dtRechitCluster_match_RPCBx_dPhi0p5[dt_fail_jet_muon_veto_cscdt_cscL1]',
                         bins=bins,
                         c=std_color_list[5])

    if is_mc:
      hhs = [hcsc_all, hcsc_sig, hcsc_bkg, hdt_all, hdt_sig, hdt_bkg]
      lls = [det + ' ' + reg for det in dets for reg in ['all', 'sig', 'bkg']]
    else:
      hhs = [hcsc_bkg, hdt_bkg]
      lls = [det + reg for det in dets for reg in ['bkg']]

    ymax = max([hh.GetMaximum() for hh in hhs])
    for hh, ll in zip(hhs, lls):
      hh.SetMaximum(ymax * 1.05)
      hh.SetMinimum(0.1)
      leg.AddEntry(hh, ll, 'L')
      hh.Draw('same')

    leg.Draw()
    lat.DrawLatexNDC(1, 1, cut_comment)
    c.Print(out_dir + 'time_' + stat + file_comment)

  ####################

  # Cluster size with exponential fits
  c = canvas()

  cut_comment = 'CSC+DT#leq2, CSC: 50 #leq Size < 100,DT: 50 #leq Size < 60'
  file_comment = '_limited' + ('_MC' if is_mc else '') + '.png'

  lat, leg = get_lat_leg(leg_coords=(0.4, 0.7, 0.5, 0.85))

  bins = (50, 50, 100)
  hcsc_all = ms.Histo1D('cscRechitClusterSize[csc_limited]', bins=bins, c=std_color_list[0])
  hdt_all = ms.Histo1D('dtRechitClusterSize[dt_limited]', bins=bins, c=std_color_list[1])

  leg.AddEntry(hcsc_all, 'CSC', 'L')
  leg.AddEntry(hdt_all, 'DT', 'L')

  ymax = max([h.GetMaximum() for h in (hcsc_all, hdt_all)])
  for h in (hcsc_all, hdt_all):
    h.SetMaximum(ymax * 1.05)
    h.Draw('same')

  npars = 2
  f1 = TF1('exp_fit_csc_limited', 'expo', 50, 100)
  f2 = TF1('exp_fit_dt_limited', 'expo', 50, 60)

  f1.SetLineColor(hcsc_all.GetLineColor())
  f2.SetLineColor(hdt_all.GetLineColor())
  f1.SetLineWidth(3)
  f2.SetLineWidth(3)

  hcsc_all.Fit(f1, 'R')
  hdt_all.Fit(f2, 'R')

  leg.Draw()
  lat.DrawLatexNDC(1, 1, cut_comment)

  pars1, errs1 = [f1.GetParameter(p) for p in range(2)], [f1.GetParError(p) for p in range(2)]
  pars2, errs2 = [f2.GetParameter(p) for p in range(2)], [f2.GetParError(p) for p in range(2)]
  lat.SetTextColor(rt.kMagenta)
  lat.DrawLatexNDC(0.95, 0.90, f'f(x) = exp[P0 + P1*x]')
  lat.SetTextColor(f1.GetLineColor())
  lat.DrawLatexNDC(0.95, 0.85, f'P0 = {pars1[0]:.4f}#pm{errs1[0]:.4f}')
  lat.DrawLatexNDC(0.95, 0.80, f'P1 = {pars1[1]:.4f}#pm{errs1[1]:.4f}')
  lat.SetTextColor(f2.GetLineColor())
  lat.DrawLatexNDC(0.95, 0.75, f'P0 = {pars2[0]:.4f}#pm{errs2[0]:.4f}')
  lat.DrawLatexNDC(0.95, 0.70, f'P1 = {pars2[1]:.4f}#pm{errs2[1]:.4f}')

  c.Print(out_dir + 'size' + file_comment)

  ###############
  # Cluster size with exponential fits CSC-DT
  # c.SetLogy()

  cut_comment = 'CSC: Size < 100, DT: 50 #leq Size < 60, Require 1 CSC + 1 DT'
  file_comment = '_limited' + ('_MC' if is_mc else '') + '.png'

  c = canvas()
  lat, leg = get_lat_leg(leg_coords=(0.2, 0.2, 0.35, 0.5))

  sel_1csc1dt = land(ms['nCscRechitClusters'] == 1, ms['nDtRechitClusters'] == 1)

  if is_mc:
    sel_1csc1dt = land(ak.sum(ms['csc_match'], axis=1) == 1, ak.sum(ms['dt_match'], axis=1) == 1)
    cut_comment = 'Signal, ' + cut_comment

  csc_size = ms['cscRechitClusterSize']
  dt_size = ms['dtRechitClusterSize']

  csc_pass = csc_size > 100
  dt_pass = land(50 <= dt_size, dt_size < 60)

  csc_pass_jet = land(csc_pass, ms['csc_pass_jet_veto'])
  dt_pass_jet = land(dt_pass, ms['dt_pass_jet_veto'])

  csc_pass_muon = land(csc_pass, ms['csc_pass_muon_veto'])
  dt_pass_muon = land(dt_pass, ms['dt_pass_muon_veto'])

  csc_pass_jet_muon = land(csc_pass, ms['csc_pass_jet_muon_veto'])
  dt_pass_jet_muon = land(dt_pass, ms['dt_pass_jet_muon_veto'])

  evt_pass = land(sel_1csc1dt, ak.sum(csc_pass, axis=1), ak.sum(dt_pass, axis=1))

  evt_pass_csc_jet = land(sel_1csc1dt, ak.sum(csc_pass_jet, axis=1), ak.sum(dt_pass, axis=1))
  evt_pass_csc_muon = land(sel_1csc1dt, ak.sum(csc_pass_muon, axis=1), ak.sum(dt_pass, axis=1))
  evt_pass_csc_jet_muon = land(sel_1csc1dt, ak.sum(csc_pass_jet_muon, axis=1), ak.sum(dt_pass, axis=1))
  evt_pass_dt_jet = land(sel_1csc1dt, ak.sum(csc_pass, axis=1), ak.sum(dt_pass_jet, axis=1))
  evt_pass_dt_muon = land(sel_1csc1dt, ak.sum(csc_pass, axis=1), ak.sum(dt_pass_muon, axis=1))
  evt_pass_dt_jet_muon = land(sel_1csc1dt, ak.sum(csc_pass, axis=1), ak.sum(dt_pass_jet_muon, axis=1))
  evt_pass_csc_jet_dt_jet = land(sel_1csc1dt, ak.sum(csc_pass_jet, axis=1), ak.sum(dt_pass_jet, axis=1))
  evt_pass_csc_muon_dt_muon = land(sel_1csc1dt, ak.sum(csc_pass_muon, axis=1), ak.sum(dt_pass_muon, axis=1))
  evt_pass_csc_jet_muon_dt_jet_muon = land(sel_1csc1dt, ak.sum(csc_pass_jet_muon, axis=1),
                                           ak.sum(dt_pass_jet_muon, axis=1))

  sel_evts = [
      evt_pass,
      evt_pass_csc_jet,
      evt_pass_csc_muon,
      evt_pass_csc_jet_muon,
      evt_pass_dt_jet,
      evt_pass_dt_muon,
      evt_pass_dt_jet_muon,
      evt_pass_csc_jet_dt_jet,
      evt_pass_csc_muon_dt_muon,
      evt_pass_csc_jet_muon_dt_jet_muon,
  ]
  labels = [
      '1 CSC + 1 DT',
      'After CSC Jet Veto',
      'After CSC Muon Veto',
      'After CSC Jet & Muon Vetoes',
      'After DT Jet Veto',
      'After DT Muon Veto',
      'After DT Jet & Muon Vetoes',
      'After CSC/DT Jet Veto',
      'After CSC/DT Muon Veto',
      'After CSC/DT Jet & Muon Vetoes',
  ]

  a, b = 100, 250
  print(f'Integrating fit from {a} to {b}')
  hfs = []
  bins = (10, 50, 60)
  for i, (sel_evt, ll) in enumerate(zip(sel_evts, labels)):
    hh = ms.Histo1D('', title=';DT Cluster Size;Count', bins=bins, c=std_color_list[i], x=dt_size[sel_evt])
    leg.AddEntry(hh, ll, 'L')

    ff = TF1('exp_fit_dt_limited', 'expo', 50, 60)
    # ff.SetTitle(hh.GetTitle())
    # ff.GetHistogram().GetXaxis().SetTitle('DT Cluster Size')
    # ff.GetHistogram().GetYaxis().SetTitle('Count')
    ff.SetLineColor(hh.GetLineColor())
    ff.SetLineStyle(hh.GetLineStyle())
    ff.SetLineWidth(3)

    hh.Fit(ff, 'QR')
    # leg.AddEntry(ff, ll, 'L')

    intgrl, interr = ff.Integral(a, b), ff.IntegralError(a, b)
    print(f'{ll} \t {intgrl:.0f} \pm {interr:.0f}')

    hfs.append((hh, ff))

  ymax = max([hh.GetMaximum() for hh, ff in hfs])
  for i, (hh, ff) in enumerate(hfs):
    hh.SetMinimum(0)  #ymax * 1.05)
    hh.SetMaximum(ymax * 1.05)
    hh.Draw('same')
  leg.Draw()
  lat.DrawLatexNDC(1, 1, cut_comment)
  c.Print(out_dir + 'size_1csc1dt_multi_with_hists' + file_comment)

  c = canvas()
  for i, (hh, ff) in enumerate(hfs):
    ff.SetMinimum(0)  #ymax * 1.05)
    ff.SetMaximum(ymax * 1.05)
    ff.Draw('same' if i else '')  #'Lsame')
  leg.Draw()
  lat.DrawLatexNDC(1, 1, cut_comment)
  c.Print(out_dir + 'size_1csc1dt_multi_wout_hists' + file_comment)

  ###############
  # Cluster size with exponential fits CSC-CSC
  # c = canvas()
  # # c.SetLogy()

  # cut_comment = 'CSC: 50 #leq Size < 100, DT: 50 #leq Size < 60, 1 CSC or 1 DT After Size Cut'
  # file_comment = '_limited' + ('_MC' if is_mc else '') + '.png'

  # sel_1csc1dt = land(ms['nCscRechitClusters'] == 1, ms['nDtRechitClusters'] == 1)
  # csc_size = ms['cscRechitClusterSize']
  # dt_size = ms['dtRechitClusterSize']

  # csc_pass = csc_size > 100
  # dt_pass = land(50 <= dt_size, dt_size < 100)

  # csc_pass_jet = land(csc_size > 100, ms['csc_pass_jet_veto'])
  # dt_pass_jet = land(50 <= dt_size, dt_size < 100, ms['dt_pass_jet_veto'])

  # csc_pass_muon = land(csc_size > 100, ms['csc_pass_muon_veto'])
  # dt_pass_muon = land(50 <= dt_size, dt_size < 100, ms['dt_pass_muon_veto'])

  # csc_pass_jet_muon = land(csc_size > 100, ms['csc_pass_jet_muon_veto'])
  # dt_pass_jet_muon = land(50 <= dt_size, dt_size < 100, ms['dt_pass_jet_muon_veto'])

  # evt_pass = land(sel_1csc1dt, ak.sum(csc_pass, axis=1), ak.sum(dt_pass, axis=1))

  # evt_pass_csc_jet = land(sel_1csc1dt, ak.sum(csc_pass_jet, axis=1), ak.sum(dt_pass, axis=1))
  # evt_pass_csc_muon = land(sel_1csc1dt, ak.sum(csc_pass_muon, axis=1), ak.sum(dt_pass, axis=1))
  # evt_pass_csc_jet_muon = land(sel_1csc1dt, ak.sum(csc_pass_jet_muon, axis=1), ak.sum(dt_pass, axis=1))
  # evt_pass_dt_jet = land(sel_1csc1dt, ak.sum(csc_pass, axis=1), ak.sum(dt_pass_jet, axis=1))
  # evt_pass_dt_muon = land(sel_1csc1dt, ak.sum(csc_pass, axis=1), ak.sum(dt_pass_muon, axis=1))
  # evt_pass_dt_jet_muon = land(sel_1csc1dt, ak.sum(csc_pass, axis=1), ak.sum(dt_pass_jet_muon, axis=1))
  # evt_pass_csc_jet_dt_jet = land(sel_1csc1dt, ak.sum(csc_pass_jet, axis=1), ak.sum(dt_pass_jet, axis=1))
  # evt_pass_csc_muon_dt_muon = land(sel_1csc1dt, ak.sum(csc_pass_muon, axis=1), ak.sum(dt_pass_muon, axis=1))
  # evt_pass_csc_jet_muon_dt_jet_muon = land(sel_1csc1dt, ak.sum(csc_pass_jet_muon, axis=1),
  #                                          ak.sum(dt_pass_jet_muon, axis=1))

  # lat, leg = get_lat_leg(leg_coords=(0.50, 0.65, 0.9, 0.90))

  # sel_evts = [
  #     evt_pass,
  #     evt_pass_csc_jet,
  #     evt_pass_csc_muon,
  #     evt_pass_csc_jet_muon,
  #     evt_pass_dt_jet,
  #     evt_pass_dt_muon,
  #     evt_pass_dt_jet_muon,
  #     evt_pass_csc_jet_dt_jet,
  #     evt_pass_csc_muon_dt_muon,
  #     evt_pass_csc_jet_muon_dt_jet_muon,
  # ]
  # labels = [
  #     '1 CSC + 1 DT',
  #     'After CSC Jet Veto',
  #     'After CSC Muon Veto',
  #     'After CSC Jet & Muon Vetoes',
  #     'After DT Jet Veto',
  #     'After DT Muon Veto',
  #     'After DT Jet & Muon Vetoes',
  #     'After CSC Jet Veto, DT Jet Veto',
  #     'After CSC Muon Veto, DT Muon Veto',
  #     'After CSC Jet & Muon Vetoes, DT Jet & Muon Vetoes',
  # ]

  # hcscs, hdts = [], []
  # bins = (10, 50, 60)
  # for i, (sel_evt, ll) in enumerate(zip(sel_evts, labels)):
  #   # hcsc = ms.Histo1D('', title=';Cluster Size;Count', bins=bins, c=std_color_list[i], x=csc_size[sel_csc][sel_evt])
  #   hdt = ms.Histo1D('', title=';Cluster Size;Count', bins=bins, c=std_color_list[i], x=dt_size[sel_evt])

  #   # hdt.SetLineStyle(rt.kDashed)
  #   # leg.AddEntry(hcsc, 'CSC' + ll, 'L')
  #   leg.AddEntry(hdt, ll, 'L')

  #   # f1 = TF1('exp_fit_csc_limited', 'expo', 50, 100)
  #   f2 = TF1('exp_fit_dt_limited', 'expo', 50, 60)
  #   # f1.SetLineColor(hcsc.GetLineColor())
  #   f2.SetLineColor(hdt.GetLineColor())
  #   # f1.SetLineStyle(hcsc.GetLineStyle())
  #   f2.SetLineStyle(hdt.GetLineStyle())
  #   # f1.SetLineWidth(3)
  #   f2.SetLineWidth(3)

  #   # hcsc.Fit(f1, 'R')
  #   hdt.Fit(f2, 'R')

  #   # hcscs.append(hcsc)
  #   hdts.append(hdt)
  #   # hhs.extend([hcsc, hdt])
  #   # hcsc.Draw('same')
  #   # hdt.Draw('same')

  # ymax = max([hh.GetMaximum() for hh in hcscs + hdts])
  # for hh in hcscs + hdts:
  #   # hh.SetMinimum(1000)  #ymax * 1.05)
  #   hh.SetMaximum(ymax * 1.05)
  #   hh.Draw('same')

  # leg.Draw()
  # lat.DrawLatexNDC(1, 1, cut_comment)

  # c.Print(out_dir + 'size_1csc1dt_multi' + file_comment)

  ####################

  # ABCD
  # 1D plots
  c = canvas(1, 2)

  cut_comment = 'CSC+DT#leq2, CSC: (Jet & P_{T,Jet}>50GeV) | Muon | (t<-12.5ns), DT: (Jet & P_{T,Jet}>50GeV) | Muon'
  file_comment = '_bkg' + ('_MC' if is_mc else '') + '.png'

  c.cd(1).SetLogy()
  lat1, leg1 = get_lat_leg(leg_coords=(0.7, 0.7, 0.8, 0.95))
  bins = (100, 0, 3000)
  hcsc_size_all = ms.Histo1D('cscRechitClusterSize', bins=bins, c=std_color_list[0])
  hcsc_size_sig = ms.Histo1D('cscRechitClusterSize[csc_sig]', bins=bins, c=std_color_list[1])
  hcsc_size_bkg = ms.Histo1D('cscRechitClusterSize[csc_bkg]', bins=bins, c=std_color_list[2])
  hdt_size_all = ms.Histo1D('dtRechitClusterSize', bins=bins, c=std_color_list[3])
  hdt_size_sig = ms.Histo1D('dtRechitClusterSize[dt_sig]', bins=bins, c=std_color_list[4])
  hdt_size_bkg = ms.Histo1D('dtRechitClusterSize[dt_bkg]', bins=bins, c=std_color_list[5])
  hsize_all = ms.Histo1D('cscRechitClusterSize', bins=bins, c=std_color_list[6])
  hsize_sig = ms.Histo1D('cscRechitClusterSize[csc_sig]', bins=bins, c=std_color_list[7])
  hsize_bkg = ms.Histo1D('cscRechitClusterSize[csc_bkg]', bins=bins, c=std_color_list[8])
  hsize_all.Add(hdt_size_all)
  hsize_sig.Add(hdt_size_sig)
  hsize_bkg.Add(hdt_size_bkg)

  if is_mc:
    hhs = [
        hcsc_size_all, hcsc_size_sig, hcsc_size_bkg, hdt_size_all, hdt_size_sig, hdt_size_bkg, hsize_all, hsize_sig,
        hsize_bkg
    ]
    lls = [det + ' ' + reg for det in dets + ['both'] for reg in ['all', 'sig', 'bkg']]
  else:
    hhs = [hcsc_size_bkg, hdt_size_bkg, hsize_bkg]
    lls = [det + ' ' + reg for det in dets + ['both'] for reg in ['bkg']]

  ymax = max([hh.GetMaximum() for hh in hhs])
  for hh, ll in zip(hhs, lls):
    hh.SetMaximum(ymax * 1.05)
    hh.SetMinimum(0.1)
    leg1.AddEntry(hh, ll, 'L')
    hh.Draw('same')

  leg1.Draw()
  lat1.DrawLatexNDC(1, 1, cut_comment)

  c.cd(2)
  lat2, leg2 = get_lat_leg(leg_coords=(0.8, 0.15, 0.9, 0.45))
  bins = (100, 0, pi)
  hcsc_dphi_all = ms.Histo1D('cscRechitClusterMet_dPhi', bins=bins, c=std_color_list[0])
  hcsc_dphi_sig = ms.Histo1D('cscRechitClusterMet_dPhi[csc_sig]', bins=bins, c=std_color_list[1])
  hcsc_dphi_bkg = ms.Histo1D('cscRechitClusterMet_dPhi[csc_bkg]', bins=bins, c=std_color_list[2])
  hdt_dphi_all = ms.Histo1D('dtRechitClusterMet_dPhi', bins=bins, c=std_color_list[3])
  hdt_dphi_sig = ms.Histo1D('dtRechitClusterMet_dPhi[dt_sig]', bins=bins, c=std_color_list[4])
  hdt_dphi_bkg = ms.Histo1D('dtRechitClusterMet_dPhi[dt_bkg]', bins=bins, c=std_color_list[5])
  hdphi_all = ms.Histo1D('cscRechitClusterMet_dPhi', bins=bins, c=std_color_list[6])
  hdphi_sig = ms.Histo1D('cscRechitClusterMet_dPhi[csc_sig]', bins=bins, c=std_color_list[7])
  hdphi_bkg = ms.Histo1D('cscRechitClusterMet_dPhi[csc_bkg]', bins=bins, c=std_color_list[8])
  hdphi_all.Add(hdt_dphi_all)
  hdphi_sig.Add(hdt_dphi_sig)
  hdphi_bkg.Add(hdt_dphi_bkg)

  if is_mc:
    hhs = [
        hcsc_dphi_all, hcsc_dphi_sig, hcsc_dphi_bkg, hdt_dphi_all, hdt_dphi_sig, hdt_dphi_bkg, hdphi_all, hdphi_sig,
        hdphi_bkg
    ]
    lls = [det + ' ' + reg for det in dets + ['both'] for reg in ['all', 'sig', 'bkg']]
  else:
    hhs = [hcsc_dphi_bkg, hdt_dphi_bkg, hdphi_bkg]
    lls = [det + ' ' + reg for det in dets + ['both'] for reg in ['bkg']]

  ymax = max([hh.GetMaximum() for hh in hhs])
  for hh, ll in zip(hhs, lls):
    hh.SetMaximum(ymax * 1.05)
    hh.SetMinimum(0)  #.1)
    leg2.AddEntry(hh, ll, 'L')
    hh.Draw('same')

  leg2.Draw()
  lat2.DrawLatexNDC(1, 1, cut_comment)
  c.Print(out_dir + 'ABCD_1D' + file_comment)

  ####################

  if is_mc:
    # MET requiring decay in muon system
    c = canvas()
    bins = (200, 0, 300)
    hmet_gLLP_2csc0dt = ms.Histo1D('met[evt_gLLP_2csc0dt]', bins=bins, c=std_color_list[0])
    hmet_gLLP_1csc1dt = ms.Histo1D('met[evt_gLLP_1csc1dt]', bins=bins, c=std_color_list[1])
    hmet_gLLP_0csc2dt = ms.Histo1D('met[evt_gLLP_0csc2dt]', bins=bins, c=std_color_list[2])
    hmet_gLLP_2csc_or_1csc1dt = ms.Histo1D('met[evt_gLLP_2csc0dt]', bins=bins, c=std_color_list[4])
    hmet_gLLP_2csc_or_2dt_or_1csc1dt = ms.Histo1D('met[evt_gLLP_2csc0dt]', bins=bins, c=std_color_list[3])

    hmet_gLLP_2csc_or_1csc1dt.Add(hmet_gLLP_1csc1dt)
    hmet_gLLP_2csc_or_2dt_or_1csc1dt.Add(hmet_gLLP_1csc1dt)
    hmet_gLLP_2csc_or_2dt_or_1csc1dt.Add(hmet_gLLP_0csc2dt)

    cut_comment = ''
    file_comment = ('_MC' if is_mc else '') + '.png'

    lat, leg = get_lat_leg(leg_coords=(0.4, 0.7, 0.5, 0.85))

    hhs = [
        hmet_gLLP_2csc0dt, hmet_gLLP_1csc1dt, hmet_gLLP_0csc2dt, hmet_gLLP_2csc_or_1csc1dt,
        hmet_gLLP_2csc_or_2dt_or_1csc1dt
    ]
    lls = ['2 CSC', '1 CSC 1 DT', '2 DT', '2 CSC or 1 CSC 1DT', 'Any']

    ymax = max([hh.GetMaximum() for hh in hhs])
    for hh, ll in zip(hhs, lls):
      hh.SetMaximum(ymax * 1.05)
      leg.AddEntry(hh, ll, 'L')
      hh.Draw('same')

    leg.Draw()
    lat.DrawLatexNDC(1, 1, cut_comment)
    c.Print(out_dir + 'met_2gLLP_in_det' + file_comment)

    ####################

    # ctau
    # remember to remove DT-DT clusters

    htots = []
    for met in [False, True]:
      gc = []
      c = canvas()
      plot_ctaus = ['10', '100', '200', '500', '700', '1000', '2000', '5000', '10000', '50000', '100000', '1000000']
      mass, ctaus = ['15'], ['1000']
      data = ms.arrs
      sel_dt = land(data['gLLP_decay_vertex_abs_z'] < 661, data['gLLP_decay_vertex_abs_r'] > 380,
                    data['gLLP_decay_vertex_r'] < 738)

      sel_csc = land(data['gLLP_abs_eta'] < 2.4, data['gLLP_decay_vertex_abs_r'] < 695.5,
                     data['gLLP_decay_vertex_z'] < 1100, data['gLLP_decay_vertex_abs_z'] > 661)

      me11 = land(data['gLLP_abs_eta'] < 2.4, data['gLLP_decay_vertex_abs_r'] < 270,
                  data['gLLP_decay_vertex_abs_z'] < 661, data['gLLP_decay_vertex_abs_z'] > 500)

      sel_csc = lor(sel_csc, me11)  # definition to make sure no overlap, BUG?
      sel_ev_cross = land(np.sum(sel_csc, axis=1) == 1, np.sum(sel_dt, axis=1) == 1)
      sel_ev_csc2 = np.sum(sel_csc, axis=1) == 2
      sel_ev_dt2 = np.sum(sel_dt, axis=1) == 2
      sel_ev_csc1 = np.sum(sel_csc, axis=1) >= 1
      sel_ev_dt1 = np.sum(sel_dt, axis=1) >= 1
      sel_ev_2tag = lor(sel_ev_cross, sel_ev_csc2, sel_ev_dt2)
      sel_ev_2tag_no_dt2 = lor(sel_ev_cross, sel_ev_csc2)

      sel_met = np.ones_like(sel_ev_cross)
      if met:
        sel_met = data['met'] > 200
      accep_dt1 = np.sum(data['weight'][land(sel_met, sel_ev_dt1)]) / np.sum(data['weight'])
      accep_csc1 = np.sum(data['weight'][land(sel_met, sel_ev_csc1)]) / np.sum(data['weight'])
      accep_dt2 = np.sum(data['weight'][land(sel_met, sel_ev_dt2)]) / np.sum(data['weight'])
      accep_csc2 = np.sum(data['weight'][land(sel_met, sel_ev_csc2)]) / np.sum(data['weight'])
      accep_cross = np.sum(data['weight'][land(sel_met, sel_ev_cross)]) / np.sum(data['weight'])

      accep_2tag = accep_dt2 + accep_csc2 + accep_cross
      accep_2tag_no_dt2 = accep_csc2 + accep_cross

      names = [
          '1 tag DT', '1 tag CSC', '2 tag total', '2 tag DT', '2 tag CSC', '1 tag DT + 1 tag CSC', 'CSC/CSC + CSC/DT'
      ]
      variables = [accep_dt1, accep_csc1, accep_2tag, accep_dt2, accep_csc2, accep_cross, accep_2tag_no_dt2]
      selections = [sel_ev_dt1, sel_ev_csc1, sel_ev_2tag, sel_ev_dt2, sel_ev_csc2, sel_ev_cross, sel_ev_2tag_no_dt2]

      for j, m in enumerate(mass):

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
                  cond = land(sel_met, selections[i])
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
              x.append(int(ct))
              y.append(sig / np.sum(data['weight']))

          x = np.array(x)
          y = np.array(y)
          cond = y > 0.0
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
          gc.append(h[k])
          if '156' == k:
            print('here')
            htots.append(h[k])

        # with rt.TFile(out_dir + 'out.root', 'CREATE') as fout:
        #   for k, v in h.items():
        #     k = k + '_multitags' + ('_MET' if met else '')
        #     fout.write(k, v)

      leg.Draw()
      c.SetLogy()
      c.SetLogx()
      c.SetTicky(1)
      c.SetTickx(1)

      c.Print(out_dir + 'ctau_accp_multitags' + ('_MET' if met else '') + file_comment)

    #
    c = canvas()
    lat, leg = get_lat_leg(leg_coords=(0.2, 0.85, 0.3, 0.92))
    for i, hh in enumerate(htots):
      hh.SetLineColor(std_color_list[i])
      hh.Draw('CA' if i == 0 else 'Csame')
      hh.SetLineColor(std_color_list[i])
      hh.SetLineWidth(3)
      hh.GetXaxis().SetTitleOffset(1)
      hh.GetYaxis().SetTitleSize(0.05)
      hh.GetYaxis().SetTitleOffset(1.5)
      hh.GetXaxis().SetLimits(0.01, 1000.0)
      hh.GetYaxis().SetRangeUser(1e-8, 1)
      leg.AddEntry(hh, 'CSC/CSC + CSC/DT' + (', MET > 200 GeV' if i else ''), "L")

    leg.Draw()
    c.SetLogy()
    c.SetLogx()
    c.SetTicky(1)
    c.SetTickx(1)

    lat.SetTextColor(rt.kBlack)
    lat.SetTextSize(0.04)
    lat.DrawLatexNDC(0.9, 0.9, 'm_{X} = 15 GeV')

    c.Print(out_dir + 'ctau_accp_met_compare' + file_comment)

    #
    c = canvas()

  ####################
