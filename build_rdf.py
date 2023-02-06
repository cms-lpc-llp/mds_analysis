'''Parses CMS TTree and reduces it into a RDF

This script allows the user to process the MuonSystem TTree and
extract the pertinent columns for the Long Lived Particle analysis
in the Muon Endcap.

'''

import ROOT as rt
import numpy as np
# from src import CMS_lumi, tdrstyle
# import numba as nb

bad_cols = ('cscRechitCluster2ZLep1LooseIso', 'dtRechitCluster2ZLep1LooseIso', 'dtRechitCluster22ZLep1LooseIso', 'dtRechitClusterJetVetoPtJESDown', 'dtRechitClusterTightJetVetoPtJESDown')

@rt.Numba.Declare(['RVec<bool>'], 'int')
def _fix_n_col(cut_col):
  return np.sum(cut_col)


@rt.Numba.Declare(['RVec<float>', 'RVec<bool>'], 'RVec<float>')
def _fix_vector_col_f(col, cut_col):
  return col[cut_col]


@rt.Numba.Declare(['RVec<int>', 'RVec<bool>'], 'RVec<int>')
def _fix_vector_col_i(col, cut_col):
  return col[cut_col]


@rt.Numba.Declare(['RVec<bool>', 'RVec<bool>'], 'RVec<bool>')
def _fix_vector_col_b(col, cut_col):
  return col[cut_col]


def update_cols(targ_rdf, cut: str, cut_col: str, col_type: str, n_cluster: int=0):
  if cut_col in targ_rdf.GetColumnNames():
    targ_rdf = targ_rdf.Redefine(cut_col, cut)
  else:
    targ_rdf = targ_rdf.Define(cut_col, cut)

  if col_type == 'jet':
    n_col_name, cut_prefix = 'nJets', 'jet'
  elif col_type == 'csc':
    n_col_name, cut_prefix = 'nCscRechitClusters', 'cscRechitCluster'
  elif col_type == 'dt':
    n_col_name, cut_prefix = 'nDtRechitClusters', 'dtRechitCluster'
  else:
    raise ValueError(col_type)

  for col in targ_rdf.GetColumnNames():
    isvec = 'RVec' in targ_rdf.GetColumnType(col)
    ctype = targ_rdf.GetColumnType(col).split('<')[-1].split('>')[0]
    if str(col)[:len(cut_prefix)] == cut_prefix and isvec and ctype in (
        'Float_t', 'Int_t', 'Bool_t') and col not in bad_cols:
      targ_rdf = targ_rdf.Redefine(col, f'Numba::_fix_vector_col_{ctype[0].lower()}({col}, {cut_col})')

  # print(cut, cut_col, col_type, cut_prefix, n_col_name, targ_rdf.Sum(n_col_name).GetValue())
  targ_rdf = targ_rdf.Redefine(n_col_name, f'Numba::_fix_n_col({cut_col})')
  # print('\t', targ_rdf.Sum(n_col_name).GetValue())
  # print('')


  if n_cluster:
    targ_rdf = targ_rdf.Filter(f'(nCscRechitClusters + nDtRechitClusters) == {n_cluster}')

  return targ_rdf


def apply_ML_cut(targ_rdf, pass_val: int=1):
  return targ_rdf.Filter(f'ML_cut == {pass_val}')

def apply_met_cut(targ_rdf, min_met: float = 200):
  return targ_rdf.Filter(f'{min_met} <= met')


def apply_lepton_cut(targ_rdf, n_lep: int = 0):
  return targ_rdf.Filter(f'nLeptons == {n_lep}')


@rt.Numba.Declare(['RVec<float>', 'RVec<float>', 'float', 'float'], 'RVec<bool>')
def _jet_cut(jet_pt, jet_eta, min_pt, max_eta):
  return (min_pt < jet_pt) & (np.abs(jet_eta) < max_eta)


def apply_jet_cut(targ_rdf, n: int = 1, min_pt: float = 50, max_eta: float = 2.4, cut_col: str = 'temp_cut_col'):
  targ_rdf = targ_rdf.Filter('nJets > 0')
  targ_rdf = update_cols(targ_rdf, f'Numba::_jet_cut(jetPt, jetEta, {min_pt}, {max_eta})', cut_col, 'jet')
  return targ_rdf.Filter(f'nJets >= {n}')


@rt.Numba.Declare(['RVec<float>', 'RVec<float>', 'RVec<int>'], 'RVec<bool>')
def _cbid_cut(csc_eta, csc_avgstation, csc_nstation):
  # if csc_nstation > 1:
  #   return np.abs(csc_eta) < 1.9
  # elif csc_avgstation >= 4:
  #   return np.abs(csc_eta) < 1.8
  # elif csc_avgstation >= 3:
  #   return np.abs(csc_eta) < 1.6
  # elif csc_avgstation >= 2:
  #   return np.abs(csc_eta) < 1.6
  # # elif csc_avgstation >= 1:
  # return np.abs(csc_eta) < 1.1
  n, a, e = csc_nstation, np.floor(csc_avgstation), csc_eta
  return (n > 1) * (e < 1.9) + \
          (n == 1) * ((a == 4) * (e < 1.8) + \
                      (a == 3) * (e < 1.6) + \
                      (a == 2) * (e < 1.6) + \
                      (a == 1) * (e < 1.1))


def apply_cbid_cut(targ_rdf, cut_col: str = 'temp_cut_col', n_cluster: int=0):
  targ_rdf = update_cols(
      targ_rdf, 'Numba::_cbid_cut(cscRechitClusterEta, cscRechitClusterAvgStation10, cscRechitClusterNStation10)',
      cut_col, 'csc', n_cluster=2)
  return targ_rdf#.Filter('nCscRechitClusters + nDtRechitClusters > 0')


@rt.Numba.Declare(['RVec<int>', 'int', 'int'], 'RVec<bool>')
def _size_cut(size, min_size, max_size):
  return (min_size <= size) & (size < max_size)


def apply_size_cut(targ_rdf, min_size: int = 50, max_size: int = 100000, col_type: str='csc', cut_col: str = 'temp_cut_col', n_cluster: int=0):
  targ_rdf = update_cols(targ_rdf, f'Numba::_size_cut({col_type}RechitClusterSize, {min_size}, {max_size})', cut_col,
                         col_type, n_cluster=2)
  #
  return targ_rdf#.Filter('nCscRechitClusters + nDtRechitClusters > 0')


@rt.Numba.Declare(['RVec<float>', 'float'], 'RVec<bool>')
def _csc_eta_cut(csc_eta, max_eta):
  return csc_eta < max_eta


def apply_csc_eta_cut(targ_rdf, max_eta: float = 2, cut_col: str = 'temp_cut_col', n_cluster: int=0):
  targ_rdf = update_cols(targ_rdf, f'Numba::_csc_eta_cut(cscRechitClusterEta, {max_eta})', cut_col, 'csc', n_cluster=2)
  return targ_rdf#.Filter('nCscRechitClusters + nDtRechitClusters > 0')


@rt.Numba.Declare(['RVec<float>', 'float', 'float'], 'RVec<bool>')
def _time_cut(csc_time, min_t, max_t):
  return (min_t < csc_time) & (csc_time < max_t)


def apply_time_cut(targ_rdf, min_t: float = -5, max_t: float = 12.5, col_type: str='csc', cut_col: str = 'temp_cut_col', n_cluster: int=0):
  # targ_rdf = update_cols(targ_rdf, f'Numba::_time_cut({col_type}RechitClusterTime, {min_t}, {max_t})', cut_col, col_type, n_cluster=2)
  return targ_rdf#.Filter('nCscRechitClusters + nDtRechitClusters > 0')


@rt.Numba.Declare(['RVec<float>', 'float'], 'RVec<bool>')
def _csc_time_spread_cut(csc_time_spread, max_spread):
  return np.abs(csc_time_spread) < max_spread


def apply_csc_time_spread_cut(targ_rdf, max_spread: float = 20, cut_col: str = 'temp_cut_col', n_cluster: int=0):
  targ_rdf = update_cols(targ_rdf, f'Numba::_csc_time_spread_cut(cscRechitClusterTimeSpread, {max_spread})', cut_col,
                         'csc', n_cluster=2)
  return targ_rdf#.Filter('nCscRechitClusters + nDtRechitClusters > 0')


@rt.Numba.Declare(['RVec<float>', 'float'], 'RVec<bool>')
def _csc_oot_cut(csc_time, max_t):
  return csc_time < max_t


def apply_csc_oot_cut(targ_rdf, max_t: float = -12.5, cut_col: str = 'temp_cut_col', n_cluster: int=0):
  targ_rdf = update_cols(targ_rdf, f'Numba::_csc_oot_cut(cscRechitClusterTime, {max_t})', cut_col, 'csc', n_cluster=2)
  return targ_rdf#.Filter('nCscRechitClusters + nDtRechitClusters > 0')


@rt.Numba.Declare(['RVec<float>', 'RVec<float>', 'RVec<float>', 'RVec<float>', 'float'], 'RVec<bool>')
def _csc_jet_cut(csc_eta, csc_phi, jet_eta, jet_phi, min_dr):
  cut = np.abs(csc_eta) > 0
  for icsc, (ceta, cphi) in enumerate(zip(csc_eta, csc_phi)):
    for jeta, jphi in zip(jet_eta, jet_phi):
      cut[icsc] = cut[icsc] & (np.sqrt((ceta - jeta)**2 + (cphi - jphi)**2) > min_dr)
  return cut


def apply_csc_jet_cut(targ_rdf, min_dr: float = 0.4, cut_col: str = 'temp_cut_col', n_cluster: int=0):
  targ_rdf = update_cols(targ_rdf,
                         f'Numba::_csc_jet_cut(cscRechitClusterEta, cscRechitClusterPhi, jetEta, jetPhi, {min_dr})',
                         cut_col, 'csc', n_cluster=2)
  return targ_rdf#.Filter('nCscRechitClusters + nDtRechitClusters > 0')


@rt.Numba.Declare(['RVec<float>', 'RVec<float>', 'RVec<float>', 'RVec<float>', 'RVec<int>', 'float'], 'RVec<bool>')
def _csc_muon_cut(csc_eta, csc_phi, lep_eta, lep_phi, lep_pid, min_dr):
  cut = np.abs(csc_eta) > 0
  for icsc, (ceta, cphi) in enumerate(zip(csc_eta, csc_phi)):
    for lpid, leta, lphi in zip(lep_pid, lep_eta, lep_phi):
      if abs(lpid) == 13:
        cut[icsc] = cut[icsc] & (np.sqrt((ceta - leta)**2 + (cphi - lphi)**2) > min_dr)
  return cut


def apply_csc_muon_cut(targ_rdf, min_dr: float = 0.4, cut_col: str = 'temp_cut_col', n_cluster: int=0):
  targ_rdf = update_cols(
      targ_rdf, f'Numba::_csc_muon_cut(cscRechitClusterEta, cscRechitClusterPhi, lepEta, lepPhi, lepPdgId, {min_dr})',
      cut_col, 'csc', n_cluster=2)
  return targ_rdf#.Filter('nCscRechitClusters + nDtRechitClusters > 0')


@rt.Numba.Declare(['RVec<int>', 'RVec<int>', 'RVec<int>', 'RVec<int>', 'RVec<int>', 'RVec<int>', 'RVec<int>'],
                  'RVec<bool>')
def _endcap_nrechit_cut(me11, me12, p11, p12, re12, mb1, rb1):
  return (me11 + me12 + p11 + p12 + re12 + mb1 + rb1) == 0


def apply_endcap_nrechit_cut(targ_rdf, cut_col: str = 'temp_cut_col', n_cluster: int=0):
  targ_rdf = update_cols(
      targ_rdf,
      'Numba::_endcap_nrechit_cut(cscRechitClusterNRechitChamberMinus11, cscRechitClusterNRechitChamberMinus12, cscRechitClusterNRechitChamberPlus11, cscRechitClusterNRechitChamberPlus12, cscRechitCluster_match_RE12_0p4, cscRechitCluster_match_MB1Seg_0p4, cscRechitCluster_match_RB1_0p4)',
      cut_col, 'csc', n_cluster=2)
  return targ_rdf#.Filter('nCscRechitClusters + nDtRechitClusters > 0')


def apply_match_cut(targ_rdf, inverse: bool=False, cut_col: str = 'temp_cut_col', n_cluster: int=0):
  if inverse:
    targ_rdf = update_cols(targ_rdf, '!cscRechitCluster_match_gLLP', cut_col, 'csc', n_cluster=2)
    targ_rdf = update_cols(targ_rdf, '!dtRechitCluster_match_gLLP', cut_col, 'dt', n_cluster=2)
  else:
    targ_rdf = update_cols(targ_rdf, 'cscRechitCluster_match_gLLP', cut_col, 'csc', n_cluster=2)
    targ_rdf = update_cols(targ_rdf, 'dtRechitCluster_match_gLLP', cut_col, 'dt', n_cluster=2)

  return targ_rdf#.Filter('nCscRechitClusters + nDtRechitClusters > 0')


def apply_all_cuts(targ_rdf, met=200, only_event_level=False, require_2_clusters=False, oot=False, n_cluster=2):
  # Event Level
  targ_rdf = apply_met_cut(targ_rdf, met)
  targ_rdf = apply_lepton_cut(targ_rdf)
  targ_rdf = apply_jet_cut(targ_rdf)

  # Cluster Level
  if not only_event_level:
    targ_rdf = apply_csc_muon_cut(targ_rdf)
    targ_rdf = apply_csc_jet_cut(targ_rdf)
    targ_rdf = apply_endcap_nrechit_cut(targ_rdf)
    targ_rdf = apply_size_cut(targ_rdf)
    targ_rdf = apply_cbid_cut(targ_rdf)
    targ_rdf = apply_csc_eta_cut(targ_rdf)
    targ_rdf = apply_csc_time_spread_cut(targ_rdf)

    if require_2_clusters:
      targ_rdf = targ_rdf.Filter('nCscRechitClusters == 2')

    if oot:  #Out of time
      targ_rdf = apply_csc_oot_cut(targ_rdf, n_cluster=n_cluster)
    else:
      targ_rdf = apply_time_cut(targ_rdf, n_cluster=n_cluster)

  return targ_rdf#.Filter('nCscRechitClusters + nDtRechitClusters > 0')


#########################################################
@rt.Numba.Declare(['float', 'RVec<float>', 'RVec<float>'], 'RVec<float>')
def add_delta_e(higgs_e, csc_e, jet_e):
  return np.array([higgs_e - ce - je for ce in csc_e for je in jet_e])


@rt.Numba.Declare(['int', 'RVec<int>'], 'int')
def add_size_col(col_idx, csc_size):
  if len(csc_size) > 1:
    # return csc_size[col_idx]
    if col_idx == 0:
      return np.min(csc_size)
    # elif col_idx == 1:
    return np.max(csc_size)
  else:
    return -1


#########################################################

# - [ ] METNoMu tiggers
# - [ ] Cluster RecHits with η-ɸ, distance parameter ΔR = 0.2
# - [ ] Merge clusters if two clusters are within ΔR < 0.6

if __name__ == '__main__':
  date = 'jan30'
  file_prefix = 'simu'

  rt.gROOT.SetBatch(True)
  gc = []
  req_2 = False

  # tdr = tdrstyle.setTDRStyle()
  # CMS_lumi.writeExtraText = True
  lat = rt.TLatex()
  lat.SetTextSize(0.05)
  lat.SetTextAlign(11)
  rt.gStyle.SetOptFit(1011)
  rt.gStyle.SetOptStat(0)

  data_year = 'full'
  years = ['Run3_Fall22']
  fpath = 'data/raw/ggH_HToSSTobbbb_MH-125_MS-15_ctau-1000_TuneCP5_13TeV-powheg-pythia8_59740pb_weighted.root'
  # fpath = 'data/raw/WITH_ML_ggH_HToSSTobbbb_MH-125_MS-15_ctau-1000_TuneCP5_13TeV-powheg-pythia8_59740pb_weighted.root'
  tree_name = 'MuonSystem'

  rdf = rt.RDataFrame(tree_name, fpath)
  # ! Always include the CSC max size 100 cut! This blinds the data
  # TODO: Do this before any analysis and save it as a separate data file!!!
  rdfs = {'raw': {'rdf': rdf, 'line_color': rt.kBlack}}  # Blinded data
  # rdfs = {'raw': {'rdf': apply_csc_size_cut(rdf, 0, 100), 'line_color': rt.kBlack}}  # Blinded data
  # rdfs = {'raw': {'rdf': apply_csc_size_cut(rdf, 100), 'line_color': rt.kBlack}}  # Signal

  rdfs['met50'] = {'rdf': apply_met_cut(rdfs['raw']['rdf'], 50), 'line_color': rt.kRed}
  rdfs['met100'] = {'rdf': apply_met_cut(rdfs['raw']['rdf'], 100), 'line_color': rt.kRed}
  rdfs['met150'] = {'rdf': apply_met_cut(rdfs['raw']['rdf'], 150), 'line_color': rt.kRed}
  rdfs['met200'] = {'rdf': apply_met_cut(rdfs['raw']['rdf'], 200), 'line_color': rt.kRed}
  rdfs['met250'] = {'rdf': apply_met_cut(rdfs['raw']['rdf'], 250), 'line_color': rt.kRed}
  rdfs['lep0'] = {'rdf': apply_lepton_cut(rdfs['raw']['rdf']), 'line_color': rt.kGreen}
  rdfs['jet'] = {'rdf': apply_jet_cut(rdfs['raw']['rdf']), 'line_color': rt.kBlue}
  rdfs['csc_jet'] = {'rdf': apply_csc_jet_cut(rdfs['raw']['rdf']), 'line_color': rt.kBlue}
  rdfs['csc_muon'] = {'rdf': apply_csc_muon_cut(rdfs['raw']['rdf']), 'line_color': rt.kBlue}
  rdfs['endcap'] = {'rdf': apply_endcap_nrechit_cut(rdfs['raw']['rdf']), 'line_color': rt.kBlue}
  rdfs['cbid'] = {'rdf': apply_cbid_cut(rdfs['raw']['rdf']), 'line_color': rt.kMagenta}
  rdfs['csc_size'] = {'rdf': apply_size_cut(rdfs['raw']['rdf']), 'line_color': rt.kMagenta}
  rdfs['csc_eta'] = {'rdf': apply_csc_eta_cut(rdfs['raw']['rdf']), 'line_color': rt.kBlack}
  rdfs['csc_time'] = {'rdf': apply_time_cut(rdfs['raw']['rdf']), 'line_color': rt.kBlack}
  rdfs['csc_time_spread'] = {'rdf': apply_csc_time_spread_cut(rdfs['raw']['rdf']), 'line_color': rt.kBlack}
  rdfs['csc_oot'] = {'rdf': apply_csc_oot_cut(rdfs['raw']['rdf']), 'line_color': rt.kPink}
  rdfs['matched'] = {'rdf': apply_match_cut(rdfs['raw']['rdf']), 'line_color': rt.kPink}

  # rdfs['pass_ML'] = {'rdf': apply_ML_cut(rdfs['raw']['rdf'], 1), 'line_color': rt.kPink}
  # rdfs['fail_ML'] = {'rdf': apply_ML_cut(rdfs['raw']['rdf'], 0), 'line_color': rt.kPink}
  # rdfs['pass_ML_matched'] = {'rdf': apply_match_cut(apply_ML_cut(rdfs['raw']['rdf'], 1)), 'line_color': rt.kPink}
  # rdfs['fail_ML_matched'] = {'rdf': apply_match_cut(apply_ML_cut(rdfs['raw']['rdf'], 0)), 'line_color': rt.kPink}

  rdfs['event_level_met50'] = {
      'rdf': apply_all_cuts(rdfs['raw']['rdf'], 50, only_event_level=True, require_2_clusters=req_2),
      'line_color': rt.kYellow - 3,
  }
  rdfs['event_level_met50_matched'] = {
      'rdf': apply_match_cut(rdfs['event_level_met50']['rdf']),
      'line_color': rt.kYellow - 3,
  }

  rdfs['ALL_met50'] = {
      'rdf': apply_all_cuts(rdfs['raw']['rdf'], 50, require_2_clusters=req_2),
      'line_color': rt.kYellow - 3,
  }
  rdfs['ALL_met50_matched'] = {
      'rdf': apply_match_cut(rdfs['ALL_met50']['rdf']),
      'line_color': rt.kYellow - 3,
  }

  rdfs['OOT_ALL_met50'] = {
      'rdf': apply_all_cuts(rdfs['raw']['rdf'], 50, require_2_clusters=req_2, oot=True),
      'line_color': rt.kYellow - 3,
  }

  ##############################################3

  vals = [
      ('gHiggsE', 'gHiggsE', 100, 100, 1600),
      ('met', 'met', 100, 0, 800),
      ('cscRechitCluster_match_gLLP_e', 'cscRechitCluster_match_gLLP_e', 100, 0, 800),
      ('jetE', 'jetE', 100, 0, 800),
      ('deltaE', 'deltaE', 100, -250, 250),
      ('cscRechitClusterTime', 'cscRechitClusterTime', 100, -15, 15),
      ('nJets', 'nJets', 8, -0.5, 7.5),
      ('nCscRechitClusters', 'nCscRechitClusters', 7, -0.5, 6.5),
      ('cscRechitClusterPhi', 'cscRechitClusterPhi', 100, -np.pi, np.pi),
      ('cscRechitCluster_match_gLLP_phi', 'cscRechitCluster_match_gLLP_phi', 100, -np.pi, np.pi),
      ('cscRechitClusterEta', 'cscRechitClusterEta', 100, -3, 3),
      ('cscRechitCluster_match_gLLP_eta', 'cscRechitCluster_match_gLLP_eta', 100, -3, 3),
      ('cscRechitClusterZ', 'cscRechitClusterZ', 100, -1000, 1000),
      ('cscRechitCluster_match_gLLP_decay_z', 'cscRechitCluster_match_gLLP_decay_z', 100, -1000, 1000),
      # ('cscRechitClusterR', 'cscRechitClusterR', 100, 0, 1000),
      ('cscRechitCluster_match_gLLP_decay_r', 'cscRechitCluster_match_gLLP_decay_r', 100, 0, 1000),
      ('gHiggsPt', 'gHiggsPt', 100, 0, 500),
      ('gHiggsEta', 'gHiggsEta', 100, -12.5, 12.5),
      ('gHiggsPhi', 'gHiggsPhi', 100, -np.pi, np.pi),
      ('jetPt', 'jetPt', 100, 0, 500),
      ('jetEta', 'jetEta', 100, -3, 3),
      ('jetPhi', 'jetPhi', 100, -np.pi, np.pi),
      ('cscRechitClusterSize', 'cscRechitClusterSize', 100, 0, 3000),
      ('cscRechitClusterNStation10', 'cscRechitClusterNStation10', 6, -0.5, 5.5),
      ('cscRechitClusterMaxStation', 'cscRechitClusterMaxStation', 12, -5.5, 5.5),
      ('cscRechitClusterAvgStation10', 'cscRechitClusterAvgStation10', 50, -4, 4),
      ('nDtRechitClusters', 'nDtRechitClusters', 7, -0.5, 6.5),
      ('nDTRechits', 'nDTRechits', 100, 0, 3000),
      ('nCscRechits', 'nCscRechits', 100, 0, 3000),
  ]

  # npx, width = 500, 4
  # c1 = rt.TCanvas('c1', 'c1', npx * width, npx * (len(vals) + width - 1) // width)
  # c1.Divide(width, (len(vals) + width - 1) // width)
  # c1.Draw()

  for k, rdf_info in rdfs.items():

    print(f'{k:>20} - {rdf_info["rdf"].Count().GetValue():,}')
    print(f'{"weight":>20} - {rdf_info["rdf"].Sum("weight").GetValue()*25/59:,.0f}')

    print('APPLYING 2 CLUSTER CUT (DT + CSC)')
    rdf_info['rdf'] = rdf_info['rdf'].Filter('(nCscRechitClusters + nDtRechitClusters) >= 2')
    # rdf_info['rdf'] = rdf_info['rdf'].Filter('nCscRechitClusters == 2')
    
    print(f'{"2 cluster":>20} - {rdf_info["rdf"].Count().GetValue():,}')
    print(f'{"weight":>20} - {rdf_info["rdf"].Sum("weight").GetValue()*25/59:,.0f}')
    print('')

    rdf_info['rdf'] = rdf_info['rdf'].Define('deltaE',
                                             'Numba::add_delta_e(gHiggsE, cscRechitCluster_match_gLLP_e, jetE)')

    rdf_info['rdf'] = rdf_info['rdf'].Define('csc_size0', 'Numba::add_size_col(0, cscRechitClusterSize)')
    rdf_info['rdf'] = rdf_info['rdf'].Define('csc_size1', 'Numba::add_size_col(1, cscRechitClusterSize)')

    rdf_info['rdf'] = rdf_info['rdf'].Define('delta_metPhi_cscRechitClusterPhi', 'metPhi - cscRechitClusterPhi')
    rdf_info['rdf'] = rdf_info['rdf'].Define('cscRechitClusterAbsEta', 'abs(cscRechitClusterEta)')
    rdf_info['rdf'] = rdf_info['rdf'].Define('cscRechitCluster_match_gLLP_decay_abs_z',
                                             'abs(cscRechitCluster_match_gLLP_decay_z)')
    rdf_info['rdf'] = rdf_info['rdf'].Define('cscRechitClusterAbsZ', 'abs(cscRechitClusterZ)')
    rdf_info['rdf'] = rdf_info['rdf'].Define('gLLP_decay_vertex_abs_z', 'abs(gLLP_decay_vertex_z)')
    rdf_info['rdf'] = rdf_info['rdf'].Define('cscRechitCluster_match_gLLP_decay_abs_r',
                                             'abs(cscRechitCluster_match_gLLP_decay_r)')
    rdf_info['rdf'] = rdf_info['rdf'].Define('gLLP_decay_vertex_abs_r', 'abs(gLLP_decay_vertex_r)')

  #   for iv, v in enumerate(vals):
  #     c1.cd(iv + 1).SetLogy()
  #     c1.cd(iv + 1).SetGrid()
  #     if len(v) > 1:
  #       hh = rdf_info['rdf'].Histo1D(v, v[0])
  #     else:
  #       hh = rdf_info['rdf'].Histo1D(v[0])
  #     hh.SetLineColor(rdf_info['line_color'])
  #     hh.SetLineWidth(2)
  #     hh.SetMinimum(1)
  #     hh.Draw('same')
  #     gc.append(hh)

  #     c1.Print('test.pdf')


  # Cluster vs Cluster
  c1 = rt.TCanvas('c1', 'c1', 800, 400)
  c1.Divide(2, 1)
  c1.Draw()
  
  c1.cd(1).SetLogz()
  c1.cd(1).SetGrid()
  h1 = rdfs['raw']['rdf'].Histo2D(('',';csc;dt',7, -0.5, 6.5, 7, -0.5, 6.5),'nCscRechitClusters','nDtRechitClusters')
  h1.Draw('coltext')

  c1.cd(2).SetLogz()
  c1.cd(2).SetGrid()
  h2 = rdfs['matched']['rdf'].Histo2D(('',';csc;dt',7, -0.5, 6.5, 7, -0.5, 6.5),'nCscRechitClusters','nDtRechitClusters')
  h2.Draw('coltext')

  c1.Print(f'reports/weekly/{date}/{file_prefix}_ncscclusters_ndtclusters.png')

  ### MET Scan
  c1 = rt.TCanvas('c1', 'c1', 800, 800)
  c1.Draw()
  c1.cd(1).SetLogy()
  c1.cd(1).SetGrid()

  template = (('delta_metPhi_cscRechitClusterPhi', 'delta_metPhi_cscRechitClusterPhi;\phi_{ME_{t}} - \phi_{CSC}', 100, -np.pi, np.pi),
              'delta_metPhi_cscRechitClusterPhi')

  hhs = [ rdfs['raw']['rdf'].Histo1D(*template),
          rdfs['met50']['rdf'].Histo1D(*template),
          rdfs['met100']['rdf'].Histo1D(*template),
          rdfs['met150']['rdf'].Histo1D(*template),
          rdfs['met200']['rdf'].Histo1D(*template),
          rdfs['met250']['rdf'].Histo1D(*template),
          # rdfs['pass_ML']['rdf'].Histo1D(*template),
          # rdfs['fail_ML']['rdf'].Histo1D(*template)
          ]

  hms = [  apply_match_cut(rdfs['raw']['rdf']).Histo1D(*template),
            apply_match_cut(rdfs['met50']['rdf']).Histo1D(*template),
            apply_match_cut(rdfs['met100']['rdf']).Histo1D(*template),
            apply_match_cut(rdfs['met150']['rdf']).Histo1D(*template),
            apply_match_cut(rdfs['met200']['rdf']).Histo1D(*template),
            apply_match_cut(rdfs['met250']['rdf']).Histo1D(*template),
            # apply_match_cut(rdfs['pass_ML']['rdf']).Histo1D(*template),
            # apply_match_cut(rdfs['fail_ML']['rdf']).Histo1D(*template) 
            ]

  ccs = [ rt.kBlack,
          rt.kRed,
          rt.kBlue,
          rt.kGreen,
          rt.kCyan,
          rt.kMagenta,
          # rt.kRed+3,
          # rt.kRed-5 
          ]

  tts = [
    'Signal (No cuts)',
    'MET >  50 GeV',
    'MET > 100 GeV',
    'MET > 150 GeV',
    'MET > 200 GeV',
    'MET > 250 GeV',
    # 'Pass ML',
    # 'Fail ML'
  ]
  
  lat.SetTextAlign(33)
  for ih, (hh, hm, cc, tt) in enumerate(zip(hhs, hms, ccs, tts)):
    hh.SetMinimum(1)
    hh.SetLineColor(cc)
    hh.SetLineWidth(2)
    hh.Draw('same')

    hm.SetLineColor(cc)
    hm.SetLineWidth(2)
    hm.SetLineStyle(2)
    hm.Draw('same')

    lat.SetTextColor(cc)
    lat.DrawLatexNDC(0.95, 0.9 - ih * 0.05, tt)

  c1.Print(f'reports/weekly/{date}/{file_prefix}_delta_metPhi_cscRechitClusterPhi_met_scan.png')


  ### deltaE
  c1 = rt.TCanvas('c1', 'c1', 800, 800)
  c1.Divide(1, 1)
  c1.Draw()

  templates = [
      (('deltaE', 'deltaE', 100, -200, 200), 'deltaE'),
  ]

  for it, template in enumerate(templates):
    c1.cd(it + 1).SetGrid()
    h0 = rdfs['raw']['rdf'].Histo1D(*template)
    h1 = rdfs['met50']['rdf'].Histo1D(*template)
    h2 = rdfs['matched']['rdf'].Histo1D(*template)
    h3 = rdfs['ALL_met50']['rdf'].Histo1D(*template)
    h4 = rdfs['ALL_met50_matched']['rdf'].Histo1D(*template)
    h5 = rdfs['csc_oot']['rdf'].Histo1D(*template)

    hhs = [h0, h1, h2, h3, h4, h5]
    ccs = [rt.kBlack, rt.kYellow - 2, rt.kGreen, rt.kRed, rt.kMagenta, rt.kBlue]
    tts = ['Raw', '50 < MET', 'Matched', 'All Cuts, 50 < Met', '^^ + matched', 'OOT (t < -12.5)']
    for ih, (hh, cc, tt) in enumerate(zip(hhs, ccs, tts)):
      hh.SetLineColor(cc)
      hh.DrawNormalized('same')
      lat.SetTextAlign(33)
      lat.SetTextColor(cc)
      lat.DrawLatexNDC(0.95, 0.9 - ih * 0.05, tt)

    gc.extend(hhs)

  c1.Print(f'reports/weekly/{date}/{file_prefix}_deltaE.png')

  ### Efficiency Plots
  c1 = rt.TCanvas('c', 'c', 2 * 1200, 2 * 800)
  c1.Divide(2, 2)
  c1.Draw()
  vcsc_pre, vgen_pre = 'cscRechitCluster_match_gLLP_decay_abs_', 'gLLP_decay_vertex_abs_'

  nb, zbins, rbins = 100, (400, 1100), (100, 750)
  for i, (vv, bins) in enumerate(zip(['z', 'r'], [zbins, rbins])):
    c1.cd(i + 1).SetGrid()
    vcsc, vgen = vcsc_pre + vv, vgen_pre + vv

    xmin, xmax = bins
    model = ('efficiency_decay_' + vv, f' CSC |{vv.upper()}| Decay Vertex Efficiency;|{vv.upper()}| [cm];efficiency',
             nb, xmin, xmax)

    hh_raw = rdfs['matched']['rdf'].Histo1D(model, vcsc)
    hh_raw.Divide(rdfs['raw']['rdf'].Histo1D(model, vgen).GetPtr())

    # hh_met = rdfs['pass_ML']['rdf'].Histo1D(model, vcsc)
    # hh_met.Divide(rdfs['raw']['rdf'].Histo1D(model, vgen).GetPtr())
    hh_met = rdfs['ALL_met50']['rdf'].Histo1D(model, vcsc)
    hh_met.Divide(rdfs['event_level_met50']['rdf'].Histo1D(model, vgen).GetPtr())

    rmax, mmax = hh_raw.GetMaximum(), hh_met.GetMaximum()

    hh_raw.SetLineColor(rt.kBlack)
    hh_met.SetLineColor(rt.kRed)
    if rmax > mmax:
      hh_raw.Draw('same')
      hh_met.Draw('same')
    else:
      hh_met.Draw('same')
      hh_raw.Draw('same')
    gc.extend([hh_raw, hh_met])

    print(hh_raw.Integral(), hh_met.Integral())

    lat.SetTextAlign(33)
    lat.SetTextColor(rt.kBlack)
    lat.DrawLatexNDC(0.95, 0.90, 'Raw')
    lat.SetTextColor(rt.kRed)
    lat.DrawLatexNDC(0.95, 0.85, 'Event Level Cuts, 50 < MET')

  lat.SetTextAlign(33)
  model_zr = ('', 'CSC |R| vs |Z| Decay Efficiency;|Z| [cm];|R| [cm]', 20, *zbins, 20, *rbins)

  # c1.cd(3).SetLogz()
  c1.cd(3).SetGrid()
  hh_eff_z_r_raw = rdfs['matched']['rdf'].Histo2D(model_zr, vcsc_pre + 'z', vcsc_pre + 'r')
  hh_eff_z_r_raw.Divide(rdfs['raw']['rdf'].Histo2D(model_zr, vgen_pre + 'z', vgen_pre + 'r').GetPtr())

  # hh_eff_z_r_met = rdfs['pass_ML']['rdf'].Histo2D(model_zr, vcsc_pre + 'z', vcsc_pre + 'r')
  # hh_eff_z_r_met.Divide(rdfs['raw']['rdf'].Histo2D(model_zr, vgen_pre + 'z', vgen_pre + 'r').GetPtr())
  hh_eff_z_r_met = rdfs['ALL_met50']['rdf'].Histo2D(model_zr, vcsc_pre + 'z', vcsc_pre + 'r')
  hh_eff_z_r_met.Divide(rdfs['event_level_met50']['rdf'].Histo2D(model_zr, vgen_pre + 'z', vgen_pre + 'r').GetPtr())

  zmax = max(hh_eff_z_r_raw.GetMaximum(), hh_eff_z_r_met.GetMaximum())
  hh_eff_z_r_met.SetMaximum(zmax)
  hh_eff_z_r_raw.SetMaximum(zmax)

  hh_eff_z_r_met.SetMinimum(0.05)
  hh_eff_z_r_raw.SetMinimum(0.05)

  hh_eff_z_r_raw.Draw('colz')
  lat.SetTextColor(rt.kBlack)
  lat.DrawLatexNDC(0.95, 1.00, 'Raw')

  # c1.cd(4).SetLogz()
  c1.cd(4).SetGrid()
  hh_eff_z_r_met.Draw('colz')
  lat.SetTextColor(rt.kRed)
  lat.DrawLatexNDC(0.95, 1.00, 'Event Level Cuts, 50 < MET')

  c1.Print(f'reports/weekly/{date}/{file_prefix}_efficiency.png')

  #### Cluster Time (EXO)
  c1 = rt.TCanvas('c1', 'c1', 800, 800)
  c1.Draw()
  c1.cd(1).SetGrid()

  template = (('cscRechitClusterTime', 'cscRechitClusterTime', 100, -80, 80), 'cscRechitClusterTime')

  h0 = rdfs['raw']['rdf'].Histo1D(*template)
  h1 = rdfs['met50']['rdf'].Histo1D(*template)
  h2 = rdfs['matched']['rdf'].Histo1D(*template)
  h3 = rdfs['ALL_met50']['rdf'].Histo1D(*template)
  h4 = rdfs['ALL_met50_matched']['rdf'].Histo1D(*template)

  hhs = [h0, h1, h2, h3, h4]
  ccs = [rt.kBlack, rt.kYellow - 2, rt.kGreen, rt.kRed, rt.kMagenta]
  tts = ['Raw', '50 < MET', 'Matched', 'All Cuts, 50 < Met', '^^ + matched']
  for i, (hh, cc, tt) in enumerate(zip(hhs, ccs, tts)):
    hh.SetLineColor(cc)
    hh.DrawNormalized('same')
    lat.SetTextAlign(33)
    lat.SetTextColor(cc)
    lat.DrawLatexNDC(0.95, 0.9 - i * 0.05, tt)

  c1.Print(f'reports/weekly/{date}/{file_prefix}_cluster_time.png')

  #### Cluster-Level ID (EXO)
  nr, nc = 2, 3
  c1 = rt.TCanvas('c1', 'c1', nc * 800, nr * 800)
  c1.Divide(nc, nr)
  c1.Draw()

  templates = [
      (('cscRechitClusterNStation10', 'cscRechitClusterNStation10', 6, -0.5, 5.5), 'cscRechitClusterNStation10'),
      (('cscRechitClusterAbsEta', 'cscRechitClusterAbsEta', 100, 0.8, 2), 'cscRechitClusterAbsEta'),
      (('cscRechitClusterAvgStation10', 'cscRechitClusterAvgStation10', 6, -0.5, 5.5), 'cscRechitClusterAvgStation10'),
      (('nCscRechits', 'nCscRechits', 100, 0, 3000), 'nCscRechits'),
      (('cscRechitClusterSize', 'cscRechitClusterSize', 100, 0, 3000), 'cscRechitClusterSize'),
  ]

  for it, template in enumerate(templates):
    c1.cd(it + 1).SetGrid()
    c1.cd(it + 1).SetLogy()
    h0 = rdfs['raw']['rdf'].Histo1D(*template)
    h1 = rdfs['met50']['rdf'].Histo1D(*template)
    h2 = rdfs['matched']['rdf'].Histo1D(*template)
    h3 = rdfs['ALL_met50']['rdf'].Histo1D(*template)
    h4 = rdfs['ALL_met50_matched']['rdf'].Histo1D(*template)
    h5 = rdfs['csc_oot']['rdf'].Histo1D(*template)

    hhs = [h0, h1, h2, h3, h4, h5]
    ymax = max([h.GetMaximum() * 1.1 for h in hhs])
    ccs = [rt.kBlack, rt.kYellow - 2, rt.kGreen, rt.kRed, rt.kMagenta, rt.kBlue]
    tts = ['Raw', '50 < MET', 'Matched', 'All Cuts, 50 < Met', '^^ + matched', 'OOT (t < -12.5)']
    for ih, (hh, cc, tt) in enumerate(zip(hhs, ccs, tts)):
      hh.SetLineColor(cc)
      hh.SetMinimum(1)
      hh.SetMaximum(ymax)
      hh.Draw('same')
      # hh.DrawNormalized('same')
      lat.SetTextAlign(33)
      lat.SetTextColor(cc)
      lat.DrawLatexNDC(0.95, 0.9 - ih * 0.05, tt)

    gc.extend(hhs)

  c1.Print(f'reports/weekly/{date}/{file_prefix}_cluster_level_id.png')

  #### CSC 2 Cluster Correlation
  nr, nc = 2, 2
  c1 = rt.TCanvas('c1', 'c1', nc * 800, nr * 800)
  c1.Divide(nc, nr)
  c1.Draw()

  template = (('csc_corr', 'csc_corr;Small Cluster Size;Big Cluster Size', 100, 0, 1000, 100, 0, 1000), 'csc_size0',
              'csc_size1')

  h0 = rdfs['raw']['rdf'].Histo2D(*template)
  h1 = rdfs['met50']['rdf'].Histo2D(*template)
  h2 = rdfs['matched']['rdf'].Histo2D(*template)
  h3 = rdfs['ALL_met50']['rdf'].Histo2D(*template)
  # h4 = rdfs['ALL_met50_matched']['rdf'].Histo2D(*template)

  hhs = [h0, h1, h2, h3]
  # ccs = [rt.kBlack, rt.kYellow - 2, rt.kGreen, rt.kRed, rt.kMagenta]
  tts = ['Raw', 'met50', 'Matched', 'All Cuts, 50 < Met']  # '^^ + matched'
  for i, (hh, tt) in enumerate(zip(hhs, tts)):
    c1.cd(i + 1).SetGrid()
    c1.cd(i + 1).SetLogz()
    # hh.SetLineColor(cc)
    hh.Draw('colz')
    lat.SetTextAlign(33)
    lat.DrawLatexNDC(0.95, 0.95, tt)

  c1.Print(f'reports/weekly/{date}/{file_prefix}_csc_size_2cluster_corr.png')
