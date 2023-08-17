"""train_bdt.py"""
import os
import sys
import pathlib
from collections import defaultdict

import numpy as np
import sklearn as skl

# import ROOT as rt
# from ROOT import gErrorIgnoreLevel
from ROOT import RDataFrame
from ROOT import EnableImplicitMT, GetThreadPoolSize
from ROOT.gROOT import SetBatch

# from src import CMS_lumi
# from src import tdrstyle

from src.muon_system import MuonSystemRDF, MuonSystemAwkward# , multi_plot

# from src.muon_system import (get_lat_leg, draw_csc_z_boxes, draw_dt_r_boxes)
from src.helper_functions import canvas, alert


def flatten_muon_system_rdf(ms: MuonSystemRDF) -> np.ndarray:
    """Flattens a MuonSystem RDF into an Numpy array.
    
    Args:
        ms (MuonSystemRDF): 

    Returns:
        a numpy array with cluster and physics information.

    Raises:
        ValueError: A row contains more than 2 clusters

    Flattens the muon system RDF into a numpy array."""
    return ms

# tag_dR
# tag_dEta
# tag_dPhi

############


if __name__ == '__main__':
    print('+-----------------------+')
    print('| Starting train_bdt.py |')
    print('+-----------------------+\n')

    #############################
    # Parameters -- make a yaml #
    #############################
    save_dstat = 'ca_0p6'
    save_date = 'aug10'
    N_THREADS = 4

    N_EVENTS = 100_000
    CUT = True

    #############################

    gc = []
    SetBatch()
    # EnableImplicitMT(N_THREADS)
    # alert(f'Setting {N_THREADS=:,} and {GetThreadPoolSize()=:,}')
    # print(f'Running with thread pool size = {GetThreadPoolSize():,}')

    #############################

    if len(sys.argv) > 1:
        N_EVENTS = int(sys.argv[1])
        if N_EVENTS > 0:
            alert(f'Setting {N_EVENTS=:,}')
        else:
            N_EVENTS = None
            alert(f'Setting {N_EVENTS=} (Loading ALL events)')

    #############################

    if 'caltech' in os.uname()[1]:
        ff_mc = '/storage/cms/store/user/christiw/displacedJetMuonAnalyzer/Run3/V1p19/MC_Summer22EE/v1/sixie/v6/normalized/ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted.root'
        ff_r3 = '/storage/cms/store/user/christiw/displacedJetMuonAnalyzer/Run3/V1p19/Data2022/v6/normalized/DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi.root'
    else:
        ff_mc = 'data/raw/ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v6.root'
        ff_r3 = 'data/raw/DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi_v6.root'


    ms_mc = MuonSystemAwkward(ff_mc, nev=N_EVENTS, is_mc=True)
    ms_r3 = MuonSystemAwkward(ff_r3, nev=N_EVENTS, is_mc=False)

    #!!! TURNING CUTS OFF !!!!#
    ms_mc.cut, ms_r3.cut = False, False
    #!!!!!!!!!!!!!!!!!!!!!!!!!#

    print('')
    print('--- Filtering MuonSystems ---')
    print(f'   In | MC : {ms_mc.count():>10,} | R3 : {ms_r3.count():>10,}')

    ms_mc.match_mc()
    print(f'Match | MC : {ms_mc.count():>10,} | R3 : {ms_r3.count():>10,}')

    ms_mc.cut_hlt()
    ms_r3.cut_hlt()
    print(f'  HLT | MC : {ms_mc.count():>10,} | R3 : {ms_r3.count():>10,}')

    # ms_mc.cut_l1() #! Broken rn
    # ms_r3.cut_l1()
    # print(f'Events In | {ms_mc.count()=:>10,} | {ms_r3.count()=:>10,}')

    ms_mc.cut_time('csc,dt', cut_csc_spread=True, cut_rpc_hits=True)
    ms_r3.cut_time('csc,dt', cut_csc_spread=True, cut_rpc_hits=True)
    print(f'   IT | MC : {ms_mc.count():>10,} | R3 : {ms_r3.count():>10,}')

    ms_mc.f((ms_mc['nCsc']==1) & (ms_mc['nDt']==1))
    ms_r3.f((ms_r3['nCsc']==1) & (ms_r3['nDt']==1))
    print(f'CSCDT | MC : {ms_mc.count():>10,} | R3 : {ms_r3.count():>10,}')

    ms_mc.tag(tags='cscdt')
    ms_r3.tag(tags='cscdt')
    ms_r3.blind('dphi')
    print(f' dPhi | MC : {ms_mc.count():>10,} | R3 : {ms_r3.count():>10,}')
    print('')

    #!!! TURNING CUTS ON !!!!#
    ms_mc.cut, ms_r3.cut = True, True
    #!!!!!!!!!!!!!!!!!!!!!!!!#

    # ms_mc = MuonSystemRDF(ff_mc, isMC=True)  # , nev=N_EVENTS)
    # ms_r3 = MuonSystemRDF(ff_r3, isMC=False)  # , nev=N_EVENTS)

    # rdf_mc = RDataFrame('MuonSystem', ff_mc)
    # rdf_r3 = RDataFrame('MuonSystem', ff_r3)
    # print(rdf_mc.GetColumnNames())

    #! Need to do
    #! weight[k] = T.array('weight')[sel_ev[k]] * lumi
    #! if 'data' in k:weight[k] = T.array('weight')[sel_ev[k]]*0.0 + 1

    # Apply cuts:
    #   - Require HLT
    #   - Require L1 plateau
    #   - If MC: Require CSC clusters to be in CSC
    #   - If MC: Require DT clusters to be in DT
    #   - Require 1 CSC & 1 DT
    #   - If BKG: Requre dPhi > 0.5

    # print(f'Events In | {ms_mc.Count()=:>9,} | {ms_r3.Count()=:>10,}')

    # ms_mc.Filter('(HLTDecision[566] > 0) || (HLTDecision[569] > 0) || (HLTDecision[570] > 0)', implicit=True)
    # ms_mc.L1_plateau(implicit=True) #TODO: Need to fix
    # ms_mc.match_in_det(implicit=True) # Sig
    # ms_mc.time_cut('it', system='cscdt', implicit=True)
    # ms_mc.Filter('dtRechitClusterNHitStation1 == 0', 'dt', implicit=True)
    # ms_mc.define_2tag_kins_and_cut(system='cscdt', implicit=True)
    # ms_mc.Define('weight', 'weight * 23.02 * 1000', implicit=True)

    # ms_r3.Filter('(HLTDecision[566] > 0) || (HLTDecision[569] > 0) || (HLTDecision[570] > 0)', implicit=True)
    # ms_r3.L1_plateau(implicit=True) #TODO: Need to fix
    # ms_r3.time_cut('it', system='cscdt', implicit=True)
    # ms_r3.Filter('dtRechitClusterNHitStation1 == 0', 'dt', implicit=True)
    # ms_r3.define_2tag_kins_and_cut(system='cscdt', implicit=True)
    # ms_r3.Filter('tag_dPhi > 0.5', implicit=True) # Bkg
    # ms_mc.Define('weight', '1', implicit=True)

    # print(f'Events Post Cuts | {ms_mc.Count()=:>9,} | {ms_r3.Count()=:>10,}')


    # Convert to Numpy array
    vals = 'tag_dR,tag_dEta,tag_dPhi'
    vals += ',cscNStation,cscAvgStation,cscSpreadR'



    # print('')
    # print(f'| ----- | {"-"*9} | {"-"*10} |')
    # print(f'| Step  | {"MC":^9} | {"R3":^10} |')
    # print(f'| ----- | {"-"*9} | {"-"*10} |')
    # print(f'| Raw   | {ms_mc.Count():>9,} | {ms_r3.Count():>10,} |')
    # print(f'| Match | {ms_mc.Count():>9,} | {ms_r3.Count():>10,} |')
    # print(f'| HLT   | {ms_mc.Count():>9,} | {ms_r3.Count():>10,} |')
    # print(f'| L1    | {ms_mc.Count():>9,} | {ms_r3.Count():>10,} |')
    # print(f'| OOT   | {ms_mc.Count():>9,} | {ms_r3.Count():>10,} |')
    # print(f'| 2tag  | {ms_mc.Count():>9,} | {ms_r3.Count():>10,} |')
    # print(f'| ----- | {"-"*9} | {"-"*10} |')

    print('+-----------------------+')
    print('| Finished train_bdt.py |')
    print('+-----------------------+')
