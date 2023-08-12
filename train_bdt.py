"""train_bdt.py"""
import os
import sys
import pathlib
from collections import defaultdict

import numpy as np
import sklearn as skl

import ROOT as rt
# from ROOT import gErrorIgnoreLevel

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
    rt.gROOT.SetBatch()
    # rt.EnableImplicitMT(N_THREADS)
    # alert(f'Setting {N_THREADS=:,} and {rt.GetThreadPoolSize()=:,}')
    # print(f'Running with thread pool size = {rt.GetThreadPoolSize():,}')

    #############################

    if len(sys.argv) > 1:
        N_EVENTS = int(sys.argv[1])
        N_EVENTS = N_EVENTS if N_EVENTS > 0 else None
        alert(f'Setting {N_EVENTS=}')

    #############################

    ff_mc = 'data/raw/ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v6.root'
    ff_r3 = 'data/raw/DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi_v6.root'


    ms_mc = MuonSystemAwkward(ff_mc, nev=N_EVENTS, is_mc=True)
    ms_r3 = MuonSystemAwkward(ff_r3, nev=N_EVENTS, is_mc=False)

    # ms_mc = MuonSystemRDF(ff_mc, isMC=True)  # , nev=N_EVENTS)
    # ms_r3 = MuonSystemRDF(ff_r3, isMC=False)  # , nev=N_EVENTS)

    # rdf_mc = rt.RDataFrame('MuonSystem', ff_mc)
    # rdf_r3 = rt.RDataFrame('MuonSystem', ff_r3)
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
