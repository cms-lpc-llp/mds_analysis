import sys
import pathlib
import numpy as np

import pickle

import ROOT as rt
from ROOT import RDataFrame

from src.histo_utilities import std_color_list as SCL

# **************************** #
LOCAL_DIR = "/home/psimmerl/mds_analysis"
OUT_DIR = f"{LOCAL_DIR}/reports/weekly/2024-04-15"

# STAT = 'raw'
# LUMI = 23.02 * 1000
# FN_MC = f'{LOCAL_DIR}/data/raw/ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v6.root'
# FN_R3 = f'{LOCAL_DIR}/data/raw/DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi_v6.root'

STAT = "pedro"
LUMI = 1.1328524540090597e-06 * 23.02 * 1000  # ???
RUN2_BR = 2.16e-03  # Adjust weights to Run 2 BR limit
FN_MC = f"{LOCAL_DIR}/data/raw/mc_pedro.root"
# FN_R3 = f'{LOCAL_DIR}/data/raw/data_pedro.root'
# FN_MC = f'{LOCAL_DIR}/data/processed/mc_pedro_hlt569.root'
FN_R3 = f"{LOCAL_DIR}/data/processed/r3_pedro_hlt569.root"


# **** #
LOW_MET_CUTOFF = 75
HIGH_MET_CUTOFF = 150

# **** #
CUTS_L1 = [
    "acceptance",
    "HLT",
    "L1",
    "MET",
    #! reset cutflow indices here
    "CSC IT",
    "DT IT",
    "MB1",
    "1 CSC-DT",
    # 'dPhi',
]

CUTS = [
    "acceptance",
    "HLT",
    "L1",
    "MET",
    #! reset cutflow indices here
    "CSC IT",
    "DT IT",
    # 'ME1',
    "MB1",
    # 'n jets',
    # 'n leptons',
    "jet veto",
    "muon veto",
    "halo veto",
    # 'DT N Station',
    # 'DT Avg Station',
    # 'DT stn',
    # 'BDT',
    "DNN",
    "1 CSC-DT",
    # 'dR',
    # 'dEta',
    "dPhi",
    # 'ABCD'
    # 'DT size' # for ABCD?
]

# **** #
CUT_VALUES = {
    "l1_lt200": {
        "MIN_CSC_TIME": -5.0,
        "MAX_CSC_TIME": 12.5,
        "MAX_CSC_TSPREAD": 20.0,
        "MAX_RPC_BX": 0,
        "MIN_RPC_HITS": 1,
        "MAX_CSC_JET": 200,
        "MAX_DT_JET": 200,
        "MAX_CSC_MUON": 200,
        "MAX_DT_MUON": 200,
        "MAX_ME1": 0,
        "MAX_MB1": 10,
        "HALO_CUTOFF": 0,
        "MIN_DPHI": 0.4,
        "MIN_DETA": 0,
        "MAX_DETA": 4,
        "MIN_CSC_DNN": 0,
        # "MIN_DT_DNN": 0,
    },
    "scs_lt200": {
        "MIN_CSC_TIME": -5.0,
        "MAX_CSC_TIME": 12.5,
        "MAX_CSC_TSPREAD": 20.0,
        "MAX_RPC_BX": 0,
        "MIN_RPC_HITS": 1,
        "MAX_CSC_JET": 30,
        "MAX_DT_JET": 50,
        "MAX_CSC_MUON": 30,
        "MAX_DT_MUON": 10,
        "MAX_ME1": 0,
        "MAX_MB1": 0,
        "HALO_CUTOFF": 0.4,
        "MIN_DPHI": 0.4,
        "MIN_DETA": 0,
        "MAX_DETA": 4,
        "MIN_CSC_DNN": 0,
        # "MIN_DT_DNN": 0,
    },
    "scs_low": {
        "MIN_CSC_TIME": -5.0,
        "MAX_CSC_TIME": 12.5,
        "MAX_CSC_TSPREAD": 20.0,
        "MAX_RPC_BX": 0,
        "MIN_RPC_HITS": 1,
        "MAX_CSC_JET": 30,
        "MAX_DT_JET": 50,
        "MAX_CSC_MUON": 30,
        "MAX_DT_MUON": 10,
        "MAX_ME1": 0,
        "MAX_MB1": 0,
        "HALO_CUTOFF": 0.4,
        "MIN_DPHI": 0.4,
        "MIN_DETA": 0,
        "MAX_DETA": 4,
        "MIN_CSC_DNN": 0,
        # "MIN_DT_DNN": 0,
    },
    "scs_high": {
        "MIN_CSC_TIME": -5.0,
        "MAX_CSC_TIME": 12.5,
        "MAX_CSC_TSPREAD": 20.0,
        "MAX_RPC_BX": 0,
        "MIN_RPC_HITS": 1,
        "MAX_CSC_JET": 30,
        "MAX_DT_JET": 50,
        "MAX_CSC_MUON": 30,
        "MAX_DT_MUON": 10,
        "MAX_ME1": 0,
        "MAX_MB1": 5,
        "HALO_CUTOFF": 0.4,
        "MIN_DPHI": 0.4,
        "MIN_DETA": 0,
        "MAX_DETA": 4,
        "MIN_CSC_DNN": 0,
        # "MIN_DT_DNN": 0,
    },
    "lopt_lt200": {
        # Default guess
        "MIN_CSC_TIME": -5.0,
        "MAX_CSC_TIME": 12.5,
        "MAX_CSC_TSPREAD": 20.0,
        "MAX_RPC_BX": 0,
        "MIN_RPC_HITS": 1,
        # 'MAX_N_JETS' : 15,
        # 'MAX_N_LEPS' : 5,
        "MAX_CSC_JET": 100,
        "MAX_DT_JET": 100,
        "MAX_CSC_MUON": 100,
        "MAX_DT_MUON": 100,
        "MAX_ME1": 0,
        "MAX_MB1": 10,
        "HALO_CUTOFF": 0.2,
        "MIN_DPHI": 0.4,
        "MIN_DETA": 0,
        "MAX_DETA": 4,
        "MIN_CSC_DNN": 0,
        # "MIN_DT_DNN": 0,
    },
    "lopt_low": {
        # Default guess
        "MIN_CSC_TIME": -5.0,
        "MAX_CSC_TIME": 12.5,
        "MAX_CSC_TSPREAD": 20.0,
        "MAX_RPC_BX": 0,
        "MIN_RPC_HITS": 1,
        # 'MAX_N_JETS' : 15,
        # 'MAX_N_LEPS' : 5,
        "MAX_CSC_JET": 100,
        "MAX_DT_JET": 100,
        "MAX_CSC_MUON": 100,
        "MAX_DT_MUON": 100,
        "MAX_ME1": 0,
        "MAX_MB1": 10,
        "HALO_CUTOFF": 0.2,
        "MIN_DPHI": 0.4,
        "MIN_DETA": 0,
        "MAX_DETA": 4,
        "MIN_CSC_DNN": 0,
        # "MIN_DT_DNN": 0,
        # #
        # 'MIN_CSC_TIME' : -5.00,
        # 'MAX_CSC_TIME' : 12.50,
        # 'MAX_CSC_TSPREAD' : 20.00,
        # 'MAX_RPC_BX' : 0.00,
        # 'MIN_RPC_HITS' : 1.00,
        # # 'MAX_N_JETS' : 15,
        # # 'MAX_N_LEPS' : 5,
        # 'MAX_CSC_JET' : 141.81,
        # 'MAX_DT_JET' : 10.01,
        # 'MAX_CSC_MUON' : 200.00,
        # 'MAX_DT_MUON' : 200.00,
        # 'MAX_ME1' : 0.00,
        # 'MAX_MB1' : 0.98,
        # 'HALO_CUTOFF' : 0.00,
        # 'MIN_DPHI' : 0.40,
        # 'MIN_DETA' : 0.00,
        # 'MAX_DETA' : 4.00,
        # 'MIN_CSC_DNN' : 0.00,
        # 'MIN_DT_DNN' : 0.00,
    },
    "lopt_high": {
        # Default guess
        "MIN_CSC_TIME": -5.0,
        "MAX_CSC_TIME": 12.5,
        "MAX_CSC_TSPREAD": 20.0,
        "MAX_RPC_BX": 0,
        "MIN_RPC_HITS": 1,
        # 'MAX_N_JETS' : 15,
        # 'MAX_N_LEPS' : 5,
        "MAX_CSC_JET": 100,
        "MAX_DT_JET": 100,
        "MAX_CSC_MUON": 100,
        "MAX_DT_MUON": 100,
        "MAX_ME1": 0,
        "MAX_MB1": 10,
        "HALO_CUTOFF": 0.2,
        "MIN_DPHI": 0.4,
        "MIN_DETA": 0,
        "MAX_DETA": 4,
        "MIN_CSC_DNN": 0,
        # "MIN_DT_DNN": 0,
        # #
        # 'MIN_CSC_TIME' : -5.0,
        # 'MAX_CSC_TIME' : 12.5,
        # 'MAX_CSC_TSPREAD' : 20.0,
        # 'MAX_RPC_BX' : 0,
        # 'MIN_RPC_HITS' : 1,
        # # 'MAX_N_JETS' : 15,
        # # 'MAX_N_LEPS' : 5,
        # 'MAX_CSC_JET' : 100,
        # 'MAX_DT_JET' : 100,
        # 'MAX_CSC_MUON' : 100,
        # 'MAX_DT_MUON' : 100,
        # 'MAX_ME1' : 0,
        # 'MAX_MB1' : 1,
        # # 10: S = 45.9, B = 9, S/sqrt[B] = 15.4
        # #  9: S = 45.3, B = 9, S/sqrt[B] = 15.1
        # #  8: S = 44.7, B = 9, S/sqrt[B] = 14.9
        # #  7: S = 44.7, B = 8, S/sqrt[B] = 15.8
        # #  6: S = 44.2, B = 8, S/sqrt[B] = 15.6
        # #  5: S = 43.5, B = 6, S/sqrt[B] = 17.8
        # #  4: S = 42.9, B = 6, S/sqrt[B] = 17.5
        # #  3: S = 40.5, B = 6, S/sqrt[B] = 16.5
        # #  2: S = 39.9, B = 5, S/sqrt[B] = 17.8
        # #  1: S = 36.9, B = 3, S/sqrt[B] = 21.3
        # #  0: S = 33.8, B = 2, S/sqrt[B] = 23.9
        # 'HALO_CUTOFF' : 0.2,
        # 'MIN_DPHI' : 0.4,
        # 'MIN_DETA' : 0,
        # 'MAX_DETA' : 4,
        # 'MIN_CSC_DNN' : 0,
        # 'MIN_DT_DNN' : 0,
    },
    "loptDNN_lt200": {
        # Default guess
        "MIN_CSC_TIME": -5.0,
        "MAX_CSC_TIME": 12.5,
        "MAX_CSC_TSPREAD": 20.0,
        "MAX_RPC_BX": 0,
        "MIN_RPC_HITS": 1,
        # 'MAX_N_JETS' : 15,
        # 'MAX_N_LEPS' : 5,
        "MAX_CSC_JET": 100,
        "MAX_DT_JET": 100,
        "MAX_CSC_MUON": 100,
        "MAX_DT_MUON": 100,
        "MAX_ME1": 0,
        "MAX_MB1": 10,
        "HALO_CUTOFF": 0.2,
        "MIN_DPHI": 0.4,
        "MIN_DETA": 0,
        "MAX_DETA": 4,
        "MIN_CSC_DNN": 0,
        # "MIN_DT_DNN": 0,
    },
    "loptDNN_low": {
        # Default guess
        "MIN_CSC_TIME": -5.0,
        "MAX_CSC_TIME": 12.5,
        "MAX_CSC_TSPREAD": 20.0,
        "MAX_RPC_BX": 0,
        "MIN_RPC_HITS": 1,
        # 'MAX_N_JETS' : 15,
        # 'MAX_N_LEPS' : 5,
        "MAX_CSC_JET": 100,
        "MAX_DT_JET": 100,
        "MAX_CSC_MUON": 100,
        "MAX_DT_MUON": 100,
        "MAX_ME1": 0,
        "MAX_MB1": 10,
        "HALO_CUTOFF": 0.2,
        "MIN_DPHI": 0.4,
        "MIN_DETA": 0,
        "MAX_DETA": 4,
        "MIN_CSC_DNN": 0,
        # "MIN_DT_DNN": 0,
        # # Optimized with s2b
        # 'MIN_CSC_TIME' : -5.00,
        # 'MAX_CSC_TIME' : 12.50,
        # 'MAX_CSC_TSPREAD' : 20.00,
        # 'MAX_RPC_BX' : 0.00,
        # 'MIN_RPC_HITS' : 1.00,
        # 'MAX_CSC_JET' : 112.38,
        # 'MAX_DT_JET' : 61.94,
        # 'MAX_CSC_MUON' : 200.00,
        # 'MAX_DT_MUON' : 198.72,
        # 'MAX_ME1' : 0.00,
        # 'MAX_MB1' : 0.36,
        # 'HALO_CUTOFF' : 0.00,
        # 'MIN_DPHI' : 0.40,
        # 'MIN_DETA' : 0.00,
        # 'MAX_DETA' : 4.00,
        # 'MIN_CSC_DNN' : 0.69,
        # 'MIN_DT_DNN' : 0.79,
    },
    "loptDNN_high": {
        # Default guess
        "MIN_CSC_TIME": -5.0,
        "MAX_CSC_TIME": 12.5,
        "MAX_CSC_TSPREAD": 20.0,
        "MAX_RPC_BX": 0,
        "MIN_RPC_HITS": 1,
        # 'MAX_N_JETS' : 15,
        # 'MAX_N_LEPS' : 5,
        "MAX_CSC_JET": 100,
        "MAX_DT_JET": 100,
        "MAX_CSC_MUON": 100,
        "MAX_DT_MUON": 100,
        "MAX_ME1": 0,
        "MAX_MB1": 10,
        "HALO_CUTOFF": 0.2,
        "MIN_DPHI": 0.4,
        "MIN_DETA": 0,
        "MAX_DETA": 4,
        "MIN_CSC_DNN": 0,
        # "MIN_DT_DNN": 0,
        # #
        # 'MIN_CSC_TIME' : -5.00,
        # 'MAX_CSC_TIME' : 12.50,
        # 'MAX_CSC_TSPREAD' : 20.00,
        # 'MAX_RPC_BX' : 0.00,
        # 'MIN_RPC_HITS' : 1.00,
        # 'MAX_CSC_JET' : 100.00,
        # 'MAX_DT_JET' : 100.00,
        # 'MAX_CSC_MUON' : 100.00,
        # 'MAX_DT_MUON' : 100.00,
        # 'MAX_ME1' : 0.00,
        # 'MAX_MB1' : 10.00,
        # 'HALO_CUTOFF' : 0.20,
        # 'MIN_DPHI' : 0.40,
        # 'MIN_DETA' : 0.00,
        # 'MAX_DETA' : 4.00,
        # 'MIN_CSC_DNN' : 0.35,
        # 'MIN_DT_DNN' : 0.00,
    },
    "ropt_lt200": {
        # Default guess
        "MIN_CSC_TIME": -5.0,
        "MAX_CSC_TIME": 12.5,
        "MAX_CSC_TSPREAD": 20.0,
        "MAX_RPC_BX": 0,
        "MIN_RPC_HITS": 1,
        # 'MAX_N_JETS' : 15,
        # 'MAX_N_LEPS' : 5,
        "MAX_CSC_JET": 100,
        "MAX_DT_JET": 100,
        "MAX_CSC_MUON": 100,
        "MAX_DT_MUON": 100,
        "MAX_ME1": 0,
        "MAX_MB1": 10,
        "HALO_CUTOFF": 0.2,
        "MIN_DPHI": 0.4,
        "MIN_DETA": 0,
        "MAX_DETA": 4,
        "MIN_CSC_DNN": 0,
        # "MIN_DT_DNN": 0,
        #
    },
    "ropt_low": {
        # # Default guess
        # "MIN_CSC_TIME": -5.0,
        # "MAX_CSC_TIME": 12.5,
        # "MAX_CSC_TSPREAD": 20.0,
        # "MAX_RPC_BX": 0,
        # "MIN_RPC_HITS": 1,
        # # 'MAX_N_JETS' : 15,
        # # 'MAX_N_LEPS' : 5,
        # "MAX_CSC_JET": 100,
        # "MAX_DT_JET": 100,
        # "MAX_CSC_MUON": 100,
        # "MAX_DT_MUON": 100,
        # "MAX_ME1": 0,
        # "MAX_MB1": 10,
        # "HALO_CUTOFF": 0.2,
        # "MIN_DPHI": 0.4,
        # "MIN_DETA": 0,
        # "MAX_DETA": 4,
        # "MIN_CSC_DNN": 0,
        # "MIN_DT_DNN": 0,
        # Optimized with s2b12 and IT data
        # "MIN_CSC_TIME": -5.00,
        # "MAX_CSC_TIME": 12.50,
        # "MAX_CSC_TSPREAD": 20.00,
        # "MAX_RPC_BX": 0,
        # "MIN_RPC_HITS": 1,
        # "MAX_CSC_JET": 200,
        # "MAX_DT_JET": 10,
        # "MAX_CSC_MUON": 171,
        # "MAX_DT_MUON": 198,
        # "MAX_ME1": 0,
        # "MAX_MB1": 0,
        # "HALO_CUTOFF": 0.00,
        # "MIN_DPHI": 0.40,
        # "MIN_DETA": 0.0,
        # "MAX_DETA": 4.0,
        # "MIN_CSC_DNN": 0.00,
        # "MIN_DT_DNN": 0.00,
        # Optimized with s1 and IT data
        "MIN_CSC_TIME": -5.00,
        "MAX_CSC_TIME": 12.50,
        "MAX_CSC_TSPREAD": 20.00,
        "MAX_RPC_BX": 0,
        "MIN_RPC_HITS": 1,
        "MAX_CSC_JET": 10,
        "MAX_DT_JET": 10,
        "MAX_CSC_MUON": 167,
        "MAX_DT_MUON": 199,
        "MAX_ME1": 0,
        "MAX_MB1": 0,
        "HALO_CUTOFF": 0.00,
        "MIN_DPHI": 0.40,
        "MIN_DETA": 0.0,
        "MAX_DETA": 4.0,
        "MIN_CSC_DNN": 0.00,
        # "MIN_DT_DNN": 0.00,
    },
    "ropt_high": {
        # # Default guess
        # "MIN_CSC_TIME": -5.0,
        # "MAX_CSC_TIME": 12.5,
        # "MAX_CSC_TSPREAD": 20.0,
        # "MAX_RPC_BX": 0,
        # "MIN_RPC_HITS": 1,
        # # 'MAX_N_JETS' : 15,
        # # 'MAX_N_LEPS' : 5,
        # "MAX_CSC_JET": 100,
        # "MAX_DT_JET": 100,
        # "MAX_CSC_MUON": 100,
        # "MAX_DT_MUON": 100,
        # "MAX_ME1": 0,
        # "MAX_MB1": 10,
        # "HALO_CUTOFF": 0.2,
        # "MIN_DPHI": 0.4,
        # "MIN_DETA": 0,
        # "MAX_DETA": 4,
        # "MIN_CSC_DNN": 0,
        # "MIN_DT_DNN": 0,
        # Found manually
        "MIN_CSC_TIME": -5.0,
        "MAX_CSC_TIME": 12.5,
        "MAX_CSC_TSPREAD": 20.0,
        "MAX_RPC_BX": 0,
        "MIN_RPC_HITS": 1,
        # 'MAX_N_JETS' : 15,
        # 'MAX_N_LEPS' : 5,
        "MAX_CSC_JET": 999,
        "MAX_DT_JET": 999,
        "MAX_CSC_MUON": 999,
        "MAX_DT_MUON": 999,
        "MAX_ME1": 0,
        "MAX_MB1": 10,
        "HALO_CUTOFF": 0.0,
        "MIN_DPHI": 0.4,
        "MIN_DETA": 0,
        "MAX_DETA": 4,
        "MIN_CSC_DNN": 0,
        # "MIN_DT_DNN": 0,
    },
    "roptDNN_lt200": {
        # Default guess
        "MIN_CSC_TIME": -5.0,
        "MAX_CSC_TIME": 12.5,
        "MAX_CSC_TSPREAD": 20.0,
        "MAX_RPC_BX": 0,
        "MIN_RPC_HITS": 1,
        # 'MAX_N_JETS' : 15,
        # 'MAX_N_LEPS' : 5,
        "MAX_CSC_JET": 100,
        "MAX_DT_JET": 100,
        "MAX_CSC_MUON": 100,
        "MAX_DT_MUON": 100,
        "MAX_ME1": 0,
        "MAX_MB1": 10,
        "HALO_CUTOFF": 0.2,
        "MIN_DPHI": 0.4,
        "MIN_DETA": 0,
        "MAX_DETA": 4,
        "MIN_CSC_DNN": 0,
        # "MIN_DT_DNN": 0,
        #
    },
    "roptDNN_low": {
        # # Default guess
        # "MIN_CSC_TIME": -5.0,
        # "MAX_CSC_TIME": 12.5,
        # "MAX_CSC_TSPREAD": 20.0,
        # "MAX_RPC_BX": 0,
        # "MIN_RPC_HITS": 1,
        # # 'MAX_N_JETS' : 15,
        # # 'MAX_N_LEPS' : 5,
        # "MAX_CSC_JET": 100,
        # "MAX_DT_JET": 100,
        # "MAX_CSC_MUON": 100,
        # "MAX_DT_MUON": 100,
        # "MAX_ME1": 0,
        # "MAX_MB1": 10,
        # "HALO_CUTOFF": 0.2,
        # "MIN_DPHI": 0.4,
        # "MIN_DETA": 0,
        # "MAX_DETA": 4,
        # "MIN_CSC_DNN": 0,
        # # "MIN_DT_DNN": 0,
        # From ropt_low with s1 opt
        "MIN_CSC_TIME": -5.00,
        "MAX_CSC_TIME": 12.50,
        "MAX_CSC_TSPREAD": 20.00,
        "MAX_RPC_BX": 0,
        "MIN_RPC_HITS": 1,
        "MAX_CSC_JET": 10,
        "MAX_DT_JET": 10,
        "MAX_CSC_MUON": 167,
        "MAX_DT_MUON": 199,
        "MAX_ME1": 0,
        "MAX_MB1": 0,
        "HALO_CUTOFF": 0.00,
        "MIN_DPHI": 0.40,
        "MIN_DETA": 0.0,
        "MAX_DETA": 4.0,
        "MIN_CSC_DNN": 0.00,
        # "MIN_CSC_DNN": 0.917,  # bkgMC
        # "MIN_CSC_DNN": 0.97,  # bkgMC_plusBeamHalo
        # "MIN_CSC_DNN": 0.97,  # bkgOOTData # or 0.89 or 0.94 or 0.96
        # "MIN_DT_DNN": 0.00,
    },
    "roptDNN_high": {
        # Default guess
        "MIN_CSC_TIME": -5.0,
        "MAX_CSC_TIME": 12.5,
        "MAX_CSC_TSPREAD": 20.0,
        "MAX_RPC_BX": 0,
        "MIN_RPC_HITS": 1,
        # 'MAX_N_JETS' : 15,
        # 'MAX_N_LEPS' : 5,
        "MAX_CSC_JET": 100,
        "MAX_DT_JET": 100,
        "MAX_CSC_MUON": 100,
        "MAX_DT_MUON": 100,
        "MAX_ME1": 0,
        "MAX_MB1": 10,
        "HALO_CUTOFF": 0.2,
        "MIN_DPHI": 0.4,
        "MIN_DETA": 0,
        "MAX_DETA": 4,
        "MIN_CSC_DNN": 0,
        # "MIN_DT_DNN": 0,
        #
    },
}

# OPT_CUTS = ['MB1','CSC DNN','DT DNN']
OPT_CUTS = [
    "MB1",
    "CSC jet veto",
    "DT jet veto",
    "CSC muon veto",
    "DT muon veto",
    "halo veto",
    "CSC DNN",
    # "DT DNN",
]  # ,'DT size','min dPhi','min dEta','max dEta','ABCD phi', 'ABCD size',]#,'ME1']
CUT_OPT_PARS = {
    "ME1": {
        "col": "cscNHitME1",
        "cut": "MAX_ME1",
        "values": np.arange(0, 11, 1)[::-1],
        "func": lambda v, vd: vd <= v,
    },
    "MB1": {
        "col": "dtNHitStation1",
        "cut": "MAX_MB1",
        "values": np.arange(0, 11, 1)[::-1],
        "func": lambda v, vd: vd <= v,
    },
    "CSC jet veto": {
        "col": "cscJetVetoPt",
        "cut": "MAX_CSC_JET",
        "values": np.arange(10, 201, 1)[::-1],
        # 'func' : lambda v, vd: (10 < vd) & (vd < v),
        "func": lambda v, vd: vd < v,
    },
    "DT jet veto": {
        "col": "dtJetVetoPt",
        "cut": "MAX_DT_JET",
        "values": np.arange(10, 201, 1)[::-1],
        # 'func' : lambda v, vd: (10 < vd) & (vd < v),
        "func": lambda v, vd: vd < v,
    },
    "CSC muon veto": {
        "col": "cscMuonVetoPt",
        "cut": "MAX_CSC_MUON",
        "values": np.arange(0, 201, 1)[::-1],
        # 'func' : lambda v, vd: (0 < vd) & (vd < v),
        "func": lambda v, vd: vd < v,
    },
    "DT muon veto": {
        "col": "dtMuonVetoPt",
        "cut": "MAX_DT_MUON",
        "values": np.arange(0, 201, 1)[::-1],
        # 'func' : lambda v, vd: (0 < vd) & (vd < v),
        "func": lambda v, vd: vd < v,
    },
    "halo veto": {
        "col": "dtPhi",
        "cut": "HALO_CUTOFF",
        "values": np.arange(0, 0.52, 0.02),
        "func": lambda v, vd: (v < np.abs(vd)) & (np.abs(vd) < PI - v),
    },
    "min dPhi": {
        "col": "tag_dPhi",
        "cut": "MIN_DPHI",
        "values": np.arange(0, 0.62, 0.02),
        "func": lambda v, vd: v < vd,
    },
    "min dEta": {
        "col": "tag_dEta",
        "cut": "MIN_DETA",
        "values": np.arange(0, 4.1, 0.1),
        "func": lambda v, vd: v < vd,
    },
    "max dEta": {
        "col": "tag_dEta",
        "cut": "MAX_DETA",
        "values": np.arange(0, 4.1, 0.1)[::-1],
        "func": lambda v, vd: vd < v,
    },
    "CSC DNN": {
        "col": "cscDNN",
        "cut": "MIN_CSC_DNN",
        "values": np.arange(0, 1.01, 0.01),
        "func": lambda v, vd: v < vd,
    },
    # "DT DNN": {
    #     "col": "dtDNN",
    #     "cut": "MIN_DT_DNN",
    #     "values": np.arange(0, 1.01, 0.01),
    #     "func": lambda v, vd: v < vd,
    # },
    "ABCD phi": {
        "col": "tag_dPhi",
        "cut": "ABCD_PHI",
        "values": np.arange(0, 3.5, 0.01),
        "func": lambda v, vd: v <= vd,
    },
    "ABCD size": {
        "col": "dtSize",
        "cut": "ABCD_SIZE",
        "values": np.arange(50, 150, 1),
        "func": lambda v, vd: v <= vd,
    },
    "DT size": {
        "col": "dtSize",
        "cut": "MIN_DT_SIZE",
        "values": np.arange(50, 150, 1),
        "func": lambda v, vd: v <= vd,
    },
}

# **** #
C, D = "cscRechitCluster", "dtRechitCluster"

COLUMNS_OUT = [
    # 'Flag2_BadChargedCandidateFilter',
    # 'Flag2_BadPFMuonFilter',
    # 'Flag2_EcalDeadCellTriggerPrimitiveFilter',
    # 'Flag2_HBHEIsoNoiseFilter',
    # 'Flag2_HBHENoiseFilter',
    # 'Flag2_all',
    # 'Flag2_ecalBadCalibFilter',
    # 'Flag2_eeBadScFilter',
    # 'Flag2_globalSuperTightHalo2016Filter',
    # 'Flag2_globalTightHalo2016Filter',
    # 'Flag_BadChargedCandidateFilter',
    # 'Flag_BadPFMuonFilter',
    # 'Flag_CSCTightHaloFilter',
    # 'Flag_HBHEIsoNoiseFilter',
    # 'Flag_HBHENoiseFilter',
    # 'Flag_all',
    # 'Flag_ecalBadCalibFilter',
    # 'Flag_eeBadScFilter',
    # 'Flag_globalSuperTightHalo2016Filter',
    # 'Flag_goodVertices',
    # 'HLTDecision',
    # 'MC_condition',
    #
    "cscRechitClusterAvgStation10",
    "cscRechitClusterEta",
    "cscRechitClusterJetVetoE",
    "cscRechitClusterJetVetoLooseId",
    "cscRechitClusterJetVetoPt",
    "cscRechitClusterJetVetoTightId",
    "cscRechitClusterMaxChamber",
    "cscRechitClusterMaxChamberRatio",
    "cscRechitClusterMaxStation",
    "cscRechitClusterMaxStationRatio",
    "cscRechitClusterMet_dPhi",
    "cscRechitClusterMuonVetoE",
    "cscRechitClusterMuonVetoGlobal",
    "cscRechitClusterMuonVetoLooseId",
    "cscRechitClusterMuonVetoPt",
    "cscRechitClusterNChamber",
    "cscRechitClusterNRechitChamberMinus11",
    "cscRechitClusterNRechitChamberMinus12",
    "cscRechitClusterNRechitChamberMinus13",
    "cscRechitClusterNRechitChamberMinus21",
    "cscRechitClusterNRechitChamberMinus22",
    "cscRechitClusterNRechitChamberMinus31",
    "cscRechitClusterNRechitChamberMinus32",
    "cscRechitClusterNRechitChamberMinus41",
    "cscRechitClusterNRechitChamberMinus42",
    "cscRechitClusterNRechitChamberPlus11",
    "cscRechitClusterNRechitChamberPlus12",
    "cscRechitClusterNRechitChamberPlus13",
    "cscRechitClusterNRechitChamberPlus21",
    "cscRechitClusterNRechitChamberPlus22",
    "cscRechitClusterNRechitChamberPlus31",
    "cscRechitClusterNRechitChamberPlus32",
    "cscRechitClusterNRechitChamberPlus41",
    "cscRechitClusterNRechitChamberPlus42",
    "cscRechitClusterNStation10",
    "cscRechitClusterPhi",
    "cscRechitClusterSize",
    "cscRechitClusterTime",
    "cscRechitClusterTimeSpread",
    "cscRechitClusterTimeSpreadWeightedAll",
    "cscRechitClusterTimeWeighted",
    "cscRechitClusterX",
    "cscRechitClusterY",
    "cscRechitClusterZ",
    "cscRechitCluster_match_MB1Seg_0p4",
    "cscRechitCluster_match_RB1_0p4",
    "cscRechitCluster_match_RE12_0p4",
    "cscRechitCluster_match_dtSeg_0p4",
    # 'cscRechitCluster_match_gLLP',
    # 'cscRechitCluster_match_gLLP_csc',
    # 'cscRechitCluster_match_gLLP_decay_r',
    # 'cscRechitCluster_match_gLLP_decay_z',
    # 'cscRechitCluster_match_gLLP_dt',
    # 'cscRechitCluster_match_gLLP_e',
    # 'cscRechitCluster_match_gLLP_eta',
    # 'cscRechitCluster_match_gLLP_index',
    # 'cscRechitCluster_match_gLLP_minDeltaR',
    # 'cscRechitCluster_match_gLLP_phi',
    #
    # 'ctau',
    #
    "dtRechitClusterAvgStation10",
    "dtRechitClusterEta",
    "dtRechitClusterJetVetoE",
    "dtRechitClusterJetVetoLooseId",
    "dtRechitClusterJetVetoPt",
    "dtRechitClusterJetVetoTightId",
    "dtRechitClusterMaxChamber",
    "dtRechitClusterMaxChamberRatio",
    "dtRechitClusterMaxStation",
    "dtRechitClusterMaxStationRatio",
    "dtRechitClusterMet_dPhi",
    "dtRechitClusterMuonVetoE",
    "dtRechitClusterMuonVetoGlobal",
    "dtRechitClusterMuonVetoLooseId",
    "dtRechitClusterMuonVetoPt",
    "dtRechitClusterMuonVetoTightId",
    "dtRechitClusterNChamber",
    "dtRechitClusterNHitStation1",
    "dtRechitClusterNHitStation2",
    "dtRechitClusterNHitStation3",
    "dtRechitClusterNHitStation4",
    "dtRechitClusterNOppositeSegStation1",
    "dtRechitClusterNOppositeSegStation2",
    "dtRechitClusterNOppositeSegStation3",
    "dtRechitClusterNOppositeSegStation4",
    "dtRechitClusterNSegStation1",
    "dtRechitClusterNSegStation2",
    "dtRechitClusterNSegStation3",
    "dtRechitClusterNSegStation4",
    "dtRechitClusterNStation10",
    "dtRechitClusterNoiseHit",
    "dtRechitClusterNoiseHitStation1",
    "dtRechitClusterNoiseHitStation2",
    "dtRechitClusterNoiseHitStation3",
    "dtRechitClusterNoiseHitStation4",
    "dtRechitClusterOverlap",
    "dtRechitClusterPhi",
    "dtRechitClusterSize",
    "dtRechitClusterWheel",
    "dtRechitClusterX",
    "dtRechitClusterY",
    "dtRechitClusterZ",
    "dtRechitCluster_match_MB1Seg_0p4",
    "dtRechitCluster_match_MB1Seg_0p5",
    "dtRechitCluster_match_MB1hits_0p4",
    "dtRechitCluster_match_MB1hits_0p5",
    "dtRechitCluster_match_MB1hits_cosmics_minus",
    "dtRechitCluster_match_MB1hits_cosmics_plus",
    "dtRechitCluster_match_RB1_0p4",
    "dtRechitCluster_match_RB1_dPhi0p5",
    "dtRechitCluster_match_RPCBx_dPhi0p5",
    "dtRechitCluster_match_RPChits_dPhi0p5",
    # 'dtRechitCluster_match_gLLP',
    # 'dtRechitCluster_match_gLLP_csc',
    # 'dtRechitCluster_match_gLLP_decay_r',
    # 'dtRechitCluster_match_gLLP_decay_z',
    # 'dtRechitCluster_match_gLLP_dt',
    # 'dtRechitCluster_match_gLLP_e',
    # 'dtRechitCluster_match_gLLP_eta',
    # 'dtRechitCluster_match_gLLP_index',
    # 'dtRechitCluster_match_gLLP_minDeltaR',
    # 'dtRechitCluster_match_gLLP_phi',
    #
    "evtNum",
    # 'gHiggsE',
    # 'gHiggsEta',
    # 'gHiggsPhi',
    # 'gHiggsPt',
    # 'gLLP_beta',
    # 'gLLP_csc',
    # 'gLLP_ctau',
    # 'gLLP_decay_vertex_r',
    # 'gLLP_decay_vertex_x',
    # 'gLLP_decay_vertex_y',
    # 'gLLP_decay_vertex_z',
    # 'gLLP_dt',
    # 'gLLP_e',
    # 'gLLP_eta',
    # 'gLLP_phi',
    # 'gLLP_pt',
    "jetE",
    "jetEta",
    "jetPhi",
    "jetPt",
    "jetTightPassId",
    "lepDZ",
    "lepE",
    "lepEta",
    "lepPassLooseIso",
    "lepPassTightIso",
    "lepPassVTightIso",
    "lepPassVVTightIso",
    "lepPdgId",
    "lepPhi",
    "lepPt",
    "lepTightId",
    # 'lumiSec',
    # 'mH',
    # 'mX',
    "met",
    "metPhi",
    # 'nCscRechitClusters',
    # 'nCscRings',
    # 'nDtRechitClusters',
    # 'nDtRings',
    # 'nGLLP',
    "nJets",
    "nLeptons",
    # 'npu',
    # 'npv',
    # 'pileupWeight',
    # 'rho',
    # 'runNum',
    "weight",
]

# **** #
pathlib.Path(OUT_DIR).mkdir(parents=True, exist_ok=True)

rng = np.random.default_rng()

rt.gErrorIgnoreLevel = 1001  # rt.kInfo + 1
rt.gROOT.SetBatch(True)
PI = rt.TMath.Pi()

gc = []  # ROOT garbage collector


# **************************** #
def optimize_cut(_rdfs, col, cut_func, cut_values=None, cut0=None):
    vmc = _rdfs["mc"].AsNumpy(["weight", col])  # ,'tag_dPhi','dtSize'
    vr3 = _rdfs["r3"].AsNumpy(["weight", col])  # ,'tag_dPhi','dtSize'
    wmc, vmc = vmc["weight"], vmc[col]
    wr3, vr3 = vr3["weight"], vr3[col]
    # wmc_tot, wr3_tot = wmc.sum(), wr3.sum() # may have issues when cut_func has preselection?
    # vmin, vmax = min(vmc.min(), vr3.min()), max(vmc.max(), vr3.max())

    # **** #
    nmc, nr3 = _rdfs["mc"].Count().GetValue(), _rdfs["r3"].Count().GetValue()
    if isinstance(OPT_FRAC, float):
        idxs_mc = np.random.randint(0, nmc, int(OPT_FRAC * nmc))
        idxs_r3 = np.random.randint(0, nr3, int(OPT_FRAC * nr3))
    elif OPT_FRAC == "sqrt":
        idxs_mc = np.random.randint(0, nmc, int(nmc**0.5))
        idxs_r3 = np.random.randint(0, nr3, int(nr3**0.5))
    wmc, vmc = wmc[idxs_mc], vmc[idxs_mc]
    wr3, vr3 = wr3[idxs_r3], vr3[idxs_r3]

    # **** #
    # if 'ABCD' in CUTS:
    #     ABCD_DPHI = CUT_VALUES['ABCD_DPHI']
    #     ABCD_SIZE = CUT_VALUES['ABCD_SIZE']
    #     a_mc, b_mc, c_mc, d_mc =

    wmc0, wr30 = wmc[cut_func(cut0, vmc)].sum(), wr3[cut_func(cut0, vr3)].sum()
    wmc = np.array([wmc[cut_func(v, vmc)].sum() for v in cut_values])
    wr3 = np.array([wr3[cut_func(v, vr3)].sum() for v in cut_values])

    wmc_tot, wr3_tot = (
        wmc.max(),
        wr3.max(),
    )  # relative to itself? sorta fixes preselection bias
    # **** #
    # Scoring Metrics
    # S/rtB, sometimes falls into minimum where B=1 and S=small
    # s2b0 = wmc0/(wr30**0.5) if wr30 else 0
    # s2b = np.divide(wmc, np.sqrt(wr3), out=np.zeros_like(wmc), where=(wr3>0))
    # # F1 to balance precision/recall (avoid small S minimum?)
    # f10 = wmc0 /(2*wmc0+wr30+(wmc_tot-wmc0)) if 2*wmc0+wr30+(wmc_tot-wmc0) else 0
    # f1 = np.divide(wmc, 2*wmc+wr3+(wmc_tot-wmc),out=np.zeros_like(wmc), where=2*wmc+wr3+(wmc_tot-wmc)>0)
    # # Modified F-score where precision = S / (S + rtB)
    # f10 = wmc0 /(2*wmc0+(wr30**0.5)+(wmc_tot-wmc0)) if 2*wmc0+(wr30**0.5)+(wmc_tot-wmc0) else 0
    # f1 = np.divide(wmc, 2*wmc+np.sqrt(wr3)+(wmc_tot-wmc),out=np.zeros_like(wmc), where=2*wmc+np.sqrt(wr3)+(wmc_tot-wmc)>0)
    # # Modified F-score where precision = S / rtB, harmonic mean of S/rtB and S recall
    # f10 = wmc0 /(wmc_tot+(wr30**0.5)) if wmc_tot+(wr30**0.5) else 0
    # f1 = np.divide(wmc, wmc_tot+np.sqrt(wr3),out=np.zeros_like(wmc), where=wmc_tot+np.sqrt(wr3)>0)

    # https://arxiv.org/pdf/hep-ph/0204326.pdf
    s2b0 = 2 * (np.sqrt(wmc0 + wr30) - np.sqrt(wr30))
    s2b = 2 * (np.sqrt(wmc + wr3) - np.sqrt(wr3))

    # **** #
    score0, score = s2b0, s2b
    # score0, score = s2b0, s2b
    # score0, score = f10, f1
    # score0, score = s2b120, s2b12

    # # Find best at this working point
    # idx = np.argmax(score)
    # Randomly select a new cut that's better than before
    idx = rng.choice(np.arange(len(score))[score >= score0]) if (score >= score0).sum() else 0

    # **** #
    if cut0 is not None:
        return cut_values[idx], s2b[idx], wmc[idx], wr3[idx], s2b0, wmc0, wr30
    return cut_values[idx], s2b[idx], wmc[idx], wr3[idx]


# **************************** #
if __name__ == "__main__":
    print("+------------------------+")
    print("| Starting skim_cscdt.py |")
    print("+------------------------+")

    rt.EnableImplicitMT(0)
    print("    Enabled ROOT's implicit multithreading (sometimes causes a crash)")

    # **************************** #
    N_ITERATIONS = 1

    OPT_SCORES = []
    OPT_FRAC = 1.0  # if 1 just bootstrap, fraction with replacement
    LOO = False
    RAND = False

    PRINT_CUTFLOW = False
    MET_CATEGORY = ("lt200", "low", "high")[0]
    CUTSET = "scs"
    OOT = False

    DNN_VERSION = (None, "bkgMC", "bkgMC_plusBeamHalo", "bkgOOTData")[0]

    args = " ".join(sys.argv[1:]) if len(sys.argv) > 1 else ""

    if "cutflow" in args:
        print("    Printing cutflow tables")
        PRINT_CUTFLOW = True

    if "l1" in args:
        print("    Using the reduced cut set (l1)")
        CUTS = CUTS_L1
        CUTSET = "l1"

    if " low" in args:  # need space otherwise cutflow breaks it
        print("    Low met category")
        CUTS = [c.replace("MET", "low MET") if "MET" == c else c for c in CUTS]
        MET_CATEGORY = "low"
    elif "high" in args:
        print("    High met category")
        print("        REMOVING HALO, JET, & MUON VETOS")
        CUTS = [c.replace("MET", "high MET") if "MET" == c else c for c in CUTS]
        CUTS = [c for c in CUTS if "halo" not in c]
        CUTS = [c for c in CUTS if "jet" not in c]
        CUTS = [c for c in CUTS if "muon" not in c]
        OPT_CUTS = [c for c in OPT_CUTS if "halo" not in c]
        OPT_CUTS = [c for c in OPT_CUTS if "jet" not in c]
        OPT_CUTS = [c for c in OPT_CUTS if "muon" not in c]
        MET_CATEGORY = "high"
    else:
        print("    No met categorization (only met<200)")
        MET_CATEGORY = "lt200"

    if "roptDNN" in args:
        print("    Using the randomly optimized cut set with DNN (roptDNN)")
        CUTSET = "roptDNN"
    elif "ropt" in args:
        print("    Using the randomly optimized cut set (ropt)")
        CUTSET = "ropt"
    elif "loptDNN" in args:
        print("    Using the leave-one-out optimized cut set with DNN (loptDNN)")
        CUTSET = "loptDNN"
    elif "lopt" in args:
        print("    Using the leave-one-out optimized cut set (lopt)")
        CUTSET = "lopt"
    else:
        print("    Using the standard cut set (scs)")
        CUTSET = "scs"

    if "oot" in args:
        print("    Using out-of-time 2nd cluster")
        CUTS = [c.replace("DT IT", "DT OOT") if "DT IT" == c else c for c in CUTS]
        OOT = True
    else:
        print("    Using in-time 2nd cluster")

    if "loo" in args:
        print("    PERFORMING LOO OPTIMIZATION")
        if OOT:
            print("        FORCING MC TO IN-TIME")
        LOO = True
        N_ITERATIONS = 501
    if "rand" in args:
        print("    PERFORMING RAND OPTIMIZATION")
        if OOT:
            print("        FORCING MC TO IN-TIME")
        RAND = True
        N_ITERATIONS = 501
    # else:
    #     print('    REMOVING DT SIZE CUT')
    #     CUTS = [c for c in CUTS if 'DT size' not in c]

    if "DNN" not in CUTSET:
        print("    Removing DNN from CUTS")
        CUTS = [c for c in CUTS if "DNN" not in c]
        OPT_CUTS = [c for c in OPT_CUTS if "DNN" not in c]
    else:
        if RAND or LOO:
            print("        OPTIMIZING JUST DNN")
            OPT_CUTS = [c for c in OPT_CUTS if "DNN" in c]

    if "bkgMC_plusBeamHalo" in args:
        print("    Using DNN version bkgMC_plusBeamHalo")
        DNN_VERSION = "bkgMC_plusBeamHalo"
    elif "bkgMC" in args:
        print("    Using DNN version bkgMC")
        DNN_VERSION = "bkgMC"
        CUT_OPT_PARS["CSC DNN"]["values"] = np.arange(0.830, 1.001, 0.001)
    elif "bkgOOTData" in args:
        print("    Using DNN version bkgOOTData")
        DNN_VERSION = "bkgOOTData"
    if DNN_VERSION is not None:
        CUT_OPT_PARS["CSC DNN"]["col"] = f"cscDNN_{DNN_VERSION}"

    if f"{CUTSET}_{MET_CATEGORY}" in CUT_VALUES:
        CUT_VALUES = CUT_VALUES[f"{CUTSET}_{MET_CATEGORY}"]
    elif f"{CUTSET}_lt200" in CUT_VALUES:
        CUT_VALUES = CUT_VALUES[f"{CUTSET}_lt200"]
        print(f"    No cutset for {MET_CATEGORY=} found, defaulting to 'lt200'")
    else:
        raise ValueError(f"no cutset associated with {CUTSET=} and {MET_CATEGORY=} found")

    print("")

    # **** #
    if "1 CSC-DT" not in CUTS:
        raise NotImplementedError("cant handle multiple pairs yet")

    # **************************** #
    for iopt in range(N_ITERATIONS):
        # if iopt % 100>0 and iopt != N_ITERATIONS-1:
        if iopt != 0 and iopt != N_ITERATIONS - 1 and (LOO or RAND):
            OPT_CUT = rng.choice([c for c in OPT_CUTS if c != OPT_CUT or len(OPT_CUTS) == 1])
        else:
            OPT_CUT = ""

        if RAND and OPT_CUT != "":
            cv_history = [x[-1] for x in OPT_SCORES]
            cn = CUT_OPT_PARS[OPT_CUT]["cut"]
            val0 = CUT_VALUES[cn]
            skip_vals = [val0]
            for i in range(len(CUT_OPT_PARS[OPT_CUT]["values"]) - 1):
                val = rng.choice([v for v in CUT_OPT_PARS[OPT_CUT]["values"] if v not in skip_vals])
                CUT_VALUES[cn] = val
                if CUT_VALUES not in cv_history:
                    break
                skip_vals.append(val)
            if len(skip_vals) == len(CUT_OPT_PARS[OPT_CUT]["values"]):
                continue

        MIN_CSC_TIME = CUT_VALUES["MIN_CSC_TIME"]
        MAX_CSC_TIME = CUT_VALUES["MAX_CSC_TIME"]
        MAX_CSC_TSPREAD = CUT_VALUES["MAX_CSC_TSPREAD"]
        MAX_RPC_BX = CUT_VALUES["MAX_RPC_BX"]
        MIN_RPC_HITS = CUT_VALUES["MIN_RPC_HITS"]
        MAX_CSC_JET = CUT_VALUES["MAX_CSC_JET"]
        MAX_DT_JET = CUT_VALUES["MAX_DT_JET"]
        MAX_CSC_MUON = CUT_VALUES["MAX_CSC_MUON"]
        MAX_DT_MUON = CUT_VALUES["MAX_DT_MUON"]
        MAX_ME1 = CUT_VALUES["MAX_ME1"]
        MAX_MB1 = CUT_VALUES["MAX_MB1"]
        HALO_CUTOFF = CUT_VALUES["HALO_CUTOFF"]
        MIN_DPHI = CUT_VALUES["MIN_DPHI"]
        MIN_DETA = CUT_VALUES["MIN_DETA"]
        MAX_DETA = CUT_VALUES["MAX_DETA"]
        MIN_CSC_DNN = CUT_VALUES["MIN_CSC_DNN"]
        # MIN_DT_DNN = CUT_VALUES["MIN_DT_DNN"]
        # ABCD_DPHI = CUT_VALUES['ABCD_DPHI']
        # ABCD_SIZE = CUT_VALUES['ABCD_SIZE']

        rdfs = {
            "mc": RDataFrame("MuonSystem", FN_MC),
            "r3": RDataFrame("MuonSystem", FN_R3),
        }

        if iopt == 0:
            print("Events in:")
        for key, rdf in rdfs.items():
            # rdf = rdf.Range(0,100_000) # Skim a subset of events for debugging

            if key == "mc":  # fix weights
                rdf = rdf.Redefine("weight", f"weight * {LUMI}")

            count, weight = rdf.Count().GetValue(), rdf.Sum("weight").GetValue()
            if iopt == 0:
                print(f"  {key} = {count:,} ({weight:,.2f}) -- read")

            # **** #
            # Create dummy columns to store what cluster indices pass our selections
            rdf = rdf.Define("evtFlag", "weight > 0")
            rdf = rdf.Define(f"{C}Flag", f"{C}Size > 0")
            rdf = rdf.Define(f"{D}Flag", f"{D}Size > 0")

            rdf = rdf.Define("evtCutFlag", "weight > 0")
            rdf = rdf.Define(f"{C}CutFlag", f"{C}Size > 0")
            rdf = rdf.Define(f"{D}CutFlag", f"{D}Size > 0")

            # rdf = rdf.Define("cscDNN", f"Take({C}DNN_{DNN_VERSION},nCscRechitClusters)")
            # min_dnn = rdf.Min(f"cscDNN").GetValue()
            # max_dnn = rdf.Max(f"cscDNN").GetValue()
            # mean_dnn = rdf.Mean(f"cscDNN").GetValue()
            # print(f"{key} : {mean_dnn:.4f}, ({min_dnn:.4f}, {max_dnn:.4f})")

            # **** #
            # Add cluster radius columns
            for p in (C, D):
                rdf = rdf.Define(f"{p}R", f"ROOT::VecOps::sqrt( {p}X*{p}X + {p}Y*{p}Y )")
                if f"{p}R" not in COLUMNS_OUT:
                    COLUMNS_OUT.append(f"{p}R")

            rdf = rdf.Define(
                f"{C}NHitME1",
                f"{C}NRechitChamberPlus11 + {C}NRechitChamberPlus12 + {C}NRechitChamberMinus11 + {C}NRechitChamberMinus12",
            )
            if f"{C}NHitME1" not in COLUMNS_OUT:
                COLUMNS_OUT.append(f"{C}NHitME1")

            # **** #
            rdfs[key] = rdf

        #! Find original weights for optimization (will be overwritten when/if I reset indices)
        wmc_nocuts = rdfs["mc"].Sum("weight").GetValue()
        wr3_nocuts = rdfs["r3"].Sum("weight").GetValue()

        if iopt == 0:
            print("")

        # **************** #
        # Apply selections from CUTS (and/or arguments?)
        for key, rdf in rdfs.items():
            if PRINT_CUTFLOW:
                print(f"{key} cutflow:")
                print(r"\begin{center}")
                print(r"\begin{tabular}{c|ccc|ccc|ccc}")
                ec0, cc0, dc0 = (
                    rdf.Filter("evtFlag").Sum("weight"),
                    rdf.Filter("evtFlag").Sum(f"{C}Flag"),
                    rdf.Filter("evtFlag").Sum(f"{D}Flag"),
                )
                ec0, cc0, dc0 = ec0.GetValue(), cc0.GetValue(), dc0.GetValue()
                print(r"    \hline")
                print(
                    f'    \\textbf{{{"Signal" if "mc" in key else "Data"}}}'
                    + r" & \multicolumn{3}{c}{CSC-DT} & \multicolumn{3}{c}{CSC} & \multicolumn{3}{c}{DT} \\"
                )
                print(r"    \hline")
                print(
                    r"    Selection & N Events & Cut Eff & Cum Eff & N Clusters & Cut Eff & Cum Eff & N Clusters & Cut Eff & Cum Eff \\"
                )
                print(r"    \hline")
                print(f"    sample & {ec0:,.0f} & -- & -- & {cc0:,} & -- & -- & {dc0:,} & -- & -- \\\\")

            for cut in CUTS:
                # **** #
                # Reset cut flags
                rdf = rdf.Redefine("evtCutFlag", "weight > 0")
                rdf = rdf.Redefine(f"{C}CutFlag", f"{C}Size > 0")
                rdf = rdf.Redefine(f"{D}CutFlag", f"{D}Size > 0")

                # **** #
                # Signal matching for MC
                if "acceptance" in cut:
                    if "mc" in key:
                        rdf = rdf.Redefine(
                            f"{C}CutFlag",
                            f"{C}CutFlag && {C}_match_gLLP && "
                            + f"(abs({C}_match_gLLP_eta) < 3) && "
                            + f"({C}_match_gLLP_decay_r < 800) && "
                            + f"(400 < abs({C}_match_gLLP_decay_z)) &&"
                            + f" (abs({C}_match_gLLP_decay_z) < 1200)",
                        )

                        rdf = rdf.Redefine(
                            f"{D}CutFlag",
                            f"{D}CutFlag && {D}_match_gLLP && "
                            + f"(200 < {D}_match_gLLP_decay_r) && "
                            + f"({D}_match_gLLP_decay_r < 800) && "
                            + f"(abs({D}_match_gLLP_decay_z) < 700)",
                        )
                    elif "r3" in key:
                        continue

                # **** #
                # Trigger selections (HLT may be wrong in MC?)
                if "HLT" in cut:
                    if "mc" in key and "hlt" not in FN_MC:
                        continue
                    elif "r3" in key and "hlt" not in FN_R3:
                        rdf = rdf.Redefine(
                            "evtCutFlag", "evtCutFlag && HLTDecision[569]"
                        )  # 569 = HLT_L1CSCCluster_DTCluster50
                    else:
                        continue

                if "L1" in cut:
                    rdf = rdf.Redefine(
                        f"{C}CutFlag",
                        f"{C}CutFlag && ( "
                        + f"((100 < {C}R) && ({C}R < 275) && (580 < abs({C}Z)) && (abs({C}Z) < 632) && (500 <= {C}Size)) || "  # ME 11
                        + f"((275 < {C}R) && ({C}R < 465) && (668 < abs({C}Z)) && (abs({C}Z) < 724) && (200 <= {C}Size)) || "  # ME 12
                        + f"((505 < {C}R) && ({C}R < 700) && (668 < abs({C}Z)) && (abs({C}Z) < 724) && (200 <= {C}Size)) || "  # ME 13
                        +
                        #
                        f"((139 < {C}R) && ({C}R < 345) && (789 < abs({C}Z)) && (abs({C}Z) < 850) && (500 <= {C}Size)) || "  # ME 21
                        + f"((357 < {C}R) && ({C}R < 700) && (791 < abs({C}Z)) && (abs({C}Z) < 850) && (200 <= {C}Size)) || "  # ME 22
                        +
                        #
                        f"((160 < {C}R) && ({C}R < 345) && (915 < abs({C}Z)) && (abs({C}Z) < 970) && (500 <= {C}Size)) || "  # ME 31
                        + f"((357 < {C}R) && ({C}R < 700) && (911 < abs({C}Z)) && (abs({C}Z) < 970) && (200 <= {C}Size)) || "  # ME 32
                        +
                        #
                        f"((178 < {C}R) && ({C}R < 345) && (1002 < abs({C}Z)) && (abs({C}Z) < 1063) && (500 <= {C}Size)) || "  # ME 41
                        + f"((357 < {C}R) && ({C}R < 700) && (1002 < abs({C}Z)) && (abs({C}Z) < 1063) && (200 <= {C}Size)) )",
                    )  # ME 42

                # **** #
                # Missing Et requirements & categorization
                if "MET" in cut:
                    rdf = rdf.Redefine("evtCutFlag", f"evtCutFlag && (met < 200)")

                if "low MET" in cut:
                    rdf = rdf.Redefine("evtCutFlag", f"evtCutFlag && (met < {LOW_MET_CUTOFF})")

                if "high MET" in cut:
                    rdf = rdf.Redefine(
                        "evtCutFlag",
                        f"evtCutFlag && (({HIGH_MET_CUTOFF} < met) && (met < 200))",
                    )

                # **** #
                # Cluster time requirements
                if (LOO or RAND) and "mc" in key and " OOT" in cut:
                    cut = cut.replace(" OOT", " IT")

                if "CSC IT" in cut:
                    rdf = rdf.Redefine(
                        f"{C}CutFlag",
                        f"{C}CutFlag && ( ({MIN_CSC_TIME} < {C}TimeWeighted) && ({C}TimeWeighted < {MAX_CSC_TIME}) && ({C}TimeSpreadWeightedAll < {MAX_CSC_TSPREAD}) )",
                    )
                elif "CSC OOT" in cut:
                    rdf = rdf.Redefine(
                        f"{C}CutFlag",
                        f"{C}CutFlag && !( ({MIN_CSC_TIME} < {C}TimeWeighted) && ({C}TimeWeighted < {MAX_CSC_TIME}) && ({C}TimeSpreadWeightedAll < {MAX_CSC_TSPREAD}) )",
                    )

                if "DT IT" in cut:
                    rdf = rdf.Redefine(
                        f"{D}CutFlag",
                        f"{D}CutFlag && ( (abs({D}_match_RPCBx_dPhi0p5) <= {MAX_RPC_BX}) && ({D}_match_RPChits_dPhi0p5 >= {MIN_RPC_HITS}) )",
                    )
                elif "DT OOT" in cut:
                    rdf = rdf.Redefine(
                        f"{D}CutFlag",
                        f"{D}CutFlag && !( (abs({D}_match_RPCBx_dPhi0p5) <= {MAX_RPC_BX}) && ({D}_match_RPChits_dPhi0p5 >= {MIN_RPC_HITS}) )",
                    )

                # **** #
                # Station requirements from L1 trigger
                if "ME1" in cut and (OPT_CUT != "ME1" or not LOO):
                    # rdf = rdf.Redefine(f'{C}CutFlag', f'{C}CutFlag && ( '+
                    #                 f'({C}NRechitChamberPlus11 <= {MAX_ME1}) &&'+
                    #                 f'({C}NRechitChamberPlus12 <= {MAX_ME1}) &&'+
                    #                 f'({C}NRechitChamberMinus11 <= {MAX_ME1}) &&'+
                    #                 f'({C}NRechitChamberMinus12 <= {MAX_ME1}) )')
                    rdf = rdf.Redefine(f"{C}CutFlag", f"{C}CutFlag && ( {C}NHitME1 <= {MAX_ME1} ) ")

                if "MB1" in cut and (OPT_CUT != "MB1" or not LOO):
                    rdf = rdf.Redefine(f"{D}CutFlag", f"{D}CutFlag && ( {D}NHitStation1 <= {MAX_MB1} )")

                if "DT stn" in cut:
                    rdf = rdf.Redefine(
                        f"{D}CutFlag",
                        f"{D}CutFlag && ( ({D}NStation10 < 3) && !(({D}NStation10 == 2) && ({D}MaxStation == 4)) )",
                    )

                # **** #
                if "jet veto" in cut:
                    if OPT_CUT != "CSC jet veto" or not LOO:
                        rdf = rdf.Redefine(
                            f"{C}CutFlag",
                            f"{C}CutFlag && ({C}JetVetoPt < {MAX_CSC_JET})",
                        )
                    if OPT_CUT != "DT jet veto" or not LOO:
                        rdf = rdf.Redefine(
                            f"{D}CutFlag",
                            f"{D}CutFlag && ({D}JetVetoPt < {MAX_DT_JET})",
                        )

                if "muon veto" in cut:
                    rdf = rdf.Redefine(f"{C}CutFlag", f"{C}CutFlag && ({C}MuonVetoGlobal == 0)")
                    rdf = rdf.Redefine(f"{D}CutFlag", f"{D}CutFlag && ({D}MuonVetoLooseId == 0)")

                    if OPT_CUT != "CSC muon veto" or not LOO:
                        rdf = rdf.Redefine(
                            f"{C}CutFlag",
                            f"{C}CutFlag && ({C}MuonVetoPt < {MAX_CSC_MUON})",
                        )
                    if OPT_CUT != "DT muon veto" or not LOO:
                        rdf = rdf.Redefine(
                            f"{D}CutFlag",
                            f"{D}CutFlag && ({D}MuonVetoPt < {MAX_DT_MUON})",
                        )

                if "halo veto" in cut and (OPT_CUT != "halo veto" or not LOO):
                    rdf = rdf.Redefine(
                        f"{D}CutFlag",
                        f"{D}CutFlag && ( ({HALO_CUTOFF} < abs({D}Phi)) && (abs({D}Phi) < {PI} - {HALO_CUTOFF}) )",
                    )

                # **** #
                if "BDT" in cut:
                    raise NotImplementedError("BDT")

                if "DNN" in cut:
                    if True:  # OPT_CUT != "CSC DNN" or not LOO:
                        rdf = rdf.Redefine(
                            f"{C}CutFlag",
                            f"{C}CutFlag && ( Take({C}DNN_{DNN_VERSION},nCscRechitClusters) > {MIN_CSC_DNN} )",
                        )
                    # if OPT_CUT != "DT DNN" or not LOO:
                    #     rdf = rdf.Redefine(
                    #         f"{D}CutFlag",
                    #         f"{D}CutFlag && ( Take({D}DNN,nDtRechitClusters) > {MIN_DT_DNN} )",
                    #     )

                # **** #
                if "1 CSC-DT" in cut:
                    # rdf = rdf.Redefine('evtCutFlag', f'evtCutFlag && evtFlag && '+
                    #                    f'(reduce({C}Flag.begin(), {C}Flag.end()) == 1) && '+
                    #                    f'(reduce({D}Flag.begin(), {D}Flag.end()) == 1)')
                    # if LOO:
                    #     icsc, idt = rng.random(2)
                    #     rdf = rdf.Redefine(f'{C}Flag', f'{C}Flag && ( {C}Size == {C}Size[int({icsc}*reduce({C}Flag.begin(),{C}Flag.end()))] )')
                    #     rdf = rdf.Redefine(f'{D}Flag', f'{D}Flag && ( {D}Size == {D}Size[int({idt}*reduce({D}Flag.begin(),{D}Flag.end()))] )')
                    # else:
                    rdf = rdf.Redefine(f"{C}Flag", f"{C}Flag && ( {C}Size == Max({C}Size*{C}Flag) )")
                    rdf = rdf.Redefine(f"{D}Flag", f"{D}Flag && ( {D}Size == Max({D}Size*{D}Flag) )")

                    # Apply our cluster level selections to relevant columns
                    columns_out = []
                    for col in COLUMNS_OUT:
                        ocol = col.replace("RechitCluster", "")
                        if C in col[: len(C)]:
                            rdf = rdf.Define(ocol, f"{col}[{C}Flag][0]")
                        elif D in col[: len(D)]:
                            rdf = rdf.Define(ocol, f"{col}[{D}Flag][0]")
                        columns_out.append(ocol)

                    # rdf = rdf.Define("cscDNN", f"Take({C}DNN,nCscRechitClusters)[{C}Flag][0]")
                    # rdf = rdf.Define("dtDNN", f"Take({D}DNN,nDtRechitClusters)[{D}Flag][0]")
                    rdf = rdf.Define("cscDNN_bkgMC", f"Take({C}DNN_bkgMC,nCscRechitClusters)[{C}Flag][0]")
                    rdf = rdf.Define(
                        "cscDNN_bkgMC_plusBeamHalo", f"Take({C}DNN_bkgMC_plusBeamHalo,nCscRechitClusters)[{C}Flag][0]"
                    )
                    rdf = rdf.Define("cscDNN_bkgOOTData", f"Take({C}DNN_bkgOOTData,nCscRechitClusters)[{C}Flag][0]")

                    # min_dnn = rdf.Min(f"cscDNN_{DNN_VERSION}").GetValue()
                    # max_dnn = rdf.Max(f"cscDNN_{DNN_VERSION}").GetValue()
                    # mean_dnn = rdf.Mean(f"cscDNN_{DNN_VERSION}").GetValue()
                    # print(f"{key} : {mean_dnn:.4f}, ({min_dnn:.4f}, {max_dnn:.4f})")

                    # rdf = rdf.Define("dtDNN_bkgMC", f"Take({D}DNN_bkgMC,nDtRechitClusters)[{D}Flag][0]")
                    # rdf = rdf.Define(
                    #     "dtDNN_bkgMC_plusBeamHalo", f"Take({D}DNN_bkgMC_plusBeamHalo,nDtRechitClusters)[{D}Flag][0]"
                    # )
                    # rdf = rdf.Define("dtDNN_bkgOOTData", f"Take({D}DNN_bkgOOTData,nDtRechitClusters)[{D}Flag][0]")

                    rdf = rdf.Define("cscCTau", f"gLLP_ctau[{C}_match_gLLP_index[{C}Flag][0]]")
                    rdf = rdf.Define("dtCTau", f"gLLP_ctau[{D}_match_gLLP_index[{D}Flag][0]]")
                    rdf = rdf.Define("tag_dEta", "abs(cscEta - dtEta)")
                    rdf = rdf.Define(
                        "tag_dPhi",
                        f"abs( abs(cscPhi - dtPhi) - (abs(cscPhi - dtPhi) > {PI} ? 2*{PI} : 0) )",
                    )
                    rdf = rdf.Define("tag_dR", "sqrt(tag_dEta*tag_dEta + tag_dPhi*tag_dPhi)")
                    columns_out.extend(
                        [
                            "tag_dEta",
                            "tag_dPhi",
                            "tag_dR",
                            "cscCTau",
                            "dtCTau",
                            # "cscDNN",
                            # "dtDNN",
                            "cscDNN_bkgMC",
                            "cscDNN_bkgMC_plusBeamHalo",
                            "cscDNN_bkgOOTData",
                            # "dtDNN_bkgMC",
                            # "dtDNN_bkgMC_plusBeamHalo",
                            # "dtDNN_bkgOOTData",
                        ]
                    )

                # **** #
                if "dR" in cut:
                    raise NotImplementedError("dR")

                if "dEta" in cut:
                    if OPT_CUT != "min dEta" or not LOO:
                        rdf = rdf.Redefine("evtCutFlag", f"evtCutFlag && (tag_dEta > {MIN_DETA})")
                    if OPT_CUT != "max dEta" or not LOO:
                        rdf = rdf.Redefine("evtCutFlag", f"evtCutFlag && (tag_dEta < {MAX_DETA})")

                if "dPhi" in cut and (OPT_CUT != "min dPhi" or not LOO):
                    rdf = rdf.Redefine("evtCutFlag", f"evtCutFlag && (tag_dPhi > {MIN_DPHI})")

                # **** #
                # Propagate to cumulative flags
                rdf = rdf.Redefine(f"{C}Flag", f"{C}Flag && {C}CutFlag")
                rdf = rdf.Redefine(f"{D}Flag", f"{D}Flag && {D}CutFlag")

                rdf = rdf.Redefine(
                    "evtFlag",
                    "evtFlag && evtCutFlag && ("
                    + f"(reduce({C}Flag.begin(), {C}Flag.end()) > 0) && "
                    + f"(reduce({D}Flag.begin(), {D}Flag.end()) > 0))",
                )

                # rdf = rdf.Redefine(f'nC{C[1:]}s', f'reduce({C}Flag.begin(), {C}Flag.end())')
                # rdf = rdf.Redefine(f'nD{D[1:]}s', f'reduce({D}Flag.begin(), {D}Flag.end())')

                # **** #
                if PRINT_CUTFLOW:
                    rdf = rdf.Redefine(
                        "evtCutFlag",
                        "evtCutFlag && ("
                        + f"(reduce({C}CutFlag.begin(), {C}CutFlag.end()) > 0) && "
                        + f"(reduce({D}CutFlag.begin(), {D}CutFlag.end()) > 0))",
                    )

                    ec, cc, dc = (
                        rdf.Filter("evtFlag").Sum("weight"),
                        rdf.Filter("evtFlag").Sum(f"{C}Flag"),
                        rdf.Filter("evtFlag").Sum(f"{D}Flag"),
                    )
                    # ecc,ccc,dcc = rdf.Filter('evtCutFlag').Sum('weight'), rdf.Filter('evtCutFlag').Sum(f'{C}CutFlag'), rdf.Filter('evtCutFlag').Sum(f'{D}CutFlag')
                    ecc, ccc, dcc = (
                        rdf.Filter("evtCutFlag").Sum("weight"),
                        rdf.Sum(f"{C}CutFlag"),
                        rdf.Sum(f"{D}CutFlag"),
                    )

                    ec, cc, dc = ec.GetValue(), cc.GetValue(), dc.GetValue()
                    ecc, ccc, dcc = ecc.GetValue(), ccc.GetValue(), dcc.GetValue()

                    c_cum_eff = f"{100*cc/cc0:.2f}" if cc != cc0 else "--"
                    c_cut_eff = f"{100*ccc/cc0:.2f}" if ccc != cc0 else "--"
                    d_cum_eff = f"{100*dc/dc0:.2f}" if dc != dc0 else "--"
                    d_cut_eff = f"{100*dcc/dc0:.2f}" if dcc != dc0 else "--"

                    print(
                        f'    {cut.replace("_"," ")} & {ec:,.0f} & {100*ecc/ec0:.2f} & {100*ec/ec0:.2f} & '
                        + f"{cc:,.0f} & {c_cut_eff} & {c_cum_eff} & "
                        + f"{dc:,.0f} & {d_cut_eff} & {d_cum_eff} \\\\"
                    )

                    # Only filter (reset indices) after the MET cut when printing cutflow
                    if "MET" in cut:
                        print(r"    \hline")
                        ec0, cc0, dc0 = ec, cc, dc

                        rdf = rdf.Filter("evtFlag")
                        rdf = rdf.Redefine(f"{C}Size", f"{C}Size * {C}Flag")
                        rdf = rdf.Redefine(f"{D}Size", f"{D}Size * {D}Flag")
                else:
                    rdf = rdf.Filter("evtFlag")
                    rdf = rdf.Redefine(f"{C}Size", f"{C}Size * {C}Flag")
                    rdf = rdf.Redefine(f"{D}Size", f"{D}Size * {D}Flag")

                #! Assumes we don't optimize all cuts with DNN, for some optimization targets such as f1
                if (LOO or RAND) and (("DT IT" in cut or "DT OOT" in cut) or ("DNN" in CUTSET and "halo" in cut)):
                    if "mc" in key:
                        wmc_nocuts = rdf.Sum("weight").GetValue()
                    else:
                        wr3_nocuts = rdf.Sum("weight").GetValue()

            rdfs[key] = rdf.Filter("evtFlag")  # Apply the final filter and update the dictionary with the new RDF

            if PRINT_CUTFLOW:
                print(r"    \hline")
                print(r"\end{tabular}")
                print(r"\end{center}")
                print("")

        if LOO:
            _c = False
            if OPT_CUT == "":
                val, wmc, wr3 = (
                    0,
                    rdfs["mc"].Sum("weight").GetValue(),
                    rdfs["r3"].Sum("weight").GetValue(),
                )
                # score = wmc/(wr3**0.5) if wr3 else 0 # if b~1 then s2b explodes and falls into a false minimum
                score = 2 * ((wmc + wr3) ** 0.5 - wr3**0.5)  # https://arxiv.org/pdf/hep-ph/0204326.pdf
                val0, score0, wmc0, wr30 = 0, score, wmc, wr3
            else:
                val0 = CUT_VALUES[CUT_OPT_PARS[OPT_CUT]["cut"]]
                val, score, wmc, wr3, score0, wmc0, wr30 = optimize_cut(
                    _rdfs=rdfs,
                    col=CUT_OPT_PARS[OPT_CUT]["col"],
                    cut_func=CUT_OPT_PARS[OPT_CUT]["func"],
                    cut_values=CUT_OPT_PARS[OPT_CUT]["values"],
                    cut0=CUT_VALUES[CUT_OPT_PARS[OPT_CUT]["cut"]],
                )
                if wmc >= 1 and wr3 >= 3 and score >= score0:
                    _c = True
                    # CUT_VALUES[CUT_OPT_PARS[OPT_CUT]['cut']] = val
                    # CUT_VALUES[CUT_OPT_PARS[OPT_CUT]['cut']] = (val*wmc + val0*wmc0)/(wmc + wmc0)
                    CUT_VALUES[CUT_OPT_PARS[OPT_CUT]["cut"]] = (
                        val * wmc / (5 * (1 - iopt / N_ITERATIONS)) + val0 * wmc0
                    ) / (wmc / (5 * (1 - iopt / N_ITERATIONS)) + wmc0)

                del rdfs["mc"]
                del rdfs["r3"]
                del rdf

            print(
                f'{iopt:>3} | {"Y" if _c else "X"} | {OPT_CUT:>13} = {val:>6.2f} ({val0:>6.2f}) | {score:>3.0f} ({score0:>3.0f}), {wmc:>4.0f} ({wmc0:>4.0f}), {wr3:>5.0f} ({wr30:>5.0f})'
            )

            OPT_SCORES.append(
                [
                    OPT_CUT,
                    val,
                    score,
                    wmc,
                    wr3,
                    val0,
                    score0,
                    wmc0,
                    wr30,
                    {k: v for k, v in CUT_VALUES.items()},
                ]
            )
            with open(
                f'loo_scores_cscdt{"OOT" if OOT else ""}_{CUTSET}_{MET_CATEGORY}.pkl',
                "wb",
            ) as fout:
                pickle.dump(OPT_SCORES, fout)
        if RAND:
            _c = False
            wmc = rdfs["mc"].Sum("weight").GetValue()
            wr3 = rdfs["r3"].Sum("weight").GetValue()
            # score = wmc / (wr3**0.5) if wr3 else 0  # if b~1 then s2b explodes and falls into a false minimum
            score = 2 * ((wmc + wr3) ** 0.5 - wr3**0.5)  # https://arxiv.org/pdf/hep-ph/0204326.pdf
            # score = 2*((wmc*2.16e-03+wr3)**0.5 - wr3**0.5)
            # score = 2*((wmc*(2.16e-03**(iopt/N_ITERATIONS))+wr3)**0.5 - wr3**0.5) # fix score0 for next iteration too
            # limit = ((1.5/2+np.sqrt(wr3))**2-wr3)/wmc if wmc else 1
            # limit = ((1.96/2+np.sqrt(wr3))**2-wr3)/wmc if wmc else 1
            # score = 2 * wmc / (wmc + (wr3 + (wmc_nocuts - wmc)) / 2)  # F1 Score, balanced precision/recall
            # score = wmc / (wmc + wr3)
            # wmc_ad = rdfs["mc"].Filter(f"dtSize > 85").Sum("weight").GetValue()
            # wr3_ad = rdfs["r3"].Filter(f"dtSize > 85").Sum("weight").GetValue()
            # score = wmc_ad / wr3_ad**0.5 if wr3_ad else 0

            if OPT_CUT == "":
                # limit0 = limit
                score0, wmc0, wr30 = score, wmc, wr3
                val, val0 = 999, 999
            else:
                # if wmc>0 and wr3>0 and limit <= limit0:# and wmc/wmc0 > 0.9:
                if wmc > 0 and wr3 > 0 and score >= score0:  # and wmc/wmc0 > 0.9:
                    _c = True
                    CUT_VALUES[CUT_OPT_PARS[OPT_CUT]["cut"]] = val
                else:
                    CUT_VALUES[CUT_OPT_PARS[OPT_CUT]["cut"]] = val0

            # print(f'{iopt:>3} | {"Y" if _c else "X"} | {OPT_CUT:>13} = {val:>6.2f} ({val0:>6.2f}) | {score:>3.0f} ({score0:>3.0f}), {wmc:>4.0f} ({wmc0:>4.0f}), {wr3:>5.0f} ({wr30:>5.0f})')
            print(
                f'{iopt:>3} | {"Y" if _c else "X"} | {OPT_CUT:>13} = {val:>7.3f} ({val0:>7.3f}) | {score:>4.2f} ({score0:>4.2f}), {wmc:>4.0f} ({wmc0:>4.0f}), {wr3:>5.0f} ({wr30:>5.0f})'
            )
            # print(f'{iopt:>3} | {"Y" if _c else "X"} | {OPT_CUT:>13} = {val:>6.2f} ({val0:>6.2f}) | {limit:>.2e} ({limit0:>.2e}), {wmc:>4.0f} ({wmc0:>4.0f}), {wr3:>5.0f} ({wr30:>5.0f})')

            OPT_SCORES.append(
                [
                    OPT_CUT,
                    val,
                    score,
                    wmc,
                    wr3,
                    val0,
                    score0,
                    wmc0,
                    wr30,
                    {k: v for k, v in CUT_VALUES.items()},
                ]
            )
            with open(
                f'rand_scores_cscdt{"OOT" if OOT else ""}_{CUTSET}_{MET_CATEGORY}.pkl',
                "wb",
            ) as fout:
                pickle.dump(OPT_SCORES, fout)

            if _c:
                score0, wmc0, wr30 = score, wmc, wr3
                # limit0 = limit
            # score0 = 2*((wmc0*(2.16e-03**((iopt+1)/N_ITERATIONS))+wr30)**0.5 - wr30**0.5)

    # **************** #
    if LOO or RAND:
        print("\nOptimized Cuts:")
        print(f"    {MIN_CSC_TIME=:.2f}")
        print(f"    {MAX_CSC_TIME=:.2f}")
        print(f"    {MAX_CSC_TSPREAD=:.2f}")
        print(f"    {MAX_RPC_BX=:.0f}")
        print(f"    {MIN_RPC_HITS=:.0f}")
        print(f"    {MAX_CSC_JET=:.0f}")
        print(f"    {MAX_DT_JET=:.0f}")
        print(f"    {MAX_CSC_MUON=:.0f}")
        print(f"    {MAX_DT_MUON=:.0f}")
        print(f"    {MAX_ME1=:.0f}")
        print(f"    {MAX_MB1=:.0f}")
        print(f"    {HALO_CUTOFF=:.2f}")
        print(f"    {MIN_DPHI=:.2f}")
        print(f"    {MIN_DETA=:.1f}")
        print(f"    {MAX_DETA=:.1f}")
        print(f"    {MIN_CSC_DNN=:.3f}")
        # print(f"    {MIN_DT_DNN=:.2f}")
        print("")

    print("Events out:")
    hists = {}
    for key, rdf in rdfs.items():
        count, weight = rdf.Count(), rdf.Sum("weight")
        count, weight = count.GetValue(), weight.GetValue()
        if count == 0:
            print(f"{key} is empty")
            continue

        name = f'{key}_cscdt{"OOT" if OOT else ""}_{CUTSET}'
        if DNN_VERSION is not None and "DNN" in CUTSET:
            name += f"_{DNN_VERSION}"
        name += f"_{MET_CATEGORY}"

        if LOO or RAND:
            name += f"_LOO" if LOO else f"_RAND"
            if key == "mc":
                name = name.replace("cscdtOOT", "cscdt")

        rdf = rdf.Snapshot("MuonSystem_flat", f"data/processed/{name}_rdf.root", columns_out)
        rdfs[key] = rdf

        # print(f'  {key} = {count:,} ({weight:,.2e})')
        print(f"  {key} = {count:,} ({weight:,.2f})")
        for xx in (
            "met",
            "cscSize",
            "cscR",
            "cscEta",
            "cscPhi",
            "dtSize",
            "dtR",
            "dtEta",
            "dtPhi",
            "tag_dEta",
            "tag_dPhi",
            "tag_dR",
        ):
            hh = rdf.Histo1D(xx, "weight")
            hh.SetName(f"{key}_{xx}")
            hh.SetTitle(f"{key};{xx};count")

            if xx not in hists:
                hists[xx] = [hh]
            else:
                hists[xx].append(hh)
    # for xx, hhs in hists.items():
    #     canvas = rt.TCanvas('','',800,800)
    #     # xmin, xmax, std = rdf.Min(xx).GetValue(), rdf.Max(xx).GetValue(), rdf.StdDev(xx).GetValue()
    #     # nbins = int((xmax-xmin)/(2.718*std*count**(-1/3))) if std else 1
    #     # hh = rdf.Histo1D((f'{key}_{xx}',f'{key};{xx};count',nbins,xmin,xmax),f'{xx}').GetValue()
    #     hhs = [h.GetValue() for h in hhs]
    #     hmax = max([h.GetMaximum() / h.Integral() if h.Integral() else 0 for h in hhs])
    #     xmin =
    #     for ih, hh in enumerate(hhs):
    #         if hmax < 1 and hh.Integral():
    #             hh.Scale(1/hh.Integral())
    #         hh.SetMinimum(0)
    #         hh.SetMaximum(hmax*1.1)
    #         hh.GetXaxis().SetRangeUser(xmin, xmax)
    #         hh.SetLineColor(SCL[ih])
    #         hh.SetLineWidth(3)
    #         hh.Draw('hist same')
    #     canvas.Draw()
    #     canvas.Print(f'{OUT_DIR}/{name}_{xx}.png')

    print("")
