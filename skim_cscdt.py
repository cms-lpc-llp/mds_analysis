import sys
import pathlib
import numpy as np

import pickle

import ROOT as rt
from ROOT import RDataFrame

from src.histo_utilities import std_color_list as SCL

# **************************** #
LOCAL_DIR = '/home/psimmerl/mds_analysis'
OUT_DIR = f'{LOCAL_DIR}/reports/weekly/2024-04-15'

# STAT = 'raw'
# LUMI = 23.02 * 1000
# FN_MC = f'{LOCAL_DIR}/data/raw/ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v6.root'
# FN_R3 = f'{LOCAL_DIR}/data/raw/DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi_v6.root'

STAT = 'pedro'
LUMI = 1.1328524540090597e-06 * 23.02 * 1000 # ???
FN_MC = f'{LOCAL_DIR}/data/raw/mc_pedro.root'
# FN_R3 = f'{LOCAL_DIR}/data/raw/data_pedro.root'
# FN_MC = f'{LOCAL_DIR}/data/processed/mc_pedro_hlt569.root'
FN_R3 = f'{LOCAL_DIR}/data/processed/r3_pedro_hlt569.root'


# **** #
LOW_MET_CUTOFF = 75
HIGH_MET_CUTOFF = 150

# **** #
CUTS_L1 = [
    'acceptance',
    'HLT',
    'L1',
    'MET',
    #! reset cutflow indices here
    'CSC IT',
    'DT IT',
    'MB1',
    '1 CSC-DT',
    # 'dPhi',
]

CUTS = [
    'acceptance',
    'HLT',
    'L1',
    'MET',
    #! reset cutflow indices here
    'CSC IT',
    'DT IT',
    # 'ME1',
    'MB1',
    # 'n jets',
    # 'n leptons',
    'jet veto',
    'muon veto',
    'halo veto',
    # 'DT N Station',
    # 'DT Avg Station',
    # 'BDT',
    'DNN',
    '1 CSC-DT',
    # 'dR',
    'dEta',
    'dPhi',
    # 'ABCD'
    'DT size' # for ABCD?
]

# **** #
CUT_VALUES = {
    'l1_lt200' : {
        'MIN_CSC_TIME' : -5.0,
        'MAX_CSC_TIME' : 12.5,
        'MAX_CSC_TSPREAD' : 20.0,
        'MAX_RPC_BX' : 0,
        'MIN_RPC_HITS' : 1,
        'MAX_CSC_JET' : 200,
        'MAX_DT_JET' : 200,
        'MAX_CSC_MUON' : 200,
        'MAX_DT_MUON' : 200,
        'MAX_ME1' : 0,
        'MAX_MB1' : 10,
        'HALO_CUTOFF' : 0,
        'MIN_DPHI' : 0.4,
        'MIN_DETA' : 0,
        'MAX_DETA' : 4,
        'MIN_CSC_DNN' : 0,
        'MIN_DT_DNN' : 0,
        'MIN_DT_SIZE' : 50,
    },
    'scs_lt200' : {
        'MIN_CSC_TIME' : -5.0,
        'MAX_CSC_TIME' :  12.5,
        'MAX_CSC_TSPREAD' :  20.0,
        'MAX_RPC_BX' : 0,
        'MIN_RPC_HITS' : 1,
        'MAX_CSC_JET' : 30,
        'MAX_DT_JET' : 50,
        'MAX_CSC_MUON' : 30,
        'MAX_DT_MUON' : 10,
        'MAX_ME1' : 0,
        'MAX_MB1' : 0,
        'HALO_CUTOFF' : 0.4,
        'MIN_DPHI' : 0.4,
        'MIN_DETA' : 0,
        'MAX_DETA' : 4,
        'MIN_CSC_DNN' : 0,
        'MIN_DT_DNN' : 0,
        'MIN_DT_SIZE' : 50,
    },
    'ropt_low' : {
        'MIN_CSC_TIME' : -5.0,
        'MAX_CSC_TIME' :  12.5,
        'MAX_CSC_TSPREAD' :  20.0,
        'MAX_RPC_BX' : 0,
        'MIN_RPC_HITS' : 1,
        'MAX_CSC_JET' : 100,
        'MAX_DT_JET' : 10,
        'MAX_CSC_MUON' : 90,
        'MAX_DT_MUON' : 60,
        'MAX_ME1' : 0,
        'MAX_MB1' : 1,
        'HALO_CUTOFF' : 0.15,
        'MIN_DPHI' : 0.4,
        'MIN_DETA' : 0,
        'MAX_DETA' : 4,
        'MIN_CSC_DNN' : 0,
        'MIN_DT_DNN' : 0,
        'MIN_DT_SIZE' : 50,
    },
    'ropt_high' : {
        'MIN_CSC_TIME' : -5.0,
        'MAX_CSC_TIME' :  12.5,
        'MAX_CSC_TSPREAD' :  20.0,
        'MAX_RPC_BX' : 0,
        'MIN_RPC_HITS' : 1,
        'MAX_CSC_JET' : 35,
        'MAX_DT_JET' : 10,
        'MAX_CSC_MUON' : 70,
        'MAX_DT_MUON' : 95,
        'MAX_ME1' : 0,
        'MAX_MB1' : 9,
        'HALO_CUTOFF' : 0.0,
        'MIN_DPHI' : 0.4,
        'MIN_DETA' : 0,
        'MAX_DETA' : 4,
        'MIN_CSC_DNN' : 0,
        'MIN_DT_DNN' : 0,
        'MIN_DT_SIZE' : 50,
    },
    'lopt_lt200' : {
        # Default guess
        'MIN_CSC_TIME' : -5.0,
        'MAX_CSC_TIME' : 12.5,
        'MAX_CSC_TSPREAD' : 20.0,
        'MAX_RPC_BX' : 0,
        'MIN_RPC_HITS' : 1,
        # 'MAX_N_JETS' : 15,
        # 'MAX_N_LEPS' : 5,
        'MAX_CSC_JET' : 100,
        'MAX_DT_JET' : 100,
        'MAX_CSC_MUON' : 100,
        'MAX_DT_MUON' : 100,
        'MAX_ME1' : 0,
        'MAX_MB1' : 10,
        'HALO_CUTOFF' : 0.2,
        'MIN_DPHI' : 0,
        'MIN_DETA' : 0,
        'MAX_DETA' : 4,
        'MIN_CSC_DNN' : 0,
        'MIN_DT_DNN' : 0,
        'MIN_DT_SIZE' : 50,
    },
    'loptDNN_lt200' : {
        # Default guess
        'MIN_CSC_TIME' : -5.0,
        'MAX_CSC_TIME' : 12.5,
        'MAX_CSC_TSPREAD' : 20.0,
        'MAX_RPC_BX' : 0,
        'MIN_RPC_HITS' : 1,
        # 'MAX_N_JETS' : 15,
        # 'MAX_N_LEPS' : 5,
        'MAX_CSC_JET' : 100,
        'MAX_DT_JET' : 100,
        'MAX_CSC_MUON' : 100,
        'MAX_DT_MUON' : 100,
        'MAX_ME1' : 0,
        'MAX_MB1' : 10,
        'HALO_CUTOFF' : 0.2,
        'MIN_DPHI' : 0,
        'MIN_DETA' : 0,
        'MAX_DETA' : 4,
        'MIN_CSC_DNN' : 0,
        'MIN_DT_DNN' : 0,
        'MIN_DT_SIZE' : 50,
    },
    'lopt_low' : {
        # Default guess
        'MIN_CSC_TIME' : -5.0,
        'MAX_CSC_TIME' : 12.5,
        'MAX_CSC_TSPREAD' : 20.0,
        'MAX_RPC_BX' : 0,
        'MIN_RPC_HITS' : 1,
        # 'MAX_N_JETS' : 15,
        # 'MAX_N_LEPS' : 5,
        'MAX_CSC_JET' : 100,
        'MAX_DT_JET' : 100,
        'MAX_CSC_MUON' : 100,
        'MAX_DT_MUON' : 100,
        'MAX_ME1' : 0,
        'MAX_MB1' : 10,
        'HALO_CUTOFF' : 0.2,
        'MIN_DPHI' : 0,
        'MIN_DETA' : 0,
        'MAX_DETA' : 4,
        'MIN_CSC_DNN' : 0,
        'MIN_DT_DNN' : 0,
        'MIN_DT_SIZE' : 50,
    },
    'loptDNN_low' : {
        # Default guess
        'MIN_CSC_TIME' : -5.0,
        'MAX_CSC_TIME' : 12.5,
        'MAX_CSC_TSPREAD' : 20.0,
        'MAX_RPC_BX' : 0,
        'MIN_RPC_HITS' : 1,
        # 'MAX_N_JETS' : 15,
        # 'MAX_N_LEPS' : 5,
        'MAX_CSC_JET' : 100,
        'MAX_DT_JET' : 100,
        'MAX_CSC_MUON' : 100,
        'MAX_DT_MUON' : 100,
        'MAX_ME1' : 0,
        'MAX_MB1' : 10,
        'HALO_CUTOFF' : 0.2,
        'MIN_DPHI' : 0,
        'MIN_DETA' : 0,
        'MAX_DETA' : 4,
        'MIN_CSC_DNN' : 0,
        'MIN_DT_DNN' : 0,
        'MIN_DT_SIZE' : 50,
        # 'ABCD_PHI' : 75,
        # 'ABCD_SIZE' : 75,
        # # Randomly choosing a value that has a better S/sqrt[B], frac=1 niter=100, S=622.8 B=4, SR S/rt[B]=558
        # 'MIN_CSC_TIME' : -5.00,
        # 'MAX_CSC_TIME' : 12.50,
        # 'MAX_CSC_TSPREAD' : 20.00,
        # 'MAX_RPC_BX' : 0,
        # 'MIN_RPC_HITS' : 1,
        # 'MAX_CSC_JET' : 72,
        # 'MAX_DT_JET' : 199,
        # 'MAX_CSC_MUON' : 146,
        # 'MAX_DT_MUON' : 200,
        # 'MAX_ME1' : 0,
        # 'MAX_MB1' : 8,
        # 'HALO_CUTOFF' : 0.08,
        # 'MIN_DPHI' : 0.40,
        # 'MIN_DETA' : 0,
        # 'MAX_DETA' : 4,
        # 'MIN_CSC_DNN' : 0.82,
        # 'MIN_DT_DNN' : 0.81,
        # 'MIN_DT_SIZE' : 50,
        # # Randomly choosing a value that has a better S/sqrt[B], frac=0.5 niter=200, S=879.0 B=11, SR S/rt[B]=320
        # 'MIN_CSC_TIME' : -5.00,
        # 'MAX_CSC_TIME' : 12.50,
        # 'MAX_CSC_TSPREAD' : 20.00,
        # 'MAX_RPC_BX' : 0,
        # 'MIN_RPC_HITS' : 1,
        # 'MAX_CSC_JET' : 101,
        # 'MAX_DT_JET' : 177,
        # 'MAX_CSC_MUON' : 200,
        # 'MAX_DT_MUON' : 200,
        # 'MAX_ME1' : 0,
        # 'MAX_MB1' : 2,
        # 'HALO_CUTOFF' : 0.04,
        # 'MIN_DPHI' : 0.40,
        # 'MIN_DETA' : 0,
        # 'MAX_DETA' : 4,
        # 'MIN_CSC_DNN' : 0.70,
        # 'MIN_DT_DNN' : 0.63,
        # 'MIN_DT_SIZE' : 50,
        # # Randomly choosing a value that has a better S/sqrt[B], wmv/4 in avg, niter=500, S=783.5 B=8, SR S/rt[B]=320
        # 'MIN_CSC_TIME' : -5.00,
        # 'MAX_CSC_TIME' : 12.50,
        # 'MAX_CSC_TSPREAD' : 20.00,
        # 'MAX_RPC_BX' : 0,
        # 'MIN_RPC_HITS' : 1,
        # 'MAX_CSC_JET' : 123,
        # 'MAX_DT_JET' : 129,
        # 'MAX_CSC_MUON' : 159,
        # 'MAX_DT_MUON' : 171,
        # 'MAX_ME1' : 0,
        # 'MAX_MB1' : 9,
        # 'HALO_CUTOFF' : 0.14,
        # 'MIN_DPHI' : 0.54,
        # 'MIN_DETA' : 0.2,
        # 'MAX_DETA' : 3.2,
        # 'MIN_CSC_DNN' : 0.72,
        # 'MIN_DT_DNN' : 0.71,
        # 'MIN_DT_SIZE' : 50,
        #
        # 'MIN_CSC_TIME' : -5.00,
        # 'MAX_CSC_TIME' : 12.50,
        # 'MAX_CSC_TSPREAD' : 20.00,
        # 'MAX_RPC_BX' : 0,
        # 'MIN_RPC_HITS' : 1,
        # 'MAX_CSC_JET' : 57,
        # 'MAX_DT_JET' : 138,
        # 'MAX_CSC_MUON' : 146,
        # 'MAX_DT_MUON' : 158,
        # 'MAX_ME1' : 0,
        # 'MAX_MB1' : 9,
        # 'HALO_CUTOFF' : 0.04,
        # 'MIN_DPHI' : 0.38,
        # 'MIN_DETA' : 0.0,
        # 'MAX_DETA' : 2.4,
        # 'MIN_CSC_DNN' : 0.69,
        # 'MIN_DT_DNN' : 0.82,
        # 'MIN_DT_SIZE' : 50,#53,
    },
    'lopt_high' : {
        'MIN_CSC_TIME' : -5.0,
        'MAX_CSC_TIME' : 12.5,
        'MAX_CSC_TSPREAD' : 20.0,
        'MAX_RPC_BX' : 0,
        'MIN_RPC_HITS' : 1,
        # 'MAX_N_JETS' : 15,
        # 'MAX_N_LEPS' : 5,
        'MAX_CSC_JET' : 100,
        'MAX_DT_JET' : 100,
        'MAX_CSC_MUON' : 100,
        'MAX_DT_MUON' : 100,
        'MAX_ME1' : 0,
        'MAX_MB1' : 10,
        'HALO_CUTOFF' : 0.2,
        'MIN_DPHI' : 0,
        'MIN_DETA' : 0,
        'MAX_DETA' : 4,
        'MIN_CSC_DNN' : 0,
        'MIN_DT_DNN' : 0,
        'MIN_DT_SIZE' : 50,
    },
    'loptDNN_high' : {
        'MIN_CSC_TIME' : -5.0,
        'MAX_CSC_TIME' : 12.5,
        'MAX_CSC_TSPREAD' : 20.0,
        'MAX_RPC_BX' : 0,
        'MIN_RPC_HITS' : 1,
        # 'MAX_N_JETS' : 15,
        # 'MAX_N_LEPS' : 5,
        'MAX_CSC_JET' : 100,
        'MAX_DT_JET' : 100,
        'MAX_CSC_MUON' : 100,
        'MAX_DT_MUON' : 100,
        'MAX_ME1' : 0,
        'MAX_MB1' : 10,
        'HALO_CUTOFF' : 0.2,
        'MIN_DPHI' : 0,
        'MIN_DETA' : 0,
        'MAX_DETA' : 4,
        'MIN_CSC_DNN' : 0,
        'MIN_DT_DNN' : 0,
        'MIN_DT_SIZE' : 50,
        # 'ABCD_PHI' : 1.5,
        # 'ABCD_SIZE' : 75,
        # # Randomly choosing a value that has a better S/sqrt[B], wmv/4 in avg, niter=500, S=783.5 B=8, SR S/rt[B]=320
        # 'MIN_CSC_TIME' : -5.00,
        # 'MAX_CSC_TIME' : 12.50,
        # 'MAX_CSC_TSPREAD' : 20.00,
        # 'MAX_RPC_BX' : 0,
        # 'MIN_RPC_HITS' : 1,
        # 'MAX_CSC_JET' : 188,
        # 'MAX_DT_JET' : 115,
        # 'MAX_CSC_MUON' : 200,
        # 'MAX_DT_MUON' : 200,
        # 'MAX_ME1' : 0,
        # 'MAX_MB1' : 8,
        # 'HALO_CUTOFF' : 0.04,
        # 'MIN_DPHI' : 0.24,
        # 'MIN_DETA' : 0.1,
        # 'MAX_DETA' : 3.0,
        # 'MIN_CSC_DNN' : 0.21,
        # 'MIN_DT_DNN' : 0.01,
        # 'MIN_DT_SIZE' : 50,
        # 'MIN_CSC_TIME' : -5.00,
        # 'MAX_CSC_TIME' : 12.50,
        # 'MAX_CSC_TSPREAD' : 20.00,
        # 'MAX_RPC_BX' : 0,
        # 'MIN_RPC_HITS' : 1,
        # 'MAX_CSC_JET' : 118,
        # 'MAX_DT_JET' : 132,
        # 'MAX_CSC_MUON' : 97,
        # 'MAX_DT_MUON' : 183,
        # 'MAX_ME1' : 0,
        # 'MAX_MB1' : 10,
        # 'HALO_CUTOFF' : 0.02,
        # 'MIN_DPHI' : 0.12,
        # 'MIN_DETA' : 0.1,
        # 'MAX_DETA' : 2.5,
        # 'MIN_CSC_DNN' : 0.33,
        # 'MIN_DT_DNN' : 0.06,
        # 'MIN_DT_SIZE' : 50,


    },
}

LOO_CUTS = ['MB1','CSC jet veto','DT jet veto','CSC muon veto','DT muon veto','halo veto','min dPhi','min dEta','max dEta','CSC DNN','DT DNN','DT size']#,'ABCD phi', 'ABCD size',]#,'ME1']
CUT_OPT_PARS = {
    'ME1' : {
        'col' : 'cscNHitME1',
        'cut' : 'MAX_ME1',
        'values' : np.arange(0, 11, 1)[::-1],
        'func' : lambda v, vd: vd <= v,
    },
    'MB1' : {
        'col' : 'dtNHitStation1',
        'cut' : 'MAX_MB1',
        'values' : np.arange(0, 11, 1)[::-1],
        'func' : lambda v, vd: vd <= v,
    },
    'CSC jet veto' : {
        'col' : 'cscJetVetoPt',
        'cut' : 'MAX_CSC_JET',
        'values' : np.arange(10, 201, 1)[::-1],
        # 'func' : lambda v, vd: (10 < vd) & (vd < v),
        'func' : lambda v, vd: vd < v,
    },
    'DT jet veto' : {
        'col' : 'dtJetVetoPt',
        'cut' : 'MAX_DT_JET',
        'values' : np.arange(10, 201, 1)[::-1],
        # 'func' : lambda v, vd: (10 < vd) & (vd < v),
        'func' : lambda v, vd: vd < v,
    },
    'CSC muon veto' : {
        'col' : 'cscMuonVetoPt',
        'cut' : 'MAX_CSC_MUON',
        'values' : np.arange(0, 201, 1)[::-1],
        # 'func' : lambda v, vd: (0 < vd) & (vd < v),
        'func' : lambda v, vd: vd < v,
    },
    'DT muon veto' : {
        'col' : 'dtMuonVetoPt',
        'cut' : 'MAX_DT_MUON',
        'values' : np.arange(0, 201, 1)[::-1],
        # 'func' : lambda v, vd: (0 < vd) & (vd < v),
        'func' : lambda v, vd: vd < v,
    },
    'halo veto' : {
        'col' : 'dtPhi',
        'cut' : 'HALO_CUTOFF',
        'values' : np.arange(0, 0.52, 0.02),
        'func' : lambda v, vd: (v < np.abs(vd)) & (np.abs(vd) < PI - v),
    },
    'min dPhi' : {
        'col' : 'tag_dPhi',
        'cut' : 'MIN_DPHI',
        'values' : np.arange(0, 0.62, 0.02),
        'func' : lambda v, vd: v < vd,
    },
    'min dEta' : {
        'col' : 'tag_dEta',
        'cut' : 'MIN_DETA',
        'values' : np.arange(0, 4.1, 0.1),
        'func' : lambda v, vd: v < vd,
    },
    'max dEta' : {
        'col' : 'tag_dEta',
        'cut' : 'MAX_DETA',
        'values' : np.arange(0, 4.1, 0.1)[::-1],
        'func' : lambda v, vd: vd < v,
    },
    'CSC DNN' : {
        'col' : 'cscDNN',
        'cut' : 'MIN_CSC_DNN',
        'values' : np.arange(0, 1.01, 0.01),
        'func' : lambda v, vd: v < vd,
    },
    'DT DNN' : {
        'col' : 'dtDNN',
        'cut' : 'MIN_DT_DNN',
        'values' : np.arange(0, 1.01, 0.01),
        'func' : lambda v, vd: v < vd,
    },
    'CSC DNN' : {
        'col' : 'cscDNN',
        'cut' : 'MIN_CSC_DNN',
        'values' : np.arange(0, 1.01, 0.01),
        'func' : lambda v, vd: v < vd,
    },
    'DT DNN' : {
        'col' : 'dtDNN',
        'cut' : 'MIN_DT_DNN',
        'values' : np.arange(0, 1.01, 0.01),
        'func' : lambda v, vd: v < vd,
    },
    'ABCD phi' : {
        'col' : 'tag_dPhi',
        'cut' : 'ABCD_PHI',
        'values' : np.arange(0, 3.5, 0.01),
        'func' : lambda v, vd: v <= vd,
    },
    'ABCD size' : {
        'col' : 'dtSize',
        'cut' : 'ABCD_SIZE',
        'values' : np.arange(50, 150, 1),
        'func' : lambda v, vd: v <= vd,
    },
    'DT size' : {
        'col' : 'dtSize',
        'cut' : 'MIN_DT_SIZE',
        'values' : np.arange(50, 150, 1),
        'func' : lambda v, vd: v <= vd,
    },
}

# **** #
C, D = 'cscRechitCluster', 'dtRechitCluster'

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
    'cscRechitClusterAvgStation10',
    'cscRechitClusterEta',
    'cscRechitClusterJetVetoE',
    'cscRechitClusterJetVetoLooseId',
    'cscRechitClusterJetVetoPt',
    'cscRechitClusterJetVetoTightId',
    'cscRechitClusterMaxChamber',
    'cscRechitClusterMaxChamberRatio',
    'cscRechitClusterMaxStation',
    'cscRechitClusterMaxStationRatio',
    'cscRechitClusterMet_dPhi',
    'cscRechitClusterMuonVetoE',
    'cscRechitClusterMuonVetoGlobal',
    'cscRechitClusterMuonVetoLooseId',
    'cscRechitClusterMuonVetoPt',
    'cscRechitClusterNChamber',
    'cscRechitClusterNRechitChamberMinus11',
    'cscRechitClusterNRechitChamberMinus12',
    'cscRechitClusterNRechitChamberMinus13',
    'cscRechitClusterNRechitChamberMinus21',
    'cscRechitClusterNRechitChamberMinus22',
    'cscRechitClusterNRechitChamberMinus31',
    'cscRechitClusterNRechitChamberMinus32',
    'cscRechitClusterNRechitChamberMinus41',
    'cscRechitClusterNRechitChamberMinus42',
    'cscRechitClusterNRechitChamberPlus11',
    'cscRechitClusterNRechitChamberPlus12',
    'cscRechitClusterNRechitChamberPlus13',
    'cscRechitClusterNRechitChamberPlus21',
    'cscRechitClusterNRechitChamberPlus22',
    'cscRechitClusterNRechitChamberPlus31',
    'cscRechitClusterNRechitChamberPlus32',
    'cscRechitClusterNRechitChamberPlus41',
    'cscRechitClusterNRechitChamberPlus42',
    'cscRechitClusterNStation10',
    'cscRechitClusterPhi',
    'cscRechitClusterSize',
    'cscRechitClusterTime',
    'cscRechitClusterTimeSpread',
    'cscRechitClusterTimeSpreadWeightedAll',
    'cscRechitClusterTimeWeighted',
    'cscRechitClusterX',
    'cscRechitClusterY',
    'cscRechitClusterZ',
    'cscRechitCluster_match_MB1Seg_0p4',
    'cscRechitCluster_match_RB1_0p4',
    'cscRechitCluster_match_RE12_0p4',
    'cscRechitCluster_match_dtSeg_0p4',
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
    'dtRechitClusterAvgStation10',
    'dtRechitClusterEta',
    'dtRechitClusterJetVetoE',
    'dtRechitClusterJetVetoLooseId',
    'dtRechitClusterJetVetoPt',
    'dtRechitClusterJetVetoTightId',
    'dtRechitClusterMaxChamber',
    'dtRechitClusterMaxChamberRatio',
    'dtRechitClusterMaxStation',
    'dtRechitClusterMaxStationRatio',
    'dtRechitClusterMet_dPhi',
    'dtRechitClusterMuonVetoE',
    'dtRechitClusterMuonVetoGlobal',
    'dtRechitClusterMuonVetoLooseId',
    'dtRechitClusterMuonVetoPt',
    'dtRechitClusterMuonVetoTightId',
    'dtRechitClusterNChamber',
    'dtRechitClusterNHitStation1',
    'dtRechitClusterNHitStation2',
    'dtRechitClusterNHitStation3',
    'dtRechitClusterNHitStation4',
    'dtRechitClusterNOppositeSegStation1',
    'dtRechitClusterNOppositeSegStation2',
    'dtRechitClusterNOppositeSegStation3',
    'dtRechitClusterNOppositeSegStation4',
    'dtRechitClusterNSegStation1',
    'dtRechitClusterNSegStation2',
    'dtRechitClusterNSegStation3',
    'dtRechitClusterNSegStation4',
    'dtRechitClusterNStation10',
    'dtRechitClusterNoiseHit',
    'dtRechitClusterNoiseHitStation1',
    'dtRechitClusterNoiseHitStation2',
    'dtRechitClusterNoiseHitStation3',
    'dtRechitClusterNoiseHitStation4',
    'dtRechitClusterOverlap',
    'dtRechitClusterPhi',
    'dtRechitClusterSize',
    'dtRechitClusterWheel',
    'dtRechitClusterX',
    'dtRechitClusterY',
    'dtRechitClusterZ',
    'dtRechitCluster_match_MB1Seg_0p4',
    'dtRechitCluster_match_MB1Seg_0p5',
    'dtRechitCluster_match_MB1hits_0p4',
    'dtRechitCluster_match_MB1hits_0p5',
    'dtRechitCluster_match_MB1hits_cosmics_minus',
    'dtRechitCluster_match_MB1hits_cosmics_plus',
    'dtRechitCluster_match_RB1_0p4',
    'dtRechitCluster_match_RB1_dPhi0p5',
    'dtRechitCluster_match_RPCBx_dPhi0p5',
    'dtRechitCluster_match_RPChits_dPhi0p5',
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
    'evtNum',
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
    'jetE',
    'jetEta',
    'jetPhi',
    'jetPt',
    'jetTightPassId',
    'lepDZ',
    'lepE',
    'lepEta',
    'lepPassLooseIso',
    'lepPassTightIso',
    'lepPassVTightIso',
    'lepPassVVTightIso',
    'lepPdgId',
    'lepPhi',
    'lepPt',
    'lepTightId',
    # 'lumiSec',
    # 'mH',
    # 'mX',
    'met',
    'metPhi',
    # 'nCscRechitClusters',
    # 'nCscRings',
    # 'nDtRechitClusters',
    # 'nDtRings',
    # 'nGLLP',
    'nJets',
    'nLeptons',
    # 'npu',
    # 'npv',
    # 'pileupWeight',
    # 'rho',
    # 'runNum',
    'weight',
]

# **** #
pathlib.Path(OUT_DIR).mkdir(parents=True, exist_ok=True)

rng = np.random.default_rng()

rt.gErrorIgnoreLevel = 1001  # rt.kInfo + 1
rt.gROOT.SetBatch(True)
PI = rt.TMath.Pi()

gc = [] # ROOT garbage collector

# **************************** #
def optimize_cut(_rdfs, col, cut_func, cut_values=None, cut0=None):
    vmc = _rdfs['mc'].AsNumpy(['weight',col])#,'tag_dPhi','dtSize'
    vr3 = _rdfs['r3'].AsNumpy(['weight',col])#,'tag_dPhi','dtSize'
    wmc, vmc = vmc['weight'], vmc[col]
    wr3, vr3 = vr3['weight'], vr3[col]
    # wmc_tot, wr3_tot = wmc.sum(), wr3.sum() # may have issues when cut_func has preselection?
    # vmin, vmax = min(vmc.min(), vr3.min()), max(vmc.max(), vr3.max())

    # **** #
    nmc, nr3 = _rdfs['mc'].Count().GetValue(), _rdfs['r3'].Count().GetValue()
    if isinstance(LOO_FRAC, float):
        idxs_mc = np.random.randint(0,nmc,int(LOO_FRAC*nmc))
        idxs_r3 = np.random.randint(0,nr3,int(LOO_FRAC*nr3))
    elif LOO_FRAC == 'sqrt':
        idxs_mc = np.random.randint(0,nmc,int(nmc**0.5))
        idxs_r3 = np.random.randint(0,nr3,int(nr3**0.5))
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

    wmc_tot, wr3_tot = wmc.max(), wr3.max() # relative to itself? sorta fixes preselection bias
    # **** #
    # Scoring Metrics
    # S/rtB, sometimes falls into minimum where B=1 and S=small
    s2b0 = wmc0/(wr30**0.5) if wr30 else 0
    s2b = np.divide(wmc, np.sqrt(wr3), out=np.zeros_like(wmc), where=(wr3>0))
    # # F1 to balance precision/recall (avoid small S minimum?)
    # f10 = wmc0 /(2*wmc0+wr30+(wmc_tot-wmc0)) if 2*wmc0+wr30+(wmc_tot-wmc0) else 0
    # f1 = np.divide(wmc, 2*wmc+wr3+(wmc_tot-wmc),out=np.zeros_like(wmc), where=2*wmc+wr3+(wmc_tot-wmc)>0)
    # # Modified F-score where precision = S / (S + rtB)
    # f10 = wmc0 /(2*wmc0+(wr30**0.5)+(wmc_tot-wmc0)) if 2*wmc0+(wr30**0.5)+(wmc_tot-wmc0) else 0
    # f1 = np.divide(wmc, 2*wmc+np.sqrt(wr3)+(wmc_tot-wmc),out=np.zeros_like(wmc), where=2*wmc+np.sqrt(wr3)+(wmc_tot-wmc)>0)
    # # Modified F-score where precision = S / rtB, harmonic mean of S/rtB and S recall
    # f10 = wmc0 /(wmc_tot+(wr30**0.5)) if wmc_tot+(wr30**0.5) else 0
    # f1 = np.divide(wmc, wmc_tot+np.sqrt(wr3),out=np.zeros_like(wmc), where=wmc_tot+np.sqrt(wr3)>0)

    # **** #
    score0, score = s2b0, s2b
    # score0, score = s2b0, s2b
    # score0, score = f10, f1

    # # Find best at this working point
    # idx = np.argmax(score)
    # Randomly select a new cut that's better than before
    idx = rng.choice(np.arange(len(score))[score >= score0]) if (score >= score0).sum() else 0
    
    # # Best signal yield while still being a better s2b
    # idx = np.argmax(wmc * (s2b >= s2b0))

    # idx = np.argmax(score * (s2b >= s2b0))
    # idx = np.argmax(s2b * (score >= score0))
    # idx = rng.choice(np.arange(len(score))[s2b >= s2b0]) if (s2b >= s2b0).sum() else 0

    # **** #
    if cut0 is not None:
        return cut_values[idx], s2b[idx], wmc[idx], wr3[idx], s2b0, wmc0, wr30
    return cut_values[idx], s2b[idx], wmc[idx], wr3[idx]

# **************************** #
if __name__ == '__main__':
    print('+------------------------+')
    print('| Starting skim_cscdt.py |')
    print('+------------------------+')

    rt.EnableImplicitMT(0)
    print('    Enabled ROOT\'s implicit multithreading (sometimes causes a crash)')

    # **************************** #
    N_ITERATIONS = 1

    LOO_SCORES = []
    LOO_FRAC = 1.0 # if 1 just bootstrap, fraction with replacement
    LOO = False
    
    PRINT_CUTFLOW = False
    MET_CATEGORY = 'lt200'
    CUTSET = 'scs'
    OOT = False

    args = ' '.join(sys.argv[1:]) if len(sys.argv) > 1 else ''

    if 'cutflow' in args:
        print('    Printing cutflow tables')
        PRINT_CUTFLOW = True

    if 'l1' in args:
        print('    Using the reduced cut set (l1)')
        CUTS = CUTS_L1
        CUTSET = 'l1'
    elif 'ropt' in args:
        print('    Using the randomly optimized cut set (ropt)')
        CUTSET = 'ropt'
    elif 'lopt' in args:
        print('    Using the leave-one-out optimized cut set (lopt)')
        CUTSET = 'lopt'
    else:
        print('    Using the standard cut set (scs)')
        CUTSET = 'scs'
    
    if 'low' in args:
        print('    Low met category')
        CUTS = [c.replace('MET', 'low MET') if 'MET' == c else c for c in CUTS]
        MET_CATEGORY = 'low'
    elif 'high' in args:
        print('    High met category')
        CUTS = [c.replace('MET', 'high MET') if 'MET' == c else c for c in CUTS]
        MET_CATEGORY = 'high'
    else:
        print('    No met categorization (only met<200)')
        MET_CATEGORY = 'lt200'

    if 'oot' in args:
        print('    Using out-of-time 2nd cluster')
        CUTS = [c.replace('DT IT', 'DT OOT') if 'DT IT' == c else c for c in CUTS]
        OOT = True
    else:
        print('    Using in-time 2nd cluster')

    if 'loo' in args:
        print('    PERFORMING LOO OPTIMIZATION')
        if OOT:
            print('        FORCING MC TO IN-TIME')
        LOO = True
        N_ITERATIONS = 5001
    else:
        CUTS = [c for c in CUTS if 'DT size' not in c]

    if 'DNN' not in CUTSET:
        print('    Removing DNN from CUTS')
        CUTS = [c for c in CUTS if 'DNN' not in c]
        LOO_CUTS = [c for c in LOO_CUTS if 'DNN' not in c]

    if f'{CUTSET}_{MET_CATEGORY}' in CUT_VALUES:
        CUT_VALUES = CUT_VALUES[f'{CUTSET}_{MET_CATEGORY}']
    elif f'{CUTSET}_lt200' in CUT_VALUES:
        CUT_VALUES = CUT_VALUES[f'{CUTSET}_lt200']
        print(f'    No cutset for {MET_CATEGORY=} found, defaulting to \'lt200\'')
    else:
        raise ValueError(f'no cutset associated with {CUTSET=} and {MET_CATEGORY=} found')

    print('')

    # **** #
    if '1 CSC-DT' not in CUTS:
        raise NotImplementedError('cant handle multiple pairs yet')

    # **************************** #
    for iloo in range(N_ITERATIONS):
        # if iloo % 100>0 and iloo != N_ITERATIONS-1:
        if iloo != 0 and iloo != N_ITERATIONS-1:
            LOO_CUT = rng.choice([c for c in LOO_CUTS if c != LOO_CUT])
        else:
            LOO_CUT = ''

        MIN_CSC_TIME = CUT_VALUES['MIN_CSC_TIME']
        MAX_CSC_TIME = CUT_VALUES['MAX_CSC_TIME']
        MAX_CSC_TSPREAD = CUT_VALUES['MAX_CSC_TSPREAD']
        MAX_RPC_BX = CUT_VALUES['MAX_RPC_BX']
        MIN_RPC_HITS = CUT_VALUES['MIN_RPC_HITS']
        MAX_CSC_JET = CUT_VALUES['MAX_CSC_JET']
        MAX_DT_JET = CUT_VALUES['MAX_DT_JET']
        MAX_CSC_MUON = CUT_VALUES['MAX_CSC_MUON']
        MAX_DT_MUON = CUT_VALUES['MAX_DT_MUON']
        MAX_ME1 = CUT_VALUES['MAX_ME1']
        MAX_MB1 = CUT_VALUES['MAX_MB1']
        HALO_CUTOFF = CUT_VALUES['HALO_CUTOFF']
        MIN_DPHI = CUT_VALUES['MIN_DPHI']
        MIN_DETA = CUT_VALUES['MIN_DETA']
        MAX_DETA = CUT_VALUES['MAX_DETA']
        MIN_CSC_DNN = CUT_VALUES['MIN_CSC_DNN']
        MIN_DT_DNN = CUT_VALUES['MIN_DT_DNN']
        # ABCD_DPHI = CUT_VALUES['ABCD_DPHI']
        # ABCD_SIZE = CUT_VALUES['ABCD_SIZE']
        MIN_DT_SIZE = CUT_VALUES['MIN_DT_SIZE']

        rdfs = {
            'mc' : RDataFrame('MuonSystem', FN_MC),
            'r3' : RDataFrame('MuonSystem', FN_R3),
        }

        if iloo == 0:
            print('Events in:')
        for key, rdf in rdfs.items():
            # rdf = rdf.Range(0,100_000) # Skim a subset of events for debugging

            if key == 'mc': # fix weights
                rdf = rdf.Redefine('weight', f'weight * {LUMI}')

            count, weight = rdf.Count().GetValue(), rdf.Sum('weight').GetValue()
            if iloo == 0:
                print(f'  {key} = {count:,} ({weight:,.2f}) -- read')

            # **** #
            # Create dummy columns to store what cluster indices pass our selections
            rdf = rdf.Define('evtFlag', 'weight > 0')
            rdf = rdf.Define(f'{C}Flag', f'{C}Size > 0')
            rdf = rdf.Define(f'{D}Flag', f'{D}Size > 0')

            rdf = rdf.Define('evtCutFlag', 'weight > 0')
            rdf = rdf.Define(f'{C}CutFlag', f'{C}Size > 0')
            rdf = rdf.Define(f'{D}CutFlag', f'{D}Size > 0')

            # **** #
            # Add cluster radius columns
            for p in (C, D):
                rdf = rdf.Define(f'{p}R', f'ROOT::VecOps::sqrt( {p}X*{p}X + {p}Y*{p}Y )')
                if f'{p}R' not in COLUMNS_OUT:
                    COLUMNS_OUT.append(f'{p}R')
            
            rdf = rdf.Define(f'{C}NHitME1', f'{C}NRechitChamberPlus11 + {C}NRechitChamberPlus12 + {C}NRechitChamberMinus11 + {C}NRechitChamberMinus12')
            if f'{C}NHitME1' not in COLUMNS_OUT:
                COLUMNS_OUT.append(f'{C}NHitME1')

            # **** #
            rdfs[key] = rdf

        if iloo == 0:
            print('')

        # **************** #
        # Apply selections from CUTS (and/or arguments?)
        for key, rdf in rdfs.items():
            if PRINT_CUTFLOW:
                print(f'{key} cutflow:')
                print(r'\begin{center}')
                print(r'\begin{tabular}{c|ccc|ccc|ccc}')
                ec0,cc0,dc0 = rdf.Filter('evtFlag').Sum('weight'), rdf.Filter('evtFlag').Sum(f'{C}Flag'), rdf.Filter('evtFlag').Sum(f'{D}Flag')
                ec0,cc0,dc0 = ec0.GetValue(),cc0.GetValue(),dc0.GetValue()
                print(r'    \hline')
                print(f'    \\textbf{{{"Signal" if "mc" in key else "Data"}}}'+r' & \multicolumn{3}{c}{CSC-DT} & \multicolumn{3}{c}{CSC} & \multicolumn{3}{c}{DT} \\')
                print(r'    \hline')
                print(r'    Selection & N Events & Cut Eff & Cum Eff & N Clusters & Cut Eff & Cum Eff & N Clusters & Cut Eff & Cum Eff \\')
                print(r'    \hline')
                print(f'    sample & {ec0:,.0f} & -- & -- & {cc0:,} & -- & -- & {dc0:,} & -- & -- \\\\')

            for cut in CUTS:
                # **** #
                # Reset cut flags
                rdf = rdf.Redefine('evtCutFlag', 'weight > 0')
                rdf = rdf.Redefine(f'{C}CutFlag', f'{C}Size > 0')
                rdf = rdf.Redefine(f'{D}CutFlag', f'{D}Size > 0')

                # **** #
                # Signal matching for MC
                if 'acceptance' in cut: 
                    if 'mc' in key:
                        rdf = rdf.Redefine(f'{C}CutFlag', f'{C}CutFlag && {C}_match_gLLP && '+
                                        f'(abs({C}_match_gLLP_eta) < 3) && '+
                                        f'({C}_match_gLLP_decay_r < 800) && '+
                                        f'(400 < abs({C}_match_gLLP_decay_z)) &&'+
                                        f' (abs({C}_match_gLLP_decay_z) < 1200)')

                        rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && {D}_match_gLLP && '+
                                        f'(200 < {D}_match_gLLP_decay_r) && '+
                                        f'({D}_match_gLLP_decay_r < 800) && '+
                                        f'(abs({D}_match_gLLP_decay_z) < 700)')
                    elif 'r3' in key:
                        continue

                # **** #
                # Trigger selections (HLT may be wrong in MC?)
                if 'HLT' in cut:
                    if 'mc' in key and 'hlt' not in FN_MC:
                        continue
                    elif 'r3' in key and 'hlt' not in FN_R3:
                        rdf = rdf.Redefine('evtCutFlag', 'evtCutFlag && HLTDecision[569]') # 569 = HLT_L1CSCCluster_DTCluster50
                    else:
                        continue

                if 'L1' in cut:
                    rdf = rdf.Redefine(f'{C}CutFlag', f'{C}CutFlag && ( '+
                                    f'((100 < {C}R) && ({C}R < 275) && (580 < abs({C}Z)) && (abs({C}Z) < 632) && (500 <= {C}Size)) || ' + # ME 11
                                    f'((275 < {C}R) && ({C}R < 465) && (668 < abs({C}Z)) && (abs({C}Z) < 724) && (200 <= {C}Size)) || ' + # ME 12
                                    f'((505 < {C}R) && ({C}R < 700) && (668 < abs({C}Z)) && (abs({C}Z) < 724) && (200 <= {C}Size)) || ' + # ME 13
                                    #
                                    f'((139 < {C}R) && ({C}R < 345) && (789 < abs({C}Z)) && (abs({C}Z) < 850) && (500 <= {C}Size)) || ' + # ME 21
                                    f'((357 < {C}R) && ({C}R < 700) && (791 < abs({C}Z)) && (abs({C}Z) < 850) && (200 <= {C}Size)) || ' + # ME 22
                                    #
                                    f'((160 < {C}R) && ({C}R < 345) && (915 < abs({C}Z)) && (abs({C}Z) < 970) && (500 <= {C}Size)) || ' + # ME 31
                                    f'((357 < {C}R) && ({C}R < 700) && (911 < abs({C}Z)) && (abs({C}Z) < 970) && (200 <= {C}Size)) || ' + # ME 32
                                    #
                                    f'((178 < {C}R) && ({C}R < 345) && (1002 < abs({C}Z)) && (abs({C}Z) < 1063) && (500 <= {C}Size)) || ' + # ME 41
                                    f'((357 < {C}R) && ({C}R < 700) && (1002 < abs({C}Z)) && (abs({C}Z) < 1063) && (200 <= {C}Size)) )') # ME 42

                # **** #
                # Missing Et requirements & categorization
                if 'MET' in cut:
                    rdf = rdf.Redefine('evtCutFlag', f'evtCutFlag && (met < 200)')

                if 'low MET' in cut:
                    rdf = rdf.Redefine('evtCutFlag', f'evtCutFlag && (met < {LOW_MET_CUTOFF})')

                if 'high MET' in cut:
                    rdf = rdf.Redefine('evtCutFlag', f'evtCutFlag && (({HIGH_MET_CUTOFF} < met) && (met < 200))')

                # **** #
                # Cluster time requirements
                if LOO and 'mc' in key and ' OOT' in cut:
                    cut = cut.replace(' OOT', ' IT')

                if 'CSC IT' in cut:
                    rdf = rdf.Redefine(f'{C}CutFlag', f'{C}CutFlag && ( ({MIN_CSC_TIME} < {C}TimeWeighted) && ({C}TimeWeighted < {MAX_CSC_TIME}) && ({C}TimeSpreadWeightedAll < {MAX_CSC_TSPREAD}) )')
                elif 'CSC OOT' in cut:
                    rdf = rdf.Redefine(f'{C}CutFlag', f'{C}CutFlag && !( ({MIN_CSC_TIME} < {C}TimeWeighted) && ({C}TimeWeighted < {MAX_CSC_TIME}) && ({C}TimeSpreadWeightedAll < {MAX_CSC_TSPREAD}) )')

                if 'DT IT' in cut:
                    rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( (abs({D}_match_RPCBx_dPhi0p5) <= {MAX_RPC_BX}) && ({D}_match_RPChits_dPhi0p5 >= {MIN_RPC_HITS}) )')
                elif 'DT OOT' in cut:
                    rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && !( (abs({D}_match_RPCBx_dPhi0p5) <= {MAX_RPC_BX}) && ({D}_match_RPChits_dPhi0p5 >= {MIN_RPC_HITS}) )')

                # **** #
                # Station requirements from L1 trigger
                if 'ME1' in cut and LOO_CUT != 'ME1':
                    # rdf = rdf.Redefine(f'{C}CutFlag', f'{C}CutFlag && ( '+
                    #                 f'({C}NRechitChamberPlus11 <= {MAX_ME1}) &&'+
                    #                 f'({C}NRechitChamberPlus12 <= {MAX_ME1}) &&'+
                    #                 f'({C}NRechitChamberMinus11 <= {MAX_ME1}) &&'+
                    #                 f'({C}NRechitChamberMinus12 <= {MAX_ME1}) )')
                    rdf = rdf.Redefine(f'{C}CutFlag', f'{C}CutFlag && ( {C}NHitME1 <= {MAX_ME1} ) ')

                if 'MB1' in cut and LOO_CUT != 'MB1':
                    rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( {D}NHitStation1 <= {MAX_MB1} )')

                if 'DT stn' in cut:
                    rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( ({D}NStation10 < 3) && !(({D}NStation10 == 2) && ({D}MaxStation == 4)) )')

                # **** #
                if 'jet veto' in cut:
                    if LOO_CUT != 'CSC jet veto':
                        rdf = rdf.Redefine(f'{C}CutFlag', f'{C}CutFlag && ({C}JetVetoPt < {MAX_CSC_JET})')
                    if LOO_CUT != 'DT jet veto':
                        rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ({D}JetVetoPt < {MAX_DT_JET})')
                
                if 'muon veto' in cut:
                    rdf = rdf.Redefine(f'{C}CutFlag', f'{C}CutFlag && ({C}MuonVetoGlobal == 0)')
                    rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ({D}MuonVetoLooseId == 0)')

                    if LOO_CUT != 'CSC muon veto':
                        rdf = rdf.Redefine(f'{C}CutFlag', f'{C}CutFlag && ({C}MuonVetoPt < {MAX_CSC_MUON})')
                    if LOO_CUT != 'DT muon veto':
                        rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ({D}MuonVetoPt < {MAX_DT_MUON})')

                if 'halo veto' in cut and LOO_CUT != 'halo veto':
                    rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( ({HALO_CUTOFF} < abs({D}Phi)) && (abs({D}Phi) < {PI} - {HALO_CUTOFF}) )')

                # **** #
                if 'BDT' in cut:
                    raise NotImplementedError('BDT')

                if 'DNN' in cut:
                    if LOO_CUT != 'CSC DNN':
                        rdf = rdf.Redefine(f'{C}CutFlag', f'{C}CutFlag && ( Take({C}DNN,nCscRechitClusters) > {MIN_CSC_DNN} )')
                    if LOO_CUT != 'DT DNN':
                        rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( Take({D}DNN,nDtRechitClusters) > {MIN_DT_DNN} )')

                # **** #
                if '1 CSC-DT' in cut:
                    # rdf = rdf.Redefine('evtCutFlag', f'evtCutFlag && evtFlag && '+
                    #                    f'(reduce({C}Flag.begin(), {C}Flag.end()) == 1) && '+
                    #                    f'(reduce({D}Flag.begin(), {D}Flag.end()) == 1)')
                    if LOO:
                        icsc, idt = rng.random(2)
                        rdf = rdf.Redefine(f'{C}Flag', f'{C}Flag && ( {C}Size == {C}Size[int({icsc}*reduce({C}Flag.begin(),{C}Flag.end()))] )')
                        rdf = rdf.Redefine(f'{D}Flag', f'{D}Flag && ( {D}Size == {D}Size[int({idt}*reduce({D}Flag.begin(),{D}Flag.end()))] )')
                    else:
                        rdf = rdf.Redefine(f'{C}Flag', f'{C}Flag && ( {C}Size == Max({C}Size*{C}Flag) )')
                        rdf = rdf.Redefine(f'{D}Flag', f'{D}Flag && ( {D}Size == Max({D}Size*{D}Flag) )')

                    # Apply our cluster level selections to relevant columns
                    columns_out = []
                    for col in COLUMNS_OUT:
                        ocol = col.replace('RechitCluster', '')
                        if C in col[:len(C)]:
                            rdf = rdf.Define(ocol, f'{col}[{C}Flag][0]')
                        elif D in col[:len(D)]:
                            rdf = rdf.Define(ocol, f'{col}[{D}Flag][0]')
                        columns_out.append(ocol)

                    rdf = rdf.Define('cscDNN', f'Take({C}DNN,nCscRechitClusters)[{C}Flag][0]')
                    rdf = rdf.Define('dtDNN', f'Take({D}DNN,nDtRechitClusters)[{D}Flag][0]')
            
                    rdf = rdf.Define('cscCTau', f'gLLP_ctau[{C}_match_gLLP_index[{C}Flag][0]]')
                    rdf = rdf.Define('dtCTau', f'gLLP_ctau[{D}_match_gLLP_index[{D}Flag][0]]')
                    rdf = rdf.Define('tag_dEta', 'abs(cscEta - dtEta)')
                    rdf = rdf.Define('tag_dPhi', f'abs( abs(cscPhi - dtPhi) - (abs(cscPhi - dtPhi) > {PI} ? 2*{PI} : 0) )')
                    rdf = rdf.Define('tag_dR', 'sqrt(tag_dEta*tag_dEta + tag_dPhi*tag_dPhi)')
                    columns_out.extend(['tag_dEta', 'tag_dPhi', 'tag_dR', 'cscCTau', 'dtCTau', 'cscDNN', 'dtDNN'])

                # **** #
                if 'dR' in cut:
                    raise NotImplementedError('dR')

                if 'dEta' in cut:
                    if LOO_CUT != 'min dEta':
                        rdf = rdf.Redefine('evtCutFlag', f'evtCutFlag && (tag_dEta > {MIN_DETA})')
                    if LOO_CUT != 'max dEta':
                        rdf = rdf.Redefine('evtCutFlag', f'evtCutFlag && (tag_dEta < {MAX_DETA})')

                if 'dPhi' in cut and LOO_CUT != 'min dPhi':
                    rdf = rdf.Redefine('evtCutFlag', f'evtCutFlag && (tag_dPhi > {MIN_DPHI})')

                # **** #
                if 'DT size' in cut and LOO_CUT != 'DT size': # for ABCD?
                    rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( {MIN_DT_SIZE} <= {D}Size )')

                # **** #
                # Propagate to cumulative flags
                rdf = rdf.Redefine(f'{C}Flag', f'{C}Flag && {C}CutFlag')
                rdf = rdf.Redefine(f'{D}Flag', f'{D}Flag && {D}CutFlag')

                rdf = rdf.Redefine('evtFlag', 'evtFlag && evtCutFlag && ('+
                                f'(reduce({C}Flag.begin(), {C}Flag.end()) > 0) && '+
                                f'(reduce({D}Flag.begin(), {D}Flag.end()) > 0))')

                # rdf = rdf.Redefine(f'nC{C[1:]}s', f'reduce({C}Flag.begin(), {C}Flag.end())')
                # rdf = rdf.Redefine(f'nD{D[1:]}s', f'reduce({D}Flag.begin(), {D}Flag.end())')

                # **** #
                if PRINT_CUTFLOW:
                    rdf = rdf.Redefine('evtCutFlag', 'evtCutFlag && ('+
                                    f'(reduce({C}CutFlag.begin(), {C}CutFlag.end()) > 0) && '+
                                    f'(reduce({D}CutFlag.begin(), {D}CutFlag.end()) > 0))')

                    ec,cc,dc = rdf.Filter('evtFlag').Sum('weight'), rdf.Filter('evtFlag').Sum(f'{C}Flag'), rdf.Filter('evtFlag').Sum(f'{D}Flag')
                    # ecc,ccc,dcc = rdf.Filter('evtCutFlag').Sum('weight'), rdf.Filter('evtCutFlag').Sum(f'{C}CutFlag'), rdf.Filter('evtCutFlag').Sum(f'{D}CutFlag')
                    ecc,ccc,dcc = rdf.Filter('evtCutFlag').Sum('weight'), rdf.Sum(f'{C}CutFlag'), rdf.Sum(f'{D}CutFlag')
                    
                    ec,cc,dc = ec.GetValue(),cc.GetValue(),dc.GetValue()
                    ecc,ccc,dcc = ecc.GetValue(),ccc.GetValue(),dcc.GetValue()

                    c_cum_eff = f'{100*cc/cc0:.2f}' if cc != cc0 else '--'
                    c_cut_eff = f'{100*ccc/cc0:.2f}' if ccc != cc0 else '--'
                    d_cum_eff = f'{100*dc/dc0:.2f}' if dc != dc0 else '--'
                    d_cut_eff = f'{100*dcc/dc0:.2f}' if dcc != dc0 else '--'

                    print(f'    {cut.replace("_"," ")} & {ec:,.0f} & {100*ecc/ec0:.2f} & {100*ec/ec0:.2f} & '+
                        f'{cc:,.0f} & {c_cut_eff} & {c_cum_eff} & '+
                        f'{dc:,.0f} & {d_cut_eff} & {d_cum_eff} \\\\')

                    # Only filter (reset indices) after the MET cut when printing cutflow
                    if 'MET' in cut:
                        print(r'    \hline')
                        ec0, cc0, dc0 = ec, cc, dc
                        
                        rdf = rdf.Filter('evtFlag')
                        rdf = rdf.Redefine(f'{C}Size', f'{C}Size * {C}Flag')
                        rdf = rdf.Redefine(f'{D}Size', f'{D}Size * {D}Flag')
                else:
                    rdf = rdf.Filter('evtFlag')
                    rdf = rdf.Redefine(f'{C}Size', f'{C}Size * {C}Flag')
                    rdf = rdf.Redefine(f'{D}Size', f'{D}Size * {D}Flag')

            rdfs[key] = rdf.Filter('evtFlag') # Apply the final filter and update the dictionary with the new RDF

            if PRINT_CUTFLOW:
                print(r'    \hline')
                print(r'\end{tabular}')
                print(r'\end{center}')
                print('')

        if LOO:
            _c = False
            if LOO_CUT == '':
                val, wmc, wr3 = 0, rdfs['mc'].Sum('weight').GetValue(), rdfs['r3'].Sum('weight').GetValue()
                s2b = wmc/(wr3**0.5) if wr3 else 0
                val0, s2b0, wmc0, wr30 = 0, s2b, wmc, wr3
            else:
                val0 = CUT_VALUES[CUT_OPT_PARS[LOO_CUT]['cut']]
                val, s2b, wmc, wr3, s2b0, wmc0, wr30 = optimize_cut(
                    _rdfs=rdfs,
                    col=CUT_OPT_PARS[LOO_CUT]['col'],
                    cut_func=CUT_OPT_PARS[LOO_CUT]['func'],
                    cut_values=CUT_OPT_PARS[LOO_CUT]['values'],
                    cut0=CUT_VALUES[CUT_OPT_PARS[LOO_CUT]['cut']],
                )
                if wmc>0 and wr3>0:# and s2b > s2b0:#and wmc > wmc0:
                    _c = True
                    # CUT_VALUES[CUT_OPT_PARS[LOO_CUT]['cut']] = val
                    # sometimes I fall into a local minimum where the nbkg~1 so I get large s2b but a low signal acceptance
                    CUT_VALUES[CUT_OPT_PARS[LOO_CUT]['cut']] = (val*wmc/(10*(1-iloo/N_ITERATIONS)) + val0*wmc0)/(wmc/(10*(1-iloo/N_ITERATIONS)) + wmc0)
                    # CUT_VALUES[CUT_OPT_PARS[LOO_CUT]['cut']] = (s2b*val + s2b0*val0)/(s2b + s2b0)

            # print(f'{iloo:>3} | {"Y" if _c else "X"} | {LOO_CUT:>13} = {val:>6.2f} | {s2b:>3.0f}, {wmc:>6.0f}, {wr3:>6.0f}')
            # print(f'{iloo:>3} | {"Y" if _c else "X"} | {LOO_CUT:>13} = {val:>6.2f} ({val0:>6.2f}) | {s2b:>3.0f} ({s2b0:>3.0f}), {wmc:>6.0f} ({wmc0:>6.0f}), {wr3:>6.0f} ({wr30:>6.0f})')
            print(f'{iloo:>3} | {LOO_CUT:>13} = {val:>6.2f} ({val0:>6.2f}) | {s2b:>3.0f} ({s2b0:>3.0f}), {wmc:>6.0f} ({wmc0:>6.0f}), {wr3:>6.0f} ({wr30:>6.0f})')

            LOO_SCORES.append([LOO_CUT, val, s2b, wmc, wr3, {k:v for k,v in CUT_VALUES.items()}])
            with open(f'loo_scores_cscdt{"OOT" if OOT else ""}_{CUTSET}_{MET_CATEGORY}.pkl', 'wb') as fout:
                pickle.dump(LOO_SCORES, fout)

    # **************** #
    if LOO:
        print('\nLeave-One-Out Optimized Cuts:')
        print(f'    {MIN_CSC_TIME=:.2f}')
        print(f'    {MAX_CSC_TIME=:.2f}')
        print(f'    {MAX_CSC_TSPREAD=:.2f}')
        print(f'    {MAX_RPC_BX=:.0f}')
        print(f'    {MIN_RPC_HITS=:.0f}')
        print(f'    {MAX_CSC_JET=:.0f}')
        print(f'    {MAX_DT_JET=:.0f}')
        print(f'    {MAX_CSC_MUON=:.0f}')
        print(f'    {MAX_DT_MUON=:.0f}')
        print(f'    {MAX_ME1=:.0f}')
        print(f'    {MAX_MB1=:.0f}')
        print(f'    {HALO_CUTOFF=:.2f}')
        print(f'    {MIN_DPHI=:.2f}')
        print(f'    {MIN_DETA=:.1f}')
        print(f'    {MAX_DETA=:.1f}')
        print(f'    {MIN_CSC_DNN=:.2f}')
        print(f'    {MIN_DT_DNN=:.2f}')
        print(f'    {MIN_DT_SIZE=:.2f}')
        print('')

    print('Events out:')
    hists = {}
    for key, rdf in rdfs.items():
        count, weight = rdf.Count(), rdf.Sum('weight')
        count, weight = count.GetValue(), weight.GetValue()
        if count == 0:
            print(f'{key} is empty')
            continue

        name = f'{key}_cscdt{"OOT" if OOT else ""}_{CUTSET}_{MET_CATEGORY}'

        if LOO: 
            name = name.replace(f'_{CUTSET}_',f'_{CUTSET}LOO_')
            if key == 'mc':
                name = name.replace('OOT','')

        rdf = rdf.Snapshot('MuonSystem_flat', f'data/processed/{name}_rdf.root', columns_out)
        rdfs[key] = rdf

        print(f'  {key} = {count:,} ({weight:,.2f})')
        for xx in ('met', 'cscSize', 'cscR', 'cscEta', 'cscPhi', 'dtSize', 'dtR', 'dtEta', 'dtPhi', 'tag_dEta','tag_dPhi','tag_dR'):
            hh = rdf.Histo1D(xx,'weight')
            hh.SetName(f'{key}_{xx}')
            hh.SetTitle(f'{key};{xx};count')

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

    print('')
