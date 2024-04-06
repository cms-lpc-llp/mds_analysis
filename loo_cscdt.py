import os
import sys
import pathlib
import argparse

# import math
import numpy as np
# import numba as nb
# from math import ceil, floor

# from src.helper_functions import alert, Table
# from src.histo_utilities import std_color_list as SCL

import ROOT as rt
from ROOT import RDataFrame
# from ROOT import TNtuple, TTree, TBranch, RDataFrame
# from ROOT import TCanvas, TLatex, TLegend, TLine, TBox
# from ROOT import TH1D, TH2D, TGraph, TGraphErrors, TGraphAsymmErrors

import pickle
# **************************** #
CUTS = [
    'acceptance',
    'HLT',
    'L1',
    'low MET',
    # 'high MET',
    #! I reset cutflow indices here
    'CSC IT',
    'DT IT',
    # 'DT OOT',
    # 'DT I/OOT',
    # 'MET',
    'DNN',
    # 'BDT',
    # 'ME1',
    'MB1',
    'jet veto',
    'muon veto',
    # 'N Lep',
    'halo veto',
    # 'DT stn',
    # 'CSCSIZE',
    '1 CSC-DT',
    # 'BLINDSR',
    # 'DR',
    'dPhi',
    # 'dPhi_0.2',
]
CUTS_LOW = [
    'acceptance',
    'HLT',
    'L1',
    'low MET',
    #! I reset cutflow indices here
    'CSC IT',
    'DT IT',
    # 'DT OOT',
    # 'DT I/OOT',
    'MB1',
    'jet veto',
    'halo veto',
    # 'DT stn',
    'DNN',
    '1 CSC-DT',
    'dPhi_0.4',
]
CUTS_HIGH = [
    'acceptance',
    'HLT',
    'L1',
    'high MET',
    #! I reset cutflow indices here
    'CSC IT',
    'DT IT',
    # 'DT OOT',
    # 'DT I/OOT',
    'MB1',
    # 'jet veto',
    'DNN',
    '1 CSC-DT',
    'dPhi_0.4',
]


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
    # 'jetE',
    # 'jetEta',
    # 'jetPhi',
    # 'jetPt',
    # 'jetTightPassId',
    # 'lepDZ',
    # 'lepE',
    # 'lepEta',
    # 'lepPassLooseIso',
    # 'lepPassTightIso',
    # 'lepPassVTightIso',
    # 'lepPassVVTightIso',
    # 'lepPdgId',
    # 'lepPhi',
    # 'lepPt',
    # 'lepTightId',
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


PRINT_CUTFLOW = False
# **************************** #
OUT_DIR = 'reports/weekly/2024-01-22'
T2_OUT_DIR = '/storage/af/user/psimmerl/LLP/mds_analysis'  # os.getcwd()
LOCAL_OUT_DIR = '/home/psimmerl/mds_analysis'  # os.getcwd()

DATA_VERSION = '6'
LUMI = 23.02 * 1000
PI = rt.TMath.Pi()

T2_DATA_DIR = '/storage/cms/store/user/christiw/displacedJetMuonAnalyzer/Run3/V1p19'
LOCAL_DATA_DIR = '/home/psimmerl/mds_analysis/data/raw'  # os.getcwd() + '/data/raw'
DATA_DIR = 'TIER2' if 'caltech' in os.uname()[1] else 'LOCAL'

FN_MC = 'ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted'
FN_R3 = 'DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi'

if 'TIER2' in DATA_DIR:
    OUT_DIR = f'{T2_OUT_DIR}/{OUT_DIR}'
    FN_MC = f'{T2_DATA_DIR}/MC_Summer22EE/v1/sixie/v{DATA_VERSION}/normalized/{FN_MC}.root'
    FN_R3 = f'{T2_DATA_DIR}/Data2022/v{DATA_VERSION}/normalized/{FN_R3}.root'
else:
    OUT_DIR = f'{LOCAL_OUT_DIR}/{OUT_DIR}'
    FN_MC = f'{LOCAL_DATA_DIR}/{FN_MC}_v{DATA_VERSION}.root'
    FN_R3 = f'{LOCAL_DATA_DIR}/{FN_R3}_v{DATA_VERSION}.root'
pathlib.Path(OUT_DIR).mkdir(parents=True, exist_ok=True)

ROOT_ERROR_LEVEL = 1001  # rt.kInfo + 1
rt.gErrorIgnoreLevel = ROOT_ERROR_LEVEL
rt.gROOT.SetBatch(True)

gc = []
# **************************** #
C, D = 'cscRechitCluster', 'dtRechitCluster'


if __name__ == '__main__':
    print('+------------------------------+')
    print('| Starting loo_optimization.py |')
    print('+ -----------------------------+')


    if len(sys.argv) > 1:
        if sys.argv[1] == 'low':
            print('    Using low met cuts')
            CUTS = CUTS_LOW
            MET_CATEGORY = 'low'
        elif sys.argv[1] == 'high':
            print('    Using high met cuts')
            CUTS = CUTS_HIGH
            MET_CATEGORY = 'high'
        else:
            raise ValueError('Unknown arguments: '+' '.join(sys.argv[1:]))

    # # All of the data
    # FN_MC, TN_MC = '/home/psimmerl/mds_analysis/data/raw/ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v6.root', 'MuonSystem'
    # # FN_R3, TN_R3 ='/home/psimmerl/mds_analysis/data/raw/DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi_v6.root', 'MuonSystem'

    # # My skim where it passes HLT_L1CSCCluster_DTCluster50 (HLT 569)
    # # FN_MC, TN_MC ='/home/psimmerl/mds_analysis/data/processed/mc_hlt569.root', 'MuonSystem_HLT569'
    # FN_R3, TN_R3 ='/home/psimmerl/mds_analysis/data/processed/r3_hlt569.root', 'MuonSystem_HLT569'

    # Pedro's DNN skim
    LUMI = 1.1328524540090597e-06 * LUMI # ???
    FN_MC, TN_MC ='/home/psimmerl/mds_analysis/data/processed/mc_pedro.root', 'MuonSystem'
    # FN_R3, TN_R3 ='/home/psimmerl/mds_analysis/data/processed/data_pedro.root', 'MuonSystem'
    FN_R3, TN_R3 ='/home/psimmerl/mds_analysis/data/processed/r3_pedro_hlt569.root', 'MuonSystem'

    if '1 CSC-DT' not in CUTS:
        raise NotImplementedError('cant handle multiple pairs yet')

    # rt.EnableImplicitMT(4)
    # print('    Enabled ROOT\'s implicit multithreading (sometimes causes a crash)')
    
    print('')

    # **************** #

    rng = np.random.default_rng()

    tests = []

    # With CSC time cut optimization
    test_cuts = [
        [000.00, 0000.00, 00.00, -5.00, 12.5, 20.00, 30.00, 50.00, 999.9, 999.9, 1.00, 0.40, 0.40, 0.00, 0.00],
        [121.95, 1104.28, 82.00, -2.27, 6.07, 13.26, 17.24, 48.54, 32.30, 45.82, 3.79, 0.01, 0.58, 0.10, 0.38],
        [112.74, 1144.15, 103.00, -7.87, 4.86, 17.71, 1.00, 62.55, 62.11, 39.70, 7.42, 0.04, 0.35, 0.30, 0.55],
        [112.71, 1032.99, 84.00, -7.27, 2.15, 34.85, 79.81, 46.13, 70.73, 18.38, 0.49, 0.07, 0.24, 0.53, 0.32],
        [107.51, 1122.40, 109.00, -6.94, 4.22, 21.09, 92.61, 52.32, 16.21, 92.34, 9.75, 0.10, 0.06, 0.48, 0.60],
        [107.48, 1080.11, 101.00, -6.82, 3.53, 14.60, 79.05, 75.57, 37.24, 95.08, 9.78, 0.30, 0.30, 0.55, 0.28],
        [106.65, 1173.14, 121.00, -4.39, 3.68, 35.45, 75.97, 1.64, 27.23, 27.19, 3.69, 0.10, 0.38, 0.53, 0.32],
        [106.07, 1011.85, 91.00, -18.50, 8.00, 41.15, 40.06, 36.02, 64.87, 74.79, 5.34, 0.01, 0.52, 0.61, 0.62],
        [104.53, 1034.81, 98.00, -1.79, 2.92, 49.34, 69.95, 14.63, 70.19, 60.26, 0.08, 0.11, 0.18, 0.20, 0.42],
        [100.79, 1144.75, 129.00, -15.99, 1.93, 49.97, 77.46, 88.45, 20.15, 63.64, 6.92, 0.01, 0.60, 0.57, 0.47],
        [100.25, 1235.97, 152.00, -11.25, 5.03, 21.32, 79.95, 33.69, 57.35, 6.05, 2.23, 0.00, 0.58, 0.36, 0.46],
    ]
    # Fixed CSC time cut
    test_cuts = [
        [000.00, 0000.00, 00.00, -5.00, 12.5, 20.00, 30.00, 50.00, 999.9, 999.9, 1.00, 0.40, 0.40, 0.00, 0.00],
        [122.31, 1023.33, 70.00, -5.00, 12.50, 20.00, 178.42, 170.92, 168.68, 199.16, 8.45, 0.05, 0.44, 0.59, 0.69],
        [118.26, 1044.47, 78.00, -5.00, 12.50, 20.00, 57.57, 145.94, 100.43, 77.14, 9.42, 0.09, 0.53, 0.41, 0.74],
        [108.54, 1035.41, 91.00, -5.00, 12.50, 20.00, 86.82, 165.34, 25.03, 135.81, 2.76, 0.12, 0.54, 0.46, 0.60],
        [105.80, 1042.06, 97.00, -5.00, 12.50, 20.00, 66.25, 68.82, 31.59, 88.50, 9.65, 0.06, 0.24, 0.31, 0.78],
        [102.56, 1055.95, 106.00, -5.00, 12.50, 20.00, 171.30, 16.31, 181.93, 24.57, 2.25, 0.07, 0.30, 0.35, 0.68],
        [102.02, 1060.18, 108.00, -5.00, 12.50, 20.00, 158.27, 118.35, 78.16, 5.67, 3.49, 0.05, 0.46, 0.72, 0.36],
        [98.92, 1097.03, 123.00, -5.00, 12.50, 20.00, 106.16, 113.13, 33.15, 176.01, 8.35, 0.11, 0.57, 0.73, 0.30],
        [98.74, 1164.08, 139.00, -5.00, 12.50, 20.00, 65.34, 181.75, 60.82, 150.19, 5.41, 0.06, 0.55, 0.68, 0.29],
        [96.00, 1086.15, 128.00, -5.00, 12.50, 20.00, 81.03, 112.16, 111.71, 58.03, 7.35, 0.01, 0.26, 0.22, 0.79],
        [95.96, 1106.69, 133.00, -5.00, 12.50, 20.00, 71.00, 2.75, 85.21, 46.74, 1.58, 0.13, 0.58, 0.48, 0.32],
    ]

    n_iterations = 100_000
    for i in range(n_iterations):
    # for i in range(len(test_cuts)):
        # SHUFFLE CUTS #
        MIN_CSC_TIME = -5.0 #rng.uniform(-25, 0) # 
        MAX_CSC_TIME =  12.5 #rng.uniform(0, 25) #
        MAX_CSC_TSPREAD =  20.0 #rng.uniform(0, 50) #

        MAX_CSC_JET = rng.uniform(0, 200)
        MAX_DT_JET = rng.uniform(0, 200)
        MAX_CSC_MUON = rng.uniform(0, 200)
        MAX_DT_MUON = rng.uniform(0, 200)

        MAX_MB1 = rng.uniform(0, 10)

        HALO_CUTOFF = rng.uniform(0, 0.5)
        MIN_DPHI = rng.uniform(0, 0.6)

        MIN_CSC_DNN = rng.uniform(0, 1)
        MIN_DT_DNN = rng.uniform(0, 1)
        if i == 0:
            print('old cuts')
            s2b, s, b, MIN_CSC_TIME, MAX_CSC_TIME, MAX_CSC_TSPREAD, MAX_CSC_JET, MAX_DT_JET, MAX_CSC_MUON, MAX_DT_MUON, MAX_MB1, HALO_CUTOFF, MIN_DPHI, MIN_CSC_DNN, MIN_DT_DNN = test_cuts[i]

        # ************ #
        
        rdfs = {
            'mc' : rt.RDataFrame(TN_MC, FN_MC),
            'r3' : rt.RDataFrame(TN_R3, FN_R3),
        }

        for key, rdf in rdfs.items():
            # rdf = rdf.Range(0,100_000) # Skim a subset of events for debugging
            count, weight = rdf.Count().GetValue(), rdf.Sum('weight').GetValue()

            # **** #
            # Create dummy columns to store what cluster indices pass our selections
            rdf = rdf.Define('evtFlag', 'weight > 0')
            rdf = rdf.Define(f'{C}Flag', f'{C}Size > 0')
            rdf = rdf.Define(f'{D}Flag', f'{D}Size > 0')

            rdf = rdf.Define('evtCutFlag', 'weight > 0')
            rdf = rdf.Define(f'{C}CutFlag', f'{C}Size > 0')
            rdf = rdf.Define(f'{D}CutFlag', f'{D}Size > 0')

            for p in (C, D):
                rdf = rdf.Define(f'{p}R', f'ROOT::VecOps::sqrt( {p}X*{p}X + {p}Y*{p}Y )')
            
            if 'mc' in key: # Fix the weight branch
                rdf = rdf.Redefine('weight', f'weight * {LUMI}')

            rdfs[key] = rdf

        # **************** #
        # Apply selections from CUTS (and/or arguments?)
        for key, rdf in rdfs.items():
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
                # Trigger selections (HLT might be wrong for MC so use the reproduced L1 trigger)
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
                if 'low MET' in cut:
                    rdf = rdf.Redefine('evtCutFlag', f'evtCutFlag && (met < 75)')

                if 'high MET' in cut:
                    rdf = rdf.Redefine('evtCutFlag', f'evtCutFlag && ((150 < met) && (met < 200))')

                # **** #
                # Cluster level selections
                if 'CSC IT' in cut:
                    rdf = rdf.Redefine(f'{C}CutFlag', f'{C}CutFlag && ( ({MIN_CSC_TIME} < {C}TimeWeighted) && ({C}TimeWeighted < {MAX_CSC_TIME}) && ({C}TimeSpreadWeightedAll < {MAX_CSC_TSPREAD}) )')
                elif 'CSC OOT' in cut:
                    rdf = rdf.Redefine(f'{C}CutFlag', f'{C}CutFlag && !( ({MIN_CSC_TIME} < {C}TimeWeighted) && ({C}TimeWeighted < {MAX_CSC_TIME}) && ({C}TimeSpreadWeightedAll < {MAX_CSC_TSPREAD}) )')

                if 'DT I/OOT' in cut:
                    if 'mc' in key:
                        cut = 'DT IT*'
                    else:
                        cut = 'DT OOT*'

                if 'DT IT' in cut:
                    rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( ({D}_match_RPCBx_dPhi0p5 == 0) && ({D}_match_RPChits_dPhi0p5 > 0) )')
                elif 'DT OOT' in cut:
                    rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && !( ({D}_match_RPCBx_dPhi0p5 == 0) && ({D}_match_RPChits_dPhi0p5 > 0) )')

                # **** #
                if 'ME1' in cut:
                    rdf = rdf.Redefine(f'{C}CutFlag', f'{C}CutFlag && ( '+
                                    f'({C}NRechitChamberPlus11 == 0) &&'+
                                    f'({C}NRechitChamberPlus12 == 0) &&'+
                                    f'({C}NRechitChamberMinus11 == 0) &&'+
                                    f'({C}NRechitChamberMinus12 == 0) )')

                if 'MB1' in cut:
                    # if 'low MET' in CUTS:
                    #     rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( {D}NHitStation1 == 0 )')
                    # else:
                    rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( {D}NHitStation1 < {MAX_MB1} )')

                if 'DT stn' in cut:
                    rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( ({D}NStation10 < 3) && !(({D}NStation10 == 2) && ({D}MaxStation == 4)) )')

                # **** #
                if 'jet veto' in cut:
                    if 'high MET' in CUTS:
                        rdf = rdf.Redefine(f'{C}CutFlag', f'{C}CutFlag && ( ({C}JetVetoLooseId == 0) || ({C}JetVetoPt < {MAX_CSC_JET}) )')
                        rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( ({D}JetVetoLooseId == 0) || ({D}JetVetoPt < {MAX_DT_JET}) )')
                    else:
                        rdf = rdf.Redefine(f'{C}CutFlag', f'{C}CutFlag && ( ({C}JetVetoLooseId == 0) || ({C}JetVetoPt < {MAX_CSC_JET}) )')
                        rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( ({D}JetVetoLooseId == 0) || ({D}JetVetoPt < {MAX_DT_JET}) )')
                
                if 'muon veto' in cut:
                    rdf = rdf.Redefine(f'{C}CutFlag', f'{C}CutFlag && ( ({C}MuonVetoGlobal == 0) || ({C}MuonVetoPt < {MAX_CSC_MUON}) )')
                    rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( ({D}MuonVetoLooseId == 0) || ({D}MuonVetoPt < {MAX_DT_MUON}) )')

                if 'halo veto' in cut:
                    rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( ({HALO_CUTOFF} < abs({D}Phi)) && (abs({D}Phi) < {PI} - {HALO_CUTOFF}) )')

                # **** #
                # if 'BDT' in cut:
                #     raise NotImplementedError('BDT')

                if 'DNN' in cut:
                    # rdf = rdf.Redefine(f'{C}CutFlag', f'{C}CutFlag && ( Take({C}DNN,nCscRechitClusters) > 0.85 )')
                    # rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( Take({D}DNN,nDtRechitClusters) > 0.85 )')
                    if 'high MET' in CUTS:
                        rdf = rdf.Redefine(f'{C}CutFlag', f'{C}CutFlag && ( Take({C}DNN,nCscRechitClusters) > {MIN_CSC_DNN} )')
                        rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( Take({D}DNN,nDtRechitClusters) > {MIN_DT_DNN} )')
                    else:
                        rdf = rdf.Redefine(f'{C}CutFlag', f'{C}CutFlag && ( Take({C}DNN,nCscRechitClusters) > {MIN_CSC_DNN} )')
                        rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( Take({D}DNN,nDtRechitClusters) > {MIN_DT_DNN} )')

                # **** #
                if '1 CSC-DT' in cut:
                    # rdf = rdf.Redefine('evtCutFlag', f'evtCutFlag && evtFlag && '+
                    #                    f'(reduce({C}Flag.begin(), {C}Flag.end()) == 1) && '+
                    #                    f'(reduce({D}Flag.begin(), {D}Flag.end()) == 1)')
                    rdf = rdf.Redefine(f'{C}Flag', f'{C}Flag && ( {C}Size == Max({C}Size*{C}Flag) )')
                    rdf = rdf.Redefine(f'{D}Flag', f'{D}Flag && ( {D}Size == Max({D}Size*{D}Flag) )')

                    for col in COLUMNS_OUT:
                        ocol = col.replace('RechitCluster', '')
                        if C in col[:len(C)]:
                            rdf = rdf.Define(ocol, f'{col}[{C}Flag][0]')
                        elif D in col[:len(D)]:
                            rdf = rdf.Define(ocol, f'{col}[{D}Flag][0]')

                    # Apply our cluster level selections to relevant columns
                    rdf = rdf.Define('cscCTau', f'gLLP_ctau[{C}_match_gLLP_index[{C}Flag][0]]')
                    rdf = rdf.Define('dtCTau', f'gLLP_ctau[{D}_match_gLLP_index[{D}Flag][0]]')
                    rdf = rdf.Define('tag_dEta', 'abs(cscEta - dtEta)')
                    rdf = rdf.Define('tag_dPhi', f'abs( abs(cscPhi - dtPhi) - (abs(cscPhi - dtPhi) > {PI} ? 2*{PI} : 0) )')
                    rdf = rdf.Define('tag_dR', 'sqrt(tag_dEta*tag_dEta + tag_dPhi*tag_dPhi)')

                # **** #
                # if 'BLINDSR' in cut:
                #     raise NotImplementedError('BLINDSR')

                # **** #
                # if 'DR' in cut:
                #     raise NotImplementedError('DR')
                # if 'DETA' in cut: #! Has not been tested -- if gg->H->ss is exclusive then conservation
                #     rdf = rdf.Redefine('evtCutFlag', '((0 < cscEta) && ( dtEta < 0)) || (( cscEta < 0) && (0 < dtEta)) ||')
                if 'dPhi' in cut:
                    rdf = rdf.Redefine('evtCutFlag', f'evtCutFlag && (tag_dPhi > {MIN_DPHI})')

                # **** #
                # Finally require that we still have just one or at least one pair of CSC+DT per event
                rdf = rdf.Redefine('evtCutFlag', 'evtCutFlag && ('+
                                f'(reduce({C}CutFlag.begin(), {C}CutFlag.end()) > 0) && '+
                                f'(reduce({D}CutFlag.begin(), {D}CutFlag.end()) > 0))')

                # **** #
                # Propagate to cumulative flags
                rdf = rdf.Redefine('evtFlag', 'evtFlag && evtCutFlag')
                rdf = rdf.Redefine(f'{C}Flag', f'{C}Flag && {C}CutFlag')
                rdf = rdf.Redefine(f'{D}Flag', f'{D}Flag && {D}CutFlag')

                rdf = rdf.Redefine('evtFlag', 'evtFlag && ('+
                                f'(reduce({C}Flag.begin(), {C}Flag.end()) > 0) && '+
                                f'(reduce({D}Flag.begin(), {D}Flag.end()) > 0))')

                # rdf = rdf.Redefine(f'nC{C[1:]}s', f'reduce({C}Flag.begin(), {C}Flag.end())')
                # rdf = rdf.Redefine(f'nD{D[1:]}s', f'reduce({D}Flag.begin(), {D}Flag.end())')

                # **** #
                rdf = rdf.Filter('evtFlag')
                rdf = rdf.Redefine(f'{C}Size', f'{C}Size * {C}Flag')
                rdf = rdf.Redefine(f'{D}Size', f'{D}Size * {D}Flag')

            rdfs[key] = rdf.Filter('evtFlag') # Apply event level cut and update the dictionary with the new RDF
            
        s = rdfs['mc'].Sum('weight').GetValue()
        b = rdfs['r3'].Sum('weight').GetValue()
        s2b = s / ( b )**(1/2) if b else 0

        # **************** #
        # if s>= 4 and b >= 3:
        #     sizes = np.linspace(50,150,50) #if ABCD_SIZE is None else np.array([ABCD_SIZE])
        #     dphis = np.linspace(0.4,np.pi,50) #if ABCD_DPHI is None else np.array([ABCD_DPHI])

        #     # if 'low MET' in CUTS:
        #     a_mc = [ [rdfs['mc'].Filter(f'(tag_dPhi >= {dp}) && (dtSize >= {sz})').Sum('weight') for sz in sizes] for dp in dphis ]
        #     b_mc = [ [rdfs['mc'].Filter(f'(tag_dPhi >= {dp}) && (dtSize < {sz})').Sum('weight') for sz in sizes] for dp in dphis ]
        #     c_mc = [ [rdfs['mc'].Filter(f'(tag_dPhi < {dp}) && (dtSize < {sz})').Sum('weight') for sz in sizes] for dp in dphis ]
        #     d_mc = [ [rdfs['mc'].Filter(f'(tag_dPhi < {dp}) && (dtSize >= {sz})').Sum('weight') for sz in sizes] for dp in dphis ]

        #     a_r3 = [ [rdfs['r3'].Filter(f'(tag_dPhi >= {dp}) && (dtSize >= {sz})').Sum('weight') for sz in sizes] for dp in dphis ]
        #     b_r3 = [ [rdfs['r3'].Filter(f'(tag_dPhi >= {dp}) && (dtSize < {sz})').Sum('weight') for sz in sizes] for dp in dphis ]
        #     c_r3 = [ [rdfs['r3'].Filter(f'(tag_dPhi < {dp}) && (dtSize < {sz})').Sum('weight') for sz in sizes] for dp in dphis ]
        #     d_r3 = [ [rdfs['r3'].Filter(f'(tag_dPhi < {dp}) && (dtSize >= {sz})').Sum('weight') for sz in sizes] for dp in dphis ]
        #     # if 'high MET' in CUTS:
        #     #     a_mc = [ [rdfs['mc'].Filter(f'(tag_dPhi < {dp}) && (dtSize >= {sz})').Sum('weight') for sz in sizes] for dp in dphis ]
        #     #     b_mc = [ [rdfs['mc'].Filter(f'(tag_dPhi < {dp}) && (dtSize < {sz})').Sum('weight') for sz in sizes] for dp in dphis ]
        #     #     c_mc = [ [rdfs['mc'].Filter(f'(tag_dPhi >= {dp}) && (dtSize < {sz})').Sum('weight') for sz in sizes] for dp in dphis ]
        #     #     d_mc = [ [rdfs['mc'].Filter(f'(tag_dPhi >= {dp}) && (dtSize >= {sz})').Sum('weight') for sz in sizes] for dp in dphis ]

        #     #     a_r3 = [ [rdfs['r3'].Filter(f'(tag_dPhi < {dp}) && (dtSize >= {sz})').Sum('weight') for sz in sizes] for dp in dphis ]
        #     #     b_r3 = [ [rdfs['r3'].Filter(f'(tag_dPhi < {dp}) && (dtSize < {sz})').Sum('weight') for sz in sizes] for dp in dphis ]
        #     #     c_r3 = [ [rdfs['r3'].Filter(f'(tag_dPhi >= {dp}) && (dtSize < {sz})').Sum('weight') for sz in sizes] for dp in dphis ]
        #     #     d_r3 = [ [rdfs['r3'].Filter(f'(tag_dPhi >= {dp}) && (dtSize >= {sz})').Sum('weight') for sz in sizes] for dp in dphis ]
        #     # else:
        #     #     raise

        #     a_mc = np.array([ [xx.GetValue() for xx in x] for x in a_mc])
        #     b_mc = np.array([ [xx.GetValue() for xx in x] for x in b_mc])
        #     c_mc = np.array([ [xx.GetValue() for xx in x] for x in c_mc])
        #     d_mc = np.array([ [xx.GetValue() for xx in x] for x in d_mc])

        #     a_r3 = np.array([ [xx.GetValue() for xx in x] for x in a_r3])
        #     b_r3 = np.array([ [xx.GetValue() for xx in x] for x in b_r3])
        #     c_r3 = np.array([ [xx.GetValue() for xx in x] for x in c_r3])
        #     d_r3 = np.array([ [xx.GetValue() for xx in x] for x in d_r3])

        #     ap_mc = np.divide(b_mc*d_mc, c_mc, where=c_mc>0, out=np.zeros_like(a_mc))
        #     ap_r3 = np.divide(b_r3*d_r3, c_r3, where=c_r3>0, out=np.zeros_like(a_r3))

        #     ####

        #     # There has to be signal and data in each of the 4 bins for the ABCD method to work
        #     cond = ((a_mc>0) & (b_mc>0) & (c_mc>0) & (d_mc>0)) & \
        #             ((ap_r3>0) & (b_r3>0) & (c_r3>0) & (d_r3>0))

        #     # Here we add a condition that the predicted value must be within 2 sigma of the actual value
        #     ae_r3 = np.sqrt(a_r3)
        #     ape_r3 = ap_r3 * np.sqrt(
        #         np.divide(1,b_r3, where=cond, out=np.zeros_like(b_r3)) +
        #         np.divide(1,c_r3, where=cond, out=np.zeros_like(c_r3)) +
        #         np.divide(1,d_r3, where=cond, out=np.zeros_like(d_r3))
        #     )

        #     cond = cond & (np.divide(np.abs(a_r3-ap_r3),np.sqrt(ae_r3**2 + ape_r3**2), where=cond, out=999*np.zeros_like(a_r3))<1)
                        
        #     ####
        #     s2bs = np.divide(a_mc, np.sqrt(ap_r3), where=cond, out=np.zeros_like(a_mc))

        #     idx_s2b = np.unravel_index(np.argmax(s2bs), s2bs.shape)

        #     s2b_a = s2bs[idx_s2b]
        #     s_a = a_mc[idx_s2b]
        #     b_a = ap_r3[idx_s2b]

        #     ABCD_DPHI = dphis[idx_s2b[0]]
        #     ABCD_SIZE = sizes[idx_s2b[1]]
        # else:
        #     s2b_a = 0
        #     s_a = 0
        #     b_a = 0
        #     ABCD_DPHI = 0
        #     ABCD_SIZE = 0
            
        # **************** #

        tests.append([
            s2b,
            s,
            b,
            # s2b_a,
            # s_a,
            # b_a,
            # ABCD_DPHI,
            # ABCD_SIZE,
            MIN_CSC_TIME,
            MAX_CSC_TIME,
            MAX_CSC_TSPREAD,
            MAX_CSC_JET,
            MAX_DT_JET,
            MAX_CSC_MUON,
            MAX_DT_MUON,
            MAX_MB1,
            HALO_CUTOFF,
            MIN_DPHI,
            MIN_CSC_DNN,
            MIN_DT_DNN,
        ])

        with open(f'loo_tests_{MET_CATEGORY}.pkl', 'wb') as fout:
            pickle.dump(tests, fout)
        
        print(f'{s2b:>7.2f}, {s:>7.2f}, {b:>7.2f} | {MIN_CSC_TIME:>6.2f}, {MAX_CSC_TIME:>5.2f}, {MAX_CSC_TSPREAD:>5.2f} | {MAX_CSC_JET:>6.2f}, {MAX_DT_JET:>6.2f} | {MAX_CSC_MUON:>6.2f}, {MAX_DT_MUON:>6.2f} | {MAX_MB1:>5.2f}, {HALO_CUTOFF:>4.2f}, {MIN_DPHI:>4.2f} | {MIN_CSC_DNN:>4.2f}, {MIN_DT_DNN:>4.2f}')
        # print(f'{s2b:>7.2f} {s:>7.2f}, {b:>7.2f} | {s2b_a:>7.2f} {s_a:>7.2f}, {b_a:>6.2f} | {ABCD_DPHI:>4.2f}, {ABCD_SIZE:>6.2f} | {MIN_CSC_TIME:>6.2f}, {MAX_CSC_TIME:>5.2f}, {MAX_CSC_TSPREAD:>5.2f} | {MAX_CSC_JET:>6.2f}, {MAX_DT_JET:>6.2f} | {MAX_CSC_MUON:>6.2f}, {MAX_DT_MUON:>6.2f} | {HALO_CUTOFF:>4.2f} | {MIN_CSC_DNN:>4.2f}, {MIN_DT_DNN:>4.2f}')
        

        if i % 50 == 0 or i+1 == n_iterations:
            _t = np.array(tests)
            # idx_best = np.argmax(_t[:,3])
            print('')
            print(f'Completed {i+1} iterations')
            print('Best S/sqrt[ B ]:')
            idxs = np.argsort(_t[:,0])[::-1]
            for iidx, idx in enumerate(idxs[:10]):
                print(f'{iidx}: ' + ', '.join([f'{x:.2f}' for x in tests[idx]]))
            
            print('Best S/sqrt[ B ], with S > 1000:')
            idxs = np.argsort(_t[:,0] * (_t[:,1]>1000))[::-1]
            for iidx, idx in enumerate(idxs[:10]):
                print(f'{iidx}: ' + ', '.join([f'{x:.2f}' for x in tests[idx]]))
            print('')

