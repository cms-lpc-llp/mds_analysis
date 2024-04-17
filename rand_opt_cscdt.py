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
CUTS_L1 = [
    'acceptance',
    'HLT',
    'L1',
    'MET',
    #! I reset cutflow indices here
    'CSC IT',
    'DT IT', # can set OOT using command line
    'MB1',
    '1 CSC-DT',
    # 'dPhi',
]

CUTS = [
    'acceptance',
    'HLT',
    'L1',
    'MET',
    #! I reset cutflow indices here
    'CSC IT',
    'DT IT', # can set OOT using command line
    # 'DT OOT',
    # 'DT I/OOT',
    # 'MET',
    # 'DNN',
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
    print('+-----------------------+')
    print('| Starting loo_cscdt.py |')
    print('+ ----------------------+')

    n_iterations = 100_000
    SAVE_BEST = False

    if len(sys.argv) > 1:
        args = ' '.join(sys.argv[1:])
        if 'l1' in args:
            print('    Using the reduced cut set')
            CUTS = CUTS_L1

        if 'low' in args:
            print('    Using low met cuts')
            CUTS = [c.replace('MET', 'low MET') if 'MET' == c else c for c in CUTS]
            # CUTS = CUTS_LOW
            MET_CATEGORY = 'low'
        if 'high' in args:
            print('    Using high met cuts')
            CUTS = [c.replace('MET', 'high MET') if 'MET' == c else c for c in CUTS]
            # CUTS = CUTS_HIGH
            MET_CATEGORY = 'high'

        if 'oot' in args:
            print('    Using out-of-time 2nd cluster, DATA ONLY')
            CUTS = [c.replace('DT IT', 'DT OOT') if 'DT IT' == c else c for c in CUTS]
        else:#if 'it' in args:
            print('    Using in-time 2nd cluster')

        if 'save' in args:
            SAVE_BEST = True

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

    rt.EnableImplicitMT(4)
    print('    Enabled ROOT\'s implicit multithreading (sometimes causes a crash)')
    
    print('')

    # **************** #

    rng = np.random.default_rng()
    tests = []

# LOW MET, Min dPhi < 1.5, 1350 iterations
# 0: 49.94, 1121.19, 504.00, -5.00, 12.50, 20.00, 15.00, 10.00, 45.00, 95.00, 0.00, 0.00, 1.35, 0.00, 0.00
# 1: 49.55, 1155.62, 544.00, -5.00, 12.50, 20.00, 55.00, 10.00, 65.00, 80.00, 0.00, 0.00, 1.50, 0.00, 0.00
# 2: 47.66, 1136.90, 569.00, -5.00, 12.50, 20.00, 35.00, 10.00, 30.00, 60.00, 0.00, 0.00, 1.40, 0.00, 0.00
# 3: 45.63, 1154.42, 640.00, -5.00, 12.50, 20.00, 80.00, 10.00, 45.00, 75.00, 0.00, 0.00, 1.20, 0.00, 0.00
# 4: 45.35, 1083.74, 571.00, -5.00, 12.50, 20.00, 10.00, 100.00, 45.00, 100.00, 0.00, 0.00, 1.25, 0.00, 0.00
# 5: 44.38, 1135.69, 655.00, -5.00, 12.50, 20.00, 20.00, 10.00, 90.00, 40.00, 0.00, 0.00, 1.15, 0.00, 0.00
# 6: 44.26, 1150.79, 676.00, -5.00, 12.50, 20.00, 10.00, 15.00, 70.00, 85.00, 1.00, 0.00, 1.30, 0.00, 0.00
# 7: 43.85, 1172.54, 715.00, -5.00, 12.50, 20.00, 75.00, 10.00, 15.00, 40.00, 1.00, 0.00, 1.50, 0.00, 0.00
# 8: 43.66, 1239.59, 806.00, -5.00, 12.50, 20.00, 85.00, 10.00, 55.00, 70.00, 1.00, 0.00, 1.35, 0.00, 0.00
# 9: 43.14, 1124.21, 679.00, -5.00, 12.50, 20.00, 10.00, 75.00, 15.00, 95.00, 1.00, 0.00, 1.45, 0.00, 0.00
# Min dPhi = 0.4
# Completed 4151 iterations
# Best S/sqrt[ B ], B >= 3:
# 0: 36.77, 1077.09, 858.00, -5.00, 12.50, 20.00, 0.00, 1.00, 85.00, 10.00, 15.00, 40.00, -1.00, 0.00, 0.05, 0.40, 0.00, 0.00
# 1: 36.51, 1036.62, 806.00, -5.00, 12.50, 20.00, 0.00, 1.00, 100.00, 10.00, 90.00, 60.00, -1.00, 0.00, 0.15, 0.40, 0.00, 0.00
# 2: 36.38, 1015.48, 779.00, -5.00, 12.50, 20.00, 0.00, 1.00, 10.00, 70.00, 55.00, 85.00, -1.00, 0.00, 0.10, 0.40, 0.00, 0.00
# 3: 36.29, 1013.66, 780.00, -5.00, 12.50, 20.00, 0.00, 1.00, 100.00, 10.00, 95.00, 95.00, -1.00, 0.00, 0.20, 0.40, 0.00, 0.00
# 4: 36.23, 1010.64, 778.00, -5.00, 12.50, 20.00, 0.00, 1.00, 25.00, 10.00, 85.00, 90.00, -1.00, 0.00, 0.20, 0.40, 0.00, 0.00
# 5: 36.13, 1005.21, 774.00, -5.00, 12.50, 20.00, 0.00, 1.00, 30.00, 10.00, 65.00, 80.00, -1.00, 0.00, 0.20, 0.40, 0.00, 0.00
# 6: 36.09, 996.14, 762.00, -5.00, 12.50, 20.00, 0.00, 1.00, 25.00, 10.00, 30.00, 70.00, -1.00, 0.00, 0.20, 0.40, 0.00, 0.00
# 7: 35.79, 1016.08, 806.00, -5.00, 12.50, 20.00, 0.00, 1.00, 90.00, 10.00, 10.00, 45.00, -1.00, 0.00, 0.10, 0.40, 0.00, 0.00
# 8: 35.67, 1178.58, 1092.00, -5.00, 12.50, 20.00, 0.00, 1.00, 15.00, 10.00, 15.00, 85.00, -1.00, 1.00, 0.00, 0.40, 0.00, 0.00
# 9: 35.53, 1045.68, 866.00, -5.00, 12.50, 20.00, 0.00, 1.00, 10.00, 85.00, 85.00, 25.00, -1.00, 0.00, 0.00, 0.40, 0.00, 0.00
# Best S/sqrt[ B ], S > 900, B >= 3:
# 0: 36.77, 1077.09, 858.00, -5.00, 12.50, 20.00, 0.00, 1.00, 85.00, 10.00, 15.00, 40.00, -1.00, 0.00, 0.05, 0.40, 0.00, 0.00
# 1: 36.51, 1036.62, 806.00, -5.00, 12.50, 20.00, 0.00, 1.00, 100.00, 10.00, 90.00, 60.00, -1.00, 0.00, 0.15, 0.40, 0.00, 0.00
# 2: 36.38, 1015.48, 779.00, -5.00, 12.50, 20.00, 0.00, 1.00, 10.00, 70.00, 55.00, 85.00, -1.00, 0.00, 0.10, 0.40, 0.00, 0.00
# 3: 36.29, 1013.66, 780.00, -5.00, 12.50, 20.00, 0.00, 1.00, 100.00, 10.00, 95.00, 95.00, -1.00, 0.00, 0.20, 0.40, 0.00, 0.00
# 4: 36.23, 1010.64, 778.00, -5.00, 12.50, 20.00, 0.00, 1.00, 25.00, 10.00, 85.00, 90.00, -1.00, 0.00, 0.20, 0.40, 0.00, 0.00
# 5: 36.13, 1005.21, 774.00, -5.00, 12.50, 20.00, 0.00, 1.00, 30.00, 10.00, 65.00, 80.00, -1.00, 0.00, 0.20, 0.40, 0.00, 0.00
# 6: 36.09, 996.14, 762.00, -5.00, 12.50, 20.00, 0.00, 1.00, 25.00, 10.00, 30.00, 70.00, -1.00, 0.00, 0.20, 0.40, 0.00, 0.00
# 7: 35.79, 1016.08, 806.00, -5.00, 12.50, 20.00, 0.00, 1.00, 90.00, 10.00, 10.00, 45.00, -1.00, 0.00, 0.10, 0.40, 0.00, 0.00
# 8: 35.67, 1178.58, 1092.00, -5.00, 12.50, 20.00, 0.00, 1.00, 15.00, 10.00, 15.00, 85.00, -1.00, 1.00, 0.00, 0.40, 0.00, 0.00
# 9: 35.53, 1045.68, 866.00, -5.00, 12.50, 20.00, 0.00, 1.00, 10.00, 85.00, 85.00, 25.00, -1.00, 0.00, 0.00, 0.40, 0.00, 0.00

# HIGH MET, Min dPhi < 1.5, 1400 iterations
# 0: 23.02, 39.87, 3.00, -5.00, 12.50, 20.00, 45.00, 10.00, 40.00, 100.00, 7.00, 0.00, 0.05, 0.00, 0.00
# 1: 22.67, 39.27, 3.00, -5.00, 12.50, 20.00, 75.00, 10.00, 15.00, 95.00, 7.00, 0.00, 0.05, 0.00, 0.00
# 2: 22.32, 38.66, 3.00, -5.00, 12.50, 20.00, 100.00, 10.00, 40.00, 85.00, 10.00, 0.00, 0.15, 0.00, 0.00
# 3: 21.62, 37.45, 3.00, -5.00, 12.50, 20.00, 60.00, 10.00, 20.00, 70.00, 9.00, 0.00, 0.10, 0.00, 0.00
# 4: 20.58, 35.64, 3.00, -5.00, 12.50, 20.00, 55.00, 90.00, 35.00, 95.00, 1.00, 0.00, 0.00, 0.00, 0.00
# 5: 20.23, 35.04, 3.00, -5.00, 12.50, 20.00, 95.00, 15.00, 20.00, 65.00, 4.00, 0.00, 0.15, 0.00, 0.00
# 6: 19.88, 34.43, 3.00, -5.00, 12.50, 20.00, 15.00, 10.00, 25.00, 60.00, 7.00, 0.00, 0.25, 0.00, 0.00
# 7: 19.63, 39.27, 4.00, -5.00, 12.50, 20.00, 80.00, 10.00, 25.00, 30.00, 9.00, 0.00, 0.00, 0.00, 0.00
# 8: 19.63, 39.27, 4.00, -5.00, 12.50, 20.00, 65.00, 10.00, 45.00, 95.00, 9.00, 0.00, 0.10, 0.00, 0.00
# 9: 19.63, 39.27, 4.00, -5.00, 12.50, 20.00, 35.00, 10.00, 70.00, 80.00, 7.00, 0.00, 0.05, 0.00, 0.00
# Min dPhi = 0.4
# Completed 2501 iterations
# Best S/sqrt[ B ], B >= 3:
# 0: 20.58, 35.64, 3.00, -5.00, 12.50, 20.00, 0.00, 1.00, 35.00, 10.00, 70.00, 95.00, -1.00, 9.00, 0.00, 0.40, 0.00, 0.00
# 1: 19.53, 33.83, 3.00, -5.00, 12.50, 20.00, 0.00, 1.00, 75.00, 15.00, 10.00, 85.00, -1.00, 8.00, 0.00, 0.40, 0.00, 0.00
# 2: 19.53, 33.83, 3.00, -5.00, 12.50, 20.00, 0.00, 1.00, 95.00, 15.00, 10.00, 75.00, -1.00, 9.00, 0.00, 0.40, 0.00, 0.00
# 3: 19.18, 33.22, 3.00, -5.00, 12.50, 20.00, 0.00, 1.00, 45.00, 15.00, 25.00, 50.00, -1.00, 10.00, 0.05, 0.40, 0.00, 0.00
# 4: 19.18, 33.22, 3.00, -5.00, 12.50, 20.00, 0.00, 1.00, 55.00, 20.00, 40.00, 75.00, -1.00, 6.00, 0.00, 0.40, 0.00, 0.00
# 5: 19.18, 33.22, 3.00, -5.00, 12.50, 20.00, 0.00, 1.00, 30.00, 65.00, 35.00, 95.00, -1.00, 4.00, 0.00, 0.40, 0.00, 0.00
# 6: 19.18, 33.22, 3.00, -5.00, 12.50, 20.00, 0.00, 1.00, 85.00, 10.00, 45.00, 50.00, -1.00, 8.00, 0.00, 0.40, 0.00, 0.00
# 7: 19.18, 33.22, 3.00, -5.00, 12.50, 20.00, 0.00, 1.00, 55.00, 75.00, 15.00, 95.00, -1.00, 5.00, 0.00, 0.40, 0.00, 0.00
# 8: 19.18, 33.22, 3.00, -5.00, 12.50, 20.00, 0.00, 1.00, 95.00, 15.00, 15.00, 60.00, -1.00, 8.00, 0.00, 0.40, 0.00, 0.00
# 9: 18.83, 32.62, 3.00, -5.00, 12.50, 20.00, 0.00, 1.00, 25.00, 75.00, 40.00, 80.00, -1.00, 6.00, 0.05, 0.40, 0.00, 0.00

    for i in range(n_iterations):
        MIN_CSC_TIME = -5.0
        MAX_CSC_TIME =  12.5
        MAX_CSC_TSPREAD =  20.0

        MAX_RPC_BX = 0
        MIN_RPC_HITS = 1

        MAX_CSC_JET = 30
        MAX_DT_JET = 50
        MAX_CSC_MUON = 30
        MAX_DT_MUON = 10

        MAX_ME1 = -1
        MAX_MB1 = 0

        HALO_CUTOFF = 0.4
        MIN_DPHI = 0.4

        MIN_CSC_DNN = 0
        MIN_DT_DNN = 0
        if i == 0:
            print('Running first iteration with default cuts')
        else:
            MIN_CSC_TIME = -5.0 # rng.uniform(-25, 0)
            MAX_CSC_TIME =  12.5 # rng.uniform(0, 25)
            MAX_CSC_TSPREAD =  20.0 # rng.uniform(0, 50)

            MAX_RPC_BX = 0
            MIN_RPC_HITS = 1

            MAX_CSC_JET = rng.choice(np.arange(10, 100+5, 5)) # rng.uniform(10, 100) # 30
            MAX_DT_JET = rng.choice(np.arange(10, 100+5, 5)) # rng.uniform(10, 100) # 50
            MAX_CSC_MUON = rng.choice(np.arange(10, 100+5, 5)) # rng.uniform(10, 100) # 30
            MAX_DT_MUON = rng.choice(np.arange(10, 100+5, 5)) # rng.uniform(10, 100) # 10

            MAX_ME1 = -1 # rng.choice(np.arange(0, 10+1, 1)) # rng.uniform(0, 50) # 0
            MAX_MB1 = rng.choice(np.arange(0, 10+1, 1)) # rng.uniform(0, 50) # 0

            HALO_CUTOFF = rng.choice(np.arange(0, 0.6+0.05, 0.05))
            MIN_DPHI = 0.4 # rng.choice(np.arange(0, 0.6+0.05, 0.05)) # rng.uniform(0, 1.5) # 1

            MIN_CSC_DNN = rng.choice(np.arange(0, 1+0.05, 0.05)) # rng.uniform(0, 1) # 0
            MIN_DT_DNN = rng.choice(np.arange(0, 1+0.05, 0.05)) # rng.uniform(0, 1) # 0

        if 'DNN' not in CUTS:
            MIN_CSC_DNN, MIN_DT_DNN = 0, 0

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
                if f'{p}R' not in COLUMNS_OUT:
                    COLUMNS_OUT.append(f'{p}R')
            
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
                if 'MET' in cut:
                    rdf = rdf.Redefine('evtCutFlag', f'evtCutFlag && (met < 200)')

                if 'low MET' in cut:
                    rdf = rdf.Redefine('evtCutFlag', f'evtCutFlag && (met < 75)')

                if 'high MET' in cut:
                    rdf = rdf.Redefine('evtCutFlag', f'evtCutFlag && ((150 < met) && (met < 200))')

                # **** #
                # Cluster level selections
                if 'mc' in key and ' OOT' in cut:
                    cut = cut.replace(' OOT', ' IT')

                if 'CSC IT' in cut or 'mc' in key:
                    rdf = rdf.Redefine(f'{C}CutFlag', f'{C}CutFlag && ( ({MIN_CSC_TIME} < {C}TimeWeighted) && ({C}TimeWeighted < {MAX_CSC_TIME}) && ({C}TimeSpreadWeightedAll < {MAX_CSC_TSPREAD}) )')
                elif 'CSC OOT' in cut:
                    rdf = rdf.Redefine(f'{C}CutFlag', f'{C}CutFlag && !( ({MIN_CSC_TIME} < {C}TimeWeighted) && ({C}TimeWeighted < {MAX_CSC_TIME}) && ({C}TimeSpreadWeightedAll < {MAX_CSC_TSPREAD}) )')

                if 'DT IT' in cut or 'mc' in key:
                    rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( (abs({D}_match_RPCBx_dPhi0p5) <= {MAX_RPC_BX}) && ({D}_match_RPChits_dPhi0p5 >= {MIN_RPC_HITS}) )')
                elif 'DT OOT' in cut:
                    rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && !( (abs({D}_match_RPCBx_dPhi0p5) <= {MAX_RPC_BX}) && ({D}_match_RPChits_dPhi0p5 >= {MIN_RPC_HITS}) )')

                # **** #
                if 'ME1' in cut:
                    rdf = rdf.Redefine(f'{C}1CutFlag', f'{C}1CutFlag && ( '+
                                    f'({C}NRechitChamberPlus11 <= {MAX_ME1}) &&'+
                                    f'({C}NRechitChamberPlus12 <= {MAX_ME1}) &&'+
                                    f'({C}NRechitChamberMinus11 <= {MAX_ME1}) &&'+
                                    f'({C}NRechitChamberMinus12 <= {MAX_ME1}) )')
                    # rdf = rdf.Redefine(f'{C}1CutFlag', f'{C}1CutFlag && ( '+
                    #                 f'({C}NRechitChamberPlus11 +'+
                    #                 f'{C}NRechitChamberPlus12 +'+
                    #                 f'{C}NRechitChamberMinus11 +'+
                    #                 f'{C}NRechitChamberMinus12 <= {MAX_ME1}) ')

                if 'MB1' in cut:
                    # rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( {D}NHitStation1 == 0 )')
                    rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( {D}NHitStation1 <= {MAX_MB1} )')

                if 'DT stn' in cut:
                    rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( ({D}NStation10 < 3) && !(({D}NStation10 == 2) && ({D}MaxStation == 4)) )')

                # **** #
                if 'jet veto' in cut:
                    rdf = rdf.Redefine(f'{C}CutFlag', f'{C}CutFlag && ( ({C}JetVetoLooseId == 0) && ({C}JetVetoPt < {MAX_CSC_JET}) )')
                    rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( ({D}JetVetoLooseId == 0) && ({D}JetVetoPt < {MAX_DT_JET}) )')
                
                if 'muon veto' in cut:
                    rdf = rdf.Redefine(f'{C}CutFlag', f'{C}CutFlag && ( ({C}MuonVetoGlobal == 0) && ({C}MuonVetoPt < {MAX_CSC_MUON}) )')
                    rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( ({D}MuonVetoLooseId == 0) && ({D}MuonVetoPt < {MAX_DT_MUON}) )')

                if 'halo veto' in cut:
                    rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( ({HALO_CUTOFF} < abs({D}Phi)) && (abs({D}Phi) < {PI} - {HALO_CUTOFF}) )')

                # **** #
                # if 'BDT' in cut:
                #     raise NotImplementedError('BDT')

                if 'DNN' in cut:
                    rdf = rdf.Redefine(f'{C}CutFlag', f'{C}CutFlag && ( Take({C}DNN,nCscRechitClusters) > {MIN_CSC_DNN} )')
                    rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( Take({D}DNN,nDtRechitClusters) > {MIN_DT_DNN} )')

                # **** #
                if '1 CSC-DT' in cut:
                    # rdf = rdf.Redefine('evtCutFlag', f'evtCutFlag && evtFlag && '+
                    #                    f'(reduce({C}Flag.begin(), {C}Flag.end()) == 1) && '+
                    #                    f'(reduce({D}Flag.begin(), {D}Flag.end()) == 1)')

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

                    # Apply our cluster level selections to relevant columns
                    rdf = rdf.Define('cscCTau', f'gLLP_ctau[{C}_match_gLLP_index[{C}Flag][0]]')
                    rdf = rdf.Define('dtCTau', f'gLLP_ctau[{D}_match_gLLP_index[{D}Flag][0]]')
                    rdf = rdf.Define('tag_dEta', 'abs(cscEta - dtEta)')
                    rdf = rdf.Define('tag_dPhi', f'abs( abs(cscPhi - dtPhi) - (abs(cscPhi - dtPhi) > {PI} ? 2*{PI} : 0) )')
                    rdf = rdf.Define('tag_dR', 'sqrt(tag_dEta*tag_dEta + tag_dPhi*tag_dPhi)')
                    columns_out.extend(['tag_dEta', 'tag_dPhi', 'tag_dR', 'cscCTau', 'dtCTau'])

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
                # Propagate to cumulative flags
                rdf = rdf.Redefine(f'{C}Flag', f'{C}Flag && {C}CutFlag')
                rdf = rdf.Redefine(f'{D}Flag', f'{D}Flag && {D}CutFlag')

                rdf = rdf.Redefine('evtFlag', 'evtFlag && evtCutFlag && ('+
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

        tests.append([
            s2b,
            s,
            b,
            MIN_CSC_TIME,
            MAX_CSC_TIME,
            MAX_CSC_TSPREAD,
            MAX_RPC_BX,
            MIN_RPC_HITS,
            MAX_CSC_JET,
            MAX_DT_JET,
            MAX_CSC_MUON,
            MAX_DT_MUON,
            MAX_ME1,
            MAX_MB1,
            HALO_CUTOFF,
            MIN_DPHI,
            MIN_CSC_DNN,
            MIN_DT_DNN,
            # MIN_CSC_TIME,
            # MAX_CSC_TIME,
            # MAX_CSC_TSPREAD,
            # MAX_CSC_JET,
            # MAX_DT_JET,
            # MAX_CSC_MUON,
            # MAX_DT_MUON,
            # MAX_MB1,
            # HALO_CUTOFF,
            # MIN_DPHI,
            # MIN_CSC_DNN,
            # MIN_DT_DNN,
        ])

        with open(f'loo_tests_cscdt_{MET_CATEGORY}.pkl', 'wb') as fout:
            pickle.dump(tests, fout)
        
        print(f'{s2b:>5.1f}, {s:>4.0f}, {b:>6.0f} | {MIN_CSC_TIME:>5.1f}, {MAX_CSC_TIME:>4.2f}, {MAX_CSC_TSPREAD:>4.1f} | {MAX_RPC_BX:>1.0f}, {MIN_RPC_HITS:>2.0f} | {MAX_CSC_JET:>3.0f}, {MAX_DT_JET:>3.0f} | {MAX_CSC_MUON:>3.0f}, {MAX_DT_MUON:>3.0f} | {MAX_ME1:>2.0F}, {MAX_MB1:>2.0f} | {HALO_CUTOFF:>4.2f}, {MIN_DPHI:>4.2f} | {MIN_CSC_DNN:>4.3f}, {MIN_DT_DNN:>4.3f}')

        if i % 50 == 0 or i+1 == n_iterations:
            _t = np.array(tests)
            # idx_best = np.argmax(_t[:,3])
            print('')
            print(f'Completed {i+1} iterations')
            print('Best S/sqrt[ B ], B >= 3:')
            idxs = np.argsort(_t[:,0] * (_t[:,2]>=3))[::-1]
            for iidx, idx in enumerate(idxs[:10]):
                print(f'{iidx}: ' + ', '.join([f'{x:.2f}' for x in tests[idx]]))
            
            print('Best S/sqrt[ B ], S > 900, B >= 3:')
            idxs = np.argsort(_t[:,0] * (_t[:,1]>900) * (_t[:,2]>=3))[::-1]
            for iidx, idx in enumerate(idxs[:10]):
                print(f'{iidx}: ' + ', '.join([f'{x:.2f}' for x in tests[idx]]))
            print('')

        if n_iterations == 1:
            print('saving loo rdfs using the best cut selection found')
            # for key, rdf in rdfs.items():
            #     fout = f'data/processed/{key}_cscdt_loo_{MET_CATEGORY}_rdf.root'
            #     if 'DT OOT' in CUTS:
            #         fout = fout.replace('cscdt_', 'cscdtOOT_')
            #     rdf = rdf.Snapshot('MuonSystem_flat', fout, columns_out)