import os
import sys
import pathlib
import argparse

# import math
# import numpy as np
# import numba as nb
# from math import ceil, floor

# from src.helper_functions import alert, Table
# from src.histo_utilities import std_color_list as SCL

import ROOT as rt
from ROOT import RDataFrame
# from ROOT import TNtuple, TTree, TBranch, RDataFrame
# from ROOT import TCanvas, TLatex, TLegend, TLine, TBox
# from ROOT import TH1D, TH2D, TGraph, TGraphErrors, TGraphAsymmErrors

# **************************** #
CUTS = [
    'acceptance',
    'HLT',
    'L1',
    # 'low MET',
    'high MET',
    #! I reset cutflow indices here
    'CSC IT',
    'DT IT',
    # 'DT OOT',
    # 'DT I/OOT',
    # 'MET',
    # 'DNN',
    # 'BDT',
    # 'ME1',
    'MB1',
    'jet veto',
    # 'muon veto',
    # 'N Lep',
    # 'halo veto',
    # 'DT stn',
    # 'CSCSIZE',
    '1 CSC-DT',
    # 'BLINDSR',
    # 'DR',
    # 'dPhi_0.4',
    'dPhi_0.2',
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
    'DT stn',
    # 'DNN',
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

if __name__ == '__main__':
    print('+------------------+')
    print('| Starting skim.py |')
    print('+ -----------------+')


    if len(sys.argv) > 1:
        if sys.argv[1] == 'low':
            print('    Using low met cuts')
            CUTS = CUTS_LOW
        elif sys.argv[1] == 'high':
            print('    Using high met cuts')
            CUTS = CUTS_HIGH
        else:
            raise ValueError('Unknown arguments: '+' '.join(sys.argv[1:]))

    # Arguments:
    #   - input path
    #   - output path
    #   - apply mc cuts
    #   - cut options (eg in-time vs out-of-time)
    #   - multithreading

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
    print('Events in:')
    
    rdfs = {
        'mc' : rt.RDataFrame(TN_MC, FN_MC),
        'r3' : rt.RDataFrame(TN_R3, FN_R3),
    }

    for key, rdf in rdfs.items():
        # rdf = rdf.Range(0,100_000) # Skim a subset of events for debugging
        count, weight = rdf.Count().GetValue(), rdf.Sum('weight').GetValue()
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

        # **** #
        # Fix the weight branch
        if 'mc' in key:
            rdf = rdf.Redefine('weight', f'weight * {LUMI}')

        # rdf = rdf.Filter(f'(reduce({C}Flag.begin(), {C}Flag.end()) > 0) && (reduce({D}Flag.begin(), {D}Flag.end()) > 0)')
        # count, weight = rdf.Count().GetValue(), rdf.Sum('weight').GetValue()
        # print(f'  {key} = {count:,} ({weight:,.2f}) -- has at least 1 csc and 1 dt')

        # **** #
        rdfs[key] = rdf

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
            print(f'    sample & {ec0:,.0f} & -- & -- & '+
                  f'{cc0:,} & -- & -- & '+
                  f'{dc0:,} & -- & -- \\\\')

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
                rdf = rdf.Redefine(f'{C}CutFlag', f'{C}CutFlag && ( (-5 < {C}TimeWeighted) && ({C}TimeWeighted < 12.5) && ({C}TimeSpreadWeightedAll < 20) )')
            elif 'CSC OOT' in cut:
                rdf = rdf.Redefine(f'{C}CutFlag', f'{C}CutFlag && !( (-5 < {C}TimeWeighted) && ({C}TimeWeighted < 12.5) && ({C}TimeSpreadWeightedAll < 20) )')

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
                if 'low MET' in CUTS:
                    rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( {D}NHitStation1 == 0 )')
                else:
                    rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( {D}NHitStation1 <= 3 )')

            if 'DT stn' in cut:
                rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( ({D}NStation10 < 3) && !(({D}NStation10 == 2) && ({D}MaxStation == 4)) )')

            # **** #
            if 'jet veto' in cut:
                if 'high MET' in CUTS:
                    rdf = rdf.Redefine(f'{C}CutFlag', f'{C}CutFlag && ( ({C}JetVetoLooseId == 0) || ({C}JetVetoPt < 30) )')
                    rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( ({D}JetVetoLooseId == 0) || ({D}JetVetoPt < 50) )')
                else:
                    rdf = rdf.Redefine(f'{C}CutFlag', f'{C}CutFlag && ( ({C}JetVetoLooseId == 0) || ({C}JetVetoPt < 10) )')
                    rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( ({D}JetVetoLooseId == 0) || ({D}JetVetoPt < 10) )')
            
            if 'muon veto' in cut:
                rdf = rdf.Redefine(f'{C}CutFlag', f'{C}CutFlag && ( ({C}MuonVetoGlobal == 0) || ({C}MuonVetoPt < 30) )')
                rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( ({D}MuonVetoLooseId == 0) || ({D}MuonVetoPt < 50) )')

            if 'halo veto' in cut:
                rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( (0.4 < abs({D}Phi)) && (abs({D}Phi) < {PI} - 0.4) )')

            # **** #
            # if 'BDT' in cut:
            #     raise NotImplementedError('BDT')

            if 'DNN' in cut:
                # rdf = rdf.Redefine(f'{C}CutFlag', f'{C}CutFlag && ( Take({C}DNN,nCscRechitClusters) > 0.85 )')
                # rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( Take({D}DNN,nDtRechitClusters) > 0.85 )')
                if 'high MET' in CUTS:
                    rdf = rdf.Redefine(f'{C}CutFlag', f'{C}CutFlag && ( Take({C}DNN,nCscRechitClusters) > 0.01 )')
                    rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( Take({D}DNN,nDtRechitClusters) > 0.01 )')
                else:
                    rdf = rdf.Redefine(f'{C}CutFlag', f'{C}CutFlag && ( Take({C}DNN,nCscRechitClusters) > 0.75 )')
                    rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( Take({D}DNN,nDtRechitClusters) > 0.75 )')

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
                if '_' in cut:
                    rdf = rdf.Redefine('evtCutFlag', f'evtCutFlag && (tag_dPhi > {cut.split("_")[-1]})')
                else:
                    rdf = rdf.Redefine('evtCutFlag', f'evtCutFlag && (tag_dPhi > 0.4)')

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
            if PRINT_CUTFLOW:
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
                
                if 'MET' in cut:
                    print(r'    \hline')
                    ec0, cc0, dc0 = ec, cc, dc
                    
                    rdf = rdf.Filter('evtFlag')
                    # for col in rdf.GetColumnNames(): #! this causes an error where the length of the RVec is wrong
                    #     col = str(col)
                    #     if C in col[:len(C)]:
                    #         rdf = rdf.Redefine(col, f'{col}[{C}Flag]')
                    #     elif D in col[:len(D)]:
                    #         rdf = rdf.Redefine(col, f'{col}[{D}Flag]')

                    # rdf = rdf.Redefine('weight', 'weight * evtFlag')
                    rdf = rdf.Redefine(f'{C}Size', f'{C}Size * {C}Flag')
                    rdf = rdf.Redefine(f'{D}Size', f'{D}Size * {D}Flag')
            else:
                rdf = rdf.Filter('evtFlag')
                rdf = rdf.Redefine(f'{C}Size', f'{C}Size * {C}Flag')
                rdf = rdf.Redefine(f'{D}Size', f'{D}Size * {D}Flag')


        rdfs[key] = rdf.Filter('evtFlag') # Apply event level cut and update the dictionary with the new RDF

        if PRINT_CUTFLOW:
            print(r'    \hline')
            print(r'\end{tabular}')
            print(r'\end{center}')
            print('')

    # **************** #
    print('Events out:')
    for key, rdf in rdfs.items():
        if rdf.Count().GetValue() == 0:
            print(f'{key} is empty')
            continue

        rdf = rdf.Snapshot('MuonSystem_flat', f'data/processed/{key}_rdf.root', columns_out)
        rdfs[key] = rdf

        count, weight = rdf.Count().GetValue(), rdf.Sum('weight').GetValue()
        print(f'  {key} = {count:,} ({weight:,.2f})')
        for xx in ('met', 'cscSize', 'cscR', 'cscEta', 'cscPhi', 'dtSize', 'dtR', 'dtEta', 'dtPhi', 'tag_dEta','tag_dPhi','tag_dR'):
            canvas = rt.TCanvas('','',800,800)
            xmin, xmax, std = rdf.Min(xx).GetValue(), rdf.Max(xx).GetValue(), rdf.StdDev(xx).GetValue()
            nbins = int((xmax-xmin)/(2.718*std*count**(-1/3)))
            hh = rdf.Histo1D((f'{key}_{xx}',f'{key};{xx};count',nbins,xmin,xmax),f'{xx}').GetValue()
            hh.SetMinimum(0)
            hh.Draw()
            canvas.Draw()
            canvas.Print(f'{OUT_DIR}/{key}_{xx}.png')
    print('')
