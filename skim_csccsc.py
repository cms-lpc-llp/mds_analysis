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
    'MET',
    #! I reset cutflow indices here
    'CSC0 IT',
    'CSC1 IT',
    # 'CSC I/OOT',
    # 'CSC OOT',
    'DT IT',
    'ME1',
    'MB1',
    'jet veto',
    'muon veto',
    # 'no leptons',
    # 'halo veto',
    # 'DT stn',
    # 'CSCSIZE',
    # 'DNN',
    '1 CSC-CSC',
    # 'BLINDSR',
    # 'DR',
    'dPhi_0.4',
]


PRINT_CUTFLOW = True
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
    'cscRechitClusterR',
    #
    # 'ctau',
    #
    # 'dtRechitClusterAvgStation10',
    # 'dtRechitClusterEta',
    # 'dtRechitClusterJetVetoE',
    # 'dtRechitClusterJetVetoLooseId',
    # 'dtRechitClusterJetVetoPt',
    # 'dtRechitClusterJetVetoTightId',
    # 'dtRechitClusterMaxChamber',
    # 'dtRechitClusterMaxChamberRatio',
    # 'dtRechitClusterMaxStation',
    # 'dtRechitClusterMaxStationRatio',
    # 'dtRechitClusterMet_dPhi',
    # 'dtRechitClusterMuonVetoE',
    # 'dtRechitClusterMuonVetoGlobal',
    # 'dtRechitClusterMuonVetoLooseId',
    # 'dtRechitClusterMuonVetoPt',
    # 'dtRechitClusterMuonVetoTightId',
    # 'dtRechitClusterNChamber',
    # 'dtRechitClusterNHitStation1',
    # 'dtRechitClusterNHitStation2',
    # 'dtRechitClusterNHitStation3',
    # 'dtRechitClusterNHitStation4',
    # 'dtRechitClusterNOppositeSegStation1',
    # 'dtRechitClusterNOppositeSegStation2',
    # 'dtRechitClusterNOppositeSegStation3',
    # 'dtRechitClusterNOppositeSegStation4',
    # 'dtRechitClusterNSegStation1',
    # 'dtRechitClusterNSegStation2',
    # 'dtRechitClusterNSegStation3',
    # 'dtRechitClusterNSegStation4',
    # 'dtRechitClusterNStation10',
    # 'dtRechitClusterNoiseHit',
    # 'dtRechitClusterNoiseHitStation1',
    # 'dtRechitClusterNoiseHitStation2',
    # 'dtRechitClusterNoiseHitStation3',
    # 'dtRechitClusterNoiseHitStation4',
    # 'dtRechitClusterOverlap',
    # 'dtRechitClusterPhi',
    # 'dtRechitClusterSize',
    # 'dtRechitClusterWheel',
    # 'dtRechitClusterX',
    # 'dtRechitClusterY',
    # 'dtRechitClusterZ',
    # 'dtRechitCluster_match_MB1Seg_0p4',
    # 'dtRechitCluster_match_MB1Seg_0p5',
    # 'dtRechitCluster_match_MB1hits_0p4',
    # 'dtRechitCluster_match_MB1hits_0p5',
    # 'dtRechitCluster_match_MB1hits_cosmics_minus',
    # 'dtRechitCluster_match_MB1hits_cosmics_plus',
    # 'dtRechitCluster_match_RB1_0p4',
    # 'dtRechitCluster_match_RB1_dPhi0p5',
    # 'dtRechitCluster_match_RPCBx_dPhi0p5',
    # 'dtRechitCluster_match_RPChits_dPhi0p5',
    # # 'dtRechitCluster_match_gLLP',
    # # 'dtRechitCluster_match_gLLP_csc',
    # # 'dtRechitCluster_match_gLLP_decay_r',
    # # 'dtRechitCluster_match_gLLP_decay_z',
    # # 'dtRechitCluster_match_gLLP_dt',
    # # 'dtRechitCluster_match_gLLP_e',
    # # 'dtRechitCluster_match_gLLP_eta',
    # # 'dtRechitCluster_match_gLLP_index',
    # # 'dtRechitCluster_match_gLLP_minDeltaR',
    # # 'dtRechitCluster_match_gLLP_phi',
    #
    # 'dtRechitClusterR',
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
        if sys.argv[1] == 'it':
            print('    Using in-time 2nd cluster')
        elif sys.argv[1] == 'oot':
            print('    Using out-of-time 2nd cluster')
            CUTS = [c.replace('CSC1 IT', 'CSC1 OOT') if 'CSC1 IT' == c else c for c in CUTS]
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
    # FN_R3, TN_R3 ='/home/psimmerl/mds_analysis/data/raw/DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi_v6.root', 'MuonSystem'

    # # My skim where it passes (HLT 566)
    # FN_MC, TN_MC ='/home/psimmerl/mds_analysis/data/processed/mc_hlt566.root', 'MuonSystem'
    # FN_R3, TN_R3 ='/home/psimmerl/mds_analysis/data/processed/r3_hlt566.root', 'MuonSystem'

    # Pedro's DNN skim
    LUMI = 1.1328524540090597e-06 * LUMI # ???
    # FN_MC, TN_MC ='/home/psimmerl/mds_analysis/data/processed/mc_pedro.root', 'MuonSystem'
    # FN_R3, TN_R3 ='/home/psimmerl/mds_analysis/data/processed/data_pedro.root', 'MuonSystem'
    # FN_R3, TN_R3 ='/home/psimmerl/mds_analysis/data/processed/r3_pedro_hlt566.root', 'MuonSystem'
    FN_MC, TN_MC ='/home/psimmerl/mds_analysis/data/processed/mc_pedro.root', 'MuonSystem'
    FN_R3, TN_R3 ='/home/psimmerl/mds_analysis/data/processed/r3_pedro_hlt566.root', 'MuonSystem'

    if '1 CSC-CSC' not in CUTS:
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
        rdf = rdf.Define(f'{C}0Size', f'{C}Size')
        rdf = rdf.Define(f'{C}1Size', f'{C}Size')

        rdf = rdf.Define('evtFlag', 'weight > 0')
        rdf = rdf.Define(f'{C}0Flag', f'{C}0Size > 0')
        rdf = rdf.Define(f'{C}1Flag', f'{C}1Size > 0')
        rdf = rdf.Define(f'{D}Flag', f'{D}Size > 0')

        rdf = rdf.Define('evtCutFlag', 'weight > 0')
        rdf = rdf.Define(f'{C}0CutFlag', f'{C}0Size > 0')
        rdf = rdf.Define(f'{C}1CutFlag', f'{C}1Size > 0')
        rdf = rdf.Define(f'{D}CutFlag', f'{D}Size > 0')

        # **** #
        # Add cluster radius columns
        for p in (C, D):
            rdf = rdf.Define(f'{p}R', f'ROOT::VecOps::sqrt( {p}X*{p}X + {p}Y*{p}Y )')

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
    columns_out = []

    for key, rdf in rdfs.items():
        if PRINT_CUTFLOW:
            print(f'{key} cutflow:')
            print(r'\begin{center}')
            print(r'\begin{tabular}{c|ccc|ccc|ccc|ccc}')
            ec0,cc00,cc10,dc0 = rdf.Filter('evtFlag').Sum('weight'), rdf.Filter('evtFlag').Sum(f'{C}0Flag'), rdf.Filter('evtFlag').Sum(f'{C}1Flag'), rdf.Filter('evtFlag').Sum(f'{D}Flag')
            ec0,cc00,cc10,dc0 = ec0.GetValue(),cc00.GetValue(),cc10.GetValue(),dc0.GetValue()
            print(r'    \hline')
            print(f'    \\textbf{{{"Signal" if "mc" in key else "Data"}}}'+r' & \multicolumn{3}{c}{CSC-CSC} & \multicolumn{3}{c}{Trigger CSC} & \multicolumn{3}{c}{Second CSC} & \multicolumn{3}{c}{DT} \\')
            print(r'    \hline')
            print(r'    Selection & N Events & Cut Eff & Cum Eff & N Clusters & Cut Eff & Cum Eff & N Clusters & Cut Eff & Cum Eff & N Clusters & Cut Eff & Cum Eff \\')
            print(r'    \hline')
            print(f'    sample & {ec0:,.0f} & -- & -- & '+
                  f'{cc00:,} & -- & -- & '+
                  f'{cc10:,} & -- & -- & '+
                  f'{dc0:,} & -- & -- \\\\')

        for cut in CUTS:
            # **** #
            # Reset cut flags
            rdf = rdf.Redefine('evtCutFlag', 'weight > 0')
            rdf = rdf.Redefine(f'{C}0CutFlag', f'{C}0Size > 0')
            rdf = rdf.Redefine(f'{C}1CutFlag', f'{C}1Size > 0')
            rdf = rdf.Redefine(f'{D}CutFlag', f'{D}Size > 0')

            # **** #
            # Signal matching for MC
            if 'acceptance' in cut: 
                if 'mc' in key:
                    rdf = rdf.Redefine(f'{C}0CutFlag', f'{C}0CutFlag && {C}_match_gLLP && '+
                                    f'(abs({C}_match_gLLP_eta) < 3) && '+
                                    f'({C}_match_gLLP_decay_r < 800) && '+
                                    f'(400 < abs({C}_match_gLLP_decay_z)) &&'+
                                    f' (abs({C}_match_gLLP_decay_z) < 1200)')
                    rdf = rdf.Redefine(f'{C}1CutFlag', f'{C}1CutFlag && {C}_match_gLLP && '+
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
                    rdf = rdf.Redefine('evtCutFlag', 'evtCutFlag && HLTDecision[566]') # 566 = CscCluster_Loose
                else:
                    continue

            if 'L1' in cut: # First passes, second fails
                rdf = rdf.Redefine(f'{C}0CutFlag', f'{C}0CutFlag && ( '+
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

                rdf = rdf.Redefine(f'{C}1CutFlag', f'{C}1CutFlag && !( '+
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
                rdf = rdf.Redefine('evtCutFlag', f'evtCutFlag && (150 < met)')

            # **** #
            # Cluster level selections
            if 'CSC I/OOT' in cut:
                rdf = rdf.Redefine(f'{C}0CutFlag', f'{C}0CutFlag && ( (-5 < {C}TimeWeighted) && ({C}TimeWeighted < 12.5) && ({C}TimeSpreadWeightedAll < 20) )')
                rdf = rdf.Redefine(f'{C}1CutFlag', f'{C}1CutFlag && !( (-5 < {C}TimeWeighted) && ({C}TimeWeighted < 12.5) && ({C}TimeSpreadWeightedAll < 20) )')

            if 'CSC0 IT' in cut:
                rdf = rdf.Redefine(f'{C}0CutFlag', f'{C}0CutFlag && ( (-5 < {C}TimeWeighted) && ({C}TimeWeighted < 12.5) && ({C}TimeSpreadWeightedAll < 20) )')
            if 'CSC1 IT' in cut:
                rdf = rdf.Redefine(f'{C}1CutFlag', f'{C}1CutFlag && ( (-5 < {C}TimeWeighted) && ({C}TimeWeighted < 12.5) && ({C}TimeSpreadWeightedAll < 20) )')
            
            if 'CSC0 OOT' in cut:
                rdf = rdf.Redefine(f'{C}0CutFlag', f'{C}0CutFlag && !( (-5 < {C}TimeWeighted) && ({C}TimeWeighted < 12.5) && ({C}TimeSpreadWeightedAll < 20) )')
            if 'CSC1 OOT' in cut:
                rdf = rdf.Redefine(f'{C}1CutFlag', f'{C}1CutFlag && !( (-5 < {C}TimeWeighted) && ({C}TimeWeighted < 12.5) && ({C}TimeSpreadWeightedAll < 20) )')

            if 'DT IT' in cut:
                rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( ({D}_match_RPCBx_dPhi0p5 == 0) && ({D}_match_RPChits_dPhi0p5 > 0) )')
            elif 'DT OOT' in cut:
                rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && !( ({D}_match_RPCBx_dPhi0p5 == 0) && ({D}_match_RPChits_dPhi0p5 > 0) )')

            # **** #
            if 'ME1' in cut:
                rdf = rdf.Redefine(f'{C}1CutFlag', f'{C}1CutFlag && ( '+
                                   f'({C}NRechitChamberPlus11 == 0) &&'+
                                   f'({C}NRechitChamberPlus12 == 0) &&'+
                                   f'({C}NRechitChamberMinus11 == 0) &&'+
                                   f'({C}NRechitChamberMinus12 == 0) )')

            if 'MB1' in cut:
                rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( {D}NHitStation1 == 0 )')
                # rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( {D}NHitStation1 <= 3 )')

            # **** #
            if 'jet veto' in cut:
                rdf = rdf.Redefine(f'{C}0CutFlag', f'{C}0CutFlag && ( ({C}JetVetoLooseId == 0) || ({C}JetVetoPt < 30) )')
                rdf = rdf.Redefine(f'{C}1CutFlag', f'{C}1CutFlag && ( ({C}JetVetoLooseId == 0) || ({C}JetVetoPt < 30) )')
                rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( ({D}JetVetoLooseId == 0) || ({D}JetVetoPt < 50) )')
            
            if 'muon veto' in cut:
                rdf = rdf.Redefine(f'{C}0CutFlag', f'{C}0CutFlag && ( ({C}MuonVetoGlobal == 0) || ({C}MuonVetoPt < 30) )')
                rdf = rdf.Redefine(f'{C}1CutFlag', f'{C}1CutFlag && ( ({C}MuonVetoGlobal == 0) || ({C}MuonVetoPt < 30) )')
                rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( ({D}MuonVetoLooseId == 0) || ({D}MuonVetoPt < 10) )')

            if 'halo veto' in cut:
                rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( (0.1 < abs({D}Phi)) && (abs({D}Phi) < {PI} - 0.1) )')

            # **** #
            # if 'BDT' in cut:
            #     raise NotImplementedError('BDT')

            if 'DNN' in cut:
                rdf = rdf.Redefine(f'{C}0CutFlag', f'{C}0CutFlag && ( Take({C}DNN,nCscRechitClusters) > 0.01 )')
                rdf = rdf.Redefine(f'{C}1CutFlag', f'{C}1CutFlag && ( Take({C}DNN,nCscRechitClusters) > 0.01 )')
                rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( Take({D}DNN,nDtRechitClusters) > 0.01 )')

            # **** #
            if '1 CSC-CSC' in cut:
                # rdf = rdf.Redefine('evtCutFlag', f'evtCutFlag && evtFlag && '+
                #                    f'(reduce({C}0Flag.begin(), {C}0Flag.end()) == 1) && '+
                #                    f'(reduce({C}1Flag.begin(), {C}1Flag.end()) == 1)')
                rdf = rdf.Redefine('evtCutFlag', f'evtCutFlag && evtFlag && (reduce({D}Flag.begin(), {D}Flag.end()) == 0)')
                rdf = rdf.Redefine(f'{C}0Flag', f'{C}0Flag && ( {C}Size == Max({C}Size*{C}0Flag) )')
                rdf = rdf.Redefine(f'{C}1Flag', f'{C}1Flag && ( {C}Size == Max({C}Size*{C}1Flag) )')
                # rdf = rdf.Redefine(f'{D}Flag', f'{D}Flag && ( {D}Size == Max({D}Size*{D}Flag) )')
                

                # Apply our cluster level selections to relevant columns
                for col in COLUMNS_OUT:
                    if C in col[:len(C)]:
                        ocol0 = col.replace('RechitCluster', '0')
                        ocol1 = col.replace('RechitCluster', '1')
                        rdf = rdf.Define(ocol0, f'{col}[{C}0Flag][0]')
                        rdf = rdf.Define(ocol1, f'{col}[{C}1Flag][0]')
                        if ocol0 not in columns_out:
                            columns_out.append(ocol0)
                        if ocol1 not in columns_out:
                            columns_out.append(ocol1)
                    elif D in col[:len(D)]:
                        ocol = col.replace('RechitCluster', '')
                        rdf = rdf.Define(ocol, f'{col}[{D}Flag][0]')
                        if ocol not in columns_out:
                            columns_out.append(ocol)
                    else:
                        if col not in columns_out:
                            columns_out.append(col)

                rdf = rdf.Define('csc0CTau', f'gLLP_ctau[{C}_match_gLLP_index[{C}0Flag][0]]')
                rdf = rdf.Define('csc1CTau', f'gLLP_ctau[{C}_match_gLLP_index[{C}1Flag][0]]')
                rdf = rdf.Define('tag_dEta', 'abs(csc0Eta - csc1Eta)')
                rdf = rdf.Define('tag_dPhi', f'abs( abs(csc0Phi - csc1Phi) - (abs(csc0Phi - csc1Phi) > {PI} ? 2*{PI} : 0) )')
                rdf = rdf.Define('tag_dR', 'sqrt(tag_dEta*tag_dEta + tag_dPhi*tag_dPhi)')
                for ocol in ['tag_dEta', 'tag_dPhi', 'tag_dR', 'csc0CTau', 'csc1CTau']:
                    if ocol not in columns_out:
                        columns_out.append(ocol)

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
            # Finally require that we still have just one or at least one pair of CSC+CSC per event
            rdf = rdf.Redefine('evtCutFlag', 'evtCutFlag && ('+
                               f'(reduce({C}0CutFlag.begin(), {C}0CutFlag.end()) > 0) && '+
                               f'(reduce({C}1CutFlag.begin(), {C}1CutFlag.end()) > 0))')

            # **** #
            # Propagate to cumulative flags
            rdf = rdf.Redefine('evtFlag', 'evtFlag && evtCutFlag')
            rdf = rdf.Redefine(f'{C}0Flag', f'{C}0Flag && {C}0CutFlag')
            rdf = rdf.Redefine(f'{C}1Flag', f'{C}1Flag && {C}1CutFlag')
            rdf = rdf.Redefine(f'{D}Flag', f'{D}Flag && {D}CutFlag')

            rdf = rdf.Redefine('evtFlag', 'evtFlag && ('+
                               f'(reduce({C}0Flag.begin(), {C}0Flag.end()) > 0) && '+
                               f'(reduce({C}1Flag.begin(), {C}1Flag.end()) > 0))')

            # rdf = rdf.Redefine(f'nC{C[1:]}s', f'reduce({C}Flag.begin(), {C}Flag.end())')
            # rdf = rdf.Redefine(f'nD{D[1:]}s', f'reduce({D}Flag.begin(), {D}Flag.end())')

            # **** #
            if PRINT_CUTFLOW:
                ec,cc0,cc1,dc = rdf.Filter('evtFlag').Sum('weight'), rdf.Filter('evtFlag').Sum(f'{C}0Flag'), rdf.Filter('evtFlag').Sum(f'{C}1Flag'), rdf.Filter('evtFlag').Sum(f'{D}Flag')
                # ecc,ccc0,ccc1,dcc = rdf.Filter('evtCutFlag').Sum('weight'), rdf.Filter('evtCutFlag').Sum(f'{C}0CutFlag'), rdf.Filter('evtCutFlag').Sum(f'{C}1CutFlag'), rdf.Filter('evtCutFlag').Sum(f'{D}CutFlag')
                ecc,ccc0,ccc1,dcc = rdf.Filter('evtCutFlag').Sum('weight'), rdf.Sum(f'{C}0CutFlag'), rdf.Sum(f'{C}1CutFlag'), rdf.Sum(f'{D}CutFlag')
                
                ec,cc0,cc1,dc = ec.GetValue(),cc0.GetValue(),cc1.GetValue(),dc.GetValue()
                ecc,ccc0,ccc1,dcc = ecc.GetValue(),ccc0.GetValue(),ccc1.GetValue(),dcc.GetValue()

                c0_cum_eff = f'{100*cc0/cc00:.2f}' if cc0 != cc00 else '--'
                c0_cut_eff = f'{100*ccc0/cc00:.2f}' if ccc0 != cc00 else '--'
                c1_cum_eff = f'{100*cc1/cc10:.2f}' if cc1 != cc10 else '--'
                c1_cut_eff = f'{100*ccc1/cc10:.2f}' if ccc1 != cc10 else '--'
                d_cum_eff = f'{100*dc/dc0:.2f}' if dc != dc0 else '--'
                d_cut_eff = f'{100*dcc/dc0:.2f}' if dcc != dc0 else '--'

                print(f'    {cut.replace("_"," ")} & {ec:,.0f} & {100*ecc/ec0:.2f} & {100*ec/ec0:.2f} & '+
                      f'{cc0:,.0f} & {c0_cut_eff} & {c0_cum_eff} & '+
                      f'{cc1:,.0f} & {c1_cut_eff} & {c1_cum_eff} & '+
                      f'{dc:,.0f} & {d_cut_eff} & {d_cum_eff} \\\\')
                
                if 'MET' in cut:
                    print(r'    \hline')
                    ec0, cc00, cc10, dc0 = ec, cc0, cc1, dc
                    
                    rdf = rdf.Filter('evtFlag')
                    # for col in rdf.GetColumnNames(): #! this causes an error where the length of the RVec is wrong
                    #     col = str(col)
                    #     if C in col[:len(C)]:
                    #         rdf = rdf.Redefine(col, f'{col}[{C}Flag]')
                    #     elif D in col[:len(D)]:
                    #         rdf = rdf.Redefine(col, f'{col}[{D}Flag]')

                    # rdf = rdf.Redefine('weight', 'weight * evtFlag')
                    rdf = rdf.Redefine(f'{C}0Size', f'{C}0Size * {C}0Flag')
                    rdf = rdf.Redefine(f'{C}1Size', f'{C}1Size * {C}1Flag')
                    rdf = rdf.Redefine(f'{D}Size', f'{D}Size * {D}Flag')
            else:
                rdf = rdf.Filter('evtFlag')
                rdf = rdf.Redefine(f'{C}0Size', f'{C}0Size * {C}0Flag')
                rdf = rdf.Redefine(f'{C}1Size', f'{C}1Size * {C}1Flag')
                rdf = rdf.Redefine(f'{D}Size', f'{D}Size * {D}Flag')


        rdfs[key] = rdf.Filter('evtFlag') # Apply event level cut and update the dictionary with the new RDF

        if PRINT_CUTFLOW:
            print(r'    \hline')
            print(r'\end{tabular}')
            print(r'\end{center}')
            print('')

        all_cols = [str(x) for x in rdf.GetColumnNames()]
        for ocol in columns_out:
            if ocol not in all_cols:
                print(ocol)
    # **************** #
    print('Events out:')
    for key, rdf in rdfs.items():
        if rdf.Count().GetValue() == 0:
            print(f'{key} is empty')
            continue

        rdf = rdf.Snapshot('MuonSystem_flat', f'data/processed/{key}_csccsc_rdf.root', columns_out)
        rdfs[key] = rdf

        count, weight = rdf.Count().GetValue(), rdf.Sum('weight').GetValue()
        print(f'  {key} = {count:,} ({weight:,.2f})')
        for xx in ('met', 'csc0Size', 'csc0R', 'csc0Eta', 'csc0Phi', 'csc1Size', 'csc1R', 'csc1Eta', 'csc1Phi', 'tag_dEta','tag_dPhi','tag_dR'):
            canvas = rt.TCanvas('','',800,800)
            xmin, xmax, std = rdf.Min(xx).GetValue(), rdf.Max(xx).GetValue(), rdf.StdDev(xx).GetValue()
            nbins = int((xmax-xmin)/(2.718*std*count**(-1/3)))
            hh = rdf.Histo1D((f'{key}_{xx}',f'{key};{xx};count',nbins,xmin,xmax),f'{xx}').GetValue()
            hh.SetMinimum(0)
            hh.Draw()
            canvas.Draw()
            canvas.Print(f'{OUT_DIR}/{key}_csccsc_{xx}.png')
    print('')
