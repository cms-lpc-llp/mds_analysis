import os
# import sys
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
    'MATCH',
    'HLT',
    # 'CSC&DT>0',
    'L1',
    'CSCIT',
    'DTIT',
    # 'MET',
    # 'ME1',
    'MB1',
    'JET',
    # 'MUON',
    # 'NLEP',
    # 'BDT',
    'HALO',
    # 'CSCSIZE',
    'DTSTN',
    '1CSC1DT',
    # 'BLINDSR',
    # 'DR',
    'DPHI_0.4',
    # 'DPHI_0.2',
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
    print('Starting skim.py')
    print('----------------')

    # Arguments:
    #   - input path
    #   - output path
    #   - apply mc cuts
    #   - cut options (eg in-time vs out-of-time)
    #   - multithreading

    # ROOT MT causes the code to crash (not enough memory?)
    rt.EnableImplicitMT()
    print('  Enabled ROOT\'s implicit multithreading (sometimes causes a crash)')

    print('')

    # **************** #

    rdfs = {
        'mc' : rt.RDataFrame('MuonSystem', FN_MC),
        # 'r3' : rt.RDataFrame('MuonSystem', FN_R3),
        # 'mc' : rt.RDataFrame('MuonSystem_HLT569', '/home/psimmerl/mds_analysis/data/processed/mc_hlt569.root'),
        'r3' : rt.RDataFrame('MuonSystem_HLT569', '/home/psimmerl/mds_analysis/data/processed/r3_hlt569.root'),
    }
    # columns_in = [n for n in rdfs['mc'].GetColumnNames()]
    # for c in columns_in:
    #     print(c)

    print('Events In:')
    for key, rdf in rdfs.items():
        # rdf = rdf.Range(0,100_000)
        # rdfs[key] = rdf

        # Fix the weight branch for MC
        if 'mc' in key:
            rdf = rdf.Redefine('weight', f'weight * {LUMI}')
            rdfs[key] = rdf

        count, weight = rdf.Count().GetValue(), rdf.Sum('weight').GetValue()
        print(f'  {key} = {count:,} ({weight:,.2f})')
    print('')


    C, D = 'cscRechitCluster', 'dtRechitCluster'

    for key, rdf in rdfs.items():
        # Create dummy columns to store what cluster indices pass our selections
        rdf = rdf.Define(f'{C}Flag', f'{C}Size > 0')
        rdf = rdf.Define(f'{D}Flag', f'{D}Size > 0')
        # ec,cc,dc = rdf.Count(), rdf.Sum(f'{C}Flag'), rdf.Sum(f'{D}Flag')
        # ec,cc,dc=ec.GetValue(),cc.GetValue(),dc.GetValue()
        # print(f' 0 | {ec=:,}, {cc=:,}, {dc=:,}')

        # Add cluster radius columns
        for p in (C, D):
            rdf = rdf.Define(f'{p}R', f'ROOT::VecOps::sqrt( {p}X*{p}X + {p}Y*{p}Y )')
            if f'{p}R' not in COLUMNS_OUT:
                COLUMNS_OUT.append(f'{p}R')

        # Apply selections from CUTS (and/or arguments?)
        if 'MATCH' in CUTS and 'mc' in key:
            # print(key, 'MATCH')
            rdf = rdf.Redefine(f'{C}Flag', f'{C}Flag && ( {C}_match_gLLP && (abs({C}_match_gLLP_eta) < 3) && ({C}_match_gLLP_decay_r < 800) && (400 < abs({C}_match_gLLP_decay_z)) && (abs({C}_match_gLLP_decay_z) < 1200) )')
            rdf = rdf.Redefine(f'{D}Flag', f'{D}Flag && ( {D}_match_gLLP && (200 < {D}_match_gLLP_decay_r) && ({D}_match_gLLP_decay_r < 800) && (abs({D}_match_gLLP_decay_z) < 700) )')

        # ec,cc,dc = rdf.Count(), rdf.Sum(f'{C}Flag'), rdf.Sum(f'{D}Flag')
        # ec,cc,dc=ec.GetValue(),cc.GetValue(),dc.GetValue()
        # print(f' 1 | {ec=:,}, {cc=:,}, {dc=:,}')

        # Trigger selections (HLT might be wrong for MC so use the reproduced L1 trigger)
        # if 'HLT' in CUTS and 'mc' not in key:
        #     rdf = rdf.Filter('HLTDecision[569]') # 569 = HLT_L1CSCCluster_DTCluster50
        if 'L1' in CUTS:
            rdf = rdf.Redefine(f'{C}Flag', f'{C}Flag && ( '+
                               f'(({C}R > 100) && ({C}R < 275) && (abs({C}Z) > 580) && (abs({C}Z) < 632) && ({C}Size >= 500)) || ' + # ME 11
                               f'(({C}R > 275) && ({C}R < 465) && (abs({C}Z) > 668) && (abs({C}Z) < 724) && ({C}Size >= 200)) || ' + # ME 12
                               f'(({C}R > 505) && ({C}R < 700) && (abs({C}Z) > 668) && (abs({C}Z) < 724) && ({C}Size >= 200)) || ' + # ME 13
                               #
                               f'(({C}R > 139) && ({C}R < 345) && (abs({C}Z) > 789) && (abs({C}Z) < 850) && ({C}Size >= 500)) || ' + # ME 21
                               f'(({C}R > 357) && ({C}R < 700) && (abs({C}Z) > 791) && (abs({C}Z) < 850) && ({C}Size >= 200)) || ' + # ME 22
                               #
                               f'(({C}R > 160) && ({C}R < 345) && (abs({C}Z) > 915) && (abs({C}Z) < 970) && ({C}Size >= 500)) || ' + # ME 31
                               f'(({C}R > 357) && ({C}R < 700) && (abs({C}Z) > 911) && (abs({C}Z) < 970) && ({C}Size >= 200)) || ' + # ME 32
                               #
                               f'(({C}R > 178) && ({C}R < 345) && (abs({C}Z) > 1002) && (abs({C}Z) < 1063) && ({C}Size >= 500)) || ' + # ME 41
                               f'(({C}R > 357) && ({C}R < 700) && (abs({C}Z) > 1002) && (abs({C}Z) < 1063) && ({C}Size >= 200)) )') # ME 42
        # ec,cc,dc = rdf.Count(), rdf.Sum(f'{C}Flag'), rdf.Sum(f'{D}Flag')
        # ec,cc,dc=ec.GetValue(),cc.GetValue(),dc.GetValue()
        # print(f' 2 | {ec=:,}, {cc=:,}, {dc=:,}')

        if 'MET' in CUTS:
            rdf.Filter('met > 200')

        # Cluster level selections
        # if key == 'mc':
        #     rdf = rdf.Redefine(f'{C}Flag', f'{C}Flag && ( (-5 < {C}TimeWeighted) && ({C}TimeWeighted < 12.5) && ({C}TimeSpreadWeightedAll < 20) )')
        # else:
        #     rdf = rdf.Redefine(f'{C}Flag', f'{C}Flag && !( (-5 < {C}TimeWeighted) && ({C}TimeWeighted < 12.5) && ({C}TimeSpreadWeightedAll < 20) )')

        if 'CSCIT' in CUTS:
            rdf = rdf.Redefine(f'{C}Flag', f'{C}Flag && ( (-5 < {C}TimeWeighted) && ({C}TimeWeighted < 12.5) && ({C}TimeSpreadWeightedAll < 20) )')
        elif 'CSCOOT' in CUTS:
            rdf = rdf.Redefine(f'{C}Flag', f'{C}Flag && !( (-5 < {C}TimeWeighted) && ({C}TimeWeighted < 12.5) && ({C}TimeSpreadWeightedAll < 20) )')
        # ec,cc,dc = rdf.Count(), rdf.Sum(f'{C}Flag'), rdf.Sum(f'{D}Flag')
        # ec,cc,dc=ec.GetValue(),cc.GetValue(),dc.GetValue()
        # print(f' 3 | {ec=:,}, {cc=:,}, {dc=:,}')
        if 'DTIT' in CUTS:
            rdf = rdf.Redefine(f'{D}Flag', f'{D}Flag && ( ({D}_match_RPCBx_dPhi0p5 == 0) && ({D}_match_RPChits_dPhi0p5 > 0) )')
        elif 'DTOOT' in CUTS:
            rdf = rdf.Redefine(f'{D}Flag', f'{D}Flag && !( ({D}_match_RPCBx_dPhi0p5 == 0) && ({D}_match_RPChits_dPhi0p5 > 0) )')
            
        # ec,cc,dc = rdf.Count(), rdf.Sum(f'{C}Flag'), rdf.Sum(f'{D}Flag')
        # ec,cc,dc=ec.GetValue(),cc.GetValue(),dc.GetValue()
        # print(f' 3 | {ec=:,}, {cc=:,}, {dc=:,}')
 
        if 'ME1' in CUTS:
            rdf = rdf.Redefine(f'{C}Flag', f'{C}Flag && ( '+
                               f'({C}NRechitChamberPlus11 == 0) &&'+
                               f'({C}NRechitChamberPlus12 == 0) &&'+
                               f'({C}NRechitChamberMinus11 == 0) &&'+
                               f'({C}NRechitChamberMinus12 == 0) )')
        if 'MB1' in CUTS:
            rdf = rdf.Redefine(f'{D}Flag', f'{D}Flag && ( {D}NHitStation1 == 0 )')
        # ec,cc,dc = rdf.Count(), rdf.Sum(f'{C}Flag'), rdf.Sum(f'{D}Flag')
        # ec,cc,dc=ec.GetValue(),cc.GetValue(),dc.GetValue()
        # print(f' 4 | {ec=:,}, {cc=:,}, {dc=:,}')
        if 'JET' in CUTS:
            rdf = rdf.Redefine(f'{C}Flag', f'{C}Flag && ( ({C}JetVetoPt < 10) )')#30) || ({C}JetVetoLooseId == 0) )')
            rdf = rdf.Redefine(f'{D}Flag', f'{D}Flag && ( ({D}JetVetoPt < 10) )')#50) || ({D}JetVetoLooseId == 0) )')
        # ec,cc,dc = rdf.Count(), rdf.Sum(f'{C}Flag'), rdf.Sum(f'{D}Flag')
        # ec,cc,dc=ec.GetValue(),cc.GetValue(),dc.GetValue()
        # print(f' 5 | {ec=:,}, {cc=:,}, {dc=:,}')
        if 'MUON' in CUTS:
            rdf = rdf.Redefine(f'{C}Flag', f'{C}Flag && ( ({C}MuonVetoPt < 30) )') # || ({C}MuonVetoGlobal == 0) )')
            rdf = rdf.Redefine(f'{D}Flag', f'{D}Flag && ( ({D}MuonVetoPt < 50) )') # || ({D}MuonVetoLooseId == 0) )')

        if 'HALO' in CUTS:
            rdf = rdf.Redefine(f'{D}Flag', f'{D}Flag && ( (0.4 < abs({D}Phi)) && (abs({D}Phi) < {PI} - 0.4) )')
        # ec,cc,dc = rdf.Count(), rdf.Sum(f'{C}Flag'), rdf.Sum(f'{D}Flag')
        # ec,cc,dc=ec.GetValue(),cc.GetValue(),dc.GetValue()
        # print(f' 6 | {ec=:,}, {cc=:,}, {dc=:,}')
        if 'DTSTN' in CUTS:
            rdf = rdf.Redefine(f'{D}Flag', f'{D}Flag && ( ({D}NStation10 < 3) && !(({D}NStation10 == 2) && ({D}MaxStation == 4)) )')
        # ec,cc,dc = rdf.Count(), rdf.Sum(f'{C}Flag'), rdf.Sum(f'{D}Flag')
        # ec,cc,dc=ec.GetValue(),cc.GetValue(),dc.GetValue()
        # print(f' 7 | {ec=:,}, {cc=:,}, {dc=:,}')

        # if 'BDT' in CUTS:
        #     raise NotImplementedError('BDT')

        # Finally require that we still have just one or at least one pair of CSC+DT per event
        # rdf = rdf.Define('nCsc', f'reduce({C}Flag.begin(), {C}Flag.end())')
        # rdf = rdf.Define('nDt', f'reduce({D}Flag.begin(), {D}Flag.end())')

        # canvas = rt.TCanvas('','',800,800)
        # canvas.SetLogz()
        # hh = rdf.Histo2D((f'{key}_nCsc_nDt',f'{key};nCsc;nDt',6,-0.5, 5.5,6,-0.5,5.5),'nCsc','nDt').GetValue()
        # hh.SetMinimum(0)
        # hh.Draw('text colz')
        # canvas.Draw()
        # canvas.Print(f'{OUT_DIR}/{key}_nCsc_nDt.png')

        # for x in (C,D):
        #     for nx in range(0, 4):
        #         print(x, nx, rdf.Filter(f'(reduce({x}Flag.begin(), {x}Flag.end()) == {nx})').Count().GetValue())

        if '1CSC1DT' in CUTS:
            rdf = rdf.Filter(f'(reduce({C}Flag.begin(), {C}Flag.end()) == 1) && (reduce({D}Flag.begin(), {D}Flag.end()) == 1)')
        else:
            rdf = rdf.Filter(f'(reduce({C}Flag.begin(), {C}Flag.end()) > 0) && (reduce({D}Flag.begin(), {D}Flag.end()) > 0)')

        # ec,cc,dc = rdf.Count(), rdf.Sum(f'{C}Flag'), rdf.Sum(f'{D}Flag')
        # ec,cc,dc=ec.GetValue(),cc.GetValue(),dc.GetValue()
        # print(f' 8 | {ec=:,}, {cc=:,}, {dc=:,}')
            # raise NotImplementedError('cant handle multiple pairs yet')

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

        # if 'BLINDSR' in CUTS:
        #     raise NotImplementedError('BLINDSR')
        # if 'DR' in CUTS:
        #     raise NotImplementedError('DR')
        # if 'DETA' in CUTS: #! Has not been tested -- if gg->H->ss is exclusive then conservation
        #     rdf = rdf.Filter('((0 < cscEta) && ( dtEta < 0)) || (( cscEta < 0) && (0 < dtEta)) ||')
        if 'DPHI_0.4' in CUTS:
            rdf = rdf.Filter('tag_dPhi > 0.4')
            # raise NotImplementedError('DPHI')
        if 'DPHI_0.2' in CUTS:
            rdf = rdf.Filter('tag_dPhi > 0.2')
        
        # Update the dictionary with the new RDF
        rdfs[key] = rdf

    print('Events Out:')
    for key, rdf in rdfs.items():
        rdf = rdf.Snapshot('MuonSystem_flat', f'data/processed/{key}_rdf.root', columns_out)
        rdfs[key] = rdf

        count, weight = rdf.Count().GetValue(), rdf.Sum('weight').GetValue()
        print(f'  {key} = {count:,} ({weight:,.2f})')
        for xx in ('cscSize', 'cscR', 'cscEta', 'cscPhi', 'dtSize', 'dtR', 'dtEta', 'dtPhi', 'tag_dEta','tag_dPhi','tag_dR'):
            canvas = rt.TCanvas('','',800,800)
            xmin, xmax, std = rdf.Min(xx).GetValue(), rdf.Max(xx).GetValue(), rdf.StdDev(xx).GetValue()
            hh = rdf.Histo1D((f'{key}_{xx}',f'{key};{xx};count',int((xmax-xmin)/(2.718*std*count**(-1/3))),xmin,xmax),f'{xx}').GetValue()
            hh.SetMinimum(0)
            hh.Draw()
            canvas.Draw()
            canvas.Print(f'{OUT_DIR}/{key}_{xx}.png')
    print('')
