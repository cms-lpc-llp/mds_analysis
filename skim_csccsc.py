import sys
import pathlib
import numpy as np

import ROOT as rt
from ROOT import RDataFrame

# **************************** #
LOCAL_DIR = '/eos/user/f/fernanpe/mds_analysis/'
OUT_DIR = f'{LOCAL_DIR}/reports/weekly/2024-04-15'

# STAT = 'raw'
# LUMI = 23.02 * 1000
# FN_MC = f'{LOCAL_DIR}/data/raw/ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v6.root'
# FN_R3 = f'{LOCAL_DIR}/data/raw/DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi_v6.root'

STAT = 'pedro'
LUMI = 1.1328524540090597e-06 * 27.82 * 1000 # ???

FN_MC = f'{LOCAL_DIR}/data/processed/mc_pedro_hlt566_2023.root'
FN_R3 = f'{LOCAL_DIR}/data/processed/r3_pedro_hlt566_2023.root'

# **** #
LOW_MET_CUTOFF = 150
HIGH_MET_CUTOFF = 150

# **** #
CUTS_L1 = [
    'acceptance',
    'HLT',
    'L1',
    'MET',
    #! I reset cutflow indices here
    'CSC0 IT',
    'CSC1 IT',
    'DT IT',
    'ME1',
    'MB1',
    '1 CSC-CSC',
    # 'dPhi',
]

CUTS = [
    'acceptance',
    'HLT',
    'L1',
    'MET',
    #! I reset cutflow indices here
    'CSC0 IT',
    'CSC1 IT',
    #'DT IT',
    'ME1',
    #'MB1',
    # 'no leptons',
    # 'halo veto',
    # 'DT stn',
    # 'BDT',
    '1 CSC-CSC',
    # 'dR',
    # 'dEta',
    'dPhi $>$ 1.8',
    'DNN $>$ 0.96',
    # TO BE REMOVED? vetoes doesn't seem to help s/sqrt(b)
    'muon veto',
    'jet veto',
    

]
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
    #'cscRechitClusterDNN_bkgMC',
    #'cscRechitClusterDNN_bkgMC_plusBeamHalo',
    #'cscRechitClusterDNN_bkgOOTData',
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

# **** #
pathlib.Path(OUT_DIR).mkdir(parents=True, exist_ok=True)

rt.gErrorIgnoreLevel = 1001  # rt.kInfo + 1
rt.gROOT.SetBatch(True)
PI = rt.TMath.Pi()

gc = []

# **************************** #

if __name__ == '__main__':
    print('+-------------------------+')
    print('| Starting skim_csccsc.py |')
    print('+ ------------------------+')

    rt.EnableImplicitMT(4)
    print('    Enabled ROOT\'s implicit multithreading (sometimes causes a crash)')

    # **************************** #
    PRINT_CUTFLOW = False
    MET_CATEGORY = 'lt200'
    CUTSET = 'scs'
    OOT = False

    args = ' '.join(sys.argv[1:]) if len(sys.argv) > 1 else ''

    if 'CUTFLOW' in args:
        print('    Printing cutflow tables')
        PRINT_CUTFLOW = True

    if 'l1' in args:
        print('    Using the reduced cut set (l1)')
        CUTS = CUTS_L1
        CUTSET = 'l1'
    elif 'ropt' in args:
        print('    Using the randomly optimized set (ropt)')
        CUTSET = 'ropt'
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
        CUTS = [c.replace('CSC1 IT', 'CSC1 OOT') if 'CSC1 IT' == c else c for c in CUTS]
        OOT = True
    else:
        print('    Using in-time 2nd cluster')

    print('')

    # **** #
    if '1 CSC-CSC' not in CUTS:
        raise NotImplementedError('cant handle multiple pairs yet')

    # **************************** #
    if CUTSET == 'l1': # Level-1 & HLT Cuts only
        MIN_CSC_TIME = -5.0
        MAX_CSC_TIME =  12.5
        MAX_CSC_TSPREAD =  20.0

        MAX_RPC_BX = 0
        MIN_RPC_HITS = 1

        MAX_CSC_JET = 999
        MAX_DT_JET = 999
        MAX_CSC_MUON = 999
        MAX_DT_MUON = 999

        MAX_ME1 = 0
        MAX_MB1 = 0

        HALO_CUTOFF = 0
        MIN_DPHI = 1.0

        MIN_CSC_DNN = 0
        MIN_DT_DNN = 0

    if CUTSET == 'scs': # Standard cut selections for CSC-DT (no met categorization)
        MIN_CSC_TIME = -5.0
        MAX_CSC_TIME =  12.5
        MAX_CSC_TSPREAD =  20.0

        MAX_RPC_BX = 0
        MIN_RPC_HITS = 1

        MAX_CSC_JET = 30
        MAX_DT_JET = 50
        MAX_CSC_MUON = 30
        MAX_DT_MUON = 10

        MAX_ME1 = 0
        MAX_MB1 = 0

        HALO_CUTOFF = 0
        MIN_DPHI = 1.0

        MIN_CSC_DNN = 0
        MIN_DT_DNN = 0

    if CUTSET == 'ropt': # Standard cut selections for CSC-DT (no met categorization)
        MIN_CSC_TIME = -5.0
        MAX_CSC_TIME =  12.5
        MAX_CSC_TSPREAD =  20.0

        #MAX_RPC_BX = 0
        #MIN_RPC_HITS = 1

        MAX_CSC_JET = 20
        #MAX_DT_JET = 50
        MAX_CSC_MUON = 10
        #MAX_DT_MUON = 100

        MAX_ME1 = 0
        #MAX_MB1 = 0

        #HALO_CUTOFF = 0
        MIN_DPHI = 1.8

        MIN_CSC_DNN = 0.96137905
        #MIN_DT_DNN = 0
    
    # **************************** #
    rdfs = {
        'mc' : RDataFrame('MuonSystem', FN_MC),
        'r3' : RDataFrame('MuonSystem', FN_R3),
    }

    print('Events in:')
    for key, rdf in rdfs.items():
        # rdf = rdf.Range(0,100_000) # Skim a subset of events for debugging

        if key == 'mc': # fix weights
            rdf = rdf.Redefine('weight', f'weight * {LUMI}')

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
            if f'{p}R' not in COLUMNS_OUT:
                COLUMNS_OUT.append(f'{p}R')

        rdf = rdf.Define(f'{C}NHitME1', f'{C}NRechitChamberPlus11 + {C}NRechitChamberPlus12 + {C}NRechitChamberMinus11 + {C}NRechitChamberMinus12')
        if f'{C}NHitME1' not in COLUMNS_OUT:
            COLUMNS_OUT.append(f'{C}NHitME1')

        # **** #
        rdfs[key] = rdf

    print('')

    # **************** #
    # Apply selections from CUTS (and/or arguments?)
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
            print(f'    sample & {ec0:,.0f} & -- & -- & {cc00:,} & -- & -- & {cc10:,} & -- & -- & {dc0:,} & -- & -- \\\\')

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
            # Trigger selections (HLT may be wrong in MC?)
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
            # Missing Et requirements & categorization
            if 'MET' in cut:
                rdf = rdf.Redefine('evtCutFlag', f'evtCutFlag && (met < 200)')

            if 'low MET' in cut:
                rdf = rdf.Redefine('evtCutFlag', f'evtCutFlag && (met < {LOW_MET_CUTOFF})')

            if 'high MET' in cut:
                rdf = rdf.Redefine('evtCutFlag', f'evtCutFlag && ((met > {HIGH_MET_CUTOFF}) && (met < 200))')

            # **** #
            # Cluster time requirements
            if 'CSC0 IT' in cut:
                rdf = rdf.Redefine(f'{C}0CutFlag', f'{C}0CutFlag && ( ({MIN_CSC_TIME} < {C}TimeWeighted) && ({C}TimeWeighted < {MAX_CSC_TIME}) && ({C}TimeSpreadWeightedAll < {MAX_CSC_TSPREAD}) )')
            if 'CSC1 IT' in cut:
                rdf = rdf.Redefine(f'{C}1CutFlag', f'{C}1CutFlag && ( ({MIN_CSC_TIME} < {C}TimeWeighted) && ({C}TimeWeighted < {MAX_CSC_TIME}) && ({C}TimeSpreadWeightedAll < {MAX_CSC_TSPREAD}) )')
            
            if 'CSC0 OOT' in cut:
                rdf = rdf.Redefine(f'{C}0CutFlag', f'{C}0CutFlag && !( ({MIN_CSC_TIME} < {C}TimeWeighted) && ({C}TimeWeighted < {MAX_CSC_TIME}) && ({C}TimeSpreadWeightedAll < {MAX_CSC_TSPREAD}) )')
            if 'CSC1 OOT' in cut:
                rdf = rdf.Redefine(f'{C}1CutFlag', f'{C}1CutFlag && !( ({MIN_CSC_TIME} < {C}TimeWeighted) && ({C}TimeWeighted < {MAX_CSC_TIME}) && ({C}TimeSpreadWeightedAll < {MAX_CSC_TSPREAD}) )')

            if 'DT IT' in cut:
                rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( (abs({D}_match_RPCBx_dPhi0p5) <= {MAX_RPC_BX}) && ({D}_match_RPChits_dPhi0p5 >= {MIN_RPC_HITS}) )')
            elif 'DT OOT' in cut:
                rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && !( (abs({D}_match_RPCBx_dPhi0p5) <= {MAX_RPC_BX}) && ({D}_match_RPChits_dPhi0p5 >= {MIN_RPC_HITS}) )')

            # **** #
            # Station requirements from L1 trigger
            if 'ME1' in cut:
                # rdf = rdf.Redefine(f'{C}1CutFlag', f'{C}1CutFlag && ( '+
                #                 f'({C}NRechitChamberPlus11 <= {MAX_ME1}) &&'+
                #                 f'({C}NRechitChamberPlus12 <= {MAX_ME1}) &&'+
                #                 f'({C}NRechitChamberMinus11 <= {MAX_ME1}) &&'+
                #                 f'({C}NRechitChamberMinus12 <= {MAX_ME1}) )')
                    rdf = rdf.Redefine(f'{C}1CutFlag', f'{C}1CutFlag && ( {C}NHitME1 <= {MAX_ME1} ) ')

            if 'MB1' in cut:
                rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( {D}NHitStation1 <= {MAX_MB1} )')

            if 'DT stn' in cut:
                rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( ({D}NStation10 < 3) && !(({D}NStation10 == 2) && ({D}MaxStation == 4)) )')

            # **** #
            if 'jet veto' in cut:
                rdf = rdf.Redefine(f'{C}0CutFlag', f'{C}0CutFlag && ({C}JetVetoPt < {MAX_CSC_JET})')
                rdf = rdf.Redefine(f'{C}1CutFlag', f'{C}1CutFlag && ({C}JetVetoPt < {MAX_CSC_JET})')

            
            if 'muon veto' in cut:
                rdf = rdf.Redefine(f'{C}0CutFlag', f'{C}0CutFlag && (({C}MuonVetoGlobal == 0) || (({C}MuonVetoGlobal == 1) && ({C}MuonVetoPt < {MAX_CSC_MUON})))')
                rdf = rdf.Redefine(f'{C}1CutFlag', f'{C}1CutFlag && (({C}MuonVetoGlobal == 0) || (({C}MuonVetoGlobal == 1) && ({C}MuonVetoPt < {MAX_CSC_MUON})))')

            if 'halo veto' in cut:
                rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( ({HALO_CUTOFF} < abs({D}Phi)) && (abs({D}Phi) < {PI} - {HALO_CUTOFF}) )')

            # **** #
            if 'BDT' in cut:
                raise NotImplementedError('BDT')

            if 'DNN' in cut:
                rdf = rdf.Redefine(f'{C}0CutFlag', f'{C}0CutFlag && ( Take({C}DNN_bkgMC_plusBeamHalo,nCscRechitClusters) > {MIN_CSC_DNN} )')
                rdf = rdf.Redefine(f'{C}1CutFlag', f'{C}1CutFlag && ( Take({C}DNN_bkgMC_plusBeamHalo,nCscRechitClusters) > {MIN_CSC_DNN} )')
                #rdf = rdf.Redefine(f'{D}CutFlag', f'{D}CutFlag && ( Take({D}DNN,nDtRechitClusters) > {MIN_DT_DNN} )')

            # **** #
            if '1 CSC-CSC' in cut:
                # rdf = rdf.Redefine('evtCutFlag', f'evtCutFlag && evtFlag && '+
                #                    f'(reduce({C}0Flag.begin(), {C}0Flag.end()) == 1) && '+
                #                    f'(reduce({C}1Flag.begin(), {C}1Flag.end()) == 1)')
                rdf = rdf.Redefine(f'{C}0Flag', f'{C}0Flag && ( {C}Size == Max({C}Size*{C}0Flag) )')
                rdf = rdf.Redefine(f'{C}1Flag', f'{C}1Flag && ( {C}Size == Max({C}Size*{C}1Flag) )')

                rdf = rdf.Redefine('evtCutFlag', f'evtCutFlag && evtFlag && (reduce({D}Flag.begin(), {D}Flag.end()) == 0)')
                # rdf = rdf.Redefine(f'{D}Flag', f'{D}Flag && ( {D}Size == Max({D}Size*{D}Flag) )')

                # Apply our cluster level selections to relevant columns
                columns_out = []
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
                columns_out.extend(['tag_dEta', 'tag_dPhi', 'tag_dR', 'csc0CTau', 'csc1CTau'])

            # **** #
            if 'dR' in cut:
                raise NotImplementedError('dR')

            if 'dEta' in cut:
                raise NotImplementedError('dEta')

            if 'dPhi' in cut:
                rdf = rdf.Redefine('evtCutFlag', f'evtCutFlag && (tag_dPhi > {MIN_DPHI})')

            # **** #
            # Propagate to cumulative flags
            rdf = rdf.Redefine(f'{C}0Flag', f'{C}0Flag && {C}0CutFlag')
            rdf = rdf.Redefine(f'{C}1Flag', f'{C}1Flag && {C}1CutFlag')
            rdf = rdf.Redefine(f'{D}Flag', f'{D}Flag && {D}CutFlag')

            rdf = rdf.Redefine('evtFlag', 'evtFlag && evtCutFlag && ('+
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

    # **************** #
    print('Events out:')
    for key, rdf in rdfs.items():
        count, weight = rdf.Count(), rdf.Sum('weight')
        count, weight = count.GetValue(), weight.GetValue()
        if count == 0:
            print(f'{key} is empty')
            continue
        
        name = f'{key}_csccsc{"OOT" if OOT else ""}_{CUTSET}'
        name += f"_{MET_CATEGORY}"
        rdf = rdf.Snapshot('MuonSystem_flat', f'data/processed/{name}_rdf.root', columns_out)
        rdfs[key] = rdf

        print(f'  {key} = {count:,} ({weight:,.2f})')
        for xx in ('met', 'csc0Size', 'csc0R', 'csc0Eta', 'csc0Phi', 'csc1Size', 'csc1R', 'csc1Eta', 'csc1Phi', 'tag_dEta','tag_dPhi','tag_dR'):
            canvas = rt.TCanvas('','',800,800)
            xmin, xmax, std = rdf.Min(xx).GetValue(), rdf.Max(xx).GetValue(), rdf.StdDev(xx).GetValue()
            nbins = int((xmax-xmin)/(2.718*std*count**(-1/3))) if std else 1
            hh = rdf.Histo1D((f'{key}_{xx}',f'{key};{xx};count',nbins,xmin,xmax),f'{xx}').GetValue()
            hh.SetMinimum(0)
            hh.Draw()
            canvas.Draw()
            canvas.Print(f'{OUT_DIR}/{name}_{xx}.png')
    print('')
