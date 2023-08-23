"""train_bdt.py

To run : 'python train_bdt.py {n_events}'
use n_events=-1 for all events
"""

__author__ = 'Paul Simmerling'
__email__ = 'psimmerl@caltech.edu'
__credits__ = [
    'Paul Simmerling', 
    'Christina Wang', 
    'Lisa Benato', 
    'Si Xie',
    'Cristian Pena',
    'Martin Kwok',
    'Pedro Fernandez Manteca',
    'Maria Spiropulu',
]

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

from ROOT import TCanvas, TLegend, TH1F, TH2F, TLatex, TBox, TLine
from ROOT import kRed, kBlue, kGreen, kCyan, kMagenta, kYellow, kBlack, kAzure
from ROOT.VecOps import RVec

from src.histo_utilities import create_TH1D, create_TGraph, std_color_list

from src.muon_system import MuonSystemRDF, MuonSystemAwkward# , multi_plot

# from src.muon_system import (get_lat_leg, draw_csc_z_boxes, draw_dt_r_boxes)
from src.helper_functions import canvas, alert


from src import CMS_lumi, tdrstyle

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.metrics import roc_curve, roc_auc_score

#############################
# Parameters -- make a yaml #
#############################
save_dstat = 'ca_0p6'
save_date = 'aug24'
OUT_DIR = f'reports/weekly/{save_date}'

LUMI = 23.02*1000
N_EVENTS = -1
CUT = True
ROOT_BATCH = True

gc = []
#############################

def hist1D(x_labels, bins, systems, out_dir=OUT_DIR):
    for ix, xl in enumerate(x_labels):
        cn = 'c'+str(np.random.randint(999999999))
        c = TCanvas(cn,cn, 800, 800)
        leg = []
        lat = TLatex()
        lat.SetTextAlign(11)
        lat.SetTextSize(0.04)
        hmin, hmax = 999, -1

        for ims, ms in enumerate(systems):
            k, var, weight = ms.name, ms[xl], ms['weight']

            h = create_TH1D(
                var,
                axis_title=[xl, 'fraction of events'],
                binning=bins[ix],
                weights=weight
            )
            h.SetLineColor(std_color_list[ims])

            # print(k, h.Integral(), np.sum(weight[np.abs(var)>1.5])/np.sum(weight))
            leg.append((k, std_color_list[ims]))

            if h.Integral()>0:
                h.Scale(1./h.Integral())

            hmax = max(hmax, h.GetMaximum())
            hmin = min(hmin, h.GetMinimum(0))

            # h.SetMinimum(1e-3)
            # h.SetMaximum(4e-1)
            h.Draw('hist same')
            gc.append(h)

        for ileg, (text, color) in enumerate(leg):
            lat.SetTextColor(color)
            lat.DrawLatexNDC(0.80-0.1*ileg,0.92,text)

        gc[-1].SetMaximum(hmax)
        gc[-2].SetMaximum(hmax)
        gc[-1].SetMinimum(hmin)
        gc[-2].SetMinimum(hmin)

        c.SetLogy()
        c.SetGrid()
        c.SetRightMargin(0.04)

        if not ROOT_BATCH:
            c.Draw()

        c.Print(f'{out_dir}/{xl}.png')


if __name__ == '__main__':
    alert('Starting train_bdt.py', '=', 'c')

    if ROOT_BATCH:
        SetBatch(ROOT_BATCH)

    #############################

    tdrstyle.setTDRStyle()
    CMS_lumi.writeExtraText = 0

    #############################

    if 'caltech' in os.uname()[1]:
        ff_mc = '/storage/cms/store/user/christiw/displacedJetMuonAnalyzer/Run3/V1p19/MC_Summer22EE/v1/sixie/v6/normalized/ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted.root'
        ff_r3 = '/storage/cms/store/user/christiw/displacedJetMuonAnalyzer/Run3/V1p19/Data2022/v6/normalized/DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi.root'
        OUT_DIR = '/storage/af/user/psimmerl/LLP/mdc_analysis/' + OUT_DIR
    else:
        ff_mc = 'data/raw/ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v6.root'
        ff_r3 = 'data/raw/DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi_v6.root'
        OUT_DIR = '/home/psimmerl/LLP/mdc_analysis/' + OUT_DIR

    #############################

    pathlib.Path(OUT_DIR).mkdir(parents=True, exist_ok=True)
    alert(f'Using output directory "{OUT_DIR}"')

    #############################


    if len(sys.argv) > 1:
        N_EVENTS = int(sys.argv[1])

    if N_EVENTS > 0:
        alert(f'Setting {N_EVENTS=:,}')
    else:
        N_EVENTS = None
        alert(f'Setting {N_EVENTS=} (Loading ALL events)')

    #############################

    ms_mc = MuonSystemAwkward(ff_mc, name='signal', nev=N_EVENTS, is_mc=True)
    ms_r3 = MuonSystemAwkward(ff_r3, name='data', nev=N_EVENTS, is_mc=False)

    #!!! TURNING CUTS OFF !!!!#
    ms_mc.cut, ms_r3.cut = False, False
    #!!!!!!!!!!!!!!!!!!!!!!!!!#

    print('')
    print('+-------------------------------------------+')
    print('|    Filtering MuonSystems (CSC>0 & DT>0)   |')
    print('+-------+-----------------+-----------------+')
    print(f'|    In | MC : {ms_mc.count():>10,} | R3 : {ms_r3.count():>10,} |')

    ms_mc.match_mc('csc,dt')
    print(f'| Match | MC : {ms_mc.count():>10,} | R3 : {ms_r3.count():>10,} |')

    # ! I cannot figure out how to load HLTDecision without overflowing memory
    # ms_mc.cut_hlt()
    # ms_r3.cut_hlt()
    # print(f'|   HLT | MC : {ms_mc.count():>10,} | R3 : {ms_r3.count():>10,} |')

    ms_mc.f(ms_mc['cscRechitClusterMe11Ratio'] + ms_mc['cscRechitClusterMe12Ratio'] == 0, 'csc')
    ms_r3.f(ms_r3['cscRechitClusterMe11Ratio'] + ms_r3['cscRechitClusterMe12Ratio'] == 0, 'csc')
    print(f'|   ME1 | MC : {ms_mc.count():>10,} | R3 : {ms_r3.count():>10,} |')

    ms_mc.f(ms_mc['dtRechitClusterNHitStation1'] == 0, 'dt')
    ms_r3.f(ms_r3['dtRechitClusterNHitStation1'] == 0, 'dt')
    print(f'|   MB1 | MC : {ms_mc.count():>10,} | R3 : {ms_r3.count():>10,} |')

    ms_mc.cut_l1()
    ms_r3.cut_l1()
    print(f'|    L1 | MC : {ms_mc.count():>10,} | R3 : {ms_r3.count():>10,} |')

    # ms_mc.cut_time('csc,dt', cut_csc_spread=True, cut_rpc_hits=True)
    # ms_r3.cut_time('csc,dt', cut_csc_spread=True, cut_rpc_hits=True)
    # print(f'|    IT | MC : {ms_mc.count():>10,} | R3 : {ms_r3.count():>10,} |')

    ms_mc.cut_time('csc', cut_csc_spread=True, cut_rpc_hits=True, invert=False)
    ms_r3.cut_time('csc', cut_csc_spread=True, cut_rpc_hits=True, invert=False)
    print(f'| CSCIT | MC : {ms_mc.count():>10,} | R3 : {ms_r3.count():>10,} |')

    ms_mc.cut_time('dt', cut_csc_spread=True, cut_rpc_hits=True, invert=False)
    ms_r3.cut_time('dt', cut_csc_spread=True, cut_rpc_hits=True, invert=True)
    print(f'| DTIOT | MC : {ms_mc.count():>10,} | R3 : {ms_r3.count():>10,} |')

    ms_mc.tag(tags='cscdt')
    ms_r3.tag(tags='cscdt')
    print(f'|  2tag | MC : {ms_mc.count():>10,} | R3 : {ms_r3.count():>10,} |')

    # ms_r3.blind('dphi')
    # print(f'|  dPhi | MC : {ms_mc.count():>10,} | R3 : {ms_r3.count():>10,} |')
    print('+-------+-----------------+-----------------+')
    print('')

    # Fix weights
    ms_mc['weight'] = ms_mc['weight'] * LUMI
    ms_r3['weight'] = ms_r3['weight'] * 0 + 1

    #!!! TURNING CUTS ON !!!!#
    ms_mc.cut, ms_r3.cut = True, True
    #!!!!!!!!!!!!!!!!!!!!!!!!#

    ######################################################
    ######################################################
    ######################################################

    alert('SKIPPING PLOTS', form='+', c='r')
    # alert('Making CSC Plots', form='-', c='g')

    # xls = [
    #     'cscMuonVetoPt',
    #     'cscJetVetoPt',
    #     'cscPhi',
    #     'cscEta',
    #     'cscSize',
    #     'cscNStation10',
    #     'cscMaxChamber',
    #     'cscAvgStation10',
    #     'cscMaxStation',
    #     'cscMe11Ratio',
    #     'cscMe12Ratio'
    # ]
    # bins = [[50,0,100],
    #         [50,0,100],
    #         [25,-np.pi,np.pi],
    #         [50,0.8,2.5],
    #         [50,0,2000],
    #         [6,-0.5,5.5],
    #         [50,0,50],
    #         [50,0,5],
    #         [6,-0.5,5.5],
    #         [50,0,1],
    #         [50,0,1]]

    # for ix, xl in enumerate(xls):
    #     cn = 'c'+str(np.random.randint(999999999))
    #     c = TCanvas(cn,cn, 800, 800)
    #     leg = []
    #     lat = TLatex()
    #     lat.SetTextAlign(11)
    #     lat.SetTextSize(0.04)
    #     hmin, hmax = 999, -1

    #     for ims, ms in enumerate([ms_r3, ms_mc]):
    #         k, var, weight = ms.name, ms[xl], ms['weight']
    #         if xl == 'cscEta':
    #             var = np.abs(var)

    #         # cond = np.abs(cscClusterSize[k][:,0])>1000
    #         # cond = np.abs(cscClusterEta[k][:,0])>1.6

    #         # cond = np.abs(cscClusterEta[k][:,0])<1.2
    #         #! Why did you do dPhi > -0.5?
    #         cond = (ms['tag_dPhi'] > 0.5) & (ms['dtSize'][:,0]<100)
    #         # cond = np.logical_and(np.abs(cscClusterEta[k][:,0])<1.6, np.abs(cscClusterEta[k][:,0])>1.2)
    #         # 0.9 1.2
    #         # 1.2 1.6
    #         # cond = np.logical_and(cond, np.abs(dPhi[k])>1.5)
    #         # cond = np.logical_and(cond, np.abs(cscClusterSize[k][:,0])>300)
    #         if not ms.is_mc:
    #             var, weight = var[cond], weight[cond]

    #         h = create_TH1D(
    #             var,
    #             axis_title=[xl, 'fraction of events'],
    #             binning=bins[ix],
    #             weights=weight
    #         )
    #         h.SetLineColor(std_color_list[ims])

    #         # print(k, h.Integral(), np.sum(weight[np.abs(var)>1.5])/np.sum(weight))
    #         leg.append((k, std_color_list[ims]))

    #         if h.Integral()>0:
    #             h.Scale(1./h.Integral())

    #         hmax = max(hmax, h.GetMaximum())
    #         hmin = min(hmin, h.GetMinimum(0))

    #         # h.SetMinimum(1e-3)
    #         # h.SetMaximum(4e-1)
    #         h.Draw('hist same')
    #         gc.append(h)

    #     for ileg, (text, color) in enumerate(leg):
    #         lat.SetTextColor(color)
    #         lat.DrawLatexNDC(0.80-0.1*ileg,0.92,text)

    #     gc[-1].SetMaximum(hmax)
    #     gc[-2].SetMaximum(hmax)
    #     gc[-1].SetMinimum(hmin)
    #     gc[-2].SetMinimum(hmin)

    #     c.SetLogy()
    #     c.SetGrid()
    #     c.SetRightMargin(0.04)

    #     if not ROOT_BATCH:
    #         c.Draw()

    #     c.Print(f'{OUT_DIR}/{xl}.png')

    # alert('Making DT Plots', form='-', c='g')

    # xls = [
    #     'dtJetVetoPt',
    #     'dtMuonVetoPt',
    #     'dtPhi',
    #     'dtEta',
    #     'dtSize',
    #     'dt_match_RPCBx_dPhi0p5',
    #     'dtNStation10',
    #     'dtMaxStation',
    #     'dtNHitStation1',
    #     'dtMb1Ratio',
    # ]
    # bins = [[50,0,100],
    #         [50,0,100],
    #         [25,-np.pi,np.pi],
    #         [50,0,1.3],
    #         [50,0,500],
    #         [20,-10,10],
    #         [6,-0.5,5.5],
    #         [6,-0.5,5.5],
    #         [50,0,50],
    #         [50,0,1]]

    # for ix, xl in enumerate(xls):
    #     cn = 'c'+str(np.random.randint(999999999))
    #     c = TCanvas(cn,cn, 800, 800)
    #     leg = []
    #     lat = TLatex()
    #     lat.SetTextAlign(11)
    #     lat.SetTextSize(0.04)
    #     hmin, hmax = 999, -1

    #     for ims, ms in enumerate([ms_r3, ms_mc]):
    #         k, var, weight = ms.name, ms[xl], ms['weight']
    #         if xl == 'cscEta':
    #             var = np.abs(var)

    #         # if k == 'data':sel_dtcluster =  dtClusterSize[k]<80
    #         # else: sel_dtcluster = dtClusterSize[k]>=50
    #         # sel_dtcluster = np.logical_and(sel_dtcluster, dtClusterTime[k] == 0)
    #         # sel_dtcluster = np.logical_and(sel_dtcluster, dtClusterJetVetoPt[k] < 10)
    #         # sel_dtcluster = np.logical_and(sel_dtcluster, dtClusterMuonVetoPt[k] < 10)

    #         cond = (ms['tag_dPhi'] > 0.5) & (ms['cscSize'][:,0]<300)
    #         if not ms.is_mc:
    #             var, weight = var[cond], weight[cond]

    #         h = create_TH1D(
    #             var,
    #             axis_title=[xl, 'fraction of events'],
    #             binning=bins[ix],
    #             weights=weight
    #         )
    #         h.SetLineColor(std_color_list[ims])

    #         # print(k, h.Integral(), np.sum(weight[np.abs(var)>1.5])/np.sum(weight))
    #         leg.append((k, std_color_list[ims]))

    #         if h.Integral()>0:
    #             h.Scale(1./h.Integral())

    #         hmax = max(hmax, h.GetMaximum())
    #         hmin = min(hmin, h.GetMinimum(0))

    #         # h.SetMinimum(1e-3)
    #         # h.SetMaximum(4e-1)
    #         h.Draw('hist same')
    #         gc.append(h)

    #     for ileg, (text, color) in enumerate(leg):
    #         lat.SetTextColor(color)
    #         lat.DrawLatexNDC(0.80-0.1*ileg,0.92,text)

    #     gc[-1].SetMaximum(hmax)
    #     gc[-2].SetMaximum(hmax)
    #     gc[-1].SetMinimum(hmin)
    #     gc[-2].SetMinimum(hmin)

    #     c.SetLogy()
    #     c.SetGrid()
    #     c.SetRightMargin(0.04)

    #     if not ROOT_BATCH:
    #         c.Draw()

    #     c.Print(f'{OUT_DIR}/{xl}.png')

    # alert('Making Event/2tag Plots', form='-', c='g')

    # xls = ['tag_dEta',
    #        'tag_dPhi',
    #        'tag_dR',
    #        'nCsc',
    #        'nDt',
    #        'met',
    #        'runNum']
    # bins = [[32,0,3.5],
    #         [32,0,np.pi],
    #         [32,0,5],
    #         [6,-0.5,5.5],
    #         [6,-0.5,5.5],
    #         [30,0,500],
    #         [2741,360019, 362760]]

    # for ix, xl in enumerate(xls):
    #     cn = 'c'+str(np.random.randint(999999999))
    #     c = TCanvas(cn,cn, 800, 800)
    #     leg = []
    #     lat = TLatex()
    #     lat.SetTextAlign(11)
    #     lat.SetTextSize(0.04)
    #     hmin, hmax = 999, -1

    #     for ims, ms in enumerate([ms_r3, ms_mc]):
    #         k, var, weight = ms.name, ms[xl], ms['weight']

    #         h = create_TH1D(
    #             var,
    #             axis_title=[xl, 'fraction of events'],
    #             binning=bins[ix],
    #             weights=weight
    #         )
    #         h.SetLineColor(std_color_list[ims])

    #         # print(k, h.Integral(), np.sum(weight[np.abs(var)>1.5])/np.sum(weight))
    #         leg.append((k, std_color_list[ims]))

    #         if h.Integral()>0:
    #             h.Scale(1./h.Integral())

    #         hmax = max(hmax, h.GetMaximum())
    #         hmin = min(hmin, h.GetMinimum(0))

    #         # h.SetMinimum(1e-3)
    #         # h.SetMaximum(4e-1)
    #         h.Draw('hist same')
    #         gc.append(h)

    #     for ileg, (text, color) in enumerate(leg):
    #         lat.SetTextColor(color)
    #         lat.DrawLatexNDC(0.80-0.1*ileg,0.92,text)

    #     gc[-1].SetMaximum(hmax)
    #     gc[-2].SetMaximum(hmax)
    #     gc[-1].SetMinimum(hmin)
    #     gc[-2].SetMinimum(hmin)

    #     c.SetLogy()
    #     c.SetGrid()
    #     c.SetRightMargin(0.04)

    #     if not ROOT_BATCH:
    #         c.Draw()

    #     c.Print(f'{OUT_DIR}/{xl}.png')

    ######################################################
    ######################################################
    ######################################################

    alert('Training BDT', form='-', c='g')

    # Convert to Numpy array
    feats_csc = [
        # 'cscSize',
        'cscPhi',
        'cscEta',
        # 'cscX',
        # 'cscY',
        'cscZ',
        'cscR',
        'cscNStation10',
        'cscAvgStation10',
        # 'cscMaxStation',
        # 'cscMe11Ratio',
        # 'cscMe12Ratio',
        'cscJetVetoPt',
        'cscMuonVetoPt',
    ]
    feats_dt = [
        'tag_dR',
        'dtSize',
        'dtPhi',
        'dtEta',
        # 'dtX',
        # 'dtY',
        'dtZ',
        'dtR',
        'dtNStation10',
        'dtAvgStation10',
        # 'dtMaxStation',
        # 'dtMb1Ratio',
        'dtJetVetoPt',
        'dtMuonVetoPt',
    ]
    feats_2tag = [
        'tag_dR',
        'tag_dEta',
        'tag_dPhi',
    ]

    #! Disbaled rn
    mc_blinded_idxs = ms_mc['tag_dPhi'] > -1# 0.5
    r3_blinded_idxs = ms_r3['tag_dPhi'] > -1#< 0.5

    X_csc = np.array([np.r_[
        np.ravel(ms_mc[feat][mc_blinded_idxs]),
        np.ravel(ms_r3[feat][r3_blinded_idxs])
        ] for feat in feats_csc]).T
    X_dt = np.array([np.r_[
        np.ravel(ms_mc[feat][mc_blinded_idxs]),
        np.ravel(ms_r3[feat][r3_blinded_idxs])
        ] for feat in feats_dt]).T
    
    y = np.r_[
        np.ones(len(ms_mc['sel_evt'][mc_blinded_idxs]), dtype=bool),
        np.zeros(len(ms_r3['sel_evt'][r3_blinded_idxs]), dtype=bool)
        ]
    w = np.r_[
        np.ravel(ms_mc['weight'][mc_blinded_idxs]),
        np.ravel(ms_r3['weight'][r3_blinded_idxs])
        ]

    X_csc, X_dt = np.abs(X_csc), np.abs(X_dt)

    print(X_csc.shape, X_dt.shape, y.shape, np.sum(y)/y.shape[0])

    # n_estimators = np.sqrt(min(len(feats_csc), len(feats_dt))*min(np.sum(y_csc), len(y_csc)-np.sum(y_csc)))
    n_estimators = np.sqrt(min(len(feats_csc), len(feats_dt))*len(y))
    # n_estimators = np.sqrt(min(len(feats_csc), len(feats_dt))*min(np.sum(y_csc), len(y_csc)-np.sum(y_csc)))
    n_estimators = 500#int(n_estimators)
    print(f'{n_estimators=}')

    # clf_csc = RandomForestClassifier(n_estimators=n_estimators, random_state=42, n_jobs=2)
    # clf_dt = RandomForestClassifier(n_estimators=n_estimators, random_state=42, n_jobs=2)
    clf_csc = GradientBoostingClassifier(n_estimators=n_estimators, random_state=42)#, max_depth=10)
    clf_dt = GradientBoostingClassifier(n_estimators=n_estimators, random_state=42)#, max_depth=10)

    X_trn_csc, X_tst_csc, X_trn_dt, X_tst_dt, y_trn, y_tst, w_trn, w_tst = train_test_split(X_csc, X_dt, y, w, test_size=2, random_state=42)

    sclr_csc, sclr_dt = StandardScaler(), StandardScaler()
    X_trn_csc = sclr_csc.fit_transform(X_trn_csc)
    X_tst_csc = sclr_csc.transform(X_tst_csc)

    X_trn_dt = sclr_dt.fit_transform(X_trn_dt)
    X_tst_dt = sclr_dt.transform(X_tst_dt)

    clf_csc.fit(X_trn_csc, y_trn, w_trn)
    clf_dt.fit(X_trn_dt, y_trn, w_trn)

    try:
        pred_csc = clf_csc.decision_function(X_csc)
        pred_dt = clf_dt.decision_function(X_dt)
    except:
        pred_csc = clf_csc.predict_proba(X_csc)[:,1]
        pred_dt = clf_dt.predict_proba(X_dt)[:,1]

    for fts, clf in [(feats_csc, clf_csc), (feats_dt, clf_dt)]:
        fts, wts = np.asarray(fts), clf.feature_importances_
        idxs = np.argsort(wts)[::-1]
        for ift, (ft, wt) in enumerate(zip(fts[idxs], wts[idxs])):
            print(f'| {ift:>2} - {ft:>16} | {wt:>6.4f} |')
        print('')

    ######################################################
    ######################################################
    ######################################################

    alert('Making ROC Plots', form='-', c='g')

    bkgs, sigs, aucs = {}, {}, {}
    bkg_effs, sig_effs = {}, {}
    names = ['tag_dPhi', 'cscSize', 'dtSize', 'cscNStation10', 'dtNStation10', 'GBT_CSC', 'GBT_DT']

    print('+-----------------------------------------------------------------------------------+')
    print('|                    Evaluating Signal/Background Discriminators                    |')
    print('+------------------+-------------+-------------+-----------+-----------+------------+')
    denom_r3 = r3_blinded_idxs #ms_r3['tag_dPhi'] > 0.5
    denom_mc = mc_blinded_idxs #ms_mc['tag_dPhi'] > 0.5
    wmc, wr3 = w[y], w[~y]
    for i, name in enumerate(names):
        bkgs[name], sigs[name], aucs[name] = [], [], 0
        bkg_effs[name], sig_effs[name] = [], []
        if 'GBT_CSC' in name:
            vmc, vr3 = pred_csc[y], pred_csc[~y]
        elif 'GBT_DT' in name:
            vmc, vr3 = pred_dt[y], pred_dt[~y]
        else:
            vmc, vr3 = ms_mc[name][denom_mc], ms_r3[name][denom_r3]

        vmc, vr3 = np.asarray(vmc), np.asarray(vr3)
        if len(vmc.shape) == 2:
            if vmc.shape[1] == 1:
                vmc, vr3 = vmc[:,0], vr3[:,0]
            else:
                print(f'too many columns, expected 1 saw {vmc.shape[1]}. (using col 0)')
                vmc, vr3 = vmc[:,0], vr3[:,0]

        #  if name == 'tag_dPhi':
        #     threshold = np.arange(0,np.pi,np.pi/50)
        # elif 'GBT' in name:
        #     _min, _max = np.min(np.r_[vmc,vr3]),np.max(np.r_[vmc,vr3])
        #     threshold = np.arange(_min,_max,(_max-_min)/100)
        # else:
        #     threshold = np.arange(0,1000,1)
        _min, _max = np.min(np.r_[vmc,vr3]),np.max(np.r_[vmc,vr3])
        threshold = np.arange(_min,_max,(_max-_min)/200)

        for th in threshold:
            cond = vr3 > th
            # cond = np.abs(vr3) > th
            bkg_effs[name].append(np.sum(wr3[cond])/np.sum(wr3))
            bkgs[name].append(np.sum(wr3[cond]))

            cond = vmc > th
            sig_effs[name].append(np.sum(wmc[cond])/np.sum(wmc))
            sigs[name].append(np.sum(wmc[cond]))

        sig_effs[name] = np.array(sig_effs[name])
        bkg_effs[name] = np.array(bkg_effs[name])
        sigs[name] = np.array(sigs[name])
        bkgs[name] = np.array(bkgs[name])
        aucs[name] = roc_auc_score(y_true=np.r_[np.ones_like(vmc), np.zeros_like(vr3)],
                                   y_score=np.r_[vmc, vr3],
                                   sample_weight=np.r_[wmc, wr3])

        sigs[name] = sigs[name][bkg_effs[name]>0]
        bkgs[name] = bkgs[name][bkg_effs[name]>0]
        sig_effs[name] = sig_effs[name][bkg_effs[name]>0]
        bkg_effs[name] = bkg_effs[name][bkg_effs[name]>0]

        sig, bkg, sig_eff, bkg_eff, auc = sigs[name], bkgs[name], sig_effs[name], bkg_effs[name], aucs[name]

        idx = np.argmax(sig / (np.sqrt(bkg)))# if bkg else 1))
        # if isinstance(idx, np.ndarray):
        #     alert(f'idx is a list ({idx=})', '#', 'r')
        #     idx = idx[0]
        # if idx >= len(sig_eff):
        #     alert(f'idx too long ({idx=}, {len(sig_eff)=})', '#', 'r')
        #     idx = len(sig_eff) - 1

        sig, bkg, sig_eff, bkg_eff = sig[idx], bkg[idx], sig_eff[idx], bkg_eff[idx]

        print(f'| {name:>16} | S: {sig:>8.2f} | B: {bkg:>8.2f} | SE: {sig_eff:>5.3f} | BE: {bkg_eff:>5.3f} | AUC: {auc:5.3f} |')
    print('+------------------+-------------+-------------+-----------+-----------+------------+')


    cn = 'c'+str(np.random.randint(999999999))
    c = TCanvas(cn,cn, 800, 800)
    leg = []
    lat = TLatex()
    lat.SetTextAlign(31)
    lat.SetTextSize(0.04)

    graph = {}
    for i, v in enumerate(names):
        graph[v] = create_TGraph(sig_effs[v], 1/bkg_effs[v], axis_title= ['signal efficiency', 'bkg rejection'])
        graph[v].SetLineWidth(5)
        graph[v].SetLineColor(std_color_list[i])
        leg.append([v, std_color_list[i]])
        # graph[v].SetMaximum(50)
        graph[v].Draw('ac' if i == 0 else 'c same')

    for ileg, (text, color) in enumerate(leg):
        lat.SetTextColor(color)
        # lat.DrawLatexNDC(0.94, 0.92 - ileg*0.04, text)
        lat.DrawLatexNDC(0.94, 0.92 - ileg*(0.92-0.08)/(len(leg)+1), text)

    # c.SetLogy()
    c.SetGrid()
    # c.SetRightMargin(0.04)

    if not ROOT_BATCH:
        c.Draw()

    c.Print(f'{OUT_DIR}/rocs.png')

    ######################################################
    ######################################################
    ######################################################

    alert('Finished train_bdt.py', '=', 'c')
