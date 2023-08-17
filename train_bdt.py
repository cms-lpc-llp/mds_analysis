"""train_bdt.py
To run with all events :
'python train_bdt.py -1'
"""
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
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import roc_curve

############

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
    save_date = 'aug17'
    OUT_DIR = f'reports/weekly/{save_date}'

    LUMI = 23.02*1000
    N_EVENTS = 100_000
    CUT = True
    ROOT_BATCH = True

    gc = []
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
    print('|           Filtering MuonSystems           |')
    print('+-------+-----------------+-----------------+')
    print(f'|    In | MC : {ms_mc.count():>10,} | R3 : {ms_r3.count():>10,} |')

    ms_mc.match_mc('csc,dt')
    print(f'| Match | MC : {ms_mc.count():>10,} | R3 : {ms_r3.count():>10,} |')

    # ! I cannot figure out how to load HLTDecision without overflowing memory
    # ms_mc.cut_hlt()
    # ms_r3.cut_hlt()
    # print(f'|   HLT | MC : {ms_mc.count():>10,} | R3 : {ms_r3.count():>10,} |')

    # ms_mc.f(ms_mc['dtRechitClusterNHitStation1'] == 0, 'dt')
    # print(f'|   MB1 | MC : {ms_mc.count():>10,} | R3 : {ms_r3.count():>10,} |')

    ms_mc.cut_l1()
    ms_r3.cut_l1()
    print(f'|    L1 | MC : {ms_mc.count():>10,} | R3 : {ms_r3.count():>10,} |')

    ms_mc.cut_time('csc,dt', cut_csc_spread=True, cut_rpc_hits=True)
    ms_r3.cut_time('csc,dt', cut_csc_spread=True, cut_rpc_hits=True)
    print(f'|    IT | MC : {ms_mc.count():>10,} | R3 : {ms_r3.count():>10,} |')

    ms_mc.tag(tags='cscdt')
    ms_r3.tag(tags='cscdt')
    print(f'| CSCDT | MC : {ms_mc.count():>10,} | R3 : {ms_r3.count():>10,} |')

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

    print('\t+--------------+')
    print('\t| Making Plots |')
    print('\t+--------------+')


    xls = ['tag_dEta',
           'tag_dPhi',
           'tag_dR',
           'nCsc',
           'nDt',
           'met',
           'runNum']
    bins = [[32,0,3.5],
            [32,0,np.pi],
            [32,0,5],
            [6,-0.5,5.5],
            [6,-0.5,5.5],
            [30,0,500],
            [2741,360019, 362760]]


    for ix, xl in enumerate(xls):
        cn = 'c'+str(np.random.randint(999999999))
        c = TCanvas(cn,cn, 800, 800)
        leg = []
        lat = TLatex()
        lat.SetTextAlign(11)
        lat.SetTextSize(0.04)
        hmin, hmax = 999, -1

        for ims, ms in enumerate([ms_r3, ms_mc]):
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

        c.Print(f'{OUT_DIR}/{xl}.png')

    print('\t+------------------+')
    print('\t| Making CSC Plots |')
    print('\t+------------------+')

    xls = [
        'cscMuonVetoPt',
        'cscJetVetoPt',
        'cscPhi',
        'cscEta',
        'cscSize',
        'cscNStation10',
        'cscMaxChamber',
        'cscAvgStation10',
        'cscMaxStation',
        'cscMe11Ratio',
        'cscMe12Ratio'
    ]

    bins = [[50,0,100],
            [50,0,100],
            [25,-np.pi,np.pi],
            [50,0.8,2.5],
            [50,0,2000],
            [6,-0.5,5.5],
            [50,0,50],
            [50,0,5],
            [6,-0.5,5.5],
            [50,0,1],
            [50,0,1]]

    for ix, xl in enumerate(xls):
        cn = 'c'+str(np.random.randint(999999999))
        c = TCanvas(cn,cn, 800, 800)
        leg = []
        lat = TLatex()
        lat.SetTextAlign(11)
        lat.SetTextSize(0.04)
        hmin, hmax = 999, -1

        for ims, ms in enumerate([ms_r3, ms_mc]):
            k, var, weight = ms.name, ms[xl], ms['weight']
            if xl == 'cscEta':
                var = np.abs(var)

            # cond = np.abs(cscClusterSize[k][:,0])>1000
            # cond = np.abs(cscClusterEta[k][:,0])>1.6

            # cond = np.abs(cscClusterEta[k][:,0])<1.2
            #! Why did you do dPhi > -0.5?
            cond = (ms['tag_dPhi'] > 0.5) & (ms['dtSize'][:,0]<100)
            # cond = np.logical_and(np.abs(cscClusterEta[k][:,0])<1.6, np.abs(cscClusterEta[k][:,0])>1.2)
            # 0.9 1.2
            # 1.2 1.6
            # cond = np.logical_and(cond, np.abs(dPhi[k])>1.5)
            # cond = np.logical_and(cond, np.abs(cscClusterSize[k][:,0])>300)
            if not ms.is_mc:
                var, weight = var[cond], weight[cond]

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

        c.Print(f'{OUT_DIR}/{xl}.png')

    print('\t+-----------------+')
    print('\t| Making DT Plots |')
    print('\t+-----------------+')

    xls = [
        'dtJetVetoPt',
        'dtMuonVetoPt',
        'dtPhi',
        'dtEta',
        'dtSize',
        'dt_match_RPCBx_dPhi0p5',
        'dtNStation10',
        'dtMaxStation',
        'dtNHitStation1',
        'dtMb1Ratio',
    ]

    bins = [[50,0,100],
            [50,0,100],
            [25,-np.pi,np.pi],
            [50,0,1.3],
            [50,0,500],
            [20,-10,10],
            [6,-0.5,5.5],
            [6,-0.5,5.5],
            [50,0,50],
            [50,0,1]]

    for ix, xl in enumerate(xls):
        cn = 'c'+str(np.random.randint(999999999))
        c = TCanvas(cn,cn, 800, 800)
        leg = []
        lat = TLatex()
        lat.SetTextAlign(11)
        lat.SetTextSize(0.04)
        hmin, hmax = 999, -1

        for ims, ms in enumerate([ms_r3, ms_mc]):
            k, var, weight = ms.name, ms[xl], ms['weight']
            if xl == 'cscEta':
                var = np.abs(var)

            # if k == 'data':sel_dtcluster =  dtClusterSize[k]<80
            # else: sel_dtcluster = dtClusterSize[k]>=50
            # sel_dtcluster = np.logical_and(sel_dtcluster, dtClusterTime[k] == 0)
            # sel_dtcluster = np.logical_and(sel_dtcluster, dtClusterJetVetoPt[k] < 10)
            # sel_dtcluster = np.logical_and(sel_dtcluster, dtClusterMuonVetoPt[k] < 10)

            cond = (ms['tag_dPhi'] > 0.5) & (ms['cscSize'][:,0]<300)
            if not ms.is_mc:
                var, weight = var[cond], weight[cond]

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

        c.Print(f'{OUT_DIR}/{xl}.png')


    ######################################################
    ######################################################
    ######################################################
    ######################################################

    print('+--------------+')
    print('| Training BDT |')
    print('+--------------+')

    # Convert to Numpy array
    feats_csc = [
        'cscPhi',
        'cscEta',
        'cscNStation10',
        'cscAvgStation10',
        'cscMaxStation',
        # 'cscMe11Ratio',
        # 'cscMe12Ratio',
        # 'cscJetVetoPt',
        # 'cscMuonVetoPt',
    ]
    feats_dt = [
        'dtPhi',
        'dtEta',
        'dtNStation10',
        'dtAvgStation10',
        'dtMaxStation',
        # 'dtMb1Ratio',
        # 'dtJetVetoPt',
        # 'dtMuonVetoPt',
    ]
    # feats_2tag = [
    #     'tag_dR',
    #     'tag_dEta',
    #     'tag_dPhi',
    # ]

    mc_blinded_idxs = ms_mc['tag_dPhi'] > 0.5
    r3_blinded_idxs = ms_r3['tag_dPhi'] > 0.5

    X_csc = np.array([np.r_[
        np.ravel(ms_mc[feat][mc_blinded_idxs]),
        np.ravel(ms_r3[feat][r3_blinded_idxs])
        ] for feat in feats_csc]).T
    X_dt = np.array([np.r_[
        np.ravel(ms_mc[feat][mc_blinded_idxs]),
        np.ravel(ms_r3[feat][r3_blinded_idxs])
        ] for feat in feats_dt]).T
    y_csc = np.r_[
        np.ones(len(ms_mc['sel_evt'][mc_blinded_idxs]), dtype=bool), 
        np.zeros(len(ms_r3['sel_evt'][r3_blinded_idxs]), dtype=bool)
        ]
    y_dt = np.r_[
        np.ones(len(ms_mc['sel_evt'][mc_blinded_idxs]), dtype=bool),
        np.zeros(len(ms_r3['sel_evt'][r3_blinded_idxs]), dtype=bool)
        ]
    print(X_csc.shape, X_dt.shape, y_csc.shape, y_dt.shape, np.sum(y_csc), np.sum(y_dt))

    sclr_csc, sclr_dt = StandardScaler(), StandardScaler()
    gbt_csc, gbt_dt = GradientBoostingClassifier(), GradientBoostingClassifier()

    X_trn_csc, X_tst_csc, y_trn_csc, y_tst_csc = train_test_split(X_csc, y_csc, test_size=0.1)
    X_trn_dt, X_tst_dt, y_trn_dt, y_tst_dt = train_test_split(X_dt, y_dt, test_size=0.1)

    X_trn_csc = sclr_csc.fit_transform(X_trn_csc)
    X_tst_csc = sclr_csc.transform(X_tst_csc)

    X_trn_dt = sclr_dt.fit_transform(X_trn_dt)
    X_tst_dt = sclr_dt.transform(X_tst_dt)

    gbt_csc.fit(X_trn_csc, y_trn_csc)
    gbt_dt.fit(X_trn_dt, y_trn_dt)

    pred_csc = gbt_csc.decision_function(X_csc)
    pred_dt = gbt_dt.decision_function(X_dt)

    print('\t+------------------+')
    print('\t| Making ROC Plots |')
    print('\t+------------------+')
    bkg = {}
    sig = {}
    bkg_eff = {}
    sig_eff = {}
    names = ['tag_dPhi', 'cscSize', 'dtSize', 'GBT_CSC', 'GBT_DT']

    for i, name in enumerate(names):
        bkg[name] = []
        sig[name] = []
        bkg_eff[name] = []
        sig_eff[name] = []
        denom_r3 = r3_blinded_idxs #ms_r3['tag_dPhi'] > 0.5
        denom_mc = mc_blinded_idxs #ms_mc['tag_dPhi'] > 0.5
        wmc, wr3 = ms_mc['weight'][denom_mc], ms_r3['weight'][denom_r3]
        if 'GBT_CSC' in name:
            vmc, vr3 = pred_csc[y_csc], pred_csc[~y_csc]
        elif 'GBT_DT' in name:
            vmc, vr3 = pred_dt[y_dt], pred_dt[~y_dt]
        else:
            vmc, vr3 = np.ravel(ms_mc[name][denom_mc]), np.ravel(ms_r3[name][denom_r3])


        if name == 'tag_dPhi':
            threshold = np.arange(0,3.2,0.01)
        elif 'GBT' in name:
            _min, _max = np.min(np.r_[vmc,vr3]),np.max(np.r_[vmc,vr3])
            threshold = np.arange(_min,_max,(_max-_min)/100)
        else:
            threshold = np.arange(0,1000,1)
        for th in threshold:
            cond = np.abs(vr3) > th
            bkg_eff[name].append(np.sum(wr3[cond])/np.sum(wr3))
            bkg[name].append(np.sum(wr3[cond]))

            cond = np.abs(vmc) > th
            sig_eff[name].append(np.sum(wmc[cond])/np.sum(wmc))
            sig[name].append(np.sum(wmc[cond]))

        sig_eff[name] = np.array(sig_eff[name])
        bkg_eff[name] = np.array(bkg_eff[name])

        sig_eff[name] = sig_eff[name][bkg_eff[name]>0]
        bkg_eff[name] = bkg_eff[name][bkg_eff[name]>0]

        idx = np.argmax(sig[name] / np.sqrt(bkg[name]))
        print(name, sig_eff[name][min(idx, len(sig_eff[name])-1)], bkg_eff[name][min(idx, len(bkg_eff[name])-1)])


    cn = 'c'+str(np.random.randint(999999999))
    c = TCanvas(cn,cn, 800, 800)
    leg = []
    lat = TLatex()
    lat.SetTextAlign(31)
    lat.SetTextSize(0.04)

    graph = {}
    for i, v in enumerate(names):
        graph[v] = create_TGraph(sig_eff[v], 1/bkg_eff[v], axis_title= ['signal efficiency', 'bkg rejection'])
        graph[v].SetLineWidth(5)
        graph[v].SetLineColor(std_color_list[i])
        leg.append([v, std_color_list[i]])
        graph[v].SetMaximum(50)
        graph[v].Draw('ac' if i == 0 else 'c same')

    for ileg, (text, color) in enumerate(leg):
        lat.SetTextColor(color)
        lat.DrawLatexNDC(0.94, 0.92 - ileg*0.04, text)

    # c.SetLogy()
    c.SetGrid()
    # c.SetRightMargin(0.04)

    if not ROOT_BATCH:
        c.Draw()

    c.Print(f'{OUT_DIR}/rocs.png')

    print(bkg_eff['GBT_DT'],     sig_eff['GBT_DT'], )


    print('+-----------------------+')
    print('| Finished train_bdt.py |')
    print('+-----------------------+')
