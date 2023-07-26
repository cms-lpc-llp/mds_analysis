"""Analyzes MuonSystem TTrees using RDataFrame using the MuonSystemRDF class
    - Generates (many...) histograms for further plotting + analysis
    - Saves data selections for faster processing"""

import os
import sys
import pathlib
from collections import defaultdict

import ROOT as rt
# from ROOT import gErrorIgnoreLevel
import numpy as np
import sklearn as skl
# import awkward as ak

from src import CMS_lumi, tdrstyle

from src.muon_system import MuonSystemRDF
from src.muon_system import (get_lat_leg, multi_plot, draw_csc_z_boxes, draw_dt_r_boxes)

# from src.muon_system import MuonSystem, MuonSystemRDF
# from src.helper_functions import lnot, land, lor, lxor, asum, ABS
# from src.muon_system import (get_lat_leg, ms.Histo1D(, ms.Histo2D(, multi_plot, make_cluster_eff_1D, draw_csc_z_boxes, draw_dt_r_boxes)

from src.helper_functions import canvas

pi = rt.TMath.Pi()
bins = {
    "met": (100, 0, 400),
    "metPhi": (100, -pi, pi),
    #
    "cscZ": (100, 400, 1100),
    "cscR": (100, 0, 800),
    "cscN": (11, -0.5, 10.5),
    "cscSize": (100, 0, 500),
    "cscTime": (100, -60, 60),
    "cscTimeSpread": (100, 0, 80),
    #'cscE' : (100, -100, 100),
    #'cscPt' : (100, -100, 100),
    "cscEta": (100, 0.5, 3),
    "cscPhi": (100, -pi, pi),
    "cscMetDPhi": (100, -pi, pi),
    "cscNSt": (9, -0.5, 8.5),
    "cscASt": (9, -0.5, 8.5),
    "cscMSt": (9, -0.5, 8.5),
    #
    "dtZ": (100, 0, 700),
    "dtR": (100, 180, 800),
    "dtN": (11, -0.5, 10.5),
    "dtSize": (100, 0, 500),
    "dtSizeSmall": (50, 40, 200),
    "dtTime": (7, -3.5, 3.5),
    #'dtTimeSpread' : (100, -100, 100),
    #'dtE' : (100, -100, 100),
    #'dtPt' : (100, -100, 100),
    "dtEta": (100, 0, 1.5),
    "dtPhi": (100, -pi, pi),
    "dtMetDPhi": (100, -pi, pi),
    "dtNSt": (9, -0.5, 8.5),
    "dtASt": (9, -0.5, 8.5),
    "dtMSt": (9, -0.5, 8.5),
    "mb1": (50, 0, 100),
    "mb1Ratio": (100, 0, 0.1),
    #
    "gllpN": (4, -0.5, 3.5),
    "gllpE": (100, 0, 200),
    "gllpPt": (100, 0, 200),
    "gllpEta": (100, -5, 5),
    "gllpPhi": (100, -pi, pi),
    #
    "lepN": (6, -0.5, 5.5),
    "lepE": (100, 0, 200),
    "lepPt": (100, 0, 150),
    "lepEta": (100, -3, 3),
    "lepPhi": (100, -pi, pi),
    #
    "jetN": (11, -0.5, 10.5),
    "jetE": (100, 0, 300),
    "jetPt": (100, 0, 200),
    "jetEta": (100, -5, 5),
    "jetPhi": (100, -pi, pi),
    #
    "dEta": (100, -3.5, 3.5),
    "dPhi": (100, -pi, pi),
    "dR": (100, 0, 5),
    #
    "AdEta": (100, 0, 3.5),
    "AdPhi": (100, 0, pi),
}

bins["cmsZ"] = (bins['cscZ'][0], min(bins['cscZ'][1], bins['dtZ'][1]), max(bins['cscZ'][2], bins['dtZ'][2]))
bins["cmsR"] = (bins['cscR'][0], min(bins['cscR'][1], bins['dtR'][1]), max(bins['cscR'][2], bins['dtR'][2]))

bins2D = {
    "met_cscMetDPhi": bins["met"] + bins["cscMetDPhi"],
    "met_dtMetDPhi": bins["met"] + bins["dtMetDPhi"],
    "dPhi_dEta": bins["dPhi"] + bins["dEta"],
    "dPhi_dR": bins["dPhi"] + bins["dR"],
    "dEta_dR": bins["dEta"] + bins["dR"],
    "dPhi_met": bins["dPhi"] + bins["met"],
    "dEta_met": bins["dEta"] + bins["met"],
    "dR_met": bins["dR"] + bins["met"],
    "AdPhi_dtSize": bins["AdPhi"] + bins["dtSize"],
    "AdPhi_dtSizeSmall": bins["AdPhi"] + bins["dtSizeSmall"],
    "AdPhi_AdEta": bins["AdPhi"] + bins["AdEta"],
    "AdPhi_dR": bins["AdPhi"] + bins["dR"],
    "AdEta_dR": bins["AdEta"] + bins["dR"],
    "AdPhi_met": bins["AdPhi"] + bins["met"],
    "AdEta_met": bins["AdEta"] + bins["met"],
    "cmsZ_cmsR": bins["cmsZ"] + bins["cmsR"],
    "jetPt_mb1": bins["jetPt"] + bins["mb1"],
    "muonPt_mb1": bins["lepPt"] + bins["mb1"],
}

bins = {**bins, **bins2D}


def build_hists(ms, hhs, comment=''):

    if comment and comment[0] != '_':
        comment = '_' + comment
    #################################
    ## System & Cluster Kinematics ##
    #################################
    # print("Generating system and cluster kinematics plots")

    # System level
    hhs["met" + comment].append(ms.Histo1D(("met", ";E_{T}^{miss} [GeV];count", *bins["met"]), "met"))
    hhs["metPhi" + comment].append(ms.Histo1D(("metPhi", ";#phi[E_{T}^{miss}];count", *bins["metPhi"]), "metPhi"))

    # kinematics: (csc, dt, gllp, jet, lep)
    # N, Size, Time, TimeSpread, NSt, ASt, MSt, MetDPhi, E, Pt, Eta, Phi
    hhs["cscN" + comment].append(
        ms.Histo1D(("nCscRechitClusters", ";N CSC Clusters;count", *bins["cscN"]), "nCscRechitClusters"))
    hhs["cscSize" + comment].append(
        ms.Histo1D(("cscRechitClusterSize", ";CSC Rechit Size;count", *bins["cscSize"]), "cscRechitClusterSize"))
    hhs["cscTime" + comment].append(
        ms.Histo1D(("cscRechitClusterTimeWeighted", ";CSC Time [ns];count", *bins["cscTime"]),
                   "cscRechitClusterTimeWeighted"))
    hhs["cscTimeSpread" + comment].append(
        ms.Histo1D(("cscTimeSpread", ";CSC Time Spread [ns];count", *bins["cscTimeSpread"]),
                   "cscRechitClusterTimeSpreadWeightedAll"))
    hhs["cscNSt" + comment].append(
        ms.Histo1D(("cscRechitClusterNStation10", ";N Stations / CSC Cluster;count", *bins["cscNSt"]),
                   "cscRechitClusterNStation10"))
    hhs["cscASt" + comment].append(
        ms.Histo1D(("cscASt", ";Avg Station / CSC Cluster;count", *bins["cscASt"]), "ABScscRechitClusterAvgStation10"))
    hhs["cscMSt" + comment].append(
        ms.Histo1D(("cscMSt", ";Max Station / CSC Cluster;count", *bins["cscMSt"]), "ABScscRechitClusterMaxStation"))
    hhs["cscMetDPhi" + comment].append(
        ms.Histo1D(("cscMetDPhi", ";#Delta#phi(CSC,E_{T}^{miss});count", *bins["cscMetDPhi"]),
                   "cscRechitClusterMet_dPhi"))
    hhs["cscEta" + comment].append(
        ms.Histo1D(("cscRechitClusterEta", ";CSC #eta;count", *bins["cscEta"]), "ABScscRechitClusterEta"))
    hhs["cscPhi" + comment].append(
        ms.Histo1D(("cscRechitClusterPhi", ";CSC #phi;count", *bins["cscPhi"]), "cscRechitClusterPhi"))

    hhs["dtN" + comment].append(
        ms.Histo1D(("nDtRechitClusters", ";N DT Clusters;count", *bins["dtN"]), "nDtRechitClusters"))
    hhs["dtSize" + comment].append(
        ms.Histo1D(("dtRechitClusterSize", ";DT Rechit Size;count", *bins["dtSize"]), "dtRechitClusterSize"))
    hhs["dtTime" + comment].append(
        ms.Histo1D(("dtRechitCluster_match_RPCBx_dPhi0p5", ";RPC_{DT,matched} B_{X};count", *bins["dtTime"]),
                   "dtRechitCluster_match_RPCBx_dPhi0p5"))
    hhs["dtNSt" + comment].append(
        ms.Histo1D(("dtRechitClusterNStation10", ";N Stations / DT Cluster;count", *bins["dtNSt"]),
                   "dtRechitClusterNStation10"))
    hhs["dtASt" + comment].append(
        ms.Histo1D(("dtRechitClusterAvgStation10", ";Avg Station / DT Cluster;count", *bins["dtASt"]),
                   "ABSdtRechitClusterAvgStation10"))
    hhs["dtMSt" + comment].append(
        ms.Histo1D(("dtRechitClusterMaxStation", ";Max Station / DT Cluster;count", *bins["dtMSt"]),
                   "ABSdtRechitClusterMaxStation"))
    hhs["dtMetDPhi" + comment].append(
        ms.Histo1D(("dtRechitClusterMet_dPhi", ";#Delta#phi(DT,E_{T}^{miss});count", *bins["dtMetDPhi"]),
                   "dtRechitClusterMet_dPhi"))
    hhs["dtEta" + comment].append(
        ms.Histo1D(("dtRechitClusterEta", ";DT #eta;count", *bins["dtEta"]), "ABSdtRechitClusterEta"))
    hhs["dtPhi" + comment].append(
        ms.Histo1D(("dtRechitClusterPhi", ";DT #phi;count", *bins["dtPhi"]), "dtRechitClusterPhi"))

    hhs["gllpN" + comment].append(ms.Histo1D(("nGLLP", ";N gLLPs;count", *bins["gllpN"]), "nGLLP"))
    hhs["gllpE" + comment].append(ms.Histo1D(("gLLP_e", ";gLLP E [GeV];count", *bins["gllpE"]), "gLLP_e"))
    hhs["gllpPt" + comment].append(ms.Histo1D(("gLLP_pt", ";gLLP P_{T} [GeV];count", *bins["gllpPt"]), "gLLP_pt"))
    hhs["gllpEta" + comment].append(ms.Histo1D(("gLLP_eta", ";gLLP #eta;count", *bins["gllpEta"]), "gLLP_eta"))
    hhs["gllpPhi" + comment].append(ms.Histo1D(("gLLP_phi", ";gLLP #phi;count", *bins["gllpPhi"]), "gLLP_phi"))

    hhs["jetN" + comment].append(ms.Histo1D(("nJets", ";N Jets;count", *bins["jetN"]), "nJets"))
    hhs["jetE" + comment].append(ms.Histo1D(("jetE", ";Jet E [GeV];count", *bins["jetE"]), "jetE"))
    hhs["jetPt" + comment].append(ms.Histo1D(("jetPt", ";Jet P_{T} [GeV];count", *bins["jetPt"]), "jetPt"))
    hhs["jetEta" + comment].append(ms.Histo1D(("jetEta", ";Jet #eta;count", *bins["jetEta"]), "jetEta"))
    hhs["jetPhi" + comment].append(ms.Histo1D(("jetPhi", ";Jet #phi;count", *bins["jetPhi"]), "jetPhi"))

    hhs["lepN" + comment].append(ms.Histo1D(("nLeptons", ";N Leptons;count", *bins["lepN"]), "nLeptons"))
    hhs["lepE" + comment].append(ms.Histo1D(("lepE", ";Lepton E [GeV];count", *bins["lepE"]), "lepE"))
    hhs["lepPt" + comment].append(ms.Histo1D(("lepPt", ";Lepton P_{T} [GeV];count", *bins["lepPt"]), "lepPt"))
    hhs["lepEta" + comment].append(ms.Histo1D(("lepEta", ";Lepton #eta;count", *bins["lepEta"]), "lepEta"))
    hhs["lepPhi" + comment].append(ms.Histo1D(("lepPhi", ";Lepton #phi;count", *bins["lepPhi"]), "lepPhi"))

    hhs["mb1" + comment].append(ms.Histo1D(("mb1", ";N MB1;count", *bins["mb1"]), "dtRechitClusterNHitStation1"))
    hhs["mb1Ratio" + comment].append(
        ms.Histo1D(("mb1Ratio", ";N MB1 / N DT;count", *bins["mb1Ratio"]), "dtRechitClusterMB1Ratio"))

    return hhs


def alert(msg: str):
    print('')
    print('!' * (len(msg) + 6))
    print('!! ' + msg + ' !!')
    print('!' * (len(msg) + 6))


if __name__ == "__main__":

    # Parameters #
    # TODO - make a control card (yaml)
    isBlind = True
    justMC = False
    isCut = True
    save_dstat = "ca_0p6"
    save_date = 'may25'
    nev = 100_000
    nThreads = 8
    ##############

    gc = []
    rt.gROOT.SetBatch()

    cur_dir = os.getcwd()
    if 'Documents' in cur_dir:
        rt.EnableImplicitMT()
    else:
        rt.EnableImplicitMT(nThreads)

    print(f"Running with thread pool size = {rt.GetThreadPoolSize():,}")

    a = tdrstyle.setTDRStyle()
    # CMS_lumi.writeExtraText = 0
    rt.gStyle.SetOptFit(0)  # 1011)
    ssOpt = rt.RDF.RSnapshotOptions()
    ssOpt.fLazy = True

    out_dir = cur_dir + f"/reports/weekly/{save_date}/"
    out_plots_dir = out_dir + "/plots/"
    in_data_dir = cur_dir + "/data/raw/"
    out_data_dir = cur_dir + "/data/processed/"
    ending = ".png"

    if len(sys.argv) > 1:
        nev = int(sys.argv[1])

    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)  # make out directory if it doesn't exist
    pathlib.Path(out_plots_dir).mkdir(parents=True, exist_ok=True)  # make plots directory if it doesn't exist
    pathlib.Path(out_data_dir).mkdir(parents=True, exist_ok=True)  # make out data directory if it doesn't exist
    print(f"Using output directory '{out_dir}'")

    fout_name = out_data_dir + f"histograms_{save_date}.root"
    fout_root = rt.TFile(fout_name, "recreate")

    mc_db_0p4 = in_data_dir + "ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v4.root"
    mc_ca_0p4 = in_data_dir + "ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v5.root"
    mc_ca_0p5 = in_data_dir + "ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v7.root"
    mc_ca_0p6 = in_data_dir + "ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v6.root"
    mc_ca_0p8 = in_data_dir + "ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v8.root"
    mc_ca_1p0 = in_data_dir + "ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v9.root"

    r3_db_0p4 = in_data_dir + "DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi_v4.root"
    r3_ca_0p4 = in_data_dir + "DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi_v5.root"
    r3_ca_0p5 = in_data_dir + "DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi_v7.root"
    r3_ca_0p6 = in_data_dir + "DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi_v6.root"
    r3_ca_0p8 = in_data_dir + "DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi_v8.root"
    r3_ca_1p0 = in_data_dir + "DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi_v9.root"

    files_mc = [
        # mc_db_0p4,
        # mc_ca_0p4,
        # mc_ca_0p5,
        mc_ca_0p6,
        # mc_ca_0p8,
        # mc_ca_1p0,
    ]
    files_r3 = [
        # r3_db_0p4,
        # r3_ca_0p4,
        # r3_ca_0p5,
        r3_ca_0p6,
        # r3_ca_0p8,
        # r3_ca_1p0,
    ]

    labels_mc = [
        # "Signal: DBSCAN 0.4",
        # "Signal: CA 0.4",
        # "Signal: CA 0.5",
        "Signal: CA 0.6",
        # "Signal: CA 0.8",
        # "Signal: CA 1.0",
    ]
    labels_r3 = [
        # "Data: DBSCAN 0.4",
        # "Data: CA 0.4",
        # "Data: CA 0.5",
        "Data: CA 0.6",
        # "Data: CA 0.8",
        # "Data: CA 1.0",
    ]

    stats_mc = [
        # "mc_db_0p4",
        # "mc_ca_0p4",
        # "mc_ca_0p5",
        "mc_ca_0p6",
        # "mc_ca_0p8",
        # "mc_ca_1p0",
    ]
    stats_r3 = [
        # "r3_db_0p4",
        # "r3_ca_0p4",
        # "r3_ca_0p5",
        "r3_ca_0p6",
        # "r3_ca_0p8",
        # "r3_ca_1p0",
    ]
    # muon_systems_mc = [MuonSystem(ff, isMC=True, nev=nev) for ff in files_mc]
    # muon_systems_r3 = [MuonSystem(ff, isMC=False, nev=nev) for ff in files_r3]

    # muon_systems = muon_systems_mc + muon_systems_r3
    files = files_mc + files_r3
    labels = labels_mc + labels_r3
    stats = stats_mc + stats_r3

    # ================================ #
    # NOTE: Using convention "x vs y"  #
    # ================================ #
    hhs = defaultdict(list)  # hist dictionary

    if justMC:
        alert("WARNING: RUNNING WITH JUST MC")

    for ff, label, stat in zip(files, labels, stats):
        isMC = "mc" in stat
        if justMC and not isMC:
            continue

        print('\n' + 2 * (80 * '=' + '\n'))
        oms = MuonSystemRDF(ff, isMC=isMC, nev=nev)
        print(f"Loaded {stat} ({oms.Count():,} events)")

        oms.Define("ABScscRechitClusterAvgStation10", "abs(cscRechitClusterAvgStation10)")
        oms.Define("ABScscRechitClusterMaxStation", "abs(cscRechitClusterMaxStation)")
        oms.Define("ABScscRechitClusterEta", "abs(cscRechitClusterEta)")
        oms.Define("ABSdtRechitClusterAvgStation10", "abs(dtRechitClusterAvgStation10)")
        oms.Define("ABSdtRechitClusterMaxStation", "abs(dtRechitClusterMaxStation)")
        oms.Define("ABSdtRechitClusterEta", "abs(dtRechitClusterEta)")
        oms.Define(
            "dtRechitClusterMB1Ratio",
            "auto _mb1Ratio = dtRechitCluster_match_MB1hits_0p5; for (int i=0; i<nDtRechitClusters; i++) { _mb1Ratio[i] = dtRechitClusterNHitStation1[i] / dtRechitClusterSize[i];} return _mb1Ratio; "
        )

        # Make copy of orignal muon system to work on
        #   - For efficient analysis I should define all my RDFs now in as few ops as possible...
        ms = oms.Filter('met >= 0', implicit=False)

        #########################
        ## MAKE CUT-FLOW TABLE ##
        #########################

        # alert("MAKING CUT-FLOW TABLE")
        # ms.print_cutflow_table(match=isMC)
        # print('repeating but with DT size < 80 (BLINDED BKG REGION)')
        # ms.Filter('dtRechitClusterSize < 80', 'dt', implicit=False).print_cutflow_table(match=isMC)

        # # ---------------------------------------- #
        # # ---------------------------------------- #
        # # ---------------------------------------- #

        # if isMC:
        #     template = (("dtAllTime", ";RPC_{DT,matched} B_{X};count", *bins["dtTime"]),
        #                 "dtRechitCluster_match_RPCBx_dPhi0p5")
        #     cdl = "(nCscRechitClusters == 1) && (nDtRechitClusters == 1) && (nLeptons == 0)"
        #     mb1 = "dtRechitClusterNSegStation1 == 0"  

        #     hRPC_all = ms.Histo1D(*template)
        #     hRPC_jmv = ms.jet_cut(implicit=False).muon_cut().Histo1D(*template)
        #     hRPC_mb1 = ms.Filter(mb1, 'dt', implicit=False).Histo1D(*template)
        #     hRPC_cut = ms.Filter(cdl, implicit=False).Histo1D(*template)
        #     hRPC_mat = ms.match_clusters('cscdt', implicit=False).Histo1D(*template)
        #     hRPC_mat_jmv = ms.match_clusters('cscdt', implicit=False).jet_cut().muon_cut().Histo1D(*template)
        #     hRPC_mat_mb1 = ms.match_clusters('cscdt', implicit=False).Filter(mb1, 'dt').Histo1D(*template)
        #     hRPC_mat_cut = ms.match_clusters('cscdt', implicit=False).Filter(cdl).Histo1D(*template)

        #     _hhs = [hRPC_all, hRPC_jmv, hRPC_mb1, hRPC_cut, hRPC_mat, hRPC_mat_mb1, hRPC_mat_cut]
        #     _lls = [
        #         'no cuts', 'Jet+Muon Veto', 'MB1 Veto', '1CSC+1DT+0Lep', 'signal', 'sig+JM Veto', 'sig+MB1 Veto',
        #         'sig+(1CSC+1DT+0Lep)'
        #     ]
        #     _lxy = (0.2, 0.7, 0.45, 0.9)

        #     c = canvas(1, 1)
        #     c.cd(1).SetLogy()
        #     _hll = multi_plot(_hhs, _lls, ymin=1, ymax='log', lw=3, legxy=_lxy)
        #     c.Print(out_dir + 'rpc_compare_' + stat + ending)

        #     c = canvas(1, 1)
        #     c.cd(1).SetLogy()
        #     _hll = multi_plot(_hhs, _lls, ymin=3e-4, ymax=1.3, lw=3, legxy=_lxy, norm=True)
        #     c.Print(out_dir + 'rpc_compare_norm_' + stat + ending)

        #     hRPC_mb1_cl = hRPC_mb1.Clone()

        #     c = canvas(1, 1)
        #     hRPC_jmv.Divide(hRPC_all.GetPtr())
        #     hRPC_mb1.Divide(hRPC_all.GetPtr())
        #     hRPC_cut.Divide(hRPC_all.GetPtr())
        #     hRPC_mat_jmv.Divide(hRPC_mat.GetPtr())
        #     hRPC_mat_mb1.Divide(hRPC_mat.GetPtr())
        #     hRPC_mat_cut.Divide(hRPC_mat.GetPtr())
        #     _hhs = [hRPC_jmv, hRPC_mb1, hRPC_cut, hRPC_mat_jmv, hRPC_mat_mb1, hRPC_mat_cut]
        #     for hh in _hhs:
        #         hh.GetYaxis().SetTitle('acceptance')
        #     _lls = [
        #         'Jet+Muon Veto / no cuts', 'MB1 Veto / no cuts', '1CSC+1DT+0Lep / no cuts', 'sig+JM Veto / sig',
        #         'sig+MB1 Veto / sig', 'sig+(1CSC+1DT+0Lep) / sig'
        #     ]
        #     _hll = multi_plot(_hhs, _lls, ymin=0, ymax=1, lw=3, legxy=_lxy)
        #     c.Print(out_dir + 'rpc_compare_eff_' + stat + ending)

        #     c = canvas(1, 1)
        #     hRPC_mb1_cl.Divide(hRPC_mat.GetPtr())
        #     hRPC_mb1_cl.GetYaxis().SetTitle('efficiency')
        #     _hhs = [hRPC_mb1_cl]
        #     _lls = ['MB1 Veto / sig']
        #     _hll = multi_plot(_hhs, _lls, ymin=0, ymax=1, lw=3, legxy=_lxy)
        #     c.Print(out_dir + 'rpc_mb1_div_sig_eff_' + stat + ending)

        # ##

        # ms_csc = ms.time_cut('oot', 'csc', implicit=False)
        # hCsc = ms_csc.Histo1D(("nCscRechitClusters", ";N CSC Clusters;count", *bins["cscN"]), "nCscRechitClusters")

        # ms_cscT = ms.time_cut('oot', 'cscT', implicit=False)
        # hCscT = ms_cscT.Histo1D(("nCscRechitClusters", ";N CSC Clusters;count", *bins["cscN"]), "nCscRechitClusters")

        # c = canvas(1, 1)
        # _hll = multi_plot([hCsc, hCscT], ['csc', 'cscT'], ymin=0)  #, norm=True)
        # c.Print(out_dir + 'csc_oot_test_' + stat + ending)

        # # ---------------------------------------- #
        # # ---------------------------------------- #
        # # ---------------------------------------- #

        # ===================================== #
        # Making time plots before cuts on time #
        # ===================================== #

        # hhs["cscAllTime"].append(
        #     ms.Histo1D(("cscAllTime", ";CSC Time [ns];count", *bins["cscTime"]), "cscRechitClusterTimeWeighted"))
        # hhs["cscAllTimeSpread"].append(
        #     ms.Histo1D(("cscAllTimeSpread", ";CSC Time Spread [ns];count", *bins["cscTimeSpread"]),
        #                "cscRechitClusterTimeSpreadWeightedAll"))
        # hhs["dtAllTime"].append(
        #     ms.Histo1D(("dtAllTime", ";RPC_{DT,matched} B_{X};count", *bins["dtTime"]),
        #                "dtRechitCluster_match_RPCBx_dPhi0p5"))

        ##############################
        ## CUT, BLIND, & MATCH DATA ##
        ##############################

        # Blind Types:
        #   - DT OOT (or not trigger CSC OOT)
        #   - CSC+DT IT, 1.0 < dPhi < 2.0
        #   - CSC+DT IT, dtSize < 80

        # Histo Cuts:
        #   - Raw
        #   - match
        #   - jet
        #   - muon
        #   - L1
        #   - jet+muon+L1
        #   - ME1
        #   - MB1
        #   - ME1+MB1
        #   - jet+muon+L1+ME1+MB1

        # More info:
        #   - 2 cluster selection after jet+muon+...+... ?

        if isMC:
            alert("MATCHING SIGNAL MC")
            print("    - Matching to gLLP")
            print("    - Requiring gLLP to decay in a MD")
            ms.match_clusters(system='cscdt', in_det=True)

        build_hists(ms, hhs, comment='raw')

        # ===================== #
        # Simple Cut Histograms #
        # ===================== #

        alert("MAKING SIMPLE CUT HISTS")
        print("    - 1CSC1DT")
        print("    - L1 Plateau")
        print("    - Jet Veto")
        print("    - Muon Veto")
        print("    - HLT")
        print("    - Jet + Muon Veto")
        print("    - L1 + Jet + Muon Veto")
        print("    - 1CSC1DT + L1 + Jet + Muon Veto")

        # build_hists(ms.L1_plateau(implicit=False), hhs, comment='l1')
        # build_hists(ms.jet_cut(implicit=False), hhs, comment='jet')
        # build_hists(ms.muon_cut(implicit=False), hhs, comment='muon')
        # build_hists(ms.jet_cut(implicit=False).muon_cut(), hhs, comment='jet_muon')
        # build_hists(ms.L1_plateau(implicit=False).jet_cut().muon_cut(), hhs, comment='l1_jet_muon')

        _ms = ms.L1_plateau(implicit=False).Filter('dtRechitClusterSize > 100', 'dt').define_2tag_kins_and_cut('cscdt')
        if not isMC:
            _ms.Filter('tag_dPhi < 2')
        build_hists(_ms, hhs, comment='formb1')
        # hhs["2D_jetPt_mb1"].append(
        #     _ms.Filter("dtRechitClusterJetVetoLooseId == 1", "dt",
        #                implicit=False).Histo2D(("2D_jetPt_mb1", ";Jet P_{T} [GeV];N MB1 Hits;", *bins["jetPt_mb1"]),
        #                                        "dtRechitClusterJetVetoPt", "dtRechitClusterNSegStation1"))
        # hhs["2D_muonPt_mb1"].append(
        #     _ms.Filter("dtRechitClusterMuonVetoLooseId == 1", "dt",
        #                implicit=False).Histo2D(("2D_muonPt_mb1", ";Muon P_{T} [GeV];N MB1 Hits;", *bins["muonPt_mb1"]),
        #                                        "dtRechitClusterMuonVetoPt", "dtRechitClusterNSegStation1"))

        _ms = ms.L1_plateau(implicit=False)
        if not isMC:
            _ms.time_cut('oot', 'dt')
        build_hists(_ms, hhs, comment='dtOOT')
        # if isCut:
        #     alert("APPLYING CUTS")
        #     print("    - Applying jet cut on CSC & DT")
        #     print("    - Applying muon cut on CSC & DT")
        #     print("    - Applying L1 plateau on event (uses csc)")

        #     ms.jet_cut()
        #     ms.muon_cut()
        #     # HLT_CscCluster_Loose or HLT_L1CSCCluster_DTCluster50
        #     ms.Filter("(HLTDecision[566] > 0) || (HLTDecision[569] > 0)")

        #     build_hists(ms, hhs, comment='jmL1')

        # if isBlind and not isMC:
        #     alert("BLINDING DATA")
        #     print("    - Requiring DT to be OOT")
        #     print("    - Requiring CSC (not trigger) to be OOT")

        #     ms.time_cut('oot', 'dt cscT')
        # build_hists(ms, hhs, comment='jmL1_dtOOT')

        # MET vs (csc, dt) MET dPhi
        # hhs["2D_met_cscMetDPhi"].append(
        #     ms.Histo2D(
        #         ("2D_met_cscMetDPhi", ";E_{T}^{miss} [GeV];#Delta#phi(csc, #phi_{T}^{miss});", *bins["met_cscMetDPhi"]),
        #         "met", "cscRechitClusterMet_dPhi"))
        # hhs["2D_met_dtMetDPhi"].append(
        #     ms.Histo2D(("2D_met_dtMetDPhi", ";E_{T}^{miss} [GeV];#Delta#phi(dt, #phi_{T}^{miss});count",
        #                 *bins["met_dtMetDPhi"]), "met", "dtRechitClusterMet_dPhi"))

        ###########################################
        ## Single Pair CSC-CSC/CSC-DT Kinematics ##
        ###########################################

        alert("Reducing to 2tag analysis")
        ms.Filter("(HLTDecision[566] > 0) || (HLTDecision[569] > 0)")
        ms.jet_cut()
        ms.muon_cut()
        ms.L1_plateau()
        if isBlind and not isMC:
            ms.Filter('dtRechitClusterSize < 80', 'dt')
        ms.define_2tag_kins_and_cut(system="csccsc,cscdt")

        for second in ('dt',):  #('csc', 'dt'):
            print(f"Generating CSC-{second.upper()} kinematics plots")
            _ms = ms.Filter(f"tag_type == {0 if second == 'csc' else 1}", implicit=False)

            if isMC:
                lpair = f"(CSC,{second.upper()})"
            else:
                lpair = f"(CSC,{second.upper()}_"  # + "{OOT})"

            hhs["csc" + second + "_dEta"].append(
                _ms.Histo1D(("csc" + second + "_dEta", ";#Delta#eta" + lpair + ";count", *bins["AdEta"]), "tag_dEta"))
            hhs["csc" + second + "_dPhi"].append(
                _ms.Histo1D(("csc" + second + "_dPhi", ";#Delta#phi" + lpair + ";count", *bins["AdPhi"]), "tag_dPhi"))
            hhs["csc" + second + "_dR"].append(
                _ms.Histo1D(("csc" + second + "_dR", ";#DeltaR" + lpair + ";count", *bins["dR"]), "tag_dR"))

            # d__ vs d__
            hhs["csc" + second + "_2D_dPhi_dtSize"].append(
                _ms.Histo2D(("csc" + second + "_2D_dPhi_dtSize", ";#Delta#phi" + lpair + ";DT Size;count",
                             *bins["AdPhi_dtSizeSmall"]), "tag_dPhi", "dtRechitClusterSize"))
            hhs["csc" + second + "_2D_dPhi_dEta"].append(
                _ms.Histo2D(("csc" + second + "_2D_dPhi_dEta", ";#Delta#phi" + lpair + ";#Delta#eta" + lpair + ";count",
                             *bins["AdPhi_AdEta"]), "tag_dPhi", "tag_dEta"))
            hhs["csc" + second + "_2D_dPhi_dR"].append(
                _ms.Histo2D(("csc" + second + "_2D_dPhi_dR", ";#Delta#phi" + lpair + ";#DeltaR" + lpair + ";count",
                             *bins["AdPhi_dR"]), "tag_dPhi", "tag_dR"))
            hhs["csc" + second + "_2D_dEta_dR"].append(
                _ms.Histo2D(("csc" + second + "_2D_dEta_dR", ";#Delta#eta" + lpair + ";#DeltaR" + lpair + ";count",
                             *bins["AdEta_dR"]), "tag_dEta", "tag_dR"))

            #
            # hhs["csc" + second + "_2D_dPhi_met"].append(
            #     _ms.Histo2D(("csc" + second + "_2D_dPhi_met", ";#Delta#phi" + lpair + ";E_{T}^{miss} [GeV];count",
            #                  *bins["AdPhi_met"]), "tag_dPhi", "met"))
            # hhs["csc" + second + "_2D_dEta_met"].append(
            #     _ms.Histo2D(("csc" + second + "_2D_dEta_met", ";#Delta#eta" + lpair + ";E_{T}^{miss} [GeV];count",
            #                  *bins["AdEta_met"]), "tag_dEta", "met"))
            # hhs["csc" + second + "_2D_dR_met"].append(
            #     _ms.Histo2D(
            #         ("csc" + second + "_2D_dR_met", ";#DeltaR" + lpair + ";E_{T}^{miss} [GeV];count", *bins["dR_met"]),
            #         "tag_dR", "met"))

        # ###############################################
        # ## Save data selections for further analysis ##
        # ###############################################
        # alert("SAVING DATA SELECTIONS")

        # #=====================================================#
        # #  SAVING JET and MUON DATA TO DEVELOP BETTER MB1 Cut #
        # #=====================================================#
        # print(f'    - dt_mb1_jet_muon_{stat}')

        # if isMC:
        #     # require DT matched and decayed DT
        #     sms = oms.match_clusters('cscdt', implicit=False).time_cut('it').L1_plateau()  #.Filter(
        #     # 'dtRechitClusterSize > 100', 'dt')
        #     sms.define_2tag_kins_and_cut('cscdt').Filter('nDtRechitClusters > 0')

        #     # sms = oms.match_clusters('dt', implicit=False).Filter('nDtRechitClusters > 0')
        # else:
        #     # require HLT DT trigger, blind using CSC+DT IT, DT>100, 1<dPhi<2
        #     sms = oms.time_cut('it', implicit=False)  #.Filter('dtRechitClusterSize > 100', 'dt')
        #     sms.define_2tag_kins_and_cut('cscdt').Filter('nDtRechitClusters > 0')
        #     sms.Filter("(0<tag_dPhi)&&(tag_dPhi<2)")  #&&(HLTDecision[569] == 1)")

        #     # sms = oms.time_cut('oot', 'dt', implicit=False).Filter('dtRechitClusterSize > 100',
        #     #                                                        'dt').Filter('nDtRechitClusters > 0')

        # build_hists(sms, hhs, 'sel_mb1')

        # col_names = [
        #     "nDtRechitClusters", "dtRechitClusterSize", "dtRechitClusterNSegStation1", "dtRechitClusterJetVetoLooseId",
        #     "dtRechitClusterJetVetoPt", "dtRechitClusterMuonVetoLooseId", "dtRechitClusterMuonVetoPt"
        # ]

        # sms.rdf.Snapshot('ReducedMuonSystem', out_data_dir + f"dt_mb1_jet_muon_{stat}_{save_date}.root",
        #                  col_names)  #, ssOpt)
        # print(f'          n events = {sms.Count():,}')
        #=====================================================#

        # if save_dstat in stat:
        #     print("Saving data to disk.")
        #     header_out = ['met','dR','dEta','dPhi', 'isMC', 'isCscCsc']
        #     data_out = np.c_[met, dR, dEta, dPhi, isMC*np.ones_like(met, dtype=bool), (second=='csc')*np.ones_like(met,dtype=bool)]
        #     np.savetxt(out_data_dir + "csc"+second+"_processed_"+stat+".csv", data_out, delimiter=',', comments='', header=','.join(header_out))

        ###############################################
        ## Save plots to disk and make simple histos ##
        ###############################################

        # TODO: save plots to a .root file
        alert("Saving plots to disk.")
        print("    - Note: this will call all queued lazy operations!")
        for k, hs in hhs.items():
            hs = [h.GetPtr() for h in hs]
            for i, ss in enumerate(stats[:len(hs)]):
                fout_root.WriteObject(hs[i], k + f"_{ss}")

            if "2D" in k:
                for i, ss in enumerate(stats[:len(hs)]):
                    c = canvas(1, 1)
                    c.cd(1).SetLogz()
                    hs[i].Draw("colz")
                    c.Print(out_dir + k + f"_{ss}" + ending)
            else:
                c = canvas(1, 1)
                _hll = multi_plot(hs[::-1], labels[::-1], ymin=0, norm=True)
                c.Print(out_dir + k + ending)

    fout_root.Close()
    # #########################
    # ## 1D Efficiency Plots ##
    # #########################
    # print("Generating efficiency plots")

    # for isCut in [False, True]:
    #     c = canvas(2, 2)
    #     for i, det in enumerate(["csc", "dt"]):
    #         for j, xl in enumerate(["z", "r"]):
    #             c.cd(2 * i + j + 1)
    #             _hhs = []
    #             for ff, ss in zip(files_mc, stats_mc):
    #                 ms = MuonSystem(ff, isMC=True, nev=nev)
    #                 _hhs.append(make_cluster_eff_1D(ms, det, xl, cuts=isCut))
    #                 # Write Hists to memory
    #                 #fout_root.WriteObject(_hhs[-1], f"eff_{det}_{xl}" + ("_CUT" if isCut else "")+f"_{ss}")

    #             if isCut:
    #                 legxy = (0.3, 0.7, 0.5, 0.85)
    #             else:
    #                 legxy = (0.5, 0.15, 0.7, 0.4)
    #             _ll = multi_plot(_hhs, labels_mc, legxy=legxy)
    #             for h in _hhs:
    #                 h.SetMinimum(0.0)
    #                 h.SetMaximum(1.0)

    #             if det == "dt" and xl == "r":
    #                 boxes = draw_dt_r_boxes(_hhs[-1])
    #                 gc.extend(boxes)
    #             if det == "csc" and xl == "z":
    #                 boxes = draw_csc_z_boxes(_hhs[-1])
    #                 gc.extend(boxes)
    #             gc.extend(_hhs + list(_ll))

    #     c.Print(out_dir + "eff_compare_DBSCAN-CA" + ("_CUT" if isCut else "") + ending)

    # print("Finished!")
