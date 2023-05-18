"""For analyzing multiple files at once (e.g. comparing clustering algorithms)"""

import os
import sys
import pathlib
from collections import defaultdict

import ROOT as rt
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
    "dtTime": (6, -3.5, 3.5),
    #'dtTimeSpread' : (100, -100, 100),
    #'dtE' : (100, -100, 100),
    #'dtPt' : (100, -100, 100),
    "dtEta": (100, 0, 1.5),
    "dtPhi": (100, -pi, pi),
    "dtMetDPhi": (100, -pi, pi),
    "dtNSt": (9, -0.5, 8.5),
    "dtASt": (9, -0.5, 8.5),
    "dtMSt": (9, -0.5, 8.5),
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
    "dEta": (100, -5, 5),
    "dPhi": (100, -pi, pi),
    "dR": (100, 0, 5),
    #
    "AdEta": (100, 0, 5),
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
    "AdPhi_AdEta": bins["AdPhi"] + bins["AdEta"],
    "AdPhi_dR": bins["AdPhi"] + bins["dR"],
    "AdEta_dR": bins["AdEta"] + bins["dR"],
    "AdPhi_met": bins["AdPhi"] + bins["met"],
    "AdEta_met": bins["AdEta"] + bins["met"],
    "cmsZ_cmsR": bins["cmsZ"] + bins["cmsR"],
}

bins = {**bins, **bins2D}

if __name__ == "__main__":
    isBlind = True

    gc = []
    rt.gROOT.SetBatch()

    cur_dir = os.getcwd()
    if 'Documents' in cur_dir:
        rt.EnableImplicitMT()
    else:
        rt.EnableImplicitMT(8)

    print(f"Running with thread pool size = {rt.GetThreadPoolSize():,}")

    a = tdrstyle.setTDRStyle()
    # CMS_lumi.writeExtraText = 0
    rt.gStyle.SetOptFit(0)  # 1011)

    out_dir = cur_dir + "/reports/weekly/may15/"
    in_data_dir = cur_dir + "/data/raw/"
    out_data_dir = cur_dir + "/data/processed/"
    ending = ".png"

    isCut = True
    save_dstat = "ca_0p6"
    if len(sys.argv) > 1:
        nev = int(sys.argv[1])
    else:
        nev = 100_000

    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)  # make out directory if it doesn't exist
    pathlib.Path(out_data_dir).mkdir(parents=True, exist_ok=True)  # make out data directory if it doesn't exist
    print(f"Using output directory '{out_dir}'")

    fout_name = out_data_dir + "processed_data.root"
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

    hhs = defaultdict(list)

    for ff, label, stat in zip(files, labels, stats):
        print(f"\nLoading {stat} (N Max Events = {nev:,})")  # ms.file_name)
        isMC = "mc" in stat
        ms = MuonSystemRDF(ff, isMC=isMC, nev=nev)
        print(f"\tLoaded {ms.Count():,} events.")

        # # ---------------------------------------- #
        # # ---------------------------------------- #
        # # ---------------------------------------- #

        # No cuts		 1118331.31	 100.000
        # 2 CSC + 0 lep	 15543.23	 1.390
        # ME1 veto [0]	 7406.14	 0.662
        # ME1 veto [1]	 4898.57	 0.438
        # ME1 veto [0&1]	 2461.67	 0.220
        # Trigger plateau	 6951.87	 0.622
        # Trigger + ME1	 1277.05	 0.114
        # Min dphi>1.8	 5470.03	 0.489
        # Delta eta<1.9	 4833.32	 0.432
        # In time cut	 4566.32	 0.408
        # Time spread	 4339.78	 0.388
        # Jet muon vetoes	 3494.06	 0.312
        # ABCD plane	 3494.06	 0.31
        # C	 822.77	 0.07
        # B	 953.86	 0.09
        # D	 767.19	 0.07
        # A	 950.23	 0.08

        cscdt = "(nCscRechitClusters == 1) && (nDtRechitClusters == 1) && (nLeptons == 0)"
        # ##
        # print("Not implicit")

        print(r"\begin{table}[]")
        print(r"\begin{tabular}{c|rr}")
        print(r"Selection & Yield & Eff. vs no cuts \\ \hline")

        hhs["met_all_ni"].append(ms.Histo1D(('meta_ni', ';met, all;', *bins['met']), 'met'))
        yd, nc = ms.Count(), ms.Count()
        print(f"All    & {yd:,} & {yd/nc*100:,.3f}\\% \\\\")

        #! THIS ONE IS IMPLICIT!!!!
        ms2 = ms.Filter("(nCscRechitClusters == 1) && (nDtRechitClusters == 1) && (nLeptons == 0)", implicit=True)
        yd = ms2.Count()
        print(f"1 CSC/DT + 0 lep  & {yd:,} & {yd/nc*100:,.3f}\\% \\\\")
        nc = yd
        #! !!!!!!!!!!!!!!!!!!!!!!!!

        ms2 = ms.Filter(
            "(cscRechitClusterNRechitChamberMinus11 + cscRechitClusterNRechitChamberMinus12 + cscRechitClusterNRechitChamberPlus11 + cscRechitClusterNRechitChamberPlus12 == 0)",
            'csc',
            implicit=False).Filter(cscdt)
        yd = ms2.Count()
        print(f"ME1 veto [CSC]  & {yd:,} & {yd/nc*100:,.3f}\\% \\\\")

        ms2 = ms.Filter("dtRechitClusterNSegStation1 == 0", 'dt', implicit=False).Filter(cscdt)
        yd = ms2.Count()
        print(f"MB1 veto [DT]  & {yd:,} & {yd/nc*100:,.3f}\\% \\\\")

        ms2 = ms.Filter(
            "(cscRechitClusterNRechitChamberMinus11 + cscRechitClusterNRechitChamberMinus12 + cscRechitClusterNRechitChamberPlus11 + cscRechitClusterNRechitChamberPlus12 == 0)",
            'csc',
            implicit=False).Filter("(dtRechitClusterNSegStation1 == 0)", 'dt').Filter(cscdt)
        yd = ms2.Count()
        print(f"ME1 + MB1 veto [CSC&DT]  & {yd:,} & {yd/nc*100:,.3f}\\% \\\\")

        ms2 = ms2.time_cut("it", "csc", implicit=False).Filter(cscdt)
        yd = ms2.Count()
        print(f"CSC In-time & {yd:,} & {yd/nc*100:,.3f}\\% \\\\")

        ms2 = ms2.time_cut("it", "dt", implicit=False).Filter(cscdt)
        yd = ms2.Count()
        print(f"DT In-time & {yd:,} & {yd/nc*100:,.3f}\\% \\\\")

        ms2 = ms2.time_cut("it", "cscdt", implicit=False).Filter(cscdt)
        yd = ms2.Count()
        print(f"CSC + DT In-time & {yd:,} & {yd/nc*100:,.3f}\\% \\\\")

        ms2 = ms.jet_cut(implicit=False).muon_cut().Filter(cscdt)
        yd = ms2.Count()
        print(f"Jet + Muon veto & {yd:,} & {yd/nc*100:,.3f}\\% \\\\")

        print(r"\end{tabular}")
        print(r"\end{table}")
        print("")

        # ##

        # print("Implicit")

        # print(r"\begin{table}[]")
        # print(r"\begin{tabular}{c|rr}")
        # print(r"Selection & Yield & Eff. vs no cuts \\ \hline")

        # hhs["met_all"].append(ms.Histo1D(('meta', ';met, all;', *bins['met']), 'met'))
        # yd, nc = ms.Count(), ms.Count()
        # print(f"All    & {yd:,} & {yd/nc*100:,.3f}\\% \\\\")

        # ms.Filter("(nCscRechitClusters == 1) & (nDtRechitClusters == 1)")
        # hhs["met_2tag"].append(ms.Histo1D(('met2tag', ';met, CSC+DT;', *bins['met']), 'met'))
        # yd = ms.Count()
        # print(f"+2 Tag  & {yd:,} & {yd/nc*100:,.3f}\\% \\\\")

        # ms.jet_cut()
        # ms.Filter("(nCscRechitClusters == 1) & (nDtRechitClusters == 1)")
        # yd = ms.Count()
        # print(f"+Jet    & {yd:,} & {yd/nc*100:,.3f}\\% \\\\")

        # ms.muon_cut()
        # ms.Filter("(nCscRechitClusters == 1) & (nDtRechitClusters == 1)")
        # yd = ms.Count()
        # print(f"+Muon   & {yd:,} & {yd/nc*100:,.3f}\\% \\\\")

        # ms.Filter("(nCscRechitClusters == 1) & (nDtRechitClusters == 1)")
        # ms.time_cut("oot", "dt")
        # hhs["met_2t_DtOot"].append(ms.Histo1D(('met2tDtOot', ';met, CSC+DT_{OOT};', *bins['met']), 'met'))
        # yd = ms.Count()
        # print(f"+DT OOT & {yd:,} & {yd/nc*100:,.3f}\\% \\\\")
        # print(r"\end{tabular}")
        # print(r"\end{table}")
        # print("")

        # ##

        # ms_csc = ms.time_cut('oot', 'csc', implicit=False).Filter("nCscRechitClusters + nDtRechitClusters > 0")
        # hCsc = ms_csc.Histo1D(("nCscRechitClusters", ";N CSC Clusters;count", *bins["cscN"]), "nCscRechitClusters")

        # ms_cscT = ms.time_cut('oot', 'cscT', implicit=False).Filter("nCscRechitClusters + nDtRechitClusters > 0")
        # hCscT = ms_cscT.Histo1D(("nCscRechitClusters", ";N CSC Clusters;count", *bins["cscN"]), "nCscRechitClusters")

        # c = canvas(1, 1)
        # _hll = multi_plot([hCsc, hCscT], ['csc', 'cscT'], ymin=0)  #, norm=True)
        # c.Print(out_dir + 'csc_oot_test_' + stat + ending)

        # ##

        ms = MuonSystemRDF(ff, isMC=isMC, nev=nev)
        # # ---------------------------------------- #
        # # ---------------------------------------- #
        # # ---------------------------------------- #

        # ===================================== #
        # Making time plots before cuts on time #
        # ===================================== #

        hhs["cscAllTime"].append(
            ms.Histo1D(("cscAllTime", ";CSC Time [ns];count", *bins["cscTime"]), "cscRechitClusterTimeWeighted"))
        hhs["cscAllTimeSpread"].append(
            ms.Histo1D(("cscAllTimeSpread", ";CSC Time Spread [ns];count", *bins["cscTimeSpread"]),
                       "cscRechitClusterTimeSpreadWeightedAll"))
        hhs["dtAllTime"].append(
            ms.Histo1D(("dtAllTime", ";DT Time [ns];count", *bins["dtTime"]), "dtRechitCluster_match_RPCBx_dPhi0p5"))

        ##############################
        ## CUT, BLIND, & MATCH DATA ##
        ##############################

        if isCut:
            ms.jet_cut()
            ms.muon_cut()

        if isMC:
            print("!!!!!!!!!!!!!!!!!!!")
            print("!! MATCHING DATA !!")
            print("!!!!!!!!!!!!!!!!!!!")
            ms.match_clusters()

        if isBlind and not isMC:
            print("!!!!!!!!!!!!!!!!!!!")
            print("!! BLINDING DATA !!")
            print("!!!!!!!!!!!!!!!!!!!")
            print("    - Requiring DT to be OOT")
            print("    - Requiring CSC (not trigger) to be OOT")
            ms.time_cut('oot', 'dt cscT')
            # ms.L1_plateau()

        ms.Filter("nCscRechitClusters + nDtRechitClusters > 0")

        # ================================ #
        # NOTE: Using convention "x vs y"  #
        # ================================ #

        #################################
        ## System & Cluster Kinematics ##
        #################################
        print("Generating system and cluster kinematics plots")

        ms.Define("ABScscRechitClusterAvgStation10", "abs(cscRechitClusterAvgStation10)")
        ms.Define("ABScscRechitClusterMaxStation", "abs(cscRechitClusterMaxStation)")
        ms.Define("ABScscRechitClusterEta", "abs(cscRechitClusterEta)")
        ms.Define("ABSdtRechitClusterAvgStation10", "abs(dtRechitClusterAvgStation10)")
        ms.Define("ABSdtRechitClusterMaxStation", "abs(dtRechitClusterMaxStation)")
        ms.Define("ABSdtRechitClusterEta", "abs(dtRechitClusterEta)")

        # System level
        hhs["met"].append(ms.Histo1D(("met", ";E_{T}^{miss} [GeV];count", *bins["met"]), "met"))
        hhs["metPhi"].append(ms.Histo1D(("metPhi", ";#phi[E_{T}^{miss}];count", *bins["metPhi"]), "metPhi"))

        # kinematics: (csc, dt, gllp, jet, lep)
        # N, Size, Time, TimeSpread, NSt, ASt, MSt, MetDPhi, E, Pt, Eta, Phi
        hhs["cscN"].append(
            ms.Histo1D(("nCscRechitClusters", ";N CSC Clusters;count", *bins["cscN"]), "nCscRechitClusters"))
        hhs["cscSize"].append(
            ms.Histo1D(("cscRechitClusterSize", ";CSC Rechit Size;count", *bins["cscSize"]), "cscRechitClusterSize"))
        hhs["cscTime"].append(
            ms.Histo1D(("cscRechitClusterTimeWeighted", ";CSC Time [ns];count", *bins["cscTime"]),
                       "cscRechitClusterTimeWeighted"))
        hhs["cscTimeSpread"].append(
            ms.Histo1D(("cscTimeSpread", ";CSC Time Spread [ns];count", *bins["cscTimeSpread"]),
                       "cscRechitClusterTimeSpreadWeightedAll"))
        hhs["cscNSt"].append(
            ms.Histo1D(("cscRechitClusterNStation10", ";N Stations / CSC Cluster;count", *bins["cscNSt"]),
                       "cscRechitClusterNStation10"))
        hhs["cscASt"].append(
            ms.Histo1D(("cscASt", ";Avg Station / CSC Cluster;count", *bins["cscASt"]),
                       "ABScscRechitClusterAvgStation10"))
        hhs["cscMSt"].append(
            ms.Histo1D(("cscMSt", ";Max Station / CSC Cluster;count", *bins["cscMSt"]),
                       "ABScscRechitClusterMaxStation"))
        hhs["cscMetDPhi"].append(
            ms.Histo1D(("cscMetDPhi", ";#Delta#phi(CSC,E_{T}^{miss});count", *bins["cscMetDPhi"]),
                       "cscRechitClusterMet_dPhi"))
        hhs["cscEta"].append(
            ms.Histo1D(("cscRechitClusterEta", ";CSC #eta;count", *bins["cscEta"]), "ABScscRechitClusterEta"))
        hhs["cscPhi"].append(
            ms.Histo1D(("cscRechitClusterPhi", ";CSC #phi;count", *bins["cscPhi"]), "cscRechitClusterPhi"))

        hhs["dtN"].append(ms.Histo1D(("nDtRechitClusters", ";N DT Clusters;count", *bins["dtN"]), "nDtRechitClusters"))
        hhs["dtSize"].append(
            ms.Histo1D(("dtRechitClusterSize", ";DT Rechit Size;count", *bins["dtSize"]), "dtRechitClusterSize"))
        hhs["dtTime"].append(
            ms.Histo1D(("dtRechitCluster_match_RPCBx_dPhi0p5", ";DT Time [ns];count", *bins["dtTime"]),
                       "dtRechitCluster_match_RPCBx_dPhi0p5"))
        hhs["dtNSt"].append(
            ms.Histo1D(("dtRechitClusterNStation10", ";N Stations / DT Cluster;count", *bins["dtNSt"]),
                       "dtRechitClusterNStation10"))
        hhs["dtASt"].append(
            ms.Histo1D(("dtRechitClusterAvgStation10", ";Avg Station / DT Cluster;count", *bins["dtASt"]),
                       "ABSdtRechitClusterAvgStation10"))
        hhs["dtMSt"].append(
            ms.Histo1D(("dtRechitClusterMaxStation", ";Max Station / DT Cluster;count", *bins["dtMSt"]),
                       "ABSdtRechitClusterMaxStation"))
        hhs["dtMetDPhi"].append(
            ms.Histo1D(("dtRechitClusterMet_dPhi", ";#Delta#phi(DT,E_{T}^{miss});count", *bins["dtMetDPhi"]),
                       "dtRechitClusterMet_dPhi"))
        hhs["dtEta"].append(
            ms.Histo1D(("dtRechitClusterEta", ";DT #eta;count", *bins["dtEta"]), "ABSdtRechitClusterEta"))
        hhs["dtPhi"].append(ms.Histo1D(("dtRechitClusterPhi", ";DT #phi;count", *bins["dtPhi"]), "dtRechitClusterPhi"))

        hhs["gllpN"].append(ms.Histo1D(("nGLLP", ";N gLLPs;count", *bins["gllpN"]), "nGLLP"))
        hhs["gllpE"].append(ms.Histo1D(("gLLP_e", ";gLLP E [GeV];count", *bins["gllpE"]), "gLLP_e"))
        hhs["gllpPt"].append(ms.Histo1D(("gLLP_pt", ";gLLP P_{T} [GeV];count", *bins["gllpPt"]), "gLLP_pt"))
        hhs["gllpEta"].append(ms.Histo1D(("gLLP_eta", ";gLLP #eta;count", *bins["gllpEta"]), "gLLP_eta"))
        hhs["gllpPhi"].append(ms.Histo1D(("gLLP_phi", ";gLLP #phi;count", *bins["gllpPhi"]), "gLLP_phi"))

        hhs["jetN"].append(ms.Histo1D(("nJets", ";N Jets;count", *bins["jetN"]), "nJets"))
        hhs["jetE"].append(ms.Histo1D(("jetE", ";Jet E [GeV];count", *bins["jetE"]), "jetE"))
        hhs["jetPt"].append(ms.Histo1D(("jetPt", ";Jet P_{T} [GeV];count", *bins["jetPt"]), "jetPt"))
        hhs["jetEta"].append(ms.Histo1D(("jetEta", ";Jet #eta;count", *bins["jetEta"]), "jetEta"))
        hhs["jetPhi"].append(ms.Histo1D(("jetPhi", ";Jet #phi;count", *bins["jetPhi"]), "jetPhi"))

        hhs["lepN"].append(ms.Histo1D(("nLeptons", ";N Leptons;count", *bins["lepN"]), "nLeptons"))
        hhs["lepE"].append(ms.Histo1D(("lepE", ";Lepton E [GeV];count", *bins["lepE"]), "lepE"))
        hhs["lepPt"].append(ms.Histo1D(("lepPt", ";Lepton P_{T} [GeV];count", *bins["lepPt"]), "lepPt"))
        hhs["lepEta"].append(ms.Histo1D(("lepEta", ";Lepton #eta;count", *bins["lepEta"]), "lepEta"))
        hhs["lepPhi"].append(ms.Histo1D(("lepPhi", ";Lepton #phi;count", *bins["lepPhi"]), "lepPhi"))

        # MET vs (csc, dt) MET dPhi
        hhs["2D_met_cscMetDPhi"].append(
            ms.Histo2D(
                ("2D_met_cscMetDPhi", ";E_{T}^{miss} [GeV];#Delta#phi(csc, #phi_{T}^{miss});", *bins["met_cscMetDPhi"]),
                "met", "cscRechitClusterMet_dPhi"))
        hhs["2D_met_dtMetDPhi"].append(
            ms.Histo2D(("2D_met_dtMetDPhi", ";E_{T}^{miss} [GeV];#Delta#phi(dt, #phi_{T}^{miss});count",
                        *bins["met_dtMetDPhi"]), "met", "dtRechitClusterMet_dPhi"))

        ###########################################
        ## Single Pair CSC-CSC/CSC-DT Kinematics ##
        ###########################################

        print("Reducing to 2tag analysis!")
        ms.define_tag_kins(system="csccsc,cscdt")

        for second in ('csc', 'dt'):
            print(f"Generating CSC-{second.upper()} kinematics plots")
            _ms = ms.Filter(f"tag_type == {0 if second == 'csc' else 1}", implicit=False)

            if isMC:
                lpair = f"(CSC,{second.upper()})"
            else:
                lpair = f"(CSC,{second.upper()}_" + "{OOT})"
                if second == 'csc':
                    _ms.Filter("HLTDecision[566] > 0 ")  # HLT_CscCluster_Loose
                elif second == 'dt':
                    _ms.Filter("HLTDecision[569] > 0 ")  # HLT_L1CSCCluster_DTCluster50

            hhs["csc" + second + "_dEta"].append(
                _ms.Histo1D(("csc" + second + "_dEta", ";#Delta#eta" + lpair + ";count", *bins["AdEta"]), "tag_dEta"))
            hhs["csc" + second + "_dPhi"].append(
                _ms.Histo1D(("csc" + second + "_dPhi", ";#Delta#phi" + lpair + ";count", *bins["AdPhi"]), "tag_dPhi"))
            hhs["csc" + second + "_dR"].append(
                _ms.Histo1D(("csc" + second + "_dR", ";#DeltaR" + lpair + ";count", *bins["dR"]), "tag_dR"))

            # d__ vs d__
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
            hhs["csc" + second + "_2D_dPhi_met"].append(
                _ms.Histo2D(("csc" + second + "_2D_dPhi_met", ";#Delta#phi" + lpair + ";E_{T}^{miss} [GeV];count",
                             *bins["AdPhi_met"]), "tag_dPhi", "met"))
            hhs["csc" + second + "_2D_dEta_met"].append(
                _ms.Histo2D(("csc" + second + "_2D_dEta_met", ";#Delta#eta" + lpair + ";E_{T}^{miss} [GeV];count",
                             *bins["AdEta_met"]), "tag_dEta", "met"))
            hhs["csc" + second + "_2D_dR_met"].append(
                _ms.Histo2D(
                    ("csc" + second + "_2D_dR_met", ";#DeltaR" + lpair + ";E_{T}^{miss} [GeV];count", *bins["dR_met"]),
                    "tag_dR", "met"))

        #     ###############################################
        #     ## Save data selections for further analysis ##
        #     ###############################################
        #     if save_dstat in stat:
        #         print("Saving data to disk.")
        #         header_out = ['met','dR','dEta','dPhi', 'isMC', 'isCscCsc']
        #         data_out = np.c_[met, dR, dEta, dPhi, isMC*np.ones_like(met, dtype=bool), (second=='csc')*np.ones_like(met,dtype=bool)]
        #         np.savetxt(out_data_dir + "csc"+second+"_processed_"+stat+".csv", data_out, delimiter=',', comments='', header=','.join(header_out))

        ###############################################
        ## Save plots to disk and make simple histos ##
        ###############################################
        # TODO: save plots to a .root file
        print("Saving plots to disk.")
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
                _hll = multi_plot(hs, labels, ymin=0, norm=True)
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
