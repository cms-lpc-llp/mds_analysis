"""For analyzing multiple files at once (e.g. comparing clustering algorithms)"""

import os
import sys
import pathlib
from collections import defaultdict

import ROOT as rt
import numpy as np
import sklearn as skl
import awkward as ak

from src import CMS_lumi, tdrstyle
from src.muon_system import MuonSystem, MuonSystemRDF
from src.helper_functions import lnot, land, lor, lxor, asum, aabs

from src.muon_system import (get_lat_leg, H1D, H2D, multi_plot, make_cluster_eff_1D, draw_csc_z_boxes, draw_dt_r_boxes)
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
    "dEta": (100, -4, 4),
    "dPhi": (100, -pi, pi),
    "dR": (100, 0, 5),
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
    "cmsZ_cmsR": bins["cmsZ"] + bins["cmsR"],
}

bins = {**bins, **bins2D}

if __name__ == "__main__":
    isBlind = True

    gc = []
    rt.gROOT.SetBatch()

    a = tdrstyle.setTDRStyle()
    # CMS_lumi.writeExtraText = 0
    rt.gStyle.SetOptFit(0)  # 1011)

    cur_dir = os.getcwd()
    out_dir = cur_dir + "/reports/weekly/may8/"
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
        print(f"\nLoading {stat}. (N Max Events = {nev:,})")  # ms.file_name)
        isMC = "mc" in stat
        ms = MuonSystemRDF(ff, isMC=isMC, nev=nev)
        print(f"\tLoaded {ms.Count():,} events.")

        ######################

        print("Not implicit")

        print(r"\begin{table}[]")
        print(r"\begin{tabular}{c|rr}")
        print(r"Selection & Yield & Eff. vs no cuts (%) \\ \hline")

        ms2 = ms.Filter("met != -3.14159")  # Useless cut just to make a copy
        hhs["met_all"].append(ms2.Histo1D(('meta', ';met, all;', *bins['met']), 'met'))
        yd, nc = ms2.Count(), ms2.Count()
        print(f"All    & {yd:,} & {yd/nc*100:,}% \\")

        ms2.Filter("(nCscRechitClusters == 1) & (nDtRechitClusters == 1)", implicit=False)
        hhs["met_2tag"].append(ms2.Histo1D(('met2tag', ';met, CSC+DT;', *bins['met']), 'met'))
        yd = ms2.Count()
        print(f"2 Tag  & {yd:,} & {yd/nc*100:,}% \\")

        ms2.jet_cut(implicit=False)
        ms2.Filter("(nCscRechitClusters == 1) & (nDtRechitClusters == 1)", implicit=False)
        yd = ms2.Count()
        print(f"Jet    & {yd:,} & {yd/nc*100:,}% \\")

        ms2.muon_cut(implicit=False)
        ms2.Filter("(nCscRechitClusters == 1) & (nDtRechitClusters == 1)", implicit=False)
        yd = ms2.Count()
        print(f"Muon   & {yd:,} & {yd/nc*100:,}% \\")

        ms2.Filter("(nCscRechitClusters == 1) & (nDtRechitClusters == 1)", implicit=False)
        ms2.time_cut("oot", "dt", implicit=False)
        hhs["met_2t_DtOot"].append(ms2.Histo1D(('met2tDtOot', ';met, CSC+DT_{OOT};', *bins['met']), 'met'))
        yd = ms2.Count()
        print(f"DT OOT & {yd:,} & {yd/nc*100:,}% \\")
        print(r"\end{tabular}")
        print(r"\end{table}")
        print("")

        ######################

        print("Implicit")

        print(r"\begin{table}[]")
        print(r"\begin{tabular}{c|rr}")
        print(r"Selection & Yield & Eff. vs no cuts (%) \\ \hline")

        hhs["met_all"].append(ms.Histo1D(('meta', ';met, all;', *bins['met']), 'met'))
        yd, nc = ms.Count(), ms.Count()
        print(f"All    & {yd:,} & {yd/nc*100:,}% \\")

        ms.Filter("(nCscRechitClusters == 1) & (nDtRechitClusters == 1)")
        hhs["met_2tag"].append(ms.Histo1D(('met2tag', ';met, CSC+DT;', *bins['met']), 'met'))
        yd = ms.Count()
        print(f"+2 Tag  & {yd:,} & {yd/nc*100:,}% \\")

        ms.jet_cut()
        ms.Filter("(nCscRechitClusters == 1) & (nDtRechitClusters == 1)")
        yd = ms.Count()
        print(f"+Jet    & {yd:,} & {yd/nc*100:,}% \\")

        ms.muon_cut()
        ms.Filter("(nCscRechitClusters == 1) & (nDtRechitClusters == 1)")
        yd = ms.Count()
        print(f"+Muon   & {yd:,} & {yd/nc*100:,}% \\")

        ms.Filter("(nCscRechitClusters == 1) & (nDtRechitClusters == 1)")
        ms.time_cut("oot", "dt")
        hhs["met_2t_DtOot"].append(ms.Histo1D(('met2tDtOot', ';met, CSC+DT_{OOT};', *bins['met']), 'met'))
        yd = ms.Count()
        print(f"+DT OOT & {yd:,} & {yd/nc*100:,}% \\")
        print(r"\end{tabular}")
        print(r"\end{table}")
        print("")

        ######################

        # # =======================================#
        # # Making time plots before cuts on time #
        # # =======================================#

        # print("here")
        # hhs["cscAllTime"].append(H1D(ms["cscRechitClusterTimeWeighted"], ";CSC Time [ns];count", bins["cscTime"]))
        # hhs["cscAllTimeSpread"].append(
        #     H1D(ms["cscRechitClusterTimeSpreadWeightedAll"], ";CSC Time Spread [ns];count", bins["cscTimeSpread"])
        # )
        # hhs["dtAllTime"].append(H1D(ms["dtRechitCluster_match_RPCBx_dPhi0p5"], ";DT Time [ns];count", bins["dtTime"]))

        # ##############################
        # ## CUT, BLIND, & MATCH DATA ##
        # ##############################

        # if isCut:
        #     sel_csc_jet, sel_dt_jet = ms.jet_veto_cut()
        #     ms.apply_cut(sel_csc_jet, system="csc")
        #     ms.apply_cut(sel_dt_jet, system="dt")
        #     ms.apply_cut(ms["nCscRechitClusters"] + ms["nDtRechitClusters"] > 0, system="event")

        #     sel_csc_muon, sel_dt_muon = ms.muon_veto_cut()
        #     ms.apply_cut(sel_csc_muon, system="csc")
        #     ms.apply_cut(sel_dt_muon, system="dt")
        #     ms.apply_cut(ms["nCscRechitClusters"] + ms["nDtRechitClusters"] > 0, system="event")
        # if isMC:
        #     print("!!!!!!!!!!!!!!!!!!!")
        #     print("!! MATCHING DATA !!")
        #     print("!!!!!!!!!!!!!!!!!!!")
        #     idx_csc_match, idx_dt_match = ms.match_cut()
        #     ms.apply_cut(idx_csc_match, system="csc")
        #     ms.apply_cut(idx_dt_match, system="dt")

        # elif isBlind:
        #     print("!!!!!!!!!!!!!!!!!!!")
        #     print("!! BLINDING DATA !!")
        #     print("!!!!!!!!!!!!!!!!!!!")

        #     # sel_csc_oot = aabs(ms["cscRechitClusterTimeWeighted"]) > 50
        #     sel_csc_oot = ms["cscRechitClusterTimeWeighted"] < -12.5
        #     sel_dt_oot = ms["dtRechitCluster_match_RPCBx_dPhi0p5"] != 0

        #     sel_bkg = lor(
        #         land(ms["nCscRechitClusters"] == 1, asum(sel_csc_oot) == 1),
        #         land(ms["nCscRechitClusters"] == 2, asum(sel_csc_oot) == 1),  # NOTE: change to req only 2nd clstr oot
        #         land(ms["nDtRechitClusters"] == 1, asum(sel_dt_oot) == 1),
        #         land(ms["nDtRechitClusters"] == 2, asum(sel_dt_oot) == 1),
        #         land(ms["nCscRechitClusters"] == 1, ms["nDtRechitClusters"] == 1, asum(sel_dt_oot) == 1),
        #     )

        #     ms.apply_cut(sel_bkg)

        # # ================================ #
        # # NOTE: Using convention "x vs y"  #
        # # ================================ #

        # #################################
        # ## System & Cluster Kinematics ##
        # #################################
        # print("Generating system and cluster kinematics plots")

        # # System level
        # hhs["met"].append(H1D(ms["met"], ";E_{T}^{miss} [GeV];count", bins["met"]))
        # hhs["metPhi"].append(H1D(ms["metPhi"], ";#phi[E_{T}^{miss}];count", bins["metPhi"]))

        # # kinematics: (csc, dt, gllp, jet, lep)
        # # N, Size, Time, TimeSpread, NSt, ASt, MSt, MetDPhi, E, Pt, Eta, Phi
        # hhs["cscN"].append(H1D(ms["nCscRechitClusters"], ";N CSC Clusters;count", bins["cscN"]))
        # hhs["cscSize"].append(H1D(ms["cscRechitClusterSize"], ";CSC Rechit Size;count", bins["cscSize"]))
        # hhs["cscTime"].append(H1D(ms["cscRechitClusterTimeWeighted"], ";CSC Time [ns];count", bins["cscTime"]))
        # hhs["cscTimeSpread"].append(
        #     H1D(ms["cscRechitClusterTimeSpreadWeightedAll"], ";CSC Time Spread [ns];count", bins["cscTimeSpread"])
        # )
        # hhs["cscNSt"].append(H1D(ms["cscRechitClusterNStation10"], ";N Stations per CSC Cluster;count", bins["cscNSt"]))
        # hhs["cscASt"].append(
        #     H1D(aabs(ms["cscRechitClusterAvgStation10"]), ";Avg Station per CSC Cluster;count", bins["cscASt"])
        # )
        # hhs["cscMSt"].append(
        #     H1D(aabs(ms["cscRechitClusterMaxStation"]), ";Max Station per CSC Cluster;count", bins["cscMSt"])
        # )
        # hhs["cscMetDPhi"].append(
        #     H1D(ms["cscRechitClusterMet_dPhi"], ";#Delta#phi(CSC,E_{T}^{miss});count", bins["cscMetDPhi"])
        # )
        # hhs["cscEta"].append(H1D(aabs(ms["cscRechitClusterEta"]), ";CSC #eta;count", bins["cscEta"]))
        # hhs["cscPhi"].append(H1D(ms["cscRechitClusterPhi"], ";CSC #phi;count", bins["cscPhi"]))

        # hhs["dtN"].append(H1D(ms["nDtRechitClusters"], ";N DT Clusters;count", bins["dtN"]))
        # hhs["dtSize"].append(H1D(ms["dtRechitClusterSize"], ";DT Rechit Size;count", bins["dtSize"]))
        # hhs["dtTime"].append(H1D(ms["dtRechitCluster_match_RPCBx_dPhi0p5"], ";DT Time [ns];count", bins["dtTime"]))
        # hhs["dtNSt"].append(H1D(ms["dtRechitClusterNStation10"], ";N Stations per DT Cluster;count", bins["dtNSt"]))
        # hhs["dtASt"].append(
        #     H1D(aabs(ms["dtRechitClusterAvgStation10"]), ";Avg Station per DT Cluster;count", bins["dtASt"])
        # )
        # hhs["dtMSt"].append(
        #     H1D(aabs(ms["dtRechitClusterMaxStation"]), ";Max Station per DT Cluster;count", bins["dtMSt"])
        # )
        # hhs["dtMetDPhi"].append(
        #     H1D(ms["dtRechitClusterMet_dPhi"], ";#Delta#phi(DT,E_{T}^{miss});count", bins["dtMetDPhi"])
        # )
        # hhs["dtEta"].append(H1D(aabs(ms["dtRechitClusterEta"]), ";DT #eta;count", bins["dtEta"]))
        # hhs["dtPhi"].append(H1D(ms["dtRechitClusterPhi"], ";DT #phi;count", bins["dtPhi"]))

        # hhs["gllpN"].append(H1D(ms["nGLLP"], ";N gLLPs;count", bins["gllpN"]))
        # hhs["gllpE"].append(H1D(ms["gLLP_e"], ";gLLP E [GeV];count", bins["gllpE"]))
        # hhs["gllpPt"].append(H1D(ms["gLLP_pt"], ";gLLP P_{T} [GeV];count", bins["gllpPt"]))
        # hhs["gllpEta"].append(H1D(ms["gLLP_eta"], ";gLLP #eta;count", bins["gllpEta"]))
        # hhs["gllpPhi"].append(H1D(ms["gLLP_phi"], ";gLLP #phi;count", bins["gllpPhi"]))

        # hhs["jetN"].append(H1D(ms["nJets"], ";N Jets;count", bins["jetN"]))
        # hhs["jetE"].append(H1D(ms["jetE"], ";Jet E [GeV];count", bins["jetE"]))
        # hhs["jetPt"].append(H1D(ms["jetPt"], ";Jet P_{T} [GeV];count", bins["jetPt"]))
        # hhs["jetEta"].append(H1D(ms["jetEta"], ";Jet #eta;count", bins["jetEta"]))
        # hhs["jetPhi"].append(H1D(ms["jetPhi"], ";Jet #phi;count", bins["jetPhi"]))

        # hhs["lepN"].append(H1D(ms["nLeptons"], ";N Leptons;count", bins["lepN"]))
        # hhs["lepE"].append(H1D(ms["lepE"], ";Lepton E [GeV];count", bins["lepE"]))
        # hhs["lepPt"].append(H1D(ms["lepPt"], ";Lepton P_{T} [GeV];count", bins["lepPt"]))
        # hhs["lepEta"].append(H1D(ms["lepEta"], ";Lepton #eta;count", bins["lepEta"]))
        # hhs["lepPhi"].append(H1D(ms["lepPhi"], ";Lepton #phi;count", bins["lepPhi"]))

        # # MET vs (csc, dt) MET dPhi
        # hhs["2D_met_cscMetDPhi"].append(
        #     H2D(
        #         ms["met"],
        #         ms["cscRechitClusterMet_dPhi"],
        #         ";E_{T}^{miss} [GeV];#Delta#phi(csc, #phi_{T}^{miss});",
        #         bins["met_cscMetDPhi"],
        #     )
        # )
        # hhs["2D_met_dtMetDPhi"].append(
        #     H2D(
        #         ms["met"],
        #         ms["dtRechitClusterMet_dPhi"],
        #         ";E_{T}^{miss} [GeV];#Delta#phi(dt, #phi_{T}^{miss});count",
        #         bins["met_dtMetDPhi"],
        #     )
        # )

        # ###################################
        # ## Single Pair CSC-DT Kinematics ##
        # ###################################

        # for second in ('csc', 'dt'):
        #     print(f"Generating CSC-{second.upper()} kinematics plots")

        #     if second == 'csc':
        #         sel_hlt = ms["HLT_CscCluster_Loose"]
        #         sel_2csc = ms["nCscRechitClusters"] == 2

        #         if not isMC:
        #             lpair = "(CSC,CSC_{OOT})"
        #             sel = land(sel_hlt, sel_2csc)
        #         else:
        #             lpair = "(CSC,CSC)"
        #             sel = sel_2csc

        #         # dEta, dPhi, dR
        #         met = ms["met"][sel]
        #         dEta = ms["cscRechitClusterEta"][sel,0] - ms["cscRechitClusterEta"][sel,1]
        #         dPhi = ms["cscRechitClusterPhi"][sel,0] - ms["cscRechitClusterPhi"][sel,1]
        #         dPhi = dPhi - (dPhi > pi) * 2 * pi
        #         dR = np.sqrt(dEta * dEta + dPhi * dPhi)
        #     if second == 'dt':
        #         sel_hlt = ms["HLT_L1CSCCluster_DTCluster50"]
        #         sel_1csc1dt = land(ms["nCscRechitClusters"] == 1, ms["nDtRechitClusters"] == 1)

        #         if not isMC:
        #             lpair = "(CSC,DT_{OOT})"
        #             sel = land(sel_hlt, sel_1csc1dt)
        #         else:
        #             lpair = "(CSC,DT)"
        #             sel = sel_1csc1dt

        #         # dEta, dPhi, dR
        #         met = ms["met"][sel]
        #         dEta = ms["cscRechitClusterEta"][sel] - ms["dtRechitClusterEta"][sel]
        #         dPhi = ms["cscRechitClusterPhi"][sel] - ms["dtRechitClusterPhi"][sel]
        #         dPhi = dPhi - (dPhi > pi) * 2 * pi
        #         dR = np.sqrt(dEta * dEta + dPhi * dPhi)

        #     hhs["csc"+second+"_dEta"].append(H1D(dEta, ";#Delta#eta"+lpair+";count", bins["dEta"]))
        #     hhs["csc"+second+"_dPhi"].append(H1D(dPhi, ";#Delta#phi"+lpair+";count", bins["dPhi"]))
        #     hhs["csc"+second+"_dR"].append(H1D(dR, ";#DeltaR"+lpair+";count", bins["dR"]))

        #     # d__ vs d__
        #     hhs["csc"+second+"_2D_dPhi_dEta"].append(
        #         H2D(dPhi, dEta, ";#Delta#phi"+lpair+";#Delta#eta"+lpair+";count", bins["dPhi_dEta"])
        #     )
        #     hhs["csc"+second+"_2D_dPhi_dR"].append(
        #         H2D(dPhi, dR, ";#Delta#phi"+lpair+";#DeltaR"+lpair+";count", bins["dPhi_dR"])
        #     )
        #     hhs["csc"+second+"_2D_dEta_dR"].append(
        #         H2D(dEta, dR, ";#Delta#eta"+lpair+";#DeltaR"+lpair+";count", bins["dEta_dR"])
        #     )

        #     #
        #     hhs["csc"+second+"_2D_dPhi_met"].append(
        #         H2D(dPhi, met, ";#Delta#phi"+lpair+";E_{T}^{miss} [GeV];count", bins["dPhi_met"])
        #     )
        #     hhs["csc"+second+"_2D_dEta_met"].append(
        #         H2D(dEta, met, ";#Delta#eta"+lpair+";E_{T}^{miss} [GeV];count", bins["dEta_met"])
        #     )
        #     hhs["csc"+second+"_2D_dR_met"].append(
        #         H2D(dR, met, ";#DeltaR"+lpair+";E_{T}^{miss} [GeV];count", bins["dPhi_met"])
        #     )

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

    #     #####################################################
    #     ## Clear memory because im having memory issues :( ##
    #     #####################################################
    #     del ms

    # fout_root.Close()
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
