"""For analyzing multiple files at once (e.g. comparing clustering algorithms)"""

import os
import pathlib
from collections import defaultdict

import ROOT as rt
import numpy as np
import sklearn as skl

from src import CMS_lumi, tdrstyle
from src.muon_system import MuonSystem
from src.helper_functions import lnot, land, lor, lxor, asum, aabs

from src.muon_system import get_lat_leg, H1D, H2D, multi_plot, make_cluster_eff_1D, draw_csc_z_boxes
from src.helper_functions import canvas

pi = rt.TMath.Pi()
bins = {
    "met": (100, 0, 400),
    "metPhi": (100, -pi, pi),
    #
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

bins2D = {
    "met_cscMetDPhi": bins["met"] + bins["cscMetDPhi"],
    "met_dtMetDPhi": bins["met"] + bins["dtMetDPhi"],
    "dPhi_dEta": bins["dPhi"] + bins["dEta"],
    "dPhi_dR": bins["dPhi"] + bins["dR"],
    "dEta_dR": bins["dEta"] + bins["dR"],
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
    out_dir = cur_dir + "/reports/weekly/apr6/"
    data_dir = cur_dir + "/data/raw/"
    # out_dir = "/home/psimmerl/Documents/CMS/LLP/reports/weekly/apr3/"
    # data_dir = "/home/psimmerl/Documents/CMS/LLP/data/raw/"
    ending = ".png"

    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)  # make out directory if it doesn't exist
    print(f"Using output directory '{out_dir}'")

    mc_db_0p4 = data_dir + "ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v4.root"
    mc_ca_0p4 = data_dir + "ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v5.root"
    mc_ca_0p5 = data_dir + "ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v7.root"
    mc_ca_0p6 = data_dir + "ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v6.root"
    mc_ca_0p8 = data_dir + "ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v8.root"
    mc_ca_1p0 = data_dir + "ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v9.root"

    r3_db_0p4 = data_dir + "DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi_v4.root"
    r3_ca_0p4 = data_dir + "DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi_v5.root"
    r3_ca_0p5 = data_dir + "DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi_v7.root"
    r3_ca_0p6 = data_dir + "DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi_v6.root"
    r3_ca_0p8 = data_dir + "DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi_v8.root"
    r3_ca_1p0 = data_dir + "DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi_v9.root"

    files_mc = [mc_db_0p4, mc_ca_0p4, mc_ca_0p5, mc_ca_0p6, mc_ca_0p8, mc_ca_1p0]
    files_r3 = [r3_db_0p4, r3_ca_0p4, r3_ca_0p5, r3_ca_0p6, r3_ca_0p8, r3_ca_1p0]

    labels_mc = [
        "Signal: DBSCAN 0.4",
        "Signal: CA 0.4",
        "Signal: CA 0.5",
        "Signal: CA 0.6",
        "Signal: CA 0.8",
        "Signal: CA 1.0",
    ]
    labels_r3 = [
        "Data: DBSCAN 0.4",
        "Data: CA 0.4",
        "Data: CA 0.5",
        "Data: CA 0.6",
        "Data: CA 0.8",
        "Data: CA 1.0",
    ]

    stats_mc = [
        "mc_db_0p4",
        "mc_ca_0p4",
        "mc_ca_0p5",
        "mc_ca_0p6",
        "mc_ca_0p8",
        "mc_ca_1p0",
    ]
    stats_r3 = [
        "r3_db_0p4",
        "r3_ca_0p4",
        "r3_ca_0p5",
        "r3_ca_0p6",
        "r3_ca_0p8",
        "r3_ca_1p0",
    ]
    # muon_systems_mc = [MuonSystem(ff, isMC=True) for ff in files_mc]
    # muon_systems_r3 = [MuonSystem(ff, isMC=False) for ff in files_r3]

    # muon_systems = muon_systems_mc + muon_systems_r3
    files = files_mc + files_r3
    labels = labels_mc + labels_r3
    stats = stats_mc + stats_r3

    hhs = defaultdict(list)

    for ff, label, stat in zip(files, labels, stats):
        print("\n", stat)#ms.file_name)
        isMC = "mc" in stat
        ms = MuonSystem(ff, isMC=isMC)
        print(f"Loaded {len(ms['met']):,} events.")
        
        #=======================================#
        # Making time plots before cuts on time #
        #=======================================#

        hhs["cscAllTime"].append(H1D(ms["cscRechitClusterTimeWeighted"], ";CSC Time [ns];count", bins["cscTime"]))
        hhs["cscAllTimeSpread"].append(
            H1D(ms["cscRechitClusterTimeSpreadWeightedAll"], ";CSC Time Spread [ns];count", bins["cscTimeSpread"])
        )
        hhs["dtAllTime"].append(
            H1D(ms["dtRechitCluster_match_RPCBx_dPhi0p5"], ";DT Time [ns];count", bins["dtTime"])
        )        

        #########################
        ## BLIND OR MATCH DATA ##
        #########################
        if isMC:
            print("!!!!!!!!!!!!!!!!!!!")
            print("!! MATCHING DATA !!")
            print("!!!!!!!!!!!!!!!!!!!")
            idx_csc_match, idx_dt_match = ms.match_cut()
            ms.apply_cut(idx_csc_match, system='csc')
            ms.apply_cut(idx_dt_match, system='dt')

        elif isBlind:
            print("!!!!!!!!!!!!!!!!!!!")
            print("!! BLINDING DATA !!")
            print("!!!!!!!!!!!!!!!!!!!")

            #sel_csc_oot = aabs(ms["cscRechitClusterTimeWeighted"]) > 50
            sel_csc_oot = ms["cscRechitClusterTimeWeighted"] < -12.5
            sel_dt_oot = ms["dtRechitCluster_match_RPCBx_dPhi0p5"] != 0

            sel_bkg = lor(
                land(ms["nCscRechitClusters"] == 1, asum(sel_csc_oot) == 1),
                land(
                    ms["nCscRechitClusters"] == 2, asum(sel_csc_oot) == 1
                ),  # NOTE: change to req only 2nd clstr oot
                land(ms["nDtRechitClusters"] == 1, asum(sel_dt_oot) == 1),
                land(ms["nDtRechitClusters"] == 2, asum(sel_dt_oot) == 1),
                land(ms["nCscRechitClusters"] == 1, ms["nDtRechitClusters"] == 1, asum(sel_dt_oot) == 1),
            )

            ms.apply_cut(sel_bkg)

        #=================================#
        # NOTE: Using convention "x vs y" #
        #=================================#

        #################################
        ## System & Cluster Kinematics ##
        #################################
        print("Generating system and cluster kinematics plots")

        # System level
        hhs["met"].append(H1D(ms["met"], ";E_{T}^{miss} [GeV];count", bins["met"]))
        hhs["metPhi"].append(H1D(ms["metPhi"], ";#phi[E_{T}^{miss}];count", bins["metPhi"]))

        # kinematics: (csc, dt, gllp, jet, lep)
        # N, Size, Time, TimeSpread, NSt, ASt, MSt, MetDPhi, E, Pt, Eta, Phi
        hhs["cscN"].append(H1D(ms["nCscRechitClusters"], ";N CSC Clusters;count", bins["cscN"]))
        hhs["cscSize"].append(H1D(ms["cscRechitClusterSize"], ";CSC Rechit Size;count", bins["cscSize"]))
        hhs["cscTime"].append(H1D(ms["cscRechitClusterTimeWeighted"], ";CSC Time [ns];count", bins["cscTime"]))
        hhs["cscTimeSpread"].append(
            H1D(ms["cscRechitClusterTimeSpreadWeightedAll"], ";CSC Time Spread [ns];count", bins["cscTimeSpread"])
        )
        hhs["cscNSt"].append(H1D(ms["cscRechitClusterNStation10"], ";N Stations per CSC Cluster;count", bins["cscNSt"]))
        hhs["cscASt"].append(
            H1D(aabs(ms["cscRechitClusterAvgStation10"]), ";Avg Station per CSC Cluster;count", bins["cscASt"])
        )
        hhs["cscMSt"].append(
            H1D(aabs(ms["cscRechitClusterMaxStation"]), ";Max Station per CSC Cluster;count", bins["cscMSt"])
        )
        hhs["cscMetDPhi"].append(
            H1D(ms["cscRechitClusterMet_dPhi"], ";#Delta#phi(CSC,E_{T}^{miss});count", bins["cscMetDPhi"])
        )
        hhs["cscEta"].append(H1D(aabs(ms["cscRechitClusterEta"]), ";CSC #eta;count", bins["cscEta"]))
        hhs["cscPhi"].append(H1D(ms["cscRechitClusterPhi"], ";CSC #phi;count", bins["cscPhi"]))

        hhs["dtN"].append(H1D(ms["nDtRechitClusters"], ";N DT Clusters;count", bins["dtN"]))
        hhs["dtSize"].append(H1D(ms["dtRechitClusterSize"], ";DT Rechit Size;count", bins["dtSize"]))
        hhs["dtTime"].append(H1D(ms["dtRechitCluster_match_RPCBx_dPhi0p5"], ";DT Time [ns];count", bins["dtTime"]))
        hhs["dtNSt"].append(H1D(ms["dtRechitClusterNStation10"], ";N Stations per DT Cluster;count", bins["dtNSt"]))
        hhs["dtASt"].append(
            H1D(aabs(ms["dtRechitClusterAvgStation10"]), ";Avg Station per DT Cluster;count", bins["dtASt"])
        )
        hhs["dtMSt"].append(
            H1D(aabs(ms["dtRechitClusterMaxStation"]), ";Max Station per DT Cluster;count", bins["dtMSt"])
        )
        hhs["dtMetDPhi"].append(
            H1D(ms["dtRechitClusterMet_dPhi"], ";#Delta#phi(DT,E_{T}^{miss});count", bins["dtMetDPhi"])
        )
        hhs["dtEta"].append(H1D(aabs(ms["dtRechitClusterEta"]), ";DT #eta;count", bins["dtEta"]))
        hhs["dtPhi"].append(H1D(ms["dtRechitClusterPhi"], ";DT #phi;count", bins["dtPhi"]))

        hhs["gllpN"].append(H1D(ms["nGLLP"], ";N gLLPs;count", bins["gllpN"]))
        hhs["gllpE"].append(H1D(ms["gLLP_e"], ";gLLP E [GeV];count", bins["gllpE"]))
        hhs["gllpPt"].append(H1D(ms["gLLP_pt"], ";gLLP P_{T} [GeV];count", bins["gllpPt"]))
        hhs["gllpEta"].append(H1D(ms["gLLP_eta"], ";gLLP #eta;count", bins["gllpEta"]))
        hhs["gllpPhi"].append(H1D(ms["gLLP_phi"], ";gLLP #phi;count", bins["gllpPhi"]))

        hhs["jetN"].append(H1D(ms["nJets"], ";N Jets;count", bins["jetN"]))
        hhs["jetE"].append(H1D(ms["jetE"], ";Jet E [GeV];count", bins["jetE"]))
        hhs["jetPt"].append(H1D(ms["jetPt"], ";Jet P_{T} [GeV];count", bins["jetPt"]))
        hhs["jetEta"].append(H1D(ms["jetEta"], ";Jet #eta;count", bins["jetEta"]))
        hhs["jetPhi"].append(H1D(ms["jetPhi"], ";Jet #phi;count", bins["jetPhi"]))

        hhs["lepN"].append(H1D(ms["nLeptons"], ";N Leptons;count", bins["lepN"]))
        hhs["lepE"].append(H1D(ms["lepE"], ";Lepton E [GeV];count", bins["lepE"]))
        hhs["lepPt"].append(H1D(ms["lepPt"], ";Lepton P_{T} [GeV];count", bins["lepPt"]))
        hhs["lepEta"].append(H1D(ms["lepEta"], ";Lepton #eta;count", bins["lepEta"]))
        hhs["lepPhi"].append(H1D(ms["lepPhi"], ";Lepton #phi;count", bins["lepPhi"]))

        # MET vs (csc, dt) MET dPhi
        hhs["2D_met_cscMetDPhi_" + stat].append(
            H2D(
                ms["met"],
                ms["cscRechitClusterMet_dPhi"],
                ";E_{T}^{miss} [GeV];#Delta#phi(csc, #phi_{T}^{miss});",
                bins["met_cscMetDPhi"],
            )
        )
        hhs["2D_met_dtMetDPhi_" + stat].append(
            H2D(
                ms["met"],
                ms["dtRechitClusterMet_dPhi"],
                ";E_{T}^{miss} [GeV];#Delta#phi(dt, #phi_{T}^{miss});count",
                bins["met_dtMetDPhi"],
            )
        )

        ##################################
        # Single Pair CSC-DT Kinematics ##
        ##################################
        print("Generating CSC-DT kinematics plots")
        sel_hlt = ms["HLTDecision"][:, 569]  # HLT_L1CSCCluster_DTCluster50
        sel_1csc1dt_orig = land(ms["nCscRechitClusters"] == 1, ms["nDtRechitClusters"] == 1)
        sel_1csc1dt = land(asum(ms['cscRechitClusterSize']>0) == 1, asum(ms['dtRechitClusterSize']>0) == 1)
        if np.sum(sel_1csc1dt_orig ^ sel_1csc1dt) != 0:
            print('UH OH SOMETHING BAD HAPPENED!')
            exit()
        if not isMC:
            sel = land(sel_hlt, sel_1csc1dt)
        else:
            sel = sel_1csc1dt
        
        # dEta, dPhi, dR
        dEta = ms["cscRechitClusterEta"][sel] - ms["dtRechitClusterEta"][sel]
        dPhi = ms["cscRechitClusterPhi"][sel] - ms["dtRechitClusterPhi"][sel]
        dPhi = dPhi - (dPhi > pi) * 2 * pi
        dR = np.sqrt(dEta * dEta + dPhi * dPhi)

        hhs["dEta"].append(H1D(dEta, ";#Delta#eta(CSC,DT_{OOT});count", bins["dEta"]))
        hhs["dPhi"].append(H1D(dPhi, ";#Delta#phi(CSC,DT_{OOT});count", bins["dPhi"]))
        hhs["dR"].append(H1D(dR, ";#DeltaR(CSC,DT_{OOT});count", bins["dR"]))

        # d__ vs d__
        hhs["2D_dPhi_dEta_"+stat].append(
            H2D(dPhi, dEta, ";#Delta#phi(CSC,DT_{OOT});#Delta#eta(CSC,DT_{OOT});count", bins["dPhi_dEta"])
        )
        hhs["2D_dPhi_dR_"+stat].append(H2D(dPhi, dR, ";#Delta#phi(CSC,DT_{OOT});#DeltaR(CSC,DT_{OOT});count", bins["dPhi_dR"]))
        hhs["2D_dEta_dR_"+stat].append(H2D(dEta, dR, ";#Delta#eta(CSC,DT_{OOT});#DeltaR(CSC,DT_{OOT});count", bins["dEta_dR"]))

        print("Saving plots to disk.")
        for k, hs in hhs.items():
            c = canvas(1, 1)
            if "2D" in k:
                c.cd(1).SetLogz()
                hs[0].Draw("colz")
            else:
                _ll = multi_plot(hs, labels, ymin=0, norm=True)
            c.Print(out_dir + k + ending)
    
    #########################
    ## 1D Efficiency Plots ##
    #########################
    print("Generating efficiency plots")

    for isCut in [False, True]:
        c = canvas(2, 2)
        for i, det in enumerate(["csc", "dt"]):
            for j, xl in enumerate(["z", "r"]):
                c.cd(2 * i + j + 1)
                _hhs = []
                for ff in files_mc:
                    ms = MuonSystem(ff, isMC=True)
                    _hhs.append(make_cluster_eff_1D(ms, det, xl, cuts=isCut))
                if isCut:
                    legxy = (0.3, 0.7, 0.5, 0.85)
                else:
                    legxy = (0.5, 0.15, 0.7, 0.4)
                _ll = multi_plot(_hhs, labels_mc, legxy=legxy)
                for h in _hhs:
                    h.SetMinimum(0.0)
                    h.SetMaximum(1.0)

                if det == "csc" and xl == "z":
                    boxes = draw_csc_z_boxes(_hhs[-1])
                    gc.extend(boxes)
                gc.extend(_hhs + list(_ll))

        c.Print(out_dir + "eff_compare_DBSCAN-CA" + ("_CUT" if isCut else "") + ending)

    print("Finished!")
