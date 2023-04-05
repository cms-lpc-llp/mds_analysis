"""For analyzing multiple files at once (e.g. comparing clustering algorithms)"""

import pathlib
import ROOT as rt
from src import CMS_lumi, tdrstyle
from src.muon_system import MuonSystem

from src.muon_system import get_lat_leg, H1D, multi_plot, make_cluster_eff_1D, draw_csc_z_boxes
from src.helper_functions import canvas


if __name__ == "__main__":
    gc = []
    rt.gROOT.SetBatch()

    a = tdrstyle.setTDRStyle()
    # CMS_lumi.writeExtraText = 0
    rt.gStyle.SetOptFit(0)  # 1011)

    out_dir = "/home/psimmerl/Documents/CMS/LLP/reports/weekly/apr3/"
    data_dir = "/home/psimmerl/Documents/CMS/LLP/data/raw/"
    ending = ".png"

    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)  # make out directory if it doesn't exist
    print(f"Using output directory '{out_dir}'")

    mc_db_0p4 = data_dir + "ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v4.root"
    mc_ca_0p4 = data_dir + "ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v5.root"
    mc_ca_0p5 = data_dir + "ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v7.root"
    mc_ca_0p6 = data_dir + "ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v6.root"
    mc_ca_0p8 = data_dir + "ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v8.root"
    mc_ca_1p0 = data_dir + "ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v9.root"
    run3_file = data_dir + "DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi.root"

    files_mc = [mc_db_0p4, mc_ca_0p4, mc_ca_0p5, mc_ca_0p6, mc_ca_0p8, mc_ca_1p0]
    labels = [
        "DBSCAN 0.4",
        "CA 0.4",
        "CA 0.5",
        "CA 0.6",
        "CA 0.8",
        "CA 1.0",
    ]
    muon_systems_mc = [MuonSystem(ff, isMC=True) for ff in files_mc]
    muon_systems = muon_systems_mc + [MuonSystem(run3_file, isMC=False)]

    #########################
    ## 1D Efficiency Plots ##
    #########################

    for isCut in [False, True]:
        c = canvas(2, 2)
        for i, det in enumerate(["csc", "dt"]):
            for j, xl in enumerate(["z", "r"]):
                c.cd(2 * i + j + 1)
                hhs = []
                for ms in muon_systems_mc:
                    print(ms.file_name)
                    hhs.append(make_cluster_eff_1D(ms, det, xl, cuts=isCut))
                if isCut:
                    legxy = (0.3, 0.7, 0.5, 0.85)
                else:
                    legxy = (0.5, 0.25, 0.7, 0.4)
                _ll = multi_plot(hhs, labels, legxy=legxy)
                for h in hhs:
                    h.SetMinimum(0.0)
                    h.SetMaximum(1.0)

                if det == "csc" and xl == "z":
                    boxes = draw_csc_z_boxes(hhs[-1])
                    gc.extend(boxes)
                gc.extend(hhs + list(_ll))

        c.Print(out_dir + "eff_compare_DBSCAN-CA" + ("_CUT" if isCut else "") + ending)

    #################################
    ## System & Cluster Kinematics ##
    #################################

    for ms in muon_systems_mc:
        # MET
        print(ms.file_name)
        hhs["met"].append(H1D(ms["met"], ";E_{T}^{miss} [GeV];count", (100, 0, 2000)))

        hhs["csc_dPhiMET"].append(
            H1D(ms["cscRechitCluster_dPhi"], f";#Delta#phi(csc," + "#phi_{T}^{miss});count", (100, 0, 4))
        )
        hhs["dt_dPhiMET"].append(
            H1D(ms["dtRechitCluster_dPhi"], f";#Delta#phi(dt," + "#phi_{T}^{miss});count", (100, 0, 4))
        )

    # (both, csc, dt) dPhiMET
    # MET vs (both, csc, dt) dPhiMET

    ##################################
    # Single Pair CSC-DT Kinematics ##
    ##################################

    # dEta, dPhi, dPt, dR
    # d__ vs d__
