"""train_bdt.py
**Note - turned into more of a 'main.py'**

To run:
For all events
    'python train_bdt.py'
For a subset
    'python train_bdt.py {N_EVENTS}'
"""

__author__ = "Paul Simmerling"
__email__ = "psimmerl@caltech.edu"
__credits__ = [
    "Paul Simmerling",
    "Christina Wang",
    "Lisa Benato",
    "Si Xie",
    "Cristian Pena",
    "Martin Kwok",
    "Pedro Fernandez Manteca",
    "Maria Spiropulu",
]

import sys
import pathlib
from os import uname
import json

import numpy as np
import ROOT as rt
from math import ceil, floor

from src.muon_system import MuonSystemAwkward
from src import CMS_lumi, tdrstyle
from src.helper_functions import alert, Table  # , canvas
from src.histo_utilities import create_TH1D, create_TH2D, create_TGraph, std_color_list

from sklearn.metrics import roc_auc_score  # , roc_curve
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import GradientBoostingClassifier  # , RandomForestClassifier

# TODO: make a yaml to hold the parameter
# **************************** #
# *        Parameters        * #
# **************************** #
SAVE_STAT = "TEST"  # "DTOOT"
BLIND_TYPE = ["DTOOT", "SIZE"][1]
# SAVE_STAT = "BLINDSR"
# BLIND_TYPE = "BLINDSR"

N_EVENTS = -1
TEST_SIZE = 0.8

ABCD_DPHI = 2.75
ABCD_DTSIZE = 90
ABCD_INVERT_DPHI = False
# **************************** #
CUTS = [
    "match",
    #!!!!!!!#
    "HLT",  # TODO: Load HLTDecision without overflowing memory
    #!!!!!!!#
    "CSC&DT>0",
    "L1",
    "CSCIT",
    "DTIT",
    # "MET",
    # "ME11/12",
    "MB1",
    "JET",
    # "MUON",
    # "BDT",
    "HALO",
    # "CSCSIZE",
    "DTSTN",
    "1CSC1DT",
    "BLINDSR",
    # "DR",
    "DPHI",
]

TRAIN_BDT = False
SECTIONS = [
    "EFFICIENCY",
    "CSC",
    "DT",
    "EVENT",
    "2D",
    "ABCD",
    "ROCS",
]
# **************************** #
#! Paths are absolute so I know *exactly* what is being read/written
OUT_DIR = "reports/weekly/2023-09-21c"
T2_OUT_DIR = "/storage/af/user/psimmerl/LLP/mdc_analysis"  # os.getcwd()
LOCAL_OUT_DIR = "/home/psimmerl/LLP/mdc_analysis"  # os.getcwd()
# **************************** #
DATA_VERSION = "6"
LUMI = 23.02 * 1000

T2_DATA_DIR = "/storage/cms/store/user/christiw/displacedJetMuonAnalyzer/Run3/V1p19"
LOCAL_DATA_DIR = "/home/psimmerl/LLP/mdc_analysis/data/raw"  # os.getcwd() + "/data/raw"
DATA_DIR = "TIER2" if "caltech" in uname()[1] else "LOCAL"

FN_MC = "ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted"
FN_R3 = "DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi"

TN_MC = "MuonSystem"
TN_R3 = "MuonSystem_HLT569" if "HLT" in CUTS else "MuonSystem"

# with open("default_bins.json", "r") as f: #TODO: this should be absolute for consitency
#     DEFAULT_BINS = json.load(f)
# **************************** #
VERBOSE = True
ROOT_BATCH = True
ROOT_ERROR_LEVEL = 1001  # rt.kInfo + 1
BOT_MARGIN, TOP_MARGIN = 0.025, 0.1

gc = []
# **************************** #


def create_hists(
    values,
    bins,
    names,
    titles,
    weights=1,
    colors=None,
    styles=None,
    log=None,
    norm=False,
    share_scale=True,
    canvas=None,
):
    if not isinstance(names, (list, tuple)):
        (
            values,
            names,
        ) = [
            values
        ], [names]

    if len(titles) == 3:  # values[0].shape[1] == 2
        hist_type = "2D"
    else:
        hist_type = "1D"
        # TODO: flatten?

    if colors is None:
        colors = std_color_list

    if not isinstance(weights, (list, tuple)):
        weights = [weights] * len(names)
    if not isinstance(colors, (list, tuple)):
        colors = [colors] * len(names)
    if not isinstance(styles, (list, tuple)):
        styles = [styles] * len(names)

    logx, logy, logz = False, False, False
    if isinstance(log, (list, tuple)):
        logx = log[0]
        logy = log[1]
        if len(log) > 2:
            logz = log[2]
            log = logz
        else:
            log = logy
    else:
        if hist_type == "1D":
            logy = log
        else:
            logz = log

    hists = []

    RAND_NUM = f".{SAVE_STAT}.{np.random.randint(999999999)}"
    if canvas is None:
        if hist_type == "2D":
            canvas = rt.TCanvas("c" + RAND_NUM, "c" + RAND_NUM, len(names) * 800, 800)
            canvas.Divide(len(names), 1)
            for ic in range(len(names)):
                canvas.cd(ic + 1).SetGrid()
                canvas.cd(ic + 1).SetRightMargin(0.04)
        else:
            canvas = rt.TCanvas("c" + RAND_NUM, "c" + RAND_NUM, 800, 800)
            canvas.SetGrid()
            canvas.SetRightMargin(0.04)

    lat = rt.TLatex()
    lat.SetTextAlign(21)  # left,cent,right=1,2,3 , bot,cent,top=1,2,3
    lat.SetTextSize(0.06)

    legend = rt.TLegend(0.17, 0.91 - 0.04 * len(names), 0.21 + 0.016 * (max([len(n) for n in names]) + 11), 0.94)
    legend.SetTextFont(62)
    legend.SetTextSize(0.03)
    legend.SetBorderSize(0)
    legend.SetFillColorAlpha(1, 0.2)
    legend.SetEntrySeparation(0.01)

    hmin, hmax = 1e9, -1

    for n, v, w, c, s in zip(names, values, weights, colors, styles):
        if hist_type == "1D":
            hist = create_TH1D(v, name=n + RAND_NUM, axis_title=titles, binning=bins, weights=w)
        else:
            hist = create_TH2D(v, name=n + RAND_NUM, axis_title=titles, binning=bins, weights=w)
        hist.SetLineColor(c)
        hist.SetMarkerColor(c)
        # hist.SetFillColorAlpha(c, 0.3)
        if s is not None:
            hist.SetLineStyle(c)

        if norm and hist.Integral() > 0:
            hist.Scale(1.0 / hist.Integral())

        hmax = max(hmax, hist.GetMaximum())
        hmin = min(hmin, hist.GetMinimum(0))
        hists.append(hist)

    # print(log, log is None, hmax / hmin)
    if share_scale:
        if log or (log is None and hmax / hmin > 1000):  # clean this up
            hmin, hmax = np.log10(hmin if hmin == 1 else 1e-1), np.log10(hmax)
            hmin, hmax = hmin - BOT_MARGIN * (hmax - hmin), hmax + TOP_MARGIN * (hmax - hmin)
            hmin, hmax = 10 ** (hmin), 10 ** (hmax)
            log = True
        else:
            hmin, hmax = hmin - BOT_MARGIN * (hmax - hmin), hmax + TOP_MARGIN * (hmax - hmin)
            hmin = 0 if hmin < 50 else hmin
            log = False

    for ih, hist in enumerate(hists):
        if share_scale:
            hist.SetMaximum(hmax)
            hist.SetMinimum(hmin)
        if hist_type == "2D":
            canvas.cd(ih + 1)
            hist.Draw("colz")
            lat.SetTextColor(hist.GetLineColor())
            lat.DrawLatexNDC(0.5, 0.96, hist.GetName().split(".")[0])
        else:
            hist.Draw("hist same")
            le = legend.AddEntry(hist, hist.GetName().split(".")[0], "PE")

    if hist_type == "1D":
        legend.Draw()

    if logx:
        canvas.SetLogx()
    if logy or (log and hist_type == "1D"):
        canvas.SetLogy()
    if logz or (log and hist_type == "2D"):
        canvas.SetLogz()

    if not ROOT_BATCH:
        canvas.Draw()

    return canvas, legend, hists


# ************************************************************ #
# ************************************************************ #
# ************************************************************ #


def main():
    print("")
    alert("Building MuonSystem(s)", form="-", c="g")
    ms_mc = MuonSystemAwkward(FN_MC, name="Signal", tree_name=TN_MC, nev=N_EVENTS, is_mc=True, lumi=LUMI)
    ms_r3 = MuonSystemAwkward(FN_R3, name="Data", tree_name=TN_R3, nev=N_EVENTS, is_mc=False, lumi=LUMI)
    mss = [ms_mc, ms_r3]

    if BLIND_TYPE != "BLINDSR":
        CUTS.remove("BLINDSR")

    if len(sys.argv) > 2:
        if "standard" in sys.argv[2]:
            pass
        elif "neither" in sys.argv[2]:
            CUTS.remove("HALO")
            CUTS.remove("DPHI")

        elif "nodtstn" in sys.argv[2]:
            CUTS.remove("DTSTN")
        elif "dtstn" in sys.argv[2]:
            ms_nocut = MuonSystemAwkward(
                FN_R3, name="NoDtStnCut", tree_name=TN_R3, nev=N_EVENTS, is_mc=False, lumi=LUMI
            )
            mss.append(ms_nocut)

        elif "nohalo" in sys.argv[2]:
            CUTS.remove("HALO")
        elif "halo" in sys.argv[2]:
            ms_nocut = MuonSystemAwkward(FN_R3, name="NoHaloCut", tree_name=TN_R3, nev=N_EVENTS, is_mc=False, lumi=LUMI)
            mss.append(ms_nocut)

        elif "nodphi" in sys.argv[2]:
            CUTS.remove("DPHI")
        elif "dphi" in sys.argv[2]:
            ms_nocut = MuonSystemAwkward(FN_R3, name="NoDPhiCut", tree_name=TN_R3, nev=N_EVENTS, is_mc=False, lumi=LUMI)
            mss.append(ms_nocut)
        elif "both" in sys.argv[2]:
            ms_cut = MuonSystemAwkward(
                FN_R3, name="NoHaloNoDPhiCuts", tree_name=TN_R3, nev=N_EVENTS, is_mc=False, lumi=LUMI
            )
            mss.append(ms_cut)

    # ************************************************************ #

    for ims, ms in enumerate(mss):
        if ms.colors is None:
            ms.colors = [std_color_list[ims], std_color_list[len(mss) + ims]]
        #!!! TURNING MuonSystem CUTS OFF !!!!#
        # cuts are not applied when you get MuonSystem data
        # but they are 'AND'ed with sel_evt, sel_csc, and sel_dt.
        ms.cut = False
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

    # ************************************************************ #

    print("")
    alert("Cutting MuonSystem(s)", form="-", c="g")

    columns = [f"{ms.name} (eff.)" for ms in mss]
    cut_table = Table(columns, "MuonSystem Cuts (req. CSC>0 & DT>0 after each cut)", 1)
    eff_decimals = 4

    raw_counts = [ms.count() for ms in mss]
    cut_table.add_row("read", raw_counts)

    for cut in CUTS:
        counts = [["", ""] for i in range(len(mss))]
        for ims, ms in enumerate(mss):
            was_cut, mc, msn = True, ms.is_mc, ms.name
            match cut:
                case "match":
                    if mc:
                        ms.match_mc("csc,dt")
                    else:
                        was_cut = False
                case "CSC&DT>0":  # just applies has_clusters=True
                    ms.f((ms["nCsc"] > 0) & (ms["nDt"] > 0), has_clusters=False)
                    raw_counts[ims] = ms.count()

                    #!!! TURNING EFFICIENCY ON !!!!#
                    if "EFFICIENCY" in SECTIONS:
                        ms.efficiency_denom = "eff"
                        ms.ms_read["eff_evt"] = ms.ms_read["sel_evt"]
                        ms.ms_read["eff_csc"] = ms.ms_read["sel_csc"]
                        ms.ms_read["eff_dt"] = ms.ms_read["sel_dt"]
                        ms.efficiency = 1
                    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

                case "HLT":
                    # if not mc:
                    #     ms.cut_hlt()
                    pass
                case "ME11/12":
                    ms.f(ms["cscMe11Ratio"] + ms["cscMe12Ratio"] == 0, "csc")
                case "MB1":
                    ms.f(ms["dtNHitStation1"] == 0, "dt")
                case "JET":
                    ms.cut_jet("csc,dt", invert=False)
                case "MUON":
                    ms.cut_muon("csc,dt", invert=False)
                case "L1":
                    ms.cut_l1()
                case "CSCIT":
                    ms.cut_time("csc", invert=False)
                case "DTIT":
                    ms.cut_time("dt", invert=(not mc and BLIND_TYPE == "DTOOT"))
                case "HALO":
                    if "NoHalo" in msn:
                        was_cut = False
                    else:
                        ms.cut_halo(invert="Halo" in msn)
                case "CSCSIZE":
                    ms.f(ms["cscSize"] > 250, "csc")
                # case "DTSIZE":
                #     ms.f(ms["dtSize"] < 200, "dt")
                # case "CSCSTN":
                #     ms.f(ms["cscSize"] > 250, "csc")
                case "DTSTN":
                    if "NoDtStn" in msn:
                        was_cut = False
                    else:
                        ms.f((ms["dtNStation10"] < 3) & ~((ms["dtNStation10"] == 2) & (ms["dtMaxStation"] == 4)), "dt")
                case "1CSC1DT":
                    ms.tag(tags="cscdt")
                case "BLINDSR":
                    if not mc:
                        ms.f((ms["dtSize"] < ABCD_DTSIZE) | (ms["tag_dPhi"] < ABCD_DPHI), "dt")
                    else:
                        was_cut = False
                case "DPHI":
                    if "NoDPhi" in msn:
                        was_cut = False
                    else:
                        ms.f(0.4 < ms["tag_dPhi"], invert="dphi" in msn.lower())
                case default:
                    raise KeyError(f"Cut '{cut}' not recognized!")
            if was_cut:
                counts[ims][0] = ms.count()
                if cut in ("match", "HLT"):  # skip efficiency
                    counts[ims] = counts[ims][0]
                elif "EFFICIENCY" in SECTIONS:
                    counts[ims][1] = f"({ms.efficiency:.{eff_decimals}f})"

        cut_label = "DTIT/OOT" if cut == "DTIT" and BLIND_TYPE == "DTOOT" else cut
        cut_table.add_row(cut_label, counts)
        if cut in ("match",):
            cut_table.add_spacer()
        if VERBOSE:
            print(f"Finished cut '{cut}'")  # , counts)
    cut_table.print()

    # ************************************************************ #

    for ms in mss:
        #!!! TURNING MuonSystem CUTS ON !!!!#
        # cuts are now applied when you MuonSystem get data
        ms.cut = True
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

        #!!! TURNING EFFICIENCY OFF !!!!#
        if ms.efficiency_denom is not None:
            del ms.ms_read[f"{ms.efficiency_denom}_evt"]
            del ms.ms_read[f"{ms.efficiency_denom}_csc"]
            del ms.ms_read[f"{ms.efficiency_denom}_dt"]
            ms.efficiency_denom = None
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
        # Some may *need* 2 tag
        if "EVENT" in SECTIONS:
            ms["cscR-dtR"] = ms["cscR"] - ms["dtR"]
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

    # ************************************************************ #
    # ************************************************************ #
    # ************************************************************ #

    if "CSC" in SECTIONS:
        print("")
        alert("Making CSC Plots", form="-", c="g")

        xls = [
            "cscMuonVetoPt",
            "cscJetVetoPt",
            "cscR",
            "cscZ.abs",
            "cscPhi.abs",
            "cscEta.abs",
            "cscSize",
            "cscNStation10",
            "cscMaxChamber",
            "cscAvgStation10",
            "cscMaxStation",
            "cscMe11Ratio",
            "cscMe12Ratio",
        ]
        bins = [
            [50, 0, 100],
            [50, 0, 100],
            [50, 100, 700],
            [50, 600, 1100],
            [32, 0, np.pi],  # [32, -np.pi, np.pi],
            [50, 0.8, 2.4],
            [50, 0, 2000],
            [6, -0.5, 5.5],
            [50, 0, 50],
            [20, -0.5, 5.5],
            [6, -0.5, 5.5],
            [20, 0, 1],
            [20, 0, 1],
        ]

        names, colors = [ms.name for ms in mss], [ms.colors[0] for ms in mss]
        for xl, _bins in zip(xls, bins):
            weights = [ms["weight"] for ms in mss]
            values = [ms[xl] for ms in mss]
            xl = xl.split(".")[0]

            if "SIZE" in BLIND_TYPE:
                conds = [wt > 0 for wt in weights]
                if "dt" in xl and "csc" not in xl:
                    conds = [c if ms.is_mc else c & (ms["cscSize"][:, 0] < 300) for c, ms in zip(conds, mss)]
                else:
                    conds = [c if ms.is_mc else c & (ms["dtSize"][:, 0] < 100) for c, ms in zip(conds, mss)]

                weights = [w[c] for c, w in zip(conds, weights)]
                values = [v[c] for c, v in zip(conds, values)]

            canvas, legend, hists = create_hists(
                values, _bins, names, [xl, "events"], weights, colors, styles=None, log=None, norm=False
            )

            canvas.Print(f"{OUT_DIR}/{xl}_{SAVE_STAT}.png")
            if VERBOSE:
                print(f"Saved '{xl}_{SAVE_STAT}.png'")

    # **************************** #

    if "DT" in SECTIONS:
        print("")
        alert("Making DT Plots", form="-", c="g")

        xls = [
            "dtJetVetoPt",
            "dtMuonVetoPt",
            "dtR",
            "dtZ.abs",
            "dtPhi.abs",
            "dtEta.abs",
            "dtSize",
            "dt_match_RPCBx_dPhi0p5",
            "dtNStation10",
            "dtAvgStation10",
            "dtMaxStation",
            "dtNHitStation1",
            "dtMb1Ratio",
        ]
        bins = [
            [50, 0, 100],
            [50, 0, 100],
            [50, 450, 800],
            [50, 0, 700],
            [31, 0, np.pi],  # [32, -np.pi, np.pi],
            [50, 0, 1.3],
            [50, 0, 500],
            [13, -6.5, 6.5],
            [5, -0.5, 4.5],
            [20, -0.5, 5],
            [6, -0.5, 5.5],
            [50, 0, 50],
            [20, 0, 1],
        ]

        names, colors = [ms.name for ms in mss], [ms.colors[0] for ms in mss]
        for xl, _bins in zip(xls, bins):
            weights = [ms["weight"] for ms in mss]
            values = [ms[xl] for ms in mss]
            xl = xl.split(".")[0]

            if "SIZE" in BLIND_TYPE:
                conds = [wt > 0 for wt in weights]
                if "dt" in xl and "csc" not in xl:
                    conds = [c if ms.is_mc else c & (ms["cscSize"][:, 0] < 300) for c, ms in zip(conds, mss)]
                else:
                    conds = [c if ms.is_mc else c & (ms["dtSize"][:, 0] < 100) for c, ms in zip(conds, mss)]

                weights = [w[c] for c, w in zip(conds, weights)]
                values = [v[c] for c, v in zip(conds, values)]

            canvas, legend, hists = create_hists(
                values, _bins, names, [xl, "events"], weights, colors, styles=None, log=None, norm=False
            )

            canvas.Print(f"{OUT_DIR}/{xl}_{SAVE_STAT}.png")
            if VERBOSE:
                print(f"Saved '{xl}_{SAVE_STAT}.png'")

    # **************************** #

    if "EVENT" in SECTIONS:
        print("")
        alert("Making Event/2tag Plots", form="-", c="g")

        xls = [
            "tag_dEta",
            "tag_dPhi",
            "tag_dR",
            "nCsc",
            "nDt",
            "met",
            "runNum",
            "cscR-dtR",
        ]
        bins = [
            [32, 0, 3],
            [31, 0, np.pi],
            [32, 0, 4.5],
            [6, -0.5, 5.5],
            [6, -0.5, 5.5],
            [30, 0, 400],
            [2741, 360019, 362760],
            [25, -650, 650],
            # [25, 0, 700],
        ]

        names, colors = [ms.name for ms in mss], [ms.colors[0] for ms in mss]
        for xl, _bins in zip(xls, bins):
            weights = [ms["weight"] for ms in mss]
            values = [ms[xl] for ms in mss]
            xl = xl.split(".")[0]

            if "SIZE" in BLIND_TYPE:
                conds = [wt > 0 for wt in weights]
                if "dt" in xl and "csc" not in xl:
                    conds = [c if ms.is_mc else c & (ms["cscSize"][:, 0] < 300) for c, ms in zip(conds, mss)]
                else:
                    conds = [c if ms.is_mc else c & (ms["dtSize"][:, 0] < 100) for c, ms in zip(conds, mss)]

                weights = [w[c] for c, w in zip(conds, weights)]
                values = [v[c] for c, v in zip(conds, values)]

            canvas, legend, hists = create_hists(
                values, _bins, names, [xl, "events"], weights, colors, styles=None, log=None, norm=False
            )

            canvas.Print(f"{OUT_DIR}/{xl}_{SAVE_STAT}.png")
            if VERBOSE:
                print(f"Saved '{xl}_{SAVE_STAT}.png'")

    # **************************** #

    if "2D" in SECTIONS:
        print("")
        alert("Making 2D Plots", form="-", c="g")

        xyls = [
            ("cscPhi.abs", "cscR"),
            ("cscPhi.abs", "cscEta.abs"),
            ("cscEta.abs", "cscR"),
            #
            ("dtPhi.abs", "dtR"),
            ("dtPhi.abs", "dtEta.abs"),
            ("dtEta.abs", "dtR"),
            #
            ("cscR", "dtR"),
            ("cscPhi.abs", "dtPhi.abs"),
            ("cscEta.abs", "dtEta.abs"),
            #
            ("tag_dPhi", "met"),
            ("tag_dPhi", "tag_dEta"),
            ("tag_dPhi", "tag_dR"),
            ("tag_dEta", "tag_dR"),
            #
            ("tag_dPhi", "cscSize"),
            ("tag_dPhi", "cscPhi.abs"),
            ("tag_dPhi", "cscEta.abs"),
            #
            ("tag_dPhi", "dtSize"),
            ("tag_dPhi", "dtPhi.abs"),
            ("tag_dPhi", "dtEta.abs"),
        ]
        bins = [
            [16, 0, np.pi, 25, 100, 700],
            [16, 0, np.pi, 25, 0.8, 2.4],
            [25, 0.8, 2.4, 25, 100, 700],
            #
            [16, 0, np.pi, 25, 450, 800],
            [16, 0, np.pi, 25, 0, 1.3],
            [25, 0, 1.3, 25, 450, 800],
            #
            [25, 100, 700, 25, 450, 800],
            [16, 0.0, np.pi, 16, 0.0, np.pi],
            [25, 0.8, 2.4, 25, 0.0, 1.3],
            #
            [16, 0, np.pi, 25, 0, 400],
            [16, 0, np.pi, 25, 0, 3],
            [16, 0, np.pi, 25, 0, 5],
            [25, 0, 3, 25, 0, 5],
            #
            [16, 0, np.pi, 25, 200, 1000],
            [16, 0, np.pi, 16, 0, np.pi],
            [16, 0, np.pi, 25, 0.8, 2.4],
            #
            [16, 0, np.pi, 25, 50, 500],
            [16, 0, np.pi, 16, 0, np.pi],
            [16, 0, np.pi, 25, 0.0, 1.3],
        ]

        names, colors = [ms.name for ms in mss], [ms.colors[0] for ms in mss]
        for (xl, yl), _bins in zip(xyls, bins):
            weights = [ms["weight"] for ms in mss]
            values = [np.c_[ms[xl], ms[yl]] for ms in mss]
            xl, yl = xl.split(".")[0], yl.split(".")[0]

            if "SIZE" in BLIND_TYPE:
                conds = [wt > 0 for wt in weights]
                if ("dt" in xl or "dt" in yl) and ("csc" not in xl and "csc" not in yl):
                    conds = [c if ms.is_mc else c & (ms["cscSize"][:, 0] < 300) for c, ms in zip(conds, mss)]
                else:
                    conds = [c if ms.is_mc else c & (ms["dtSize"][:, 0] < 100) for c, ms in zip(conds, mss)]

                weights = [w[c] for c, w in zip(conds, weights)]
                values = [v[c] for c, v in zip(conds, values)]

            canvas, legend, hists = create_hists(
                values, _bins, names, [xl, yl, "events"], weights, colors, styles=None, log=None, norm=False
            )

            canvas.Print(f"{OUT_DIR}/{xl}_{yl}_{SAVE_STAT}.png")
            if VERBOSE:
                print(f"Saved '{xl}_{yl}_{SAVE_STAT}.png'")

    # **************************** #
    # **************************** #
    # **************************** #

    if "ABCD" in SECTIONS:
        print("")
        alert("Making ABCD Plots", form="-", c="g")

        n_abcd = 0
        xyls, bins = [], []
        for th in np.arange(0, 1.0 + 0.1, step=0.1):
            if "DPHI" in CUTS and th < 0.4:
                continue
            xyls.append(("tag_dPhi", "dtSize"))
            bins.append([th, ABCD_DPHI, np.pi, 50, ABCD_DTSIZE, 500, 3, 3])
            n_abcd += 1

        abcds = {ms.name: [] for ms in mss}
        names, colors = [ms.name for ms in mss], [ms.colors[0] for ms in mss]
        for (xl, yl), _bins in zip(xyls, bins):
            weights = [np.ravel(ms["weight"]) for ms in mss]
            values = [np.c_[np.ravel(ms[xl]), np.ravel(ms[yl])] for ms in mss]
            xl, yl = xl.split(".")[0], yl.split(".")[0]

            if "SIZE" in BLIND_TYPE:
                conds = [wt > 0 for wt in weights]
                conds = [
                    c if ms.is_mc else c & ((ms["dtSize"][:, 0] < ABCD_DTSIZE) | (ms["tag_dPhi"] < ABCD_DPHI))
                    for c, ms in zip(conds, mss)
                ]
                weights = [w[c] for c, w in zip(conds, weights)]
                values = [v[c] for c, v in zip(conds, values)]

            canvas, legend, hists = create_hists(
                values, _bins, names, [xl, yl, "events"], weights, colors, styles=None, log=None, norm=False
            )
            for ms, hh in zip(mss, hists):
                abcds[ms.name].append(
                    [
                        _bins[0],
                        hh.GetBinContent(2, 2),  # A
                        hh.GetBinError(2, 2),
                        hh.GetBinContent(2, 1),  # B
                        hh.GetBinError(2, 1),
                        hh.GetBinContent(1, 1),  # C
                        hh.GetBinError(1, 1),
                        hh.GetBinContent(1, 2),  # D
                        hh.GetBinError(1, 2),
                    ]
                )

            canvas.Print(f"{OUT_DIR}/ABCD_{_bins[0]*100:03.0f}_{SAVE_STAT}.png")
            if VERBOSE:
                print(f"Saved 'ABCD_{_bins[0]*100:03.0f}_{SAVE_STAT}.png'")

        print("")
        alert("Performing ABCD Calculations", form="-", c="g")

        tables = []

        abcd_table_sig = Table(
            cols=["Pred A", "A", "B", "C", "D", "X2"],
            header="ABCD (X: dPhi, Y: dtSize)",
            rv_decimals=2,
            rv_join=" ± ",
            rl_label=["", "dPhi"],
        )
        abcd_table_bkg = Table(
            cols=["Pred A", "A", "B", "C", "D", "X2"],
            header="ABCD (X: dPhi, Y: dtSize)",
            rv_decimals=2,
            rv_join=" ± ",
            rl_label=["", "dPhi"],
        )

        table_rows_sig = [["", None]] * n_abcd * len(mss)
        table_rows_bkg = [["", None]] * n_abcd * len(mss)

        abcd_graphs = {}
        graph_axis_labels = ["Minimum |#Delta#phi_{CSC,DT}|", "events in SR (A)"]  # "#scale[0.5]{#int} L"]
        hmin_bkg, hmax_bkg = 1e9, -1
        hmin_sig, hmax_sig = 1e9, -1

        for ims, ms in enumerate(mss):
            k = ms.name
            abcd = abcds[k]

            ths, _as, aes, pas, paes = [], [], [], [], []
            for ith, (th, a, ae, b, be, _c, ce, d, de) in enumerate(abcd):
                pa, pae = 0, 0
                if b > 0 and _c > 0 and d > 0:  # TODO: Clean
                    pa = d * b / _c
                    pae = pa * ((be / b) ** 2 + (ce / _c) ** 2 + (de / d) ** 2) ** (1 / 2)
                    cl = (a - pa) ** 2 / (ae**2 + pae**2)  # chisq?
                    # cl = abs((a - pa) / (ae**2 + pae**2) ** (1 / 2))
                    if ("SIZE" in BLIND_TYPE or "BLINDSR" in BLIND_TYPE) and not ms.is_mc:
                        vals = [[pa, pae], "   ", [b, be], [_c, ce], [d, de], ""]
                    else:
                        vals = [[pa, pae], [a, ae], [b, be], [_c, ce], [d, de], cl]
                else:
                    if ("SIZE" in BLIND_TYPE or "BLINDSR" in BLIND_TYPE) and not ms.is_mc:
                        vals = ["", "   ", [b, be], [_c, ce], [d, de], ""]
                    else:
                        vals = ["", [a, ae], [b, be], [_c, ce], [d, de], ""]

                am, pam = (a - ae if a - ae > 0 else 1), (pa - pae if pa - pae > 0 else 1)
                ap, pap = a + ae, pa + pae

                hmin_sig = min(hmin_sig, am, pam)
                hmax_sig = max(hmax_sig, ap, pap)
                if ims > 0:
                    hmin_bkg = min(hmin_bkg, am, pam)
                    hmax_bkg = max(hmax_bkg, ap, pap)

                table_rows_sig[ith * len(mss) + ims] = [[k, th if k == "Signal" else ""], vals]
                if not ms.is_mc:
                    table_rows_bkg[ith * len(mss) + ims] = [[k, th if k == "Data" else ""], vals]

                ths.append(th - (0.04 * (len(mss) - 1 / 2) / 2) + 0.04 * ims)
                _as.append(a)
                aes.append(ae)
                pas.append(pa)
                paes.append(pae)

            abcd_graphs[k] = [
                create_TGraph(ths, _as, [0] * len(ths), aes, graph_axis_labels),
                create_TGraph([th + 0.02 for th in ths], pas, [0] * len(ths), paes, graph_axis_labels),
            ]
            gc.extend(abcd_graphs[k])

        # * Graphs * #

        for imss, _mss in enumerate([mss, mss[1:]]):
            cn = "c" + str(np.random.randint(999999999))
            c = rt.TCanvas(cn, cn, 800, 800)
            c.SetGrid()
            c.SetRightMargin(0.04)

            # leg = []
            # lat = rt.TLatex()
            # lat.SetTextAlign(11)
            # lat.SetTextSize(0.04)

            msn_max = max([len(ms.name) for ms in _mss])
            # leg = rt.TLegend((msn_max+8)*0.015,(2*len(_mss))*0.5)
            leg = rt.TLegend(0.17, 0.91 - 0.04 * (2 * len(_mss)), 0.21 + 0.016 * (msn_max + 11), 0.94)
            leg.SetTextFont(62)
            leg.SetTextSize(0.03)
            leg.SetBorderSize(0)
            leg.SetFillColorAlpha(1, 0.3)
            leg.SetEntrySeparation(0.01)

            if "Signal" in [ms.name for ms in _mss]:
                _min, _max = hmin_sig, hmax_sig
                label = "_sig"
            else:
                _min, _max = hmin_bkg, hmax_bkg
                label = "_bkg"

            if _max / _min > 1000:
                c.SetLogy()
                _min, _max = np.log10(_min if _min == 1 else 1e-1), np.log10(_max)
                _min, _max = _min - BOT_MARGIN * (_max - _min), _max + TOP_MARGIN * (_max - _min)
                _min, _max = 10 ** (_min), 10 ** (_max)
            else:
                _min, _max = _min - BOT_MARGIN * (_max - _min), _max + TOP_MARGIN * (_max - _min)
                _min = 0 if _min < 50 else _min

            for ims, ms in enumerate(_mss):
                for igr, gr in enumerate(abcd_graphs[ms.name]):
                    gr.SetMinimum(_min)
                    gr.SetMaximum(_max)

                    gr.SetLineColor(ms.colors[igr])
                    gr.SetLineWidth(4)
                    gr.SetLineStyle(9 if igr else 1)  # '--' if meas else '-'

                    gr.SetMarkerColor(ms.colors[igr])
                    gr.SetMarkerSize(2.5 if igr else 2)
                    gr.SetMarkerStyle(47 if igr else 20)  # x if meas else o

                    lab = f"{ms.name:>{msn_max}}: " + ("predicted" if igr else "observed")
                    le = leg.AddEntry(gr, lab, "PE")
                    # le.SetTextColor(ms.colors[igr])

                    gr.Draw(("AF" if igr + ims == 0 else "") + "PE same")

            leg.Draw()
            if not ROOT_BATCH:
                c.Draw()

            c.Print(f"{OUT_DIR}/closure{label}_{SAVE_STAT}.png")
            if VERBOSE:
                print(f"Saved 'closure{label}_{SAVE_STAT}.png'")

        # * Tables * #

        for table, rows in [(abcd_table_sig, table_rows_sig), (abcd_table_bkg, table_rows_bkg)]:
            last_spacer = 0
            for irow, (label, vals) in enumerate([row for row in rows if row[1] is not None]):
                if label[1] != "" and irow:  # and rows[irow-1][0][0] == "":
                    if irow - last_spacer > 1:
                        table.add_spacer()
                    last_spacer = irow
                table.add_row(label, vals)
            table.print()

        # ********** #

        wt = np.asarray(ms_r3["weight"])
        dtsize = np.asarray(ms_r3["dtSize"][:, 0])
        dphi = np.asarray(ms_r3["tag_dPhi"])

        for label in ("dphi", "dtsize", "min_dphi", "min_cscsize", "dtastn_lower", "dtastn_upper"):
            if label == "dphi":
                graph_axis_labels = ["|#Delta#phi_{CSC,DT}| ABCD Boundary", "events in SR"]  # "#scale[0.5]{#int} L"]
                ths = np.linspace(1.3, 1.7, 100)
                a_cond = lambda _th: (dtsize > ABCD_DTSIZE) & (dphi > _th)
                b_cond = lambda _th: (dtsize < ABCD_DTSIZE) & (dphi > _th)
                c_cond = lambda _th: (dtsize < ABCD_DTSIZE) & (dphi < _th)
                d_cond = lambda _th: (dtsize > ABCD_DTSIZE) & (dphi < _th)
            elif label == "dtsize":
                graph_axis_labels = ["DT_{SIZE} ABCD Boundary", "events in SR"]  # "#scale[0.5]{#int} L"]
                ths = np.linspace(100, 200, 100)
                a_cond = lambda _th: (dtsize > _th) & (dphi > ABCD_DPHI)
                b_cond = lambda _th: (dtsize < _th) & (dphi > ABCD_DPHI)
                c_cond = lambda _th: (dtsize < _th) & (dphi < ABCD_DPHI)
                d_cond = lambda _th: (dtsize > _th) & (dphi < ABCD_DPHI)
            elif label == "min_dphi":
                graph_axis_labels = ["|#Delta#phi_{CSC,DT}| Lower Cutoff", "events in SR"]  # "#scale[0.5]{#int} L"]
                if "DPHI" in CUTS:
                    ths = np.linspace(0.4, 1.2, 100)
                else:
                    ths = np.linspace(0.0, 1.2, 100)
                a_cond = lambda _th: (dtsize > ABCD_DTSIZE) & (dphi > ABCD_DPHI) & (_th < dphi)
                b_cond = lambda _th: (dtsize < ABCD_DTSIZE) & (dphi > ABCD_DPHI) & (_th < dphi)
                c_cond = lambda _th: (dtsize < ABCD_DTSIZE) & (dphi < ABCD_DPHI) & (_th < dphi)
                d_cond = lambda _th: (dtsize > ABCD_DTSIZE) & (dphi < ABCD_DPHI) & (_th < dphi)
            elif label == "min_cscsize":
                graph_axis_labels = ["CSC_{SIZE} Lower Cutoff", "events in SR"]  # "#scale[0.5]{#int} L"]
                var = np.asarray(ms_r3["cscSize"][:, 0])
                ths = np.linspace(100, 300, 100)
                a_cond = lambda _th: (dtsize > ABCD_DTSIZE) & (dphi > ABCD_DPHI) & (var > _th)
                b_cond = lambda _th: (dtsize < ABCD_DTSIZE) & (dphi > ABCD_DPHI) & (var > _th)
                c_cond = lambda _th: (dtsize < ABCD_DTSIZE) & (dphi < ABCD_DPHI) & (var > _th)
                d_cond = lambda _th: (dtsize > ABCD_DTSIZE) & (dphi < ABCD_DPHI) & (var > _th)
            elif label == "dtastn_lower":
                graph_axis_labels = ["DT Avg Station Lower Cutoff", "events in SR"]  # "#scale[0.5]{#int} L"]
                var = np.asarray(ms_r3["dtAvgStation10"][:, 0])
                ths = np.linspace(2, 4, 100)
                a_cond = lambda _th: (dtsize > ABCD_DTSIZE) & (dphi > ABCD_DPHI) & (var > _th)
                b_cond = lambda _th: (dtsize < ABCD_DTSIZE) & (dphi > ABCD_DPHI) & (var > _th)
                c_cond = lambda _th: (dtsize < ABCD_DTSIZE) & (dphi < ABCD_DPHI) & (var > _th)
                d_cond = lambda _th: (dtsize > ABCD_DTSIZE) & (dphi < ABCD_DPHI) & (var > _th)
            elif label == "dtastn_upper":
                graph_axis_labels = ["DT Avg Station Upper Cutoff", "events in SR"]  # "#scale[0.5]{#int} L"]
                var = np.asarray(ms_r3["dtAvgStation10"][:, 0])
                ths = np.linspace(2, 4, 100)
                a_cond = lambda _th: (dtsize > ABCD_DTSIZE) & (dphi > ABCD_DPHI) & (var < _th)
                b_cond = lambda _th: (dtsize < ABCD_DTSIZE) & (dphi > ABCD_DPHI) & (var < _th)
                c_cond = lambda _th: (dtsize < ABCD_DTSIZE) & (dphi < ABCD_DPHI) & (var < _th)
                d_cond = lambda _th: (dtsize > ABCD_DTSIZE) & (dphi < ABCD_DPHI) & (var < _th)

            _a, _ae, _pa, _pae = [], [], [], []
            for th in ths:
                a = np.sum(wt[a_cond(th)])
                b = np.sum(wt[b_cond(th)])
                c = np.sum(wt[c_cond(th)])
                d = np.sum(wt[d_cond(th)])
                pa, pae, ae, be, ce, de = 0, 0, np.sqrt(a), np.sqrt(b), np.sqrt(c), np.sqrt(d)
                # a/b = d/c, a = bd/c
                if b * c * d > 0:
                    pa = b * d / c
                    pae = pa * ((be / b) ** 2 + (ce / c) ** 2 + (de / d) ** 2) ** (1 / 2)
                _a.append(a)
                _ae.append(ae)
                _pa.append(pa)
                _pae.append(pae)
            _a, _ae, _pa, _pae = np.asarray(_a), np.asarray(_ae), np.asarray(_pa), np.asarray(_pae)

            cn = "c" + str(np.random.randint(999999999))
            canvas = rt.TCanvas(cn, cn, 800, 800)
            canvas.SetGrid()
            canvas.SetRightMargin(0.04)

            if BLIND_TYPE == "BLINDSR":
                leg = rt.TLegend(0.17, 0.89, 0.35, 0.94)
            else:
                leg = rt.TLegend(0.17, 0.84, 0.35, 0.94)

            leg.SetTextFont(62)
            leg.SetTextSize(0.03)
            leg.SetBorderSize(0)
            leg.SetFillColorAlpha(1, 0.2)
            leg.SetEntrySeparation(0.01)

            gr_m = create_TGraph(ths, _a, [ths[1] - ths[0]] * len(ths), _ae, graph_axis_labels)
            gr_p = create_TGraph(ths, _pa, [ths[1] - ths[0]] * len(ths), _pae, graph_axis_labels)

            for igr, gr in enumerate([gr_m, gr_p]):
                gr.SetMinimum(0)
                gr.SetMaximum(max(np.max(_a + _ae), np.max(_pa + _pae)) * (1 + TOP_MARGIN))

                gr.SetLineColor(rt.kRed if igr else rt.kBlack)
                gr.SetMarkerColor(rt.kRed if igr else rt.kBlack)
                gr.SetFillColorAlpha(rt.kRed if igr else rt.kBlack, 0.3)

                gr.SetLineWidth(4)
                gr.SetMarkerSize(2.5 if igr else 2)
                # gr.SetLineStyle(9 if igr else 1)  # '--' if meas else '-'
                gr.SetMarkerStyle(47 if igr else 20)  # x if meas else o

                if BLIND_TYPE == "BLINDSR":
                    if igr:
                        gr.Draw("A CE4")
                        leg.AddEntry(gr, "predicted", "PE")
                else:
                    gr.Draw(("same" if igr else "A") + " CE4")
                    leg.AddEntry(gr, "predicted" if igr else "observed", "PE")

            leg.Draw()
            if not ROOT_BATCH:
                c.Draw()

            canvas.Print(f"{OUT_DIR}/closure_scan_{label}_{SAVE_STAT}.png")
            if VERBOSE:
                print(f"Saved 'closure_scan_{label}_{SAVE_STAT}.png'")

    # ************************************************************ #
    # ************************************************************ #
    # ************************************************************ #

    if TRAIN_BDT:
        print("")
        alert("Training BDTs", form="-", c="g")

        # Convert to Numpy array
        feats_csc = [
            # 'cscSize',
            "cscPhi",
            "cscEta",
            # 'cscX',
            # 'cscY',
            "cscZ",
            "cscR",
            "cscNStation10",
            "cscAvgStation10",
            # 'cscMaxStation',
            # 'cscMe11Ratio',
            # 'cscMe12Ratio',
            "cscJetVetoPt",
            "cscMuonVetoPt",
        ]
        feats_dt = [
            # 'dtSize',
            "dtPhi",
            "dtEta",
            # 'dtX',
            # 'dtY',
            "dtZ",
            "dtR",
            "dtNStation10",
            "dtAvgStation10",
            # 'dtMaxStation',
            # 'dtMb1Ratio',
            "dtJetVetoPt",
            "dtMuonVetoPt",
        ]
        feats_2tag = [
            "tag_dR",
            "tag_dEta",
            "tag_dPhi",
        ]

        X_csc = np.array([np.r_[np.ravel(ms_mc[feat]), np.ravel(ms_r3[feat])] for feat in feats_csc]).T
        X_dt = np.array([np.r_[np.ravel(ms_mc[feat]), np.ravel(ms_r3[feat])] for feat in feats_dt]).T

        y = np.r_[
            np.ones(len(ms_mc["weight"]), dtype=bool),
            np.zeros(len(ms_r3["weight"]), dtype=bool),
        ]
        w = np.r_[np.ravel(ms_mc["weight"]), np.ravel(ms_r3["weight"])]

        X_csc, X_dt = np.abs(X_csc), np.abs(X_dt)

        X_trn_csc, X_tst_csc, X_trn_dt, X_tst_dt, y_trn, y_tst, w_trn, w_tst = train_test_split(
            X_csc, X_dt, y, w, test_size=TEST_SIZE, random_state=42
        )
        nmc_trn, nr3_trn = np.sum(y_trn), len(y_trn) - np.sum(y_trn)
        # n_estimators = np.sqrt(min(len(feats_csc), len(feats_dt))*min(np.sum(y_csc), len(y_csc)-np.sum(y_csc)))
        n_estimators_csc = np.sqrt(len(feats_csc) * len(y_trn))
        n_estimators_dt = np.sqrt(len(feats_dt) * len(y_trn))
        # n_estimators = np.sqrt(min(len(feats_csc), len(feats_dt))*min(np.sum(y_csc), len(y_csc)-np.sum(y_csc)))
        n_estimators_csc, n_estimators_dt = int(n_estimators_csc), int(n_estimators_dt)
        if VERBOSE:
            print(f"n feats: csc = {len(feats_csc):,}, dt = {len(feats_dt):,}")
            print(f"n trn: mc = {nmc_trn:,}, r3 = {nr3_trn:,} (mc/r3 = {nmc_trn/(nmc_trn+nr3_trn):.3f})")
            print(f"n_estimators: csc = {n_estimators_csc:,}, dt = {n_estimators_dt:,}")

        # clf_csc = RandomForestClassifier(n_estimators=n_estimators, random_state=42, n_jobs=2)
        # clf_dt = RandomForestClassifier(n_estimators=n_estimators, random_state=42, n_jobs=2)
        clf_csc = GradientBoostingClassifier(n_estimators=n_estimators_csc, random_state=42)  # , max_depth=10)
        clf_dt = GradientBoostingClassifier(n_estimators=n_estimators_dt, random_state=42)  # , max_depth=10)

        sclr_csc, sclr_dt = StandardScaler(), StandardScaler()
        X_trn_csc = sclr_csc.fit_transform(X_trn_csc)
        X_tst_csc = sclr_csc.transform(X_tst_csc)
        X_csc = sclr_csc.transform(X_csc)

        X_trn_dt = sclr_dt.fit_transform(X_trn_dt)
        X_tst_dt = sclr_dt.transform(X_tst_dt)
        X_dt = sclr_dt.transform(X_dt)

        clf_csc.fit(X_trn_csc, y_trn, w_trn)
        clf_dt.fit(X_trn_dt, y_trn, w_trn)

        try:
            pred_csc = clf_csc.decision_function(X_csc)
            pred_dt = clf_dt.decision_function(X_dt)
        except AttributeError as e:
            pred_csc = clf_csc.predict_proba(X_csc)[:, 0]
            pred_dt = clf_dt.predict_proba(X_dt)[:, 0]

        ms_mc["GBT_CSC"] = pred_csc[y]
        ms_mc["GBT_DT"] = pred_dt[y]

        ms_r3["GBT_CSC"] = pred_csc[~y]
        ms_r3["GBT_DT"] = pred_dt[~y]

        for clfn, fts, clf in [("CSC", feats_csc, clf_csc), ("DT", feats_dt, clf_dt)]:
            feat_table = Table(["Importance"], f"GBT {clfn}", 4, rl_label=["", "Feature"])
            fts, wts = np.asarray(fts), clf.feature_importances_
            idxs = np.argsort(wts)[::-1]
            for ift, (ft, wt) in enumerate(zip(fts[idxs], wts[idxs])):
                feat_table.add_row([ift + 1, ft], [wt])
            feat_table.print()

    # **************************** #

    if "ROCS" in SECTIONS:
        print("")
        alert("Making ROC Plots", form="-", c="g")

        bkgs, sigs, aucs = {}, {}, {}
        bkg_effs, sig_effs = {}, {}
        names = [
            "met",
            "tag_dPhi",
            "tag_dEta",
            "tag_dR",
            "dtPhi.abs",
            "cscSize",
            "dtSize",
            "cscNStation10",
            "dtNStation10",
        ]

        if TRAIN_BDT:
            names.extend(["GBT_CSC", "GBT_DT"])

        clf_table = Table(
            ["AUC", "S/√[B]", "Thresh", "Sig", "Sig Eff", "Bkg", "Bkg Eff"], "Evaluating Discriminators", 3
        )

        wmc, wr3 = ms_mc["weight"], ms_r3["weight"]
        for i, name in enumerate(names):
            bkgs[name], sigs[name], aucs[name] = [], [], 0
            bkg_effs[name], sig_effs[name] = [], []
            vmc, vr3 = np.asarray(ms_mc[name]), np.asarray(ms_r3[name])
            if len(vmc.shape) == 2:
                vmc, vr3 = vmc[:, 0], vr3[:, 0]

            _min, _max = np.min(np.r_[vmc, vr3]), np.max(np.r_[vmc, vr3])
            threshold = np.arange(_min, _max, (_max - _min) / 100)

            for th in threshold:
                cond = vr3 > th
                bkg_effs[name].append(np.sum(wr3[cond]) / np.sum(wr3))
                bkgs[name].append(np.sum(wr3[cond]))

                cond = vmc > th
                sig_effs[name].append(np.sum(wmc[cond]) / np.sum(wmc))
                sigs[name].append(np.sum(wmc[cond]))

            sig_effs[name] = np.array(sig_effs[name])
            bkg_effs[name] = np.array(bkg_effs[name])
            sigs[name] = np.array(sigs[name])
            bkgs[name] = np.array(bkgs[name])
            aucs[name] = roc_auc_score(
                y_true=np.r_[np.ones_like(vmc), np.zeros_like(vr3)],
                y_score=np.r_[vmc, vr3],
                sample_weight=np.r_[wmc, wr3],
            )

            sigs[name] = sigs[name][bkg_effs[name] > 0]
            bkgs[name] = bkgs[name][bkg_effs[name] > 0]
            sig_effs[name] = sig_effs[name][bkg_effs[name] > 0]
            bkg_effs[name] = bkg_effs[name][bkg_effs[name] > 0]

            # sig, bkg, sig_eff, bkg_eff, auc = bkgs[name], sig_effs[name], bkg_effs[name], aucs[name]

            auc, idx = aucs[name], np.argmax(sigs[name] / (np.sqrt(bkgs[name])))
            th = threshold[idx]
            sig, sig_eff = sigs[name][idx], sig_effs[name][idx]
            bkg, bkg_eff = bkgs[name][idx], bkg_effs[name][idx]

            clf_table.add_row(name, [auc, sig / np.sqrt(bkg), th, sig, sig_eff, bkg, bkg_eff])

        clf_table.sort(lambda x: x[1].replace(",", ""))
        clf_table.print()
        # **************************** #

        cn = "c" + str(np.random.randint(999999999))
        c = rt.TCanvas(cn, cn, 800, 800)
        leg = []
        lat = rt.TLatex()
        lat.SetTextAlign(31)
        lat.SetTextSize(0.04)

        graph, ymax = {}, -1
        for i, v in enumerate(names):
            graph[v] = create_TGraph(sig_effs[v], 1 / bkg_effs[v], axis_title=["signal efficiency", "bkg rejection"])
            graph[v].SetLineWidth(5)
            graph[v].SetLineColor(std_color_list[i])
            leg.append([v, std_color_list[i]])
            # graph[v].SetMaximum(50)
            graph[v].Draw("ac" if i == 0 else "c same")
            ymax = max(ymax, np.max(1 / bkg_effs[v]))

        for ileg, (text, color) in enumerate(leg):
            lat.SetTextColor(color)
            # lat.DrawLatexNDC(0.94, 0.92 - ileg*0.04, text)
            # lat.DrawLatexNDC(0.94, 0.92 - ileg*(0.92-0.08)/(len(leg)+1), text)
            lat.DrawLatexNDC(0.95, 0.92 - ileg * (0.92 - 0.40) / (len(leg) + 1), text)

        if ymax > 100:
            c.SetLogy()
        c.SetGrid()
        # c.SetRightMargin(0.04)

        if not ROOT_BATCH:
            c.Draw()

        c.Print(f"{OUT_DIR}/rocs_{SAVE_STAT}.png")
        if VERBOSE:
            print(f"Saved 'rocs_{SAVE_STAT}.png'")

    # ************************************************************ #
    # ************************************************************ #
    # ************************************************************ #
    # print("")
    # alert("Getting run number of high hit clusters", form="-", c="g")

    # ************************************************************ #


if __name__ == "__main__":
    alert("Starting train_bdt.py", form="=", c="c")

    # **************************** #
    if len(sys.argv) > 1:
        N_EVENTS = int(sys.argv[1])
    if len(sys.argv) > 2:
        OUT_DIR = OUT_DIR + f"/{sys.argv[2]}"
    if len(sys.argv) > 3:
        BLIND_TYPE = sys.argv[3]
        SAVE_STAT = BLIND_TYPE
    # SAVE_STAT = f"{SAVE_STAT}_N{N_EVENTS // 1_000_000:2.0f}"

    alert(f"  VERBOSE      = {VERBOSE}", c="c")
    alert(f"  SAVE_STAT    = '{SAVE_STAT}'", c="c")
    alert(f"  OUT_DIR      = '{OUT_DIR}'", c="c")
    alert(f"  BLIND_TYPE   = '{BLIND_TYPE}'", c="c")
    alert(f"  TRAIN_BDT    = {TRAIN_BDT}", c="c")
    alert(f"  TEST_SIZE    = {TEST_SIZE}", c="c")
    # alert(f"  TEST_SIZE    = {TEST_SIZE:.2f}", c="c")
    alert(f"  N_EVENTS     = {N_EVENTS:,}", c="c")
    alert(f"  DATA_VERSION = {DATA_VERSION}", c="c")
    alert(f"  LUMI         = {LUMI:,}", c="c")
    # alert(f"  LUMI         = {LUMI:,.2f}", c="c")
    alert(f"  DATA_DIR     = '{DATA_DIR}'", c="c")
    # **************************** #
    if "TIER2" in DATA_DIR:
        OUT_DIR = f"{T2_OUT_DIR}/{OUT_DIR}"
        FN_MC = f"{T2_DATA_DIR}/MC_Summer22EE/v1/sixie/v{DATA_VERSION}/normalized/{FN_R3}.root"
        FN_R3 = f"{T2_DATA_DIR}/Data2022/v{DATA_VERSION}/normalized/{FN_R3}.root"

        FN_HLT_MC = f"{LOCAL_DATA_DIR}/../processed/mc_hlt569.rdf"  #! broken
        FN_HLT_R3 = f"{LOCAL_DATA_DIR}/../processed/r3_hlt569.rdf"  #! broken
    else:
        OUT_DIR = f"{LOCAL_OUT_DIR}/{OUT_DIR}"
        FN_MC = f"{LOCAL_DATA_DIR}/{FN_MC}_v{DATA_VERSION}.root"
        FN_R3 = f"{LOCAL_DATA_DIR}/{FN_R3}_v{DATA_VERSION}.root"

        FN_HLT_MC = f"{LOCAL_DATA_DIR}/../processed/mc_hlt569.rdf"
        FN_HLT_R3 = f"{LOCAL_DATA_DIR}/../processed/r3_hlt569.rdf"

    if "HLT" in CUTS:
        FN_MC, FN_R3 = FN_MC, FN_HLT_R3
    else:
        FN_MC, FN_R3 = FN_MC, FN_R3
    pathlib.Path(OUT_DIR).mkdir(parents=True, exist_ok=True)
    # **************************** #
    rt.gErrorIgnoreLevel = ROOT_ERROR_LEVEL
    if ROOT_BATCH:
        rt.gROOT.SetBatch(ROOT_BATCH)
    tdrstyle.setTDRStyle()
    CMS_lumi.writeExtraText = 0

    # **************************** #
    # **************************** #
    # **************************** #

    main()

    print("")
    alert("Finished train_bdt.py", form="=", c="c")
