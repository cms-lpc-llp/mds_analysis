"""train_bdt.py

To run :
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
################################
##         Parameters         ##
################################
SAVE_STAT = "dtOOT"
BLIND_TYPE = "DTOOT"  # , "DTSIZE", "DYNAMIC", "NONE"
N_EVENTS = -1

TRAIN_BDT = True
TEST_SIZE = 0.8

CUTS = [
    "match",
    # "HLT", # TODO: Load HLTDecision without overflowing memory
    "CSC&DT>0",
    "L1",
    "CSCIT",
    "DTIT",
    # "MET",
    "ME1",
    "MB1",
    # "JET",
    # "MUON",
    # "BDT",
    "HALO",
    "1CSC1DT",
    "DPHI",
]
################################
ABCD_DPHI = 1.5
ABCD_DTSIZE = 200
################################
OUT_DIR = "reports/weekly/2023-09-07"
T2_OUT_DIR = "/storage/af/user/psimmerl/LLP/mdc_analysis"  # os.getcwd()
LOCAL_OUT_DIR = "/home/psimmerl/LLP/mdc_analysis"  # os.getcwd()
################################
DATA_VERSION = "6"
LUMI = 23.02 * 1000

T2_DATA_DIR = "/storage/cms/store/user/christiw/displacedJetMuonAnalyzer/Run3/V1p19"
LOCAL_DATA_DIR = "/home/psimmerl/LLP/mdc_analysis/data/raw"
DATA_DIR = "TIER2" if "caltech" in uname()[1] else "LOCAL"

FN_MC = "ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted"
FN_R3 = "DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi"
################################
VERBOSE = True
ROOT_BATCH = True
ROOT_ERROR_LEVEL = 1001 #rt.kInfo + 1
gc = []

################################

# def hist(samples, sample_labels, plot_labels, weights=None, bins=None, log=False, colors="std_color_list", abs=False, canvas=None, save=True):
#     cn = "c" + str(np.random.randint(999999999))
#     c = rt.TCanvas(cn, cn, 800, 800)

#     leg = []
#     lat = rt.TLatex()
#     lat.SetTextAlign(11)
#     lat.SetTextSize(0.04)

#     hmin, hmax = 1e9 - 1, -1

#     for ims, ms in enumerate(mss):
#         k, var, weight = ms.name, ms[xl], ms["weight"]
#         if abs:
#             var = np.abs(var)

#         h = create_TH1D(var, axis_title=[xl, "fraction of events"], binning=bins[ix], weights=weight)
#         h.SetLineColor(std_color_list[ims])

#         # print(k, h.Integral(), np.sum(weight[np.abs(var)>1.5])/np.sum(weight))
#         leg.append((k, std_color_list[ims]))

#         if h.Integral() > 0:
#             h.Scale(1.0 / h.Integral())

#         hmax = max(hmax, h.GetMaximum())
#         hmin = min(hmin, h.GetMinimum(0))
#         if xl in ("cscPhi",):
#             hmin = 0

#         # h.SetMinimum(1e-3)
#         # h.SetMaximum(4e-1)
#         h.Draw("hist same")
#         gc.append(h)

#     for ileg, (text, color) in enumerate(leg):
#         gc[-(ileg + 1)].SetMaximum(hmax)
#         gc[-(ileg + 1)].SetMinimum(hmin)
#         lat.SetTextColor(color)
#         lat.DrawLatexNDC(0.80 - 0.1 * ileg, 0.92, text)

#     if hmax/hmin >= 100:
#         c.SetLogy()
#     c.SetGrid()
#     c.SetRightMargin(0.04)

#     if not ROOT_BATCH:
#         c.Draw()

#     c.Print(f"{OUT_DIR}/{xl}_{SAVE_STAT}.png")
#     if VERBOSE:
#         print(f"Saved '{xl}_{SAVE_STAT}.png'")

# def save(canvas, plots, file, label):
#     canvas.save

################################

if __name__ == "__main__":
    alert("Starting train_bdt.py", form="=", c="c")

    if len(sys.argv) > 1:
        N_EVENTS = int(sys.argv[1])
    if len(sys.argv) > 2:
        OUT_DIR = OUT_DIR + f"/{sys.argv[2]}"

    # SAVE_STAT = f"{SAVE_STAT}_N{N_EVENTS // 1_000_000}"

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

    ################################
    if "TIER2" in DATA_DIR:
        OUT_DIR = f"{T2_OUT_DIR}/{OUT_DIR}"
        ff_mc = f"{T2_DATA_DIR}/MC_Summer22EE/v1/sixie/v{DATA_VERSION}/normalized/{FN_R3}.root"
        ff_r3 = f"{T2_DATA_DIR}/Data2022/v{DATA_VERSION}/normalized/{FN_R3}.root"
    else:
        OUT_DIR = f"{LOCAL_OUT_DIR}/{OUT_DIR}"
        ff_mc = f"{LOCAL_DATA_DIR}/{FN_MC}_v{DATA_VERSION}.root"
        ff_r3 = f"{LOCAL_DATA_DIR}/{FN_R3}_v{DATA_VERSION}.root"
    pathlib.Path(OUT_DIR).mkdir(parents=True, exist_ok=True)
    ################################
    rt.gErrorIgnoreLevel = ROOT_ERROR_LEVEL
    if ROOT_BATCH:
        rt.gROOT.SetBatch(ROOT_BATCH)
    tdrstyle.setTDRStyle()
    CMS_lumi.writeExtraText = 0
    ################################################################
    ################################################################
    ################################################################

    print("")
    alert("Building MuonSystem(s)", form="-", c="g")
    ms_mc = MuonSystemAwkward(ff_mc, name="signal", nev=N_EVENTS, is_mc=True, lumi=LUMI)
    ms_r3 = MuonSystemAwkward(ff_r3, name="data", nev=N_EVENTS, is_mc=False, lumi=LUMI)
    mss, ms_snames = [ms_mc, ms_r3], ["Sig", "Bkg"]

    if "standard" in sys.argv[2]:
        CUTS.remove("DPHI")
    elif "neither" in sys.argv[2]:
        CUTS.remove("HALO")
        CUTS.remove("DPHI")
    elif "halo" in sys.argv[2]:
        CUTS.remove("DPHI")
        ms_cut = MuonSystemAwkward(ff_r3, name="halo", nev=N_EVENTS, is_mc=False, lumi=LUMI)
        mss.append(ms_cut)
        ms_snames.append("Halo")
    elif "dphi" in sys.argv[2]:
        CUTS.remove("HALO")
        ms_cut = MuonSystemAwkward(ff_r3, name="dphi", nev=N_EVENTS, is_mc=False, lumi=LUMI)
        mss.append(ms_cut)
        ms_snames.append("dPhi")
    elif "both" in sys.argv[2]:
        pass
        # ms_cut = MuonSystemAwkward(ff_r3, name="halo&dphi", nev=N_EVENTS, is_mc=False, lumi=LUMI)
        # mss.append(ms_cut)
        # ms_snames.append("HO & dF")
    else:
        pass

    ms_names = [ms.name for ms in mss]
    ################################################################

    for ms in mss:
        #!!! TURNING MuonSystem CUTS OFF !!!!#
        # cuts are not applied when you get MuonSystem data
        # Note, they are 'AND'ed with sel_evt, sel_csc, and sel_dt.
        ms.cut = False
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

    ################################################################

    print("")
    alert("Cutting MuonSystem(s)", form="-", c="g")

    columns = [f"{ms.name} (eff.)" for ms in mss]
    cut_table = Table(columns, "MuonSystem Cuts (req. CSC>0 & DT>0 after each cut)", 1, rv_join=" ")
    eff_decimals = 4

    raw_counts = [ms.count() for ms in mss]
    cut_table.add_row("read", raw_counts)

    for cut in CUTS:
        counts = [ ["", ""] for i in range(len(mss))]
        for ims, ms in enumerate(mss):
            was_cut, mc, msn = True, ms.is_mc, ms.name
            match cut:
                case "match":
                    if mc:
                        ms.match_mc("csc,dt", has_clusters=True)
                    else:
                        was_cut = False
                case "CSC&DT>0":
                    ms.f((ms["nCsc"]>0) & (ms["nDt"]>0), has_clusters=False)
                    raw_counts[ims] = ms.count()

                    ms.efficiency_denom = "eff"
                    ms.ms_read["eff_evt"] = ms.ms_read["sel_evt"]
                    ms.ms_read["eff_csc"] = ms.ms_read["sel_csc"]
                    ms.ms_read["eff_dt"] = ms.ms_read["sel_dt"]
                    ms.efficiency = 1

                case "HLT":
                    ms.cut_hlt()
                case "ME1":
                    ms.f(ms["cscMe11Ratio"] + ms["cscMe12Ratio"] == 0, "csc")
                case "MB1":
                    ms.f(ms["dtNHitStation1"] == 0, "dt")
                case "L1":
                    ms.cut_l1()
                case "CSCIT":
                    ms.cut_time("csc", invert=False)
                case "DTIT":
                    ms.cut_time("dt", invert=(not mc and BLIND_TYPE=="DTOOT"))
                case "HALO":
                    ms.cut_halo(invert="halo" in msn)
                case "1CSC1DT":
                    ms.tag(tags="cscdt")
                case "DPHI":
                    ms.f(0.4 < ms["tag_dPhi"], invert="dphi" in msn)
                case default:
                    raise ValueError(f"Cut '{cut}' not recognized!", c="r")
            if was_cut:
                counts[ims][0] = ms.count()
                if cut in ("match", "HLT"): # skip efficiency
                    counts[ims] = counts[ims][0]
                else:
                    counts[ims][1] = f"({ms.efficiency:.{eff_decimals}f})"
        cut_label = ("DTIT/OOT" if cut =="DTIT" and BLIND_TYPE=="DTOOT" else cut)
        cut_table.add_row(cut_label, counts)
        if cut in ("match", ):
            cut_table.add_spacer()
        if VERBOSE:
            print(f"Finished cut '{cut}'")
    cut_table.print()

    ################################################################

    for ms in mss:
        #!!! TURNING MuonSystem CUTS ON !!!!#
        # cuts are now applied when you MuonSystem get data
        ms.cut = True
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

        #!!! TURNING EFFICIENCY OFF !!!!#
        del ms.ms_read[f"{ms.efficiency_denom}_evt"]
        del ms.ms_read[f"{ms.efficiency_denom}_csc"]
        del ms.ms_read[f"{ms.efficiency_denom}_dt"]
        ms.efficiency_denom = None
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

        # Adding new variables (some require 2tag)
        ms["cscR-dtR"] = ms["cscR"] - ms["dtR"]

    ################################################################
    ################################################################
    ################################################################

    print("")
    alert("Making CSC Plots", form="-", c="g")

    xls = [
        "cscMuonVetoPt",
        "cscJetVetoPt",
        "cscR",
        "cscZ",
        "cscPhi",
        "cscEta",
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

    for ix, xl in enumerate(xls):
        cn = "c" + str(np.random.randint(999999999))
        c = rt.TCanvas(cn, cn, 800, 800)
        leg = []
        lat = rt.TLatex()
        lat.SetTextAlign(11)
        lat.SetTextSize(0.04)
        hmin, hmax = 1e9 - 1, -1

        for ims, ms in enumerate(mss):
            k, var, weight = ms.name, ms[xl], ms["weight"]
            if xl in ("cscEta", "cscPhi", "cscZ"):
                var = np.abs(var)

            # if "tag_dPhi" in xl or "halo" in k:
            #     cond = ms["tag_dPhi"] > 0.0
            # else:
            #     cond = ms["tag_dPhi"] > 0.5
            cond = ms["tag_dPhi"] > -1

            if "size" in BLIND_TYPE:
                cond = cond & (ms["dtSize"][:, 0] < 100)

            if not ms.is_mc:
                var, weight = var[cond], weight[cond]

            h = create_TH1D(var, axis_title=[xl, "fraction of events"], binning=bins[ix], weights=weight)
            h.SetLineColor(std_color_list[ims])

            # print(k, h.Integral(), np.sum(weight[np.abs(var)>1.5])/np.sum(weight))
            leg.append((k, std_color_list[ims]))

            if h.Integral() > 0:
                h.Scale(1.0 / h.Integral())

            hmax = max(hmax, h.GetMaximum())
            hmin = min(hmin, h.GetMinimum(0))
            if xl in ("cscPhi",):
                hmin = 0

            # h.SetMinimum(1e-3)
            # h.SetMaximum(4e-1)
            h.Draw("hist same")
            gc.append(h)

        for ileg, (text, color) in enumerate(leg):
            gc[-(ileg + 1)].SetMaximum(hmax)
            gc[-(ileg + 1)].SetMinimum(hmin)
            lat.SetTextColor(color)
            lat.DrawLatexNDC(0.80 - 0.1 * ileg, 0.92, text)

        if xl not in ("cscPhi",):
            c.SetLogy()
        c.SetGrid()
        c.SetRightMargin(0.04)

        if not ROOT_BATCH:
            c.Draw()

        c.Print(f"{OUT_DIR}/{xl}_{SAVE_STAT}.png")
        if VERBOSE:
            print(f"Saved '{xl}_{SAVE_STAT}.png'")

    ################################

    print("")
    alert("Making DT Plots", form="-", c="g")

    xls = [
        "dtJetVetoPt",
        "dtMuonVetoPt",
        "dtR",
        "dtZ",
        "dtPhi",
        "dtEta",
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

    for ix, xl in enumerate(xls):
        cn = "c" + str(np.random.randint(999999999))
        c = rt.TCanvas(cn, cn, 800, 800)
        leg = []
        lat = rt.TLatex()
        lat.SetTextAlign(11)
        lat.SetTextSize(0.04)
        hmin, hmax = 1e9 - 1, -1

        for ims, ms in enumerate(mss):
            k, var, weight = ms.name, ms[xl], ms["weight"]
            if xl in ("dtEta", "dtPhi", "dtZ"):
                var = np.abs(var)

            # if "tag_dPhi" in xl or "halo" in k:
            #     cond = ms["tag_dPhi"] > 0.0
            # else:
            #     cond = ms["tag_dPhi"] > 0.5
            cond = ms["tag_dPhi"] > -1

            if "size" in BLIND_TYPE:
                cond = cond & (ms["cscSize"][:, 0] < 300)

            if not ms.is_mc:
                var, weight = var[cond], weight[cond]

            h = create_TH1D(var, axis_title=[xl, "fraction of events"], binning=bins[ix], weights=weight)
            h.SetLineColor(std_color_list[ims])

            # print(k, h.Integral(), np.sum(weight[np.abs(var)>1.5])/np.sum(weight))
            leg.append((k, std_color_list[ims]))

            if h.Integral() > 0:
                h.Scale(1.0 / h.Integral())

            hmax = max(hmax, h.GetMaximum())
            hmin = min(hmin, h.GetMinimum(0))
            if xl in ("dtPhi",):
                hmin = 0

            # h.SetMinimum(1e-3)
            # h.SetMaximum(4e-1)
            h.Draw("hist same")
            gc.append(h)

        for ileg, (text, color) in enumerate(leg):
            gc[-(ileg + 1)].SetMaximum(hmax)
            gc[-(ileg + 1)].SetMinimum(hmin)
            lat.SetTextColor(color)
            lat.DrawLatexNDC(0.80 - 0.1 * ileg, 0.92, text)

        if xl not in ("dtPhi",):
            c.SetLogy()
        c.SetGrid()
        c.SetRightMargin(0.04)

        if not ROOT_BATCH:
            c.Draw()

        c.Print(f"{OUT_DIR}/{xl}_{SAVE_STAT}.png")
        if VERBOSE:
            print(f"Saved '{xl}_{SAVE_STAT}.png'")

    ################################

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

    for ix, xl in enumerate(xls):
        cn = "c" + str(np.random.randint(999999999))
        c = rt.TCanvas(cn, cn, 800, 800)
        leg = []
        lat = rt.TLatex()
        lat.SetTextAlign(11)
        lat.SetTextSize(0.04)
        hmin, hmax = 1e9 - 1, -1

        for ims, ms in enumerate(mss):
            k, var, weight = ms.name, ms[xl], ms["weight"]
            # if xl in ("cscR-dtR", ""):
            #     var = np.abs(var)

            # if "tag_dPhi" in xl or "halo" in k:
            #     cond = ms["tag_dPhi"] > 0.0
            # else:
            #     cond = ms["tag_dPhi"] > 0.5
            cond = ms["tag_dPhi"] > -1

            if "size" in BLIND_TYPE:
                if "dt" in xl.lower():
                    cond = cond & (ms["cscSize"][:, 0] < 300)
                else:
                    cond = cond & (ms["dtSize"][:, 0] < 100)
            if not ms.is_mc:
                var, weight = var[cond], weight[cond]

            h = create_TH1D(var, axis_title=[xl, "fraction of events"], binning=bins[ix], weights=weight)
            h.SetLineColor(std_color_list[ims])

            # print(k, h.Integral(), np.sum(weight[np.abs(var)>1.5])/np.sum(weight))
            leg.append((k, std_color_list[ims]))

            if h.Integral() > 0:
                h.Scale(1.0 / h.Integral())

            hmax = max(hmax, h.GetMaximum())
            hmin = min(hmin, h.GetMinimum(0))
            # if xl in ("cscPhi",):
            #     hmin = 0

            # h.SetMinimum(1e-3)
            # h.SetMaximum(4e-1)
            h.Draw("hist same")
            gc.append(h)

        for ileg, (text, color) in enumerate(leg):
            gc[-(ileg + 1)].SetMaximum(hmax)
            gc[-(ileg + 1)].SetMinimum(hmin)
            lat.SetTextColor(color)
            lat.DrawLatexNDC(0.80 - 0.1 * ileg, 0.92, text)

        if xl not in ("dtPhi",):
            c.SetLogy()
        c.SetGrid()
        c.SetRightMargin(0.04)

        if not ROOT_BATCH:
            c.Draw()

        c.Print(f"{OUT_DIR}/{xl}_{SAVE_STAT}.png")
        if VERBOSE:
            print(f"Saved '{xl}_{SAVE_STAT}.png'")

    ################################

    print("")
    alert("Making 2D Plots", form="-", c="g")

    kabcds = []

    xyls = [
        ("cscPhi", "cscR"),
        ("cscPhi", "cscEta"),
        ("cscEta", "cscR"),
        #
        ("dtPhi", "dtR"),
        ("dtPhi", "dtEta"),
        ("dtEta", "dtR"),
        #
        ("cscR", "dtR"),
        ("cscPhi", "dtPhi"),
        ("cscEta", "dtEta"),
        #
        ("tag_dPhi", "met"),
        ("tag_dPhi", "tag_dEta"),
        ("tag_dPhi", "tag_dR"),
        ("tag_dEta", "tag_dR"),
        #
        ("tag_dPhi", "cscSize"),
        ("tag_dPhi", "cscPhi"),
        ("tag_dPhi", "cscEta"),
        #
        ("tag_dPhi", "dtSize"),
        ("tag_dPhi", "dtPhi"),
        ("tag_dPhi", "dtEta"),
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

    n_abcd = 0
    for th in np.arange(0, 1.0+0.1, step=0.1):
        if "DPHI" in CUTS and th < 0.4:
            continue
        xyls.append(("tag_dPhi", "dtSize"))
        bins.append([th, ABCD_DPHI, np.pi, 50, ABCD_DTSIZE, 500, 3, 3])
        n_abcd += 1

    abcds = {ms.name : [] for ms in mss}
    for ii, (xl, yl) in enumerate(xyls):
        cn = "c" + str(np.random.randint(999999999))
        c = rt.TCanvas(cn, cn, 800 * len(mss), 800)
        c.Divide(len(mss), 1)
        leg = []
        lat = rt.TLatex()
        lat.SetTextAlign(21)  # left,cent,right=1,2,3 , bot,cent,top=1,2,3
        lat.SetTextSize(0.06)
        hmin, hmax = 1e10 - 1, -1

        if len(bins[ii]) == 6:
            plot_name = f"{xl}_{yl}"
        else:  # custom bins -> abcd
            plot_name = f"ABCD_{bins[ii][0]*100:03.0f}"

        for ims, ms in enumerate(mss):
            c.cd(ims + 1)
            k, xvar, yvar, weight = ms.name, ms[xl], ms[yl], ms["weight"]

            if xl in ("cscEta", "cscPhi", "cscZ", "dtPhi", "dtEta", "dtZ"):
                xvar = np.abs(xvar)
            if yl in ("cscEta", "cscPhi", "cscZ", "dtPhi", "dtEta", "dtZ"):
                yvar = np.abs(yvar)

            # if "tag_dPhi" in xl or "tag_dPhi" in yl or "halo" in k:
            #     cond = ms["tag_dPhi"] > 0.0
            # else:
            #     cond = ms["tag_dPhi"] > 0.5
            cond = ms["tag_dPhi"] > -1

            if "size" in BLIND_TYPE:
                if ("dt" in xl or "dt" in yl) and not ("csc" in xl or "csc" in yl):
                    cond = cond & (ms["cscSize"][:, 0] < 300)
                else:
                    cond = cond & (ms["dtSize"][:, 0] < 100)

            if not ms.is_mc:
                xvar, yvar, weight = xvar[cond], yvar[cond], weight[cond]

            var = np.c_[np.asarray(xvar).flat, np.asarray(yvar).flat]

            h = create_TH2D(var, axis_title=[xl, yl, "fraction of events"], binning=bins[ii], weights=weight)
            # h.SetLineColor(std_color_list[ims])

            # print(k, h.Integral(), np.sum(weight[np.abs(var)>1.5])/np.sum(weight))
            leg.append((k, std_color_list[ims]))

            if "ABCD" in plot_name:
                a, ae = h.GetBinContent(2, 2), h.GetBinError(2, 2)
                b, be = h.GetBinContent(2, 1), h.GetBinError(2, 1)
                _c, ce = h.GetBinContent(1, 1), h.GetBinError(1, 1)
                d, de = h.GetBinContent(1, 2), h.GetBinError(1, 2)
                abcds[k].append([bins[ii][0], a, ae, b, be, _c, ce, d, de])

            if h.Integral() > 0:
                h.Scale(1.0 / h.Integral())

            hmax = max(hmax, h.GetMaximum())
            hmin = min(hmin, h.GetMinimum(0))

            h.Draw("colz" + (" text" if "ABCD" in plot_name else ""))
            gc.append(h)

        for ileg, (text, color) in enumerate(leg):
            gc[-(ileg + 1)].SetMaximum(hmax)
            gc[-(ileg + 1)].SetMinimum(hmin)
            c.cd(ileg + 1).SetLogz()
            c.cd(ileg + 1).SetGrid()
            # c.cd(ileg + 1).SetRightMargin(0.04)
            lat.SetTextColor(color)
            lat.DrawLatexNDC(0.5, 0.96, text)

        if not ROOT_BATCH:
            c.Draw()

        c.Print(f"{OUT_DIR}/{plot_name}_{SAVE_STAT}.png")
        if VERBOSE:
            print(f"Saved '{plot_name}_{SAVE_STAT}.png'")

    ################################
    ################################
    ################################

    print("")
    alert("Performing ABCD Calculations", form="-", c="g")

    cn = "c" + str(np.random.randint(999999999))
    c = rt.TCanvas(cn, cn, 800, 800)
    # leg = []
    # lat = rt.TLatex()
    # lat.SetTextAlign(11)
    # lat.SetTextSize(0.04)
    leg = rt.TLegend(0.12+0.4,0.7,0.37+0.5,0.87)
    leg.SetTextSize(0.03)
    leg.SetBorderSize(0)
    leg.SetEntrySeparation(0.01)
    # hmin, hmax = 1e9 - 1, -1

    table_format = [ None ] * n_abcd * len(mss)

    for ims, (ms, sn) in enumerate(zip(mss, ms_snames)):
        k = ms.name
        abcd = abcds[k]

        ths, _as, aes, pas, paes = [], [], [], [], []
        for ith, (th, a, ae, b, be, _c, ce, d, de) in enumerate(abcd):
            pa, pae = 0, 0
            if a > 0 and b>0 and _c > 0 and d>0:
                pa = d * b / _c
                pae = pa * ( (be/b)**2 + (ce/_c)**2 + (de/d)**2 )**(1/2)
                cl = abs((a - pa) / (ae**2 + pae**2)**(1/2))
                vals = [[d,de], [_c,ce], [b,be], [a,ae], [pa,pae], cl]
            else:
                vals = [[d,de], [_c,ce], [b,be], [a,ae], "", ""]
            
            hmax = max(hmax, a+ae, pa+pae)
        
            table_format[ith*len(mss) + ims] = [[sn, th if k=="signal" else ""], vals]
            
            ths.append(th - ( 0.02 * (len(mss)-1/2)/2) + 0.02 * ims)
            _as.append(a)
            aes.append(ae)
            pas.append(pa)
            paes.append(pae)

        gr_m = create_TGraph(ths, _as, [0]*len(ths), aes, ["Minimum dPhi","events"])
        gr_p = create_TGraph([th+0.01 for th in ths], pas, [0]*len(ths), paes, ["Minimum dPhi","events"])
        
        gr_m.SetLineColor(std_color_list[ims])
        gr_p.SetLineColor(std_color_list[2*len(mss) + ims])

        gr_m.SetMarkerColor(std_color_list[ims])
        gr_p.SetMarkerColor(std_color_list[2*len(mss) + ims])

        gr_m.SetMarkerSize(1.4)
        gr_p.SetMarkerSize(1.4)

        gr_m.SetLineWidth(3)#2*(2 * len(mss) - ims) + 2)
        gr_p.SetLineWidth(3)#2*(2 * len(mss) - ims - 1) + 2)

        gr_m.SetMarkerStyle(20) # filled o
        gr_p.SetMarkerStyle(47) # filled x

        # if not ms.is_mc:
        leg.AddEntry(gr_m,f"{ms.name}: A, measured","PE")
        leg.AddEntry(gr_p,f"           A, predicted","PE")
        # leg.append((k, std_color_list[ims]))
        
        # hmax = max(hmax, gr.GetMaximum())
        # hmin = min(hmin, #gr.GetMinimum(0))
        gr_m.SetMaximum(hmax*1.1)
        gr_m.SetMinimum(0)
        gr_p.SetMaximum(hmax*1.1)
        gr_p.SetMinimum(0)

        gr_m.Draw(("" if ims else "A") + "PE same")
        gr_p.Draw("PE same")
        gc.append(gr_m)
        gc.append(gr_p)

    leg.Draw()
    c.SetGrid()
    c.SetRightMargin(0.04)
    if not ROOT_BATCH:
        c.Draw()
    
    c.Print(f"{OUT_DIR}/closure.png")
    if VERBOSE:
        print(f"Saved 'closure.png'")

    abcd_table = Table(["D","C","B","A","Pred A", "CL"], "ABCD (X: dPhi, Y: dtSize)", rl_label=["","dPhi"])
    for label, vals in table_format:
        if label[1] != "":
            abcd_table.add_spacer()
        abcd_table.add_row(label, vals)
    abcd_table.print()

    ################################################################
    ################################################################
    ################################################################

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

        X_csc = np.array(
            [np.r_[np.ravel(ms_mc[feat]), np.ravel(ms_r3[feat])] for feat in feats_csc]
        ).T
        X_dt = np.array(
            [np.r_[np.ravel(ms_mc[feat]), np.ravel(ms_r3[feat])] for feat in feats_dt]
        ).T

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
            feat_table = Table(["Importance"], f"GBT {clfn}", 4, rl_label=["","Feature"])
            fts, wts = np.asarray(fts), clf.feature_importances_
            idxs = np.argsort(wts)[::-1]
            for ift, (ft, wt) in enumerate(zip(fts[idxs], wts[idxs])):
                feat_table.add_row([ift+1, ft], [wt])
            feat_table.print()

    ################################

    print("")
    alert("Making ROC Plots", form="-", c="g")

    bkgs, sigs, aucs = {}, {}, {}
    bkg_effs, sig_effs = {}, {}
    names = ["tag_dPhi", "tag_dEta", "tag_dR", "cscSize", "dtSize", "cscNStation10", "dtNStation10"]

    if TRAIN_BDT:
        names.extend(["GBT_CSC", "GBT_DT"])

    clf_table = Table(["Sig", "Sig Eff", "Bkg", "Bkg Eff", "S/√[B]", "AUC"], "Evaluating Discriminators", 3)

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
            y_true=np.r_[np.ones_like(vmc), np.zeros_like(vr3)], y_score=np.r_[vmc, vr3], sample_weight=np.r_[wmc, wr3]
        )

        sigs[name] = sigs[name][bkg_effs[name] > 0]
        bkgs[name] = bkgs[name][bkg_effs[name] > 0]
        sig_effs[name] = sig_effs[name][bkg_effs[name] > 0]
        bkg_effs[name] = bkg_effs[name][bkg_effs[name] > 0]

        # sig, bkg, sig_eff, bkg_eff, auc = bkgs[name], sig_effs[name], bkg_effs[name], aucs[name]

        auc, idx = aucs[name], np.argmax(sigs[name] / (np.sqrt(bkgs[name])))
        sig, sig_eff = sigs[name][idx], sig_effs[name][idx]
        bkg, bkg_eff = bkgs[name][idx], bkg_effs[name][idx]

        clf_table.add_row(name, [sig, sig_eff, bkg, bkg_eff, sig / np.sqrt(bkg), auc])

    clf_table.sort(lambda x: x[-2])
    clf_table.print()
    ################################

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

    ################################################################
    ################################################################
    ################################################################

    print("")
    alert("Finished train_bdt.py", form="=", c="c")
