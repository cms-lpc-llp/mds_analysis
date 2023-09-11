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
################################
##         Parameters         ##
################################
SAVE_STAT = "dtOOT"
BLIND_TYPE = "DTOOT"  # , "DTSIZE", "DYNAMIC", "NONE"

N_EVENTS = -1
TEST_SIZE = 0.8

ABCD_DPHI = 1.5
ABCD_DTSIZE = 200

################################
TRAIN_BDT = False

CUTS = [
    "match",
    #!!!!!!!#
    # "HLT", # TODO: Load HLTDecision without overflowing memory
    #!!!!!!!#
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

SECTIONS = [
    "EFFICIENCY",
    "CSC",
    "DT",
    "EVENT",
    "2D",
    "ABCD",
    "ROC",
]
################################
#! Paths are absolute so I know *exactly* what is being read/written
OUT_DIR = "reports/weekly/2023-09-07"
T2_OUT_DIR = "/storage/af/user/psimmerl/LLP/mdc_analysis"  # os.getcwd()
LOCAL_OUT_DIR = "/home/psimmerl/LLP/mdc_analysis"  # os.getcwd()
################################
DATA_VERSION = "6"
LUMI = 23.02 * 1000

T2_DATA_DIR = "/storage/cms/store/user/christiw/displacedJetMuonAnalyzer/Run3/V1p19"
LOCAL_DATA_DIR = "/home/psimmerl/LLP/mdc_analysis/data/raw" # os.getcwd() + "/data/raw"
DATA_DIR = "TIER2" if "caltech" in uname()[1] else "LOCAL"

FN_MC = "ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted"
FN_R3 = "DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi"

# with open("default_bins.json", "r") as f: #TODO: this should be absolute for consitency
#     DEFAULT_BINS = json.load(f)
################################
VERBOSE = True
ROOT_BATCH = True
ROOT_ERROR_LEVEL = 1001 #rt.kInfo + 1
BOT_MARGIN, TOP_MARGIN = 0.025, 0.1

gc = []
################################

def create_hists(values, bins, names, titles, weights=1, colors=None, styles=None, log=None, norm=False, canvas=None):
    if not isinstance(names, (list, tuple)):
        values, names, = [ values ], [ names ]
    if isinstance(values[0], (list, tuple)):
        hist_type = "2D"
    else:
        hist_type = "1D"

    if colors is None:
        colors = std_color_list

    if not isinstance(weights, (list, tuple)):
        weights = [ weights ] * len(names)
    if not isinstance(colors, (list, tuple)):
        colors = [ colors ] * len(names)
    if not isinstance(styles, (list, tuple)):
        styles = [ styles ] * len(names)

    logx, logy, logz = False, False, False
    if isinstance(log, (list, tuple)):
        logx = log[0]
        if len(log) > 1:
            logy = log[1]
        if len(log) > 2:
            logz = log[2]
    else:
        if hist_type == "1D":
            logy = log
        else:
            logz = log

    hists = []

    if canvas is None:
        cn = "c" + str(np.random.randint(999999999))
        canvas = rt.TCanvas(cn, cn, 800, 800)

    lat = rt.TLatex()
    lat.SetTextAlign(11)
    lat.SetTextSize(0.04)
    hmin, hmax = 1e9, -1

    for n, v, w, c, s in zip(names, values, weights, colors, styles):
        if hist_type == "1D":
            hist = create_TH1D(v, name=n, axis_title=titles, binning=bins, weights=w)
        else:
            hist = create_TH2D(v, name=n, axis_title=titles, binning=bins, weights=w)
        # lat.SetTextColor(color)
        # lat.DrawLatexNDC((0.94 - 0.17)/max(5, len(names))*(len(names) - i - 1), 0.92, text)
        hist.SetLineColor(c)
        if s is not None:
            hist.SetLineStyle(c)


        if norm and hist.Integral() > 0:
            hist.Scale(1.0 / hist.Integral())

        hmax = max(hmax, hist.GetMaximum())
        hmin = min(hmin, hist.GetMinimum(0))
        hists.append(hist)

    if hmax/hmin > 1000:# clean this up
        hmin, hmax = np.log10( hmin if hmin==1 else 1e-1), np.log10(hmax)
        hmin, hmax = hmin - BOT_MARGIN*(hmax-hmin), hmax + TOP_MARGIN*(hmax-hmin)
        hmin, hmax = 10**(hmin), 10**(hmax)
        log=True
    else:
        hmin, hmax = hmin - BOT_MARGIN*(hmax-hmin), hmax + TOP_MARGIN*(hmax-hmin)
        hmin = 0 if hmin < 50 else hmin
        log=False

    if hist_type == "1D": # clean this up
        logy = log
    else:
        logz = log

    for hist in hists:
        hist.SetMaximum(hmax)
        hist.SetMinimum(hmin)
        hist.Draw("hist same")

    if logx:
        canvas.SetLogx()
    if logy:
        canvas.SetLogx()
    if logz:
        canvas.SetLogx()

    canvas.SetGrid()
    canvas.SetRightMargin(0.04)

    if not ROOT_BATCH:
        canvas.Draw()

    return canvas, hists


################################################################
################################################################
################################################################

def main():
    print("")
    alert("Building MuonSystem(s)", form="-", c="g")
    ms_mc = MuonSystemAwkward(FN_MC, name="Signal", nev=N_EVENTS, is_mc=True, lumi=LUMI)
    ms_r3 = MuonSystemAwkward(FN_R3, name="Data", nev=N_EVENTS, is_mc=False, lumi=LUMI)
    mss = [ms_mc, ms_r3]

    if "standard" in sys.argv[2]:
        CUTS.remove("HALO")
    elif "neither" in sys.argv[2]:
        CUTS.remove("HALO")
        CUTS.remove("DPHI")
    elif "halo" in sys.argv[2]:
        CUTS.remove("DPHI")
        ms_nocut = MuonSystemAwkward(FN_R3, name="NoHaloCut", nev=N_EVENTS, is_mc=False, lumi=LUMI)
        mss.append(ms_nocut)
    elif "dphi" in sys.argv[2]:
        ms_nocut = MuonSystemAwkward(FN_R3, name="NoDPhiCut", nev=N_EVENTS, is_mc=False, lumi=LUMI)
        mss.append(ms_nocut)
    elif "both" in sys.argv[2]:
        pass
        # ms_cut = MuonSystemAwkward(ff_r3, name="halo&dphi", nev=N_EVENTS, is_mc=False, lumi=LUMI)
        # mss.append(ms_cut)
    else:
        pass

    ################################################################

    for ims, ms in enumerate(mss):
        if ms.colors is None:
            ms.colors = [std_color_list[ims], std_color_list[len(mss) + ims]]
        #!!! TURNING MuonSystem CUTS OFF !!!!#
        # cuts are not applied when you get MuonSystem data
        # but they are 'AND'ed with sel_evt, sel_csc, and sel_dt.
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

                    #!!! TURNING EFFICIENCY ON !!!!#
                    if "EFFICIENCY" in SECTIONS:
                        ms.efficiency_denom = "eff"
                        ms.ms_read["eff_evt"] = ms.ms_read["sel_evt"]
                        ms.ms_read["eff_csc"] = ms.ms_read["sel_csc"]
                        ms.ms_read["eff_dt"] = ms.ms_read["sel_dt"]
                        ms.efficiency = 1
                    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

                case "HLT":
                    #!!!!!!!!!!!#
                    ms.cut_hlt()
                    #!!!!!!!!!!!#
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
                    if "NoHaloCut" in msn:
                        was_cut = False
                    else:
                        ms.cut_halo(invert="Halo" in msn)
                case "1CSC1DT":
                    ms.tag(tags="cscdt")
                case "DPHI":
                    if "NoDPhiCut" in msn:
                        was_cut = False
                    else:
                        ms.f(0.4 < ms["tag_dPhi"], invert="dphi" in msn.lower())
                case default:
                    raise ValueError(f"Cut '{cut}' not recognized!", c="r")
            if was_cut:
                counts[ims][0] = ms.count()
                if cut in ("match", "HLT"): # skip efficiency
                    counts[ims] = counts[ims][0]
                elif "EFFICIENCY" in SECTIONS:
                    counts[ims][1] = f"({ms.efficiency:.{eff_decimals}f})"

        cut_label = ("DTIT/OOT" if cut =="DTIT" and BLIND_TYPE=="DTOOT" else cut)
        cut_table.add_row(cut_label, counts)
        if cut in ("match", ):
            cut_table.add_spacer()
        if VERBOSE:
            print(f"Finished cut '{cut}'")#, counts)
    cut_table.print()

    ################################################################

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
        print(ms.ms_read.keys())
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
        # Some may *need* 2 tag
        if "EVENT" in SECTIONS:
            ms["cscR-dtR"] = ms["cscR"] - ms["dtR"]
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#


    # ################################################################
    # ################################################################
    # ################################################################

    if "CSC" in SECTIONS:
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
        # names, colors, weights = [(ms.name, ms.color[0], ms["weight"]) for ms in mss]
        # for xl, _bins in zip(xls, bins):
        #     values = [ms[xl] for ms in mss]
        #     if xl in ("cscEta", "cscPhi", "cscZ"):
        #         values = [np.abs(v) for v in values]

        #     canvas, hists = create_hists(values, _bins, names, ["xl","events"], weights, colors, styles=None, log=None, norm=False, canvas=None)
        #     canvas.Print(f"{OUT_DIR}/{xl}_{SAVE_STAT}.png")
        #     if VERBOSE:
        #         print(f"Saved '{xl}_{SAVE_STAT}.png'")

        # for ix, xl in enumerate(xls):
        #     cn = "c" + str(np.random.randint(999999999))
        #     c = rt.TCanvas(cn, cn, 800, 800)
        #     leg = []
        #     lat = rt.TLatex()
        #     lat.SetTextAlign(11)
        #     lat.SetTextSize(0.04)
        #     hmin, hmax = 1e9, -1

        #     for ims, ms in enumerate(mss):
        #         k, var, weight = ms.name, ms[xl], ms["weight"]
        #         if xl in ("cscEta", "cscPhi", "cscZ"):
        #             var = np.abs(var)

        #         # if "tag_dPhi" in xl or "halo" in k:
        #         #     cond = ms["tag_dPhi"] > 0.0
        #         # else:
        #         #     cond = ms["tag_dPhi"] > 0.5
        #         cond = ms["tag_dPhi"] > -1

        #         if "size" in BLIND_TYPE:
        #             cond = cond & (ms["dtSize"][:, 0] < 100)

        #         if not ms.is_mc:
        #             var, weight = var[cond], weight[cond]

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

        #     if xl not in ("cscPhi",):
        #         c.SetLogy()
        #     c.SetGrid()
        #     c.SetRightMargin(0.04)

        #     if not ROOT_BATCH:
        #         c.Draw()

        #     c.Print(f"{OUT_DIR}/{xl}_{SAVE_STAT}.png")
        #     if VERBOSE:
        #         print(f"Saved '{xl}_{SAVE_STAT}.png'")

    ################################

    if "DT" in SECTIONS:
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
            hmin, hmax = 1e9, -1

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

        for ix, xl in enumerate(xls):
            cn = "c" + str(np.random.randint(999999999))
            c = rt.TCanvas(cn, cn, 800, 800)
            leg = []
            lat = rt.TLatex()
            lat.SetTextAlign(11)
            lat.SetTextSize(0.04)
            hmin, hmax = 1e9, -1

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

    if "2D" in SECTIONS or "ABCD" in SECTIONS:
        print("")
        alert("Making 2D Plots", form="-", c="g")

        xyls, bins = [], []
        if "2D" in SECTIONS:
            xyls.extend([
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
            ])
            bins.extend([
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
            ])

        if "ABCD" in SECTIONS:
            n_abcd = 0
            for th in np.arange(0, 1.0+0.1, step=0.1):
                if "DPHI" in CUTS and th < 0.4:
                    continue
                xyls.append(("tag_dPhi", "dtSize"))
                bins.append([th, ABCD_DPHI, np.pi, 50, ABCD_DTSIZE, 500, 3, 3])
                n_abcd += 1

            abcds = {ms.name : [] for ms in mss}

        if "2D" in SECTIONS or "ABCD" in SECTIONS:
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

    if "ABCD" in SECTIONS:
        print("")
        alert("Performing ABCD Calculations", form="-", c="g")

        abcd_table_sig = Table(
            cols=["D","C","B","A","Pred A", "CL"],
            header="ABCD (X: dPhi, Y: dtSize)",
            rv_decimals=2,
            rl_label=["","dPhi"]
        )
        abcd_table_bkg = Table(
            cols=["D","C","B","A","Pred A", "CL"],
            header="ABCD (X: dPhi, Y: dtSize)",
            rv_decimals=2,
            rl_label=["","dPhi"]
        )

        table_rows_sig = [ ["", None] ] * n_abcd * len(mss)
        table_rows_bkg = [ ["", None] ] * n_abcd * len(mss)

        abcd_graphs = {}
        graph_axis_labels = ["Minimum |#Delta#phi_{CSC,DT}|","#scale[0.5]{#int} L"]
        hmin_bkg, hmax_bkg = 1e9, -1
        hmin_sig, hmax_sig = 1e9, -1

        for ims, ms in enumerate(mss):
            k = ms.name
            abcd = abcds[k]

            ths, _as, aes, pas, paes = [], [], [], [], []
            for ith, (th, a, ae, b, be, _c, ce, d, de) in enumerate(abcd):
                pa, pae = 0, 0
                if a > 0 and b > 0 and _c > 0 and d > 0:
                    pa = d * b / _c
                    pae = pa * ( (be/b)**2 + (ce/_c)**2 + (de/d)**2 )**(1/2)
                    cl = abs((a - pa) / (ae**2 + pae**2)**(1/2))
                    vals = [[d,de], [_c,ce], [b,be], [a,ae], [pa,pae], cl]
                else:
                    vals = [[d,de], [_c,ce], [b,be], [a,ae], "", ""]

                am, pam = (a-ae if a-ae>0 else 1), (pa-pae if pa-pae>0 else 1)
                ap, pap = a+ae, pa+pae

                hmin_sig = min(hmin_sig, am, pam)
                hmax_sig = max(hmax_sig, ap, pap)
                if ims > 0:
                    hmin_bkg = min(hmin_bkg, am, pam)
                    hmax_bkg = max(hmax_bkg, ap, pap)

                table_rows_sig[ith*len(mss) + ims] = [[k, th if k=="Signal" else ""], vals]
                table_rows_bkg[ith*len(mss) + ims] = [[k, th if k=="Data" else ""], vals]

                ths.append(th - ( 0.04 * (len(mss)-1/2)/2) + 0.04 * ims)
                _as.append(a)
                aes.append(ae)
                pas.append(pa)
                paes.append(pae)

            abcd_graphs[k] = [
                create_TGraph(ths, _as, [0.04]*len(ths), aes, graph_axis_labels),
                create_TGraph([th+0.02 for th in ths], pas, [0]*len(ths), paes, graph_axis_labels)
            ]
            gc.extend(abcd_graphs[k])

        ### Graphs ##

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
            leg = rt.TLegend(0.17, 0.91 - 0.04*(2*len(_mss)), 0.21 + 0.016*(msn_max+11), 0.94)
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

            if _max/_min > 1000:
                c.SetLogy()
                _min, _max = np.log10( _min if _min==1 else 1e-1), np.log10(_max)
                _min, _max = _min - BOT_MARGIN*(_max-_min), _max + TOP_MARGIN*(_max-_min)
                _min, _max = 10**(_min), 10**(_max)
            else:
                _min, _max = _min - BOT_MARGIN*(_max-_min), _max + TOP_MARGIN*(_max-_min)
                _min = 0 if _min < 50 else _min

            for ims, ms in enumerate(_mss):
                for igr, gr in enumerate(abcd_graphs[ms.name]):
                    gr.SetMinimum(_min)
                    gr.SetMaximum(_max)

                    gr.SetLineColor(ms.colors[igr])
                    gr.SetLineWidth(4)
                    gr.SetLineStyle(9 if igr else 1) # '--' if meas else '-'

                    gr.SetMarkerColor(ms.colors[igr])
                    gr.SetMarkerSize(2.5 if igr else 2)
                    gr.SetMarkerStyle(47 if igr else 20) # x if meas else o

                    lab = f"{ms.name:>{msn_max}}: "+("predicted" if igr else "measured")
                    le = leg.AddEntry(gr,lab,"PE")
                    # le.SetTextColor(ms.colors[igr])

                    gr.Draw(("AF" if igr+ims == 0 else "") + "PE same")


            leg.Draw()
            if not ROOT_BATCH:
                c.Draw()

            c.Print(f"{OUT_DIR}/closure{label}.png")
            if VERBOSE:
                print(f"Saved 'closure{label}.png'")

        ### Tables ##

        for table, row in [(abcd_table_sig, table_rows_sig), (abcd_table_bkg, table_rows_bkg)]:
            for irow, (label, vals) in enumerate(row):
                if label[1] != "" and irow:
                    table.add_spacer()
                table.add_row(label, vals)
            table.print()

        #############

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

    if "ROCS" in SECTIONS:
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
                y_true=np.r_[np.ones_like(vmc), np.zeros_like(vr3)],
                y_score=np.r_[vmc, vr3], sample_weight=np.r_[wmc, wr3]
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

if __name__ == "__main__":
    alert("Starting train_bdt.py", form="=", c="c")

    if len(sys.argv) > 1:
        N_EVENTS = int(sys.argv[1])
    if len(sys.argv) > 2:
        OUT_DIR = OUT_DIR + f"/{sys.argv[2]}"
    if len(sys.argv) > 3:
        BLIND_TYPE = sys.argv[3]

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
        FN_MC = f"{T2_DATA_DIR}/MC_Summer22EE/v1/sixie/v{DATA_VERSION}/normalized/{FN_R3}.root"
        FN_R3 = f"{T2_DATA_DIR}/Data2022/v{DATA_VERSION}/normalized/{FN_R3}.root"
    else:
        OUT_DIR = f"{LOCAL_OUT_DIR}/{OUT_DIR}"
        FN_MC = f"{LOCAL_DATA_DIR}/{FN_MC}_v{DATA_VERSION}.root"
        FN_R3 = f"{LOCAL_DATA_DIR}/{FN_R3}_v{DATA_VERSION}.root"
    pathlib.Path(OUT_DIR).mkdir(parents=True, exist_ok=True)
    ################################
    rt.gErrorIgnoreLevel = ROOT_ERROR_LEVEL
    if ROOT_BATCH:
        rt.gROOT.SetBatch(ROOT_BATCH)
    tdrstyle.setTDRStyle()
    CMS_lumi.writeExtraText = 0
    ################################

    main()

    print("")
    alert("Finished train_bdt.py", form="=", c="c")
