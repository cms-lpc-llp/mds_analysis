import sys
import pathlib
from os import uname
import json

import numpy as np
import ROOT as rt
from math import ceil, floor
import numba as nb

from src.muon_system import MuonSystemAwkward
from src import CMS_lumi, tdrstyle
from src.helper_functions import alert, Table  # , canvas
from src.histo_utilities import create_TH1D, create_TH2D, create_TGraph, std_color_list

from histo_skim import (
    create_hists,
    ROOT_ERROR_LEVEL,
    DATA_VERSION,
    LUMI,
    FN_MC,
    FN_R3,
    T2_OUT_DIR,
    T2_DATA_DIR,
    LOCAL_OUT_DIR,
    LOCAL_DATA_DIR,
    OUT_DIR,
    DATA_DIR,
    gc,
    ABCD_DPHI,
    ABCD_DTSIZE,
    ABCD_INVERT_DPHI,
    TOP_MARGIN,
    BOT_MARGIN,
    CUTS,
)

# from sklearn.metrics import roc_auc_score  # , roc_curve
# from sklearn.model_selection import train_test_split
# from sklearn.preprocessing import StandardScaler
# from sklearn.ensemble import GradientBoostingClassifier  # , RandomForestClassifier


# **************************** #
if "TIER2" in DATA_DIR:
    OUT_DIR = f"{T2_OUT_DIR}/{OUT_DIR}"

    FN_MC = f"{T2_DATA_DIR}/MC_Summer22EE/v1/sixie/v{DATA_VERSION}/normalized/{FN_R3}.root"
    FN_R3 = f"{T2_DATA_DIR}/Data2022/v{DATA_VERSION}/normalized/{FN_R3}.root"

    FN_HLT_MC = f"{LOCAL_DATA_DIR}/../processed/mc_hlt569.root"  #! BROKEN
    FN_HLT_R3 = f"{LOCAL_DATA_DIR}/../processed/r3_hlt569.root"  #! BROKEN
else:
    OUT_DIR = f"{LOCAL_OUT_DIR}/{OUT_DIR}"

    FN_MC = f"{LOCAL_DATA_DIR}/{FN_MC}_v{DATA_VERSION}.root"
    FN_R3 = f"{LOCAL_DATA_DIR}/{FN_R3}_v{DATA_VERSION}.root"

    FN_HLT_MC = f"{LOCAL_DATA_DIR}/../processed/mc_hlt569.root"
    FN_HLT_R3 = f"{LOCAL_DATA_DIR}/../processed/r3_hlt569.root"

pathlib.Path(OUT_DIR).mkdir(parents=True, exist_ok=True)
# **************************** #
rt.gErrorIgnoreLevel = ROOT_ERROR_LEVEL
rt.gROOT.SetBatch(True)
tdrstyle.setTDRStyle()
CMS_lumi.writeExtraText = 0
# **************************** #

# ABCD_INVERT_DPHI = False
# ABCD_DPHI = 2.75
# ABCD_DTSIZE = 90

print(f"{ABCD_DPHI=}\n{ABCD_DTSIZE=}\n{ABCD_INVERT_DPHI=}")

N_EVENTS = -1
CUTS = [
    "match",
    "HLT",  #! Loading TBranch HLTDecision causes memory overflow so I use second precut TTree
    # "CSC&DT>0",
    "L1",
    "CSCIT",
    "DTIT",
    "MET",
    # "ME11/12",
    "MB1",
    "JET",
    # "MUON",
    # "BDT",
    "HALO",
    # "CSCSIZE",
    "DTSTN",
    "1CSC1DT",
    # "BLINDSR",
    # "DR",
    "DPHI",
]

if "HLT" in CUTS:
    ff_mc, ff_r3 = FN_MC, FN_HLT_R3
else:
    ff_mc, ff_r3 = FN_MC, FN_R3

# **** #


def weight_calc(llp_ct, new_ctau, old_ctau):
    if llp_ct.ndim > 1:
        llp_ct = np.array(np.sum(llp_ct, axis=1))
    source = np.exp(-1.0 * llp_ct / old_ctau) / old_ctau**2
    weight = 1.0 / new_ctau**2 * np.exp(-1.0 * llp_ct / new_ctau) / source
    return weight


from scipy.optimize import curve_fit, minimize, least_squares


def find_CI_limit(signal, observed, Z=1.65, do_verbose=False):
    """find_CI_limit
    Z=1.65 for 1-sided 95% CI"""

    def flimit(x):  # , xe):
        norm, x[0] = x[0], x[1] * x[3] / x[2]
        # norme, xe[0] = x[0], x[0] * np.sqrt(np.sum((xe[1:]/x[1:])**2))
        xe = np.maximum(np.sqrt(x), 1)
        norme, xe[0] = 1, x[0] * np.sqrt(np.sum((xe / x)[1:] ** 2))
        loss = ((norm * sig + x) - obs) / np.sqrt(obse**2 + (norm * sige) ** 2 + xe**2)
        # loss = ((norm * sig + x) - obs)/np.sqrt(obse**2)
        # loss[0] = ((norm * sig)/np.sqrt((norm*sige)**2 + xe**2) - Z)[0]
        # loss[0] = ((norm * sig)/np.sqrt(xe**2 + (norm*sige)**2) - Z)[0]
        loss_CI = ((norm * sig) / np.sqrt(2 * (xe**2) + (norm * sige) ** 2) - Z)[0]
        # loss_CI = ((norm * sig)/np.sqrt(xe**2) - Z)[0]
        return np.sum(loss[1:] ** 2) + loss_CI**2

    ##########
    sig, obs = np.asarray(signal).astype(float), np.asarray(observed).astype(float)
    sige, obse = np.sqrt(sig), np.sqrt(obs)
    obs[0] = np.divide(obs[1] * obs[3], obs[2], where=obs[2] > 0, out=np.zeros_like(obs[2]))
    obse[0] = obs[0] * np.sqrt(np.sum((np.divide(obse, obs, where=obs > 0, out=np.zeros_like(obs)))[1:]))
    # sige, obse = sig * 0.2, np.sqrt(obs)
    ##########
    p0 = np.array([np.sum(obs) / np.sum(sig), *obs[1:]])
    bounds = [(1e-21, max(7 * o, 1)) for o in obs]
    bounds[0] = (1e-12, np.inf)
    ##
    ftol = 1e-15
    res = minimize(flimit, p0, bounds=bounds, tol=ftol)
    # fit_func = lambda x: flimit(x, xe=np.maximum(1,np.sqrt(p0)))
    # for i in range(20):
    #     res = minimize(fit_func, p0, bounds=bounds, tol=ftol)

    p0 = res.x
    pe = np.diag(res.hess_inv.todense())  # np.ones_like(pars)
    pe = np.sqrt(max(1, abs(res.fun)) * ftol * pe)
    # fit_func = lambda x: flimit(x, xe=pe)

    norm, pbkg = res.x[0], res.x
    pbkge = np.diag(res.hess_inv.todense())  # np.ones_like(pars)
    pbkge = np.sqrt(max(1, abs(res.fun)) * ftol * pbkge)
    norme = pbkge[0]

    pbkge.setflags(write=True)
    pbkg[0] = pbkg[1] * pbkg[3] / pbkg[2]
    pbkge[0] = pbkg[0] * np.sqrt(np.sum((pbkge / pbkg)[1:] ** 2))

    if do_verbose:
        print("|   |    Obs | FitBkg |    Sig |")
        print("|---+--------+--------+--------+")
        print(f"| A | {obs[0]:6.2f} | {pbkg[0]:6.2f} | {norm*sig[0]:6.2f} |")
        print(f"| B | {obs[1]:6.2f} | {pbkg[1]:6.2f} | {norm*sig[1]:6.2f} |")
        print(f"| C | {obs[2]:6.2f} | {pbkg[2]:6.2f} | {norm*sig[2]:6.2f} |")
        print(f"| D | {obs[3]:6.2f} | {pbkg[3]:6.2f} | {norm*sig[3]:6.2f} |")
        print(f"Limit = {100*norm:.3f}%")
    # print(norm, norme, pbkg, pbkge)
    return norm, pbkg, norme, pbkge


def calc_abcd(weight, a_cond, b_cond, c_cond, d_cond, blind=False, g_cond=None, invert_dphi=None):
    if g_cond is not None:
        a_cond, b_cond = a_cond & g_cond, b_cond & g_cond
        c_cond, d_cond = c_cond & g_cond, d_cond & g_cond

    # if invert_dphi is None:
    #     flip = np.sum(a_cond) < a_cond.shape[0] / 2
    if invert_dphi:
        a_cond, d_cond = d_cond, a_cond
        b_cond, c_cond = c_cond, b_cond

    a = 0 if blind else np.sum(weight[a_cond])
    b = np.sum(weight[b_cond])
    c = np.sum(weight[c_cond])
    d = np.sum(weight[d_cond])
    ae, be, ce, de = np.sqrt(a), np.sqrt(b), np.sqrt(c), np.sqrt(d)
    pa, pae = 0, 0

    if b * c * d > 0:
        pa, pae = b * d / c, b * d / c * ((be / b) ** 2 + (ce / c) ** 2 + (de / d) ** 2) ** 0.5
    return (pa, a, b, c, d), (pae, ae, be, ce, de)


def smear(*args, fit=False, bootstrap=True, size=None):
    nev = args[0].shape[0]
    # if not bootstrap and size is None:
    #     return [arg if i != 2 else np.random.normal(arg, np.sqrt(arg)) for i, arg in enumerate(args)]

    if size is None:
        size = nev
    elif size < 1000:
        size = size * nev
    size = int(size)

    if bootstrap:
        idxs = np.random.randint(0, nev, size)
    else:
        idxs = np.arange(0, size) % nev

    out = []
    for i, arg in enumerate(args):
        # nbin = max(10, int(np.sqrt(arg.shape[0])))
        # nbin = floor(np.sqrt(arg.shape[0]))
        # _min, _max = np.quantile(arg, [0.1, 0.9])
        # _min, _max = np.quantile(arg, [0.2, 0.8])
        # _min, _max = np.quantile(arg, [0.16, 0.84]) # ~ +/- 1 sigma
        # _min, _max = np.quantile(arg, [0.05, 0.95]) # 90%
        _min, _max = np.quantile(arg, [0.00, 1.00])  # 90%
        _max = min(_max, 200)
        nbin = int(np.sqrt(np.sum((_min <= arg) & (arg <= _max))))
        # print(i, nbin, _min, _max)

        RAND_NUM = f"{np.random.randint(999999999)}"
        ft = rt.TF1("fit_" + RAND_NUM, "expo", _min, _max)
        hh = rt.TH1D("hh_" + RAND_NUM, "", nbin, _min, _max)

        if i == 0:
            out.append(arg[idxs] * nev / size)
        elif i == 1:
            _out = np.array([])
            if fit:
                # ft = rt.TF1("fit_df", "expo", 0.4, np.pi)
                # hh = rt.TH1D("", "", nbin, 0.4, np.pi)
                flip = np.mean(arg) > (np.pi - 0.4) * 3 / 4
                # print(flip, np.mean(arg), _min, _max)
                for x in arg:
                    hh.Fill(x)
                # hh.Fit("fit_df", "RL")
                hh.Fit("fit_" + RAND_NUM, "RL")
                ww = np.abs(ft.GetParameter(1))
                # print(i, ww)
                while len(_out) < size:
                    if flip:
                        _temp = np.pi - np.random.exponential(1 / ww, size)
                    else:
                        _temp = np.random.exponential(1 / ww, size) + 0.4
                    _out = np.append(_out, _temp[(0.4 < _temp) & (_temp < np.pi)])
            else:
                if bootstrap:
                    idxs = np.random.randint(0, nev, size)
                    # out.append(arg[idxs])
                    _out = arg[idxs]
                    # while len(_out) < size:
                    #     _temp = np.pi - np.random.exponential(1/1.76736e+00, size)
                    #     _out = np.append(_out, _temp[_temp >= 0.4])
            out.append(_out[:size])
        elif i == 2:
            if fit:
                # ft = rt.TF1("fit_sz", "expo", 50, _max)
                # hh = rt.TH1D("", "", nbin, 50, _max)
                for x in arg:
                    hh.Fill(x)
                # hh.Fit("fit_sz", "RL")
                hh.Fit("fit_" + RAND_NUM, "RL")
                ww = np.abs(ft.GetParameter(1))
                # print(i, ww)
                out.append(np.random.exponential(1 / ww, size) + 50)
            else:
                if bootstrap:
                    idxs = np.random.randint(0, nev, size)
                    # out.append(np.random.normal(arg[idxs], np.sqrt(arg[idxs])))
                    out.append(arg[idxs])
                    # out.append(np.random.exponential(1/1.21694e-02, size) + 50)
        else:
            out.append(arg[idxs])
        # fits.append(ft)

    return out


# **** #

ctaus = np.logspace(1, 7, 55)
# ctaus = np.logspace(np.log10(10), np.log10(30000), 38)
dmet = 1
max_met = 99999999
# for min_met in range(0, 300+dmet, dmet):
# mets = [ 0, 22, 34, 46, 1e99]

old_sbs = {}
bounds = {}
# mid_mets = [25, 50, 75, 100, 150, 200, 250, 300]
mid_mets = [x for x in range(20, 150, 1)]
# mid_mets = [50]
min_mets = [0] + [0] * len(mid_mets) + mid_mets
max_mets = [999] + mid_mets + [999] * len(mid_mets)

split_met_last_best = min(mid_mets)

for min_met, max_met in zip(min_mets, max_mets):  # [0] + mets[:-1], [1e99] + mets[1:]):
    print("\n\n" + 3 * (80 * "#" + "\n"))
    print(f"{min_met=}, {max_met=}")

    ms_mc = MuonSystemAwkward(ff_mc, tree_name="MuonSystem", name="Signal", nev=N_EVENTS, is_mc=True, lumi=LUMI)
    ms_r3_dtit = MuonSystemAwkward(
        ff_r3, tree_name="MuonSystem_HLT569", name="Data, DT IT", nev=N_EVENTS, is_mc=False, lumi=LUMI
    )
    # ms_r3_dtoot = MuonSystemAwkward(ff_r3, tree_name="MuonSystem_HLT569", name="Data, DT OOT", nev=N_EVENTS, is_mc=False, lumi=LUMI)

    # ms_mc = MuonSystemAwkward(FN_MC, name="Signal, Muon Vetoes", nev=N_EVENTS, is_mc=True, lumi=LUMI)
    # ms_r3_other = MuonSystemAwkward(FN_R3, name="Data, DT IT + Muon Vetoes", nev=N_EVENTS, is_mc=False, lumi=LUMI)
    # ms_mc = MuonSystemAwkward(FN_MC, name="Signal, Jet Vetoes", nev=N_EVENTS, is_mc=True, lumi=LUMI)
    # ms_r3_other = MuonSystemAwkward(FN_R3, name="Data, DT IT + Jet Vetoes", nev=N_EVENTS, is_mc=False, lumi=LUMI)

    mss = [ms_mc, ms_r3_dtit]
    # mss = [ms_mc, ms_r3_dtoot]
    # mss = [ms_mc, ms_r3_dtit, ms_r3_dtoot]

    # **** #
    print("")
    names = [ms.name for ms in mss]
    for ims, ms in enumerate(mss):
        # !!! #
        ms.cut = False
        # !!! #
        msn, is_mc = ms.name.lower(), ms.is_mc
        for cut in CUTS:
            match cut:
                case "match":
                    if is_mc:
                        ms.match_mc("csc,dt", has_clusters=True)
                case "HLT":
                    # if not is_mc:
                    #     ms.cut_hlt()
                    pass  #! Loading TBranch HLTDecision causes memory overflow so I use a precut TTree
                case "CSC&DT>0":
                    ms.f((ms["nCsc"] > 0) & (ms["nDt"] > 0))
                case "L1":
                    ms.cut_l1()
                case "CSCIT":
                    ms.cut_time("csc", invert="csc oot" in msn)
                case "DTIT":
                    ms.cut_time("dt", invert="dt oot" in msn)
                case "MET":
                    ms.f((min_met <= ms["met"]) & (ms["met"] < max_met))
                # case "ME11/12":
                #     ms.f(
                #         ms["cscNRechitChamberPlus11"]
                #         + ms["cscNRechitChamberMinus11"]
                #         + ms["cscNRechitChamberPlus12"]
                #         + ms["cscNRechitChamberMinus12"]
                #         == 0,
                #         "csc",
                #     )
                case "MB1":
                    ms.f(ms["dtNHitStation1"] == 0, "dt")
                case "JET":
                    # ms.cut_jet("csc,dt")
                    ms.cut_jet("csc,dt")  # , csc_pt=150, dt_pt=10)
                    # ms.cut_jet("dt")
                # case "MUON":
                #     ms.cut_muon("csc,dt")
                # case "BDT":
                #     pass
                case "HALO":
                    ms.cut_halo(invert=False)  #! HALO CUT
                # case "CSCSIZE":
                #     # ms.f(ms["cscSize"] > 250, "csc")
                #     ms.f(ms["cscSize"] > 200, "csc")
                # case "DTSIZE":
                #     ms.f(ms["dtSize"] < 200, "dt")
                # case "CSCSTN":
                #     ms.f(ms["cscSize"] > 250, "csc")
                case "DTSTN":
                    ms.f((ms["dtNStation10"] < 3) & ~((ms["dtNStation10"] == 2) & (ms["dtMaxStation"] == 4)), "dt")
                case "1CSC1DT":
                    ms.tag(tags="cscdt")
                # case "BLINDSR":
                #     if not is_mc:
                #         ms.f((ms["dtSize"] < ABCD_DTSIZE) | ((ms["tag_dPhi"] < ABCD_DPHI) ^ ABCD_INVERT_DPHI), "dt")
                # case "DR":
                #     ms.f(ms["tag_dR"] > 0.5)
                case "DPHI":
                    ms.f(ms["tag_dPhi"] > 0.4, invert=False)

        # !!! #
        ms.cut = True
        # !!! #

        ms.colors = [std_color_list[ims], std_color_list[len(mss) + ims]]
        print(f"{ms.name:>{max([len(n) for n in names])}} : {ms.count():,.1f} events")

    # **** #
    weights = [np.asarray(ms["weight"]) for ms in mss]
    colors = [ms.colors[0] for ms in mss]
    blind_conds = [
        np.ones_like(ms["tag_dPhi"], dtype=bool)
        if ms.is_mc
        else (ms["dtSize"][:, 0] < ABCD_DTSIZE) | ((ms["tag_dPhi"] < ABCD_DPHI) ^ ABCD_INVERT_DPHI)
        for ms in mss
    ]

    # if sum([ms.count() < 3 for ms in mss]):
    if sum([ms.count() < 1 for ms in mss]):
        print(f"LESS THAN 1 EVENTS FOUND -- ({min_met=})")
        # print(f"LESS THAN 1 EVENTS FOUND -- EXITING ({min_met=})")

    # *** #
    split_met = min_met if min_met else max_met
    if sum([ms.count(use_weight=False) < 3 for ms in mss]):
        if split_met in old_sbs:  # np.sqrt(np.sum(weights[0])) * np.sum(weights[1]) < old_sbs[split_met]:
            alert(f"NOT ENOUGH in {min_met}, {max_met} USING BOUNDS FOUND FROM {split_met_last_best}", "r")
            ABCD_DPHI_S2B, ABCD_DTSIZE_S2B, ABCD_FLIP_S2B = bounds[split_met_last_best]
            if np.median(mss[0]["tag_dPhi"]) < np.pi / 2:
                ABCD_FLIP_S2B = True
                alert(f"MOST SIG NEAR dPHI=0 -- SETTING {ABCD_FLIP_S2B=}", "g")
            # elif min_met >= 100:
            elif min_met >= 75:
                ABCD_FLIP_S2B = True
                alert(f"FORCING ABCD FLIP TO {ABCD_FLIP_S2B=}", "g")
            else:
                alert(f"KEEPING LAST ABCD FLIP STATE {ABCD_FLIP_S2B=}", "g")
        else:
            alert(f"SOMETHING WRONG {min_met}, {max_met}", "r")
            # alert(f"TOO FEW EVENTS AND NOT IN SBS {min_met}, {max_met}", "r")
            break

    else:
        split_met_last_best = split_met
        # if split_met == min_met:
        #     min_mets.append(0)
        #     max_mets.append(min_met)
        old_sbs[split_met] = {"sqrt[|s|.|b|]": np.sqrt(np.sum(weights[0])) * np.sum(weights[1])}

        NEV_SMEAR = max([10 * w.shape[0] for w in weights])
        # NEV_SMEAR = max([100_000] + [w.shape[0] for w in weights])
        print(f"{NEV_SMEAR=:,.0f}")

        nth_sz, nth_df = (200 - 50) // 1 + 1, (315 - 40) // 1 + 1  # 31, 19
        size_ths, dphi_ths = np.linspace(50, 200, nth_sz), np.linspace(0.40, 3.15, nth_df)
        dth_sz, dth_df = size_ths[1] - size_ths[0], dphi_ths[1] - dphi_ths[0]
        thresh_info = {}  # {ms.name : [] for ms in mss}

        # sz_conds, df_conds, fp_conds = [], [], [] # what is this for???

        print(np.median(mss[0]["tag_dPhi"]), np.median(mss[1]["tag_dPhi"]))
        flip_dfc = np.median(mss[0]["tag_dPhi"]) < np.median(mss[1]["tag_dPhi"])
        for ims, ms in enumerate(mss):
            wt, df, sz = np.asarray(ms["weight"]), np.asarray(ms["tag_dPhi"]), np.asarray(ms["dtSize"][:, 0])
            # if not ms.is_mc:
            #     blind_idx =  (sz < np.max(size_ths)) | (df < np.max(dphi_ths))
            #     wt, df, sz = wt[blind_idx], df[blind_idx], sz[blind_idx]
            # wt = wt[np.random.randint(0, sz.shape[0], NEV_SMEAR)] * sz.shape[0] / NEV_SMEAR
            # df = df[np.random.randint(0, sz.shape[0], NEV_SMEAR)]
            # sz = sz[np.random.randint(0, sz.shape[0], NEV_SMEAR)]
            # sz = np.random.normal(sz, np.sqrt(sz))
            df_cond, sz_cond = df[None, :] > dphi_ths[:, None], sz[None, :] > size_ths[:, None]

            df_cond = ~df_cond if flip_dfc else df_cond

            _has_data = np.array(
                [
                    (
                        ((df_cond[idf]) & (~sz_cond[isz])).any()
                        and ((~df_cond[idf]) & (~sz_cond[isz])).any()
                        and ((~df_cond[idf]) & (sz_cond[isz])).any()
                    )
                    for idf in range(nth_df)
                    for isz in range(nth_sz)
                ]
            )

            # _has_data = np.array(
            #     [
            #         (~sz_cond[isz] | (~df_cond[idf]) & (sz_cond[isz])).any()
            #         for idf in range(nth_df)
            #         for isz in range(nth_sz)
            #     ]
            # )

            if NEV_SMEAR > 0:
                wt, df, sz = smear(wt, df, sz, fit=ims > 0, bootstrap=True, size=NEV_SMEAR)

            df_cond, sz_cond = df[None, :] > dphi_ths[:, None], sz[None, :] > size_ths[:, None]
            df_cond = ~df_cond if flip_dfc else df_cond

            _abcd = np.array(
                [
                    [
                        np.sum(wt[(df_cond[idf]) & (sz_cond[isz])]),
                        np.sum(wt[(df_cond[idf]) & (~sz_cond[isz])]),
                        np.sum(wt[(~df_cond[idf]) & (~sz_cond[isz])]),
                        np.sum(wt[(~df_cond[idf]) & (sz_cond[isz])]),
                    ]
                    for idf in range(nth_df)
                    for isz in range(nth_sz)
                ]
            ).T
            _err = np.sqrt(_abcd[0])

            # sz_conds.append(np.sum(sz_cond, 1).astype(float))
            # df_conds.append(np.sum(df_cond, 1).astype(float))
            # fp_conds.append(fp_cond.astype(float)) # what are these lines for???

            # _abcd[_abcd<0.5] = 0
            # _abcd[_abcd<1] = 0
            # _abcd[_abcd<2] = 0
            # _abcd = np.round(_abcd)
            # _abcd[_abcd<1] = np.round(_abcd[_abcd<1])
            # _abcd[_abcd<2] = np.ceil(_abcd[_abcd<2])
            # _abcd[_abcd<1] = np.sqrt(_abcd[_abcd<1])
            # _abcd[_abcd<1] = np.round(np.sqrt(_abcd[_abcd<1]))

            if not ms.is_mc:
                _abcd[0] = np.divide(_abcd[1] * _abcd[3], _abcd[2], out=np.zeros_like(_abcd[0]), where=_abcd[2] != 0)

                _idxs = (_abcd[1] * _abcd[2] * _abcd[3]) > 0
                _err[_idxs] = _abcd[0, _idxs] * np.sqrt(np.sum(1 / _abcd[1:, _idxs]))
                _err[~_idxs] = 1

            thresh_info[ms.name] = np.array(
                [
                    [
                        dphi_ths[idf],
                        size_ths[isz],
                        flip_dfc,
                        *_abcd[:, isz + idf * nth_sz],
                        _err[isz + idf * nth_sz],
                        _has_data[isz + idf * nth_sz],
                    ]
                    for idf in range(nth_df)
                    for isz in range(nth_sz)
                ]
            )

        # *** #
        # S / Sqrt[ B ]
        bins, log, norm = (
            [
                nth_df,
                dphi_ths[0] - dth_df / 2,
                dphi_ths[-1] + dth_df / 2,
                nth_sz,
                size_ths[0] - dth_sz / 2,
                size_ths[-1] + dth_sz / 2,
            ],
            False,
            False,
        )
        axis_titles = [
            "|#Delta#phi_{CSC,DT}| Boundary",
            "DT_{NRechit} Boundary",
            f"{'fraction of ' if norm else ''}" + "S / #sqrt{B}",
        ]

        # *** #
        _evts = [thresh_info[ms.name][:, 3] for ms in mss]
        _hds = [thresh_info[ms.name][:, 8] for ms in mss]
        # *** #
        idxs = (_evts[0] > 0) & (_evts[1] > 0) & (_hds[0] > 0) & (_hds[1] > 0)
        # idxs = _evts[1] > 0
        if np.sum(idxs) == 0:
            print(f"NO THRESHOLDS FOR 1 EVENT IN B, C, & D -- EXITING ({min_met=})")
            break
        _evts = [evts[idxs] for evts in _evts]
        _values = [thresh_info[ms.name][idxs, :2] for ms in mss]
        _flips = [thresh_info[ms.name][idxs, 2] for ms in mss]
        _abcd = [thresh_info[ms.name][idxs, 3:7] for ms in mss]
        _valid_abcd = [(thresh_info[ms.name][idxs, 4:7] > 0).all() for ms in mss]
        _errs = [thresh_info[ms.name][idxs, 7] for ms in mss]
        # *** #
        _s2b = _evts[0] / np.sqrt(_evts[1])

        if np.sum(_s2b) == 0:
            print(f"NO _S2B -- EXITING ({min_met=})")
            continue
        _lim = [[1.0, 1.0, 1.0, 1.0] for sig, obs in zip(_abcd[0], _abcd[1])]
        # _lim = [find_CI_limit(sig, obs, do_verbose=False) for sig, obs in zip(_abcd[0], _abcd[1])]
        _lim, _lime = np.asarray([l[0] for l in _lim]), np.asarray([l[2] for l in _lim])
        # *** #
        # idxs = _lim > 1e-4
        # _evts = [x[idxs] for x in _evts]
        # _values = [x[idxs] for x in _values]
        # _abcd = [x[idxs] for x in _abcd]
        # _errs = [x[idxs] for x in _errs]
        # _s2b = _s2b[idxs]
        # _lim = _lim[idxs]
        # *** #

        _v, _n, _c, _w = (
            [_values[0], _values[0]],
            ["S2B", "BR Limit (CI 95%)"],
            [std_color_list[0], std_color_list[1]],
            [_s2b, _lim],
        )

        lls = []

        for ih in range(len(_n)):  # , hh in enumerate(hhs):
            v, sig, bkg = _v[0], _evts[0], _evts[1]
            n, c, w = _n[ih], _c[ih], _w[ih]
            vv = _valid_abcd[0] * _valid_abcd[1]

            if ih == 0:
                idx_s2b = np.argmax(w * ((sig >= 1) & (bkg >= 1)))
                # if ((sig >= 1) & (bkg >= 1)).any():
                #     idx_s2b = np.argmax(w * ((sig >= 1) & (bkg >= 1)))
                # else:  # use old lim bnds
                #     alert(f"USING PREVIOUS {idx_s2b=}, (right now {min_met=}, {max_met=})")
                #     _values[1][idx_s2b] = [ABCD_DPHI_S2B, ABCD_DTSIZE_S2B]
                #     _flips[0][idx_s2b] = ABCD_FLIP_S2B
                #     pass
                idx = idx_s2b
                m = "X"
            elif ih == 1:
                if ((sig >= 1) & (bkg >= 1)).any():
                    idx_lim = np.argmin(w + np.max(w) * ((sig < 1) | (bkg < 1)))
                else:  # use old lim bnds
                    pass
                idx = idx_lim
                m = "O"

            df, sz, s, b = v[idx][0], v[idx][1], sig[idx], bkg[idx]
            fp = _flips[0][idx]
            s2b, lim = _s2b[idx], _lim[idx]
            lls.append(
                (
                    c,
                    df,
                    sz,
                    f"{m} | dPhi = {df:4.2f}, dtSize = {sz:5.1f}, dfFlip = {fp} -- Limit = {100*lim:5.3f}%, S/Sqrt[B] = {s2b:6.2f}, S = {s:6.2f}, B = {b:5.2f}",
                )
            )

        for ih, (n, c) in enumerate(zip(_n, _c)):
            for cl, df, sz, ll in lls:
                if ih == 0:
                    print(ll)

        ABCD_DPHI_S2B, ABCD_DTSIZE_S2B = _values[1][idx_s2b]
        ABCD_FLIP_S2B = _flips[0][idx_s2b]
        # ABCD_DPHI_S2B, ABCD_DTSIZE_S2B = _values[1][idx_lim]
        bounds[split_met] = (ABCD_DPHI_S2B, ABCD_DTSIZE_S2B, ABCD_FLIP_S2B)

    # *** #

    limits, elimits = [], []
    s2bs = []
    for ct in ctaus:
        # if ct != 1000:
        #     continue
        abcd_sig, abcd_obs = [], []
        for ims, ms in enumerate(mss):
            dtsize, dphi = ms["dtSize"][:, 0], ms["tag_dPhi"]
            sz_cond, df_cond = dtsize > ABCD_DTSIZE_S2B, dphi > ABCD_DPHI_S2B
            df_cond = ~df_cond if ABCD_FLIP_S2B else df_cond
            a_cond = sz_cond & df_cond
            b_cond = ~sz_cond & df_cond
            c_cond = ~sz_cond & ~df_cond
            d_cond = sz_cond & ~df_cond

            wt = np.asarray(ms["weight"])
            if ims == 0 and ct != 1000:  # ms["dt_match_gLLP_index"][:,0]
                wt *= weight_calc(np.asarray(ms["gLLP_ctau"][:, 0]), ct / 10, 1000 / 10)

            vals, errs = calc_abcd(wt, a_cond, b_cond, c_cond, d_cond, blind="dt it" in ms.name.lower())

            if ims == 0:
                abcd_sig = list(vals[1:5])
                asig = vals[1]
            elif ims == 1:
                abcd_obs = list(vals[1:5])
                abcd_obs[0] = vals[0]
                s2bs.append(asig / np.sqrt(vals[0]) if vals[0] else 999)

                if ct == 1000:
                    tsig, tbkg = np.sum(mss[0]["weight"]), np.sum(mss[1]["weight"])
                    asig, abkg = abcd_sig[0], abcd_obs[0]
                    ts2b = tsig / np.sqrt(tbkg) if tbkg else 999
                    as2b = asig / np.sqrt(abkg) if abkg else 999
                    with open("info.txt", "a") as ff_info:
                        print(f"min_met = {min_met}, max_met = {max_met}")
                        print(f"DPHI = {ABCD_DPHI_S2B}")
                        print(f"SIZE = {ABCD_DTSIZE_S2B}")
                        print(f"FLIP = {ABCD_FLIP_S2B}")
                        print(f"Total: Sig = {tsig}, Bkg = {np.sum(mss[1]['weight'])}, S/Sqrt[B] = {ts2b}")
                        print(f" SR A: Sig = {asig}, Bkg = {abkg}, S/Sqrt[B] = {as2b}")

                        ff_info.write(f"min_met = {min_met}, max_met = {max_met}\n")
                        ff_info.write(f"DPHI = {ABCD_DPHI_S2B}\n")
                        ff_info.write(f"SIZE = {ABCD_DTSIZE_S2B}\n")
                        ff_info.write(f"FLIP = {ABCD_FLIP_S2B}\n")
                        ff_info.write(f"Total: Sig = {tsig}, Bkg = {tbkg}, S/Sqrt[B] = {ts2b}\n")
                        ff_info.write(f" SR A: Sig = {asig}, Bkg = {asig}, S/Sqrt[B] = {as2b}\n\n")

                # print(f"\ta = " +(f"{vals[1]:7.3f} ± {errs[1]:7.3f}" if ms.is_mc else f"{vals[0]:7.3f} ± {errs[0]:7.3f} (b*d/c)"))
                # print(f"\tb = {vals[2]:7.3f} ± {errs[2]:7.3f}")
                # print(f"\tc = {vals[3]:7.3f} ± {errs[3]:7.3f}")
                # print(f"\td = {vals[4]:7.3f} ± {errs[4]:7.3f}")
        if (
            len(abcd_sig) == 0 or len(abcd_obs) == 0
        ):  # or (np.asarray(abcd_sig) == 0).any() or (np.asarray(abcd_obs) == 0).any():
            print("???")
            continue

        s2b = s2bs[-1]
        # lv, pv, le, pe = find_CI_limit(abcd_sig, abcd_obs, do_verbose=False)
        lv, pv, le, pe = 0, 0, 0, 0
        limits.append(lv)
        elimits.append(le)

        # print(f"{ct:0>7.0f} - {asig:7.2f}, {vals[0]:5.2f}, {asig/np.sqrt(vals[0]):7.2f}")

        norm = 10 ** int(np.log10(5 / sum(abcd_sig)))
        fdc = f"datacard_ct{ct:0>8.0f}_minMet{min_met:0>3.0f}_maxMet{max_met:0>3.0f}.txt"

        # print(f"{min_met=}, {(abcd_sig/np.sqrt(abcd_obs))[0]}")
        # print(f"combine -M AsymptoticLimits --freezeParam norm --setParameters norm={norm} {fdc}")
        str_dc = f"""# norm {norm}
# CTAU {ct}
# ABCD_DPHI {ABCD_DPHI_S2B}
# ABCD_SIZE {ABCD_DTSIZE_S2B}
# ABCD_FLIP {ABCD_FLIP_S2B}
# S2B,LIMIT {s2b} {lv} {le}
imax 4
jmax 1
kmax *
shapes * * FAKE
---------------
---------------
bin \t chA \t chB \t chC \t chD
observation \t {abcd_obs[0]:.3f} \t {abcd_obs[1]:.3f} \t {abcd_obs[2]:.3f} \t {abcd_obs[3]:.3f}
------------------------------
bin \t chA \t chA \t chB \t chB \t chC \t chC \t chD \t chD
process \t sig \t bkg \t sig \t bkg \t sig \t bkg \t sig \t bkg
process \t 0 \t 1 \t 0 \t 1 \t 0 \t 1 \t 0 \t 1
rate \t {abcd_sig[0]:.6e} \t 1 \t {abcd_sig[1]:.6e} \t 1 \t {abcd_sig[2]:.6e} \t 1 \t {abcd_sig[3]:.6e} \t 1
------------------------------
single_A \t rateParam \t chA \t bkg \t (@0*@2/@1) \t single_B,single_C,single_D
single_B \t rateParam \t chB \t bkg \t {abcd_obs[1]:.3f} \t [0,{max(7, 7*abcd_obs[1]):.3f}]
single_C \t rateParam \t chC \t bkg \t {abcd_obs[2]:.3f} \t [0,{max(7, 7*abcd_obs[2]):.3f}]
single_D \t rateParam \t chD \t bkg \t {abcd_obs[3]:.3f} \t [0,{max(7, 7*abcd_obs[3]):.3f}]
norm rateParam * sig 1
single_lumi \t lnN \t 1.2 \t - \t 1.2 \t - \t 1.2 \t - \t 1.2 \t -
"""
        with open(f"data/datacard_scan/{fdc}", "w") as fdc:
            fdc.write(str_dc)
