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

from train_bdt import (
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
    ABCD_DTSIZE,
    ABCD_DPHI,
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

    FN_HLT_MC = f"{LOCAL_DATA_DIR}/../processed/mc_hlt569.rdf" #! BROKEN
    FN_HLT_R3 = f"{LOCAL_DATA_DIR}/../processed/r3_hlt569.rdf" #! BROKEN
else:
    OUT_DIR = f"{LOCAL_OUT_DIR}/{OUT_DIR}"

    FN_MC = f"{LOCAL_DATA_DIR}/{FN_MC}_v{DATA_VERSION}.root"
    FN_R3 = f"{LOCAL_DATA_DIR}/{FN_R3}_v{DATA_VERSION}.root"

    FN_HLT_MC = f"{LOCAL_DATA_DIR}/../processed/mc_hlt569.rdf"
    FN_HLT_R3 = f"{LOCAL_DATA_DIR}/../processed/r3_hlt569.rdf"

pathlib.Path(OUT_DIR).mkdir(parents=True, exist_ok=True)
# **************************** #
rt.gErrorIgnoreLevel = ROOT_ERROR_LEVEL
rt.gROOT.SetBatch(True)
tdrstyle.setTDRStyle()
CMS_lumi.writeExtraText = 0
# **************************** #

def smear(*args, fit=False, bootstrap=True, size=None):
    nev = args[0].shape[0]
    # if not bootstrap and size is None:
    #     return [arg if i != 2 else np.random.normal(arg, np.sqrt(arg)) for i, arg in enumerate(args)]

    if size is None:
        size = nev
    elif size < 1000:
        size = size*nev
    size = int(size)

    if bootstrap:
        idxs = np.random.randint(0, nev, size)
    else:
        idxs = np.arange(0, size) % nev

    out = []
    for i, arg in enumerate(args):
        nbin = max(10, int(np.sqrt(arg.shape[0])))
        if i == 0:
            out.append(arg[idxs] * nev/size)
        elif i == 1:
            _out = np.array([])
            if fit:
                ff = rt.TF1("fit_df","expo",0.4,np.pi)
                hh = rt.TH1D("","",nbin,0.4,np.pi)
                for x in arg:
                  hh.Fill(x)
                hh.Fit("fit_df", "RL")
                ww = np.abs(ff.GetParameter(1))
                print(ww)
                # out.append(np.random.uniform(0.4, np.pi, size))
                while len(_out) < size:
                    _temp = np.random.exponential(1/ww, size) + 0.4
                    _out = np.append(_out, _temp[_temp <= np.pi])
            else:
                # out.append(arg[idxs])
                _out = arg[idxs]
                # while len(_out) < size:
                #     _temp = np.pi - np.random.exponential(1/1.76736e+00, size)
                #     _out = np.append(_out, _temp[_temp >= 0.4])
            out.append(_out[:size])
        elif i == 2:
            if fit:
                _max_sz = np.quantile(arg, 0.95)
                ff = rt.TF1("fit_sz","expo",50,_max_sz)
                hh = rt.TH1D("","",nbin,50,_max_sz)
                for x in arg:   
                    hh.Fill(x)
                hh.Fit("fit_sz", "RL")
                ww = np.abs(ff.GetParameter(1))
                print(ww)
                # out.append(np.random.normal(arg[idxs], np.sqrt(arg[idxs])))
                out.append(np.random.exponential(1/ww, size) + 50)
            else:
                out.append(np.random.normal(arg[idxs], np.sqrt(arg[idxs])))
                # out.append(arg[idxs])
                # out.append(np.random.exponential(1/1.21694e-02, size) + 50)
        else:
            out.append(arg[idxs])

    return out
def calc_abcd(weight, a_cond, b_cond, c_cond, d_cond, blind=False, g_cond=None):
    if g_cond is not None:
        a_cond, b_cond = a_cond & g_cond, b_cond & g_cond
        c_cond, d_cond = c_cond & g_cond, d_cond & g_cond
    a = 0 if blind else np.sum(weight[a_cond])
    b = np.sum(weight[b_cond])
    c = np.sum(weight[c_cond])
    d = np.sum(weight[d_cond])
    ae, be, ce, de = np.sqrt(a), np.sqrt(b), np.sqrt(c), np.sqrt(d)
    pa, pae = 0, 0

    if b * c * d > 0:
        pa, pae = b * d / c, b * d / c * ((be / b) ** 2 + (ce / c) ** 2 + (de / d) ** 2) ** 0.5
    return (pa,a,b,c,d), (pae,ae,be,ce,de)

# *** #

N_EVENTS = -1
CUTS = [
    "match",
    "HLT", #! Loading TBranch HLTDecision causes memory overflow so I use second precut TTree
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
    "CSCSIZE",
    "DTSTN",
    "1CSC1DT",
    # "BLINDSR",
    # "DR",
    "DPHI",
]

ff_mc, ff_r3 = FN_MC, FN_HLT_R3

NEV_SMEAR = 100_000  # max([w.shape[0] for w in weights])
nth_sz, nth_df = (200-60)//1 + 1, (300-150)//5 + 1 # 31, 19
size_ths, dphi_ths = np.linspace(60, 200, nth_sz), np.linspace(1.5, 3.0, nth_df)
dth_sz, dth_df = size_ths[1]-size_ths[0], dphi_ths[1]-dphi_ths[0]

# **** #

header = ["s2b", "s", "b", "df", "sz", "met", "csz", "cjp", "djp", "pid"]
with open("data/datacard_scan.csv", "w") as fout:
    fout.write(" ".join(header)+"\n")

params = []
param_ids = set()

# for i in range(1000):
while True:
    irun = len(param_ids)
    print("-"*50)
    _met = max(0, np.random.uniform(-50, 200))
    # _halo
    _csc_size = np.random.uniform(150, 500)
    _csc_jet_pt = np.random.uniform(0.1, 175)
    _dt_jet_pt = np.random.uniform(0.1, 175)
    _csc_size = _csc_size if _csc_size>200 else 0
    _csc_jet_pt = _csc_jet_pt if _csc_jet_pt>10 else 10
    _dt_jet_pt = _dt_jet_pt if _dt_jet_pt>10 else 10
    _csc_jet_pt = _csc_jet_pt if _csc_jet_pt<150 else -1
    _dt_jet_pt = _dt_jet_pt if _dt_jet_pt<150 else -1

    if irun == 0:
        _met = 0
        _csc_size = 250
        _csc_jet_pt = 10
        _dt_jet_pt = 10

    # _min_dphi
    print(_met, _csc_size, _csc_jet_pt, _dt_jet_pt)


    ms_mc = MuonSystemAwkward(ff_mc, tree_name="MuonSystem", name="Signal", nev=N_EVENTS, is_mc=True, lumi=LUMI)
    ms_r3_dtit = MuonSystemAwkward(ff_r3, tree_name="MuonSystem_HLT569", name="Data, DT IT", nev=N_EVENTS, is_mc=False, lumi=LUMI)
    mss = [ms_mc, ms_r3_dtit]

    # **** #
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
                    pass  #! Loading TBranch HLTDecision causes memory overflow so I use second precut TTree
                case "CSC&DT>0":
                    ms.f((ms["nCsc"] > 0) & (ms["nDt"] > 0))
                case "L1":
                    ms.cut_l1()
                case "CSCIT":
                    ms.cut_time("csc", invert="csc oot" in msn)
                case "DTIT":
                    ms.cut_time("dt", invert="dt oot" in msn)
                case "MET":
                    ms.f(ms["met"] > _met)
                case "ME11/12":
                    ms.f(
                        ms["cscNRechitChamberPlus11"]
                        + ms["cscNRechitChamberMinus11"]
                        + ms["cscNRechitChamberPlus12"]
                        + ms["cscNRechitChamberMinus12"]
                        == 0,
                        "csc",
                    )
                case "MB1":
                    ms.f(ms["dtNHitStation1"] == 0, "dt")
                case "JET":
                    cats = []
                    if _csc_jet_pt>0:
                        cats.append("csc")
                    if _dt_jet_pt>0:
                        cats.append("dt")
                    ms.cut_jet(",".join(cats), csc_pt=_csc_jet_pt, dt_pt=_dt_jet_pt)
                case "MUON":
                    ms.cut_muon("csc,dt")
                case "BDT":
                    pass
                case "HALO":
                    ms.cut_halo(invert=False)  #! HALO CUT
                case "CSCSIZE":
                    ms.f(ms["cscSize"] > _csc_size, "csc")
                case "DTSIZE":
                    ms.f(ms["dtSize"] < 200, "dt")
                case "CSCSTN":
                    pass
                case "DTSTN":
                    ms.f((ms["dtNStation10"] < 3) & ~((ms["dtNStation10"] == 2) & (ms["dtMaxStation"] == 4)), "dt")
                case "1CSC1DT":
                    ms.tag(tags="cscdt")
                case "BLINDSR":
                    if not is_mc:
                        ms.f((ms["dtSize"] < ABCD_DTSIZE) | (ms["tag_dPhi"] < ABCD_DPHI), "dt")
                case "DR":
                    pass
                case "DPHI":
                    ms.f(0.4 < ms["tag_dPhi"], invert=False)

        # !!! #
        ms.cut = True
        # !!! #

    # **** #
    print([ms.count() for ms in mss])
    if sum([ms.count() < 3 for ms in mss]):
        continue

    # **** #
    try:
        thresh_info = {}#{ms.name : [] for ms in mss}
        for ims, ms in enumerate(mss):
            wt, df, sz = np.asarray(ms["weight"]), np.asarray(ms["tag_dPhi"]), np.asarray(ms["dtSize"][:,0])
            # if NEV_SMEAR>0:
            #     wt, df, sz = smear(wt, df, sz, fit=ims>0, bootstrap=True, size=NEV_SMEAR)
            
            df_cond, sz_cond = df[None,:]>dphi_ths[:,None], sz[None,:]>size_ths[:,None]
            _abcd = np.array([[np.sum(wt[(df_cond[idf])&(sz_cond[isz])]),np.sum(wt[(df_cond[idf])&(~sz_cond[isz])]),np.sum(wt[(~df_cond[idf])&(~sz_cond[isz])]),np.sum(wt[(~df_cond[idf])&(sz_cond[isz])])] for idf in range(nth_df) for isz in range(nth_sz)]).T
            
            # _abcd = np.round(_abcd) 
            # _abcd[_abcd<1] = np.round(np.sqrt(_abcd[_abcd<1])) 
            # _abcd[_abcd<1] = np.round(_abcd[_abcd<1]) 
            _abcd[_abcd<1] = 0

            if not ms.is_mc:
                _abcd[0] = np.divide(_abcd[1]*_abcd[3], _abcd[2], out=np.zeros_like(_abcd[0]), where=_abcd[2]!=0)
                abcd_obs =  np.array([ _abcd[:,isz + idf*nth_sz] for idf in range(nth_df) for isz in range(nth_sz) ])
            else:
                abcd_sig =  np.array([ _abcd[:,isz + idf*nth_sz] for idf in range(nth_df) for isz in range(nth_sz) ])

            thresh_info[ms.name] = np.array([ [dphi_ths[idf], size_ths[isz], _abcd[0,isz + idf*nth_sz]] for idf in range(nth_df) for isz in range(nth_sz) ])


        # S / Sqrt[ B ]
        _evts = [thresh_info[ms.name][:,2] for ms in mss]
        _values = [thresh_info[ms.name][evts>0,:2] for ms, evts in zip(mss, _evts)]
        _s2bs = [_evts[0][evts>0] / np.sqrt(evts[evts>0]) for evts in _evts]

        sig, bkg = _evts[0][_evts[1]>0], _evts[1][_evts[1]>0]
        idx_s2b = np.argmax(_s2bs[1])# * (sig > 300))# * (bkg < 1))

        df, sz = _values[1][idx_s2b]
        s2b, s, b = _s2bs[1][idx_s2b], sig[idx_s2b], bkg[idx_s2b]
        print(f"dPhi = {df:4.2f}, dtSize = {sz:5.1f} -- S/Sqrt[B] = {s2b:6.2f}, S = {s:6.2f}, B = {b:5.2f}")

        par_id = 0
        while par_id in param_ids:
            par_id = np.random.randint(0, 1_000_000)
        params.append([s2b, s, b, df, sz, _met, _csc_size, _csc_jet_pt, _dt_jet_pt, par_id])
        param_ids.add(par_id)

        with open("data/datacard_scan.csv", "a") as fout:
            fout.write(" ".join(str(x) for x in params[-1])+"\n")

        abcd_sig = abcd_sig[_evts[1]>0][idx_s2b]
        abcd_obs = abcd_obs[_evts[1]>0][idx_s2b]

        with open(f"data/datacard_scan/datacard_{par_id:0>6}.txt", "w") as fout:
            fout.write(f"""# norm {10**int(np.log10(5/sum(abcd_sig)))}
# ABCD_DPHI {df}
# ABCD_SIZE {sz}
# params {' '.join(str(x) for x in params[-1])}
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
rate \t {abcd_sig[0]:.3f} \t 1 \t {abcd_sig[1]:.3f} \t 1 \t {abcd_sig[2]:.3f} \t 1 \t {abcd_sig[3]:.3f} \t 1
------------------------------
single_A \t rateParam \t chA \t bkg \t (@0*@2/@1) \t single_B,single_C,single_D
single_B \t rateParam \t chB \t bkg \t {abcd_obs[1]:.3f} \t [0,{7*abcd_obs[1]:.3f}]
single_C \t rateParam \t chC \t bkg \t {abcd_obs[2]:.3f} \t [0,{7*abcd_obs[2]:.3f}]
single_D \t rateParam \t chD \t bkg \t {abcd_obs[3]:.3f} \t [0,{7*abcd_obs[3]:.3f}]
norm rateParam * sig 1
single_lumi \t lnN \t 1.2 \t - \t 1.2 \t - \t 1.2 \t - \t 1.2 \t -
""")

        print("")
        params.sort(key=lambda x:x[0])
        for s2b, s, b, df, sz, met, csz, cjp, djp, pid in params[-10:]:
            print(f"S/Sqrt[B]={s2b:6.2f}, S={s:6.2f}, B={b:5.2f} | dPhi={df:4.2f}, dtSize={sz:3.0f} | met={met:6.2f}, cscSize={csz:6.2f}, cscJetPt={cjp:6.2f}, dtJetPt={djp:6.2f}")
        print("\n")
    except:
        continue
