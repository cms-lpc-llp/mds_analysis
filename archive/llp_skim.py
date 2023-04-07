# llp_skim.py

import ROOT as rt

# import sys

import numpy as np

# import pandas as pd

rt.gInterpreter.ProcessLine('#include "src/LLPEvent.h"')

data_year = "full"
years = ["Run3_Fall22"]
fname = "data/raw/ggH_HToSSTobbbb_MH-125_MS-15_ctau-1000_TuneCP5_13TeV-powheg-pythia8_59740pb_weighted.root"
tname = "MuonSystem"
pi = rt.TMath.Pi()
ff = rt.TFile(fname, "READ")
tt = ff.Get(tname)

# ntuple = rt.TNtuple('llp_skim', 'met;deltaPhi')

gc = []
evs_dict = {
    "raw": {},
    "rawM": {},
}

value_templates = {
    "abs_cscRechitCluster_match_gLLP_decay_z": (100, 0, 1000),
    "abs_cscRechitCluster_match_gLLP_decay_r": (100, 0, 1000),
    "abs_gLLP_decay_vertex_z": (100, 0, 1000),
    "abs_gLLP_decay_vertex_r": (100, 0, 1000),
    "nDTRechits": (100, 0, 1000),
    "nCscRechits": (100, 0, 1000),
    "nDtRechitClusters": (6, -0.5, 5.5),
    "nCscRechitClusters": (6, -0.5, 5.5),
    "cscRechitCluster_match_gLLP_AtLeastOne": (6, -0.5, 5.5),
    "met": (100, 0, 800),
    "delta_metPhi_cscRechitClusterPhi_fixed": (100, -np.pi, np.pi),
    "cscRechitClusterEta": (100, -2, 2),
    "cscRechitClusterSize": (100, 0, 3000),
    "cscRechitClusterTime": (100, -10, 10),
    "cscRechitClusterNStation10": (4, 0.5, 4.5),
    "cscRechitClusterMaxStation": (9, -4.5, 4.5),
    "cscRechitClusterAvgStation10": (100, -4.5, 4.5),
    "cscRechitCluster_match_gLLP_e": (100, -10, 10),
    "cscRechitCluster_match_gLLP_pt": (100, -10, 10),
    "cscRechitCluster_match_gLLP_eta": (100, -10, 10),
    "cscRechitCluster_match_gLLP_phi": (100, -10, 10),
    "cscRechitCluster_match_gLLP_ctau": (100, 0, 3.25),
    "cscRechitCluster_match_gLLP_decay_z": (100, -10, 10),
    "cscRechitCluster_match_gLLP_decay_r": (100, -10, 10),
    "gLLP_decay_vertex_z": (100, -10, 10),
    "gLLP_decay_vertex_r": (100, -10, 10),
    "gLLP_decay_z_efficiency": (100, -10, 10),
    "gLLP_decay_r_efficiency": (100, -10, 10),
    "gHiggsE": (100, -10, 10),
    "gHiggsPt": (100, -10, 10),
    "gHiggsEta": (100, -10, 10),
    "gHiggsPhi": (100, -10, 10),
}
comment_templates = {}
histos = {}


def fill_histo(comment, xname, xval, yname=None, yval=[], rep="c", xbins=None, ybins=None):
    xval, yval = [x for x in xval], [y for y in yval]

    if yname:  # 2D
        if xbins is None and ybins is None:
            if comment in comment_templates:
                templates = comment_templates[comment]
                if (xname, yname) in templates:
                    xbins, ybins = templates[(xname, yname)]
                elif (yname, xname) in templates:
                    xbins, ybins = templates[(yname, xname)]
                else:
                    if xname in templates:
                        xbins = templates[xbins]
                    if yname in templates:
                        ybins = templates[ybins]
            if xbins is None:
                xbins = value_templates[xname]
            if ybins is None:
                ybins = value_templates[yname]

            if (xname, yname) in value_templates:
                xbins, ybins = value_templates[(xname, yname)]
            elif (yname, xname) in value_templates:
                xbins, ybins = value_templates[(yname, xname)]

        name = f"{xname}_{yname}_{comment}"
        title = f"{name};{comment} (x:{xname}, y:{yname});{xname};{yname};count"
        if name not in histos:
            histos[name] = rt.TH2F(name, title, *xbins, *ybins)

        if rep == "c":
            xval2, yval2 = [], []
            while xval:
                xv = xval.pop(0)
                if isinstance(xv, (list, np.ndarray)):
                    if xv.size:
                        xval.extend(xv)
                else:
                    xval2.append(xv)

            while yval:
                yv = yval.pop(0)
                if isinstance(yv, (list, np.ndarray)):
                    if yv.size:
                        yval.extend(yv)
                else:
                    yval2.append(yv)

            for xv in xval2:
                for yv in yval2:
                    histos[name].Fill(xv, yv)
        # elif rep == 'o':
        #   for xv, yv in xval, yv:
        #       histos[name].Fill(xv, yv)

    else:  # 1D
        if xbins is None and ybins is None:
            if comment in comment_templates:
                templates = comment_templates[comment]
                if xname in templates:
                    xbins = templates[xbins]
            if xbins is None:
                xbins = value_templates[xname]

        name = f"{xname}_{comment}"
        title = f"{name};{comment} (x:{xname});{xname};count"

        if name not in histos:
            histos[name] = rt.TH1F(name, title, *xbins)

        if rep == "c":
            xval2 = []
            while xval:
                xv = xval.pop(0)
                if isinstance(xv, (list, np.ndarray)):
                    if xv.size:
                        xval.extend(xv)
                else:
                    xval2.append(xv)

            for xv in xval2:
                histos[name].Fill(xv)
        # elif rep == 'o':
        #   for xv, yv in xval, yv:
        #       histos[name].Fill(xv, yv)
        #       filled = True


def fill_all_histos(comment, values):
    # 1D Histograms
    # hists = [
    # ]

    for valx in value_templates:
        if valx in values:
            fill_histo(comment, valx, values[valx])
            for valy in value_templates:
                if valy in values:
                    fill_histo(comment, valx, values[valx], valy, values[valy])


############
### CUTS ###
############

met_min, cscRechitClusterSize_min = 200, 200
jetPt_max, jetEta_max = 50, 2.4
cbid_eta_max = [1.9, 1.1, 1.6, 1.6, 1.8]

############
############
############


def build_ntuple_row(ev, cut_val=None, col=None, tdict=None, more_cols_dict=None):
    cut_val, col = cut_val if cut_val else [], col if col else []

    if not isinstance(col, list):
        cut_val, col = [cut_val], [col]
    met = ev.met

    nCscRechits = ev.nCscRechits
    nCscRechitClusters = ev.nCscRechitClusters

    cscRechitClusterSize = [x for x in ev.cscRechitClusterSize]
    cscRechitClusterTime = [x for x in ev.cscRechitClusterTime]
    cscRechitClusterEta = [x for x in ev.cscRechitClusterEta]
    cscRechitClusterNStation10 = [x for x in ev.cscRechitClusterNStation10]
    cscRechitClusterMaxStation = [x for x in ev.cscRechitClusterMaxStation]
    cscRechitClusterAvgStation10 = [x for x in ev.cscRechitClusterAvgStation10]
    cscRechitCluster_match_gLLP_e = [x for x in ev.cscRechitCluster_match_gLLP_e]
    cscRechitCluster_match_gLLP_pt = [x for x in ev.cscRechitCluster_match_gLLP_pt]
    cscRechitCluster_match_gLLP_eta = [x for x in ev.cscRechitCluster_match_gLLP_eta]
    cscRechitCluster_match_gLLP_phi = [x for x in ev.cscRechitCluster_match_gLLP_phi]

    delta_metPhi_cscRechitClusterPhi = [ev.metPhi - x for x in ev.cscRechitClusterPhi]
    delta_metPhi_cscRechitClusterPhi_fixed = [x + 2 * pi * (-1) ** (x > pi) for x in delta_metPhi_cscRechitClusterPhi]

    abs_cscRechitCluster_match_gLLP_decay_z = [abs(x) for x in ev.cscRechitCluster_match_gLLP_decay_z]
    abs_cscRechitCluster_match_gLLP_decay_r = [abs(x) for x in ev.cscRechitCluster_match_gLLP_decay_r]
    abs_gLLP_decay_vertex_z = [abs(x) for x in ev.gLLP_decay_vertex_z]
    abs_gLLP_decay_vertex_r = [abs(x) for x in ev.gLLP_decay_vertex_r]

    nDTRechits = ev.nDTRechits
    nDtRechitClusters = ev.nDtRechitClusters

    gHiggsE = ev.gHiggsE
    gHiggsPt = ev.gHiggsPt
    gHiggsEta = ev.gHiggsEta
    gHiggsPhi = ev.gHiggsPhi

    tdict_temp = {
        "met": met,
        "nCscRechits": nCscRechits,
        "nCscRechitClusters": nCscRechitClusters,
        "cscRechitClusterSize": cscRechitClusterSize,
        "cscRechitClusterTime": cscRechitClusterTime,
        "cscRechitClusterEta": cscRechitClusterEta,
        "cscRechitClusterNStation10": cscRechitClusterNStation10,
        "cscRechitClusterMaxStation": cscRechitClusterMaxStation,
        "cscRechitClusterAvgStation10": cscRechitClusterAvgStation10,
        "cscRechitCluster_match_gLLP_e": cscRechitCluster_match_gLLP_e,
        "cscRechitCluster_match_gLLP_pt": cscRechitCluster_match_gLLP_pt,
        "cscRechitCluster_match_gLLP_eta": cscRechitCluster_match_gLLP_eta,
        "cscRechitCluster_match_gLLP_phi": cscRechitCluster_match_gLLP_phi,
        "delta_metPhi_cscRechitClusterPhi_fixed": delta_metPhi_cscRechitClusterPhi_fixed,
        "abs_cscRechitCluster_match_gLLP_decay_z": abs_cscRechitCluster_match_gLLP_decay_z,
        "abs_cscRechitCluster_match_gLLP_decay_r": abs_cscRechitCluster_match_gLLP_decay_r,
        "abs_gLLP_decay_vertex_z": abs_gLLP_decay_vertex_z,
        "abs_gLLP_decay_vertex_r": abs_gLLP_decay_vertex_r,
        "nDTRechits": nDTRechits,
        "nDtRechitClusters": nDtRechitClusters,
        "gHiggsE": gHiggsE,
        "gHiggsPt": gHiggsPt,
        "gHiggsEta": gHiggsEta,
        "gHiggsPhi": gHiggsPhi,
    }

    if isinstance(more_cols_dict, dict):
        tdict_temp.update(more_cols_dict)

    for k in tdict_temp:
        if not isinstance(tdict_temp[k], (list, np.ndarray)):
            tdict_temp[k] = [tdict_temp[k]]
        tdict_temp[k] = np.asarray(tdict_temp[k])

    for cv, cc in zip(cut_val, col):
        for k in tdict_temp:
            if cv and cc in k and isinstance(tdict_temp[k], list):
                tdict_temp[k] = tdict_temp[k][cv]

    # if len(cscRechitClusterSize) == 0:
    #   return tdict

    for k, v in tdict_temp.items():
        if k in tdict:
            tdict[k].append(v)
        else:
            tdict[k] = [v]

    return tdict


# def fill_histograms(ev)

############
############
############

tot = 0

for iev, ev in enumerate(tt):
    # if iev > 100:
    #   break
    if iev % 10 == 0 and iev != 0:
        print(iev, len(evs_dict["raw"]["met"]))  # , len(evs_dict['rawM']['met']))

    tot += len(ev.cscRechitClusterSize) > 0

    abs_cscRechitCluster_match_gLLP_decay_z = [abs(x) for x in ev.cscRechitCluster_match_gLLP_decay_z]
    abs_cscRechitCluster_match_gLLP_decay_r = [abs(x) for x in ev.cscRechitCluster_match_gLLP_decay_r]
    abs_gLLP_decay_vertex_z = [abs(x) for x in ev.gLLP_decay_vertex_z]
    abs_gLLP_decay_vertex_r = [abs(x) for x in ev.gLLP_decay_vertex_r]

    delta_metPhi_cscRechitClusterPhi = [ev.metPhi - x for x in ev.cscRechitClusterPhi]
    delta_metPhi_cscRechitClusterPhi_fixed = [x + 2 * pi * (-1) ** (x > pi) for x in delta_metPhi_cscRechitClusterPhi]

    ################
    ### Matching ###
    ################
    pass_csc_match = [x for x in ev.cscRechitCluster_match_gLLP]

    ###########
    ### MET ###
    ###########

    pass_met = ev.met > met_min
    # if not pass_met:
    #   continue

    pass_met50 = ev.met > 50
    pass_met100 = ev.met > 100
    pass_met150 = ev.met > 150

    #########################
    ### Rechit Size > 200 ###
    #########################

    pass_cscRechitClusterSize_idxs = [x > cscRechitClusterSize_min for x in ev.cscRechitClusterSize]
    pass_cscRechitClusterSize_AtLeastOneAbove200 = sum([x > 200 for x in ev.cscRechitClusterSize]) > 0
    # if sum(pass_cscRechitClusterSize_idxs) < 1:
    #   continue

    #############################
    ### 1 jet and lepton cuts ###
    #############################

    pass_jet = 0
    pass_lep = ev.nLeptons == 0
    for jetPt, jetEta in zip(ev.jetPt, ev.jetEta):
        pass_jet += (jetPt > 50) * (abs(jetEta) < 2.4)
    # if not (pass_jet and pass_lep):
    #   continue

    ############
    ### CBID ###
    ############

    pass_CBID_idxs = []

    for eta, nstation10, avgstation10 in zip(
        ev.cscRechitClusterEta, ev.cscRechitClusterNStation10, ev.cscRechitClusterAvgStation10
    ):
        pass_CBID_idxs.append(abs(eta) < cbid_eta_max[int(avgstation10) if nstation10 > 1 else 0])
    # if sum(pass_CBID_idxs) < 1:
    #   continue

    ####################################
    ### dR Cluster / Jet / Muon Cuts ###
    ####################################

    # ! Not quite?
    pass_CC_dR_idxs = [np.inf] * ev.nCscRechitClusters
    pass_CJ_dR_idxs = [np.inf] * ev.nCscRechitClusters
    pass_CM_dR_idxs = [np.inf] * ev.nCscRechitClusters
    for i, (csc_eta, csc_phi) in enumerate(zip(ev.cscRechitClusterEta, ev.cscRechitClusterPhi)):
        for j, (csc2_eta, csc2_phi) in enumerate(zip(ev.cscRechitClusterEta, ev.cscRechitClusterPhi)):
            pass_CC_dR_idxs[i] = min(pass_CC_dR_idxs[i], np.sqrt((csc_eta - csc2_eta) ** 2 + (csc_phi - csc2_phi) ** 2))

        for j, (jet_eta, jet_phi) in enumerate(zip(ev.cscRechitClusterEta, ev.cscRechitClusterPhi)):
            pass_CJ_dR_idxs[i] = min(pass_CJ_dR_idxs[i], np.sqrt((csc_eta - jet_eta) ** 2 + (csc_phi - jet_phi) ** 2))

        for j, (lep_eta, lep_phi, lep_pdg) in enumerate(zip(ev.lepEta, ev.lepPhi, ev.lepPdgId)):
            if abs(lep_pdg) == 13:
                pass_CM_dR_idxs[i] = min(
                    pass_CM_dR_idxs[i], np.sqrt((csc_eta - lep_eta) ** 2 + (csc_phi - lep_phi) ** 2)
                )

    pass_CC_dR_idxs = [x > 0.6 for x in pass_CC_dR_idxs]
    pass_CJ_dR_idxs = [x > 0.4 for x in pass_CJ_dR_idxs]
    pass_CM_dR_idxs = [x > 0.4 for x in pass_CM_dR_idxs]

    #########################
    ### Muon Chamber cuts ###
    #########################

    pass_cluster_rechits_ME11_ME12 = [
        (M11 + M12 + P11 + P12) == 0
        for M11, M12, P11, P12 in zip(
            ev.cscRechitClusterNRechitChamberMinus11,
            ev.cscRechitClusterNRechitChamberMinus12,
            ev.cscRechitClusterNRechitChamberPlus11,
            ev.cscRechitClusterNRechitChamberPlus12,
        )
    ]

    pass_cluster_rechits_RE12 = [x == 0 for x in ev.cscRechitCluster_match_RE12_0p4]
    pass_cluster_rechits_MB1_RB1 = [
        (MB1 + RB1) == 0 for MB1, RB1 in zip(ev.cscRechitCluster_match_MB1Seg_0p4, ev.cscRechitCluster_match_RB1_0p4)
    ]
    pass_cluster_eta_cut = [abs(x) < 2.0 for x in ev.cscRechitClusterEta]

    #################
    ### Time Cuts ###
    #################

    pass_cluster_time_cut = [(-5.0 < x) and (x < 12.5) for x in ev.cscRechitClusterTime]
    pass_cluster_time_spread_cut = [x < 20 for x in ev.cscRechitClusterTimeSpread]
    pass_BKG_cluster_OOT_cut = [x < -12.5 for x in ev.cscRechitClusterTime]  # Background

    ############################
    ############################
    ############################
    ############################

    more_cols_dict = {
        "pass_csc_match": pass_csc_match,  #!
        "pass_met": pass_met,
        "pass_met50": pass_met50,
        "pass_met100": pass_met100,
        "pass_met150": pass_met150,
        "pass_cscRechitClusterSize_idxs": pass_cscRechitClusterSize_idxs,  #!
        "pass_cscRechitClusterSize_AtLeastOneAbove200": pass_cscRechitClusterSize_AtLeastOneAbove200,
        "pass_jet": pass_jet,
        "pass_lep": pass_lep,
        "pass_cscRechitCluster_CBID_idxs": pass_CBID_idxs,  #!
        "pass_cscRechitCluster_CC_dR_idxs": pass_CC_dR_idxs,  #!
        "pass_cscRechitCluster_CJ_dR_idxs": pass_CJ_dR_idxs,  #!
        "pass_cscRechitCluster_CM_dR_idxs": pass_CM_dR_idxs,  #!
        "pass_cscRechitCluster_rechits_ME11_ME12": pass_cluster_rechits_ME11_ME12,  #!
        "pass_cscRechitCluster_rechits_RE12": pass_cluster_rechits_RE12,  #!
        "pass_cscRechitCluster_rechits_MB1_RB1": pass_cluster_rechits_MB1_RB1,  #!
        "pass_cscRechitCluster_eta_cut": pass_cluster_eta_cut,  #!
        "pass_cscRechitCluster_time_cut": pass_cluster_time_cut,  #!
        "pass_cscRechitCluster_time_spread_cut": pass_cluster_time_spread_cut,  #!
        "pass_BKG_cscRechitCluster_OOT_cut": pass_BKG_cluster_OOT_cut,  #!
    }

    evs_dict["raw"] = build_ntuple_row(ev, tdict=evs_dict["raw"], more_cols_dict=more_cols_dict)
    evs_dict["rawM"] = build_ntuple_row(
        ev, cut_val=pass_csc_match, col="cscRechitCluster", tdict=evs_dict["rawM"], more_cols_dict=more_cols_dict
    )

    # Rechit, CBID, dR, Chamber, Time

    ###############################

    fill_all_histos("raw", evs_dict["raw"])

    cut_dict = {
        "raw": ([], ""),
        "rechit": ([], ""),
        "cbid": ([], ""),
        "dR": ([], ""),
        "chamber": ([], ""),
        "time": ([], ""),
    }

    # for cut_name, cut in

# print(evs_dict)
print(iev, len(evs_dict["raw"]["met"]))  # , len(evs_dict['rawM']['met']))

fout = rt.TFile(f'{".".join(fname.split(".")[:-1])}_histos.root', "RECREATE")
for name, histo in histos.items():
    fout.Add(histo)
fout.Write()
fout.Close()

ff.Close()
