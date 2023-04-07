"""Class holding functions for processing muon system ntuples"""
import numpy as np

import ROOT as rt
import uproot as upr
from root_numpy import fill_hist

from itertools import product as combs
from src.helper_functions import lnot, land, lor, lxor, asum, aabs

############
from src.histo_utilities import std_color_list
import awkward as ak


def get_lat_leg(leg_coords=(0.6, 0.7, 0.95, 0.8)):
    lat = rt.TLatex()
    lat.SetTextColor(rt.kRed)
    lat.SetTextSize(0.03)
    lat.SetTextAlign(33)

    leg = rt.TLegend(*leg_coords)
    leg.SetTextSize(0.025)
    leg.SetBorderSize(0)
    leg.SetEntrySeparation(0.01)

    return lat, leg


def H1D(x, title, bins, **kwargs):
    name = title.split(";")[0]
    if name == "":
        name = str(np.random.randint(0, 10000))

    hh = rt.TH1F(name, title, *bins)
    hh.SetLineWidth(4)

    if isinstance(x, ak.Array):
        x = ak.flatten(x, -1)

    #weights = None
    #if "norm" in kwargs and kwargs['norm']:
    #    weights = np.ones_like(x) / len(x)
    #fill_hist(hh, x, weights)
    for xx in x:
        hh.Fill(xx, 1 / len(x) if "norm" in kwargs and kwargs["norm"] == True else 1)

    if "lw" in kwargs:
        hh.SetLineWidth(kwargs["lw"])
    if "c" in kwargs:
        hh.SetLineColor(kwargs["c"])
    if "ymin" in kwargs:
        hh.SetMinimum(kwargs["ymin"])
    if "ymax" in kwargs:
        hh.SetMaximum(kwargs["ymax"])

    return hh


def H2D(x, y, title, bins, method="c", **kwargs):
    name = title.split(";")[0]
    if name == "":
        name = str(np.random.randint(0, 10000))

    hh = rt.TH2F(name, title, *bins)

    if method == "s":  # SEQUENTIAL
        if isinstance(x, ak.Array):
            x = ak.flatten(x, -1)
        if isinstance(y, ak.Array):
            y = ak.flatten(y, -1)

        for xx, yy in zip(x, y):
            hh.Fill(xx, yy)  # , 1 / len(x) if "norm" in kwargs and kwargs["norm"] == True else 1)
    elif method == "c":  # COMBINATORIAL
        if len(x) != len(y):
            raise ValueError(f"x and y don't have the same length! ({len(x):=}, {len(y):=}")

        for xr, yr in zip(x, y):
            # print(xr, yr)
            if isinstance(xr, ak.Array):
                xr = ak.to_list(xr)
            else:
                xr = [xr]
            if isinstance(yr, ak.Array):
                yr = ak.to_list(yr)
            else:
                yr = [yr]

            if len(xr) < 1 or len(yr) < 1:
                continue
            # if isinstance(xr, ak.Array):
            #     if len(xr) == 0:
            #         continue
            #     xr = ak.flatten(xr, -1)
            #     print("here1")
            # else:
            #     xr = np.array(xr)
            #     print("here2")

            # if isinstance(yr, ak.Array):
            #     if len(yr) == 0:
            #         continue
            #     yr = ak.flatten(yr, -1)
            # else:
            #     yr = np.array(yr)

            # print('\t', xr, yr)
            for xx, yy in combs(xr, yr):
                hh.Fill(xx, yy)  # , 1 / len(x) if "norm" in kwargs and kwargs["norm"] == True else 1)

    return hh


def multi_plot(hhs, tts, **kwargs):
    ccs = std_color_list
    if "ccs" in kwargs:
        ccs = kwargs["ccs"]

    lat, leg = get_lat_leg(kwargs["legxy"] if "legxy" in kwargs else (0.6, 0.7, 0.85, 0.95))
    if 'norm' in kwargs and kwargs['norm']:
        for hh in hhs:
            hh.Scale(1/(hh.Integral() if hh.Integral() else 1))
    ymax = max([hh.GetMaximum() * (kwargs["ymax_mult"] if "ymax_mult" in kwargs else 1.05) for hh in hhs])
    for hh, tt, cc in zip(hhs, tts, ccs):
        hh.SetMaximum(ymax)
        if "ymin" in kwargs:
            hh.SetMinimum(kwargs["ymin"])
        hh.SetLineColor(cc)
        hh.Draw("same hist")

        leg.AddEntry(hh, tt, "L")

    leg.Draw()
    return lat, leg


############
def draw_csc_z_boxes(hh):
    ymin, ymax = hh.GetMinimum(), hh.GetMaximum()
    xmin, xmax = hh.GetXaxis().GetXmin(), hh.GetXaxis().GetXmax()
    boxes = []

    boxes.append(rt.TBox(xmin, ymin, 568, ymax))  # in front of ME11
    boxes.append(rt.TBox(632, ymin, 671, ymax))  # between ME11 and ME12
    boxes.append(rt.TBox(724, ymin, 789, ymax))  # between ME12 and station2
    boxes.append(rt.TBox(849, ymin, 911, ymax))  # between station2 and station3
    boxes.append(rt.TBox(970, ymin, 1002, ymax))  # between station3 and station4
    boxes.append(rt.TBox(1073, ymin, xmax, ymax))  # beyond CMS
    for b in boxes:
        b.SetFillColor(15)
        b.SetFillStyle(3001)
        b.Draw("same")

    l = rt.TLatex()
    l.SetTextSize(0.08)
    l.SetTextColor(12)
    l.SetTextAngle(90)
    l.DrawLatex(xmin + 80, ymax * 0.4, "Steel")

    l2 = rt.TLatex()
    l2.SetTextSize(0.06)
    l2.SetTextColor(13)
    l2.SetTextAngle(90)
    l2.DrawLatex(1110, ymax * 0.5, "Beyond CMS")
    text = rt.TLatex()
    text.SetTextSize(0.04)
    text.DrawLatex(570, ymax * 1.01, "ME1/1")
    text.DrawLatex(660, ymax * 1.01, "ME1/2-3")
    text.DrawLatex(795, ymax * 1.01, "ME2")
    text.DrawLatex(920, ymax * 1.01, "ME3")
    text.DrawLatex(1015, ymax * 1.01, "ME4")

    return boxes


def make_cluster_eff_1D(ms, det, xl="z", cuts=False):
    lcltr, lgllp = det + "RechitCluster_match_gLLP_decay_", "gLLP_decay_vertex_"
    title = f";gLLP {xl.upper()} Decay Vertex [cm];{det.upper()} efficiency"

    if xl == "z":
        if det == "csc":
            bins = (100, 400, 1100)
        elif det == "dt":
            bins = (100, 0, 700)
    elif xl == "r":
        if det == "csc":
            bins = (100, 0, 800)
        elif det == "dt":
            bins = (100, 0, 800)

    sel_gllp_in_det = ms.in_det_cut(det, "gllp")
    sel_cltr_in_det = ms.in_det_cut(det, "rechit")

    sel_match = ms.match_cut(det)
    sel_num = land(asum(sel_match) > 0, asum(sel_cltr_in_det) == asum(sel_match))
    if cuts:
        sel_jets = ms.jet_veto_cut(det)
        sel_muon = ms.muon_veto_cut(det)
        sel_time = ms.time_cut(det)
        sel_match = land(sel_match, sel_jets, sel_muon, sel_time)

    hnum = H1D(aabs(ms[lcltr + xl])[sel_match][sel_num], title, bins=bins)
    hden = H1D(aabs(ms[lgllp + xl])[sel_gllp_in_det], title, bins=bins)

    hnum.Divide(hden)

    return hnum


############


class MuonSystem:
    """Handler for working with Muon System ntuples"""

    def __init__(self, file_name, tree_name="MuonSystem", isMC=False) -> None:
        self.isMC = isMC
        self.file_name = file_name
        self.tree_name = tree_name
        self.fuproot = upr.open(file_name + ":" + tree_name)
        self.ms = self.fuproot.arrays(entry_stop=5_000_000)#5_000_000)

    def __getitem__(self, key):
        return self.ms[key]

    def __setitem__(self, key, value):
        self.ms[key] = value

    def get(self, name):
        return self.ms[name]

    def apply_cut(self, idxs, system: str = "event"):
        system = system.lower()

        if system == "event":
            self.ms = self.ms[idxs]
        else:
            pre = {
                "csc": "cscRechitCluster",
                "dt": "dtRechitCluster",
                "gllp": "gLLP",
                "lep": "lep",
                "jet": "jet",
            }

            if system not in pre:
                raise ValueError(f"Invaid system {system}.")
            
            pre = pre[system]
            for k in self.ms.fields:
                if pre == k[: len(pre)]:
                    self.ms[k] = self.ms[k][idxs]
        self.fix_nbranch()

    def fix_nbranch(self):
        self.ms["nCscRechitClusters"] = asum(self["cscRechitClusterSize"] > 0)
        self.ms["nDtRechitClusters"] = asum(self["dtRechitClusterSize"] > 0)
        self.ms["nGLLP"] = asum(self["gLLP_pt"] > 0)
        self.ms["nLeptons"] = asum(self["lepPt"] > 0)
        self.ms["nJets"] = asum(self["jetPt"] > 0)

    #################################################
    ## CUTS, returns idxs, implicit OFF by default ##
    #################################################

    def match_cut(self, dets=["csc", "dt"], implicit=False):
        if not self.isMC:
            raise ValueError("Not Monte Carlo data, cannot match clusters!")

        if not isinstance(dets, (tuple, list)):
            dets = [dets]
        cuts = []

        for det in dets:
            cuts.append(self[det + "RechitCluster_match_gLLP"])

        return cuts[0] if len(cuts) == 1 else cuts

    def in_det_cut(self, dets=["csc", "dt"], system="rechit", implicit=False):
        max_csc_eta, max_csc_r, min_csc_z, max_csc_z = 3, 800, 400, 1200
        min_dt_r, max_dt_r, max_dt_z = 200, 800, 700
        if not isinstance(dets, (tuple, list)):
            dets = [dets]
        cuts = []

        for det in dets:
            if system == "rechit":
                pree = det + "RechitCluster_match_gLLP_"
                prev = pree + "decay_"
            elif system == "gllp":
                pree = "gLLP_"
                prev = pree + "decay_vertex_"
            if "csc" == det:
                cuts.append(
                    land(
                        aabs(self[pree + "eta"]) < max_csc_eta,
                        aabs(self[prev + "r"]) < max_csc_r,
                        aabs(self[prev + "z"]) > min_csc_z,
                        aabs(self[prev + "z"]) < max_csc_z,
                    )
                )
            if "dt" == det:
                cuts.append(
                    land(
                        aabs(self[prev + "r"]) > min_dt_r,
                        aabs(self[prev + "r"]) < max_dt_r,
                        aabs(self[prev + "z"]) < max_dt_z,
                    )
                )

        return cuts[0] if len(cuts) == 1 else cuts

    def hlt_cut(self, implicit=False):
        """
        562  HLT_CaloMET60_DTCluster50
        563  HLT_CaloMET60_DTClusterNoMB1S50
        564  HLT_L1MET_DTCluster50
        565  HLT_L1MET_DTClusterNoMB1S50
        566  HLT_CscCluster_Loose
        567  HLT_CscCluster_Medium
        568  HLT_CscCluster_Tight
        569  HLT_L1CSCShower_DTCluster50
        570  HLT_L1CSCShower_DTCluster75"""
        return asum(self["HLTDecision"][:, 566:571]) > 0

    def ndet_cut(self, ncsc=1, ndt=1, op="&", implicit=False):
        """if ncsc or ndt is an iterable -> [inclusive, exclusive)
        if you're not cutting on one of the vars LEAVE OP as default"""
        lcsc, ldt = "nCscRechitClusters", "nDtRechitClusters"
        if isinstance(ncsc, (tuple, list)):
            idx_csc = land(self[lcsc] <= ncsc[0], self[lcsc] <= ncsc[1])
        elif isinstance(ncsc, int):
            idx_csc = self[lcsc] == ncsc
        else:
            idx_csc = np.ones_like(self["met"], dtype=bool)

        if isinstance(ndt, (tuple, list)):
            idx_dt = land(self[ldt] <= ndt[0], self[ldt] <= ndt[1])
        elif isinstance(ndt, int):
            idx_dt = self[ldt] == ndt
        else:
            idx_dt = np.ones_like(self["met"], dtype=bool)

        if op == "&":
            return land(idx_csc, idx_dt)
        elif op == "|":
            return lor(idx_csc, idx_dt)
        elif op == "^":
            return lxor(idx_csc, idx_dt)

    def met_cut(self, met_min=50, met_max=None, implicit=False):
        if met_min is None:
            return self["met"] < met_max
        if met_max is None:
            return met_min <= self["met"]

        return land(met_min <= self["met"], self["met"] < met_max)

    def muon_veto_cut(self, dets=["csc", "dt"], implicit=False):
        if not isinstance(dets, (tuple, list)):
            dets = [dets]
        cuts = []
        for det in dets:
            if "csc" == det:
                cuts.append(
                    lnot(
                        lor(
                            land(
                                self["cscRechitClusterMuonVetoLooseId"],
                                self["cscRechitClusterMuonVetoPt"] > 30,
                                aabs(self["cscRechitClusterEta"]) < 2.4,
                            ),
                            self["cscRechitClusterNRechitChamberMinus11"]
                            + self["cscRechitClusterNRechitChamberMinus12"]
                            + self["cscRechitClusterNRechitChamberPlus11"]
                            + self["cscRechitClusterNRechitChamberPlus12"]
                            > 0,
                        )
                    )
                )  # AN-19-154
            if "dt" == det:
                cuts.append(
                    lnot(
                        lor(
                            land(
                                self["dtRechitClusterMuonVetoLooseId"],
                                self["dtRechitClusterMuonVetoPt"] > 10,
                                aabs(self["dtRechitClusterEta"]) < 2.4,
                            ),
                            self["dtRechitClusterNSegStation1"] > 0,
                        )
                    )
                )  # AN-19-154
        return cuts[0] if len(cuts) == 1 else cuts

    def jet_veto_cut(self, dets=["csc", "dt"], implicit=False):
        if not isinstance(dets, (tuple, list)):
            dets = [dets]

        cuts = []
        for det in dets:
            if det == "csc":
                cuts.append(
                    lnot(
                        land(
                            self["cscRechitClusterJetVetoLooseId"],
                            self["cscRechitClusterJetVetoPt"] > 30,
                            aabs(self["cscRechitClusterEta"]) < 2.4,
                        )
                    )
                )  # AN-21-124
            if det == "dt":
                cuts.append(
                    lnot(
                        land(
                            self["dtRechitClusterJetVetoLooseId"],
                            self["dtRechitClusterJetVetoPt"] > 50,
                            aabs(self["dtRechitClusterEta"]) < 2.4,
                        )
                    )
                )  # AN-21-124

        return cuts[0] if len(cuts) == 1 else cuts

    def time_cut(self, dets=["csc", "dt"], implicit=False):
        if not isinstance(dets, (tuple, list)):
            dets = [dets]

        cuts = []
        for det in dets:
            if "csc" == det:
                cuts.append(
                    land(
                        -5 < self["cscRechitClusterTimeWeighted"],
                        self["cscRechitClusterTimeWeighted"] < 12.5,
                    )
                )
            if "dt" == det:
                cuts.append(self["dtRechitCluster_match_RPCBx_dPhi0p5"] == 0)

        return cuts[0] if len(cuts) == 1 else cuts
