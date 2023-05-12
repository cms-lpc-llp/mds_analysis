"""Class holding functions for processing muon system ntuples"""
import numpy as np

import ROOT as rt
import uproot as upr

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
        name = str(np.random.randint(0, 1e9))

    hh = rt.TH1F(name, title, *bins)
    hh.SetLineWidth(4)

    if isinstance(x, ak.Array):
        x = ak.flatten(x, -1)

    norm = 1 / len(x) if "norm" in kwargs and kwargs["norm"] == True else 1
    [hh.Fill(xx, norm) for xx in x]

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
        name = str(np.random.randint(0, 1e9))

    hh = rt.TH2F(name, title, *bins)

    if method == "s":  # SEQUENTIAL
        if isinstance(x, ak.Array):
            x = ak.flatten(x, -1)
        if isinstance(y, ak.Array):
            y = ak.flatten(y, -1)

        #rtnp.fill_hist(hh, x, y)
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
    hhs = [hh.Clone() for hh in hhs]
    if "ccs" in kwargs:
        ccs = kwargs["ccs"]

    lat, leg = get_lat_leg(kwargs["legxy"] if "legxy" in kwargs else (0.6, 0.7, 0.85, 0.95))
    if "norm" in kwargs and kwargs["norm"]:
        for hh in hhs:
            hh.Scale(1 / (hh.Integral() if hh.Integral() else 1))

    ymax = max([hh.GetMaximum() * (kwargs["ymax_mult"] if "ymax_mult" in kwargs else 1.05) for hh in hhs])
    for hh, tt, cc in zip(hhs, tts, ccs):
        hh.SetMaximum(ymax)
        if "ymin" in kwargs:
            hh.SetMinimum(kwargs["ymin"])
        hh.SetLineColor(cc)
        hh.Draw("same hist")

        leg.AddEntry(hh, tt, "L")

    leg.Draw()
    return hhs, lat, leg


############
def draw_dt_r_boxes(hh):
    ymin, ymax = hh.GetMinimum(), hh.GetMaximum()
    xmin, xmax = hh.GetXaxis().GetXmin(), hh.GetXaxis().GetXmax()
    boxes = []

    # if xmin < 181.1:
    #     boxes.append(rt.TBox(max(xmin, 181.1), ymin, 286.4, ymax))  # HCAL barrel
    # if xmin < 295:
    #     boxes.append(rt.TBox(max(xmin, 295), ymin, 377.5, ymax))  # solenoid
    #     #boxes.append(rt.TBox(max(xmin, 295), ymin, 380, ymax))  # solenoid
    # boxes.append(rt.TBox(max(xmin, 380), ymin, 402, ymax))  # b4 MB1
    boxes.append(rt.TBox(max(xmin, 180), ymin, 402, ymax))  # b4 MB1
    boxes.append(rt.TBox(449, ymin, 490, ymax))  # b4 MB2
    boxes.append(rt.TBox(533, ymin, 597, ymax))  # b4 MB3
    boxes.append(rt.TBox(636, ymin, 700, ymax))  # b4 MB4
    boxes.append(rt.TBox(738, ymin, xmax, ymax))  # beyond CMS

    for b in boxes:
        b.SetFillColor(15)
        b.SetFillStyle(3001)
        b.Draw("same")

    l = rt.TLatex()
    l.SetTextSize(0.08)
    l.SetTextAlign(12)
    l.SetTextColor(12)
    l.SetTextAngle(90)
    l.DrawLatex(230, ymax * 0.3, "Steel")
    #l.DrawLatex(230, ymax * 0.3, "HCAL")
    #l.DrawLatex(335, ymax * 0.3, "Solenoid")
    #l.DrawLatex(390, ymax * 0.3, "Steel")

    l2 = rt.TLatex()
    l2.SetTextSize(0.06)
    l2.SetTextColor(13)
    l2.SetTextAngle(90)
    l2.DrawLatex(780, ymax * 0.5, "Beyond CMS")
    text = rt.TLatex()
    text.SetTextSize(0.04)
    text.DrawLatex(400, ymax * 1.01, "MB1")
    text.DrawLatex(490, ymax * 1.01, "MB2")
    text.DrawLatex(600, ymax * 1.01, "MB3")
    text.DrawLatex(700, ymax * 1.01, "MB4")

    return boxes


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
            bins = (100, 180, 800)

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


##################################################
##################################################
##################################################


class MuonSystemRDF:
    """Handler for working with muon system ntuples using an RDataFrame, works with 2tag CSC-DT (CSC triggers) only"""

    def __init__(self, file_name, tree_name="MuonSystem", isMC=False, nev=None, rdf=None) -> None:
        self.isMC = isMC
        self.file_name = file_name
        self.tree_name = tree_name
        self.nev = nev

        if rdf is None:
            self.rdf = rt.RDataFrame(tree_name, file_name)
            if rt.IsImplicitMTEnabled():
                print("Disabling 'Range' for IMT -- will load ALL events!")
            else:
                if self.nev is not None:
                    self.rdf = self.rdf.Range(self.nev)
        else:
            self.rdf = rdf

        self.nev = self.rdf.Count().GetValue()

    def get(self, key):
        ctype = self.rdf.GetColumnType(key)
        print(ctype)
        if ctype == "int":
            return self.rdf.Take[int](key)
        elif ctype == "float":
            return self.rdf.Take[float](key)
        elif ctype == "bool":
            return self.rdf.Take[bool](key)
        elif ctype == "RVec<int>":
            return self.rdf.Take[rt.VecOps.RVec[int]](key)
        elif ctype == "RVec<float>":
            return self.rdf.Take[rt.VecOps.RVec[float]](key)
        elif ctype == "RVec<bool>":
            return self.rdf.Take[rt.VecOps.RVec[bool]](key)
        else:
            raise ValueError(f"Column type '{ctype}' not known!")
        #return self.rdf.AsNumpy(key)

    def __getitem__(self, key):
        return self.get(key)

    def __setitem__(self, key, value):
        self.Define(key, value)

    def __copy__(self):
        pass

    def Define(self, key, value, implicit=True):
        if key in self.rdf.GetColumnNames():
            print("In define - redefining -", key, value)
            rdf = self.rdf.Redefine(key, value)
            print("In define - finished redefining")
        else:
            print("In define - defining -", key, value)
            rdf = self.rdf.Define(key, value)
            print("In define - finished defining")

        if implicit:
            self.rdf = rdf
        else:
            self = MuonSystemRDF(self.file_name, self.tree_name, self.isMC, self.nev, rdf=rdf)

        return self

    def Filter(self, f, system="event", implicit=True):
        if system == "event":
            print("In filter - event -", f)
            rdf = self.rdf.Filter(f)
            print("In filter - finished event filter")
            if implicit:
                self.rdf = rdf
            else:
                self = MuonSystemRDF(self.file_name, self.tree_name, self.isMC, self.nev, rdf=rdf)

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

            system = pre[system]
            self = self.Define("cut", f, implicit=implicit)
            for k in self.rdf.GetColumnNames():
                k = str(k)
                if system == k[:len(system)]:
                    self = self.Define(k, f"{k}[cut]", implicit=implicit)

            self = self.fix_nbranch(implicit=implicit)

        return self

    def Histo1D(self, *args):
        return self.rdf.Histo1D(*args).GetPtr()

    def Histo2D(self, *args):
        #Add in behavior for RVec vs RVec, not sure how...
        return self.rdf.Histo2D(*args).GetPtr()

    def Count(self):
        return self.rdf.Count().GetValue()

    def fix_nbranch(self, implicit=True):
        self = self.Define("nCscRechitClusters", "Sum(cscRechitClusterSize>0)", implicit=implicit)
        self = self.Define("nDtRechitClusters", "Sum(dtRechitClusterSize>0)", implicit=implicit)
        self = self.Define("nJets", "Sum(jetPt>0)", implicit=implicit)
        self = self.Define("nLeptons", "Sum(lepPt>0)", implicit=implicit)
        return self

    def find_2tag_pair(self, implicit=True):
        return self
        pass

    def L1_plateau(self):
        return self
        pass

    def jet_cut(self, implicit=True):  # AN-21-124
        self = self.Filter(
            "!(cscRechitClusterJetVetoLooseId & (cscRechitClusterJetVetoPt > 30.) & (abs(cscRechitClusterEta) < 2.4))",
            system="csc",
            implicit=implicit)
        self = self.Filter(
            "!(dtRechitClusterJetVetoLooseId & (dtRechitClusterJetVetoPt > 50.) & (abs(dtRechitClusterEta) < 2.4))",
            system="dt",
            implicit=implicit)
        return self

    def muon_cut(self, implicit=True):  # AN-19-154
        self = self.Filter(
            "!( (cscRechitClusterMuonVetoLooseId & (cscRechitClusterMuonVetoPt > 30.) & (abs(cscRechitClusterEta) < 2.4)) | (cscRechitClusterNRechitChamberMinus11 + cscRechitClusterNRechitChamberMinus12 + cscRechitClusterNRechitChamberPlus11 + cscRechitClusterNRechitChamberPlus12 > 0) )",
            system="csc",
            implicit=implicit)
        self = self.Filter(
            "!( (dtRechitClusterMuonVetoLooseId & (dtRechitClusterMuonVetoPt > 10.) & (abs(dtRechitClusterEta) < 2.4)) | (dtRechitClusterNSegStation1 > 0) )",
            system="dt",
            implicit=implicit)
        return self

    def time_cut(self, time="it", system='cscdt', implicit=True):
        if time == "oot":
            if 'csc' in system:
                self = self.Filter("(cscRechitClusterTimeWeighted  < -12.5) | (cscRechitClusterTimeWeighted > 50)",
                                   system="csc",
                                   implicit=implicit)
            if 'dt' in system:
                self = self.Filter("dtRechitCluster_match_RPCBx_dPhi0p5 != 0", system="dt", implicit=implicit)
        elif time == "it":
            if 'csc' in system:
                self = self.Filter("abs(cscRechitClusterTimeWeighted) < 12.5", system="csc", implicit=implicit)
            if 'dt' in system:
                self = self.Filter("dtRechitCluster_match_RPCBx_dPhi0p5 == 0", system="dt", implicit=implicit)
        return self


##################################################
##################################################
##################################################


class MuonSystem:
    """Handler for working with Muon System ntuples using uproot"""
    _hlt_columns = {
        "HLT_CaloMET60_DTCluster50": 562,
        "HLT_CaloMET60_DTClusterNoMB1S50": 563,
        "HLT_L1MET_DTCluster50": 564,
        "HLT_L1MET_DTClusterNoMB1S50": 565,
        "HLT_CscCluster_Loose": 566,
        "HLT_CscCluster_Medium": 567,
        "HLT_CscCluster_Tight": 568,
        "HLT_L1CSCCluster_DTCluster50": 569,
        "HLT_L1CSCCluser_DTCluster75": 570,
    }

    def __init__(self, file_name, tree_name="MuonSystem", isMC=False, nev=None) -> None:
        self.isMC = isMC
        self.file_name = file_name
        self.tree_name = tree_name
        self.nev = nev
        self.fupr = upr.open(file_name + ":" + tree_name)

        # self.keys = {"met", "runNum"}
        # self.ms = self.fupr.arrays(("met", "runNum"), entry_stop=self.nev)
        self.ms = {"met": self.fupr["met"].array(entry_stop=self.nev)}
        self.cuts = []

        #self.ms = self.fupr.arrays(array_cache="inherit",entry_stop=self.nev)

    def __getitem__(self, key):
        return self.get(key)

    def __setitem__(self, key, value):
        if isinstance(key, str):
            self.ms[key] = value
        else:
            for k in self.ms:  # what was this supposed to be? set a slice?
                self.ms[k] = self.ms[k][value]
        return self

    def get(self, key):
        # if key not in self.ms.fields:
        #    dd = self.fupr.arrays(key, entry_stop=self.nev)[key]
        if key not in self.ms:
            if 'HLT' in key[:3]:
                if key == 'HLT' or key == 'HLTDecision':
                    raise ValueError("Memory overflow if you try to load all of HLT.")
                if key[4:].isdigit():
                    idx = int(key[4:])
                else:
                    idx = self._hlt_columns[key]

                print(idx)

                dd = self.fupr['HLTDecision'].array(entry_stop=self.nev)
                #dd = self['HLTDecision']
                dd = dd[:, idx]
                #dd = self.fupr['HLTDecision'].array(filter_branch=lambda b: b[:,idx], entry_stop=self.nev)
            else:
                dd = self.fupr[key].array(entry_stop=self.nev)

            for system, idxs in self.cuts:
                if system == "event":
                    dd = dd[idxs]
                else:
                    if system == key[:len(system)]:
                        dd = dd[idxs]

            self.ms[key] = dd

            self.fix_nbranch()
        return self.ms[key]

    def apply_cut(self, idxs, system: str = "event"):
        system = system.lower()

        if system == "event":
            for k in self.ms:
                self[k] = self[k][idxs]
            # self.ms = self.ms[idxs]
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

            system = pre[system]
            # for k in self.ms.fields:
            for k in self.ms:
                if system == k[:len(system)]:
                    self.ms[k] = self.ms[k][idxs]

        self.cuts.append((system, idxs))
        self.fix_nbranch()

    def fix_nbranch(self):
        self["nCscRechitClusters"] = asum(self.get("cscRechitClusterSize") > 0)
        self["nDtRechitClusters"] = asum(self.get("dtRechitClusterSize") > 0)
        self["nGLLP"] = asum(self.get("gLLP_pt") > 0)
        self["nLeptons"] = asum(self.get("lepPt") > 0)
        self["nJets"] = asum(self.get("jetPt") > 0)

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
            cuts.append(self.get(det + "RechitCluster_match_gLLP"))

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
                        aabs(self.get(pree + "eta")) < max_csc_eta,
                        aabs(self.get(prev + "r")) < max_csc_r,
                        aabs(self.get(prev + "z")) > min_csc_z,
                        aabs(self.get(prev + "z")) < max_csc_z,
                    ))
            if "dt" == det:
                cuts.append(
                    land(
                        aabs(self.get(prev + "r")) > min_dt_r,
                        aabs(self.get(prev + "r")) < max_dt_r,
                        aabs(self.get(prev + "z")) < max_dt_z,
                    ))

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
        return asum(self.get("HLTDecision")[:, 566:571]) > 0

    def ndet_cut(self, ncsc=1, ndt=1, op="&", implicit=False):
        """if ncsc or ndt is an iterable -> [inclusive, exclusive)
        if you're not cutting on one of the vars LEAVE OP as default"""
        lcsc, ldt = "nCscRechitClusters", "nDtRechitClusters"
        if isinstance(ncsc, (tuple, list)):
            idx_csc = land(self.get(lcsc) <= ncsc[0], self.get(lcsc) <= ncsc[1])
        elif isinstance(ncsc, int):
            idx_csc = self.get(lcsc) == ncsc
        else:
            idx_csc = np.ones_like(self.get("met"), dtype=bool)

        if isinstance(ndt, (tuple, list)):
            idx_dt = land(self.get(ldt) <= ndt[0], self.get(ldt) <= ndt[1])
        elif isinstance(ndt, int):
            idx_dt = self.get(ldt) == ndt
        else:
            idx_dt = np.ones_like(self.get("met"), dtype=bool)

        if op == "&":
            return land(idx_csc, idx_dt)
        elif op == "|":
            return lor(idx_csc, idx_dt)
        elif op == "^":
            return lxor(idx_csc, idx_dt)

    def met_cut(self, met_min=50, met_max=None, implicit=False):
        if met_min is None:
            return self.get("met") < met_max
        if met_max is None:
            return met_min <= self.get("met")

        return land(met_min <= self.get("met"), self.get("met") < met_max)

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
                                self.get("cscRechitClusterMuonVetoLooseId"),
                                self.get("cscRechitClusterMuonVetoPt") > 30,
                                aabs(self.get("cscRechitClusterEta")) < 2.4,
                            ),
                            self.get("cscRechitClusterNRechitChamberMinus11") +
                            self.get("cscRechitClusterNRechitChamberMinus12") +
                            self.get("cscRechitClusterNRechitChamberPlus11") +
                            self.get("cscRechitClusterNRechitChamberPlus12") > 0,
                        )))  # AN-19-154
            if "dt" == det:
                cuts.append(
                    lnot(
                        lor(
                            land(
                                self.get("dtRechitClusterMuonVetoLooseId"),
                                self.get("dtRechitClusterMuonVetoPt") > 10,
                                aabs(self.get("dtRechitClusterEta")) < 2.4,
                            ),
                            self.get("dtRechitClusterNSegStation1") > 0,
                        )))  # AN-19-154
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
                            self.get("cscRechitClusterJetVetoLooseId"),
                            self.get("cscRechitClusterJetVetoPt") > 30,
                            aabs(self.get("cscRechitClusterEta")) < 2.4,
                        )))  # AN-21-124
            if det == "dt":
                cuts.append(
                    lnot(
                        land(
                            self.get("dtRechitClusterJetVetoLooseId"),
                            self.get("dtRechitClusterJetVetoPt") > 50,
                            aabs(self.get("dtRechitClusterEta")) < 2.4,
                        )))  # AN-21-124

        return cuts[0] if len(cuts) == 1 else cuts

    def time_cut(self, dets=["csc", "dt"], implicit=False):
        if not isinstance(dets, (tuple, list)):
            dets = [dets]

        cuts = []
        for det in dets:
            if "csc" == det:
                cuts.append(
                    land(
                        -5 < self.get("cscRechitClusterTimeWeighted"),
                        self.get("cscRechitClusterTimeWeighted") < 12.5,
                    ))
            if "dt" == det:
                cuts.append(self.get("dtRechitCluster_match_RPCBx_dPhi0p5") == 0)

        return cuts[0] if len(cuts) == 1 else cuts
