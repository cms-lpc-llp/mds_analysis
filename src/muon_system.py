"""Class holding functions for processing muon system ntuples

Tasks:
    - MuonSystemAwkward
        - continue developing the basic structure
        - L1 cuts
        - Test HLT cut 
        - fix match
        - 2tag
    
    - MuonSystemRDF
        - I messed up 'self_out = self'
        - Won't propagate cuts in some functions
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

import numpy as np
import numba as nb

import ROOT as rt
from ROOT import TLatex, TLegend, TBox, TCanvas, TH1D, TH2D
# from ROOT import kRed, kBlue, kGreen, kCyan, kMagenta, kYellow, kBlack, kAzure

# from ROOT.VecOps import RVec

import uproot as upr
import awkward as ak
import awkward.numba
from src.helper_functions import alert, lnot, land, lor, lxor, asum, aabs

# from itertools import product as combs
from src.histo_utilities import std_color_list
SCL = std_color_list


################################################################
################################################################
################################################################


class MuonSystemAwkward:
    """Handler for working with muon system ntuples using an Awkward arrays"""

    def __init__(
        self,
        file_name: str,
        name: str,
        tree_name: str="MuonSystem",
        nev: int=None,
        is_mc: bool=False,
        lumi: float = None,
        cache=False,
        implicit: bool = True,
        verbose: bool = True
    ) -> None:
        """Initialize a new instance of MuonSystemAwkward"""

        if verbose:
            print(f"Building MuonSystemAwkward '{name}' -")
            print(f"  is_mc  = {is_mc}")
            print(f"  events = {nev}")
            print(f"  tree   = '{tree_name}'")
            print(f"  file   = '{file_name}'")

        ##########

        self.file_name = file_name
        self.name = name
        self.tree_name = tree_name
        self.is_mc = is_mc
        self.lumi = lumi
        self.nev = nev
        self.cache = cache
        self.implicit = implicit
        self.verbose = verbose

        ##########

        self.cut = True
        self.efficiency_denom = None
        self.efficiency = None
        self.colors = None

        ##########

        # From uproot documentation:
        # "entry_stop (None or int) –
        # The first entry to exclude (i.e. one greater than the last entry to include).
        # If None, stop at num_entries. If negative, count from the end, like a Python slice."
        if self.nev < 0:
            self.nev = None

        if not self.implicit:
            raise NotImplementedError("Non-implicit cuts do not work yet.")

        ##########

        self.ms_read = {}
        self.ms = upr.open(path=self.file_name + ":" + self.tree_name)  # , array_cache='100 kB')
        self.ms_read["sel_evt"] = ak.ones_like(self.ms["met"].array(entry_stop=self.nev), dtype=bool)
        self.ms_read["sel_csc"] = ak.ones_like(self.ms["cscRechitClusterSize"].array(entry_stop=self.nev), dtype=bool)
        self.ms_read["sel_dt"] = ak.ones_like(self.ms["dtRechitClusterSize"].array(entry_stop=self.nev), dtype=bool)

        if len(self.ms_read["sel_evt"]) != self.nev:
            self.nev = len(self.ms_read["sel_evt"])
        print(f"  Extracted {self.nev:,} events")

        ##########

        # self.aliases = {
        #         "cscRechitClusterMe11Ratio" : "(cscRechitClusterNRechitChamberPlus11 + cscRechitClusterNRechitChamberMinus11) / cscRechitClusterSize",
        #         "cscRechitClusterMe12Ratio" : "(cscRechitClusterNRechitChamberPlus12 + cscRechitClusterNRechitChamberMinus12) / cscRechitClusterSize",
        #         "dtRechitClusterMb1Ratio" : "dtRechitClusterNHitStation1 / dtRechitClusterSize",
        #         "cscRechitClusterR" : "sqrt(cscRechitClusterX**2 + cscRechitClusterY**2)",
        #         "dtRechitClusterR" : "sqrt(dtRechitClusterX**2 + dtRechitClusterY**2)",
        #         "tag_dPhi": ,
        #         "tag_dEta": ,
        #         "tag_dR": ,
        #     }
        
        # hlt_info = {
        #     'HLT_CaloMET60_DTCluster50': 562,
        #     'HLT_CaloMET60_DTClusterNoMB1S50': 563,
        #     'HLT_L1MET_DTCluster50': 564,
        #     'HLT_L1MET_DTClusterNoMB1S50': 565,
        #     'HLT_CscCluster_Loose': 566,
        #     'HLT_CscCluster_Medium': 567,
        #     'HLT_CscCluster_Tight': 568,
        #     "HLT_L1CSCCluster_DTCluster50": 569,
        #     'HLT_L1CSCCluser_DTCluster75': 570,
        # }

        # self.aliases = self.aliases | {k: f"HLTDecision[:{self.nev},c]" for k,c in hlt_info.items()}

        ##########

        # ! I cannot figure out how to load HLTDecision without overflowing memory
        dev = 100_000
        hlt = np.array([], dtype=bool)
        for i in range(0, nev//dev + 1):
            # hlt = np.r_[hlt, ms["HLTDecision"].array(entry_start=i*dev, entry_stop=(i+1)*dev)[:,569]]
            hlt = np.r_[hlt, np.take(self.ms["HLTDecision"].array(entry_start=i*dev, entry_stop=(i+1)*dev), 569, 1)]
        self.ms_read["HLT_L1CSCCluster_DTCluster50"] = hlt
        # self.ms_read["HLT_L1CSCCluster_DTCluster50"] = self.ms["HLTDecision"].array(entry_stop=self.nev)[:,569]
        # self.ms_read["HLT_L1CSCCluster_DTCluster50"] = self.ms.arrays(
        #     "HLT_L1CSCCluster_DTCluster50", aliases=self.hlt_aliases, entry_stop=self.nev
        # )["HLT_L1CSCCluster_DTCluster50"]

        # self.ms_read["nCscRechitClusters"] = np.sum(self.ms_read["sel_csc"], axis=1)
        # self.ms_read["nDtRechitClusters"] = np.sum(self.ms_read["sel_dt"], axis=1)

    def __getitem__(self, key):
        # hlt_aa = self.ms.arrays(list(self.hlt_aliases.keys()), , entry_stop=self.nev)
        # for k, v in self.hlt_info.items():
        #     self.ms_read[k] = hlt_aa[k]

        mods = [""]
        if isinstance(key, str):
            key = self._fix_key(key)

            if "." in key:
                key, mods = key.split(".")[0], key.split(".")[1:]

            if key in self.ms_read:
                arr = self.ms_read[key]
            elif "nCsc" in key:
                arr = np.sum(self.ms_read["sel_csc"], axis=1)
            elif "nDt" in key:
                arr = np.sum(self.ms_read["sel_dt"], axis=1)
            elif key in self.ms:# or key in self.aliases:
                arr = self.ms[key].array(entry_stop=self.nev)
                # arr = self.ms.arrays(key, aliases=self.aliases, entry_stop=self.nev)[key]
            # else:
            #     raise ValueError(f"Invalid key '{key}'.")
            # elif key in self.aliases:
            #     arr
            else:
                # self.cut, pcut = False, self.cut
                if key == "cscRechitClusterMe11Ratio":
                    arr = (
                        self["cscRechitClusterNRechitChamberPlus11"] + self["cscRechitClusterNRechitChamberMinus11"]
                    ) / self["cscRechitClusterSize"]
                elif key == "cscRechitClusterMe12Ratio":
                    arr = (
                        self["cscRechitClusterNRechitChamberPlus12"] + self["cscRechitClusterNRechitChamberMinus12"]
                    ) / self["cscRechitClusterSize"]
                elif key == "dtRechitClusterMb1Ratio":
                    arr = self["dtRechitClusterNHitStation1"] / self["dtRechitClusterSize"]
                elif key == "cscRechitClusterR":
                    arr = np.sqrt(self["cscRechitClusterX"] ** 2 + self["cscRechitClusterY"] ** 2)
                elif key == "dtRechitClusterR":
                    arr = np.sqrt(self["dtRechitClusterX"] ** 2 + self["dtRechitClusterY"] ** 2)
                else:
                    raise ValueError(f"Invalid key '{key}'.")

                # self.cut = pcut

            if self.cut and key not in ("sel_evt","sel_csc","sel_dt") and len(arr) == self.nev:
                if "csc" in key[:3]:
                    arr = arr[self.ms_read["sel_csc"]]
                elif "dt" in key[:2]:
                    arr = arr[self.ms_read["sel_dt"]]
                arr = arr[self.ms_read["sel_evt"]]

            if self.lumi is not None and key == "weight":
                if self.is_mc:
                    arr = arr * self.lumi

            if self.cache and key not in self.ms_read:
                self.ms_read[key] = arr

            #TODO: for loop through mods so it's in order and chainable (move to front?)
            if "abs" in mods: 
                arr = np.abs(arr)
            if "log" in mods:
                arr = np.log(arr)
            if "log10" in mods:
                arr = np.log10(arr)
            if "sqrt" in mods:
                arr = np.log10(arr)

            return arr

        else:
            raise NotImplementedError("BROKEN implicit setting")
        #     imp = self.implicit
        #     self.set_implicit(False)
        #     sout = self.filter(sel=key, system="evt")
        #     self.set_implicit(imp)
        #     return sout

    def __setitem__(self, key, value):
        # print(key, value)
        key = self._fix_key(key)
        self.ms_read[key] = value

    def __copy__(self):
        pass

    def _fix_key(self, key):
        if "csc" in key[:3] and "RechitCluster" not in key:
            key = "cscRechitCluster" + key[3:]
        if "dt" in key[:2] and "RechitCluster" not in key:
            key = "dtRechitCluster" + key[2:]

        if "nCsc" == key:
            key = "nCscRechitClusters"
        if "nDt" == key:
            key = "nDtRechitClusters"
        if "nLep" == key:
            key = "nLeptons"
        if "nJet" == key:
            key = "nJets"

        if "dPhi" == key:
            key = "tag_dPhi"
        if "dEta" == key:
            key = "tag_dEta"
        if "dR" == key:
            key = "tag_dR"

        return key

    # def _fix_n_column(self, system):
    #     if system == "evt":
    #         return

    #     n_columns = {
    #         "csc": ("nCscRechitClusters", "sel_csc"),
    #         "dt": ("nDtRechitClusters", "sel_dt"),
    #         "gllp": ("nGLLP", "sel_gllp"),
    #         "jet": ("nJets", "sel_jet"),
    #         "lep": ("nLeptons", "sel_lep"),
    #     }
    #     n_key, test_col = n_columns[system]
    #     self[n_key] = np.sum(self[test_col], axis=1)

    def count(self, use_weight=True):
        if use_weight:
            if self.cut:
                return np.sum(self["weight"])
            else:
                return np.sum(self["weight"][self["sel_evt"]])
        return np.sum(self["sel_evt"])

    # def set_implicit(self, implicit: bool = True) -> bool:
    #     self.implicit = implicit
    #     return self.implicit

    def filter(self, sel, system="evt", invert=False, has_clusters=True):
        if invert:
            sel = ~ sel

        if self.efficiency_denom is not None:
            ed = self.efficiency_denom
            self.efficiency = np.sum(sel & self[f"{ed}_{system}"]) / np.sum(self[f"{ed}_{system}"])

        self.ms_read[f"sel_{system}"] = self.ms_read[f"sel_{system}"] & sel

        if has_clusters and system != "evt":
            # self.filter(self['nCscRechitClusters'] + self['nDtRechitClusters'] > 0, system='evt')
            self.ms_read["sel_evt"] = self.ms_read["sel_evt"] & (np.sum(self.ms_read[f"sel_{system}"], axis=1) > 0)

        return self

    def f(self, *args, **kwargs):
        return self.filter(*args, **kwargs)
    
    # def hist(name, title)

    def match_mc(self, system="csc,dt", check_decay=True, **kwargs):
        self.cut, pcut = False, self.cut

        max_csc_eta, max_csc_r, min_csc_z, max_csc_z = 3, 800, 400, 1200
        min_dt_r, max_dt_r, max_dt_z = 200, 800, 700
        if "csc" in system:
            self.filter(self["cscRechitCluster_match_gLLP"], system="csc")
            if check_decay:
                sel_decay = (
                    (self["cscRechitCluster_match_gLLP_eta"] < max_csc_eta)
                    & (self["cscRechitCluster_match_gLLP_decay_r"] < max_csc_r)
                    & (min_csc_z < np.abs(self["cscRechitCluster_match_gLLP_decay_z"]))
                    & (np.abs(self["cscRechitCluster_match_gLLP_decay_z"]) < max_csc_z)
                )
                self.filter(sel_decay, system="csc", **kwargs)

        if "dt" in system:
            self.filter(self["dtRechitCluster_match_gLLP"], system="dt")
            if check_decay:
                sel_decay = (
                    (min_dt_r < self["dtRechitCluster_match_gLLP_decay_r"])
                    & (self["dtRechitCluster_match_gLLP_decay_r"] < max_dt_r)
                    & (np.abs(self["dtRechitCluster_match_gLLP_decay_z"]) < max_dt_z)
                )
                self.filter(sel_decay, system="dt", **kwargs)

        self.cut = pcut
        return self

    def blind(self, region, **kwargs):
        self.cut, pcut = False, self.cut
        if region == "dtSize":
            pass
        elif region == "cscSize":
            pass
        elif region == "dtOOT":
            pass
        else:
            raise ValueError(f"No blind cut with region='{region}'.")
        self.cut = pcut
        return self

    def tag(self, tags="cscdt"): # only works with cscdt rn
        self.cut, pcut = False, self.cut

        @nb.njit
        def _delta_kin(csc_kin, dt_kin, sel_csc, sel_dt, tag, kin):
            dkin = np.zeros(len(tag))
            for i in range(len(tag)):
                if tag[i] == 0:
                    continue

                k0, k1 = 0, 0
                i0, i1 = -1, -1
                # eh, (np.argwhere(sel_csc[i]) not compatible with compiled ArrayView(?))
                for j, v in enumerate(sel_csc[i]):
                    if v:
                        if i0 == -1:
                            i0 = j
                            k0 = csc_kin[i][j]
                        else:
                            raise NotImplementedError("More than one CSC clusters found")
                for j, v in enumerate(sel_dt[i]):
                    if v:
                        if i1 == -1:
                            i1 = j
                            k1 = dt_kin[i][j]
                        else:
                            raise NotImplementedError("More than one DT clusters found")
                if i0 == -1 or i1 == -1:
                    raise NotImplementedError("Less than two clusters found.")

                if kin == "dphi":
                    dkin[i] = min(2 * np.pi - np.abs(k0 - k1), np.abs(k0 - k1))
                elif kin == "deta":
                    dkin[i] = np.abs(k0 - k1)
            return dkin
        
        @nb.njit
        def _delta_kins(csc_phi, csc_eta, dt_phi, dt_eta, met_phi, sel_csc, sel_dt, tag):
            dphi, deta, dr = (
                np.zeros(len(tag)),
                np.zeros(len(tag)),
                np.zeros(len(tag)),
            )
            for i in range(len(csc_eta)):
                if tag[i] == 0:
                    continue

                p0, p1, e0, e1 = 0, met_phi[i], 0, 0
                c0, c1, d0, d1 = -1, -1, -1, -1
                # eh, (np.argwhere(sel_csc[i]) not compatible with compiled ArrayView(?))
                for j, v in enumerate(sel_csc[i]):
                    if v:
                        if c0 == -1:
                            c0 = j
                        elif c1 == -1:
                            c1 = j
                        else:
                            raise NotImplementedError("More than two CSC clusters found")
                for j, v in enumerate(sel_dt[i]):
                    if v:
                        if d0 == -1:
                            d0 = j
                        elif d1 == -1:
                            d1 = j
                        else:
                            raise NotImplementedError("More than two DT clusters found")

                if sum([cd_id == -1 for cd_id in [d0, d1, c0, c1]]) > 2:
                    raise NotImplementedError("2+1 (or 1+2) clusters found")

                if tag[i] == 20:
                    p0, e0 = csc_phi[i][c0], csc_eta[i][c0]
                    p1, e1 = csc_phi[i][c1], csc_eta[i][c1]
                if tag[i] == 2:
                    p0, e0 = dt_phi[i][d0], dt_eta[i][d0]
                    p1, e1 = dt_phi[i][d1], dt_eta[i][d1]
                if tag[i] == 11:
                    p0, e0 = csc_phi[i][c0], csc_eta[i][c0]
                    p1, e1 = dt_phi[i][d0], dt_eta[i][d0]
                if tag[i] == 10:
                    p0, e0 = csc_phi[i][c0], csc_eta[i][c0]
                if tag[i] == 1:
                    p0, e0 = dt_phi[i][d0], dt_eta[i][d0]

                dphi[i] = min(2 * np.pi - np.abs(p0 - p1), np.abs(p0 - p1))
                deta[i] = np.abs(e0 - e1)
                dr[i] = np.sqrt(dphi[i] * dphi[i] + deta[i] * deta[i])
            return dphi, deta, dr

        tags = tags.split(",")
        ncsc, ndt = self["nCscRechitClusters"], self["nDtRechitClusters"]
        self["tag"] = (
            10 * ((ncsc == 1) & (ndt == 0)) * ("csc" in tags)
            + 1 * ((ncsc == 0) & (ndt == 1)) * ("dt" in tags)
            + 20 * ((ncsc == 2) & (ndt == 0)) * ("csccsc" in tags)
            + 2 * ((ncsc == 0) & (ndt == 2)) * ("dtdt" in tags)
            + 11 * ((ncsc == 1) & (ndt == 1)) * ("cscdt" in tags)
        )
        self.filter(self["tag"] > 0, system="evt", has_clusters=False)

        self["tag_dPhi"] = _delta_kin(self["cscRechitClusterPhi"], self["dtRechitClusterPhi"], self["sel_csc"], self["sel_dt"], self["tag"], "dphi")
        self["tag_dEta"] = _delta_kin(self["cscRechitClusterEta"], self["dtRechitClusterEta"], self["sel_csc"], self["sel_dt"], self["tag"], "deta")
        self["tag_dR"] = np.sqrt(self["tag_dPhi"]**2 + self["tag_dEta"]**2)

        self.cut = pcut
        return self

    ################################
    #             Cuts             #
    ################################

    def cut_hlt(self, invert=False, **kwargs):
        self.cut, pcut = False, self.cut
        # sel = self['HLT_CscCluster_Loose'] | self['HLT_L1CSCCluster_DTCluster50']
        sel = self["HLT_L1CSCCluster_DTCluster50"]
        if invert:
            sel = ~sel
        self.filter(sel, system="evt", **kwargs)
        self.cut = pcut
        return self

    def cut_jet(self, system="csc,dt", csc_pt=10, dt_pt=10, invert=False, **kwargs):
        self.cut, pcut = False, self.cut
        if "csc" in system:
            sel = self["cscRechitClusterJetVetoPt"] < csc_pt
            if invert:
                sel_t = ~sel_t
            self = self.filter(sel, system="csc", **kwargs)
        if "dt" in system:
            sel = self["dtRechitClusterJetVetoPt"] < dt_pt
            if invert:
                sel_t = ~sel_t
            self = self.filter(sel, system="dt", **kwargs)
        self.cut = pcut
        return self
    
    def cut_muon(self, system="csc,dt", csc_pt=30, dt_pt=50, invert=False, **kwargs):
        self.cut, pcut = False, self.cut
        if "csc" in system:
            sel = self["cscRechitClusterMuonVetoPt"] < csc_pt
            # sel = sel | (self["cscRechitClusterMuonVetoGlobal"] == 0)
            if invert:
                sel_t = ~sel_t
            self = self.filter(sel, system="csc", **kwargs)
            
        if "dt" in system:
            sel = self["dtRechitClusterMuonVetoPt"] < dt_pt
            # sel = sel | (self["dtRechitClusterMuonVetoLooseId"] == 0)
            if invert:
                sel_t = ~sel_t
            self = self.filter(sel, system="dt", **kwargs)
        self.cut = pcut
        return self

    def cut_time(self, system="csc,dt", invert=False, cut_csc_spread=True, cut_rpc_hits=True, **kwargs):
        self.cut, pcut = False, self.cut
        min_csc_t, max_csc_t = -5, 12.5
        if "csc" in system:
            sel_t = (min_csc_t < self["cscRechitClusterTimeWeighted"]) & (
                self["cscRechitClusterTimeWeighted"] < max_csc_t
            )
            if cut_csc_spread:
                sel_t = sel_t & (self["cscRechitClusterTimeSpreadWeightedAll"] < 20)
            if invert:
                sel_t = ~sel_t
            self.filter(sel_t, system="csc", **kwargs)

        if "dt" in system:
            sel_t = self["dtRechitCluster_match_RPCBx_dPhi0p5"] == 0
            if cut_rpc_hits:
                sel_t = sel_t & (self["dtRechitCluster_match_RPChits_dPhi0p5"] > 0)
            if invert:
                sel_t = ~sel_t
            self.filter(sel_t, system="dt", **kwargs)

        self.cut = pcut
        return self

    def cut_l1(self, invert=False, **kwargs):
        self.cut, pcut = False, self.cut

        csc_z, csc_size = np.abs(self["cscRechitClusterZ"]), self["cscRechitClusterSize"]
        csc_r = np.sqrt(self["cscRechitClusterX"] ** 2 + self["cscRechitClusterY"] ** 2)
        first_in_plateau = (
              ((100 < csc_r) & (csc_r < 275) & (580 < csc_z) & (csc_z < 632) & (500 <= csc_size)) # ME 11
            | ((275 < csc_r) & (csc_r < 465) & (668 < csc_z) & (csc_z < 724) & (200 <= csc_size)) # ME 12
            | ((505 < csc_r) & (csc_r < 700) & (668 < csc_z) & (csc_z < 724) & (200 <= csc_size)) # ME 13
            #
            | ((139 < csc_r) & (csc_r < 345) & (789 < csc_z) & (csc_z < 850) & (500 <= csc_size)) # ME 21
            | ((357 < csc_r) & (csc_r < 700) & (791 < csc_z) & (csc_z < 850) & (200 <= csc_size)) # ME 22
            #
            | ((160 < csc_r) & (csc_r < 345) & (915 < csc_z) & (csc_z < 970) & (500 <= csc_size)) # ME 31
            | ((357 < csc_r) & (csc_r < 700) & (911 < csc_z) & (csc_z < 970) & (200 <= csc_size)) # ME 32
            #
            | ((178 < csc_r) & (csc_r < 345) & (1002 < csc_z) & (csc_z < 1063) & (500 <= csc_size)) # ME 41
            | ((357 < csc_r) & (csc_r < 700) & (1002 < csc_z) & (csc_z < 1063) & (200 <= csc_size)) # ME 42
        )
        self.filter(first_in_plateau, system="csc", **kwargs)

        self.cut = pcut
        return self

    def cut_halo(self, version="dtPhi", invert=False, **kwargs):
        self.cut, pcut = False, self.cut

        if version == "dtPhi":
            abs_dtphi = np.abs(self["dtRechitClusterPhi"])
            sel, system = ((abs_dtphi > 0.4) & (abs_dtphi < np.pi - 0.4)), "dt"
        elif version == "dPhi":
            sel, system = self["tag_dPhi"] > 0.5, "evt"
        else:
            raise ValueError(f"No halo cut of version='{version}'.")

        if invert:
            sel = ~sel

        self.filter(sel, system=system, **kwargs)

        self.cut = pcut
        return self

################################################################
################################################################
################################################################

class MuonSystemsAnalyzer:
    """Analyzer for the MuonSystem.
    Holds methods for plotting 1D and 2D histograms as well as
    important tables, graphs, and measurements."""

    def __init__(self, muon_systems) -> None:
        pass

    # title, xlabel, ylabel, binning, log=False
    def hist_1d(self, key, func):
        pass

    def hist_2d(self, xkey, ykey, xfunc):
        pass

    def plot_efficiency_1d(self):
        pass

    def plot_efficiency_2d(self):
        pass
