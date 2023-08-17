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

__author__ = 'Paul Simmerling'
__email__ = 'psimmerl@caltech.edu'
__credits__ = [
    'Paul Simmerling', 
    'Christina Wang', 
    'Lisa Benato', 
    'Si Xie', 
    'Cristian Pena', 
    'Martin Kwok', 
    'Pedro ',
]

import sys
import numpy as np
import numba as nb

# import ROOT as rt
from ROOT import IsImplicitMTEnabled, RDataFrame, TH1F, TH2F
from ROOT import TLatex, TLegend, TBox
from ROOT import (kRed, kBlue, kGreen, kCyan, kMagenta, kYellow, kBlack, kAzure)
from ROOT.VecOps import RVec

import uproot as upr
import awkward as ak
import awkward.numba

from itertools import product as combs
from src.helper_functions import lnot, land, lor, lxor, asum, aabs, alert
from src.histo_utilities import std_color_list

###############################

def get_lat_leg(leg_coords=(0.6, 0.7, 0.95, 0.8)):
    lat = TLatex()
    lat.SetTextColor(kRed)
    lat.SetTextSize(0.03)
    lat.SetTextAlign(33)

    leg = TLegend(*leg_coords)
    leg.SetTextSize(0.025)
    leg.SetBorderSize(0)
    leg.SetEntrySeparation(0.01)

    return lat, leg

def H1D(x, title, bins, **kwargs):
    name = title.split(';')[0]
    if name == '':
        name = str(np.random.randint(0, 1e9))

    hh = TH1F(name, title, *bins)
    hh.SetLineWidth(4)

    if isinstance(x, ak.Array):
        x = ak.flatten(x, -1)

    norm = 1 / len(x) if 'norm' in kwargs and kwargs['norm'] else 1
    # [hh.Fill(xx, norm) for xx in x]

    if 'lw' in kwargs:
        hh.SetLineWidth(kwargs['lw'])
    if 'c' in kwargs:
        hh.SetLineColor(kwargs['c'])
    if 'ymin' in kwargs:
        hh.SetMinimum(kwargs['ymin'])
    if 'ymax' in kwargs:
        hh.SetMaximum(kwargs['ymax'])

    return hh

def H2D(x, y, title, bins, method='c', **kwargs):
    name = title.split(';')[0]
    if name == '':
        name = str(np.random.randint(0, 1e9))

    hh = TH2F(name, title, *bins)

    if method == 's':  # SEQUENTIAL
        if isinstance(x, ak.Array):
            x = ak.flatten(x, -1)
        if isinstance(y, ak.Array):
            y = ak.flatten(y, -1)

        #rtnp.fill_hist(hh, x, y)
        for xx, yy in zip(x, y):
            hh.Fill(xx, yy)  # , 1 / len(x) if 'norm' in kwargs and kwargs['norm'] == True else 1)
    elif method == 'c':  # COMBINATORIAL
        if len(x) != len(y):
            raise ValueError(f'x and y don\'t have the same length! ({len(x):=}, {len(y):=}')

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
            #     print('here1')
            # else:
            #     xr = np.array(xr)
            #     print('here2')

            # if isinstance(yr, ak.Array):
            #     if len(yr) == 0:
            #         continue
            #     yr = ak.flatten(yr, -1)
            # else:
            #     yr = np.array(yr)

            # print('\t', xr, yr)
            for xx, yy in combs(xr, yr):
                hh.Fill(xx, yy)  # , 1 / len(x) if 'norm' in kwargs and kwargs['norm'] == True else 1)

    return hh

def multi_plot(hhs, tts, **kwargs):
    ccs = std_color_list
    hhs = [hh.Clone() for hh in hhs]
    if 'ccs' in kwargs:
        ccs = kwargs['ccs']

    lat, leg = get_lat_leg(kwargs['legxy'] if 'legxy' in kwargs else (0.6, 0.7, 0.85, 0.95))
    if 'norm' in kwargs and kwargs['norm']:
        for hh in hhs:
            hh.Scale(1 / (hh.Integral() if hh.Integral() else 1))

    ymax = max([hh.GetMaximum() * (kwargs['ymax_mult'] if 'ymax_mult' in kwargs else 1.05) for hh in hhs])
    # ymin = max([hh.GetMinimum() * (kwargs['ymin_mult'] if 'ymin_mult' in kwargs else 1) for hh in hhs])

    if 'ymax' in kwargs and isinstance(kwargs['ymax'], str):
        if 'log' in kwargs['ymax']:
            ymax = 10**(np.ceil(np.log10(ymax)))
    elif 'ymax' in kwargs:
        ymax = kwargs['ymax']

    # if 'ymin' in kwargs and isinstance(kwargs['ymin'], str):
    #     if 'log' in kwargs['ymin']:
    #         kwargs['ymin'] = 10**(np.floor(np.log10(ymin)))

    for hh, tt, cc in zip(hhs, tts, ccs):
        hh.SetMaximum(ymax)
        if 'ymin' in kwargs:
            hh.SetMinimum(kwargs['ymin'])

        hh.SetLineColor(cc)

        if 'lw' in kwargs:
            hh.SetLineWidth(kwargs['lw'])
        hh.Draw('same hist')

        leg.AddEntry(hh, tt, 'L')

    leg.Draw()
    return hhs, lat, leg

###############################

def draw_dt_r_boxes(hh):
    ymin, ymax = hh.GetMinimum(), hh.GetMaximum()
    xmin, xmax = hh.GetXaxis().GetXmin(), hh.GetXaxis().GetXmax()
    boxes = []

    # if xmin < 181.1:
    #     boxes.append(TBox(max(xmin, 181.1), ymin, 286.4, ymax))  # HCAL barrel
    # if xmin < 295:
    #     boxes.append(TBox(max(xmin, 295), ymin, 377.5, ymax))  # solenoid
    #     #boxes.append(TBox(max(xmin, 295), ymin, 380, ymax))  # solenoid
    # boxes.append(TBox(max(xmin, 380), ymin, 402, ymax))  # b4 MB1
    boxes.append(TBox(max(xmin, 180), ymin, 402, ymax))  # b4 MB1
    boxes.append(TBox(449, ymin, 490, ymax))  # b4 MB2
    boxes.append(TBox(533, ymin, 597, ymax))  # b4 MB3
    boxes.append(TBox(636, ymin, 700, ymax))  # b4 MB4
    boxes.append(TBox(738, ymin, xmax, ymax))  # beyond CMS

    for b in boxes:
        b.SetFillColor(15)
        b.SetFillStyle(3001)
        b.Draw('same')

    l = TLatex()
    l.SetTextSize(0.08)
    l.SetTextAlign(12)
    l.SetTextColor(12)
    l.SetTextAngle(90)
    l.DrawLatex(230, ymax * 0.3, 'Steel')
    #l.DrawLatex(230, ymax * 0.3, 'HCAL')
    #l.DrawLatex(335, ymax * 0.3, 'Solenoid')
    #l.DrawLatex(390, ymax * 0.3, 'Steel')

    l2 = TLatex()
    l2.SetTextSize(0.06)
    l2.SetTextColor(13)
    l2.SetTextAngle(90)
    l2.DrawLatex(780, ymax * 0.5, 'Beyond CMS')
    text = TLatex()
    text.SetTextSize(0.04)
    text.DrawLatex(400, ymax * 1.01, 'MB1')
    text.DrawLatex(490, ymax * 1.01, 'MB2')
    text.DrawLatex(600, ymax * 1.01, 'MB3')
    text.DrawLatex(700, ymax * 1.01, 'MB4')

    return boxes

def draw_csc_z_boxes(hh):
    ymin, ymax = hh.GetMinimum(), hh.GetMaximum()
    xmin, xmax = hh.GetXaxis().GetXmin(), hh.GetXaxis().GetXmax()
    boxes = []

    boxes.append(TBox(xmin, ymin, 568, ymax))  # in front of ME11
    boxes.append(TBox(632, ymin, 671, ymax))  # between ME11 and ME12
    boxes.append(TBox(724, ymin, 789, ymax))  # between ME12 and station2
    boxes.append(TBox(849, ymin, 911, ymax))  # between station2 and station3
    boxes.append(TBox(970, ymin, 1002, ymax))  # between station3 and station4
    boxes.append(TBox(1073, ymin, xmax, ymax))  # beyond CMS
    for b in boxes:
        b.SetFillColor(15)
        b.SetFillStyle(3001)
        b.Draw('same')

    l = TLatex()
    l.SetTextSize(0.08)
    l.SetTextColor(12)
    l.SetTextAngle(90)
    l.DrawLatex(xmin + 80, ymax * 0.4, 'Steel')

    l2 = TLatex()
    l2.SetTextSize(0.06)
    l2.SetTextColor(13)
    l2.SetTextAngle(90)
    l2.DrawLatex(1110, ymax * 0.5, 'Beyond CMS')
    text = TLatex()
    text.SetTextSize(0.04)
    text.DrawLatex(570, ymax * 1.01, 'ME1/1')
    text.DrawLatex(660, ymax * 1.01, 'ME1/2-3')
    text.DrawLatex(795, ymax * 1.01, 'ME2')
    text.DrawLatex(920, ymax * 1.01, 'ME3')
    text.DrawLatex(1015, ymax * 1.01, 'ME4')

    return boxes

def make_cluster_eff_1D(ms, det, xl='z', cuts=False):
    lcltr, lgllp = det + 'RechitCluster_match_gLLP_decay_', 'gLLP_decay_vertex_'
    title = f';gLLP {xl.upper()} Decay Vertex [cm];{det.upper()} efficiency'

    if xl == 'z':
        if det == 'csc':
            bins = (100, 400, 1100)
        elif det == 'dt':
            bins = (100, 0, 700)
    elif xl == 'r':
        if det == 'csc':
            bins = (100, 0, 800)
        elif det == 'dt':
            bins = (100, 180, 800)

    sel_gllp_in_det = ms.in_det_cut(det, 'gllp')
    sel_cltr_in_det = ms.in_det_cut(det, 'rechit')

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

###############################################################
###############################################################

class MuonSystemAwkward:
    """Handler for working with muon system ntuples using an Awkward arrays"""

    def __init__(self, file_name, nev=None, is_mc=False, tree_name='MuonSystem', implicit: bool=True) -> None:
        """Initialize a new instance of MuonSystemAwkward"""
        
        print(f'Building MuonSystemAwkward (\'{tree_name}\') -')
        print(f'  \'{file_name}\'')
        
        ##########

        self.file_name  = file_name
        self.tree_name = tree_name
        self.is_mc = is_mc
        self.nev = nev

        ##########
        
        self.implicit = implicit
        self.cut = True

        ##########
        
        self.hlt_info = {
            # 'HLT_CaloMET60_DTCluster50': 562,
            # 'HLT_CaloMET60_DTClusterNoMB1S50': 563,
            # 'HLT_L1MET_DTCluster50': 564,
            # 'HLT_L1MET_DTClusterNoMB1S50': 565,
            # 'HLT_CscCluster_Loose': 566,
            # 'HLT_CscCluster_Medium': 567,
            # 'HLT_CscCluster_Tight': 568,
            'HLT_L1CSCCluster_DTCluster50': 569,
            # 'HLT_L1CSCCluser_DTCluster75': 570,
        }

        ##########

        self.ms_read = { }
        self.ms = upr.open(path=self.file_name + ':' + self.tree_name)#, array_cache='100 kB')
        self.ms_read['sel_evt'] = self.ms['met'].array(entry_stop=self.nev) > 0

        if len(self.ms_read['sel_evt']) != self.nev:
            self.nev = len(self.ms_read['sel_evt'])
            alert(f'Extracted {self.nev=:,}', '!')

        self.hlt_aliases = {k : f'HLTDecision[:{self.nev if self.nev else ""},{v}]' for k, v in self.hlt_info.items()}
        self.ms_read['HLT_L1CSCCluster_DTCluster50'] = self.ms.arrays('HLT_L1CSCCluster_DTCluster50', aliases=self.hlt_aliases, entry_stop=self.nev)['HLT_L1CSCCluster_DTCluster50'] > 0

        self.ms_read['sel_csc'] = self.ms['cscRechitClusterSize'].array(entry_stop=self.nev) > 0
        self.ms_read['sel_dt'] = self.ms['dtRechitClusterSize'].array(entry_stop=self.nev) > 0
        self.ms_read['nCscRechitClusters'] = np.sum(self.ms_read['sel_csc'], axis=1)
        self.ms_read['nDtRechitClusters'] = np.sum(self.ms_read['sel_dt'], axis=1)

    def __getitem__(self, key):
        # hlt_aa = self.ms.arrays(list(self.hlt_aliases.keys()), , entry_stop=self.nev)
        # for k, v in self.hlt_info.items():
        #     self.ms_read[k] = hlt_aa[k]

        if isinstance(key, str):
            key = self._fix_key(key)

            if self.cut:
                sel = self.ms_read['sel_evt']
                if 'csc' in key[:3]:
                    sel = (sel & self.ms_read['sel_csc'])
                if 'dt' in key[:2]:
                    sel = (sel & self.ms_read['sel_dt'])

                if key in self.ms_read:
                    arr = self.ms_read[sel]
                else:
                    arr = self.ms[key].array(entry_stop=self.nev)
                    # arr = self.ms.arrays(key, cut=sel, aliases=self.hlt_aliases, entry_stop=self.nev)[key]
            else:
                if key in self.ms_read:
                    arr = self.ms_read[key]
                else:
                    arr = self.ms[key].array(entry_stop=self.nev)
                    # arr = self.ms.arrays(key, aliases=self.hlt_aliases, entry_stop=self.nev)[key]
            return arr

        else:
            alert('BROKEN - Filters MuonSystem', '!', 'r')
            imp = self.implicit
            self.set_implicit(False)
            sout = self.filter(sel=key, system='evt')
            self.set_implicit(imp)
            return sout

    def __setitem__(self, key, value):
        key = self._fix_key(key)
        self.ms_read[key] = value

    def __copy__(self):
        pass

    def _fix_key(self, key):
        if 'csc' in key[:3] and 'RechitCluster' not in key:
            key = 'cscRechitCluster' + key[3:]
        if 'dt' in key[:2] and 'RechitCluster' not in key:
            key = 'dtRechitCluster' + key[2:]

        if 'nCsc' == key:
            key = 'nCscRechitClusters'
        if 'nDt' == key:
            key = 'nDtRechitClusters'
        if 'nLep' == key:
            key = 'nLeptons'
        if 'nJet' == key:
            key = 'nJets'

        return key

    def _fix_n_column(self, system):
        if system == 'evt':
            return

        n_columns = {
            'csc' : ('nCscRechitClusters', 'sel_csc'),
            'dt' : ('nDtRechitClusters', 'sel_dt'),
            'gllp' : ('nGLLP', 'sel_gllp'),
            'jet' : ('nJets', 'sel_jet'),
            'lep' : ('nLeptons', 'sel_lep'),
        }
        n_key, test_col = n_columns[system]
        self.ms_read[n_key] = np.sum(self.ms_read[test_col], axis=1)

    def count(self):
        return np.sum(self.ms_read['sel_evt'])

    def set_implicit(self, implicit: bool=True) -> bool:
        """Set's the cut/filter method for the class.
        
        args:
            implicit : set's filter to cut implicitly or return a new filtered copy of the class
            
        returns:
            class filther method"""
        if not implicit:
            alert('Not loading HLT', '!', 'r')

        self.implicit = implicit
        return self.implicit

    def filter(self, sel, system='evt', has_clusters=True):
        self.ms_read[f'sel_{system}'] = (self.ms_read[f'sel_{system}'] & sel)
        # for k, v in self.ms_dict.items():
        #     if 'sel_' in k:
        #         continue
 
        #     if self.implicit:
        #         if system == 'evt' or system == k[:len(system)]:
        #             self.ms_dict[k] = v[sel]
        #     else: #TODO
        #         alert('Non-implicit behavior does not work (yet)')
        #         # make new ms_dict -> copy MuonSystemAwkward into new object w/ new dict
        #         pass

        if system != 'evt':
            self._fix_n_column(system=system)
            if has_clusters:
                self.filter(self.ms_read['nCscRechitClusters'] + self.ms_read['nCscRechitClusters'] > 0, system='evt')

        return self

    def f(self, sel, system='evt'):
        return self.filter(sel=sel, system=system)

    def match_mc(self, system='csc,dt', check_decay=True):
        self.cut, pcut = False, self.cut

        max_csc_eta, max_csc_r, min_csc_z, max_csc_z = 3, 800, 400, 1200
        min_dt_r, max_dt_r, max_dt_z = 200, 800, 700
        if 'csc' in system:
            self.filter(self['cscRechitCluster_match_gLLP'], system='csc')
            if check_decay:
                sel_decay = (
                    ( self['cscRechitCluster_match_gLLP_eta'] < max_csc_eta ) &
                    ( self['cscRechitCluster_match_gLLP_decay_r'] < max_csc_r ) &
                    ( min_csc_z < np.abs(self['cscRechitCluster_match_gLLP_decay_z']) ) &
                    ( np.abs(self['cscRechitCluster_match_gLLP_decay_z']) < max_csc_z )
                )
                self.filter(sel_decay, system='csc')

        if 'dt' in system:
            self.filter(self['dtRechitCluster_match_gLLP'], system='dt')
            if check_decay:
                sel_decay = (
                    ( min_dt_r < self['dtRechitCluster_match_gLLP_decay_r'] ) &
                    ( self['dtRechitCluster_match_gLLP_decay_r'] < max_dt_r ) &
                    ( np.abs(self['dtRechitCluster_match_gLLP_decay_z']) < max_dt_z )
                )
                self.filter(sel_decay, system='dt')

        self.cut = pcut
        return self

    def blind(self, region):
        self.cut, pcut = False, self.cut
        if region == 'dphi':
            self.filter(self['tag_dPhi'] > 0.5, system='evt')

        self.cut = pcut
        return self

    def tag(self, tags='csccsc,cscdt'):
        self.cut, pcut = False, self.cut
        @nb.njit
        def _delta_kins(csc_phi, csc_eta, dt_phi, dt_eta, met_phi, sel_csc, sel_dt, tag):
            dphi, deta, dr = np.zeros(len(tag)), np.zeros(len(tag)), np.zeros(len(tag)),
            for i in range(len(csc_eta)):
                if tag[i] == 0:
                    continue

                p0, p1, e0, e1 = 0, met_phi[i], 0, 0
                c0, c1, d0, d1 = -1, -1, -1, -1
                # Ew but works, (np.argwhere(sel_csc[i]) not compatible with compiled ArrayView(?))
                for j, v in enumerate(sel_csc[i]):
                    if v:
                        if c0 == -1:
                            c0 = j
                        else:
                            c1 = j
                # Ew but works, (np.argwhere(sel_dt[i]) not compatible with compiled ArrayView(?))
                for j, v in enumerate(sel_dt[i]):
                    if v:
                        if c0 == -1:
                            d0 = j
                        else:
                            d1 = j

                # if len(sel_csc[i]):
                #     cidxs = np.arange(len(sel_csc[i])) * sel_csc[i]
                # if len(sel_dt[i]):
                #     didxs = np.arange(len(sel_dt[i])) * sel_dt[i]

                if tag[i] == 20:
                    p0, e0 = csc_phi[i][c0], csc_eta[i][c0]
                    p1, e1 = csc_phi[i][c1], csc_eta[i][c1]
                if tag[i] ==  2:
                    p0, e0 = dt_phi[i][d0], dt_eta[i][d0]
                    p1, e1 = dt_phi[i][d1], dt_eta[i][d1]
                if tag[i] == 11:
                    p0, e0 = csc_phi[i][c0], csc_eta[i][c0]
                    p1, e1 = dt_phi[i][d0], dt_eta[i][d0]
                if tag[i] == 10:
                    p0, e0 = csc_phi[i][c0], csc_eta[i][c0]
                if tag[i] ==  1:
                    p0, e0 = dt_phi[i][d0], dt_eta[i][d0]

                dphi[i] = np.abs(p0 - p1)
                dphi[i] -= 2 * np.pi * (dphi[i] > 2 * np.pi)
                deta[i] = np.abs(e0 - e1)
                dr[i] = np.sqrt(dphi[i]*dphi[i] + deta[i]*deta[i])
            return dphi, deta, dr
            # return ak.JaggedArray(dPhi, dEta, dR)

        tags = tags.split(',')

        ncsc, ndt = self['nCscRechitClusters'], self['nDtRechitClusters']
        tag = (
            10 * ((ncsc == 1) & (ndt == 0)) * ('csc' in tags) +
             1 * ((ncsc == 0) & (ndt == 1)) * ('dt' in tags) +
            20 * ((ncsc == 2) & (ndt == 0)) * ('csccsc' in tags) +
             2 * ((ncsc == 0) & (ndt == 2)) * ('dtdt' in tags) +
            11 * ((ncsc == 1) & (ndt == 1)) * ('cscdt' in tags)
        )
        self.filter(tag > 0, system='evt')

        csc_phi, csc_eta = self['cscRechitClusterPhi'], self['cscRechitClusterEta']
        dt_phi, dt_eta = self['dtRechitClusterPhi'], self['dtRechitClusterEta']
        sel_csc, sel_dt, met_phi = self['sel_csc'], self['sel_dt'], self['metPhi']

        tag_dphi, tag_deta, tag_dr = _delta_kins(csc_phi, csc_eta, dt_phi, dt_eta, met_phi, sel_csc, sel_dt, tag)

        self['tag'] = tag
        self['tag_dPhi'] = tag_dphi
        self['tag_dEta'] = tag_deta
        self['tag_dR'] = tag_dr

        self.cut = pcut
        return self

    def cut_hlt(self, invert=False):
        self.cut, pcut = False, self.cut
        # self.filter(self['HLT_CscCluster_Loose'] | self['HLT_L1CSCCluster_DTCluster50'], system='evt')
        self.filter(self['HLT_L1CSCCluster_DTCluster50'], system='evt')
        self.cut = pcut
        return self

    def cut_time(self, system='csc,dt', invert=False, cut_csc_spread=True, cut_rpc_hits=True):
        self.cut, pcut = False, self.cut
        min_csc_t, max_csc_t = -5, 12.5
        if 'csc' in system:
            sel_t = (
                (min_csc_t < self['cscRechitClusterTimeWeighted']) &
                (self['cscRechitClusterTimeWeighted'] < max_csc_t)
            )
            if cut_csc_spread:
                sel_t = (sel_t & (self['cscRechitClusterTimeSpreadWeightedAll'] < 20))
            if invert:
                sel_t = not sel_t
            self.filter(sel_t, system='csc')

        if 'dt' in system:
            sel_t = (self['dtRechitCluster_match_RPCBx_dPhi0p5'] == 0)
            if cut_rpc_hits:
                sel_t = (sel_t & (self['dtRechitCluster_match_RPChits_dPhi0p5'] > 0))
            if invert:
                sel_t = not sel_t
            self.filter(sel_t, system='dt')

        self.cut = pcut
        return self

    def cut_l1(self, invert=False):
        self.cut, pcut = False, self.cut

        csc_z, csc_size =  self['cscRechitClusterZ'], self['cscRechitClusterSize']
        csc_r = np.sqrt(self['cscRechitClusterX']**2 + self['cscRechitClusterY']**2)
        first_in_plateau = (
            ( (csc_r>100) & (csc_r<275) & (np.abs(csc_z)>580) & (np.abs(csc_z)<632) & (csc_size>=500) ) |
            ( (csc_r>139) & (csc_r<345) & (np.abs(csc_z)>789) & (np.abs(csc_z)<850) & (csc_size>=500) ) |
            ( (csc_r>160) & (csc_r<345) & (np.abs(csc_z)>915) & (np.abs(csc_z)<970) & (csc_size>=500) ) |
            ( (csc_r>178) & (csc_r<345) & (np.abs(csc_z)>1002) & (np.abs(csc_z)<1063) & (csc_size>=500) ) |
            ( (csc_r>275) & (csc_r<465) & (np.abs(csc_z)>668) & (np.abs(csc_z)<724) & (csc_size>=200) ) |
            ( (csc_r>505) & (csc_r<700) & (np.abs(csc_z)>668) & (np.abs(csc_z)<724) & (csc_size>=200) ) |
            ( (csc_r>357) & (csc_r<700) & (np.abs(csc_z)>791) & (np.abs(csc_z)<850) & (csc_size>=200) ) |
            ( (csc_r>357) & (csc_r<700) & (np.abs(csc_z)>911) & (np.abs(csc_z)<970) & (csc_size>=200) ) |
            ( (csc_r>357) & (csc_r<700) & (np.abs(csc_z)>1002) & (np.abs(csc_z)<1063) & (csc_size>=200) )
        )
        self.filter(first_in_plateau, system='csc')

        self.cut = pcut
        return self




# class MuonSystemAwkward:
#     """Handler for working with muon system ntuples using an Awkward arrays"""

#     def __init__(self, file_name, nev=None, is_mc=False, tree_name='MuonSystem', implicit: bool=True) -> None:
#         """Initialize a new instance of MuonSystemAwkward"""
#         print(f'Building MuonSystemAwkward (\'{tree_name}\') -')
#         print(f'  \'{file_name}\'')
#         self.file_name  = file_name
#         self.tree_name = tree_name
#         self.is_mc = is_mc
#         self.nev = nev
#         self.implicit = implicit

#         self.ms_dict = {
#             'met' : ak.Array,
#             'metPhi' : ak.Array,
#             'runNum' : ak.Array,
#             'weight' : ak.Array,
#             ################################
#             # 'HLTDecision' : ak.Array, #! Needs to be cut first
#             'HLT_CaloMET60_DTCluster50' : ak.Array,
#             'HLT_CaloMET60_DTClusterNoMB1S50' : ak.Array,
#             'HLT_L1MET_DTCluster50' : ak.Array,
#             'HLT_L1MET_DTClusterNoMB1S50' : ak.Array,
#             'HLT_CscCluster_Loose' : ak.Array,
#             'HLT_CscCluster_Medium' : ak.Array,
#             'HLT_CscCluster_Tight' : ak.Array,
#             'HLT_L1CSCCluster_DTCluster50' : ak.Array,
#             'HLT_L1CSCCluser_DTCluster75' : ak.Array,
#             ################################
#             'nCscRechitClusters' : ak.Array,
#             'nCscRings' : ak.Array,
#             #
#             'cscRechitClusterX' : ak.Array,
#             'cscRechitClusterY' : ak.Array,
#             'cscRechitClusterZ' : ak.Array,
#             'cscRechitClusterSize' : ak.Array,
#             'cscRechitClusterEta' : ak.Array,
#             'cscRechitClusterPhi' : ak.Array,
#             'cscRechitClusterTime' : ak.Array,
#             'cscRechitClusterTimeSpread' : ak.Array,
#             'cscRechitClusterTimeSpreadWeightedAll' : ak.Array,
#             'cscRechitClusterTimeWeighted' : ak.Array,
#             'cscRechitClusterMet_dPhi' : ak.Array,
#             #
#             'cscRechitClusterJetVetoLooseId' : ak.Array,
#             'cscRechitClusterJetVetoPt' : ak.Array,
#             'cscRechitClusterMuonVetoLooseId' : ak.Array,
#             'cscRechitClusterMuonVetoPt' : ak.Array,
#             #
#             'cscRechitClusterNStation10' : ak.Array,
#             'cscRechitClusterAvgStation10' : ak.Array,
#             'cscRechitClusterMaxStation' : ak.Array,
#             'cscRechitClusterMaxStationRatio' : ak.Array,
#             #
#             'cscRechitClusterNChamber' : ak.Array,
#             'cscRechitClusterMaxChamber' : ak.Array,
#             'cscRechitClusterMaxChamberRatio' : ak.Array,
#             #
#             # 'cscRechitClusterNRechitChamberMinus11' : ak.Array,
#             # 'cscRechitClusterNRechitChamberMinus12' : ak.Array,
#             # 'cscRechitClusterNRechitChamberMinus13' : ak.Array,
#             # 'cscRechitClusterNRechitChamberMinus21' : ak.Array,
#             # 'cscRechitClusterNRechitChamberMinus22' : ak.Array,
#             # 'cscRechitClusterNRechitChamberMinus31' : ak.Array,
#             # 'cscRechitClusterNRechitChamberMinus32' : ak.Array,
#             # 'cscRechitClusterNRechitChamberMinus41' : ak.Array,
#             # 'cscRechitClusterNRechitChamberMinus42' : ak.Array,
#             # 'cscRechitClusterNRechitChamberPlus11' : ak.Array,
#             # 'cscRechitClusterNRechitChamberPlus12' : ak.Array,
#             # 'cscRechitClusterNRechitChamberPlus13' : ak.Array,
#             # 'cscRechitClusterNRechitChamberPlus21' : ak.Array,
#             # 'cscRechitClusterNRechitChamberPlus22' : ak.Array,
#             # 'cscRechitClusterNRechitChamberPlus31' : ak.Array,
#             # 'cscRechitClusterNRechitChamberPlus32' : ak.Array,
#             # 'cscRechitClusterNRechitChamberPlus41' : ak.Array,
#             # 'cscRechitClusterNRechitChamberPlus42' : ak.Array,
#             #
#             'cscRechitCluster_match_MB1Seg_0p4' : ak.Array,
#             'cscRechitCluster_match_RB1_0p4' : ak.Array,
#             'cscRechitCluster_match_RE12_0p4' : ak.Array,
#             #
#             'cscRechitCluster_match_gLLP' : ak.Array,
#             'cscRechitCluster_match_gLLP_csc' : ak.Array,
#             'cscRechitCluster_match_gLLP_e' : ak.Array,
#             'cscRechitCluster_match_gLLP_eta' : ak.Array,
#             'cscRechitCluster_match_gLLP_phi' : ak.Array,
#             'cscRechitCluster_match_gLLP_decay_r' : ak.Array,
#             'cscRechitCluster_match_gLLP_decay_z' : ak.Array,
#             ################################
#             'nDtRechitClusters' : ak.Array,
#             'nDtRings' : ak.Array,
#             #
#             'dtRechitClusterX' : ak.Array,
#             'dtRechitClusterY' : ak.Array,
#             'dtRechitClusterZ' : ak.Array,
#             'dtRechitClusterSize' : ak.Array,
#             'dtRechitClusterEta' : ak.Array,
#             'dtRechitClusterPhi' : ak.Array,
#             'dtRechitCluster_match_RPCBx_dPhi0p5' : ak.Array,
#             'dtRechitCluster_match_RPChits_dPhi0p5' : ak.Array,
#             'dtRechitClusterMet_dPhi' : ak.Array,
#             #
#             'dtRechitClusterJetVetoLooseId' : ak.Array,
#             'dtRechitClusterJetVetoPt' : ak.Array,
#             'dtRechitClusterMuonVetoLooseId' : ak.Array,
#             'dtRechitClusterMuonVetoPt' : ak.Array,
#             #
#             'dtRechitClusterNStation10' : ak.Array,
#             'dtRechitClusterAvgStation10' : ak.Array,
#             'dtRechitClusterMaxStation' : ak.Array,
#             'dtRechitClusterMaxStationRatio' : ak.Array,
#             #
#             'dtRechitClusterNChamber' : ak.Array,
#             'dtRechitClusterMaxChamber' : ak.Array,
#             'dtRechitClusterMaxChamberRatio' : ak.Array,
#             #
#             # 'dtRechitClusterNHitStation1' : ak.Array,
#             # 'dtRechitClusterNHitStation2' : ak.Array,
#             # 'dtRechitClusterNHitStation3' : ak.Array,
#             # 'dtRechitClusterNHitStation4' : ak.Array,
#             # 'dtRechitClusterNOppositeSegStation1' : ak.Array,
#             # 'dtRechitClusterNOppositeSegStation2' : ak.Array,
#             # 'dtRechitClusterNOppositeSegStation3' : ak.Array,
#             # 'dtRechitClusterNOppositeSegStation4' : ak.Array,
#             # 'dtRechitClusterNSegStation1' : ak.Array,
#             # 'dtRechitClusterNSegStation2' : ak.Array,
#             # 'dtRechitClusterNSegStation3' : ak.Array,
#             # 'dtRechitClusterNSegStation4' : ak.Array,
#             #
#             'dtRechitCluster_match_MB1Seg_0p4' : ak.Array,
#             'dtRechitCluster_match_MB1Seg_0p5' : ak.Array,
#             'dtRechitCluster_match_MB1hits_0p4' : ak.Array,
#             'dtRechitCluster_match_MB1hits_0p5' : ak.Array,
#             #
#             'dtRechitCluster_match_gLLP' : ak.Array,
#             'dtRechitCluster_match_gLLP_dt' : ak.Array,
#             'dtRechitCluster_match_gLLP_e' : ak.Array,
#             'dtRechitCluster_match_gLLP_eta' : ak.Array,
#             'dtRechitCluster_match_gLLP_phi' : ak.Array,
#             'dtRechitCluster_match_gLLP_decay_r' : ak.Array,
#             'dtRechitCluster_match_gLLP_decay_z' : ak.Array,
#             ################################
#             'nGLLP' : ak.Array,
#             'gLLP_e' : ak.Array,
#             'gLLP_pt' : ak.Array,
#             'gLLP_eta' : ak.Array,
#             'gLLP_phi' : ak.Array,
#             'gLLP_decay_vertex_r' : ak.Array,
#             'gLLP_decay_vertex_x' : ak.Array,
#             'gLLP_decay_vertex_y' : ak.Array,
#             'gLLP_decay_vertex_z' : ak.Array,
#             ################################
#             'nJets' : ak.Array,
#             'jetE' : ak.Array,
#             'jetPt' : ak.Array,
#             'jetEta' : ak.Array,
#             'jetPhi' : ak.Array,
#             ################################
#             'nLeptons' : ak.Array,
#             'lepE' : ak.Array,
#             'lepPt' : ak.Array,
#             'lepPhi' : ak.Array,
#             'lepEta' : ak.Array,
#             'lepPdgId' : ak.Array,
#         }

#         self.column_names = list(self.ms_dict.keys())
#         # self.column_names.remove('HLTDecision')

#         hlt_columns = {
#             'HLT_CaloMET60_DTCluster50': 562,
#             'HLT_CaloMET60_DTClusterNoMB1S50': 563,
#             'HLT_L1MET_DTCluster50': 564,
#             'HLT_L1MET_DTClusterNoMB1S50': 565,
#             'HLT_CscCluster_Loose': 566,
#             'HLT_CscCluster_Medium': 567,
#             'HLT_CscCluster_Tight': 568,
#             'HLT_L1CSCCluster_DTCluster50': 569,
#             'HLT_L1CSCCluser_DTCluster75': 570,
#         }
#         hlt_aliases = {k : f'HLTDecision[:,{v}]' for k,v in hlt_columns.items()}

#         pft = self.file_name + ':' + self.tree_name
#         with upr.open(path=pft, array_cache='100 kB') as fms:
#             ms_ak = fms.arrays(
#                 self.column_names,
#                 entry_stop=self.nev,
#                 aliases=hlt_aliases,
#             )

#             # ms_ak = ms.arrays(self.column_names, entry_stop=self.nev)
#             for k in self.ms_dict:
#                 # if k == 'HLTDecision': #TODO: apply a cut here to reduce array size
#                 #     alert('Not loading HLT', '!', 'r')
#                 #     # for k_hlt, i_hlt in hlt_columns.items():
#                 #     #     self.ms_dict[k_hlt] = ms[k][:,i_hlt].array(entry_stop=self.nev)

#                 #     # self.ms_dict[k] = fms.arrays(k+'[:,566:571]', entry_stop=self.nev)
#                 # else:
#                     # self.ms_dict[k] = ms_ak[k]
#                 # self.ms_dict[k] = fms[k].array(entry_stop=self.nev)
#                 self.ms_dict[k] = ms_ak[k]

#         self.ms_dict['HLT'] = self.ms_dict['HLT_CscCluster_Loose'] | self.ms_dict['HLT_L1CSCCluster_DTCluster50']

#         if len(self.ms_dict['weight']) != self.nev:
#             self.nev = len(self.ms_dict['weight'])
#             alert(f'ExtAracted {self.nev=:,}', '!')

#     def __getitem__(self, key):
#         if isinstance(key, str):
#             key = self._fix_key(key)
#             return self.ms_dict[key]
#         else:
#             alert('BROKEN - Filters MuonSystem', '!', 'r')
#             imp = self.implicit
#             self.set_implicit(False)
#             sout = self.filter(sel=key, system='evt')
#             self.set_implicit(imp)
#             return sout

#     def __setitem__(self, key, value):
#         key = self._fix_key(key)
#         self.ms_dict[key] = value

#     def __copy__(self):
#         pass

#     def _fix_key(self, key):
#         if 'csc' in key[:3] and 'RechitCluster' not in key:
#             key = 'cscRechitCluster' + key[3:]
#         if 'dt' in key[:2] and 'RechitCluster' not in key:
#             key = 'dtRechitCluster' + key[2:]

#         if 'nCsc' == key:
#             key = 'nCscRechitClusters'
#         if 'nDt' == key:
#             key = 'nDtRechitClusters'
#         if 'nLep' == key:
#             key = 'nLeptons'
#         if 'nJet' == key:
#             key = 'nJets'

#         return key

#     def _fix_n_column(self, system):
#         if system == 'evt':
#             return

#         n_columns = {
#             'csc' : ('nCscRechitClusters', 'cscRechitClusterSize'),
#             'dt' : ('nDtRechitClusters', 'dtRechitClusterSize'),
#             'gllp' : ('nGLLP', 'gllp_pt'),
#             'jet' : ('nJets', 'jetPt'),
#             'lep' : ('nLeptons', 'lepPt'),
#         }
#         n_key, test_col = n_columns[system]
#         self.ms_dict[n_key] = np.sum(self.ms_dict[test_col] > 0, axis=1)

#     def count(self):
#         return len(self.ms_dict['weight'])

#     def set_implicit(self, implicit: bool=True) -> bool:
#         """Set's the cut/filter method for the class.
        
#         args:
#             implicit : set's filter to cut implicitly or return a new filtered copy of the class
            
#         returns:
#             class filther method"""
#         if not implicit:
#             alert('Not loading HLT', '!', 'r')

#         self.implicit = implicit
#         return self.implicit

#     def filter(self, sel, system='evt', has_clusters=True):
#         for k, v in self.ms_dict.items():
#             if self.implicit:
#                 if system == 'evt' or system == k[:len(system)]:
#                     self.ms_dict[k] = v[sel]
#             else: #TODO
#                 alert('Non-implicit behavior does not work (yet)')
#                 # make new ms_dict -> copy MuonSystemAwkward into new object w/ new dict
#                 pass

#         if system != 'evt':
#             self._fix_n_column(system=system)
#             if has_clusters:
#                 self.filter(self.ms_dict['nCscRechitClusters'] + self.ms_dict['nCscRechitClusters'] > 0, system='evt')

#         return self

#     def f(self, sel, system='evt'):
#         return self.filter(sel=sel, system=system)

#     def match_mc(self, system='csc,dt', check_decay=True):
#         max_csc_eta, max_csc_r, min_csc_z, max_csc_z = 3, 800, 400, 1200
#         min_dt_r, max_dt_r, max_dt_z = 200, 800, 700
#         if 'csc' in system:
#             self.filter(self.ms_dict['cscRechitCluster_match_gLLP'], system='csc')
#             if check_decay:
#                 sel_decay = (
#                     ( self.ms_dict['cscRechitCluster_match_gLLP_eta'] < max_csc_eta ) &
#                     ( self.ms_dict['cscRechitCluster_match_gLLP_decay_r'] < max_csc_r ) &
#                     ( min_csc_z < np.abs(self.ms_dict['cscRechitCluster_match_gLLP_decay_z']) ) &
#                     ( np.abs(self.ms_dict['cscRechitCluster_match_gLLP_decay_z']) < max_csc_z )
#                 )
#                 self.filter(sel_decay, system='csc')

#         if 'dt' in system:
#             self.filter(self.ms_dict['dtRechitCluster_match_gLLP'], system='dt')
#             if check_decay:
#                 sel_decay = (
#                     ( min_dt_r < self.ms_dict['dtRechitCluster_match_gLLP_decay_r'] ) &
#                     ( self.ms_dict['dtRechitCluster_match_gLLP_decay_r'] < max_dt_r ) &
#                     ( np.abs(self.ms_dict['dtRechitCluster_match_gLLP_decay_z']) < max_dt_z )
#                 )
#                 self.filter(sel_decay, system='dt')

#         return self

#     def blind(self, region):
#         if region == 'dphi':
#             self.filter(self.ms_dict['tag_dPhi'] > 0.5, system='evt')

#         return self

#     def tag(self, tags='csccsc,cscdt'):
#         @nb.njit
#         def _delta_kins(cscPhi, cscEta, dtPhi, dtEta, metPhi, tag):
#             dPhi, dEta, dR = np.zeros((len(tag),1)), np.zeros((len(tag),1)), np.zeros((len(tag),1))
#             for i in range(len(cscEta)):
#                 p0, p1, e0, e1 = 0, metPhi[i], 0, 0
#                 if tag[i] == 20:
#                     p0, e0 = cscPhi[i][0], cscEta[i][0]
#                     p1, e1 = cscPhi[i][1], cscEta[i][1]
#                 if tag[i] ==  2:
#                     p0, e0 = dtPhi[i][0], dtEta[i][0]
#                     p1, e1 = dtPhi[i][1], dtEta[i][1]
#                 if tag[i] == 11:
#                     p0, e0 = cscPhi[i][0], cscEta[i][0]
#                     p1, e1 = dtPhi[i][0], dtEta[i][0]
#                 if tag[i] == 10:
#                     p0, e0 = cscPhi[i][0], cscEta[i][0]
#                 if tag[i] ==  1:
#                     p0, e0 = dtPhi[i][0], dtEta[i][0]
#                 dPhi[i] = np.abs(p0 - p1)
#                 dPhi[i] -= 2 * np.pi * (dPhi[i] > 2 * np.pi)
#                 dEta[i] = np.abs(e0 - e1)
#                 dR[i] = np.sqrt(dPhi[i]*dPhi[i] + dEta[i]*dEta[i])
#             return dPhi, dEta, dR
#             # return ak.JaggedArray(dPhi, dEta, dR)

#         tags = tags.split(',')

#         ncsc, ndt = self.ms_dict['nCscRechitClusters'], self.ms_dict['nDtRechitClusters']
#         tag = ( #! Remember to filter no cluster and many cluster events out... for now :)
#             10 * ((ncsc == 1) & (ndt == 0)) * ('csc' in tags) +
#              1 * ((ncsc == 0) & (ndt == 1)) * ('dt' in tags) +
#             20 * ((ncsc == 2) & (ndt == 0)) * ('csccsc' in tags) +
#              2 * ((ncsc == 0) & (ndt == 2)) * ('dtdt' in tags) +
#             11 * ((ncsc == 1) & (ndt == 1)) * ('cscdt' in tags)
#         )

#         self.filter(tag > 0, system='evt')
#         ncsc, ndt, tag = ncsc[tag>0], ndt[tag>0], tag[tag>0]

#         cscPhi, cscEta = self.ms_dict['cscRechitClusterPhi'], self.ms_dict['cscRechitClusterEta']
#         dtPhi, dtEta = self.ms_dict['dtRechitClusterPhi'], self.ms_dict['dtRechitClusterEta']
#         metPhi = self.ms_dict['metPhi']

#         tag_dphi, tag_deta, tag_dr = _delta_kins(cscPhi, cscEta, dtPhi, dtEta, metPhi, tag)
#         self.ms_dict['tag'] = tag
#         self.ms_dict['tag_dPhi'] = tag_dphi
#         self.ms_dict['tag_dEta'] = tag_deta
#         self.ms_dict['tag_dR'] = tag_dr

#         return self

#     def cut_time(self, system='csc,dt', invert=False, cut_csc_spread=True, cut_rpc_hits=True):
#         min_csc_t, max_csc_t = -5, 12.5
#         if 'csc' in system:
#             sel_t = (
#                 (min_csc_t < self.ms_dict['cscRechitClusterTimeWeighted']) &
#                 (self.ms_dict['cscRechitClusterTimeWeighted'] < max_csc_t)
#             )
#             if cut_csc_spread:
#                 sel_t = (sel_t & (self.ms_dict['cscRechitClusterTimeSpreadWeightedAll'] < 20))
#             if invert:
#                 sel_t = not sel_t
#             self.filter(sel_t, system='csc')

#         if 'dt' in system:
#             sel_t = (self.ms_dict['dtRechitCluster_match_RPCBx_dPhi0p5'] == 1)
#             if cut_rpc_hits:
#                 sel_t = (sel_t & (self.ms_dict['dtRechitCluster_match_RPChits_dPhi0p5'] > 0))
#             if invert:
#                 sel_t = not sel_t
#             self.filter(sel_t, system='dt')
        
#         return self

#     def cut_l1(self, invert=False):
#         return self


###############################################################
###############################################################

#! BROKEN
class MuonSystemRDF:
    """Handler for working with muon system ntuples using an RDataFrame, works with 2tag CSC-DT (CSC triggers) only"""

    def __init__(self, file_name, tree_name='MuonSystem', is_mc=False, nev=None, rdf=None) -> None:
        if tree_name == 'MuonSystem':
            raise SystemError('I broke this class -- need to fix all of the \'self_out = self\' ')

        self.is_mc = is_mc
        self.file_name = file_name
        self.tree_name = tree_name
        self.nev = nev

        if rdf is None:
            self.rdf = RDataFrame(tree_name, file_name)
            if IsImplicitMTEnabled():
                print('Disabling \'Range\' for IMT -- will load ALL events!')
            else:
                if self.nev is not None:
                    self.rdf = self.rdf.Range(self.nev)
            self.Count()
        else:
            self.rdf = rdf

    def get(self, key):
        ctype = self.rdf.GetColumnType(key)
        # print(ctype)
        if ctype == 'int':
            return self.rdf.Take[int](key)
        elif ctype == 'float':
            return self.rdf.Take[float](key)
        elif ctype == 'bool':
            return self.rdf.Take[bool](key)
        elif ctype == 'RVec<int>':
            return self.rdf.Take[RVec[int]](key)
        elif ctype == 'RVec<float>':
            return self.rdf.Take[RVec[float]](key)
        elif ctype == 'RVec<bool>':
            return self.rdf.Take[RVec[bool]](key)
        else:
            raise ValueError(f'Column type \'{ctype}\' not known!')
        #return self.rdf.AsNumpy(key)

    def __getitem__(self, key):
        return self.get(key)

    def __setitem__(self, key, value):
        self.Define(key, value)

    def __copy__(self):
        pass

    def Define(self, key, value, implicit=True):
        if key in self.rdf.GetColumnNames():
            # print('In define - redefining -', key, value)
            rdf = self.rdf.Redefine(key, value)
            # print('In define - finished redefining')
        else:
            # print('In define - defining -', key, value)
            rdf = self.rdf.Define(key, value)
            # print('In define - finished defining')

        if implicit:
            self.rdf = rdf
        else:
            self_out = MuonSystemRDF(self.file_name, self.tree_name, self.is_mc, self.nev, rdf=rdf)

        return self_out

    def Filter(self, f, system='event', implicit=True):
        if system == 'event':
            # print('In filter - event -', f)
            rdf = self.rdf.Filter(f)
            # print('In filter - finished event filter')
            if implicit:
                self.rdf = rdf
                self_out = self
            else:
                self_out = MuonSystemRDF(self.file_name, self.tree_name, self.is_mc, self.nev, rdf=rdf)

        else:
            pre = {
                'csc': 'cscRechitCluster',
                'dt': 'dtRechitCluster',
                'gllp': 'gLLP',
                'lep': 'lep',
                'jet': 'jet',
            }

            if system not in pre:
                raise ValueError(f'Invaid system {system}.')

            system = pre[system]
            self_out = self.Define('cut', f, implicit=implicit)
            for k in self.rdf.GetColumnNames():
                k = str(k)
                if system == k[:len(system)]:
                    self_out = self.Define(k, f'{k}[cut]', implicit=True)

            self_out = self.fix_nbranch(implicit=True)
            self.rdf = self.rdf.Filter('nCscRechitClusters + nDtRechitClusters > 0')

        return self_out

    def Histo1D(self, *args):
        return self.rdf.Histo1D(*args)  #.GetPtr()

    def Histo2D(self, *args):
        #Add in behavior for RVec vs RVec, not sure how...
        xval, yval = args[-2:]
        if 'RVec' in self.rdf.GetColumnType(xval) and 'RVec' in self.rdf.GetColumnType(xval):
            print(f'Both \'{xval}\' and \'{yval}\' are RVectors! Will fail :(')
            # raise TypeError(f'Both '{xval}' and '{yval}' are RVectors! Will fail :(')
        return self.rdf.Histo2D(*args)  #.GetPtr()

    def Count(self):
        self.nev = self.rdf.Count().GetValue()
        return self.nev

    def fix_nbranch(self, implicit=True):
        self_out = self.Define('nCscRechitClusters', 'Sum(cscRechitClusterSize>0)', implicit=implicit)
        self_out = self.Define('nDtRechitClusters', 'Sum(dtRechitClusterSize>0)', implicit=True)
        self_out = self.Define('nJets', 'Sum(jetPt>0)', implicit=True)
        self_out = self.Define('nLeptons', 'Sum(lepPt>0)', implicit=True)
        return self_out

    def match_clusters(self, system='cscdt', in_det=True, implicit=True):
        if 'csc' in system:
            # self_out = self.Filter('(cscRechitCluster_match_gLLP == 1) && (cscRechitCluster_match_gLLP_minDeltaR < 0.4)', system='csc', implicit=implicit)
            self_out = self.Filter('cscRechitCluster_match_gLLP == 1', system='csc', implicit=implicit)
        if 'dt' in system:
            # self_out = self.Filter('(dtRechitCluster_match_gLLP == 1) && (dtRechitCluster_match_gLLP_minDeltaR < 0.4)', system='dt', implicit=implicit)
            self_out = self.Filter('dtRechitCluster_match_gLLP == 1', system='dt', implicit=implicit)

        if in_det:
            self_out = self.match_in_det(system=system, implicit=implicit)

        return self_out

    def match_in_det(self, system='cscdt', implicit=False):
        max_csc_eta, max_csc_r, min_csc_z, max_csc_z = 3, 800, 400, 1200
        min_dt_r, max_dt_r, max_dt_z = 200, 800, 700

        if 'csc' in system:
            self_out = self.Filter( \
                f'''(abs(cscRechitCluster_match_gLLP_eta) < {max_csc_eta}) && 
                    (abs(cscRechitCluster_match_gLLP_decay_r) < {max_csc_r}) && 
                    (abs(cscRechitCluster_match_gLLP_decay_z) > {min_csc_z}) && 
                    (abs(cscRechitCluster_match_gLLP_decay_z) < {max_csc_z})''',
                            system='csc',
                            implicit=implicit)

        if 'dt' in system:
            self_out = self.Filter( \
                f'''(abs(dtRechitCluster_match_gLLP_decay_r) > {min_dt_r}) && 
                    (abs(dtRechitCluster_match_gLLP_decay_r) < {max_dt_r}) && 
                    (abs(dtRechitCluster_match_gLLP_decay_z) < {max_dt_z})''',
                            system='dt',
                            implicit=implicit)

        return self_out

    # def blind(self, sector, implicit=True):
    #     if sector == 'dPhi':
    #         self_out = self.Filter()
    #     pass

    def find_csc_trigger(self, implicit=True):
        """Adds column with idx of the CSC cluster that may have triggered the event
            - Not in ME11

        """

        self_out = self.Define('idxCscRechitClusterTrigger',
                           '''
        int iTrig = -1;
        for (int iCsc = 0; iCsc < nCscRechitClusters; iCsc++) {
        }
        ''',
                           implicit=implicit)
        return self_out

# for k in tree_keys:
#     #OR over the two clusters
#     sel_trgCluster_tr1[k] = np.logical_and( cscClusterSize[k] >= 100, np.logical_and(cscClusterNStation[k]>=2, np.abs(cscClusterEta[k])<1.9))
#     sel_trgCluster_tr2[k] = np.logical_and( cscClusterSize[k] >= 200, np.logical_and(cscClusterNStation[k]==1, np.abs(cscClusterEta[k])<1.9))
#     sel_trgCluster_tr3[k] = np.logical_and( cscClusterSize[k] >= 500, np.abs(cscClusterEta[k])>=1.9)

#     sel_trgCluster_tr1_minus_size[k] = np.logical_and(cscClusterNStation[k]>=2, np.abs(cscClusterEta[k])<1.9)
#     sel_trgCluster_tr2_minus_size[k] = np.logical_and(cscClusterNStation[k]==1, np.abs(cscClusterEta[k])<1.9)
#     sel_trgCluster_tr3_minus_size[k] = np.abs(cscClusterEta[k])>=1.9

#     sel_HLT_OR[k] = np.logical_or(sel_trgCluster_tr1[k],np.logical_or(sel_trgCluster_tr2[k],sel_trgCluster_tr3[k]))

#     #Event level
#     L1_plateau[k] = ((cscClusterSize[k] >= 200).any()==True)
#     HLT_plateau[k] = np.logical_or( sel_trgCluster_tr1[k] , np.logical_or(sel_trgCluster_tr2[k],sel_trgCluster_tr3[k])  ).any()==True

#     #First cluster specific
#     first_in_ME11[k] = (cscClusterR[k][:,0]>100)&(cscClusterR[k][:,0]<275) &(abs(cscClusterZ[k][:,0])>580)&(abs(cscClusterZ[k][:,0])<632)
#     first_in_ME12[k] = (cscClusterR[k][:,0]>275)&(cscClusterR[k][:,0]<465) &(abs(cscClusterZ[k][:,0])>668)&(abs(cscClusterZ[k][:,0])<724)
#     first_in_ME13[k] = (cscClusterR[k][:,0]>505)&(cscClusterR[k][:,0]<700) &(abs(cscClusterZ[k][:,0])>668)&(abs(cscClusterZ[k][:,0])<724)

#     first_in_ME21[k] = (cscClusterR[k][:,0]>139)&(cscClusterR[k][:,0]<345) &(abs(cscClusterZ[k][:,0])>789)&(abs(cscClusterZ[k][:,0])<850)
#     first_in_ME22[k] = (cscClusterR[k][:,0]>357)&(cscClusterR[k][:,0]<700) &(abs(cscClusterZ[k][:,0])>791)&(abs(cscClusterZ[k][:,0])<850)

#     first_in_ME31[k] = (cscClusterR[k][:,0]>160)&(cscClusterR[k][:,0]<345) &(abs(cscClusterZ[k][:,0])>915)&(abs(cscClusterZ[k][:,0])<970)
#     first_in_ME32[k] = (cscClusterR[k][:,0]>357)&(cscClusterR[k][:,0]<700) &(abs(cscClusterZ[k][:,0])>911)&(abs(cscClusterZ[k][:,0])<970)

#     first_in_ME41[k] = (cscClusterR[k][:,0]>178)&(cscClusterR[k][:,0]<345) &(abs(cscClusterZ[k][:,0])>1002)&(abs(cscClusterZ[k][:,0])<1063)
#     first_in_ME42[k] = (cscClusterR[k][:,0]>357)&(cscClusterR[k][:,0]<700) &(abs(cscClusterZ[k][:,0])>1002)&(abs(cscClusterZ[k][:,0])<1063)

#     first_in_plateau_ME11[k] = first_in_ME11[k] & (cscClusterSize[k][:,0]>=500)
#     first_in_plateau_ME21[k] = first_in_ME21[k] & (cscClusterSize[k][:,0]>=500)
#     first_in_plateau_ME31[k] = first_in_ME31[k] & (cscClusterSize[k][:,0]>=500)
#     first_in_plateau_ME41[k] = first_in_ME41[k] & (cscClusterSize[k][:,0]>=500)

#     first_in_plateau_ME12[k] = first_in_ME12[k] & (cscClusterSize[k][:,0]>=200)
#     first_in_plateau_ME13[k] = first_in_ME13[k] & (cscClusterSize[k][:,0]>=200)
#     first_in_plateau_ME22[k] = (first_in_ME22[k]) & (cscClusterSize[k][:,0]>=200)
#     first_in_plateau_ME32[k] = first_in_ME32[k] & (cscClusterSize[k][:,0]>=200)
#     first_in_plateau_ME42[k] = first_in_ME42[k] & (cscClusterSize[k][:,0]>=200)

#     first_in_plateau[k] = first_in_plateau_ME11[k] | first_in_plateau_ME12[k] | first_in_plateau_ME13[k] | first_in_plateau_ME21[k] | first_in_plateau_ME22[k] | first_in_plateau_ME31[k] | first_in_plateau_ME32[k] | first_in_plateau_ME41[k] | first_in_plateau_ME42[k]

    def L1_plateau(self, implicit=True):
        """https://github.com/cms-lpc-llp/run3_muon_system_analysis/blob/main/study_triggered_events_CSC_v2.ipynb"""
        self_out = self.Filter('''
        auto sel_trgCluster_tr1 = (cscRechitClusterSize >= 100) && (cscRechitClusterNStation10 >= 2) && (abs(cscRechitClusterEta) < 1.9);
        auto sel_trgCluster_tr2 = (cscRechitClusterSize >= 200) && (cscRechitClusterNStation10 == 1) && (abs(cscRechitClusterEta) < 1.9);
        auto sel_trgCluster_tr3 = (cscRechitClusterSize >= 500) && (abs(cscRechitClusterEta) >= 1.9);

        // Event level
        auto L1_plateau = Sum(cscRechitClusterSize >= 200) > 0;
        auto HLT_plateau = Sum(sel_trgCluster_tr1 || sel_trgCluster_tr2 || sel_trgCluster_tr3) > 0;
        return (L1_plateau && HLT_plateau);
        ''',
                           implicit=implicit)

        return self_out

    def jet_cut(self, system='cscdt', implicit=True):  # AN-21-124
        if 'csc' in system:
            self_out = self.Filter(
                '!(cscRechitClusterJetVetoLooseId && (cscRechitClusterJetVetoPt > 30.) && (abs(cscRechitClusterEta) < 2.4))',
                system='csc',
                implicit=implicit)
        if 'dt' in system:
            self_out = self.Filter(
                '!(dtRechitClusterJetVetoLooseId && (dtRechitClusterJetVetoPt > 50.) && (abs(dtRechitClusterEta) < 2.4))',
                system='dt',
                implicit=implicit)
        return self_out

    def muon_cut(self, system='cscdt', implicit=True):  # AN-19-154
        if 'csc' in system:
            self_out = self.Filter(
                '!( (cscRechitClusterMuonVetoLooseId && (cscRechitClusterMuonVetoPt > 30.) && (abs(cscRechitClusterEta) < 2.4)) )',
                # '!( (cscRechitClusterMuonVetoLooseId && (cscRechitClusterMuonVetoPt > 30.) && (abs(cscRechitClusterEta) < 2.4)) || (cscRechitClusterNRechitChamberMinus11 + cscRechitClusterNRechitChamberMinus12 + cscRechitClusterNRechitChamberPlus11 + cscRechitClusterNRechitChamberPlus12 > 0) )',
                system='csc',
                implicit=implicit)
        if 'dt' in system:
            self_out = self.Filter(
                '!( (dtRechitClusterMuonVetoLooseId && (dtRechitClusterMuonVetoPt > 10.) && (abs(dtRechitClusterEta) < 2.4)) )',
                # '!( (dtRechitClusterMuonVetoLooseId && (dtRechitClusterMuonVetoPt > 10.) && (abs(dtRechitClusterEta) < 2.4)) || (dtRechitClusterNSegStation1 > 0) )',
                system='dt',
                implicit=True)
        return self_out

    def time_cut(self, time='it', system='cscdt', implicit=True):
        if time == 'oot':
            if 'cscT' in system:
                self_out = self.Filter(
                    'auto tcut = (cscRechitClusterTimeWeighted  < -12.5) || (cscRechitClusterTimeWeighted > 50); tcut[0] = 1; return tcut',
                    system='csc',
                    implicit=implicit)
            elif 'csc' in system:
                self_out = self.Filter('(cscRechitClusterTimeWeighted  < -12.5) || (cscRechitClusterTimeWeighted > 50)',
                                   system='csc',
                                   implicit=implicit)
            if 'dt' in system:
                self_out = self.Filter('dtRechitCluster_match_RPCBx_dPhi0p5 != 0', system='dt', implicit=implicit)
        elif time == 'it':
            if 'csc' in system:
                self_out = self.Filter('-5 < cscRechitClusterTimeWeighted && cscRechitClusterTimeWeighted < 12.5 && cscRechitClusterTimeSpreadWeightedAll < 20',
                                   system='csc',
                                   implicit=implicit)
            if 'dt' in system:
                self_out = self.Filter('(dtRechitCluster_match_RPCBx_dPhi0p5 == 0) && (dtRechitCluster_match_RPChits_dPhi0p5 > 0)', system='dt', implicit=implicit)
        return self_out

    def define_2tag_kins_and_cut(self, system='csccsc,cscdt', implicit=True):
        system_cut = []
        if 'csccsc' in system:
            system_cut.append('(nCscRechitClusters == 2)')
        if 'cscdt' in system:
            system_cut.append('(nCscRechitClusters == 1 && nDtRechitClusters == 1)')
        if '1csc' in system:
            system_cut.append('(nCscRechitClusters == 1)')
        if '1dt' in system:
            system_cut.append('(nDtRechitClusters == 1)')

        system_cut = ' || '.join(system_cut)
        self_out = self.Filter(system_cut, implicit=implicit)

        self_out = self.Define('vals',
                           '''
        double ttype = -999;
        double dEta = -999;
        double dPhi = -999;
        double dR = -999;
        if (nCscRechitClusters == 2) {
            ttype = 0;
            dEta = abs(cscRechitClusterEta[0] - cscRechitClusterEta[1]);
            dPhi = abs(cscRechitClusterPhi[0] - cscRechitClusterPhi[1]);
            dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);
        } else {
            ttype = 1;
            dEta = abs(cscRechitClusterEta[0] - dtRechitClusterEta[0]);
            dPhi = abs(cscRechitClusterPhi[0] - dtRechitClusterPhi[0]);
            dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);
        }
        return std::vector<double>{ttype,dEta,dPhi,dR}; 
        ''',
                           implicit=True)
        self_out = self.Define('tag_type', 'vals[0]', implicit=True)
        self_out = self.Define('tag_dEta', 'vals[1]', implicit=True)
        self_out = self.Define('tag_dPhi', 'vals[2]', implicit=True)
        self_out = self.Define('tag_dR', 'vals[3]', implicit=True)

        return self_out

    def print_cutflow_table(self, match=False):
        ms = self.Filter('met >= 0', implicit=False)  # Make a copy of the MuonSystem
        cscdt = '(nCscRechitClusters == 1) && (nDtRechitClusters == 1) && (nLeptons == 0)'
        me1veto = '(cscRechitClusterNRechitChamberMinus11 + cscRechitClusterNRechitChamberMinus12 + cscRechitClusterNRechitChamberPlus11 + cscRechitClusterNRechitChamberPlus12 == 0)'
        mb1veto = 'dtRechitClusterNSegStation1 == 0'
        # ##
        # print('Not implicit')

        print(r'\begin{table}[]')
        print(r'\begin{tabular}{c|rrr}')
        print(r'Selection & Yield & Eff. vs ' + ('matched' if match else 'no cut') + r' & Eff. vs 1CSC+1DT \\ \hline')

        yd = ms.Count()
        nall = yd
        print(f'All    & {yd:,} & ' + ('--' if match else f'{yd/nall*100:,.3f}\\%') + ' & -- \\\\')

        if match:
            #! THIS ONE IS IMPLICIT!!!!
            yd = ms.match_clusters(implicit=True).Filter('nCscRechitClusters + nDtRechitClusters > 0').Count()
            nall = yd
            print(f'Matched*  & {yd:,} & {yd/nall*100:,.3f}\\% & -- \\\\')
            #! !!!!!!!!!!!!!!!!!!!!!!!!

        #! THIS ONE IS IMPLICIT!!!!
        yd = ms.Filter(cscdt, implicit=True).Count()
        ntag = yd
        print(f'1CSC + 1DT + 0Lep*  & {yd:,} & {yd/nall*100:,.3f}\\% & {yd/ntag*100:,.3f}\\% \\\\')
        #! !!!!!!!!!!!!!!!!!!!!!!!!

        yd = ms.Filter(me1veto, 'csc', implicit=False).Filter(cscdt).Count()
        print(f'ME1 veto [CSC]  & {yd:,} & {yd/nall*100:,.3f}\\% & {yd/ntag*100:,.3f}\\% \\\\')

        yd = ms.Filter(mb1veto, 'dt', implicit=False).Filter(cscdt).Count()
        print(f'MB1 veto [DT]  & {yd:,} & {yd/nall*100:,.3f}\\% & {yd/ntag*100:,.3f}\\% \\\\')

        yd = ms.Filter(me1veto, 'csc', implicit=False).Filter(mb1veto, 'dt').Filter(cscdt).Count()
        print(f'ME1 + MB1 veto [CSC\\&DT]  & {yd:,} & {yd/nall*100:,.3f}\\% & {yd/ntag*100:,.3f}\\% \\\\')

        yd = ms.L1_plateau(implicit=False).Filter(cscdt).Count()
        print(f'L1 Plateau & {yd:,} & {yd/nall*100:,.3f}\\% & {yd/ntag*100:,.3f}\\% \\\\')

        yd = ms.L1_plateau(implicit=False).Filter(me1veto, 'csc').Filter(mb1veto, 'dt').Filter(cscdt).Count()
        print(f'L1 + ME1 + MB1 veto  & {yd:,} & {yd/nall*100:,.3f}\\% & {yd/ntag*100:,.3f}\\% \\\\')

        yd = ms.time_cut('it', 'csc', implicit=False).Filter(cscdt).Count()
        print(f'CSC In-time & {yd:,} & {yd/nall*100:,.3f}\\% & {yd/ntag*100:,.3f}\\% \\\\')

        yd = ms.time_cut('it', 'dt', implicit=False).Filter(cscdt).Count()
        print(f'DT In-time & {yd:,} & {yd/nall*100:,.3f}\\% & {yd/ntag*100:,.3f}\\% \\\\')

        yd = ms.time_cut('it', 'cscdt', implicit=False).Filter(cscdt).Count()
        print(f'CSC + DT In-time & {yd:,} & {yd/nall*100:,.3f}\\% & {yd/ntag*100:,.3f}\\% \\\\')

        yd = ms.jet_cut(implicit=False).Filter(cscdt).Count()
        print(f'Jet veto & {yd:,} & {yd/nall*100:,.3f}\\% & {yd/ntag*100:,.3f}\\% \\\\')

        yd = ms.muon_cut(implicit=False).Filter(cscdt).Count()
        print(f'Muon veto & {yd:,} & {yd/nall*100:,.3f}\\% & {yd/ntag*100:,.3f}\\% \\\\')

        yd = ms.jet_cut(implicit=False).muon_cut().Filter(cscdt).Count()
        print(f'Jet + Muon veto & {yd:,} & {yd/nall*100:,.3f}\\% & {yd/ntag*100:,.3f}\\% \\\\')

        print(r'\end{tabular}')
        # print(r'\caption{*This cut is applied to all consecutive rows}')
        print(r'\end{table}')
        print('')


###############################################################
###############################################################

# class MuonSystem:
#     """Handler for working with Muon System ntuples using uproot"""
#     _hlt_columns = {
#         'HLT_CaloMET60_DTCluster50': 562,
#         'HLT_CaloMET60_DTClusterNoMB1S50': 563,
#         'HLT_L1MET_DTCluster50': 564,
#         'HLT_L1MET_DTClusterNoMB1S50': 565,
#         'HLT_CscCluster_Loose': 566,
#         'HLT_CscCluster_Medium': 567,
#         'HLT_CscCluster_Tight': 568,
#         'HLT_L1CSCCluster_DTCluster50': 569,
#         'HLT_L1CSCCluser_DTCluster75': 570,
#     }

#     def __init__(self, file_name, tree_name='MuonSystem', isMC=False, nev=None) -> None:
#         self.isMC = isMC
#         self.file_name = file_name
#         self.tree_name = tree_name
#         self.nev = nev
#         self.fupr = upr.open(file_name + ':' + tree_name)

#         # self.keys = {'met', 'runNum'}
#         # self.ms = self.fupr.arrays(('met', 'runNum'), entry_stop=self.nev)
#         self.ms = {'met': self.fupr['met'].array(entry_stop=self.nev)}
#         self.cuts = []

#         #self.ms = self.fupr.arrays(array_cache='inherit',entry_stop=self.nev)

#     def __getitem__(self, key):
#         return self.get(key)

#     def __setitem__(self, key, value):
#         if isinstance(key, str):
#             self.ms[key] = value
#         else:
#             for k in self.ms:  # what was this supposed to be? set a slice?
#                 self.ms[k] = self.ms[k][value]
#         return self

#     def get(self, key):
#         # if key not in self.ms.fields:
#         #    dd = self.fupr.arrays(key, entry_stop=self.nev)[key]
#         if key not in self.ms:
#             if 'HLT' in key[:3]:
#                 if key == 'HLT' or key == 'HLTDecision':
#                     raise ValueError('Memory overflow if you try to load all of HLT.')
#                 if key[4:].isdigit():
#                     idx = int(key[4:])
#                 else:
#                     idx = self._hlt_columns[key]

#                 print(idx)

#                 dd = self.fupr['HLTDecision'].array(entry_stop=self.nev)
#                 #dd = self['HLTDecision']
#                 dd = dd[:, idx]
#                 #dd = self.fupr['HLTDecision'].array(filter_branch=lambda b: b[:,idx], entry_stop=self.nev)
#             else:
#                 dd = self.fupr[key].array(entry_stop=self.nev)

#             for system, idxs in self.cuts:
#                 if system == 'event':
#                     dd = dd[idxs]
#                 else:
#                     if system == key[:len(system)]:
#                         dd = dd[idxs]

#             self.ms[key] = dd

#             self.fix_nbranch()
#         return self.ms[key]

#     def apply_cut(self, idxs, system: str = 'event'):
#         system = system.lower()

#         if system == 'event':
#             for k in self.ms:
#                 self[k] = self[k][idxs]
#             # self.ms = self.ms[idxs]
#         else:
#             pre = {
#                 'csc': 'cscRechitCluster',
#                 'dt': 'dtRechitCluster',
#                 'gllp': 'gLLP',
#                 'lep': 'lep',
#                 'jet': 'jet',
#             }

#             if system not in pre:
#                 raise ValueError(f'Invaid system {system}.')

#             system = pre[system]
#             # for k in self.ms.fields:
#             for k in self.ms:
#                 if system == k[:len(system)]:
#                     self.ms[k] = self.ms[k][idxs]

#         self.cuts.append((system, idxs))
#         self.fix_nbranch()

#     def fix_nbranch(self):
#         self['nCscRechitClusters'] = asum(self.get('cscRechitClusterSize') > 0)
#         self['nDtRechitClusters'] = asum(self.get('dtRechitClusterSize') > 0)
#         self['nGLLP'] = asum(self.get('gLLP_pt') > 0)
#         self['nLeptons'] = asum(self.get('lepPt') > 0)
#         self['nJets'] = asum(self.get('jetPt') > 0)

#     #################################################
#     ## CUTS, returns idxs, implicit OFF by default ##
#     #################################################

#     def match_cut(self, dets=['csc', 'dt'], implicit=False):
#         if not self.isMC:
#             raise ValueError('Not Monte Carlo data, cannot match clusters!')

#         if not isinstance(dets, (tuple, list)):
#             dets = [dets]
#         cuts = []

#         for det in dets:
#             cuts.append(self.get(det + 'RechitCluster_match_gLLP'))

#         return cuts[0] if len(cuts) == 1 else cuts

#     def in_det_cut(self, dets=['csc', 'dt'], system='rechit', implicit=False):
#         max_csc_eta, max_csc_r, min_csc_z, max_csc_z = 3, 800, 400, 1200
#         min_dt_r, max_dt_r, max_dt_z = 200, 800, 700
#         if not isinstance(dets, (tuple, list)):
#             dets = [dets]
#         cuts = []

#         for det in dets:
#             if system == 'rechit':
#                 pree = det + 'RechitCluster_match_gLLP_'
#                 prev = pree + 'decay_'
#             elif system == 'gllp':
#                 pree = 'gLLP_'
#                 prev = pree + 'decay_vertex_'
#             if 'csc' == det:
#                 cuts.append(
#                     land(
#                         aabs(self.get(pree + 'eta')) < max_csc_eta,
#                         aabs(self.get(prev + 'r')) < max_csc_r,
#                         aabs(self.get(prev + 'z')) > min_csc_z,
#                         aabs(self.get(prev + 'z')) < max_csc_z,
#                     ))
#             if 'dt' == det:
#                 cuts.append(
#                     land(
#                         aabs(self.get(prev + 'r')) > min_dt_r,
#                         aabs(self.get(prev + 'r')) < max_dt_r,
#                         aabs(self.get(prev + 'z')) < max_dt_z,
#                     ))

#         return cuts[0] if len(cuts) == 1 else cuts

#     def hlt_cut(self, implicit=False):
#         """
#         562  HLT_CaloMET60_DTCluster50
#         563  HLT_CaloMET60_DTClusterNoMB1S50
#         564  HLT_L1MET_DTCluster50
#         565  HLT_L1MET_DTClusterNoMB1S50
#         566  HLT_CscCluster_Loose
#         567  HLT_CscCluster_Medium
#         568  HLT_CscCluster_Tight
#         569  HLT_L1CSCShower_DTCluster50
#         570  HLT_L1CSCShower_DTCluster75"""
#         return asum(self.get('HLTDecision')[:, 566:571]) > 0

#     def ndet_cut(self, ncsc=1, ndt=1, op='&', implicit=False):
#         """if ncsc or ndt is an iterable -> [inclusive, exclusive)
#         if you're not cutting on one of the vars LEAVE OP as default"""
#         lcsc, ldt = 'nCscRechitClusters', 'nDtRechitClusters'
#         if isinstance(ncsc, (tuple, list)):
#             idx_csc = land(self.get(lcsc) <= ncsc[0], self.get(lcsc) <= ncsc[1])
#         elif isinstance(ncsc, int):
#             idx_csc = self.get(lcsc) == ncsc
#         else:
#             idx_csc = np.ones_like(self.get('met'), dtype=bool)

#         if isinstance(ndt, (tuple, list)):
#             idx_dt = land(self.get(ldt) <= ndt[0], self.get(ldt) <= ndt[1])
#         elif isinstance(ndt, int):
#             idx_dt = self.get(ldt) == ndt
#         else:
#             idx_dt = np.ones_like(self.get('met'), dtype=bool)

#         if op == '&':
#             return land(idx_csc, idx_dt)
#         elif op == '|':
#             return lor(idx_csc, idx_dt)
#         elif op == '^':
#             return lxor(idx_csc, idx_dt)

#     def met_cut(self, met_min=50, met_max=None, implicit=False):
#         if met_min is None:
#             return self.get('met') < met_max
#         if met_max is None:
#             return met_min <= self.get('met')

#         return land(met_min <= self.get('met'), self.get('met') < met_max)

#     def muon_veto_cut(self, dets=['csc', 'dt'], implicit=False):
#         if not isinstance(dets, (tuple, list)):
#             dets = [dets]
#         cuts = []
#         for det in dets:
#             if 'csc' == det:
#                 cuts.append(
#                     lnot(
#                         lor(
#                             land(
#                                 self.get('cscRechitClusterMuonVetoLooseId'),
#                                 self.get('cscRechitClusterMuonVetoPt') > 30,
#                                 aabs(self.get('cscRechitClusterEta')) < 2.4,
#                             ),
#                             self.get('cscRechitClusterNRechitChamberMinus11') +
#                             self.get('cscRechitClusterNRechitChamberMinus12') +
#                             self.get('cscRechitClusterNRechitChamberPlus11') +
#                             self.get('cscRechitClusterNRechitChamberPlus12') > 0,
#                         )))  # AN-19-154
#             if 'dt' == det:
#                 cuts.append(
#                     lnot(
#                         lor(
#                             land(
#                                 self.get('dtRechitClusterMuonVetoLooseId'),
#                                 self.get('dtRechitClusterMuonVetoPt') > 10,
#                                 aabs(self.get('dtRechitClusterEta')) < 2.4,
#                             ),
#                             self.get('dtRechitClusterNSegStation1') > 0,
#                         )))  # AN-19-154
#         return cuts[0] if len(cuts) == 1 else cuts

#     def jet_veto_cut(self, dets=['csc', 'dt'], implicit=False):
#         if not isinstance(dets, (tuple, list)):
#             dets = [dets]

#         cuts = []
#         for det in dets:
#             if det == 'csc':
#                 cuts.append(
#                     lnot(
#                         land(
#                             self.get('cscRechitClusterJetVetoLooseId'),
#                             self.get('cscRechitClusterJetVetoPt') > 30,
#                             aabs(self.get('cscRechitClusterEta')) < 2.4,
#                         )))  # AN-21-124
#             if det == 'dt':
#                 cuts.append(
#                     lnot(
#                         land(
#                             self.get('dtRechitClusterJetVetoLooseId'),
#                             self.get('dtRechitClusterJetVetoPt') > 50,
#                             aabs(self.get('dtRechitClusterEta')) < 2.4,
#                         )))  # AN-21-124

#         return cuts[0] if len(cuts) == 1 else cuts

#     def time_cut(self, dets=['csc', 'dt'], implicit=False):
#         if not isinstance(dets, (tuple, list)):
#             dets = [dets]

#         cuts = []
#         for det in dets:
#             if 'csc' == det:
#                 cuts.append(
#                     land(
#                         -5 < self.get('cscRechitClusterTimeWeighted'),
#                         self.get('cscRechitClusterTimeWeighted') < 12.5,
#                     ))
#             if 'dt' == det:
#                 cuts.append(self.get('dtRechitCluster_match_RPCBx_dPhi0p5') == 0)

#         return cuts[0] if len(cuts) == 1 else cuts
