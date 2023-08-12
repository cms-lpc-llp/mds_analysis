"""Class holding functions for processing muon system ntuples
Tasks:
    - MuonSystemAwkward
        - continue developing the basic structure
    
    - MuonSystemRDF
        - I messed up 'self_out = self'
        - Won't propagate cuts in some functions
"""

__author__ = 'Paul Simmerling'
__email__ = 'psimmerl@caltech.edu'
__credits__ = ['Paul Simmerling', 'Christina Wang', 'Lisa Benato', 'Si Xie', 'Cristian Pena', 'Martin Kwok',]

import numpy as np

import ROOT as rt

import uproot as upr
import awkward as ak

from itertools import product as combs
from src.helper_functions import lnot, land, lor, lxor, asum, aabs, alert
from src.histo_utilities import std_color_list

###############################

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
    name = title.split(';')[0]
    if name == '':
        name = str(np.random.randint(0, 1e9))

    hh = rt.TH1F(name, title, *bins)
    hh.SetLineWidth(4)

    if isinstance(x, ak.Array):
        x = ak.flatten(x, -1)

    norm = 1 / len(x) if 'norm' in kwargs and kwargs['norm'] == True else 1
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

    hh = rt.TH2F(name, title, *bins)

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
        b.Draw('same')

    l = rt.TLatex()
    l.SetTextSize(0.08)
    l.SetTextAlign(12)
    l.SetTextColor(12)
    l.SetTextAngle(90)
    l.DrawLatex(230, ymax * 0.3, 'Steel')
    #l.DrawLatex(230, ymax * 0.3, 'HCAL')
    #l.DrawLatex(335, ymax * 0.3, 'Solenoid')
    #l.DrawLatex(390, ymax * 0.3, 'Steel')

    l2 = rt.TLatex()
    l2.SetTextSize(0.06)
    l2.SetTextColor(13)
    l2.SetTextAngle(90)
    l2.DrawLatex(780, ymax * 0.5, 'Beyond CMS')
    text = rt.TLatex()
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

    boxes.append(rt.TBox(xmin, ymin, 568, ymax))  # in front of ME11
    boxes.append(rt.TBox(632, ymin, 671, ymax))  # between ME11 and ME12
    boxes.append(rt.TBox(724, ymin, 789, ymax))  # between ME12 and station2
    boxes.append(rt.TBox(849, ymin, 911, ymax))  # between station2 and station3
    boxes.append(rt.TBox(970, ymin, 1002, ymax))  # between station3 and station4
    boxes.append(rt.TBox(1073, ymin, xmax, ymax))  # beyond CMS
    for b in boxes:
        b.SetFillColor(15)
        b.SetFillStyle(3001)
        b.Draw('same')

    l = rt.TLatex()
    l.SetTextSize(0.08)
    l.SetTextColor(12)
    l.SetTextAngle(90)
    l.DrawLatex(xmin + 80, ymax * 0.4, 'Steel')

    l2 = rt.TLatex()
    l2.SetTextSize(0.06)
    l2.SetTextColor(13)
    l2.SetTextAngle(90)
    l2.DrawLatex(1110, ymax * 0.5, 'Beyond CMS')
    text = rt.TLatex()
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

###############################

class MuonSystemAwkward:
    """Handler for working with muon system ntuples using an Awkward arrays"""

    def __init__(self, file_name, nev=None, is_mc=False, tree_name='MuonSystem', implicit: bool=True) -> None:
        """Initialize a new instance of MuonSystemAwkward"""
        print('Building MuonSystemAwkward')
        print('\t', file_name)
        self.file_name  = file_name
        self.tree_name = tree_name
        self.is_mc = is_mc
        self.nev = nev
        self.implicit = implicit

        self.ms_dict = {
            'HLTDecision' : [], #! Needs to be cut before loading
            'met' : [],
            'metPhi' : [],
            'weight' : [],
            ################################
            'nCscRechitClusters' : [],
            'nCscRings' : [],
            #
            'cscRechitClusterX' : [],
            'cscRechitClusterY' : [],
            'cscRechitClusterZ' : [],
            'cscRechitClusterSize' : [],
            'cscRechitClusterEta' : [],
            'cscRechitClusterPhi' : [],
            'cscRechitClusterTime' : [],
            'cscRechitClusterTimeSpread' : [],
            'cscRechitClusterTimeSpreadWeightedAll' : [],
            'cscRechitClusterTimeWeighted' : [],
            'cscRechitClusterMet_dPhi' : [],
            #
            'cscRechitClusterJetVetoLooseId' : [],
            'cscRechitClusterJetVetoPt' : [],
            'cscRechitClusterMuonVetoLooseId' : [],
            'cscRechitClusterMuonVetoPt' : [],
            #
            'cscRechitClusterNStation10' : [],
            'cscRechitClusterAvgStation10' : [],
            'cscRechitClusterMaxStation' : [],
            'cscRechitClusterMaxStationRatio' : [],
            #
            'cscRechitClusterNChamber' : [],
            'cscRechitClusterMaxChamber' : [],
            'cscRechitClusterMaxChamberRatio' : [],
            #
            # 'cscRechitClusterNRechitChamberMinus11' : [],
            # 'cscRechitClusterNRechitChamberMinus12' : [],
            # 'cscRechitClusterNRechitChamberMinus13' : [],
            # 'cscRechitClusterNRechitChamberMinus21' : [],
            # 'cscRechitClusterNRechitChamberMinus22' : [],
            # 'cscRechitClusterNRechitChamberMinus31' : [],
            # 'cscRechitClusterNRechitChamberMinus32' : [],
            # 'cscRechitClusterNRechitChamberMinus41' : [],
            # 'cscRechitClusterNRechitChamberMinus42' : [],
            # 'cscRechitClusterNRechitChamberPlus11' : [],
            # 'cscRechitClusterNRechitChamberPlus12' : [],
            # 'cscRechitClusterNRechitChamberPlus13' : [],
            # 'cscRechitClusterNRechitChamberPlus21' : [],
            # 'cscRechitClusterNRechitChamberPlus22' : [],
            # 'cscRechitClusterNRechitChamberPlus31' : [],
            # 'cscRechitClusterNRechitChamberPlus32' : [],
            # 'cscRechitClusterNRechitChamberPlus41' : [],
            # 'cscRechitClusterNRechitChamberPlus42' : [],
            #
            'cscRechitCluster_match_MB1Seg_0p4' : [],
            'cscRechitCluster_match_RB1_0p4' : [],
            'cscRechitCluster_match_RE12_0p4' : [],
            #
            'cscRechitCluster_match_gLLP' : [],
            'cscRechitCluster_match_gLLP_csc' : [],
            'cscRechitCluster_match_gLLP_e' : [],
            'cscRechitCluster_match_gLLP_eta' : [],
            'cscRechitCluster_match_gLLP_phi' : [],
            'cscRechitCluster_match_gLLP_decay_r' : [],
            'cscRechitCluster_match_gLLP_decay_z' : [],
            ################################
            'nDtRechitClusters' : [],
            'nDtRings' : [],
            #
            'dtRechitClusterX' : [],
            'dtRechitClusterY' : [],
            'dtRechitClusterZ' : [],
            'dtRechitClusterSize' : [],
            'dtRechitClusterEta' : [],
            'dtRechitClusterPhi' : [],
            'dtRechitCluster_match_RPCBx_dPhi0p5' : [],
            'dtRechitCluster_match_RPChits_dPhi0p5' : [],
            'dtRechitClusterMet_dPhi' : [],
            #
            'dtRechitClusterJetVetoLooseId' : [],
            'dtRechitClusterJetVetoPt' : [],
            'dtRechitClusterMuonVetoLooseId' : [],
            'dtRechitClusterMuonVetoPt' : [],
            #
            'dtRechitClusterNStation10' : [],
            'dtRechitClusterAvgStation10' : [],
            'dtRechitClusterMaxStation' : [],
            'dtRechitClusterMaxStationRatio' : [],
            #
            'dtRechitClusterNChamber' : [],
            'dtRechitClusterMaxChamber' : [],
            'dtRechitClusterMaxChamberRatio' : [],
            #
            # 'dtRechitClusterNHitStation1' : [],
            # 'dtRechitClusterNHitStation2' : [],
            # 'dtRechitClusterNHitStation3' : [],
            # 'dtRechitClusterNHitStation4' : [],
            # 'dtRechitClusterNOppositeSegStation1' : [],
            # 'dtRechitClusterNOppositeSegStation2' : [],
            # 'dtRechitClusterNOppositeSegStation3' : [],
            # 'dtRechitClusterNOppositeSegStation4' : [],
            # 'dtRechitClusterNSegStation1' : [],
            # 'dtRechitClusterNSegStation2' : [],
            # 'dtRechitClusterNSegStation3' : [],
            # 'dtRechitClusterNSegStation4' : [],
            #
            'dtRechitCluster_match_MB1Seg_0p4' : [],
            'dtRechitCluster_match_MB1Seg_0p5' : [],
            'dtRechitCluster_match_MB1hits_0p4' : [],
            'dtRechitCluster_match_MB1hits_0p5' : [],
            #
            'dtRechitCluster_match_gLLP' : [],
            'dtRechitCluster_match_gLLP_dt' : [],
            'dtRechitCluster_match_gLLP_e' : [],
            'dtRechitCluster_match_gLLP_eta' : [],
            'dtRechitCluster_match_gLLP_phi' : [],
            'dtRechitCluster_match_gLLP_decay_r' : [],
            'dtRechitCluster_match_gLLP_decay_z' : [],
            ################################
            'nGLLP' : [],
            'gLLP_e' : [],
            'gLLP_pt' : [],
            'gLLP_eta' : [],
            'gLLP_phi' : [],
            'gLLP_decay_vertex_r' : [],
            'gLLP_decay_vertex_x' : [],
            'gLLP_decay_vertex_y' : [],
            'gLLP_decay_vertex_z' : [],
            ################################
            'nJets' : [],
            'jetE' : [],
            'jetPt' : [],
            'jetEta' : [],
            'jetPhi' : [],
            ################################
            'nLeptons' : [],
            'lepE' : [],
            'lepPt' : [],
            'lepPhi' : [],
            'lepEta' : [],
            'lepPdgId' : [],
        }

        self.column_names = list(self.ms_dict.keys())
        self.column_names.remove('HLTDecision')

        hlt_columns = {
            'HLT_CaloMET60_DTCluster50': 562,
            'HLT_CaloMET60_DTClusterNoMB1S50': 563,
            'HLT_L1MET_DTCluster50': 564,
            'HLT_L1MET_DTClusterNoMB1S50': 565,
            'HLT_CscCluster_Loose': 566,
            'HLT_CscCluster_Medium': 567,
            'HLT_CscCluster_Tight': 568,
            'HLT_L1CSCCluster_DTCluster50': 569,
            'HLT_L1CSCCluser_DTCluster75': 570,
        }
        # self.ff_upr = upr.open(
        #     path=self.file_name + ':' + self.tree_name,
        #     array_cache='100 kB',
        # )
        # ftpath =  + ':' + self.tree_name
        with upr.open(path=self.file_name, array_cache='100 kB') as fupr:
            ms = fupr[tree_name]
            print('HERE1')
            # ms_ak = ms.arrays(self.column_names, entry_stop=self.nev)
            print('HERE2')
            for k in self.ms_dict:
                print('HERE3', k)
                if k == 'HLTDecision': #TODO: apply a cut here to reduce array size
                    continue
                    for k_hlt, i_hlt in hlt_columns.items():
                        self.ms_dict[k_hlt] = ms[k][:,i_hlt].array()

                    # self.ms_dict[k] = ms.arrays(k+'[:,566:571]', entry_stop=self.nev)
                else:
                    self.ms_dict[k] = ms[k].array()
                    # self.ms_dict[k] = ms_ak[k]

        # del ms_ak # Clean up a potentially (really) big array

        if len(self.ms_dict['weight']) != self.nev:
            self.nev = len(self.ms_dict['weight'])
            print(self.nev)
            alert(f'Extracted {self.nev=:,}', '!')


    def __getitem__(self, key):
        if ('csc' in key[:3] or 'dt' in key[:2]) and 'RechitCluster' not in key:
            key.replace('csc', 'cscRechitCluster')
            key.replace('dt', 'dtRechitCluster')

        for s in ('csc','dt','gllp','jet','lep'):
            if s == key[:len(s)]:
                return self.ms_dict[s][key]
        return self.ms_dict[key]

    def __setitem__(self, key, value):
        if ('csc' in key[:3] or 'dt' in key[:2]) and 'RechitCluster' not in key:
            key.replace('csc', 'cscRechitCluster')
            key.replace('dt', 'dtRechitCluster')

        self.ms_dict[key] = value

    def __copy__(self):
        pass

    def set_implicit(self, implicit: bool=True) -> bool:
        """Set's the cut/filter method for the class.
        
        args:
            implicit : set's filter to cut implicitly or return a new filtered copy of the class
            
        returns:
            class filther method"""
        self.implicit = implicit
        return self.implicit

    def filter(self, sel, system='event'):
        if system == 'event':
            flat_sel = sel
        else:
            flat_sel = np.sum(sel, axis=1) > 0

        for k, v in self.ms_dict.items():
            if self.implicit:
                self.ms_dict[k] = v[sel if k[:len(system)] == system else flat_sel]
            else: #TODO: Set up 
                # make new ms_dict -> copy MuonSystemAwkward into new object w/ new dict
                pass
        n_columns = {

        }
        if system == 'csc':
            self.ms_dict['nCscRechitClusters'] = np.sum(self.ms_dict['cscRechitClusterSize'] > 0, axis=1)
        elif system == 'dt':
            self.ms_dict['nDtRechitClusters'] = np.sum(self.ms_dict['dtRechitClusterSize'] > 0, axis=1)
        elif system == 'gllp':
            self.ms_dict['nGLLP'] = np.sum(self.ms_dict['gllp_pt'] > 0, axis=1)
        elif system == 'jet':
            self.ms_dict['nJets'] = np.sum(self.ms_dict['jetPt'] > 0, axis=1)
        elif system == 'lep':
            self.ms_dict['nLeptons'] = np.sum(self.ms_dict['lepPt'] > 0, axis=1)

    def f(self, sel, system='event'):
        self.filter(sel=sel, system=system)

##################################################
##################################################
##################################################

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
            self.rdf = rt.RDataFrame(tree_name, file_name)
            if rt.IsImplicitMTEnabled():
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
            return self.rdf.Take[rt.VecOps.RVec[int]](key)
        elif ctype == 'RVec<float>':
            return self.rdf.Take[rt.VecOps.RVec[float]](key)
        elif ctype == 'RVec<bool>':
            return self.rdf.Take[rt.VecOps.RVec[bool]](key)
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


##################################################
##################################################
##################################################

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
