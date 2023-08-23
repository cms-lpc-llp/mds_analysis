"""helper_functions.py

Various helper functions to make MDC LLP analysis a little
easier to understand."""

import numpy as np
import math
import ROOT
import awkward as ak

###########################
####                   ####
####   original code   ####
####                   ####
###########################

def getRecoTime(algorithm, rechit_cut, rechit_time, rechit_energy):
    # 0 is energy weighted, 1 is energy squared weighted, 2 is median
    rechit_energy = rechit_energy[np.logical_not(rechit_time == -666)]
    rechit_time = rechit_time[np.logical_not(rechit_time == -666)]
    rechit_time = rechit_time[rechit_energy > rechit_cut]
    rechit_energy = rechit_energy[rechit_energy > rechit_cut]
    assert len(rechit_time) == len(rechit_energy)
    if np.sum(rechit_energy) > 0.0 and len(rechit_time) > 0:
        if algorithm == 0:
            return np.sum(np.multiply(rechit_time, rechit_energy) / np.sum(rechit_energy))
        elif algorithm == 1:
            return np.sum(
                np.multiply(rechit_time, rechit_energy * rechit_energy) / np.sum(rechit_energy * rechit_energy))
        elif algorithm == 2:
            return np.median(rechit_time)
    else:
        return None

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]

def deltaPhi(p1, p2):
    """Computes delta phi, handling periodic limit conditions."""
    res = abs(p1 - p2)
    if res > math.pi:
        res -= 2 * math.pi
    return res

def deltaR(e1, p1, e2, p2):
    de = e1 - e2
    dp = deltaPhi(p1, p2)
    return math.sqrt(de * de + dp * dp)


################################
####                        ####
####   psimmerl additions   ####
####                        ####
################################

def weight_calc(llp_ct, new_ctau, old_ctau):
    if llp_ct.ndim > 1:
        llp_ct = np.array(np.sum(llp_ct, axis=1))
    source = np.exp(-1.0 * llp_ct / old_ctau) / old_ctau**2
    weight = 1.0 / new_ctau**2 * np.exp(-1.0 * llp_ct / new_ctau) / source
    return weight

def land(*args):
    out = args[0]
    for v in args[1:]:
        out = np.logical_and(out, v)
    return out

def lor(*args):
    out = args[0]
    for v in args[1:]:
        out = np.logical_or(out, v)
    return out

def lxor(*args):
    out = args[0]
    for v in args[1:]:
        out = np.logical_xor(out, v)
    return out

def lnot(arg):
    return np.logical_not(arg)

def asum(arg):
    return ak.sum(arg, axis=1)

def aabs(arg):
    return np.abs(arg)

def canvas(nrows=1, ncols=1, width=None, height=None, name=None, grid=True):
    if width is None:
        width = ncols * 1000
    if height is None:
        height = nrows * 800
    if name is None:
        name = 'c' + '_' + str(np.random.randint(0, 10000))
    c = ROOT.TCanvas(name, name, width, height)
    c.Clear()
    if ncols * nrows > 1:
        c.Divide(ncols, nrows)
        if grid:
            for i in range(nrows * ncols):
                c.cd(i + 1).SetGrid()
        c.cd(1)
    else:
        if grid:
            c.SetGrid()

    c.Draw()
    return c


def alert(msg: str, form: str='-', c: str='w'):
    # print('')
    color_prefix = {
        'w' : '', # White
        'r' : '\033[91m', # Red
        'g' : '\033[92m', # Green
        'y' : '\033[93m', # Yellow
        'b' : '\033[94m', # Blue
        'm' : '\033[95m', # Purple/Magenta
        'c' : '\033[96m', # Cyan
        'lg' : '\033[97m', # LightGray
        'k' : '\033[98m', # Black
    }[c]
    print(color_prefix, end='')

    if form == '!':
        print('!' * (len(msg) + 6))
        print('!! ' + msg + ' !!')
        print('!' * (len(msg) + 6) + ('\033[0m' if c!='w' else ''))
    elif form == '#':
        print('#' * (len(msg) + 6))
        print('#' * (len(msg) + 6))
        print('## ' + msg + ' ##')
        print('#' * (len(msg) + 6))
        print('#' * (len(msg) + 6) + ('\033[0m' if c!='w' else ''))
    elif form == '=':
        print('#' + ('=' * (len(msg) + 2)) + '#')
        print('# ' + msg + ' #')
        print('#' + ('=' * (len(msg) + 2)) + '#' + ('\033[0m' if c!='w' else ''))
    else:# form == '-':
        print('+' + ('-' * (len(msg) + 2)) + '+')
        print('| ' + msg + ' |')
        print('+' + ('-' * (len(msg) + 2)) + '+' + ('\033[0m' if c!='w' else ''))
