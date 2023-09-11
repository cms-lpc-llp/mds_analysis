"""helper_functions.py

Various helper functions to make MDC LLP analysis a little
easier to understand."""

import numpy as np
import math
import ROOT as rt
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
                np.multiply(rechit_time, rechit_energy * rechit_energy) / np.sum(rechit_energy * rechit_energy)
            )
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
        name = "c" + "_" + str(np.random.randint(0, 10000))
    c = rt.TCanvas(name, name, width, height)
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


def alert(msg: str, c: str = "w", form: str = ""):
    # print('')
    color_prefix = {
        "w": "",  # White
        "r": "\033[91m",  # Red
        "g": "\033[92m",  # Green
        "y": "\033[93m",  # Yellow
        "b": "\033[94m",  # Blue
        "m": "\033[95m",  # Purple/Magenta
        "c": "\033[96m",  # Cyan
        "lg": "\033[97m",  # LightGray
        "k": "\033[98m",  # Black
    }[c]
    print(color_prefix, end="")

    if len(form) > 1:
        raise NotImplementedError(f"Form length was too long, single character only! (form='{form}')")

    if form == "":
        print(msg, end="")
    elif form == "-":
        print("+" + ("-" * (len(msg) + 2)) + "+")
        print("| " + msg + " |")
        print("+" + ("-" * (len(msg) + 2)) + "+", end="")
    elif form == "=":
        print("#" + ("=" * (len(msg) + 2)) + "#")
        print("# " + msg + " #")
        print("#" + ("=" * (len(msg) + 2)) + "#", end="")
    elif form == "#":
        print("#" * (len(msg) + 4))
        print("# " + msg + " #")
        print("#" * (len(msg) + 4), end="")
    elif form == "!":
        print("!" * (len(msg) + 4))
        print("! " + msg + " !")
        print("!" * (len(msg) + 4), end="")
    else:  # form == '-':
        print(form * (len(msg) + 4))
        print(form + " " + msg + " " + form)
        print(form * (len(msg) + 4), end="")

    print("\033[0m" if c != "w" else "")


class Table:
    """Table
    TODO: Add multiple row labels, add double col sep for row labels"""
    def __init__(self, cols, header=None, rv_decimals=2, **kwargs) -> None:
        self.header = header

        self.row_spacers = False
        if "row_spacers" in kwargs:
            self.row_spacers = kwargs["row_spacers"]
        self.spacers = []

        ###

        self.rl_join = " "
        if "rl_join" in kwargs:
            self.rl_join = kwargs["rl_join"]

        self.same_rlp = False
        if "same_rlp" in kwargs:
            self.same_rlp = kwargs["same_rlp"]

        self.rl_decimals = 2
        if "rl_decimals" in kwargs:
            self.rl_decimals = kwargs["rl_decimals"]

        self.row_label_labels = [""]
        if "rl_label" in kwargs:
            self.row_label_labels = kwargs["rl_label"]
            if not isinstance(self.row_label_labels, (list, tuple)):
                self.row_label_labels = [self.row_label_labels]
            else:
                self.row_label_labels = [ rll for rll in self.row_label_labels ] # arrays get passed as pointers
        self.row_labels = []
        self.row_label_pads = [len(l) for l in self.row_label_labels]

        ###

        self.rv_join = " +/- "
        if "rv_join" in kwargs:
            self.rv_join = kwargs["rv_join"]

        self.same_rvp = False
        if "same_rvp" in kwargs:
            self.same_rvp = kwargs["same_rvp"]

        self.rv_decimals = rv_decimals

        self.row_value_labels = [rvl for rvl in cols] # arrays get passed as pointers

        self.row_values = []
        self.row_value_pads = [len(l) for l in self.row_value_labels]

    def _to_str(self, val, decimals, join):
        if isinstance(val, (list, tuple, np.ndarray)):
            for iv, v in enumerate(val):
                val[iv] = self._to_str(val=v, decimals=decimals, join=join)
            val = join.join(val)
        elif isinstance(val, (int, np.integer)):
            val = f"{val:,}"
        elif isinstance(val, (float, np.floating)):
            val = f"{val:,.{decimals}f}"
        else:
            val = f"{val}"
        return val

    def add_row(self, label, values):
        values = [ v for v in values ] # apparently arrays get passed as pointers
        if not isinstance(label, (list, tuple, np.ndarray)):
            label = [ label ]
        else:
            label = [ l for l in label ] # apparently arrays get passed as pointers

        for il, lab in enumerate(label):
            label[il] = self._to_str(lab, self.rl_decimals, self.rl_join)
            self.row_label_pads[il] = max(len(label[il]), self.row_label_pads[il])
        self.row_labels.append(label)

        if self.same_rvp:
            self.row_label_pads = len(self.row_label_pads) * [ max(self.row_label_pads) ]

        for iv, val in enumerate(values):
            values[iv] = self._to_str(val, self.rv_decimals, self.rv_join)
            self.row_value_pads[iv] = max(len(values[iv]), self.row_value_pads[iv])
        self.row_values.append(values)

        if self.same_rvp:
            self.row_value_pads = len(self.row_value_pads) * [ max(self.row_value_pads) ]

    def add_spacer(self):
        self.spacers.append(len(self.row_labels))

    def sort(self, func): # TODO: do this before we change type to str
        table_idxs = np.argsort([float(func(rvs)) for rvs in self.row_values]).astype(int)
        self.row_labels = [self.row_labels[idx] for idx in table_idxs]
        self.row_values = [self.row_values[idx] for idx in table_idxs]

    def print_spacer(self, split_columns=True):
        lps, vps = self.row_label_pads, self.row_value_pads
        tp = (sum(lps) + 3 * max(len(lps)-1,0)) + 4 + (sum(vps) + 3 * max(len(vps)-1, 0))
        lps[0] += max(len(self.header) - tp, 0)
        if split_columns:
            print("+-" + "-+-".join([p*"-" for p in lps]) + "-++-" + "-+-".join([p*"-" for p in vps]) + "-+")
        else:
            print("+-" + "---".join([p*"-" for p in lps]) + "----" + "---".join([p*"-" for p in vps]) + "-+")

    def print_header(self):
        lps, vps = self.row_label_pads, self.row_value_pads
        lls, vls = self.row_label_labels, self.row_value_labels
        tp = (sum(lps) + 3 * max(len(lps)-1,0)) + 4 + (sum(vps) + 3 * max(len(vps)-1, 0))
        lps[0] += max(len(self.header) - tp, 0)
        if self.header is not None:
            self.print_spacer(split_columns=False)
            print(f"| {self.header:^{tp}} |")

        self.print_spacer(split_columns=True)
        lls, vls = [f"{l:^{p}}" for l, p in zip(lls, lps)], [f"{l:^{p}}" for l, p in zip(vls, vps)]
        print("| " + " | ".join(lls) + " || " + " | ".join(vls) + " |")
        self.print_spacer(split_columns=True)

    def print_row(self, index):
        if isinstance(index, str):
            index = self.row_labels.index(index)

        if index in self.spacers and index:
            self.print_spacer(split_columns=True)

        lps, vps = self.row_label_pads, self.row_value_pads
        rls, rvs = self.row_labels[index], self.row_values[index]
        tp = (sum(lps) + 3 * max(len(lps)-1,0)) + 4 + (sum(vps) + 3 * max(len(vps)-1, 0))
        lps[0] += max(len(self.header) - tp, 0)

        rls, rvs = [f"{l:<{p}}" for l, p in zip(rls, lps)], [f"{l:>{p}}" for l, p in zip(rvs, vps)]
        print("| " + " | ".join(rls) + " || " + " | ".join(rvs) + " |")


    def print(self):
        self.print_header()
        for ir in range(len(self.row_labels)):
            self.print_row(ir)
            if self.row_spacers:
                self.print_spacer(split_columns=True)

        if not self.row_spacers:
            self.print_spacer(split_columns=True)