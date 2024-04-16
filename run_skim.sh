#!/bin/bash

# pyroot=/home/psimmerl/.conda/envs/pyroot/bin/python
# pyroot=/home/psimmerl/.miniconda3/envs/pyroot/bin/python
pyroot=/home/psimmerl/mambaforge/envs/pyroot/bin/python

# ###########
# # CSC-CSC #
# ###########
# # L1
# echo "################################"
# echo "CSC-CSC - L1 Selections"
# $pyroot skim_csccsc.py l1 it lt200
# $pyroot skim_csccsc.py l1 it low
# $pyroot skim_csccsc.py l1 it high

# $pyroot skim_csccsc.py l1 oot lt200
# $pyroot skim_csccsc.py l1 oot low
# $pyroot skim_csccsc.py l1 oot high
# echo ""

# # Standard Cuts
# echo "################################"
# echo "CSC-CSC - Standard Selections"
# $pyroot skim_csccsc.py scs it lt200
# $pyroot skim_csccsc.py scs it low
# $pyroot skim_csccsc.py scs it high

# $pyroot skim_csccsc.py scs oot lt200
# $pyroot skim_csccsc.py scs oot low
# $pyroot skim_csccsc.py scs oot high
# echo ""

# # Randomly Optimized Cuts
# echo "################################"
# echo "CSC-CSC - Randomly Optimized Selections"
# $pyroot skim_csccsc.py ropt it lt200
# $pyroot skim_csccsc.py ropt it low
# $pyroot skim_csccsc.py ropt it high

# $pyroot skim_csccsc.py ropt oot lt200
# $pyroot skim_csccsc.py ropt oot low
# $pyroot skim_csccsc.py ropt oot high
# echo ""

##########
# CSC-DT #
##########
# L1
echo "################################"
echo "CSC-DT - L1 Selections"
$pyroot skim_cscdt.py l1 it lt200
$pyroot skim_cscdt.py l1 it low
$pyroot skim_cscdt.py l1 it high

$pyroot skim_cscdt.py l1 oot lt200
$pyroot skim_cscdt.py l1 oot low
$pyroot skim_cscdt.py l1 oot high
echo ""

# Standard Cuts
echo "################################"
echo "CSC-DT - Standard Selections"
$pyroot skim_cscdt.py scs it lt200
$pyroot skim_cscdt.py scs it low
$pyroot skim_cscdt.py scs it high

$pyroot skim_cscdt.py scs oot lt200
$pyroot skim_cscdt.py scs oot low
$pyroot skim_cscdt.py scs oot high
echo ""

# Randomly Optimized Cuts
echo "################################"
echo "CSC-DT - Randomly Optimized Selections"
$pyroot skim_cscdt.py ropt it lt200
$pyroot skim_cscdt.py ropt it low
$pyroot skim_cscdt.py ropt it high

$pyroot skim_cscdt.py ropt oot lt200
$pyroot skim_cscdt.py ropt oot low
$pyroot skim_cscdt.py ropt oot high
echo ""


