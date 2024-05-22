#!/bin/bash

# pyroot=/home/psimmerl/.conda/envs/pyroot/bin/python
# pyroot=/home/psimmerl/.miniconda3/envs/pyroot/bin/python
# pyroot=/home/psimmerl/mambaforge/envs/pyroot/bin/python
pyroot=/home/psimmerl/.miniforge3/envs/pyroot/bin/python

echo "###########"
echo "# CSC-CSC #"
echo "###########"
# echo ""
# echo "TIGHT CUTS"
# echo ""
# # $pyroot skim_csccsc.py tight it lt200 cutflow
# $pyroot skim_csccsc.py tight it low cutflow
# $pyroot skim_csccsc.py tight it high cutflow

# # $pyroot skim_csccsc.py tight oot lt200 cutflow
# $pyroot skim_csccsc.py tight oot low cutflow
# $pyroot skim_csccsc.py tight oot high cutflow


# # $pyroot skim_csccsc.py tightDNN it lt200 bkgMC_plusBeamHalo cutflow
# $pyroot skim_csccsc.py tightDNN it low bkgMC_plusBeamHalo cutflow
# $pyroot skim_csccsc.py tightDNN it high bkgMC_plusBeamHalo cutflow

# # $pyroot skim_csccsc.py tightDNN oot lt200 bkgMC_plusBeamHalo cutflow
# $pyroot skim_csccsc.py tightDNN oot low bkgMC_plusBeamHalo cutflow
# $pyroot skim_csccsc.py tightDNN oot high bkgMC_plusBeamHalo cutflow

echo ""
echo "ROPT CUTS"
echo ""
# $pyroot skim_csccsc.py ropt it lt200 cutflow
$pyroot skim_csccsc.py ropt it low cutflow
$pyroot skim_csccsc.py ropt it high cutflow

# $pyroot skim_csccsc.py ropt oot lt200 cutflow
$pyroot skim_csccsc.py ropt oot low cutflow
$pyroot skim_csccsc.py ropt oot high cutflow


# $pyroot skim_csccsc.py roptDNN it lt200 bkgMC_plusBeamHalo cutflow
$pyroot skim_csccsc.py roptDNN it low bkgMC_plusBeamHalo cutflow
$pyroot skim_csccsc.py roptDNN it high bkgMC_plusBeamHalo cutflow

# $pyroot skim_csccsc.py roptDNN oot lt200 bkgMC_plusBeamHalo cutflow
$pyroot skim_csccsc.py roptDNN oot low bkgMC_plusBeamHalo cutflow
$pyroot skim_csccsc.py roptDNN oot high bkgMC_plusBeamHalo cutflow

# echo ""
# echo ""
# echo "##########"
# echo "# CSC-DT #"
# echo "##########"
# echo ""
# echo "TIGHT CUTS"
# echo ""
# # $pyroot skim_cscdt.py tight it lt200 cutflow
# $pyroot skim_cscdt.py tight it low cutflow
# $pyroot skim_cscdt.py tight it high cutflow

# # $pyroot skim_cscdt.py tight oot lt200 cutflow
# $pyroot skim_cscdt.py tight oot low cutflow
# $pyroot skim_cscdt.py tight oot high cutflow

# # $pyroot skim_cscdt.py tightDNN it lt200 bkgMC_plusBeamHalo cutflow
# $pyroot skim_cscdt.py tightDNN it low bkgMC_plusBeamHalo cutflow
# $pyroot skim_cscdt.py tightDNN it high bkgMC_plusBeamHalo cutflow

# # $pyroot skim_cscdt.py tightDNN oot lt200 bkgMC_plusBeamHalo cutflow
# $pyroot skim_cscdt.py tightDNN oot low bkgMC_plusBeamHalo cutflow
# $pyroot skim_cscdt.py tightDNN oot high bkgMC_plusBeamHalo cutflow

# echo ""
# echo "ROPT CUTS"
# echo ""
# # $pyroot skim_cscdt.py ropt it lt200 cutflow
# $pyroot skim_cscdt.py ropt it low cutflow
# $pyroot skim_cscdt.py ropt it high cutflow

# # $pyroot skim_cscdt.py ropt oot lt200 cutflow
# $pyroot skim_cscdt.py ropt oot low cutflow
# $pyroot skim_cscdt.py ropt oot high cutflow

# # $pyroot skim_cscdt.py roptDNN it lt200 bkgMC_plusBeamHalo cutflow
# $pyroot skim_cscdt.py roptDNN it low bkgMC_plusBeamHalo cutflow
# $pyroot skim_cscdt.py roptDNN it high bkgMC_plusBeamHalo cutflow

# # $pyroot skim_cscdt.py roptDNN oot lt200 bkgMC_plusBeamHalo cutflow
# $pyroot skim_cscdt.py roptDNN oot low bkgMC_plusBeamHalo cutflow
# $pyroot skim_cscdt.py roptDNN oot high bkgMC_plusBeamHalo cutflow