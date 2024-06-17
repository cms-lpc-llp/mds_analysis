#!/bin/bash

# pyroot=/home/psimmerl/.conda/envs/pyroot/bin/python
# pyroot=/home/psimmerl/.miniconda3/envs/pyroot/bin/python
pyroot=/home/psimmerl/mambaforge/envs/pyroot/bin/python
# pyroot=/home/psimmerl/.miniforge3/envs/pyroot/bin/python

echo "###########"
echo "# CSC-CSC #"
echo "###########"
echo ""
echo "+============================+"
echo "| CSC-CSC: TIGHT CUTS (2022) |"
echo "+============================+"
echo ""
# # $pyroot skim_csccsc.py tight it lt200 2022 #cutflow
# $pyroot skim_csccsc.py tight it low 2022 #cutflow
# $pyroot skim_csccsc.py tight it high 2022 #cutflow

# # $pyroot skim_csccsc.py tight oot lt200 2022 #cutflow
# $pyroot skim_csccsc.py tight oot low 2022 #cutflow
# $pyroot skim_csccsc.py tight oot high 2022 #cutflow

# $pyroot skim_csccsc.py tightDNN it lt200 bkgMC_plusBeamHalo 2022 #cutflow
$pyroot skim_csccsc.py tightDNN it low bkgMC_plusBeamHalo 2022 #cutflow
$pyroot skim_csccsc.py tightDNN it high bkgMC_plusBeamHalo 2022 #cutflow

# $pyroot skim_csccsc.py tightDNN oot lt200 bkgMC_plusBeamHalo 2022 #cutflow
$pyroot skim_csccsc.py tightDNN oot low bkgMC_plusBeamHalo 2022 #cutflow
$pyroot skim_csccsc.py tightDNN oot high bkgMC_plusBeamHalo 2022 #cutflow

echo ""
echo "+============================+"
echo "| CSC-CSC: TIGHT CUTS (2023) |"
echo "+============================+"
echo ""
# # $pyroot skim_csccsc.py tight it lt200 2023 #cutflow
# $pyroot skim_csccsc.py tight it low 2023 #cutflow
# $pyroot skim_csccsc.py tight it high 2023 #cutflow

# # $pyroot skim_csccsc.py tight oot lt200 2023 #cutflow
# $pyroot skim_csccsc.py tight oot low 2023 #cutflow
# $pyroot skim_csccsc.py tight oot high 2023 #cutflow

# $pyroot skim_csccsc.py tightDNN it lt200 bkgMC_plusBeamHalo 2023 #cutflow
$pyroot skim_csccsc.py tightDNN it low bkgMC_plusBeamHalo 2023 #cutflow
$pyroot skim_csccsc.py tightDNN it high bkgMC_plusBeamHalo 2023 #cutflow

# $pyroot skim_csccsc.py tightDNN oot lt200 bkgMC_plusBeamHalo 2023 #cutflow
$pyroot skim_csccsc.py tightDNN oot low bkgMC_plusBeamHalo 2023 #cutflow
$pyroot skim_csccsc.py tightDNN oot high bkgMC_plusBeamHalo 2023 #cutflow



echo ""
echo ""
echo "##########"
echo "# CSC-DT #"
echo "##########"
echo ""
echo "+===========================+"
echo "| CSC-DT: TIGHT CUTS (2022) |"
echo "+===========================+"
echo ""
# # $pyroot skim_cscdt.py tight it lt200 2022 #cutflow
# $pyroot skim_cscdt.py tight it low 2022 #cutflow
# $pyroot skim_cscdt.py tight it high 2022 #cutflow

# # $pyroot skim_cscdt.py tight oot lt200 2022 #cutflow
# $pyroot skim_cscdt.py tight oot low 2022 #cutflow
# $pyroot skim_cscdt.py tight oot high 2022 #cutflow

# $pyroot skim_cscdt.py tightDNN it lt200 bkgMC_plusBeamHalo 2022 #cutflow
$pyroot skim_cscdt.py tightDNN it low bkgMC_plusBeamHalo 2022 #cutflow
$pyroot skim_cscdt.py tightDNN it high bkgMC_plusBeamHalo 2022 #cutflow

# $pyroot skim_cscdt.py tightDNN oot lt200 bkgMC_plusBeamHalo 2022 #cutflow
$pyroot skim_cscdt.py tightDNN oot low bkgMC_plusBeamHalo 2022 #cutflow
$pyroot skim_cscdt.py tightDNN oot high bkgMC_plusBeamHalo 2022 #cutflow

echo ""
echo "+===========================+"
echo "| CSC-DT: TIGHT CUTS (2023) |"
echo "+===========================+"
echo ""
# # $pyroot skim_cscdt.py tight it lt200 2023 #cutflow
# $pyroot skim_cscdt.py tight it low 2023 #cutflow
# $pyroot skim_cscdt.py tight it high 2023 #cutflow

# # $pyroot skim_cscdt.py tight oot lt200 2023 #cutflow
# $pyroot skim_cscdt.py tight oot low 2023 #cutflow
# $pyroot skim_cscdt.py tight oot high 2023 #cutflow

# $pyroot skim_cscdt.py tightDNN it lt200 bkgMC_plusBeamHalo 2023 #cutflow
$pyroot skim_cscdt.py tightDNN it low bkgMC_plusBeamHalo 2023 #cutflow
$pyroot skim_cscdt.py tightDNN it high bkgMC_plusBeamHalo 2023 #cutflow

# $pyroot skim_cscdt.py tightDNN oot lt200 bkgMC_plusBeamHalo 2023 #cutflow
$pyroot skim_cscdt.py tightDNN oot low bkgMC_plusBeamHalo 2023 #cutflow
$pyroot skim_cscdt.py tightDNN oot high bkgMC_plusBeamHalo 2023 #cutflow
