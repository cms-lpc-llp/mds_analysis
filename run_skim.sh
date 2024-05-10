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

# ##########
# # CSC-DT #
# ##########
# # L1
# echo "################################"
# echo "CSC-DT - L1 Selections"
# $pyroot skim_cscdt.py l1 it lt200
# $pyroot skim_cscdt.py l1 it low
# $pyroot skim_cscdt.py l1 it high

# $pyroot skim_cscdt.py l1 oot lt200
# $pyroot skim_cscdt.py l1 oot low
# $pyroot skim_cscdt.py l1 oot high
# echo ""

# # Standard Cuts
# echo "################################"
# echo "CSC-DT - Standard Selections"
# $pyroot skim_cscdt.py scs it lt200
# $pyroot skim_cscdt.py scs it low
# $pyroot skim_cscdt.py scs it high

# $pyroot skim_cscdt.py scs oot lt200
# $pyroot skim_cscdt.py scs oot low
# $pyroot skim_cscdt.py scs oot high
# echo ""

# # Randomly Optimized Cuts
# echo "################################"
# echo "CSC-DT - Randomly Optimized Selections"
# $pyroot skim_cscdt.py ropt it lt200
# $pyroot skim_cscdt.py ropt it low
# $pyroot skim_cscdt.py ropt it high

# # $pyroot skim_cscdt.py ropt oot lt200
# $pyroot skim_cscdt.py ropt oot low
# $pyroot skim_cscdt.py ropt oot high
# echo ""

# ##############
# # CSC-DT OPT #
# ##############
# # LOO without DNN
# echo "################################"
# echo "CSC-DT - Opt without DNN"
# $pyroot skim_cscdt.py ropt oot lt200 rand
# $pyroot skim_cscdt.py ropt oot low rand
# $pyroot skim_cscdt.py ropt oot high rand

# $pyroot skim_cscdt.py ropt it lt200 rand
# $pyroot skim_cscdt.py ropt it low rand
# $pyroot skim_cscdt.py ropt it high rand
# echo ""


# echo "################################"
# echo "CSC-DT - Opt with DNN"
# $pyroot skim_cscdt.py roptDNN oot lt200 rand
# $pyroot skim_cscdt.py roptDNN oot low rand
# $pyroot skim_cscdt.py roptDNN oot high rand

# $pyroot skim_cscdt.py roptDNN it lt200 rand
# $pyroot skim_cscdt.py roptDNN it low rand
# $pyroot skim_cscdt.py roptDNN it high rand
# echo ""

#
$pyroot skim_cscdt.py l1 it lt200
$pyroot skim_cscdt.py l1 it low
$pyroot skim_cscdt.py l1 it high

$pyroot skim_cscdt.py l1 oot lt200
$pyroot skim_cscdt.py l1 oot low
$pyroot skim_cscdt.py l1 oot high

#
$pyroot skim_cscdt.py scs it lt200

$pyroot skim_cscdt.py scs it low noMB1Veto
mv data/processed/mc_cscdt_scs_low_rdf.root data/processed/mc_cscdt_scs_low_noMB1Veto_rdf.root
mv data/processed/r3_cscdt_scs_low_rdf.root data/processed/r3_cscdt_scs_low_noMB1Veto_rdf.root
$pyroot skim_cscdt.py scs it low noCSCJetVeto
mv data/processed/mc_cscdt_scs_low_rdf.root data/processed/mc_cscdt_scs_low_noCSCJetVeto_rdf.root
mv data/processed/r3_cscdt_scs_low_rdf.root data/processed/r3_cscdt_scs_low_noCSCJetVeto_rdf.root
$pyroot skim_cscdt.py scs it low noDTJetVeto
mv data/processed/mc_cscdt_scs_low_rdf.root data/processed/mc_cscdt_scs_low_noDTJetVeto_rdf.root
mv data/processed/r3_cscdt_scs_low_rdf.root data/processed/r3_cscdt_scs_low_noDTJetVeto_rdf.root
$pyroot skim_cscdt.py scs it low noCSCMuonVeto
mv data/processed/mc_cscdt_scs_low_rdf.root data/processed/mc_cscdt_scs_low_noCSCMuonVeto_rdf.root
mv data/processed/r3_cscdt_scs_low_rdf.root data/processed/r3_cscdt_scs_low_noCSCMuonVeto_rdf.root
$pyroot skim_cscdt.py scs it low noDTMuonVeto
mv data/processed/mc_cscdt_scs_low_rdf.root data/processed/mc_cscdt_scs_low_noDTMuonVeto_rdf.root
mv data/processed/r3_cscdt_scs_low_rdf.root data/processed/r3_cscdt_scs_low_noDTMuonVeto_rdf.root
$pyroot skim_cscdt.py scs it low noHaloVeto
mv data/processed/mc_cscdt_scs_low_rdf.root data/processed/mc_cscdt_scs_low_noHaloVeto_rdf.root
mv data/processed/r3_cscdt_scs_low_rdf.root data/processed/r3_cscdt_scs_low_noHaloVeto_rdf.root
$pyroot skim_cscdt.py scs it low

$pyroot skim_cscdt.py scs it high noMB1Veto
mv data/processed/mc_cscdt_scs_high_rdf.root data/processed/mc_cscdt_scs_high_noMB1Veto_rdf.root
mv data/processed/r3_cscdt_scs_high_rdf.root data/processed/r3_cscdt_scs_high_noMB1Veto_rdf.root
$pyroot skim_cscdt.py scs it high noCSCJetVeto
mv data/processed/mc_cscdt_scs_high_rdf.root data/processed/mc_cscdt_scs_high_noCSCJetVeto_rdf.root
mv data/processed/r3_cscdt_scs_high_rdf.root data/processed/r3_cscdt_scs_high_noCSCJetVeto_rdf.root
$pyroot skim_cscdt.py scs it high noDTJetVeto
mv data/processed/mc_cscdt_scs_high_rdf.root data/processed/mc_cscdt_scs_high_noDTJetVeto_rdf.root
mv data/processed/r3_cscdt_scs_high_rdf.root data/processed/r3_cscdt_scs_high_noDTJetVeto_rdf.root
$pyroot skim_cscdt.py scs it high noCSCMuonVeto
mv data/processed/mc_cscdt_scs_high_rdf.root data/processed/mc_cscdt_scs_high_noCSCMuonVeto_rdf.root
mv data/processed/r3_cscdt_scs_high_rdf.root data/processed/r3_cscdt_scs_high_noCSCMuonVeto_rdf.root
$pyroot skim_cscdt.py scs it high noDTMuonVeto
mv data/processed/mc_cscdt_scs_high_rdf.root data/processed/mc_cscdt_scs_high_noDTMuonVeto_rdf.root
mv data/processed/r3_cscdt_scs_high_rdf.root data/processed/r3_cscdt_scs_high_noDTMuonVeto_rdf.root
$pyroot skim_cscdt.py scs it high noHaloVeto
mv data/processed/mc_cscdt_scs_high_rdf.root data/processed/mc_cscdt_scs_high_noHaloVeto_rdf.root
mv data/processed/r3_cscdt_scs_high_rdf.root data/processed/r3_cscdt_scs_high_noHaloVeto_rdf.root
$pyroot skim_cscdt.py scs it high

$pyroot skim_cscdt.py scs oot lt200

$pyroot skim_cscdt.py scs oot low noMB1Veto
mv data/processed/mc_cscdtOOT_scs_low_rdf.root data/processed/mc_cscdtOOT_scs_low_noMB1Veto_rdf.root
mv data/processed/r3_cscdtOOT_scs_low_rdf.root data/processed/r3_cscdtOOT_scs_low_noMB1Veto_rdf.root
$pyroot skim_cscdt.py scs oot low noCSCJetVeto
mv data/processed/mc_cscdtOOT_scs_low_rdf.root data/processed/mc_cscdtOOT_scs_low_noCSCJetVeto_rdf.root
mv data/processed/r3_cscdtOOT_scs_low_rdf.root data/processed/r3_cscdtOOT_scs_low_noCSCJetVeto_rdf.root
$pyroot skim_cscdt.py scs oot low noDTJetVeto
mv data/processed/mc_cscdtOOT_scs_low_rdf.root data/processed/mc_cscdtOOT_scs_low_noDTJetVeto_rdf.root
mv data/processed/r3_cscdtOOT_scs_low_rdf.root data/processed/r3_cscdtOOT_scs_low_noDTJetVeto_rdf.root
$pyroot skim_cscdt.py scs oot low noCSCMuonVeto
mv data/processed/mc_cscdtOOT_scs_low_rdf.root data/processed/mc_cscdtOOT_scs_low_noCSCMuonVeto_rdf.root
mv data/processed/r3_cscdtOOT_scs_low_rdf.root data/processed/r3_cscdtOOT_scs_low_noCSCMuonVeto_rdf.root
$pyroot skim_cscdt.py scs oot low noDTMuonVeto
mv data/processed/mc_cscdtOOT_scs_low_rdf.root data/processed/mc_cscdtOOT_scs_low_noDTMuonVeto_rdf.root
mv data/processed/r3_cscdtOOT_scs_low_rdf.root data/processed/r3_cscdtOOT_scs_low_noDTMuonVeto_rdf.root
$pyroot skim_cscdt.py scs oot low noHaloVeto
mv data/processed/mc_cscdtOOT_scs_low_rdf.root data/processed/mc_cscdtOOT_scs_low_noHaloVeto_rdf.root
mv data/processed/r3_cscdtOOT_scs_low_rdf.root data/processed/r3_cscdtOOT_scs_low_noHaloVeto_rdf.root
$pyroot skim_cscdt.py scs oot low

$pyroot skim_cscdt.py scs oot high noMB1Veto
mv data/processed/mc_cscdtOOT_scs_high_rdf.root data/processed/mc_cscdtOOT_scs_high_noMB1Veto_rdf.root
mv data/processed/r3_cscdtOOT_scs_high_rdf.root data/processed/r3_cscdtOOT_scs_high_noMB1Veto_rdf.root
$pyroot skim_cscdt.py scs oot high noCSCJetVeto
mv data/processed/mc_cscdtOOT_scs_high_rdf.root data/processed/mc_cscdtOOT_scs_high_noCSCJetVeto_rdf.root
mv data/processed/r3_cscdtOOT_scs_high_rdf.root data/processed/r3_cscdtOOT_scs_high_noCSCJetVeto_rdf.root
$pyroot skim_cscdt.py scs oot high noDTJetVeto
mv data/processed/mc_cscdtOOT_scs_high_rdf.root data/processed/mc_cscdtOOT_scs_high_noDTJetVeto_rdf.root
mv data/processed/r3_cscdtOOT_scs_high_rdf.root data/processed/r3_cscdtOOT_scs_high_noDTJetVeto_rdf.root
$pyroot skim_cscdt.py scs oot high noCSCMuonVeto
mv data/processed/mc_cscdtOOT_scs_high_rdf.root data/processed/mc_cscdtOOT_scs_high_noCSCMuonVeto_rdf.root
mv data/processed/r3_cscdtOOT_scs_high_rdf.root data/processed/r3_cscdtOOT_scs_high_noCSCMuonVeto_rdf.root
$pyroot skim_cscdt.py scs oot high noDTMuonVeto
mv data/processed/mc_cscdtOOT_scs_high_rdf.root data/processed/mc_cscdtOOT_scs_high_noDTMuonVeto_rdf.root
mv data/processed/r3_cscdtOOT_scs_high_rdf.root data/processed/r3_cscdtOOT_scs_high_noDTMuonVeto_rdf.root
$pyroot skim_cscdt.py scs oot high noHaloVeto
mv data/processed/mc_cscdtOOT_scs_high_rdf.root data/processed/mc_cscdtOOT_scs_high_noHaloVeto_rdf.root
mv data/processed/r3_cscdtOOT_scs_high_rdf.root data/processed/r3_cscdtOOT_scs_high_noHaloVeto_rdf.root
$pyroot skim_cscdt.py scs oot high

#
$pyroot skim_cscdt.py ropt it low noMB1Veto
mv data/processed/mc_cscdt_ropt_low_rdf.root data/processed/mc_cscdt_ropt_low_noMB1Veto_rdf.root
mv data/processed/r3_cscdt_ropt_low_rdf.root data/processed/r3_cscdt_ropt_low_noMB1Veto_rdf.root
$pyroot skim_cscdt.py ropt it low noCSCJetVeto
mv data/processed/mc_cscdt_ropt_low_rdf.root data/processed/mc_cscdt_ropt_low_noCSCJetVeto_rdf.root
mv data/processed/r3_cscdt_ropt_low_rdf.root data/processed/r3_cscdt_ropt_low_noCSCJetVeto_rdf.root
$pyroot skim_cscdt.py ropt it low noDTJetVeto
mv data/processed/mc_cscdt_ropt_low_rdf.root data/processed/mc_cscdt_ropt_low_noDTJetVeto_rdf.root
mv data/processed/r3_cscdt_ropt_low_rdf.root data/processed/r3_cscdt_ropt_low_noDTJetVeto_rdf.root
$pyroot skim_cscdt.py ropt it low noCSCMuonVeto
mv data/processed/mc_cscdt_ropt_low_rdf.root data/processed/mc_cscdt_ropt_low_noCSCMuonVeto_rdf.root
mv data/processed/r3_cscdt_ropt_low_rdf.root data/processed/r3_cscdt_ropt_low_noCSCMuonVeto_rdf.root
$pyroot skim_cscdt.py ropt it low noDTMuonVeto
mv data/processed/mc_cscdt_ropt_low_rdf.root data/processed/mc_cscdt_ropt_low_noDTMuonVeto_rdf.root
mv data/processed/r3_cscdt_ropt_low_rdf.root data/processed/r3_cscdt_ropt_low_noDTMuonVeto_rdf.root
$pyroot skim_cscdt.py ropt it low noHaloVeto
mv data/processed/mc_cscdt_ropt_low_rdf.root data/processed/mc_cscdt_ropt_low_noHaloVeto_rdf.root
mv data/processed/r3_cscdt_ropt_low_rdf.root data/processed/r3_cscdt_ropt_low_noHaloVeto_rdf.root
$pyroot skim_cscdt.py ropt it low

$pyroot skim_cscdt.py ropt it high noMB1Veto
mv data/processed/mc_cscdt_ropt_high_rdf.root data/processed/mc_cscdt_ropt_high_noMB1Veto_rdf.root
mv data/processed/r3_cscdt_ropt_high_rdf.root data/processed/r3_cscdt_ropt_high_noMB1Veto_rdf.root
$pyroot skim_cscdt.py ropt it high noCSCJetVeto
mv data/processed/mc_cscdt_ropt_high_rdf.root data/processed/mc_cscdt_ropt_high_noCSCJetVeto_rdf.root
mv data/processed/r3_cscdt_ropt_high_rdf.root data/processed/r3_cscdt_ropt_high_noCSCJetVeto_rdf.root
$pyroot skim_cscdt.py ropt it high noDTJetVeto
mv data/processed/mc_cscdt_ropt_high_rdf.root data/processed/mc_cscdt_ropt_high_noDTJetVeto_rdf.root
mv data/processed/r3_cscdt_ropt_high_rdf.root data/processed/r3_cscdt_ropt_high_noDTJetVeto_rdf.root
$pyroot skim_cscdt.py ropt it high noCSCMuonVeto
mv data/processed/mc_cscdt_ropt_high_rdf.root data/processed/mc_cscdt_ropt_high_noCSCMuonVeto_rdf.root
mv data/processed/r3_cscdt_ropt_high_rdf.root data/processed/r3_cscdt_ropt_high_noCSCMuonVeto_rdf.root
$pyroot skim_cscdt.py ropt it high noDTMuonVeto
mv data/processed/mc_cscdt_ropt_high_rdf.root data/processed/mc_cscdt_ropt_high_noDTMuonVeto_rdf.root
mv data/processed/r3_cscdt_ropt_high_rdf.root data/processed/r3_cscdt_ropt_high_noDTMuonVeto_rdf.root
$pyroot skim_cscdt.py ropt it high noHaloVeto
mv data/processed/mc_cscdt_ropt_high_rdf.root data/processed/mc_cscdt_ropt_high_noHaloVeto_rdf.root
mv data/processed/r3_cscdt_ropt_high_rdf.root data/processed/r3_cscdt_ropt_high_noHaloVeto_rdf.root
$pyroot skim_cscdt.py ropt it high

$pyroot skim_cscdt.py ropt oot low noMB1Veto
mv data/processed/mc_cscdtOOT_ropt_low_rdf.root data/processed/mc_cscdtOOT_ropt_low_noMB1Veto_rdf.root
mv data/processed/r3_cscdtOOT_ropt_low_rdf.root data/processed/r3_cscdtOOT_ropt_low_noMB1Veto_rdf.root
$pyroot skim_cscdt.py ropt oot low noCSCJetVeto
mv data/processed/mc_cscdtOOT_ropt_low_rdf.root data/processed/mc_cscdtOOT_ropt_low_noCSCJetVeto_rdf.root
mv data/processed/r3_cscdtOOT_ropt_low_rdf.root data/processed/r3_cscdtOOT_ropt_low_noCSCJetVeto_rdf.root
$pyroot skim_cscdt.py ropt oot low noDTJetVeto
mv data/processed/mc_cscdtOOT_ropt_low_rdf.root data/processed/mc_cscdtOOT_ropt_low_noDTJetVeto_rdf.root
mv data/processed/r3_cscdtOOT_ropt_low_rdf.root data/processed/r3_cscdtOOT_ropt_low_noDTJetVeto_rdf.root
$pyroot skim_cscdt.py ropt oot low noCSCMuonVeto
mv data/processed/mc_cscdtOOT_ropt_low_rdf.root data/processed/mc_cscdtOOT_ropt_low_noCSCMuonVeto_rdf.root
mv data/processed/r3_cscdtOOT_ropt_low_rdf.root data/processed/r3_cscdtOOT_ropt_low_noCSCMuonVeto_rdf.root
$pyroot skim_cscdt.py ropt oot low noDTMuonVeto
mv data/processed/mc_cscdtOOT_ropt_low_rdf.root data/processed/mc_cscdtOOT_ropt_low_noDTMuonVeto_rdf.root
mv data/processed/r3_cscdtOOT_ropt_low_rdf.root data/processed/r3_cscdtOOT_ropt_low_noDTMuonVeto_rdf.root
$pyroot skim_cscdt.py ropt oot low noHaloVeto
mv data/processed/mc_cscdtOOT_ropt_low_rdf.root data/processed/mc_cscdtOOT_ropt_low_noHaloVeto_rdf.root
mv data/processed/r3_cscdtOOT_ropt_low_rdf.root data/processed/r3_cscdtOOT_ropt_low_noHaloVeto_rdf.root
$pyroot skim_cscdt.py ropt oot low

$pyroot skim_cscdt.py ropt oot high noMB1Veto
mv data/processed/mc_cscdtOOT_ropt_high_rdf.root data/processed/mc_cscdtOOT_ropt_high_noMB1Veto_rdf.root
mv data/processed/r3_cscdtOOT_ropt_high_rdf.root data/processed/r3_cscdtOOT_ropt_high_noMB1Veto_rdf.root
$pyroot skim_cscdt.py ropt oot high noCSCJetVeto
mv data/processed/mc_cscdtOOT_ropt_high_rdf.root data/processed/mc_cscdtOOT_ropt_high_noCSCJetVeto_rdf.root
mv data/processed/r3_cscdtOOT_ropt_high_rdf.root data/processed/r3_cscdtOOT_ropt_high_noCSCJetVeto_rdf.root
$pyroot skim_cscdt.py ropt oot high noDTJetVeto
mv data/processed/mc_cscdtOOT_ropt_high_rdf.root data/processed/mc_cscdtOOT_ropt_high_noDTJetVeto_rdf.root
mv data/processed/r3_cscdtOOT_ropt_high_rdf.root data/processed/r3_cscdtOOT_ropt_high_noDTJetVeto_rdf.root
$pyroot skim_cscdt.py ropt oot high noCSCMuonVeto
mv data/processed/mc_cscdtOOT_ropt_high_rdf.root data/processed/mc_cscdtOOT_ropt_high_noCSCMuonVeto_rdf.root
mv data/processed/r3_cscdtOOT_ropt_high_rdf.root data/processed/r3_cscdtOOT_ropt_high_noCSCMuonVeto_rdf.root
$pyroot skim_cscdt.py ropt oot high noDTMuonVeto
mv data/processed/mc_cscdtOOT_ropt_high_rdf.root data/processed/mc_cscdtOOT_ropt_high_noDTMuonVeto_rdf.root
mv data/processed/r3_cscdtOOT_ropt_high_rdf.root data/processed/r3_cscdtOOT_ropt_high_noDTMuonVeto_rdf.root
$pyroot skim_cscdt.py ropt oot high noHaloVeto
mv data/processed/mc_cscdtOOT_ropt_high_rdf.root data/processed/mc_cscdtOOT_ropt_high_noHaloVeto_rdf.root
mv data/processed/r3_cscdtOOT_ropt_high_rdf.root data/processed/r3_cscdtOOT_ropt_high_noHaloVeto_rdf.root
$pyroot skim_cscdt.py ropt oot high