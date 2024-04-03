#!/bin/bash

pyroot=/home/psimmerl/.miniconda3/envs/pyroot/bin/python
# pyroot=/home/psimmerl/mambaforge/envs/pyroot/bin/python

###########
# CSC-CSC #
###########
# Reduced cut set
ll="_l1"
time $pyroot skim_csccsc.py l1 it
mv data/processed/mc_csccsc_rdf.root data/processed/mc_csccsc${ll}_rdf.root
mv data/processed/r3_csccsc_rdf.root data/processed/r3_csccsc${ll}_rdf.root

ll="OOT_l1"
time $pyroot skim_csccsc.py l1 oot
mv data/processed/mc_csccsc_rdf.root data/processed/mc_csccsc${ll}_rdf.root
mv data/processed/r3_csccsc_rdf.root data/processed/r3_csccsc${ll}_rdf.root

# # Standard selections until I LOO
# ll="_scs"
# time $pyroot skim_csccsc.py it low
# mv data/processed/mc_csccsc_rdf.root data/processed/mc_csccsc${ll}_low_rdf.root
# mv data/processed/r3_csccsc_rdf.root data/processed/r3_csccsc${ll}_low_rdf.root

# time $pyroot skim_csccsc.py it high
# mv data/processed/mc_csccsc_rdf.root data/processed/mc_csccsc${ll}_high_rdf.root
# mv data/processed/r3_csccsc_rdf.root data/processed/r3_csccsc${ll}_high_rdf.root

# ll="OOT_scs"
# time $pyroot skim_csccsc.py oot low
# mv data/processed/mc_csccsc_rdf.root data/processed/mc_csccsc${ll}_low_rdf.root
# mv data/processed/r3_csccsc_rdf.root data/processed/r3_csccsc${ll}_low_rdf.root

# time $pyroot skim_csccsc.py oot high
# mv data/processed/mc_csccsc_rdf.root data/processed/mc_csccsc${ll}_high_rdf.root
# mv data/processed/r3_csccsc_rdf.root data/processed/r3_csccsc${ll}_high_rdf.root


# ##########
# # CSC-DT #
# ##########
# ll="_loo"
# time $pyroot skim_cscdt.py it low
# mv data/processed/mc_cscdt_rdf.root data/processed/mc_cscdt${ll}_low_rdf.root
# mv data/processed/r3_cscdt_rdf.root data/processed/r3_cscdt${ll}_low_rdf.root

# time $pyroot skim_cscdt.py it high
# mv data/processed/mc_cscdt_rdf.root data/processed/mc_cscdt${ll}_high_rdf.root
# mv data/processed/r3_cscdt_rdf.root data/processed/r3_cscdt${ll}_high_rdf.root

# ll="OOT_loo"
# time $pyroot skim_cscdt.py oot low
# mv data/processed/mc_cscdt_rdf.root data/processed/mc_cscdt${ll}_low_rdf.root
# mv data/processed/r3_cscdt_rdf.root data/processed/r3_cscdt${ll}_low_rdf.root

# time $pyroot skim_cscdt.py oot high
# mv data/processed/mc_cscdt_rdf.root data/processed/mc_cscdt${ll}_high_rdf.root
# mv data/processed/r3_cscdt_rdf.root data/processed/r3_cscdt${ll}_high_rdf.root

