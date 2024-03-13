#!/bin/bash

pyroot=/home/psimmerl/.miniconda3/envs/pyroot/bin/python

# CSC-CSC
ll="_pedro"
time $pyroot skim_csccsc.py it
mv data/processed/mc_csccsc_rdf.root data/processed/mc_csccsc${ll}_rdf.root
mv data/processed/r3_csccsc_rdf.root data/processed/r3_csccsc${ll}_rdf.root

ll="_pedro_cscOOT"
time $pyroot skim_csccsc.py oot
mv data/processed/mc_csccsc_rdf.root data/processed/mc_csccsc${ll}_rdf.root
mv data/processed/r3_csccsc_rdf.root data/processed/r3_csccsc${ll}_rdf.root

# CSC-DT
ll="_loo"
time $pyroot skim_cscdt.py it low
mv data/processed/mc_cscdt_rdf.root data/processed/mc_cscdt${ll}_low_rdf.root
mv data/processed/r3_cscdt_rdf.root data/processed/r3_cscdt${ll}_low_rdf.root

time $pyroot skim_cscdt.py it high
mv data/processed/mc_cscdt_rdf.root data/processed/mc_cscdt${ll}_high_rdf.root
mv data/processed/r3_cscdt_rdf.root data/processed/r3_cscdt${ll}_high_rdf.root

ll="_loo_dtOOT"
time $pyroot skim_cscdt.py oot low
mv data/processed/mc_cscdt_rdf.root data/processed/mc_cscdt${ll}_low_rdf.root
mv data/processed/r3_cscdt_rdf.root data/processed/r3_cscdt${ll}_low_rdf.root

time $pyroot skim_cscdt.py oot high
mv data/processed/mc_cscdt_rdf.root data/processed/mc_cscdt${ll}_high_rdf.root
mv data/processed/r3_cscdt_rdf.root data/processed/r3_cscdt${ll}_high_rdf.root

