"""For analyzing single files (data + MC), main analysis"""
import ROOT as rt
from src.muon_system import MuonSystem


def main(**kwargs):
    pass


if __name__ == "__main__":
    out_dir = "/home/psimmerl/Documents/CMS/LLP/reports/weekly/mar2/"
    data_dir = "/home/psimmerl/Documents/CMS/LLP/data/raw/"

    mc_db_0p4 = data_dir + "ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v4.root"
    mc_ca_0p4 = data_dir + "ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v5.root"
    mc_ca_0p5 = data_dir + "ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v7.root"
    mc_ca_0p6 = data_dir + "ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v6.root"
    mc_ca_0p8 = data_dir + "ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v8.root"
    mc_ca_1p0 = data_dir + "ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v9.root"
    run3_file = data_dir + "DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi.root"

    sig = MuonSystem(mc_ca_0p6, isMC=True)
    dat = MuonSystem(run3_file, isMC=False)
    main()
