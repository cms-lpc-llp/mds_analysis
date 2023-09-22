import pathlib
from os import uname

import ROOT as rt

OUT_DIR = "reports/weekly/2023-09-14"
T2_OUT_DIR = "/storage/af/user/psimmerl/LLP/mdc_analysis"  # os.getcwd()
LOCAL_OUT_DIR = "/home/psimmerl/LLP/mdc_analysis"  # os.getcwd()

DATA_VERSION = "6"

T2_DATA_DIR = "/storage/cms/store/user/christiw/displacedJetMuonAnalyzer/Run3/V1p19"
LOCAL_DATA_DIR = "/home/psimmerl/LLP/mdc_analysis/data/raw"  # os.getcwd() + "/data/raw"
DATA_DIR = "TIER2" if "caltech" in uname()[1] else "LOCAL"

FN_MC = "ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted"
FN_R3 = "DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi"

# **************************** #
if "TIER2" in DATA_DIR:
    OUT_DIR = f"{T2_OUT_DIR}/{OUT_DIR}"
    FN_MC = f"{T2_DATA_DIR}/MC_Summer22EE/v1/sixie/v{DATA_VERSION}/normalized/{FN_R3}.root"
    FN_R3 = f"{T2_DATA_DIR}/Data2022/v{DATA_VERSION}/normalized/{FN_R3}.root"
else:
    OUT_DIR = f"{LOCAL_OUT_DIR}/{OUT_DIR}"
    FN_MC = f"{LOCAL_DATA_DIR}/{FN_MC}_v{DATA_VERSION}.root"
    FN_R3 = f"{LOCAL_DATA_DIR}/{FN_R3}_v{DATA_VERSION}.root"
# pathlib.Path(OUT_DIR).mkdir(parents=True, exist_ok=True)
# **************************** #

# rt.EnableImplicitMT()

rdfn = "r3"
rdf = rt.RDataFrame("MuonSystem", FN_MC if rdfn == "mc" else FN_R3)
columns = [n for n in rdf.GetColumnNames()]
columns.remove("HLTDecision")
rdf = rdf.Filter("HLTDecision[569] && (nCscRechitClusters > 0 || nDtRechitClusters > 0)")

rdf.Snapshot("MuonSystem_HLT569", f"data/processed/{rdfn}_hlt569.rdf", columns)