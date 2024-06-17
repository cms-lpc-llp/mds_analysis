import ROOT as rt

# **************************** #
# STAT = 'raw'
# LUMI = 23.02 * 1000
# FN_MC = f'/home/psimmerl/mds_analysis/data/raw/ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted_v6.root'
# FN_R3 = f'/home/psimmerl/mds_analysis/data/raw/DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi_v6.root'

ORIGINAL_WEIGHT_FROM_NTUPLES = 48.580953026775774  #! idk why the weights are different in pedro's ntuples and if this scales correctly between years

YEAR, LUMI = 2022, 23.02 * 1000 * ORIGINAL_WEIGHT_FROM_NTUPLES  # 1.1328524540090597e-06 *  # ???
# YEAR, LUMI = 2023, 27.82 * 1000 * ORIGINAL_WEIGHT_FROM_NTUPLES  # 1.1328524540090597e-06 *  # ???

STAT = f"{YEAR}" #"pedro"

FN_MC = f"/home/psimmerl/mds_analysis/data/raw/mc_{YEAR}.root"
FN_R3 = f"/home/psimmerl/mds_analysis/data/raw/data_{YEAR}.root"

# LUMI = 27.82 * 1000 # 1.1328524540090597e-06 *  # ???
# FN_MC = '/home/psimmerl/mds_analysis/data/raw/mc_2023.root'
# FN_R3 = '/home/psimmerl/mds_analysis/data/raw/data_2023.root'

# **************************** #
print("Starting rdf_hlt_builder.py")
rt.EnableImplicitMT(4)
print("  Enabled ROOT's implicit multithreading (sometimes causes a crash)")
print(f"  {YEAR=:.0f}")
print(f"  {LUMI=:.0f}")
print("")

for rdfn in ("mc", "r3"):
    print("--------------------------------")
    print(f"Loading RDF: {rdfn}")
    rdf = rt.RDataFrame("MuonSystem", FN_MC if rdfn == "mc" else FN_R3)
    count_raw, wtsum_raw = rdf.Count(), rdf.Sum("weight")

    if rdfn == "mc":
        print(f"  Reweighting to LUMI (weight raw = {wtsum_raw.GetValue():,.3f})")
        rdf = rdf.Redefine("weight", f"weight * {LUMI} / {wtsum_raw.GetValue()}")
        wtsum_raw = rdf.Sum("weight")

    print(f"  Count  Raw = {count_raw.GetValue():,.0f}")
    print(f"  Weight Raw = {wtsum_raw.GetValue():,.3f}")
    print("")

    columns = [n for n in rdf.GetColumnNames()]
    columns.remove("HLTDecision")

    ################################################################
    # CSC-CSC Category, CscCluster_Loose
    print("Filtering for CSC-CSC Category (HLTDecision = CscCluster_Loose)")
    rdf566 = rdf.Filter("HLTDecision[566] && (nCscRechitClusters > 1)")  # At least 2 CSC clusters

    fn_out = f"data/processed/{rdfn}_hlt566_{STAT}.root"
    print(f"  Output Path = `{fn_out}`")
    rdf566 = rdf566.Snapshot("MuonSystem", fn_out, columns)

    count_cut, wtsum_cut = rdf566.Count(), rdf566.Sum("weight")
    print(f"  Count  Cut = {count_cut.GetValue():,.0f}")
    print(f"  Weight Cut = {wtsum_cut.GetValue():,.3f}")
    print("")

    ################################################################
    # CSC-DT Category, HLT_L1CSCCluster_DTCluster50
    print("Filtering for CSC-DT Category (HLTDecision = HLT_L1CSCCluster_DTCluster50)")
    rdf569 = rdf.Filter(
        "HLTDecision[569] && ( (nCscRechitClusters > 0) && (nDtRechitClusters > 0) )"
    )  # At least 1 CSC and 1 DT

    fn_out = f"data/processed/{rdfn}_hlt569_{STAT}.root"
    print(f"  Output Path = `{fn_out}`")
    rdf569 = rdf569.Snapshot("MuonSystem", fn_out, columns)

    count_cut, wtsum_cut = rdf569.Count(), rdf569.Sum("weight")
    print(f"  Count  Cut = {count_cut.GetValue():,.0f}")
    print(f"  Weight Cut = {wtsum_cut.GetValue():,.3f}")
    print("")

    ################################################################
    print(f"Finished RDF: {rdfn}")
    print("")

print("--------------------------------")
print("Finished building HLT reduced RDFs")
