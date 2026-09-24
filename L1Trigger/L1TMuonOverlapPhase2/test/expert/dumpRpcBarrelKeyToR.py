import FWCore.ParameterSet.Config as cms

process = cms.Process("RpcBarrelKeyToRDump")

# Matches docs/wp8_p1scale_swnew_run/firmware_events_gen_cmssw20_swnew_p1scale.py's
# own geometry/GlobalTag EXACTLY (not the endcap dumpRpcKeyToR.py's D121/
# 131X pair) -- a real cross-check against that committed golden dataset
# (see reference_manifest.yaml's WP13 barrel entry) found a small,
# reproducible eta/phi discrepancy traced directly to this geometry-tag
# mismatch, not a bug in the LUT/index logic; re-dumping with the matching
# tag removes it.
process.load('Configuration.Geometry.GeometryExtendedRun4D110Reco_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T35', '')

process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1))

process.dumper = cms.EDAnalyzer(
    "RpcBarrelKeyToRDumper",
    outPath=cms.string("rpc_barrel_key_to_r.csv"),
)

process.p = cms.Path(process.dumper)
