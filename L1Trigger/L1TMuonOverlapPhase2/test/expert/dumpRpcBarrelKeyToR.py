import FWCore.ParameterSet.Config as cms

process = cms.Process("RpcBarrelKeyToRDump")

# Must match firmware_events_gen_cmssw20_swnew_p1scale.py's geometry/
# GlobalTag (not endcap dumpRpcKeyToR.py's D121/131X pair) -- a mismatch
# here was previously traced to a small eta/phi discrepancy, not a LUT bug.
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
