import FWCore.ParameterSet.Config as cms

process = cms.Process("RpcKeyToRDump")

process.load('Configuration.Geometry.GeometryExtendedRun4D121Reco_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '131X_mcRun4_realistic_v9', '')

process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1))

process.dumper = cms.EDAnalyzer(
    "RpcKeyToRDumper",
    outPath=cms.string("rpc_key_to_r.csv"),
)

process.p = cms.Path(process.dumper)
