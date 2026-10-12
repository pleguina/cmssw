import FWCore.ParameterSet.Config as cms

process = cms.Process("RpcKeyToRDump")

process.load('Configuration.Geometry.GeometryExtendedRun4D110Reco_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T35', '')

process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1))

process.dumper = cms.EDAnalyzer(
    "RpcKeyToRDumper",
    outPath=cms.string("rpc_endcap_key_to_r_v2.csv"),
)

process.p = cms.Path(process.dumper)
