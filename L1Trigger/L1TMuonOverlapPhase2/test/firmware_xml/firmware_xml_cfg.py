# Golden XML for the OMTF Phase-2 firmware tests: re-runs the emulator (step 2) on a
# GEN-SIM-DIGI-L1 file and dumps TestEvents XML, plus DetailedSingleEv XML for debugEvents.
# Run through generate.sh, which records what produced the files.
#
#   cmsRun firmware_xml_cfg.py inputFiles=file:<file.root> outputName=<name> [noRPC=1] [debugEvents=1,2,...]
# -*- coding: utf-8 -*-
# WP1A: generate the blessed SW-NEW golden XML dataset (0x0210/t35/nn435
# profile) from CMSSW_20_0_0, using the same native input sample and ported
# generator already verified against the legacy 0x0209/t20/bits profile in
# firmware_events_gen_cmssw20.py. Unlike that script, this one does NOT
# override configXMLFile/extrapolFactorsFilename/patternsXMLFile/
# stubEtaEncoding/useStubQualInExtr/useEndcapStubsRInExtr -- CMSSW_20_0_0's
# own simOmtfPhase2Digis_cfi.py defaults already ARE the SW-NEW profile
# (0x0210/t35/valueP2Scale). Only the NN checkpoint file and the same
# generator/debug plumbing are added on top.
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing("analysis")
options.register("outputName", "TestEvents", VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string, "TestEvents XML name, without .xml")
options.register("noRPC", 0, VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int, "1: drop RPC from the algorithm input and the XML")
options.register("debugEvents", [], VarParsing.VarParsing.multiplicity.list,
                 VarParsing.VarParsing.varType.int, "events to dump as DetailedSingleEv XML")
options.parseArguments()

process = cms.Process("L1TMuonEmulation")

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.MessageLogger = cms.Service("MessageLogger",
   destinations   = cms.untracked.vstring(
                                           'omtfEventPrint'
                ),
   categories        = cms.untracked.vstring('l1tOmtfEventPrint', 'OMTFReconstruction'),
   omtfEventPrint = cms.untracked.PSet(
                     filename  = cms.untracked.string('log_CSV_Export_Test'),
                     extension = cms.untracked.string('.txt'),
                     threshold = cms.untracked.string('DEBUG'),
                     default = cms.untracked.PSet( limit = cms.untracked.int32(0) ),
                     l1tOmtfEventPrint = cms.untracked.PSet( limit = cms.untracked.int32(1000000000) ),
                     OMTFReconstruction = cms.untracked.PSet( limit = cms.untracked.int32(1000000000) )
                   ),
   debugModules = cms.untracked.vstring('simOmtfPhase2Digis')
)

# Same CMSSW_20_0_0-native input sample used for the legacy-profile
# verification run (see reference_manifest.yaml's 2026-09-08 WP1A entries) --
# the profile (0x0209/t20 vs 0x0210/t35/nn435) only changes how the OMTF
# algorithm interprets the same raw DT/CSC/RPC digis, not what digis exist.
process.source = cms.Source('PoolSource',
  fileNames = cms.untracked.vstring(options.inputFiles),
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryExtendedRun4D110Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T35', '')

process.load("L1Trigger.DTTriggerPhase2.CalibratedDigis_cfi")
process.CalibratedDigis.dtDigiTag = "simMuonDTDigis"
process.CalibratedDigis.scenario = 0

process.load("L1Trigger.DTTriggerPhase2.dtTriggerPhase2PrimitiveDigis_cfi")
process.dtTriggerPhase2PrimitiveDigis.debug = True
process.dtTriggerPhase2PrimitiveDigis.dump = False
process.dtTriggerPhase2PrimitiveDigis.scenario = 0

process.load('L1Trigger.L1TMuonOverlapPhase2.simOmtfPhase2Digis_cfi')

# nTestRefHits investigation (2026-09-21): real CMSSW caps reference-hit
# testing at configXMLFile's own <GlobalData nTestRefHits="..."> (see
# OMTFProcessor.cc:699-714) -- a real, config-driven THROUGHPUT ceiling
# inherited from the old 160 MHz OMTF clock (160MHz*25ns = 4 cycles/BX),
# not a physics/algorithm constant. This project's own new firmware clock
# is 360 MHz, where the same 25ns BX window is 9 cycles -- matching
# GB_PAYLOAD=9. Override to a local copy of this area's own
# hwToLogicLayer_0x0210.xml with nTestRefHits bumped 4->9, everything else
# byte-identical, to regenerate golden truth reflecting the new hardware's
# real intended per-BX capacity instead of the old clock's throughput
# limit. See reference_manifest.yaml's 2026-09-21 entries.
process.simOmtfPhase2Digis.configXMLFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/hwToLogicLayer_0x0210_ntest9.xml")

process.simOmtfPhase2Digis.dumpResultToXML = cms.bool(True)
process.simOmtfPhase2Digis.XMLDumpFileName = cms.string(options.outputName + ".xml")

# TEST (per WP8 finding, 2026-09-10): the cfi default stubEtaEncoding=
# valueP2Scale carries a #TODO in CMSSW's own source ("change to
# valueP1Scale when InputMakerPhase2 is modified") and produces eta values
# (~up to 530) far outside the firmware's established range (~73-130,
# matching the legacy 'bits' encoding) -- overriding here to test whether
# InputMakerPhase2 already handles valueP1Scale correctly despite the TODO.
process.simOmtfPhase2Digis.stubEtaEncoding = cms.string("valueP1Scale")
process.simOmtfPhase2Digis.eventCaptureDebug = cms.bool(True)
process.simOmtfPhase2Digis.xmlEventsPerFile = cms.int32(500)
process.simOmtfPhase2Digis.dumpDigisToCSV = cms.bool(False)
process.simOmtfPhase2Digis.csvOutputDir = cms.string("./csv_output")

process.simOmtfPhase2Digis.dumpDetailedDebug = cms.bool(True)
process.simOmtfPhase2Digis.debugEventNumbers = cms.vint32(options.debugEvents)
process.simOmtfPhase2Digis.debugOutputDir = cms.string("./")

# NN regression checkpoint (nn435) -- the "near-free backport" confirmed in
# WP1A Phase C: NNRegression.cc:460-483 already exports regressionNN to XML,
# gated only on dumpResultToXML/getDumpResultToXML(), already set above.
process.simOmtfPhase2Digis.neuralNetworkFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/lutNN_omtfRegression_v435_FP.xml")

process.simOmtfPhase2Digis.lctCentralBx = cms.int32(8)

process.simOmtfPhase2Digis.dropRPCPrimitives = cms.bool(options.noRPC == 1)
process.simOmtfPhase2Digis.dropCSCPrimitives = cms.bool(False)

# Firmware export mode (RTL-faithful fixed-point values), same rationale as
# the legacy-profile run -- orthogonal to which pattern/NN profile is used.
process.simOmtfPhase2Digis.dtFixedPointEtaForFirmware = cms.bool(True)
process.simOmtfPhase2Digis.cscFixedPointPhiForFirmware = cms.bool(True)

# WP13 (2026-09-13): turned on to get real per-cluster RPC ground truth
# (<digi type="RPCCluster">, added by docs/wp13_rpc_cluster_export's CMSSW
# patch -- must be applied to MuonStubMakerBase.cc first) for rpc_interface
# verification. Purely additive to <inputDigis>: confirmed by diffing the
# regenerated XML against the previously-committed files with the new
# digi lines stripped out -- byte-identical, zero effect on any existing
# <inputStubs>/<referenceHits> content or on region_filter's own real
# LogicRegion-12 result (which doesn't read <inputDigis> at all).
process.simOmtfPhase2Digis.dumpRPCDigis = cms.bool(options.noRPC == 0)

# CSC r and eta from the firmware key-wire-group tables, as the firmware computes them
process.simOmtfPhase2Digis.cscFixedPointRForFirmware = cms.bool(True)
process.simOmtfPhase2Digis.cscFixedPointEtaForFirmware = cms.bool(True)

process.L1TMuonSeq = cms.Sequence(process.simOmtfPhase2Digis)

process.L1TMuonPath = cms.Path(process.L1TMuonSeq)
