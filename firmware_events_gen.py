# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms
process = cms.Process("L1TMuonEmulation")
import os
import sys

process.load("FWCore.MessageLogger.MessageLogger_cfi")

# Configure logging for debugging
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

# Use the same input file as runMuonOverlap.py (with file: prefix for local access)
process.source = cms.Source('PoolSource',  
  fileNames = cms.untracked.vstring('file:/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/OMTF/PrivateProductionForOMTFStudy/13_1_0_03_04_2024/SingleMu_ch0_OneOverPt_Run2029_13_1_0_03_04_2024/13_1_0_03_04_2024/240403_080928/0000/SingleMu_OneOverPt_1_100_m_1.root')                  
)
	                    
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100))

# import of standard configurations - same as runMuonOverlap.py
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryExtended2026D95Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '') 

# Calibrate Digis - same as runMuonOverlap.py
process.load("L1Trigger.DTTriggerPhase2.CalibratedDigis_cfi")
process.CalibratedDigis.dtDigiTag = "simMuonDTDigis" 
process.CalibratedDigis.scenario = 0

# DTTriggerPhase2 - same as runMuonOverlap.py
process.load("L1Trigger.DTTriggerPhase2.dtTriggerPhase2PrimitiveDigis_cfi")
process.dtTriggerPhase2PrimitiveDigis.debug = True
process.dtTriggerPhase2PrimitiveDigis.dump = False
process.dtTriggerPhase2PrimitiveDigis.scenario = 0

# OMTF Emulator with CSV export enabled
process.load('L1Trigger.L1TMuonOverlapPhase2.simOmtfPhase2Digis_cfi')

process.simOmtfPhase2Digis.dumpResultToXML = cms.bool(True)
process.simOmtfPhase2Digis.eventCaptureDebug = cms.bool(True)
# === ENABLE HLS CSV EXPORT ===
process.simOmtfPhase2Digis.dumpDigisToCSV = cms.bool(False)
process.simOmtfPhase2Digis.csvOutputDir = cms.string("./csv_output")

# === ENABLE DETAILED DEBUG EXPORT (for one specific event) ===
# This generates a SEPARATE detailed XML for ONLY the specified event ID
# debugEventNumber = the actual EVENT ID from the data file
# To find event IDs, check TestEvents.xml after running once
process.simOmtfPhase2Digis.dumpDetailedDebug = cms.bool(True)
process.simOmtfPhase2Digis.debugEventNumber = cms.int32(55)  # Exact event ID to debug
process.simOmtfPhase2Digis.debugOutputDir = cms.string("./")

# === USE LOOKUP TABLE FOR ETA CALCULATION ===
process.simOmtfPhase2Digis.stubEtaEncoding = cms.string("bits")

# === ENABLE STUB QUALITY IN EXTRAPOLATION ===
# This must be True to match the ExtrapolationFactors XML which uses KeyType="quality"
process.simOmtfPhase2Digis.useStubQualInExtr = cms.bool(True)

# === USE ETA (NOT R) FOR ENDCAP STUB EXTRAPOLATION ===
# Set to False to use eta as the key for CSC/RPC endcap layers (matches XML KeyType="eta")
process.simOmtfPhase2Digis.useEndcapStubsRInExtr = cms.bool(False)

# Keep your specific XML configuration files
process.simOmtfPhase2Digis.extrapolFactorsFilename = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/ExtrapolationFactors_ExtraplMB1nadMB2DTQual_ValueP1Scale_t20.xml")
process.simOmtfPhase2Digis.patternsXMLFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/Patterns_ExtraplMB1nadMB2DTQualAndEtaFixedP_ValueP1Scale_t20_v1_SingleMu_iPt_and_OneOverPt_classProb17_recalib2_minDP0.xml")

# Use same configuration as default simOmtfPhase2Digis_cfi.py
process.simOmtfPhase2Digis.lctCentralBx = cms.int32(8)

# Test the Phase2 fix: enable all primitives to see hwName in all digi types
process.simOmtfPhase2Digis.dropRPCPrimitives = cms.bool(False)
process.simOmtfPhase2Digis.dropCSCPrimitives = cms.bool(False)

# === FIRMWARE ETA MODE ===
# Use fixed mid-chamber DT eta values to match RTL constants (firmware export only)
process.simOmtfPhase2Digis.dtFixedPointEtaForFirmware = cms.bool(True)

# === FIRMWARE CSC PHI MODE ===
# Use Q2.8 fixed-point CSC phi arithmetic to match RTL CSC interface (firmware export only)
process.simOmtfPhase2Digis.cscFixedPointPhiForFirmware = cms.bool(True)

# === DISABLE RPC DIGI EXPORT (keep RPC stubs) ===
process.simOmtfPhase2Digis.dumpRPCDigis = cms.bool(False)

process.L1TMuonSeq = cms.Sequence(process.simOmtfPhase2Digis)

process.L1TMuonPath = cms.Path(process.L1TMuonSeq)