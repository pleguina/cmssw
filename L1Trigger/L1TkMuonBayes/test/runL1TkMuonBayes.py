# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms
process = cms.Process("L1TMuonEmulation")
import os
import sys
import commands

from Configuration.Eras.Era_Phase2C9_cff import Phase2C9

process.load("FWCore.MessageLogger.MessageLogger_cfi")

verbose = True

if verbose: 
    process.MessageLogger = cms.Service("MessageLogger",
       #suppressInfo       = cms.untracked.vstring('AfterSource', 'PostModule'),
       destinations   = cms.untracked.vstring(
                                               #'detailedInfo',
                                               #'critical',
                                               #'cout',
                                               #'cerr',
                                               'muCorrelatorEventPrint'
                    ),
       categories        = cms.untracked.vstring('l1tOmtfEventPrint'),
       muCorrelatorEventPrint = cms.untracked.PSet(    
                         extension = cms.untracked.string('.txt'),                
                         threshold = cms.untracked.string('DEBUG'),
                         default = cms.untracked.PSet( limit = cms.untracked.int32(0) ), 
                         #INFO   =  cms.untracked.int32(0),
                         #DEBUG   = cms.untracked.int32(0),
                         l1tOmtfEventPrint = cms.untracked.PSet( limit = cms.untracked.int32(100000000) )
                       ),
       debugModules = cms.untracked.vstring('L1TkMuonBayesTrackProducer', 'simOmtfDigis', 'simL1TkMuonBayesTrackProducer') 
       #debugModules = cms.untracked.vstring('*')
    )

    #process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)
if not verbose:
    process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)
    process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False), 
                                         #SkipEvent = cms.untracked.vstring('ProductNotFound') 
                                     )


#######################################TTTracks################################################

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
#process.load('Configuration.Geometry.GeometryExtended2026D41Reco_cff')
#process.load('Configuration.Geometry.GeometryExtended2026D41_cff')

process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff') #TODO!!!!!!!!!!!!! this geometry is required for Phase2HLTTDRWinter20

process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '111X_mcRun4_realistic_T15_v3', '')


############################################################
# input and output
############################################################

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(20))

Source_Files = cms.untracked.vstring(
          #'file:///eos/user/k/kbunkow/cms_data/mc/Phase2HLTTDRWinter20/Phase2HLTTDRWinter20DIGI__Muminus_Pt10-gun_NoPU_E6F1BC5E-BD51-A948-ADDC-8D84EFF14174_dump100Ev.root'
          'file:///eos/user/k/kbunkow/cms_data/mc/Phase2HLTTDRWinter20/Phase2HLTTDRWinter20DIGI_DoubleMuon_gun_FlatPt-1To100_NoPU_3FD40D17-5C29-804C-B49A-029CC02B63DC_dump1016Ev.root'
)

process.source = cms.Source("PoolSource", fileNames = Source_Files,
        inputCommands=cms.untracked.vstring(
        'keep *',
        'drop l1tEMTFHit2016Extras_simEmtfDigis_CSC_HLT',
        'drop l1tEMTFHit2016Extras_simEmtfDigis_RPC_HLT',
        'drop l1tEMTFHit2016s_simEmtfDigis__HLT',
        'drop l1tEMTFTrack2016Extras_simEmtfDigis__HLT',
        'drop l1tEMTFTrack2016s_simEmtfDigis__HLT')
)

#process.TFileService = cms.Service("TFileService", fileName = cms.string('muCorrelatorTTAnalysis1.root'), closeFileFast = cms.untracked.bool(True))


############################################################
# remake L1 stubs and/or cluster/stub truth ??
############################################################

process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
from L1Trigger.TrackTrigger.TTStubAlgorithmRegister_cfi import *
process.load("SimTracker.TrackTriggerAssociation.TrackTriggerAssociator_cff")

from SimTracker.TrackTriggerAssociation.TrackTriggerAssociator_cff import *
TTClusterAssociatorFromPixelDigis.digiSimLinks = cms.InputTag("simSiPixelDigis","Tracker")

process.TTClusterStub = cms.Path(process.TrackTriggerClustersStubs)
process.TTClusterStubTruth = cms.Path(process.TrackTriggerAssociatorClustersStubs)


from L1Trigger.TrackFindingTracklet.L1HybridEmulationTracks_cff import *


############################################################
# L1 tracking
############################################################

#from L1Trigger.TrackFindingTracklet.Tracklet_cfi import *
#if GEOMETRY == "D10": 
#    TTTracksFromTracklet.trackerGeometry = cms.untracked.string("flat")
#TTTracksFromTracklet.asciiFileName = cms.untracked.string("evlist.txt")

process.L1TrackTrigger_step = cms.Path(process.L1TrackTrigger)

#process.TTTracksWithTruth = cms.Path(process.L1TrackletTracksWithAssociators)


#######################################TTTracks################################################

##This overrides the tracker geometry and the TTTriger does not work!!!!!!!!!!!!
# # PostLS1 geometry used
# process.load('Configuration.Geometry.GeometryExtended2015Reco_cff')
# process.load('Configuration.Geometry.GeometryExtended2015_cff')
# ############################
# process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
# from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
# process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')


####Event Setup Producer
# process.load('L1Trigger.L1TMuonOverlap.fakeOmtfParams_cff')
# process.esProd = cms.EDAnalyzer("EventSetupRecordDataGetter",
#    toGet = cms.VPSet(
#       cms.PSet(record = cms.string('L1TMuonOverlapParamsRcd'),
#                data = cms.vstriL1TMuonBayesrlapParams'))
#                    ),
#    verbose = cms.untracked.bool(False)
# )


####OMTF Emulator
process.load('L1Trigger.L1TkMuonBayes.simL1TkMuonBayesTrackProducer_cfi')
process.simL1TkMuonBayesTrackProducer.usePhase2DTPrimitives = cms.bool(False)

process.simL1TkMuonBayesTrackProducer.L1TrackInputTag = cms.InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks")  
#process.simL1TkMuonBayesTrackProducer.ttTracksSource = cms.string("TRACKING_PARTICLES") #
process.simL1TkMuonBayesTrackProducer.ttTracksSource = cms.string("L1_TRACKER") #

process.simL1TkMuonBayesTrackProducer.TrackingParticleInputTag= cms.InputTag("mix", "MergedTrackTruth") #trackingParticleTag

process.simL1TkMuonBayesTrackProducer.lctCentralBx = cms.int32(8)#<<<<<<<<<<<<<<<<!!!!!!!!!!!!!!!!!!!!TODO this was changed in CMSSW 10(?) to 8. if the data were generated with the previous CMSSW then you have to use 6

#process.simL1TkMuonBayesTrackProducer.pdfModuleFile = cms.FileInPath("L1Trigger/L1TMuonOverlapPhase1/test/pdfModule.xml")

process.dumpED = cms.EDAnalyzer("EventContentAnalyzer")
process.dumpES = cms.EDAnalyzer("PrintEventSetupContent")

process.L1TMuonSeq = cms.Sequence( #process.esProd +    
                                   process.L1TrackTrigger + L1HybridTracksWithAssociators#+     
                                   + process.simL1TkMuonBayesTrackProducer 
                                   #+ process.dumpED
                                   #+ process.dumpES
)



process.L1TMuonPath = cms.Path(process.L1TMuonSeq)

process.out = cms.OutputModule("PoolOutputModule", 
    fileName = cms.untracked.string("outCollections.root"),
    outputCommands=cms.untracked.vstring(
        #'drop *',
        'keep *',
         'drop l1tEMTFHit2016Extras_simEmtfDigis_CSC_HLT',
        'drop l1tEMTFHit2016Extras_simEmtfDigis_RPC_HLT',
        'drop l1tEMTFHit2016s_simEmtfDigis__HLT',
        'drop l1tEMTFTrack2016Extras_simEmtfDigis__HLT',
        'drop l1tEMTFTrack2016s_simEmtfDigis__HLT',
        'drop *HGCal*_*_*_*'
        )
)
process.output_step = cms.EndPath(process.out)


############################################################

# use this if you want to re-run the stub making
#process.schedule = cms.Schedule(process.TTClusterStub,process.TTClusterStubTruth,process.TTTracksWithTruth,process.ana)

# use this if cluster/stub associators not available 
#process.schedule = cms.Schedule(process.TTClusterStubTruth, process.TTTracksWithTruth, process.L1TMuonPath)
#process.schedule = cms.Schedule(process.TTTracksWithTruth, process.L1TMuonPath)
#process.schedule = cms.Schedule(process.L1TMuonPath)
process.schedule = cms.Schedule(#process.TTClusterStub, 
                                #process.TTClusterStubTruth, 
                                process.L1TMuonPath)

# use this to only run tracking + track associator
#process.schedule = cms.Schedule(process.TTTracksWithTruth,process.ana)

process.schedule.extend([process.output_step])
