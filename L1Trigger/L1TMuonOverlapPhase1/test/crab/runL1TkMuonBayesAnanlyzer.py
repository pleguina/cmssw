# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms
process = cms.Process("L1TMuonEmulation")
import os
import sys
import commands

process.load("FWCore.MessageLogger.MessageLogger_cfi")

verbose = True

if verbose: 
    process.MessageLogger = cms.Service("MessageLogger",
       #suppressInfo       = cms.untracked.vstring('AfterSource', 'PostModule'),
       destinations   = cms.untracked.vstring(
                                               #'detailedInfo',
                                               #'critical',
                                               'cout',
                                               #'cerr',
                                               'l1TkMuonBayesLog'
                    ),
       categories        = cms.untracked.vstring('l1tOmtfEventPrint'),
       l1TkMuonBayesLog = cms.untracked.PSet(    
                         extension = cms.untracked.string('.txt'),                
                         threshold = cms.untracked.string('INFO'),
                         default = cms.untracked.PSet( limit = cms.untracked.int32(0) ), 
                         #INFO   =  cms.untracked.int32(0),
                         #DEBUG   = cms.untracked.int32(0),
                         l1tOmtfEventPrint = cms.untracked.PSet( limit = cms.untracked.int32(100000000) )
                       ),
              cout = cms.untracked.PSet(    
                         threshold = cms.untracked.string('INFO'),
                         default = cms.untracked.PSet( limit = cms.untracked.int32(0) ), 
                         #INFO   =  cms.untracked.int32(0),
                         #DEBUG   = cms.untracked.int32(0),
                         l1tOmtfEventPrint = cms.untracked.PSet( limit = cms.untracked.int32(10000000) ),
                         #OMTFReconstruction = cms.untracked.PSet( limit = cms.untracked.int32(10000000) )
                         FwkReport = cms.untracked.PSet(
                             reportEvery = cms.untracked.int32(1000),
                             optionalPSet = cms.untracked.bool(True),
                             limit = cms.untracked.int32(10000000)
                             ),
                       ),
       debugModules = cms.untracked.vstring('L1TkMuonBayesTrackProducer', 'muCorrelatorAnalyzer', 'simOmtfDigis', 'muCorrelatorAnalyzer', 'simL1TkMuonBayesTrackProducer') 
       #debugModules = cms.untracked.vstring('*')
    )

    #process.MessageLogger.cout.FwkReport.reportEvery = cms.untracked.int32(100)
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
#process.load('Configuration.Geometry.GeometryExtended2023D41Reco_cff')
#process.load('Configuration.Geometry.GeometryExtended2023D41_cff')

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

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

Source_Files = cms.untracked.vstring(
 #       "/store/relval/CMSSW_10_0_0_pre1/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/94X_upgrade2023_realistic_v2_2023D17noPU-v2/10000/06C888F3-CFCE-E711-8928-0CC47A4D764C.root"
         #"/store/relval/CMSSW_9_3_2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/93X_upgrade2023_realistic_v2_2023D17noPU-v1/10000/0681719F-AFA6-E711-87C9-0CC47A4C8E14.root"
         #"file:///eos/user/k/kbunkow/cms_data/0681719F-AFA6-E711-87C9-0CC47A4C8E14.root"
         #'file:///afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonBayes/test/F4EEAE55-C937-E811-8C29-48FD8EE739D1_dump1000Events.root'
         #'file:///eos/user/k/kbunkow/cms_data/mc/PhaseIIFall17D/SingleMu_PU200_32DF01CC-A342-E811-9FE7-48D539F3863E_dump500Events.root' #not works in 10_1_7
         #'file:///eos/user/k/kbunkow/cms_data/mc/PhaseIIFall17D/ZMM_EE29AF8E-51AF-E811-A2BD-484D7E8DF0D3_dump1000Events.root'
         #"file:///eos/cms/store/group/upgrade/sandhya/SMP-PhaseIIFall17D-00001.root"
         #'file:///afs/cern.ch/work/k/kbunkow/private/omtf_data/SingleMu_15_p_1_1_qtl.root' 
         #'file:///eos/user/k/kbunkow/cms_data/mc/Phase2HLTTDRWinter20/Phase2HLTTDRWinter20DIGI__Muminus_Pt10-gun_NoPU_E6F1BC5E-BD51-A948-ADDC-8D84EFF14174_dump100Ev.root'
         #'file:///eos/user/k/kbunkow/cms_data/mc/Phase2HLTTDRWinter20/Phase2HLTTDRWinter20DIGI_DoubleMuon_gun_FlatPt-1To100_NoPU_3FD40D17-5C29-804C-B49A-029CC02B63DC_dump100Ev.root'
         'file:///eos/user/k/kbunkow/cms_data/mc/Phase2HLTTDRWinter20/Phase2HLTTDRWinter20DIGI_DoubleMuon_gun_FlatPt-1To100_NoPU_3FD40D17-5C29-804C-B49A-029CC02B63DC_dump1016Ev.root'
)


process.source = cms.Source("PoolSource", fileNames = Source_Files,
        inputCommands=cms.untracked.vstring(
        'keep *',
        'drop l1tEMTFHit2016Extras_simEmtfDigis_CSC_HLT',
        'drop l1tEMTFHit2016Extras_simEmtfDigis_RPC_HLT',
        'drop l1tEMTFHit2016s_simEmtfDigis__HLT',
        'drop l1tEMTFTrack2016Extras_simEmtfDigis__HLT',
        'drop l1tEMTFTrack2016s_simEmtfDigis__HLT',
#           'drop l1tHGCalClusterBXVector_hgcalTriggerPrimitiveDigiProducer_cluster2D_HLT',
#           'drop *HGCal*_*_*_*',
#           'drop *hgcal*_*_*_*',
#           'drop *Ecal*_*_*_*',
#           'drop *Hcal*_*_*_*',
#           'drop *Calo*_*_*_*',
#           
#           'drop *_*HGCal*_*_*',
#           'drop *_*hgcal*_*_*',
#           'drop *_*Ecal*_*_*',
#           'drop *_*Hcal*_*_*',
#           'drop *_*Calo*_*_*',
#           
#           'drop *_*_*HGCal*_*',
#           'drop *_*_*hgcal*_*',
#           'drop *_*_*Ecal*_*',
#           'drop *_*_*Hcal*_*',
#           'drop *_*_*Calo*_*'
          )
)

process.TFileService = cms.Service("TFileService", fileName = cms.string('muCorrelatorTTAnalysis1.root'), closeFileFast = cms.untracked.bool(True))


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


#######################################TTTracks################################################

####L1TkMuonBayes
process.load('L1Trigger.L1TkMuonBayes.simL1TkMuonBayesTrackProducer_cfi')
process.simL1TkMuonBayesTrackProducer.usePhase2DTPrimitives = cms.bool(False)

#process.TFileService = cms.Service("TFileService", fileName = cms.string('muCorrelatorHists.root'), closeFileFast = cms.untracked.bool(True))

process.simL1TkMuonBayesTrackProducer.L1TrackInputTag = cms.InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks")  
process.simL1TkMuonBayesTrackProducer.ttTracksSource = cms.string("L1_TRACKER")
#process.simL1TkMuonBayesTrackProducer.ttTracksSource = cms.string("SIM_TRACKS") #TODO !!!!!!!

#process.simL1TkMuonBayesTrackProducer.TrackingParticleInputTag= cms.InputTag("mix", "MergedTrackTruth") #

process.simL1TkMuonBayesTrackProducer.lctCentralBx = cms.int32(8)#<process.simL1TkMuonBayesTrackProducer.lctCentralBx = cms.int32(8)#<

process.simL1TkMuonBayesTrackProducer.pdfModuleType = cms.string("PdfModuleWithStats") #TODO
#process.simL1TkMuonBayesTrackProducer.pdfModuleFile = cms.FileInPath("L1Trigger/L1TMuonBayes/test/pdfModule.xml") #TODO!!!!!!!!!!!!!!!!!!!!!!!!!!11
#process.simL1TkMuonBayesTrackProducer.pdfModuleFile = cms.FileInPath("L1Trigger/L1TMuonBayes/test/pdfModuleSimTracks100FilesWithiRPC.xml")
#process.simL1TkMuonBayesTrackProducer.timingModuleFile  = cms.FileInPath("L1Trigger/L1TMuonBayes/test/muTimingModule100FilesWithiRPC.xml")
#process.simL1TkMuonBayesTrackProducer.timingModuleFile  = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/muTimingModuleTest.xml")
#process.simL1TkMuonBayesTrackProducer.pdfModuleFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/pdfModuleSimTracks100FilesSigma1p3.xml")  



process.simL1TkMuonBayesTrackProducer.generateTiming = cms.bool(False)
process.simL1TkMuonBayesTrackProducer.useStubsFromAdditionalBxs = cms.int32(3)

process.dumpED = cms.EDAnalyzer("EventContentAnalyzer")
process.dumpES = cms.EDAnalyzer("PrintEventSetupContent")

process.L1TMuonSeq = cms.Sequence( #process.esProd +         
                                   process.L1TrackTrigger + process.L1HybridTracksWithAssociators#+ 
                                   + process.simL1TkMuonBayesTrackProducer 
                                   #+ process.dumpED
                                   #+ process.dumpES
)

process.L1TMuonPath = cms.Path(process.L1TMuonSeq)

# process.out = cms.OutputModule("PoolOutputModule", 
#    fileName = cms.untracked.string("l1tomtf_superprimitives1.root")
# )
#process.output_step = cms.EndPath(process.out)

############################################################

analysisType = "efficiency" # or rate
  
for a in sys.argv :
    if a == "efficiency" or a ==  "rate" or a == "withTrackPart" :
        analysisType = a
        break;
    
print "analysisType=" + analysisType

process.muCorrelatorAnalyzer= cms.EDAnalyzer("MuCorrelatorAnalyzer", 
                                 outRootFile = cms.string("muCorrelatorTTAnalysis1.root"),
                                 etaCutFrom = cms.double(0.), #OMTF eta range
                                 etaCutTo = cms.double(2.4),
                                          
                                       MyProcess = cms.int32(1),
                                       DebugMode = cms.bool(verbose),      # printout lots of debug statements
                                       SaveAllTracks = cms.bool(True),   # save *all* L1 tracks, not just truth matched to primary particle
                                       SaveStubs = cms.bool(False),      # save some info for *all* stubs
                                       LooseMatch = cms.bool(True),     # turn on to use "loose" MC truth association
                                       L1Tk_minNStub = cms.int32(4),     # L1 tracks with >= 4 stubs
                                       TP_minNStub = cms.int32(4),       # require TP to have >= X number of stubs associated with it
                                       TP_minNStubLayer = cms.int32(4),  # require TP to have stubs in >= X layers/disks
                                       TP_minPt = cms.double(2.0),       # only save TPs with pt > X GeV
                                       TP_maxEta = cms.double(2.4),      # only save TPs with |eta| < X
                                       TP_maxZ0 = cms.double(30.0),      # only save TPs with |z0| < X cm
                                       TP_maxRho = cms.double(30.0),     # for efficiency analysis, to not inlude the muons from the far decays 
                                       L1TrackInputTag = cms.InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks") ,               ## TTTrack input
                                       MCTruthTrackInputTag = cms.InputTag("TTTrackAssociatorFromPixelDigis", "Level1TTTracks"), ## MCTruth input 
                                       # other input collections
                                       L1StubInputTag = cms.InputTag("TTStubsFromPhase2TrackerDigis","StubAccepted"),
                                       MCTruthClusterInputTag = cms.InputTag("TTClusterAssociatorFromPixelDigis", "ClusterAccepted"),
                                       MCTruthStubInputTag = cms.InputTag("TTStubAssociatorFromPixelDigis", "StubAccepted"),
                                       TrackingParticleInputTag = cms.InputTag("mix", "MergedTrackTruth"),
                                       TrackingVertexInputTag = cms.InputTag("mix", "MergedTrackTruth"),
                                       
                                       muCandQualityCut = cms.int32(12),
                                       analysisType = cms.string(analysisType)
                                        )
process.muCorrelatorAnalyzerPath = cms.Path(process.muCorrelatorAnalyzer)

###########################################################

############################################################

# use this if you want to re-run the stub making
#process.schedule = cms.Schedule(process.TTClusterStub,process.TTClusterStubTruth,process.TTTracksWithTruth,process.ana)

# use this if cluster/stub associators not available
#process.TTClusterStubTruth, 
process.schedule = cms.Schedule(process.L1TMuonPath, process.muCorrelatorAnalyzerPath)

# use this to only run tracking + track associator
#process.schedule = cms.Schedule(process.TTTracksWithTruth,process.ana)

#process.schedule.extend([process.output_step])
