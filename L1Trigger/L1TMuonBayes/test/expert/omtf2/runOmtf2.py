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
                                               #'cout',
                                               #'cerr',
                                               'omtf2EventPrint'
                    ),
       categories        = cms.untracked.vstring('l1tMuBayesEventPrint'),
       omtf2EventPrint = cms.untracked.PSet(    
                         extension = cms.untracked.string('.txt'),                
                         threshold = cms.untracked.string('DEBUG'),
                         default = cms.untracked.PSet( limit = cms.untracked.int32(0) ), 
                         #INFO   =  cms.untracked.int32(0),
                         #DEBUG   = cms.untracked.int32(0),
                         l1tMuBayesEventPrint = cms.untracked.PSet( limit = cms.untracked.int32(10000000) )
                       ),
       debugModules = cms.untracked.vstring('L1TMuonBayesMuCorrelatorTrackProducer', 'OmtfTTAnalyzer', 'simOmt2fDigis', 'omtfTTAnalyzer', 'simBayesMuCorrelatorTrackProducer') 
       #debugModules = cms.untracked.vstring('*')
    )

    #process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)
if not verbose:
    process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
    process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False), 
                                         #SkipEvent = cms.untracked.vstring('ProductNotFound') 
                                        )
    
    
    
path = '/eos/user/a/akalinow/Data/SingleMu/9_3_14_FullEta_v2/'
#path = '/eos/user/a/akalinow/Data/9_3_14_Displaced_v4/'

process.source = cms.Source('PoolSource',
 #fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/g/gflouris/public/SingleMuPt6180_noanti_10k_eta1.root')
#fileNames = cms.untracked.vstring('file:///eos/user/a/akalinow/Data/SingleMu/9_3_14_FullEta_v2/SingleMu_10_m_1.root'),      
#fileNames = cms.untracked.vstring('file:///eos/user/a/akalinow/Data/9_3_14_Displaced_v3/SingleMuFlatPt_50GeVto100GeV_cfi_py_GEN_SIM_DIGI_L1_L1TrackTrigger_DIGI2RAW_HLT_500.root'),                    
fileNames = cms.untracked.vstring(#'file://' + path + 'SingleMu_18_p_1.root', 
                                  'file://' + path + 'SingleMu_10_p_1.root',
                                  #'file://' + '/eos/user/k/kbunkow/cms_data/RelValDisplacedMuonGun/' + 'DisplacedMuonGun_Pt30To100_Dxy_0_1000_C2D5C228-7F62-E911-AAA6-0CC47A78A42C_dumpAllEv.root'
                                  ),                           


skipEvents =  cms.untracked.uint32(0),
        inputCommands=cms.untracked.vstring(
        'keep *',
        'drop l1tEMTFHit2016Extras_simEmtfDigis_CSC_HLT',
        'drop l1tEMTFHit2016Extras_simEmtfDigis_RPC_HLT',
        'drop l1tEMTFHit2016s_simEmtfDigis__HLT',
        'drop l1tEMTFTrack2016Extras_simEmtfDigis__HLT',
        'drop l1tEMTFTrack2016s_simEmtfDigis__HLT')
)
	               

                        
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500))

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
#process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D17_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '103X_upgrade2023_realistic_v2', '') 

process.TFileService = cms.Service("TFileService", fileName = cms.string('omtf2.root'), closeFileFast = cms.untracked.bool(True))

####Event Setup Producer
process.load('L1Trigger.L1TMuonBayes.fakeOmtfParams_cff')
process.esProd = cms.EDAnalyzer("EventSetupRecordDataGetter",
   toGet = cms.VPSet(
      cms.PSet(record = cms.string('L1TMuonOverlapParamsRcd'),
               data = cms.vstring('L1TMuonOverlapParams'))
                   ),
   verbose = cms.untracked.bool(False)
)


####OMTF Emulator
process.load('L1Trigger.L1TMuonBayes.simOmtf2Digis_cfi')

process.simOmt2fDigis.dumpResultToXML = cms.bool(True)

process.simOmt2fDigis.ptModuleType = cms.string("PtModuleLut2D")
process.simOmt2fDigis.ptModuleFile = cms.FileInPath("L1Trigger/L1TMuonBayes/test/expert/omtf2/ptModuleLut2D_GenStdDevs_muons.xml")
#process.simOmt2fDigis.ptModuleFile = cms.FileInPath("L1Trigger/L1TMuonBayes/test/expert/omtf2/ptModuleLut2D_GenStdDevs_displacedMuons.xml")

process.dumpED = cms.EDAnalyzer("EventContentAnalyzer")
process.dumpES = cms.EDAnalyzer("PrintEventSetupContent")

process.L1TMuonSeq = cms.Sequence( process.esProd          
                                   + process.simOmt2fDigis 
                                   #+ process.dumpED
                                   #+ process.dumpES
)

process.L1TMuonPath = cms.Path(process.L1TMuonSeq)

# process.out = cms.OutputModule("PoolOutputModule", 
#    fileName = cms.untracked.string("l1tomtf_superprimitives1.root")
# )

#process.output_step = cms.EndPath(process.out)
process.schedule = cms.Schedule(process.L1TMuonPath)
#process.schedule.extend([process.output_step])
