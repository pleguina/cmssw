# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms
process = cms.Process("L1TMuonEmulation")
import os
import sys
import commands
import re
from os import listdir
from os.path import isfile, join

process.load("FWCore.MessageLogger.MessageLogger_cfi")

verbose = True

print("verbose", verbose) 

if verbose: 
    process.MessageLogger = cms.Service("MessageLogger",
       #suppressInfo       = cms.untracked.vstring('AfterSource', 'PostModule'),
       destinations   = cms.untracked.vstring(
                                               #'detailedInfo',
                                               #'critical',
                                               #'cout',
                                               #'cerr',
                                               'omtf2EventPrintLutGen'
                    ),
       categories        = cms.untracked.vstring('l1tMuBayesEventPrint'),
       omtf2EventPrintLutGen = cms.untracked.PSet(    
                         extension = cms.untracked.string('.txt'),                
                         threshold = cms.untracked.string('INFO'), #DEBUG INFO
                         default = cms.untracked.PSet( limit = cms.untracked.int32(0) ), 
                         #INFO   =  cms.untracked.int32(0),
                         #DEBUG   = cms.untracked.int32(0),
                         l1tMuBayesEventPrint = cms.untracked.PSet( limit = cms.untracked.int32(100000) )
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
        
    
#path = '/eos/user/a/akalinow/Data/SingleMu/9_3_14_FullEta_v2/'
#path = '/eos/user/a/akalinow/Data/9_3_14_Displaced_v4/'


path = '/eos/user/a/akalinow/Data/SingleMu/9_3_14_FullEta_v2/'

filesNameLike = sys.argv[2]

sampleType = "muons" #"muons" displacedMuons
genMode = "calcualteStdDevs" #calcualteStdDevs calcualteMeans

print("sampleType " + sampleType + ", genMode " + genMode) 

chosenFiles = []

filesPerPtBin = 100

if (sampleType == "muons") : 
    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    
    if filesNameLike == 'allPt' :
        for ptCode in range(31, 3, -1) :
            for sign in ['_m', '_p'] : #, m
                selFilesPerPtBin = 0
                for i in range(1, 101, 1): #TODO
                    for f in onlyfiles:
                        #if (( '_' + str(ptCode) + sign + '_' + str(i) + '_') in f): 
                        if (( '_' + str(ptCode) + sign + '_' + str(i) + ".") in f): 
                            #print f
                            chosenFiles.append('file://' + path + f) 
                            selFilesPerPtBin += 1
                    if(selFilesPerPtBin >= filesPerPtBin):
                        break
                            
    else :
        for i in range(1, filesPerPtBin +1, 1):
            for f in onlyfiles:
                #if (( filesNameLike + '_' + str(i) + '_') in f): 
                if (( filesNameLike + '_' + str(i) + '.') in f): 
                    print f
                    chosenFiles.append('file://' + path + f) 
                    
    fileNames = chosenFiles

# chosenFiles = []
# filesCnt = 100
# path = '/eos/user/a/akalinow/Data/9_3_14_Displaced_v3/'
# for i in range(1, filesCnt+1, 1):
#     chosenFiles.append('file://' + path + 'SingleMuFlatPt_50GeVto100GeV_cfi_py_GEN_SIM_DIGI_L1_L1TrackTrigger_DIGI2RAW_HLT_' + str(i) + '.root')
  
# path = '/eos/user/a/akalinow/Data/9_3_14_Displaced_v4/'
# for i in range(1, filesCnt+1, 1):
#     chosenFiles.append('file://' + path + 'SingleMuFlatPt_50GeVto100GeV_cfi_py_GEN_SIM_DIGI_L1_L1TrackTrigger_DIGI2RAW_HLT_' + str(i) + '.root')

if (sampleType == "displacedMuons") : 
    chosenFiles = []
    #chosenFiles.append('file:///afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_0_pre4/src/L1Trigger/L1TMuonBayes/test/expert/DisplacedMuonGun_Pt30To100_Dxy_0_1000_E68C6334-7F62-E911-8AA5-0025905B8610_dump2000Ev.root')
     
    path = '/eos/user/k/kbunkow/cms_data/RelValDisplacedMuonGun/'
    chosenFiles.append('file://' + path + 'DisplacedMuonGun_Pt30To100_Dxy_0_1000_C2D5C228-7F62-E911-AAA6-0CC47A78A42C_dumpAllEv.root')
    chosenFiles.append('file://' + path + 'DisplacedMuonGun_Pt30To100_Dxy_0_1000_E68C6334-7F62-E911-8AA5-0025905B8610_dumpAllEv.root')
    chosenFiles.append('file://' + path + 'DisplacedMuonGun_Pt30To100_Dxy_0_1000_E68C6334-7F62-E911-8AA5-0025905B8610_dumpAllEv.root')
    chosenFiles.append('file://' + path + 'DisplacedMuonGun_Pt30To100_Dxy_0_1000_3C065D44-7E62-E911-9F3F-0CC47A4C8F12_dumpAllEv.root')
    chosenFiles.append('file://' + path + 'DisplacedMuonGun_Pt30To100_Dxy_0_1000_34016F34-7F62-E911-AAB8-0025905AA9CC_dumpAllEv.root')

print "chosenFiles"
for chFile in chosenFiles:
    print chFile

if len(chosenFiles) == 0 :
    print "no files selected!!!!!!!!!!!!!!!"
    exit

process.source = cms.Source('PoolSource',
fileNames = cms.untracked.vstring(list(chosenFiles) ),
skipEvents =  cms.untracked.uint32(0),

        inputCommands=cms.untracked.vstring(
        'keep *',
        'drop l1tEMTFHit2016Extras_simEmtfDigis_CSC_HLT',
        'drop l1tEMTFHit2016Extras_simEmtfDigis_RPC_HLT',
        'drop l1tEMTFHit2016s_simEmtfDigis__HLT',
        'drop l1tEMTFTrack2016Extras_simEmtfDigis__HLT',
        'drop l1tEMTFTrack2016s_simEmtfDigis__HLT')
)



                        
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))

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

process.TFileService = cms.Service("TFileService", fileName = cms.string('omtf2LutGen_' + sampleType + '_' + genMode + '.root'), closeFileFast = cms.untracked.bool(True))

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

process.simOmt2fDigis.dumpResultToXML = cms.bool(False)

process.simOmt2fDigis.ptModuleType = cms.string("PtModuleLut2DGen")

if (genMode == "calcualteStdDevs"):
    process.simOmt2fDigis.ptModuleFile = cms.FileInPath("L1Trigger/L1TMuonBayes/test/expert/omtf2/ptModuleLut2D_GenMean_" + sampleType +  ".xml")

process.simOmt2fDigis.genMode = cms.string(genMode)
process.simOmt2fDigis.sampleType = cms.string(sampleType)

process.simOmt2fDigis.etaCutFrom = cms.double(0.8)
process.simOmt2fDigis.etaCutTo = cms.double(1.24)


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
