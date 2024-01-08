# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms
process = cms.Process("L1TMuonEmulation")
import os
import sys
import re
from os import listdir
from os.path import isfile, join
import fnmatch


process.load("FWCore.MessageLogger.MessageLogger_cfi")

verbose = True

useExtraploationAlgo = True

version = 't22__'

if useExtraploationAlgo :
    version = version + 'Patterns_ExtraplMB1nadMB2SimplifiedFP_t17_classProb17_recalib2_minDP0_v3_gpFinalize10'
    #version = version + 'Patterns_ExtraplMB1nadMB2SimplifiedFP_t17_classProb17_recalib2_gpFinalize10'
    #version = version + 'Patterns_ExtraplMB1nadMB2FullAlgo_t16_classProb17_recalib2_gpFinalize10'
else :
    version = version + 'Patterns_0x00012'

runDebug = "INFO" # or "INFO" DEBUG
#useExtraploationAlgo = True


analysisType = "efficiency" # or rate
  
for a in sys.argv :
    if a == "efficiency" or a ==  "rate" or a == "withTrackPart" :
        analysisType = a
        break;
    
filesNameLike = sys.argv[2]
    
outFilesName = 'omtfAnalysis2_' 
if analysisType == "efficiency" :
    outFilesName = outFilesName + "eff_"
elif analysisType == "rate" :
    outFilesName = outFilesName + "rate_"    
    
outFilesName = outFilesName + version + "__" + filesNameLike

if(runDebug == "DEBUG") :
    outFilesName = outFilesName + "_test6"

if verbose: 
    process.MessageLogger = cms.Service("MessageLogger",
       #suppressInfo       = cms.untracked.vstring('AfterSource', 'PostModule'),
       destinations   = cms.untracked.vstring(
                                               #'detailedInfo',
                                               #'critical',
                                               #'cout',
                                               #'cerr',
                                               'omtfEventPrint'
                    ),
       categories        = cms.untracked.vstring( 'OMTFReconstruction', 'l1tOmtfEventPrint', 'l1MuonAnalyzerOmtf'), #'l1tOmtfEventPrint', 'l1MuonAnalyzerOmtf'
       omtfEventPrint = cms.untracked.PSet(    
                         filename  = cms.untracked.string(outFilesName),
                         extension = cms.untracked.string('.txt'),                
                         threshold = cms.untracked.string("INFO"), #DEBUG
                         default = cms.untracked.PSet( limit = cms.untracked.int32(0) ), 
                         #INFO   =  cms.untracked.int32(0),
                         #DEBUG   = cms.untracked.int32(0),
                         l1tOmtfEventPrint = cms.untracked.PSet( limit = cms.untracked.int32(1000000000) ),
                         OMTFReconstruction = cms.untracked.PSet( limit = cms.untracked.int32(1000000000) ),
                         l1MuonAnalyzerOmtf = cms.untracked.PSet( limit = cms.untracked.int32(1000000000) ),
                       ),
       debugModules = cms.untracked.vstring('L1MuonAnalyzerOmtf', 'simOmtfDigis', 'omtfParamsSource', 'omtfParams', "esProd", 'L1TMuonOverlapPhase1ParamsESProducer') #'L1MuonAnalyzerOmtf',
       #debugModules = cms.untracked.vstring('*')
    )

    #process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)
if not verbose:
    process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
    process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False), 
                                         #SkipEvent = cms.untracked.vstring('ProductNotFound') 
                                     )
    
    
# PostLS1 geometry used
process.load('Configuration.Geometry.GeometryExtended2015Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2015_cff')
############################
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
#from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')    
    
    
# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
#process.load('Configuration.Geometry.GeometryExtended2026D41Reco_cff')
#process.load('Configuration.Geometry.GeometryExtended2026D41_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
#process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '103X_upgrade2023_realistic_v2', '') 

chosenFiles = []

fileCnt = 100000 #1000 

 
       
if filesNameLike == "SingleMu_9_3_14_FullEta_v2" :    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cscBx = 6
    matchUsingPropagation  = False 
    paths = [
        '/eos/user/a/akalinow/Data/SingleMu/9_3_14_FullEta_v2/'
        ]  
    # fileCnt = 10 #<<<<<<!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if filesNameLike == 'mcWaw_2024_01_03_OneOverPt' :
    cscBx = 8
    matchUsingPropagation  = False 
    paths = [    
             {"path": "/eos/user/a/akalinow/Data/SingleMu/13_1_0_03_01_2024/SingleMu_ch0_OneOverPt_Run2029_13_1_0_03_01_2024/", "fileCnt" : 500}, #1000 files
             {"path": "/eos/user/a/akalinow/Data/SingleMu/13_1_0_03_01_2024/SingleMu_ch2_OneOverPt_Run2029_13_1_0_03_01_2024/", "fileCnt" : 500} #1000 files
             ]
    #fileCnt = 10 #<<<<<<!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if filesNameLike == 'mcWaw2023_OneOverPt_and_iPt2':
    cscBx = 8
    matchUsingPropagation  = False 
    paths = [
             # "/eos/user/a/akalinow/Data/SingleMu/12_5_2_p1_20_04_2023/SingleMu_ch0_OneOverPt_12_5_2_p1_20_04_2023/", #500 files
             # "/eos/user/a/akalinow/Data/SingleMu/12_5_2_p1_20_04_2023/SingleMu_ch2_OneOverPt_12_5_2_p1_20_04_2023/", #500 files
             #
             # "/eos/user/a/akalinow/Data/SingleMu/12_5_2_p1_14_04_2023/SingleMu_ch0_OneOverPt_12_5_2_p1_14_04_2023/", #500 files
             # "/eos/user/a/akalinow/Data/SingleMu/12_5_2_p1_14_04_2023/SingleMu_ch2_OneOverPt_12_5_2_p1_14_04_2023/", #500 files
             #
             {"path": "/eos/user/a/akalinow/Data/SingleMu/12_5_2_p1_04_04_2023/SingleMu_ch0_OneOverPt_12_5_2_p1_04_04_2023/", "fileCnt" : 300}, #500 files
             {"path": "/eos/user/a/akalinow/Data/SingleMu/12_5_2_p1_04_04_2023/SingleMu_ch2_OneOverPt_12_5_2_p1_04_04_2023/", "fileCnt" : 300}, #500 files
             #
             # "/eos/user/a/akalinow/Data/SingleMu/12_5_2_p1_22_02_2023/SingleMu_ch0_OneOverPt_12_5_2_p1_22_02_2023/", #200 files
             # "/eos/user/a/akalinow/Data/SingleMu/12_5_2_p1_22_02_2023/SingleMu_ch2_OneOverPt_12_5_2_p1_22_02_2023/", #200 files
             #
             # "/eos/user/a/akalinow/Data/SingleMu/12_5_2_p1_15_02_2023/SingleMu_ch0_OneOverPt_12_5_2_p1_15_02_2023/", ##100 files
             # "/eos/user/a/akalinow/Data/SingleMu/12_5_2_p1_15_02_2023/SingleMu_ch2_OneOverPt_12_5_2_p1_15_02_2023/", ##100 files
              {"path": "/eos/user/a/akalinow/Data/SingleMu/12_5_2_p1_04_04_2023/SingleMu_ch0_iPt2_12_5_2_p1_04_04_2023/", "fileCnt" : 10000}, #500 files
              {"path": "/eos/user/a/akalinow/Data/SingleMu/12_5_2_p1_04_04_2023/SingleMu_ch2_iPt2_12_5_2_p1_04_04_2023/", "fileCnt" : 10000}, #500 files
             ]

if filesNameLike == "EfeMC_HTo2LongLivedTo2mu2jets" :    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cscBx = 8
    matchUsingPropagation  = True 
    paths = [
        {"path": '/eos/cms/store/user/eyigitba/dispDiMu/crabOut/CRAB_PrivateMC/', "fileCnt" : 10000},
        ]   
    
if filesNameLike == "Displaced_Dxy5m_pT0To1000_condRun3" :    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cscBx = 8
    matchUsingPropagation  = False 
    paths = [
        {"path": '/eos/user/a/akalinow/Data/SingleMu/Displaced_Dxy5m_pT0To1000_condRun3_131X_mcRun3_2023_realistic_v10/DisplacedMu_ch0_iPt0_Run2023_13_1_0_23_11_2023', "fileCnt" : 100},
        {"path": '/eos/user/a/akalinow/Data/SingleMu/Displaced_Dxy5m_pT0To1000_condRun3_131X_mcRun3_2023_realistic_v10/DisplacedMu_ch0_iPt1_Run2023_13_1_0_23_11_2023', "fileCnt" : 100},
        {"path": '/eos/user/a/akalinow/Data/SingleMu/Displaced_Dxy5m_pT0To1000_condRun3_131X_mcRun3_2023_realistic_v10/DisplacedMu_ch0_iPt2_Run2023_13_1_0_23_11_2023', "fileCnt" : 100},
        {"path": '/eos/user/a/akalinow/Data/SingleMu/Displaced_Dxy5m_pT0To1000_condRun3_131X_mcRun3_2023_realistic_v10/DisplacedMu_ch2_iPt0_Run2023_13_1_0_23_11_2023', "fileCnt" : 100},
        {"path": '/eos/user/a/akalinow/Data/SingleMu/Displaced_Dxy5m_pT0To1000_condRun3_131X_mcRun3_2023_realistic_v10/DisplacedMu_ch2_iPt1_Run2023_13_1_0_23_11_2023', "fileCnt" : 100},
        {"path": '/eos/user/a/akalinow/Data/SingleMu/Displaced_Dxy5m_pT0To1000_condRun3_131X_mcRun3_2023_realistic_v10/DisplacedMu_ch2_iPt2_Run2023_13_1_0_23_11_2023', "fileCnt" : 100},
        ]   
    

if filesNameLike == "NeutrinoGun_PU200_Alibordi" :   
    cscBx = 8 
    matchUsingPropagation  = False 
    paths = [
        {"path": '/eos/user/a/almuhamm/ZMu_Test/simPrivateProduction/NeutrinoGun_PU200_ForRateEstimation/', "fileCnt" : 10000},
        ]   
    analysisType = "rate"
    
        
print("input data paths", paths)        

if(runDebug == "DEBUG") :
    fileCnt = 2;
        
for path in paths :
    root_files = []
    for root, dirs, files in os.walk(path["path"]):
        for file in fnmatch.filter(files, '*.root'):
            root_files.append(os.path.join(root, file))  
            
    file_num = 0    
    for root_file in root_files :
        if isfile(root_file) :
            chosenFiles.append('file://' + root_file)
            file_num += 1
        else :
            print("file not found!!!!!!!: " + root_file)   
            
        if file_num >= path["fileCnt"] :
            break         
        if file_num >= fileCnt :
            break            

print("chosenFiles")
for chFile in chosenFiles:
    print(chFile)


print("chosen file count", len(chosenFiles) )

if len(chosenFiles) == 0 :
    print("no files selected!!!!!!!!!!!!!!!")
    exit

print("running version", version)
print("analysisType", analysisType)
print("outFilesName", outFilesName)


firstEv = 0#40000
#nEvents = 1000

# input files (up to 255 files accepted)
process.source = cms.Source('PoolSource',
fileNames = cms.untracked.vstring( 
    #'file:/eos/user/k/kbunkow/cms_data/SingleMuFullEta/721_FullEta_v4/SingleMu_16_p_1_1_xTE.root',
    #'file:/afs/cern.ch/user/k/kpijanow/Neutrino_Pt-2to20_gun_50.root',
    list(chosenFiles), ),
    skipEvents =  cms.untracked.uint32(0),
    inputCommands=cms.untracked.vstring(
        'keep *',
        'drop l1tEMTFHit2016Extras_simEmtfDigis_CSC_HLT',
        'drop l1tEMTFHit2016Extras_simEmtfDigis_RPC_HLT',
        'drop l1tEMTFHit2016s_simEmtfDigis__HLT',
        'drop l1tEMTFTrack2016Extras_simEmtfDigis__HLT',
        'drop l1tEMTFTrack2016s_simEmtfDigis__HLT')
)
	                    
if(runDebug == "DEBUG") :
    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20000))
else :
    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))


####Event Setup Producer
process.load('L1Trigger.L1TMuonOverlapPhase1.fakeOmtfParams_cff')
if useExtraploationAlgo :
    process.omtfParams.configXMLFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/hwToLogicLayer_0x0009.xml")

process.esProd = cms.EDAnalyzer("EventSetupRecordDataGetter",
   toGet = cms.VPSet(
      cms.PSet(record = cms.string('L1TMuonOverlapParamsRcd'),
               data = cms.vstring('L1TMuonOverlapParams'))
                   ),
   verbose = cms.untracked.bool(False)
)

process.TFileService = cms.Service("TFileService", fileName = cms.string(outFilesName + '.root'), closeFileFast = cms.untracked.bool(True) )
                                   
####OMTF Emulator
if useExtraploationAlgo :
    process.load('L1Trigger.L1TMuonOverlapPhase1.simOmtfDigis_extrapolSimple_cfi')
else :
    process.load('L1Trigger.L1TMuonOverlapPhase1.simOmtfDigis_cfi') 

if(runDebug == "DEBUG") :
    process.simOmtfDigis.dumpResultToXML = cms.bool(True)
    process.simOmtfDigis.XMLDumpFileName = cms.string("TestEvents__" + outFilesName + ".xml")
else :
    process.simOmtfDigis.dumpResultToXML = cms.bool(False)


if(runDebug == "DEBUG") :
    process.simOmtfDigis.eventCaptureDebug = cms.bool(True)
else :
    process.simOmtfDigis.eventCaptureDebug = cms.bool(False)    
#process.simOmtfDigis.simTracksTag = cms.InputTag('g4SimHits')

#needed only for the hits dumper
#process.simOmtfDigis.simTracksTag = cms.InputTag('g4SimHits')
#process.simOmtfDigis.simVertexesTag = cms.InputTag('g4SimHits')
#process.simOmtfDigis.muonMatcherFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/muonMatcherHists_100files_smoothStdDev_withOvf.root")


process.simOmtfDigis.dumpHitsToROOT = cms.bool(False)
process.simOmtfDigis.candidateSimMuonMatcher = cms.bool(False)


#process.simOmtfDigis.sorterType = cms.string("byLLH")
#process.simOmtfDigis.ghostBusterType = cms.string("byRefLayer") # byLLH byRefLayer GhostBusterPreferRefDt


#process.simOmtfDigis.patternsXMLFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/Patterns_0x0009_oldSample_3_10Files.xml")
#process.simOmtfDigis.patternsXMLFile = cms.FileInPath("L1Trigger/L1TMuonOverlapPhase1/test/expert/omtf/Patterns_0x0009_oldSample_3_10Files_classProb1.xml")
#process.simOmtfDigis.patternsXMLFile = cms.FileInPath("L1Trigger/L1TMuonOverlapPhase1/test/expert/omtf/GPs_parametrised_v1_classProb3.xml")
#process.simOmtfDigis.patternsXMLFile = cms.FileInPath("L1Trigger/L1TMuonOverlapPhase1/test/expert/omtf/Patterns_0x00012_oldSample_3_30Files_grouped1_classProb1_recalib.xml")
#process.simOmtfDigis.patternsXMLFile = cms.FileInPath("L1Trigger/L1TMuonOverlapPhase1/test/expert/omtf/Patterns_0x00012_oldSample_3_30Files_grouped1_classProb11_recalib2.xml")
#process.simOmtfDigis.patternsXMLFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/Patterns_0x00012_oldSample_3_30Files_grouped1_classProb17_recalib2.xml")

#process.simOmtfDigis.patternsXMLFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/Patterns_0x0003.xml")
#process.simOmtfDigis.patternsXMLFiles = cms.VPSet(cms.PSet(patternsXMLFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/GPs_parametrised_plus_v1.xml")),
#                                                       cms.PSet(patternsXMLFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/GPs_parametrised_minus_v1.xml"))
#)

if useExtraploationAlgo :
    #process.simOmtfDigis.patternsXMLFile = cms.FileInPath("L1Trigger/L1TMuonOverlapPhase1/test/expert/omtf/Patterns_layerStat_ExtraplMB1nadMB2_t10_classProb17_recalib2_test.xml")
    #process.simOmtfDigis.patternsXMLFile = cms.FileInPath("L1Trigger/L1TMuonOverlapPhase1/test/expert/omtf/Patterns_ExtraplMB1nadMB2Simplified_t14_classProb17_recalib2.xml")
    #process.simOmtfDigis.patternsXMLFile = cms.FileInPath("L1Trigger/L1TMuonOverlapPhase1/test/expert/omtf/Patterns_ExtraplMB1nadMB2FullAlgo_t16_classProb17_recalib2.xml")
    #process.simOmtfDigis.patternsXMLFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/Patterns_ExtraplMB1nadMB2SimplifiedFP_t17_classProb17_recalib2.xml")
    process.simOmtfDigis.patternsXMLFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/Patterns_ExtraplMB1nadMB2SimplifiedFP_t17_classProb17_recalib2_minDP0_v3.xml")
else :
    process.simOmtfDigis.patternsXMLFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/Patterns_0x00012_oldSample_3_30Files_grouped1_classProb17_recalib2.xml")

  
process.simOmtfDigis.rpcMaxClusterSize = cms.int32(3)
process.simOmtfDigis.rpcMaxClusterCnt = cms.int32(2)
process.simOmtfDigis.rpcDropAllClustersIfMoreThanMax = cms.bool(True)


process.simOmtfDigis.noHitValueInPdf = cms.bool(True)

process.simOmtfDigis.lctCentralBx = cms.int32(cscBx);#<<<<<<<<<<<<<<<<!!!!!!!!!!!!!!!!!!!!TODO this was changed in CMSSW 10(?) to 8. if the data were generated with the previous CMSSW then you have to use 6

if useExtraploationAlgo :
    process.simOmtfDigis.dtRefHitMinQuality =  cms.int32(4)

    process.simOmtfDigis.usePhiBExtrapolationFromMB1 = cms.bool(True)
    process.simOmtfDigis.usePhiBExtrapolationFromMB2 = cms.bool(True)
    
    process.simOmtfDigis.goldenPatternResultFinalizeFunction = cms.int32(10) #valid values are 0, 1, 2, 3, 5
    
    process.simOmtfDigis.minDtPhiQuality = cms.int32(2)
    process.simOmtfDigis.minDtPhiBQuality = cms.int32(4) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!!!!!!!!!!!!!!!!!!
    
    #process.simOmtfDigis.useEndcapStubsRInExtr  = cms.bool(True)   #TODO REMOVE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    #process.simOmtfDigis.useFloatingPointExtrapolation  = cms.bool(False)
    #process.simOmtfDigis.extrapolFactorsFilename = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/ExtrapolationFactors_withQAndEta.xml")
else :
    process.simOmtfDigis.minDtPhiQuality = cms.int32(2)
    process.simOmtfDigis.minDtPhiBQuality = cms.int32(2) #in 2023 it was 2, but 4 reduces the rate  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!!!!!!!!!!!!!!!!!!
         
#process.simOmtfDigis.stubEtaEncoding = cms.string("valueP1Scale")  
#process.simOmtfDigis.stubEtaEncoding = cms.string("bits")   

#nn_pThresholds = [0.36, 0.38, 0.40, 0.42, 0.44, 0.46, 0.48, 0.50, 0.52, 0.54 ]
#nn_pThresholds = [0.40, 0.50] 
#nn_pThresholds = [0.35, 0.40, 0.45, 0.50, 0.55] 
 
#process.simOmtfDigis.neuralNetworkFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/omtfClassifier_withPtBins_v34.txt")
#process.simOmtfDigis.ptCalibrationFileName = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/PtCalibration_v34.root")

#process.simOmtfDigis.nn_pThresholds = cms.vdouble(nn_pThresholds)


#process.dumpED = cms.EDAnalyzer("EventContentAnalyzer")
#process.dumpES = cms.EDAnalyzer("PrintEventSetupContent")

#process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
#process.load("Configuration.StandardSequences.MagneticField_38T_cff")

process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi")
#process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi")
#process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")


if matchUsingPropagation :
    process.L1MuonAnalyzerOmtf= cms.EDAnalyzer("L1MuonAnalyzerOmtf", 
                                 etaCutFrom = cms.double(0.82), #OMTF eta range
                                 etaCutTo = cms.double(1.24),
                                 L1OMTFInputTag  = cms.InputTag("simOmtfDigis","OMTF"),
                                 #nn_pThresholds = cms.vdouble(nn_pThresholds), 
                                 analysisType = cms.string(analysisType),
                                 
                                 simTracksTag = cms.InputTag('g4SimHits'),
                                 simVertexesTag = cms.InputTag('g4SimHits'),
                                 
                                 matchUsingPropagation = cms.bool(matchUsingPropagation),
                                 muonMatcherFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/muonMatcherHists_100files_smoothStdDev_withOvf.root") #if you want to make this file, remove this entry#if you want to make this file, remove this entry
                                 #muonMatcherFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/muonMatcherHists_noPropagation_t74.root")
                                 )
else :
    process.L1MuonAnalyzerOmtf= cms.EDAnalyzer("L1MuonAnalyzerOmtf", 
                                 etaCutFrom = cms.double(0.82), #OMTF eta range
                                 etaCutTo = cms.double(1.24),
                                 L1OMTFInputTag  = cms.InputTag("simOmtfDigis","OMTF"),
                                 #nn_pThresholds = cms.vdouble(nn_pThresholds), 
                                 analysisType = cms.string(analysisType),
                                                                  
                                 simTracksTag = cms.InputTag('g4SimHits'),
                                 simVertexesTag = cms.InputTag('g4SimHits'),
                                 matchUsingPropagation = cms.bool(matchUsingPropagation), #if this is defined, useMatcher is true, for rate analysis this mus be removed, but for efficiency is needed
                                 muonMatcherFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/muonMatcherHists_100files_smoothStdDev_withOvf.root")                                     
                                )

process.l1MuonAnalyzerOmtfPath = cms.Path(process.L1MuonAnalyzerOmtf)


process.L1TMuonSeq = cms.Sequence( process.esProd          
                                   + process.simOmtfDigis 
                                   #+ process.dumpED
                                   #+ process.dumpES
)

process.L1TMuonPath = cms.Path(process.L1TMuonSeq)

process.schedule = cms.Schedule(process.L1TMuonPath, process.l1MuonAnalyzerOmtfPath)

#process.out = cms.OutputModule("PoolOutputModule", 
#   fileName = cms.untracked.string("l1tomtf_superprimitives1.root")
#)

#process.output_step = cms.EndPath(process.out)
#process.schedule = cms.Schedule(process.L1TMuonPath)
#process.schedule.extend([process.output_step])
