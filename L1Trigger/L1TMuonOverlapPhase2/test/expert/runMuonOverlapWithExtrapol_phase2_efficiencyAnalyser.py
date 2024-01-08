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
    version = version + 'Patterns_ExtraplMB1nadMB2DTQualAndEtaFixedP_ValueP1Scale_t20_v1_SingleMu_iPt_and_OneOverPt_classProb17_recalib2_minDP0'
else :
    version = version + 'Patterns_0x00012'

runDebug = "DEBUG" # or "INFO" DEBUG
#useExtraploationAlgo = True

regeneratedL1DT = False

analysisType = "efficiency" # or rate

useNN = True

if useNN :
   version = version + "_NN_FP_v217"
   
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
    
if "NeutrinoGun" in filesNameLike : 
    outFilesName = 'omtfAnalysis2_'  + "rate_"
    
outFilesName = outFilesName + version + "__" + filesNameLike

if(runDebug == "DEBUG") :
    outFilesName = outFilesName + "_test4"

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
                         threshold = cms.untracked.string("DEBUG"), #DEBUG
                         default = cms.untracked.PSet( limit = cms.untracked.int32(0) ), 
                         #INFO   =  cms.untracked.int32(0),
                         #DEBUG   = cms.untracked.int32(0),
                         l1tOmtfEventPrint = cms.untracked.PSet( limit = cms.untracked.int32(1000000000) ),
                         OMTFReconstruction = cms.untracked.PSet( limit = cms.untracked.int32(1000000000) ),
                         l1MuonAnalyzerOmtf = cms.untracked.PSet( limit = cms.untracked.int32(1000000000) ),
                       ),
       debugModules = cms.untracked.vstring('L1MuonAnalyzerOmtf', 'simOmtfPhase2Digis', 'omtfParamsSource', 'omtfParams', "esProd", 'L1TMuonOverlapPhase1ParamsESProducer') #'L1MuonAnalyzerOmtf',
       #debugModules = cms.untracked.vstring('*')
    )

    #process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)
if not verbose:
    process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
    process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False), 
                                         #SkipEvent = cms.untracked.vstring('ProductNotFound') 
                                     )
    
    
# PostLS1 geometry used
process.load('Configuration.Geometry.GeometryExtended2026D86Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D86_cff') 
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

if filesNameLike == 'Displaced_cTau5m_XTo2LLTo4Mu':
    matchUsingPropagationInAnlyzer = True 
    matchUsingPropagationInDumper = True 
    paths = [    
             "/eos/user/a/almuhamm/ZMu_Test/simPrivateProduction/Displaced_cTau5m_XTo2LLTo4Mu_condPhase2_realistic/XTo2LLPTo4Mu_CTau5m_Phase2Exotic/231203_175643/0000/",  # 500 files
             ]    

if filesNameLike == 'Displaced_Dxy3m_pT0To1000_condPhase2_realistic':
    matchUsingPropagationInAnlyzer = True 
    matchUsingPropagationInDumper = False 
    paths = [    
             "/eos/user/a/almuhamm/ZMu_Test/simPrivateProduction/Displaced_Dxy3m_pT0To1000_condPhase2_realistic/DisplacedMu_ch0_iPt0_Run2029_13_1_0_01_12_2023",  # 500 files
             "/eos/user/a/almuhamm/ZMu_Test/simPrivateProduction/Displaced_Dxy3m_pT0To1000_condPhase2_realistic/DisplacedMu_ch2_iPt0_Run2029_13_1_0_01_12_2023",  # 500 files
             "/eos/user/a/almuhamm/ZMu_Test/simPrivateProduction/Displaced_Dxy3m_pT0To1000_condPhase2_realistic/DisplacedMu_ch0_iPt1_Run2029_13_1_0_01_12_2023",  # 500 files
             "/eos/user/a/almuhamm/ZMu_Test/simPrivateProduction/Displaced_Dxy3m_pT0To1000_condPhase2_realistic/DisplacedMu_ch2_iPt1_Run2029_13_1_0_01_12_2023",  # 500 files
             "/eos/user/a/almuhamm/ZMu_Test/simPrivateProduction/Displaced_Dxy3m_pT0To1000_condPhase2_realistic/DisplacedMu_ch0_iPt2_Run2029_13_1_0_01_12_2023",  # 500 files
             "/eos/user/a/almuhamm/ZMu_Test/simPrivateProduction/Displaced_Dxy3m_pT0To1000_condPhase2_realistic/DisplacedMu_ch2_iPt2_Run2029_13_1_0_01_12_2023",  # 500 files
             ]   

if filesNameLike == "NeutrinoGun_PU200_Alibordi" :   
    cscBx = 8 
    matchUsingPropagation  = False 
    paths = [
        {"path": '/eos/user/a/almuhamm/MuSampleSharedDirectory/simPrivateProduction/NeutrinoGun_PU200_ForRateEstimation/', "fileCnt" : 10000},
        ]   
    analysisType = "rate"

if filesNameLike == "MinBias_Phase2Spring23_PU140" :   
    cscBx = 8 
    matchUsingPropagation  = False 
    regeneratedL1DT = True
    paths = [
        {"path": '/eos/user/a/akalinow/Data/SingleMu/MinBias_TuneCP5_14TeV-pythia8/crab_MinBias_TuneCP5_14TeV-pythia8_Phase2Spring23DIGIRECOMiniAOD-PU140/', "fileCnt" : 10000},
        ]   
    analysisType = "rate"    
        
print("input data paths", paths)        

if(runDebug == "DEBUG") :
    fileCnt = 1;
        
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
    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500))
else :
    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))

# Calibrate Digis
process.load("L1Trigger.DTTriggerPhase2.CalibratedDigis_cfi")
process.CalibratedDigis.dtDigiTag = "simMuonDTDigis" 
process.CalibratedDigis.scenario = 0

# DTTriggerPhase2
process.load("L1Trigger.DTTriggerPhase2.dtTriggerPhase2PrimitiveDigis_cfi")
process.dtTriggerPhase2PrimitiveDigis.debug = False
process.dtTriggerPhase2PrimitiveDigis.dump = False
process.dtTriggerPhase2PrimitiveDigis.scenario = 0

####Event Setup Producer
process.load('L1Trigger.L1TMuonOverlapPhase1.fakeOmtfParams_cff')
process.omtfParams.configXMLFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/hwToLogicLayer_0x0209.xml")
process.omtfParams.patternsXMLFiles = cms.VPSet(
        cms.PSet(patternsXMLFile=cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/Patterns_ExtraplMB1nadMB2DTQualAndEtaFixedP_ValueP1Scale_t20_v1_SingleMu_iPt_and_OneOverPt_classProb17_recalib2_minDP0.xml")),)

process.esProd = cms.EDAnalyzer("EventSetupRecordDataGetter",
   toGet=cms.VPSet(
      cms.PSet(record=cms.string('L1TMuonOverlapParamsRcd'),
               data=cms.vstring('L1TMuonOverlapParams'))
                   ),
   verbose=cms.untracked.bool(False)
)

process.TFileService = cms.Service("TFileService", fileName = cms.string(outFilesName + '.root'), closeFileFast = cms.untracked.bool(True) )
                                   
####OMTF Emulator
if useExtraploationAlgo :
    process.load('L1Trigger.L1TMuonOverlapPhase2.simOmtfPhase2Digis_extrapol_cfi')
else :
    process.load('L1Trigger.L1TMuonOverlapPhase2.simOmtfPhase2Digis_cfi')

if(runDebug == "DEBUG") :
    process.simOmtfPhase2Digis.dumpResultToXML = cms.bool(True)
    process.simOmtfPhase2Digis.XMLDumpFileName = cms.string("TestEvents__" + outFilesName + ".xml")
else :
    process.simOmtfPhase2Digis.dumpResultToXML = cms.bool(False)


if(runDebug == "DEBUG") :
    process.simOmtfPhase2Digis.eventCaptureDebug = cms.bool(True)
else :
    process.simOmtfPhase2Digis.eventCaptureDebug = cms.bool(False)    
#process.simOmtfPhase2Digis.simTracksTag = cms.InputTag('g4SimHits')

#needed only for the hits dumper
#process.simOmtfPhase2Digis.simTracksTag = cms.InputTag('g4SimHits')
#process.simOmtfPhase2Digis.simVertexesTag = cms.InputTag('g4SimHits')
#process.simOmtfPhase2Digis.muonMatcherFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/muonMatcherHists_100files_smoothStdDev_withOvf.root")


process.simOmtfPhase2Digis.dumpHitsToROOT = cms.bool(False)
process.simOmtfPhase2Digis.candidateSimMuonMatcher = cms.bool(False)


#process.simOmtfPhase2Digis.sorterType = cms.string("byLLH")
#process.simOmtfPhase2Digis.ghostBusterType = cms.string("byRefLayer") # byLLH byRefLayer GhostBusterPreferRefDt


if useExtraploationAlgo :
    process.simOmtfPhase2Digis.patternsXMLFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/Patterns_ExtraplMB1nadMB2DTQualAndEtaFixedP_ValueP1Scale_t20_v1_SingleMu_iPt_and_OneOverPt_classProb17_recalib2_minDP0.xml")
else :
    process.simOmtfPhase2Digis.patternsXMLFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/Patterns_0x00012_oldSample_3_30Files_grouped1_classProb17_recalib2.xml") ##todo

  
process.simOmtfPhase2Digis.rpcMaxClusterSize = cms.int32(3)
process.simOmtfPhase2Digis.rpcMaxClusterCnt = cms.int32(2)
process.simOmtfPhase2Digis.rpcDropAllClustersIfMoreThanMax = cms.bool(True)


process.simOmtfPhase2Digis.noHitValueInPdf = cms.bool(True)

process.simOmtfPhase2Digis.lctCentralBx = cms.int32(cscBx);#<<<<<<<<<<<<<<<<!!!!!!!!!!!!!!!!!!!!TODO this was changed in CMSSW 10(?) to 8. if the data were generated with the previous CMSSW then you have to use 6

if useExtraploationAlgo :
    process.simOmtfPhase2Digis.dtRefHitMinQuality =  cms.int32(4)

    process.simOmtfPhase2Digis.usePhiBExtrapolationFromMB1 = cms.bool(True)
    process.simOmtfPhase2Digis.usePhiBExtrapolationFromMB2 = cms.bool(True)
    
    #process.simOmtfPhase2Digis.goldenPatternResultFinalizeFunction = cms.int32(10) #valid values are 0, 1, 2, 3, 5
    
    process.simOmtfPhase2Digis.minDtPhiQuality = cms.int32(2)
    process.simOmtfPhase2Digis.minDtPhiBQuality = cms.int32(2) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!!!!!!!!!!!!!!!!!!
    
    #process.simOmtfPhase2Digis.useEndcapStubsRInExtr  = cms.bool(True)   #TODO REMOVE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    #process.simOmtfPhase2Digis.useFloatingPointExtrapolation  = cms.bool(False)
    #process.simOmtfPhase2Digis.extrapolFactorsFilename = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/ExtrapolationFactors_withQAndEta.xml")
else :
    process.simOmtfPhase2Digis.minDtPhiQuality = cms.int32(2)
    process.simOmtfPhase2Digis.minDtPhiBQuality = cms.int32(2) #in 2023 it was 2, but 4 reduces the rate  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!!!!!!!!!!!!!!!!!!
         

if useNN :
    process.simOmtfPhase2Digis.neuralNetworkFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/lutNN_omtfRegression_FP_v217.xml")


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
                                 L1OMTFInputTag  = cms.InputTag("simOmtfPhase2Digis","OMTF"),
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
                                 L1OMTFInputTag  = cms.InputTag("simOmtfPhase2Digis","OMTF"),
                                 #nn_pThresholds = cms.vdouble(nn_pThresholds), 
                                 analysisType = cms.string(analysisType),
                                                                  
                                 simTracksTag = cms.InputTag('g4SimHits'),
                                 simVertexesTag = cms.InputTag('g4SimHits'),
                                 #matchUsingPropagation = cms.bool(matchUsingPropagation), #if this is defined, useMatcher is true, for rate analysis this mus be removed, but for efficiency is needed
                                 muonMatcherFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/muonMatcherHists_100files_smoothStdDev_withOvf.root")                                     
                                )

process.l1MuonAnalyzerOmtfPath = cms.Path(process.L1MuonAnalyzerOmtf)


process.L1TMuonSeq = cms.Sequence( process.esProd          
                                   + process.simOmtfPhase2Digis 
                                   #+ process.dumpED
                                   #+ process.dumpES
)

#process.L1TMuonPath = cms.Path(process.L1TMuonSeq) ########################################<<<<<<!!!!!!!!!!!!!!!!!!!!!!!!!!!
process.L1TMuonPath = cms.Path(process.CalibratedDigis * 
                            process.dtTriggerPhase2PrimitiveDigis * process.L1TMuonSeq)

process.schedule = cms.Schedule(process.L1TMuonPath, process.l1MuonAnalyzerOmtfPath)

#process.out = cms.OutputModule("PoolOutputModule", 
#   fileName = cms.untracked.string("l1tomtf_superprimitives1.root")
#)

#process.output_step = cms.EndPath(process.out)
#process.schedule = cms.Schedule(process.L1TMuonPath)
#process.schedule.extend([process.output_step])
