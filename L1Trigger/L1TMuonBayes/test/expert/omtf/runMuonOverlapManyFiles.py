# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms
process = cms.Process("L1TMuonEmulation")
import os
import sys
import commands
from os import listdir
from os.path import isfile, join

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))

path = '/afs/cern.ch/work/k/kbunkow/public/data/SingleMuFullEta/721_FullEta_v4/'

onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
chosenFiles = ['file://' + path + f for f in onlyfiles if (('_p_10_' in f) or ('_m_10_' in f))]
#chosenFiles = ['file://' + path + f for f in onlyfiles if (('_5_p_10_' in f))]
#chosenFiles = ['file://' + path + f for f in onlyfiles if (re.match('.*_._p_10.*', f))]
firstEv = 0#40000
nEvents = 10000

# input files (up to 255 files accepted)
process.source = cms.Source('PoolSource',
fileNames = cms.untracked.vstring( 
    #'file:/eos/user/k/kbunkow/cms_data/SingleMuFullEta/721_FullEta_v4/SingleMu_16_p_1_1_xTE.root',
    #'file:/afs/cern.ch/user/k/kpijanow/Neutrino_Pt-2to20_gun_50.root',
    set(chosenFiles),
                                  ),
eventsToProcess = cms.untracked.VEventRange(
 '3:' + str(firstEv) + '-3:' +   str(firstEv + nEvents),
 '4:' + str(firstEv) + '-4:' +   str(firstEv + nEvents),
 '5:' + str(firstEv) + '-5:' +   str(firstEv + nEvents),
 '6:' + str(firstEv) + '-6:' +   str(firstEv + nEvents),
 '7:' + str(firstEv) + '-7:' +   str(firstEv + nEvents),
 '8:' + str(firstEv) + '-8:' +   str(firstEv + nEvents),
 '9:' + str(firstEv) + '-9:' +   str(firstEv + nEvents),
'10:' + str(firstEv) + '-10:' +  str(firstEv + nEvents),
'11:' + str(firstEv) + '-11:' +  str(firstEv + nEvents),
'12:' + str(firstEv) + '-12:' +  str(firstEv + nEvents),
'13:' + str(firstEv) + '-13:' +  str(firstEv + nEvents),
'14:' + str(firstEv) + '-14:' +  str(firstEv + nEvents),
'15:' + str(firstEv) + '-15:' +  str(firstEv + nEvents),
'16:' + str(firstEv) + '-16:' +  str(firstEv + nEvents),
'17:' + str(firstEv) + '-17:' +  str(firstEv + nEvents),
'18:' + str(firstEv) + '-18:' +  str(firstEv + nEvents),
'19:' + str(firstEv) + '-19:' +  str(firstEv + nEvents),
'20:' + str(firstEv) + '-20:' +  str(firstEv + nEvents),
'21:' + str(firstEv) + '-21:' +  str(firstEv + nEvents),
'22:' + str(firstEv) + '-22:' +  str(firstEv + nEvents),
'23:' + str(firstEv) + '-23:' +  str(firstEv + nEvents),
'24:' + str(firstEv) + '-24:' +  str(firstEv + nEvents),
'25:' + str(firstEv) + '-25:' +  str(firstEv + nEvents),
'26:' + str(firstEv) + '-26:' +  str(firstEv + nEvents),
'27:' + str(firstEv) + '-27:' +  str(firstEv + nEvents),
'28:' + str(firstEv) + '-28:' +  str(firstEv + nEvents),
'29:' + str(firstEv) + '-29:' +  str(firstEv + nEvents),
'30:' + str(firstEv) + '-30:' +  str(firstEv + nEvents),
'31:' + str(firstEv) + '-31:' +  str(firstEv + nEvents)),
skipEvents =  cms.untracked.uint32(0)
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))

# PostLS1 geometry used
process.load('Configuration.Geometry.GeometryExtended2015Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2015_cff')
############################
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')


####Event Setup Producer
process.load('L1Trigger.L1TMuonOverlap.fakeOmtfParams_cff')
process.esProd = cms.EDAnalyzer("EventSetupRecordDataGetter",
   toGet = cms.VPSet(
      cms.PSet(record = cms.string('L1TMuonOverlapParamsRcd'),
               data = cms.vstring('L1TMuonOverlapParams'))
                   ),
   verbose = cms.untracked.bool(False)
)


####OMTF Emulator
process.load('L1Trigger.L1TMuonOverlap.simOmtfDigis_cfi')

process.dumpED = cms.EDAnalyzer("EventContentAnalyzer")
process.dumpES = cms.EDAnalyzer("PrintEventSetupContent")

process.cwiczenie= cms.EDAnalyzer("Cwiczenie", 
                                  outRootFile = cms.string("omtfAnalysis_oldPatsOldGb.root") )

process.L1TMuonSeq = cms.Sequence( process.esProd          
                                   + process.simOmtfDigis 
                                   + process.cwiczenie
                                   #+ process.dumpED
                                   #+ process.dumpES
)

process.L1TMuonPath = cms.Path(process.L1TMuonSeq)

process.out = cms.OutputModule("PoolOutputModule", 
   fileName = cms.untracked.string("l1tomtf_superprimitives1.root")
)

process.output_step = cms.EndPath(process.out)
process.schedule = cms.Schedule(process.L1TMuonPath)
process.schedule.extend([process.output_step])
