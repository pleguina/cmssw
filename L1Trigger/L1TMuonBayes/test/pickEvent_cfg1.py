import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('analysis')
# add a list of strings for events to process
# options.register ('eventsToProcess',
#                                   '',
#                                   VarParsing.multiplicity.list,
#                                   VarParsing.varType.string,
#                                   "Events to process")
# options.parseArguments()


Source_Files = cms.untracked.vstring(
#        "/store/relval/CMSSW_10_0_0_pre1/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/94X_upgrade2023_realistic_v2_2023D17noPU-v2/10000/06C888F3-CFCE-E711-8928-0CC47A4D764C.root"
         #"/store/relval/CMSSW_9_3_2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/93X_upgrade2023_realistic_v2_2023D17noPU-v1/10000/0681719F-AFA6-E711-87C9-0CC47A4C8E14.root"
         #"file:///eos/user/k/kbunkow/cms_data/0681719F-AFA6-E711-87C9-0CC47A4C8E14.root"
         #"file:///eos/cms/store/group/upgrade/sandhya/SMP-PhaseIIFall17D-00001.root"
         #'file:///afs/cern.ch/work/k/kbunkow/private/omtf_data/SingleMu_15_p_1_1_qtl.root' no high eta in tis file
         #'file:///eos/user/k/kbunkow/cms_data/SingleMuFullEta/721_FullEta_v4/SingleMu_12_p_10_1_mro.root' ,
         #'file:///eos/user/k/kbunkow/cms_data/SingleMuFullEta/721_FullEta_v4/SingleMu_20_p_118_1_sTk.root' ,
         #'file:///eos/user/k/kbunkow/cms_data/SingleMuFullEta/721_FullEta_v4/SingleMu_5_p_81_1_Ql3.root',
         #'file:///eos/user/k/kbunkow/cms_data/SingleMuFullEta/721_FullEta_v4/SingleMu_31_p_89_2_MJS.root',
         "/store/mc/PhaseIIFall17D/SingleMu_FlatPt-2to100/GEN-SIM-DIGI-RAW/L1TnoPU_93X_upgrade2023_realistic_v5-v1/00000/F4EEAE55-C937-E811-8C29-48FD8EE739D1.root"
)

process = cms.Process("PickEvent")
process.source = cms.Source ("PoolSource",
          fileNames = Source_Files, #cms.untracked.vstring (options.inputFiles),
          #eventsToProcess = cms.untracked.VEventRange (options.eventsToProcess),
          inputCommands=cms.untracked.vstring(
          'keep *',
          'drop l1tEMTFHit2016Extras_simEmtfDigis_CSC_HLT',
          'drop l1tEMTFHit2016Extras_simEmtfDigis_RPC_HLT',
          'drop l1tEMTFHit2016s_simEmtfDigis__HLT',
          'drop l1tEMTFTrack2016Extras_simEmtfDigis__HLT',
          'drop l1tEMTFTrack2016s_simEmtfDigis__HLT')                               
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))

#outputFileNme = '/eos/user/k/kbunkow/cms_data/mc/PhaseIIFall17D/SingleMu_FlatPt-2to100/GEN-SIM-DIGI-RAW/L1TnoPU_93X_upgrade2023_realistic_v5-v1/00000/F4EEAE55-C937-E811-8C29-48FD8EE739D1_dump1000Events.root'
outputFileNme = 'F4EEAE55-C937-E811-8C29-48FD8EE739D1_dump1000Events.root'

process.Out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string (outputFileNme)
)

process.end = cms.EndPath(process.Out)
