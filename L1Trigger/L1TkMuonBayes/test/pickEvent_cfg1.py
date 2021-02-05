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
         #"/store/mc/PhaseIIFall17D/SingleMu_FlatPt-2to100/GEN-SIM-DIGI-RAW/L1TnoPU_93X_upgrade2023_realistic_v5-v1/00000/F4EEAE55-C937-E811-8C29-48FD8EE739D1.root"
        #"/store/mc/PhaseIIFall17D/HSCPppstau_M_871_TuneCUETP8M1_14TeV_pythia8/GEN-SIM-DIGI-RAW/L1TnoPU_93X_upgrade2023_realistic_v5-v1/120000/18156A80-66EC-E811-AE02-0CC47AFCC62A.root"
        #"/store/mc/PhaseIIFall17D/HSCPppstau_M_200_TuneCUETP8M1_14TeV_pythia8/GEN-SIM-DIGI-RAW/L1TnoPU_93X_upgrade2023_realistic_v5-v1/120000/FE3D8AD6-B6D0-E811-8FBD-141877412793.root"
        #'/store/mc/RunIISummer16DR80Premix/HSCPppstau_M-651_TuneZ2star_13TeV_pythia6/AODSIM/PUMoriond17_HSCP_customise_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/0E0D542C-A9C8-E611-981C-A0000420FE80.root'
        #'/store/mc/RunIISummer16DR80Premix/HSCPppstau_M-308_TuneZ2star_13TeV-pythia6/AODSIM/PUMoriond17_HSCP_customise_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/3EAEE028-AED2-E611-83E8-002590E7DFEE.root'
        #'/store/mc/RunIISummer16DR80Premix/HSCPppstau_M-308_TuneZ2star_13TeV-pythia6/AODSIM/PUMoriond17_HSCP_customise_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/9A3F1E44-AED2-E611-9AD6-1CC1DE18CFDE.root'
        #'/store/mc/PhaseIIFall17D/ZMM_14TeV_TuneCUETP8M1_Pythia8/GEN-SIM-DIGI-RAW/L1TnoPU_93X_upgrade2023_realistic_v5-v1/60000/EE29AF8E-51AF-E811-A2BD-484D7E8DF0D3.root'
        #'/store/mc/PhaseIITDRSpring17GS/HSCPppstau_M_1599_TuneCUETP8M1_14TeV_pythia8_Customised/GEN-SIM/Customised_91X_upgrade2023_realistic_v3-v2/00000/CC5458A2-5BA9-E711-86E1-0025905D1D7A.root'
        #'/store/mc/PhaseIIFall17D/SingleMu_FlatPt-2to100/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/00000/32DF01CC-A342-E811-9FE7-48D539F3863E.root'
        #'/store/relval/CMSSW_9_3_7/RelValDisplacedMuonGun_Pt30To100_Dxy_0_1000/GEN-SIM-DIGI-RAW/93X_upgrade2023_realistic_v5_2023D17noPULxyMx1k-v1/10000/E68C6334-7F62-E911-8AA5-0025905B8610.root'
         #'/store/relval/CMSSW_9_3_7/RelValDisplacedMuonGun_Pt30To100_Dxy_0_1000/GEN-SIM-DIGI-RAW/93X_upgrade2023_realistic_v5_2023D17noPULxyMx1k-v1/10000/C2D5C228-7F62-E911-AAA6-0CC47A78A42C.root'
         #'/store/relval/CMSSW_9_3_7/RelValDisplacedMuonGun_Pt30To100_Dxy_0_1000/GEN-SIM-DIGI-RAW/93X_upgrade2023_realistic_v5_2023D17noPULxyMx1k-v1/10000/A2885F26-7F62-E911-84CF-0CC47A4D7666.root'
         #'/store/relval/CMSSW_9_3_7/RelValDisplacedMuonGun_Pt30To100_Dxy_0_1000/GEN-SIM-DIGI-RAW/93X_upgrade2023_realistic_v5_2023D17noPULxyMx1k-v1/10000/3C065D44-7E62-E911-9F3F-0CC47A4C8F12.root'
         #'/store/relval/CMSSW_9_3_7/RelValDisplacedMuonGun_Pt30To100_Dxy_0_1000/GEN-SIM-DIGI-RAW/93X_upgrade2023_realistic_v5_2023D17noPULxyMx1k-v1/10000/34016F34-7F62-E911-AAB8-0025905AA9CC.root'
         #'/store/mc/PhaseIITDRSpring19DR/GluGluHToZZTo4L_M125_14TeV_powheg2_JHUgenV724_pythia8/GEN-SIM-DIGI-RAW/NoPU_106X_upgrade2023_realistic_v3-v2/260000/FB98030B-16C5-9842-9698-8371EB8D8B01.root'
         #'/store/mc/PhaseIITDRSpring19DR/HSCPppstau_M_871_TuneCUETP8M1_14TeV_pythia8/GEN-SIM-DIGI-RAW/NoPU_106X_upgrade2023_realistic_v3-v1/40000/F7B16025-96C6-7F49-97CD-E2650DDA67C1.root'
         #'/store/mc/PhaseIITDRSpring19DR/HSCPppstau_M_200_TuneCUETP8M1_14TeV_pythia8/GEN-SIM-DIGI-RAW/NoPU_106X_upgrade2023_realistic_v3-v2/50000/F9EFA962-CDD8-C643-8F62-8F75384875F0.root'
         #'/store/mc/PhaseIITDRSpring19DR/JPsiToMuMu_Pt0to100-pythia8-gun/GEN-SIM-DIGI-RAW/NoPU_106X_upgrade2023_realistic_v3-v1/40000/FDA71CB6-4C3B-4540-99EB-803077C6EC2D.root'
         #'/store/mc/PhaseIIMTDTDRAutumn18DR/BsToMuMu_SoftQCDnonD_TuneCP5_14TeV-pythia8-evtgen/FEVT/NoPU_103X_upgrade2023_realistic_v2-v1/40000/F8FB78DD-B8E0-834B-9C0B-853BB0B9C13B.root'
         #"/store/mc/PhaseIITDRSpring19DR/HSCPppstau_M_871_TuneCP5_14TeV_pythia8/GEN-SIM-DIGI-RAW/PU200_106X_upgrade2023_realistic_v3-v2/60000/FFCC029A-296A-3547-BFF3-23D06CAC0E78.root"
         #"/store/mc/PhaseIITDRSpring19DR/HSCPppstau_M_871_TuneCP5_14TeV_pythia8/GEN-SIM-DIGI-RAW/NoPU_106X_upgrade2023_realistic_v3-v2/60000/F9357CE3-E1BD-C64C-8F43-895CFA3A0AFC.root"
         #"/store/mc/PhaseIITDRSpring19DR/HSCPppstau_M_200_TuneCP5_14TeV_pythia8/GEN-SIM-DIGI-RAW/NoPU_106X_upgrade2023_realistic_v3-v2/60000/A91AA4D8-5187-5544-8304-365404899406.root"
         #"/store/mc/PhaseIITDRSpring19DR/Mu_FlatPt2to100-pythia8-gun/GEN-SIM-DIGI-RAW/NoPU_106X_upgrade2023_realistic_v3-v1/60000/E0D5C6A5-B855-D14F-9124-0B2C9B28D0EA.root"
         #"/store/mc/PhaseIITDRSpring19DR/HSCPppstau_M_200_TuneCP5_14TeV_pythia8/GEN-SIM-DIGI-RAW/PU200_HSCP_106X_upgrade2023_realistic_v3-v1/40000/ACF9C8E2-0570-6A4A-983A-E2B230F6FCAA.root"
         #"/store/mc/PhaseIITDRSpring19DR/HSCPppstau_M_871_TuneCP5_14TeV_pythia8/GEN-SIM-DIGI-RAW/PU200_HSCP_106X_upgrade2023_realistic_v3-v2/260000/1ADE9D9E-8C0C-1948-A405-5DFDA1AF5172.root"
         #"/store/mc/PhaseIITDRSpring19DR/Nu_E10-pythia8-gun/GEN-SIM-DIGI-RAW/PU200_106X_upgrade2023_realistic_v3-v3/70001/FFB3195D-E113-3744-877D-44E21C060358.root"
         #"/store/mc/PhaseIITDRSpring19DR/Nu_E10-pythia8-gun/GEN-SIM-DIGI-RAW/PU250_106X_upgrade2023_realistic_v3_ext2-v1/50002/FFE07316-3810-6E44-97A1-5753A3070D12.root"
         #"/store/mc/Phase2HLTTDRWinter20DIGI/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW/NoPU_110X_mcRun4_realistic_v3-v2/50000/CFBA4AE8-00DE-1743-9B83-7582C69FC7F7.root"
         #"/store/mc/Phase2HLTTDRWinter20DIGI/Muminus_Pt10-gun/GEN-SIM-DIGI-RAW/NoPU_110X_mcRun4_realistic_v3-v2/50000/E6F1BC5E-BD51-A948-ADDC-8D84EFF14174.root"
         #'/store/mc/Phase2HLTTDRWinter20DIGI/DoubleMuon_gun_FlatPt-1To100/GEN-SIM-DIGI-RAW/NoPU_110X_mcRun4_realistic_v3-v2/30000/3FD40D17-5C29-804C-B49A-029CC02B63DC.root'
         '/store/mc/Phase2HLTTDRWinter20DIGI/DoubleMuon_gun_FlatPt-1To100/GEN-SIM-DIGI-RAW/NoPU_110X_mcRun4_realistic_v3-v2/30000/D19F8A25-CD3A-674B-B856-DFBD7B445493.root'
         
)

process = cms.Process("PickEvent")
process.source = cms.Source ("PoolSource",
          fileNames = Source_Files, #cms.untracked.vstring (options.inputFiles),
          #eventsToProcess = cms.untracked.VEventRange (options.eventsToProcess),
          dropDescendantsOfDroppedBranches=cms.untracked.bool(False),
          inputCommands=cms.untracked.vstring(
          'keep *',
          'drop l1tEMTFHit2016Extras_simEmtfDigis_CSC_HLT',
          'drop l1tEMTFHit2016Extras_simEmtfDigis_RPC_HLT',
          'drop l1tEMTFHit2016s_simEmtfDigis__HLT',
          'drop l1tEMTFTrack2016Extras_simEmtfDigis__HLT',
          'drop l1tEMTFTrack2016s_simEmtfDigis__HLT',
          'drop l1tHGCalClusterBXVector_hgcalTriggerPrimitiveDigiProducer_cluster2D_HLT',
          'drop *HGCal*_*_*_*',
          'drop *hgcal*_*_*_*',
          'drop *Ecal*_*_*_*',
          'drop *Hcal*_*_*_*',
          'drop *Calo*_*_*_*',
          
          'drop *_*HGCal*_*_*',
          'drop *_*hgcal*_*_*',
          'drop *_*Ecal*_*_*',
          'drop *_*Hcal*_*_*',
          'drop *_*Calo*_*_*',
          
          'drop *_*_*HGCal*_*',
          'drop *_*_*hgcal*_*',
          'drop *_*_*Ecal*_*',
          'drop *_*_*Hcal*_*',
          'drop *_*_*Calo*_*',
          )                               
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

#outputFileNme = '/eos/user/k/kbunkow/cms_data/mc/PhaseIIFall17D/SingleMu_FlatPt-2to100/GEN-SIM-DIGI-RAW/L1TnoPU_93X_upgrade2023_realistic_v5-v1/00000/F4EEAE55-C937-E811-8C29-48FD8EE739D1_dump1000Events.root'
#outputFileNme = 'HSCPppstau_M-651_TuneZ2star_13TeV_0E0D542C-A9C8-E611-981C-A0000420FE80_dump100Events.root'
#outputFileNme = 'SingleMu_PU200_32DF01CC-A342-E811-9FE7-48D539F3863E_dump500Events.root'
#outputFileNme = 'HSCPppstau_M_871_PU200_v3-v2_1ADE9D9E-8C0C-1948-A405-5DFDA1AF5172_dump100Ev.root'
#outputFileNme = 'Nu_E10-pythia8-gun_PU200_v3-v3_FFB3195D-E113-3744-877D-44E21C060358_dump100Ev.root'
#outputFileNme = 'Phase2HLTTDRWinter20DIGI__ZprimeToMuMu_M-6000_NoPU_CFBA4AE8-00DE-1743-9B83-7582C69FC7F7_dump100Ev.root'
outputFileNme = 'Phase2HLTTDRWinter20DIGI_DoubleMuon_gun_FlatPt-1To100_NoPU_3FD40D17-5C29-804C-B49A-029CC02B63DC_dump100Ev.root'

process.Out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string (outputFileNme)
)

process.end = cms.EndPath(process.Out)
