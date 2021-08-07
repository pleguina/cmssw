from CRABClient.UserUtilities import config #, getUsernameFromSiteDB
config = config()

type = 'omtf_nn'
#type = 'omtf_0x0006'

config.General.requestName = type + '_MC_analysis_ZprimeToMuMu_NoPU_v3_t100_1'
#config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
#config.JobType.psetName = 'runMuonOverlapTTMergerAnalyzerCrab.py'

if type == 'omtf_nn' :
    config.JobType.psetName = 'runMuonOverlap_nn_phase2.py'
else :
    config.JobType.psetName = 'runMuonOverlap_0x0006.py'

config.JobType.pyCfgParams = ['efficiency']

config.Data.inputDataset = '/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/PhaseIITDRSpring19DR-NoPU_106X_upgrade2023_realistic_v3-v2/GEN-SIM-DIGI-RAW' 
#'/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/Run3Winter20DRMiniAOD-FlatPU0to80_110X_mcRun3_2021_realistic_v6-v3/GEN-SIM-RAW'  <-- no DIGI here
#'/SingleMu_FlatPt-2to100/PhaseIIFall17D-L1TPU200_93X_upgrade2023_realistic_v5-v1/GEN-SIM-DIGI-RAW'
#config.Data.inputDataset = '/SingleMu_FlatPt-2to100/PhaseIIFall17D-L1TnoPU_93X_upgrade2023_realistic_v5-v1/GEN-SIM-DIGI-RAW'

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
#config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = 'CRAB_' + type + '_MC_analysis_ZprimeToMuMu_NoPu_v3_t100_1'
config.Data.totalUnits = 50
config.Data.ignoreLocality = False

config.section_("Debug")
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']

#config.Site.storageSite = 'T2_PL_Swierk'
config.Site.storageSite = 'T2_CH_CERN'

