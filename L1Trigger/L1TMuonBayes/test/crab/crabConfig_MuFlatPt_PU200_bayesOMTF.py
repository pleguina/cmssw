from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'omtf_MC_analysis_MuFlatPt_PU200_v1_t16'
#config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
#config.JobType.psetName = 'runMuonOverlapTTMergerAnalyzerCrab.py'
config.JobType.psetName = 'runMuonOverlap.py'
config.JobType.pyCfgParams = ['efficiency']

config.Data.inputDataset = '/Mu_FlatPt2to100-pythia8-gun/PhaseIITDRSpring19DR-PU200_106X_upgrade2023_realistic_v3-v2/GEN-SIM-DIGI-RAW' 
#'/SingleMu_FlatPt-2to100/PhaseIIFall17D-L1TPU200_93X_upgrade2023_realistic_v5-v1/GEN-SIM-DIGI-RAW'
#config.Data.inputDataset = '/SingleMu_FlatPt-2to100/PhaseIIFall17D-L1TnoPU_93X_upgrade2023_realistic_v5-v1/GEN-SIM-DIGI-RAW'

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = 'CRAB3_omtf_MC_analysis_MuFlatPt_PU200_v1_t16'
config.Data.totalUnits = 50 #321
config.Data.ignoreLocality = False

config.section_("Debug")
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']

config.Site.storageSite = 'T2_PL_Swierk'
