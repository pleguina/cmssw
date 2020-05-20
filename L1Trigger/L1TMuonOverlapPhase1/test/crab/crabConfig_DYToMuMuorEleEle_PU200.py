from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'muCorr_MC_analysis_DYToMuMuorEleEle_v1_t18'
#config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
#config.JobType.psetName = 'runMuonOverlapTTMergerAnalyzerCrab.py'
config.JobType.psetName = 'runMuonCorrelatorTTTracksAnanlyzer.py'
config.JobType.pyCfgParams = ['efficiency']

config.Data.inputDataset = '/DYToMuMuorEleEle_M-20_14TeV_pythia8/PhaseIITDRSpring19DR-PU200_106X_upgrade2023_realistic_v3-v1/GEN-SIM-DIGI-RAW'
#'/SingleMu_FlatPt-2to100/PhaseIIFall17D-L1TPU200_93X_upgrade2023_realistic_v5-v1/GEN-SIM-DIGI-RAW'
#config.Data.inputDataset = '/SingleMu_FlatPt-2to100/PhaseIIFall17D-L1TnoPU_93X_upgrade2023_realistic_v5-v1/GEN-SIM-DIGI-RAW'

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = 'CRAB3_muCorr_MC_analysis_DYToMuMuorEleEle_v1_t18'
config.Data.totalUnits = 155 #321

config.Site.storageSite = 'T2_PL_Swierk'
