from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'muCorr_MC_analysis_DYToLL_mtd5_v1_t11_1'
#config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
#config.JobType.psetName = 'runMuonOverlapTTMergerAnalyzerCrab.py'
config.JobType.psetName = 'runMuonCorrelatorTTTracksAnanlyzer_mtd5.py'
config.JobType.pyCfgParams = ['efficiency']

config.Data.inputDataset = '/DYToLL_M-50_14TeV_TuneCP5_pythia8/PhaseIIMTDTDRAutumn18DR-PU200_103X_upgrade2023_realistic_v2-v2/FEVT'
#'/SingleMu_FlatPt-2to100/PhaseIIFall17D-L1TPU200_93X_upgrade2023_realistic_v5-v1/GEN-SIM-DIGI-RAW'
#config.Data.inputDataset = '/SingleMu_FlatPt-2to100/PhaseIIFall17D-L1TnoPU_93X_upgrade2023_realistic_v5-v1/GEN-SIM-DIGI-RAW'

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = 'CRAB3_muCorr_MC_analysis_DYToLL_mtd5_v1_t11_1'
config.Data.totalUnits = 800 #321 total files is 3952 

config.Site.storageSite = 'T2_PL_Swierk'
