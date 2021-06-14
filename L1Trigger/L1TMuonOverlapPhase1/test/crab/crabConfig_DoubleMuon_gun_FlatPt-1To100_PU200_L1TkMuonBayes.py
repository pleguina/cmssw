from CRABClient.UserUtilities import config #, getUsernameFromSiteDB
config = config()

config.General.requestName = 'L1TkMuonBayes_MC_analysis_DoubleMuon_gun_FlatPt-1To100_PU200_t114'
#config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'

config.JobType.psetName = 'runL1TkMuonBayesAnanlyzer.py'


config.JobType.pyCfgParams = ['efficiency']

config.Data.inputDataset = '/DoubleMuon_gun_FlatPt-1To100/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v3/GEN-SIM-DIGI-RAW' 
#'/SingleMu_FlatPt-2to100/PhaseIIFall17D-L1TPU200_93X_upgrade2023_realistic_v5-v1/GEN-SIM-DIGI-RAW'
#config.Data.inputDataset = '/SingleMu_FlatPt-2to100/PhaseIIFall17D-L1TnoPU_93X_upgrade2023_realistic_v5-v1/GEN-SIM-DIGI-RAW'

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
#config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = 'CRAB3_L1TkMuonBayes_MC_analysis_DoubleMuon_gun_FlatPt-1To100_PU200_t114'
config.Data.totalUnits = 1118 #1118
config.Data.ignoreLocality = False

config.section_("Debug")
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']

#config.Site.storageSite = 'T2_PL_Swierk'
config.Site.storageSite = 'T2_CH_CERN'

