from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'omtf_run3_ZeroBias_Run2018D_t127'
#config.General.workArea = 'jobs_SM_Run2017E-ZMu-17Nov2017'
#config.General.workArea = 'jobs_JHT_2018D'
config.General.transferLogs = True 
config.General.transferOutputs = True 

config.section_("Data")
#config.Data.inputDataset = '/SingleMuon/Run2016H-ZMu-PromptReco-v2/RAW-RECO'
#config.Data.inputDataset = '/SingleMuon/Run2016H-PromptReco-v2/AOD'
#config.Data.inputDataset = '/ExpressPhysics/Run2017A-Express-v1/FEVT'
#config.Data.inputDataset = '/SingleMuon/Run2017F-ZMu-PromptReco-v1/RAW-RECO'
#config.Data.inputDataset = '/SingleMuon/Run2018D-ZMu-PromptReco-v2/RAW-RECO'
#config.Data.inputDataset = '/JetHT/Run2018D-JetHTJetPlusHOFilter-PromptReco-v2/RAW-RECO'
#config.Data.inputDataset = '/JetHT/Run2017F-JetHTJetPlusHOFilter-09May2018-v1/RAW-RECO'
#config.Data.inputDataset = '/SingleMuon/Run2017E-ZMu-17Nov2017-v1/RAW-RECO'
config.Data.inputDataset = '/ZeroBias/Run2018D-v1/RAW'


#config.Data.lumiMask='Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
config.Data.lumiMask='Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt'

#config.Data.runRange = '282000-283000'
#config.Data.runRange = '295606'

config.Data.useParent = False 
config.Data.inputDBS = 'global'
#config.Data.splitting = 'LumiBased'
config.Data.splitting = 'Automatic'
#config.Data.unitsPerJob = 200 #number of files per jobs
config.Data.totalUnits =  20000000 #number of event
#config.Data.outLFNDirBase = '/store/user/konec/test/'
config.Data.publication = False 
#config.Data.outputDatasetTag = 'CRAB3_tutorial_May2015_Data_analysis'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runMuonOverlap_run3_realData.py'
config.JobType.pyCfgParams = ['rate']
#config.JobType.disableAutomaticOutputCollection = True
#config.JobType.outputFiles = ['omtfTree.root', 'omtfHelper.root']

config.section_("Site")
#config.Site.whitelist = ['T3_CH_CERNCAF']
#config.Site.whitelist = ['T2_CH_CERN']
config.Site.storageSite = 'T2_PL_Swierk'
#config.Site.blacklist = ['T2_KR_*','T2_CN_*','T2_BR_*','T2_US_Florida','T2_US_UCSD']
