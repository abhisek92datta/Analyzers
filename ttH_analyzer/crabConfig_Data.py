from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'ttH_Analyzer_SingleElectron_Run2016B'
#config.General.requestName = 'ttH_Analyzer_SingleMuon_Run2016B'
#config.General.requestName = 'ttH_Analyzer_DoubleElectron_Run2016B'
#config.General.requestName = 'ttH_Analyzer_ElectronMuon_Run2016B'
#config.General.requestName = 'ttH_Analyzer_DoubleMuon_Run2016B'

#config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'test/CU_ttH_data_EDA_cfg.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.inputFiles = ['data']
config.JobType.maxJobRuntimeMin = 3600
#config.JobType.outputFiles = ['Cornell_tth_80X.csv']

config.section_("Data")

#Single Electron
#for B
config.Data.inputDataset = '/SingleElectron/Run2016B-03Feb2017_ver2-v2/MINIAOD'
#for C,D,E,F,G
#config.Data.inputDataset = '/SingleElectron/Run2016C-03Feb2017-v1/MINIAOD'
#for H
#config.Data.inputDataset = '/SingleElectron/Run2016H-03Feb2017_ver2-v1/MINIAOD'
#config.Data.inputDataset = '/SingleElectron/Run2016H-03Feb2017_ver3-v1/MINIAOD'

#Single Muon
#for B
#config.Data.inputDataset = '/SingleMuon/Run2016B-03Feb2017_ver2-v2/MINIAOD'
#for C,D,E,F,G
#config.Data.inputDataset = '/SingleMuon/Run2016C-03Feb2017-v1/MINIAOD'
#for H
#config.Data.inputDataset = '/SingleMuon/Run2016H-03Feb2017_ver2-v1/MINIAOD'
#config.Data.inputDataset = '/SingleMuon/Run2016H-03Feb2017_ver3-v1/MINIAOD'

#Double Electron
#for B
#config.Data.inputDataset = '/DoubleEG/Run2016B-03Feb2017_ver2-v2/MINIAOD'
#for C,D,E,F,G
#config.Data.inputDataset = '/DoubleEG/Run2016C-03Feb2017-v1/MINIAOD'
#for H
#config.Data.inputDataset = '/DoubleEG/Run2016H-03Feb2017_ver2-v1/MINIAOD'
#config.Data.inputDataset = '/DoubleEG/Run2016H-03Feb2017_ver3-v1/MINIAOD'

#Electron Muon
#for B
#config.Data.inputDataset = '/MuonEG/Run2016B-03Feb2017_ver2-v2/MINIAOD'
#for C,D,E,F,G
#config.Data.inputDataset = '/MuonEG/Run2016C-03Feb2017-v1/MINIAOD'
#for H
#config.Data.inputDataset = '/MuonEG/Run2016H-03Feb2017_ver2-v1/MINIAOD'
#config.Data.inputDataset = '/MuonEG/Run2016H-03Feb2017_ver3-v1/MINIAOD'

#Double Muon
#for B
#config.Data.inputDataset = '/DoubleMuon/Run2016B-03Feb2017_ver2-v2/MINIAOD'
#for C,D,E,F,G
#config.Data.inputDataset = '/DoubleMuon/Run2016C-03Feb2017-v1/MINIAOD'
#for H
#config.Data.inputDataset = '/DoubleMuon/Run2016H-03Feb2017_ver2-v1/MINIAOD'
#config.Data.inputDataset = '/DoubleMuon/Run2016H-03Feb2017_ver3-v1/MINIAOD'

config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/'
config.Data.allowNonValidInputDataset = True
config.Data.splitting = 'LumiBased'
config.Data.lumiMask = 'data/JSON/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.Data.unitsPerJob = 20
config.Data.publication = True
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/'
config.Data.outputDatasetTag = 'ttH_Analyzer_SingleElectron_Run2016B'
#config.Data.outputDatasetTag = 'ttH_Analyzer_SingleMuon_Run2016B'
#config.Data.outputDatasetTag = 'ttH_Analyzer_DoubleElectron_Run2016B'
#config.Data.outputDatasetTag = 'ttH_Analyzer_ElectronMuon_Run2016B'
#config.Data.outputDatasetTag = 'ttH_Analyzer_DoubleMuon_Run2016B'
config.Data.outLFNDirBase = '/store/user/abdatta/Trigger_Analysis/'
config.Data.ignoreLocality = True

config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'
