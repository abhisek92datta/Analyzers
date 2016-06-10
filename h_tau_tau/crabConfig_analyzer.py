from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'ttH_Analyzer_ntuple_tau_tau'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ttH_analyzer/test/CU_ttH_EDA_cfg.py'

config.Data.userInputFiles = ['/store/mc/RunIIFall15MiniAODv2/ttHJetToTT_M125_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v3/60000/0C6DA13E-38C8-E511-8F6E-00259055220A.root']
config.Data.outputPrimaryDataset = 'ttH_Analyzer'
#config.Data.inputDataset = '/store/mc/RunIIFall15MiniAODv2/ttHJetToTT_M125_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v3/60000/0C6DA13E-38C8-E511-8F6E-00259055220A.root'
#config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/abdatta/' 
config.Data.publication = True
config.Data.outputDatasetTag = 'ttH_Analyzer_ntuple_tau_tau'

config.Site.storageSite = 'T3_US_Cornell'
config.Site.whitelist = ['T3_US_Cornell']
