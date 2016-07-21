from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'ttH_Analyzer_hbb'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'test/CU_ttH_mc_EDA_cfg.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.inputFiles = ['data']
config.JobType.outputFiles = ['Cornell_tth_80X.csv']

config.Data.userInputFiles = ['/store/mc/RunIISpring16MiniAODv2/ttHTobb_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/40000/0089CC67-6338-E611-947D-0025904C4E2A.root']
config.Data.outputPrimaryDataset = 'ttH_Analyzer_hbb'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/abdatta/' 
config.Data.publication = True
config.Data.outputDatasetTag = 'ttH_Analyzer_hbb'

config.Site.storageSite = 'T3_US_FNALLPC'
#config.Site.storageSite = 'T3_US_Cornell'
#config.Site.whitelist = ['T3_US_Cornell']
