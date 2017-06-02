from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'ttH_Analyzer_tthbb'
#config.General.requestName = 'ttH_Analyzer_tthnonbb'
#config.General.requestName = 'ttH_Analyzer_ttjets_ol'
#config.General.requestName = 'ttH_Analyzer_ttjets_pp'
#config.General.requestName = 'ttH_Analyzer_ttjets_pp_sl'
#config.General.requestName = 'ttH_Analyzer_ttjets_pp_dl'
#config.General.requestName = 'ttH_Analyzer_ttjets_pp_sl_hf'
#config.General.requestName = 'ttH_Analyzer_ttjets_pp_dl_hf'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'test/CU_ttH_mc_EDA_cfg.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.inputFiles = ['data']
config.JobType.maxJobRuntimeMin = 3600
#config.JobType.outputFiles = ['Cornell_tth_80X.csv']

config.section_("Data")
config.Data.inputDataset = '/ttHTobb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
#config.Data.inputDataset = '/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
#config.Data.inputDataset = '/ttbb_4FS_OpenLoops_13TeV-sherpa/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
#config.Data.inputDataset = '/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
#config.Data.inputDataset = '/TTToSemilepton_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
#config.Data.inputDataset = '/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
#config.Data.inputDataset = '/TTToSemilepton_ttbbFilter_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_TTbbWithttHFGenFilter_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
#config.Data.inputDataset = '/TTTo2L2Nu_ttbbFilter_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_TTbbWithttHFGenFilter_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'

config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/'
config.Data.allowNonValidInputDataset = True
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.publication = True
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/'
config.Data.outputDatasetTag = 'ttH_Analyzer_tthbb'
#config.Data.outputDatasetTag = 'ttH_Analyzer_tthnonbb'
#config.Data.outputDatasetTag = 'ttH_Analyzer_ttjets_ol'
#config.Data.outputDatasetTag = 'ttH_Analyzer_ttjets_pp'
#config.Data.outputDatasetTag = 'ttH_Analyzer_ttjets_pp_sl'
#config.Data.outputDatasetTag = 'ttH_Analyzer_ttjets_pp_dl'
#config.Data.outputDatasetTag = 'ttH_Analyzer_ttjets_pp_sl_hf'
#config.Data.outputDatasetTag = 'ttH_Analyzer_ttjets_pp_dl_hf'
config.Data.outLFNDirBase = '/store/user/abdatta/ttH_Analysis/'
config.Data.ignoreLocality = True

config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'
