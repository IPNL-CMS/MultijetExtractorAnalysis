from CRABClient.UserUtilities import config
config = config()

# General section
config.General.requestName = "@uiworkingdir@"
config.General.workArea = "tasks"
config.General.transferOutputs = True
config.General.transferLogs = True

# JobType section
config.JobType.pluginName = "Analysis"
config.JobType.psetName = "Extractor_MULTIJET_MC.py"
config.JobType.pyCfgParams = []
config.JobType.outputFiles = ['extracted_mc.root']
config.JobType.allowUndistributedCMSSW = True

# Data section
config.Data.inputDataset = "@datasetname@"
config.Data.inputDBS = "@dbs_url@"
config.Data.outLFNDirBase = "@remote_dir@"
config.Data.splitting = "EventAwareLumiBased"
config.Data.unitsPerJob = 50000

# Site section
config.Site.storageSite = "T3_FR_IPNL"
