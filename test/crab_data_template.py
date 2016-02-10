from CRABClient.UserUtilities import config
config = config()

# General section
config.General.requestName = "@uiworkingdir@"
config.General.workArea = "tasks"
config.General.transferOutputs = True
config.General.transferLogs = True

# JobType section
config.JobType.pluginName = "Analysis"
config.JobType.psetName = "Extractor_MULTIJET_data.py"
config.JobType.pyCfgParams = ['globalTag=@globaltag@']
config.JobType.outputFiles = ['extracted.root']
config.JobType.allowUndistributedCMSSW = True

# Data section
config.Data.inputDataset = "@datasetname@"
config.Data.inputDBS = "@dbs_url@"
config.Data.outLFNDirBase = "@remote_dir@"
config.Data.splitting = "LumiBased"
config.Data.unitsPerJob = 15
config.Data.lumiMask = "@golden_json_file@"
# crab bug: remove ignoreLocality and whitelist ; see https://hypernews.cern.ch/HyperNews/CMS/get/computing-tools/958/2/1/1/1/1/1/1/1/1/1/1/1/1/1/1.html
#config.Data.ignoreLocality = True

# Site section
config.Site.storageSite = "T3_FR_IPNL"
#config.Site.whitelist = ['T0_CH_CERN_Disk', 'T1_IT_CNAF_Buffer', 'T1_IT_CNAF_Disk', 'T1_IT_CNAF_MSS', 'T2_BE_UCL', 'T2_CH_CERN', 'T2_DE_DESY', 'T2_FR_IPHC', 'T2_IT_Legnaro', 'T2_UK_London_IC']
#config.Site.whitelist = ['T2_BE_UCL', 'T2_CH_CERN', 'T2_DE_DESY', 'T2_FR_IPHC', 'T2_IT_Legnaro', 'T2_UK_London_IC']

