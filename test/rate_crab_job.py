import sys
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'Run3_Hbb_rate'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.allowUndistributedCMSSW = True
config.JobType.maxJobRuntimeMin = 300
config.JobType.numCores = 8
config.JobType.maxMemoryMB = 9000
config.JobType.psetName = '/afs/hep.wisc.edu/home/jorgeeh/CMSSW_14_0_14/src/L1Trigger/Run3Ntuplizer/test/testL1TCaloSummary-ZeroBias.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['l1TNtuple-ZeroBias.root']

config.Data.inputDataset = '/ZeroBias/Run2024F-PromptReco-v1/MINIAOD' 
config.Data.secondaryInputDataset = '/ZeroBias/Run2024F-v1/RAW'
config.Data.lumiMask = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions24/Cert_Collisions2024_378981_384052_Golden.json'


config.Data.unitsPerJob = 1
config.Data.splitting = 'FileBased'
config.Data.outLFNDirBase = '/store/user/joeherna'


config.Site.storageSite = 'T2_US_Wisconsin'
