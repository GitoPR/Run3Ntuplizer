import sys
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'l1TNtuple-ggHBB'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.allowUndistributedCMSSW = True
config.JobType.maxJobRuntimeMin = 300
config.JobType.numCores = 8
config.JobType.maxMemoryMB = 9000
config.JobType.psetName = '/afs/hep.wisc.edu/home/jorgeeh/CMSSW_14_0_14/src/L1Trigger/Run3Ntuplizer/test/testL1TCaloSummary-ggHBB.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['l1TNtuple-ggHBB.root']

config.Data.inputDataset = '/Boosted_ggHbb_M-125_Pt-250_GENSIM/joeherna-Boosted_ggHbb_M-125_Pt-250_MINIAOD-42ad603b91d99d3cc77a3641b840e853/USER'
config.Data.secondaryInputDataset =  '/Boosted_ggHbb_M-125_Pt-250_GENSIM/joeherna-Boosted_ggHbb_M-125_Pt-250_DIGIRAW-905680e847f46ad1a422d47384efadca/USER'
config.Data.inputDBS = 'phys03'


config.Data.unitsPerJob = 1
config.Data.splitting = 'FileBased'
config.Data.outLFNDirBase = '/store/user/joeherna'


config.Site.storageSite = 'T2_US_Wisconsin'

