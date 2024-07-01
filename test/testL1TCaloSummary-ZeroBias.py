import os
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_2023_cff import Run3_2023
process = cms.Process("L1TCaloSummaryTest", Run3_2023)

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag

process.GlobalTag = GlobalTag(process.GlobalTag, '130X_dataRun3_Prompt_v4', '')

process.load('L1Trigger.Configuration.SimL1Emulator_cff')
process.load('L1Trigger.Configuration.CaloTriggerPrimitives_cff')
process.load('EventFilter.L1TXRawToDigi.caloLayer1Stage2Digis_cfi')
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.load('CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi')

process.raw2digi_step = cms.Path(process.RawToDigi)
process.endjob_step = cms.EndPath(process.endOfProcess)

process.schedule = cms.Schedule(process.raw2digi_step, process.endjob_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

from L1Trigger.Configuration.customiseReEmul import L1TReEmulFromRAW
process = L1TReEmulFromRAW(process)


process.load("L1Trigger.Run3Ntuplizer.l1BoostedJetStudies_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                'root://cms-xrd-global.cern.ch//store/data/Run2023C/ZeroBias/MINIAOD/PromptReco-v1/000/367/229/00000/9dce1b71-bb1a-43e7-9855-b7f754e26c0b.root'
                            ),
                            secondaryFileNames = cms.untracked.vstring(
                                'root://cms-xrd-global.cern.ch://store/data/Run2023C/ZeroBias/RAW/v1/000/367/229/00000/153e5368-5613-4c28-abb4-9d2ce62f7100.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2023C/ZeroBias/RAW/v1/000/367/229/00000/271c6501-3d4b-4ee3-b9bc-b61bfe92eb3b.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2023C/ZeroBias/RAW/v1/000/367/229/00000/36c68528-7e81-413a-91a4-c6328ce86b1b.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2023C/ZeroBias/RAW/v1/000/367/229/00000/bba5dd27-5cf3-4bfa-9a68-781b821c3ff7.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2023C/ZeroBias/RAW/v1/000/367/229/00000/d5754704-c739-4511-a6a2-3ea707f69771.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2023C/ZeroBias/RAW/v1/000/367/229/00000/e3e15114-bf4b-4c96-89f8-c9f8fa329781.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2023C/ZeroBias/RAW/v1/000/367/229/00000/f852ae46-5926-4967-ad65-d13fd523c3ef.root',
                            )
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step2 nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)


process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("l1TFullEvent.root"),
    outputCommands = cms.untracked.vstring('keep *')
    #outputCommands = cms.untracked.vstring('drop *') #'keep *_*_*_L1TCaloSummaryTest')
    #outputCommands = cms.untracked.vstring('drop *', 'keep *_l1tCaloLayer1Digis_*_*, keep *_*_*_L1TCaloSummaryTest' )
)


#Output
process.TFileService = cms.Service(
	"TFileService",
	fileName = cms.string("l1TNtuple-ZeroBias.root")
)

process.p = cms.Path(process.l1tCaloLayer1Digis*process.simCaloStage2Layer1Digis*process.l1NtupleProducer)
process.schedule.append(process.p)

process.e = cms.EndPath(process.out)
#process.schedule.append(process.e)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)

# Multi-threading
process.options.numberOfThreads=cms.untracked.uint32(8)
process.options.numberOfStreams=cms.untracked.uint32(0)

# End adding early deletion

#dump_file = open('dump.py','w')
#dump_file.write(process.dumpPython())
