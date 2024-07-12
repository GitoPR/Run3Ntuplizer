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
                                'root://cms-xrd-global.cern.ch://store/data/Run2024F/ZeroBias/MINIAOD/PromptReco-v1/000/382/725/00000/25e880d8-5f24-45cf-8cf5-2dff21a81c88.root'
                            ),
                            secondaryFileNames = cms.untracked.vstring(
                                'root://cms-xrd-global.cern.ch://store/data/Run2024F/ZeroBias/RAW/v1/000/382/725/00000/039feb9f-dc91-48c0-bcc4-022ec918e558.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2024F/ZeroBias/RAW/v1/000/382/725/00000/086eb686-6d06-4289-8ea2-0191721908b1.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2024F/ZeroBias/RAW/v1/000/382/725/00000/10142e41-f20d-40d1-9546-0ad7a99f0305.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2024F/ZeroBias/RAW/v1/000/382/725/00000/1ca5be12-bbff-4a2d-923f-7f3850cb2a14.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2024F/ZeroBias/RAW/v1/000/382/725/00000/21dbcd9b-e8d0-49bd-b41c-dee21dfe9e90.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2024F/ZeroBias/RAW/v1/000/382/725/00000/2d9d18fe-15bc-4b1e-a192-156e71232fb6.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2024F/ZeroBias/RAW/v1/000/382/725/00000/301caaf8-fba4-4394-b60b-3e42e98cadaa.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2024F/ZeroBias/RAW/v1/000/382/725/00000/371ab874-20ed-4b5d-b743-d63eb3358fe6.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2024F/ZeroBias/RAW/v1/000/382/725/00000/3ae63922-2f44-4e12-ac85-c533f52d4b37.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2024F/ZeroBias/RAW/v1/000/382/725/00000/57f090b4-b045-40c1-85d4-c9bfb4049a09.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2024F/ZeroBias/RAW/v1/000/382/725/00000/5c4560fb-ef2a-4729-bcab-7f7fc0b58f2d.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2024F/ZeroBias/RAW/v1/000/382/725/00000/6d801b8a-2165-415f-85f4-2dbaebaf7bef.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2024F/ZeroBias/RAW/v1/000/382/725/00000/7bf39f26-5d71-4062-a667-2226108bec43.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2024F/ZeroBias/RAW/v1/000/382/725/00000/90da2f6a-8bff-4822-aa54-373df6f5ff6b.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2024F/ZeroBias/RAW/v1/000/382/725/00000/930e53b9-590f-499b-bed5-b8af5971162a.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2024F/ZeroBias/RAW/v1/000/382/725/00000/a5dc515f-ab1c-4e2b-9a29-33853441faba.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2024F/ZeroBias/RAW/v1/000/382/725/00000/a98610c2-ed72-4c47-9fe8-b7eed8facad7.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2024F/ZeroBias/RAW/v1/000/382/725/00000/aa8fdb2a-193e-474b-aa31-28f1f8a51eec.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2024F/ZeroBias/RAW/v1/000/382/725/00000/b450e2e9-4b45-45f6-a9c8-68352973a9c4.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2024F/ZeroBias/RAW/v1/000/382/725/00000/b78de906-a9da-404e-9a4d-9cf6fd33031a.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2024F/ZeroBias/RAW/v1/000/382/725/00000/bde23107-8d32-4abe-bfbe-c88fbea27fe3.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2024F/ZeroBias/RAW/v1/000/382/725/00000/bf80e503-4e1b-45fc-99f1-bd20ff09bbdd.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2024F/ZeroBias/RAW/v1/000/382/725/00000/c2a185b9-4e0e-409c-9b39-be64e54b416e.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2024F/ZeroBias/RAW/v1/000/382/725/00000/caefad16-87dd-4311-93ab-03d427fc801e.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2024F/ZeroBias/RAW/v1/000/382/725/00000/cfa1a3b9-b91a-4cfd-a113-0567287a27d0.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2024F/ZeroBias/RAW/v1/000/382/725/00000/dff32c2c-8d17-4431-ae47-d566339adfe3.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2024F/ZeroBias/RAW/v1/000/382/725/00000/e80b58c6-8cdb-469b-bedc-8a3e68f48c7c.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2024F/ZeroBias/RAW/v1/000/382/725/00000/ec076523-4201-4495-9f58-6804b8b9d578.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2024F/ZeroBias/RAW/v1/000/382/725/00000/f316ba1c-3547-4594-80bf-22fb0e944798.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2024F/ZeroBias/RAW/v1/000/382/725/00000/f3e27b91-5849-477d-b451-ae9f78442e53.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2024F/ZeroBias/RAW/v1/000/382/725/00000/fde673cb-6980-4e13-a488-47139a1db039.root',
                                'root://cms-xrd-global.cern.ch://store/data/Run2024F/ZeroBias/RAW/v1/000/382/725/00000/57f090b4-b045-40c1-85d4-c9bfb4049a09.root'
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
