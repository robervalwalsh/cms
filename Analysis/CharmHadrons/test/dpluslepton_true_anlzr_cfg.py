import FWCore.ParameterSet.Config as cms

process = cms.Process("CHARMPLOTS")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#process.GlobalTag.globaltag = "START36_V10::All"
process.GlobalTag.globaltag = "DESIGN_36_V10::All"

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.TFileService = cms.Service("TFileService", 
   fileName = cms.string("Dplus_TrueDplusLepton.root"),
   closeFileFast = cms.untracked.bool(True)
)

process.DplusLeptonTrueAnalyzer = cms.EDProducer('DplusLeptonTrueAnalyzer'
   , GeneratedParticles = cms.InputTag("genParticles")
   , PtMin = cms.double(3.0)
   , EtaRange = cms.vdouble(-2.5,2.5)
)


process.p = cms.Path(
   process.DplusLeptonTrueAnalyzer
)

process.out = cms.OutputModule('PoolOutputModule',
# Data
   fileName = cms.untracked.string("test.root"),
   SelectEvents = cms.untracked.PSet(
      SelectEvents = cms.vstring("p")
   ),
   outputCommands = cms.untracked.vstring(
      "keep *",
   )
)

#process.endPath = cms.EndPath(process.out)


readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
process.source = cms.Source ('PoolSource',
         fileNames = readFiles,
         secondaryFileNames = secFiles)

readFiles.extend( [
      '/store/user/rwalsh/MinBias_TuneD6T_7TeV-pythia6/Summer10-START36_V10_SP10-v1_DplusTrueFilter/ea5be3b5e51b9702e615335dd9650f27/Dplus_10_1_kcv.root',
      '/store/user/rwalsh/MinBias_TuneD6T_7TeV-pythia6/Summer10-START36_V10_SP10-v1_DplusTrueFilter/ea5be3b5e51b9702e615335dd9650f27/Dplus_11_1_oJq.root',
      '/store/user/rwalsh/MinBias_TuneD6T_7TeV-pythia6/Summer10-START36_V10_SP10-v1_DplusTrueFilter/ea5be3b5e51b9702e615335dd9650f27/Dplus_12_1_VzW.root',
      '/store/user/rwalsh/MinBias_TuneD6T_7TeV-pythia6/Summer10-START36_V10_SP10-v1_DplusTrueFilter/ea5be3b5e51b9702e615335dd9650f27/Dplus_13_1_Lcq.root',
      '/store/user/rwalsh/MinBias_TuneD6T_7TeV-pythia6/Summer10-START36_V10_SP10-v1_DplusTrueFilter/ea5be3b5e51b9702e615335dd9650f27/Dplus_14_1_Hjg.root',
      '/store/user/rwalsh/MinBias_TuneD6T_7TeV-pythia6/Summer10-START36_V10_SP10-v1_DplusTrueFilter/ea5be3b5e51b9702e615335dd9650f27/Dplus_15_1_EBs.root'
] );



secFiles.extend( [
               ] )
