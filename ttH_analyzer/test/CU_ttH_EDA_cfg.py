import FWCore.ParameterSet.Config as cms

from RecoJets.Configuration.RecoJets_cff import *
from RecoJets.Configuration.RecoPFJets_cff import *
from JetMETCorrections.Configuration.JetCorrectionProducersAllAlgos_cff import *
from JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff import *
from JetMETCorrections.Configuration.JetCorrectionServices_cff import *

process = cms.Process("MAOD")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load( "Configuration.StandardSequences.FrontierConditions_GlobalTag_cff" )
#process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_v12'
#process.GlobalTag.globaltag = '74X_mcRun2_asymptotic_v4'
process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_miniAODv2'

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet(
	input = cms.untracked.int32(-1)
)

#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(
        #htobb
        '/store/mc/RunIISpring16MiniAODv2/ttHTobb_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v2/90000/0204C78B-9434-E611-B3B0-02163E014C6F.root'
        #'/store/mc/RunIISpring15MiniAODv2/ttHTobb_M125_13TeV_powheg_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/30000/DA1B6FD6-C46D-E511-9C7B-00A0D1EE29B8.root'
        # tt+jet 
        #'/store/mc/RunIISpring15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/0EE7E064-BE6F-E511-BB41-E4115BB4C4BC.root'
       	)
)

# new electron MVA developed by the EGamma POG 
process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")
# load the analysis:
process.load("Analyzers.ttH_analyzer.ttHbb_cfi")
#ttHf categorization
process.load("JetMCAlgos.GenHFHadronMatcher_cfi");
    
process.TFileService = cms.Service("TFileService",
	fileName = cms.string('ttHbbNtuple.root')
)

process.p = cms.Path(
    process.electronMVAValueMapProducer
    * process.matchGenHFHadron
    * process.ttHbb
)
