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

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True), SkipEvent = cms.untracked.vstring('ProductNotFound') )

process.maxEvents = cms.untracked.PSet(
	input = cms.untracked.int32(-1)
)

#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(
        #htobb
        #'/store/mc/RunIISpring16MiniAODv2/ttHTobb_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v2/90000/0204C78B-9434-E611-B3B0-02163E014C6F.root'
        # tt+jet 
        '/store/mc/RunIISpring16MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext3-v1/00000/000B9244-4B27-E611-91D2-7845C4FC3C6B.root'
        )
)

#ttHf categorization
#process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
#process.load("PhysicsTools.JetMCAlgos.GenHFHadronMatcher_cfi")
# new electron MVA developed by the EGamma POG 
process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")
# load the analysis:
process.load("Analyzers.ttH_analyzer.ttHbb_cfi")

#process.matchGenHFHadron.genParticles = cms.InputTag('prunedGenParticles')
#genParticleCollection = 'prunedGenParticles'
#genJetCollection = 'slimmedGenJets'

#from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartons
#process.selectedHadronsAndPartons = selectedHadronsAndPartons.clone(
#    particles = genParticleCollection
#)

#from PhysicsTools.JetMCAlgos.AK4PFJetsMCFlavourInfos_cfi import ak4JetFlavourInfos
#process.genJetFlavourInfos = ak4JetFlavourInfos.clone(
#    jets = genJetCollection
#)

# Plugin for analysing B hadrons
# MUST use the same particle collection as in selectedHadronsAndPartons
#from PhysicsTools.JetMCAlgos.GenHFHadronMatcher_cff import matchGenBHadron
#process.matchGenBHadron = matchGenBHadron.clone(
#    genParticles = genParticleCollection,
#    jetFlavourInfos = "genJetFlavourInfos"
#)

# Plugin for analysing C hadrons
# MUST use the same particle collection as in selectedHadronsAndPartons
#from PhysicsTools.JetMCAlgos.GenHFHadronMatcher_cff import matchGenCHadron
#process.matchGenCHadron = matchGenCHadron.clone(
#    genParticles = genParticleCollection,
#    jetFlavourInfos = "genJetFlavourInfos"
#)

process.TFileService = cms.Service("TFileService",
	fileName = cms.string('ttHbbNtuple.root')
)

process.p = cms.Path(
    process.electronMVAValueMapProducer
    #* process.selectedHadronsAndPartons
    #* process.genJetFlavourInfos
    #* process.matchGenCHadron
    #* process.matchGenBHadron
    * process.ttHbb
)
