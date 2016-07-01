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
#process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_miniAODv2'
process.GlobalTag.globaltag = '80X_dataRun2_Prompt_v8'

#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True), SkipEvent = cms.untracked.vstring('ProductNotFound') )
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.maxEvents = cms.untracked.PSet(
	input = cms.untracked.int32(-1)
)

#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(
        ## Single Electron
        '/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/150/00000/0A6284C7-D719-E611-93E6-02163E01421D.root',
        '/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/158/00000/06277EC1-181A-E611-870F-02163E0145E5.root',
        '/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/158/00000/0A7BD549-131A-E611-8287-02163E0134FC.root',
        '/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/158/00000/0C65F947-0F1A-E611-A1E6-02163E0144EA.root',
        '/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/158/00000/1CCC1100-0E1A-E611-98C7-02163E014332.root',
        '/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/158/00000/2209DFC5-191A-E611-9809-02163E014272.root',
        '/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/158/00000/26868D6F-4A1A-E611-8916-02163E011D33.root',
        '/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/158/00000/361F045C-111A-E611-AEF0-02163E01457C.root',
        '/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/158/00000/3A8562C4-2B1A-E611-96AC-02163E0141D2.root',
        '/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/158/00000/429D71B2-1D1A-E611-A5A9-02163E013926.root'
        
        ## Single Muon
        #'/store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/150/00000/34A57FB8-D819-E611-B0A4-02163E0144EE.root',
        #'/store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/158/00000/02D9C19F-571A-E611-AD8E-02163E013732.root',
        #'/store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/158/00000/18383F36-2E1A-E611-8C57-02163E014186.root',
        #'/store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/158/00000/1AF0711C-241A-E611-AC07-02163E0141E1.root',
        #'/store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/158/00000/2096A6AC-261A-E611-99BD-02163E01355E.root',
        #'/store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/158/00000/2226C1A1-571A-E611-B1E2-02163E013441.root',
        #'/store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/158/00000/422892A9-571A-E611-9BB8-02163E011D78.root',
        #'/store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/158/00000/66E2AFB2-571A-E611-BF2A-02163E013674.root',
        #'/store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/158/00000/701EA1A3-571A-E611-95AC-02163E011D78.root',
        #'/store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/158/00000/8A14EDAB-571A-E611-96F5-02163E0123DD.root'
        
        ## Double Electron
        #'/store/data/Run2016B/DoubleEG/MINIAOD/PromptReco-v2/000/273/150/00000/66C653C9-D719-E611-A715-02163E0142D9.root',
        #'/store/data/Run2016B/DoubleEG/MINIAOD/PromptReco-v2/000/273/158/00000/0227DB1C-E719-E611-872C-02163E0141F9.root',
        #'/store/data/Run2016B/DoubleEG/MINIAOD/PromptReco-v2/000/273/158/00000/1E5ABF54-E019-E611-AAED-02163E01293F.root',
        #'/store/data/Run2016B/DoubleEG/MINIAOD/PromptReco-v2/000/273/158/00000/2227B110-EE19-E611-9B2C-02163E0145A3.root',
        #'/store/data/Run2016B/DoubleEG/MINIAOD/PromptReco-v2/000/273/158/00000/2255BB59-EA19-E611-8305-02163E011D10.root',
        #'/store/data/Run2016B/DoubleEG/MINIAOD/PromptReco-v2/000/273/158/00000/26F9B42D-E319-E611-BBEB-02163E013665.root',
        #'/store/data/Run2016B/DoubleEG/MINIAOD/PromptReco-v2/000/273/158/00000/2881ECF1-EE19-E611-9B73-02163E0125FE.root',
        #'/store/data/Run2016B/DoubleEG/MINIAOD/PromptReco-v2/000/273/158/00000/2EEDB2BB-DD19-E611-8C26-02163E011E0E.root',
        #'/store/data/Run2016B/DoubleEG/MINIAOD/PromptReco-v2/000/273/158/00000/308C0A95-E819-E611-9413-02163E01386F.root',
        #'/store/data/Run2016B/DoubleEG/MINIAOD/PromptReco-v2/000/273/158/00000/32EBF0E5-E619-E611-9903-02163E014118.root'
        
        ## Electron Muon
        #'/store/data/Run2016B/MuonEG/MINIAOD/PromptReco-v2/000/273/150/00000/7EEDA70D-D819-E611-8491-02163E01349D.root',
        #'/store/data/Run2016B/MuonEG/MINIAOD/PromptReco-v2/000/273/158/00000/26281378-291A-E611-AE69-02163E011E9B.root',
        #'/store/data/Run2016B/MuonEG/MINIAOD/PromptReco-v2/000/273/158/00000/A2E8DB05-371A-E611-9B55-02163E01413C.root',
        #'/store/data/Run2016B/MuonEG/MINIAOD/PromptReco-v2/000/273/158/00000/A407F64D-221A-E611-8060-02163E011998.root',
        #'/store/data/Run2016B/MuonEG/MINIAOD/PromptReco-v2/000/273/158/00000/CC82FDA4-201A-E611-BC48-02163E013552.root',
        #'/store/data/Run2016B/MuonEG/MINIAOD/PromptReco-v2/000/273/290/00000/0A065771-161A-E611-8712-02163E0144F0.root',
        #'/store/data/Run2016B/MuonEG/MINIAOD/PromptReco-v2/000/273/291/00000/3A2AFF67-411A-E611-B4BB-02163E0139C8.root',
        #'/store/data/Run2016B/MuonEG/MINIAOD/PromptReco-v2/000/273/292/00000/C8EF3916-231A-E611-AC8E-02163E01411A.root',
        #'/store/data/Run2016B/MuonEG/MINIAOD/PromptReco-v2/000/273/294/00000/C2E47BF3-231A-E611-BC1D-02163E0137D0.root',
        #'/store/data/Run2016B/MuonEG/MINIAOD/PromptReco-v2/000/273/295/00000/A835EC00-241A-E611-A677-02163E0144FB.root'
        
        ## Double Muon
        #'/store/data/Run2016B/DoubleMuon/MINIAOD/PromptReco-v2/000/273/150/00000/680BED0F-D919-E611-85E6-02163E01424F.root',
        #'/store/data/Run2016B/DoubleMuon/MINIAOD/PromptReco-v2/000/273/158/00000/2C8772DF-F319-E611-AEC1-02163E014122.root',
        #'/store/data/Run2016B/DoubleMuon/MINIAOD/PromptReco-v2/000/273/158/00000/6AB163DF-121A-E611-9F00-02163E011CDC.root',
        #'/store/data/Run2016B/DoubleMuon/MINIAOD/PromptReco-v2/000/273/158/00000/80DEDA0B-131A-E611-8960-02163E011C3C.root',
        #'/store/data/Run2016B/DoubleMuon/MINIAOD/PromptReco-v2/000/273/158/00000/A6AC80E5-121A-E611-A689-02163E01439E.root',
        #'/store/data/Run2016B/DoubleMuon/MINIAOD/PromptReco-v2/000/273/158/00000/AE16A7ED-121A-E611-982A-02163E011AB8.root',
        #'/store/data/Run2016B/DoubleMuon/MINIAOD/PromptReco-v2/000/273/158/00000/B2738FB4-E919-E611-8FA9-02163E012338.root',
        #'/store/data/Run2016B/DoubleMuon/MINIAOD/PromptReco-v2/000/273/158/00000/C0DF4FD8-121A-E611-8E96-02163E011E3E.root',
        #'/store/data/Run2016B/DoubleMuon/MINIAOD/PromptReco-v2/000/273/158/00000/DE48F8C5-EB19-E611-A747-02163E014430.root',
        #'/store/data/Run2016B/DoubleMuon/MINIAOD/PromptReco-v2/000/273/158/00000/E25D967A-EE19-E611-92DB-02163E01380C.root'
        
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
