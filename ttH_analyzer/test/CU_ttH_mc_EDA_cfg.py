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

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
#process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")

process.load( "Configuration.StandardSequences.FrontierConditions_GlobalTag_cff" )
#process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_miniAODv2_v1'
process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v7'

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

process.maxEvents = cms.untracked.PSet(
	input = cms.untracked.int32(-1)
)

process.ak4PFCHSL1Fastjet = cms.ESProducer(
    'L1FastjetCorrectionESProducer',
    level       = cms.string('L1FastJet'),
    algorithm   = cms.string('AK4PFchs'),
    srcRho      = cms.InputTag( 'fixedGridRhoFastjetAll' )
    )

process.ak4PFchsL2Relative = ak4CaloL2Relative.clone( algorithm = 'AK4PFchs' )
process.ak4PFchsL3Absolute = ak4CaloL3Absolute.clone( algorithm = 'AK4PFchs' )

process.ak4PFchsL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
	'ak4PFCHSL1Fastjet', 
    'ak4PFchsL2Relative', 
    'ak4PFchsL3Absolute')
)

# HIP mitigted B-tagger

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(
  process,
  jetSource = cms.InputTag('slimmedJets'),
  jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None' ),  
  btagDiscriminators = ['pfCombinedInclusiveSecondaryVertexV2BJetTags'],
  runIVF=True,
  #btagPrefix = 'new' # optional, in case interested in accessing both the old and new discriminator values
)

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
# If you only want to re-correct and get the proper uncertainties
runMetCorAndUncFromMiniAOD(process,
  isData=False,
  #jecUncFile = "data/JEC/Spring16_25nsV6_MC_Uncertainty_AK4PFchs.txt"
)

process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")

#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(
        #htobb	
 	##reHLT
        ##'/store/mc/RunIISpring16MiniAODv2/ttHTobb_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/40000/0089CC67-6338-E611-947D-0025904C4E2A.root'
        #withHLT tranche 3
        '/store/mc/RunIISpring16MiniAODv2/ttHTobb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/MINIAODSIM/premix_withHLT_80X_mcRun2_asymptotic_v14-v1/100000/086ED46A-1E76-E611-ABE5-180373FF8456.root'	
	#htoNonbb
        #witHLT tranche 3
        #'/store/mc/RunIISpring16MiniAODv2/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/MINIAODSIM/premix_withHLT_80X_mcRun2_asymptotic_v14-v1/90000/00212379-8674-E611-956C-0CC47A4C8E22.root'
        #tt+jets 
        ##reHLT
        ##'/store/mc/RunIISpring16MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext3-v1/00000/0064B539-803A-E611-BDEA-002590D0B060.root'
        #withHLT tranche 3
        #semileptonic 
        #'/store/mc/RunIISpring16MiniAODv2/TTToSemilepton_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/MINIAODSIM/premix_withHLT_80X_mcRun2_asymptotic_v14-v1/100000/00B939E4-F982-E611-AB1C-F46D042E833B.root'
        #dileptonic
        #'/store/mc/RunIISpring16MiniAODv2/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/MINIAODSIM/premix_withHLT_80X_mcRun2_asymptotic_v14-v1/100000/005054EA-CB80-E611-9836-14187733AD89.root'
        )
)

# new electron MVA developed by the EGamma POG 
#process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff']
for idmod in my_id_modules:  
	setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

# load the analysis:
process.load("Analyzers.ttH_analyzer.ttHbb_cfi")

#ttHf categorization

# Setting input particle collections to be used by the tools
genJetCollection = 'ak4GenJetsCustom'
genParticleCollection = 'prunedGenParticles'
genJetInputParticleCollection = 'packedGenParticles'
    
## producing a subset of particles to be used for jet clustering
from RecoJets.Configuration.GenJetParticles_cff import genParticlesForJetsNoNu
process.genParticlesForJetsNoNu = genParticlesForJetsNoNu.clone(
	src = genJetInputParticleCollection
)
    
# Supplies PDG ID to real name resolution of MC particles
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
    
# Producing own jets for testing purposes
from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
process.ak4GenJetsCustom = ak4GenJets.clone(
        src = 'genParticlesForJetsNoNu',
        #    src = genJetInputParticleCollection,
        rParam = cms.double(0.4),
        jetAlgorithm = cms.string("AntiKt")
)
    
# Ghost particle collection used for Hadron-Jet association 
# MUST use proper input particle collection
from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartons
process.selectedHadronsAndPartons = selectedHadronsAndPartons.clone(
        particles = genParticleCollection
)
    
# Input particle collection for matching to gen jets (partons + leptons) 
# MUST use use proper input jet collection: the jets to which hadrons should be associated
# rParam and jetAlgorithm MUST match those used for jets to be associated with hadrons
# More details on the tool: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools#New_jet_flavour_definition

from PhysicsTools.JetMCAlgos.AK4PFJetsMCFlavourInfos_cfi import ak4JetFlavourInfos
process.genJetFlavourInfos = ak4JetFlavourInfos.clone(
        jets = genJetCollection,
        rParam = cms.double(0.4),
        jetAlgorithm = cms.string("AntiKt")
)
    
from PhysicsTools.JetMCAlgos.GenHFHadronMatcher_cff import *
# Plugin for analysing B hadrons
# MUST use the same particle collection as in selectedHadronsAndPartons
process.matchGenBHadron = matchGenBHadron.clone(
        genParticles = genParticleCollection,
        jetFlavourInfos = "genJetFlavourInfos"
)
    
# Plugin for analysing C hadrons
# MUST use the same particle collection as in selectedHadronsAndPartons
process.matchGenCHadron = matchGenCHadron.clone(
        genParticles = genParticleCollection,
        jetFlavourInfos = "genJetFlavourInfos"
)
    
## Producer for ttbar categorisation ID
# MUST use same genJetCollection as used for tools above
from PhysicsTools.JetMCAlgos.GenTtbarCategorizer_cfi import categorizeGenTtbar
process.categorizeGenTtbar = categorizeGenTtbar.clone(
        genJetPtMin = 20.,
        genJetAbsEtaMax = 2.4,
        genJets = genJetCollection,
)

process.TFileService = cms.Service("TFileService",
	fileName = cms.string('ttHbbNtuple.root')
)

process.p = cms.Path(
    #process.electronMVAValueMapProducer
    process.egmGsfElectronIDSequence	
    * process.genParticlesForJetsNoNu
    * process.ak4GenJetsCustom
    * process.selectedHadronsAndPartons
    * process.genJetFlavourInfos
    * process.matchGenBHadron
    * process.matchGenCHadron
    * process.categorizeGenTtbar
    * process.fullPatMetSequence
    * process.BadPFMuonFilter
    * process.BadChargedCandidateFilter
    * process.ttHbb
)
