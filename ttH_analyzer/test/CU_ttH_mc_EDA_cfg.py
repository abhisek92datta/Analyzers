import os
import FWCore.ParameterSet.Config as cms
import PhysicsTools.PythonAnalysis.LumiList as LumiList
from FWCore.ParameterSet.VarParsing import VarParsing

from RecoJets.Configuration.RecoJets_cff import *
from RecoJets.Configuration.RecoPFJets_cff import *
from JetMETCorrections.Configuration.JetCorrectionProducersAllAlgos_cff import *
from JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff import *
from JetMETCorrections.Configuration.JetCorrectionServices_cff import *

#
# options
#
options = VarParsing("python")

# add custom options
options.register("realData",
                 False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "input dataset contains real data"
                 )
options.register("dataEra",
                 "",
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "the era of the data taking period, e.g. '2016B', empty for MC"
                 )
options.register("deterministicSeeds",
                 True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "create collections with deterministic seeds"
                 )
options.register("electronRegression",
                 "GT",
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "'GT' or an absolute path to a sqlite file for electron energy regression"
                 )
options.register("electronSmearing",
                 "Moriond17_23Jan",
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "correction type for electron energy smearing"
                 )
options.register("recorrectMET",
                 True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "recorrect MET using latest JES and e/g corrections"
                 )
options.register("updatePUJetId",
                 True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "update the PUJetId values"
                 )
options.register("isTtbar",
                 False, # set to True for all ttbar datasets
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "creates the ttbar gen id and performs ttbar heavy flavour tagging"
                 )
options.parseArguments()

#
# collection placeholders
#
electronCollection = cms.InputTag("slimmedElectrons", "", "PAT")
muonCollection     = cms.InputTag("slimmedMuons", "", "PAT")
tauCollection      = cms.InputTag("slimmedTaus", "", "PAT")
photonCollection   = cms.InputTag("slimmedPhotons", "", "PAT")
METCollection      = cms.InputTag("slimmedMETs", "", "PAT")
jetCollection      = cms.InputTag("slimmedJets", "", "PAT")


process = cms.Process("MAOD")

#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    #htobb
    #tranch IV
    '/store/mc/RunIISummer16MiniAODv2/ttHTobb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/44949CF4-96C6-E611-B9A0-0025905A6122.root'
    #tt+jets
    #tranche IV
    #'/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/0693E0E7-97BE-E611-B32F-0CC47A78A3D8.root'
    )
)


# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
#process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")

# Supplies PDG ID to real name resolution of MC particles
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.load( "Configuration.StandardSequences.FrontierConditions_GlobalTag_cff" )
process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

process.maxEvents = cms.untracked.PSet(
	input = cms.untracked.int32(100)
)

seq = cms.Sequence()

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


#
# deterministic seed producer
#


if options.deterministicSeeds:
    process.load("PhysicsTools.PatUtils.deterministicSeeds_cfi")
    process.deterministicSeeds.produceCollections = cms.bool(True)
    process.deterministicSeeds.produceValueMaps   = cms.bool(False)
    process.deterministicSeeds.electronCollection = electronCollection
    process.deterministicSeeds.muonCollection     = muonCollection
    process.deterministicSeeds.tauCollection      = tauCollection
    process.deterministicSeeds.photonCollection   = photonCollection
    process.deterministicSeeds.jetCollection      = jetCollection
    process.deterministicSeeds.METCollection      = METCollection

    # overwrite output collections
    electronCollection = cms.InputTag("deterministicSeeds", "electronsWithSeed", process.name_())
    muonCollection     = cms.InputTag("deterministicSeeds", "muonsWithSeed", process.name_())
    tauCollection      = cms.InputTag("deterministicSeeds", "tausWithSeed", process.name_())
    photonCollection   = cms.InputTag("deterministicSeeds", "photonsWithSeed", process.name_())
    jetCollection      = cms.InputTag("deterministicSeeds", "jetsWithSeed", process.name_())
    METCollection      = cms.InputTag("deterministicSeeds", "METsWithSeed", process.name_())


#
# electron energy regression
#


if options.electronRegression:
    if options.electronRegression == "GT":
        from EgammaAnalysis.ElectronTools.regressionWeights_cfi import regressionWeights
        process = regressionWeights(process)
    else:
        from EgammaAnalysis.ElectronTools.regressionWeights_local_cfi import GBRDWrapperRcd
        GBRDWrapperRcd.connect = cms.string("sqlite_file:" + options.electronRegression)
        process.regressions = GBRDWrapperRcd
        process.regressions.DumpStat = cms.untracked.bool(False)
        process.es_prefer_regressions = cms.ESPrefer("PoolDBESSource", "regressions")
    process.load("EgammaAnalysis.ElectronTools.regressionApplication_cff")
    seq += process.regressionApplication

    # set the electron and photon sources
    process.slimmedElectrons.src = electronCollection
    process.slimmedPhotons.src = photonCollection

    # overwrite output collections
    electronCollection = cms.InputTag("slimmedElectrons", "", process.name_())
    photonCollection = cms.InputTag("slimmedPhotons", "", process.name_())



#
# electron energy smearing
#


# only allowed when regression was used
if options.electronSmearing and options.electronRegression:
    # the smearing procedure requires a preselection
    process.selectedElectrons = cms.EDFilter("PATElectronSelector",
        src = electronCollection,
        cut = cms.string("pt>5 && abs(superCluster.eta)<2.5")
    )
    electronCollection = cms.InputTag("selectedElectrons", "", process.name_())

    # setup the smearing
    process.load("EgammaAnalysis.ElectronTools.calibratedPatElectronsRun2_cfi")
    from EgammaAnalysis.ElectronTools.calibratedPatElectronsRun2_cfi import files
    process.calibratedPatElectrons.isMC           = cms.bool(not options.realData)
    process.calibratedPatElectrons.correctionFile = cms.string(files[options.electronSmearing])
    process.calibratedPatElectrons.electrons      = electronCollection
    seq += process.calibratedPatElectrons

    # use our deterministic seeds or a random generator service
    if options.deterministicSeeds:
       process.calibratedPatElectrons.seedUserInt = process.deterministicSeeds.seedUserInt
    else:
       process.load("Configuration.StandardSequences.Services_cff")
       process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
            calibratedPatElectrons = cms.PSet(
                 initialSeed = cms.untracked.uint32(81),
                 engineName  = cms.untracked.string("TRandom3")
            )
       )
                                                         
    # overwrite output collections
    electronCollection = cms.InputTag("calibratedPatElectrons", "", process.name_())



#
# electron VIDs
#

from PhysicsTools.SelectorUtils.tools.vid_id_tools import DataFormat, \
    switchOnVIDElectronIdProducer, setupAllVIDIdsInModule, setupVIDElectronSelection

eleVIDModules = [
                 "RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff",
                 "RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff",
                 "RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff"
                 ]

switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)

for mod in eleVIDModules:
    setupAllVIDIdsInModule(process, mod, setupVIDElectronSelection)

# update some VID modules to work with potentially changed electron collections
process.egmGsfElectronIDs.physicsObjectSrc = electronCollection
process.electronRegressionValueMapProducer.srcMiniAOD = electronCollection
process.electronMVAValueMapProducer.srcMiniAOD = electronCollection


#
# MET corrections and uncertainties
#

if options.recorrectMET:
    # use the standard tool
    from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
    # do not use a postfix here!
    runMetCorAndUncFromMiniAOD(process,
       isData           = options.realData,
       electronColl     = electronCollection.value(),
       muonColl         = muonCollection.value(),
       tauColl          = tauCollection.value(),
       photonColl       = photonCollection.value(),
       jetCollUnskimmed = jetCollection.value(),
       recoMetFromPFCs  = True
   )

    # overwrite output collections
    METCollection = cms.InputTag("slimmedMETs", "", process.name_())

    # also add MET corrections due to e/g corrections, such as the slew rate fix in reMiniAOD
    if options.realData:
        from PhysicsTools.PatUtils.tools.corMETFromMuonAndEG import corMETFromMuonAndEG
        corMETFromMuonAndEG(process,
           pfCandCollection      = "",
           electronCollection    = "slimmedElectronsBeforeGSFix",
           photonCollection      = "slimmedPhotonsBeforeGSFix",
           corElectronCollection = electronCollection.value(),
           corPhotonCollection   = photonCollection.value(),
           allMETEGCorrected     = True,
           muCorrection          = False,
           eGCorrection          = True,
           runOnMiniAOD          = True,
           postfix               = "MuEGClean"
        )
        process.slimmedMETsMuEGClean = process.slimmedMETsRecorrected.clone(
           src             = cms.InputTag("patPFMetT1MuEGClean"),
           rawVariation    = cms.InputTag("patPFMetRawMuEGClean"),
           t1Uncertainties = cms.InputTag("patPFMetT1%sMuEGClean")
        )
        del process.slimmedMETsMuEGClean.caloMET

        # overwrite output collections
        METCollection = cms.InputTag("slimmedMETsMuEGClean", "", process.name_())


    # patch the phi correction parameter sets that are used in runMetCorAndUncFromMiniAOD,
    # we only need to overwrite patMultPhiCorrParams_T1Txy_25ns with the new one
    if options.realData:
        if options.dataEra in ("2016B", "2016C", "2016D", "2016E", "2016F"):
            from MetTools.MetPhiCorrections.tools.multPhiCorr_ReMiniAOD_Data_BCDEF_80X_sumPt_cfi \
                import multPhiCorr_Data_BCDEF_80X as metPhiCorrParams
        else: # "2016G", "2016Hv2", "2016Hv3"
            from MetTools.MetPhiCorrections.tools.multPhiCorr_ReMiniAOD_Data_GH_80X_sumPt_cfi \
                import multPhiCorr_Data_GH_80X as metPhiCorrParams
    else:
        from MetTools.MetPhiCorrections.tools.multPhiCorr_Summer16_MC_DY_80X_sumPt_cfi \
            import multPhiCorr_MC_DY_sumPT_80X as metPhiCorrParams
    # actual patch
    getattr(process, "patPFMetTxyCorr").parameters = cms.VPSet(pset for pset in metPhiCorrParams)

#
# custom MET filters
#

process.load("RecoMET.METFilters.BadPFMuonFilter_cfi")
process.BadPFMuonFilter.muons        = muonCollection
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

process.load("RecoMET.METFilters.BadChargedCandidateFilter_cfi")
process.BadChargedCandidateFilter.muons        = muonCollection
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")

process.load("RecoMET.METFilters.badGlobalMuonTaggersMiniAOD_cff")
process.badGlobalMuonTaggerMAOD.muons         = muonCollection
process.badGlobalMuonTaggerMAOD.taggingMode   = cms.bool(True)
process.cloneGlobalMuonTaggerMAOD.muons       = muonCollection
process.cloneGlobalMuonTaggerMAOD.taggingMode = cms.bool(True)

seq += process.BadPFMuonFilter + process.BadChargedCandidateFilter + process.badGlobalMuonTaggerMAOD + process.cloneGlobalMuonTaggerMAOD

#
# update PUJetId values
#
if options.updatePUJetId:
    process.load("RecoJets.JetProducers.PileupJetID_cfi")
    process.pileupJetIdUpdated = process.pileupJetId.clone(
      jets             = jetCollection,
      vertexes         = cms.InputTag("offlineSlimmedPrimaryVertices"),
      inputIsCorrected = cms.bool(True),
      applyJec         = cms.bool(True)
    )
    seq += process.pileupJetIdUpdated

    process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")
    process.updatedPatJets.jetSource         = jetCollection
    process.updatedPatJets.addJetCorrFactors = cms.bool(False)
    process.updatedPatJets.userData.userFloats.src.append("pileupJetIdUpdated:fullDiscriminant")
    process.updatedPatJets.userData.userInts.src.append("pileupJetIdUpdated:fullId")

    # overwrite output collections
    jetCollection = cms.InputTag("updatedPatJets", "", process.name_())


#
# ttbar related setup
#

#if options.isTtbar:
#    process.load("FOO.BAR.ttbarSequence_cff") # replace FOO.BAR with your subystem.module



# new electron MVA developed by the EGamma POG
#process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")
#from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
#switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
#my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff']
#for idmod in my_id_modules:
#	setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)



# load the analysis:
process.load("Analyzers.ttH_analyzer.ttHbb_cfi")

# pat object collections
process.ttHbb.input_tags.electrons = electronCollection
process.ttHbb.input_tags.muons     = muonCollection
process.ttHbb.input_tags.mets      = METCollection
process.ttHbb.input_tags.jets     = jetCollection

#electron VID collections
#process.ttHbb.electronVIDCollections = cms.VInputTag(
#     "egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp80",
#     "egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp90",
#     "egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp80",
#     "egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp90",
#     "egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto",
#     "egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose",
#     "egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium",
#     "egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"
#)

# additional MET filters
process.ttHbb.additionalMETFilterCollections = cms.VInputTag(
     "BadPFMuonFilter",
     "BadChargedCandidateFilter",
     "badGlobalMuonTaggerMAOD",
     "cloneGlobalMuonTaggerMAOD"
)

# ttbar-related collections
process.ttHbb.ttHFGenFilterCollection = cms.InputTag("ttHFGenFilter")



#ttHf categorization

# Setting input particle collections to be used by the tools
genJetCollection = 'ak4GenJetsCustom'
#genJetCollection = 'slimmedGenJets'
genParticleCollection = 'prunedGenParticles'
genJetInputParticleCollection = 'packedGenParticles'
    
## producing a subset of particles to be used for jet clustering
from RecoJets.Configuration.GenJetParticles_cff import genParticlesForJetsNoNu
process.genParticlesForJetsNoNu = genParticlesForJetsNoNu.clone(
	src = genJetInputParticleCollection
)
seq += process.genParticlesForJetsNoNu

    
# Producing own jets for testing purposes
from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
process.ak4GenJetsCustom = ak4GenJets.clone(
        src = 'genParticlesForJetsNoNu',
        #    src = genJetInputParticleCollection,
        rParam = cms.double(0.4),
        jetAlgorithm = cms.string("AntiKt")
)
seq += process.ak4GenJetsCustom
    
# Ghost particle collection used for Hadron-Jet association 
# MUST use proper input particle collection
from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartons
process.selectedHadronsAndPartons = selectedHadronsAndPartons.clone(
        particles = genParticleCollection
)
seq += process.selectedHadronsAndPartons
    
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
seq += process.genJetFlavourInfos
    
from PhysicsTools.JetMCAlgos.GenHFHadronMatcher_cff import *
# Plugin for analysing B hadrons
# MUST use the same particle collection as in selectedHadronsAndPartons
process.matchGenBHadron = matchGenBHadron.clone(
        genParticles = genParticleCollection,
        jetFlavourInfos = "genJetFlavourInfos"
)
seq += process.matchGenBHadron
    
# Plugin for analysing C hadrons
# MUST use the same particle collection as in selectedHadronsAndPartons
process.matchGenCHadron = matchGenCHadron.clone(
        genParticles = genParticleCollection,
        jetFlavourInfos = "genJetFlavourInfos"
)
seq += process.matchGenCHadron
    
## Producer for ttbar categorisation ID
# MUST use same genJetCollection as used for tools above
from PhysicsTools.JetMCAlgos.GenTtbarCategorizer_cfi import categorizeGenTtbar
process.categorizeGenTtbar = categorizeGenTtbar.clone(
        genJetPtMin = 20.,
        genJetAbsEtaMax = 2.4,
        genJets = genJetCollection,
)

seq += process.categorizeGenTtbar

process.TFileService = cms.Service("TFileService",
	fileName = cms.string("ttHbbNtuple.root")
)

process.p = cms.Path(
    #process.electronMVAValueMapProducer
    #process.egmGsfElectronIDSequence
    #* process.genParticlesForJetsNoNu
    #* process.ak4GenJetsCustom
    #* process.selectedHadronsAndPartons
    #* process.genJetFlavourInfos
    #* process.matchGenBHadron
    #* process.matchGenCHadron
    #* process.categorizeGenTtbar
    #process.fullPatMetSequence
    #process.BadPFMuonFilter
    #process.BadChargedCandidateFilter
    seq + process.ttHbb
)
