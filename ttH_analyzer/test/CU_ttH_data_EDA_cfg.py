import os
import FWCore.ParameterSet.Config as cms
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
                 True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "input dataset contains real data"
                 )
options.register("dataEra",
                 "2016B",
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
   ## Single Electron
   '/store/data/Run2016B/SingleElectron/MINIAOD/03Feb2017_ver2-v2/110000/003B2C1F-50EB-E611-A8F1-002590E2D9FE.root',
   '/store/data/Run2016B/SingleElectron/MINIAOD/03Feb2017_ver2-v2/110000/0043D133-4DEB-E611-BFE4-003048F5ADEC.root',
   #'/store/data/Run2016C/SingleElectron/MINIAOD/03Feb2017-v1/100000/02169BE7-81EB-E611-BB99-02163E0137CD.root',
   #'/store/data/Run2016C/SingleElectron/MINIAOD/03Feb2017-v1/100000/0244B3B4-78EB-E611-9D42-0CC47A009E24.root',
   #'/store/data/Run2016D/SingleElectron/MINIAOD/03Feb2017-v1/110000/001A5DEC-77EB-E611-95CD-0CC47A4C8EE8.root',
   #'/store/data/Run2016D/SingleElectron/MINIAOD/03Feb2017-v1/110000/001F5992-6DEA-E611-AF62-A4BF01013DD5.root',
   #'/store/data/Run2016E/SingleElectron/MINIAOD/03Feb2017-v1/110000/00022327-8BEA-E611-86CB-0025905B8566.root',
   #'/store/data/Run2016E/SingleElectron/MINIAOD/03Feb2017-v1/110000/002070E4-1BEB-E611-BA80-0025905A608E.root',
   #'/store/data/Run2016F/SingleElectron/MINIAOD/03Feb2017-v1/100000/00B336D6-6AEC-E611-8581-E0071B7AC7B0.root',
   #'/store/data/Run2016F/SingleElectron/MINIAOD/03Feb2017-v1/100000/0229A4E9-5DEC-E611-B931-A0000420FE80.root',
   #'/store/data/Run2016G/SingleElectron/MINIAOD/03Feb2017-v1/50000/004A75AB-B2EA-E611-B000-24BE05CEFDF1.root',
   #'/store/data/Run2016G/SingleElectron/MINIAOD/03Feb2017-v1/50000/00539EDD-F8EA-E611-AC1D-0025900B5648.root',',
   #'/store/data/Run2016H/SingleElectron/MINIAOD/03Feb2017_ver2-v1/100000/00553E5F-29EC-E611-ADB0-00259074AE8C.root',
   #'/store/data/Run2016H/SingleElectron/MINIAOD/03Feb2017_ver2-v1/100000/006EC85C-29EC-E611-B2E1-20CF3019DF17.root',
   #'/store/data/Run2016H/SingleElectron/MINIAOD/03Feb2017_ver3-v1/110000/02973E99-69EC-E611-9913-5065F381A2F1.root',
   #'/store/data/Run2016H/SingleElectron/MINIAOD/03Feb2017_ver3-v1/110000/04C8C9AF-62EC-E611-AB90-A0000420FE80.root'

   ## Single Muon
   #'/store/data/Run2016B/SingleMuon/MINIAOD/03Feb2017_ver2-v2/100000/000C6E52-8BEC-E611-B3FF-0025905C42FE.root',
   #'/store/data/Run2016B/SingleMuon/MINIAOD/03Feb2017_ver2-v2/100000/001E3E7D-57EB-E611-8469-0CC47A7C35D2.root',
   #'/store/data/Run2016C/SingleMuon/MINIAOD/03Feb2017-v1/50000/001CF316-1AEB-E611-BBBD-0CC47A4C8EE2.root',
   #'/store/data/Run2016C/SingleMuon/MINIAOD/03Feb2017-v1/50000/0022D65B-05EB-E611-84E9-0025905A6104.root',
   #'/store/data/Run2016D/SingleMuon/MINIAOD/03Feb2017-v1/100000/00622F98-20EB-E611-A0A4-28924A33AFF6.root',
   #'/store/data/Run2016D/SingleMuon/MINIAOD/03Feb2017-v1/100000/023071E7-97EA-E611-A89A-0025904C67B6.root',
   #'/store/data/Run2016E/SingleMuon/MINIAOD/03Feb2017-v1/110000/001D567A-0CEB-E611-A438-D8D385AE8848.root',
   #'/store/data/Run2016E/SingleMuon/MINIAOD/03Feb2017-v1/110000/004309C7-E6EA-E611-92B5-0025905A60DA.root',
   #'/store/data/Run2016F/SingleMuon/MINIAOD/03Feb2017-v1/100000/00E6AB6D-BCEC-E611-8F6E-0025905C3D98.root',
   #'/store/data/Run2016F/SingleMuon/MINIAOD/03Feb2017-v1/100000/040B13AE-C3EC-E611-8082-0025904C6564.root',
   #'/store/data/Run2016G/SingleMuon/MINIAOD/03Feb2017-v1/100000/00E6DF50-70EA-E611-ACC4-0CC47A1E089C.root',
   #'/store/data/Run2016G/SingleMuon/MINIAOD/03Feb2017-v1/100000/02382B19-D1EA-E611-B2F9-0CC47ABAC11C.root',
   #'/store/data/Run2016H/SingleMuon/MINIAOD/03Feb2017_ver2-v1/110000/00633FF0-85EA-E611-811C-001E674FB25C.root',
   #'/store/data/Run2016H/SingleMuon/MINIAOD/03Feb2017_ver2-v1/110000/006517CB-8AEA-E611-8CF6-0CC47AC08BD4.root',
   #'/store/data/Run2016H/SingleMuon/MINIAOD/03Feb2017_ver3-v1/80000/0040ECBB-76EA-E611-8FE7-A0000420FE80.root',
   #'/store/data/Run2016H/SingleMuon/MINIAOD/03Feb2017_ver3-v1/80000/16F28614-84EA-E611-8083-A0369F310374.root'

   ## Double Electron
   #'/store/data/Run2016B/DoubleEG/MINIAOD/03Feb2017_ver2-v2/50000/00054938-CEEA-E611-889E-0CC47A4D7650.root',
   #'/store/data/Run2016B/DoubleEG/MINIAOD/03Feb2017_ver2-v2/50000/001BC16F-90EA-E611-87E2-F04DA27540BB.root',
   #'/store/data/Run2016C/DoubleEG/MINIAOD/03Feb2017-v1/80000/00371362-6AEC-E611-9845-842B2B758BAA.root',
   #'/store/data/Run2016C/DoubleEG/MINIAOD/03Feb2017-v1/80000/006E86B9-77EC-E611-BA8F-02163E019CE7.root',
   #'/store/data/Run2016D/DoubleEG/MINIAOD/03Feb2017-v1/100000/002CE21C-0BEB-E611-8597-001E67E6F8E6.root',
   #'/store/data/Run2016D/DoubleEG/MINIAOD/03Feb2017-v1/100000/006BB826-05EB-E611-97F4-0025904A8ED2.root',
   #'/store/data/Run2016E/DoubleEG/MINIAOD/03Feb2017-v1/110000/003AF399-ABEA-E611-92CF-002590E2DA08.root',
   #'/store/data/Run2016E/DoubleEG/MINIAOD/03Feb2017-v1/110000/048E7FCD-98EA-E611-AB8A-002590D9D966.root',
   #'/store/data/Run2016F/DoubleEG/MINIAOD/03Feb2017-v1/80000/0006AFD8-F8EA-E611-9F9D-0CC47A13D09C.root',
   #'/store/data/Run2016F/DoubleEG/MINIAOD/03Feb2017-v1/80000/00A7E4D9-F8EA-E611-A62B-002590E3A004.root',
   #'/store/data/Run2016G/DoubleEG/MINIAOD/03Feb2017-v1/100000/002F14FF-D0EA-E611-952E-008CFA197AF4.root',
   #'/store/data/Run2016G/DoubleEG/MINIAOD/03Feb2017-v1/100000/02642443-F0EA-E611-9D24-008CFA197D74.root',
   #'/store/data/Run2016H/DoubleEG/MINIAOD/03Feb2017_ver2-v1/100000/023E858B-F7EC-E611-889C-047D7BD6DDF2.root',
   #'/store/data/Run2016H/DoubleEG/MINIAOD/03Feb2017_ver2-v1/100000/02AD337A-95ED-E611-B1AE-047D7B881D90.root',
   #'/store/data/Run2016H/DoubleEG/MINIAOD/03Feb2017_ver3-v1/1030000/D41C6358-4DF0-E611-BBAC-002590DB927A.root',
   #'/store/data/Run2016H/DoubleEG/MINIAOD/03Feb2017_ver3-v1/110000/78C7FC4D-7CED-E611-870F-0CC47A7FC858.root'


   ## Electron Muon
   #'/store/data/Run2016B/MuonEG/MINIAOD/03Feb2017_ver2-v2/100000/0017ADAB-F6EC-E611-8C01-0090FAA57EA4.root',
   #'/store/data/Run2016B/MuonEG/MINIAOD/03Feb2017_ver2-v2/100000/00680B67-A0EC-E611-ABC5-00259073E4E2.root',
   #'/store/data/Run2016C/MuonEG/MINIAOD/03Feb2017-v1/110000/0649FFBD-50EC-E611-9875-24BE05CEEDE1.root',
   #'/store/data/Run2016C/MuonEG/MINIAOD/03Feb2017-v1/110000/0AEA5448-35EC-E611-B1AF-A0369F310374.root',
   #'/store/data/Run2016D/MuonEG/MINIAOD/03Feb2017-v1/80000/02264DFC-6EEB-E611-95AF-0090FAA572B0.root',
   #'/store/data/Run2016D/MuonEG/MINIAOD/03Feb2017-v1/80000/0261593F-D6EA-E611-8577-0090FAA57630.root',
   #'/store/data/Run2016E/MuonEG/MINIAOD/03Feb2017-v1/110000/00003F8D-05EB-E611-92D9-02163E019DC6.root',
   #'/store/data/Run2016E/MuonEG/MINIAOD/03Feb2017-v1/110000/005F2FB7-7CEB-E611-9454-0CC47A745298.root',
   #'/store/data/Run2016F/MuonEG/MINIAOD/03Feb2017-v1/50000/0496325A-05EB-E611-953B-0025905A60DE.root',
   #'/store/data/Run2016F/MuonEG/MINIAOD/03Feb2017-v1/50000/0C20A32A-A6EA-E611-ADFB-0025905A60F8.root',
   #'/store/data/Run2016G/MuonEG/MINIAOD/03Feb2017-v1/100000/08127890-99EC-E611-82BE-0CC47A7E6A6C.root',
   #'/store/data/Run2016G/MuonEG/MINIAOD/03Feb2017-v1/100000/0ABB558F-99EC-E611-8484-002590DE6E6E.root',
   #'/store/data/Run2016H/MuonEG/MINIAOD/03Feb2017_ver2-v1/100000/044366C7-4AEE-E611-8CF7-0025905B856E.root',
   #'/store/data/Run2016H/MuonEG/MINIAOD/03Feb2017_ver2-v1/100000/04BB134D-3BEE-E611-B80E-0025905A48D6.root',
   #'/store/data/Run2016H/MuonEG/MINIAOD/03Feb2017_ver3-v1/110000/10D7E08C-96EA-E611-9119-0025905B858E.root',
   #'/store/data/Run2016H/MuonEG/MINIAOD/03Feb2017_ver3-v1/110000/1468B227-F4EA-E611-BF5A-0025905A60B6.root'

  ## Double Muon
  #'/store/data/Run2016B/DoubleMuon/MINIAOD/03Feb2017_ver2-v2/100000/0005AD9F-64ED-E611-A952-0CC47A78A42C.root',
   #'/store/data/Run2016B/DoubleMuon/MINIAOD/03Feb2017_ver2-v2/100000/00415FAC-B5EC-E611-A1C9-00266CF3E130.root',
   #'/store/data/Run2016C/DoubleMuon/MINIAOD/03Feb2017-v1/110000/02292378-59ED-E611-BA0B-0CC47A4D768E.root',
   #'/store/data/Run2016C/DoubleMuon/MINIAOD/03Feb2017-v1/110000/022E014E-5FED-E611-BDCC-0025905B858E.root',
   #'/store/data/Run2016D/DoubleMuon/MINIAOD/03Feb2017-v1/100000/00E0F123-F7ED-E611-9F34-0CC47A7FC736.root',
   #'/store/data/Run2016D/DoubleMuon/MINIAOD/03Feb2017-v1/100000/025BE3ED-EBED-E611-95E9-047D7BD6DDA4.root',
   #'/store/data/Run2016E/DoubleMuon/MINIAOD/03Feb2017-v1/100000/022FEC03-9AED-E611-9AE9-0025905A60B2.root',
   #'/store/data/Run2016E/DoubleMuon/MINIAOD/03Feb2017-v1/100000/04328A2E-13ED-E611-BC83-0CC47A4D76D2.root',
   #'/store/data/Run2016F/DoubleMuon/MINIAOD/03Feb2017-v1/100000/0055F5B5-76EB-E611-B313-002590D9D9FC.root',
   #'/store/data/Run2016F/DoubleMuon/MINIAOD/03Feb2017-v1/100000/04B25C93-58EB-E611-86E4-002590FD5694.root',
   #'/store/data/Run2016G/DoubleMuon/MINIAOD/03Feb2017-v1/100000/00182C13-EEEA-E611-8897-001E675A6C2A.root',
   #'/store/data/Run2016G/DoubleMuon/MINIAOD/03Feb2017-v1/100000/007796A5-78EB-E611-8EFA-A4BF01011FD0.root',
   #'/store/data/Run2016H/DoubleMuon/MINIAOD/03Feb2017_ver2-v1/110000/08ADA6AA-D3EC-E611-AF17-B083FED42488.root',
   #'/store/data/Run2016H/DoubleMuon/MINIAOD/03Feb2017_ver2-v1/110000/0A1CCEC7-78EC-E611-863C-141877411FED.root',
   #'/store/data/Run2016H/DoubleMuon/MINIAOD/03Feb2017_ver3-v1/50000/36C4C2FB-4AEB-E611-ADD7-008CFA580778.root',
   #'/store/data/Run2016H/DoubleMuon/MINIAOD/03Feb2017_ver3-v1/50000/44C1C1FB-4AEB-E611-9597-008CFA1113F4.root'
   )
)

import FWCore.PythonUtilities.LumiList as LumiList
process.source.lumisToProcess = LumiList.LumiList(filename = 'data/JSON/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt').getVLuminosityBlockRange()

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
#process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")

process.load( "Configuration.StandardSequences.FrontierConditions_GlobalTag_cff" )
#process.GlobalTag.globaltag = '80X_dataRun2_Prompt_ICHEP16JEC_v0'


# Combined B-H
process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v7'

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))
process.options.allowUnscheduled = cms.untracked.bool(True)

process.maxEvents = cms.untracked.PSet(
	input = cms.untracked.int32(-1)
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
process.ak4PFchsResidual = ak4CaloResidual.clone( algorithm = 'AK4PFchs' )

process.ak4PFchsL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
    'ak4PFCHSL1Fastjet', 
    'ak4PFchsL2Relative', 
    'ak4PFchsL3Absolute',
    'ak4PFchsResidual')
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
      process.slimmedMETsMuEGClean = getattr(process, METCollection.getModuleLabel()).clone(
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

#seq += process.BadPFMuonFilter + process.BadChargedCandidateFilter + process.badGlobalMuonTaggerMAOD + process.cloneGlobalMuonTaggerMAOD


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
#    "egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp80",
#    "egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp90",
#    "egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp80",
#    "egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp90",
#    "egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto",
#    "egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose",
#    "egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium",
#    "egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"
#)

# additional MET filters
process.ttHbb.additionalMETFilterCollections = cms.VInputTag(
    "BadPFMuonFilter",
    "BadChargedCandidateFilter",
    "badGlobalMuonTaggerMAOD",
    "cloneGlobalMuonTaggerMAOD"
)

process.TFileService = cms.Service("TFileService",
	fileName = cms.string('ttHbbNtuple.root')
)

process.p = cms.Path(
    #process.electronMVAValueMapProducer
    #process.egmGsfElectronIDSequence
    #process.fullPatMetSequence
    #process.BadPFMuonFilter
    #process.BadChargedCandidateFilter
    seq + process.ttHbb
)
