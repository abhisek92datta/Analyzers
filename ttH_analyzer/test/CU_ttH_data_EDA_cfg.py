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
#process.GlobalTag.globaltag = '80X_dataRun2_Prompt_ICHEP16JEC_v0'

#For Datasets B-G
process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v6'
#For Dataset H
#process.GlobalTag.globaltag = '80X_dataRun2_Prompt_v15'

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))
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
process.ak4PFchsResidual = ak4CaloResidual.clone( algorithm = 'AK4PFchs' )

process.ak4PFchsL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
    'ak4PFCHSL1Fastjet', 
    'ak4PFchsL2Relative', 
    'ak4PFchsL3Absolute',
    'ak4PFchsResidual')
)


#from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
#updateJetCollection(
#  process,
#  jetSource = cms.InputTag('slimmedJets'),
#  jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None' ),  
#  btagDiscriminators = ['pfCombinedInclusiveSecondaryVertexV2BJetTags'],
#  runIVF=True,
#  #btagPrefix = 'new' # optional, in case interested in accessing both the old and new discriminator values
#)

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
# If you only want to re-correct and get the proper uncertainties
runMetCorAndUncFromMiniAOD(process,
 isData=True,
 #jecUncFile = "data/JEC/Spring16_25nsV6_DATA_Uncertainty_AK4PFchs.txt"
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
        ## Single Electron
        '/store/data/Run2016B/SingleElectron/MINIAOD/23Sep2016-v3/00000/00099863-E799-E611-A876-141877343E6D.root',
	'/store/data/Run2016B/SingleElectron/MINIAOD/23Sep2016-v3/00000/0077B9D8-1E9A-E611-963A-0CC47A57CCE8.root',
        '/store/data/Run2016C/SingleElectron/MINIAOD/23Sep2016-v1/50000/001B31C0-248C-E611-B7BD-0025905B860E.root',
        '/store/data/Run2016C/SingleElectron/MINIAOD/23Sep2016-v1/50000/001C375D-098B-E611-9453-0025905B85E8.root',
        '/store/data/Run2016D/SingleElectron/MINIAOD/23Sep2016-v1/70000/04E8F72C-AF89-E611-9D2F-FA163E1D7951.root',
        '/store/data/Run2016D/SingleElectron/MINIAOD/23Sep2016-v1/70000/0654EC37-FA8A-E611-820F-FA163E066046.root',
        '/store/data/Run2016E/SingleElectron/MINIAOD/23Sep2016-v1/100000/00827A71-F98C-E611-9639-0CC47A4D7650.root',
	'/store/data/Run2016E/SingleElectron/MINIAOD/23Sep2016-v1/100000/008A51FD-AC8C-E611-87AE-008CFA197E84.root',
	'/store/data/Run2016F/SingleElectron/MINIAOD/23Sep2016-v1/100000/00FFEF1F-EA96-E611-9A24-0090FAA583C4.root',
	'/store/data/Run2016F/SingleElectron/MINIAOD/23Sep2016-v1/100000/04173AAE-F893-E611-8B11-002590E7D5AE.root',
	'/store/data/Run2016G/SingleElectron/MINIAOD/23Sep2016-v1/100000/004A7893-A990-E611-B29F-002590E7DE36.root',
	'/store/data/Run2016G/SingleElectron/MINIAOD/23Sep2016-v1/100000/0083C5D5-968F-E611-B630-7845C4FC35CF.root',
	'/store/data/Run2016H/SingleElectron/MINIAOD/PromptReco-v2/000/281/207/00000/989A10BD-6282-E611-975A-02163E011C96.root',
	'/store/data/Run2016H/SingleElectron/MINIAOD/PromptReco-v2/000/281/209/00000/884F1086-6382-E611-98D8-02163E0143CB.root',
	'/store/data/Run2016H/SingleElectron/MINIAOD/PromptReco-v3/000/284/036/00000/1257A268-599F-E611-A437-02163E011E7A.root',
	'/store/data/Run2016H/SingleElectron/MINIAOD/PromptReco-v3/000/284/036/00000/1CBE1DEB-589F-E611-ABBB-02163E0143B5.root'
        
        ## Single Muon
        #'/store/data/Run2016B/SingleMuon/MINIAOD/23Sep2016-v3/00000/00AE0629-1F98-E611-921A-008CFA1112CC.root',
        #'/store/data/Run2016B/SingleMuon/MINIAOD/23Sep2016-v3/00000/041C2515-1698-E611-B909-901B0E542804.root',
        #'/store/data/Run2016C/SingleMuon/MINIAOD/23Sep2016-v1/70000/001F13A2-7E8D-E611-B910-FA163E782438.root',
        #'/store/data/Run2016C/SingleMuon/MINIAOD/23Sep2016-v1/70000/00DFF8E4-828C-E611-AF50-0CC47A0AD6AA.root',
        #'/store/data/Run2016D/SingleMuon/MINIAOD/23Sep2016-v1/010000/0201B90C-C79B-E611-8763-00266CFFBF80.root',
        #'/store/data/Run2016D/SingleMuon/MINIAOD/23Sep2016-v1/010000/0CEAC70F-C79B-E611-BE1E-00266CFFC7CC.root',
        #'/store/data/Run2016E/SingleMuon/MINIAOD/23Sep2016-v1/50000/00CFC689-8D8D-E611-9F90-0CC47A13D16E.root',
        #'/store/data/Run2016E/SingleMuon/MINIAOD/23Sep2016-v1/50000/0230DB91-868D-E611-A532-0025904A96BC.root',
        #'/store/data/Run2016F/SingleMuon/MINIAOD/23Sep2016-v1/100000/0A06DDB0-0A8D-E611-9714-0242AC130002.root',
        #'/store/data/Run2016F/SingleMuon/MINIAOD/23Sep2016-v1/100000/0CFCD6D1-FC8C-E611-A7C5-A0369F6369D2.root'
        #'/store/data/Run2016G/SingleMuon/MINIAOD/23Sep2016-v1/1110000/72446D9C-D89C-E611-9060-002590A3C984.root',
        #'/store/data/Run2016G/SingleMuon/MINIAOD/23Sep2016-v1/1110000/9039FEA1-B19C-E611-9CD9-7845C4FC361A.root',
        #'/store/data/Run2016H/SingleMuon/MINIAOD/PromptReco-v2/000/281/207/00000/C01B8838-6282-E611-9884-02163E01414B.root',
        #'/store/data/Run2016H/SingleMuon/MINIAOD/PromptReco-v2/000/281/209/00000/CC401A44-6682-E611-A6DE-02163E0141D8.root',
        #'/store/data/Run2016H/SingleMuon/MINIAOD/PromptReco-v3/000/284/036/00000/0E02D50E-989F-E611-A962-FA163EE15C80.root',
        #'/store/data/Run2016H/SingleMuon/MINIAOD/PromptReco-v3/000/284/036/00000/129CD4B5-5D9F-E611-A9AB-02163E014220.root'
        
        ## Double Electron
        #'/store/data/Run2016B/DoubleEG/MINIAOD/23Sep2016-v3/00000/0060C751-C097-E611-9FE6-FA163EFD4308.root',
	#'/store/data/Run2016B/DoubleEG/MINIAOD/23Sep2016-v3/00000/026637B9-C797-E611-B157-14187741212B.root',
	#'/store/data/Run2016C/DoubleEG/MINIAOD/23Sep2016-v1/100000/00831315-BA89-E611-80F4-0CC47A7C3412.root',
	#'/store/data/Run2016C/DoubleEG/MINIAOD/23Sep2016-v1/100000/1669252F-BA89-E611-87AF-0025905A48F0.root',
	#'/store/data/Run2016D/DoubleEG/MINIAOD/23Sep2016-v1/100000/206CD6B5-AE87-E611-8B2B-0CC47A4D769A.root',
	#'/store/data/Run2016D/DoubleEG/MINIAOD/23Sep2016-v1/100000/2620E901-AF87-E611-B7A5-0CC47A4D76B8.root',
	#'/store/data/Run2016E/DoubleEG/MINIAOD/23Sep2016-v1/80000/00A347CF-E289-E611-974D-1CC1DE192766.root',
	#'/store/data/Run2016E/DoubleEG/MINIAOD/23Sep2016-v1/80000/022F96B2-1989-E611-BAFB-0025905746AA.root',
	#'/store/data/Run2016F/DoubleEG/MINIAOD/23Sep2016-v1/50000/0068A385-D989-E611-A227-008CFA166188.root',
	#'/store/data/Run2016F/DoubleEG/MINIAOD/23Sep2016-v1/50000/00B858F2-8188-E611-AC6A-0CC47A4C8E46.root',
	#'/store/data/Run2016G/DoubleEG/MINIAOD/23Sep2016-v1/100000/0608426D-3F8E-E611-A52D-00237DF28460.root',
	#'/store/data/Run2016G/DoubleEG/MINIAOD/23Sep2016-v1/100000/06FC0D6B-F88A-E611-A24D-00215E2EB6EE.root',
	#'/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v2/000/281/207/00000/FAC8E8C0-6482-E611-B1B3-02163E013671.root',
	#'/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v2/000/281/209/00000/3285CB51-6682-E611-9D1E-02163E014751.root',
	#'/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v3/000/284/036/00000/1878DF24-619F-E611-A962-02163E0146C8.root',
	#'/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v3/000/284/036/00000/2AC0860C-5F9F-E611-9235-FA163EA23E61.root'
        
        ## Electron Muon
        #'/store/data/Run2016B/MuonEG/MINIAOD/23Sep2016-v3/00000/024ADA16-1F98-E611-AD32-0242AC130005.root',
	#'/store/data/Run2016B/MuonEG/MINIAOD/23Sep2016-v3/00000/02DEC0BD-2098-E611-99BD-0242AC130003.root',
	#'/store/data/Run2016C/MuonEG/MINIAOD/23Sep2016-v1/70000/003F4D0A-5E87-E611-9105-0242AC130002.root',
	#'/store/data/Run2016C/MuonEG/MINIAOD/23Sep2016-v1/70000/0059AB0D-2289-E611-B1A2-00269E95B1BC.root',
	#'/store/data/Run2016D/MuonEG/MINIAOD/23Sep2016-v1/100000/02140D15-8D89-E611-8930-FA163E65CBE5.root',
	#'/store/data/Run2016D/MuonEG/MINIAOD/23Sep2016-v1/100000/04C348B2-2689-E611-BF3B-02163E013EFC.root',
	#'/store/data/Run2016E/MuonEG/MINIAOD/23Sep2016-v1/100000/023BB6F0-588F-E611-9F5A-0CC47A4C8E7E.root',
	#'/store/data/Run2016E/MuonEG/MINIAOD/23Sep2016-v1/100000/08F0E9FD-258F-E611-B593-0242AC130004.root',
	#'/store/data/Run2016F/MuonEG/MINIAOD/23Sep2016-v1/100000/026529CA-F891-E611-A364-02163E00E615.root',
	#'/store/data/Run2016F/MuonEG/MINIAOD/23Sep2016-v1/100000/0C56508A-D492-E611-BC09-0090FAA58204.root',
	#'/store/data/Run2016G/MuonEG/MINIAOD/23Sep2016-v1/100000/005AB7E9-0B93-E611-AC81-848F69FD2925.root',
	#'/store/data/Run2016G/MuonEG/MINIAOD/23Sep2016-v1/100000/005FB5C1-ED8F-E611-BAE7-0025905A607E.root',
	#'/store/data/Run2016H/MuonEG/MINIAOD/PromptReco-v2/000/281/207/00000/3AB357FC-6282-E611-B7EE-02163E014751.root',
	#'/store/data/Run2016H/MuonEG/MINIAOD/PromptReco-v2/000/281/209/00000/5645FE65-6382-E611-995D-FA163E823565.root',
	#'/store/data/Run2016H/MuonEG/MINIAOD/PromptReco-v3/000/284/036/00000/2C4B1CAE-A99F-E611-8739-02163E0141DA.root',
	#'/store/data/Run2016H/MuonEG/MINIAOD/PromptReco-v3/000/284/036/00000/4A12488F-C19F-E611-BFF9-02163E01456E.root'
        
        ## Double Muon
        #'/store/data/Run2016B/DoubleMuon/MINIAOD/23Sep2016-v3/00000/00269AA1-479B-E611-A359-0025905A6082.root',
        #'/store/data/Run2016B/DoubleMuon/MINIAOD/23Sep2016-v3/00000/002F0F59-3F9B-E611-95B9-0CC47A7C354A.root',
        #'/store/data/Run2016C/DoubleMuon/MINIAOD/23Sep2016-v1/80000/0044DA75-708B-E611-8F8B-008CFA1974A4.root',
        #'/store/data/Run2016C/DoubleMuon/MINIAOD/23Sep2016-v1/80000/005599F4-5787-E611-A034-0025905C54C6.root',
        #'/store/data/Run2016D/DoubleMuon/MINIAOD/23Sep2016-v1/100000/022B1BC6-A789-E611-808B-B499BAAB427C.root',
        #'/store/data/Run2016D/DoubleMuon/MINIAOD/23Sep2016-v1/100000/0274567D-078A-E611-91AF-00266CFCC988.root',
        #'/store/data/Run2016E/DoubleMuon/MINIAOD/23Sep2016-v1/100000/0049B6D2-278C-E611-AEB4-0025902D944E.root',
        #'/store/data/Run2016E/DoubleMuon/MINIAOD/23Sep2016-v1/100000/0694CCB8-578C-E611-AF97-00304867FDAB.root',
        #'/store/data/Run2016F/DoubleMuon/MINIAOD/23Sep2016-v1/50000/003E1D68-618E-E611-A816-0CC47AD9908C.root',
        #'/store/data/Run2016F/DoubleMuon/MINIAOD/23Sep2016-v1/50000/040EDEBA-0490-E611-A424-008CFA110C68.root',
        #'/tore/data/Run2016G/DoubleMuon/MINIAOD/23Sep2016-v1/100000/00993A51-DF90-E611-A4EE-7845C4FC3650.root',
        #'/tore/data/Run2016G/DoubleMuon/MINIAOD/23Sep2016-v1/100000/00DD00F8-008C-E611-8CD0-00266CFFC9C4.root',
        #'/store/data/Run2016H/DoubleMuon/MINIAOD/PromptReco-v2/000/281/207/00000/2CC1A774-6382-E611-AD53-FA163E7F172A.root',
        #'/store/data/Run2016H/DoubleMuon/MINIAOD/PromptReco-v2/000/281/209/00000/6AC26751-6482-E611-B678-FA163E0B3058.root',
        #'/store/data/Run2016H/DoubleMuon/MINIAOD/PromptReco-v3/000/284/036/00000/04DC0281-C89F-E611-81C6-02163E0141E6.root',
        #'/store/data/Run2016H/DoubleMuon/MINIAOD/PromptReco-v3/000/284/036/00000/5AE9F885-A19F-E611-989B-02163E014742.root' 
        )
)

import FWCore.PythonUtilities.LumiList as LumiList
#process.source.lumisToProcess = LumiList.LumiList(filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.txt').getVLuminosityBlockRange()
#process.source.lumisToProcess = LumiList.LumiList(filename = 'data/JSON/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON_unprescaled.txt').getVLuminosityBlockRange()
process.source.lumisToProcess = LumiList.LumiList(filename = 'data/JSON/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON_unprescaled.txt').getVLuminosityBlockRange()


# new electron MVA developed by the EGamma POG 
#process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff']
for idmod in my_id_modules:  
	setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

# load the analysis:
process.load("Analyzers.ttH_analyzer.ttHbb_cfi")

process.TFileService = cms.Service("TFileService",
	fileName = cms.string('ttHbbNtuple.root')
)

process.p = cms.Path(
    #process.electronMVAValueMapProducer
    process.egmGsfElectronIDSequence
    * process.fullPatMetSequence
    * process.BadPFMuonFilter
    * process.BadChargedCandidateFilter
    * process.ttHbb
)
