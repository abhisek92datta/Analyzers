# Configurations for ttH bb analysis
#
# Abhisek Datta
#

import FWCore.ParameterSet.Config as cms


ttHbb =  cms.EDAnalyzer('CU_ttH_EDA',
        # Analysis type choice
        analysis_type = cms.string("dilepton"),
        # Generic
        verbosity = cms.bool(False),
        print_HLT_event_path = cms.bool(False),
        HLT_config_tag = cms.string('HLT2'),
        filter_config_tag = cms.string('PAT'),
        # Triggers
        collect_trigger_stats = cms.bool(False),
        ## Single lepton triggers:
        HLT_electron_triggers = cms.vstring(['HLT_Ele27_eta2p1_WPTight_Gsf_v']),
        HLT_muon_triggers = cms.vstring(['HLT_IsoMu22_v','HLT_IsoTkMu22_v']),
        ## Dilepton triggers:
        HLT_electron_electron_triggers = cms.vstring(['HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v']),
        HLT_electron_muon_triggers = cms.vstring( ['HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v','HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v']),
        HLT_muon_muon_triggers = cms.vstring(['HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v','HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v']),
        # Cuts
        min_ele_pT = cms.double(30.0),
        min_mu_pT = cms.double(25.0),
        min_veto_ele_pT = cms.double(15.0),
        min_veto_mu_pT = cms.double(15.0),
        min_di_ele1_pT = cms.double(25.0),
        min_di_ele2_pT = cms.double(15.0),
        min_di_mu1_pT = cms.double(25.0),
        min_di_mu2_pT = cms.double(15.0),
        min_jet_pT = cms.double(30.0),
        min_jet2_pT = cms.double(20.0),
        min_bjet_pT = cms.double(20.0),
        max_ele_eta = cms.double(2.1),
        max_mu_eta = cms.double(2.1),
        max_veto_ele_eta = cms.double(2.4),
        max_veto_mu_eta = cms.double(2.4),
        max_di_ele1_eta = cms.double(2.4),
        max_di_ele2_eta = cms.double(2.4),
        max_di_mu1_eta = cms.double(2.4),
        max_di_mu2_eta = cms.double(2.4),
        max_jet_eta = cms.double(2.4),
        max_bjet_eta = cms.double(2.4),
        min_njets = cms.int32(4),
        min_di_njets = cms.int32(2),
        min_nbtags = cms.int32(2),
        min_di_nbtags = cms.int32(1),
        min_di_mll = cms.double(20.0),
        min_di_met = cms.double(40.0),
        # MiniAODhelper
        using_real_data = cms.bool(False),
        ## available choices '-': none, 'L': loose, 'M': medium, 'T': tight
        b_tag_strength = cms.string('M'),
                          
        # InputTags
        input_tags = cms.PSet(
            pv = cms.InputTag("offlineSlimmedPrimaryVertices"),
            sv = cms.InputTag("slimmedSecondaryVertices"),
            pileup = cms.InputTag("addPileupInfo"),
            rho = cms.InputTag("fixedGridRhoFastjetAll"),
            electrons = cms.InputTag("slimmedElectrons"),
            muons = cms.InputTag("slimmedMuons"),
            jets = cms.InputTag("slimmedJets"),
            mets = cms.InputTag("slimmedMETs"),
            genjets = cms.InputTag("slimmedGenJets"),
            pfcand = cms.InputTag("packedPFCandidates"),
            beamspot = cms.InputTag("offlineBeamSpot"),
            mvaValues = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Values"),
            mvaCategories = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Categories"),
            genTtbarId = cms.InputTag("categorizeGenTtbar", "genTtbarId"),
            pileupinfo = cms.InputTag("slimmedAddPileupInfo"),
            lhepprod = cms.InputTag("externalLHEProducer")
        )
)
