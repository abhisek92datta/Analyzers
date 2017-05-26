# Configurations for ttH bb analysis
#
# Abhisek Datta
#

import FWCore.ParameterSet.Config as cms

from Analyzers.ttH_analyzer.CU_ttH_read_yaml import *

ttHbb =  cms.EDAnalyzer('CU_ttH_EDA',
        # Analysis type choice
        analysis_type = cms.string(analysis_type),
        # Generic
        verbosity = cms.bool(verbosity),
        print_HLT_event_path = cms.bool(print_HLT_event_path),
        HLT_config_tag = cms.string(HLT_config_tag),
        filter_config_tag = cms.string(filter_config_tag),
        # Triggers
        collect_trigger_stats = cms.bool(collect_trigger_stats),
        ## Single lepton triggers:
        HLT_electron_triggers = cms.vstring(HLT_electron_triggers),
        HLT_muon_triggers = cms.vstring(HLT_muon_triggers),
        ## Dilepton triggers:
        HLT_electron_electron_triggers = cms.vstring(HLT_electron_electron_triggers),
        HLT_electron_muon_triggers = cms.vstring(HLT_electron_muon_triggers),
        HLT_muon_muon_triggers = cms.vstring(HLT_muon_muon_triggers),
        MET_filter_names = cms.vstring(filter_names),
        # Cuts
        min_ele_pT = cms.double(min_ele_pT),
        min_mu_pT = cms.double(min_mu_pT),
        min_veto_ele_pT = cms.double(min_veto_ele_pT),
        min_veto_mu_pT = cms.double(min_veto_mu_pT),
        min_di_ele1_pT = cms.double(min_di_ele1_pT),
        min_di_ele2_pT = cms.double(min_di_ele2_pT),
        min_di_mu1_pT = cms.double(min_di_mu1_pT),
        min_di_mu2_pT = cms.double(min_di_mu2_pT),
        min_jet_pT = cms.double(min_jet_pT),
        min_jet2_pT = cms.double(min_jet2_pT),
        min_bjet_pT = cms.double(min_bjet_pT),
        max_ele_eta = cms.double(max_ele_eta),
        max_mu_eta = cms.double(max_mu_eta),
        max_veto_ele_eta = cms.double(max_veto_ele_eta),
        max_veto_mu_eta = cms.double(max_veto_mu_eta),
        max_di_ele1_eta = cms.double(max_di_ele1_eta),
        max_di_ele2_eta = cms.double(max_di_ele2_eta),
        max_di_mu1_eta = cms.double(max_di_mu1_eta),
        max_di_mu2_eta = cms.double(max_di_mu2_eta),
        max_jet_eta = cms.double(max_jet_eta),
        max_bjet_eta = cms.double(max_bjet_eta),
        min_njets = cms.int32(min_njets),
        min_di_njets = cms.int32(min_di_njets),
        min_nbtags = cms.int32(min_nbtags),
        min_di_nbtags = cms.int32(min_di_nbtags),
        min_di_mll = cms.double(min_di_mll),
        min_di_met = cms.double(min_di_met),
        # MiniAODhelper
        using_real_data = cms.bool(using_real_data),
        dataset = cms.int32(dataset),
        ## available choices '-': none, 'L': loose, 'M': medium, 'T': tight
        b_tag_strength = cms.string(b_tag_strength),
                          
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
            #mets = cms.InputTag("slimmedMETs","","MAOD"),
            genjets = cms.InputTag("slimmedGenJets"),
            pfcand = cms.InputTag("packedPFCandidates"),
            beamspot = cms.InputTag("offlineBeamSpot"),
            #mvaValues = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
            #mvaCategories = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),
            #eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
            #eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp80"),
            eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"),
            #mvaValues = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"),
            #mvaCategories = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Categories"),
            genTtbarId = cms.InputTag("categorizeGenTtbar", "genTtbarId"),
            pileupinfo = cms.InputTag("slimmedAddPileupInfo"),
            lhepprod = cms.InputTag("externalLHEProducer"),
            badchcandfilter = cms.InputTag("BadChargedCandidateFilter"),
            badpfmufilter = cms.InputTag("BadPFMuonFilter"),
            badglobalmutagger = cms.InputTag("badGlobalMuonTaggerMAOD"),
            cloneglobalmutagger = cms.InputTag("cloneGlobalMuonTaggerMAOD")
        )
)
