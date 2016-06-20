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
        # Cuts
        min_tight_lepton_pT = cms.double(min_tight_lepton_pT),
        min_ele_pT = cms.double(min_ele_pT),
        min_mu_pT = cms.double(min_mu_pT),
        min_veto_ele_pT = cms.double(min_veto_ele_pT),
        min_veto_mu_pT = cms.double(min_veto_mu_pT),
        #min_tau_pT = cms.double(min_tau_pT),
        min_jet_pT = cms.double(min_jet_pT),
        min_bjet_pT = cms.double(min_bjet_pT),
        max_ele_eta = cms.double(max_ele_eta),
        max_mu_eta = cms.double(max_mu_eta),
        max_veto_ele_eta = cms.double(max_veto_ele_eta),
        max_veto_mu_eta = cms.double(max_veto_mu_eta),
        max_jet_eta = cms.double(max_jet_eta),
        max_bjet_eta = cms.double(max_bjet_eta),
        min_njets = cms.int32(min_njets),
        min_nbtags = cms.int32(min_nbtags),
        # Jets
        jet_corrector = cms.InputTag(jet_corrector),
        # MiniAODhelper
        using_real_data = cms.bool(using_real_data),
        ## available choices '-': none, 'L': loose, 'M': medium, 'T': tight
        b_tag_strength = cms.string(b_tag_strength),
                          
        # InputTags
        input_tags = cms.PSet(
            pv = cms.InputTag("offlineSlimmedPrimaryVertices"),
            sv = cms.InputTag("slimmedSecondaryVertices"),
            pileup = cms.InputTag("addPileupInfo"),
            rho = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
            #"fixedGridRhoFastjetAll"
            electrons = cms.InputTag("slimmedElectrons"),
            muons = cms.InputTag("slimmedMuons"),
            #taus = cms.InputTag("slimmedTaus"),
            jets = cms.InputTag("slimmedJets"),
            mets = cms.InputTag("slimmedMETs"),
            pfcand = cms.InputTag("packedPFCandidates"),
            beamspot = cms.InputTag("offlineBeamSpot"),
            packedgen = cms.InputTag("packedGenParticles"),
            prunedgen = cms.InputTag("prunedGenParticles"),
            mvaValues = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Values"),
            mvaCategories = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Categories")
            #mvaValues = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig25nsV1Values",""),
            #mvaCategories = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig25nsV1Categories","")
        )
)
