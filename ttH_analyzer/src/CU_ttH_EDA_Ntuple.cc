#ifndef CU_ttH_EDA_Ntuple_cc
#define CU_ttH_EDA_Ntuple_cc

#include "Analyzers/ttH_analyzer/interface/CU_ttH_EDA_Ntuple.h"
#include "Analyzers/ttH_analyzer/interface/CU_ttH_EDA_event_vars.h"
#include <cmath>

// Constructor
CU_ttH_EDA_Ntuple::CU_ttH_EDA_Ntuple() {}

// Destructor
CU_ttH_EDA_Ntuple::~CU_ttH_EDA_Ntuple() {}

// Member functions

void CU_ttH_EDA_Ntuple::write_ntuple_SL(const CU_ttH_EDA_event_vars &local,
                                        const MiniAODHelper &miniAODhelper) {
  // Event variales
  nEvent = local.event_nr;
  ls = local.lumisection_nr;
  run = local.run_nr;

  npv = local.truenpv;

  n_ele = local.n_electrons;
  n_mu = local.n_muons;
  n_lep = local.n_leptons;

  fill_ntuple_muons(local.mu_selected, local.corr_mu_sl, miniAODhelper);

  fill_ntuple_electrons(local.e_selected, miniAODhelper);

  fill_ntuple_jets(local.jets_sl_selected_sorted, miniAODhelper);
  fill_ntuple_bjets(local.jets_sl_selected_tag_sorted, miniAODhelper);

  n_jets = local.n_sl_jets;
  n_btags = local.n_sl_btags;

  // MET/MHT
  PFMETpt = local.met_pt_phi_corr; // sqrt(local.pfMET.px()*local.pfMET.px()+local.pfMET.py()*local.pfMET.py());
  PFMETphi = local.met_phi_phi_corr;

  if(!local.isdata)
      fill_ntuple_gen(local);

  fill_SF_SL(local);

  ht = local.ht_sl;
}

void CU_ttH_EDA_Ntuple::write_ntuple_DL(const CU_ttH_EDA_event_vars &local,
                                        const MiniAODHelper &miniAODhelper) {
  // Event variables
  nEvent = local.event_nr;
  ls = local.lumisection_nr;
  run = local.run_nr;

  n_ele = local.n_di_electrons;
  n_mu = local.n_di_muons;
  n_lep = local.n_di_leptons;

  npv = local.truenpv;

  // Muons
  fill_ntuple_muons(local.mu_di_selected, local.corr_mu_di, miniAODhelper);

  // Electrons
  fill_ntuple_electrons(local.e_di_selected, miniAODhelper);

  fill_ntuple_jets(local.jets_di_selected_sorted, miniAODhelper);
  fill_ntuple_bjets(local.jets_di_selected_tag_sorted, miniAODhelper);

  n_jets = local.n_di_jets;
  n_btags = local.n_di_btags;

  // MET/MHT
  PFMETpt = local.met_pt_phi_corr; // sqrt(local.pfMET.px()*local.pfMET.px()+local.pfMET.py()*local.pfMET.py());
  PFMETphi = local.met_phi_phi_corr;

  if(!local.isdata)
    fill_ntuple_gen(local);

  fill_SF_DL(local);

  mll = local.mll;
  ht = local.ht_di;
}


inline void
CU_ttH_EDA_Ntuple::fill_ntuple_muons(const std::vector<pat::Muon> &muons, const std::vector<TLorentzVector> &corr_mu,
                                     const MiniAODHelper &miniAODhelper) {

  if (muons.size() == 1) {
    mu_pt.push_back(corr_mu[0].Pt());
    mu_eta.push_back(muons[0].eta());
    mu_phi.push_back(muons[0].phi());
    mu_E.push_back(muons[0].energy());
    mu_charge.push_back(muons[0].charge());
    mu_iso.push_back(((miniAODhelper.GetMuonRelIso(muons[0], coneSize::R04,
                       corrType::deltaBeta)) * (muons[0].pt()/corr_mu[0].Pt())));
    /*
    mu0_jetNDauChargedMVASel = muons[0].userFloat("nearestJetNDauCharged");
    mu0_miniRelIso = muons[0].userFloat("miniIso");
    mu0_miniIsoCharged = muons[0].userFloat("miniAbsIsoCharged");
    mu0_miniIsoNeutral = muons[0].userFloat("miniAbsIsoNeutralcorr");
    mu0_jetPtRel = muons[0].userFloat("nearestJetPtRel");
    mu0_jetPtRatio = muons[0].userFloat("nearestJetPtRatio");
    mu0_jetCSV = muons[0].userFloat("nearestJetCsv");
    mu0_sip3D = muons[0].userFloat("sip3D");
    mu0_dxy = muons[0].userFloat("dxy");
    mu0_dz = muons[0].userFloat("dz");
    mu0_segmentCompatibility = muons[0].segmentCompatibility();
    mu0_leptonMVA = muons[0].userFloat("leptonMVA");
    mu0_mediumID = muons[0].isMediumMuon();
    mu0_dpt_div_pt =
    muons[0].muonBestTrack()->ptError()/muons[0].muonBestTrack()->pt();
    mu0_isfakeablesel = muons[0].userFloat("idFakeable") > 0.5;
    mu0_iscutsel = muons[0].userFloat("idCutBased") > 0.5;
    mu0_ismvasel = muons[0].userFloat("idMVABased") > 0.5;
    */
  }

  else if (muons.size() == 2) {

    int lead = 0;
    int sublead = 1;

    if(corr_mu[0].Pt() < corr_mu[1].Pt()){
        lead = 1;
        sublead = 0;
    }

    mu_pt.push_back(corr_mu[lead].Pt());
    mu_eta.push_back(muons[lead].eta());
    mu_phi.push_back(muons[lead].phi());
    mu_E.push_back(muons[lead].energy());
    mu_charge.push_back(muons[lead].charge());
    mu_iso.push_back(((miniAODhelper.GetMuonRelIso(muons[lead], coneSize::R04,
                       corrType::deltaBeta)) * (muons[lead].pt()/corr_mu[lead].Pt())));

    mu_pt.push_back(corr_mu[sublead].Pt());
    mu_eta.push_back(muons[sublead].eta());
    mu_phi.push_back(muons[sublead].phi());
    mu_E.push_back(muons[sublead].energy());
    mu_charge.push_back(muons[sublead].charge());
    mu_iso.push_back(((miniAODhelper.GetMuonRelIso(muons[sublead], coneSize::R04,
                        corrType::deltaBeta)) * (muons[sublead].pt()/corr_mu[sublead].Pt())));
    /*
    mu1_jetNDauChargedMVASel = muons[1].userFloat("nearestJetNDauCharged");
    mu1_miniRelIso = muons[1].userFloat("miniIso");
    mu1_miniIsoCharged = muons[1].userFloat("miniAbsIsoCharged");
    mu1_miniIsoNeutral = muons[1].userFloat("miniAbsIsoNeutralcorr");
    mu1_jetPtRel = muons[1].userFloat("nearestJetPtRel");
    mu1_jetPtRatio = muons[1].userFloat("nearestJetPtRatio");
    mu1_jetCSV = muons[1].userFloat("nearestJetCsv");
    mu1_sip3D = muons[1].userFloat("sip3D");
    mu1_dxy = muons[1].userFloat("dxy");
    mu1_dz = muons[1].userFloat("dz");
    mu1_segmentCompatibility = muons[1].segmentCompatibility();
    mu1_leptonMVA = muons[1].userFloat("leptonMVA");
    mu1_mediumID = muons[1].isMediumMuon();
    mu1_dpt_div_pt =
    muons[1].muonBestTrack()->ptError()/muons[1].muonBestTrack()->pt();
    mu1_isfakeablesel = muons[1].userFloat("idFakeable") > 0.5;
    mu1_iscutsel = muons[1].userFloat("idCutBased") > 0.5;
    mu1_ismvasel = muons[1].userFloat("idMVABased") > 0.5;
    */
  }
}

inline void CU_ttH_EDA_Ntuple::fill_ntuple_electrons(
    const std::vector<pat::Electron> &electrons,
    const MiniAODHelper &miniAODhelper) {

  if (electrons.size() == 1) {
    ele_pt.push_back(electrons[0].pt());
    ele_eta.push_back(electrons[0].eta());
    ele_phi.push_back(electrons[0].phi());
    ele_E.push_back(electrons[0].energy());
    ele_charge.push_back(electrons[0].charge());
    ele_iso.push_back(miniAODhelper.GetElectronRelIso(electrons[0], coneSize::R03, corrType::rhoEA, effAreaType::spring15));
    /*
    ele0_jetNDauChargedMVASel = electrons[0].userFloat("nearestJetNDauCharged");
    ele0_miniRelIso = electrons[0].userFloat("miniIso");
    ele0_miniIsoCharged = electrons[0].userFloat("miniAbsIsoCharged");
    ele0_miniIsoNeutral = electrons[0].userFloat("miniAbsIsoNeutralcorr");
    ele0_jetPtRel = electrons[0].userFloat("nearestJetPtRel");
    ele0_jetPtRatio = electrons[0].userFloat("nearestJetPtRatio");
    ele0_jetCSV = electrons[0].userFloat("nearestJetCsv");
    ele0_sip3D = electrons[0].userFloat("sip3D");
    ele0_dxy = electrons[0].userFloat("dxy");
    ele0_dz = electrons[0].userFloat("dz");
    ele0_ntMVAeleID = electrons[0].userFloat("eleMvaId"); // non-triggering mva
    ID
    ele0_leptonMVA = electrons[0].userFloat("leptonMVA");
    ele0_isChargeConsistent = electrons[0].isGsfCtfScPixChargeConsistent();
    ele0_passesConversionVeto = electrons[0].passConversionVeto();
    ele0_nMissingHits = electrons[0].userFloat("numMissingHits");
    ele0_iscutsel = electrons[0].userFloat("idCutBased") > 0.5;
    ele0_ismvasel = electrons[0].userFloat("idMVABased") > 0.5;
    ele0_isfakeablesel = electrons[0].userFloat("idFakeable") > 0.5;
    */
  }

  else if (electrons.size() == 2 ) {

    int lead = 0;
    int sublead = 1;

    if(electrons[0].pt() < electrons[1].pt()){
        lead = 1;
        sublead = 0;
    }

    ele_pt.push_back(electrons[lead].pt());
    ele_eta.push_back(electrons[lead].eta());
    ele_phi.push_back(electrons[lead].phi());
    ele_E.push_back(electrons[lead].energy());
    ele_charge.push_back(electrons[lead].charge());
    ele_iso.push_back(miniAODhelper.GetElectronRelIso(electrons[lead], coneSize::R03, corrType::rhoEA, effAreaType::spring15));
      
    ele_pt.push_back(electrons[sublead].pt());
    ele_eta.push_back(electrons[sublead].eta());
    ele_phi.push_back(electrons[sublead].phi());
    ele_E.push_back(electrons[sublead].energy());
    ele_charge.push_back(electrons[sublead].charge());
    ele_iso.push_back(miniAODhelper.GetElectronRelIso(electrons[sublead], coneSize::R03, corrType::rhoEA, effAreaType::spring15));
    /*
    ele1_jetNDauChargedMVASel = electrons[1].userFloat("nearestJetNDauCharged");
    ele1_miniRelIso = electrons[1].userFloat("miniIso");
    ele1_miniIsoCharged = electrons[1].userFloat("miniAbsIsoCharged");
    ele1_miniIsoNeutral = electrons[1].userFloat("miniAbsIsoNeutralcorr");
    ele1_jetPtRel = electrons[1].userFloat("nearestJetPtRel");
    ele1_jetPtRatio = electrons[1].userFloat("nearestJetPtRatio");
    ele1_jetCSV = electrons[1].userFloat("nearestJetCsv");
    ele1_sip3D = electrons[1].userFloat("sip3D");
    ele1_dxy = electrons[1].userFloat("dxy");
    ele1_dz = electrons[1].userFloat("dz");
    ele1_ntMVAeleID = electrons[1].userFloat("eleMvaId"); // non-triggering mva
    ID
    ele1_leptonMVA = electrons[1].userFloat("leptonMVA");
    ele1_isChargeConsistent = electrons[1].isGsfCtfScPixChargeConsistent();
    ele1_passesConversionVeto = electrons[1].passConversionVeto();
    ele1_nMissingHits = electrons[1].userFloat("numMissingHits");
    ele1_iscutsel = electrons[1].userFloat("idCutBased") > 0.5;
    ele1_ismvasel = electrons[1].userFloat("idMVABased") > 0.5;
    ele1_isfakeablesel = electrons[1].userFloat("idFakeable") > 0.5;
    */
  }
}

inline void
CU_ttH_EDA_Ntuple::fill_ntuple_jets(const std::vector<pat::Jet> &jets,
                                    const MiniAODHelper &miniAODhelper) {

  for(unsigned int i=0; i<jets.size(); i++){
      jet_pt.push_back(jets[i].pt());
      jet_eta.push_back(jets[i].eta());
      jet_phi.push_back(jets[i].phi());
      jet_E.push_back(jets[i].energy());
      jet_CSV.push_back(miniAODhelper.GetJetCSV(jets[i], "pfCombinedInclusiveSecondaryVertexV2BJetTags"));
  }

}

inline void
CU_ttH_EDA_Ntuple::fill_ntuple_bjets(const std::vector<pat::Jet> &bjets,
                                    const MiniAODHelper &miniAODhelper) {

    for(unsigned int i=0; i<bjets.size(); i++){
        bjet_pt.push_back(bjets[i].pt());
        bjet_eta.push_back(bjets[i].eta());
        bjet_phi.push_back(bjets[i].phi());
        bjet_E.push_back(bjets[i].energy());
        bjet_CSV.push_back(miniAODhelper.GetJetCSV(bjets[i], "pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    }
    
}

inline void
CU_ttH_EDA_Ntuple::fill_ntuple_gen(const CU_ttH_EDA_event_vars &local){

    /*
    reco::GenParticle lead_ele;
    reco::GenParticle sublead_ele;
    reco::GenParticle lead_mu;
    reco::GenParticle sublead_mu;
    double temp1, temp2;

    temp1 = temp2 = -9999;
    if(local.genelectrons_selected.size() == 1)
        lead_ele = local.genelectrons_selected[0];
    else {
        for( unsigned int i=0; i<local.genelectrons_selected.size(); i++ ){
            if(local.genelectrons_selected[i].pt() > temp1){
                lead_ele = local.genelectrons_selected[i];
                temp1 = local.genelectrons_selected[i].pt();
            }
        }
        for( unsigned int i=0; i<local.genelectrons_selected.size(); i++ ){
            if( local.genelectrons_selected[i].pt() > temp2  && local.genelectrons_selected[i].pt() < temp1  ){
                sublead_ele = local.genelectrons_selected[i];
                temp2 = local.genelectrons_selected[i].pt();
            }
        }
    }

    temp1 = temp2 = -9999;
    if(local.genmuons_selected.size() == 1)
        lead_mu = local.genmuons_selected[0];
    else {
        for( unsigned int i=0; i<local.genmuons_selected.size(); i++ ){
            if(local.genmuons_selected[i].pt() > temp1){
                lead_mu = local.genmuons_selected[i];
                temp1 = local.genmuons_selected[i].pt();
            }
        }
        for( unsigned int i=0; i<local.genmuons_selected.size(); i++ ){
            if( local.genmuons_selected[i].pt() > temp2  && local.genmuons_selected[i].pt() < temp1  ){
                sublead_mu = local.genmuons_selected[i];
                temp2 = local.genmuons_selected[i].pt();
            }
        }
    }

    if( local.is_e == 1 ){
        if( lead_ele.pt() > 30 && fabs(lead_ele.eta()) < 2.1 ){
            genele_pt.push_back(lead_ele.pt());
            genele_eta.push_back(lead_ele.eta());
            genele_phi.push_back(lead_ele.phi());
            genele_E.push_back(lead_ele.energy());
            genele_charge.push_back(lead_ele.charge());
        }
        for( unsigned int i=0; i<local.genjets_selected.size(); i++ ){
            if(local.genjets_selected[i].pt() > 30){
                genjet_pt.push_back(local.genjets_selected[i].pt());
                genjet_eta.push_back(local.genjets_selected[i].eta());
                genjet_phi.push_back(local.genjets_selected[i].phi());
                genjet_E.push_back(local.genjets_selected[i].energy());
            }
        }
    }
    if( local.is_mu == 1 ){
        if( lead_mu.pt() > 26 && fabs(lead_mu.eta()) < 2.1 ){
            genmu_pt.push_back(lead_mu.pt());
            genmu_eta.push_back(lead_mu.eta());
            genmu_phi.push_back(lead_mu.phi());
            genmu_E.push_back(lead_mu.energy());
            genmu_charge.push_back(lead_mu.charge());
        }
        for( unsigned int i=0; i<local.genjets_selected.size(); i++ ){
            if(local.genjets_selected[i].pt() > 30){
                genjet_pt.push_back(local.genjets_selected[i].pt());
                genjet_eta.push_back(local.genjets_selected[i].eta());
                genjet_phi.push_back(local.genjets_selected[i].phi());
                genjet_E.push_back(local.genjets_selected[i].energy());
            }
        }
    }
    */

    for( unsigned int i=0; i<local.genelectrons_selected.size(); i++ ){
        genele_pt.push_back(local.genelectrons_selected[i].pt());
        genele_eta.push_back(local.genelectrons_selected[i].eta());
        genele_phi.push_back(local.genelectrons_selected[i].phi());
        genele_E.push_back(local.genelectrons_selected[i].energy());
        genele_charge.push_back(local.genelectrons_selected[i].charge());
    }
    for( unsigned int i=0; i<local.genmuons_selected.size(); i++ ){
        genmu_pt.push_back(local.genmuons_selected[i].pt());
        genmu_eta.push_back(local.genmuons_selected[i].eta());
        genmu_phi.push_back(local.genmuons_selected[i].phi());
        genmu_E.push_back(local.genmuons_selected[i].energy());
        genmu_charge.push_back(local.genmuons_selected[i].charge());
    }
    for( unsigned int i=0; i<local.genjets_selected.size(); i++ ){
        genjet_pt.push_back(local.genjets_selected[i].pt());
        genjet_eta.push_back(local.genjets_selected[i].eta());
        genjet_phi.push_back(local.genjets_selected[i].phi());
        genjet_E.push_back(local.genjets_selected[i].energy());
    }

}

inline void
CU_ttH_EDA_Ntuple::fill_SF_SL(const CU_ttH_EDA_event_vars &local) {

    if (local.isdata == 0){
        lep_sf_id.push_back(local.lep_sf_id_sl[0]);
        lep_sf_id.push_back(1.00);
        lep_sf_iso.push_back(local.lep_sf_iso_sl[0]);
        lep_sf_iso.push_back(1.00);
        lep_sf_trig = local.lep_sf_trig_sl;

        b_weight = local.b_weight_sl;
        gen_weight = local.gen_weight;
        PU_weight = local.PU_weight;
        pdf_weight_up = local.pdf_weight_up;
        pdf_weight_down = local.pdf_weight_down;
        q2_weight_up = local.q2_weight_up;
        q2_weight_down = local.q2_weight_down;
    }
    else {
        lep_sf_id.push_back(1.00);
        lep_sf_id.push_back(1.00);
        lep_sf_iso.push_back(1.00);
        lep_sf_iso.push_back(1.00);
        lep_sf_trig = 1;

        b_weight = 1;
        gen_weight = 1;
        PU_weight = 1;
        pdf_weight_up = 1;
        pdf_weight_down = 1;
        q2_weight_up = 1;
        q2_weight_down = 1;
    }
}

inline void
CU_ttH_EDA_Ntuple::fill_SF_DL(const CU_ttH_EDA_event_vars &local) {

    if (local.isdata == 0){
        lep_sf_id.push_back(local.lep_sf_id_di[0]);
        lep_sf_id.push_back(local.lep_sf_id_di[1]);
        lep_sf_iso.push_back(local.lep_sf_iso_di[0]);
        lep_sf_iso.push_back(local.lep_sf_iso_di[1]);
        lep_sf_trig = local.lep_sf_trig_di;

        b_weight = local.b_weight_di;
        gen_weight = local.gen_weight;
        PU_weight = local.PU_weight;
        pdf_weight_up = local.pdf_weight_up;
        pdf_weight_down = local.pdf_weight_down;
        q2_weight_up = local.q2_weight_up;
        q2_weight_down = local.q2_weight_down;
    }
    else {
        lep_sf_id.push_back(1.00);
        lep_sf_id.push_back(1.00);
        lep_sf_iso.push_back(1.00);
        lep_sf_iso.push_back(1.00);
        lep_sf_trig = 1;

        b_weight = 1;
        gen_weight = 1;
        PU_weight = 1;
        pdf_weight_up = 1;
        pdf_weight_down = 1;
        q2_weight_up = 1;
        q2_weight_down = 1;
    }
}

void CU_ttH_EDA_Ntuple::initialize() {
  // event variables

  nEvent = -9999;
  ls = -9999;
  run = -9999;

  npv = -9999;
  ttHf_cat = -9999;
  ttHFGenFilter = -9999;
  SL_tag = -9999;
  DL_tag = -9999;
  FH_tag = -9999;

  n_ele = -9999;
  n_mu = -9999;
  n_lep = -9999;

  mll = -9999;
  ht = -9999;

  // muons
  mu_pt.clear();
  mu_eta.clear();
  mu_phi.clear();
  mu_E.clear();
  mu_charge.clear();
  mu_iso.clear();
  /*
  mu0_jetNDauChargedMVASel = -9999;
  mu0_miniRelIso = -9999.;
  mu0_miniIsoCharged = -9999.;
  mu0_miniIsoNeutral = -9999.;
  mu0_jetPtRel = -9999.;
  mu0_jetPtRatio = -9999.;
  mu0_jetCSV = -9999.;
  mu0_sip3D = -9999.;
  mu0_dxy = -9999.;
  mu0_dz = -9999.;
  mu0_segmentCompatibility = -9999.;
  mu0_leptonMVA = -9999.;
  mu0_mediumID = -9999.;
  mu0_dpt_div_pt = -9999.;
  mu0_iscutsel = -9999;
  mu0_ismvasel = -9999;
  mu0_isfakeablesel = -9999;
  mu1_jetNDauChargedMVASel = -9999;
  mu1_miniRelIso = -9999.;
  mu1_miniIsoCharged = -9999.;
  mu1_miniIsoNeutral = -9999.;
  mu1_jetPtRel = -9999.;
  mu1_jetPtRatio = -9999.;
  mu1_jetCSV = -9999.;
  mu1_sip3D = -9999.;
  mu1_dxy = -9999.;
  mu1_dz = -9999.;
  mu1_segmentCompatibility = -9999.;
  mu1_leptonMVA = -9999.;
  mu1_mediumID = -9999.;
  mu1_dpt_div_pt = -9999.;
  mu1_iscutsel = -9999;
  mu1_ismvasel = -9999;
  mu1_isfakeablesel = -9999;
  */

  // electrons
  ele_pt.clear();
  ele_eta.clear();
  ele_phi.clear();
  ele_E.clear();
  ele_charge.clear();
  ele_iso.clear();
  /*
  ele0_jetNDauChargedMVASel = -9999;
  ele0_miniRelIso = -9999.;
  ele0_miniIsoCharged = -9999.;
  ele0_miniIsoNeutral = -9999.;
  ele0_jetPtRel = -9999.;
  ele0_jetPtRatio = -9999.;
  ele0_jetCSV = -9999.;
  ele0_sip3D = -9999.;
  ele0_dxy = -9999.;
  ele0_dz = -9999.;
  ele0_ntMVAeleID = -9999.;
  ele0_leptonMVA = -9999.;
  ele0_isChargeConsistent = -9999;
  ele0_passesConversionVeto = -9999;
  ele0_nMissingHits = -9999;
  ele0_iscutsel = -9999;
  ele0_ismvasel = -9999;
  ele0_isfakeablesel = -9999;
  ele1_jetNDauChargedMVASel = -9999;
  ele1_miniRelIso = -9999.;
  ele1_miniIsoCharged = -9999.;
  ele1_miniIsoNeutral = -9999.;
  ele1_jetPtRel = -9999.;
  ele1_jetPtRatio = -9999.;
  ele1_jetCSV = -9999.;
  ele1_sip3D = -9999.;
  ele1_dxy = -9999.;
  ele1_dz = -9999.;
  ele1_ntMVAeleID = -9999.;
  ele1_leptonMVA = -9999.;
  ele1_isChargeConsistent = -9999;
  ele1_passesConversionVeto = -9999;
  ele1_nMissingHits = -9999;
  ele1_iscutsel = -9999;
  ele1_ismvasel = -9999;
  ele1_isfakeablesel = -9999;
  */

  n_jets = -9999;
  n_btags = -9999;

  // jets
  jet_pt.clear();
  jet_eta.clear();
  jet_phi.clear();
  jet_E.clear();
  jet_CSV.clear();

  bjet_pt.clear();
  bjet_eta.clear();
  bjet_phi.clear();
  bjet_E.clear();
  bjet_CSV.clear();

  // MET
  PFMETpt = -9999.;
  PFMETphi = -9999.;

  // Gen-Level Info

  genmu_pt.clear();
  genmu_eta.clear();
  genmu_phi.clear();
  genmu_E.clear();
  genmu_charge.clear();

  genele_pt.clear();
  genele_eta.clear();
  genele_phi.clear();
  genele_E.clear();
  genele_charge.clear();

  genjet_pt.clear();
  genjet_eta.clear();
  genjet_phi.clear();
  genjet_E.clear();

  // SF and event_weight

  lep_sf_id.clear();
  lep_sf_iso.clear();
  lep_sf_trig = -9999;

  b_weight = -9999;
  gen_weight = -9999;
  PU_weight = -9999;
  pdf_weight_up = -9999;
  pdf_weight_down = -9999;
  q2_weight_up = -9999;
  q2_weight_down = -9999;

}


void CU_ttH_EDA_Ntuple::set_up_branches(TTree *tree) {
  // initialize ntuple
  this->initialize();

  // set up tree branches
  tree->Branch("nEvent", &nEvent);
  tree->Branch("ls", &ls);
  tree->Branch("run", &run);

  tree->Branch("npv", &npv);
  tree->Branch("ttHf_cat", &ttHf_cat);
  tree->Branch("ttHFGenFilter", &ttHFGenFilter);
  tree->Branch("SL_tag", &SL_tag);
  tree->Branch("DL_tag", &DL_tag);
  tree->Branch("FH_tag", &FH_tag);


  tree->Branch("n_ele", &n_ele);
  tree->Branch("n_mu", &n_mu);
  tree->Branch("n_lep", &n_lep);

  tree->Branch("mll", &mll);
  tree->Branch("ht", &ht);

  // muons
  tree->Branch("mu_pt", &mu_pt);
  tree->Branch("mu_eta", &mu_eta);
  tree->Branch("mu_phi", &mu_phi);
  tree->Branch("mu_E", &mu_E);
  tree->Branch("mu_charge", &mu_charge);
  tree->Branch("mu_iso", &mu_iso);
  
  
  /*tree->Branch("mu0_jetNDauChargedMVASel", &mu0_jetNDauChargedMVASel);
  tree->Branch("mu0_miniRelIso",           &mu0_miniRelIso);
  tree->Branch("mu0_miniIsoCharged",       &mu0_miniIsoCharged);
  tree->Branch("mu0_miniIsoNeutral",       &mu0_miniIsoNeutral);
  tree->Branch("mu0_jetPtRel",             &mu0_jetPtRel);
  tree->Branch("mu0_jetPtRatio",           &mu0_jetPtRatio);
  tree->Branch("mu0_jetCSV",               &mu0_jetCSV);
  tree->Branch("mu0_sip3D",                &mu0_sip3D);
  tree->Branch("mu0_dxy",                  &mu0_dxy);
  tree->Branch("mu0_dz",                   &mu0_dz);
  tree->Branch("mu0_segmentCompatibility", &mu0_segmentCompatibility);
  tree->Branch("mu0_leptonMVA",            &mu0_leptonMVA);
  tree->Branch("mu0_mediumID", &mu0_mediumID);
  tree->Branch("mu0_dpt_div_pt", &mu0_dpt_div_pt);
  tree->Branch("mu0_iscutsel", &mu0_iscutsel);
  tree->Branch("mu0_ismvasel", &mu0_ismvasel);
  tree->Branch("mu0_isfakeablesel", &mu0_isfakeablesel);
  tree->Branch("mu1_jetNDauChargedMVASel", &mu1_jetNDauChargedMVASel);
  tree->Branch("mu1_miniRelIso",           &mu1_miniRelIso);
  tree->Branch("mu1_miniIsoCharged",       &mu1_miniIsoCharged);
  tree->Branch("mu1_miniIsoNeutral",       &mu1_miniIsoNeutral);
  tree->Branch("mu1_jetPtRel",             &mu1_jetPtRel);
  tree->Branch("mu1_jetPtRatio",           &mu1_jetPtRatio);
  tree->Branch("mu1_jetCSV",               &mu1_jetCSV);
  tree->Branch("mu1_sip3D",                &mu1_sip3D);
  tree->Branch("mu1_dxy",                  &mu1_dxy);
  tree->Branch("mu1_dz",                   &mu1_dz);
  tree->Branch("mu1_segmentCompatibility", &mu1_segmentCompatibility);
  tree->Branch("mu1_leptonMVA",            &mu1_leptonMVA);
  tree->Branch("mu1_mediumID", &mu1_mediumID);
  tree->Branch("mu1_dpt_div_pt", &mu1_dpt_div_pt);
  tree->Branch("mu1_iscutsel", &mu1_iscutsel);
  tree->Branch("mu1_ismvasel", &mu1_ismvasel);
  tree->Branch("mu1_isfakeablesel", &mu1_isfakeablesel);
  */
  // electrons
  
  tree->Branch("ele_pt", &ele_pt);
  tree->Branch("ele_eta", &ele_eta);
  tree->Branch("ele_phi", &ele_phi);
  tree->Branch("ele_E", &ele_E);
  tree->Branch("ele_charge", &ele_charge);
  tree->Branch("ele_iso", &ele_iso);
  
  /*
  tree->Branch("ele0_jetNDauChargedMVASel", &ele0_jetNDauChargedMVASel);
  tree->Branch("ele0_miniRelIso",           &ele0_miniRelIso);
  tree->Branch("ele0_miniIsoCharged",       &ele0_miniIsoCharged);
  tree->Branch("ele0_miniIsoNeutral",       &ele0_miniIsoNeutral);
  tree->Branch("ele0_jetPtRel",             &ele0_jetPtRel);
  tree->Branch("ele0_jetPtRatio",           &ele0_jetPtRatio);
  tree->Branch("ele0_jetCSV",               &ele0_jetCSV);
  tree->Branch("ele0_sip3D",                &ele0_sip3D);
  tree->Branch("ele0_dxy",                  &ele0_dxy);
  tree->Branch("ele0_dz",                   &ele0_dz);
  tree->Branch("ele0_ntMVAeleID",           &ele0_ntMVAeleID);
  tree->Branch("ele0_leptonMVA",            &ele0_leptonMVA);
  tree->Branch("ele0_isChargeConsistent", &ele0_isChargeConsistent);
  tree->Branch("ele0_passesConversionVeto", &ele0_passesConversionVeto);
  tree->Branch("ele0_nMissingHits", &ele0_nMissingHits);
  tree->Branch("ele0_iscutsel", &ele0_iscutsel);
  tree->Branch("ele0_ismvasel", &ele0_ismvasel);
  tree->Branch("ele0_isfakeablesel", &ele0_isfakeablesel);
  tree->Branch("ele1_jetNDauChargedMVASel", &ele1_jetNDauChargedMVASel);
  tree->Branch("ele1_miniRelIso",           &ele1_miniRelIso);
  tree->Branch("ele1_miniIsoCharged",       &ele1_miniIsoCharged);
  tree->Branch("ele1_miniIsoNeutral",       &ele1_miniIsoNeutral);
  tree->Branch("ele1_jetPtRel",             &ele1_jetPtRel);
  tree->Branch("ele1_jetPtRatio",           &ele1_jetPtRatio);
  tree->Branch("ele1_jetCSV",               &ele1_jetCSV);
  tree->Branch("ele1_sip3D",                &ele1_sip3D);
  tree->Branch("ele1_dxy",                  &ele1_dxy);
  tree->Branch("ele1_dz",                   &ele1_dz);
  tree->Branch("ele1_ntMVAeleID",           &ele1_ntMVAeleID);
  tree->Branch("ele1_leptonMVA",            &ele1_leptonMVA);
  tree->Branch("ele1_isChargeConsistent", &ele1_isChargeConsistent);
  tree->Branch("ele1_passesConversionVeto", &ele1_passesConversionVeto);
  tree->Branch("ele1_nMissingHits", &ele1_nMissingHits);
  tree->Branch("ele1_iscutsel", &ele1_iscutsel);
  tree->Branch("ele1_ismvasel", &ele1_ismvasel);
  tree->Branch("ele1_isfakeablesel", &ele1_isfakeablesel);
  */
  
  tree->Branch("n_jets", &n_jets);
  tree->Branch("n_btags", &n_btags);
  // jets
  tree->Branch("jet_pt", &jet_pt);
  tree->Branch("jet_eta", &jet_eta);
  tree->Branch("jet_phi", &jet_phi);
  tree->Branch("jet_E", &jet_E);
  tree->Branch("jet_CSV", &jet_CSV);

  tree->Branch("bjet_pt", &bjet_pt);
  tree->Branch("bjet_eta", &bjet_eta);
  tree->Branch("bjet_phi", &bjet_phi);
  tree->Branch("bjet_E", &bjet_E);
  tree->Branch("bjet_CSV", &bjet_CSV);

  // MET
  tree->Branch("PFMETpt", &PFMETpt);
  tree->Branch("PFMETphi", &PFMETphi);

  // Gen-Level Info
  tree->Branch("genmu_pt", &genmu_pt);
  tree->Branch("genmu_eta", &genmu_eta);
  tree->Branch("genmu_phi", &genmu_phi);
  tree->Branch("genmu_E", &genmu_E);
  tree->Branch("genmu_charge", &genmu_charge);

  tree->Branch("genele_pt", &genele_pt);
  tree->Branch("genele_eta", &genele_eta);
  tree->Branch("genele_phi", &genele_phi);
  tree->Branch("genele_E", &genele_E);
  tree->Branch("genele_charge", &genele_charge);

  tree->Branch("genjet_pt", &genjet_pt);
  tree->Branch("genjet_eta", &genjet_eta);
  tree->Branch("genjet_phi", &genjet_phi);
  tree->Branch("genjet_E", &genjet_E);

  // SF and event-weight
  tree->Branch("lep_sf_id", &lep_sf_id);
  tree->Branch("lep_sf_iso", &lep_sf_iso);
  tree->Branch("lep_sf_trig", &lep_sf_trig);

  tree->Branch("b_weight", &b_weight);
  tree->Branch("gen_weight", &gen_weight);
  tree->Branch("PU_weight", &PU_weight);
  tree->Branch("pdf_weight_up", &pdf_weight_up);
  tree->Branch("pdf_weight_down", &pdf_weight_down);
  tree->Branch("q2_weight_up", &q2_weight_up);
  tree->Branch("q2_weight_down", &q2_weight_down);

}

#endif
