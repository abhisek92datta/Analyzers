#ifndef CU_ttH_EDA_event_vars_h
#define CU_ttH_EDA_event_vars_h

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include <vector>

/*
 *
 * struct for per-event variables used in analyze(...)
 *
 */

struct CU_ttH_EDA_event_vars {
  double weight; // total event weight (there are partial weights)

  /// Common, run parameters
  int run_nr;
  int event_nr;
  int lumisection_nr;

  bool isdata;

  /// Number of tags per event
  int n_electrons;
  int n_veto_electrons;
  int n_muons;
  int n_veto_muons;
  int n_leptons;
  int n_sl_jets;
  int n_sl_btags;
  int n_di_jets;
  int n_di_btags;
  int n_ttags;
  int n_Htags;

  int n_di_electrons;
  int n_di_muons;
  int n_di_leptons;

  /// Passing-trigger flags
  bool pass_single_e;
  bool pass_single_mu;
  bool pass_double_mu;
  bool pass_double_e;
  bool pass_elemu;

  /// MET Filters
  bool MET_filters;
  bool filterbadChCandidate;
  bool filterbadPFMuon;
  bool badGlobalMuonTagger;
  bool cloneGlobalMuonTagger;

  /// Particle container vectors
  std::vector<pat::Electron> e_with_id;
  std::vector<pat::Electron> e_selected;
  std::vector<pat::Electron> e_veto_selected;
  std::vector<pat::Muon> mu_selected;
  std::vector<pat::Muon> mu_veto_selected;

  std::vector<TLorentzVector> corr_mu_sl;
  std::vector<TLorentzVector> corr_mu_sl_veto;
  std::vector<TLorentzVector> corr_mu_di;

  std::vector<unsigned int> ele_sl_seeds;
  std::vector<unsigned int> ele_sl_veto_seeds;
  std::vector<unsigned int> ele_di_seeds;
  std::vector<unsigned int> mu_sl_seeds;
  std::vector<unsigned int> mu_sl_veto_seeds;
  std::vector<unsigned int> mu_di_seeds;

  std::vector<pat::Electron> e_di_selected;
  std::vector<pat::Muon> mu_di_selected;

  std::vector<pat::Jet> jets_raw;
  std::vector<pat::Jet> jets_raw_puid;
  std::vector<pat::Jet> jets_corrected;

  std::vector<pat::Jet> jets_sl_raw;
  std::vector<pat::Jet> jets_sl_corrected;
  std::vector<pat::Jet> jets_sl_selected_raw;
  std::vector<pat::Jet> jets_sl_selected;
  std::vector<pat::Jet> jets_sl_selected_sorted;
  std::vector<pat::Jet> jets_sl_selected_tag_old;
  std::vector<pat::Jet> jets_sl_selected_tag;
  std::vector<pat::Jet> jets_sl_selected_tag_sorted;

  std::vector<pat::Jet> jets_di_raw;
  std::vector<pat::Jet> jets_di_corrected;
  std::vector<pat::Jet> jets_di_selected_raw;
  std::vector<pat::Jet> jets_di_selected;
  std::vector<pat::Jet> jets_di_selected_sorted;
  std::vector<pat::Jet> jets_di_selected_tag_old;
  std::vector<pat::Jet> jets_di_selected_tag;
  std::vector<pat::Jet> jets_di_selected_tag_sorted;

  std::vector<unsigned int> jet_sl_seeds;
  std::vector<unsigned int> jet_di_seeds;
  std::vector<int> jet_sl_puid;
  std::vector<int> jet_di_puid;
  std::vector<double> jet_sl_pudisc;
  std::vector<double> jet_di_pudisc;
  
  std::vector<pat::Jet> jets_inp_bjetness;

  std::vector<double> vec_jet_pt;
  std::vector<double> vec_jet_eta;
  std::vector<double> vec_jet_csv;
  std::vector<int> vec_jet_hadronFlavour;
  std::vector<int> iSys;

  /// Other quantities
  pat::MET pfMET;
  pat::MET MET_corrected;
  double dimuon_mass;
  double dielectron_mass;
  double dilepton_mass;
  double mll;
  double ht_sl;
  double ht_di;
  int mll_passed;
  int ttHf_cat;
  int ttHFGenFilter;
  int truenpv;
  int SL_tag;
  int DL_tag;
  int FH_tag;

  double MHT;
  double metLD;
  double met_pt_uncorr, met_phi_uncorr;
  double met_pt_corr, met_phi_corr;
  double met_pt_phi_corr, met_phi_phi_corr;
  int met_passed;
  double b_weight_sl, b_weight_sl_lfup, b_weight_sl_hfdown, b_weight_sl_cErr1_down;
  double b_weight_di, b_weight_di_lfup, b_weight_di_hfdown, b_weight_di_cErr1_down;
  double jet1SF_sl;
  double jet1SF_up_sl;
  double jet1SF_down_sl;
  double jet1_jesSF_PileUpDataMC_down_sl;
  double jet1_jesSF_RelativeFSR_up_sl;
  double jet1_jerSF_nominal_sl;
  double jet1SF_di;
  double jet1SF_up_di;
  double jet1SF_down_di;
  double jet1_jesSF_PileUpDataMC_down_di;
  double jet1_jesSF_RelativeFSR_up_di;
  double jet1_jerSF_nominal_di;
  double jet2SF_sl;
  double jet2SF_up_sl;
  double jet2SF_down_sl;
  double jet2_jesSF_PileUpDataMC_down_sl;
  double jet2_jesSF_RelativeFSR_up_sl;
  double jet2_jerSF_nominal_sl;
  double jet2SF_di;
  double jet2SF_up_di;
  double jet2SF_down_di;
  double jet2_jesSF_PileUpDataMC_down_di;
  double jet2_jesSF_RelativeFSR_up_di;
  double jet2_jerSF_nominal_di;

  std::vector<double> lep_sf_id_sl;
  std::vector<double> lep_sf_iso_sl;
  double lep_sf_gsf_sl;
  double lep_sf_hip_sl;
  double lep_sf_trig_sl;
  std::vector<double> lep_sf_id_di;
  std::vector<double> lep_sf_iso_di;
  double lep_sf_gsf_di;
  double lep_sf_hip_di;
  double lep_sf_trig_di;

  double gen_weight;
  double PU_weight;
  double pdf_weight_up;
  double pdf_weight_down;
  double q2_weight_up;
  double q2_weight_down;
  
  reco::Vertex PV;
  
  double bjetnessFV_num_leps; 
  double bjetnessFV_npvTrkOVcollTrk;
  double bjetnessFV_avip3d_val;
  double bjetnessFV_avip3d_sig; 
  double bjetnessFV_avsip3d_sig; 
  double bjetnessFV_avip1d_sig;

  bool is_e;
  bool is_mu;
  bool is_ee;
  bool is_emu;
  bool is_mumu;

  int n_prim_V;
  bool event_selection_SL;
  bool event_selection_DL;
};

#endif
