#ifndef CU_ttH_EDA_event_vars_h
#define CU_ttH_EDA_event_vars_h

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

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

	/// Particle container vectors
	std::vector<pat::Electron> e_with_id;
	std::vector<pat::Electron> e_selected;
	std::vector<pat::Electron> e_veto_selected;
	std::vector<pat::Electron> e_selected_sorted;
	std::vector<pat::Electron> e_veto_selected_sorted;
	std::vector<pat::Muon> mu_selected;
	std::vector<pat::Muon> mu_veto_selected;
	std::vector<pat::Muon> mu_selected_sorted;
	std::vector<pat::Muon> mu_veto_selected_sorted;
	
	std::vector<pat::Electron> e_di_selected;
	std::vector<pat::Muon> mu_di_selected;
	std::vector<pat::Electron> e_di_selected_sorted;
	std::vector<pat::Muon> mu_di_selected_sorted;
	
	std::vector<pat::Jet> jets_raw;
	std::vector<pat::Jet> jets_no_mu;
	std::vector<pat::Jet> jets_no_mu_e;
	
	std::vector<pat::Jet> jets_sl_raw;
	std::vector<pat::Jet> jets_sl_corrected_JEC;
	std::vector<pat::Jet> jets_sl_corrected;
	std::vector<pat::Jet> jets_sl_selected_JEC;
	std::vector<pat::Jet> jets_sl_selected;
	std::vector<pat::Jet> jets_sl_selected_sorted;
	std::vector<pat::Jet> jets_sl_selected_tag_old;
	std::vector<pat::Jet> jets_sl_selected_tag;
	std::vector<pat::Jet> jets_sl_selected_tag_sorted;
	
	std::vector<pat::Jet> jets_di_raw;
	std::vector<pat::Jet> jets_di_corrected_JEC;
	std::vector<pat::Jet> jets_di_corrected;
	std::vector<pat::Jet> jets_di_selected_JEC;
	std::vector<pat::Jet> jets_di_selected;
	std::vector<pat::Jet> jets_di_selected_sorted;
	std::vector<pat::Jet> jets_di_selected_tag_old;
	std::vector<pat::Jet> jets_di_selected_tag;
	std::vector<pat::Jet> jets_di_selected_tag_sorted;
	
	std::vector<double> vec_jet_pt;
  	std::vector<double> vec_jet_eta;
 	std::vector<double> vec_jet_csv;
 	std::vector<int> vec_jet_hadronFlavour;
 	int iSys;

	/// Other quantities
	pat::MET pfMET;
	pat::MET MET_corrected;
	double dimuon_mass;
	double dielectron_mass;
	double dilepton_mass;
	double mll;
	int mll_passed;
	int ttHf_cat;

	double MHT;
	double metLD;
	double met_pt, met_phi;
	int met_passed;
	double b_weight_sl;
	double b_weight_di;
	double jet1SF_sl;
	double jet1SF_up_sl;
	double jet1SF_down_sl;
	double jet1SF_di;
	double jet1SF_up_di;
	double jet1SF_down_di;
	
	double lep_sf_id_sl;
	double lep_sf_iso_sl;
	double lep_sf_trig_sl;
	double lep_sf_id_di;
	double lep_sf_iso_di;
	double lep_sf_trig_di;
	
	double PU_weight;
	double pdf_weight_up;
	double pdf_weight_down;
	double local.q2_weight_up;
	double local.q2_weight_down;

	bool is_e;
	bool is_mu;
	bool is_ee;
	bool is_emu;
	bool is_mumu;
	
	int n_prim_V;
	bool event_selection_SL;
	bool event_selection_DL;
	
	// Gen Particles
	std::vector<reco::GenParticle> genHiggs;
	std::vector<reco::GenParticle> genTops;
	reco::CandidateCollection genHiggs_daughters; // or edm::OwnVector<reco::Candidate>
	reco::CandidateCollection genTop_daughters;
	reco::CandidateCollection genW_daughters;
};

#endif
