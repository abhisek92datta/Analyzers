#ifndef CU_ttH_EDA_cc
#define CU_ttH_EDA_cc

// -*- C++ -*-
//
// Package:    ttH-LeptonPlusJets/CU_ttH_EDA
// Class:      CU_ttH_EDA
//
/**\class CU_ttH_EDA CU_ttH_EDA.cc
 ttH-LeptonPlusJets/AnalysisCode/plugins/CU_ttH_EDA.cc

 Description: [one line class summary]

 Implementation:
		 [Notes on implementation]
*/

/// Includes
#include "CU_ttH_EDA.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

/*
 * Function/method section
*/

/// Constructor
CU_ttH_EDA::CU_ttH_EDA(const edm::ParameterSet &iConfig):
	// Analysis type
	config_analysis_type (iConfig.getParameter<string>("analysis_type")),
	// Generic
	verbose_ (iConfig.getParameter<bool>("verbosity")),
	dumpHLT_ (iConfig.getParameter<bool>("print_HLT_event_path")),
	hltTag (iConfig.getParameter<string>("HLT_config_tag")),
	filterTag (iConfig.getParameter<string>("filter_config_tag")),
	// Triggers
	trigger_stats (iConfig.getParameter<bool>("collect_trigger_stats")),
	trigger_on_HLT_e (iConfig.getParameter<std::vector<string>>("HLT_electron_triggers")),
	trigger_on_HLT_mu (iConfig.getParameter<std::vector<string>>("HLT_muon_triggers")),
	trigger_on_HLT_ee (iConfig.getParameter<std::vector<string>>("HLT_electron_electron_triggers")),
	trigger_on_HLT_emu (iConfig.getParameter<std::vector<string>>("HLT_electron_muon_triggers")),
	trigger_on_HLT_mumu (iConfig.getParameter<std::vector<string>>("HLT_muon_muon_triggers")),
	// Cuts
	min_ele_pT (iConfig.getParameter<double>("min_ele_pT")),
	min_mu_pT (iConfig.getParameter<double>("min_mu_pT")),
	min_veto_ele_pT (iConfig.getParameter<double>("min_veto_ele_pT")),
	min_veto_mu_pT (iConfig.getParameter<double>("min_veto_mu_pT")),
	min_di_ele1_pT (iConfig.getParameter<double>("min_di_ele1_pT")),
	min_di_ele2_pT (iConfig.getParameter<double>("min_di_ele2_pT")),
	min_di_mu1_pT (iConfig.getParameter<double>("min_di_mu1_pT")),
	min_di_mu2_pT (iConfig.getParameter<double>("min_di_mu2_pT")),
	min_jet_pT (iConfig.getParameter<double>("min_jet_pT")),
	min_jet2_pT (iConfig.getParameter<double>("min_jet2_pT")),
	min_bjet_pT (iConfig.getParameter<double>("min_bjet_pT")),
	max_ele_eta (iConfig.getParameter<double>("max_ele_eta")),
	max_mu_eta (iConfig.getParameter<double>("max_mu_eta")),
	max_veto_ele_eta (iConfig.getParameter<double>("max_veto_ele_eta")),
	max_veto_mu_eta (iConfig.getParameter<double>("max_veto_mu_eta")),
	max_di_ele1_eta (iConfig.getParameter<double>("max_di_ele1_eta")),
	max_di_ele2_eta (iConfig.getParameter<double>("max_di_ele2_eta")),
	max_di_mu1_eta (iConfig.getParameter<double>("max_di_mu1_eta")),
	max_di_mu2_eta (iConfig.getParameter<double>("max_di_mu2_eta")),
	max_jet_eta (iConfig.getParameter<double>("max_jet_eta")),
	max_bjet_eta (iConfig.getParameter<double>("max_bjet_eta")),
	min_njets (iConfig.getParameter<int>("min_njets")),
	min_di_njets (iConfig.getParameter<int>("min_di_njets")),
	min_nbtags (iConfig.getParameter<int>("min_nbtags")),
	min_di_nbtags (iConfig.getParameter<int>("min_di_nbtags")),
	min_di_mll (iConfig.getParameter<double>("min_di_mll")),
	min_di_met (iConfig.getParameter<double>("min_di_met")),
	// miniAODhelper
	isdata (iConfig.getParameter<bool>("using_real_data")),
	MAODHelper_b_tag_strength (iConfig.getParameter<string>("b_tag_strength")[0])
{
	/*
	 * now do whatever initialization is needed
	*/

	/// temporary mock-up parameters
	MAODHelper_era = "2015_74x";
	MAODHelper_sample_nr = 2500;

	total_xs = 831.76;
	sample_n = 25446993;
	int_lumi = 10000;
	weight_sample = int_lumi * total_xs / sample_n;

	Load_configuration_set_type(config_analysis_type);
	Load_configuration_MAODH(isdata);

	Set_up_tokens(iConfig.getParameter<edm::ParameterSet>("input_tags"));
	Set_up_histograms();
	Set_up_output_files();

	Set_up_Tree();
	
	Set_up_b_weights();
	
	init_PU_weight();
	init_PDF_weight();
}

/// Destructor
CU_ttH_EDA::~CU_ttH_EDA()
{
	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

	Close_output_files();
	
	delete f_CSVwgt_HF;
	delete f_CSVwgt_LF;
}

// ------------ method called for each event  ------------
void CU_ttH_EDA::analyze(const edm::Event &iEvent,
						 const edm::EventSetup &iSetup)
{
	using namespace edm;
	++event_count;

	/// Declaring local struct for data readout and manipulations
	CU_ttH_EDA_event_vars local;
	
	/// Triggers have not fired yet. Check_triggers, Check_filters will adjust
	local.pass_single_e = false;
	local.pass_single_mu = false;
	local.pass_double_mu = false;
	local.pass_double_e = false;
	local.pass_elemu = false;
	Update_common_vars(iEvent, local);
	
	//if (local.event_nr != 1410845 )
	//	return;
	
	
	/// Create and set up edm:Handles in stack mem.
	edm_Handles handle;
	Set_up_handles(iEvent, handle, token);
	
	// for PDF weight
	if(!isdata)
		iEvent.getByToken(token.event_gen_info, handle.event_gen_info);
	
	// for Q2 weight
	if(!isdata)
		iEvent.getByToken( token.lheptoken, handle.EvtHandle );
	
	/// Run checks on event containers via their handles
	Check_triggers(handle.triggerResults, local);
	//Check_filters(handle.filterResults);
	Check_vertices_set_MAODhelper(handle.vertices);

	local.n_prim_V = Check_PV(handle.vertices);

	// Setting rho
	auto rho = handle.srcRho;
	miniAODhelper.SetRho(*rho);

	// weight_gen = event_gen_info.product()->weight();
	//local.weight = weight_sample * (handle.event_gen_info.product()->weight());
	local.weight = 0;

	if (trigger_stats) {
		h_hlt->Fill(0., 1);
		h_flt->Fill(0., 1);
	}

	/// Lepton selection
	local.e_with_id = miniAODhelper.GetElectronsWithMVAid(handle.electrons_for_mva, handle.mvaValues, handle.mvaCategories);
	// Single Lepton
	local.mu_selected = miniAODhelper.GetSelectedMuons(
		*(handle.muons), min_mu_pT, muonID::muonTight, coneSize::R04, corrType::deltaBeta, max_mu_eta);
	local.mu_veto_selected = miniAODhelper.GetSelectedMuons(
		*(handle.muons), min_veto_mu_pT, muonID::muonTightDL, coneSize::R04, corrType::deltaBeta, max_veto_mu_eta);
	local.e_selected = miniAODhelper.GetSelectedElectrons(
		local.e_with_id, min_ele_pT, electronID::electronEndOf15MVA80iso0p15, max_ele_eta);
	local.e_veto_selected = miniAODhelper.GetSelectedElectrons(
		local.e_with_id, min_veto_ele_pT, electronID::electronEndOf15MVA80iso0p15, max_veto_ele_eta);
	local.n_electrons = static_cast<int>(local.e_selected.size());
	local.n_veto_electrons = static_cast<int>(local.e_veto_selected.size());
	local.n_muons = static_cast<int>(local.mu_selected.size());
	local.n_veto_muons = static_cast<int>(local.mu_veto_selected.size());
	local.n_leptons = local.n_electrons + local.n_muons;
	
	// Dilepton
	local.mu_di_selected = miniAODhelper.GetSelectedMuons(
		*(handle.muons), min_di_mu2_pT, muonID::muonTightDL, coneSize::R04, corrType::deltaBeta, max_di_mu2_eta);
	local.e_di_selected = miniAODhelper.GetSelectedElectrons(
		local.e_with_id, min_di_ele2_pT, electronID::electronEndOf15MVA80iso0p15, max_di_ele2_eta);
	local.n_di_electrons = static_cast<int>(local.e_di_selected.size());
	local.n_di_muons = static_cast<int>(local.mu_di_selected.size());
	/// Sort leptons by pT
	local.mu_di_selected_sorted = miniAODhelper.GetSortedByPt(local.mu_di_selected);
	local.e_di_selected_sorted = miniAODhelper.GetSortedByPt(local.e_di_selected);
	local.n_di_leptons = local.n_di_electrons + local.n_di_muons;
	

	/// Jet selection
	
	// Uncorrected jets
	local.jets_raw = miniAODhelper.GetUncorrectedJets(*(handle.jets));
	// Jet Energy Correction
	SetFactorizedJetCorrector();
	local.jets_sl_corrected_JEC = GetCorrectedJets_JEC(local.jets_raw, *rho);
	local.jets_di_corrected_JEC = GetCorrectedJets_JEC(local.jets_raw, *rho);
	if(isdata) {
		local.jets_sl_corrected = local.jets_sl_corrected_JEC;
		local.jets_di_corrected = local.jets_di_corrected_JEC;
	}
	else {
		local.jets_sl_corrected = GetCorrectedJets_JER(local.jets_sl_corrected_JEC, *rho);
		local.jets_di_corrected = GetCorrectedJets_JER(local.jets_di_corrected_JEC, *rho);
	}
	// ID selection
	local.jets_sl_corrected = CheckJetID(local.jets_sl_corrected, *(handle.jets));
	local.jets_di_corrected = CheckJetID(local.jets_di_corrected, *(handle.jets));
	// Remove Overlap
	local.jets_sl_corrected = miniAODhelper.GetDeltaRCleanedJets(local.jets_sl_corrected, local.mu_veto_selected, local.e_veto_selected, 0.4);
	local.jets_di_corrected = miniAODhelper.GetDeltaRCleanedJets(local.jets_di_corrected, local.mu_di_selected, local.e_di_selected, 0.4);

	// for b-weight
	local.iSys = 0; // none - 0,  JESUp - 7 , JESDown - 8		
	
	// Jet selection
	index.clear();
	i=0;
	for ( auto& jet : local.jets_sl_corrected) {
		if ( (jet.pt() > min_jet_pT ) && (fabs(jet.eta()) < max_jet_eta ) ) {
			local.jets_sl_selected.push_back(jet);	
			index.push_back(i);
		}
		i++;
	}
	if(index.size()!=0) {
		for (int j : index) {
			local.jets_sl_selected_JEC.push_back(local.jets_sl_corrected_JEC[j]);
		}
	}
	index.clear();
	i=0;
	for ( auto& jet : local.jets_di_corrected) {
		if ( (jet.pt() > min_jet2_pT ) && (fabs(jet.eta()) < max_jet_eta ) ) {
			local.jets_di_selected.push_back(jet);	
			index.push_back(i);
		}
		i++;
	}
	if(index.size()!=0) {
		for (int j : index) {
			local.jets_di_selected_JEC.push_back(local.jets_di_corrected_JEC[j]);
		}
	}

	// b-tagged jet selection
	local.jets_sl_selected_tag = miniAODhelper.GetSelectedJets(
		local.jets_sl_selected, min_bjet_pT, max_bjet_eta, jetID::none,
		MAODHelper_b_tag_strength);
	local.jets_di_selected_tag = miniAODhelper.GetSelectedJets(
		local.jets_di_selected, min_bjet_pT, max_bjet_eta, jetID::none,
		MAODHelper_b_tag_strength);
	
	/*
	for (const auto& jet : local.jets_sl_selected_tag_old) {
		if (miniAODhelper.GetJetCSV(jet,"pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.89) {
			local.jets_sl_selected_tag.push_back(jet);
		}
	}
	*/
	
	local.n_sl_jets = static_cast<int>(local.jets_sl_selected.size());
	local.n_sl_btags = static_cast<int>(local.jets_sl_selected_tag.size());
	local.n_di_jets = static_cast<int>(local.jets_di_selected.size());
	local.n_di_btags = static_cast<int>(local.jets_di_selected_tag.size());
	
	/// Sort jets by pT
	local.jets_sl_selected_sorted =
		miniAODhelper.GetSortedByPt(local.jets_sl_selected);
	local.jets_sl_selected_tag_sorted =
		miniAODhelper.GetSortedByPt(local.jets_sl_selected_tag);
	local.jets_sl_selected_JEC =
		miniAODhelper.GetSortedByPt(local.jets_sl_selected_JEC);
		
	local.jets_di_selected_sorted =
		miniAODhelper.GetSortedByPt(local.jets_di_selected);
	local.jets_di_selected_tag_sorted =
		miniAODhelper.GetSortedByPt(local.jets_di_selected_tag);
	local.jets_di_selected_JEC =
		miniAODhelper.GetSortedByPt(local.jets_di_selected_JEC);

	/// Top and Higgs tagging using collections through handles. adjusts
	/// local.<tag>
	//Top_tagger(handle.top_jets, local);
	//Higgs_tagger(handle.subfilter_jets, local);

	/// MET
	local.pfMET = handle.METs->front();
	local.met_pt = local.pfMET.pt();
        local.met_phi = local.pfMET.phi();
	local.met_passed = 0;
	local.mll_passed = 0;

	local.is_e = false;
	local.is_mu = false;
	local.is_ee = false;
	local.is_emu = false;
	local.is_mumu = false;

	// flag for determining whether to select an event for writing
	local.event_selection_SL = false;
	local.event_selection_DL = false;

	/*
	if(!isdata) {
		local.pass_single_e = 1;
		local.pass_single_mu = 1;
		local.pass_double_e = 1;
		local.pass_double_mu = 1;
		local.pass_elemu = 1;
	}
	*/
	
	local.PU_weight = -1;
	local.pdf_weight_up = -1;
	local.pdf_weight_down = -1;
	local.q2_weight_up = -1;
	local.q2_weight_down = -1;
	
	// Event selection criteria for single lepton events
	Check_SL_Event_Selection(local);
	
	// Event selection criteria for dilepton events
	Check_DL_Event_Selection(local);
	
	if (local.event_selection_SL!=0 || local.event_selection_DL!=0){
		selection_count++;
	}
	
	
	/// Check tags, fill hists, print events
	if (local.event_selection_SL!=0) {
		Fill_addn_quant(local, *rho, handle);
		Check_Fill_Print_single_lepton(local);
	}
	else if (local.event_selection_DL!=0) {
		Fill_addn_quant(local, *rho, handle);
		Check_Fill_Print_di_lepton(local);
	}

	eventTree->Fill();
}

// ------------ method called once each job just before starting event loop
// ------------
void CU_ttH_EDA::beginJob()
{
	TH1::SetDefaultSumw2(true);

	event_count = 0;
	selection_count = 0;
	sl_e = sl_mu = dl_ee = dl_emu = dl_mumu = 0;
}

// ------------ method called once each job just after ending the event loop
// ------------
void CU_ttH_EDA::endJob() { return; }

// ------------ method called when starting to processes a run  ------------
void CU_ttH_EDA::beginRun(const edm::Run &iRun, const edm::EventSetup &iSetup)
{
	/// Update HLTConfigProvider(s) for the new run
	bool hlt_config_changed = true; // init() updates this one
	if (!hlt_config.init(iRun, iSetup, hltTag, hlt_config_changed)) {
		std::cerr << "Warning, didn't find trigger process HLT,\t" << hltTag
				  << std::endl;
		return;
	}

	if (hlt_config_changed)
		std::cout << "New " << hltTag << " config has been loaded.\n";
	
	/*
	bool filter_config_changed = true; // init() updates this one
	if (!filter_config.init(iRun, iSetup, filterTag, filter_config_changed)) {
		std::cerr << "Warning, didn't find filter process HLT,\t" << filterTag
				  << std::endl;
		return;
	}
	
	if (filter_config_changed)
		std::cout << "New " << filterTag << " config has been loaded.\n";
	*/
	/// Set up filter and trigger name vectors and maps
	if (trigger_stats) {
		Set_up_trigger_name_vectors();

		if (Set_up_Run_histograms_triggers() != 0) {
			std::cerr << "Setting up histograms for trigger/filter counts has "
					  << "failed\n";
			return;
		}
	}
	
}

// ------------ method called when ending the processing of a run  ------------
void CU_ttH_EDA::endRun(const edm::Run &, const edm::EventSetup &)
{
	// report results of sync exercises
	if (analysis_type == Analyze_lepton_jet) {
	
	}

	if (analysis_type == Analyze_dilepton) {
	
	}

	if (trigger_stats)
		End_Run_hist_fill_triggers();

	// report on triggers fired
	if (trigger_stats && dumpHLT_) {
		std::cout
			<< "***************************************************************"
			<< std::endl;
		std::cout << "  Summary for HLT: Total number of events = "
				  << event_count << std::endl;
		for (std::map<std::string, unsigned long>::const_iterator iter =
				 n_trigger_fired.begin();
			 iter != n_trigger_fired.end(); ++iter) {
			std::string name = iter->first;
			double eff = 100 * double(iter->second) / double(event_count);
			printf("\t %s,\t %lu,\t %.1f \n", name.c_str(), iter->second, eff);
		}
		std::cout
			<< "***************************************************************"
			<< std::endl;
		std::cout << "  Summary for Filters: Total number of events = "
				  << event_count << std::endl;
		for (std::map<std::string, unsigned long>::const_iterator iter =
				 n_filter_fired.begin();
			 iter != n_filter_fired.end(); ++iter) {
			std::string name = iter->first;
			double eff = 100 * double(iter->second) / double(event_count);
			printf("\t %s,\t %lu,\t %.1f \n", name.c_str(), iter->second, eff);
		}
		std::cout
			<< "***************************************************************"
			<< std::endl;
	}


	std::cout
		<< "***************************************************************"
		<< std::endl;
	std::cout <<"SL e = "<<sl_e<< std::endl;
	std::cout <<"SL mu = "<<sl_mu<< std::endl;
	std::cout <<"DL ee = "<<dl_ee<< std::endl;
	std::cout <<"DL emu = "<<dl_emu<< std::endl;
	std::cout <<"DL mumu = "<<dl_mumu<< std::endl;
	std::cout
		<< "***************************************************************"
		<< std::endl;
	std::cout
		<< "***************************************************************"
		<< std::endl;
	std::cout <<"Number of events selected = "<<selection_count<< std::endl;
	std::cout
		<< "***************************************************************"
		<< std::endl;
	std::cout
		<< "***************************************************************"
		<< std::endl;
	std::cout << "  Total number of events = " << event_count << std::endl;
	std::cout
		<< "***************************************************************"
		<< std::endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the
// module  ------------
void CU_ttH_EDA::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
	// The following says we do not know what parameters are allowed so do no
	// validation
	// Please change this to state exactly what you do use, even if it is no
	// parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

// define this as a CMSSW plugin
DEFINE_FWK_MODULE(CU_ttH_EDA);

#endif
