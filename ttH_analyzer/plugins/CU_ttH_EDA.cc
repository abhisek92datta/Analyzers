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
CU_ttH_EDA::CU_ttH_EDA(const edm::ParameterSet &iConfig)
    : // Analysis type
      config_analysis_type(iConfig.getParameter<string>("analysis_type")),
      // Generic
      verbose_(iConfig.getParameter<bool>("verbosity")),
      dumpHLT_(iConfig.getParameter<bool>("print_HLT_event_path")),
      hltTag(iConfig.getParameter<string>("HLT_config_tag")),
      filterTag(iConfig.getParameter<string>("filter_config_tag")),
      // Triggers
      trigger_stats(iConfig.getParameter<bool>("collect_trigger_stats")),
      trigger_on_HLT_e(
          iConfig.getParameter<std::vector<string>>("HLT_electron_triggers")),
      trigger_on_HLT_mu(
          iConfig.getParameter<std::vector<string>>("HLT_muon_triggers")),
      trigger_on_HLT_ee(iConfig.getParameter<std::vector<string>>(
          "HLT_electron_electron_triggers")),
      trigger_on_HLT_emu(iConfig.getParameter<std::vector<string>>(
          "HLT_electron_muon_triggers")),
      trigger_on_HLT_mumu(
          iConfig.getParameter<std::vector<string>>("HLT_muon_muon_triggers")),
      MET_filter_names(
          iConfig.getParameter<std::vector<string>>("MET_filter_names")),
      // Cuts
      min_ele_pT(iConfig.getParameter<double>("min_ele_pT")),
      min_mu_pT(iConfig.getParameter<double>("min_mu_pT")),
      min_veto_ele_pT(iConfig.getParameter<double>("min_veto_ele_pT")),
      min_veto_mu_pT(iConfig.getParameter<double>("min_veto_mu_pT")),
      min_di_ele1_pT(iConfig.getParameter<double>("min_di_ele1_pT")),
      min_di_ele2_pT(iConfig.getParameter<double>("min_di_ele2_pT")),
      min_di_mu1_pT(iConfig.getParameter<double>("min_di_mu1_pT")),
      min_di_mu2_pT(iConfig.getParameter<double>("min_di_mu2_pT")),
      min_jet_pT(iConfig.getParameter<double>("min_jet_pT")),
      min_jet2_pT(iConfig.getParameter<double>("min_jet2_pT")),
      min_bjet_pT(iConfig.getParameter<double>("min_bjet_pT")),
      max_ele_eta(iConfig.getParameter<double>("max_ele_eta")),
      max_mu_eta(iConfig.getParameter<double>("max_mu_eta")),
      max_veto_ele_eta(iConfig.getParameter<double>("max_veto_ele_eta")),
      max_veto_mu_eta(iConfig.getParameter<double>("max_veto_mu_eta")),
      max_di_ele1_eta(iConfig.getParameter<double>("max_di_ele1_eta")),
      max_di_ele2_eta(iConfig.getParameter<double>("max_di_ele2_eta")),
      max_di_mu1_eta(iConfig.getParameter<double>("max_di_mu1_eta")),
      max_di_mu2_eta(iConfig.getParameter<double>("max_di_mu2_eta")),
      max_jet_eta(iConfig.getParameter<double>("max_jet_eta")),
      max_bjet_eta(iConfig.getParameter<double>("max_bjet_eta")),
      min_njets(iConfig.getParameter<int>("min_njets")),
      min_di_njets(iConfig.getParameter<int>("min_di_njets")),
      min_nbtags(iConfig.getParameter<int>("min_nbtags")),
      min_di_nbtags(iConfig.getParameter<int>("min_di_nbtags")),
      min_di_mll(iConfig.getParameter<double>("min_di_mll")),
      min_di_met(iConfig.getParameter<double>("min_di_met")),
      // miniAODhelper
      isdata(iConfig.getParameter<bool>("using_real_data")),
      dataset(iConfig.getParameter<int>("dataset")),
      write_csv(iConfig.getParameter<bool>("write_csv")),
      is_OLS(iConfig.getParameter<bool>("is_OLS")),
      MAODHelper_b_tag_strength(
          iConfig.getParameter<string>("b_tag_strength")[0])
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
    //Set_up_histograms();
    if(write_csv)
        Set_up_output_files();

    // For Jet Correction
    SetFactorizedJetCorrector();
    // miniAODhelper.SetJetCorrectorUncertainty();

    SetpT_ResFile();

    Set_up_Tree();
    Set_up_weights();

    // Rochester Correction
    rc.init("data/rcdata.2016.v3");
    //rc = new RoccoR("data/rcdata.2016.v3");

    //r = new TRandom3(0);
}

/// Destructor
CU_ttH_EDA::~CU_ttH_EDA()
{
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    /*
    delete _jetCorrector_MC;
    delete _jetCorrector_BCD;
    delete _jetCorrector_EF;
    delete _jetCorrector_G;
    delete _jetCorrector_H;
    delete _jetCorrectorUnc;
    */
    
    rnd.SetSeed(0);
    if(write_csv)
        Close_output_files();

    //delete r;
    //delete events_combined;
    //delete eventTree;

    f_CSVwgt_HF->Close();
    f_CSVwgt_LF->Close();
    delete f_CSVwgt_HF;
    delete f_CSVwgt_LF;
    delete pdfSet;
    //delete h_hlt;
    //delete h_flt;
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

    local.isdata = isdata;


    //if (local.event_nr != 315602 && local.event_nr != 315735)
    //	return;
    //std::cout<<local.event_nr<<"    "<<event_count<<"\n\n";


    /// Create and set up edm:Handles in stack mem.
    edm_Handles handle;
    Set_up_handles(iEvent, iSetup, handle, token, isdata);

    // for JEC
    /*
    const JetCorrector *corrector =
        JetCorrector::getJetCorrector("ak4PFchsL1L2L3", iSetup);
    miniAODhelper.SetJetCorrector(corrector);
    miniAODhelper.SetJetCorrectorUncertainty(iSetup);
    */

    // for JER (from GT)
    //resolution = JME::JetResolution::get(iSetup, "AK4PFchs_pt");
    //resolution_sf = JME::JetResolutionScaleFactor::get(iSetup, "AK4PFchs");

    // for PDF weight
    if (!isdata)
        iEvent.getByToken(token.event_gen_info, handle.event_gen_info);
    // for Q2 weight
    if (!isdata && !is_OLS)
        iEvent.getByToken(token.lheptoken, handle.EvtHandle);
    // for ttHf categorization
    if (!isdata) {
        iEvent.getByToken(token.genTtbarIdToken_, handle.genTtbarId);
        iEvent.getByToken(token.ttHFGenFilterToken_, handle.ttHFGenFilter);
    }

    /// Run checks on event containers via their handles
    Check_triggers(handle.triggerResults, local);
    Check_filters(handle.filterResults, local);


    /// MET filters
    local.filterbadChCandidate = *handle.ifilterbadChCand;
    local.filterbadPFMuon = *handle.ifilterbadPFMuon;
    local.badGlobalMuonTagger = *handle.ibadGlobalMuonTagger;
    local.cloneGlobalMuonTagger = *handle.icloneGlobalMuonTagger;

    Check_vertices_set_MAODhelper(handle.vertices);
   
    local.n_prim_V = Check_PV(handle.vertices);
    prim_vertex = *(handle.vertices->begin());

    // Setting rho
    auto rho = handle.srcRho;
    miniAODhelper.SetRho(*rho);

    // weight_gen = event_gen_info.product()->weight();
    // local.weight = weight_sample *
    // (handle.event_gen_info.product()->weight());
    local.weight = 0;

	/*
    if (trigger_stats) {
        h_hlt->Fill(0., 1);
        h_flt->Fill(0., 1);
    }
	*/
	
    /// Lepton selection
    Select_Leptons(local, handle);

    /// Jet selection
    //Select_Jets(local, iEvent, iSetup, handle, *rho, resolution);
    Select_Jets(local, iEvent, iSetup, handle, *rho);

    /// MET
    Init_Mets(local, handle);

    init_flags(local);
    init_weights(local);
    init_SFs(local);
    init_bjetness_var(local);

    // Event selection criteria for single lepton events
    Check_SL_Event_Selection(local);

    // Event selection criteria for dilepton events
    Check_DL_Event_Selection(local);

    if (local.event_selection_SL != 0 || local.event_selection_DL != 0) {
        ++selection_count;
    }

    /*
    std::cout<<"\n"<<local.event_nr<<"  "<<local.pass_single_mu<<"  "<<local.MET_filters<<"  "<<local.filterbadChCandidate<<"  "<<local.filterbadPFMuon<<"  "<<local.badGlobalMuonTagger<<"  "<<local.cloneGlobalMuonTagger<<"\n";
    std::cout<<local.n_muons<<"  "<<local.n_veto_electrons<<"  "<<local.n_veto_muons<<"\n";
    std::cout<<local.n_sl_jets<<"  "<<local.n_sl_btags<<"\n";
    std::cout<<"\n";
    for(unsigned int i=0; i<local.corr_mu_sl.size(); i++)
        std::cout<<local.corr_mu_sl[i].Pt()<<"  "<<local.corr_mu_sl[i].Eta()<<"\n";
    std::cout<<"\n";


    for(unsigned int i=0; i<local.jets_raw.size(); i++)
        std::cout<<local.jets_raw[i].pt()<<" "<<miniAODhelper.GetJetCSV(local.jets_raw[i],"pfCombinedInclusiveSecondaryVertexV2BJetTags")<<"\n";
    std::cout<<"\n";

    for(unsigned int i=0; i<local.jets_raw_puid.size(); i++)
        std::cout<<local.jets_raw_puid[i].pt()<<" "<<miniAODhelper.GetJetCSV(local.jets_raw_puid[i],"pfCombinedInclusiveSecondaryVertexV2BJetTags")<<"\n";
    std::cout<<"\n";

    for(unsigned int i=0; i<local.jets_sl_raw.size(); i++)
        std::cout<<local.jets_sl_raw[i].pt()<<" "<<miniAODhelper.GetJetCSV(local.jets_sl_raw[i],"pfCombinedInclusiveSecondaryVertexV2BJetTags")<<"\n";
    std::cout<<"\n";

    for(unsigned int i=0; i<local.jets_sl_corrected.size(); i++)
        std::cout<<local.jets_sl_corrected[i].pt()<<" "<<local.jets_sl_corrected[i].eta()<<"  "<<
        miniAODhelper.GetJetCSV(local.jets_sl_corrected[i],"pfCombinedInclusiveSecondaryVertexV2BJetTags")<<"\n";
    std::cout<<"\n";

    for(int i=0; i<local.n_sl_jets; i++)
        std::cout<<local.jets_sl_selected_sorted[i].pt()<<" "<<miniAODhelper.GetJetCSV(local.jets_sl_selected_sorted[i],"pfCombinedInclusiveSecondaryVertexV2BJetTags")<<"\n";
    std::cout<<"\n";

    for(int i=0; i<local.n_sl_btags; i++)
        std::cout<<local.jets_sl_selected_tag_sorted[i].pt()<<"  "<<miniAODhelper.GetJetCSV(local.jets_sl_selected_tag_sorted[i],"pfCombinedInclusiveSecondaryVertexV2BJetTags")<<"\n";
     */


    // Initialize Ntuple and filling some general info
    hbbNtuple.initialize();

	// Event variables
	hbbNtuple.nEvent = local.event_nr;
	hbbNtuple.ls = local.lumisection_nr;
	hbbNtuple.run = local.run_nr;

    // generator weight
    local.gen_weight = 1;
    if (!isdata)
        local.gen_weight = handle.event_gen_info->weight();
    hbbNtuple.gen_weight = local.gen_weight;

    //ttHF categorization
    local.ttHf_cat = -1;
    if (!isdata && handle.genTtbarId.isValid())
        local.ttHf_cat = *handle.genTtbarId;
    hbbNtuple.ttHf_cat = local.ttHf_cat;

    // ttHFGenFilter
    local.ttHFGenFilter = false;
    if(!isdata) {
        local.ttHFGenFilter = *(handle.ttHFGenFilter);
		if(local.ttHFGenFilter == true)
			hbbNtuple.ttHFGenFilter = 1;
		else
			hbbNtuple.ttHFGenFilter = 0;
    }
	else
		hbbNtuple.ttHFGenFilter = -1;

    //SL, DL and FH tagger
    if (!isdata)
        Lepton_Tag(*(handle.genparticles), local);
    else {
        local.SL_tag = -1;
        local.DL_tag = -1;
        local.FH_tag = -1;
    }
    hbbNtuple.SL_tag = local.SL_tag;
    hbbNtuple.DL_tag = local.DL_tag;
    hbbNtuple.FH_tag = local.FH_tag;
    gen_SL_count = gen_SL_count + local.SL_tag;
    gen_DL_count = gen_DL_count + local.DL_tag;
    gen_FH_count = gen_FH_count + local.FH_tag;
	if(local.ttHFGenFilter){
		gen_SL_hf_count = gen_SL_hf_count + local.SL_tag;
		gen_DL_hf_count = gen_DL_hf_count + local.DL_tag;
		gen_FH_hf_count = gen_FH_hf_count + local.FH_tag;
	}
	else{
		gen_SL_nonhf_count = gen_SL_nonhf_count + local.SL_tag;
		gen_DL_nonhf_count = gen_DL_nonhf_count + local.DL_tag;
		gen_FH_nonhf_count = gen_FH_nonhf_count + local.FH_tag;
	}
	gen_tot_hf_count = gen_SL_hf_count + gen_DL_hf_count + gen_FH_hf_count;
	gen_tot_nonhf_count = gen_SL_nonhf_count + gen_DL_nonhf_count + gen_FH_nonhf_count;
    gen_tot_count = gen_SL_count + gen_DL_count + gen_FH_count;

    // Generator Level b-quark info
    if (!isdata){
        Fill_Gen_b_info(*(handle.genparticles), local);
        hbbNtuple.fill_ntuple_gen_b(local);
    }


    /// Check tags, fill ntuple, print events
    if (local.event_selection_SL != 0) {
        Fill_addn_quant(local, iEvent, iSetup, *rho, handle);
        if(write_csv)
            Check_Fill_Print_single_lepton(local);
        hbbNtuple.write_ntuple_SL(local, miniAODhelper);
    } else if (local.event_selection_DL != 0) {
        Fill_addn_quant(local, iEvent, iSetup, *rho, handle);
        if(write_csv)
            Check_Fill_Print_di_lepton(local);
        hbbNtuple.write_ntuple_DL(local, miniAODhelper);
    }

    // Fill Ntuple
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
    gen_SL_count = gen_DL_count = gen_FH_count = gen_tot_count = gen_SL_hf_count = gen_SL_nonhf_count = gen_DL_hf_count = gen_DL_nonhf_count = gen_FH_hf_count = gen_FH_nonhf_count = gen_tot_hf_count = gen_tot_nonhf_count = 0;
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

    bool filter_config_changed = true; // init() updates this one
    if (!filter_config.init(iRun, iSetup, filterTag, filter_config_changed)) {
        std::cerr << "Warning, didn't find filter process HLT,\t" << filterTag
                  << std::endl;
        return;
    }

    if (filter_config_changed)
        std::cout << "New " << filterTag << " config has been loaded.\n";
    
    /// Set up filter and trigger name vectors and maps
    /*
    if (trigger_stats) {
        Set_up_trigger_name_vectors();

        if (Set_up_Run_histograms_triggers() != 0) {
            std::cerr << "Setting up histograms for trigger/filter counts has "
                      << "failed\n";
            return;
        }
    }
    */
}

// ------------ method called when ending the processing of a run  ------------
void CU_ttH_EDA::endRun(const edm::Run &, const edm::EventSetup &)
{
    // report results of sync exercises
    /*
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
	*/
    std::cout
        << "***************************************************************"
        << std::endl;
    std::cout << "SL e = " << sl_e << std::endl;
    std::cout << "SL mu = " << sl_mu << std::endl;
    std::cout << "DL ee = " << dl_ee << std::endl;
    std::cout << "DL emu = " << dl_emu << std::endl;
    std::cout << "DL mumu = " << dl_mumu << std::endl;
    std::cout
        << "***************************************************************"
        << std::endl;
    std::cout
        << "***************************************************************"
        << std::endl;
    std::cout << "Number of events selected = " << selection_count << std::endl;
    std::cout
        << "***************************************************************"
        << std::endl;
    std::cout
        << "***************************************************************"
        << std::endl;
    std::cout << "Total number of events = " << event_count << std::endl;
    std::cout
        << "***************************************************************"
        << std::endl;
    std::cout
        << "***************************************************************"
        << std::endl;
    std::cout << "Total number of SL tagged events = " << gen_SL_count << ", SL hf : "<< gen_SL_hf_count<<", SL nonhf : "<< gen_SL_nonhf_count<< std::endl;
    std::cout << "Total number of DL tagged events = " << gen_DL_count << ", DL hf : "<< gen_DL_hf_count<<", DL nonhf : "<< gen_DL_nonhf_count<< std::endl;
    std::cout << "Total number of FH tagged events = " << gen_FH_count << ", FH hf : "<< gen_FH_hf_count<<", FH nonhf : "<< gen_FH_nonhf_count<< std::endl;
    std::cout << "Total number of tagged events = " << gen_tot_count << ", Total hf : "<< gen_tot_hf_count<<", Total nonhf : "<< gen_tot_nonhf_count<< std::endl;
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
