#ifndef CU_ttH_EDA_Setups_cc
#define CU_ttH_EDA_Setups_cc

/// Includes
#include "CU_ttH_EDA.h"

void CU_ttH_EDA::init_PU_weight()
{
	fin.open("ttH_analyzer/data/PU_weights.txt");
	for(int i=0; i<50; i++) {
		fin>>PU_x[i]>>PU_y[i];
	}
	fin.close();
}

void CU_ttH_EDA::init_PDF_weight()
{
	CT14nlo_PDFSet = new LHAPDF::PDFSet("CT14nlo");
	_systPDFs = CT14nlo_PDFSet->mkPDFs();
	//return CT14nlo_PDFSet;
}

void CU_ttH_EDA::Close_output_files()
{
	fclose(events_combined);
	//if (analysis_type == Analyze_lepton_jet) {
		//fclose(events_single_lepton);
		/*
		fclose(events_e_cut1);
		fclose(events_e_cut2);
		fclose(events_e_cut3);
		fclose(events_e_cut4);
		fclose(events_e_cut5);
		fclose(events_e_cut6);
		fclose(events_e_cut7);

		fclose(events_mu_cut1);
		fclose(events_mu_cut2);
		fclose(events_mu_cut3);
		fclose(events_mu_cut4);
		fclose(events_mu_cut5);
		fclose(events_mu_cut6);
		fclose(events_mu_cut7);
		*/
	//}

	//if (analysis_type == Analyze_dilepton) {
		//fclose(events_di_lepton);
		/*
		fclose(events_dimu_cut1);
		fclose(events_dimu_cut2);
		fclose(events_dimu_cut3);
		fclose(events_dimu_cut4);
		fclose(events_dimu_cut5);
		fclose(events_dimu_cut6);
		fclose(events_dimu_cut7);

		fclose(events_diele_cut1);
		fclose(events_diele_cut2);
		fclose(events_diele_cut3);
		fclose(events_diele_cut4);
		fclose(events_diele_cut5);
		fclose(events_diele_cut6);
		fclose(events_diele_cut7);

		fclose(events_elemu_cut1);
		fclose(events_elemu_cut2);
		fclose(events_elemu_cut3);
		fclose(events_elemu_cut4);
		fclose(events_elemu_cut5);
		*/
	//}
}

void CU_ttH_EDA::Set_up_histograms()
{
	// 	h_electron_selection = fs_->make<TH1D>("h_electron_selection",
	// ";electron cut", 12, 0 , 12 );
	// 	h_muon_selection = fs_->make<TH1D>("h_muon_selection", ";muon cut", 12,
	// 0 , 12 );
	//
	// 	h_electron_selection->GetXaxis()->SetBinLabel(1, "All");
	// 	h_electron_selection->GetXaxis()->SetBinLabel(2, "p_{T}>20, |#eta|<2.4,
	// !inCrack");
	// 	h_electron_selection->GetXaxis()->SetBinLabel(3, "full5x5 #sigma_{i#eta
	// i#eta}");
	// 	h_electron_selection->GetXaxis()->SetBinLabel(4, "|#Delta #eta_{in}|");
	// 	h_electron_selection->GetXaxis()->SetBinLabel(5, "|#Delta #phi_{in}|");
	// 	h_electron_selection->GetXaxis()->SetBinLabel(6, "hOverE");
	// 	h_electron_selection->GetXaxis()->SetBinLabel(7, "1/E - 1/p");
	// 	h_electron_selection->GetXaxis()->SetBinLabel(8, "d0");
	// 	h_electron_selection->GetXaxis()->SetBinLabel(9, "dZ");
	// 	h_electron_selection->GetXaxis()->SetBinLabel(10,
	// "expectedMissingInnerHits");
	// 	h_electron_selection->GetXaxis()->SetBinLabel(11, "passConversionVeto");
	// 	h_electron_selection->GetXaxis()->SetBinLabel(12, "relIso (#Delta Beta,
	// 0.3)");
	//
	// 	h_muon_selection->GetXaxis()->SetBinLabel(1, "All");
	// 	h_muon_selection->GetXaxis()->SetBinLabel(2, "p_{T}>20, |#eta|<2.4");
	// 	h_muon_selection->GetXaxis()->SetBinLabel(3, "GlobalMuon ||
	// TrackerMuon");
	// 	h_muon_selection->GetXaxis()->SetBinLabel(4, "PFMuon");
	// 	h_muon_selection->GetXaxis()->SetBinLabel(5, "#Chi^{2}");
	// 	h_muon_selection->GetXaxis()->SetBinLabel(6, "validMuonHit");
	// 	h_muon_selection->GetXaxis()->SetBinLabel(7, "validPixelHit");
	// 	h_muon_selection->GetXaxis()->SetBinLabel(8, "trk layers w/meas");
	// 	h_muon_selection->GetXaxis()->SetBinLabel(9, "matched stations");
	// 	h_muon_selection->GetXaxis()->SetBinLabel(10, "d0");
	// 	h_muon_selection->GetXaxis()->SetBinLabel(11, "dZ");
	// 	h_muon_selection->GetXaxis()->SetBinLabel(12, "relIso < 0.1");

	if (analysis_type == Analyze_lepton_jet) {
		/*
		h_tth_syncex1_ele =
			fs_->make<TH1D>("h_tth_syncex1_ele", ";cut", 8, 0, 8);
		h_tth_syncex1_ele->GetXaxis()->SetBinLabel(1, "All events");
		h_tth_syncex1_ele->GetXaxis()->SetBinLabel(2, "Single ele trig");
		h_tth_syncex1_ele->GetXaxis()->SetBinLabel(3, "==1 electron");
		h_tth_syncex1_ele->GetXaxis()->SetBinLabel(4, "==0 muons");
		h_tth_syncex1_ele->GetXaxis()->SetBinLabel(5, ">=4 jets");
		h_tth_syncex1_ele->GetXaxis()->SetBinLabel(6, ">=2 b-tags");
		h_tth_syncex1_ele->GetXaxis()->SetBinLabel(7, ">=1 top-tags");
		h_tth_syncex1_ele->GetXaxis()->SetBinLabel(8, ">=1 Higgs-tags");

		h_tth_syncex1_mu = fs_->make<TH1D>("h_tth_syncex1_mu", ";cut", 8, 0, 8);
		h_tth_syncex1_mu->GetXaxis()->SetBinLabel(1, "All events");
		h_tth_syncex1_mu->GetXaxis()->SetBinLabel(2, "Single mu trig");
		h_tth_syncex1_mu->GetXaxis()->SetBinLabel(3, "==1 muon");
		h_tth_syncex1_mu->GetXaxis()->SetBinLabel(4, "==0 electrons");
		h_tth_syncex1_mu->GetXaxis()->SetBinLabel(5, ">=4 jets");
		h_tth_syncex1_mu->GetXaxis()->SetBinLabel(6, ">=2 b-tags");
		h_tth_syncex1_mu->GetXaxis()->SetBinLabel(7, ">=1 top-tags");
		h_tth_syncex1_mu->GetXaxis()->SetBinLabel(8, ">=1 Higgs-tags");
		*/
	}

	if (analysis_type == Analyze_dilepton) {
		/*
		h_tth_syncex1_dimu =
			fs_->make<TH1D>("h_tth_syncex1_dimu", ";cut", 8, 0, 8);
		h_tth_syncex1_dimu->GetXaxis()->SetBinLabel(1, "All events");
		h_tth_syncex1_dimu->GetXaxis()->SetBinLabel(2, "Double mu trig");
		h_tth_syncex1_dimu->GetXaxis()->SetBinLabel(3, ">=2 muons");
		h_tth_syncex1_dimu->GetXaxis()->SetBinLabel(4, "Mll > 20");
		h_tth_syncex1_dimu->GetXaxis()->SetBinLabel(5, "Z Veto   ");
		h_tth_syncex1_dimu->GetXaxis()->SetBinLabel(6, ">=2 jets");
		h_tth_syncex1_dimu->GetXaxis()->SetBinLabel(7, "MET > 40");
		h_tth_syncex1_dimu->GetXaxis()->SetBinLabel(8, ">=1 b-tags");

		h_tth_syncex1_diele =
			fs_->make<TH1D>("h_tth_syncex1_diele", ";cut", 8, 0, 8);
		h_tth_syncex1_diele->GetXaxis()->SetBinLabel(1, "All events");
		h_tth_syncex1_diele->GetXaxis()->SetBinLabel(2, "Double ele trig");
		h_tth_syncex1_diele->GetXaxis()->SetBinLabel(3, ">=2 electrons");
		h_tth_syncex1_diele->GetXaxis()->SetBinLabel(4, "Mll > 20");
		h_tth_syncex1_diele->GetXaxis()->SetBinLabel(5, "Z Veto   ");
		h_tth_syncex1_diele->GetXaxis()->SetBinLabel(6, ">=2 jets");
		h_tth_syncex1_diele->GetXaxis()->SetBinLabel(7, "MET > 40");
		h_tth_syncex1_diele->GetXaxis()->SetBinLabel(8, ">=1 b-tags");

		h_tth_syncex1_elemu =
			fs_->make<TH1D>("h_tth_syncex1_elemu", ";cut", 6, 0, 6);
		h_tth_syncex1_elemu->GetXaxis()->SetBinLabel(1, "All events");
		h_tth_syncex1_elemu->GetXaxis()->SetBinLabel(2, "Ele-mu trig");
		h_tth_syncex1_elemu->GetXaxis()->SetBinLabel(3, ">=2 leptons");
		h_tth_syncex1_elemu->GetXaxis()->SetBinLabel(4, "Mll > 20");
		h_tth_syncex1_elemu->GetXaxis()->SetBinLabel(5, ">=2 jets");
		h_tth_syncex1_elemu->GetXaxis()->SetBinLabel(6, ">=1 b-tags");
		*/
	}

}

/// Prepare histograms for trigger/filter counts
int CU_ttH_EDA::Set_up_Run_histograms_triggers()
{
	unsigned int numHLT = trigger_names_no_ver.size();
	h_hlt = fs_->make<TH1D>("h_hlt", ";HLT path", numHLT, 0, numHLT);
	if (!h_hlt)
		return 1;

	TAxis *axis = h_hlt->GetXaxis();
	if (!axis)
		return 1;

	for (unsigned int i = 0; i < numHLT; ++i)
		axis->SetBinLabel(i + 1, trigger_names_no_ver[i].c_str());

	unsigned int numFLT = filter_names_no_ver.size();
	h_flt = fs_->make<TH1D>("h_flt", ";Filter path", numFLT, 0, numFLT);
	if (!h_flt)
		return 1;

	axis = h_flt->GetXaxis();
	if (!axis)
		return 1;

	for (unsigned int i = 0; i < numFLT; ++i)
		axis->SetBinLabel(i + 1, filter_names_no_ver[i].c_str());

	return 0;
}

void CU_ttH_EDA::Set_up_trigger_name_vectors()
{
	/// Fill trigger name vectors and counters
	trigger_names = hlt_config.triggerNames();

	trigger_names_no_ver.clear();
	trigger_names_no_ver.push_back("All");
	std::string prefix = "HLT_";
	for (unsigned int i = 0; i < trigger_names.size(); ++i) {
		std::string pathNameNoVer = hlt_config.removeVersion(trigger_names[i]);

		n_trigger_fired[pathNameNoVer] = 0;

		if (trigger_names[i].compare(0, prefix.length(), prefix) == 0)
			trigger_names_no_ver.push_back(pathNameNoVer);
	}

	/// Fill filter name vectors and counters
	filter_names = filter_config.triggerNames();

	filter_names_no_ver.clear();
	filter_names_no_ver.push_back("All");
	for (unsigned int i = 0; i < filter_names.size(); ++i) {
		std::string pathNameNoVer =
			filter_config.removeVersion(filter_names[i]);

		n_filter_fired[pathNameNoVer] = 0;

		filter_names_no_ver.push_back(pathNameNoVer);
	}
}

/// Make and open write-out files
void CU_ttH_EDA::Set_up_output_files()
{
	events_combined = fopen("CU_tth_80X.csv", "w");
	fprintf(events_combined, "run,lumi,event,is_e,is_mu,is_ee,is_emu,is_mumu,n_jets,n_btags,lep1_pt,lep1_iso,lep1_pdgId,lep2_pt,lep2_iso,lep2_pdgId,jet1_pt,jet2_pt,jet1_CSVv2,jet2_CSVv2,jet1_JecSF,jet1_JecSF_up,jet1_JecSF_down,MET_pt,MET_phi,mll,ttHFCategory,PUWeight,bWeight,triggerSF,lepSF,lepISOSF,Q2_upup,Q2_downdown,pdf_up,pdf_down\n");
	/*
	if (analysis_type == Analyze_lepton_jet) {
		events_single_lepton = fopen("CU_events_single_lepton_cuts.csv", "w");
		fprintf(events_single_lepton, "run, lumi, event, is_SL,	is_DL, lep1_pt, lep1_eta, lep1_phi, lep1_iso, lep1_pdgId, lep2_pt, lep2_eta, lep2_phi, lep2_iso, lep2_pdgId, mll,  mll_passed, jet1_pt,	jet2_pt, jet3_pt, jet4_pt, jet1_CSVv2, jet2_CSVv2, jet3_CSVv2, jet4_CSVv2, MET_pt, MET_phi, met_passed,	n_jets,	n_btags, bWeight, ttHFCategory,	final_discriminant1, final_discriminant2, n_fatjets, pt_fatjet_1, pt_fatjet_2, pt_nonW_1, pt_nonW_2, pt_W1_1, pt_W1_2, pt_W2_1,	pt_W2_2, pt_top_1, pt_top_2, m_top_1, m_top_2, higgstag_fatjet_1, higgstag_fatjet_2, csv2_fatjet_1, csv2_fatjet_2 \n");
		
		events_e_cut1 = fopen("Outputs/CU_events_e_cut1.dat", "w");
		events_e_cut2 = fopen("Outputs/CU_events_e_cut2.dat", "w");
		events_e_cut3 = fopen("Outputs/CU_events_e_cut3.dat", "w");
		events_e_cut4 = fopen("Outputs/CU_events_e_cut4.dat", "w");
		events_e_cut5 = fopen("Outputs/CU_events_e_cut5.dat", "w");
		events_e_cut6 = fopen("Outputs/CU_events_e_cut6.dat", "w");
		events_e_cut7 = fopen("Outputs/CU_events_e_cut7.dat", "w");

		events_mu_cut1 = fopen("Outputs/CU_events_mu_cut1.dat", "w");
		events_mu_cut2 = fopen("Outputs/CU_events_mu_cut2.dat", "w");
		events_mu_cut3 = fopen("Outputs/CU_events_mu_cut3.dat", "w");
		events_mu_cut4 = fopen("Outputs/CU_events_mu_cut4.dat", "w");
		events_mu_cut5 = fopen("Outputs/CU_events_mu_cut5.dat", "w");
		events_mu_cut6 = fopen("Outputs/CU_events_mu_cut6.dat", "w");
		events_mu_cut7 = fopen("Outputs/CU_events_mu_cut7.dat", "w");
		
	}

	if (analysis_type == Analyze_dilepton) {
		events_di_lepton = fopen("CU_events_di_lepton_cuts.csv", "w");
		fprintf(events_di_lepton, "run, lumi, event, is_SL,	is_DL, lep1_pt, lep1_eta, lep1_phi, lep1_iso, lep1_pdgId, lep2_pt, lep2_eta, lep2_phi, lep2_iso, lep2_pdgId, mll,  mll_passed, jet1_pt,	jet2_pt, jet3_pt, jet4_pt, jet1_CSVv2, jet2_CSVv2, jet3_CSVv2, jet4_CSVv2, MET_pt, MET_phi, met_passed,	n_jets,	n_btags, bWeight, ttHFCategory,	final_discriminant1, final_discriminant2, n_fatjets, pt_fatjet_1, pt_fatjet_2, pt_nonW_1, pt_nonW_2, pt_W1_1, pt_W1_2, pt_W2_1,	pt_W2_2, pt_top_1, pt_top_2, m_top_1, m_top_2, higgstag_fatjet_1, higgstag_fatjet_2, csv2_fatjet_1, csv2_fatjet_2 \n");
		
		events_dimu_cut1 = fopen("Outputs/CU_events_dimu_cut1.dat", "w");
		events_dimu_cut2 = fopen("Outputs/CU_events_dimu_cut2.dat", "w");
		events_dimu_cut3 = fopen("Outputs/CU_events_dimu_cut3.dat", "w");
		events_dimu_cut4 = fopen("Outputs/CU_events_dimu_cut4.dat", "w");
		events_dimu_cut5 = fopen("Outputs/CU_events_dimu_cut5.dat", "w");
		events_dimu_cut6 = fopen("Outputs/CU_events_dimu_cut6.dat", "w");
		events_dimu_cut7 = fopen("Outputs/CU_events_dimu_cut7.dat", "w");

		events_diele_cut1 = fopen("Outputs/CU_events_diele_cut1.dat", "w");
		events_diele_cut2 = fopen("Outputs/CU_events_diele_cut2.dat", "w");
		events_diele_cut3 = fopen("Outputs/CU_events_diele_cut3.dat", "w");
		events_diele_cut4 = fopen("Outputs/CU_events_diele_cut4.dat", "w");
		events_diele_cut5 = fopen("Outputs/CU_events_diele_cut5.dat", "w");
		events_diele_cut6 = fopen("Outputs/CU_events_diele_cut6.dat", "w");
		events_diele_cut7 = fopen("Outputs/CU_events_diele_cut7.dat", "w");

		events_elemu_cut1 = fopen("Outputs/CU_events_elemu_cut1.dat", "w");
		events_elemu_cut2 = fopen("Outputs/CU_events_elemu_cut2.dat", "w");
		events_elemu_cut3 = fopen("Outputs/CU_events_elemu_cut3.dat", "w");
		events_elemu_cut4 = fopen("Outputs/CU_events_elemu_cut4.dat", "w");
		events_elemu_cut5 = fopen("Outputs/CU_events_elemu_cut5.dat", "w");
		
	}
	*/
}

void CU_ttH_EDA::Set_up_tokens(const edm::ParameterSet &config)
{
	if(!isdata)
		token.event_gen_info =
			consumes<GenEventInfoProduct>(edm::InputTag(std::string("generator")));
	token.triggerResults = consumes<edm::TriggerResults>(
		edm::InputTag(std::string("TriggerResults"), std::string(""), hltTag));
	token.filterResults = consumes<edm::TriggerResults>(edm::InputTag(
		std::string("TriggerResults"), std::string(""), filterTag));

	token.vertices = consumes<reco::VertexCollection>(
	    config.getParameter<edm::InputTag>("pv"));
	token.sec_vertices = consumes<reco::VertexCompositePtrCandidateCollection>(
	    config.getParameter<edm::InputTag>("sv"));
	token.PU_info = consumes<std::vector<PileupSummaryInfo>>(
	    config.getParameter<edm::InputTag>("pileup"));
	token.srcRho = consumes<double>(
	    config.getParameter<edm::InputTag>("rho"));
	token.electrons = consumes<pat::ElectronCollection>(
	    config.getParameter<edm::InputTag>("electrons"));
	token.muons = consumes<pat::MuonCollection>(
		config.getParameter<edm::InputTag>("muons"));
	token.jets = consumes<pat::JetCollection>(
	    config.getParameter<edm::InputTag>("jets"));
	token.METs = consumes<pat::METCollection>(
	    config.getParameter<edm::InputTag>("mets"));
	token.PF_candidates = consumes<pat::PackedCandidateCollection>(
	    config.getParameter<edm::InputTag>("pfcand"));
	token.BS = consumes<reco::BeamSpot>(
	    config.getParameter<edm::InputTag>("beamspot"));
	//token.MC_particles = consumes<reco::GenParticleCollection>(
	//    config.getParameter<edm::InputTag>("prunedgen"));
	//token.MC_packed = consumes<pat::PackedGenParticleCollection>(
	//    config.getParameter<edm::InputTag>("packedgen"));
	token.mvaValuesMapToken_ = consumes<edm::ValueMap<float>>(
	    config.getParameter<edm::InputTag>("mvaValues"));
   	token.mvaCategoriesMapToken_ = consumes<edm::ValueMap<int>>(
   	    config.getParameter<edm::InputTag>("mvaCategories"));
   	token.electrons_for_mva_token = consumes<edm::View<pat::Electron>>(
   	    config.getParameter<edm::InputTag>("electrons"));
   	token.genTtbarIdToken_ = consumes<int>(
   	    config.getParameter<edm::InputTag>("genTtbarId"));
	token.puInfoToken = consumes<std::vector< PileupSummaryInfo > >(
   	    config.getParameter<edm::InputTag>("pileupinfo"));
   	token.lheptoken = consumes<LHEEventProduct>(
   	    config.getParameter<edm::InputTag>("lhepprod"));
   	
	
}

void CU_ttH_EDA::Set_up_Tree()
{
	eventTree = fs_->make<TTree>("eventTree", "Event tree");
	
	/*
	// If ntuple class is inherited from Root TClass (ToDo)
	//
	//eventTree -> Branch("ntuple_", "CU_ttH_EDA_Ntuple", &ntuple);
	//std::cout << "IsTObject :" << ntuple->IsTObject() <<  std::endl;
	//std::cout << "GetNdata() :" << ntuple->GetNdata() << std::endl;
	//std::cout << "CanSplit() :" << ntuple->CanSplit() << std::endl;
	//ntuple->Dump();
	*/
	//tauNtuple.set_up_branches(eventTree);
}

void CU_ttH_EDA::Set_up_b_weights(){
	//inputFileHF = "MiniAOD/MiniAODHelper/data/csv_rwt_fit_hf_2015_11_20.root";
  	//inputFileLF = "MiniAOD/MiniAODHelper/data/csv_rwt_fit_lf_2015_11_20.root";
  	inputFileHF = "Analyzers/data/csv_rwt_fit_hf_v2_final_2016_06_30test.root";
  	inputFileLF = "Analyzers/data/csv_rwt_fit_lf_v2_final_2016_06_30test.root";
	f_CSVwgt_HF = new TFile ((std::string(getenv("CMSSW_BASE")) + "/src/" + inputFileHF).c_str());
	f_CSVwgt_LF = new TFile ((std::string(getenv("CMSSW_BASE")) + "/src/" + inputFileLF).c_str());
	fillCSVHistos(f_CSVwgt_HF, f_CSVwgt_LF);
}

#endif
