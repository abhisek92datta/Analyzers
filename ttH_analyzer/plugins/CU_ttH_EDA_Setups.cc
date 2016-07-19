#ifndef CU_ttH_EDA_Setups_cc
#define CU_ttH_EDA_Setups_cc

/// Includes
#include "CU_ttH_EDA.h"

void CU_ttH_EDA::init_PU_weight()
{
	ifstream fin;
	//fin.open("ttH_analyzer/data/PU_weights.txt");
	fin.open(string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/PU_weights.txt");
	for(int i=0; i<50; i++) {
		fin>>PU_x[i]>>PU_y[i];
	}
	fin.close();
}

void CU_ttH_EDA::init_PDF_weight()
{
	NNPDF30_nlo_as_0118_PDFSet = new LHAPDF::PDFSet("NNPDF30_nlo_as_0118");
	_systPDFs = NNPDF30_nlo_as_0118_PDFSet->mkPDFs();
}

void CU_ttH_EDA::Close_output_files()
{
	fclose(events_combined);
	}

void CU_ttH_EDA::Set_up_histograms()
{
	if (analysis_type == Analyze_lepton_jet) {
	}
	if (analysis_type == Analyze_dilepton) {
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
	fprintf(events_combined, "run,lumi,event,is_e,is_mu,is_ee,is_emu,is_mumu,n_jets,n_btags,lep1_pt,lep1_iso,lep1_pdgId,lep2_pt,lep2_iso,lep2_pdgId,jet1_pt,jet2_pt,jet1_CSVv2,jet2_CSVv2,jet1_JecSF,jet1_JecSF_up,jet1_JecSF_down,MET_pt,MET_phi,mll,ttHFCategory,PUWeight,bWeight,triggerSF,lepIDSF,lepISOSF,Q2_upup,Q2_downdown,pdf_up,pdf_down\n");
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
}

void CU_ttH_EDA::Set_up_b_weights(){
	inputFileHF = "Analyzers/ttH_analyzer/data/csv_rwt_fit_hf_v2_final_2016_06_30test.root";
  	inputFileLF = "Analyzers/ttH_analyzer/data/csv_rwt_fit_lf_v2_final_2016_06_30test.root";
	f_CSVwgt_HF = new TFile ((std::string(getenv("CMSSW_BASE")) + "/src/" + inputFileHF).c_str());
	f_CSVwgt_LF = new TFile ((std::string(getenv("CMSSW_BASE")) + "/src/" + inputFileLF).c_str());
	fillCSVHistos(f_CSVwgt_HF, f_CSVwgt_LF);
}

#endif
