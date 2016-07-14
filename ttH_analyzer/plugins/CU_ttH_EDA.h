#ifndef CU_ttH_EDA_h
#define CU_ttH_EDA_h

/// Core libraries
#include <memory>
#include <cstdio>	// printf, fprintf
#include <stdexcept> // standard exceptions
#include <cmath>      // arctan
#include <fstream>

/// CMSSW user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

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
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMapFwd.h"
#include "DataFormats/L1Trigger/interface/L1HFRings.h"
#include "DataFormats/L1Trigger/interface/L1HFRingsFwd.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/strbitset.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

#include "MiniAOD/MiniAODHelper/interface/MiniAODHelper.h"
#include "MiniAOD/MiniAODHelper/interface/LeptonSFHelper.h"

#include "LHAPDF/LHAPDF.h"

/// ROOT includes
#include "TH1.h"
#include "TTree.h"
#include "TRandom3.h"

/// Higgs and top tagger
//#include "MiniAOD/BoostedObjects/interface/HTTTopJet.h"
//#include "MiniAOD/BoostedObjects/interface/SubFilterJet.h"
//#include "BoostedTTH/BoostedAnalyzer/interface/BoostedUtils.hpp"
#include "MiniAOD/MiniAODHelper/interface/TopTagger.h"
#include "MiniAOD/MiniAODHelper/interface/HiggsTagger.h"

/// structs for holding multiple edm::Handle and EDGetTokenT
#include "CU_ttH_EDA_Handles.h"

#include "Analyzers/ttH_analyzer/interface/CU_ttH_EDA_Ntuple.h"
#include "Analyzers/ttH_analyzer/interface/LHEEventProduct.h"

/// Configuration reader
#include "yaml-cpp/yaml.h"

/// TMVA
#include "TMVA/Reader.h"

/*
 *
 * enum for analysis type.
 * Purpose: allows faster in-function comparisons
 *
 */

enum analysis_types {
	Analyze_lepton_jet,
	Analyze_dilepton
};

 /*
 *
 * Analyzer class
 *
 */

class CU_ttH_EDA : public edm::EDAnalyzer
{
  public:
	explicit CU_ttH_EDA(const edm::ParameterSet &);
	~CU_ttH_EDA();

	static void fillDescriptions(edm::ConfigurationDescriptions &);

  private:
	/*
	 * Function/method section
	*/

	/// Standard EDAnalyzer functions
	void analyze(const edm::Event &, const edm::EventSetup &) override;

	void beginJob() override;
	void endJob() override;

	void beginRun(const edm::Run &, const edm::EventSetup &) override;
	void endRun(const edm::Run &, const edm::EventSetup &) override;

	/// One-time-run functions
	void Close_output_files();		 // at ~CU_ttH_EDA()
	void Load_configuration(string); // at CU_ttH_EDA(), runs _MAODH()
	void Load_configuration_set_type(const string &); // sets analysis_type
	void Load_configuration_MAODH(bool); // runs miniAODhelper.SetUp
	void Set_up_histograms();			 // at CU_ttH_EDA()
	void Set_up_output_files();			 // at CU_ttH_EDA()
	void Set_up_tokens(const edm::ParameterSet &);
	void Set_up_b_weights();			// at CU_ttH_EDA()
	void Set_up_Tree();

	int Set_up_Run_histograms_triggers(); // at beginRun(), after
										  // Set_up_name_vectors()

	void Set_up_trigger_name_vectors(); // at beginRun()

	int End_Run_hist_fill_triggers(); // fill histograms at endRun()

	/// Per-event functions
	void Update_common_vars(const edm::Event &, CU_ttH_EDA_event_vars &);

	/// Object checks. Returns 1 in case of an error
	int Check_beam_spot(edm::Handle<reco::BeamSpot>);
	int Check_triggers(edm::Handle<edm::TriggerResults>,
					   CU_ttH_EDA_event_vars &); // adjusts event variables
	int Check_filters(edm::Handle<edm::TriggerResults>);
	int Check_vertices_set_MAODhelper(edm::Handle<reco::VertexCollection>);
	int Check_PV(edm::Handle<reco::VertexCollection>);   // check primary vertex

	// trigger iterator, part of Check_triggers()
	bool Check_triggers_iterator(const vector<string> &, edm::Handle<edm::TriggerResults>);
	
	//std::vector<pat::Jet> CheckJetID (const std::vector<pat::Jet>&);
	std::vector<pat::Jet> CheckJetID (const std::vector<pat::Jet>&, const std::vector<pat::Jet>&);
	std::vector<pat::Jet> GetCorrectedJets_JEC(const std::vector<pat::Jet>&, double, const sysType::sysType iSysType=sysType::NA, const float& corrFactor = 1, const float& uncFactor = 1);
	std::vector<pat::Jet> GetCorrectedJets_JER(const std::vector<pat::Jet>&, double, JME::JetResolution, const sysType::sysType iSysType=sysType::NA, const float& corrFactor = 1, const float& uncFactor = 1);
	void SetFactorizedJetCorrector(const sysType::sysType iSysType=sysType::NA);
	double getJERfactor( const int, const double, const double, const double );
	double GetJetSF( pat::Jet, const sysType::sysType, double);
	void getbweight (CU_ttH_EDA_event_vars &);
	double getPUweight ( edm::Handle<std::vector< PileupSummaryInfo > >  );
	
	
	double getCSVWeight(std::vector<double> jetPts, std::vector<double> jetEtas, std::vector<double> jetCSVs,
                       std::vector<int> jetFlavors, int iSys, double &csvWgtHF, double &csvWgtLF, double &csvWgtCF);
        void fillCSVHistos(TFile *fileHF, TFile *fileLF);
        
        //void Fill_addn_quant(CU_ttH_EDA_event_vars &, double, edm::Handle<int>, edm::Handle<GenEventInfoProduct> );
        void Fill_addn_quant(CU_ttH_EDA_event_vars &, double, edm_Handles );
	void getPDFweight(CU_ttH_EDA_event_vars &, edm::Handle<GenEventInfoProduct> );
	double getQ2weight( edm::Handle<GenEventInfoProduct>, edm::Handle<LHEEventProduct>, string);
	
	void Check_SL_Event_Selection(CU_ttH_EDA_event_vars &);
	void Check_DL_Event_Selection(CU_ttH_EDA_event_vars &);
	void init_PU_weight();
	void init_PDF_weight();
	
	ifstream fin;
	
	/// Taggers. Returns 1 in case of an error
	//int Higgs_tagger(Handle<boosted::SubFilterJetCollection>,
	//				 CU_ttH_EDA_event_vars &); // FIXME: uses b-tag medium WP
	//int Top_tagger(Handle<boosted::HTTTopJetCollection>,
	//				   CU_ttH_EDA_event_vars &);
	//TopTagger toptagger;

	/// Other functions
	void Check_Fill_Print_single_lepton(CU_ttH_EDA_event_vars &);
	void Check_Fill_Print_di_lepton(CU_ttH_EDA_event_vars &);
	void Check_Fill_Print_ej(CU_ttH_EDA_event_vars &);
	void Check_Fill_Print_muj(CU_ttH_EDA_event_vars &);
	void Check_Fill_Print_dimuj(CU_ttH_EDA_event_vars &);
	void Check_Fill_Print_dielej(CU_ttH_EDA_event_vars &);
	void Check_Fill_Print_elemuj(CU_ttH_EDA_event_vars &);

	void printDecayChain(const reco::Candidate &p, int &index, int mother_index,
						 bool details);
	
	template <class lepton>
	int Print_event_in_file1(FILE *, lepton &, std::vector<pat::Jet> &,
							 CU_ttH_EDA_event_vars &);
	template <class lep1, class lep2>
	// template<class lepton2>
	int Print_event_in_file1_dilepton(FILE *, lep1 &, lep2 &, double,
									  std::vector<pat::Jet> &,
									  CU_ttH_EDA_event_vars &);

	template <typename T1, typename T2>
		std::vector<T1>
		removeOverlapdR(const std::vector<T1>& v1, const std::vector<T2>& v2, double dR = 0.02);

	//float getMHT(CU_ttH_EDA_event_vars &);
	
	/*
	* Variable section
	*/

	// Analysis type
	analysis_types analysis_type;
	std::string config_analysis_type;
	
	/// debug flags
	bool verbose_;
	bool dumpHLT_;
	
	edm_Tokens token; // common tokens for all events

	/// Triggers, paths: configs filled/updated via run
	HLTConfigProvider hlt_config;
	HLTConfigProvider filter_config;

	/// Triggers, paths.
	// Used for trigger statistics, filled via run (trigger_stats = true)
	std::string hltTag;
	std::string filterTag;

	bool trigger_stats;
	
	// counters (trigger_stats = true)
	std::map<std::string, unsigned long> n_trigger_fired; // HLT
	std::map<std::string, unsigned long> n_filter_fired;

	std::vector<std::string> trigger_names;
	std::vector<std::string> trigger_names_no_ver;

	std::vector<std::string> filter_names;
	std::vector<std::string> filter_names_no_ver;

	// triggers of interest. Provided by config
	std::vector<std::string> trigger_on_HLT_e;	// single electron trigger
	std::vector<std::string> trigger_on_HLT_mu;   // single muon trigger
	std::vector<std::string> trigger_on_HLT_ee;   // dielectron tigger
	std::vector<std::string> trigger_on_HLT_emu;  // electron+muon trigger
	std::vector<std::string> trigger_on_HLT_mumu; // dimuon trigger

	/// Output file is opened/closed through CMS py config
	edm::Service<TFileService> fs_;

	/// Common sample parameters
	unsigned long event_count; // running event counter
	unsigned long selection_count; // counting event selections
	int sl_e;
	int sl_mu;
	int dl_ee;
	int dl_emu;
	int dl_mumu;

	double total_xs;	  // total cross section
	double int_lumi;	  // integrated luminosity
	double sample_n;	  // total nr of events. Should be long if compatible
	double weight_sample; // int lumi * xs / sample_n
	// double weight_gen;

	/// Cuts
	float min_ele_pT;	
	float min_mu_pT;
	float min_veto_ele_pT;
	float min_veto_mu_pT;
	float min_di_ele1_pT;
	float min_di_ele2_pT;
	float min_di_mu1_pT;
	float min_di_mu2_pT;
	float min_jet_pT;
	float min_jet2_pT;
	float min_bjet_pT;
	float max_ele_eta;
	float max_mu_eta;
	float max_veto_ele_eta;
	float max_veto_mu_eta;
	float max_di_ele1_eta;
	float max_di_ele2_eta;
	float max_di_mu1_eta;
	float max_di_mu2_eta;
	float max_jet_eta;
	float max_bjet_eta;
	int min_njets;
	int min_di_njets;
	int min_nbtags;
	int min_di_nbtags;
	float min_di_mll;
	float min_di_met;
	
	// for JEC
	FactorizedJetCorrector* _jetCorrector;
	JetCorrectionUncertainty* _jetCorrectorUnc;
	
	// for JER
	//JME::JetResolution resolution;
	
	int i;
	std::vector<int> index;
	
	double E;
	double p;
	double pz;
	double py;
	double px;
	
	// for b-weights
	std::string inputFileHF;
  	std::string inputFileLF;
  	TFile *f_CSVwgt_HF;
  	TFile *f_CSVwgt_LF;
  	TH1D *h_csv_wgt_hf[9][6];
    	TH1D *c_csv_wgt_hf[9][6];
    	TH1D *h_csv_wgt_lf[9][4][3];
	
	/// Selection helper
	MiniAODHelper miniAODhelper;
	LeptonSFHelper leptonSFhelper;
	
	// for PDF weight
	LHAPDF::PDFSet *NNPDF30_nlo_as_0118_PDFSet;
	std::vector<LHAPDF::PDF*> _systPDFs;

	bool isdata;
	char MAODHelper_b_tag_strength;
	int MAODHelper_sample_nr; // past insample_, in-development var. for
							  // MAODHelper?
	std::string MAODHelper_era;

	int PU_x[100];
	double PU_y[100];
	
	 TRandom3 *r;

	/// Histograms
	TH1D *h_tth_syncex1_ele;
	TH1D *h_tth_syncex1_mu;
	TH1D *h_tth_syncex1_dimu;
	TH1D *h_tth_syncex1_diele;
	TH1D *h_tth_syncex1_elemu;

	TH1D *h_tth_syncex_dimutauh;
	TH1D *h_tth_syncex_dieletauh;
	TH1D *h_tth_syncex_elemutauh;
	
	TH1D *h_tth_syncex_eleditauh;
	TH1D *h_tth_syncex_muditauh;

	TH1D *h_hlt;
	TH1D *h_flt;

	/// Write-out files
	FILE *events_single_lepton;
	FILE *events_di_lepton;
	FILE *events_combined;
	
	FILE *events_e_cut1, *events_e_cut2, *events_e_cut3, *events_e_cut4,
		*events_e_cut5, *events_e_cut6, *events_e_cut7;

	FILE *events_mu_cut1, *events_mu_cut2, *events_mu_cut3, *events_mu_cut4,
		*events_mu_cut5, *events_mu_cut6, *events_mu_cut7;

	FILE *events_dimu_cut1, *events_dimu_cut2, *events_dimu_cut3,
		*events_dimu_cut4, *events_dimu_cut5, *events_dimu_cut6,
		*events_dimu_cut7;

	FILE *events_diele_cut1, *events_diele_cut2, *events_diele_cut3,
		*events_diele_cut4, *events_diele_cut5, *events_diele_cut6,
		*events_diele_cut7;

	FILE *events_elemu_cut1, *events_elemu_cut2, *events_elemu_cut3,
		*events_elemu_cut4, *events_elemu_cut5;

	// tree and ntuple
	TTree *eventTree;
	CU_ttH_EDA_Ntuple bbNtuple;

};

template <typename T1, typename T2>
std::vector<T1>
CU_ttH_EDA::removeOverlapdR(const std::vector<T1> &v1, const std::vector<T2> &v2, double dR)
{
	std::vector<T1> res;
	for (const auto& o1: v1) {
		bool keep = true;
		for (const auto& o2: v2)
			if (miniAODhelper.DeltaR(&o1, &o2) < dR)
				keep = false;
		if (keep)
			res.push_back(o1);
	}
	return res;
}

#endif
