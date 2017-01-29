#ifndef CU_ttH_EDA_Setups_cc
#define CU_ttH_EDA_Setups_cc

/// Includes
#include "CU_ttH_EDA.h"

void CU_ttH_EDA::Set_up_weights()
{
    Set_up_b_weights();
    init_PU_weight();
    init_PDF_weight();
}

void CU_ttH_EDA::init_weights(CU_ttH_EDA_event_vars &local)
{
    local.gen_weight = -1;
    local.PU_weight = -1;
    local.pdf_weight_up = -1;
    local.pdf_weight_down = -1;
    local.q2_weight_up = -1;
    local.q2_weight_down = -1;
}

void CU_ttH_EDA::init_bjetness_var(CU_ttH_EDA_event_vars &local)
{
  local.bjetnessFV_num_leps = -1; 
  local.bjetnessFV_npvTrkOVcollTrk =-1;
  local.bjetnessFV_avip3d_val =-1;
  local.bjetnessFV_avip3d_sig =-1; 
  local.bjetnessFV_avsip3d_sig =-1; 
  local.bjetnessFV_avip1d_sig =-1;
} 

void CU_ttH_EDA::init_flags(CU_ttH_EDA_event_vars &local)
{
    local.is_e = false;
    local.is_mu = false;
    local.is_ee = false;
    local.is_emu = false;
    local.is_mumu = false;

    // flag for determining whether to select an event for writing
    local.event_selection_SL = false;
    local.event_selection_DL = false;
}

void CU_ttH_EDA::init_PU_weight()
{
    ifstream fin;
    fin.open("data/PU_weight/PU_weights.txt");
    // fin.open("/afs/cern.ch/user/a/abdatta/public/PU_weight/PU_weights.txt");
    for (int i = 0; i < 75; ++i) {
        fin >> PU_x[i] >> PU_y[i];
    }
    fin.close();
}

void CU_ttH_EDA::init_PDF_weight()
{
    NNPDF30_nlo_as_0118_PDFSet = new LHAPDF::PDFSet("NNPDF30_nlo_as_0118");
    _systPDFs = NNPDF30_nlo_as_0118_PDFSet->mkPDFs();
}

void CU_ttH_EDA::Close_output_files() { fclose(events_combined); }

/*
void CU_ttH_EDA::Set_up_histograms()
{
    if (analysis_type == Analyze_lepton_jet) {
    }
    if (analysis_type == Analyze_dilepton) {
    }
}
*/

/*
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
*/

/*
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
*/

/// Make and open write-out files
void CU_ttH_EDA::Set_up_output_files()
{
    events_combined = fopen("Cornell_tth_80X.csv", "w");
    fprintf(events_combined,
            "run,lumi,event,is_e,is_mu,is_ee,is_emu,is_mumu,n_jets,n_btags,"
            "lep1_pt,lep1_iso,lep1_pdgId,lep2_pt,lep2_iso,lep2_pdgId,jet1_pt,"
            "jet1_eta,jet1_phi,jet1_jesSF,jet1_jesSF_up,jet1_jesSF_down,jet1_csv,"
            "jet2_pt,jet2_eta,jet2_phi,jet2_jesSF,jet2_jesSF_up,jet2_jesSF_down,jet2_csv,"
            "MET_pt,MET_phi,mll,ttHFCategory,n_interactions,PUWeight,csvSF,"
            "csvSF_lf_up,csvSF_hf_down,csvSF_cErr1_down," 				
            "pdf_up,pdf_down,me_up,me_down\n");
}

//bjetnessFV_num_leps,bjetnessFV_npvTrkOVcollTrk,bjetnessFV_avip3d_val,bjetnessFV_avip3d_sig,bjetnessFV_avsip3d_sig,bjetnessFV_avip1d_sig

void CU_ttH_EDA::Set_up_tokens(const edm::ParameterSet &config)
{
    if (!isdata)
        token.event_gen_info = consumes<GenEventInfoProduct>(
            edm::InputTag(std::string("generator")));
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
    token.srcRho = consumes<double>(config.getParameter<edm::InputTag>("rho"));
    token.electrons = consumes<pat::ElectronCollection>(
        config.getParameter<edm::InputTag>("electrons"));
    token.muons = consumes<pat::MuonCollection>(
        config.getParameter<edm::InputTag>("muons"));
    token.jets = consumes<pat::JetCollection>(
        config.getParameter<edm::InputTag>("jets"));
    token.METs = consumes<pat::METCollection>(
        config.getParameter<edm::InputTag>("mets"));
    token.genjets = consumes<reco::GenJetCollection>(
        config.getParameter<edm::InputTag>("genjets"));
    token.BadChCandFilterToken_ = consumes<bool>(
    	config.getParameter<edm::InputTag>("badchcandfilter"));
    token.BadPFMuonFilterToken_ = consumes<bool>(
    	config.getParameter<edm::InputTag>("badpfmufilter"));
    //token.BadChCandFilterToken_ = consumes<bool>(
    //    edm::InputTag(std::string("BadChargedCandidateFilter")));
    //token.BadPFMuonFilterToken_ = consumes<bool>(
    //    edm::InputTag(std::string("BadPFMuonFilter")));
    token.PF_candidates = consumes<pat::PackedCandidateCollection>(
        config.getParameter<edm::InputTag>("pfcand"));
    token.BS = consumes<reco::BeamSpot>(
        config.getParameter<edm::InputTag>("beamspot"));
    token.eleMediumIdMapToken_ = consumes<edm::ValueMap<bool> >(
        config.getParameter<edm::InputTag>("eleMediumIdMap")),
    token.mvaValuesMapToken_ = consumes<edm::ValueMap<float>>(
        config.getParameter<edm::InputTag>("mvaValues"));
    token.mvaCategoriesMapToken_ = consumes<edm::ValueMap<int>>(
        config.getParameter<edm::InputTag>("mvaCategories"));
    //token.electrons_for_mva_token = mayConsume<edm::View<reco::GsfElectron>>(
    //    config.getParameter<edm::InputTag>("electrons"));
    token.electrons_for_mva_token = consumes<edm::View<pat::Electron>>(
        config.getParameter<edm::InputTag>("electrons"));
    token.muon_h_token = consumes<edm::View<pat::Muon>>(
        config.getParameter<edm::InputTag>("muons"));
    token.genTtbarIdToken_ =
        consumes<int>(config.getParameter<edm::InputTag>("genTtbarId"));
    token.puInfoToken = consumes<std::vector<PileupSummaryInfo>>(
        config.getParameter<edm::InputTag>("pileupinfo"));
    token.lheptoken = consumes<LHEEventProduct>(
        config.getParameter<edm::InputTag>("lhepprod"));
}

void CU_ttH_EDA::Set_up_Tree()
{
    eventTree = fs_->make<TTree>("eventTree", "Event tree");
    hbbNtuple.set_up_branches(eventTree);
    //hbbNtuple = 0;
    //eventTree->Branch("hbbNtuple", "CU_ttH_EDA_Ntuple", &hbbNtuple, 8000, 1);
}

void CU_ttH_EDA::Set_up_b_weights()
{
    inputFileHF =
    //    "data/csv_weights/csv_rwt_fit_hf_v2_final_2016_06_30test.root";
    //	  "data/csv_weights/csv_rwt_fit_hf_v2_final_2016_09_7test.root";
          "data/csv_weights/csv_rwt_fit_hf_v2_final_2017_1_10test.root";
    inputFileLF =
    //    "data/csv_weights/csv_rwt_fit_lf_v2_final_2016_06_30test.root";
    //	  "data/csv_weights/csv_rwt_fit_lf_v2_final_2016_09_7test.root";
          "data/csv_weights/csv_rwt_fit_lf_v2_final_2017_1_10test.root";
    f_CSVwgt_HF = new TFile((inputFileHF).c_str());
    f_CSVwgt_LF = new TFile((inputFileLF).c_str());
    fillCSVHistos(f_CSVwgt_HF, f_CSVwgt_LF);
}

#endif
