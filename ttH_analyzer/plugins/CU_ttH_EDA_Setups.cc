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

void CU_ttH_EDA::init_SFs(CU_ttH_EDA_event_vars &local)
{
    local.b_weight_sl = -1;
    local.b_weight_sl_lfup = -1;
    local.b_weight_sl_hfdown = -1;
    local.b_weight_sl_cErr1_down = -1;
    local.b_weight_di = -1;
    local.b_weight_di_lfup = -1;
    local.b_weight_di_hfdown = -1;
    local.b_weight_di_cErr1_down = -1;

    local.lep_sf_id_sl.clear();
    local.lep_sf_iso_sl.clear();
    local.lep_sf_gsf_sl = -1;
    local.lep_sf_hip_sl = -1;
    local.lep_sf_trig_sl = -1;
    local.lep_sf_id_di.clear();
    local.lep_sf_iso_di.clear();
    local.lep_sf_gsf_di = -1;
    local.lep_sf_hip_di = -1;
    local.lep_sf_trig_di = -1;
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
    for (int i = 0; i < 75; ++i) {
        fin >> PU_x[i] >> PU_y[i];
    }
    fin.close();
}

void CU_ttH_EDA::init_PDF_weight()
{
    pdfSet = new LHAPDF::PDFSet("PDF4LHC15_nlo_30");
    _systPDFs = pdfSet->mkPDFs();
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
            "lep1_pt,lep1_eta,lep1_iso,lep1_pdgId,lep1_idSF,lep1_isoSF,lep1_seed,"
            "lep2_pt,lep2_eta,lep2_iso,lep2_pdgId,lep2_idSF,lep2_isoSF,lep2_seed,"
            "jet1_pt,jet1_eta,jet1_phi,jet1_jesSF,jet1_jesSF_up,jet1_jesSF_down,jet1_jesSF_PileUpDataMC_down,jet1_jesSF_RelativeFSR_up,"
            "jet1_jerSF_nominal,jet1_csv,jet1_PUJetId,jet1_PUJetDiscriminant,jet1_seed,"
            "jet2_pt,jet2_eta,jet2_phi,jet2_jesSF,jet2_jesSF_up,jet2_jesSF_down,jet2_jesSF_PileUpDataMC_down,jet2_jesSF_RelativeFSR_up,"
            "jet2_jerSF_nominal,jet2_csv,jet2_PUJetId,jet2_PUJetDiscriminant,jet2_seed,"
            "MET_pt,MET_phi,MET_pt_phiCor,MET_phi_phiCor,mll,ttHFCategory,ttHFGenFilterTag,n_interactions,puWeight,csvSF,"
            "csvSF_lf_up,csvSF_hf_down,csvSF_cErr1_down," 				
            "pdf_up,pdf_down,me_up,me_down,triggerSF,top_pt_weight,bdt_output,dnn_ttH_output,dnn_ttbb_output\n");
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
    token.BadGlobalMuonTaggerToken_ = consumes<bool>(
        config.getParameter<edm::InputTag>("badglobalmutagger"));
    token.CloneGlobalMuonTaggerToken_ = consumes<bool>(
        config.getParameter<edm::InputTag>("cloneglobalmutagger"));
    token.PF_candidates = consumes<pat::PackedCandidateCollection>(
        config.getParameter<edm::InputTag>("pfcand"));
    token.BS = consumes<reco::BeamSpot>(
        config.getParameter<edm::InputTag>("beamspot"));
    token.eleTightIdMapToken_ = consumes<edm::ValueMap<bool> >(
        config.getParameter<edm::InputTag>("eleTightIdMap")),
    //token.mvaValuesMapToken_ = consumes<edm::ValueMap<float>>(
    //    config.getParameter<edm::InputTag>("mvaValues"));
    //token.mvaCategoriesMapToken_ = consumes<edm::ValueMap<int>>(
    //    config.getParameter<edm::InputTag>("mvaCategories"));
    //token.electrons_for_mva_token = mayConsume<edm::View<reco::GsfElectron>>(
    //    config.getParameter<edm::InputTag>("electrons"));
    token.electrons_for_mva_token = consumes<edm::View<pat::Electron>>(
        config.getParameter<edm::InputTag>("electrons"));
    token.muon_h_token = consumes<edm::View<pat::Muon>>(
        config.getParameter<edm::InputTag>("muons"));
    if (!isdata) {
        token.genTtbarIdToken_ =
            consumes<int>(config.getParameter<edm::InputTag>("genTtbarId"));
        token.ttHFGenFilterToken_ =
            consumes<bool>(config.getParameter<edm::InputTag>("ttHFGenFilter"));
    }
    token.puInfoToken = consumes<std::vector<PileupSummaryInfo>>(
        config.getParameter<edm::InputTag>("pileupinfo"));
    if (!isdata && !is_OLS) {
    token.lheptoken = consumes<LHEEventProduct>(
        config.getParameter<edm::InputTag>("lhepprod"));
    }
}

void CU_ttH_EDA::Set_up_Tree()
{
    //eventTree = fs_->make<TTree>();
    eventTree = fs_->make<TTree>("eventTree", "Event tree");
    hbbNtuple.set_up_branches(eventTree);
    
    //hbbNtuple = 0;
    //eventTree->Branch("hbbNtuple", "CU_ttH_EDA_Ntuple", &hbbNtuple, 8000, 1);
}

void CU_ttH_EDA::Set_up_b_weights()
{
    inputFileHF =
    //    "data/csv_weights/csv_rwt_fit_hf_v2_final_2017_1_10test.root";
          "data/csv_weights/csv_rwt_fit_hf_v2_final_2017_3_29test.root";
    inputFileLF =
    //    "data/csv_weights/csv_rwt_fit_lf_v2_final_2017_1_10test.root";
          "data/csv_weights/csv_rwt_fit_lf_v2_final_2017_3_29test.root";

    f_CSVwgt_HF = new TFile((inputFileHF).c_str());
    f_CSVwgt_LF = new TFile((inputFileLF).c_str());
    fillCSVHistos(f_CSVwgt_HF, f_CSVwgt_LF);
}

#endif
