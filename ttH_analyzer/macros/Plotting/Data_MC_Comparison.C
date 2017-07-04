#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH3.h"
#include "TH2F.h"
#include "TF1.h"
#include "TF2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TAxis.h"
#include "TMath.h"
#include "TRandom3.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <exception>
#include <cmath> 
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include "TGraphAsymmErrors.h"
#include "TVector.h"
#include "TLorentzVector.h"

void Data_MC_Comparison( int maxNentries=-1, int Njobs=1, int jobN=1 ) {

    int N_total = 0;
    int N_4b_pass = 0;
    int N_ttbb_pass = 0;
    int N_sel = 0;
    int N_SL = 0;
    int N_DL = 0;
	int N_e = 0;
	int N_mu = 0;
    int N_ee = 0;
    int N_emu = 0;
    int N_mumu = 0;

    double sum_gen_weights = 0;

    int isttbar = 0;  // Set 1 for ttbar datasets

    // input filename
    std::string treefilename = "ttHbbNtuple.root";
    //std::string treefilename = "/eos/uscms/store/user/adatta/ttH_Analysis/v1/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/ttH_Analyzer_ttjets_pp/170603_213946/0000/ttHbbNtuple*.root";
    //std::string treefilename = "/eos/uscms/store/user/adatta/ttH_Analysis/v1/ttbb_4FS_OpenLoops_13TeV-sherpa/ttH_Analyzer_ttjets_ol/170606_013446/0000/ttHbbNtuple*.root";

    std::string histofilename = "Distribution.root";
	//std::string histofilename = "Distribution_pp.root";
    //std::string histofilename = "Distribution_ol.root";
	
	TChain *chain = new TChain("ttHbb/eventTree");

    chain->Add(treefilename.c_str());

    TFile histofile(histofilename.c_str(),"recreate");
    histofile.cd();
    
	TH1::SetDefaultSumw2();

    // Single Electron Histograms
    TH1D* pt_e1_sl_e = new TH1D("pt_e1_sl_e","Leading Electron pT Distribution (Single Electron) ;Leading Electron pT (GeV);Nr. of Events",50,0,300);
    TH1D* eta_e1_sl_e = new TH1D("eta_e1_sl_e","Leading Electron eta Distribution (Single Electron) ;Leading Electron #eta; Nr. of Events",15,-3,3);
    TH1D* phi_e1_sl_e = new TH1D("phi_e1_sl_e","Leading Electron phi Distribution (Single Electron) ;Leading Electron #phi; Nr. of Events",10,-3,3);
    TH1D* pt_jet1_sl_e = new TH1D("pt_jet1_sl_e","Leading Jet pT Distribution (Single Electron) ;Jet 1 pT (GeV);Nr. of Events",50,0,300);
    TH1D* pt_jet2_sl_e = new TH1D("pt_jet2_sl_e","2nd Leading Jet pT Distribution (Single Electron) ;Jet 2 pT (GeV);Nr. of Events",50,0,300);
    TH1D* pt_jet3_sl_e = new TH1D("pt_jet3_sl_e","3rd Leading Jet pT Distribution (Single Electron) ;Jet 3 pT (GeV);Nr. of Events",50,0,300);
    TH1D* pt_jet4_sl_e = new TH1D("pt_jet4_sl_e","4th Leading Jet pT Distribution (Single Electron) ;Jet 4 pT (GeV);Nr. of Events",50,0,300);
    TH1D* pt_jet5_sl_e = new TH1D("pt_jet5_sl_e","5th Leading Jet pT Distribution (Single Electron) ;Jet 5 pT (GeV);Nr. of Events",50,0,300);
    TH1D* pt_jet6_sl_e = new TH1D("pt_jet6_sl_e","6th Leading Jet pT Distribution (Single Electron) ;Jet 6 pT (GeV);Nr. of Events",50,0,300);
    TH1D* pt_jet_all_sl_e = new TH1D("pt_jet_all_sl_e","Jet (all) pT Distribution (Single Electron) ;Jet pT (GeV);Nr. of Events",50,0,300);
    TH1D* csv_jet_all_sl_e = new TH1D("csv_jet_all_sl_e","Jet (all) CSV Distribution (Single Electron) ;Jet CSV (GeV);Nr. of Events",50,0,1);
    TH1D* ht_sl_e = new TH1D("ht_sl_e","HT Distribution (Single Electron) ;HT (GeV);Nr. of Events",50,0,1000);
    TH2D* njets_ht_sl_e = new TH2D("njets_ht_2D_sl_e","Jet Multiplicity and HT 2D Distribution (Single Electron) ;Nr. of jets ;HT (GeV)",15,0,15,50,0,1000);
    TH1D* njets_sl_e = new TH1D("njets_sl_e","Jet Multiplicity (Single Electron) ;Nr. of jets; Nr. of Events",15,0,15);
    TH1D* nbtags_sl_e = new TH1D("nbtags_sl_e","B-jet Multiplicity (Single Electron) ;Nr. of b-jets; Nr. of Events",8,0,8);
    TH1D* npv_sl_e = new TH1D("npv_sl_e","Nr. of Primary Vertices (Single Electron)  ;NPV; Nr. of Events",20,0,75);

    // Single Muon Histograms
    TH1D* pt_mu1_sl_mu = new TH1D("pt_mu1_sl_mu","Leading Muon pT Distribution (Single Muon) ;Leading Muon pT (GeV);Nr. of Events",50,0,300);
    TH1D* eta_mu1_sl_mu = new TH1D("eta_mu1_sl_mu","Leading Muon eta Distribution (Single Muon) ;Leading Muon #eta; Nr. of Events",15,-3,3);
    TH1D* phi_mu1_sl_mu = new TH1D("phi_mu1_sl_mu","Leading Muon phi Distribution (Single Muon) ;Leading Muon #phi; Nr. of Events",10,-3,3);
    TH1D* pt_jet1_sl_mu = new TH1D("pt_jet1_sl_mu","Leading Jet pT Distribution (Single Muon) ;Jet 1 pT (GeV);Nr. of Events",50,0,300);
    TH1D* pt_jet2_sl_mu = new TH1D("pt_jet2_sl_mu","2nd Leading Jet pT Distribution (Single Muon) ;Jet 2 pT (GeV);Nr. of Events",50,0,300);
    TH1D* pt_jet3_sl_mu = new TH1D("pt_jet3_sl_mu","3rd Leading Jet pT Distribution (Single Muon) ;Jet 3 pT (GeV);Nr. of Events",50,0,300);
    TH1D* pt_jet4_sl_mu = new TH1D("pt_jet4_sl_mu","4th Leading Jet pT Distribution (Single Muon) ;Jet 4 pT (GeV);Nr. of Events",50,0,300);
    TH1D* pt_jet5_sl_mu = new TH1D("pt_jet5_sl_mu","5th Leading Jet pT Distribution (Single Muon) ;Jet 5 pT (GeV);Nr. of Events",50,0,300);
    TH1D* pt_jet6_sl_mu = new TH1D("pt_jet6_sl_mu","6th Leading Jet pT Distribution (Single Muon) ;Jet 6 pT (GeV);Nr. of Events",50,0,300);
    TH1D* pt_jet_all_sl_mu = new TH1D("pt_jet_all_sl_mu","Jet (all) pT Distribution (Single Muon) ;Jet pT (GeV);Nr. of Events",50,0,300);
    TH1D* csv_jet_all_sl_mu = new TH1D("csv_jet_all_sl_mu","Jet (all) CSV Distribution (Single Muon) ;Jet CSV (GeV);Nr. of Events",50,0,1);
    TH1D* ht_sl_mu = new TH1D("ht_sl_mu","HT Distribution (Single Muon) ;HT (GeV);Nr. of Events",50,0,1000);
    TH2D* njets_ht_sl_mu = new TH2D("njets_ht_2D_sl_mu","Jet Multiplicity and HT 2D Distribution (Single Muon) ;Nr. of jets ;HT (GeV)",15,0,15,50,0,1000);
    TH1D* njets_sl_mu = new TH1D("njets_sl_mu","Jet Multiplicity (Single Muon) ;Nr. of jets; Nr. of Events",15,0,15);
    TH1D* nbtags_sl_mu = new TH1D("nbtags_sl_mu","B-jet Multiplicity (Single Muon) ;Nr. of b-jets; Nr. of Events",8,0,8);
    TH1D* npv_sl_mu = new TH1D("npv_sl_mu","Nr. of Primary Vertices (Single Muon)  ;NPV; Nr. of Events",20,0,75);

    // Double Electron Histograms
    TH1D* pt_e1_di_ee = new TH1D("pt_e1_di_ee","Leading Electron pT Distribution (Double Electron) ;Leading Electron pT (GeV);Nr. of Events",50,0,300);
    TH1D* eta_e1_di_ee = new TH1D("eta_e1_di_ee","Leading Electron eta Distribution (Double Electron) ;Leading Electron #eta; Nr. of Events",15,-3,3);
    TH1D* phi_e1_di_ee = new TH1D("phi_e1_di_ee","Leading Electron phi Distribution (Double Electron) ;Leading Electron #phi; Nr. of Events",10,-3,3);
    TH1D* pt_e2_di_ee = new TH1D("pt_e2_di_ee","Sub-leading Electron pT Distribution (Double Electron) ;Sub-leading Electron pT (GeV);Nr. of Events",50,0,300);
    TH1D* eta_e2_di_ee = new TH1D("eta_e2_di_ee","Sub-leading Electron eta Distribution (Double Electron) ;Sub-leading Electron #eta; Nr. of Events",10,-3,3);
    TH1D* phi_e2_di_ee = new TH1D("phi_e2_di_ee","Sub-leading Electron phi Distribution (Double Electron) ;Sub-leading Electron #phi; Nr. of Events",10,-3,3);
    TH1D* pt_jet1_di_ee = new TH1D("pt_jet1_di_ee","Leading Jet pT Distribution (Double Electron) ;Jet 1 pT (GeV);Nr. of Events",50,0,300);
    TH1D* pt_jet2_di_ee = new TH1D("pt_jet2_di_ee","2nd Leading Jet pT Distribution (Double Electron) ;Jet 2 pT (GeV);Nr. of Events",50,0,300);
    TH1D* pt_jet3_di_ee = new TH1D("pt_jet3_di_ee","3rd Leading Jet pT Distribution (Double Electron) ;Jet 3 pT (GeV);Nr. of Events",50,0,300);
    TH1D* pt_jet4_di_ee = new TH1D("pt_jet4_di_ee","4th Leading Jet pT Distribution (Double Electron) ;Jet 4 pT (GeV);Nr. of Events",50,0,300);
    TH1D* pt_jet5_di_ee = new TH1D("pt_jet5_di_ee","5th Leading Jet pT Distribution (Double Electron) ;Jet 5 pT (GeV);Nr. of Events",50,0,300);
    TH1D* pt_jet6_di_ee = new TH1D("pt_jet6_di_ee","6th Leading Jet pT Distribution (Double Electron) ;Jet 6 pT (GeV);Nr. of Events",50,0,300);
    TH1D* pt_jet_all_di_ee = new TH1D("pt_jet_all_di_ee","Jet (all) pT Distribution (Double Electron) ;Jet pT (GeV);Nr. of Events",50,0,300);
    TH1D* csv_jet_all_di_ee = new TH1D("csv_jet_all_di_ee","Jet (all) CSV Distribution (Double Electron) ;Jet CSV (GeV);Nr. of Events",50,0,1);
    TH1D* ht_di_ee = new TH1D("ht_di_ee","HT Distribution (Double Electron) ;HT (GeV);Nr. of Events",50,0,1000);
    TH2D* njets_ht_di_ee = new TH2D("njets_ht_2D_di_ee","Jet Multiplicity and HT 2D Distribution (Double Electron) ;Nr. of jets ;HT (GeV)",15,0,15,50,0,1000);
    TH1D* njets_di_ee = new TH1D("njets_di_ee","Jet Multiplicity (Double Electron) ;Nr. of jets; Nr. of Events",15,0,15);
    TH1D* nbtags_di_ee = new TH1D("nbtags_di_ee","B-jet Multiplicity (Double Electron) ;Nr. of b-jets; Nr. of Events",8,0,8);
    TH1D* npv_di_ee = new TH1D("npv_di_ee","Nr. of Primary Vertices (Double Electron)  ;NPV; Nr. of Events",20,0,75);

    // Electron Muon Histograms
    TH1D* pt_e1_di_emu = new TH1D("pt_e1_di_emu","Electron pT Distribution (Electron Muon) ;Electron pT (GeV);Nr. of Events",50,0,300);
    TH1D* eta_e1_di_emu = new TH1D("eta_e1_di_emu","Electron eta Distribution (Electron Muon) ;Electron #eta; Nr. of Events",15,-3,3);
    TH1D* phi_e1_di_emu = new TH1D("phi_e1_di_emu","Electron phi Distribution (Electron Muon) ;Electron #phi; Nr. of Events",10,-3,3);
    TH1D* pt_mu1_di_emu = new TH1D("pt_mu1_di_emu","Muon pT Distribution (Electron Muon) ;Muon pT (GeV);Nr. of Events",50,0,300);
    TH1D* eta_mu1_di_emu = new TH1D("eta_mu1_di_emu","Muon eta Distribution (Electron Muon) ;Muon #eta; Nr. of Events",10,-3,3);
    TH1D* phi_mu1_di_emu = new TH1D("phi_mu1_di_emu","Muon phi Distribution (Electron Muon) ;Muon #phi; Nr. of Events",10,-3,3);
    TH1D* pt_jet1_di_emu = new TH1D("pt_jet1_di_emu","Leading Jet pT Distribution (Electron Muon) ;Jet 1 pT (GeV);Nr. of Events",50,0,300);
    TH1D* pt_jet2_di_emu = new TH1D("pt_jet2_di_emu","2nd Leading Jet pT Distribution (Electron Muon) ;Jet 2 pT (GeV);Nr. of Events",50,0,300);
    TH1D* pt_jet3_di_emu = new TH1D("pt_jet3_di_emu","3rd Leading Jet pT Distribution (Electron Muon) ;Jet 3 pT (GeV);Nr. of Events",50,0,300);
    TH1D* pt_jet4_di_emu = new TH1D("pt_jet4_di_emu","4th Leading Jet pT Distribution (Electron Muon) ;Jet 4 pT (GeV);Nr. of Events",50,0,300);
    TH1D* pt_jet5_di_emu = new TH1D("pt_jet5_di_emu","5th Leading Jet pT Distribution (Electron Muon) ;Jet 5 pT (GeV);Nr. of Events",50,0,300);
    TH1D* pt_jet6_di_emu = new TH1D("pt_jet6_di_emu","6th Leading Jet pT Distribution (Electron Muon) ;Jet 6 pT (GeV);Nr. of Events",50,0,300);
    TH1D* pt_jet_all_di_emu = new TH1D("pt_jet_all_di_emu","Jet (all) pT Distribution (Electron Muon) ;Jet pT (GeV);Nr. of Events",50,0,300);
    TH1D* csv_jet_all_di_emu = new TH1D("csv_jet_all_di_emu","Jet (all) CSV Distribution (Electron Muon) ;Jet CSV (GeV);Nr. of Events",50,0,1);
    TH1D* ht_di_emu = new TH1D("ht_di_emu","HT Distribution (Electron Muon) ;HT (GeV);Nr. of Events",50,0,1000);
    TH2D* njets_ht_di_emu = new TH2D("njets_ht_2D_di_emu","Jet Multiplicity and HT 2D Distribution (Electron Muon) ;Nr. of jets ;HT (GeV)",15,0,15,50,0,1000);
    TH1D* njets_di_emu = new TH1D("njets_di_emu","Jet Multiplicity (Electron Muon) ;Nr. of jets; Nr. of Events",15,0,15);
    TH1D* nbtags_di_emu = new TH1D("nbtags_di_emu","B-jet Multiplicity (Electron Muon) ;Nr. of b-jets; Nr. of Events",8,0,8);
    TH1D* npv_di_emu = new TH1D("npv_di_emu","Nr. of Primary Vertices (Electron Muon)  ;NPV; Nr. of Events",20,0,75);

    // Double Muon Histograms
    TH1D* pt_mu1_di_mumu = new TH1D("pt_mu1_di_mumu","Leading Muon pT Distribution (Double Muon) ;Leading Muon pT (GeV);Nr. of Events",50,0,300);
    TH1D* eta_mu1_di_mumu = new TH1D("eta_mu1_di_mumu","Leading Muon eta Distribution (Double Muon) ;Leading Muon #eta; Nr. of Events",15,-3,3);
    TH1D* phi_mu1_di_mumu = new TH1D("phi_mu1_di_mumu","Leading Muon phi Distribution (Double Muon) ;Leading Muon #phi; Nr. of Events",10,-3,3);
    TH1D* pt_mu2_di_mumu = new TH1D("pt_mu2_di_mumu","Sub-leading Muon pT Distribution (Double Muon) ;Sub-leading Muon pT (GeV);Nr. of Events",50,0,300);
    TH1D* eta_mu2_di_mumu = new TH1D("eta_mu2_di_mumu","Sub-leading Muon eta Distribution (Double Muon) ;Sub-leading Muon #eta; Nr. of Events",10,-3,3);
    TH1D* phi_mu2_di_mumu = new TH1D("phi_mu2_di_mumu","Sub-leading Muon phi Distribution (Double Muon) ;Sub-leading Muon #phi; Nr. of Events",10,-3,3);
    TH1D* pt_jet1_di_mumu = new TH1D("pt_jet1_di_mumu","Leading Jet pT Distribution (Double Muon) ;Jet 1 pT (GeV);Nr. of Events",50,0,300);
    TH1D* pt_jet2_di_mumu = new TH1D("pt_jet2_di_mumu","2nd Leading Jet pT Distribution (Double Muon) ;Jet 2 pT (GeV);Nr. of Events",50,0,300);
    TH1D* pt_jet3_di_mumu = new TH1D("pt_jet3_di_mumu","3rd Leading Jet pT Distribution (Double Muon) ;Jet 3 pT (GeV);Nr. of Events",50,0,300);
    TH1D* pt_jet4_di_mumu = new TH1D("pt_jet4_di_mumu","4th Leading Jet pT Distribution (Double Muon) ;Jet 4 pT (GeV);Nr. of Events",50,0,300);
    TH1D* pt_jet5_di_mumu = new TH1D("pt_jet5_di_mumu","5th Leading Jet pT Distribution (Double Muon) ;Jet 5 pT (GeV);Nr. of Events",50,0,300);
    TH1D* pt_jet6_di_mumu = new TH1D("pt_jet6_di_mumu","6th Leading Jet pT Distribution (Double Muon) ;Jet 6 pT (GeV);Nr. of Events",50,0,300);
    TH1D* pt_jet_all_di_mumu = new TH1D("pt_jet_all_di_mumu","Jet (all) pT Distribution (Double Muon) ;Jet pT (GeV);Nr. of Events",50,0,300);
    TH1D* csv_jet_all_di_mumu = new TH1D("csv_jet_all_di_mumu","Jet (all) CSV Distribution (Double Muon) ;Jet CSV (GeV);Nr. of Events",50,0,1);
    TH1D* ht_di_mumu = new TH1D("ht_di_mumu","HT Distribution (Double Muon) ;HT (GeV);Nr. of Events",50,0,1000);
    TH2D* njets_ht_di_mumu = new TH2D("njets_ht_2D_di_mumu","Jet Multiplicity and HT 2D Distribution (Double Muon) ;Nr. of jets ;HT (GeV)",15,0,15,50,0,1000);
    TH1D* njets_di_mumu = new TH1D("njets_di_mumu","Jet Multiplicity (Double Muon) ;Nr. of jets; Nr. of Events",15,0,15);
    TH1D* nbtags_di_mumu = new TH1D("nbtags_di_mumu","B-jet Multiplicity (Double Muon) ;Nr. of b-jets; Nr. of Events",8,0,8);
    TH1D* npv_di_mumu = new TH1D("npv_di_mumu","Nr. of Primary Vertices (Double Muon)  ;NPV; Nr. of Events",20,0,75);
    

	int nentries = chain->GetEntries();
	int NeventsPerJob = int( double(nentries)/double(Njobs) + 0.000001 ) + 1;

  	int firstEvent = (jobN-1)*NeventsPerJob + 1;
  	int lastEvent  = firstEvent + NeventsPerJob;
 	if( jobN==Njobs ) lastEvent = -1;
  	if( jobN==1 ) firstEvent = 0;
 
	int nEvent;
 	int n_lep, n_ele, n_mu;
	int n_jets, n_btags;

    int npv;
    int ttHf_cat;
    int ttHFGenFilter;

    std::vector<double> * mu_pt = 0;
    std::vector<double> * mu_eta = 0;
    std::vector<double> * mu_phi = 0;
    std::vector<double> * mu_E = 0;
    std::vector<int> * mu_charge = 0;
    std::vector<double> * mu_iso = 0;

    std::vector<double> * ele_pt = 0;
    std::vector<double> * ele_eta = 0;
    std::vector<double> * ele_phi = 0;
    std::vector<double> * ele_E = 0;
    std::vector<int> * ele_charge = 0;
    std::vector<double> * ele_iso = 0;

    std::vector<double> * jet_pt = 0;
    std::vector<double> * jet_eta = 0;
    std::vector<double> * jet_phi = 0;
    std::vector<double> * jet_E = 0;
    std::vector<double> * jet_CSV = 0;

    std::vector<double> * bjet_pt = 0;
    std::vector<double> * bjet_eta = 0;
    std::vector<double> * bjet_phi = 0;
    std::vector<double> * bjet_E = 0;
    std::vector<double> * bjet_CSV = 0;

    double PFMETpt;
    double PFMETphi;
    double mll;
    double ht;

    std::vector<double> * lep_sf_id = 0;
    std::vector<double> * lep_sf_iso = 0;
    double lep_sf_trig;

    double b_weight;
    double gen_weight;
    double PU_weight;
    double pdf_weight_up;
    double pdf_weight_down;
    double q2_weight_up;
    double q2_weight_down;

    std::vector<double> * genmu_pt = 0;
    std::vector<double> * genmu_eta = 0;
    std::vector<double> * genmu_phi = 0;
    std::vector<double> * genmu_E = 0;
    std::vector<int> * genmu_charge = 0;
    std::vector<int> * genmu_genid = 0;
    std::vector<int> * genmu_parentid = 0;
    std::vector<int> * genmu_grandparentid = 0;

    std::vector<double> * genele_pt = 0;
    std::vector<double> * genele_eta = 0;
    std::vector<double> * genele_phi = 0;
    std::vector<double> * genele_E = 0;
    std::vector<int> * genele_charge = 0;
    std::vector<int> * genele_genid = 0;
    std::vector<int> * genele_parentid = 0;
    std::vector<int> * genele_grandparentid = 0;

    std::vector<double> * genjet_pt = 0;
    std::vector<double> * genjet_eta = 0;
    std::vector<double> * genjet_phi = 0;
    std::vector<double> * genjet_E = 0;

    std::vector<double> * genbquarks_pt = 0;
    std::vector<double> * genbquarks_eta = 0;
    std::vector<double> * genbquarks_phi = 0;
    std::vector<int> * genbquarks_genid = 0;
    std::vector<int> * genbquarks_imm_parentid = 0;
    std::vector<int> * genbquarks_imm_daughterid = 0;
    std::vector<int> * genbquarks_parentid = 0;
    std::vector<int> * genbquarks_grandparentid = 0;


    chain->SetBranchAddress("nEvent", &nEvent );
    chain->SetBranchAddress("n_lep", &n_lep );
    chain->SetBranchAddress("n_ele", &n_ele );
    chain->SetBranchAddress("n_mu", &n_mu );
    chain->SetBranchAddress("npv", &npv );
    chain->SetBranchAddress("ttHf_cat", &ttHf_cat );
    chain->SetBranchAddress("ttHFGenFilter", &ttHFGenFilter );
    chain->SetBranchAddress("ele_pt", &ele_pt );
    chain->SetBranchAddress("ele_eta", &ele_eta );
    chain->SetBranchAddress("ele_phi", &ele_phi );
    chain->SetBranchAddress("ele_charge", &ele_charge );
    chain->SetBranchAddress("ele_iso", &ele_iso );
    chain->SetBranchAddress("mu_pt", &mu_pt );
    chain->SetBranchAddress("mu_eta", &mu_eta );
    chain->SetBranchAddress("mu_phi", &mu_phi );
    chain->SetBranchAddress("mu_charge", &mu_charge );
    chain->SetBranchAddress("mu_iso", &mu_iso );
    chain->SetBranchAddress("jet_pt", &jet_pt );
    chain->SetBranchAddress("jet_eta", &jet_eta );
    chain->SetBranchAddress("jet_phi", &jet_phi );
    chain->SetBranchAddress("jet_CSV", &jet_CSV );
    chain->SetBranchAddress("n_jets", &n_jets );
    chain->SetBranchAddress("n_btags", &n_btags );
    chain->SetBranchAddress("PFMETpt", &PFMETpt );
    chain->SetBranchAddress("PFMETphi", &PFMETphi );
    chain->SetBranchAddress("mll", &mll );
    chain->SetBranchAddress("ht", &ht );
    chain->SetBranchAddress("lep_sf_id", &lep_sf_id );
    chain->SetBranchAddress("lep_sf_iso", &lep_sf_iso );
    chain->SetBranchAddress("lep_sf_trig", &lep_sf_trig );
    chain->SetBranchAddress("b_weight", &b_weight );
    chain->SetBranchAddress("gen_weight", &gen_weight );
    chain->SetBranchAddress("PU_weight", &PU_weight );
    chain->SetBranchAddress("pdf_weight_up", &pdf_weight_up );
    chain->SetBranchAddress("pdf_weight_down", &pdf_weight_down );
    chain->SetBranchAddress("q2_weight_up", &q2_weight_up );
    chain->SetBranchAddress("q2_weight_down", &q2_weight_down );
    chain->SetBranchAddress("genmu_pt", &genmu_pt );
    chain->SetBranchAddress("genmu_eta", &genmu_eta );
    chain->SetBranchAddress("genmu_phi", &genmu_phi );
    chain->SetBranchAddress("genmu_E", &genmu_E );
    chain->SetBranchAddress("genmu_charge", &genmu_charge );
    chain->SetBranchAddress("genmu_genid", &genmu_genid );
    chain->SetBranchAddress("genmu_parentid", &genmu_parentid );
    chain->SetBranchAddress("genmu_grandparentid", &genmu_grandparentid );
    chain->SetBranchAddress("genele_pt", &genele_pt );
    chain->SetBranchAddress("genele_eta", &genele_eta );
    chain->SetBranchAddress("genele_phi", &genele_phi );
    chain->SetBranchAddress("genele_E", &genele_E );
    chain->SetBranchAddress("genele_charge", &genele_charge );
    chain->SetBranchAddress("genele_genid", &genele_genid );
    chain->SetBranchAddress("genele_parentid", &genele_parentid );
    chain->SetBranchAddress("genele_grandparentid", &genele_grandparentid );
    chain->SetBranchAddress("genjet_pt", &genjet_pt );
    chain->SetBranchAddress("genjet_eta", &genjet_eta );
    chain->SetBranchAddress("genjet_phi", &genjet_phi );
    chain->SetBranchAddress("genjet_E", &genjet_E );
    chain->SetBranchAddress("genbquarks_pt", &genbquarks_pt );
    chain->SetBranchAddress("genbquarks_eta", &genbquarks_eta );
    chain->SetBranchAddress("genbquarks_phi", &genbquarks_phi );
    chain->SetBranchAddress("genbquarks_genid", &genbquarks_genid );
    chain->SetBranchAddress("genbquarks_imm_parentid", &genbquarks_imm_parentid );
    chain->SetBranchAddress("genbquarks_imm_daughterid", &genbquarks_imm_daughterid );
    chain->SetBranchAddress("genbquarks_parentid", &genbquarks_parentid );
    chain->SetBranchAddress("genbquarks_grandparentid", &genbquarks_grandparentid );


	std::cout << "========  Starting Event Loop  ========" << std::endl;
 	for (Long64_t ievt=0; ievt<chain->GetEntries();ievt++) {   
  
  		++N_total;

	    if( ievt<firstEvent ) continue;
	    if( ievt==lastEvent ) break;

      	if( ievt==0 )        std::cout << "     Event " << ievt+1 << std::endl;
      	if( ievt%10000==0 && ievt!=1) std::cout << "           " << ievt << "\t" 
  					     << int(double(ievt-1)/double(nentries)*100) << "% done" << std::endl;
	  
	 	if( ievt==(maxNentries+1) && ievt!=0 ) break;
		
		chain->GetEntry(ievt);

        /*
        if(nEvent == 31533){
            std::cout<<nEvent<<","<<n_ele<<","<<n_mu<<","<<n_jets<<","<<n_btags<<","<<(*ele_pt)[0]<<","<<(*ele_eta)[0]<<","<<(*ele_iso)[0]<<","<<(*ele_charge)[0]<<",";
            std::cout<<(*lep_sf_id)[0]<<","<<(*lep_sf_iso)[0]<<","<<(*jet_pt)[0]<<","<<(*jet_eta)[0]<<","<<(*jet_phi)[0]<<","<<(*jet_CSV)[0]<<",";
            std::cout<<(*jet_pt)[1]<<","<<(*jet_eta)[1]<<","<<(*jet_phi)[1]<<","<<(*jet_CSV)[1]<<","<<PFMETpt<<","<<PFMETphi<<","<<mll<<",";
            std::cout<<ttHf_cat<<","<<npv<<","<<PU_weight<<","<<b_weight<<","<<pdf_weight_up<<","<<pdf_weight_down<<","<<lep_sf_trig<<"\n\n";
        }

        if(nEvent == 31544){
            std::cout<<nEvent<<","<<n_ele<<","<<n_mu<<","<<n_jets<<","<<n_btags<<","<<(*mu_pt)[0]<<","<<(*mu_eta)[0]<<","<<(*mu_iso)[0]<<","<<(*mu_charge)[0]<<",";
            std::cout<<(*lep_sf_id)[0]<<","<<(*lep_sf_iso)[0]<<","<<(*jet_pt)[0]<<","<<(*jet_eta)[0]<<","<<(*jet_phi)[0]<<","<<(*jet_CSV)[0]<<",";
            std::cout<<(*jet_pt)[1]<<","<<(*jet_eta)[1]<<","<<(*jet_phi)[1]<<","<<(*jet_CSV)[1]<<","<<PFMETpt<<","<<PFMETphi<<","<<mll<<",";
            std::cout<<ttHf_cat<<","<<npv<<","<<PU_weight<<","<<b_weight<<","<<pdf_weight_up<<","<<pdf_weight_down<<","<<lep_sf_trig<<"\n\n";
        }

        if (nEvent == 31602){
            std::cout<<nEvent<<","<<n_ele<<","<<n_mu<<","<<n_jets<<","<<n_btags<<","<<(*mu_pt)[0]<<","<<(*mu_eta)[0]<<","<<(*mu_iso)[0]<<","<<(*mu_charge)[0]<<",";
            std::cout<<(*lep_sf_id)[0]<<","<<(*lep_sf_iso)[0]<<","<<(*mu_pt)[1]<<","<<(*mu_eta)[1]<<","<<(*mu_iso)[1]<<","<<(*mu_charge)[1]<<",";
            std::cout<<(*lep_sf_id)[1]<<","<<(*lep_sf_iso)[1]<<",";
            std::cout<<(*jet_pt)[0]<<","<<(*jet_eta)[0]<<","<<(*jet_phi)[0]<<","<<(*jet_CSV)[0]<<",";
            std::cout<<(*jet_pt)[1]<<","<<(*jet_eta)[1]<<","<<(*jet_phi)[1]<<","<<(*jet_CSV)[1]<<","<<PFMETpt<<","<<PFMETphi<<","<<mll<<",";
            std::cout<<ttHf_cat<<","<<npv<<","<<PU_weight<<","<<b_weight<<","<<pdf_weight_up<<","<<pdf_weight_down<<","<<lep_sf_trig<<"\n\n";
        }

        if(nEvent == 31650){
            std::cout<<nEvent<<","<<n_ele<<","<<n_mu<<","<<n_jets<<","<<n_btags<<","<<(*ele_pt)[0]<<","<<(*ele_eta)[0]<<","<<(*ele_iso)[0]<<","<<(*ele_charge)[0]<<",";
            std::cout<<(*lep_sf_id)[0]<<","<<(*lep_sf_iso)[0]<<","<<(*mu_pt)[0]<<","<<(*mu_eta)[0]<<","<<(*mu_iso)[0]<<","<<(*mu_charge)[0]<<",";
            std::cout<<(*lep_sf_id)[1]<<","<<(*lep_sf_iso)[1]<<",";
            std::cout<<(*jet_pt)[0]<<","<<(*jet_eta)[0]<<","<<(*jet_phi)[0]<<","<<(*jet_CSV)[0]<<",";
            std::cout<<(*jet_pt)[1]<<","<<(*jet_eta)[1]<<","<<(*jet_phi)[1]<<","<<(*jet_CSV)[1]<<","<<PFMETpt<<","<<PFMETphi<<","<<mll<<",";
            std::cout<<ttHf_cat<<","<<npv<<","<<PU_weight<<","<<b_weight<<","<<pdf_weight_up<<","<<pdf_weight_down<<","<<lep_sf_trig<<"\n\n";
        }

        if(nEvent == 1970314){
            std::cout<<nEvent<<","<<n_ele<<","<<n_mu<<","<<n_jets<<","<<n_btags<<","<<(*ele_pt)[0]<<","<<(*ele_eta)[0]<<","<<(*ele_iso)[0]<<","<<(*ele_charge)[0]<<",";
            std::cout<<(*lep_sf_id)[0]<<","<<(*lep_sf_iso)[0]<<","<<(*ele_pt)[1]<<","<<(*ele_eta)[1]<<","<<(*ele_iso)[1]<<","<<(*ele_charge)[1]<<",";
            std::cout<<(*lep_sf_id)[1]<<","<<(*lep_sf_iso)[1]<<",";
            std::cout<<(*jet_pt)[0]<<","<<(*jet_eta)[0]<<","<<(*jet_phi)[0]<<","<<(*jet_CSV)[0]<<",";
            std::cout<<(*jet_pt)[1]<<","<<(*jet_eta)[1]<<","<<(*jet_phi)[1]<<","<<(*jet_CSV)[1]<<","<<PFMETpt<<","<<PFMETphi<<","<<mll<<",";
            std::cout<<ttHf_cat<<","<<npv<<","<<PU_weight<<","<<b_weight<<","<<pdf_weight_up<<","<<pdf_weight_down<<","<<lep_sf_trig<<"\n\n";
        }
        */

        
        sum_gen_weights = sum_gen_weights + gen_weight;

        // for taking tt jets hf filtered events only
        if(ttHFGenFilter!=1)
              continue;

        //if(n_btags<4)
        //    continue;
        

        //for(unsigned int i=0; i< (*genbquarks_pt).size(); i++){
        //    std::cout<<(*genbquarks_genid)[i]<<"  "<<(*genbquarks_imm_parentid)[i]<<"  "<<(*genbquarks_imm_daughterid)[i]<<"  "<<(*genbquarks_parentid)[i]<<"  "<<(*genbquarks_grandparentid)[i]<<"  "<<(*genbquarks_pt)[i]<<"  "<<(*genbquarks_eta)[i]<<"  "<<(*genbquarks_phi)[i]<<"\n";
        //}

        // Checking at least 4 b-quarks at Gen Level
        int n_bquarks = 0;
        for(unsigned int i=0; i< (*genbquarks_pt).size(); i++){
            if((*genbquarks_imm_parentid)[i] != -99 && (*genbquarks_imm_daughterid)[i] != -99){
                if((*genbquarks_imm_daughterid)[i] != 5 && (*genbquarks_imm_daughterid)[i] != -5){
                    if( (*genbquarks_pt)[i] > 30 )
                        n_bquarks++;
                }
            }
        }
        if(n_bquarks < 4)
            continue;

         ++N_4b_pass;

        // Checking ttjets category = ttbb
        double mod_ttHFcat = ttHf_cat%100;
        if( (mod_ttHFcat < 53) || (mod_ttHFcat > 55) )
            continue;

        ++N_ttbb_pass;

        // Event weight and SFs
        double tot_weight = 1;
        if (n_lep == 1 || n_lep == 2){
            tot_weight = tot_weight * (*lep_sf_id)[0] * (*lep_sf_id)[1] * (*lep_sf_iso)[0] * (*lep_sf_iso)[1];
            tot_weight = tot_weight * lep_sf_trig * gen_weight * PU_weight * b_weight;
        }
        

        // Single Lepton events
		if (n_lep == 1) {
			++N_SL;
			if (n_ele == 1) {
				++N_e;
				pt_e1_sl_e->Fill((*ele_pt)[0], tot_weight);
				eta_e1_sl_e->Fill((*ele_eta)[0], tot_weight);
                phi_e1_sl_e->Fill((*ele_phi)[0], tot_weight);
                pt_jet1_sl_e->Fill((*jet_pt)[0], tot_weight);
                pt_jet2_sl_e->Fill((*jet_pt)[1], tot_weight);
                pt_jet3_sl_e->Fill((*jet_pt)[2], tot_weight);
                pt_jet4_sl_e->Fill((*jet_pt)[3], tot_weight);
                pt_jet5_sl_e->Fill((*jet_pt)[4], tot_weight);
                pt_jet6_sl_e->Fill((*jet_pt)[5], tot_weight);

                for(int i=0; i<n_jets; i++){
                    pt_jet_all_sl_e->Fill((*jet_pt)[i], tot_weight);
                    csv_jet_all_sl_e->Fill((*jet_CSV)[i], tot_weight);
                }

                ht_sl_e->Fill(ht, tot_weight);
                njets_ht_sl_e->Fill(n_jets, ht, tot_weight);
                njets_sl_e->Fill(n_jets, tot_weight);
                nbtags_sl_e->Fill(n_btags, tot_weight);
                npv_sl_e->Fill(npv, tot_weight);
			}
			else if (n_mu == 1) {
				++N_mu;
                pt_mu1_sl_mu->Fill((*mu_pt)[0], tot_weight);
                eta_mu1_sl_mu->Fill((*mu_eta)[0], tot_weight);
                phi_mu1_sl_mu->Fill((*mu_phi)[0], tot_weight);
                pt_jet1_sl_mu->Fill((*jet_pt)[0], tot_weight);
                pt_jet2_sl_mu->Fill((*jet_pt)[1], tot_weight);
                pt_jet3_sl_mu->Fill((*jet_pt)[2], tot_weight);
                pt_jet4_sl_mu->Fill((*jet_pt)[3], tot_weight);
                pt_jet5_sl_mu->Fill((*jet_pt)[4], tot_weight);
                pt_jet6_sl_mu->Fill((*jet_pt)[5], tot_weight);

                for(int i=0; i<n_jets; i++){
                    pt_jet_all_sl_mu->Fill((*jet_pt)[i], tot_weight);
                    csv_jet_all_sl_mu->Fill((*jet_CSV)[i], tot_weight);
                }

                ht_sl_mu->Fill(ht, tot_weight);
                njets_ht_sl_mu->Fill(n_jets, ht, tot_weight);
                njets_sl_mu->Fill(n_jets, tot_weight);
                nbtags_sl_mu->Fill(n_btags, tot_weight);
                npv_sl_mu->Fill(npv, tot_weight);
			}
		}

        // Dilepton events
        else if (n_lep == 2) {
            ++N_DL;
            if (n_ele == 2) {
                ++N_ee;
                pt_e1_di_ee->Fill((*ele_pt)[0], tot_weight);
                eta_e1_di_ee->Fill((*ele_eta)[0], tot_weight);
                phi_e1_di_ee->Fill((*ele_phi)[0], tot_weight);
                pt_e2_di_ee->Fill((*ele_pt)[1], tot_weight);
                eta_e2_di_ee->Fill((*ele_eta)[1], tot_weight);
                phi_e2_di_ee->Fill((*ele_phi)[1], tot_weight);
                pt_jet1_di_ee->Fill((*jet_pt)[0], tot_weight);
                pt_jet2_di_ee->Fill((*jet_pt)[1], tot_weight);
                pt_jet3_di_ee->Fill((*jet_pt)[2], tot_weight);
                pt_jet4_di_ee->Fill((*jet_pt)[3], tot_weight);
                pt_jet5_di_ee->Fill((*jet_pt)[4], tot_weight);
                pt_jet6_di_ee->Fill((*jet_pt)[5], tot_weight);

                for(int i=0; i<n_jets; i++){
                    pt_jet_all_di_ee->Fill((*jet_pt)[i], tot_weight);
                    csv_jet_all_di_ee->Fill((*jet_CSV)[i], tot_weight);
                }

                ht_di_ee->Fill(ht, tot_weight);
                njets_ht_di_ee->Fill(n_jets, ht, tot_weight);
                njets_di_ee->Fill(n_jets, tot_weight);
                nbtags_di_ee->Fill(n_btags, tot_weight);
                npv_di_ee->Fill(npv, tot_weight);
            }
            else if (n_mu == 2){
                ++N_mumu;
                pt_mu1_di_mumu->Fill((*mu_pt)[0], tot_weight);
                eta_mu1_di_mumu->Fill((*mu_eta)[0], tot_weight);
                phi_mu1_di_mumu->Fill((*mu_phi)[0], tot_weight);
                pt_mu2_di_mumu->Fill((*mu_pt)[1], tot_weight);
                eta_mu2_di_mumu->Fill((*mu_eta)[1], tot_weight);
                phi_mu2_di_mumu->Fill((*mu_phi)[1], tot_weight);
                pt_jet1_di_mumu->Fill((*jet_pt)[0], tot_weight);
                pt_jet2_di_mumu->Fill((*jet_pt)[1], tot_weight);
                pt_jet3_di_mumu->Fill((*jet_pt)[2], tot_weight);
                pt_jet4_di_mumu->Fill((*jet_pt)[3], tot_weight);
                pt_jet5_di_mumu->Fill((*jet_pt)[4], tot_weight);
                pt_jet6_di_mumu->Fill((*jet_pt)[5], tot_weight);

                for(int i=0; i<n_jets; i++){
                    pt_jet_all_di_mumu->Fill((*jet_pt)[i], tot_weight);
                    csv_jet_all_di_mumu->Fill((*jet_CSV)[i], tot_weight);
                }

                ht_di_mumu->Fill(ht, tot_weight);
                njets_ht_di_mumu->Fill(n_jets, ht, tot_weight);
                njets_di_mumu->Fill(n_jets, tot_weight);
                nbtags_di_mumu->Fill(n_btags, tot_weight);
                npv_di_mumu->Fill(npv, tot_weight);

            }
            else {
                ++N_emu;
                pt_e1_di_emu->Fill((*ele_pt)[0], tot_weight);
                eta_e1_di_emu->Fill((*ele_eta)[0], tot_weight);
                phi_e1_di_emu->Fill((*ele_phi)[0], tot_weight);
                pt_mu1_di_emu->Fill((*mu_pt)[0], tot_weight);
                eta_mu1_di_emu->Fill((*mu_eta)[0], tot_weight);
                phi_mu1_di_emu->Fill((*mu_phi)[0], tot_weight);
                pt_jet1_di_emu->Fill((*jet_pt)[0], tot_weight);
                pt_jet2_di_emu->Fill((*jet_pt)[1], tot_weight);
                pt_jet3_di_emu->Fill((*jet_pt)[2], tot_weight);
                pt_jet4_di_emu->Fill((*jet_pt)[3], tot_weight);
                pt_jet5_di_emu->Fill((*jet_pt)[4], tot_weight);
                pt_jet6_di_emu->Fill((*jet_pt)[5], tot_weight);

                for(int i=0; i<n_jets; i++){
                    pt_jet_all_di_emu->Fill((*jet_pt)[i], tot_weight);
                    csv_jet_all_di_emu->Fill((*jet_CSV)[i], tot_weight);
                }

                ht_di_emu->Fill(ht, tot_weight);
                njets_ht_di_emu->Fill(n_jets, ht, tot_weight);
                njets_di_emu->Fill(n_jets, tot_weight);
                nbtags_di_emu->Fill(n_btags, tot_weight);
                npv_di_emu->Fill(npv, tot_weight);
            }
        }
		
	}   

    N_sel = N_SL + N_DL;

    std::cout << " Done! " << std::endl;
 	std::cout<<"**********************************************************************************************\n";
    std::cout<<"Total No. of Events : "<<N_total<<"\n";
    std::cout<<"Total No. of Events having more than 4 b-quarks : "<<N_4b_pass<<" , "<<((double)N_4b_pass/N_total)*100<<"\n";
    std::cout<<"Total No. of Events having more than 4 b-quarks + in ttbb category : "<<N_ttbb_pass<<" , "<<((double)N_ttbb_pass/N_total)*100<<"\n";
    std::cout<<"Total No. of Events having more than 4 b-quarks + in ttbb category + passing Selection : "<<N_sel<<" , "<<((double)N_sel/N_total)*100<<"\n";
    std::cout<<"**********************************************************************************************\n";
    std::cout<<"No. of Single Electron Events: "<<N_e<<"\n";
    std::cout<<"No. of Single Muon Events: "<<N_mu<<"\n";
    std::cout<<"No. of Double Electron Events: "<<N_ee<<"\n";
    std::cout<<"No. of Electron Muon Events: "<<N_emu<<"\n";
    std::cout<<"No. of Double Muon Events: "<<N_mumu<<"\n";
    std::cout<<"Total No. of Single Lepton Events: "<<N_SL<<"\n";
    std::cout<<"Total No. of Dilepton Events: "<<N_DL<<"\n";
    std::cout<<"**********************************************************************************************\n";
    std::cout<<"Sum of Generator weights of all Events: "<<sum_gen_weights<<"\n";
    std::cout<<"**********************************************************************************************\n";
    histofile.Write();
    histofile.Close();
    
    return;
}


