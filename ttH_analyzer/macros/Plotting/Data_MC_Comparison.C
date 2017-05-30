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
    int N_SL = 0;
    int N_DL = 0;
	int N_e = 0;
	int N_mu = 0;
    int N_ee = 0;
    int N_emu = 0;
    int N_mumu = 0;

	std::string treefilename = "ttHbbNtuple.root";
	
	std::string histofilename = "Data_MC_Comparison.root";
	
	TChain *chain = new TChain("ttHbb/eventTree");

    chain->Add(treefilename.c_str());

    TFile histofile(histofilename.c_str(),"recreate");
    histofile.cd();
    
	TH1::SetDefaultSumw2();

    // Single Electron Histograms
	TH1D* pt_e1_e = new TH1D("pt_e1_e","Leading Electron pT Distribution (Single Electron) ;Leading Electron pT (GeV);Nr. of Events",150,0,300);
	TH1D* eta_e1_e = new TH1D("eta_e1_e","Leading Electron eta Distribution (Single Electron) ;Leading Electron #eta; Nr. of Events",30,-3,3);
	TH1D* phi_e1_e = new TH1D("phi_e1_e","Leading Electron phi Distribution (Single Electron) ;Leading Electron #phi; Nr. of Events",30,-3,3);
    TH1D* pt_jet1_e = new TH1D("pt_jet1_e","Leading Jet pT Distribution (Single Electron) ;Jet 1 pT (GeV);Nr. of Events",150,0,300);
    TH1D* pt_jet2_e = new TH1D("pt_jet2_e","2nd Leading Jet pT Distribution (Single Electron) ;Jet 2 pT (GeV);Nr. of Events",150,0,300);
    TH1D* pt_jet3_e = new TH1D("pt_jet3_e","3rd Leading Jet pT Distribution (Single Electron) ;Jet 3 pT (GeV);Nr. of Events",150,0,300);
    TH1D* pt_jet4_e = new TH1D("pt_jet4_e","4th Leading Jet pT Distribution (Single Electron) ;Jet 4 pT (GeV);Nr. of Events",150,0,300);
    TH1D* pt_jet5_e = new TH1D("pt_jet5_e","5th Leading Jet pT Distribution (Single Electron) ;Jet 5 pT (GeV);Nr. of Events",150,0,300);
    TH1D* pt_jet6_e = new TH1D("pt_jet6_e","6th Leading Jet pT Distribution (Single Electron) ;Jet 6 pT (GeV);Nr. of Events",150,0,300);
    TH1D* pt_jet_all_e = new TH1D("pt_jet_all_e","Jet (all) pT Distribution (Single Electron) ;Jet pT (GeV);Nr. of Events",150,0,300);
    TH1D* csv_jet_all_e = new TH1D("csv_jet_all_e","Jet (all) CSV Distribution (Single Electron) ;Jet CSV (GeV);Nr. of Events",100,0,1);
    TH1D* ht_e = new TH1D("ht_e","HT Distribution (Single Electron) ;HT (GeV);Nr. of Events",200,0,1000);
    TH1D* njets_e = new TH1D("njets_e","Jet Multiplicity (Single Electron) ;Nr. of jets; Nr. of Events",15,0,15);
    TH1D* nbtags_e = new TH1D("nbtags_e","B-jet Multiplicity (Single Electron) ;Nr. of b-jets; Nr. of Events",8,0,8);
    TH1D* npv_e = new TH1D("npv_e","Nr. of Primary Vertices (Single Electron)  ;NPV; Nr. of Events",75,0,75);

    // Single Muon Histograms
	TH1D* pt_mu1_mu = new TH1D("pt_mu1_mu","Leading Muon pT Distribution (Single Muon) ;Leading Muon pT (GeV);Nr. of Events",150,0,300);
	TH1D* eta_mu1_mu = new TH1D("eta_mu1_mu","Leading Muon eta Distribution (Single Muon) ;Leading Muon #eta; Nr. of Events",30,-3,3);
	TH1D* phi_mu1_mu = new TH1D("phi_mu1_mu","Leading Muon phi Distribution (Single Muon) ;Leading Muon #phi; Nr. of Events",30,-3,3);
    TH1D* pt_jet1_mu = new TH1D("pt_jet1_mu","Leading Jet pT Distribution (Single Muon) ;Jet 1 pT (GeV);Nr. of Events",150,0,300);
    TH1D* pt_jet2__mu = new TH1D("pt_jet2_mu","2nd Leading Jet pT Distribution (Single Muon) ;Jet 2 pT (GeV);Nr. of Events",150,0,300);
    TH1D* pt_jet3_mu = new TH1D("pt_jet3_mu","3rd Leading Jet pT Distribution (Single Muon) ;Jet 3 pT (GeV);Nr. of Events",150,0,300);
    TH1D* pt_jet4_mu = new TH1D("pt_jet4_mu","4th Leading Jet pT Distribution (Single Muon) ;Jet 4 pT (GeV);Nr. of Events",150,0,300);
    TH1D* pt_jet5_mu = new TH1D("pt_jet5_mu","5th Leading Jet pT Distribution (Single Muon) ;Jet 5 pT (GeV);Nr. of Events",150,0,300);
    TH1D* pt_jet6_mu = new TH1D("pt_jet6_mu","6th Leading Jet pT Distribution (Single Muon) ;Jet 6 pT (GeV);Nr. of Events",150,0,300);
    TH1D* pt_jet_all_mu = new TH1D("pt_jet_all_mu","Jet (all) pT Distribution (Single Muon) ;Jet pT (GeV);Nr. of Events",150,0,300);
    TH1D* csv_jet_all_mu = new TH1D("csv_jet_all_mu","Jet (all) CSV Distribution (Single Muon) ;Jet CSV (GeV);Nr. of Events",100,0,1);
    TH1D* ht_mu = new TH1D("ht_mu","HT Distribution (Single Muon) ;HT (GeV);Nr. of Events",200,0,1000);
    TH1D* njets_mu = new TH1D("njets_mu","Jet Multiplicity (Single Muon) ;Nr. of jets; Nr. of Events",15,0,15);
    TH1D* nbtags_mu = new TH1D("nbtags_mu","B-jet Multiplicity (Single Muon) ;Nr. of b-jets; Nr. of Events",8,0,8);
    TH1D* npv_mu = new TH1D("npv_mu","Nr. of Primary Vertices (Single Muon)  ;NPV; Nr. of Events",75,0,75);

    // Double Electron Histograms
    TH1D* pt_e1_ee = new TH1D("pt_e1_ee","Leading Electron pT Distribution (Double Electron) ;Leading Electron pT (GeV);Nr. of Events",150,0,300);
    TH1D* eta_e1_ee = new TH1D("eta_e1_ee","Leading Electron eta Distribution (Double Electron) ;Leading Electron #eta; Nr. of Events",30,-3,3);
    TH1D* phi_e1_ee = new TH1D("phi_e1_ee","Leading Electron phi Distribution (Double Electron) ;Leading Electron #phi; Nr. of Events",30,-3,3);
    TH1D* pt_e2_ee = new TH1D("pt_e2_ee","Sub-leading Electron pT Distribution (Double Electron) ;Sub-leading Electron pT (GeV);Nr. of Events",150,0,300);
    TH1D* eta_e2_ee = new TH1D("eta_e2_ee","Sub-leading Electron eta Distribution (Double Electron) ;Sub-leading Electron #eta; Nr. of Events",30,-3,3);
    TH1D* phi_e2_ee = new TH1D("phi_e2_ee","Sub-leading Electron phi Distribution (Double Electron) ;Sub-leading Electron #phi; Nr. of Events",30,-3,3);
    TH1D* pt_jet1_ee = new TH1D("pt_jet1_ee","Leading Jet pT Distribution (Double Electron) ;Jet 1 pT (GeV);Nr. of Events",150,0,300);
    TH1D* pt_jet2_ee = new TH1D("pt_jet2_ee","2nd Leading Jet pT Distribution (Double Electron) ;Jet 2 pT (GeV);Nr. of Events",150,0,300);
    TH1D* pt_jet3_ee = new TH1D("pt_jet3_ee","3rd Leading Jet pT Distribution (Double Electron) ;Jet 3 pT (GeV);Nr. of Events",150,0,300);
    TH1D* pt_jet4_ee = new TH1D("pt_jet4_ee","4th Leading Jet pT Distribution (Double Electron) ;Jet 4 pT (GeV);Nr. of Events",150,0,300);
    TH1D* pt_jet5_ee = new TH1D("pt_jet5_ee","5th Leading Jet pT Distribution (Double Electron) ;Jet 5 pT (GeV);Nr. of Events",150,0,300);
    TH1D* pt_jet6_ee = new TH1D("pt_jet6_ee","6th Leading Jet pT Distribution (Double Electron) ;Jet 6 pT (GeV);Nr. of Events",150,0,300);
    TH1D* pt_jet_all_ee = new TH1D("pt_jet_all_ee","Jet (all) pT Distribution (Double Electron) ;Jet pT (GeV);Nr. of Events",150,0,300);
    TH1D* csv_jet_all_ee = new TH1D("csv_jet_all_ee","Jet (all) CSV Distribution (Double Electron) ;Jet CSV (GeV);Nr. of Events",100,0,1);
    TH1D* ht_ee = new TH1D("ht_ee","HT Distribution (Double Electron) ;HT (GeV);Nr. of Events",200,0,1000);
    TH1D* njets_ee = new TH1D("njets_ee","Jet Multiplicity (Double Electron) ;Nr. of jets; Nr. of Events",15,0,15);
    TH1D* nbtags_ee = new TH1D("nbtags_ee","B-jet Multiplicity (Double Electron) ;Nr. of b-jets; Nr. of Events",8,0,8);
    TH1D* npv_ee = new TH1D("npv_ee","Nr. of Primary Vertices (Double Electron)  ;NPV; Nr. of Events",75,0,75);

    // Electron Muon Histograms
    TH1D* pt_e1_emu = new TH1D("pt_e1_emu","Electron pT Distribution (Electron Muon) ;Electron pT (GeV);Nr. of Events",150,0,300);
    TH1D* eta_e1_emu = new TH1D("eta_e1_emu","Electron eta Distribution (Electron Muon) ;Electron #eta; Nr. of Events",30,-3,3);
    TH1D* phi_e1_emu = new TH1D("phi_e1_emu","Electron phi Distribution (Electron Muon) ;Electron #phi; Nr. of Events",30,-3,3);
    TH1D* pt_mu1_emu = new TH1D("pt_mu1_emu","Muon pT Distribution (Electron Muon) ;Muon pT (GeV);Nr. of Events",150,0,300);
    TH1D* eta_mu1_emu = new TH1D("eta_mu1_emu","Muon eta Distribution (Electron Muon) ;Muon #eta; Nr. of Events",30,-3,3);
    TH1D* phi_mu1_emu = new TH1D("phi_mu1_emu","Muon phi Distribution (Electron Muon) ;Muon #phi; Nr. of Events",30,-3,3);
    TH1D* pt_jet1_emu = new TH1D("pt_jet1_emu","Leading Jet pT Distribution (Electron Muon) ;Jet 1 pT (GeV);Nr. of Events",150,0,300);
    TH1D* pt_jet2_emu = new TH1D("pt_jet2_emu","2nd Leading Jet pT Distribution (Electron Muon) ;Jet 2 pT (GeV);Nr. of Events",150,0,300);
    TH1D* pt_jet3_emu = new TH1D("pt_jet3_emu","3rd Leading Jet pT Distribution (Electron Muon) ;Jet 3 pT (GeV);Nr. of Events",150,0,300);
    TH1D* pt_jet4_emu = new TH1D("pt_jet4_emu","4th Leading Jet pT Distribution (Electron Muon) ;Jet 4 pT (GeV);Nr. of Events",150,0,300);
    TH1D* pt_jet5_emu = new TH1D("pt_jet5_emu","5th Leading Jet pT Distribution (Electron Muon) ;Jet 5 pT (GeV);Nr. of Events",150,0,300);
    TH1D* pt_jet6_emu = new TH1D("pt_jet6_emu","6th Leading Jet pT Distribution (Electron Muon) ;Jet 6 pT (GeV);Nr. of Events",150,0,300);
    TH1D* pt_jet_all_emu = new TH1D("pt_jet_all_emu","Jet (all) pT Distribution (Electron Muon) ;Jet pT (GeV);Nr. of Events",150,0,300);
    TH1D* csv_jet_all_emu = new TH1D("csv_jet_all_emu","Jet (all) CSV Distribution (Electron Muon) ;Jet CSV (GeV);Nr. of Events",100,0,1);
    TH1D* ht_emu = new TH1D("ht_emu","HT Distribution (Electron Muon) ;HT (GeV);Nr. of Events",200,0,1000);
    TH1D* njets_emu = new TH1D("njets_emu","Jet Multiplicity (Electron Muon) ;Nr. of jets; Nr. of Events",15,0,15);
    TH1D* nbtags_emu = new TH1D("nbtags_emu","B-jet Multiplicity (Electron Muon) ;Nr. of b-jets; Nr. of Events",8,0,8);
    TH1D* npv_emu = new TH1D("npv_emu","Nr. of Primary Vertices (Electron Muon)  ;NPV; Nr. of Events",75,0,75);

    // Double Muon Histograms
    TH1D* pt_mu1_mumu = new TH1D("pt_mu1_mumu","Leading Muon pT Distribution (Double Muon) ;Leading Muon pT (GeV);Nr. of Events",150,0,300);
    TH1D* eta_mu1_mumu = new TH1D("eta_mu1_mumu","Leading Muon eta Distribution (Double Muon) ;Leading Muon #eta; Nr. of Events",30,-3,3);
    TH1D* phi_mu1_mumu = new TH1D("phi_mu1_mumu","Leading Muon phi Distribution (Double Muon) ;Leading Muon #phi; Nr. of Events",30,-3,3);
    TH1D* pt_mu2_mumu = new TH1D("pt_mu2_mumu","Sub-leading Muon pT Distribution (Double Muon) ;Sub-leading Muon pT (GeV);Nr. of Events",150,0,300);
    TH1D* eta_mu2_mumu = new TH1D("eta_mu2_mumu","Sub-leading Muon eta Distribution (Double Muon) ;Sub-leading Muon #eta; Nr. of Events",30,-3,3);
    TH1D* phi_mu2_mumu = new TH1D("phi_mu2_mumu","Sub-leading Muon phi Distribution (Double Muon) ;Sub-leading Muon #phi; Nr. of Events",30,-3,3);
    TH1D* pt_jet1_mumu = new TH1D("pt_jet1_mumu","Leading Jet pT Distribution (Double Muon) ;Jet 1 pT (GeV);Nr. of Events",150,0,300);
    TH1D* pt_jet2_mumu = new TH1D("pt_jet2_mumu","2nd Leading Jet pT Distribution (Double Muon) ;Jet 2 pT (GeV);Nr. of Events",150,0,300);
    TH1D* pt_jet3_mumu = new TH1D("pt_jet3_mumu","3rd Leading Jet pT Distribution (Double Muon) ;Jet 3 pT (GeV);Nr. of Events",150,0,300);
    TH1D* pt_jet4_mumu = new TH1D("pt_jet4_mumu","4th Leading Jet pT Distribution (Double Muon) ;Jet 4 pT (GeV);Nr. of Events",150,0,300);
    TH1D* pt_jet5_mumu = new TH1D("pt_jet5_mumu","5th Leading Jet pT Distribution (Double Muon) ;Jet 5 pT (GeV);Nr. of Events",150,0,300);
    TH1D* pt_jet6_mumu = new TH1D("pt_jet6_mumu","6th Leading Jet pT Distribution (Double Muon) ;Jet 6 pT (GeV);Nr. of Events",150,0,300);
    TH1D* pt_jet_all_mumu = new TH1D("pt_jet_all_mumu","Jet (all) pT Distribution (Double Muon) ;Jet pT (GeV);Nr. of Events",150,0,300);
    TH1D* csv_jet_all_mumu = new TH1D("csv_jet_all_mumu","Jet (all) CSV Distribution (Double Muon) ;Jet CSV (GeV);Nr. of Events",100,0,1);
    TH1D* ht_mumu = new TH1D("ht_mumu","HT Distribution (Double Muon) ;HT (GeV);Nr. of Events",200,0,1000);
    TH1D* njets_mumu = new TH1D("njets_mumu","Jet Multiplicity (Double Muon) ;Nr. of jets; Nr. of Events",15,0,15);
    TH1D* nbtags_mumu = new TH1D("nbtags_mumu","B-jet Multiplicity (Double Muon) ;Nr. of b-jets; Nr. of Events",8,0,8);
    TH1D* npv_mumu = new TH1D("npv_mumu","Nr. of Primary Vertices (Double Muon)  ;NPV; Nr. of Events",75,0,75);


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


    chain->SetBranchAddress("nEvent", &nEvent );
    chain->SetBranchAddress("n_lep", &n_lep );
    chain->SetBranchAddress("n_ele", &n_ele );
    chain->SetBranchAddress("n_mu", &n_mu );
    chain->SetBranchAddress("npv", &npv );
    chain->SetBranchAddress("ttHf_cat", &ttHf_cat );
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

		// Single Lepton events
		if (n_lep == 1) {
			++N_SL;	
			//N_jets->Fill(n_jets);
			//N_btags->Fill(n_btags);
			if (n_ele == 1) {
				++N_e;
				//pt_e1->Fill(ele_pt[0]);
				//eta_e1->Fill(ele_eta[0]);
				//phi_e1 ->Fill(ele_phi[0]);
			}
			else if (n_mu == 1) {
				++N_mu;
				//pt_mu1->Fill(mu_pt[0]);
				//eta_mu1->Fill(mu_eta[0]);
				//phi_mu1 ->Fill(mu_phi[0]);
			}
		}

        else if (n_lep == 2) {
            if (n_ele == 2)
                ++N_ee;
            else if (n_mu == 2)
                ++N_mumu;
            else
                ++N_emu;
        }
		
	}   

    std::cout << " Done! " << std::endl;
 	std::cout<<"**********************************************************************************************\n";
    std::cout<<"Total No. of events : "<<N_total<<"\n";   
    std::cout<<"No. of Single Electron Events: "<<N_e<<"\n";
    std::cout<<"No. of Single Muon Events: "<<N_mu<<"\n";
    std::cout<<"No. of Double Electron Events: "<<N_ee<<"\n";
    std::cout<<"No. of Electron Muon Events: "<<N_emu<<"\n";
    std::cout<<"No. of Double Muon Events: "<<N_mumu<<"\n";
    std::cout<<"Total No. of Single Lepton Events: "<<N_SL<<"\n";
    std::cout<<"**********************************************************************************************\n";
    histofile.Write();
    histofile.Close();
    
    return;
}


