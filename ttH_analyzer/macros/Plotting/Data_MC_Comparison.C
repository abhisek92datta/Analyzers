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
#include "Math/Interpolator.h"

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

	TH1D* pt_e1 = new TH1D("pt_e1","Leading Electron pT Distribution ;Electron pT (GeV);Nr. of Events",150,0,300);
	TH1D* eta_e1 = new TH1D("eta_e1","Leading Electron eta Distribution ;Electron #eta; Nr. of Events",30,-3,3);
	TH1D* phi_e1 = new TH1D("phi_e1","Leading Electron phi Distribution ;Electron #phi; Nr. of Events",30,-3,3);
	TH1D* pt_mu1 = new TH1D("pt_mu1","Leading Muon pT Distribution ;Muon pT (GeV);Nr. of Events",150,0,300);
	TH1D* eta_mu1 = new TH1D("eta_mu1","Leading Muon eta Distribution ;Muon #eta; Nr. of Events",30,-3,3);
	TH1D* phi_mu1 = new TH1D("phi_mu1","Leading Muon phi Distribution ;Muon #phi; Nr. of Events",30,-3,3);
    TH1D* pt_jet1 = new TH1D("pt_jet1","Leading Jet pT Distribution ;Jet pT (GeV);Nr. of Events",150,0,300);
	TH1D* N_jets = new TH1D("N_jets","Jet Multiplicity ;Nr. of jets; Nr. of Events",15,0,15);
	TH1D* N_btags = new TH1D("N_btags","B-jet Multiplicity ;Nr. of b-jets; Nr. of Events",8,0,8);

	int nentries = chain->GetEntries();
	int NeventsPerJob = int( double(nentries)/double(Njobs) + 0.000001 ) + 1;

  	int firstEvent = (jobN-1)*NeventsPerJob + 1;
  	int lastEvent  = firstEvent + NeventsPerJob;
 	if( jobN==Njobs ) lastEvent = -1;
  	if( jobN==1 ) firstEvent = 0;
 
	int nEvent;
 	int n_lep, n_ele, n_mu;
	int n_jets, n_btags;

    std::vector<double> mu_pt;
    std::vector<double> mu_eta;
    std::vector<double> mu_phi;
    std::vector<double> mu_E;
    std::vector<int> mu_charge;
    std::vector<double> mu_iso;

    std::vector<double> ele_pt;
    std::vector<double> ele_eta;
    std::vector<double> ele_phi;
    std::vector<double> ele_E;
    std::vector<int> ele_charge;
    std::vector<double> ele_iso;

    std::vector<double> jet_pt;
    std::vector<double> jet_eta;
    std::vector<double> jet_phi;
    std::vector<double> jet_E;
    std::vector<double> jet_CSV;

    std::vector<double> bjet_pt;
    std::vector<double> bjet_eta;
    std::vector<double> bjet_phi;
    std::vector<double> bjet_E;
    std::vector<double> bjet_CSV;

    double PFMETpt;
    double PFMETphi;

    std::vector<double> lep_sf_id;
    std::vector<double> lep_sf_iso;
    double lep_sf_trig;

    double b_weight;
    double gen_weight;
    double PU_weight;
    double pdf_weight_up;
    double pdf_weight_down;
    double q2_weight_up;
    double q2_weight_down;



	// to ensure proper reading from root file in the event loop	
	chain->GetEntry(0);
    chain->SetBranchAddress("nEvent", &nEvent );
    chain->SetBranchAddress("n_lep", &n_lep );
    chain->SetBranchAddress("n_ele", &n_ele );
    chain->SetBranchAddress("n_mu", &n_mu );
    chain->SetBranchAddress("ele_pt", &ele_pt );
    chain->SetBranchAddress("ele_eta", &ele_eta );
    chain->SetBranchAddress("ele_phi", &ele_phi );
    chain->SetBranchAddress("mu_pt", &mu_pt );
    chain->SetBranchAddress("mu_eta", &mu_eta );
    chain->SetBranchAddress("mu_phi", &mu_phi );
    chain->SetBranchAddress("jet_pt", &jet_pt );
    chain->SetBranchAddress("n_jets", &n_jets );
    chain->SetBranchAddress("n_btags", &n_btags );



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
		
		chain->SetBranchAddress("nEvent", &nEvent );
		chain->SetBranchAddress("n_lep", &n_lep );
		chain->SetBranchAddress("n_ele", &n_ele );
		chain->SetBranchAddress("n_mu", &n_mu );	
		chain->SetBranchAddress("ele_pt", &ele_pt );
		chain->SetBranchAddress("ele_eta", &ele_eta );
		chain->SetBranchAddress("ele_phi", &ele_phi );
		chain->SetBranchAddress("mu_pt", &mu_pt );
        chain->SetBranchAddress("mu_eta", &mu_eta );
        chain->SetBranchAddress("mu_phi", &mu_phi );
        chain->SetBranchAddress("jet_pt", &jet_pt );
		chain->SetBranchAddress("n_jets", &n_jets );
        chain->SetBranchAddress("n_btags", &n_btags );

		// Single Lepton events
		if (n_lep == 1) {
			++N_SL;	
			N_jets->Fill(n_jets);
			N_btags->Fill(n_btags);
			if (n_ele == 1) {
				++N_e;
				pt_e1->Fill(ele_pt[0]);
				eta_e1->Fill(ele_eta[0]);
				phi_e1 ->Fill(ele_phi[0]);
			}
			else if (n_mu == 1) {
				++N_mu;
				pt_mu1->Fill(mu_pt[0]);
				eta_mu1->Fill(mu_eta[0]);
				phi_mu1 ->Fill(mu_phi[0]);
			}
		}

        else if (n_lep == 2) {
            if (n_ele == 2)
                ++N_ee;
            else if (n_mu == 2)
                ++N_mu;
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


