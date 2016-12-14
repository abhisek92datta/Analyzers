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

//#include "Analyzers/ttH_analyzer/interface/CU_ttH_EDA_Ntuple.h"

void Control_Plots_tool( int maxNentries=-1, int Njobs=1, int jobN=1 ) {
	
	int N_total = 0;
	int N_SE = 0;
	int N_SMu = 0;
	int N_SL = 0;
	std::string treefilename = "ttHbbNtuple_tth.root";	
	
	std::string histofilename = "Control_Plots_tth.root";
	
	TChain *chain = new TChain("ttHbb/eventTree");
    chain->Add(treefilename.c_str());

	//CU_ttH_EDA_Ntuple *hbbNtuple = 0;
    //chain->SetBranchAddress("hbbNtuple", &hbbNtuple );

    TFile histofile(histofilename.c_str(),"recreate");
    histofile.cd();
    
	TH1::SetDefaultSumw2();

	TH1D* pt_e0 = new TH1D("pt_e0","Leading Electron pT Distribution ;pT (GeV);Nr. of Events",150,0,300);
	TH1D* eta_e0 = new TH1D("eta_e0","Leading Electron eta Distribution ;#eta; Nr. of Events",30,-3,3);
	TH1D* phi_e0 = new TH1D("phi_e0","Leading Electron phi Distribution ;#phi; Nr. of Events",30,-3,3);
	TH1D* pt_mu0 = new TH1D("pt_mu0","Leading Muon pT Distribution ;pT (GeV);Nr. of Events",150,0,300);
	TH1D* eta_mu0 = new TH1D("eta_mu0","Leading Muon eta Distribution ;#eta; Nr. of Events",30,-3,3);
	TH1D* phi_mu0 = new TH1D("phi_mu0","Leading Muon phi Distribution ;#phi; Nr. of Events",30,-3,3);
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
	double ele0_pt,ele0_eta,ele0_phi,mu0_pt,mu0_eta,mu0_phi;
	int n_jets, n_btags;

	// to ensure proper reading from root file in the event loop	
	chain->GetEntry(0);
        chain->SetBranchAddress("nEvent", &nEvent );
        chain->SetBranchAddress("n_lep", &n_lep );
        chain->SetBranchAddress("n_ele", &n_ele );
        chain->SetBranchAddress("n_mu", &n_mu );
        chain->SetBranchAddress("ele0_pt", &ele0_pt );
        chain->SetBranchAddress("ele0_eta", &ele0_eta );
        chain->SetBranchAddress("ele0_phi", &ele0_phi );
        chain->SetBranchAddress("mu0_pt", &mu0_pt );
        chain->SetBranchAddress("mu0_eta", &mu0_eta );
        chain->SetBranchAddress("mu0_phi", &mu0_phi );
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
		chain->SetBranchAddress("ele0_pt", &ele0_pt );
		chain->SetBranchAddress("ele0_eta", &ele0_eta );
		chain->SetBranchAddress("ele0_phi", &ele0_phi );
		chain->SetBranchAddress("mu0_pt", &mu0_pt );
        chain->SetBranchAddress("mu0_eta", &mu0_eta );
        chain->SetBranchAddress("mu0_phi", &mu0_phi );
		chain->SetBranchAddress("n_jets", &n_jets );
        chain->SetBranchAddress("n_btags", &n_btags );

		// Single Lepton events
		if (n_lep == 1) {
			++N_SL;	
			N_jets->Fill(n_jets);
			N_btags->Fill(n_btags);
			if (n_ele == 1) {
				++N_SE;
				pt_e0->Fill(ele0_pt);
				eta_e0->Fill(ele0_eta);
				phi_e0 ->Fill(ele0_phi);
			}
			else if (n_mu == 1) {
				++N_SMu;
				pt_mu0->Fill(mu0_pt);
				eta_mu0->Fill(mu0_eta);
				phi_mu0 ->Fill(mu0_phi);
			}
		}
		
		
	}   
       
    std::cout << " Done! " << std::endl;
 	std::cout<<"**********************************************************************************************\n";
    std::cout<<"Total No. of events : "<<N_total<<"\n";   
    std::cout<<"No. of Single Electron Events: "<<N_SE<<"\n";
    std::cout<<"No. of Single Muon Events: "<<N_SMu<<"\n";
    std::cout<<"Total No. of Single Lepton Events: "<<N_SL<<"\n";
    std::cout<<"**********************************************************************************************\n";
    histofile.Write();
    histofile.Close();
    
    return;
}


