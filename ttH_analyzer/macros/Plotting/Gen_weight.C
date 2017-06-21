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

void Gen_weight( int maxNentries=-1, int Njobs=1, int jobN=1 ) {

    int N_total = 0;
    int N_sel = 0;

    double sum_gen_weights = 0;

    // input filename
    //std::string treefilename = "/eos/uscms/store/user/adatta/ttH_Analysis/v1/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/ttH_Analyzer_ttjets_pp/170603_213946/0000/ttHbbNtuple*.root";
    std::string treefilename = "/eos/uscms/store/user/adatta/ttH_Analysis/v1/ttbb_4FS_OpenLoops_13TeV-sherpa/ttH_Analyzer_ttjets_ol/170606_013446/0000/ttHbbNtuple*.root";
	
	TChain *chain = new TChain("ttHbb/eventTree");

    chain->Add(treefilename.c_str());

	int nentries = chain->GetEntries();
	int NeventsPerJob = int( double(nentries)/double(Njobs) + 0.000001 ) + 1;

  	int firstEvent = (jobN-1)*NeventsPerJob + 1;
  	int lastEvent  = firstEvent + NeventsPerJob;
 	if( jobN==Njobs ) lastEvent = -1;
  	if( jobN==1 ) firstEvent = 0;
 
    double gen_weight;
    int ttHFGenFilter;
    int n_lep;

    chain->SetBranchAddress("n_lep", &n_lep );
    chain->SetBranchAddress("gen_weight", &gen_weight );
    chain->SetBranchAddress("ttHFGenFilter", &ttHFGenFilter );

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

        sum_gen_weights = sum_gen_weights + gen_weight;

        if(ttHFGenFilter!=1)
            continue;

        if(n_lep==1 || n_lep==2)
            N_sel++;


	}   

    std::cout << " Done! " << std::endl;
 	std::cout<<"**********************************************************************************************\n";
    std::cout<<"Total No. of events : "<<N_total<<"\n";
    std::cout<<"Total No. of selected events : "<<N_sel<<"\n";
    std::cout<<"Sum of Generator weights of all Events: "<<sum_gen_weights<<"\n";
    std::cout<<"**********************************************************************************************\n";
    
    return;
}


