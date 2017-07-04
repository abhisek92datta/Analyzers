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
    int N_1_check = 0;
    int N_2_check = 0;

    int N_ttHFGenFilter = 0;
    int N_ttbb = 0;
    int N_ttb = 0;
    int N_tt2b = 0;
    int N_ttcc = 0;
    int N_ttlf = 0;
    int N_all_cat = 0;
    int N_1xx = 0;
    int N_2xx = 0;


    double per_ttHFGenFilter = 0;
    double per_ttbb = 0;
    double per_ttb = 0;
    double per_tt2b = 0;
    double per_ttcc = 0;
    double per_ttlf = 0;

    double sum_gen_weights = 0;

    // input filename
    //std::string treefilename = "ttHbbNtuple.root";

    std::string treefilename = "/eos/uscms/store/user/adatta/ttH_Analysis/v1/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/ttH_Analyzer_ttjets_pp/170603_213946/0000/ttHbbNtuple*.root";
    //std::string treefilename = "/eos/uscms/store/user/adatta/ttH_Analysis/v1/ttbb_4FS_OpenLoops_13TeV-sherpa/ttH_Analyzer_ttjets_ol/170606_013446/0000/ttHbbNtuple*.root";

    //std::string treefilename = "/eos/uscms/store/user/adatta/ttH_Analysis/v1/TTToSemilepton_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/ttH_Analyzer_ttjets_pp_sl/170603_214133/0000/ttHbbNtuple*.root";
    //std::string treefilename = "/eos/uscms/store/user/adatta/ttH_Analysis/v1/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/ttH_Analyzer_ttjets_pp_dl/170603_214359/0000/ttHbbNtuple*.root";
    //std::string treefilename = "/eos/uscms/store/user/adatta/ttH_Analysis/v1/TTToSemilepton_ttbbFilter_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/ttH_Analyzer_ttjets_pp_sl_hf/170603_214816/0000/ttHbbNtuple*.root";
    //std::string treefilename = "/eos/uscms/store/user/adatta/ttH_Analysis/v1/TTTo2L2Nu_ttbbFilter_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/ttH_Analyzer_ttjets_pp_dl_hf/170603_214912/0000/ttHbbNtuple*.root";


	
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
    int ttHf_cat;
    int n_lep, n_btags;

    int mod_ttHFcat;

    chain->SetBranchAddress("n_lep", &n_lep );
    chain->SetBranchAddress("gen_weight", &gen_weight );
    chain->SetBranchAddress("ttHFGenFilter", &ttHFGenFilter );
    chain->SetBranchAddress("ttHf_cat", &ttHf_cat );
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

        sum_gen_weights = sum_gen_weights + gen_weight;

        if(ttHFGenFilter!=-9999)
            N_1_check++;
        if(ttHf_cat!=-9999)
            N_2_check++;

        if(n_lep != 1 && n_lep != 2)
            continue;

        //if(n_btags<4)
        //    continue;

        N_sel++;

        if(ttHFGenFilter == 1)
            N_ttHFGenFilter++;

        //if(ttHFGenFilter != 1)
        //    continue;
        //N_ttHFGenFilter++;

        mod_ttHFcat = ttHf_cat%100;
        //mod_ttHFcat = ttHf_cat;

        if(ttHf_cat>=100 && ttHf_cat<200 )
            N_1xx++;
        else if(ttHf_cat >= 200 && ttHf_cat < 300)
            N_2xx++;

        if( (mod_ttHFcat >= 53) && (mod_ttHFcat <= 55) ){
            N_ttbb++;
            //std::cout<<"ttbb : "<<ttHf_cat<<"  "<<mod_ttHFcat<<"\n";
        }

        if( mod_ttHFcat == 51 ){
            N_ttb++;
            //std::cout<<"ttb : "<<ttHf_cat<<"  "<<mod_ttHFcat<<"\n";
        }

        if( mod_ttHFcat == 52 ){
            N_tt2b++;
            //std::cout<<"tt2b : "<<ttHf_cat<<"  "<<mod_ttHFcat<<"\n";
        }

        if( (mod_ttHFcat >= 41) && (mod_ttHFcat <= 45) ){
            N_ttcc++;
            //std::cout<<"ttcc : "<<ttHf_cat<<"  "<<mod_ttHFcat<<"\n";
        }

        if( mod_ttHFcat == 0 ){
            N_ttlf++;
            //std::cout<<"ttlf : "<<ttHf_cat<<"  "<<mod_ttHFcat<<"\n\n";
        }


	}

    N_all_cat = N_ttbb + N_ttb + N_tt2b + N_ttcc + N_ttlf;

    per_ttHFGenFilter = (double(N_ttHFGenFilter)/N_sel)*100;
    per_ttbb = (double(N_ttbb)/N_sel)*100;
    per_ttb = (double(N_ttb)/N_sel)*100;
    per_tt2b = (double(N_tt2b)/N_sel)*100;
    per_ttcc = (double(N_ttcc)/N_sel)*100;
    per_ttlf = (double(N_ttlf)/N_sel)*100;

    /*
    per_ttHFGenFilter = (double(N_ttHFGenFilter)/N_sel)*100;
    per_ttbb = (double(N_ttbb)/N_ttHFGenFilter)*100;
    per_ttb = (double(N_ttb)/N_ttHFGenFilter)*100;
    per_tt2b = (double(N_tt2b)/N_ttHFGenFilter)*100;
    per_ttcc = (double(N_ttcc)/N_ttHFGenFilter)*100;
    per_ttlf = (double(N_ttlf)/N_ttHFGenFilter)*100;
     */


    std::cout << " Done! " << std::endl;
 	std::cout<<"**********************************************************************************************\n";
    std::cout<<"Total No. of events : "<<N_total<<"\n";
    std::cout<<"Total No. of selected events : "<<N_sel<<"\n";
    std::cout<<"ttHFGenFilter check : "<<N_1_check<<"\n";
    std::cout<<"ttHFcat check : "<<N_2_check<<"\n\n";
    std::cout<<"Ratio of N_1xx / N_2xx : "<<(double)N_1xx/N_2xx<<"\n";
    std::cout<<"% of selected events passing ttHFGenFilter : "<<per_ttHFGenFilter<<"\n";
    std::cout<<"% of selected events which are ttbb : "<<per_ttbb<<"\n";
    std::cout<<"% of selected events which are ttb : "<<per_ttb<<"\n";
    std::cout<<"% of selected events which are tt2b : "<<per_tt2b<<"\n";
    std::cout<<"% of selected events which are ttcc : "<<per_ttcc<<"\n";
    std::cout<<"% of selected events which are ttlf : "<<per_ttlf<<"\n";
    std::cout<<"Total No. of selected events in all categories : "<<N_all_cat<<"\n";
    std::cout<<"Sum of Generator weights of all Events: "<<sum_gen_weights<<"\n";
    std::cout<<"**********************************************************************************************\n";
    
    return;
}


