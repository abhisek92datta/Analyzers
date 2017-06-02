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

    double sum_gen_weights = 0;

    // input filename
	std::string treefilename = "ttHbbNtuple.root";
	
	TChain *chain = new TChain("ttHbb/eventTree");

    chain->Add(treefilename.c_str());

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
                njets_di_emu->Fill(n_jets, tot_weight);
                nbtags_di_emu->Fill(n_btags, tot_weight);
                npv_di_emu->Fill(npv, tot_weight);
            }
        }
		
	}   

    N_sel = N_SL + N_DL;

    std::cout << " Done! " << std::endl;
 	std::cout<<"**********************************************************************************************\n";
    std::cout<<"Total No. of events : "<<N_total<<"\n";   
    std::cout<<"No. of Single Electron Events: "<<N_e<<"\n";
    std::cout<<"No. of Single Muon Events: "<<N_mu<<"\n";
    std::cout<<"No. of Double Electron Events: "<<N_ee<<"\n";
    std::cout<<"No. of Electron Muon Events: "<<N_emu<<"\n";
    std::cout<<"No. of Double Muon Events: "<<N_mumu<<"\n";
    std::cout<<"Total No. of Single Lepton Events: "<<N_SL<<"\n";
    std::cout<<"Total No. of Dilepton Events: "<<N_DL<<"\n";
    std::cout<<"Total No. of Selected Events: "<<N_sel<<"\n";
    std::cout<<"Sum of Generator weights of all Events: "<<sum_gen_weights<<"\n";
    std::cout<<"**********************************************************************************************\n";
    histofile.Write();
    histofile.Close();
    
    return;
}


