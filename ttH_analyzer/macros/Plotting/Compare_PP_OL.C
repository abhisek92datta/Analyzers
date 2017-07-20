#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH1D.h"
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
#include "TEfficiency.h"
#include "THStack.h"


void Compare_PP_OL() {

    TH1::SetDefaultSumw2();

    // Set Integrated Luminosity

    double L = 35.9182; // 1/fb

    // List of Histograms

    ifstream fin;
    fin.open("dist_histo_names.txt");
    char histonames[200][200];
    int nhistos = 0;

    while(!fin.eof()){
        fin>>histonames[nhistos];
        nhistos++;
    }
    fin.close();

    // Output Parameters
    std::vector<std::string> out_file_e;
    std::vector<std::string> out_file_mu;
    std::vector<std::string> out_file_ee;
    std::vector<std::string> out_file_emu;
    std::vector<std::string> out_file_mumu;
    std::vector<std::string> out_file_e_2D;
    std::vector<std::string> out_file_mu_2D;
    std::vector<std::string> out_file_ee_2D;
    std::vector<std::string> out_file_emu_2D;
    std::vector<std::string> out_file_mumu_2D;


    ///////////////////////////////////////////////////////////////////////
    // Reading the Histograms from the different files
    ///////////////////////////////////////////////////////////////////////

    TFile *f1 = new TFile("Distribution_cut20_pp_noGenHF_2.root");
    TFile *f2 = new TFile("Distribution_cut20_ol_noGenHF_2.root");
    //TFile *f3 = new TFile("Distribution_data.root");

    std::vector<TH1D*> ttjets_pp_e;
    std::vector<TH1D*> ttjets_pp_mu;
    std::vector<TH1D*> ttjets_pp_ee;
    std::vector<TH1D*> ttjets_pp_emu;
    std::vector<TH1D*> ttjets_pp_mumu;
    std::vector<TH2D*> ttjets_pp_e_2D;
    std::vector<TH2D*> ttjets_pp_mu_2D;
    std::vector<TH2D*> ttjets_pp_ee_2D;
    std::vector<TH2D*> ttjets_pp_emu_2D;
    std::vector<TH2D*> ttjets_pp_mumu_2D;

    std::vector<TH1D*> ttjets_pp_e_norm;
    std::vector<TH1D*> ttjets_pp_mu_norm;
    std::vector<TH1D*> ttjets_pp_ee_norm;
    std::vector<TH1D*> ttjets_pp_emu_norm;
    std::vector<TH1D*> ttjets_pp_mumu_norm;

    std::vector<TH1D*> ttjets_ol_e;
    std::vector<TH1D*> ttjets_ol_mu;
    std::vector<TH1D*> ttjets_ol_ee;
    std::vector<TH1D*> ttjets_ol_emu;
    std::vector<TH1D*> ttjets_ol_mumu;
    std::vector<TH2D*> ttjets_ol_e_2D;
    std::vector<TH2D*> ttjets_ol_mu_2D;
    std::vector<TH2D*> ttjets_ol_ee_2D;
    std::vector<TH2D*> ttjets_ol_emu_2D;
    std::vector<TH2D*> ttjets_ol_mumu_2D;

    std::vector<TH1D*> ttjets_ol_e_norm;
    std::vector<TH1D*> ttjets_ol_mu_norm;
    std::vector<TH1D*> ttjets_ol_ee_norm;
    std::vector<TH1D*> ttjets_ol_emu_norm;
    std::vector<TH1D*> ttjets_ol_mumu_norm;

    /*
    std::vector<TH1D*> data_e;
    std::vector<TH1D*> data_mu;
    std::vector<TH1D*> data_ee;
    std::vector<TH1D*> data_emu;
    std::vector<TH1D*> data_mumu;
    */


    int n_e, n_mu, n_ee, n_emu, n_mumu, n_e_2D, n_mu_2D, n_ee_2D, n_emu_2D, n_mumu_2D;
    n_e = n_mu = n_ee = n_emu = n_mumu = n_e_2D = n_mu_2D = n_ee_2D = n_emu_2D = n_mumu_2D = 0;

    for(int i=0; i<nhistos; i++){

        std::string histo = histonames[i];
        TH1D *h1, *h11;
        TH2D *h1_2D;
        TH1D *h2, *h22;
        TH2D *h2_2D;
        //TH1D *h3;
        //TH2D *h3_2D;

        // ttjets_pp
        if(histo.find("2D") != std::string::npos){
            h1_2D = (TH2D*)f1->Get((histo).c_str());
            h1_2D->SetName(("ttjets_pp_" + histo).c_str());
            h1_2D->SetTitle(("MC ttjets_pp : " + histo).c_str());
        }
        else {
            h1 = (TH1D*)f1->Get((histo).c_str());
            h1->SetName(("ttjets_pp_" + histo).c_str());
            h1->SetTitle(("MC ttjets_pp : " + histo).c_str());
			h11 = (TH1D*)f1->Get((histo).c_str());
			h11->SetName(("norm_ttjets_pp_" + histo).c_str());
			h11->SetTitle(("MC norm_ttjets_pp : " + histo).c_str());
        }

        // ttjets_ol
        if(histo.find("2D") != std::string::npos){
            h2_2D = (TH2D*)f2->Get((histo).c_str());
            h2_2D->SetName(("ttjets_ol_" + histo).c_str());
            h2_2D->SetTitle(("MC ttjets_ol : " + histo).c_str());
        }
        else {
            h2 = (TH1D*)f2->Get((histo).c_str());
            h2->SetName(("ttjets_ol_" + histo).c_str());
            h2->SetTitle(("MC ttjets_ol : " + histo).c_str());
			h22 = (TH1D*)f2->Get((histo).c_str());
			h22->SetName(("norm_ttjets_ol_" + histo).c_str());
			h22->SetTitle(("MC norm_ttjets_ol : " + histo).c_str());

        }

        /*
        // data
        if(histo.find("2D") != std::string::npos){
            h3_2D = (TH2D*)f3->Get((histo).c_str());
            h3_2D->SetName(("data_" + histo).c_str());
            h3_2D->SetTitle(("DATA : " + histo).c_str());
        }
        else {
            h3 = (TH1D*)f3->Get((histo).c_str());
            h3->SetName(("data_" + histo).c_str());
            h3->SetTitle(("DATA : " + histo).c_str());
        }
        */

        if(histo.find("_sl_e") != std::string::npos){
            if(histo.find("2D") != std::string::npos){
                ttjets_pp_e_2D.push_back(h1_2D);
                ttjets_ol_e_2D.push_back(h2_2D);
                //data_e_2D.push_back(h3_2D);

                out_file_e_2D.push_back((histo + ".png").c_str());
                n_e_2D++;
            }
            else{
                ttjets_pp_e.push_back(h1);
                ttjets_ol_e.push_back(h2);
				ttjets_pp_e_norm.push_back(h11);
				ttjets_ol_e_norm.push_back(h22);
                //data_e.push_back(h3);

                out_file_e.push_back((histo + ".png").c_str());
                n_e++;
            }
        }

        else if(histo.find("_sl_mu") != std::string::npos){
            if(histo.find("2D") != std::string::npos){
                ttjets_pp_mu_2D.push_back(h1_2D);
                ttjets_ol_mu_2D.push_back(h2_2D);
                //data_mu_2D.push_back(h3_2D);

                out_file_mu_2D.push_back((histo + ".png").c_str());
                n_mu_2D++;
            }
            else{
                ttjets_pp_mu.push_back(h1);
                ttjets_ol_mu.push_back(h2);
				ttjets_pp_mu_norm.push_back(h11);
				ttjets_ol_mu_norm.push_back(h22);
                //data_mu.push_back(h3);

                out_file_mu.push_back((histo + ".png").c_str());
                n_mu++;
            }
        }

        else if(histo.find("_di_ee") != std::string::npos){
            if(histo.find("2D") != std::string::npos){
                ttjets_pp_ee_2D.push_back(h1_2D);
                ttjets_ol_ee_2D.push_back(h2_2D);
                //data_ee_2D.push_back(h3_2D);

                out_file_ee_2D.push_back((histo + ".png").c_str());
                n_ee_2D++;
            }
            else{
                ttjets_pp_ee.push_back(h1);
                ttjets_ol_ee.push_back(h2);
				ttjets_pp_ee_norm.push_back(h11);
				ttjets_ol_ee_norm.push_back(h22);
                //data_ee.push_back(h3);

                out_file_ee.push_back((histo + ".png").c_str());
                n_ee++;
            }
        }

        else if(histo.find("_di_emu") != std::string::npos){
            if(histo.find("2D") != std::string::npos){
                ttjets_pp_emu_2D.push_back(h1_2D);
                ttjets_ol_emu_2D.push_back(h2_2D);
                //data_emu_2D.push_back(h3_2D);

                out_file_emu_2D.push_back((histo + ".png").c_str());
                n_emu_2D++;
            }
            else{
                ttjets_pp_emu.push_back(h1);
                ttjets_ol_emu.push_back(h2);
				ttjets_pp_emu_norm.push_back(h11);
				ttjets_ol_emu_norm.push_back(h22);
                //data_emu.push_back(h3);

                out_file_emu.push_back((histo + ".png").c_str());
                n_emu++;
            }
        }

        else if(histo.find("_di_mumu") != std::string::npos){
            if(histo.find("2D") != std::string::npos){
                ttjets_pp_mumu_2D.push_back(h1_2D);
                ttjets_ol_mumu_2D.push_back(h2_2D);
                //data_mumu_2D.push_back(h3_2D);

                out_file_mumu_2D.push_back((histo + ".png").c_str());
                n_mumu_2D++;
            }
            else{
                ttjets_pp_mumu.push_back(h1);
                ttjets_ol_mumu.push_back(h2);
				ttjets_pp_mumu_norm.push_back(h11);
				ttjets_ol_mumu_norm.push_back(h22);
                //data_mumu.push_back(h3);

                out_file_mumu.push_back((histo + ".png").c_str());
                n_mumu++;
            }
        }

    }




    ///////////////////////////////////////////////////////////////////////
    // Scaling and Normlization of MC signal and background
    ///////////////////////////////////////////////////////////////////////

    //ttjets_pp

    double factor_ttjets_pp = 1.0;
	double sigma_ttjets_pp = 730.6; // pb
    double N_total_ttjets_pp = 77078997;
    double sum_gen_weight_ttjets_pp = 7.7079e+07;
    double norm_ttjets_pp = (L*sigma_ttjets_pp*1000)/(sum_gen_weight_ttjets_pp);
    double scale_ttjets_pp = factor_ttjets_pp*norm_ttjets_pp;

    //ttjets_ol
    double factor_ttjets_ol = 1.0;
    double sigma_ttjets_ol = 25.9; // pb
    double N_total_ttjets_ol = 19659611;
    double sum_gen_weight_ttjets_ol = 1.76729e+12;
    double norm_ttjets_ol = (L*sigma_ttjets_ol*1000)/(sum_gen_weight_ttjets_ol);
    double scale_ttjets_ol = factor_ttjets_ol*norm_ttjets_ol;

    double integral;

	for(int i=0; i<n_e; i++){
        ttjets_pp_e[i]->Scale(scale_ttjets_pp);
        integral = ttjets_pp_e[i]->Integral();
		ttjets_pp_e_norm[i]->Scale(scale_ttjets_pp/integral);
        ttjets_ol_e[i]->Scale(scale_ttjets_ol);
        integral = ttjets_ol_e[i]->Integral();
		ttjets_ol_e_norm[i]->Scale(scale_ttjets_ol/integral);
    }

    for(int i=0; i<n_e_2D; i++){
        ttjets_pp_e_2D[i]->Scale(scale_ttjets_pp);
        ttjets_ol_e_2D[i]->Scale(scale_ttjets_ol);
    }
    for(int i=0; i<n_mu; i++){
        ttjets_pp_mu[i]->Scale(scale_ttjets_pp);
        integral = ttjets_pp_mu[i]->Integral();
        ttjets_pp_mu_norm[i]->Scale(scale_ttjets_pp/integral);
        ttjets_ol_mu[i]->Scale(scale_ttjets_ol);
        integral = ttjets_ol_mu[i]->Integral();
        ttjets_ol_mu_norm[i]->Scale(scale_ttjets_ol/integral);
    }
    for(int i=0; i<n_mu_2D; i++){
        ttjets_pp_mu_2D[i]->Scale(scale_ttjets_pp);
        ttjets_ol_mu_2D[i]->Scale(scale_ttjets_ol);
    }
    for(int i=0; i<n_ee; i++){
        ttjets_pp_ee[i]->Scale(scale_ttjets_pp);
        integral = ttjets_pp_ee[i]->Integral();
        ttjets_pp_ee_norm[i]->Scale(scale_ttjets_pp/integral);
        ttjets_ol_ee[i]->Scale(scale_ttjets_ol);
        integral = ttjets_ol_ee[i]->Integral();
        ttjets_ol_ee_norm[i]->Scale(scale_ttjets_ol/integral);
    }
    for(int i=0; i<n_ee_2D; i++){
        ttjets_pp_ee_2D[i]->Scale(scale_ttjets_pp);
        ttjets_ol_ee_2D[i]->Scale(scale_ttjets_ol);
    }
    for(int i=0; i<n_emu; i++){
        ttjets_pp_emu[i]->Scale(scale_ttjets_pp);
        integral = ttjets_pp_emu[i]->Integral();
        ttjets_pp_emu_norm[i]->Scale(scale_ttjets_pp/integral);
        ttjets_ol_emu[i]->Scale(scale_ttjets_ol);
        integral = ttjets_ol_emu[i]->Integral();
        ttjets_ol_emu_norm[i]->Scale(scale_ttjets_ol/integral);
    }
    for(int i=0; i<n_emu_2D; i++){
        ttjets_pp_emu_2D[i]->Scale(scale_ttjets_pp);
        ttjets_ol_emu_2D[i]->Scale(scale_ttjets_ol);
    }
    for(int i=0; i<n_mumu; i++){
        ttjets_pp_mumu[i]->Scale(scale_ttjets_pp);
        integral = ttjets_pp_mumu[i]->Integral();
        ttjets_pp_mumu_norm[i]->Scale(scale_ttjets_pp/integral);
        ttjets_ol_mumu[i]->Scale(scale_ttjets_ol);
        integral = ttjets_ol_mumu[i]->Integral();
        ttjets_ol_mumu_norm[i]->Scale(scale_ttjets_ol/integral);
    }
    for(int i=0; i<n_mumu_2D; i++){
        ttjets_pp_mumu_2D[i]->Scale(scale_ttjets_pp);
        ttjets_ol_mumu_2D[i]->Scale(scale_ttjets_ol);
    }


    ///////////////////////////////////////////////////////////////////////
    // Plotting
    ///////////////////////////////////////////////////////////////////////

    // setting plot parameters
    std::vector<double> x_l_sl, x_u_sl, y_l_sl, y_u_sl, y_n_l_sl, y_n_u_sl;
    std::vector<double> x_l_di, x_u_di, y_l_di, y_u_di, y_n_l_di, y_n_u_di;
    std::vector<double> x_l_sl_2D, x_u_sl_2D, y_l_sl_2D, y_u_sl_2D;
    std::vector<double> x_l_di_2D, x_u_di_2D, y_l_di_2D, y_u_di_2D;
    std::vector<std::string> canvas_e;
    std::vector<std::string> canvas_mu;
    std::vector<std::string> canvas_ee;
    std::vector<std::string> canvas_emu;
    std::vector<std::string> canvas_mumu;
    std::vector<std::string> canvas_e_ratio;
    std::vector<std::string> canvas_mu_ratio;
    std::vector<std::string> canvas_ee_ratio;
    std::vector<std::string> canvas_emu_ratio;
    std::vector<std::string> canvas_mumu_ratio;
    std::vector<std::string> canvas_e_2D;
    std::vector<std::string> canvas_mu_2D;
    std::vector<std::string> canvas_ee_2D;
    std::vector<std::string> canvas_emu_2D;
    std::vector<std::string> canvas_mumu_2D;

    canvas_e.push_back("Electron pT (for Single Electron Channel) ; pT (Electron) [GeV] ; Nr. of Events");
    canvas_e.push_back("Electron eta (for Single Electron Channel) ; eta (Electron) ; Nr. of Events");
    canvas_e.push_back("Electron phi (for Single Electron Channel) ; phi (Electron) ; Nr. of Events");
    canvas_e.push_back("Leading Jet pT (for Single Electron Channel) ; pT (Leading Jet) [GeV] ; Nr. of Events");
    canvas_e.push_back("2nd Leading Jet pT (for Single Electron Channel) ; pT (2nd Leading Jet) [GeV] ; Nr. of Events");
    canvas_e.push_back("3rd Leading Jet pT (for Single Electron Channel) ; pT (3rd Leading Jet) [GeV] ; Nr. of Events");
    canvas_e.push_back("4th Leading Jet pT (for Single Electron Channel) ; pT (4th Leading Jet) [GeV] ; Nr. of Events");
    canvas_e.push_back("5th Leading Jet pT (for Single Electron Channel) ; pT (5th Leading Jet) [GeV] ; Nr. of Events");
    canvas_e.push_back("6th Leading Jet pT (for Single Electron Channel) ; pT (6th Leading Jet) [GeV] ; Nr. of Events");
    canvas_e.push_back("Jet pT (for Single Electron Channel) ; pT (Jet) [GeV] ; Nr. of Events");
    canvas_e.push_back("Jet CSV (for Single Electron Channel) ; CSV (Jet) ; Nr. of Events");
    canvas_e.push_back("HT (for Single Electron Channel) ; HT [GeV] ; Nr. of Events");
    canvas_e.push_back("Nr. of Jets (for Single Electron Channel) ; Nr. of jets ; Nr. of Events");
    canvas_e.push_back("Nr. of B-tags (for Single Electron Channel) ; Nr. of b-tags ; Nr. of Events");
    canvas_e.push_back("Nr. of Primary Vertices (for Single Electron Channel) ; Nr. of primary vertices ; Nr. of Events");
    canvas_e.push_back("Leading Gen Electron pT (for Single Electron Channel) ; pT (Leading Gen Electron) [GeV] ; Nr. of Events");
    canvas_e.push_back("Leading Gen Electron eta (for Single Electron Channel) ; eta (Leading Gen Electron) ; Nr. of Events");
    canvas_e.push_back("Leading Gen Jet pT (for Single Electron Channel) ; pT (Leading Gen Jet) [GeV] ; Nr. of Events");
    canvas_e.push_back("2nd Leading Gen Jet pT (for Single Electron Channel) ; pT (2nd Leading Gen Jet) [GeV] ; Nr. of Events");
    canvas_e.push_back("Gen HT (for Single Electron Channel) ; Gen HT [GeV] ; Nr. of Events");
    canvas_e.push_back("Nr. of Gen Jets(for Single Electron Channel) ; Nr. of genjets ; Nr. of Events");
    canvas_e_ratio.push_back("; pT (Electron) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_e_ratio.push_back("; eta (Electron) ; Powheg-Pythia/Sherpa+OL");
    canvas_e_ratio.push_back("; phi (Electron) ; Powheg-Pythia/Sherpa+OL");
    canvas_e_ratio.push_back("; pT (Leading Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_e_ratio.push_back("; pT (2nd Leading Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_e_ratio.push_back("; pT (3rd Leading Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_e_ratio.push_back("; pT (4th Leading Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_e_ratio.push_back("; pT (5th Leading Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_e_ratio.push_back("; pT (6th Leading Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_e_ratio.push_back("; pT (Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_e_ratio.push_back("; CSV (Jet) ; Powheg-Pythia/Sherpa+OL");
    canvas_e_ratio.push_back("; HT [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_e_ratio.push_back("; Nr. of jets ; Powheg-Pythia/Sherpa+OL");
    canvas_e_ratio.push_back("; Nr. of b-tags ; Powheg-Pythia/Sherpa+OL");
    canvas_e_ratio.push_back("; Nr. of primary vertices ; Powheg-Pythia/Sherpa+OL");
    canvas_e_ratio.push_back("; pT (Leading Gen Electron) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_e_ratio.push_back("; eta (Leading Gen Electron) ; Powheg-Pythia/Sherpa+OL");
    canvas_e_ratio.push_back("; pT (Leading Gen Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_e_ratio.push_back("; pT (2nd Leading Gen Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_e_ratio.push_back("; HT [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_e_ratio.push_back("; Nr. of genjets ; Powheg-Pythia/Sherpa+OL");
    canvas_e_2D.push_back("2D : Njets & HT (for Single Electron Channel) ; Nr. of jets ; HT [GeV]");

    canvas_mu.push_back("Muon pT (for Single Muon Channel) ; pT (Muon) [GeV] ; Nr. of Events");
    canvas_mu.push_back("Muon eta (for Single Muon Channel) ; eta (Muon) ; Nr. of Events");
    canvas_mu.push_back("Muon phi (for Single Muon Channel) ; phi (Muon) ; Nr. of Events");
    canvas_mu.push_back("Leading Jet pT (for Single Muon Channel) ; pT (Leading Jet) [GeV] ; Nr. of Events");
    canvas_mu.push_back("2nd Leading Jet pT (for Single Muon Channel) ; pT (2nd Leading Jet) [GeV] ; Nr. of Events");
    canvas_mu.push_back("3rd Leading Jet pT (for Single Muon Channel) ; pT (3rd Leading Jet) [GeV] ; Nr. of Events");
    canvas_mu.push_back("4th Leading Jet pT (for Single Muon Channel) ; pT (4th Leading Jet) [GeV] ; Nr. of Events");
    canvas_mu.push_back("5th Leading Jet pT (for Single Muon Channel) ; pT (5th Leading Jet) [GeV] ; Nr. of Events");
    canvas_mu.push_back("6th Leading Jet pT (for Single Muon Channel) ; pT (6th Leading Jet) [GeV] ; Nr. of Events");
    canvas_mu.push_back("Jet pT (for Single Muon Channel) ; pT (Jet) [GeV] ; Nr. of Events");
    canvas_mu.push_back("Jet CSV (for Single Muon Channel) ; CSV (Jet) ; Nr. of Events");
    canvas_mu.push_back("HT (for Single Muon Channel) ; HT [GeV] ; Nr. of Events");
    canvas_mu.push_back("Nr. of Jets (for Single Muon Channel) ; Nr. of jets ; Nr. of Events");
    canvas_mu.push_back("Nr. of B-tags (for Single Muon Channel) ; Nr. of b-tags ; Nr. of Events");
    canvas_mu.push_back("Nr. of Primary Vertices (for Single Muon Channel) ; Nr. of primary vertices ; Nr. of Events");
    canvas_mu.push_back("Leading Gen Muon pT (for Single Muon Channel) ; pT (Leading Gen Muon) [GeV] ; Nr. of Events");
    canvas_mu.push_back("Leading Gen Muon eta (for Single Muon Channel) ; eta (Leading Gen Muon) ; Nr. of Events");
    canvas_mu.push_back("Leading Gen Jet pT (for Single Muon Channel) ; pT (Leading Gen Jet) [GeV] ; Nr. of Events");
    canvas_mu.push_back("2nd Leading Gen Jet pT (for Single Muon Channel) ; pT (2nd Leading Gen Jet) [GeV] ; Nr. of Events");
    canvas_mu.push_back("Gen HT (for Single Muon Channel) ; Gen HT [GeV] ; Nr. of Events");
    canvas_mu.push_back("Nr. of Gen Jets(for Single Muon Channel) ; Nr. of genjets ; Nr. of Events");
    canvas_mu_ratio.push_back("; pT (Muon) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_mu_ratio.push_back("; eta (Muon) ; Powheg-Pythia/Sherpa+OL");
    canvas_mu_ratio.push_back("; phi (Muon) ; Powheg-Pythia/Sherpa+OL");
    canvas_mu_ratio.push_back("; pT (Leading Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_mu_ratio.push_back("; pT (2nd Leading Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_mu_ratio.push_back("; pT (3rd Leading Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_mu_ratio.push_back("; pT (4th Leading Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_mu_ratio.push_back("; pT (5th Leading Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_mu_ratio.push_back("; pT (6th Leading Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_mu_ratio.push_back("; pT (Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_mu_ratio.push_back("; CSV (Jet) ; Powheg-Pythia/Sherpa+OL");
    canvas_mu_ratio.push_back("; HT [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_mu_ratio.push_back("; Nr. of jets ; Powheg-Pythia/Sherpa+OL");
    canvas_mu_ratio.push_back("; Nr. of b-tags ; Powheg-Pythia/Sherpa+OL");
    canvas_mu_ratio.push_back("; Nr. of primary vertices ; Powheg-Pythia/Sherpa+OL");
    canvas_mu_ratio.push_back("; pT (Leading Gen Muon) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_mu_ratio.push_back("; eta (Leading Gen Muon) ; Powheg-Pythia/Sherpa+OL");
    canvas_mu_ratio.push_back("; pT (Leading Gen Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_mu_ratio.push_back("; pT (2nd Leading Gen Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_mu_ratio.push_back("; HT [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_mu_ratio.push_back("; Nr. of genjets ; Powheg-Pythia/Sherpa+OL");
    canvas_mu_2D.push_back("2D : Njets & HT (for Single Muon Channel) ; Nr. of jets ; HT [GeV]");

    canvas_ee.push_back("Leading Electron pT (for Double Electron Channel) ; pT (Leading Electron) [GeV] ; Nr. of Events");
    canvas_ee.push_back("Leading Electron eta (for Double Electron Channel) ; eta (Leading Electron) ; Nr. of Events");
    canvas_ee.push_back("Leading Electron phi (for Double Electron Channel) ; phi (Leading Electron) ; Nr. of Events");
    canvas_ee.push_back("Sub-leading Electron pT (for Double Electron Channel) ; pT (Sub-leading Electron) [GeV] ; Nr. of Events");
    canvas_ee.push_back("Sub-leading Electron eta (for Double Electron Channel) ; eta (Sub-leading Electron) ; Nr. of Events");
    canvas_ee.push_back("Sub-leading Electron phi (for Double Electron Channel) ; phi (Sub-leading Electron) ; Nr. of Events");
    canvas_ee.push_back("Leading Jet pT (for Double Electron Channel) ; pT (Leading Jet) [GeV] ; Nr. of Events");
    canvas_ee.push_back("2nd Leading Jet pT (for Double Electron Channel) ; pT (2nd Leading Jet) [GeV] ; Nr. of Events");
    canvas_ee.push_back("3rd Leading Jet pT (for Double Electron Channel) ; pT (3rd Leading Jet) [GeV] ; Nr. of Events");
    canvas_ee.push_back("4th Leading Jet pT (for Double Electron Channel) ; pT (4th Leading Jet) [GeV] ; Nr. of Events");
    canvas_ee.push_back("5th Leading Jet pT (for Double Electron Channel) ; pT (5th Leading Jet) [GeV] ; Nr. of Events");
    canvas_ee.push_back("6th Leading Jet pT (for Double Electron Channel) ; pT (6th Leading Jet) [GeV] ; Nr. of Events");
    canvas_ee.push_back("Jet pT (for Double Electron Channel) ; pT (Jet) [GeV] ; Nr. of Events");
    canvas_ee.push_back("Jet CSV (for Double Electron Channel) ; CSV (Jet) ; Nr. of Events");
    canvas_ee.push_back("HT (for Double Electron Channel) ; HT [GeV] ; Nr. of Events");
    canvas_ee.push_back("Nr. of Jets (for Double Electron Channel) ; Nr. of jets ; Nr. of Events");
    canvas_ee.push_back("Nr. of B-tags (for Double Electron Channel) ; Nr. of b-tags ; Nr. of Events");
    canvas_ee.push_back("Nr. of Primary Vertices (for Double Electron Channel) ; Nr. of primary vertices ; Nr. of Events");
    canvas_ee.push_back("Leading Gen Electron pT (for Double Electron Channel) ; pT (Leading Gen Electron) [GeV] ; Nr. of Events");
    canvas_ee.push_back("Leading Gen Electron eta (for Double Electron Channel) ; eta (Leading Gen Electron) ; Nr. of Events");
    canvas_ee.push_back("Sub-leading Gen Electron pT (for Double Electron Channel) ; pT (Sub-leading Gen Electron) [GeV] ; Nr. of Events");
    canvas_ee.push_back("Sub-leading Gen Electron eta (for Double Electron Channel) ; eta (Sub-leading Gen Electron) ; Nr. of Events");
    canvas_ee.push_back("Leading Gen Jet pT (for Double Electron Channel) ; pT (Leading Gen Jet) [GeV] ; Nr. of Events");
    canvas_ee.push_back("2nd Leading Gen Jet pT (for Double Electron Channel) ; pT (2nd Leading Gen Jet) [GeV] ; Nr. of Events");
    canvas_ee.push_back("Gen HT (for Double Electron Channel) ; Gen HT [GeV] ; Nr. of Events");
    canvas_ee.push_back("Nr. of Gen Jets(for Double Electron Channel) ; Nr. of genjets ; Nr. of Events");
    canvas_ee_ratio.push_back("; pT (Leading Electron) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_ee_ratio.push_back("; eta (Leading Electron) ; Powheg-Pythia/Sherpa+OL");
    canvas_ee_ratio.push_back("; phi (Leading Electron) ; Powheg-Pythia/Sherpa+OL");
    canvas_ee_ratio.push_back("; pT (Sub-Leading Electron) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_ee_ratio.push_back("; eta (Sub-Leading Electron) ; Powheg-Pythia/Sherpa+OL");
    canvas_ee_ratio.push_back("; phi (Sub-Leading Electron) ; Powheg-Pythia/Sherpa+OL");
    canvas_ee_ratio.push_back("; pT (Leading Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_ee_ratio.push_back("; pT (2nd Leading Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_ee_ratio.push_back("; pT (3rd Leading Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_ee_ratio.push_back("; pT (4th Leading Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_ee_ratio.push_back("; pT (5th Leading Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_ee_ratio.push_back("; pT (6th Leading Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_ee_ratio.push_back("; pT (Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_ee_ratio.push_back("; CSV (Jet) ; Powheg-Pythia/Sherpa+OL");
    canvas_ee_ratio.push_back("; HT [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_ee_ratio.push_back("; Nr. of jets ; Powheg-Pythia/Sherpa+OL");
    canvas_ee_ratio.push_back("; Nr. of b-tags ; Powheg-Pythia/Sherpa+OL");
    canvas_ee_ratio.push_back("; Nr. of primary vertices ; Powheg-Pythia/Sherpa+OL");
    canvas_ee_ratio.push_back("; pT (Leading Gen Electron) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_ee_ratio.push_back("; eta (Leading Gen Electron) ; Powheg-Pythia/Sherpa+OL");
    canvas_ee_ratio.push_back("; pT (Sub-leading Gen Electron) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_ee_ratio.push_back("; eta (Sub-leading Gen Electron) ; Powheg-Pythia/Sherpa+OL");
    canvas_ee_ratio.push_back("; pT (Leading Gen Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_ee_ratio.push_back("; pT (2nd Leading Gen Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_ee_ratio.push_back("; HT [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_ee_ratio.push_back("; Nr. of genjets ; Powheg-Pythia/Sherpa+OL");
    canvas_ee_2D.push_back("2D : Njets & HT (for Double Electron Channel) ; Nr. of jets ; HT [GeV]");

    canvas_emu.push_back("Electron pT (for Electron Muon Channel) ; pT (Electron) [GeV] ; Nr. of Events");
    canvas_emu.push_back("Electron eta (for Electron Muon Channel) ; eta (Electron) ; Nr. of Events");
    canvas_emu.push_back("Electron phi (for Electron Muon Channel) ; phi (Electron) ; Nr. of Events");
    canvas_emu.push_back("Muon pT (for Electron Muon Channel) ; pT (Muon) [GeV] ; Nr. of Events");
    canvas_emu.push_back("Muon eta (for Electron Muon Channel) ; eta (Muon) ; Nr. of Events");
    canvas_emu.push_back("Muon phi (for Electron Muon Channel) ; phi (Muon) ; Nr. of Events");
    canvas_emu.push_back("Leading Jet pT (for Electron Muon Channel) ; pT (Leading Jet) [GeV] ; Nr. of Events");
    canvas_emu.push_back("2nd Leading Jet pT (for Electron Muon Channel) ; pT (2nd Leading Jet) [GeV] ; Nr. of Events");
    canvas_emu.push_back("3rd Leading Jet pT (for Electron Muon Channel) ; pT (3rd Leading Jet) [GeV] ; Nr. of Events");
    canvas_emu.push_back("4th Leading Jet pT (for Electron Muon Channel) ; pT (4th Leading Jet) [GeV] ; Nr. of Events");
    canvas_emu.push_back("5th Leading Jet pT (for Electron Muon Channel) ; pT (5th Leading Jet) [GeV] ; Nr. of Events");
    canvas_emu.push_back("6th Leading Jet pT (for Electron Muon Channel) ; pT (6th Leading Jet) [GeV] ; Nr. of Events");
    canvas_emu.push_back("Jet pT (for Electron Muon Channel) ; pT (Jet) [GeV] ; Nr. of Events");
    canvas_emu.push_back("Jet CSV (for Electron Muon Channel) ; CSV (Jet) ; Nr. of Events");
    canvas_emu.push_back("HT (for Electron Muon Channel) ; HT [GeV] ; Nr. of Events");
    canvas_emu.push_back("Nr. of Jets (for Electron Muon Channel) ; Nr. of jets ; Nr. of Events");
    canvas_emu.push_back("Nr. of B-tags (for Electron Muon Channel) ; Nr. of b-tags ; Nr. of Events");
    canvas_emu.push_back("Nr. of Primary Vertices (for Electron Muon Channel) ; Nr. of primary vertices ; Nr. of Events");
    canvas_emu.push_back("Leading Gen Electron pT (for Electron Muon Channel) ; pT (Leading Gen Electron) [GeV] ; Nr. of Events");
    canvas_emu.push_back("Leading Gen Electron eta (for Electron Muon Channel) ; eta (Leading Gen Electron) ; Nr. of Events");
    canvas_emu.push_back("Leading Gen Muon pT (for Electron Muon Channel) ; pT (Leading Gen Muon) [GeV] ; Nr. of Events");
    canvas_emu.push_back("Leading Gen Muon eta (for Electron Muon Channel) ; eta (Leading Gen Muon) ; Nr. of Events");
    canvas_emu.push_back("Leading Gen Jet pT (for Electron Muon Channel) ; pT (Leading Gen Jet) [GeV] ; Nr. of Events");
    canvas_emu.push_back("2nd Leading Gen Jet pT (for Electron Muon Channel) ; pT (2nd Leading Gen Jet) [GeV] ; Nr. of Events");
    canvas_emu.push_back("Gen HT (for Electron Muon Channel) ; Gen HT [GeV] ; Nr. of Events");
    canvas_emu.push_back("Nr. of Gen Jets(for Electron Muon Channel) ; Nr. of genjets ; Nr. of Events");
    canvas_emu_ratio.push_back("; pT (Electron) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_emu_ratio.push_back("; eta (Electron) ; Powheg-Pythia/Sherpa+OL");
    canvas_emu_ratio.push_back("; phi (Electron) ; Powheg-Pythia/Sherpa+OL");
    canvas_emu_ratio.push_back("; pT (Muon) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_emu_ratio.push_back("; eta (Muon) ; Powheg-Pythia/Sherpa+OL");
    canvas_emu_ratio.push_back("; phi (Muon) ; Powheg-Pythia/Sherpa+OL");
    canvas_emu_ratio.push_back("; pT (Leading Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_emu_ratio.push_back("; pT (2nd Leading Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_emu_ratio.push_back("; pT (3rd Leading Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_emu_ratio.push_back("; pT (4th Leading Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_emu_ratio.push_back("; pT (5th Leading Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_emu_ratio.push_back("; pT (6th Leading Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_emu_ratio.push_back("; pT (Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_emu_ratio.push_back("; CSV (Jet) ; Powheg-Pythia/Sherpa+OL");
    canvas_emu_ratio.push_back("; HT [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_emu_ratio.push_back("; Nr. of jets ; Powheg-Pythia/Sherpa+OL");
    canvas_emu_ratio.push_back("; Nr. of b-tags ; Powheg-Pythia/Sherpa+OL");
    canvas_emu_ratio.push_back("; Nr. of primary vertices ; Powheg-Pythia/Sherpa+OL");
    canvas_emu_ratio.push_back("; pT (Leading Gen Electron) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_emu_ratio.push_back("; eta (Leading Gen Electron) ; Powheg-Pythia/Sherpa+OL");
    canvas_emu_ratio.push_back("; pT (Leading Gen Muon) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_emu_ratio.push_back("; eta (Leading Gen Muon) ; Powheg-Pythia/Sherpa+OL");
    canvas_emu_ratio.push_back("; pT (Leading Gen Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_emu_ratio.push_back("; pT (2nd Leading Gen Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_emu_ratio.push_back("; HT [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_emu_ratio.push_back("; Nr. of genjets ; Powheg-Pythia/Sherpa+OL");
    canvas_emu_2D.push_back("2D : Njets & HT (for Electron Muon Channel) ; Nr. of jets ; HT [GeV]");

    canvas_mumu.push_back("Leading Muon pT (for Double Muon Channel) ; pT (Leading Muon) [GeV] ; Nr. of Events");
    canvas_mumu.push_back("Leading Muon eta (for Double Muon Channel) ; eta (Leading Muon) ; Nr. of Events");
    canvas_mumu.push_back("Leading Muon phi (for Double Muon Channel) ; phi (Leading Muon) ; Nr. of Events");
    canvas_mumu.push_back("Sub-leading Muon pT (for Double Muon Channel) ; pT (Sub-leading Muon) [GeV] ; Nr. of Events");
    canvas_mumu.push_back("Sub-leading Muon eta (for Double Muon Channel) ; eta (Sub-leading Muon) ; Nr. of Events");
    canvas_mumu.push_back("Sub-leading Muon phi (for Double Muon Channel) ; phi (Sub-leading Muon) ; Nr. of Events");
    canvas_mumu.push_back("Leading Jet pT (for Double Muon Channel) ; pT (Leading Jet) [GeV] ; Nr. of Events");
    canvas_mumu.push_back("2nd Leading Jet pT (for Double Muon Channel) ; pT (2nd Leading Jet) [GeV] ; Nr. of Events");
    canvas_mumu.push_back("3rd Leading Jet pT (for Double Muon Channel) ; pT (3rd Leading Jet) [GeV] ; Nr. of Events");
    canvas_mumu.push_back("4th Leading Jet pT (for Double Muon Channel) ; pT (4th Leading Jet) [GeV] ; Nr. of Events");
    canvas_mumu.push_back("5th Leading Jet pT (for Double Muon Channel) ; pT (5th Leading Jet) [GeV] ; Nr. of Events");
    canvas_mumu.push_back("6th Leading Jet pT (for Double Muon Channel) ; pT (6th Leading Jet) [GeV] ; Nr. of Events");
    canvas_mumu.push_back("Jet pT (for Double Muon Channel) ; pT (Jet) [GeV] ; Nr. of Events");
    canvas_mumu.push_back("Jet CSV (for Double Muon Channel) ; CSV (Jet) ; Nr. of Events");
    canvas_mumu.push_back("HT (for Double Muon Channel) ; HT [GeV] ; Nr. of Events");
    canvas_mumu.push_back("Nr. of Jets (for Double Muon Channel) ; Nr. of jets ; Nr. of Events");
    canvas_mumu.push_back("Nr. of B-tags (for Double Muon Channel) ; Nr. of b-tags ; Nr. of Events");
    canvas_mumu.push_back("Nr. of Primary Vertices (for Double Muon Channel) ; Nr. of primary vertices ; Nr. of Events");
    canvas_mumu.push_back("Leading Gen Muon pT (for Double Muon Channel) ; pT (Leading Gen Muon) [GeV] ; Nr. of Events");
    canvas_mumu.push_back("Leading Gen Muon eta (for Double Muon Channel) ; eta (Leading Gen Muon) ; Nr. of Events");
    canvas_mumu.push_back("Sub-leading Gen Muon pT (for Double Muon Channel) ; pT (Sub-leading Gen Muon) [GeV] ; Nr. of Events");
    canvas_mumu.push_back("Sub-leading Gen Muon eta (for Double Muon Channel) ; eta (Sub-leading Gen Muon) ; Nr. of Events");
    canvas_mumu.push_back("Leading Gen Jet pT (for Double Muon Channel) ; pT (Leading Gen Jet) [GeV] ; Nr. of Events");
    canvas_mumu.push_back("2nd Leading Gen Jet pT (for Double Muon Channel) ; pT (2nd Leading Gen Jet) [GeV] ; Nr. of Events");
    canvas_mumu.push_back("Gen HT (for Double Muon Channel) ; Gen HT [GeV] ; Nr. of Events");
    canvas_mumu.push_back("Nr. of Gen Jets(for Double Muon Channel) ; Nr. of genjets ; Nr. of Events");
    canvas_mumu_ratio.push_back("; pT (Leading Muon) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_mumu_ratio.push_back("; eta (Leading Muon) ; Powheg-Pythia/Sherpa+OL");
    canvas_mumu_ratio.push_back("; phi (Leading Muon) ; Powheg-Pythia/Sherpa+OL");
    canvas_mumu_ratio.push_back("; pT (Sub-leading Muon) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_mumu_ratio.push_back("; eta (Sub-leading Muon) ; Powheg-Pythia/Sherpa+OL");
    canvas_mumu_ratio.push_back("; phi (Sub-leading Muon) ; Powheg-Pythia/Sherpa+OL");
    canvas_mumu_ratio.push_back("; pT (Leading Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_mumu_ratio.push_back("; pT (2nd Leading Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_mumu_ratio.push_back("; pT (3rd Leading Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_mumu_ratio.push_back("; pT (4th Leading Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_mumu_ratio.push_back("; pT (5th Leading Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_mumu_ratio.push_back("; pT (6th Leading Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_mumu_ratio.push_back("; pT (Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_mumu_ratio.push_back("; CSV (Jet) ; Powheg-Pythia/Sherpa+OL");
    canvas_mumu_ratio.push_back("; HT [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_mumu_ratio.push_back("; Nr. of jets ; Powheg-Pythia/Sherpa+OL");
    canvas_mumu_ratio.push_back("; Nr. of b-tags ; Powheg-Pythia/Sherpa+OL");
    canvas_mumu_ratio.push_back("; Nr. of primary vertices ; Powheg-Pythia/Sherpa+OL");
    canvas_mumu_ratio.push_back("; pT (Leading Gen Muon) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_mumu_ratio.push_back("; eta (Leading Gen Muon) ; Powheg-Pythia/Sherpa+OL");
    canvas_mumu_ratio.push_back("; pT (Sub-leading Gen Muon) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_mumu_ratio.push_back("; eta (Sub-leading Gen Muon) ; Powheg-Pythia/Sherpa+OL");
    canvas_mumu_ratio.push_back("; pT (Leading Gen Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_mumu_ratio.push_back("; pT (2nd Leading Gen Jet) [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_mumu_ratio.push_back("; HT [GeV] ; Powheg-Pythia/Sherpa+OL");
    canvas_mumu_ratio.push_back("; Nr. of genjets ; Powheg-Pythia/Sherpa+OL");
    canvas_mumu_2D.push_back("2D : Njets & HT (for Double Muon Channel) ; Nr. of jets ; HT [GeV]");

    // Lepton pT
    x_l_sl.push_back(0);
    y_l_sl.push_back(0);
    x_u_sl.push_back(300);
    //y_u_sl.push_back(1000);
    y_u_sl.push_back(400);
    y_n_l_sl.push_back(0);
    y_n_u_sl.push_back(0.1);
    // Lepton eta
    x_l_sl.push_back(-3);
    y_l_sl.push_back(0);
    x_u_sl.push_back(3);
    //y_u_sl.push_back(2000);
    y_u_sl.push_back(600);
    y_n_l_sl.push_back(0);
    y_n_u_sl.push_back(0.15);
    // Lepton phi
    x_l_sl.push_back(-3);
    y_l_sl.push_back(0);
    x_u_sl.push_back(3);
    //y_u_sl.push_back(2000);
    y_u_sl.push_back(500);
    y_n_l_sl.push_back(0);
    y_n_u_sl.push_back(0.1);
    // Jet 1 pT
    x_l_sl.push_back(0);
    y_l_sl.push_back(0);
    x_u_sl.push_back(300);
    //y_u_sl.push_back(700);
    y_u_sl.push_back(250);
    y_n_l_sl.push_back(0);
    y_n_u_sl.push_back(0.1);
    // Jet 2 pT
    x_l_sl.push_back(0);
    y_l_sl.push_back(0);
    x_u_sl.push_back(300);
    //y_u_sl.push_back(700);
    y_u_sl.push_back(250);
    y_n_l_sl.push_back(0);
    y_n_u_sl.push_back(0.1);
    // Jet 3 pT
    x_l_sl.push_back(0);
    y_l_sl.push_back(0);
    x_u_sl.push_back(300);
    //y_u_sl.push_back(700);
    y_u_sl.push_back(250);
    y_n_l_sl.push_back(0);
    y_n_u_sl.push_back(0.1);
    // Jet 4 pT
    x_l_sl.push_back(0);
    y_l_sl.push_back(0);
    x_u_sl.push_back(300);
    //y_u_sl.push_back(700);
    y_u_sl.push_back(250);
    y_n_l_sl.push_back(0);
    y_n_u_sl.push_back(0.1);
    // Jet 5 pT
    x_l_sl.push_back(0);
    y_l_sl.push_back(0);
    x_u_sl.push_back(300);
    //y_u_sl.push_back(700);
    y_u_sl.push_back(250);
    y_n_l_sl.push_back(0);
    y_n_u_sl.push_back(0.1);
    // Jet 6 pT
    x_l_sl.push_back(0);
    y_l_sl.push_back(0);
    x_u_sl.push_back(300);
    //y_u_sl.push_back(700);
    y_u_sl.push_back(250);
    y_n_l_sl.push_back(0);
    y_n_u_sl.push_back(0.1);
    // Jet (all) pT
    x_l_sl.push_back(0);
    y_l_sl.push_back(0);
    x_u_sl.push_back(300);
    //y_u_sl.push_back(10000);
    y_u_sl.push_back(2500);
    y_n_l_sl.push_back(0);
    y_n_u_sl.push_back(0.1);
    // Jet (all) CSV
    x_l_sl.push_back(0);
    y_l_sl.push_back(0);
    x_u_sl.push_back(1);
    //y_u_sl.push_back(1000);
    y_u_sl.push_back(2500);
    y_n_l_sl.push_back(0);
    y_n_u_sl.push_back(0.16);
    // HT
    x_l_sl.push_back(0);
    y_l_sl.push_back(0);
    x_u_sl.push_back(1000);
    //y_u_sl.push_back(1500);
    y_u_sl.push_back(250);
    y_n_l_sl.push_back(0);
    y_n_u_sl.push_back(0.1);
    // Nr. of jets
    x_l_sl.push_back(4);
    y_l_sl.push_back(0);
    x_u_sl.push_back(13);
    //y_u_sl.push_back(10000);
    y_u_sl.push_back(1500);
    y_n_l_sl.push_back(0);
    y_n_u_sl.push_back(0.7);
    // Nr. of btags
    x_l_sl.push_back(2);
    y_l_sl.push_back(0);
    x_u_sl.push_back(9);
    //y_u_sl.push_back(10000);
    y_u_sl.push_back(4000);
    y_n_l_sl.push_back(0);
    y_n_u_sl.push_back(1);
    // Nr. of primary vertices
    x_l_sl.push_back(0);
    y_l_sl.push_back(0);
    x_u_sl.push_back(80);
    //y_u_sl.push_back(1000);
    y_u_sl.push_back(1000);
    y_n_l_sl.push_back(0);
    y_n_u_sl.push_back(0.2);
    // Nr. of jets & HT 2D
    x_l_sl_2D.push_back(4);
    y_l_sl_2D.push_back(0);
    x_u_sl_2D.push_back(13);
    y_u_sl_2D.push_back(1000);
    // Gen Lepton pT
    x_l_sl.push_back(0);
    y_l_sl.push_back(0);
    x_u_sl.push_back(300);
    //y_u_sl.push_back(1000);
    y_u_sl.push_back(500);
    y_n_l_sl.push_back(0);
    y_n_u_sl.push_back(0.1);
    // Gen Lepton eta
    x_l_sl.push_back(-3);
    y_l_sl.push_back(0);
    x_u_sl.push_back(3);
    //y_u_sl.push_back(2000);
    y_u_sl.push_back(6000);
    y_n_l_sl.push_back(0);
    y_n_u_sl.push_back(0.15);
    // Gen Jet 1 pT
    x_l_sl.push_back(0);
    y_l_sl.push_back(0);
    x_u_sl.push_back(300);
    //y_u_sl.push_back(700);
    y_u_sl.push_back(250);
    y_n_l_sl.push_back(0);
    y_n_u_sl.push_back(0.1);
    // Gen Jet 2 pT
    x_l_sl.push_back(0);
    y_l_sl.push_back(0);
    x_u_sl.push_back(300);
    //y_u_sl.push_back(700);
    y_u_sl.push_back(250);
    y_n_l_sl.push_back(0);
    y_n_u_sl.push_back(0.1);
    // Gen HT
    x_l_sl.push_back(0);
    y_l_sl.push_back(0);
    x_u_sl.push_back(1000);
    //y_u_sl.push_back(1500);
    y_u_sl.push_back(250);
    y_n_l_sl.push_back(0);
    y_n_u_sl.push_back(0.1);
    // Nr. of genjets
    x_l_sl.push_back(0);
    y_l_sl.push_back(0);
    x_u_sl.push_back(13);
    //y_u_sl.push_back(10000);
    y_u_sl.push_back(1500);
    y_n_l_sl.push_back(0);
    y_n_u_sl.push_back(0.7);

    // Lepton 1 pT
    x_l_di.push_back(0);
    y_l_di.push_back(0);
    x_u_di.push_back(300);
    //y_u_di.push_back(1000);
    y_u_di.push_back(75);
    y_n_l_di.push_back(0);
    y_n_u_di.push_back(0.1);
    // Lepton 1 eta
    x_l_di.push_back(-3);
    y_l_di.push_back(0);
    x_u_di.push_back(3);
    //y_u_di.push_back(2000);
    y_u_di.push_back(50);
    y_n_l_di.push_back(0);
    y_n_u_di.push_back(0.15);
    // Lepton 1 phi
    x_l_di.push_back(-3);
    y_l_di.push_back(0);
    x_u_di.push_back(3);
    //y_u_di.push_back(2000);
    y_u_di.push_back(100);
    y_n_l_di.push_back(0);
    y_n_u_di.push_back(0.1);
    // Lepton 2 pT
    x_l_di.push_back(0);
    y_l_di.push_back(0);
    x_u_di.push_back(300);
    //y_u_di.push_back(1000);
    y_u_di.push_back(75);
    y_n_l_di.push_back(0);
    y_n_u_di.push_back(0.1);
    // Lepton 2 eta
    x_l_di.push_back(-3);
    y_l_di.push_back(0);
    x_u_di.push_back(3);
    //y_u_di.push_back(2000);
    y_u_di.push_back(50);
    y_n_l_di.push_back(0);
    y_n_u_di.push_back(0.15);
    // Lepton 2 phi
    x_l_di.push_back(-3);
    y_l_di.push_back(0);
    x_u_di.push_back(3);
    //y_u_di.push_back(2000);
    y_u_di.push_back(100);
    y_n_l_di.push_back(0);
    y_n_u_di.push_back(0.1);
    // Jet 1 pT
    x_l_di.push_back(0);
    y_l_di.push_back(0);
    x_u_di.push_back(300);
    //y_u_di.push_back(700);
    y_u_di.push_back(75);
    y_n_l_di.push_back(0);
    y_n_u_di.push_back(0.1);
    // Jet 2 pT
    x_l_di.push_back(0);
    y_l_di.push_back(0);
    x_u_di.push_back(300);
    //y_u_di.push_back(700);
    y_u_di.push_back(75);
    y_n_l_di.push_back(0);
    y_n_u_di.push_back(0.1);
    // Jet 3 pT
    x_l_di.push_back(0);
    y_l_di.push_back(0);
    x_u_di.push_back(300);
    //y_u_di.push_back(700);
    y_u_di.push_back(75);
    y_n_l_di.push_back(0);
    y_n_u_di.push_back(0.1);
    // Jet 4 pT
    x_l_di.push_back(0);
    y_l_di.push_back(0);
    x_u_di.push_back(300);
    //y_u_di.push_back(700);
    y_u_di.push_back(75);
    y_n_l_di.push_back(0);
    y_n_u_di.push_back(0.1);
    // Jet 5 pT
    x_l_di.push_back(0);
    y_l_di.push_back(0);
    x_u_di.push_back(300);
    //y_u_di.push_back(700);
    y_u_di.push_back(75);
    y_n_l_di.push_back(0);
    y_n_u_di.push_back(0.1);
    // Jet 6 pT
    x_l_di.push_back(0);
    y_l_di.push_back(0);
    x_u_di.push_back(300);
    //y_u_di.push_back(700);
    y_u_di.push_back(75);
    y_n_l_di.push_back(0);
    y_n_u_di.push_back(0.1);
    // Jet (all) pT
    x_l_di.push_back(0);
    y_l_di.push_back(0);
    x_u_di.push_back(300);
    //y_u_di.push_back(10000);
    y_u_di.push_back(500);
    y_n_l_di.push_back(0);
    y_n_u_di.push_back(0.1);
    // Jet (all) CSV
    x_l_di.push_back(0);
    y_l_di.push_back(0);
    x_u_di.push_back(1);
    //y_u_di.push_back(1000);
    y_u_di.push_back(250);
    y_n_l_di.push_back(0);
    y_n_u_di.push_back(0.16);
    // HT
    x_l_di.push_back(0);
    y_l_di.push_back(0);
    x_u_di.push_back(1000);
    //y_u_di.push_back(1500);
    y_u_di.push_back(50);
    y_n_l_di.push_back(0);
    y_n_u_di.push_back(0.1);
    // Nr. of jets
    x_l_di.push_back(2);
    y_l_di.push_back(0);
    x_u_di.push_back(13);
    //y_u_di.push_back(10000);
    y_u_di.push_back(500);
    y_n_l_di.push_back(0);
    y_n_u_di.push_back(0.7);
    // Nr. of btags
    x_l_di.push_back(1);
    y_l_di.push_back(0);
    x_u_di.push_back(9);
    //y_u_di.push_back(10000);
    y_u_di.push_back(500);
    y_n_l_di.push_back(0);
    y_n_u_di.push_back(1);
    // Nr. of primary vertices
    x_l_di.push_back(0);
    y_l_di.push_back(0);
    x_u_di.push_back(80);
    //y_u_di.push_back(1000);
    y_u_di.push_back(100);
    y_n_l_di.push_back(0);
    y_n_u_di.push_back(0.2);
    // Nr. of jets & HT 2D
    x_l_di_2D.push_back(2);
    y_l_di_2D.push_back(0);
    x_u_di_2D.push_back(13);
    y_u_di_2D.push_back(1000);
    // Gen Lepton 1 pT
    x_l_di.push_back(0);
    y_l_di.push_back(0);
    x_u_di.push_back(300);
    //y_u_di.push_back(1000);
    y_u_di.push_back(75);
    y_n_l_di.push_back(0);
    y_n_u_di.push_back(0.1);
    // Gen Lepton 1 eta
    x_l_di.push_back(-3);
    y_l_di.push_back(0);
    x_u_di.push_back(3);
    //y_u_di.push_back(2000);
    y_u_di.push_back(50);
    y_n_l_di.push_back(0);
    y_n_u_di.push_back(0.15);
    // Gen Lepton 2 pT
    x_l_di.push_back(0);
    y_l_di.push_back(0);
    x_u_di.push_back(300);
    //y_u_di.push_back(1000);
    y_u_di.push_back(75);
    y_n_l_di.push_back(0);
    y_n_u_di.push_back(0.1);
    // Gen Lepton 2 eta
    x_l_di.push_back(-3);
    y_l_di.push_back(0);
    x_u_di.push_back(3);
    //y_u_di.push_back(2000);
    y_u_di.push_back(50);
    y_n_l_di.push_back(0);
    y_n_u_di.push_back(0.15);
    // Gen Jet 1 pT
    x_l_di.push_back(0);
    y_l_di.push_back(0);
    x_u_di.push_back(300);
    //y_u_di.push_back(700);
    y_u_di.push_back(75);
    y_n_l_di.push_back(0);
    y_n_u_di.push_back(0.1);
    // Gen Jet 2 pT
    x_l_di.push_back(0);
    y_l_di.push_back(0);
    x_u_di.push_back(300);
    //y_u_di.push_back(700);
    y_u_di.push_back(75);
    y_n_l_di.push_back(0);
    y_n_u_di.push_back(0.1);
    // Gen HT
    x_l_di.push_back(0);
    y_l_di.push_back(0);
    x_u_di.push_back(1000);
    //y_u_di.push_back(1500);
    y_u_di.push_back(50);
    y_n_l_di.push_back(0);
    y_n_u_di.push_back(0.1);
    // Nr. of genjets
    x_l_di.push_back(0);
    y_l_di.push_back(0);
    x_u_di.push_back(13);
    //y_u_di.push_back(10000);
    y_u_di.push_back(500);
    y_n_l_di.push_back(0);
    y_n_u_di.push_back(0.7);


    gStyle->SetLegendTextSize(0.033);
    gStyle->SetLabelSize(0.045,"X");
    gStyle->SetLabelSize(0.045,"Y");
    gStyle->SetTitleXSize(0.045);
    gStyle->SetTitleYSize(0.045);

    // Plotting : Not Normalized to Unity

    // Single Electron

    for(int i=0; i<n_e; i++){

        TCanvas *c1 = new TCanvas("c1","test",1100,650);
        //c1->DrawFrame(x_l_sl[i],y_l_sl[i],x_u_sl[i],y_u_sl[i],(canvas_e[i]).c_str());
        TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
        pad1->Draw();
        pad1->cd();
        pad1->DrawFrame(x_l_sl[i],y_l_sl[i],x_u_sl[i],y_u_sl[i],(canvas_e[i]).c_str());
        pad1->SetBottomMargin(0);

        TLegend* leg1 = new TLegend(0.55,0.65,0.85,0.85);
        leg1->SetFillColor(kWhite);
        leg1->SetFillStyle(1001);
        leg1->AddEntry(ttjets_pp_e[i],"MC : Background Powheg-Pythia","L");
        leg1->AddEntry(ttjets_ol_e[i],"MC : Background OpenLoops-Sherpa","L");
        //leg1->AddEntry(data_e[i],"Data","L");
        ttjets_pp_e[i]->SetLineColor(kRed);
        ttjets_pp_e[i]->SetLineWidth(2);
        ttjets_pp_e[i]->Draw("same");

        ttjets_ol_e[i]->SetLineColor(kGreen+3);
        ttjets_ol_e[i]->SetLineWidth(3);
        ttjets_ol_e[i]->Draw("same");

        //data_e[i]->SetLineColor(kBlue);
        //data_e[i]->SetLineWidth(2);
        //data_e[i]->Draw("same");

        leg1->Draw("same");

        c1->cd();
        TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.25);
        pad2->Draw();
        pad2->cd();
        pad2->DrawFrame(x_l_sl[i],0,x_u_sl[i],2,(canvas_e_ratio[i]).c_str());
        pad2->SetTopMargin(0);
        pad2->SetBottomMargin(0.3);
        TLine *line1 = new TLine(x_l_sl[i],1,x_u_sl[i],1);
        line1->SetLineStyle(kDashed);
        line1->SetLineWidth(1);
        TH1D *h = (TH1D*)ttjets_pp_e[i]->Clone("h");
        h->Divide(ttjets_ol_e[i]);
        h->SetStats(0);
        h->SetTitle("");
        h->SetLineColor(kBlack);
        h->SetMinimum(0);
        h->SetMaximum(2);
        h->GetXaxis()->SetRangeUser(x_l_sl[i],x_u_sl[i]);
        //h->GetXaxis()->SetTitle("X");
        //h->GetYaxis()->SetTitle("Powheg-Pythia/Sherpa+OL");
        h->GetYaxis()->SetNdivisions(4);
        h->GetXaxis()->SetTitleSize(0.15);
        h->GetYaxis()->SetTitleSize(0.12);
        h->GetXaxis()->SetLabelSize(0.13);
        h->GetYaxis()->SetLabelSize(0.13);
        h->Draw("ep");
        line1->Draw("same");
        c1->cd();

        c1->Print((out_file_e[i]).c_str());
        delete c1;
        delete leg1;
        delete line1;

    }

    for(int i=0; i<n_e_2D; i++){

        TCanvas *c1 = new TCanvas("c1","test",1100,650);
        c1->DrawFrame(x_l_sl_2D[i],y_l_sl_2D[i],x_u_sl_2D[i],y_u_sl_2D[i],(canvas_e_2D[i]).c_str());
        TLegend* leg1 = new TLegend(0.55,0.65,0.85,0.85);
        leg1->SetFillColor(kWhite);
        leg1->SetFillStyle(1001);
        leg1->AddEntry(ttjets_pp_e_2D[i],"MC : Background Powheg-Pythia","L");
        leg1->AddEntry(ttjets_ol_e_2D[i],"MC : Background OpenLoops-Sherpa","L");
        //leg1->AddEntry(data_e_2D[i],"Data","L");
        ttjets_pp_e_2D[i]->SetLineColor(kRed);
        ttjets_pp_e_2D[i]->SetLineWidth(2);
        ttjets_pp_e_2D[i]->Draw("same");

        ttjets_ol_e_2D[i]->SetLineColor(kGreen+3);
        ttjets_ol_e_2D[i]->SetLineWidth(3);
        ttjets_ol_e_2D[i]->Draw("same");

        //data_e_2D[i]->SetLineColor(kBlue);
        //data_e_2D[i]->SetLineWidth(2);
        //data_e_2D[i]->Draw("same");

        leg1->Draw("same");
        c1->Print((out_file_e_2D[i]).c_str());
        delete c1;
        delete leg1;
        
    }

    // Single Muon
    for(int i=0; i<n_mu; i++){

        TCanvas *c1 = new TCanvas("c1","test",1100,650);
        //c1->DrawFrame(x_l_sl[i],y_l_sl[i],x_u_sl[i],y_u_sl[i],(canvas_mu[i]).c_str());
        TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
        pad1->Draw();
        pad1->cd();
        pad1->DrawFrame(x_l_sl[i],y_l_sl[i],x_u_sl[i],y_u_sl[i],(canvas_mu[i]).c_str());
        pad1->SetBottomMargin(0);


        TLegend* leg1 = new TLegend(0.55,0.65,0.85,0.85);
        leg1->SetFillColor(kWhite);
        leg1->SetFillStyle(1001);
        leg1->AddEntry(ttjets_pp_mu[i],"MC : Background Powheg-Pythia","L");
        leg1->AddEntry(ttjets_ol_mu[i],"MC : Background OpenLoops-Sherpa","L");
        //leg1->AddEntry(data_mu[i],"Data","L");
        ttjets_pp_mu[i]->SetLineColor(kRed);
        ttjets_pp_mu[i]->SetLineWidth(2);
        ttjets_pp_mu[i]->Draw("same");

        ttjets_ol_mu[i]->SetLineColor(kGreen+3);
        ttjets_ol_mu[i]->SetLineWidth(3);
        ttjets_ol_mu[i]->Draw("same");

        //data_mu[i]->SetLineColor(kBlue);
        //data_mu[i]->SetLineWidth(2);
        //data_mu[i]->Draw("same");

        leg1->Draw("same");

        c1->cd();
        TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.25);
        pad2->Draw();
        pad2->cd();
        pad2->DrawFrame(x_l_sl[i],0,x_u_sl[i],2,(canvas_mu_ratio[i]).c_str());
        pad2->SetTopMargin(0);
        pad2->SetBottomMargin(0.3);
        TLine *line1 = new TLine(x_l_sl[i],1,x_u_sl[i],1);
        line1->SetLineStyle(kDashed);
        line1->SetLineWidth(1);
        TH1D *h = (TH1D*)ttjets_pp_mu[i]->Clone("h");
        h->Divide(ttjets_ol_mu[i]);
        h->SetStats(0);
        h->SetTitle("");
        h->SetLineColor(kBlack);
        h->SetMinimum(0);
        h->SetMaximum(2);
        h->GetXaxis()->SetRangeUser(x_l_sl[i],x_u_sl[i]);
        //h->GetXaxis()->SetTitle("X");
        //h->GetYaxis()->SetTitle("Powheg-Pythia/Sherpa+OL");
        h->GetYaxis()->SetNdivisions(4);
        h->GetXaxis()->SetTitleSize(0.15);
        h->GetYaxis()->SetTitleSize(0.12);
        h->GetXaxis()->SetLabelSize(0.13);
        h->GetYaxis()->SetLabelSize(0.13);
        h->Draw("ep");
        line1->Draw("same");
        c1->cd();

        c1->Print((out_file_mu[i]).c_str());
        delete c1;
        delete leg1;
        delete line1;
        
    }

    for(int i=0; i<n_mu_2D; i++){

        TCanvas *c1 = new TCanvas("c1","test",1100,650);
        c1->DrawFrame(x_l_sl_2D[i],y_l_sl_2D[i],x_u_sl_2D[i],y_u_sl_2D[i],(canvas_mu_2D[i]).c_str());
        TLegend* leg1 = new TLegend(0.55,0.65,0.85,0.85);
        leg1->SetFillColor(kWhite);
        leg1->SetFillStyle(1001);
        leg1->AddEntry(ttjets_pp_mu_2D[i],"MC : Background Powheg-Pythia","L");
        leg1->AddEntry(ttjets_ol_mu_2D[i],"MC : Background OpenLoops-Sherpa","L");
        //leg1->AddEntry(data_mu_2D[i],"Data","L");
        ttjets_pp_mu_2D[i]->SetLineColor(kRed);
        ttjets_pp_mu_2D[i]->SetLineWidth(2);
        ttjets_pp_mu_2D[i]->Draw("same");

        ttjets_ol_mu_2D[i]->SetLineColor(kGreen+3);
        ttjets_ol_mu_2D[i]->SetLineWidth(3);
        ttjets_ol_mu_2D[i]->Draw("same");

        //data_mu_2D[i]->SetLineColor(kBlue);
        //data_mu_2D[i]->SetLineWidth(2);
        //data_mu_2D[i]->Draw("same");

        leg1->Draw("same");
        c1->Print((out_file_mu_2D[i]).c_str());
        delete c1;
        delete leg1;
        
    }

    // Double Electron
    for(int i=0; i<n_ee; i++){

        TCanvas *c1 = new TCanvas("c1","test",1100,650);
        //c1->DrawFrame(x_l_di[i],y_l_di[i],x_u_di[i],y_u_di[i],(canvas_ee[i]).c_str());
        TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
        pad1->Draw();
        pad1->cd();
        pad1->DrawFrame(x_l_di[i],y_l_di[i],x_u_di[i],y_u_di[i],(canvas_ee[i]).c_str());
        pad1->SetBottomMargin(0);

        TLegend* leg1 = new TLegend(0.55,0.65,0.85,0.85);
        leg1->SetFillColor(kWhite);
        leg1->SetFillStyle(1001);
        leg1->AddEntry(ttjets_pp_ee[i],"MC : Background Powheg-Pythia","L");
        leg1->AddEntry(ttjets_ol_ee[i],"MC : Background OpenLoops-Sherpa","L");
        //leg1->AddEntry(data_ee[i],"Data","L");
        ttjets_pp_ee[i]->SetLineColor(kRed);
        ttjets_pp_ee[i]->SetLineWidth(2);
        ttjets_pp_ee[i]->Draw("same");

        ttjets_ol_ee[i]->SetLineColor(kGreen+3);
        ttjets_ol_ee[i]->SetLineWidth(3);
        ttjets_ol_ee[i]->Draw("same");

        //data_ee[i]->SetLineColor(kBlue);
        //data_ee[i]->SetLineWidth(2);
        //data_ee[i]->Draw("same");

        leg1->Draw("same");

        c1->cd();
        TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.25);
        pad2->Draw();
        pad2->cd();
        pad2->DrawFrame(x_l_di[i],0,x_u_di[i],2,(canvas_ee_ratio[i]).c_str());
        pad2->SetTopMargin(0);
        pad2->SetBottomMargin(0.3);
        TLine *line1 = new TLine(x_l_di[i],1,x_u_di[i],1);
        line1->SetLineStyle(kDashed);
        line1->SetLineWidth(1);
        TH1D *h = (TH1D*)ttjets_pp_ee[i]->Clone("h");
        h->Divide(ttjets_ol_ee[i]);
        h->SetStats(0);
        h->SetTitle("");
        h->SetLineColor(kBlack);
        h->SetMinimum(0);
        h->SetMaximum(2);
        h->GetXaxis()->SetRangeUser(x_l_di[i],x_u_di[i]);
        //h->GetXaxis()->SetTitle("X");
        //h->GetYaxis()->SetTitle("Powheg-Pythia/Sherpa+OL");
        h->GetYaxis()->SetNdivisions(4);
        h->GetXaxis()->SetTitleSize(0.15);
        h->GetYaxis()->SetTitleSize(0.12);
        h->GetXaxis()->SetLabelSize(0.13);
        h->GetYaxis()->SetLabelSize(0.13);
        h->Draw("ep");
        line1->Draw("same");
        c1->cd();

        c1->Print((out_file_ee[i]).c_str());
        delete c1;
        delete leg1;
        delete line1;
        
    }

    for(int i=0; i<n_ee_2D; i++){

        TCanvas *c1 = new TCanvas("c1","test",1100,650);
        c1->DrawFrame(x_l_di_2D[i],y_l_di_2D[i],x_u_di_2D[i],y_u_di_2D[i],(canvas_ee_2D[i]).c_str());
        TLegend* leg1 = new TLegend(0.55,0.65,0.85,0.85);
        leg1->SetFillColor(kWhite);
        leg1->SetFillStyle(1001);
        leg1->AddEntry(ttjets_pp_ee_2D[i],"MC : Background Powheg-Pythia","L");
        leg1->AddEntry(ttjets_ol_ee_2D[i],"MC : Background OpenLoops-Sherpa","L");
        //leg1->AddEntry(data_ee_2D[i],"Data","L");
        ttjets_pp_ee_2D[i]->SetLineColor(kRed);
        ttjets_pp_ee_2D[i]->SetLineWidth(2);
        ttjets_pp_ee_2D[i]->Draw("same");

        ttjets_ol_ee_2D[i]->SetLineColor(kGreen+3);
        ttjets_ol_ee_2D[i]->SetLineWidth(3);
        ttjets_ol_ee_2D[i]->Draw("same");

        //data_ee_2D[i]->SetLineColor(kBlue);
        //data_ee_2D[i]->SetLineWidth(2);
        //data_ee_2D[i]->Draw("same");

        leg1->Draw("same");
        c1->Print((out_file_ee_2D[i]).c_str());
        delete c1;
        delete leg1;
        
    }

    // Electron Muon
    for(int i=0; i<n_emu; i++){

        TCanvas *c1 = new TCanvas("c1","test",1100,650);
        //c1->DrawFrame(x_l_di[i],y_l_di[i],x_u_di[i],y_u_di[i],(canvas_emu[i]).c_str());
        TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
        pad1->Draw();
        pad1->cd();
        pad1->DrawFrame(x_l_di[i],y_l_di[i],x_u_di[i],y_u_di[i],(canvas_emu[i]).c_str());
        pad1->SetBottomMargin(0);

        TLegend* leg1 = new TLegend(0.55,0.65,0.85,0.85);
        leg1->SetFillColor(kWhite);
        leg1->SetFillStyle(1001);
        leg1->AddEntry(ttjets_pp_emu[i],"MC : Background Powheg-Pythia","L");
        leg1->AddEntry(ttjets_ol_emu[i],"MC : Background OpenLoops-Sherpa","L");
        //leg1->AddEntry(data_emu[i],"Data","L");
        ttjets_pp_emu[i]->SetLineColor(kRed);
        ttjets_pp_emu[i]->SetLineWidth(2);
        ttjets_pp_emu[i]->Draw("same");

        ttjets_ol_emu[i]->SetLineColor(kGreen+3);
        ttjets_ol_emu[i]->SetLineWidth(3);
        ttjets_ol_emu[i]->Draw("same");

        //data_emu[i]->SetLineColor(kBlue);
        //data_emu[i]->SetLineWidth(2);
        //data_emu[i]->Draw("same");

        leg1->Draw("same");

        c1->cd();
        TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.25);
        pad2->Draw();
        pad2->cd();
        pad2->DrawFrame(x_l_di[i],0,x_u_di[i],2,(canvas_emu_ratio[i]).c_str());
        pad2->SetTopMargin(0);
        pad2->SetBottomMargin(0.3);
        TLine *line1 = new TLine(x_l_di[i],1,x_u_di[i],1);
        line1->SetLineStyle(kDashed);
        line1->SetLineWidth(1);
        TH1D *h = (TH1D*)ttjets_pp_emu[i]->Clone("h");
        h->Divide(ttjets_ol_emu[i]);
        h->SetStats(0);
        h->SetTitle("");
        h->SetLineColor(kBlack);
        h->SetMinimum(0);
        h->SetMaximum(2);
        h->GetXaxis()->SetRangeUser(x_l_di[i],x_u_di[i]);
        //h->GetXaxis()->SetTitle("X");
        //h->GetYaxis()->SetTitle("Powheg-Pythia/Sherpa+OL");
        h->GetYaxis()->SetNdivisions(4);
        h->GetXaxis()->SetTitleSize(0.15);
        h->GetYaxis()->SetTitleSize(0.12);
        h->GetXaxis()->SetLabelSize(0.13);
        h->GetYaxis()->SetLabelSize(0.13);
        h->Draw("ep");
        line1->Draw("same");
        c1->cd();

        c1->Print((out_file_emu[i]).c_str());
        delete c1;
        delete leg1;
        delete line1;
        
    }

    for(int i=0; i<n_emu_2D; i++){

        TCanvas *c1 = new TCanvas("c1","test",1100,650);
        c1->DrawFrame(x_l_di_2D[i],y_l_di_2D[i],x_u_di_2D[i],y_u_di_2D[i],(canvas_emu_2D[i]).c_str());
        TLegend* leg1 = new TLegend(0.55,0.65,0.85,0.85);
        leg1->SetFillColor(kWhite);
        leg1->SetFillStyle(1001);
        leg1->AddEntry(ttjets_pp_emu_2D[i],"MC : Background Powheg-Pythia","L");
        leg1->AddEntry(ttjets_ol_emu_2D[i],"MC : Background OpenLoops-Sherpa","L");
        //leg1->AddEntry(data_emu_2D[i],"Data","L");
        ttjets_pp_emu_2D[i]->SetLineColor(kRed);
        ttjets_pp_emu_2D[i]->SetLineWidth(2);
        ttjets_pp_emu_2D[i]->Draw("same");

        ttjets_ol_emu_2D[i]->SetLineColor(kGreen+3);
        ttjets_ol_emu_2D[i]->SetLineWidth(3);
        ttjets_ol_emu_2D[i]->Draw("same");

        //data_emu_2D[i]->SetLineColor(kBlue);
        //data_emu_2D[i]->SetLineWidth(2);
        //data_emu_2D[i]->Draw("same");

        leg1->Draw("same");
        c1->Print((out_file_emu_2D[i]).c_str());
        delete c1;
        delete leg1;
        
    }

    // Double Muon
    for(int i=0; i<n_mumu; i++){

        TCanvas *c1 = new TCanvas("c1","test",1100,650);
        //c1->DrawFrame(x_l_di[i],y_l_di[i],x_u_di[i],y_u_di[i],(canvas_mumu[i]).c_str());
        TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
        pad1->Draw();
        pad1->cd();
        pad1->DrawFrame(x_l_di[i],y_l_di[i],x_u_di[i],y_u_di[i],(canvas_mumu[i]).c_str());
        pad1->SetBottomMargin(0);

        TLegend* leg1 = new TLegend(0.55,0.65,0.85,0.85);
        leg1->SetFillColor(kWhite);
        leg1->SetFillStyle(1001);
        leg1->AddEntry(ttjets_pp_mumu[i],"MC : Background Powheg-Pythia","L");
        leg1->AddEntry(ttjets_ol_mumu[i],"MC : Background OpenLoops-Sherpa","L");
        //leg1->AddEntry(data_mumu[i],"Data","L");
        ttjets_pp_mumu[i]->SetLineColor(kRed);
        ttjets_pp_mumu[i]->SetLineWidth(2);
        ttjets_pp_mumu[i]->Draw("same");

        ttjets_ol_mumu[i]->SetLineColor(kGreen+3);
        ttjets_ol_mumu[i]->SetLineWidth(3);
        ttjets_ol_mumu[i]->Draw("same");

        //data_mumu[i]->SetLineColor(kBlue);
        //data_mumu[i]->SetLineWidth(2);
        //data_mumu[i]->Draw("same");

        leg1->Draw("same");

        c1->cd();
        TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.25);
        pad2->Draw();
        pad2->cd();
        pad2->DrawFrame(x_l_di[i],0,x_u_di[i],2,(canvas_mumu_ratio[i]).c_str());
        pad2->SetTopMargin(0);
        pad2->SetBottomMargin(0.3);
        TLine *line1 = new TLine(x_l_di[i],1,x_u_di[i],1);
        line1->SetLineStyle(kDashed);
        line1->SetLineWidth(1);
        TH1D *h = (TH1D*)ttjets_pp_mumu[i]->Clone("h");
        h->Divide(ttjets_ol_mumu[i]);
        h->SetStats(0);
        h->SetTitle("");
        h->SetLineColor(kBlack);
        h->SetMinimum(0);
        h->SetMaximum(2);
        h->GetXaxis()->SetRangeUser(x_l_di[i],x_u_di[i]);
        //h->GetXaxis()->SetTitle("X");
        //h->GetYaxis()->SetTitle("Powheg-Pythia/Sherpa+OL");
        h->GetYaxis()->SetNdivisions(4);
        h->GetXaxis()->SetTitleSize(0.15);
        h->GetYaxis()->SetTitleSize(0.12);
        h->GetXaxis()->SetLabelSize(0.13);
        h->GetYaxis()->SetLabelSize(0.13);
        h->Draw("ep");
        line1->Draw("same");
        c1->cd();

        c1->Print((out_file_mumu[i]).c_str());
        delete c1;
        delete leg1;
        delete line1;
        
    }

    for(int i=0; i<n_mumu_2D; i++){

        TCanvas *c1 = new TCanvas("c1","test",1100,650);
        c1->DrawFrame(x_l_di_2D[i],y_l_di_2D[i],x_u_di_2D[i],y_u_di_2D[i],(canvas_mumu_2D[i]).c_str());
        TLegend* leg1 = new TLegend(0.55,0.65,0.85,0.85);
        leg1->SetFillColor(kWhite);
        leg1->SetFillStyle(1001);
        leg1->AddEntry(ttjets_pp_mumu_2D[i],"MC : Background Powheg-Pythia","L");
        leg1->AddEntry(ttjets_ol_mumu_2D[i],"MC : Background OpenLoops-Sherpa","L");
        //leg1->AddEntry(data_mumu_2D[i],"Data","L");
        ttjets_pp_mumu_2D[i]->SetLineColor(kRed);
        ttjets_pp_mumu_2D[i]->SetLineWidth(2);
        ttjets_pp_mumu_2D[i]->Draw("same");

        ttjets_ol_mumu_2D[i]->SetLineColor(kGreen+3);
        ttjets_ol_mumu_2D[i]->SetLineWidth(3);
        ttjets_ol_mumu_2D[i]->Draw("same");

        //data_mumu_2D[i]->SetLineColor(kBlue);
        //data_mumu_2D[i]->SetLineWidth(2);
        //data_mumu_2D[i]->Draw("same");

        leg1->Draw("same");
        c1->Print((out_file_mumu_2D[i]).c_str());
        delete c1;
        delete leg1;
        
    }


    // Plotting : Normalized to unity

    // Single Electron
    for(int i=0; i<n_e; i++){

        TCanvas *c1 = new TCanvas("c1","test",1100,650);
        //c1->DrawFrame(x_l_sl[i],y_n_l_sl[i],x_u_sl[i],y_n_u_sl[i],(canvas_e[i]).c_str());
        TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
        pad1->Draw();
        pad1->cd();
        pad1->DrawFrame(x_l_sl[i],y_n_l_sl[i],x_u_sl[i],y_n_u_sl[i],(canvas_e[i]).c_str());
        pad1->SetBottomMargin(0);

        TLegend* leg1 = new TLegend(0.55,0.65,0.85,0.85);
        leg1->SetFillColor(kWhite);
        leg1->SetFillStyle(1001);
        leg1->AddEntry(ttjets_pp_e_norm[i],"MC : Background Powheg-Pythia","L");
        leg1->AddEntry(ttjets_ol_e_norm[i],"MC : Background OpenLoops-Sherpa","L");
        //leg1->AddEntry(data_e[i],"Data","L");
        ttjets_pp_e_norm[i]->SetLineColor(kRed);
        ttjets_pp_e_norm[i]->SetLineWidth(2);
        ttjets_pp_e_norm[i]->Draw("same");

        ttjets_ol_e_norm[i]->SetLineColor(kGreen+3);
        ttjets_ol_e_norm[i]->SetLineWidth(3);
        ttjets_ol_e_norm[i]->Draw("same");

        //data_e[i]->SetLineColor(kBlue);
        //data_e[i]->SetLineWidth(2);
        //data_e[i]->Draw("same");

        leg1->Draw("same");

        c1->cd();
        TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.25);
        pad2->Draw();
        pad2->cd();
        pad2->DrawFrame(x_l_sl[i],0,x_u_sl[i],2,(canvas_e_ratio[i]).c_str());
        pad2->SetTopMargin(0);
        pad2->SetBottomMargin(0.3);
        TLine *line1 = new TLine(x_l_sl[i],1,x_u_sl[i],1);
        line1->SetLineStyle(kDashed);
        line1->SetLineWidth(1);
        TH1D *h = (TH1D*)ttjets_pp_e_norm[i]->Clone("h");
        h->Divide(ttjets_ol_e_norm[i]);
        h->SetStats(0);
        h->SetTitle("");
        h->SetLineColor(kBlack);
        h->SetMinimum(0);
        h->SetMaximum(2);
        h->GetXaxis()->SetRangeUser(x_l_sl[i],x_u_sl[i]);
        //h->GetXaxis()->SetTitle("X");
        //h->GetYaxis()->SetTitle("Powheg-Pythia/Sherpa+OL");
        h->GetYaxis()->SetNdivisions(4);
        h->GetXaxis()->SetTitleSize(0.15);
        h->GetYaxis()->SetTitleSize(0.12);
        h->GetXaxis()->SetLabelSize(0.13);
        h->GetYaxis()->SetLabelSize(0.13);
        h->Draw("ep");
        line1->Draw("same");
        c1->cd();

        c1->Print(("norm_unity_" + out_file_e[i]).c_str());
        delete c1;
        delete leg1;
        delete line1;

    }

    // Single Muon
    for(int i=0; i<n_mu; i++){

        TCanvas *c1 = new TCanvas("c1","test",1100,650);
        //c1->DrawFrame(x_l_sl[i],y_n_l_sl[i],x_u_sl[i],y_n_u_sl[i],(canvas_mu[i]).c_str());
        TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
        pad1->Draw();
        pad1->cd();
        pad1->DrawFrame(x_l_sl[i],y_n_l_sl[i],x_u_sl[i],y_n_u_sl[i],(canvas_mu[i]).c_str());
        pad1->SetBottomMargin(0);


        TLegend* leg1 = new TLegend(0.55,0.65,0.85,0.85);
        leg1->SetFillColor(kWhite);
        leg1->SetFillStyle(1001);
        leg1->AddEntry(ttjets_pp_mu_norm[i],"MC : Background Powheg-Pythia","L");
        leg1->AddEntry(ttjets_ol_mu_norm[i],"MC : Background OpenLoops-Sherpa","L");
        //leg1->AddEntry(data_mu[i],"Data","L");
        ttjets_pp_mu_norm[i]->SetLineColor(kRed);
        ttjets_pp_mu_norm[i]->SetLineWidth(2);
        ttjets_pp_mu_norm[i]->Draw("same");

        ttjets_ol_mu_norm[i]->SetLineColor(kGreen+3);
        ttjets_ol_mu_norm[i]->SetLineWidth(3);
        ttjets_ol_mu_norm[i]->Draw("same");

        //data_mu[i]->SetLineColor(kBlue);
        //data_mu[i]->SetLineWidth(2);
        //data_mu[i]->Draw("same");

        leg1->Draw("same");

        c1->cd();
        TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.25);
        pad2->Draw();
        pad2->cd();
        pad2->DrawFrame(x_l_sl[i],0,x_u_sl[i],2,(canvas_mu_ratio[i]).c_str());
        pad2->SetTopMargin(0);
        pad2->SetBottomMargin(0.3);
        TLine *line1 = new TLine(x_l_sl[i],1,x_u_sl[i],1);
        line1->SetLineStyle(kDashed);
        line1->SetLineWidth(1);
        TH1D *h = (TH1D*)ttjets_pp_mu_norm[i]->Clone("h");
        h->Divide(ttjets_ol_mu_norm[i]);
        h->SetStats(0);
        h->SetTitle("");
        h->SetLineColor(kBlack);
        h->SetMinimum(0);
        h->SetMaximum(2);
        h->GetXaxis()->SetRangeUser(x_l_sl[i],x_u_sl[i]);
        //h->GetXaxis()->SetTitle("X");
        //h->GetYaxis()->SetTitle("Powheg-Pythia/Sherpa+OL");
        h->GetYaxis()->SetNdivisions(4);
        h->GetXaxis()->SetTitleSize(0.15);
        h->GetYaxis()->SetTitleSize(0.12);
        h->GetXaxis()->SetLabelSize(0.13);
        h->GetYaxis()->SetLabelSize(0.13);
        h->Draw("ep");
        line1->Draw("same");
        c1->cd();

        c1->Print(("norm_unity_" + out_file_mu[i]).c_str());
        delete c1;
        delete leg1;
        delete line1;

    }

    // Double Electron
    for(int i=0; i<n_ee; i++){

        TCanvas *c1 = new TCanvas("c1","test",1100,650);
        //c1->DrawFrame(x_l_di[i],y_n_l_di[i],x_u_di[i],y_n_u_di[i],(canvas_ee[i]).c_str());
        TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
        pad1->Draw();
        pad1->cd();
        pad1->DrawFrame(x_l_di[i],y_n_l_di[i],x_u_di[i],y_n_u_di[i],(canvas_ee[i]).c_str());
        pad1->SetBottomMargin(0);

        TLegend* leg1 = new TLegend(0.55,0.65,0.85,0.85);
        leg1->SetFillColor(kWhite);
        leg1->SetFillStyle(1001);
        leg1->AddEntry(ttjets_pp_ee_norm[i],"MC : Background Powheg-Pythia","L");
        leg1->AddEntry(ttjets_ol_ee_norm[i],"MC : Background OpenLoops-Sherpa","L");
        //leg1->AddEntry(data_ee[i],"Data","L");
        ttjets_pp_ee_norm[i]->SetLineColor(kRed);
        ttjets_pp_ee_norm[i]->SetLineWidth(2);
        ttjets_pp_ee_norm[i]->Draw("same");

        ttjets_ol_ee_norm[i]->SetLineColor(kGreen+3);
        ttjets_ol_ee_norm[i]->SetLineWidth(3);
        ttjets_ol_ee_norm[i]->Draw("same");

        //data_ee[i]->SetLineColor(kBlue);
        //data_ee[i]->SetLineWidth(2);
        //data_ee[i]->Draw("same");

        leg1->Draw("same");

        c1->cd();
        TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.25);
        pad2->Draw();
        pad2->cd();
        pad2->DrawFrame(x_l_di[i],0,x_u_di[i],2,(canvas_ee_ratio[i]).c_str());
        pad2->SetTopMargin(0);
        pad2->SetBottomMargin(0.3);
        TLine *line1 = new TLine(x_l_di[i],1,x_u_di[i],1);
        line1->SetLineStyle(kDashed);
        line1->SetLineWidth(1);
        TH1D *h = (TH1D*)ttjets_pp_ee_norm[i]->Clone("h");
        h->Divide(ttjets_ol_ee_norm[i]);
        h->SetStats(0);
        h->SetTitle("");
        h->SetLineColor(kBlack);
        h->SetMinimum(0);
        h->SetMaximum(2);
        h->GetXaxis()->SetRangeUser(x_l_di[i],x_u_di[i]);
        //h->GetXaxis()->SetTitle("X");
        //h->GetYaxis()->SetTitle("Powheg-Pythia/Sherpa+OL");
        h->GetYaxis()->SetNdivisions(4);
        h->GetXaxis()->SetTitleSize(0.15);
        h->GetYaxis()->SetTitleSize(0.12);
        h->GetXaxis()->SetLabelSize(0.13);
        h->GetYaxis()->SetLabelSize(0.13);
        h->Draw("ep");
        line1->Draw("same");
        c1->cd();

        c1->Print(("norm_unity_" + out_file_ee[i]).c_str());
        delete c1;
        delete leg1;
        delete line1;

    }


    // Electron Muon
    for(int i=0; i<n_emu; i++){

        TCanvas *c1 = new TCanvas("c1","test",1100,650);
        //c1->DrawFrame(x_l_di[i],y_n_l_di[i],x_u_di[i],y_n_u_di[i],(canvas_emu[i]).c_str());
        TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
        pad1->Draw();
        pad1->cd();
        pad1->DrawFrame(x_l_di[i],y_n_l_di[i],x_u_di[i],y_n_u_di[i],(canvas_emu[i]).c_str());
        pad1->SetBottomMargin(0);

        TLegend* leg1 = new TLegend(0.55,0.65,0.85,0.85);
        leg1->SetFillColor(kWhite);
        leg1->SetFillStyle(1001);
        leg1->AddEntry(ttjets_pp_emu_norm[i],"MC : Background Powheg-Pythia","L");
        leg1->AddEntry(ttjets_ol_emu_norm[i],"MC : Background OpenLoops-Sherpa","L");
        //leg1->AddEntry(data_emu[i],"Data","L");
        ttjets_pp_emu_norm[i]->SetLineColor(kRed);
        ttjets_pp_emu_norm[i]->SetLineWidth(2);
        ttjets_pp_emu_norm[i]->Draw("same");

        ttjets_ol_emu_norm[i]->SetLineColor(kGreen+3);
        ttjets_ol_emu_norm[i]->SetLineWidth(3);
        ttjets_ol_emu_norm[i]->Draw("same");

        //data_emu[i]->SetLineColor(kBlue);
        //data_emu[i]->SetLineWidth(2);
        //data_emu[i]->Draw("same");

        leg1->Draw("same");

        c1->cd();
        TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.25);
        pad2->Draw();
        pad2->cd();
        pad2->DrawFrame(x_l_di[i],0,x_u_di[i],2,(canvas_emu_ratio[i]).c_str());
        pad2->SetTopMargin(0);
        pad2->SetBottomMargin(0.3);
        TLine *line1 = new TLine(x_l_di[i],1,x_u_di[i],1);
        line1->SetLineStyle(kDashed);
        line1->SetLineWidth(1);
        TH1D *h = (TH1D*)ttjets_pp_emu_norm[i]->Clone("h");
        h->Divide(ttjets_ol_emu_norm[i]);
        h->SetStats(0);
        h->SetTitle("");
        h->SetLineColor(kBlack);
        h->SetMinimum(0);
        h->SetMaximum(2);
        h->GetXaxis()->SetRangeUser(x_l_di[i],x_u_di[i]);
        //h->GetXaxis()->SetTitle("X");
        //h->GetYaxis()->SetTitle("Powheg-Pythia/Sherpa+OL");
        h->GetYaxis()->SetNdivisions(4);
        h->GetXaxis()->SetTitleSize(0.15);
        h->GetYaxis()->SetTitleSize(0.12);
        h->GetXaxis()->SetLabelSize(0.13);
        h->GetYaxis()->SetLabelSize(0.13);
        h->Draw("ep");
        line1->Draw("same");
        c1->cd();

        c1->Print(("norm_unity_" + out_file_emu[i]).c_str());
        delete c1;
        delete leg1;
        delete line1;

    }

    // Double Muon
    for(int i=0; i<n_mumu; i++){

        TCanvas *c1 = new TCanvas("c1","test",1100,650);
        //c1->DrawFrame(x_l_di[i],y_n_l_di[i],x_u_di[i],y_n_u_di[i],(canvas_mumu[i]).c_str());
        TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
        pad1->Draw();
        pad1->cd();
        pad1->DrawFrame(x_l_di[i],y_n_l_di[i],x_u_di[i],y_n_u_di[i],(canvas_mumu[i]).c_str());
        pad1->SetBottomMargin(0);

        TLegend* leg1 = new TLegend(0.55,0.65,0.85,0.85);
        leg1->SetFillColor(kWhite);
        leg1->SetFillStyle(1001);
        leg1->AddEntry(ttjets_pp_mumu_norm[i],"MC : Background Powheg-Pythia","L");
        leg1->AddEntry(ttjets_ol_mumu_norm[i],"MC : Background OpenLoops-Sherpa","L");
        //leg1->AddEntry(data_mumu[i],"Data","L");
        ttjets_pp_mumu_norm[i]->SetLineColor(kRed);
        ttjets_pp_mumu_norm[i]->SetLineWidth(2);
        ttjets_pp_mumu_norm[i]->Draw("same");

        ttjets_ol_mumu_norm[i]->SetLineColor(kGreen+3);
        ttjets_ol_mumu_norm[i]->SetLineWidth(3);
        ttjets_ol_mumu_norm[i]->Draw("same");

        //data_mumu[i]->SetLineColor(kBlue);
        //data_mumu[i]->SetLineWidth(2);
        //data_mumu[i]->Draw("same");

        leg1->Draw("same");

        c1->cd();
        TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.25);
        pad2->Draw();
        pad2->cd();
        pad2->DrawFrame(x_l_di[i],0,x_u_di[i],2,(canvas_mumu_ratio[i]).c_str());
        pad2->SetTopMargin(0);
        pad2->SetBottomMargin(0.3);
        TLine *line1 = new TLine(x_l_di[i],1,x_u_di[i],1);
        line1->SetLineStyle(kDashed);
        line1->SetLineWidth(1);
        TH1D *h = (TH1D*)ttjets_pp_mumu_norm[i]->Clone("h");
        h->Divide(ttjets_ol_mumu_norm[i]);
        h->SetStats(0);
        h->SetTitle("");
        h->SetLineColor(kBlack);
        h->SetMinimum(0);
        h->SetMaximum(2);
        h->GetXaxis()->SetRangeUser(x_l_di[i],x_u_di[i]);
        //h->GetXaxis()->SetTitle("X");
        //h->GetYaxis()->SetTitle("Powheg-Pythia/Sherpa+OL");
        h->GetYaxis()->SetNdivisions(4);
        h->GetXaxis()->SetTitleSize(0.15);
        h->GetYaxis()->SetTitleSize(0.12);
        h->GetXaxis()->SetLabelSize(0.13);
        h->GetYaxis()->SetLabelSize(0.13);
        h->Draw("ep");
        line1->Draw("same");
        c1->cd();

        c1->Print(("norm_unity_" + out_file_mumu[i]).c_str());
        delete c1;
        delete leg1;
        delete line1;

    }


    f1->Close();
    f2->Close();
    //f3->Close();
	delete f1;
    delete f2;
    //delete f3;

	return;
}


