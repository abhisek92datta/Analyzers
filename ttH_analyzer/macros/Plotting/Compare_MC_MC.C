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


void Compare_Data_MC() {

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

    ///////////////////////////////////////////////////////////////////////
    // Reading the Histograms from the different files
    ///////////////////////////////////////////////////////////////////////

	TFile *f1 = new TFile("Distribution_mc_ttjets_pp.root");
    //TFile *f2 = new TFile("Distribution_mc_ttjets_ol.root");
    //TFile *f3 = new TFile("Distribution_data_2016.root");

    std::vector<TH1D*> ttjets_pp_e;
    std::vector<TH1D*> ttjets_pp_mu;
    std::vector<TH1D*> ttjets_pp_ee;
    std::vector<TH1D*> ttjets_pp_emu;
    std::vector<TH1D*> ttjets_pp_mumu;
    /*
    std::vector<TH1D*> ttjets_ol_e;
    std::vector<TH1D*> ttjets_ol_mu;
    std::vector<TH1D*> ttjets_ol_ee;
    std::vector<TH1D*> ttjets_ol_emu;
    std::vector<TH1D*> ttjets_ol_mumu;
    std::vector<TH1D*> data_e;
    std::vector<TH1D*> data_mu;
    std::vector<TH1D*> data_ee;
    std::vector<TH1D*> data_emu;
    std::vector<TH1D*> data_mumu;
    */

    int n_e, n_mu, n_ee, n_emu, n_mumu;
    n_e = n_mu = n_ee = n_emu = n_mumu = 0;

    for(int i=0; i<nhistos; i++){

        std::string histo = histonames[i];
        TH1D *h1;
        //TH1D *h2;
        //TH1D *h3;

        // ttjets_pp
        h1 = (TH1D*)f1->Get((histo).c_str());
        h1->SetName(("ttjets_pp_" + histo).c_str());
        h1->SetTitle(("MC ttjets_pp : " + histo).c_str());

        /*
        // ttjets_ol
        h2 = (TH1D*)f2->Get((histo).c_str());
        h2->SetName(("ttjets_ol_" + histo).c_str());
        h2->SetTitle(("MC ttjets_ol : " + histo).c_str());

        // data
        h3 = (TH1D*)f3->Get((histo).c_str());
        h3->SetName(("data_" + histo).c_str());
        h3->SetTitle(("DATA : " + histo).c_str());
        */

        if(histo.find("_sl_e") != std::string::npos){
            ttjets_pp_e.push_back(h1);
            //ttjets_ol_e.push_back(h2);
            //data_e.push_back(h3);

            out_file_e.push_back((histo + ".png").c_str());
            n_e++;
        }

        else if(histo.find("_sl_mu") != std::string::npos){
            ttjets_pp_mu.push_back(h1);
            //ttjets_ol_mu.push_back(h2);
            //data_mu.push_back(h3);

            out_file_mu.push_back((histo + ".png").c_str());
            n_mu++;
        }

        else if(histo.find("_di_ee") != std::string::npos){
            ttjets_pp_ee.push_back(h1);
            //ttjets_ol_ee.push_back(h2);
            //data_ee.push_back(h3);

            out_file_ee.push_back((histo + ".png").c_str());
            n_ee++;
        }

        else if(histo.find("_di_emu") != std::string::npos){
            ttjets_pp_emu.push_back(h1);
            //ttjets_ol_emu.push_back(h2);
            //data_emu.push_back(h3);

            out_file_emu.push_back((histo + ".png").c_str());
            n_emu++;
        }

        else if(histo.find("_di_mumu") != std::string::npos){
            ttjets_pp_mumu.push_back(h1);
            //ttjets_ol_mumu.push_back(h2);
            //data_mumu.push_back(h3);

            out_file_mumu.push_back((histo + ".png").c_str());
            n_mumu++;
        }

    }


    ///////////////////////////////////////////////////////////////////////
    // Scaling and Normlization of MC signal and background
    ///////////////////////////////////////////////////////////////////////

    //ttjets_pp
    double factor_ttjets_pp = 1.0;
    double sigma_ttjets_pp = ; // pb
    double N_total_ttjets_pp = ;
    double sum_gen_weight_ttjets_pp = ;
    double norm_ttjets_pp = (L*sigma_ttjets_pp*1000)/(sum_gen_weight_ttjets_pp);
    double scale_ttjets_pp = factor_ttjets_pp*norm_ttjets_pp;

    /*
    //ttjets_ol
    double factor_ttjets_ol = 1.0;
    double sigma_ttjets_ol = 831.76; // pb
    double N_total_ttjets_ol = ;
    double sum_gen_weight_ttjets_ol = ;
    double norm_ttjets_ol = (L*sigma_ttjets_ol*1000)/(sum_gen_weight_ttjets_ol);
    double scale_ttjets_ol = factor_ttjets_ol*norm_ttjets_ol;
    */

    for(int i=0; i<n_e; i++){
        ttjets_pp_e[i]->Scale(scale_ttjets_pp);
        //ttjets_ol_e[i]->Scale(scale_ttjets_ol);
    }
    for(int i=0; i<n_mu; i++){
        ttjets_pp_mu[i]->Scale(scale_ttjets_pp);
        //ttjets_ol_mu[i]->Scale(scale_ttjets_ol);
    }
    for(int i=0; i<n_ee; i++){
        ttjets_pp_ee[i]->Scale(scale_ttjets_pp);
        //ttjets_ol_ee[i]->Scale(scale_ttjets_ol);
    }
    for(int i=0; i<n_emu; i++){
        ttjets_pp_emu[i]->Scale(scale_ttjets_pp);
        //ttjets_ol_emu[i]->Scale(scale_ttjets_ol);
    }
    for(int i=0; i<n_mumu; i++){
        ttjets_pp_mumu[i]->Scale(scale_ttjets_pp);
        //ttjets_ol_mumu[i]->Scale(scale_ttjets_ol);
    }


    ///////////////////////////////////////////////////////////////////////
    // Plotting
    ///////////////////////////////////////////////////////////////////////

    // setting plot parameters
    std::vector<double> x_l_sl, x_u_sl, y_l_sl, y_u_sl;
    std::vector<double> x_l_di, x_u_di, y_l_di, y_u_di;
    std::vector<std::string> canvas_e;
    std::vector<std::string> canvas_mu;
    std::vector<std::string> canvas_ee;
    std::vector<std::string> canvas_emu;
    std::vector<std::string> canvas_mumu;

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

    canvas_emu.push_back("Electron pT (for Double Electron Channel) ; pT (Electron) [GeV] ; Nr. of Events");
    canvas_emu.push_back("Electron eta (for Double Electron Channel) ; eta (Electron) ; Nr. of Events");
    canvas_emu.push_back("Electron phi (for Double Electron Channel) ; phi (Electron) ; Nr. of Events");
    canvas_emu.push_back("Muon pT (for Double Electron Channel) ; pT (Muon) [GeV] ; Nr. of Events");
    canvas_emu.push_back("Muon eta (for Double Electron Channel) ; eta (Muon) ; Nr. of Events");
    canvas_emu.push_back("Muon phi (for Double Electron Channel) ; phi (Muon) ; Nr. of Events");
    canvas_emu.push_back("Leading Jet pT (for Double Electron Channel) ; pT (Leading Jet) [GeV] ; Nr. of Events");
    canvas_emu.push_back("2nd Leading Jet pT (for Double Electron Channel) ; pT (2nd Leading Jet) [GeV] ; Nr. of Events");
    canvas_emu.push_back("3rd Leading Jet pT (for Double Electron Channel) ; pT (3rd Leading Jet) [GeV] ; Nr. of Events");
    canvas_emu.push_back("4th Leading Jet pT (for Double Electron Channel) ; pT (4th Leading Jet) [GeV] ; Nr. of Events");
    canvas_emu.push_back("5th Leading Jet pT (for Double Electron Channel) ; pT (5th Leading Jet) [GeV] ; Nr. of Events");
    canvas_emu.push_back("6th Leading Jet pT (for Double Electron Channel) ; pT (6th Leading Jet) [GeV] ; Nr. of Events");
    canvas_emu.push_back("Jet pT (for Double Electron Channel) ; pT (Jet) [GeV] ; Nr. of Events");
    canvas_emu.push_back("Jet CSV (for Double Electron Channel) ; CSV (Jet) ; Nr. of Events");
    canvas_emu.push_back("HT (for Double Electron Channel) ; HT [GeV] ; Nr. of Events");
    canvas_emu.push_back("Nr. of Jets (for Double Electron Channel) ; Nr. of jets ; Nr. of Events");
    canvas_emu.push_back("Nr. of B-tags (for Double Electron Channel) ; Nr. of b-tags ; Nr. of Events");
    canvas_emu.push_back("Nr. of Primary Vertices (for Double Electron Channel) ; Nr. of primary vertices ; Nr. of Events");

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

    // Lepton pT
    x_l_sl.push_back(0);
    y_l_sl.push_back(0);
    x_u_sl.push_back(300);
    y_u_sl.push_back(13000);
    // Lepton eta
    x_l_sl.push_back(-3);
    y_l_sl.push_back(0);
    x_u_sl.push_back(3);
    y_u_sl.push_back(30000);
    // Lepton phi
    x_l_sl.push_back(-3);
    y_l_sl.push_back(4000);
    x_u_sl.push_back(3);
    y_u_sl.push_back(12000);
    // Jet 1 pT
    x_l_sl.push_back(0);
    y_l_sl.push_back(0);
    x_u_sl.push_back(300);
    y_u_sl.push_back(7000);
    // Jet 2 pT
    x_l_sl.push_back(0);
    y_l_sl.push_back(0);
    x_u_sl.push_back(300);
    y_u_sl.push_back(7000);
    // Jet 3 pT
    x_l_sl.push_back(0);
    y_l_sl.push_back(0);
    x_u_sl.push_back(300);
    y_u_sl.push_back(7000);
    // Jet 4 pT
    x_l_sl.push_back(0);
    y_l_sl.push_back(0);
    x_u_sl.push_back(300);
    y_u_sl.push_back(7000);
    // Jet 5 pT
    x_l_sl.push_back(0);
    y_l_sl.push_back(0);
    x_u_sl.push_back(300);
    y_u_sl.push_back(7000);
    // Jet 6 pT
    x_l_sl.push_back(0);
    y_l_sl.push_back(0);
    x_u_sl.push_back(300);
    y_u_sl.push_back(7000);
    // Jet (all) pT
    x_l_sl.push_back(0);
    y_l_sl.push_back(0);
    x_u_sl.push_back(300);
    y_u_sl.push_back(50000);
    // Jet (all) CSV
    x_l_sl.push_back(0);
    y_l_sl.push_back(0);
    x_u_sl.push_back(1);
    y_u_sl.push_back(30000);
    // HT
    x_l_sl.push_back(0);
    y_l_sl.push_back(0);
    x_u_sl.push_back(1000);
    y_u_sl.push_back(7000);
    // Nr. of jets
    x_l_sl.push_back(4);
    y_l_sl.push_back(0);
    x_u_sl.push_back(13);
    y_u_sl.push_back(200000);
    // Nr. of btags
    x_l_sl.push_back(2);
    y_l_sl.push_back(0);
    x_u_sl.push_back(9);
    y_u_sl.push_back(350000);
    // Nr. of primary vertices
    x_l_sl.push_back(0);
    y_l_sl.push_back(0);
    x_u_sl.push_back(80);
    y_u_sl.push_back(50000);

    // Lepton 1 pT
    x_l_di.push_back(0);
    y_l_di.push_back(0);
    x_u_di.push_back(300);
    y_u_di.push_back(13000);
    // Lepton 1 eta
    x_l_di.push_back(-3);
    y_l_di.push_back(0);
    x_u_di.push_back(3);
    y_u_di.push_back(30000);
    // Lepton 1 phi
    x_l_di.push_back(-3);
    y_l_di.push_back(4000);
    x_u_di.push_back(3);
    y_u_di.push_back(12000);
    // Lepton 2 pT
    x_l_di.push_back(0);
    y_l_di.push_back(0);
    x_u_di.push_back(300);
    y_u_di.push_back(13000);
    // Lepton 2 eta
    x_l_di.push_back(-3);
    y_l_di.push_back(0);
    x_u_di.push_back(3);
    y_u_di.push_back(30000);
    // Lepton 2 phi
    x_l_di.push_back(-3);
    y_l_di.push_back(4000);
    x_u_di.push_back(3);
    y_u_di.push_back(12000);
    // Jet 1 pT
    x_l_di.push_back(0);
    y_l_di.push_back(0);
    x_u_di.push_back(300);
    y_u_di.push_back(7000);
    // Jet 2 pT
    x_l_di.push_back(0);
    y_l_di.push_back(0);
    x_u_di.push_back(300);
    y_u_di.push_back(7000);
    // Jet 3 pT
    x_l_di.push_back(0);
    y_l_di.push_back(0);
    x_u_di.push_back(300);
    y_u_di.push_back(7000);
    // Jet 4 pT
    x_l_di.push_back(0);
    y_l_di.push_back(0);
    x_u_di.push_back(300);
    y_u_di.push_back(7000);
    // Jet 5 pT
    x_l_di.push_back(0);
    y_l_di.push_back(0);
    x_u_di.push_back(300);
    y_u_di.push_back(7000);
    // Jet 6 pT
    x_l_di.push_back(0);
    y_l_di.push_back(0);
    x_u_di.push_back(300);
    y_u_di.push_back(7000);
    // Jet (all) pT
    x_l_di.push_back(0);
    y_l_di.push_back(0);
    x_u_di.push_back(300);
    y_u_di.push_back(50000);
    // Jet (all) CSV
    x_l_di.push_back(0);
    y_l_di.push_back(0);
    x_u_di.push_back(1);
    y_u_di.push_back(30000);
    // HT
    x_l_di.push_back(0);
    y_l_di.push_back(0);
    x_u_di.push_back(1000);
    y_u_di.push_back(7000);
    // Nr. of jets
    x_l_di.push_back(4);
    y_l_di.push_back(0);
    x_u_di.push_back(13);
    y_u_di.push_back(200000);
    // Nr. of btags
    x_l_di.push_back(2);
    y_l_di.push_back(0);
    x_u_di.push_back(9);
    y_u_di.push_back(350000);
    // Nr. of primary vertices
    x_l_di.push_back(0);
    y_l_di.push_back(0);
    x_u_di.push_back(80);
    y_u_di.push_back(50000);


    // Plotting

    // Single Electron
    for(int i=0; i<n_e; i++){

        TCanvas *c1 = new TCanvas("c1","test",1100,650);
        c1->DrawFrame(x_l_sl[i],y_l_sl[i],x_u_sl[i],y_u_sl[i],(canvas_e[i]).c_str());
        TLegend* leg1 = new TLegend(0.65,0.70,0.85,0.85);
        leg1->SetFillColor(kWhite);
        leg1->SetFillStyle(1001);
        leg1->AddEntry(ttjets_pp_e[i],"MC : ttH x 600","L");
        //leg1->AddEntry(ttjets_ol_e[i],"MC : Background","L");
        //leg1->AddEntry(data_e[i],"Data","L");
        ttjets_pp_e[i]->SetLineColor(kRed);
        ttjets_pp_e[i]->SetLineWidth(2);
        ttjets_pp_e[i]->Draw("same");
        /*
        ttjets_ol_e[i]->SetLineColor(kGreen+3);
        ttjets_ol_e[i]->SetLineWidth(3);
        ttjets_ol_e[i]->Draw("same");
        data_e[i]->SetLineColor(kBlue);
        data_e[i]->SetLineWidth(2);
        data_e[i]->Draw("same");
        */
        leg1->Draw("same");
        c1->Print((out_file_e[i]).c_str());
        delete c1;
        delete leg1;

    }

    // Single Muon
    for(int i=0; i<n_mu; i++){

        TCanvas *c1 = new TCanvas("c1","test",1100,650);
        c1->DrawFrame(x_l_sl[i],y_l_sl[i],x_u_sl[i],y_u_sl[i],(canvas_mu[i]).c_str());
        TLegend* leg1 = new TLegend(0.65,0.70,0.85,0.85);
        leg1->SetFillColor(kWhite);
        leg1->SetFillStyle(1001);
        leg1->AddEntry(ttjets_pp_mu[i],"MC : ttH x 600","L");
        //leg1->AddEntry(ttjets_ol_mu[i],"MC : Background","L");
        //leg1->AddEntry(data_mu[i],"Data","L");
        ttjets_pp_mu[i]->SetLineColor(kRed);
        ttjets_pp_mu[i]->SetLineWidth(2);
        ttjets_pp_mu[i]->Draw("same");
        /*
        ttjets_ol_mu[i]->SetLineColor(kGreen+3);
        ttjets_ol_mu[i]->SetLineWidth(3);
        ttjets_ol_mu[i]->Draw("same");
        data_mu[i]->SetLineColor(kBlue);
        data_mu[i]->SetLineWidth(2);
        data_mu[i]->Draw("same");
        */
        leg1->Draw("same");
        c1->Print((out_file_mu[i]).c_str());
        delete c1;
        delete leg1;
        
    }

    // Double Electron
    for(int i=0; i<n_ee; i++){

        TCanvas *c1 = new TCanvas("c1","test",1100,650);
        c1->DrawFrame(x_l_di[i],y_l_di[i],x_u_di[i],y_u_di[i],(canvas_ee[i]).c_str());
        TLegend* leg1 = new TLegend(0.65,0.70,0.85,0.85);
        leg1->SetFillColor(kWhite);
        leg1->SetFillStyle(1001);
        leg1->AddEntry(ttjets_pp_ee[i],"MC : ttH x 600","L");
        //leg1->AddEntry(ttjets_ol_ee[i],"MC : Background","L");
        //leg1->AddEntry(data_ee[i],"Data","L");
        ttjets_pp_ee[i]->SetLineColor(kRed);
        ttjets_pp_ee[i]->SetLineWidth(2);
        ttjets_pp_ee[i]->Draw("same");
        /*
        ttjets_ol_ee[i]->SetLineColor(kGreen+3);
        ttjets_ol_ee[i]->SetLineWidth(3);
        ttjets_ol_ee[i]->Draw("same");
        data_ee[i]->SetLineColor(kBlue);
        data_ee[i]->SetLineWidth(2);
        data_ee[i]->Draw("same");
        */
        leg1->Draw("same");
        c1->Print((out_file_ee[i]).c_str());
        delete c1;
        delete leg1;
        
    }

    // Electron Muon
    for(int i=0; i<n_emu; i++){

        TCanvas *c1 = new TCanvas("c1","test",1100,650);
        c1->DrawFrame(x_l_di[i],y_l_di[i],x_u_di[i],y_u_di[i],(canvas_emu[i]).c_str());
        TLegend* leg1 = new TLegend(0.65,0.70,0.85,0.85);
        leg1->SetFillColor(kWhite);
        leg1->SetFillStyle(1001);
        leg1->AddEntry(ttjets_pp_emu[i],"MC : ttH x 600","L");
        //leg1->AddEntry(ttjets_ol_emu[i],"MC : Background","L");
        //leg1->AddEntry(data_emu[i],"Data","L");
        ttjets_pp_emu[i]->SetLineColor(kRed);
        ttjets_pp_emu[i]->SetLineWidth(2);
        ttjets_pp_emu[i]->Draw("same");
        /*
        ttjets_ol_emu[i]->SetLineColor(kGreen+3);
        ttjets_ol_emu[i]->SetLineWidth(3);
        ttjets_ol_emu[i]->Draw("same");
        data_emu[i]->SetLineColor(kBlue);
        data_emu[i]->SetLineWidth(2);
        data_emu[i]->Draw("same");
        */
        leg1->Draw("same");
        c1->Print((out_file_emu[i]).c_str());
        delete c1;
        delete leg1;
        
    }

    // Double Muon
    for(int i=0; i<n_mumu; i++){

        TCanvas *c1 = new TCanvas("c1","test",1100,650);
        c1->DrawFrame(x_l_di[i],y_l_di[i],x_u_di[i],y_u_di[i],(canvas_mumu[i]).c_str());
        TLegend* leg1 = new TLegend(0.65,0.70,0.85,0.85);
        leg1->SetFillColor(kWhite);
        leg1->SetFillStyle(1001);
        leg1->AddEntry(ttjets_pp_mumu[i],"MC : ttH x 600","L");
        //leg1->AddEntry(ttjets_ol_mumu[i],"MC : Background","L");
        //leg1->AddEntry(data_mumu[i],"Data","L");
        ttjets_pp_mumu[i]->SetLineColor(kRed);
        ttjets_pp_mumu[i]->SetLineWidth(2);
        ttjets_pp_mumu[i]->Draw("same");
        /*
        ttjets_ol_mumu[i]->SetLineColor(kGreen+3);
        ttjets_ol_mumu[i]->SetLineWidth(3);
        ttjets_ol_mumu[i]->Draw("same");
        data_mumu[i]->SetLineColor(kBlue);
        data_mumu[i]->SetLineWidth(2);
        data_mumu[i]->Draw("same");
        */
        leg1->Draw("same");
        c1->Print((out_file_mumu[i]).c_str());
        delete c1;
        delete leg1;
        
    }

	
    f1->Close();
    //f2->Close();
    //f3->Close();
	delete f1;
    //delete f2;
    //delete f3;

	return;
}


