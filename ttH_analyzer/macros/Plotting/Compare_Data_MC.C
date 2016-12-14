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


void Compare_Data_MC() {
	
	double norm_tth, norm_ttjets;
	
	TFile *f1 = new TFile("Control_Plots_tth.root");
	TH1D *tth_h1 = (TH1D*)f1->Get("pt_e0");	
	double tth_h1_n = tth_h1->GetEntries(); 
	
	TFile *f2 = new TFile("Control_Plots_ttjets.root");
	TH1D *ttjets_h1 = (TH1D*)f2->Get("pt_e0");
	double ttjets_h1_n = ttjets_h1->GetEntries(); 
	
	TFile *f3 = new TFile("Control_Plots_data_SE.root");
	TH1D *data_SE_h1 = (TH1D*)f3->Get("pt_e0");
	double data_SE_h1_n = data_SE_h1->GetEntries(); 
	
	TFile *f4 = new TFile("Control_Plots_data_SMu.root");
	TH1D *data_SMu_h1 = (TH1D*)f4->Get("pt_mu0");
	double data_SMu_h1_n = data_SMu_h1->GetEntries(); 
	
	TH1D* mc_h1 = new TH1D("pt_tth_ttjets","Leading Electron pT Distribution ;pT (GeV);No. of Events",150,0,300);
	norm_tth = data_h1_n/tth_h1_n;
	norm_ttjets = data_h1_n/ttjets_h1_n;
	mc_h1->Add(tth_h1,ttjets_h1,norm_tth,norm_ttjets);   // weighting can be done here
	
	
	// Leading electron pt 
    
	TCanvas *c1 = new TCanvas("c1","test",1100,650);
	c1->DrawFrame(0,0,330,15,"Leading Electron pT Distribution ;pT (GeV);No. of Events");
	TLegend* leg1 = new TLegend(0.65,0.65,0.85,0.85);
    leg1->SetFillColor(kWhite);
    leg1->SetFillStyle(1001);
    leg1->AddEntry(mc_h1,"MC","L");
    leg1->AddEntry(data_SE_h1,"Data","L");	
	mc_h1->SetLineColor(kRed);
	mc_h1->SetLineWidth(2);
	mc_h1->Draw("b same");
	data_SE_h1->SetLineColor(kBlue);
	data_SE_h1->SetLineWidth(2);
	data_SE_h1->Draw("same");
	leg1->Draw("same");
	c1->Print("ele0_pt_SE.png");
	delete c1;
	delete leg1;
	
	delete f1;
	delete f2;
	delete f3;
	delete f4;
	
	return;
}


