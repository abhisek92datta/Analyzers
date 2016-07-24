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


void PU_data_hist_prod()
{
	TFile *f = new TFile("output.root");
	TH1D *h = (TH1D*)f->Get("pileup");

	for(int i=0; i<50; i++)
	{
		std::cout<<h->GetBinContent(i+1)<<"\n";
	}
	return;
}
