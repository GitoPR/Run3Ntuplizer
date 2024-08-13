#define controlplot_cxx
#include "controlplot.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

int main(int argc, char *argv[])
{

  if(argc > 1)
    {
      controlplot t(argv[1], argv[2], argv[3], argv[4]);
      t.Loop(argv[3], argv[4]);
    }
  return 0;
}

using namespace std;

void controlplot::Loop(const char* recoeta, const char* l1pt)
{
//   In a ROOT session, you can do:
//      root> .L controlplot.C
//      root> controlplot t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   std::stringstream sstrm1(recoeta);
   std::string recoetarange;
   sstrm1 >> recoetarange;
   bool recoetacut;
   std::stringstream sstrm2(l1pt);
   double l1ptcut;
   sstrm2 >> l1ptcut;
   double SF = 1.0; // naive scale factor
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
	Long64_t ientry = LoadTree(jentry);
	if (ientry < 0) break;
	nb = fChain->GetEntry(jentry);   nbytes += nb;

	//efficiencies
	jetClusterPt = SF*jetClusterPt; 
	recoetacut = false;
	if(recoetarange == "barrel") recoetacut = (abs(recoEta_1) <= 1.474);
	if(recoetarange == "endcap") recoetacut = (abs(recoEta_1) > 1.474 && abs(recoEta_1) <= 3);
	if(recoetarange == "all") recoetacut = true;
        if(recoPt_1 >= 200. && genId == 25 && genDR < 0.4) {
		recojetpt_eff_den1->Fill(recoPt_1); 
		recojeteta_eff_den1->Fill(recoEta_1); 
		recojetphi_eff_den1->Fill(recoPhi_1); 
		if (jetClusterPt >= l1ptcut) {
			recojetpt_eff_num1->Fill(recoPt_1);  
			recojeteta_eff_num1->Fill(recoEta_1); 
			recojetphi_eff_num1->Fill(recoPhi_1); 
		} 
	}

	//rates
	l1jetpt->Fill(jetClusterPt*SF);
	if(l1Jets->size() > 0) {
		double fillpt = -99;
		for(size_t i = 0; i < l1Jets->size(); i++){
			if(l1Jets->at(i).Pt() > fillpt) fillpt = l1Jets->at(i).Pt();
		}
		if(fillpt > -99) l1jetpt_rate1->Fill(fillpt*SF);
	}	
   }
}

void controlplot::BookHistos(const char* file2){
        fileName = new TFile(file2, "RECREATE");
        fileName->cd();
        char name[100];

        sprintf(name, "l1jetpt");
        l1jetpt = new TH1F (name,"l1jetpt", 100, 0, 1100);
        l1jetpt->GetXaxis()->SetTitle("p_{T} [GeV]");

        sprintf(name, "l1jetpt_rate1");
        l1jetpt_rate1 = new TH1F (name,"l1jetpt_rate1", 100, 0, 1100);
        l1jetpt_rate1->GetXaxis()->SetTitle("p_{T} [GeV]");

        sprintf(name, "recojetpt_eff_den1");
        recojetpt_eff_den1 = new TH1F (name,"recojetpt_eff_den1", 40, 200, 500);
        recojetpt_eff_den1->GetXaxis()->SetTitle("p_{T} [GeV]");

        sprintf(name, "recojetpt_eff_num1");
        recojetpt_eff_num1 = new TH1F (name,"recojetpt_eff_num1", 40, 200, 500);
        recojetpt_eff_num1->GetXaxis()->SetTitle("p_{T} [GeV]");

        sprintf(name, "recojeteta_eff_den1");
        recojeteta_eff_den1 = new TH1F (name,"recojeteta_eff_den1", 40, -5, 5);
        recojeteta_eff_den1->GetXaxis()->SetTitle("#eta");

        sprintf(name, "recojeteta_eff_num1");
        recojeteta_eff_num1 = new TH1F (name,"recojeteta_eff_num1", 40, -5, 5);
        recojeteta_eff_num1->GetXaxis()->SetTitle("#eta");

        sprintf(name, "recojetphi_eff_den1");
        recojetphi_eff_den1 = new TH1F (name,"recojetphi_eff_den1", 40, -M_PI, M_PI);
        recojetphi_eff_den1->GetXaxis()->SetTitle("#phi");

        sprintf(name, "recojetphi_eff_num1");
        recojetphi_eff_num1 = new TH1F (name,"recojetphi_eff_num1", 40, -M_PI, M_PI);
        recojetphi_eff_num1->GetXaxis()->SetTitle("#phi");
}
