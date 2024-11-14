
#include "TH2F.h"
#include "TEllipse.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
#include "ROOT/RDataFrame.hxx"

long long int entry = 0;
long long int kEntry = 0;
int nEtaBins = 14 ;


namespace {
  constexpr std::array<double, 15> fillEtaValues() {
    std::array<double, 15> EtaValues = {{0}};
    
    EtaValues[0]= -2.5; //ieta 4
    EtaValues[1]= -2.172; //ieta 5
    EtaValues[2]= -1.740;  //ieta 6
    EtaValues[3]= -1.392; //ieta 7
    EtaValues[4]= -1.044;
    EtaValues[5]= -0.696;
    EtaValues[6]= -0.348;
    EtaValues[7]= 0.001;
    EtaValues[8]= 0.348; //eta 10
    EtaValues[9]= 0.696;  //ieta 11
    EtaValues[10]=1.044;
    EtaValues[11]= 1.392;
    EtaValues[12]= 1.74;
    EtaValues[13]=2.172;
    EtaValues[14]=2.5;
 
    return EtaValues;
  }
  constexpr std::array<double, 15> EtaValues = fillEtaValues();
}


// namespace
//std::vector<double> createArray(double start, double end, int numPoints) {
//  std::vector<double> array(numPoints);
//  double step = (end - start) / (numPoints - 1);
//  for(int i = 0; i < numPoints; ++i) {
//      array[i] = start + i * step;
//  }
//  return array;
//
//  double start = -2.4;
//  double end = 2.4 ;
//  numPoints = 14

//phi values
    
void plotEventDisplay2(long long int _kEntry) {
  TFile* f = TFile::Open("../l1TNtuple-ggHBB_10.root");
  kEntry = _kEntry;
  auto makeHisto = [](ROOT::RVec<UShort_t> arr) {
    entry++;
    int cols = 14;
    int rows = 18;
    

    TH2F *h = new TH2F(Form("myHisto_%lld",entry),Form("Histo #%lld",entry), nEtaBins ,EtaValues.data(),rows,0,17);
    
    if(entry == kEntry){
      std::cout << arr << std::endl;
      std::cout << arr.size() << std::endl;
    }
    for (int j = 0; j < rows; j++){
      for(int i = 0; i < cols; i++) {
    if(entry == kEntry) {
      std::cout << arr[i + j * cols] << " ";
    }
    h->SetBinContent(i+1,j+1,arr[i + j * cols]);
      }
      if(entry == kEntry) std::cout << std::endl;
    }
    if (entry == kEntry)
      h->Draw();
  };
  
  ROOT::RDataFrame df("l1NtupleProducer/efficiencyTree",f,{"cregions"});
  df.Foreach(makeHisto,{"cregions"});
}

void plotEventDisplay2() {
  TFile* f = TFile::Open("../l1TNtuple-ggHBB-maxEtfinder-test.root");
  ROOT::RDataFrame df("l1NtupleProducer/efficiencyTree",f);
  
  auto makeHisto = [](int genId , double genEta_1, double  genPhi_1,ROOT::RVec<UShort_t> arr, ROOT::RVec<double> dgt_phi, ROOT::RVec<double> dgt_et, ROOT::RVec<int> dgt_id, ROOT::RVec<double> dgt_eta, ROOT::RVec<double> clusterCord) {
    entry++;
    int cols = 14;
    int rows = 18;

    // Where does the coordinate for the cluster algo come from

    //What is the coordinate with the most et. 

  // +++++++++++ for testing the results of maximum et finder  +++++++++++++++++++   
    /*    int size_arr = 252;

    double maxValue = arr[0];
    int index = 0;

    for (int i = 1; i < size_arr; ++i) {
        if (arr[i] > maxValue) {
            maxValue = arr[i];
            index = i;
        }
    }

    int test_ieta = index % 14;
    int test_iphi =static_cast<int>( index / 14);
    std::cout << "test max Et: " << maxValue << std::endl;
    std::cout << "test et Index: " << index << std::endl;
    std::cout << "test ieta: " << test_ieta << std::endl; 
    std::cout << "test  Iphi: " << test_iphi << "\n" << std::endl; */ 

    
    TCanvas *canvas = new TCanvas(Form("canvas_%lld", entry), Form("Canvas for Event %lld", entry), 800, 600);
    gStyle->SetOptStat(0); // Turn off statistics box
    TH2F *h = new TH2F(Form("myHisto_%lld",entry),Form("Histo #%lld",entry),nEtaBins,EtaValues.data() ,rows,-3.4, 3.4);
    for (int j = 0; j < rows; j++){
      for(int i = 0; i < cols; i++) {
    h->SetBinContent(i+1,j+1,arr[i + j * cols]);
    //std::cout << "this is the x axis: " << i+ 1 <<  "\t" <<  "this is the y : " << j+1 << "\t" << "regiod id "std::endl;  
      }
    }
    h->Draw("COLZ"); // Use "COLZ" to use color for bin content

    //Draw a circle at the position of a Higgs.

    
      TEllipse *circ = new TEllipse(genEta_1,genPhi_1,.1 , .1);
      circ->SetFillStyle(1001);
      circ->SetFillColor(kRed);
      circ->SetLineColor(kRed);
      circ->Draw("SAME");

      
      
      auto size = dgt_eta.size(); 
      for (size_t i = 0; i < size ; i++){
	/*std::cout << "dau_eta: " << dgt_eta[i] << std::endl;
	std::cout << "dau_phi: " << dgt_phi[i] << std::endl;
	std::cout << "dau_et: " << dgt_et[i] << std::endl;
	std::cout << "dau_id: " << dgt_id[i] << std::endl;  */ 
	

	circ ->SetFillStyle(1001);
	circ-> SetFillColor(kOrange + 2);
	circ -> SetLineColor(kOrange +2);
	circ -> Draw("SAME");
	
	
      }


      /*      
      double max_et_phi = 1.2;  
      double max_et_eta = 1.5; 
      TDiamond *diamond = new TDiamond (max_et_eta+2 , max_et_phi, max_et_eta, max_et_phi+2);
      diamond -> SetFillStyle(1001);
      diamond -> SetFillColor(kRed +2);
      diamond -> SetLineColor(kRed);
      diamond -> Draw("SAME");*/

      //Draw a square on th are of a jet region.

      std::cout << clusterCord[0] << "\t" << clusterCord[1] << "\t" << std::endl; 
      double seed_eta = clusterCord[0];
      double seed_phi = clusterCord[1];
      
      double right_corner_eta = clusterCord[2];
      double right_corner_phi = clusterCord[3];
      
      TBox *box = new TBox(seed_eta - 0.1 , seed_phi - 0.1, seed_eta + 0.1, seed_phi + 0.1) ;
      box -> SetFillStyle(1001);
      box -> SetFillColor(kRed+2);
      box->SetLineColor(kRed);
      box -> Draw("SAME"); 
      
    canvas->SaveAs(Form("histogram_event_%lld.png", entry)); // Save the canvas as an image file
  };
  //ROOT::RDataFrame df("l1NtupleProducer/efficiencyTree",f,{"cregions"});
  df.Foreach(makeHisto,{"genId" ,"genEta_1" , "genPhi_1", "cregions", "dgt_phi" , "dgt_et" , "dgt_id" , "dgt_eta","clusterCord"});
}
