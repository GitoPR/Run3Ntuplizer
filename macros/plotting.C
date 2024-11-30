#include "TH2F.h"
#include "TEllipse.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
#include "ROOT/RDataFrame.hxx"

void plot(){
  TFile* f = TFile::Open("gghbb_genmatching.root");
  ROOT::RDataFrame df("l1NtupleProducer/efficiencyTree",f);

  auto pt_resolution = [](double recoPt_1 , double jetClusterPt)  {


    if(recoPt_1 > 0 && jetClusterPt > 0) { 
      return (recoPt_1 - jetClusterPt)/recoPt_1;
    }

    return -99.0;

  };

  auto eta_resolution = [](double recoEta_1 , double jetClusterEta)  {


    if(recoEta_1 > -99 && jetClusterEta > -99) {
      return (recoEta_1 - jetClusterEta)/recoEta_1;
    }

    return -99.0;

  };

  auto phi_resolution = [](double recoPhi_1 , double jetClusterPhi)  {


    if(recoPhi_1 > -99 && jetClusterPhi> -99) {
      return (recoPhi_1 - jetClusterPhi)/recoPhi_1;
    }

    return -99.0;

  };
 
 

    TCanvas *c1 = new TCanvas();
    gStyle->SetOptStat(0); // Turn off statistics box                                                                                                                                                       
    int nbinsx = 50;
    float xmax = 10;
    float xmin = -10;

    auto reso  = df.Define("resolution", pt_resolution, {"recoPt_1", "jetClusterPt"}).Histo1D({"resolution", "Jet Clustering Pt_Resolution", nbinsx, xmin, xmax}, "resolution");

    reso->Draw();
    c1->SaveAs("clustering_pt_resolution.png");

    TCanvas *c2 = new TCanvas ();
    gStyle->SetOptStat(0);

    auto eta_reso  = df.Define("Eta Resolution", eta_resolution, {"recoEta_1", "jetClusterEta"}).Histo1D({"Eta Resolution", "Jet Clustering Eta Resolution", nbinsx, xmin, xmax}, "Eta Resolution");

    eta_reso->Draw();
    c2->SaveAs("clustering_eta_resolution.png");


    TCanvas *c3 = new TCanvas ();
    gStyle->SetOptStat(0);


    auto phi_reso  = df.Define("Phi Resolution", eta_resolution, {"recoPhi_1", "jetClusterPhi"}).Histo1D({"Phi Resolution", "Jet Clustering Phi Resolution", nbinsx, xmin, xmax}, "Phi Resolution");

    phi_reso->Draw();
    c3->SaveAs("clustering_phi_resolution.png");

    
    
  //resolution for eta


  //resolution for phi

  
}
