//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Feb 18 15:14:21 2020 by ROOT version 6.14/09
// from TTree efficiencyTree/Reco Matched Jet Tree 
// found on file: l1TNtuple-ZeroBias_test.root
//////////////////////////////////////////////////////////

#ifndef controlplot_h
#define controlplot_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <string>
#include <iostream>
#include <fstream>
#include <TMap.h>
#include <vector>
//#ifdef __MAKECINT__
//#pragma link C++ class vector<std::vector<int> >+;
//#pragma link C++ class vector<std::vector<std::string> >+;
//#pragma link C++ class vector<std::vector<float> >+;
//#pragma link C++ class vector<std::vector<bool> >+;
//#pragma extra_include "vector";
//#endif

using namespace std;

// Header file for the classes stored in the TTree if any.

class controlplot {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
 
   TFile *fileName;
   TH1F *recojetpt_eff_den1, *recojetpt_eff_num1, *recojeteta_eff_den1, *recojeteta_eff_num1, *recojetphi_eff_den1, *recojetphi_eff_num1, *l1jetpt, *l1jetpt_rate1;

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           lumi;
   Int_t           event;
   Double_t        genDR;
   Int_t           genId;
   Double_t        recoPt_1;
   Double_t        recoEta_1;
   Double_t        recoPhi_1;
   Int_t           recoNthJet_1;
   Double_t        recoPt_2;
   Double_t        recoEta_2;
   Double_t        recoPhi_2;
   Int_t           recoNthJet_2;
   Double_t        recoDeltaEta;
   Double_t        recoDeltaPhi;
   Double_t        recoDeltaR;
   Double_t        recoMass;
   Double_t        l1Pt_1;
   Double_t        l1Eta_1;
   Double_t        l1Phi_1;
   Double_t        jetClusterPt;
   Double_t        jetClusterEta;
   Double_t        jetClusterPhi;
   Int_t           l1NthJet_1;
   Int_t           l1NTau_1;
   Double_t        l1Pt_2;
   Double_t        l1Eta_2;
   Double_t        l1Phi_2;
   Int_t           l1NthJet_2;
   Int_t           l1NTau_2;
   Double_t        l1DeltaEta;
   Double_t        l1DeltaPhi;
   Double_t        l1DeltaR;
   Double_t        l1Mass;
   Int_t           l1Matched_1;
   Int_t           l1Matched_2;
   Int_t           nRecoJets;
   Int_t           nL1Jets;
   Double_t        vbfBDT;
   vector<float>   *tau1;
   vector<float>   *tau2;
   vector<float>   *tau3;
   vector<int>     *nSubJets;
   vector<vector<int> > *subJetHFlav;
   vector<int>     *nBHadrons;
   vector<int>     *nL1Taus;
   vector<TLorentzVector> *allRegions;
   vector<TLorentzVector> *hcalTPGs;
   vector<TLorentzVector> *ecalTPGs;
   vector<TLorentzVector> *caloClusters;
   vector<TLorentzVector> *l1Jets;
   vector<TLorentzVector> *ak8Jets;
   vector<TLorentzVector> *subJets;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_event;   //!
   TBranch        *b_genDR;   //!
   TBranch        *b_genId;   //!
   TBranch        *b_recoPt_1;   //!
   TBranch        *b_recoEta_1;   //!
   TBranch        *b_recoPhi_1;   //!
   TBranch        *b_recoNthJet_1;   //!
   TBranch        *b_recoPt_2;   //!
   TBranch        *b_recoEta_2;   //!
   TBranch        *b_recoPhi_2;   //!
   TBranch        *b_recoNthJet_2;   //!
   TBranch        *b_recoDeltaEta;   //!
   TBranch        *b_recoDeltaPhi;   //!
   TBranch        *b_recoDeltaR;   //!
   TBranch        *b_recoMass;   //!
   TBranch        *b_l1Pt_1;   //!
   TBranch        *b_l1Eta_1;   //!
   TBranch        *b_l1Phi_1;
   TBranch       *jetClusterPt;
  TBranch        *jetClusterEta;
  TBranch        *jetClusterPhi;//!
   TBranch        *b_l1NthJet_1;   //!
   TBranch        *b_l1NTau_1;   //!
   TBranch        *b_l1Pt_2;   //!
   TBranch        *b_l1Eta_2;   //!
   TBranch        *b_l1Phi_2;   //!
   TBranch        *b_l1NthJet_2;   //!
   TBranch        *b_l1NTau_2;   //!
   TBranch        *b_l1DeltaEta;   //!
   TBranch        *b_l1DeltaPhi;   //!
   TBranch        *b_l1DeltaR;   //!
   TBranch        *b_l1Mass;   //!
   TBranch        *b_l1Matched_1;   //!
   TBranch        *b_l1Matched_2;   //!
   TBranch        *b_nRecoJets;   //!
   TBranch        *b_nL1Jets;   //!
   TBranch        *b_vbfBDT;   //!
   TBranch        *b_tau1;   //!
   TBranch        *b_tau2;   //!
   TBranch        *b_tau3;   //!
   TBranch        *b_nSubJets;   //!
   TBranch        *b_subJetHFlav;   //!
   TBranch        *b_nBHadrons;   //!
   TBranch        *b_nL1Taus;   //!
   TBranch        *b_allRegions;   //!
   TBranch        *b_hcalTPGs;   //!
   TBranch        *b_ecalTPGs;   //!
   TBranch        *b_caloClusters;   //!
   TBranch        *b_l1Jets;   //!
   TBranch        *b_ak8Jets;   //!
   TBranch        *b_subJets;   //!

   controlplot(const char* file1, const char* file2, const char* recoeta, const char* l1pt);
   virtual ~controlplot();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(const char* recoeta, const char* l1pt);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     BookHistos(const char* file2);
};

#endif

#ifdef controlplot_cxx
controlplot::controlplot(const char* file1, const char* file2, const char* recoeta, const char* l1pt)
{
   BookHistos(file2);
   TChain *chain = new TChain("l1NtupleProducer/efficiencyTree");
   ifstream file;
   file.open(file1, ifstream::in );
   char filename[2000];
   while (true) {
      file >> filename;
      if( file.eof() ) break;
         chain->Add(filename);
         cout<<"Added "<<filename<<endl;
   }//loop over while

   //BookHistos(file2);
   //TChain *chain = new TChain("l1NtupleProducer/efficiencyTree");
   //chain->Add(file1);

   Init(chain);
}

controlplot::~controlplot()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   fileName->cd();
   fileName->Write();
   fileName->Close();
}

Int_t controlplot::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t controlplot::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void controlplot::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   tau1 = 0;
   tau2 = 0;
   tau3 = 0;
   nSubJets = 0;
   subJetHFlav = 0;
   nBHadrons = 0;
   nL1Taus = 0;
   allRegions = 0;
   hcalTPGs = 0;
   ecalTPGs = 0;
   caloClusters = 0;
   l1Jets = 0;
   ak8Jets = 0;
   subJets = 0;

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("genDR", &genDR, &b_genDR);
   fChain->SetBranchAddress("genId", &genId, &b_genId);
   fChain->SetBranchAddress("recoPt_1", &recoPt_1, &b_recoPt_1);
   fChain->SetBranchAddress("recoEta_1", &recoEta_1, &b_recoEta_1);
   fChain->SetBranchAddress("recoPhi_1", &recoPhi_1, &b_recoPhi_1);
   fChain->SetBranchAddress("recoNthJet_1", &recoNthJet_1, &b_recoNthJet_1);
   fChain->SetBranchAddress("recoPt_2", &recoPt_2, &b_recoPt_2);
   fChain->SetBranchAddress("recoEta_2", &recoEta_2, &b_recoEta_2);
   fChain->SetBranchAddress("recoPhi_2", &recoPhi_2, &b_recoPhi_2);
   fChain->SetBranchAddress("recoNthJet_2", &recoNthJet_2, &b_recoNthJet_2);
   fChain->SetBranchAddress("recoDeltaEta", &recoDeltaEta, &b_recoDeltaEta);
   fChain->SetBranchAddress("recoDeltaPhi", &recoDeltaPhi, &b_recoDeltaPhi);
   fChain->SetBranchAddress("recoDeltaR", &recoDeltaR, &b_recoDeltaR);
   fChain->SetBranchAddress("recoMass", &recoMass, &b_recoMass);
   fChain->SetBranchAddress("l1Pt_1", &l1Pt_1, &b_l1Pt_1);
   fChain->SetBranchAddress("l1Eta_1", &l1Eta_1, &b_l1Eta_1);
   fChain->SetBranchAddress("l1Phi_1", &l1Phi_1, &b_l1Phi_1);
   fChain->SetBranchAddress("l1NthJet_1", &l1NthJet_1, &b_l1NthJet_1);
   fChain->SetBranchAddress("l1NTau_1", &l1NTau_1, &b_l1NTau_1);
   fChain->SetBranchAddress("l1Pt_2", &l1Pt_2, &b_l1Pt_2);
   fChain->SetBranchAddress("l1Eta_2", &l1Eta_2, &b_l1Eta_2);
   fChain->SetBranchAddress("l1Phi_2", &l1Phi_2, &b_l1Phi_2);
   fChain->SetBranchAddress("l1NthJet_2", &l1NthJet_2, &b_l1NthJet_2);
   fChain->SetBranchAddress("l1NTau_2", &l1NTau_2, &b_l1NTau_2);
   fChain->SetBranchAddress("l1DeltaEta", &l1DeltaEta, &b_l1DeltaEta);
   fChain->SetBranchAddress("l1DeltaPhi", &l1DeltaPhi, &b_l1DeltaPhi);
   fChain->SetBranchAddress("l1DeltaR", &l1DeltaR, &b_l1DeltaR);
   fChain->SetBranchAddress("l1Mass", &l1Mass, &b_l1Mass);
   fChain->SetBranchAddress("l1Matched_1", &l1Matched_1, &b_l1Matched_1);
   fChain->SetBranchAddress("l1Matched_2", &l1Matched_2, &b_l1Matched_2);
   fChain->SetBranchAddress("nRecoJets", &nRecoJets, &b_nRecoJets);
   fChain->SetBranchAddress("nL1Jets", &nL1Jets, &b_nL1Jets);
   fChain->SetBranchAddress("vbfBDT", &vbfBDT, &b_vbfBDT);
   fChain->SetBranchAddress("tau1", &tau1, &b_tau1);
   fChain->SetBranchAddress("tau2", &tau2, &b_tau2);
   fChain->SetBranchAddress("tau3", &tau3, &b_tau3);
   fChain->SetBranchAddress("nSubJets", &nSubJets, &b_nSubJets);
   fChain->SetBranchAddress("subJetHFlav", &subJetHFlav, &b_subJetHFlav);
   fChain->SetBranchAddress("nBHadrons", &nBHadrons, &b_nBHadrons);
   fChain->SetBranchAddress("nL1Taus", &nL1Taus, &b_nL1Taus);
   fChain->SetBranchAddress("allRegions", &allRegions, &b_allRegions);
   fChain->SetBranchAddress("hcalTPGs", &hcalTPGs, &b_hcalTPGs);
   fChain->SetBranchAddress("ecalTPGs", &ecalTPGs, &b_ecalTPGs);
   fChain->SetBranchAddress("caloClusters", &caloClusters, &b_caloClusters);
   fChain->SetBranchAddress("l1Jets", &l1Jets, &b_l1Jets);
   fChain->SetBranchAddress("ak8Jets", &ak8Jets, &b_ak8Jets);
   fChain->SetBranchAddress("subJets", &subJets, &b_subJets);

   Notify();
}

Bool_t controlplot::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void controlplot::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t controlplot::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef controlplot_cxx
