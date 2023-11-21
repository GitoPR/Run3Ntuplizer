// system include files
#include <memory>

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

#include "L1Trigger/L1TCaloLayer1/interface/UCTLayer1.hh"
#include "L1Trigger/L1TCaloLayer1/interface/UCTCrate.hh"
#include "L1Trigger/L1TCaloLayer1/interface/UCTCard.hh"
#include "L1Trigger/L1TCaloLayer1/interface/UCTRegion.hh"
#include "L1Trigger/L1TCaloLayer1/interface/UCTTower.hh"
#include "L1Trigger/L1TCaloLayer1/interface/UCTGeometry.hh"
#include "L1Trigger/L1TCaloLayer1/interface/UCTObject.hh"
#include "L1Trigger/L1TCaloLayer1/interface/UCTSummaryCard.hh"
#include "L1Trigger/L1TCaloLayer1/interface/UCTGeometryExtended.hh"
#include "L1Trigger/L1TCaloLayer1/interface/UCTParameters.hh"
//#include "L1Trigger/L1TCaloLayer1/interface/L1TCaloLayer1FetchLUTs.hh"

#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticleFwd.h"

#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegion.h"

#include "DataFormats/Math/interface/LorentzVector.h"

//#include "L1Trigger/Run3Ntuplizer/plugins/helpers.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

// GCT and RCT data formats
#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctCollections.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1Trigger/interface/Muon.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

using namespace l1t;
using namespace l1tcalo;
using namespace l1extra;
using namespace std;

bool compareByPt (l1extra::L1JetParticle i, l1extra::L1JetParticle j) { return(i.pt()>j.pt()); };
bool compareByObjectEt (UCTObject* i, UCTObject* j) { return(i->et() > j->et()); };

float towerEtaMap[28]= {
    0.0435,
    0.1305, 0.2175, 0.3045, 0.3915, 0.4785,
    0.5655, 0.6525, 0.7395, 0.8265, 0.9135,
    1.0005, 1.0875, 1.1745, 1.2615, 1.3485,
    1.4355, 1.5225, 1.6095, 1.6965, 1.7835,
    1.8705, 1.9575, 2.0445, 2.217, 2.391,
    2.565, //2.739,
    2.913,
    };

float towerPhiMap[72]=
    {-0.131, -0.044, 0.044, 0.131, 0.218, 0.305, 0.393, 0.480, 0.567, 0.654, 0.742, 0.829, 0.916, 1.004, 1.091, 1.178, 1.265, 1.353, 1.440, 1.527, 1.614, 1.702, 1.789, 1.876, 1.963, 2.051, 2.138, 2.225, 2.313, 2.400, 2.487, 2.574, 2.662, 2.749, 2.836, 2.923, 3.011, 3.098,
      -3.098, -3.011, -2.923, -2.836, -2.749, -2.662, -2.574, -2.487, -2.400, -2.313, -2.225, -2.138, -2.051, -1.963, -1.876, -1.789, -1.702, -1.614, -1.527, -1.440, -1.353, -1.265, -1.178, -1.091, -1.004, -0.916, -0.829, -0.742, -0.654, -0.567, -0.480, -0.393, -0.305, -0.218};

float getRecoEta(int ieta, short zside){
  float eta = -999;
  if(ieta<0 || ieta>(28*2)){
    std::cout<<"Error!!! towereta out of bounds in triggerGeometryTools.h "<<std::endl;
    std::cout<<"ieta "<<ieta<<std::endl;
    exit(0);
  }
  if(zside == 1)
    eta = towerEtaMap[ieta];
  else if(zside == -1)
    eta = towerEtaMap[ieta];
  else{
    std::cout<<"Error!!! zside out of bounds in triggerGeometryTools.h "<<std::endl;
    std::cout<<"zside "<<zside<<std::endl;
    exit(0);
  }
  return eta;
};

float getRecoEtaNew(int caloEta){
  float eta = -999.;
  static bool first = true;
  static double twrEtaValues[42];
  if(first) {
    twrEtaValues[0] = 0;
    for(unsigned int i = 0; i < 20; i++) {
      twrEtaValues[i + 1] = 0.0436 + i * 0.0872;
    }
    twrEtaValues[21] = 1.785;
    twrEtaValues[22] = 1.880;
    twrEtaValues[23] = 1.9865;
    twrEtaValues[24] = 2.1075;
    twrEtaValues[25] = 2.247;
    twrEtaValues[26] = 2.411;
    twrEtaValues[27] = 2.575;
    twrEtaValues[28] = 2.825;
    twrEtaValues[29] = 999.;
    twrEtaValues[30] = (3.15+2.98)/2.;
    twrEtaValues[31] = (3.33+3.15)/2.;
    twrEtaValues[32] = (3.50+3.33)/2.;
    twrEtaValues[33] = (3.68+3.50)/2.;
    twrEtaValues[34] = (3.68+3.85)/2.;
    twrEtaValues[35] = (3.85+4.03)/2.;
    twrEtaValues[36] = (4.03+4.20)/2.;
    twrEtaValues[37] = (4.20+4.38)/2.;
    twrEtaValues[38] = (4.74+4.38*3)/4.;
    twrEtaValues[39] = (4.38+4.74*3)/4.;
    twrEtaValues[40] = (5.21+4.74*3)/4.;
    twrEtaValues[41] = (4.74+5.21*3)/4.;
    first = false;
  }
  uint32_t absCaloEta = abs(caloEta);
  if(absCaloEta <= 41) {
    if(caloEta < 0)
      eta =  -twrEtaValues[absCaloEta];
    else
      eta = +twrEtaValues[absCaloEta];
  }
  return eta;
};

float getRecoPhi(int iphi){
  return towerPhiMap[iphi-1];
};

float getRecoPhiNew(int caloPhi){
  float phi = -999.;
  if(caloPhi > 72) phi = +999.;
  uint32_t absCaloPhi = std::abs(caloPhi) - 1;
  if(absCaloPhi < 36)
    phi = (((double) absCaloPhi + 0.5) * 0.0872);
  else
    phi = (-(71.5 - (double) absCaloPhi) * 0.0872);
  return phi;
};


//
// class declaration
//

class BoostedJetStudies : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit BoostedJetStudies(const edm::ParameterSet&);
  ~BoostedJetStudies();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  void zeroOutAllVariables();

private:
  void analyze(const edm::Event& evt, const edm::EventSetup& es);      
  virtual void beginJob() override;
  virtual void endJob() override;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);

  // ----------member data ---------------------------
  edm::EDGetTokenT<EcalTrigPrimDigiCollection> ecalTPSource;
  edm::EDGetTokenT<HcalTrigPrimDigiCollection> hcalTPSource;

  edm::EDGetTokenT<vector<reco::CaloJet> > jetSrc_;
  edm::EDGetTokenT<vector<pat::Jet> > jetSrcAK8_;
  edm::EDGetTokenT<reco::GenParticleCollection> genSrc_;

  edm::EDGetTokenT<BXVector<l1t::Jet>> stage2JetToken_;
  edm::EDGetTokenT<BXVector<l1t::Tau>> stage2TauToken_;
  edm::EDGetTokenT<vector<l1extra::L1JetParticle>> l1BoostedToken_;
  edm::EDGetTokenT<L1CaloRegionCollection> regionToken_;
  //const L1TCaloLayer1FetchLUTsTokens lutsTokens;

  std::vector<std::array<std::array<std::array<uint32_t, nEtBins>, nCalSideBins>, nCalEtaBins> > ecalLUT;
  std::vector<std::array<std::array<std::array<uint32_t, nEtBins>, nCalSideBins>, nCalEtaBins> > hcalLUT;
  std::vector<std::array<std::array<uint32_t, nEtBins>, nHfEtaBins> > hfLUT;

  std::vector<unsigned int> ePhiMap;
  std::vector<unsigned int> hPhiMap;
  std::vector<unsigned int> hfPhiMap;

  TH1F* nEvents;

  int run, lumi, event;

  double genPt_1, genEta_1, genPhi_1, genM_1, genDR;
  int genId, genMother;
  double recoPt_1, recoEta_1, recoPhi_1;
  double l1Pt_1, l1Eta_1, l1Phi_1;
  double seedPt_1, seedEta_1, seedPhi_1;
  
  int l1NthJet_1;
  int recoNthJet_1;
  int seedNthJet_1;

  double recoPt_;
  std::vector<int> nSubJets, nBHadrons, HFlav;
  std::vector<std::vector<int>> subJetHFlav;
  std::vector<float> tau1, tau2, tau3;

  std::vector<TLorentzVector> *l1Jets  = new std::vector<TLorentzVector>;
  std::vector<TLorentzVector> *seed180  = new std::vector<TLorentzVector>;
  std::vector<TLorentzVector> *tauseed  = new std::vector<TLorentzVector>;
  std::vector<TLorentzVector> *ak8Jets  = new std::vector<TLorentzVector>;
  std::vector<TLorentzVector> *subJets  = new std::vector<TLorentzVector>;

  uint32_t nPumBins;

  std::vector< std::vector< std::vector < uint32_t > > > pumLUT;

  std::vector< UCTTower* > twrList;

  double caloScaleFactor;

  uint32_t jetSeed;
  uint32_t tauSeed;
  double tauIsolationFactor;
  uint32_t eGammaSeed;
  double eGammaIsolationFactor;
  double boostedJetPtFactor;
  bool verbose;
  int fwVersion;
  double activityFraction12;

  //UCTSummaryCard *summaryCard;

  std::vector<TLorentzVector> *allRegions  = new std::vector<TLorentzVector>;
  std::vector<TLorentzVector> *allL1Jets = new std::vector<TLorentzVector>;
  std::vector<bool> allL1Signals;
  std::vector<float> allL1DR_withReco;
  std::vector<std::vector<int>> jetRegionEt;
  std::vector<std::vector<int>> jetRegionEGVeto;
  std::vector<std::vector<int>> jetRegionTauVeto;
  std::vector<string> regionEta, regionPhi;

  //std::vector<string> etaBits, phiBits, etaBits12, phiBits12, regionEta, regionPhi;

  void createBranches(TTree *tree);
  TTree* efficiencyTree;
  edm::Service<TFileService> tfs_;  

};

BoostedJetStudies::BoostedJetStudies(const edm::ParameterSet& iConfig) :
  jetSrc_(    consumes<vector<reco::CaloJet> >(iConfig.getParameter<edm::InputTag>("recoJets"))),
  jetSrcAK8_( consumes<vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("recoJetsAK8"))),
  genSrc_( consumes<reco::GenParticleCollection> (iConfig.getParameter<edm::InputTag>( "genParticles"))),
  stage2JetToken_(consumes<BXVector<l1t::Jet>>( edm::InputTag("caloStage2Digis","Jet","RECO"))),
  stage2TauToken_(consumes<BXVector<l1t::Tau>>( edm::InputTag("caloStage2Digis","Tau","RECO"))),
  l1BoostedToken_(consumes<vector<l1extra::L1JetParticle>>( edm::InputTag("uct2016EmulatorDigis","Boosted",""))),
  regionToken_(consumes<L1CaloRegionCollection>(edm::InputTag("simCaloStage2Layer1Digis"))),
  nPumBins(iConfig.getParameter<unsigned int>("nPumBins")),
  pumLUT(nPumBins, std::vector<std::vector<uint32_t>>(2, std::vector<uint32_t>(13))),
  caloScaleFactor(iConfig.getParameter<double>("caloScaleFactor")),
  jetSeed(iConfig.getParameter<unsigned int>("jetSeed")),
  tauSeed(iConfig.getParameter<unsigned int>("tauSeed")),
  tauIsolationFactor(iConfig.getParameter<double>("tauIsolationFactor")),
  eGammaSeed(iConfig.getParameter<unsigned int>("eGammaSeed")),
  eGammaIsolationFactor(iConfig.getParameter<double>("eGammaIsolationFactor")),
  boostedJetPtFactor(iConfig.getParameter<double>("boostedJetPtFactor")),
  verbose(iConfig.getParameter<bool>("verbose")),
  fwVersion(iConfig.getParameter<int>("firmwareVersion"))
{

  std::vector<double> pumLUTData;
  char pumLUTString[10];
  for (uint32_t pumBin = 0; pumBin < nPumBins; pumBin++) {
    for (uint32_t side = 0; side < 2; side++) {
      if (side == 0)
        sprintf(pumLUTString, "pumLUT%2.2dp", pumBin);
      else
        sprintf(pumLUTString, "pumLUT%2.2dn", pumBin);
      pumLUTData = iConfig.getParameter<std::vector<double>>(pumLUTString);
      for (uint32_t iEta = 0; iEta < std::max((uint32_t)pumLUTData.size(), MaxUCTRegionsEta); iEta++) {
        pumLUT[pumBin][side][iEta] = (uint32_t)round(pumLUTData[iEta] / caloScaleFactor);
      }
      if (pumLUTData.size() != (MaxUCTRegionsEta))
        edm::LogError("L1TCaloSummary") << "PUM LUT Data size integrity check failed; Expected size = "
                                        << MaxUCTRegionsEta << "; Provided size = " << pumLUTData.size()
                                        << "; Will use what is provided :(" << std::endl;
    }
  }

  //summaryCard = new UCTSummaryCard(&pumLUT, jetSeed, tauSeed, tauIsolationFactor, eGammaSeed, eGammaIsolationFactor);

  // Initialize the Tree

  recoPt_      = iConfig.getParameter<double>("recoPtCut");
  nEvents      = tfs_->make<TH1F>( "nEvents"  , "nEvents", 2,  0., 1. );
  efficiencyTree = tfs_->make<TTree>("efficiencyTree", "Gen Matched Jet Tree ");
  createBranches(efficiencyTree);
}

BoostedJetStudies::~BoostedJetStudies() {
  delete allRegions;
  delete allL1Jets;
  delete l1Jets;
  delete seed180;
  delete tauseed;
  delete ak8Jets;
  delete subJets;
  //if(summaryCard != 0) delete summaryCard;
}

//
// member functions
//

// ------------ method called to produce the data  ------------

void BoostedJetStudies::analyze( const edm::Event& evt, const edm::EventSetup& es )
{
  using namespace edm;

  nEvents->Fill(1);
  run = evt.id().run();
  lumi = evt.id().luminosityBlock();
  event = evt.id().event();
  Handle<L1CaloRegionCollection> regions;
   
  std::vector<reco::CaloJet> goodJets;
  std::vector<pat::Jet> goodJetsAK8;
  std::vector<l1t::Jet> seeds;

  l1Jets->clear();
  seed180->clear();
  tauseed->clear();
  ak8Jets->clear();
  subJets->clear();
  nSubJets.clear();
  nBHadrons.clear();
  subJetHFlav.clear();
  tau1.clear();
  tau2.clear();
  tau3.clear();

  allRegions->clear();
  allL1Jets->clear();
  allL1Signals.clear();
  allL1DR_withReco.clear();
  jetRegionEt.clear();
  jetRegionEGVeto.clear();
  jetRegionTauVeto.clear();
  regionEta.clear();
  regionPhi.clear();

  UCTGeometry g;
  //summaryCard->clearRegions();
  UCTSummaryCard summaryCard = UCTSummaryCard(&pumLUT, jetSeed, tauSeed, tauIsolationFactor, eGammaSeed, eGammaIsolationFactor);
  std::vector<UCTRegion*> inputRegions;
  inputRegions.clear();
  edm::Handle<std::vector<L1CaloRegion>> regionCollection;
  if (!evt.getByToken(regionToken_, regionCollection))
    edm::LogError("L1TCaloSummary") << "UCT: Failed to get regions from region collection!";
  evt.getByToken(regionToken_, regionCollection);
  for (const L1CaloRegion& i : *regionCollection) {  //for (auto i : rgnCollection) {
    UCTRegionIndex r = g.getUCTRegionIndexFromL1CaloRegion(i.gctEta(), i.gctPhi());
    UCTTowerIndex t = g.getUCTTowerIndexFromL1CaloRegion(r, i.raw());
    uint32_t absCaloEta = std::abs(t.first);
    uint32_t absCaloPhi = std::abs(t.second);
    bool negativeEta = false;
    if (t.first < 0)
      negativeEta = true;
    uint32_t crate = g.getCrate(t.first, t.second);
    uint32_t card = g.getCard(t.first, t.second);
    uint32_t region = g.getRegion(absCaloEta, absCaloPhi);
    UCTRegion* test = new UCTRegion(crate, card, negativeEta, region, fwVersion);
    //std::cout<<"processing...."<<std::endl;
    if(!test->process()) std::cout<<"failed to process region"<<std::endl;
    test->setRegionSummary(i.raw());
    inputRegions.push_back(test);

    uint16_t test_raw = i.raw();
    uint32_t test_et = i.et();
    //std::cout<<test_et<<std::endl;
    uint32_t test_rEta = i.id().ieta();
    uint32_t test_rPhi = i.id().iphi();
    UCTRegionIndex test_rIndex = g.getUCTRegionIndexFromL1CaloRegion(test_rEta, test_rPhi);
    UCTTowerIndex test_tIndex = g.getUCTTowerIndexFromL1CaloRegion(test_rIndex, test_raw);
    int test_cEta = test_tIndex.first;
    int test_cPhi = test_tIndex.second;
    float pt = test_et;
    float eta = getRecoEtaNew(test_cEta);;
    float phi = getRecoPhiNew(test_cPhi);
    TLorentzVector temp ;
    temp.SetPtEtaPhiE(pt,eta,phi,pt);
    allRegions->push_back(temp);
  }
  //summaryCard->setRegionData(inputRegions);
  summaryCard.setRegionData(inputRegions);

  //if (!summaryCard->process()) {
  if (!summaryCard.process()) {
    edm::LogError("L1TCaloSummary") << "UCT: Failed to process summary card" << std::endl;
    exit(1);
  }

  Handle<vector<reco::CaloJet> > jets;
  if(evt.getByToken(jetSrc_, jets)){//Begin Getting Reco Jets
    for (const reco::CaloJet &jet : *jets) {
      //std::cout<<jet.pt()<<std::endl;
      if(jet.pt() > recoPt_ ) {
        goodJets.push_back(jet);
      }
    }
  }
  else
    cout<<"Error getting calo jets"<<std::endl;

  Handle<vector<pat::Jet> > jetsAK8;

  if(evt.getByToken(jetSrcAK8_, jetsAK8)){//Begin Getting AK8 Jets
    for (const pat::Jet &jetAK8 : *jetsAK8) {
      if(jetAK8.pt() > recoPt_ ) {
        nSubJets.push_back(jetAK8.subjets("SoftDropPuppi").size());
        nBHadrons.push_back(jetAK8.jetFlavourInfo().getbHadrons().size());
        TLorentzVector temp ;
        temp.SetPtEtaPhiE(jetAK8.pt(),jetAK8.eta(),jetAK8.phi(),jetAK8.et());
        ak8Jets->push_back(temp);
        if(jetAK8.subjets("SoftDropPuppi").size() ==  2 && jetAK8.jetFlavourInfo().getbHadrons().size() > 1){
          goodJetsAK8.push_back(jetAK8);
        }
      }
    }
  }
  else
    cout<<"Error getting AK8 jets"<<std::endl;

  double pt = 0;
  double eta = -999.;
  double phi = -999.;
  double caloScaleFactor = 0.5;
  std::vector<int> region_ets;
  std::vector<int> region_eg;
  std::vector<int> region_tau;

  // Sort boosted objects by Et
  std::list<UCTObject*> boostedJetObjs = summaryCard.getBoostedJetObjs();
  vector<UCTObject*> BoostedSorted;
  for(std::list<UCTObject*>::const_iterator i = boostedJetObjs.begin(); i != boostedJetObjs.end(); i++) {
    BoostedSorted.push_back(*i);
  }
  std::sort(BoostedSorted.begin(), BoostedSorted.end(), compareByObjectEt);

  int count = 0; // only consider up to 20 jets per event, match boosted objects to both ak4 and ak8 jets within dR < 0.2 (0.4 selects too many boosted)

  vector<int> indices_ak4;
  indices_ak4.clear();
  count = 0;
  for(auto object : BoostedSorted){
    float jetDR = 0.2;
    int temp_index = -99;
    for(auto jet:goodJets){
      if(reco::deltaR(jet.eta(), jet.phi(), g.getUCTTowerEta(object->iEta()), g.getUCTTowerPhi(object->iPhi())) < 0.2){
        jetDR = reco::deltaR(jet.eta(), jet.phi(), g.getUCTTowerEta(object->iEta()), g.getUCTTowerPhi(object->iPhi()));
        temp_index = count;
      }
    }
    if(jetDR < 0.2) indices_ak4.push_back(temp_index);
    count++;
  }

  vector<int> indices_ak8;
  indices_ak8.clear();
  count = 0;
  for(auto object : BoostedSorted){
    float jetDR = 0.2;
    int temp_index = -99;
    for(auto jet:goodJetsAK8){
      if(reco::deltaR(jet.eta(), jet.phi(), g.getUCTTowerEta(object->iEta()), g.getUCTTowerPhi(object->iPhi())) < 0.2){
        jetDR = reco::deltaR(jet.eta(), jet.phi(), g.getUCTTowerEta(object->iEta()), g.getUCTTowerPhi(object->iPhi()));
        temp_index = count;
      } 
    }
    if(jetDR < 0.2) indices_ak8.push_back(temp_index);
    count++;
  }   

  count = -1;

  //for(std::list<UCTObject*>::const_iterator i = boostedJetObjs.begin(); i != boostedJetObjs.end(); i++) {
  for(auto object : BoostedSorted){
    count++;
    if(!( std::find(indices_ak4.begin(), indices_ak4.end(), count) != indices_ak4.end() || std::find(indices_ak8.begin(), indices_ak8.end(), count) != indices_ak8.end())) continue; // only consider matched boosted objects
    region_ets.clear();
    region_eg.clear();
    region_tau.clear();
    //const UCTObject* object = *i;
    pt = ((double) object->et()) * caloScaleFactor * boostedJetPtFactor;
    eta = g.getUCTTowerEta(object->iEta());
    phi = g.getUCTTowerPhi(object->iPhi());
    bitset<3> activeRegionEtaPattern = 0;
    for (uint32_t iEta = 0; iEta < 3; iEta++) {
      bool activeStrip = false;
      for (uint32_t iPhi = 0; iPhi < 3; iPhi++) {
        region_ets.push_back(object->boostedJetRegionET()[3 * iEta + iPhi]);
        region_eg.push_back(object->boostedJetRegionEGammaVeto()[3 * iEta + iPhi]);
        region_tau.push_back(object->boostedJetRegionTauVeto()[3 * iEta + iPhi]);
        if (object->boostedJetRegionET()[3 * iEta + iPhi] > 30 &&
            object->boostedJetRegionET()[3 * iEta + iPhi] > object->et() * 0.0625)
          activeStrip = true;
      }
      if (activeStrip)
        activeRegionEtaPattern |= (0x1 << iEta);
    }
    bitset<3> activeRegionPhiPattern = 0;
    for (uint32_t iPhi = 0; iPhi < 3; iPhi++) {
      bool activeStrip = false;
      for (uint32_t iEta = 0; iEta < 3; iEta++) {
        if (object->boostedJetRegionET()[3 * iEta + iPhi] > 30 &&
            object->boostedJetRegionET()[3 * iEta + iPhi] > object->et() * 0.0625)
          activeStrip = true;
      }
      if (activeStrip)
        activeRegionPhiPattern |= (0x1 << iPhi);
    }
    regionEta.push_back(activeRegionEtaPattern.to_string<char,std::string::traits_type,std::string::allocator_type>());
    regionPhi.push_back(activeRegionPhiPattern.to_string<char,std::string::traits_type,std::string::allocator_type>());

    //std::cout<<pt<<"\t"<<eta<<"\t"<<phi<<"\t"<<object->boostedJetRegionET()[0]<<"\t"<<object->boostedJetRegionET()[1]<<"\t"<<object->boostedJetRegionET()[2]<<"\t"<<object->boostedJetRegionET()[3]<<std::endl;
    TLorentzVector temp;
    temp.SetPtEtaPhiE(pt,eta,phi,pt);
    allL1Jets->push_back(temp);

    //jetRegionEt.push_back(object->boostedJetRegionET());
    //jetRegionEGVeto.push_back(object->boostedJetRegionEGammaVeto());
    //jetRegionTauVeto.push_back(object->boostedJetRegionTauVeto());
    //std::cout<<region_ets.size()<<"\t"<<region_eg.size()<<"\t"<<region_tau.size()<<std::endl;
    jetRegionEt.push_back(region_ets);
    jetRegionEGVeto.push_back(region_eg);
    jetRegionTauVeto.push_back(region_tau);
    //count++;
    if (count == 19) break;

    //bitset<12> eta_in = object->activeTowerEta(); 
    //bitset<12> phi_in = object->activeTowerPhi();
    //std::cout<<"eta_in: "<<eta_in<<"\t"<<"phi_in: "<<phi_in<<std::endl;
    //etaBits.push_back(eta_in.to_string<char,std::string::traits_type,std::string::allocator_type>());
    //phiBits.push_back(phi_in.to_string<char,std::string::traits_type,std::string::allocator_type>());

    //bool activeTower[12][12];
    //uint32_t activityLevel = object->et()*activityFraction12;
    //for(uint32_t iPhi = 0; iPhi < 12; iPhi++){
    //  for(uint32_t iEta = 0; iEta < 12; iEta++){
    //    uint32_t towerET = object->boostedJetTowers()[iEta*12+iPhi];
    //    //std::cout<<towerET<<std::endl;
    //    if(towerET > activityLevel) {
    //      activeTower[iEta][iPhi] = true;
    //    }
    //    else activeTower[iEta][iPhi] = false;
    //  }
    //}
    //bitset<12> activeTowerEtaPattern = 0;
    //for(uint32_t iEta = 0; iEta < 12; iEta++){
    //  bool activeStrip = false;
    //  for(uint32_t iPhi = 0; iPhi < 12; iPhi++){
    //    if(activeTower[iEta][iPhi]) activeStrip = true;
    //  }
    //  if(activeStrip) activeTowerEtaPattern |= (0x1 << iEta);
    //}
    //bitset<12> activeTowerPhiPattern = 0;
    //for(uint32_t iPhi = 0; iPhi < 12; iPhi++){
    //  bool activeStrip = false;
    //  for(uint32_t iEta = 0; iEta < 12; iEta++){
    //    if(activeTower[iEta][iPhi]) activeStrip = true;
    //  }
    //  if(activeStrip) activeTowerPhiPattern |= (0x1 << iPhi);
    //}
    //etaBits12.push_back(activeTowerEtaPattern.to_string<char,std::string::traits_type,std::string::allocator_type>());
    //phiBits12.push_back(activeTowerPhiPattern.to_string<char,std::string::traits_type,std::string::allocator_type>());
    //int nActiveRegion = 0;
    //for(int i = 0; i < 9; i++){
    //  if(object->boostedJetRegionET()[i] > object->et()*0.0625 ) nActiveRegion++;
    //}
   // std::cout<<"nActiveRegion: "<<nActiveRegion<<std::endl;
  }


  // Accessing existing L1 seed stored in MINIAOD
  edm::Handle<BXVector<l1t::Jet>> stage2Jets;
  if(!evt.getByToken(stage2JetToken_, stage2Jets)) cout<<"ERROR GETTING THE STAGE 2 JETS"<<std::endl;
  evt.getByToken(stage2JetToken_, stage2Jets);
  const BXVector<l1t::Jet> &s2j = *stage2Jets;
  for(auto obj : s2j) {
    seeds.push_back(obj);
    TLorentzVector temp;
    temp.SetPtEtaPhiE(obj.pt(), obj.eta(), obj.phi(), obj.pt());
    seed180->push_back(temp);
  }

  edm::Handle<BXVector<l1t::Tau>> stage2Taus;
  if(!evt.getByToken(stage2TauToken_, stage2Taus)) cout<<"ERROR GETTING THE STAGE 2 TAUS"<<std::endl;
  evt.getByToken(stage2TauToken_, stage2Taus);
  const BXVector<l1t::Tau> &s2t = *stage2Taus;
  for(auto obj : s2t) {
    TLorentzVector temp;
    temp.SetPtEtaPhiE(obj.pt(), obj.eta(), obj.phi(), obj.pt());
    tauseed->push_back(temp);
  }

  // Accessing L1boosted collection
  edm::Handle<vector<l1extra::L1JetParticle>> l1Boosted;
  if(!evt.getByToken(l1BoostedToken_, l1Boosted)) cout<<"ERROR GETTING THE L1BOOSTED JETS"<<std::endl;
  evt.getByToken(l1BoostedToken_, l1Boosted);
  const vector<l1extra::L1JetParticle> &l1B = *l1Boosted;
  for(auto obj : l1B) {
    TLorentzVector temp;
    temp.SetPtEtaPhiE(obj.pt(), obj.eta(), obj.phi(), obj.pt());
    l1Jets->push_back(temp);
  }

  // Start Running Analysis
  //Handle<vector<reco::CaloJet> > jets;
  //if(evt.getByToken(jetSrc_, jets)){//Begin Getting Reco Jets
  //  for (const reco::CaloJet &jet : *jets) {
  //    if(jet.pt() > recoPt_ ) {
  //      goodJets.push_back(jet);
  //    }
  //  }
  //}
  //else
  //  cout<<"Error getting calo jets"<<std::endl;

  //Handle<vector<pat::Jet> > jetsAK8;

  //if(evt.getByToken(jetSrcAK8_, jetsAK8)){//Begin Getting AK8 Jets
  //  for (const pat::Jet &jetAK8 : *jetsAK8) {
  //    if(jetAK8.pt() > recoPt_ ) {
  //      nSubJets.push_back(jetAK8.subjets("SoftDropPuppi").size());
  //      nBHadrons.push_back(jetAK8.jetFlavourInfo().getbHadrons().size());
  //      TLorentzVector temp ;
  //      temp.SetPtEtaPhiE(jetAK8.pt(),jetAK8.eta(),jetAK8.phi(),jetAK8.et());
  //      ak8Jets->push_back(temp);
  //      if(jetAK8.subjets("SoftDropPuppi").size() ==  2 && jetAK8.jetFlavourInfo().getbHadrons().size() > 1){
  //        goodJetsAK8.push_back(jetAK8);
  //      }
  //    }
  //  }
  //}
  //else
  //  cout<<"Error getting AK8 jets"<<std::endl;

  zeroOutAllVariables();
  if(goodJetsAK8.size()>0){

    for(auto jet:goodJetsAK8){
      tau1.push_back(jet.userFloat("NjettinessAK8Puppi:tau1"));
      tau2.push_back(jet.userFloat("NjettinessAK8Puppi:tau2"));
      tau3.push_back(jet.userFloat("NjettinessAK8Puppi:tau3"));
      HFlav.clear();
      for(unsigned int isub=0; isub<((jet.subjets("SoftDropPuppi")).size()); isub++){
        HFlav.push_back(jet.subjets("SoftDropPuppi")[isub]->hadronFlavour());
        TLorentzVector temp;
        temp.SetPtEtaPhiE(jet.subjets("SoftDropPuppi")[isub]->pt(),jet.subjets("SoftDropPuppi")[isub]->eta(),jet.subjets("SoftDropPuppi")[isub]->phi(),jet.subjets("SoftDropPuppi")[isub]->et());
        subJets->push_back(temp);
      }
      subJetHFlav.push_back(HFlav);
      // take more variables from here: https://github.com/gouskos/HiggsToBBNtupleProducerTool/blob/opendata_80X/NtupleAK8/src/FatJetInfoFiller.cc#L215-L217
      // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools
    }

    // Match to boosted jets and see if we can match subjettiness functions...
    vector<l1extra::L1JetParticle> l1JetsSorted;
    for( vector<l1extra::L1JetParticle>::const_iterator l1Jet = l1Boosted->begin(); l1Jet != l1Boosted->end(); l1Jet++ ){
      l1JetsSorted.push_back(*l1Jet);
    }
    if(l1JetsSorted.size() > 1){  std::sort(l1JetsSorted.begin(),l1JetsSorted.end(),compareByPt);}

    pat::Jet recoJet_1;

    recoPt_1  = goodJetsAK8.at(0).pt();
    recoEta_1 = goodJetsAK8.at(0).eta();
    recoPhi_1 = goodJetsAK8.at(0).phi();
    recoJet_1 = goodJetsAK8.at(0);

    int i = 0;
    int foundL1Jet_1 = 0;
    l1extra::L1JetParticle l1Jet_1;
    if(l1JetsSorted.size() > 0){
      for(auto jet : l1JetsSorted){
        if(reco::deltaR(jet, recoJet_1)<0.4 && foundL1Jet_1 == 0 ){
          l1Jet_1 = jet;
          l1Pt_1  = jet.pt();
          l1Eta_1 = jet.eta();
          l1Phi_1 = jet.phi();
          l1NthJet_1 = i;
          foundL1Jet_1 = 1;
        }
        i++;
      }
    }

    int j = 0;
    int foundSeed_1 = 0;
    if(seeds.size() > 0){
      for(auto seed : seeds){
        if(reco::deltaR(seed, recoJet_1)<0.4 && foundSeed_1 == 0 ){
          seedPt_1  = seed.pt();
          seedEta_1 = seed.eta();
          seedPhi_1 = seed.phi();
          seedNthJet_1 = j;
          foundSeed_1 = 1;
        }
        j++;
      }
    }
  }

  edm::Handle<reco::GenParticleCollection> genParticles;
  if(evt.getByToken(genSrc_, genParticles)){//Begin Getting Gen Particles
    for (reco::GenParticleCollection::const_iterator genparticle = genParticles->begin(); genparticle != genParticles->end(); genparticle++){
      double DR = reco::deltaR(recoEta_1, recoPhi_1, genparticle->eta(), genparticle->phi());
      //if (DR < genDR && genparticle->status() > 21 && genparticle->status() < 41){
      if(genparticle->pdgId() == 25){
        genDR = DR;
        genId = genparticle->pdgId();
        genMother = genparticle->motherRef(0)->pdgId();
        genPt_1 = genparticle->pt();
        genEta_1 = genparticle->eta();
        genPhi_1 = genparticle->phi();
        genM_1 = genparticle->mass();
      }
    }
  }

  count = 0;
  double L1DR = 0.4;
  int L1Match = 99;
  for( vector<TLorentzVector>::const_iterator l1Jet = allL1Jets->begin(); l1Jet != allL1Jets->end(); l1Jet++ ){
    allL1Signals.push_back(false);
    allL1DR_withReco.push_back(reco::deltaR(recoEta_1, recoPhi_1, l1Jet->Eta(), l1Jet->Phi())); // store distance from the AK8 jet matched to the Higgs
    double DR = reco::deltaR(genEta_1, genPhi_1, l1Jet->Eta(), l1Jet->Phi()); // consider distance from the gen Higgs
    if(DR < L1DR) { L1Match = count; L1DR = DR; }
    count++;
  }
  if(L1Match < 99) allL1Signals[L1Match] = true; // only one boosted Higgs in each event!

  efficiencyTree->Fill();
  inputRegions.clear();
}

void BoostedJetStudies::zeroOutAllVariables(){
  genPt_1=-99; genEta_1=-99; genPhi_1=-99; genM_1=-99; genDR=99; genId=-99; genMother=-99;
  seedPt_1=-99; seedEta_1=-99; seedPhi_1=-99; seedNthJet_1=-99;
  recoPt_1=-99; recoEta_1=-99; recoPhi_1=-99; recoNthJet_1=-99;
  l1Pt_1=-99; l1Eta_1=-99; l1Phi_1=-99; l1NthJet_1=-99;
  recoPt_=-99;
}

void BoostedJetStudies::createBranches(TTree *tree){
    tree->Branch("run",     &run,     "run/I");
    tree->Branch("lumi",    &lumi,    "lumi/I");
    tree->Branch("event",   &event,   "event/I");
    tree->Branch("genPt_1",       &genPt_1,     "genPt_1/D");
    tree->Branch("genEta_1",      &genEta_1,    "genEta_1/D");
    tree->Branch("genPhi_1",      &genPhi_1,    "genPhi_1/D");
    tree->Branch("genM_1",        &genM_1,      "genM_1/D");
    tree->Branch("genDR",         &genDR,       "genDR/D");
    tree->Branch("genId",         &genId,       "genId/I");
    tree->Branch("genMother",     &genMother,   "genMother/I");
    tree->Branch("seedPt_1",      &seedPt_1,     "seedPt_1/D");
    tree->Branch("seedEta_1",     &seedEta_1,    "seedEta_1/D");
    tree->Branch("seedPhi_1",     &seedPhi_1,    "seedPhi_1/D");
    tree->Branch("seedNthJet_1",  &seedNthJet_1, "seedNthJet_1/I");
    tree->Branch("recoPt_1",      &recoPt_1,     "recoPt_1/D");
    tree->Branch("recoEta_1",     &recoEta_1,    "recoEta_1/D");
    tree->Branch("recoPhi_1",     &recoPhi_1,    "recoPhi_1/D");
    tree->Branch("recoNthJet_1",  &recoNthJet_1, "recoNthJet_1/I");
    tree->Branch("l1Pt_1",        &l1Pt_1,       "l1Pt_1/D"); 
    tree->Branch("l1Eta_1",       &l1Eta_1,      "l1Eta_1/D");
    tree->Branch("l1Phi_1",       &l1Phi_1,      "l1Phi_1/D");
    tree->Branch("l1NthJet_1",    &l1NthJet_1,   "l1NthJet_1/I");
    tree->Branch("tau1",          &tau1);
    tree->Branch("tau2",          &tau2);
    tree->Branch("tau3",          &tau3);
    tree->Branch("nSubJets",      &nSubJets);
    tree->Branch("subJetHFlav",   &subJetHFlav);
    tree->Branch("nBHadrons",     &nBHadrons);
    tree->Branch("bp_jets", "vector<TLorentzVector>", &l1Jets, 32000, 0);
    tree->Branch("seed180", "vector<TLorentzVector>", &seed180, 32000, 0);
    tree->Branch("tauseed", "vector<TLorentzVector>", &tauseed, 32000, 0);
    tree->Branch("ak8Jets", "vector<TLorentzVector>", &ak8Jets, 32000, 0);
    tree->Branch("subJets", "vector<TLorentzVector>", &subJets, 32000, 0);
    tree->Branch("allRegions", "vector<TLorentzVector>", &allRegions, 32000, 0);
    tree->Branch("allL1Jets", "vector<TLorentzVector>", &allL1Jets, 32000, 0);
    tree->Branch("signals",     &allL1Signals);
    tree->Branch("allL1DR_withReco", &allL1DR_withReco);
    tree->Branch("regionEta",        &regionEta);
    tree->Branch("regionPhi",        &regionPhi);
    tree->Branch("jetRegionEt",      &jetRegionEt);
    tree->Branch("jetRegionEGVeto",  &jetRegionEGVeto);
    tree->Branch("jetRegionTauVeto", &jetRegionTauVeto );

  }


// ------------ method called once each job just before starting event loop  ------------
void 
BoostedJetStudies::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
BoostedJetStudies::endJob() {
}

// ------------ method called when starting to processes a run  ------------

void
BoostedJetStudies::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BoostedJetStudies::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BoostedJetStudies);
