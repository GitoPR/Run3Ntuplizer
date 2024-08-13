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
#include "DataFormats/L1Trigger/interface/EtSum.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"


using namespace l1extra;
using namespace std;

bool compareByPt (l1extra::L1JetParticle i, l1extra::L1JetParticle j) { return(i.pt()>j.pt()); };

  constexpr std::array<double, 42> fillTwrEtaValues() {
    std::array<double, 42> twrEtaValues = {{0}};
    twrEtaValues[0] = 0;
    for (unsigned int i = 0; i < 20; i++) {
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
    twrEtaValues[30] = (3.15 + 2.98) / 2.;
    twrEtaValues[31] = (3.33 + 3.15) / 2.;
    twrEtaValues[32] = (3.50 + 3.33) / 2.;
    twrEtaValues[33] = (3.68 + 3.50) / 2.;
    twrEtaValues[34] = (3.68 + 3.85) / 2.;
    twrEtaValues[35] = (3.85 + 4.03) / 2.;
    twrEtaValues[36] = (4.03 + 4.20) / 2.;
    twrEtaValues[37] = (4.20 + 4.38) / 2.;
    twrEtaValues[38] = (4.74 + 4.38 * 3) / 4.;
    twrEtaValues[39] = (4.38 + 4.74 * 3) / 4.;
    twrEtaValues[40] = (5.21 + 4.74 * 3) / 4.;
    twrEtaValues[41] = (4.74 + 5.21 * 3) / 4.;
    return twrEtaValues;
  }
  constexpr std::array<double, 42> twrEtaValues = fillTwrEtaValues();


double getUCTTowerEta(int caloEta) {
  uint32_t absCaloEta = std::abs(caloEta);
  if (absCaloEta <= 41) {
    if (caloEta < 0)
      return -twrEtaValues[absCaloEta];
    else
      return +twrEtaValues[absCaloEta];
  } else
    return -999.;
}

double getUCTTowerPhi(int caloPhi) {
  if (caloPhi < 1)
    return -999.;
  else if (caloPhi > 72)
    return +999.;
  uint32_t absCaloPhi = std::abs(caloPhi) - 1;
  if (absCaloPhi < 36)
    return (((double)absCaloPhi + 0.5) * 0.0872);
  else
    return (-(71.5 - (double)absCaloPhi) * 0.0872);
}

namespace l1tcalo {
  constexpr uint32_t RegionETMask{0x000003FF};
  constexpr uint32_t RegionEGVeto{0x00000400};
  constexpr uint32_t RegionTauVeto{0x00000800};
  constexpr uint32_t HitTowerBits{0x0000F000};
  constexpr uint32_t RegionNoBits{0x000F0000};
  constexpr uint32_t CardNoBits{0x00700000};
  constexpr uint32_t CrateNoBits{0x01800000};
  constexpr uint32_t NegEtaBit{0x80000000};
  constexpr uint32_t LocationBits{0xFFFFF000};
  constexpr uint32_t LocationShift{12};
  constexpr uint32_t RegionNoShift{16};
  constexpr uint32_t CardNoShift{20};
  constexpr uint32_t CrateNoShift{23};
}  // namespace l1tcalo

typedef int loop; //loop type(i guess)                                                                                                                                                               
//Jet class                                                                                                                                                                                                 

class jetInfo{
public:
  int seedEnergy;
  int energy;
  int phiMax;
  int etaMax;

  jetInfo(){
    seedEnergy = 0;
    energy = 0;
    phiMax = 0;
    etaMax = 0;
  }

  jetInfo& operator=(const jetInfo& rhs){
    seedEnergy = rhs.seedEnergy;
    energy = rhs.energy;
    phiMax = rhs.phiMax;
    etaMax = rhs.etaMax;
    return *this;
  }
};

//
// class declaration
//

static constexpr int nSTPhi = 18;
static constexpr int nSTEta = 14;

namespace gctobj {

class towerMax {
  public:
  int energy;
  int iphi;
  int ieta;
  int towerEta;
  int towerPhi; 
  double eta;
  double phi; 
    towerMax() {
      energy = 0;
      phi = 0;
      eta = 0;
      towerEta = 0 ;
      towerPhi = 0;
      ieta = 0;
      iphi = 0;
    }
  };

  typedef struct {
    int et;
    int ieta;
    int iphi;
    int towerEta; 
    int towerPhi;
    double eta;
    double phi; 
  } GCTsupertower_t;


  typedef struct {
    GCTsupertower_t cr[nSTPhi];
  } etaStrip_t;

  typedef struct {
    GCTsupertower_t pk[nSTEta];
  } etaStripPeak_t;

  inline GCTsupertower_t bestOf2(const GCTsupertower_t& calotp0, const GCTsupertower_t& calotp1) {
   GCTsupertower_t x;
    x = (calotp0.et > calotp1.et) ? calotp0 : calotp1;
    return x;
  }


  inline GCTsupertower_t getPeakBin18N(const etaStrip_t& etaStrip) {
    GCTsupertower_t best01 = bestOf2(etaStrip.cr[0], etaStrip.cr[1]);
    GCTsupertower_t best23 = bestOf2(etaStrip.cr[2], etaStrip.cr[3]);
    GCTsupertower_t best45 = bestOf2(etaStrip.cr[4], etaStrip.cr[5]);
    GCTsupertower_t best67 = bestOf2(etaStrip.cr[6], etaStrip.cr[7]);
    GCTsupertower_t best89 = bestOf2(etaStrip.cr[8], etaStrip.cr[9]);
    GCTsupertower_t best1011 = bestOf2(etaStrip.cr[10], etaStrip.cr[11]);
    GCTsupertower_t best1213 = bestOf2(etaStrip.cr[12], etaStrip.cr[13]);
    GCTsupertower_t best1415 = bestOf2(etaStrip.cr[14], etaStrip.cr[15]);
    GCTsupertower_t best1617 = bestOf2(etaStrip.cr[16], etaStrip.cr[17]);



    GCTsupertower_t best0123 = bestOf2(best01, best23);
    GCTsupertower_t best4567 = bestOf2(best45, best67);
    GCTsupertower_t best891011 = bestOf2(best89, best1011);
    GCTsupertower_t best12131415 = bestOf2(best1213, best1415);

    GCTsupertower_t best0to7 = bestOf2(best0123, best4567);
    GCTsupertower_t best8to15 = bestOf2(best12131415, best891011);

    GCTsupertower_t best8to17 = bestOf2(best8to15, best1617);

    GCTsupertower_t bestOf18 = bestOf2(best0to7, best8to17);

    return bestOf18;
 }

   inline towerMax getPeakBin14N(const etaStripPeak_t& etaStrip) {
    towerMax x;

    GCTsupertower_t best01 = bestOf2(etaStrip.pk[0], etaStrip.pk[1]);
    GCTsupertower_t best23 = bestOf2(etaStrip.pk[2], etaStrip.pk[3]);
    GCTsupertower_t best45 = bestOf2(etaStrip.pk[4], etaStrip.pk[5]);
    GCTsupertower_t best67 = bestOf2(etaStrip.pk[6], etaStrip.pk[7]);
    GCTsupertower_t best89 = bestOf2(etaStrip.pk[8], etaStrip.pk[9]);
    GCTsupertower_t best1011 = bestOf2(etaStrip.pk[10], etaStrip.pk[11]);
    GCTsupertower_t best1213 = bestOf2(etaStrip.pk[12], etaStrip.pk[13]);

    GCTsupertower_t best0123 = bestOf2(best01, best23);
    GCTsupertower_t best4567 = bestOf2(best45, best67);
    GCTsupertower_t best891011 = bestOf2(best89, best1011);

    GCTsupertower_t best8to13 = bestOf2(best891011, best1213);
    GCTsupertower_t best0to7 = bestOf2(best0123, best4567);

    GCTsupertower_t bestOf14 = bestOf2(best0to7, best8to13);

    x.energy = bestOf14.et;
    x.iphi = bestOf14.iphi;
    x.ieta = bestOf14.ieta;
    x.eta = bestOf14.eta;
    x.phi = bestOf14.phi;
    x.towerEta = bestOf14.towerEta;
    x.towerPhi = bestOf14.towerPhi;
    return x;
  }

    inline towerMax getTowerMax(GCTsupertower_t temp[nSTEta][nSTPhi]) {
    etaStripPeak_t etaStripPeak;

    for (int i = 0; i < nSTEta; i++) {
      etaStrip_t test;
      for (int j = 0; j < nSTPhi; j++) {
        test.cr[j] = temp[i][j];
      }
      etaStripPeak.pk[i] = getPeakBin18N(test);
    }

    towerMax peakIn14;
    peakIn14 = getPeakBin14N(etaStripPeak);
    return peakIn14;
  }


} // namespace gctobj

jetInfo getJetValues(gctobj::GCTsupertower_t tempX[nSTEta][nSTPhi], int seed_eta, int seed_phi ){

  int temp[nSTEta+2][nSTPhi+2] ;

  int eta_slice[3] ;

  jetInfo jet_tmp;



  for(loop i=0; i<nSTEta+2; i++){
    for(loop k=0; k<nSTPhi+2; k++){
      temp[i][k] = 0 ;
    }
  }

  for(loop i=0; i<nSTEta; i++){

    //    std::cout<< "tempX[i][17].et : " << tempX[i][17].et << "\n"<< "tempX[i][0].et : " << tempX[i][0].et << "\n"<<std::endl;
    
    temp[i+1][0] = tempX[i][17].et;
    temp[i+1][19]= tempX[i][0].et;

    
    for(loop k=0; k<nSTPhi; k++){
      temp[i+1][k+1] = tempX[i][k].et ;
      //      std::cout << "tempX[i][k].et : " << tempX[i][k].et << std::endl; 
    }
  }


  int seed_eta1,  seed_phi1 ;

  seed_eta1 = seed_eta ; //to start from corner
  seed_phi1 = seed_phi ; //to start from corner
  int tmp1, tmp2, tmp3 ;

  for(loop j=0; j<nSTEta; j++){
    for(loop k=0; k<nSTPhi; k++){
      if(j== seed_eta1 && k == seed_phi1){
	std::cout << "seed_eta1 : " << j  << "\t" << "seed_phi1 : " << k << std::endl;
        for(loop m=0; m<3 ; m++){
          tmp1 = temp[j+m][k] ;
          tmp2 = temp[j+m][k+1] ;
          tmp3 = temp[j+m][k+2] ;
	  //std::cout << "tmp1 : " << tmp1 << "\t" << "tmp2 : " << tmp2 << "\t" << "tmp3 : " << tmp3 << "\n" << std::endl;
          eta_slice[m] = tmp1 + tmp2 + tmp3 ; // Sum the energies of 3 3x1 adjacent slices to make the 3x3.
	} 
      }
    }
  }
  
 jet_tmp.energy=eta_slice[0] + eta_slice[1] + eta_slice[2];
 //   std::cout << "eta_slice[0] : " << eta_slice[0] << std::endl;
 //std::cout << "eta_slice[1] : " << eta_slice[1] << std::endl;
 //std::cout << "eta_slice[2] : " << eta_slice[2] << std::endl;
 //    std::cout << " jet_tmp.energy : " <<  jet_tmp.energy <<"\n" <<std::endl; 
  

  for(loop i=0; i<nSTEta; i++){
    if(i+1>=seed_eta && i<=seed_eta+1){
      for(loop k=0; k<nSTPhi; k++){
        if(k+1>=seed_phi && k<=seed_phi+1)  tempX[i][k].et = 0 ; // set the 3x3 energies to 0
      } 
    }
  }



  return jet_tmp ;
} //end of the getJetValues function

typedef struct
{
	int iphi; // -41 to 41
	int ieta; // -28 to 28
	int side; // 1: negative side, 0: positive side

} calo_coor_t;

//for converting ieta to tower eta:
const calo_coor_t calo_coor[252] =
{
		{ 71 , 25 , 1 },
		{ 71 , 21 , 1 },
		{ 71 , 17 , 1 },
		{ 71 , 13 , 1 },
		{ 71 , 9 , 1 },
		{ 71 , 5 , 1 },
		{ 71 , 1 , 1 },
		{ 71 , 1 , 0 },
		{ 71 , 5 , 0 },
		{ 71 , 9 , 0 },
		{ 71 , 13 , 0 },
		{ 71 , 17 , 0 },
		{ 71 , 21 , 0 },
		{ 71 , 25 , 0 },
		{ 3 , 25 , 1 },
		{ 3 , 21 , 1 },
		{ 3 , 17 , 1 },
		{ 3 , 13 , 1 },
		{ 3 , 9 , 1 },
		{ 3 , 5 , 1 },
		{ 3 , 1 , 1 },
		{ 3 , 1 , 0 },
		{ 3 , 5 , 0 },
		{ 3 , 9 , 0 },
		{ 3 , 13 , 0 },
		{ 3 , 17 , 0 },
		{ 3 , 21 , 0 },
		{ 3 , 25 , 0 },
		{ 7 , 25 , 1 },
		{ 7 , 21 , 1 },
		{ 7 , 17 , 1 },
		{ 7 , 13 , 1 },
		{ 7 , 9 , 1 },
		{ 7 , 5 , 1 },
		{ 7 , 1 , 1 },
		{ 7 , 1 , 0 },
		{ 7 , 5 , 0 },
		{ 7 , 9 , 0 },
		{ 7 , 13 , 0 },
		{ 7 , 17 , 0 },
		{ 7 , 21 , 0 },
		{ 7 , 25 , 0 },
		{ 11 , 25 , 1 },
		{ 11 , 21 , 1 },
		{ 11 , 17 , 1 },
		{ 11 , 13 , 1 },
		{ 11 , 9 , 1 },
		{ 11 , 5 , 1 },
		{ 11 , 1 , 1 },
		{ 11 , 1 , 0 },
		{ 11 , 5 , 0 },
		{ 11 , 9 , 0 },
		{ 11 , 13 , 0 },
		{ 11 , 17 , 0 },
		{ 11 , 21 , 0 },
		{ 11 , 25 , 0 },
		{ 15 , 25 , 1 },
		{ 15 , 21 , 1 },
		{ 15 , 17 , 1 },
		{ 15 , 13 , 1 },
		{ 15 , 9 , 1 },
		{ 15 , 5 , 1 },
		{ 15 , 1 , 1 },
		{ 15 , 1 , 0 },
		{ 15 , 5 , 0 },
		{ 15 , 9 , 0 },
		{ 15 , 13 , 0 },
		{ 15 , 17 , 0 },
		{ 15 , 21 , 0 },
		{ 15 , 25 , 0 },
		{ 19 , 25 , 1 },
		{ 19 , 21 , 1 },
		{ 19 , 17 , 1 },
		{ 19 , 13 , 1 },
		{ 19 , 9 , 1 },
		{ 19 , 5 , 1 },
		{ 19 , 1 , 1 },
		{ 19 , 1 , 0 },
		{ 19 , 5 , 0 },
		{ 19 , 9 , 0 },
		{ 19 , 13 , 0 },
		{ 19 , 17 , 0 },
		{ 19 , 21 , 0 },
		{ 19 , 25 , 0 },
		{ 23 , 25 , 1 },
		{ 23 , 21 , 1 },
		{ 23 , 17 , 1 },
		{ 23 , 13 , 1 },
		{ 23 , 9 , 1 },
		{ 23 , 5 , 1 },
		{ 23 , 1 , 1 },
		{ 23 , 1 , 0 },
		{ 23 , 5 , 0 },
		{ 23 , 9 , 0 },
		{ 23 , 13 , 0 },
		{ 23 , 17 , 0 },
		{ 23 , 21 , 0 },
		{ 23 , 25 , 0 },
		{ 27 , 25 , 1 },
		{ 27 , 21 , 1 },
		{ 27 , 17 , 1 },
		{ 27 , 13 , 1 },
		{ 27 , 9 , 1 },
		{ 27 , 5 , 1 },
		{ 27 , 1 , 1 },
		{ 27 , 1 , 0 },
		{ 27 , 5 , 0 },
		{ 27 , 9 , 0 },
		{ 27 , 13 , 0 },
		{ 27 , 17 , 0 },
		{ 27 , 21 , 0 },
		{ 27 , 25 , 0 },
		{ 31 , 25 , 1 },
		{ 31 , 21 , 1 },
		{ 31 , 17 , 1 },
		{ 31 , 13 , 1 },
		{ 31 , 9 , 1 },
		{ 31 , 5 , 1 },
		{ 31 , 1 , 1 },
		{ 31 , 1 , 0 },
		{ 31 , 5 , 0 },
		{ 31 , 9 , 0 },
		{ 31 , 13 , 0 },
		{ 31 , 17 , 0 },
		{ 31 , 21 , 0 },
		{ 31 , 25 , 0 },
		{ 35 , 25 , 1 },
		{ 35 , 21 , 1 },
		{ 35 , 17 , 1 },
		{ 35 , 13 , 1 },
		{ 35 , 9 , 1 },
		{ 35 , 5 , 1 },
		{ 35 , 1 , 1 },
		{ 35 , 1 , 0 },
		{ 35 , 5 , 0 },
		{ 35 , 9 , 0 },
		{ 35 , 13 , 0 },
		{ 35 , 17 , 0 },
		{ 35 , 21 , 0 },
		{ 35 , 25 , 0 },
		{ 39 , 25 , 1 },
		{ 39 , 21 , 1 },
		{ 39 , 17 , 1 },
		{ 39 , 13 , 1 },
		{ 39 , 9 , 1 },
		{ 39 , 5 , 1 },
		{ 39 , 1 , 1 },
		{ 39 , 1 , 0 },
		{ 39 , 5 , 0 },
		{ 39 , 9 , 0 },
		{ 39 , 13 , 0 },
		{ 39 , 17 , 0 },
		{ 39 , 21 , 0 },
		{ 39 , 25 , 0 },
		{ 43 , 25 , 1 },
		{ 43 , 21 , 1 },
		{ 43 , 17 , 1 },
		{ 43 , 13 , 1 },
		{ 43 , 9 , 1 },
		{ 43 , 5 , 1 },
		{ 43 , 1 , 1 },
		{ 43 , 1 , 0 },
		{ 43 , 5 , 0 },
		{ 43 , 9 , 0 },
		{ 43 , 13 , 0 },
		{ 43 , 17 , 0 },
		{ 43 , 21 , 0 },
		{ 43 , 25 , 0 },
		{ 47 , 25 , 1 },
		{ 47 , 21 , 1 },
		{ 47 , 17 , 1 },
		{ 47 , 13 , 1 },
		{ 47 , 9 , 1 },
		{ 47 , 5 , 1 },
		{ 47 , 1 , 1 },
		{ 47 , 1 , 0 },
		{ 47 , 5 , 0 },
		{ 47 , 9 , 0 },
		{ 47 , 13 , 0 },
		{ 47 , 17 , 0 },
		{ 47 , 21 , 0 },
		{ 47 , 25 , 0 },
		{ 51 , 25 , 1 },
		{ 51 , 21 , 1 },
		{ 51 , 17 , 1 },
		{ 51 , 13 , 1 },
		{ 51 , 9 , 1 },
		{ 51 , 5 , 1 },
		{ 51 , 1 , 1 },
		{ 51 , 1 , 0 },
		{ 51 , 5 , 0 },
		{ 51 , 9 , 0 },
		{ 51 , 13 , 0 },
		{ 51 , 17 , 0 },
		{ 51 , 21 , 0 },
		{ 51 , 25 , 0 },
		{ 55 , 25 , 1 },
		{ 55 , 21 , 1 },
		{ 55 , 17 , 1 },
		{ 55 , 13 , 1 },
		{ 55 , 9 , 1 },
		{ 55 , 5 , 1 },
		{ 55 , 1 , 1 },
		{ 55 , 1 , 0 },
		{ 55 , 5 , 0 },
		{ 55 , 9 , 0 },
		{ 55 , 13 , 0 },
		{ 55 , 17 , 0 },
		{ 55 , 21 , 0 },
		{ 55 , 25 , 0 },
		{ 59 , 25 , 1 },
		{ 59 , 21 , 1 },
		{ 59 , 17 , 1 },
		{ 59 , 13 , 1 },
		{ 59 , 9 , 1 },
		{ 59 , 5 , 1 },
		{ 59 , 1 , 1 },
		{ 59 , 1 , 0 },
		{ 59 , 5 , 0 },
		{ 59 , 9 , 0 },
		{ 59 , 13 , 0 },
		{ 59 , 17 , 0 },
		{ 59 , 21 , 0 },
		{ 59 , 25 , 0 },
		{ 63 , 25 , 1 },
		{ 63 , 21 , 1 },
		{ 63 , 17 , 1 },
		{ 63 , 13 , 1 },
		{ 63 , 9 , 1 },
		{ 63 , 5 , 1 },
		{ 63 , 1 , 1 },
		{ 63 , 1 , 0 },
		{ 63 , 5 , 0 },
		{ 63 , 9 , 0 },
		{ 63 , 13 , 0 },
		{ 63 , 17 , 0 },
		{ 63 , 21 , 0 },
		{ 63 , 25 , 0 },
		{ 67 , 25 , 1 },
		{ 67 , 21 , 1 },
		{ 67 , 17 , 1 },
		{ 67 , 13 , 1 },
		{ 67 , 9 , 1 },
		{ 67 , 5 , 1 },
		{ 67 , 1 , 1 },
		{ 67 , 1 , 0 },
		{ 67 , 5 , 0 },
		{ 67 , 9 , 0 },
		{ 67 , 13 , 0 },
		{ 67 , 17 , 0 },
		{ 67 , 21 , 0 },
		{ 67 , 25 , 0 }

};



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

  //  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<vector<reco::CaloJet> > jetSrc_;
  edm::EDGetTokenT<vector<pat::Jet> > jetSrcAK8_;
  edm::EDGetTokenT<reco::GenParticleCollection> genSrc_;

  edm::EDGetTokenT<BXVector<l1t::Jet>> stage2JetToken_;
  edm::EDGetTokenT<BXVector<l1t::Tau>> stage2TauToken_;
  edm::EDGetTokenT<l1t::EtSumBxCollection> stage2EtSumToken_;
  edm::EDGetTokenT<vector<l1extra::L1JetParticle>> l1BoostedToken_;
  //add calo region token 
  edm::EDGetTokenT<std::vector <L1CaloRegion> > regionsToken_;

  TH1F* nEvents;

  int run, lumi, event;

  double genPt_1, genEta_1, genPhi_1, genM_1, genDR, dau_eta, dauET, dau_phi;
  int genId, genMother, dauID;
  double recoPt_1, recoEta_1, recoPhi_1;
  double l1Pt_1, l1Eta_1, l1Phi_1;
  double jetClusterPt, jetClusterEta, jetClusterPhi ;
  double seedPt_1, seedEta_1, seedPhi_1;
  std::vector<double> clusterCord; 
  
  int l1NthJet_1;
  int recoNthJet_1;
  int seedNthJet_1;

  // add array of cregions
  std::vector<uint16_t> cregions;

  // Daugther particles info
   std::vector<int> dgt_id;
  std::vector<double>dgt_et;
  std::vector<double>dgt_eta;
  std::vector<double>dgt_phi;

    
  double recoPt_;
  std::vector<int> nSubJets, nBHadrons, HFlav;
  std::vector<std::vector<int>> subJetHFlav;
  std::vector<float> tau1, tau2, tau3;

  std::vector<TLorentzVector> *l1Jets  = new std::vector<TLorentzVector>;
  std::vector<TLorentzVector> *seed180  = new std::vector<TLorentzVector>;
  std::vector<TLorentzVector> *tauseed  = new std::vector<TLorentzVector>;
  std::vector<TLorentzVector> *ak8Jets  = new std::vector<TLorentzVector>;
  std::vector<TLorentzVector> *subJets  = new std::vector<TLorentzVector>;

  void createBranches(TTree *tree);
  TTree* efficiencyTree;
  TH1F* leadL1Pt;
  TH1F* leadL1Eta;
  TH1F* leadL1Phi;
  TH1F* leadSingleJet;
  TH1F* EtSum_HT;
  TH1F* EtSum_ETMHF;
  edm::Service<TFileService> tfs_;  

};

BoostedJetStudies::BoostedJetStudies(const edm::ParameterSet& iConfig) :
  jetSrc_(    consumes<vector<reco::CaloJet> >(iConfig.getParameter<edm::InputTag>("recoJets"))),
  jetSrcAK8_( consumes<vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("recoJetsAK8"))),
  genSrc_( consumes<reco::GenParticleCollection> (iConfig.getParameter<edm::InputTag>( "genParticles"))),
  stage2JetToken_(consumes<BXVector<l1t::Jet>>( edm::InputTag("caloStage2Digis","Jet","RECO"))),
  stage2TauToken_(consumes<BXVector<l1t::Tau>>( edm::InputTag("caloStage2Digis","Tau","RECO"))),
  stage2EtSumToken_(consumes<l1t::EtSumBxCollection>( edm::InputTag("caloStage2Digis","EtSum","RECO"))),
  l1BoostedToken_(consumes<vector<l1extra::L1JetParticle>>( edm::InputTag("simCaloStage2Layer1Summary","Boosted",""))),
  regionsToken_(consumes<std::vector <L1CaloRegion> >(iConfig.getUntrackedParameter<edm::InputTag>("UCTRegion")))

{
  // Initialize the Tree

  recoPt_      = iConfig.getParameter<double>("recoPtCut");
  nEvents      = tfs_->make<TH1F>( "nEvents"  , "nEvents", 2,  0., 1. );
  efficiencyTree = tfs_->make<TTree>("efficiencyTree", "Gen Matched Jet Tree ");
  createBranches(efficiencyTree);
  leadL1Pt     = tfs_->make<TH1F>( "leadL1Pt", "leadL1Pt", 100, 0., 1100.);
  leadL1Eta    = tfs_->make<TH1F>( "leadL1Eta", "leadL1Eta", 100, -5., 5.);
  leadL1Phi    = tfs_->make<TH1F>( "leadL1Phi", "leadL1Phi", 100, -M_PI, M_PI);
  leadSingleJet= tfs_->make<TH1F>( "leadSingleJet", "leadSingleJet", 100, 0., 1100.);
  EtSum_HT     = tfs_->make<TH1F>( "EtSum_HT", "EtSum_HT", 100, 0., 1100.);
  EtSum_ETMHF  = tfs_->make<TH1F>( "EtSum_ETMHF", "EtSum_ETMHF", 100, 0., 1100.);
}

BoostedJetStudies::~BoostedJetStudies() {
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

  
  uint16_t regionColl[252];
  uint16_t regionColl_input[14][18];
  
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
  cregions.clear();
  dgt_id.clear();
  dgt_eta.clear();
  dgt_phi.clear();
  dgt_et.clear();
  clusterCord.clear(); 
  //  maxEta.clear();
  //maxPhi.clear(); 

  gctobj::GCTsupertower_t temp[nSTEta][nSTPhi]; 
  for (const auto& region : evt.get(regionsToken_)){

    uint32_t ieta = region.id().ieta() - 4; // Subtract off the offset for HF
    uint32_t iphi = region.id().iphi();
    double  et = region.et();
    //  std::cout << "et : " << et << std::endl;
    
    uint16_t regionSummary = region.raw();
    
    uint16_t rloc_eta  = ((0xFFFF & regionSummary) >> 14);
    uint16_t rloc_phi = ((0x3FFF & regionSummary) >> 12); 
    
    int calo_index = 14*iphi + ieta; // if we need to we can make another index table to reduce the cost of resources of this operation.

    calo_coor_t calo_coor_event = calo_coor[calo_index];
    //    std::cout << "calo_coor_event.side : " << calo_coor_event.side << std::endl; 

    int towerEta = calo_coor_event.ieta + rloc_eta;
    int towerPhi = calo_coor_event.iphi + rloc_phi;
      if (calo_coor_event.side > 0){ 
	towerEta = -towerEta; 
	
      }

      //    std::cout << "towerEta : " << towerEta << "\t" << "towerPhi : " << towerPhi << "\n" << std::endl;
    
    double eta =  getUCTTowerEta(towerEta);
    double phi = getUCTTowerPhi(towerPhi); 
      
    regionColl_input[ieta][iphi] = et;  //regionSummary;

    temp[ieta][iphi].ieta  = ieta;
    temp[ieta][iphi].iphi = iphi;
    temp[ieta][iphi].et = et;
    temp[ieta][iphi].towerEta = towerEta;
    temp[ieta][iphi].towerPhi = towerPhi;
    temp[ieta][iphi].eta = eta;
    temp[ieta][iphi].phi = phi; 
  }

  for (unsigned int phi = 0; phi < 18; phi++){
    for (int eta = 0; eta < 14; eta++){
      //      std::cout << "regionColl_input[eta][phi] : "<<regionColl_input[eta][phi] << std::endl; 
      cregions.push_back(regionColl_input[eta][phi]);
      
    }
  }

  
  // ++++++++++ For testing out regionSummary outputs +++++++++++++ 
  /*  for(int idx = 0 ; idx < 252 ; idx++){
    std::cout << idx << std::endl; 
    //std::cout << "rloc_eta: "<<((0xFFFF & cregions[idx]) >> 14)<<"\t"<<"rloc_phi: "<<((0x3FFF & cregions[idx]) >> 12)<< std::endl;
    uint32_t location = ((cregions[idx] & l1tcalo::LocationBits) >> l1tcalo::LocationShift);
    bool eleBit = !((l1tcalo::RegionEGVeto & cregions[idx]) == l1tcalo::RegionEGVeto);
    bool tauBit = !((l1tcalo::RegionTauVeto & cregions[idx]) == l1tcalo::RegionTauVeto);
    uint32_t hitTowerLocation = (location & 0xF); 
    //    std::cout << "Hello" << std::endl;
    //    value  = cregions[idx] ;

    uint16_t rloc_eta  = ((0xFFFF & cregions[idx]) >> 14);
    uint16_t rloc_phi = ((0x3FFF & cregions[idx]) >> 12);
    
    towerEta = 4*
    //std::cout<< "\t"<<"et: "<<cregions[idx]<<"\t"<<"ieta : "<<idx%14<<"\t"<<"iphi: "<< idx/14 <<"\t"<<"rloc_eta: "<<((0xFFFF & cregions[idx]) >> 14)<<"\t"<<"rloc_phi: "<<((0x3FFF & cregions[idx]) >> 12) << "\t"<<"location: "<<hitTowerLocation<<"\t"<<"eleBit: "<<eleBit<<"\t"<<"tauBit: "<<tauBit<<std::endl;
    //std::cout<<"count: " <<count<< "\t"<< std::endl; // "et: "<<cregions[idx]<<"\t" << std::endl; 

    }*/ 

  

  

  // +++++++++++ for testing the results of maximum et finder  +++++++++++++++++++ 
  /*int size = 252; 

    double maxValue = cregions[0];
    int index = 0;

    for (int i = 1; i < size; ++i) {
        if (cregions[i] > maxValue) {
            maxValue = cregions[i];
            index = i;
        }
    }
    
    int test_ieta = index % 14; 
    int test_iphi =static_cast<int>( index / 14); 
    std::cout << "test max Et: " << maxValue << std::endl;
    std::cout << "test et Index: " << index << std::endl;
    std::cout << "test ieta: " << test_ieta << std::endl;
    std::cout << "test  Iphi: " << test_iphi << "\n" << std::endl; */


    // testing results
    gctobj::towerMax maxTower  = gctobj::getTowerMax(temp);
    /* std::cout << "algorithm max et = " << maxTower.energy <<std::endl;
    std::cout << "algorithm ieta = " << maxTower.eta << std::endl;
    std::cout <<"algorithm iphi = " << maxTower.phi << "\n" << std::endl; 
    std::cout << "checking the value of gctobj at the indeces given by test code: " << temp[test_ieta][test_iphi].et << "\n" << std::endl;  */

    jetInfo test_jet;

    test_jet.seedEnergy = maxTower.energy;
    test_jet.etaMax = maxTower.eta;
    test_jet.phiMax = maxTower.phi;

    std::cout << "maxTower.towerEta" << maxTower.towerEta<<"\t"  << "maxTower.towerPhi : " << maxTower.towerPhi << "\n" <<std::endl; 
    /*std::cout << "maxTower.ieta : " << maxTower.ieta << std::endl;
    std::cout << "maxTower.iphi : " << maxTower.iphi << std::endl;
    std::cout << "maxTower.eta : " << maxTower.eta << std::endl;
    std::cout << "maxTower.phi : " << maxTower.phi << "\n" << std::endl;*/
    clusterCord.insert(clusterCord.end(), {maxTower.eta ,maxTower.phi});
    jetClusterPt = test_jet.seedEnergy;
    jetClusterEta = test_jet.etaMax;
    jetClusterPhi = test_jet.phiMax;



     
     
    //    maxTower.seedEta*

    jetInfo tmp_jet;
    
    
    tmp_jet = getJetValues(temp,maxTower.ieta, maxTower.iphi);
    test_jet.energy = tmp_jet.energy;

    /*std::cout << "test_jet.energy : " << test_jet.energy << std::endl;
    std::cout << "test_jet.phi : " << test_jet.phiMax << std::endl;
    std::cout << "test_jet.eta : "<< test_jet.etaMax <<  "\n"<<std:: endl;*/
 
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
  if(seed180->size() > 0) leadSingleJet->Fill(seed180->at(0).Pt());

  edm::Handle<BXVector<l1t::Tau>> stage2Taus;
  if(!evt.getByToken(stage2TauToken_, stage2Taus)) cout<<"ERROR GETTING THE STAGE 2 TAUS"<<std::endl;
  evt.getByToken(stage2TauToken_, stage2Taus);
  const BXVector<l1t::Tau> &s2t = *stage2Taus;
  for(auto obj : s2t) {
    TLorentzVector temp;
    temp.SetPtEtaPhiE(obj.pt(), obj.eta(), obj.phi(), obj.pt());
    tauseed->push_back(temp);
  }

  edm::Handle<l1t::EtSumBxCollection> stage2EtSum;
  if(!evt.getByToken(stage2EtSumToken_, stage2EtSum)) cout<<"ERROR GETTING THE STAGE 2 ETSUM"<<endl;
  evt.getByToken(stage2EtSumToken_, stage2EtSum);
  for (l1t::EtSumBxCollection::const_iterator obj = stage2EtSum->begin(0); obj != stage2EtSum->end(0); obj++) {
    if(obj->getType() == l1t::EtSum::kTotalHt){
      EtSum_HT->Fill(obj->pt());
    }
    if(obj->getType() == l1t::EtSum::kMissingHtHF){
      EtSum_ETMHF->Fill(obj->pt());
    }
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
  if(l1Jets->size() > 0) {
    leadL1Pt->Fill(l1Jets->at(0).Pt());
    leadL1Eta->Fill(l1Jets->at(0).Eta());
    leadL1Phi->Fill(l1Jets->at(0).Phi());
  }

  // Start Runing Analysis
  Handle<vector<reco::CaloJet> > jets;
  if(evt.getByToken(jetSrc_, jets)){//Begin Getting Reco Jets
    for (const reco::CaloJet &jet : *jets) {
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
      //take more variables from here: https://github.com/gouskos/HiggsToBBNtupleProducerTool/blob/opendata_80X/NtupleAK8/src/FatJetInfoFiller.cc#L215-L217
      // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools
    }
    
    
    //Match to boosted jets and see if we can match subjettiness functions...
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
      if ( genparticle->pdgId() == 25 && genparticle->status() > 21 && genparticle->status() < 41){

	//	cout<< "check1"  << std::endl;
	genDR = DR;
	genId = genparticle->pdgId();
	genPt_1 = genparticle->pt();
	genEta_1 = genparticle->eta();
	genPhi_1 = genparticle->phi();
	genM_1 = genparticle->mass();
	/*   cout<< "genEta: " << genEta_1 << std::endl;
	      cout<<"genId :" << genId <<std::endl;
	      cout<< "genEta: " << genEta_1 << std::endl;
	      cout << "\n" << std::endl;  */ 

      
      }

      const reco::Candidate *m =      genparticle -> mother();
      //      cout<< "check2"  << std::endl;
      if (m != nullptr && m->pdgId() == 25){
	if(genparticle -> pdgId() != 25 && m -> pdgId()  == 25) {

	  //	  cout<< "check3"  << std::endl;
	  dauID  = genparticle -> pdgId();
	  dauET = genparticle -> pt();
	  dau_phi = genparticle -> phi();
	  dau_eta = genparticle -> eta();
	  
	  dgt_id.push_back(dauID);  // daughter pdgid                                                                                                                                                                               
	  dgt_eta.push_back(dau_eta);  // eta
	  dgt_phi.push_back(dau_phi); // phi                                                                                                                                                                                 
	  dgt_et.push_back(dauET); //phi                                                                                                                                                                                                         


      
	}
      }
      //      cout<< "check4"  << std::endl;

    }
  }
  if (abs(genEta_1) < 2.5) {
    efficiencyTree->Fill();
  }

  //  efficiencyTree->Fill();
  //  cout<< "check5"  << std::endl;
}


  void BoostedJetStudies::zeroOutAllVariables(){
    genPt_1=-99; genEta_1=-99; genPhi_1=-99; genM_1=-99; genDR=99; genId=-99; genMother=-99;
    seedPt_1=-99; seedEta_1=-99; seedPhi_1=-99; seedNthJet_1=-99;
    recoPt_1=-99; recoEta_1=-99; recoPhi_1=-99; recoNthJet_1=-99;
    l1Pt_1=-99; l1Eta_1=-99; l1Phi_1=-99; l1NthJet_1=-99;
    jetClusterPt = -99 ; jetClusterEta = -99;  jetClusterPhi = - 99; 
    recoPt_=-99; dauID = -99;  dauET = -99; dau_phi =-99; dau_eta = -99; // maxEta = -99 ;maxPhi = -99; 
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
    tree->Branch("jetClusterPt",  &jetClusterPt,       "jetClusterPt/D");
    tree->Branch("jetClusterEta",  &jetClusterEta,       "jetClusterEta/D");
    tree->Branch("jetClusterPhi",  &jetClusterPhi,       "jetClusterPhi/D");
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
    tree->Branch("l1Jets", "vector<TLorentzVector>", &l1Jets, 32000, 0);
    tree->Branch("seed180", "vector<TLorentzVector>", &seed180, 32000, 0);
    tree->Branch("tauseed", "vector<TLorentzVector>", &tauseed, 32000, 0);
    tree->Branch("ak8Jets", "vector<TLorentzVector>", &ak8Jets, 32000, 0);
    tree->Branch("subJets", "vector<TLorentzVector>", &subJets, 32000, 0);
    tree->Branch("cregions",     &cregions);
    tree->Branch("dgt_id",     &dgt_id);
    tree->Branch("dgt_et" , &dgt_et);
    tree->Branch("dgt_eta" , &dgt_eta);
    tree->Branch("dgt_phi" , &dgt_phi);
    tree->Branch("clusterCord" , &clusterCord);
    
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

  //void
  //BoostedJetStudies::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
  //}
 
  // ------------ method called when ending the processing of a run  ------------
  /*
    void
    BoostedJetStudies::endRun(edm::Run const&, edm::EventSetup const&)
    {
    }
  */
 
  // ------------ method called when starting to processes a luminosity block  ------------
  /*
    void
    BoostedJetStudies::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
    {
    }
  */
 
  // ------------ method called when ending the processing of a luminosity block  ------------
  /*
    void
    BoostedJetStudies::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
    {
    }
  */
 
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
 
