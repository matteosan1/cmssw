// -*- C++ -*-
//
// Package:    PileupWeightProducer
// Class:      PileupWeightProducer
// 
/**\class PileupWeightProducer PileupWeightProducer.cc 

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ricardo Vasquez Sierra,6 R-025,+41227672274,
//         Created:  Mon Nov 21 15:05:26 CET 2011
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/Utilities/interface/Lumi3DReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
#include <vector>
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include <iostream>
//
// class declaration
//

class PileupWeightProducer : public edm::EDProducer {
   public:
      explicit PileupWeightProducer(const edm::ParameterSet&);
      ~PileupWeightProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

  //bool firsttime_;
  bool hardcodedWeights_;
  std::string pileupMC_;
  std::string pileupData_;
  edm::LumiReWeighting LumiWeights_;
  edm::Lumi3DReWeighting LumiWeightsNominal_;
  edm::Lumi3DReWeighting LumiWeightsUp_;
  edm::Lumi3DReWeighting LumiWeightsDown_;
  std::vector< float > Data_;
  std::vector<float> MC_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
PileupWeightProducer::PileupWeightProducer(const edm::ParameterSet& iConfig) {

  //firsttime_ =  iConfig.getUntrackedParameter<bool>("FirstTime");
  hardcodedWeights_ = iConfig.getUntrackedParameter<bool>("hardcodedWeights");
  pileupMC_ = iConfig.existsAs<std::string>("PileupMCFile") ? iConfig.getParameter<std::string>("PileupMCFile") : "PUMC_dist.root" ;
  pileupData_ = iConfig.existsAs<std::string>("PileupDataFile") ? iConfig.getParameter<std::string>("PileupDataFile") : "PUData_dist.root" ;

  //register your products
  produces<std::vector<float> >( "pileupWeights" ).setBranchAlias( "pileupWeights" );

  if (hardcodedWeights_) {
    double Asymp50ns[60]= {0.000000,  0.000000,  0.000034,  0.000219,  0.001040,  0.002635,  0.006730,  0.012016,  0.021400,  0.031514,  0.043007,  0.055493,  0.064754,  0.074051,  0.077450,  0.079318,  0.076225,  0.072531,  0.065585,  0.058338,  0.050649,  0.042521,  0.035597,  0.029352,  0.023825,  0.017925,  0.014861,  0.011396,  0.008582,  0.006260,  0.004950,  0.003400,  0.002607,  0.001880,  0.001231,  0.000849,  0.000620,  0.000426,  0.000254,  0.000169,  0.000081,  0.000056,  0.000060,  0.000047,  0.000009,  0.000034,  0.000000,  0.000003,  0.000009,  0.000000,  0.000003,  0.000003,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000};
    
    double Data2015_50ns[60] = {0.000000,  0.000000,  18.000000,  28.000000,  55.000000,  102.000000,  274.000000,  480.000000,  711.000000,  1183.000000,  1523.000000,  1859.000000,  2229.000000,  2511.000000,  2605.000000,  2587.000000,  2243.000000,  2030.000000,  1761.000000,  1492.000000,  1129.000000,  870.000000,  632.000000,  402.000000,  285.000000,  200.000000,  136.000000,  98.000000,  65.000000,  28.000000,  20.000000,  18.000000,  5.000000,  1.000000,  2.000000,  1.000000,  2.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000};
    
    
    for( int i=0; i<60; ++i) {
      Data_.push_back(Data2015_50ns[i]);
      MC_.push_back(Asymp50ns[i]);
    }
    LumiWeights_ = edm::LumiReWeighting(MC_,Data_);
  } else {
    LumiWeights_ = edm::LumiReWeighting(pileupMC_, pileupData_, "pileup", "pileup");
  }
}

PileupWeightProducer::~PileupWeightProducer()
{}

// ------------ method called to produce the data  ------------
void
PileupWeightProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  using namespace edm;
  std::auto_ptr<std::vector<float> > pileupWeights( new std::vector<float> );

  Handle<std::vector< PileupSummaryInfo > >  PupInfo;
  iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);
  
  std::vector<PileupSummaryInfo>::const_iterator PVI;
  
  float Tnpv = -1;
  for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
    
    int BX = PVI->getBunchCrossing();
    
    if(BX == 0) { 
      Tnpv = PVI->getTrueNumInteractions();
      continue;
    }
  }
  
  double MyWeight = LumiWeights_.weight( Tnpv ); 
  pileupWeights->push_back( MyWeight );
  iEvent.put(pileupWeights, "pileupWeights"); 
}

// ------------ method called once each job just before starting event loop  ------------
void 
PileupWeightProducer::beginJob()
{}

// ------------ method called once each job just after ending the event loop  ------------
void 
PileupWeightProducer::endJob() 
{}

// ------------ method called when starting to processes a run  ------------
void 
PileupWeightProducer::beginRun(edm::Run&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a run  ------------
void 
PileupWeightProducer::endRun(edm::Run&, edm::EventSetup const&)
{}

// ------------ method called when starting to processes a luminosity block  ------------
void 
PileupWeightProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
PileupWeightProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PileupWeightProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PileupWeightProducer);
