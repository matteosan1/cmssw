#include "PhysicsTools/TagAndProbe/plugins/SampleInfoTree.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <iostream>

tnp::SampleInfoTree::SampleInfoTree(const edm::ParameterSet& iConfig) : 
  weightSrcToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genInfo")))  {
  
  // make trees as requested
  edm::Service<TFileService> fs;
  addTree_ = fs->make<TTree>("sampleInfo", "sampleInfo");
    
  addTree_->Branch("sumWeight", &totSumWeight_, "sumWeight/D");
  addTree_->Branch("nEvents", &totNEvents_, "nEvents/D");

  totSumWeight_ = 0.0;
  totNEvents_ = 0.0;

  produces<edm::MergeableDouble, edm::InLumi> ("totalGenWeight"); 
  produces<edm::MergeableDouble, edm::InLumi> ("totalEvent");
}
    
void tnp::SampleInfoTree::beginLuminosityBlock(const edm::LuminosityBlock & theLuminosityBlock, const edm::EventSetup & theSetup) {
  sumWeight_ = 0.0;
  nEvents_ = 0.0;
}

void tnp::SampleInfoTree::endLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iSetup) 
{}

void tnp::SampleInfoTree::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<GenEventInfoProduct> genInfo;
  
  if(!iEvent.isRealData()) {
    iEvent.getByToken(weightSrcToken_, genInfo);
    
    const auto & weights = genInfo->weights(); 
    if(!weights.empty()) {
      sumWeight_ += weights[0];
      totSumWeight_ += weights[0];
    }
  } else {
    sumWeight_ += 1.;
    totSumWeight_ += 1.;
  }

  totNEvents_ += 1.0;
  nEvents_ += 1.;
}  

void tnp::SampleInfoTree::endJob() {
  //std::cout << totSumWeight_ << " " << totNEvents_ << std::endl;
  addTree_->Fill();
}

void tnp::SampleInfoTree::endLuminosityBlockProduce(edm::LuminosityBlock & theLuminosityBlock, const edm::EventSetup & theSetup) {
  
  std::auto_ptr<edm::MergeableDouble> numWeightsPtr(new edm::MergeableDouble);
  numWeightsPtr->value = sumWeight_;
  theLuminosityBlock.put(numWeightsPtr, "totalGenWeight");
  
  std::auto_ptr<edm::MergeableDouble> numEventsPtr(new edm::MergeableDouble);
  numEventsPtr->value = nEvents_;
  theLuminosityBlock.put(numEventsPtr, "totalEvent");
}

DEFINE_FWK_MODULE(tnp::SampleInfoTree);
