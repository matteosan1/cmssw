#include "PhysicsTools/TagAndProbe/interface/ElectronMatchedCandidateProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/MakerMacros.h"

//#include "DataFormats/Math/interface/deltaR.h"


ElectronMatchedCandidateProducer::ElectronMatchedCandidateProducer(const edm::ParameterSet &params):
  electronCollectionToken_(consumes<edm::RefVector<pat::ElectronCollection> >(params.getUntrackedParameter<edm::InputTag>("ReferenceElectronCollection"))),
  scCollectionToken_(consumes<reco::RecoEcalCandidateCollection> (params.getParameter<edm::InputTag>("src"))),
  candSelector_(params.getParameter<std::string>("cut")) {
  //delRMatchingCut_(params.getUntrackedParameter<double>("deltaR", 0.30)) {
  produces<edm::RefVector<reco::RecoEcalCandidateCollection> >();
}

ElectronMatchedCandidateProducer::~ElectronMatchedCandidateProducer()
{}

void ElectronMatchedCandidateProducer::produce(edm::Event &event,
					       const edm::EventSetup &eventSetup) {

  std::auto_ptr<edm::RefVector<reco::RecoEcalCandidateCollection> > outCol (new edm::RefVector<reco::RecoEcalCandidateCollection>);

  // Read electrons
  edm::Handle<edm::RefVector<pat::ElectronCollection> > electrons;
  event.getByToken(electronCollectionToken_, electrons);

  //Read candidates
  edm::Handle<reco::RecoEcalCandidateCollection> recoCandColl;
  event.getByToken(scCollectionToken_ , recoCandColl);

  for (size_t sc=0; sc<recoCandColl->size(); sc++) {
    
    reco::RecoEcalCandidateRef ref(recoCandColl, sc);
    if (candSelector_(*ref)) {
      for (size_t elec=0; elec<electrons->size(); elec++) {
	if ((*electrons)[elec]->superCluster() == ref->superCluster())
	  outCol->push_back(ref);

	//reco::SuperClusterRef eSC = (*electrons)[elec]->superCluster();
	//
	//double dRval = reco::deltaR((float)ref->eta(), (float)ref->phi(),
	//			    eSC->eta(), eSC->phi());
	//
	//if(dRval < delRMatchingCut_) {
	//  
	//  outCol->push_back(ref);
	//  //outColRef->push_back( recoCandColl->refAt(counter) );
	//  //outColPtr->push_back( recoCandColl->ptrAt(counter) );
	//  
	//} // end if loop
      } 
    } 
  }

  event.put(outCol);
}

void ElectronMatchedCandidateProducer::beginJob() 
{}

void ElectronMatchedCandidateProducer::endJob() 
{}

DEFINE_FWK_MODULE(ElectronMatchedCandidateProducer);

