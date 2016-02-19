#include "PhysicsTools/TagAndProbe/plugins/ElectronMatchedCandidateProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/MakerMacros.h"

ElectronMatchedCandidateProducer::ElectronMatchedCandidateProducer(const edm::ParameterSet &params):
  electronCollectionToken_(consumes<edm::RefVector<pat::ElectronCollection> >(params.getUntrackedParameter<edm::InputTag>("ReferenceElectronCollection"))),
  scCollectionToken_(consumes<reco::RecoEcalCandidateCollection> (params.getParameter<edm::InputTag>("src"))),
  candSelector_(params.getParameter<std::string>("cut")) 
{
  produces<edm::RefVector<reco::RecoEcalCandidateCollection> >("superclusters");
  produces<edm::RefVector<std::vector<pat::Electron> > >("electrons");
}

ElectronMatchedCandidateProducer::~ElectronMatchedCandidateProducer()
{}

void ElectronMatchedCandidateProducer::produce(edm::Event &event,
					       const edm::EventSetup &eventSetup) {

  std::auto_ptr<edm::RefVector<reco::RecoEcalCandidateCollection> > outCol (new edm::RefVector<reco::RecoEcalCandidateCollection>);
  std::auto_ptr<edm::RefVector<std::vector<pat::Electron> > > outCol2 (new edm::RefVector<std::vector<pat::Electron> >);

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
	if ((*electrons)[elec]->superCluster() == ref->superCluster()) {
	  outCol->push_back(ref);
	  outCol2->push_back((*electrons)[elec]);
	}
      } 
    } 
  }

  event.put(outCol, "superclusters");
  event.put(outCol2, "electrons");
}

DEFINE_FWK_MODULE(ElectronMatchedCandidateProducer);

