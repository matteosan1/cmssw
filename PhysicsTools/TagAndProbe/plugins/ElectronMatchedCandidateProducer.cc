#include "PhysicsTools/TagAndProbe/interface/ElectronMatchedCandidateProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Math/interface/deltaR.h"


ElectronMatchedCandidateProducer::ElectronMatchedCandidateProducer(const edm::ParameterSet &params):
  electronCollectionToken_(consumes<edm::RefVector<pat::ElectronCollection> >(params.getUntrackedParameter<edm::InputTag>("ReferenceElectronCollection"))),
  scCollectionToken_(consumes<reco::RecoEcalCandidateCollection> (params.getParameter<edm::InputTag>("src"))),
  candSelector_(params.getParameter<std::string>("cut")),
  delRMatchingCut_(params.getUntrackedParameter<double>("deltaR", 0.30)) {
  
  //const edm::InputTag allelectrons("gsfElectrons");
  //  electronCollectionToken_ =
  //  consumes<edm::View<reco::GsfElectron> >(params.getUntrackedParameter<edm::InputTag>("ReferenceElectronCollection",
  //											allelectrons));
  
  
  //scCollectionToken_ =
  //  consumes<edm::View<reco::Candidate> >(params.getParameter<edm::InputTag>("src"));

  
  //produces< edm::PtrVector<reco::Candidate> >();
  //produces< edm::RefToBaseVector<reco::Candidate> >();
  produces<edm::RefVector<reco::RecoEcalCandidateCollection> >();
}

ElectronMatchedCandidateProducer::~ElectronMatchedCandidateProducer()
{}

void ElectronMatchedCandidateProducer::produce(edm::Event &event,
					       const edm::EventSetup &eventSetup) {

  // Create the output collection
  //std::auto_ptr< edm::RefToBaseVector<reco::Candidate> >
  //  outColRef( new edm::RefToBaseVector<reco::Candidate> );
  //std::auto_ptr< edm::PtrVector<reco::Candidate> >
  //  outColPtr( new edm::PtrVector<reco::Candidate> );
  std::auto_ptr<edm::RefVector<reco::RecoEcalCandidateCollection> > outCol (new edm::RefVector<reco::RecoEcalCandidateCollection>);
  
  // Read electrons
  //  edm::Handle<edm::View<reco::GsfElectron> > electrons;
  //event.getByToken(electronCollectionToken_, electrons);

  //Read candidates
  //edm::Handle<edm::View<reco::Candidate> > recoCandColl;
  //event.getByToken( scCollectionToken_ , recoCandColl);

  // Read electrons
  edm::Handle<edm::RefVector<pat::ElectronCollection> > electrons;
  event.getByToken(electronCollectionToken_, electrons);

  //Read candidates
  edm::Handle<reco::RecoEcalCandidateCollection> recoCandColl;
  event.getByToken(scCollectionToken_ , recoCandColl);

  //unsigned int counter=0;

  // Loop over candidates
  //for(edm::View<reco::Candidate>::const_iterator scIt = recoCandColl->begin();
  //    scIt != recoCandColl->end(); ++scIt, ++counter){
  for (size_t sc=0; sc<recoCandColl->size(); sc++) {
    
    reco::RecoEcalCandidateRef ref(recoCandColl, sc);
    // Now loop over electrons

    if (candSelector_(*ref)) {
      //for(edm::reco::GsfElectronCollection::const_iterator elec = electrons->begin();
      //elec != electrons->end(); ++elec) {
      for (size_t elec=0; elec<electrons->size(); elec++) {
	
	reco::SuperClusterRef eSC = (*electrons)[elec]->superCluster();
	
	double dRval = reco::deltaR((float)ref->eta(), (float)ref->phi(),
				    eSC->eta(), eSC->phi());
	
	if(dRval < delRMatchingCut_) {
	  
	  outCol->push_back(ref);
	  //outColRef->push_back( recoCandColl->refAt(counter) );
	  //outColPtr->push_back( recoCandColl->ptrAt(counter) );
	  
	} // end if loop
      } // end electron loop
    } // end candidate loop
  }
  //event.put(outColRef);
  //event.put(outColPtr);
  event.put(outCol);
}

void ElectronMatchedCandidateProducer::beginJob() 
{}

void ElectronMatchedCandidateProducer::endJob() 
{}

DEFINE_FWK_MODULE(ElectronMatchedCandidateProducer);

