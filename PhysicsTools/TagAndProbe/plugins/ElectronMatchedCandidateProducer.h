#ifndef PhysicsTools_TagAndProbe_ElectronMatchedCandidateProducer_h
#define PhysicsTools_TagAndProbe_ElectronMatchedCandidateProducer_h

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Electron.h"

#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

class ElectronMatchedCandidateProducer : public edm::EDProducer {
 public:
  explicit ElectronMatchedCandidateProducer(const edm::ParameterSet&);
  ~ElectronMatchedCandidateProducer();

 private:
  virtual void produce(edm::Event&, const edm::EventSetup&) override;

  edm::EDGetTokenT<edm::RefVector<pat::ElectronCollection> > electronCollectionToken_;
  edm::EDGetTokenT<reco::RecoEcalCandidateCollection> scCollectionToken_;
  StringCutObjectSelector<reco::RecoEcalCandidate> candSelector_;
};

#endif
