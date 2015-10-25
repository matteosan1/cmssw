#include "MiniAODTriggerCandProducer.h"
  
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

template <>
bool MiniAODTriggerCandProducer<pat::Electron>::onlineOfflineMatching(const edm::TriggerNames & triggerNames, 
								      edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects, 
								      pat::ElectronRef ref, std::string filterLabel, float dRmin) {
  
  for (pat::TriggerObjectStandAlone obj : *triggerObjects) { 
    obj.unpackPathNames(triggerNames); 
    if (obj.hasFilterLabel(filterLabel)) {
      float dR = deltaR(ref->superCluster()->position(), obj.p4());
      if (dR < dRmin)
	return true;
    }
  }

  return false;
}

template <>
bool MiniAODTriggerCandProducer<pat::Photon>::onlineOfflineMatching(const edm::TriggerNames & triggerNames, 
								    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects, 
								    pat::PhotonRef ref, std::string filterLabel, float dRmin) {
  
  for (pat::TriggerObjectStandAlone obj : *triggerObjects) { 
    obj.unpackPathNames(triggerNames); 
    if (obj.hasFilterLabel(filterLabel)) {
      float dR = deltaR(ref->superCluster()->position(), obj.p4());
      if (dR < dRmin)
	return true;
    }
  }

  return false;
}

template <>
bool MiniAODTriggerCandProducer<reco::RecoEcalCandidate>::onlineOfflineMatching(const edm::TriggerNames & triggerNames, 
										edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects, 
										edm::Ref<std::vector<reco::RecoEcalCandidate>> ref, 
										std::string filterLabel, float dRmin) {

  for (pat::TriggerObjectStandAlone obj : *triggerObjects) { 
    obj.unpackPathNames(triggerNames); 
    if (obj.hasFilterLabel(filterLabel)) {
      float dR = deltaR(ref->superCluster()->position(), obj.p4());
      if (dR < dRmin)
	return true;
    }
  }

  return false;
}

typedef MiniAODTriggerCandProducer<pat::Electron> PatElectronTriggerCandProducer;
DEFINE_FWK_MODULE(PatElectronTriggerCandProducer);

typedef MiniAODTriggerCandProducer<pat::Photon> PatPhotonTriggerCandProducer;
DEFINE_FWK_MODULE(PatPhotonTriggerCandProducer);

typedef MiniAODTriggerCandProducer<reco::RecoEcalCandidate> RecoEcalCandidateTriggerCandProducer;
DEFINE_FWK_MODULE(RecoEcalCandidateTriggerCandProducer);
