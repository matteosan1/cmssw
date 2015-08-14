#ifndef SelectorByValueMap_h
#define SelectorByValueMap_h

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/Handle.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include <vector>
#include <string>

template <typename C>
struct CompatibleConfigurationType {
  typedef C type;
};
    
// "float" is not allowed, as it conflicts with "double"
template <>
struct CompatibleConfigurationType<float> {
  typedef double type;
};
  
template <typename T, typename C>
class SelectorByValueMap : public edm::EDProducer {
public:
  explicit SelectorByValueMap(edm::ParameterSet const & config);
  
private:
  typedef T candidate_type;
  typedef C selection_type;
  typedef typename CompatibleConfigurationType<selection_type>::type cut_type;
  typedef std::vector<candidate_type> candidateCollection;
  typedef edm::Ref<candidateCollection> candidateRef;
  typedef edm::RefVector<candidateCollection> candidateRefVector;

  void produce(edm::Event & event, edm::EventSetup const & setup);
    
  //edm::EDGetTokenT<candidateCollection>    token_inputs;
  edm::EDGetTokenT<candidateRefVector>    token_inputs;
  StringCutObjectSelector<T> candSelector;
  edm::EDGetTokenT<edm::ValueMap<selection_type> >  token_selection;
  cut_type m_cut;
  bool isGreaterThan_;
};
  
template <typename T, typename C>
SelectorByValueMap<T,C>::SelectorByValueMap(edm::ParameterSet const & config) :
  //token_inputs(consumes<candidateCollection>(config.getParameter<edm::InputTag>("input"))),
  token_inputs(consumes<candidateRefVector>(config.getParameter<edm::InputTag>("input"))),
  candSelector(config.getParameter<std::string>("cut")),
  token_selection(consumes<edm::ValueMap<selection_type> >(config.getParameter<edm::InputTag>("selection"))),
  m_cut(config.getParameter<cut_type>("id_cut")) {

  isGreaterThan_ = true;
  if (config.existsAs<edm::ParameterSet>("isGreaterThan")) 
    isGreaterThan_ = config.getParameter<bool>("isGreaterThan");

  produces<candidateRefVector>();
}
  
template <typename T, typename C>
void SelectorByValueMap<T, C>::produce(edm::Event & event, const edm::EventSetup & setup) {
  std::auto_ptr<candidateRefVector> candidates(new candidateRefVector());

  edm::Handle<candidateRefVector> h_inputs;
  event.getByToken(token_inputs, h_inputs);

  // read the selection map from the Event
  edm::Handle<edm::ValueMap<selection_type> > h_selection;
  event.getByToken(token_selection, h_selection);
  edm::ValueMap<selection_type> const & selectionMap = * h_selection;
  
  for (unsigned int i = 0; i < h_inputs->size(); ++i) {

    candidateRef ptr = (*h_inputs)[i];
    if (candSelector(*ptr)) {
      if (isGreaterThan_) {
	if (selectionMap[ptr] >= m_cut)
	  candidates->push_back(ptr);
      } else {
	if (selectionMap[ptr] < m_cut)
	  candidates->push_back(ptr);
      }
    }
  }
  
  event.put(candidates);
}

template<>
void SelectorByValueMap<pat::Electron, unsigned int>::produce(edm::Event & event, const edm::EventSetup & setup) {
  std::auto_ptr<candidateRefVector> candidates(new candidateRefVector());

  edm::Handle<candidateRefVector> h_inputs;
  event.getByToken(token_inputs, h_inputs);

  // read the selection map from the Event
  edm::Handle<edm::ValueMap<selection_type> > h_selection;
  event.getByToken(token_selection, h_selection);
  edm::ValueMap<selection_type> const & selectionMap = * h_selection;
  
  for (unsigned int i = 0; i < h_inputs->size(); ++i) {

    candidateRef ptr = (*h_inputs)[i];
    if (candSelector(*ptr)) {
      if ((selectionMap[ptr]&m_cut) == m_cut)
	candidates->push_back(ptr);
    }
  }
  
  event.put(candidates);
}

template<>
void SelectorByValueMap<pat::Photon, unsigned int>::produce(edm::Event & event, const edm::EventSetup & setup) {
  std::auto_ptr<candidateRefVector> candidates(new candidateRefVector());

  edm::Handle<candidateRefVector> h_inputs;
  event.getByToken(token_inputs, h_inputs);

  // read the selection map from the Event
  edm::Handle<edm::ValueMap<selection_type> > h_selection;
  event.getByToken(token_selection, h_selection);
  edm::ValueMap<selection_type> const & selectionMap = * h_selection;
  
  for (unsigned int i = 0; i < h_inputs->size(); ++i) {

    candidateRef ptr = (*h_inputs)[i];
    if (candSelector(*ptr)) {
      if ((selectionMap[ptr]&m_cut) == m_cut)
	candidates->push_back(ptr);
    }
  }
  
  event.put(candidates);
}
#endif

typedef SelectorByValueMap<pat::Electron, bool> PatElectronSelectorByValueMap;
DEFINE_FWK_MODULE(PatElectronSelectorByValueMap);

typedef SelectorByValueMap<pat::Photon, bool> PatPhotonSelectorByValueMap;
DEFINE_FWK_MODULE(PatPhotonSelectorByValueMap);

typedef SelectorByValueMap<pat::Electron, unsigned int> PatElectronSelectorByBitMap;
DEFINE_FWK_MODULE(PatElectronSelectorByBitMap);

typedef SelectorByValueMap<pat::Photon, unsigned int> PatPhotonSelectorByBitMap;
DEFINE_FWK_MODULE(PatPhotonSelectorByBitMap);
