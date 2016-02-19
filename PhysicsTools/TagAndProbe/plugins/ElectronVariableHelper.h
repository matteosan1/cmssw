#ifndef _ELECTRONVARIABLEHELPER_H
#define _ELECTRONVARIABLEHELPER_H

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"

#include "DataFormats/Math/interface/deltaR.h"

template <class T>
class ElectronVariableHelper : public edm::EDProducer {
 public:
  explicit ElectronVariableHelper(const edm::ParameterSet & iConfig);
  virtual ~ElectronVariableHelper() ;
  
  virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) override;
  
private:
  edm::EDGetTokenT<std::vector<T> > probesToken_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<l1extra::L1EmParticleCollection> l1NonIsoToken_;
  edm::EDGetTokenT<l1extra::L1EmParticleCollection> l1IsoToken_;
};

template<class T>
ElectronVariableHelper<T>::ElectronVariableHelper(const edm::ParameterSet & iConfig) :
  probesToken_(consumes<std::vector<T> >(iConfig.getParameter<edm::InputTag>("probes"))),
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"))),
  l1NonIsoToken_(consumes<l1extra::L1EmParticleCollection>(edm::InputTag("l1extraParticles:NonIsolated"))),
  l1IsoToken_(consumes<l1extra::L1EmParticleCollection>(edm::InputTag("l1extraParticles:Isolated"))) {

  produces<edm::ValueMap<float> >("dz");
  produces<edm::ValueMap<float> >("dxy");
  produces<edm::ValueMap<float> >("missinghits");
  produces<edm::ValueMap<float> >("l1et");
}

template<class T>
ElectronVariableHelper<T>::~ElectronVariableHelper()
{}

template<class T>
void ElectronVariableHelper<T>::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {

  // read input
  edm::Handle<std::vector<T> > probes;
  edm::Handle<reco::VertexCollection> vtxH;
  edm::Handle<std::vector<l1extra::L1EmParticle> > l1NonIsoH, l1IsoH;
  
  iEvent.getByToken(probesToken_,  probes);
  iEvent.getByToken(vtxToken_, vtxH);
  const reco::VertexRef vtx(vtxH, 0);
  iEvent.getByToken(l1NonIsoToken_, l1NonIsoH);
  const l1extra::L1EmParticleCollection* l1NonIso = l1NonIsoH.product();
  iEvent.getByToken(l1IsoToken_, l1IsoH);
  const l1extra::L1EmParticleCollection* l1Iso = l1IsoH.product();

  // prepare vector for output
  std::vector<float> dzValues;
  std::vector<float> dxyValues;
  std::vector<float> mhValues;
  std::vector<float> l1Values;

  typename std::vector<T>::const_iterator probe, endprobes = probes->end();

  for (probe = probes->begin(); probe != endprobes; ++probe) {
    
    dzValues.push_back(probe->gsfTrack()->dz(vtx->position()));
    dxyValues.push_back(probe->gsfTrack()->dxy(vtx->position()));
    mhValues.push_back(float(probe->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS)));
    
    float l1et = -1;
    float dRmin = 99999;
    for (l1extra::L1EmParticleCollection::const_iterator l1It = l1Iso->begin(); l1It != l1Iso->end(); l1It++) {
      float dR = deltaR(*l1It, *(probe->superCluster()));
      if (dR < dRmin) {
	dRmin = dR;
	l1et = l1It->et();
      }
    }
     
    if (l1et < 0) {
      for (l1extra::L1EmParticleCollection::const_iterator l1It = l1NonIso->begin(); l1It != l1NonIso->end(); l1It++) {
	float dR = deltaR(*l1It, *(probe->superCluster()));
	if (dR < dRmin) {
	  dRmin = dR;
	  l1et = l1It->et();
	}
      }
    }
    
    l1Values.push_back(l1et);
  }

  
  // convert into ValueMap and store
  std::auto_ptr<edm::ValueMap<float> > dzValMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler dzFiller(*dzValMap);
  dzFiller.insert(probes, dzValues.begin(), dzValues.end());
  dzFiller.fill();
  iEvent.put(dzValMap, "dz");

  std::auto_ptr<edm::ValueMap<float> > dxyValMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler dxyFiller(*dxyValMap);
  dxyFiller.insert(probes, dxyValues.begin(), dxyValues.end());
  dxyFiller.fill();
  iEvent.put(dxyValMap, "dxy");

  std::auto_ptr<edm::ValueMap<float> > mhValMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler mhFiller(*mhValMap);
  mhFiller.insert(probes, mhValues.begin(), mhValues.end());
  mhFiller.fill();
  iEvent.put(mhValMap, "missinghits");

  std::auto_ptr<edm::ValueMap<float> > l1ValMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler l1Filler(*l1ValMap);
  l1Filler.insert(probes, l1Values.begin(), l1Values.end());
  l1Filler.fill();
  iEvent.put(l1ValMap, "l1et");
}

#endif
