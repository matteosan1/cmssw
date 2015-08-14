#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

class ElectronVariableHelper : public edm::EDProducer {
public:
  explicit ElectronVariableHelper(const edm::ParameterSet & iConfig);
  virtual ~ElectronVariableHelper() ;
  
  virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) override;
  
private:
  edm::EDGetTokenT<std::vector<pat::Electron> > probesToken_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
};

ElectronVariableHelper::ElectronVariableHelper(const edm::ParameterSet & iConfig) :
  probesToken_(consumes<std::vector<pat::Electron> >(iConfig.getParameter<edm::InputTag>("probes"))),
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"))) {

  produces<edm::ValueMap<float> >("dz");
  produces<edm::ValueMap<float> >("dxy");
  produces<edm::ValueMap<float> >("missinghits");
}

ElectronVariableHelper::~ElectronVariableHelper()
{}

void ElectronVariableHelper::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {

  // read input
  edm::Handle<std::vector<pat::Electron> > probes;
  edm::Handle<reco::VertexCollection> vtxH;

  iEvent.getByToken(probesToken_,  probes);
  iEvent.getByToken(vtxToken_, vtxH);
  const reco::VertexRef vtx(vtxH, 0);

  // prepare vector for output
  std::vector<float> dzValues;
  std::vector<float> dxyValues;
  std::vector<float> mhValues;

  std::vector<pat::Electron>::const_iterator probe, endprobes = probes->end();
  //const std::vector<pat::Electron>* electronCollection = elHandle.product();
  //reco::GsfElectronCollection::const_iterator eleIt = electronCollection->begin();

  for (probe = probes->begin(); probe != endprobes; ++probe) {
    
    dzValues.push_back(probe->gsfTrack()->dz(vtx->position()));
    dxyValues.push_back(probe->gsfTrack()->dxy(vtx->position()));
    mhValues.push_back(float(probe->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS)));
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
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ElectronVariableHelper);
