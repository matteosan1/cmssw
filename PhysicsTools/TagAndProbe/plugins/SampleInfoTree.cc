#include "PhysicsTools/TagAndProbe/plugins/SampleInfoTree.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

tnp::SampleInfoTree::SampleInfoTree(const edm::ParameterSet& iConfig) {
  // make trees as requested
  edm::Service<TFileService> fs;
  addTree_ = fs->make<TTree>("sampleInfo", "sampleInfo");
  hNvtx = fs->make<TH1F>("nvtx", "nvtx", 100, 0., 100.);

  sumWeight_ = 0.0;
  nEvents_ = 0;
 
  weightSrcToken_ = consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genInfo"));
  recVtxsToken_   = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"));
  
  addTree_->Branch("sumWeight", &sumWeight_, "sumWeight/D");
  addTree_->Branch("nEvents", &nEvents_, "nEvents/D");
}

tnp::SampleInfoTree::~SampleInfoTree() {}

void tnp::SampleInfoTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  edm::Handle<GenEventInfoProduct> weight;
  iEvent.getByToken(weightSrcToken_, weight);
  
  if (iEvent.eventAuxiliary().isRealData())
    sumWeight_ += 1.;
  else
    sumWeight_ += weight->weight();

  nEvents_ += 1.;

  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByToken(recVtxsToken_,recVtxs); 

  if (iEvent.eventAuxiliary().isRealData())
    hNvtx->Fill(recVtxs->size(), 1);
  else
    hNvtx->Fill(recVtxs->size(), weight->weight());
}

void tnp::SampleInfoTree::endJob() {
  addTree_->Fill();
}

DEFINE_FWK_MODULE(tnp::SampleInfoTree);
