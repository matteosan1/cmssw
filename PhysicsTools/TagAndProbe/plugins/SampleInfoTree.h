#ifndef PhysicsTools_TagAndProbe_SampleInfoTree_h
#define PhysicsTools_TagAndProbe_SampleInfoTree_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include <TTree.h>
#include <TH1F.h>

namespace tnp {

  class SampleInfoTree : public edm::EDAnalyzer, boost::noncopyable {
  public:
    explicit SampleInfoTree(const edm::ParameterSet& config);
    ~SampleInfoTree();

  private:
    //void beginJob();
    void analyze(const edm::Event&, const edm::EventSetup&);
    //void endLuminosityBlock(const edm::LuminosityBlock &, const edm::EventSetup &);
    void endJob();
    
    edm::EDGetTokenT<GenEventInfoProduct> weightSrcToken_;
    edm::EDGetTokenT<reco::VertexCollection> recVtxsToken_;
    
    mutable TTree * addTree_;
    mutable double sumWeight_;
    mutable double nEvents_; 
    mutable TH1F* hNvtx;
  };
}

#endif
