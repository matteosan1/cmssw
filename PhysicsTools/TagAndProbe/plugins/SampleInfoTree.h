#ifndef PhysicsTools_TagAndProbe_SampleInfoTree_h
#define PhysicsTools_TagAndProbe_SampleInfoTree_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "DataFormats/Common/interface/MergeableDouble.h"

#include "DataFormats/Common/interface/Handle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include <TTree.h>
#include <TH1F.h>

namespace tnp {
  class SampleInfoTree : public edm::one::EDProducer<edm::one::WatchLuminosityBlocks, 
    edm::EndLuminosityBlockProducer> {

  public:
    explicit SampleInfoTree(const edm::ParameterSet& config);
    ~SampleInfoTree() {};

  private:
    virtual void beginLuminosityBlock(const edm::LuminosityBlock &, const edm::EventSetup&);
    virtual void endLuminosityBlock(edm::LuminosityBlock const&, const edm::EventSetup&) override;
    virtual void endLuminosityBlockProduce(edm::LuminosityBlock &, const edm::EventSetup&) override;
    
    //void beginJob(); 
    virtual void produce(edm::Event &, const edm::EventSetup&) override;
    //void produce(const edm::Event&, const edm::EventSetup&);
    //void analyze(const edm::Event&, const edm::EventSetup&);
    //void endLuminosityBlock(const edm::LuminosityBlock &, const edm::EventSetup &);
    void endJob();
    
    edm::EDGetTokenT<GenEventInfoProduct> weightSrcToken_;
    //edm::EDGetTokenT<reco::VertexCollection> recVtxsToken_;
    
    mutable TTree * addTree_;
    mutable double sumWeight_;
    mutable double nEvents_; 

    mutable double totSumWeight_;
    mutable double totNEvents_; 
  };
}

#endif
