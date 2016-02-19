// -*- C++ -*-
//
// Package:    TagProbeFitTreeProducer
// Class:      TagProbeFitTreeProducer
//
/**\class TagProbeFitTreeProducer TagProbeFitTreeProducer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Sep 15 09:45
//         Created:  Mon Sep 15 09:49:08 CEST 2008
//
//


// system include files
#include <memory>
#include <ctype.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/MergeableDouble.h"

#include "PhysicsTools/TagAndProbe/interface/TPTreeFiller.h"
#include "PhysicsTools/TagAndProbe/interface/TagProbePairMaker.h"

#include <set>

#include "FWCore/ParameterSet/interface/Registry.h"

#include <DataFormats/Math/interface/deltaR.h>

class TagProbeFitTreeProducer : public edm::EDAnalyzer {
public:
  explicit TagProbeFitTreeProducer(const edm::ParameterSet&);
  ~TagProbeFitTreeProducer();
    
private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override ;
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void beginLuminosityBlock(const edm::LuminosityBlock &, const edm::EventSetup&) override;

  //---- MC truth information
  bool isMC_;
  /// Token for an edm::Association<reco::GenParticle> from tags & probes to MC truth
  //edm::EDGetTokenT<edm::Association<std::vector<reco::GenParticle> > > tagMatchesToken_, probeMatchesToken_;
  edm::EDGetTokenT<std::vector<reco::GenParticle> > genParticleToken_;

  /// Possible pdgids for the mother. If empty, any truth-matched mu will be considered good
  std::set<int32_t> motherPdgId_;
  bool useTauDecays_;
  int pdgId_;
  bool checkCharge_;

  /// Return true if ref is not null and has an ancestor with pdgId inside 'motherPdgId_'
  //bool checkMother(const reco::GenParticleRef &ref) const ;
  std::pair<const reco::GenParticleRef, bool> checkStatus(const reco::CandidateBaseRef& ref, edm::Handle<std::vector<reco::GenParticle> > genP);
  
  //---- Unbiased MC truth information
  /// Do we have to compute this
  //bool makeMCUnbiasTree_;
  /// Check mother pdgId in unbiased inefficiency measurement
  //bool checkMotherInUnbiasEff_;
  
  /// InputTag to the collection of all probes
  edm::EDGetTokenT<reco::CandidateView> allProbesToken_;
  edm::EDGetTokenT<edm::MergeableDouble> totGenWeightToken_;
  edm::EDGetTokenT<edm::MergeableDouble> totEventsToken_;

  /// The object that produces pairs of tags and probes, making any arbitration needed
  tnp::TagProbePairMaker tagProbePairMaker_;

  /// The object that actually computes variables and fills the tree for T&P
  std::auto_ptr<tnp::TPTreeFiller> treeFiller_;
  std::auto_ptr<tnp::BaseTreeFiller> mcUnbiasFiller_;
  std::auto_ptr<tnp::BaseTreeFiller> oldTagFiller_;
  std::auto_ptr<tnp::BaseTreeFiller> tagFiller_;
  std::auto_ptr<tnp::BaseTreeFiller> pairFiller_;
  std::auto_ptr<tnp::BaseTreeFiller> mcFiller_;

  double totGenWeight_;
  double totEvents_;
};

TagProbeFitTreeProducer::TagProbeFitTreeProducer(const edm::ParameterSet& iConfig) :
    isMC_(iConfig.getParameter<bool>("isMC")),
    //makeMCUnbiasTree_(isMC_ ? iConfig.getParameter<bool>("makeMCUnbiasTree") : false),
    //checkMotherInUnbiasEff_(makeMCUnbiasTree_ ? iConfig.getParameter<bool>("checkMotherInUnbiasEff") : false),
    tagProbePairMaker_(iConfig, consumesCollector()),
    treeFiller_(new tnp::TPTreeFiller(iConfig, consumesCollector())),
    oldTagFiller_((iConfig.existsAs<bool>("fillTagTree") && iConfig.getParameter<bool>("fillTagTree")) ? new tnp::BaseTreeFiller("tag_tree",iConfig, consumesCollector()) : 0)
{
  if (isMC_) {
    // For mc efficiency we need the MC matches for tags & probes
    //tagMatchesToken_ = consumes<edm::Association<std::vector<reco::GenParticle> > >(iConfig.getParameter<edm::InputTag>("tagMatches"));
    //probeMatchesToken_ = consumes<edm::Association<std::vector<reco::GenParticle> > >(iConfig.getParameter<edm::InputTag>("probeMatches"));
    totGenWeightToken_ = consumes<edm::MergeableDouble, edm::InLumi>(edm::InputTag("sampleInfo:totalGenWeight"));
    totEventsToken_ = consumes<edm::MergeableDouble, edm::InLumi>(edm::InputTag("sampleInfo:totalEvent"));

    genParticleToken_ = consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"));
    useTauDecays_ = iConfig.getParameter<bool>("useTauDecays");
    pdgId_ = iConfig.getParameter<int>("pdgId");
    checkCharge_ = iConfig.getParameter<bool>("checkCharge");

    //.. and the pdgids of the possible mothers
    if (iConfig.existsAs<int32_t>("motherPdgId")) {
      motherPdgId_.insert(iConfig.getParameter<int32_t>("motherPdgId"));
    } else {
      std::vector<int32_t> motherIds = iConfig.getParameter<std::vector<int32_t> >("motherPdgId");
      motherPdgId_.insert(motherIds.begin(), motherIds.end());
    }
    
    // For unbiased efficiency we also need the collection of all probes
    //if (makeMCUnbiasTree_) {
    //  allProbesToken_ = consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("allProbes"));
    //  mcUnbiasFiller_.reset(new tnp::BaseTreeFiller("mcUnbias_tree",iConfig, consumesCollector()));
    //} 

    totGenWeight_ = 0.;
    totEvents_ = 0.;
  }
  
  edm::ParameterSet tagPSet;
  if (iConfig.existsAs<edm::ParameterSet>("tagVariables")) tagPSet.addParameter<edm::ParameterSet>("variables", iConfig.getParameter<edm::ParameterSet>("tagVariables"));
  if (iConfig.existsAs<edm::ParameterSet>("tagFlags"    )) tagPSet.addParameter<edm::ParameterSet>("flags",     iConfig.getParameter<edm::ParameterSet>("tagFlags"));
  if (!tagPSet.empty()) { tagFiller_.reset(new tnp::BaseTreeFiller(*treeFiller_, tagPSet, consumesCollector(), "tag_")); }
  edm::ParameterSet mcPSet;
  if (iConfig.existsAs<edm::ParameterSet>("mcVariables")) mcPSet.addParameter<edm::ParameterSet>("variables", iConfig.getParameter<edm::ParameterSet>("mcVariables"));
  if (iConfig.existsAs<edm::ParameterSet>("mcFlags"    )) mcPSet.addParameter<edm::ParameterSet>("flags",     iConfig.getParameter<edm::ParameterSet>("mcFlags"));
  if (!mcPSet.empty()) { mcFiller_.reset(new tnp::BaseTreeFiller(*treeFiller_, mcPSet, consumesCollector(), "mc_")); }
  edm::ParameterSet pairPSet;
  if (iConfig.existsAs<edm::ParameterSet>("pairVariables")) pairPSet.addParameter<edm::ParameterSet>("variables", iConfig.getParameter<edm::ParameterSet>("pairVariables"));
  if (iConfig.existsAs<edm::ParameterSet>("pairFlags"    )) pairPSet.addParameter<edm::ParameterSet>("flags",     iConfig.getParameter<edm::ParameterSet>("pairFlags"));
  if (!pairPSet.empty()) { pairFiller_.reset(new tnp::BaseTreeFiller(*treeFiller_, pairPSet, consumesCollector(), "pair_")); }
}


TagProbeFitTreeProducer::~TagProbeFitTreeProducer()
{}

// ------------ method called to for each event  ------------
void TagProbeFitTreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm; using namespace std;
  Handle<reco::CandidateView> src, allProbes;
  //Handle<Association<vector<reco::GenParticle> > > tagMatches, probeMatches;
  Handle<vector<reco::GenParticle> > genParticlesH;

  treeFiller_->init(iEvent); // read out info from the event if needed (external vars, list of passing probes, ...)
  if (oldTagFiller_.get()) oldTagFiller_->init(iEvent);
  if (tagFiller_.get())    tagFiller_->initPerObject(iEvent);
  if (pairFiller_.get())   pairFiller_->initPerObject(iEvent);
  if (mcFiller_.get())     mcFiller_->initPerObject(iEvent);

  // get the list of (tag+probe) pairs, performing arbitration
  tnp::TagProbePairs pairs = tagProbePairMaker_.run(iEvent);
  
  // on mc we want to load also the MC match info
  if (isMC_) {
    iEvent.getByToken(genParticleToken_, genParticlesH);
    //Ievent.getByToken(tagMatchesToken_, tagMatches);
    //iEvent.getByToken(probeMatchesToken_, probeMatches);
  }

  // loop on them to fill the tree
  for (tnp::TagProbePairs::const_iterator it = pairs.begin(), ed = pairs.end(); it != ed; ++it) {
    // on mc, fill mc info (on non-mc, let it to 'true', the treeFiller will ignore it anyway
    bool mcTrue = false;

    if (isMC_) {
      //reco::GenParticleRef mtag = (*tagMatches)[it->tag], mprobe = (*probeMatches)[it->probe];
      std::pair<reco::GenParticleRef,bool> mtag = checkStatus(it->tag, genParticlesH);
      std::pair<reco::GenParticleRef,bool> mprobe = checkStatus(it->probe, genParticlesH);
      //mcTrue = checkMother(mtag) && checkMother(mprobe);
      //std::cout << mtag.second  << " " << mprobe.second << std::endl;
      mcTrue = mtag.second && mprobe.second;
      if (mcTrue && mcFiller_.get()) mcFiller_->fill(reco::CandidateBaseRef(mprobe.first));
    }

    // fill in the variables for this t+p pair
    if (tagFiller_.get())    tagFiller_->fill(it->tag);
    if (oldTagFiller_.get()) oldTagFiller_->fill(it->tag);
    if (pairFiller_.get())   pairFiller_->fill(it->pair);
    treeFiller_->fill(it->probe, it->mass, mcTrue);
  }
  
  //if (isMC_ && makeMCUnbiasTree_) {
  //  // read full collection of probes
  //  iEvent.getByToken(allProbesToken_, allProbes);
  //  // init the tree filler
  //  mcUnbiasFiller_->init(iEvent);
  //  // loop on probes
  //  for (size_t i = 0, n = allProbes->size(); i < n; ++i) {
  //    const reco::CandidateBaseRef & probe = allProbes->refAt(i);
  //    // check mc match, and possibly mother match
  //    reco::GenParticleRef probeMatch = (*probeMatches)[probe];
  //    bool probeOk = checkMotherInUnbiasEff_ ? checkMother(probeMatch) : probeMatch.isNonnull();
  //    // fill the tree only for good ones
  //    if (probeOk) mcUnbiasFiller_->fill(probe);
  //  }
  //}  
}

std::pair<const reco::GenParticleRef, bool> TagProbeFitTreeProducer::checkStatus(const reco::CandidateBaseRef& ref, edm::Handle<std::vector<reco::GenParticle> > genParticlesH) {
  
  if (ref.isNonnull()) { 
    for (size_t gp=0; gp<genParticlesH->size(); gp++) {
      const reco::GenParticleRef p(genParticlesH, gp);
      if ((abs(p->pdgId()) == pdgId_ and !checkCharge_) or
	  (p->pdgId() == pdgId_ and checkCharge_)) {

	if (p->isPromptFinalState() or (p->isDirectPromptTauDecayProductFinalState() && useTauDecays_) or
	    (p->status() == 23)) { // THIS ONE IS FOR FLASHGG

	  if (motherPdgId_.empty()) { //or checkMother(p)
	    float dR = deltaR(p->p4(), ref->p4());
	    
	    if (dR < 0.2) // FIXME MAKE IT CONFIGURABLE
	      return std::pair<const reco::GenParticleRef, bool>(p, true);
	  }
	}
      }
    }
  }
  
  return std::pair<const reco::GenParticleRef, bool>(reco::GenParticleRef(), false);
}
  
//bool TagProbeFitTreeProducer::checkMother(const reco::GenParticleRef &ref) const {
//
//  if (ref.isNull()) {
//    if(ref->isPromptFinalState())
//      return true;
//    else
//      return false;
//  }
//
//
//  if (motherPdgId_.find(abs(ref->pdgId())) != motherPdgId_.end()) 
//    return true;
//
//  reco::GenParticle::mothers m = ref->motherRefVector();
//  for (reco::GenParticle::mothers::const_iterator it = m.begin(), e = m.end(); it != e; ++it) {
//    if (checkMother(*it))
//      return true;
//  }
//
//  return false;
//}

void TagProbeFitTreeProducer::beginLuminosityBlock(const edm::LuminosityBlock& iLumi, const edm::EventSetup& iEventS) 
{}
  
void TagProbeFitTreeProducer::endLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iEventS) {

  if (isMC_) {
    edm::Handle<edm::MergeableDouble> totWeight;
    iLumi.getByToken(totGenWeightToken_, totWeight);
    if (totWeight.isValid()) 
      totGenWeight_ += totWeight->value; 
    
    edm::Handle<edm::MergeableDouble> totEvents;
    iLumi.getByToken(totEventsToken_, totEvents);
    
    if (totEvents.isValid())
      totEvents_ += totEvents->value;
  }
}

void TagProbeFitTreeProducer::endJob() {
  // ask to write the current PSet info into the TTree header
  // FIXME !
  //treeFiller_->writeProvenance(edm::getProcessParameterSet());
  if (isMC_) {    
    treeFiller_->addTotWeightBranch(totGenWeight_, totEvents_);
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(TagProbeFitTreeProducer);
