import FWCore.ParameterSet.Config as cms

def setModules(process, options):

    process.sampleInfo = cms.EDProducer("tnp::SampleInfoTree",
                                        #vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                        genInfo = cms.InputTag("generator")
                                        )

    process.eleVarHelper = cms.EDProducer("GsfElectronVariableHelper",
                                          probes = cms.InputTag(options['ELECTRON_COLL']),
                                          vertexCollection = cms.InputTag("offlinePrimaryVertices")
                                          )

    process.hltVarHelper = cms.EDProducer("GsfElectronHLTVariableHelper",
                                          probes = cms.InputTag(options['ELECTRON_COLL']),
                                          hltCandidateCollection = cms.InputTag("hltEgammaCandidates"),
                                          mapOutputNames = cms.vstring("hltsieie",
                                                                       "hltecaliso",
                                                                       "hlthcaliso",
                                                                       "hlthoe",
                                                                       "hlttkiso",
                                                                       "hltdeta",
                                                                       "hltdetaseed",
                                                                       "hltdphi",
                                                                       "hlteop",
                                                                       "hltmishits"),
                                          mapInputTags = cms.VInputTag("hltEgammaClusterShape:sigmaIEtaIEta5x5",
                                                                       "hltEgammaEcalPFClusterIso",
                                                                       "hltEgammaHcalPFClusterIso",
                                                                       "hltEgammaHoverE", 
                                                                       "hltEgammaEleGsfTrackIso",
                                                                       "hltEgammaGsfTrackVars:Deta",
                                                                       "hltEgammaGsfTrackVars:DetaSeed",
                                                                       "hltEgammaGsfTrackVars:Dphi",
                                                                       "hltEgammaGsfTrackVars:OneOESuperMinusOneOP",
                                                                       "hltEgammaGsfTrackVars:MissingHits")
                                          )

    from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
    process.hltFilter = hltHighLevel.clone()
    process.hltFilter.throw = cms.bool(True)
    process.hltFilter.HLTPaths = options['TnPPATHS']

    from PhysicsTools.TagAndProbe.pileupConfiguration_cfi import pileupProducer
    process.pileupReweightingProducer = pileupProducer.clone()

###################################################################                                                                               
## ELECTRON MODULES                                                                                                                                    
###################################################################                                    
    
    process.goodElectrons = cms.EDFilter("GsfElectronRefSelector",
                                         src = cms.InputTag(options['ELECTRON_COLL']),
                                         cut = cms.string(options['ELECTRON_CUTS'])
                                         )
    
###################################################################                                                                     
## SUPERCLUSTER MODULES                                                     
###################################################################         
    
    process.superClusterCands = cms.EDProducer("ConcreteEcalCandidateProducer",
                                               src = cms.InputTag(options['SUPERCLUSTER_COLL']),
                                               particleType = cms.int32(11),
                                               )
    
    process.goodSuperClusters = cms.EDFilter("RecoEcalCandidateRefSelector",
                                             src = cms.InputTag("superClusterCands"),
                                             cut = cms.string(options['SUPERCLUSTER_CUTS']),
                                             filter = cms.bool(True)
                                             )
    
    process.GsfMatchedSuperClusterCands = cms.EDProducer("ElectronMatchedCandidateProducer",
                                                         src     = cms.InputTag("superClusterCands"),
                                                         ReferenceElectronCollection = cms.untracked.InputTag("goodElectrons"),
                                                         cut = cms.string(options['SUPERCLUSTER_CUTS'])
                                                         )
    
###################################################################
## TRIGGER MATCHING
###################################################################
    
    process.goodElectronsTagHLT = cms.EDProducer("GsfElectronTriggerCandProducer",
                                                 filterNames = cms.vstring(options['TnPHLTTagFilters']),
                                                 inputs      = cms.InputTag("goodElectronsTAGCutBasedTight"),
                                                 bits        = cms.InputTag('TriggerResults::HLT'),
                                                 objects     = cms.InputTag('hltTriggerSummaryAOD'),
                                                 dR          = cms.double(0.3),
                                                 isAND       = cms.bool(True)
                                                 )
    
    process.goodElectronsProbeHLT = cms.EDProducer("GsfElectronTriggerCandProducer",
                                                   filterNames = cms.vstring(options['TnPHLTProbeFilters']),
                                                   inputs      = cms.InputTag("goodElectrons"),
                                                   bits        = cms.InputTag('TriggerResults::HLT'),
                                                   objects     = cms.InputTag('hltTriggerSummaryAOD'),
                                                   dR          = cms.double(0.3),
                                                   isAND       = cms.bool(True)
                                                   )
    
    process.goodElectronsProbeMeasureHLT = cms.EDProducer("GsfElectronTriggerCandProducer",
                                                          filterNames = cms.vstring(options['TnPHLTProbeFilters']),
                                                          inputs      = cms.InputTag("goodElectrons"),
                                                          bits        = cms.InputTag('TriggerResults::HLT'),
                                                          objects     = cms.InputTag('hltTriggerSummaryAOD'),
                                                          dR          = cms.double(0.3),
                                                          isAND       = cms.bool(True)
                                                          )
    
    process.goodElectronsMeasureHLT = cms.EDProducer("GsfElectronTriggerCandProducer",
                                                     filterNames = cms.vstring(options['HLTFILTERTOMEASURE']),
                                                     inputs      = cms.InputTag("goodElectronsProbeMeasureHLT"),
                                                     bits        = cms.InputTag('TriggerResults::HLT'),
                                                     objects     = cms.InputTag('hltTriggerSummaryAOD'),
                                                     dR          = cms.double(0.3),
                                                     isAND       = cms.bool(False)
                                                     )
    
    process.goodSuperClustersHLT = cms.EDProducer("RecoEcalCandidateTriggerCandProducer",
                                                  filterNames  = cms.vstring(options['TnPHLTProbeFilters']),
                                                  inputs       = cms.InputTag("goodSuperClusters"),
                                                  bits         = cms.InputTag('TriggerResults::HLT'),
                                                  objects      = cms.InputTag('hltTriggerSummaryAOD'),
                                                  dR           = cms.double(0.3),
                                                  isAND        = cms.bool(True)
                                                  )
    
###################################################################
## TnP PAIRS
###################################################################
    
    process.tagTightHLT = cms.EDProducer("CandViewShallowCloneCombiner",
                                         decay = cms.string("goodElectronsTagHLT@+ goodElectronsProbeMeasureHLT@-"), 
                                         checkCharge = cms.bool(True),
                                         cut = cms.string("60<mass<120"),
                                         )
    
    process.tagTightSC = cms.EDProducer("CandViewShallowCloneCombiner",
                                        decay = cms.string("goodElectronsTagHLT goodSuperClustersHLT"), 
                                        checkCharge = cms.bool(False),
                                        cut = cms.string("60<mass<120"),
                                        )
    
    process.tagTightRECO = cms.EDProducer("CandViewShallowCloneCombiner",
                                          decay = cms.string("goodElectronsTagHLT@+ goodElectronsProbeHLT@-"), 
                                          checkCharge = cms.bool(True),
                                          cut = cms.string("60<mass<120"),
                                          )
    
