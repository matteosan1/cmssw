import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import sys

process = cms.Process("tnp")

###################################################################
varOptions = VarParsing('analysis')
varOptions.register(
    "isMC",
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Compute MC efficiencies"
    )

varOptions.parseArguments()
if (varOptions.isMC == None):
    raise Exception("Must select either data or MC")

if (varOptions.isMC):
    from PhysicsTools.TagAndProbe.treeMakerOptionsMC_cfi import options
    process.pileupReweightingProducer = cms.EDProducer("PileupWeightProducer",
                                                   hardcodedWeights = cms.untracked.bool(True)
                                                   )

    process.GsfDRToNearestTauProbe = cms.EDProducer("DeltaRNearestGenPComputer",
                                                    probes = cms.InputTag(options['ELECTRON_COLL']),
                                                    objects = cms.InputTag('prunedGenParticles'),
                                                    objectSelection = cms.string("abs(pdgId)==15"),
                                                    )
    
    process.GsfDRToNearestTauSC = cms.EDProducer("DeltaRNearestGenPComputer",
                                                 probes = cms.InputTag("superClusterCands"),
                                                 objects = cms.InputTag('prunedGenParticles'),
                                                 objectSelection = cms.string("abs(pdgId)==15"),
                                                 )
    
    process.GsfDRToNearestTauTag = cms.EDProducer("DeltaRNearestGenPComputer",
                                                  probes = cms.InputTag(options['ELECTRON_COLL']),
                                                  objects = cms.InputTag('prunedGenParticles'),
                                                  objectSelection = cms.string("abs(pdgId)==15"),
                                                  )
else:
    from PhysicsTools.TagAndProbe.treeMakerOptionsData_cfi import options

###################################################################

process.eleVarHelper = cms.EDProducer("ElectronVariableHelper",
                                      probes = cms.InputTag(options['ELECTRON_COLL']),
                                      vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices")
                                      )

process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
process.hltHighLevel.throw = cms.bool(True)
process.hltHighLevel.HLTPaths = options['TnPPATHS']

process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = options['GLOBALTAG']

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(options['INPUT_FILE_NAME']),
                            eventsToProcess = options['EVENTSToPROCESS']
                            )

process.maxEvents = cms.untracked.PSet( input = options['MAXEVENTS'])
    
###################################################################
## ELECTRON MODULES
###################################################################

process.goodElectrons = cms.EDFilter("PATElectronRefSelector",
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
## ID MODULES
###################################################################

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

dataFormat = DataFormat.MiniAOD
if (options['useAOD']):
    dataFormat = DataFormat.AOD

switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff']
                 
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process, idmod, setupVIDElectronSelection)

# List of cuts
#0 MinPtCut
#1 GsfEleSCEtaMultiRangeCut
#2 GsfEleDEtaInCut
#3 GsfEleDPhiInCut
#4 GsfEleFull5x5SigmaIEtaIEtaCut
#5 GsfEleHadronicOverEMCut
#6 GsfEleDxyCut 
#7 GsfEleDzCut
#8 GsfEleEInverseMinusPInverseCut
#9 GsfEleEffAreaPFIsoCut
#10 GsfEleConversionVetoCut
#11 GsfEleMissingHitsCut

process.goodElectronsPROBECutBasedVeto = cms.EDProducer("PatElectronNm1Selector",
                                                        input     = cms.InputTag("goodElectrons"),
                                                        cut       = cms.string(options['ELECTRON_CUTS']),
                                                        selection = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
                                                        #cutIndicesToMask = cms.vuint32(2, 3),
                                                        cutNamesToMask = cms.vstring("GsfEleDEtaInCut_0", "GsfEleDPhiInCut_0")
                                                        )

process.goodElectronsPROBECutBasedLoose = process.goodElectronsPROBECutBasedVeto.clone()
process.goodElectronsPROBECutBasedLoose.selection = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose")
process.goodElectronsPROBECutBasedMedium = process.goodElectronsPROBECutBasedVeto.clone()
process.goodElectronsPROBECutBasedMedium.selection = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium")
process.goodElectronsPROBECutBasedTight = process.goodElectronsPROBECutBasedVeto.clone()
process.goodElectronsPROBECutBasedTight.selection = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight")


process.goodElectronsTAGCutBasedVeto = cms.EDProducer("PatElectronSelectorByValueMap",
                                                      input     = cms.InputTag("goodElectrons"),
                                                      cut       = cms.string(options['ELECTRON_TAG_CUTS']),
                                                      selection = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
                                                      id_cut    = cms.bool(True),
                                                      )

process.goodElectronsTAGCutBasedLoose = process.goodElectronsTAGCutBasedVeto.clone()
process.goodElectronsTAGCutBasedLoose.selection = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose")
process.goodElectronsTAGCutBasedMedium = process.goodElectronsTAGCutBasedVeto.clone()
process.goodElectronsTAGCutBasedMedium.selection = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium")
process.goodElectronsTAGCutBasedTight = process.goodElectronsTAGCutBasedVeto.clone()
process.goodElectronsTAGCutBasedTight.selection = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight")

###################################################################
## TRIGGER MATCHING
###################################################################

process.goodElectronsTagHLT = cms.EDProducer("PatElectronTriggerCandProducer",
                                             filterNames = cms.vstring(options['TnPHLTTagFilters']),
                                             inputs      = cms.InputTag("goodElectronsTAGCutBasedTight"),
                                             bits        = cms.InputTag('TriggerResults::HLT'),
                                             objects     = cms.InputTag('selectedPatTrigger'),
                                             dR          = cms.double(0.3),
                                             isAND       = cms.bool(True)
                                             )

process.goodElectronsProbeHLT = cms.EDProducer("PatElectronTriggerCandProducer",
                                               filterNames = cms.vstring(options['TnPHLTProbeFilters']),
                                               inputs      = cms.InputTag("goodElectrons"),
                                               bits        = cms.InputTag('TriggerResults::HLT'),
                                               objects     = cms.InputTag('selectedPatTrigger'),
                                               dR          = cms.double(0.3),
                                               isAND       = cms.bool(True)
                                               )

process.goodElectronsProbeMeasureHLT = cms.EDProducer("PatElectronTriggerCandProducer",
                                                      filterNames = cms.vstring(options['TnPHLTProbeFilters']),
                                                      inputs      = cms.InputTag("goodElectrons"),
                                                      bits        = cms.InputTag('TriggerResults::HLT'),
                                                      objects     = cms.InputTag('selectedPatTrigger'),
                                                      dR          = cms.double(0.3),
                                                      isAND       = cms.bool(True)
                                                      )

process.goodElectronsMeasureHLT = cms.EDProducer("PatElectronTriggerCandProducer",
                                                 filterNames = cms.vstring(options['HLTFILTERTOMEASURE']),
                                                 inputs      = cms.InputTag("goodElectronsProbeMeasureHLT"),
                                                 bits        = cms.InputTag('TriggerResults::HLT'),
                                                 objects     = cms.InputTag('selectedPatTrigger'),
                                                 dR          = cms.double(0.3),
                                                 isAND       = cms.bool(False)
                                                 )

process.goodSuperClustersHLT = cms.EDProducer("RecoEcalCandidateTriggerCandProducer",
                                              filterNames  = cms.vstring(options['TnPHLTProbeFilters']),
                                              inputs       = cms.InputTag("goodSuperClusters"),
                                              bits         = cms.InputTag('TriggerResults::HLT'),
                                              objects      = cms.InputTag('selectedPatTrigger'),
                                              dR           = cms.double(0.3),
                                              isAND        = cms.bool(True)
                                              )

###################################################################
## SEQUENCES
###################################################################

process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag(options['ELECTRON_COLL'])
process.ele_sequence = cms.Sequence(
    process.goodElectrons +
    process.egmGsfElectronIDSequence +
    process.goodElectronsPROBECutBasedVeto +
    process.goodElectronsPROBECutBasedLoose +
    process.goodElectronsPROBECutBasedMedium +
    process.goodElectronsPROBECutBasedTight +
    process.goodElectronsTAGCutBasedVeto +
    process.goodElectronsTAGCutBasedLoose +
    process.goodElectronsTAGCutBasedMedium +
    process.goodElectronsTAGCutBasedTight +
    process.goodElectronsTagHLT +
    process.goodElectronsProbeHLT +
    process.goodElectronsProbeMeasureHLT +
    process.goodElectronsMeasureHLT
    )

process.sc_sequence = cms.Sequence(process.superClusterCands +
                                   process.goodSuperClusters +
                                   process.goodSuperClustersHLT +
                                   process.GsfMatchedSuperClusterCands
                                   )

###################################################################
## TnP PAIRS
###################################################################

process.tagTightHLT = cms.EDProducer("CandViewShallowCloneCombiner",
                                     decay = cms.string("goodElectronsTagHLT@+ goodElectronsProbeMeasureHLT@-"), 
                                     checkCharge = cms.bool(True),
                                     cut = cms.string("40<mass<1000"),
                                     )

process.tagTightSC = cms.EDProducer("CandViewShallowCloneCombiner",
                                    decay = cms.string("goodElectronsTagHLT goodSuperClustersHLT"), 
                                    checkCharge = cms.bool(False),
                                    cut = cms.string("40<mass<1000"),
                                    )

process.tagTightRECO = cms.EDProducer("CandViewShallowCloneCombiner",
                                      decay = cms.string("goodElectronsTagHLT@+ goodElectronsProbeHLT@-"), 
                                      checkCharge = cms.bool(True),
                                      cut = cms.string("40<mass<1000"),
                                    )


process.allTagsAndProbes = cms.Sequence()

if (options['DOTRIGGER']):
    process.allTagsAndProbes *= process.tagTightHLT

if (options['DORECO']):
    process.allTagsAndProbes *= process.tagTightSC

if (options['DOID']):
    process.allTagsAndProbes *= process.tagTightRECO

###################################################################
## MC MATCHES
###################################################################

process.McMatchHLT = cms.EDProducer("MCTruthDeltaRMatcherNew",
                                    matchPDGId = cms.vint32(11),
                                    src = cms.InputTag("goodElectrons"),
                                    distMin = cms.double(0.3),
                                    matched = cms.InputTag("prunedGenParticles"),
                                    checkCharge = cms.bool(True)
                                    )

process.McMatchSC = cms.EDProducer("MCTruthDeltaRMatcherNew",
                                   matchPDGId = cms.vint32(11),
                                   src = cms.InputTag("goodSuperClusters"),
                                   distMin = cms.double(0.3),
                                   matched = cms.InputTag("prunedGenParticles"),
                                   checkCharge = cms.bool(False)
                                   )
                     
process.McMatchTag = cms.EDProducer("MCTruthDeltaRMatcherNew",
                                    matchPDGId = cms.vint32(11),
                                    src = cms.InputTag("goodElectronsTAGCutBasedTight"),
                                    distMin = cms.double(0.2),
                                    matched = cms.InputTag("prunedGenParticles"),
                                    checkCharge = cms.bool(True)
                                    )

process.McMatchRECO = cms.EDProducer("MCTruthDeltaRMatcherNew",
                                     matchPDGId = cms.vint32(11),
                                     src = cms.InputTag("goodElectrons"),
                                     distMin = cms.double(0.2),
                                     matched = cms.InputTag("prunedGenParticles"),
                                     checkCharge = cms.bool(True)
                                    )

process.mc_sequence = cms.Sequence()

if (varOptions.isMC):
    process.mc_sequence *= (process.McMatchHLT + process.McMatchTag + process.McMatchSC + process.McMatchRECO)

##########################################################################
## TREE CONTENT
#########################################################################

ZVariablesToStore = cms.PSet(
    eta = cms.string("eta"),
    abseta = cms.string("abs(eta)"),
    pt  = cms.string("pt"),
    mass  = cms.string("mass"),
    )   

SCProbeVariablesToStore = cms.PSet(
    probe_sc_eta    = cms.string("eta"),
    probe_sc_abseta = cms.string("abs(eta)"),
    probe_sc_pt     = cms.string("pt"),
    probe_sc_et     = cms.string("et"),
    probe_sc_e      = cms.string("energy"),
    
    probe_dRTau  = cms.InputTag("GsfDRToNearestTauSC")
)

ProbeVariablesToStore = cms.PSet(
    probe_Ele_eta    = cms.string("eta"),
    probe_Ele_abseta = cms.string("abs(eta)"),
    probe_Ele_pt     = cms.string("pt"),
    probe_Ele_et     = cms.string("et"),
    probe_Ele_e      = cms.string("energy"),
    probe_Ele_q      = cms.string("charge"),

## super cluster quantities
    probe_sc_energy = cms.string("superCluster.energy"),
    probe_sc_et     = cms.string("superCluster.energy*sin(superClusterPosition.theta)"),    
    probe_sc_eta    = cms.string("superCluster.eta"),
    probe_sc_abseta = cms.string("abs(superCluster.eta)"),

#id based
    probe_Ele_dEtaIn        = cms.string("deltaEtaSuperClusterTrackAtVtx"),
    probe_Ele_dPhiIn        = cms.string("deltaPhiSuperClusterTrackAtVtx"),
    probe_Ele_sigmaIEtaIEta = cms.string("sigmaIetaIeta"),
    probe_Ele_hoe           = cms.string("hadronicOverEm"),
    probe_Ele_ooemoop       = cms.string("(1.0/ecalEnergy - eSuperClusterOverP/ecalEnergy)"),
    probe_Ele_mHits         = cms.InputTag("eleVarHelper:missinghits"),
    probe_Ele_dz            = cms.InputTag("eleVarHelper:dz"),
    probe_Ele_dxy           = cms.InputTag("eleVarHelper:dxy"),
    probe_Ele_mva           = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Phys14NonTrigValues"),
    probe_dRTau             = cms.InputTag("GsfDRToNearestTauProbe"),

#isolation
    probe_Ele_chIso         = cms.string("pfIsolationVariables().sumChargedHadronPt"),
    probe_Ele_phoIso        = cms.string("pfIsolationVariables().sumPhotonEt"),
    probe_Ele_neuIso        = cms.string("pfIsolationVariables().sumNeutralHadronEt"),
)

TagVariablesToStore = cms.PSet(
    Ele_eta    = cms.string("eta"),
    Ele_abseta = cms.string("abs(eta)"),
    Ele_pt     = cms.string("pt"),
    Ele_et     = cms.string("et"),
    Ele_e      = cms.string("energy"),
    Ele_q      = cms.string("charge"),
    Ele_dRTau = cms.InputTag("GsfDRToNearestTauTag"),
    
    ## super cluster quantities
    sc_energy = cms.string("superCluster.energy"),
    sc_et     = cms.string("superCluster.energy*sin(superClusterPosition.theta)"),    
    sc_eta    = cms.string("superCluster.eta"),
    sc_abseta = cms.string("abs(superCluster.eta)"),
)

CommonStuffForGsfElectronProbe = cms.PSet(
    variables = cms.PSet(ProbeVariablesToStore),
    ignoreExceptions =  cms.bool (True),
    addRunLumiInfo   =  cms.bool (True),
    addEventVariablesInfo   =  cms.bool(True),
    vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    #pfMet = cms.InputTag(""),
    pairVariables =  cms.PSet(ZVariablesToStore),
    pairFlags     =  cms.PSet(
        mass60to120 = cms.string("60 < mass < 120")
        ),
    tagVariables   =  cms.PSet(TagVariablesToStore),
    tagFlags       =  cms.PSet(),    
    )

CommonStuffForSuperClusterProbe = CommonStuffForGsfElectronProbe.clone()
CommonStuffForSuperClusterProbe.variables = cms.PSet(SCProbeVariablesToStore)

if varOptions.isMC:
    mcTruthCommonStuff = cms.PSet(
        isMC = cms.bool(True),
        tagMatches = cms.InputTag("McMatchTag"),
        motherPdgId = cms.vint32(22,23),
        #motherPdgId = cms.vint32(443), # JPsi
        #motherPdgId = cms.vint32(553), # Yupsilon
        makeMCUnbiasTree = cms.bool(False),
        checkMotherInUnbiasEff = cms.bool(False),
        mcVariables = cms.PSet(
            probe_eta = cms.string("eta"),
            probe_abseta = cms.string("abs(eta)"),
            probe_et  = cms.string("et"),
            probe_e  = cms.string("energy"),
            ),
        mcFlags     =  cms.PSet(
            probe_flag = cms.string("pt>0")
            ),      
        )
else:
    mcTruthCommonStuff = cms.PSet(
        isMC = cms.bool(False)
        )

##########################################################################
## TREE MAKER OPTIONS
##########################################################################

process.GsfElectronToTrigger = cms.EDAnalyzer("TagProbeFitTreeProducer",
                                              CommonStuffForSuperClusterProbe, mcTruthCommonStuff,
                                              tagProbePairs = cms.InputTag("tagTightHLT"),
                                              arbitration   = cms.string("None"),
                                              flags         = cms.PSet(passingHLT    = cms.InputTag("goodElectronsMeasureHLT")
                                                                       ),                                               
                                              allProbes     = cms.InputTag("goodElectronsProbeMeasureHLT"),
                                              )

if (varOptions.isMC):
    process.GsfElectronToTrigger.probeMatches  = cms.InputTag("McMatchHLT")
    process.GsfElectronToTrigger.eventWeight   = cms.InputTag("generator")
    process.GsfElectronToTrigger.PUWeightSrc   = cms.InputTag("pileupReweightingProducer","pileupWeights")

process.GsfElectronToSC = cms.EDAnalyzer("TagProbeFitTreeProducer",
                                         CommonStuffForSuperClusterProbe, mcTruthCommonStuff,
                                         tagProbePairs = cms.InputTag("tagTightSC"),
                                         arbitration   = cms.string("None"),
                                         flags         = cms.PSet(passingRECO   = cms.InputTag("GsfMatchedSuperClusterCands", "superclusters"),         
                                                                  ),                                               
                                         allProbes     = cms.InputTag("goodSuperClustersHLT"),
                                         )

if (varOptions.isMC):
    process.GsfElectronToSC.probeMatches  = cms.InputTag("McMatchSC")
    process.GsfElectronToSC.eventWeight   = cms.InputTag("generator")
    process.GsfElectronToSC.PUWeightSrc   = cms.InputTag("pileupReweightingProducer","pileupWeights")

process.GsfElectronToRECO = cms.EDAnalyzer("TagProbeFitTreeProducer",
                                           mcTruthCommonStuff, CommonStuffForGsfElectronProbe,
                                           tagProbePairs = cms.InputTag("tagTightRECO"),
                                           arbitration   = cms.string("None"),
                                           flags         = cms.PSet(passingVeto   = cms.InputTag("goodElectronsPROBECutBasedVeto"),
                                                                    passingLoose  = cms.InputTag("goodElectronsPROBECutBasedLoose"),
                                                                    passingMedium = cms.InputTag("goodElectronsPROBECutBasedMedium"),
                                                                    passingTight  = cms.InputTag("goodElectronsPROBECutBasedTight"),
                                                                    ),                                               
                                           allProbes     = cms.InputTag("goodElectronsProbeHLT"),
                                           )

if (varOptions.isMC):
    process.GsfElectronToRECO.probeMatches  = cms.InputTag("McMatchRECO")
    process.GsfElectronToRECO.eventWeight   = cms.InputTag("generator")
    process.GsfElectronToRECO.PUWeightSrc   = cms.InputTag("pileupReweightingProducer","pileupWeights")

process.tree_sequence = cms.Sequence()
if (options['DOTRIGGER']):
    process.tree_sequence *= process.GsfElectronToTrigger

if (options['DORECO']):
    process.tree_sequence *= process.GsfElectronToSC

if (options['DOID']):
    process.tree_sequence *= process.GsfElectronToRECO

##########################################################################
##    ____       _   _     
##   |  _ \ __ _| |_| |__  
##   | |_) / _` | __| '_ \ 
##   |  __/ (_| | |_| | | |
##   |_|   \__,_|\__|_| |_|
##
##########################################################################

process.out = cms.OutputModule("PoolOutputModule", 
                               fileName = cms.untracked.string(options['OUTPUTEDMFILENAME']),
                               SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("p"))
                               )
process.outpath = cms.EndPath(process.out)
if (not options['DEBUG']):
    process.outpath.remove(process.out)

if (varOptions.isMC):
    process.p = cms.Path(
        process.hltHighLevel +
        process.ele_sequence + 
        process.sc_sequence +
        process.allTagsAndProbes +
        process.pileupReweightingProducer +
        process.mc_sequence +
        process.eleVarHelper +
        process.GsfDRToNearestTauProbe + 
        process.GsfDRToNearestTauTag + 
        process.GsfDRToNearestTauSC + 
        process.tree_sequence
        )
else:
    process.p = cms.Path(
        process.hltHighLevel +
        process.ele_sequence + 
        process.sc_sequence +
        process.allTagsAndProbes +
        process.mc_sequence +
        process.eleVarHelper +
        process.tree_sequence
        )

process.TFileService = cms.Service(
    "TFileService", fileName = cms.string(options['OUTPUT_FILE_NAME'])
    )
