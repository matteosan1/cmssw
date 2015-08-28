import FWCore.ParameterSet.Config as cms
import sys

from PhysicsTools.TagAndProbe.treeMakerOptionsMC_cfi import options

process = cms.Process("tnp")

process.pileupReweightingProducer = cms.EDProducer("PileupWeightProducer",
                                                   hardcodedWeights = cms.untracked.bool(True)
                                                   )

process.eleVarHelper = cms.EDProducer("ElectronVariableHelper",
                                      probes = cms.InputTag(options['ELECTRON_COLL']),
                                      vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices")
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

process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
process.hltHighLevel.throw = cms.bool(True)
process.hltHighLevel.HLTPaths = options['TnPPATHS']

###################################################################
##    ___            _           _      
##   |_ _|_ __   ___| |_   _  __| | ___ 
##    | || '_ \ / __| | | | |/ _` |/ _ \
##    | || | | | (__| | |_| | (_| |  __/
##   |___|_| |_|\___|_|\__,_|\__,_|\___|
##
###################################################################

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

###################################################################
##   ____             _ ____ 
##  |  _ \ ___   ___ | / ___|  ___  _   _ _ __ ___ ___ 
##  | |_) / _ \ / _ \| \___ \ / _ \| | | | '__/ __/ _ \
##  |  __/ (_) | (_) | |___) | (_) | |_| | | | (_|  __/
##  |_|   \___/ \___/|_|____/ \___/ \__,_|_|  \___\___|
##  
###################################################################

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(options['INPUT_FILE_NAME']),
                            eventsToProcess = options['EVENTSToPROCESS']
                            )

process.maxEvents = cms.untracked.PSet( input = options['MAXEVENTS'])
    
###################################################################
##    ____      __ _____ _           _                   
##   / ___|___ / _| ____| | ___  ___| |_ _ __ ___  _ __  
##  | |  _/ __| |_|  _| | |/ _ \/ __| __| '__/ _ \| '_ \ 
##  | |_| \__ \  _| |___| |  __/ (__| |_| | | (_) | | | |
##   \____|___/_| |_____|_|\___|\___|\__|_|  \___/|_| |_|
##  
###################################################################

process.goodElectrons = cms.EDFilter("PATElectronRefSelector",
                                     src = cms.InputTag(options['ELECTRON_COLL']),
                                     cut = cms.string(options['ELECTRON_CUTS'])    
                                     )

###################################################################
##   ____                         ____ _           _            
##  / ___| _   _ _ __   ___ _ __ / ___| |_   _ ___| |_ ___ _ __ 
##  \___ \| | | | '_ \ / _ \ '__| |   | | | | / __| __/ _ \ '__|
##   ___) | |_| | |_) |  __/ |  | |___| | |_| \__ \ ||  __/ |   
##  |____/ \__,_| .__/ \___|_|   \____|_|\__,_|___/\__\___|_|   
##  
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
##    _____ _           _                     ___    _ 
##   | ____| | ___  ___| |_ _ __ ___  _ __   |_ _|__| |
##   |  _| | |/ _ \/ __| __| '__/ _ \| '_ \   | |/ _` |
##   | |___| |  __/ (__| |_| | | (_) | | | |  | | (_| |
##   |_____|_|\___|\___|\__|_|  \___/|_| |_| |___\__,_|
##   
###################################################################

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

dataFormat = DataFormat.MiniAOD
if (options['useAOD']):
    dataFormat = DataFormat.AOD

switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V2_cff',
                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV51_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_PHYS14_PU20bx25_nonTrig_V1_cff']
                 

for idmod in my_id_modules:
    setupAllVIDIdsInModule(process, idmod, setupVIDElectronSelection)

process.goodElectronsPROBECutBasedVeto = cms.EDProducer("PatElectronSelectorByValueMap",
                                                        input     = cms.InputTag("goodElectrons"), #options['ELECTRON_COLL']),
                                                        cut       = cms.string(options['ELECTRON_CUTS']),
                                                        selection = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto"),
                                                        id_cut    = cms.bool(True)
                                                        )

process.goodElectronsPROBECutBasedLoose = process.goodElectronsPROBECutBasedVeto.clone()
process.goodElectronsPROBECutBasedLoose.selection = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-loose")
process.goodElectronsPROBECutBasedMedium = process.goodElectronsPROBECutBasedVeto.clone()
process.goodElectronsPROBECutBasedMedium.selection = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium")
process.goodElectronsPROBECutBasedTight = process.goodElectronsPROBECutBasedVeto.clone()
process.goodElectronsPROBECutBasedTight.selection = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-tight")

process.goodElectronsTAGCutBasedVeto = cms.EDProducer("PatElectronSelectorByValueMap",
                                                   input     = cms.InputTag("goodElectrons"), #options['ELECTRON_COLL']),
                                                   cut       = cms.string(options['ELECTRON_TAG_CUTS']),
                                                   selection = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto"),
                                                   id_cut    = cms.bool(True)
                                                   )

process.goodElectronsTAGCutBasedLoose = process.goodElectronsTAGCutBasedVeto.clone()
process.goodElectronsTAGCutBasedLoose.selection = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-loose")
process.goodElectronsTAGCutBasedMedium = process.goodElectronsTAGCutBasedVeto.clone()
process.goodElectronsTAGCutBasedMedium.selection = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium")
process.goodElectronsTAGCutBasedTight = process.goodElectronsTAGCutBasedVeto.clone()
process.goodElectronsTAGCutBasedTight.selection = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-tight")


###################################################################
##    _____     _                         __  __       _       _     _             
##   |_   _| __(_) __ _  __ _  ___ _ __  |  \/  | __ _| |_ ___| |__ (_)_ __   __ _ 
##     | || '__| |/ _` |/ _` |/ _ \ '__| | |\/| |/ _` | __/ __| '_ \| | '_ \ / _` |
##     | || |  | | (_| | (_| |  __/ |    | |  | | (_| | || (__| | | | | | | | (_| |
##     |_||_|  |_|\__, |\__, |\___|_|    |_|  |_|\__,_|\__\___|_| |_|_|_| |_|\__, |
##                |___/ |___/                                                |___/ 
###################################################################

#if (len(options['TnPHLTTagFilters']) != len(options['TnPHLTProbeFilters'])):
#    print "ERROR: different number of tag and probe filters, please fix it"
#    sys.exit(-1)

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
##    _____ ___   ____    ____       _          
##   |_   _( _ ) |  _ \  |  _ \ __ _(_)_ __ ___ 
##     | | / _ \/\ |_) | | |_) / _` | | '__/ __|
##     | || (_>  <  __/  |  __/ (_| | | |  \__ \
##     |_| \___/\/_|     |_|   \__,_|_|_|  |___/
##   
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
##    __  __  ____   __  __       _       _               
##   |  \/  |/ ___| |  \/  | __ _| |_ ___| |__   ___  ___ 
##   | |\/| | |     | |\/| |/ _` | __/ __| '_ \ / _ \/ __|
##   | |  | | |___  | |  | | (_| | || (__| | | |  __/\__ \
##   |_|  |_|\____| |_|  |_|\__,_|\__\___|_| |_|\___||___/
##                                   
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

if (options['MC_FLAG']):
    process.mc_sequence *= (process.McMatchHLT + process.McMatchTag + process.McMatchSC + process.McMatchRECO)

##########################################################################
##    _____           _       _ ____            _            _   _  ____  
##   |_   _|_ _  __ _( )_ __ ( )  _ \ _ __ ___ | |__   ___  | \ | |/ ___| 
##     | |/ _` |/ _` |/| '_ \|/| |_) | '__/ _ \| '_ \ / _ \ |  \| | |  _  
##     | | (_| | (_| | | | | | |  __/| | | (_) | |_) |  __/ | |\  | |_| | 
##     |_|\__,_|\__, | |_| |_| |_|   |_|  \___/|_.__/ \___| |_| \_|\____| 
##              |___/                                                     
##                                                                        
##########################################################################
##    ____                      _     _           
##   |  _ \ ___ _   _ ___  __ _| |__ | | ___  ___ 
##   | |_) / _ \ | | / __|/ _` | '_ \| |/ _ \/ __|
##   |  _ <  __/ |_| \__ \ (_| | |_) | |  __/\__ \
##   |_| \_\___|\__,_|___/\__,_|_.__/|_|\___||___/
##
## I define some common variables for re-use later.
## This will save us repeating the same code for each efficiency category
#########################################################################

ZVariablesToStore = cms.PSet(
    eta = cms.string("eta"),
    abseta = cms.string("abs(eta)"),
    pt  = cms.string("pt"),
    mass  = cms.string("mass"),
    )   

SCProbeVariablesToStore = cms.PSet(
    probe_eta    = cms.string("eta"),
    probe_abseta = cms.string("abs(eta)"),
    probe_pt     = cms.string("pt"),
    probe_et     = cms.string("et"),
    probe_e      = cms.string("energy"),
    
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
    pairVariables =  cms.PSet(ZVariablesToStore),
    pairFlags     =  cms.PSet(
        mass60to120 = cms.string("60 < mass < 120")
        ),
    tagVariables   =  cms.PSet(TagVariablesToStore),
    tagFlags       =  cms.PSet(),    
    )

CommonStuffForSuperClusterProbe = CommonStuffForGsfElectronProbe.clone()
CommonStuffForSuperClusterProbe.variables = cms.PSet(SCProbeVariablesToStore)

if options['MC_FLAG']:
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
##   ____      __       __    ___                 ___    _ 
##  / ___|___ / _|      \ \  |_ _|___  ___       |_ _|__| |
## | |  _/ __| |_   _____\ \  | |/ __|/ _ \       | |/ _` |
## | |_| \__ \  _| |_____/ /  | |\__ \ (_) |  _   | | (_| |
##  \____|___/_|        /_/  |___|___/\___/  ( ) |___\__,_|
##                                           |/            
##########################################################################

process.GsfElectronToTrigger = cms.EDAnalyzer("TagProbeFitTreeProducer",
                                              CommonStuffForSuperClusterProbe, mcTruthCommonStuff,
                                              tagProbePairs = cms.InputTag("tagTightHLT"),
                                              arbitration   = cms.string("None"),
                                              flags         = cms.PSet(passingHLT    = cms.InputTag("goodElectronsMeasureHLT")
                                                                       ),                                               
                                              allProbes     = cms.InputTag("goodElectronsProbeMeasureHLT"),
                                              )

if (options['MC_FLAG']):
    process.GsfElectronToTrigger.probeMatches  = cms.InputTag("McMatchHLT")
    process.GsfElectronToTrigger.eventWeight   = cms.InputTag("generator")
    process.GsfElectronToTrigger.PUWeightSrc   = cms.InputTag("pileupReweightingProducer","pileupWeights")

process.GsfElectronToSC = cms.EDAnalyzer("TagProbeFitTreeProducer",
                                         CommonStuffForSuperClusterProbe, mcTruthCommonStuff,
                                         tagProbePairs = cms.InputTag("tagTightSC"),
                                         arbitration   = cms.string("None"),
                                         flags         = cms.PSet(passingRECO   = cms.InputTag("GsfMatchedSuperClusterCands"),                                                                  
                                                                  ),                                               
                                         allProbes     = cms.InputTag("goodSuperClustersHLT"),
                                         )

if (options['MC_FLAG']):
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

if (options['MC_FLAG']):
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

if (options['MC_FLAG']):
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
