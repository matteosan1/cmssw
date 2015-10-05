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
    from PhysicsTools.TagAndProbe.treeMakerOptionsPhotonsMC_cfi import options
    process.pileupReweightingProducer = cms.EDProducer("PileupWeightProducer",
                                                       hardcodedWeights = cms.untracked.bool(True)
                                                       )

    process.GsfDRToNearestTauProbe = cms.EDProducer("DeltaRNearestGenPComputer",
                                                    probes = cms.InputTag(options['PHOTON_COLL']),
                                                    objects = cms.InputTag('prunedGenParticles'),
                                                    objectSelection = cms.string("abs(pdgId)==15"),
                                                    )

    process.GsfDRToNearestTauTag = cms.EDProducer("DeltaRNearestGenPComputer",
                                                  probes = cms.InputTag(options['PHOTON_COLL']),
                                                  objects = cms.InputTag('prunedGenParticles'),
                                                  objectSelection = cms.string("abs(pdgId)==15"),
                                                  )
else:
    from PhysicsTools.TagAndProbe.treeMakerOptionsPhotonsData_cfi import options

###################################################################

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

process.goodPhotons = cms.EDFilter("PATPhotonRefSelector",
                                    src = cms.InputTag(options['PHOTON_COLL']),
                                    cut = cms.string(options['PHOTON_CUTS'])    
                                    )

###################################################################
## IDs
###################################################################

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

dataFormat = DataFormat.MiniAOD
if (options['useAOD']):
    dataFormat = DataFormat.AOD

switchOnVIDPhotonIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring15_50ns_V1_cff',
                 'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring15_50ns_nonTrig_V2_cff']

for idmod in my_id_modules:
    setupAllVIDIdsInModule(process, idmod, setupVIDPhotonSelection)

process.goodPhotonsPROBECutBasedLoose = cms.EDProducer("PatPhotonSelectorByValueMap",
                                                       input     = cms.InputTag("goodPhotons"),
                                                       cut       = cms.string(options['PHOTON_CUTS']),
                                                       selection = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-50ns-V1-standalone-loose"),
                                                       id_cut    = cms.bool(True)
                                                       )

process.goodPhotonsPROBECutBasedMedium = process.goodPhotonsPROBECutBasedLoose.clone()
process.goodPhotonsPROBECutBasedMedium.selection = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-50ns-V1-standalone-medium")
process.goodPhotonsPROBECutBasedTight = process.goodPhotonsPROBECutBasedLoose.clone()
process.goodPhotonsPROBECutBasedTight.selection = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-50ns-V1-standalone-tight")

process.goodPhotonsTAGCutBasedLoose = cms.EDProducer("PatPhotonSelectorByValueMap",
                                                     input     = cms.InputTag("goodPhotons"),
                                                     cut       = cms.string(options['PHOTON_TAG_CUTS']),
                                                     selection = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-50ns-V1-standalone-loose"),
                                                     id_cut    = cms.bool(True)
                                                     )

process.goodPhotonsTAGCutBasedMedium = process.goodPhotonsTAGCutBasedLoose.clone()
process.goodPhotonsTAGCutBasedMedium.selection = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-50ns-V1-standalone-medium")
process.goodPhotonsTAGCutBasedTight = process.goodPhotonsTAGCutBasedLoose.clone()
process.goodPhotonsTAGCutBasedTight.selection = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-50ns-V1-standalone-tight")


###################################################################
## PHOTON ISOLATION
###################################################################
process.load("RecoEgamma/PhotonIdentification/PhotonIDValueMapProducer_cfi")

###################################################################
## HLT MATCHING
###################################################################

process.goodPhotonsTagHLT = cms.EDProducer("PatPhotonTriggerCandProducer",
                                           filterNames = options['TnPHLTTagFilters'],
                                           inputs      = cms.InputTag("goodPhotonsTAGCutBasedTight"),
                                           bits        = cms.InputTag('TriggerResults::HLT'),
                                           objects     = cms.InputTag('selectedPatTrigger'),
                                           dR          = cms.double(0.3),
                                           isAND       = cms.bool(True)
                                           )

process.goodPhotonsProbeHLT = cms.EDProducer("PatPhotonTriggerCandProducer",
                                             filterNames = options['TnPHLTProbeFilters'],
                                             inputs      = cms.InputTag("goodPhotons"),
                                             bits        = cms.InputTag('TriggerResults::HLT'),
                                             objects     = cms.InputTag('selectedPatTrigger'),
                                             dR          = cms.double(0.3),
                                             isAND       = cms.bool(True)
                                             )

process.egmPhotonIDs.physicsObjectSrc = cms.InputTag(options['PHOTON_COLL'])
process.pho_sequence = cms.Sequence(
    process.goodPhotons +
    process.egmPhotonIDSequence +
    process.photonIDValueMapProducer +
    process.goodPhotonsPROBECutBasedLoose +
    process.goodPhotonsPROBECutBasedMedium +
    process.goodPhotonsPROBECutBasedTight +
    process.goodPhotonsTAGCutBasedLoose +
    process.goodPhotonsTAGCutBasedMedium +
    process.goodPhotonsTAGCutBasedTight +
    process.goodPhotonsTagHLT +
    process.goodPhotonsProbeHLT 
    )

###################################################################
## TnP PAIRS
###################################################################

process.tagTightRECO = cms.EDProducer("CandViewShallowCloneCombiner",
                                      decay = cms.string("goodPhotonsTagHLT@+ goodPhotonsProbeHLT@-"), 
                                      checkCharge = cms.bool(False),
                                      cut = cms.string("40<mass<1000"),
                                      )

process.allTagsAndProbes = cms.Sequence(process.tagTightRECO)

###################################################################
## MC MATCHING
###################################################################
                     
process.McMatchTag = cms.EDProducer("MCTruthDeltaRMatcherNew",
                                    matchPDGId = cms.vint32(11),
                                    src = cms.InputTag("goodPhotonsTAGCutBasedTight"),
                                    distMin = cms.double(0.2),
                                    matched = cms.InputTag("prunedGenParticles"),
                                    checkCharge = cms.bool(False)
                                    )

process.McMatchRECO = cms.EDProducer("MCTruthDeltaRMatcherNew",
                                     matchPDGId = cms.vint32(11),
                                     src = cms.InputTag("goodPhotons"),
                                     distMin = cms.double(0.2),
                                     matched = cms.InputTag("prunedGenParticles"),
                                     checkCharge = cms.bool(False)
                                    )

process.mc_sequence = cms.Sequence()

if (varOptions.isMC):
    process.mc_sequence *= (process.McMatchTag + process.McMatchRECO)

##########################################################################
## TREE CONTENT
#########################################################################

ZVariablesToStore = cms.PSet(
    eta    = cms.string("eta"),
    abseta = cms.string("abs(eta)"),
    pt     = cms.string("pt"),
    mass   = cms.string("mass")
    )   

ProbeVariablesToStore = cms.PSet(
    probe_Pho_eta    = cms.string("eta"),
    probe_Pho_abseta = cms.string("abs(eta)"),
    probe_Pho_et     = cms.string("et"),
    probe_Pho_e      = cms.string("energy"),
    probe_Pho_sieie  = cms.string("full5x5_sigmaIetaIeta"),
## super cluster quantities
    probe_sc_energy = cms.string("superCluster.energy"),
    probe_sc_et     = cms.string("superCluster.energy*sin(superCluster.position.theta)"),    
    probe_sc_eta    = cms.string("superCluster.eta"),
    probe_sc_abseta = cms.string("abs(superCluster.eta)"),

#id based
    probe_Pho_sigmaIEtaIEta = cms.string("full5x5_sigmaIetaIeta"),
    probe_Pho_ESsigma       = cms.InputTag("photonIDValueMapProducer:phoESEffSigmaRR"),
    probe_Pho_sigmaIEtaIPhi = cms.InputTag("photonIDValueMapProducer:phoFull5x5SigmaIEtaIPhi"),
    probe_Pho_hoe           = cms.string("hadronicOverEm"),

#iso
    probe_Pho_chIso    = cms.InputTag("photonIDValueMapProducer:phoChargedIsolation"),
    probe_Pho_neuIso   = cms.InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation"),
    probe_Pho_phoIso   = cms.InputTag("photonIDValueMapProducer:phoPhotonIsolation"),
    probe_Pho_chWorIso = cms.InputTag("photonIDValueMapProducer:phoWorstChargedIsolation"), 
    
    probe_dRTau    = cms.InputTag("GsfDRToNearestTauProbe"),
)

TagVariablesToStore = cms.PSet(
    Pho_eta    = cms.string("eta"),
    Pho_abseta = cms.string("abs(eta)"),
    Pho_pt     = cms.string("pt"),
    Pho_et     = cms.string("et"),
    Pho_e      = cms.string("energy"),
    Pho_dRTau  = cms.InputTag("GsfDRToNearestTauProbe"),
    
    ## super cluster quantities
    sc_energy = cms.string("superCluster.energy"),
    sc_et     = cms.string("superCluster.energy*sin(superCluster.position.theta)"),    
    sc_eta    = cms.string("superCluster.eta"),
    sc_abseta = cms.string("abs(superCluster.eta)"),
)

CommonStuffForPhotonProbe = cms.PSet(
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

if varOptions.isMC:
    mcTruthCommonStuff = cms.PSet(
        isMC        = cms.bool(True),
        tagMatches  = cms.InputTag("McMatchTag"),
        motherPdgId = cms.vint32(22,23),
        #motherPdgId = cms.vint32(443), # JPsi
        #motherPdgId = cms.vint32(553), # Yupsilon
        makeMCUnbiasTree       = cms.bool(False),
        checkMotherInUnbiasEff = cms.bool(False),
        mcVariables = cms.PSet(
            probe_eta    = cms.string("eta"),
            probe_abseta = cms.string("abs(eta)"),
            probe_et     = cms.string("et"),
            probe_e      = cms.string("energy"),
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
## TREE MAKER
##########################################################################

process.PhotonToRECO = cms.EDAnalyzer("TagProbeFitTreeProducer",
                                      mcTruthCommonStuff, CommonStuffForPhotonProbe,
                                      tagProbePairs = cms.InputTag("tagTightRECO"),
                                      arbitration   = cms.string("None"),
                                      flags         = cms.PSet(passingLoose  = cms.InputTag("goodPhotonsPROBECutBasedLoose"),
                                                               passingMedium = cms.InputTag("goodPhotonsPROBECutBasedMedium"),
                                                               passingTight  = cms.InputTag("goodPhotonsPROBECutBasedTight"),
                                                               ),                                               
                                      allProbes     = cms.InputTag("goodPhotonsProbeHLT"),
                                      )

if (varOptions.isMC):
    process.PhotonToRECO.probeMatches  = cms.InputTag("McMatchRECO")
    process.PhotonToRECO.eventWeight   = cms.InputTag("generator")
    process.PhotonToRECO.PUWeightSrc   = cms.InputTag("pileupReweightingProducer","pileupWeights")

process.tree_sequence = cms.Sequence(process.PhotonToRECO)

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
        process.pho_sequence + 
        process.allTagsAndProbes +
        process.pileupReweightingProducer +
        process.mc_sequence + 
        process.GsfDRToNearestTauProbe + 
        process.GsfDRToNearestTauTag + 
        process.tree_sequence
        )
else:
    process.p = cms.Path(
        process.hltHighLevel +
        process.pho_sequence + 
        process.allTagsAndProbes +
        process.mc_sequence +
        process.tree_sequence
        )

process.TFileService = cms.Service(
    "TFileService", fileName = cms.string(options['OUTPUT_FILE_NAME'])
    )
