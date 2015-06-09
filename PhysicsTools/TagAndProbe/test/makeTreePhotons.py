import FWCore.ParameterSet.Config as cms
import sys

from PhysicsTools.TagAndProbe.treeMakerOptionsPhotonsMC_cfi import options
#from PhysicsTools.TagAndProbe.treeMakerOptionsData_cfi import options

process = cms.Process("tnp")

process.pileupReweightingProducer = cms.EDProducer("PileupWeightProducer",
                                                   FirstTime = cms.untracked.bool(True)
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
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
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

process.goodPhotons = cms.EDFilter("PATPhotonRefSelector",
                                    src = cms.InputTag(options['PHOTON_COLL']),
                                    cut = cms.string(options['PHOTON_CUTS'])    
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

switchOnVIDPhotonIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_PHYS14_PU20bx25_V2_cff',]

for idmod in my_id_modules:
    setupAllVIDIdsInModule(process, idmod, setupVIDPhotonSelection)

process.goodPhotonsPROBECutBasedLoose = cms.EDProducer("PatPhotonSelectorByValueMap",
                                                       input     = cms.InputTag("goodPhotons"),
                                                       cut       = cms.string(options['PHOTON_CUTS']),
                                                       selection = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-loose"),
                                                       id_cut    = cms.bool(True)
                                                       )

process.goodPhotonsPROBECutBasedMedium = process.goodPhotonsPROBECutBasedLoose.clone()
process.goodPhotonsPROBECutBasedMedium.selection = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-medium")
process.goodPhotonsPROBECutBasedTight = process.goodPhotonsPROBECutBasedLoose.clone()
process.goodPhotonsPROBECutBasedTight.selection = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-tight")

process.goodPhotonsTAGCutBasedLoose = cms.EDProducer("PatPhotonSelectorByValueMap",
                                                     input     = cms.InputTag("goodPhotons"),
                                                     cut       = cms.string(options['PHOTON_TAG_CUTS']),
                                                     selection = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-loose"),
                                                     id_cut    = cms.bool(True)
                                                     )

process.goodPhotonsTAGCutBasedMedium = process.goodPhotonsTAGCutBasedLoose.clone()
process.goodPhotonsTAGCutBasedMedium.selection = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-medium")
process.goodPhotonsTAGCutBasedTight = process.goodPhotonsTAGCutBasedLoose.clone()
process.goodPhotonsTAGCutBasedTight.selection = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-tight")

###################################################################
##    _____     _                         __  __       _       _     _             
##   |_   _| __(_) __ _  __ _  ___ _ __  |  \/  | __ _| |_ ___| |__ (_)_ __   __ _ 
##     | || '__| |/ _` |/ _` |/ _ \ '__| | |\/| |/ _` | __/ __| '_ \| | '_ \ / _` |
##     | || |  | | (_| | (_| |  __/ |    | |  | | (_| | || (__| | | | | | | | (_| |
##     |_||_|  |_|\__, |\__, |\___|_|    |_|  |_|\__,_|\__\___|_| |_|_|_| |_|\__, |
##                |___/ |___/                                                |___/ 
###################################################################

if (len(options['TnPHLTTagFilters']) != len(options['TnPHLTProbeFilters'])):
    print "ERROR: different number of tag and probe filters, please fix it"
    sys.exit(-1)

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
##    _____ ___   ____    ____       _          
##   |_   _( _ ) |  _ \  |  _ \ __ _(_)_ __ ___ 
##     | | / _ \/\ |_) | | |_) / _` | | '__/ __|
##     | || (_>  <  __/  |  __/ (_| | | |  \__ \
##     |_| \___/\/_|     |_|   \__,_|_|_|  |___/
##   
###################################################################

process.tagTightRECO = cms.EDProducer("CandViewShallowCloneCombiner",
                                      decay = cms.string("goodPhotonsTagHLT@+ goodPhotonsProbeHLT@-"), 
                                      checkCharge = cms.bool(False),
                                      cut = cms.string("40<mass<1000"),
                                      )

process.allTagsAndProbes = cms.Sequence(process.tagTightRECO)

###################################################################
##    __  __  ____   __  __       _       _               
##   |  \/  |/ ___| |  \/  | __ _| |_ ___| |__   ___  ___ 
##   | |\/| | |     | |\/| |/ _` | __/ __| '_ \ / _ \/ __|
##   | |  | | |___  | |  | | (_| | || (__| | | |  __/\__ \
##   |_|  |_|\____| |_|  |_|\__,_|\__\___|_| |_|\___||___/
##                                   
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

if (options['MC_FLAG']):
    process.mc_sequence *= (process.McMatchTag + process.McMatchRECO)

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
    probe_Pho_sigmaIEtaIEta = cms.string("sigmaIetaIeta"),
    probe_Pho_hoe           = cms.string("hadronicOverEm"),
)

TagVariablesToStore = cms.PSet(
    Pho_eta    = cms.string("eta"),
    Pho_abseta = cms.string("abs(eta)"),
    Pho_pt     = cms.string("pt"),
    Pho_et     = cms.string("et"),
    Pho_e      = cms.string("energy"),
    
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
    pairVariables =  cms.PSet(ZVariablesToStore),
    pairFlags     =  cms.PSet(
        mass60to120 = cms.string("60 < mass < 120")
        ),
    tagVariables   =  cms.PSet(TagVariablesToStore),
    tagFlags       =  cms.PSet(),    
    )

if options['MC_FLAG']:
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
##   ____      __       __    ___                 ___    _ 
##  / ___|___ / _|      \ \  |_ _|___  ___       |_ _|__| |
## | |  _/ __| |_   _____\ \  | |/ __|/ _ \       | |/ _` |
## | |_| \__ \  _| |_____/ /  | |\__ \ (_) |  _   | | (_| |
##  \____|___/_|        /_/  |___|___/\___/  ( ) |___\__,_|
##                                           |/            
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
                                      PUWeightSrc   = cms.InputTag("pileupReweightingProducer","pileupWeights")
                                      )

if (options['MC_FLAG']):
    process.PhotonToRECO.probeMatches  = cms.InputTag("McMatchRECO")

process.tree_sequence = cms.Sequence(process.PhotonToRECO)

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

process.p = cms.Path(
    process.hltHighLevel +
    process.pho_sequence + 
    ####process.GsfDRToNearestTau+
    process.allTagsAndProbes +
    process.pileupReweightingProducer +
    process.mc_sequence +
    process.tree_sequence
    )

process.TFileService = cms.Service(
    "TFileService", fileName = cms.string(options['OUTPUT_FILE_NAME'])
    )
