import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import sys

process = cms.Process("tnp")

###################################################################
options = dict()
varOptions = VarParsing('analysis')
varOptions.register(
    "isMC",
    True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Compute MC efficiencies"
    )

varOptions.parseArguments()

options['HLTProcessName']          = "HLT"
options['ELECTRON_COLL']           = "slimmedElectrons"
options['ELECTRON_CUTS']           = "(abs(superCluster.eta)<2.5) && (ecalEnergy*sin(superClusterPosition.theta)>10.0)"
options['ELECTRON_TAG_CUTS']       = "(abs(superCluster.eta)<=2.5) && !(1.4442<=abs(superCluster.eta)<=1.566) && pt >= 25.0"
options['SUPERCLUSTER_COLL']       = "reducedEgamma:reducedSuperClusters"
options['SUPERCLUSTER_CUTS']       = "abs(eta)<2.5 && !(1.4442< abs(eta) <1.566) && et>10.0"
options['MAXEVENTS']               = cms.untracked.int32(-1) 
options['useAOD']                  = cms.bool(False)
options['DOTRIGGER']               = cms.bool(True)
options['DORECO']                  = cms.bool(True)
options['DOID']                    = cms.bool(True)
options['OUTPUTEDMFILENAME']       = 'edmFile.root'
options['DEBUG']                   = cms.bool(False)

from PhysicsTools.TagAndProbe.treeMakerOptions_cfi import *

if (varOptions.isMC):
    options['INPUT_FILE_NAME']     = "/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/00C4781D-6B08-E511-8A0A-0025905A6084.root"
    options['OUTPUT_FILE_NAME']    = "TnPTree_mc.root"
    options['TnPPATHS']            = cms.vstring("HLT_Ele23_WP75_Gsf_v*")
    options['TnPHLTTagFilters']    = cms.vstring("hltEle23WP75GsfTrackIsoFilter")
    options['TnPHLTProbeFilters']  = cms.vstring()
    options['HLTFILTERTOMEASURE']  = cms.vstring("")
    options['GLOBALTAG']           = 'MCRUN2_74_V9A'
    options['EVENTSToPROCESS']     = cms.untracked.VEventRange()
else:
    options['INPUT_FILE_NAME']     = "/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/244/00000/12EE24E2-8F27-E511-80D1-02163E013793.root"
    options['OUTPUT_FILE_NAME']    = "TnPTree_data.root"
    options['TnPPATHS']            = ["HLT_Ele23_WPLoose_Gsf_v1",]
    options['TnPHLTTagFilters']    = ["hltEle23WPLooseGsfTrackIsoFilter"]
    options['TnPHLTProbeFilters']  = cms.vstring()
    options['HLTFILTERTOMEASURE']  = cms.vstring("")
    options['GLOBALTAG']           = 'GR_P_V56'
    options['EVENTSToPROCESS']     = cms.untracked.VEventRange()

###################################################################

setModules(process, options)
from PhysicsTools.TagAndProbe.treeContent_cfi import *

process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = options['GLOBALTAG']

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
## ID
###################################################################

from PhysicsTools.TagAndProbe.electronIDModules_cfi import *
setIDs(process, options)
process.goodElectronsPROBECutBasedVeto.saveSCRef = cms.bool(True)
process.goodElectronsPROBECutBasedVeto.recoEcalCandidates = cms.InputTag("GsfMatchedSuperClusterCands","superclusters")
process.goodElectronsPROBECutBasedLoose.saveSCRef = cms.bool(True)
process.goodElectronsPROBECutBasedLoose.recoEcalCandidates = cms.InputTag("GsfMatchedSuperClusterCands","superclusters")
process.goodElectronsPROBECutBasedMedium.saveSCRef = cms.bool(True)
process.goodElectronsPROBECutBasedMedium.recoEcalCandidates = cms.InputTag("GsfMatchedSuperClusterCands","superclusters")
process.goodElectronsPROBECutBasedTight.saveSCRef = cms.bool(True)
process.goodElectronsPROBECutBasedTight.recoEcalCandidates = cms.InputTag("GsfMatchedSuperClusterCands","superclusters")

###################################################################
## SEQUENCES
###################################################################

process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag(options['ELECTRON_COLL'])
process.ele_sequence = cms.Sequence(
    process.goodElectrons +
    process.egmGsfElectronIDSequence +
    process.goodElectronsTAGCutBasedVeto +
    process.goodElectronsTAGCutBasedLoose +
    process.goodElectronsTAGCutBasedMedium +
    process.goodElectronsTAGCutBasedTight +
    process.goodElectronsTagHLT
    )

process.sc_sequence = cms.Sequence(process.superClusterCands +
                                   process.goodSuperClusters +
                                   process.goodSuperClustersHLT +
                                   process.GsfMatchedSuperClusterCands +
                                   process.goodElectronsPROBECutBasedVeto +
                                   process.goodElectronsPROBECutBasedLoose +
                                   process.goodElectronsPROBECutBasedMedium +
                                   process.goodElectronsPROBECutBasedTight 
                                   )

###################################################################
## TnP PAIRS
###################################################################
process.tagTightRECO.decay = cms.string("goodElectronsTagHLT goodSuperClustersHLT") 
process.allTagsAndProbes = cms.Sequence()
process.allTagsAndProbes *= process.tagTightRECO

process.mc_sequence = cms.Sequence()

#if (varOptions.isMC):
#    process.mc_sequence *= (process.McMatchTag + process.McMatchSC)

##########################################################################
## TREE MAKER OPTIONS
##########################################################################
if (not varOptions.isMC):
    mcTruthCommonStuff = cms.PSet(
        isMC = cms.bool(False)
        )

process.GsfElectronToSC = cms.EDAnalyzer("TagProbeFitTreeProducer",
                                         CommonStuffForSuperClusterProbe, mcTruthCommonStuff,
                                         tagProbePairs = cms.InputTag("tagTightRECO"), 
                                         arbitration   = cms.string("Random2"),
                                         flags         = cms.PSet(passingVeto   = cms.InputTag("goodElectronsPROBECutBasedVeto", "superclusters"),
                                                                  passingLoose  = cms.InputTag("goodElectronsPROBECutBasedLoose","superclusters"),
                                                                  passingMedium = cms.InputTag("goodElectronsPROBECutBasedMedium","superclusters"),
                                                                  passingTight  = cms.InputTag("goodElectronsPROBECutBasedTight","superclusters"),
                                                                  passingRECO   = cms.InputTag("GsfMatchedSuperClusterCands","superclusters")
                                                                  ),                                               
                                         allProbes     = cms.InputTag("goodSuperClustersHLT"),
                                         )

if (varOptions.isMC):
    #process.GsfElectronToSC.probeMatches  = cms.InputTag("McMatchSC")
    process.GsfElectronToSC.eventWeight   = cms.InputTag("generator")
    process.GsfElectronToSC.PUWeightSrc   = cms.InputTag("pileupReweightingProducer","pileupWeights")
    #process.GsfElectronToSC.variables.probe_dRTau = cms.InputTag("GsfDRToNearestTauSC")
    #process.GsfElectronToSC.tagVariables.Ele_dRTau = cms.InputTag("GsfDRToNearestTauTag")

##########################################################################
## PATHS
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
        process.sampleInfo +
        process.hltFilter +
        process.ele_sequence + 
        process.sc_sequence +
        process.allTagsAndProbes +
        process.pileupReweightingProducer +
        process.mc_sequence +
        process.eleVarHelper +
        #process.GsfDRToNearestTauProbe + 
        #process.GsfDRToNearestTauTag + 
        #process.GsfDRToNearestTauSC + 
        process.GsfElectronToSC
        )
else:
    process.p = cms.Path(
        process.sampleInfo +
        process.hltFilter +
        process.ele_sequence + 
        process.sc_sequence +
        process.allTagsAndProbes +
        process.mc_sequence +
        process.eleVarHelper +
        process.GsfElectronToSC
        )

process.TFileService = cms.Service(
    "TFileService", fileName = cms.string(options['OUTPUT_FILE_NAME']),
    closeFileFast = cms.untracked.bool(True)
    )
