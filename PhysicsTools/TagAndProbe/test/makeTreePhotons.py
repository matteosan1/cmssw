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

options['HLTProcessName']        = "HLT"
options['PHOTON_COLL']           = "slimmedPhotons"
options['PHOTON_CUTS']           = "(abs(superCluster.eta)<2.5) && ((superCluster.energy*sin(superCluster.position.theta))>15.0)"
options['PHOTON_TAG_CUTS']       = "(abs(superCluster.eta)<=2.5) && !(1.4442<=abs(superCluster.eta)<=1.566) && (superCluster.energy*sin(superCluster.position.theta))>25.0"
options['MAXEVENTS']             = cms.untracked.int32(1000) 
options['useAOD']                = cms.bool(False)
options['OUTPUTEDMFILENAME']     = 'edmFile.root'
options['DEBUG']                 = cms.bool(False)

from PhysicsTools.TagAndProbe.treeMakerOptionsPhotons_cfi import *

if (varOptions.isMC):
    options['INPUT_FILE_NAME']       = ("/store/relval/CMSSW_7_4_1/RelValZEE_13/MINIAODSIM/MCRUN2_74_V9_gensim_740pre7-v1/00000/1E35CCF8-32EC-E411-8F29-0025905A48D0.root",)
    options['OUTPUT_FILE_NAME']      = "TnPTree_mc.root"
    options['TnPPATHS']            = cms.vstring("HLT_Ele23_WP75_Gsf_v*")
    options['TnPHLTTagFilters']    = cms.vstring("hltEle23WP75GsfTrackIsoFilter")
    options['TnPHLTProbeFilters']  = cms.vstring()
    options['HLTFILTERTOMEASURE']  = cms.vstring("")
     options['GLOBALTAG']             = 'MCRUN2_74_V9'
    options['EVENTSToPROCESS']       = cms.untracked.VEventRange()
else:
    options['INPUT_FILE_NAME']       = ("/store/relval/CMSSW_7_4_1/RelValZEE_13/MINIAODSIM/MCRUN2_74_V9_gensim_740pre7-v1/00000/1E35CCF8-32EC-E411-8F29-0025905A48D0.root")
    options['OUTPUT_FILE_NAME']      = "TnPTree_data.root"
    options['TnPPATHS']            = ["HLT_Ele23_WPLoose_Gsf_v1",]
    options['TnPHLTTagFilters']    = ["hltEle23WPLooseGsfTrackIsoFilter"]
    options['TnPHLTProbeFilters']  = cms.vstring()
    options['HLTFILTERTOMEASURE']  = cms.vstring("")
    options['GLOBALTAG']             = 'MCRUN2_74_V9'
    options['EVENTSToPROCESS']       = cms.untracked.VEventRange()

###################################################################

setModules(process, options)
from PhysicsTools.TagAndProbe.treeContentPhotons_cfi import *

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

from PhysicsTools.TagAndProbe.photonIDModules_cfi import *
setIDs(process, options)

###################################################################
## SEQUENCES
###################################################################

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

process.allTagsAndProbes = cms.Sequence()
process.allTagsAndProbes *= process.tagTightRECO

process.mc_sequence = cms.Sequence()

if (varOptions.isMC):
    process.mc_sequence *= (process.McMatchTag + process.McMatchRECO)

##########################################################################
## TREE MAKER OPTIONS
##########################################################################
if (not varOptions.isMC):
    mcTruthCommonStuff = cms.PSet(
        isMC = cms.bool(False)
        )

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
    process.PhotonToRECO.variables.Pho_dRTau  = cms.InputTag("GsfDRToNearestTauProbe")
    process.PhotonToRECO.tagVariables.probe_dRTau    = cms.InputTag("GsfDRToNearestTauProbe")

process.tree_sequence = cms.Sequence(process.PhotonToRECO)

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
        process.sampleInfo +
        process.hltFilter +
        process.pho_sequence + 
        process.allTagsAndProbes +
        process.mc_sequence +
        process.tree_sequence
        )

process.TFileService = cms.Service(
    "TFileService", fileName = cms.string(options['OUTPUT_FILE_NAME']),
    closeFileFast = cms.untracked.bool(True)
    )
