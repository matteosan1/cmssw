import FWCore.ParameterSet.Config as cms

###################################################################
## ID MODULES
###################################################################

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

def setIDs(process, options):

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

    process.goodElectronsPROBECutBasedVeto = cms.EDProducer("GsfElectronSelectorByValueMap",
                                                            input     = cms.InputTag("goodElectrons"),
                                                            cut       = cms.string(options['ELECTRON_CUTS']),
                                                            selection = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
                                                            id_cut    = cms.bool(True)
                                                            )

    ######################################################################################
    ## MODULES FOR N-1 CUT BASED STUDIES
    ######################################################################################
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
    
    #process.goodElectronsPROBECutBasedVeto = cms.EDProducer("PatElectronNm1Selector",
    #                                                    input     = cms.InputTag("goodElectrons"),
    #                                                    cut       = cms.string(options['ELECTRON_CUTS']),
    #                                                    selection = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
    #                                                    #cutIndicesToMask = cms.vuint32(2, 3),
    #                                                    cutNamesToMask = cms.vstring("GsfEleDEtaInCut_0", "GsfEleDPhiInCut_0")
    #                                                    )
    
    process.goodElectronsPROBECutBasedLoose = process.goodElectronsPROBECutBasedVeto.clone()
    process.goodElectronsPROBECutBasedLoose.selection = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose")
    process.goodElectronsPROBECutBasedMedium = process.goodElectronsPROBECutBasedVeto.clone()
    process.goodElectronsPROBECutBasedMedium.selection = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium")
    process.goodElectronsPROBECutBasedTight = process.goodElectronsPROBECutBasedVeto.clone()
    process.goodElectronsPROBECutBasedTight.selection = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight")
    
    process.goodElectronsTAGCutBasedVeto = cms.EDProducer("GsfElectronSelectorByValueMap",
                                                          input     = cms.InputTag("goodElectrons"),
                                                          cut       = cms.string(options['ELECTRON_TAG_CUTS']),
                                                          selection = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
                                                          id_cut    = cms.bool(True)
                                                          )
    
    process.goodElectronsTAGCutBasedLoose = process.goodElectronsTAGCutBasedVeto.clone()
    process.goodElectronsTAGCutBasedLoose.selection = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose")
    process.goodElectronsTAGCutBasedMedium = process.goodElectronsTAGCutBasedVeto.clone()
    process.goodElectronsTAGCutBasedMedium.selection = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium")
    process.goodElectronsTAGCutBasedTight = process.goodElectronsTAGCutBasedVeto.clone()
    process.goodElectronsTAGCutBasedTight.selection = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight")
