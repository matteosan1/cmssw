import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

################################################
##                      _              _       
##   ___ ___  _ __  ___| |_ __ _ _ __ | |_ ___ 
##  / __/ _ \| '_ \/ __| __/ _` | '_ \| __/ __|
## | (_| (_) | | | \__ \ || (_| | | | | |_\__ \
##  \___\___/|_| |_|___/\__\__,_|_| |_|\__|___/
##                                              
################################################

isMC = True
InputFileName = "prova.root"
OutputFilePrefix = "efficiency-data-"
PDFName = "pdfSignalPlusBackground"

if isMC:
    InputFileName = "TnPTree_mc.root"
    PDFName = "pdfSignalPlusBackground"
    OutputFilePrefix = "efficiency-mc-pt-high-"

################################################
#specifies the binning of parameters
EfficiencyBins = cms.PSet(probe_et = cms.vdouble(20, 1000 ),
                          probe_abseta = cms.vdouble(0, 1.5, 2.5)
                          )

#### For data: except for HLT step
EfficiencyBinningSpecification = cms.PSet(
    #specifies what unbinned variables to include in the dataset, the mass is needed for the fit
    UnbinnedVariables = cms.vstring("mass", "totWeight", "Ele_dRTau", "probe_dRTau"),
    #specifies the binning of parameters
    BinnedVariables = cms.PSet(EfficiencyBins),
    #first string is the default followed by binRegExp - PDFname pairs
    BinToPDFmap = cms.vstring(PDFName)
)

#### For MC truth: do truth matching
EfficiencyBinningSpecificationMC = cms.PSet(
    UnbinnedVariables = cms.vstring("mass", "totWeight", "Ele_dRTau", "probe_dRTau"),
    BinnedVariables = cms.PSet(EfficiencyBins,
                               mcTrue = cms.vstring("true")
                               ),
    BinToPDFmap = cms.vstring(PDFName)  
)

############################################################################################

if isMC:
    mcTruthModules = cms.PSet(
        MCtruth_Reco = cms.PSet(EfficiencyBinningSpecificationMC,
                                EfficiencyCategoryAndState = cms.vstring("passingRECO", "pass"),
                                ),
        )
else:
    mcTruthModules = cms.PSet()

############################################################################################
############################################################################################
####### GsfElectron->Id / selection efficiency 
############################################################################################
############################################################################################

process.GsfElectronToRECO = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
                                           InputFileNames = cms.vstring(InputFileName),
                                           InputDirectoryName = cms.string("GsfElectronToSC"),
                                           InputTreeName = cms.string("fitter_tree"),
                                           OutputFileName = cms.string(OutputFilePrefix+"GsfElectronToRECO.root"),
                                           NumCPU = cms.uint32(8),
                                           SaveWorkspace = cms.bool(False),
                                           doCutAndCount = cms.bool(False),
                                           floatShapeParameters = cms.bool(True),
                                           binnedFit = cms.bool(True),
                                           binsForFit = cms.uint32(60),
                                           WeightVariable = cms.string("totWeight"),
                                           #fixVars = cms.vstring("mean"),
                                           
                                           # defines all the real variables of the probes available in the input tree and intended for use in the efficiencies
                                           Variables = cms.PSet(mass = cms.vstring("Tag-Probe Mass", "60.0", "120.0", "GeV/c^{2}"),
                                                                probe_et = cms.vstring("Probe E_{T}", "0", "1000", "GeV/c"),
                                                                probe_abseta = cms.vstring("Probe #eta", "0", "2.5", ""),                
                                                                totWeight = cms.vstring("totWeight", "0", "10000000", ""),         
                                                                Ele_dRTau = cms.vstring("Ele_dRTau", "0.2", "100000", ""),
                                                                probe_dRTau = cms.vstring("probe_dRTau", "0.2", "100000", ""),       
                                                                ),
                                           
                                           # defines all the discrete variables of the probes available in the input tree and intended for use in the efficiency calculations
                                           Categories = cms.PSet(mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
                                                                 passingRECO = cms.vstring("passingRECO", "dummy[pass=1,fail=0]"),
                                                                 ),
                                           
                                           # defines all the PDFs that will be available for the efficiency calculations; 
                                           # uses RooFit's "factory" syntax;
                                           # each pdf needs to define "signal", "backgroundPass", "backgroundFail" pdfs, "efficiency[0.9,0,1]" 
                                           # and "signalFractionInPassing[0.9]" are used for initial values  
                                           PDFs = cms.PSet(pdfSignalPlusBackground = cms.vstring(
            "RooCBExGaussShape::signalResPass(mass,meanP[-0.0,-5.000,5.000],sigmaP[0.956,0.001,5.000],alphaP[0.999,0.000,5.00],nP[1.405,0.000,5.000],sigmaP_2[2.000,0.00,5.00])",
            "RooCBExGaussShape::signalResFail(mass,meanF[-0.0,-5.000,5.000],sigmaF[3.331,0.001,5.000],alphaF[1.586,0.000,5.0],nF[2.464,0.000,5.0],sigmaF_2[1.675,0.000,5.0])",
            "ZGeneratorLineShape::signalPhy(mass)", ### NLO line shape
            "RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,70.], betaPass[0.001, 0.,0.1], gammaPass[0.001, 0.,0.1], peakPass[90.0])",
            "RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,70.], betaFail[0.001, 0.,0.1], gammaFail[0.001, 0.,0.1], peakFail[90.0])",
            "FCONV::signalPass(mass, signalPhy, signalResPass)",
            "FCONV::signalFail(mass, signalPhy, signalResFail)",     
            "efficiency[0.5,0,1]",
            "signalFractionInPassing[1.0]"     
            ),
                                                           ),
                                           
                                           # defines a set of efficiency calculations, what PDF to use for fitting and how to bin the data;
                                           # there will be a separate output directory for each calculation that includes a simultaneous fit, side band subtraction and counting. 
                                           Efficiencies = cms.PSet(mcTruthModules,
                                                                   #the name of the parameter set becomes the name of the directory
                                                                   #Reco = cms.PSet(EfficiencyBinningSpecification,
                                                                   #                EfficiencyCategoryAndState = cms.vstring("passingRECO", "pass"),
                                                                   #                ),
                                                                   )
                                           )

process.fit = cms.Path(
    process.GsfElectronToRECO
    )
