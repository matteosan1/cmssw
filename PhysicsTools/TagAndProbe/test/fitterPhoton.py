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
    InputFileName = "TnPTree.root"
    PDFName = "pdfSignalPlusBackground"
    OutputFilePrefix = "efficiency-mc-"

################################################
#specifies the binning of parameters
EfficiencyBins = cms.PSet(probe_sc_et = cms.vdouble( 25, 40 ),
                          probe_sc_abseta = cms.vdouble( 0.0, 1.5, 2.5 )
                          )

#### For data: except for HLT step
EfficiencyBinningSpecification = cms.PSet(
    #specifies what unbinned variables to include in the dataset, the mass is needed for the fit
    UnbinnedVariables = cms.vstring("mass"),
    #specifies the binning of parameters
    BinnedVariables = cms.PSet(EfficiencyBins),
    #first string is the default followed by binRegExp - PDFname pairs
    BinToPDFmap = cms.vstring(PDFName)
)

#### For MC truth: do truth matching
EfficiencyBinningSpecificationMC = cms.PSet(
    UnbinnedVariables = cms.vstring("mass"),
    BinnedVariables = cms.PSet(EfficiencyBins,
                               mcTrue = cms.vstring("true")
                               ),
    BinToPDFmap = cms.vstring(PDFName)  
)

############################################################################################

if isMC:
    mcTruthModules = cms.PSet(
        MCtruth_Medium = cms.PSet(EfficiencyBinningSpecificationMC,
                                  EfficiencyCategoryAndState = cms.vstring("passingMedium", "pass"),
                                  ),
        )
else:
    mcTruthModules = cms.PSet()

############################################################################################
############################################################################################
####### GsfElectron->Id / selection efficiency 
############################################################################################
############################################################################################

process.PhotonToId = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
                                    InputFileNames = cms.vstring(InputFileName),
                                    InputDirectoryName = cms.string("PhotonToRECO"),
                                    InputTreeName = cms.string("fitter_tree"),
                                    OutputFileName = cms.string(OutputFilePrefix+"PhotonToId.root"),
                                    NumCPU = cms.uint32(8),
                                    SaveWorkspace = cms.bool(False),
                                    doCutAndCount = cms.bool(False),
                                    floatShapeParameters = cms.bool(True),
                                    binnedFit = cms.bool(False),
                                    binsForFit = cms.uint32(60),
                                    #WeightVariable = cms.string("PUweight"),
                                    #fixVars = cms.vstring("mean"),
                                    
                                    # defines all the real variables of the probes available in the input tree and intended for use in the efficiencies
                                    Variables = cms.PSet(mass = cms.vstring("Tag-Probe Mass", "60.0", "120.0", "GeV/c^{2}"),
                                                         probe_sc_et = cms.vstring("Probe E_{T}", "0", "1000", "GeV/c"),
                                                         probe_sc_abseta = cms.vstring("Probe #eta", "0", "2.5", ""),                
                                                         ),
                                    
                                    # defines all the discrete variables of the probes available in the input tree and intended for use in the efficiency calculations
                                    Categories = cms.PSet(mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
                                                          #probe_passConvRej = cms.vstring("probe_passConvRej", "dummy[pass=1,fail=0]"), 
                                                          passingMedium = cms.vstring("passingMedium", "dummy[pass=1,fail=0]"),
                                                          ),
                                    
                                    # defines all the PDFs that will be available for the efficiency calculations; 
                                    # uses RooFit's "factory" syntax;
                                    # each pdf needs to define "signal", "backgroundPass", "backgroundFail" pdfs, "efficiency[0.9,0,1]" 
                                    # and "signalFractionInPassing[0.9]" are used for initial values  
                                    PDFs = cms.PSet(pdfSignalPlusBackground = cms.vstring(
            "RooCBExGaussShape::signalResPass(mass, meanP[0, -5., 5.], sigmaP[1.5, 1., 50.],alphaP[0.01, 0, 5], nP[.6, 0, 20], sigmaP_2[2, 1., 40.], fracP[0.6,0, 1])",
            "RooCBExGaussShape::signalResFail(mass, meanF[0., -5., 5.], sigmaF[1.5, 1., 50.],alphaF[0.01, 0, 5], nF[.6, 0, 20], sigmaF_2[2, 1., 40.], fracF[0.6, 0, 1])",
            #"RooCBExGaussShape::signalPass(mass, meanP[90, 85., 95.], sigmaP[1.5, 1.4, 50.],alphaP[0.01, 0, 5], nP[.6, 0, 20], sigmaP_2[2, 1., 40.], fracP[0.6,0, 1])",
            #"RooCBExGaussShape::signalFail(mass, meanF[90., 85., 95.], sigmaF[1.5, 1.4, 50.],alphaF[0.01, 0, 5], nF[.6, 0, 20], sigmaF_2[2, 1., 40.], fracF[0.6, 0, 1])",
            "ZGeneratorLineShape::signalPhy(mass)", ### NLO line shape
            "RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,70.], betaPass[0.001, 0.,0.1], betaPass, peakPass[90.0])",
            "RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,70.], betaFail[0.001, 0.,0.1], betaFail, peakFail[90.0])",
            "FCONV::signalPass(mass, signalPhy, signalResPass)",
            "FCONV::signalFail(mass, signalPhy, signalResFail)",     
            "efficiency[0.5,0,1]",
            "signalFractionInPassing[1.0]"     
            ),
                                                    ),
                                    
                                    # defines a set of efficiency calculations, what PDF to use for fitting and how to bin the data;
                                    # there will be a separate output directory for each calculation that includes a simultaneous fit, side band subtraction and counting. 
                                    Efficiencies = cms.PSet(mcTruthModules,
                                                            #Medium = cms.PSet(EfficiencyBinningSpecification,
                                                            #                  EfficiencyCategoryAndState = cms.vstring("passingMedium", "pass"),
                                                            #                  ),
                                                            )
                                    )


process.fit = cms.Path(
    process.PhotonToId  
    )
