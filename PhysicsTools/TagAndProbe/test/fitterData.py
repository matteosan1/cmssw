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

isMC = False
InputFileName = "TnPTree_data.root"
OutputFilePrefix = "efficiency-data-bin0-tight"
PDFName = "pdfSignalPlusBackground"

if isMC:
    InputFileName = "TnPTree_mc_11.root"
    PDFName = "pdfSignalPlusBackground"
    OutputFilePrefix = "efficiency-mc-"

################################################
#specifies the binning of parameters
EfficiencyBins = cms.PSet(probe_sc_et = cms.vdouble( 20, 1000 ),
                          probe_sc_abseta = cms.vdouble( 0., 1.5 ),
                          #probe_sc_abseta = cms.vdouble( 1.5, 2.5 ),
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

process.GsfElectronToId = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
                                         InputFileNames = cms.vstring(InputFileName),
                                         InputDirectoryName = cms.string("GsfElectronToRECO"),
                                         InputTreeName = cms.string("fitter_tree"), 
                                         OutputFileName = cms.string(OutputFilePrefix+"GsfElectronToId.root"),
                                         NumCPU = cms.uint32(1),
                                         SaveWorkspace = cms.bool(False),
                                         doCutAndCount = cms.bool(False),
                                         floatShapeParameters = cms.bool(True),
                                         binnedFit = cms.bool(True),
                                         binsForFit = cms.uint32(60),
                                         #fixVars = cms.vstring("meanP", "meanF", "sigmaP", "sigmaF", "sigmaP_2", "sigmaF_2"),
                                         
                                         # defines all the real variables of the probes available in the input tree and intended for use in the efficiencies
                                         Variables = cms.PSet(mass = cms.vstring("Tag-Probe Mass", "60.0", "120.0", "GeV/c^{2}"),
                                                              probe_sc_et = cms.vstring("Probe E_{T}", "0", "1000", "GeV/c"),
                                                              probe_sc_abseta = cms.vstring("Probe #eta", "0", "2.5", ""),                
                                                              tag_sc_et = cms.vstring("Tag et", "0", "10000", ""),                
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
            #"RooCBExGaussShape::signalResPass(mass,meanP[0,-5.000,5.000],sigmaP[1.419,1.00,15.000],alphaP[0.142,0, 20],nP[50.08],sigmaP_2[1, 1, 15.000])",
            #"RooCBExGaussShape::signalResFail(mass,meanF[1.1,-5.000,5.000],sigmaF[1.331,0.100,15.000],alphaF[0.210, 0, 20], nF[4.242],sigmaF_2[1, 1, 15.000])",
            #"RooCBExGaussShape::signalResPass(mass,meanP[0,-5.000,5.000],sigmaP[1.419,1.00,15.000],alphaP[21.558],nP[37.08],sigmaP_2[1, 1, 15.000])",
            #"RooCBExGaussShape::signalResFail(mass,meanF[1.1,-5.000,5.000],sigmaF[1.331,0.100,15.000],alphaF[0.915], nF[7.178],sigmaF_2[1, 1, 15.000])",
            
            "RooCBExGaussShape::signalResPass(mass,meanP[0,-5.000,5.000],sigmaP[1.419,1.00,15.000],alphaP[0.130],nP[50.08],sigmaP_2[1, 1, 15.000])",
            "RooCBExGaussShape::signalResFail(mass,meanF[1.1,-5.000,5.000],sigmaF[1.331,0.100,15.000],alphaF[0.284], nF[3.350],sigmaF_2[1, 1, 15.000])",
            "ZGeneratorLineShape::signalPhy(mass)", ### NLO line shape
            #"RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,120.], betaPass[0.001, 0.,1.0], gammaPass[0.001, 0.,1.0], peakPass[90.0])",
            #"RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,120.], betaFail[0.001, 0.,1.0], gammaFail[0.001, 0.,1.0], peakFail[90.0])",
            "RooExponential::backgroundPass(mass, aPass[-0.1, -1., 0])", 
            "RooExponential::backgroundFail(mass, aFail[-0.1, -1., 0])", 
            "FCONV::signalPass(mass, signalPhy, signalResPass)",
            "FCONV::signalFail(mass, signalPhy, signalResFail)",     
            "efficiency[0.5,0,1]",
            "signalFractionInPassing[0.9]"     
            ),
                                                         ),

                                         # defines a set of efficiency calculations, what PDF to use for fitting and how to bin the data;
                                         # there will be a separate output directory for each calculation that includes a simultaneous fit, side band subtraction and counting. 
                                         Efficiencies = cms.PSet(
                                                                 #the name of the parameter set becomes the name of the directory
                                                                 Medium = cms.PSet(EfficiencyBinningSpecification,
                                                                                  EfficiencyCategoryAndState = cms.vstring("passingMedium", "pass"),
                                                                                   ),
                                                                 )
                                         )

process.fit = cms.Path(
    process.GsfElectronToId  
    )
