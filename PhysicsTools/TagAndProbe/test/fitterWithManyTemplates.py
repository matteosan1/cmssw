import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import PhysicsTools.TagAndProbe.parametricTemplatesWP80MC as common

options = VarParsing('analysis')
options.register(
    "isMC",
    True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Compute efficiency for MC"
    )

options.register(
    "inputFileName",
    "/afs/cern.ch/work/i/ishvetso/public/for_Matteo/TnPTree_mc-powheg.root",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Input filename"
    )

options.register(
    "outputFileName",
    "",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output filename"
    )

options.register(
    "idName",
    "passingTight",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "ID variable name as in the fitter_tree"
    )

options.register(
    "dirName",
    "GsfElectronToRECO",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Folder name containing the fitter_tree"
    )

options.register(
    "doCutAndCount",
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Do not compute fitting, just cut and count"
    )
options.parseArguments()

process = cms.Process("TagProbe")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

################################################

InputFileName = options.inputFileName
OutputFile = "efficiency-mc-"+options.idName
if (not options.isMC):
    OutputFile = "efficiency-data-"+options.idName

if (options.outputFileName != ""):
    OutputFile = OutputFile+"-"+options.outputFileName+".root"
else:
    OutputFile = OutputFile+".root"

################################################

#specifies the binning of parameters
EfficiencyBins = cms.PSet(
    probe_Ele_et = cms.vdouble(20. ,40. ,60. ,100.),
    probe_sc_eta = cms.vdouble(-2.5, -1.0, 0.0, 1.0, 2.5),
    )

DataBinningSpecification = cms.PSet(
    UnbinnedVariables = cms.vstring("mass"),
    BinnedVariables = cms.PSet(EfficiencyBins),
    BinToPDFmap = cms.vstring(
        "tight_20p0To40p0_0p0To1p5", 
        "*et_bin0*eta_bin0*","tight_20p0To40p0_m2p5Tom1p0",
        "*et_bin0*eta_bin1*","tight_20p0To40p0_m1p0To0p0",
        "*et_bin0*eta_bin2*","tight_20p0To40p0_0p0To1p0",
        "*et_bin0*eta_bin3*","tight_20p0To40p0_1p0To2p5",

        "*et_bin1*eta_bin0*","tight_40p0To60p0_m2p5Tom1p0",
        "*et_bin1*eta_bin1*","tight_40p0To60p0_m1p0To0p0",
        "*et_bin1*eta_bin2*","tight_40p0To60p0_0p0To1p0",
        "*et_bin1*eta_bin3*","tight_40p0To60p0_1p0To2p5",

        "*et_bin2*eta_bin0*","tight_60p0To100p0_m2p5Tom1p0",
        "*et_bin2*eta_bin1*","tight_60p0To100p0_m1p0To0p0",
        "*et_bin2*eta_bin2*","tight_60p0To100p0_0p0To1p0",
        "*et_bin2*eta_bin3*","tight_60p0To100p0_1p0To2p5",
        )
    )

McBinningSpecification = cms.PSet(
    UnbinnedVariables = cms.vstring("mass", "totWeight"),
    BinnedVariables = cms.PSet(EfficiencyBins, mcTrue = cms.vstring("true")),
    BinToPDFmap = cms.vstring(
        "tight_20p0To40p0_0p0To1p5", 
        "*et_bin0*eta_bin0*","tight_20p0To40p0_0p0To1p5",
        "*et_bin1*eta_bin0*","tight_40p0To60p0_0p0To1p5",
        "*et_bin2*eta_bin0*","tight_60p0To100p0_0p0To1p5",
        "*et_bin0*eta_bin1*","tight_20p0To40p0_1p5To2p5",
        "*et_bin1*eta_bin1*","tight_40p0To60p0_1p5To2p5",
        "*et_bin2*eta_bin1*","tight_60p0To100p0_1p5To2p5",
        )
)

########################

process.TnPMeasurement = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
                                        InputFileNames = cms.vstring(InputFileName),
                                        InputDirectoryName = cms.string(options.dirName),
                                        InputTreeName = cms.string("fitter_tree"), 
                                        OutputFileName = cms.string(OutputFile),
                                        NumCPU = cms.uint32(1),
                                        SaveWorkspace = cms.bool(False), #VERY TIME CONSUMING FOR MC
                                        doCutAndCount = cms.bool(options.doCutAndCount),
                                        floatShapeParameters = cms.bool(True),
                                        binnedFit = cms.bool(True),
                                        binsForFit = cms.uint32(60),
                                        WeightVariable = cms.string("totWeight"),
                                        # defines all the real variables of the probes available in the input tree and intended for use in the efficiencies
                                        Variables = cms.PSet(
        mass = cms.vstring("Tag-Probe Mass", "60.0", "120.0", "GeV/c^{2}"),
        probe_Ele_et = cms.vstring("Probe E_{T}", "0", "100", "GeV/c"),
        probe_sc_eta = cms.vstring("Probe #eta", "-2.5", "2.5", ""), 
        totWeight = cms.vstring("totWeight", "-1000000", "100000000", ""),
        ),
                                        
                                        # defines all the discrete variables of the probes available in the input tree and intended for use in the efficiency calculation
                                        Expressions = cms.PSet(),
                                        Categories = cms.PSet(),
                                        PDFs = common.all_pdfs,
                                        Efficiencies = cms.PSet()
                                        )

setattr(process.TnPMeasurement.Categories, options.idName, cms.vstring(options.idName, "dummy[pass=1,fail=0]"))
setattr(process.TnPMeasurement.Categories, "mcTrue", cms.vstring("MC true", "dummy[true=1,false=0]"))

if (not options.isMC):
    delattr(process.TnPMeasurement, "WeightVariable")
    process.TnPMeasurement.Variables = cms.PSet(
        mass = cms.vstring("Tag-Probe Mass", "60.0", "120.0", "GeV/c^{2}"),
        probe_Ele_et = cms.vstring("Probe E_{T}", "20", "100", "GeV/c"),
        probe_sc_eta = cms.vstring("Probe #eta", "-2.5", "2.5", ""), 
        )
    for pdf in process.TnPMeasurement.PDFs.__dict__:
        param =  process.TnPMeasurement.PDFs.getParameter(pdf)
        if (type(param) is not cms.vstring):
            continue
        for i, l in enumerate(getattr(process.TnPMeasurement.PDFs, pdf)):
            if l.find("signalFractionInPassing") != -1:
                getattr(process.TnPMeasurement.PDFs, pdf)[i] = l.replace("[1.0]","[0.5,0.,1.]")

    setattr(process.TnPMeasurement.Efficiencies, options.idName, DataBinningSpecification)    
    setattr(getattr(process.TnPMeasurement.Efficiencies, options.idName) , "EfficiencyCategoryAndState", cms.vstring(options.idName, "pass"))
else:   
    setattr(process.TnPMeasurement.Efficiencies, "MCtruth_" + options.idName, McBinningSpecification)    
    setattr(getattr(process.TnPMeasurement.Efficiencies, "MCtruth_" + options.idName), "EfficiencyCategoryAndState", cms.vstring(options.idName, "pass"))

    for pdf in process.TnPMeasurement.PDFs.__dict__:
        param =  process.TnPMeasurement.PDFs.getParameter(pdf)
        if (type(param) is not cms.vstring):
            continue
        for i, l in enumerate(getattr(process.TnPMeasurement.PDFs, pdf)):
            if l.find("backgroundPass") != -1:
                getattr(process.TnPMeasurement.PDFs, pdf)[i] = "RooPolynomial::backgroundPass(mass, a[0.0])"
            if l.find("backgroundFail") != -1:
                getattr(process.TnPMeasurement.PDFs, pdf)[i] = "RooPolynomial::backgroundFail(mass, a[0.0])"

process.fit = cms.Path(
    process.TnPMeasurement  
    )


