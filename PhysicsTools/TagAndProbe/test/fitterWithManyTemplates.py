import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import PhysicsTools.TagAndProbe.commonFit as common

options = VarParsing('analysis')

options.register(
    "isData",
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Compute data efficiencies"
    )

options.register(
    "isMC",
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Compute MC efficiencies"
    )

options.register(
    "inputFileName",
    "TnP_Run2015D_30Sep.root",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Don't compute MC efficiencies"
    )

options.register(
    "outputFileName",
    "TnP_Run2015D_30Sep.root",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Don't compute MC efficiencies"
    )

options.parseArguments()

if (not options.isMC) and (not options.isData):
    raise Exception("Must select either data or MC")

if (options.isData):
    for pdf in common.all_pdfs.__dict__:
        param = common.all_pdfs.getParameter(pdf)
        if type(param) is not cms.vstring:
            continue
        i = 0
        for line in getattr(common.all_pdfs, pdf):
            if line.find("signalFractionInPassing") != -1:
                getattr(common.all_pdfs, pdf)[i] = line.replace("[1.]","[0.5,0.,1.]")
            i = i + 1


process = cms.Process("TagProbe")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

InputFileName = options.inputFileName
OutputFile = options.outputFileName
PDFName = "pdfSignalPlusBackground"

################################################

#specifies the binning of parameters
EfficiencyBins = cms.PSet(
    probe_Pho_et = cms.vdouble(20. ,40. ,60. ,100.),
    probe_sc_abseta = cms.vdouble(0.0, 1.0, 1.5, 2.0, 2.5),
    )

DataBinningSpecification = cms.PSet(
    UnbinnedVariables = cms.vstring("mass"), #"totWeight"), #, "Ele_dRTau", "probe_dRTau"),
    BinnedVariables = cms.PSet(EfficiencyBins),
    BinToPDFmap = cms.vstring(
        "id_20.0To40.0_1.0To1.5", # FIXME ADDING A DEFAULT CONFIG
        "*et_bin0*abseta_bin0*","id_20p0To40p0_0p0To1p0",
        "*et_bin1*abseta_bin0*","id_40p0To60p0_0p0To1p0",
        "*et_bin2*abseta_bin0*","id_60p0To100p0_0p0To1p0",
        "*et_bin0*abseta_bin1*","id_20p0To40p0_1p0To1p5",
        "*et_bin1*abseta_bin1*","id_40p0To60p0_1p0To1p5",
        "*et_bin2*abseta_bin1*","id_60p0To100p0_1p0To1p5",
        "*et_bin0*abseta_bin2*","id_20p0To40p0_1p5To2p0",
        "*et_bin1*abseta_bin2*","id_40p0To60p0_1p5To2p0",
        "*et_bin2*abseta_bin2*","id_60p0To100p0_1p5To2p0",
        "*et_bin0*abseta_bin3*","id_20p0To40p0_2p0To2p5",
        "*et_bin1*abseta_bin3*","id_40p0To60p0_2p0To2p5",
        "*et_bin2*abseta_bin3*","id_60p0To100p0_2p0To2p5",
        )
    )

McBinningSpecification = cms.PSet(
    UnbinnedVariables = cms.vstring("mass", "totWeight"), #"Ele_dRTau", "probe_dRTau"),
    BinnedVariables = cms.PSet(EfficiencyBins, mcTrue = cms.vstring("true")),
    BinToPDFmap = cms.vstring(
        "id_20.0To40.0_1.0To1.5", # FIXME ADDING A DEFAULT CONFIG
        "*et_bin0*abseta_bin0*","id_20p0To40p0_0p0To1p0",
        "*et_bin1*abseta_bin0*","id_40p0To60p0_0p0To1p0",
        "*et_bin2*abseta_bin0*","id_60p0To100p0_0p0To1p0",
        "*et_bin0*abseta_bin1*","id_20p0To40p0_1p0To1p5",
        "*et_bin1*abseta_bin1*","id_40p0To60p0_1p0To1p5",
        "*et_bin2*abseta_bin1*","id_60p0To100p0_1p0To1p5",
        "*et_bin0*abseta_bin2*","id_20p0To40p0_1p5To2p0",
        "*et_bin1*abseta_bin2*","id_40p0To60p0_1p5To2p0",
        "*et_bin2*abseta_bin2*","id_60p0To100p0_1p5To2p0",
        "*et_bin0*abseta_bin3*","id_20p0To40p0_2p0To2p5",
        "*et_bin1*abseta_bin3*","id_40p0To60p0_2p0To2p5",
        "*et_bin2*abseta_bin3*","id_60p0To100p0_2p0To2p5",
        )
)

########################
########################
####### MCEfficiency
########################
########################

process.McEfficiency = cms.EDAnalyzer(
    "TagProbeFitTreeAnalyzer",
    InputFileNames = cms.vstring(InputFileName),
    InputDirectoryName = cms.string("PhotonToRECO"),
    InputTreeName = cms.string("fitter_tree"), 
    OutputFileName = cms.string(OutputFile),
    NumCPU = cms.uint32(6),
    SaveWorkspace = cms.bool(False), #VERY TIME CONSUMING FOR MC
    doCutAndCount = cms.bool(False),
    floatShapeParameters = cms.bool(True),
    binnedFit = cms.bool(True),
    binsForFit = cms.uint32(60),
    WeightVariable = cms.string("totWeight"),
    # defines all the real variables of the probes available in the input tree and intended for use in the efficiencies
    Variables = cms.PSet(
        mass = cms.vstring("Tag-Probe Mass", "60.0", "120.0", "GeV/c^{2}"),
        probe_Pho_et = cms.vstring("Probe E_{T}", "20", "100", "GeV/c"),
        probe_sc_abseta = cms.vstring("Probe #eta", "0", "2.5", ""), 
        totWeight = cms.vstring("totWeight", "-1000000", "100000000", ""), 
        #Ele_dRTau = cms.vstring("Ele_dRTau", "0.2", "100000", ""),
        #probe_dRTau = cms.vstring("probe_dRTau", "0.2", "100000", ""),
        ),

    # defines all the discrete variables of the probes available in the input tree and intended for use in the efficiency calculation
    Categories = cms.PSet(
        mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
        passingID = cms.vstring("passingID", "dummy[pass=1,fail=0]"),
        ),
    PDFs = common.all_pdfs,
    # defines a set of efficiency calculations, what PDF to use for fitting and how to bin the data;
    # there will be a separate output directory for each calculation that includes a
    # simultaneous fit, side band subtraction and counting. 
    Efficiencies = cms.PSet(
        MCtruth_ID = cms.PSet(
            McBinningSpecification,
            EfficiencyCategoryAndState = cms.vstring("passingID", "pass"),
            ),
        )
    )


process.DataEfficiency = process.McEfficiency.clone()
process.DataEfficiency.InputFileNames = cms.vstring(InputFileName)
process.DataEfficiency.OutputFileName = cms.string(OutputFile)
process.DataEfficiency.doCutAndCount = cms.bool(False)
delattr(process.DataEfficiency, "WeightVariable")
process.DataEfficiency.Variables = cms.PSet(
    mass = cms.vstring("Tag-Probe Mass", "60.0", "120.0", "GeV/c^{2}"),
    probe_Pho_et = cms.vstring("Probe E_{T}", "20", "100", "GeV/c"),
    probe_sc_abseta = cms.vstring("Probe #eta", "0", "2.5", ""), 
    )
process.DataEfficiency.Categories = cms.PSet(passingID = cms.vstring("passingID", "dummy[pass=1,fail=0]"))
process.DataEfficiency.Efficiencies = cms.PSet(
    ID = cms.PSet(
        DataBinningSpecification,
        EfficiencyCategoryAndState = cms.vstring("passingID", "pass"),
        ),
    )

process.seq = cms.Sequence()

if (not options.noMC):
    process.seq += process.McEfficiency

if (not options.noData):
    process.seq += process.DataEfficiency

process.fit = cms.Path(process.seq)
