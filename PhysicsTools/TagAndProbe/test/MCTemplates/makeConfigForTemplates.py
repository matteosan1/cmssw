#FIXME ADD TEMPLATE FOR BG PDF CHOICE
import ROOT
from optparse import OptionParser
from getTemplatesFromMC import findBins

def main(options):
    var1s = []
    for v in options.var1Bins.split(","):
        var1s.append(float(v))
    var2s = []
    for v in options.var2Bins.split(","):
        var2s.append(float(v))
    
    outputFile = file(options.outputFile, "w")
    
    outputFile.write("import FWCore.ParameterSet.Config as cms\n\n")
    outputFile.write("all_pdfs = cms.PSet(\n")

    for binVar1 in xrange(len(var1s)-1):
        for binVar2 in xrange(len(var2s)-1):
            psetName = options.idLabel+"_"+str(var1s[binVar1])+"To"+str(var1s[binVar1+1])+"_"+str(var2s[binVar2])+"To"+str(var2s[binVar2+1])
            psetName = psetName.replace(".", "p")
            psetName = psetName.replace("-", "m")
            psetName = psetName + " = cms.vstring(\n"
            outputFile.write(psetName)

            outputFile.write("\"RooGaussian::signalResPass(mass, meanP[1.0,-5.000,5.000],sigmaP[0.5,0.001,5.000])\",\n")
            outputFile.write("\"RooGaussian::signalResFail(mass, meanF[1.0,-5.000,5.000],sigmaF[0.5,0.001,5.000])\",\n")
            histNameSt = "hMass_"+str(var1s[binVar1])+"To"+str(var1s[binVar1+1])+"_"+str(var2s[binVar2])+"To"+str(var2s[binVar2+1])
            outputFile.write("\"ZGeneratorLineShape::signalPhyPass(mass,\\\""+options.templateFile+"\\\", \\\""+histNameSt+"_Pass\\\")\",\n"),
            outputFile.write("\"ZGeneratorLineShape::signalPhyFail(mass,\\\""+options.templateFile+"\\\", \\\""+histNameSt+"_Fail\\\")\",\n"),
            outputFile.write("\""+options.passBkgPdf+"\",\n")
            outputFile.write("\""+options.failBkgPdf+"\",\n")
            outputFile.write("\"FCONV::signalPass(mass, signalPhyPass, signalResPass)\",\n")
            outputFile.write("\"FCONV::signalFail(mass, signalPhyFail, signalResFail)\",\n")     
            outputFile.write("\"efficiency[0.5,0,1]\",\n")
            outputFile.write("\"signalFractionInPassing[1.0]\"\n")     
            outputFile.write("),\n")
            outputFile.write("\n")
    outputFile.write(")")
    outputFile.close()

if __name__ == "__main__":  
    parser = OptionParser()

    parser.add_option("-i", "--idLabel", default="pdf", help="Prefix for block of PDFs")
    parser.add_option("-o", "--outputFile", default="../../python/commonFit.py", help="Output filename")
    parser.add_option("-t", "--templateFile", default="templatesID.root", help="Output filename")
    parser.add_option("", "--var1Bins", default="20,30,40,50,200", help="Binning to use in var1")
    parser.add_option("", "--var2Bins", default="0.0,1.0,1.4442,1.566,2.0,2.5", help="Binning to use in var2")
    parser.add_option("", "--failBkgPdf", default="RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,70.], betaFail[0.001, 0.,0.1], gammaFail[0.1, 0, 1], peakFail[90.0, 70, 80])", help="Background PDF for passing probes")
    parser.add_option("", "--passBkgPdf", default="RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,70.], betaPass[0.001, 0.,0.1], gammaPass[0.1, 0, 1], peakPass[90.0, 70, 80])", help="Background PDF for failing probes")

    (options, arg) = parser.parse_args()
     
    main(options)
