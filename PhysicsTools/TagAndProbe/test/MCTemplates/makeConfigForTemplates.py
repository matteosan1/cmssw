#FIXME ADD TEMPLATE FOR BG PDF CHOICE
import ROOT
from optparse import OptionParser
from getTemplatesFromMC import findBins

def main(options):
    pts = []
    for v in options.ptbins.split(","):
        pts.append(float(v))
    etas = []
    for v in options.etabins.split(","):
        etas.append(float(v))
    
    outputFile = file(options.outputFile, "w")
    
    outputFile.write("import FWCore.ParameterSet.Config as cms\n\n")
    outputFile.write("all_pdfs = cms.PSet(\n")

    for binPt in xrange(len(pts)-1):
        for binEta in xrange(len(etas)-1):
            psetName = options.idLabel+"_"+str(pts[binPt])+"To"+str(pts[binPt+1])+"_"+str(etas[binEta])+"To"+str(etas[binEta+1])
            psetName = psetName.replace(".", "p") + " = cms.vstring(\n"
            outputFile.write(psetName)
            #outputFile.write(options.idLabel+"_ptBin"+str(binPt)+"_etaBin"+str(binEta)+" = cms.vstring(\n")
            outputFile.write("\"RooGaussian::signalResPass(mass, meanP[.0,-5.000,5.000],sigmaP[0.956,0.00,5.000])\",\n")
            outputFile.write("\"RooGaussian::signalResFail(mass, meanF[.0,-5.000,5.000],sigmaF[0.956,0.00,5.000])\",\n")
            histNameSt = "hMass_"+str(pts[binPt])+"To"+str(pts[binPt+1])+"_"+str(etas[binEta])+"To"+str(etas[binEta+1])
            outputFile.write("\"ZGeneratorLineShape::signalPhyPass(mass,\\\""+options.templateFile+"\\\", \\\""+histNameSt+"_Pass\\\")\",\n"),
            outputFile.write("\"ZGeneratorLineShape::signalPhyFail(mass,\\\""+options.templateFile+"\\\", \\\""+histNameSt+"_Fail\\\")\",\n"),
            outputFile.write("\"RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,70.], betaPass[0.001, 0.,0.1], gammaPass[0.1, 0, 1], peakPass[90.0])\",\n")
            outputFile.write("\"RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,70.], betaFail[0.001, 0.,0.1], gammaFail[0.1, 0, 1], peakFail[90.0])\",\n")
            outputFile.write("\"FCONV::signalPass(mass, signalPhyPass, signalResPass)\",\n")
            outputFile.write("\"FCONV::signalFail(mass, signalPhyFail, signalResFail)\",\n")     
            outputFile.write("\"efficiency[0.5,0,1]\",\n")
            outputFile.write("\"signalFractionInPassing[1.]\"\n")     
            outputFile.write("),\n")
            outputFile.write("\n")
    outputFile.write(")")
    outputFile.close()



if __name__ == "__main__":  
    parser = OptionParser()

    parser.add_option("-i", "--idLabel", default="pdf", help="Prefix for block of PDFs")
    parser.add_option("-o", "--outputFile", default="commonFit.py", help="Output filename")
    parser.add_option("-t", "--templateFile", default="templates.root", help="Output filename")
    parser.add_option("", "--ptbins", default="20,30,40,50,200", help="Binning to use in pT")
    parser.add_option("", "--etabins", default="0.0,1.0,1.4442,1.566,2.0,2.5", help="Binning to use in eta")

    (options, arg) = parser.parse_args()
     
    main(options)
