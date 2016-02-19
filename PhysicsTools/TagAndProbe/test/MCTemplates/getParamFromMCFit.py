import ROOT
from optparse import OptionParser

def main(options):
    
    var1s = []
    for v in options.var1Bins.split(","):
        var1s.append(float(v))
    var2s = []
    for v in options.var2Bins.split(","):
        var2s.append(float(v))
    
    listToConst = [p for p in options.params.split(",")]
    #print listToConst
    f = ROOT.TFile(options.input)

    pdfDict = {}
    directory = f.Get(options.directory)
    directory.cd()
    keys = directory.GetListOfKeys()
    for k in keys:
        if ("pdfSignalPlusBackground" in k.GetName()):
            directory = f.Get(options.directory+"/"+k.GetName())
            directory.cd()
            results = directory.Get("fitresults")
            params = results.floatParsFinal()

            pdfDict[k.GetName()] = {}
            varString = ""
            for p in xrange(params.getSize()):
                myPar = params.at(p)
                if (myPar.GetName() in listToConst): 
                    varString = "%s[%.3f]"%(myPar.GetName(), myPar.getVal())
                else:
                    varString = "%s[%.3f,%.3f,%.3f]"%(myPar.GetName(), myPar.getVal(), myPar.getMin(), myPar.getMax())
                
                pdfDict[k.GetName()][myPar.GetName()] = varString

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

            binName = options.var1Name +"_bin"+str(binVar1)+"__"+options.var2Name+"_bin"+str(binVar2)+"__mcTrue_true__pdfSignalPlusBackground"

            line = "\"RooCBExGaussShape::signalResPass(mass," + pdfDict[binName]['meanP']+","+ pdfDict[binName]['sigmaP']+","+ pdfDict[binName]['alphaP']+","+ pdfDict[binName]['nP']+","+ pdfDict[binName]['sigmaP_2']+")\",\n"
           
            outputFile.write(line)
            line = "\"RooCBExGaussShape::signalResFail(mass," + pdfDict[binName]['meanF']+","+ pdfDict[binName]['sigmaF']+","+ pdfDict[binName]['alphaF']+","+ pdfDict[binName]['nF']+","+ pdfDict[binName]['sigmaF_2']+")\",\n"
            outputFile.write(line)
            histNameSt = "hMass_"+str(var1s[binVar1])+"To"+str(var1s[binVar1+1])+"_"+str(var2s[binVar2])+"To"+str(var2s[binVar2+1])
            outputFile.write("\"ZGeneratorLineShape::signalPhyPass(mass)\",\n"),
            outputFile.write("\"ZGeneratorLineShape::signalPhyFail(mass)\",\n"),
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
            #w = directory.Get("w")
            #print "GsfElectronToRECO/Medium/"+k.GetName()
            #w.Print()
            #for s in options.signal.split(","):
            #    print s
            #    continue
            #    pdf = w.pdf("signal")
            #
            #    argset = pdf.getParameters(ROOT.RooArgSet(w.var("mass")))
            #    it = argset.createIterator()
            #    var = it.Next()
            #
            #    funcString = pdf.ClassName()+"::"+pdf.GetName()+"(mass"
            #    while (var):
            #        if (var.GetName() in listToConst):
            #            funcString += ", %s[%.3f]"%(var.GetName(), var.getVal())
            #        else:
            #            funcString += ", %s[%.3f,%.3f,%.3f]"%(var.GetName(), var.getVal(), var.getMax(), var.getMin())
            #
            #        var = it.Next()
            #        funcString += ")"
            #    print funcString
        

if __name__ == "__main__":  
    parser = OptionParser()
    parser.add_option("-i", "--input", default="test.root", help="Input filename")
    parser.add_option("-o", "--outputFile", default="../python/parametricTemplates.py", help="Output filename")
    parser.add_option("-t", "--templateFile", default="templatesID.root", help="Output filename")
    parser.add_option("-d", "--directory", default="GsfElectronToRECO/Medium", help="Directory with workspace")
    parser.add_option("-p", "--params", default="", help="List of parameters to set to constant")
    parser.add_option("", "--var1Bins", default="20,30,40,50,200", help="Binning to use in var1")
    parser.add_option("", "--var2Bins", default="0.0,1.0,1.4442,1.566,2.0,2.5", help="Binning to use in var2")
    parser.add_option("", "--var1Name", default="probe_sc_eta", help="Variable1 branch name")
    parser.add_option("", "--var2Name", default="probe_sc_et", help="Variable2 branch name")
    parser.add_option("", "--failBkgPdf", default="RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,70.], betaFail[0.001, 0.,0.1], gammaFail[0.1, 0, 1], peakFail[90.0])", help="Background PDF for passing probes")
    parser.add_option("", "--passBkgPdf", default="RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,70.], betaPass[0.001, 0.,0.1], gammaPass[0.1, 0, 1], peakPass[90.0])", help="Background PDF for failing probes")
    parser.add_option("", "--idLabel", default="passingPresel", help="String identifying ID WP to measure")

    (options, arg) = parser.parse_args()
     
    main(options)
