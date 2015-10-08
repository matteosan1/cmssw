import ROOT
from optparse import OptionParser

def main(options):
    
    listToConst = [options.params.split(",")]

    f = ROOT.TFile(options.input)
    directory = f.Get(options.directory)
    directory.cd()
    keys = directory.GetListOfKeys()
    for k in keys:
        if ("pdfSignalPlusBackground" in k.GetName()):
            print k.GetName()
            directory = f.Get(options.directory+"/"+k.GetName())
            directory.cd()
            results = directory.Get("fitresults")
            params = results.floatParsFinal()

            for p in xrange(params.getSize()):
                myPar = params.at(p)
                if (myPar.GetName() in listToConst): 
                    print "%s[%.3f]"%(myPar.GetName(), myPar.getVal())
                else:
                    print "%s[%.3f,%.3f,%.3f]"%(myPar.GetName(), myPar.getVal(), myPar.getMin(), myPar.getMax())
            print
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
    parser.add_option("-d", "--directory", default="GsfElectronToRECO/Medium", help="Directory with workspace")
    parser.add_option("-p", "--params", default="", help="List of parameters to set to constant")
    parser.add_option("-s", "--signal", default="", help="List of parameters to set to constant")

    (options, arg) = parser.parse_args()
     
    main(options)
