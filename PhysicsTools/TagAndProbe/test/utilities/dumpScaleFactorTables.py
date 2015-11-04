import ROOT
import math
from optparse import OptionParser

def makeTable(hnum, hden, tablefilename):
    nX = hnum.GetNbinsX()
    nY = hnum.GetNbinsY()
  
    f = open(tablefilename, "w+")
  
    for i in xrange(1, nX+1):
        pT0 = hnum.GetXaxis().GetBinLowEdge(i)
        pT1 = hnum.GetXaxis().GetBinLowEdge(i+1)
    
        for j in xrange(1, nY+1):
            x = hnum.GetBinContent(i,j)/hden.GetBinContent(i,j)
            dx1 = hnum.GetBinError(i,j)/hnum.GetBinContent(i,j)
            dx2 = hden.GetBinError(i,j)/hden.GetBinContent(i,j)
            dx = math.sqrt(dx1*dx1+dx2*dx2)*x
            print hnum.GetBinError(i,j), hnum.GetBinContent(i,j)
            print dx1, dx2, dx
            eta0 = hnum.GetYaxis().GetBinLowEdge(j)
            eta1 = hnum.GetYaxis().GetBinLowEdge(j+1)
            
            f.write("%4.1f  %4.1f   %+6.4f   %+6.4f  %6.4f   %6.4f \n"%(pT0, pT1, eta0, eta1, x, dx))
            
    f.close()

def main(options):
    fData = ROOT.TFile(options.data)
    fMC   = ROOT.TFile(options.mc)
    hData = ""
    hMC = ""

    temp = "%s/%s/fit_eff_plots/" % (options.directory, options.name)
    #if("ToHLT" in temp):
    #    temp = "%s/%s/cnt_eff_plots/" %(options.directory, options.name)
    fData.cd(temp)

    keyList = [key.GetName() for key in ROOT.gDirectory.GetListOfKeys()]
    for k in  keyList:
        obj = ROOT.gDirectory.GetKey(k).ReadObj();
        innername = obj.GetName()
        if (obj.ClassName() == "TCanvas"):
            for p in obj.GetListOfPrimitives():
                if (p.ClassName() == "TH2F"):
                    hData = p

    temp = "%s/MCtruth_%s/fit_eff_plots/" % (options.directory, options.name)
    fMC.cd(temp)
    keyList = [key.GetName() for key in ROOT.gDirectory.GetListOfKeys()]
    for k in  keyList:
        obj = ROOT.gDirectory.GetKey(k).ReadObj();
        innername = obj.GetName()
        if (obj.ClassName() == "TCanvas"):
            for p in obj.GetListOfPrimitives():
                if (p.ClassName() == "TH2F"):
                    hMC = p

    temp = "ScaleFactor_%s_%s.txt"%(options.directory, options.name)
    #hData.Divide(hMC)
    makeTable(hData, hMC, temp)
    
    fData.Close()
    fMC.Close()

if (__name__ == "__main__"):
    parser = OptionParser()
    parser.add_option("", "--mc", default="efficiency-mc-GsfElectronToId.root", help="Input filename for MC")
    parser.add_option("", "--data", default="efficiency-data-GsfElectronToId.root", help="Input filename for data")
    parser.add_option("-d", "--directory", default="GsfElectronToRECO", help="Directory with workspace")
    parser.add_option("-n", "--name", default="Medium", help="Subdirectory with results")
    parser.add_option("-b", dest="batch", action="store_true", help="ROOT batch mode", default=False)
    
    (options, arg) = parser.parse_args()

    if (options.batch):
        ROOT.gROOT.SetBatch(True)

    main(options)
