import ROOT
from optparse import OptionParser
import sys

def makeTable(h, tablefilename):
    nX = h.GetNbinsX()
    nY = h.GetNbinsY()
  
    print "Writing...", tablefilename
    f = open(tablefilename, "w+")

    for i in xrange(1, nX+1):
    
        pT0 = h.GetXaxis().GetBinLowEdge(i)
        pT1 = h.GetXaxis().GetBinLowEdge(i+1)
    
        for j in xrange(1, nY+1):
            x    = h.GetBinContent(i,j)
            dx   = h.GetBinError(i,j)
            eta0 = h.GetYaxis().GetBinLowEdge(j)
            eta1 = h.GetYaxis().GetBinLowEdge(j+1)
      
            f.write("%4.2f  %4.2f   %+6.4f   %+6.4f  %6.4f   %6.4f \n" %(pT0, pT1, eta0, eta1, x, dx))
  
    f.close()


def main(options):
    ROOT.gStyle.SetPaintTextFormat("1.4f")
    print "##################################################   "
    f = ROOT.TFile(options.input)
    f.cd(options.directory+"/"+options.name)
    
    keyList = [key.GetName() for key in ROOT.gDirectory.GetListOfKeys()]

    # EFFICIENCY PLOTS + TABLE
    dirName = "cnt_eff_plots"
    if(not options.cc):
        dirName = "fit_eff_plots"
    ROOT.gDirectory.cd(dirName)

    keyList2 = [key.GetName() for key in ROOT.gDirectory.GetListOfKeys()]
    tableDone = False

    for k in  keyList2:
        obj = ROOT.gDirectory.GetKey(k).ReadObj();
        innername = obj.GetName()
        if (obj.ClassName() == "TCanvas"):
            for p in obj.GetListOfPrimitives():
                if (p.ClassName() == "TH2F" and not tableDone):
                    makeTable(p, options.name+".txt")
                    tableDone = True
            obj.Draw()
            innername = innername.replace("&", "")            
            plotname = innername + ".png"
            obj.SaveAs(plotname)

    if (options.dumpPlots and not options.cc):
        ROOT.gDirectory.cd("../")
        keyList = [key.GetName() for key in ROOT.gDirectory.GetListOfKeys()]
        for k in  keyList:
            obj = ROOT.gDirectory.GetKey(k).ReadObj();
            innername = obj.GetName()
            if (not obj.IsA().InheritsFrom("TDirectory") or not "_bin" in innername):
                continue
            ROOT.gDirectory.cd(innername)
            c = ROOT.gDirectory.Get("fit_canvas")
            c.Draw()
            plotname = "fit" + options.name + "_" + innername + ".png"
            #plotname = plotname.replace("probe_sc_", "")
            plotname = plotname.replace("&", "")
            plotname = plotname.replace("__pdfSignalPlusBackground", "")
            c.SaveAs(plotname)
            ROOT.gDirectory.cd("../")

if (__name__ == "__main__"):
    parser = OptionParser()
    parser.add_option("-i", "--input", default="efficiency-mc-GsfElectronToId.root", help="Input filename")
    parser.add_option("-d", "--directory", default="GsfElectronToRECO", help="Directory with workspace")
    parser.add_option("-c", dest="cc", action="store_true", help="Is simple Cut and Count", default=False)
    parser.add_option("-b", dest="batch", action="store_true", help="ROOT batch mode", default=False)
    parser.add_option("-m", "--name", help="Subdirectory with results", default="xxx")
    parser.add_option("-p", "--dump-plot", dest="dumpPlots", action="store_true", help="Dump fit plots", default=False)

    (options, arg) = parser.parse_args()

    if (options.batch):
        ROOT.gROOT.SetBatch(True)

    main(options)


