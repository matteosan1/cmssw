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
      
            f.write("%4.1f  %4.1f   %+6.4f   %+6.4f  %6.4f   %6.4f \n" %(pT0, pT1, eta0, eta1, x, dx))
  
    f.close()


def main(options):

    print "##################################################   "
    if(options.cc):
        print "Plotting efficiency from cut & count. No background subtraction performed !"
        print "If you want to plot MC truth efficiency, please set: isMCTruth = true."
    else:
        print "Plotting efficiency from simultaneous fit."
        print "##################################################   "
  
        
    f = ROOT.TFile(options.input)
    f.cd(options.directory)

    keyList = [key.GetName() for key in ROOT.gDirectory.GetListOfKeys()]
    for k in  keyList:
        obj = ROOT.gDirectory.GetKey(k).ReadObj();
        name = obj.GetName()

        if (not obj.IsA().InheritsFrom("TDirectory")):
            continue
        if (not options.isMCTruth and "MCtruth" in name):
            continue
        if (options.isMCTruth and not "MCtruth" in name):
            continue

        print  "   ==================================================   "
        dirName = "%s/cnt_eff_plots/"%(name)
        if(not options.cc):
            dirName = "%s/fit_eff_plots/"%(name)

        print "****************************dirName = ", dirName
        ROOT.gDirectory.cd(dirName)
        keyList2 = [key.GetName() for key in ROOT.gDirectory.GetListOfKeys()]
        tableDone = False
        for k in  keyList2:
            obj = ROOT.gDirectory.GetKey(k).ReadObj();
            innername = obj.GetName()
            if (obj.ClassName() == "TCanvas"):
                for p in obj.GetListOfPrimitives():
                    if (p.ClassName() == "TH2F" and not tableDone):
                        makeTable(p, name+".txt")
                        tableDone = True
                obj.Draw()
                plotname = innername + ".png"
                obj.SaveAs(plotname)

    if(not options.cc):
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
            plotname = "fit" + name + "_" + innername + ".png"
            plotname.replace("probe_sc_", "")
            plotname.replace("__pdfSignalPlusBackground", "")
            c.SaveAs(plotname)
            ROOT.gDirectory.cd("../")

    if(not options.cc): 
        ROOT.gDirectory.cd("../")
    else:
        ROOT.gDirectory.cd("../../")

    print "   ==================================================   "


if (__name__ == "__main__"):
    parser = OptionParser()
    parser.add_option("-i", "--input", default="efficiency-mc-GsfElectronToId.root", help="Input filename")
    parser.add_option("-d", "--directory", default="GsfElectronToRECO", help="Directory with workspace")
    parser.add_option("-c", dest="cc", action="store_true", help="Is simple Cut and Count", default=False)
    parser.add_option("-b", dest="batch", action="store_true", help="ROOT batch mode", default=False)
    parser.add_option("-m", "--mc", dest="isMCTruth", action="store_true", default="Medium", help="Subdirectory with results", default=False)

    (options, arg) = parser.parse_args()

    if (options.batch):
        ROOT.gROOT.SetBatch(True)

    main(options)


