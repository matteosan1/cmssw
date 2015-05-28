import ROOT

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


def main():
    inputFileName = "efficiency-mc-GsfElectronToId.root"
    dir="GsfElectronToRECO"
    isCutAndCount = False 
    isMCTruth = False   
    
    print "   ##################################################   "

    if(isCutAndCount):
        print "Plotting efficiency from cut & count. No background subtraction performed !"
        print "If you want to plot MC truth efficiency, please set: isMCTruth = true."
    else:
        print "Plotting efficiency from simultaneous fit."
        print "   ##################################################   "
  
        
    f = ROOT.TFile(inputFileName)
    f.cd(dir)

    keyList = [key.GetName() for key in ROOT.gDirectory.GetListOfKeys()]
    for k in  keyList:
        obj = ROOT.gDirectory.GetKey(k).ReadObj();
        name = obj.GetName()
           
        if (not obj.IsA().InheritsFrom("TDirectory")):
            continue
        if (not isMCTruth and "MCtruth" in name):
            continue
        if (isMCTruth and not "MCtruth" in name):
            continue

        print  "   ==================================================   "
        dirName = "%s/cnt_eff_plots/"%(name)
        if(not isCutAndCount):
            dirName = "%s/fit_eff_plots/"%(name)

        print "****************************dirName = ", dirName
        ROOT.gDirectory.cd(dirName)
        keyList2 = [key.GetName() for key in ROOT.gDirectory.GetListOfKeys()]
        for k in  keyList2:
            obj = ROOT.gDirectory.GetKey(k).ReadObj();
            innername = obj.GetName()
            if (obj.ClassName() == "TCanvas"):
                for p in obj.GetListOfPrimitives():
                    if (p.ClassName() == "TH2F"):
                        makeTable(p, name+".txt")
                        

    if(not isCutAndCount):
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

    if(not isCutAndCount): 
        ROOT.gDirectory.cd("../")
    else:
        ROOT.gDirectory.cd("../../")

    print "   ==================================================   "


if (__name__ == "__main__"):
    main()


