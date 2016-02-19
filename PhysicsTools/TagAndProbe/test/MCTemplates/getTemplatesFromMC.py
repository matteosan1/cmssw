import ROOT
from optparse import OptionParser

def removeNegativeBins(h):
    for i in xrange(h.GetNbinsX()):
        if (h.GetBinContent(i) < 0):
            h.SetBinContent(i, 0)
    
    

def findBins(array, var):
    size = len(array)
    bin = "dump";
    for i in xrange(size-1):
        low = array[i]
        hi  = array[i+1]
        if (low <= var and hi > var):
            bin = str(low)+"To"+str(hi)
  
    return bin;

def main(options):
    
    var1s  = []
    for v in options.var1Bins.split(","):
        var1s.append(float(v))
    var2s = []
    for v in options.var2Bins.split(","):
        var2s.append(float(v))

    inFile = ROOT.TFile(options.input)
    inFile.cd(options.directory)
    fDir = inFile.Get(options.directory)
    fChain = fDir.Get("fitter_tree")

    histos = dict()

    for binVar1 in xrange(len(var1s)-1):
        for binVar2 in xrange(len(var2s)-1):
            histNameSt = "hMass_"+str(var1s[binVar1])+"To"+str(var1s[binVar1+1])+"_"+str(var2s[binVar2])+"To"+str(var2s[binVar2+1])
            #histNameSt = histNameSt.replace("-", "m")
            print "Doing templates for "+histNameSt,
            hp = histNameSt+"_Pass"
            hf = histNameSt+"_Fail"
            histos[hp] = ROOT.TH1D(hp, hp, 120, 60, 120)
            histos[hf] = ROOT.TH1D(hf, hf, 120, 60, 120)
            
            #binning = options.tagTauVarName+" > 0.2 && "+options.probeTauVarName+" > 0.2 && mcTrue == 1 && pair_mass60to120 && "+options.etVarName +">"+str(pts[binVar1])+" && "+options.etVarName +"<"+str(pts[binVar1+1])+" && "+options.etaVarName +">"+str(etas[binVar2])+" && "+options.etaVarName +"<"+str(etas[binVar2+1])            
            binning = "mcTrue == 1 && pair_mass60to120 && "+options.var1Name +">"+str(var1s[binVar1])+" && "+options.var1Name +"<"+str(var1s[binVar1+1])+" && "+options.var2Name +">"+str(var2s[binVar2])+" && "+options.var2Name +"<"+str(var2s[binVar2+1])            

            cuts = "("
            if (options.addProbeCut != ""):
                cuts = cuts + options.addProbeCut + " && "
            cuts = cuts + binning + " && "+options.idprobe+"==1"+")*"+options.weightVarName
            fChain.Draw("mass>>"+histos[hp].GetName(), cuts, "goff")
            #print cuts
            
            cuts = "("
            if (options.addProbeCut != ""):
                cuts = cuts + options.addProbeCut + " && "
            cuts = cuts + binning + " && "+options.idprobe+"==0"+")*"+options.weightVarName
            fChain.Draw("mass>>"+histos[hf].GetName(), cuts, "goff")
            #print cuts

            removeNegativeBins(histos[hp])
            removeNegativeBins(histos[hf])

            hpassInt = histos[hp].Integral()
            hfailInt = histos[hf].Integral()
            
            print "%.2f %.2f %1.4f"%(hpassInt, hfailInt, hpassInt/(hpassInt+hfailInt))
    
    outFile = ROOT.TFile(options.output, "RECREATE")
    for k in histos:
        histos[k].Write()
    outFile.Close()


if __name__ == "__main__":  
    parser = OptionParser()
    parser.add_option("-i", "--input", default="../TnPTree_mc.root", help="Input filename")
    parser.add_option("-o", "--output", default="mc_templates.root", help="Output filename")
    parser.add_option("-d", "--directory", default="GsfElectronToRECO", help="Directory with fitter_tree")
    parser.add_option("", "--idprobe", default="passingMedium", help="String identifying ID WP to measure")
    parser.add_option("", "--var1Bins", default="20,30,40,50,200", help="Binning to use in var1")
    parser.add_option("", "--var2Bins", default="0.0,1.0,1.4442,1.566,2.0,2.5", help="Binning to use in var2")
    parser.add_option("", "--var1Name", default="probe_sc_eta", help="Variable1 branch name")
    parser.add_option("", "--var2Name", default="probe_sc_et", help="Variable2 branch name")
    parser.add_option("", "--addProbeCut", default="", help="Additional cut on the probe")
    parser.add_option("", "--weightVarName", default="totWeight", help="Weight variable branch name")
    #parser.add_option("", "--tagTauVarName", default="", help="Tag to tau dr variable branch name")
    #parser.add_option("", "--probeTauVarName", default="", help="Tag to tau dr variable branch name")

    (options, arg) = parser.parse_args()
     
    main(options)
