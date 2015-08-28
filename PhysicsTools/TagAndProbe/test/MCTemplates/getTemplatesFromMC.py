import ROOT
from optparse import OptionParser

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
    
    pts  = []
    for v in options.ptbins.split(","):
        pts.append(float(v))
    etas = []
    for v in options.etabins.split(","):
        etas.append(float(v))

    inFile = ROOT.TFile(options.input)
    inFile.cd(options.directory)
    fDir = inFile.Get(options.directory)
    fChain = fDir.Get("fitter_tree")

    histos = dict()

    for binPt in xrange(len(pts)-1):
        for binEta in xrange(len(etas)-1):
            print "Doing templates for "+str(pts[binPt])+"To"+str(pts[binPt+1])+"_"+str(etas[binEta])+"To"+str(etas[binEta+1])
            histNameSt = "hMass_"+str(pts[binPt])+"To"+str(pts[binPt+1])+"_"+str(etas[binEta])+"To"+str(etas[binEta+1])
            hp = histNameSt+"_Pass"
            hf = histNameSt+"_Fail"
            histos[hp] = ROOT.TH1D(hp, hp, 120, 60, 120)
            histos[hf] = ROOT.TH1D(hf, hf, 120, 60, 120)
            
            binning = options.tagTauVarName+" > 0.2 && "+options.probeTauVarName+" > 0.2 && mcTrue == 1 && pair_mass60to120 && "+options.etVarName +">"+str(pts[binPt])+" && "+options.etVarName +"<"+str(pts[binPt+1])+" && "+options.etaVarName +">"+str(etas[binEta])+" && "+options.etaVarName +"<"+str(etas[binEta+1])            
            cuts = "("+binning + " && passing"+options.idprobe+"==1"+")*"+options.weightVarName
            fChain.Draw("mass>>"+histos[hp].GetName(), cuts, "goff")
            cuts = binning + " && passing"+options.idprobe+"==0"
            fChain.Draw("mass>>"+histos[hf].GetName(), cuts, "goff")
    
    outFile = ROOT.TFile(options.output, "RECREATE")
    for k in histos:
        histos[k].Write()
    outFile.Close()


if __name__ == "__main__":  
    parser = OptionParser()
    parser.add_option("-i", "--input", default="../TnPTree_MC.root", help="Input filename")
    parser.add_option("-o", "--output", default="mc_templates.root", help="Output filename")
    parser.add_option("-d", "--directory", default="GsfElectronToRECO", help="Directory with fitter_tree")
    parser.add_option("", "--idprobe", default="Medium", help="String identifying ID WP to measure")
    parser.add_option("", "--ptbins", default="20,30,40,50,200", help="Binning to use in pT")
    parser.add_option("", "--etabins", default="0.0,1.0,1.4442,1.566,2.0,2.5", help="Binning to use in eta")
    parser.add_option("", "--etaVarName", default="probe_sc_eta", help="Eta variable branch name")
    parser.add_option("", "--etVarName", default="probe_sc_et", help="Et variable branch name")
    parser.add_option("", "--weightVarName", default="totWeight", help="Weight variable branch name")
    parser.add_option("", "--tagTauVarName", default="Ele_dRTau", help="Tag to tau dr variable branch name")
    parser.add_option("", "--probeTauVarName", default="probe_dRTau", help="Tag to tau dr variable branch name")

    (options, arg) = parser.parse_args()
     
    main(options)
