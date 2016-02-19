import ROOT
import array
from optparse import OptionParser
import sys

def updatePUWeights(hmc, hdata):
    hmc.Scale(1/hmc.Integral())
    hdata.Scale(1/hdata.Integral())

    hdata.Divide(hmc)
    
    print "You can take the following weights to update the PileupWeightProducer plugin."
    s = ""
    for i in xrange(60):
        print "%1.3f, "%(hdata.GetBinContent(i+1)), 
    print

    return hdata

def main(options):
    fmc = ROOT.TFile(options.mc)
    hmc = fmc.Get("sampleInfo/nvtx")
    mcSampleTree = fmc.Get("sampleInfo/sampleInfo")
    mcSampleTree.GetEntry(0)

    fdata = ROOT.TFile(options.data)
    hdata = fdata.Get("sampleInfo/nvtx")
    
    ratio = updatePUWeights(hmc, hdata)
    
    fmc.cd()
    directories = fmc.GetListOfKeys()

    outFile = options.mc.split(".root")[0]+"_norm.root"    
    newFile = ROOT.TFile(outFile, "RECREATE")

    for d in directories:
        if (d.GetName() == "sampleInfo"):
            continue
        print "Reducing tree in directory: ", d.GetName()
        directory = fmc.Get(d.GetName())

        tree = directory.Get("fitter_tree")
        entries = tree.GetEntries()  

        #--- Write to new file
        tree.SetBranchStatus("*", 1)
        tree.SetBranchStatus("totWeight", 0) 
        tree.SetBranchStatus("PUweight", 0) 

        directory_new = newFile.mkdir(d.GetName())
        directory_new.cd()
        tree_new = tree.CloneTree(0)

        tree.SetBranchStatus("totWeight", 1)
        tree.SetBranchStatus("PUweight", 1)

        # FIXME TO NORMALIZE TO LUMINOSITY
        normalization = mcSampleTree.nEvents/mcSampleTree.sumWeight
        if (options.noNegative):
            sumW = 0.
            newNEvents = 0.
            for z in range(entries):
                tree.GetEntry(z)
                if (tree.weight >=0):
                    sumW = sumW + tree.weight
                    newNEvents = newNEvents + 1.
    
            normalization = newNEvents/sumW
    
        b_totWeight = array.array('f', [0])
        b_PUweight = array.array('f', [0])
        tree_new.Branch("totWeight", b_totWeight, "totWeight/F")
        tree_new.Branch("PUweight", b_PUweight, "PUweight/F")

        for z in range(entries):
            tree.GetEntry(z)
            if (tree.weight >=0 or not options.noNegative):
                b_PUweight[0] = ratio.GetBinContent(tree.event_nPV)      
                b_totWeight[0] = tree.weight*b_PUweight[0]*normalization
                tree_new.Fill()
        tree_new.GetCurrentFile().Write()

    tree_new.GetCurrentFile().Close() 


if __name__ == "__main__":  
    parser = OptionParser()
    parser.add_option("-m", "--mc",   default="TnPTree_mc.root",   help="MC input filename")
    parser.add_option("-d", "--data", default="TnPTree_data.root", help="Data input filename")
    parser.add_option("-n", "--no-negative-weights", dest="noNegative", action="store_true", help="Remove negative weighted events", default=False)

    (options, arg) = parser.parse_args()
    
    main(options)

