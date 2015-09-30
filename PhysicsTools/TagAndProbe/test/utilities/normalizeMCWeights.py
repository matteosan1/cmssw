import ROOT
import array
from optparse import OptionParser

def main(options):
    file = ROOT.TFile(options.input)
    directory = file.Get(options.directory)
    directory.cd()
    tree = directory.Get("fitter_tree")
    entries = tree.GetEntries()  

    #--- Compute sum of MC weights
    tree.SetBranchStatus("*", 0)
    tree.SetBranchStatus("weight", 1)
    sumWeight = 0
    weightedEntries = 0
    for z in range(entries):
        tree.GetEntry(z)
        if (tree.weight >= 0):
            weightedEntries = weightedEntries + 1.
            sumWeight = sumWeight + tree.weight
            
    sumWeight = sumWeight/weightedEntries
    #--- Write to new file
    tree.SetBranchStatus("*", 1)
    outFile = options.input.split(".root")[0]+"_norm.root"
    newFile = ROOT.TFile(outFile, "RECREATE")
    directory_new = newFile.mkdir(options.directory)
    directory_new.cd()
    tree_new = tree.CloneTree(0)
    b_totWeight = array.array('f', [0])
    tree_new.Branch("totWeight", b_totWeight, "totWeight/F")

    for z in range(entries):
        tree.GetEntry(z)
        
        if (tree.weight>= 0):
            b_totWeight[0] = tree.weight/sumWeight*tree.PUweight
            tree_new.Fill()
    
    tree_new.GetCurrentFile().Write()
    tree_new.GetCurrentFile().Close() 


if __name__ == "__main__":  
    parser = OptionParser()
    parser.add_option("-i", "--input",     default="TnPTree_mc.root",           help="Input filename")
    parser.add_option("-d", "--directory", default="GsfElectronToRECO",         help="Directory with tree")

    (options, arg) = parser.parse_args()
    
    main(options)

