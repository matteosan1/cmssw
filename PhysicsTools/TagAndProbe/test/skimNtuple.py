import ROOT
import json 
from optparse import OptionParser

data = ""

def readJSON(jsonName):
    global data
    json_data=open(jsonName).read()
    data = json.loads(json_data)

def main(options):
    global data
    file = ROOT.TFile(options.input)
    directory = file.Get(options.directory)
    directory.cd()
    tree = directory.Get("fitter_tree")
    entries = tree.GetEntries()  

    #--- Write to new file
    outFile = options.input.split(".root")[0]+"_slim.root"
    newFile = ROOT.TFile(outFile, "RECREATE")
    directory_new = newFile.mkdir(options.directory)
    directory_new.cd()
    tree_new = tree.CloneTree(0)

    for z in range(entries):
        tree.GetEntry(z)
    
        passCut = False
        if (options.isMC):
            #--- Only write out certain events that pass some cut
            if (tree.event%10 == 0):
                passCut = True
        else:
            if (unicode(tree.run) in data):
                for r in data[unicode(tree.run)]:
                    if (tree.lumi in xrange(r[0], r[1])):
                        passCut = True
                        break
        
        if (passCut):
            tree_new.Fill()
    
    tree_new.GetCurrentFile().Write()
    tree_new.GetCurrentFile().Close() 


if __name__ == "__main__":  
    parser = OptionParser()
    parser.add_option("-i", "--input",     default="TnPTree_mc.root",           help="Input filename")
    parser.add_option("-d", "--directory", default="GsfElectronToRECO",         help="Directory with tree")
    parser.add_option("--mc", dest="isMC", action="store_true", help="MC file or not", default=False)
    parser.add_option("--json",            default="json_DCSONLY_Run2015B.txt", help="Name of JSON file")

    (options, arg) = parser.parse_args()
     
    if (not options.isMC):
        readJSON(options.json)
    
    main(options)

