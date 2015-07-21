import ROOT
import json 

json_data=open("json_DCSONLY_Run2015B.txt").read()
data = json.loads(json_data)

for k in data.keys():
    print type(k)

file = ROOT.TFile("TnPTree_mc.root")
#directory = file.Get("GsfElectronToSC")
directory = file.Get("GsfElectronToRECO")
directory.cd()
tree = directory.Get("fitter_tree")
entries = tree.GetEntries()  

#--- Write to new file
outFile = "TnPTree_mc_slim.root"
newFile = ROOT.TFile(outFile, "RECREATE")
#directory_new = newFile.mkdir("GsfElectronToSC")
directory_new = newFile.mkdir("GsfElectronToRECO")
directory_new.cd()
tree_new = tree.CloneTree(0)


for z in range(entries):
    tree.GetEntry(z)
    #--- Only write out certain events that pass some cut
    if (tree.event%10 == 0):
        tree_new.Fill()
    
    isGood = False
    #if (unicode(tree.run) in data):
    #    for r in data[unicode(tree.run)]:
    #        if (tree.lumi in xrange(r[0], r[1])):
    #            isGood = True
    #            break
    #
    #if (isGood):
    #tree_new.Fill()
    
# use GetCurrentFile just in case we went over the
# (customizable) maximum file size
tree_new.GetCurrentFile().Write()
tree_new.GetCurrentFile().Close() 

