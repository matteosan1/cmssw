import ROOT
from optparse import OptionParser
import sys

def psetDumpFromTTree(options):
    
    f = ROOT.TFile(options.input)
    t = f.Get(options.directory+"/fitter_tree")
    md = t.GetUserInfo()
    md.Print()

if (__name__ == "__main__"):
    parser = OptionParser()
    parser.add_option("-i", "--input", default="TnPTree_mc.root", help="Input filename")
    parser.add_option("-d", "--directory", default="GsfElectronToRECO", help="Directory with workspace")

    (options, arg) = parser.parse_args()

    psetDumpFromTTree(options)
