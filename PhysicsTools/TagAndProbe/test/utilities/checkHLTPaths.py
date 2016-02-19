# import ROOT in batch mode
import sys
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

# load FWLite C++ libraries
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.AutoLibraryLoader.enable()

# load FWlite python libraries
from DataFormats.FWLite import Handle, Events

triggerBits, triggerBitLabel = Handle("edm::TriggerResults"), ("TriggerResults","","HLT")
triggerObjects, triggerObjectLabel  = Handle("std::vector<pat::TriggerObjectStandAlone>"), "selectedPatTrigger"
triggerPrescales, triggerPrescaleLabel  = Handle("pat::PackedTriggerPrescales"), "patTrigger"

# open file (you can use 'edmFileUtil -d /store/whatever.root' to get the physical file name)
events = Events(sys.argv[1])

for iev,event in enumerate(events):
    event.getByLabel(triggerBitLabel, triggerBits)
    event.getByLabel(triggerObjectLabel, triggerObjects)
    #event.getByLabel(triggerPrescaleLabel, triggerPrescales)

    print "\nEvent %d: run %6d, lumi %4d, event %12d" % (iev,event.eventAuxiliary().run(), event.eventAuxiliary().luminosityBlock(),event.eventAuxiliary().event())
    print "\n === TRIGGER PATHS ==="
    names = event.object().triggerNames(triggerBits.product())
    for i in xrange(triggerBits.product().size()):
        #print "Trigger ", names.triggerName(i), ", prescale ", triggerPrescales.product().getPrescaleForIndex(i), ": ", ("PASS" if triggerBits.product().accept(i) else "fail (or not run)") 
        print "Trigger ", names.triggerName(i), ": ", ("PASS" if triggerBits.product().accept(i) else "fail (or not run)") 

    print "\n === TRIGGER OBJECTS ==="
    for j,to in enumerate(triggerObjects.product()):
        to.unpackPathNames(names);
        print "Trigger object pt %6.2f eta %+5.3f phi %+5.3f  " % (to.pt(),to.eta(),to.phi())
        print "         collection: ", to.collection()
        print "         type ids: ", ", ".join([str(f) for f in to.filterIds()])
        print "         filters: ", ", ".join([str(f) for f in to.filterLabels()])
        pathslast = set(to.pathNames(True))
        print "         paths:   ", ", ".join([("%s*" if f in pathslast else "%s")%f for f in to.pathNames()]) 
    if iev > 10: break
