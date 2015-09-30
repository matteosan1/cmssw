import FWCore.ParameterSet.Config as cms

options = dict()

options['HLTProcessName']          = "HLT"
options['INPUT_FILE_NAME']         = "/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/00C4781D-6B08-E511-8A0A-0025905A6084.root"
options['OUTPUT_FILE_NAME']        = "TnPTree_mc.root"
options['ELECTRON_COLL']           = "slimmedElectrons"
options['ELECTRON_CUTS']           = "(abs(superCluster.eta)<2.5) && (ecalEnergy*sin(superClusterPosition.theta)>10.0)"
options['ELECTRON_TAG_CUTS']       = "(abs(superCluster.eta)<=2.5) && !(1.4442<=abs(superCluster.eta)<=1.566) && pt >= 25.0"
options['SUPERCLUSTER_COLL']       = "reducedEgamma:reducedSuperClusters"
options['SUPERCLUSTER_CUTS']       = "abs(eta)<2.5 && !(1.4442< abs(eta) <1.566) && et>10.0"
options['TnPPATHS']                = cms.vstring("HLT_Ele23_WP75_Gsf_v*")
options['TnPHLTTagFilters']        = cms.vstring("hltEle23WP75GsfTrackIsoFilter")
options['TnPHLTProbeFilters']      = cms.vstring()
options['HLTFILTERTOMEASURE']      = cms.vstring("hltEle27WP75GsfTrackIsoFilter")
options['GLOBALTAG']               = 'MCRUN2_74_V9A'
options['EVENTSToPROCESS']         = cms.untracked.VEventRange()
options['MAXEVENTS']               = cms.untracked.int32(1000) 
options['useAOD']                  = cms.bool(False)
options['DOTRIGGER']               = cms.bool(True)
options['DORECO']                  = cms.bool(True)
options['DOID']                    = cms.bool(True)
options['OUTPUTEDMFILENAME']       = 'edmFile.root'
options['DEBUG']                   = cms.bool(False)
