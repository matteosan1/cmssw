import FWCore.ParameterSet.Config as cms

options = dict()

options['MC_FLAG']             = cms.bool(False)
options['HLTProcessName']      = "HLT"
options['INPUT_FILE_NAME']     = "/store/relval/CMSSW_7_4_1/DoubleElectron/MINIAOD/GR_R_74_V12A_RelVal_zEl2012D-v1/00000/F4073AF3-85EC-E411-8117-0025905B8596.root"
options['OUTPUT_FILE_NAME']    = "TnPTree.root"
options['ELECTRON_COLL']       = "slimmedElectrons"
options['ELECTRON_CUTS']       = "(abs(superCluster.eta)<2.5) && (ecalEnergy*sin(superClusterPosition.theta)>10.0)"
options['ELECTRON_TAG_CUTS']   = "(abs(superCluster.eta)<=2.5) && !(1.4442<=abs(superCluster.eta)<=1.566) && pt >= 25.0"
options['SUPERCLUSTER_COLL']   = "reducedEgamma:reducedSuperClusters"
options['SUPERCLUSTER_CUTS']   = "abs(eta)<2.5 && !(1.4442< abs(eta) <1.566) && et>10.0"
options['TnPPATHS']            = []
options['TnPHLTTagFilters']    = ["hltEle20WP60Ele8TrackIsoFilter", "hltEle25WP60SC4TrackIsoFilter"]
options['TnPHLTProbeFilters']  = ["hltEle20WP60Ele8PixelMatchUnseededFilter", "hltEle25WP60SC4EtUnseededFilter"]
options['HLTPathToMeasure']    = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v15"
options['HLTFILTERTOMEASURE']  = "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter"
options['GLOBALTAG']           = 'GR_R_74_V8'
options['EVENTSToPROCESS']     = cms.untracked.VEventRange('1:28065077','1:28065427','1:28033364','1:54509262')
options['MAXEVENTS']           = cms.untracked.int32(-1) 
options['useAOD']              = cms.bool(False)
options['DOTRIGGER']           = cms.bool(False)
options['DORECO']              = cms.bool(False)
options['DOID']                = cms.bool(True)
options['OUTPUTEDMFILENAME']   = 'edmFile.root'
options['DEBUG']               = cms.bool(False)
