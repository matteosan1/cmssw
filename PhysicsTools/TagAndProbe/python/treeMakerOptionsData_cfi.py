import FWCore.ParameterSet.Config as cms

options = dict()

options['MC_FLAG']             = cms.bool(False)
options['HLTProcessName']      = "HLT"
options['INPUT_FILE_NAME']     = "/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/244/00000/12EE24E2-8F27-E511-80D1-02163E013793.root"
options['OUTPUT_FILE_NAME']    = "TnPTree_data.root"
options['ELECTRON_COLL']       = "slimmedElectrons"
options['ELECTRON_CUTS']       = "(abs(superCluster.eta)<2.5) && (ecalEnergy*sin(superClusterPosition.theta)>10.0)"
options['ELECTRON_TAG_CUTS']   = "(abs(superCluster.eta)<=2.5) && !(1.4442<=abs(superCluster.eta)<=1.566) && pt >= 25.0"
options['SUPERCLUSTER_COLL']   = "reducedEgamma:reducedSuperClusters"
options['SUPERCLUSTER_CUTS']   = "abs(eta)<2.5 && !(1.4442< abs(eta) <1.566) && et>10.0"
options['TnPPATHS']            = ["HLT_Ele27_eta2p1_WPLoose_Gsf_v1",]
options['TnPHLTTagFilters']    = ["hltEle27WPLooseGsfTrackIsoFilter",]
options['TnPHLTProbeFilters']  = ["*"]
options['HLTPathToMeasure']    = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v15"
options['HLTFILTERTOMEASURE']  = "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter"
options['GLOBALTAG']           = 'GR_P_V56'
options['EVENTSToPROCESS']     = cms.untracked.VEventRange()
options['MAXEVENTS']           = cms.untracked.int32(-1) 
options['useAOD']              = cms.bool(False)
options['DOTRIGGER']           = cms.bool(False)
options['DORECO']              = cms.bool(True)
options['DOID']                = cms.bool(True)
options['OUTPUTEDMFILENAME']   = 'edmFile.root'
options['DEBUG']               = cms.bool(False)
