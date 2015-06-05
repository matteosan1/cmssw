import FWCore.ParameterSet.Config as cms

options = dict()

options['MC_FLAG']               = cms.bool(True)
options['HLTProcessName']        = "HLT"
options['INPUT_FILE_NAME']       = "/store/relval/CMSSW_7_4_1/RelValZEE_13/MINIAODSIM/MCRUN2_74_V9_gensim_740pre7-v1/00000/1E35CCF8-32EC-E411-8F29-0025905A48D0.root"
options['OUTPUT_FILE_NAME']      = "TnPTree.root"
options['PHOTON_COLL']           = "slimmedPhotons"
options['PHOTON_CUTS']           = "(abs(superCluster.eta)<2.5) && ((superCluster.energy*sin(superCluster.position.theta))>15.0) && !passElectronVeto"
options['PHOTON_TAG_CUTS']       = "(abs(superCluster.eta)<=2.5) && !(1.4442<=abs(superCluster.eta)<=1.566) && (superCluster.energy*sin(superCluster.position.theta))>25.0 && !passElectronVeto"
options['TnPPATHS']              = cms.vstring("HLT_Ele20WP60_Ele8_Mass55_v*", "HLT_Ele25WP60_SC4_Mass55_v*")
options['TnPHLTTagFilters']      = cms.vstring("hltEle20WP60Ele8TrackIsoFilter", "hltEle25WP60SC4TrackIsoFilter")
options['TnPHLTProbeFilters']    = cms.vstring("hltEle20WP60Ele8PixelMatchUnseededFilter", "hltEle25WP60SC4EtUnseededFilter")
options['HLTPathToMeasure']      = cms.vstring("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v15")
options['HLTFILTERTOMEASURE']    = cms.vstring("hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter")
options['GLOBALTAG']             = 'PHYS14_25_V1::All'
options['EVENTSToPROCESS']       = cms.untracked.VEventRange()
options['MAXEVENTS']             = cms.untracked.int32(-1) 
options['useAOD']                = cms.bool(False)
options['OUTPUTEDMFILENAME']     = 'edmFile.root'
options['DEBUG']                 = cms.bool(False)
