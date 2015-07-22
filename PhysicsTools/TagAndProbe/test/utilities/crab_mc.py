from WMCore.Configuration import Configuration

config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.workArea = 'data'
config.General.transferLogs=True

config.section_('JobType')
config.JobType.psetName = 'makeTree.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['TnPTree_mc.root']
config.section_('Data')

config.Data.inputDataset= '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/MINIAODSIM'
config.Data.unitsPerJob = 10
config.Data.splitting = 'FileBased'
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/DCSOnly/json_DCSONLY.txt'
#config.Data.runRange = '251562,251561,251559'
config.Data.publication = False
config.Data.useParent = False
config.Data.outLFNDirBase = '/store/user/sani/mc/2015/tnp'

config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_US_UCSD'
    
