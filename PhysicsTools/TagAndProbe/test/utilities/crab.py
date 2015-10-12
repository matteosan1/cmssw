from WMCore.Configuration import Configuration
    
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.workArea = 'data'
config.General.transferLogs=True

config.section_('JobType')
config.JobType.psetName = 'makeTreeData.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['TnPTree_data.root']
config.section_('Data')

config.Data.inputDataset= '/SingleElectron/Run2015B-PromptReco-v1/MINIAOD'
config.Data.unitsPerJob = 10
config.Data.splitting = 'FileBased'
#config.Data.lumiMask = "https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/DCSOnly/json_DCSONLY.txt"
#config.Data.runRange = '251562,251561,251559'
config.Data.publication = False
config.Data.useParent = False
config.Data.outLFNDirBase = '/store/user/sani/data/2015/tnp'

config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_US_UCSD'
