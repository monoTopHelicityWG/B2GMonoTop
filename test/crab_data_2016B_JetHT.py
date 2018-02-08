from CRABClient.UserUtilities import config, getUsernameFromSiteDB
import os

config = config()
#config.General.requestName = 'monotopTreeV6_ST_t-channel_top_4f_inclusiveDecays'
config.General.requestName = 'monotopTree_JetHT_Run2016B-23Sep2016-v3_v9'
config.General.transferOutputs = True
config.General.transferLogs = True
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_B2GTTbarTreeMaker_MonoTopNew_Toolbox_data.py'
#config.JobType.maxJobRuntimeMin = 2750
config.JobType.maxMemoryMB = 2500

#config.Data.userInputFiles = open('mt_mu_2tev_RH.txt').readlines()
config.Data.inputDBS = 'global'
config.JobType.inputFiles = [

'PUweight_FinalJSON2016_PileupJSONFeb2017_Xsec69200nominal_MCRunIISummer16MiniAODv2_PUMoriond17.root',
os.environ['CMSSW_BASE'] + '/src/JMEAnalysis/JECDatabase/textFiles/Summer16_23Sep2016BCDV6_DATA/Summer16_23Sep2016BCDV6_DATA_L1FastJet_AK8PFchs.txt',
os.environ['CMSSW_BASE'] + '/src/JMEAnalysis/JECDatabase/textFiles/Summer16_23Sep2016BCDV6_DATA/Summer16_23Sep2016BCDV6_DATA_L2Relative_AK8PFchs.txt',
os.environ['CMSSW_BASE'] + '/src/JMEAnalysis/JECDatabase/textFiles/Summer16_23Sep2016BCDV6_DATA/Summer16_23Sep2016BCDV6_DATA_L3Absolute_AK8PFchs.txt',
os.environ['CMSSW_BASE'] + '/src/JMEAnalysis/JECDatabase/textFiles/Summer16_23Sep2016BCDV6_DATA/Summer16_23Sep2016BCDV6_DATA_L2L3Residual_AK8PFchs.txt',
os.environ['CMSSW_BASE'] + '/src/JMEAnalysis/JECDatabase/textFiles/Summer16_23Sep2016BCDV6_DATA/Summer16_23Sep2016BCDV6_DATA_Uncertainty_AK8PFchs.txt',
os.environ['CMSSW_BASE'] + '/src/JMEAnalysis/JECDatabase/textFiles/Summer16_23Sep2016BCDV6_DATA/Summer16_23Sep2016BCDV6_DATA_L1FastJet_AK4PFchs.txt',
os.environ['CMSSW_BASE'] + '/src/JMEAnalysis/JECDatabase/textFiles/Summer16_23Sep2016BCDV6_DATA/Summer16_23Sep2016BCDV6_DATA_L2Relative_AK4PFchs.txt',
os.environ['CMSSW_BASE'] + '/src/JMEAnalysis/JECDatabase/textFiles/Summer16_23Sep2016BCDV6_DATA/Summer16_23Sep2016BCDV6_DATA_L3Absolute_AK4PFchs.txt',
os.environ['CMSSW_BASE'] + '/src/JMEAnalysis/JECDatabase/textFiles/Summer16_23Sep2016BCDV6_DATA/Summer16_23Sep2016BCDV6_DATA_L2L3Residual_AK4PFchs.txt',
os.environ['CMSSW_BASE'] + '/src/JMEAnalysis/JECDatabase/textFiles/Summer16_23Sep2016BCDV6_DATA/Summer16_23Sep2016BCDV6_DATA_Uncertainty_AK4PFchs.txt',
os.environ['CMSSW_BASE'] + '/src/JMEAnalysis/JECDatabase/textFiles/Summer16_23Sep2016BCDV6_DATA/Summer16_23Sep2016BCDV6_DATA_L1FastJet_AK8PFPuppi.txt',
os.environ['CMSSW_BASE'] + '/src/JMEAnalysis/JECDatabase/textFiles/Summer16_23Sep2016BCDV6_DATA/Summer16_23Sep2016BCDV6_DATA_L2Relative_AK8PFPuppi.txt',
os.environ['CMSSW_BASE'] + '/src/JMEAnalysis/JECDatabase/textFiles/Summer16_23Sep2016BCDV6_DATA/Summer16_23Sep2016BCDV6_DATA_L3Absolute_AK8PFPuppi.txt',
os.environ['CMSSW_BASE'] + '/src/JMEAnalysis/JECDatabase/textFiles/Summer16_23Sep2016BCDV6_DATA/Summer16_23Sep2016BCDV6_DATA_L2L3Residual_AK8PFPuppi.txt',
os.environ['CMSSW_BASE'] + '/src/JMEAnalysis/JECDatabase/textFiles/Summer16_23Sep2016BCDV6_DATA/Summer16_23Sep2016BCDV6_DATA_Uncertainty_AK8PFPuppi.txt',
os.environ['CMSSW_BASE'] + '/src/JMEAnalysis/JECDatabase/textFiles/Summer16_23Sep2016BCDV6_DATA/Summer16_23Sep2016BCDV6_DATA_L1FastJet_AK4PFPuppi.txt',
os.environ['CMSSW_BASE'] + '/src/JMEAnalysis/JECDatabase/textFiles/Summer16_23Sep2016BCDV6_DATA/Summer16_23Sep2016BCDV6_DATA_L2Relative_AK4PFPuppi.txt',
os.environ['CMSSW_BASE'] + '/src/JMEAnalysis/JECDatabase/textFiles/Summer16_23Sep2016BCDV6_DATA/Summer16_23Sep2016BCDV6_DATA_L3Absolute_AK4PFPuppi.txt',
os.environ['CMSSW_BASE'] + '/src/JMEAnalysis/JECDatabase/textFiles/Summer16_23Sep2016BCDV6_DATA/Summer16_23Sep2016BCDV6_DATA_L2L3Residual_AK4PFPuppi.txt',
os.environ['CMSSW_BASE'] + '/src/JMEAnalysis/JECDatabase/textFiles/Summer16_23Sep2016BCDV6_DATA/Summer16_23Sep2016BCDV6_DATA_Uncertainty_AK4PFPuppi.txt',
os.environ['CMSSW_BASE'] + '/src/JMEAnalysis/JRDatabase/textFiles/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA_SF_AK8PFchs.txt',

os.environ['CMSSW_BASE'] + '/src/JMEAnalysis/JRDatabase/textFiles/Spring16_25nsV10_MC/Spring16_25nsV10_MC_PtResolution_AK8PFchs.txt',
os.environ['CMSSW_BASE'] + '/src/JMEAnalysis/JRDatabase/textFiles/Spring16_25nsV10_MC/Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt'
]
config.Data.splitting = 'FileBased'
#config.Data.inputDataset = '/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM'
config.Data.inputDataset = '/JetHT/Run2016B-23Sep2016-v3/MINIAOD'
config.Data.inputDBS = 'global'

config.Data.unitsPerJob = 1
config.Data.totalUnits = 100
config.Data.outLFNDirBase = '/store/user/rymuelle/MonoTop/Data/'
config.Data.publication = False
config.Site.storageSite = 'T3_US_FNALLPC'

config.Data.ignoreLocality = True
config.Site.whitelist = ['T3_US_FNALLPC']
config.Site.ignoreGlobalBlacklist = True
