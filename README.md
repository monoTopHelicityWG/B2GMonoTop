# B2GMonoTop


##B2GMonoTop2016 recipe (ReReco data + Summer16 MC):
```
setenv SCRAM_ARCH slc6_amd64_gcc530
cmsrel CMSSW_8_0_26
cd CMSSW_8_0_26/src/
cmsenv
git cms-init
git cms-merge-topic gkasieczka:test-httv2-8014 #adding HEP TT v2 https://twiki.cern.ch/twiki/bin/view/CMS/JetTopTagging
git cms-merge-topic -u cms-met:fromCMSSW_8_0_20_postICHEPfilter
git cms-merge-topic cms-met:METRecipe_8020
git cms-merge-topic Sam-Harper:HEEPV70VID_8010_ReducedCheckout 
git cms-merge-topic ikrav:egm_id_80X_v3 
mkdir -p ../external/slc6_amd64_gcc530/data/RecoEgamma/ElectronIdentification/ 
git clone git@github.com:cms-data/RecoEgamma-ElectronIdentification ../external/slc6_amd64_gcc530/data/RecoEgamma/ElectronIdentification/data 
git cms-addpkg CondFormats/BTauObjects
git clone https://github.com/rappoccio/PredictedDistribution.git Analysis/PredictedDistribution
git clone https://github.com/monoTopHelicityWG/B2GMonoTop.git Analysis/B2GMonoTop
git clone git@github.com:cms-jet/JetToolbox.git JMEAnalysis/JetToolbox -b jetToolbox_80X
git clone git@github.com:cms-jet/JECDatabase.git JMEAnalysis/JECDatabase
git clone git@github.com:cms-jet/JRDatabase.git JMEAnalysis/JRDatabase
git clone git@github.com:thaarres/PuppiSoftdropMassCorr.git JMEAnalysis/PuppiSoftdropMassCorr
scramv1 b -j 16
cd Analysis/B2GTTbar/test/
cmsRun run_B2GTTbarTreeMaker_MC_testLocal_ZP2000w200.py


**For analysis with Loop tree, you must also checkout the btag package and remove one line:**
cd CMSSW_8_0_26/src/
*Edit CondFormats/BTauObjects/src/classes.h  and comment out line 30 (it says BTagCalibration btc1;)* 
cd CondFormats/BTauObjects 
scramv1 b
cd ../../Analysis/B2GTTbar/test/
*run loop tree*
```

https://github.com/cmsb2g/B2GTTbar/blob/master/plugins/B2GTTbarTreeMaker.cc
