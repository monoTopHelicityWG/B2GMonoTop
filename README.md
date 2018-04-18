# B2GMonoTop


## B2GMonoTop2016 recipe (ReReco data + Summer16 MC):
```bash
setenv SCRAM_ARCH slc6_amd64_gcc530
cmsrel CMSSW_8_0_26
cd CMSSW_8_0_26/src/
cmsenv
git cms-init

git cms-merge-topic -u cms-met:fromCMSSW_8_0_20_postICHEPfilter
git cms-merge-topic cms-met:METRecipe_8020
git clone https://github.com/monoTopHelicityWG/B2GMonoTop.git Analysis/B2GMonoTop
git clone git@github.com:cms-jet/JetToolbox.git JMEAnalysis/JetToolbox -b jetToolbox_80X
git clone git@github.com:cms-jet/JECDatabase.git JMEAnalysis/JECDatabase
git clone git@github.com:cms-jet/JRDatabase.git JMEAnalysis/JRDatabase
git clone git@github.com:thaarres/PuppiSoftdropMassCorr.git JMEAnalysis/PuppiSoftdropMassCorr
git cms-merge-topic Sam-Harper:HEEPV70VID_8010_ReducedCheckout 
scramv1 b -j 16
cd Analysis/B2GTTbar/test/
cmsRun run_B2GTTbarTreeMaker_MonoTopNew_Toolbox.py
```
Based on:
https://github.com/cmsb2g/B2GTTbar/blob/master/plugins/B2GTTbarTreeMaker.cc
