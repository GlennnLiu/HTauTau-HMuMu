# HTauTauHMuMu

Measure Y(tautau)/Y(mumu) to constrain lepton universality violation. More is being developed.

The CMSSW version is **CMSSW_10_2_22**.
Do the following:
```
cmsrel CMSSW_10_2_22
cd CMSSW_10_2_22/src/
cmsenv

git cms-init

#Preliminary electron scale and smearing corrections according to https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPostRecoRecipes#2018_Preliminary_Energy_Correcti
#We need the ElectronTools package to calculate smear and scale uncertainties so just download the ScaleAndSmearing files manualy 
git cms-merge-topic cms-egamma:EgammaPostRecoTools
git cms-merge-topic cms-egamma:PhotonIDValueMapSpeedup1029
git cms-merge-topic cms-egamma:slava77-btvDictFix_10210
git cms-addpkg EgammaAnalysis/ElectronTools
(rm -rf EgammaAnalysis/ElectronTools/data;git clone https://github.com/cms-data/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data;)

# New Jet PU ID: dedicated training for each year
git cms-addpkg  RecoJets/JetProducers

# 2016 and 2018 retraining for electron BDT
git cms-merge-topic mkovac:Electron_XGBoost_MVA_2016_and_2018_CMSSW_10_2_15

#MET corrections according to https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription#Instructions_for_9_4_X_X_0_for_M
git cms-merge-topic cms-met:METFixEE2017_949_v2_backport_to_102X

# Muon MVA
git clone https://github.com/bonanomi/MuonMVAReader.git MuonMVAReader

#kinematic refitting
git clone https://github.com/mhl0116/KinZfitter-1.git KinZfitter
(cd KinZfitter ; git checkout -b from-27daebb 27daebb)

# SVfit
git clone https://github.com/LLRCMS/ClassicSVfit.git TauAnalysis/ClassicSVfit -b bbtautau_LegacyRun2
git clone https://github.com/svfit/SVfitTF TauAnalysis/SVfitTF

#Add TauPOG corrections (TES and EES)
git clone https://github.com/cms-tau-pog/TauIDSFs TauPOG/TauIDSFs

#Add DeepTau code from Tau POG repository (note "-u" option preventing checkout of unnecessary stuff)
git cms-merge-topic -u cms-tau-pog:CMSSW_10_2_X_tau-pog_DeepTau2017v2p1_nanoAOD

# KinFit
git clone git@github.com:LLRCMS/HHKinFit2.git -b bbtautau_LegacyRun2

scramv1 b -j 8

cd HHKinFit2/
ln -ns interface include
source setup.sh
./compile.sh
cd ..

# This package
git clone git@github.com:GlennnLiu/HTauTauHMuMu.git -b Run2_Legacy
```

