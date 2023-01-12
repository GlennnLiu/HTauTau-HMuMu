# HTauTau-HMuMu

Measure Y(tautau)/Y(mumu) to constrain lepton universality violation, with the Run2 UltraLegacy samples. More is being developed.

The CMSSW version is **CMSSW_10_6_26**.
Do the following:
```
cmsrel CMSSW_10_6_26
cd CMSSW_10_6_26/src/
cmsenv

git cms-init

# New Jet PU ID: dedicated training for each year
git cms-addpkg  RecoJets/JetProducers

# STXS Categorisation: now directly implemented in CMSSW
git cms-addpkg GeneratorInterface/RivetInterface
git cms-addpkg SimDataFormats/HTXS


# Updated for UL. See: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaUL2016To2018
git cms-addpkg RecoEgamma/EgammaTools
git clone https://github.com/cms-egamma/EgammaPostRecoTools.git
mv EgammaPostRecoTools/python/EgammaPostRecoTools.py RecoEgamma/EgammaTools/python/.
git cms-addpkg EgammaAnalysis/ElectronTools
(rm -rf EgammaAnalysis/ElectronTools/data;git clone https://github.com/jainshilpi/EgammaAnalysis-ElectronTools.git -b ULSSfiles_correctScaleSysMC EgammaAnalysis/ElectronTools/data;)

#### Please do not add any custom (non-CMSSW) package before this line ####
git clone git@github.com:GlennnLiu/HTauTauHMuMu.git -b Run2_UL

