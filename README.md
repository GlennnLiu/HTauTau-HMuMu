# HTauTau-HMuMu

Measure Y(tautau)/Y(mumu) to constrain lepton universality violation, with the Run2 UltraLegacy samples. In the beginning we would like to use trigger-and-probe method, which trigger on associated objects of the Higgs boson (VBF or weak boson). The branch Run2_UL_ZZ is the framework for a side-studies: ZZ process, with one Z to be triggered on decaying to ee or mumu, the other to be probed decaying to tautau or mumu.
While following the discussions, since there's no inclusive VBF triggers for Run2, purly triggering on VBF jets is not possible; triggering on leptonically decaying weak bosons would lead to a too small total branching fraction. We decide to do the following:

Four processes: H->tautau, H->mumu, Z->tautau, Z->mumu. H/Z->tautau are triggered on a series of HLT paths but they should be the same between H and Z; the same for H/Z->mumu. The final POI is $\frac{BR(H->mumu)/BR(Z->mumu)}{BR(H->tautau)/BR(H->tautau)}$. In this way, we can cancel some of the trigger efficiency systematics and well as reco/selection systematics.

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

#Add TauPOG corrections (TES and EES)
git clone https://github.com/cms-tau-pog/TauIDSFs TauPOG/TauIDSFs

# SVfit
git clone https://github.com/LLRCMS/ClassicSVfit.git TauAnalysis/ClassicSVfit -b bbtautau_LegacyRun2
git clone https://github.com/svfit/SVfitTF TauAnalysis/SVfitTF

# This package
git clone git@github.com:GlennnLiu/HTauTauHMuMu.git -b Run2_UL

# Muon MVA
git cms-addpkg CondFormats/EgammaObjects
git cms-addpkg CommonTools/MVAUtils
git clone https://github.com/bonanomi/MuonMVAReader.git MuonMVAReader
(cd MuonMVAReader; git checkout 3d53269)

#NanoAODTools
git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
(cd PhysicsTools/NanoAODTools ; git checkout -b from-c32f055 c32f055)

# compile
scramv1 b -j 8
```
