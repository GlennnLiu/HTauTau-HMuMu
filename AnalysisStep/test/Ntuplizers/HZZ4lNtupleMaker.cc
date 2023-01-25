// -*- C++ -*-
//
//
// Fill a tree for selected candidates.
//


// system include files
#include <cassert>
#include <memory>
#include <cmath>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Framework/interface/LuminosityBlock.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <DataFormats/Common/interface/TriggerResults.h>
#include <FWCore/Common/interface/TriggerNames.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>

#include <DataFormats/Common/interface/View.h>
#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h" //Atbbf
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Tau.h>
#include <DataFormats/PatCandidates/interface/Photon.h>
#include <DataFormats/PatCandidates/interface/Jet.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/METReco/interface/PFMET.h>
#include <DataFormats/METReco/interface/PFMETCollection.h>
#include <DataFormats/JetReco/interface/PFJet.h>
#include <DataFormats/JetReco/interface/PFJetCollection.h>
#include <DataFormats/Math/interface/LorentzVector.h>
#include <CommonTools/UtilAlgos/interface/TFileService.h>
#include <DataFormats/Common/interface/MergeableCounter.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include <SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h>
#include <SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h>

#include <HTauTauHMuMu/AnalysisStep/interface/DaughterDataHelpers.h>
#include <HTauTauHMuMu/AnalysisStep/interface/FinalStates.h>
#include <HTauTauHMuMu/AnalysisStep/interface/MCHistoryTools.h>
#include <HTauTauHMuMu/AnalysisStep/interface/PileUpWeight.h>
#include "SimDataFormats/HTXS/interface/HiggsTemplateCrossSections.h"

//ATjets Additional libraries for GenJet variables
#include <DataFormats/PatCandidates/interface/Jet.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectorParameters.h>
#include <JetMETCorrections/Objects/interface/JetCorrectionsRecord.h>
#include <JetMETCorrections/Modules/interface/JetResolution.h>

#include "HTauTauHMuMu/AnalysisStep/interface/EwkCorrections.h"
#include "HTauTauHMuMu/AnalysisStep/src/kFactors.C"
#include <HTauTauHMuMu/AnalysisStep/interface/bitops.h>
#include <HTauTauHMuMu/AnalysisStep/interface/LeptonIsoHelper.h>
#include <HTauTauHMuMu/AnalysisStep/interface/PhotonIDHelper.h>
#include <HTauTauHMuMu/AnalysisStep/interface/JetCleaner.h>
#include <HTauTauHMuMu/AnalysisStep/interface/utils.h>
#include <HTauTauHMuMu/AnalysisStep/interface/miscenums.h>
#include <HTauTauHMuMu/AnalysisStep/interface/ggF_qcd_uncertainty_2017.h>
#include <HTauTauHMuMu/AnalysisStep/interface/LeptonSFHelper.h>

#include <HTauTauHMuMu/AnalysisStep/interface/SVfit.h>
#include <TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h>
#include <TauAnalysis/ClassicSVfit/interface/ClassicSVfitIntegrand.h>
#include <TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h>
#include <TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h>

#include "ZZ4lConfigHelper.h"
#include "HZZ4lNtupleFactory.h"

#include <TRandom3.h>
#include <TH2D.h>
#include "TLorentzVector.h"
#include "TSpline.h"
#include "TGraphErrors.h"

#include <string>

bool verbose = false; //ATbbf

namespace {
  bool writeJets = false;     // Write jets in the tree. FIXME: make this configurable
  bool writePhotons = false; // Write photons in the tree. FIXME: make this configurable
  bool addKinRefit = false;
  bool addZZKinfit = false;
  bool addSVfit = true;
  bool addVtxFit = false;
  bool addFSRDetails = false;
  bool addQGLInputs = true;
  bool addGenAngles = false;
  bool skipMuDataMCWeight = false; // skip computation of data/MC weight for mu 
  bool skipEleDataMCWeight = false; // skip computation of data/MC weight for ele
  bool skipFakeWeight = true;   // skip computation of fake rate weight for CRs
  bool skipHqTWeight = true;    // skip computation of hQT weight

  //List of variables with default values
  Int_t RunNumber  = 0;
  Long64_t EventNumber  = 0;
  Int_t LumiNumber  = 0;
  Short_t NRecoMu  = 0;
  Short_t NRecoEle  = 0;
  Short_t NRecoTau  = 0;
  Short_t Nvtx  = 0;
  Short_t NObsInt  = 0;
  Float_t NTrueInt  = 0;
  Float_t PUWeight  = 0;
  Float_t PUWeight_Up  = 0;
  Float_t PUWeight_Dn  = 0;

  Float_t KFactor_QCD_ggZZ_Nominal = 0;
  Float_t KFactor_QCD_ggZZ_PDFScaleDn = 0;
  Float_t KFactor_QCD_ggZZ_PDFScaleUp = 0;
  Float_t KFactor_QCD_ggZZ_QCDScaleDn = 0;
  Float_t KFactor_QCD_ggZZ_QCDScaleUp = 0;
  Float_t KFactor_QCD_ggZZ_AsDn = 0;
  Float_t KFactor_QCD_ggZZ_AsUp = 0;
  Float_t KFactor_QCD_ggZZ_PDFReplicaDn = 0;
  Float_t KFactor_QCD_ggZZ_PDFReplicaUp = 0;
  Float_t KFactor_EW_qqZZ = 0;
  Float_t KFactor_EW_qqZZ_unc = 0;
  Float_t KFactor_QCD_qqZZ_dPhi = 0;
  Float_t KFactor_QCD_qqZZ_M = 0;
  Float_t KFactor_QCD_qqZZ_Pt = 0;
  
  //-------------------------------
  //-------------MET object--------
  //-------------------------------
  Float_t GenMET = -99;
  Float_t GenMETPhi = -99;
  Float_t PFMET = -99;
  Float_t PFMETPhi = -99;
  Float_t METx = -99;
  Float_t METy = -99;
  TMatrixD covMET(2, 2);
  Float_t METxUPTES = -99;
  Float_t METyUPTES = -99;
  Float_t METxDOWNTES = -99;
  Float_t METyDOWNTES = -99;
  Float_t METxUPEES = -99;
  Float_t METyUPEES = -99;
  Float_t METxDOWNEES = -99;
  Float_t METyDOWNEES = -99;
  Float_t METxUPMES = -99;
  Float_t METyUPMES = -99;
  Float_t METxDOWNMES = -99;
  Float_t METyDOWNMES = -99;
  Float_t METxUPJES = -99;
  Float_t METyUPJES = -99;
  Float_t METxDOWNJES = -99;
  Float_t METyDOWNJES = -99;
  Float_t METxUPJER = -99;
  Float_t METyUPJER = -99;
  Float_t METxDOWNJER = -99;
  Float_t METyDOWNJER = -99;
    
  // MET with no HF
  //Float_t PFMETNoHF  =  -99;
  //Float_t PFMETNoHFPhi  =  -99;
  Short_t nCleanedJets  =  0;
  Short_t nCleanedJetsPt30  = 0;
  Short_t nCleanedJetsPt30_jesUp  = 0;
  Short_t nCleanedJetsPt30_jesUp_Total           = 0;
  Short_t nCleanedJetsPt30_jesUp_Abs             = 0;
  Short_t nCleanedJetsPt30_jesUp_Abs_year        = 0;
  Short_t nCleanedJetsPt30_jesUp_BBEC1           = 0;
  Short_t nCleanedJetsPt30_jesUp_BBEC1_year      = 0;
  Short_t nCleanedJetsPt30_jesUp_EC2             = 0;
  Short_t nCleanedJetsPt30_jesUp_EC2_year        = 0;
  Short_t nCleanedJetsPt30_jesUp_FlavQCD         = 0;
  Short_t nCleanedJetsPt30_jesUp_HF              = 0;
  Short_t nCleanedJetsPt30_jesUp_HF_year         = 0;
  Short_t nCleanedJetsPt30_jesUp_RelBal          = 0;
  Short_t nCleanedJetsPt30_jesUp_RelSample_year  = 0;
  Short_t nCleanedJetsPt30_jesDn  = 0;
  Short_t nCleanedJetsPt30_jesDn_Total           = 0;
  Short_t nCleanedJetsPt30_jesDn_Abs             = 0;
  Short_t nCleanedJetsPt30_jesDn_Abs_year        = 0;
  Short_t nCleanedJetsPt30_jesDn_BBEC1           = 0;
  Short_t nCleanedJetsPt30_jesDn_BBEC1_year      = 0;
  Short_t nCleanedJetsPt30_jesDn_EC2             = 0;
  Short_t nCleanedJetsPt30_jesDn_EC2_year        = 0;
  Short_t nCleanedJetsPt30_jesDn_FlavQCD         = 0;
  Short_t nCleanedJetsPt30_jesDn_HF              = 0;
  Short_t nCleanedJetsPt30_jesDn_HF_year         = 0;
  Short_t nCleanedJetsPt30_jesDn_RelBal          = 0;
  Short_t nCleanedJetsPt30_jesDn_RelSample_year  = 0;
  Short_t nCleanedJetsPt30_jerUp  = 0;
  Short_t nCleanedJetsPt30_jerDn  = 0;
  Short_t nCleanedJetsPt30BTagged  = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF  = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesUp  = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesUp_Total           = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesUp_Abs             = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesUp_Abs_year        = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesUp_BBEC1           = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesUp_BBEC1_year      = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesUp_EC2             = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesUp_EC2_year        = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesUp_FlavQCD         = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesUp_HF              = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesUp_HF_year         = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesUp_RelBal          = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesUp_RelSample_year  = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesDn  = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesDn_Total           = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesDn_Abs             = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesDn_Abs_year        = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesDn_BBEC1           = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesDn_BBEC1_year      = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesDn_EC2             = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesDn_EC2_year        = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesDn_FlavQCD         = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesDn_HF              = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesDn_HF_year         = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesDn_RelBal          = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesDn_RelSample_year  = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jerUp  = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jerDn  = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSFUp  = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSFDn  = 0;
  Short_t trigWord  = 0;
 
  Float_t ZZMass  = 0;
  Float_t ZZMassErr  = 0;
  Float_t ZZMassErrCorr  = 0;
  Float_t ZZMassPreFSR  = 0;
  Short_t ZZsel  = 0;
  Float_t ZZPt  = 0;
  Float_t ZZEta  = 0;
  Float_t ZZPhi  = 0;
  Float_t ZZjjPt = 0;
 
//   Float_t ZZGoodMass  = 0;
  Bool_t muHLTMatch  = false;
  Bool_t eleHLTMatch  = false;
 
  Int_t CRflag  = 0;
  Float_t Z1Mass  = 0;
  Float_t Z1Pt  = 0;
  Float_t Z1Eta  = 0;
  Float_t Z1Phi  = 0;
  Short_t Z1Flav  = 0;
  Float_t ZZMassRefit  = 0;
  Float_t ZZMassRefitErr  = 0;
  Float_t ZZMassUnrefitErr  = 0;
  Float_t ZZMassCFit  = 0;
  Float_t ZZChi2CFit  = 0;
  Float_t Z2Mass  = 0;
  Float_t Z2Pt  = 0;
  Float_t Z2Eta  = 0;
  Float_t Z2Phi  = 0;
  Short_t Z2Flav  = 0;
  Float_t costhetastar  = 0;
  Float_t helphi  = 0;
  Float_t helcosthetaZ1  = 0;
  Float_t helcosthetaZ2  = 0;
  Float_t phistarZ1  = 0;
  Float_t phistarZ2  = 0;
  Float_t xi  = 0;
  Float_t xistar  = 0;
  Float_t TLE_dR_Z = -1; // Delta-R between a TLE and the Z it does not belong to.
  Float_t TLE_min_dR_3l = 999; // Minimum DR between a TLE and any of the other leptons
  Short_t evtPassMETTrigger = 0;

  Bool_t Z1muHLTMatch1 = false;
  Bool_t Z1muHLTMatch2 = false;
  Bool_t Z1muHLTMatch = false;
  Bool_t Z1eleHLTMatch1 = false;
  Bool_t Z1eleHLTMatch2 = false;
  Bool_t Z1eleHLTMatch = false;
  Bool_t Z2muHLTMatch1 = false;
  Bool_t Z2muHLTMatch2 = false;
  Bool_t Z2muHLTMatch = false;
  Bool_t Z2eleHLTMatch1 = false;
  Bool_t Z2eleHLTMatch2 = false;
  Bool_t Z2eleHLTMatch = false;
  Short_t Z2isGoodTau = 0;
  Short_t Z2isGoodTau_Ele = 0;
  Short_t Z2isGoodTau_Mu = 0;
  Short_t Z2isGoodTau_Tau = 0;


//-------------------------------------------------------------------
//---------------------SV Fit variables------------------------------
//-------------------------------------------------------------------

  Float_t Z2SVMass = -999;
  Float_t Z2SVPt = -999;
  Float_t Z2SVEta = -999;
  Float_t Z2SVPhi = -999;
  Float_t Z2GoodMass = -999;
    
  Float_t Z2Mass_TESup = -999;
  Float_t Z2SVMass_TESup = -999;
  Float_t Z2SVPt_TESup = -999;
  Float_t Z2SVEta_TESup = -999;
  Float_t Z2SVPhi_TESup = -999;
  Float_t Z2GoodMass_TESup = -999;
    
  Float_t Z2Mass_EESup = -999;
  Float_t Z2SVMass_EESup = -999;
  Float_t Z2SVPt_EESup = -999;
  Float_t Z2SVEta_EESup = -999;
  Float_t Z2SVPhi_EESup = -999;
  Float_t Z2GoodMass_EESup = -999;

  Float_t Z2Mass_MESup = -999;
  Float_t Z2SVMass_MESup = -999;
  Float_t Z2SVPt_MESup = -999;
  Float_t Z2SVEta_MESup = -999;
  Float_t Z2SVPhi_MESup = -999;
  Float_t Z2GoodMass_MESup = -999;
    
  Float_t Z2Mass_JESup = -999;
  Float_t Z2SVMass_JESup = -999;
  Float_t Z2SVPt_JESup = -999;
  Float_t Z2SVEta_JESup = -999;
  Float_t Z2SVPhi_JESup = -999;
  Float_t Z2GoodMass_JESup = -999;
    
  Float_t Z2Mass_JERup = -999;
  Float_t Z2SVMass_JERup = -999;
  Float_t Z2SVPt_JERup = -999;
  Float_t Z2SVEta_JERup = -999;
  Float_t Z2SVPhi_JERup = -999;
  Float_t Z2GoodMass_JERup = -999;
      
  Float_t Z2Mass_TESdn = -999;
  Float_t Z2SVMass_TESdn = -999;
  Float_t Z2SVPt_TESdn = -999;
  Float_t Z2SVEta_TESdn = -999;
  Float_t Z2SVPhi_TESdn = -999;
  Float_t Z2GoodMass_TESdn = -999;
    
  Float_t Z2Mass_EESdn = -999;
  Float_t Z2SVMass_EESdn = -999;
  Float_t Z2SVPt_EESdn = -999;
  Float_t Z2SVEta_EESdn = -999;
  Float_t Z2SVPhi_EESdn = -999;
  Float_t Z2GoodMass_EESdn = -999;

  Float_t Z2Mass_MESdn = -999;
  Float_t Z2SVMass_MESdn = -999;
  Float_t Z2SVPt_MESdn = -999;
  Float_t Z2SVEta_MESdn = -999;
  Float_t Z2SVPhi_MESdn = -999;
  Float_t Z2GoodMass_MESdn = -999;
    
  Float_t Z2Mass_JESdn = -999;
  Float_t Z2SVMass_JESdn = -999;
  Float_t Z2SVPt_JESdn = -999;
  Float_t Z2SVEta_JESdn = -999;
  Float_t Z2SVPhi_JESdn = -999;
  Float_t Z2GoodMass_JESdn = -999;
    
  Float_t Z2Mass_JERdn = -999;
  Float_t Z2SVMass_JERdn = -999;
  Float_t Z2SVPt_JERdn = -999;
  Float_t Z2SVEta_JERdn = -999;
  Float_t Z2SVPhi_JERdn = -999;
  Float_t Z2GoodMass_JERdn = -999;
  //------------------------------End SV fit variables--------------------------------
    
  //MET
  Float_t MET = 0;
  Float_t METPhi = 0;

  std::vector<float> LepPt;
  std::vector<float> LepEta;
  std::vector<float> LepPhi;
  std::vector<float> LepM;
  std::vector<float> LepSCEta;
  std::vector<short> LepLepId;
  std::vector<float> LepSIP;
  std::vector<float> Lepdxy;
  std::vector<float> Lepdz;
  std::vector<float> LepTime;
  std::vector<bool> LepisID;
  std::vector<float> LepBDT;
  std::vector<bool> LepisCrack;
  std::vector<char> LepMissingHit;
  std::vector<float> LepChargedHadIso;
  std::vector<float> LepNeutralHadIso;
  std::vector<float> LepPhotonIso;
  std::vector<float> LepPUIsoComponent;
  std::vector<float> LepCombRelIsoPF;
  std::vector<short> LepisLoose;
  std::vector<float> LepSF;
  std::vector<float> LepSF_Unc;
  std::vector<float> LepScale_Total_Up;
  std::vector<float> LepScale_Total_Dn;
  std::vector<float> LepScale_Stat_Up;
  std::vector<float> LepScale_Stat_Dn;
  std::vector<float> LepScale_Syst_Up;
  std::vector<float> LepScale_Syst_Dn;
  std::vector<float> LepScale_Gain_Up;
  std::vector<float> LepScale_Gain_Dn;
  std::vector<float> LepSigma_Total_Up;
  std::vector<float> LepSigma_Total_Dn;
  std::vector<float> LepSigma_Rho_Up;
  std::vector<float> LepSigma_Rho_Dn;
  std::vector<float> LepSigma_Phi_Up;
  std::vector<float> LepSigma_Phi_Dn;

//tau specified
  std::vector<short> TauVSmu;
  std::vector<short> TauVSe;
  std::vector<short> TauVSjet;
  std::vector<float> TauDecayMode;
  std::vector<short> TauGenMatch;
  std::vector<float> TauTES_p_Up;
  std::vector<float> TauTES_p_Dn;
  std::vector<float> TauTES_m_Up;
  std::vector<float> TauTES_m_Dn;
  std::vector<float> TauTES_e_Up;
  std::vector<float> TauTES_e_Dn;
  std::vector<float> TauFES_p_Up;
  std::vector<float> TauFES_p_Dn;
  std::vector<float> TauFES_m_Up;
  std::vector<float> TauFES_m_Dn;
  std::vector<float> TauFES_e_Up;
  std::vector<float> TauFES_e_Dn;

//HLT match
  std::vector<short> HLTMatch1;
//  std::vector<short> HLTMatch2;

  std::vector<float> fsrPt;
  std::vector<float> fsrEta;
  std::vector<float> fsrPhi;
  std::vector<float> fsrDR;
  std::vector<short> fsrLept;
  std::vector<short> fsrLeptID;
  std::vector<float> fsrGenPt;
  Bool_t passIsoPreFSR = 0;

  std::vector<float> JetPt ;
  std::vector<float> JetEta ;
  std::vector<float> JetPhi ;
  std::vector<float> JetMass ;
  std::vector<float> JetEnergy ;
  std::vector<float> JetBTagger ;
  std::vector<float> JetIsBtagged;
  std::vector<float> JetIsBtaggedWithSF;
  std::vector<float> JetIsBtaggedWithSFUp;
  std::vector<float> JetIsBtaggedWithSFDn;
  std::vector<float> JetQGLikelihood;
  std::vector<float> JetAxis2;
  std::vector<float> JetMult;
  std::vector<float> JetPtD;
  std::vector<float> JetSigma ;
  std::vector<float> JetSigma_Total ;
  std::vector<float> JetSigma_Abs ;
  std::vector<float> JetSigma_Abs_year ;
  std::vector<float> JetSigma_BBEC1 ;
  std::vector<float> JetSigma_BBEC1_year ;
  std::vector<float> JetSigma_EC2 ;
  std::vector<float> JetSigma_EC2_year ;
  std::vector<float> JetSigma_FlavQCD ;
  std::vector<float> JetSigma_HF ;
  std::vector<float> JetSigma_HF_year ;
  std::vector<float> JetSigma_RelBal ;
  std::vector<float> JetSigma_RelSample_year ;
  std::vector<short> JetHadronFlavour;
  std::vector<short> JetPartonFlavour;
   
  std::vector<float> JetPtJEC_noJER;
  std::vector<float> JetRawPt;

  std::vector<float> JetPUValue;
  std::vector<short> JetPUID;
  std::vector<float> JetPUID_score;
    
  std::vector<short> JetID;
   
  std::vector<float> JetJESUp ;
  std::vector<float> JetJESUp_Total ;
  std::vector<float> JetJESUp_Abs ;
  std::vector<float> JetJESUp_Abs_year ;
  std::vector<float> JetJESUp_BBEC1 ;
  std::vector<float> JetJESUp_BBEC1_year ;
  std::vector<float> JetJESUp_EC2 ;
  std::vector<float> JetJESUp_EC2_year ;
  std::vector<float> JetJESUp_FlavQCD ;
  std::vector<float> JetJESUp_HF ;
  std::vector<float> JetJESUp_HF_year ;
  std::vector<float> JetJESUp_RelBal ;
  std::vector<float> JetJESUp_RelSample_year ;
  std::vector<float> JetJESDown ;
  std::vector<float> JetJESDown_Total ;
  std::vector<float> JetJESDown_Abs ;
  std::vector<float> JetJESDown_Abs_year ;
  std::vector<float> JetJESDown_BBEC1 ;
  std::vector<float> JetJESDown_BBEC1_year ;
  std::vector<float> JetJESDown_EC2 ;
  std::vector<float> JetJESDown_EC2_year ;
  std::vector<float> JetJESDown_FlavQCD ;
  std::vector<float> JetJESDown_HF ;
  std::vector<float> JetJESDown_HF_year ;
  std::vector<float> JetJESDown_RelBal ;
  std::vector<float> JetJESDown_RelSample_year ;
   
  std::vector<float> JetJERUp ;
  std::vector<float> JetJERDown ;

  Float_t DiJetMass  = -99;
  Float_t DiJetDEta  = -99;

  Short_t nExtraLep  = 0;
  Short_t nExtraZ  = 0;
  std::vector<float> ExtraLepPt;
  std::vector<float> ExtraLepEta;
  std::vector<float> ExtraLepPhi ;
  std::vector<short> ExtraLepLepId;
  
  // Photon info
  std::vector<float> PhotonPt ;
  std::vector<float> PhotonEta ;
  std::vector<float> PhotonPhi ;
  std::vector<bool> PhotonIsCutBasedLooseID;
   
   
  Short_t genFinalState  = 0;
  Int_t genProcessId  = 0;
  Float_t genHEPMCweight  = 0;
  Float_t genHEPMCweight_NNLO  = 0;

  Float_t PythiaWeight_isr_muRoneoversqrt2 = 0;
  Float_t PythiaWeight_fsr_muRoneoversqrt2 = 0;
  Float_t PythiaWeight_isr_muRsqrt2 = 0;
  Float_t PythiaWeight_fsr_muRsqrt2 = 0;
  Float_t PythiaWeight_isr_muR0p5 = 0;
  Float_t PythiaWeight_fsr_muR0p5 = 0;
  Float_t PythiaWeight_isr_muR2 = 0;
  Float_t PythiaWeight_fsr_muR2 = 0;
  Float_t PythiaWeight_isr_muR0p25 = 0;
  Float_t PythiaWeight_fsr_muR0p25 = 0;
  Float_t PythiaWeight_isr_muR4 = 0;
  Float_t PythiaWeight_fsr_muR4 = 0;


  Short_t genExtInfo  = 0;
  Float_t xsection  = 0;
  Float_t genxsection = 0;
  Float_t genbranchingratio = 0;
  Float_t dataMCWeight  = 0;
  Float_t trigEffWeight  = 0;
  Float_t HqTMCweight  = 0;
  Float_t ZXFakeweight  = 0;
  Float_t overallEventWeight  = 0;
  Float_t L1prefiringWeight = 0;
  Float_t L1prefiringWeightUp = 0;
  Float_t L1prefiringWeightDn = 0;
  Float_t L1prefiringWeight_ECAL = 0;
  Float_t L1prefiringWeightUp_ECAL = 0;
  Float_t L1prefiringWeightDn_ECAL = 0;
  Float_t L1prefiringWeight_Mu = 0;
  Float_t L1prefiringWeightUp_Mu = 0;
  Float_t L1prefiringWeightDn_Mu = 0;
  Float_t GenHMass  = 0;
  Float_t GenHPt  = 0;
  Float_t GenHRapidity  = 0;
  Float_t GenZ1Mass  = 0;
  Float_t GenZ1Eta  = 0;
  Float_t GenZ1Pt  = 0;
  Float_t GenZ1Phi  = 0;
  Float_t GenZ1Flav  = 0;
  Float_t GenZ2Mass  = 0;
  Float_t GenZ2Eta  = 0;
  Float_t GenZ2Pt  = 0;
  Float_t GenZ2Phi  = 0;
  Float_t GenZ2Flav  = 0;
  Float_t GenLep1Pt  = 0;
  Float_t GenLep1Eta  = 0;
  Float_t GenLep1Phi  = 0;
  Float_t GenLep1M  = 0;
  Short_t GenLep1Id  = 0;
  Float_t GenLep2Pt  = 0;
  Float_t GenLep2Eta  = 0;
  Float_t GenLep2Phi  = 0;
  Float_t GenLep2M  = 0;
  Short_t GenLep2Id  = 0;
  Float_t GenLep3Pt  = 0;
  Float_t GenLep3Eta  = 0;
  Float_t GenLep3Phi  = 0;
  Float_t GenLep3M  = 0;
  Short_t GenLep3Id  = 0;
  Float_t GenLep4Pt  = 0;
  Float_t GenLep4Eta  = 0;
  Float_t GenLep4Phi  = 0;
  Float_t GenLep4M  = 0;
  Short_t GenLep4Id  = 0;
  //Visible information
  Float_t GenVisZ1Mass  = 0;
  Float_t GenVisZ1Eta  = 0;
  Float_t GenVisZ1Pt  = 0;
  Float_t GenVisZ1Phi  = 0;
  Float_t GenVisZ1Flav  = 0;
  Float_t GenVisZ2Mass  = 0;
  Float_t GenVisZ2Eta  = 0;
  Float_t GenVisZ2Pt  = 0;
  Float_t GenVisZ2Phi  = 0;
  Float_t GenVisZ2Flav  = 0;
  Float_t GenVisLep1Pt  = 0;
  Float_t GenVisLep1Eta  = 0;
  Float_t GenVisLep1Phi  = 0;
  Float_t GenVisLep1M  = 0;
  Short_t GenVisLep1Id  = 0;
  Float_t GenVisLep2Pt  = 0;
  Float_t GenVisLep2Eta  = 0;
  Float_t GenVisLep2Phi  = 0;
  Float_t GenVisLep2M  = 0;
  Short_t GenVisLep2Id  = 0;
  Float_t GenVisLep3Pt  = 0;
  Float_t GenVisLep3Eta  = 0;
  Float_t GenVisLep3Phi  = 0;
  Float_t GenVisLep3M  = 0;
  Short_t GenVisLep3Id  = 0;
  Float_t GenVisLep4Pt  = 0;
  Float_t GenVisLep4Eta  = 0;
  Float_t GenVisLep4Phi  = 0;
  Float_t GenVisLep4M  = 0;
  Short_t GenVisLep4Id  = 0;

  Float_t GenAssocLep1Pt  = 0;
  Float_t GenAssocLep1Eta  = 0;
  Float_t GenAssocLep1Phi  = 0;
  Short_t GenAssocLep1Id  = 0;
  Float_t GenAssocLep2Pt  = 0;
  Float_t GenAssocLep2Eta  = 0;
  Float_t GenAssocLep2Phi  = 0;
  Short_t GenAssocLep2Id  = 0;

  Float_t GenLep1Iso  = 0; //AT
  Float_t GenLep2Iso  = 0; //AT
  Float_t GenLep3Iso  = 0; //AT
  Float_t GenLep4Iso  = 0; //AT
  Float_t Gencosthetastar  = 0; //AT
  Float_t GenhelcosthetaZ1  = 0; //AT
  Float_t GenhelcosthetaZ2  = 0; //AT
  Float_t Genhelphi  = 0; //AT
  Float_t GenphistarZ1 = 0; //AT
  std::vector<float> GenJetPt; //ATjets
  std::vector<float> GenJetMass; //ATjets
  std::vector<float> GenJetEta; //ATjets
  std::vector<float> GenJetPhi; //ATjets
  std::vector<float> GenJetRapidity; //ATjets
  Int_t nGenJet = 0; //ATjets
  std::vector<float> GenCleanedJetPt; //ATjets
  std::vector<float> GenCleanedJetMass; //ATjets
  std::vector<float> GenCleanedJetEta; //ATjets
  std::vector<float> GenCleanedJetPhi; //ATjets
  std::vector<float> GenCleanedJetRapidity; //ATjets
  std::vector<float> GenCleanedJetHadronFlavour; //ATjets
  Int_t nCleanedGenJet = 0; //ATjets

  Int_t   htxsNJets = -1;
  Float_t htxsHPt = 0;
  Int_t   htxs_errorCode=-1;
  Int_t   htxs_prodMode=-1;
  Int_t   htxs_stage0_cat = -1;
  Int_t   htxs_stage1p1_cat = -1;
  Int_t   htxs_stage1p0_cat = -1;
  Int_t   htxs_stage1p2_cat = -1;
  Float_t ggH_NNLOPS_weight = 0;
  Float_t ggH_NNLOPS_weight_unc = 0;
  std::vector<float> qcd_ggF_uncertSF;
    
  //ATbbf
  //Event variables
  // Short_t GENfinalState;
  bool passedFiducialSelection_bbf;
  // lepton variables
  std::vector<float> GENlep_pt; std::vector<float> GENlep_eta; std::vector<float> GENlep_phi; std::vector<float> GENlep_mass;
  std::vector<TLorentzVector> GENlep_tetra;
  std::vector<short> GENlep_id; std::vector<short> GENlep_status;
  std::vector<short> GENlep_MomId; std::vector<short> GENlep_MomMomId;
  // Short_t GENlep_Hindex[4];//position of Higgs candidate leptons in lep_p4: 0 = Z1 lead, 1 = Z1 sub, 2 = Z2 lead, 3 = Z3 sub
  std::vector<Short_t> GENlep_Hindex;
  // std::vector<float> GENlep_isoCH; std::vector<float> GENlep_isoNH; std::vector<float> GENlep_isoPhot;
  std::vector<float> GENlep_RelIso;
  // Higgs candidate variables (calculated using selected gen leptons)
  std::vector<float> GENH_pt; std::vector<float> GENH_eta; std::vector<float> GENH_phi; std::vector<float> GENH_mass;
  float_t GENmass4l, GENpT4l, GENeta4l, GENphi4l, GENrapidity4l;
  // GENmass4e, GENmass4mu, GENmass2e2mu, GENpT4l, GENeta4l, GENrapidity4l;
  // float_t GENMH; //mass directly from gen particle with id==25
  float_t GENcosTheta1, GENcosTheta2, GENcosThetaStar, GENPhi, GENPhi1;
  // Z candidate variables
  std::vector<float> GENZ_pt; std::vector<float> GENZ_eta; std::vector<float> GENZ_phi; std::vector<float> GENZ_mass;
  std::vector<short> GENZ_DaughtersId; std::vector<short> GENZ_MomId;
  float_t  GENmassZ1, GENmassZ2;
  // GENpTZ1, GENpTZ2,
  float_t GENdPhiZZ, GENmassZZ, GENpTZZ;
  // Jets
  // std::vector<float> GENjet_pt; std::vector<float> GENjet_eta; std::vector<float> GENjet_phi; std::vector<float> GENjet_mass;
  std::vector<float> GENjetsPt_pt30_eta4p7; std::vector<float> GENjetsEta_pt30_eta4p7; std::vector<float> GENjetsPhi_pt30_eta4p7; std::vector<float> GENjetsMass_pt30_eta4p7;
  std::vector<float> GENjetsPt_pt30_eta2p5; std::vector<float> GENjetsEta_pt30_eta2p5; std::vector<float> GENjetsPhi_pt30_eta2p5; std::vector<float> GENjetsMass_pt30_eta2p5;

  Short_t GENnjets_pt30_eta4p7;
  Short_t GENnjets_pt30_eta2p5;


//FIXME: temporary fix to the mismatch of charge() and sign(pdgId()) for muons with BTT=4
  int getPdgId(const reco::Candidate* p) {
    int id = p->pdgId();
    if (id!=22 && //for TLEs
	signbit(id) && p->charge()<0) id*=-1; // negative pdgId must be positive charge
    return id;
  }

}

using namespace std;
using namespace edm;
//
// class declaration
//
class HZZ4lNtupleMaker : public edm::EDAnalyzer {
public:
  explicit HZZ4lNtupleMaker(const edm::ParameterSet&);
  ~HZZ4lNtupleMaker();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  static float EvalSpline(TSpline3* const& sp, float xval);

  static void addweight(float &weight, float weighttoadd);

private:
  virtual void beginJob() ;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  void BookAllBranches();
  virtual void FillKFactors(edm::Handle<GenEventInfoProduct>& genInfo, std::vector<const reco::Candidate *>& genZLeps);
  virtual void FillCandidate(const pat::CompositeCandidate& higgs, bool evtPass, const edm::Event&, const Int_t CRflag);
  virtual void FillJet(const pat::Jet& jet);
  virtual void FillPhoton(int year, const pat::Photon& photon);
  virtual void endJob() ;
    
  virtual void FillGENCandidate(const pat::CompositeCandidate& cand); //AT

  void FillHGenInfo(const math::XYZTLorentzVector Hp, float w);
  void FillZGenInfo(Short_t Z1Id, Short_t Z2Id,
                    const math::XYZTLorentzVector pZ1, const math::XYZTLorentzVector pZ2);
  void FillLepGenInfo(Short_t Lep1Id, Short_t Lep2Id, Short_t Lep3Id, Short_t Lep4Id,
    const math::XYZTLorentzVector Lep1, const math::XYZTLorentzVector Lep2, const math::XYZTLorentzVector Lep3, const math::XYZTLorentzVector Lep4);
  void FillVisZGenInfo(Short_t Z1Id, Short_t Z2Id,
                    const math::XYZTLorentzVector pZ1, const math::XYZTLorentzVector pZ2);
  void FillVisLepGenInfo(Short_t Lep1Id, Short_t Lep2Id, Short_t Lep3Id, Short_t Lep4Id,
    const math::XYZTLorentzVector Lep1, const math::XYZTLorentzVector Lep2, const math::XYZTLorentzVector Lep3, const math::XYZTLorentzVector Lep4);
  void FillAssocLepGenInfo(std::vector<const reco::Candidate *>& AssocLeps);
  void FillLepGenIso(float_t Lep1Iso, float_t Lep2Iso, float_t Lep3Iso, float_t Lep4Iso); //AT
  // void FillJetGen(float_t Lep1Iso, float_t Lep2Iso, float_t Lep3Iso, float_t Lep4Iso); //AT
  
  Float_t getAllWeight(const vector<const reco::Candidate*>& leptons);
  Float_t getHqTWeight(double mH, double genPt) const;
  Float_t getFakeWeight(Float_t LepPt, Float_t LepEta, Int_t LepID, Int_t LepZ1ID);
  Int_t FindBinValue(TGraphErrors *tgraph, double value);

  void getCheckedUserFloat(const pat::CompositeCandidate& cand, const std::string& strval, Float_t& setval, Float_t defaultval=0);
	
	

  // ----------member data ---------------------------
  ZZ4lConfigHelper myHelper;
  int theChannel;
  std::string theCandLabel;
  TString theFileName;

  HZZ4lNtupleFactory *myTree;
  TH1F *hCounter;

  bool isMC;
  bool preVFP=false;
  bool is_loose_ele_selection; // Collection includes candidates with loose electrons/TLEs
  bool applySkim;       //   "     "      "         skim (if skipEmptyEvents=true)
  bool skipEmptyEvents; // Skip events whith no selected candidate (otherwise, gen info is preserved for all events; candidates not passing trigger&&skim are flagged with negative ZZsel)
  FailedTreeLevel failedTreeLevel;  //if/how events with no selected candidate are written to a separate tree (see miscenums.h for details)
  edm::InputTag metTag;
  bool applyTrigger;    // Keep only events passing trigger (overriden if skipEmptyEvents=False)
  bool applyTrigEffWeight;// apply trigger efficiency weight (concerns samples where trigger is not applied)
  Float_t xsec;
  Float_t genxsec;
  Float_t genbr;
  int year;
  double sqrts;
  double Hmass;

  int apply_K_NNLOQCD_ZZGG; // 0: Do not; 1: NNLO/LO; 2: NNLO/NLO; 3: NLO/LO
  bool apply_K_NNLOQCD_ZZQQB;
  bool apply_K_NLOEW_ZZQQB;
  bool apply_QCD_GGF_UNCERT;

  edm::EDGetTokenT<edm::View<reco::Candidate> > genParticleToken;
  edm::Handle<edm::View<reco::Candidate> > genParticles;
    
  edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_bbf;//ATbbf
  edm::Handle<reco::GenParticleCollection> genParticles_bbf;//ATbbf
  edm::Handle<edm::View<pat::PackedGenParticle> > packedgenParticles; //ATbbf
  edm::Handle<edm::View<reco::GenJet> > genJets; //ATjets
  edm::EDGetTokenT<edm::View<reco::GenJet> > genJetsToken; //ATjets
  edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedgenParticlesToken; //ATbbf
    
  edm::EDGetTokenT<GenEventInfoProduct> genInfoToken;
  edm::EDGetTokenT<edm::View<pat::CompositeCandidate> > candToken;
  edm::EDGetTokenT<edm::TriggerResults> triggerResultToken;
  edm::EDGetTokenT<vector<reco::Vertex> > vtxToken;
  edm::EDGetTokenT<edm::View<pat::Jet> > jetToken;
  edm::EDGetTokenT<pat::PhotonCollection> photonToken;
  
  edm::EDGetTokenT<pat::METCollection> metToken;
  edm::EDGetTokenT<math::Error<2>::type> theCovTag;
  edm::EDGetTokenT<double> theMETdxUPTESTag;
  edm::EDGetTokenT<double> theMETdyUPTESTag;
  edm::EDGetTokenT<double> theMETdxDOWNTESTag;
  edm::EDGetTokenT<double> theMETdyDOWNTESTag;
  edm::EDGetTokenT<double> theMETdxUPEESTag;
  edm::EDGetTokenT<double> theMETdyUPEESTag;
  edm::EDGetTokenT<double> theMETdxDOWNEESTag;
  edm::EDGetTokenT<double> theMETdyDOWNEESTag;
  edm::EDGetTokenT<double> theMETdxUPMESTag;
  edm::EDGetTokenT<double> theMETdyUPMESTag;
  edm::EDGetTokenT<double> theMETdxDOWNMESTag;
  edm::EDGetTokenT<double> theMETdyDOWNMESTag;
  edm::EDGetTokenT<double> theMETdxUPJESTag;
  edm::EDGetTokenT<double> theMETdyUPJESTag;
  edm::EDGetTokenT<double> theMETdxDOWNJESTag;
  edm::EDGetTokenT<double> theMETdyDOWNJESTag;
  edm::EDGetTokenT<double> theMETdxUPJERTag;
  edm::EDGetTokenT<double> theMETdyUPJERTag;
  edm::EDGetTokenT<double> theMETdxDOWNJERTag;
  edm::EDGetTokenT<double> theMETdyDOWNJERTag;
    
  edm::EDGetTokenT<pat::MuonCollection> muonToken;
  edm::EDGetTokenT<pat::ElectronCollection> electronToken;
  edm::EDGetTokenT<HTXS::HiggsClassification> htxsToken;
  edm::EDGetTokenT<edm::MergeableCounter> preSkimToken;
   
  edm::EDGetTokenT<edm::View<pat::CompositeCandidate>> GENCandidatesToken; //ATMELA
  edm::Handle<edm::View<pat::CompositeCandidate> > GENCandidates; //ATMELA
    
  edm::EDGetTokenT< double > prefweight_token;
  edm::EDGetTokenT< double > prefweightup_token;
  edm::EDGetTokenT< double > prefweightdown_token;
    
  edm::EDGetTokenT< double > prefweightECAL_token;
  edm::EDGetTokenT< double > prefweightupECAL_token;
  edm::EDGetTokenT< double > prefweightdownECAL_token;

  edm::EDGetTokenT< double > prefweightMuon_token;
  edm::EDGetTokenT< double > prefweightupMuon_token;
  edm::EDGetTokenT< double > prefweightdownMuon_token;

  PileUpWeight* pileUpReweight;

  //counters
  Float_t Nevt_Gen;
  Float_t Nevt_Gen_lumiBlock;

  Float_t gen_ZZ4mu;
  Float_t gen_ZZ4e;
  Float_t gen_ZZ2mu2e;
  Float_t gen_ZZ2l2tau;
  Float_t gen_ZZ2emu2tau;
  Float_t gen_ZZ2mu2tau;
  Float_t gen_ZZ2e2tau;
  Float_t gen_ZZ4tau;
  Float_t gen_ZZ4mu_EtaAcceptance;
  Float_t gen_ZZ4mu_LeptonAcceptance;
  Float_t gen_ZZ4e_EtaAcceptance;
  Float_t gen_ZZ4e_LeptonAcceptance;
  Float_t gen_ZZ2mu2e_EtaAcceptance;
  Float_t gen_ZZ2mu2e_LeptonAcceptance;
  Float_t gen_BUGGY;
  Float_t gen_Unknown;

  Float_t gen_sumPUWeight;
  Float_t gen_sumGenMCWeight;
  Float_t gen_sumWeights;

  string sampleName;
    
  string dataTag;

  std::vector<const reco::Candidate *> genFSR;

  std::vector<std::vector<float> > ewkTable;
  TSpline3* spkfactor_ggzz_nnlo[9]; // Nominal, PDFScaleDn, PDFScaleUp, QCDScaleDn, QCDScaleUp, AsDn, AsUp, PDFReplicaDn, PDFReplicaUp
  TSpline3* spkfactor_ggzz_nlo[9]; // Nominal, PDFScaleDn, PDFScaleUp, QCDScaleDn, QCDScaleUp, AsDn, AsUp, PDFReplicaDn, PDFReplicaUp

  LeptonSFHelper *lepSFHelper;

  TH2D* h_weight; //HqT weights
  //TH2F *h_ZXWeightMuo;
  //TH2F *h_ZXWeightEle;
  TH2D* h_ZXWeight[4];
	
  TGraphErrors *gr_NNLOPSratio_pt_powheg_0jet;
  TGraphErrors *gr_NNLOPSratio_pt_powheg_1jet;
  TGraphErrors *gr_NNLOPSratio_pt_powheg_2jet;
  TGraphErrors *gr_NNLOPSratio_pt_powheg_3jet;

  //tau discriminants
  std::vector<std::string> VSmuDisc;
  std::vector<std::string> VSeDisc;
  std::vector<std::string> VSjetDisc;
    
  bool firstRun;

};

//
// constructors and destructor
//
HZZ4lNtupleMaker::HZZ4lNtupleMaker(const edm::ParameterSet& pset) :
  myHelper(pset),
  theChannel(myHelper.channel()), // Valid options: ZZ, ZLL, ZL
  theCandLabel(pset.getUntrackedParameter<string>("CandCollection")), // Name of input ZZ collection
  theFileName(pset.getUntrackedParameter<string>("fileName")),
  myTree(nullptr),
  skipEmptyEvents(pset.getParameter<bool>("skipEmptyEvents")), // Do not store events with no selected candidate (normally: true)
  failedTreeLevel(FailedTreeLevel(pset.getParameter<int>("failedTreeLevel"))),

  metTag(pset.getParameter<edm::InputTag>("metSrc")),
  
  applyTrigger(pset.getParameter<bool>("applyTrigger")), // Reject events that do not pass trigger (normally: true)
  applyTrigEffWeight(pset.getParameter<bool>("applyTrigEff")), //Apply an additional efficiency weights for MC samples where triggers are not present (normally: false)
  xsec(pset.getParameter<double>("xsec")),
  genxsec(pset.getParameter<double>("GenXSEC")),
  genbr(pset.getParameter<double>("GenBR")),
  year(pset.getParameter<int>("setup")),
  sqrts(SetupToSqrts(year)),
  //recoMElist(pset.getParameter<std::vector<std::string>>("recoProbabilities")),

  //lheHandler(nullptr),
  apply_K_NNLOQCD_ZZGG(pset.getParameter<int>("Apply_K_NNLOQCD_ZZGG")),
  apply_K_NNLOQCD_ZZQQB(pset.getParameter<bool>("Apply_K_NNLOQCD_ZZQQB")),
  apply_K_NLOEW_ZZQQB(pset.getParameter<bool>("Apply_K_NLOEW_ZZQQB")),
  apply_QCD_GGF_UNCERT(pset.getParameter<bool>("Apply_QCD_GGF_UNCERT")),

  pileUpReweight(nullptr),
  sampleName(pset.getParameter<string>("sampleName")),
  dataTag(pset.getParameter<string>("dataTag")),
  h_weight(0),

  firstRun(true)
{
  //cout<< "Beginning Constructor\n\n\n" <<endl;
  consumesMany<std::vector< PileupSummaryInfo > >();
  genParticleToken = consumes<edm::View<reco::Candidate> >(edm::InputTag("prunedGenParticles"));
  genParticleToken_bbf = consumes<reco::GenParticleCollection>(edm::InputTag("prunedGenParticles"));
  packedgenParticlesToken = consumes<edm::View<pat::PackedGenParticle> > (edm::InputTag("packedGenParticles")); //ATbbf
  genInfoToken = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
  genJetsToken = consumes<edm::View<reco::GenJet> >(edm::InputTag("slimmedGenJets")); //ATjets
  GENCandidatesToken = consumes<edm::View<pat::CompositeCandidate> >(edm::InputTag("GENLevel"));
  consumesMany<LHEEventProduct>();
  candToken = consumes<edm::View<pat::CompositeCandidate> >(edm::InputTag(theCandLabel));

  is_loose_ele_selection = false;
  if(pset.exists("is_loose_ele_selection")) { 
    is_loose_ele_selection = pset.getParameter<bool>("is_loose_ele_selection");
  }
  triggerResultToken = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults"));
  vtxToken = consumes<vector<reco::Vertex> >(edm::InputTag("goodPrimaryVertices"));
  jetToken = consumes<edm::View<pat::Jet> >(edm::InputTag("dressedJets"));
  photonToken = consumes<pat::PhotonCollection>(edm::InputTag("slimmedPhotons"));
      
  metToken = consumes<pat::METCollection>(pset.getParameter<edm::InputTag>("metSrc"));
  theCovTag = consumes<math::Error<2>::type>(pset.getParameter<edm::InputTag>("covSrc"));
  theMETdxUPTESTag = consumes<double>(pset.getParameter<edm::InputTag>("METdxUPTES"));
  theMETdyUPTESTag = consumes<double>(pset.getParameter<edm::InputTag>("METdyUPTES"));
  theMETdxDOWNTESTag = consumes<double>(pset.getParameter<edm::InputTag>("METdxDOWNTES"));
  theMETdyDOWNTESTag = consumes<double>(pset.getParameter<edm::InputTag>("METdyDOWNTES"));
  theMETdxUPEESTag = consumes<double>(pset.getParameter<edm::InputTag>("METdxUPEES"));
  theMETdyUPEESTag = consumes<double>(pset.getParameter<edm::InputTag>("METdyUPEES"));
  theMETdxDOWNEESTag = consumes<double>(pset.getParameter<edm::InputTag>("METdxDOWNEES"));
  theMETdyDOWNEESTag = consumes<double>(pset.getParameter<edm::InputTag>("METdyDOWNEES"));
  theMETdxUPMESTag = consumes<double>(pset.getParameter<edm::InputTag>("METdxUPMES"));
  theMETdyUPMESTag = consumes<double>(pset.getParameter<edm::InputTag>("METdyUPMES"));
  theMETdxDOWNMESTag = consumes<double>(pset.getParameter<edm::InputTag>("METdxDOWNMES"));
  theMETdyDOWNMESTag = consumes<double>(pset.getParameter<edm::InputTag>("METdyDOWNMES"));
  theMETdxUPJESTag = consumes<double>(pset.getParameter<edm::InputTag>("METdxUPJES"));
  theMETdyUPJESTag = consumes<double>(pset.getParameter<edm::InputTag>("METdyUPJES"));
  theMETdxDOWNJESTag = consumes<double>(pset.getParameter<edm::InputTag>("METdxDOWNJES"));
  theMETdyDOWNJESTag = consumes<double>(pset.getParameter<edm::InputTag>("METdyDOWNJES"));
  theMETdxUPJERTag = consumes<double>(pset.getParameter<edm::InputTag>("METdxUPJER"));
  theMETdyUPJERTag = consumes<double>(pset.getParameter<edm::InputTag>("METdyUPJER"));
  theMETdxDOWNJERTag = consumes<double>(pset.getParameter<edm::InputTag>("METdxDOWNJER"));
  theMETdyDOWNJERTag = consumes<double>(pset.getParameter<edm::InputTag>("METdyDOWNJER"));
      
  //metNoHFToken = consumes<pat::METCollection>(edm::InputTag("slimmedMETsNoHF"));
  muonToken = consumes<pat::MuonCollection>(edm::InputTag("slimmedMuons"));
  electronToken = consumes<pat::ElectronCollection>(edm::InputTag("slimmedElectrons"));
  preSkimToken = consumes<edm::MergeableCounter,edm::InLumi>(edm::InputTag("preSkimCounter"));
   
  if (skipEmptyEvents) {
    applySkim=true;
  } else {
    applyTrigger=false; // This overrides the card applyTrigger
    applySkim=false;
    failedTreeLevel=noFailedTree; // This overrides the card failedTreeLevel
  }

  if (applyTrigEffWeight&&applyTrigger) {
    cout << "ERROR: cannot have applyTrigEffWeight == applyTrigger == true" << endl;
  }

  isMC = myHelper.isMC();
   
  if( isMC )
  {
     prefweight_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProb"));
     prefweightup_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbUp"));
     prefweightdown_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbDown"));
      
     // CMSSW_10_6_26 stores ECAL and muon weights separately
     // Access them here but not used in the analysis since prefiringweight = prefiringweightECAL*prefiringweightMuon
     prefweightECAL_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbECAL"));
     prefweightupECAL_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbECALUp"));
     prefweightdownECAL_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbECALDown"));

     prefweightMuon_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbMuon"));
     prefweightupMuon_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbMuonUp"));
     prefweightdownMuon_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbMuonDown"));
  }
   
   
  if (isMC){
    //    lheHandler = new LHEHandler(
    //      ((MELAEvent::CandidateVVMode)(pset.getParameter<int>("VVMode")+1)), // FIXME: Need to pass strings and interpret them instead!
    //      pset.getParameter<int>("VVDecayMode"),
    //      (addLHEKinematics ? LHEHandler::doHiggsKinematics : LHEHandler::noKinematics),
    //      year, LHEHandler::tryNNPDF30, LHEHandler::tryNLO, LHEHandler::CMS_Run2_UL
    //    );
    htxsToken = consumes<HTXS::HiggsClassification>(edm::InputTag("rivetProducerHTXS","HiggsClassification"));
    pileUpReweight = new PileUpWeight(myHelper.sampleType(), myHelper.setup());
  }

  Nevt_Gen = 0;
  Nevt_Gen_lumiBlock = 0;

  //For Efficiency studies
  gen_ZZ4mu = 0;
  gen_ZZ4e = 0;
  gen_ZZ2mu2e = 0;
  gen_ZZ2l2tau = 0;
  gen_ZZ2emu2tau = 0;
  gen_ZZ2mu2tau = 0;
  gen_ZZ2e2tau = 0;
  gen_ZZ4tau = 0;
  gen_ZZ4mu_EtaAcceptance = 0;
  gen_ZZ4mu_LeptonAcceptance = 0;
  gen_ZZ4e_EtaAcceptance = 0;
  gen_ZZ4e_LeptonAcceptance = 0;
  gen_ZZ2mu2e_EtaAcceptance = 0;
  gen_ZZ2mu2e_LeptonAcceptance = 0;
  gen_BUGGY = 0;
  gen_Unknown = 0;

  gen_sumPUWeight = 0.f;
  gen_sumGenMCWeight = 0.f;
  gen_sumWeights =0.f;

  std::string fipPath;

  // Read EWK K-factor table from file
  edm::FileInPath ewkFIP("HTauTauHMuMu/AnalysisStep/data/kfactors/ZZ_EwkCorrections.dat");
  fipPath=ewkFIP.fullPath();
  ewkTable = EwkCorrections::readFile_and_loadEwkTable(fipPath.data());

  // Read the ggZZ k-factor shape from file
  TString strZZGGKFVar[9]={
    "Nominal", "PDFScaleDn", "PDFScaleUp", "QCDScaleDn", "QCDScaleUp", "AsDn", "AsUp", "PDFReplicaDn", "PDFReplicaUp"
  };
  edm::FileInPath ggzzFIP_NNLO("HTauTauHMuMu/AnalysisStep/data/kfactors/Kfactor_Collected_ggHZZ_2l2l_NNLO_NNPDF_NarrowWidth_13TeV.root");
  fipPath=ggzzFIP_NNLO.fullPath();
  TFile* ggZZKFactorFile = TFile::Open(fipPath.data());
  for (unsigned int ikf=0; ikf<9; ikf++) spkfactor_ggzz_nnlo[ikf] = (TSpline3*)ggZZKFactorFile->Get(Form("sp_kfactor_%s", strZZGGKFVar[ikf].Data()))->Clone(Form("sp_kfactor_%s_NNLO", strZZGGKFVar[ikf].Data()));
  ggZZKFactorFile->Close();
  edm::FileInPath ggzzFIP_NLO("HTauTauHMuMu/AnalysisStep/data/kfactors/Kfactor_Collected_ggHZZ_2l2l_NLO_NNPDF_NarrowWidth_13TeV.root");
  fipPath=ggzzFIP_NLO.fullPath();
  ggZZKFactorFile = TFile::Open(fipPath.data());
  for (unsigned int ikf=0; ikf<9; ikf++) spkfactor_ggzz_nlo[ikf] = (TSpline3*)ggZZKFactorFile->Get(Form("sp_kfactor_%s", strZZGGKFVar[ikf].Data()))->Clone(Form("sp_kfactor_%s_NLO", strZZGGKFVar[ikf].Data()));
  ggZZKFactorFile->Close();

 
  edm::FileInPath NNLOPS_weight_path("HTauTauHMuMu/AnalysisStep/data/ggH_NNLOPS_Weights/NNLOPS_reweight.root");
  fipPath=NNLOPS_weight_path.fullPath();
  TFile* NNLOPS_weight_file = TFile::Open(fipPath.data());
  gr_NNLOPSratio_pt_powheg_0jet = (TGraphErrors*)NNLOPS_weight_file->Get("gr_NNLOPSratio_pt_powheg_0jet");
  gr_NNLOPSratio_pt_powheg_1jet = (TGraphErrors*)NNLOPS_weight_file->Get("gr_NNLOPSratio_pt_powheg_1jet");
  gr_NNLOPSratio_pt_powheg_2jet = (TGraphErrors*)NNLOPS_weight_file->Get("gr_NNLOPSratio_pt_powheg_2jet");
  gr_NNLOPSratio_pt_powheg_3jet = (TGraphErrors*)NNLOPS_weight_file->Get("gr_NNLOPSratio_pt_powheg_3jet");
      
  if(dataTag=="ULAPV"){
    preVFP=true;
  }

  //Scale factors for data/MC efficiency
  if (!skipEleDataMCWeight && isMC) { lepSFHelper = new LeptonSFHelper(preVFP); }

  if (!skipHqTWeight) {
    //HqT weights
    edm::FileInPath HqTfip("HTauTauHMuMu/AnalysisStep/test/Macros/HqTWeights.root");
    std::string fipPath=HqTfip.fullPath();
    TFile *fHqt = TFile::Open(fipPath.data(),"READ");
    h_weight = (TH2D*)fHqt->Get("wH_400")->Clone();//FIXME: Ask simon to provide the 2D histo
    fHqt->Close();
  }

  if (!skipFakeWeight) {
    //CR fake rate weight
    TString filename;filename.Form("HTauTauHMuMu/AnalysisStep/test/Macros/FR2_2011_AA_electron.root");
    if(year==2015)filename.Form("HTauTauHMuMu/AnalysisStep/test/Macros/FR2_AA_ControlSample_ABCD.root");
    edm::FileInPath fipEleZX(filename.Data());
    std::string fipPath=fipEleZX.fullPath();
    TFile *FileZXWeightEle = TFile::Open(fipPath.data(),"READ");

    filename.Form("HTauTauHMuMu/AnalysisStep/test/Macros/FR2_2011_AA_muon.root");
    if(year==2015)filename.Form("HTauTauHMuMu/AnalysisStep/test/Macros/FR2_AA_muon.root");
    edm::FileInPath fipMuZX(filename.Data());
    fipPath=fipMuZX.fullPath();
    TFile *FileZXWeightMuo = TFile::Open(fipPath.data(),"READ");

    h_ZXWeight[0]=(TH2D*)FileZXWeightEle->Get("eff_Z1ee_plus_electron")->Clone();
    h_ZXWeight[1]=(TH2D*)FileZXWeightEle->Get("eff_Z1mumu_plus_electron")->Clone();
    h_ZXWeight[2]=(TH2D*)FileZXWeightMuo->Get("eff_Z1ee_plus_muon")->Clone();
    h_ZXWeight[3]=(TH2D*)FileZXWeightMuo->Get("eff_Z1mumu_plus_muon")->Clone();

    FileZXWeightEle->Close();
    FileZXWeightMuo->Close();
  }
  
  VSmuDisc = {
    "byVLooseDeepTau2017v2p1VSmu",
    "byLooseDeepTau2017v2p1VSmu", 
    "byMediumDeepTau2017v2p1VSmu",
    "byTightDeepTau2017v2p1VSmu"
  };
  VSeDisc = {
    "byVVVLooseDeepTau2017v2p1VSe",  
    "byVVLooseDeepTau2017v2p1VSe", 
    "byVLooseDeepTau2017v2p1VSe",   
    "byLooseDeepTau2017v2p1VSe",   
    "byMediumDeepTau2017v2p1VSe",   
    "byTightDeepTau2017v2p1VSe",   
    "byVTightDeepTau2017v2p1VSe",   
    "byVVTightDeepTau2017v2p1VSe"
  };
  VSjetDisc = {
    "byVVVLooseDeepTau2017v2p1VSjet",
    "byVVLooseDeepTau2017v2p1VSjet", 
    "byVLooseDeepTau2017v2p1VSjet",  
    "byLooseDeepTau2017v2p1VSjet",   
    "byMediumDeepTau2017v2p1VSjet",  
    "byTightDeepTau2017v2p1VSjet",   
    "byVTightDeepTau2017v2p1VSjet",  
    "byVVTightDeepTau2017v2p1VSjet"
  };

}

HZZ4lNtupleMaker::~HZZ4lNtupleMaker()
{
  //clearMELABranches(); // Cleans LHE branches
  //delete lheHandler;
  delete pileUpReweight;
}


// ------------ method called for each event  ------------
void HZZ4lNtupleMaker::analyze(const edm::Event& event, const edm::EventSetup& eSetup)
{
  myTree->InitializeVariables();
  //cout<<"HZZ4lNtupleMaker:"<<theChannel<<endl;
  //----------------------------------------------------------------------
  // Analyze MC truth; collect MC weights and update counters (this is done for all generated events,
  // including those that do not pass skim, trigger etc!)

  bool gen_ZZ4lInEtaAcceptance = false;   // All 4 gen leptons in eta acceptance
  bool gen_ZZ4lInEtaPtAcceptance = false; // All 4 gen leptons in eta,pT acceptance

  const reco::Candidate * genH = 0;
  std::vector<const reco::Candidate *> genZLeps;
  std::vector<const reco::Candidate *> genAssocLeps;
  std::vector<const reco::Candidate *> genVisZLeps;
  std::vector<const reco::Candidate *> genTauNus;
    
  std::vector<float> genIso; //AT
  std::vector<const reco::GenJet *> genJet; //ATjets
  std::vector<const reco::GenJet *> genCleanedJet; //ATjets

  edm::Handle<GenEventInfoProduct> genInfo;

  if (isMC) {
    // get PU weights
    vector<Handle<std::vector< PileupSummaryInfo > > >  PupInfos; //FIXME support for miniAOD v1/v2 where name changed; catch does not work...
    event.getManyByType(PupInfos);
    Handle<std::vector< PileupSummaryInfo > > PupInfo = PupInfos.front();

    std::vector<PileupSummaryInfo>::const_iterator PVI;
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      if(PVI->getBunchCrossing() == 0) {
        NObsInt  = PVI->getPU_NumInteractions();
        NTrueInt = PVI->getTrueNumInteractions();
        break;
      }
    }

    // get PU weight
    PUWeight = pileUpReweight->weight(NTrueInt);
    PUWeight_Up = pileUpReweight->weight(NTrueInt, PileUpWeight::PUvar::VARUP);
    PUWeight_Dn = pileUpReweight->weight(NTrueInt, PileUpWeight::PUvar::VARDOWN);
     
    // L1 prefiring weights
    // From CMSSW_10_6_26 available for all the years
    // ECAL and muon weights accessible separately
    if( year == 2016 || year == 2017 || year == 2018 )
    {
       edm::Handle< double > theprefweight;
       event.getByToken(prefweight_token, theprefweight ) ;
       L1prefiringWeight =(*theprefweight);
        
       edm::Handle< double > theprefweightup;
       event.getByToken(prefweightup_token, theprefweightup ) ;
       L1prefiringWeightUp =(*theprefweightup);
        
       edm::Handle< double > theprefweightdown;
       event.getByToken(prefweightdown_token, theprefweightdown ) ;
       L1prefiringWeightDn =(*theprefweightdown);

       // ECAL and Mu weights (access for redundancy but not used)
       edm::Handle< double > theprefweightECAL;
       event.getByToken(prefweightECAL_token, theprefweightECAL ) ;
       L1prefiringWeight_ECAL =(*theprefweightECAL);

       edm::Handle< double > theprefweightupECAL;
       event.getByToken(prefweightupECAL_token, theprefweightupECAL ) ;
       L1prefiringWeightUp_ECAL =(*theprefweightupECAL);

       edm::Handle< double > theprefweightdownECAL;
       event.getByToken(prefweightdownECAL_token, theprefweightdownECAL ) ;
       L1prefiringWeightDn_ECAL =(*theprefweightdownECAL);

       edm::Handle< double > theprefweightMuon;
       event.getByToken(prefweightMuon_token, theprefweightMuon ) ;
       L1prefiringWeight_Mu =(*theprefweightMuon);

       edm::Handle< double > theprefweightupMuon;
       event.getByToken(prefweightupMuon_token, theprefweightupMuon ) ;
       L1prefiringWeightUp_Mu =(*theprefweightupMuon);

       edm::Handle< double > theprefweightdownMuon;
       event.getByToken(prefweightdownMuon_token, theprefweightdownMuon ) ;
       L1prefiringWeightDn_Mu =(*theprefweightdownMuon);
    }

    event.getByToken(genParticleToken, genParticles);
    event.getByToken(genInfoToken, genInfo);

    event.getByToken(genJetsToken, genJets); //ATjets
    event.getByToken(packedgenParticlesToken, packedgenParticles); //ATbbf
    event.getByToken(genParticleToken_bbf, genParticles_bbf); //ATbbf
    event.getByToken(GENCandidatesToken, GENCandidates);//ATMELA      

    edm::Handle<HTXS::HiggsClassification> htxs;
    event.getByToken(htxsToken,htxs);

    MCHistoryTools mch(event, sampleName, genParticles, genInfo, genJets, packedgenParticles);
    genFinalState = mch.genFinalState();
    genProcessId = mch.getProcessID();
    genHEPMCweight_NNLO = genHEPMCweight = mch.gethepMCweight(); // Overridden by LHEHandler if genHEPMCweight==1.
                                                                 // For 2017 MC, genHEPMCweight is reweighted later from NNLO to NLO
    const auto& genweights = genInfo->weights();
    if (genweights.size() > 1){
      if ((genweights.size() != 14 && genweights.size() != 46) || genweights[0] != genweights[1]){
        cms::Exception e("GenWeights");
        e << "Expected to find 1 gen weight, or 14 or 46 with the first two the same, found " << genweights.size() << ":\n";
        for (auto w : genweights) e << w << " ";
        throw e;
      }
      auto nominal = genweights[0];
      // see twiki here with definitions and order : https://twiki.cern.ch/twiki/bin/view/CMS/HowToPDF#Parton_shower_weights
      PythiaWeight_isr_muRoneoversqrt2 = genweights[24] / nominal;
      PythiaWeight_fsr_muRoneoversqrt2 = genweights[2] / nominal;
      PythiaWeight_isr_muRsqrt2 = genweights[25] / nominal;
      PythiaWeight_fsr_muRsqrt2 = genweights[3] / nominal;

      PythiaWeight_isr_muR0p5 = genweights[26] / nominal;
      PythiaWeight_fsr_muR0p5 = genweights[4] / nominal;
      PythiaWeight_isr_muR2 = genweights[27] / nominal;
      PythiaWeight_fsr_muR2 = genweights[5] / nominal;

      PythiaWeight_isr_muR0p25 = genweights[28] / nominal;
      PythiaWeight_fsr_muR0p25 = genweights[6] / nominal;
      PythiaWeight_isr_muR4 = genweights[29] / nominal;
      PythiaWeight_fsr_muR4 = genweights[7] / nominal;
    } else {
      PythiaWeight_isr_muRsqrt2 = PythiaWeight_isr_muRoneoversqrt2 = PythiaWeight_isr_muR2 =
      PythiaWeight_isr_muR0p5 = PythiaWeight_isr_muR4 = PythiaWeight_isr_muR0p25 =
      PythiaWeight_fsr_muRsqrt2 = PythiaWeight_fsr_muRoneoversqrt2 = PythiaWeight_fsr_muR2 =
      PythiaWeight_fsr_muR0p5 = PythiaWeight_fsr_muR4 = PythiaWeight_fsr_muR0p25 = 1;
    }

   htxsNJets = htxs->jets30.size();
   htxsHPt = htxs->higgs.Pt();
   htxs_stage0_cat = htxs->stage0_cat;
   htxs_stage1p0_cat = htxs->stage1_cat_pTjet30GeV;
   htxs_stage1p1_cat = htxs->stage1_1_cat_pTjet30GeV;
   htxs_stage1p2_cat = htxs->stage1_2_cat_pTjet30GeV;
   htxs_errorCode=htxs->errorCode;
   htxs_prodMode= htxs->prodMode;

   genExtInfo = mch.genAssociatedFS();

   //Information on generated candidates, will be used later
   genH = mch.genH();
   genZLeps     = mch.sortedGenZZLeps();
   genAssocLeps = mch.genAssociatedLeps();
   genFSR       = mch.genFSR();
   genVisZLeps	= mch.sortedVisGenZZLeps();
   genTauNus	= mch.genTauNus();
      
   genIso       = mch.genIso(); //AT
   genJet      = mch.GenJets(); //ATjets
   genCleanedJet = mch.GenCleanedJets(); //ATjets


    if(genH != 0){
      FillHGenInfo(genH->p4(),getHqTWeight(genH->p4().M(),genH->p4().Pt()));
    }
    else if(genZLeps.size()==4){ // for 4l events take the mass of the ZZ(4l) system
      FillHGenInfo((genZLeps.at(0)->p4()+genZLeps.at(1)->p4()+genZLeps.at(2)->p4()+genZLeps.at(3)->p4()),0);
    }

    // ATjets
    if (genJet.size() != 0) {
      for (unsigned int i = 0; i<genJet.size(); ++i){
        GenJetPt.push_back(genJet[i]->pt());
        GenJetMass.push_back(genJet[i]->mass());
        GenJetEta.push_back(genJet[i]->eta());
        GenJetPhi.push_back(genJet[i]->phi());
        GenJetRapidity.push_back(genJet[i]->rapidity());
      }
      nGenJet = genJet.size();

      for (unsigned int i = 0; i<genCleanedJet.size(); ++i){
        GenCleanedJetPt.push_back(genCleanedJet[i]->pt());
        GenCleanedJetMass.push_back(genCleanedJet[i]->mass());
        GenCleanedJetEta.push_back(genCleanedJet[i]->eta());
        GenCleanedJetPhi.push_back(genCleanedJet[i]->phi());
        GenCleanedJetRapidity.push_back(genCleanedJet[i]->rapidity());
      }
      nCleanedGenJet = genCleanedJet.size();

    }
      
    if (genFinalState!=BUGGY) {

      if (genZLeps.size()==4) {

        // "generated Zs" defined with standard pairing applied on gen leptons (genZLeps is sorted by MCHistoryTools)
        FillZGenInfo(genZLeps.at(0)->pdgId()*genZLeps.at(1)->pdgId(),
                     genZLeps.at(2)->pdgId()*genZLeps.at(3)->pdgId(),
                     genZLeps.at(0)->p4()+genZLeps.at(1)->p4(),
                     genZLeps.at(2)->p4()+genZLeps.at(3)->p4());

        math::XYZTLorentzVector genVisLep1p4,genVisLep2p4,genVisLep3p4,genVisLep4p4;
	if (abs(genVisZLeps.at(0)->pdgId())!=15 || genTauNus.at(0)==0) { genVisLep1p4=genVisZLeps.at(0)->p4(); }
	else { genVisLep1p4=genVisZLeps.at(0)->p4()-genTauNus.at(0)->p4(); }
        if (abs(genVisZLeps.at(1)->pdgId())!=15 || genTauNus.at(1)==0) { genVisLep2p4=genVisZLeps.at(1)->p4(); }
        else { genVisLep2p4=genVisZLeps.at(1)->p4()-genTauNus.at(1)->p4(); }
        if (abs(genVisZLeps.at(2)->pdgId())!=15 || genTauNus.at(2)==0) { genVisLep3p4=genVisZLeps.at(2)->p4(); }
        else { genVisLep3p4=genVisZLeps.at(2)->p4()-genTauNus.at(2)->p4(); }
        if (abs(genVisZLeps.at(3)->pdgId())!=15 || genTauNus.at(3)==0) { genVisLep4p4=genVisZLeps.at(3)->p4(); }
        else { genVisLep4p4=genVisZLeps.at(3)->p4()-genTauNus.at(3)->p4(); }

	FillVisZGenInfo(genVisZLeps.at(0)->pdgId()*genVisZLeps.at(1)->pdgId(),
		        genVisZLeps.at(2)->pdgId()*genVisZLeps.at(3)->pdgId(),
		        genVisLep1p4+genVisLep2p4, 
			genVisLep3p4+genVisLep4p4);
        // Gen leptons
        FillLepGenInfo(genZLeps.at(0)->pdgId(), genZLeps.at(1)->pdgId(), genZLeps.at(2)->pdgId(), genZLeps.at(3)->pdgId(),
           genZLeps.at(0)->p4(), genZLeps.at(1)->p4(), genZLeps.at(2)->p4(), genZLeps.at(3)->p4());

	FillVisLepGenInfo(genVisZLeps.at(0)->pdgId(), genVisZLeps.at(1)->pdgId(), genVisZLeps.at(2)->pdgId(), genVisZLeps.at(3)->pdgId(),
	   genVisLep1p4, genVisLep2p4, genVisLep3p4, genVisLep4p4);

      }

      if (genZLeps.size()==3) {

        math::XYZTLorentzVector genVisLep1p4,genVisLep2p4,genVisLep3p4;
        if (abs(genVisZLeps.at(0)->pdgId())!=15 || genTauNus.at(0)==0) { genVisLep1p4=genVisZLeps.at(0)->p4(); }
        else { genVisLep1p4=genVisZLeps.at(0)->p4()-genTauNus.at(0)->p4(); }
        if (abs(genVisZLeps.at(1)->pdgId())!=15 || genTauNus.at(1)==0) { genVisLep2p4=genVisZLeps.at(1)->p4(); }
        else { genVisLep2p4=genVisZLeps.at(1)->p4()-genTauNus.at(1)->p4(); }
        if (abs(genVisZLeps.at(2)->pdgId())!=15 || genTauNus.at(2)==0) { genVisLep3p4=genVisZLeps.at(2)->p4(); }
        else { genVisLep3p4=genVisZLeps.at(2)->p4()-genTauNus.at(2)->p4(); }

        FillLepGenInfo(genZLeps.at(0)->pdgId(), genZLeps.at(1)->pdgId(), genZLeps.at(2)->pdgId(), 0,
                       genZLeps.at(0)->p4(), genZLeps.at(1)->p4(), genZLeps.at(2)->p4(), *(new math::XYZTLorentzVector));
	FillVisLepGenInfo(genVisZLeps.at(0)->pdgId(), genVisZLeps.at(1)->pdgId(), genVisZLeps.at(2)->pdgId(), 0,
           genVisLep1p4, genVisLep2p4, genVisLep3p4, *(new math::XYZTLorentzVector));
      }
      if (genZLeps.size()==2) {

        math::XYZTLorentzVector genVisLep1p4,genVisLep2p4;
        if (abs(genVisZLeps.at(0)->pdgId())!=15 || genTauNus.at(0)==0) { genVisLep1p4=genVisZLeps.at(0)->p4(); }
        else { genVisLep1p4=genVisZLeps.at(0)->p4()-genTauNus.at(0)->p4(); }
        if (abs(genVisZLeps.at(1)->pdgId())!=15 || genTauNus.at(1)==0) { genVisLep2p4=genVisZLeps.at(1)->p4(); }
        else { genVisLep2p4=genVisZLeps.at(1)->p4()-genTauNus.at(1)->p4(); }
	
	FillLepGenInfo(genZLeps.at(0)->pdgId(), genZLeps.at(1)->pdgId(), 0, 0,
                       genZLeps.at(0)->p4(), genZLeps.at(1)->p4(), *(new math::XYZTLorentzVector), *(new math::XYZTLorentzVector));
	FillVisLepGenInfo(genVisZLeps.at(0)->pdgId(), genVisZLeps.at(1)->pdgId(), 0, 0,
           genVisLep1p4, genVisLep2p4, *(new math::XYZTLorentzVector), *(new math::XYZTLorentzVector));
      }

      if (genAssocLeps.size()==1 || genAssocLeps.size()==2) {
        FillAssocLepGenInfo(genAssocLeps);
      }
        
      if(genIso.size()==4){
        FillLepGenIso(genIso.at(0), genIso.at(1), genIso.at(2), genIso.at(3)); //AT
      }
      else if(genIso.size()==3){
        FillLepGenIso(genIso.at(0), genIso.at(1), genIso.at(2), -1);
      }
      else if(genIso.size()==2){
        FillLepGenIso(genIso.at(0), genIso.at(1), -1, -1);
      }
      else if(genIso.size()==1){
        FillLepGenIso(genIso.at(0), -1, -1, -1);
      }
      else if(genIso.size()==0){
        FillLepGenIso(-1, -1, -1, -1);
      }

      // keep track of sum of weights
      addweight(gen_sumPUWeight, PUWeight);
      addweight(gen_sumGenMCWeight, genHEPMCweight);
      addweight(gen_sumWeights, PUWeight*genHEPMCweight);

      mch.genAcceptance(gen_ZZ4lInEtaAcceptance, gen_ZZ4lInEtaPtAcceptance);
    }

    addweight(Nevt_Gen_lumiBlock, 1); // Needs to be outside the if-block

    if (genFinalState == EEEE) {
      addweight(gen_ZZ4e, 1);
      if (gen_ZZ4lInEtaAcceptance) addweight(gen_ZZ4e_EtaAcceptance, 1);
      if (gen_ZZ4lInEtaPtAcceptance) addweight(gen_ZZ4e_LeptonAcceptance, 1);
    } else if (genFinalState == MMMM) {
      addweight(gen_ZZ4mu, 1);
      if (gen_ZZ4lInEtaAcceptance) addweight(gen_ZZ4mu_EtaAcceptance, 1);
      if (gen_ZZ4lInEtaPtAcceptance) addweight(gen_ZZ4mu_LeptonAcceptance, 1);
    } else if (genFinalState == EEMM) {
      addweight(gen_ZZ2mu2e, 1);
      if (gen_ZZ4lInEtaAcceptance) addweight(gen_ZZ2mu2e_EtaAcceptance, 1);
      if (gen_ZZ4lInEtaPtAcceptance) addweight(gen_ZZ2mu2e_LeptonAcceptance, 1);
    } else if (genFinalState == MMTT){
      addweight(gen_ZZ2mu2tau, 1);
      addweight(gen_ZZ2emu2tau, 1);
      addweight(gen_ZZ2l2tau, 1);
    } else if (genFinalState == EETT){
      addweight(gen_ZZ2e2tau, 1);
      addweight(gen_ZZ2emu2tau, 1);
      addweight(gen_ZZ2l2tau, 1);
    } else if (genFinalState == TTTT){
      addweight(gen_ZZ4tau, 1);
      addweight(gen_ZZ2l2tau, 1);
    } else if (genFinalState == BUGGY){ // handle MCFM ZZ->4tau mZ<2mtau bug
      addweight(gen_BUGGY, 1);
      return; // BUGGY events are skipped
    } else {
      addweight(gen_Unknown, 1);
    }

// End of MC history analysis ------------------------------------------
  } else {
    ++Nevt_Gen_lumiBlock; // keep track of # events for data as well
  }



  // Get candidate collection
  edm::Handle<edm::View<pat::CompositeCandidate> > candHandle;
  event.getByToken(candToken, candHandle);
  if(candHandle.failedToGet()) {
    if(is_loose_ele_selection) return; // The collection can be missing in this case since we have a filter to skip the module when a regular candidate is present.
    else edm::LogError("") << "ZZ collection not found in non-loose electron flow. This should never happen";
  }
  const edm::View<pat::CompositeCandidate>* cands = candHandle.product();

  // For Z+L CRs, we want only events with exactly 1 Z+l candidate. FIXME: this has to be reviewed.
  if (theChannel==ZL && cands->size() != 1) return;


  // Retrieve trigger results
  Handle<edm::TriggerResults> triggerResults;
  event.getByToken(triggerResultToken, triggerResults);

  bool failed = false;

  // Apply MC filter (skip event)
  // Heshy note: I'm not turning return into failed = true because it looks like it's applied even if !skipEmptyEvents.
  //             It only does anything if the MCFILTER variable is set in the csv file, which is not currently the case.
  if (isMC && !(myHelper.passMCFilter(event,triggerResults))) return;

  // Apply skim
  bool evtPassSkim = myHelper.passSkim(event,triggerResults,trigWord);
  if (applySkim && !evtPassSkim) failed = true;       //but gen information will still be recorded if failedTreeLevel != 0

  // Apply trigger request (skip event)
  bool evtPassTrigger = myHelper.passTrigger(event,triggerResults,trigWord);
  if (applyTrigger && !evtPassTrigger) failed = true; //but gen information will still be recorded if failedTreeLevel != 0
	
  // Apply MET trigger request (skip event)
  evtPassMETTrigger = myHelper.passMETTrigger(event,triggerResults);

  if (skipEmptyEvents && !failedTreeLevel && (cands->size() == 0 || failed)) return; // Skip events with no candidate, unless skipEmptyEvents = false or failedTreeLevel != 0

  //Fill MC truth information
  if (isMC) FillKFactors(genInfo, genZLeps);

  // General event information
  RunNumber=event.id().run();
  LumiNumber=event.luminosityBlock();
  EventNumber=event.id().event();
  xsection=xsec;
  genxsection=genxsec;
  genbranchingratio=genbr;

  // General event information
  RunNumber=event.id().run();
  LumiNumber=event.luminosityBlock();
  EventNumber=event.id().event();
  xsection=xsec;
  genxsection=genxsec;
  genbranchingratio=genbr;

  // Primary vertices
  Handle<vector<reco::Vertex> > vertices;
  event.getByToken(vtxToken,vertices);
  Nvtx=vertices->size();

  // Jets (cleaned wrt all tight isolated leptons)
  Handle<edm::View<pat::Jet> > CleanedJets;
  event.getByToken(jetToken, CleanedJets);
  vector<const pat::Jet*> cleanedJets;
  for(edm::View<pat::Jet>::const_iterator jet = CleanedJets->begin(); jet != CleanedJets->end(); ++jet){
    cleanedJets.push_back(&*jet);
  }
   
   // Photons
   Handle<pat::PhotonCollection> photonCands;
   event.getByToken(photonToken, photonCands);
   vector<const pat::Photon*> photons;
   
   for(unsigned int i = 0; i< photonCands->size(); ++i){
      const pat::Photon* photon = &((*photonCands)[i]);
      photons.push_back(&*photon);
   }

   
   if (writePhotons){
      for (unsigned i=0; i<photons.size(); ++i) {
            FillPhoton(year, *(photons.at(i)));
         }
   }
      
  // MET
  Handle<pat::METCollection> metHandle;
  event.getByToken(metToken, metHandle);
  Handle<math::Error<2>::type> covHandle;
  event.getByToken(theCovTag, covHandle);
  Handle<double> METdxUPTESHandle, METdyUPTESHandle, METdxDOWNTESHandle, METdyDOWNTESHandle, METdxUPEESHandle, METdyUPEESHandle, METdxDOWNEESHandle, METdyDOWNEESHandle, METdxUPMESHandle, METdyUPMESHandle, METdxDOWNMESHandle, METdyDOWNMESHandle, METdxUPJESHandle, METdyUPJESHandle, METdxDOWNJESHandle, METdyDOWNJESHandle, METdxUPJERHandle, METdyUPJERHandle, METdxDOWNJERHandle, METdyDOWNJERHandle;
  
  event.getByToken(theMETdxUPTESTag, METdxUPTESHandle);
  event.getByToken(theMETdyUPTESTag, METdyUPTESHandle);
  event.getByToken(theMETdxDOWNTESTag, METdxDOWNTESHandle);
  event.getByToken(theMETdyDOWNTESTag, METdyDOWNTESHandle);
  event.getByToken(theMETdxUPEESTag, METdxUPEESHandle);
  event.getByToken(theMETdyUPEESTag, METdyUPEESHandle);
  event.getByToken(theMETdxDOWNEESTag, METdxDOWNEESHandle);
  event.getByToken(theMETdyDOWNEESTag, METdyDOWNEESHandle);
  event.getByToken(theMETdxUPMESTag, METdxUPMESHandle);
  event.getByToken(theMETdyUPMESTag, METdyUPMESHandle);
  event.getByToken(theMETdxDOWNMESTag, METdxDOWNMESHandle);
  event.getByToken(theMETdyDOWNMESTag, METdyDOWNMESHandle);
  event.getByToken(theMETdxUPJESTag, METdxUPJESHandle);
  event.getByToken(theMETdyUPJESTag, METdyUPJESHandle);
  event.getByToken(theMETdxDOWNJESTag, METdxDOWNJESHandle);
  event.getByToken(theMETdyDOWNJESTag, METdyDOWNJESHandle);
  event.getByToken(theMETdxUPJERTag, METdxUPJERHandle);
  event.getByToken(theMETdyUPJERTag, METdyUPJERHandle);
  event.getByToken(theMETdxDOWNJERTag, METdxDOWNJERHandle);
  event.getByToken(theMETdyDOWNJERTag, METdyDOWNJERHandle);
    
  const pat::MET &met = metHandle->front();
  PFMET = met.pt();
  PFMETPhi = met.phi();
  METx = met.px();
  METy = met.py();
  covMET[0][0] = (*covHandle)(0,0);
  covMET[1][0] = (*covHandle)(1,0);
  covMET[0][1] = covMET[1][0]; // (1,0) is the only one saved
  covMET[1][1] = (*covHandle)(1,1);
  METxUPTES=*METdxUPTESHandle;
  METyUPTES=*METdyUPTESHandle;
  METxDOWNTES=*METdxDOWNTESHandle;
  METyDOWNTES=*METdyDOWNTESHandle;
  METxUPEES=*METdxUPEESHandle;
  METyUPEES=*METdyUPEESHandle;
  METxDOWNEES=*METdxDOWNEESHandle;
  METyDOWNEES=*METdyDOWNEESHandle;
  METxUPMES=*METdxUPMESHandle;
  METyUPMES=*METdyUPMESHandle;
  METxDOWNMES=*METdxDOWNMESHandle;
  METyDOWNMES=*METdyDOWNMESHandle;
  METxUPJES=*METdxUPJESHandle;
  METyUPJES=*METdyUPJESHandle;
  METxDOWNJES=*METdxDOWNJESHandle;
  METyDOWNJES=*METdyDOWNJESHandle;
  METxUPJER=*METdxUPJERHandle;
  METyUPJER=*METdyUPJERHandle;
  METxDOWNJER=*METdxDOWNJERHandle;
  METyDOWNJER=*METdyDOWNJERHandle;
    
  GenMET=GenMETPhi=-99;
  if (isMC && met.genMET()){
      GenMET = met.genMET()->pt();
      GenMETPhi = met.genMET()->phi();
  }
  else if (isMC){
      cms::Exception e("GenMET");
      e << "No met.genMET!";
      throw e;
    }
//   if (metHandle.isValid()){
//     const pat::MET &met = metHandle->front();
//   }

  // number of reconstructed leptons
  edm::Handle<pat::MuonCollection> muonHandle;
  event.getByToken(muonToken, muonHandle);
  for(unsigned int i = 0; i< muonHandle->size(); ++i){
    const pat::Muon* m = &((*muonHandle)[i]);
    if(m->pt()>5 && m->isPFMuon()) // these cuts are implicit in miniAOD
      NRecoMu++;
  }
  edm::Handle<pat::ElectronCollection> electronHandle;
  event.getByToken(electronToken, electronHandle);
  for(unsigned int i = 0; i< electronHandle->size(); ++i){
    const pat::Electron* e = &((*electronHandle)[i]);
    if(e->pt()>5) // this cut is implicit in miniAOD
      NRecoEle++;
  }

  if(isMC && apply_QCD_GGF_UNCERT)
  {

	  

    if (htxsNJets==0)
    {
      ggH_NNLOPS_weight = gr_NNLOPSratio_pt_powheg_0jet->Eval(min((double) htxsHPt, 125.0));
      ggH_NNLOPS_weight_unc=(gr_NNLOPSratio_pt_powheg_0jet->GetErrorY(FindBinValue(gr_NNLOPSratio_pt_powheg_0jet, min((double) htxsHPt, 125.0))))/ggH_NNLOPS_weight;

    }
	  else if (htxsNJets==1)
	  {
		  ggH_NNLOPS_weight = gr_NNLOPSratio_pt_powheg_1jet->Eval(min((double)htxsHPt,625.0));
		  ggH_NNLOPS_weight_unc=(gr_NNLOPSratio_pt_powheg_1jet->GetErrorY(FindBinValue(gr_NNLOPSratio_pt_powheg_1jet,min((double)htxsHPt,125.0))))/ggH_NNLOPS_weight;
	  }
	  else if (htxsNJets==2)
	  {
		  ggH_NNLOPS_weight = gr_NNLOPSratio_pt_powheg_2jet->Eval(min((double)htxsHPt,800.0));
		  ggH_NNLOPS_weight_unc=(gr_NNLOPSratio_pt_powheg_2jet->GetErrorY(FindBinValue(gr_NNLOPSratio_pt_powheg_2jet,min((double)htxsHPt,125.0))))/ggH_NNLOPS_weight;
	  }
	  else if (htxsNJets>=3)
	  {
		  ggH_NNLOPS_weight = gr_NNLOPSratio_pt_powheg_3jet->Eval(min((double)htxsHPt,925.0));
		  ggH_NNLOPS_weight_unc=(gr_NNLOPSratio_pt_powheg_3jet->GetErrorY(FindBinValue(gr_NNLOPSratio_pt_powheg_3jet,min((double)htxsHPt,125.0))))/ggH_NNLOPS_weight;
	  }
	  else
	  {
		  ggH_NNLOPS_weight = 1.0;
		  ggH_NNLOPS_weight_unc = 0.0;
	  }
	  std::vector<double> qcd_ggF_uncertSF_tmp;
	  qcd_ggF_uncertSF.clear();
     ////////////////////////////////////////////////////////////
    //////////////////     CHECK THIS!!!!!    //////////////////
    // Why is this done with the STXS 1.0 bins uncertainties? //
    //////////////////     CHECK THIS!!!!!    //////////////////
    ////////////////////////////////////////////////////////////

	  qcd_ggF_uncertSF_tmp = qcd_ggF_uncertSF_2017(htxsNJets, htxsHPt, htxs_stage1p0_cat);
	  qcd_ggF_uncertSF = std::vector<float>(qcd_ggF_uncertSF_tmp.begin(),qcd_ggF_uncertSF_tmp.end());


  }
   
  //Loop on the candidates
  vector<Int_t> CRFLAG(cands->size());
  for( edm::View<pat::CompositeCandidate>::const_iterator cand = cands->begin(); cand != cands->end(); ++cand) {
    if (failed) break; //don't waste time on this
    size_t icand= cand-cands->begin();

    //    int candChannel = cand->userFloat("candChannel"); // This is currently the product of pdgId of leptons (eg 14641, 28561, 20449)

    if (theChannel==ZLL) {
      // Cross check region for Z + 1 loose electron + 1 loose TLE (defined only in loose_ele paths) 
      if (is_loose_ele_selection) {
        if (cand->userFloat("isBestCRZLL")&&cand->userFloat("CRZLL")) set_bit(CRFLAG[icand],ZLL);
      }

      // AA CRs
      if (cand->userFloat("isBestCRZLLss")&&cand->userFloat("CRZLLss")) set_bit(CRFLAG[icand],CRZLLss);

      // A CRs
      if (cand->userFloat("isBestCRZLLos_2P2F")&&cand->userFloat("CRZLLos_2P2F")) set_bit(CRFLAG[icand],CRZLLos_2P2F);
      if (cand->userFloat("isBestCRZLLos_3P1F")&&cand->userFloat("CRZLLos_3P1F")) set_bit(CRFLAG[icand],CRZLLos_3P1F);

      if (CRFLAG[icand]) { // This candidate belongs to one of the CRs: perform additional jet cleaning.
        // Note that this is (somewhat incorrectly) done per-event, so there could be some over-cleaning in events with >1 CR candidate.
        for (unsigned i=0; i<cleanedJets.size(); ++i) {
          if (cleanedJets[i]!=0  && (!jetCleaner::isGood(*cand, *(cleanedJets[i])))) {
            cleanedJets[i]=0;
          }
        }
      }
    }
  }

  // Count and store jets, after additional cleaning for CRs...
  for (unsigned i=0; i<cleanedJets.size(); ++i) {
    if (cleanedJets[i]==0) {
      continue; // Jet has been suppressed by additional cleaning
    }

    ++nCleanedJets;

    // count jes up/down njets pt30
    float jes_unc = cleanedJets[i]->userFloat("jes_unc");
    float jes_unc_Total = cleanedJets[i]->userFloat("jes_unc_split_Total");
    float jes_unc_Abs = cleanedJets[i]->userFloat("jes_unc_split_Abs");
    float jes_unc_Abs_year = cleanedJets[i]->userFloat("jes_unc_split_Abs_year");
    float jes_unc_BBEC1 = cleanedJets[i]->userFloat("jes_unc_split_BBEC1");
    float jes_unc_BBEC1_year = cleanedJets[i]->userFloat("jes_unc_split_BBEC1_year");
    float jes_unc_EC2 = cleanedJets[i]->userFloat("jes_unc_split_EC2");
    float jes_unc_EC2_year = cleanedJets[i]->userFloat("jes_unc_split_EC2_year");
    float jes_unc_FlavQCD = cleanedJets[i]->userFloat("jes_unc_split_FlavQCD");
    float jes_unc_HF = cleanedJets[i]->userFloat("jes_unc_split_HF");
    float jes_unc_HF_year = cleanedJets[i]->userFloat("jes_unc_split_HF_year");
    float jes_unc_RelBal = cleanedJets[i]->userFloat("jes_unc_split_RelBal");
    float jes_unc_RelSample_year = cleanedJets[i]->userFloat("jes_unc_split_RelSample_year");
      
    float pt_nominal = cleanedJets[i]->pt();
    float pt_jes_up = pt_nominal * (1.0 + jes_unc);
    float pt_jes_up_Total = pt_nominal * (1.0 + jes_unc_Total);
    float pt_jes_up_Abs = pt_nominal * (1.0 + jes_unc_Abs);
    float pt_jes_up_Abs_year = pt_nominal * (1.0 + jes_unc_Abs_year);
    float pt_jes_up_BBEC1 = pt_nominal * (1.0 + jes_unc_BBEC1);
    float pt_jes_up_BBEC1_year = pt_nominal * (1.0 + jes_unc_BBEC1_year);
    float pt_jes_up_EC2 = pt_nominal * (1.0 + jes_unc_EC2);
    float pt_jes_up_EC2_year = pt_nominal * (1.0 + jes_unc_EC2_year);
    float pt_jes_up_FlavQCD = pt_nominal * (1.0 + jes_unc_FlavQCD);
    float pt_jes_up_HF = pt_nominal * (1.0 + jes_unc_HF);
    float pt_jes_up_HF_year = pt_nominal * (1.0 + jes_unc_HF_year);
    float pt_jes_up_RelBal = pt_nominal * (1.0 + jes_unc_RelBal);
    float pt_jes_up_RelSample_year = pt_nominal * (1.0 + jes_unc_RelSample_year);
    float pt_jes_dn = pt_nominal * (1.0 - jes_unc);
    float pt_jes_dn_Total = pt_nominal * (1.0 - jes_unc_Total);
    float pt_jes_dn_Abs = pt_nominal * (1.0 - jes_unc_Abs);
    float pt_jes_dn_Abs_year = pt_nominal * (1.0 - jes_unc_Abs_year);
    float pt_jes_dn_BBEC1 = pt_nominal * (1.0 - jes_unc_BBEC1);
    float pt_jes_dn_BBEC1_year = pt_nominal * (1.0 - jes_unc_BBEC1_year);
    float pt_jes_dn_EC2 = pt_nominal * (1.0 - jes_unc_EC2);
    float pt_jes_dn_EC2_year = pt_nominal * (1.0 - jes_unc_EC2_year);
    float pt_jes_dn_FlavQCD = pt_nominal * (1.0 - jes_unc_FlavQCD);
    float pt_jes_dn_HF = pt_nominal * (1.0 - jes_unc_HF);
    float pt_jes_dn_HF_year = pt_nominal * (1.0 - jes_unc_HF_year);
    float pt_jes_dn_RelBal = pt_nominal * (1.0 - jes_unc_RelBal);
    float pt_jes_dn_RelSample_year = pt_nominal * (1.0 - jes_unc_RelSample_year);

    if(pt_nominal>30){
      ++nCleanedJetsPt30;
      if(cleanedJets[i]->userFloat("isBtagged")) ++nCleanedJetsPt30BTagged;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF_Up")) ++nCleanedJetsPt30BTagged_bTagSFUp;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF_Dn")) ++nCleanedJetsPt30BTagged_bTagSFDn;
    }
    if(pt_jes_up>30){
      ++nCleanedJetsPt30_jesUp;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesUp;
    }
    if(pt_jes_up_Total>30){
      ++nCleanedJetsPt30_jesUp_Total;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesUp_Total;
    }
    if(pt_jes_up_Abs>30){
      ++nCleanedJetsPt30_jesUp_Abs;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesUp_Abs;
    }
    if(pt_jes_up_Abs_year>30){
      ++nCleanedJetsPt30_jesUp_Abs_year;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesUp_Abs_year;
    }
    if(pt_jes_up_BBEC1>30){
      ++nCleanedJetsPt30_jesUp_BBEC1;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesUp_BBEC1;
    }
    if(pt_jes_up_BBEC1_year>30){
      ++nCleanedJetsPt30_jesUp_BBEC1_year;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesUp_BBEC1_year;
    }
    if(pt_jes_up_EC2>30){
      ++nCleanedJetsPt30_jesUp_EC2;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesUp_EC2;
    }
    if(pt_jes_up_EC2_year>30){
      ++nCleanedJetsPt30_jesUp_EC2_year;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesUp_EC2_year;
    }
    if(pt_jes_up_FlavQCD>30){
      ++nCleanedJetsPt30_jesUp_FlavQCD;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesUp_FlavQCD;
    }
    if(pt_jes_up_HF>30){
      ++nCleanedJetsPt30_jesUp_HF;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesUp_HF;
    }
    if(pt_jes_up_HF_year>30){
      ++nCleanedJetsPt30_jesUp_HF_year;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesUp_HF_year;
    }
    if(pt_jes_up_RelBal>30){
      ++nCleanedJetsPt30_jesUp_RelBal;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesUp_RelBal;
    }
    if(pt_jes_up_RelSample_year>30){
      ++nCleanedJetsPt30_jesUp_RelSample_year;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesUp_RelSample_year;
    }
    if(pt_jes_dn>30){
      ++nCleanedJetsPt30_jesDn;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesDn;
    }
    if(pt_jes_dn_Total>30){
      ++nCleanedJetsPt30_jesDn_Total;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesDn_Total;
    }
    if(pt_jes_dn_Abs>30){
      ++nCleanedJetsPt30_jesDn_Abs;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesDn_Abs;
    }
    if(pt_jes_dn_Abs_year>30){
      ++nCleanedJetsPt30_jesDn_Abs_year;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesDn_Abs_year;
    }
    if(pt_jes_dn_BBEC1>30){
      ++nCleanedJetsPt30_jesDn_BBEC1;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesDn_BBEC1;
    }
    if(pt_jes_dn_BBEC1_year>30){
      ++nCleanedJetsPt30_jesDn_BBEC1_year;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesDn_BBEC1_year;
    }
    if(pt_jes_dn_EC2>30){
      ++nCleanedJetsPt30_jesDn_EC2;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesDn_EC2;
    }
    if(pt_jes_dn_EC2_year>30){
      ++nCleanedJetsPt30_jesDn_EC2_year;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesDn_EC2_year;
    }
    if(pt_jes_dn_FlavQCD>30){
      ++nCleanedJetsPt30_jesDn_FlavQCD;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesDn_FlavQCD;
    }
    if(pt_jes_dn_HF>30){
      ++nCleanedJetsPt30_jesDn_HF;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesDn_HF;
    }
    if(pt_jes_dn_HF_year>30){
      ++nCleanedJetsPt30_jesDn_HF_year;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesDn_HF_year;
    }
    if(pt_jes_dn_RelBal>30){
      ++nCleanedJetsPt30_jesDn_RelBal;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesDn_RelBal;
    }
    if(pt_jes_dn_RelSample_year>30){
      ++nCleanedJetsPt30_jesDn_RelSample_year;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesDn_RelSample_year;
    }

    // count jer up/down njets pt30
    float pt_jer_up = cleanedJets[i]->userFloat("pt_jerup");
    float pt_jer_dn = cleanedJets[i]->userFloat("pt_jerdn");
     
    if(pt_jer_up>30){
       ++nCleanedJetsPt30_jerUp;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jerUp;
    }
    if(pt_jer_dn>30){
       ++nCleanedJetsPt30_jerDn;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jerDn;
    }
     
    if (writeJets) FillJet(*(cleanedJets.at(i))); // No additional pT cut (for JEC studies)
  }

  // Now we can write the variables for candidates
  int nFilled=0;
  for( edm::View<pat::CompositeCandidate>::const_iterator cand = cands->begin(); cand != cands->end(); ++cand) {
    if (failed) break; //don't waste time on this
    size_t icand= cand-cands->begin();

    if (!( theChannel==ZL || CRFLAG[icand] || (bool)(cand->userFloat("isBestCand")) )) continue; // Skip events other than the best cand (or CR candidates in the CR)

    //For the SR, also fold information about acceptance in CRflag.
    if (isMC && (theChannel==ZZ)) {
      if (gen_ZZ4lInEtaAcceptance)   set_bit(CRFLAG[icand],28);
      if (gen_ZZ4lInEtaPtAcceptance) set_bit(CRFLAG[icand],29);
    }
    FillCandidate(*cand, evtPassTrigger&&evtPassSkim, event, CRFLAG[icand]);

    // Fill the candidate as one entry in the tree. Do not reinitialize the event variables, as in CRs
    // there could be several candidates per event.
    myTree->FillCurrentTree(true);
    ++nFilled;
  }

  // If no candidate was filled but we still want to keep gen-level and weights, we need to fill one entry anyhow.
  if (nFilled==0) {
    if (skipEmptyEvents==false)
      myTree->FillCurrentTree(true);
    else
      myTree->FillCurrentTree(false); //puts it in the failed tree if there is one
  }
}


void HZZ4lNtupleMaker::FillJet(const pat::Jet& jet)
{
   JetPt  .push_back( jet.pt());
   JetEta .push_back( jet.eta());
   JetPhi .push_back( jet.phi());
   JetMass .push_back( jet.p4().M());
   JetEnergy .push_back( jet.p4().energy());
   JetBTagger .push_back( jet.userFloat("bTagger"));
   JetIsBtagged .push_back( jet.userFloat("isBtagged"));
   JetIsBtaggedWithSF .push_back( jet.userFloat("isBtaggedWithSF"));
   JetIsBtaggedWithSFUp .push_back( jet.userFloat("isBtaggedWithSF_Up"));
   JetIsBtaggedWithSFDn .push_back( jet.userFloat("isBtaggedWithSF_Dn"));
   JetQGLikelihood .push_back( jet.userFloat("qgLikelihood"));
   if(addQGLInputs){
     JetAxis2 .push_back( jet.userFloat("axis2"));
     JetMult .push_back( jet.userFloat("mult"));
     JetPtD .push_back( jet.userFloat("ptD"));
   }
   JetSigma .push_back(jet.userFloat("jes_unc"));
   JetSigma_Total .push_back(jet.userFloat("jes_unc_split_Total"));
   JetSigma_Abs .push_back(jet.userFloat("jes_unc_split_Abs"));
   JetSigma_Abs_year .push_back(jet.userFloat("jes_unc_split_Abs_year"));
   JetSigma_BBEC1 .push_back(jet.userFloat("jes_unc_split_BBEC1"));
   JetSigma_BBEC1_year .push_back(jet.userFloat("jes_unc_split_BBEC1_year"));
   JetSigma_EC2 .push_back(jet.userFloat("jes_unc_split_EC2"));
   JetSigma_EC2_year .push_back(jet.userFloat("jes_unc_split_EC2_year"));
   JetSigma_FlavQCD .push_back(jet.userFloat("jes_unc_split_FlavQCD"));
   JetSigma_HF .push_back(jet.userFloat("jes_unc_split_HF"));
   JetSigma_HF_year .push_back(jet.userFloat("jes_unc_split_HF_year"));
   JetSigma_RelBal .push_back(jet.userFloat("jes_unc_split_RelBal"));
   JetSigma_RelSample_year .push_back(jet.userFloat("jes_unc_split_RelSample_year"));
    
   JetRawPt  .push_back( jet.userFloat("RawPt"));
   JetPtJEC_noJER .push_back( jet.userFloat("pt_JEC_noJER"));
   
   JetJESUp .push_back(jet.userFloat("pt_jesup"));
   JetJESUp_Total .push_back(jet.userFloat("pt_jesup_split_Total"));
   JetJESUp_Abs .push_back(jet.userFloat("pt_jesup_split_Abs"));
   JetJESUp_Abs_year .push_back(jet.userFloat("pt_jesup_split_Abs_year"));
   JetJESUp_BBEC1 .push_back(jet.userFloat("pt_jesup_split_BBEC1"));
   JetJESUp_BBEC1_year .push_back(jet.userFloat("pt_jesup_split_BBEC1_year"));
   JetJESUp_EC2 .push_back(jet.userFloat("pt_jesup_split_EC2"));
   JetJESUp_EC2_year .push_back(jet.userFloat("pt_jesup_split_EC2_year"));
   JetJESUp_FlavQCD .push_back(jet.userFloat("pt_jesup_split_FlavQCD"));
   JetJESUp_HF .push_back(jet.userFloat("pt_jesup_split_HF"));
   JetJESUp_HF_year .push_back(jet.userFloat("pt_jesup_split_HF_year"));
   JetJESUp_RelBal .push_back(jet.userFloat("pt_jesup_split_RelBal"));
   JetJESUp_RelSample_year .push_back(jet.userFloat("pt_jesup_split_RelSample_year"));
   JetJESDown .push_back(jet.userFloat("pt_jesdn"));
   JetJESDown_Total .push_back(jet.userFloat("pt_jesdn_split_Total"));
   JetJESDown_Abs .push_back(jet.userFloat("pt_jesdn_split_Abs"));
   JetJESDown_Abs_year .push_back(jet.userFloat("pt_jesdn_split_Abs_year"));
   JetJESDown_BBEC1 .push_back(jet.userFloat("pt_jesdn_split_BBEC1"));
   JetJESDown_BBEC1_year .push_back(jet.userFloat("pt_jesdn_split_BBEC1_year"));
   JetJESDown_EC2 .push_back(jet.userFloat("pt_jesdn_split_EC2"));
   JetJESDown_EC2_year .push_back(jet.userFloat("pt_jesdn_split_EC2_year"));
   JetJESDown_FlavQCD .push_back(jet.userFloat("pt_jesdn_split_FlavQCD"));
   JetJESDown_HF .push_back(jet.userFloat("pt_jesdn_split_HF"));
   JetJESDown_HF_year .push_back(jet.userFloat("pt_jesdn_split_HF_year"));
   JetJESDown_RelBal .push_back(jet.userFloat("pt_jesdn_split_RelBal"));
   JetJESDown_RelSample_year .push_back(jet.userFloat("pt_jesdn_split_RelSample_year"));

   JetJERUp .push_back(jet.userFloat("pt_jerup"));
   JetJERDown .push_back(jet.userFloat("pt_jerdn"));
    
   JetID.push_back(jet.userFloat("JetID"));
   JetPUID.push_back(jet.userFloat("PUjetID"));
   JetPUID_score.push_back(jet.userFloat("PUjetID_score"));

   if (jet.hasUserFloat("pileupJetIdUpdated:fullDiscriminant")) { // if JEC is reapplied, we set this
     JetPUValue.push_back(jet.userFloat("pileupJetIdUpdated:fullDiscriminant"));
   } else {
     JetPUValue.push_back(jet.userFloat("pileupJetId:fullDiscriminant"));
   }
   

   JetHadronFlavour .push_back(jet.hadronFlavour());
   JetPartonFlavour .push_back(jet.partonFlavour());
}

void HZZ4lNtupleMaker::FillPhoton(int year, const pat::Photon& photon)
{
   PhotonPt  .push_back( photon.pt());
   PhotonEta .push_back( photon.eta());
   PhotonPhi .push_back( photon.phi());
   
   PhotonIsCutBasedLooseID .push_back( PhotonIDHelper::isCutBasedID_Loose(year, photon) );
}

float HZZ4lNtupleMaker::EvalSpline(TSpline3* const& sp, float xval){
  double xmin = sp->GetXmin();
  double xmax = sp->GetXmax();
  double res=0;
  if (xval<xmin){
    res=sp->Eval(xmin);
    double deriv=sp->Derivative(xmin);
    res += deriv*(xval-xmin);
  }
  else if (xval>xmax){
    res=sp->Eval(xmax);
    double deriv=sp->Derivative(xmax);
    res += deriv*(xval-xmax);
  }
  else res=sp->Eval(xval);
  return res;
}

void HZZ4lNtupleMaker::FillKFactors(edm::Handle<GenEventInfoProduct>& genInfo, std::vector<const reco::Candidate *>& genZLeps){
  KFactor_QCD_ggZZ_Nominal=1;
  KFactor_QCD_ggZZ_PDFScaleDn=1;
  KFactor_QCD_ggZZ_PDFScaleUp=1;
  KFactor_QCD_ggZZ_QCDScaleDn=1;
  KFactor_QCD_ggZZ_QCDScaleUp=1;
  KFactor_QCD_ggZZ_AsDn=1;
  KFactor_QCD_ggZZ_AsUp=1;
  KFactor_QCD_ggZZ_PDFReplicaDn=1;
  KFactor_QCD_ggZZ_PDFReplicaUp=1;
  KFactor_QCD_qqZZ_dPhi=1;
  KFactor_QCD_qqZZ_M=1;
  KFactor_QCD_qqZZ_Pt=1;
  KFactor_EW_qqZZ=1;
  KFactor_EW_qqZZ_unc=0;

  if (isMC){
    GenEventInfoProduct  genInfoP = *(genInfo.product());
    if (apply_K_NNLOQCD_ZZGG>0 && apply_K_NNLOQCD_ZZGG!=3){
      if (spkfactor_ggzz_nnlo[0]!=0) KFactor_QCD_ggZZ_Nominal = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nnlo[0], GenHMass);
      if (spkfactor_ggzz_nnlo[1]!=0) KFactor_QCD_ggZZ_PDFScaleDn = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nnlo[1], GenHMass);
      if (spkfactor_ggzz_nnlo[2]!=0) KFactor_QCD_ggZZ_PDFScaleUp = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nnlo[2], GenHMass);
      if (spkfactor_ggzz_nnlo[3]!=0) KFactor_QCD_ggZZ_QCDScaleDn = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nnlo[3], GenHMass);
      if (spkfactor_ggzz_nnlo[4]!=0) KFactor_QCD_ggZZ_QCDScaleUp = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nnlo[4], GenHMass);
      if (spkfactor_ggzz_nnlo[5]!=0) KFactor_QCD_ggZZ_AsDn = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nnlo[5], GenHMass);
      if (spkfactor_ggzz_nnlo[6]!=0) KFactor_QCD_ggZZ_AsUp = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nnlo[6], GenHMass);
      if (spkfactor_ggzz_nnlo[7]!=0) KFactor_QCD_ggZZ_PDFReplicaDn = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nnlo[7], GenHMass);
      if (spkfactor_ggzz_nnlo[8]!=0) KFactor_QCD_ggZZ_PDFReplicaUp = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nnlo[8], GenHMass);
      if (apply_K_NNLOQCD_ZZGG==2){
        if (spkfactor_ggzz_nlo[0]!=0){
          float divisor = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nlo[0], GenHMass);
          KFactor_QCD_ggZZ_Nominal /= divisor;
          KFactor_QCD_ggZZ_PDFScaleDn /= divisor;
          KFactor_QCD_ggZZ_PDFScaleUp /= divisor;
          KFactor_QCD_ggZZ_QCDScaleDn /= divisor;
          KFactor_QCD_ggZZ_QCDScaleUp /= divisor;
          KFactor_QCD_ggZZ_AsDn /= divisor;
          KFactor_QCD_ggZZ_AsUp /= divisor;
          KFactor_QCD_ggZZ_PDFReplicaDn /= divisor;
          KFactor_QCD_ggZZ_PDFReplicaUp /= divisor;
        }
        else{
          KFactor_QCD_ggZZ_Nominal=0;
          KFactor_QCD_ggZZ_PDFScaleDn=0;
          KFactor_QCD_ggZZ_PDFScaleUp=0;
          KFactor_QCD_ggZZ_QCDScaleDn=0;
          KFactor_QCD_ggZZ_QCDScaleUp=0;
          KFactor_QCD_ggZZ_AsDn=0;
          KFactor_QCD_ggZZ_AsUp=0;
          KFactor_QCD_ggZZ_PDFReplicaDn=0;
          KFactor_QCD_ggZZ_PDFReplicaUp=0;
        }
      }
    }
    else if (apply_K_NNLOQCD_ZZGG==3){
      if (spkfactor_ggzz_nlo[0]!=0) KFactor_QCD_ggZZ_Nominal = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nlo[0], GenHMass);
      if (spkfactor_ggzz_nlo[1]!=0) KFactor_QCD_ggZZ_PDFScaleDn = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nlo[1], GenHMass);
      if (spkfactor_ggzz_nlo[2]!=0) KFactor_QCD_ggZZ_PDFScaleUp = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nlo[2], GenHMass);
      if (spkfactor_ggzz_nlo[3]!=0) KFactor_QCD_ggZZ_QCDScaleDn = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nlo[3], GenHMass);
      if (spkfactor_ggzz_nlo[4]!=0) KFactor_QCD_ggZZ_QCDScaleUp = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nlo[4], GenHMass);
      if (spkfactor_ggzz_nlo[5]!=0) KFactor_QCD_ggZZ_AsDn = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nlo[5], GenHMass);
      if (spkfactor_ggzz_nlo[6]!=0) KFactor_QCD_ggZZ_AsUp = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nlo[6], GenHMass);
      if (spkfactor_ggzz_nlo[7]!=0) KFactor_QCD_ggZZ_PDFReplicaDn = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nlo[7], GenHMass);
      if (spkfactor_ggzz_nlo[8]!=0) KFactor_QCD_ggZZ_PDFReplicaUp = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nlo[8], GenHMass);
    }

    if (genFinalState!=BUGGY){
      if (genZLeps.size()==4) {
        // Calculate NNLO/NLO QCD K factors for qqZZ
        if (apply_K_NNLOQCD_ZZQQB){
          bool sameflavor=(genZLeps.at(0)->pdgId()*genZLeps.at(1)->pdgId() == genZLeps.at(2)->pdgId()*genZLeps.at(3)->pdgId());
          KFactor_QCD_qqZZ_dPhi = kfactor_qqZZ_qcd_dPhi(fabs(GenZ1Phi-GenZ2Phi), (sameflavor) ? 1 : 2);
          KFactor_QCD_qqZZ_M    = kfactor_qqZZ_qcd_M(GenHMass, (sameflavor) ? 1 : 2, 2) / kfactor_qqZZ_qcd_M(GenHMass, (sameflavor) ? 1 : 2, 1);
          KFactor_QCD_qqZZ_Pt   = kfactor_qqZZ_qcd_Pt(GenHPt, (sameflavor) ? 1 : 2);
        }
        // Calculate NLO EWK K factors for qqZZ
        if (apply_K_NLOEW_ZZQQB){
          TLorentzVector GENZ1Vec, GENZ2Vec, GENZZVec;
          GENZ1Vec.SetPtEtaPhiM(GenZ1Pt, GenZ1Eta, GenZ1Phi, GenZ1Mass);
          GENZ2Vec.SetPtEtaPhiM(GenZ2Pt, GenZ2Eta, GenZ2Phi, GenZ2Mass);
          GENZZVec = GENZ1Vec + GENZ2Vec;

          KFactor_EW_qqZZ = EwkCorrections::getEwkCorrections(genParticles, ewkTable, genInfoP, GENZ1Vec, GENZ2Vec);

          bool sameflavor=(genZLeps.at(0)->pdgId()*genZLeps.at(1)->pdgId() == genZLeps.at(2)->pdgId()*genZLeps.at(3)->pdgId());
          float K_NNLO_LO = kfactor_qqZZ_qcd_M(GenHMass, (sameflavor) ? 1 : 2, 2);
          float rho = GENZZVec.Pt()/(GenLep1Pt+GenLep2Pt+GenLep3Pt+GenLep4Pt);
          if (rho<0.3) KFactor_EW_qqZZ_unc = fabs((K_NNLO_LO-1.)*(1.-KFactor_EW_qqZZ));
          else KFactor_EW_qqZZ_unc = fabs(1.-KFactor_EW_qqZZ);
        }
      }
    }
  }

}


void HZZ4lNtupleMaker::FillCandidate(const pat::CompositeCandidate& cand, bool evtPass, const edm::Event& event, Int_t CRFLAG)
{
  //Initialize a new candidate into the tree
  //myTree->createNewCandidate(); // this doesn't do anything anymore

  //Reinitialize the per-candidate vectors (necessary because in CRs we can store more than 1 candidate per event)
  LepPt.clear();
  LepEta.clear();
  LepPhi.clear();
  LepM.clear();
  LepSCEta.clear();
  LepLepId.clear();
  LepSIP.clear();
  Lepdxy.clear();
  Lepdz.clear();
  LepTime.clear();
  LepisID.clear();
  LepBDT.clear();
  LepisCrack.clear();
  LepMissingHit.clear();
  LepChargedHadIso.clear();
  LepNeutralHadIso.clear();
  LepPhotonIso.clear();
  LepPUIsoComponent.clear();
  LepCombRelIsoPF.clear();

  LepSF.clear();
  LepSF_Unc.clear();
	
  LepScale_Total_Up.clear();
  LepScale_Total_Dn.clear();
  LepScale_Stat_Up.clear();
  LepScale_Stat_Dn.clear();
  LepScale_Syst_Up.clear();
  LepScale_Syst_Dn.clear();
  LepScale_Gain_Up.clear();
  LepScale_Gain_Dn.clear();
  LepSigma_Total_Up.clear();
  LepSigma_Total_Dn.clear();
  LepSigma_Rho_Up.clear();
  LepSigma_Rho_Dn.clear();
  LepSigma_Phi_Up.clear();
  LepSigma_Phi_Dn.clear();

  HLTMatch1.clear();
  //HLTMatch2.clear();

  TauVSmu.clear();
  TauVSe.clear();
  TauVSjet.clear();
  TauDecayMode.clear();
  TauGenMatch.clear();
  TauTES_p_Up.clear();
  TauTES_p_Dn.clear();
  TauTES_m_Up.clear();
  TauTES_m_Dn.clear();
  TauTES_e_Up.clear();
  TauTES_e_Dn.clear();
  TauFES_p_Up.clear();
  TauFES_p_Dn.clear();
  TauFES_m_Up.clear();
  TauFES_m_Dn.clear();
  TauFES_e_Up.clear();
  TauFES_e_Dn.clear();

 
  TLE_dR_Z = -1;
  fsrPt.clear();
  fsrEta.clear();
  fsrPhi.clear();
  fsrLept.clear();
  fsrLeptID.clear();
  fsrDR.clear();
  fsrGenPt.clear();
  ExtraLepPt.clear();
  ExtraLepEta.clear();
  ExtraLepPhi.clear();
  ExtraLepLepId.clear();

  CRflag = CRFLAG;

  if(theChannel!=ZL){
    //Fill the info on the Higgs candidate
    ZZMass = cand.p4().mass();
    ZZMassErr = cand.userFloat("massError");
    ZZMassErrCorr = cand.userFloat("massErrorCorr");
    ZZMassPreFSR = cand.userFloat("m4l");

    ZZPt  = cand.p4().pt();
    ZZEta = cand.p4().eta();
    ZZPhi = cand.p4().phi();
    
    ZZjjPt = cand.userFloat("ZZjjPt");
    
//     ZZGoodMass = cand.userFloat("goodMass");
    eleHLTMatch = cand.userFloat("eleHLTMatch");
    muHLTMatch  = cand.userFloat("muHLTMatch");

    MET = cand.userFloat("MET");
    METPhi = cand.userFloat("METPhi");

    if(addVtxFit){
      ZZMassCFit = cand.userFloat("CFitM");
      ZZChi2CFit = cand.userFloat("CFitChi2");
    }

    DiJetMass  = cand.userFloat("DiJetMass");
    DiJetDEta  = cand.userFloat("DiJetDEta");

    //Fill the angular variables
    helcosthetaZ1 = cand.userFloat("costheta1");
    helcosthetaZ2 = cand.userFloat("costheta2");
    helphi       = cand.userFloat("phi");
    costhetastar = cand.userFloat("costhetastar");
    phistarZ1      = cand.userFloat("phistar1");
    //phistarZ2      = cand.userFloat("phistar2");
    xi            = cand.userFloat("xi");
    xistar        = cand.userFloat("xistar");

  }

  //Z1 and Z2 variables
  const reco::Candidate* Z1;
  const reco::Candidate* Z2;
  vector<const reco::Candidate*> leptons;
  vector<const reco::Candidate*> fsrPhot;
  vector<short> fsrIndex;
  vector<string> labels;

  if (theChannel!=ZL) { // Regular 4l candidates
    Z1   = cand.daughter("Z1");
    Z2   = cand.daughter("Z2");
    userdatahelpers::getSortedLeptons(cand, leptons, labels, fsrPhot, fsrIndex);
  } else {              // Special handling of Z+l candidates
    Z1   = cand.daughter(0); // the Z
    Z2   = cand.daughter(1); // This is actually the additional lepton!
    userdatahelpers::getSortedLeptons(cand, leptons, labels, fsrPhot, fsrIndex, false); // note: we get just 3 leptons in this case.
  }
  Z1Mass = Z1->mass();
  Z1Pt =   Z1->pt();
  Z1Eta =  Z1->eta();
  Z1Phi =  Z1->phi();
  Z1Flav =  getPdgId(Z1->daughter(0)) * getPdgId(Z1->daughter(1));
//   Z1GoodMass = userdatahelpers::getUserFloat(Z1,"goodMass");
  Z1muHLTMatch1 = userdatahelpers::getUserFloat(Z1,"muHLTMatch1");
  Z1muHLTMatch2 = userdatahelpers::getUserFloat(Z1,"muHLTMatch2");
  Z1muHLTMatch = userdatahelpers::getUserFloat(Z1,"muHLTMatch");
  Z1eleHLTMatch1 = userdatahelpers::getUserFloat(Z1,"eleHLTMatch1");
  Z1eleHLTMatch2 = userdatahelpers::getUserFloat(Z1,"eleHLTMatch2");
  Z1eleHLTMatch = userdatahelpers::getUserFloat(Z1,"eleHLTMatch");

  Z2Mass = Z2->mass();
  Z2Pt =   Z2->pt();
  Z2Eta  = Z2->eta();
  Z2Phi  = Z2->phi();
  Z2Flav = theChannel==ZL ? getPdgId(Z2) : getPdgId(Z2->daughter(0)) * getPdgId(Z2->daughter(1));
//   if(userdatahelpers::hasUserFloat(Z2,"goodMass")){
//     Z2GoodMass = userdatahelpers::getUserFloat(Z2,"goodMass");
//   }
  if(userdatahelpers::hasUserFloat(Z2,"muHLTMatch")){
    Z2muHLTMatch1 = userdatahelpers::getUserFloat(Z2,"muHLTMatch1");
    Z2muHLTMatch2 = userdatahelpers::getUserFloat(Z2,"muHLTMatch2");
    Z2muHLTMatch = userdatahelpers::getUserFloat(Z2,"muHLTMatch");
  }
  if(userdatahelpers::hasUserFloat(Z2,"eleHLTMatch")){
    Z2eleHLTMatch1 = userdatahelpers::getUserFloat(Z2,"eleHLTMatch1");
    Z2eleHLTMatch2 = userdatahelpers::getUserFloat(Z2,"eleHLTMatch2");
    Z2eleHLTMatch = userdatahelpers::getUserFloat(Z2,"eleHLTMatch");
  }
  if (theChannel!=ZL) 
    Z2isGoodTau=userdatahelpers::getUserFloat(Z2,"isGoodTau");
  else if (abs(Z2Flav)==15) {
    Z2isGoodTau_Ele=userdatahelpers::getUserFloat(Z2,"isGood_Ele");
    Z2isGoodTau_Mu=userdatahelpers::getUserFloat(Z2,"isGood_Mu");
    Z2isGoodTau_Tau=userdatahelpers::getUserFloat(Z2,"isGood_Tau");
  }

  if (theChannel!=ZL) {
//-------------------------------------------------------------------------------------------------------
//----------------------------------------SV FIT---------------------------------------------------------
//-------------------------------------------------------------------------------------------------------
  
  bool doSVFit=false;
  if (theChannel!=ZL && (abs(Z2Flav)==165 || abs(Z2Flav)==195 || abs(Z2Flav)==225 || abs(Z2Flav)==143)) doSVFit=true;
  bool swi;
  if (abs(leptons[2]->pdgId()) > abs(leptons[3]->pdgId())) swi=true;
  else if (abs(leptons[2]->pdgId()) < abs(leptons[3]->pdgId())) swi=false;
  else if (leptons[2]->pt() < leptons[3]->pt()) swi=true;
  else swi=false;
  int idx1,idx2;
  if (swi) {idx1=3;idx2=2;}
  else {idx1=2;idx2=3;}

  //central
  TLorentzVector Z2p4;
  Z2p4.SetPtEtaPhiM(Z2Pt,Z2Eta,Z2Phi,Z2Mass);
  TLorentzVector tau1,tau2,met;
  int pairType;
  int dm1,dm2;
  tau1.SetPxPyPzE(leptons[idx1]->px(),leptons[idx1]->py(),leptons[idx1]->pz(),leptons[idx1]->energy());
  tau2.SetPxPyPzE(leptons[idx2]->px(),leptons[idx2]->py(),leptons[idx2]->pz(),leptons[idx2]->energy());
  met.SetPxPyPzE(METx,METy,0,std::hypot(METx,METy));
  if (abs(Z2Flav)==195) pairType=0;
  else if (abs(Z2Flav)==165) pairType=1;
  else if (abs(Z2Flav)==225) pairType=2;
  else pairType=3;
  dm1=abs(leptons[idx1]->pdgId())==15?userdatahelpers::getUserFloat(leptons[idx1],"decayMode"):-1;
  dm2=abs(leptons[idx2]->pdgId())==15?userdatahelpers::getUserFloat(leptons[idx2],"decayMode"):-1;
  if (doSVFit) {
      SVfit algo_central(0,tau1,tau2,met,covMET,pairType,dm1,dm2);
      std::vector<double> results_central=algo_central.FitAndGetResult();
      Z2SVPt=results_central.at(0);
      Z2SVEta=results_central.at(1);
      Z2SVPhi=results_central.at(2);
      Z2SVMass=results_central.at(3);
  }

  //TES UP/DOWN
  TLorentzVector tau1_tesup,tau2_tesup,met_tesup;
  if (abs(leptons[idx1]->pdgId())==15 && userdatahelpers::getUserInt(leptons[idx1],"isTESShifted")) tau1_tesup.SetPxPyPzE(
      userdatahelpers::getUserFloat(leptons[idx1],"px_TauUp"),userdatahelpers::getUserFloat(leptons[idx1],"py_TauUp"),
      userdatahelpers::getUserFloat(leptons[idx1],"pz_TauUp"),userdatahelpers::getUserFloat(leptons[idx1],"e_TauUp"));
  else tau1_tesup.SetPxPyPzE(leptons[idx1]->px(),leptons[idx1]->py(),leptons[idx1]->pz(),leptons[idx1]->energy());
  if (abs(leptons[idx2]->pdgId())==15 && userdatahelpers::getUserInt(leptons[idx2],"isTESShifted")) tau2_tesup.SetPxPyPzE(
      userdatahelpers::getUserFloat(leptons[idx2],"px_TauUp"),userdatahelpers::getUserFloat(leptons[idx2],"py_TauUp"),
      userdatahelpers::getUserFloat(leptons[idx2],"pz_TauUp"),userdatahelpers::getUserFloat(leptons[idx2],"e_TauUp"));
  else tau2_tesup.SetPxPyPzE(leptons[idx2]->px(),leptons[idx2]->py(),leptons[idx2]->pz(),leptons[idx2]->energy());
  Z2Mass_TESup=(Z2p4-tau1+tau1_tesup-tau2+tau2_tesup).M();
  met_tesup.SetPxPyPzE(METxUPTES,METyUPTES,0,std::hypot(METxUPTES,METyUPTES));
  if (doSVFit) {
      if (METxUPTES==METx && METyUPTES==METy) {
          Z2SVPt_TESup=Z2SVPt;
          Z2SVEta_TESup=Z2SVEta;
          Z2SVPhi_TESup=Z2SVPhi;
          Z2SVMass_TESup=Z2SVMass;
      }
      else {
          SVfit algo_tesup(0,tau1_tesup,tau2_tesup,met_tesup,covMET,pairType,dm1,dm2);
          std::vector<double> results_tesup=algo_tesup.FitAndGetResult();
          Z2SVPt_TESup=results_tesup.at(0);
          Z2SVEta_TESup=results_tesup.at(1);
          Z2SVPhi_TESup=results_tesup.at(2);
          Z2SVMass_TESup=results_tesup.at(3);
      }
  }

  TLorentzVector tau1_tesdn,tau2_tesdn,met_tesdn;
  if (abs(leptons[idx1]->pdgId())==15 && userdatahelpers::getUserInt(leptons[idx1],"isTESShifted")) tau1_tesdn.SetPxPyPzE(
      userdatahelpers::getUserFloat(leptons[idx1],"px_TauDown"),userdatahelpers::getUserFloat(leptons[idx1],"py_TauDown"),
      userdatahelpers::getUserFloat(leptons[idx1],"pz_TauDown"),userdatahelpers::getUserFloat(leptons[idx1],"e_TauDown"));
  else tau1_tesdn.SetPxPyPzE(leptons[idx1]->px(),leptons[idx1]->py(),leptons[idx1]->pz(),leptons[idx1]->energy());
  if (abs(leptons[idx2]->pdgId())==15 && userdatahelpers::getUserInt(leptons[idx2],"isTESShifted")) tau2_tesdn.SetPxPyPzE(
      userdatahelpers::getUserFloat(leptons[idx2],"px_TauDown"),userdatahelpers::getUserFloat(leptons[idx2],"py_TauDown"),
      userdatahelpers::getUserFloat(leptons[idx2],"pz_TauDown"),userdatahelpers::getUserFloat(leptons[idx2],"e_TauDown"));
  else tau2_tesdn.SetPxPyPzE(leptons[idx2]->px(),leptons[idx2]->py(),leptons[idx2]->pz(),leptons[idx2]->energy());
  Z2Mass_TESdn=(Z2p4-tau1+tau1_tesdn-tau2+tau2_tesdn).M();
  met_tesdn.SetPxPyPzE(METxDOWNTES,METyDOWNTES,0,std::hypot(METxDOWNTES,METyDOWNTES));
  if (doSVFit) {
      if (METxDOWNTES==METx && METyDOWNTES==METy) {
          Z2SVPt_TESdn=Z2SVPt;
          Z2SVEta_TESdn=Z2SVEta;
          Z2SVPhi_TESdn=Z2SVPhi;
          Z2SVMass_TESdn=Z2SVMass;
      }
      else {
          SVfit algo_tesdn(0,tau1_tesdn,tau2_tesdn,met_tesdn,covMET,pairType,dm1,dm2);
          std::vector<double> results_tesdn=algo_tesdn.FitAndGetResult();
          Z2SVPt_TESdn=results_tesdn.at(0);
          Z2SVEta_TESdn=results_tesdn.at(1);
          Z2SVPhi_TESdn=results_tesdn.at(2);
          Z2SVMass_TESdn=results_tesdn.at(3);
      }
  }

  //EES UP/DOWN
  TLorentzVector tau1_eesup,tau2_eesup,met_eesup;
  if (abs(leptons[idx1]->pdgId())==15 && userdatahelpers::getUserInt(leptons[idx1],"isEESShifted")) tau1_eesup.SetPxPyPzE(
      userdatahelpers::getUserFloat(leptons[idx1],"px_EleUp"),userdatahelpers::getUserFloat(leptons[idx1],"py_EleUp"),
      userdatahelpers::getUserFloat(leptons[idx1],"pz_EleUp"),userdatahelpers::getUserFloat(leptons[idx1],"e_EleUp"));
  else tau1_eesup.SetPxPyPzE(leptons[idx1]->px(),leptons[idx1]->py(),leptons[idx1]->pz(),leptons[idx1]->energy());
  if (abs(leptons[idx2]->pdgId())==15 && userdatahelpers::getUserInt(leptons[idx2],"isEESShifted")) tau2_eesup.SetPxPyPzE(
      userdatahelpers::getUserFloat(leptons[idx2],"px_EleUp"),userdatahelpers::getUserFloat(leptons[idx2],"py_EleUp"),
      userdatahelpers::getUserFloat(leptons[idx2],"pz_EleUp"),userdatahelpers::getUserFloat(leptons[idx2],"e_EleUp"));
  else tau2_eesup.SetPxPyPzE(leptons[idx2]->px(),leptons[idx2]->py(),leptons[idx2]->pz(),leptons[idx2]->energy());
  Z2Mass_EESup=(Z2p4-tau1+tau1_eesup-tau2+tau2_eesup).M();
  met_eesup.SetPxPyPzE(METxUPEES,METyUPEES,0,std::hypot(METxUPEES,METyUPEES));
  if (doSVFit) {
      if (METxUPEES==METx && METyUPEES==METy) {
          Z2SVPt_EESup=Z2SVPt;
          Z2SVEta_EESup=Z2SVEta;
          Z2SVPhi_EESup=Z2SVPhi;
          Z2SVMass_EESup=Z2SVMass;
      }
      else {
          SVfit algo_eesup(0,tau1_eesup,tau2_eesup,met_eesup,covMET,pairType,dm1,dm2);
          std::vector<double> results_eesup=algo_eesup.FitAndGetResult();
          Z2SVPt_EESup=results_eesup.at(0);
          Z2SVEta_EESup=results_eesup.at(1);
          Z2SVPhi_EESup=results_eesup.at(2);
          Z2SVMass_EESup=results_eesup.at(3);
      }
  }

  TLorentzVector tau1_eesdn,tau2_eesdn,met_eesdn;
  if (abs(leptons[idx1]->pdgId())==15 && userdatahelpers::getUserInt(leptons[idx1],"isEESShifted")) tau1_eesdn.SetPxPyPzE(
      userdatahelpers::getUserFloat(leptons[idx1],"px_EleDown"),userdatahelpers::getUserFloat(leptons[idx1],"py_EleDown"),
      userdatahelpers::getUserFloat(leptons[idx1],"pz_EleDown"),userdatahelpers::getUserFloat(leptons[idx1],"e_EleDown"));
  else tau1_eesdn.SetPxPyPzE(leptons[idx1]->px(),leptons[idx1]->py(),leptons[idx1]->pz(),leptons[idx1]->energy());
  if (abs(leptons[idx2]->pdgId())==15 && userdatahelpers::getUserInt(leptons[idx2],"isEESShifted")) tau2_eesdn.SetPxPyPzE(
      userdatahelpers::getUserFloat(leptons[idx2],"px_EleDown"),userdatahelpers::getUserFloat(leptons[idx2],"py_EleDown"),
      userdatahelpers::getUserFloat(leptons[idx2],"pz_EleDown"),userdatahelpers::getUserFloat(leptons[idx2],"e_EleDown"));
  else tau2_eesdn.SetPxPyPzE(leptons[idx2]->px(),leptons[idx2]->py(),leptons[idx2]->pz(),leptons[idx2]->energy());
  Z2Mass_EESdn=(Z2p4-tau1+tau1_eesdn-tau2+tau2_eesdn).M();
  met_eesdn.SetPxPyPzE(METxDOWNEES,METyDOWNEES,0,std::hypot(METxDOWNEES,METyDOWNEES));
  if (doSVFit) {
      if (METxDOWNEES==METx && METyDOWNEES==METy) {
          Z2SVPt_EESdn=Z2SVPt;
          Z2SVEta_EESdn=Z2SVEta;
          Z2SVPhi_EESdn=Z2SVPhi;
          Z2SVMass_EESdn=Z2SVMass;
      }
      else {
          SVfit algo_eesdn(0,tau1_eesdn,tau2_eesdn,met_eesdn,covMET,pairType,dm1,dm2);
          std::vector<double> results_eesdn=algo_eesdn.FitAndGetResult();
          Z2SVPt_EESdn=results_eesdn.at(0);
          Z2SVEta_EESdn=results_eesdn.at(1);
          Z2SVPhi_EESdn=results_eesdn.at(2);
          Z2SVMass_EESdn=results_eesdn.at(3);
      }
  }

  //MES UP/DOWN
  TLorentzVector tau1_mesup,tau2_mesup,met_mesup;
  if (abs(leptons[idx1]->pdgId())==15 && (userdatahelpers::getUserFloat(leptons[idx1],"genmatch")==2 || userdatahelpers::getUserFloat(leptons[idx1],"genmatch")==4)) tau1_mesup.SetPxPyPzE(leptons[idx1]->px()*1.01,leptons[idx1]->py()*1.01,leptons[idx1]->pz()*1.01,leptons[idx1]->energy()*1.01);
  else tau1_mesup.SetPxPyPzE(leptons[idx1]->px(),leptons[idx1]->py(),leptons[idx1]->pz(),leptons[idx1]->energy());
  if (abs(leptons[idx2]->pdgId())==15 && (userdatahelpers::getUserFloat(leptons[idx2],"genmatch")==2 || userdatahelpers::getUserFloat(leptons[idx2],"genmatch")==4)) tau2_mesup.SetPxPyPzE(leptons[idx2]->px()*1.01,leptons[idx2]->py()*1.01,leptons[idx2]->pz()*1.01,leptons[idx2]->energy()*1.01);
  else tau2_mesup.SetPxPyPzE(leptons[idx2]->px(),leptons[idx2]->py(),leptons[idx2]->pz(),leptons[idx2]->energy());
  Z2Mass_MESup=(Z2p4-tau1+tau1_mesup-tau2+tau2_mesup).M();
  met_mesup.SetPxPyPzE(METxUPMES,METyUPMES,0,std::hypot(METxUPMES,METyUPMES));
  if (doSVFit) {
      if (METxUPMES==METx && METyUPMES==METy) {
          Z2SVPt_MESup=Z2SVPt;
          Z2SVEta_MESup=Z2SVEta;
          Z2SVPhi_MESup=Z2SVPhi;
          Z2SVMass_MESup=Z2SVMass;
      }
      else {
          SVfit algo_mesup(0,tau1_mesup,tau2_mesup,met_mesup,covMET,pairType,dm1,dm2);
          std::vector<double> results_mesup=algo_mesup.FitAndGetResult();
          Z2SVPt_MESup=results_mesup.at(0);
          Z2SVEta_MESup=results_mesup.at(1);
          Z2SVPhi_MESup=results_mesup.at(2);
          Z2SVMass_MESup=results_mesup.at(3);
      }
  }

  TLorentzVector tau1_mesdn,tau2_mesdn,met_mesdn;
  if (abs(leptons[idx1]->pdgId())==15 && (userdatahelpers::getUserFloat(leptons[idx1],"genmatch")==2 || userdatahelpers::getUserFloat(leptons[idx1],"genmatch")==4)) tau1_mesdn.SetPxPyPzE(leptons[idx1]->px()*0.99,leptons[idx1]->py()*0.99,leptons[idx1]->pz()*0.99,leptons[idx1]->energy()*0.99);
  else tau1_mesdn.SetPxPyPzE(leptons[idx1]->px(),leptons[idx1]->py(),leptons[idx1]->pz(),leptons[idx1]->energy());
  if (abs(leptons[idx2]->pdgId())==15 && (userdatahelpers::getUserFloat(leptons[idx2],"genmatch")==2 || userdatahelpers::getUserFloat(leptons[idx2],"genmatch")==4)) tau2_mesdn.SetPxPyPzE(leptons[idx2]->px()*0.99,leptons[idx2]->py()*0.99,leptons[idx2]->pz()*0.99,leptons[idx2]->energy()*0.99);
  else tau2_mesdn.SetPxPyPzE(leptons[idx2]->px(),leptons[idx2]->py(),leptons[idx2]->pz(),leptons[idx2]->energy());
  Z2Mass_MESdn=(Z2p4-tau1+tau1_mesdn-tau2+tau2_mesdn).M();
  met_mesdn.SetPxPyPzE(METxDOWNMES,METyDOWNMES,0,std::hypot(METxDOWNMES,METyDOWNMES));
  if (doSVFit) {
      if (METxDOWNMES==METx && METyDOWNMES==METy) {
          Z2SVPt_MESdn=Z2SVPt;
          Z2SVEta_MESdn=Z2SVEta;
          Z2SVPhi_MESdn=Z2SVPhi;
          Z2SVMass_MESdn=Z2SVMass;
      }
      else {
          SVfit algo_mesdn(0,tau1_mesdn,tau2_mesdn,met_mesdn,covMET,pairType,dm1,dm2);
          std::vector<double> results_mesdn=algo_mesdn.FitAndGetResult();
          Z2SVPt_MESdn=results_mesdn.at(0);
          Z2SVEta_MESdn=results_mesdn.at(1);
          Z2SVPhi_MESdn=results_mesdn.at(2);
          Z2SVMass_MESdn=results_mesdn.at(3);
      }
  }

  Z2Mass_JESup=Z2Mass_JESdn=Z2Mass_JERup=Z2Mass_JERdn=Z2Mass;
  if (doSVFit) {
      //JES UP/DOWN
      TLorentzVector tau1_jesup,tau2_jesup,met_jesup;
      tau1_jesup.SetPxPyPzE(leptons[idx1]->px(),leptons[idx1]->py(),leptons[idx1]->pz(),leptons[idx1]->energy());
      tau2_jesup.SetPxPyPzE(leptons[idx2]->px(),leptons[idx2]->py(),leptons[idx2]->pz(),leptons[idx2]->energy());
      met_jesup.SetPxPyPzE(METxUPJES,METyUPJES,0,std::hypot(METxUPJES,METyUPJES));
      if (METxUPJES==METx && METyUPJES==METy) {
          Z2SVPt_JESup=Z2SVPt;
          Z2SVEta_JESup=Z2SVEta;
          Z2SVPhi_JESup=Z2SVPhi;
          Z2SVMass_JESup=Z2SVMass;
      }
      else {
          SVfit algo_jesup(0,tau1_jesup,tau2_jesup,met_jesup,covMET,pairType,dm1,dm2);
          std::vector<double> results_jesup=algo_jesup.FitAndGetResult();
          Z2SVPt_JESup=results_jesup.at(0);
          Z2SVEta_JESup=results_jesup.at(1);
          Z2SVPhi_JESup=results_jesup.at(2);
          Z2SVMass_JESup=results_jesup.at(3);
      }

      TLorentzVector tau1_jesdn,tau2_jesdn,met_jesdn;
      tau1_jesdn.SetPxPyPzE(leptons[idx1]->px(),leptons[idx1]->py(),leptons[idx1]->pz(),leptons[idx1]->energy());
      tau2_jesdn.SetPxPyPzE(leptons[idx2]->px(),leptons[idx2]->py(),leptons[idx2]->pz(),leptons[idx2]->energy());
      met_jesdn.SetPxPyPzE(METxDOWNJES,METyDOWNJES,0,std::hypot(METxDOWNJES,METyDOWNJES));
      if (METxDOWNJES==METx && METyDOWNJES==METy) {
          Z2SVPt_JESdn=Z2SVPt;
          Z2SVEta_JESdn=Z2SVEta;
          Z2SVPhi_JESdn=Z2SVPhi;
          Z2SVMass_JESdn=Z2SVMass;
      }
      else {
          SVfit algo_jesdn(0,tau1_jesdn,tau2_jesdn,met_jesdn,covMET,pairType,dm1,dm2);
          std::vector<double> results_jesdn=algo_jesdn.FitAndGetResult();
          Z2SVPt_JESdn=results_jesdn.at(0);
          Z2SVEta_JESdn=results_jesdn.at(1);
          Z2SVPhi_JESdn=results_jesdn.at(2);
          Z2SVMass_JESdn=results_jesdn.at(3);
      }

      //JER UP/DOWN
      TLorentzVector tau1_jerup,tau2_jerup,met_jerup;
      tau1_jerup.SetPxPyPzE(leptons[idx1]->px(),leptons[idx1]->py(),leptons[idx1]->pz(),leptons[idx1]->energy());
      tau2_jerup.SetPxPyPzE(leptons[idx2]->px(),leptons[idx2]->py(),leptons[idx2]->pz(),leptons[idx2]->energy());
      met_jerup.SetPxPyPzE(METxUPJER,METyUPJER,0,std::hypot(METxUPJER,METyUPJER));
      if (METxUPJER==METx && METyUPJER==METy) {
          Z2SVPt_JERup=Z2SVPt;
          Z2SVEta_JERup=Z2SVEta;
          Z2SVPhi_JERup=Z2SVPhi;
          Z2SVMass_JERup=Z2SVMass;
      }
      else {
          SVfit algo_jerup(0,tau1_jerup,tau2_jerup,met_jerup,covMET,pairType,dm1,dm2);
          std::vector<double> results_jerup=algo_jerup.FitAndGetResult();
          Z2SVPt_JERup=results_jerup.at(0);
          Z2SVEta_JERup=results_jerup.at(1);
          Z2SVPhi_JERup=results_jerup.at(2);
          Z2SVMass_JERup=results_jerup.at(3);
      }

      TLorentzVector tau1_jerdn,tau2_jerdn,met_jerdn;
      tau1_jerdn.SetPxPyPzE(leptons[idx1]->px(),leptons[idx1]->py(),leptons[idx1]->pz(),leptons[idx1]->energy());
      tau2_jerdn.SetPxPyPzE(leptons[idx2]->px(),leptons[idx2]->py(),leptons[idx2]->pz(),leptons[idx2]->energy());
      met_jerdn.SetPxPyPzE(METxDOWNJER,METyDOWNJER,0,std::hypot(METxDOWNJER,METyDOWNJER));
      if (METxDOWNJER==METx && METyDOWNJER==METy) {
          Z2SVPt_JERdn=Z2SVPt;
          Z2SVEta_JERdn=Z2SVEta;
          Z2SVPhi_JERdn=Z2SVPhi;
          Z2SVMass_JERdn=Z2SVMass;
      }
      else {
          SVfit algo_jerdn(0,tau1_jerdn,tau2_jerdn,met_jerdn,covMET,pairType,dm1,dm2);
          std::vector<double> results_jerdn=algo_jerdn.FitAndGetResult();
          Z2SVPt_JERdn=results_jerdn.at(0);
          Z2SVEta_JERdn=results_jerdn.at(1);
          Z2SVPhi_JERdn=results_jerdn.at(2);
          Z2SVMass_JERdn=results_jerdn.at(3);
      }
  }
  
  Z2GoodMass=(Z2SVMass>0?Z2SVMass:Z2Mass);
  Z2GoodMass_TESup=(Z2SVMass_TESup>0?Z2SVMass_TESup:Z2Mass_TESup);
  Z2GoodMass_TESdn=(Z2SVMass_TESdn>0?Z2SVMass_TESdn:Z2Mass_TESdn);
  Z2GoodMass_EESup=(Z2SVMass_EESup>0?Z2SVMass_EESup:Z2Mass_EESup);
  Z2GoodMass_EESdn=(Z2SVMass_EESdn>0?Z2SVMass_EESdn:Z2Mass_EESdn);
  Z2GoodMass_MESup=(Z2SVMass_MESup>0?Z2SVMass_MESup:Z2Mass_MESup);
  Z2GoodMass_MESdn=(Z2SVMass_MESdn>0?Z2SVMass_MESdn:Z2Mass_MESdn);
  Z2GoodMass_JESup=(Z2SVMass_JESup>0?Z2SVMass_JESup:Z2Mass_JESup);
  Z2GoodMass_JESdn=(Z2SVMass_JESdn>0?Z2SVMass_JESdn:Z2Mass_JESdn);
  Z2GoodMass_JERup=(Z2SVMass_JERup>0?Z2SVMass_JERup:Z2Mass_JERup);
  Z2GoodMass_JERdn=(Z2SVMass_JERdn>0?Z2SVMass_JERdn:Z2Mass_JERdn);

//-------------------------------------------------------------------------------------------------------
//----------------------------------------END SV FIT-----------------------------------------------------
//-------------------------------------------------------------------------------------------------------
  }
    
  const reco::Candidate* non_TLE_Z = nullptr;
  size_t TLE_index = 999;
  if(abs(Z1Flav) == 11*11 || abs(Z1Flav) == 13*13) non_TLE_Z = Z1;
  if(abs(Z2Flav) == 11*11 || abs(Z2Flav) == 13*13) non_TLE_Z = Z2;
  for (size_t i=0; i<leptons.size(); ++i){
    if(abs(leptons[i]->pdgId()) == 22) TLE_index = i;
  }
  if(TLE_index < 999 && non_TLE_Z != nullptr) {
    TLE_dR_Z = reco::deltaR(non_TLE_Z->p4(), leptons[TLE_index]->p4()); 
  }

  Int_t sel = 0;
  if(theChannel==ZZ){

    // Precomputed selections
    bool candPass70Z2Loose = cand.userFloat("Z2Mass") &&
                             cand.userFloat("MAllComb") &&
                             cand.userFloat("pt1")>20 && cand.userFloat("pt2")>10. &&
                             ZZMass>70.;
    bool candPassFullSel70 = cand.userFloat("SR");
    bool candPassFullSel   = cand.userFloat("FullSel");
    bool candIsBest = cand.userFloat("isBestCand");
    bool passMz_zz = (Z1Mass>60. && Z1Mass<120. && Z2Mass>60. && Z2Mass<120.);   //FIXME hardcoded cut

    if (candIsBest) {
      //    sel = 10; //FIXME see above
      if (candPass70Z2Loose) sel=70;
      if (candPassFullSel70){ // includes MZ2 > 12
        sel = 90;
        if (candPassFullSel){
          sel=100;
          if (passMz_zz) sel = 120;
        }
      }
    }
  } else if(theChannel==ZLL) sel = 20; 
  else if(theChannel==ZL) sel = 10;
 

  if (!(evtPass)) {sel = -sel;} // avoid confusion when we write events which do not pass trigger/skim

  // Retrieve the userFloat of the leptons in vectors ordered in the same way.
  vector<float> SIP(4);
  vector<float> combRelIsoPF(4);
  passIsoPreFSR = true;
  for (unsigned int i=0; i<leptons.size(); ++i){
    float curr_dR = 999;
    if(i != TLE_index && TLE_index < 999)
      curr_dR = reco::deltaR(leptons[i]->p4(), leptons[TLE_index]->p4());
    if(curr_dR < TLE_min_dR_3l) TLE_min_dR_3l = curr_dR;

    short lepFlav = std::abs(leptons[i]->pdgId());

    SIP[i]             = userdatahelpers::getUserFloat(leptons[i],"SIP");
    passIsoPreFSR      = passIsoPreFSR&&(userdatahelpers::getUserFloat(leptons[i],"combRelIsoPF")<LeptonIsoHelper::isoCut(leptons[i]));

    //in the Legacy approach,  FSR-corrected iso is attached to the Z, not to the lepton!
    if (theChannel!=ZL) {
      combRelIsoPF[i]    = cand.userFloat(labels[i]+"combRelIsoPFFSRCorr"); // Note: the
      assert(SIP[i] == cand.userFloat(labels[i]+"SIP")); // Check that I don't mess up with labels[] and leptons[]
    } else {
      //FIXME: at the moment,  FSR-corrected iso is not computed for Z+L CRs
      combRelIsoPF[i]    = userdatahelpers::getUserFloat(leptons[i],"combRelIsoPF");
    }

    //Fill the info on the lepton candidates
    LepPt .push_back( leptons[i]->pt() );
    LepEta.push_back( leptons[i]->eta() );
    LepPhi.push_back( leptons[i]->phi() );
    LepM.push_back( leptons[i]->mass() );
    LepSCEta.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"SCeta") : -99. );
    int id =  leptons[i]->pdgId();
    if(id == 22 && (i == 1 || i == 3)) id=-22; //FIXME this assumes a standard ordering of leptons.
    LepLepId.push_back( id );
    LepSIP  .push_back( SIP[i] );
    Lepdxy  .push_back( userdatahelpers::getUserFloat(leptons[i],"dxy") );
    Lepdz   .push_back( userdatahelpers::getUserFloat(leptons[i],"dz") );
    LepTime .push_back( lepFlav==13 ? userdatahelpers::getUserFloat(leptons[i],"time") : 0. );
    LepisID .push_back( userdatahelpers::getUserFloat(leptons[i],"ID") );
    LepBDT  .push_back( (lepFlav==13 || lepFlav==11) ? userdatahelpers::getUserFloat(leptons[i],"BDT") : -99. );
    LepisCrack.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"isCrack") : 0 );
    LepMissingHit.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"missingHit") : 0 );
    if (userdatahelpers::hasUserFloat(leptons[i],"scale_total_up")) { // These are not set when APPLYMUCORR=false
      LepScale_Total_Up.push_back( userdatahelpers::getUserFloat(leptons[i],"scale_total_up") );
      LepScale_Total_Dn.push_back( userdatahelpers::getUserFloat(leptons[i],"scale_total_dn") );
      LepSigma_Total_Up.push_back( userdatahelpers::getUserFloat(leptons[i],"sigma_total_up") );
      LepSigma_Total_Dn.push_back( userdatahelpers::getUserFloat(leptons[i],"sigma_total_dn") );
    }
    LepScale_Stat_Up.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"scale_stat_up") : -99. );
    LepScale_Stat_Dn.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"scale_stat_dn") : -99. );
    LepScale_Syst_Up.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"scale_syst_up") : -99. );
    LepScale_Syst_Dn.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"scale_syst_dn") : -99. );
    LepScale_Gain_Up.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"scale_gain_up") : -99. );
    LepScale_Gain_Dn.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"scale_gain_dn") : -99. );
    LepSigma_Rho_Up.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"sigma_rho_up") : -99. );
    LepSigma_Rho_Dn.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"sigma_rho_dn") : -99. );
    LepSigma_Phi_Up.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"sigma_phi_up") : -99. );
    LepSigma_Phi_Dn.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"sigma_phi_dn") : -99. );
    LepChargedHadIso.push_back( userdatahelpers::getUserFloat(leptons[i],"PFChargedHadIso") );
    LepNeutralHadIso.push_back( userdatahelpers::getUserFloat(leptons[i],"PFNeutralHadIso") );
    LepPhotonIso.push_back( (lepFlav==13 || lepFlav==11) ? userdatahelpers::getUserFloat(leptons[i],"PFPhotonIso") : -99. );
    LepPUIsoComponent.push_back( lepFlav==13 ? userdatahelpers::getUserFloat(leptons[i],"PFPUChargedHadIso") : 0. );
    LepCombRelIsoPF.push_back( combRelIsoPF[i] );
    LepisLoose.push_back(userdatahelpers::hasUserFloat(leptons[i],"isLoose") == 1 ? userdatahelpers::getUserFloat(leptons[i],"isLoose") : -2);

    if(lepFlav==15){
	short tauVSmu=0,tauVSe=0,tauVSjet=0;
	for (size_t iDisc = 0; iDisc<VSmuDisc.size(); ++iDisc) 
	    if (userdatahelpers::getUserInt(leptons[i],VSmuDisc[iDisc].c_str())==1) tauVSmu=iDisc+1;
	for (size_t iDisc = 0; iDisc<VSeDisc.size(); ++iDisc)
	    if (userdatahelpers::getUserInt(leptons[i],VSeDisc[iDisc].c_str())==1) tauVSe=iDisc+1;
	for (size_t iDisc = 0; iDisc<VSjetDisc.size(); ++iDisc)
            if (userdatahelpers::getUserInt(leptons[i],VSjetDisc[iDisc].c_str())==1) tauVSjet=iDisc+1;
	TauVSmu.push_back(tauVSmu);
	TauVSe.push_back(tauVSe);
	TauVSjet.push_back(tauVSjet);

	TauDecayMode.push_back( userdatahelpers::getUserFloat(leptons[i],"decayMode") );
	TauGenMatch.push_back( userdatahelpers::getUserFloat(leptons[i],"genmatch") );
	if (userdatahelpers::getUserInt(leptons[i],"isTESShifted")){
	    TauTES_p_Up.push_back( userdatahelpers::getUserFloat(leptons[i],"px_TauUp")/leptons[i]->px() );
	    TauTES_p_Dn.push_back( userdatahelpers::getUserFloat(leptons[i],"px_TauDown")/leptons[i]->px() );
            TauTES_m_Up.push_back( userdatahelpers::getUserFloat(leptons[i],"m_TauUp")/leptons[i]->mass() );
            TauTES_m_Dn.push_back( userdatahelpers::getUserFloat(leptons[i],"m_TauDown")/leptons[i]->mass() );
            TauTES_e_Up.push_back( userdatahelpers::getUserFloat(leptons[i],"e_TauUp")/leptons[i]->energy() );
            TauTES_e_Dn.push_back( userdatahelpers::getUserFloat(leptons[i],"e_TauDown")/leptons[i]->energy() );
	} else{
	    TauTES_p_Up.push_back(0.);
            TauTES_p_Dn.push_back(0.);
            TauTES_m_Up.push_back(0.);
            TauTES_m_Dn.push_back(0.);
            TauTES_e_Up.push_back(0.);
            TauTES_e_Dn.push_back(0.);
	}
	if (userdatahelpers::getUserInt(leptons[i],"isEESShifted")){
	    TauFES_p_Up.push_back( userdatahelpers::getUserFloat(leptons[i],"px_EleUp")/leptons[i]->px() );
            TauFES_p_Dn.push_back( userdatahelpers::getUserFloat(leptons[i],"px_EleDown")/leptons[i]->px() );
            TauFES_m_Up.push_back( userdatahelpers::getUserFloat(leptons[i],"m_EleUp")/leptons[i]->mass() );
            TauFES_m_Dn.push_back( userdatahelpers::getUserFloat(leptons[i],"m_EleDown")/leptons[i]->mass() );
            TauFES_e_Up.push_back( userdatahelpers::getUserFloat(leptons[i],"e_EleUp")/leptons[i]->energy() );
            TauFES_e_Dn.push_back( userdatahelpers::getUserFloat(leptons[i],"e_EleDown")/leptons[i]->energy() );
	} else{
	    TauFES_p_Up.push_back(0.);
            TauFES_p_Dn.push_back(0.);
            TauFES_m_Up.push_back(0.);
            TauFES_m_Dn.push_back(0.);
            TauFES_e_Up.push_back(0.);
            TauFES_e_Dn.push_back(0.);
	}
    } else{
	TauVSmu.push_back(-1);
	TauVSe.push_back(-1);
	TauVSjet.push_back(-1);
	TauDecayMode.push_back(-1);
	TauGenMatch.push_back(-1);
	TauTES_p_Up.push_back(0.);
        TauTES_p_Dn.push_back(0.);
        TauTES_m_Up.push_back(0.);
        TauTES_m_Dn.push_back(0.);
        TauTES_e_Up.push_back(0.);
        TauTES_e_Dn.push_back(0.);
	TauFES_p_Up.push_back(0.);
        TauFES_p_Dn.push_back(0.);
        TauFES_m_Up.push_back(0.);
        TauFES_m_Dn.push_back(0.);
        TauFES_e_Up.push_back(0.);
        TauFES_e_Dn.push_back(0.);
    }

    HLTMatch1.push_back( userdatahelpers::hasUserFloat(leptons[i],"HLTMatch1") ? userdatahelpers::getUserFloat(leptons[i],"HLTMatch1") : -1 );
    //HLTMatch2.push_back( userdatahelpers::hasUserFloat(leptons[i],"HLTMatch2") ? userdatahelpers::getUserFloat(leptons[i],"HLTMatch2") : -1 );

  }

  // FSR
  for (unsigned i=0; i<fsrPhot.size(); ++i) {
    math::XYZTLorentzVector fsr = fsrPhot[i]->p4();
    fsrPt.push_back(fsr.pt());
    fsrEta.push_back(fsr.eta());
    fsrPhi.push_back(fsr.phi());
    fsrLept.push_back(fsrIndex[i]+1);
    fsrLeptID.push_back(leptons[fsrIndex[i]]->pdgId());
    fsrDR.push_back(ROOT::Math::VectorUtil::DeltaR(leptons[fsrIndex[i]]->momentum(), fsrPhot[i]->momentum()));
    int igen = MCHistoryTools::fsrMatch(fsrPhot[i], genFSR);
    double dRGenVsReco = -1.;
    double genpT = -1.;

    if (igen>=0) {
      dRGenVsReco = ROOT::Math::VectorUtil::DeltaR(genFSR[igen]->momentum(), fsrPhot[i]->momentum());
//       pTGen = genFSR[igen]->pt();
//       etaGen = genFSR[igen]->eta();
//       phiGen = genFSR[igen]->phi();
      if (dRGenVsReco<0.3) {// matching cut - FIXME
        genpT=genFSR[igen]->pt();
      }
    }
    fsrGenPt.push_back(genpT);


  }

  //convention: 0 -> 4mu   1 -> 4e   2 -> 2mu2e
  if(CRFLAG){
    ZXFakeweight = 1.;
    for(int izx=0;izx<2;izx++)
      ZXFakeweight *= getFakeWeight(Z2->daughter(izx)->pt(),Z2->daughter(izx)->eta(),Z2->daughter(izx)->pdgId(),Z1->daughter(0)->pdgId());
  }
  ZZsel = sel;

  if (theChannel!=ZL) {

    //Fill the info on categorization
    nExtraLep = cand.userFloat("nExtraLep"); // Why is this still a float at this point?
    nExtraZ=cand.userFloat("nExtraZ");

    //Fill the info on the extra leptons
    for (int iExtraLep=1; iExtraLep<=(int)nExtraLep; iExtraLep++){
      TString extraString;
      extraString.Form("ExtraLep%d",iExtraLep);
      if (cand.hasUserCand(extraString.Data())){
        //for(int iextra=0;iextra<4;iextra++)varExtra[iextra].Prepend(extraString.Data());
        reco::CandidatePtr candPtr=cand.userCand(extraString.Data());
        ExtraLepPt.push_back(candPtr->pt());
        ExtraLepEta.push_back(candPtr->eta());
        ExtraLepPhi.push_back(candPtr->phi());
        ExtraLepLepId.push_back(candPtr->pdgId());
      }
    }

  }

  //Compute the data/MC weight
  dataMCWeight = 1.;
  //When the trigger is not applied in the MC, apply a trigger efficiency factor instead (FIXME: here hardcoding the efficiencies computed for ICHEP2016)
  trigEffWeight = 1.;
  if(isMC) {
    
    dataMCWeight = getAllWeight(leptons);
    
    if (applyTrigEffWeight){
      Int_t ZZFlav = abs(Z1Flav*Z2Flav);
      if(ZZFlav==121*121 || ZZFlav==121*242) //4e
	trigEffWeight = 0.992;
      if(ZZFlav==169*169) //4mu
	trigEffWeight = 0.996;
      if(ZZFlav==169*121 || ZZFlav==169*242) //2e2mu
	trigEffWeight = 0.995;
    }
  }
  //Store an overall event weight (which is normalized by gen_sumWeights)
  overallEventWeight = PUWeight * genHEPMCweight * dataMCWeight * trigEffWeight;

  /* // debug printout for weights
  cout<<"Event "<<event.id().run()<<":"<<event.luminosityBlock()<<":"<<event.id().event()<<endl;
  cout<<" pileup weight =         "<<PUWeight<<endl;
  cout<<" sign of gen. weight =   "<<genHEPMCweight/fabs(genHEPMCweight)<<endl;
  cout<<" lepton data/MC weight = "<<dataMCWeight<<endl;
  for(unsigned int i=0; i<leptons.size(); ++i)
    cout<<"   lepton ID="<<leptons[i]->pdgId()<<", pT="<<leptons[i]->pt()<<", weight="<<getAllWeight(leptons[i])<<endl;
  cout<<" trigger eff. weight =   "<<trigEffWeight<<endl;
  cout<<"product of all =         "<<overallEventWeight/fabs(genHEPMCweight)<<endl;
  cout<<endl;
  //*/

}


void HZZ4lNtupleMaker::FillGENCandidate(const pat::CompositeCandidate& cand){ //AT
  passedFiducialSelection_bbf = cand.userData<bool>("passedFiducial");
  if(cand.userFloat("event") == 91257) {
    cout << passedFiducialSelection_bbf << endl;
  }
}


void HZZ4lNtupleMaker::getCheckedUserFloat(const pat::CompositeCandidate& cand, const std::string& strval, Float_t& setval, Float_t defaultval){
  if (cand.hasUserFloat(strval)) setval = cand.userFloat(strval);
  else setval = defaultval;
}



// ------------ method called once each job just before starting event loop  ------------
void HZZ4lNtupleMaker::beginJob()
{
  edm::Service<TFileService> fs;
  TTree *candTree = fs->make<TTree>(theFileName,"Event Summary");
  TTree *candTree_failed = 0;
  if (failedTreeLevel)
    candTree_failed = fs->make<TTree>(theFileName+"_failed","Event Summary");
  myTree = new HZZ4lNtupleFactory(candTree, candTree_failed);
  const int nbins = 45;
  hCounter = fs->make<TH1F>("Counters", "Counters", nbins, 0., nbins);
  BookAllBranches();
}

// ------------ method called once each job just after ending the event loop  ------------
void HZZ4lNtupleMaker::endJob()
{
  hCounter->SetBinContent(0 ,gen_sumWeights); // also stored in bin 40
  hCounter->SetBinContent(1 ,Nevt_Gen-gen_BUGGY);
  hCounter->SetBinContent(2 ,gen_ZZ4mu);
  hCounter->SetBinContent(3 ,gen_ZZ4e);
  hCounter->SetBinContent(4 ,gen_ZZ2mu2e);
  hCounter->SetBinContent(5 ,gen_ZZ2l2tau);
  hCounter->SetBinContent(6 ,gen_ZZ4mu_EtaAcceptance);
  hCounter->SetBinContent(7 ,gen_ZZ4mu_LeptonAcceptance);
  hCounter->SetBinContent(8 ,gen_ZZ2emu2tau);
  hCounter->SetBinContent(30,gen_ZZ2mu2tau);
  hCounter->SetBinContent(31,gen_ZZ2e2tau);
  hCounter->SetBinContent(9 ,gen_ZZ4tau);
  hCounter->SetBinContent(10,gen_ZZ4e_EtaAcceptance);
  hCounter->SetBinContent(11,gen_ZZ4e_LeptonAcceptance);
  hCounter->SetBinContent(14,gen_ZZ2mu2e_EtaAcceptance);
  hCounter->SetBinContent(15,gen_ZZ2mu2e_LeptonAcceptance);
  hCounter->SetBinContent(19,gen_BUGGY);
  hCounter->SetBinContent(20,gen_Unknown);

  hCounter->SetBinContent(40,gen_sumWeights); // Also stored in underflow bin; added here for convenience
  hCounter->SetBinContent(41,gen_sumGenMCWeight);
  hCounter->SetBinContent(42,gen_sumPUWeight);

  TH1 *h[1] ={ hCounter };
  for (int i = 0; i < 1; i++) {
    h[i]->GetXaxis()->SetBinLabel(1 ,"Nevt_Gen");
    h[i]->GetXaxis()->SetBinLabel(2 ,"gen_ZZ4mu");
    h[i]->GetXaxis()->SetBinLabel(3 ,"gen_ZZ4e");
    h[i]->GetXaxis()->SetBinLabel(4 ,"gen_ZZ2mu2e");
    h[i]->GetXaxis()->SetBinLabel(5 ,"gen_ZZ2l2tau");
    h[i]->GetXaxis()->SetBinLabel(6 ,"gen_ZZ4mu_EtaAcceptance");
    h[i]->GetXaxis()->SetBinLabel(7 ,"gen_ZZ4mu_LeptonAcceptance");
    h[i]->GetXaxis()->SetBinLabel(8 ,"gen_ZZ2emu2tau");
    h[i]->GetXaxis()->SetBinLabel(30,"gen_ZZ2mu2tau");
    h[i]->GetXaxis()->SetBinLabel(31,"gen_ZZ2e2tau");
    h[i]->GetXaxis()->SetBinLabel(9 ,"gen_ZZ4tau");
    h[i]->GetXaxis()->SetBinLabel(10,"gen_ZZ4e_EtaAcceptance");
    h[i]->GetXaxis()->SetBinLabel(11,"gen_ZZ4e_LeptonAcceptance");
    h[i]->GetXaxis()->SetBinLabel(14,"gen_ZZ2mu2e_EtaAcceptance");
    h[i]->GetXaxis()->SetBinLabel(15,"gen_ZZ2mu2e_LeptonAcceptance");
    h[i]->GetXaxis()->SetBinLabel(19,"gen_BUGGY");
    h[i]->GetXaxis()->SetBinLabel(20,"gen_Unknown");

    h[i]->GetXaxis()->SetBinLabel(40,"gen_sumWeights");
    h[i]->GetXaxis()->SetBinLabel(41,"gen_sumGenMCWeight");
    h[i]->GetXaxis()->SetBinLabel(42,"gen_sumPUWeight");
  }

  delete myTree;

  return;
}

// ------------ method called when starting to processes a run  ------------
void HZZ4lNtupleMaker::beginRun(edm::Run const& iRun, edm::EventSetup const&)
{
  if (firstRun){
    //if (lheHandler){
    //  edm::Handle<LHERunInfoProduct> lhe_runinfo;
    //  iRun.getByLabel(edm::InputTag("externalLHEProducer"), lhe_runinfo);
    //  lheHandler->setHeaderFromRunInfo(&lhe_runinfo);
    //}
    firstRun=false;
  }
}

// ------------ method called when ending the processing of a run  ------------
void HZZ4lNtupleMaker::endRun(edm::Run const& iRun, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void HZZ4lNtupleMaker::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
  Nevt_Gen_lumiBlock = 0;
}

// ------------ method called when ending the processing of a luminosity block  ------------
void HZZ4lNtupleMaker::endLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iSetup)
{
  Float_t Nevt_preskim = -1.;
  edm::Handle<edm::MergeableCounter> preSkimCounter;
  if (iLumi.getByToken(preSkimToken, preSkimCounter)) { // Counter before skim. Does not exist for non-skimmed samples.
    Nevt_preskim = preSkimCounter->value;
    // We do not use a filtering skim for the time being; so this is just left as a check in case we need it again in the future.
    if (!std::uncaught_exception() && Nevt_preskim>=0.) assert(Nevt_preskim == Nevt_Gen_lumiBlock);
  }

  Nevt_Gen += Nevt_Gen_lumiBlock;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HZZ4lNtupleMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


Float_t HZZ4lNtupleMaker::getAllWeight(const vector<const reco::Candidate*>& leptons) 
{
  Float_t totWeight = 1.;
	
  for(unsigned int i=0; i<leptons.size(); ++i){ 
    Int_t   myLepID = abs(leptons[i]->pdgId());
    if (skipMuDataMCWeight&& myLepID==13) return 1.;
    if (skipEleDataMCWeight&& myLepID==11) return 1.;
    
    float SF = 1.0;
    float SF_Unc = 0.0;
  
    Float_t myLepPt = leptons[i]->pt();
    Float_t myLepEta = leptons[i]->eta();
     
    Float_t SCeta;
    if (myLepID == 11) SCeta = userdatahelpers::getUserFloat(leptons[i],"SCeta");
    else SCeta = myLepEta;
    
    Float_t mySCeta;
     
    // Deal with very rare cases when SCeta is out of 2.5 bonds
    if ( myLepEta <= 2.5 && SCeta >= 2.5) mySCeta = 2.49;
    else if ( myLepEta >= -2.5 && SCeta <= -2.5) mySCeta = -2.49;
    else mySCeta = SCeta;
     
    bool isCrack;
    if (myLepID == 11) isCrack = userdatahelpers::getUserFloat(leptons[i],"isCrack");
    else isCrack = false;
     
 
    SF = lepSFHelper->getSF(year,myLepID,myLepPt,myLepEta, mySCeta, isCrack);
    SF_Unc = lepSFHelper->getSFError(year,myLepID,myLepPt,myLepEta, mySCeta, isCrack);

    LepSF.push_back(SF);
    LepSF_Unc.push_back(SF_Unc);

    totWeight *= SF;
  } 

  return totWeight;
}


Float_t HZZ4lNtupleMaker::getHqTWeight(double mH, double genPt) const
{
  if (skipHqTWeight) return 1.;

  //cout<<"mH = "<<mH<<", genPt = "<<genPt<<endl;
  if (mH<400 || genPt>250) return 1.;

  double weight = 1.;

  const int masses[4] = {400,600,800,1000};
  double massDiff = 1000;
  int iMass = -1;
  for (int i=0; i<4; ++i){
    double massDiffTmp = std::fabs(mH-masses[i]);
    if (massDiffTmp<massDiff){
      massDiff = massDiffTmp;
      iMass = i;
    }
  }

  if (iMass>=0) {
    weight = h_weight->GetBinContent(h_weight->FindBin(genPt));
  }
  return weight;
}


// Added by CO
Float_t HZZ4lNtupleMaker::getFakeWeight(Float_t LepPt, Float_t LepEta, Int_t LepID, Int_t LepZ1ID)
{
  if (skipFakeWeight) return 1.;

  // year 0 = 2011
  // year 1 = 2012

  Float_t weight  = 1.;

  Int_t nHisto=0;

  Float_t myLepPt   = LepPt;
  Float_t myLepEta  = fabs(LepEta);
  Int_t   myLepID   = abs(LepID);
  Int_t   myZ1LepID = abs(LepZ1ID);

  //cout << " pt = " << myLepPt << " eta = " << myLepEta << " ZID = " << myZ1LepID << " LepID = " << myLepID << endl;

  //avoid to go out of the TH boundary
  if(myLepPt > 79.) myLepPt = 79.;
  if(myLepID==13)nHisto+=2;
  if(myZ1LepID==13)nHisto+=1;


  TString Z1flavor = "Z1ee";     if(myZ1LepID==13) Z1flavor = "Z1mumu";
  TString Z2flavor = "electron"; if(myLepID==13)   Z2flavor = "muon";
  TString histo_name = "eff_"+Z1flavor+"_plus_"+Z2flavor;

  //cout << " histo = " << histo_name << endl;

  weight = h_ZXWeight[nHisto]->GetBinContent(h_ZXWeight[nHisto]->GetXaxis()->FindBin(myLepPt), h_ZXWeight[nHisto]->GetYaxis()->FindBin(myLepEta));
  /*
  h_ZXWeight[0]=FileZXWeightEle->Get("eff_Z1ee_plus_electron")->Clone();
  h_ZXWeight[2]=FileZXWeightEle->Get("eff_Z1mumu_plus_electron")->Clone();

  h_ZXWeight[5]=FileZXWeightMuo->Get("eff_Z1ee_plus_muon")->Clone();
  h_ZXWeight[7]=FileZXWeightMuo->Get("eff_Z1mumu_plus_muon")->Clone();
*/
  return weight;

} // end of getFakeWeight
void HZZ4lNtupleMaker::FillZGenInfo(Short_t Z1Id, Short_t Z2Id,
                                    const math::XYZTLorentzVector pZ1, const math::XYZTLorentzVector pZ2)
{
  GenZ1Mass= pZ1.M();
  GenZ1Pt= pZ1.Pt();
  GenZ1Eta= pZ1.Eta();
  GenZ1Phi= pZ1.Phi();
  GenZ1Flav= Z1Id;

  GenZ2Mass= pZ2.M();
  GenZ2Pt= pZ2.Pt();
  GenZ2Eta= pZ2.Eta();
  GenZ2Phi= pZ2.Phi();
  GenZ2Flav= Z2Id;

  return;
}

void HZZ4lNtupleMaker::FillVisZGenInfo(Short_t Z1Id, Short_t Z2Id,
                                    const math::XYZTLorentzVector pZ1, const math::XYZTLorentzVector pZ2)
{
  GenVisZ1Mass= pZ1.M();
  GenVisZ1Pt= pZ1.Pt();
  GenVisZ1Eta= pZ1.Eta();
  GenVisZ1Phi= pZ1.Phi();
  GenVisZ1Flav= Z1Id;

  GenVisZ2Mass= pZ2.M();
  GenVisZ2Pt= pZ2.Pt();
  GenVisZ2Eta= pZ2.Eta();
  GenVisZ2Phi= pZ2.Phi();
  GenVisZ2Flav= Z2Id;

  return;
}


void HZZ4lNtupleMaker::FillLepGenInfo(Short_t Lep1Id, Short_t Lep2Id, Short_t Lep3Id, Short_t Lep4Id,
                                      const math::XYZTLorentzVector Lep1, const math::XYZTLorentzVector Lep2,
                                      const math::XYZTLorentzVector Lep3, const math::XYZTLorentzVector Lep4)
{
  GenLep1Pt=Lep1.Pt();
  GenLep1Eta=Lep1.Eta();
  GenLep1Phi=Lep1.Phi();
  GenLep1M=Lep1.M();
  GenLep1Id=Lep1Id;

  GenLep2Pt=Lep2.Pt();
  GenLep2Eta=Lep2.Eta();
  GenLep2Phi=Lep2.Phi();
  GenLep2M=Lep2.M();
  GenLep2Id=Lep2Id;

  GenLep3Pt=Lep3.Pt();
  GenLep3Eta=Lep3.Eta();
  GenLep3Phi=Lep3.Phi();
  GenLep3M=Lep3.M();
  GenLep3Id=Lep3Id;

  GenLep4Pt=Lep4.Pt();
  GenLep4Eta=Lep4.Eta();
  GenLep4Phi=Lep4.Phi();
  GenLep4M=Lep4.M();
  GenLep4Id=Lep4Id;

//   if (addGenAngles) {
//     TUtil::computeAngles(Gencosthetastar, GenhelcosthetaZ1, GenhelcosthetaZ2, Genhelphi, GenphistarZ1, zzanalysis::tlv(Lep1), Lep1Id, zzanalysis::tlv(Lep2), Lep2Id, zzanalysis::tlv(Lep3), Lep3Id, zzanalysis::tlv(Lep4), Lep4Id);
//   }
    
  return;
}

void HZZ4lNtupleMaker::FillVisLepGenInfo(Short_t Lep1Id, Short_t Lep2Id, Short_t Lep3Id, Short_t Lep4Id,
                                      const math::XYZTLorentzVector Lep1, const math::XYZTLorentzVector Lep2,
                                      const math::XYZTLorentzVector Lep3, const math::XYZTLorentzVector Lep4)
{
  GenVisLep1Pt=Lep1.Pt();
  GenVisLep1Eta=Lep1.Eta();
  GenVisLep1Phi=Lep1.Phi();
  GenVisLep1M=Lep1.M();
  GenVisLep1Id=Lep1Id;

  GenVisLep2Pt=Lep2.Pt();
  GenVisLep2Eta=Lep2.Eta();
  GenVisLep2Phi=Lep2.Phi();
  GenVisLep2M=Lep2.M();
  GenVisLep2Id=Lep2Id;

  GenVisLep3Pt=Lep3.Pt();
  GenVisLep3Eta=Lep3.Eta();
  GenVisLep3Phi=Lep3.Phi();
  GenVisLep3M=Lep3.M();
  GenVisLep3Id=Lep3Id;

  GenVisLep4Pt=Lep4.Pt();
  GenVisLep4Eta=Lep4.Eta();
  GenVisLep4Phi=Lep4.Phi();
  GenVisLep4M=Lep4.M();
  GenVisLep4Id=Lep4Id;

  return;
}

void HZZ4lNtupleMaker::FillAssocLepGenInfo(std::vector<const reco::Candidate *>& AssocLeps)
{
  if (AssocLeps.size() >= 1) {
    GenAssocLep1Pt =AssocLeps.at(0)->p4().Pt();
    GenAssocLep1Eta=AssocLeps.at(0)->p4().Eta();
    GenAssocLep1Phi=AssocLeps.at(0)->p4().Phi();
    GenAssocLep1Id =AssocLeps.at(0)->pdgId();
  }
  if (AssocLeps.size() >= 2) {
    GenAssocLep2Pt =AssocLeps.at(1)->p4().Pt();
    GenAssocLep2Eta=AssocLeps.at(1)->p4().Eta();
    GenAssocLep2Phi=AssocLeps.at(1)->p4().Phi();
    GenAssocLep2Id =AssocLeps.at(1)->pdgId();
  }

  return;
}


void HZZ4lNtupleMaker::FillHGenInfo(const math::XYZTLorentzVector pH, float w)
{
  GenHMass=pH.M();
  GenHPt=pH.Pt();
  GenHRapidity=pH.Rapidity();

  HqTMCweight=w;

  return;
}


void HZZ4lNtupleMaker::FillLepGenIso(float_t Lep1Iso, float_t Lep2Iso, float_t Lep3Iso, float_t Lep4Iso)
{
  GenLep1Iso = Lep1Iso;
  GenLep2Iso = Lep2Iso;
  GenLep3Iso = Lep3Iso;
  GenLep4Iso = Lep4Iso;

  return;
}//AT

void HZZ4lNtupleMaker::BookAllBranches(){
   //Event variables
  myTree->Book("RunNumber",RunNumber, failedTreeLevel >= minimalFailedTree);
  myTree->Book("EventNumber",EventNumber, failedTreeLevel >= minimalFailedTree);
  myTree->Book("LumiNumber",LumiNumber, failedTreeLevel >= minimalFailedTree);
  myTree->Book("NRecoMu",NRecoMu, failedTreeLevel >= fullFailedTree);
  myTree->Book("NRecoEle",NRecoEle, failedTreeLevel >= fullFailedTree);
  myTree->Book("Nvtx",Nvtx, failedTreeLevel >= fullFailedTree);
  myTree->Book("NObsInt",NObsInt, failedTreeLevel >= fullFailedTree);
  myTree->Book("NTrueInt",NTrueInt, failedTreeLevel >= fullFailedTree);

  myTree->Book("GenMET", GenMET, failedTreeLevel >= minimalFailedTree);
  myTree->Book("GenMETPhi", GenMETPhi, failedTreeLevel >= minimalFailedTree);
  myTree->Book("PFMET", PFMET, failedTreeLevel >= fullFailedTree);
  myTree->Book("PFMETPhi", PFMETPhi, failedTreeLevel >= fullFailedTree);
  myTree->Book("METx", METx, failedTreeLevel >= fullFailedTree);
  myTree->Book("METy", METy, failedTreeLevel >= fullFailedTree);
  myTree->Book("METxUPTES", METxUPTES, failedTreeLevel >= fullFailedTree);
  myTree->Book("METyUPTES", METyUPTES, failedTreeLevel >= fullFailedTree);
  myTree->Book("METxDOWNTES", METxDOWNTES, failedTreeLevel >= fullFailedTree);
  myTree->Book("METyDOWNTES", METyDOWNTES, failedTreeLevel >= fullFailedTree);
  myTree->Book("METxUPEES", METxUPEES, failedTreeLevel >= fullFailedTree);
  myTree->Book("METyUPEES", METyUPEES, failedTreeLevel >= fullFailedTree);
  myTree->Book("METxDOWNEES", METxDOWNEES, failedTreeLevel >= fullFailedTree);
  myTree->Book("METyDOWNEES", METyDOWNEES, failedTreeLevel >= fullFailedTree);
  myTree->Book("METxUPMES", METxUPMES, failedTreeLevel >= fullFailedTree);
  myTree->Book("METyUPMES", METyUPMES, failedTreeLevel >= fullFailedTree);
  myTree->Book("METxDOWNMES", METxDOWNMES, failedTreeLevel >= fullFailedTree);
  myTree->Book("METyDOWNMES", METyDOWNMES, failedTreeLevel >= fullFailedTree);
  myTree->Book("METxUPJES", METxUPJES, failedTreeLevel >= fullFailedTree);
  myTree->Book("METyUPJES", METyUPJES, failedTreeLevel >= fullFailedTree);
  myTree->Book("METxDOWNJES", METxDOWNJES, failedTreeLevel >= fullFailedTree);
  myTree->Book("METyDOWNJES", METyDOWNJES, failedTreeLevel >= fullFailedTree);
  myTree->Book("METxUPJER", METxUPJER, failedTreeLevel >= fullFailedTree);
  myTree->Book("METyUPJER", METyUPJER, failedTreeLevel >= fullFailedTree);
  myTree->Book("METxDOWNJER", METxDOWNJER, failedTreeLevel >= fullFailedTree);
  myTree->Book("METyDOWNJER", METyDOWNJER, failedTreeLevel >= fullFailedTree);
    
  myTree->Book("nCleanedJets",nCleanedJets, failedTreeLevel >= minimalFailedTree);
  myTree->Book("nCleanedJetsPt30",nCleanedJetsPt30, failedTreeLevel >= minimalFailedTree);
  myTree->Book("nCleanedJetsPt30_jesUp",nCleanedJetsPt30_jesUp, failedTreeLevel >= minimalFailedTree);
  myTree->Book("nCleanedJetsPt30_jesUp_Total",nCleanedJetsPt30_jesUp_Total, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30_jesUp_Abs",nCleanedJetsPt30_jesUp_Abs, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30_jesUp_Abs_year",nCleanedJetsPt30_jesUp_Abs_year, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30_jesUp_BBEC1",nCleanedJetsPt30_jesUp_BBEC1, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30_jesUp_BBEC1_year",nCleanedJetsPt30_jesUp_BBEC1_year, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30_jesUp_EC2",nCleanedJetsPt30_jesUp_EC2, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30_jesUp_EC2_year",nCleanedJetsPt30_jesUp_EC2_year, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30_jesUp_FlavQCD",nCleanedJetsPt30_jesUp_FlavQCD, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30_jesUp_HF",nCleanedJetsPt30_jesUp_HF, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30_jesUp_HF_year",nCleanedJetsPt30_jesUp_HF_year, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30_jesUp_RelBal",nCleanedJetsPt30_jesUp_RelBal, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30_jesUp_RelSample_year",nCleanedJetsPt30_jesUp_RelSample_year, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30_jesDn",nCleanedJetsPt30_jesDn, failedTreeLevel >= minimalFailedTree);
  myTree->Book("nCleanedJetsPt30_jesDn_Total",nCleanedJetsPt30_jesDn_Total, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30_jesDn_Abs",nCleanedJetsPt30_jesDn_Abs, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30_jesDn_Abs_year",nCleanedJetsPt30_jesDn_Abs_year, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30_jesDn_BBEC1",nCleanedJetsPt30_jesDn_BBEC1, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30_jesDn_BBEC1_year",nCleanedJetsPt30_jesDn_BBEC1_year, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30_jesDn_EC2",nCleanedJetsPt30_jesDn_EC2, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30_jesDn_EC2_year",nCleanedJetsPt30_jesDn_EC2_year, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30_jesDn_FlavQCD",nCleanedJetsPt30_jesDn_FlavQCD, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30_jesDn_HF",nCleanedJetsPt30_jesDn_HF, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30_jesDn_HF_year",nCleanedJetsPt30_jesDn_HF_year, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30_jesDn_RelBal",nCleanedJetsPt30_jesDn_RelBal, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30_jesDn_RelSample_year",nCleanedJetsPt30_jesDn_RelSample_year, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30_jerUp",nCleanedJetsPt30_jerUp, failedTreeLevel >= minimalFailedTree);
  myTree->Book("nCleanedJetsPt30_jerDn",nCleanedJetsPt30_jerDn, failedTreeLevel >= minimalFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged",nCleanedJetsPt30BTagged, failedTreeLevel >= minimalFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF",nCleanedJetsPt30BTagged_bTagSF, failedTreeLevel >= minimalFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesUp",nCleanedJetsPt30BTagged_bTagSF_jesUp, failedTreeLevel >= minimalFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesUp_Total",nCleanedJetsPt30BTagged_bTagSF_jesUp_Total, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesUp_Abs",nCleanedJetsPt30BTagged_bTagSF_jesUp_Abs, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesUp_Abs_year",nCleanedJetsPt30BTagged_bTagSF_jesUp_Abs_year, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesUp_BBEC1",nCleanedJetsPt30BTagged_bTagSF_jesUp_BBEC1, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesUp_BBEC1_year",nCleanedJetsPt30BTagged_bTagSF_jesUp_BBEC1_year, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesUp_EC2",nCleanedJetsPt30BTagged_bTagSF_jesUp_EC2, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesUp_EC2_year",nCleanedJetsPt30BTagged_bTagSF_jesUp_EC2_year, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesUp_FlavQCD",nCleanedJetsPt30BTagged_bTagSF_jesUp_FlavQCD, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesUp_HF",nCleanedJetsPt30BTagged_bTagSF_jesUp_HF, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesUp_HF_year",nCleanedJetsPt30BTagged_bTagSF_jesUp_HF_year, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesUp_RelBal",nCleanedJetsPt30BTagged_bTagSF_jesUp_RelBal, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesUp_RelSample_year",nCleanedJetsPt30BTagged_bTagSF_jesUp_RelSample_year, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesDn",nCleanedJetsPt30BTagged_bTagSF_jesDn, failedTreeLevel >= minimalFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesDn_Total",nCleanedJetsPt30BTagged_bTagSF_jesDn_Total, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesDn_Abs",nCleanedJetsPt30BTagged_bTagSF_jesDn_Abs, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesDn_Abs_year",nCleanedJetsPt30BTagged_bTagSF_jesDn_Abs_year, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesDn_BBEC1",nCleanedJetsPt30BTagged_bTagSF_jesDn_BBEC1, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesDn_BBEC1_year",nCleanedJetsPt30BTagged_bTagSF_jesDn_BBEC1_year, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesDn_EC2",nCleanedJetsPt30BTagged_bTagSF_jesDn_EC2, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesDn_EC2_year",nCleanedJetsPt30BTagged_bTagSF_jesDn_EC2_year, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesDn_FlavQCD",nCleanedJetsPt30BTagged_bTagSF_jesDn_FlavQCD, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesDn_HF",nCleanedJetsPt30BTagged_bTagSF_jesDn_HF, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesDn_HF_year",nCleanedJetsPt30BTagged_bTagSF_jesDn_HF_year, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesDn_RelBal",nCleanedJetsPt30BTagged_bTagSF_jesDn_RelBal, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesDn_RelSample_year",nCleanedJetsPt30BTagged_bTagSF_jesDn_RelSample_year, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jerUp",nCleanedJetsPt30BTagged_bTagSF_jerUp, failedTreeLevel >= minimalFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jerDn",nCleanedJetsPt30BTagged_bTagSF_jerDn, failedTreeLevel >= minimalFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSFUp",nCleanedJetsPt30BTagged_bTagSFUp, failedTreeLevel >= minimalFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSFDn",nCleanedJetsPt30BTagged_bTagSFDn, failedTreeLevel >= minimalFailedTree);
    
  myTree->Book("trigWord",trigWord, failedTreeLevel >= minimalFailedTree);
  myTree->Book("evtPassMETFilter",evtPassMETTrigger, failedTreeLevel >= minimalFailedTree);
  myTree->Book("ZZMass",ZZMass, false);
  myTree->Book("ZZMassErr",ZZMassErr, false);
  myTree->Book("ZZMassErrCorr",ZZMassErrCorr, false);
  myTree->Book("ZZMassPreFSR",ZZMassPreFSR, false);
  myTree->Book("ZZsel",ZZsel, false);
  myTree->Book("ZZPt",ZZPt, false);
  myTree->Book("ZZEta",ZZEta, false);
  myTree->Book("ZZPhi",ZZPhi, false);
  myTree->Book("ZZjjPt",ZZjjPt, false);
//   myTree->Book("ZZGoodMass",ZZGoodMass, false);
  myTree->Book("eleHLTMatch",eleHLTMatch, false);
  myTree->Book("muHLTMatch",muHLTMatch, false);
  myTree->Book("CRflag",CRflag, false);
  myTree->Book("MET",MET, false);
  myTree->Book("METPhi",METPhi, false);

  //Kin refitted info
  if (addKinRefit) {
    myTree->Book("ZZMassRefit",ZZMassRefit, false);
    myTree->Book("ZZMassRefitErr",ZZMassRefitErr, false);
    myTree->Book("ZZMassUnrefitErr",ZZMassUnrefitErr, false);
  }
  if (addVtxFit){
    myTree->Book("ZZMassCFit",ZZMassCFit, false);
    myTree->Book("ZZChi2CFit",ZZChi2CFit, false);
  }
//   if (addZZKinfit){
//     myTree->Book("ZZKMass",ZZKMass, false);
//     myTree->Book("ZZKChi2",ZZKChi2, false);
//   }
//   if (addSVfit){
//     myTree->Book("ZZSVMass",ZZSVMass, false);
//     myTree->Book("ZZSVPt",ZZSVPt, false);
//     myTree->Book("ZZSVEta",ZZSVEta, false);
//     myTree->Book("ZZSVPhi",ZZSVPhi, false);
//   }

  //Z1 variables
  myTree->Book("Z1Mass",Z1Mass, false);
//   myTree->Book("Z1GoodMass",Z1GoodMass, false);
  myTree->Book("Z1Pt",Z1Pt, false);
  myTree->Book("Z1Eta",Z1Eta, false);
  myTree->Book("Z1Phi",Z1Phi, false);
  myTree->Book("Z1Flav",Z1Flav, false);
  myTree->Book("Z1muHLTMatch1",Z1muHLTMatch1, false);
  myTree->Book("Z1muHLTMatch2",Z1muHLTMatch2, false);
  myTree->Book("Z1muHLTMatch",Z1muHLTMatch, false);
  myTree->Book("Z1eleHLTMatch1",Z1eleHLTMatch1, false);
  myTree->Book("Z1eleHLTMatch2",Z1eleHLTMatch2, false);
  myTree->Book("Z1eleHLTMatch",Z1eleHLTMatch, false);

  //Z2 variables
  myTree->Book("Z2Mass",Z2Mass, false);
//   myTree->Book("Z2GoodMass",Z2GoodMass, false);
  myTree->Book("Z2Pt",Z2Pt, false);
  myTree->Book("Z2Eta",Z2Eta, false);
  myTree->Book("Z2Phi",Z2Phi, false);
  myTree->Book("Z2Flav",Z2Flav, false);
  if (theChannel!=ZL)
    myTree->Book("Z2isGoodTau",Z2isGoodTau, false);
  else if (abs(Z2Flav)==15) {
    myTree->Book("Z2isGoodTau_Ele",Z2isGoodTau_Ele, false);
    myTree->Book("Z2isGoodTau_Mu",Z2isGoodTau_Mu, false);
    myTree->Book("Z2isGoodTau_Tau",Z2isGoodTau_Tau, false);
  }
  if (theChannel!=ZL && addSVfit){
    myTree->Book("Z2SVMass",Z2SVMass, false);
    myTree->Book("Z2SVPt",Z2SVPt, false);
    myTree->Book("Z2SVEta",Z2SVEta, false);
    myTree->Book("Z2SVPhi",Z2SVPhi, false);
    myTree->Book("Z2GoodMass",Z2GoodMass, false);
      
    myTree->Book("Z2Mass_TESup",Z2Mass_TESup, false);
    myTree->Book("Z2SVMass_TESup",Z2SVMass_TESup, false);
    myTree->Book("Z2SVPt_TESup",Z2SVPt_TESup, false);
    myTree->Book("Z2SVEta_TESup",Z2SVEta_TESup, false);
    myTree->Book("Z2SVPhi_TESup",Z2SVPhi_TESup, false);
    myTree->Book("Z2GoodMass_TESup",Z2GoodMass_TESup, false);

    myTree->Book("Z2Mass_TESdn",Z2Mass_TESdn, false);
    myTree->Book("Z2SVMass_TESdn",Z2SVMass_TESdn, false);
    myTree->Book("Z2SVPt_TESdn",Z2SVPt_TESdn, false);
    myTree->Book("Z2SVEta_TESdn",Z2SVEta_TESdn, false);
    myTree->Book("Z2SVPhi_TESdn",Z2SVPhi_TESdn, false);
    myTree->Book("Z2GoodMass_TESdn",Z2GoodMass_TESdn, false);

    myTree->Book("Z2Mass_EESup",Z2Mass_EESup, false);
    myTree->Book("Z2SVMass_EESup",Z2SVMass_EESup, false);
    myTree->Book("Z2SVPt_EESup",Z2SVPt_EESup, false);
    myTree->Book("Z2SVEta_EESup",Z2SVEta_EESup, false);
    myTree->Book("Z2SVPhi_EESup",Z2SVPhi_EESup, false);
    myTree->Book("Z2GoodMass_EESup",Z2GoodMass_EESup, false);

    myTree->Book("Z2Mass_EESdn",Z2Mass_EESdn, false);
    myTree->Book("Z2SVMass_EESdn",Z2SVMass_EESdn, false);
    myTree->Book("Z2SVPt_EESdn",Z2SVPt_EESdn, false);
    myTree->Book("Z2SVEta_EESdn",Z2SVEta_EESdn, false);
    myTree->Book("Z2SVPhi_EESdn",Z2SVPhi_EESdn, false);
    myTree->Book("Z2GoodMass_EESdn",Z2GoodMass_EESdn, false);

    myTree->Book("Z2Mass_MESup",Z2Mass_MESup, false);
    myTree->Book("Z2SVMass_MESup",Z2SVMass_MESup, false);
    myTree->Book("Z2SVPt_MESup",Z2SVPt_MESup, false);
    myTree->Book("Z2SVEta_MESup",Z2SVEta_MESup, false);
    myTree->Book("Z2SVPhi_MESup",Z2SVPhi_MESup, false);
    myTree->Book("Z2GoodMass_MESup",Z2GoodMass_MESup, false);

    myTree->Book("Z2Mass_MESdn",Z2Mass_MESdn, false);
    myTree->Book("Z2SVMass_MESdn",Z2SVMass_MESdn, false);
    myTree->Book("Z2SVPt_MESdn",Z2SVPt_MESdn, false);
    myTree->Book("Z2SVEta_MESdn",Z2SVEta_MESdn, false);
    myTree->Book("Z2SVPhi_MESdn",Z2SVPhi_MESdn, false);
    myTree->Book("Z2GoodMass_MESdn",Z2GoodMass_MESdn, false);

    myTree->Book("Z2Mass_JESup",Z2Mass_JESup, false);
    myTree->Book("Z2SVMass_JESup",Z2SVMass_JESup, false);
    myTree->Book("Z2SVPt_JESup",Z2SVPt_JESup, false);
    myTree->Book("Z2SVEta_JESup",Z2SVEta_JESup, false);
    myTree->Book("Z2SVPhi_JESup",Z2SVPhi_JESup, false);
    myTree->Book("Z2GoodMass_JESup",Z2GoodMass_JESup, false);

    myTree->Book("Z2Mass_JESdn",Z2Mass_JESdn, false);
    myTree->Book("Z2SVMass_JESdn",Z2SVMass_JESdn, false);
    myTree->Book("Z2SVPt_JESdn",Z2SVPt_JESdn, false);
    myTree->Book("Z2SVEta_JESdn",Z2SVEta_JESdn, false);
    myTree->Book("Z2SVPhi_JESdn",Z2SVPhi_JESdn, false);
    myTree->Book("Z2GoodMass_JESdn",Z2GoodMass_JESdn, false);

    myTree->Book("Z2Mass_JERup",Z2Mass_JERup, false);
    myTree->Book("Z2SVMass_JERup",Z2SVMass_JERup, false);
    myTree->Book("Z2SVPt_JERup",Z2SVPt_JERup, false);
    myTree->Book("Z2SVEta_JERup",Z2SVEta_JERup, false);
    myTree->Book("Z2SVPhi_JERup",Z2SVPhi_JERup, false);
    myTree->Book("Z2GoodMass_JERup",Z2GoodMass_JERup, false);

    myTree->Book("Z2Mass_JERdn",Z2Mass_JERdn, false);
    myTree->Book("Z2SVMass_JERdn",Z2SVMass_JERdn, false);
    myTree->Book("Z2SVPt_JERdn",Z2SVPt_JERdn, false);
    myTree->Book("Z2SVEta_JERdn",Z2SVEta_JERdn, false);
    myTree->Book("Z2SVPhi_JERdn",Z2SVPhi_JERdn, false);
    myTree->Book("Z2GoodMass_JERdn",Z2GoodMass_JERdn, false);
  }
  myTree->Book("Z2muHLTMatch1",Z2muHLTMatch1, false);
  myTree->Book("Z2muHLTMatch2",Z2muHLTMatch2, false);
  myTree->Book("Z2muHLTMatch",Z2muHLTMatch, false);
  myTree->Book("Z2eleHLTMatch1",Z2eleHLTMatch1, false);
  myTree->Book("Z2eleHLTMatch2",Z2eleHLTMatch2, false);
  myTree->Book("Z2eleHLTMatch",Z2eleHLTMatch, false);

//   myTree->Book("costhetastar",costhetastar, false);
//   myTree->Book("helphi",helphi, false);
//   myTree->Book("helcosthetaZ1",helcosthetaZ1, false);
//   myTree->Book("helcosthetaZ2",helcosthetaZ2, false);
//   myTree->Book("phistarZ1",phistarZ1, false);
//   myTree->Book("phistarZ2",phistarZ2, false);
//   myTree->Book("xi",xi, false);
//   myTree->Book("xistar",xistar, false);

  if (is_loose_ele_selection) {
    myTree->Book("TLE_dR_Z",TLE_dR_Z, false);
    myTree->Book("TLE_min_dR_3l",TLE_min_dR_3l, false);
  }

  myTree->Book("LepPt",LepPt, false);
  myTree->Book("LepEta",LepEta, false);
  myTree->Book("LepPhi",LepPhi, false);
  myTree->Book("LepM",LepM, false);
  myTree->Book("LepSCEta",LepSCEta, false);
  myTree->Book("LepLepId",LepLepId, false);
  myTree->Book("LepSIP",LepSIP, false);
  myTree->Book("Lepdxy",Lepdxy, false);
  myTree->Book("Lepdz",Lepdz, false);
  myTree->Book("LepTime",LepTime, false);
  myTree->Book("LepisID",LepisID, false);
  myTree->Book("LepisLoose",LepisLoose, false);
  myTree->Book("LepBDT",LepBDT, false);
  myTree->Book("LepisCrack",LepisCrack, false);
  myTree->Book("LepMissingHit",LepMissingHit, false);
  myTree->Book("LepChargedHadIso",LepChargedHadIso, false);
  myTree->Book("LepNeutralHadIso",LepNeutralHadIso, false);
  myTree->Book("LepPhotonIso",LepPhotonIso, false);
  myTree->Book("LepPUIsoComponent",LepPUIsoComponent, false);
  myTree->Book("LepCombRelIsoPF",LepCombRelIsoPF, false);
  myTree->Book("LepSF",LepSF, false);
  myTree->Book("LepSF_Unc",LepSF_Unc, false);
  myTree->Book("LepScale_Total_Up",LepScale_Total_Up, false);
  myTree->Book("LepScale_Total_Dn",LepScale_Total_Dn, false);
  myTree->Book("LepScale_Stat_Up",LepScale_Stat_Up, false);
  myTree->Book("LepScale_Stat_Dn",LepScale_Stat_Dn, false);
  myTree->Book("LepScale_Syst_Up",LepScale_Syst_Up, false);
  myTree->Book("LepScale_Syst_Dn",LepScale_Syst_Dn, false);
  myTree->Book("LepScale_Gain_Up",LepScale_Gain_Up, false);
  myTree->Book("LepScale_Gain_Dn",LepScale_Gain_Dn, false);
  myTree->Book("LepSigma_Total_Up",LepSigma_Total_Up, false);
  myTree->Book("LepSigma_Total_Dn",LepSigma_Total_Dn, false);
  myTree->Book("LepSigma_Rho_Up",LepSigma_Rho_Up, false);
  myTree->Book("LepSigma_Rho_Dn",LepSigma_Rho_Dn, false);
  myTree->Book("LepSigma_Phi_Up",LepSigma_Phi_Up, false);
  myTree->Book("LepSigma_Phi_Dn",LepSigma_Phi_Up, false);

  myTree->Book("TauVSmu",TauVSmu, false);
  myTree->Book("TauVSe",TauVSe, false);
  myTree->Book("TauVSjet",TauVSjet, false);
  myTree->Book("TauDecayMode",TauDecayMode, false);
  myTree->Book("TauGenMatch",TauGenMatch, false);
  myTree->Book("TauTES_p_Up",TauTES_p_Up, false);
  myTree->Book("TauTES_p_Dn",TauTES_p_Dn, false);
  myTree->Book("TauTES_m_Up",TauTES_m_Up, false);
  myTree->Book("TauTES_m_Dn",TauTES_m_Dn, false);
  myTree->Book("TauTES_e_Up",TauTES_e_Up, false);
  myTree->Book("TauTES_e_Dn",TauTES_e_Dn, false);
  myTree->Book("TauFES_p_Up",TauFES_p_Up, false);
  myTree->Book("TauFES_p_Dn",TauFES_p_Dn, false);
  myTree->Book("TauFES_m_Up",TauFES_m_Up, false);
  myTree->Book("TauFES_m_Dn",TauFES_m_Dn, false);
  myTree->Book("TauFES_e_Up",TauFES_e_Up, false);
  myTree->Book("TauFES_e_Dn",TauFES_e_Dn, false);

  myTree->Book("HLTMatch1",HLTMatch1, false);
//  myTree->Book("HLTMatch2",HLTMatch2, false);

  myTree->Book("fsrPt",fsrPt, false);
  myTree->Book("fsrEta",fsrEta, false);
  myTree->Book("fsrPhi",fsrPhi, false);
  myTree->Book("fsrLept",fsrLept, false);
  myTree->Book("passIsoPreFSR",passIsoPreFSR, false);
  if (addFSRDetails) {
    myTree->Book("fsrDR",fsrDR, false);
    myTree->Book("fsrLeptId",fsrLeptID, false);
    myTree->Book("fsrGenPt",fsrGenPt, false);
  }

  //Jet variables
  myTree->Book("JetPt",JetPt, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetEta",JetEta, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPhi",JetPhi, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetMass",JetMass, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetEnergy",JetEnergy, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetBTagger",JetBTagger, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetIsBtagged",JetIsBtagged, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetIsBtaggedWithSF",JetIsBtaggedWithSF, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetIsBtaggedWithSFUp",JetIsBtaggedWithSFUp, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetIsBtaggedWithSFDn",JetIsBtaggedWithSFDn, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetQGLikelihood",JetQGLikelihood, failedTreeLevel >= fullFailedTree);
  if(addQGLInputs){
    myTree->Book("JetAxis2",JetAxis2, failedTreeLevel >= fullFailedTree);
    myTree->Book("JetMult",JetMult, failedTreeLevel >= fullFailedTree);
    myTree->Book("JetPtD",JetPtD, failedTreeLevel >= fullFailedTree);
  }
  myTree->Book("JetSigma",JetSigma, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetSigma_Total",JetSigma_Total, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetSigma_Abs",JetSigma_Abs, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetSigma_Abs_year",JetSigma_Abs_year, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetSigma_BBEC1",JetSigma_BBEC1, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetSigma_BBEC1_year",JetSigma_BBEC1_year, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetSigma_EC2",JetSigma_EC2, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetSigma_EC2_year",JetSigma_EC2_year, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetSigma_FlavQCD",JetSigma_FlavQCD, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetSigma_HF",JetSigma_HF, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetSigma_HF_year",JetSigma_HF_year, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetSigma_RelBal",JetSigma_RelBal, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetSigma_RelSample_year",JetSigma_RelSample_year, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetHadronFlavour",JetHadronFlavour, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPartonFlavour",JetPartonFlavour, failedTreeLevel >= fullFailedTree);

  myTree->Book("JetRawPt",JetRawPt, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPtJEC_noJER",JetPtJEC_noJER, failedTreeLevel >= fullFailedTree);
  

  myTree->Book("JetPt_JESUp",JetJESUp, failedTreeLevel >= minimalFailedTree);
  myTree->Book("JetPt_JESUp_Total",JetJESUp_Total, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPt_JESUp_Abs",JetJESUp_Abs, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPt_JESUp_Abs_year",JetJESUp_Abs_year, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPt_JESUp_BBEC1",JetJESUp_BBEC1, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPt_JESUp_BBEC1_year",JetJESUp_BBEC1_year, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPt_JESUp_EC2",JetJESUp_EC2, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPt_JESUp_EC2_year",JetJESUp_EC2_year, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPt_JESUp_FlavQCD",JetJESUp_FlavQCD, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPt_JESUp_HF",JetJESUp_HF, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPt_JESUp_HF_year",JetJESUp_HF_year, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPt_JESUp_RelBal",JetJESUp_RelBal, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPt_JESUp_RelSample_year",JetJESUp_RelSample_year, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPt_JESDown",JetJESDown, failedTreeLevel >= minimalFailedTree);
  myTree->Book("JetPt_JESDown_Total",JetJESDown_Total, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPt_JESDown_Abs",JetJESDown_Abs, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPt_JESDown_Abs_year",JetJESDown_Abs_year, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPt_JESDown_BBEC1",JetJESDown_BBEC1, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPt_JESDown_BBEC1_year",JetJESDown_BBEC1_year, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPt_JESDown_EC2",JetJESDown_EC2, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPt_JESDown_EC2_year",JetJESDown_EC2_year, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPt_JESDown_FlavQCD",JetJESDown_FlavQCD, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPt_JESDown_HF",JetJESDown_HF, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPt_JESDown_HF_year",JetJESDown_HF_year, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPt_JESDown_RelBal",JetJESDown_RelBal, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPt_JESDown_RelSample_year",JetJESDown_RelSample_year, failedTreeLevel >= fullFailedTree);
   
  myTree->Book("JetPt_JERUp",JetJERUp, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPt_JERDown",JetJERDown, failedTreeLevel >= fullFailedTree);

  myTree->Book("JetID", JetID, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPUID", JetPUID, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPUID_score", JetPUID_score, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPUValue", JetPUValue, failedTreeLevel >= fullFailedTree);

  myTree->Book("DiJetMass",DiJetMass, false);
  myTree->Book("DiJetDEta",DiJetDEta, false);
  
  //Photon variables
  myTree->Book("PhotonPt",PhotonPt, failedTreeLevel >= fullFailedTree);
  myTree->Book("PhotonEta",PhotonEta, failedTreeLevel >= fullFailedTree);
  myTree->Book("PhotonPhi",PhotonPhi, failedTreeLevel >= fullFailedTree);
  myTree->Book("PhotonIsCutBasedLooseID",PhotonIsCutBasedLooseID, failedTreeLevel >= fullFailedTree);
   
  myTree->Book("nExtraLep",nExtraLep, false);
  myTree->Book("nExtraZ",nExtraZ, false);
  myTree->Book("ExtraLepPt",ExtraLepPt, false);
  myTree->Book("ExtraLepEta",ExtraLepEta, false);
  myTree->Book("ExtraLepPhi",ExtraLepPhi, false);
  myTree->Book("ExtraLepLepId",ExtraLepLepId, false);

  myTree->Book("ZXFakeweight", ZXFakeweight, false);

  if (isMC){
    if (apply_K_NNLOQCD_ZZGG>0){
      myTree->Book("KFactor_QCD_ggZZ_Nominal", KFactor_QCD_ggZZ_Nominal, failedTreeLevel >= minimalFailedTree);
      myTree->Book("KFactor_QCD_ggZZ_PDFScaleDn", KFactor_QCD_ggZZ_PDFScaleDn, failedTreeLevel >= minimalFailedTree);
      myTree->Book("KFactor_QCD_ggZZ_PDFScaleUp", KFactor_QCD_ggZZ_PDFScaleUp, failedTreeLevel >= minimalFailedTree);
      myTree->Book("KFactor_QCD_ggZZ_QCDScaleDn", KFactor_QCD_ggZZ_QCDScaleDn, failedTreeLevel >= minimalFailedTree);
      myTree->Book("KFactor_QCD_ggZZ_QCDScaleUp", KFactor_QCD_ggZZ_QCDScaleUp, failedTreeLevel >= minimalFailedTree);
      myTree->Book("KFactor_QCD_ggZZ_AsDn", KFactor_QCD_ggZZ_AsDn, failedTreeLevel >= minimalFailedTree);
      myTree->Book("KFactor_QCD_ggZZ_AsUp", KFactor_QCD_ggZZ_AsUp, failedTreeLevel >= minimalFailedTree);
      myTree->Book("KFactor_QCD_ggZZ_PDFReplicaDn", KFactor_QCD_ggZZ_PDFReplicaDn, failedTreeLevel >= minimalFailedTree);
      myTree->Book("KFactor_QCD_ggZZ_PDFReplicaUp", KFactor_QCD_ggZZ_PDFReplicaUp, failedTreeLevel >= minimalFailedTree);
    }
    if (apply_K_NLOEW_ZZQQB){
      myTree->Book("KFactor_EW_qqZZ", KFactor_EW_qqZZ, failedTreeLevel >= minimalFailedTree);
      myTree->Book("KFactor_EW_qqZZ_unc", KFactor_EW_qqZZ_unc, failedTreeLevel >= minimalFailedTree);
    }
    if (apply_K_NNLOQCD_ZZQQB){
      myTree->Book("KFactor_QCD_qqZZ_dPhi", KFactor_QCD_qqZZ_dPhi, failedTreeLevel >= minimalFailedTree);
      myTree->Book("KFactor_QCD_qqZZ_M", KFactor_QCD_qqZZ_M, failedTreeLevel >= minimalFailedTree);
      myTree->Book("KFactor_QCD_qqZZ_Pt", KFactor_QCD_qqZZ_Pt, failedTreeLevel >= minimalFailedTree);
    }

    myTree->Book("genFinalState", genFinalState, failedTreeLevel >= minimalFailedTree);
    myTree->Book("genProcessId", genProcessId, failedTreeLevel >= minimalFailedTree);
    myTree->Book("genHEPMCweight", genHEPMCweight, failedTreeLevel >= minimalFailedTree);
    if (year == 2017 || year == 2018) myTree->Book("genHEPMCweight_NNLO", genHEPMCweight_NNLO, failedTreeLevel >= minimalFailedTree);
    //myTree->Book("genHEPMCweight_POWHEGonly", genHEPMCweight_POWHEGonly, failedTreeLevel >= minimalFailedTree);
    myTree->Book("PUWeight", PUWeight, failedTreeLevel >= minimalFailedTree);
    myTree->Book("PUWeight_Dn", PUWeight_Dn, failedTreeLevel >= minimalFailedTree);
    myTree->Book("PUWeight_Up", PUWeight_Up, failedTreeLevel >= minimalFailedTree);
    myTree->Book("dataMCWeight", dataMCWeight, false);
    myTree->Book("trigEffWeight", trigEffWeight, false);
    myTree->Book("overallEventWeight", overallEventWeight, false);
    myTree->Book("L1prefiringWeight", L1prefiringWeight, false);
    myTree->Book("L1prefiringWeightUp", L1prefiringWeightUp, false);
    myTree->Book("L1prefiringWeightDn", L1prefiringWeightDn, false);
    myTree->Book("HqTMCweight", HqTMCweight, failedTreeLevel >= minimalFailedTree);
    myTree->Book("xsec", xsection, failedTreeLevel >= minimalFailedTree);
    myTree->Book("genxsec", genxsection, failedTreeLevel >= minimalFailedTree);
    myTree->Book("genBR", genbranchingratio, failedTreeLevel >= minimalFailedTree);
    myTree->Book("genExtInfo", genExtInfo, failedTreeLevel >= minimalFailedTree);
    myTree->Book("GenHMass", GenHMass, failedTreeLevel >= minimalFailedTree);
    myTree->Book("GenHPt", GenHPt, failedTreeLevel >= minimalFailedTree);
    myTree->Book("GenHRapidity", GenHRapidity, failedTreeLevel >= minimalFailedTree);
    myTree->Book("GenZ1Mass", GenZ1Mass, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenZ1Pt", GenZ1Pt, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenZ1Phi", GenZ1Phi, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenZ1Eta", GenZ1Eta, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenZ1Flav", GenZ1Flav, failedTreeLevel >= minimalFailedTree);
    myTree->Book("GenZ2Mass", GenZ2Mass, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenZ2Pt", GenZ2Pt, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenZ2Phi", GenZ2Phi, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenZ2Eta", GenZ2Eta, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenZ2Flav", GenZ2Flav, failedTreeLevel >= minimalFailedTree);
    myTree->Book("GenLep1Pt", GenLep1Pt, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLep1Eta", GenLep1Eta, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLep1Phi", GenLep1Phi, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLep1M", GenLep1M, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLep1Id", GenLep1Id, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLep2Pt", GenLep2Pt, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLep2Eta", GenLep2Eta, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLep2Phi", GenLep2Phi, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLep2M", GenLep2M, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLep2Id", GenLep2Id, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLep3Pt", GenLep3Pt, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLep3Eta", GenLep3Eta, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLep3Phi", GenLep3Phi, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLep3M", GenLep3M, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLep3Id", GenLep3Id, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLep4Pt", GenLep4Pt, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLep4Eta", GenLep4Eta, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLep4Phi", GenLep4Phi, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLep4M", GenLep4M, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLep4Id", GenLep4Id, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisZ1Mass", GenVisZ1Mass, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisZ1Pt", GenVisZ1Pt, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisZ1Eta", GenVisZ1Eta, failedTreeLevel >= minimalFailedTree);
    myTree->Book("GenVisZ1Phi", GenVisZ1Phi, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisZ1Flav", GenVisZ1Flav, failedTreeLevel >= minimalFailedTree);
    myTree->Book("GenVisZ2Mass", GenVisZ2Mass, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisZ2Pt", GenVisZ2Pt, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisZ2Eta", GenVisZ2Eta, failedTreeLevel >= minimalFailedTree);
    myTree->Book("GenVisZ2Phi", GenVisZ2Phi, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisZ2Flav", GenVisZ2Flav, failedTreeLevel >= minimalFailedTree);
    myTree->Book("GenVisLep1Pt", GenVisLep1Pt, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisLep1Eta", GenVisLep1Eta, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisLep1Phi", GenVisLep1Phi, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisLep1M", GenVisLep1M, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisLep1Id", GenVisLep1Id, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisLep2Pt", GenVisLep2Pt, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisLep2Eta", GenVisLep2Eta, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisLep2Phi", GenVisLep2Phi, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisLep2M", GenVisLep2M, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisLep2Id", GenVisLep2Id, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisLep3Pt", GenVisLep3Pt, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisLep3Eta", GenVisLep3Eta, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisLep3Phi", GenVisLep3Phi, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisLep3M", GenVisLep3M, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisLep3Id", GenVisLep3Id, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisLep4Pt", GenVisLep4Pt, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisLep4Eta", GenVisLep4Eta, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisLep4Phi", GenVisLep4Phi, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisLep4M", GenVisLep4M, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisLep4Id", GenVisLep4Id, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenAssocLep1Pt", GenAssocLep1Pt, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenAssocLep1Eta", GenAssocLep1Eta, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenAssocLep1Phi", GenAssocLep1Phi, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenAssocLep1Id", GenAssocLep1Id, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenAssocLep2Pt", GenAssocLep2Pt, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenAssocLep2Eta", GenAssocLep2Eta, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenAssocLep2Phi", GenAssocLep2Phi, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenAssocLep2Id", GenAssocLep2Id, failedTreeLevel >= fullFailedTree);

    myTree->Book("GenLep1Iso", GenLep1Iso, failedTreeLevel >= minimalFailedTree); //AT
    myTree->Book("GenLep2Iso", GenLep2Iso, failedTreeLevel >= minimalFailedTree); //AT
    myTree->Book("GenLep3Iso", GenLep3Iso, failedTreeLevel >= minimalFailedTree); //AT
    myTree->Book("GenLep4Iso", GenLep4Iso, failedTreeLevel >= minimalFailedTree); //AT
//     if (addGenAngles) {
//       myTree->Book("Gencosthetastar", Gencosthetastar, failedTreeLevel >= minimalFailedTree); //AT
//       myTree->Book("GenhelcosthetaZ1", GenhelcosthetaZ1, failedTreeLevel >= minimalFailedTree); //AT
//       myTree->Book("GenhelcosthetaZ2", GenhelcosthetaZ2, failedTreeLevel >= minimalFailedTree); //AT
//       myTree->Book("Genhelphi", Genhelphi, failedTreeLevel >= minimalFailedTree); //AT
//       myTree->Book("GenphistarZ1", Genhelphi, failedTreeLevel >= minimalFailedTree); //AT
//     }
    myTree->Book("GenJetPt", GenJetPt, failedTreeLevel >= minimalFailedTree); //ATjets
    myTree->Book("GenJetMass", GenJetMass, failedTreeLevel >= minimalFailedTree); //ATjets
    myTree->Book("GenJetEta", GenJetEta, failedTreeLevel >= minimalFailedTree); //ATjets
    myTree->Book("GenJetPhi", GenJetPhi, failedTreeLevel >= minimalFailedTree); //ATjets
    myTree->Book("GenJetRapidity", GenJetRapidity, failedTreeLevel >= minimalFailedTree); //ATjets
    myTree->Book("nGenJet", nGenJet, failedTreeLevel >= minimalFailedTree); //ATjets
    myTree->Book("GenCleanedJetPt", GenCleanedJetPt, failedTreeLevel >= minimalFailedTree); //ATjets
    myTree->Book("GenCleanedJetMass", GenCleanedJetMass, failedTreeLevel >= minimalFailedTree); //ATjets
    myTree->Book("GenCleanedJetEta", GenCleanedJetEta, failedTreeLevel >= minimalFailedTree); //ATjets
    myTree->Book("GenCleanedJetPhi", GenCleanedJetPhi, failedTreeLevel >= minimalFailedTree); //ATjets
    myTree->Book("GenCleanedJetRapidity", GenCleanedJetRapidity, failedTreeLevel >= minimalFailedTree); //ATjets
    myTree->Book("GenCleanedJetHadronFlavour", GenCleanedJetHadronFlavour, failedTreeLevel >= minimalFailedTree); //ATjets
    myTree->Book("nCleanedGenJet", nCleanedGenJet, failedTreeLevel >= minimalFailedTree); //ATjets
      
    myTree->Book("htxs_errorCode", htxs_errorCode, failedTreeLevel >= minimalFailedTree);
    myTree->Book("htxs_prodMode", htxs_prodMode, failedTreeLevel >= minimalFailedTree);
    myTree->Book("htxsNJets", htxsNJets, failedTreeLevel >= minimalFailedTree);
    myTree->Book("htxsHPt", htxsHPt, failedTreeLevel >= minimalFailedTree);
    myTree->Book("htxs_stage0_cat", htxs_stage0_cat, failedTreeLevel >= minimalFailedTree);
    myTree->Book("htxs_stage1p1_cat", htxs_stage1p1_cat, failedTreeLevel >= minimalFailedTree);
    myTree->Book("htxs_stage1p2_cat", htxs_stage1p2_cat, failedTreeLevel >= minimalFailedTree);

    if(apply_QCD_GGF_UNCERT)
      {
	myTree->Book("ggH_NNLOPS_weight", ggH_NNLOPS_weight, failedTreeLevel >= minimalFailedTree);
	myTree->Book("ggH_NNLOPS_weight_unc", ggH_NNLOPS_weight_unc, failedTreeLevel >= minimalFailedTree);
	myTree->Book("qcd_ggF_uncertSF", qcd_ggF_uncertSF, failedTreeLevel >= minimalFailedTree);
      }	  

    myTree->Book("PythiaWeight_isr_muR4", PythiaWeight_isr_muR4, failedTreeLevel >= minimalFailedTree);
    myTree->Book("PythiaWeight_isr_muR2", PythiaWeight_isr_muR2, failedTreeLevel >= minimalFailedTree);
    myTree->Book("PythiaWeight_isr_muRsqrt2", PythiaWeight_isr_muRsqrt2, failedTreeLevel >= minimalFailedTree);
    myTree->Book("PythiaWeight_isr_muRoneoversqrt2", PythiaWeight_isr_muRoneoversqrt2, failedTreeLevel >= minimalFailedTree);
    myTree->Book("PythiaWeight_isr_muR0p5", PythiaWeight_isr_muR0p5, failedTreeLevel >= minimalFailedTree);
    myTree->Book("PythiaWeight_isr_muR0p25", PythiaWeight_isr_muR0p25, failedTreeLevel >= minimalFailedTree);
    myTree->Book("PythiaWeight_fsr_muR4", PythiaWeight_fsr_muR4, failedTreeLevel >= minimalFailedTree);
    myTree->Book("PythiaWeight_fsr_muR2", PythiaWeight_fsr_muR2, failedTreeLevel >= minimalFailedTree);
    myTree->Book("PythiaWeight_fsr_muRsqrt2", PythiaWeight_fsr_muRsqrt2, failedTreeLevel >= minimalFailedTree);
    myTree->Book("PythiaWeight_fsr_muRoneoversqrt2", PythiaWeight_fsr_muRoneoversqrt2, failedTreeLevel >= minimalFailedTree);
    myTree->Book("PythiaWeight_fsr_muR0p5", PythiaWeight_fsr_muR0p5, failedTreeLevel >= minimalFailedTree);
    myTree->Book("PythiaWeight_fsr_muR0p25", PythiaWeight_fsr_muR0p25, failedTreeLevel >= minimalFailedTree);
  }
// MELA branches are booked under buildMELA
}

void HZZ4lNtupleMaker::addweight(float &weight, float weighttoadd) {
  weight += weighttoadd;
}


Int_t HZZ4lNtupleMaker::FindBinValue(TGraphErrors *tgraph, double value)
{
   Double_t x_prev,x,y;
   Int_t bin = 0;
   x_prev = 0.;
   for(int i=0;i<tgraph->GetN();i++){
      tgraph->GetPoint(i,x,y);
      if(value > x_prev && value < x){
         bin = i;
         break;
      }
      else x_prev = x;
   }
   if (bin == 0) bin = 1;
   return bin-1;
}

//define this as a plug-in
DEFINE_FWK_MODULE(HZZ4lNtupleMaker);
