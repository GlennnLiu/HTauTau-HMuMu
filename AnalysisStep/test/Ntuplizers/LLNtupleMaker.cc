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
#include <CommonTools/Utils/interface/StringCutObjectSelector.h>
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

#include <TauPOG/TauIDSFs/interface/TauIDSFTool.h>

#include <TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h>
#include <TauAnalysis/ClassicSVfit/interface/ClassicSVfitIntegrand.h>
#include <TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h>
#include <TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h>

#include <HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h>
#include <HTT-utilities/RecoilCorrections/interface/MEtSys.h>

#include "LLConfigHelper.h"
#include "LLNtupleFactory.h"

#include <TRandom3.h>
#include <TH2D.h>
#include "TLorentzVector.h"
#include "TSpline.h"
#include "TGraphErrors.h"

#include <string>

bool verbose = false; //ATbbf

namespace {
  bool writeJets = true;     // Write jets in the tree. FIXME: make this configurable
  bool writePhotons = false; // Write photons in the tree. FIXME: make this configurable
  bool addSVfit = true;
  bool addFSRDetails = false;
  bool addQGLInputs = true;
  bool skipMuDataMCWeight = false; // skip computation of data/MC weight for mu 
  bool skipEleDataMCWeight = false; // skip computation of data/MC weight for ele
  bool skipTauDataMCWeight = false;

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
  Float_t GenMETx = -99;
  Float_t GenMETy = -99;

  Float_t PFMETRecoil = -99;
  Float_t PFMETPhiRecoil = -99;
  Float_t METxRecoil = -99;
  Float_t METyRecoil = -99;
  Float_t PFMET = -99;
  Float_t PFMETPhi = -99;
  Float_t METx = -99;
  Float_t METy = -99;

  TMatrixD covMET(2, 2);

  Float_t Pzeta1;
  Float_t Pzeta2;
  Float_t MtLMET;
//   Float_t METxUPTES = -99;
//   Float_t METyUPTES = -99;
//   Float_t METxDOWNTES = -99;
//   Float_t METyDOWNTES = -99;
//   Float_t METxUPEES = -99;
//   Float_t METyUPEES = -99;
//   Float_t METxDOWNEES = -99;
//   Float_t METyDOWNEES = -99;
//   Float_t METxUPMES = -99;
//   Float_t METyUPMES = -99;
//   Float_t METxDOWNMES = -99;
//   Float_t METyDOWNMES = -99;
//   Float_t METxUPJES = -99;
//   Float_t METyUPJES = -99;
//   Float_t METxDOWNJES = -99;
//   Float_t METyDOWNJES = -99;
//   Float_t METxUPJER = -99;
//   Float_t METyUPJER = -99;
//   Float_t METxDOWNJER = -99;
//   Float_t METyDOWNJER = -99;
    
  // MET with no HF
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
  Short_t nCleanedJetsPt25BTagged_bTagSF  = 0;
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
  Short_t metWord = 0;
 
  Float_t LLMass  = 0;
  Float_t LLPt  = 0;
  Float_t LLEta  = 0;
  Float_t LLPhi  = 0;
  Float_t LLDR = 0;
  Short_t LLFlav  = 0;
  Short_t LLisGoodTau = 0;
  Short_t LLisGoodTau_Ele = 0;
  Short_t LLisGoodTau_Mu = 0;
  Short_t LLisGoodTau_Tau = 0;


//-------------------------------------------------------------------
//---------------------SV Fit variables------------------------------
//-------------------------------------------------------------------

  Float_t LLSVMass = -99;
  Float_t LLSVPt = -99;
  Float_t LLSVEta = -99;
  Float_t LLSVPhi = -99;
  Float_t LLGoodMass = -99;

  std::vector<float> LLSVMass_up;
  std::vector<float> LLSVPt_up;
  std::vector<float> LLSVEta_up;
  std::vector<float> LLSVPhi_up;
  std::vector<float> LLMass_up;
  std::vector<float> LLGoodMass_up;
  std::vector<float> LLSVMass_dn;
  std::vector<float> LLSVPt_dn;
  std::vector<float> LLSVEta_dn;
  std::vector<float> LLSVPhi_dn;
  std::vector<float> LLMass_dn;
  std::vector<float> LLGoodMass_dn;
  // std::vector<std::string> Names_shift;
    
  //------------------------------End SV fit variables--------------------------------
    

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
  // std::vector<float> LepBDT;
  std::vector<bool> LepisCrack;
  std::vector<short> LepMissingHit;
  std::vector<bool> LepConversionVeto;
  std::vector<float> LepChargedHadIso;
  std::vector<float> LepNeutralHadIso;
  std::vector<float> LepPhotonIso;
  std::vector<float> LepPUIsoComponent;
  std::vector<float> LepCombRelIsoPF;
  std::vector<short> LepisLoose;

  std::vector<float> LepSF;
  std::vector<float> LepSF_UncUp;
  std::vector<float> LepSF_UncDn;
  std::vector<float> LepSF_UncUp_RECO_syst;
  std::vector<float> LepSF_UncUp_RECO_stat;
  std::vector<float> LepSF_UncUp_ID_syst;
  std::vector<float> LepSF_UncUp_ID_stat;
  std::vector<float> LepSF_UncUp_ISO_syst;
  std::vector<float> LepSF_UncUp_ISO_stat;
  std::vector<float> LepSF_UncDn_RECO_syst;
  std::vector<float> LepSF_UncDn_RECO_stat;
  std::vector<float> LepSF_UncDn_ID_syst;
  std::vector<float> LepSF_UncDn_ID_stat;
  std::vector<float> LepSF_UncDn_ISO_syst;
  std::vector<float> LepSF_UncDn_ISO_stat;
  std::vector<float> LepSF_UncUp_uncert0;
  std::vector<float> LepSF_UncUp_uncert1;
  std::vector<float> LepSF_UncUp_syst_alleras;
  std::vector<float> LepSF_UncUp_syst_year;
  std::vector<float> LepSF_UncUp_syst_dm_year;
  std::vector<float> LepSF_UncUp_fakeEle;
  std::vector<float> LepSF_UncUp_fakeMu;
  std::vector<float> LepSF_UncDn_uncert0;
  std::vector<float> LepSF_UncDn_uncert1;
  std::vector<float> LepSF_UncDn_syst_alleras;
  std::vector<float> LepSF_UncDn_syst_year;
  std::vector<float> LepSF_UncDn_syst_dm_year;
  std::vector<float> LepSF_UncDn_fakeEle;
  std::vector<float> LepSF_UncDn_fakeMu;

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

//HLT trigger match
  Short_t pass_SingleTrigger = 0;
  Short_t pass_CrossTrigger = 0;
  Short_t pass_Trigger = 0;

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

  //VBF jets
  Float_t DeltaEtaJJ = 0;
  Float_t DiJetMass = 0;
  Short_t VBFJetIdx1 = -1;
  Short_t VBFJetIdx2 = -1;
  
  // Photon info
  std::vector<float> PhotonPt ;
  std::vector<float> PhotonEta ;
  std::vector<float> PhotonPhi ;
  std::vector<bool> PhotonIsCutBasedLooseID;
   
  Short_t genFinalState  = 0;
  Float_t genHEPMCweight  = 0;

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
  Float_t trigEffWeight  = 1;
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

  Float_t GenLLMass  = 0;
  Float_t GenLLEta  = 0;
  Float_t GenLLPt  = 0;
  Float_t GenLLPhi  = 0;
  Float_t GenLLFlav  = 0;
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
  //Visible information
  Float_t GenVisLLMass  = 0;
  Float_t GenVisLLEta  = 0;
  Float_t GenVisLLPt  = 0;
  Float_t GenVisLLPhi  = 0;
  Float_t GenVisLLFlav  = 0;
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
  Float_t GenLep1Iso  = 0; //AT
  Float_t GenLep2Iso  = 0; //AT
    
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
class LLNtupleMaker : public edm::EDAnalyzer {
public:
  explicit LLNtupleMaker(const edm::ParameterSet&);
  ~LLNtupleMaker();

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
  virtual void FillKFactors(edm::Handle<GenEventInfoProduct>& genInfo, std::vector<const reco::Candidate *>& genZ, std::vector<const reco::Candidate *>& genLeps);
  virtual void FillCandidate(const pat::CompositeCandidate& higgs, bool evtPass, const edm::Event&);//, const Int_t CRflag);
  virtual void FillJet(const pat::Jet& jet);
  virtual void FillPhoton(int year, const pat::Photon& photon);
  virtual void endJob() ;
    

  void FillLLGenInfo(Short_t ZId, const math::XYZTLorentzVector pZ);
  void FillLepGenInfo(Short_t Lep1Id, Short_t Lep2Id, const math::XYZTLorentzVector Lep1, const math::XYZTLorentzVector Lep2);
  void FillVisLLGenInfo(Short_t ZId, const math::XYZTLorentzVector pZ);
  void FillVisLepGenInfo(Short_t Lep1Id, Short_t Lep2Id, const math::XYZTLorentzVector Lep1, const math::XYZTLorentzVector Lep2);
  void FillLepGenIso(float_t Lep1Iso, float_t Lep2Iso); //AT
  
  Float_t getAllWeight(const vector<const reco::Candidate*>& leptons);
  // Float_t getHqTWeight(double mH, double genPt) const;
  Int_t FindBinValue(TGraphErrors *tgraph, double value);

  void getCheckedUserFloat(const pat::CompositeCandidate& cand, const std::string& strval, Float_t& setval, Float_t defaultval=0);
	
	

  // ----------member data ---------------------------
  LLConfigHelper myHelper;
  int theChannel;
  std::string theCandLabel;
  TString theFileName;

  LLNtupleFactory *myTree;
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
  
  const StringCutObjectSelector<pat::CompositeCandidate, true> cut;


  int apply_K_NNLOQCD_ZZGG; // 0: Do not; 1: NNLO/LO; 2: NNLO/NLO; 3: NLO/LO
  bool apply_K_NNLOQCD_ZZQQB;
  bool apply_K_NLOEW_ZZQQB;
  bool apply_QCD_GGF_UNCERT;

  bool do_MET_Recoil;

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
  edm::EDGetTokenT<edm::TriggerResults> triggerResultPATToken;
  edm::EDGetTokenT<vector<reco::Vertex> > vtxToken;
  edm::EDGetTokenT<edm::View<pat::Jet> > jetToken;
  edm::EDGetTokenT<pat::PhotonCollection> photonToken;
  
  edm::EDGetTokenT<pat::METCollection> metToken;
  edm::EDGetTokenT<math::Error<2>::type> theCovTag;
  edm::EDGetTokenT<pat::METCollection> metJESupToken;
  edm::EDGetTokenT<pat::METCollection> metJESdnToken;
//   edm::EDGetTokenT<double> theMETdxUPTESTag;
//   edm::EDGetTokenT<double> theMETdyUPTESTag;
//   edm::EDGetTokenT<double> theMETdxDOWNTESTag;
//   edm::EDGetTokenT<double> theMETdyDOWNTESTag;
//   edm::EDGetTokenT<double> theMETdxUPEESTag;
//   edm::EDGetTokenT<double> theMETdyUPEESTag;
//   edm::EDGetTokenT<double> theMETdxDOWNEESTag;
//   edm::EDGetTokenT<double> theMETdyDOWNEESTag;
//   edm::EDGetTokenT<double> theMETdxUPMESTag;
//   edm::EDGetTokenT<double> theMETdyUPMESTag;
//   edm::EDGetTokenT<double> theMETdxDOWNMESTag;
//   edm::EDGetTokenT<double> theMETdyDOWNMESTag;
//   edm::EDGetTokenT<double> theMETdxUPJESTag;
//   edm::EDGetTokenT<double> theMETdyUPJESTag;
//   edm::EDGetTokenT<double> theMETdxDOWNJESTag;
//   edm::EDGetTokenT<double> theMETdyDOWNJESTag;
//   edm::EDGetTokenT<double> theMETdxUPJERTag;
//   edm::EDGetTokenT<double> theMETdyUPJERTag;
//   edm::EDGetTokenT<double> theMETdxDOWNJERTag;
//   edm::EDGetTokenT<double> theMETdyDOWNJERTag;
    
  edm::EDGetTokenT<pat::MuonCollection> muonToken;
  edm::EDGetTokenT<pat::ElectronCollection> electronToken;
  edm::EDGetTokenT<HTXS::HiggsClassification> htxsToken;
  edm::EDGetTokenT<edm::MergeableCounter> preSkimToken;
    
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

  Float_t gen_mumu,gen_mumu_LeptonAcceptance,gen_mumu_EtaAcceptance;
  Float_t gen_etau,gen_etau_LeptonAcceptance,gen_etau_EtaAcceptance;
  Float_t gen_mutau,gen_mutau_LeptonAcceptance,gen_mutau_EtaAcceptance;
  Float_t gen_tautau,gen_tautau_LeptonAcceptance,gen_tautau_EtaAcceptance;
  Float_t gen_emu,gen_emu_LeptonAcceptance,gen_emu_EtaAcceptance;
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
  TauIDSFTool *DeepTauSF_VSe_ETau, *DeepTauSF_VSmu_ETau, *DeepTauSF_VSjet_ETau, *DeepTauSF_VSe_MuTau, *DeepTauSF_VSmu_MuTau, *DeepTauSF_VSjet_MuTau, *DeepTauSF_VSe_TauTau, *DeepTauSF_VSmu_TauTau, *DeepTauSF_VSjet_TauTau;
  RecoilCorrector *recoilPFMetCorrector;
  MEtSys *recoilPFMetSyst;
  std::vector<string> uncSources {};


  TGraphErrors *gr_NNLOPSratio_pt_powheg_0jet;
  TGraphErrors *gr_NNLOPSratio_pt_powheg_1jet;
  TGraphErrors *gr_NNLOPSratio_pt_powheg_2jet;
  TGraphErrors *gr_NNLOPSratio_pt_powheg_3jet;
    
  bool firstRun;

};

//
// constructors and destructor
//
LLNtupleMaker::LLNtupleMaker(const edm::ParameterSet& pset) :
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
  cut(pset.getParameter<std::string>("cut")),

  //lheHandler(nullptr),
  apply_K_NNLOQCD_ZZGG(pset.getParameter<int>("Apply_K_NNLOQCD_ZZGG")),
  apply_K_NNLOQCD_ZZQQB(pset.getParameter<bool>("Apply_K_NNLOQCD_ZZQQB")),
  apply_K_NLOEW_ZZQQB(pset.getParameter<bool>("Apply_K_NLOEW_ZZQQB")),
  apply_QCD_GGF_UNCERT(pset.getParameter<bool>("Apply_QCD_GGF_UNCERT")),

  do_MET_Recoil(pset.getParameter<bool>("doMETRecoil")),

  pileUpReweight(nullptr),
  sampleName(pset.getParameter<string>("sampleName")),
  dataTag(pset.getParameter<string>("dataTag")),

  firstRun(true)
{
  //cout<< "Beginning Constructor\n\n\n" <<endl;
  consumesMany<std::vector< PileupSummaryInfo > >();
  genParticleToken = consumes<edm::View<reco::Candidate> >(edm::InputTag("prunedGenParticles"));
  genParticleToken_bbf = consumes<reco::GenParticleCollection>(edm::InputTag("prunedGenParticles"));
  packedgenParticlesToken = consumes<edm::View<pat::PackedGenParticle> > (edm::InputTag("packedGenParticles")); //ATbbf
  genInfoToken = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
  genJetsToken = consumes<edm::View<reco::GenJet> >(edm::InputTag("slimmedGenJets")); //ATjets
  // GENCandidatesToken = consumes<edm::View<pat::CompositeCandidate> >(edm::InputTag("GENLevel"));
  consumesMany<LHEEventProduct>();
  candToken = consumes<edm::View<pat::CompositeCandidate> >(edm::InputTag(theCandLabel));

  triggerResultToken = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults"));
  triggerResultPATToken = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","PAT"));
  vtxToken = consumes<vector<reco::Vertex> >(edm::InputTag("goodPrimaryVertices"));
  jetToken = consumes<edm::View<pat::Jet> >(edm::InputTag("cleanJets"));
  photonToken = consumes<pat::PhotonCollection>(edm::InputTag("slimmedPhotons"));
      
  metToken = consumes<pat::METCollection>(pset.getParameter<edm::InputTag>("metSrc"));
  theCovTag = consumes<math::Error<2>::type>(pset.getParameter<edm::InputTag>("covSrc"));
  metJESupToken = consumes<pat::METCollection>(pset.getParameter<edm::InputTag>("metJESup"));
  metJESdnToken = consumes<pat::METCollection>(pset.getParameter<edm::InputTag>("metJESdown"));
      
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
  if (theChannel!=SR) failedTreeLevel=noFailedTree;

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
    htxsToken = consumes<HTXS::HiggsClassification>(edm::InputTag("rivetProducerHTXS","HiggsClassification"));
    pileUpReweight = new PileUpWeight(myHelper.sampleType(), myHelper.setup());
  }

  Nevt_Gen = 0;
  Nevt_Gen_lumiBlock = 0;

  //For Efficiency studies
  gen_mumu = 0;gen_mumu_LeptonAcceptance = 0;gen_mumu_EtaAcceptance = 0;
  gen_etau = 0;gen_etau_LeptonAcceptance = 0;gen_etau_EtaAcceptance = 0;
  gen_mutau = 0;gen_mutau_LeptonAcceptance = 0;gen_mutau_EtaAcceptance = 0;
  gen_tautau = 0;gen_tautau_LeptonAcceptance = 0;gen_tautau_EtaAcceptance = 0;
  gen_emu = 0;gen_emu_LeptonAcceptance = 0;gen_emu_EtaAcceptance = 0;
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
  if (!skipTauDataMCWeight && isMC) {
    string period;
    if (year==2016 && preVFP) period="UL2016_preVFP";
    else if (year==2016 && !preVFP) period="UL2016_postVFP";
    else if (year==2017) period="UL2017";
    else period="UL2018";

    DeepTauSF_VSe_ETau = new TauIDSFTool(period,"DeepTau2017v2p1VSe","Tight");
    DeepTauSF_VSmu_ETau = new TauIDSFTool(period,"DeepTau2017v2p1VSmu","VLoose");
    DeepTauSF_VSjet_ETau = new TauIDSFTool(period,"DeepTau2017v2p1VSjet","Medium", "Tight", false, true);
    DeepTauSF_VSe_MuTau = new TauIDSFTool(period,"DeepTau2017v2p1VSe","VVLoose");
    DeepTauSF_VSmu_MuTau = new TauIDSFTool(period,"DeepTau2017v2p1VSmu","Tight");
    DeepTauSF_VSjet_MuTau = new TauIDSFTool(period,"DeepTau2017v2p1VSjet","Medium", "VVLoose", false, true);
    DeepTauSF_VSe_TauTau = new TauIDSFTool(period,"DeepTau2017v2p1VSe","VVLoose");
    DeepTauSF_VSmu_TauTau = new TauIDSFTool(period,"DeepTau2017v2p1VSmu","VLoose");
    DeepTauSF_VSjet_TauTau = new TauIDSFTool(period,"DeepTau2017v2p1VSjet","Medium", "VVLoose", false, true);

    if (do_MET_Recoil) {
      TString recoilFile;
      if (year == 2016) recoilFile="HTT-utilities/RecoilCorrections/data/TypeI-PFMet_Run2016_legacy.root";
      else if (year == 2017) recoilFile="HTT-utilities/RecoilCorrections/data/Type1_PFMET_2017.root";
      else recoilFile="HTT-utilities/RecoilCorrections/data/TypeI-PFMet_Run2018.root";

      recoilPFMetCorrector = new RecoilCorrector(recoilFile);

      if (year == 2016) recoilFile="HTT-utilities/RecoilCorrections/data/PFMEtSys_2016.root";
      else if (year == 2017) recoilFile="HTT-utilities/RecoilCorrections/data/PFMEtSys_2017.root";
      else recoilFile="HTT-utilities/RecoilCorrections/data/PFMEtSys_2018.root";

      recoilPFMetSyst = new MEtSys(recoilFile);
    }

    uncSources.push_back("Total");
    uncSources.push_back("Abs");
    uncSources.push_back("Abs_year");
    uncSources.push_back("BBEC1");
    uncSources.push_back("BBEC1_year");
    uncSources.push_back("EC2");
    uncSources.push_back("EC2_year");
    uncSources.push_back("FlavQCD");
    uncSources.push_back("HF");
    uncSources.push_back("HF_year");
    uncSources.push_back("RelBal");
    uncSources.push_back("RelSample_year");
  }
}

LLNtupleMaker::~LLNtupleMaker()
{
  delete pileUpReweight;
}


// ------------ method called for each event  ------------
void LLNtupleMaker::analyze(const edm::Event& event, const edm::EventSetup& eSetup)
{
  myTree->InitializeVariables();
  //cout<<"LLNtupleMaker:"<<theChannel<<endl;
  //----------------------------------------------------------------------
  // Analyze MC truth; collect MC weights and update counters (this is done for all generated events,
  // including those that do not pass skim, trigger etc!)

  bool InEtaAcceptance = false;   // All 4 gen leptons in eta acceptance
  bool InEtaPtAcceptance = false; // All 4 gen leptons in eta,pT acceptance

  const reco::Candidate * genH = 0;
  std::vector<const reco::Candidate *> genV;
  std::vector<const reco::Candidate *> genLeps;
  std::vector<const reco::Candidate *> genAssocLeps;
  std::vector<const reco::Candidate *> genVisLeps;
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
    // event.getByToken(GENCandidatesToken, GENCandidates);//ATMELA      

    edm::Handle<HTXS::HiggsClassification> htxs;
    event.getByToken(htxsToken,htxs);

    // cout<<"begin MC history"<<endl;
    MCHistoryTools mch(event, sampleName, genParticles, genInfo, genJets, packedgenParticles);
    genFinalState = mch.genFinalState();
    // cout<<"initiate"<<endl;
    // genProcessId = mch.getProcessID();
    genHEPMCweight = mch.gethepMCweight(); // Overridden by LHEHandler if genHEPMCweight==1.
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
    // cout<<"weights"<<endl;
    htxsNJets = htxs->jets30.size();
    htxsHPt = htxs->higgs.Pt();
    htxs_stage0_cat = htxs->stage0_cat;
    htxs_stage1p0_cat = htxs->stage1_cat_pTjet30GeV;
    htxs_stage1p1_cat = htxs->stage1_1_cat_pTjet30GeV;
    htxs_stage1p2_cat = htxs->stage1_2_cat_pTjet30GeV;
    htxs_errorCode=htxs->errorCode;
    htxs_prodMode= htxs->prodMode;
    // cout<<"htxs"<<endl;
    genExtInfo = mch.genAssociatedFS();

    //Information on generated candidates, will be used later
    genH = mch.genH();
    genV = mch.genVs();
    genLeps     = mch.sortedGenZZLeps();
    genAssocLeps = mch.genAssociatedLeps();
    genFSR       = mch.genFSR();
    genVisLeps	= mch.sortedVisGenZZLeps();
    genTauNus  	= mch.genTauNus();
        
    genIso       = mch.genIso(); //AT
    genJet       = mch.GenJets(); //ATjets
    genCleanedJet= mch.GenCleanedJets(); //ATjets

    // cout<<genZ.size()<<","<<genLeps.size()<<endl;



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
      
    // if (genFinalState!=BUGGY && genFinalState!=NONE) {
    if (genLeps.size()>=2) {
      if (genH!=0) {
        FillLLGenInfo(genLeps.at(0)->pdgId()*genLeps.at(1)->pdgId(), genH->p4());
      }
      else if (genV.size()>0) {
        FillLLGenInfo(genLeps.at(0)->pdgId()*genLeps.at(1)->pdgId(), genV.at(0)->p4());
      }
      else {
        FillLLGenInfo(genLeps.at(0)->pdgId()*genLeps.at(1)->pdgId(), genLeps.at(0)->p4()+genLeps.at(1)->p4());
      }
      FillLepGenInfo(genLeps.at(0)->pdgId(),genLeps.at(1)->pdgId(),genLeps.at(0)->p4(),genLeps.at(1)->p4());
      math::XYZTLorentzVector genVisLep1p4,genVisLep2p4;
      if (abs(genVisLeps.at(0)->pdgId())!=15 || genTauNus.at(0)==0) genVisLep1p4=genVisLeps.at(0)->p4();
      else genVisLep1p4=genVisLeps.at(0)->p4()-genTauNus.at(0)->p4();
      if (abs(genVisLeps.at(1)->pdgId())!=15 || genTauNus.at(1)==0) genVisLep2p4=genVisLeps.at(1)->p4();
      else genVisLep2p4=genVisLeps.at(1)->p4()-genTauNus.at(1)->p4();
      FillVisLLGenInfo(genVisLeps.at(0)->pdgId()*genVisLeps.at(1)->pdgId(), genVisLep1p4+genVisLep2p4);
      FillVisLepGenInfo(genVisLeps.at(0)->pdgId(),genVisLeps.at(1)->pdgId(), genVisLep1p4,genVisLep2p4);
      FillLepGenIso(genIso.at(0), genIso.at(1));
    }
    else {
      if (genH!=0) {
        FillLLGenInfo(0, genH->p4());
      }
      else if (genV.size()>0) {
        FillLLGenInfo(0, genV.at(0)->p4());
      }
    }
    // if (genLeps.size()!=2) cout<<"[WARNING] Number of leptons = "<<genLeps.size()<<"!"<<endl;
    // else {
    //   if (genH!=0) {
    //     FillLLGenInfo(genLeps.at(0)->pdgId()*genLeps.at(1)->pdgId(), genH->p4());
    //   }
    //   else {
    //     FillLLGenInfo(genLeps.at(0)->pdgId()*genLeps.at(1)->pdgId(), genLeps.at(0)->mother()->p4());
    //   }
    //   FillLepGenInfo(genLeps.at(0)->pdgId(),genLeps.at(1)->pdgId(),genLeps.at(0)->p4(),genLeps.at(1)->p4());
    //   math::XYZTLorentzVector genVisLep1p4,genVisLep2p4;
    //   if (abs(genVisLeps.at(0)->pdgId())!=15 || genTauNus.at(0)==0) genVisLep1p4=genVisLeps.at(0)->p4();
    //   else genVisLep1p4=genVisLeps.at(0)->p4()-genTauNus.at(0)->p4();
    //   if (abs(genVisLeps.at(1)->pdgId())!=15 || genTauNus.at(1)==0) genVisLep2p4=genVisLeps.at(1)->p4();
    //   else genVisLep2p4=genVisLeps.at(1)->p4()-genTauNus.at(1)->p4();
    //   FillVisLLGenInfo(genVisLeps.at(0)->pdgId()*genVisLeps.at(1)->pdgId(), genVisLep1p4+genVisLep2p4);
    //   FillVisLepGenInfo(genVisLeps.at(0)->pdgId(),genVisLeps.at(1)->pdgId(), genVisLep1p4,genVisLep2p4);
    //   FillLepGenIso(genIso.at(0), genIso.at(1));
    // }
    // }

      // keep track of sum of weights
      addweight(gen_sumPUWeight, PUWeight);
      addweight(gen_sumGenMCWeight, genHEPMCweight);
      addweight(gen_sumWeights, PUWeight*genHEPMCweight);

      mch.genAcceptance(InEtaAcceptance, InEtaPtAcceptance);

    addweight(Nevt_Gen_lumiBlock, 1); // Needs to be outside the if-block

    if (genFinalState == mumu) {
      addweight(gen_mumu, 1);
      if (InEtaAcceptance) addweight(gen_mumu_EtaAcceptance, 1);
      if (InEtaPtAcceptance) addweight(gen_mumu_LeptonAcceptance, 1);
    } else if (genFinalState == etau) {
      addweight(gen_etau, 1);
      if (InEtaAcceptance) addweight(gen_etau_EtaAcceptance, 1);
      if (InEtaPtAcceptance) addweight(gen_etau_LeptonAcceptance, 1);
    } else if (genFinalState == mutau) {
      addweight(gen_mutau, 1);
      if (InEtaAcceptance) addweight(gen_mutau_EtaAcceptance, 1);
      if (InEtaPtAcceptance) addweight(gen_mutau_LeptonAcceptance, 1);
    } else if (genFinalState == tautau) {
      addweight(gen_tautau, 1);
      if (InEtaAcceptance) addweight(gen_tautau_EtaAcceptance, 1);
      if (InEtaPtAcceptance) addweight(gen_tautau_LeptonAcceptance, 1);
    } else if (genFinalState == emu) {
      addweight(gen_emu, 1);
      if (InEtaAcceptance) addweight(gen_emu_EtaAcceptance, 1);
      if (InEtaPtAcceptance) addweight(gen_emu_LeptonAcceptance, 1);
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
    else edm::LogError("") << "LL collection not found in non-loose electron flow. This should never happen";
  }
  const edm::View<pat::CompositeCandidate>* cands = candHandle.product();

  // For Z+L CRs, we want only events with exactly 1 Z+l candidate. FIXME: this has to be reviewed.

  // Retrieve trigger results
  Handle<edm::TriggerResults> triggerResults, triggerResultsPAT;
  event.getByToken(triggerResultToken, triggerResults);
  event.getByToken(triggerResultPATToken, triggerResultsPAT);

  bool failed = false;

  // Apply MC filter (skip event)
  // Heshy note: I'm not turning return into failed = true because it looks like it's applied even if !skipEmptyEvents.
  //             It only does anything if the MCFILTER variable is set in the csv file, which is not currently the case.
  if (isMC && !(myHelper.passMCFilter(event,triggerResults))) return;

  if (!myHelper.passMETFilter(event,triggerResultsPAT,metWord)) return;

  // Apply skim
  bool evtPassSkim = myHelper.passSkim(event,triggerResults,trigWord);
  if (applySkim && !evtPassSkim) failed = true;       //but gen information will still be recorded if failedTreeLevel != 0

  // Apply trigger request (skip event)
  bool evtPassTrigger = myHelper.passTrigger(event,triggerResults,trigWord);
  if (applyTrigger && !evtPassTrigger) failed = true; //but gen information will still be recorded if failedTreeLevel != 0
	
  // Apply MET trigger request (skip event)
  // evtPassMETTrigger = myHelper.passMETTrigger(event,triggerResults);

  if (skipEmptyEvents && !failedTreeLevel && (cands->size() == 0 || failed)) return; // Skip events with no candidate, unless skipEmptyEvents = false or failedTreeLevel != 0

  //Fill MC truth information
  if (isMC) FillKFactors(genInfo, genV, genLeps);

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
    
  const pat::MET &met = metHandle->front();
  PFMET = met.pt();
  PFMETPhi = met.phi();
  METx = met.px();
  METy = met.py();
  covMET[0][0] = (*covHandle)(0,0);
  covMET[1][0] = (*covHandle)(1,0);
  covMET[0][1] = covMET[1][0]; // (1,0) is the only one saved
  covMET[1][1] = (*covHandle)(1,1);
    
  GenMET=GenMETPhi=-99;
  if (isMC && met.genMET()){
      GenMET = met.genMET()->pt();
      GenMETPhi = met.genMET()->phi();
      GenMETx = met.genMET()->px();
      GenMETy = met.genMET()->py();
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
  // vector<Int_t> CRFLAG(cands->size());
  // for( edm::View<pat::CompositeCandidate>::const_iterator cand = cands->begin(); cand != cands->end(); ++cand) {
  //   if (failed) break; //don't waste time on this
  //   size_t icand= cand-cands->begin();

  //   //    int candChannel = cand->userFloat("candChannel"); // This is currently the product of pdgId of leptons (eg 14641, 28561, 20449)

  //   if (theChannel==ZLL) {
  //     // Cross check region for Z + 1 loose electron + 1 loose TLE (defined only in loose_ele paths) 
  //     if (is_loose_ele_selection) {
  //       if (cand->userFloat("isBestCRZLL")&&cand->userFloat("CRZLL")) set_bit(CRFLAG[icand],ZLL);
  //     }

  //     // AA CRs
  //     if (cand->userFloat("isBestCRZLLss")&&cand->userFloat("CRZLLss")) set_bit(CRFLAG[icand],CRZLLss);

  //     // A CRs
  //     if (cand->userFloat("isBestCRZLLos_2P2F")&&cand->userFloat("CRZLLos_2P2F")) set_bit(CRFLAG[icand],CRZLLos_2P2F);
  //     if (cand->userFloat("isBestCRZLLos_3P1F")&&cand->userFloat("CRZLLos_3P1F")) set_bit(CRFLAG[icand],CRZLLos_3P1F);

  //     if (CRFLAG[icand]) { // This candidate belongs to one of the CRs: perform additional jet cleaning.
  //       // Note that this is (somewhat incorrectly) done per-event, so there could be some over-cleaning in events with >1 CR candidate.
  //       for (unsigned i=0; i<cleanedJets.size(); ++i) {
  //         if (cleanedJets[i]!=0  && (!jetCleaner::isGood(*cand, *(cleanedJets[i])))) {
  //           cleanedJets[i]=0;
  //         }
  //       }
  //     }
  //   }
  // }

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

    float abseta = fabs(cleanedJets[i]->eta());

    if(pt_nominal>25 && abseta<2.4) {
      if (cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt25BTagged_bTagSF;
    }
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

  //VBF jets
  Float_t best_DiJetMass=0;
  for (unsigned i=0; i<cleanedJets.size(); ++i) {
    for (unsigned j=i+1; j<cleanedJets.size(); ++j) {
      Float_t _DeltaEtaJJ, _DiJetMass;
      _DeltaEtaJJ=fabs(cleanedJets[i]->eta()-cleanedJets[j]->eta());
      _DiJetMass=(cleanedJets[i]->p4()+cleanedJets[j]->p4()).M();
      if (_DeltaEtaJJ<2.5 || _DiJetMass<350) continue;
      if (_DiJetMass>DiJetMass) {
        DiJetMass=_DiJetMass;
        DeltaEtaJJ=_DeltaEtaJJ;
        VBFJetIdx1=i;
        VBFJetIdx2=j;
      }
    }
  }

  //MET recoil
  if (do_MET_Recoil) {
    recoilPFMetCorrector->CorrectByMeanResolution(
      METx,
      METy,
      GenLLPt*cos(GenLLPhi),
      GenLLPt*sin(GenLLPhi),
      GenVisLLPt*cos(GenVisLLPhi),
      GenVisLLPt*sin(GenVisLLPhi),
      nCleanedJetsPt30,
      METxRecoil,
      METyRecoil
    );
  }
  else {
    METxRecoil=METx;
    METyRecoil=METy;
  }
  PFMETRecoil=TMath::Sqrt(METxRecoil*METxRecoil+METyRecoil*METyRecoil);
  PFMETPhiRecoil=TMath::ATan2(METyRecoil,METxRecoil);

  // Now we can write the variables for candidates
  int nFilled=0;
  for(unsigned int i = 0; i < cands->size(); ++i) {
    const pat::CompositeCandidate& c = (*cands)[i];
    pat::CompositeCandidate cand(c);
    if (failed) break; //don't waste time on this
    // size_t icand= cand-cands->begin();

    if (!(bool)(cand.userFloat("isBestCand")) ) continue; // Skip events other than the best cand

    if (abs(cand.daughter(0)->pdgId()*cand.daughter(1)->pdgId())==195 || abs(cand.daughter(0)->pdgId()*cand.daughter(1)->pdgId())==165) {
      TLorentzVector tau1,METRecoil;
      tau1.SetPxPyPzE(cand.daughter(0)->px(),cand.daughter(0)->py(),cand.daughter(0)->pz(),cand.daughter(0)->energy());
      METRecoil.SetPxPyPzE(METxRecoil,METyRecoil,0,std::hypot(METxRecoil,METyRecoil));
      MtLMET=(tau1+METRecoil).Mt();
    }
    else {
      MtLMET=-99;
    }

    cand.addUserFloat("NmedB",nCleanedJetsPt25BTagged_bTagSF);
    cand.addUserFloat("MtLMET",MtLMET);
    if (!cut(cand)) continue;

    //For the SR, also fold information about acceptance in CRflag.
    // if (isMC && (theChannel==ZZ)) {
    //   if (InEtaAcceptance)   set_bit(CRFLAG[icand],28);
    //   if (InEtaPtAcceptance) set_bit(CRFLAG[icand],29);
    // }
    FillCandidate(cand, evtPassTrigger&&evtPassSkim, event);

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


void LLNtupleMaker::FillJet(const pat::Jet& jet)
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

void LLNtupleMaker::FillPhoton(int year, const pat::Photon& photon)
{
   PhotonPt  .push_back( photon.pt());
   PhotonEta .push_back( photon.eta());
   PhotonPhi .push_back( photon.phi());
   
   PhotonIsCutBasedLooseID .push_back( PhotonIDHelper::isCutBasedID_Loose(year, photon) );
}

float LLNtupleMaker::EvalSpline(TSpline3* const& sp, float xval){
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

void LLNtupleMaker::FillKFactors(edm::Handle<GenEventInfoProduct>& genInfo, std::vector<const reco::Candidate *>& genZ, std::vector<const reco::Candidate *>& genLeps){
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
    Float_t GenHMass=0,GenHPt=0;
    if (genZ.size()>=2) {
      GenHMass=(genZ.at(0)->p4()+genZ.at(1)->p4()).M();
      GenHPt=(genZ.at(0)->p4()+genZ.at(1)->p4()).Pt();
    }
    if (apply_K_NNLOQCD_ZZGG>0 && apply_K_NNLOQCD_ZZGG!=3){
      if (spkfactor_ggzz_nnlo[0]!=0) KFactor_QCD_ggZZ_Nominal = LLNtupleMaker::EvalSpline(spkfactor_ggzz_nnlo[0], GenHMass);
      if (spkfactor_ggzz_nnlo[1]!=0) KFactor_QCD_ggZZ_PDFScaleDn = LLNtupleMaker::EvalSpline(spkfactor_ggzz_nnlo[1], GenHMass);
      if (spkfactor_ggzz_nnlo[2]!=0) KFactor_QCD_ggZZ_PDFScaleUp = LLNtupleMaker::EvalSpline(spkfactor_ggzz_nnlo[2], GenHMass);
      if (spkfactor_ggzz_nnlo[3]!=0) KFactor_QCD_ggZZ_QCDScaleDn = LLNtupleMaker::EvalSpline(spkfactor_ggzz_nnlo[3], GenHMass);
      if (spkfactor_ggzz_nnlo[4]!=0) KFactor_QCD_ggZZ_QCDScaleUp = LLNtupleMaker::EvalSpline(spkfactor_ggzz_nnlo[4], GenHMass);
      if (spkfactor_ggzz_nnlo[5]!=0) KFactor_QCD_ggZZ_AsDn = LLNtupleMaker::EvalSpline(spkfactor_ggzz_nnlo[5], GenHMass);
      if (spkfactor_ggzz_nnlo[6]!=0) KFactor_QCD_ggZZ_AsUp = LLNtupleMaker::EvalSpline(spkfactor_ggzz_nnlo[6], GenHMass);
      if (spkfactor_ggzz_nnlo[7]!=0) KFactor_QCD_ggZZ_PDFReplicaDn = LLNtupleMaker::EvalSpline(spkfactor_ggzz_nnlo[7], GenHMass);
      if (spkfactor_ggzz_nnlo[8]!=0) KFactor_QCD_ggZZ_PDFReplicaUp = LLNtupleMaker::EvalSpline(spkfactor_ggzz_nnlo[8], GenHMass);
      if (apply_K_NNLOQCD_ZZGG==2){
        if (spkfactor_ggzz_nlo[0]!=0){
          float divisor = LLNtupleMaker::EvalSpline(spkfactor_ggzz_nlo[0], GenHMass);
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
      if (spkfactor_ggzz_nlo[0]!=0) KFactor_QCD_ggZZ_Nominal = LLNtupleMaker::EvalSpline(spkfactor_ggzz_nlo[0], GenHMass);
      if (spkfactor_ggzz_nlo[1]!=0) KFactor_QCD_ggZZ_PDFScaleDn = LLNtupleMaker::EvalSpline(spkfactor_ggzz_nlo[1], GenHMass);
      if (spkfactor_ggzz_nlo[2]!=0) KFactor_QCD_ggZZ_PDFScaleUp = LLNtupleMaker::EvalSpline(spkfactor_ggzz_nlo[2], GenHMass);
      if (spkfactor_ggzz_nlo[3]!=0) KFactor_QCD_ggZZ_QCDScaleDn = LLNtupleMaker::EvalSpline(spkfactor_ggzz_nlo[3], GenHMass);
      if (spkfactor_ggzz_nlo[4]!=0) KFactor_QCD_ggZZ_QCDScaleUp = LLNtupleMaker::EvalSpline(spkfactor_ggzz_nlo[4], GenHMass);
      if (spkfactor_ggzz_nlo[5]!=0) KFactor_QCD_ggZZ_AsDn = LLNtupleMaker::EvalSpline(spkfactor_ggzz_nlo[5], GenHMass);
      if (spkfactor_ggzz_nlo[6]!=0) KFactor_QCD_ggZZ_AsUp = LLNtupleMaker::EvalSpline(spkfactor_ggzz_nlo[6], GenHMass);
      if (spkfactor_ggzz_nlo[7]!=0) KFactor_QCD_ggZZ_PDFReplicaDn = LLNtupleMaker::EvalSpline(spkfactor_ggzz_nlo[7], GenHMass);
      if (spkfactor_ggzz_nlo[8]!=0) KFactor_QCD_ggZZ_PDFReplicaUp = LLNtupleMaker::EvalSpline(spkfactor_ggzz_nlo[8], GenHMass);
    }

    if (genFinalState!=BUGGY){
      if (genLeps.size()==4) {
        // Calculate NNLO/NLO QCD K factors for qqZZ
        if (apply_K_NNLOQCD_ZZQQB){
          bool sameflavor=(genLeps.at(0)->pdgId()*genLeps.at(1)->pdgId() == genLeps.at(2)->pdgId()*genLeps.at(3)->pdgId());
          Float_t GenZ1Phi=0,GenZ2Phi=0;
          if (genZ.size()>=2) {
            GenZ1Phi=genZ.at(0)->phi();
            GenZ2Phi=genZ.at(1)->phi();
          }
          KFactor_QCD_qqZZ_dPhi = kfactor_qqZZ_qcd_dPhi(fabs(GenZ1Phi-GenZ2Phi), (sameflavor) ? 1 : 2);
          KFactor_QCD_qqZZ_M    = kfactor_qqZZ_qcd_M(GenHMass, (sameflavor) ? 1 : 2, 2) / kfactor_qqZZ_qcd_M(GenHMass, (sameflavor) ? 1 : 2, 1);
          KFactor_QCD_qqZZ_Pt   = kfactor_qqZZ_qcd_Pt(GenHPt, (sameflavor) ? 1 : 2);
        }
        // Calculate NLO EWK K factors for qqZZ
        if (apply_K_NLOEW_ZZQQB){
          TLorentzVector GENZ1Vec, GENZ2Vec, GENZZVec;
          if (genZ.size()>=2) {
            GENZ1Vec.SetPtEtaPhiM(genZ.at(0)->pt(),genZ.at(0)->eta(),genZ.at(0)->phi(),genZ.at(0)->mass());
            GENZ2Vec.SetPtEtaPhiM(genZ.at(1)->pt(),genZ.at(1)->eta(),genZ.at(1)->phi(),genZ.at(1)->mass());
            GENZZVec = GENZ1Vec + GENZ2Vec;
          }
          KFactor_EW_qqZZ = EwkCorrections::getEwkCorrections(genParticles, ewkTable, genInfoP, GENZ1Vec, GENZ2Vec);

          bool sameflavor=(genLeps.at(0)->pdgId()*genLeps.at(1)->pdgId() == genLeps.at(2)->pdgId()*genLeps.at(3)->pdgId());
          float K_NNLO_LO = kfactor_qqZZ_qcd_M(GenHMass, (sameflavor) ? 1 : 2, 2);
          float rho = genLeps.size()>=4 ? GENZZVec.Pt()/(genLeps.at(0)->pt()+genLeps.at(1)->pt()+genLeps.at(2)->pt()+genLeps.at(3)->pt()) : 0.;
          if (rho<0.3) KFactor_EW_qqZZ_unc = fabs((K_NNLO_LO-1.)*(1.-KFactor_EW_qqZZ));
          else KFactor_EW_qqZZ_unc = fabs(1.-KFactor_EW_qqZZ);
        }
      }
    }
  }

}


void LLNtupleMaker::FillCandidate(const pat::CompositeCandidate& cand, bool evtPass, const edm::Event& event)
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
  // LepBDT.clear();
  LepisCrack.clear();
  LepMissingHit.clear();
  LepConversionVeto.clear();
  LepChargedHadIso.clear();
  LepNeutralHadIso.clear();
  LepPhotonIso.clear();
  LepPUIsoComponent.clear();
  LepCombRelIsoPF.clear();

  LepSF.clear();
  LepSF_UncUp.clear();
	LepSF_UncDn.clear();
  LepSF_UncUp_RECO_syst.clear();
  LepSF_UncUp_RECO_stat.clear();
  LepSF_UncUp_ID_syst.clear();
  LepSF_UncUp_ID_stat.clear();
  LepSF_UncUp_ISO_syst.clear();
  LepSF_UncUp_ISO_stat.clear();
  LepSF_UncDn_RECO_syst.clear();
  LepSF_UncDn_RECO_stat.clear();
  LepSF_UncDn_ID_syst.clear();
  LepSF_UncDn_ID_stat.clear();
  LepSF_UncDn_ISO_syst.clear();
  LepSF_UncDn_ISO_stat.clear();
  LepSF_UncUp_uncert0.clear();
  LepSF_UncUp_uncert1.clear();
  LepSF_UncUp_syst_alleras.clear();
  LepSF_UncUp_syst_year.clear();
  LepSF_UncUp_syst_dm_year.clear();
  LepSF_UncUp_fakeEle.clear();
  LepSF_UncUp_fakeMu.clear();
  LepSF_UncDn_uncert0.clear();
  LepSF_UncDn_uncert1.clear();
  LepSF_UncDn_syst_alleras.clear();
  LepSF_UncDn_syst_year.clear();
  LepSF_UncDn_syst_dm_year.clear();
  LepSF_UncDn_fakeEle.clear();
  LepSF_UncDn_fakeMu.clear();


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

  // HLTMatch1.clear();
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
 
  fsrPt.clear();
  fsrEta.clear();
  fsrPhi.clear();
  fsrLept.clear();
  fsrLeptID.clear();
  fsrDR.clear();
  fsrGenPt.clear();

  LLSVPhi_up.clear();
  LLSVEta_up.clear();
  LLSVPt_up.clear();
  LLSVMass_up.clear();
  LLMass_up.clear();
  LLGoodMass_up.clear();
  LLSVPhi_dn.clear();
  LLSVEta_dn.clear();
  LLSVPt_dn.clear();
  LLSVMass_dn.clear();
  LLMass_dn.clear();
  LLGoodMass_dn.clear();
  // Names_shift.clear();

  vector<const reco::Candidate*> leptons;
  vector<const reco::Candidate*> fsrPhot;
  vector<short> fsrIndex;
  vector<string> labels;

  userdatahelpers::getSortedLeptons(cand, leptons, labels, fsrPhot, fsrIndex);

  LLMass = cand.mass();
  LLPt = cand.pt();
  LLEta = cand.eta();
  LLPhi = cand.phi();
  LLDR = deltaR(cand.daughter(0)->p4(),cand.daughter(1)->p4());
  TLorentzVector LLP4;
  LLP4.SetPtEtaPhiM(LLPt,LLEta,LLPhi,LLMass);

  LLFlav = getPdgId(cand.daughter(0)) * getPdgId(cand.daughter(1));

  LLisGoodTau = cand.userFloat("isGoodTau");

  pass_SingleTrigger = cand.userInt("pass_SingleTrigger");
  pass_CrossTrigger = cand.userInt("pass_CrossTrigger");
  pass_Trigger = cand.userInt("pass_Trigger");


  //-------------------------------------------------------------------------------------------------------
  //----------------------------------------SV FIT---------------------------------------------------------
  //-------------------------------------------------------------------------------------------------------
  // cout<<"Begin begin SVFit: "<<LLFlav<<endl;
  bool doSVFit=false;
  if (abs(LLFlav)==165 || abs(LLFlav)==195 || abs(LLFlav)==225 || abs(LLFlav)==143) doSVFit=true;

  bool swi;
  if (abs(leptons[0]->pdgId()) > abs(leptons[1]->pdgId())) swi=true;
  else if (abs(leptons[0]->pdgId()) < abs(leptons[1]->pdgId())) swi=false;
  else if (leptons[0]->pt() < leptons[1]->pt()) swi=true;
  else swi=false;
  int idx1,idx2;
  if (swi) {idx1=1;idx2=0;}
  else {idx1=0;idx2=1;}

  TVector3 tau1_3,tau2_3,met_3;
  tau1_3.SetXYZ(leptons[idx1]->px(),leptons[idx1]->py(),leptons[idx1]->pz());
  tau2_3.SetXYZ(leptons[idx2]->px(),leptons[idx2]->py(),leptons[idx2]->pz());
  met_3.SetXYZ(METx,METy,0);
  TVector3 zeta(leptons[idx1]->px()+leptons[idx2]->px(),leptons[idx1]->py()+leptons[idx2]->py(),0);
  if (zeta.Mag()==0) {
    Pzeta1=-999;
    Pzeta2=-999;
    }
  else {
    Float_t Pall1=met_3.Dot(zeta)/zeta.Mag();
    Float_t Pall2=(tau1_3+tau2_3+met_3).Dot(zeta)/zeta.Mag();
    Float_t Pvis=(tau1_3+tau2_3).Dot(zeta)/zeta.Mag();
    Pzeta1=Pall1-0.85*Pvis;
    Pzeta2=Pall2-0.85*Pvis;
  }

  //central
  TLorentzVector LLp4;
  LLp4.SetPtEtaPhiM(LLPt,LLEta,LLPhi,LLMass);
  TLorentzVector tau1,tau2,MET,METTau,METJet,METRaw,METRecoil;
  int pairType;
  int dm1,dm2;
  tau1.SetPxPyPzE(leptons[idx1]->px(),leptons[idx1]->py(),leptons[idx1]->pz(),leptons[idx1]->energy());
  tau2.SetPxPyPzE(leptons[idx2]->px(),leptons[idx2]->py(),leptons[idx2]->pz(),leptons[idx2]->energy());
  // MET.SetPxPyPzE(METx,METy,0,std::hypot(METx,METy));
  // cout<<"Fill other METs"<<endl;
  // METTau.SetPxPyPzE(METxTau,METyTau,0,std::hypot(METxTau,METyTau));
  // METJet.SetPxPyPzE(METxJet,METyJet,0,std::hypot(METxJet,METyJet));
  // METRaw.SetPxPyPzE(METxRaw,METyRaw,0,std::hypot(METxRaw,METyRaw));
  METRecoil.SetPxPyPzE(METxRecoil,METyRecoil,0,std::hypot(METxRecoil,METyRecoil));

  if (abs(LLFlav)==195) pairType=0;
  else if (abs(LLFlav)==165) pairType=1;
  else if (abs(LLFlav)==225) pairType=2;
  else pairType=3;
  dm1=abs(leptons[idx1]->pdgId())==15?userdatahelpers::getUserFloat(leptons[idx1],"decayMode"):-1;
  dm2=abs(leptons[idx2]->pdgId())==15?userdatahelpers::getUserFloat(leptons[idx2],"decayMode"):-1;
  if (doSVFit) {
    // cout<<"Begin SVFit"<<endl;
    SVfit algo_central_Recoil(0,tau1,tau2,METRecoil,covMET,pairType,dm1,dm2);
    std::vector<double> results_central_Recoil=algo_central_Recoil.FitAndGetResult();
    LLSVPt=results_central_Recoil.at(0);
    LLSVEta=results_central_Recoil.at(1);
    LLSVPhi=results_central_Recoil.at(2);
    LLSVMass=results_central_Recoil.at(3);
  }
  else {
    LLSVPt=-99;
    LLSVEta=-99;
    LLSVPhi=-99;
    LLSVMass=-99;
  }

  //-------------------------------------------------------------------------------
  //-----------SYSTEMATICS AFFECTING SVFIT MASS AND VISIBLE MASS-------------------
  //-------------------------------------------------------------------------------

  // Electron energy corrections
  std::vector<std::string> correctionNames={"scale_stat","scale_syst","scale_gain","sigma_rho","sigma_phi"};
  std::vector<std::string> NPNames={"CMS_scale_e_stat_year","CMS_scale_e_syst","CMS_scale_e_gain_year","CMS_res_e_rho","CMS_res_e_rho"};
  if (theChannel==SR) {
    if (doSVFit && (abs(LLFlav)==165 || abs(LLFlav)==143)) {
      TLorentzVector tau1_up, tau1_dn;
      for (size_t iNP=0;iNP<correctionNames.size();iNP++) {
        tau1_up=tau1*userdatahelpers::getUserFloat(leptons[idx1],(correctionNames[iNP]+"_up").c_str());
        SVfit algo_up(0,tau1_up,tau2,METRecoil,covMET,pairType,dm1,dm2);
        std::vector<double> results_up=algo_up.FitAndGetResult();
        LLSVPt_up.push_back(results_up.at(0));
        LLSVEta_up.push_back(results_up.at(1));
        LLSVPhi_up.push_back(results_up.at(2));
        LLSVMass_up.push_back(results_up.at(3));
        LLMass_up.push_back((LLP4-tau1+tau1_up).M());
        tau1_dn=tau1*userdatahelpers::getUserFloat(leptons[idx1],(correctionNames[iNP]+"_dn").c_str());
        SVfit algo_dn(0,tau1_dn,tau2,METRecoil,covMET,pairType,dm1,dm2);
        std::vector<double> results_dn=algo_dn.FitAndGetResult();
        LLSVPt_dn.push_back(results_dn.at(0));
        LLSVEta_dn.push_back(results_dn.at(1));
        LLSVPhi_dn.push_back(results_dn.at(2));
        LLSVMass_dn.push_back(results_dn.at(3));
        LLMass_dn.push_back((LLP4-tau1+tau1_dn).M());
        // Names_shift.push_back(NPNames[iNP]);
      }
    }
    else {
      for (size_t iNP=0;iNP<correctionNames.size();iNP++) {
        LLSVPt_up.push_back(LLSVPt);
        LLSVEta_up.push_back(LLSVEta);
        LLSVPhi_up.push_back(LLSVPhi);
        LLSVMass_up.push_back(LLSVMass);
        LLMass_up.push_back(LLMass);
        LLSVPt_dn.push_back(LLSVPt);
        LLSVEta_dn.push_back(LLSVEta);
        LLSVPhi_dn.push_back(LLSVPhi);
        LLSVMass_dn.push_back(LLSVMass);
        LLMass_dn.push_back(LLMass);
        // Names_shift.push_back(NPNames[iNP]);
      }
    }
  }

  // Muon energy corrections
  correctionNames.clear();
  correctionNames={"scale_total","sigma_total"};
  NPNames.clear();
  NPNames={"CMS_scale_m","CMS_res_m"};
  if (theChannel==SR) {
    TLorentzVector tau1_up, tau2_up, tau1_dn, tau2_dn;
    for (size_t iNP=0;iNP<correctionNames.size();iNP++) {
      if (abs(leptons[idx1]->pdgId())==13) {
        tau1_up=tau1*userdatahelpers::getUserFloat(leptons[idx1],(correctionNames[iNP]+"_up").c_str());
        tau1_dn=tau1*userdatahelpers::getUserFloat(leptons[idx1],(correctionNames[iNP]+"_dn").c_str());
      }
      else {
        tau1_up=tau1;
        tau1_dn=tau1;
      }
      if (abs(leptons[idx2]->pdgId())==13) {
        tau2_up=tau2*userdatahelpers::getUserFloat(leptons[idx2],(correctionNames[iNP]+"_up").c_str());
        tau2_dn=tau2*userdatahelpers::getUserFloat(leptons[idx2],(correctionNames[iNP]+"_dn").c_str());
      }
      else {
        tau2_up=tau2;
        tau2_dn=tau2;
      }
      LLMass_up.push_back((LLP4-tau1-tau2+tau1_up+tau2_up).M());
      LLMass_dn.push_back((LLP4-tau1-tau2+tau1_dn+tau2_dn).M());
      if (doSVFit && (abs(leptons[idx1]->pdgId())==13 || abs(leptons[idx2]->pdgId())==13)) {
        SVfit algo_up(0,tau1_up,tau2_up,METRecoil,covMET,pairType,dm1,dm2);
        std::vector<double> results_up=algo_up.FitAndGetResult();
        LLSVPt_up.push_back(results_up.at(0));
        LLSVEta_up.push_back(results_up.at(1));
        LLSVPhi_up.push_back(results_up.at(2));
        LLSVMass_up.push_back(results_up.at(3));
        SVfit algo_dn(0,tau1_dn,tau2_dn,METRecoil,covMET,pairType,dm1,dm2);
        std::vector<double> results_dn=algo_dn.FitAndGetResult();
        LLSVPt_dn.push_back(results_dn.at(0));
        LLSVEta_dn.push_back(results_dn.at(1));
        LLSVPhi_dn.push_back(results_dn.at(2));
        LLSVMass_dn.push_back(results_dn.at(3));
      }
      else {
        LLSVPt_up.push_back(LLSVPt);
        LLSVEta_up.push_back(LLSVEta);
        LLSVPhi_up.push_back(LLSVPhi);
        LLSVMass_up.push_back(LLSVMass);
        LLSVPt_dn.push_back(LLSVPt);
        LLSVEta_dn.push_back(LLSVEta);
        LLSVPhi_dn.push_back(LLSVPhi);
        LLSVMass_dn.push_back(LLSVMass);
      }
    }
  }

  // Tau lepton energy corrections
  correctionNames.clear();
  correctionNames={"Tau","Ele","Mu"};
  std::vector<std::string> correctionNames1={"isTESShifted","isEESShifted","isMESShifted"};
  NPNames.clear();
  NPNames={"CMS_scale_t_year","CMS_scale_efaket_year","CMS_scale_mfaket_year"};
  if (theChannel==SR) {
    TLorentzVector tau1_up, tau2_up, tau1_dn, tau2_dn;
    for (size_t iNP=0;iNP<correctionNames.size();iNP++) {
      bool changed=false;
      if (iNP<correctionNames.size()-1) {
        if (abs(leptons[idx1]->pdgId())==15 && userdatahelpers::getUserInt(leptons[idx1],correctionNames1[iNP].c_str())) {
          changed=true;
          tau1_up.SetPxPyPzE(userdatahelpers::getUserFloat(leptons[idx1],("px_"+correctionNames[iNP]+"Up").c_str()),userdatahelpers::getUserFloat(leptons[idx1],("py_"+correctionNames[iNP]+"Up").c_str()),userdatahelpers::getUserFloat(leptons[idx1],("pz_"+correctionNames[iNP]+"Up").c_str()),userdatahelpers::getUserFloat(leptons[idx1],("e_"+correctionNames[iNP]+"Up").c_str()));
          tau1_dn.SetPxPyPzE(userdatahelpers::getUserFloat(leptons[idx1],("px_"+correctionNames[iNP]+"Down").c_str()),userdatahelpers::getUserFloat(leptons[idx1],("py_"+correctionNames[iNP]+"Down").c_str()),userdatahelpers::getUserFloat(leptons[idx1],("pz_"+correctionNames[iNP]+"Down").c_str()),userdatahelpers::getUserFloat(leptons[idx1],("e_"+correctionNames[iNP]+"Down").c_str()));
        }
        else {
          tau1_up=tau1;
          tau1_dn=tau1;
        }
        if (abs(leptons[idx2]->pdgId())==15 && userdatahelpers::getUserInt(leptons[idx2],correctionNames1[iNP].c_str())) {
          changed=true;
          tau2_up.SetPxPyPzE(userdatahelpers::getUserFloat(leptons[idx2],("px_"+correctionNames[iNP]+"Up").c_str()),userdatahelpers::getUserFloat(leptons[idx2],("py_"+correctionNames[iNP]+"Up").c_str()),userdatahelpers::getUserFloat(leptons[idx2],("pz_"+correctionNames[iNP]+"Up").c_str()),userdatahelpers::getUserFloat(leptons[idx2],("e_"+correctionNames[iNP]+"Up").c_str()));
          tau2_dn.SetPxPyPzE(userdatahelpers::getUserFloat(leptons[idx2],("px_"+correctionNames[iNP]+"Down").c_str()),userdatahelpers::getUserFloat(leptons[idx2],("py_"+correctionNames[iNP]+"Down").c_str()),userdatahelpers::getUserFloat(leptons[idx2],("pz_"+correctionNames[iNP]+"Down").c_str()),userdatahelpers::getUserFloat(leptons[idx2],("e_"+correctionNames[iNP]+"Down").c_str()));
        }
        else {
          tau2_up=tau2;
          tau2_dn=tau2;
        }
      }
      else {
        if (abs(leptons[idx1]->pdgId())==15 && (userdatahelpers::getUserFloat(leptons[idx1],"genmatch")==2 || userdatahelpers::getUserFloat(leptons[idx1],"genmatch")==4)) {
          changed=true;
          tau1_up=1.01*tau1;
          tau1_dn=0.99*tau1;
        }
        else {
          tau1_up=tau1;
          tau1_dn=tau1;
        }
        if (abs(leptons[idx2]->pdgId())==15 && (userdatahelpers::getUserFloat(leptons[idx2],"genmatch")==2 || userdatahelpers::getUserFloat(leptons[idx2],"genmatch")==4)) {
          changed=true;
          tau2_up=1.01*tau2;
          tau2_dn=0.99*tau2;
        }
        else {
          tau2_up=tau2;
          tau2_dn=tau2;
        }
      }
      LLMass_up.push_back((LLP4-tau1-tau2+tau1_up+tau2_up).M());
      LLMass_dn.push_back((LLP4-tau1-tau2+tau1_dn+tau2_dn).M());
      if (doSVFit && changed) {
        SVfit algo_up(0,tau1_up,tau2_up,METRecoil,covMET,pairType,dm1,dm2);
        std::vector<double> results_up=algo_up.FitAndGetResult();
        LLSVPt_up.push_back(results_up.at(0));
        LLSVEta_up.push_back(results_up.at(1));
        LLSVPhi_up.push_back(results_up.at(2));
        LLSVMass_up.push_back(results_up.at(3));
        SVfit algo_dn(0,tau1_dn,tau2_dn,METRecoil,covMET,pairType,dm1,dm2);
        std::vector<double> results_dn=algo_dn.FitAndGetResult();
        LLSVPt_dn.push_back(results_dn.at(0));
        LLSVEta_dn.push_back(results_dn.at(1));
        LLSVPhi_dn.push_back(results_dn.at(2));
        LLSVMass_dn.push_back(results_dn.at(3));
      }
      else {
        LLSVPt_up.push_back(LLSVPt);
        LLSVEta_up.push_back(LLSVEta);
        LLSVPhi_up.push_back(LLSVPhi);
        LLSVMass_up.push_back(LLSVMass);
        LLSVPt_dn.push_back(LLSVPt);
        LLSVEta_dn.push_back(LLSVEta);
        LLSVPhi_dn.push_back(LLSVPhi);
        LLSVMass_dn.push_back(LLSVMass);
      }
    }
  }

  // Jet energy corrections
  correctionNames.clear();
  correctionNames=uncSources;
  NPNames.clear();
  NPNames=uncSources;
  
  Handle<pat::METCollection> metJESupHandle;
  event.getByToken(metJESupToken, metJESupHandle);
  Handle<pat::METCollection> metJESdnHandle;
  event.getByToken(metJESdnToken, metJESdnHandle); 
  
  if (theChannel==SR) {
    TLorentzVector MET_up, MET_dn;
    for (size_t iNP=0;iNP<correctionNames.size();iNP++) {
      const pat::MET &met_up = metJESupHandle->at(iNP);
      Float_t metx_up, mety_up;
      if (do_MET_Recoil) {
        recoilPFMetCorrector->CorrectByMeanResolution(met_up.px(),met_up.py(),GenLLPt*cos(GenLLPhi),GenLLPt*sin(GenLLPhi),GenVisLLPt*cos(GenVisLLPhi),GenVisLLPt*sin(GenVisLLPhi),nCleanedJetsPt30,metx_up,mety_up);
      }
      else {
        metx_up=met_up.px();
        mety_up=met_up.py();
      }
      MET_up.SetPxPyPzE(metx_up,mety_up,0,std::hypot(metx_up,mety_up));

      const pat::MET &met_dn = metJESdnHandle->at(iNP);
      Float_t metx_dn, mety_dn;
      if (do_MET_Recoil) {
        recoilPFMetCorrector->CorrectByMeanResolution(met_dn.px(),met_dn.py(),GenLLPt*cos(GenLLPhi),GenLLPt*sin(GenLLPhi),GenVisLLPt*cos(GenVisLLPhi),GenVisLLPt*sin(GenVisLLPhi),nCleanedJetsPt30,metx_dn,mety_dn);
      }
      else {
        metx_dn=met_dn.px();
        mety_dn=met_dn.py();
      }
      MET_dn.SetPxPyPzE(metx_dn,mety_dn,0,std::hypot(metx_dn,mety_dn));

      LLMass_up.push_back(LLMass);
      LLMass_dn.push_back(LLMass);

      if (doSVFit) {
        SVfit algo_up(0,tau1,tau2,MET_up,covMET,pairType,dm1,dm2);
        std::vector<double> results_up=algo_up.FitAndGetResult();
        LLSVPt_up.push_back(results_up.at(0));
        LLSVEta_up.push_back(results_up.at(1));
        LLSVPhi_up.push_back(results_up.at(2));
        LLSVMass_up.push_back(results_up.at(3));
        SVfit algo_dn(0,tau1,tau2,MET_dn,covMET,pairType,dm1,dm2);
        std::vector<double> results_dn=algo_dn.FitAndGetResult();
        LLSVPt_dn.push_back(results_dn.at(0));
        LLSVEta_dn.push_back(results_dn.at(1));
        LLSVPhi_dn.push_back(results_dn.at(2));
        LLSVMass_dn.push_back(results_dn.at(3));
      }
      else {
        LLSVPt_up.push_back(-99);
        LLSVEta_up.push_back(-99);
        LLSVPhi_up.push_back(-99);
        LLSVMass_up.push_back(-99);
        LLSVPt_dn.push_back(-99);
        LLSVEta_dn.push_back(-99);
        LLSVPhi_dn.push_back(-99);
        LLSVMass_dn.push_back(-99);
      }
    }
  }

  // MET recoil
  correctionNames.clear();
  correctionNames={"CMS_scale_met","CMS_res_met"};
  NPNames.clear();
  NPNames={"CMS_scale_met","CMS_res_met"};
  if (theChannel==SR) {
    if (do_MET_Recoil) {
      TLorentzVector MET_up, MET_dn;

      Float_t metx_up, mety_up;
      Float_t metx_dn, mety_dn;

      LLMass_up.push_back(LLMass);
      LLMass_dn.push_back(LLMass);

      Float_t tmpGenLLPt=GenLLPt>99?99:GenLLPt;

      if (doSVFit) {

        recoilPFMetSyst->ApplyMEtSys(METxRecoil,METyRecoil,tmpGenLLPt*cos(GenLLPhi),tmpGenLLPt*sin(GenLLPhi),GenVisLLPt*cos(GenVisLLPhi),GenVisLLPt*sin(GenVisLLPhi),nCleanedJetsPt30,0,0,0,metx_up, mety_up);
        MET_up.SetPxPyPzE(metx_up,mety_up,0,std::hypot(metx_up,mety_up));

        recoilPFMetSyst->ApplyMEtSys(METxRecoil,METyRecoil,tmpGenLLPt*cos(GenLLPhi),tmpGenLLPt*sin(GenLLPhi),GenVisLLPt*cos(GenVisLLPhi),GenVisLLPt*sin(GenVisLLPhi),nCleanedJetsPt30,0,0,1,metx_up, mety_up);
        MET_dn.SetPxPyPzE(metx_dn,mety_dn,0,std::hypot(metx_dn,mety_dn));
        
        SVfit algo_up(0,tau1,tau2,MET_up,covMET,pairType,dm1,dm2);
        std::vector<double> results_up=algo_up.FitAndGetResult();
        LLSVPt_up.push_back(results_up.at(0));
        LLSVEta_up.push_back(results_up.at(1));
        LLSVPhi_up.push_back(results_up.at(2));
        LLSVMass_up.push_back(results_up.at(3));
        SVfit algo_dn(0,tau1,tau2,MET_dn,covMET,pairType,dm1,dm2);
        std::vector<double> results_dn=algo_dn.FitAndGetResult();
        LLSVPt_dn.push_back(results_dn.at(0));
        LLSVEta_dn.push_back(results_dn.at(1));
        LLSVPhi_dn.push_back(results_dn.at(2));
        LLSVMass_dn.push_back(results_dn.at(3));
      }
      else {
        LLSVPt_up.push_back(-99);
        LLSVEta_up.push_back(-99);
        LLSVPhi_up.push_back(-99);
        LLSVMass_up.push_back(-99);
        LLSVPt_dn.push_back(-99);
        LLSVEta_dn.push_back(-99);
        LLSVPhi_dn.push_back(-99);
        LLSVMass_dn.push_back(-99);
      }

      LLMass_up.push_back(LLMass);
      LLMass_dn.push_back(LLMass);

      if (doSVFit) {
        recoilPFMetSyst->ApplyMEtSys(METxRecoil,METyRecoil,tmpGenLLPt*cos(GenLLPhi),tmpGenLLPt*sin(GenLLPhi),GenVisLLPt*cos(GenVisLLPhi),GenVisLLPt*sin(GenVisLLPhi),nCleanedJetsPt30,0,1,0,metx_up, mety_up);
        MET_up.SetPxPyPzE(metx_up,mety_up,0,std::hypot(metx_up,mety_up));

        recoilPFMetSyst->ApplyMEtSys(METxRecoil,METyRecoil,tmpGenLLPt*cos(GenLLPhi),tmpGenLLPt*sin(GenLLPhi),GenVisLLPt*cos(GenVisLLPhi),GenVisLLPt*sin(GenVisLLPhi),nCleanedJetsPt30,0,1,1,metx_up, mety_up);
        MET_dn.SetPxPyPzE(metx_dn,mety_dn,0,std::hypot(metx_dn,mety_dn));

        SVfit algo_up(0,tau1,tau2,MET_up,covMET,pairType,dm1,dm2);
        std::vector<double> results_up=algo_up.FitAndGetResult();
        LLSVPt_up.push_back(results_up.at(0));
        LLSVEta_up.push_back(results_up.at(1));
        LLSVPhi_up.push_back(results_up.at(2));
        LLSVMass_up.push_back(results_up.at(3));
        SVfit algo_dn(0,tau1,tau2,MET_dn,covMET,pairType,dm1,dm2);
        std::vector<double> results_dn=algo_dn.FitAndGetResult();
        LLSVPt_dn.push_back(results_dn.at(0));
        LLSVEta_dn.push_back(results_dn.at(1));
        LLSVPhi_dn.push_back(results_dn.at(2));
        LLSVMass_dn.push_back(results_dn.at(3));
      }
      else {
        LLSVPt_up.push_back(-99);
        LLSVEta_up.push_back(-99);
        LLSVPhi_up.push_back(-99);
        LLSVMass_up.push_back(-99);
        LLSVPt_dn.push_back(-99);
        LLSVEta_dn.push_back(-99);
        LLSVPhi_dn.push_back(-99);
        LLSVMass_dn.push_back(-99);
      }
    }
    else {
      LLSVPt_up.push_back(LLSVPt);
      LLSVEta_up.push_back(LLSVEta);
      LLSVPhi_up.push_back(LLSVPhi);
      LLSVMass_up.push_back(LLSVMass);
      LLMass_up.push_back(LLMass);
      LLSVPt_dn.push_back(LLSVPt);
      LLSVEta_dn.push_back(LLSVEta);
      LLSVPhi_dn.push_back(LLSVPhi);
      LLSVMass_dn.push_back(LLSVMass);
      LLMass_dn.push_back(LLMass);
      LLSVPt_up.push_back(LLSVPt);
      LLSVEta_up.push_back(LLSVEta);
      LLSVPhi_up.push_back(LLSVPhi);
      LLSVMass_up.push_back(LLSVMass);
      LLMass_up.push_back(LLMass);
      LLSVPt_dn.push_back(LLSVPt);
      LLSVEta_dn.push_back(LLSVEta);
      LLSVPhi_dn.push_back(LLSVPhi);
      LLSVMass_dn.push_back(LLSVMass);
      LLMass_dn.push_back(LLMass);
    }
  }

  
  // //TES UP/DOWN
  // TLorentzVector tau1_tesup,tau2_tesup,met_tesup;
  // if (abs(leptons[idx1]->pdgId())==15 && userdatahelpers::getUserInt(leptons[idx1],"isTESShifted")) tau1_tesup.SetPxPyPzE(
  //     userdatahelpers::getUserFloat(leptons[idx1],"px_TauUp"),userdatahelpers::getUserFloat(leptons[idx1],"py_TauUp"),
  //     userdatahelpers::getUserFloat(leptons[idx1],"pz_TauUp"),userdatahelpers::getUserFloat(leptons[idx1],"e_TauUp"));
  // else tau1_tesup.SetPxPyPzE(leptons[idx1]->px(),leptons[idx1]->py(),leptons[idx1]->pz(),leptons[idx1]->energy());
  // if (abs(leptons[idx2]->pdgId())==15 && userdatahelpers::getUserInt(leptons[idx2],"isTESShifted")) tau2_tesup.SetPxPyPzE(
  //     userdatahelpers::getUserFloat(leptons[idx2],"px_TauUp"),userdatahelpers::getUserFloat(leptons[idx2],"py_TauUp"),
  //     userdatahelpers::getUserFloat(leptons[idx2],"pz_TauUp"),userdatahelpers::getUserFloat(leptons[idx2],"e_TauUp"));
  // else tau2_tesup.SetPxPyPzE(leptons[idx2]->px(),leptons[idx2]->py(),leptons[idx2]->pz(),leptons[idx2]->energy());
  // Z2Mass_TESup=(Z2p4-tau1+tau1_tesup-tau2+tau2_tesup).M();
  // met_tesup.SetPxPyPzE(METxUPTES,METyUPTES,0,std::hypot(METxUPTES,METyUPTES));
  // if (doSVFit) {
  //     if (METxUPTES==METx && METyUPTES==METy) {
  //         Z2SVPt_TESup=Z2SVPt;
  //         Z2SVEta_TESup=Z2SVEta;
  //         Z2SVPhi_TESup=Z2SVPhi;
  //         Z2SVMass_TESup=Z2SVMass;
  //     }
  //     else {
  //         SVfit algo_tesup(0,tau1_tesup,tau2_tesup,met_tesup,covMET,pairType,dm1,dm2);
  //         std::vector<double> results_tesup=algo_tesup.FitAndGetResult();
  //         Z2SVPt_TESup=results_tesup.at(0);
  //         Z2SVEta_TESup=results_tesup.at(1);
  //         Z2SVPhi_TESup=results_tesup.at(2);
  //         Z2SVMass_TESup=results_tesup.at(3);
  //     }
  // }

  // TLorentzVector tau1_tesdn,tau2_tesdn,met_tesdn;
  // if (abs(leptons[idx1]->pdgId())==15 && userdatahelpers::getUserInt(leptons[idx1],"isTESShifted")) tau1_tesdn.SetPxPyPzE(
  //     userdatahelpers::getUserFloat(leptons[idx1],"px_TauDown"),userdatahelpers::getUserFloat(leptons[idx1],"py_TauDown"),
  //     userdatahelpers::getUserFloat(leptons[idx1],"pz_TauDown"),userdatahelpers::getUserFloat(leptons[idx1],"e_TauDown"));
  // else tau1_tesdn.SetPxPyPzE(leptons[idx1]->px(),leptons[idx1]->py(),leptons[idx1]->pz(),leptons[idx1]->energy());
  // if (abs(leptons[idx2]->pdgId())==15 && userdatahelpers::getUserInt(leptons[idx2],"isTESShifted")) tau2_tesdn.SetPxPyPzE(
  //     userdatahelpers::getUserFloat(leptons[idx2],"px_TauDown"),userdatahelpers::getUserFloat(leptons[idx2],"py_TauDown"),
  //     userdatahelpers::getUserFloat(leptons[idx2],"pz_TauDown"),userdatahelpers::getUserFloat(leptons[idx2],"e_TauDown"));
  // else tau2_tesdn.SetPxPyPzE(leptons[idx2]->px(),leptons[idx2]->py(),leptons[idx2]->pz(),leptons[idx2]->energy());
  // Z2Mass_TESdn=(Z2p4-tau1+tau1_tesdn-tau2+tau2_tesdn).M();
  // met_tesdn.SetPxPyPzE(METxDOWNTES,METyDOWNTES,0,std::hypot(METxDOWNTES,METyDOWNTES));
  // if (doSVFit) {
  //     if (METxDOWNTES==METx && METyDOWNTES==METy) {
  //         Z2SVPt_TESdn=Z2SVPt;
  //         Z2SVEta_TESdn=Z2SVEta;
  //         Z2SVPhi_TESdn=Z2SVPhi;
  //         Z2SVMass_TESdn=Z2SVMass;
  //     }
  //     else {
  //         SVfit algo_tesdn(0,tau1_tesdn,tau2_tesdn,met_tesdn,covMET,pairType,dm1,dm2);
  //         std::vector<double> results_tesdn=algo_tesdn.FitAndGetResult();
  //         Z2SVPt_TESdn=results_tesdn.at(0);
  //         Z2SVEta_TESdn=results_tesdn.at(1);
  //         Z2SVPhi_TESdn=results_tesdn.at(2);
  //         Z2SVMass_TESdn=results_tesdn.at(3);
  //     }
  // }

  // //EES UP/DOWN
  // TLorentzVector tau1_eesup,tau2_eesup,met_eesup;
  // if (abs(leptons[idx1]->pdgId())==15 && userdatahelpers::getUserInt(leptons[idx1],"isEESShifted")) tau1_eesup.SetPxPyPzE(
  //     userdatahelpers::getUserFloat(leptons[idx1],"px_EleUp"),userdatahelpers::getUserFloat(leptons[idx1],"py_EleUp"),
  //     userdatahelpers::getUserFloat(leptons[idx1],"pz_EleUp"),userdatahelpers::getUserFloat(leptons[idx1],"e_EleUp"));
  // else tau1_eesup.SetPxPyPzE(leptons[idx1]->px(),leptons[idx1]->py(),leptons[idx1]->pz(),leptons[idx1]->energy());
  // if (abs(leptons[idx2]->pdgId())==15 && userdatahelpers::getUserInt(leptons[idx2],"isEESShifted")) tau2_eesup.SetPxPyPzE(
  //     userdatahelpers::getUserFloat(leptons[idx2],"px_EleUp"),userdatahelpers::getUserFloat(leptons[idx2],"py_EleUp"),
  //     userdatahelpers::getUserFloat(leptons[idx2],"pz_EleUp"),userdatahelpers::getUserFloat(leptons[idx2],"e_EleUp"));
  // else tau2_eesup.SetPxPyPzE(leptons[idx2]->px(),leptons[idx2]->py(),leptons[idx2]->pz(),leptons[idx2]->energy());
  // Z2Mass_EESup=(Z2p4-tau1+tau1_eesup-tau2+tau2_eesup).M();
  // met_eesup.SetPxPyPzE(METxUPEES,METyUPEES,0,std::hypot(METxUPEES,METyUPEES));
  // if (doSVFit) {
  //     if (METxUPEES==METx && METyUPEES==METy) {
  //         Z2SVPt_EESup=Z2SVPt;
  //         Z2SVEta_EESup=Z2SVEta;
  //         Z2SVPhi_EESup=Z2SVPhi;
  //         Z2SVMass_EESup=Z2SVMass;
  //     }
  //     else {
  //         SVfit algo_eesup(0,tau1_eesup,tau2_eesup,met_eesup,covMET,pairType,dm1,dm2);
  //         std::vector<double> results_eesup=algo_eesup.FitAndGetResult();
  //         Z2SVPt_EESup=results_eesup.at(0);
  //         Z2SVEta_EESup=results_eesup.at(1);
  //         Z2SVPhi_EESup=results_eesup.at(2);
  //         Z2SVMass_EESup=results_eesup.at(3);
  //     }
  // }

  // TLorentzVector tau1_eesdn,tau2_eesdn,met_eesdn;
  // if (abs(leptons[idx1]->pdgId())==15 && userdatahelpers::getUserInt(leptons[idx1],"isEESShifted")) tau1_eesdn.SetPxPyPzE(
  //     userdatahelpers::getUserFloat(leptons[idx1],"px_EleDown"),userdatahelpers::getUserFloat(leptons[idx1],"py_EleDown"),
  //     userdatahelpers::getUserFloat(leptons[idx1],"pz_EleDown"),userdatahelpers::getUserFloat(leptons[idx1],"e_EleDown"));
  // else tau1_eesdn.SetPxPyPzE(leptons[idx1]->px(),leptons[idx1]->py(),leptons[idx1]->pz(),leptons[idx1]->energy());
  // if (abs(leptons[idx2]->pdgId())==15 && userdatahelpers::getUserInt(leptons[idx2],"isEESShifted")) tau2_eesdn.SetPxPyPzE(
  //     userdatahelpers::getUserFloat(leptons[idx2],"px_EleDown"),userdatahelpers::getUserFloat(leptons[idx2],"py_EleDown"),
  //     userdatahelpers::getUserFloat(leptons[idx2],"pz_EleDown"),userdatahelpers::getUserFloat(leptons[idx2],"e_EleDown"));
  // else tau2_eesdn.SetPxPyPzE(leptons[idx2]->px(),leptons[idx2]->py(),leptons[idx2]->pz(),leptons[idx2]->energy());
  // Z2Mass_EESdn=(Z2p4-tau1+tau1_eesdn-tau2+tau2_eesdn).M();
  // met_eesdn.SetPxPyPzE(METxDOWNEES,METyDOWNEES,0,std::hypot(METxDOWNEES,METyDOWNEES));
  // if (doSVFit) {
  //     if (METxDOWNEES==METx && METyDOWNEES==METy) {
  //         Z2SVPt_EESdn=Z2SVPt;
  //         Z2SVEta_EESdn=Z2SVEta;
  //         Z2SVPhi_EESdn=Z2SVPhi;
  //         Z2SVMass_EESdn=Z2SVMass;
  //     }
  //     else {
  //         SVfit algo_eesdn(0,tau1_eesdn,tau2_eesdn,met_eesdn,covMET,pairType,dm1,dm2);
  //         std::vector<double> results_eesdn=algo_eesdn.FitAndGetResult();
  //         Z2SVPt_EESdn=results_eesdn.at(0);
  //         Z2SVEta_EESdn=results_eesdn.at(1);
  //         Z2SVPhi_EESdn=results_eesdn.at(2);
  //         Z2SVMass_EESdn=results_eesdn.at(3);
  //     }
  // }

  // //MES UP/DOWN
  // TLorentzVector tau1_mesup,tau2_mesup,met_mesup;
  // if (abs(leptons[idx1]->pdgId())==15 && (userdatahelpers::getUserFloat(leptons[idx1],"genmatch")==2 || userdatahelpers::getUserFloat(leptons[idx1],"genmatch")==4)) tau1_mesup.SetPxPyPzE(leptons[idx1]->px()*1.01,leptons[idx1]->py()*1.01,leptons[idx1]->pz()*1.01,leptons[idx1]->energy()*1.01);
  // else tau1_mesup.SetPxPyPzE(leptons[idx1]->px(),leptons[idx1]->py(),leptons[idx1]->pz(),leptons[idx1]->energy());
  // if (abs(leptons[idx2]->pdgId())==15 && (userdatahelpers::getUserFloat(leptons[idx2],"genmatch")==2 || userdatahelpers::getUserFloat(leptons[idx2],"genmatch")==4)) tau2_mesup.SetPxPyPzE(leptons[idx2]->px()*1.01,leptons[idx2]->py()*1.01,leptons[idx2]->pz()*1.01,leptons[idx2]->energy()*1.01);
  // else tau2_mesup.SetPxPyPzE(leptons[idx2]->px(),leptons[idx2]->py(),leptons[idx2]->pz(),leptons[idx2]->energy());
  // Z2Mass_MESup=(Z2p4-tau1+tau1_mesup-tau2+tau2_mesup).M();
  // met_mesup.SetPxPyPzE(METxUPMES,METyUPMES,0,std::hypot(METxUPMES,METyUPMES));
  // if (doSVFit) {
  //     if (METxUPMES==METx && METyUPMES==METy) {
  //         Z2SVPt_MESup=Z2SVPt;
  //         Z2SVEta_MESup=Z2SVEta;
  //         Z2SVPhi_MESup=Z2SVPhi;
  //         Z2SVMass_MESup=Z2SVMass;
  //     }
  //     else {
  //         SVfit algo_mesup(0,tau1_mesup,tau2_mesup,met_mesup,covMET,pairType,dm1,dm2);
  //         std::vector<double> results_mesup=algo_mesup.FitAndGetResult();
  //         Z2SVPt_MESup=results_mesup.at(0);
  //         Z2SVEta_MESup=results_mesup.at(1);
  //         Z2SVPhi_MESup=results_mesup.at(2);
  //         Z2SVMass_MESup=results_mesup.at(3);
  //     }
  // }

  // TLorentzVector tau1_mesdn,tau2_mesdn,met_mesdn;
  // if (abs(leptons[idx1]->pdgId())==15 && (userdatahelpers::getUserFloat(leptons[idx1],"genmatch")==2 || userdatahelpers::getUserFloat(leptons[idx1],"genmatch")==4)) tau1_mesdn.SetPxPyPzE(leptons[idx1]->px()*0.99,leptons[idx1]->py()*0.99,leptons[idx1]->pz()*0.99,leptons[idx1]->energy()*0.99);
  // else tau1_mesdn.SetPxPyPzE(leptons[idx1]->px(),leptons[idx1]->py(),leptons[idx1]->pz(),leptons[idx1]->energy());
  // if (abs(leptons[idx2]->pdgId())==15 && (userdatahelpers::getUserFloat(leptons[idx2],"genmatch")==2 || userdatahelpers::getUserFloat(leptons[idx2],"genmatch")==4)) tau2_mesdn.SetPxPyPzE(leptons[idx2]->px()*0.99,leptons[idx2]->py()*0.99,leptons[idx2]->pz()*0.99,leptons[idx2]->energy()*0.99);
  // else tau2_mesdn.SetPxPyPzE(leptons[idx2]->px(),leptons[idx2]->py(),leptons[idx2]->pz(),leptons[idx2]->energy());
  // Z2Mass_MESdn=(Z2p4-tau1+tau1_mesdn-tau2+tau2_mesdn).M();
  // met_mesdn.SetPxPyPzE(METxDOWNMES,METyDOWNMES,0,std::hypot(METxDOWNMES,METyDOWNMES));
  // if (doSVFit) {
  //     if (METxDOWNMES==METx && METyDOWNMES==METy) {
  //         Z2SVPt_MESdn=Z2SVPt;
  //         Z2SVEta_MESdn=Z2SVEta;
  //         Z2SVPhi_MESdn=Z2SVPhi;
  //         Z2SVMass_MESdn=Z2SVMass;
  //     }
  //     else {
  //         SVfit algo_mesdn(0,tau1_mesdn,tau2_mesdn,met_mesdn,covMET,pairType,dm1,dm2);
  //         std::vector<double> results_mesdn=algo_mesdn.FitAndGetResult();
  //         Z2SVPt_MESdn=results_mesdn.at(0);
  //         Z2SVEta_MESdn=results_mesdn.at(1);
  //         Z2SVPhi_MESdn=results_mesdn.at(2);
  //         Z2SVMass_MESdn=results_mesdn.at(3);
  //     }
  // }

  // Z2Mass_JESup=Z2Mass_JESdn=Z2Mass_JERup=Z2Mass_JERdn=Z2Mass;
  // if (doSVFit) {
  //     //JES UP/DOWN
  //     TLorentzVector tau1_jesup,tau2_jesup,met_jesup;
  //     tau1_jesup.SetPxPyPzE(leptons[idx1]->px(),leptons[idx1]->py(),leptons[idx1]->pz(),leptons[idx1]->energy());
  //     tau2_jesup.SetPxPyPzE(leptons[idx2]->px(),leptons[idx2]->py(),leptons[idx2]->pz(),leptons[idx2]->energy());
  //     met_jesup.SetPxPyPzE(METxUPJES,METyUPJES,0,std::hypot(METxUPJES,METyUPJES));
  //     if (METxUPJES==METx && METyUPJES==METy) {
  //         Z2SVPt_JESup=Z2SVPt;
  //         Z2SVEta_JESup=Z2SVEta;
  //         Z2SVPhi_JESup=Z2SVPhi;
  //         Z2SVMass_JESup=Z2SVMass;
  //     }
  //     else {
  //         SVfit algo_jesup(0,tau1_jesup,tau2_jesup,met_jesup,covMET,pairType,dm1,dm2);
  //         std::vector<double> results_jesup=algo_jesup.FitAndGetResult();
  //         Z2SVPt_JESup=results_jesup.at(0);
  //         Z2SVEta_JESup=results_jesup.at(1);
  //         Z2SVPhi_JESup=results_jesup.at(2);
  //         Z2SVMass_JESup=results_jesup.at(3);
  //     }

  //     TLorentzVector tau1_jesdn,tau2_jesdn,met_jesdn;
  //     tau1_jesdn.SetPxPyPzE(leptons[idx1]->px(),leptons[idx1]->py(),leptons[idx1]->pz(),leptons[idx1]->energy());
  //     tau2_jesdn.SetPxPyPzE(leptons[idx2]->px(),leptons[idx2]->py(),leptons[idx2]->pz(),leptons[idx2]->energy());
  //     met_jesdn.SetPxPyPzE(METxDOWNJES,METyDOWNJES,0,std::hypot(METxDOWNJES,METyDOWNJES));
  //     if (METxDOWNJES==METx && METyDOWNJES==METy) {
  //         Z2SVPt_JESdn=Z2SVPt;
  //         Z2SVEta_JESdn=Z2SVEta;
  //         Z2SVPhi_JESdn=Z2SVPhi;
  //         Z2SVMass_JESdn=Z2SVMass;
  //     }
  //     else {
  //         SVfit algo_jesdn(0,tau1_jesdn,tau2_jesdn,met_jesdn,covMET,pairType,dm1,dm2);
  //         std::vector<double> results_jesdn=algo_jesdn.FitAndGetResult();
  //         Z2SVPt_JESdn=results_jesdn.at(0);
  //         Z2SVEta_JESdn=results_jesdn.at(1);
  //         Z2SVPhi_JESdn=results_jesdn.at(2);
  //         Z2SVMass_JESdn=results_jesdn.at(3);
  //     }

  //     //JER UP/DOWN
  //     TLorentzVector tau1_jerup,tau2_jerup,met_jerup;
  //     tau1_jerup.SetPxPyPzE(leptons[idx1]->px(),leptons[idx1]->py(),leptons[idx1]->pz(),leptons[idx1]->energy());
  //     tau2_jerup.SetPxPyPzE(leptons[idx2]->px(),leptons[idx2]->py(),leptons[idx2]->pz(),leptons[idx2]->energy());
  //     met_jerup.SetPxPyPzE(METxUPJER,METyUPJER,0,std::hypot(METxUPJER,METyUPJER));
  //     if (METxUPJER==METx && METyUPJER==METy) {
  //         Z2SVPt_JERup=Z2SVPt;
  //         Z2SVEta_JERup=Z2SVEta;
  //         Z2SVPhi_JERup=Z2SVPhi;
  //         Z2SVMass_JERup=Z2SVMass;
  //     }
  //     else {
  //         SVfit algo_jerup(0,tau1_jerup,tau2_jerup,met_jerup,covMET,pairType,dm1,dm2);
  //         std::vector<double> results_jerup=algo_jerup.FitAndGetResult();
  //         Z2SVPt_JERup=results_jerup.at(0);
  //         Z2SVEta_JERup=results_jerup.at(1);
  //         Z2SVPhi_JERup=results_jerup.at(2);
  //         Z2SVMass_JERup=results_jerup.at(3);
  //     }

  //     TLorentzVector tau1_jerdn,tau2_jerdn,met_jerdn;
  //     tau1_jerdn.SetPxPyPzE(leptons[idx1]->px(),leptons[idx1]->py(),leptons[idx1]->pz(),leptons[idx1]->energy());
  //     tau2_jerdn.SetPxPyPzE(leptons[idx2]->px(),leptons[idx2]->py(),leptons[idx2]->pz(),leptons[idx2]->energy());
  //     met_jerdn.SetPxPyPzE(METxDOWNJER,METyDOWNJER,0,std::hypot(METxDOWNJER,METyDOWNJER));
  //     if (METxDOWNJER==METx && METyDOWNJER==METy) {
  //         Z2SVPt_JERdn=Z2SVPt;
  //         Z2SVEta_JERdn=Z2SVEta;
  //         Z2SVPhi_JERdn=Z2SVPhi;
  //         Z2SVMass_JERdn=Z2SVMass;
  //     }
  //     else {
  //         SVfit algo_jerdn(0,tau1_jerdn,tau2_jerdn,met_jerdn,covMET,pairType,dm1,dm2);
  //         std::vector<double> results_jerdn=algo_jerdn.FitAndGetResult();
  //         Z2SVPt_JERdn=results_jerdn.at(0);
  //         Z2SVEta_JERdn=results_jerdn.at(1);
  //         Z2SVPhi_JERdn=results_jerdn.at(2);
  //         Z2SVMass_JERdn=results_jerdn.at(3);
  //     }
  // }
  // for (size_t i=0;i<LLSVMass.size();i++)
  //   LLGoodMass.push_back(LLSVMass[i]>0?LLSVMass[i]:LLMass);
  LLGoodMass=LLSVMass>0?LLSVMass:LLMass;
  for (size_t i=0;i<LLSVMass_up.size();i++) {
    LLGoodMass_up.push_back(LLSVMass_up[i]>0?LLSVMass_up[i]:LLMass_up[i]);
    LLGoodMass_dn.push_back(LLSVMass_dn[i]>0?LLSVMass_dn[i]:LLMass_dn[i]);
  }
  // Z2GoodMass_TESup=(Z2SVMass_TESup>0?Z2SVMass_TESup:Z2Mass_TESup);
  // Z2GoodMass_TESdn=(Z2SVMass_TESdn>0?Z2SVMass_TESdn:Z2Mass_TESdn);
  // Z2GoodMass_EESup=(Z2SVMass_EESup>0?Z2SVMass_EESup:Z2Mass_EESup);
  // Z2GoodMass_EESdn=(Z2SVMass_EESdn>0?Z2SVMass_EESdn:Z2Mass_EESdn);
  // Z2GoodMass_MESup=(Z2SVMass_MESup>0?Z2SVMass_MESup:Z2Mass_MESup);
  // Z2GoodMass_MESdn=(Z2SVMass_MESdn>0?Z2SVMass_MESdn:Z2Mass_MESdn);
  // Z2GoodMass_JESup=(Z2SVMass_JESup>0?Z2SVMass_JESup:Z2Mass_JESup);
  // Z2GoodMass_JESdn=(Z2SVMass_JESdn>0?Z2SVMass_JESdn:Z2Mass_JESdn);
  // Z2GoodMass_JERup=(Z2SVMass_JERup>0?Z2SVMass_JERup:Z2Mass_JERup);
  // Z2GoodMass_JERdn=(Z2SVMass_JERdn>0?Z2SVMass_JERdn:Z2Mass_JERdn);

//-------------------------------------------------------------------------------------------------------
//----------------------------------------END SV FIT-----------------------------------------------------
//-------------------------------------------------------------------------------------------------------
  // }
    
  // const reco::Candidate* non_TLE_Z = nullptr;
  // size_t TLE_index = 999;
  // if(abs(LLFlav) == 11*11 || abs(LLFlav) == 13*13) non_TLE_Z = (const reco::Candidate) cand;
  // for (size_t i=0; i<leptons.size(); ++i){
  //   if(abs(leptons[i]->pdgId()) == 22) TLE_index = i;
  // }
  // if(TLE_index < 999 && non_TLE_Z != nullptr) {
  //   TLE_dR_Z = reco::deltaR(non_TLE_Z->p4(), leptons[TLE_index]->p4()); 
  // }

    // Precomputed selections
  //   bool candPass70Z2Loose = cand.userFloat("Z2Mass") &&
  //                            cand.userFloat("MAllComb") &&
  //                            cand.userFloat("pt1")>20 && cand.userFloat("pt2")>10. &&
  //                            ZZMass>70.;
  //   bool candPassFullSel70 = cand.userFloat("SR");
  //   bool candPassFullSel   = cand.userFloat("FullSel");
  //   bool candIsBest = cand.userFloat("isBestCand");
  //   bool passMz_zz = (Z1Mass>60. && Z1Mass<120. && Z2Mass>60. && Z2Mass<120.);   //FIXME hardcoded cut

  //   if (candIsBest) {
  //     //    sel = 10; //FIXME see above
  //     if (candPass70Z2Loose) sel=70;
  //     if (candPassFullSel70){ // includes MZ2 > 12
  //       sel = 90;
  //       if (candPassFullSel){
  //         sel=100;
  //         if (passMz_zz) sel = 120;
  //       }
  //     }
  //   }
  // } else if(theChannel==ZLL) sel = 20; 
  // else if(theChannel==ZL) sel = 10;
 

  int sel = 1;
  if (!(evtPass)) {sel = -sel;} // avoid confusion when we write events which do not pass trigger/skim

  // Retrieve the userFloat of the leptons in vectors ordered in the same way.
  vector<float> SIP(2);
  vector<float> combRelIsoPF(2);
  passIsoPreFSR = true;

  for (unsigned int i=0; i<leptons.size(); ++i){
    // float curr_dR = 999;
    // if(i != TLE_index && TLE_index < 999)
    //   curr_dR = reco::deltaR(leptons[i]->p4(), leptons[TLE_index]->p4());
    // if(curr_dR < TLE_min_dR_3l) TLE_min_dR_3l = curr_dR;

    short lepFlav = std::abs(leptons[i]->pdgId());

    SIP[i]             = userdatahelpers::getUserFloat(leptons[i],"SIP");
    passIsoPreFSR      = passIsoPreFSR&&(userdatahelpers::getUserFloat(leptons[i],"combRelIsoPF")<LeptonIsoHelper::isoCut(leptons[i]));

    //in the Legacy approach,  FSR-corrected iso is attached to the Z, not to the lepton!
    combRelIsoPF[i] = cand.userFloat(labels[i]+"combRelIsoPFFSRCorr"); // Note: the
    assert(SIP[i] == cand.userFloat(labels[i]+"SIP")); // Check that I don't mess up with labels[] and leptons[]

    //Fill the info on the lepton candidates
    LepPt .push_back( leptons[i]->pt() );
    LepEta.push_back( leptons[i]->eta() );
    LepPhi.push_back( leptons[i]->phi() );
    LepM.push_back( leptons[i]->mass() );
    LepSCEta.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"SCeta") : -99. );
    int id =  leptons[i]->pdgId();
    if(id == 22 && i == 1) id=-22; //FIXME this assumes a standard ordering of leptons.
    LepLepId.push_back( id );
    LepSIP  .push_back( SIP[i] );
    Lepdxy  .push_back( userdatahelpers::getUserFloat(leptons[i],"dxy") );
    Lepdz   .push_back( userdatahelpers::getUserFloat(leptons[i],"dz") );
    LepTime .push_back( lepFlav==13 ? userdatahelpers::getUserFloat(leptons[i],"time") : 0. );
    LepisID .push_back( userdatahelpers::getUserFloat(leptons[i],"isGood") );
    // LepBDT  .push_back( (lepFlav==13 || lepFlav==11) ? userdatahelpers::getUserFloat(leptons[i],"BDT") : -99. );
    LepisCrack.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"isCrack") : 0 );
    LepMissingHit.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"missingHit") : 0 );
    LepConversionVeto.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"ConversionVeto") : true);
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
      TauVSmu.push_back(userdatahelpers::getUserInt(leptons[i],"tauVSmu"));
      TauVSe.push_back(userdatahelpers::getUserInt(leptons[i],"tauVSe"));
      TauVSjet.push_back(userdatahelpers::getUserInt(leptons[i],"tauVSjet"));

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
    } else {
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

    // HLTMatch1.push_back( userdatahelpers::hasUserFloat(leptons[i],"HLTMatch1") ? userdatahelpers::getUserFloat(leptons[i],"HLTMatch1") : -1 );
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

  //Compute the data/MC weight
  dataMCWeight = 1.;
  //When the trigger is not applied in the MC, apply a trigger efficiency factor instead (FIXME: here hardcoding the efficiencies computed for ICHEP2016)
  if(isMC) {
    dataMCWeight = getAllWeight(leptons);
  }
  //Store an overall event weight (which is normalized by gen_sumWeights)
  overallEventWeight = PUWeight * genHEPMCweight * dataMCWeight * trigEffWeight;
}

// void LLNtupleMaker::FillGENCandidate(const pat::CompositeCandidate& cand){ //AT
//   passedFiducialSelection_bbf = cand.userData<bool>("passedFiducial");
//   if(cand.userFloat("event") == 91257) {
//     cout << passedFiducialSelection_bbf << endl;
//   }
// }


void LLNtupleMaker::getCheckedUserFloat(const pat::CompositeCandidate& cand, const std::string& strval, Float_t& setval, Float_t defaultval){
  if (cand.hasUserFloat(strval)) setval = cand.userFloat(strval);
  else setval = defaultval;
}



// ------------ method called once each job just before starting event loop  ------------
void LLNtupleMaker::beginJob()
{
  edm::Service<TFileService> fs;
  TTree *candTree = fs->make<TTree>(theFileName,"Event Summary");
  TTree *candTree_failed = 0;
  if (failedTreeLevel)
    candTree_failed = fs->make<TTree>(theFileName+"_failed","Event Summary");
  myTree = new LLNtupleFactory(candTree, candTree_failed);
  const int nbins = 45;
  hCounter = fs->make<TH1F>("Counters", "Counters", nbins, 0., nbins);
  BookAllBranches();
}

// ------------ method called once each job just after ending the event loop  ------------
void LLNtupleMaker::endJob()
{
  hCounter->SetBinContent(0 ,gen_sumWeights); // also stored in bin 40
  hCounter->SetBinContent(1 ,Nevt_Gen-gen_BUGGY);
  hCounter->SetBinContent(2 ,gen_mumu);
  hCounter->SetBinContent(3 ,gen_mumu_LeptonAcceptance);
  hCounter->SetBinContent(4 ,gen_mumu_EtaAcceptance);
  hCounter->SetBinContent(5 ,gen_etau);
  hCounter->SetBinContent(6 ,gen_etau_LeptonAcceptance);
  hCounter->SetBinContent(7 ,gen_etau_EtaAcceptance);
  hCounter->SetBinContent(8 ,gen_mutau);
  hCounter->SetBinContent(9 ,gen_mutau_LeptonAcceptance);
  hCounter->SetBinContent(10,gen_mutau_EtaAcceptance);
  hCounter->SetBinContent(11,gen_tautau);
  hCounter->SetBinContent(12,gen_tautau_LeptonAcceptance);
  hCounter->SetBinContent(13,gen_tautau_EtaAcceptance);
  hCounter->SetBinContent(14,gen_emu);
  hCounter->SetBinContent(15,gen_emu_EtaAcceptance);
  hCounter->SetBinContent(16,gen_emu_LeptonAcceptance);
  hCounter->SetBinContent(19,gen_BUGGY);
  hCounter->SetBinContent(20,gen_Unknown);

  hCounter->SetBinContent(40,gen_sumWeights);
  hCounter->SetBinContent(41,gen_sumGenMCWeight);
  hCounter->SetBinContent(42,gen_sumPUWeight);

  TH1 *h[1] ={ hCounter };
  for (int i = 0; i < 1; i++) {
    h[i]->GetXaxis()->SetBinLabel(1 ,"Nevt_Gen");
    h[i]->GetXaxis()->SetBinLabel(2 ,"gen_mumu");
    h[i]->GetXaxis()->SetBinLabel(3 ,"gen_mumu_LeptonAcceptance");
    h[i]->GetXaxis()->SetBinLabel(4 ,"gen_mumu_EtaAcceptance");
    h[i]->GetXaxis()->SetBinLabel(5 ,"gen_etau");
    h[i]->GetXaxis()->SetBinLabel(6 ,"gen_etau_LeptonAcceptance");
    h[i]->GetXaxis()->SetBinLabel(7 ,"gen_etau_EtaAcceptance");
    h[i]->GetXaxis()->SetBinLabel(8 ,"gen_mutau");
    h[i]->GetXaxis()->SetBinLabel(9 ,"gen_mutau_LeptonAcceptance");
    h[i]->GetXaxis()->SetBinLabel(10,"gen_mutau_EtaAcceptance");
    h[i]->GetXaxis()->SetBinLabel(11,"gen_tautau");
    h[i]->GetXaxis()->SetBinLabel(12,"gen_tautau_LeptonAcceptance");
    h[i]->GetXaxis()->SetBinLabel(13,"gen_tautau_EtaAcceptance");
    h[i]->GetXaxis()->SetBinLabel(14,"gen_emu");
    h[i]->GetXaxis()->SetBinLabel(15,"gen_emu_EtaAcceptance");
    h[i]->GetXaxis()->SetBinLabel(16,"gen_emu_LeptonAcceptance");
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
void LLNtupleMaker::beginRun(edm::Run const& iRun, edm::EventSetup const&)
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
void LLNtupleMaker::endRun(edm::Run const& iRun, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void LLNtupleMaker::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
  Nevt_Gen_lumiBlock = 0;
}

// ------------ method called when ending the processing of a luminosity block  ------------
void LLNtupleMaker::endLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iSetup)
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
void LLNtupleMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

Float_t LLNtupleMaker::getAllWeight(const vector<const reco::Candidate*>& leptons) 
{
  Float_t totWeight = 1.;
  int flav=abs(leptons[0]->pdgId()*leptons[1]->pdgId());
	
  for(unsigned int i=0; i<leptons.size(); ++i){
    Int_t   myLepID = abs(leptons[i]->pdgId());
    if (skipMuDataMCWeight && myLepID==13) return 1.;
    if (skipEleDataMCWeight && myLepID==11) return 1.;
    if (skipTauDataMCWeight && myLepID==15) return 1.;
    
    float SF = 1.0;
    //ele
    float SF_UncUp = 1.0;
    float SF_UncDn = 1.0;
    //mu
    float SF_UncUp_RECO_syst = 1.0;
    float SF_UncUp_RECO_stat = 1.0;
    float SF_UncUp_ID_syst = 1.0;
    float SF_UncUp_ID_stat = 1.0;
    float SF_UncUp_ISO_syst = 1.0;
    float SF_UncUp_ISO_stat = 1.0;
    float SF_UncDn_RECO_syst = 1.0;
    float SF_UncDn_RECO_stat = 1.0;
    float SF_UncDn_ID_syst = 1.0;
    float SF_UncDn_ID_stat = 1.0;
    float SF_UncDn_ISO_syst = 1.0;
    float SF_UncDn_ISO_stat = 1.0;
    //tau
    float SF_UncUp_uncert0 = 1.0;
    float SF_UncUp_uncert1 = 1.0;
    float SF_UncUp_syst_alleras = 1.0;
    float SF_UncUp_syst_year = 1.0;
    float SF_UncUp_syst_dm_year = 1.0;
    float SF_UncUp_fakeEle = 1.0;
    float SF_UncUp_fakeMu = 1.0;
    float SF_UncDn_uncert0 = 1.0;
    float SF_UncDn_uncert1 = 1.0;
    float SF_UncDn_syst_alleras = 1.0;
    float SF_UncDn_syst_year = 1.0;
    float SF_UncDn_syst_dm_year = 1.0;
    float SF_UncDn_fakeEle = 1.0;
    float SF_UncDn_fakeMu = 1.0;


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
     
    // bool isCrack;
    // if (myLepID == 11) isCrack = userdatahelpers::getUserFloat(leptons[i],"isCrack");
    // else isCrack = false;

    if (myLepID == 11) {
      SF = lepSFHelper->getSF(year, myLepID, myLepPt, myLepEta, mySCeta, "", "");
      SF_UncUp = SF + lepSFHelper->getSF(year, myLepID, myLepPt, myLepEta, mySCeta, "unc", "");
      SF_UncDn = 2*SF - SF_UncUp;
    }
    else if (myLepID == 13) {
      SF = lepSFHelper->getSF(year, myLepID, myLepPt, myLepEta, mySCeta, "", "");
      SF_UncUp_RECO_syst = SF + lepSFHelper->getSF(year, myLepID, myLepPt, myLepEta, mySCeta, "unc", "RECO_syst");
      SF_UncDn_RECO_syst = 2*SF - SF_UncUp_RECO_syst;
      SF_UncUp_RECO_stat = SF + lepSFHelper->getSF(year, myLepID, myLepPt, myLepEta, mySCeta, "unc", "RECO_stat");
      SF_UncDn_RECO_stat = 2*SF - SF_UncUp_RECO_stat;
      SF_UncUp_ID_syst = SF + lepSFHelper->getSF(year, myLepID, myLepPt, myLepEta, mySCeta, "unc", "ID_syst");
      SF_UncDn_ID_syst = 2*SF - SF_UncUp_ID_syst;
      SF_UncUp_ID_stat = SF + lepSFHelper->getSF(year, myLepID, myLepPt, myLepEta, mySCeta, "unc", "ID_stat");
      SF_UncDn_ID_stat = 2*SF - SF_UncUp_ID_stat;
      SF_UncUp_ISO_syst = SF + lepSFHelper->getSF(year, myLepID, myLepPt, myLepEta, mySCeta, "unc", "ISO_syst");
      SF_UncDn_ISO_syst = 2*SF - SF_UncUp_ISO_syst;
      SF_UncUp_ISO_stat = SF + lepSFHelper->getSF(year, myLepID, myLepPt, myLepEta, mySCeta, "unc", "ISO_stat");
      SF_UncDn_ISO_stat = 2*SF - SF_UncUp_ISO_stat;
    }
    else {
      int gm=userdatahelpers::getUserFloat(leptons[i],"genmatch");
      int dm= userdatahelpers::getUserFloat(leptons[i],"decayMode");


      if (gm==5) {
        if (flav==165) {
          SF=DeepTauSF_VSjet_ETau->getSFvsDMandPT(myLepPt,dm,gm);
          SF_UncUp_uncert0=DeepTauSF_VSjet_ETau->getSFvsDMandPT(myLepPt,dm,gm,"uncert0_up");
          SF_UncUp_uncert1=DeepTauSF_VSjet_ETau->getSFvsDMandPT(myLepPt,dm,gm,"uncert1_up");
          SF_UncUp_syst_alleras=DeepTauSF_VSjet_ETau->getSFvsDMandPT(myLepPt,dm,gm,"syst_alleras_up");
          SF_UncUp_syst_year=DeepTauSF_VSjet_ETau->getSFvsDMandPT(myLepPt,dm,gm,"syst_"+std::to_string(year)+(year==2016?(preVFP?"_preVFP":"_postVFP"):"")+"_up");
          SF_UncUp_syst_dm_year=DeepTauSF_VSjet_ETau->getSFvsDMandPT(myLepPt,dm,gm,"syst_dm"+std::to_string(dm)+"_"+std::to_string(year)+(year==2016?(preVFP?"_preVFP":"_postVFP"):"")+"_up");
          SF_UncDn_uncert0=DeepTauSF_VSjet_ETau->getSFvsDMandPT(myLepPt,dm,gm,"uncert0_down");
          SF_UncDn_uncert1=DeepTauSF_VSjet_ETau->getSFvsDMandPT(myLepPt,dm,gm,"uncert1_down");
          SF_UncDn_syst_alleras=DeepTauSF_VSjet_ETau->getSFvsDMandPT(myLepPt,dm,gm,"syst_alleras_down");
          SF_UncDn_syst_year=DeepTauSF_VSjet_ETau->getSFvsDMandPT(myLepPt,dm,gm,"syst_"+std::to_string(year)+(year==2016?(preVFP?"_preVFP":"_postVFP"):"")+"_down");
          SF_UncDn_syst_dm_year=DeepTauSF_VSjet_ETau->getSFvsDMandPT(myLepPt,dm,gm,"syst_dm"+std::to_string(dm)+"_"+std::to_string(year)+(year==2016?(preVFP?"_preVFP":"_postVFP"):"")+"_down");

        }
        else if (flav==195) {
          SF=DeepTauSF_VSjet_MuTau->getSFvsDMandPT(myLepPt,dm,gm);
          SF_UncUp_uncert0=DeepTauSF_VSjet_MuTau->getSFvsDMandPT(myLepPt,dm,gm,"uncert0_up");
          SF_UncUp_uncert1=DeepTauSF_VSjet_MuTau->getSFvsDMandPT(myLepPt,dm,gm,"uncert1_up");
          SF_UncUp_syst_alleras=DeepTauSF_VSjet_MuTau->getSFvsDMandPT(myLepPt,dm,gm,"syst_alleras_up");
          SF_UncUp_syst_year=DeepTauSF_VSjet_MuTau->getSFvsDMandPT(myLepPt,dm,gm,"syst_"+std::to_string(year)+(year==2016?(preVFP?"_preVFP":"_postVFP"):"")+"_up");
          SF_UncUp_syst_dm_year=DeepTauSF_VSjet_MuTau->getSFvsDMandPT(myLepPt,dm,gm,"syst_dm"+std::to_string(dm)+"_"+std::to_string(year)+(year==2016?(preVFP?"_preVFP":"_postVFP"):"")+"_up");
          SF_UncDn_uncert0=DeepTauSF_VSjet_MuTau->getSFvsDMandPT(myLepPt,dm,gm,"uncert0_down");
          SF_UncDn_uncert1=DeepTauSF_VSjet_MuTau->getSFvsDMandPT(myLepPt,dm,gm,"uncert1_down");
          SF_UncDn_syst_alleras=DeepTauSF_VSjet_MuTau->getSFvsDMandPT(myLepPt,dm,gm,"syst_alleras_down");
          SF_UncDn_syst_year=DeepTauSF_VSjet_MuTau->getSFvsDMandPT(myLepPt,dm,gm,"syst_"+std::to_string(year)+(year==2016?(preVFP?"_preVFP":"_postVFP"):"")+"_down");
          SF_UncDn_syst_dm_year=DeepTauSF_VSjet_MuTau->getSFvsDMandPT(myLepPt,dm,gm,"syst_dm"+std::to_string(dm)+"_"+std::to_string(year)+(year==2016?(preVFP?"_preVFP":"_postVFP"):"")+"_down");
        }
        else if (flav==225) {
          SF=DeepTauSF_VSjet_TauTau->getSFvsDMandPT(myLepPt,dm,gm);
          SF_UncUp_uncert0=DeepTauSF_VSjet_TauTau->getSFvsDMandPT(myLepPt,dm,gm,"uncert0_up");
          SF_UncUp_uncert1=DeepTauSF_VSjet_TauTau->getSFvsDMandPT(myLepPt,dm,gm,"uncert1_up");
          SF_UncUp_syst_alleras=DeepTauSF_VSjet_TauTau->getSFvsDMandPT(myLepPt,dm,gm,"syst_alleras_up");
          SF_UncUp_syst_year=DeepTauSF_VSjet_TauTau->getSFvsDMandPT(myLepPt,dm,gm,"syst_"+std::to_string(year)+(year==2016?(preVFP?"_preVFP":"_postVFP"):"")+"_up");
          SF_UncUp_syst_dm_year=DeepTauSF_VSjet_TauTau->getSFvsDMandPT(myLepPt,dm,gm,"syst_dm"+std::to_string(dm)+"_"+std::to_string(year)+(year==2016?(preVFP?"_preVFP":"_postVFP"):"")+"_up");
          SF_UncDn_uncert0=DeepTauSF_VSjet_TauTau->getSFvsDMandPT(myLepPt,dm,gm,"uncert0_down");
          SF_UncDn_uncert1=DeepTauSF_VSjet_TauTau->getSFvsDMandPT(myLepPt,dm,gm,"uncert1_down");
          SF_UncDn_syst_alleras=DeepTauSF_VSjet_TauTau->getSFvsDMandPT(myLepPt,dm,gm,"syst_alleras_down");
          SF_UncDn_syst_year=DeepTauSF_VSjet_TauTau->getSFvsDMandPT(myLepPt,dm,gm,"syst_"+std::to_string(year)+(year==2016?(preVFP?"_preVFP":"_postVFP"):"")+"_down");
          SF_UncDn_syst_dm_year=DeepTauSF_VSjet_TauTau->getSFvsDMandPT(myLepPt,dm,gm,"syst_dm"+std::to_string(dm)+"_"+std::to_string(year)+(year==2016?(preVFP?"_preVFP":"_postVFP"):"")+"_down");
        }
      }
      else if (gm==1 || gm==3) {
        if (flav==165) {
          SF=DeepTauSF_VSe_ETau->getSFvsEta(myLepEta,gm);
          SF_UncUp_fakeEle=DeepTauSF_VSe_ETau->getSFvsEta(myLepEta,gm,"Up");
          SF_UncDn_fakeEle=DeepTauSF_VSe_ETau->getSFvsEta(myLepEta,gm,"Down");
        }
        else if (flav==195) {
          SF=DeepTauSF_VSe_MuTau->getSFvsEta(myLepEta,gm);
          SF_UncUp_fakeEle=DeepTauSF_VSe_MuTau->getSFvsEta(myLepEta,gm,"Up");
          SF_UncDn_fakeEle=DeepTauSF_VSe_MuTau->getSFvsEta(myLepEta,gm,"Dn");
        }
        else if (flav==225) {
          SF=DeepTauSF_VSe_TauTau->getSFvsEta(myLepEta,gm);
          SF_UncUp_fakeEle=DeepTauSF_VSe_TauTau->getSFvsEta(myLepEta,gm,"Up");
          SF_UncDn_fakeEle=DeepTauSF_VSe_TauTau->getSFvsEta(myLepEta,gm,"Dn");
        }
      }
      else if (gm==2 || gm==4) {
        if (flav==165) { 
          SF=DeepTauSF_VSmu_ETau->getSFvsEta(myLepEta,gm);
          SF_UncUp_fakeMu=DeepTauSF_VSmu_ETau->getSFvsEta(myLepEta,gm,"Up");
          SF_UncDn_fakeMu=DeepTauSF_VSmu_ETau->getSFvsEta(myLepEta,gm,"Dn");
        }
        else if (flav==195) {
          SF=DeepTauSF_VSmu_MuTau->getSFvsEta(myLepEta,gm);
          SF_UncUp_fakeMu=DeepTauSF_VSmu_MuTau->getSFvsEta(myLepEta,gm,"Up");
          SF_UncDn_fakeMu=DeepTauSF_VSmu_MuTau->getSFvsEta(myLepEta,gm,"Dn");
        }
        else if (flav==225) {
          SF=DeepTauSF_VSmu_TauTau->getSFvsEta(myLepEta,gm);
          SF_UncUp_fakeMu=DeepTauSF_VSmu_TauTau->getSFvsEta(myLepEta,gm,"Up");
          SF_UncDn_fakeMu=DeepTauSF_VSmu_TauTau->getSFvsEta(myLepEta,gm,"Dn");
        }
      }

    }
    LepSF.push_back(SF);
    LepSF_UncUp.push_back(SF_UncUp);
    LepSF_UncDn.push_back(SF_UncDn);

    LepSF_UncUp_RECO_syst.push_back(SF_UncUp_RECO_syst);
    LepSF_UncDn_RECO_syst.push_back(SF_UncDn_RECO_syst);
    LepSF_UncUp_RECO_stat.push_back(SF_UncUp_RECO_stat);
    LepSF_UncDn_RECO_stat.push_back(SF_UncDn_RECO_stat);
    LepSF_UncUp_ID_syst.push_back(SF_UncUp_ID_syst);
    LepSF_UncDn_ID_syst.push_back(SF_UncDn_ID_syst);
    LepSF_UncUp_ID_stat.push_back(SF_UncUp_ID_stat);
    LepSF_UncDn_ID_stat.push_back(SF_UncDn_ID_stat);
    LepSF_UncUp_ISO_syst.push_back(SF_UncUp_ISO_syst);
    LepSF_UncDn_ISO_syst.push_back(SF_UncDn_ISO_syst);
    LepSF_UncUp_ISO_stat.push_back(SF_UncUp_ISO_stat);
    LepSF_UncDn_ISO_stat.push_back(SF_UncDn_ISO_stat);

    LepSF_UncUp_uncert0.push_back(SF_UncUp_uncert0);
    LepSF_UncUp_uncert1.push_back(SF_UncUp_uncert1);
    LepSF_UncUp_syst_alleras.push_back(SF_UncUp_syst_alleras);
    LepSF_UncUp_syst_year.push_back(SF_UncUp_syst_year);
    LepSF_UncUp_syst_dm_year.push_back(SF_UncUp_syst_dm_year);
    LepSF_UncUp_fakeEle.push_back(SF_UncUp_fakeEle);
    LepSF_UncUp_fakeMu.push_back(SF_UncUp_fakeMu);
    LepSF_UncDn_uncert0.push_back(SF_UncDn_uncert0);
    LepSF_UncDn_uncert1.push_back(SF_UncDn_uncert1);
    LepSF_UncDn_syst_alleras.push_back(SF_UncDn_syst_alleras);
    LepSF_UncDn_syst_year.push_back(SF_UncDn_syst_year);
    LepSF_UncDn_syst_dm_year.push_back(SF_UncDn_syst_dm_year);
    LepSF_UncDn_fakeEle.push_back(SF_UncDn_fakeEle);
    LepSF_UncDn_fakeMu.push_back(SF_UncDn_fakeMu);
    
    totWeight *= SF;
  } 

  return totWeight;
}


// Float_t LLNtupleMaker::getHqTWeight(double mH, double genPt) const
// {
//   if (skipHqTWeight) return 1.;

//   //cout<<"mH = "<<mH<<", genPt = "<<genPt<<endl;
//   if (mH<400 || genPt>250) return 1.;

//   double weight = 1.;

//   const int masses[4] = {400,600,800,1000};
//   double massDiff = 1000;
//   int iMass = -1;
//   for (int i=0; i<4; ++i){
//     double massDiffTmp = std::fabs(mH-masses[i]);
//     if (massDiffTmp<massDiff){
//       massDiff = massDiffTmp;
//       iMass = i;
//     }
//   }

//   if (iMass>=0) {
//     weight = h_weight->GetBinContent(h_weight->FindBin(genPt));
//   }
//   return weight;
// }


void LLNtupleMaker::FillLLGenInfo(Short_t ZId, const math::XYZTLorentzVector pZ)
{
  GenLLMass= pZ.M();
  GenLLPt= pZ.Pt();
  GenLLEta= pZ.Eta();
  GenLLPhi= pZ.Phi();
  GenLLFlav= ZId;

  return;
}

void LLNtupleMaker::FillVisLLGenInfo(Short_t ZId, const math::XYZTLorentzVector pZ)
{
  GenVisLLMass= pZ.M();
  GenVisLLPt= pZ.Pt();
  GenVisLLEta= pZ.Eta();
  GenVisLLPhi= pZ.Phi();
  GenVisLLFlav= ZId;

  return;
}

void LLNtupleMaker::FillLepGenInfo(Short_t Lep1Id, Short_t Lep2Id, const math::XYZTLorentzVector Lep1, const math::XYZTLorentzVector Lep2)
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
    
  return;
}

void LLNtupleMaker::FillVisLepGenInfo(Short_t Lep1Id, Short_t Lep2Id, const math::XYZTLorentzVector Lep1, const math::XYZTLorentzVector Lep2)
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

  return;
}



void LLNtupleMaker::FillLepGenIso(float_t Lep1Iso, float_t Lep2Iso)
{
  GenLep1Iso = Lep1Iso;
  GenLep2Iso = Lep2Iso;

  return;
}//AT

void LLNtupleMaker::BookAllBranches(){
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
  myTree->Book("GenMETx", GenMETx, failedTreeLevel >= minimalFailedTree);
  myTree->Book("GenMETy", GenMETy, failedTreeLevel >= minimalFailedTree);
  myTree->Book("PFMET", PFMET, failedTreeLevel >= fullFailedTree);
  myTree->Book("PFMETPhi", PFMETPhi, failedTreeLevel >= fullFailedTree);
  myTree->Book("METx", METx, failedTreeLevel >= fullFailedTree);
  myTree->Book("METy", METy, failedTreeLevel >= fullFailedTree);
  // myTree->Book("PFMETTau", PFMETTau, failedTreeLevel >= fullFailedTree);
  // myTree->Book("PFMETPhiTau", PFMETPhiTau, failedTreeLevel >= fullFailedTree);
  // myTree->Book("METxTau", METxTau, failedTreeLevel >= fullFailedTree);
  // myTree->Book("METyTau", METyTau, failedTreeLevel >= fullFailedTree);
  // myTree->Book("PFMETJet", PFMETJet, failedTreeLevel >= fullFailedTree);
  // myTree->Book("PFMETPhiJet", PFMETPhiJet, failedTreeLevel >= fullFailedTree);
  // myTree->Book("METxJet", METxJet, failedTreeLevel >= fullFailedTree);
  // myTree->Book("METyJet", METyJet, failedTreeLevel >= fullFailedTree);
  // myTree->Book("PFMETRaw", PFMETRaw, failedTreeLevel >= fullFailedTree);
  // myTree->Book("PFMETPhiRaw", PFMETPhiRaw, failedTreeLevel >= fullFailedTree);
  // myTree->Book("METxRaw", METxRaw, failedTreeLevel >= fullFailedTree);
  // myTree->Book("METyRaw", METyRaw, failedTreeLevel >= fullFailedTree);
  myTree->Book("PFMETRecoil", PFMETRecoil, failedTreeLevel >= fullFailedTree);
  myTree->Book("PFMETPhiRecoil", PFMETPhiRecoil, failedTreeLevel >= fullFailedTree);
  myTree->Book("METxRecoil", METxRecoil, failedTreeLevel >= fullFailedTree);
  myTree->Book("METyRecoil", METyRecoil, failedTreeLevel >= fullFailedTree);
  myTree->Book("Pzeta1", Pzeta1, failedTreeLevel >= fullFailedTree);
  myTree->Book("Pzeta2", Pzeta2, failedTreeLevel >= fullFailedTree);
  myTree->Book("MtLMET", MtLMET, failedTreeLevel >= fullFailedTree);

  myTree->Book("DeltaEtaJJ", DeltaEtaJJ, failedTreeLevel >= fullFailedTree);
  myTree->Book("DiJetMass", DiJetMass, failedTreeLevel >= fullFailedTree);
  myTree->Book("VBFJetIdx1", VBFJetIdx1, failedTreeLevel >= fullFailedTree);
  myTree->Book("VBFJetIdx2", VBFJetIdx2, failedTreeLevel >= fullFailedTree);

  // myTree->Book("METxUPTES", METxUPTES, failedTreeLevel >= fullFailedTree);
  // myTree->Book("METyUPTES", METyUPTES, failedTreeLevel >= fullFailedTree);
  // myTree->Book("METxDOWNTES", METxDOWNTES, failedTreeLevel >= fullFailedTree);
  // myTree->Book("METyDOWNTES", METyDOWNTES, failedTreeLevel >= fullFailedTree);
  // myTree->Book("METxUPEES", METxUPEES, failedTreeLevel >= fullFailedTree);
  // myTree->Book("METyUPEES", METyUPEES, failedTreeLevel >= fullFailedTree);
  // myTree->Book("METxDOWNEES", METxDOWNEES, failedTreeLevel >= fullFailedTree);
  // myTree->Book("METyDOWNEES", METyDOWNEES, failedTreeLevel >= fullFailedTree);
  // myTree->Book("METxUPMES", METxUPMES, failedTreeLevel >= fullFailedTree);
  // myTree->Book("METyUPMES", METyUPMES, failedTreeLevel >= fullFailedTree);
  // myTree->Book("METxDOWNMES", METxDOWNMES, failedTreeLevel >= fullFailedTree);
  // myTree->Book("METyDOWNMES", METyDOWNMES, failedTreeLevel >= fullFailedTree);
  // myTree->Book("METxUPJES", METxUPJES, failedTreeLevel >= fullFailedTree);
  // myTree->Book("METyUPJES", METyUPJES, failedTreeLevel >= fullFailedTree);
  // myTree->Book("METxDOWNJES", METxDOWNJES, failedTreeLevel >= fullFailedTree);
  // myTree->Book("METyDOWNJES", METyDOWNJES, failedTreeLevel >= fullFailedTree);
  // myTree->Book("METxUPJER", METxUPJER, failedTreeLevel >= fullFailedTree);
  // myTree->Book("METyUPJER", METyUPJER, failedTreeLevel >= fullFailedTree);
  // myTree->Book("METxDOWNJER", METxDOWNJER, failedTreeLevel >= fullFailedTree);
  // myTree->Book("METyDOWNJER", METyDOWNJER, failedTreeLevel >= fullFailedTree);
    
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
  myTree->Book("nCleanedJetsPt25BTagged_bTagSF",nCleanedJetsPt25BTagged_bTagSF, failedTreeLevel >= minimalFailedTree);
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
  myTree->Book("metWord",metWord, failedTreeLevel >= minimalFailedTree);

  // myTree->Book("evtPassMETFilter",evtPassMETTrigger, failedTreeLevel >= minimalFailedTree);
//   myTree->Book("ZZGoodMass",ZZGoodMass, false);
  myTree->Book("pass_SingleTrigger",pass_SingleTrigger, false);
  myTree->Book("pass_CrossTrigger",pass_CrossTrigger, false);
  myTree->Book("pass_Trigger",pass_Trigger, false);

  myTree->Book("LLMass",LLMass, false);
//   myTree->Book("Z2GoodMass",Z2GoodMass, false);
  myTree->Book("LLPt",LLPt, false);
  myTree->Book("LLEta",LLEta, false);
  myTree->Book("LLPhi",LLPhi, false);
  myTree->Book("LLDR",LLDR, false);
  myTree->Book("LLFlav",LLFlav, false);
  myTree->Book("LLisGoodTau",LLisGoodTau, false);

  if (addSVfit){
    myTree->Book("LLSVMass",LLSVMass, false);
    myTree->Book("LLSVPt",LLSVPt, false);
    myTree->Book("LLSVEta",LLSVEta, false);
    myTree->Book("LLSVPhi",LLSVPhi, false);
    myTree->Book("LLGoodMass",LLGoodMass, false);

    if (theChannel==SR) {
      myTree->Book("LLSVMass_up",LLSVMass_up, false);
      myTree->Book("LLSVPt_up",LLSVPt_up, false);
      myTree->Book("LLSVEta_up",LLSVEta_up, false);
      myTree->Book("LLSVPhi_up",LLSVPhi_up, false);
      myTree->Book("LLGoodMass_up",LLGoodMass_up, false);
      myTree->Book("LLSVMass_dn",LLSVMass_dn, false);
      myTree->Book("LLSVPt_dn",LLSVPt_dn, false);
      myTree->Book("LLSVEta_dn",LLSVEta_dn, false);
      myTree->Book("LLSVPhi_dn",LLSVPhi_dn, false);
      myTree->Book("LLGoodMass_dn",LLGoodMass_dn, false);
      // myTree->Book("Names_shift",Names_shift, false);
    }
      
    // myTree->Book("Z2Mass_TESup",Z2Mass_TESup, false);
    // myTree->Book("Z2SVMass_TESup",Z2SVMass_TESup, false);
    // myTree->Book("Z2SVPt_TESup",Z2SVPt_TESup, false);
    // myTree->Book("Z2SVEta_TESup",Z2SVEta_TESup, false);
    // myTree->Book("Z2SVPhi_TESup",Z2SVPhi_TESup, false);
    // myTree->Book("Z2GoodMass_TESup",Z2GoodMass_TESup, false);

    // myTree->Book("Z2Mass_TESdn",Z2Mass_TESdn, false);
    // myTree->Book("Z2SVMass_TESdn",Z2SVMass_TESdn, false);
    // myTree->Book("Z2SVPt_TESdn",Z2SVPt_TESdn, false);
    // myTree->Book("Z2SVEta_TESdn",Z2SVEta_TESdn, false);
    // myTree->Book("Z2SVPhi_TESdn",Z2SVPhi_TESdn, false);
    // myTree->Book("Z2GoodMass_TESdn",Z2GoodMass_TESdn, false);

    // myTree->Book("Z2Mass_EESup",Z2Mass_EESup, false);
    // myTree->Book("Z2SVMass_EESup",Z2SVMass_EESup, false);
    // myTree->Book("Z2SVPt_EESup",Z2SVPt_EESup, false);
    // myTree->Book("Z2SVEta_EESup",Z2SVEta_EESup, false);
    // myTree->Book("Z2SVPhi_EESup",Z2SVPhi_EESup, false);
    // myTree->Book("Z2GoodMass_EESup",Z2GoodMass_EESup, false);

    // myTree->Book("Z2Mass_EESdn",Z2Mass_EESdn, false);
    // myTree->Book("Z2SVMass_EESdn",Z2SVMass_EESdn, false);
    // myTree->Book("Z2SVPt_EESdn",Z2SVPt_EESdn, false);
    // myTree->Book("Z2SVEta_EESdn",Z2SVEta_EESdn, false);
    // myTree->Book("Z2SVPhi_EESdn",Z2SVPhi_EESdn, false);
    // myTree->Book("Z2GoodMass_EESdn",Z2GoodMass_EESdn, false);

    // myTree->Book("Z2Mass_MESup",Z2Mass_MESup, false);
    // myTree->Book("Z2SVMass_MESup",Z2SVMass_MESup, false);
    // myTree->Book("Z2SVPt_MESup",Z2SVPt_MESup, false);
    // myTree->Book("Z2SVEta_MESup",Z2SVEta_MESup, false);
    // myTree->Book("Z2SVPhi_MESup",Z2SVPhi_MESup, false);
    // myTree->Book("Z2GoodMass_MESup",Z2GoodMass_MESup, false);

    // myTree->Book("Z2Mass_MESdn",Z2Mass_MESdn, false);
    // myTree->Book("Z2SVMass_MESdn",Z2SVMass_MESdn, false);
    // myTree->Book("Z2SVPt_MESdn",Z2SVPt_MESdn, false);
    // myTree->Book("Z2SVEta_MESdn",Z2SVEta_MESdn, false);
    // myTree->Book("Z2SVPhi_MESdn",Z2SVPhi_MESdn, false);
    // myTree->Book("Z2GoodMass_MESdn",Z2GoodMass_MESdn, false);

    // myTree->Book("Z2Mass_JESup",Z2Mass_JESup, false);
    // myTree->Book("Z2SVMass_JESup",Z2SVMass_JESup, false);
    // myTree->Book("Z2SVPt_JESup",Z2SVPt_JESup, false);
    // myTree->Book("Z2SVEta_JESup",Z2SVEta_JESup, false);
    // myTree->Book("Z2SVPhi_JESup",Z2SVPhi_JESup, false);
    // myTree->Book("Z2GoodMass_JESup",Z2GoodMass_JESup, false);

    // myTree->Book("Z2Mass_JESdn",Z2Mass_JESdn, false);
    // myTree->Book("Z2SVMass_JESdn",Z2SVMass_JESdn, false);
    // myTree->Book("Z2SVPt_JESdn",Z2SVPt_JESdn, false);
    // myTree->Book("Z2SVEta_JESdn",Z2SVEta_JESdn, false);
    // myTree->Book("Z2SVPhi_JESdn",Z2SVPhi_JESdn, false);
    // myTree->Book("Z2GoodMass_JESdn",Z2GoodMass_JESdn, false);

    // myTree->Book("Z2Mass_JERup",Z2Mass_JERup, false);
    // myTree->Book("Z2SVMass_JERup",Z2SVMass_JERup, false);
    // myTree->Book("Z2SVPt_JERup",Z2SVPt_JERup, false);
    // myTree->Book("Z2SVEta_JERup",Z2SVEta_JERup, false);
    // myTree->Book("Z2SVPhi_JERup",Z2SVPhi_JERup, false);
    // myTree->Book("Z2GoodMass_JERup",Z2GoodMass_JERup, false);

    // myTree->Book("Z2Mass_JERdn",Z2Mass_JERdn, false);
    // myTree->Book("Z2SVMass_JERdn",Z2SVMass_JERdn, false);
    // myTree->Book("Z2SVPt_JERdn",Z2SVPt_JERdn, false);
    // myTree->Book("Z2SVEta_JERdn",Z2SVEta_JERdn, false);
    // myTree->Book("Z2SVPhi_JERdn",Z2SVPhi_JERdn, false);
    // myTree->Book("Z2GoodMass_JERdn",Z2GoodMass_JERdn, false);
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
  // myTree->Book("LepBDT",LepBDT, false);
  myTree->Book("LepisCrack",LepisCrack, false);
  myTree->Book("LepMissingHit",LepMissingHit, false);  
  myTree->Book("LepConversionVeto",LepConversionVeto, false);
  myTree->Book("LepChargedHadIso",LepChargedHadIso, false);
  myTree->Book("LepNeutralHadIso",LepNeutralHadIso, false);
  myTree->Book("LepPhotonIso",LepPhotonIso, false);
  myTree->Book("LepPUIsoComponent",LepPUIsoComponent, false);
  myTree->Book("LepCombRelIsoPF",LepCombRelIsoPF, false);

  myTree->Book("LepSF",LepSF, false);
  myTree->Book("LepSF_UncUp",LepSF_UncUp, false);
  myTree->Book("LepSF_UncDn",LepSF_UncDn, false);
  myTree->Book("LepSF_UncUp_RECO_syst",LepSF_UncUp_RECO_syst, false);
  myTree->Book("LepSF_UncUp_RECO_stat",LepSF_UncUp_RECO_stat, false);
  myTree->Book("LepSF_UncUp_ID_syst",LepSF_UncUp_ID_syst, false);
  myTree->Book("LepSF_UncUp_ID_stat",LepSF_UncUp_ID_stat, false);
  myTree->Book("LepSF_UncUp_ISO_syst",LepSF_UncUp_ISO_syst, false);
  myTree->Book("LepSF_UncUp_ISO_stat",LepSF_UncUp_ISO_stat, false);
  myTree->Book("LepSF_UncDn_RECO_syst",LepSF_UncDn_RECO_syst, false);
  myTree->Book("LepSF_UncDn_RECO_stat",LepSF_UncDn_RECO_stat, false);
  myTree->Book("LepSF_UncDn_ID_syst",LepSF_UncDn_ID_syst, false);
  myTree->Book("LepSF_UncDn_ID_stat",LepSF_UncDn_ID_stat, false);
  myTree->Book("LepSF_UncDn_ISO_syst",LepSF_UncDn_ISO_syst, false);
  myTree->Book("LepSF_UncDn_ISO_stat",LepSF_UncDn_ISO_stat, false);

  myTree->Book("LepSF_UncUp_uncert0",LepSF_UncUp_uncert0, false);
  myTree->Book("LepSF_UncUp_uncert1",LepSF_UncUp_uncert1, false);
  myTree->Book("LepSF_UncUp_syst_alleras",LepSF_UncUp_syst_alleras, false);
  myTree->Book("LepSF_UncUp_syst_year",LepSF_UncUp_syst_year, false);
  myTree->Book("LepSF_UncUp_syst_dm_year",LepSF_UncUp_syst_dm_year, false);
  myTree->Book("LepSF_UncUp_fakeEle",LepSF_UncUp_fakeEle, false);
  myTree->Book("LepSF_UncUp_fakeMu",LepSF_UncUp_fakeMu, false);
  myTree->Book("LepSF_UncDn_uncert0",LepSF_UncDn_uncert0, false);
  myTree->Book("LepSF_UncDn_uncert1",LepSF_UncDn_uncert1, false);
  myTree->Book("LepSF_UncDn_syst_alleras",LepSF_UncDn_syst_alleras, false);
  myTree->Book("LepSF_UncDn_syst_year",LepSF_UncDn_syst_year, false);
  myTree->Book("LepSF_UncDn_syst_dm_year",LepSF_UncDn_syst_dm_year, false);
  myTree->Book("LepSF_UncDn_fakeEle",LepSF_UncDn_fakeEle, false);
  myTree->Book("LepSF_UncDn_fakeMu",LepSF_UncDn_fakeMu, false);


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

  
  //Photon variables
  myTree->Book("PhotonPt",PhotonPt, failedTreeLevel >= fullFailedTree);
  myTree->Book("PhotonEta",PhotonEta, failedTreeLevel >= fullFailedTree);
  myTree->Book("PhotonPhi",PhotonPhi, failedTreeLevel >= fullFailedTree);
  myTree->Book("PhotonIsCutBasedLooseID",PhotonIsCutBasedLooseID, failedTreeLevel >= fullFailedTree);
   
  // myTree->Book("nExtraLep",nExtraLep, false);
  // myTree->Book("nExtraZ",nExtraZ, false);
  // myTree->Book("ExtraLepPt",ExtraLepPt, false);
  // myTree->Book("ExtraLepEta",ExtraLepEta, false);
  // myTree->Book("ExtraLepPhi",ExtraLepPhi, false);
  // myTree->Book("ExtraLepLepId",ExtraLepLepId, false);

  // myTree->Book("ZXFakeweight", ZXFakeweight, false);

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
    // myTree->Book("genProcessId", genProcessId, failedTreeLevel >= minimalFailedTree);
    myTree->Book("genHEPMCweight", genHEPMCweight, failedTreeLevel >= minimalFailedTree);
    myTree->Book("PUWeight", PUWeight, failedTreeLevel >= minimalFailedTree);
    myTree->Book("PUWeight_Dn", PUWeight_Dn, failedTreeLevel >= minimalFailedTree);
    myTree->Book("PUWeight_Up", PUWeight_Up, failedTreeLevel >= minimalFailedTree);
    myTree->Book("dataMCWeight", dataMCWeight, false);
    myTree->Book("trigEffWeight", trigEffWeight, false);
    myTree->Book("overallEventWeight", overallEventWeight, false);
    myTree->Book("L1prefiringWeight", L1prefiringWeight, false);
    myTree->Book("L1prefiringWeightUp", L1prefiringWeightUp, false);
    myTree->Book("L1prefiringWeightDn", L1prefiringWeightDn, false);
    myTree->Book("L1prefiringWeight_ECAL", L1prefiringWeight_ECAL, false);
    myTree->Book("L1prefiringWeightUp_ECAL", L1prefiringWeightUp_ECAL, false);
    myTree->Book("L1prefiringWeightDn_ECAL", L1prefiringWeightDn_ECAL, false);
    myTree->Book("L1prefiringWeight_Mu", L1prefiringWeight_Mu, false);
    myTree->Book("L1prefiringWeightUp_Mu", L1prefiringWeightUp_Mu, false);
    myTree->Book("L1prefiringWeightDn_Mu", L1prefiringWeightDn_Mu, false);

    myTree->Book("xsec", xsection, failedTreeLevel >= minimalFailedTree);
    myTree->Book("genxsec", genxsection, failedTreeLevel >= minimalFailedTree);
    myTree->Book("genBR", genbranchingratio, failedTreeLevel >= minimalFailedTree);
    myTree->Book("genExtInfo", genExtInfo, failedTreeLevel >= minimalFailedTree);
    myTree->Book("GenLLMass", GenLLMass, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLLPt", GenLLPt, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLLPhi", GenLLPhi, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLLEta", GenLLEta, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLLFlav", GenLLFlav, failedTreeLevel >= minimalFailedTree);
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
    myTree->Book("GenVisLLMass", GenVisLLMass, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisLLPt", GenVisLLPt, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisLLEta", GenVisLLEta, failedTreeLevel >= minimalFailedTree);
    myTree->Book("GenVisLLPhi", GenVisLLPhi, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisLLFlav", GenVisLLFlav, failedTreeLevel >= minimalFailedTree);
    myTree->Book("GenVisLep1Pt", GenVisLep1Pt, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisLep1Eta", GenVisLep1Eta, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisLep1Phi", GenVisLep1Phi, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisLep1M", GenVisLep1M, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisLep1Id", GenVisLep1Id, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisLep2Pt", GenVisLep2Pt, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisLep2Eta", GenVisLep2Eta, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisLep2Phi", GenVisLep2Phi, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenVisLep2M", GenVisLep2M, failedTreeLevel >= fullFailedTree);
    // myTree->Book("GenVisLep2Id", GenVisLep2Id, failedTreeLevel >= fullFailedTree);
    // myTree->Book("GenAssocLep1Pt", GenAssocLep1Pt, failedTreeLevel >= fullFailedTree);
    // myTree->Book("GenAssocLep1Eta", GenAssocLep1Eta, failedTreeLevel >= fullFailedTree);
    // myTree->Book("GenAssocLep1Phi", GenAssocLep1Phi, failedTreeLevel >= fullFailedTree);
    // myTree->Book("GenAssocLep1Id", GenAssocLep1Id, failedTreeLevel >= fullFailedTree);
    // myTree->Book("GenAssocLep2Pt", GenAssocLep2Pt, failedTreeLevel >= fullFailedTree);
    // myTree->Book("GenAssocLep2Eta", GenAssocLep2Eta, failedTreeLevel >= fullFailedTree);
    // myTree->Book("GenAssocLep2Phi", GenAssocLep2Phi, failedTreeLevel >= fullFailedTree);
    // myTree->Book("GenAssocLep2Id", GenAssocLep2Id, failedTreeLevel >= fullFailedTree);

    myTree->Book("GenLep1Iso", GenLep1Iso, failedTreeLevel >= minimalFailedTree); //AT
    myTree->Book("GenLep2Iso", GenLep2Iso, failedTreeLevel >= minimalFailedTree); //AT
    
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

void LLNtupleMaker::addweight(float &weight, float weighttoadd) {
  weight += weighttoadd;
}


Int_t LLNtupleMaker::FindBinValue(TGraphErrors *tgraph, double value)
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
DEFINE_FWK_MODULE(LLNtupleMaker);
