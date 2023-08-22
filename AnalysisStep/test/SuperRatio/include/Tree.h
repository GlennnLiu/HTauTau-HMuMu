//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jul 21 14:52:28 2016 by ROOT version 6.02/05
// from TTree candTree/Event Summary
// found on file: ICHEP_2016/ggH125/ZZ4lAnalysis.root
//////////////////////////////////////////////////////////

#ifndef Tree_h
#define Tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

using namespace std;

class Tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           RunNumber;
   Long64_t        EventNumber;
   Int_t           LumiNumber;
   Short_t         Nvtx;
   Short_t         NObsInt;
   Float_t         NTrueInt;

   Float_t         PFMET;
   Float_t         PFMETPhi;
   Float_t         PFMETRecoil;
   Float_t         PFMETPhiRecoil;

   Short_t          pass_Trigger;
   Short_t          pass_SingleTrigger;
   Short_t          pass_CrossTrigger;
//   Float_t         PFMETNoHF;
//   Float_t         PFMETNoHFPhi;

   Short_t         nCleanedJets;
   Short_t         nCleanedJetsPt30;
   Short_t         nCleanedJetsPt30BTagged;
   Short_t         nCleanedJetsPt25BTagged_bTagSF;
   Short_t         nCleanedJetsPt30BTagged_bTagSF;
   Short_t         trigWord;
   Float_t         LLMass;
   Float_t         LLGoodMass;
   Float_t         LLPt;
   Float_t         LLEta;
   Float_t         LLPhi;
   Short_t         LLFlav;
   Float_t         LLDR;
   Float_t         LLSVPt;


   Float_t         MtLMET;
   Float_t         Pzeta1;
   Float_t         Pzeta2;

   
   Float_t         DeltaEtaJJ;
   Float_t         DiJetMass;
   Short_t         VBFJetIdx1;
   Short_t         VBFJetIdx2;

   vector<float>   *LepPt;
   vector<float>   *LepEta;
   vector<float>   *LepPhi;
   vector<short>   *LepLepId;
   vector<float>   *LepSIP;
   vector<float>   *Lepdxy;
   vector<float>   *Lepdz;
   vector<float>   *LepTime;
   vector<bool>    *LepisID;
   // vector<short>   *LepisLoose;
   // vector<float>   *LepBDT;
   vector<char>    *LepMissingHit;
   vector<float>   *LepCombRelIsoPF;
//tau
   vector<short>   *TauVSmu;
   vector<short>   *TauVSe;
   vector<short>   *TauVSjet;
   vector<float>   *TauDecayMode;
   vector<short>   *TauGenMatch;

   vector<float>   *fsrPt;
   vector<float>   *fsrEta;
   vector<float>   *fsrPhi;
   vector<short>   *fsrLept;
   Bool_t          passIsoPreFSR;
   
   vector<float>   *JetPt;
   vector<float>   *JetEta;
   vector<float>   *JetPhi;
   vector<float>   *JetMass;
   vector<float>   *JetBTagger;
   vector<float>   *JetIsBtagged;
   vector<float>   *JetIsBtaggedWithSF;
   vector<float>   *JetIsBtaggedWithSFUp;
   vector<float>   *JetIsBtaggedWithSFDn;
   vector<float>   *JetQGLikelihood;
   vector<float>   *JetAxis2;
   vector<float>   *JetMult;
   vector<float>   *JetPtD;
   vector<float>   *JetSigma;

   Float_t         ggH_NNLOPS_weight;
   Float_t         KFactor_QCD_ggZZ_Nominal;
   Float_t         KFactor_QCD_ggZZ_PDFScaleDn;
   Float_t         KFactor_QCD_ggZZ_PDFScaleUp;
   Float_t         KFactor_QCD_ggZZ_QCDScaleDn;
   Float_t         KFactor_QCD_ggZZ_QCDScaleUp;
   Float_t         KFactor_QCD_ggZZ_AsDn;
   Float_t         KFactor_QCD_ggZZ_AsUp;
   Float_t         KFactor_QCD_ggZZ_PDFReplicaDn;
   Float_t         KFactor_QCD_ggZZ_PDFReplicaUp;
   Float_t         KFactor_EW_qqZZ;
   Float_t         KFactor_EW_qqZZ_unc;
   Float_t         KFactor_QCD_qqZZ_dPhi;
   Float_t         KFactor_QCD_qqZZ_M;
   Float_t         KFactor_QCD_qqZZ_Pt;
   Short_t         genFinalState;
   // Int_t           genProcessId;
   Float_t         genHEPMCweight;
   Float_t         PUWeight;
   Float_t         dataMCWeight;
   Float_t         overallEventWeight;
   Float_t         L1prefiringWeight;
   Float_t         xsec;
   Short_t         genExtInfo;
   Float_t         GenLLMass;
   Float_t         GenLLPt;
   Float_t         GenLLEta;
   Float_t         GenLLPhi;
   Float_t         GenLLFlav;
   Float_t         GenLep1Pt;
   Float_t         GenLep1Eta;
   Float_t         GenLep1Phi;
   Short_t         GenLep1Id;
   Float_t         GenLep2Pt;
   Float_t         GenLep2Eta;
   Float_t         GenLep2Phi;
   Short_t         GenLep2Id;
   //Float_t         LHEweight_QCDscale_muR1_muF1;
   //Float_t         LHEweight_QCDscale_muR1_muF2;
   //Float_t         LHEweight_QCDscale_muR1_muF0p5;
   //Float_t         LHEweight_QCDscale_muR2_muF1;
   //Float_t         LHEweight_QCDscale_muR2_muF2;
   //Float_t         LHEweight_QCDscale_muR2_muF0p5;
   //Float_t         LHEweight_QCDscale_muR0p5_muF1;
   //Float_t         LHEweight_QCDscale_muR0p5_muF2;
   //Float_t         LHEweight_QCDscale_muR0p5_muF0p5;

   // List of branches

   TBranch        *b_RunNumber;
   TBranch        *b_EventNumber;
   TBranch        *b_LumiNumber;
   TBranch        *b_Nvtx;
   TBranch        *b_NObsInt;
   TBranch        *b_NTrueInt;
   TBranch        *b_PFMET;
   TBranch        *b_PFMETPhi;
   TBranch        *b_PFMETRecoil;
   TBranch        *b_PFMETPhiRecoil;
   TBranch        *b_pass_Trigger;
   TBranch        *b_pass_SingleTrigger;
   TBranch        *b_pass_CrossTrigger;
//   TBranch        *b_PFMETNoHF;
//   TBranch        *b_PFMETNoHFPhi;
   TBranch        *b_nCleanedJets;
   TBranch        *b_nCleanedJetsPt30;
   TBranch        *b_nCleanedJetsPt30BTagged;
   TBranch        *b_nCleanedJetsPt25BTagged_bTagSF;
   TBranch        *b_nCleanedJetsPt30BTagged_bTagSF;
   TBranch        *b_trigWord;
   TBranch        *b_LLMass;
   TBranch        *b_LLGoodMass;
   TBranch        *b_LLPt;
   TBranch        *b_LLEta;
   TBranch        *b_LLPhi;
   TBranch        *b_LLFlav;
   TBranch        *b_LLDR;
   TBranch        *b_LLSVPt;

   TBranch        *b_MtLMET;
   TBranch        *b_Pzeta1;
   TBranch        *b_Pzeta2;

   TBranch        *b_DeltaEtaJJ;
   TBranch        *b_DiJetMass;
   TBranch        *b_VBFJetIdx1;
   TBranch        *b_VBFJetIdx2;

   TBranch        *b_LepPt;
   TBranch        *b_LepEta;
   TBranch        *b_LepPhi;
   TBranch        *b_LepLepId;
   TBranch        *b_LepSIP;
   TBranch        *b_Lepdxy;
   TBranch        *b_Lepdz;
   TBranch        *b_LepTime;
   TBranch        *b_LepisID;
   TBranch        *b_LepisLoose;
   TBranch        *b_LepBDT;
   TBranch        *b_LepMissingHit;
   TBranch        *b_LepCombRelIsoPF;
//tau
   TBranch        *b_TauVSmu;
   TBranch        *b_TauVSe;
   TBranch        *b_TauVSjet;
   TBranch        *b_TauDecayMode;
   TBranch        *b_TauGenMatch;

   TBranch        *b_fsrPt;
   TBranch        *b_fsrEta;
   TBranch        *b_fsrPhi;
   TBranch        *b_fsrLept;
   TBranch        *b_passIsoPreFSR;
   
   TBranch        *b_JetPt;
   TBranch        *b_JetEta;
   TBranch        *b_JetPhi;
   TBranch        *b_JetMass;
   TBranch        *b_JetBTagger;
   TBranch        *b_JetIsBtagged;
   TBranch        *b_JetIsBtaggedWithSF;
   TBranch        *b_JetIsBtaggedWithSFUp;
   TBranch        *b_JetIsBtaggedWithSFDn;
   TBranch        *b_JetQGLikelihood;
   TBranch        *b_JetAxis2;
   TBranch        *b_JetMult;
   TBranch        *b_JetPtD;
   TBranch        *b_JetSigma;

   TBranch        *b_ggH_NNLOPS_weight;
   TBranch        *b_KFactor_QCD_ggZZ_Nominal;
   TBranch        *b_KFactor_QCD_ggZZ_PDFScaleDn;
   TBranch        *b_KFactor_QCD_ggZZ_PDFScaleUp;
   TBranch        *b_KFactor_QCD_ggZZ_QCDScaleDn;
   TBranch        *b_KFactor_QCD_ggZZ_QCDScaleUp;
   TBranch        *b_KFactor_QCD_ggZZ_AsDn;
   TBranch        *b_KFactor_QCD_ggZZ_AsUp;
   TBranch        *b_KFactor_QCD_ggZZ_PDFReplicaDn;
   TBranch        *b_KFactor_QCD_ggZZ_PDFReplicaUp;
   TBranch        *b_KFactor_EW_qqZZ;
   TBranch        *b_KFactor_EW_qqZZ_unc;
   TBranch        *b_KFactor_QCD_qqZZ_dPhi;
   TBranch        *b_KFactor_QCD_qqZZ_M;
   TBranch        *b_KFactor_QCD_qqZZ_Pt;
   TBranch        *b_genFinalState;
   // TBranch        *b_genProcessId;
   TBranch        *b_genHEPMCweight;
   TBranch        *b_PUWeight;
   TBranch        *b_dataMCWeight;
   TBranch        *b_overallEventWeight;
   TBranch        *b_L1prefiringWeight;
   TBranch        *b_xsec;
   TBranch        *b_genExtInfo;
   TBranch        *b_GenLLMass;
   TBranch        *b_GenLLPt;
   TBranch        *b_GenLLEta;
   TBranch        *b_GenLLPhi;
   TBranch        *b_GenLLFlav;
   TBranch        *b_GenLep1Pt;
   TBranch        *b_GenLep1Eta;
   TBranch        *b_GenLep1Phi;
   TBranch        *b_GenLep1Id;
   TBranch        *b_GenLep2Pt;
   TBranch        *b_GenLep2Eta;
   TBranch        *b_GenLep2Phi;
   TBranch        *b_GenLep2Id;
   //TBranch        *b_LHEweight_QCDscale_muR1_muF1;   //!
   //TBranch        *b_LHEweight_QCDscale_muR1_muF2;   //!
   //TBranch        *b_LHEweight_QCDscale_muR1_muF0p5;   //!
   //TBranch        *b_LHEweight_QCDscale_muR2_muF1;   //!
   //TBranch        *b_LHEweight_QCDscale_muR2_muF2;   //!
   //TBranch        *b_LHEweight_QCDscale_muR2_muF0p5;   //!
   //TBranch        *b_LHEweight_QCDscale_muR0p5_muF1;   //!
   //TBranch        *b_LHEweight_QCDscale_muR0p5_muF2;   //!
   //TBranch        *b_LHEweight_QCDscale_muR0p5_muF0p5;   //!

   Tree(TTree *tree=0);
   virtual ~Tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree, TString input_file_name, bool notZLregion);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Tree_cxx
Tree::Tree(TTree *tree) : fChain(0) 
{
}

Tree::~Tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Tree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Tree::Init(TTree *tree, TString input_file_name, bool notZLregion)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   LepPt = 0;
   LepEta = 0;
   LepPhi = 0;
   LepLepId = 0;
   LepSIP = 0;
   Lepdxy = 0;
   Lepdz = 0;
   LepTime = 0;
   LepisID = 0;
   // LepisLoose = 0;
   // LepBDT = 0;
   LepMissingHit = 0;
   LepCombRelIsoPF = 0;
//tau
   TauVSmu = 0;
   TauVSe = 0;
   TauVSjet = 0;
   TauDecayMode = 0;
   TauGenMatch = 0;

   fsrPt = 0;
   fsrEta = 0;
   fsrPhi = 0;
   fsrLept = 0;
   
   JetPt = 0;
   JetEta = 0;
   JetPhi = 0;
   JetMass = 0;
   JetBTagger = 0;
   JetIsBtagged = 0;
   JetIsBtaggedWithSF = 0;
   JetIsBtaggedWithSFUp = 0;
   JetIsBtaggedWithSFDn = 0;
   JetQGLikelihood = 0;
   JetAxis2 = 0;
   JetMult = 0;
   JetPtD = 0;
   JetSigma = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("LumiNumber", &LumiNumber, &b_LumiNumber);
   fChain->SetBranchAddress("Nvtx", &Nvtx, &b_Nvtx);
   fChain->SetBranchAddress("NObsInt", &NObsInt, &b_NObsInt);
   fChain->SetBranchAddress("NTrueInt", &NTrueInt, &b_NTrueInt);

   fChain->SetBranchAddress("PFMET", &PFMET, &b_PFMET);
   fChain->SetBranchAddress("PFMETPhi", &PFMETPhi, &b_PFMETPhi);
   fChain->SetBranchAddress("PFMETRecoil", &PFMETRecoil, &b_PFMETRecoil);
   fChain->SetBranchAddress("PFMETPhiRecoil", &PFMETPhiRecoil, &b_PFMETPhiRecoil);
   fChain->SetBranchAddress("pass_Trigger", &pass_Trigger, &b_pass_Trigger);
   fChain->SetBranchAddress("pass_SingleTrigger", &pass_SingleTrigger, &b_pass_SingleTrigger);
   fChain->SetBranchAddress("pass_CrossTrigger", &pass_CrossTrigger, &b_pass_CrossTrigger);

   fChain->SetBranchAddress("nCleanedJets", &nCleanedJets, &b_nCleanedJets);
   fChain->SetBranchAddress("nCleanedJetsPt30", &nCleanedJetsPt30, &b_nCleanedJetsPt30);
   fChain->SetBranchAddress("nCleanedJetsPt30BTagged", &nCleanedJetsPt30BTagged, &b_nCleanedJetsPt30BTagged);
   fChain->SetBranchAddress("nCleanedJetsPt25BTagged_bTagSF", &nCleanedJetsPt25BTagged_bTagSF, &b_nCleanedJetsPt25BTagged_bTagSF);
   fChain->SetBranchAddress("nCleanedJetsPt30BTagged_bTagSF", &nCleanedJetsPt30BTagged_bTagSF, &b_nCleanedJetsPt30BTagged_bTagSF);
   fChain->SetBranchAddress("trigWord", &trigWord, &b_trigWord);

   fChain->SetBranchAddress("LLMass", &LLMass, &b_LLMass);
   fChain->SetBranchAddress("LLGoodMass", &LLGoodMass, &b_LLGoodMass);
   fChain->SetBranchAddress("LLPt", &LLPt, &b_LLPt);
   fChain->SetBranchAddress("LLEta", &LLEta, &b_LLEta);
   fChain->SetBranchAddress("LLPhi", &LLPhi, &b_LLPhi);
   fChain->SetBranchAddress("LLFlav", &LLFlav, &b_LLFlav);
   fChain->SetBranchAddress("LLDR", &LLDR, &b_LLDR);
   fChain->SetBranchAddress("LLSVPt", &LLSVPt, &b_LLSVPt);

   fChain->SetBranchAddress("MtLMET", &MtLMET, &b_MtLMET);
   fChain->SetBranchAddress("Pzeta1", &Pzeta1, &b_Pzeta1);
   fChain->SetBranchAddress("Pzeta2", &Pzeta2, &b_Pzeta2);

   fChain->SetBranchAddress("DeltaEtaJJ", &DeltaEtaJJ, &b_DeltaEtaJJ);
   fChain->SetBranchAddress("DiJetMass", &DiJetMass, &b_DiJetMass);
   fChain->SetBranchAddress("VBFJetIdx1", &VBFJetIdx1, &b_VBFJetIdx1);
   fChain->SetBranchAddress("VBFJetIdx2", &VBFJetIdx2, &b_VBFJetIdx2);

   fChain->SetBranchAddress("LepPt", &LepPt, &b_LepPt);
   fChain->SetBranchAddress("LepEta", &LepEta, &b_LepEta);
   fChain->SetBranchAddress("LepPhi", &LepPhi, &b_LepPhi);
   fChain->SetBranchAddress("LepLepId", &LepLepId, &b_LepLepId);
   fChain->SetBranchAddress("LepSIP", &LepSIP, &b_LepSIP);
   fChain->SetBranchAddress("Lepdxy", &Lepdxy, &b_Lepdxy);
   fChain->SetBranchAddress("Lepdz", &Lepdz, &b_Lepdz);
   fChain->SetBranchAddress("LepTime", &LepTime, &b_LepTime);
   fChain->SetBranchAddress("LepisID", &LepisID, &b_LepisID);
   // fChain->SetBranchAddress("LepisLoose", &LepisLoose, &b_LepisLoose);
   // fChain->SetBranchAddress("LepBDT", &LepBDT, &b_LepBDT);
   fChain->SetBranchAddress("LepMissingHit", &LepMissingHit, &b_LepMissingHit);
   fChain->SetBranchAddress("LepCombRelIsoPF", &LepCombRelIsoPF, &b_LepCombRelIsoPF);
//tau
   fChain->SetBranchAddress("TauVSmu", &TauVSmu, &b_TauVSmu);
   fChain->SetBranchAddress("TauVSe", &TauVSe, &b_TauVSe);
   fChain->SetBranchAddress("TauVSjet", &TauVSjet, &b_TauVSjet);
   fChain->SetBranchAddress("TauDecayMode", &TauDecayMode, &b_TauDecayMode);
   fChain->SetBranchAddress("TauGenMatch", &TauGenMatch, &b_TauGenMatch);

   fChain->SetBranchAddress("fsrPt", &fsrPt, &b_fsrPt);
   fChain->SetBranchAddress("fsrEta", &fsrEta, &b_fsrEta);
   fChain->SetBranchAddress("fsrPhi", &fsrPhi, &b_fsrPhi);
   fChain->SetBranchAddress("fsrLept", &fsrLept, &b_fsrLept);
   fChain->SetBranchAddress("passIsoPreFSR", &passIsoPreFSR, &b_passIsoPreFSR);
   
   fChain->SetBranchAddress("JetPt", &JetPt, &b_JetPt);
   fChain->SetBranchAddress("JetEta", &JetEta, &b_JetEta);
   fChain->SetBranchAddress("JetPhi", &JetPhi, &b_JetPhi);
   fChain->SetBranchAddress("JetMass", &JetMass, &b_JetMass);
   fChain->SetBranchAddress("JetBTagger", &JetBTagger, &b_JetBTagger);
   fChain->SetBranchAddress("JetIsBtagged", &JetIsBtagged, &b_JetIsBtagged);
   fChain->SetBranchAddress("JetIsBtaggedWithSF", &JetIsBtaggedWithSF, &b_JetIsBtaggedWithSF);
   fChain->SetBranchAddress("JetIsBtaggedWithSFUp", &JetIsBtaggedWithSFUp, &b_JetIsBtaggedWithSFUp);
   fChain->SetBranchAddress("JetIsBtaggedWithSFDn", &JetIsBtaggedWithSFDn, &b_JetIsBtaggedWithSFDn);
   fChain->SetBranchAddress("JetQGLikelihood", &JetQGLikelihood, &b_JetQGLikelihood);
   fChain->SetBranchAddress("JetAxis2", &JetAxis2, &b_JetAxis2);
   fChain->SetBranchAddress("JetMult", &JetMult, &b_JetMult);
   fChain->SetBranchAddress("JetPtD", &JetPtD, &b_JetPtD);
   fChain->SetBranchAddress("JetSigma", &JetSigma, &b_JetSigma);
   
   if ( !(input_file_name.Contains("Data")) )
   {
      if ( input_file_name.Contains("ggH") )
      {
         fChain->SetBranchAddress("ggH_NNLOPS_weight", &ggH_NNLOPS_weight, &b_ggH_NNLOPS_weight);
      }
      if ( input_file_name.Contains("ggT") )
      {
         fChain->SetBranchAddress("KFactor_QCD_ggZZ_Nominal", &KFactor_QCD_ggZZ_Nominal, &b_KFactor_QCD_ggZZ_Nominal);
         fChain->SetBranchAddress("KFactor_QCD_ggZZ_PDFScaleDn", &KFactor_QCD_ggZZ_PDFScaleDn, &b_KFactor_QCD_ggZZ_PDFScaleDn);
         fChain->SetBranchAddress("KFactor_QCD_ggZZ_PDFScaleUp", &KFactor_QCD_ggZZ_PDFScaleUp, &b_KFactor_QCD_ggZZ_PDFScaleUp);
         fChain->SetBranchAddress("KFactor_QCD_ggZZ_QCDScaleDn", &KFactor_QCD_ggZZ_QCDScaleDn, &b_KFactor_QCD_ggZZ_QCDScaleDn);
         fChain->SetBranchAddress("KFactor_QCD_ggZZ_QCDScaleUp", &KFactor_QCD_ggZZ_QCDScaleUp, &b_KFactor_QCD_ggZZ_QCDScaleUp);
         fChain->SetBranchAddress("KFactor_QCD_ggZZ_AsDn", &KFactor_QCD_ggZZ_AsDn, &b_KFactor_QCD_ggZZ_AsDn);
         fChain->SetBranchAddress("KFactor_QCD_ggZZ_AsUp", &KFactor_QCD_ggZZ_AsUp, &b_KFactor_QCD_ggZZ_AsUp);
         fChain->SetBranchAddress("KFactor_QCD_ggZZ_PDFReplicaDn", &KFactor_QCD_ggZZ_PDFReplicaDn, &b_KFactor_QCD_ggZZ_PDFReplicaDn);
         fChain->SetBranchAddress("KFactor_QCD_ggZZ_PDFReplicaUp", &KFactor_QCD_ggZZ_PDFReplicaUp, &b_KFactor_QCD_ggZZ_PDFReplicaUp);
      }
      
      if ( input_file_name.Contains("ZZTo") )
      {
         fChain->SetBranchAddress("KFactor_EW_qqZZ", &KFactor_EW_qqZZ, &b_KFactor_EW_qqZZ);
         fChain->SetBranchAddress("KFactor_EW_qqZZ_unc", &KFactor_EW_qqZZ_unc, &b_KFactor_EW_qqZZ_unc);
         fChain->SetBranchAddress("KFactor_QCD_qqZZ_dPhi", &KFactor_QCD_qqZZ_dPhi, &b_KFactor_QCD_qqZZ_dPhi);
         fChain->SetBranchAddress("KFactor_QCD_qqZZ_M", &KFactor_QCD_qqZZ_M, &b_KFactor_QCD_qqZZ_M);
         fChain->SetBranchAddress("KFactor_QCD_qqZZ_Pt", &KFactor_QCD_qqZZ_Pt, &b_KFactor_QCD_qqZZ_Pt);
      }
      
      fChain->SetBranchAddress("genFinalState", &genFinalState, &b_genFinalState);
      // fChain->SetBranchAddress("genProcessId", &genProcessId, &b_genProcessId);
      fChain->SetBranchAddress("genHEPMCweight", &genHEPMCweight, &b_genHEPMCweight);
      fChain->SetBranchAddress("PUWeight", &PUWeight, &b_PUWeight);
      fChain->SetBranchAddress("dataMCWeight", &dataMCWeight, &b_dataMCWeight);
      fChain->SetBranchAddress("overallEventWeight", &overallEventWeight, &b_overallEventWeight);
      fChain->SetBranchAddress("L1prefiringWeight", &L1prefiringWeight, &b_L1prefiringWeight);
      fChain->SetBranchAddress("xsec", &xsec, &b_xsec);
      fChain->SetBranchAddress("genExtInfo", &genExtInfo, &b_genExtInfo);
         
      fChain->SetBranchAddress("GenLLMass", &GenLLMass, &b_GenLLMass);
      fChain->SetBranchAddress("GenLLPt", &GenLLPt, &b_GenLLPt);
      fChain->SetBranchAddress("GenLLEta", &GenLLEta, &b_GenLLEta);
      fChain->SetBranchAddress("GenLLPhi", &GenLLPhi, &b_GenLLPhi);
      fChain->SetBranchAddress("GenLLFlav", &GenLLFlav, &b_GenLLFlav);
      
      fChain->SetBranchAddress("GenLep1Pt", &GenLep1Pt, &b_GenLep1Pt);
      fChain->SetBranchAddress("GenLep1Eta", &GenLep1Eta, &b_GenLep1Eta);
      fChain->SetBranchAddress("GenLep1Phi", &GenLep1Phi, &b_GenLep1Phi);
      fChain->SetBranchAddress("GenLep1Id", &GenLep1Id, &b_GenLep1Id);
      fChain->SetBranchAddress("GenLep2Pt", &GenLep2Pt, &b_GenLep2Pt);
      fChain->SetBranchAddress("GenLep2Eta", &GenLep2Eta, &b_GenLep2Eta);
      fChain->SetBranchAddress("GenLep2Phi", &GenLep2Phi, &b_GenLep2Phi);
      fChain->SetBranchAddress("GenLep2Id", &GenLep2Id, &b_GenLep2Id);
      //   fChain->SetBranchAddress("reweightingweights", &reweightingweights, &b_reweightingweights);
      //fChain->SetBranchAddress("LHEPDFScale", &LHEPDFScale, &b_LHEPDFScale);
      //fChain->SetBranchAddress("LHEweight_QCDscale_muR1_muF1", &LHEweight_QCDscale_muR1_muF1, &b_LHEweight_QCDscale_muR1_muF1);
      //fChain->SetBranchAddress("LHEweight_QCDscale_muR1_muF2", &LHEweight_QCDscale_muR1_muF2, &b_LHEweight_QCDscale_muR1_muF2);
      //fChain->SetBranchAddress("LHEweight_QCDscale_muR1_muF0p5", &LHEweight_QCDscale_muR1_muF0p5, &b_LHEweight_QCDscale_muR1_muF0p5);
      //fChain->SetBranchAddress("LHEweight_QCDscale_muR2_muF1", &LHEweight_QCDscale_muR2_muF1, &b_LHEweight_QCDscale_muR2_muF1);
      //fChain->SetBranchAddress("LHEweight_QCDscale_muR2_muF2", &LHEweight_QCDscale_muR2_muF2, &b_LHEweight_QCDscale_muR2_muF2);
      //fChain->SetBranchAddress("LHEweight_QCDscale_muR2_muF0p5", &LHEweight_QCDscale_muR2_muF0p5, &b_LHEweight_QCDscale_muR2_muF0p5);
      //fChain->SetBranchAddress("LHEweight_QCDscale_muR0p5_muF1", &LHEweight_QCDscale_muR0p5_muF1, &b_LHEweight_QCDscale_muR0p5_muF1);
      //fChain->SetBranchAddress("LHEweight_QCDscale_muR0p5_muF2", &LHEweight_QCDscale_muR0p5_muF2, &b_LHEweight_QCDscale_muR0p5_muF2);
      //fChain->SetBranchAddress("LHEweight_QCDscale_muR0p5_muF0p5", &LHEweight_QCDscale_muR0p5_muF0p5, &b_LHEweight_QCDscale_muR0p5_muF0p5);
   }
   
   Notify();
}

Bool_t Tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Tree_cxx
