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
    Int_t RunNumber;
    Long64_t EventNumber;
    Int_t LumiNumber;

    Float_t GenMET;
    Float_t GenMETPhi;
    Float_t PFMET;
    Float_t PFMETPhi;
    Float_t METx;
    Float_t METy;
    Float_t METxUPTES;
    Float_t METyUPTES;
    Float_t METxDOWNTES;
    Float_t METyDOWNTES;
    Float_t METxUPEES;
    Float_t METyUPEES;
    Float_t METxDOWNEES;
    Float_t METyDOWNEES;
    Float_t METxUPMES;
    Float_t METyUPMES;
    Float_t METxDOWNMES;
    Float_t METyDOWNMES;
    Float_t METxUPJES;
    Float_t METyUPJES;
    Float_t METxDOWNJES;
    Float_t METyDOWNJES;
    Float_t METxUPJER;
    Float_t METyUPJER;
    Float_t METxDOWNJER;
    Float_t METyDOWNJER;
    
    Short_t nCleanedJets;
    Float_t ZZMass;
    Short_t ZZsel;
    Float_t ZZPt;
    Float_t ZZjjPt;
    Float_t ZZEta;
    Float_t ZZPhi;
    Int_t   CRflag;
    
    Float_t Z1Mass;
    Float_t Z1Pt;
    Short_t Z1Flav;
    
    Float_t Z2Mass;
    Float_t Z2Pt;
    Short_t Z2Flav;
    Short_t Z2isGoodTau;
    
    Float_t Z2SVMass;
    Float_t Z2SVPt;
    Float_t Z2SVEta;
    Float_t Z2SVPhi;
    Float_t Z2GoodMass;

    Float_t Z2Mass_TESup;
    Float_t Z2SVMass_TESup;
    Float_t Z2SVPt_TESup;
    Float_t Z2SVEta_TESup;
    Float_t Z2SVPhi_TESup;
    Float_t Z2GoodMass_TESup;

    Float_t Z2Mass_EESup;
    Float_t Z2SVMass_EESup;
    Float_t Z2SVPt_EESup;
    Float_t Z2SVEta_EESup;
    Float_t Z2SVPhi_EESup;
    Float_t Z2GoodMass_EESup;

    Float_t Z2Mass_MESup;
    Float_t Z2SVMass_MESup;
    Float_t Z2SVPt_MESup;
    Float_t Z2SVEta_MESup;
    Float_t Z2SVPhi_MESup;
    Float_t Z2GoodMass_MESup;

    Float_t Z2Mass_JESup;
    Float_t Z2SVMass_JESup;
    Float_t Z2SVPt_JESup;
    Float_t Z2SVEta_JESup;
    Float_t Z2SVPhi_JESup;
    Float_t Z2GoodMass_JESup;

    Float_t Z2Mass_JERup;
    Float_t Z2SVMass_JERup;
    Float_t Z2SVPt_JERup;
    Float_t Z2SVEta_JERup;
    Float_t Z2SVPhi_JERup;
    Float_t Z2GoodMass_JERup;

    Float_t Z2Mass_TESdn;
    Float_t Z2SVMass_TESdn;
    Float_t Z2SVPt_TESdn;
    Float_t Z2SVEta_TESdn;
    Float_t Z2SVPhi_TESdn;
    Float_t Z2GoodMass_TESdn;

    Float_t Z2Mass_EESdn;
    Float_t Z2SVMass_EESdn;
    Float_t Z2SVPt_EESdn;
    Float_t Z2SVEta_EESdn;
    Float_t Z2SVPhi_EESdn;
    Float_t Z2GoodMass_EESdn;

    Float_t Z2Mass_MESdn;
    Float_t Z2SVMass_MESdn;
    Float_t Z2SVPt_MESdn;
    Float_t Z2SVEta_MESdn;
    Float_t Z2SVPhi_MESdn;
    Float_t Z2GoodMass_MESdn;

    Float_t Z2Mass_JESdn;
    Float_t Z2SVMass_JESdn;
    Float_t Z2SVPt_JESdn;
    Float_t Z2SVEta_JESdn;
    Float_t Z2SVPhi_JESdn;
    Float_t Z2GoodMass_JESdn;

    Float_t Z2Mass_JERdn;
    Float_t Z2SVMass_JERdn;
    Float_t Z2SVPt_JERdn;
    Float_t Z2SVEta_JERdn;
    Float_t Z2SVPhi_JERdn;
    Float_t Z2GoodMass_JERdn;
    
    Float_t costhetastar;
    Float_t helphi;
    Float_t helcosthetaZ1;
    Float_t helcosthetaZ2;
    Float_t phistarZ1;
    Float_t phistarZ2;
    Float_t xi;
    Float_t xistar;
    
    vector<float>   *LepPt;
    vector<float>   *LepEta;
    vector<float>   *LepPhi;
    vector<short>   *LepLepId;
    vector<float>   *LepSIP;
    vector<float>   *Lepdxy;
    vector<float>   *Lepdz;
    vector<float>   *LepTime;
    vector<bool>    *LepisID;
    vector<short>   *LepisLoose;
    vector<float>   *LepBDT;
    vector<char>    *LepMissingHit;
    vector<float>   *LepCombRelIsoPF;
    vector<float>   *LepSF;
    vector<float>   *LepSF_Unc;
    //tau
    vector<short>   *TauVSmu;
    vector<short>   *TauVSe;
    vector<short>   *TauVSjet;
    vector<float>   *TauDecayMode;
    vector<short>   *TauGenMatch;
    vector<float>   *TauTES_p_Up;
    vector<float>   *TauFES_p_Up;

    vector<float>   *fsrPt;
    vector<float>   *fsrEta;
    vector<float>   *fsrPhi;
    vector<short>   *fsrLept;
    Bool_t  passIsoPreFSR;

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
    Float_t DiJetMass;
    Float_t DiJetDEta;
    Float_t DiJetFisher;
    Short_t nExtraLep;
    Short_t nExtraZ;
    vector<float>   *ExtraLepPt;
    vector<float>   *ExtraLepEta;
    vector<float>   *ExtraLepPhi;
    vector<short>   *ExtraLepLepId;
    Float_t ZXFakeweight;
    Float_t ggH_NNLOPS_weight;
    Float_t KFactor_QCD_ggZZ_Nominal;
    Float_t KFactor_QCD_ggZZ_PDFScaleDn;
    Float_t KFactor_QCD_ggZZ_PDFScaleUp;
    Float_t KFactor_QCD_ggZZ_QCDScaleDn;
    Float_t KFactor_QCD_ggZZ_QCDScaleUp;
    Float_t KFactor_QCD_ggZZ_AsDn;
    Float_t KFactor_QCD_ggZZ_AsUp;
    Float_t KFactor_QCD_ggZZ_PDFReplicaDn;
    Float_t KFactor_QCD_ggZZ_PDFReplicaUp;
    Float_t KFactor_EW_qqZZ;
    Float_t KFactor_EW_qqZZ_unc;
    Float_t KFactor_QCD_qqZZ_dPhi;
    Float_t KFactor_QCD_qqZZ_M;
    Float_t KFactor_QCD_qqZZ_Pt;
    
    Float_t PythiaWeight_isr_muR4;
    Float_t PythiaWeight_isr_muR0p25;
    Float_t PythiaWeight_fsr_muR4;
    Float_t PythiaWeight_fsr_muR0p25;

    Short_t genFinalState;
    Int_t   genProcessId;
    Float_t genHEPMCweight;
    Float_t PUWeight;
    Float_t dataMCWeight;
    Float_t trigEffWeight;
    Float_t overallEventWeight;
    Float_t L1prefiringWeight;
    Float_t HqTMCweight;
    Float_t xsec;
    Short_t genExtInfo;
    Float_t GenHMass;
    Float_t GenHPt;
    Float_t GenHRapidity;
    Float_t GenZ1Mass;
    Float_t GenZ1Pt;
    Float_t GenZ1Phi;
    Float_t GenZ1Flav;
    Float_t GenZ2Mass;
    Float_t GenZ2Pt;
    Float_t GenZ2Phi;
    Float_t GenZ2Flav;
    Float_t GenLep1Pt;
    Float_t GenLep1Eta;
    Float_t GenLep1Phi;
    Short_t GenLep1Id;
    Float_t GenLep2Pt;
    Float_t GenLep2Eta;
    Float_t GenLep2Phi;
    Short_t GenLep2Id;
    Float_t GenLep3Pt;
    Float_t GenLep3Eta;
    Float_t GenLep3Phi;
    Short_t GenLep3Id;
    Float_t GenLep4Pt;
    Float_t GenLep4Eta;
    Float_t GenLep4Phi;
    Short_t GenLep4Id;
    Float_t GenAssocLep1Pt;
    Float_t GenAssocLep1Eta;
    Float_t GenAssocLep1Phi;
    Short_t GenAssocLep1Id;
    Float_t GenAssocLep2Pt;
    Float_t GenAssocLep2Eta;
    Float_t GenAssocLep2Phi;
    Short_t GenAssocLep2Id;

    // List of branches
    TBranch *b_RunNumber;
    TBranch *b_EventNumber;
    TBranch *b_LumiNumber;

    TBranch *b_GenMET;
    TBranch *b_GenMETPhi;
    TBranch *b_PFMET;
    TBranch *b_PFMETPhi;
    TBranch *b_METx;
    TBranch *b_METy;
    TBranch *b_METxUPTES;
    TBranch *b_METyUPTES;
    TBranch *b_METxDOWNTES;
    TBranch *b_METyDOWNTES;
    TBranch *b_METxUPEES;
    TBranch *b_METyUPEES;
    TBranch *b_METxDOWNEES;
    TBranch *b_METyDOWNEES;
    TBranch *b_METxUPMES;
    TBranch *b_METyUPMES;
    TBranch *b_METxDOWNMES;
    TBranch *b_METyDOWNMES;
    TBranch *b_METxUPJES;
    TBranch *b_METyUPJES;
    TBranch *b_METxDOWNJES;
    TBranch *b_METyDOWNJES;
    TBranch *b_METxUPJER;
    TBranch *b_METyUPJER;
    TBranch *b_METxDOWNJER;
    TBranch *b_METyDOWNJER;
    
    TBranch *b_nCleanedJets;
    TBranch *b_ZZMass;
    TBranch *b_ZZsel;
    TBranch *b_ZZPt;
    TBranch *b_ZZjjPt;
    TBranch *b_ZZEta;
    TBranch *b_ZZPhi;
    TBranch *b_CRflag;
    
    TBranch *b_Z1Mass;
    TBranch *b_Z1Pt;
    TBranch *b_Z1Flav;
    
    TBranch *b_Z2Mass;
    TBranch *b_Z2Pt;
    TBranch *b_Z2Flav;
    TBranch *b_Z2isGoodTau;
    
    TBranch *b_Z2SVMass;
    TBranch *b_Z2SVPt;
    TBranch *b_Z2SVEta;
    TBranch *b_Z2SVPhi;
    TBranch *b_Z2GoodMass;

    TBranch *b_Z2Mass_TESup;
    TBranch *b_Z2SVMass_TESup;
    TBranch *b_Z2SVPt_TESup;
    TBranch *b_Z2SVEta_TESup;
    TBranch *b_Z2SVPhi_TESup;
    TBranch *b_Z2GoodMass_TESup;

    TBranch *b_Z2Mass_EESup;
    TBranch *b_Z2SVMass_EESup;
    TBranch *b_Z2SVPt_EESup;
    TBranch *b_Z2SVEta_EESup;
    TBranch *b_Z2SVPhi_EESup;
    TBranch *b_Z2GoodMass_EESup;

    TBranch *b_Z2Mass_MESup;
    TBranch *b_Z2SVMass_MESup;
    TBranch *b_Z2SVPt_MESup;
    TBranch *b_Z2SVEta_MESup;
    TBranch *b_Z2SVPhi_MESup;
    TBranch *b_Z2GoodMass_MESup;

    TBranch *b_Z2Mass_JESup;
    TBranch *b_Z2SVMass_JESup;
    TBranch *b_Z2SVPt_JESup;
    TBranch *b_Z2SVEta_JESup;
    TBranch *b_Z2SVPhi_JESup;
    TBranch *b_Z2GoodMass_JESup;

    TBranch *b_Z2Mass_JERup;
    TBranch *b_Z2SVMass_JERup;
    TBranch *b_Z2SVPt_JERup;
    TBranch *b_Z2SVEta_JERup;
    TBranch *b_Z2SVPhi_JERup;
    TBranch *b_Z2GoodMass_JERup;

    TBranch *b_Z2Mass_TESdn;
    TBranch *b_Z2SVMass_TESdn;
    TBranch *b_Z2SVPt_TESdn;
    TBranch *b_Z2SVEta_TESdn;
    TBranch *b_Z2SVPhi_TESdn;
    TBranch *b_Z2GoodMass_TESdn;

    TBranch *b_Z2Mass_EESdn;
    TBranch *b_Z2SVMass_EESdn;
    TBranch *b_Z2SVPt_EESdn;
    TBranch *b_Z2SVEta_EESdn;
    TBranch *b_Z2SVPhi_EESdn;
    TBranch *b_Z2GoodMass_EESdn;

    TBranch *b_Z2Mass_MESdn;
    TBranch *b_Z2SVMass_MESdn;
    TBranch *b_Z2SVPt_MESdn;
    TBranch *b_Z2SVEta_MESdn;
    TBranch *b_Z2SVPhi_MESdn;
    TBranch *b_Z2GoodMass_MESdn;

    TBranch *b_Z2Mass_JESdn;
    TBranch *b_Z2SVMass_JESdn;
    TBranch *b_Z2SVPt_JESdn;
    TBranch *b_Z2SVEta_JESdn;
    TBranch *b_Z2SVPhi_JESdn;
    TBranch *b_Z2GoodMass_JESdn;

    TBranch *b_Z2Mass_JERdn;
    TBranch *b_Z2SVMass_JERdn;
    TBranch *b_Z2SVPt_JERdn;
    TBranch *b_Z2SVEta_JERdn;
    TBranch *b_Z2SVPhi_JERdn;
    TBranch *b_Z2GoodMass_JERdn;
    
    TBranch *b_costhetastar;
    TBranch *b_helphi;
    TBranch *b_helcosthetaZ1;
    TBranch *b_helcosthetaZ2;
    TBranch *b_phistarZ1;
    TBranch *b_phistarZ2;
    TBranch *b_xi;
    TBranch *b_xistar;
    
    
    TBranch *b_LepPt;
    TBranch *b_LepEta;
    TBranch *b_LepPhi;
    TBranch *b_LepLepId;
    TBranch *b_LepSIP;
    TBranch *b_Lepdxy;
    TBranch *b_Lepdz;
    TBranch *b_LepTime;
    TBranch *b_LepisID;
    TBranch *b_LepisLoose;
    TBranch *b_LepBDT;
    TBranch *b_LepMissingHit;
    TBranch *b_LepCombRelIsoPF;
    TBranch *b_LepSF;
    TBranch *b_LepSF_Unc;
    //tau
    TBranch *b_TauVSmu;
    TBranch *b_TauVSe;
    TBranch *b_TauVSjet;
    TBranch *b_TauDecayMode;
    TBranch *b_TauGenMatch;
    TBranch *b_TauTES_p_Up;
    TBranch *b_TauFES_p_Up;

    TBranch *b_fsrPt;
    TBranch *b_fsrEta;
    TBranch *b_fsrPhi;
    TBranch *b_fsrLept;
    TBranch *b_passIsoPreFSR;

    TBranch *b_JetPt;
    TBranch *b_JetEta;
    TBranch *b_JetPhi;
    TBranch *b_JetMass;
    TBranch *b_JetBTagger;
    TBranch *b_JetIsBtagged;
    TBranch *b_JetIsBtaggedWithSF;
    TBranch *b_JetIsBtaggedWithSFUp;
    TBranch *b_JetIsBtaggedWithSFDn;
    TBranch *b_JetQGLikelihood;
    TBranch *b_JetAxis2;
    TBranch *b_JetMult;
    TBranch *b_JetPtD;
    TBranch *b_JetSigma;
    TBranch *b_DiJetMass;
    TBranch *b_DiJetDEta;
    TBranch *b_DiJetFisher;
    TBranch *b_nExtraLep;
    TBranch *b_nExtraZ;
    TBranch *b_ExtraLepPt;
    TBranch *b_ExtraLepEta;
    TBranch *b_ExtraLepPhi;
    TBranch *b_ExtraLepLepId;
    TBranch *b_ZXFakeweight;
    TBranch *b_ggH_NNLOPS_weight;
    TBranch *b_KFactor_QCD_ggZZ_Nominal;
    TBranch *b_KFactor_QCD_ggZZ_PDFScaleDn;
    TBranch *b_KFactor_QCD_ggZZ_PDFScaleUp;
    TBranch *b_KFactor_QCD_ggZZ_QCDScaleDn;
    TBranch *b_KFactor_QCD_ggZZ_QCDScaleUp;
    TBranch *b_KFactor_QCD_ggZZ_AsDn;
    TBranch *b_KFactor_QCD_ggZZ_AsUp;
    TBranch *b_KFactor_QCD_ggZZ_PDFReplicaDn;
    TBranch *b_KFactor_QCD_ggZZ_PDFReplicaUp;
    TBranch *b_KFactor_EW_qqZZ;
    TBranch *b_KFactor_EW_qqZZ_unc;
    TBranch *b_KFactor_QCD_qqZZ_dPhi;
    TBranch *b_KFactor_QCD_qqZZ_M;
    TBranch *b_KFactor_QCD_qqZZ_Pt;
    
    TBranch *b_PythiaWeight_isr_muR4;
    TBranch *b_PythiaWeight_isr_muR0p25;
    TBranch *b_PythiaWeight_fsr_muR4;
    TBranch *b_PythiaWeight_fsr_muR0p25;

    TBranch *b_genFinalState;
    TBranch *b_genProcessId;
    TBranch *b_genHEPMCweight;
    TBranch *b_PUWeight;
    TBranch *b_dataMCWeight;
    TBranch *b_trigEffWeight;
    TBranch *b_overallEventWeight;
    TBranch *b_L1prefiringWeight;
    TBranch *b_HqTMCweight;
    TBranch *b_xsec;
    TBranch *b_genExtInfo;
    TBranch *b_GenHMass;
    TBranch *b_GenHPt;
    TBranch *b_GenHRapidity;
    TBranch *b_GenZ1Mass;
    TBranch *b_GenZ1Pt;
    TBranch *b_GenZ1Phi;
    TBranch *b_GenZ1Flav;
    TBranch *b_GenZ2Mass;
    TBranch *b_GenZ2Pt;
    TBranch *b_GenZ2Phi;
    TBranch *b_GenZ2Flav;
    TBranch *b_GenLep1Pt;
    TBranch *b_GenLep1Eta;
    TBranch *b_GenLep1Phi;
    TBranch *b_GenLep1Id;
    TBranch *b_GenLep2Pt;
    TBranch *b_GenLep2Eta;
    TBranch *b_GenLep2Phi;
    TBranch *b_GenLep2Id;
    TBranch *b_GenLep3Pt;
    TBranch *b_GenLep3Eta;
    TBranch *b_GenLep3Phi;
    TBranch *b_GenLep3Id;
    TBranch *b_GenLep4Pt;
    TBranch *b_GenLep4Eta;
    TBranch *b_GenLep4Phi;
    TBranch *b_GenLep4Id;
    TBranch *b_GenAssocLep1Pt;
    TBranch *b_GenAssocLep1Eta;
    TBranch *b_GenAssocLep1Phi;
    TBranch *b_GenAssocLep1Id;
    TBranch *b_GenAssocLep2Pt;
    TBranch *b_GenAssocLep2Eta;
    TBranch *b_GenAssocLep2Phi;
    TBranch *b_GenAssocLep2Id;


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
    LepisLoose = 0;
    LepBDT = 0;
    LepMissingHit = 0;
    LepCombRelIsoPF = 0;
    LepSF = 0;
    LepSF_Unc = 0;
    
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

    ExtraLepPt = 0;
    ExtraLepEta = 0;
    ExtraLepPhi = 0;
    ExtraLepLepId = 0;
    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
    fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
    fChain->SetBranchAddress("LumiNumber", &LumiNumber, &b_LumiNumber);
    fChain->SetBranchAddress("GenMET", &GenMET, &b_GenMET);
    fChain->SetBranchAddress("GenMETPhi", &GenMETPhi, &b_GenMETPhi);
    fChain->SetBranchAddress("PFMET", &PFMET, &b_PFMET);
    fChain->SetBranchAddress("PFMETPhi", &PFMETPhi, &b_PFMETPhi);
    fChain->SetBranchAddress("METx", &METx, &b_METx);
    fChain->SetBranchAddress("METy", &METy, &b_METy);
    fChain->SetBranchAddress("METxUPTES", &METxUPTES, &b_METxUPTES);
    fChain->SetBranchAddress("METyUPTES", &METyUPTES, &b_METyUPTES);
    fChain->SetBranchAddress("METxDOWNTES", &METxDOWNTES, &b_METxDOWNTES);
    fChain->SetBranchAddress("METyDOWNTES", &METyDOWNTES, &b_METyDOWNTES);
    fChain->SetBranchAddress("METxUPEES", &METxUPEES, &b_METxUPEES);
    fChain->SetBranchAddress("METyUPEES", &METyUPEES, &b_METyUPEES);
    fChain->SetBranchAddress("METxDOWNEES", &METxDOWNEES, &b_METxDOWNEES);
    fChain->SetBranchAddress("METyDOWNEES", &METyDOWNEES, &b_METyDOWNEES);
    fChain->SetBranchAddress("METxUPMES", &METxUPMES, &b_METxUPMES);
    fChain->SetBranchAddress("METyUPMES", &METyUPMES, &b_METyUPMES);
    fChain->SetBranchAddress("METxDOWNMES", &METxDOWNMES, &b_METxDOWNMES);
    fChain->SetBranchAddress("METyDOWNMES", &METyDOWNMES, &b_METyDOWNMES);
    fChain->SetBranchAddress("METxUPJES", &METxUPJES, &b_METxUPJES);
    fChain->SetBranchAddress("METyUPJES", &METyUPJES, &b_METyUPJES);
    fChain->SetBranchAddress("METxDOWNJES", &METxDOWNJES, &b_METxDOWNJES);
    fChain->SetBranchAddress("METyDOWNJES", &METyDOWNJES, &b_METyDOWNJES);
    fChain->SetBranchAddress("METxUPJER", &METxUPJER, &b_METxUPJER);
    fChain->SetBranchAddress("METyUPJER", &METyUPJER, &b_METyUPJER);
    fChain->SetBranchAddress("METxDOWNJER", &METxDOWNJER, &b_METxDOWNJER);
    fChain->SetBranchAddress("METyDOWNJER", &METyDOWNJER, &b_METyDOWNJER);
    fChain->SetBranchAddress("nCleanedJets", &nCleanedJets, &b_nCleanedJets);
    fChain->SetBranchAddress("ZZMass", &ZZMass, &b_ZZMass);
    fChain->SetBranchAddress("ZZsel", &ZZsel, &b_ZZsel);
    fChain->SetBranchAddress("ZZPt", &ZZPt, &b_ZZPt);
    fChain->SetBranchAddress("ZZjjPt", &ZZjjPt, &b_ZZjjPt);
    fChain->SetBranchAddress("ZZEta", &ZZEta, &b_ZZEta);
    fChain->SetBranchAddress("ZZPhi", &ZZPhi, &b_ZZPhi);
    fChain->SetBranchAddress("CRflag", &CRflag, &b_CRflag);
    fChain->SetBranchAddress("Z1Mass", &Z1Mass, &b_Z1Mass);
    fChain->SetBranchAddress("Z1Pt", &Z1Pt, &b_Z1Pt);
    fChain->SetBranchAddress("Z1Flav", &Z1Flav, &b_Z1Flav);
    fChain->SetBranchAddress("Z2Mass", &Z2Mass, &b_Z2Mass);
    fChain->SetBranchAddress("Z2Pt", &Z2Pt, &b_Z2Pt);
    fChain->SetBranchAddress("Z2Flav", &Z2Flav, &b_Z2Flav);
    fChain->SetBranchAddress("Z2isGoodTau", &Z2isGoodTau, &b_Z2isGoodTau);
    fChain->SetBranchAddress("Z2SVMass", &Z2SVMass, &b_Z2SVMass);
    fChain->SetBranchAddress("Z2SVPt", &Z2SVPt, &b_Z2SVPt);
    fChain->SetBranchAddress("Z2SVEta", &Z2SVEta, &b_Z2SVEta);
    fChain->SetBranchAddress("Z2SVPhi", &Z2SVPhi, &b_Z2SVPhi);
    fChain->SetBranchAddress("Z2GoodMass", &Z2GoodMass, &b_Z2GoodMass);
    fChain->SetBranchAddress("Z2Mass_TESup", &Z2Mass_TESup, &b_Z2Mass_TESup);
    fChain->SetBranchAddress("Z2SVMass_TESup", &Z2SVMass_TESup, &b_Z2SVMass_TESup);
    fChain->SetBranchAddress("Z2SVPt_TESup", &Z2SVPt_TESup, &b_Z2SVPt_TESup);
    fChain->SetBranchAddress("Z2SVEta_TESup", &Z2SVEta_TESup, &b_Z2SVEta_TESup);
    fChain->SetBranchAddress("Z2SVPhi_TESup", &Z2SVPhi_TESup, &b_Z2SVPhi_TESup);
    fChain->SetBranchAddress("Z2GoodMass_TESup", &Z2GoodMass_TESup, &b_Z2GoodMass_TESup);
    fChain->SetBranchAddress("Z2Mass_EESup", &Z2Mass_EESup, &b_Z2Mass_EESup);
    fChain->SetBranchAddress("Z2SVMass_EESup", &Z2SVMass_EESup, &b_Z2SVMass_EESup);
    fChain->SetBranchAddress("Z2SVPt_EESup", &Z2SVPt_EESup, &b_Z2SVPt_EESup);
    fChain->SetBranchAddress("Z2SVEta_EESup", &Z2SVEta_EESup, &b_Z2SVEta_EESup);
    fChain->SetBranchAddress("Z2SVPhi_EESup", &Z2SVPhi_EESup, &b_Z2SVPhi_EESup);
    fChain->SetBranchAddress("Z2GoodMass_EESup", &Z2GoodMass_EESup, &b_Z2GoodMass_EESup);
    fChain->SetBranchAddress("Z2Mass_MESup", &Z2Mass_MESup, &b_Z2Mass_MESup);
    fChain->SetBranchAddress("Z2SVMass_MESup", &Z2SVMass_MESup, &b_Z2SVMass_MESup);
    fChain->SetBranchAddress("Z2SVPt_MESup", &Z2SVPt_MESup, &b_Z2SVPt_MESup);
    fChain->SetBranchAddress("Z2SVEta_MESup", &Z2SVEta_MESup, &b_Z2SVEta_MESup);
    fChain->SetBranchAddress("Z2SVPhi_MESup", &Z2SVPhi_MESup, &b_Z2SVPhi_MESup);
    fChain->SetBranchAddress("Z2GoodMass_MESup", &Z2GoodMass_MESup, &b_Z2GoodMass_MESup);
    fChain->SetBranchAddress("Z2Mass_JESup", &Z2Mass_JESup, &b_Z2Mass_JESup);
    fChain->SetBranchAddress("Z2SVMass_JESup", &Z2SVMass_JESup, &b_Z2SVMass_JESup);
    fChain->SetBranchAddress("Z2SVPt_JESup", &Z2SVPt_JESup, &b_Z2SVPt_JESup);
    fChain->SetBranchAddress("Z2SVEta_JESup", &Z2SVEta_JESup, &b_Z2SVEta_JESup);
    fChain->SetBranchAddress("Z2SVPhi_JESup", &Z2SVPhi_JESup, &b_Z2SVPhi_JESup);
    fChain->SetBranchAddress("Z2GoodMass_JESup", &Z2GoodMass_JESup, &b_Z2GoodMass_JESup);
    fChain->SetBranchAddress("Z2Mass_JERup", &Z2Mass_JERup, &b_Z2Mass_JERup);
    fChain->SetBranchAddress("Z2SVMass_JERup", &Z2SVMass_JERup, &b_Z2SVMass_JERup);
    fChain->SetBranchAddress("Z2SVPt_JERup", &Z2SVPt_JERup, &b_Z2SVPt_JERup);
    fChain->SetBranchAddress("Z2SVEta_JERup", &Z2SVEta_JERup, &b_Z2SVEta_JERup);
    fChain->SetBranchAddress("Z2SVPhi_JERup", &Z2SVPhi_JERup, &b_Z2SVPhi_JERup);
    fChain->SetBranchAddress("Z2GoodMass_JERup", &Z2GoodMass_JERup, &b_Z2GoodMass_JERup);
    fChain->SetBranchAddress("Z2Mass_TESdn", &Z2Mass_TESdn, &b_Z2Mass_TESdn);
    fChain->SetBranchAddress("Z2SVMass_TESdn", &Z2SVMass_TESdn, &b_Z2SVMass_TESdn);
    fChain->SetBranchAddress("Z2SVPt_TESdn", &Z2SVPt_TESdn, &b_Z2SVPt_TESdn);
    fChain->SetBranchAddress("Z2SVEta_TESdn", &Z2SVEta_TESdn, &b_Z2SVEta_TESdn);
    fChain->SetBranchAddress("Z2SVPhi_TESdn", &Z2SVPhi_TESdn, &b_Z2SVPhi_TESdn);
    fChain->SetBranchAddress("Z2GoodMass_TESdn", &Z2GoodMass_TESdn, &b_Z2GoodMass_TESdn);
    fChain->SetBranchAddress("Z2Mass_EESdn", &Z2Mass_EESdn, &b_Z2Mass_EESdn);
    fChain->SetBranchAddress("Z2SVMass_EESdn", &Z2SVMass_EESdn, &b_Z2SVMass_EESdn);
    fChain->SetBranchAddress("Z2SVPt_EESdn", &Z2SVPt_EESdn, &b_Z2SVPt_EESdn);
    fChain->SetBranchAddress("Z2SVEta_EESdn", &Z2SVEta_EESdn, &b_Z2SVEta_EESdn);
    fChain->SetBranchAddress("Z2SVPhi_EESdn", &Z2SVPhi_EESdn, &b_Z2SVPhi_EESdn);
    fChain->SetBranchAddress("Z2GoodMass_EESdn", &Z2GoodMass_EESdn, &b_Z2GoodMass_EESdn);
    fChain->SetBranchAddress("Z2Mass_MESdn", &Z2Mass_MESdn, &b_Z2Mass_MESdn);
    fChain->SetBranchAddress("Z2SVMass_MESdn", &Z2SVMass_MESdn, &b_Z2SVMass_MESdn);
    fChain->SetBranchAddress("Z2SVPt_MESdn", &Z2SVPt_MESdn, &b_Z2SVPt_MESdn);
    fChain->SetBranchAddress("Z2SVEta_MESdn", &Z2SVEta_MESdn, &b_Z2SVEta_MESdn);
    fChain->SetBranchAddress("Z2SVPhi_MESdn", &Z2SVPhi_MESdn, &b_Z2SVPhi_MESdn);
    fChain->SetBranchAddress("Z2GoodMass_MESdn", &Z2GoodMass_MESdn, &b_Z2GoodMass_MESdn);
    fChain->SetBranchAddress("Z2Mass_JESdn", &Z2Mass_JESdn, &b_Z2Mass_JESdn);
    fChain->SetBranchAddress("Z2SVMass_JESdn", &Z2SVMass_JESdn, &b_Z2SVMass_JESdn);
    fChain->SetBranchAddress("Z2SVPt_JESdn", &Z2SVPt_JESdn, &b_Z2SVPt_JESdn);
    fChain->SetBranchAddress("Z2SVEta_JESdn", &Z2SVEta_JESdn, &b_Z2SVEta_JESdn);
    fChain->SetBranchAddress("Z2SVPhi_JESdn", &Z2SVPhi_JESdn, &b_Z2SVPhi_JESdn);
    fChain->SetBranchAddress("Z2GoodMass_JESdn", &Z2GoodMass_JESdn, &b_Z2GoodMass_JESdn);
    fChain->SetBranchAddress("Z2Mass_JERdn", &Z2Mass_JERdn, &b_Z2Mass_JERdn);
    fChain->SetBranchAddress("Z2SVMass_JERdn", &Z2SVMass_JERdn, &b_Z2SVMass_JERdn);
    fChain->SetBranchAddress("Z2SVPt_JERdn", &Z2SVPt_JERdn, &b_Z2SVPt_JERdn);
    fChain->SetBranchAddress("Z2SVEta_JERdn", &Z2SVEta_JERdn, &b_Z2SVEta_JERdn);
    fChain->SetBranchAddress("Z2SVPhi_JERdn", &Z2SVPhi_JERdn, &b_Z2SVPhi_JERdn);
    fChain->SetBranchAddress("Z2GoodMass_JERdn", &Z2GoodMass_JERdn, &b_Z2GoodMass_JERdn);
    fChain->SetBranchAddress("costhetastar", &costhetastar, &b_costhetastar);
    fChain->SetBranchAddress("helphi", &helphi, &b_helphi);
    fChain->SetBranchAddress("helcosthetaZ1", &helcosthetaZ1, &b_helcosthetaZ1);
    fChain->SetBranchAddress("helcosthetaZ2", &helcosthetaZ2, &b_helcosthetaZ2);
    fChain->SetBranchAddress("phistarZ1", &phistarZ1, &b_phistarZ1);
    fChain->SetBranchAddress("phistarZ2", &phistarZ2, &b_phistarZ2);
    fChain->SetBranchAddress("xi", &xi, &b_xi);
    fChain->SetBranchAddress("xistar", &xistar, &b_xistar);
    fChain->SetBranchAddress("LepPt", &LepPt, &b_LepPt);
    fChain->SetBranchAddress("LepEta", &LepEta, &b_LepEta);
    fChain->SetBranchAddress("LepPhi", &LepPhi, &b_LepPhi);
    fChain->SetBranchAddress("LepLepId", &LepLepId, &b_LepLepId);
    fChain->SetBranchAddress("LepSIP", &LepSIP, &b_LepSIP);
    fChain->SetBranchAddress("Lepdxy", &Lepdxy, &b_Lepdxy);
    fChain->SetBranchAddress("Lepdz", &Lepdz, &b_Lepdz);
    fChain->SetBranchAddress("LepTime", &LepTime, &b_LepTime);
    fChain->SetBranchAddress("LepisID", &LepisID, &b_LepisID);
    fChain->SetBranchAddress("LepisLoose", &LepisLoose, &b_LepisLoose);
    fChain->SetBranchAddress("LepBDT", &LepBDT, &b_LepBDT);
    fChain->SetBranchAddress("LepMissingHit", &LepMissingHit, &b_LepMissingHit);
    fChain->SetBranchAddress("LepCombRelIsoPF", &LepCombRelIsoPF, &b_LepCombRelIsoPF);
    fChain->SetBranchAddress("LepSF", &LepSF, &b_LepSF);
    fChain->SetBranchAddress("LepSF_Unc", &LepSF_Unc, &b_LepSF_Unc);
    fChain->SetBranchAddress("TauVSmu", &TauVSmu, &b_TauVSmu);
    fChain->SetBranchAddress("TauVSe", &TauVSe, &b_TauVSe);
    fChain->SetBranchAddress("TauVSjet", &TauVSjet, &b_TauVSjet);
    fChain->SetBranchAddress("TauDecayMode", &TauDecayMode, &b_TauDecayMode);
    fChain->SetBranchAddress("TauGenMatch", &TauGenMatch, &b_TauGenMatch);
    fChain->SetBranchAddress("TauTES_p_Up", &TauTES_p_Up, &b_TauTES_p_Up);
    fChain->SetBranchAddress("TauFES_p_Up", &TauFES_p_Up, &b_TauFES_p_Up);
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
    fChain->SetBranchAddress("DiJetMass", &DiJetMass, &b_DiJetMass);
    fChain->SetBranchAddress("DiJetDEta", &DiJetDEta, &b_DiJetDEta);
    fChain->SetBranchAddress("DiJetFisher", &DiJetFisher, &b_DiJetFisher);
    fChain->SetBranchAddress("nExtraLep", &nExtraLep, &b_nExtraLep);
    fChain->SetBranchAddress("nExtraZ", &nExtraZ, &b_nExtraZ);
    fChain->SetBranchAddress("ExtraLepPt", &ExtraLepPt, &b_ExtraLepPt);
    fChain->SetBranchAddress("ExtraLepEta", &ExtraLepEta, &b_ExtraLepEta);
    fChain->SetBranchAddress("ExtraLepPhi", &ExtraLepPhi, &b_ExtraLepPhi);
    fChain->SetBranchAddress("ExtraLepLepId", &ExtraLepLepId, &b_ExtraLepLepId);
    fChain->SetBranchAddress("ZXFakeweight", &ZXFakeweight, &b_ZXFakeweight);
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
      
      if ( input_file_name.Contains("ZZTo4l") )
      {
         fChain->SetBranchAddress("KFactor_EW_qqZZ", &KFactor_EW_qqZZ, &b_KFactor_EW_qqZZ);
         fChain->SetBranchAddress("KFactor_EW_qqZZ_unc", &KFactor_EW_qqZZ_unc, &b_KFactor_EW_qqZZ_unc);
         fChain->SetBranchAddress("KFactor_QCD_qqZZ_dPhi", &KFactor_QCD_qqZZ_dPhi, &b_KFactor_QCD_qqZZ_dPhi);
         fChain->SetBranchAddress("KFactor_QCD_qqZZ_M", &KFactor_QCD_qqZZ_M, &b_KFactor_QCD_qqZZ_M);
         fChain->SetBranchAddress("KFactor_QCD_qqZZ_Pt", &KFactor_QCD_qqZZ_Pt, &b_KFactor_QCD_qqZZ_Pt);
      }
      
    fChain->SetBranchAddress("PythiaWeight_isr_muR4", &PythiaWeight_isr_muR4, &b_PythiaWeight_isr_muR4);
    fChain->SetBranchAddress("PythiaWeight_isr_muR0p25", &PythiaWeight_isr_muR0p25, &b_PythiaWeight_isr_muR0p25);
    fChain->SetBranchAddress("PythiaWeight_fsr_muR4", &PythiaWeight_fsr_muR4, &b_PythiaWeight_fsr_muR4);
    fChain->SetBranchAddress("PythiaWeight_fsr_muR0p25", &PythiaWeight_fsr_muR0p25, &b_PythiaWeight_fsr_muR0p25);
    fChain->SetBranchAddress("genFinalState", &genFinalState, &b_genFinalState);
    fChain->SetBranchAddress("genProcessId", &genProcessId, &b_genProcessId);
    fChain->SetBranchAddress("genHEPMCweight", &genHEPMCweight, &b_genHEPMCweight);
    fChain->SetBranchAddress("PUWeight", &PUWeight, &b_PUWeight);
    fChain->SetBranchAddress("dataMCWeight", &dataMCWeight, &b_dataMCWeight);
    fChain->SetBranchAddress("trigEffWeight", &trigEffWeight, &b_trigEffWeight);
    fChain->SetBranchAddress("overallEventWeight", &overallEventWeight, &b_overallEventWeight);
    fChain->SetBranchAddress("L1prefiringWeight", &L1prefiringWeight, &b_L1prefiringWeight);
    fChain->SetBranchAddress("HqTMCweight", &HqTMCweight, &b_HqTMCweight);
    fChain->SetBranchAddress("xsec", &xsec, &b_xsec);
    fChain->SetBranchAddress("genExtInfo", &genExtInfo, &b_genExtInfo);
    fChain->SetBranchAddress("GenHMass", &GenHMass, &b_GenHMass);
    fChain->SetBranchAddress("GenHPt", &GenHPt, &b_GenHPt);
    fChain->SetBranchAddress("GenHRapidity", &GenHRapidity, &b_GenHRapidity);
    fChain->SetBranchAddress("GenZ1Mass", &GenZ1Mass, &b_GenZ1Mass);
    fChain->SetBranchAddress("GenZ1Pt", &GenZ1Pt, &b_GenZ1Pt);
    fChain->SetBranchAddress("GenZ1Phi", &GenZ1Phi, &b_GenZ1Phi);
    fChain->SetBranchAddress("GenZ1Flav", &GenZ1Flav, &b_GenZ1Flav);
    fChain->SetBranchAddress("GenZ2Mass", &GenZ2Mass, &b_GenZ2Mass);
    fChain->SetBranchAddress("GenZ2Pt", &GenZ2Pt, &b_GenZ2Pt);
    fChain->SetBranchAddress("GenZ2Phi", &GenZ2Phi, &b_GenZ2Phi);
    fChain->SetBranchAddress("GenZ2Flav", &GenZ2Flav, &b_GenZ2Flav);
    fChain->SetBranchAddress("GenLep1Pt", &GenLep1Pt, &b_GenLep1Pt);
    fChain->SetBranchAddress("GenLep1Eta", &GenLep1Eta, &b_GenLep1Eta);
    fChain->SetBranchAddress("GenLep1Phi", &GenLep1Phi, &b_GenLep1Phi);
    fChain->SetBranchAddress("GenLep1Id", &GenLep1Id, &b_GenLep1Id);
    fChain->SetBranchAddress("GenLep2Pt", &GenLep2Pt, &b_GenLep2Pt);
    fChain->SetBranchAddress("GenLep2Eta", &GenLep2Eta, &b_GenLep2Eta);
    fChain->SetBranchAddress("GenLep2Phi", &GenLep2Phi, &b_GenLep2Phi);
    fChain->SetBranchAddress("GenLep2Id", &GenLep2Id, &b_GenLep2Id);
    fChain->SetBranchAddress("GenLep3Pt", &GenLep3Pt, &b_GenLep3Pt);
    fChain->SetBranchAddress("GenLep3Eta", &GenLep3Eta, &b_GenLep3Eta);
    fChain->SetBranchAddress("GenLep3Phi", &GenLep3Phi, &b_GenLep3Phi);
    fChain->SetBranchAddress("GenLep3Id", &GenLep3Id, &b_GenLep3Id);
    fChain->SetBranchAddress("GenLep4Pt", &GenLep4Pt, &b_GenLep4Pt);
    fChain->SetBranchAddress("GenLep4Eta", &GenLep4Eta, &b_GenLep4Eta);
    fChain->SetBranchAddress("GenLep4Phi", &GenLep4Phi, &b_GenLep4Phi);
    fChain->SetBranchAddress("GenLep4Id", &GenLep4Id, &b_GenLep4Id);
    fChain->SetBranchAddress("GenAssocLep1Pt", &GenAssocLep1Pt, &b_GenAssocLep1Pt);
    fChain->SetBranchAddress("GenAssocLep1Eta", &GenAssocLep1Eta, &b_GenAssocLep1Eta);
    fChain->SetBranchAddress("GenAssocLep1Phi", &GenAssocLep1Phi, &b_GenAssocLep1Phi);
    fChain->SetBranchAddress("GenAssocLep1Id", &GenAssocLep1Id, &b_GenAssocLep1Id);
    fChain->SetBranchAddress("GenAssocLep2Pt", &GenAssocLep2Pt, &b_GenAssocLep2Pt);
    fChain->SetBranchAddress("GenAssocLep2Eta", &GenAssocLep2Eta, &b_GenAssocLep2Eta);
    fChain->SetBranchAddress("GenAssocLep2Phi", &GenAssocLep2Phi, &b_GenAssocLep2Phi);
    fChain->SetBranchAddress("GenAssocLep2Id", &GenAssocLep2Id, &b_GenAssocLep2Id);

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
