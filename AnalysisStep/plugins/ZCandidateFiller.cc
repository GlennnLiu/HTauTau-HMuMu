/** \class ZCandidateFiller
 *
 *  Create a collection of Zs including FSR and additional variables stored as userFloats.
 *
 *  \author N. Amapane (Torino)
 *  \author S. Bolognesi (JHU)
 *  \author C. Botta (CERN)
 *  \author S. Casasso (Torino)
 */

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h>
#include <FWCore/Utilities/interface/InputTag.h>

#include <CommonTools/Utils/interface/StringCutObjectSelector.h>
#include <CommonTools/Utils/interface/StringObjectFunction.h>

#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/Candidate/interface/ShallowCloneCandidate.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Tau.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/METReco/interface/PFMET.h>
#include <DataFormats/METReco/interface/PFMETCollection.h>
#include <DataFormats/METReco/interface/CommonMETData.h>

#include <TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h>
#include <TauAnalysis/ClassicSVfit/interface/ClassicSVfitIntegrand.h>
#include <TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h>
#include <TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h>

#include <HTauTauHMuMu/AnalysisStep/interface/PhotonFwd.h>
#include <HTauTauHMuMu/AnalysisStep/interface/CutSet.h>
#include <HTauTauHMuMu/AnalysisStep/interface/LeptonIsoHelper.h>
#include <HTauTauHMuMu/AnalysisStep/interface/DaughterDataHelpers.h>

#include <TLorentzVector.h>
#include <string>
#include <Math/VectorUtil.h>
#include <vector>
#include <string>
#include <cmath>

using namespace edm;
using namespace std;
using namespace reco;
using namespace classic_svFit;
using METUncertainty = pat::MET::METUncertainty;


class ZCandidateFiller : public edm::EDProducer {
 public:
  /// Constructor
  explicit ZCandidateFiller(const edm::ParameterSet&);
    
  /// Destructor
  ~ZCandidateFiller(){};  

 private:
  virtual void beginJob(){};  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};

  classic_svFit::MeasuredTauLepton::kDecayType GetDecayTypeFlag (int pdgId);
  bool Switch (classic_svFit::MeasuredTauLepton::kDecayType type1, double pt1, float l1_iso, classic_svFit::MeasuredTauLepton::kDecayType type2, double pt2, float l2_iso);
  double GetMass (classic_svFit::MeasuredTauLepton::kDecayType type, double candMass);
  bool IsInteresting (const reco::Candidate *l1, const reco::Candidate *l2); // if true, compute SVFit
  
  //  void ComputeSVFit(pat::CompositeCandidate * candPtr,
  //		    double METx, double METy, double uncorrMETx, double uncorrMETy, TMatrixD covMET, float significance,
  //		    double METx_UP_JES, double METy_UP_JES, double METx_DOWN_JES, double METy_DOWN_JES,
  //		    double METx_UP_TES, double METy_UP_TES, double METx_DOWN_TES, double METy_DOWN_TES,
  //		    double METx_UP_EES, double METy_UP_EES, double METx_DOWN_EES, double METy_DOWN_EES
  //		    );

  edm::EDGetTokenT<edm::View<reco::CompositeCandidate> > candidateToken;
  const StringCutObjectSelector<pat::CompositeCandidate, true> preBestZSelection;
  int sampleType;
  int setup;
  const CutSet<pat::CompositeCandidate> cuts;
  int FSRMode;
  bool embedDaughterFloats;

  //For SV fit
  edm::EDGetTokenT<View<pat::MET> > theMETTag;
  edm::EDGetTokenT<double> theSigTag;
  edm::EDGetTokenT<math::Error<2>::type> theCovTag;
  bool _usePairMET;
  bool _computeForUpDownTES;
  bool _computeForUpDownMET;
  edm::EDGetTokenT<double> theMETdxUPTag;
  edm::EDGetTokenT<double> theMETdyUPTag;
  edm::EDGetTokenT<double> theMETdxDOWNTag;
  edm::EDGetTokenT<double> theMETdyDOWNTag;
  edm::EDGetTokenT<double> theMETdxUPEESTag;
  edm::EDGetTokenT<double> theMETdyUPEESTag;
  edm::EDGetTokenT<double> theMETdxDOWNEESTag;
  edm::EDGetTokenT<double> theMETdyDOWNEESTag;
  
  //For HLT match
  vector<string> muHLTPaths2_;
  vector<string> eleHLTPaths2_;

  enum pairType {
    kMuHad  = 0,
    kEHad   = 1,
    kHadHad = 2,
    kMuMu   = 3,
    kEE     = 4,
    kEMu    = 5,
    kEEPrompt = 6, // prompt Z->ee/mumu decays
    kMuMuPrompt = 7,
    kOther  = 8 // for e.g. h->bb
  };


};


ZCandidateFiller::ZCandidateFiller(const edm::ParameterSet& iConfig) :
  candidateToken(consumes<edm::View<reco::CompositeCandidate> >(iConfig.getParameter<edm::InputTag>("src"))),
  preBestZSelection(iConfig.getParameter<std::string>("bestZAmong")),
  sampleType(iConfig.getParameter<int>("sampleType")),
  setup(iConfig.getParameter<int>("setup")),
  cuts(iConfig.getParameter<edm::ParameterSet>("flags")),
  embedDaughterFloats(iConfig.getUntrackedParameter<bool>("embedDaughterFloats",true)),
  theMETTag(consumes<View<pat::MET>>(iConfig.getParameter<edm::InputTag>("srcMET"))),
  theSigTag(consumes<double>(iConfig.getParameter<edm::InputTag>("srcSig"))),
  theCovTag(consumes<math::Error<2>::type>(iConfig.getParameter<edm::InputTag>("srcCov"))),
  theMETdxUPTag(consumes<double>(iConfig.getParameter<edm::InputTag>("METdxUP"))),
  theMETdyUPTag(consumes<double>(iConfig.getParameter<edm::InputTag>("METdyUP"))),
  theMETdxDOWNTag(consumes<double>(iConfig.getParameter<edm::InputTag>("METdxDOWN"))),
  theMETdyDOWNTag(consumes<double>(iConfig.getParameter<edm::InputTag>("METdyDOWN"))),
  theMETdxUPEESTag(consumes<double>(iConfig.getParameter<edm::InputTag>("METdxUP_EES"))),
  theMETdyUPEESTag(consumes<double>(iConfig.getParameter<edm::InputTag>("METdyUP_EES"))),
  theMETdxDOWNEESTag(consumes<double>(iConfig.getParameter<edm::InputTag>("METdxDOWN_EES"))),
  theMETdyDOWNEESTag(consumes<double>(iConfig.getParameter<edm::InputTag>("METdyDOWN_EES")))
{
  
  string mode = iConfig.getParameter<string>("FSRMode");
  if      (mode == "skip")   FSRMode = 0;
  else if (mode == "Legacy") FSRMode = 1;
  else if (mode == "RunII")  FSRMode = 2;
  else {
    cout << "ZCandidateFiller: FSRMode " << FSRMode << " not supported" << endl;
    abort();
  }
  
  _usePairMET = iConfig.getParameter<bool>("usePairMET");
  _computeForUpDownTES = false;
  _computeForUpDownMET = false;

  if (sampleType == 2016)
  {
	muHLTPaths2_ = 
	{
	"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",//DiMu
	"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*",
	"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*",
	"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*",
	};
	eleHLTPaths2_ = 
	{
	"HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
	"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
	"HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*",
	};
  }
  else if (sampleType == 2017)
  {
	muHLTPaths2_ =
	{
	"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*",
	"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*",
	};
	eleHLTPaths2_ = 
	{
	"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
	"HLT_DoubleEle33_CaloIdL_MW_v*",
	};
  }
  else if (sampleType == 2018)
  {
	muHLTPaths2_ = 
	{
	"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*",
	};
	eleHLTPaths2_ = 
	{
	"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
	"HLT_DoubleEle25_CaloIdL_MW_v*",
	};
  }
  
  produces<pat::CompositeCandidateCollection>();
}


void
ZCandidateFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{  
  using namespace edm;
  using namespace std;
  using namespace reco;
  
  auto result = std::make_unique<pat::CompositeCandidateCollection>();

  //-- Get LL candidates
  Handle<View<reco::CompositeCandidate> > LLCands;
  iEvent.getByToken(candidateToken, LLCands);

  //SV fit info
  Handle<View<pat::MET> > METHandle;


  // intialize MET
  double METx = 0.;
  double METy = 0.; 
  double uncorrMETx = -999.;
  double uncorrMETy = -999.; 
  TMatrixD covMET(2, 2);
  float significance = -999.;

  double METx_UP_JES   = 0.;
  double METy_UP_JES   = 0.;
  double METx_DOWN_JES = 0.;
  double METy_DOWN_JES = 0.;

  double METx_UP_TES = 0.;
  double METy_UP_TES = 0.;
  double METx_DOWN_TES = 0.;
  double METy_DOWN_TES = 0.;

  double METx_UP_EES = 0.;
  double METy_UP_EES = 0.;
  double METx_DOWN_EES = 0.;
  double METy_DOWN_EES = 0.;

  iEvent.getByToken(theMETTag, METHandle);

  if (!_usePairMET)
  {   
     unsigned int metNumber = METHandle->size();
     if (metNumber != 1)     
        edm::LogWarning("pfMetHasNotSizeOne") << "(ClassicSVfitInterface) Warning! Using single pf MEt, but input MEt collection size is different from 1"
                                                           << "   --> using MET entry num. 0";
     const pat::MET& patMET = (*METHandle)[0];
     METx = patMET.px();
     METy = patMET.py();

     Handle<double> significanceHandle;
     Handle<math::Error<2>::type> covHandle;
     Handle<double> METxUPHandle;
     Handle<double> METyUPHandle;
     Handle<double> METxDOWNHandle;
     Handle<double> METyDOWNHandle;
     Handle<double> METxUPEESHandle;
     Handle<double> METyUPEESHandle;
     Handle<double> METxDOWNEESHandle;
     Handle<double> METyDOWNEESHandle;

     iEvent.getByToken (theSigTag, significanceHandle);
     iEvent.getByToken (theCovTag, covHandle);
     
     iEvent.getByToken (theMETdxUPTag, METxUPHandle);
     iEvent.getByToken (theMETdyUPTag, METyUPHandle);
     iEvent.getByToken (theMETdxDOWNTag, METxDOWNHandle);
     iEvent.getByToken (theMETdyDOWNTag, METyDOWNHandle);
     iEvent.getByToken (theMETdxUPEESTag, METxUPEESHandle);
     iEvent.getByToken (theMETdyUPEESTag, METyUPEESHandle);
     iEvent.getByToken (theMETdxDOWNEESTag, METxDOWNEESHandle);
     iEvent.getByToken (theMETdyDOWNEESTag, METyDOWNEESHandle);

     covMET[0][0] = (*covHandle)(0,0);
     covMET[1][0] = (*covHandle)(1,0);
     covMET[0][1] = covMET[1][0]; // (1,0) is the only one saved
     covMET[1][1] = (*covHandle)(1,1);

     significance = (float) (*significanceHandle);
     
     METx_UP_TES   = *METxUPHandle;
     METy_UP_TES   = *METyUPHandle;
     METx_DOWN_TES = *METxDOWNHandle;
     METy_DOWN_TES = *METyDOWNHandle;
     METx_UP_EES   = *METxUPEESHandle;
     METy_UP_EES   = *METyUPEESHandle;
     METx_DOWN_EES = *METxDOWNEESHandle;
     METy_DOWN_EES = *METyDOWNEESHandle;

     // protection against singular matrices
     if (covMET[0][0] == 0 && covMET[1][0] == 0 && covMET[0][1] == 0 && covMET[1][1] == 0)
        edm::LogWarning("SingularCovarianceMatrix") << "(ClassicSVfitInterface) Warning! Input covariance matrix is singular"
                                                    << "   --> SVfit algorithm will probably crash...";

     // now get the MET shifted UP/DOWN by the JES
     LorentzVector patMET_UP_JES   = patMET.shiftedP4(METUncertainty::JetEnUp);
     METx_UP_JES = patMET_UP_JES.px();
     METy_UP_JES = patMET_UP_JES.py();
     LorentzVector patMET_DOWN_JES = patMET.shiftedP4(METUncertainty::JetEnDown);
     METx_DOWN_JES = patMET_DOWN_JES.px();
     METy_DOWN_JES = patMET_DOWN_JES.py();
  }

  //--- Fill user info
  const float ZmassValue = 91.1876;

  float closestZeeMassDiff = 99999.;
  float closestZmmMassDiff = 99999.;
  float closestZttMassDiff = 99999.;
  float closestZMassDiff   = 99999.;
  float closestLLMassDiff = 99999.;

  int bestZeeidx = -1;
  int bestZmmidx = -1;
  int bestZttidx = -1;
  int bestZidx = -1;
  int bestLLidx = -1;

  //--- Loop over LL Candidates
  cout<<LLCands->size()<<" bare Z candidates"<<endl;
  for(unsigned int i = 0; i < LLCands->size(); ++i) {
    const CompositeCandidate& c = (*LLCands)[i];
    pat::CompositeCandidate myCand(c); 

    if (embedDaughterFloats){  
      userdatahelpers::embedDaughterData(myCand);
    }

    int id0 = myCand.daughter(0)->pdgId();
    int id1 = myCand.daughter(1)->pdgId();
    bool OS = (id0*id1)<0;
    bool SF = ( abs(id0)!=15 && abs(id1)!=15 && abs(id0)==abs(id1) ) || abs(id0)==15 || abs(id0)==15;//ee, mumu, etau, mutau, tautau

    // ------------------------------
    // FSR recovery
    // ------------------------------

    if (FSRMode==1) { // Legacy
      //loop on the 2 daughters; apply mass cuts on (llg) and store 
      // the highest-pT and the lowest-DR assocated gamma.
      double    maxPT = -1.;
      double    minDR = 9999.;
      const pat::PFParticle* maxPTg=0;
      const pat::PFParticle* minDRg=0;
      int maxPTgLep=-1; // Index of daughter to which the above photons
      int minDRgLep=-1; // are associated

      for (int dauIdx=0; dauIdx<2; ++dauIdx) { 
        const Candidate* d = myCand.daughter(dauIdx);
        const PhotonPtrVector* gammas = userdatahelpers::getUserPhotons(d);
        if (gammas==0) continue;
        for (PhotonPtrVector::const_iterator g = gammas->begin();
             g!= gammas->end(); ++g) {
          const pat::PFParticle* gamma = g->get();
          reco::Candidate::LorentzVector p4G = gamma->p4();
          reco::Candidate::LorentzVector p4LL = myCand.p4();
          double mLLG = (p4LL + p4G).M();
          bool movesToZPeak = (fabs(mLLG-ZmassValue) < fabs(myCand.mass()-ZmassValue));
          if (movesToZPeak && mLLG<100. && mLLG>4) { // Mass cuts (4 is implicit)

            double pt = gamma->pt();
            if (pt>maxPT) {
              maxPT  = pt;
              maxPTg = gamma;
              maxPTgLep = dauIdx;
            }

            double dR = ROOT::Math::VectorUtil::DeltaR(gamma->momentum(),d->momentum());
            if (dR<minDR) {
              minDR  = dR;
              minDRg = gamma;
              minDRgLep = dauIdx;
            }
          }
        } // end loop on photons
      } // end loop on daughters (leptons)
    
      // Define the selected FSR photon.
      const pat::PFParticle* fsr=0;    
      int lepWithFsr=-1; 
      if (maxPTg!=0) { // at least 1 photon selected
        if (maxPT>4) { // First case: take highest-pT
          fsr=maxPTg;
          lepWithFsr=maxPTgLep;
        } else {
          fsr=minDRg;
          lepWithFsr=minDRgLep;
        }
      }

      myCand.addUserFloat("dauWithFSR",lepWithFsr); // Index of the cand daughter with associated FSR photon //FIXME must be removed

      if (fsr!=0) {
        // Add daughter and set p4.
        myCand.addUserFloat("mll",myCand.mass()); // for debug purposes
        myCand.setP4(myCand.p4()+fsr->p4());
        //      myCand.addDaughter(reco::ShallowCloneCandidate(fsr->masterClone()),"FSR"); //FIXME: fsr does not have a masterClone
        pat::PFParticle myFsr(*fsr);
        myFsr.setPdgId(22); // Fix: photons that are isFromMu have abs(pdgId)=13!!!
        myFsr.addUserFloat("leptIdx",lepWithFsr);
        myCand.addDaughter(myFsr,"FSR");
      }

    } else if (FSRMode==2) { // Run II
      float mll = myCand.mass(); // pre-FSR mass
      for (int dauIdx=0; dauIdx<2; ++dauIdx) { 
        const Candidate* d = myCand.daughter(dauIdx);
        const PhotonPtrVector* gammas = userdatahelpers::getUserPhotons(d);
        if (gammas==0) continue;
        assert(gammas->size()<=1); // Must have already been preselected.
        if (gammas->size()==1){
          const pat::PFParticle* fsr = gammas->begin()->get();
          pat::PFParticle myFsr(*fsr);
          myCand.setP4(myCand.p4()+fsr->p4());
          myFsr.setPdgId(22);
          myFsr.addUserFloat("leptIdx",dauIdx);
          //        myFsr.addUserFloat("gRelIso",0.);
          myCand.addDaughter(myFsr);
        }
      }
      
      if (myCand.numberOfDaughters()>2) {
        myCand.addUserFloat("mll",mll);
      }

    } else if (FSRMode==0) { // no FSR
      myCand.addUserFloat("dauWithFSR",-1);
      myCand.addUserFloat("d0.combRelIsoPFFSRCorr",myCand.userFloat("d0.combRelIsoPF"));
      myCand.addUserFloat("d1.combRelIsoPFFSRCorr",myCand.userFloat("d1.combRelIsoPF"));
    }
    
    //    pat::CompositeCandidate * candPtr = &myCand;
    //    ComputeSVFit(candPtr,
    //		 METx,METy,uncorrMETx,uncorrMETy,covMET,significance,
    //		 METx_UP_JES,METy_UP_JES,METx_DOWN_JES,METy_DOWN_JES,
    //		 METx_UP_TES,METy_UP_TES,METx_DOWN_TES,METy_DOWN_TES,
    //		 METx_UP_EES,METy_UP_EES,METx_DOWN_EES,METy_DOWN_EES
    //		 );
    
    // ---------------------------------------------------------------------------
    // ---------------------------------------------------------------------------
    // ------------------------- Begin SV Fit ------------------------------------
    // ---------------------------------------------------------------------------
    // ---------------------------------------------------------------------------
    const Candidate *l1 = myCand.daughter(0);
    const Candidate *l2 = myCand.daughter(1);

    classic_svFit::MeasuredTauLepton::kDecayType l1Type = GetDecayTypeFlag (l1->pdgId());
    classic_svFit::MeasuredTauLepton::kDecayType l2Type = GetDecayTypeFlag (l2->pdgId());
    double mass1 = GetMass (l1Type, l1->mass());
    double mass2 = GetMass (l2Type, l2->mass());
    double mass1_tauUP   = -1.;
    double mass1_tauDOWN = -1.;
    double mass2_tauUP   = -1.;
    double mass2_tauDOWN = -1.;
    double mass1_eleUP   = -1.;
    double mass1_eleDOWN = -1.;
    double mass2_eleUP   = -1.;
    double mass2_eleDOWN = -1.;
   
    int decay1 = -1;
    int decay2 = -1;
    if (l1Type == classic_svFit::MeasuredTauLepton::kTauToHadDecay) decay1 = (int)(userdatahelpers::getUserFloat(l1,"decayMode"));
    if (l2Type == classic_svFit::MeasuredTauLepton::kTauToHadDecay) decay2 = (int)(userdatahelpers::getUserFloat(l2,"decayMode"));
    
    float l1_iso = -1.;
    float l2_iso = -1.;
    if (l1Type == classic_svFit::MeasuredTauLepton::kTauToHadDecay) l1_iso = userdatahelpers::getUserFloat(l1,"byDeepTau2017v2p1VSjetraw");
    if (l2Type == classic_svFit::MeasuredTauLepton::kTauToHadDecay) l2_iso = userdatahelpers::getUserFloat(l2,"byDeepTau2017v2p1VSjetraw");

    bool swi = Switch (l1Type, l1->pt(), l1_iso, l2Type, l2->pt(), l2_iso);

    //Prepare tau nominal, up and down candidates
    std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;
    std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsTauUp;
    std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsTauDown;
    std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsEleUp;
    std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsEleDown;

    // init shifted version to nominal
    TLorentzVector l1_tauUP (l1->px(), l1->py(), l1->pz(), l1->energy());
    TLorentzVector l1_tauDOWN = l1_tauUP;
    TLorentzVector l2_tauUP (l2->px(), l2->py(), l2->pz(), l2->energy());
    TLorentzVector l2_tauDOWN = l2_tauUP;
    TLorentzVector l1_eleUP (l1->px(), l1->py(), l1->pz(), l1->energy());
    TLorentzVector l1_eleDOWN = l1_eleUP;
    TLorentzVector l2_eleUP (l2->px(), l2->py(), l2->pz(), l2->energy());
    TLorentzVector l2_eleDOWN = l2_eleUP;

    bool l1shifted_tau = false;
    bool l2shifted_tau = false;
    bool l1shifted_ele = false;
    bool l2shifted_ele = false;

    if (userdatahelpers::hasUserInt(l1,"isTESShifted")) l1shifted_tau = (userdatahelpers::getUserInt(l1,"isTESShifted") == 1 ? true : false);
    if (userdatahelpers::hasUserInt(l2,"isTESShifted")) l2shifted_tau = (userdatahelpers::getUserInt(l2,"isTESShifted") == 1 ? true : false);
    if (userdatahelpers::hasUserInt(l1,"isEESShifted")) l1shifted_ele = (userdatahelpers::getUserInt(l1,"isEESShifted") == 1 ? true : false);
    if (userdatahelpers::hasUserInt(l2,"isEESShifted")) l2shifted_ele = (userdatahelpers::getUserInt(l2,"isEESShifted") == 1 ? true : false);

    if (l1shifted_tau)
    {
      float pxUp = userdatahelpers::getUserFloat(l1,"px_TauUp");
      float pyUp = userdatahelpers::getUserFloat(l1,"py_TauUp");
      float pzUp = userdatahelpers::getUserFloat(l1,"pz_TauUp");
      float eUp  = userdatahelpers::getUserFloat(l1,"e_TauUp");
      float mUp  = userdatahelpers::getUserFloat(l1,"m_TauUp");
      l1_tauUP.SetPxPyPzE (pxUp, pyUp, pzUp, eUp);

      float pxDown = userdatahelpers::getUserFloat(l1,"px_TauDown");
      float pyDown = userdatahelpers::getUserFloat(l1,"py_TauDown");
      float pzDown = userdatahelpers::getUserFloat(l1,"pz_TauDown");
      float eDown  = userdatahelpers::getUserFloat(l1,"e_TauDown");
      float mDown  = userdatahelpers::getUserFloat(l1,"m_TauDown");
      l1_tauDOWN.SetPxPyPzE (pxDown, pyDown, pzDown, eDown);

      mass1_tauUP   = GetMass (l1Type, mUp);
      mass1_tauDOWN = GetMass (l1Type, mDown);
    }

    if (l2shifted_tau)
    {
      float pxUp = userdatahelpers::getUserFloat(l2,"px_TauUp");
      float pyUp = userdatahelpers::getUserFloat(l2,"py_TauUp");
      float pzUp = userdatahelpers::getUserFloat(l2,"pz_TauUp");
      float eUp  = userdatahelpers::getUserFloat(l2,"e_TauUp");
      float mUp  = userdatahelpers::getUserFloat(l2,"m_TauUp");
      l2_tauUP.SetPxPyPzE (pxUp, pyUp, pzUp, eUp);

      float pxDown = userdatahelpers::getUserFloat(l2,"px_TauDown");
      float pyDown = userdatahelpers::getUserFloat(l2,"py_TauDown");
      float pzDown = userdatahelpers::getUserFloat(l2,"pz_TauDown");
      float eDown  = userdatahelpers::getUserFloat(l2,"e_TauDown");
      float mDown  = userdatahelpers::getUserFloat(l2,"m_TauDown");
      l2_tauDOWN.SetPxPyPzE (pxDown, pyDown, pzDown, eDown);

      mass2_tauUP   = GetMass (l2Type, mUp);
      mass2_tauDOWN = GetMass (l2Type, mDown);
    }

    if (l1shifted_ele)
    {
      float pxUp = userdatahelpers::getUserFloat(l1,"px_EleUp");
      float pyUp = userdatahelpers::getUserFloat(l1,"py_EleUp");
      float pzUp = userdatahelpers::getUserFloat(l1,"pz_EleUp");
      float eUp  = userdatahelpers::getUserFloat(l1,"e_EleUp");
      float mUp  = userdatahelpers::getUserFloat(l1,"m_EleUp");
      l1_eleUP.SetPxPyPzE (pxUp, pyUp, pzUp, eUp);

      float pxDown = userdatahelpers::getUserFloat(l1,"px_EleDown");
      float pyDown = userdatahelpers::getUserFloat(l1,"py_EleDown");
      float pzDown = userdatahelpers::getUserFloat(l1,"pz_EleDown");
      float eDown  = userdatahelpers::getUserFloat(l1,"e_EleDown");
      float mDown  = userdatahelpers::getUserFloat(l1,"m_EleDown");
      l1_eleDOWN.SetPxPyPzE (pxDown, pyDown, pzDown, eDown);

      mass1_eleUP   = GetMass (l1Type, mUp);
      mass1_eleDOWN = GetMass (l1Type, mDown);
    }

    if (l2shifted_ele)
    {
      float pxUp = userdatahelpers::getUserFloat(l2,"px_EleUp");
      float pyUp = userdatahelpers::getUserFloat(l2,"py_EleUp");
      float pzUp = userdatahelpers::getUserFloat(l2,"pz_EleUp");
      float eUp  = userdatahelpers::getUserFloat(l2,"e_EleUp");
      float mUp  = userdatahelpers::getUserFloat(l2,"m_EleUp");
      l2_eleUP.SetPxPyPzE (pxUp, pyUp, pzUp, eUp);

      float pxDown = userdatahelpers::getUserFloat(l2,"px_EleDown");
      float pyDown = userdatahelpers::getUserFloat(l2,"py_EleDown");
      float pzDown = userdatahelpers::getUserFloat(l2,"pz_EleDown");
      float eDown  = userdatahelpers::getUserFloat(l2,"e_EleDown");
      float mDown  = userdatahelpers::getUserFloat(l2,"m_EleDown");
      l2_eleDOWN.SetPxPyPzE (pxDown, pyDown, pzDown, eDown);

      mass2_eleUP   = GetMass (l2Type, mUp);
      mass2_eleDOWN = GetMass (l2Type, mDown);
    }

    // set lepton vector, ordered for SVfit
    if (swi)  // 2 first, 1 second (switch)
    {
      measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(l2Type, l2->pt(), l2->eta(), l2->phi(), mass2, decay2 ));
      measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(l1Type, l1->pt(), l1->eta(), l1->phi(), mass1, decay1 ));
      
      if (l1shifted_tau && l2shifted_tau){
      measuredTauLeptonsTauUp.push_back(classic_svFit::MeasuredTauLepton(l2Type, l2_tauUP.Pt(), l2_tauUP.Eta(), l2_tauUP.Phi(), mass2_tauUP, decay2 ));
      measuredTauLeptonsTauUp.push_back(classic_svFit::MeasuredTauLepton(l1Type, l1_tauUP.Pt(), l1_tauUP.Eta(), l1_tauUP.Phi(), mass1_tauUP, decay1 ));

      measuredTauLeptonsTauDown.push_back(classic_svFit::MeasuredTauLepton(l2Type, l2_tauDOWN.Pt(), l2_tauDOWN.Eta(), l2_tauDOWN.Phi(), mass2_tauDOWN, decay2 ));
      measuredTauLeptonsTauDown.push_back(classic_svFit::MeasuredTauLepton(l1Type, l1_tauDOWN.Pt(), l1_tauDOWN.Eta(), l1_tauDOWN.Phi(), mass1_tauDOWN, decay1 ));
      }

      if (l1shifted_ele && l2shifted_ele){
      measuredTauLeptonsEleUp.push_back(classic_svFit::MeasuredTauLepton(l2Type, l2_eleUP.Pt(), l2_eleUP.Eta(), l2_eleUP.Phi(), mass2_eleUP, decay2 ));
      measuredTauLeptonsEleUp.push_back(classic_svFit::MeasuredTauLepton(l1Type, l1_eleUP.Pt(), l1_eleUP.Eta(), l1_eleUP.Phi(), mass1_eleUP, decay1 ));

      measuredTauLeptonsEleDown.push_back(classic_svFit::MeasuredTauLepton(l2Type, l2_eleDOWN.Pt(), l2_eleDOWN.Eta(), l2_eleDOWN.Phi(), mass2_eleDOWN, decay2 ));
      measuredTauLeptonsEleDown.push_back(classic_svFit::MeasuredTauLepton(l1Type, l1_eleDOWN.Pt(), l1_eleDOWN.Eta(), l1_eleDOWN.Phi(), mass1_eleDOWN, decay1 ));
      }
    }

    else
    {
      measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(l1Type, l1->pt(), l1->eta(), l1->phi(), mass1, decay1 ));
      measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(l2Type, l2->pt(), l2->eta(), l2->phi(), mass2, decay2 ));

      if (l1shifted_tau && l2shifted_tau){
      measuredTauLeptonsTauUp.push_back(classic_svFit::MeasuredTauLepton(l1Type, l1_tauUP.Pt(), l1_tauUP.Eta(), l1_tauUP.Phi(), mass1_tauUP, decay1 ));
      measuredTauLeptonsTauUp.push_back(classic_svFit::MeasuredTauLepton(l2Type, l2_tauUP.Pt(), l2_tauUP.Eta(), l2_tauUP.Phi(), mass2_tauUP, decay2 ));

      measuredTauLeptonsTauDown.push_back(classic_svFit::MeasuredTauLepton(l1Type, l1_tauDOWN.Pt(), l1_tauDOWN.Eta(), l1_tauDOWN.Phi(), mass1_tauDOWN, decay1 ));
      measuredTauLeptonsTauDown.push_back(classic_svFit::MeasuredTauLepton(l2Type, l2_tauDOWN.Pt(), l2_tauDOWN.Eta(), l2_tauDOWN.Phi(), mass2_tauDOWN, decay2 ));
      }

      if (l1shifted_ele && l2shifted_ele){
      measuredTauLeptonsEleUp.push_back(classic_svFit::MeasuredTauLepton(l1Type, l1_eleUP.Pt(), l1_eleUP.Eta(), l1_eleUP.Phi(), mass1_eleUP, decay1 ));
      measuredTauLeptonsEleUp.push_back(classic_svFit::MeasuredTauLepton(l2Type, l2_eleUP.Pt(), l2_eleUP.Eta(), l2_eleUP.Phi(), mass2_eleUP, decay2 ));

      measuredTauLeptonsEleDown.push_back(classic_svFit::MeasuredTauLepton(l1Type, l1_eleDOWN.Pt(), l1_eleDOWN.Eta(), l1_eleDOWN.Phi(), mass1_eleDOWN, decay1 ));
      measuredTauLeptonsEleDown.push_back(classic_svFit::MeasuredTauLepton(l2Type, l2_eleDOWN.Pt(), l2_eleDOWN.Eta(), l2_eleDOWN.Phi(), mass2_eleDOWN, decay2 ));
      }
    }


    // define algorithm (set the debug level to 3 for testing)
    unsigned int verbosity = 0;
    
    bool ComputeSV = false;

    double SVfitMass = -999.;
    double SVfitMassTauUp = -999.;
    double SVfitMassTauDown = -999.;
    double SVfitMassMETUp = -999.;
    double SVfitMassMETDown = -999.;
    double SVfitMassEleUp = -999.;
    double SVfitMassEleDown = -999.;

    double SVfitTransverseMass = -999.;
    double SVfitTransverseMassTauUp = -999.;
    double SVfitTransverseMassTauDown = -999.;
    double SVfitTransverseMassMETUp = -999.;
    double SVfitTransverseMassMETDown = -999.;
    double SVfitTransverseMassEleUp = -999.;
    double SVfitTransverseMassEleDown = -999.;

    double SVpt = -999.;
    double SVptTauUp = -999.;
    double SVptTauDown = -999.;
    double SVptMETUp = -999.;
    double SVptMETDown = -999.;
    double SVptEleUp = -999.;
    double SVptEleDown = -999.;

    double SVeta = -999.;
    double SVetaTauUp = -999.;
    double SVetaTauDown = -999.;
    double SVetaMETUp = -999.;
    double SVetaMETDown = -999.;
    double SVetaEleUp = -999.;
    double SVetaEleDown = -999.;

    double SVphi = -999.;
    double SVphiTauUp = -999.;
    double SVphiTauDown = -999.;
    double SVphiMETUp = -999.;
    double SVphiMETDown = -999.;
    double SVphiEleUp = -999.;
    double SVphiEleDown = -999.;

    double SVfitMassUnc = -999.;
    double SVfitMassUncTauUp = -999.;
    double SVfitMassUncTauDown = -999.;
    double SVfitMassUncMETUp = -999.;
    double SVfitMassUncMETDown = -999.;
    double SVfitMassUncEleUp = -999.;
    double SVfitMassUncEleDown = -999.;

    double SVfitTransverseMassUnc = -999.;
    double SVfitTransverseMassUncTauUp = -999.;
    double SVfitTransverseMassUncTauDown = -999.;
    double SVfitTransverseMassUncMETUp = -999.;
    double SVfitTransverseMassUncMETDown = -999.;
    double SVfitTransverseMassUncEleUp = -999.;
    double SVfitTransverseMassUncEleDown = -999.;

    double SVptUnc = -999.;
    double SVptUncTauUp = -999.;
    double SVptUncTauDown = -999.;
    double SVptUncMETUp = -999.;
    double SVptUncMETDown = -999.;
    double SVptUncEleUp = -999.;
    double SVptUncEleDown = -999.;

    double SVetaUnc = -999.;
    double SVetaUncTauUp = -999.;
    double SVetaUncTauDown = -999.;
    double SVetaUncMETUp = -999.;
    double SVetaUncMETDown = -999.;
    double SVetaUncEleUp = -999.;
    double SVetaUncEleDown = -999.;

    double SVphiUnc = -999.;
    double SVphiUncTauUp = -999.;
    double SVphiUncTauDown = -999.;
    double SVphiUncMETUp = -999.;
    double SVphiUncMETDown = -999.;
    double SVphiUncEleUp = -999.;
    double SVphiUncEleDown = -999.;

    double SVMETRho = -999.;        // fitted MET
    double SVMETRhoTauUp = -999.;   // fitted MET
    double SVMETRhoTauDown = -999.; // fitted MET
    double SVMETRhoMETUp = -999.;   // fitted MET
    double SVMETRhoMETDown = -999.; // fitted MET
    double SVMETRhoEleUp = -999.;   // fitted MET
    double SVMETRhoEleDown = -999.; // fitted MET

    double SVMETPhi = -999.;        // fitted MET
    double SVMETPhiTauUp = -999.;   // fitted MET
    double SVMETPhiTauDown = -999.; // fitted MET
    double SVMETPhiMETUp = -999.;   // fitted MET
    double SVMETPhiMETDown = -999.; // fitted MET
    double SVMETPhiEleUp = -999.;   // fitted MET
    double SVMETPhiEleDown = -999.; // fitted MET

    // assessing pair type
    int apdg1 = abs(l1->pdgId());
    int apdg2 = abs(l2->pdgId());

    int nmu = 0;
    int nele = 0;
    int ntau = 0;

    if (apdg1 == 13) nmu++;
    if (apdg1 == 11) nele++;
    if (apdg1 == 15) ntau++;

    if (apdg2 == 13) nmu++;
    if (apdg2 == 11) nele++;
    if (apdg2 == 15) ntau++;

    pairType pType = kOther;
    if (nmu == 1 && nele == 0 && ntau == 1) pType = kMuHad;
    if (nmu == 0 && nele == 1 && ntau == 1) pType = kEHad;
    if (nmu == 0 && nele == 0 && ntau == 2) pType = kHadHad;
    if (nmu == 2 && nele == 0 && ntau == 0) pType = kMuMu;
    if (nmu == 0 && nele == 2 && ntau == 0) pType = kEE;
    if (nmu == 1 && nele == 1 && ntau == 0) pType = kEMu;

    // Define the k factor
    double kappa; // use 3 for emu, 4 for etau and mutau, 5 for tautau channel
    if      (pType == kMuHad ) kappa = 4.;  // mutau
    else if (pType == kEHad  ) kappa = 4.;  // etau
    else if (pType == kHadHad) kappa = 5.;  // tautau
    else                       kappa = 3.;  // ee, emu, mumu
    

    // only run SVfit if taus are passing discriminator, skip mumu and ee pairs, apply very loose quality cuts on objects
    if (IsInteresting(l1, l2))
    {
      ClassicSVfit algo(verbosity);
      algo.addLogM_fixed(false, kappa);
      algo.addLogM_dynamic(false);
      algo.integrate(measuredTauLeptons, METx, METy, covMET);

      if ( algo.isValidSolution() )
      {
        ComputeSV	       = true;
	SVfitMass              = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getMass(); 
        SVfitTransverseMass    = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getTransverseMass();
        SVpt                   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getPt();
        SVeta                  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getEta();
        SVphi                  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getPhi();
        SVfitMassUnc           = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getMassErr();
        SVfitTransverseMassUnc = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getTransverseMassErr();
        SVptUnc                = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getPtErr();
        SVetaUnc               = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getEtaErr();
        SVphiUnc               = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getPhiErr();
     
        cout << "--- SVFit Output Debug ---" << endl;
        cout << "SVfitMass           = " << SVfitMass << endl;
        //cout << "SVfitTransverseMass = " << SVfitTransverseMass << endl;
        //cout << "SVpt 	             = " << SVpt << endl;
        //cout << "SVeta	             = " << SVeta << endl;
        //cout << "SVphi	             = " << SVphi << endl;
        
        ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> measuredTau1(l1->pt(), l1->eta(), l1->phi(), mass1);
        ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> measuredTau2(l2->pt(), l2->eta(), l2->phi(), mass2);
        ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> measuredDiTauSystem = measuredTau1 + measuredTau2;
        ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> fittedDiTauSystem(SVpt, SVeta, SVphi, SVfitMass);
        Vector fittedMET = (fittedDiTauSystem.Vect() - measuredDiTauSystem.Vect());
        SVMETRho = fittedMET.Rho();
        SVMETPhi = fittedMET.Phi();
      }
      else
        SVfitMass = -111; // -111: SVfit failed (cfr: -999: SVfit not computed)

      // ------------------------------------------------------------------
      // compute up/down SVfit variations
      // ------------------------------------------------------------------
      if ( (l1shifted_tau || l2shifted_tau) && _computeForUpDownTES)
      {

        // UP
        ClassicSVfit algoTauUp(verbosity);
        algoTauUp.addLogM_fixed(false, kappa);
        algoTauUp.addLogM_dynamic(false);
        algoTauUp.integrate(measuredTauLeptonsTauUp, METx_UP_TES, METy_UP_TES, covMET);
        
        if ( algoTauUp.isValidSolution() )
        {
          SVfitMassTauUp              = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoTauUp.getHistogramAdapter())->getMass();
          SVfitTransverseMassTauUp    = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoTauUp.getHistogramAdapter())->getTransverseMass();
          SVptTauUp                   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoTauUp.getHistogramAdapter())->getPt();
          SVetaTauUp                  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoTauUp.getHistogramAdapter())->getEta();
          SVphiTauUp                  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoTauUp.getHistogramAdapter())->getPhi();
          SVfitMassUncTauUp           = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoTauUp.getHistogramAdapter())->getMassErr();
          SVfitTransverseMassUncTauUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoTauUp.getHistogramAdapter())->getTransverseMassErr();
          SVptUncTauUp                = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoTauUp.getHistogramAdapter())->getPtErr();
          SVetaUncTauUp               = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoTauUp.getHistogramAdapter())->getEtaErr();
          SVphiUncTauUp               = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoTauUp.getHistogramAdapter())->getPhiErr();
          
          ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> measuredTau1Up(l1_tauUP.Pt(), l1_tauUP.Eta(), l1_tauUP.Phi(), mass1_tauUP);
          ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> measuredTau2Up(l2_tauUP.Pt(), l2_tauUP.Eta(), l2_tauUP.Phi(), mass2_tauUP);
          ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> measuredDiTauSystemUp = measuredTau1Up + measuredTau2Up;
          ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> fittedDiTauSystemUp(SVptTauUp, SVetaTauUp, SVphiTauUp, SVfitMassTauUp);
          Vector fittedMETUp = (fittedDiTauSystemUp.Vect() - measuredDiTauSystemUp.Vect());
          SVMETRhoTauUp = fittedMETUp.Rho();
          SVMETPhiTauUp = fittedMETUp.Phi();
        }
        else
          SVfitMassTauUp = -111; // -111: SVfit failed (cfr: -999: SVfit not computed)

        // DOWN
        ClassicSVfit algoTauDown(verbosity);
        algoTauDown.addLogM_fixed(false, kappa);
        algoTauDown.addLogM_dynamic(false);
        algoTauDown.integrate(measuredTauLeptonsTauDown, METx_DOWN_TES, METy_DOWN_TES, covMET);

        if ( algoTauDown.isValidSolution() )
        {
        
          SVfitMassTauDown              = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoTauDown.getHistogramAdapter())->getMass();
          SVfitTransverseMassTauDown    = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoTauDown.getHistogramAdapter())->getTransverseMass();
          SVptTauDown                   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoTauDown.getHistogramAdapter())->getPt();
          SVetaTauDown                  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoTauDown.getHistogramAdapter())->getEta();
          SVphiTauDown                  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoTauDown.getHistogramAdapter())->getPhi();
          SVfitMassUncTauDown           = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoTauDown.getHistogramAdapter())->getMassErr();
          SVfitTransverseMassUncTauDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoTauDown.getHistogramAdapter())->getTransverseMassErr();
          SVptUncTauDown                = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoTauDown.getHistogramAdapter())->getPtErr();
          SVetaUncTauDown               = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoTauDown.getHistogramAdapter())->getEtaErr();
          SVphiUncTauDown               = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoTauDown.getHistogramAdapter())->getPhiErr();

          ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> measuredTau1Down(l1_tauDOWN.Pt(), l1_tauDOWN.Eta(), l1_tauDOWN.Phi(), mass1_tauDOWN);
          ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> measuredTau2Down(l2_tauDOWN.Pt(), l2_tauDOWN.Eta(), l2_tauDOWN.Phi(), mass2_tauDOWN);
          ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> measuredDiTauSystemDown = measuredTau1Down + measuredTau2Down;
          ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> fittedDiTauSystemDown(SVptTauDown, SVetaTauDown, SVphiTauDown, SVfitMassTauDown);
          Vector fittedMETDown = (fittedDiTauSystemDown.Vect() - measuredDiTauSystemDown.Vect());
          SVMETRhoTauDown = fittedMETDown.Rho();
          SVMETPhiTauDown = fittedMETDown.Phi();
        }
        else
          SVfitMassTauDown = -111; // -111: SVfit failed (cfr: -999: SVfit not computed)

      }
      else if (_computeForUpDownTES) // if I asked to have UP/DOWN variation, but this pair has not tau shifted, simply put central value. 
      {                              // instead, if I dindn't ask for up/down, I get -999 everywhere to remember my mistakes
          SVfitMassTauUp = SVfitMassTauDown = SVfitMass ;
          SVfitTransverseMassTauUp = SVfitTransverseMassTauDown = SVfitTransverseMass ;
          SVptTauUp = SVptTauDown = SVpt ;
          SVetaTauUp = SVetaTauDown = SVeta ;
          SVphiTauUp = SVphiTauDown = SVphi ;
          SVfitMassUncTauUp = SVfitMassUncTauDown = SVfitMassUnc ;
          SVfitTransverseMassUncTauUp = SVfitTransverseMassUncTauDown = SVfitTransverseMassUnc ;
          SVptUncTauUp = SVptUncTauDown = SVptUnc ;
          SVetaUncTauUp = SVetaUncTauDown = SVetaUnc ;
          SVphiUncTauUp = SVphiUncTauDown = SVphiUnc ;
          SVMETRhoTauUp = SVMETRhoTauDown = SVMETRho ;
          SVMETPhiTauUp = SVMETPhiTauDown = SVMETPhi ;
      }

      // ------------------------------------------------------------------
      // compute up/down SVfit variations due to e->tau ES
      // ------------------------------------------------------------------
      if ( (l1shifted_ele || l2shifted_ele) && _computeForUpDownTES)
      {
        // UP E->tau ES
        ClassicSVfit algoEleUp(verbosity);
        algoEleUp.addLogM_fixed(false, kappa);
        algoEleUp.addLogM_dynamic(false);
        algoEleUp.integrate(measuredTauLeptonsEleUp, METx_UP_EES, METy_UP_EES, covMET);

        if ( algoEleUp.isValidSolution() )
        {
          SVfitMassEleUp              = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoEleUp.getHistogramAdapter())->getMass();
          SVfitTransverseMassEleUp    = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoEleUp.getHistogramAdapter())->getTransverseMass();
          SVptEleUp                   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoEleUp.getHistogramAdapter())->getPt();
          SVetaEleUp                  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoEleUp.getHistogramAdapter())->getEta();
          SVphiEleUp                  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoEleUp.getHistogramAdapter())->getPhi();
          SVfitMassUncEleUp           = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoEleUp.getHistogramAdapter())->getMassErr();
          SVfitTransverseMassUncEleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoEleUp.getHistogramAdapter())->getTransverseMassErr();
          SVptUncEleUp                = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoEleUp.getHistogramAdapter())->getPtErr();
          SVetaUncEleUp               = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoEleUp.getHistogramAdapter())->getEtaErr();
          SVphiUncEleUp               = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoEleUp.getHistogramAdapter())->getPhiErr();

          ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> measuredTau1Up(l1_eleUP.Pt(), l1_eleUP.Eta(), l1_eleUP.Phi(), mass1_eleUP);
          ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> measuredTau2Up(l2_eleUP.Pt(), l2_eleUP.Eta(), l2_eleUP.Phi(), mass2_eleUP);
          ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> measuredDiTauSystemUp = measuredTau1Up + measuredTau2Up;
          ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> fittedDiTauSystemUp(SVptEleUp, SVetaEleUp, SVphiEleUp, SVfitMassEleUp);
          Vector fittedMETUp = (fittedDiTauSystemUp.Vect() - measuredDiTauSystemUp.Vect());
          SVMETRhoEleUp = fittedMETUp.Rho();
          SVMETPhiEleUp = fittedMETUp.Phi();
        }
        else
          SVfitMassEleUp = -111; // -111: SVfit failed (cfr: -999: SVfit not computed)

        // DOWN E->tau ES
        ClassicSVfit algoEleDown(verbosity);
        algoEleDown.addLogM_fixed(false, kappa);
        algoEleDown.addLogM_dynamic(false);
        algoEleDown.integrate(measuredTauLeptonsEleDown, METx_DOWN_EES, METy_DOWN_EES, covMET);

        if ( algoEleDown.isValidSolution() )
        {
          SVfitMassEleDown              = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoEleDown.getHistogramAdapter())->getMass();
          SVfitTransverseMassEleDown    = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoEleDown.getHistogramAdapter())->getTransverseMass();
          SVptEleDown                   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoEleDown.getHistogramAdapter())->getPt();
          SVetaEleDown                  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoEleDown.getHistogramAdapter())->getEta();
          SVphiEleDown                  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoEleDown.getHistogramAdapter())->getPhi();
          SVfitMassUncEleDown           = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoEleDown.getHistogramAdapter())->getMassErr();
          SVfitTransverseMassUncEleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoEleDown.getHistogramAdapter())->getTransverseMassErr();
          SVptUncEleDown                = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoEleDown.getHistogramAdapter())->getPtErr();
          SVetaUncEleDown               = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoEleDown.getHistogramAdapter())->getEtaErr();
          SVphiUncEleDown               = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoEleDown.getHistogramAdapter())->getPhiErr();

          ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> measuredTau1Down(l1_eleDOWN.Pt(), l1_eleDOWN.Eta(), l1_eleDOWN.Phi(), mass1_eleDOWN);
          ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> measuredTau2Down(l2_eleDOWN.Pt(), l2_eleDOWN.Eta(), l2_eleDOWN.Phi(), mass2_eleDOWN);
          ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> measuredDiTauSystemDown = measuredTau1Down + measuredTau2Down;
          ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> fittedDiTauSystemDown(SVptEleDown, SVetaEleDown, SVphiEleDown, SVfitMassEleDown);
          Vector fittedMETDown = (fittedDiTauSystemDown.Vect() - measuredDiTauSystemDown.Vect());
          SVMETRhoEleDown = fittedMETDown.Rho();
          SVMETPhiEleDown = fittedMETDown.Phi();
        }
        else
          SVfitMassEleDown = -111; // -111: SVfit failed (cfr: -999: SVfit not computed)

      }
      else if (_computeForUpDownTES) // if I asked to have UP/DOWN variation, but this pair has not tau shifted, simply put central value.
      {                              // instead, if I dindn't ask for up/down, I get -999 everywhere to remember my mistakes
          SVfitMassEleUp = SVfitMassEleDown = SVfitMass ;
          SVfitTransverseMassEleUp = SVfitTransverseMassEleDown = SVfitTransverseMass ;
          SVptEleUp = SVptEleDown = SVpt ;
          SVetaEleUp = SVetaEleDown = SVeta ;
          SVphiEleUp = SVphiEleDown = SVphi ;
          SVfitMassUncEleUp = SVfitMassUncEleDown = SVfitMassUnc ;
          SVfitTransverseMassUncEleUp = SVfitTransverseMassUncEleDown = SVfitTransverseMassUnc ;
          SVptUncEleUp = SVptUncEleDown = SVptUnc ;
          SVetaUncEleUp = SVetaUncEleDown = SVetaUnc ;
          SVphiUncEleUp = SVphiUncEleDown = SVphiUnc ;
          SVMETRhoEleUp = SVMETRhoEleDown = SVMETRho ;
          SVMETPhiEleUp = SVMETPhiEleDown = SVMETPhi ;
      }

      // ----------------------------------------------------------------
      // compute UP/DOWN due to the MET JES shift
      // ----------------------------------------------------------------
      if (_computeForUpDownMET)
      {
        // UP MET
        ClassicSVfit algoMETUp(verbosity);
        algoMETUp.addLogM_fixed(false, kappa);
        algoMETUp.addLogM_dynamic(false);
        //algoTauUp.shiftVisPt(true, inputFile_visPtResolution_); //not in Classic_svFit
        algoMETUp.integrate(measuredTauLeptons, METx_UP_JES, METy_UP_JES, covMET);

        if ( algoMETUp.isValidSolution() )
        {
          SVfitMassMETUp              = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoMETUp.getHistogramAdapter())->getMass();
          SVfitTransverseMassMETUp    = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoMETUp.getHistogramAdapter())->getTransverseMass();
          SVptMETUp                   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoMETUp.getHistogramAdapter())->getPt();
          SVetaMETUp                  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoMETUp.getHistogramAdapter())->getEta();
          SVphiMETUp                  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoMETUp.getHistogramAdapter())->getPhi();
          SVfitMassUncMETUp           = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoMETUp.getHistogramAdapter())->getMassErr();
          SVfitTransverseMassUncMETUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoMETUp.getHistogramAdapter())->getTransverseMassErr();
          SVptUncMETUp                = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoMETUp.getHistogramAdapter())->getPtErr();
          SVetaUncMETUp               = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoMETUp.getHistogramAdapter())->getEtaErr();
          SVphiUncMETUp               = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoMETUp.getHistogramAdapter())->getPhiErr();

          ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> measuredTau1(l1->pt(), l1->eta(), l1->phi(), mass1);
          ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> measuredTau2(l2->pt(), l2->eta(), l2->phi(), mass2);
          ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> measuredDiTauSystem = measuredTau1 + measuredTau2;
          ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> fittedDiTauSystem(SVptMETUp, SVetaMETUp, SVphiMETUp, SVfitMassMETUp);
          Vector fittedMET = (fittedDiTauSystem.Vect() - measuredDiTauSystem.Vect());
          SVMETRhoMETUp = fittedMET.Rho();
          SVMETPhiMETUp = fittedMET.Phi();
        }
        else
          SVfitMassMETUp = -111; // -111: SVfit failed (cfr: -999: SVfit not computed)

        // DOWN MET
        ClassicSVfit algoMETDown(verbosity);
        algoMETDown.addLogM_fixed(false, kappa);
        algoMETDown.addLogM_dynamic(false);
        //algoTauUp.shiftVisPt(true, inputFile_visPtResolution_); //not in Classic_svFit
        algoMETDown.integrate(measuredTauLeptons, METx_DOWN_JES, METy_DOWN_JES, covMET);

        if ( algoMETDown.isValidSolution() )
        {
          SVfitMassMETDown              = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoMETDown.getHistogramAdapter())->getMass();
          SVfitTransverseMassMETDown    = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoMETDown.getHistogramAdapter())->getTransverseMass();
          SVptMETDown                   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoMETDown.getHistogramAdapter())->getPt();
          SVetaMETDown                  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoMETDown.getHistogramAdapter())->getEta();
          SVphiMETDown                  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoMETDown.getHistogramAdapter())->getPhi();
          SVfitMassUncMETDown           = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoMETDown.getHistogramAdapter())->getMassErr();
          SVfitTransverseMassUncMETDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoMETDown.getHistogramAdapter())->getTransverseMassErr();
          SVptUncMETDown                = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoMETDown.getHistogramAdapter())->getPtErr();
          SVetaUncMETDown               = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoMETDown.getHistogramAdapter())->getEtaErr();
          SVphiUncMETDown               = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algoMETDown.getHistogramAdapter())->getPhiErr();

          ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> measuredTau1(l1->pt(), l1->eta(), l1->phi(), mass1);
          ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> measuredTau2(l2->pt(), l2->eta(), l2->phi(), mass2);
          ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> measuredDiTauSystem = measuredTau1 + measuredTau2;
          ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> fittedDiTauSystem(SVptMETDown, SVetaMETDown, SVphiMETDown, SVfitMassMETDown);
          Vector fittedMET = (fittedDiTauSystem.Vect() - measuredDiTauSystem.Vect());
          SVMETRhoMETDown = fittedMET.Rho();
          SVMETPhiMETDown = fittedMET.Phi();
        }
        else
          SVfitMassMETDown = -111; // -111: SVfit failed (cfr: -999: SVfit not computed)
      }

    } // end of quality checks IF

    // add user floats: SVfit mass, met properties, etc..  
    myCand.addUserFloat("ComputeSV", ComputeSV);
    myCand.addUserFloat("SVfitMass", (float) SVfitMass);
    myCand.addUserFloat("SVfitTransverseMass", (float) SVfitTransverseMass);
    myCand.addUserFloat("SVfit_pt", (float) SVpt);
    myCand.addUserFloat("SVfit_eta", (float) SVeta);
    myCand.addUserFloat("SVfit_phi", (float) SVphi);
    myCand.addUserFloat("SVfitMassUnc", (float) SVfitMassUnc);
    myCand.addUserFloat("SVfitTransverseMassUnc", (float) SVfitTransverseMassUnc);
    myCand.addUserFloat("SVfit_ptUnc", (float) SVptUnc);
    myCand.addUserFloat("SVfit_etaUnc", (float) SVetaUnc);
    myCand.addUserFloat("SVfit_phiUnc", (float) SVphiUnc);
    myCand.addUserFloat("SVfit_METRho", (float) SVMETRho);
    myCand.addUserFloat("SVfit_METPhi", (float) SVMETPhi);

    if (_computeForUpDownTES)
    {
      myCand.addUserFloat("SVfitMassTauUp", (float) SVfitMassTauUp);
      myCand.addUserFloat("SVfitMassTauDown", (float) SVfitMassTauDown);
      myCand.addUserFloat("SVfitTransverseMassTauUp", (float) SVfitTransverseMassTauUp);
      myCand.addUserFloat("SVfitTransverseMassTauDown", (float) SVfitTransverseMassTauDown);
      myCand.addUserFloat("SVfit_ptTauUp", (float) SVptTauUp);
      myCand.addUserFloat("SVfit_ptTauDown", (float) SVptTauDown);
      myCand.addUserFloat("SVfit_etaTauUp", (float) SVetaTauUp);
      myCand.addUserFloat("SVfit_etaTauDown", (float) SVetaTauDown);
      myCand.addUserFloat("SVfit_phiTauUp", (float) SVphiTauUp);
      myCand.addUserFloat("SVfit_phiTauDown", (float) SVphiTauDown);
      myCand.addUserFloat("SVfitMassUncTauUp", (float) SVfitMassUncTauUp);
      myCand.addUserFloat("SVfitMassUncTauDown", (float) SVfitMassUncTauDown);
      myCand.addUserFloat("SVfitTransverseMassUncTauUp", (float) SVfitTransverseMassUncTauUp);
      myCand.addUserFloat("SVfitTransverseMassUncTauDown", (float) SVfitTransverseMassUncTauDown);
      myCand.addUserFloat("SVfit_ptUncTauUp", (float) SVptUncTauUp);
      myCand.addUserFloat("SVfit_ptUncTauDown", (float) SVptUncTauDown);
      myCand.addUserFloat("SVfit_etaUncTauUp", (float) SVetaUncTauUp);
      myCand.addUserFloat("SVfit_etaUncTauDown", (float) SVetaUncTauDown);
      myCand.addUserFloat("SVfit_phiUncTauUp", (float) SVphiUncTauUp);
      myCand.addUserFloat("SVfit_phiUncTauDown", (float) SVphiUncTauDown);
      myCand.addUserFloat("SVfit_METRhoTauUp", (float) SVMETRhoTauUp);
      myCand.addUserFloat("SVfit_METRhoTauDown", (float) SVMETRhoTauDown);
      myCand.addUserFloat("SVfit_METPhiTauUp", (float) SVMETPhiTauUp);
      myCand.addUserFloat("SVfit_METPhiTauDown", (float) SVMETPhiTauDown);

      myCand.addUserFloat("SVfitMassEleUp", (float) SVfitMassEleUp);
      myCand.addUserFloat("SVfitMassEleDown", (float) SVfitMassEleDown);
      myCand.addUserFloat("SVfitTransverseMassEleUp", (float) SVfitTransverseMassEleUp);
      myCand.addUserFloat("SVfitTransverseMassEleDown", (float) SVfitTransverseMassEleDown);
      myCand.addUserFloat("SVfit_ptEleUp", (float) SVptEleUp);
      myCand.addUserFloat("SVfit_ptEleDown", (float) SVptEleDown);
      myCand.addUserFloat("SVfit_etaEleUp", (float) SVetaEleUp);
      myCand.addUserFloat("SVfit_etaEleDown", (float) SVetaEleDown);
      myCand.addUserFloat("SVfit_phiEleUp", (float) SVphiEleUp);
      myCand.addUserFloat("SVfit_phiEleDown", (float) SVphiEleDown);
      myCand.addUserFloat("SVfitMassUncEleUp", (float) SVfitMassUncEleUp);
      myCand.addUserFloat("SVfitMassUncEleDown", (float) SVfitMassUncEleDown);
      myCand.addUserFloat("SVfitTransverseMassUncEleUp", (float) SVfitTransverseMassUncEleUp);
      myCand.addUserFloat("SVfitTransverseMassUncEleDown", (float) SVfitTransverseMassUncEleDown);
      myCand.addUserFloat("SVfit_ptUncEleUp", (float) SVptUncEleUp);
      myCand.addUserFloat("SVfit_ptUncEleDown", (float) SVptUncEleDown);
      myCand.addUserFloat("SVfit_etaUncEleUp", (float) SVetaUncEleUp);
      myCand.addUserFloat("SVfit_etaUncEleDown", (float) SVetaUncEleDown);
      myCand.addUserFloat("SVfit_phiUncEleUp", (float) SVphiUncEleUp);
      myCand.addUserFloat("SVfit_phiUncEleDown", (float) SVphiUncEleDown);
      myCand.addUserFloat("SVfit_METRhoEleUp", (float) SVMETRhoEleUp);
      myCand.addUserFloat("SVfit_METRhoEleDown", (float) SVMETRhoEleDown);
      myCand.addUserFloat("SVfit_METPhiEleUp", (float) SVMETPhiEleUp);
      myCand.addUserFloat("SVfit_METPhiEleDown", (float) SVMETPhiEleDown);
    }

    if (_computeForUpDownMET)
    {
      myCand.addUserFloat("SVfitMassMETUp", (float) SVfitMassMETUp);
      myCand.addUserFloat("SVfitMassMETDown", (float) SVfitMassMETDown);
      myCand.addUserFloat("SVfitTransverseMassMETUp", (float) SVfitTransverseMassMETUp);
      myCand.addUserFloat("SVfitTransverseMassMETDown", (float) SVfitTransverseMassMETDown);
      myCand.addUserFloat("SVfit_ptMETUp", (float) SVptMETUp);
      myCand.addUserFloat("SVfit_ptMETDown", (float) SVptMETDown);
      myCand.addUserFloat("SVfit_etaMETUp", (float) SVetaMETUp);
      myCand.addUserFloat("SVfit_etaMETDown", (float) SVetaMETDown);
      myCand.addUserFloat("SVfit_phiMETUp", (float) SVphiMETUp);
      myCand.addUserFloat("SVfit_phiMETDown", (float) SVphiMETDown);
      myCand.addUserFloat("SVfitMassUncMETUp", (float) SVfitMassUncMETUp);
      myCand.addUserFloat("SVfitMassUncMETDown", (float) SVfitMassUncMETDown);
      myCand.addUserFloat("SVfitTransverseMassUncMETUp", (float) SVfitTransverseMassUncMETUp);
      myCand.addUserFloat("SVfitTransverseMassUncMETDown", (float) SVfitTransverseMassUncMETDown);
      myCand.addUserFloat("SVfit_ptUncMETUp", (float) SVptUncMETUp);
      myCand.addUserFloat("SVfit_ptUncMETDown", (float) SVptUncMETDown);
      myCand.addUserFloat("SVfit_etaUncMETUp", (float) SVetaUncMETUp);
      myCand.addUserFloat("SVfit_etaUncMETDown", (float) SVetaUncMETDown);
      myCand.addUserFloat("SVfit_phiUncMETUp", (float) SVphiUncMETUp);
      myCand.addUserFloat("SVfit_phiUncMETDown", (float) SVphiUncMETDown);
      myCand.addUserFloat("SVfit_METRhoMETUp", (float) SVMETRhoMETUp);
      myCand.addUserFloat("SVfit_METRhoMETDown", (float) SVMETRhoMETDown);
      myCand.addUserFloat("SVfit_METPhiMETUp", (float) SVMETPhiMETUp);
      myCand.addUserFloat("SVfit_METPhiMETDown", (float) SVMETPhiMETDown);
    }

    myCand.addUserFloat("MEt_px", (float) METx);
    myCand.addUserFloat("MEt_py", (float) METy);
    myCand.addUserFloat("uncorrMEt_px", (float) uncorrMETx);
    myCand.addUserFloat("uncorrMEt_py", (float) uncorrMETy);
    myCand.addUserFloat("MEt_cov00", (float) covMET[0][0]);
    myCand.addUserFloat("MEt_cov01", (float) covMET[0][1]);
    myCand.addUserFloat("MEt_cov10", (float) covMET[1][0]);
    myCand.addUserFloat("MEt_cov11", (float) covMET[1][1]);
    myCand.addUserFloat("MEt_significance", significance);
    myCand.addUserFloat("MEt_px_UP_JES", (float) METx_UP_JES);
    myCand.addUserFloat("MEt_py_UP_JES", (float) METy_UP_JES);
    myCand.addUserFloat("MEt_px_DOWN_JES", (float) METx_DOWN_JES);
    myCand.addUserFloat("MEt_py_DOWN_JES", (float) METy_DOWN_JES);
    myCand.addUserFloat("MEt_px_UP_TES", (float) METx_UP_TES);
    myCand.addUserFloat("MEt_py_UP_TES", (float) METy_UP_TES);
    myCand.addUserFloat("MEt_px_DOWN_TES", (float) METx_DOWN_TES);
    myCand.addUserFloat("MEt_py_DOWN_TES", (float) METy_DOWN_TES);
    myCand.addUserFloat("MEt_px_UP_EES", (float) METx_UP_EES);
    myCand.addUserFloat("MEt_py_UP_EES", (float) METy_UP_EES);
    myCand.addUserFloat("MEt_px_DOWN_EES", (float) METx_DOWN_EES);
    myCand.addUserFloat("MEt_py_DOWN_EES", (float) METy_DOWN_EES);

    // ---------------------------------------------------------------------------
    // ---------------------------------------------------------------------------
    // ------------------------- End SV Fit --------------------------------------
    // ---------------------------------------------------------------------------
    // ---------------------------------------------------------------------------


    //if (false && myCand.hasUserFloat("ComputeSV")){
      if (myCand.userFloat("ComputeSV")){
	myCand.addUserFloat("goodMass",myCand.userFloat("SVfitMass"));
      }
    //}
    else
	myCand.addUserFloat("goodMass",myCand.mass());
    
    //int id0 = myCand.daughter(0)->pdgId();
    //int id1 = myCand.daughter(1)->pdgId();

    //In terms of etau, mutau, tautau, special good and iso requirements are needed
    bool goodTau = true;
    if (abs(id0)==15 && abs(id1)==11 && !myCand.userFloat("d0.isGood_Ele"))
        goodTau=false;
    if (abs(id0)==15 && abs(id1)==13 && !myCand.userFloat("d0.isGood_Mu"))
        goodTau=false;
    if (abs(id0)==15 && abs(id1)==15 && ( !myCand.userFloat("d0.isGood_Tau") || !myCand.userFloat("d1.isGood_Tau") ))
        goodTau=false;
    if (abs(id0)==11 && abs(id1)==15 && !myCand.userFloat("d1.isGood_Ele"))
        goodTau=false;
    if (abs(id0)==13 && abs(id1)==15 && !myCand.userFloat("d1.isGood_Mu"))
        goodTau=false;
    myCand.addUserFloat("isGoodTau",goodTau);

    //HLTMatch
    //The candidate passes HLTMatch if at least one lepton daughter matches the single trigger or both matches to the double trigger
    bool muHLTMatch;
    bool eleHLTMatch;
    muHLTMatch=false;
    eleHLTMatch=false;
    if (abs(id0)==13 && abs(id1)==13) {
        if (myCand.userFloat("d0.HLTMatch1") || myCand.userFloat("d1.HLTMatch1"))
            muHLTMatch=true;
        for (size_t j=0; j<muHLTPaths2_.size(); ++j)
            if (myCand.userFloat("d0."+muHLTPaths2_[j]) && myCand.userFloat("d1."+muHLTPaths2_[j]))
                muHLTMatch=true;
    }
    if (abs(id0)==11 && abs(id1)==11) {
        if (myCand.userFloat("d0.HLTMatch1") || myCand.userFloat("d1.HLTMatch1"))
            eleHLTMatch=true;
        for (size_t j=0; j<eleHLTPaths2_.size(); ++j)
            if (myCand.userFloat("d0."+eleHLTPaths2_[j]) && myCand.userFloat("d1."+eleHLTPaths2_[j]))
                eleHLTMatch=true;
    }
    myCand.addUserFloat("muHLTMatch",muHLTMatch);
    myCand.addUserFloat("eleHLTMatch",eleHLTMatch);


    myCand.addUserFloat("OSSF",OS && SF);
    //--- Find "best Z" (closest to mZ) among those passing the "bestZAmong" selection (2011 PRL logic), now deprecated!!!
    if (preBestZSelection(myCand)) {
      float diffZmass = fabs(ZmassValue - myCand.userFloat("goodMass"));
      if (diffZmass < closestLLMassDiff) { // Best among any ll in the collection
        bestLLidx = i;
        closestLLMassDiff = diffZmass;
      }
      if (OS&&SF) {
        if (diffZmass < closestZMassDiff) { // Best among all OSSF pairs in the collection
          bestZidx = i;
          closestZMassDiff = diffZmass;
        }
	if (abs(id0) == 15 || abs(id1) == 15) {
	  if (diffZmass < closestZttMassDiff) {
	    bestZttidx = i;
	    closestZttMassDiff = diffZmass;
	  }
        } else if (abs(id0) == 13) { 
          if (diffZmass < closestZmmMassDiff) { // Best among all mu+mu- pairs in the collection
            bestZmmidx = i;
            closestZmmMassDiff = diffZmass;
          }
        } else if (abs(id0) == 11) {
          if (diffZmass < closestZeeMassDiff) { // Best among all e+e- pairs in the collection
            bestZeeidx = i;
            closestZeeMassDiff = diffZmass;
          }
        }
      }
    }
    result->push_back(myCand);
  }

  //--- Embed isBestZ flag (must be done in a separate loop)
  for (int i = 0; i< (int)result->size(); ++i) {
    pat::CompositeCandidate& myCand = (*result)[i];    
    myCand.addUserFloat("isBestZ",  (i==bestZidx));
    myCand.addUserFloat("isBestZmm",(i==bestZmmidx));
    myCand.addUserFloat("isBestZee",(i==bestZeeidx));
    myCand.addUserFloat("isBestZtt",(i==bestZttidx));
    myCand.addUserFloat("isBestInColl", (i==bestLLidx));

    //--- Embed flags (ie cuts specified in the "flags" pset)
    //    We do this here so that isBestZ is available within the cuts
    for(CutSet<pat::CompositeCandidate>::const_iterator cut = cuts.begin(); cut != cuts.end(); ++cut) {
      myCand.addUserFloat(cut->first,int((*(cut->second))(myCand)));
    }    

  }
  cout<<result->size()<<" Z candidates"<<endl;
  iEvent.put(std::move(result));
}

//void ZCandidateFiller::ComputeSVFit(pat::CompositeCandidate * candPtr,
//		    double METx, double METy, double uncorrMETx, double uncorrMETy, TMatrixD covMET, float significance,
//		    double METx_UP_JES, double METy_UP_JES, double METx_DOWN_JES, double METy_DOWN_JES,
//		    double METx_UP_TES, double METy_UP_TES, double METx_DOWN_TES, double METy_DOWN_TES,
//           	    double METx_UP_EES, double METy_UP_EES, double METx_DOWN_EES, double METy_DOWN_EES)
//{

// }


classic_svFit::MeasuredTauLepton::kDecayType ZCandidateFiller::GetDecayTypeFlag (int pdgId)
{
    if (abs(pdgId) == 11) return classic_svFit::MeasuredTauLepton::kTauToElecDecay;
    if (abs(pdgId) == 13) return classic_svFit::MeasuredTauLepton::kTauToMuDecay;
    if (abs(pdgId) == 15) return classic_svFit::MeasuredTauLepton::kTauToHadDecay;
    
    edm::LogWarning("WrongDecayModePdgID")
       << "(ClassicSVfitInterface): unable to identify decay type from pdgId"
       << "     ---> Decay will be treated as an hadronic decay";
    return classic_svFit::MeasuredTauLepton::kTauToHadDecay;
}

// decide if leptons 1 and 2 must be switched to respect SVfit conventions
bool ZCandidateFiller::Switch (classic_svFit::MeasuredTauLepton::kDecayType type1, double pt1, float l1_iso, classic_svFit::MeasuredTauLepton::kDecayType type2, double pt2, float l2_iso)
{
    // e e, mu mu, tau tau
    if (type1 == type2)
    {
      if (type1 == classic_svFit::MeasuredTauLepton::kTauToHadDecay && type2 == classic_svFit::MeasuredTauLepton::kTauToHadDecay)
        return (l1_iso < l2_iso);
      else
        return (pt1 < pt2);
    }
    
    // e tau, mu tau
    if ( (type1 == classic_svFit::MeasuredTauLepton::kTauToElecDecay || type1 == classic_svFit::MeasuredTauLepton::kTauToMuDecay) &&
         type2 == classic_svFit::MeasuredTauLepton::kTauToHadDecay ) {return false;}
    if ( (type2 == classic_svFit::MeasuredTauLepton::kTauToElecDecay || type2 == classic_svFit::MeasuredTauLepton::kTauToMuDecay) &&
         type1 == classic_svFit::MeasuredTauLepton::kTauToHadDecay ) {return true;}

    // e mu
    if (type1 == classic_svFit::MeasuredTauLepton::kTauToElecDecay && type2 == classic_svFit::MeasuredTauLepton::kTauToMuDecay) {return false;}
    if (type2 == classic_svFit::MeasuredTauLepton::kTauToElecDecay && type1 == classic_svFit::MeasuredTauLepton::kTauToMuDecay) {return true;}
    
    cout << "SVfit Standalone: ordering not done (should never happen)" << endl;
    return false;
}

// set mass (pdg ele/mu or leave cand mass for tauh
double ZCandidateFiller::GetMass (classic_svFit::MeasuredTauLepton::kDecayType type, double candMass)
{
    if (type == classic_svFit::MeasuredTauLepton::kTauToElecDecay) return 0.51100e-3;
    if (type == classic_svFit::MeasuredTauLepton::kTauToMuDecay)   return 105.658e-3;

    return candMass; // for tauh and all exceptions return cand mass
}

bool ZCandidateFiller::IsInteresting (const reco::Candidate *l1, const reco::Candidate *l2)
{
  int apdg1 = abs(l1->pdgId());
  int apdg2 = abs(l2->pdgId());

  int nmu = 0;
  int nele = 0;
  int ntau = 0;

  if (apdg1 == 13) nmu++;
  if (apdg1 == 11) nele++;
  if (apdg1 == 15) ntau++;

  if (apdg2 == 13) nmu++;
  if (apdg2 == 11) nele++;
  if (apdg2 == 15) ntau++;

  pairType pType = kOther;
  if (nmu == 1 && nele == 0 && ntau == 1) pType = kMuHad;
  if (nmu == 0 && nele == 1 && ntau == 1) pType = kEHad;
  if (nmu == 0 && nele == 0 && ntau == 2) pType = kHadHad;
  if (nmu == 2 && nele == 0 && ntau == 0) pType = kMuMu;
  if (nmu == 0 && nele == 2 && ntau == 0) pType = kEE;
  if (nmu == 1 && nele == 1 && ntau == 0) pType = kEMu;

  ///////

  // switch to apply different requirements to the objects
  if (deltaR(l1->p4(), l2->p4()) < 0.1)
  //if (deltaR(l1->p4(), l2->p4()) < 0.4)
    return false; // for overlap removal

  ///////

  // create pointers with usual pair ordering -- easier to apply cuts later
  const reco::Candidate* dau1;
  const reco::Candidate* dau2;

  if (pType == kMuHad)
  {
    dau1 = (apdg1 == 13 ? l1 : l2);
    dau2 = (apdg1 == 13 ? l2 : l1);

    //if (dau1->pt() < 17.)
    if (dau1->pt() < 10.)
      return false;

    if (dau2->pt() < 10.)
      return false;

    if (userdatahelpers::getUserInt(l2,"decayModeFindingNewDMs") != 1) // decayModeFinding == decayModeFindingOldDMs
      return false;

    //bool iso1 = (userdatahelpers::getUserFloat(l1,"combRelIsoPF") < 0.3);
    //bool iso1 = (userdatahelpers::getUserFloat(l1,"combRelIsoPF") < 0.2); //Commented during March 2020 sync: inconsistency with KLUB
    //bool iso2 = (userdatahelpers::getUserInt(l2,"byVLooseIsolationMVArun2v1DBoldDMwLT") == 1);
    //bool iso2 = (userdatahelpers::getUserInt(l2,"byVVLooseIsolationMVArun2017v2DBoldDMwLT2017") == 1); //FRA 2017
    //bool iso2 = (userdatahelpers::getUserInt(l2,"byVVVLooseDeepTau2017v2p1VSjet") == 1);
    bool iso2 = (userdatahelpers::getUserInt(l2,"byVVVLooseDeepTau2017v2p1VSjet") == 1);

    //if (!iso1 || !iso2)
    if (!iso2)
      return false;

    return true; // passed all requirements
  }

  else if (pType == kEHad)
  {
    dau1 = (apdg1 == 11 ? l1 : l2);
    dau2 = (apdg1 == 11 ? l2 : l1);

    //if (dau1->pt() < 19.)
    if (dau1->pt() < 10.)
      return false;

    if (dau2->pt() < 10.)
      return false;

    if (userdatahelpers::getUserInt(l2,"decayModeFindingNewDMs") != 1)  // decayModeFinding == decayModeFindingOldDMs
      return false;

    //bool iso1 = (userdatahelpers::getUserFloat(l1,"combRelIsoPF") < 0.3);
    //bool iso1 = (userdatahelpers::getUserFloat(l1,"combRelIsoPF") < 0.2); //Commented during March 2020 sync: inconsistency with KLUB
    //bool iso2 = (userdatahelpers::getUserInt(l2,"byVLooseIsolationMVArun2v1DBoldDMwLT") == 1);
    //bool iso2 = (userdatahelpers::getUserInt(l2,"byVVLooseIsolationMVArun2017v2DBoldDMwLT2017") == 1); //FRA 2017
    //bool iso2 = (userdatahelpers::getUserInt(l2,"byVVVLooseDeepTau2017v2p1VSjet") == 1);
    bool iso2 = (userdatahelpers::getUserInt(l2,"byVVVLooseDeepTau2017v2p1VSjet") == 1);

    //if (!iso1 || !iso2)
    if (!iso2)
      return false;

    return true; // passed all requirements
  }

  else if (pType == kHadHad)
  {
    dau1 = ((l1->pt() > l2->pt()) ? l1 : l2);
    dau2 = ((l1->pt() > l2->pt()) ? l2 : l1);

    //if (dau1->pt() < 30.)
    if (dau1->pt() < 10.)
      return false;
    
    //if (dau2->pt() < 30.)
    if (dau2->pt() < 10.)
      return false;
    
    if (userdatahelpers::getUserInt(l1,"decayModeFindingNewDMs") != 1)  // decayModeFinding == decayModeFindingOldDMs
      return false;
    
    if (userdatahelpers::getUserInt(l2,"decayModeFindingNewDMs") != 1)  // decayModeFinding == decayModeFindingOldDMs
      return false;

    //bool iso1 = (userdatahelpers::getUserInt(l1,"byVLooseIsolationMVArun2v1DBoldDMwLT") == 1);
    //bool iso2 = (userdatahelpers::getUserInt(l2,"byVLooseIsolationMVArun2v1DBoldDMwLT") == 1);
    //bool iso1 = (userdatahelpers::getUserInt(l1,"byVVLooseIsolationMVArun2017v2DBoldDMwLT2017") == 1); //FRA 2017
    //bool iso2 = (userdatahelpers::getUserInt(l2,"byVVLooseIsolationMVArun2017v2DBoldDMwLT2017") == 1); //FRA 2017
    //bool iso1 = (userdatahelpers::getUserInt(l1,"byVVVLooseDeepTau2017v2p1VSjet") == 1);
    //bool iso2 = (userdatahelpers::getUserInt(l2,"byVVVLooseDeepTau2017v2p1VSjet") == 1);
    bool iso1 = (userdatahelpers::getUserInt(l1,"byVVVLooseDeepTau2017v2p1VSjet") == 1);
    bool iso2 = (userdatahelpers::getUserInt(l2,"byVVVLooseDeepTau2017v2p1VSjet") == 1);

    if (!iso1 || !iso2)
      return false;

    return true; // passed all requirements
  }

  else if (pType == kMuMu)
    return false;
  
  else if (pType == kEE)
    return false;
  
  else if (pType == kEMu)
    return false;
  
  else
  {
    // should never happen
    edm::LogWarning("Unrecognised pair") << "(ClassicSVfitInterface) Warning! could not assess the pair type, won't compute SVFit";
    return false;
  }
}

#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(ZCandidateFiller);
