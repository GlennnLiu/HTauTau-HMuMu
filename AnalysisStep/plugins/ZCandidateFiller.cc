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
#include "DataFormats/Math/interface/deltaR.h"

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
// using namespace classic_svFit;
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

  edm::EDGetTokenT<edm::View<reco::CompositeCandidate> > candidateToken;
  const StringCutObjectSelector<pat::CompositeCandidate, true> preBestZSelection;
  const edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  const edm::EDGetTokenT< edm::TriggerResults > triggerResultsToken_;

  int sampleType;
  int setup;
  const CutSet<pat::CompositeCandidate> cuts;
  int FSRMode;
  bool embedDaughterFloats;
  
  //For HLT match
  vector<string> muHLTPaths2_;
  vector<string> muHLTFilters2_;
  vector<string> muHLTFilters2_leg1;
  vector<string> muHLTFilters2_leg2;
  vector<float> muHLTFilters2_DZ;
  vector<float> muHLTFilters2_DR;

  vector<string> eleHLTPaths2_;
  vector<string> eleHLTFilters2_;
  vector<string> eleHLTFilters2_leg1;
  vector<string> eleHLTFilters2_leg2;
  vector<float> eleHLTFilters2_DZ;
  vector<float> eleHLTFilters2_DR;

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
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection> (iConfig.getParameter<edm::InputTag>("TriggerSet"))),
  triggerResultsToken_( consumes< edm::TriggerResults >( iConfig.getParameter< edm::InputTag >( "TriggerResults" ) ) ),
  sampleType(iConfig.getParameter<int>("sampleType")),
  setup(iConfig.getParameter<int>("setup")),
  cuts(iConfig.getParameter<edm::ParameterSet>("flags")),
  embedDaughterFloats(iConfig.getUntrackedParameter<bool>("embedDaughterFloats",true))
{
  
  string mode = iConfig.getParameter<string>("FSRMode");
  if      (mode == "skip")   FSRMode = 0;
  else if (mode == "Legacy") FSRMode = 1;
  else if (mode == "RunII")  FSRMode = 2;
  else {
    cout << "ZCandidateFiller: FSRMode " << FSRMode << " not supported" << endl;
    abort();
  }

  if (sampleType == 2016)
  {
	muHLTPaths2_ = 
	{
	"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",//DiMu
	"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*",
	"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*",
	"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*",
	};
	muHLTFilters2_ = 
	{
	"hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4DzFiltered0p2",
	"hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4DzFiltered0p2",
	"hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4",
	"hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4",
	};
	muHLTFilters2_leg1 = 
	{
	"hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4",
	"hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4",
	"hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4",
	"hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4",
	};
        muHLTFilters2_leg2 =
        {
	"hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4",
        "hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4",
	"hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4",
	"hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4",
	};
	muHLTFilters2_DZ = 
	{
	0.2,0.2,-1,-1,
	};
        muHLTFilters2_DR =
        {
	0.001,0.001,-1,-1,
        };

	eleHLTPaths2_ = 
	{
	"HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
	"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
	"HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*",
	};
        eleHLTFilters2_ =
        {
	"hltEle17Ele12CaloIdLTrackIdLIsoVLDZFilter",
        "hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter",
	"hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter",
	};
        eleHLTFilters2_leg1 =
        {
	"hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter",
	"hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter",
	"hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter",
        };
        eleHLTFilters2_leg2 =
        {
	"hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter",
        "hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter",
	"hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter",
	};
        eleHLTFilters2_DZ =
        {
	0.2,0.2,-1,
        };
        eleHLTFilters2_DR =
        {
	-1,-1,-1,
        };

  }
  else if (sampleType == 2017)
  {
	muHLTPaths2_ =
	{
	"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*",
	"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*",
	};
        muHLTFilters2_ =
        {
	"hltDiMuon178RelTrkIsoFiltered0p4DzFiltered0p2",//hltDiMuon178Mass3p8Filtered
        "hltDiMuon178RelTrkIsoFiltered0p4DzFiltered0p2",//hltDiMuon178Mass8Filtered
	};
        muHLTFilters2_leg1 =
        {
	"hltDiMuon178RelTrkIsoFiltered0p4",
        "hltDiMuon178RelTrkIsoFiltered0p4",
	};
        muHLTFilters2_leg2 =
        {
	"hltDiMuon178RelTrkIsoFiltered0p4",
        "hltDiMuon178RelTrkIsoFiltered0p4",
	};
        muHLTFilters2_DZ =
        {
	0.2,0.2,
        };
        muHLTFilters2_DR =
        {
	0.001,0.001,
        };

	eleHLTPaths2_ = 
	{
	"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
	"HLT_DoubleEle33_CaloIdL_MW_v*",
	};
        eleHLTFilters2_ =
        {
	"pass",
	"hltDiEle33CaloIdLMWPMS2UnseededFilter",
        };
        eleHLTFilters2_leg1 =
        {
	"hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter",
        "hltDiEle33CaloIdLMWPMS2UnseededFilter",
	};
        eleHLTFilters2_leg2 =
        {
	"hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter",
        "hltDiEle33CaloIdLMWPMS2UnseededFilter",
	};
        eleHLTFilters2_DZ =
        {
	-1,-1,
        };
        eleHLTFilters2_DR =
        {
	-1,-1,
        };

  }
  else if (sampleType == 2018)
  {
	muHLTPaths2_ = 
	{
	"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*",
	};
        muHLTFilters2_ =
        {
	"hltDiMuon178RelTrkIsoFiltered0p4DzFiltered0p2",//hltDiMuon178Mass3p8Filtered
        };
        muHLTFilters2_leg1 =
        {
	"hltDiMuon178RelTrkIsoFiltered0p4",
        };
        muHLTFilters2_leg2 =
        {
	"hltDiMuon178RelTrkIsoFiltered0p4",
        };
        muHLTFilters2_DZ =
        {
	0.2,
        };
        muHLTFilters2_DR =
        {
	0.001,
        };

	eleHLTPaths2_ = 
	{
	"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
	"HLT_DoubleEle25_CaloIdL_MW_v*",
	};
        eleHLTFilters2_ =
        {
	"pass",
	"hltDiEle25CaloIdLMWPMS2UnseededFilter",
        };
        eleHLTFilters2_leg1 =
        {
	"hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter",
	"hltDiEle25CaloIdLMWPMS2UnseededFilter",
        };
        eleHLTFilters2_leg2 =
        {
	"hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter",
	"hltDiEle25CaloIdLMWPMS2UnseededFilter",
        };
        eleHLTFilters2_DZ =
        {
	-1,-1,
        };
        eleHLTFilters2_DR =
        {
	-1,-1,
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

  edm::Handle< edm::TriggerResults > triggerResults;
  iEvent.getByToken( triggerResultsToken_, triggerResults );

  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);


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
//   cout<<LLCands->size()<<" bare Z candidates"<<endl;
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
    
    const Candidate *l1 = myCand.daughter(0);
    const Candidate *l2 = myCand.daughter(1);

    //pt sum, for sorting in ZZ cand
    float ptSum=l1->pt()+l2->pt();
    myCand.addUserFloat("ptSum",ptSum);

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
    if (abs(id0)*abs(id1)!=225) {
        myCand.addUserFloat("d0.isGoodTau",abs(id0)==15?goodTau:true);
        myCand.addUserFloat("d1.isGoodTau",abs(id1)==15?goodTau:true);
    }
    else {
	myCand.addUserFloat("d0.isGoodTau",myCand.userFloat("d0.isGood_Tau"));
	myCand.addUserFloat("d1.isGoodTau",myCand.userFloat("d1.isGood_Tau"));
    }

    //HLTMatch
    //The candidate passes HLTMatch if at least one lepton daughter matches the single trigger or both matches to the double trigger
    bool muHLTMatch1=false,muHLTMatch2=false,muHLTMatch=false;
    bool eleHLTMatch1=false,eleHLTMatch2=false,eleHLTMatch=false;
    
    if (abs(id0)==13 && abs(id1)==13) {
        if (myCand.userFloat("d0.HLTMatch1") || myCand.userFloat("d1.HLTMatch1"))
            muHLTMatch1=true;
	
	for (size_t idxto = 0; idxto < triggerObjects->size(); ++idxto) {
	    for (size_t jdxto = 0; jdxto < triggerObjects->size(); ++jdxto) {
		if (!myCand.hasUserInt("d0.TrgObj"+std::to_string(idxto)) || !myCand.hasUserInt("d1.TrgObj"+std::to_string(jdxto)))
		    continue;
	  	pat::TriggerObjectStandAlone obj1 = triggerObjects->at(idxto);
	        obj1.unpackFilterLabels(iEvent,*triggerResults );
                pat::TriggerObjectStandAlone obj2 = triggerObjects->at(jdxto);
                obj2.unpackFilterLabels(iEvent,*triggerResults );
		for (size_t j=0; j<muHLTPaths2_.size(); ++j) {
		    if (!( obj1.hasFilterLabel(muHLTFilters2_[j]) && obj2.hasFilterLabel(muHLTFilters2_[j]) ))
			continue;
		    if (!( obj1.hasFilterLabel(muHLTFilters2_leg1[j]) && obj2.hasFilterLabel(muHLTFilters2_leg2[j]) ) && !( obj1.hasFilterLabel(muHLTFilters2_leg2[j]) && obj2.hasFilterLabel(muHLTFilters2_leg1[j]) ))
			continue;
		    float dr=deltaR(obj1,obj2);
		    float dz=abs(obj1.vz()-obj2.vz());
		    if ( muHLTFilters2_DR[j] >= 0 && dr<muHLTFilters2_DR[j])
			continue;
		    if ( muHLTFilters2_DZ[j] >= 0 && dz>muHLTFilters2_DZ[j])
			continue;
		    muHLTMatch2=true;
		    break;
		}
	    }
	}
    }
    muHLTMatch=muHLTMatch1 || muHLTMatch2;

    if (abs(id0)==11 && abs(id1)==11) {
        if (myCand.userFloat("d0.HLTMatch1") || myCand.userFloat("d1.HLTMatch1"))
            eleHLTMatch1=true;

        for (size_t idxto = 0; idxto < triggerObjects->size(); ++idxto) {
            for (size_t jdxto = 0; jdxto < triggerObjects->size(); ++jdxto) {
                if (!myCand.hasUserInt("d0.TrgObj"+std::to_string(idxto)) || !myCand.hasUserInt("d1.TrgObj"+std::to_string(jdxto)))
                    continue;
                pat::TriggerObjectStandAlone obj1 = triggerObjects->at(idxto);
                obj1.unpackFilterLabels(iEvent,*triggerResults );
                pat::TriggerObjectStandAlone obj2 = triggerObjects->at(jdxto);
                obj2.unpackFilterLabels(iEvent,*triggerResults );
                for (size_t j=0; j<eleHLTPaths2_.size(); ++j) {
                    if (!( obj1.hasFilterLabel(eleHLTFilters2_[j]) && obj2.hasFilterLabel(eleHLTFilters2_[j]) ))
                        continue;
                    if (!( obj1.hasFilterLabel(eleHLTFilters2_leg1[j]) && obj2.hasFilterLabel(eleHLTFilters2_leg2[j]) ) && !( obj1.hasFilterLabel(eleHLTFilters2_leg2[j]) && obj2.hasFilterLabel(eleHLTFilters2_leg1[j]) ))
                        continue;
                    float dr=deltaR(obj1,obj2);
                    float dz=abs(obj1.vz()-obj2.vz());
                    if ( eleHLTFilters2_DR[j] >= 0 && dr<eleHLTFilters2_DR[j])
                        continue;
                    if ( eleHLTFilters2_DZ[j] >= 0 && dz>eleHLTFilters2_DZ[j])
                        continue;
		    if (eleHLTPaths2_[j].find("Mass3p8")!=string::npos) {
			if (!( obj1.hasFilterLabel("hltDiMuon178Mass3p8Filtered") && obj2.hasFilterLabel("hltDiMuon178Mass3p8Filtered") ))
			    continue;
			float invMass=(obj1.p4()+obj2.p4()).M();
			if (invMass<3.8)
			    continue;
			}
		    if (eleHLTPaths2_[j].find("Mass8")!=string::npos) {
                        if (!( obj1.hasFilterLabel("hltDiMuon178Mass8Filtered") && obj2.hasFilterLabel("hltDiMuon178Mass8Filtered") ))
                            continue;
                        float invMass=(obj1.p4()+obj2.p4()).M();
                        if (invMass<8.)
                            continue;
			}
                    eleHLTMatch2=true;
                    break;
                }
            }
        }
    }
    eleHLTMatch=eleHLTMatch1 || eleHLTMatch2;

    myCand.addUserFloat("muHLTMatch1",muHLTMatch1);
    myCand.addUserFloat("muHLTMatch2",muHLTMatch2);
    myCand.addUserFloat("muHLTMatch",muHLTMatch);
    myCand.addUserFloat("eleHLTMatch1",eleHLTMatch1);
    myCand.addUserFloat("eleHLTMatch2",eleHLTMatch2);
    myCand.addUserFloat("eleHLTMatch",eleHLTMatch);


    myCand.addUserFloat("OSSF",OS && SF);
    //--- Find "best Z" (closest to mZ) among those passing the "bestZAmong" selection (2011 PRL logic), now deprecated!!!
    if (preBestZSelection(myCand)) {
      float diffZmass = fabs(ZmassValue - myCand.mass());//userFloat("goodMass"));
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
//   cout<<result->size()<<" Z candidates"<<endl;
  iEvent.put(std::move(result));
}

#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(ZCandidateFiller);
