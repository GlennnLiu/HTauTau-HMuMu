/** \class MuFiller
 *
 *  No description available.
 *
 *  $Date: 2013/05/14 10:08:19 $
 *  $Revision: 1.18 $
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

#include <DataFormats/PatCandidates/interface/Muon.h>
#include "DataFormats/VertexReco/interface/Vertex.h"
// #include <DataFormats/Common/interface/TriggerResults.h>
#include "DataFormats/Math/interface/deltaR.h"

#include <HTauTauHMuMu/AnalysisStep/interface/CutSet.h>
#include <HTauTauHMuMu/AnalysisStep/interface/LeptonIsoHelper.h>

#include "MuonMVAReader/Reader/interface/MuonGBRForestReader.hpp"

#include <vector>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;


class MuFiller : public edm::EDProducer {
   public:
   /// Constructor
   explicit MuFiller(const edm::ParameterSet&);
   
   /// Destructor
   ~MuFiller();
   
   private:
   virtual void beginJob(){};
   virtual void produce(edm::Event&, const edm::EventSetup&);
   virtual void endJob(){};
   
   edm::EDGetTokenT<pat::MuonRefVector> muonToken;
   //edm::EDGetTokenT<vector<pat::Muon> > muonToken;

   int sampleType;
   int setup;
   const StringCutObjectSelector<pat::Muon, true> cut;
   // const edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
   // const edm::EDGetTokenT< edm::TriggerResults > triggerResultsToken_;
   const CutSet<pat::Muon> flags;
   edm::EDGetTokenT<double> rhoToken;
   edm::EDGetTokenT<vector<Vertex> > vtxToken;
   
   // MVA Reader
   MuonGBRForestReader *r;
};


MuFiller::MuFiller(const edm::ParameterSet& iConfig) :
//muonToken(consumes<vector<pat::Muon> >(iConfig.getParameter<edm::InputTag>("src"))),
muonToken(consumes<pat::MuonRefVector>(iConfig.getParameter<edm::InputTag>("src"))),
sampleType(iConfig.getParameter<int>("sampleType")),
setup(iConfig.getParameter<int>("setup")),
cut(iConfig.getParameter<std::string>("cut")),
// triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection> (iConfig.getParameter<edm::InputTag>("TriggerSet"))),
// triggerResultsToken_( consumes< edm::TriggerResults >( iConfig.getParameter< edm::InputTag >( "TriggerResultsLabel" ) ) ),
flags(iConfig.getParameter<edm::ParameterSet>("flags"))
{
   rhoToken = consumes<double>(LeptonIsoHelper::getMuRhoTag(sampleType, setup));
   vtxToken = consumes<vector<Vertex> >(edm::InputTag("goodPrimaryVertices"));
   produces<pat::MuonCollection>();
   
   // MVA Reader
   r = new MuonGBRForestReader(setup, 2);
	
}

MuFiller::~MuFiller(){}


void
MuFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   
   //--- Get leptons and rho
   //edm::Handle<vector<pat::Muon> > muonHandle;
   edm::Handle<pat::MuonRefVector> muonHandle;
   iEvent.getByToken(muonToken, muonHandle);
   //const vector<pat::Muon> *inputMuons = muonHandle.product();  

 
   // edm::Handle< edm::TriggerResults > triggerResults;
   // iEvent.getByToken( triggerResultsToken_, triggerResults );
   //const edm::TriggerNames &names = iEvent.triggerNames(*triggerResults);
  
   // edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
   // iEvent.getByToken(triggerObjects_, triggerObjects);
 
   edm::Handle<double> rhoHandle;
   iEvent.getByToken(rhoToken, rhoHandle);
   double rho = *rhoHandle;
   
   edm::Handle<vector<Vertex> > vertices;
   iEvent.getByToken(vtxToken,vertices);
   
   
   // Output collection
   auto result = std::make_unique<pat::MuonCollection>();
   

   //for (unsigned i=0; i<inputMuons->size(); ++i) {
      //pat::Muon l = inputMuons->at(i);
   for (unsigned int i = 0; i< muonHandle->size(); ++i){
      //---Clone the pat::Muon
      pat::Muon l(*((*muonHandle)[i].get()));
      
      //--- PF ISO
      // for cone size R=0.3 :
      float PFChargedHadIso   = l.pfIsolationR03().sumChargedHadronPt;
      float PFNeutralHadIso   = l.pfIsolationR03().sumNeutralHadronEt;
      float PFPhotonIso       = l.pfIsolationR03().sumPhotonEt;
      float PFPUChargedHadIso = l.pfIsolationR03().sumPUPt;
      
      float combRelIsoPF = LeptonIsoHelper::combRelIsoPF(sampleType, setup, rho, l);
      
      //--- SIP, dxy, dz
      float IP      = std::abs(l.dB(pat::Muon::PV3D));
      float IPError = l.edB(pat::Muon::PV3D);
      float SIP     = IP/IPError;
      float dxy = std::abs(l.dB(pat::Muon::PV2D));
      float dz  = std::abs(l.dB(pat::Muon::PVDZ));
      
      float mediumID;
      
      bool goodGlob = l.isGlobalMuon() && 
         l.globalTrack()->normalizedChi2() < 3 && 
         l.combinedQuality().chi2LocalPosition < 12 && 
         l.combinedQuality().trkKink < 20; 
      mediumID = muon::isLooseMuon(l) && 
         l.innerTrack()->validFraction() > 0.49 && 
         muon::segmentCompatibility(l) > (goodGlob ? 0.303 : 0.451); 

      l.addUserFloat("mediumID",mediumID);

      if(!l.hasUserFloat("isPFMuon")) {
         l.addUserFloat("isPFMuon",l.isPFMuon());
      }
      if(!l.hasUserFloat("isGlobalMuon")) {
         l.addUserFloat("isGlobalMuon",l.isGlobalMuon());
      }
      if(!l.hasUserFloat("isTrackerMuon")) {
         l.addUserFloat("isTrackerMuon",l.isTrackerMuon());
      }

      //--- Muon Timing
      float muontime = 0;
      if (l.time().nDof>4) muontime= l.time().timeAtIpInOut;
      
      
      //--- Embed user variables
      l.addUserFloat("PFChargedHadIso",PFChargedHadIso);
      l.addUserFloat("PFNeutralHadIso",PFNeutralHadIso);
      l.addUserFloat("PFPhotonIso",PFPhotonIso);
      l.addUserFloat("PFPUChargedHadIso",PFPUChargedHadIso);
      l.addUserFloat("combRelIsoPF",combRelIsoPF);
      l.addUserFloat("rho",rho);
      l.addUserFloat("SIP",SIP);
      l.addUserFloat("dxy",dxy);
      l.addUserFloat("dz",dz);
      // l.addUserFloat("BDT",BDT);
      // l.addUserFloat("isBDT",isBDT);
      // l.addUserCand("MCMatch",genMatch); // FIXME
      l.addUserFloat("time",muontime);
      
      if (!l.hasUserFloat("correctedPtError")) {
         l.addUserFloat("correctedPtError",l.muonBestTrack()->ptError()); //This is expected by the kin fitter
      }
      
      //--- MC parent code
      //     MCHistoryTools mch(iEvent);
      //     if (mch.isMC()) {
      //       int MCParentCode = 0;//FIXME: does not work on cmg mch.getParentCode((l.genParticleRef()).get());
      //       l.addUserFloat("MCParentCode",MCParentCode);
      //     }
      
      //--- Check selection cut. Being done here, flags are not available; but this way we
      //    avoid wasting time on rejected leptons.
      if (!cut(l)) continue;
      
      //--- Embed flags (ie flags specified in the "flags" pset)
      for(CutSet<pat::Muon>::const_iterator flag = flags.begin(); flag != flags.end(); ++flag) {
         l.addUserFloat(flag->first,int((*(flag->second))(l)));
      }
      
      result->push_back(l);
   }
   //cout<<result->size()<<" soft muons"<<endl;
   iEvent.put(std::move(result));
}


#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(MuFiller);

