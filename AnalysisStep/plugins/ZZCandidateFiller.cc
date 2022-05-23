/** \class ZZCandidateFiller
 *
 *
 *  \author N. Amapane - Torino
 *  \author C. Botta - Torino
 *  \author G. Ortona - LLR
 */

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/ParameterSet/interface/FileInPath.h>

#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/Jet.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/GeometryVector/interface/Point3DBase.h>
#include <DataFormats/METReco/interface/PFMET.h>
#include <DataFormats/METReco/interface/PFMETCollection.h>

#include <SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h>

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h>
#include <RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h>
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>

#include <HTauTauHMuMu/AnalysisStep/interface/CutSet.h>
#include <HTauTauHMuMu/AnalysisStep/interface/DaughterDataHelpers.h>
#include <HTauTauHMuMu/AnalysisStep/interface/FinalStates.h>
#include <HTauTauHMuMu/AnalysisStep/interface/CompositeCandMassResolution.h>

#include <HTauTauHMuMu/AnalysisStep/interface/Fisher.h>
#include <HTauTauHMuMu/AnalysisStep/interface/Comparators.h>
#include <HTauTauHMuMu/AnalysisStep/interface/utils.h>
#include <HTauTauHMuMu/AnalysisStep/interface/LeptonIsoHelper.h>
#include <HTauTauHMuMu/AnalysisStep/interface/JetCleaner.h>

#include <KinZfitter/KinZfitter/interface/KinZfitter.h>

#include <HHKinFit2/interface/HHKinFitMasterHeavyHiggs.h>
#include <HHKinFit2/interface/exceptions/HHInvMConstraintException.h>
#include <HHKinFit2/interface/exceptions/HHEnergyRangeException.h>
#include <HHKinFit2/interface/exceptions/HHEnergyConstraintException.h>

#include <HHKinFit2/src/HHKinFitMasterHeavyHiggs.cpp>
#include <HHKinFit2/src/HHKinFit.cpp>
#include <HHKinFit2/src/HHFitObjectEConstM.cpp>
#include <HHKinFit2/src/HHFitObjectEConstBeta.cpp>
#include <HHKinFit2/src/HHFitObjectE.cpp>
#include <HHKinFit2/src/HHFitObjectMET.cpp>
#include <HHKinFit2/src/HHFitObject.cpp>
#include <HHKinFit2/src/HHFitObjectComposite.cpp>
#include <HHKinFit2/src/HHFitConstraint4Vector.cpp>
#include <HHKinFit2/src/HHFitConstraint4VectorBJet.cpp>
#include <HHKinFit2/src/HHFitConstraint.cpp>
#include <HHKinFit2/src/HHFitConstraintEHardM.cpp>
#include <HHKinFit2/src/HHFitConstraintSoftBoundary.cpp>
#include <HHKinFit2/src/HHLorentzVector.cpp>
#include <HHKinFit2/src/PSMath.cpp>


#include "TH2F.h"
#include "TFile.h"
#include "TLorentzVector.h"

#include <string>
#include <vector>

using namespace htautauhmumu;
//using namespace BranchHelpers;

typedef std::pair<int, TLorentzVector> SimpleParticle_t;
typedef std::vector<SimpleParticle_t> SimpleParticleCollection_t;

bool doVtxFit = false;

class ZZCandidateFiller : public edm::EDProducer {
public:
  /// Constructor
  explicit ZZCandidateFiller(const edm::ParameterSet&);

  /// Destructor
  ~ZZCandidateFiller();

private:
  typedef map<const reco::Candidate*, const pat::PFParticle*> FSRToLepMap;

  virtual void beginJob(){};
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};

  void getPairMass(const reco::Candidate* lp, const reco::Candidate* lm, FSRToLepMap& photons, float& mass, int& ID);
//  void getSVMass(const reco::Candidate* l1, const reco::Candidate* l2, 

  edm::EDGetTokenT<edm::View<reco::CompositeCandidate> > candidateToken;
  edm::EDGetTokenT<edm::View<pat::MET> > metToken;
  edm::EDGetTokenT<math::Error<2>::type> metCovToken;

  const CutSet<pat::CompositeCandidate> preBestCandSelection;
  const CutSet<pat::CompositeCandidate> cuts;
  int sampleType;
  int setup;

  bool embedDaughterFloats;
  bool ZRolesByMass;
  reco::CompositeCandidate::role_collection rolesZ1Z2;
  reco::CompositeCandidate::role_collection rolesZ2Z1;
  bool isMC;
  bool doKinFit,doKinFitOld;
  bool debug;

  // float muon_iso_cut, electron_iso_cut;
  TH2F* corrSigmaMu;
  TH2F* corrSigmaEle;
  Comparators::ComparatorTypes bestCandType;
  KinZfitter *kinZfitter;
  edm::EDGetTokenT<edm::View<pat::Jet> > jetToken;
//  edm::EDGetTokenT<pat::METCollection> metToken;
  edm::EDGetTokenT<edm::View<reco::Candidate> > softLeptonToken;
  edm::EDGetTokenT<edm::View<reco::CompositeCandidate> > ZCandToken;
   
  int _num_of_JEC_variations;
};


ZZCandidateFiller::ZZCandidateFiller(const edm::ParameterSet& iConfig) :
  candidateToken(consumes<edm::View<reco::CompositeCandidate> >(iConfig.getParameter<edm::InputTag>("src"))),
  preBestCandSelection(iConfig.getParameter<edm::ParameterSet>("bestCandAmong")),
  cuts(iConfig.getParameter<edm::ParameterSet>("flags")),
  sampleType(iConfig.getParameter<int>("sampleType")),
  setup(iConfig.getParameter<int>("setup")),
  embedDaughterFloats(iConfig.getUntrackedParameter<bool>("embedDaughterFloats", true)),
  ZRolesByMass(iConfig.getParameter<bool>("ZRolesByMass")),
  isMC(iConfig.getParameter<bool>("isMC")),
  doKinFit(iConfig.getParameter<bool>("doKinFit")),
  doKinFitOld(iConfig.getParameter<bool>("doKinFitOld")),
  debug(iConfig.getParameter<bool>("debug")),
  corrSigmaMu(0),
  corrSigmaEle(0),
  kinZfitter(0)
{
  produces<pat::CompositeCandidateCollection>();

  jetToken = consumes<edm::View<pat::Jet> >(edm::InputTag("cleanJets"));
  metToken = consumes<edm::View<pat::MET> >(iConfig.getParameter<edm::InputTag>("srcMET"));//("slimmedMETs"));
  metCovToken = consumes<math::Error<2>::type>(iConfig.getParameter<edm::InputTag>("srcCov"));
  softLeptonToken = consumes<edm::View<reco::Candidate> >(edm::InputTag("softLeptons"));
  ZCandToken = consumes<edm::View<reco::CompositeCandidate> >(edm::InputTag("ZCand"));

  rolesZ1Z2 = {"Z1", "Z2"};
  rolesZ2Z1 = {"Z2", "Z1"};

  _num_of_JEC_variations = 5; // Define number of total JEC variations 4 (JESUp, JESDn, JERUp, JERDn) + 1 for Nominal

  if (setup < 2015) {// FIXME:  EbE corrections to be updated for Run II
    // Run I ebe corrections; obsolete
    edm::FileInPath fip("HTauTauHMuMu/AnalysisStep/data/ebeOverallCorrections.Legacy2013.v0.root");
    std::string ebePath=fip.fullPath();

    // EbE corrections
    TFile* fCorrSigma = new TFile(ebePath.data()); // FIXME: is leaked
    std::string sigmaCorrType = (isMC?"mc":"reco");
    std::string sigmaCorrYear = "";
    if (setup==2011) sigmaCorrYear = "42x";
    else if (setup==2012) sigmaCorrYear = "53x";

    corrSigmaMu=  (TH2F*)fCorrSigma->Get(("mu_"+sigmaCorrType+sigmaCorrYear).data());
    corrSigmaEle= (TH2F*)fCorrSigma->Get(("el_"+sigmaCorrType+sigmaCorrYear).data());
  }

  string cmp=iConfig.getParameter<string>("bestCandComparator");
  if      (cmp=="byBestZ1bestZ2") bestCandType=Comparators::byBestZ1bestZ2;
  else if (cmp=="byBestKD")       bestCandType=Comparators::byBestKD;
  else if (cmp=="byBestKD_VH")    bestCandType=Comparators::byBestKD_VH;
  else if (cmp=="byBestPsig")    bestCandType=Comparators::byBestPsig;
  else if (cmp=="byMHWindow")    bestCandType=Comparators::byMHWindow;
  else abort();
  //cout<<iConfig.getParameter<edm::InputTag>("src").label()<<endl;
  //-- kinematic refitter
  kinZfitter = new KinZfitter(!isMC);
  // No longer used, but keept for future needs
//   muon_iso_cut = iConfig.getParameter<double>("muon_iso_cut");
//   electron_iso_cut = iConfig.getParameter<double>("electron_iso_cut");
}

ZZCandidateFiller::~ZZCandidateFiller(){
  delete kinZfitter;
}


void ZZCandidateFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;
  using namespace std;
  using namespace reco;

  auto result = std::make_unique<pat::CompositeCandidateCollection>();

  const float ZmassValue = 91.1876;//PDGHelpers::Zmass;

  // Get LLLL candidates
  if (debug) cout<<"Get LLLL candidates"<<endl;
  Handle<edm::View<CompositeCandidate> > LLLLCands;
  iEvent.getByToken(candidateToken, LLLLCands);

  // Get jets
  if (debug) cout<<"Get jets"<<endl;
  Handle<edm::View<pat::Jet> > CleanJets;
  iEvent.getByToken(jetToken, CleanJets);

  // Get MET
  if (debug) cout<<"Get MET"<<endl;
  float PFMET = 0.;
  float PFMETPhi = 0.;
  Handle<edm::View<pat::MET> > metHandle;
  iEvent.getByToken(metToken, metHandle);
  if(metHandle.isValid()){
    PFMET = metHandle->front().pt();
    PFMETPhi = metHandle->front().phi();
  }
  TVector2 ptmiss(metHandle->front().px(),metHandle->front().py());
  
  Handle<math::Error<2>::type> covHandle;
  iEvent.getByToken (metCovToken, covHandle);
  TMatrixD covMET(2, 2);
  covMET[0][0] = (*covHandle)(0,0);
  covMET[1][0] = (*covHandle)(1,0);
  covMET[0][1] = covMET[1][0]; // (1,0) is the only one saved
  covMET[1][1] = (*covHandle)(1,1);
  
  // FIXME: May need to correct MET in the MC and use the corrected MET to calculate the VH MEs

  // Get leptons (in order to store extra leptons)
  if (debug) cout<<"Get leptons"<<endl;
  Handle<View<reco::Candidate> > softleptoncoll;
  iEvent.getByToken(softLeptonToken, softleptoncoll);
  vector<reco::CandidatePtr> goodisoleptonPtrs;
  for( View<reco::Candidate>::const_iterator lep = softleptoncoll->begin(); lep != softleptoncoll->end(); ++ lep ){
    if((bool)userdatahelpers::getUserFloat(&*lep,"isGood")
       //       && (bool)userdatahelpers::getUserFloat(&*lep,"isIsoFSRUncorr") // with old FSR strategy
       && (bool)userdatahelpers::getUserFloat(&*lep,"passCombRelIsoPFFSRCorr") // with new FSR strategy
       ){
      const reco::CandidatePtr lepPtr(softleptoncoll,lep-softleptoncoll->begin());
      goodisoleptonPtrs.push_back(lepPtr);
    }
  }

  // Get Z Candidates
  if (debug) cout<<"Get Z Candidates"<<endl;
  Handle<View<CompositeCandidate> > ZCands;
  iEvent.getByToken(ZCandToken, ZCands);

  // Get processID
//   edm::Handle<GenEventInfoProduct> gen;
//   iEvent.getByLabel( "generator", gen );
//   int processID = gen->signalProcessID();

  // to calculate mass resolution
  if (debug) cout<<"Calculate mass resolution"<<endl;
  CompositeCandMassResolution errorBuilder;
  errorBuilder.init(iSetup);

  vector<int> bestCandIdx(preBestCandSelection.size(),-1);
  vector<float> maxPtSum(preBestCandSelection.size(),-1);
  vector< vector<int> > preSelCands(preBestCandSelection.size());
  //----------------------------------------------------------------------
  //--- Loop over input candidates
  if (debug) cout<<"Loop over candidates"<<endl;
  cout<<LLLLCands->size()<<" bare ZZ candidates"<<endl;
  for( View<CompositeCandidate>::const_iterator cand = LLLLCands->begin(); cand != LLLLCands->end(); ++ cand ) {
    int icand = distance(LLLLCands->begin(),cand);

    pat::CompositeCandidate myCand(*cand);

    if (embedDaughterFloats){
      userdatahelpers::embedDaughterData(myCand);
    }

    //--- Set id of the Z1 and "Z1"/"Z2" labels. This allows to call e.g. aHiggs->daughter("Z1").
    // if ZRolesByMass is true, 'Z1' and iZ1 refer to the  Z closest to mZ; otherwise Z1 = daughter(0). The latter is used for control regions.
    if (debug) cout<<"Role Z by mass"<<endl;
    const reco::CompositeCandidate::role_collection* ZRoles = &rolesZ1Z2;
    int iZ1 = 0;
    int iZ2 = 1;
    if (ZRolesByMass) {
      if(std::abs(myCand.userFloat("d0.goodMass")-ZmassValue)>=std::abs(myCand.userFloat("d1.goodMass")-ZmassValue)){
        swap(iZ1,iZ2);
        ZRoles = &rolesZ2Z1;
      }
    }
    myCand.setRoles(*ZRoles);
    myCand.applyRoles();

    //--- Z pointers
    if (debug) cout<<"Z pointers"<<endl;
    const reco::Candidate* Z1= myCand.daughter(iZ1);
    const reco::Candidate* Z2= myCand.daughter(iZ2);
    vector<const reco::Candidate*> Zs = {Z1, Z2}; // in the original order

    //--- Lepton pointers in the original order
    if (debug) cout<<"Lepton pointers"<<endl;
    const reco::Candidate* Z1L1= Z1->daughter(0);
    const reco::Candidate* Z1L2= Z1->daughter(1);
    const reco::Candidate* Z2L1= Z2->daughter(0);
    const reco::Candidate* Z2L2= Z2->daughter(1);
    vector<const reco::Candidate*> ZZLeps = {Z1L1,Z1L2,Z2L1,Z2L2}; // array, in the original order

    // Create corresponding array of fourmomenta; will add FSR (below)
    vector<math::XYZTLorentzVector> pij(4);
    std::transform(ZZLeps.begin(), ZZLeps.end(),pij.begin(), [](const reco::Candidate* c){return c->p4();});

    //--- Collect FSR photons and map them to the corresponding leptons
    if (debug) cout<<"FSR photons"<<endl;
    FSRToLepMap FSRMap;
    for (unsigned iZ=0; iZ<2; ++iZ) {
      for (unsigned ifsr=2; ifsr<Zs[iZ]->numberOfDaughters(); ++ifsr) {
    const pat::PFParticle* fsr = static_cast<const pat::PFParticle*>(Zs[iZ]->daughter(ifsr));
    int ilep = iZ*2+fsr->userFloat("leptIdx");
    FSRMap[ZZLeps[ilep]]= fsr;
    pij[ilep]+=fsr->p4();
      }
    }

    //--- Lepton four-vectors in the original order; with FSR added
    math::XYZTLorentzVector p11 = pij[0];
    math::XYZTLorentzVector p12 = pij[1];
    math::XYZTLorentzVector p21 = pij[2];
    math::XYZTLorentzVector p22 = pij[3];
    int id11 = Z1L1->pdgId();
    int id12 = Z1L2->pdgId();
    int id21 = Z2L1->pdgId();
    int id22 = Z2L2->pdgId();
    int candChannel = id11*id12*id21*id22;

    if((id11 == 22 && id12 == 22) || (id21 == 22 && id22 == 22)) LogError("Z with 2 tle") << "Found a Z candidate made up of 2 trackless electrons";

    if(id11 == 22) id11 = -1 * id12;
    if(id12 == 22) id12 = -1 * id11;
    if(id21 == 22) id21 = -1 * id22;
    if(id22 == 22) id22 = -1 * id21;


    // Compute worst-lepton isolation
    if (debug) cout<<"Worst lepton isolation"<<endl;
    for (int zIdx=0; zIdx<2; ++zIdx) {
      float worstMuIso=0;
      float worstEleIso=0;
      for (int dauIdx=0; dauIdx<2; ++dauIdx) {
    const reco::Candidate* z = myCand.daughter(zIdx);
    const reco::Candidate* d = z->daughter(dauIdx);
    float combRelIsoPFCorr = userdatahelpers::getUserFloat(d,"combRelIsoPFFSRCorr");
    if (d->isMuon())		worstMuIso  = max(worstMuIso,  combRelIsoPFCorr);
    else if (d->isElectron())	worstEleIso = max(worstEleIso, combRelIsoPFCorr);
      }
      string base = (zIdx==0?"d0.":"d1.");
      myCand.addUserFloat(base+"worstMuIso",worstMuIso);
      myCand.addUserFloat(base+"worstEleIso",worstEleIso);
    }


    //----------------------------------------------------------------------
    //--- Alternative lepton pairings: "smart cut" and QCD suppression and

    //--- Sign-ordered leptons and leptopn four-vectors (without FSR), to be used to compute mZa, mZb, mZalpha, mZbeta
    if (debug) cout<<"Order by charge"<<endl;
    const reco::Candidate* Z1Lp(Z1L1);
    const reco::Candidate* Z1Lm(Z1L2);
    const reco::Candidate* Z2Lp(Z2L1);
    const reco::Candidate* Z2Lm(Z2L2);

    // Sort leptons for OS Z candidates; no sorting for the same-sign collections used for CRs
    if (Z1Lp->charge() < 0 && Z1Lp->charge()*Z1Lm->charge()<0) {
      swap(Z1Lp,Z1Lm);
    }
    if (Z2Lp->charge() < 0 && Z2Lp->charge()*Z2Lm->charge()<0) {
      swap(Z2Lp,Z2Lm);
    }

    math::XYZTLorentzVector p1p(Z1Lp->p4());
    math::XYZTLorentzVector p1m(Z1Lm->p4());
    math::XYZTLorentzVector p2p(Z2Lp->p4());
    math::XYZTLorentzVector p2m(Z2Lm->p4());

    if (debug) cout<<"Other combination options and smart cut"<<endl;
    // Build the other SF/OS combination
    float mZ1= userdatahelpers::getUserFloat(Z1,"goodMass");
    float mZa, mZb;
    int ZaID, ZbID;
    getPairMass(Z1Lp,Z2Lm,FSRMap,mZa,ZaID);
    getPairMass(Z1Lm,Z2Lp,FSRMap,mZb,ZbID);

    // For same-sign CRs, the Z2 leptons are same sign, so we need to check also the other combination.
    float mZalpha, mZbeta;
    int ZalphaID, ZbetaID;
    getPairMass(Z1Lp,Z2Lp,FSRMap,mZalpha,ZalphaID);
    getPairMass(Z1Lm,Z2Lm,FSRMap,mZbeta,ZbetaID);

    // Sort (mZa,mZb and) (mZalpha,mZbeta) so that a and alpha are the ones closest to mZ
    if (std::abs(mZa-ZmassValue)>=std::abs(mZb-ZmassValue)) {
      swap(mZa,mZb);
      swap(ZaID,ZbID);
    }
    if (std::abs(mZalpha-ZmassValue)>=std::abs(mZbeta-ZmassValue)) {
      swap(mZalpha,mZbeta);
      swap(ZalphaID,ZbetaID);
    }

    // "smart cut" mll logic: veto the candidate if by swapping leptons we find a better Z1 and the Z2 is below 12 GeV.
    // To handle same-sign CRs, we have to check both alternate pairings, and consider those that have a SF/OS Z1.
    bool passSmartMLL = true;
    if (((ZaID==-121||ZaID==-169) && std::abs(mZa-ZmassValue)<std::abs(mZ1-ZmassValue) && mZb<12) ||
        ((ZalphaID==-121||ZalphaID==-169) && std::abs(mZalpha-ZmassValue)<std::abs(mZ1-ZmassValue) && mZbeta<12)) passSmartMLL = false;


    //--- QCD suppression cut
    if (debug) cout<<"QCD suppression cut"<<endl;
    vector<const reco::Candidate*> lep;
    lep.push_back(Z1Lm);
    lep.push_back(Z1Lp);
    lep.push_back(Z2Lm);
    lep.push_back(Z2Lp);

    float mll6 = 9999;
    float mll4 = 9999;
    for (int i=0;i<4;++i) {
      for (int j=i+1;j<4;++j) {
        float mll = (lep[i]->p4()+lep[j]->p4()).mass();
        mll6 = min(mll, mll6);
        if (lep[i]->charge()*lep[j]->charge()<0) { //OS
          mll4 = min (mll,mll4);
        }
      }
    }


    //--- worst SIP value
    if (debug) cout<<"Worst SIP"<<endl;
    vector<double> SIPS = {myCand.userFloat("d0.d0.SIP"), myCand.userFloat("d0.d1.SIP"), myCand.userFloat("d1.d0.SIP"), myCand.userFloat("d1.d1.SIP")};
    sort(SIPS.begin(),SIPS.end());
    float SIP4 = SIPS[3];

    //--- Sorted pTs
    if (debug) cout<<"Sort pt"<<endl;
    vector<pair<double, int>> ptS;
    ptS.push_back(make_pair(Z1Lm->pt(), Z1Lm->pdgId()));
    ptS.push_back(make_pair(Z1Lp->pt(), Z1Lp->pdgId()));
    ptS.push_back(make_pair(Z2Lm->pt(), Z2Lm->pdgId()));
    ptS.push_back(make_pair(Z2Lp->pt(), Z2Lp->pdgId()));
    sort(ptS.begin(),ptS.end());

    //--- Mass and Lepton uncertainties
    if (debug) cout<<"Mass and lepton uncertainties"<<endl;
    std::vector<double> errs;
    float massError;
    if (abs(Z1Lm->pdgId())==15 || abs(Z1Lp->pdgId())==15 || abs(Z2Lm->pdgId())==15 || abs(Z2Lp->pdgId())==15 ){
	massError=0;
	for (size_t i = 0; i<myCand.daughter(0)->numberOfDaughters()+myCand.daughter(1)->numberOfDaughters(); ++i)
	    errs.push_back(0.);
    }
    else{
    	massError = errorBuilder.getMassResolutionWithComponents(myCand, errs);
    }
    int offset =0;
    float sigma[2][3] = {{0,0,0}, {0,0,0}};

    myCand.addUserFloat("massError",      massError);
    myCand.addUserFloat("massError11",    errs[0]);
    sigma[0][0] = errs[0];
    myCand.addUserFloat("massError12",    errs[1]);
    sigma[0][1] = errs[1];
    if (myCand.daughter(0)->numberOfDaughters()==3){
      myCand.addUserFloat("massError13",    errs[2]);
      sigma[0][2] = errs[2];
      offset = 1;
    }
    myCand.addUserFloat("massError21",    errs[2+offset]);
    sigma[1][0] = errs[2+offset];
    myCand.addUserFloat("massError22",    errs[3+offset]);
    sigma[1][1] = errs[3+offset];
    if (myCand.daughter(1)->numberOfDaughters()==3){
      myCand.addUserFloat("massError23",    errs[4+offset]);
      sigma[1][2]=errs[4+offset];
    }

    float massErrorCorr=0;
    for (int i=0; i<2; ++i) {
      for (int j=0; j<2; ++j) {
        const reco::Candidate* l=cand->daughter(i)->daughter(j);
        const TH2F* h;
        if (l->isMuon()) h = corrSigmaMu;
        else             h = corrSigmaEle;
    float ebecorr=1.;
    if (h!=0) {
      int ptBin  = min(max(1,h->GetXaxis()->FindBin(l->pt())), h->GetNbinsX());
      int etaBin = min(max(1,h->GetYaxis()->FindBin(fabs(l->eta()))), h->GetNbinsY());
      ebecorr = h->GetBinContent(ptBin, etaBin);
    }
        massErrorCorr+= (sigma[i][j]*ebecorr)*(sigma[i][j]*ebecorr);
      }
    }
    massErrorCorr += (sigma[0][2])*(sigma[0][2]);
    massErrorCorr += (sigma[1][2])*(sigma[1][2]);
    massErrorCorr = sqrt(massErrorCorr);
    myCand.addUserFloat("massErrorCorr",      massErrorCorr);


    //--- store good isolated leptons that are not involved in the current ZZ candidate
    if (debug) cout<<"Extra leptons"<<endl;
    int nExtraLep = 0;
    SimpleParticleCollection_t associatedLeptons;
    for (vector<reco::CandidatePtr>::const_iterator lepPtr = goodisoleptonPtrs.begin(); lepPtr != goodisoleptonPtrs.end(); ++lepPtr){
      const reco::Candidate* lep = lepPtr->get();
      if (
        reco::deltaR(lep->p4(), Z1L1->p4()) > 0.02 &&
        reco::deltaR(lep->p4(), Z1L2->p4()) > 0.02 &&
        reco::deltaR(lep->p4(), Z2L1->p4()) > 0.02 &&
        reco::deltaR(lep->p4(), Z2L2->p4()) > 0.02
        ){
        nExtraLep++;
        myCand.addUserCand("ExtraLep"+to_string(nExtraLep), *lepPtr);

        SimpleParticle_t theLepton(
        lep->pdgId(),
        TLorentzVector(lep->p4().x(), lep->p4().y(), lep->p4().z(), lep->p4().t())
        );
        bool inserted=false;
        for (SimpleParticleCollection_t::iterator ielo=associatedLeptons.begin(); ielo<associatedLeptons.end(); ielo++){
          if (lep->pt()>(*ielo).second.Pt()){
            inserted=true;
            associatedLeptons.insert(ielo, theLepton);
            break;
          }
        }
        if (!inserted) associatedLeptons.push_back(theLepton);
      }
    }
    myCand.addUserFloat("nExtraLep",nExtraLep);

    // Leptonically decaying WH
    if (debug) cout<<"Extra W"<<endl;
    if (nExtraLep>=1){
      // Take leading-pT lepton to compute fake neutrino
      int nuid = -associatedLeptons.at(0).first + (associatedLeptons.at(0).first>0 ? -1 : +1);

      // Take neutrino momentum from the MET, using a W mass constraint to solve for the z component
      float a = associatedLeptons.at(0).second.X();
      float b = associatedLeptons.at(0).second.Y();
      float c = associatedLeptons.at(0).second.Z();
      float f = associatedLeptons.at(0).second.T();
      TLorentzVector myLep(a, b, c, f);
      float x = PFMET*cos(PFMETPhi);
      float y = PFMET*sin(PFMETPhi);
      float m = 80.399;//PDGHelpers::Wmass;
      float delta = pow(c*(-4*a*a - 4*b*b - 4*c*c + 4*f*f - 4*m*m - 8*a*x - 8*b*y), 2) - 4*(-4*c*c + 4*f*f)*(-a*a*a*a - 2*a*a*b*b - b*b*b*b - 2*a*a*c*c - 2*b*b*c*c - c*c*c*c + 2*a*a*f*f + 2*b*b*f*f + 2*c*c*f*f - f*f*f*f - 2*a*a*m*m - 2*b*b*m*m - 2*c*c*m*m + 2*f*f*m*m - m*m*m*m - 4*a*a*a*x - 4*a*b*b*x - 4*a*c*c*x + 4*a*f*f*x - 4*a*m*m*x - 4*a*a*x*x + 4*f*f*x*x - 4*a*a*b*y - 4*b*b*b*y - 4*b*c*c*y + 4*b*f*f*y - 4*b*m*m*y - 8*a*b*x*y - 4*b*b*y*y + 4*f*f*y*y);

      if (delta>=0.){
        float z1 = (-c*(-4*a*a - 4*b*b - 4*c*c + 4*f*f - 4*m*m - 8*a*x - 8*b*y) - sqrt(delta)) / (2*(-4*c*c + 4*f*f));
        float z2 = (-c*(-4*a*a - 4*b*b - 4*c*c + 4*f*f - 4*m*m - 8*a*x - 8*b*y) + sqrt(delta)) / (2*(-4*c*c + 4*f*f));
        TLorentzVector myNu0(x, y, z1, sqrt(x*x+y*y+z1*z1));
        TLorentzVector myNu1(x, y, z2, sqrt(x*x+y*y+z2*z2));
        associatedLeptons.push_back(SimpleParticle_t(nuid, myNu0));
        associatedLeptons.push_back(SimpleParticle_t(nuid, myNu1));
      }
      else{
        TLorentzVector myNu(x, y, 0, TMath::Sqrt(x*x+y*y));
        associatedLeptons.push_back(SimpleParticle_t(nuid, myNu));
      }
    }


    //--- store Z candidates whose leptons are not involved in the current ZZ candidate
    if (debug) cout<<"Extra Z"<<endl;
    int nExtraZ = 0;
    vector<const CompositeCandidate*> extraZs;
    for( View<CompositeCandidate>::const_iterator zcand = ZCands->begin(); zcand != ZCands->end(); ++ zcand ) {
      if( reco::deltaR( zcand->daughter(0)->p4(), Z1L1->p4() ) > 0.02 &&
      reco::deltaR( zcand->daughter(0)->p4(), Z1L2->p4() ) > 0.02 &&
      reco::deltaR( zcand->daughter(0)->p4(), Z2L1->p4() ) > 0.02 &&
      reco::deltaR( zcand->daughter(0)->p4(), Z2L2->p4() ) > 0.02 &&
      reco::deltaR( zcand->daughter(1)->p4(), Z1L1->p4() ) > 0.02 &&
      reco::deltaR( zcand->daughter(1)->p4(), Z1L2->p4() ) > 0.02 &&
      reco::deltaR( zcand->daughter(1)->p4(), Z2L1->p4() ) > 0.02 &&
      reco::deltaR( zcand->daughter(1)->p4(), Z2L2->p4() ) > 0.02    ){
    const reco::CandidatePtr myZCand(ZCands,zcand-ZCands->begin());
    if((bool)userdatahelpers::getUserFloat(&*myZCand,"GoodIsoLeptons")){
      nExtraZ++;
      extraZs.push_back(&*zcand);
      myCand.addUserCand("assocZ"+to_string(nExtraZ),myZCand);
    }
      }
    }
    myCand.addUserFloat("nExtraZ",nExtraZ);


    // Lepton TLorentzVectors, including FSR
    SimpleParticleCollection_t daughters;
    daughters.push_back(SimpleParticle_t(id11, TLorentzVector(p11.x(), p11.y(), p11.z(), p11.t())));
    daughters.push_back(SimpleParticle_t(id12, TLorentzVector(p12.x(), p12.y(), p12.z(), p12.t())));
    daughters.push_back(SimpleParticle_t(id21, TLorentzVector(p21.x(), p21.y(), p21.z(), p21.t())));
    daughters.push_back(SimpleParticle_t(id22, TLorentzVector(p22.x(), p22.y(), p22.z(), p22.t())));

    if (debug) cout<<"Angles, not computed"<<endl;
    //--- Compute angles, better done here
    float costheta1=0, costheta2=0, phi=0, costhetastar=0, phistar1=0;
    //TUtil::computeAngles(
    //  costhetastar, costheta1, costheta2, phi, phistar1,
    //  daughters.at(0).second, daughters.at(0).first,
    //  daughters.at(1).second, daughters.at(1).first,
    //  daughters.at(2).second, daughters.at(2).first,
    //  daughters.at(3).second, daughters.at(3).first
    //);
    //--- compute higgs azimuthal angles, xi
    TLorentzVector Z14vec = daughters.at(0).second + daughters.at(1).second;
    TLorentzVector higgs = Z14vec + daughters.at(2).second + daughters.at(3).second;
    TVector3 Xaxis(1, 0, 0);
    float xi = higgs.Phi();
    // boost Z1 into rest frame of higgs
    // xistar is the angle between Z decay plane and x-axis
    Z14vec.Boost(-higgs.BoostVector());
    float xistar = Z14vec.Phi();

    // detaJJ, Mjj and Fisher. These are per-event variables in the SR, but not necessarily in the CR as we clean jets also
    // for loose-but-not-tight leptons.
    if (debug) cout<<"Dijet stuff"<<endl;
    float DiJetMass  = -99;
    float DiJetDEta  = -99;
    float DiJetFisher  = -99;
    float ZZjjPt     = -99;

    unsigned int nCandidates=0; // Should equal jecnum after the loop below
    
    // Loop over following JEC variations:
    // JES
    // JER
    if (debug) cout<<"JES, JER"<<endl;
    for (int jecnum = 0; jecnum < _num_of_JEC_variations; jecnum++){
      SimpleParticleCollection_t associated;

      vector<const pat::Jet*> cleanedJetsPt30Jec;
      vector<float> jec_ratio;
      for (edm::View<pat::Jet>::const_iterator jet = CleanJets->begin(); jet != CleanJets->end(); ++jet){
        
        // Nominal jet
        float ratio = 1.;
        float newPt = jet->pt();
        
        // calculate all JEC uncertainty up/down
        
        //JES Up uncertainty
        if (jecnum == 1 )      {
           ratio = 1. + jet->userFloat("jes_unc");
           newPt = jet->pt() * ratio;
        }
        
        //JES Down uncertainty
        else if (jecnum == 2 ) {
           ratio = 1. - jet->userFloat("jes_unc");
           newPt = jet->pt() * ratio;
        }
         
        //JER Up uncertainty
        else if (jecnum == 3 ) {
           ratio = jet->userFloat("pt_jerup") / jet->pt();
           newPt = jet->userFloat("pt_jerup");
        }
         
        //JER Down uncertainty
        else if (jecnum == 4 ) {
           ratio = jet->userFloat("pt_jerdn") / jet->pt();
           newPt = jet->userFloat("pt_jerdn");
        }

        if (newPt<=30.) continue;
        // additional jets cleaning for loose leptons belonging to this candidate (for CRs only;
        // does nothing for the SR as jets are already cleaned with all tight isolated leptons )
        if (!jetCleaner::isGood(myCand, *jet)) continue;
        // store jets and up/down ratio
        cleanedJetsPt30Jec.push_back(&*jet);
        jec_ratio.push_back(ratio);
      }
      if (jecnum==0 && cleanedJetsPt30Jec.size()>1){
        const pat::Jet& jet1 = *(cleanedJetsPt30Jec.at(0));
        const pat::Jet& jet2 = *(cleanedJetsPt30Jec.at(1));
        DiJetDEta = jet1.eta()-jet2.eta();
        DiJetMass = (jet1.p4()+jet2.p4()).M();
        DiJetFisher = fisher(DiJetMass, DiJetDEta);
      }
      
      vector<const pat::Jet*> cleanedJetsPt30;
      for (edm::View<pat::Jet>::const_iterator jet = CleanJets->begin(); jet != CleanJets->end(); ++jet){
        if (jet->pt() > 30.) cleanedJetsPt30.push_back(&*jet);
      }
      
      if(cleanedJetsPt30.size() > 1)
      {
        const pat::Jet& jet1 = *(cleanedJetsPt30.at(0));
        const pat::Jet& jet2 = *(cleanedJetsPt30.at(1));
        ZZjjPt = (Z1Lm->p4()+Z1Lp->p4()+Z2Lm->p4()+Z2Lp->p4()+jet1.p4()+jet2.p4()).pt();
      }
      
      for (unsigned int ijet = 0; ijet<cleanedJetsPt30Jec.size(); ijet++){
        TLorentzVector jet(
          cleanedJetsPt30Jec[ijet]->p4().x()*jec_ratio.at(ijet),
          cleanedJetsPt30Jec[ijet]->p4().y()*jec_ratio.at(ijet),
          cleanedJetsPt30Jec[ijet]->p4().z()*jec_ratio.at(ijet),
          cleanedJetsPt30Jec[ijet]->p4().t()*jec_ratio.at(ijet)
          );
        associated.push_back(SimpleParticle_t(0, jet));
      }
      for (unsigned int ilep=0; ilep<associatedLeptons.size(); ilep++) associated.push_back(associatedLeptons.at(ilep));
      for (unsigned int ijet = 0; ijet<cleanedJetsPt30Jec.size(); ijet++){
        TLorentzVector jet(
          cleanedJetsPt30Jec[ijet]->p4().x()*jec_ratio.at(ijet),
          cleanedJetsPt30Jec[ijet]->p4().y()*jec_ratio.at(ijet),
          cleanedJetsPt30Jec[ijet]->p4().z()*jec_ratio.at(ijet),
	  cleanedJetsPt30Jec[ijet]->p4().t()*jec_ratio.at(ijet)
 	  );
        SimpleParticleCollection_t stableTopDaughters; // Just a collection with one jet as the top
        stableTopDaughters.push_back(SimpleParticle_t(0, jet));
      }
      nCandidates++;
    }

    //----------------------------------------------------------------------
    //--- kinematic fit for ZZ containing taus
    if (debug) cout<<"Kinematic fit"<<endl;
    float ZZKMass = -999;
    float ZZKChi2 = -999;
    if (doKinFit){
    bool wrongFit=false;
    bool checkFlavor=false;
    TLorentzVector p_L1, p_L2, p_Tau1, p_Tau2;
    if ( (abs(id11)==15 or abs(id12)==15) && abs(id21)!=15 && abs(id22)!=15 ){
	checkFlavor=true;
	p_L1.SetPxPyPzE(p21.x(),p21.y(),p21.z(),p21.t());
	p_L2.SetPxPyPzE(p22.x(),p22.y(),p22.z(),p22.t());
	p_Tau1.SetPxPyPzE(p11.x(),p11.y(),p11.z(),p11.t());
	p_Tau2.SetPxPyPzE(p12.x(),p12.y(),p12.z(),p12.t());
    }
    else if ( (abs(id21)==15 or abs(id22)==15) && abs(id11)!=15 && abs(id12)!=15 ){
	checkFlavor=true;
	p_L1.SetPxPyPzE(p11.x(),p11.y(),p11.z(),p11.t());
        p_L2.SetPxPyPzE(p12.x(),p12.y(),p12.z(),p12.t());
        p_Tau1.SetPxPyPzE(p21.x(),p21.y(),p21.z(),p21.t());
        p_Tau2.SetPxPyPzE(p22.x(),p22.y(),p22.z(),p22.t());
    }
    if (checkFlavor){
	HHKinFit2::HHKinFitMasterHeavyHiggs kinFits = HHKinFit2::HHKinFitMasterHeavyHiggs(p_L1, p_L2, p_Tau1, p_Tau2, ptmiss, covMET, 0., 0.);
	for (int i = 0; i<10; i++){
	    double mZ2=(i+1)*0.1*ZmassValue;
	    kinFits.addHypo(ZmassValue,mZ2);
	    kinFits.addHypo(mZ2,ZmassValue);
	}
	try{ kinFits.fit();}
	catch(HHKinFit2::HHInvMConstraintException &e){
	    cout<<"INVME THIS EVENT WAS WRONG, INV MASS CONSTRAIN EXCEPTION"<<endl;
            wrongFit=true;
	}
	catch (HHKinFit2::HHEnergyRangeException &e){
	    cout<<"ERANGE THIS EVENT WAS WRONG, ENERGY RANGE EXCEPTION"<<endl;
	    wrongFit=true;
	}
	catch(HHKinFit2::HHEnergyConstraintException &e){
	    cout<<"ECON THIS EVENT WAS WRONG, ENERGY CONSTRAIN EXCEPTION"<<endl;
	    wrongFit=true;
	}
	if (!wrongFit){
	    ZZKMass=kinFits.getMH();
	    ZZKChi2=kinFits.getChi2();
	}
	else{
	    ZZKMass=-333.;
	    ZZKChi2=-333.;
	}
    }
    else{
        ZZKMass=-666.;
	ZZKChi2=-666.;
    }
    }
    myCand.addUserFloat("ZZKMass",ZZKMass);
    myCand.addUserFloat("ZZKChi2",ZZKChi2);


    //----------------------------------------------------------------------
    //--- Get the mass by combining l1l2 + SVFit momentum
    if (debug) cout<<"SV fit propagation"<<endl;
    bool checkFlavor=false;
    float SVfitMass=-999., SVpt=-999., SVeta=-999., SVphi=-999.;
    TLorentzVector pZtt,pZll,pZZ;
    if ( (abs(id11)==15 or abs(id12)==15) && abs(id21)!=15 && abs(id22)!=15 ){
	if (userdatahelpers::hasUserFloat(Z1,"ComputeSV") && userdatahelpers::getUserFloat(Z1,"ComputeSV")){
	    pZtt.SetPtEtaPhiM(userdatahelpers::getUserFloat(Z1,"SVpt"),userdatahelpers::getUserFloat(Z1,"SVeta"),userdatahelpers::getUserFloat(Z1,"SVphi"),userdatahelpers::getUserFloat(Z1,"SVfitMass"));
	    pZll.SetPtEtaPhiM(Z2->pt(),Z2->eta(),Z2->phi(),Z2->mass());
	    checkFlavor=true;
	}
	else{
	    SVfitMass=-666.;
	    SVpt=-666.;
	    SVeta=-666.;
	    SVphi=-666.;
	}
    }
    else if ( (abs(id21)==15 or abs(id22)==15) && abs(id11)!=15 && abs(id12)!=15 ){
	if (userdatahelpers::hasUserFloat(Z2,"ComputeSV") && userdatahelpers::getUserFloat(Z2,"ComputeSV")){
            pZtt.SetPtEtaPhiM(userdatahelpers::getUserFloat(Z2,"SVpt"),userdatahelpers::getUserFloat(Z2,"SVeta"),userdatahelpers::getUserFloat(Z2,"SVphi"),userdatahelpers::getUserFloat(Z2,"SVfitMass"));
	    pZll.SetPtEtaPhiM(Z1->pt(),Z1->eta(),Z1->phi(),Z1->mass());
	    checkFlavor=true;
        }
        else{
            SVfitMass=-666.;
            SVpt=-666.;
            SVeta=-666.;
            SVphi=-666.;
        }
    }
    else{
	SVfitMass=-333.;
        SVpt=-333.;
        SVeta=-333.;
        SVphi=-333.;
    }
    if (checkFlavor){
	pZZ=pZtt+pZll;
	SVfitMass=pZZ.M();
	SVpt=pZZ.Pt();
	SVeta=pZZ.Eta();
	SVphi=pZZ.Phi();
    }
    
    myCand.addUserFloat("SVfitMass",SVfitMass);
    myCand.addUserFloat("SVpt",SVpt);
    myCand.addUserFloat("SVeta",SVeta);
    myCand.addUserFloat("SVphi",SVphi);


    //----------------------------------------------------------------------
    //--- Decide which 4l mass to be used
    if (debug) cout<<"Good mass"<<endl;
    float goodMass;
    if (SVfitMass < 0)
	goodMass=myCand.mass();
    else
	goodMass=SVfitMass;//ZZKMass;
    myCand.addUserFloat("goodMass",goodMass);
	
    //----------------------------------------------------------------------
    //--- kinematic refitting using Z mass constraint
    if (debug) cout<<"Kinamtic refitting"<<endl;
    float ZZMassRefit = 0.;
    float ZZMassRefitErr = 0.;
    float ZZMassUnrefitErr = 0.;

    if(doKinFitOld){

      vector<reco::Candidate *> selectedLeptons;
      std::map<unsigned int, TLorentzVector> selectedFsrMap;

      for(unsigned ilep=0; ilep<4; ilep++){

    selectedLeptons.push_back((reco::Candidate*)(ZZLeps[ilep]->masterClone().get()));

    if(FSRMap.find(ZZLeps[ilep])!=FSRMap.end()){
      pat::PFParticle fsr = *(FSRMap[ZZLeps[ilep]]);
      TLorentzVector p4;
      p4.SetPxPyPzE(fsr.px(),fsr.py(),fsr.pz(),fsr.energy());
      selectedFsrMap[ilep] = p4;
    }

      }

      kinZfitter->Setup(selectedLeptons, selectedFsrMap);
      kinZfitter->KinRefitZ();

      ZZMassRefit = kinZfitter->GetRefitM4l();
      ZZMassRefitErr = kinZfitter->GetRefitM4lErrFullCov();
      ZZMassUnrefitErr = kinZfitter->GetM4lErr();

      // four 4-vectors after refitting order by Z1_1,Z1_2,Z2_1,Z2_2
      //vector<TLorentzVector> p4 = kinZfitter->GetRefitP4s();

    }


    //----------------------------------------------------------------------
    //--- 4l vertex fits (experimental)

    //CandConstraintFit::fit(&myCand, iSetup);
    if (debug) cout<<"Vertex fit"<<endl;
    if (doVtxFit && abs(id11)==13 && abs(id12)==13 && abs(id21)==13 && abs(id22)==13) {

      edm::ESHandle<TransientTrackBuilder> theTTBuilder;
      iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);

      //Creating a KinematicParticleFactory
      KinematicParticleFactoryFromTransientTrack factory;

      const ParticleMass muon_mass = 0.1056583;
      const ParticleMass electron_mass = 0.0005;
      float muon_sigma = 0.0000000001;
      float electron_sigma = 0.0000000001;

      //initial chi2 and ndf before kinematic fits. The chi2 of the reconstruction is not considered
      float chi;
      float ndof;

      vector<RefCountedKinematicParticle> Particles;
      //vector<TransientTrack> t_tks;

      for (unsigned k = 0; k < myCand.numberOfDaughters(); ++k ) {
        const reco::Candidate* Z = myCand.daughter(k);
        for (unsigned l = 0; l < Z->numberOfDaughters(); ++l ) {
          chi = 0.; ndof = 0.;

          const reco::Candidate* lepton= Z->daughter(l);

          if (lepton->isGlobalMuon() || lepton->isTrackerMuon()){
            TransientTrack tt = theTTBuilder->build(lepton->get<TrackRef>());
            Particles.push_back(factory.particle (tt,muon_mass,chi,ndof,muon_sigma));
            //t_tks.push_back(tt);
          }
          else if (lepton->isElectron()){
            TransientTrack tt = theTTBuilder->build(lepton->get<GsfTrackRef>());
            Particles.push_back(factory.particle (tt,electron_mass,chi,ndof,electron_sigma));
            //t_tks.push_back(tt);
          }
        }
      }

      //cout << "Number of particle for constrain fitter= " << Particles.size()<< endl;

      if (Particles.size()>=4){
        KinematicParticleVertexFitter fitter;
        RefCountedKinematicTree myTree = fitter.fit(Particles);

        if ( !myTree->isEmpty()) {
          //accessing the tree components
          myTree->movePointerToTheTop();

          RefCountedKinematicParticle allLeptonsCand     = myTree->currentParticle();
          RefCountedKinematicVertex allLeptonsVertex     = myTree->currentDecayVertex();

          // if(dbg) cout << "m(" << myCand->numberOfDaughters() << "l): " << allLeptonsCand->currentState().mass() << " +- "
          //                   << sqrt(allLeptonsCand->currentState().kinematicParametersError().matrix()(6,6) ) << endl;

          reco::Vertex constrainedVertex(reco::Vertex::Point(allLeptonsVertex->position()),
                                         allLeptonsVertex->error().matrix(),
                                         allLeptonsVertex->chiSquared(),
                                         allLeptonsVertex->degreesOfFreedom(),0);

          //      if(dbg) cout << "kinematicFit vertex, ndof, chi2, prob: "
          //                   << allLeptonsVertex->position() << " , "
          //                   << allLeptonsVertex->degreesOfFreedom() << " , "
          //                   << allLeptonsVertex->chiSquared()   << " , "
          //                   << TMath::Prob(allLeptonsVertex->chiSquared(),allLeptonsVertex->degreesOfFreedom()) << endl;

          //myCand->addUserData("ConstrainedCandVtx",constrainedVertex);
          myCand.addUserFloat("CFitM",allLeptonsCand->currentState().mass());
          myCand.addUserFloat("CFitSigmaM",allLeptonsCand->currentState().kinematicParametersError().matrix()(6,6));
          myCand.addUserFloat("CFitNdof",allLeptonsVertex->degreesOfFreedom());
          myCand.addUserFloat("CFitChi2",allLeptonsVertex->chiSquared());

        } else {
          cout << " ERROR CandConstraintFit: KinematicParticleVertexFitter failed " << endl;
        }
      }
    }
    
    if (debug) cout<<"Good tau and HLT match"<<endl;
    bool goodTaus=false;
    if (myCand.userFloat("d0.isGoodTau") && myCand.userFloat("d1.isGoodTau"))
	goodTaus=true;
    
    bool muHLTMatch=false;
    bool eleHLTMatch=false;
    if (myCand.userFloat("d0.muHLTMatch") || myCand.userFloat("d1.muHLTMatch"))
	muHLTMatch=true;
    if (myCand.userFloat("d0.eleHLTMatch") || myCand.userFloat("d1.eleHLTMatch"))
	eleHLTMatch=true;

    
    //----------------------------------------------------------------------
    //--- Embed variables
    if (debug) cout<<"Embed variables"<<endl;
    myCand.addUserFloat("candChannel",    candChannel);
    myCand.addUserFloat("SIP4",           SIP4);
    myCand.addUserFloat("pt1",            ptS.at(3).first); // leading-pT
    myCand.addUserFloat("pdgId1",	  ptS.at(3).second);
    myCand.addUserFloat("pt2",            ptS.at(2).first); // sub-leading pT
    myCand.addUserFloat("pdgId2",         ptS.at(2).second); // sub-leading pT
    myCand.addUserFloat("mZa",            mZa);
    myCand.addUserFloat("mZb",            mZb);
    myCand.addUserFloat("ZaID",           ZaID);
    myCand.addUserFloat("ZbID",           ZbID);
    myCand.addUserFloat("mZalpha",        mZalpha);
    myCand.addUserFloat("mZbeta",         mZbeta);
    myCand.addUserFloat("ZalphaID",       ZalphaID);
    myCand.addUserFloat("ZbetaID",        ZbetaID);
    myCand.addUserFloat("mLL4",           mll4); // smallest mass of any AF/OS pair
    myCand.addUserFloat("mLL6",           mll6);   // smallest mass of any AF/AS pair
    myCand.addUserFloat("passSmartMLL",   passSmartMLL);
    myCand.addUserFloat("costheta1",      costheta1);
    myCand.addUserFloat("costheta2",      costheta2);
    myCand.addUserFloat("phi",            phi);
    myCand.addUserFloat("costhetastar",   costhetastar);
    myCand.addUserFloat("phistar1",       phistar1);
    myCand.addUserFloat("xistar",         xistar);  //azimuthal angle of higgs in rest frame of higgs
    myCand.addUserFloat("xi",             xi);      //azimuthal angle of higgs in lab frame

    myCand.addUserFloat("m4l",            (Z1Lm->p4()+Z1Lp->p4()+Z2Lm->p4()+Z2Lp->p4()).mass()); // mass without FSR
    if(doKinFit) {
      myCand.addUserFloat("ZZMassRefit"   , ZZMassRefit);
      myCand.addUserFloat("ZZMassRefitErr", ZZMassRefitErr);
      myCand.addUserFloat("ZZMassUnrefitErr", ZZMassUnrefitErr);
    }

    myCand.addUserFloat("goodTaus", goodTaus);
    myCand.addUserFloat("muHLTMatch", muHLTMatch);
    myCand.addUserFloat("eleHLTMatch", eleHLTMatch);

    // Jet quantities
    myCand.addUserFloat("DiJetMass", DiJetMass);
    myCand.addUserFloat("DiJetDEta", DiJetDEta);
    myCand.addUserFloat("DiJetFisher", DiJetFisher);
    
    myCand.addUserFloat("ZZjjPt", ZZjjPt);
    
    // MET
    myCand.addUserFloat("MET", PFMET);
    myCand.addUserFloat("METPhi", PFMETPhi);


    //--- MC matching. To be revised, cf. MuFiller, EleFiller
//     if (isMC) {
//       int refID = 25; // FIXME: handle ZZ (sigId = 23)
//       bool MC_isRight = (myCand.userFloat("d0.d0.MCParentCode")==refID &&
//                       myCand.userFloat("d0.d1.MCParentCode")==refID &&
//                       myCand.userFloat("d1.d0.MCParentCode")==refID &&
//                       myCand.userFloat("d1.d1.MCParentCode")==refID);
//       bool MC_isRightPair = false; //FIXME to be

//       myCand.addUserFloat("MC_isRight",     MC_isRight);
//       myCand.addUserFloat("MC_isRightPair", MC_isRightPair);
//     }



    //----------------------------------------------------------------------
    //--- Check if candedate passes the "bestCandAmong" selections (2011 PRL logic)
    if (debug) cout<<"Best cand among"<<endl;
    int iCRname=0;
    for(CutSet<pat::CompositeCandidate>::const_iterator bca = preBestCandSelection.begin(); bca != preBestCandSelection.end(); ++bca){
      if (debug) cout<<bca->first<<endl;
      int preBestCandResult= int((*(bca->second))(myCand));

      if (preBestCandResult){
        // Fill preSelCands matrix
        preSelCands[iCRname].push_back(icand);
      }
      iCRname++;
    }
    result->push_back(myCand);

  } // End of loop over input candidates


  //--- For each of the bestCandAmong preselections, find the best candidate and store its index (bestCandIdx)
  if (debug) cout<<"Best candidate"<<endl;
  Comparators::BestCandComparator myComp(*result, bestCandType);
  for (int iCRname=0; iCRname<(int)preSelCands.size(); ++iCRname) {
    cout<<preSelCands[iCRname].size()<<" bestCandAmong"<<endl;
    if (preSelCands[iCRname].size() > 0) {
      bestCandIdx[iCRname] = *std::min_element( preSelCands[iCRname].begin(), preSelCands[iCRname].end(), myComp);
      cout<<"iCRname "<<iCRname<<", bestCandIdx "<<bestCandIdx[iCRname]<<endl;
    }
  }

  //--- Embed best candidate flag (must be done in a separate loop)
  if (debug) cout<<"Embed flags"<<endl;
  for (int i = 0; i< (int)result->size(); ++i) {
    pat::CompositeCandidate& myCand = (*result)[i];

    int iCRname=0;
    for(CutSet<pat::CompositeCandidate>::const_iterator bca = preBestCandSelection.begin(); bca != preBestCandSelection.end(); ++bca){
      bool isBestCand = (i==bestCandIdx[iCRname]);
      myCand.addUserFloat(bca->first,isBestCand);
      iCRname++;
    }

    //--- Embed flags (ie cuts specified in the "flags" pset).
    //    We do this here so that isBestCand is available within the cuts.
    for(CutSet<pat::CompositeCandidate>::const_iterator cut = cuts.begin(); cut != cuts.end(); ++cut) {
      myCand.addUserFloat(cut->first,int((*(cut->second))(myCand)));
    }
  }
  cout<<result->size()<<" ZZ candidates"<<endl;
  iEvent.put(std::move(result));

}


void
ZZCandidateFiller::getPairMass(const reco::Candidate* lp, const reco::Candidate* lm, ZZCandidateFiller::FSRToLepMap& photons, float& mass, int& ID){
  math::XYZTLorentzVector llp4 = lp->p4()+lm->p4();
  auto lpp = photons.find(lp);
  auto lmp = photons.find(lm);
  if (lpp!=photons.end()) llp4+=lpp->second->p4();
  if (lmp!=photons.end()) llp4+=lmp->second->p4();
  mass=llp4.mass();
  ID=lp->pdgId()*lm->pdgId();
}


#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(ZZCandidateFiller);

