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
  bool triggerFlavour(pat::TriggerObjectStandAlone OBJ, int flavour);

  edm::EDGetTokenT<edm::View<reco::CompositeCandidate> > candidateToken;
  const StringCutObjectSelector<pat::CompositeCandidate, true> preBestZSelection;
  const edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  const edm::EDGetTokenT< edm::TriggerResults > triggerResultsToken_;
  edm::EDGetTokenT<edm::View<reco::Candidate> > softLeptonToken;

  int sampleType;
  int setup;
  const CutSet<pat::CompositeCandidate> cuts;
  int FSRMode;
  bool embedDaughterFloats;
  
  //For HLT match
  vector<string> SingleMuPaths, SingleMuFilters;
  vector<string> SingleElePaths, SingleEleFilters;
    
  vector<string> DiMuPaths, DiMuFilters, DiMuFiltersLeg1, DiMuFiltersLeg2; vector<float> DiMuFiltersDZ, DiMuFiltersDR;
  vector<string> DiTauPaths, DiTauFilters, DiTauFiltersLeg1, DiTauFiltersLeg2; vector<float> DiTauFiltersDZ, DiTauFiltersDR;
  vector<string> MuTauPaths, MuTauFilters, MuTauFiltersLeg1, MuTauFiltersLeg2; vector<float> MuTauFiltersDZ, MuTauFiltersDR;
  vector<string> EleTauPaths, EleTauFilters, EleTauFiltersLeg1, EleTauFiltersLeg2; vector<float> EleTauFiltersDZ, EleTauFiltersDR;
  vector<string> MuElePaths, MuEleFilters, MuEleFiltersLeg1, MuEleFiltersLeg2; vector<float> MuEleFiltersDZ, MuEleFiltersDR;


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
  triggerResultsToken_( consumes< edm::TriggerResults >( iConfig.getParameter< edm::InputTag >( "TriggerResultsLabel" ) ) ),
  sampleType(iConfig.getParameter<int>("sampleType")),
  setup(iConfig.getParameter<int>("setup")),
  cuts(iConfig.getParameter<edm::ParameterSet>("flags")),
  embedDaughterFloats(iConfig.getUntrackedParameter<bool>("embedDaughterFloats",true))
{
  
  string mode = iConfig.getParameter<string>("FSRMode");
  if      (mode == "skip")   FSRMode = 0;
  // "Legacy" Run I mode (1) is no longer supported.
  else if (mode == "RunII")  FSRMode = 2;
  else {
    cout << "ZCandidateFiller: FSRMode " << FSRMode << " not supported" << endl;
    abort();
  }
      
  softLeptonToken = consumes<edm::View<reco::Candidate> >(edm::InputTag("softLeptons"));

  if (sampleType == 2016)
  {
      //SingleMu
      SingleMuPaths = {
          "HLT_IsoMu22_v*",
          "HLT_IsoMu22_eta2p1_v*",
          "HLT_IsoTkMu22_v*",
          "HLT_IsoTkMu22_eta2p1_v*",
          "HLT_IsoMu24_v*",
          "HLT_IsoTkMu24_v*"
      };
      SingleMuFilters = {
          "hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09",
          "hltL3crIsoL1sSingleMu20erL1f0L2f10QL3f22QL3trkIsoFiltered0p09",
          "hltL3fL1sMu20L1f0Tkf22QL3trkIsoFiltered0p09",
          "hltL3fL1sMu20erL1f0Tkf22QL3trkIsoFiltered0p09",
          "hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09",
          "hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09"
      };
      //SingleEle
      SingleElePaths = {
          "HLT_Ele25_eta2p1_WPTight_Gsf_v*",
          "HLT_Ele27_WPTight_Gsf_v*",
          "HLT_Ele27_eta2p1_WPLoose_Gsf_v*",
          "HLT_Ele32_eta2p1_WPTight_Gsf_v*"
      };
      SingleEleFilters = {
          "hltEle25erWPTightGsfTrackIsoFilter",
          "hltEle27WPTightGsfTrackIsoFilter",
          "hltEle27erWPLooseGsfTrackIsoFilter",
          "hltEle32WPTightGsfTrackIsoFilter",
      };
      //DiMu
      DiMuPaths = {
          "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",
          "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*",
          "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*",
          "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*"
      };
      DiMuFilters = {
          "hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4DzFiltered0p2",
          "hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4DzFiltered0p2",
          "hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4",
          "hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4"
      };
      DiMuFiltersLeg1 = {
          "hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4",
          "hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4",
          "hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4",
          "hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4"

      };
      DiMuFiltersLeg2 = {
          "hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4",
          "hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4",
          "hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4",
          "hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4"
      };
      DiMuFiltersDZ = {0.2,0.2,-1,-1};
      DiMuFiltersDR = {0.001,0.001,-1,-1};
      //DiTau
      DiTauPaths = {
          "HLT_DoubleMediumIsoPFTau32_Trk1_eta2p1_Reg_v*",
          "HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v*",
          "HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v*",
          "HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v*",
          "HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v*",
          "HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v*"
      };
      DiTauFilters = {
          "hltDoublePFTau32TrackPt1MediumIsolationDz02Reg",
          "hltDoublePFTau35TrackPt1MediumIsolationDz02Reg",
          "hltDoublePFTau40TrackPt1MediumIsolationDz02Reg",
          "hltDoublePFTau35TrackPt1MediumCombinedIsolationDz02Reg",
          "hltDoublePFTau40TrackPt1MediumCombinedIsolationDz02Reg",
          "hltDoublePFTau40TrackPt1MediumCombinedIsolationDz02"
      };
      DiTauFiltersLeg1 = {
          "hltDoublePFTau32TrackPt1MediumIsolationDz02Reg",
          "hltDoublePFTau35TrackPt1MediumIsolationDz02Reg",
          "hltDoublePFTau40TrackPt1MediumIsolationDz02Reg",
          "hltDoublePFTau35TrackPt1MediumCombinedIsolationDz02Reg",
          "hltDoublePFTau40TrackPt1MediumCombinedIsolationDz02Reg",
          "hltDoublePFTau40TrackPt1MediumCombinedIsolationDz02"
      };
      DiTauFiltersLeg2 = {
          "hltDoublePFTau32TrackPt1MediumIsolationDz02Reg",
          "hltDoublePFTau35TrackPt1MediumIsolationDz02Reg",
          "hltDoublePFTau40TrackPt1MediumIsolationDz02Reg",
          "hltDoublePFTau35TrackPt1MediumCombinedIsolationDz02Reg",
          "hltDoublePFTau40TrackPt1MediumCombinedIsolationDz02Reg",
          "hltDoublePFTau40TrackPt1MediumCombinedIsolationDz02"
      };
      DiTauFiltersDZ = {0.2,0.2,0.2,0.2,0.2,0.2};
      DiTauFiltersDR = {0.5,0.5,0.5,0.5,0.5,0.5};
      //MuTau
      MuTauPaths = {
          "HLT_IsoMu19_eta2p1_LooseIsoPFTau20_v*",
          "HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_v*"
      };
      MuTauFilters = {
          "hltOverlapFilterIsoMu19LooseIsoPFTau20",
          "hltOverlapFilterSingleIsoMu19LooseIsoPFTau20"
      };
      MuTauFiltersLeg1 = {//Mu
          "hltL3crIsoL1sMu18erTauJet20erL1f0L2f10QL3f19QL3trkIsoFiltered0p09",
          "hltL3crIsoL1sSingleMu18erIorSingleMu20erL1f0L2f10QL3f19QL3trkIsoFiltered0p09",
      };
      MuTauFiltersLeg2 = {//Tau
          "hltPFTau20TrackLooseIsoAgainstMuon",
          "hltPFTau20TrackLooseIsoAgainstMuon"
      };
      MuTauFiltersDZ = {-1,-1};
      MuTauFiltersDR = {0.3,0.3};
      //EleTau
      EleTauPaths = {};
      EleTauFilters = {};
      EleTauFiltersLeg1 = {};
      EleTauFiltersLeg2 = {};
      EleTauFiltersDZ = {};
      EleTauFiltersDR = {};
      //MuEle
      MuElePaths = {
          "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*",
          "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*",
          "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
          "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*"
      };
      MuEleFilters = {
          "hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter",
          "hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLDZFilter",
          "hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter",
          "hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLDZFilter"
      };
      MuEleFiltersLeg1 = {//Ele
          "hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegDphiFilter",
          "hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter",
          "hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegDphiFilter",
          "hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"
      };
      MuEleFiltersLeg2 = {//Mu
          "hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegDphiFilter",
          "hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8",
          "hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegDphiFilter",
          "hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23"
      };
      MuEleFiltersDZ = {-1,0.2,-1,0.2};
      MuEleFiltersDR = {-1,-1,-1,-1};
  }
  else if (sampleType == 2017)
  {
      //SingleMu
      SingleMuPaths = {
          "HLT_IsoMu24_v*",
          "HLT_IsoMu27_v*"
      };
      SingleMuFilters = {
          "hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07",
          "hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07"
      };
      //SingleEle
      SingleElePaths = {
          "HLT_Ele32_WPTight_Gsf_v*",
          "HLT_Ele35_WPTight_Gsf_v*"
      };
      SingleEleFilters = {
          "hltEle32WPTightGsfTrackIsoFilter",
          "hltEle35noerWPTightGsfTrackIsoFilter"
      };
      //DiMu
      DiMuPaths = {
          "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*",
          "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*"
      };
      DiMuFilters = {
          "hltDiMuon178Mass3p8Filtered",
          "hltDiMuon178Mass8Filtered",
      };
      DiMuFiltersLeg1 = {
          "hltDiMuon178RelTrkIsoFiltered0p4DzFiltered0p2",
          "hltDiMuon178RelTrkIsoFiltered0p4DzFiltered0p2"
      };
      DiMuFiltersLeg2 = {
          "hltDiMuon178RelTrkIsoFiltered0p4DzFiltered0p2",
          "hltDiMuon178RelTrkIsoFiltered0p4DzFiltered0p2"
      };
      DiMuFiltersDZ = {0.2,0.2};
      DiMuFiltersDR = {0.001,0.001};
      //DiTau
      DiTauPaths = {
          "HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v*",
          "HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v*",
          "HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v*"
      };
      DiTauFilters = {
          "hltDoublePFTau40TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg ",
          "hltDoublePFTau40TrackPt1TightChargedIsolationDz02Reg",
          "hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg "
      };
      DiTauFiltersLeg1 = {
          "hltL1JetsHLTDoublePFTauTrackPt1MediumChargedIsolationAndTightOOSCPhotonsMatchReg",
          "hltL1JetsHLTDoublePFTauTrackPt1TightChargedIsolationMatchReg",
          "hltL1JetsHLTDoublePFTauTrackPt1TightChargedIsolationAndTightOOSCPhotonsMatchReg"
      };
      DiTauFiltersLeg2 = {
          "hltL1JetsHLTDoublePFTauTrackPt1MediumChargedIsolationAndTightOOSCPhotonsMatchReg",
          "hltL1JetsHLTDoublePFTauTrackPt1TightChargedIsolationMatchReg",
          "hltL1JetsHLTDoublePFTauTrackPt1TightChargedIsolationAndTightOOSCPhotonsMatchReg"
      };
      DiTauFiltersDZ = {0.2,0.2,0.2};
      DiTauFiltersDR = {0.5,0.5,0.5};
      //MuTau
      MuTauPaths = {
          "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v*"
      };
      MuTauFilters = {
          "hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded"
      };
      MuTauFiltersLeg1 = {//Mu
          "hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07"
      };
      MuTauFiltersLeg2 = {//Tau
          "hltSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched"
      };
      MuTauFiltersDZ = {-1};
      MuTauFiltersDR = {0.3};
      //EleTau
      EleTauPaths = {
          "HLT_Ele24_eta2p1_WPTight_Gsf_Loose_ChargedIsoPFTau30_eta2p1_CrossL1_v*"
      };
      EleTauFilters = {
          "hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30"
      };
      EleTauFiltersLeg1 = {//Ele
          "hltEle24erWPTightGsfTrackIsoFilterForTau"
      };
      EleTauFiltersLeg2 = {
          "hltSelectedPFTau30LooseChargedIsolationL1HLTMatched"
      };
      EleTauFiltersDZ = {-1};
      EleTauFiltersDR = {0.3};
      //MuEle
      MuElePaths = {
          "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*",
          "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*"
      };
      MuEleFilters = {
          "hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLDZFilter",
          "hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLDZFilter"
      };
      MuEleFiltersLeg1 = {//Ele
          "hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter",
          "hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"
      };
      MuEleFiltersLeg2 = {//Mu
          "hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8",
          "hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23"
      };
      MuEleFiltersDZ = {0.2,0.2};
      MuEleFiltersDR = {-1,-1};
  }
  else if (sampleType == 2018)
  {
      //SingleMu
      SingleMuPaths = {
          "HLT_IsoMu24_v*",
          "HLT_IsoMu27_v*"
      };
      SingleMuFilters = {
          "hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07",
          "hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07"
      };
      //SingleEle
      SingleElePaths = {
          "HLT_Ele32_WPTight_Gsf_v*",
          "HLT_Ele35_WPTight_Gsf_v*"
      };
      SingleEleFilters = {
          "hltEle32WPTightGsfTrackIsoFilter",
          "hltEle35noerWPTightGsfTrackIsoFilter"
      };
      //DiMu
      DiMuPaths = {
          "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*"
      };
      DiMuFilters = {
          "hltDiMuon178Mass3p8Filtered"
      };
      DiMuFiltersLeg1 = {
          "hltDiMuon178RelTrkIsoFiltered0p4DzFiltered0p2"
      };
      DiMuFiltersLeg2 = {
          "hltDiMuon178RelTrkIsoFiltered0p4DzFiltered0p2"
      };
      DiMuFiltersDZ = {0.2};
      DiMuFiltersDR = {0.001};
      //DiTau
      DiTauPaths = {
          "HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v*",
          "HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v*",
          "HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v*",
          "HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_eta2p1_Reg_v*"
      };
      DiTauFilters = {
          "hltDoublePFTau40TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg",
          "hltHpsDoublePFTau40TrackPt1TightChargedIsolationDz02Reg",
          "hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg",
          "hltHpsDoublePFTau40TrackPt1TightChargedIsolationDz02Reg"
      };
      DiTauFiltersLeg1 = {
          "hltL1JetsHLTDoublePFTauTrackPt1MediumChargedIsolationAndTightOOSCPhotonsMatchReg",
          "hltHpsL1JetsHLTDoublePFTauTrackPt1TightChargedIsolationMatchReg",
          "hltL1JetsHLTDoublePFTauTrackPt1TightChargedIsolationAndTightOOSCPhotonsMatchReg",
          "hltHpsL1JetsHLTDoublePFTauTrackPt1TightChargedIsolationMatchReg"
      };
      DiTauFiltersLeg2 = {
          "hltL1JetsHLTDoublePFTauTrackPt1MediumChargedIsolationAndTightOOSCPhotonsMatchReg",
          "hltHpsL1JetsHLTDoublePFTauTrackPt1TightChargedIsolationMatchReg",
          "hltL1JetsHLTDoublePFTauTrackPt1TightChargedIsolationAndTightOOSCPhotonsMatchReg",
          "hltHpsL1JetsHLTDoublePFTauTrackPt1TightChargedIsolationMatchReg"
      };
      DiTauFiltersDZ = {0.2,0.2,0.2,0.2};
      DiTauFiltersDR = {0.5,0.5,0.5,0.5};
      //MuTau
      MuTauPaths = {
          "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_v*",
          "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v*"
      };
      MuTauFilters = {
          "hltHpsOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded",
          "hltHpsOverlapFilterIsoMu20LooseChargedIsoTightOOSCPhotonsPFTau27L1Seeded"
      };
      MuTauFiltersLeg1 = {//Mu
          "hltL3crIsoBigORMu18erTauXXer2p1L1f0L2f10QL3f20QL3trkIsoFiltered0p07",
          "hltL3crIsoBigORMu18erTauXXer2p1L1f0L2f10QL3f20QL3trkIsoFiltered0p07"
      };
      MuTauFiltersLeg2 = {//Tau
          "hltHpsSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched",
          "hltHpsSelectedPFTau27LooseChargedIsolationTightOOSCPhotonsAgainstMuonL1HLTMatched"
      };
      MuTauFiltersDZ = {-1,-1};
      MuTauFiltersDR = {0.3,0.3};
      //EleTau
      EleTauPaths = {
          "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_v*",
          "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_v*"
      };
      EleTauFilters = {
          "hltHpsOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30",
          "hltHpsOverlapFilterIsoEle24WPTightGsfLooseIsoTightOOSCPhotonsPFTau30"
      };
      EleTauFiltersLeg1 = {//Ele
          "hltEle24erWPTightGsfTrackIsoFilterForTau",
          "hltEle24erWPTightGsfTrackIsoFilterForTau"
      };
      EleTauFiltersLeg2 = {
          "hltHpsSelectedPFTau30LooseChargedIsolationL1HLTMatched",
          "hltHpsSelectedPFTau30LooseChargedIsolationTightOOSCPhotonsL1HLTMatched"
      };
      EleTauFiltersDZ = {-1,-1};
      EleTauFiltersDR = {0.3,0.3};
      //MuEle
      MuElePaths = {
          "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
          "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*"
      };
      MuEleFilters = {
          "hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLDZFilter",
          "hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLDZFilter"
      };
      MuEleFiltersLeg1 = {//Ele
          "hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter",
          "hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"
      };
      MuEleFiltersLeg2 = {//Mu
          "hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23",
          "hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8"
      };
      MuEleFiltersDZ = {0.2,0.2};
      MuEleFiltersDR = {-1,-1};
  }
  
  produces<pat::CompositeCandidateCollection>();
}


void ZCandidateFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  Handle<View<reco::Candidate> > softleptoncoll;
  iEvent.getByToken(softLeptonToken, softleptoncoll);
  vector<reco::CandidatePtr> goodisoleptonPtrs;
  for( View<reco::Candidate>::const_iterator lep = softleptoncoll->begin(); lep != softleptoncoll->end(); ++ lep ){
    if((bool)userdatahelpers::getUserFloat(&*lep,"isGood")
       && (bool)userdatahelpers::getUserFloat(&*lep,"passCombRelIsoPFFSRCorr") // FSR-corrected iso
       && abs(lep->pdgId())!=15
       ){
      const reco::CandidatePtr lepPtr(softleptoncoll,lep-softleptoncoll->begin());
      goodisoleptonPtrs.push_back(lepPtr);
    }
  }
    
  int bestZidx = -1;
  int tau1Iso = -1;
  float tau1Pt = -1;
  int tau2Iso = -1;
  float tau2Pt = -1;

  //--- Loop over LL Candidates
  for(unsigned int i = 0; i < LLCands->size(); ++i) {
    const CompositeCandidate& c = (*LLCands)[i];
    pat::CompositeCandidate myCand(c); 

    if (embedDaughterFloats){  
      userdatahelpers::embedDaughterData(myCand);
    }

    int id0 = myCand.daughter(0)->pdgId();
    int id1 = myCand.daughter(1)->pdgId();

    // ------------------------------
    // FSR recovery
    // ------------------------------
    if (FSRMode==2) { // Run II
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

    myCand.addUserFloat("DR",deltaR(*l1,*l2));

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

    // cout<<"HLTMatch: "<<id0*id1<<endl;
    //--------------------------------------------------------------
    //-------------------------HLTMatch-----------------------------
    //--------------------------------------------------------------

    // const edm::TriggerNames &names = iEvent.triggerNames(*triggerResults);

    vector<bool> HLTMatch_singleLeg1, HLTMatch_singleLeg2, HLTMatch_cross;
    vector<string> SinglePaths,SingleFilters,DiPaths,DiFilters,DiFiltersLeg1,DiFiltersLeg2;
    vector<float> DZ,DR;
    
    if (abs(id0*id1)==169) {//mumu
        SinglePaths=SingleMuPaths; SingleFilters=SingleMuFilters;
        DiPaths=DiMuPaths; DiFilters=DiMuFilters; DiFiltersLeg1=DiMuFiltersLeg1; DiFiltersLeg2=DiMuFiltersLeg2; 
        DZ=DiMuFiltersDZ; DR=DiMuFiltersDR;
    }
    else if (abs(id0*id1)==165) {//etau
        SinglePaths=SingleElePaths; SingleFilters=SingleEleFilters;
        DiPaths=EleTauPaths; DiFilters=EleTauFilters; DiFiltersLeg1=EleTauFiltersLeg1; DiFiltersLeg2=EleTauFiltersLeg2; 
        DZ=EleTauFiltersDZ; DR=EleTauFiltersDR;
    }
    else if (abs(id0*id1)==195) {//mutau
        SinglePaths=SingleMuPaths; SingleFilters=SingleMuFilters;
        DiPaths=MuTauPaths; DiFilters=MuTauFilters; DiFiltersLeg1=MuTauFiltersLeg1; DiFiltersLeg2=MuTauFiltersLeg2; 
        DZ=MuTauFiltersDZ; DR=MuTauFiltersDR;
    }
    else if (abs(id0*id1)==225) {//tautau
        DiPaths=DiTauPaths; DiFilters=DiTauFilters; DiFiltersLeg1=DiTauFiltersLeg1; DiFiltersLeg2=DiTauFiltersLeg2; 
        DZ=DiTauFiltersDZ; DR=DiTauFiltersDR;
    }
    else if (abs(id0*id1)==143) {//emu
        SinglePaths.clear();
        SinglePaths.insert(SinglePaths.end(),SingleMuPaths.begin(),SingleMuPaths.end());
        SinglePaths.insert(SinglePaths.end(),SingleElePaths.begin(),SingleElePaths.end());
        SingleFilters.clear();
        SingleFilters.insert(SingleFilters.end(),SingleMuFilters.begin(),SingleMuFilters.end());
        SingleFilters.insert(SingleFilters.end(),SingleEleFilters.begin(),SingleEleFilters.end());
        DiPaths=MuElePaths; DiFilters=MuEleFilters; DiFiltersLeg1=MuEleFiltersLeg1; DiFiltersLeg2=MuEleFiltersLeg2; 
        DZ=MuEleFiltersDZ; DR=MuEleFiltersDR;
    }
    else {
        cout<<"Wrong flavour: "<<id0*id1<<endl;
        continue;
    }
    
    for (size_t itrg=0; itrg<SinglePaths.size();itrg++) {
        HLTMatch_singleLeg1.push_back(false);
        HLTMatch_singleLeg2.push_back(false);
    }
    
    // cout<<"step 1"<<endl;
    for (size_t idxto = 0; idxto < triggerObjects->size(); ++idxto) {//single triggers
        pat::TriggerObjectStandAlone obj = triggerObjects->at(idxto);
        obj.unpackFilterLabels(iEvent,*triggerResults);
        // obj.unpackPathNames(names);
        // if (id0*id1==-165) cout<<SinglePaths.size()<<endl;
        if (deltaR2(*l1,obj)<0.25 && triggerFlavour(obj,id0) && SinglePaths.size() > 0) {
            // if (id0*id1==-165) cout<<SinglePaths.size()<<endl;
            for (size_t itrg=0; itrg<SinglePaths.size();itrg++) {
                // if (id0*id1==-165) cout<<SingleFilters[itrg]<<endl;
                if (obj.hasFilterLabel(SingleFilters[itrg])) HLTMatch_singleLeg1[itrg]=true;
            }
        }
        if (deltaR2(*l2,obj)<0.25 && triggerFlavour(obj,id1) && SinglePaths.size() > 0) {
            // if (id0*id1==-165) cout<<SinglePaths.size()<<endl;
            for (size_t itrg=0; itrg<SinglePaths.size();itrg++) {
                // if (id0*id1==-165) cout<<SingleFilters[itrg]<<endl;
                if (obj.hasFilterLabel(SingleFilters[itrg])) HLTMatch_singleLeg2[itrg]=true;
            }
        }
    }
    // cout<<"step 2"<<endl;
    for (size_t itrg=0; itrg<DiPaths.size();itrg++) {
        HLTMatch_cross.push_back(false);
    }
    for (size_t idxto = 0; idxto < triggerObjects->size(); ++idxto) {//cross triggers
        pat::TriggerObjectStandAlone obj1 = triggerObjects->at(idxto);
        obj1.unpackFilterLabels(iEvent,*triggerResults );
        // obj1.unpackPathNames(names);
        for (size_t jdxto = 0; jdxto < triggerObjects->size(); ++jdxto) {
            pat::TriggerObjectStandAlone obj2 = triggerObjects->at(jdxto);
            obj2.unpackFilterLabels(iEvent,*triggerResults );
            // obj2.unpackPathNames(names);
            if (deltaR(obj1,obj2)<0.02) continue;
            if ( (deltaR2(*l1,obj1)<0.25 && triggerFlavour(obj1,id0) && deltaR2(*l2,obj2)<0.25 && triggerFlavour(obj2,id1) ) || (deltaR2(*l1,obj2)<0.25 && triggerFlavour(obj2,id0) && deltaR2(*l2,obj1)<0.25 && triggerFlavour(obj1,id1) )) {
                for (size_t itrg=0;itrg<DiPaths.size();itrg++) {
                    if (obj1.hasFilterLabel(DiFilters[itrg]) && obj2.hasFilterLabel(DiFilters[itrg])) {
                        if ( ( obj1.hasFilterLabel(DiFiltersLeg1[itrg]) && obj2.hasFilterLabel(DiFiltersLeg2[itrg]) ) || ( obj1.hasFilterLabel(DiFiltersLeg2[itrg]) && obj2.hasFilterLabel(DiFiltersLeg1[itrg]) ) ) {
                            float dr=deltaR(obj1,obj2);
                            float dz=abs(obj1.vz()-obj2.vz());
                            if ((DZ[itrg]<0 || dz<=DZ[itrg]) && (DR[itrg]<0 || dr>=DR[itrg])) {
                                if ( (DiPaths[itrg].find("Mass3p8")==string::npos || (obj1.p4()+obj2.p4()).M() > 3.8) && (DiPaths[itrg].find("Mass8")==string::npos || (obj1.p4()+obj2.p4()).M() > 8) ) {
                                    HLTMatch_cross[itrg]=true;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
      
    // cout<<"step 3"<<endl;
    bool pass_SingleTrigger1=false, pass_SingleTrigger2=false, pass_SingleTrigger=false, pass_CrossTrigger=false, pass_Trigger=false;
    for (size_t itrg=0;itrg<SinglePaths.size();itrg++) {
        myCand.addUserInt("d0.pass_"+SinglePaths[itrg],HLTMatch_singleLeg1[itrg]);
        if (HLTMatch_singleLeg1[itrg]) pass_SingleTrigger1=true;
        myCand.addUserInt("d1.pass_"+SinglePaths[itrg],HLTMatch_singleLeg2[itrg]);
        if (HLTMatch_singleLeg2[itrg]) pass_SingleTrigger2=true;
    }
    for (size_t itrg=0;itrg<DiPaths.size();itrg++) {
        myCand.addUserInt("pass_"+DiPaths[itrg],HLTMatch_cross[itrg]);
        if (HLTMatch_cross[itrg]) pass_CrossTrigger=true;
    }
    // cout<<"step 4"<<endl;
    pass_SingleTrigger=pass_SingleTrigger1 || pass_SingleTrigger2;
    pass_Trigger=pass_SingleTrigger || pass_CrossTrigger;
    myCand.addUserInt("d0.pass_SingleTrigger",pass_SingleTrigger1);
    myCand.addUserInt("d1.pass_SingleTrigger",pass_SingleTrigger2);
    myCand.addUserInt("pass_SingleTrigger",pass_SingleTrigger);
    myCand.addUserInt("pass_CrossTrigger",pass_CrossTrigger);
    myCand.addUserInt("pass_Trigger",pass_Trigger);
    
    //--------------------------------------------------------------
    //-----------------------END HLTMatch---------------------------
    //--------------------------------------------------------------
      
    //veto extra electrons or muons
    int nExtraLep = 0;
    for (vector<reco::CandidatePtr>::const_iterator lepPtr = goodisoleptonPtrs.begin(); lepPtr != goodisoleptonPtrs.end(); ++lepPtr){
      const reco::Candidate* lep = lepPtr->get();
      if (
        deltaR(lep->p4(), l1->p4()) > 0.02 &&
        deltaR(lep->p4(), l2->p4()) > 0.02
        ){
        nExtraLep++;
      }
    }
    if (nExtraLep > 0) continue;
      

    //--- Find "best Z" (closest to mZ) among those passing the "bestZAmong" selection
    if (preBestZSelection(myCand)) {
      if (abs(id0*id1)==143 || abs(id0*id1)==169) {
          bestZidx=i;
      }
      else if (abs(id0*id1)==165 || abs(id0*id1)==195) {
          if (abs(id0)==15) {
              if (userdatahelpers::getUserInt(l1,"tauVSjet") > tau1Iso || (userdatahelpers::getUserInt(l1,"tauVSjet") == tau1Iso && l1->pt() > tau1Pt)) {
                  bestZidx=i;
                  tau1Iso=userdatahelpers::getUserInt(l1,"tauVSjet");
                  tau1Pt=l1->pt();
              }
          }
          else {
              if (userdatahelpers::getUserInt(l2,"tauVSjet") > tau1Iso || (userdatahelpers::getUserInt(l2,"tauVSjet") == tau1Iso && l2->pt() > tau1Pt)) {
                  bestZidx=i;
                  tau1Iso=userdatahelpers::getUserInt(l2,"tauVSjet");
                  tau1Pt=l2->pt();
              }
          }
      }
      else if (abs(id0*id1)==225) {
          if (userdatahelpers::getUserInt(l1,"tauVSjet") > userdatahelpers::getUserInt(l2,"tauVSjet") || (userdatahelpers::getUserInt(l1,"tauVSjet") == userdatahelpers::getUserInt(l2,"tauVSjet") && l1->pt() > l2->pt() )) {
              if ((userdatahelpers::getUserInt(l1,"tauVSjet") > tau1Iso) || (userdatahelpers::getUserInt(l1,"tauVSjet") == tau1Iso && l1->pt() > tau1Pt) || (userdatahelpers::getUserInt(l1,"tauVSjet") == tau1Iso && l1->pt() == tau1Pt && userdatahelpers::getUserInt(l2,"tauVSjet") > tau2Iso) || (userdatahelpers::getUserInt(l1,"tauVSjet") == tau1Iso && l1->pt() == tau1Pt && userdatahelpers::getUserInt(l2,"tauVSjet") == tau2Iso && l2->pt() > tau2Pt)) {
                  bestZidx=i;
                  tau1Iso=userdatahelpers::getUserInt(l1,"tauVSjet");
                  tau1Pt=l1->pt();
                  tau2Iso=userdatahelpers::getUserInt(l2,"tauVSjet");
                  tau2Pt=l2->pt();
              }
          }
          else {
              if ((userdatahelpers::getUserInt(l2,"tauVSjet") > tau1Iso) || (userdatahelpers::getUserInt(l2,"tauVSjet") == tau1Iso && l2->pt() > tau1Pt) || (userdatahelpers::getUserInt(l2,"tauVSjet") == tau1Iso && l2->pt() == tau1Pt && userdatahelpers::getUserInt(l1,"tauVSjet") > tau2Iso) || (userdatahelpers::getUserInt(l2,"tauVSjet") == tau1Iso && l2->pt() == tau1Pt && userdatahelpers::getUserInt(l1,"tauVSjet") == tau2Iso && l1->pt() > tau2Pt)) {
                  bestZidx=i;
                  tau1Iso=userdatahelpers::getUserInt(l2,"tauVSjet");
                  tau1Pt=l2->pt();
                  tau2Iso=userdatahelpers::getUserInt(l1,"tauVSjet");
                  tau2Pt=l1->pt();
              }
          }
      }
    }
    result->push_back(myCand);
  }

  //--- Embed isBestZ flag (must be done in a separate loop)
  for (int i = 0; i< (int)result->size(); ++i) {
      pat::CompositeCandidate& myCand = (*result)[i];
      myCand.addUserFloat("isBestCand",  (i==bestZidx));

    //--- Embed flags (ie cuts specified in the "flags" pset)
    //    We do this here so that isBestZ is available within the cuts
    for(CutSet<pat::CompositeCandidate>::const_iterator cut = cuts.begin(); cut != cuts.end(); ++cut) {
      myCand.addUserFloat(cut->first,int((*(cut->second))(myCand)));
    }    

  }
//   cout<<result->size()<<" Z candidates"<<endl;
  iEvent.put(std::move(result));
}

bool ZCandidateFiller::triggerFlavour(pat::TriggerObjectStandAlone OBJ, int flavour)
{
    if (abs(flavour)==11)
        return ( OBJ.hasTriggerObjectType(trigger::TriggerElectron) || OBJ.hasTriggerObjectType(trigger::TriggerPhoton) );
    else if (abs(flavour)==13)
        return OBJ.hasTriggerObjectType(trigger::TriggerMuon);
    else if (abs(flavour)==15) {
        // cout<<"Tau trigger: "<<(OBJ.hasTriggerObjectType(trigger::TriggerTau) || OBJ.hasTriggerObjectType(trigger::TriggerL1TauJet))<<endl;
        return ( OBJ.hasTriggerObjectType(trigger::TriggerTau) || OBJ.hasTriggerObjectType(trigger::TriggerL1TauJet) || OBJ.hasTriggerObjectType(trigger::TriggerL1Tau));
    }
    else return false;
}

#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(ZCandidateFiller);
