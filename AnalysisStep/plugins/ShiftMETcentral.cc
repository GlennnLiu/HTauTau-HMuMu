/*
** class  : ShiftMETcentral
** author : F. Brivio (MIB)
** date   : 12 March 2020
** brief  : takes in input the met, the uncorrected taus (bareTaus) and the corrected taus (softTaus) 
**          and produces a pat::METCollection shifted due to TES and EES central corrections of the taus
*/

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Utilities/interface/InputTag.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/METReco/interface/PFMET.h>
#include <DataFormats/METReco/interface/PFMETCollection.h>
#include <DataFormats/METReco/interface/CommonMETData.h>
#include <DataFormats/PatCandidates/interface/Tau.h>
#include <DataFormats/PatCandidates/interface/Jet.h>
#include "TLorentzVector.h"
#include <HTauTauHMuMu/AnalysisStep/interface/DaughterDataHelpers.h>
#include <DataFormats/Candidate/interface/Candidate.h>

#include <vector>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;

using LorentzVectorE = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double>>;

class ShiftMETcentral : public edm::EDProducer {
    public: 
        /// Constructor
        explicit ShiftMETcentral(const edm::ParameterSet&);
        /// Destructor
        ~ShiftMETcentral(){};

    private:
        virtual void beginJob(){};  
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob(){};

        edm::EDGetTokenT<View<pat::MET>> theMETTag;
        edm::EDGetTokenT<pat::TauRefVector> theTauUncorrectedTag;
        edm::EDGetTokenT<pat::TauCollection> theTauCorrectedTag;
        edm::EDGetTokenT<View<pat::Jet>> theJetTag;
};

ShiftMETcentral::ShiftMETcentral(const edm::ParameterSet& iConfig) :
theMETTag(consumes<View<pat::MET>>(iConfig.getParameter<edm::InputTag>("srcMET"))),
theTauUncorrectedTag(consumes<pat::TauRefVector>(iConfig.getParameter<edm::InputTag>("tauUncorrected"))),
theTauCorrectedTag(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("tauCorrected"))),
theJetTag(consumes<View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jetCollection")))
{
    produces<pat::METCollection>();
}

void ShiftMETcentral::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    // Declare ptrs to save MET variations
    std::unique_ptr<pat::METCollection> out_MET_ptr(new pat::METCollection());

    // Get the MET
    Handle<View<pat::MET> > METHandle;
    iEvent.getByToken(theMETTag, METHandle);
    const pat::MET& patMET = (*METHandle)[0];
   
    // Get the Uncorrected Taus
    Handle<pat::TauRefVector> tauUncorrectedHandle;
    iEvent.getByToken(theTauUncorrectedTag, tauUncorrectedHandle);

    // Get the Corrected Taus
    Handle<pat::TauCollection> tauCorrectedHandle;
    iEvent.getByToken(theTauCorrectedTag, tauCorrectedHandle);
    
    // Get Jets
    Handle<View<pat::Jet>> jetHandle;
    iEvent.getByToken(theJetTag, jetHandle);

    // Define the correction of the met
    TLorentzVector deltaTaus, deltaJets;
    //cout << "-- Delta : " << deltaTaus.Pt() << " / " << deltaTaus.Eta() << endl;
    //cout << "-- MET   : " << patMET.px() << " / " << patMET.py() << endl;
    //cout << "tauUncorrected size: " << tauUncorrectedHandle->size() << endl;
    //cout << "tauCorrected size  : " << tauCorrectedHandle->size() << endl;
    
    // Loop on taus
    //cout << "--------- *** BEGIN LOOP" << endl;
    int itau=0;
    for(pat::TauCollection::const_iterator inputTau = tauCorrectedHandle->begin(); inputTau != tauCorrectedHandle->end(); ++inputTau, ++itau)
    {
      pat::Tau tauCorrected(*inputTau);
      pat::Tau tauUncorrected(*((*tauUncorrectedHandle)[itau].get()));

      // cast to reco::Candidate*
      const reco::Candidate* l_Uncorrected = (const reco::Candidate*)&tauUncorrected;
      const reco::Candidate* l_Corrected   = (const reco::Candidate*)&tauCorrected;
      // Unshifted tau
      TLorentzVector pfour;
      pfour.SetPxPyPzE(l_Uncorrected->px(), l_Uncorrected->py(), l_Uncorrected->pz(), l_Uncorrected->energy());
      
      // Shifted tau
      TLorentzVector pfourShifted;
      pfourShifted.SetPxPyPzE(l_Corrected->px(), l_Corrected->py(), l_Corrected->pz(), l_Corrected->energy());

      // Compute delta
      deltaTaus += ( pfourShifted - pfour );
    } // end loop on taus

    for(auto jet = jetHandle->begin(); jet != jetHandle->end(); ++jet) {

      pat::Jet j(*jet);
      
      TLorentzVector pfour, pfourShifted;
      pfourShifted.SetPxPyPzE(j.px(),j.py(),j.pz(),j.energy());
      
      if (j.pt()>0) deltaJets += pfourShifted*(1-j.userFloat("pt_JEC_noJER")/j.pt());
    }
    
    // //MET recoil correction
    // TString fileName;
    // if (setup==2016) fileName="RecoilCorrections/data/TypeI-PFMet_Run2016_legacy.root";
    // else if (setup==2017) fileName="RecoilCorrections/data/Type1_PFMET_2017.root";
    // else fileName="RecoilCorrections/data/TypeI-PFMet_Run2018.root";

    // raw one
    out_MET_ptr->push_back(patMET);
    
    // Calculate the correction
    float shiftMetPx = patMET.px() - deltaTaus.Px() - deltaJets.Px();
    float shiftMetPy = patMET.py() - deltaTaus.Py() - deltaJets.Py();
    reco::Candidate::LorentzVector shiftedMetP4(shiftMetPx, shiftMetPy, 0., sqrt(shiftMetPx*shiftMetPx + shiftMetPy*shiftMetPy));
    pat::MET corrMEt(patMET);
    corrMEt.setP4(shiftedMetP4);
    corrMEt.setSignificanceMatrix(patMET.getSignificanceMatrix());
    out_MET_ptr->push_back(corrMEt);

    // only correct tau
    float shiftMetPxTau = patMET.px() - deltaTaus.Px();
    float shiftMetPyTau = patMET.py() - deltaTaus.Py();
    reco::Candidate::LorentzVector shiftedMetP4Tau(shiftMetPxTau, shiftMetPyTau, 0., sqrt(shiftMetPxTau*shiftMetPxTau + shiftMetPyTau*shiftMetPyTau));
    pat::MET corrMEtTau(patMET);
    corrMEtTau.setP4(shiftedMetP4Tau);
    corrMEtTau.setSignificanceMatrix(patMET.getSignificanceMatrix());
    out_MET_ptr->push_back(corrMEtTau);

    // only correct jet
    float shiftMetPxJet = patMET.px() - deltaJets.Px();
    float shiftMetPyJet = patMET.py() - deltaJets.Py();
    reco::Candidate::LorentzVector shiftedMetP4Jet(shiftMetPxJet, shiftMetPyJet, 0., sqrt(shiftMetPxJet*shiftMetPxJet + shiftMetPyJet*shiftMetPyJet));
    pat::MET corrMEtJet(patMET);
    corrMEtJet.setP4(shiftedMetP4Jet);
    corrMEtJet.setSignificanceMatrix(patMET.getSignificanceMatrix());
    out_MET_ptr->push_back(corrMEtJet);

    iEvent.put(std::move(out_MET_ptr));
}

#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(ShiftMETcentral);
