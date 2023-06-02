#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Utilities/interface/InputTag.h>
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

class ShiftMETforJES : public edm::EDProducer {
    public: 
        /// Constructor
        explicit ShiftMETforJES(const edm::ParameterSet&);
        /// Destructor
        ShiftMETforJES(){};

    private:
        virtual void beginJob(){};  
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob(){};

        edm::EDGetTokenT<View<pat::MET>> theMETTag;
        edm::EDGetTokenT<View<pat::Jet>> theJetTag;
        std::vector<string> uncSources {};
};

ShiftMETforJES::ShiftMETforJES(const edm::ParameterSet& iConfig) :
theMETTag(consumes<View<pat::MET>>(iConfig.getParameter<edm::InputTag>("srcMET"))),
theJetTag(consumes<View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jetCollection")))
{
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
    
    produces<pat::METCollection>("up");
    produces<pat::METCollection>("down");
}

void ShiftMETforJES::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    // Declare ptrs to save MET variations
    std::unique_ptr<pat::METCollection> out_MET_up(new pat::METCollection());
    std::unique_ptr<pat::METCollection> out_MET_down(new pat::METCollection());

    // Get the MET
    Handle<View<pat::MET> > METHandle;
    iEvent.getByToken(theMETTag, METHandle);
    const pat::MET& patMET = (*METHandle)[0];
//     cout<<patMET.px()<<patMET.py()<<endl;
   
    // Get the Jets
    Handle<View<pat::Jet>> jetHandle;
    iEvent.getByToken(theJetTag, jetHandle);
    
    for (unsigned s_unc = 0; s_unc < uncSources.size(); s_unc++) {
        // Define the correction of the met
        TLorentzVector deltaJets_UP  ;
        TLorentzVector deltaJets_DOWN;
        
        //cout << "--------- *** SHIFTED BEGIN LOOP" << endl;
        
        // Loop on taus
        //for (unsigned int itau = 0; itau < tauHandle->size(); ++itau)
        for (auto jet = jetHandle->begin(); jet != jetHandle->end(); ++jet) {
            //cout << "-----> iTau " << itau << endl;

            //---Clone the pat::Tau
            pat::Jet j(*jet);
                
            // Unshifted tau
            TLorentzVector pfour;
            pfour.SetPxPyPzE (j.px(), j.py(), j.pz(), j.energy());
            
            // Shifted taus
            TLorentzVector pfourJESUp;
            TLorentzVector pfourJESDown;

            pfourJESUp=pfour*(j.userFloat("pt_jesup_split_"+uncSources[s_unc])/j.pt());
            pfourJESDown=pfour*(j.userFloat("pt_jesdn_split_"+uncSources[s_unc])/j.pt());
            
            deltaJets_UP += ( pfourJESUp - pfour );
            deltaJets_DOWN += ( pfourJESDown - pfour );
            
        } // end loop on taus
        float shiftMetPxUp = patMET.px() - deltaJets_UP.Px();
        float shiftMetPyUp = patMET.py() - deltaJets_UP.Py();
        reco::Candidate::LorentzVector shiftedMetP4Up(shiftMetPxUp, shiftMetPyUp, 0., sqrt(shiftMetPxUp*shiftMetPxUp + shiftMetPyUp*shiftMetPyUp));
        pat::MET corrMEtUp(patMET);
        corrMEtUp.setP4(shiftedMetP4Up);
        corrMEtUp.setSignificanceMatrix(patMET.getSignificanceMatrix());
        out_MET_up->push_back(corrMEtUp);

        float shiftMetPxDn = patMET.px() - deltaJets_DOWN.Px();
        float shiftMetPyDn = patMET.py() - deltaJets_DOWN.Py();
        reco::Candidate::LorentzVector shiftedMetP4Dn(shiftMetPxDn, shiftMetPyDn, 0., sqrt(shiftMetPxDn*shiftMetPxDn + shiftMetPyDn*shiftMetPyDn));
        pat::MET corrMEtDn(patMET);
        corrMEtDn.setP4(shiftedMetP4Dn);
        corrMEtDn.setSignificanceMatrix(patMET.getSignificanceMatrix());
        out_MET_down->push_back(corrMEtDn);
    }
    
    //cout << "SHIFTED UP  : " << *dx_UP_ptr << " / " << *dy_UP_ptr << endl;
    //cout << "SHIFTED DOWN: " << *dx_DOWN_ptr << " / " << *dy_DOWN_ptr << endl;

    iEvent.put( std::move(out_MET_up)  , "up"   );
    iEvent.put( std::move(out_MET_down)  , "down"   );
    
}

#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(ShiftMETforJES);
