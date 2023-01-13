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
};

ShiftMETforJES::ShiftMETforJES(const edm::ParameterSet& iConfig) :
theMETTag(consumes<View<pat::MET>>(iConfig.getParameter<edm::InputTag>("srcMET"))),
theJetTag(consumes<View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jetCollection")))
{
    produces<double>("METdxUPJES");
    produces<double>("METdyUPJES");
    produces<double>("METdxDOWNJES");
    produces<double>("METdyDOWNJES");
    
}

void ShiftMETforJES::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    // Declare ptrs to save MET variations
    std::unique_ptr<double> dx_UP_ptr   (new double);
    std::unique_ptr<double> dy_UP_ptr   (new double);
    std::unique_ptr<double> dx_DOWN_ptr (new double);
    std::unique_ptr<double> dy_DOWN_ptr (new double);

    // Get the MET
    Handle<View<pat::MET> > METHandle;
    iEvent.getByToken(theMETTag, METHandle);
    const pat::MET& patMET = (*METHandle)[0];
//     cout<<patMET.px()<<patMET.py()<<endl;
   
    // Get the Jets
    Handle<View<pat::Jet>> jetHandle;
    iEvent.getByToken(theJetTag, jetHandle);
    
    // Define the correction of the met
    TLorentzVector deltaJets_UP  ;
    TLorentzVector deltaJets_DOWN;
    
    //cout << "--------- *** SHIFTED BEGIN LOOP" << endl;
    
    // Loop on taus
    //for (unsigned int itau = 0; itau < tauHandle->size(); ++itau)
    for(auto jet = jetHandle->begin(); jet != jetHandle->end(); ++jet)
    {
      //cout << "-----> iTau " << itau << endl;

      //---Clone the pat::Tau
      pat::Jet j(*jet);
        
      // Unshifted tau
      TLorentzVector pfour;
      pfour.SetPxPyPzE (j.px(), j.py(), j.pz(), j.energy());
      
      // Shifted taus
      TLorentzVector pfourJESUp;
      TLorentzVector pfourJESDown;

      pfourJESUp=pfour*(j.userFloat("pt_jesup")/j.pt());
      pfourJESDown=pfour*(j.userFloat("pt_jesdn")/j.pt());
      
      deltaJets_UP += ( pfourJESUp - pfour );
      deltaJets_DOWN += ( pfourJESDown - pfour );
    
    } // end loop on taus
    
    // Calculate the correction
    (*dx_UP_ptr) = patMET.px() - deltaJets_UP.Px();
    (*dy_UP_ptr) = patMET.py() - deltaJets_UP.Py();

    (*dx_DOWN_ptr) = patMET.px() - deltaJets_DOWN.Px();
    (*dy_DOWN_ptr) = patMET.py() - deltaJets_DOWN.Py();
    
    //cout << "SHIFTED UP  : " << *dx_UP_ptr << " / " << *dy_UP_ptr << endl;
    //cout << "SHIFTED DOWN: " << *dx_DOWN_ptr << " / " << *dy_DOWN_ptr << endl;

    iEvent.put( std::move(dx_UP_ptr)  , "METdxUPJES"   );
    iEvent.put( std::move(dy_UP_ptr)  , "METdyUPJES"   );
    iEvent.put( std::move(dx_DOWN_ptr), "METdxDOWNJES" );
    iEvent.put( std::move(dy_DOWN_ptr), "METdyDOWNJES" );
    
}

#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(ShiftMETforJES);
