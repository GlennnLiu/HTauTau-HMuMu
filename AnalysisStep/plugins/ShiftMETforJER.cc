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

class ShiftMETforJER : public edm::EDProducer {
    public: 
        /// Constructor
        explicit ShiftMETforJER(const edm::ParameterSet&);
        /// Destructor
        ShiftMETforJER(){};

    private:
        virtual void beginJob(){};  
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob(){};

        edm::EDGetTokenT<View<pat::MET>> theMETTag;
        edm::EDGetTokenT<View<pat::Jet>> theJetTag;
};

ShiftMETforJER::ShiftMETforJER(const edm::ParameterSet& iConfig) :
theMETTag(consumes<View<pat::MET>>(iConfig.getParameter<edm::InputTag>("srcMET"))),
theJetTag(consumes<View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jetCollection")))
{
    produces<double>("METdxUPJER");
    produces<double>("METdyUPJER");
    produces<double>("METdxDOWNJER");
    produces<double>("METdyDOWNJER");
    
}

void ShiftMETforJER::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
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
    int itau=0;
    for(auto jet = jetHandle->begin(); jet != jetHandle->end(); ++jet)
    {
      //cout << "-----> iTau " << itau << endl;

      //---Clone the pat::Tau
      pat::Jet j(*jet);
        
      // Unshifted tau
      TLorentzVector pfour;
      pfour.SetPxPyPzE (j.px(), j.py(), j.pz(), j.energy());
      
      // Shifted taus
      TLorentzVector pfourJERUp;
      TLorentzVector pfourJERDown;

      pfourJERUp=pfour*(j.userFloat("pt_jerup")/j.pt());
      pfourJERDown=pfour*(j.userFloat("pt_jerdn")/j.pt());
      
      deltaJets_UP += ( pfourJERUp - pfour );
      deltaJets_DOWN += ( pfourJERDown - pfour );
    
    } // end loop on taus
    
    // Calculate the correction
    (*dx_UP_ptr) = patMET.px() - deltaJets_UP.Px();
    (*dy_UP_ptr) = patMET.py() - deltaJets_UP.Py();

    (*dx_DOWN_ptr) = patMET.px() - deltaJets_DOWN.Px();
    (*dy_DOWN_ptr) = patMET.py() - deltaJets_DOWN.Py();
    
    //cout << "SHIFTED UP  : " << *dx_UP_ptr << " / " << *dy_UP_ptr << endl;
    //cout << "SHIFTED DOWN: " << *dx_DOWN_ptr << " / " << *dy_DOWN_ptr << endl;

    iEvent.put( std::move(dx_UP_ptr)  , "METdxUPJER"   );
    iEvent.put( std::move(dy_UP_ptr)  , "METdyUPJER"   );
    iEvent.put( std::move(dx_DOWN_ptr), "METdxDOWNJER" );
    iEvent.put( std::move(dy_DOWN_ptr), "METdyDOWNJER" );
    
}

#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(ShiftMETforJER);
