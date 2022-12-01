// Include classes
#include <HTauTauHMuMu/AnalysisStep/test/ZpXEstimation/include/FakeRates.h>

using namespace std;

// Constructor
//===============================================
FakeRates::FakeRates( TString input_file_FR_name )
{
 
   input_file_FR = new TFile(input_file_FR_name);

   // for OS FR files
   if (input_file_FR_name.Contains("OS"))
   {
      g_FR_mu_EB = (TH1F*)input_file_FR->Get("FR_OS_muon_EB");
      g_FR_mu_EE = (TH1F*)input_file_FR->Get("FR_OS_muon_EE");
      g_FR_e_EB  = (TH1F*)input_file_FR->Get("FR_OS_electron_EB");
      g_FR_e_EE  = (TH1F*)input_file_FR->Get("FR_OS_electron_EE");
      g_FR_tauE_Decay0 = (TH2F*)input_file_FR->Get("FR_OS_tauE_Decay0");
      g_FR_tauE_Decay1 = (TH2F*)input_file_FR->Get("FR_OS_tauE_Decay1");
      g_FR_tauE_Decay10 = (TH2F*)input_file_FR->Get("FR_OS_tauE_Decay10");
      g_FR_tauE_Decay11 = (TH2F*)input_file_FR->Get("FR_OS_tauE_Decay11");
      g_FR_tauMu_Decay0 = (TH2F*)input_file_FR->Get("FR_OS_tauMu_Decay0");
      g_FR_tauMu_Decay1 = (TH2F*)input_file_FR->Get("FR_OS_tauMu_Decay1");
      g_FR_tauMu_Decay10 = (TH2F*)input_file_FR->Get("FR_OS_tauMu_Decay10");
      g_FR_tauMu_Decay11 = (TH2F*)input_file_FR->Get("FR_OS_tauMu_Decay11");
      g_FR_tauTau_Decay0 = (TH2F*)input_file_FR->Get("FR_OS_tauTau_Decay0");
      g_FR_tauTau_Decay1 = (TH2F*)input_file_FR->Get("FR_OS_tauTau_Decay1");
      g_FR_tauTau_Decay10 = (TH2F*)input_file_FR->Get("FR_OS_tauTau_Decay10");
      g_FR_tauTau_Decay11 = (TH2F*)input_file_FR->Get("FR_OS_tauTau_Decay11");
   }
   
   // for SS FR files
   if (input_file_FR_name.Contains("SS"))
   {
      g_FR_mu_EB = (TH1F*)input_file_FR->Get("FR_SS_muon_EB");
      g_FR_mu_EE = (TH1F*)input_file_FR->Get("FR_SS_muon_EE");
      g_FR_e_EB  = (TH1F*)input_file_FR->Get("FR_SS_electron_EB");
      g_FR_e_EE  = (TH1F*)input_file_FR->Get("FR_SS_electron_EE");
      g_FR_tauE_Decay0 = (TH2F*)input_file_FR->Get("FR_SS_tauE_Decay0");
      g_FR_tauE_Decay1 = (TH2F*)input_file_FR->Get("FR_SS_tauE_Decay1");
      g_FR_tauE_Decay10 = (TH2F*)input_file_FR->Get("FR_SS_tauE_Decay10");
      g_FR_tauE_Decay11 = (TH2F*)input_file_FR->Get("FR_SS_tauE_Decay11");
      g_FR_tauMu_Decay0 = (TH2F*)input_file_FR->Get("FR_SS_tauMu_Decay0");
      g_FR_tauMu_Decay1 = (TH2F*)input_file_FR->Get("FR_SS_tauMu_Decay1");
      g_FR_tauMu_Decay10 = (TH2F*)input_file_FR->Get("FR_SS_tauMu_Decay10");
      g_FR_tauMu_Decay11 = (TH2F*)input_file_FR->Get("FR_SS_tauMu_Decay11");
      g_FR_tauTau_Decay0 = (TH2F*)input_file_FR->Get("FR_SS_tauTau_Decay0");
      g_FR_tauTau_Decay1 = (TH2F*)input_file_FR->Get("FR_SS_tauTau_Decay1");
      g_FR_tauTau_Decay10 = (TH2F*)input_file_FR->Get("FR_SS_tauTau_Decay10");
      g_FR_tauTau_Decay11 = (TH2F*)input_file_FR->Get("FR_SS_tauTau_Decay11");
   }

}
//===============================================



//======================
FakeRates::~FakeRates() {}
//======================



//==================================================================
float FakeRates::GetFakeRate(float lep_Pt, float MET, float lep_eta, int decayMode, int lep_ID, int tauChannel)
{
    float my_lep_Pt = lep_Pt >= 80. ? 79. : lep_Pt;
    float my_MET = MET >= 100. ? 99. : MET;
    int   my_lep_ID = abs(lep_ID);

    if ( my_lep_ID == 11 ) {
        if ( fabs(lep_eta) < 1.479 ) return g_FR_e_EB->GetBinContent(g_FR_e_EB->FindBin(my_lep_Pt));
        else return g_FR_e_EE->GetBinContent(g_FR_e_EE->FindBin(my_lep_Pt));
    }
    else if ( my_lep_ID == 13 ) {
        if ( fabs(lep_eta) < 1.2 ) return g_FR_e_EB->GetBinContent(g_FR_e_EB->FindBin(my_lep_Pt));
        else return g_FR_e_EE->GetBinContent(g_FR_e_EE->FindBin(my_lep_Pt));
    }
    else if ( my_lep_ID == 15 ) {
        if ( tauChannel == 0 ) {
            if (decayMode==0) return g_FR_tauE_Decay0->GetBinContent(g_FR_tauE_Decay0->FindBin(my_lep_Pt,my_MET));
            else if (decayMode==1) return g_FR_tauE_Decay1->GetBinContent(g_FR_tauE_Decay1->FindBin(my_lep_Pt,my_MET));
            else if (decayMode==10) return g_FR_tauE_Decay10->GetBinContent(g_FR_tauE_Decay10->FindBin(my_lep_Pt,my_MET));
            else if (decayMode==11) return g_FR_tauE_Decay11->GetBinContent(g_FR_tauE_Decay11->FindBin(my_lep_Pt,my_MET));
            else { cout<<"[ERROR] Decay mode "<<decayMode<<" not correct!"<<endl; return 0;}
        }
        else if ( tauChannel == 1 ) {
            if (decayMode==0) return g_FR_tauMu_Decay0->GetBinContent(g_FR_tauMu_Decay0->FindBin(my_lep_Pt,my_MET));
            else if (decayMode==1) return g_FR_tauMu_Decay1->GetBinContent(g_FR_tauMu_Decay1->FindBin(my_lep_Pt,my_MET));
            else if (decayMode==10) return g_FR_tauMu_Decay10->GetBinContent(g_FR_tauMu_Decay10->FindBin(my_lep_Pt,my_MET));
            else if (decayMode==11) return g_FR_tauMu_Decay11->GetBinContent(g_FR_tauMu_Decay11->FindBin(my_lep_Pt,my_MET));
            else { cout<<"[ERROR] Decay mode "<<decayMode<<" not correct!"<<endl; return 0;}
        }
        else if ( tauChannel == 2 ) {
            if (decayMode==0) return g_FR_tauTau_Decay0->GetBinContent(g_FR_tauTau_Decay0->FindBin(my_lep_Pt,my_MET));
            else if (decayMode==1) return g_FR_tauTau_Decay1->GetBinContent(g_FR_tauTau_Decay1->FindBin(my_lep_Pt,my_MET));
            else if (decayMode==10) return g_FR_tauTau_Decay10->GetBinContent(g_FR_tauTau_Decay10->FindBin(my_lep_Pt,my_MET));
            else if (decayMode==11) return g_FR_tauTau_Decay11->GetBinContent(g_FR_tauTau_Decay11->FindBin(my_lep_Pt,my_MET));
            else { cout<<"[ERROR] Decay mode "<<decayMode<<" not correct!"<<endl; return 0;}
        }
        else {
            cout << "[ERROR] Wrong tau channel: " << tauChannel << endl;
            return 0;
        }
    }
    else {
        cout << "[ERROR] Wrong lepton ID: " << my_lep_ID << endl;
        return 0;
    }
}
//==================================================================

//==================================================================
float FakeRates::GetFakeRate_Dn(float lep_Pt, float MET, float lep_eta, int decayMode, int lep_ID, int tauChannel)
{
    float my_lep_Pt = lep_Pt >= 80. ? 79. : lep_Pt;
    float my_MET = MET >= 100. ? 99. : MET;
    int   my_lep_ID = abs(lep_ID);

    if ( my_lep_ID == 11 ) {
        if ( fabs(lep_eta) < 1.479 ) return g_FR_e_EB->GetBinContent(g_FR_e_EB->FindBin(my_lep_Pt)) - g_FR_e_EB->GetBinError(g_FR_e_EB->FindBin(my_lep_Pt));
        else return g_FR_e_EE->GetBinContent(g_FR_e_EE->FindBin(my_lep_Pt)) - g_FR_e_EE->GetBinError(g_FR_e_EE->FindBin(my_lep_Pt));
    }
    else if ( my_lep_ID == 13 ) {
        if ( fabs(lep_eta) < 1.2 ) return g_FR_e_EB->GetBinContent(g_FR_e_EB->FindBin(my_lep_Pt)) - g_FR_e_EB->GetBinError(g_FR_e_EB->FindBin(my_lep_Pt));
        else return g_FR_e_EE->GetBinContent(g_FR_e_EE->FindBin(my_lep_Pt)) - g_FR_e_EE->GetBinError(g_FR_e_EE->FindBin(my_lep_Pt));
    }
    else if ( my_lep_ID == 15 ) {
        if ( tauChannel == 0 ) {
            if (decayMode==0) return g_FR_tauE_Decay0->GetBinContent(g_FR_tauE_Decay0->FindBin(my_lep_Pt,my_MET)) - g_FR_tauE_Decay0->GetBinError(g_FR_tauE_Decay0->FindBin(my_lep_Pt,my_MET));
            else if (decayMode==1) return g_FR_tauE_Decay1->GetBinContent(g_FR_tauE_Decay1->FindBin(my_lep_Pt,my_MET)) - g_FR_tauE_Decay1->GetBinError(g_FR_tauE_Decay1->FindBin(my_lep_Pt,my_MET));
            else if (decayMode==10) return g_FR_tauE_Decay10->GetBinContent(g_FR_tauE_Decay10->FindBin(my_lep_Pt,my_MET)) - g_FR_tauE_Decay10->GetBinError(g_FR_tauE_Decay10->FindBin(my_lep_Pt,my_MET));
            else if (decayMode==11) return g_FR_tauE_Decay11->GetBinContent(g_FR_tauE_Decay11->FindBin(my_lep_Pt,my_MET)) - g_FR_tauE_Decay11->GetBinError(g_FR_tauE_Decay11->FindBin(my_lep_Pt,my_MET));
            else { cout<<"[ERROR] Decay mode "<<decayMode<<" not correct!"<<endl; return 0;}
        }
        else if ( tauChannel == 1 ) {
            if (decayMode==0) return g_FR_tauMu_Decay0->GetBinContent(g_FR_tauMu_Decay0->FindBin(my_lep_Pt,my_MET)) - g_FR_tauMu_Decay0->GetBinError(g_FR_tauMu_Decay0->FindBin(my_lep_Pt,my_MET));
            else if (decayMode==1) return g_FR_tauMu_Decay1->GetBinContent(g_FR_tauMu_Decay1->FindBin(my_lep_Pt,my_MET)) - g_FR_tauMu_Decay1->GetBinError(g_FR_tauMu_Decay1->FindBin(my_lep_Pt,my_MET));
            else if (decayMode==10) return g_FR_tauMu_Decay10->GetBinContent(g_FR_tauMu_Decay10->FindBin(my_lep_Pt,my_MET)) - g_FR_tauMu_Decay10->GetBinError(g_FR_tauMu_Decay10->FindBin(my_lep_Pt,my_MET));
            else if (decayMode==11) return g_FR_tauMu_Decay11->GetBinContent(g_FR_tauMu_Decay11->FindBin(my_lep_Pt,my_MET)) - g_FR_tauMu_Decay11->GetBinError(g_FR_tauMu_Decay11->FindBin(my_lep_Pt,my_MET));
            else { cout<<"[ERROR] Decay mode "<<decayMode<<" not correct!"<<endl; return 0;}
        }
        else if ( tauChannel == 2 ) {
            if (decayMode==0) return g_FR_tauTau_Decay0->GetBinContent(g_FR_tauTau_Decay0->FindBin(my_lep_Pt,my_MET)) - g_FR_tauTau_Decay0->GetBinError(g_FR_tauTau_Decay0->FindBin(my_lep_Pt,my_MET));
            else if (decayMode==1) return g_FR_tauTau_Decay1->GetBinContent(g_FR_tauTau_Decay1->FindBin(my_lep_Pt,my_MET)) - g_FR_tauTau_Decay1->GetBinError(g_FR_tauTau_Decay1->FindBin(my_lep_Pt,my_MET));
            else if (decayMode==10) return g_FR_tauTau_Decay10->GetBinContent(g_FR_tauTau_Decay10->FindBin(my_lep_Pt,my_MET)) - g_FR_tauTau_Decay10->GetBinError(g_FR_tauTau_Decay10->FindBin(my_lep_Pt,my_MET));
            else if (decayMode==11) return g_FR_tauTau_Decay11->GetBinContent(g_FR_tauTau_Decay11->FindBin(my_lep_Pt,my_MET)) - g_FR_tauTau_Decay11->GetBinError(g_FR_tauTau_Decay11->FindBin(my_lep_Pt,my_MET));
            else { cout<<"[ERROR] Decay mode "<<decayMode<<" not correct!"<<endl; return 0;}
        }
        else {
            cout << "[ERROR] Wrong tau channel: " << tauChannel << endl;
            return 0;
        }
    }
    else {
        cout << "[ERROR] Wrong lepton ID: " << my_lep_ID << endl;
        return 0;
    }
}
//==================================================================

//==================================================================
float FakeRates::GetFakeRate_Up(float lep_Pt, float MET, float lep_eta, int decayMode, int lep_ID, int tauChannel)
{
    float my_lep_Pt = lep_Pt >= 80. ? 79. : lep_Pt;
    float my_MET = MET >= 100. ? 99. : MET;
    int   my_lep_ID = abs(lep_ID);

    if ( my_lep_ID == 11 ) {
        if ( fabs(lep_eta) < 1.479 ) return g_FR_e_EB->GetBinContent(g_FR_e_EB->FindBin(my_lep_Pt)) + g_FR_e_EB->GetBinError(g_FR_e_EB->FindBin(my_lep_Pt));
        else return g_FR_e_EE->GetBinContent(g_FR_e_EE->FindBin(my_lep_Pt)) + g_FR_e_EE->GetBinError(g_FR_e_EE->FindBin(my_lep_Pt));
    }
    else if ( my_lep_ID == 13 ) {
        if ( fabs(lep_eta) < 1.2 ) return g_FR_e_EB->GetBinContent(g_FR_e_EB->FindBin(my_lep_Pt)) + g_FR_e_EB->GetBinError(g_FR_e_EB->FindBin(my_lep_Pt));
        else return g_FR_e_EE->GetBinContent(g_FR_e_EE->FindBin(my_lep_Pt)) + g_FR_e_EE->GetBinError(g_FR_e_EE->FindBin(my_lep_Pt));
    }
    else if ( my_lep_ID == 15 ) {
        if ( tauChannel == 0 ) {
            if (decayMode==0) return g_FR_tauE_Decay0->GetBinContent(g_FR_tauE_Decay0->FindBin(my_lep_Pt,my_MET)) + g_FR_tauE_Decay0->GetBinError(g_FR_tauE_Decay0->FindBin(my_lep_Pt,my_MET));
            else if (decayMode==1) return g_FR_tauE_Decay1->GetBinContent(g_FR_tauE_Decay1->FindBin(my_lep_Pt,my_MET)) + g_FR_tauE_Decay1->GetBinError(g_FR_tauE_Decay1->FindBin(my_lep_Pt,my_MET));
            else if (decayMode==10) return g_FR_tauE_Decay10->GetBinContent(g_FR_tauE_Decay10->FindBin(my_lep_Pt,my_MET)) + g_FR_tauE_Decay10->GetBinError(g_FR_tauE_Decay10->FindBin(my_lep_Pt,my_MET));
            else if (decayMode==11) return g_FR_tauE_Decay11->GetBinContent(g_FR_tauE_Decay11->FindBin(my_lep_Pt,my_MET)) + g_FR_tauE_Decay11->GetBinError(g_FR_tauE_Decay11->FindBin(my_lep_Pt,my_MET));
            else { cout<<"[ERROR] Decay mode "<<decayMode<<" not correct!"<<endl; return 0;}
        }
        else if ( tauChannel == 1 ) {
            if (decayMode==0) return g_FR_tauMu_Decay0->GetBinContent(g_FR_tauMu_Decay0->FindBin(my_lep_Pt,my_MET)) + g_FR_tauMu_Decay0->GetBinError(g_FR_tauMu_Decay0->FindBin(my_lep_Pt,my_MET));
            else if (decayMode==1) return g_FR_tauMu_Decay1->GetBinContent(g_FR_tauMu_Decay1->FindBin(my_lep_Pt,my_MET)) + g_FR_tauMu_Decay1->GetBinError(g_FR_tauMu_Decay1->FindBin(my_lep_Pt,my_MET));
            else if (decayMode==10) return g_FR_tauMu_Decay10->GetBinContent(g_FR_tauMu_Decay10->FindBin(my_lep_Pt,my_MET)) + g_FR_tauMu_Decay10->GetBinError(g_FR_tauMu_Decay10->FindBin(my_lep_Pt,my_MET));
            else if (decayMode==11) return g_FR_tauMu_Decay11->GetBinContent(g_FR_tauMu_Decay11->FindBin(my_lep_Pt,my_MET)) + g_FR_tauMu_Decay11->GetBinError(g_FR_tauMu_Decay11->FindBin(my_lep_Pt,my_MET));
            else { cout<<"[ERROR] Decay mode "<<decayMode<<" not correct!"<<endl; return 0;}
        }
        else if ( tauChannel == 2 ) {
            if (decayMode==0) return g_FR_tauTau_Decay0->GetBinContent(g_FR_tauTau_Decay0->FindBin(my_lep_Pt,my_MET)) + g_FR_tauTau_Decay0->GetBinError(g_FR_tauTau_Decay0->FindBin(my_lep_Pt,my_MET));
            else if (decayMode==1) return g_FR_tauTau_Decay1->GetBinContent(g_FR_tauTau_Decay1->FindBin(my_lep_Pt,my_MET)) + g_FR_tauTau_Decay1->GetBinError(g_FR_tauTau_Decay1->FindBin(my_lep_Pt,my_MET));
            else if (decayMode==10) return g_FR_tauTau_Decay10->GetBinContent(g_FR_tauTau_Decay10->FindBin(my_lep_Pt,my_MET)) + g_FR_tauTau_Decay10->GetBinError(g_FR_tauTau_Decay10->FindBin(my_lep_Pt,my_MET));
            else if (decayMode==11) return g_FR_tauTau_Decay11->GetBinContent(g_FR_tauTau_Decay11->FindBin(my_lep_Pt,my_MET)) + g_FR_tauTau_Decay11->GetBinError(g_FR_tauTau_Decay11->FindBin(my_lep_Pt,my_MET));
            else { cout<<"[ERROR] Decay mode "<<decayMode<<" not correct!"<<endl; return 0;}
        }
        else {
            cout << "[ERROR] Wrong tau channel: " << tauChannel << endl;
            return 0;
        }
    }
    else {
        cout << "[ERROR] Wrong lepton ID: " << my_lep_ID << endl;
        return 0;
    }
}
//==================================================================
