// Include classes
#include <HTauTauHMuMu/AnalysisStep/test/Plotter/include/FakeRates.h>

using namespace std;

// Constructor
//===============================================
FakeRates::FakeRates( TString input_file_FR_name )
{
 
   input_file_FR = new TFile(input_file_FR_name);

   // for OS FR files
   if (input_file_FR_name.Contains("OS"))
   {
      g_FR_mu_EB = (TGraph*)input_file_FR->Get("FR_OS_muon_EB");
      g_FR_mu_EE = (TGraph*)input_file_FR->Get("FR_OS_muon_EE");
      g_FR_e_EB  = (TGraph*)input_file_FR->Get("FR_OS_electron_EB");
      g_FR_e_EE  = (TGraph*)input_file_FR->Get("FR_OS_electron_EE");
      g_FR_tauE_EB = (TGraph*)input_file_FR->Get("FR_OS_tauE_EB");
      g_FR_tauE_EE = (TGraph*)input_file_FR->Get("FR_OS_tauE_EE");
      g_FR_tauMu_EB = (TGraph*)input_file_FR->Get("FR_OS_tauMu_EB");
      g_FR_tauMu_EE = (TGraph*)input_file_FR->Get("FR_OS_tauMu_EE");
      g_FR_tauTau_EB = (TGraph*)input_file_FR->Get("FR_OS_tauTau_EB");
      g_FR_tauTau_EE = (TGraph*)input_file_FR->Get("FR_OS_tauTau_EE");
   }
   
   // for SS FR files
   if (input_file_FR_name.Contains("SS"))
   {
      g_FR_mu_EB = (TGraph*)input_file_FR->Get("FR_SS_muon_EB");
      g_FR_mu_EE = (TGraph*)input_file_FR->Get("FR_SS_muon_EE");
      g_FR_e_EB  = (TGraph*)input_file_FR->Get("FR_SS_electron_EB");
      g_FR_e_EE  = (TGraph*)input_file_FR->Get("FR_SS_electron_EE");
      g_FR_tauE_EB = (TGraph*)input_file_FR->Get("FR_SS_tauE_EB");
      g_FR_tauE_EE = (TGraph*)input_file_FR->Get("FR_SS_tauE_EE");
      g_FR_tauMu_EB = (TGraph*)input_file_FR->Get("FR_SS_tauMu_EB");
      g_FR_tauMu_EE = (TGraph*)input_file_FR->Get("FR_SS_tauMu_EE");
      g_FR_tauTau_EB = (TGraph*)input_file_FR->Get("FR_SS_tauTau_EB");
      g_FR_tauTau_EE = (TGraph*)input_file_FR->Get("FR_SS_tauTau_EE");
   }

}
//===============================================



//======================
FakeRates::~FakeRates() {}
//======================



//==================================================================
float FakeRates::GetFakeRate(float lep_Pt, float lep_eta, int lep_ID, int tauChannel)
{

   float my_lep_Pt = lep_Pt >= 80. ? 79. : lep_Pt;
   int   my_lep_ID = abs(lep_ID);
   //cout<<lep_Pt<<","<<lep_ID<<","<<tauChannel<<endl;
   int bin = 0;
   if ( my_lep_Pt > 5 && my_lep_Pt <= 7 ) bin = 0;
   else if ( my_lep_Pt >  7 && my_lep_Pt <= 10 ) bin = 1;
   else if ( my_lep_Pt > 10 && my_lep_Pt <= 20 ) bin = 2;
   else if ( my_lep_Pt > 20 && my_lep_Pt <= 30 ) bin = 3;
   else if ( my_lep_Pt > 30 && my_lep_Pt <= 40 ) bin = 4;
   else if ( my_lep_Pt > 40 && my_lep_Pt <= 50 ) bin = 5;
   else if ( my_lep_Pt > 50 && my_lep_Pt <= 80 ) bin = 6;
   
   if ( fabs(my_lep_ID) == 11 ) bin = bin-1; // there is no [5, 7] bin in the electron fake rate
   if ( fabs(my_lep_ID) == 15 ) bin = bin-3; // there is no [5, 20] bin in the hadronic tau fake rate
   //cout<<bin<<endl;
   if ( my_lep_ID == 11 ) {
      if ( fabs(lep_eta) < 1.479 ) return (g_FR_e_EB->GetY())[bin];
      else return (g_FR_e_EE->GetY())[bin];
   }
   else if ( my_lep_ID == 13 ) {
      if ( fabs(lep_eta) < 1.2 ) return (g_FR_mu_EB->GetY())[bin];
      else return (g_FR_mu_EE->GetY())[bin];
   }
   else if ( my_lep_ID == 15 ) {
      if ( tauChannel == 0 ) {
	 if ( fabs(lep_eta) < 1.479 ) return (g_FR_tauE_EB->GetY())[bin];
	 else return (g_FR_tauE_EE->GetY())[bin];
      }
      else if ( tauChannel == 1 ) {
         if ( fabs(lep_eta) < 1.479 ) return (g_FR_tauMu_EB->GetY())[bin];
         else return (g_FR_tauMu_EE->GetY())[bin];
      }
      else if ( tauChannel == 2 ) {
         if ( fabs(lep_eta) < 1.479 ) return (g_FR_tauTau_EB->GetY())[bin];
         else return (g_FR_tauTau_EE->GetY())[bin];
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
float FakeRates::GetFakeRate_Up(float lep_Pt, float lep_eta, int lep_ID, int tauChannel)
{

   float my_lep_Pt = lep_Pt >= 80. ? 79. : lep_Pt;
   int   my_lep_ID = abs(lep_ID);

   int bin = 0;
   if ( my_lep_Pt > 5 && my_lep_Pt <= 7 ) bin = 0;
   else if ( my_lep_Pt >  7 && my_lep_Pt <= 10 ) bin = 1;
   else if ( my_lep_Pt > 10 && my_lep_Pt <= 20 ) bin = 2;
   else if ( my_lep_Pt > 20 && my_lep_Pt <= 30 ) bin = 3;
   else if ( my_lep_Pt > 30 && my_lep_Pt <= 40 ) bin = 4;
   else if ( my_lep_Pt > 40 && my_lep_Pt <= 50 ) bin = 5;
   else if ( my_lep_Pt > 50 && my_lep_Pt <= 80 ) bin = 6;
	
   if ( fabs(my_lep_ID) == 11 ) bin = bin-1; // there is no [5, 7] bin in the electron fake rate
   if ( fabs(my_lep_ID) == 15 ) bin = bin-3; // there is no [5, 20] bin in the hadronic tau fake rate

   if ( my_lep_ID == 11 ) {
      if ( fabs(lep_eta) < 1.479 )
         return (g_FR_e_EB->GetY())[bin] + g_FR_e_EB->GetErrorY(bin);
      else
         return (g_FR_e_EE->GetY())[bin] + g_FR_e_EE->GetErrorY(bin);
   }
   else if ( my_lep_ID == 13 ) {
      if ( fabs(lep_eta) < 1.2 )
         return (g_FR_mu_EB->GetY())[bin] + g_FR_mu_EB->GetErrorY(bin);
      else
         return (g_FR_mu_EE->GetY())[bin] + g_FR_mu_EE->GetErrorY(bin);
   }
   else if ( my_lep_ID == 15 ) {
      if ( tauChannel == 0 ) {
         if ( fabs(lep_eta) < 1.479 ) return (g_FR_tauE_EB->GetY())[bin] + g_FR_tauE_EB->GetErrorY(bin);
         else return (g_FR_tauE_EE->GetY())[bin] + g_FR_tauE_EE->GetErrorY(bin);
      }
      else if ( tauChannel == 1 ) {
         if ( fabs(lep_eta) < 1.479 ) return (g_FR_tauMu_EB->GetY())[bin] + g_FR_tauMu_EB->GetErrorY(bin);
         else return (g_FR_tauMu_EE->GetY())[bin] + g_FR_tauMu_EE->GetErrorY(bin);
      }
      else if ( tauChannel == 2 ) {
         if ( fabs(lep_eta) < 1.479 ) return (g_FR_tauTau_EB->GetY())[bin] + g_FR_tauTau_EB->GetErrorY(bin);
         else return (g_FR_tauTau_EE->GetY())[bin] + g_FR_tauTau_EE->GetErrorY(bin);
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
float FakeRates::GetFakeRate_Dn(float lep_Pt, float lep_eta, int lep_ID, int tauChannel)
{

   float my_lep_Pt = lep_Pt >= 80. ? 79. : lep_Pt;
   int   my_lep_ID = abs(lep_ID);

   int bin = 0;
   if ( my_lep_Pt > 5 && my_lep_Pt <= 7 ) bin = 0;
   else if ( my_lep_Pt >  7 && my_lep_Pt <= 10 ) bin = 1;
   else if ( my_lep_Pt > 10 && my_lep_Pt <= 20 ) bin = 2;
   else if ( my_lep_Pt > 20 && my_lep_Pt <= 30 ) bin = 3;
   else if ( my_lep_Pt > 30 && my_lep_Pt <= 40 ) bin = 4;
   else if ( my_lep_Pt > 40 && my_lep_Pt <= 50 ) bin = 5;
   else if ( my_lep_Pt > 50 && my_lep_Pt <= 80 ) bin = 6;
	
   if ( fabs(my_lep_ID) == 11 ) bin = bin-1; // there is no [5, 7] bin in the electron fake rate
   if ( fabs(my_lep_ID) == 15 ) bin = bin-3; // there is no [5, 20] bin in the hadronic tau fake rate

   if ( my_lep_ID == 11 ) {
      if ( fabs(lep_eta) < 1.479 )
         return (g_FR_e_EB->GetY())[bin] - g_FR_e_EB->GetErrorY(bin);
      else
         return (g_FR_e_EE->GetY())[bin] - g_FR_e_EE->GetErrorY(bin);
   }
   else if ( my_lep_ID == 13 ) {
      if ( fabs(lep_eta) < 1.2 )
         return (g_FR_mu_EB->GetY())[bin] - g_FR_mu_EB->GetErrorY(bin);
      else
         return (g_FR_mu_EE->GetY())[bin] - g_FR_mu_EE->GetErrorY(bin);
   }
   else if ( my_lep_ID == 15 ) {
      if ( tauChannel == 0 ) {
         if ( fabs(lep_eta) < 1.479 ) return (g_FR_tauE_EB->GetY())[bin] - g_FR_tauE_EB->GetErrorY(bin);
         else return (g_FR_tauE_EE->GetY())[bin] - g_FR_tauE_EE->GetErrorY(bin);
      }
      else if ( tauChannel == 1 ) {
         if ( fabs(lep_eta) < 1.479 ) return (g_FR_tauMu_EB->GetY())[bin] - g_FR_tauMu_EB->GetErrorY(bin);
         else return (g_FR_tauMu_EE->GetY())[bin] - g_FR_tauMu_EE->GetErrorY(bin);
      }
      else if ( tauChannel == 2 ) {
         if ( fabs(lep_eta) < 1.479 ) return (g_FR_tauTau_EB->GetY())[bin] - g_FR_tauTau_EB->GetErrorY(bin);
         else return (g_FR_tauTau_EE->GetY())[bin] - g_FR_tauTau_EE->GetErrorY(bin);
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
