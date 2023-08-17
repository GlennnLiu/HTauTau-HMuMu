#ifndef Plots_h
#define Plots_h

// C++
#include <iostream>
#include <vector>

// ROOT
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"


using namespace std;

class Plots
{
   
public:
   Plots ();
   ~Plots();

   struct LLMass
   {
      TString var_X_label = "m_{LLvis} (GeV)";
      TString var_Y_label = "Events / 4 GeV";
      TString var_cut_label = "";
      Int_t var_N_bin = 30;
      Float_t var_min = 50;
      Float_t var_max = 170;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
      Int_t rebinningDYTTbar = 1;
   };
   
   struct LLGoodMass
   {
      TString var_X_label = "m_{LL} (GeV)";
      TString var_Y_label = "Events / 6 GeV";
      TString var_cut_label = "";
      Int_t var_N_bin = 40;
      Float_t var_min = 50;
      Float_t var_max = 290;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
      Int_t rebinningDYTTbar = 1;
   };

   struct Pzeta
   {
      TString var_X_label = "P_{#zeta} (GeV)";
      TString var_Y_label = "Events / 8 GeV";
      TString var_cut_label = "";
      Int_t var_N_bin = 50;
      Float_t var_min = -200;
      Float_t var_max = 200;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
      Int_t rebinningDYTTbar = 1;
   };

   struct MtLMET
   {
      TString var_X_label = "M_{T,l+MET} (GeV)";
      TString var_Y_label = "Events / 4 GeV";
      TString var_cut_label = "";
      Int_t var_N_bin = 45;
      Float_t var_min = 20;
      Float_t var_max = 200;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
      Int_t rebinningDYTTbar = 1;
   };

   struct MtLLMET
   {
      TString var_X_label = "M_{T,ll+MET} (GeV)";
      TString var_Y_label = "Events / 4 GeV";
      TString var_cut_label = "";
      Int_t var_N_bin = 50;
      Float_t var_min = 40;
      Float_t var_max = 240;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
      Int_t rebinningDYTTbar = 1;
   };

   struct LLPT
   {
      TString var_X_label = "P_{T,LL} (GeV)";
      TString var_Y_label = "Events / 5 GeV";
      TString var_cut_label = "";
      Int_t var_N_bin = 50;
      Float_t var_min = 15;
      Float_t var_max = 265;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
      Int_t rebinningDYTTbar = 1;
   };

private:

};

#endif
