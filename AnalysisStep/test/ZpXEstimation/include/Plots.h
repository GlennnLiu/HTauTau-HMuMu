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
   
   //=============
   // M4l
   //=============
   struct M4l
   {
      TString var_X_label = "m_{4#font[12]{l}} (GeV)";
      TString var_Y_label = "Events / 20 GeV";
      TString var_cut_label = "";
      Int_t var_N_bin = 40;
      Float_t var_min = 70;
      Float_t var_max = 870;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
      Int_t rebinningDYTTbar = 1;
   };

   struct ZtagMass
   {
      TString var_X_label = "m_{Z_{#font[12]{tag}}} (GeV)";
      TString var_Y_label = "Events / GeV";
      TString var_cut_label = "";
      Int_t var_N_bin = 20;
      Float_t var_min = 81.1876;
      Float_t var_max = 101.1876;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
      Int_t rebinningDYTTbar = 1;
   };
   
   struct ZprobMass
   {
      TString var_X_label = "m_{Z_{#font[12]{prob}}} (GeV)";
      TString var_Y_label = "Events / 4 GeV";
      TString var_cut_label = "";
      Int_t var_N_bin = 35;
      Float_t var_min = 0;
      Float_t var_max = 140;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
      Int_t rebinningDYTTbar = 1;
   };

   struct ZprobPt
   {
      TString var_X_label = "P_{t,Z_{#font[12]{prob}}} (GeV)";
      TString var_Y_label = "Events / 10 GeV";
      TString var_cut_label = "";
      Int_t var_N_bin = 30;
      Float_t var_min = 0;
      Float_t var_max = 300;
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
