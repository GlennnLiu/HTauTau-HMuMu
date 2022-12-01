#ifndef FAKERATES_H
#define FAKERATES_H

// C++
#include <iostream>
#include <fstream>

// ROOT
#include "TFile.h"
#include "TString.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace std;

class FakeRates
{

public:
	
	FakeRates( TString );
	~FakeRates();
   float GetFakeRate( float, float, float, int, int, int );
   float GetFakeRate_Up( float, float, float, int, int, int );
   float GetFakeRate_Dn( float, float, float, int, int, int );
   
   private:
   
   TFile *input_file_FR;
   
   TH1F *g_FR_mu_EB, *g_FR_mu_EE, *g_FR_e_EB, *g_FR_e_EE;
   TH2F *g_FR_tauE_Decay0, *g_FR_tauE_Decay1, *g_FR_tauE_Decay10, *g_FR_tauE_Decay11, *g_FR_tauMu_Decay0, *g_FR_tauMu_Decay1, *g_FR_tauMu_Decay10, *g_FR_tauMu_Decay11, *g_FR_tauTau_Decay0, *g_FR_tauTau_Decay1, *g_FR_tauTau_Decay10, *g_FR_tauTau_Decay11;

};

#endif
