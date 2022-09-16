#ifndef FAKERATES_H
#define FAKERATES_H

// C++
#include <iostream>
#include <fstream>

// ROOT
#include "TFile.h"
#include "TString.h"
#include "TGraph.h"

using namespace std;

class FakeRates
{

public:
	
	FakeRates( TString );
	~FakeRates();
   float GetFakeRate( float, float, int, int );
   float GetFakeRate_Up( float, float, int, int );
   float GetFakeRate_Dn( float, float, int, int );
   
   private:
   
   TFile *input_file_FR;
   
   TGraph *g_FR_mu_EB, *g_FR_mu_EE, *g_FR_e_EB, *g_FR_e_EE, *g_FR_tauE_EB, *g_FR_tauE_EE, *g_FR_tauMu_EB, *g_FR_tauMu_EE, *g_FR_tauTau_EB, *g_FR_tauTaub;

};

#endif
