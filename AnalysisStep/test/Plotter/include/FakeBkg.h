#ifndef FakeBkg_h
#define FakeBkg_h

// C++
#include <iostream>
#include <fstream>
#include <iomanip> // For setprecision
#include <vector>
#include <map>

// ROOT
#include <TFile.h>   // TFile
#include <TH1.h>     // TH1
#include <TF1.h>     // TF1
#include <TH2.h>
#include <TString.h> // Form
#include <string>    // std::string
#include <vector>    // std::vector
#include <map>       // std::map
#include <stdlib.h>  // getenv
#include <functional>

using namespace std;

class FakeBkg {

    protected:

    TFile *f_LTau_FR, *f_LTau_Closure_FR, *f_LTau_QCDvsSR_FR, *f_LTau_WJvsSR, *f_LTau_Frac;
    const TF1 *LTau_f_FR[3/*#fakeBkg*/][2/*#finalState*/][3/*#njetBin*/];
    TH1F *LTau_hClosure_FR[3/*#fakeBkg*/][2/*#finalState*/][3/*#ntauptBin*/];
    const TF1 *LTau_fClosure_FR[3/*#fakeBkg*/][2/*#finalState*/][3/*#ntauptBin*/];
    const TF1 *LTau_fQCDvsSR_FR[2/*#finalState*/];
    const TF1 *LTau_fWJvsSR_FR[2/*#finalState*/];
    TH1F *LTau_hFrac_FR[3/*#fakeBkg*/][2/*#finalState*/][6/*#category*/];

    TFile *f_TauTau_FR, *f_TauTau_Closure_FR, *f_TauTau_QCDvsSR_FR;
    const TF1 *TauTau_f_FR[3/*#njetBin*/];
    const TF1 *TauTau_fClosure_FR;
    const TF1 *TauTau_fQCDvsSR_FR;

    TFile *f_LL_FR, *f_LL_Closure_FR, *f_LL_QCDvsSR_FR, *f_LL_Frac;
    const TF1 *LL_f_FR[2/*#fakeBkg*/];
    TH2F *LL_hClosure_FR[2/*#fakeBkg*/];
    TH2F *LL_hQCDvsSR_FR;
    TH1F *LL_hFrac_FR[2/*#fakeBkg*/][6/*#category*/];

    void disabled() const{
      std::cerr << std::endl << "ERROR! Method has been disabled!"<< std::endl;
      assert(0);
    }

    public:

    FakeBkg(const std::string& year, const std::string& path);
    ~FakeBkg() {}

    float get_LTau_FR(int fs, int cat, int njet, float lpt, float taupt, float DR, float MT, float Mvis);
    float get_TauTau_FR(int njet, float tau1pt, float tau2pt, float Mvis);
    float get_LL_FR(int cat, float DR, float l1pt, float l2pt, float Mvis);

    private:
    vector<string>  _s_final_state, _s_categories, _s_fake_bkg;
    string _histo_name;


};

#endif

