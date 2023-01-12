#ifndef Datacards_h
#define Datacards_h

// C++
#include <iostream>
#include <fstream>
#include <iomanip> // For setprecision
#include <vector>
#include <map>

// ROOT
#include "TApplication.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TLorentzVector.h"
#include "THStack.h"

#include <HTauTauHMuMu/AnalysisStep/interface/FinalStates.h>
#include <HTauTauHMuMu/AnalysisStep/interface/bitops.h>
#include <HTauTauHMuMu/AnalysisStep/test/Datacards/include/FakeRates.h>
#include <HTauTauHMuMu/AnalysisStep/interface/TauIDSFTool.h>
#include <HTauTauHMuMu/AnalysisStep/test/Datacards/include/Tree.h>
#include <HTauTauHMuMu/AnalysisStep/test/Datacards/include/Settings.h>

using namespace std;

const int num_of_processes         = Settings::num_of_processes;
const int num_of_final_states      = Settings::num_of_final_states;
const int num_of_gen_final_states  = Settings::num_of_gen_final_states;
const int num_of_norm_unc          = 5;
const int num_of_shape_unc         = 8;
const int num_of_unc               = 13;

class Datacards: public Tree
{

public:
    Datacards();
    ~Datacards();
    
    void SetPaths(string,string);
    void ToHistos(string*,int,int);
    void ToZXHistos(string*,int,TString);
    void ComputeUncValues();
    void PrintUncValues();
    void ToDatacards();
    
private:
    void DeclareHistos(int);
    int FindFinalState();
    int FindGenFinalState();
    float calculate_K_factor(TString, const string&);
    float calculate_TauIDSF_OneLeg(short, float, float, int, int, const string&);
    float calculate_TauIDSF(TString, const string&);
    int GENMatch(int,TString);
    
    string normUnc[num_of_norm_unc]={"CMS_lumi","FR","qqZZ_K","ggZZ_K","pythia_scale"}, shapeUnc[num_of_shape_unc]={"eff_e","eff_mu","ID_tau","tes","ees","mes","jes","jer"}, Unc[num_of_unc]={"CMS_lumi","FR","qqZZ_K","ggZZ_K","pythia_scale","eff_e","eff_mu","ID_tau","tes","ees","mes","jes","jer"};
    float UncValues[num_of_processes-1][num_of_gen_final_states-1][num_of_final_states-1][num_of_unc][2];
    
    string _path, _savepath;
    vector<string> _s_process, _s_final_state, _s_gen_final_state;
    int _current_process, _current_final_state, _current_gen_final_state;
    float _lumi, _k_factor, _TauIDSF, _yield_SR;
    double gen_sum_weights, _event_weight, _event_weight_up, _event_weight_dn;
    vector<float> _fs_ROS_SS, _cb_ss, _FR_uncUp, _FR_uncDn;
    
    TFile *input_file, *output_file;
    TTree *input_tree;
    TH1F *hCounters;
    
    TString _histo_name, _histo_labels;
    TH1F *histos_nom[num_of_processes][num_of_gen_final_states][num_of_final_states];
    TH1F *histos_unc[num_of_processes-1][num_of_gen_final_states][num_of_final_states][num_of_unc][2/*0:up,1:down*/];
    TH1F *histos_unc_ggzz[num_of_gen_final_states][num_of_final_states][4][2];
    
    TauIDSFTool *DeepTauSF_VSe_ETau, *DeepTauSF_VSmu_ETau, *DeepTauSF_VSjet_ETau, *DeepTauSF_VSe_MuTau, *DeepTauSF_VSmu_MuTau, *DeepTauSF_VSjet_MuTau, *DeepTauSF_VSe_TauTau, *DeepTauSF_VSmu_TauTau, *DeepTauSF_VSjet_TauTau;
};
#endif