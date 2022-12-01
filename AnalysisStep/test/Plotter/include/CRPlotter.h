#ifndef CRPlotter_h
#define CRPlotter_h

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

#include <HTauTauHMuMu/AnalysisStep/test/Plotter/include/Settings.h>
#include <HTauTauHMuMu/AnalysisStep/interface/FinalStates.h>
#include <HTauTauHMuMu/AnalysisStep/interface/bitops.h>
#include <HTauTauHMuMu/AnalysisStep/test/Plotter/include/FakeRates.h>
#include <HTauTauHMuMu/AnalysisStep/test/Plotter/include/Plots.h>
#include <HTauTauHMuMu/AnalysisStep/interface/TauIDSFTool.h>
#include <HTauTauHMuMu/AnalysisStep/test/Plotter/include/Tree.h>
#include <HTauTauHMuMu/AnalysisStep/test/Plotter/include/CMS_lumi.h>

using namespace std;

const int num_of_final_states      = Settings::num_of_final_states;
const int num_of_zl_final_states   = Settings::num_of_zl_final_states;

class CRPlotter: public Tree
{

public:

    CRPlotter();
    ~CRPlotter();
    
    void SetVariable(string, string, string, string, int, float, float, int);
    void SetRegion(string);
    void SetData(string, string);
    void SetSim(string, string);
    void DeclareHistos();
    void FillHistos();
    void Plotter(bool);
    
private:
    
    double gen_sum_weights, _event_weight;
    float _lumi, _k_factor, _TauIDSF, _yield_SR;
    int _current_process, _current_final_state, _current_region;
    vector<string> _s_final_state;
    vector<string> _s_zl_final_state;
    vector<string> _s_region;
    
    Float_t varFloat;
    Short_t varShort;
    TBranch *b_var;
    TString var, varName, xlabel, ylabel;
    int varType;
    int nbins, range_min, range_max;
    TString region;
    TString pathData, nameData;
    vector<TString> pathSim, nameSim;
    TH1F * hData[3][num_of_final_states];
    TH1F * hSim[10][3][num_of_final_states];
    TString _histo_name, _histo_labels;
    
    TFile *input_file, *output_file;
    TTree *input_tree;
    TH1F *hCounters;
    
    TauIDSFTool *DeepTauSF_VSe_ETau, *DeepTauSF_VSmu_ETau, *DeepTauSF_VSjet_ETau, *DeepTauSF_VSe_MuTau, *DeepTauSF_VSmu_MuTau, *DeepTauSF_VSjet_MuTau, *DeepTauSF_VSe_TauTau, *DeepTauSF_VSmu_TauTau, *DeepTauSF_VSjet_TauTau, *DeepTauSF_VSe_bare, *DeepTauSF_VSmu_bare, *DeepTauSF_VSjet_bare;
    
    void DeclareHistosData();
    void DeclareHistosSim(size_t);
    void FillHistosData();
    void FillHistosSim(size_t);
    void SumHistos();
    void SaveHistos();
    void GetHistos();
    int FindFinalState();
    float calculate_K_factor(TString);
    int GENMatch(int, TString);
    float calculate_TauIDSF_OneLeg(short, float, float, int, int, bool);
    float calculate_TauIDSF(TString, bool);
    void SetColor(TH1F *,int);
    TString ToFSName(string);
    
};
#endif
    