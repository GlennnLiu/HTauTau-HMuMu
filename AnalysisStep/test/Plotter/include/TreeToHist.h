#ifndef TreeToHist_h
#define TreeToHist_h

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

const int num_of_processes         = Settings::num_of_processes;
const int num_of_final_states      = Settings::num_of_final_states;
const int num_of_gen_final_states      = Settings::num_of_gen_final_states;

class TreeToHist: public Tree
{

public:
    TreeToHist();
    ~TreeToHist();
    
    void SetPaths(string,string,string);
    void ToHistos(string*,int,int,bool);
    void ToHistosZX(string*,int,bool,TString);
    void GetHistos(string*,int,int,bool);
    void SumTotalMC();
    void ToPlots(bool);
    void ToPlotsMC(bool);
    
private:
    void DeclareHistos(string*,int,int);
    int FindFinalState();
    int FindGenFinalState();
    float calculate_K_factor(TString);
    int GENMatch(int, TString );
    float calculate_TauIDSF_OneLeg(short, float, float, int, int);
    float calculate_TauIDSF(TString);
    void SumGroups(int,int);
    void saveHistos(int,int);
    void SetColor(TH1F*,int);
    TString ToFSName(string);
    TString ToGFSName(string);

    string _path, _file_name, _savepath;
    vector<string> _pathList[num_of_processes];
    vector<string> _s_process, _s_final_state, _s_gen_final_state;
    int _current_process, _current_final_state, _current_gen_final_state;
    float _lumi, _k_factor, _TauIDSF, _yield_SR;
    double gen_sum_weights, _event_weight;
    vector<float> _fs_ROS_SS, _cb_ss;
    
    
    TFile *input_file, *output_file;
    TTree *input_tree;
    TH1F *hCounters;
    
    TH1F *histos_1D[num_of_processes][num_of_final_states][num_of_gen_final_states][99][4];
    TString _histo_name, _histo_labels;
    
    TauIDSFTool *DeepTauSF_VSe_ETau, *DeepTauSF_VSmu_ETau, *DeepTauSF_VSjet_ETau, *DeepTauSF_VSe_MuTau, *DeepTauSF_VSmu_MuTau, *DeepTauSF_VSjet_MuTau, *DeepTauSF_VSe_TauTau, *DeepTauSF_VSmu_TauTau, *DeepTauSF_VSjet_TauTau;
    
};
#endif