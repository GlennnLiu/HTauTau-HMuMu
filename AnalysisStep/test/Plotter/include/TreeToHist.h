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
// #include <HTauTauHMuMu/AnalysisStep/test/Plotter/include/FakeRates.h>
#include <HTauTauHMuMu/AnalysisStep/test/Plotter/include/Plots.h>
// #include <HTauTauHMuMu/AnalysisStep/interface/TauIDSFTool.h>
#include <HTauTauHMuMu/AnalysisStep/test/Plotter/include/Tree.h>
#include <HTauTauHMuMu/AnalysisStep/test/Plotter/include/CMS_lumi.h>
#include <HTauTauHMuMu/AnalysisStep/test/Plotter/include/FakeBkg.h>

using namespace std;

const int num_of_processes         = Settings::num_of_processes;
const int num_of_final_states      = Settings::num_of_final_states;
const int num_of_categories      = Settings::num_of_categories;
const int num_of_variables       = 6;

class TreeToHist: public Tree
{

public:
    TreeToHist();
    ~TreeToHist();
    
    void SetPaths(string,string,string,string);
    void ToHistos(string*,int,int);
    void ToHistosFake(string*, int, string*, int, int);
    void GetHistos(int);
    void SumTotalMC(bool);
    void ToPlots(bool,bool,bool);
    void ToPlotsRatio(bool,bool);
    
private:
    void DeclareHistos(int);
    int FindFinalState();
    bool pass_Extra_Cuts(int,float,float);
    int FindProcess(int);
    int FindCategory(int,bool);
    float calculate_K_factor(TString);
    void SumGroups(int);
    void saveHistos(int);
    void SetColor(TH1F*,int);
    TString ToFSName(string);
    TString ToProcessName(string);
    int NumberOfJets();

    string _path, _file_name, _savepath;
    TString _directory;
    vector<string> _pathList[num_of_processes];
    vector<string> _s_process, _s_final_state, _s_category;
    int _current_process, _current_final_state, _current_category;
    float _lumi, _k_factor, _TauIDSF;
    double gen_sum_weights, _event_weight;
    
    TFile *input_file, *output_file;
    TTree *input_tree;
    TH1F *hCounters;
    
    TH1F *histos_1D[num_of_processes][num_of_final_states][num_of_categories][num_of_variables];
    TString _histo_name, _histo_labels;

    FakeBkg* fakeBkg;
};
#endif