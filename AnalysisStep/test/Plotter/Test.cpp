// C++
#include <iostream>
#include <fstream>
#include <string>

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

#include <HTauTauHMuMu/AnalysisStep/test/ZpXEstimation/include/FakeRates.h>
#include <HTauTauHMuMu/AnalysisStep/interface/TauIDSFTool.h>
#include <HTauTauHMuMu/AnalysisStep/test/Plotter/include/Settings.h>
#include <HTauTauHMuMu/AnalysisStep/test/Plotter/include/TreeToHist.h>

using namespace std;

int main( int argc, char *argv[] )

{
    string path = "/eos/home-g/geliu/LepUni/BigTrees/Simp/";
    string file_name = "/HTauTauHMuMu_new.root";
    
    string allProc[]={"AllData","ggH125","VBFH125","WminusH125","WplusH125","ZH125","bbH125","ggZH125","tHW125","tqH125","ttH125","ZZTo4l","ggTo2e2mu_Contin_MCFM701","ggTo2e2tau_Contin_MCFM701","ggTo2mu2tau_Contin_MCFM701","ggTo4e_Contin_MCFM701","ggTo4mu_Contin_MCFM701","ggTo4tau_Contin_MCFM701","TTWW","TTZZ","WWZ","WZZ","ZZZ","TTZJets_M10_MLM","TTZToLLNuNu_M10","TTZToLL_M1to1O_MLM","Z+X"};
    string Data[]={"DataF"};
    string H[]={"ggH125","VBFH125","WminusH125","WplusH125","ZH125","bbH125","ggZH125","tHW125","tqH125","ttH125"};
    string qqZZ[]={"ZZTo4l"};
    string ggZZ[]={"ggTo2e2mu_Contin_MCFM701","ggTo2e2tau_Contin_MCFM701","ggTo2mu2tau_Contin_MCFM701","ggTo4e_Contin_MCFM701","ggTo4mu_Contin_MCFM701","ggTo4tau_Contin_MCFM701"};
    string ZX[]={"DYJetsToLL_M50","TTTo2L2Nu","WZTo2L2Q","ZZTo2L2Q","WJetsToLNu"};
    string WZ[]={"WZTo3LNu"};
    string rare[]={"TTWW","TTZZ","WWZ","WZZ","ZZZ","TTZJets_M10_MLM","TTZToLLNuNu_M10","TTZToLL_M1to1O_MLM"};
    //string ZX[]={"Z+X"};
    
    string savepath="/eos/home-g/geliu/LepUni/Histos/Simp/Test/";
    
    TreeToHist *TH = new TreeToHist();
    
    TH->SetPaths(path,file_name,savepath);
    
    TH->ToHistosCR(Data,Settings::ZX,false,"/eos/home-g/geliu/LepUni/FakeRates/Test/FakeRates_SS.root");
}