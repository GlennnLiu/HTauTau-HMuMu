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

#include <HTauTauHMuMu/AnalysisStep/test/Datacards/include/FakeRates.h>
#include <HTauTauHMuMu/AnalysisStep/interface/TauIDSFTool.h>
#include <HTauTauHMuMu/AnalysisStep/test/Datacards/include/Settings.h>
#include <HTauTauHMuMu/AnalysisStep/test/Datacards/include/Datacards.h>

using namespace std;

int main( int argc, char *argv[] )

{
    string path = "/eos/home-g/geliu/LepUni/BigTrees/Simp/";
    
    string allProc[]={"AllData","ggH125","VBFH125","WminusH125","WplusH125","ZH125","bbH125","ggZH125","tHW125","tqH125","ttH125","ZZTo4l","ggTo2e2mu_Contin_MCFM701","ggTo2e2tau_Contin_MCFM701","ggTo2mu2tau_Contin_MCFM701","ggTo4e_Contin_MCFM701","ggTo4mu_Contin_MCFM701","ggTo4tau_Contin_MCFM701","TTWW","TTZZ","WWZ","WZZ","ZZZ","TTZJets_M10_MLM","TTZToLLNuNu_M10","TTZToLL_M1to1O_MLM","Z+X"};
    string Data[]={"AllData"};
    string H[]={"ggH125","VBFH125","WminusH125","WplusH125","ZH125","bbH125","ggZH125","tHW125","tqH125","ttH125"};
    string qqZZ[]={"ZZTo4l"};
    string ggZZ[]={"ggTo2e2mu_Contin_MCFM701","ggTo2e2tau_Contin_MCFM701","ggTo2mu2tau_Contin_MCFM701","ggTo4e_Contin_MCFM701","ggTo4mu_Contin_MCFM701","ggTo4tau_Contin_MCFM701"};
    string ZX[]={"ZX"};
    string rare[]={"TTWW","TTZZ","WWZ","WZZ","ZZZ","TTZJets_M10_MLM","TTZToLLNuNu_M10","TTZToLL_M1to1O_MLM"};
    
    string savepath="/eos/home-g/geliu/LepUni/Datacards/Simp/";

    Datacards *DC = new Datacards();
    
    DC->SetPaths(path,savepath);
    
    DC->ToHistos(H,size(H),Settings::H);
    DC->ToHistos(qqZZ,size(qqZZ),Settings::qqZZ);
    DC->ToHistos(ggZZ,size(ggZZ),Settings::ggZZ);
    DC->ToHistos(rare,size(rare),Settings::rare);
    DC->ToZXHistos(Data,Settings::ZX,"/eos/home-g/geliu/LepUni/FakeRates/withDM_MET/FakeRates_SS.root");
    DC->ToHistos(Data,size(Data),Settings::Data);
    
    DC->ComputeUncValues();
    DC->PrintUncValues();
    DC->ToDatacards();
    
    delete DC;
}