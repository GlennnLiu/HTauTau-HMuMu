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
#include <HTauTauHMuMu/AnalysisStep/test/Plotter/include/CRPlotter.h>

using namespace std;

int main( int argc, char *argv[] )

{
    CRPlotter *Plotter=new CRPlotter();
    
//     Plotter->SetVariable("ZZMass","M_{ZZ}","M_{ZZ} (GeV)","Events / 20 GeV", 14, 70, 350, 2);
//     Plotter->SetVariable("Lep3Pt","p_{t,l3}","p_{t,l3} (GeV)","Events / 10 GeV", 20, 0, 200, 2);
    Plotter->SetVariable("Z2Mass","M_{Z_{prob}}","M_{Z_{prob}} (GeV)","Events / 10 GeV", 14, 0, 140, 2);
    Plotter->SetRegion("CRZLLTree");
    
    string folder="/eos/home-g/geliu/LepUni/BigTrees/Simp/";
    Plotter->SetData(folder+"AllData/HTauTauHMuMu.root","AllData");
    Plotter->SetSim(folder+"DYJetsToLL_M50/HTauTauHMuMu.root","DYJetsToLL");
    Plotter->SetSim(folder+"WZTo3LNu/HTauTauHMuMu.root","WZTo3LNu");
    Plotter->SetSim(folder+"ZZTo4l/HTauTauHMuMu.root","ZZTo4l");
    Plotter->SetSim(folder+"TTTo2L2Nu/HTauTauHMuMu.root","TTTo2L2Nu");

    Plotter->DeclareHistos();
//     Plotter->FillHistos();
    Plotter->Plotter(true);
}