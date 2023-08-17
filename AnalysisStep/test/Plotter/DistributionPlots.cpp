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

// #include <HTauTauHMuMu/AnalysisStep/test/Plotter/include/FakeRates.h>
// #include <HTauTauHMuMu/AnalysisStep/interface/TauIDSFTool.h>
#include <HTauTauHMuMu/AnalysisStep/test/Plotter/include/Settings.h>
#include <HTauTauHMuMu/AnalysisStep/test/Plotter/include/TreeToHist.h>

using namespace std;

int main( int argc, char *argv[] )

{
    string path = "/eos/home-g/geliu/LepUni/BigTrees/UL/2016/";
    string file_name = "/HTauTauHMuMu_new.root";
    
    string directory;
    if (argc>=3 && std::string(argv[1])=="--region") directory = argv[2];
    else directory = "SRTree";

    bool addData;
    if (argc>=4 && std::string(argv[3])=="--addData") addData=true;

    bool addSignal=false;
    bool setLog=false;
    if (directory=="SRTree") {addSignal=true; setLog=true;}
    
    // string allProc[]={"AllData","ggH125","VBFH125","WminusH125","WplusH125","ZH125","bbH125","ggZH125","tHW125","tqH125","ttH125","ZZTo4l","ggTo2e2mu_Contin_MCFM701","ggTo2e2tau_Contin_MCFM701","ggTo2mu2tau_Contin_MCFM701","ggTo4e_Contin_MCFM701","ggTo4mu_Contin_MCFM701","ggTo4tau_Contin_MCFM701","TTWW","TTZZ","WWZ","WZZ","ZZZ","TTZJets_M10_MLM","TTZToLLNuNu_M10","TTZToLL_M1to1O_MLM","Z+X"};
    string data[]={"SingleEle2016","SingleMuon2016","DoubleMu2016","MuonEG2016","Tau2016"};
    string Htt[]={"ggH125ToTauTau","VBFH125ToTauTau","WplusH125ToTauTau","WminusH125ToTauTau","ZH125ToTauTau","ttH125ToTauTau"};
    string Hmm[]={"ggH125ToMuMu","VBFH125ToMuMu","WplusH125ToMuMu","WminusH125ToMuMu","ZH125ToMuMu","ttH125ToMuMu"};
    string Hww[]={"ggH125ToWW","VBFH125ToWW","WplusH125ToWW","WminusH125ToWW","ZH125ToWW","ggZH125ToWW"};
    string diBoson[]={"ZZTo4l","ZZTo2L2Q","ZZTo2L2Nu","WZTo2L2Q","WZTo3LNu","WWTo2L2Nu","ggTo2e2mu","ggTo2e2tau","ggTo2mu2tau","ggTo4e","ggTo4mu","ggTo4tau","ggTo2e2nu","ggTo2mu2nu"};
    string triBoson[]={"ZZZ","WZZ","WWZ","WWW"};
    string top[]={"tW_top","tW_antitop","tZq","TTTo2L2Nu","TTW","TTZZ","TTWW","TTZToLLNuNu_M10","TTZToLL_M1to10","ST_antitop","ST_s","ST_top"};
    string EWKZ[]={"EWKZToLL"};
    string DY[]={"DYJetsToLL"};
    string fake[]={"WJetsToLNu","TTToHadronic1","TTToHadronic2","TTToSemiLeptonic1","TTToSemiLeptonic2","TTToSemiLeptonic3"};
    string bkgSub[]={"DYJetsToLL","EWKZToLL","TTTo2L2Nu","TTW","TTWW","TTZToLLNuNu_M10","TTZToLL_M1to10","TTZZ","WWTo2L2Nu","WWW","WWZ","WZTo2L2Q","WZTo3LNu","WZZ","ZZTo2L2Nu","ZZTo2L2Q","ZZTo4l","ZZZ","ggTo2e2mu","ggTo2e2nu","ggTo2e2tau","ggTo2mu2nu","ggTo2mu2tau","ggTo4e","ggTo4mu","ggTo4tau","tW_antitop","tW_top","tZq","ST_antitop","ST_s","ST_top"};
    
    string savepath="/eos/home-g/geliu/LepUni/Histos/UL2016/DistributionPlots/";
    
    TreeToHist *TH = new TreeToHist();
    
    TH->SetPaths(path,file_name,savepath,directory);
    
    if (addSignal) {
        TH->ToHistos(Htt,size(Htt),Settings::Htt);
        TH->ToHistos(Hmm,size(Hmm),Settings::Hmm);
        TH->ToHistos(Hww,size(Hww),Settings::Hww);
    }
    TH->ToHistos(diBoson,size(diBoson),Settings::diBoson);
    TH->ToHistos(triBoson,size(triBoson),Settings::triBoson);
    TH->ToHistos(top,size(top),Settings::top);
    TH->ToHistos(EWKZ,size(EWKZ),Settings::EWKZee);
    if (directory=="SRTree") {
        TH->ToHistosFake(data,size(data),bkgSub,size(bkgSub),Settings::fake);
    }
    else {
        TH->ToHistos(fake,size(fake),Settings::fake);
    }
    TH->ToHistos(DY,size(DY),Settings::DYee);

    if (addData) {
        TH->ToHistos(data,size(data),Settings::Data);
    }


    //TH->SumTotalMC();

    if (addSignal) {
        TH->GetHistos(Settings::Htt);
        TH->GetHistos(Settings::Hmm);
        TH->GetHistos(Settings::Hww);
    }
    TH->GetHistos(Settings::diBoson);
    TH->GetHistos(Settings::triBoson);
    TH->GetHistos(Settings::top);
    TH->GetHistos(Settings::EWKZee);
    TH->GetHistos(Settings::fake);
    TH->GetHistos(Settings::DYee);
    if (addData) {
        TH->GetHistos(Settings::Data);
    }

    TH->SumTotalMC(addSignal);
    TH->ToPlotsRatio(addSignal,setLog);
}
