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
#include <HTauTauHMuMu/AnalysisStep/test/SuperRatio/include/Settings.h>
#include <HTauTauHMuMu/AnalysisStep/test/SuperRatio/include/LL.h>

using namespace std;

int main( int argc, char *argv[] )

{
    string path = "/eos/home-g/geliu/LepUni/BigTrees/UL/2016/";
    string file_name = "/HTauTauHMuMu_new.root";
    
    string data[]={"SingleEle2016","SingleMuon2016","DoubleMu2016","MuonEG2016","Tau2016"};
    string bkg[]={"DYJetsToLL","EWKZToLL","TTTo2L2Nu","TTW","TTWW","TTZToLLNuNu_M10","TTZToLL_M1to10","TTZZ","WWTo2L2Nu","WWW","WWZ","WZTo2L2Q","WZTo3LNu","WZZ","ZZTo2L2Nu","ZZTo2L2Q","ZZTo4l","ZZZ","ggTo2e2mu","ggTo2e2nu","ggTo2e2tau","ggTo2mu2nu","ggTo2mu2tau","ggTo4e","ggTo4mu","ggTo4tau","tW_antitop","tW_top","tZq","TTToHadronic1","TTToHadronic2","TTToSemiLeptonic1","TTToSemiLeptonic2","TTToSemiLeptonic3","ST_antitop","ST_s","ST_top"};
    string tt[]={"TTToHadronic1","TTToHadronic2","TTToSemiLeptonic1","TTToSemiLeptonic2","TTToSemiLeptonic3"};
    
    string savepath="/eos/home-g/geliu/LepUni/SuperRatio/";

    float l1pt_bins[]={15,20,30,50,80,150};
    float l2pt_bins[]={15,20,30,50,80,150};
    float Mvis_bins[]={20,60,90,120,150,200,300};
    float DR_bins[]={0,1,1.5,2,2.5,3,3.5,4,4.5,5,6};
    
    LL *ll = new LL();
    
    ll->SetPaths(path,file_name,savepath);

    ll->SetFileList(data,size(data),bkg,size(bkg),tt,size(tt));

    ll->Set_lpt_bin(size(l1pt_bins),l1pt_bins,1);
    ll->Set_lpt_bin(size(l2pt_bins),l2pt_bins,2);
    ll->Set_Mvis_bin(size(Mvis_bins),Mvis_bins);
    ll->Set_DR_bin(size(DR_bins),DR_bins);

    cout<<"---------------------------------------------------"<<endl;
    cout<<"----------Step 1: compute fake rates---------------"<<endl;
    cout<<"---------------------------------------------------"<<endl;
    // ll->Step1_FakeRate_DeclareHistos();
    // ll->Step1_FakeRate_FillHistos();
    // ll->Step1_FakeRate_Compute();
    // ll->Step1_FakeRate_SaveHistos();
    // ll->Step1_FakeRate_GetHistos();

    cout<<"---------------------------------------------------"<<endl;
    cout<<"----------Step 2: compute closure factors----------"<<endl;
    cout<<"---------------------------------------------------"<<endl;
    // ll->Step2_Closure_DeclareHistos();
    // ll->Step2_Closure_FillHistos();
    // ll->Step2_Closure_Compute();
    // ll->Step2_Closure_SaveHistos();
    // ll->Step2_Closure_GetHistos();

    cout<<"---------------------------------------------------"<<endl;
    cout<<"----------Step 3.1: Diff vs SR for QCD-------------"<<endl;
    cout<<"---------------------------------------------------"<<endl;
    // ll->Step3_QCD_FakeRate_DeclareHistos();
    // ll->Step3_QCD_FakeRate_FillHistos();
    // ll->Step3_QCD_FakeRate_Compute();
    // ll->Step3_QCD_FakeRate_SaveHistos();
    // ll->Step3_QCD_FakeRate_GetHistos();

    // ll->Step3_QCD_Closure_DeclareHistos();
    // ll->Step3_QCD_Closure_FillHistos();
    // ll->Step3_QCD_Closure_Compute();
    // ll->Step3_QCD_Closure_SaveHistos();
    // ll->Step3_QCD_Closure_GetHistos();

    // ll->Step3_QCD_vsSR_DeclareHistos();
    // ll->Step3_QCD_vsSR_FillHistos();
    // ll->Step3_QCD_vsSR_Compute();
    // ll->Step3_QCD_vsSR_SaveHistos();
    // ll->Step3_QCD_vsSR_GetHistos();

    cout<<"---------------------------------------------------"<<endl;
    cout<<"----------Step 4: Compute the fractions------------"<<endl;
    cout<<"---------------------------------------------------"<<endl;
    ll->Step4_Fraction_DeclareHistos();
    ll->Step4_Fraction_FillHistos();
    ll->Step4_Fraction_Compute();
    ll->Step4_Fraction_SaveHistos();
    ll->Step4_Fraction_GetHistos();
}
