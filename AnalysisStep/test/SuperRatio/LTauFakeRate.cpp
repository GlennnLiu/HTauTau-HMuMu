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
#include <HTauTauHMuMu/AnalysisStep/test/SuperRatio/include/LTau.h>

using namespace std;

int main( int argc, char *argv[] )

{
    string path = "/eos/home-g/geliu/LepUni/BigTrees/UL/2016/";
    string file_name = "/HTauTauHMuMu_new.root";
    
    string data[]={"SingleEle2016","SingleMuon2016","DoubleMu2016","MuonEG2016","Tau2016"};
    string bkg[]={"DYJetsToLL","EWKZToLL","TTTo2L2Nu","TTW","TTWW","TTZToLLNuNu_M10","TTZToLL_M1to10","TTZZ","WJetsToLNu","WWTo2L2Nu","WWW","WWZ","WZTo2L2Q","WZTo3LNu","WZZ","ZZTo2L2Nu","ZZTo2L2Q","ZZTo4l","ZZZ","ggTo2e2mu","ggTo2e2nu","ggTo2e2tau","ggTo2mu2nu","ggTo2mu2tau","ggTo4e","ggTo4mu","ggTo4tau","tW_antitop","tW_top","tZq","TTToHadronic1","TTToHadronic2","TTToSemiLeptonic1","TTToSemiLeptonic2","TTToSemiLeptonic3","ST_antitop","ST_s","ST_top"};
    string tt[]={"TTToHadronic1","TTToHadronic2","TTToSemiLeptonic1","TTToSemiLeptonic2","TTToSemiLeptonic3"};
    string Wjet[]={"WJetsToLNu"};
    
    string savepath="/eos/home-g/geliu/LepUni/SuperRatio/";

    float taupt_bins[]={30,35,40,45,50,55,60,70,80,100,150};
    float lpt_bins[]={20,23,26,30,35,40,45,50,55,60,70,80,90,100,150};
    float Mvis_bins[]={20,60,90,120,150,200,300};
    float MT_bins[]={20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180};
    float DR_bins[]={0,1,1.5,2,2.5,3,3.5,4,4.5,5,6};
    
    LTau *ltau = new LTau();
    
    ltau->SetPaths(path,file_name,savepath);

    ltau->SetFileList(data,size(data),bkg,size(bkg),Wjet,size(Wjet),tt,size(tt));

    ltau->Set_taupt_bin(size(taupt_bins),taupt_bins);
    ltau->Set_lpt_bin(size(lpt_bins),lpt_bins);
    ltau->Set_Mvis_bin(size(Mvis_bins),Mvis_bins);
    ltau->Set_MT_bin(size(MT_bins),MT_bins);
    ltau->Set_DR_bin(size(DR_bins),DR_bins);

    cout<<"---------------------------------------------------"<<endl;
    cout<<"----------Step 1: compute fake rates---------------"<<endl;
    cout<<"---------------------------------------------------"<<endl;
    ltau->Step1_FakeRate_DeclareHistos();
    ltau->Step1_FakeRate_FillHistos();
    ltau->Step1_FakeRate_Compute();
    ltau->Step1_FakeRate_SaveHistos();
    ltau->Step1_FakeRate_GetHistos();

    cout<<"---------------------------------------------------"<<endl;
    cout<<"----------Step 2: compute closure factors----------"<<endl;
    cout<<"---------------------------------------------------"<<endl;
    ltau->Step2_Closure_DeclareHistos();
    ltau->Step2_Closure_FillHistos();
    ltau->Step2_Closure_Compute();
    ltau->Step2_Closure_SaveHistos();

    cout<<"---------------------------------------------------"<<endl;
    cout<<"----------Step 3.1: Diff vs SR for QCD-------------"<<endl;
    cout<<"---------------------------------------------------"<<endl;
    // ltau->Step3_QCD_FakeRate_DeclareHistos();
    // ltau->Step3_QCD_FakeRate_FillHistos();
    // ltau->Step3_QCD_FakeRate_Compute();
    // ltau->Step3_QCD_FakeRate_SaveHistos();
    // ltau->Step3_QCD_FakeRate_GetHistos();

    // ltau->Step3_QCD_Closure_DeclareHistos();
    // ltau->Step3_QCD_Closure_FillHistos();
    // ltau->Step3_QCD_Closure_Compute();
    // ltau->Step3_QCD_Closure_SaveHistos();
    // ltau->Step3_QCD_Closure_GetHistos();

    // ltau->Step3_QCD_vsSR_DeclareHistos();
    // ltau->Step3_QCD_vsSR_FillHistos();
    // ltau->Step3_QCD_vsSR_Compute();
    // ltau->Step3_QCD_vsSR_SaveHistos();
    // ltau->Step3_QCD_vsSR_GetHistos();

    cout<<"---------------------------------------------------"<<endl;
    cout<<"----------Step 3.1: Diff vs SR for WJ-------------"<<endl;
    cout<<"---------------------------------------------------"<<endl;
    // ltau->Step3_WJ_FakeRate_DeclareHistos();
    // ltau->Step3_WJ_FakeRate_FillHistos();
    // ltau->Step3_WJ_FakeRate_Compute();
    // ltau->Step3_WJ_FakeRate_SaveHistos();
    // ltau->Step3_WJ_FakeRate_GetHistos();

    // ltau->Step3_WJ_Closure_DeclareHistos();
    // ltau->Step3_WJ_Closure_FillHistos();
    // ltau->Step3_WJ_Closure_Compute();
    // ltau->Step3_WJ_Closure_SaveHistos();
    // ltau->Step3_WJ_Closure_GetHistos();

    // ltau->Step3_WJ_vsSR_DeclareHistos();
    // ltau->Step3_WJ_vsSR_FillHistos();
    // ltau->Step3_WJ_vsSR_Compute();
    // ltau->Step3_WJ_vsSR_SaveHistos();
    // ltau->Step3_WJ_vsSR_GetHistos();

    cout<<"---------------------------------------------------"<<endl;
    cout<<"----------Step 4: Compute the fractions------------"<<endl;
    cout<<"---------------------------------------------------"<<endl;
    ltau->Step4_Fraction_DeclareHistos();
    ltau->Step4_Fraction_FillHistos();
    ltau->Step4_Fraction_Compute();
    ltau->Step4_Fraction_SaveHistos();
    ltau->Step4_Fraction_GetHistos();
}
