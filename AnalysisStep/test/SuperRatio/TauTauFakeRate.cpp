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

#include <HTauTauHMuMu/AnalysisStep/test/SuperRatio/include/Settings.h>
#include <HTauTauHMuMu/AnalysisStep/test/SuperRatio/include/TauTau.h>

using namespace std;

int main( int argc, char *argv[] )

{
    string path = "/eos/home-g/geliu/LepUni/BigTrees/UL/2016/";
    string file_name = "/HTauTauHMuMu_new.root";
    
    string data[]={"SingleEle2016","SingleMuon2016","DoubleMu2016","MuonEG2016","Tau2016"};
    string bkg[]={"DYJetsToLL","EWKZToLL","TTTo2L2Nu","TTW","TTWW","TTZToLLNuNu_M10","TTZToLL_M1to10","TTZZ","WJetsToLNu","WWTo2L2Nu","WWW","WWZ","WZTo2L2Q","WZTo3LNu","WZZ","ZZTo2L2Nu","ZZTo2L2Q","ZZTo4l","ZZZ","ggTo2e2mu","ggTo2e2nu","ggTo2e2tau","ggTo2mu2nu","ggTo2mu2tau","ggTo4e","ggTo4mu","ggTo4tau","tW_antitop","tW_top","tZq","TTToHadronic1","TTToHadronic2","TTToSemiLeptonic1","TTToSemiLeptonic2","TTToSemiLeptonic3","ST_antitop","ST_s","ST_top"};
    
    string savepath="/eos/home-g/geliu/LepUni/SuperRatio/";

    float taupt_bins[]={40,42,44,46,48,50,52,54,56,58,60,65,70,75,80,90,100,110,120,150};
    float lpt_bins[]={40,42,44,46,48,50,52,54,56,58,60,65,70,75,80,90,100,110,120,150};
    float Mvis_fine_bins[]={20,40,60,70,80,90,100,110,120,130,140,150,160,180,200,230,260,300};
    float Mvis_bins[]={20,60,90,120,150,200,300};
    
    TauTau *tautau = new TauTau();
    
    tautau->SetPaths(path,file_name,savepath);

    tautau->SetFileList(data,size(data),bkg,size(bkg));

    tautau->Set_taupt_bin(size(taupt_bins),taupt_bins);
    tautau->Set_lpt_bin(size(lpt_bins),lpt_bins);
    tautau->Set_Mvis_bin(size(Mvis_bins),Mvis_bins,false);
    tautau->Set_Mvis_bin(size(Mvis_fine_bins),Mvis_fine_bins,true);

    cout<<"---------------------------------------------------"<<endl;
    cout<<"----------Step 1: compute fake rates---------------"<<endl;
    cout<<"---------------------------------------------------"<<endl;
    tautau->Step1_FakeRate_DeclareHistos();
    tautau->Step1_FakeRate_FillHistos();
    tautau->Step1_FakeRate_Compute();
    tautau->Step1_FakeRate_SaveHistos();
    tautau->Step1_FakeRate_GetHistos();

    cout<<"---------------------------------------------------"<<endl;
    cout<<"----------Step 2: compute closure factors----------"<<endl;
    cout<<"---------------------------------------------------"<<endl;
    tautau->Step2_Closure_DeclareHistos();
    tautau->Step2_Closure_FillHistos();
    tautau->Step2_Closure_Compute();
    tautau->Step2_Closure_SaveHistos();
    tautau->Step2_Closure_GetHistos();

    cout<<"---------------------------------------------------"<<endl;
    cout<<"----------Step 3.1: Diff vs SR for QCD-------------"<<endl;
    cout<<"---------------------------------------------------"<<endl;
    tautau->Step3_QCD_FakeRate_DeclareHistos();
    tautau->Step3_QCD_FakeRate_FillHistos();
    tautau->Step3_QCD_FakeRate_Compute();
    tautau->Step3_QCD_FakeRate_SaveHistos();
    tautau->Step3_QCD_FakeRate_GetHistos();

    tautau->Step3_QCD_Closure_DeclareHistos();
    tautau->Step3_QCD_Closure_FillHistos();
    tautau->Step3_QCD_Closure_Compute();
    tautau->Step3_QCD_Closure_SaveHistos();
    tautau->Step3_QCD_Closure_GetHistos();

    tautau->Step3_QCD_vsSR_DeclareHistos();
    tautau->Step3_QCD_vsSR_FillHistos();
    tautau->Step3_QCD_vsSR_Compute();
    tautau->Step3_QCD_vsSR_SaveHistos();
    tautau->Step3_QCD_vsSR_GetHistos();
}
