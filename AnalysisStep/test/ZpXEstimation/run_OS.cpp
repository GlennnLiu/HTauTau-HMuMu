// C++
#include <iostream>
#include <fstream>
#include <string>

// ROOT
#include "TApplication.h"
#include <TROOT.h>
#include "TFile.h"
#include "TString.h"
#include "TStyle.h"

// My own files
#include <HTauTauHMuMu/AnalysisStep/test/ZpXEstimation/include/OSmethod.h>
#include <HTauTauHMuMu/AnalysisStep/test/ZpXEstimation/src/setTDRStyle.cpp>


using namespace std;

int main( int argc, char *argv[] )
{
   setTDRStyle();
   
   TString path = "/eos/home-g/geliu/LepUni/BigTrees/Simp/";
   TString file_name = "/HTauTauHMuMu.root";
   
   TString Data    = path + "AllData" + file_name;
   TString WZ      = path + "WZTo3LNu/HTauTauHMuMu.root";//       + file_name;
   TString ZZ      = path + "ZZTo4l/HTauTauHMuMu.root";//      + file_name;
   TString ttbar   = path + "TTTo2L2Nu"      + file_name;
   TString DY      = path + "DYJetsToLL_M50/HTauTauHMuMu.root";// + file_name;
	
   bool SubtractWZ = true;
   bool Remove_NegBins_FR = true;
   bool Remove_NegBins_ZX = true;
	
   
   float pT_bins[] = {5, 7, 10, 20, 25, 30, 40, 50, 80};
   float MET_bins[] = {0, 10, 18, 25, 32, 40, 50, 60, 80, 100};
   
   
   OSmethod *os = new OSmethod();
   
   os->SetLumi(35.92); //2016
   //os->SetLumi(41.53);   //2017
   //os->SetLumi(59.74); //2018
   //cout<<"begin"<<endl;
   ///////////////////////////////////
   // Fill control histos           //
   ///////////////////////////////////
   os->FillDataMCPlots(Data);
   os->FillDataMCPlots(WZ);
   os->FillDataMCPlots(ZZ);
   os->FillDataMCPlots(ttbar);
   os->FillDataMCPlots(DY);
   os->SaveDataMCHistos("DataMC_OS.root");

   ///////////////////////////////////
   // Fill passing/failling histos  //
   ///////////////////////////////////
   os->FillFRHistos(Data);
   os->FillFRHistos(WZ);
   os->FillFRHistos(ZZ);
   os->SaveFRHistos("Histos_OS.root", SubtractWZ, Remove_NegBins_FR);
   ///////////////////////////////////
   // Calculate fake rates          //
   ///////////////////////////////////
   os->GetFRHistos("Histos_OS.root");
   os->Set_pT_binning(9, pT_bins);
   os->Set_MET_binning(10, MET_bins);
   os->ProduceFakeRates("FakeRates_OS.root");

   ///////////////////////////////////
   // Fill ZX contributions histos  //
   ///////////////////////////////////
   os->MakeHistogramsZX(Data, "FakeRates_OS.root");
   os->MakeZXMCContribution(ZZ, "FakeRates_OS.root");
   os->SaveZXHistos("ZXHistos_OS.root", Remove_NegBins_ZX);

   ///////////////////////////////////
   // Plot control plots            //
   ///////////////////////////////////
   os->GetZXHistos("ZXHistos_OS.root");
   os->PrintZXYields();
   os->GetDataMCHistos("DataMC_OS.root");
   os->PlotDataMC("M4l", "Plots");
   os->PlotDataMC_2P2F( "M4l", "Plots" );
   os->PlotDataMC_3P1F( "M4l", "Plots" );

   ///////////////////////////////////
   // Plot Z+X plots                //
   ///////////////////////////////////
   os->GetZXHistos("ZXHistos_OS.root");
   os->PlotZXContributions("Plots");
   os->FitZX("Plots");
	
   delete os;
}
