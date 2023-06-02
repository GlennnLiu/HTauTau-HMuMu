#ifndef LEPTONSFHELPER_H
#define LEPTONSFHELPER_H

#include <string>
#include <iostream>
#include <vector>
#include <utility>

#include <cmath>
#include "TString.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH2D.h"

#include <FWCore/ParameterSet/interface/FileInPath.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h>


enum SFsyst {central = 0, up = 1, down = 2};

class LeptonSFHelper
{

 public:

  LeptonSFHelper(bool preVFP);
  ~LeptonSFHelper();
  
  float getSF (int year, int flav, float pt, float eta, float SCeta, const std::string& unc="", const std::string& level="") const;
   
 private:
   TFile *root_file;
   
   // Electron SF map histograms
   TH2F *h_Ele_2016, *h_Ele_2017, *h_Ele_2018;
   
   // Muons SF map histograms
  //  TH2D *h_Mu_TRG_2016, *h_Mu_TRG_2017, *h_Mu_TRG_2018;
   TH2D *h_Mu_RECO_syst_2016, *h_Mu_RECO_syst_2017, *h_Mu_RECO_syst_2018;
   TH2D *h_Mu_ID_syst_2016, *h_Mu_ID_syst_2017, *h_Mu_ID_syst_2018;
   TH2D *h_Mu_ISO_syst_2016, *h_Mu_ISO_syst_2017, *h_Mu_ISO_syst_2018;
   TH2D *h_Mu_RECO_stat_2016, *h_Mu_RECO_stat_2017, *h_Mu_RECO_stat_2018;
   TH2D *h_Mu_ID_stat_2016, *h_Mu_ID_stat_2017, *h_Mu_ID_stat_2018;
   TH2D *h_Mu_ISO_stat_2016, *h_Mu_ISO_stat_2017, *h_Mu_ISO_stat_2018;

};

#endif
