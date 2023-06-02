#include <HTauTauHMuMu/AnalysisStep/interface/LeptonSFHelper.h>

using namespace std;
//using namespace edm;

LeptonSFHelper::LeptonSFHelper(bool preVFP)
{
   // 2016 preVFP Electrons
   if(preVFP)
   {  
      TString fipEle_2016 = Form("$CMSSW_BASE/src/HTauTauHMuMu/AnalysisStep/data/LeptonEffScaleFactors/egammaEffi.txt_Ele_wp90noiso_preVFP_EGM2D.root");
      root_file = TFile::Open(fipEle_2016.Data(),"READ");
      h_Ele_2016 = (TH2F*) root_file->Get("EGamma_SF2D")->Clone();
   }
   // 2016 postVFP Electrons
   else
   {  
      TString fipEle_2016 = Form("$CMSSW_BASE/src/HTauTauHMuMu/AnalysisStep/data/LeptonEffScaleFactors/egammaEffi.txt_Ele_wp90noiso_postVFP_EGM2D.root");
      root_file = TFile::Open(fipEle_2016.Data(),"READ");
      h_Ele_2016 = (TH2F*) root_file->Get("EGamma_SF2D")->Clone();
   }

   // 2017 Electrons
   TString fipEle_2017 = Form("$CMSSW_BASE/src/HTauTauHMuMu/AnalysisStep/data/LeptonEffScaleFactors/egammaEffi.txt_EGM2D_MVA90noIso_UL17.root");
   root_file = TFile::Open(fipEle_2017.Data(),"READ");
   h_Ele_2017 = (TH2F*) root_file->Get("EGamma_SF2D")->Clone();

   // 2018 Electrons
   TString fipEle_2018 = Form("$CMSSW_BASE/src/HTauTauHMuMu/AnalysisStep/data/LeptonEffScaleFactors/egammaEffi.txt_Ele_wp90noiso_EGM2D.root");
   root_file = TFile::Open(fipEle_2018.Data(),"READ");
   h_Ele_2018 = (TH2F*) root_file->Get("EGamma_SF2D")->Clone();

   // 2016 preVFP Muons
   if (preVFP) {
      TString fipMu_RECO_2016 = Form("$CMSSW_BASE/src/HTauTauHMuMu/AnalysisStep/data/LeptonEffScaleFactors/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_RECO.root");
      root_file = TFile::Open(fipMu_RECO_2016.Data(),"READ");
      h_Mu_RECO_syst_2016 = (TH2D*)root_file->Get("NUM_TrackerMuons_DEN_genTracks_abseta_pt_syst");
      h_Mu_RECO_stat_2016 = (TH2D*)root_file->Get("NUM_TrackerMuons_DEN_genTracks_abseta_pt_stat");

      TString fipMu_ID_2016 = Form("$CMSSW_BASE/src/HTauTauHMuMu/AnalysisStep/data/LeptonEffScaleFactors/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ID.root");
      root_file = TFile::Open(fipMu_ID_2016.Data(),"READ");
      h_Mu_ID_syst_2016 = (TH2D*)root_file->Get("NUM_MediumID_DEN_TrackerMuons_abseta_pt_syst");
      h_Mu_ID_stat_2016 = (TH2D*)root_file->Get("NUM_MediumID_DEN_TrackerMuons_abseta_pt_stat");

      TString fipMu_ISO_2016 = Form("$CMSSW_BASE/src/HTauTauHMuMu/AnalysisStep/data/LeptonEffScaleFactors/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ISO.root");
      root_file = TFile::Open(fipMu_ISO_2016.Data(),"READ");
      h_Mu_ISO_syst_2016 = (TH2D*)root_file->Get("NUM_TightRelIso_DEN_MediumID_abseta_pt_syst");
      h_Mu_ISO_stat_2016 = (TH2D*)root_file->Get("NUM_TightRelIso_DEN_MediumID_abseta_pt_stat");
   }
   // 2016 postVFP Muons
   else {
      TString fipMu_RECO_2016 = Form("$CMSSW_BASE/src/HTauTauHMuMu/AnalysisStep/data/LeptonEffScaleFactors/Efficiencies_muon_generalTracks_Z_Run2016_UL_RECO.root");
      root_file = TFile::Open(fipMu_RECO_2016.Data(),"READ");
      h_Mu_RECO_syst_2016 = (TH2D*)root_file->Get("NUM_TrackerMuons_DEN_genTracks_abseta_pt_syst");
      h_Mu_RECO_stat_2016 = (TH2D*)root_file->Get("NUM_TrackerMuons_DEN_genTracks_abseta_pt_stat");

      TString fipMu_ID_2016 = Form("$CMSSW_BASE/src/HTauTauHMuMu/AnalysisStep/data/LeptonEffScaleFactors/Efficiencies_muon_generalTracks_Z_Run2016_UL_ID.root");
      root_file = TFile::Open(fipMu_ID_2016.Data(),"READ");
      h_Mu_ID_syst_2016 = (TH2D*)root_file->Get("NUM_MediumID_DEN_TrackerMuons_abseta_pt_syst");
      h_Mu_ID_stat_2016 = (TH2D*)root_file->Get("NUM_MediumID_DEN_TrackerMuons_abseta_pt_stat");

      TString fipMu_ISO_2016 = Form("$CMSSW_BASE/src/HTauTauHMuMu/AnalysisStep/data/LeptonEffScaleFactors/Efficiencies_muon_generalTracks_Z_Run2016_UL_ISO.root");
      root_file = TFile::Open(fipMu_ISO_2016.Data(),"READ");
      h_Mu_ISO_syst_2016 = (TH2D*)root_file->Get("NUM_TightRelIso_DEN_MediumID_abseta_pt_syst");
      h_Mu_ISO_stat_2016 = (TH2D*)root_file->Get("NUM_TightRelIso_DEN_MediumID_abseta_pt_stat");
   }
   // 2017 Muons
   TString fipMu_RECO_2017 = Form("$CMSSW_BASE/src/HTauTauHMuMu/AnalysisStep/data/LeptonEffScaleFactors/Efficiencies_muon_generalTracks_Z_Run2017_UL_RECO.root");
   root_file = TFile::Open(fipMu_RECO_2017.Data(),"READ");
   h_Mu_RECO_syst_2017 = (TH2D*)root_file->Get("NUM_TrackerMuons_DEN_genTracks_abseta_pt_syst");
   h_Mu_RECO_stat_2017 = (TH2D*)root_file->Get("NUM_TrackerMuons_DEN_genTracks_abseta_pt_stat");

   TString fipMu_ID_2017 = Form("$CMSSW_BASE/src/HTauTauHMuMu/AnalysisStep/data/LeptonEffScaleFactors/Efficiencies_muon_generalTracks_Z_Run2017_UL_ID.root");
   root_file = TFile::Open(fipMu_ID_2017.Data(),"READ");
   h_Mu_ID_syst_2017 = (TH2D*)root_file->Get("NUM_MediumID_DEN_TrackerMuons_abseta_pt_syst");
   h_Mu_ID_stat_2017 = (TH2D*)root_file->Get("NUM_MediumID_DEN_TrackerMuons_abseta_pt_stat");

   TString fipMu_ISO_2017 = Form("$CMSSW_BASE/src/HTauTauHMuMu/AnalysisStep/data/LeptonEffScaleFactors/Efficiencies_muon_generalTracks_Z_Run2017_UL_ISO.root");
   root_file = TFile::Open(fipMu_ISO_2017.Data(),"READ");
   h_Mu_ISO_syst_2017 = (TH2D*)root_file->Get("NUM_TightRelIso_DEN_MediumID_abseta_pt_syst");
   h_Mu_ISO_stat_2017 = (TH2D*)root_file->Get("NUM_TightRelIso_DEN_MediumID_abseta_pt_stat");

   // 2018 Muons
   TString fipMu_RECO_2018 = Form("$CMSSW_BASE/src/HTauTauHMuMu/AnalysisStep/data/LeptonEffScaleFactors/Efficiencies_muon_generalTracks_Z_Run2018_UL_RECO.root");
   root_file = TFile::Open(fipMu_RECO_2018.Data(),"READ");
   h_Mu_RECO_syst_2018 = (TH2D*)root_file->Get("NUM_TrackerMuons_DEN_genTracks_abseta_pt_syst");
   h_Mu_RECO_stat_2018 = (TH2D*)root_file->Get("NUM_TrackerMuons_DEN_genTracks_abseta_pt_stat");

   TString fipMu_ID_2018 = Form("$CMSSW_BASE/src/HTauTauHMuMu/AnalysisStep/data/LeptonEffScaleFactors/Efficiencies_muon_generalTracks_Z_Run2018_UL_ID.root");
   root_file = TFile::Open(fipMu_ID_2018.Data(),"READ");
   h_Mu_ID_syst_2018 = (TH2D*)root_file->Get("NUM_MediumID_DEN_TrackerMuons_abseta_pt_syst");
   h_Mu_ID_stat_2018 = (TH2D*)root_file->Get("NUM_MediumID_DEN_TrackerMuons_abseta_pt_stat");

   TString fipMu_ISO_2018 = Form("$CMSSW_BASE/src/HTauTauHMuMu/AnalysisStep/data/LeptonEffScaleFactors/Efficiencies_muon_generalTracks_Z_Run2018_UL_ISO.root");
   root_file = TFile::Open(fipMu_ISO_2018.Data(),"READ");
   h_Mu_ISO_syst_2018 = (TH2D*)root_file->Get("NUM_TightRelIso_DEN_MediumID_abseta_pt_syst");
   h_Mu_ISO_stat_2018 = (TH2D*)root_file->Get("NUM_TightRelIso_DEN_MediumID_abseta_pt_stat");

   cout << "[LeptonSFHelper] SF maps opened from root files." << endl;
}

LeptonSFHelper::~LeptonSFHelper()
{
}

float LeptonSFHelper::getSF(int year, int flav, float pt, float eta, float SCeta, const std::string& unc, const std::string& level) const
{
   float SF = 1.0;
   float RECOSF = 1.0;
   float IDSF = 1.0;
   float ISOSF = 1.0;

   // Electron SFs
   if(abs(flav) == 11) {
      if(year == 2016)
      {
         if (unc=="unc") SF = h_Ele_2016->GetBinError(h_Ele_2016->GetXaxis()->FindBin(SCeta),h_Ele_2016->GetYaxis()->FindBin(std::min(pt,499.f)));
         else SF = h_Ele_2016->GetBinContent(h_Ele_2016->GetXaxis()->FindBin(SCeta),h_Ele_2016->GetYaxis()->FindBin(std::min(pt,499.f)));
      }
      else if(year == 2017)
      {
         if (unc=="unc") SF = h_Ele_2017->GetBinContent(h_Ele_2017->GetXaxis()->FindBin(SCeta),h_Ele_2017->GetYaxis()->FindBin(std::min(pt,499.f)));
         else SF = h_Ele_2017->GetBinError(h_Ele_2017->GetXaxis()->FindBin(SCeta),h_Ele_2017->GetYaxis()->FindBin(std::min(pt,499.f)));
      }
      else if(year == 2018)
      {
         if (unc=="unc") SF = h_Ele_2018->GetBinContent(h_Ele_2018->GetXaxis()->FindBin(SCeta),h_Ele_2018->GetYaxis()->FindBin(std::min(pt,499.f)));
         else SF = h_Ele_2018->GetBinError(h_Ele_2018->GetXaxis()->FindBin(SCeta),h_Ele_2018->GetYaxis()->FindBin(std::min(pt,499.f)));
      }
      else {
         edm::LogError("LeptonSFHelper::") << "Ele SFs for " << year << " is not supported!";
         abort();
      }
   }

   //Muon SF
   else if(abs(flav) == 13 )
   {
      if(year == 2016)
      {
         bool RECO_unc=true, ID_unc=true, ISO_unc=true;
         if (unc=="unc") {
            if (level=="RECO_syst") {
               RECOSF = h_Mu_RECO_syst_2016->GetBinError(h_Mu_RECO_syst_2016->GetXaxis()->FindBin(fabs(eta)),h_Mu_RECO_syst_2016->GetYaxis()->FindBin(50.));
            }
            else if (level=="RECO_stat") {
               RECOSF = h_Mu_RECO_stat_2016->GetBinError(h_Mu_RECO_stat_2016->GetXaxis()->FindBin(fabs(eta)),h_Mu_RECO_stat_2016->GetYaxis()->FindBin(50.));
            }
            else {
               RECO_unc=false;
               RECOSF = h_Mu_RECO_syst_2016->GetBinContent(h_Mu_RECO_syst_2016->GetXaxis()->FindBin(fabs(eta)),h_Mu_RECO_syst_2016->GetYaxis()->FindBin(50.));
            }

            if (level=="ID_syst") {
               IDSF = h_Mu_ID_syst_2016->GetBinError(h_Mu_ID_syst_2016->GetXaxis()->FindBin(fabs(eta)),h_Mu_ID_syst_2016->GetYaxis()->FindBin(std::min(pt,119.f)));
            }
            else if (level=="ID_stat") {
               IDSF = h_Mu_ID_stat_2016->GetBinError(h_Mu_ID_stat_2016->GetXaxis()->FindBin(fabs(eta)),h_Mu_ID_stat_2016->GetYaxis()->FindBin(std::min(pt,119.f)));
            }
            else {
               ID_unc=false;
               IDSF = h_Mu_ID_syst_2016->GetBinContent(h_Mu_ID_syst_2016->GetXaxis()->FindBin(fabs(eta)),h_Mu_ID_syst_2016->GetYaxis()->FindBin(std::min(pt,119.f)));
            }

            if (level=="ISO_syst") {
               ISOSF = h_Mu_ISO_syst_2016->GetBinError(h_Mu_ISO_syst_2016->GetXaxis()->FindBin(fabs(eta)),h_Mu_ISO_syst_2016->GetYaxis()->FindBin(std::min(pt,119.f)));
            }
            else if (level=="ISO_stat") {
               ISOSF = h_Mu_ISO_stat_2016->GetBinError(h_Mu_ISO_stat_2016->GetXaxis()->FindBin(fabs(eta)),h_Mu_ISO_stat_2016->GetYaxis()->FindBin(std::min(pt,119.f)));
            }
            else {
               ISO_unc=false;
               ISOSF = h_Mu_ISO_syst_2016->GetBinContent(h_Mu_ISO_syst_2016->GetXaxis()->FindBin(fabs(eta)),h_Mu_ISO_syst_2016->GetYaxis()->FindBin(std::min(pt,119.f)));
            }
            if (RECO_unc==false && ID_unc==false && ISO_unc==false) {
               cout<<"[ERROR] Require uncertainties but no sources are required!!!";
               abort();
            }
         }
         else {
            RECOSF = h_Mu_RECO_syst_2016->GetBinContent(h_Mu_RECO_syst_2016->GetXaxis()->FindBin(fabs(eta)),h_Mu_RECO_syst_2016->GetYaxis()->FindBin(50.));
            IDSF = h_Mu_ID_syst_2016->GetBinContent(h_Mu_ID_syst_2016->GetXaxis()->FindBin(fabs(eta)),h_Mu_ID_syst_2016->GetYaxis()->FindBin(std::min(pt,119.f)));
            ISOSF = h_Mu_ISO_syst_2016->GetBinContent(h_Mu_ISO_syst_2016->GetXaxis()->FindBin(fabs(eta)),h_Mu_ISO_syst_2016->GetYaxis()->FindBin(std::min(pt,119.f)));
         }
      }
      else if(year == 2017)
      {
         bool RECO_unc=true, ID_unc=true, ISO_unc=true;
         if (unc=="error") {
            if (level=="RECO_syst") {
               RECOSF = h_Mu_RECO_syst_2017->GetBinError(h_Mu_RECO_syst_2017->GetXaxis()->FindBin(fabs(eta)),h_Mu_RECO_syst_2017->GetYaxis()->FindBin(50.));
            }
            else if (level=="RECO_stat") {
               RECOSF = h_Mu_RECO_stat_2017->GetBinError(h_Mu_RECO_stat_2017->GetXaxis()->FindBin(fabs(eta)),h_Mu_RECO_stat_2017->GetYaxis()->FindBin(50.));
            }
            else {
               RECO_unc=false;
               RECOSF = h_Mu_RECO_syst_2017->GetBinContent(h_Mu_RECO_syst_2017->GetXaxis()->FindBin(fabs(eta)),h_Mu_RECO_syst_2017->GetYaxis()->FindBin(50.));
            }

            if (level=="ID_syst") {
               IDSF = h_Mu_ID_syst_2017->GetBinError(h_Mu_ID_syst_2017->GetXaxis()->FindBin(fabs(eta)),h_Mu_ID_syst_2017->GetYaxis()->FindBin(std::min(pt,119.f)));
            }
            else if (level=="ID_stat") {
               IDSF = h_Mu_ID_stat_2017->GetBinError(h_Mu_ID_stat_2017->GetXaxis()->FindBin(fabs(eta)),h_Mu_ID_stat_2017->GetYaxis()->FindBin(std::min(pt,119.f)));
            }
            else {
               ID_unc=false;
               IDSF = h_Mu_ID_syst_2017->GetBinContent(h_Mu_ID_syst_2017->GetXaxis()->FindBin(fabs(eta)),h_Mu_ID_syst_2017->GetYaxis()->FindBin(std::min(pt,119.f)));
            }

            if (level=="ISO_syst") {
               ISOSF = h_Mu_RECO_syst_2017->GetBinError(h_Mu_RECO_syst_2017->GetXaxis()->FindBin(fabs(eta)),h_Mu_RECO_syst_2017->GetYaxis()->FindBin(std::min(pt,119.f)));
            }
            else if (level=="ISO_stat") {
               ISOSF = h_Mu_RECO_stat_2017->GetBinError(h_Mu_RECO_stat_2017->GetXaxis()->FindBin(fabs(eta)),h_Mu_RECO_stat_2017->GetYaxis()->FindBin(std::min(pt,119.f)));
            }
            else {
               ISO_unc=false;
               ISOSF = h_Mu_RECO_syst_2017->GetBinContent(h_Mu_RECO_syst_2017->GetXaxis()->FindBin(fabs(eta)),h_Mu_RECO_syst_2017->GetYaxis()->FindBin(std::min(pt,119.f)));
            }
            if (RECO_unc==false && ID_unc==false && ISO_unc==false) {
               cout<<"[ERROR] Require uncertainties but no sources are required!!!";
               abort();
            }
         }
         else {
            RECOSF = h_Mu_RECO_syst_2017->GetBinContent(h_Mu_RECO_syst_2017->GetXaxis()->FindBin(fabs(eta)),h_Mu_RECO_syst_2017->GetYaxis()->FindBin(50.));
            IDSF = h_Mu_ID_syst_2017->GetBinContent(h_Mu_ID_syst_2017->GetXaxis()->FindBin(fabs(eta)),h_Mu_ID_syst_2017->GetYaxis()->FindBin(std::min(pt,119.f)));
            ISOSF = h_Mu_RECO_syst_2017->GetBinContent(h_Mu_RECO_syst_2017->GetXaxis()->FindBin(fabs(eta)),h_Mu_RECO_syst_2017->GetYaxis()->FindBin(std::min(pt,119.f)));
         }
      }
      else if(year == 2018)
      {
         bool RECO_unc=true, ID_unc=true, ISO_unc=true;
         if (unc=="error") {
            if (level=="RECO_syst") {
               RECOSF = h_Mu_RECO_syst_2018->GetBinError(h_Mu_RECO_syst_2018->GetXaxis()->FindBin(fabs(eta)),h_Mu_RECO_syst_2018->GetYaxis()->FindBin(50.));
            }
            else if (level=="RECO_stat") {
               RECOSF = h_Mu_RECO_stat_2018->GetBinError(h_Mu_RECO_stat_2018->GetXaxis()->FindBin(fabs(eta)),h_Mu_RECO_stat_2018->GetYaxis()->FindBin(50.));
            }
            else {
               RECO_unc=false;
               RECOSF = h_Mu_RECO_syst_2018->GetBinContent(h_Mu_RECO_syst_2018->GetXaxis()->FindBin(fabs(eta)),h_Mu_RECO_syst_2018->GetYaxis()->FindBin(50.));
            }

            if (level=="ID_syst") {
               IDSF = h_Mu_ID_syst_2018->GetBinError(h_Mu_ID_syst_2018->GetXaxis()->FindBin(fabs(eta)),h_Mu_ID_syst_2018->GetYaxis()->FindBin(std::min(pt,119.f)));
            }
            else if (level=="ID_stat") {
               IDSF = h_Mu_ID_stat_2018->GetBinError(h_Mu_ID_stat_2018->GetXaxis()->FindBin(fabs(eta)),h_Mu_ID_stat_2018->GetYaxis()->FindBin(std::min(pt,119.f)));
            }
            else {
               ID_unc=false;
               IDSF = h_Mu_ID_syst_2018->GetBinContent(h_Mu_ID_syst_2018->GetXaxis()->FindBin(fabs(eta)),h_Mu_ID_syst_2018->GetYaxis()->FindBin(std::min(pt,119.f)));
            }

            if (level=="ISO_syst") {
               ISOSF = h_Mu_RECO_syst_2018->GetBinError(h_Mu_RECO_syst_2018->GetXaxis()->FindBin(fabs(eta)),h_Mu_RECO_syst_2018->GetYaxis()->FindBin(std::min(pt,119.f)));
            }
            else if (level=="ISO_stat") {
               ISOSF = h_Mu_RECO_stat_2018->GetBinError(h_Mu_RECO_stat_2018->GetXaxis()->FindBin(fabs(eta)),h_Mu_RECO_stat_2018->GetYaxis()->FindBin(std::min(pt,119.f)));
            }
            else {
               ISO_unc=false;
               ISOSF = h_Mu_RECO_syst_2018->GetBinContent(h_Mu_RECO_syst_2018->GetXaxis()->FindBin(fabs(eta)),h_Mu_RECO_syst_2018->GetYaxis()->FindBin(std::min(pt,119.f)));
            }
            if (RECO_unc==false && ID_unc==false && ISO_unc==false) {
               cout<<"[ERROR] Require uncertainties but no sources are required!!!";
               abort();
            }
         }
         else {
            RECOSF = h_Mu_RECO_syst_2018->GetBinContent(h_Mu_RECO_syst_2018->GetXaxis()->FindBin(fabs(eta)),h_Mu_RECO_syst_2018->GetYaxis()->FindBin(50.));
            IDSF = h_Mu_ID_syst_2018->GetBinContent(h_Mu_ID_syst_2018->GetXaxis()->FindBin(fabs(eta)),h_Mu_ID_syst_2018->GetYaxis()->FindBin(std::min(pt,119.f)));
            ISOSF = h_Mu_RECO_syst_2018->GetBinContent(h_Mu_RECO_syst_2018->GetXaxis()->FindBin(fabs(eta)),h_Mu_RECO_syst_2018->GetYaxis()->FindBin(std::min(pt,119.f)));
         }  
      }
      else {
         edm::LogError("LeptonSFHelper::") << "Ele SFs for " << year << " is not supported!";
         abort();
      }

      SF = RECOSF * IDSF * ISOSF;
   }

    return SF;
}