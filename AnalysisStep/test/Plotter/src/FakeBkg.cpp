#include <HTauTauHMuMu/AnalysisStep/test/Plotter/include/FakeBkg.h>
#include <iostream> // std::cerr, std::endl
#include <iomanip> 
#include <assert.h> // assert

TFile* ensureTFile(const TString filename, bool verbose=false){
  if(verbose)
    std::cout << "Opening " << filename << std::endl;
  TFile* file = new TFile(filename);
  if(!file or file->IsZombie()) {
    std::cerr << std::endl << "ERROR! Failed to open input file = '" << filename << "'!" << std::endl;
    assert(0);
  }
  return file;
}

TH1F* extractTH1(const TFile* file, const std::string& histname){
  TH1F* hist = dynamic_cast<TH1F*>((const_cast<TFile*>(file))->Get(histname.data()));
  if(!hist){
    std::cerr << std::endl << "ERROR! Failed to load histogram = '" << histname << "' from input file!" << std::endl;
    assert(0);
  }
  return hist;
}

TH2F* extractTH2(const TFile* file, const std::string& histname){
  TH2F* hist = dynamic_cast<TH2F*>((const_cast<TFile*>(file))->Get(histname.data()));
  if(!hist){
    std::cerr << std::endl << "ERROR! Failed to load histogram = '" << histname << "' from input file!" << std::endl;
    assert(0);
  }
  return hist;
}

const TF1* extractTF1(const TFile* file, const std::string& funcname){
  const TF1* function = dynamic_cast<TF1*>((const_cast<TFile*>(file))->Get(funcname.data()));
  if(!function){
    std::cerr << std::endl << "ERROR! Failed to load function = '" << funcname << "' from input file!" << std::endl;
    assert(0);
  }
  return function;
}

FakeBkg::FakeBkg(const std::string& year, const std::string& path) {

    bool verbose = false;

    _s_final_state.push_back("etau");
    _s_final_state.push_back("mutau");
    _s_final_state.push_back("tautau");
    _s_final_state.push_back("emu");

    _s_categories.push_back("GGH");
    _s_categories.push_back("VBF_ptHl200");
    _s_categories.push_back("VBF_ptHg200");
    _s_categories.push_back("Boost_1j");
    _s_categories.push_back("Boost_2j");
    _s_categories.push_back("All");

    _s_fake_bkg.push_back("QCD");
    _s_fake_bkg.push_back("Wjet");
    _s_fake_bkg.push_back("TT");

    f_LTau_FR=ensureTFile(path+"LTau.Step1_FakeRate.root");
    for (int i_bkg=0;i_bkg<3;i_bkg++) {
        for (int i_fs=0;i_fs<2;i_fs++) {
            for (int i_njet=0;i_njet<3;i_njet++) {
                _histo_name="f_FR_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs)+"_Njet"+std::to_string(i_njet);
                LTau_f_FR[i_bkg][i_fs][i_njet]=extractTF1(f_LTau_FR,_histo_name);
            }
        }
    }

    f_LTau_Closure_FR=ensureTFile(path+"LTau.Step2_Closure.root");
    for (int i_bkg=0;i_bkg<3;i_bkg++) {
        for (int i_fs=0;i_fs<2;i_fs++) {
            for (int i_taupt=0;i_taupt<3;i_taupt++) {
                _histo_name="hClosure_FR_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs)+"_Ntaupt"+std::to_string(i_taupt);
                LTau_hClosure_FR[i_bkg][i_fs][i_taupt]=extractTH1(f_LTau_Closure_FR,_histo_name);
                _histo_name="fClosure_FR_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs)+"_Ntaupt"+std::to_string(i_taupt);
                LTau_fClosure_FR[i_bkg][i_fs][i_taupt]=extractTF1(f_LTau_Closure_FR,_histo_name);
            }
        }
    }

    f_LTau_QCDvsSR_FR=ensureTFile(path+"LTau.Step3_QCD_vsSR.root");
    for (int i_fs=0;i_fs<2;i_fs++) {
        _histo_name="fQCDvsSR_FR_"+_s_final_state.at(i_fs);
        LTau_fQCDvsSR_FR[i_fs]=extractTF1(f_LTau_QCDvsSR_FR,_histo_name);
    }

    f_LTau_WJvsSR=ensureTFile(path+"LTau.Step3_WJ_vsSR.root");
    for (int i_fs=0;i_fs<2;i_fs++) {
        _histo_name="fWJvsSR_FR_"+_s_final_state.at(i_fs);
        LTau_fWJvsSR_FR[i_fs]=extractTF1(f_LTau_WJvsSR,_histo_name);
    }

    f_LTau_Frac=ensureTFile(path+"LTau.Step4_Fraction.root");
    for (int i_bkg=0;i_bkg<2+1;i_bkg++) {
        for (int i_fs=0;i_fs<2;i_fs++) {
            for (int i_cat=0;i_cat<5;i_cat++) {
                _histo_name="hFrac_FR_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs)+"_"+_s_categories.at(i_cat);
                LTau_hFrac_FR[i_bkg][i_fs][i_cat]=extractTH1(f_LTau_Frac,_histo_name);
            }
        }
    }

    f_TauTau_FR=ensureTFile(path+"TauTau.Step1_FakeRate.root");
    for (int i_njet=0;i_njet<3;i_njet++) {
        _histo_name="f_FR_"+_s_fake_bkg.at(0)+"_"+_s_final_state.at(2)+"_Njet"+std::to_string(i_njet);
        TauTau_f_FR[i_njet]=extractTF1(f_TauTau_FR,_histo_name);
    }

    f_TauTau_Closure_FR=ensureTFile(path+"TauTau.Step2_Closure.root");
    _histo_name="fClosure_FR_"+_s_fake_bkg.at(0)+"_"+_s_final_state.at(2);
    TauTau_fClosure_FR=extractTF1(f_TauTau_Closure_FR,_histo_name);

    f_TauTau_QCDvsSR_FR=ensureTFile(path+"TauTau.Step3_QCD_vsSR.root");
    _histo_name="fQCDvsSR_FR_"+_s_final_state.at(2);
    TauTau_fQCDvsSR_FR=extractTF1(f_TauTau_QCDvsSR_FR,_histo_name);

    f_LL_FR=ensureTFile(path+"LL.Step1_FakeRate.root");
    for (int i_bkg=0;i_bkg<2;i_bkg++) {
        _histo_name="f_FR_"+_s_fake_bkg.at(i_bkg==0?0:2)+"_"+_s_final_state.at(3);
        LL_f_FR[i_bkg]=extractTF1(f_LL_FR,_histo_name);
    }

    f_LL_Closure_FR=ensureTFile(path+"LL.Step2_Closure.root");
    for (int i_bkg=0;i_bkg<2;i_bkg++) {
        _histo_name="hClosure_FR_"+_s_fake_bkg.at(i_bkg==0?0:2)+"_"+_s_final_state.at(3);
        LL_hClosure_FR[i_bkg]=extractTH2(f_LL_Closure_FR,_histo_name);
    }

    f_LL_QCDvsSR_FR=ensureTFile(path+"LL.Step3_QCD_vsSR.root");
    _histo_name="hQCDvsSR_FR_"+_s_final_state.at(3);
    LL_hQCDvsSR_FR=extractTH2(f_LL_QCDvsSR_FR,_histo_name);

    f_LL_Frac=ensureTFile(path+"LL.Step4_Fraction.root");
    for (int i_bkg=0;i_bkg<2;i_bkg++) {
        for (int i_cat=0;i_cat<5;i_cat++) {
            _histo_name="hFrac_FR_"+_s_fake_bkg.at(i_bkg==0?0:2)+"_"+_s_final_state.at(3)+"_"+_s_categories.at(i_cat);
            LL_hFrac_FR[i_bkg][i_cat]=extractTH1(f_LL_Frac,_histo_name);
        }
    }

}

float FakeBkg::get_LTau_FR(int fs, int cat, int njet, float lpt, float taupt, float DR, float MT, float Mvis)
{
    int taupt_bin;
    if (taupt<=40) taupt_bin=0;
    else if (taupt<=50) taupt_bin=1;
    else taupt_bin=2;

    if (taupt>150) taupt=150;
    if (lpt>150) lpt=150;
    if (DR>6) DR=6;
    if (Mvis<20) Mvis=21;
    if (Mvis>300) Mvis=299;

    float FR[3],CL[3],vsSR[3],frac[3];
    for (int i_bkg=0;i_bkg<3;i_bkg++) {
        FR[i_bkg]=LTau_f_FR[i_bkg][fs][njet]->Eval(taupt);
        if (lpt<=23 && LTau_hClosure_FR[i_bkg][fs][taupt_bin]->GetBinContent(1)>0) {
            CL[i_bkg]=LTau_hClosure_FR[i_bkg][fs][taupt_bin]->GetBinContent(1);
        }
        else {
            CL[i_bkg]=LTau_fClosure_FR[i_bkg][fs][taupt_bin]->Eval(lpt);
        }
        frac[i_bkg]=LTau_hFrac_FR[i_bkg][fs][cat]->GetBinContent(LTau_hFrac_FR[i_bkg][fs][cat]->FindBin(Mvis));
    }

    vsSR[0]=LTau_fQCDvsSR_FR[fs]->Eval(DR);
    vsSR[1]=LTau_fWJvsSR_FR[fs]->Eval(MT);
    vsSR[2]=1.;

    float tot=0;
    for (int i_bkg=0;i_bkg<3;i_bkg++) {
        tot+=FR[i_bkg]/(1.-FR[i_bkg])*CL[i_bkg]*vsSR[i_bkg]*frac[i_bkg];
    }
    return tot;
}

float FakeBkg::get_TauTau_FR(int njet, float tau1pt, float tau2pt, float Mvis)
{
    if (tau2pt>150) tau2pt=150;
    if (tau1pt>150) tau1pt=150;
    if (Mvis<20) Mvis=21;
    if (Mvis>300) Mvis=299;

    float FR=TauTau_f_FR[njet]->Eval(tau2pt);
    float CL=TauTau_fClosure_FR->Eval(tau1pt);
    float vsSR=TauTau_fQCDvsSR_FR->Eval(Mvis);
    float tot=FR/(1.-FR)*CL*vsSR;
    return tot;
}
float FakeBkg::get_LL_FR(int cat, float DR, float l1pt, float l2pt, float Mvis)
{
    if (l1pt>150) l1pt=150;
    if (l2pt>150) l2pt=150;
    if (DR>6) DR=6;
    if (Mvis<20) Mvis=21;
    if (Mvis>300) Mvis=299;

    float FR[2],CL[2],vsSR[2],frac[2];
    for (int i_bkg=0;i_bkg<2;i_bkg++) {
        FR[i_bkg]=LL_f_FR[i_bkg]->Eval(DR);
        CL[i_bkg]=LL_hClosure_FR[i_bkg]->GetBinContent(LL_hClosure_FR[i_bkg]->FindBin(l1pt,l2pt));
        if (CL[i_bkg]<=0) CL[i_bkg]=1.;
        frac[i_bkg]=LL_hFrac_FR[i_bkg][cat]->GetBinContent(LL_hFrac_FR[i_bkg][cat]->FindBin(Mvis));
    }
    vsSR[0]=LL_hQCDvsSR_FR->GetBinContent(LL_hQCDvsSR_FR->FindBin(l1pt,l2pt));
    if (vsSR[0]<=0) vsSR[0]=1.;
    vsSR[1]=1.;

    float tot=0;
    for (int i_bkg=0;i_bkg<2;i_bkg++) {
        tot+=FR[i_bkg]*CL[i_bkg]*vsSR[i_bkg]*frac[i_bkg];
    }
    return tot;
}