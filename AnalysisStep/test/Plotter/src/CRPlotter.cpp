// Include classes
#include <HTauTauHMuMu/AnalysisStep/test/Plotter/include/CRPlotter.h>

// Constructor
//============================================================
CRPlotter::CRPlotter():Tree()
{    
    _current_process = -999;
    _current_final_state = -999;
    
    _s_final_state.push_back("4mu");
    _s_final_state.push_back("2muetau");
    _s_final_state.push_back("2mumutau");
    _s_final_state.push_back("2mutautau");
    _s_final_state.push_back("2e2mu");
    _s_final_state.push_back("2eetau");
    _s_final_state.push_back("2emutau");
    _s_final_state.push_back("2etautau");
    _s_final_state.push_back("4l");
    
    _s_zl_final_state.push_back("2mue");
    _s_zl_final_state.push_back("2mumu");
    _s_zl_final_state.push_back("2mutau");
    _s_zl_final_state.push_back("2ee");
    _s_zl_final_state.push_back("2emu");
    _s_zl_final_state.push_back("2etau");
    _s_zl_final_state.push_back("3l");
    
    _s_region.push_back("2P2F");
    _s_region.push_back("3P1F");
    _s_region.push_back("SS");
    
    _lumi=35.92;
    string year="2016Legacy";

    DeepTauSF_VSe_ETau = new TauIDSFTool(year,"DeepTau2017v2p1VSe","Medium");
    DeepTauSF_VSmu_ETau = new TauIDSFTool(year,"DeepTau2017v2p1VSmu","Tight");
    DeepTauSF_VSjet_ETau = new TauIDSFTool(year,"DeepTau2017v2p1VSjet","Medium");
    DeepTauSF_VSe_MuTau = new TauIDSFTool(year,"DeepTau2017v2p1VSe","VLoose");
    DeepTauSF_VSmu_MuTau = new TauIDSFTool(year,"DeepTau2017v2p1VSmu","Tight");
    DeepTauSF_VSjet_MuTau = new TauIDSFTool(year,"DeepTau2017v2p1VSjet","Tight");
    DeepTauSF_VSe_TauTau = new TauIDSFTool(year,"DeepTau2017v2p1VSe","VLoose");
    DeepTauSF_VSmu_TauTau = new TauIDSFTool(year,"DeepTau2017v2p1VSmu","Tight");
    DeepTauSF_VSjet_TauTau = new TauIDSFTool(year,"DeepTau2017v2p1VSjet","Tight");
    
//    DeepTauSF_VSe_bare = new TauIDSFTool(year,"DeepTau2017v2p1VSe","VVVLoose");
   DeepTauSF_VSmu_bare = new TauIDSFTool(year,"DeepTau2017v2p1VSmu","VLoose");
   DeepTauSF_VSjet_bare = new TauIDSFTool(year,"DeepTau2017v2p1VSjet","VVVLoose");
}
//--------------------------------------------------------------------------------------

// Destructor
//====================
CRPlotter::~CRPlotter()
{
}
//====================

//====================
void CRPlotter::SetVariable(string Var, string VarName, string X, string Y, int Nbins, float Range_min, float Range_max, int Type)
{
    var=Var;
    varName=VarName;
    xlabel=X;
    ylabel=Y;
    nbins=Nbins;
    range_min=Range_min;
    range_max=Range_max;
    varType=Type;
}

//====================
void CRPlotter::SetRegion(string Region)
{
    region=Region;
}

//====================
void CRPlotter::SetData(string path, string name)
{
    pathData=path;
    nameData=name;
}

//====================
void CRPlotter::SetSim(string path, string name)
{
    pathSim.push_back(path);
    nameSim.push_back(name);
}

//====================
void CRPlotter::DeclareHistos()
{
    DeclareHistosData();
    for (size_t i=0;i<=pathSim.size();i++) {
        DeclareHistosSim(i);
    }
    cout<<"[INFO] Simulation Histograms declared."<<endl;
}

//====================
void CRPlotter::FillHistos()
{
    FillHistosData();
    for (size_t i=0;i<pathSim.size();i++) {
        FillHistosSim(i);
    }
    SumHistos();
    SaveHistos();
}

//====================
void CRPlotter::Plotter(bool addData)
{
    GetHistos();
    int colors[]={kGreen,kMagenta,kBlue,kCyan,kRed,kYellow};
    for (int i_re=0;i_re<(region=="CRZLLTree"?3:1);i_re++) {
        for (int i_fs=0;i_fs<=(region=="CRZLLTree"?(Settings::fs4l):(Settings::fs3l));i_fs++) {
            TString name=var+"_"+(region=="CRZLLTree"?_s_region.at(i_re)+"_":"")+(region=="CRZLLTree"?_s_final_state.at(i_fs):_s_zl_final_state.at(i_fs));
            TCanvas *c=new TCanvas(name,name,0,0,1000,800);

            THStack *mc=new THStack(name,name);
            float x=0.55;//(i_fs==0 || i_fs==4)?0.2:0.55);
            TLegend *l=new TLegend( x, 0.4, x+0.35, 0.9 ,ToFSName(region=="CRZLLTree"?_s_final_state.at(i_fs):_s_zl_final_state.at(i_fs)));
            l->SetFillStyle(0);
            l->SetBorderSize(0);
            l->SetTextSize(0.05);
            vector<TH1F *> HList;

            for (size_t i_proc=0;i_proc<pathSim.size();i_proc++) {
                if (hSim[i_proc][i_re][i_fs]->Integral()>0) {
                    hSim[i_proc][i_re][i_fs]->SetStats(0);
                    SetColor(hSim[i_proc][i_re][i_fs],colors[i_proc]);
                    HList.push_back(hSim[i_proc][i_re][i_fs]);
                    TString labelName;
                    labelName=nameSim.at(i_proc);
                    l->AddEntry(hSim[i_proc][i_re][i_fs],labelName);
                }
            }


            for (int i=HList.size()-1;i>=0;i--)
                mc->Add(HList.at((size_t)i));

            if (addData) {
                hData[i_re][i_fs]->SetBinErrorOption(TH1::kPoisson);
                hData[i_re][i_fs]->SetStats(0);
                hData[i_re][i_fs]->SetMarkerSize(1.2);
                hData[i_re][i_fs]->SetMarkerStyle(8);
                SetColor(hData[i_re][i_fs],kBlack);
                l->AddEntry(hData[i_re][i_fs],"Data","lep");
            }

            TH1F *base=(TH1F *)hSim[0][i_re][i_fs]->Clone();
            base->Reset();
            base->SetStats(0);
            base->GetXaxis()->SetTitleSize(0.05);
            base->GetXaxis()->SetLabelSize(0.045);
            base->GetYaxis()->SetTitleSize(0.05);
            base->GetYaxis()->SetTitleOffset(0.9);
            base->GetYaxis()->SetLabelSize(0.045);
            if (addData)
                base->SetMaximum(max(hData[i_re][i_fs]->GetMaximum(),mc->GetMaximum())*1.3);
            else
                base->SetMaximum(mc->GetMaximum()*1.3);

            base->Draw();
            mc->Draw("histsame");
            if (addData)
                hData[i_re][i_fs]->Draw("Esame");
            l->Draw();

            CMS_lumi *lumi = new CMS_lumi;
            lumi->set_lumi(c, _lumi, 0);

            c->SetBottomMargin(0.12);
            c->SetGrid();
            c->SaveAs("Plots/"+region+"_"+name+(addData?"_data.png":".png"));
            c->SaveAs("Plots/"+region+"_"+name+(addData?"_data.pdf":".pdf"));
        }
    }
}


//=====================================private==========================================
//====================
void CRPlotter::DeclareHistosData()
{
    if (region=="CRZLLTree") {
        for (int i_re=0;i_re<3;i_re++) {
            for (int i_fs=0;i_fs<=Settings::fs4l;i_fs++) {
                _histo_name=var+"_"+nameData+"_"+(region=="CRZLLTree"?_s_region.at(i_re)+"_":"")+_s_final_state.at(i_fs);
                _histo_labels=";"+xlabel+ ";" + ylabel;
                hData[i_re][i_fs]=new TH1F(_histo_name, _histo_labels, nbins, range_min, range_max);
            }
        }
    }
    else {
        for (int i_fs=0;i_fs<=Settings::fs3l;i_fs++) {
            _histo_name=var+"_"+nameData+"_"+_s_zl_final_state.at(i_fs);
            _histo_labels=";"+xlabel+ ";" + ylabel;
            hData[0][i_fs]=new TH1F(_histo_name, _histo_labels, nbins, range_min, range_max);
        }
    }
    cout<<"[INFO] Data Histograms declared."<<endl;
}

//====================
void CRPlotter::DeclareHistosSim(size_t i)
{
    if (region=="CRZLLTree") {
        for (int i_re=0;i_re<3;i_re++) {
            for (int i_fs=0;i_fs<=Settings::fs4l;i_fs++) {
                _histo_name=var+"_"+(i<pathSim.size()? nameSim.at(i):"TotalSim")+"_"+(region=="CRZLLTree"?_s_region.at(i_re)+"_":"")+_s_final_state.at(i_fs);
                _histo_labels=";"+xlabel+ ";" + ylabel;
                hSim[i][i_re][i_fs]=new TH1F(_histo_name, _histo_labels, nbins, range_min, range_max);
            }
        }
    }
    else {
        for (int i_fs=0;i_fs<=Settings::fs3l;i_fs++) {
            _histo_name=var+"_"+(i<pathSim.size()? nameSim.at(i):"TotalSim")+"_"+_s_zl_final_state.at(i_fs);
            _histo_labels=";"+xlabel+ ";" + ylabel;
            hSim[i][0][i_fs]=new TH1F(_histo_name, _histo_labels, nbins, range_min, range_max);
        }
    }
}
//====================
void CRPlotter::FillHistosData()
{
    cout<<pathData<<endl;
    input_file = TFile::Open(pathData,"read");

    hCounters = (TH1F*)input_file->Get(region+"/Counters");
    gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);

    input_tree = (TTree*)input_file->Get(region+"/candTree");
    Init( input_tree, pathData, true);

    if (fChain == 0) {return;}

//     if (varType==1) fChain->SetBranchAddress(var, &varShort, &b_var);
//     else if (varType==2) fChain->SetBranchAddress(var, &varFloat, &b_var);
//     else {cout<<"[ERROR] variable type not correct: "<<varType<<endl; return;}
    
    Long64_t nentries = fChain->GetEntriesFast();
    cout<<"[INFO] Processing "<<pathData<<", total event: "<<nentries<<endl;
    Long64_t nbytes = 0, nb = 0;

    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry%50000==0) cout<<ientry<<endl;
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;
        
//         cout<<"step1"<<endl;
        if ( fabs(LepEta->at(2)) > 2.5) {continue;}
        if ( (abs(LepLepId->at(2))==15 && fabs(LepEta->at(2)) > 2.3)) {continue;}
        if ( (LepSIP->at(2) > 4. || Lepdxy->at(2) > 0.5 || Lepdz->at(2) > 1.0) && (abs(LepLepId->at(2)) == 11 || abs(LepLepId->at(2)) == 13) ) {continue;} // Included dxy/dz cuts for 3rd LEPTON (ele and mu)             
        if ( ( Lepdz->at(2) > 10 || LepPt->at(2) < 20 ||!LepisID->at(2)) && abs(LepLepId->at(2))==15 ) {continue;}

//         cout<<"step2"<<endl;
        if (region!="CRZLTree") {
//             cout<<"???"<<endl;
            if ( fabs(LepEta->at(3)) > 2.5) {continue;}
            if ( (abs(LepLepId->at(3))==15 && fabs(LepEta->at(3)) > 2.3)) {continue;}
            if ( (LepSIP->at(3) > 4. || Lepdxy->at(3) > 0.5 || Lepdz->at(3) > 1.0) && (abs(LepLepId->at(2)) == 11 || abs(LepLepId->at(2)) == 13) ) {continue;} // Included dxy/dz cuts for 4th LEPTON (ele and mu)
            if ( ( Lepdz->at(3) > 10 || LepPt->at(3) < 20 ||!LepisID->at(3)) && abs(LepLepId->at(3))==15 ) {continue;}
            if ( ZZMass < 70. ) continue;
        }
        
//         cout<<"step3"<<endl;
        if (region=="CRZLLTree") {
            if (test_bit(CRflag, CRZLLos_2P2F)) _current_region=0;
            else if (test_bit(CRflag, CRZLLos_3P1F)) _current_region=1;
            else if (test_bit(CRflag, CRZLLss)) _current_region=2;
            else continue;
        }
        else _current_region=0;
        
        _current_final_state = FindFinalState();
        if (_current_final_state<0) {continue;}
        
//         cout<<"Fill "<<LepPt->at(2)<<","<<_current_region<<","<<_current_final_state<<endl;
        hData[_current_region][_current_final_state]->Fill(LepPt->at(2));
//         if (varType==1) hData[_current_final_state]->Fill(varShort);
//         else if (varType==2) hData[_current_final_state]->Fill(varFloat);
    }
    input_file->Close();
    cout<<"[INFO] Processing "<<pathData<<" done"<<endl;
}

//====================
void CRPlotter::FillHistosSim(size_t i)
{
    cout<<i<<","<<pathSim.at(i)<<endl;
    input_file = TFile::Open(pathSim.at(i),"read");

    hCounters = (TH1F*)input_file->Get(region+"/Counters");
    gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);

    input_tree = (TTree*)input_file->Get(region+"/candTree");
    Init( input_tree, pathSim.at(i), true);

    if (fChain == 0) {return;}

//     if (varType==1) fChain->SetBranchAddress(var, &varShort, &b_var);
//     else if (varType==2) fChain->SetBranchAddress(var, &varFloat, &b_var);
//     else {cout<<"[ERROR] variable type not correct: "<<varType<<endl; return;}
    
    Long64_t nentries = fChain->GetEntriesFast();
    cout<<"[INFO] Processing "<<pathSim.at(i)<<", total event: "<<nentries<<endl;
    Long64_t nbytes = 0, nb = 0;

    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry%50000==0) cout<<ientry<<endl;
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;
               
        if ( fabs(LepEta->at(2)) > 2.5) {continue;}
        if ( (abs(LepLepId->at(2))==15 && fabs(LepEta->at(2)) > 2.3)) {continue;}
        if ( (LepSIP->at(2) > 4. || Lepdxy->at(2) > 0.5 || Lepdz->at(2) > 1.0) && (abs(LepLepId->at(2)) == 11 || abs(LepLepId->at(2)) == 13) ) {continue;} // Included dxy/dz cuts for 3rd LEPTON (ele and mu)             
        if ( ( Lepdz->at(2) > 10 || LepPt->at(2) < 20 ||!LepisID->at(2)) && abs(LepLepId->at(2))==15 ) {continue;}
        
        if (region!="CRZLTree") {
            if ( fabs(LepEta->at(3)) > 2.5) {continue;}
            if ( (abs(LepLepId->at(3))==15 && fabs(LepEta->at(3)) > 2.3)) {continue;}
            if ( (LepSIP->at(3) > 4. || Lepdxy->at(3) > 0.5 || Lepdz->at(3) > 1.0) && (abs(LepLepId->at(2)) == 11 || abs(LepLepId->at(2)) == 13) ) {continue;} // Included dxy/dz cuts for 4th LEPTON (ele and mu)
            if ( ( Lepdz->at(3) > 10 || LepPt->at(3) < 20 ||!LepisID->at(3)) && abs(LepLepId->at(3))==15 ) {continue;}
            if ( ZZMass < 70. ) continue;
        }
        
        if (region=="CRZLLTree") {
            if (test_bit(CRflag, CRZLLos_2P2F)) _current_region=0;
            else if (test_bit(CRflag, CRZLLos_3P1F)) _current_region=1;
            else if (test_bit(CRflag, CRZLLss)) _current_region=2;
            else continue;
        }
        _current_final_state = FindFinalState();
        if (_current_final_state<0) {continue;}
        
        _k_factor = calculate_K_factor(pathSim.at(i));
        _TauIDSF = (region!="CRZLTree"? calculate_TauIDSF(nameSim.at(i), false) : calculate_TauIDSF_OneLeg(TauDecayMode->at(2), LepPt->at(2), LepEta->at(2), GENMatch(2,nameSim.at(i)), 165, false) );
//         cout<<_TauIDSF<<endl;
        _event_weight = (_lumi * 1000 * xsec * _k_factor * _TauIDSF * overallEventWeight * L1prefiringWeight) / gen_sum_weights;
//         if (varType==1) hSim[i][_current_final_state]->Fill(varShort, _event_weight);
//         else if (varType==2) hSim[i][_current_final_state]->Fill(varFloat, _event_weight);
        hSim[i][_current_region][_current_final_state]->Fill(LepPt->at(2), _event_weight);
    }
    input_file->Close();
    cout<<"[INFO] Processing "<<pathSim.at(i)<<" done"<<endl;
}

//====================
void CRPlotter::SumHistos()
{
    cout<<"step0"<<endl;
    for (int i_re=0;i_re<(region=="CRZLLTree"?3:1);i_re++) {
        for (int i_fs=0;i_fs<(region=="CRZLLTree"?(Settings::fs4l):(Settings::fs3l));i_fs++) {
            hData[i_re][region=="CRZLLTree"?(Settings::fs4l):(Settings::fs3l)]->Add(hData[i_re][i_fs]);
        }
        cout<<"step1"<<endl;
        for (size_t i=0;i<pathSim.size();i++) {
            for (int i_fs=0;i_fs<(region=="CRZLLTree"?(Settings::fs4l):(Settings::fs3l));i_fs++) {
                hSim[i][i_re][region=="CRZLLTree"?(Settings::fs4l):(Settings::fs3l)]->Add(hSim[i][i_re][i_fs]);
            }
        }
        cout<<"step2"<<endl;
        for (size_t i=0;i<pathSim.size();i++) {
            for (int i_fs=0;i_fs<=(region=="CRZLLTree"?(Settings::fs4l):(Settings::fs3l));i_fs++) {
            hSim[pathSim.size()][i_re][i_fs]->Add(hSim[i][i_re][i_fs]);
            }
        }
    }
    cout<<"[INFO] Histos summed."<<endl;
}

//===============================================================================
void CRPlotter::SaveHistos()
{
    TString output_file_name="Histos"+region+".root";
    output_file = TFile::Open(output_file_name, "recreate");
    output_file->cd();
    for (int i_re=0;i_re<(region=="CRZLLTree"?3:1);i_re++) {
        for (int i_fs=0;i_fs<=(region=="CRZLLTree"?(Settings::fs4l):(Settings::fs3l));i_fs++) {
            hData[i_re][i_fs]->Write();
        }
        for (size_t i_proc=0;i_proc<=pathSim.size();i_proc++) {
            for (int i_fs=0;i_fs<=(region=="CRZLLTree"?(Settings::fs4l):(Settings::fs3l));i_fs++) {
                hSim[i_proc][i_re][i_fs]->Write();
            }
        }
    }
    output_file->Close();
    cout<<"[INFO] Saved Histos to "<<output_file_name.Data()<<endl;
}

//===============================================================================
void CRPlotter::GetHistos()
{
    TString output_file_name="Histos"+region+".root";
    output_file = TFile::Open(output_file_name, "read");
    if (region=="CRZLLTree") {
        for (int i_re=0;i_re<3;i_re++) {
            for (int i_fs=0;i_fs<=Settings::fs4l;i_fs++) {
                _histo_name=var+"_"+nameData+"_"+(region=="CRZLLTree"?_s_region.at(i_re)+"_":"")+_s_final_state.at(i_fs);
                hData[i_re][i_fs]=(TH1F*)output_file->Get(_histo_name);
            }
            for (size_t i_proc=0;i_proc<=pathSim.size();i_proc++) {
                for (int i_fs=0;i_fs<=Settings::fs4l;i_fs++) {
                    _histo_name=var+"_"+(i_proc<pathSim.size()? nameSim.at(i_proc):"TotalSim")+"_"+(region=="CRZLLTree"?_s_region.at(i_re)+"_":"")+_s_final_state.at(i_fs);
                    hSim[i_proc][i_re][i_fs]=(TH1F*)output_file->Get(_histo_name);
                }
            }
        }
    }
    else {
        for (int i_fs=0;i_fs<=Settings::fs3l;i_fs++) {
            _histo_name=var+"_"+nameData+"_"+_s_zl_final_state.at(i_fs);
            hData[0][i_fs]=(TH1F*)output_file->Get(_histo_name);
        }
        for (size_t i_proc=0;i_proc<=pathSim.size();i_proc++) {
            for (int i_fs=0;i_fs<=Settings::fs3l;i_fs++) {
                _histo_name=var+"_"+(i_proc<pathSim.size()? nameSim.at(i_proc):"TotalSim")+"_"+_s_zl_final_state.at(i_fs);
                hSim[i_proc][0][i_fs]=(TH1F*)output_file->Get(_histo_name);
            }
        }
    }
    cout<<"[INFO] Get histos from "<<output_file_name<<endl;
}

//====================
int CRPlotter::FindFinalState()
{
    int final_state = -999;

    if ( Z1Flav == -121 )
    {
        if ( abs(Z2Flav) == 169 ) final_state = Settings::fs2e2mu;
        else if ( abs(Z2Flav) == 165) final_state = Settings::fs2eetau;
        else if ( abs(Z2Flav) == 195) final_state = Settings::fs2emutau;
        else if ( abs(Z2Flav) == 225) final_state = Settings::fs2etautau;
        else if ( abs(Z2Flav) == 11)  final_state = Settings::fs2ee;
        else if ( abs(Z2Flav) == 13)  final_state = Settings::fs2emu;
        else if ( abs(Z2Flav) == 15)  final_state = Settings::fs2etau;
        else final_state = -1;
    }
    else if ( Z1Flav == -169 )
    {
        if ( abs(Z2Flav) == 169 ) final_state = Settings::fs4mu;
        else if ( abs(Z2Flav) == 165) final_state = Settings::fs2muetau;
        else if ( abs(Z2Flav) == 195) final_state = Settings::fs2mumutau;
        else if ( abs(Z2Flav) == 225) final_state = Settings::fs2mutautau;
        else if ( abs(Z2Flav) == 11)  final_state = Settings::fs2mue;
        else if ( abs(Z2Flav) == 13)  final_state = Settings::fs2mumu;
        else if ( abs(Z2Flav) == 15)  final_state = Settings::fs2mutau;
        else final_state = -1;
    }
    else
    {
        final_state = -1;
    }
   
    return final_state;
}


//=================================
float CRPlotter::calculate_K_factor(TString input_file_name)
{
   
   float k_factor = 1;
   
   if ( input_file_name.Contains("ZZTo4l"))
   {
      k_factor = KFactor_EW_qqZZ * KFactor_QCD_qqZZ_M; // As of Moriond2016
   }
   else if ( input_file_name.Contains("ggTo"))
   {
      k_factor = KFactor_QCD_ggZZ_Nominal; // as of Moriond2016
   }
   return k_factor;
}

//=================================
//DeepTau ID scale factors
int CRPlotter::GENMatch(int ileg,TString input_file_name)
{
   if (abs(LepLepId->at(ileg))!=15) {
//       cout<<"Only applicable to pdgId = 15, not "<<LepLepId->at(ileg)<<endl;
      return 6;
   }
   if ( input_file_name.Contains("WZTo") || input_file_name.Contains("ZZ") || input_file_name.Contains("H") || input_file_name.Contains("ggTo") ) {
      if (TauTES_p_Up->at(ileg)>0)
	 return 5;
      else if (TauFES_p_Up->at(ileg)>0)
	 return 1;
      else
	 return 2;
   }
   else {// if ( input_file_name.Contains("WW") || input_file_name.Contains("TTZ") ) {
      if (TauTES_p_Up->at(ileg)>0)
         return 5;
      else if (TauFES_p_Up->at(ileg)>0)
         return 1;
      else
         return 6;
   }
}

float CRPlotter::calculate_TauIDSF_OneLeg(short DM, float pt, float eta, int genmatch, int flav, bool tight)
{
    eta=fabs(eta);
    if (tight) {
        if (genmatch==6)
          return 1.;
        else if (genmatch==5) {
          //TauIDSFTool *DeepTauSF = new TauIDSFTool("2016Legacy","DeepTau2017v2p1VSjet",WP[0].Data());
          if (flav==165) return DeepTauSF_VSjet_ETau->getSFvsPT(pt,genmatch);
          else if (flav==195) return DeepTauSF_VSjet_MuTau->getSFvsPT(pt,genmatch);
          else if (flav==225) return DeepTauSF_VSjet_TauTau->getSFvsPT(pt,genmatch);
          else return 1.;
        }
        else if (genmatch==1 || genmatch==3) {
          //TauIDSFTool *DeepTauSF = new TauIDSFTool("2016Legacy","DeepTau2017v2p1VSe",WP[1].Data());
          if (flav==165) return DeepTauSF_VSe_ETau->getSFvsEta(eta,genmatch);
          else if (flav==195) return DeepTauSF_VSe_MuTau->getSFvsEta(eta,genmatch);
          else if (flav==225) return DeepTauSF_VSe_TauTau->getSFvsEta(eta,genmatch);
          else return 1.;
        }
        else if (genmatch==2 || genmatch==4) {
          //TauIDSFTool *DeepTauSF = new TauIDSFTool("2016Legacy","DeepTau2017v2p1VSmu",WP[2].Data());
          if (flav==165) return DeepTauSF_VSmu_ETau->getSFvsEta(eta,genmatch);
          else if (flav==195) return DeepTauSF_VSmu_MuTau->getSFvsEta(eta,genmatch);
          else if (flav==225) return DeepTauSF_VSmu_TauTau->getSFvsEta(eta,genmatch);
          else return 1.;
        }
        else {
          cout<<"The genmatch "<<genmatch<<" is wrong"<<endl;
          return 1.;
        }
    }
    else {
        if (genmatch==6)
          return 1.;
        else if (genmatch==5) return DeepTauSF_VSjet_bare->getSFvsPT(pt,genmatch);
        else if (genmatch==1 || genmatch==3) return 1.;//DeepTauSF_VSe_bare->getSFvsEta(eta,genmatch);
        else if (genmatch==2 || genmatch==4) return DeepTauSF_VSmu_bare->getSFvsEta(eta,genmatch);
        else {
          cout<<"The genmatch "<<genmatch<<" is wrong"<<endl;
          return 1.;
        }
    }
}

float CRPlotter::calculate_TauIDSF(TString input_file_name,bool tight)
{
    float SF = 1;
    for (size_t ileg=2;ileg<=3;ileg++) {
        if (abs(LepLepId->at(ileg))!=15)
            continue;
        int genmatch=GENMatch(ileg,input_file_name);
        SF=SF*calculate_TauIDSF_OneLeg(TauDecayMode->at(ileg),LepPt->at(ileg),LepEta->at(ileg),genmatch,abs(Z2Flav),tight);
    }
    return SF;
}

//===============================================================================
void CRPlotter::SetColor(TH1F *h,int color)
{
    h->SetMarkerColor(color);
    h->SetLineColor(color);
    h->SetFillColor(color);
}

//===============================================================================
TString CRPlotter::ToFSName(string name)
{
    TString n=name;
    n.ReplaceAll("mu","#mu");
    n.ReplaceAll("tau","#tau_{h}");
//     n.ReplaceAll("tautau","#tau_{h}#tau_{h}");
//     n.ReplaceAll("etau","#tau_{e}#tau_{h}");
//     n.ReplaceAll("#mutau","#tau_{#mu}#tau_{h}");
    return n;
}