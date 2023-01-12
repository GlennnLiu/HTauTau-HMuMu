// Include classes
#include <HTauTauHMuMu/AnalysisStep/test/Datacards/include/Datacards.h>

// Constructor
//============================================================
Datacards::Datacards():Tree()
{    
    _current_process = -999;
    _current_final_state = -999;
    _current_gen_final_state = -999;
    
    _fs_ROS_SS.push_back(0.954003);
    _fs_ROS_SS.push_back(1.05664);
    _fs_ROS_SS.push_back(1.03268);
    _fs_ROS_SS.push_back(0.983579);
    _fs_ROS_SS.push_back(0.985223);
    _fs_ROS_SS.push_back(1.04922);
    _fs_ROS_SS.push_back(1.02563);
    _fs_ROS_SS.push_back(0.995229);
    
    _cb_ss.push_back(1.11134383);
    _cb_ss.push_back(1.11641482);
    _cb_ss.push_back(1.138206556);
    _cb_ss.push_back(1.07851374);
    _cb_ss.push_back(1.157347331);
    _cb_ss.push_back(1.136168688);
    _cb_ss.push_back(1.126838945);
    _cb_ss.push_back(1.035982102);
    
    _FR_uncUp.push_back(0.273855443);
    _FR_uncUp.push_back(0.275892969);
    _FR_uncUp.push_back(0.26910266);
    _FR_uncUp.push_back(0.284858774);
    _FR_uncUp.push_back(0.263630405);
    _FR_uncUp.push_back(0.271129198);
    _FR_uncUp.push_back(0.271688517);
    _FR_uncUp.push_back(0.297240493);
    
    _FR_uncDn.push_back(0.274119559);
    _FR_uncDn.push_back(0.276396085);
    _FR_uncDn.push_back(0.269425182);
    _FR_uncDn.push_back(0.285388729);
    _FR_uncDn.push_back(0.263906581);
    _FR_uncDn.push_back(0.271628339);
    _FR_uncDn.push_back(0.271983693);
    _FR_uncDn.push_back(0.29789368);
        
    
    _s_process.push_back("H");
    _s_process.push_back("qqZZ");
    _s_process.push_back("ggZZ");
    _s_process.push_back("rare");
    _s_process.push_back("ZX");
    _s_process.push_back("data_obs");
    
    _s_final_state.push_back("4mu");
    _s_final_state.push_back("2muetau");
    _s_final_state.push_back("2mumutau");
    _s_final_state.push_back("2mutautau");
    _s_final_state.push_back("2e2mu");
    _s_final_state.push_back("2eetau");
    _s_final_state.push_back("2emutau");
    _s_final_state.push_back("2etautau");
    _s_final_state.push_back("4l");
    
    _s_gen_final_state.push_back("4mu");
    _s_gen_final_state.push_back("2mu2tau");
    _s_gen_final_state.push_back("2e2mu");
    _s_gen_final_state.push_back("2e2tau");
    _s_gen_final_state.push_back("bkg");
    _s_gen_final_state.push_back("all");
    
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
}
//--------------------------------------------------------------------------------------

// Destructor
//====================
Datacards::~Datacards()
{
}
//====================

//===============================================================================
void Datacards::SetPaths(string path, string savepath)
{
    _path=path;
    _savepath=savepath;
}
//===============================================================================
void Datacards::ToHistos(string* process, int n_proc, int i_proc)
{
    DeclareHistos(i_proc);
    
    for (int i=0;i<n_proc;i++) {
        TString input_file_name=_path+process[i]+(i_proc<=3?"/HTauTauHMuMu_unc.root":"/HTauTauHMuMu.root");
        cout<<i<<","<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");
        
        hCounters = (TH1F*)input_file->Get("ZZTree/Counters");
        gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);
        
        input_tree = (TTree*)input_file->Get("ZZTree/candTree");
        Init( input_tree, input_file_name, true);
        
        if (fChain == 0) {return;}
        
        Long64_t nentries = fChain->GetEntriesFast();
        cout<<"[INFO] Processing "<<input_file_name<<", total event: "<<nentries<<endl;
        Long64_t nbytes = 0, nb = 0;
        
        for (Long64_t jentry=0; jentry<nentries;jentry++) {
            Long64_t ientry = LoadTree(jentry);
            if (ientry%50000==0) cout<<ientry<<endl;
            if (ientry < 0) break;
            nb = fChain->GetEntry(jentry);
            nbytes += nb;
            
            if ( fabs(LepEta->at(2)) > 2.5 || fabs(LepEta->at(3)) > 2.5) {continue;}
            if ( (abs(LepLepId->at(2))==15 && fabs(LepEta->at(2)) > 2.3) || (abs(LepLepId->at(3))==15 && fabs(LepEta->at(3)) > 2.3) ) {continue;}
            if ( (LepSIP->at(2) > 4. || Lepdxy->at(2) > 0.5 || Lepdz->at(2) > 1.0) && (abs(LepLepId->at(2)) == 11 || abs(LepLepId->at(2)) == 13) ) {continue;} // Included dxy/dz cuts for 3rd LEPTON (ele and mu)             
            if ( (LepSIP->at(3) > 4. || Lepdxy->at(3) > 0.5 || Lepdz->at(3) > 1.0) && (abs(LepLepId->at(2)) == 11 || abs(LepLepId->at(2)) == 13) ) {continue;} // Included dxy/dz cuts for 4th LEPTON (ele and mu)
            if ( ( Lepdz->at(2) > 10 || LepPt->at(2) < 20 || TauDecayMode->at(2)==5 || TauDecayMode->at(2)==6 || !LepisID->at(2)) && abs(LepLepId->at(2))==15 ) {continue;}
            if ( ( Lepdz->at(3) > 10 || LepPt->at(3) < 20 || TauDecayMode->at(3)==5 || TauDecayMode->at(3)==6 || !LepisID->at(3)) && abs(LepLepId->at(3))==15 ) {continue;}
            if ( ZZMass < 70. ) continue;
            
            _current_final_state = FindFinalState();
            if (_current_final_state<0) {continue;}
            _current_gen_final_state = (i_proc>=3?(Settings::gfsbkg):FindGenFinalState());
            if (i_proc==Settings::Data) _event_weight = 1.;
            else {
                _k_factor = calculate_K_factor(input_file_name,"");
                _TauIDSF = calculate_TauIDSF(input_file_name,"");
                _event_weight = (_lumi * 1000 * xsec * _k_factor * _TauIDSF * overallEventWeight * L1prefiringWeight) / gen_sum_weights;
            }
            histos_nom[i_proc][_current_gen_final_state][_current_final_state]->Fill(Z2GoodMass, _event_weight);
            
            //-----------------------------------------------------------------
            //------------------------uncertainties----------------------------
            //-----------------------------------------------------------------
            if (i_proc!=Settings::Data) {
                //CMS_lumi2
                histos_unc[i_proc][_current_gen_final_state][_current_final_state][0][0]->Fill(Z2GoodMass, _event_weight*1.026);
                histos_unc[i_proc][_current_gen_final_state][_current_final_state][0][1]->Fill(Z2GoodMass, _event_weight*0.974);

                //FR: since FR is not filled here, all these are nominal
                histos_unc[i_proc][_current_gen_final_state][_current_final_state][1][0]->Fill(Z2GoodMass, _event_weight);
                histos_unc[i_proc][_current_gen_final_state][_current_final_state][1][1]->Fill(Z2GoodMass, _event_weight);
                
                //qqZZ_K
                if (i_proc==Settings::qqZZ) {
                    _event_weight_up = (_lumi * 1000 * xsec * calculate_K_factor(input_file_name,"Up") * _TauIDSF * overallEventWeight * L1prefiringWeight) / gen_sum_weights;
                    _event_weight_dn = (_lumi * 1000 * xsec * calculate_K_factor(input_file_name,"Dn") * _TauIDSF * overallEventWeight * L1prefiringWeight) / gen_sum_weights;
                }
                else _event_weight_up=_event_weight_dn=_event_weight;
                histos_unc[i_proc][_current_gen_final_state][_current_final_state][2][0]->Fill(Z2GoodMass, _event_weight_up);
                histos_unc[i_proc][_current_gen_final_state][_current_final_state][2][1]->Fill(Z2GoodMass, _event_weight_dn);
                
                //ggZZ_K
                if (i_proc==Settings::ggZZ) {
                    histos_unc_ggzz[_current_gen_final_state][_current_final_state][0][0]->Fill(Z2GoodMass, (_lumi * 1000 * xsec * calculate_K_factor(input_file_name,"UpPDF") * _TauIDSF * overallEventWeight * L1prefiringWeight) / gen_sum_weights );
                    histos_unc_ggzz[_current_gen_final_state][_current_final_state][1][0]->Fill(Z2GoodMass, (_lumi * 1000 * xsec * calculate_K_factor(input_file_name,"UpQCD") * _TauIDSF * overallEventWeight * L1prefiringWeight) / gen_sum_weights );
                    histos_unc_ggzz[_current_gen_final_state][_current_final_state][2][0]->Fill(Z2GoodMass, (_lumi * 1000 * xsec * calculate_K_factor(input_file_name,"UpAs") * _TauIDSF * overallEventWeight * L1prefiringWeight) / gen_sum_weights );
                    histos_unc_ggzz[_current_gen_final_state][_current_final_state][3][0]->Fill(Z2GoodMass, (_lumi * 1000 * xsec * calculate_K_factor(input_file_name,"UpPDFReplica") * _TauIDSF * overallEventWeight * L1prefiringWeight) / gen_sum_weights );
                    histos_unc_ggzz[_current_gen_final_state][_current_final_state][0][1]->Fill(Z2GoodMass, (_lumi * 1000 * xsec * calculate_K_factor(input_file_name,"DnPDF") * _TauIDSF * overallEventWeight * L1prefiringWeight) / gen_sum_weights );
                    histos_unc_ggzz[_current_gen_final_state][_current_final_state][1][1]->Fill(Z2GoodMass, (_lumi * 1000 * xsec * calculate_K_factor(input_file_name,"DnQCD") * _TauIDSF * overallEventWeight * L1prefiringWeight) / gen_sum_weights );
                    histos_unc_ggzz[_current_gen_final_state][_current_final_state][2][1]->Fill(Z2GoodMass, (_lumi * 1000 * xsec * calculate_K_factor(input_file_name,"DnAs") * _TauIDSF * overallEventWeight * L1prefiringWeight) / gen_sum_weights );
                    histos_unc_ggzz[_current_gen_final_state][_current_final_state][3][1]->Fill(Z2GoodMass, (_lumi * 1000 * xsec * calculate_K_factor(input_file_name,"DnPDFReplica") * _TauIDSF * overallEventWeight * L1prefiringWeight) / gen_sum_weights );
                }
                
                //pythia_scale
                histos_unc[i_proc][_current_gen_final_state][_current_final_state][4][0]->Fill(Z2GoodMass, _event_weight*PythiaWeight_isr_muR4*PythiaWeight_fsr_muR4);
                histos_unc[i_proc][_current_gen_final_state][_current_final_state][4][1]->Fill(Z2GoodMass, _event_weight*PythiaWeight_isr_muR0p25*PythiaWeight_fsr_muR0p25);
                
                //eff_e
                _event_weight_up=_event_weight_dn=_event_weight;
                for (size_t ilep=0;ilep<4;ilep++) {
                    if (abs(LepLepId->at(ilep))==11) {
                        _event_weight_up*=1+LepSF_Unc->at(ilep);
                        _event_weight_dn*=1-LepSF_Unc->at(ilep);
                    }
                }
                histos_unc[i_proc][_current_gen_final_state][_current_final_state][5][0]->Fill(Z2GoodMass, _event_weight_up);
                histos_unc[i_proc][_current_gen_final_state][_current_final_state][5][1]->Fill(Z2GoodMass, _event_weight_dn);
                
                //eff_mu
                _event_weight_up=_event_weight_dn=_event_weight;
                for (size_t ilep=0;ilep<4;ilep++) {
                    if (abs(LepLepId->at(ilep))==13) {
                        _event_weight_up*=1+LepSF_Unc->at(ilep);
                        _event_weight_dn*=1-LepSF_Unc->at(ilep);
                    }
                }
                histos_unc[i_proc][_current_gen_final_state][_current_final_state][6][0]->Fill(Z2GoodMass, _event_weight_up);
                histos_unc[i_proc][_current_gen_final_state][_current_final_state][6][1]->Fill(Z2GoodMass, _event_weight_dn);
                
                //ID_tau
                _event_weight_up = (_lumi * 1000 * xsec * _k_factor * calculate_TauIDSF(input_file_name,"Up") * overallEventWeight * L1prefiringWeight) / gen_sum_weights;
                _event_weight_dn = (_lumi * 1000 * xsec * _k_factor * calculate_TauIDSF(input_file_name,"Dn") * overallEventWeight * L1prefiringWeight) / gen_sum_weights;
                histos_unc[i_proc][_current_gen_final_state][_current_final_state][7][0]->Fill(Z2GoodMass, _event_weight_up);
                histos_unc[i_proc][_current_gen_final_state][_current_final_state][7][1]->Fill(Z2GoodMass, _event_weight_dn);
                
                if (process[i]=="VBFH125" || process[i]=="WplusH125")
                    for (int i_unc=8;i_unc<13;i_unc++) {
                        histos_unc[i_proc][_current_gen_final_state][_current_final_state][i_unc][0]->Fill(Z2GoodMass, _event_weight);
                        histos_unc[i_proc][_current_gen_final_state][_current_final_state][i_unc][1]->Fill(Z2GoodMass, _event_weight);
                    }
                else {
                //tes
                histos_unc[i_proc][_current_gen_final_state][_current_final_state][8][0]->Fill(Z2GoodMass_TESup, _event_weight);
                histos_unc[i_proc][_current_gen_final_state][_current_final_state][8][1]->Fill(Z2GoodMass_TESdn, _event_weight);

                //ees
                histos_unc[i_proc][_current_gen_final_state][_current_final_state][9][0]->Fill(Z2GoodMass_EESup, _event_weight);
                histos_unc[i_proc][_current_gen_final_state][_current_final_state][9][1]->Fill(Z2GoodMass_EESdn, _event_weight);

                //mes
                histos_unc[i_proc][_current_gen_final_state][_current_final_state][10][0]->Fill(Z2GoodMass_MESup, _event_weight);
                histos_unc[i_proc][_current_gen_final_state][_current_final_state][10][1]->Fill(Z2GoodMass_MESdn, _event_weight);
      
                //jes
                histos_unc[i_proc][_current_gen_final_state][_current_final_state][11][0]->Fill(Z2GoodMass_JESup, _event_weight);
                histos_unc[i_proc][_current_gen_final_state][_current_final_state][11][1]->Fill(Z2GoodMass_JESdn, _event_weight);

                //jer
                histos_unc[i_proc][_current_gen_final_state][_current_final_state][12][0]->Fill(Z2GoodMass_JERup, _event_weight);
                histos_unc[i_proc][_current_gen_final_state][_current_final_state][12][1]->Fill(Z2GoodMass_JERdn, _event_weight);
                }
                
            }
            //-----------------------------------------------------------------
            //----------------------END uncertainties--------------------------
            //-----------------------------------------------------------------
        }
        input_file->Close();
        cout<<"[INFO] Processing "<<input_file_name<<" done"<<endl;
    }
}

void Datacards::ToZXHistos(string* process, int i_proc, TString input_file_FR_name)
{
    if (i_proc!=Settings::ZX) return;
    
    DeclareHistos(i_proc);
    
    TString input_file_name=_path+process[0]+"/HTauTauHMuMu.root";
    
    FakeRates *FR = new FakeRates( input_file_FR_name );
    
    input_file = TFile::Open( input_file_name);
    input_tree = (TTree*)input_file->Get("CRZLLTree/candTree");
    Init( input_tree, input_file_name , true);
    
    if (fChain == 0) return;
    
    Long64_t nentries = fChain->GetEntriesFast();
    cout<<"[INFO] Processing "<<input_file_name<<", total event: "<<nentries<<endl;
    Long64_t nbytes = 0, nb = 0;
    
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        
        Long64_t ientry = LoadTree(jentry);
        if (ientry%50000==0) cout<<ientry<<endl;
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;
        
        if ( !CRflag ) continue;
        if ( !test_bit(CRflag, CRZLLss) ) continue;
        
        if ( fabs(LepEta->at(2)) > 2.5 || fabs(LepEta->at(3)) > 2.5) {continue;}
        if ( (abs(LepLepId->at(2))==15 && fabs(LepEta->at(2)) > 2.3) || (abs(LepLepId->at(3))==15 && fabs(LepEta->at(3)) > 2.3) ) {continue;}
        if ( (LepSIP->at(2) > 4. || Lepdxy->at(2) > 0.5 || Lepdz->at(2) > 1.0) && (abs(LepLepId->at(2)) == 11 || abs(LepLepId->at(2)) == 13) ) {continue;} // Included dxy/dz cuts for 3rd LEPTON (ele and mu)             
        if ( (LepSIP->at(3) > 4. || Lepdxy->at(3) > 0.5 || Lepdz->at(3) > 1.0) && (abs(LepLepId->at(2)) == 11 || abs(LepLepId->at(2)) == 13) ) {continue;} // Included dxy/dz cuts for 4th LEPTON (ele and mu)
        if ( ( Lepdz->at(2) > 10 || LepPt->at(2) < 20 || !LepisID->at(2)) && abs(LepLepId->at(2))==15 ) {continue;}
        if ( ( Lepdz->at(3) > 10 || LepPt->at(3) < 20 || !LepisID->at(3)) && abs(LepLepId->at(3))==15 ) {continue;}
        if ( ZZMass < 70. ) continue;
        
        _current_final_state = FindFinalState();
        if (_current_final_state<0) {continue;}
        _current_gen_final_state = Settings::gfsbkg;
    
        int tauChannel=-1;
        if (abs(Z2Flav)==165) tauChannel=0;
        else if (abs(Z2Flav)==195) tauChannel=1;
        else if (abs(Z2Flav)==225) tauChannel=2;
        
        _yield_SR    = _fs_ROS_SS.at(_current_final_state)*_cb_ss.at(_current_final_state)*FR->GetFakeRate(LepPt->at(2),PFMET,LepEta->at(2),TauDecayMode->at(2),LepLepId->at(2),tauChannel)*FR->GetFakeRate(LepPt->at(3),PFMET,LepEta->at(3),TauDecayMode->at(3),LepLepId->at(3),tauChannel);
        
        histos_nom[i_proc][_current_gen_final_state][_current_final_state]->Fill(Z2GoodMass, _yield_SR);
        
    }
    for (int i_fs=0;i_fs<Settings::fs4l;i_fs++) {
        for (int i_unc=0;i_unc<num_of_unc;i_unc++) {
            if (i_unc==1) {
                histos_unc[i_proc][Settings::gfsbkg][i_fs][i_unc][0]->Add(histos_nom[i_proc][Settings::gfsbkg][i_fs],1+_FR_uncUp.at(i_fs));
                histos_unc[i_proc][Settings::gfsbkg][i_fs][i_unc][1]->Add(histos_nom[i_proc][Settings::gfsbkg][i_fs],1-_FR_uncDn.at(i_fs));
            }
            else {
                histos_unc[i_proc][Settings::gfsbkg][i_fs][i_unc][0]->Add(histos_nom[i_proc][Settings::gfsbkg][i_fs]);
                histos_unc[i_proc][Settings::gfsbkg][i_fs][i_unc][1]->Add(histos_nom[i_proc][Settings::gfsbkg][i_fs]);
            }
        }
    }
    input_file->Close();
    
    cout<<"[INFO] Processing "<<input_file_name<<" done"<<endl;
        
}
//===============================================================================
void Datacards::ComputeUncValues()
{
    for (int i_fs=0;i_fs<Settings::fs4l;i_fs++) {
        for (int i_proc=0;i_proc<Settings::Data;i_proc++) {
            if (i_proc<3)
                for (int i_gfs=0;i_gfs<=Settings::gfsbkg;i_gfs++) {
                    float nom=histos_nom[i_proc][i_gfs][i_fs]->Integral(), nom_mean=histos_nom[i_proc][i_gfs][i_fs]->GetMean();
                    if (nom<=0) {cout<<_s_process.at(i_proc)+"_"+_s_gen_final_state.at(i_gfs)+", "+_s_final_state.at(i_fs)<<endl;continue;}
                    for (int i_unc=0;i_unc<num_of_unc;i_unc++) {
                        if (i_unc==3) {
                            if (i_proc==Settings::ggZZ) {
                                float up, dn;
                                float tot_up=0, tot_dn=0;
                                for (int i=0;i<4;i++) {
                                    up=histos_unc_ggzz[i_gfs][i_fs][i][0]->Integral();
                                    dn=histos_unc_ggzz[i_gfs][i_fs][i][1]->Integral();
                                    tot_up+=(nom-up)*(nom-up);
                                    tot_dn+=(nom-dn)*(nom-dn);
                                }//end loop on ggZZ uncertainties
                                UncValues[i_proc][i_gfs][i_fs][i_unc][0]=1.+sqrt(tot_up)/nom;
                                UncValues[i_proc][i_gfs][i_fs][i_unc][1]=1.-sqrt(tot_dn)/nom;
                            }//end if
                            else {
                                UncValues[i_proc][i_gfs][i_fs][i_unc][0]=1.;
                                UncValues[i_proc][i_gfs][i_fs][i_unc][1]=1.;
                            }
                        }//end i_unc==3
                        else if (i_unc<8) {
                            UncValues[i_proc][i_gfs][i_fs][i_unc][0]=histos_unc[i_proc][i_gfs][i_fs][i_unc][0]->Integral()/nom;
                            UncValues[i_proc][i_gfs][i_fs][i_unc][1]=histos_unc[i_proc][i_gfs][i_fs][i_unc][1]->Integral()/nom;
                        }
                        else {
                            UncValues[i_proc][i_gfs][i_fs][i_unc][0]=histos_unc[i_proc][i_gfs][i_fs][i_unc][0]->GetMean()/nom_mean;
                            UncValues[i_proc][i_gfs][i_fs][i_unc][1]=histos_unc[i_proc][i_gfs][i_fs][i_unc][1]->GetMean()/nom_mean;
                        }
                    }//end loop on uncertainties
                }//end loop on gen final states and i_proc<3
            else {
                int i_gfs=Settings::gfsbkg;
                float nom=histos_nom[i_proc][i_gfs][i_fs]->Integral(), nom_mean=histos_nom[i_proc][i_gfs][i_fs]->GetMean();
                if (nom<=0) continue;
                for (int i_unc=0;i_unc<num_of_unc;i_unc++) {
                    if (i_unc==3) {
                        UncValues[i_proc][i_gfs][i_fs][i_unc][0]=1.;
                        UncValues[i_proc][i_gfs][i_fs][i_unc][1]=1.;
                    }
                    else if (i_unc<8) {
                        UncValues[i_proc][i_gfs][i_fs][i_unc][0]=histos_unc[i_proc][i_gfs][i_fs][i_unc][0]->Integral()/nom;
                        UncValues[i_proc][i_gfs][i_fs][i_unc][1]=histos_unc[i_proc][i_gfs][i_fs][i_unc][1]->Integral()/nom;
                    }
                    else {
                        UncValues[i_proc][i_gfs][i_fs][i_unc][0]=histos_unc[i_proc][i_gfs][i_fs][i_unc][0]->GetMean()/nom_mean;
                        UncValues[i_proc][i_gfs][i_fs][i_unc][1]=histos_unc[i_proc][i_gfs][i_fs][i_unc][1]->GetMean()/nom_mean;
                    }
                }//end loop on uncertainties
            }//end else
        }//end loop on processes
    }//end loop on final states
}

//===============================================================================
void Datacards::PrintUncValues()
{
    for (int i_fs=0;i_fs<Settings::fs4l;i_fs++) {
        cout<<"--------------------------------------"<<endl<<"Channel: "<<_s_final_state.at(i_fs)<<endl<<"--------------------------------------"<<endl;
        
        float tot_yield=0;
        for (int i_proc=0;i_proc<Settings::Data;i_proc++) {
            if (i_proc<3) {
                for (int i_gfs=0;i_gfs<=Settings::gfsbkg;i_gfs++) 
                    tot_yield+=histos_nom[i_proc][i_gfs][i_fs]->Integral();
            }
            else tot_yield+=histos_nom[i_proc][Settings::gfsbkg][i_fs]->Integral();
        }
        
        cout<<"Process\t";
        for (int i_proc=0;i_proc<Settings::Data;i_proc++) {
            if (i_proc<3) {
                for (int i_gfs=0;i_gfs<=Settings::gfsbkg;i_gfs++) {
                    if (histos_nom[i_proc][i_gfs][i_fs]->Integral()>0.001*tot_yield) 
                        cout<<_s_process.at(i_proc)+"_"+_s_gen_final_state.at(i_gfs)+"\t";
                }
            }
            else {
                if (histos_nom[i_proc][Settings::gfsbkg][i_fs]->Integral()>0.001*tot_yield) 
                    cout<<_s_process.at(i_proc)+"\t";
            }
        }
        cout<<endl;
        for (int i_unc=0;i_unc<num_of_unc;i_unc++) {
            cout<<Unc[i_unc]<<"\t";
            for (int i_proc=0;i_proc<Settings::Data;i_proc++) {
                if (i_proc<3) {
                    for (int i_gfs=0;i_gfs<=Settings::gfsbkg;i_gfs++) {
                        if (histos_nom[i_proc][i_gfs][i_fs]->Integral()>0.001*tot_yield) 
                            cout<<Form("%4.3f/%4.3f\t",UncValues[i_proc][i_gfs][i_fs][i_unc][0],UncValues[i_proc][i_gfs][i_fs][i_unc][1]);
                    }
                }
                else {
                    if (histos_nom[i_proc][Settings::gfsbkg][i_fs]->Integral()>0.001*tot_yield) 
                        cout<<Form("%4.3f/%4.3f\t",UncValues[i_proc][Settings::gfsbkg][i_fs][i_unc][0],UncValues[i_proc][Settings::gfsbkg][i_fs][i_unc][1]);
                }
            }
            cout<<endl;
        }
        cout<<"--------------------------------------"<<endl<<endl;
    }
}

//===============================================================================
void Datacards::ToDatacards()
{
    for (int i_fs=0;i_fs<Settings::fs4l;i_fs++) {
        
        ofstream card(_savepath+_s_final_state.at(i_fs)+"2016.txt");
        
        card<<"Datacard for channel: ch_"<<_s_final_state.at(i_fs)<<endl;
        card<<"imax 1 number of bins"<<endl;
        card<<"jmax * number of processes minus 1"<<endl;
        card<<"kmax * number of nuisance parameters"<<endl;
        card<<"-----------------------------------------"<<endl<<endl;
        
        card<<"shapes * ch_"<<_s_final_state.at(i_fs)<<" "<<_s_final_state.at(i_fs)<<"2016.root $PROCESS $PROCESS_$SYSTEMATIC"<<endl;
        card<<"-----------------------------------------"<<endl;
        card<<"bin ch_"<<_s_final_state.at(i_fs)<<endl;
        card<<"observation -1"<<endl;
        card<<"-----------------------------------------"<<endl<<endl;
        
        float tot_yield=0;
        for (int i_proc=0;i_proc<Settings::Data;i_proc++) {
            if (i_proc<3) {
                for (int i_gfs=0;i_gfs<=Settings::gfsbkg;i_gfs++) 
                    tot_yield+=histos_nom[i_proc][i_gfs][i_fs]->Integral();
            }
            else tot_yield+=histos_nom[i_proc][Settings::gfsbkg][i_fs]->Integral();
        }
        
        card<<"bin\t";
        for (int i_proc=0;i_proc<Settings::Data;i_proc++) {
            if (i_proc<3) {
                for (int i_gfs=0;i_gfs<=Settings::gfsbkg;i_gfs++) {
                    if (histos_nom[i_proc][i_gfs][i_fs]->Integral()>0.001*tot_yield) 
                        card<<"ch_"+_s_final_state.at(i_fs)<<"\t";
                }
            }
            else {
                if (histos_nom[i_proc][Settings::gfsbkg][i_fs]->Integral()>0.001*tot_yield) 
                    card<<"ch_"+_s_final_state.at(i_fs)<<"\t";
            }
        }//only write the processes taking up at least 0.1%
        card<<endl;
        
        card<<"process\t";
        for (int i_proc=0;i_proc<Settings::Data;i_proc++) {
            if (i_proc<3) {
                for (int i_gfs=0;i_gfs<=Settings::gfsbkg;i_gfs++) {
                    if (histos_nom[i_proc][i_gfs][i_fs]->Integral()>0.001*tot_yield) 
                        card<<_s_process.at(i_proc)+"_"+_s_gen_final_state.at(i_gfs)<<"\t";
                }
            }
            else {
                if (histos_nom[i_proc][Settings::gfsbkg][i_fs]->Integral()>0.001*tot_yield) 
                    card<<_s_process.at(i_proc)<<"\t";
            }
        }
        card<<endl;
        
        card<<"process\t";
        for (int i_proc=0;i_proc<Settings::Data;i_proc++) {
            if (i_proc<3) {
                for (int i_gfs=0;i_gfs<=Settings::gfsbkg;i_gfs++) {
                    if (histos_nom[i_proc][i_gfs][i_fs]->Integral()>0.001*tot_yield) 
                        card<<Form("%d\t",i_gfs==Settings::gfsbkg? (i_proc*5+1):(-i_proc*5-i_gfs-1));
                }
            }
            else {
                if (histos_nom[i_proc][Settings::gfsbkg][i_fs]->Integral()>0.001*tot_yield) 
                    card<<Form("%d\t",i_proc*5+1);
            }
        }
        card<<endl;
        
        card<<"rate\t";
        for (int i_proc=0;i_proc<Settings::Data;i_proc++) {
            if (i_proc<3) {
                for (int i_gfs=0;i_gfs<=Settings::gfsbkg;i_gfs++) {
                    if (histos_nom[i_proc][i_gfs][i_fs]->Integral()>0.001*tot_yield) 
                        card<<Form("%.3f\t",histos_nom[i_proc][i_gfs][i_fs]->Integral());
                }
            }
            else {
                if (histos_nom[i_proc][Settings::gfsbkg][i_fs]->Integral()>0.001*tot_yield) 
                    card<<Form("%.3f\t",histos_nom[i_proc][Settings::gfsbkg][i_fs]->Integral());
            }
        }
        card<<endl;
        
        for (int i_unc=0;i_unc<5;i_unc++) {//norm uncertainties
            if (i_unc==1) card<<Unc[i_unc]+"_"+_s_final_state.at(i_fs)+"\tlnN\t";
            else card<<Unc[i_unc]<<"\tlnN\t";
            for (int i_proc=0;i_proc<Settings::Data;i_proc++) {
                if (i_proc<3)
                    for (int i_gfs=0;i_gfs<=Settings::gfsbkg;i_gfs++) {
                        if (histos_nom[i_proc][i_gfs][i_fs]->Integral()>0.001*tot_yield) {
                            if (fabs(1-UncValues[i_proc][i_gfs][i_fs][i_unc][0])>0.00045 || fabs(1-UncValues[i_proc][i_gfs][i_fs][i_unc][1])>0.00045) 
                                card<<Form("%4.3f/%4.3f\t",UncValues[i_proc][i_gfs][i_fs][i_unc][0],UncValues[i_proc][i_gfs][i_fs][i_unc][1]);
                            else card<<"-\t";
                        }
                    }
                else if (histos_nom[i_proc][Settings::gfsbkg][i_fs]->Integral()>0.001*tot_yield) {
                    if (fabs(1-UncValues[i_proc][Settings::gfsbkg][i_fs][i_unc][0])>0.00045 || fabs(1-UncValues[i_proc][Settings::gfsbkg][i_fs][i_unc][1])>0.00045) 
                        card<<Form("%4.3f/%4.3f\t",UncValues[i_proc][Settings::gfsbkg][i_fs][i_unc][0],UncValues[i_proc][Settings::gfsbkg][i_fs][i_unc][1]);
                    else card<<"-\t";
                }
            }
            card<<endl;
        }
        
        for (int i_unc=5;i_unc<13;i_unc++) {
            card<<Unc[i_unc]<<"\tshape\t";
            for (int i_proc=0;i_proc<Settings::Data;i_proc++) {
                if (i_proc<3) {
                    for (int i_gfs=0;i_gfs<=Settings::gfsbkg;i_gfs++) {
                        if (histos_nom[i_proc][i_gfs][i_fs]->Integral()>0.001*tot_yield) card<<"1.000\t";
                    }
                }
                else {
                    if (histos_nom[i_proc][Settings::gfsbkg][i_fs]->Integral()>0.001*tot_yield) card<<"1.000\t";
                }
            }
            card<<endl;
        }
        
        card.close();
        
        TString output_file_name=_savepath+_s_final_state.at(i_fs)+"2016.root";
        output_file=TFile::Open(output_file_name,"recreate");
        output_file->cd();
        
        for (int i_proc=0;i_proc<=Settings::Data;i_proc++) {
            if (i_proc<3) {
                for (int i_gfs=0;i_gfs<=Settings::gfsbkg;i_gfs++) {
                    if (histos_nom[i_proc][i_gfs][i_fs]->Integral()>0.001*tot_yield) {
                        histos_nom[i_proc][i_gfs][i_fs]->Write((_s_process.at(i_proc)+"_"+_s_gen_final_state.at(i_gfs)).c_str());
                        for (int i_unc=5;i_unc<13;i_unc++) {
                            histos_unc[i_proc][i_gfs][i_fs][i_unc][0]->Write((_s_process.at(i_proc)+"_"+_s_gen_final_state.at(i_gfs)+"_"+Unc[i_unc]+"Up").c_str());
                            histos_unc[i_proc][i_gfs][i_fs][i_unc][1]->Write((_s_process.at(i_proc)+"_"+_s_gen_final_state.at(i_gfs)+"_"+Unc[i_unc]+"Down").c_str());
                        }
                    }
                }
            }
            else if (i_proc<Settings::Data) {
                if (histos_nom[i_proc][Settings::gfsbkg][i_fs]->Integral()>0.001*tot_yield) {
                    int i_gfs=Settings::gfsbkg;
                    histos_nom[i_proc][i_gfs][i_fs]->Write((_s_process.at(i_proc)).c_str());
                    for (int i_unc=5;i_unc<13;i_unc++) {
                        histos_unc[i_proc][i_gfs][i_fs][i_unc][0]->Write((_s_process.at(i_proc)+"_"+Unc[i_unc]+"Up").c_str());
                        histos_unc[i_proc][i_gfs][i_fs][i_unc][1]->Write((_s_process.at(i_proc)+"_"+Unc[i_unc]+"Down").c_str());
                    }
                }
            }
            else {
                int i_gfs=Settings::gfsbkg;
                histos_nom[i_proc][i_gfs][i_fs]->Write((_s_process.at(i_proc)).c_str());
            }
        }
        
        output_file->Close();
    }
}

//===============================================================================
void Datacards::DeclareHistos(int i_proc)
{
    for (int i_fs=0;i_fs<Settings::fs4l;i_fs++) {
        if (i_proc<3)
            for (int i_gfs=0;i_gfs<=Settings::gfsbkg;i_gfs++) {
                _histo_name=_s_process.at(i_proc)+"_"+_s_gen_final_state.at(i_gfs)+"_"+_s_final_state.at(i_fs);
                _histo_labels=";M_{Z_{prob}} (GeV);Weighted events / 10 GeV";
                histos_nom[i_proc][i_gfs][i_fs]=new TH1F(_histo_name,_histo_labels,14,0,140);
                for (int i_unc=0;i_unc<num_of_unc;i_unc++) {
                    _histo_name=_s_process.at(i_proc)+"_"+_s_gen_final_state.at(i_gfs)+"_"+Unc[i_unc]+"Up"+"_"+_s_final_state.at(i_fs);
                    _histo_labels=";M_{Z_{prob}} (GeV);Weighted events / 10 GeV";
                    histos_unc[i_proc][i_gfs][i_fs][i_unc][0]=new TH1F(_histo_name,_histo_labels,14,0,140);
                    _histo_name=_s_process.at(i_proc)+"_"+_s_gen_final_state.at(i_gfs)+"_"+Unc[i_unc]+"Down"+"_"+_s_final_state.at(i_fs);
                    _histo_labels=";M_{Z_{prob}} (GeV);Weighted events / 10 GeV";
                    histos_unc[i_proc][i_gfs][i_fs][i_unc][1]=new TH1F(_histo_name,_histo_labels,14,0,140);
                }//end loop on uncertainties
                if (i_proc==Settings::ggZZ)
                    for (int i_unc=0;i_unc<4;i_unc++) {
                        _histo_name=_s_process.at(i_proc)+"_"+_s_gen_final_state.at(i_gfs)+"_"+Form("%d",i_unc)+"Up"+"_"+_s_final_state.at(i_fs);
                        _histo_labels=";M_{Z_{prob}} (GeV);Weighted events / 10 GeV";
                        histos_unc_ggzz[i_gfs][i_fs][i_unc][0]=new TH1F(_histo_name,_histo_labels,14,0,140);
                        _histo_name=_s_process.at(i_proc)+"_"+_s_gen_final_state.at(i_gfs)+"_"+Form("%d",i_unc)+"Down"+"_"+_s_final_state.at(i_fs);
                        _histo_labels=";M_{Z_{prob}} (GeV);Weighted events / 10 GeV";
                        histos_unc_ggzz[i_gfs][i_fs][i_unc][1]=new TH1F(_histo_name,_histo_labels,14,0,140);
                    }//end loop on ggZZ uncertainties
            }//end loop on gen final states
        else {
            int i_gfs=Settings::gfsbkg;
            _histo_name=_s_process.at(i_proc)+"_"+_s_final_state.at(i_fs);
            _histo_labels=";M_{Z_{prob}} (GeV);Weighted events / 10 GeV";
            histos_nom[i_proc][i_gfs][i_fs]=new TH1F(_histo_name,_histo_labels,14,0,140);
            for (int i_unc=0;i_unc<num_of_unc;i_unc++) {
                _histo_name=_s_process.at(i_proc)+"_"+Unc[i_unc]+"Up"+"_"+_s_final_state.at(i_fs);
                _histo_labels=";M_{Z_{prob}} (GeV);Weighted events / 10 GeV";
                histos_unc[i_proc][i_gfs][i_fs][i_unc][0]=new TH1F(_histo_name,_histo_labels,14,0,140);
                _histo_name=_s_process.at(i_proc)+"_"+Unc[i_unc]+"Down"+"_"+_s_final_state.at(i_fs);
                _histo_labels=";M_{Z_{prob}} (GeV);Weighted events / 10 GeV";
                histos_unc[i_proc][i_gfs][i_fs][i_unc][1]=new TH1F(_histo_name,_histo_labels,14,0,140);
            }//end loop on uncertainties
        }//end else
    }//end loop on final states
}

//===============================================================================
int Datacards::FindFinalState()
{
    int final_state = -999;
    if ( Z1Flav == -121 )
    {
        if ( abs(Z2Flav) == 169 ) final_state = Settings::fs2e2mu;
        else if (abs(Z2Flav) == 165) final_state = Settings::fs2eetau;
        else if (abs(Z2Flav) == 195) final_state = Settings::fs2emutau;
        else if (abs(Z2Flav) == 225) final_state = Settings::fs2etautau;
        else final_state = -2;
    }
    else if ( Z1Flav == -169 )
    {
        if ( abs(Z2Flav) == 169 ) final_state = Settings::fs4mu;
        else if (abs(Z2Flav) == 165) final_state = Settings::fs2muetau;
        else if (abs(Z2Flav) == 195) final_state = Settings::fs2mumutau;
        else if (abs(Z2Flav) == 225) final_state = Settings::fs2mutautau;
        else final_state = -2;
    }
    else
    {
      final_state = -1;
    }
    return final_state;
}

//===============================================================================
int Datacards::FindGenFinalState()
{
    int gen_final_state = -999;
    if (GenZ1Flav*GenZ2Flav==169*169) gen_final_state=Settings::gfs4mu;
    else if (GenZ1Flav*GenZ2Flav==169*225) gen_final_state=Settings::gfs2mu2tau;
    else if (GenZ1Flav*GenZ2Flav==121*169) gen_final_state=Settings::gfs2e2mu;
    else if (GenZ1Flav*GenZ2Flav==121*225) gen_final_state=Settings::gfs2e2tau;
    else gen_final_state=Settings::gfsbkg;
    return gen_final_state;
}

//===============================================================================
float Datacards::calculate_K_factor(TString input_file_name, const string& unc="")
{
    float k_factor = 1;
    if ( input_file_name.Contains("ZZTo4l")) {
        if (unc=="Up") k_factor = KFactor_EW_qqZZ * KFactor_QCD_qqZZ_M*1.001;
        else if (unc=="Dn") k_factor = KFactor_EW_qqZZ * KFactor_QCD_qqZZ_M*0.999;
        else k_factor = KFactor_EW_qqZZ * KFactor_QCD_qqZZ_M;
    }
    else if ( input_file_name.Contains("ggTo")) {
        if (unc=="UpPDF") k_factor = KFactor_QCD_ggZZ_PDFScaleUp;
        else if (unc=="UpQCD") k_factor = KFactor_QCD_ggZZ_QCDScaleUp;
        else if (unc=="UpAs") k_factor = KFactor_QCD_ggZZ_AsUp;
        else if (unc=="UpPDFReplica") k_factor = KFactor_QCD_ggZZ_PDFReplicaUp;
        else if (unc=="DnPDF") k_factor = KFactor_QCD_ggZZ_PDFScaleDn;
        else if (unc=="DnQCD") k_factor = KFactor_QCD_ggZZ_QCDScaleDn;
        else if (unc=="DnAs") k_factor = KFactor_QCD_ggZZ_AsDn;
        else if (unc=="DnPDFReplica") k_factor = KFactor_QCD_ggZZ_PDFReplicaDn;
        else k_factor = KFactor_QCD_ggZZ_Nominal;
    }
    else if ( input_file_name.Contains("ggH")) k_factor = ggH_NNLOPS_weight;
    return k_factor;
}

//===============================================================================
float Datacards::calculate_TauIDSF_OneLeg(short DM, float pt, float eta, int genmatch, int flav, const string& unc)
{
    eta=fabs(eta);
    if (genmatch==6)
        return 1.;
    else if (genmatch==5) {
        if (flav==165) return DeepTauSF_VSjet_ETau->getSFvsPT(pt,genmatch,unc);
        else if (flav==195) return DeepTauSF_VSjet_MuTau->getSFvsPT(pt,genmatch,unc);
        else if (flav==225) return DeepTauSF_VSjet_TauTau->getSFvsPT(pt,genmatch,unc);
        else return 1.;
    }
    else if (genmatch==1 || genmatch==3) {
        if (flav==165) return DeepTauSF_VSe_ETau->getSFvsEta(eta,genmatch,unc);
        else if (flav==195) return DeepTauSF_VSe_MuTau->getSFvsEta(eta,genmatch,unc);
        else if (flav==225) return DeepTauSF_VSe_TauTau->getSFvsEta(eta,genmatch,unc);
        else return 1.;
    }
    else if (genmatch==2 || genmatch==4) {
        if (flav==165) return DeepTauSF_VSmu_ETau->getSFvsEta(eta,genmatch,unc);
        else if (flav==195) return DeepTauSF_VSmu_MuTau->getSFvsEta(eta,genmatch,unc);
        else if (flav==225) return DeepTauSF_VSmu_TauTau->getSFvsEta(eta,genmatch,unc);
        else return 1.;
    }
    else {
        cout<<"The genmatch "<<genmatch<<" is wrong"<<endl;
        return 1.;
    }
}

//===============================================================================
float Datacards::calculate_TauIDSF(TString input_file_name, const string& unc="")
{
    float SF = 1;
    for (size_t ileg=2;ileg<=3;ileg++) {
        if (abs(LepLepId->at(ileg))!=15)
            continue;
        int genmatch;
        if (input_file_name.Contains("VBFH125") || input_file_name.Contains("WplusH125"))
            genmatch=GENMatch(ileg,input_file_name);
        else
            genmatch=TauGenMatch->at(ileg);
        SF=SF*calculate_TauIDSF_OneLeg(TauDecayMode->at(ileg),LepPt->at(ileg),LepEta->at(ileg),genmatch,abs(Z2Flav),unc);
    }
    return SF;
}

//===============================================================================
int Datacards::GENMatch(int ileg,TString input_file_name)
{
    if (abs(LepLepId->at(ileg))!=15) {
        cout<<"Only applicable to pdgId = 15, not "<<ileg<<endl;
        return 6;
    }
    if (input_file_name.Contains("WZTo") || input_file_name.Contains("ZZ") || input_file_name.Contains("H") || input_file_name.Contains("ggTo")) {
        if (TauTES_p_Up->at(ileg)>0) return 5;
        else if (TauFES_p_Up->at(ileg)>0) return 1;
        else return 2;
    }
    else {
        if (TauTES_p_Up->at(ileg)>0) return 5;
        else if (TauFES_p_Up->at(ileg)>0) return 1;
        else return 6;
    }
}