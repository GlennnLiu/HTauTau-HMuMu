// Include classes
#include <HTauTauHMuMu/AnalysisStep/test/SuperRatio/include/LTau.h>

// Constructor
//============================================================
LTau::LTau():Tree()
{    
    _current_process = -999;
    _current_final_state = -999;
    _current_category = -999;
    _current_njet_bin = -999;
    _current_taupt_bin = -999;

    _s_process.clear();
    _s_process.push_back("Htt");
    _s_process.push_back("Hmm");
    _s_process.push_back("Hww");
    _s_process.push_back("diBoson");
    _s_process.push_back("triBoson");
    _s_process.push_back("top");
    _s_process.push_back("fake");
    _s_process.push_back("EWKZee");
    _s_process.push_back("EWKZtt");
    _s_process.push_back("EWKZmm");
    _s_process.push_back("DYee");
    _s_process.push_back("DYtt");
    _s_process.push_back("DYmm");
    _s_process.push_back("TotalMC");
    _s_process.push_back("Data");
    
    // _s_final_state.push_back("mumu");
    _s_final_state.push_back("etau");
    _s_final_state.push_back("mutau");
    // _s_final_state.push_back("tautau");
    // _s_final_state.push_back("emu");
    // _s_final_state.push_back("2l");
    
    _s_categories.push_back("GGH_ptHl10");
    _s_categories.push_back("GGH_ptHg10");
    _s_categories.push_back("VBF_ptHl200");
    _s_categories.push_back("VBF_ptHg200");
    _s_categories.push_back("Boost_1j");
    _s_categories.push_back("Boost_2j");
    _s_categories.push_back("All");

    _s_fake_bkg.push_back("QCD");
    _s_fake_bkg.push_back("Wjet");
    _s_fake_bkg.push_back("TT");
    _s_fake_bkg.push_back("All");

    _s_variables.push_back("Mvis");
    _s_variables.push_back("Mgood");
    
    _lumi=16.8;
    string year="UL2016_postVFP";
}
//--------------------------------------------------------------------------------------

// Destructor
//====================
LTau::~LTau()
{
}
//====================

//===============================================================================
void LTau::SetPaths(string path, string file_name, string savepath)
{
    _path=path;
    _file_name=file_name;
    _savepath=savepath;
}

//===============================================================================
void LTau::SetFileList(string* data, int ndata, string* bkg, int nbkg, string* wj, int nwj, string* tt, int ntt)
{
    for (int i=0;i<ndata;i++) _dataFiles.push_back(data[i]);
    for (int i=0;i<nbkg;i++) _MCBkgFiles.push_back(bkg[i]);
    for (int i=0;i<nwj;i++) _WJFiles.push_back(wj[i]);
    for (int i=0;i<ntt;i++) _TTFiles.push_back(tt[i]);
    
}

//===============================================================================
void LTau::Set_taupt_bin(int n, float* bins)
{
    _n_taupt_bins=n;
    for (int i=0;i<n;i++) _taupt_bins[i]=bins[i];
}

//===============================================================================
void LTau::Set_lpt_bin(int n, float* bins)
{
    _n_lpt_bins=n;
    for (int i=0;i<n;i++) _lpt_bins[i]=bins[i];
}

//===============================================================================
void LTau::Set_Mvis_bin(int n, float* bins)
{
    _n_Mvis_bins=n;
    for (int i=0;i<n;i++) _Mvis_bins[i]=bins[i];
}

//===============================================================================
void LTau::Set_MT_bin(int n, float* bins)
{
    _n_MT_bins=n;
    for (int i=0;i<n;i++) _MT_bins[i]=bins[i];
}

//===============================================================================
void LTau::Set_DR_bin(int n, float* bins)
{
    _n_DR_bins=n;
    for (int i=0;i<n;i++) _DR_bins[i]=bins[i];
}

//===============================================================================
void LTau::Step1_FakeRate_DeclareHistos()
{
    cout<<"[INFO] Declare histograms and functions"<<endl;
    for (int i_bkg=0;i_bkg<num_of_fake_bkg;i_bkg++) {
        for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
            for (int i_njet=0;i_njet<num_of_njet_bins;i_njet++) {
                for (int i=0;i<2;i++) {
                    _histo_name="h_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs)+"_Njet"+std::to_string(i_njet)+(i==1?"_pass":"_fail");
                    _histo_labels=";#tau_{h} p_{T};Weighted events";
                    h[i_bkg][i_fs][i_njet][i]=new TH1F(_histo_name, _histo_labels, _n_taupt_bins-1, _taupt_bins);
                }
                _histo_name="h_FR_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs)+"_Njet"+std::to_string(i_njet);
                _histo_labels=";#tau_{h} p_{T};Weighted events";
                h_FR[i_bkg][i_fs][i_njet]=new TH1F(_histo_name, _histo_labels, _n_taupt_bins-1, _taupt_bins);
                _histo_name="f_FR_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs)+"_Njet"+std::to_string(i_njet);
                f_FR[i_bkg][i_fs][i_njet]=new TF1(_histo_name, "pol1", _taupt_bins[0], _taupt_bins[_n_taupt_bins-1]);
            }
        }
    }
}

//===============================================================================
void LTau::Step1_FakeRate_FillHistos()
{
    cout<<"[INFO] Fill in the histograms"<<endl;
    for (size_t i_proc=0;i_proc<_dataFiles.size();i_proc++) {
        TString input_file_name=_path+_dataFiles[i_proc]+_file_name;
        cout<<"[INFO] Processing "<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");

        for (int i_bkg=0;i_bkg<num_of_fake_bkg-1;i_bkg++) {
            TString tree_name=GetTreeName(i_bkg);

            // hCounters = (TH1F*)input_file->Get(tree_name+"/Counters");
            // gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);
            
            input_tree = (TTree*)input_file->Get(tree_name+"/candTree");
            Init( input_tree, input_file_name, true);

            if (fChain == 0) {return;}

            Long64_t nentries = fChain->GetEntriesFast();
            cout<<"[INFO] Processing "<<input_file_name<<" "<<tree_name<<", total event: "<<nentries<<endl;
            Long64_t nbytes = 0, nb = 0;

            for (Long64_t jentry=0; jentry<nentries;jentry++) {
                Long64_t ientry = LoadTree(jentry);
                if (ientry%50000==0) cout<<ientry<<endl;
                if (ientry < 0) break;
                nb = fChain->GetEntry(jentry);
                nbytes += nb;
                
                if (!pass_Trigger) continue;

                _current_final_state = FindFinalState();
                if (_current_final_state<0) continue;

                // if (nCleanedJetsPt30==0) _current_njet_bin=0;
                // else if (nCleanedJetsPt30==1) _current_njet_bin=1;
                // else _current_njet_bin=2;
                _current_njet_bin=NumberOfJets();

                // _k_factor = calculate_K_factor(input_file_name);

                // _event_weight = (_lumi * 1000 * xsec * _k_factor * PUWeight * genHEPMCweight * dataMCWeight * L1prefiringWeight) / gen_sum_weights;

                if (TauVSjet->at(1)>=5) _pass=1;
                else _pass=0;

                h[i_bkg][_current_final_state][_current_njet_bin][_pass]->Fill(LepPt->at(1));
            }
        }
    }

    for (size_t i_proc=0;i_proc<_MCBkgFiles.size();i_proc++) {
        TString input_file_name=_path+_MCBkgFiles[i_proc]+_file_name;
        cout<<"[INFO] Processing "<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");

        for (int i_bkg=0;i_bkg<num_of_fake_bkg-1;i_bkg++) {
            if (i_bkg==1 && _MCBkgFiles[i_proc]==_WJFiles[0]) continue;

            TString tree_name=GetTreeName(i_bkg);

            hCounters = (TH1F*)input_file->Get(tree_name+"/Counters");
            gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);
            
            input_tree = (TTree*)input_file->Get(tree_name+"/candTree");
            Init( input_tree, input_file_name, true);

            if (fChain == 0) {return;}

            Long64_t nentries = fChain->GetEntriesFast();
            cout<<"[INFO] Processing "<<input_file_name<<" "<<tree_name<<", total event: "<<nentries<<endl;
            Long64_t nbytes = 0, nb = 0;

            for (Long64_t jentry=0; jentry<nentries;jentry++) {
                Long64_t ientry = LoadTree(jentry);
                if (ientry%50000==0) cout<<ientry<<endl;
                if (ientry < 0) break;
                nb = fChain->GetEntry(jentry);
                nbytes += nb;
                
                if (!pass_Trigger) continue;

                _current_final_state = FindFinalState();
                if (_current_final_state<0) continue;

                // if (nCleanedJetsPt30==0) _current_njet_bin=0;
                // else if (nCleanedJetsPt30==1) _current_njet_bin=1;
                // else _current_njet_bin=2;
                _current_njet_bin=NumberOfJets();

                if (h[i_bkg][_current_final_state][_current_njet_bin][_pass]->GetBinContent(h[i_bkg][_current_final_state][_current_njet_bin][_pass]->FindBin(LepPt->at(1)))>0) {
                    _k_factor = calculate_K_factor(input_file_name);

                    _event_weight = (_lumi * 1000 * xsec * _k_factor * PUWeight * genHEPMCweight * dataMCWeight * L1prefiringWeight) / gen_sum_weights;

                    if (TauVSjet->at(1)>=5) _pass=1;
                    else _pass=0;

                    h[i_bkg][_current_final_state][_current_njet_bin][_pass]->Fill(LepPt->at(1),-1.*_event_weight);
                }
            }
        }
    }

    for (size_t i_proc=0;i_proc<_TTFiles.size();i_proc++) {
        TString input_file_name=_path+_TTFiles[i_proc]+_file_name;
        cout<<"[INFO] Processing "<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");

        size_t i_bkg=2;
        TString tree_name=GetTreeName(i_bkg);

        hCounters = (TH1F*)input_file->Get(tree_name+"/Counters");
        gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);
        
        input_tree = (TTree*)input_file->Get(tree_name+"/candTree");
        Init( input_tree, input_file_name, true);

        if (fChain == 0) {return;}

        Long64_t nentries = fChain->GetEntriesFast();
        cout<<"[INFO] Processing "<<input_file_name<<" "<<tree_name<<", total event: "<<nentries<<endl;
        Long64_t nbytes = 0, nb = 0;

        for (Long64_t jentry=0; jentry<nentries;jentry++) {
            Long64_t ientry = LoadTree(jentry);
            if (ientry%50000==0) cout<<ientry<<endl;
            if (ientry < 0) break;
            nb = fChain->GetEntry(jentry);
            nbytes += nb;
            
            if (!pass_Trigger) continue;

            _current_final_state = FindFinalState();
            if (_current_final_state<0) continue;

            // if (nCleanedJetsPt30==0) _current_njet_bin=0;
            // else if (nCleanedJetsPt30==1) _current_njet_bin=1;
            // else _current_njet_bin=2;
            _current_njet_bin=NumberOfJets();

            _k_factor = calculate_K_factor(input_file_name);

            _event_weight = (_lumi * 1000 * xsec * _k_factor * PUWeight * genHEPMCweight * dataMCWeight * L1prefiringWeight) / gen_sum_weights;

            if (TauVSjet->at(1)>=5) _pass=1;
            else _pass=0;

            h[i_bkg][_current_final_state][_current_njet_bin][_pass]->Fill(LepPt->at(1),_event_weight);
        }
    }
}

//===============================================================================
void LTau::Step1_FakeRate_Compute()
{
    cout<<"[INFO] Computing and fit fake rates"<<endl;
    for (int i_bkg=0;i_bkg<num_of_fake_bkg;i_bkg++) {
        for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
            for (int i_njet=0;i_njet<num_of_njet_bins;i_njet++) {
                h_tmp=(TH1F*)h[i_bkg][i_fs][i_njet][1]->Clone();
                h_tmp->Add(h[i_bkg][i_fs][i_njet][0]);
                DoDivision(h_FR[i_bkg][i_fs][i_njet],h[i_bkg][i_fs][i_njet][1],h_tmp,true);
                h_FR[i_bkg][i_fs][i_njet]->Fit(f_FR[i_bkg][i_fs][i_njet],"R");
            }
        }
    }
}

//===============================================================================
void LTau::Step1_FakeRate_SaveHistos()
{
    cout<<"[INFO] Save histograms and functions"<<endl;
    TString output_file_name=_savepath+"LTau.Step1_FakeRate.root";
    output_file = TFile::Open(output_file_name, "recreate");
    output_file->cd();
    for (int i_bkg=0;i_bkg<num_of_fake_bkg;i_bkg++) {
        for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
            for (int i_njet=0;i_njet<num_of_njet_bins;i_njet++) {
                for (int i=0;i<2;i++) {
                    h[i_bkg][i_fs][i_njet][i]->Write();
                }
                h_FR[i_bkg][i_fs][i_njet]->Write();
                f_FR[i_bkg][i_fs][i_njet]->Write();
            }
        }
    }
    output_file->Close();
}

//===============================================================================
void LTau::Step1_FakeRate_GetHistos()
{
    cout<<"[INFO] Get the fake rate files"<<endl;
    TString FR_file_name=_savepath+"LTau.Step1_FakeRate.root";
    TFile *FR_file = TFile::Open(FR_file_name, "read");
    for (int i_bkg=0;i_bkg<num_of_fake_bkg;i_bkg++) {
        for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
            for (int i_njet=0;i_njet<num_of_njet_bins;i_njet++) {
                for (int i=0;i<2;i++) {
                    _histo_name="h_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs)+"_Njet"+std::to_string(i_njet)+(i==1?"_pass":"_fail");
                    h[i_bkg][i_fs][i_njet][i]=(TH1F*)FR_file->Get(_histo_name);
                }
                _histo_name="h_FR_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs)+"_Njet"+std::to_string(i_njet);
                h_FR[i_bkg][i_fs][i_njet]=(TH1F*)FR_file->Get(_histo_name);
                _histo_name="f_FR_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs)+"_Njet"+std::to_string(i_njet);
                f_FR[i_bkg][i_fs][i_njet]=(TF1*)FR_file->Get(_histo_name);
            }
        }
    }
}

//===============================================================================
void LTau::Step2_Closure_DeclareHistos()
{
    cout<<"[INFO] Declare histograms and functions"<<endl;
    for (int i_bkg=0;i_bkg<num_of_fake_bkg;i_bkg++) {
        for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
            for (int i_taupt=0;i_taupt<num_of_taupt_bins;i_taupt++) {
                for (int i=0;i<2;i++) {
                    _histo_name="hClosure_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs)+"_Ntaupt"+std::to_string(i_taupt)+(i==1?"_pass":"_fail");
                    _histo_labels=";#l p_{T};Weighted events";
                    hClosure[i_bkg][i_fs][i_taupt][i]=new TH1F(_histo_name, _histo_labels, _n_lpt_bins-1, _lpt_bins);
                }
                _histo_name="hClosure_FR_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs)+"_Ntaupt"+std::to_string(i_taupt);
                _histo_labels=";#l p_{T};Weighted events";
                hClosure_FR[i_bkg][i_fs][i_taupt]=new TH1F(_histo_name, _histo_labels, _n_lpt_bins-1, _lpt_bins);
                _histo_name="fClosure_FR_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs)+"_Ntaupt"+std::to_string(i_taupt);
                fClosure_FR[i_bkg][i_fs][i_taupt]=new TF1(_histo_name, "pol3", _lpt_bins[1], _lpt_bins[_n_lpt_bins-1]);
            }
        }
    }
}

//===============================================================================
void LTau::Step2_Closure_FillHistos()
{
    cout<<"[INFO] Fill in the histograms and apply fake rates"<<endl;    
    for (size_t i_proc=0;i_proc<_dataFiles.size();i_proc++) {
        TString input_file_name=_path+_dataFiles[i_proc]+_file_name;
        cout<<"[INFO] Processing "<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");

        for (int i_bkg=0;i_bkg<num_of_fake_bkg-1;i_bkg++) {
            TString tree_name=GetTreeName(i_bkg);

            // hCounters = (TH1F*)input_file->Get(tree_name+"/Counters");
            // gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);
            
            input_tree = (TTree*)input_file->Get(tree_name+"/candTree");
            Init( input_tree, input_file_name, true);

            if (fChain == 0) {return;}

            Long64_t nentries = fChain->GetEntriesFast();
            cout<<"[INFO] Processing "<<input_file_name<<" "<<tree_name<<", total event: "<<nentries<<endl;
            Long64_t nbytes = 0, nb = 0;

            for (Long64_t jentry=0; jentry<nentries;jentry++) {
                Long64_t ientry = LoadTree(jentry);
                if (ientry%50000==0) cout<<ientry<<endl;
                if (ientry < 0) break;
                nb = fChain->GetEntry(jentry);
                nbytes += nb;
                
                if (!pass_Trigger) continue;

                _current_final_state = FindFinalState();
                if (_current_final_state<0) continue;

                // if (nCleanedJetsPt30==0) _current_njet_bin=0;
                // else if (nCleanedJetsPt30==1) _current_njet_bin=1;
                // else _current_njet_bin=2;
                _current_njet_bin=NumberOfJets();

                if (LepPt->at(1)<=40) _current_taupt_bin=0;
                else if (LepPt->at(1)<=50) _current_taupt_bin=1;
                else _current_taupt_bin=2;

                // _k_factor = calculate_K_factor(input_file_name);

                // _event_weight = (_lumi * 1000 * xsec * _k_factor * PUWeight * genHEPMCweight * dataMCWeight * L1prefiringWeight) / gen_sum_weights;

                if (TauVSjet->at(1)>=5) {
                    _pass=1;
                    _event_weight=1;
                    }
                else {
                    _pass=0;
                    float FR=f_FR[i_bkg][_current_final_state][_current_njet_bin]->Eval(LepPt->at(1)<=150?LepPt->at(1):150);
                    _event_weight=FR/(1.-FR);
                }

                hClosure[i_bkg][_current_final_state][_current_taupt_bin][_pass]->Fill(LepPt->at(0),_event_weight);
            }
        }
    }

    for (size_t i_proc=0;i_proc<_MCBkgFiles.size();i_proc++) {
        TString input_file_name=_path+_MCBkgFiles[i_proc]+_file_name;
        cout<<"[INFO] Processing "<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");

        for (int i_bkg=0;i_bkg<num_of_fake_bkg-1;i_bkg++) {
            if (i_bkg==1 && _MCBkgFiles[i_proc]==_WJFiles[0]) continue;

            TString tree_name=GetTreeName(i_bkg);

            hCounters = (TH1F*)input_file->Get(tree_name+"/Counters");
            gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);
            
            input_tree = (TTree*)input_file->Get(tree_name+"/candTree");
            Init( input_tree, input_file_name, true);

            if (fChain == 0) {return;}

            Long64_t nentries = fChain->GetEntriesFast();
            cout<<"[INFO] Processing "<<input_file_name<<" "<<tree_name<<", total event: "<<nentries<<endl;
            Long64_t nbytes = 0, nb = 0;

            for (Long64_t jentry=0; jentry<nentries;jentry++) {
                Long64_t ientry = LoadTree(jentry);
                if (ientry%50000==0) cout<<ientry<<endl;
                if (ientry < 0) break;
                nb = fChain->GetEntry(jentry);
                nbytes += nb;
                
                if (!pass_Trigger) continue;

                _current_final_state = FindFinalState();
                if (_current_final_state<0) continue;

                // if (nCleanedJetsPt30==0) _current_njet_bin=0;
                // else if (nCleanedJetsPt30==1) _current_njet_bin=1;
                // else _current_njet_bin=2;
                _current_njet_bin=NumberOfJets();

                if (LepPt->at(1)<=40) _current_taupt_bin=0;
                else if (LepPt->at(1)<=50) _current_taupt_bin=1;
                else _current_taupt_bin=2;

                if (hClosure[i_bkg][_current_final_state][_current_taupt_bin][_pass]->GetBinContent(hClosure[i_bkg][_current_final_state][_current_taupt_bin][_pass]->FindBin(LepPt->at(0)))>0) {
                    _k_factor = calculate_K_factor(input_file_name);

                    _event_weight = (_lumi * 1000 * xsec * _k_factor * PUWeight * genHEPMCweight * dataMCWeight * L1prefiringWeight) / gen_sum_weights;

                    if (TauVSjet->at(1)>=5) {
                        _pass=1;
                        _event_weight*=-1;
                        }
                    else {
                        _pass=0;
                        float FR=f_FR[i_bkg][_current_final_state][_current_njet_bin]->Eval(LepPt->at(1)<=150?LepPt->at(1):150);
                        _event_weight*=-1.*FR/(1.-FR);
                    }
                    hClosure[i_bkg][_current_final_state][_current_taupt_bin][_pass]->Fill(LepPt->at(0),_event_weight);
                }
            }
        }
    }

    for (size_t i_proc=0;i_proc<_TTFiles.size();i_proc++) {
        TString input_file_name=_path+_TTFiles[i_proc]+_file_name;
        cout<<"[INFO] Processing "<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");

        size_t i_bkg=2;
        TString tree_name=GetTreeName(i_bkg);

        hCounters = (TH1F*)input_file->Get(tree_name+"/Counters");
        gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);
        
        input_tree = (TTree*)input_file->Get(tree_name+"/candTree");
        Init( input_tree, input_file_name, true);

        if (fChain == 0) {return;}

        Long64_t nentries = fChain->GetEntriesFast();
        cout<<"[INFO] Processing "<<input_file_name<<" "<<tree_name<<", total event: "<<nentries<<endl;
        Long64_t nbytes = 0, nb = 0;

        for (Long64_t jentry=0; jentry<nentries;jentry++) {
            Long64_t ientry = LoadTree(jentry);
            if (ientry%50000==0) cout<<ientry<<endl;
            if (ientry < 0) break;
            nb = fChain->GetEntry(jentry);
            nbytes += nb;
            
            if (!pass_Trigger) continue;

            _current_final_state = FindFinalState();
            if (_current_final_state<0) continue;

            // if (nCleanedJetsPt30==0) _current_njet_bin=0;
            // else if (nCleanedJetsPt30==1) _current_njet_bin=1;
            // else _current_njet_bin=2;
            _current_njet_bin=NumberOfJets();

            if (LepPt->at(1)<=40) _current_taupt_bin=0;
            else if (LepPt->at(1)<=50) _current_taupt_bin=1;
            else _current_taupt_bin=2;

            _k_factor = calculate_K_factor(input_file_name);

            _event_weight = (_lumi * 1000 * xsec * _k_factor * PUWeight * genHEPMCweight * dataMCWeight * L1prefiringWeight) / gen_sum_weights;

            if (TauVSjet->at(1)>=5) {
                _pass=1;
                }
            else {
                _pass=0;
                float FR=f_FR[i_bkg][_current_final_state][_current_njet_bin]->Eval(LepPt->at(1)<=150?LepPt->at(1):150);
                _event_weight*=FR/(1.-FR);
            }

            hClosure[i_bkg][_current_final_state][_current_taupt_bin][_pass]->Fill(LepPt->at(0),_event_weight);
        }
    }
}

//===============================================================================
void LTau::Step2_Closure_Compute()
{
    cout<<"[INFO] Computing and fit closure"<<endl;
    for (int i_bkg=0;i_bkg<num_of_fake_bkg;i_bkg++) {
        for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
            for (int i_taupt=0;i_taupt<num_of_taupt_bins;i_taupt++) {
                DoDivision(hClosure_FR[i_bkg][i_fs][i_taupt],hClosure[i_bkg][i_fs][i_taupt][1],hClosure[i_bkg][i_fs][i_taupt][0],false);
                hClosure_FR[i_bkg][i_fs][i_taupt]->Fit(fClosure_FR[i_bkg][i_fs][i_taupt],"R");
            }
        }
    }
}

//===============================================================================
void LTau::Step2_Closure_SaveHistos()
{
    cout<<"[INFO] Save histograms and functions"<<endl;
    TString output_file_name=_savepath+"LTau.Step2_Closure.root";
    output_file = TFile::Open(output_file_name, "recreate");
    output_file->cd();
    for (int i_bkg=0;i_bkg<num_of_fake_bkg;i_bkg++) {
        for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
            for (int i_taupt=0;i_taupt<num_of_taupt_bins;i_taupt++) {
                for (int i=0;i<2;i++) {
                    hClosure[i_bkg][i_fs][i_taupt][i]->Write();
                }
                hClosure_FR[i_bkg][i_fs][i_taupt]->Write();
                fClosure_FR[i_bkg][i_fs][i_taupt]->Write();
            }
        }
    }
    output_file->Close();
}

//===============================================================================
void LTau::Step2_Closure_GetHistos()
{
    cout<<"[INFO] Get the closure files"<<endl;
    TString CL_file_name=_savepath+"LTau.Step2_Closure.root";
    TFile *CL_file = TFile::Open(CL_file_name, "read");
    for (int i_bkg=0;i_bkg<num_of_fake_bkg;i_bkg++) {
        for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
            for (int i_taupt=0;i_taupt<num_of_taupt_bins;i_taupt++) {
                for (int i=0;i<2;i++) {
                    _histo_name="hClosure_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs)+"_Ntaupt"+std::to_string(i_taupt)+(i==1?"_pass":"_fail");
                    hClosure[i_bkg][i_fs][i_taupt][i]=(TH1F*)CL_file->Get(_histo_name);
                }
                _histo_name="hClosure_FR_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs)+"_Ntaupt"+std::to_string(i_taupt);
                hClosure_FR[i_bkg][i_fs][i_taupt]=(TH1F*)CL_file->Get(_histo_name);
                _histo_name="fClosure_FR_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs)+"_Ntaupt"+std::to_string(i_taupt);
                fClosure_FR[i_bkg][i_fs][i_taupt]=(TF1*)CL_file->Get(_histo_name);
            }
        }
    }
}

//===============================================================================
void LTau::Step3_QCD_FakeRate_DeclareHistos()
{
    cout<<"[INFO] Declare histograms and functions"<<endl;
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i_njet=0;i_njet<num_of_njet_bins;i_njet++) {
            for (int i=0;i<2;i++) {
                _histo_name="hQCD_"+_s_final_state.at(i_fs)+"_Njet"+std::to_string(i_njet)+(i==1?"_pass":"_fail");
                _histo_labels=";#tau_{h} p_{T};Weighted events";
                hQCD[i_fs][i_njet][i]=new TH1F(_histo_name, _histo_labels, _n_taupt_bins-1, _taupt_bins);
            }
            _histo_name="hQCD_FR_"+_s_final_state.at(i_fs)+"_Njet"+std::to_string(i_njet);
            _histo_labels=";#tau_{h} p_{T};Weighted events";
            hQCD_FR[i_fs][i_njet]=new TH1F(_histo_name, _histo_labels, _n_taupt_bins-1, _taupt_bins);
            _histo_name="fQCD_FR_"+_s_final_state.at(i_fs)+"_Njet"+std::to_string(i_njet);
            fQCD_FR[i_fs][i_njet]=new TF1(_histo_name, "pol1", _taupt_bins[0], _taupt_bins[_n_taupt_bins-1]);
        }
    }
}

//===============================================================================
void LTau::Step3_QCD_FakeRate_FillHistos()
{
    cout<<"[INFO] Fill in the histograms"<<endl;
    for (size_t i_proc=0;i_proc<_dataFiles.size();i_proc++) {
        TString input_file_name=_path+_dataFiles[i_proc]+_file_name;
        cout<<"[INFO] Processing "<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");

        TString tree_name="CRQCDvSRTree";
        
        input_tree = (TTree*)input_file->Get(tree_name+"/candTree");
        Init( input_tree, input_file_name, true);

        if (fChain == 0) {return;}

        Long64_t nentries = fChain->GetEntriesFast();
        cout<<"[INFO] Processing "<<input_file_name<<" "<<tree_name<<", total event: "<<nentries<<endl;
        Long64_t nbytes = 0, nb = 0;

        for (Long64_t jentry=0; jentry<nentries;jentry++) {
            Long64_t ientry = LoadTree(jentry);
            if (ientry%50000==0) cout<<ientry<<endl;
            if (ientry < 0) break;
            nb = fChain->GetEntry(jentry);
            nbytes += nb;
            
            if (!pass_Trigger) continue;

            _current_final_state = FindFinalState();
            if (_current_final_state<0) continue;

            _current_njet_bin=NumberOfJets();

            if (TauVSjet->at(1)>=5) _pass=1;
            else _pass=0;

            hQCD[_current_final_state][_current_njet_bin][_pass]->Fill(LepPt->at(1));
        }
    }

    for (size_t i_proc=0;i_proc<_MCBkgFiles.size();i_proc++) {
        TString input_file_name=_path+_MCBkgFiles[i_proc]+_file_name;
        cout<<"[INFO] Processing "<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");

        TString tree_name="CRQCDvSRTree";

        hCounters = (TH1F*)input_file->Get(tree_name+"/Counters");
        gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);
        
        input_tree = (TTree*)input_file->Get(tree_name+"/candTree");
        Init( input_tree, input_file_name, true);

        if (fChain == 0) {return;}

        Long64_t nentries = fChain->GetEntriesFast();
        cout<<"[INFO] Processing "<<input_file_name<<" "<<tree_name<<", total event: "<<nentries<<endl;
        Long64_t nbytes = 0, nb = 0;

        for (Long64_t jentry=0; jentry<nentries;jentry++) {
            Long64_t ientry = LoadTree(jentry);
            if (ientry%50000==0) cout<<ientry<<endl;
            if (ientry < 0) break;
            nb = fChain->GetEntry(jentry);
            nbytes += nb;
            
            if (!pass_Trigger) continue;

            _current_final_state = FindFinalState();
            if (_current_final_state<0) continue;

            _current_njet_bin=NumberOfJets();

            if (hQCD[_current_final_state][_current_njet_bin][_pass]->GetBinContent(hQCD[_current_final_state][_current_njet_bin][_pass]->FindBin(LepPt->at(1)))>0) {
                _k_factor = calculate_K_factor(input_file_name);

                _event_weight = (_lumi * 1000 * xsec * _k_factor * PUWeight * genHEPMCweight * dataMCWeight * L1prefiringWeight) / gen_sum_weights;

                if (TauVSjet->at(1)>=5) _pass=1;
                else _pass=0;

                hQCD[_current_final_state][_current_njet_bin][_pass]->Fill(LepPt->at(1),-1.*_event_weight);
            }
        }
    }
}

//===============================================================================
void LTau::Step3_QCD_FakeRate_Compute()
{
    cout<<"[INFO] Computing and fit fake rates"<<endl;
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i_njet=0;i_njet<num_of_njet_bins;i_njet++) {
            h_tmp=(TH1F*)hQCD[i_fs][i_njet][1]->Clone();
            h_tmp->Add(hQCD[i_fs][i_njet][0]);
            DoDivision(hQCD_FR[i_fs][i_njet],hQCD[i_fs][i_njet][1],h_tmp,true);
            hQCD_FR[i_fs][i_njet]->Fit(fQCD_FR[i_fs][i_njet],"R");
        }
    }
}

//===============================================================================
void LTau::Step3_QCD_FakeRate_SaveHistos()
{
    cout<<"[INFO] Save histograms and functions"<<endl;
    TString output_file_name=_savepath+"LTau.Step3_QCD_FakeRate.root";
    output_file = TFile::Open(output_file_name, "recreate");
    output_file->cd();
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i_njet=0;i_njet<num_of_njet_bins;i_njet++) {
            for (int i=0;i<2;i++) {
                hQCD[i_fs][i_njet][i]->Write();
            }
            hQCD_FR[i_fs][i_njet]->Write();
            fQCD_FR[i_fs][i_njet]->Write();
        }
    }
    output_file->Close();
}

//===============================================================================
void LTau::Step3_QCD_FakeRate_GetHistos()
{
    cout<<"[INFO] Get the fake rate files"<<endl;
    TString QCDFR_file_name=_savepath+"LTau.Step3_QCD_FakeRate.root";
    TFile *QCDFR_file = TFile::Open(QCDFR_file_name, "read");
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i_njet=0;i_njet<num_of_njet_bins;i_njet++) {
            for (int i=0;i<2;i++) {
                _histo_name="hQCD_"+_s_final_state.at(i_fs)+"_Njet"+std::to_string(i_njet)+(i==1?"_pass":"_fail");
                hQCD[i_fs][i_njet][i]=(TH1F*)QCDFR_file->Get(_histo_name);
            }
            _histo_name="hQCD_FR_"+_s_final_state.at(i_fs)+"_Njet"+std::to_string(i_njet);
            hQCD_FR[i_fs][i_njet]=(TH1F*)QCDFR_file->Get(_histo_name);
            _histo_name="fQCD_FR_"+_s_final_state.at(i_fs)+"_Njet"+std::to_string(i_njet);
            fQCD_FR[i_fs][i_njet]=(TF1*)QCDFR_file->Get(_histo_name);
        }
    }
}

//===============================================================================
void LTau::Step3_QCD_Closure_DeclareHistos()
{
    cout<<"[INFO] Declare histograms and functions"<<endl;
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i_taupt=0;i_taupt<num_of_taupt_bins;i_taupt++) {
            for (int i=0;i<2;i++) {
                _histo_name="hQCDClosure_"+_s_final_state.at(i_fs)+"_Ntaupt"+std::to_string(i_taupt)+(i==1?"_pass":"_fail");
                _histo_labels=";#l p_{T};Weighted events";
                hQCDClosure[i_fs][i_taupt][i]=new TH1F(_histo_name, _histo_labels, _n_lpt_bins-1, _lpt_bins);
            }
            _histo_name="hQCDClosure_FR_"+_s_final_state.at(i_fs)+"_Ntaupt"+std::to_string(i_taupt);
            _histo_labels=";#l p_{T};Weighted events";
            hQCDClosure_FR[i_fs][i_taupt]=new TH1F(_histo_name, _histo_labels, _n_lpt_bins-1, _lpt_bins);
            _histo_name="fQCDClosure_FR_"+_s_final_state.at(i_fs)+"_Ntaupt"+std::to_string(i_taupt);
            fQCDClosure_FR[i_fs][i_taupt]=new TF1(_histo_name, "pol3", _lpt_bins[1], _lpt_bins[_n_lpt_bins-1]);
        }
    }
}

//===============================================================================
void LTau::Step3_QCD_Closure_FillHistos()
{
    cout<<"[INFO] Fill in the histograms and apply fake rates"<<endl;    
    for (size_t i_proc=0;i_proc<_dataFiles.size();i_proc++) {
        TString input_file_name=_path+_dataFiles[i_proc]+_file_name;
        cout<<"[INFO] Processing "<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");

        TString tree_name="CRQCDvSRTree";
        
        input_tree = (TTree*)input_file->Get(tree_name+"/candTree");
        Init( input_tree, input_file_name, true);

        if (fChain == 0) {return;}

        Long64_t nentries = fChain->GetEntriesFast();
        cout<<"[INFO] Processing "<<input_file_name<<" "<<tree_name<<", total event: "<<nentries<<endl;
        Long64_t nbytes = 0, nb = 0;

        for (Long64_t jentry=0; jentry<nentries;jentry++) {
            Long64_t ientry = LoadTree(jentry);
            if (ientry%50000==0) cout<<ientry<<endl;
            if (ientry < 0) break;
            nb = fChain->GetEntry(jentry);
            nbytes += nb;
            
            if (!pass_Trigger) continue;

            _current_final_state = FindFinalState();
            if (_current_final_state<0) continue;

            _current_njet_bin=NumberOfJets();

            if (LepPt->at(1)<=40) _current_taupt_bin=0;
            else if (LepPt->at(1)<=50) _current_taupt_bin=1;
            else _current_taupt_bin=2;

            if (TauVSjet->at(1)>=5) {
                _pass=1;
                _event_weight=1;
                }
            else {
                _pass=0;
                float FR=fQCD_FR[_current_final_state][_current_njet_bin]->Eval(LepPt->at(1)<=150?LepPt->at(1):150);
                _event_weight=FR/(1.-FR);
            }

            hQCDClosure[_current_final_state][_current_taupt_bin][_pass]->Fill(LepPt->at(0),_event_weight);
        }
    }

    for (size_t i_proc=0;i_proc<_MCBkgFiles.size();i_proc++) {
        TString input_file_name=_path+_MCBkgFiles[i_proc]+_file_name;
        cout<<"[INFO] Processing "<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");

        TString tree_name="CRQCDvSRTree";

        hCounters = (TH1F*)input_file->Get(tree_name+"/Counters");
        gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);
        
        input_tree = (TTree*)input_file->Get(tree_name+"/candTree");
        Init( input_tree, input_file_name, true);

        if (fChain == 0) {return;}

        Long64_t nentries = fChain->GetEntriesFast();
        cout<<"[INFO] Processing "<<input_file_name<<" "<<tree_name<<", total event: "<<nentries<<endl;
        Long64_t nbytes = 0, nb = 0;

        for (Long64_t jentry=0; jentry<nentries;jentry++) {
            Long64_t ientry = LoadTree(jentry);
            if (ientry%50000==0) cout<<ientry<<endl;
            if (ientry < 0) break;
            nb = fChain->GetEntry(jentry);
            nbytes += nb;
            
            if (!pass_Trigger) continue;

            _current_final_state = FindFinalState();
            if (_current_final_state<0) continue;

            _current_njet_bin=NumberOfJets();

            if (LepPt->at(1)<=40) _current_taupt_bin=0;
            else if (LepPt->at(1)<=50) _current_taupt_bin=1;
            else _current_taupt_bin=2;

            if (hQCDClosure[_current_final_state][_current_taupt_bin][_pass]->GetBinContent(hQCDClosure[_current_final_state][_current_taupt_bin][_pass]->FindBin(LepPt->at(0)))>0) {
                _k_factor = calculate_K_factor(input_file_name);

                _event_weight = (_lumi * 1000 * xsec * _k_factor * PUWeight * genHEPMCweight * dataMCWeight * L1prefiringWeight) / gen_sum_weights;

                if (TauVSjet->at(1)>=5) {
                    _pass=1;
                    _event_weight*=-1;
                    }
                else {
                    _pass=0;
                    float FR=fQCD_FR[_current_final_state][_current_njet_bin]->Eval(LepPt->at(1)<=150?LepPt->at(1):150);
                    _event_weight*=-1.*FR/(1.-FR);
                }

                hQCDClosure[_current_final_state][_current_taupt_bin][_pass]->Fill(LepPt->at(0),_event_weight);
            }
        }
    }
}

//===============================================================================
void LTau::Step3_QCD_Closure_Compute()
{
    cout<<"[INFO] Computing and fit closure"<<endl;
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i_taupt=0;i_taupt<num_of_taupt_bins;i_taupt++) {
            DoDivision(hQCDClosure_FR[i_fs][i_taupt],hQCDClosure[i_fs][i_taupt][1],hQCDClosure[i_fs][i_taupt][0],false);
            hQCDClosure_FR[i_fs][i_taupt]->Fit(fQCDClosure_FR[i_fs][i_taupt],"R");
        }
    }
}

//===============================================================================
void LTau::Step3_QCD_Closure_SaveHistos()
{
    cout<<"[INFO] Save histograms and functions"<<endl;
    TString output_file_name=_savepath+"LTau.Step3_QCD_Closure.root";
    output_file = TFile::Open(output_file_name, "recreate");
    output_file->cd();
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i_taupt=0;i_taupt<num_of_taupt_bins;i_taupt++) {
            for (int i=0;i<2;i++) {
                hQCDClosure[i_fs][i_taupt][i]->Write();
            }
            hQCDClosure_FR[i_fs][i_taupt]->Write();
            fQCDClosure_FR[i_fs][i_taupt]->Write();
        }
    }
    output_file->Close();
}

//===============================================================================
void LTau::Step3_QCD_Closure_GetHistos()
{
    cout<<"[INFO] Get the closure files"<<endl;
    TString QCDCL_file_name=_savepath+"LTau.Step3_QCD_Closure.root";
    TFile *QCDCL_file = TFile::Open(QCDCL_file_name, "read");
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i_taupt=0;i_taupt<num_of_taupt_bins;i_taupt++) {
            for (int i=0;i<2;i++) {
                _histo_name="hQCDClosure_"+_s_final_state.at(i_fs)+"_Ntaupt"+std::to_string(i_taupt)+(i==1?"_pass":"_fail");
                hQCDClosure[i_fs][i_taupt][i]=(TH1F*)QCDCL_file->Get(_histo_name);
            }
            _histo_name="hQCDClosure_FR_"+_s_final_state.at(i_fs)+"_Ntaupt"+std::to_string(i_taupt);
            hQCDClosure_FR[i_fs][i_taupt]=(TH1F*)QCDCL_file->Get(_histo_name);
            _histo_name="fQCDClosure_FR_"+_s_final_state.at(i_fs)+"_Ntaupt"+std::to_string(i_taupt);
            fQCDClosure_FR[i_fs][i_taupt]=(TF1*)QCDCL_file->Get(_histo_name);
        }
    }
}

//===============================================================================
void LTau::Step3_QCD_vsSR_DeclareHistos()
{
    cout<<"[INFO] Declare histograms and functions"<<endl;
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i=0;i<2;i++) {
            _histo_name="hQCDvsSR_"+_s_final_state.at(i_fs)+(i==1?"_pass":"_fail");
            _histo_labels=";#Delta R(l,#tau_{h});Weighted events";
            hQCDvsSR[i_fs][i]=new TH1F(_histo_name, _histo_labels, _n_DR_bins-1, _DR_bins);
        }
        _histo_name="hQCDvsSR_FR_"+_s_final_state.at(i_fs);
        _histo_labels=";#Delta R(l,#tau_{h});Weighted events";
        hQCDvsSR_FR[i_fs]=new TH1F(_histo_name, _histo_labels, _n_DR_bins-1, _DR_bins);
        _histo_name="fQCDvsSR_FR_"+_s_final_state.at(i_fs);
        fQCDvsSR_FR[i_fs]=new TF1(_histo_name, "pol2", _DR_bins[0], _DR_bins[_n_DR_bins-1]);
    }
}

//===============================================================================
void LTau::Step3_QCD_vsSR_FillHistos()
{
    cout<<"[INFO] Fill in the histograms and apply fake rates"<<endl;    
    for (size_t i_proc=0;i_proc<_dataFiles.size();i_proc++) {
        TString input_file_name=_path+_dataFiles[i_proc]+_file_name;
        cout<<"[INFO] Processing "<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");

        for (_pass=0;_pass<2;_pass++) {

            TString tree_name=_pass==0?"CRQCDvSROSTree":"CRQCDvSROS1Tree";
            
            input_tree = (TTree*)input_file->Get(tree_name+"/candTree");
            if (!input_tree) {return;}

            Init( input_tree, input_file_name, true);
            if (fChain == 0) {return;}

            Long64_t nentries = fChain->GetEntriesFast();
            cout<<"[INFO] Processing "<<input_file_name<<" "<<tree_name<<", total event: "<<nentries<<endl;
            Long64_t nbytes = 0, nb = 0;

            for (Long64_t jentry=0; jentry<nentries;jentry++) {
                Long64_t ientry = LoadTree(jentry);
                if (ientry%50000==0) cout<<ientry<<endl;
                if (ientry < 0) break;
                nb = fChain->GetEntry(jentry);
                nbytes += nb;
                
                if (!pass_Trigger) continue;

                float DR=std::sqrt((LepEta->at(0)-LepEta->at(1))*(LepEta->at(0)-LepEta->at(1))+(LepPhi->at(0)-LepPhi->at(1))*(LepPhi->at(0)-LepPhi->at(1)));

                _current_final_state = FindFinalState();
                if (_current_final_state<0) continue;

                _current_njet_bin=NumberOfJets();

                if (LepPt->at(1)<=40) _current_taupt_bin=0;
                else if (LepPt->at(1)<=50) _current_taupt_bin=1;
                else _current_taupt_bin=2;

                if (_pass==0) {
                    float FR=fQCD_FR[_current_final_state][_current_njet_bin]->Eval(LepPt->at(1)<=150?LepPt->at(1):150);
                    float CL;
                    if (LepPt->at(0)<=23 && hQCDClosure_FR[_current_final_state][_current_taupt_bin]->GetBinContent(1)>0) 
                        CL=hQCDClosure_FR[_current_final_state][_current_taupt_bin]->GetBinContent(1);
                    else
                        CL=fQCDClosure_FR[_current_final_state][_current_taupt_bin]->Eval(LepPt->at(0)<=150?LepPt->at(0):150);
                    _event_weight=FR/(1.-FR)*CL;
                }
                else {
                    _event_weight=1;
                }

                hQCDvsSR[_current_final_state][_pass]->Fill(DR,_event_weight);
            }
        }
    }

    for (size_t i_proc=0;i_proc<_MCBkgFiles.size();i_proc++) {
        for (_pass=0;_pass<2;_pass++) {
            TString input_file_name=_path+_MCBkgFiles[i_proc]+_file_name;//(_pass==0?_file_name:"HTauTauHMuMu_left.root");
            cout<<"[INFO] Processing "<<input_file_name<<endl;
            input_file = TFile::Open(input_file_name,"read");

            TString tree_name=_pass==0?"CRQCDvSROSTree":"CRQCDvSROS1Tree";

            input_tree = (TTree*)input_file->Get(tree_name+"/candTree");
            if (!input_tree) {
                input_file_name=_path+_MCBkgFiles[i_proc]+"/HTauTauHMuMu_left.root";
                input_file = TFile::Open(input_file_name,"read");
                input_tree = (TTree*)input_file->Get(tree_name+"/candTree");
            }
            
            Init( input_tree, input_file_name, true);
            if (fChain == 0) {return;}

            hCounters = (TH1F*)input_file->Get(tree_name+"/Counters");
            gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);

            Long64_t nentries = fChain->GetEntriesFast();
            cout<<"[INFO] Processing "<<input_file_name<<" "<<tree_name<<", total event: "<<nentries<<endl;
            Long64_t nbytes = 0, nb = 0;

            for (Long64_t jentry=0; jentry<nentries;jentry++) {
                Long64_t ientry = LoadTree(jentry);
                if (ientry%50000==0) cout<<ientry<<endl;
                if (ientry < 0) break;
                nb = fChain->GetEntry(jentry);
                nbytes += nb;

                if (!pass_Trigger) continue;

                float DR=std::sqrt((LepEta->at(0)-LepEta->at(1))*(LepEta->at(0)-LepEta->at(1))+(LepPhi->at(0)-LepPhi->at(1))*(LepPhi->at(0)-LepPhi->at(1)));

                _current_final_state = FindFinalState();
                if (_current_final_state<0) continue;

                _current_njet_bin=NumberOfJets();

                if (LepPt->at(1)<=40) _current_taupt_bin=0;
                else if (LepPt->at(1)<=50) _current_taupt_bin=1;
                else _current_taupt_bin=2;


                if (hQCDvsSR[_current_final_state][_pass]->GetBinContent(hQCDvsSR[_current_final_state][_pass]->FindBin(DR))>0) {
                    _k_factor = calculate_K_factor(input_file_name);

                    _event_weight = (_lumi * 1000 * xsec * _k_factor * PUWeight * genHEPMCweight * dataMCWeight * L1prefiringWeight) / gen_sum_weights;

                    if (_pass==0) {
                        float FR=fQCD_FR[_current_final_state][_current_njet_bin]->Eval(LepPt->at(1)<=150?LepPt->at(1):150);
                        float CL;
                        if (LepPt->at(0)<=23 && hQCDClosure_FR[_current_final_state][_current_taupt_bin]->GetBinContent(1)>0) 
                            CL=hQCDClosure_FR[_current_final_state][_current_taupt_bin]->GetBinContent(1);
                        else
                            CL=fQCDClosure_FR[_current_final_state][_current_taupt_bin]->Eval(LepPt->at(0)<=150?LepPt->at(0):150);
                        _event_weight*=-1.*FR/(1.-FR)*CL;
                    }
                    else {
                        _event_weight*=-1.;
                    }

                    hQCDvsSR[_current_final_state][_pass]->Fill(DR,_event_weight);
                }
            }
        }
    }
}

//===============================================================================
void LTau::Step3_QCD_vsSR_Compute()
{
    cout<<"[INFO] Computing and fit closure"<<endl;
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        DoDivision(hQCDvsSR_FR[i_fs],hQCDvsSR[i_fs][1],hQCDvsSR[i_fs][0],false);
        hQCDvsSR_FR[i_fs]->Fit(fQCDvsSR_FR[i_fs],"R");
    }
}

//===============================================================================
void LTau::Step3_QCD_vsSR_SaveHistos()
{
    cout<<"[INFO] Save histograms and functions"<<endl;
    TString output_file_name=_savepath+"LTau.Step3_QCD_vsSR.root";
    output_file = TFile::Open(output_file_name, "recreate");
    output_file->cd();
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i=0;i<2;i++) {
            hQCDvsSR[i_fs][i]->Write();
        }
        hQCDvsSR_FR[i_fs]->Write();
        fQCDvsSR_FR[i_fs]->Write();
    }
    output_file->Close();
}

//===============================================================================
void LTau::Step3_QCD_vsSR_GetHistos()
{
    cout<<"[INFO] Get the closure files"<<endl;
    TString QCDvsSR_file_name=_savepath+"LTau.Step3_QCD_vsSR.root";
    TFile *QCDvsSR_file = TFile::Open(QCDvsSR_file_name, "read");
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i=0;i<2;i++) {
            _histo_name="hQCDvsSR_"+_s_final_state.at(i_fs)+(i==1?"_pass":"_fail");
            hQCDvsSR[i_fs][i]=(TH1F*)QCDvsSR_file->Get(_histo_name);
        }
        _histo_name="hQCDvsSR_FR_"+_s_final_state.at(i_fs);
        hQCDvsSR_FR[i_fs]=(TH1F*)QCDvsSR_file->Get(_histo_name);
        _histo_name="fQCDvsSR_FR_"+_s_final_state.at(i_fs);
        fQCDvsSR_FR[i_fs]=(TF1*)QCDvsSR_file->Get(_histo_name);
    }
}

//===============================================================================
void LTau::Step3_WJ_FakeRate_DeclareHistos()
{
    cout<<"[INFO] Declare histograms and functions"<<endl;
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i_njet=0;i_njet<num_of_njet_bins;i_njet++) {
            for (int i=0;i<2;i++) {
                _histo_name="hWJ_"+_s_final_state.at(i_fs)+"_Njet"+std::to_string(i_njet)+(i==1?"_pass":"_fail");
                _histo_labels=";#tau_{h} p_{T};Weighted events";
                hWJ[i_fs][i_njet][i]=new TH1F(_histo_name, _histo_labels, _n_taupt_bins-1, _taupt_bins);
            }
            _histo_name="hWJ_FR_"+_s_final_state.at(i_fs)+"_Njet"+std::to_string(i_njet);
            _histo_labels=";#tau_{h} p_{T};Weighted events";
            hWJ_FR[i_fs][i_njet]=new TH1F(_histo_name, _histo_labels, _n_taupt_bins-1, _taupt_bins);
            _histo_name="fWJ_FR_"+_s_final_state.at(i_fs)+"_Njet"+std::to_string(i_njet);
            fWJ_FR[i_fs][i_njet]=new TF1(_histo_name, "pol1", _taupt_bins[0], _taupt_bins[_n_taupt_bins-1]);
        }
    }
}

//===============================================================================
void LTau::Step3_WJ_FakeRate_FillHistos()
{
    cout<<"[INFO] Fill in the histograms"<<endl;

    for (size_t i_proc=0;i_proc<_WJFiles.size();i_proc++) {
        TString input_file_name=_path+_WJFiles[i_proc]+_file_name;
        cout<<"[INFO] Processing "<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");

        TString tree_name="CRWJvSRTree";

        hCounters = (TH1F*)input_file->Get(tree_name+"/Counters");
        gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);
        
        input_tree = (TTree*)input_file->Get(tree_name+"/candTree");
        Init( input_tree, input_file_name, true);

        if (fChain == 0) {return;}

        Long64_t nentries = fChain->GetEntriesFast();
        cout<<"[INFO] Processing "<<input_file_name<<" "<<tree_name<<", total event: "<<nentries<<endl;
        Long64_t nbytes = 0, nb = 0;

        for (Long64_t jentry=0; jentry<nentries;jentry++) {
            Long64_t ientry = LoadTree(jentry);
            if (ientry%50000==0) cout<<ientry<<endl;
            if (ientry < 0) break;
            nb = fChain->GetEntry(jentry);
            nbytes += nb;
            
            if (!pass_Trigger) continue;
            if (MtLMET<=70) continue;

            _current_final_state = FindFinalState();
            if (_current_final_state<0) continue;

            _current_njet_bin=NumberOfJets();

            _k_factor = calculate_K_factor(input_file_name);

            _event_weight = (_lumi * 1000 * xsec * _k_factor * PUWeight * genHEPMCweight * dataMCWeight * L1prefiringWeight) / gen_sum_weights;

            if (TauVSjet->at(1)>=5) _pass=1;
            else _pass=0;

            hWJ[_current_final_state][_current_njet_bin][_pass]->Fill(LepPt->at(1),_event_weight);
        }
    }
}

//===============================================================================
void LTau::Step3_WJ_FakeRate_Compute()
{
    cout<<"[INFO] Computing and fit fake rates"<<endl;
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i_njet=0;i_njet<num_of_njet_bins;i_njet++) {
            h_tmp=(TH1F*)hWJ[i_fs][i_njet][1]->Clone();
            h_tmp->Add(hWJ[i_fs][i_njet][0]);
            DoDivision(hWJ_FR[i_fs][i_njet],hWJ[i_fs][i_njet][1],h_tmp,true);
            hWJ_FR[i_fs][i_njet]->Fit(fWJ_FR[i_fs][i_njet],"R");
        }
    }
}

//===============================================================================
void LTau::Step3_WJ_FakeRate_SaveHistos()
{
    cout<<"[INFO] Save histograms and functions"<<endl;
    TString output_file_name=_savepath+"LTau.Step3_WJ_FakeRate.root";
    output_file = TFile::Open(output_file_name, "recreate");
    output_file->cd();
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i_njet=0;i_njet<num_of_njet_bins;i_njet++) {
            for (int i=0;i<2;i++) {
                hWJ[i_fs][i_njet][i]->Write();
            }
            hWJ_FR[i_fs][i_njet]->Write();
            fWJ_FR[i_fs][i_njet]->Write();
        }
    }
    output_file->Close();
}

//===============================================================================
void LTau::Step3_WJ_FakeRate_GetHistos()
{
    cout<<"[INFO] Get the fake rate files"<<endl;
    TString WJFR_file_name=_savepath+"LTau.Step3_WJ_FakeRate.root";
    TFile *WJFR_file = TFile::Open(WJFR_file_name, "read");
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i_njet=0;i_njet<num_of_njet_bins;i_njet++) {
            for (int i=0;i<2;i++) {
                _histo_name="hWJ_"+_s_final_state.at(i_fs)+"_Njet"+std::to_string(i_njet)+(i==1?"_pass":"_fail");
                hWJ[i_fs][i_njet][i]=(TH1F*)WJFR_file->Get(_histo_name);
            }
            _histo_name="hWJ_FR_"+_s_final_state.at(i_fs)+"_Njet"+std::to_string(i_njet);
            hWJ_FR[i_fs][i_njet]=(TH1F*)WJFR_file->Get(_histo_name);
            _histo_name="fWJ_FR_"+_s_final_state.at(i_fs)+"_Njet"+std::to_string(i_njet);
            fWJ_FR[i_fs][i_njet]=(TF1*)WJFR_file->Get(_histo_name);
        }
    }
}

//===============================================================================
void LTau::Step3_WJ_Closure_DeclareHistos()
{
    cout<<"[INFO] Declare histograms and functions"<<endl;
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i_taupt=0;i_taupt<num_of_taupt_bins;i_taupt++) {
            for (int i=0;i<2;i++) {
                _histo_name="hWJClosure_"+_s_final_state.at(i_fs)+"_Ntaupt"+std::to_string(i_taupt)+(i==1?"_pass":"_fail");
                _histo_labels=";#l p_{T};Weighted events";
                hWJClosure[i_fs][i_taupt][i]=new TH1F(_histo_name, _histo_labels, _n_lpt_bins-1, _lpt_bins);
            }
            _histo_name="hWJClosure_FR_"+_s_final_state.at(i_fs)+"_Ntaupt"+std::to_string(i_taupt);
            _histo_labels=";#l p_{T};Weighted events";
            hWJClosure_FR[i_fs][i_taupt]=new TH1F(_histo_name, _histo_labels, _n_lpt_bins-1, _lpt_bins);
            _histo_name="fWJClosure_FR_"+_s_final_state.at(i_fs)+"_Ntaupt"+std::to_string(i_taupt);
            fWJClosure_FR[i_fs][i_taupt]=new TF1(_histo_name, "pol3", _lpt_bins[1], _lpt_bins[_n_lpt_bins-1]);
        }
    }
}

//===============================================================================
void LTau::Step3_WJ_Closure_FillHistos()
{
    cout<<"[INFO] Fill in the histograms and apply fake rates"<<endl;    

    for (size_t i_proc=0;i_proc<_WJFiles.size();i_proc++) {
        TString input_file_name=_path+_WJFiles[i_proc]+_file_name;
        cout<<"[INFO] Processing "<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");

        TString tree_name="CRWJvSRTree";

        hCounters = (TH1F*)input_file->Get(tree_name+"/Counters");
        gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);
        
        input_tree = (TTree*)input_file->Get(tree_name+"/candTree");
        Init( input_tree, input_file_name, true);

        if (fChain == 0) {return;}

        Long64_t nentries = fChain->GetEntriesFast();
        cout<<"[INFO] Processing "<<input_file_name<<" "<<tree_name<<", total event: "<<nentries<<endl;
        Long64_t nbytes = 0, nb = 0;

        for (Long64_t jentry=0; jentry<nentries;jentry++) {
            Long64_t ientry = LoadTree(jentry);
            if (ientry%50000==0) cout<<ientry<<endl;
            if (ientry < 0) break;
            nb = fChain->GetEntry(jentry);
            nbytes += nb;
            
            if (!pass_Trigger) continue;
            if (MtLMET<=70) continue;

            _current_final_state = FindFinalState();
            if (_current_final_state<0) continue;

            _current_njet_bin=NumberOfJets();

            if (LepPt->at(1)<=40) _current_taupt_bin=0;
            else if (LepPt->at(1)<=50) _current_taupt_bin=1;
            else _current_taupt_bin=2;

            _k_factor = calculate_K_factor(input_file_name);

            _event_weight = (_lumi * 1000 * xsec * _k_factor * PUWeight * genHEPMCweight * dataMCWeight * L1prefiringWeight) / gen_sum_weights;

            if (TauVSjet->at(1)>=5) {
                _pass=1;
                }
            else {
                _pass=0;
                float FR=fWJ_FR[_current_final_state][_current_njet_bin]->Eval(LepPt->at(1)<=150?LepPt->at(1):150);
                _event_weight*=FR/(1.-FR);
            }

            hWJClosure[_current_final_state][_current_taupt_bin][_pass]->Fill(LepPt->at(0),_event_weight);
        }
    }
}

//===============================================================================
void LTau::Step3_WJ_Closure_Compute()
{
    cout<<"[INFO] Computing and fit closure"<<endl;
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i_taupt=0;i_taupt<num_of_taupt_bins;i_taupt++) {
            DoDivision(hWJClosure_FR[i_fs][i_taupt],hWJClosure[i_fs][i_taupt][1],hWJClosure[i_fs][i_taupt][0],false);
            hWJClosure_FR[i_fs][i_taupt]->Fit(fWJClosure_FR[i_fs][i_taupt],"R");
        }
    }
}

//===============================================================================
void LTau::Step3_WJ_Closure_SaveHistos()
{
    cout<<"[INFO] Save histograms and functions"<<endl;
    TString output_file_name=_savepath+"LTau.Step3_WJ_Closure.root";
    output_file = TFile::Open(output_file_name, "recreate");
    output_file->cd();
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i_taupt=0;i_taupt<num_of_taupt_bins;i_taupt++) {
            for (int i=0;i<2;i++) {
                hWJClosure[i_fs][i_taupt][i]->Write();
            }
            hWJClosure_FR[i_fs][i_taupt]->Write();
            fWJClosure_FR[i_fs][i_taupt]->Write();
        }
    }
    output_file->Close();
}

//===============================================================================
void LTau::Step3_WJ_Closure_GetHistos()
{
    cout<<"[INFO] Get the closure files"<<endl;
    TString WJCL_file_name=_savepath+"LTau.Step3_WJ_Closure.root";
    TFile *WJCL_file = TFile::Open(WJCL_file_name, "read");
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i_taupt=0;i_taupt<num_of_taupt_bins;i_taupt++) {
            for (int i=0;i<2;i++) {
                _histo_name="hWJClosure_"+_s_final_state.at(i_fs)+"_Ntaupt"+std::to_string(i_taupt)+(i==1?"_pass":"_fail");
                hWJClosure[i_fs][i_taupt][i]=(TH1F*)WJCL_file->Get(_histo_name);
            }
            _histo_name="hWJClosure_FR_"+_s_final_state.at(i_fs)+"_Ntaupt"+std::to_string(i_taupt);
            hWJClosure_FR[i_fs][i_taupt]=(TH1F*)WJCL_file->Get(_histo_name);
            _histo_name="fWJClosure_FR_"+_s_final_state.at(i_fs)+"_Ntaupt"+std::to_string(i_taupt);
            fWJClosure_FR[i_fs][i_taupt]=(TF1*)WJCL_file->Get(_histo_name);
        }
    }
}

//===============================================================================
void LTau::Step3_WJ_vsSR_DeclareHistos()
{
    cout<<"[INFO] Declare histograms and functions"<<endl;
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i=0;i<2;i++) {
            _histo_name="hWJvsSR_"+_s_final_state.at(i_fs)+(i==1?"_pass":"_fail");
            _histo_labels=";#M_{T}(l,MET);Weighted events";
            hWJvsSR[i_fs][i]=new TH1F(_histo_name, _histo_labels, _n_MT_bins-1, _MT_bins);
        }
        _histo_name="hWJvsSR_FR_"+_s_final_state.at(i_fs);
        _histo_labels=";#M_{T}(l,MET);Weighted events";
        hWJvsSR_FR[i_fs]=new TH1F(_histo_name, _histo_labels, _n_MT_bins-1, _MT_bins);
        _histo_name="fWJvsSR_FR_"+_s_final_state.at(i_fs);
        fWJvsSR_FR[i_fs]=new TF1(_histo_name, "pol1", _MT_bins[0], _MT_bins[_n_MT_bins-1]);
    }
}

//===============================================================================
void LTau::Step3_WJ_vsSR_FillHistos()
{
    cout<<"[INFO] Fill in the histograms and apply fake rates"<<endl;
    for (size_t i_proc=0;i_proc<_WJFiles.size();i_proc++) {
        TString input_file_name=_path+_WJFiles[i_proc]+_file_name;
        cout<<"[INFO] Processing "<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");

        TString tree_name="CRWJvSRTree";

        input_tree = (TTree*)input_file->Get(tree_name+"/candTree");
        
        Init( input_tree, input_file_name, true);
        if (fChain == 0) {return;}

        hCounters = (TH1F*)input_file->Get(tree_name+"/Counters");
        gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);

        Long64_t nentries = fChain->GetEntriesFast();
        cout<<"[INFO] Processing "<<input_file_name<<" "<<tree_name<<", total event: "<<nentries<<endl;
        Long64_t nbytes = 0, nb = 0;

        for (Long64_t jentry=0; jentry<nentries;jentry++) {
            Long64_t ientry = LoadTree(jentry);
            if (ientry%50000==0) cout<<ientry<<endl;
            if (ientry < 0) break;
            nb = fChain->GetEntry(jentry);
            nbytes += nb;

            if (!pass_Trigger) continue;

            _current_final_state = FindFinalState();
            if (_current_final_state<0) continue;

            _current_njet_bin=NumberOfJets();

            if (LepPt->at(1)<=40) _current_taupt_bin=0;
            else if (LepPt->at(1)<=50) _current_taupt_bin=1;
            else _current_taupt_bin=2;

            _k_factor = calculate_K_factor(input_file_name);

            _event_weight = (_lumi * 1000 * xsec * _k_factor * PUWeight * genHEPMCweight * dataMCWeight * L1prefiringWeight) / gen_sum_weights;

            if (TauVSjet->at(1)>=5) {
                _pass=1;
                }
            else {
                _pass=0;
                float FR=fWJ_FR[_current_final_state][_current_njet_bin]->Eval(LepPt->at(1)<=150?LepPt->at(1):150);
                float CL;
                if (LepPt->at(0)<=23 && hWJClosure_FR[_current_final_state][_current_taupt_bin]->GetBinContent(1)>0)
                    CL=hWJClosure_FR[_current_final_state][_current_taupt_bin]->GetBinContent(1);
                else
                    CL=fWJClosure_FR[_current_final_state][_current_taupt_bin]->Eval(LepPt->at(0)<=150?LepPt->at(0):150);
                _event_weight*=FR/(1.-FR)*CL;
            }

            hWJvsSR[_current_final_state][_pass]->Fill(MtLMET,_event_weight);
        }
    }
}

//===============================================================================
void LTau::Step3_WJ_vsSR_Compute()
{
    cout<<"[INFO] Computing and fit closure"<<endl;
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        DoDivision(hWJvsSR_FR[i_fs],hWJvsSR[i_fs][1],hWJvsSR[i_fs][0],false);
        hWJvsSR_FR[i_fs]->Fit(fWJvsSR_FR[i_fs],"R");
    }
}

void LTau::Step3_WJ_vsSR_SaveHistos()
{
    cout<<"[INFO] Save histograms and functions"<<endl;
    TString output_file_name=_savepath+"LTau.Step3_WJ_vsSR.root";
    output_file = TFile::Open(output_file_name, "recreate");
    output_file->cd();
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i=0;i<2;i++) {
            hWJvsSR[i_fs][i]->Write();
        }
        hWJvsSR_FR[i_fs]->Write();
        fWJvsSR_FR[i_fs]->Write();
    }
    output_file->Close();
}

//===============================================================================
void LTau::Step3_WJ_vsSR_GetHistos()
{
    cout<<"[INFO] Get the closure files"<<endl;
    TString WJvsSR_file_name=_savepath+"LTau.Step3_WJ_vsSR.root";
    TFile *WJvsSR_file = TFile::Open(WJvsSR_file_name, "read");
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i=0;i<2;i++) {
            _histo_name="hWJvsSR_"+_s_final_state.at(i_fs)+(i==1?"_pass":"_fail");
            hWJvsSR[i_fs][i]=(TH1F*)WJvsSR_file->Get(_histo_name);
        }
        _histo_name="hWJvsSR_FR_"+_s_final_state.at(i_fs);
        hWJvsSR_FR[i_fs]=(TH1F*)WJvsSR_file->Get(_histo_name);
        _histo_name="fWJvsSR_FR_"+_s_final_state.at(i_fs);
        fWJvsSR_FR[i_fs]=(TF1*)WJvsSR_file->Get(_histo_name);
    }
}

//===============================================================================
void LTau::Step4_Fraction_DeclareHistos()
{
    cout<<"[INFO] Declare histograms and functions"<<endl;
    for (int i_bkg=0;i_bkg<num_of_fake_bkg+1;i_bkg++) {
        for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
            for (int i_cat=0;i_cat<num_of_categories;i_cat++) {
                _histo_name="hFrac_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs)+"_"+_s_categories.at(i_cat);
                _histo_labels=";M_{vis};Weighted events";
                hFrac[i_bkg][i_fs][i_cat]=new TH1F(_histo_name, _histo_labels, _n_Mvis_bins-1, _Mvis_bins);
                _histo_name="hFrac_FR_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs)+"_"+_s_categories.at(i_cat);
                hFrac_FR[i_bkg][i_fs][i_cat]=new TH1F(_histo_name, _histo_labels, _n_Mvis_bins-1, _Mvis_bins);
            }
        }
    }
}

//===============================================================================
void LTau::Step4_Fraction_FillHistos()
{
    cout<<"[INFO] Fill in the histograms"<<endl;
    for (size_t i_proc=0;i_proc<_dataFiles.size();i_proc++) {
        TString input_file_name=_path+_dataFiles[i_proc]+_file_name;
        cout<<"[INFO] Processing "<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");

        TString tree_name="CRAPPOSTree";

        input_tree = (TTree*)input_file->Get(tree_name+"/candTree");
        Init( input_tree, input_file_name, true);

        if (fChain == 0) {return;}

        Long64_t nentries = fChain->GetEntriesFast();
        cout<<"[INFO] Processing "<<input_file_name<<" "<<tree_name<<", total event: "<<nentries<<endl;
        Long64_t nbytes = 0, nb = 0;

        for (Long64_t jentry=0; jentry<nentries;jentry++) {
            Long64_t ientry = LoadTree(jentry);
            if (ientry%50000==0) cout<<ientry<<endl;
            if (ientry < 0) break;
            nb = fChain->GetEntry(jentry);
            nbytes += nb;
            
            if (!pass_Trigger) continue;

            _current_final_state = FindFinalState();
            if (_current_final_state<0) continue;

            _current_category = FindCategory();
            if (_current_category<0) continue;

            hFrac[3][_current_final_state][_current_category]->Fill(LLMass);
        }
    }

    for (size_t i_proc=0;i_proc<_MCBkgFiles.size();i_proc++) {

        if (std::find(_WJFiles.begin(),_WJFiles.end(),_MCBkgFiles[i_proc])!=_WJFiles.end()) {
            _current_process=1;
            }
        else if (std::find(_TTFiles.begin(),_TTFiles.end(),_MCBkgFiles[i_proc])!=_TTFiles.end()) {
            _current_process=2;
        }
        else {
            _current_process=3;
        }

        TString input_file_name=_path+_MCBkgFiles[i_proc]+_file_name;
        cout<<"[INFO] Processing "<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");

        TString tree_name="CRAPPOSTree";

        hCounters = (TH1F*)input_file->Get(tree_name+"/Counters");
        gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);
        
        input_tree = (TTree*)input_file->Get(tree_name+"/candTree");
        Init( input_tree, input_file_name, true);

        if (fChain == 0) {return;}

        Long64_t nentries = fChain->GetEntriesFast();
        cout<<"[INFO] Processing "<<input_file_name<<" "<<tree_name<<", total event: "<<nentries<<endl;
        Long64_t nbytes = 0, nb = 0;

        for (Long64_t jentry=0; jentry<nentries;jentry++) {
            Long64_t ientry = LoadTree(jentry);
            if (ientry%50000==0) cout<<ientry<<endl;
            if (ientry < 0) break;
            nb = fChain->GetEntry(jentry);
            nbytes += nb;
            
            if (!pass_Trigger) continue;

            _current_final_state = FindFinalState();
            if (_current_final_state<0) continue;

            _current_category = FindCategory();
            if (_current_category<0) continue;

            if (_current_process==3) {
                if (hFrac[3][_current_final_state][_current_category]->GetBinContent(hFrac[3][_current_final_state][_current_category]->FindBin(LLMass))>0) {
                    _k_factor = calculate_K_factor(input_file_name);

                    _event_weight = (_lumi * 1000 * xsec * _k_factor * PUWeight * genHEPMCweight * dataMCWeight * L1prefiringWeight) / gen_sum_weights;

                    hFrac[3][_current_final_state][_current_category]->Fill(LLMass,-1.*_event_weight);
                }
            }
            else {
                _k_factor = calculate_K_factor(input_file_name);

                _event_weight = (_lumi * 1000 * xsec * _k_factor * PUWeight * genHEPMCweight * dataMCWeight * L1prefiringWeight) / gen_sum_weights;

                hFrac[_current_process][_current_final_state][_current_category]->Fill(LLMass,_event_weight);
            }
        }
    }
}

//===============================================================================
void LTau::Step4_Fraction_Compute()
{
    cout<<"[INFO] Computing fractions"<<endl;
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i_cat=0;i_cat<num_of_categories;i_cat++) {
            hFrac[0][i_fs][i_cat]->Add(hFrac[3][i_fs][i_cat]);
            hFrac[0][i_fs][i_cat]->Add(hFrac[1][i_fs][i_cat],-1.);
            hFrac[0][i_fs][i_cat]->Add(hFrac[2][i_fs][i_cat],-1.);
            float integral=hFrac[0][i_fs][i_cat]->Integral();
            for (int ibin=1;ibin<=hFrac[0][i_fs][i_cat]->GetNbinsX();ibin++) {
                if (hFrac[0][i_fs][i_cat]->GetBinContent(ibin)<0) {
                    hFrac[0][i_fs][i_cat]->SetBinContent(ibin,0);
                }
            }
            if (integral>0) 
                hFrac[0][i_fs][i_cat]->Scale(integral/hFrac[0][i_fs][i_cat]->Integral());
            hFrac[3][i_fs][i_cat]->Reset();
            for (int i_bkg=0;i_bkg<3;i_bkg++) hFrac[3][i_fs][i_cat]->Add(hFrac[i_bkg][i_fs][i_cat]);
        }
    }
    for (int i_bkg=0;i_bkg<num_of_fake_bkg+1;i_bkg++) {
        for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
            for (int i_cat=0;i_cat<num_of_categories;i_cat++) {
                DoDivision(hFrac_FR[i_bkg][i_fs][i_cat],hFrac[i_bkg][i_fs][i_cat],hFrac[3][i_fs][i_cat],true);
            }
        }
    }
}

//===============================================================================
void LTau::Step4_Fraction_SaveHistos()
{
    cout<<"[INFO] Save histograms and functions"<<endl;
    TString output_file_name=_savepath+"LTau.Step4_Fraction.root";
    output_file = TFile::Open(output_file_name, "recreate");
    output_file->cd();
    for (int i_bkg=0;i_bkg<num_of_fake_bkg+1;i_bkg++) {
        for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
            for (int i_cat=0;i_cat<num_of_categories;i_cat++) {
                hFrac[i_bkg][i_fs][i_cat]->Write();
                hFrac_FR[i_bkg][i_fs][i_cat]->Write();
            }
        }
    }
    output_file->Close();
}

//===============================================================================
void LTau::Step4_Fraction_GetHistos()
{
    cout<<"[INFO] Get the fake rate files"<<endl;
    TString Frac_file_name=_savepath+"LTau.Step4_Fraction.root";
    TFile *Frac_file = TFile::Open(Frac_file_name, "read");
    for (int i_bkg=0;i_bkg<num_of_fake_bkg+1;i_bkg++) {
        for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
            for (int i_cat=0;i_cat<num_of_categories;i_cat++) {
                _histo_name="hFrac_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs)+"_"+_s_categories.at(i_cat);
                hFrac[i_bkg][i_fs][i_cat]=(TH1F*)Frac_file->Get(_histo_name);
                _histo_name="hFrac_FR_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs)+"_"+_s_categories.at(i_cat);
                hFrac_FR[i_bkg][i_fs][i_cat]=(TH1F*)Frac_file->Get(_histo_name);
            }
        }
    }
}

//===============================================================================
int LTau::FindFinalState()
{
    int final_state = -999;
    if (abs(LLFlav)==165) final_state=0;
    else if (abs(LLFlav)==195) final_state=1;
    else final_state=-1;
    return final_state;
}

//===============================================================================
int LTau::FindCategory()
{
    int category=-1;
    int njet = NumberOfJets();
    if (njet == 0) {
        if (LLPt<=10) category = 0;
        else category = 1;
    }
    else if (VBFJetIdx1>=0) {
        if (LLPt<=200) category = 2;
        else category = 3;
    }
    else {
        if (njet == 1) category = 4;
        else category = 5;
    }
    return category;
}

//===============================================================================
TString LTau::GetTreeName(int i_bkg)
{
    if (i_bkg==0) return "CRQCDTree";
    else if (i_bkg==1) return "CRWJTree";
    else return "CRTTTree";
}

//===============================================================================
float LTau::calculate_K_factor(TString input_file_name)
{
    float k_factor = 1;
    if ( input_file_name.Contains("ZZTo")) k_factor = KFactor_EW_qqZZ * KFactor_QCD_qqZZ_M;
    else if ( input_file_name.Contains("ggTo")) k_factor = KFactor_QCD_ggZZ_Nominal;
    else if ( input_file_name.Contains("ggH")) k_factor = ggH_NNLOPS_weight;
    return k_factor;
}

//===============================================================================
void LTau::SetColor(TH1F *h,int color)
{
    h->SetMarkerColor(color);
    h->SetLineColor(color);
    h->SetFillColor(color);
}

//===============================================================================
TString LTau::ToFSName(string name)
{
    TString n=name;
    n.ReplaceAll("mu","#mu");
    n.ReplaceAll("tautau","#tau_{h}#tau_{h}");
    n.ReplaceAll("etau","#tau_{e}#tau_{h}");
    n.ReplaceAll("#mutau","#tau_{#mu}#tau_{h}");
    n.ReplaceAll("emu","#tau_{e}#tau_{#mu}");
    return n;
}

//===============================================================================
TString LTau::ToProcessName(string name)
{
    TString n=name;
    n.ReplaceAll("mu","#mu");
    n.ReplaceAll("tau","#tau");
    n.ReplaceAll("H","H#rightarrow");
    n.ReplaceAll("DY","DYZ#rightarrow");
    n.ReplaceAll("EWKZ","EWKZ#rightarrow");
    return n;
}

//===============================================================================
void LTau::DoDivision(TH1F* h, TH1F* num, TH1F* den, bool corr)
{
    // h->Divide(h0);
    for (int i=0;i<h->GetNbinsX();i++) {
        double x=num->GetBinContent(i+1);
        double x0=den->GetBinContent(i+1);
        double e=num->GetBinError(i+1);
        double e0=den->GetBinError(i+1);
        cout<<"x:"<<x<<",e:"<<e<<",x0:"<<x0<<",e0:"<<e0<<endl;
        if (x0<=0) continue;
        h->SetBinContent(i+1,x/x0);
        if (corr) h->SetBinError(i+1,std::sqrt(e*e*(x0-x)*(x0-x)+(e0*e0-e*e)*x*x)/x0/x0);
        else h->SetBinError(i+1,std::sqrt(e*e*x0*x0+e0*e0*x*x)/x0/x0);
    }
}

//===============================================================================
int LTau::NumberOfJets()
{
    int N=0;
    float l1eta=LepEta->at(0);
    float l1phi=LepPhi->at(0);
    float l2eta=LepEta->at(1);
    float l2phi=LepPhi->at(1);
    for (size_t ijet=0;ijet<JetPt->size();ijet++) {
        if (JetPt->at(ijet)<30) continue;
        float jeta=JetEta->at(ijet);
        float jphi=JetPhi->at(ijet);
        if ((jeta-l1eta)*(jeta-l1eta)+(jphi-l1phi)*(jphi-l1phi)<0.25) continue;
        if ((jeta-l2eta)*(jeta-l2eta)+(jphi-l2phi)*(jphi-l2phi)<0.25) continue;
        N++;
        if (N==2) break;
    }
    return N;
}