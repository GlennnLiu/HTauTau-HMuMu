// Include classes
#include <HTauTauHMuMu/AnalysisStep/test/SuperRatio/include/LL.h>

// Constructor
//============================================================
LL::LL():Tree()
{    
    _current_process = -999;
    _current_final_state = -999;
    _current_category = -999;
    _current_njet_bin = -999;

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
    _s_final_state.push_back("emu");
    // _s_final_state.push_back("mumu");
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
LL::~LL()
{
}
//====================

//===============================================================================
void LL::SetPaths(string path, string file_name, string savepath)
{
    _path=path;
    _file_name=file_name;
    _savepath=savepath;
}

//===============================================================================
void LL::SetFileList(string* data, int ndata, string* bkg, int nbkg, string* tt, int ntt)
{
    for (int i=0;i<ndata;i++) _dataFiles.push_back(data[i]);
    for (int i=0;i<nbkg;i++) _MCBkgFiles.push_back(bkg[i]);
    for (int i=0;i<ntt;i++) _TTFiles.push_back(tt[i]);
}

//===============================================================================
void LL::Set_lpt_bin(int n, float* bins, int i)
{
    if (i==1) {
        _n_l1pt_bins=n;
        for (int i=0;i<n;i++) _l1pt_bins[i]=bins[i];
    }
    else if (i==2) {
        _n_l2pt_bins=n;
        for (int i=0;i<n;i++) _l2pt_bins[i]=bins[i];
    }
    else {
        cout<<"[ERROR] Wrong lepton index: "<<i<<"!!!"<<endl;
    }
}

//===============================================================================
void LL::Set_Mvis_bin(int n, float* bins)
{
    _n_Mvis_bins=n;
    for (int i=0;i<n;i++) _Mvis_bins[i]=bins[i];
}

//===============================================================================
void LL::Set_DR_bin(int n, float* bins)
{
    _n_DR_bins=n;
    for (int i=0;i<n;i++) _DR_bins[i]=bins[i];
}

//===============================================================================
void LL::Step1_FakeRate_DeclareHistos()
{
    cout<<"[INFO] Declare histograms and functions"<<endl;
    for (int i_bkg=0;i_bkg<num_of_fake_bkg;i_bkg++) {
        for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
            // for (int i_njet=0;i_njet<num_of_njet_bins;i_njet++) {
            for (int i=0;i<2;i++) {
                _histo_name="h_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs)+(i==1?"_pass":"_fail");
                _histo_labels=";#Delta R(l_{1},l_{2});Weighted events";
                h[i_bkg][i_fs][i]=new TH1F(_histo_name, _histo_labels, _n_DR_bins-1, _DR_bins);
            }
            _histo_name="h_FR_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs);
            _histo_labels=";#Delta R(l_{1},l_{2});Weighted events";
            h_FR[i_bkg][i_fs]=new TH1F(_histo_name, _histo_labels, _n_DR_bins-1, _DR_bins);
            _histo_name="f_FR_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs);
            f_FR[i_bkg][i_fs]=new TF1(_histo_name, i_bkg==0?"pol1":"pol4", _DR_bins[0], _DR_bins[_n_DR_bins-1]);
            // }
        }
    }
}

//===============================================================================
void LL::Step1_FakeRate_FillHistos()
{
    cout<<"[INFO] Fill in the histograms"<<endl;
    for (size_t i_proc=0;i_proc<_dataFiles.size();i_proc++) {
        TString input_file_name=_path+_dataFiles[i_proc]+_file_name;
        cout<<"[INFO] Processing "<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");

        int i_bkg=0;

        for (_pass=0;_pass<2;_pass++) {
            TString tree_name=(_pass==0?"CRAPPSSTree":"CRAPPOSTree");
            
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

                // _current_njet_bin=nCleanedJetsPt30>2?2:nCleanedJetsPt30;//NumberOfJets();

                float DR=std::sqrt((LepEta->at(0)-LepEta->at(1))*(LepEta->at(0)-LepEta->at(1))+(LepPhi->at(0)-LepPhi->at(1))*(LepPhi->at(0)-LepPhi->at(1)));

                h[i_bkg][_current_final_state][_pass]->Fill(DR);
            }
        }
    } 

    for (size_t i_proc=0;i_proc<_MCBkgFiles.size();i_proc++) {
        TString input_file_name=_path+_MCBkgFiles[i_proc]+_file_name;
        cout<<"[INFO] Processing "<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");
        
        int i_bkg=0;

        for (_pass=0;_pass<2;_pass++) {
            TString tree_name=(_pass==0?"CRAPPSSTree":"CRAPPOSTree");

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

                // _current_njet_bin=nCleanedJetsPt30>2?2:nCleanedJetsPt30;//NumberOfJets();

                float DR=std::sqrt((LepEta->at(0)-LepEta->at(1))*(LepEta->at(0)-LepEta->at(1))+(LepPhi->at(0)-LepPhi->at(1))*(LepPhi->at(0)-LepPhi->at(1)));

                // if (h[i_bkg][_current_final_state][_pass]->GetBinContent(h[i_bkg][_current_final_state][_pass]->FindBin(DR))>0) {
                _k_factor = calculate_K_factor(input_file_name);

                _event_weight = (_lumi * 1000 * xsec * _k_factor * PUWeight * genHEPMCweight * dataMCWeight * L1prefiringWeight) / gen_sum_weights;

                h[i_bkg][_current_final_state][_pass]->Fill(DR,-1.*_event_weight);
                // }
            }
        }
    }

    for (size_t i_proc=0;i_proc<_TTFiles.size();i_proc++) {
        TString input_file_name=_path+_TTFiles[i_proc]+_file_name;
        cout<<"[INFO] Processing "<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");
        
        int i_bkg=1;
        for (_pass=0;_pass<2;_pass++) {
            TString tree_name=(_pass==0?"CRQCDTree":"SRTree");

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
                if (LepCombRelIsoPF->at(1)>0.15) continue;
                if (nCleanedJetsPt25BTagged_bTagSF>=1) continue;

                _current_final_state = FindFinalState();
                if (_current_final_state<0) continue;

                float DR=std::sqrt((LepEta->at(0)-LepEta->at(1))*(LepEta->at(0)-LepEta->at(1))+(LepPhi->at(0)-LepPhi->at(1))*(LepPhi->at(0)-LepPhi->at(1)));

                // _current_njet_bin=nCleanedJetsPt30>2?2:nCleanedJetsPt30;//NumberOfJets();

                _k_factor = calculate_K_factor(input_file_name);

                _event_weight = (_lumi * 1000 * xsec * _k_factor * PUWeight * genHEPMCweight * dataMCWeight * L1prefiringWeight) / gen_sum_weights;

                h[i_bkg][_current_final_state][_pass]->Fill(DR,_event_weight);
            }
        }
    }
}

//===============================================================================
void LL::Step1_FakeRate_Compute()
{
    cout<<"[INFO] Computing and fit fake rates"<<endl;
    for (int i_bkg=0;i_bkg<num_of_fake_bkg;i_bkg++) {
        for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
            // for (int i_njet=0;i_njet<num_of_njet_bins;i_njet++) {
            DoDivision(h_FR[i_bkg][i_fs],h[i_bkg][i_fs][1],h[i_bkg][i_fs][0],false);
            h_FR[i_bkg][i_fs]->Fit(f_FR[i_bkg][i_fs],"R");
            // }
        }
    }
}

//===============================================================================
void LL::Step1_FakeRate_SaveHistos()
{
    cout<<"[INFO] Save histograms and functions"<<endl;
    TString output_file_name=_savepath+"LL.Step1_FakeRate.root";
    output_file = TFile::Open(output_file_name, "recreate");
    output_file->cd();
    for (int i_bkg=0;i_bkg<num_of_fake_bkg;i_bkg++) {
        for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
            // for (int i_njet=0;i_njet<num_of_njet_bins;i_njet++) {
            for (int i=0;i<2;i++) {
                h[i_bkg][i_fs][i]->Write();
            }
            h_FR[i_bkg][i_fs]->Write();
            f_FR[i_bkg][i_fs]->Write();
            // }
        }
    }
    output_file->Close();
}

//===============================================================================
void LL::Step1_FakeRate_GetHistos()
{
    cout<<"[INFO] Get the fake rate files"<<endl;
    TString FR_file_name=_savepath+"LL.Step1_FakeRate.root";
    TFile *FR_file = TFile::Open(FR_file_name, "read");
    for (int i_bkg=0;i_bkg<num_of_fake_bkg;i_bkg++) {
        for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
            // for (int i_njet=0;i_njet<num_of_njet_bins;i_njet++) {
            for (int i=0;i<2;i++) {
                _histo_name="h_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs)+(i==1?"_pass":"_fail");
                h[i_bkg][i_fs][i]=(TH1F*)FR_file->Get(_histo_name);
            }
            _histo_name="h_FR_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs);
            h_FR[i_bkg][i_fs]=(TH1F*)FR_file->Get(_histo_name);
            _histo_name="f_FR_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs);
            f_FR[i_bkg][i_fs]=(TF1*)FR_file->Get(_histo_name);
            // }
        }
    }
}

//===============================================================================
void LL::Step2_Closure_DeclareHistos()
{
    cout<<"[INFO] Declare histograms and functions"<<endl;
    for (int i_bkg=0;i_bkg<num_of_fake_bkg;i_bkg++) {
        for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
            for (int i=0;i<2;i++) {
                _histo_name="hClosure_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs)+(i==1?"_pass":"_fail");
                _histo_labels=";#l_{1} p_{T};l_{2} p_{T}";
                hClosure[i_bkg][i_fs][i]=new TH2F(_histo_name, _histo_labels, _n_l1pt_bins-1, _l1pt_bins,_n_l2pt_bins-1, _l2pt_bins);
            }
            _histo_name="hClosure_FR_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs);
            _histo_labels=";#l_{1} p_{T};l_{2} p_{T}";
            hClosure_FR[i_bkg][i_fs]=new TH2F(_histo_name, _histo_labels, _n_l1pt_bins-1, _l1pt_bins,_n_l2pt_bins-1, _l2pt_bins);
        }
    }
}

//===============================================================================
void LL::Step2_Closure_FillHistos()
{
    cout<<"[INFO] Fill in the histograms and apply fake rates"<<endl;    
    for (size_t i_proc=0;i_proc<_dataFiles.size();i_proc++) {
        TString input_file_name=_path+_dataFiles[i_proc]+_file_name;
        cout<<"[INFO] Processing "<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");

        int i_bkg=0;

        for (_pass=0;_pass<2;_pass++) {
            TString tree_name=(_pass==0?"CRAPPSSTree":"CRAPPOSTree");
            
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

                // _current_njet_bin=NumberOfJets();

                float DR=std::sqrt((LepEta->at(0)-LepEta->at(1))*(LepEta->at(0)-LepEta->at(1))+(LepPhi->at(0)-LepPhi->at(1))*(LepPhi->at(0)-LepPhi->at(1)));
                
                _event_weight=1;
                if (_pass==0) {
                    float FR=f_FR[i_bkg][_current_final_state]->Eval(DR>6?6:DR);
                    _event_weight=FR;
                }

                hClosure[i_bkg][_current_final_state][_pass]->Fill(LepPt->at(0),LepPt->at(1),_event_weight);
            }
        }
    }

    for (size_t i_proc=0;i_proc<_MCBkgFiles.size();i_proc++) {
        TString input_file_name=_path+_MCBkgFiles[i_proc]+_file_name;
        cout<<"[INFO] Processing "<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");

        int i_bkg=0;

        for (_pass=0;_pass<2;_pass++) {
            TString tree_name=(_pass==0?"CRAPPSSTree":"CRAPPOSTree");

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

                // _current_njet_bin=NumberOfJets();

                float DR=std::sqrt((LepEta->at(0)-LepEta->at(1))*(LepEta->at(0)-LepEta->at(1))+(LepPhi->at(0)-LepPhi->at(1))*(LepPhi->at(0)-LepPhi->at(1)));

                // if (hClosure[i_bkg][_current_final_state][_pass]->GetBinContent(hClosure[i_bkg][_current_final_state][_pass]->FindBin(LepPt->at(0),LepPt->at(1)))>0) {
                _k_factor = calculate_K_factor(input_file_name);

                _event_weight = (_lumi * 1000 * xsec * _k_factor * PUWeight * genHEPMCweight * dataMCWeight * L1prefiringWeight) / gen_sum_weights;

                if (_pass==0) {
                    float FR=f_FR[i_bkg][_current_final_state]->Eval(DR>6?6:DR);
                    _event_weight*=FR;
                }

                hClosure[i_bkg][_current_final_state][_pass]->Fill(LepPt->at(0),LepPt->at(1),-1.*_event_weight);
                // }
            }
        }
    }

    for (size_t i_proc=0;i_proc<_TTFiles.size();i_proc++) {
        TString input_file_name=_path+_TTFiles[i_proc]+_file_name;
        cout<<"[INFO] Processing "<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");

        if (std::find(_TTFiles.begin(),_TTFiles.end(),_MCBkgFiles[i_proc])!=_TTFiles.end()) {
            continue;
        }

        int i_bkg=1;

        for (_pass=0;_pass<2;_pass++) {
            TString tree_name=(_pass==0?"CRQCDTree":"SRTree");

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
                if (LepCombRelIsoPF->at(1)>0.15) continue;
                if (nCleanedJetsPt25BTagged_bTagSF>=1) continue;

                _current_final_state = FindFinalState();
                if (_current_final_state<0) continue;

                // _current_njet_bin=NumberOfJets();

                float DR=std::sqrt((LepEta->at(0)-LepEta->at(1))*(LepEta->at(0)-LepEta->at(1))+(LepPhi->at(0)-LepPhi->at(1))*(LepPhi->at(0)-LepPhi->at(1)));

                _k_factor = calculate_K_factor(input_file_name);

                _event_weight = (_lumi * 1000 * xsec * _k_factor * PUWeight * genHEPMCweight * dataMCWeight * L1prefiringWeight) / gen_sum_weights;

                if (_pass==0) {
                    float FR=f_FR[i_bkg][_current_final_state]->Eval(DR>6?6:DR);
                    _event_weight*=FR;
                }

                hClosure[i_bkg][_current_final_state][_pass]->Fill(LepPt->at(0),LepPt->at(1),_event_weight);
            }
        }
    }
}

//===============================================================================
void LL::Step2_Closure_Compute()
{
    cout<<"[INFO] Computing and fit closure"<<endl;
    for (int i_bkg=0;i_bkg<num_of_fake_bkg;i_bkg++) {
        for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
            DoDivision(hClosure_FR[i_bkg][i_fs],hClosure[i_bkg][i_fs][1],hClosure[i_bkg][i_fs][0],false);
        }
    }
}

//===============================================================================
void LL::Step2_Closure_SaveHistos()
{
    cout<<"[INFO] Save histograms and functions"<<endl;
    TString output_file_name=_savepath+"LL.Step2_Closure.root";
    output_file = TFile::Open(output_file_name, "recreate");
    output_file->cd();
    for (int i_bkg=0;i_bkg<num_of_fake_bkg;i_bkg++) {
        for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
            for (int i=0;i<2;i++) {
                hClosure[i_bkg][i_fs][i]->Write();
            }
            hClosure_FR[i_bkg][i_fs]->Write();
        }
    }
    output_file->Close();
}

//===============================================================================
void LL::Step2_Closure_GetHistos()
{
    cout<<"[INFO] Get the closure files"<<endl;
    TString CL_file_name=_savepath+"LL.Step2_Closure.root";
    TFile *CL_file = TFile::Open(CL_file_name, "read");
    for (int i_bkg=0;i_bkg<num_of_fake_bkg;i_bkg++) {
        for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
            for (int i=0;i<2;i++) {
                _histo_name="hClosure_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs)+(i==1?"_pass":"_fail");
                hClosure[i_bkg][i_fs][i]=(TH2F*)CL_file->Get(_histo_name);
            }
            _histo_name="hClosure_FR_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs);
            hClosure_FR[i_bkg][i_fs]=(TH2F*)CL_file->Get(_histo_name);
        }
    }
}

//===============================================================================
void LL::Step3_QCD_FakeRate_DeclareHistos()
{
    cout<<"[INFO] Declare histograms and functions"<<endl;
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        // for (int i_njet=0;i_njet<num_of_njet_bins;i_njet++) {
        for (int i=0;i<2;i++) {
            _histo_name="hQCD_"+_s_final_state.at(i_fs)+(i==1?"_pass":"_fail");
            _histo_labels=";#Delta R(l_{1},l_{2});Weighted events";
            hQCD[i_fs][i]=new TH1F(_histo_name, _histo_labels, _n_DR_bins-1, _DR_bins);
        }
        _histo_name="hQCD_FR_"+_s_final_state.at(i_fs);
        _histo_labels=";#Delta R(l_{1},l_{2});Weighted events";
        hQCD_FR[i_fs]=new TH1F(_histo_name, _histo_labels, _n_DR_bins-1, _DR_bins);
        _histo_name="fQCD_FR_"+_s_final_state.at(i_fs);
        fQCD_FR[i_fs]=new TF1(_histo_name, "pol1", _DR_bins[0], _DR_bins[_n_DR_bins-1]);
        // }
    }
}

//===============================================================================
void LL::Step3_QCD_FakeRate_FillHistos()
{
    cout<<"[INFO] Fill in the histograms"<<endl;
    for (size_t i_proc=0;i_proc<_dataFiles.size();i_proc++) {
        TString input_file_name=_path+_dataFiles[i_proc]+_file_name;
        cout<<"[INFO] Processing "<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");

        for (_pass=0;_pass<2;_pass++) {
            TString tree_name=(_pass==0?"CRQCDvSRSSTree":"CRQCDvSROSTree");
            
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

                // _current_njet_bin=NumberOfJets();

                float DR=std::sqrt((LepEta->at(0)-LepEta->at(1))*(LepEta->at(0)-LepEta->at(1))+(LepPhi->at(0)-LepPhi->at(1))*(LepPhi->at(0)-LepPhi->at(1)));

                hQCD[_current_final_state][_pass]->Fill(DR);
            }
        }
    }

    for (size_t i_proc=0;i_proc<_MCBkgFiles.size();i_proc++) {
        TString input_file_name=_path+_MCBkgFiles[i_proc]+_file_name;
        cout<<"[INFO] Processing "<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");
        
        for (_pass=0;_pass<2;_pass++) {
            TString tree_name=(_pass==0?"CRQCDvSRSSTree":"CRQCDvSROSTree");

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

                // _current_njet_bin=NumberOfJets();

                float DR=std::sqrt((LepEta->at(0)-LepEta->at(1))*(LepEta->at(0)-LepEta->at(1))+(LepPhi->at(0)-LepPhi->at(1))*(LepPhi->at(0)-LepPhi->at(1)));

                // if (hQCD[_current_final_state][_pass]->GetBinContent(hQCD[_current_final_state][_pass]->FindBin(DR)>0)) {
                _k_factor = calculate_K_factor(input_file_name);

                _event_weight = (_lumi * 1000 * xsec * _k_factor * PUWeight * genHEPMCweight * dataMCWeight * L1prefiringWeight) / gen_sum_weights;

                hQCD[_current_final_state][_pass]->Fill(DR,-1.*_event_weight);
                // }
            }
        }
    }
}

//===============================================================================
void LL::Step3_QCD_FakeRate_Compute()
{
    cout<<"[INFO] Computing and fit fake rates"<<endl;
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        // for (int i_njet=0;i_njet<num_of_njet_bins;i_njet++) {
            DoDivision(hQCD_FR[i_fs],hQCD[i_fs][1],hQCD[i_fs][0],false);
            hQCD_FR[i_fs]->Fit(fQCD_FR[i_fs],"R");
        // }
    }
}

//===============================================================================
void LL::Step3_QCD_FakeRate_SaveHistos()
{
    cout<<"[INFO] Save histograms and functions"<<endl;
    TString output_file_name=_savepath+"LL.Step3_QCD_FakeRate.root";
    output_file = TFile::Open(output_file_name, "recreate");
    output_file->cd();
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        // for (int i_njet=0;i_njet<num_of_njet_bins;i_njet++) {
        for (int i=0;i<2;i++) {
            hQCD[i_fs][i]->Write();
        }
        hQCD_FR[i_fs]->Write();
        fQCD_FR[i_fs]->Write();
        // }
    }
    output_file->Close();
}

//===============================================================================
void LL::Step3_QCD_FakeRate_GetHistos()
{
    cout<<"[INFO] Get the fake rate files"<<endl;
    TString FR_file_name=_savepath+"LL.Step3_QCD_FakeRate.root";
    TFile *FR_file = TFile::Open(FR_file_name, "read");
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        // for (int i_njet=0;i_njet<num_of_njet_bins;i_njet++) {
        for (int i=0;i<2;i++) {
            _histo_name="hQCD_"+_s_final_state.at(i_fs)+(i==1?"_pass":"_fail");
            hQCD[i_fs][i]=(TH1F*)FR_file->Get(_histo_name);
        }
        _histo_name="hQCD_FR_"+_s_final_state.at(i_fs);
        hQCD_FR[i_fs]=(TH1F*)FR_file->Get(_histo_name);
        _histo_name="fQCD_FR_"+_s_final_state.at(i_fs);
        fQCD_FR[i_fs]=(TF1*)FR_file->Get(_histo_name);
        // }
    }
}

//===============================================================================
void LL::Step3_QCD_Closure_DeclareHistos()
{
    cout<<"[INFO] Declare histograms and functions"<<endl;
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i=0;i<2;i++) {
            _histo_name="hQCDClosure_"+_s_final_state.at(i_fs)+(i==1?"_pass":"_fail");
            _histo_labels=";#l_{1} p_{T};l_{2} p_{T}";
            hQCDClosure[i_fs][i]=new TH2F(_histo_name, _histo_labels, _n_l1pt_bins-1, _l1pt_bins,_n_l2pt_bins-1, _l2pt_bins);
        }
        _histo_name="hQCDClosure_FR_"+_s_final_state.at(i_fs);
        _histo_labels=";#l_{1} p_{T};l_{2} p_{T}";
        hQCDClosure_FR[i_fs]=new TH2F(_histo_name, _histo_labels, _n_l1pt_bins-1, _l1pt_bins,_n_l2pt_bins-1, _l2pt_bins);
    }
}

//===============================================================================
void LL::Step3_QCD_Closure_FillHistos()
{
    cout<<"[INFO] Fill in the histograms and apply fake rates"<<endl;    
    for (size_t i_proc=0;i_proc<_dataFiles.size();i_proc++) {
        TString input_file_name=_path+_dataFiles[i_proc]+_file_name;
        cout<<"[INFO] Processing "<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");

        for (_pass=0;_pass<2;_pass++) {
            TString tree_name=(_pass==0?"CRQCDvSRSSTree":"CRQCDvSROSTree");
            
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

                float DR=std::sqrt((LepEta->at(0)-LepEta->at(1))*(LepEta->at(0)-LepEta->at(1))+(LepPhi->at(0)-LepPhi->at(1))*(LepPhi->at(0)-LepPhi->at(1)));

                if (_pass==0) {
                    float FR=fQCD_FR[_current_final_state]->Eval(DR>6?6:DR);
                    _event_weight=FR;
                }
                else {
                    _event_weight=1;
                }

                hQCDClosure[_current_final_state][_pass]->Fill(LepPt->at(0),LepPt->at(1),_event_weight);
            }
        }
    }

    for (size_t i_proc=0;i_proc<_MCBkgFiles.size();i_proc++) {
        TString input_file_name=_path+_MCBkgFiles[i_proc]+_file_name;
        cout<<"[INFO] Processing "<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");

        for (_pass=0;_pass<2;_pass++) {
            TString tree_name=(_pass==0?"CRQCDvSRSSTree":"CRQCDvSROSTree");

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

                float DR=std::sqrt((LepEta->at(0)-LepEta->at(1))*(LepEta->at(0)-LepEta->at(1))+(LepPhi->at(0)-LepPhi->at(1))*(LepPhi->at(0)-LepPhi->at(1)));

                // if (hQCDClosure[_current_final_state][_pass]->GetBinContent(hQCDClosure[_current_final_state][_pass]->FindBin(LepPt->at(0),LepPt->at(1)))>0) {
                _k_factor = calculate_K_factor(input_file_name);

                _event_weight = (_lumi * 1000 * xsec * _k_factor * PUWeight * genHEPMCweight * dataMCWeight * L1prefiringWeight) / gen_sum_weights;

                if (_pass==0) {
                    float FR=fQCD_FR[_current_final_state]->Eval(DR>6?6:DR);
                    _event_weight*=-1.*FR;
                }
                else {
                    _event_weight*=-1.;
                }

                hQCDClosure[_current_final_state][_pass]->Fill(LepPt->at(0),LepPt->at(1),_event_weight);
                // }
            }
        }
    }
}

//===============================================================================
void LL::Step3_QCD_Closure_Compute()
{
    cout<<"[INFO] Computing and fit closure"<<endl;
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        DoDivision(hQCDClosure_FR[i_fs],hQCDClosure[i_fs][1],hQCDClosure[i_fs][0],false);
    }
}

//===============================================================================
void LL::Step3_QCD_Closure_SaveHistos()
{
    cout<<"[INFO] Save histograms and functions"<<endl;
    TString output_file_name=_savepath+"LL.Step3_QCD_Closure.root";
    output_file = TFile::Open(output_file_name, "recreate");
    output_file->cd();
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i=0;i<2;i++) {
            hQCDClosure[i_fs][i]->Write();
        }
        hQCDClosure_FR[i_fs]->Write();
    }
    output_file->Close();
}

//===============================================================================
void LL::Step3_QCD_Closure_GetHistos()
{
    cout<<"[INFO] Get the closure files"<<endl;
    TString CL_file_name=_savepath+"LL.Step3_QCD_Closure.root";
    TFile *CL_file = TFile::Open(CL_file_name, "read");
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i=0;i<2;i++) {
            _histo_name="hQCDClosure_"+_s_final_state.at(i_fs)+(i==1?"_pass":"_fail");
            hQCDClosure[i_fs][i]=(TH2F*)CL_file->Get(_histo_name);
        }
        _histo_name="hQCDClosure_FR_"+_s_final_state.at(i_fs);
        hQCDClosure_FR[i_fs]=(TH2F*)CL_file->Get(_histo_name);
    }
}

//===============================================================================
void LL::Step3_QCD_vsSR_DeclareHistos()
{
    cout<<"[INFO] Declare histograms and functions"<<endl;
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i=0;i<2;i++) {
            _histo_name="hQCDvsSR_"+_s_final_state.at(i_fs)+(i==1?"_pass":"_fail");
            _histo_labels=";#l_{1} p_{T};l_{2} p_{T}";
            hQCDvsSR[i_fs][i]=new TH2F(_histo_name, _histo_labels, _n_l1pt_bins-1, _l1pt_bins,_n_l2pt_bins-1, _l2pt_bins);
        }
        _histo_name="hQCDvsSR_FR_"+_s_final_state.at(i_fs);
        _histo_labels=";#l_{1} p_{T};l_{2} p_{T}";
        hQCDvsSR_FR[i_fs]=new TH2F(_histo_name, _histo_labels, _n_l1pt_bins-1, _l1pt_bins,_n_l2pt_bins-1, _l2pt_bins);
    }
}

//===============================================================================
void LL::Step3_QCD_vsSR_FillHistos()
{
    cout<<"[INFO] Fill in the histograms and apply fake rates"<<endl;    
    for (size_t i_proc=0;i_proc<_dataFiles.size();i_proc++) {
        TString input_file_name=_path+_dataFiles[i_proc]+_file_name;
        cout<<"[INFO] Processing "<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");

        for (_pass=0;_pass<2;_pass++) {

            TString tree_name=_pass==0?"CRQCDvSRTree":"CRQCDvSROS1Tree";
            
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
                if (LepCombRelIsoPF->at(1)>0.15) continue;

                float DR=std::sqrt((LepEta->at(0)-LepEta->at(1))*(LepEta->at(0)-LepEta->at(1))+(LepPhi->at(0)-LepPhi->at(1))*(LepPhi->at(0)-LepPhi->at(1)));

                _current_final_state = FindFinalState();
                if (_current_final_state<0) continue;

                if (_pass==0) {
                    float FR=fQCD_FR[_current_final_state]->Eval(DR>6?6:DR);
                    float CL=hQCDClosure_FR[_current_final_state]->GetBinContent(hQCDClosure_FR[_current_final_state]->FindBin(LepPt->at(0)>150?150:LepPt->at(0),LepPt->at(1)>150?150:LepPt->at(1)));
                    if (CL<=0) CL=1;
                    _event_weight=FR*CL;
                }
                else {
                    _event_weight=1;
                }

                hQCDvsSR[_current_final_state][_pass]->Fill(LepPt->at(0),LepPt->at(1),_event_weight);
            }
        }
    }

    for (size_t i_proc=0;i_proc<_MCBkgFiles.size();i_proc++) {
        for (_pass=0;_pass<2;_pass++) {
            TString input_file_name=_path+_MCBkgFiles[i_proc]+_file_name;//(_pass==0?_file_name:"HTauTauHMuMu_left.root");
            cout<<"[INFO] Processing "<<input_file_name<<endl;
            input_file = TFile::Open(input_file_name,"read");

            TString tree_name=_pass==0?"CRQCDvSRTree":"CRQCDvSROS1Tree";

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
                if (LepCombRelIsoPF->at(1)>0.15) continue;

                float DR=std::sqrt((LepEta->at(0)-LepEta->at(1))*(LepEta->at(0)-LepEta->at(1))+(LepPhi->at(0)-LepPhi->at(1))*(LepPhi->at(0)-LepPhi->at(1)));

                _current_final_state = FindFinalState();
                if (_current_final_state<0) continue;

                // if (hQCDvsSR[_current_final_state][_pass]->GetBinContent(hQCDvsSR[_current_final_state][_pass]->FindBin(DR))>0) {
                _k_factor = calculate_K_factor(input_file_name);

                _event_weight = (_lumi * 1000 * xsec * _k_factor * PUWeight * genHEPMCweight * dataMCWeight * L1prefiringWeight) / gen_sum_weights;

                if (_pass==0) {
                    float FR=fQCD_FR[_current_final_state]->Eval(DR>6?6:DR);
                    float CL=hQCDClosure_FR[_current_final_state]->GetBinContent(hQCDClosure_FR[_current_final_state]->FindBin(LepPt->at(0)>150?150:LepPt->at(0),LepPt->at(1)>150?150:LepPt->at(1)));
                    if (CL<=0) CL=1;
                    _event_weight*=-1.*FR*CL;
                }
                else {
                    _event_weight*=-1.;
                }

                hQCDvsSR[_current_final_state][_pass]->Fill(DR,_event_weight);
                // }
            }
        }
    }
}

//===============================================================================
void LL::Step3_QCD_vsSR_Compute()
{
    cout<<"[INFO] Computing and fit closure"<<endl;
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        DoDivision(hQCDvsSR_FR[i_fs],hQCDvsSR[i_fs][1],hQCDvsSR[i_fs][0],false);
    }
}

//===============================================================================
void LL::Step3_QCD_vsSR_SaveHistos()
{
    cout<<"[INFO] Save histograms and functions"<<endl;
    TString output_file_name=_savepath+"LL.Step3_QCD_vsSR.root";
    output_file = TFile::Open(output_file_name, "recreate");
    output_file->cd();
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i=0;i<2;i++) {
            hQCDvsSR[i_fs][i]->Write();
        }
        hQCDvsSR_FR[i_fs]->Write();
    }
    output_file->Close();
}

//===============================================================================
void LL::Step3_QCD_vsSR_GetHistos()
{
    cout<<"[INFO] Get the closure files"<<endl;
    TString QCDvsSR_file_name=_savepath+"LL.Step3_QCD_vsSR.root";
    TFile *QCDvsSR_file = TFile::Open(QCDvsSR_file_name, "read");
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i=0;i<2;i++) {
            _histo_name="hQCDvsSR_"+_s_final_state.at(i_fs)+(i==1?"_pass":"_fail");
            hQCDvsSR[i_fs][i]=(TH2F*)QCDvsSR_file->Get(_histo_name);
        }
        _histo_name="hQCDvsSR_FR_"+_s_final_state.at(i_fs);
        hQCDvsSR_FR[i_fs]=(TH2F*)QCDvsSR_file->Get(_histo_name);
    }
}

//===============================================================================
void LL::Step4_Fraction_DeclareHistos()
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
void LL::Step4_Fraction_FillHistos()
{
    cout<<"[INFO] Fill in the histograms"<<endl;
    for (size_t i_proc=0;i_proc<_dataFiles.size();i_proc++) {
        TString input_file_name=_path+_dataFiles[i_proc]+_file_name;
        cout<<"[INFO] Processing "<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");

        TString tree_name="CRQCDTree";

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
            // if (LepCombRelIsoPF->at(1)>0.15) continue;

            _current_final_state = FindFinalState();
            if (_current_final_state<0) continue;

            _current_category = FindCategory();
            if (_current_category<0) continue;

            hFrac[2][_current_final_state][_current_category]->Fill(LLMass);
        }
    }

    for (size_t i_proc=0;i_proc<_MCBkgFiles.size();i_proc++) {

        if (std::find(_TTFiles.begin(),_TTFiles.end(),_MCBkgFiles[i_proc])!=_TTFiles.end()) {
            _current_process=1;
        }
        else {
            _current_process=2;
        }

        TString input_file_name=_path+_MCBkgFiles[i_proc]+_file_name;
        cout<<"[INFO] Processing "<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");

        TString tree_name="CRQCDTree";

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
            // if (LepCombRelIsoPF->at(1)>0.15) continue;

            _current_final_state = FindFinalState();
            if (_current_final_state<0) continue;

            _current_category = FindCategory();
            if (_current_category<0) continue;

            if (_current_process==2) {
                if (hFrac[2][_current_final_state][_current_category]->GetBinContent(hFrac[2][_current_final_state][_current_category]->FindBin(LLMass))>0) {
                    _k_factor = calculate_K_factor(input_file_name);

                    _event_weight = (_lumi * 1000 * xsec * _k_factor * PUWeight * genHEPMCweight * dataMCWeight * L1prefiringWeight) / gen_sum_weights;

                    hFrac[2][_current_final_state][_current_category]->Fill(LLMass,-1.*_event_weight);
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
void LL::Step4_Fraction_Compute()
{
    cout<<"[INFO] Computing fractions"<<endl;
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i_cat=0;i_cat<num_of_categories;i_cat++) {
            hFrac[0][i_fs][i_cat]->Add(hFrac[2][i_fs][i_cat]);
            hFrac[0][i_fs][i_cat]->Add(hFrac[1][i_fs][i_cat],-1.);
            float integral=hFrac[0][i_fs][i_cat]->Integral();
            for (int ibin=1;ibin<=hFrac[0][i_fs][i_cat]->GetNbinsX();ibin++) {
                if (hFrac[0][i_fs][i_cat]->GetBinContent(ibin)<0) {
                    cout<<hFrac[0][i_fs][i_cat]->GetBinContent(ibin)<<" removed"<<endl;
                    hFrac[0][i_fs][i_cat]->SetBinContent(ibin,0);
                }
            }
            if (integral>0)
                hFrac[0][i_fs][i_cat]->Scale(integral/hFrac[0][i_fs][i_cat]->Integral());
            hFrac[2][i_fs][i_cat]->Reset();
            for (int i_bkg=0;i_bkg<2;i_bkg++) hFrac[2][i_fs][i_cat]->Add(hFrac[i_bkg][i_fs][i_cat]);
        }
    }
    for (int i_bkg=0;i_bkg<num_of_fake_bkg+1;i_bkg++) {
        for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
            for (int i_cat=0;i_cat<num_of_categories;i_cat++) {
                DoDivision(hFrac_FR[i_bkg][i_fs][i_cat],hFrac[i_bkg][i_fs][i_cat],hFrac[2][i_fs][i_cat],true);
            }
        }
    }
}

//===============================================================================
void LL::Step4_Fraction_SaveHistos()
{
    cout<<"[INFO] Save histograms and functions"<<endl;
    TString output_file_name=_savepath+"LL.Step4_Fraction.root";
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
void LL::Step4_Fraction_GetHistos()
{
    cout<<"[INFO] Get the fake rate files"<<endl;
    TString Frac_file_name=_savepath+"LL.Step4_Fraction.root";
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
int LL::FindFinalState()
{
    int final_state = -999;
    if (abs(LLFlav)==143) final_state=0;
    // else if (abs(LLFlav)==169) final_state=1;
    else final_state=-1;
    return final_state;
}

//===============================================================================
int LL::FindCategory()
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
TString LL::GetTreeName(int i_bkg)
{
    if (i_bkg==0) return "CRQCDTree";
    else if (i_bkg==1) return "CRWJTree";
    else return "CRTTTree";
}

//===============================================================================
float LL::calculate_K_factor(TString input_file_name)
{
    float k_factor = 1;
    if ( input_file_name.Contains("ZZTo")) k_factor = KFactor_EW_qqZZ * KFactor_QCD_qqZZ_M;
    else if ( input_file_name.Contains("ggTo")) k_factor = KFactor_QCD_ggZZ_Nominal;
    else if ( input_file_name.Contains("ggH")) k_factor = ggH_NNLOPS_weight;
    return k_factor;
}

//===============================================================================
void LL::SetColor(TH1F *h,int color)
{
    h->SetMarkerColor(color);
    h->SetLineColor(color);
    h->SetFillColor(color);
}

//===============================================================================
TString LL::ToFSName(string name)
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
TString LL::ToProcessName(string name)
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
void LL::DoDivision(TH1F* h, TH1F* num, TH1F* den, bool corr)
{
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
void LL::DoDivision(TH2F* h, TH2F* num, TH2F* den, bool corr)
{
    for (int i=0;i<h->GetNbinsX();i++) {
        for (int j=0;j<h->GetNbinsY();j++) {
            double x=num->GetBinContent(i+1,j+1);
            double x0=den->GetBinContent(i+1,j+1);
            double e=num->GetBinError(i+1,j+1);
            double e0=den->GetBinError(i+1,j+1);
            cout<<"x:"<<x<<",e:"<<e<<",x0:"<<x0<<",e0:"<<e0<<endl;
            if (x0<=0) continue;
            h->SetBinContent(i+1,j+1,x/x0);
            if (corr) h->SetBinError(i+1,j+1,std::sqrt(e*e*(x0-x)*(x0-x)+(e0*e0-e*e)*x*x)/x0/x0);
            else h->SetBinError(i+1,j+1,std::sqrt(e*e*x0*x0+e0*e0*x*x)/x0/x0);
        }
    }
}

//===============================================================================
int LL::NumberOfJets()
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