// Include classes
#include <HTauTauHMuMu/AnalysisStep/test/SuperRatio/include/TauTau.h>

// Constructor
//============================================================
TauTau::TauTau():Tree()
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
    // _s_final_state.push_back("etau");
    // _s_final_state.push_back("mutau");
    _s_final_state.push_back("tautau");
    // _s_final_state.push_back("emu");
    // _s_final_state.push_back("2l");
    
    _s_categories.push_back("GGH");
    _s_categories.push_back("VBF");
    _s_categories.push_back("Boost");
    _s_categories.push_back("All");

    _s_fake_bkg.push_back("QCD");
    _s_fake_bkg.push_back("All");

    _s_variables.push_back("Mvis");
    _s_variables.push_back("Mgood");
    
    _lumi=16.8;
    string year="UL2016_postVFP";
}
//--------------------------------------------------------------------------------------

// Destructor
//====================
TauTau::~TauTau()
{
}
//====================

//===============================================================================
void TauTau::SetPaths(string path, string file_name, string savepath)
{
    _path=path;
    _file_name=file_name;
    _savepath=savepath;
}

//===============================================================================
void TauTau::SetFileList(string* data, int ndata, string* bkg, int nbkg)
{
    for (int i=0;i<ndata;i++) _dataFiles.push_back(data[i]);
    for (int i=0;i<nbkg;i++) _MCBkgFiles.push_back(bkg[i]);
}

//===============================================================================
void TauTau::Set_taupt_bin(int n, float* bins)
{
    _n_taupt_bins=n;
    for (int i=0;i<n;i++) _taupt_bins[i]=bins[i];
}

//===============================================================================
void TauTau::Set_lpt_bin(int n, float* bins)
{
    _n_lpt_bins=n;
    for (int i=0;i<n;i++) _lpt_bins[i]=bins[i];
}

//===============================================================================
void TauTau::Set_Mvis_bin(int n, float* bins, bool fine)
{
    if (fine) {
        _n_Mvis_fine_bins=n;
        for (int i=0;i<n;i++) _Mvis_fine_bins[i]=bins[i];
    }
    else {
        _n_Mvis_bins=n;
        for (int i=0;i<n;i++) _Mvis_bins[i]=bins[i];
    }
}

//===============================================================================
void TauTau::Step1_FakeRate_DeclareHistos()
{
    cout<<"[INFO] Declare histograms and functions"<<endl;
    for (int i_bkg=0;i_bkg<num_of_fake_bkg;i_bkg++) {
        for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
            for (int i_njet=0;i_njet<num_of_njet_bins;i_njet++) {
                for (int i=0;i<2;i++) {
                    _histo_name="h_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs)+"_Njet"+std::to_string(i_njet)+(i==1?"_pass":"_fail");
                    _histo_labels=";Leading #tau_{h} p_{T};Weighted events";
                    h[i_bkg][i_fs][i_njet][i]=new TH1F(_histo_name, _histo_labels, _n_taupt_bins-1, _taupt_bins);
                }
                _histo_name="h_FR_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs)+"_Njet"+std::to_string(i_njet);
                _histo_labels=";Leading #tau_{h} p_{T};Weighted events";
                h_FR[i_bkg][i_fs][i_njet]=new TH1F(_histo_name, _histo_labels, _n_taupt_bins-1, _taupt_bins);
                _histo_name="f_FR_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs)+"_Njet"+std::to_string(i_njet);
                f_FR[i_bkg][i_fs][i_njet]=new TF1(_histo_name, "[0]+x*[1]+[2]*TMath::Landau(x,[3],[4])", _taupt_bins[0], _taupt_bins[_n_taupt_bins-1]);
            }
        }
    }
}

//===============================================================================
void TauTau::Step1_FakeRate_FillHistos()
{
    cout<<"[INFO] Fill in the histograms"<<endl;
    for (size_t i_proc=0;i_proc<_dataFiles.size();i_proc++) {
        TString input_file_name=_path+_dataFiles[i_proc]+_file_name;
        cout<<"[INFO] Processing "<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");

        for (int i_bkg=0;i_bkg<num_of_fake_bkg;i_bkg++) {
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

                _current_final_state = FindFinalState();
                if (_current_final_state<0) continue;

                _current_njet_bin=NumberOfJets();

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

        for (int i_bkg=0;i_bkg<num_of_fake_bkg;i_bkg++) {
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

                _current_final_state = FindFinalState();
                if (_current_final_state<0) continue;

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
}

//===============================================================================
void TauTau::Step1_FakeRate_Compute()
{
    cout<<"[INFO] Computing and fit fake rates"<<endl;
    for (int i_bkg=0;i_bkg<num_of_fake_bkg;i_bkg++) {
        for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
            for (int i_njet=0;i_njet<num_of_njet_bins;i_njet++) {
                h_tmp=(TH1F*)h[i_bkg][i_fs][i_njet][1]->Clone();
                h_tmp->Add(h[i_bkg][i_fs][i_njet][0]);
                DoDivision(h_FR[i_bkg][i_fs][i_njet],h[i_bkg][i_fs][i_njet][1],h_tmp,true);
                f_FR[i_bkg][i_fs][i_njet]->SetParameter(2,0.1);
                f_FR[i_bkg][i_fs][i_njet]->SetParameter(3,40);
                f_FR[i_bkg][i_fs][i_njet]->SetParameter(4,5);
                h_FR[i_bkg][i_fs][i_njet]->Fit(f_FR[i_bkg][i_fs][i_njet],"R");
            }
        }
    }
}

//===============================================================================
void TauTau::Step1_FakeRate_SaveHistos()
{
    cout<<"[INFO] Save histograms and functions"<<endl;
    TString output_file_name=_savepath+"TauTau.Step1_FakeRate.root";
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
void TauTau::Step1_FakeRate_GetHistos()
{
    cout<<"[INFO] Get the fake rate files"<<endl;
    TString FR_file_name=_savepath+"TauTau.Step1_FakeRate.root";
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
void TauTau::Step2_Closure_DeclareHistos()
{
    cout<<"[INFO] Declare histograms and functions"<<endl;
    for (int i_bkg=0;i_bkg<num_of_fake_bkg;i_bkg++) {
        for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
            for (int i=0;i<2;i++) {
                _histo_name="hClosure_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs)+(i==1?"_pass":"_fail");
                _histo_labels=";Subleading #tau_{h} p_{T};Weighted events";
                hClosure[i_bkg][i_fs][i]=new TH1F(_histo_name, _histo_labels, _n_lpt_bins-1, _lpt_bins);
            }
            _histo_name="hClosure_FR_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs);
            _histo_labels=";Subleading #tau_{h} p_{T};Weighted events";
            hClosure_FR[i_bkg][i_fs]=new TH1F(_histo_name, _histo_labels, _n_lpt_bins-1, _lpt_bins);
            _histo_name="fClosure_FR_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs);
            fClosure_FR[i_bkg][i_fs]=new TF1(_histo_name, "pol2", _lpt_bins[0], _lpt_bins[_n_lpt_bins-1]);
        }
    }
}

//===============================================================================
void TauTau::Step2_Closure_FillHistos()
{
    cout<<"[INFO] Fill in the histograms and apply fake rates"<<endl;    
    for (size_t i_proc=0;i_proc<_dataFiles.size();i_proc++) {
        TString input_file_name=_path+_dataFiles[i_proc]+_file_name;
        cout<<"[INFO] Processing "<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");

        for (int i_bkg=0;i_bkg<num_of_fake_bkg;i_bkg++) {
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

                _current_final_state = FindFinalState();
                if (_current_final_state<0) continue;

                _current_njet_bin=NumberOfJets();

                if (TauVSjet->at(1)>=5) {
                    _pass=1;
                    _event_weight=1;
                    }
                else {
                    _pass=0;
                    float FR=f_FR[i_bkg][_current_final_state][_current_njet_bin]->Eval(LepPt->at(1)<=150?LepPt->at(1):150);
                    _event_weight=FR/(1.-FR);
                }

                hClosure[i_bkg][_current_final_state][_pass]->Fill(LepPt->at(0),_event_weight);
            }
        }
    }

    for (size_t i_proc=0;i_proc<_MCBkgFiles.size();i_proc++) {
        TString input_file_name=_path+_MCBkgFiles[i_proc]+_file_name;
        cout<<"[INFO] Processing "<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");

        for (int i_bkg=0;i_bkg<num_of_fake_bkg;i_bkg++) {
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

                _current_final_state = FindFinalState();
                if (_current_final_state<0) continue;

                _current_njet_bin=NumberOfJets();

                if (hClosure[i_bkg][_current_final_state][_pass]->GetBinContent(hClosure[i_bkg][_current_final_state][_pass]->FindBin(LepPt->at(0)))>0) {
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
                    hClosure[i_bkg][_current_final_state][_pass]->Fill(LepPt->at(0),_event_weight);
                }
            }
        }
    }
}

//===============================================================================
void TauTau::Step2_Closure_Compute()
{
    cout<<"[INFO] Computing and fit closure"<<endl;
    for (int i_bkg=0;i_bkg<num_of_fake_bkg;i_bkg++) {
        for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
            DoDivision(hClosure_FR[i_bkg][i_fs],hClosure[i_bkg][i_fs][1],hClosure[i_bkg][i_fs][0],false);
            hClosure_FR[i_bkg][i_fs]->Fit(fClosure_FR[i_bkg][i_fs],"R");
        }
    }
}

//===============================================================================
void TauTau::Step2_Closure_SaveHistos()
{
    cout<<"[INFO] Save histograms and functions"<<endl;
    TString output_file_name=_savepath+"TauTau.Step2_Closure.root";
    output_file = TFile::Open(output_file_name, "recreate");
    output_file->cd();
    for (int i_bkg=0;i_bkg<num_of_fake_bkg;i_bkg++) {
        for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
            for (int i=0;i<2;i++) {
                hClosure[i_bkg][i_fs][i]->Write();
            }
            hClosure_FR[i_bkg][i_fs]->Write();
            fClosure_FR[i_bkg][i_fs]->Write();
        }
    }
    output_file->Close();
}

//===============================================================================
void TauTau::Step2_Closure_GetHistos()
{
    cout<<"[INFO] Get the closure files"<<endl;
    TString CL_file_name=_savepath+"TauTau.Step2_Closure.root";
    TFile *CL_file = TFile::Open(CL_file_name, "read");
    for (int i_bkg=0;i_bkg<num_of_fake_bkg;i_bkg++) {
        for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
            for (int i=0;i<2;i++) {
                _histo_name="hClosure_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs)+(i==1?"_pass":"_fail");
                hClosure[i_bkg][i_fs][i]=(TH1F*)CL_file->Get(_histo_name);
            }
            _histo_name="hClosure_FR_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs);
            hClosure_FR[i_bkg][i_fs]=(TH1F*)CL_file->Get(_histo_name);
            _histo_name="fClosure_FR_"+_s_fake_bkg.at(i_bkg)+"_"+_s_final_state.at(i_fs);
            fClosure_FR[i_bkg][i_fs]=(TF1*)CL_file->Get(_histo_name);
        }
    }
}

//===============================================================================
void TauTau::Step3_QCD_FakeRate_DeclareHistos()
{
    cout<<"[INFO] Declare histograms and functions"<<endl;
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i_njet=0;i_njet<num_of_njet_bins;i_njet++) {
            for (int i=0;i<2;i++) {
                _histo_name="hQCD_"+_s_final_state.at(i_fs)+"_Njet"+std::to_string(i_njet)+(i==1?"_pass":"_fail");
                _histo_labels=";Leading #tau_{h} p_{T};Weighted events";
                hQCD[i_fs][i_njet][i]=new TH1F(_histo_name, _histo_labels, _n_taupt_bins-1, _taupt_bins);
            }
            _histo_name="hQCD_FR_"+_s_final_state.at(i_fs)+"_Njet"+std::to_string(i_njet);
            _histo_labels=";Leading #tau_{h} p_{T};Weighted events";
            hQCD_FR[i_fs][i_njet]=new TH1F(_histo_name, _histo_labels, _n_taupt_bins-1, _taupt_bins);
            _histo_name="fQCD_FR_"+_s_final_state.at(i_fs)+"_Njet"+std::to_string(i_njet);
            fQCD_FR[i_fs][i_njet]=new TF1(_histo_name, "[0]+x*[1]+[2]*TMath::Landau(x,[3],[4])", _taupt_bins[0], _taupt_bins[_n_taupt_bins-1]);
        }
    }
}

//===============================================================================
void TauTau::Step3_QCD_FakeRate_FillHistos()
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
void TauTau::Step3_QCD_FakeRate_Compute()
{
    cout<<"[INFO] Computing and fit fake rates"<<endl;
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i_njet=0;i_njet<num_of_njet_bins;i_njet++) {
            h_tmp=(TH1F*)hQCD[i_fs][i_njet][1]->Clone();
            h_tmp->Add(hQCD[i_fs][i_njet][0]);
            DoDivision(hQCD_FR[i_fs][i_njet],hQCD[i_fs][i_njet][1],h_tmp,true);
            fQCD_FR[i_fs][i_njet]->SetParameter(2,0.1);
            fQCD_FR[i_fs][i_njet]->SetParameter(3,40);
            fQCD_FR[i_fs][i_njet]->SetParameter(4,5);
            hQCD_FR[i_fs][i_njet]->Fit(fQCD_FR[i_fs][i_njet],"R");
        }
    }
}

//===============================================================================
void TauTau::Step3_QCD_FakeRate_SaveHistos()
{
    cout<<"[INFO] Save histograms and functions"<<endl;
    TString output_file_name=_savepath+"TauTau.Step3_QCD_FakeRate.root";
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
void TauTau::Step3_QCD_FakeRate_GetHistos()
{
    cout<<"[INFO] Get the fake rate files"<<endl;
    TString QCDFR_file_name=_savepath+"TauTau.Step3_QCD_FakeRate.root";
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
void TauTau::Step3_QCD_Closure_DeclareHistos()
{
    cout<<"[INFO] Declare histograms and functions"<<endl;
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i=0;i<2;i++) {
            _histo_name="hQCDClosure_"+_s_final_state.at(i_fs)+(i==1?"_pass":"_fail");
            _histo_labels=";Subleading #tau p_{T};Weighted events";
            hQCDClosure[i_fs][i]=new TH1F(_histo_name, _histo_labels, _n_lpt_bins-1, _lpt_bins);
        }
        _histo_name="hQCDClosure_FR_"+_s_final_state.at(i_fs);
        _histo_labels=";Subleading #tau p_{T};Weighted events";
        hQCDClosure_FR[i_fs]=new TH1F(_histo_name, _histo_labels, _n_lpt_bins-1, _lpt_bins);
        _histo_name="fQCDClosure_FR_"+_s_final_state.at(i_fs);
        fQCDClosure_FR[i_fs]=new TF1(_histo_name, "pol2", _lpt_bins[0], _lpt_bins[_n_lpt_bins-1]);
    }
}

//===============================================================================
void TauTau::Step3_QCD_Closure_FillHistos()
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

            if (TauVSjet->at(1)>=5) {
                _pass=1;
                _event_weight=1;
                }
            else {
                _pass=0;
                float FR=fQCD_FR[_current_final_state][_current_njet_bin]->Eval(LepPt->at(1)<=150?LepPt->at(1):150);
                _event_weight=FR/(1.-FR);
            }

            hQCDClosure[_current_final_state][_pass]->Fill(LepPt->at(0),_event_weight);
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

            if (hQCDClosure[_current_final_state][_pass]->GetBinContent(hQCDClosure[_current_final_state][_pass]->FindBin(LepPt->at(0)))>0) {
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

                hQCDClosure[_current_final_state][_pass]->Fill(LepPt->at(0),_event_weight);
            }
        }
    }
}

//===============================================================================
void TauTau::Step3_QCD_Closure_Compute()
{
    cout<<"[INFO] Computing and fit closure"<<endl;
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        DoDivision(hQCDClosure_FR[i_fs],hQCDClosure[i_fs][1],hQCDClosure[i_fs][0],false);
        hQCDClosure_FR[i_fs]->Fit(fQCDClosure_FR[i_fs],"R");
    }
}

//===============================================================================
void TauTau::Step3_QCD_Closure_SaveHistos()
{
    cout<<"[INFO] Save histograms and functions"<<endl;
    TString output_file_name=_savepath+"TauTau.Step3_QCD_Closure.root";
    output_file = TFile::Open(output_file_name, "recreate");
    output_file->cd();
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i=0;i<2;i++) {
            hQCDClosure[i_fs][i]->Write();
        }
        hQCDClosure_FR[i_fs]->Write();
        fQCDClosure_FR[i_fs]->Write();
    }
    output_file->Close();
}

//===============================================================================
void TauTau::Step3_QCD_Closure_GetHistos()
{
    cout<<"[INFO] Get the closure files"<<endl;
    TString QCDCL_file_name=_savepath+"TauTau.Step3_QCD_Closure.root";
    TFile *QCDCL_file = TFile::Open(QCDCL_file_name, "read");
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i=0;i<2;i++) {
            _histo_name="hQCDClosure_"+_s_final_state.at(i_fs)+(i==1?"_pass":"_fail");
            hQCDClosure[i_fs][i]=(TH1F*)QCDCL_file->Get(_histo_name);
        }
        _histo_name="hQCDClosure_FR_"+_s_final_state.at(i_fs);
        hQCDClosure_FR[i_fs]=(TH1F*)QCDCL_file->Get(_histo_name);
        _histo_name="fQCDClosure_FR_"+_s_final_state.at(i_fs);
        fQCDClosure_FR[i_fs]=(TF1*)QCDCL_file->Get(_histo_name);
    }
}

//===============================================================================
void TauTau::Step3_QCD_vsSR_DeclareHistos()
{
    cout<<"[INFO] Declare histograms and functions"<<endl;
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        for (int i=0;i<2;i++) {
            _histo_name="hQCDvsSR_"+_s_final_state.at(i_fs)+(i==1?"_pass":"_fail");
            _histo_labels=";M_{vis};Weighted events";
            hQCDvsSR[i_fs][i]=new TH1F(_histo_name, _histo_labels, _n_Mvis_fine_bins-1, _Mvis_fine_bins);
        }
        _histo_name="hQCDvsSR_FR_"+_s_final_state.at(i_fs);
        _histo_labels=";M_{vis};Weighted events";
        hQCDvsSR_FR[i_fs]=new TH1F(_histo_name, _histo_labels, _n_Mvis_fine_bins-1, _Mvis_fine_bins);
        _histo_name="fQCDvsSR_FR_"+_s_final_state.at(i_fs);
        fQCDvsSR_FR[i_fs]=new TF1(_histo_name, "([0]+[1]*x+[2]*x*x)*(x<90)+([3]*x+[0]+8100*[2]+90*[1]-90*[3])*(x>=90)", _Mvis_fine_bins[0], _Mvis_fine_bins[_n_Mvis_fine_bins-1]);
    }
}

//===============================================================================
void TauTau::Step3_QCD_vsSR_FillHistos()
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

                _current_final_state = FindFinalState();
                if (_current_final_state<0) continue;

                _current_njet_bin=NumberOfJets();

                if (_pass==0) {
                    float FR=fQCD_FR[_current_final_state][_current_njet_bin]->Eval(LepPt->at(1)<=150?LepPt->at(1):150);
                    float CL=fQCDClosure_FR[_current_final_state]->Eval(LepPt->at(0)<=150?LepPt->at(0):150);
                    // if (LepPt->at(0)<=23) 
                    //     CL=hQCDClosure_FR[_current_final_state][_current_taupt_bin]->GetBinContent(1);
                    // else
                    //     CL=fQCDClosure_FR[_current_final_state]->Eval(LepPt->at(0)<=150?LepPt->at(0):150);
                    _event_weight=FR/(1.-FR)*CL;
                }
                else {
                    _event_weight=1;
                }

                hQCDvsSR[_current_final_state][_pass]->Fill(LLMass,_event_weight);
            }
        }
    }

    for (size_t i_proc=0;i_proc<_MCBkgFiles.size();i_proc++) {
        for (_pass=0;_pass<2;_pass++) {
            TString input_file_name=_path+_MCBkgFiles[i_proc]+_file_name;
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

                _current_final_state = FindFinalState();
                if (_current_final_state<0) continue;

                _current_njet_bin=NumberOfJets();

                if (hQCDvsSR[_current_final_state][_pass]->GetBinContent(hQCDvsSR[_current_final_state][_pass]->FindBin(LLMass))>0) {
                    _k_factor = calculate_K_factor(input_file_name);

                    _event_weight = (_lumi * 1000 * xsec * _k_factor * PUWeight * genHEPMCweight * dataMCWeight * L1prefiringWeight) / gen_sum_weights;

                    if (_pass==0) {
                        float FR=fQCD_FR[_current_final_state][_current_njet_bin]->Eval(LepPt->at(1)<=150?LepPt->at(1):150);
                        float CL=fQCDClosure_FR[_current_final_state]->Eval(LepPt->at(0)<=150?LepPt->at(0):150);
                        // if (LepPt->at(0)<=23) 
                        //     CL=hQCDClosure_FR[_current_final_state][_current_taupt_bin]->GetBinContent(1);
                        // else
                        //     CL=fQCDClosure_FR[_current_final_state][_current_taupt_bin]->Eval(LepPt->at(0)<=150?LepPt->at(0):150);
                        _event_weight*=-1.*FR/(1.-FR)*CL;
                    }
                    else {
                        _event_weight*=-1.;
                    }

                    hQCDvsSR[_current_final_state][_pass]->Fill(LLMass,_event_weight);
                }
            }
        }
    }
}

//===============================================================================
void TauTau::Step3_QCD_vsSR_Compute()
{
    cout<<"[INFO] Computing and fit closure"<<endl;
    for (int i_fs=0;i_fs<num_of_final_states;i_fs++) {
        DoDivision(hQCDvsSR_FR[i_fs],hQCDvsSR[i_fs][1],hQCDvsSR[i_fs][0],false);
        hQCDvsSR_FR[i_fs]->Fit(fQCDvsSR_FR[i_fs],"R");
    }
}

void TauTau::Step3_QCD_vsSR_SaveHistos()
{
    cout<<"[INFO] Save histograms and functions"<<endl;
    TString output_file_name=_savepath+"TauTau.Step3_QCD_vsSR.root";
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
void TauTau::Step3_QCD_vsSR_GetHistos()
{
    cout<<"[INFO] Get the closure files"<<endl;
    TString QCDvsSR_file_name=_savepath+"TauTau.Step3_QCD_vsSR.root";
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
int TauTau::FindFinalState()
{
    int final_state = -999;
    if (abs(LLFlav)==225) final_state=0;
    else final_state=-1;
    return final_state;
}

//===============================================================================
int TauTau::FindCategory()
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
float TauTau::calculate_K_factor(TString input_file_name)
{
    float k_factor = 1;
    if ( input_file_name.Contains("ZZTo")) k_factor = KFactor_EW_qqZZ * KFactor_QCD_qqZZ_M;
    else if ( input_file_name.Contains("ggTo")) k_factor = KFactor_QCD_ggZZ_Nominal;
    else if ( input_file_name.Contains("ggH")) k_factor = ggH_NNLOPS_weight;
    return k_factor;
}

//===============================================================================
void TauTau::SetColor(TH1F *h,int color)
{
    h->SetMarkerColor(color);
    h->SetLineColor(color);
    h->SetFillColor(color);
}

//===============================================================================
TString TauTau::ToFSName(string name)
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
TString TauTau::ToProcessName(string name)
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
void TauTau::DoDivision(TH1F* h, TH1F* num, TH1F* den, bool corr)
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
int TauTau::NumberOfJets()
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