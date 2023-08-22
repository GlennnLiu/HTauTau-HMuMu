// Include classes
#include <HTauTauHMuMu/AnalysisStep/test/Plotter/include/TreeToHist.h>

// Constructor
//============================================================
TreeToHist::TreeToHist():Tree()
{    
    _current_process = -999;
    _current_final_state = -999;
    _current_category = -999;

    _s_process.clear();
    _s_process.push_back("Htt");
    _s_process.push_back("Hmm");
    _s_process.push_back("Hww");
    _s_process.push_back("diBoson");
    _s_process.push_back("triBoson");
    _s_process.push_back("top");
    _s_process.push_back("EWKZee");
    _s_process.push_back("EWKZtt");
    _s_process.push_back("EWKZmm");
    _s_process.push_back("fake");
    _s_process.push_back("DYee");
    _s_process.push_back("DYtt");
    _s_process.push_back("DYmm");
    _s_process.push_back("TotalMC");
    _s_process.push_back("Data");
    
    _s_final_state.push_back("mumu");
    _s_final_state.push_back("etau");
    _s_final_state.push_back("mutau");
    _s_final_state.push_back("tautau");
    _s_final_state.push_back("emu");
    _s_final_state.push_back("tautaucomb");
    _s_final_state.push_back("2l");
    
    _s_category.push_back("GGH");
    _s_category.push_back("VBF_ptHl200");
    _s_category.push_back("VBF_ptHg200");
    _s_category.push_back("Boost_1j");
    _s_category.push_back("Boost_2j");
    _s_category.push_back("All");
    
    _lumi=16.8;
    string year="UL2016_postVFP";

    fakeBkg=new FakeBkg(year,"/eos/home-g/geliu/LepUni/SuperRatio/");
}
//--------------------------------------------------------------------------------------

// Destructor
//====================
TreeToHist::~TreeToHist()
{
}
//====================

//===============================================================================
void TreeToHist::SetPaths(string path, string file_name, string savepath, string directory)
{
    _path=path;
    _file_name=file_name;
    _savepath=savepath;
    _directory=directory;
}
//===============================================================================
void TreeToHist::ToHistos(string* process, int n_proc, int i_proc)
{
    DeclareHistos(i_proc);
    
    for (int i=0;i<n_proc;i++) {
        TString input_file_name=_path+process[i]+_file_name;
        cout<<i<<","<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");
        
        hCounters = (TH1F*)input_file->Get(_directory+"/Counters");
        gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);
        
        input_tree = (TTree*)input_file->Get(_directory+"/candTree");
        Init( input_tree, input_file_name, true);
        
        if (fChain == 0) {return;}
        
        Long64_t nentries = fChain->GetEntriesFast();
        cout<<"[INFO] Processing "<<input_file_name<<", "<<_directory<<", total event: "<<nentries<<endl;
        Long64_t nbytes = 0, nb = 0;
        
        for (Long64_t jentry=0; jentry<nentries;jentry++) {
            Long64_t ientry = LoadTree(jentry);
            if (ientry%50000==0) cout<<ientry<<endl;
            if (ientry < 0) break;
            nb = fChain->GetEntry(jentry);
            nbytes += nb;

            if (!pass_Trigger) continue;
            if (nCleanedJetsPt25BTagged_bTagSF>0) continue;

            _current_final_state = FindFinalState();
            if (_current_final_state<0) {continue;}
            if (i_proc!=Settings::DYee && i_proc!=Settings::EWKZee) {_current_process=i_proc;}
            else {
                _current_process = FindProcess(i_proc);
            } 

            _current_category=FindCategory(NumberOfJets(),true);

            TLorentzVector l1,l2,MET;
            l1.SetPtEtaPhiM(LepPt->at(0),LepEta->at(0),LepPhi->at(0),LepM->at(0));
            l2.SetPtEtaPhiM(LepPt->at(1),LepEta->at(1),LepPhi->at(1),LepM->at(1));
            MET.SetPtEtaPhiM(PFMET,0,PFMETPhi,0);
            float MtlMET=(l1+MET).Mt();
            float MtllMET=(l1+l2+MET).Mt();

            if (!pass_Extra_Cuts(_current_final_state,Pzeta1,MtllMET)) continue;

            if (i_proc != Settings::Data) {
                _k_factor = calculate_K_factor(input_file_name);
                _event_weight = (_lumi * 1000 * xsec * _k_factor * PUWeight * genHEPMCweight * dataMCWeight * L1prefiringWeight) / gen_sum_weights;
            }
            else {
                _event_weight = 1;
            }

            histos_1D[_current_process][_current_final_state][_current_category][0]->Fill(LLMass,_event_weight);
            histos_1D[_current_process][_current_final_state][_current_category][1]->Fill(LLGoodMass,_event_weight);
            histos_1D[_current_process][_current_final_state][_current_category][2]->Fill(Pzeta1,_event_weight);
            histos_1D[_current_process][_current_final_state][_current_category][3]->Fill(MtlMET,_event_weight);
            histos_1D[_current_process][_current_final_state][_current_category][4]->Fill(MtllMET,_event_weight);
            histos_1D[_current_process][_current_final_state][_current_category][5]->Fill(LLPt,_event_weight);
        }
        input_file->Close();
        cout<<"[INFO] Processing "<<input_file_name<<" done"<<endl;
    }

    SumGroups(i_proc);
    saveHistos(i_proc);
}

//===============================================================================
void TreeToHist::ToHistosFake(string* data, int n_data, string* bkg, int n_bkg, int i_proc)
{
    DeclareHistos(i_proc);
    
    for (int i=0;i<n_data;i++) {
        TString input_file_name=_path+data[i]+_file_name;
        cout<<i<<","<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");
        
        input_tree = (TTree*)input_file->Get("CRAPPOSTree/candTree");
        Init( input_tree, input_file_name, true);
        
        if (fChain == 0) {return;}
        
        Long64_t nentries = fChain->GetEntriesFast();
        cout<<"[INFO] Processing "<<input_file_name<<", "<<_directory<<", total event: "<<nentries<<endl;
        Long64_t nbytes = 0, nb = 0;
        
        for (Long64_t jentry=0; jentry<nentries;jentry++) {
            Long64_t ientry = LoadTree(jentry);
            if (ientry%50000==0) cout<<ientry<<endl;
            if (ientry < 0) break;
            nb = fChain->GetEntry(jentry);
            nbytes += nb;
            
            if (!pass_Trigger) continue;
            if (nCleanedJetsPt25BTagged_bTagSF>0) continue;

            _current_final_state = FindFinalState();
            int njet = NumberOfJets();
            _current_category=FindCategory(njet,true);

            TLorentzVector l1,l2,MET;
            l1.SetPtEtaPhiM(LepPt->at(0),LepEta->at(0),LepPhi->at(0),LepM->at(0));
            l2.SetPtEtaPhiM(LepPt->at(1),LepEta->at(1),LepPhi->at(1),LepM->at(1));
            MET.SetPtEtaPhiM(PFMET,0,PFMETPhi,0);
            float MtlMET=(l1+MET).Mt();
            float MtllMET=(l1+l2+MET).Mt();

            if (_current_final_state==Settings::fsetau || _current_final_state==Settings::fsmutau) {
                int fs=_current_final_state==Settings::fsetau?0:1;
                int cat=FindCategory(njet,true);
                float DR=std::sqrt((LepEta->at(0)-LepEta->at(1))*(LepEta->at(0)-LepEta->at(1))+(LepPhi->at(0)-LepPhi->at(1))*(LepPhi->at(0)-LepPhi->at(1)));

                _event_weight=fakeBkg->get_LTau_FR(fs,cat,njet,LepPt->at(0),LepPt->at(1),DR,MtLMET,LLMass);
            }
            else if (_current_final_state==Settings::fstautau) {
                _event_weight=fakeBkg->get_TauTau_FR(njet,LepPt->at(0),LepPt->at(1),LLMass);
            }
            else {
                continue;
            }

            if (!pass_Extra_Cuts(_current_final_state,Pzeta1,MtllMET)) continue;

            histos_1D[Settings::fake][_current_final_state][_current_category][0]->Fill(LLMass,_event_weight);
            histos_1D[Settings::fake][_current_final_state][_current_category][1]->Fill(LLGoodMass,_event_weight);
            histos_1D[Settings::fake][_current_final_state][_current_category][2]->Fill(Pzeta1,_event_weight);
            histos_1D[Settings::fake][_current_final_state][_current_category][3]->Fill(MtlMET,_event_weight);
            histos_1D[Settings::fake][_current_final_state][_current_category][4]->Fill(MtllMET,_event_weight);
            histos_1D[Settings::fake][_current_final_state][_current_category][5]->Fill(LLPt,_event_weight);
        }
        input_file->Close();
        cout<<"[INFO] Processing "<<input_file_name<<" done"<<endl;
    }

    for (int i=0;i<n_bkg;i++) {
        TString input_file_name=_path+bkg[i]+_file_name;
        cout<<i<<","<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");
        
        hCounters = (TH1F*)input_file->Get("CRAPPOSTree/Counters");
        gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);
        
        input_tree = (TTree*)input_file->Get("CRAPPOSTree/candTree");
        Init( input_tree, input_file_name, true);
        
        if (fChain == 0) {return;}
        
        Long64_t nentries = fChain->GetEntriesFast();
        cout<<"[INFO] Processing "<<input_file_name<<", "<<_directory<<", total event: "<<nentries<<endl;
        Long64_t nbytes = 0, nb = 0;
        
        for (Long64_t jentry=0; jentry<nentries;jentry++) {
            Long64_t ientry = LoadTree(jentry);
            if (ientry%50000==0) cout<<ientry<<endl;
            if (ientry < 0) break;
            nb = fChain->GetEntry(jentry);
            nbytes += nb;
            
            if (!pass_Trigger) continue;
            if (nCleanedJetsPt25BTagged_bTagSF>0) continue;

            _current_final_state = FindFinalState();
            int njet = NumberOfJets();
            _current_category=FindCategory(njet,true);

            TLorentzVector l1,l2,MET;
            l1.SetPtEtaPhiM(LepPt->at(0),LepEta->at(0),LepPhi->at(0),LepM->at(0));
            l2.SetPtEtaPhiM(LepPt->at(1),LepEta->at(1),LepPhi->at(1),LepM->at(1));
            MET.SetPtEtaPhiM(PFMET,0,PFMETPhi,0);
            float MtlMET=(l1+MET).Mt();
            float MtllMET=(l1+l2+MET).Mt();

            if (_current_final_state==Settings::fsetau || _current_final_state==Settings::fsmutau) {
                int fs=_current_final_state==Settings::fsetau?0:1;
                int cat=FindCategory(njet,true);
                float DR=std::sqrt((LepEta->at(0)-LepEta->at(1))*(LepEta->at(0)-LepEta->at(1))+(LepPhi->at(0)-LepPhi->at(1))*(LepPhi->at(0)-LepPhi->at(1)));

                _event_weight=fakeBkg->get_LTau_FR(fs,cat,njet,LepPt->at(0),LepPt->at(1),DR,MtLMET,LLMass);
            }
            else if (_current_final_state==Settings::fstautau) {
                _event_weight=fakeBkg->get_TauTau_FR(njet,LepPt->at(0),LepPt->at(1),LLMass);
            }
            else {
                continue;
            }

            if (!pass_Extra_Cuts(_current_final_state,Pzeta1,MtllMET)) continue;

            _k_factor = calculate_K_factor(input_file_name);
            _event_weight *= (_lumi * 1000 * xsec * _k_factor * PUWeight * genHEPMCweight * dataMCWeight * L1prefiringWeight) / gen_sum_weights;

            histos_1D[Settings::fake][_current_final_state][_current_category][0]->Fill(LLMass,-1.*_event_weight);
            histos_1D[Settings::fake][_current_final_state][_current_category][1]->Fill(LLGoodMass,-1.*_event_weight);
            histos_1D[Settings::fake][_current_final_state][_current_category][2]->Fill(Pzeta1,-1.*_event_weight);
            histos_1D[Settings::fake][_current_final_state][_current_category][3]->Fill(MtlMET,-1.*_event_weight);
            histos_1D[Settings::fake][_current_final_state][_current_category][4]->Fill(MtllMET,-1.*_event_weight);
            histos_1D[Settings::fake][_current_final_state][_current_category][5]->Fill(LLPt,-1.*_event_weight);
        }
        input_file->Close();
        cout<<"[INFO] Processing "<<input_file_name<<" done"<<endl;
    }

    float tot_integral = 0.;
    float tot_integral_iso = 0.;

    for (int i=0;i<n_data;i++) {
        TString input_file_name=_path+data[i]+_file_name;
        cout<<i<<","<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");
        
        input_tree = (TTree*)input_file->Get("CRQCDTree/candTree");
        Init( input_tree, input_file_name, true);
        
        if (fChain == 0) {return;}
        
        Long64_t nentries = fChain->GetEntriesFast();
        cout<<"[INFO] Processing "<<input_file_name<<", "<<_directory<<", total event: "<<nentries<<endl;
        Long64_t nbytes = 0, nb = 0;
        
        for (Long64_t jentry=0; jentry<nentries;jentry++) {
            Long64_t ientry = LoadTree(jentry);
            if (ientry%50000==0) cout<<ientry<<endl;
            if (ientry < 0) break;
            nb = fChain->GetEntry(jentry);
            nbytes += nb;
            
            if (!pass_Trigger) continue;
            if (nCleanedJetsPt25BTagged_bTagSF>0) continue;

            _current_final_state = FindFinalState();
            int njet = NumberOfJets();
            _current_category=FindCategory(njet,true);

            if (_current_final_state==Settings::fsemu) {
                int cat=FindCategory(njet,true);
                float DR=std::sqrt((LepEta->at(0)-LepEta->at(1))*(LepEta->at(0)-LepEta->at(1))+(LepPhi->at(0)-LepPhi->at(1))*(LepPhi->at(0)-LepPhi->at(1)));

                _event_weight=fakeBkg->get_LL_FR(cat,DR,LepPt->at(0),LepPt->at(1),LLMass);
            }
            else {
                continue;
            }

            tot_integral+=1;
            if (LepCombRelIsoPF->at(1)>0.15) tot_integral_iso+=1;

            TLorentzVector l1,l2,MET;
            l1.SetPtEtaPhiM(LepPt->at(0),LepEta->at(0),LepPhi->at(0),LepM->at(0));
            l2.SetPtEtaPhiM(LepPt->at(1),LepEta->at(1),LepPhi->at(1),LepM->at(1));
            MET.SetPtEtaPhiM(PFMET,0,PFMETPhi,0);
            float MtlMET=(l1+MET).Mt();
            float MtllMET=(l1+l2+MET).Mt();

            if (!pass_Extra_Cuts(_current_final_state,Pzeta1,MtllMET)) continue;

            histos_1D[Settings::fake][_current_final_state][_current_category][0]->Fill(LLMass,_event_weight);
            histos_1D[Settings::fake][_current_final_state][_current_category][1]->Fill(LLGoodMass,_event_weight);
            histos_1D[Settings::fake][_current_final_state][_current_category][2]->Fill(Pzeta1,_event_weight);
            histos_1D[Settings::fake][_current_final_state][_current_category][3]->Fill(MtlMET,_event_weight);
            histos_1D[Settings::fake][_current_final_state][_current_category][4]->Fill(MtllMET, _event_weight);
            histos_1D[Settings::fake][_current_final_state][_current_category][5]->Fill(LLPt, _event_weight);
        }
        input_file->Close();
        cout<<"[INFO] Processing "<<input_file_name<<" done"<<endl;
    }

    for (int i=0;i<n_bkg;i++) {
        TString input_file_name=_path+bkg[i]+_file_name;
        cout<<i<<","<<input_file_name<<endl;
        input_file = TFile::Open(input_file_name,"read");
        
        hCounters = (TH1F*)input_file->Get("CRQCDTree/Counters");
        gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);
        
        input_tree = (TTree*)input_file->Get("CRQCDTree/candTree");
        Init( input_tree, input_file_name, true);
        
        if (fChain == 0) {return;}
        
        Long64_t nentries = fChain->GetEntriesFast();
        cout<<"[INFO] Processing "<<input_file_name<<", "<<_directory<<", total event: "<<nentries<<endl;
        Long64_t nbytes = 0, nb = 0;
        
        for (Long64_t jentry=0; jentry<nentries;jentry++) {
            Long64_t ientry = LoadTree(jentry);
            if (ientry%50000==0) cout<<ientry<<endl;
            if (ientry < 0) break;
            nb = fChain->GetEntry(jentry);
            nbytes += nb;
            
            if (!pass_Trigger) continue;
            if (nCleanedJetsPt25BTagged_bTagSF>0) continue;

            _current_final_state = FindFinalState();
            int njet = NumberOfJets();
            _current_category=FindCategory(njet,true);

            if (_current_final_state==Settings::fsemu) {
                int cat=FindCategory(njet,true);
                float DR=std::sqrt((LepEta->at(0)-LepEta->at(1))*(LepEta->at(0)-LepEta->at(1))+(LepPhi->at(0)-LepPhi->at(1))*(LepPhi->at(0)-LepPhi->at(1)));

                _event_weight=fakeBkg->get_LL_FR(cat,DR,LepPt->at(0),LepPt->at(1),LLMass);
            }
            else {
                continue;
            }

            TLorentzVector l1,l2,MET;
            l1.SetPtEtaPhiM(LepPt->at(0),LepEta->at(0),LepPhi->at(0),LepM->at(0));
            l2.SetPtEtaPhiM(LepPt->at(1),LepEta->at(1),LepPhi->at(1),LepM->at(1));
            MET.SetPtEtaPhiM(PFMET,0,PFMETPhi,0);
            float MtlMET=(l1+MET).Mt();
            float MtllMET=(l1+l2+MET).Mt();

            if (!pass_Extra_Cuts(_current_final_state,Pzeta1,MtllMET)) continue;

            _k_factor = calculate_K_factor(input_file_name);
            float _weight = (_lumi * 1000 * xsec * _k_factor * PUWeight * genHEPMCweight * dataMCWeight * L1prefiringWeight) / gen_sum_weights;
            tot_integral-=_weight;
            if (LepCombRelIsoPF->at(1)>0.15) tot_integral_iso-=_weight;
            _event_weight *= _weight;

            histos_1D[Settings::fake][_current_final_state][_current_category][0]->Fill(LLMass,-1.*_event_weight);
            histos_1D[Settings::fake][_current_final_state][_current_category][1]->Fill(LLGoodMass,-1.*_event_weight);
            histos_1D[Settings::fake][_current_final_state][_current_category][2]->Fill(Pzeta1,-1.*_event_weight);
            histos_1D[Settings::fake][_current_final_state][_current_category][3]->Fill(MtlMET,-1.*_event_weight);
            histos_1D[Settings::fake][_current_final_state][_current_category][4]->Fill(MtllMET, -1.*_event_weight);
            histos_1D[Settings::fake][_current_final_state][_current_category][5]->Fill(LLPt, -1.*_event_weight);
        }
        input_file->Close();
        cout<<"[INFO] Processing "<<input_file_name<<" done"<<endl;
    }

    for (int i_cat=0;i_cat<Settings::catAll;i_cat++) {
        for (int i_var=0;i_var<num_of_variables;i_var++) {
            histos_1D[Settings::fake][Settings::fsemu][i_cat][i_var]->Scale(tot_integral_iso/tot_integral);
        }
    }

    SumGroups(i_proc);
    saveHistos(i_proc);
}

//===============================================================================
void TreeToHist::DeclareHistos(int i_proc)
{
    if (i_proc!=Settings::DYee && i_proc!=Settings::EWKZee) {
        for (int i_fs=0;i_fs<=Settings::fs2l;i_fs++) {
            for (int i_cat=0;i_cat<=Settings::catAll;i_cat++) {
                _histo_name="LLMass_"+_s_process.at(i_proc)+"_"+_s_final_state.at(i_fs)+"_"+_s_category.at(i_cat);
                _histo_labels=";"+Plots::LLMass().var_X_label + ";" + Plots::LLMass().var_Y_label;
                histos_1D[i_proc][i_fs][i_cat][0]=new TH1F(_histo_name, _histo_labels, Plots::LLMass().var_N_bin, Plots::LLMass().var_min, Plots::LLMass().var_max);

                _histo_name="LLGoodMass_"+_s_process.at(i_proc)+"_"+_s_final_state.at(i_fs)+"_"+_s_category.at(i_cat);
                _histo_labels=";"+Plots::LLGoodMass().var_X_label + ";" + Plots::LLGoodMass().var_Y_label;
                histos_1D[i_proc][i_fs][i_cat][1]=new TH1F(_histo_name, _histo_labels, Plots::LLGoodMass().var_N_bin, Plots::LLGoodMass().var_min, Plots::LLGoodMass().var_max);
                
                _histo_name="Pzeta_"+_s_process.at(i_proc)+"_"+_s_final_state.at(i_fs)+"_"+_s_category.at(i_cat);
                _histo_labels=";"+Plots::Pzeta().var_X_label + ";" + Plots::Pzeta().var_Y_label;
                histos_1D[i_proc][i_fs][i_cat][2]=new TH1F(_histo_name, _histo_labels, Plots::Pzeta().var_N_bin, Plots::Pzeta().var_min, Plots::Pzeta().var_max);

                _histo_name="MtLMET_"+_s_process.at(i_proc)+"_"+_s_final_state.at(i_fs)+"_"+_s_category.at(i_cat);
                _histo_labels=";"+Plots::MtLMET().var_X_label + ";" + Plots::MtLMET().var_Y_label;
                histos_1D[i_proc][i_fs][i_cat][3]=new TH1F(_histo_name, _histo_labels, Plots::MtLMET().var_N_bin, Plots::MtLMET().var_min, Plots::MtLMET().var_max);

                _histo_name="MtLLMET_"+_s_process.at(i_proc)+"_"+_s_final_state.at(i_fs)+"_"+_s_category.at(i_cat);
                _histo_labels=";"+Plots::MtLLMET().var_X_label + ";" + Plots::MtLLMET().var_Y_label;
                histos_1D[i_proc][i_fs][i_cat][4]=new TH1F(_histo_name, _histo_labels, Plots::MtLLMET().var_N_bin, Plots::MtLLMET().var_min, Plots::MtLLMET().var_max);

                _histo_name="LLPT_"+_s_process.at(i_proc)+"_"+_s_final_state.at(i_fs)+"_"+_s_category.at(i_cat);
                _histo_labels=";"+Plots::LLPT().var_X_label + ";" + Plots::LLPT().var_Y_label;
                histos_1D[i_proc][i_fs][i_cat][5]=new TH1F(_histo_name, _histo_labels, Plots::LLPT().var_N_bin, Plots::LLPT().var_min, Plots::LLPT().var_max);

            }
        }
    }
    else {
        for (int i=i_proc;i<=i_proc+2;i++) {
            for (int i_fs=0;i_fs<=Settings::fs2l;i_fs++) {
                for (int i_cat=0;i_cat<=Settings::catAll;i_cat++) {
                    _histo_name="LLMass_"+_s_process.at(i)+"_"+_s_final_state.at(i_fs)+"_"+_s_category.at(i_cat);
                    _histo_labels=";"+Plots::LLMass().var_X_label + ";" + Plots::LLMass().var_Y_label;
                    histos_1D[i][i_fs][i_cat][0]=new TH1F(_histo_name, _histo_labels, Plots::LLMass().var_N_bin, Plots::LLMass().var_min, Plots::LLMass().var_max);

                    _histo_name="LLGoodMass_"+_s_process.at(i)+"_"+_s_final_state.at(i_fs)+"_"+_s_category.at(i_cat);
                    _histo_labels=";"+Plots::LLGoodMass().var_X_label + ";" + Plots::LLGoodMass().var_Y_label;
                    histos_1D[i][i_fs][i_cat][1]=new TH1F(_histo_name, _histo_labels, Plots::LLGoodMass().var_N_bin, Plots::LLGoodMass().var_min, Plots::LLGoodMass().var_max);

                    _histo_name="Pzeta_"+_s_process.at(i)+"_"+_s_final_state.at(i_fs)+"_"+_s_category.at(i_cat);
                    _histo_labels=";"+Plots::Pzeta().var_X_label + ";" + Plots::Pzeta().var_Y_label;
                    histos_1D[i][i_fs][i_cat][2]=new TH1F(_histo_name, _histo_labels, Plots::Pzeta().var_N_bin, Plots::Pzeta().var_min, Plots::Pzeta().var_max);

                    _histo_name="MtLMET_"+_s_process.at(i)+"_"+_s_final_state.at(i_fs)+"_"+_s_category.at(i_cat);
                    _histo_labels=";"+Plots::MtLMET().var_X_label + ";" + Plots::MtLMET().var_Y_label;
                    histos_1D[i][i_fs][i_cat][3]=new TH1F(_histo_name, _histo_labels, Plots::MtLMET().var_N_bin, Plots::MtLMET().var_min, Plots::MtLMET().var_max);

                    _histo_name="MtLLMET_"+_s_process.at(i)+"_"+_s_final_state.at(i_fs)+"_"+_s_category.at(i_cat);
                    _histo_labels=";"+Plots::MtLLMET().var_X_label + ";" + Plots::MtLLMET().var_Y_label;
                    histos_1D[i][i_fs][i_cat][4]=new TH1F(_histo_name, _histo_labels, Plots::MtLLMET().var_N_bin, Plots::MtLLMET().var_min, Plots::MtLLMET().var_max);
                    
                    _histo_name="LLPT_"+_s_process.at(i)+"_"+_s_final_state.at(i_fs)+"_"+_s_category.at(i_cat);
                    _histo_labels=";"+Plots::LLPT().var_X_label + ";" + Plots::LLPT().var_Y_label;
                    histos_1D[i][i_fs][i_cat][5]=new TH1F(_histo_name, _histo_labels, Plots::LLPT().var_N_bin, Plots::LLPT().var_min, Plots::LLPT().var_max);
                }
            }
        }
    }
}

//===============================================================================
void TreeToHist::SumGroups(int i_proc)
{
    if (i_proc!=Settings::DYee && i_proc!=Settings::EWKZee) {
        for (int i_fs=0;i_fs<Settings::fs2l-1;i_fs++) {
            for (int i_var=0;i_var<num_of_variables;i_var++) {
                for (int i_cat=0;i_cat<Settings::catAll;i_cat++) {
                    histos_1D[i_proc][Settings::fs2l][i_cat][i_var]->Add(histos_1D[i_proc][i_fs][i_cat][i_var]);
                }
            }
        }
        for (int i_fs=1;i_fs<Settings::fs2l-1;i_fs++) {
            for (int i_var=0;i_var<num_of_variables;i_var++) {
                for (int i_cat=0;i_cat<Settings::catAll;i_cat++) {
                    histos_1D[i_proc][Settings::fstautaucomb][i_cat][i_var]->Add(histos_1D[i_proc][i_fs][i_cat][i_var]);
                }
            }
        }
        for (int i_fs=0;i_fs<=Settings::fs2l;i_fs++) {
            for (int i_var=0;i_var<num_of_variables;i_var++) {
                for (int i_cat=0;i_cat<Settings::catAll;i_cat++) {
                    histos_1D[i_proc][i_fs][Settings::catAll][i_var]->Add(histos_1D[i_proc][i_fs][i_cat][i_var]);
                }
            }
        }
    }
    else {
        for (int i=i_proc;i<=i_proc+2;i++) {
            for (int i_fs=0;i_fs<Settings::fs2l-1;i_fs++) {
                for (int i_var=0;i_var<num_of_variables;i_var++) {
                    for (int i_cat=0;i_cat<Settings::catAll;i_cat++) {
                        histos_1D[i][Settings::fs2l][i_cat][i_var]->Add(histos_1D[i][i_fs][i_cat][i_var]);
                    }
                }
            }
            for (int i_fs=1;i_fs<Settings::fs2l-1;i_fs++) {
                for (int i_var=0;i_var<num_of_variables;i_var++) {
                    for (int i_cat=0;i_cat<Settings::catAll;i_cat++) {
                        histos_1D[i][Settings::fstautaucomb][i_cat][i_var]->Add(histos_1D[i][i_fs][i_cat][i_var]);
                    }
                }
            }
            for (int i_fs=0;i_fs<=Settings::fs2l;i_fs++) {
                for (int i_var=0;i_var<num_of_variables;i_var++) {
                    for (int i_cat=0;i_cat<Settings::catAll;i_cat++) {
                        histos_1D[i][i_fs][Settings::catAll][i_var]->Add(histos_1D[i][i_fs][i_cat][i_var]);
                    }
                }
            }
        }
    }
    cout<<"[INFO] Sum groups done"<<endl;
}

//===============================================================================
void TreeToHist::saveHistos(int i_proc)
{
    if (i_proc!=Settings::DYee && i_proc!=Settings::EWKZee) {
        TString output_file_name=_savepath+_directory+"."+_s_process.at(i_proc)+".root";
        output_file = TFile::Open(output_file_name, "recreate");
        output_file->cd();
        for (int i_fs=0;i_fs<=Settings::fs2l;i_fs++) {
            for (int i_var=0;i_var<num_of_variables;i_var++) {
                for (int i_cat=0;i_cat<=Settings::catAll;i_cat++) {
                    histos_1D[i_proc][i_fs][i_cat][i_var]->Write();
                }
            }
        }
        output_file->Close();
        cout<<"[INFO] Saved Histos to "<<output_file_name.Data()<<endl;
    }
    else {
        for (int i=i_proc;i<=i_proc+2;i++) {
            TString output_file_name=_savepath+_directory+"."+_s_process.at(i)+".root";
            output_file = TFile::Open(output_file_name, "recreate");
            output_file->cd();
            for (int i_fs=0;i_fs<=Settings::fs2l;i_fs++) {
                for (int i_var=0;i_var<num_of_variables;i_var++) {
                    for (int i_cat=0;i_cat<=Settings::catAll;i_cat++) {
                        histos_1D[i][i_fs][i_cat][i_var]->Write();
                    }
                }
            }
            output_file->Close();
            cout<<"[INFO] Saved Histos to "<<output_file_name.Data()<<endl;
        }
    }
}

//===============================================================================
void TreeToHist::GetHistos(int i_proc)
{
    string vars[]={"LLMass","LLGoodMass","Pzeta","MtLMET","MtLLMET","LLPT"};
    if (i_proc!=Settings::DYee && i_proc!=Settings::EWKZee) {
        TString output_file_name=_savepath+_directory+"."+_s_process.at(i_proc)+".root";
        output_file = TFile::Open(output_file_name, "read");
        for (int i_fs=0;i_fs<=Settings::fs2l;i_fs++) {
            for (int i_var=0;i_var<num_of_variables;i_var++) {
                for (int i_cat=0;i_cat<=Settings::catAll;i_cat++) {
                    _histo_name=vars[i_var]+"_"+_s_process.at(i_proc)+"_"+_s_final_state.at(i_fs)+"_"+_s_category.at(i_cat);
                    histos_1D[i_proc][i_fs][i_cat][i_var]=(TH1F *)output_file->Get(_histo_name);
                }
            }
        }
        cout<<"[INFO] Get histos from "<<output_file_name<<endl;
    }
    else {
        for (int i=i_proc;i<=i_proc+2;i++) {
            TString output_file_name=_savepath+_directory+"."+_s_process.at(i)+".root";
            output_file = TFile::Open(output_file_name, "read");
            for (int i_fs=0;i_fs<=Settings::fs2l;i_fs++) {
                for (int i_var=0;i_var<num_of_variables;i_var++) {
                    for (int i_cat=0;i_cat<=Settings::catAll;i_cat++) {
                        _histo_name=vars[i_var]+"_"+_s_process.at(i)+"_"+_s_final_state.at(i_fs)+"_"+_s_category.at(i_cat);
                        histos_1D[i][i_fs][i_cat][i_var]=(TH1F *)output_file->Get(_histo_name);
                    }
                }
            }
            cout<<"[INFO] Get histos from "<<output_file_name<<endl;
        }
    }
    //output_file->Close();
} 
//===============================================================================
void TreeToHist::SumTotalMC(bool addSignal)
{
    // string tmp[]={""};
    DeclareHistos(Settings::TotalMC);
    for (int i_fs=0;i_fs<=Settings::fs2l;i_fs++) 
        for (int i_var=0;i_var<num_of_variables;i_var++) 
            for (int i_cat=0;i_cat<=Settings::catAll;i_cat++) 
                for (int i_proc=(addSignal?0:3);i_proc<Settings::TotalMC;i_proc++) {
                    histos_1D[Settings::TotalMC][i_fs][i_cat][i_var]->Add(histos_1D[i_proc][i_fs][i_cat][i_var]);
                }
    cout<<"[INFO] Total MC summed"<<endl;
    saveHistos(Settings::TotalMC);
}

//===============================================================================
void TreeToHist::ToPlots(bool addData, bool addSignal, bool setLog)
{
    string vars[]={"LLMass","LLGoodMass","Pzeta","MtLMET","MtLLMET","LLPT"};
    int colors[]={kRed,kMagenta,kOrange,kCyan,kCyan+1,kYellow,kBlue+2,kBlue-9,kBlue,12,kGreen+3,kGreen-6,kGreen};
    for (int i_fs=0;i_fs<=Settings::fs2l;i_fs++) {
        for (int i_var=0;i_var<num_of_variables;i_var++) {
            for (int i_cat=0;i_cat<=Settings::catAll;i_cat++) {
                TString name=_s_final_state.at(i_fs)+"_"+_s_category.at(i_cat)+"_"+vars[i_var];
                TCanvas *c=new TCanvas(name,name,0,0,1000,800);
                
                c->cd();

                THStack *mc=new THStack(name,name);
                float x=(i_fs==Settings::fsmumu?0.6:0.13);
                TLegend *l=new TLegend( x, 0.5, x+0.35, 0.9 ,ToFSName(_s_final_state.at(i_fs))+", "+_s_category.at(i_cat));
                l->SetFillStyle(0);
                l->SetBorderSize(0);
                l->SetTextSize(0.05);

                for (int i_proc=(addSignal?0:3);i_proc<Settings::TotalMC;i_proc++) {
                    if (histos_1D[i_proc][i_fs][i_cat][i_var]->Integral()>0) {
                        histos_1D[i_proc][i_fs][i_cat][i_var]->SetStats(0);
                        SetColor(histos_1D[i_proc][i_fs][i_cat][i_var],colors[i_proc]);
                        mc->Add(histos_1D[i_proc][i_fs][i_cat][i_var]);
                        if (histos_1D[i_proc][i_fs][i_cat][i_var]->Integral()>histos_1D[Settings::TotalMC][i_fs][i_cat][i_var]->Integral()*0.001 || ((i_proc==Settings::Hmm || i_proc==Settings::Htt) && i_fs==Settings::fsmumu) || (i_proc==Settings::Hww && (i_fs!=Settings::fsetau && i_fs!=Settings::fstautau))) {
                            TString labelName;
                            labelName=ToProcessName(_s_process.at(i_proc));
                            l->AddEntry(histos_1D[i_proc][i_fs][i_cat][i_var],labelName);
                        }
                    }
                }

                if (addData) {
                    histos_1D[Settings::Data][i_fs][i_cat][i_var]->SetBinErrorOption(TH1::kPoisson);
                    histos_1D[Settings::Data][i_fs][i_cat][i_var]->SetStats(0);
                    histos_1D[Settings::Data][i_fs][i_cat][i_var]->SetMarkerSize(1.2);
                    histos_1D[Settings::Data][i_fs][i_cat][i_var]->SetMarkerStyle(8);
                    SetColor(histos_1D[Settings::Data][i_fs][i_cat][i_var],kBlack);
                    l->AddEntry(histos_1D[Settings::Data][i_fs][i_cat][i_var],"Data","lep");
                }

                TH1F *base=(TH1F *)histos_1D[3][i_fs][0][i_var]->Clone();
                base->Reset();
                base->SetStats(0);
                base->GetXaxis()->SetTitleSize(0.05);
                base->GetXaxis()->SetLabelSize(0.045);
                base->GetYaxis()->SetTitleSize(0.05);
                base->GetYaxis()->SetTitleOffset(0.9);
                base->GetYaxis()->SetLabelSize(0.045);
                // base->GetYaxis()->SetTitle("Weighted events");
                base->SetMaximum(histos_1D[addData?(Settings::Data):(Settings::TotalMC)][i_fs][i_cat][i_var]->GetMaximum()*(setLog?5:1.2));
                if (addSignal) {
                    base->SetMinimum(histos_1D[i_fs==Settings::fsmumu?Settings::Hmm : Settings::Htt][i_fs][i_cat][i_var]->GetMaximum()*0.1);
                }

                base->Draw();
                mc->Draw("histsame");
                if (addData)
                    histos_1D[Settings::Data][i_fs][i_cat][i_var]->Draw("Esame");
                l->Draw();
                
                CMS_lumi *lumi = new CMS_lumi;
                lumi->set_lumi(c, _lumi, 0);
                
                c->SetBottomMargin(0.12);
                c->SetGrid();
                if (setLog) c->SetLogy();
                c->SaveAs(_savepath+_directory+".Plots/"+name+".png");
                c->SaveAs(_savepath+_directory+".Plots/"+name+".pdf");
            }
        }
    }
}

//===============================================================================
void TreeToHist::ToPlotsRatio(bool addSignal, bool setLog)
{
    string vars[]={"LLMass","LLGoodMass","Pzeta","MtLMET","MtLLMET","LLPT"};
    int colors[]={kRed,kMagenta,kOrange,kCyan,kCyan+1,kYellow,kBlue+2,kBlue-9,kBlue,12,kGreen+3,kGreen-6,kGreen};
    for (int i_fs=0;i_fs<=Settings::fs2l;i_fs++) {
        for (int i_var=0;i_var<num_of_variables;i_var++) {
            for (int i_cat=0;i_cat<=Settings::catAll;i_cat++) {
                TString name=_s_final_state.at(i_fs)+"_"+_s_category.at(i_cat)+"_"+vars[i_var];
                TCanvas *c=new TCanvas(name,name,0,0,1000,1100);
                c->Divide(1,2);
                TPad *pad1=(TPad*)c->cd(1);
                TPad *pad2=(TPad*)c->cd(2);
                pad1->SetPad(0,0.3,1,1);
                pad2->SetPad(0,0,1,0.3);
                pad1->SetBottomMargin(0);
                pad2->SetTopMargin(0);
                
                pad1->cd();

                THStack *mc=new THStack(name,name);
                float x=(i_fs==Settings::fsmumu?0.6:0.13);
                TLegend *l=new TLegend( x, 0.47, x+0.35, 0.89 ,ToFSName(_s_final_state.at(i_fs))+", "+_s_category.at(i_cat));
                l->SetFillStyle(0);
                l->SetBorderSize(0);
                l->SetTextSize(0.05);

                for (int i_proc=(addSignal?0:3);i_proc<Settings::TotalMC;i_proc++) {
                    if (histos_1D[i_proc][i_fs][i_cat][i_var]->Integral()>0) {
                        histos_1D[i_proc][i_fs][i_cat][i_var]->SetStats(0);
                        SetColor(histos_1D[i_proc][i_fs][i_cat][i_var],colors[i_proc]);
                        mc->Add(histos_1D[i_proc][i_fs][i_cat][i_var]);
                        if (histos_1D[i_proc][i_fs][i_cat][i_var]->Integral()>histos_1D[Settings::TotalMC][i_fs][i_cat][i_var]->Integral()*0.001 || ((i_proc==Settings::Hmm || i_proc==Settings::Htt) && i_fs==Settings::fsmumu) || (i_proc==Settings::Hww && (i_fs!=Settings::fsetau && i_fs!=Settings::fstautau))) {
                            TString labelName;
                            labelName=ToProcessName(_s_process.at(i_proc));
                            l->AddEntry(histos_1D[i_proc][i_fs][i_cat][i_var],labelName);
                        }
                    }
                }

                histos_1D[Settings::Data][i_fs][i_cat][i_var]->SetBinErrorOption(TH1::kPoisson);
                histos_1D[Settings::Data][i_fs][i_cat][i_var]->SetStats(0);
                histos_1D[Settings::Data][i_fs][i_cat][i_var]->SetMarkerSize(1.2);
                histos_1D[Settings::Data][i_fs][i_cat][i_var]->SetMarkerStyle(8);
                SetColor(histos_1D[Settings::Data][i_fs][i_cat][i_var],kBlack);
                l->AddEntry(histos_1D[Settings::Data][i_fs][i_cat][i_var],"Data","lep");

                TH1F *base=(TH1F *)histos_1D[3][i_fs][0][i_var]->Clone();
                base->Reset();
                base->SetStats(0);
                base->GetXaxis()->SetTitleSize(0.05);
                base->GetXaxis()->SetLabelSize(0.045);
                base->GetYaxis()->SetTitleSize(0.05);
                base->GetYaxis()->SetTitleOffset(0.9);
                base->GetYaxis()->SetLabelSize(0.045);
                // base->GetYaxis()->SetTitle("Weighted events");
                base->SetMaximum(histos_1D[Settings::Data][i_fs][i_cat][i_var]->GetMaximum()*(setLog?10:1.2));
                if (addSignal) {
                    base->SetMinimum(histos_1D[i_fs==Settings::fsmumu?Settings::Hmm : Settings::Htt][i_fs][i_cat][i_var]->GetMaximum()*0.1);
                }

                base->Draw();
                mc->Draw("histsame");
                histos_1D[Settings::Data][i_fs][i_cat][i_var]->Draw("Esame");
                l->Draw();

                TH1F *ratio_data=(TH1F*)histos_1D[Settings::Data][i_fs][i_cat][i_var]->Clone();
                ratio_data->Divide(histos_1D[Settings::TotalMC][i_fs][i_cat][i_var]);
                for (int i_bin=1;i_bin<=ratio_data->GetNbinsX();i_bin++) {
                    float err_data=histos_1D[Settings::Data][i_fs][i_cat][i_var]->GetBinError(i_bin);
                    float n_sim=histos_1D[Settings::TotalMC][i_fs][i_cat][i_var]->GetBinContent(i_bin);
                    if (n_sim>0) ratio_data->SetBinError(i_bin,err_data/n_sim);
                }

                TH1F *ratio_sim=(TH1F*)histos_1D[Settings::TotalMC][i_fs][i_cat][i_var]->Clone();
                ratio_sim->Divide(histos_1D[Settings::TotalMC][i_fs][i_cat][i_var]);
                for (int i_bin=1;i_bin<=ratio_sim->GetNbinsX();i_bin++) {
                    float err_sim=histos_1D[Settings::TotalMC][i_fs][i_cat][i_var]->GetBinError(i_bin);
                    float n_sim=histos_1D[Settings::TotalMC][i_fs][i_cat][i_var]->GetBinContent(i_bin);
                    if (n_sim>0) ratio_sim->SetBinError(i_bin,err_sim/n_sim);
                }
                ratio_sim->SetFillColor(kGray);
                ratio_sim->SetLineColor(kGray);
                ratio_sim->SetMarkerColor(kGray);

                TLegend *l_ratio=new TLegend( 0.6, 0.75, 0.9, 0.95);
                l_ratio->AddEntry(ratio_data,"Data");
                l_ratio->AddEntry(ratio_sim,"S+B, stat unc");
                l_ratio->SetFillStyle(0);
                l_ratio->SetBorderSize(0);
                l_ratio->SetTextSize(0.11);

                TH1F *base_ratio=(TH1F *)histos_1D[3][i_fs][0][i_var]->Clone();
                base_ratio->Reset();
                base_ratio->SetStats(0);
                base_ratio->GetXaxis()->SetTitleSize(0.05*7/3);
                base_ratio->GetXaxis()->SetLabelSize(0.045*7/3);
                base_ratio->GetYaxis()->SetTitle("Data/Modeling");
                base_ratio->GetYaxis()->SetTitleSize(0.05*7/3);
                base_ratio->GetYaxis()->SetTitleOffset(0.4);
                base_ratio->GetYaxis()->SetLabelSize(0.045*7/3);
                base_ratio->SetMaximum(1.29);
                base_ratio->SetMinimum(0.71);

                pad2->cd();
                base_ratio->Draw();
                ratio_sim->Draw("E2Same");
                ratio_data->Draw("ESame");
                l_ratio->Draw();
                
                CMS_lumi *lumi = new CMS_lumi;
                lumi->set_lumi(pad1, _lumi, 0);
                
                pad1->SetGrid();
                pad2->SetGrid();
                pad1->SetRightMargin(0.07);
                pad2->SetRightMargin(0.07);
                pad2->SetBottomMargin(0.25);
                if (setLog) pad1->SetLogy();
                c->SetLeftMargin(0.12);
                c->SetBottomMargin(0.25);
                c->SetRightMargin(0);
                gPad->RedrawAxis();
                c->SaveAs(_savepath+_directory+".Plots/"+name+".png");
                // c->SaveAs(_savepath+_directory+".Plots/"+name+".pdf","pdf");
            }
        }
    }
}

//===============================================================================
int TreeToHist::FindFinalState()
{
    int final_state = -999;
    if (abs(LLFlav)==169) final_state=Settings::fsmumu;
    else if (abs(LLFlav)==165) final_state=Settings::fsetau;
    else if (abs(LLFlav)==195) final_state=Settings::fsmutau;
    else if (abs(LLFlav)==225) final_state=Settings::fstautau;
    else if (abs(LLFlav)==143) final_state=Settings::fsemu;
    else final_state=-1;
    return final_state;
}

//-------------------------------------------------------------------------------
bool TreeToHist::pass_Extra_Cuts(int i_fs,float Pzeta,float MtllMET)
{
    bool pass=true;
    if (i_fs==Settings::fsmumu) {
        if (MtllMET>120) pass=true;
    }
    else if (i_fs==Settings::fsetau || i_fs==Settings::fsmutau) {
        if (Pzeta>-70 && Pzeta<40 && MtllMET>80) pass=true;
    }
    else if (i_fs==Settings::fstautau) {
        if (Pzeta>-100 && Pzeta<30 && MtllMET>100) pass=true;
    }
    else if (i_fs==Settings::fsemu) {
        if (Pzeta>-60 && Pzeta<50 && MtllMET>60) pass=true;
    }
    return pass;
}

//===============================================================================
int TreeToHist::FindProcess(int i_proc)
{
    //cout<<GenLLFlav<<endl;
    if (abs(GenLLFlav)==225) return i_proc+1;
    else if (abs(GenLLFlav)==169) return i_proc+2;
    else return i_proc;
}
//===============================================================================
float TreeToHist::calculate_K_factor(TString input_file_name)
{
    float k_factor = 1;
    if ( input_file_name.Contains("ZZTo")) k_factor = KFactor_EW_qqZZ * KFactor_QCD_qqZZ_M;
    else if ( input_file_name.Contains("ggTo")) k_factor = KFactor_QCD_ggZZ_Nominal;
    else if ( input_file_name.Contains("ggH")) k_factor = ggH_NNLOPS_weight;
    return k_factor;
}

//===============================================================================
void TreeToHist::SetColor(TH1F *h,int color)
{
    h->SetMarkerColor(color);
    h->SetLineColor(color);
    h->SetFillColor(color);
}

//===============================================================================
TString TreeToHist::ToFSName(string name)
{
    TString n=name;
    if (n.Contains("comb")) return "#tau#tau, all";
    n.ReplaceAll("mu","#mu");
    n.ReplaceAll("tautau","#tau_{h}#tau_{h}");
    n.ReplaceAll("etau","#tau_{e}#tau_{h}");
    n.ReplaceAll("#mutau","#tau_{#mu}#tau_{h}");
    n.ReplaceAll("emu","#tau_{e}#tau_{#mu}");
    return n;
}

//===============================================================================
TString TreeToHist::ToProcessName(string name)
{
    TString n=name;
    if (n.Contains("comb")) return "#tau#tau, all";
    n.ReplaceAll("mu","#mu");
    n.ReplaceAll("tau","#tau");
    n.ReplaceAll("H","H#rightarrow");
    n.ReplaceAll("DY","DYZ#rightarrow");
    n.ReplaceAll("EWKZ","EWKZ#rightarrow");
    return n;
}

//===============================================================================
int TreeToHist::FindCategory(int njet, bool fine)
{
    int category=-1;
    if (fine) {
        if (njet == 0) {
            category = 0;
        }
        else if (VBFJetIdx1>=0) {
            if ((LLSVPt>0?LLSVPt:LLPt)<=200) category = 1;
            else category = 2;
        }
        else {
            if (njet == 1) category = 3;
            else category = 4;
        }
    }
    else {
        if (njet == 0) category = 0;
        else if (VBFJetIdx1>=0) category = 1;
        else category = 2;
    }
    return category;
}

//===============================================================================
int TreeToHist::NumberOfJets()
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