// Include classes
#include <HTauTauHMuMu/AnalysisStep/test/Plotter/include/TreeToHist.h>

// Constructor
//============================================================
TreeToHist::TreeToHist():Tree()
{    
    _current_process = -999;
    _current_final_state = -999;
    _current_gen_final_state = -999;
    
    _fs_ROS_SS.push_back(0.954003);
    _fs_ROS_SS.push_back(1.05603);
    _fs_ROS_SS.push_back(1.03091);
    _fs_ROS_SS.push_back(0.98326);
    _fs_ROS_SS.push_back(0.985223);
    _fs_ROS_SS.push_back(1.04867);
    _fs_ROS_SS.push_back(1.024);
    _fs_ROS_SS.push_back(0.994852);

    _cb_ss.push_back(1);//(1.112709813);
    _cb_ss.push_back(1);//(1.017112365);
    _cb_ss.push_back(1);//(1.045235335);
    _cb_ss.push_back(1);//(1.021858387);
    _cb_ss.push_back(1);//(1.158447633);
    _cb_ss.push_back(1);//(1.053897151);
    _cb_ss.push_back(1);//(1.054209183);
    _cb_ss.push_back(1);//(1.016327376);
    
    _s_process.push_back("H");
    _s_process.push_back("qqZZ");
    _s_process.push_back("ggZZ");
    _s_process.push_back("rare");
    _s_process.push_back("ZX");
    _s_process.push_back("WZ");
    _s_process.push_back("TotalMC");
    _s_process.push_back("Data");
    
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
    _s_gen_final_state.push_back("background");
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
TreeToHist::~TreeToHist()
{
}
//====================

//===============================================================================
void TreeToHist::SetPaths(string path, string file_name, string savepath)
{
    _path=path;
    _file_name=file_name;
    _savepath=savepath;
}
//===============================================================================
void TreeToHist::ToHistos(string* process, int n_proc, int i_proc, bool Signal)
{
    DeclareHistos(process,n_proc,i_proc);
    
    for (int i=0;i<n_proc;i++) {
        TString input_file_name=_path+process[i]+_file_name;
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
            
            /*bool cut=true;
            for (size_t ilep=0;ilep<4;ilep++)
                if (abs(LepLepId->at(ilep))==15 && (TauDecayMode->at(ilep)==5 || TauDecayMode->at(ilep)==6 || fabs(LepEta->at(ilep))>2.3 || LepPt->at(ilep)<20 || !LepisID->at(ilep))) cut=false;
	    if (!cut) continue;*/
            _current_final_state = FindFinalState();
            if (_current_final_state<0) {continue;}
            if (Signal) _current_gen_final_state = FindGenFinalState();
            else _current_gen_final_state = Settings::gfsbkg;
            _k_factor = calculate_K_factor(input_file_name);
            _TauIDSF = calculate_TauIDSF(input_file_name);
            _event_weight = (_lumi * 1000 * xsec * _k_factor * _TauIDSF * overallEventWeight * L1prefiringWeight) / gen_sum_weights;
            histos_1D[i_proc][_current_final_state][_current_gen_final_state][n_proc>1?i+1:i][0]->Fill(ZZMass,(i_proc == Settings::Data) ? 1 :  _event_weight);
            histos_1D[i_proc][_current_final_state][_current_gen_final_state][n_proc>1?i+1:i][1]->Fill(Z1Mass,(i_proc == Settings::Data) ? 1 :  _event_weight);
            histos_1D[i_proc][_current_final_state][_current_gen_final_state][n_proc>1?i+1:i][2]->Fill(Z2GoodMass,(i_proc == Settings::Data) ? 1 :  _event_weight);
            histos_1D[i_proc][_current_final_state][_current_gen_final_state][n_proc>1?i+1:i][3]->Fill(Z2Pt,(i_proc == Settings::Data) ? 1 :  _event_weight);
        }
        input_file->Close();
        cout<<"[INFO] Processing "<<input_file_name<<" done"<<endl;
    }

    SumGroups(n_proc,i_proc);
    saveHistos(n_proc,i_proc);
}

//===============================================================================
void TreeToHist::ToHistosCR(string* process, int i_proc, bool Signal,TString input_file_FR_name)
{
    if (i_proc!=Settings::ZX) return;
    if (Signal) return;
    
    DeclareHistos(process,1,i_proc);
    
    TString input_file_name=_path+process[0]+_file_name;
    
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
        if ( (!LepisID->at(2)) && abs(LepLepId->at(2))==11 ) {continue;}
        if ( (!LepisID->at(3)) && abs(LepLepId->at(3))==11 ) {continue;}
        if ( (!LepisID->at(2) || LepCombRelIsoPF->at(2) > 0.35) && abs(LepLepId->at(2))==13 ) {continue;}
        if ( (!LepisID->at(3) || LepCombRelIsoPF->at(3) > 0.35) && abs(LepLepId->at(3))==13 ) {continue;}
        if ( ( Lepdz->at(2) > 10 || LepPt->at(2) < 20 || TauDecayMode->at(2)==5 || TauDecayMode->at(2)==6 || !LepisID->at(2)) && abs(LepLepId->at(2))==15 ) {continue;}
        if ( ( Lepdz->at(3) > 10 || LepPt->at(3) < 20 || TauDecayMode->at(3)==5 || TauDecayMode->at(3)==6 || !LepisID->at(3)) && abs(LepLepId->at(3))==15 ) {continue;}
        if ( ZZMass < 70. ) continue;
        
        _current_final_state = FindFinalState();
        if (_current_final_state<0) {continue;}
        _current_gen_final_state = Settings::gfsbkg;
        
        if (abs(Z2Flav)==165) {
            if (abs(LepLepId->at(2))==15 && !(TauVSmu->at(2)>=4 && TauVSe->at(2)>=5 && TauVSjet->at(2)>=5)) continue;
            if (abs(LepLepId->at(3))==15 && !(TauVSmu->at(3)>=4 && TauVSe->at(3)>=5 && TauVSjet->at(3)>=5)) continue;
        }
        else if (abs(Z2Flav)==195 || abs(Z2Flav)==225) {
            if (abs(LepLepId->at(2))==15 && !(TauVSmu->at(2)>=4 && TauVSe->at(2)>=3 && TauVSjet->at(2)>=6)) continue;
            if (abs(LepLepId->at(3))==15 && !(TauVSmu->at(3)>=4 && TauVSe->at(3)>=3 && TauVSjet->at(3)>=6)) continue;
        }
        _yield_SR    = _fs_ROS_SS.at(_current_final_state);
//         _yield_SR    = _fs_ROS_SS.at(_current_final_state)*_cb_ss.at(_current_final_state)*(abs(LepLepId->at(2))==15?1.:FR->GetFakeRate(LepPt->at(2),LepEta->at(2),TauDecayMode->at(2),LepLepId->at(2),0))*(abs(LepLepId->at(3))==15?1.:FR->GetFakeRate(LepPt->at(3),LepEta->at(3),TauDecayMode->at(3),LepLepId->at(3),0));
            
        histos_1D[i_proc][_current_final_state][_current_gen_final_state][0][0]->Fill(ZZMass,_yield_SR);
        histos_1D[i_proc][_current_final_state][_current_gen_final_state][0][1]->Fill(Z1Mass,_yield_SR);
        histos_1D[i_proc][_current_final_state][_current_gen_final_state][0][2]->Fill(Z2GoodMass,_yield_SR);
        histos_1D[i_proc][_current_final_state][_current_gen_final_state][0][3]->Fill(Z2Pt,_yield_SR);
    }
    input_file->Close();
    
    cout<<"[INFO] Processing "<<input_file_name<<" done"<<endl;
    SumGroups(1,i_proc);
    saveHistos(1,i_proc);
}

//===============================================================================
void TreeToHist::ToHistosZXSS(string* process, int i_proc, bool Signal,TString input_file_FR_name)
{
    if (i_proc!=Settings::ZX) return;
    if (Signal) return;
    
    DeclareHistos(process,1,i_proc);
    
    TString input_file_name=_path+process[0]+_file_name;
    
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
        
        histos_1D[i_proc][_current_final_state][_current_gen_final_state][0][0]->Fill(ZZMass,_yield_SR);
        histos_1D[i_proc][_current_final_state][_current_gen_final_state][0][1]->Fill(Z1Mass,_yield_SR);
        histos_1D[i_proc][_current_final_state][_current_gen_final_state][0][2]->Fill(Z2GoodMass,_yield_SR);
        histos_1D[i_proc][_current_final_state][_current_gen_final_state][0][3]->Fill(Z2Pt,_yield_SR);
    }
    input_file->Close();
    
    cout<<"[INFO] Processing "<<input_file_name<<" done"<<endl;
    SumGroups(1,i_proc);
    saveHistos(1,i_proc);
}

//===============================================================================
void TreeToHist::ToHistosZXOS(string* process, int i_proc, bool Signal,TString input_file_FR_name)
{
    if (i_proc!=Settings::ZX) return;
    if (Signal) return;
    
    DeclareHistos(process,1,i_proc);
    
    TString input_file_name=_path+process[0]+_file_name;
    
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
        if (!(test_bit(CRflag, CRZLLos_2P2F)) && !(test_bit(CRflag, CRZLLos_3P1F))) continue;
        
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
        
        if (test_bit(CRflag, CRZLLos_2P2F)) {
            float fr3=FR->GetFakeRate(LepPt->at(2),PFMET,LepEta->at(2),TauDecayMode->at(2),LepLepId->at(2),tauChannel);
            float fr4=FR->GetFakeRate(LepPt->at(3),PFMET,LepEta->at(3),TauDecayMode->at(3),LepLepId->at(3),tauChannel);
            _yield_SR=-1.*fr3*fr4/(1-fr3)/(1-fr4);
        }
        else {
            bool tightLep3=false;
            if (abs(LepLepId->at(3))==11) tightLep3=LepisID->at(3);
            else if (abs(LepLepId->at(3))==13) tightLep3=LepisID->at(3) && LepCombRelIsoPF->at(3) < 0.35;
            else if (abs(LepLepId->at(3))==15) {
                if (tauChannel==0) tightLep3=(LepisID->at(3) && TauVSmu->at(3)>=4 && TauVSe->at(3)>=5 && TauVSjet->at(3)>=5);
                else tightLep3=(LepisID->at(3) && TauVSmu->at(3)>=4 && TauVSe->at(3)>=3 && TauVSjet->at(3)>=6);
            }
            float fr;
            if (tightLep3) fr=FR->GetFakeRate(LepPt->at(2),PFMET,LepEta->at(2),TauDecayMode->at(2),LepLepId->at(2),tauChannel);
            else fr=FR->GetFakeRate(LepPt->at(3),PFMET,LepEta->at(3),TauDecayMode->at(3),LepLepId->at(3),tauChannel);
            _yield_SR=fr/(1-fr);
        }
       
        histos_1D[i_proc][_current_final_state][_current_gen_final_state][0][0]->Fill(ZZMass,_yield_SR);
        histos_1D[i_proc][_current_final_state][_current_gen_final_state][0][1]->Fill(Z1Mass,_yield_SR);
        histos_1D[i_proc][_current_final_state][_current_gen_final_state][0][2]->Fill(Z2GoodMass,_yield_SR);
        histos_1D[i_proc][_current_final_state][_current_gen_final_state][0][3]->Fill(Z2Pt,_yield_SR);
    }
    input_file->Close();
    
    cout<<"[INFO] Processing "<<input_file_name<<" done"<<endl;
    SumGroups(1,i_proc);
    saveHistos(1,i_proc);
}

//===============================================================================
void TreeToHist::DeclareHistos(string* process,int n_proc,int i_proc)
{
    for (int i_fs=0;i_fs<=Settings::fs4l;i_fs++) {
        for (int i_gfs=0;i_gfs<=Settings::gfsall;i_gfs++) {
            if (n_proc>1) {
                for (int i=0;i<n_proc;i++) {
                    _histo_name="MZZ_"+_s_process.at(i_proc)+"_"+process[i]+"_"+_s_final_state.at(i_fs)+"_"+_s_gen_final_state.at(i_gfs);
                    _histo_labels=";"+Plots::M4l().var_X_label + ";" + Plots::M4l().var_Y_label;
                    histos_1D[i_proc][i_fs][i_gfs][i+1][0]=new TH1F(_histo_name, _histo_labels, Plots::M4l().var_N_bin, Plots::M4l().var_min, Plots::M4l().var_max);
                    _histo_name="MZ1_"+_s_process.at(i_proc)+"_"+process[i]+"_"+_s_final_state.at(i_fs)+"_"+_s_gen_final_state.at(i_gfs);
                    _histo_labels=";"+Plots::Z1Mass().var_X_label + ";" + Plots::Z1Mass().var_Y_label;
                    histos_1D[i_proc][i_fs][i_gfs][i+1][1]=new TH1F(_histo_name, _histo_labels, Plots::Z1Mass().var_N_bin, Plots::Z1Mass().var_min, Plots::Z1Mass().var_max);
                    _histo_name="MZ2_"+_s_process.at(i_proc)+"_"+process[i]+"_"+_s_final_state.at(i_fs)+"_"+_s_gen_final_state.at(i_gfs);
                    _histo_labels=";"+Plots::Z2Mass().var_X_label + ";" + Plots::Z2Mass().var_Y_label;
                    histos_1D[i_proc][i_fs][i_gfs][i+1][2]=new TH1F(_histo_name, _histo_labels, Plots::Z2Mass().var_N_bin, Plots::Z2Mass().var_min, Plots::Z2Mass().var_max);
                    _histo_name="PtZ2_"+_s_process.at(i_proc)+"_"+process[i]+"_"+_s_final_state.at(i_fs)+"_"+_s_gen_final_state.at(i_gfs);
                    _histo_labels=";"+Plots::Z2Pt().var_X_label + ";" + Plots::Z2Pt().var_Y_label;
                    histos_1D[i_proc][i_fs][i_gfs][i+1][3]=new TH1F(_histo_name, _histo_labels, Plots::Z2Pt().var_N_bin, Plots::Z2Pt().var_min, Plots::Z2Pt().var_max);
                }
            }
            _histo_name="MZZ_"+_s_process.at(i_proc)+"_"+_s_final_state.at(i_fs)+"_"+_s_gen_final_state.at(i_gfs);
            _histo_labels=";"+Plots::M4l().var_X_label + ";" + Plots::M4l().var_Y_label;
            histos_1D[i_proc][i_fs][i_gfs][0][0]=new TH1F(_histo_name, _histo_labels, Plots::M4l().var_N_bin, Plots::M4l().var_min, Plots::M4l().var_max);
            _histo_name="MZ1_"+_s_process.at(i_proc)+"_"+_s_final_state.at(i_fs)+"_"+_s_gen_final_state.at(i_gfs);
            _histo_labels=";"+Plots::Z1Mass().var_X_label + ";" + Plots::Z1Mass().var_Y_label;
            histos_1D[i_proc][i_fs][i_gfs][0][1]=new TH1F(_histo_name, _histo_labels, Plots::Z1Mass().var_N_bin, Plots::Z1Mass().var_min, Plots::Z1Mass().var_max);
            _histo_name="MZ2_"+_s_process.at(i_proc)+"_"+_s_final_state.at(i_fs)+"_"+_s_gen_final_state.at(i_gfs);
            _histo_labels=";"+Plots::Z2Mass().var_X_label + ";" + Plots::Z2Mass().var_Y_label;
            histos_1D[i_proc][i_fs][i_gfs][0][2]=new TH1F(_histo_name, _histo_labels, Plots::Z2Mass().var_N_bin, Plots::Z2Mass().var_min, Plots::Z2Mass().var_max);
            _histo_name="PtZ2_"+_s_process.at(i_proc)+"_"+_s_final_state.at(i_fs)+"_"+_s_gen_final_state.at(i_gfs);
            _histo_labels=";"+Plots::Z2Pt().var_X_label + ";" + Plots::Z2Pt().var_Y_label;
            histos_1D[i_proc][i_fs][i_gfs][0][3]=new TH1F(_histo_name, _histo_labels, Plots::Z2Pt().var_N_bin, Plots::Z2Pt().var_min, Plots::Z2Pt().var_max);
        }
    }
}

//===============================================================================
void TreeToHist::SumGroups(int n_proc,int i_proc)
{
    if (n_proc>1) {
        for (int i_fs=0;i_fs<Settings::fs4l;i_fs++)
            for (int i_gfs=0;i_gfs<Settings::gfsall;i_gfs++)
                for (int i=0;i<n_proc;i++)
                    for (int i_var=0;i_var<4;i_var++)
                        histos_1D[i_proc][i_fs][i_gfs][0][i_var]->Add(histos_1D[i_proc][i_fs][i_gfs][i+1][i_var]); 
    }
    
    for (int i_fs=0;i_fs<Settings::fs4l;i_fs++)
        for (int i_gfs=0;i_gfs<Settings::gfsall;i_gfs++)
            for (int i_var=0;i_var<4;i_var++) {
                histos_1D[i_proc][i_fs][Settings::gfsall][0][i_var]->Add(histos_1D[i_proc][i_fs][i_gfs][0][i_var]);
                histos_1D[i_proc][Settings::fs4l][i_gfs][0][i_var]->Add(histos_1D[i_proc][i_fs][i_gfs][0][i_var]);
                histos_1D[i_proc][Settings::fs4l][Settings::gfsall][0][i_var]->Add(histos_1D[i_proc][i_fs][i_gfs][0][i_var]);
            }
    cout<<"[INFO] Sum groups done"<<endl;
}

//===============================================================================
void TreeToHist::saveHistos(int n_proc,int i_proc)
{
    TString output_file_name=_savepath+_s_process.at(i_proc)+".root";
    output_file = TFile::Open(output_file_name, "recreate");
    output_file->cd();
    if (n_proc>1) {
        for (int i_fs=0;i_fs<Settings::fs4l;i_fs++)
            for (int i_gfs=0;i_gfs<Settings::gfsall;i_gfs++)
                for (int i=0;i<n_proc;i++)
                    for (int i_var=0;i_var<4;i_var++)
                        histos_1D[i_proc][i_fs][i_gfs][i+1][i_var]->Write();
    }
    for (int i_fs=0;i_fs<=Settings::fs4l;i_fs++) 
        for (int i_gfs=0;i_gfs<=Settings::gfsall;i_gfs++) 
            for (int i_var=0;i_var<4;i_var++) 
                histos_1D[i_proc][i_fs][i_gfs][0][i_var]->Write();
    output_file->Close();
    cout<<"[INFO] Saved Histos to "<<output_file_name.Data()<<endl;
}

//===============================================================================
void TreeToHist::GetHistos(string* process, int n_proc, int i_proc, bool Signal)
{
    TString output_file_name=_savepath+_s_process.at(i_proc)+".root";
    output_file = TFile::Open(output_file_name, "read");
    string vars[]={"MZZ","MZ1","MZ2","PtZ2"};
    if (n_proc>1) {
        for (int i_fs=0;i_fs<Settings::fs4l;i_fs++)
            for (int i_gfs=0;i_gfs<Settings::gfsall;i_gfs++)
                for (int i=0;i<n_proc;i++)
                    for (int i_var=0;i_var<4;i_var++) {
                        _histo_name=vars[i_var]+"_"+_s_process.at(i_proc)+"_"+process[i]+"_"+_s_final_state.at(i_fs)+"_"+_s_gen_final_state.at(i_gfs);
                        histos_1D[i_proc][i_fs][i_gfs][i+1][i_var]=(TH1F *)output_file->Get(_histo_name);
                    }
    }
    for (int i_fs=0;i_fs<=Settings::fs4l;i_fs++) 
        for (int i_gfs=0;i_gfs<=Settings::gfsall;i_gfs++) 
            for (int i_var=0;i_var<4;i_var++) {
                _histo_name=vars[i_var]+"_"+_s_process.at(i_proc)+"_"+_s_final_state.at(i_fs)+"_"+_s_gen_final_state.at(i_gfs);
                histos_1D[i_proc][i_fs][i_gfs][0][i_var]=(TH1F *)output_file->Get(_histo_name);
            }
    cout<<"[INFO] Get histos from "<<output_file_name<<endl;
    //output_file->Close();
} 
//===============================================================================
void TreeToHist::SumTotalMC()
{
    string tmp[]={""};
    DeclareHistos(tmp,1,Settings::TotalMC);
    for (int i_fs=0;i_fs<=Settings::fs4l;i_fs++) 
        for (int i_gfs=0;i_gfs<=Settings::gfsall;i_gfs++)
            for (int i_var=0;i_var<4;i_var++) 
                for (int i_proc=0;i_proc<Settings::TotalMC;i_proc++)
                    histos_1D[Settings::TotalMC][i_fs][i_gfs][0][i_var]->Add(histos_1D[i_proc][i_fs][i_gfs][0][i_var]);
    cout<<"[INFO] Total MC summed"<<endl;
    saveHistos(1,Settings::TotalMC);
}

//===============================================================================
void TreeToHist::ToPlots(bool addData)
{
    string vars[]={"MZZ","MZ1","MZ2","PtZ2"};
    int colors[]={kRed,kMagenta,kBlue,kCyan,kGreen};
    for (int i_fs=0;i_fs<=Settings::fs4l;i_fs++) {
        for (int i_var=0;i_var<4;i_var++) {
            TString name=_s_final_state.at(i_fs)+"_"+vars[i_var];
            TCanvas *c=new TCanvas(name,name,0,0,1000,800);
            
            THStack *mc=new THStack(name,name);
            float x=((i_fs==0 || i_fs==4)?0.2:0.55);
            TLegend *l=new TLegend( x, 0.4, x+0.35, 0.9 ,ToFSName(_s_final_state.at(i_fs)));
            l->SetFillStyle(0);
            l->SetBorderSize(0);
            l->SetTextSize(0.05);
            vector<TH1F *> HList;
            
            for (int i_proc=0;i_proc<Settings::rare;i_proc++)
                for (int i_gfs=0;i_gfs<=Settings::gfsbkg;i_gfs++)
                    if (histos_1D[i_proc][i_fs][i_gfs][0][i_var]->Integral()>0) {
                        histos_1D[i_proc][i_fs][i_gfs][0][i_var]->SetStats(0);
                        SetColor(histos_1D[i_proc][i_fs][i_gfs][0][i_var],colors[i_proc]-i_gfs*2);
                        HList.push_back(histos_1D[i_proc][i_fs][i_gfs][0][i_var]);
                        if (histos_1D[i_proc][i_fs][i_gfs][0][i_var]->Integral()>histos_1D[Settings::Data][i_fs][Settings::gfsall][0][i_var]->Integral()*0.001) {
                            TString labelName;
                            if (i_gfs==Settings::gfsbkg) labelName=_s_process.at(i_proc)+"#rightarrowOthers";
                            else labelName=_s_process.at(i_proc)+"#rightarrow"+ToGFSName(_s_gen_final_state.at(i_gfs));
                            l->AddEntry(histos_1D[i_proc][i_fs][i_gfs][0][i_var],labelName);
                        }
                    }

            for (int i_proc=Settings::rare;i_proc<=Settings::ZX;i_proc++)
                if (histos_1D[i_proc][i_fs][Settings::gfsbkg][0][i_var]->Integral()>0) {
                    SetColor(histos_1D[i_proc][i_fs][Settings::gfsbkg][0][i_var],colors[i_proc]);
                    HList.push_back(histos_1D[i_proc][i_fs][Settings::gfsbkg][0][i_var]);
                    TString labelName=_s_process.at(i_proc);
                    l->AddEntry(histos_1D[i_proc][i_fs][Settings::gfsbkg][0][i_var],labelName);
                }

            for (int i=HList.size()-1;i>=0;i--)
                mc->Add(HList.at((size_t)i));
            
            if (addData) {
                histos_1D[Settings::Data][i_fs][Settings::gfsall][0][i_var]->SetBinErrorOption(TH1::kPoisson);
                histos_1D[Settings::Data][i_fs][Settings::gfsall][0][i_var]->SetBinErrorOption(TH1::kPoisson);
                histos_1D[Settings::Data][i_fs][Settings::gfsall][0][i_var]->SetStats(0);
                histos_1D[Settings::Data][i_fs][Settings::gfsall][0][i_var]->SetMarkerSize(1.2);
                histos_1D[Settings::Data][i_fs][Settings::gfsall][0][i_var]->SetMarkerStyle(8);
                SetColor(histos_1D[Settings::Data][i_fs][Settings::gfsall][0][i_var],kBlack);
                l->AddEntry(histos_1D[Settings::Data][i_fs][Settings::gfsall][0][i_var],"Data","lep");
            }

            TH1F *base=(TH1F *)histos_1D[Settings::H][i_fs][0][0][i_var]->Clone();
            base->Reset();
            base->SetStats(0);
            base->GetXaxis()->SetTitleSize(0.05);
            base->GetXaxis()->SetLabelSize(0.045);
            base->GetYaxis()->SetTitleSize(0.05);
            base->GetYaxis()->SetTitleOffset(0.9);
            base->GetYaxis()->SetLabelSize(0.045);
            if (addData)
                base->SetMaximum(max(histos_1D[Settings::Data][i_fs][Settings::gfsall][0][i_var]->GetMaximum(),histos_1D[Settings::TotalMC][i_fs][Settings::gfsall][0][i_var]->GetMaximum())*1.2);
            else
                base->SetMaximum(histos_1D[Settings::TotalMC][i_fs][Settings::gfsall][0][i_var]->GetMaximum()*1.3);

            base->Draw();
            mc->Draw("histsame");
            if (addData)
                histos_1D[Settings::Data][i_fs][Settings::gfsall][0][i_var]->Draw("Esame");
            l->Draw();
            
            CMS_lumi *lumi = new CMS_lumi;
            lumi->set_lumi(c, _lumi, 0);
            
            c->SetBottomMargin(0.12);
            c->SetGrid();
            c->SaveAs(_savepath+"Plots/"+name+(addData?"data.png":".png"));
            c->SaveAs(_savepath+"Plots/"+name+(addData?"data.png":".pdf"));
        }
    }
}

//===============================================================================
void TreeToHist::ToPlotsMC(bool addData)
{
    string vars[]={"MZZ","MZ1","MZ2","PtZ2"};
    int colors[]={kRed,kMagenta,kBlue,kCyan,kGreen,kYellow};
    for (int i_fs=0;i_fs<=Settings::fs4l;i_fs++) {
        for (int i_var=0;i_var<4;i_var++) {
            TString name=_s_final_state.at(i_fs)+"_"+vars[i_var];
            TCanvas *c=new TCanvas(name,name,0,0,1000,800);
            
            THStack *mc=new THStack(name,name);
            float x=((i_fs==0 || i_fs==4)?0.2:0.55);
            TLegend *l=new TLegend( x, 0.4, x+0.35, 0.9 ,ToFSName(_s_final_state.at(i_fs)));
            l->SetFillStyle(0);
            l->SetBorderSize(0);
            l->SetTextSize(0.05);
            vector<TH1F *> HList;
            
            for (int i_proc=0;i_proc<Settings::rare;i_proc++)
                for (int i_gfs=0;i_gfs<=Settings::gfsbkg;i_gfs++)
                    if (histos_1D[i_proc][i_fs][i_gfs][0][i_var]->Integral()>0) {
                        histos_1D[i_proc][i_fs][i_gfs][0][i_var]->SetStats(0);
                        SetColor(histos_1D[i_proc][i_fs][i_gfs][0][i_var],colors[i_proc]-i_gfs*2);
                        HList.push_back(histos_1D[i_proc][i_fs][i_gfs][0][i_var]);
                        if (histos_1D[i_proc][i_fs][i_gfs][0][i_var]->Integral()>histos_1D[Settings::Data][i_fs][Settings::gfsall][0][i_var]->Integral()*0.001) {
                            TString labelName;
                            if (i_gfs==Settings::gfsbkg) labelName=_s_process.at(i_proc)+"#rightarrowOthers";
                            else labelName=_s_process.at(i_proc)+"#rightarrow"+ToGFSName(_s_gen_final_state.at(i_gfs));
                            l->AddEntry(histos_1D[i_proc][i_fs][i_gfs][0][i_var],labelName);
                        }
                    }

            for (int i_proc=Settings::rare;i_proc<=Settings::WZ;i_proc++)
                if (histos_1D[i_proc][i_fs][Settings::gfsbkg][0][i_var]->Integral()>0) {
                    SetColor(histos_1D[i_proc][i_fs][Settings::gfsbkg][0][i_var],colors[i_proc]);
                    HList.push_back(histos_1D[i_proc][i_fs][Settings::gfsbkg][0][i_var]);
                    TString labelName=_s_process.at(i_proc);
                    l->AddEntry(histos_1D[i_proc][i_fs][Settings::gfsbkg][0][i_var],labelName);
                }

            for (int i=HList.size()-1;i>=0;i--)
                mc->Add(HList.at((size_t)i));
            
            if (addData) {
                histos_1D[Settings::Data][i_fs][Settings::gfsall][0][i_var]->SetBinErrorOption(TH1::kPoisson);
                histos_1D[Settings::Data][i_fs][Settings::gfsall][0][i_var]->SetStats(0);
                histos_1D[Settings::Data][i_fs][Settings::gfsall][0][i_var]->SetMarkerSize(1.2);
                histos_1D[Settings::Data][i_fs][Settings::gfsall][0][i_var]->SetMarkerStyle(8);
                SetColor(histos_1D[Settings::Data][i_fs][Settings::gfsall][0][i_var],kBlack);
                l->AddEntry(histos_1D[Settings::Data][i_fs][Settings::gfsall][0][i_var],"Data","lep");
            }

            TH1F *base=(TH1F *)histos_1D[Settings::H][i_fs][0][0][i_var]->Clone();
            base->Reset();
            base->SetStats(0);
            base->GetXaxis()->SetTitleSize(0.05);
            base->GetXaxis()->SetLabelSize(0.045);
            base->GetYaxis()->SetTitleSize(0.05);
            base->GetYaxis()->SetTitleOffset(0.9);
            base->GetYaxis()->SetLabelSize(0.045);
            if (addData)
                base->SetMaximum(max(histos_1D[Settings::Data][i_fs][Settings::gfsall][0][i_var]->GetMaximum(),mc->GetMaximum())*1.3);
            else
                base->SetMaximum(histos_1D[Settings::TotalMC][i_fs][Settings::gfsall][0][i_var]->GetMaximum()*1.3);

            base->Draw();
            mc->Draw("histsame");
            if (addData)
                histos_1D[Settings::Data][i_fs][Settings::gfsall][0][i_var]->Draw("Esame");
            l->Draw();
            
            CMS_lumi *lumi = new CMS_lumi;
            lumi->set_lumi(c, _lumi, 0);
            
            c->SetBottomMargin(0.12);
            c->SetGrid();
            c->SaveAs(_savepath+"Plots/"+name+(addData?"data.png":".png"));
            c->SaveAs(_savepath+"Plots/"+name+(addData?"data.png":".pdf"));
        }
    }
}

//===============================================================================
int TreeToHist::FindFinalState()
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
int TreeToHist::FindGenFinalState()
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
float TreeToHist::calculate_K_factor(TString input_file_name)
{
    float k_factor = 1;
    if ( input_file_name.Contains("ZZTo4l")) k_factor = KFactor_EW_qqZZ * KFactor_QCD_qqZZ_M;
    else if ( input_file_name.Contains("ggTo")) k_factor = KFactor_QCD_ggZZ_Nominal;
    else if ( input_file_name.Contains("ggH")) k_factor = ggH_NNLOPS_weight;
    return k_factor;

}

//===============================================================================
int TreeToHist::GENMatch(int ileg,TString input_file_name)
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

//===============================================================================
float TreeToHist::calculate_TauIDSF_OneLeg(short DM, float pt, float eta, int genmatch, int flav)
{
    eta=fabs(eta);
    if (genmatch==6)
        return 1.;
    else if (genmatch==5) {
        if (flav==165) return DeepTauSF_VSjet_ETau->getSFvsPT(pt,genmatch);
        else if (flav==195) return DeepTauSF_VSjet_MuTau->getSFvsPT(pt,genmatch);
        else if (flav==225) return DeepTauSF_VSjet_TauTau->getSFvsPT(pt,genmatch);
        else return 1.;
    }
    else if (genmatch==1 || genmatch==3) {
        if (flav==165) return DeepTauSF_VSe_ETau->getSFvsEta(eta,genmatch);
        else if (flav==195) return DeepTauSF_VSe_MuTau->getSFvsEta(eta,genmatch);
        else if (flav==225) return DeepTauSF_VSe_TauTau->getSFvsEta(eta,genmatch);
        else return 1.;
    }
    else if (genmatch==2 || genmatch==4) {
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

//===============================================================================
float TreeToHist::calculate_TauIDSF(TString input_file_name)
{
    float SF = 1;
    for (size_t ileg=2;ileg<=3;ileg++) {
        if (abs(LepLepId->at(ileg))!=15)
            continue;
        int genmatch=GENMatch(ileg,input_file_name);
        SF=SF*calculate_TauIDSF_OneLeg(TauDecayMode->at(ileg),LepPt->at(ileg),LepEta->at(ileg),genmatch,abs(Z2Flav));
    }
    return SF;
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
    n.ReplaceAll("mu","#mu");
    n.ReplaceAll("tautau","#tau_{h}#tau_{h}");
    n.ReplaceAll("etau","#tau_{e}#tau_{h}");
    n.ReplaceAll("#mutau","#tau_{#mu}#tau_{h}");
    return n;
}

//===============================================================================
TString TreeToHist::ToGFSName(string name)
{
    TString n=name;
    n.ReplaceAll("mu","#mu");
    n.ReplaceAll("tau","#tau");
    return n;
}
