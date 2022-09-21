// Include classes
#include <HTauTauHMuMu/AnalysisStep/test/ZpXEstimation/include/OSmethod.h>

// Constructor
//============================================================
OSmethod::OSmethod():Tree()
{
   _current_process = -999;
   _current_final_state = -999;
   _current_category = -999;
   _current_category_stxs = -999;
   
   _s_process.push_back("Data");
   _s_process.push_back("WZ");
   _s_process.push_back("qqZZ");
   _s_process.push_back("DY");
   _s_process.push_back("ttbar");
   
   _s_flavour.push_back("ele");
   _s_flavour.push_back("mu");
   _s_flavour.push_back("tauE");
   _s_flavour.push_back("tauMu");
   _s_flavour.push_back("tauTau");
   
   _s_final_state.push_back("4mu");
   _s_final_state.push_back("2muetau");
   _s_final_state.push_back("2mumutau");
   _s_final_state.push_back("2mutautau");
   _s_final_state.push_back("2e2mu");
   _s_final_state.push_back("2eetau");
   _s_final_state.push_back("2emutau");
   _s_final_state.push_back("2etautau");
   _s_final_state.push_back("4l");
  /* 
   _s_category.push_back("UnTagged");
   _s_category.push_back("VBF1jTagged");
   _s_category.push_back("VBF2jTagged");
   _s_category.push_back("VHLeptTagged");
   _s_category.push_back("VHHadrTagged");
   _s_category.push_back("ttHLeptTagged");
   _s_category.push_back("ttHHadrTagged");
   _s_category.push_back("VHMETTagged");
   _s_category.push_back("Inclusive");
   
   _s_category_stxs.push_back("ggH_0J_PTH_0_10");
   _s_category_stxs.push_back("ggH_0J_PTH_10_200");
   _s_category_stxs.push_back("ggH_1J_PTH_0_60");
   _s_category_stxs.push_back("ggH_1J_PTH_60_120");
   _s_category_stxs.push_back("ggH_1J_PTH_120_200");
   _s_category_stxs.push_back("ggH_2J_PTH_0_60");
   _s_category_stxs.push_back("ggH_2J_PTH_60_120");
   _s_category_stxs.push_back("ggH_2J_PTH_120_200");
   _s_category_stxs.push_back("ggH_PTH_200");
   _s_category_stxs.push_back("ggH_VBF");
   _s_category_stxs.push_back("VBF_1j");
   _s_category_stxs.push_back("VBF_2j");
   _s_category_stxs.push_back("VBF_2j_mjj_350_700_2j");
   _s_category_stxs.push_back("VBF_2j_mjj_GT700_2j");
   _s_category_stxs.push_back("VBF_2j_mjj_GT350_3j");
   _s_category_stxs.push_back("VBF_GT200_2J");
   _s_category_stxs.push_back("VH_Had");
   _s_category_stxs.push_back("VBF_rest_VH");
   _s_category_stxs.push_back("VH_lep_0_150");
   _s_category_stxs.push_back("VH_Lep_GT150");
   _s_category_stxs.push_back("ttH_Lep");
   _s_category_stxs.push_back("ttH_Had");
   _s_category_stxs.push_back("Inclusive");
   */
   _s_region.push_back("2P2F");
   _s_region.push_back("3P1F");
   _s_region.push_back("OS");
	
   _s_variation.push_back("nominal");
   _s_variation.push_back("Up");
   _s_variation.push_back("Dn");
  
   std::string year="2016Legacy";

   DeepTauSF_VSe_ETau = new TauIDSFTool(year,"DeepTau2017v2p1VSe","Medium");
   DeepTauSF_VSmu_ETau = new TauIDSFTool(year,"DeepTau2017v2p1VSmu","Tight");
   DeepTauSF_VSjet_ETau = new TauIDSFTool(year,"DeepTau2017v2p1VSjet","Medium");
   DeepTauSF_VSe_MuTau = new TauIDSFTool(year,"DeepTau2017v2p1VSe","VLoose");
   DeepTauSF_VSmu_MuTau = new TauIDSFTool(year,"DeepTau2017v2p1VSmu","Tight");
   DeepTauSF_VSjet_MuTau = new TauIDSFTool(year,"DeepTau2017v2p1VSjet","Tight");
   DeepTauSF_VSe_TauTau = new TauIDSFTool(year,"DeepTau2017v2p1VSe","VLoose");
   DeepTauSF_VSmu_TauTau = new TauIDSFTool(year,"DeepTau2017v2p1VSmu","Tight");
   DeepTauSF_VSjet_TauTau = new TauIDSFTool(year,"DeepTau2017v2p1VSjet","Tight");
 
   DeclareFRHistos();
   DeclareDataMCHistos();
   DeclareZXHistos();
}
//============================================================



// Destructor
//====================
OSmethod::~OSmethod()
{
}
//====================


//===============================================================================
void OSmethod::FillFRHistos( TString input_file_data_name )
{
   input_file_data = TFile::Open( input_file_data_name);
   
   hCounters = (TH1F*)input_file_data->Get("CRZLTree/Counters");
   gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);
   
   input_tree_data = (TTree*)input_file_data->Get("CRZLTree/candTree");
   Init( input_tree_data, input_file_data_name , false);
   
   _current_process = find_current_process(input_file_data_name);
   
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   cout<<"[INFO] Processing "<<input_file_data_name<<", total event: "<<nentries<<endl;
   Long64_t nbytes = 0, nb = 0;
	
	// Define some counters for control print out
	Int_t _total_events = 0;
	Int_t _failZ1MassCut = 0;
	Int_t _failLepPtCut = 0;
	Int_t _failEtaCut = 0;
	Int_t _failSipVtxCut = 0;
	Int_t _failMETCut = 0;
	Int_t _passingSelection = 0;
	Int_t _faillingSelection = 0;
	Int_t _faillingJPsiMassCut = 0;
//tau
	Int_t _passingSelectionMu = 0;
        Int_t _passingSelectionE = 0;
        Int_t _passingSelectionTau = 0;
        Int_t _faillingSelectionMu = 0;
        Int_t _faillingSelectionE = 0;
        Int_t _faillingSelectionTau = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry%50000==0) cout<<ientry<<endl;
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
			   
      _total_events++;
		
      TLorentzVector p1,p2,p3;
      p1.SetPtEtaPhiM(LepPt->at(0), LepEta->at(0), LepPhi->at(0), 0.);
      p2.SetPtEtaPhiM(LepPt->at(1), LepEta->at(1), LepPhi->at(1), 0.);
      p3.SetPtEtaPhiM(LepPt->at(2), LepEta->at(2), LepPhi->at(2), 0.);
	   
      if ( fabs(Z1Mass - 91.2) > 7. ) {_failZ1MassCut++; continue;}
      if ( (LepPt->at(0) > LepPt->at(1)) && (LepPt->at(0) < 20. || LepPt->at(1) < 10.) ) {_failLepPtCut++; continue;}
      if ( (LepPt->at(1) > LepPt->at(0)) && (LepPt->at(1) < 20. || LepPt->at(0) < 10.) ) {_failLepPtCut++; continue;}
      if ( (LepPt->at(2) < 20. && abs(LepLepId->at(2)) == 15 ) ) {_failLepPtCut++; continue;}
      if ( (fabs(LepEta->at(2)) > 2.5 ) && ( abs(LepLepId->at(2)) == 11 || abs(LepLepId->at(2)) == 13 )) {_failEtaCut++; continue;}
      if ( (fabs(LepEta->at(2)) > 2.3 ) && ( abs(LepLepId->at(2)) == 15 )) {_failEtaCut++; continue;}
      if ( (LepSIP->at(2) > 4. || Lepdxy->at(2) > 0.5 || Lepdz->at(2) > 1.0) && (abs(LepLepId->at(2)) == 11)) { _failSipVtxCut++; continue;} // Included dxy/dz cuts for ele       
      if ( (LepSIP->at(2) > 4. || Lepdxy->at(2) > 0.5 || Lepdz->at(2) > 1.0) && (abs(LepLepId->at(2)) == 13)) { _failSipVtxCut++; continue;} // Included dxy/dz cuts for mu   
      if ( (Lepdz->at(2) > 10) && (abs(LepLepId->at(2)) == 15)) { _failSipVtxCut++; continue;}
      // NB: Included SIP cut on muons that was removed when it was included in the muon BDT
      if ( PFMET > 25. ) {_failMETCut++; continue;}
      if ( (LepLepId->at(2) < 0 && LepLepId->at(0) > 0 && (p1+p3).M() < 4.) || (LepLepId->at(2) < 0 && LepLepId->at(1) > 0 && (p2+p3).M() < 4.) ) {_faillingJPsiMassCut++; continue;}
      if ( (LepLepId->at(2) > 0 && LepLepId->at(0) < 0 && (p1+p3).M() < 4.) || (LepLepId->at(2) > 0 && LepLepId->at(1) < 0 && (p2+p3).M() < 4.) ) {_faillingJPsiMassCut++; continue;}
      else
      {
         // Final event weight
         _k_factor = calculate_K_factor(input_file_data_name);
	 //_TauIDSF = calculate_TauIDSF(input_file_data_name);
         _event_weight = (_lumi * 1000 * xsec * _k_factor * overallEventWeight * L1prefiringWeight) / gen_sum_weights;

         //if( LepisID->at(2) ) // Changed because we are not using BDT-based muon ID but PF+ISO 
         if(abs(LepLepId->at(2)) == 11 || abs(LepLepId->at(2)) == 13) {           
         if(LepisID->at(2) && ((abs(LepLepId->at(2)) == 13) ? LepCombRelIsoPF->at(2) < 0.35 : LepCombRelIsoPF->at(2) < 999999.))
         {
	   _passingSelection++;
	   if(abs(LepLepId->at(2)) == 11 ) passing[_current_process][Settings::ele]->Fill(LepPt->at(2), (fabs(LepEta->at(2)) < 1.479) ? 0.5 : 1.5 , (_current_process == Settings::Data) ? 1 :  _event_weight);
	   else if(abs(LepLepId->at(2)) == 13 ) passing[_current_process][Settings::mu]->Fill(LepPt->at(2), (fabs(LepEta->at(2)) < 1.2) ? 0.5 : 1.5 , (_current_process == Settings::Data) ? 1 :  _event_weight);
         }
         else
	   {
	     _faillingSelection++;
	     if(abs(LepLepId->at(2)) == 11 ) failing[_current_process][Settings::ele]->Fill(LepPt->at(2), (fabs(LepEta->at(2)) < 1.479) ? 0.5 : 1.5 , (_current_process == Settings::Data) ? 1 :  _event_weight);
	     else if(abs(LepLepId->at(2)) == 13 ) failing[_current_process][Settings::mu]->Fill(LepPt->at(2), (fabs(LepEta->at(2)) < 1.2) ? 0.5 : 1.5 , (_current_process == Settings::Data) ? 1 :  _event_weight);
	   }
	 }
	 else {
	 if (LepisID->at(2) && TauVSmu->at(2)>=4 && TauVSe->at(2)>=3 && TauVSjet->at(2)>=6) {
	     _passingSelectionMu++;
	     _passingSelectionTau++;
	     //TString WP[3]={"Tight","VLoose","Tight"};
	     _TauIDSF=calculate_TauIDSF_OneLeg(TauDecayMode->at(2), LepPt->at(2), LepEta->at(2), GENMatch(2,input_file_data_name), abs(Z2Flav));
	     passing[_current_process][Settings::tauMu]->Fill(LepPt->at(2), (fabs(LepEta->at(2)) < 1.479) ? 0.5 : 1.5 , (_current_process == Settings::Data) ? 1 :  _event_weight*_TauIDSF);
	     passing[_current_process][Settings::tauTau]->Fill(LepPt->at(2), (fabs(LepEta->at(2)) < 1.479) ? 0.5 : 1.5 , (_current_process == Settings::Data) ? 1 :  _event_weight*_TauIDSF);
	 }
	 else {
	     _faillingSelectionMu++;
	     _faillingSelectionTau++;
	     //TString WP[3]={"Tight","VLoose","Tight"};
             _TauIDSF=calculate_TauIDSF_OneLeg(TauDecayMode->at(2), LepPt->at(2), LepEta->at(2), GENMatch(2,input_file_data_name), abs(Z2Flav));
	     failing[_current_process][Settings::tauMu]->Fill(LepPt->at(2), (fabs(LepEta->at(2)) < 1.479) ? 0.5 : 1.5 , (_current_process == Settings::Data) ? 1 :  _event_weight*_TauIDSF);
	     failing[_current_process][Settings::tauTau]->Fill(LepPt->at(2), (fabs(LepEta->at(2)) < 1.479) ? 0.5 : 1.5 , (_current_process == Settings::Data) ? 1 :  _event_weight*_TauIDSF);
	 }
	 if (LepisID->at(2) && TauVSmu->at(2)>=4 && TauVSe->at(2)>=5 && TauVSjet->at(2)>=5) {
	    _passingSelectionE++;
	    //TString WP[3]={"Medium","Medium","Tight"};
	    _TauIDSF=calculate_TauIDSF_OneLeg(TauDecayMode->at(2), LepPt->at(2), LepEta->at(2), GENMatch(2,input_file_data_name), abs(Z2Flav));
	    passing[_current_process][Settings::tauE]->Fill(LepPt->at(2), (fabs(LepEta->at(2)) < 1.479) ? 0.5 : 1.5 , (_current_process == Settings::Data) ? 1 :  _event_weight*_TauIDSF);
	 }
	 else {
	    _faillingSelectionE++;
	    //TString WP[3]={"Medium","Medium","Tight"};
            _TauIDSF=calculate_TauIDSF_OneLeg(TauDecayMode->at(2), LepPt->at(2), LepEta->at(2), GENMatch(2,input_file_data_name), abs(Z2Flav));
	    failing[_current_process][Settings::tauE]->Fill(LepPt->at(2), (fabs(LepEta->at(2)) < 1.479) ? 0.5 : 1.5 , (_current_process == Settings::Data) ? 1 :  _event_weight*_TauIDSF);
	 }
	 }
      }
	
   } // END events loop
	
	// OS method: control printout for ele/mu in Z+L CR 
	if( _current_process == Settings::Data)
	{
		cout << endl;
		cout << "========================================================================================" << endl;
		cout << "[INFO] Control printout for Z+L control region." << endl;
		cout << "========================================================================================" << endl;
		cout << "[INFO] Total number of events in Z+L control region = " << _total_events << endl;
		cout << "[INFO] Events lost after  abs(Z1 - Z) < 7 GeV cut  = " << _failZ1MassCut << endl;
		cout << "[INFO] Events lost after LepPt > 20,10 GeV cut  = " << _failLepPtCut << endl;
		cout << "[INFO] Events lost after Eta cut  = " << _failEtaCut << endl;
		cout << "[INFO] Events lost after SIP < 4 cut and dxy-dz cuts  = " << _failSipVtxCut << endl;
		cout << "[INFO] Events lost after MET < 25 cut  = " << _failMETCut << endl;
		cout << "[INFO] Events lost after m_ll > 4 cut  = " << _faillingJPsiMassCut << endl;
		cout << "[INFO] Total events left = " << _passingSelection + _faillingSelection << endl;
		cout << "[INFO] Passing selection = " << _passingSelection  << endl;
		cout << "[INFO] Failling selection = " << _faillingSelection << endl;
		cout << "========================================================================================" << endl;
		cout << endl;
	}
   
   cout << "[INFO] Processing of " << input_file_data_name << " done." << endl;
}
//===============================================================================



//===============================================================================
void OSmethod::FillDataMCPlots( TString input_file_data_name )
{
   input_file_data = TFile::Open( input_file_data_name);
   
   hCounters = (TH1F*)input_file_data->Get("CRZLLTree/Counters");
   gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);
   
   input_tree_data = (TTree*)input_file_data->Get("CRZLLTree/candTree");
   Init( input_tree_data, input_file_data_name , true);
   
   _current_process = find_current_process(input_file_data_name);
   
   if (fChain == 0) return;
   
   Long64_t nentries = fChain->GetEntriesFast();
   cout<<"[INFO] Processing "<<input_file_data_name<<", total event: "<<nentries<<endl;

   Long64_t nbytes = 0, nb = 0;
   
   for (Long64_t jentry=0; jentry<nentries;jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry%50000==0) cout<<ientry<<endl;
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      //cout<<"step1"<<endl;
      if (!(test_bit(CRflag, CRZLLos_2P2F)) && !(test_bit(CRflag, CRZLLos_3P1F))) continue;
      //cout<<"step2"<<endl;
      _current_final_state = FindFinalState();
      if (_current_final_state<0) {continue;}
      //cout<<"step3"<<endl;
/*    
      for ( int j = 0; j < nCleanedJetsPt30; j++)
      {
         jetPt[j] = JetPt->at(j);
         jetEta[j] = JetEta->at(j);
         jetPhi[j] = JetPhi->at(j);
         jetMass[j] = JetMass->at(j);
         jetQGL[j] = JetQGLikelihood->at(j);
         jetPgOverPq[j] = 1./JetQGLikelihood->at(j) - 1.;
      }
      
      _current_category = categoryMor18(  nExtraLep,
					  nExtraZ,
					  nCleanedJetsPt30,
					  nCleanedJetsPt30BTagged_bTagSF,
					  jetQGL,
					  p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
					  p_JQCD_SIG_ghg2_1_JHUGen_JECNominal,
					  p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
					  p_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
					  pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
					  p_HadWH_SIG_ghw1_1_JHUGen_JECNominal,
					  p_HadZH_SIG_ghz1_1_JHUGen_JECNominal,
					  p_HadWH_mavjj_JECNominal,
					  p_HadWH_mavjj_true_JECNominal,
					  p_HadZH_mavjj_JECNominal,
					  p_HadZH_mavjj_true_JECNominal,
					  jetPhi,
					  ZZMass,
					  PFMET,
					  false,// Use VHMET category
					  false);// Use QG tagging
      
      _current_category_stxs = stage1_reco_1p1 ( nCleanedJetsPt30,
                                                 DiJetMass,
                                                 ZZPt,
                                                 _current_category,
                                                 ZZjjPt);
*/    //cout<<"step4"<<endl;
      _k_factor = calculate_K_factor(input_file_data_name);
      //cout<<"step5"<<endl;
      _TauIDSF = calculate_TauIDSF(input_file_data_name);
      //cout<<"step6"<<endl;
      _event_weight = (_lumi * 1000 * xsec * _k_factor * _TauIDSF * overallEventWeight * L1prefiringWeight) / gen_sum_weights;
      //cout<<test_bit(CRflag, CRZLLos_2P2F)<<","<<_current_process<<","<<_current_final_state<<endl; 
      if ( test_bit(CRflag, CRZLLos_2P2F) ) histos_1D[Settings::reg2P2F][_current_process][_current_final_state]->Fill(ZZMass, (_current_process == Settings::Data) ? 1 :  _event_weight);
      //cout<<"step7"<<endl;
      if ( test_bit(CRflag, CRZLLos_3P1F) ) histos_1D[Settings::reg3P1F][_current_process][_current_final_state]->Fill(ZZMass, (_current_process == Settings::Data) ? 1 :  _event_weight);
      //cout<<"step8"<<endl;
      if ( Z1Flav < 0 && Z2Flav < 0 )       histos_1D[Settings::regOS][_current_process][_current_final_state]->Fill(ZZMass, (_current_process == Settings::Data) ? 1 :  _event_weight);
      //cout<<"step9"<<endl;   
   } // END events loop
   
   cout << "[INFO] Processing of " << input_file_data_name << " done." << endl;
}
//===============================================================================



//===============================================================================
void OSmethod::MakeHistogramsZX( TString input_file_data_name, TString  input_file_FR_name )
{
   
   FakeRates *FR = new FakeRates( input_file_FR_name );
   
   input_file_data = TFile::Open( input_file_data_name);
   input_tree_data = (TTree*)input_file_data->Get("CRZLLTree/candTree");
   Init( input_tree_data, input_file_data_name , true);
   
   
   if (fChain == 0) return;
   
   Long64_t nentries = fChain->GetEntriesFast();
   cout<<"[INFO] Processing "<<input_file_data_name<<", total event: "<<nentries<<endl;
   
   Long64_t nbytes = 0, nb = 0;
   // FOR DEBUG
   //Int_t nevents_CRLLos      = 0;
   //Int_t nevents_CRLLos_2P2F = 0;
   //Int_t nevents_CRLLos_3P1F = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++)
   {
      
      Long64_t ientry = LoadTree(jentry);
      if (ientry%50000==0) cout<<ientry<<endl;
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      
      if (!(test_bit(CRflag, CRZLLos_2P2F)) && !(test_bit(CRflag, CRZLLos_3P1F))) continue;
      //nevents_CRLLos += 1;	

      // Included SIP and dxy/dz cuts for 3rd and 4th lepton                                                                                                                         
      if ( (fabs(LepEta->at(2)) > 2.5) || (fabs(LepEta->at(3)) > 2.5) ) {continue;}
      if ( (abs(LepLepId->at(2))==15 && fabs(LepEta->at(2)) > 2.3) || (abs(LepLepId->at(3))==15 && fabs(LepEta->at(3)) > 2.3) ) {continue;}
      if ( ( LepSIP->at(2) > 4. || Lepdxy->at(2) > 0.5 || Lepdz->at(2) > 1.0) && (abs(LepLepId->at(2)) == 11 || abs(LepLepId->at(2)) == 13)) {continue;}
      if ( ( LepSIP->at(3) > 4. || Lepdxy->at(3) > 0.5 || Lepdz->at(3) > 1.0) && (abs(LepLepId->at(3)) == 11 || abs(LepLepId->at(3)) == 13)) {continue;}
      if ( ( Lepdz->at(2) > 10 || LepPt->at(2) < 20) && abs(LepLepId->at(2))==15 ) {continue;}
      if ( ( Lepdz->at(3) > 10 || LepPt->at(3) < 20) && abs(LepLepId->at(3))==15 ) {continue;}
      if ( ZZMass < 70. ) continue;
      _current_final_state = FindFinalState();
      if (_current_final_state<0) {continue;}
      
      for ( int j = 0; j < nCleanedJetsPt30; j++)
      {
         jetPt[j] = JetPt->at(j);
         jetEta[j] = JetEta->at(j);
         jetPhi[j] = JetPhi->at(j);
         jetMass[j] = JetMass->at(j);
         jetQGL[j] = JetQGLikelihood->at(j);
         jetPgOverPq[j] = 1./JetQGLikelihood->at(j) - 1.;
      }
/* 
      _current_category = categoryMor18(  nExtraLep,
					  nExtraZ,
					  nCleanedJetsPt30,
					  nCleanedJetsPt30BTagged_bTagSF,
					  jetQGL,
					  p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
					  p_JQCD_SIG_ghg2_1_JHUGen_JECNominal,
					  p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
					  p_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
					  pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
					  p_HadWH_SIG_ghw1_1_JHUGen_JECNominal,
					  p_HadZH_SIG_ghz1_1_JHUGen_JECNominal,
					  p_HadWH_mavjj_JECNominal,
					  p_HadWH_mavjj_true_JECNominal,
					  p_HadZH_mavjj_JECNominal,
					  p_HadZH_mavjj_true_JECNominal,
					  jetPhi,
					  ZZMass,
					  PFMET,
					  false,// Use VHMET category
					  false);// Use QG tagging
      
      _current_category_stxs = stage1_reco_1p1 ( nCleanedJetsPt30,
                                                 DiJetMass,
                                                 ZZPt,
                                                 _current_category,
                                                 ZZjjPt);
*/
      int tauChannel=-1;
      if (Z2Flav==-165) tauChannel=0;
      else if (Z2Flav==-195) tauChannel=1;
      else if (Z2Flav==-225) tauChannel=2;
      //cout<<"step1"<<endl;
      _f3    = FR->GetFakeRate(LepPt->at(2),LepEta->at(2),LepLepId->at(2),tauChannel);
      //cout<<"step2"<<endl;
      _f3_Up = FR->GetFakeRate_Up(LepPt->at(2),LepEta->at(2),LepLepId->at(2),tauChannel);
      _f3_Dn = FR->GetFakeRate_Dn(LepPt->at(2),LepEta->at(2),LepLepId->at(2),tauChannel);
      _f4    = FR->GetFakeRate(LepPt->at(3),LepEta->at(3),LepLepId->at(3),tauChannel);
      //cout<<"step2.05"<<endl;
      _f4_Up = FR->GetFakeRate_Up(LepPt->at(3),LepEta->at(3),LepLepId->at(3),tauChannel);
      //cout<<"step2.1"<<endl;
      _f4_Dn = FR->GetFakeRate_Dn(LepPt->at(3),LepEta->at(3),LepLepId->at(3),tauChannel);
      //cout<<"step2.2, "<<test_bit(CRflag, CRZLLos_2P2F)<<test_bit(CRflag, CRZLLos_3P1F)<<endl;
      if ( test_bit(CRflag, CRZLLos_2P2F) )
      {
	//nevents_CRLLos_2P2F += 1;	
	//cout << "===============" << endl;
	//cout << "f3 = " << _f3 << endl;
	//cout << "f4 = " << _f4 << endl;
	//cout << "weight = " << (_f3/(1-_f3))*(_f4/(1-_f4)) << endl;
	//cout << "weight_up = " << (_f3_Up/(1-_f3_Up))*(_f4_Up/(1-_f4_Up)) << endl;
	//cout << "weight_dn = " << (_f3_Dn/(1-_f3_Dn))*(_f4_Dn/(1-_f4_Dn)) << endl;
	 //cout<<"step3"<<endl;		
         h_from2P2F_SR[Settings::nominal][_current_final_state]->Fill(ZZMass, (_f3/(1-_f3))*(_f4/(1-_f4)) );
	 //cout<<"step4"<<endl;
         h_from2P2F_3P1F[Settings::nominal][_current_final_state]->Fill(ZZMass, (_f3/(1-_f3))+(_f4/(1-_f4)) );
			
         h_from2P2F_SR[Settings::Up][_current_final_state]->Fill(ZZMass, (_f3_Up/(1-_f3_Up))*(_f4_Up/(1-_f4_Up)) );
         h_from2P2F_3P1F[Settings::Up][_current_final_state]->Fill(ZZMass, (_f3_Up/(1-_f3_Up))+(_f4_Up/(1-_f4_Up)) );
			
         h_from2P2F_SR[Settings::Dn][_current_final_state]->Fill(ZZMass, (_f3_Dn/(1-_f3_Dn))*(_f4_Dn/(1-_f4_Dn)) );
         h_from2P2F_3P1F[Settings::Dn][_current_final_state]->Fill(ZZMass, (_f3_Dn/(1-_f3_Dn))+(_f4_Dn/(1-_f4_Dn)) );
      }
      if ( test_bit(CRflag, CRZLLos_3P1F) )
      {
	//nevents_CRLLos_3P1F += 1;
	bool tightLep3=false;
	if (abs(LepLepId->at(3))==11) tightLep3=LepisID->at(3);
	else if (abs(LepLepId->at(3))==13) tightLep3=LepisID->at(3) && LepCombRelIsoPF->at(3) < 0.35;
	else if (abs(LepLepId->at(3))==15) {
	  if (tauChannel==0) tightLep3=(LepisID->at(3) && TauVSmu->at(3)>=4 && TauVSe->at(3)>=5 && TauVSjet->at(3)>=5);
	  else tightLep3=(LepisID->at(3) && TauVSmu->at(3)>=4 && TauVSe->at(3)>=3 && TauVSjet->at(3)>=6);
	}
	//cout<<"tightLep3"<<tightLep3<<endl;
	if(tightLep3)//(LepisID->at(3))
	  {
	    _f4    = FR->GetFakeRate(LepPt->at(2),LepEta->at(2),LepLepId->at(2),tauChannel);
	    _f4_Up = FR->GetFakeRate_Up(LepPt->at(2),LepEta->at(2),LepLepId->at(2),tauChannel);
	    _f4_Dn = FR->GetFakeRate_Dn(LepPt->at(2),LepEta->at(2),LepLepId->at(2),tauChannel);
	  }
	else
	  {
	    _f4    = FR->GetFakeRate(LepPt->at(3),LepEta->at(3),LepLepId->at(3),tauChannel);
	    _f4_Up = FR->GetFakeRate_Up(LepPt->at(3),LepEta->at(3),LepLepId->at(3),tauChannel);
	    _f4_Dn = FR->GetFakeRate_Dn(LepPt->at(3),LepEta->at(3),LepLepId->at(3),tauChannel);
	  }
	
	h_from3P1F_SR[Settings::nominal][_current_final_state]->Fill(ZZMass, (_f4/(1-_f4)) );
	h_from3P1F_SR[Settings::Up][_current_final_state]->Fill(ZZMass, (_f4_Up/(1-_f4_Up)) );
	h_from3P1F_SR[Settings::Dn][_current_final_state]->Fill(ZZMass, (_f4_Dn/(1-_f4_Dn)) );
      }
      
   }
   
   //std::cout << "####################################################\n";
   //std::cout << "# events CRLLos      = " << nevents_CRLLos      << '\n';
   //std::cout << "# events CRLLos_3P1F = " << nevents_CRLLos_3P1F << '\n';
   //std::cout << "# events CRLLos_2P2F = " << nevents_CRLLos_2P2F << '\n';
   //std::cout << "####################################################\n";

   cout << "[INFO] Processing of " << input_file_data_name << " done." << endl;
}
//===============================================================================

//===============================================================================
void OSmethod::MakeZXMCContribution( TString input_file_data_name, TString  input_file_FR_name )
{
   
   FakeRates *FR = new FakeRates( input_file_FR_name );
   input_file_data = TFile::Open( input_file_data_name);
   
   hCounters = (TH1F*)input_file_data->Get("CRZLLTree/Counters");
   gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);
   
   input_tree_data = (TTree*)input_file_data->Get("CRZLLTree/candTree");
   Init( input_tree_data, input_file_data_name , true);
   
   
   if (fChain == 0) return;
   
   Long64_t nentries = fChain->GetEntriesFast();
   cout<<"[INFO] Processing "<<input_file_data_name<<", total event: "<<nentries<<endl;
   
   Long64_t nbytes = 0, nb = 0;
   
   for (Long64_t jentry=0; jentry<nentries;jentry++)
   {
      
      Long64_t ientry = LoadTree(jentry);
      if (ientry%50000==0) cout<<ientry<<endl;
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      
      if (!(test_bit(CRflag, CRZLLos_3P1F))) continue;
      
      _current_final_state = FindFinalState();
      if (_current_final_state<0) {continue;}
      
      for ( int j = 0; j < nCleanedJetsPt30; j++)
      {
         jetPt[j] = JetPt->at(j);
         jetEta[j] = JetEta->at(j);
         jetPhi[j] = JetPhi->at(j);
         jetMass[j] = JetMass->at(j);
         jetQGL[j] = JetQGLikelihood->at(j);
         jetPgOverPq[j] = 1./JetQGLikelihood->at(j) - 1.;
      }
 /*     
      _current_category = categoryMor18(  nExtraLep,
					  nExtraZ,
					  nCleanedJetsPt30,
					  nCleanedJetsPt30BTagged_bTagSF,
					  jetQGL,
					  p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
					  p_JQCD_SIG_ghg2_1_JHUGen_JECNominal,
					  p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
					  p_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
					  pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
					  p_HadWH_SIG_ghw1_1_JHUGen_JECNominal,
					  p_HadZH_SIG_ghz1_1_JHUGen_JECNominal,
					  p_HadWH_mavjj_JECNominal,
					  p_HadWH_mavjj_true_JECNominal,
					  p_HadZH_mavjj_JECNominal,
					  p_HadZH_mavjj_true_JECNominal,
					  jetPhi,
					  ZZMass,
					  PFMET,
					  false,// Use VHMET category
					  false);// Use QG tagging
      
      _current_category_stxs = stage1_reco_1p1 ( nCleanedJetsPt30,
                                                 DiJetMass,
                                                 ZZPt,
                                                 _current_category,
                                                 ZZjjPt);
   */   
      _k_factor = calculate_K_factor(input_file_data_name);
      _TauIDSF = calculate_TauIDSF(input_file_data_name);
      _event_weight = (_lumi * 1000 * xsec * _k_factor * _TauIDSF * overallEventWeight * L1prefiringWeight) / gen_sum_weights;

      int tauChannel=-1;
      if (Z2Flav==-165) tauChannel=0;
      else if (Z2Flav==-195) tauChannel=1;
      else if (Z2Flav==-225) tauChannel=2;

      bool tightLep3;
      if (abs(LepLepId->at(3))==11) tightLep3=LepisID->at(3);
      else if (abs(LepLepId->at(3))==13) tightLep3=LepisID->at(3) && LepCombRelIsoPF->at(3) < 0.35;
      else if (abs(LepLepId->at(3))==15) {
        if (Z2Flav==-165) tightLep3=LepisID->at(3) && TauVSmu->at(3)>=4 && TauVSe->at(3)>=5 && TauVSjet->at(3)>=5;
        else tightLep3=LepisID->at(3) && TauVSmu->at(3)>=4 && TauVSe->at(3)>=3 && TauVSjet->at(3)>=6;
      }
      
      if(tightLep3)//( LepisID->at(3) )
      {
	_f4    = FR->GetFakeRate(LepPt->at(2),LepEta->at(2),LepLepId->at(2),tauChannel);
	_f4_Up = FR->GetFakeRate_Up(LepPt->at(2),LepEta->at(2),LepLepId->at(2),tauChannel);
	_f4_Dn = FR->GetFakeRate_Dn(LepPt->at(2),LepEta->at(2),LepLepId->at(2),tauChannel);
      }
      else
      {
	_f4    = FR->GetFakeRate(LepPt->at(3),LepEta->at(3),LepLepId->at(3),tauChannel);
	_f4_Up = FR->GetFakeRate_Up(LepPt->at(3),LepEta->at(3),LepLepId->at(3),tauChannel);
	_f4_Dn = FR->GetFakeRate_Dn(LepPt->at(3),LepEta->at(3),LepLepId->at(3),tauChannel);
      }
      
      h_from3P1F_SR_ZZonly[Settings::nominal][_current_final_state]->Fill(ZZMass, _event_weight * (_f4/(1-_f4)) );
      h_from3P1F_SR_ZZonly[Settings::Up][_current_final_state]->Fill(ZZMass, _event_weight * (_f4_Up/(1-_f4_Up)) );
      h_from3P1F_SR_ZZonly[Settings::Dn][_current_final_state]->Fill(ZZMass, _event_weight * (_f4_Dn/(1-_f4_Dn)) );
      
   }
   
   cout << "[INFO] Processing of " << input_file_data_name << " done." << endl;
}
//===============================================================================




//===============================================================
void OSmethod::DeclareFRHistos()
{
   for (int i_flav = 0; i_flav < num_of_flavours; i_flav++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
	_histo_name = "Passing_" + _s_process.at(i_proc) + "_" + _s_flavour.at(i_flav);
	passing[i_proc][i_flav] = new TH2F(_histo_name,"", 80, 0, 80, 2, 0, 2);
        
	_histo_name = "Failing_" + _s_process.at(i_proc) + "_" + _s_flavour.at(i_flav);
	failing[i_proc][i_flav] = new TH2F(_histo_name,"", 80, 0, 80, 2, 0, 2);
	
      }
      
      _histo_name = "Passing_Total_" + _s_flavour.at(i_flav);
      passing[Settings::Total][i_flav] = new TH2F(_histo_name,"", 80, 0, 80, 2, 0, 2);
      _histo_name = "Failing_Total_" + _s_flavour.at(i_flav);
      failing[Settings::Total][i_flav] = new TH2F(_histo_name,"", 80, 0, 80, 2, 0, 2);

   }

}
//===============================================================

//===============================================================
void OSmethod::DeclareDataMCHistos()
{
   for (int i_reg = 0; i_reg < num_of_regions_os; i_reg ++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         for (int i_fs = 0; i_fs < num_of_final_states; i_fs++)
         {
            //for (int i_cat = 0; i_cat < num_of_categories_stxs; i_cat++)
            //{
               _histo_name = "M4l_" + _s_region.at(i_reg) + "_" + _s_process.at(i_proc) + "_" + _s_final_state.at(i_fs);// + "_" + _s_category_stxs.at(i_cat);
               _histo_labels = ";" + Plots::M4l().var_X_label + ";" + Plots::M4l().var_Y_label;
               histos_1D[i_reg][i_proc][i_fs] = new TH1F(_histo_name, _histo_labels, Plots::M4l().var_N_bin, Plots::M4l().var_min, Plots::M4l().var_max);
            //}
         }
      }
   }
   
}
//===============================================================

//===============================================================
void OSmethod::DeclareZXHistos()
{
   for (int i_fs = 0; i_fs < num_of_final_states; i_fs++)
   {
      //for (int i_cat = 0; i_cat < num_of_categories_stxs; i_cat++)
      //{
      	for(int i_var = 0; i_var < num_of_fr_variations; i_var++)
      	{
				_histo_name = "h_from2P2F_SR_" + _s_variation.at(i_var) + "_" + _s_final_state.at(i_fs);// + "_" + _s_category_stxs.at(i_cat);
				_histo_labels = ";" + Plots::M4l().var_X_label + ";" + Plots::M4l().var_Y_label;
				h_from2P2F_SR[i_var][i_fs] = new TH1F(_histo_name, _histo_labels, Plots::M4l().var_N_bin, Plots::M4l().var_min, Plots::M4l().var_max);
				
				_histo_name = "h_from2P2F_3P1F_" + _s_variation.at(i_var) + "_" + _s_final_state.at(i_fs);// + "_" + _s_category_stxs.at(i_cat);
				_histo_labels = ";" + Plots::M4l().var_X_label + ";" + Plots::M4l().var_Y_label;
				h_from2P2F_3P1F[i_var][i_fs] = new TH1F(_histo_name, _histo_labels, Plots::M4l().var_N_bin, Plots::M4l().var_min, Plots::M4l().var_max);
				
				_histo_name = "h_from3P1F_SR_final_" + _s_variation.at(i_var) + "_" + _s_final_state.at(i_fs);// + "_" + _s_category_stxs.at(i_cat);
				_histo_labels = ";" + Plots::M4l().var_X_label + ";" + Plots::M4l().var_Y_label;
				h_from3P1F_SR_final[i_var][i_fs] = new TH1F(_histo_name, _histo_labels, Plots::M4l().var_N_bin, Plots::M4l().var_min, Plots::M4l().var_max);
				
				_histo_name = "h_from3P1F_SR_" + _s_variation.at(i_var) + "_" + _s_final_state.at(i_fs);// + "_" + _s_category_stxs.at(i_cat);
				_histo_labels = ";" + Plots::M4l().var_X_label + ";" + Plots::M4l().var_Y_label;
				h_from3P1F_SR[i_var][i_fs] = new TH1F(_histo_name, _histo_labels, Plots::M4l().var_N_bin, Plots::M4l().var_min, Plots::M4l().var_max);
				
				_histo_name = "h_from3P1F_SR_ZZonly_" + _s_variation.at(i_var) + "_" + _s_final_state.at(i_fs);// + "_" + _s_category_stxs.at(i_cat);
				_histo_labels = ";" + Plots::M4l().var_X_label + ";" + Plots::M4l().var_Y_label;
				h_from3P1F_SR_ZZonly[i_var][i_fs] = new TH1F(_histo_name, _histo_labels, Plots::M4l().var_N_bin, Plots::M4l().var_min, Plots::M4l().var_max);
				
				_histo_name = "ZX_" + _s_variation.at(i_var) + "_" + _s_final_state.at(i_fs);// + "_" + _s_category_stxs.at(i_cat);
				_histo_labels = ";" + Plots::M4l().var_X_label + ";" + Plots::M4l().var_Y_label;
				histos_ZX[i_var][i_fs] = new TH1F(_histo_name, _histo_labels, Plots::M4l().var_N_bin, Plots::M4l().var_min, Plots::M4l().var_max);
				
			}
			
      //}
   }
}
//===============================================================

//===============================================================
void OSmethod::SaveFRHistos( TString file_name,  bool subtractWZ, bool remove_negative_bins)
{
   TFile* fOutHistos = TFile::Open(file_name, "recreate");
   fOutHistos->cd();
   
   // Copy data histos to total histos, if there is no WZ subtraction this is the final histo for fake rate calculation
   for (int i_flav = 0; i_flav < num_of_flavours; i_flav++)
   {
      passing[Settings::Total][i_flav]->Add(passing[Settings::Data][i_flav], 1.);
      failing[Settings::Total][i_flav]->Add(failing[Settings::Data][i_flav], 1.);
   }
   
   if (subtractWZ ) SubtractWZ(); // Subtract WZ contribution from MC estimate
   
   if ( remove_negative_bins ) // Set negative bins to zero
   {
      for (int i_flav = 0; i_flav < num_of_flavours; i_flav++)
      {
         RemoveNegativeBins2D( passing[Settings::Total][i_flav] );
         RemoveNegativeBins2D( failing[Settings::Total][i_flav] );
      }
      cout << "[INFO] Negative bins removed." << endl;
   }

   for (int i_flav = 0; i_flav < num_of_flavours; i_flav++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         passing[i_proc][i_flav]->Write();
         failing[i_proc][i_flav]->Write();
      }
      
   passing[Settings::Total][i_flav]->Write();
   failing[Settings::Total][i_flav]->Write();
      
   }
   
   fOutHistos->Close();
   delete fOutHistos;
   
   cout << "[INFO] All FakeRate histograms saved." << endl;
}
//===============================================================

//===============================================================
void OSmethod::SaveDataMCHistos( TString file_name )
{
   FillDataMCInclusive();
   
   TFile* fOutHistos = TFile::Open(file_name, "recreate");
   fOutHistos->cd();
   
   for (int i_reg = 0; i_reg < num_of_regions_os; i_reg ++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         for (int i_fs = 0; i_fs < num_of_final_states; i_fs++)
         {
            //for (int i_cat = 0; i_cat < num_of_categories_stxs; i_cat++)
            //{
               histos_1D[i_reg][i_proc][i_fs]->Write();
            //}
         }
      }
   }
   
   fOutHistos->Close();
   delete fOutHistos;
   
   cout << "[INFO] All Data/MC histograms saved." << endl;
}
//===============================================================

//===============================================================
void OSmethod::FillDataMCInclusive( )
{
/*   for (int i_reg = 0; i_reg < num_of_regions_os; i_reg ++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         for (int i_fs = 0; i_fs < Settings::fs4l; i_fs++)
         {
            //for (int i_cat = 0; i_cat < Settings::inclusive_stxs; i_cat++)
            //{
               histos_1D[i_reg][i_proc][i_fs]->Add(histos_1D[i_reg][i_proc][i_fs][i_cat]);
               histos_1D[i_reg][i_proc][Settings::fs4l][i_cat]    ->Add(histos_1D[i_reg][i_proc][i_fs][i_cat]);
            //}
         }
      }
   }
*/   
   for (int i_reg = 0; i_reg < num_of_regions_os; i_reg ++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         for (int i_fs = 0; i_fs < Settings::fs4l; i_fs++)
         {
            histos_1D[i_reg][i_proc][Settings::fs4l]->Add(histos_1D[i_reg][i_proc][i_fs]);
         }
      }
   }
   
   cout << "[INFO] All Data/MC histograms summed." << endl;
}
//===============================================================

//===============================================================
void OSmethod::SaveZXHistos( TString file_name , bool remove_negative_bins)
{
   FillZXInclusive(remove_negative_bins);
   
   TFile* fOutHistos = TFile::Open(file_name, "recreate");
   fOutHistos->cd();
   
   for (int i_fs = 0; i_fs < num_of_final_states; i_fs++)
   {
      //for (int i_cat = 0; i_cat < num_of_categories_stxs; i_cat++)
      //{
			for(int i_var = 0; i_var < num_of_fr_variations; i_var++)
      	{
	        h_from2P2F_SR[i_var][i_fs]->Write();
         	h_from2P2F_3P1F[i_var][i_fs]->Write();
         	h_from3P1F_SR_final[i_var][i_fs]->Write();
         	h_from3P1F_SR[i_var][i_fs]->Write();
		h_from3P1F_SR_ZZonly[i_var][i_fs]->Write();
         	histos_ZX[i_var][i_fs]->Write();
      	}

      //}
   }
   
   fOutHistos->Close();
   delete fOutHistos;
   
   cout << "[INFO] All Z+X histograms saved." << endl;
}
//===============================================================

//===============================================================
void OSmethod::FillZXInclusive( bool remove_negative_bins )
{
  if ( remove_negative_bins )
    {
      for (int i_fs = 0; i_fs <= Settings::fs4l; i_fs++)
	{
	  //for (int i_cat = 0; i_cat <= Settings::inclusive_stxs; i_cat++)
	    //{
	      for(int i_var = 0; i_var < num_of_fr_variations; i_var++)
      		{
		  RemoveNegativeBins1D(h_from2P2F_SR[i_var][i_fs]);
		  RemoveNegativeBins1D(h_from2P2F_3P1F[i_var][i_fs]);
		  RemoveNegativeBins1D(h_from3P1F_SR_final[i_var][i_fs]);
		  RemoveNegativeBins1D(h_from3P1F_SR[i_var][i_fs]);
		  RemoveNegativeBins1D(h_from3P1F_SR_ZZonly[i_var][i_fs]);
		  
		}
	    //}
	}
    }
	
  for (int i_fs = 0; i_fs < Settings::fs4l; i_fs++)
    {
      //for (int i_cat = 0; i_cat < Settings::inclusive_stxs; i_cat++)
	//{
	  for(int i_var = 0; i_var < num_of_fr_variations; i_var++)
	    {
	      //h_from2P2F_SR[i_var][i_fs][Settings::inclusive_stxs]->Add(h_from2P2F_SR[i_var][i_fs]);
	      h_from2P2F_SR[i_var][Settings::fs4l]    ->Add(h_from2P2F_SR[i_var][i_fs]);
	      
	      //h_from2P2F_3P1F[i_var][i_fs][Settings::inclusive_stxs]->Add(h_from2P2F_3P1F[i_var][i_fs][i_cat]);
	      h_from2P2F_3P1F[i_var][Settings::fs4l]    ->Add(h_from2P2F_3P1F[i_var][i_fs]);
	      
	      //h_from3P1F_SR_final[i_var][i_fs][Settings::inclusive_stxs]->Add(h_from3P1F_SR_final[i_var][i_fs][i_cat]);
	      h_from3P1F_SR_final[i_var][Settings::fs4l]    ->Add(h_from3P1F_SR_final[i_var][i_fs]);
	      
	      //h_from3P1F_SR[i_var][i_fs][Settings::inclusive_stxs]->Add(h_from3P1F_SR[i_var][i_fs][i_cat]);
	      h_from3P1F_SR[i_var][Settings::fs4l]    ->Add(h_from3P1F_SR[i_var][i_fs]);
	      
	      //h_from3P1F_SR_ZZonly[i_var][i_fs][Settings::inclusive_stxs]->Add(h_from3P1F_SR_ZZonly[i_var][i_fs][i_cat]);
	      h_from3P1F_SR_ZZonly[i_var][Settings::fs4l]    ->Add(h_from3P1F_SR_ZZonly[i_var][i_fs]);
	    }
	//}
    }
  
  for (int i_fs = 0; i_fs <= Settings::fs4l; i_fs++)
    {
      //for (int i_cat = 0; i_cat <= Settings::inclusive_stxs; i_cat++)
	//{
	  for(int i_var = 0; i_var < num_of_fr_variations; i_var++)
	    {
	      h_from3P1F_SR_final[i_var][i_fs]->Add(h_from3P1F_SR[i_var][i_fs], 1.);
	      h_from3P1F_SR_final[i_var][i_fs]->Add(h_from3P1F_SR_ZZonly[i_var][i_fs], -1.);
	      h_from3P1F_SR_final[i_var][i_fs]->Add(h_from2P2F_SR[i_var][i_fs], -2.);
	    }
	//}
    }
  
  if ( remove_negative_bins )
    {
      for (int i_fs = 0; i_fs <= Settings::fs4l; i_fs++)
	{
	  //for (int i_cat = 0; i_cat <= Settings::inclusive_stxs; i_cat++)
	    //{
	      for(int i_var = 0; i_var < num_of_fr_variations; i_var++)
      		{
		  RemoveNegativeBins1D(h_from3P1F_SR_final[i_var][i_fs]);
		}
	    //}
	}
    }
  
  for (int i_fs = 0; i_fs < Settings::fs4l; i_fs++)
    {
      for(int i_var = 0; i_var < num_of_fr_variations; i_var++)
	{
	  h_from2P2F_SR[i_var][Settings::fs4l]->Add(h_from2P2F_SR[i_var][i_fs]);
	  h_from2P2F_3P1F[i_var][Settings::fs4l]->Add(h_from2P2F_3P1F[i_var][i_fs]);
	  h_from3P1F_SR_final[i_var][Settings::fs4l]->Add(h_from3P1F_SR_final[i_var][i_fs]);
	  h_from3P1F_SR[i_var][Settings::fs4l]->Add(h_from3P1F_SR[i_var][i_fs]);
	  h_from3P1F_SR_ZZonly[i_var][Settings::fs4l]->Add(h_from3P1F_SR_ZZonly[i_var][i_fs]);
	}
      
    }
  
  for (int i_fs = 0; i_fs <= Settings::fs4l; i_fs++)
    {
      //for (int i_cat = 0; i_cat <= Settings::inclusive_stxs; i_cat++)
	//{
	  for(int i_var = 0; i_var < num_of_fr_variations; i_var++)
	    {
	      histos_ZX[i_var][i_fs]->Add(h_from3P1F_SR_final[i_var][i_fs], 1.);
	      histos_ZX[i_var][i_fs]->Add(h_from2P2F_SR[i_var][i_fs], 1.);
	    }
	//}
    }
  
  
  cout << "[INFO] All Z+X histograms summed." << endl;
}
//===============================================================

//===============================================================
void OSmethod::GetFRHistos( TString file_name)
{
   TFile* histo_file = TFile::Open(file_name);
   
   for (int i_flav = 0; i_flav < num_of_flavours; i_flav++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         _histo_name = "Passing_" + _s_process.at(i_proc) + "_" + _s_flavour.at(i_flav);
         passing[i_proc][i_flav] = (TH2F*)histo_file->Get(_histo_name);
         
         _histo_name = "Failing_" + _s_process.at(i_proc) + "_" + _s_flavour.at(i_flav);
         failing[i_proc][i_flav] = (TH2F*)histo_file->Get(_histo_name);
         
      }
      
      _histo_name = "Passing_Total_" + _s_flavour.at(i_flav);
      passing[Settings::Total][i_flav] = (TH2F*)histo_file->Get(_histo_name);
      _histo_name = "Failing_Total_" + _s_flavour.at(i_flav);
      failing[Settings::Total][i_flav] = (TH2F*)histo_file->Get(_histo_name);
      
   }
   
   cout << "[INFO] All FakeRate histograms retrieved from file." << endl;
}
//===============================================================

//===============================================================
void OSmethod::GetDataMCHistos( TString file_name)
{
   TFile* histo_file = TFile::Open(file_name);
   
   for (int i_reg = 0; i_reg < num_of_regions_os; i_reg ++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         for (int i_fs = 0; i_fs < num_of_final_states; i_fs++)
         {
            //for (int i_cat = 0; i_cat < num_of_categories_stxs; i_cat++)
            //{
               _histo_name = "M4l_" + _s_region.at(i_reg) + "_" + _s_process.at(i_proc) + "_" + _s_final_state.at(i_fs);// + "_" + _s_category_stxs.at(i_cat);
               histos_1D[i_reg][i_proc][i_fs] = (TH1F*)histo_file->Get(_histo_name);
            //}
         }
      }
   }
   
   cout << "[INFO] All Data/MC histograms retrieved from file." << endl;
}

//===============================================================

//===============================================================
void OSmethod::GetZXHistos( TString file_name)
{
  TFile* histo_file = TFile::Open(file_name);
  
  for (int i_fs = 0; i_fs < num_of_final_states; i_fs++)
    {
      //for (int i_cat = 0; i_cat < num_of_categories_stxs; i_cat++)
	//{
	  for(int i_var = 0; i_var < num_of_fr_variations; i_var++)
	    {
	      _histo_name = "h_from2P2F_SR_" + _s_variation.at(i_var )+ "_" + _s_final_state.at(i_fs);// + "_" + _s_category_stxs.at(i_cat);
	      h_from2P2F_SR[i_var][i_fs] = (TH1F*)histo_file->Get(_histo_name);
	      
	      _histo_name = "h_from2P2F_3P1F_" + _s_variation.at(i_var )+ "_" + _s_final_state.at(i_fs);// + "_" + _s_category_stxs.at(i_cat);
	      h_from2P2F_3P1F[i_var][i_fs] = (TH1F*)histo_file->Get(_histo_name);
	      
	      _histo_name = "h_from3P1F_SR_final_" + _s_variation.at(i_var )+ "_" + _s_final_state.at(i_fs);// + "_" + _s_category_stxs.at(i_cat);
	      h_from3P1F_SR_final[i_var][i_fs] = (TH1F*)histo_file->Get(_histo_name);
	      
	      _histo_name = "h_from3P1F_SR_" + _s_variation.at(i_var )+ "_" + _s_final_state.at(i_fs);// + "_" + _s_category_stxs.at(i_cat);
	      h_from3P1F_SR[i_var][i_fs] = (TH1F*)histo_file->Get(_histo_name);
	      
	      _histo_name = "h_from3P1F_SR_ZZonly_" + _s_variation.at(i_var )+ "_" + _s_final_state.at(i_fs);// + "_" + _s_category_stxs.at(i_cat);
	      h_from3P1F_SR_ZZonly[i_var][i_fs] = (TH1F*)histo_file->Get(_histo_name);
	      
	      _histo_name = "ZX_" + _s_variation.at(i_var )+ "_" + _s_final_state.at(i_fs);// + "_" + _s_category_stxs.at(i_cat);
	      histos_ZX[i_var][i_fs] = (TH1F*)histo_file->Get(_histo_name);
	    }
	//}
    }
  
  cout << "[INFO] All Z+X histograms retrieved from file." << endl;
}

//===============================================================



void OSmethod::PrintZXYields()
{
  double stat;
  double syst;
  double comb;
  double syst_comp;
  double yield, yield_up;
	
  cout << endl;
  cout << "==============================================================" << endl;
  cout << "[INFO] Control printout for OS Z+X yields in final states "<< endl;
  cout << "==============================================================" << endl;
  for( int i_fs = 0; i_fs < Settings::fs4l ; i_fs++ )
    {
      //for ( int i_cat = 0; i_cat <= Settings::inclusive_stxs; i_cat++)
	//{
	  
	  yield = histos_ZX[Settings::nominal][i_fs]->IntegralAndError(0,histos_ZX[Settings::nominal][i_fs]->GetSize() - 2,stat); //statistical uncertainty
	  yield_up = histos_ZX[Settings::Up][i_fs]->Integral();
	  syst = ((yield_up/yield) - 1.) * yield; //systematical uncertainty due to fake rate variation
	  syst_comp = yield*0.3; //background composition uncertainty of 30% measured in Run I
	  comb = sqrt(stat*stat + syst*syst + syst_comp*syst_comp);
	  
	  cout << "Final state: " << _s_final_state.at(i_fs) << endl;
	  cout << yield << " +/- " << comb << " (total.)   - " << stat << " (stat.)   - " << syst << " (syst.)" << endl;
	//}
      
      cout << "============================================================" << endl;
    }
}



//========================================================================================================
void OSmethod::PlotDataMC_2P2F( TString variable_name, TString folder )
{
   TCanvas *c;
   c = new TCanvas("2P2F", variable_name, 600, 600);
   
   if ( GetVarLogX( variable_name) ) c->SetLogx();
   if ( GetVarLogY( variable_name) ) c->SetLogy();
   
   for( int i_fs = 0; i_fs < Settings::fs4l ; i_fs++ )
   {
      histos_1D[Settings::reg2P2F][Settings::WZ][i_fs]   ->SetFillColor(kMagenta-7);
      histos_1D[Settings::reg2P2F][Settings::qqZZ][i_fs] ->SetFillColor(kCyan+1);
      histos_1D[Settings::reg2P2F][Settings::DY][i_fs]   ->SetFillColor(kGreen+2);
      histos_1D[Settings::reg2P2F][Settings::ttbar][i_fs]->SetFillColor(kBlue-4);
      
      histos_1D[Settings::reg2P2F][Settings::WZ][i_fs]   ->SetLineColor(kMagenta-7);
      histos_1D[Settings::reg2P2F][Settings::qqZZ][i_fs] ->SetLineColor(kCyan+1);
      histos_1D[Settings::reg2P2F][Settings::DY][i_fs]   ->SetLineColor(kGreen+2);
      histos_1D[Settings::reg2P2F][Settings::ttbar][i_fs]->SetLineColor(kBlue-4);
      
      histos_1D[Settings::reg2P2F][Settings::Data][i_fs]->SetMarkerSize(0.8);
      histos_1D[Settings::reg2P2F][Settings::Data][i_fs]->SetMarkerStyle(20);
      histos_1D[Settings::reg2P2F][Settings::Data][i_fs]->SetBinErrorOption(TH1::kPoisson);
      histos_1D[Settings::reg2P2F][Settings::Data][i_fs]->SetLineColor(kBlack);
      
      THStack *stack = new THStack( "stack", "stack" );
		stack->Add(histos_1D[Settings::reg2P2F][Settings::qqZZ][i_fs]);
      stack->Add(histos_1D[Settings::reg2P2F][Settings::WZ][i_fs]);
		stack->Add(histos_1D[Settings::reg2P2F][Settings::ttbar][i_fs]);
      stack->Add(histos_1D[Settings::reg2P2F][Settings::DY][i_fs]);
   
      stack->Draw("HIST");
      
      float data_max = histos_1D[Settings::reg2P2F][Settings::Data][i_fs]->GetBinContent(histos_1D[Settings::reg2P2F][Settings::Data][i_fs]->GetMaximumBin());
      float data_max_error = histos_1D[Settings::reg2P2F][Settings::Data][i_fs]->GetBinErrorUp(histos_1D[Settings::reg2P2F][Settings::Data][i_fs]->GetMaximumBin());
      
      stack->SetMinimum(1e-5);
      stack->SetMaximum((data_max + data_max_error)*1.1);
      
      TString _fs_label;
      if ( i_fs == Settings::fs4mu) _fs_label = "m_{4#font[12]{#mu}} (GeV)";
      if ( i_fs == Settings::fs2muetau) _fs_label = "m_{2#font[12]{#mu}#font[12]{#tau_{e}#tau_{h}}} (GeV)";
      if ( i_fs == Settings::fs2mumutau) _fs_label = "m_{2#font[12]{#mu}#font[12]{#tau_{#mu}#tau_{h}}} (GeV)";
      if ( i_fs == Settings::fs2mutautau) _fs_label = "m_{2#font[12]{#mu}#font[12]{#tau_{h}#tau_{h}}} (GeV)";
      if ( i_fs == Settings::fs2e2mu) _fs_label = "m_{2#font[12]{e}2#font[12]{#mu}} (GeV)";
      if ( i_fs == Settings::fs2eetau) _fs_label = "m_{2#font[12]{e}#font[12]{#tau_{e}#tau_{h}}} (GeV)";
      if ( i_fs == Settings::fs2emutau) _fs_label = "m_{2#font[12]{e}#font[12]{#tau_{#mu}#tau_{h}}} (GeV)";
      if ( i_fs == Settings::fs2etautau) _fs_label = "m_{2#font[12]{e}#font[12]{#tau_{h}#tau_{h}}} (GeV)";
      stack->GetXaxis()->SetTitle(_fs_label);
      stack->GetXaxis()->SetTitleSize(0.04);
      stack->GetXaxis()->SetLabelSize(0.04);
      stack->GetYaxis()->SetTitle(histos_1D[Settings::reg2P2F][Settings::Data][i_fs]->GetYaxis()->GetTitle());
      stack->GetYaxis()->SetTitleSize(0.04);
      stack->GetYaxis()->SetLabelSize(0.04);
      
      stack->GetXaxis()->SetTitleOffset(1.2);
      stack->GetYaxis()->SetTitleOffset(1.25);
      
      histos_1D[Settings::reg2P2F][Settings::Data][i_fs]->Draw("SAME p E1 X0");
      
      TLegend *legend;
      legend  = CreateLegend_2P2F("right",histos_1D[Settings::reg2P2F][Settings::Data][i_fs],histos_1D[Settings::reg2P2F][Settings::WZ][i_fs],histos_1D[Settings::reg2P2F][Settings::qqZZ][i_fs],histos_1D[Settings::reg2P2F][Settings::DY][i_fs],histos_1D[Settings::reg2P2F][Settings::ttbar][i_fs]);
      legend->Draw();

      // Draw lumi
      CMS_lumi *lumi = new CMS_lumi;
      lumi->set_lumi(c, _lumi, 0);
      
      TString _out_file_name;
      _out_file_name = folder + "/" + variable_name + "_OS_2P2F_" + _s_final_state.at(i_fs);// + "_" + _s_category_stxs.at(Settings::inclusive_stxs);
      SavePlots(c, _out_file_name);

   }
}
//========================================================================================================


//========================================================================================================
void OSmethod::PlotDataMC_3P1F( TString variable_name, TString folder )
{
   TCanvas *c;
   c = new TCanvas("3P2F", variable_name, 600, 600);
   
   if ( GetVarLogX( variable_name) ) c->SetLogx();
   if ( GetVarLogY( variable_name) ) c->SetLogy();
   
   for( int i_fs = 0; i_fs < Settings::fs4l ; i_fs++ )
   {
      histos_1D[Settings::reg3P1F][Settings::WZ][i_fs]   ->SetFillColor(kMagenta-7);
      histos_1D[Settings::reg3P1F][Settings::qqZZ][i_fs] ->SetFillColor(kCyan+1);
      histos_1D[Settings::reg3P1F][Settings::DY][i_fs]   ->SetFillColor(kGreen+2);
      histos_1D[Settings::reg3P1F][Settings::ttbar][i_fs]->SetFillColor(kBlue-4);
      
      histos_1D[Settings::reg3P1F][Settings::WZ][i_fs]   ->SetLineColor(kMagenta-7);
      histos_1D[Settings::reg3P1F][Settings::qqZZ][i_fs] ->SetLineColor(kCyan+1);
      histos_1D[Settings::reg3P1F][Settings::DY][i_fs]   ->SetLineColor(kGreen+2);
      histos_1D[Settings::reg3P1F][Settings::ttbar][i_fs]->SetLineColor(kBlue-4);
      
      h_from2P2F_3P1F[Settings::nominal][i_fs]->SetLineColor(kRed);
      
      histos_1D[Settings::reg3P1F][Settings::Data][i_fs]->SetMarkerSize(0.8);
      histos_1D[Settings::reg3P1F][Settings::Data][i_fs]->SetMarkerStyle(20);
      histos_1D[Settings::reg3P1F][Settings::Data][i_fs]->SetBinErrorOption(TH1::kPoisson);
      histos_1D[Settings::reg3P1F][Settings::Data][i_fs]->SetLineColor(kBlack);
      
      THStack *stack = new THStack( "stack", "stack" );
      stack->Add(histos_1D[Settings::reg3P1F][Settings::qqZZ][i_fs]);
      stack->Add(histos_1D[Settings::reg3P1F][Settings::WZ][i_fs]);
      stack->Add(histos_1D[Settings::reg3P1F][Settings::ttbar][i_fs]);
      stack->Add(histos_1D[Settings::reg3P1F][Settings::DY][i_fs]);
		
      stack->Draw("HIST");
      
      h_from2P2F_3P1F[Settings::nominal][i_fs]->Draw("HIST SAME");
      
      float data_max = histos_1D[Settings::reg3P1F][Settings::Data][i_fs]->GetBinContent(histos_1D[Settings::reg3P1F][Settings::Data][i_fs]->GetMaximumBin());
      float data_max_error = histos_1D[Settings::reg3P1F][Settings::Data][i_fs]->GetBinErrorUp(histos_1D[Settings::reg3P1F][Settings::Data][i_fs]->GetMaximumBin());
      
      stack->SetMinimum(1e-5);
      stack->SetMaximum((data_max + data_max_error)*1.35);
      
      TString _fs_label;
      if ( i_fs == Settings::fs4mu) _fs_label = "m_{4#font[12]{#mu}} (GeV)";
      if ( i_fs == Settings::fs2muetau) _fs_label = "m_{2#font[12]{#mu}#font[12]{#tau_{e}#tau_{h}}} (GeV)";
      if ( i_fs == Settings::fs2mumutau) _fs_label = "m_{2#font[12]{#mu}#font[12]{#tau_{#mu}#tau_{h}}} (GeV)";
      if ( i_fs == Settings::fs2mutautau) _fs_label = "m_{2#font[12]{#mu}#font[12]{#tau_{h}#tau_{h}}} (GeV)";
      if ( i_fs == Settings::fs2e2mu) _fs_label = "m_{2#font[12]{e}2#font[12]{#mu}} (GeV)";
      if ( i_fs == Settings::fs2eetau) _fs_label = "m_{2#font[12]{e}#font[12]{#tau_{e}#tau_{h}}} (GeV)";
      if ( i_fs == Settings::fs2emutau) _fs_label = "m_{2#font[12]{e}#font[12]{#tau_{#mu}#tau_{h}}} (GeV)";
      if ( i_fs == Settings::fs2etautau) _fs_label = "m_{2#font[12]{e}#font[12]{#tau_{h}#tau_{h}}} (GeV)";
      stack->GetXaxis()->SetTitle(_fs_label);
      stack->GetXaxis()->SetTitleSize(0.04);
      stack->GetXaxis()->SetLabelSize(0.04);
      stack->GetYaxis()->SetTitle(histos_1D[Settings::reg3P1F][Settings::Data][i_fs]->GetYaxis()->GetTitle());
      stack->GetYaxis()->SetTitleSize(0.04);
      stack->GetYaxis()->SetLabelSize(0.04);
      
      stack->GetXaxis()->SetTitleOffset(1.2);
      stack->GetYaxis()->SetTitleOffset(1.25);
      
      histos_1D[Settings::reg3P1F][Settings::Data][i_fs]->Draw("SAME p E1 X0");
      
      TLegend *legend;
      legend  = CreateLegend_3P1F("right",histos_1D[Settings::reg3P1F][Settings::Data][i_fs],h_from2P2F_3P1F[Settings::nominal][i_fs],histos_1D[Settings::reg3P1F][Settings::WZ][i_fs],histos_1D[Settings::reg3P1F][Settings::qqZZ][i_fs],histos_1D[Settings::reg3P1F][Settings::DY][i_fs],histos_1D[Settings::reg3P1F][Settings::ttbar][i_fs]);
      legend->Draw();
      
      // Draw lumi
      CMS_lumi *lumi = new CMS_lumi;
      lumi->set_lumi(c, _lumi, 0);
      
      TString _out_file_name;
      _out_file_name = folder + "/" + variable_name + "_OS_3P1F_" + _s_final_state.at(i_fs);// + "_" + _s_category_stxs.at(Settings::inclusive_stxs);
      SavePlots(c, _out_file_name);
      
   }
}
//========================================================================================================

//========================================================================================================
void OSmethod::PlotDataMC( TString variable_name, TString folder )
{
   TCanvas *c;
   c = new TCanvas("OS", variable_name, 600, 600);
	
   if ( GetVarLogX( variable_name) ) c->SetLogx();
   if ( GetVarLogY( variable_name) ) c->SetLogy();
	
   for( int i_fs = 0; i_fs < Settings::fs4l ; i_fs++ )
   {
      histos_1D[Settings::regOS][Settings::WZ][i_fs]   ->SetFillColor(kMagenta-7);
      histos_1D[Settings::regOS][Settings::qqZZ][i_fs] ->SetFillColor(kCyan+1);
      histos_1D[Settings::regOS][Settings::DY][i_fs]   ->SetFillColor(kGreen+2);
      histos_1D[Settings::regOS][Settings::ttbar][i_fs]->SetFillColor(kBlue-4);
		
      histos_1D[Settings::regOS][Settings::WZ][i_fs]   ->SetLineColor(kMagenta-7);
      histos_1D[Settings::regOS][Settings::qqZZ][i_fs] ->SetLineColor(kCyan+1);
      histos_1D[Settings::regOS][Settings::DY][i_fs]   ->SetLineColor(kGreen+2);
      histos_1D[Settings::regOS][Settings::ttbar][i_fs]->SetLineColor(kBlue-4);
		
      histos_1D[Settings::regOS][Settings::Data][i_fs]->SetMarkerSize(0.8);
      histos_1D[Settings::regOS][Settings::Data][i_fs]->SetMarkerStyle(20);
      histos_1D[Settings::regOS][Settings::Data][i_fs]->SetBinErrorOption(TH1::kPoisson);
      histos_1D[Settings::regOS][Settings::Data][i_fs]->SetLineColor(kBlack);
		
      THStack *stack = new THStack( "stack", "stack" );
      stack->Add(histos_1D[Settings::regOS][Settings::qqZZ][i_fs]);
      stack->Add(histos_1D[Settings::regOS][Settings::WZ][i_fs]);
      stack->Add(histos_1D[Settings::regOS][Settings::ttbar][i_fs]);
      stack->Add(histos_1D[Settings::regOS][Settings::DY][i_fs]);
		
      stack->Draw("HIST");
		
      float data_max = histos_1D[Settings::regOS][Settings::Data][i_fs]->GetBinContent(histos_1D[Settings::regOS][Settings::Data][i_fs]->GetMaximumBin());
      float data_max_error = histos_1D[Settings::regOS][Settings::Data][i_fs]->GetBinErrorUp(histos_1D[Settings::regOS][Settings::Data][i_fs]->GetMaximumBin());
		
      stack->SetMinimum(1e-5);
      stack->SetMaximum((data_max + data_max_error)*1.35);
		
      TString _fs_label;
      if ( i_fs == Settings::fs4mu) _fs_label = "m_{4#font[12]{#mu}} (GeV)";
      if ( i_fs == Settings::fs2muetau) _fs_label = "m_{2#font[12]{#mu}#font[12]{#tau_{e}#tau_{h}}} (GeV)";
      if ( i_fs == Settings::fs2mumutau) _fs_label = "m_{2#font[12]{#mu}#font[12]{#tau_{#mu}#tau_{h}}} (GeV)";
      if ( i_fs == Settings::fs2mutautau) _fs_label = "m_{2#font[12]{#mu}#font[12]{#tau_{h}#tau_{h}}} (GeV)";
      if ( i_fs == Settings::fs2e2mu) _fs_label = "m_{2#font[12]{e}2#font[12]{#mu}} (GeV)";
      if ( i_fs == Settings::fs2eetau) _fs_label = "m_{2#font[12]{e}#font[12]{#tau_{e}#tau_{h}}} (GeV)";
      if ( i_fs == Settings::fs2emutau) _fs_label = "m_{2#font[12]{e}#font[12]{#tau_{#mu}#tau_{h}}} (GeV)";
      if ( i_fs == Settings::fs2etautau) _fs_label = "m_{2#font[12]{e}#font[12]{#tau_{h}#tau_{h}}} (GeV)";
      stack->GetXaxis()->SetTitle(_fs_label);
      stack->GetXaxis()->SetTitleSize(0.04);
      stack->GetXaxis()->SetLabelSize(0.04);
      stack->GetYaxis()->SetTitle(histos_1D[Settings::regOS][Settings::Data][i_fs]->GetYaxis()->GetTitle());
      stack->GetYaxis()->SetTitleSize(0.04);
      stack->GetYaxis()->SetLabelSize(0.04);
		
      stack->GetXaxis()->SetTitleOffset(1.2);
      stack->GetYaxis()->SetTitleOffset(1.25);
		
      histos_1D[Settings::regOS][Settings::Data][i_fs]->Draw("SAME p E1 X0");
		
      TLegend *legend;
      legend  = CreateLegend_2P2F("right",histos_1D[Settings::regOS][Settings::Data][i_fs],histos_1D[Settings::regOS][Settings::WZ][i_fs],histos_1D[Settings::regOS][Settings::qqZZ][i_fs],histos_1D[Settings::regOS][Settings::DY][i_fs],histos_1D[Settings::regOS][Settings::ttbar][i_fs]);
      legend->Draw();
		
      // Draw lumi
      CMS_lumi *lumi = new CMS_lumi;
      lumi->set_lumi(c, _lumi, 0);
		
      TString _out_file_name;
      _out_file_name = folder + "/" + variable_name + "_OS_" + _s_final_state.at(i_fs);// + "_" + _s_category_stxs.at(Settings::inclusive_stxs);
      SavePlots(c, _out_file_name);
		
   }
}
//========================================================================================================


//========================================================================================================
void OSmethod::PlotZXContributions( TString folder )
{
   TCanvas *c, *c_zx;
   TString _out_file_name;
   CMS_lumi *lumi = new CMS_lumi;

   c    = new TCanvas("c", "c", 600, 600);
   c_zx = new TCanvas("c_zx", "c_zx", 600, 600);
	
	for( int i_fs = 0; i_fs <= Settings::fs4l ; i_fs++ )
   {
		//for ( int i_cat = 0; i_cat <= Settings::inclusive_stxs; i_cat++ )
      //{
			c->cd();
			
			h_from3P1F_SR[Settings::nominal][i_fs]       ->SetLineColor(kBlue);
			h_from2P2F_SR[Settings::nominal][i_fs]       ->SetLineColor(kYellow);
			h_from3P1F_SR_final[Settings::nominal][i_fs] ->SetLineColor(kBlack);
			h_from3P1F_SR_ZZonly[Settings::nominal][i_fs]->SetLineColor(kRed);
			histos_ZX[Settings::nominal][i_fs]           ->SetLineColor(kGreen);
			
			h_from3P1F_SR[Settings::nominal][i_fs]->SetMinimum(0.0);
			
			h_from3P1F_SR[Settings::nominal][i_fs]       ->Draw("HIST");
			h_from2P2F_SR[Settings::nominal][i_fs]       ->Draw("HIST SAME");
			h_from3P1F_SR_final[Settings::nominal][i_fs] ->Draw("HIST SAME");
			h_from3P1F_SR_ZZonly[Settings::nominal][i_fs]->Draw("HIST SAME");
			histos_ZX[Settings::nominal][i_fs]           ->Draw("HIST SAME");
			
			TString _fs_label;
		        if ( i_fs == Settings::fs4mu) _fs_label = "m_{4#font[12]{#mu}} (GeV)";
		        if ( i_fs == Settings::fs2muetau) _fs_label = "m_{2#font[12]{#mu}#font[12]{#tau_{e}#tau_{h}}} (GeV)";
		        if ( i_fs == Settings::fs2mumutau) _fs_label = "m_{2#font[12]{#mu}#font[12]{#tau_{#mu}#tau_{h}}} (GeV)";
		        if ( i_fs == Settings::fs2mutautau) _fs_label = "m_{2#font[12]{#mu}#font[12]{#tau_{h}#tau_{h}}} (GeV)";
		        if ( i_fs == Settings::fs2e2mu) _fs_label = "m_{2#font[12]{e}2#font[12]{#mu}} (GeV)";
		        if ( i_fs == Settings::fs2eetau) _fs_label = "m_{2#font[12]{e}#font[12]{#tau_{e}#tau_{h}}} (GeV)";
		        if ( i_fs == Settings::fs2emutau) _fs_label = "m_{2#font[12]{e}#font[12]{#tau_{#mu}#tau_{h}}} (GeV)";
		        if ( i_fs == Settings::fs2etautau) _fs_label = "m_{2#font[12]{e}#font[12]{#tau_{h}#tau_{h}}} (GeV)";
			if ( i_fs == Settings::fs4l)    _fs_label = "m_{4#font[12]{l}} (GeV)";
			h_from3P1F_SR[Settings::nominal][i_fs]->GetXaxis()->SetTitle(_fs_label);
			h_from3P1F_SR[Settings::nominal][i_fs]->GetXaxis()->SetTitleSize(0.04);
			h_from3P1F_SR[Settings::nominal][i_fs]->GetXaxis()->SetLabelSize(0.04);
			h_from3P1F_SR[Settings::nominal][i_fs]->GetYaxis()->SetTitle(h_from2P2F_SR[Settings::nominal][i_fs]->GetYaxis()->GetTitle());
			h_from3P1F_SR[Settings::nominal][i_fs]->GetYaxis()->SetTitleSize(0.04);
			h_from3P1F_SR[Settings::nominal][i_fs]->GetYaxis()->SetLabelSize(0.04);
			
			h_from3P1F_SR[Settings::nominal][i_fs]->GetXaxis()->SetTitleOffset(1.2);
			h_from3P1F_SR[Settings::nominal][i_fs]->GetYaxis()->SetTitleOffset(1.25);
			
			TLegend *legend;
			legend  = CreateLegend_ZXcontr( "right", h_from2P2F_SR[Settings::nominal][i_fs], h_from3P1F_SR[Settings::nominal][i_fs],h_from3P1F_SR_ZZonly[Settings::nominal][i_fs],h_from3P1F_SR_final[Settings::nominal][i_fs],histos_ZX[Settings::nominal][i_fs] );
			legend->Draw();
			
			// Draw lumi
			lumi->set_lumi(c, _lumi, 0);
			
			_out_file_name = folder + "/" + "ZX_Contributions_OS_" + _s_final_state.at(i_fs);// + "_" + _s_category_stxs.at(i_cat);
			SavePlots(c, _out_file_name);
			
			c_zx->cd();
			
			histos_ZX[Settings::nominal][i_fs]->SetLineColor(kGreen+2);
			histos_ZX[Settings::nominal][i_fs]->SetFillColor(kGreen+2);
			histos_ZX[Settings::nominal][i_fs]->Draw("HIST");
			lumi->set_lumi(c_zx, _lumi, 0);
			
			_out_file_name = folder + "/" + "ZX_OS_" + _s_final_state.at(i_fs);// + "_" + _s_category_stxs.at(i_cat);
			SavePlots(c_zx, _out_file_name);
      //}
   }
	
	
}
//========================================================================================================


//========================================================================================================
void OSmethod::FitZX( TString folder )
{
   TCanvas *c_zx;
   CMS_lumi *lumi = new CMS_lumi;
   TF1  *fit_function;
   TString _out_file_name;
   c_zx = new TCanvas("c_zx", "c_zx", 600, 600);
	
	for( int i_fs = 0; i_fs <= Settings::fs4l ; i_fs++ )
   {
		//for ( int i_cat = 0; i_cat <= Settings::inclusive_stxs; i_cat++ )
      //{
			c_zx->cd();
			
			gStyle->SetOptFit();
			gStyle->SetStatY(0.85);
			gStyle->SetStatX(0.95);
			gStyle->SetStatW(0.2);
			gStyle->SetStatH(0.2);

                        fit_function = new TF1("fit_function","[0]*TMath::Landau(x, [1], [2])",70,800);
                        fit_function->SetParNames("Constant","MPV","#sigma");
                        fit_function->SetParameter(0,1.);
                        fit_function->SetParameter(1,100.);
                        fit_function->SetParameter(2,10.);
			
			histos_ZX[Settings::nominal][i_fs]->Fit("fit_function");
			histos_ZX[Settings::nominal][i_fs]->Draw("");
			
			// Draw lumi
			lumi->set_lumi(c_zx, _lumi, 0);
						
			_out_file_name = folder + "/" + "ZX_OS_fit_" + _s_final_state.at(i_fs);// + "_" + _s_category_stxs.at(i_cat);
			SavePlots(c_zx, _out_file_name);
      //}
   }
	gStyle->SetOptFit(0);
}
//========================================================================================================




//===============================================================
void OSmethod::SubtractWZ()
{
   for (int i_flav = 0; i_flav < num_of_flavours; i_flav++)
   {
      passing[Settings::Total][i_flav]->Add(passing[Settings::WZ][i_flav], -1.);
      failing[Settings::Total][i_flav]->Add(failing[Settings::WZ][i_flav], -1.);
   }

   cout << "[INFO] WZ contribution subtracted." << endl;
   
}
//===============================================================

//===============================================================
void OSmethod::ProduceFakeRates( TString file_name )
{
   for(int i_pT_bin = 0; i_pT_bin < _n_pT_bins - 1; i_pT_bin++ )
   {
      double temp_NP = 0;
      double temp_NF = 0;

      double temp_error_NP = 0;
      double temp_error_NF = 0;
      
      for (int i_flav = 0; i_flav < num_of_flavours; i_flav++)
      {
         if ( i_flav == Settings::ele && i_pT_bin == 0) continue; // electrons do not have 5 - 7 GeV bin
	 if ( (i_flav == Settings::tauE || i_flav == Settings::tauMu || i_flav == Settings::tauTau) && i_pT_bin <= 2) continue; //hadronic taus do not have 5 - 20 GeV bin
         temp_NP = passing[Settings::Total][i_flav]->IntegralAndError(passing[Settings::Total][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin]),passing[Settings::Total][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin+1]) - 1, 1, 1, temp_error_NP);
         temp_NF = failing[Settings::Total][i_flav]->IntegralAndError(failing[Settings::Total][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin]),failing[Settings::Total][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin+1]) - 1, 1, 1, temp_error_NF);
         
//         cout << "========================================" << endl;
//         cout << "pT bin = " << _pT_bins[i_pT_bin] << endl;
//         cout << "NP = " << temp_NP << endl;
//         cout << "error NP = " << temp_error_NP << endl;
//         cout << "NF = " << temp_NF << endl;
//         cout << "error NF = " << temp_error_NF << endl;
//         cout << "X = " << (_pT_bins[i_pT_bin] + _pT_bins[i_pT_bin + 1])/2 << endl;
//         cout << "error X = " << (_pT_bins[i_pT_bin + 1] - _pT_bins[i_pT_bin])/2 << endl;
//         cout << "Y = " << temp_NP/(temp_NP+temp_NF) << endl;
//         cout << "error Y = " << sqrt(pow((temp_NF/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NP,2) + pow((temp_NP/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NF,2)) << endl;
         
         vector_X[Settings::corrected][Settings::EB][i_flav].push_back((_pT_bins[i_pT_bin] + _pT_bins[i_pT_bin + 1])/2);
         vector_Y[Settings::corrected][Settings::EB][i_flav].push_back(temp_NP/(temp_NP+temp_NF));
         
         vector_EX[Settings::corrected][Settings::EB][i_flav].push_back((_pT_bins[i_pT_bin + 1] - _pT_bins[i_pT_bin])/2);
         vector_EY[Settings::corrected][Settings::EB][i_flav].push_back(sqrt(pow((temp_NF/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NP,2) + pow((temp_NP/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NF,2)));
         
         temp_NP = passing[Settings::Total][i_flav]->IntegralAndError(passing[Settings::Total][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin]),passing[Settings::Total][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin+1]) - 1, 2, 2, temp_error_NP);
         temp_NF = failing[Settings::Total][i_flav]->IntegralAndError(failing[Settings::Total][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin]),failing[Settings::Total][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin+1]) - 1, 2, 2, temp_error_NF);
         
         vector_X[Settings::corrected][Settings::EE][i_flav].push_back((_pT_bins[i_pT_bin] + _pT_bins[i_pT_bin + 1])/2);
         vector_Y[Settings::corrected][Settings::EE][i_flav].push_back(temp_NP/(temp_NP+temp_NF));
         
         vector_EX[Settings::corrected][Settings::EE][i_flav].push_back((_pT_bins[i_pT_bin + 1] - _pT_bins[i_pT_bin])/2);
         vector_EY[Settings::corrected][Settings::EE][i_flav].push_back(sqrt(pow((temp_NF/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NP,2) + pow((temp_NP/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NF,2)));
         
         // Just for fake rate plots calculate the same for histograms without WZ subtraction
         temp_NP = passing[Settings::Data][i_flav]->IntegralAndError(passing[Settings::Data][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin]),passing[Settings::Data][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin+1]) - 1, 1, 1, temp_error_NP);
         temp_NF = failing[Settings::Data][i_flav]->IntegralAndError(failing[Settings::Data][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin]),failing[Settings::Data][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin+1]) - 1, 1, 1, temp_error_NF);
         
         //         cout << "========================================" << endl;
         //         cout << "pT bin = " << _pT_bins[i_pT_bin] << endl;
         //         cout << "NP = " << temp_NP << endl;
         //         cout << "error NP = " << temp_error_NP << endl;
         //         cout << "NF = " << temp_NF << endl;
         //         cout << "error NF = " << temp_error_NF << endl;
         //         cout << "X = " << (_pT_bins[i_pT_bin] + _pT_bins[i_pT_bin + 1])/2 << endl;
         //         cout << "error X = " << (_pT_bins[i_pT_bin + 1] - _pT_bins[i_pT_bin])/2 << endl;
         //         cout << "Y = " << temp_NP/(temp_NP+temp_NF) << endl;
         //         cout << "error Y = " << sqrt(pow((temp_NF/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NP,2) + pow((temp_NP/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NF,2)) << endl;
         
         vector_X[Settings::uncorrected][Settings::EB][i_flav].push_back((_pT_bins[i_pT_bin] + _pT_bins[i_pT_bin + 1])/2);
         vector_Y[Settings::uncorrected][Settings::EB][i_flav].push_back(temp_NP/(temp_NP+temp_NF));
         
         vector_EX[Settings::uncorrected][Settings::EB][i_flav].push_back((_pT_bins[i_pT_bin + 1] - _pT_bins[i_pT_bin])/2);
         vector_EY[Settings::uncorrected][Settings::EB][i_flav].push_back(sqrt(pow((temp_NF/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NP,2) + pow((temp_NP/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NF,2)));
         
         temp_NP = passing[Settings::Data][i_flav]->IntegralAndError(passing[Settings::Data][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin]),passing[Settings::Data][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin+1]) - 1, 2, 2, temp_error_NP);
         temp_NF = failing[Settings::Data][i_flav]->IntegralAndError(failing[Settings::Data][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin]),failing[Settings::Data][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin+1]) - 1, 2, 2, temp_error_NF);
         
         vector_X[Settings::uncorrected][Settings::EE][i_flav].push_back((_pT_bins[i_pT_bin] + _pT_bins[i_pT_bin + 1])/2);
         vector_Y[Settings::uncorrected][Settings::EE][i_flav].push_back(temp_NP/(temp_NP+temp_NF));
         
         vector_EX[Settings::uncorrected][Settings::EE][i_flav].push_back((_pT_bins[i_pT_bin + 1] - _pT_bins[i_pT_bin])/2);
         vector_EY[Settings::uncorrected][Settings::EE][i_flav].push_back(sqrt(pow((temp_NF/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NP,2) + pow((temp_NP/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NF,2)));

      }
   }
//corrected
   FR_OS_electron_EB = new TGraphErrors (vector_X[Settings::corrected][Settings::EB][Settings::ele].size(),
                                         &(vector_X[Settings::corrected][Settings::EB][Settings::ele][0]),
                                         &(vector_Y[Settings::corrected][Settings::EB][Settings::ele][0]),
                                         &(vector_EX[Settings::corrected][Settings::EB][Settings::ele][0]),
                                         &(vector_EY[Settings::corrected][Settings::EB][Settings::ele][0]));
   FR_OS_electron_EB->SetName("FR_OS_electron_EB");
   
   FR_OS_electron_EE = new TGraphErrors (vector_X[Settings::corrected][Settings::EE][Settings::ele].size(),
                                         &(vector_X[Settings::corrected][Settings::EE][Settings::ele][0]),
                                         &(vector_Y[Settings::corrected][Settings::EE][Settings::ele][0]),
                                         &(vector_EX[Settings::corrected][Settings::EE][Settings::ele][0]),
                                         &(vector_EY[Settings::corrected][Settings::EE][Settings::ele][0]));
   FR_OS_electron_EE->SetName("FR_OS_electron_EE");
   
   FR_OS_muon_EB = new TGraphErrors (vector_X[Settings::corrected][Settings::EB][Settings::mu].size(),
                                     &(vector_X[Settings::corrected][Settings::EB][Settings::mu][0]),
                                     &(vector_Y[Settings::corrected][Settings::EB][Settings::mu][0]),
                                     &(vector_EX[Settings::corrected][Settings::EB][Settings::mu][0]),
                                     &(vector_EY[Settings::corrected][Settings::EB][Settings::mu][0]));
   FR_OS_muon_EB->SetName("FR_OS_muon_EB");
   
   FR_OS_muon_EE = new TGraphErrors (vector_X[Settings::corrected][Settings::EE][Settings::mu].size(),
                                     &(vector_X[Settings::corrected][Settings::EE][Settings::mu][0]),
                                     &(vector_Y[Settings::corrected][Settings::EE][Settings::mu][0]),
                                     &(vector_EX[Settings::corrected][Settings::EE][Settings::mu][0]),
                                     &(vector_EY[Settings::corrected][Settings::EE][Settings::mu][0]));
   FR_OS_muon_EE->SetName("FR_OS_muon_EE");

   FR_OS_tauE_EB = new TGraphErrors (vector_X[Settings::corrected][Settings::EB][Settings::tauE].size(),
                                     &(vector_X[Settings::corrected][Settings::EB][Settings::tauE][0]),
                                     &(vector_Y[Settings::corrected][Settings::EB][Settings::tauE][0]),
                                     &(vector_EX[Settings::corrected][Settings::EB][Settings::tauE][0]),
                                     &(vector_EY[Settings::corrected][Settings::EB][Settings::tauE][0]));
   FR_OS_tauE_EB->SetName("FR_OS_tauE_EB");

   FR_OS_tauE_EE = new TGraphErrors (vector_X[Settings::corrected][Settings::EE][Settings::tauE].size(),
                                     &(vector_X[Settings::corrected][Settings::EE][Settings::tauE][0]),
                                     &(vector_Y[Settings::corrected][Settings::EE][Settings::tauE][0]),
                                     &(vector_EX[Settings::corrected][Settings::EE][Settings::tauE][0]),
                                     &(vector_EY[Settings::corrected][Settings::EE][Settings::tauE][0]));
   FR_OS_tauE_EE->SetName("FR_OS_tauE_EE");

   FR_OS_tauMu_EB = new TGraphErrors (vector_X[Settings::corrected][Settings::EB][Settings::tauMu].size(),
                                     &(vector_X[Settings::corrected][Settings::EB][Settings::tauMu][0]),
                                     &(vector_Y[Settings::corrected][Settings::EB][Settings::tauMu][0]),
                                     &(vector_EX[Settings::corrected][Settings::EB][Settings::tauMu][0]),
                                     &(vector_EY[Settings::corrected][Settings::EB][Settings::tauMu][0]));
   FR_OS_tauMu_EB->SetName("FR_OS_tauMu_EB");

   FR_OS_tauMu_EE = new TGraphErrors (vector_X[Settings::corrected][Settings::EE][Settings::tauMu].size(),
                                     &(vector_X[Settings::corrected][Settings::EE][Settings::tauMu][0]),
                                     &(vector_Y[Settings::corrected][Settings::EE][Settings::tauMu][0]),
                                     &(vector_EX[Settings::corrected][Settings::EE][Settings::tauMu][0]),
                                     &(vector_EY[Settings::corrected][Settings::EE][Settings::tauMu][0]));
   FR_OS_tauMu_EE->SetName("FR_OS_tauMu_EE");

   FR_OS_tauTau_EB = new TGraphErrors (vector_X[Settings::corrected][Settings::EB][Settings::tauTau].size(),
                                     &(vector_X[Settings::corrected][Settings::EB][Settings::tauTau][0]),
                                     &(vector_Y[Settings::corrected][Settings::EB][Settings::tauTau][0]),
                                     &(vector_EX[Settings::corrected][Settings::EB][Settings::tauTau][0]),
                                     &(vector_EY[Settings::corrected][Settings::EB][Settings::tauTau][0]));
   FR_OS_tauTau_EB->SetName("FR_OS_tauTau_EB");

   FR_OS_tauTau_EE = new TGraphErrors (vector_X[Settings::corrected][Settings::EE][Settings::tauTau].size(),
                                     &(vector_X[Settings::corrected][Settings::EE][Settings::tauTau][0]),
                                     &(vector_Y[Settings::corrected][Settings::EE][Settings::tauTau][0]),
                                     &(vector_EX[Settings::corrected][Settings::EE][Settings::tauTau][0]),
                                     &(vector_EY[Settings::corrected][Settings::EE][Settings::tauTau][0]));
   FR_OS_tauTau_EE->SetName("FR_OS_tauTau_EE");   

//uncorrected
   FR_OS_electron_EB_unc = new TGraphErrors (vector_X[Settings::uncorrected][Settings::EB][Settings::ele].size(),
                                         &(vector_X[Settings::uncorrected][Settings::EB][Settings::ele][0]),
                                         &(vector_Y[Settings::uncorrected][Settings::EB][Settings::ele][0]),
                                         &(vector_EX[Settings::uncorrected][Settings::EB][Settings::ele][0]),
                                         &(vector_EY[Settings::uncorrected][Settings::EB][Settings::ele][0]));
   FR_OS_electron_EB_unc->SetName("FR_OS_electron_EB_unc");
   
   FR_OS_electron_EE_unc = new TGraphErrors (vector_X[Settings::uncorrected][Settings::EE][Settings::ele].size(),
                                         &(vector_X[Settings::uncorrected][Settings::EE][Settings::ele][0]),
                                         &(vector_Y[Settings::uncorrected][Settings::EE][Settings::ele][0]),
                                         &(vector_EX[Settings::uncorrected][Settings::EE][Settings::ele][0]),
                                         &(vector_EY[Settings::uncorrected][Settings::EE][Settings::ele][0]));
   FR_OS_electron_EE_unc->SetName("FR_OS_electron_EE_unc");
   
   FR_OS_muon_EB_unc = new TGraphErrors (vector_X[Settings::uncorrected][Settings::EB][Settings::mu].size(),
                                     &(vector_X[Settings::uncorrected][Settings::EB][Settings::mu][0]),
                                     &(vector_Y[Settings::uncorrected][Settings::EB][Settings::mu][0]),
                                     &(vector_EX[Settings::uncorrected][Settings::EB][Settings::mu][0]),
                                     &(vector_EY[Settings::uncorrected][Settings::EB][Settings::mu][0]));
   FR_OS_muon_EB_unc->SetName("FR_OS_muon_EB_unc");
   
   FR_OS_muon_EE_unc = new TGraphErrors (vector_X[Settings::uncorrected][Settings::EE][Settings::mu].size(),
                                     &(vector_X[Settings::uncorrected][Settings::EE][Settings::mu][0]),
                                     &(vector_Y[Settings::uncorrected][Settings::EE][Settings::mu][0]),
                                     &(vector_EX[Settings::uncorrected][Settings::EE][Settings::mu][0]),
                                     &(vector_EY[Settings::uncorrected][Settings::EE][Settings::mu][0]));
   FR_OS_muon_EE_unc->SetName("FR_OS_muon_EE_unc");

   FR_OS_tauE_EB_unc = new TGraphErrors (vector_X[Settings::uncorrected][Settings::EB][Settings::tauE].size(),
                                     &(vector_X[Settings::uncorrected][Settings::EB][Settings::tauE][0]),
                                     &(vector_Y[Settings::uncorrected][Settings::EB][Settings::tauE][0]),
                                     &(vector_EX[Settings::uncorrected][Settings::EB][Settings::tauE][0]),
                                     &(vector_EY[Settings::uncorrected][Settings::EB][Settings::tauE][0]));
   FR_OS_tauE_EB_unc->SetName("FR_OS_tauE_EB_unc");

   FR_OS_tauE_EE_unc = new TGraphErrors (vector_X[Settings::uncorrected][Settings::EE][Settings::tauE].size(),
                                     &(vector_X[Settings::uncorrected][Settings::EE][Settings::tauE][0]),
                                     &(vector_Y[Settings::uncorrected][Settings::EE][Settings::tauE][0]),
                                     &(vector_EX[Settings::uncorrected][Settings::EE][Settings::tauE][0]),
                                     &(vector_EY[Settings::uncorrected][Settings::EE][Settings::tauE][0]));
   FR_OS_tauE_EE_unc->SetName("FR_OS_tauE_EE_unc");

   FR_OS_tauMu_EB_unc = new TGraphErrors (vector_X[Settings::uncorrected][Settings::EB][Settings::tauMu].size(),
                                     &(vector_X[Settings::uncorrected][Settings::EB][Settings::tauMu][0]),
                                     &(vector_Y[Settings::uncorrected][Settings::EB][Settings::tauMu][0]),
                                     &(vector_EX[Settings::uncorrected][Settings::EB][Settings::tauMu][0]),
                                     &(vector_EY[Settings::uncorrected][Settings::EB][Settings::tauMu][0]));
   FR_OS_tauMu_EB_unc->SetName("FR_OS_tauMu_EB_unc");

   FR_OS_tauMu_EE_unc = new TGraphErrors (vector_X[Settings::uncorrected][Settings::EE][Settings::tauMu].size(),
                                     &(vector_X[Settings::uncorrected][Settings::EE][Settings::tauMu][0]),
                                     &(vector_Y[Settings::uncorrected][Settings::EE][Settings::tauMu][0]),
                                     &(vector_EX[Settings::uncorrected][Settings::EE][Settings::tauMu][0]),
                                     &(vector_EY[Settings::uncorrected][Settings::EE][Settings::tauMu][0]));
   FR_OS_tauMu_EE_unc->SetName("FR_OS_tauMu_EE_unc");

   FR_OS_tauTau_EB_unc = new TGraphErrors (vector_X[Settings::uncorrected][Settings::EB][Settings::tauTau].size(),
                                     &(vector_X[Settings::uncorrected][Settings::EB][Settings::tauTau][0]),
                                     &(vector_Y[Settings::uncorrected][Settings::EB][Settings::tauTau][0]),
                                     &(vector_EX[Settings::uncorrected][Settings::EB][Settings::tauTau][0]),
                                     &(vector_EY[Settings::uncorrected][Settings::EB][Settings::tauTau][0]));
   FR_OS_tauTau_EB_unc->SetName("FR_OS_tauTau_EB_unc");

   FR_OS_tauTau_EE_unc = new TGraphErrors (vector_X[Settings::uncorrected][Settings::EE][Settings::tauTau].size(),
                                     &(vector_X[Settings::uncorrected][Settings::EE][Settings::tauTau][0]),
                                     &(vector_Y[Settings::uncorrected][Settings::EE][Settings::tauTau][0]),
                                     &(vector_EX[Settings::uncorrected][Settings::EE][Settings::tauTau][0]),
                                     &(vector_EY[Settings::uncorrected][Settings::EE][Settings::tauTau][0]));
   FR_OS_tauTau_EE_unc->SetName("FR_OS_tauTau_EE_unc");
   

   PlotFR();
   
   TFile* fOutHistos = TFile::Open(file_name, "recreate");
   fOutHistos->cd();
   
   FR_OS_electron_EB->Write();
   FR_OS_electron_EE->Write();
   FR_OS_muon_EB->Write();
   FR_OS_muon_EE->Write();
   FR_OS_tauE_EB->Write();
   FR_OS_tauE_EE->Write();
   FR_OS_tauMu_EB->Write();
   FR_OS_tauMu_EE->Write();
   FR_OS_tauTau_EB->Write();
   FR_OS_tauTau_EE->Write();

   fOutHistos->Close();
   delete fOutHistos;
   
   cout << "[INFO] Fake rates produced and stored in a file." << endl;
}
//===============================================================

//===============================================================
void OSmethod::PlotFR()
{
   TCanvas *c_ele, *c_mu, *c_tauE, *c_tauMu, *c_tauTau;
   c_ele = new TCanvas("FR_ele", "FR_ele", 600, 600);
   c_mu  = new TCanvas("FR_mu", "FR_mu", 600, 600);
   c_tauE = new TCanvas("FR_tauE", "FR_tauMu", 600, 600);
   c_tauMu = new TCanvas("FR_tauMu", "FR_tauMu", 600, 600);
   c_tauTau = new TCanvas("FR_tauTau", "FR_tauTau", 600, 600);
   
   mg_electrons = new TMultiGraph();
   mg_muons = new TMultiGraph();
   mg_tauE = new TMultiGraph();
   mg_tauMu = new TMultiGraph();
   mg_tauTau = new TMultiGraph();

   mg_electrons->Add(FR_OS_electron_EB);
   FR_OS_electron_EB->SetLineColor(kBlue);
   FR_OS_electron_EB->SetLineStyle(2);
   FR_OS_electron_EB->SetMarkerSize(0);
   FR_OS_electron_EB->SetTitle("barel corrected");
   mg_electrons->Add(FR_OS_electron_EE);
   FR_OS_electron_EE->SetLineColor(kRed);
   FR_OS_electron_EE->SetLineStyle(2);
   FR_OS_electron_EE->SetMarkerSize(0);
   FR_OS_electron_EE->SetTitle("endcap corrected");
   mg_electrons->Add(FR_OS_electron_EB_unc);
   FR_OS_electron_EB_unc->SetLineColor(kBlue);
   FR_OS_electron_EB_unc->SetLineStyle(1);
   FR_OS_electron_EB_unc->SetMarkerSize(0);
   FR_OS_electron_EB_unc->SetTitle("barel uncorrected");
   mg_electrons->Add(FR_OS_electron_EE_unc);
   FR_OS_electron_EE_unc->SetLineColor(kRed);
   FR_OS_electron_EE_unc->SetLineStyle(1);
   FR_OS_electron_EE_unc->SetMarkerSize(0);
   FR_OS_electron_EE_unc->SetTitle("endcap uncorrected");
   
   mg_muons->Add(FR_OS_muon_EB);
   FR_OS_muon_EB->SetLineColor(kBlue);
   FR_OS_muon_EB->SetLineStyle(2);
   FR_OS_muon_EB->SetMarkerSize(0);
   FR_OS_muon_EB->SetTitle("barel corrected");
   mg_muons->Add(FR_OS_muon_EE);
   FR_OS_muon_EE->SetLineColor(kRed);
   FR_OS_muon_EE->SetLineStyle(2);
   FR_OS_muon_EE->SetMarkerSize(0);
   FR_OS_muon_EE->SetTitle("endcap corrected");
   mg_muons->Add(FR_OS_muon_EB_unc);
   FR_OS_muon_EB_unc->SetLineColor(kBlue);
   FR_OS_muon_EB_unc->SetLineStyle(1);
   FR_OS_muon_EB_unc->SetMarkerSize(0);
   FR_OS_muon_EB_unc->SetTitle("barel uncorrected");
   mg_muons->Add(FR_OS_muon_EE_unc);
   FR_OS_muon_EE_unc->SetLineColor(kRed);
   FR_OS_muon_EE_unc->SetLineStyle(1);
   FR_OS_muon_EE_unc->SetMarkerSize(0);
   FR_OS_muon_EE_unc->SetTitle("endcap uncorrected");
   
   mg_tauE->Add(FR_OS_tauE_EB);
   FR_OS_tauE_EB->SetLineColor(kBlue);
   FR_OS_tauE_EB->SetLineStyle(2);
   FR_OS_tauE_EB->SetMarkerSize(0);
   FR_OS_tauE_EB->SetTitle("barel corrected");
   mg_tauE->Add(FR_OS_tauE_EE);
   FR_OS_tauE_EE->SetLineColor(kRed);
   FR_OS_tauE_EE->SetLineStyle(2);
   FR_OS_tauE_EE->SetMarkerSize(0);
   FR_OS_tauE_EE->SetTitle("endcap corrected");
   mg_tauE->Add(FR_OS_tauE_EB_unc);
   FR_OS_tauE_EB_unc->SetLineColor(kBlue);
   FR_OS_tauE_EB_unc->SetLineStyle(1);
   FR_OS_tauE_EB_unc->SetMarkerSize(0);
   FR_OS_tauE_EB_unc->SetTitle("barel uncorrected");
   mg_tauE->Add(FR_OS_tauE_EE_unc);
   FR_OS_tauE_EE_unc->SetLineColor(kRed);
   FR_OS_tauE_EE_unc->SetLineStyle(1);
   FR_OS_tauE_EE_unc->SetMarkerSize(0);
   FR_OS_tauE_EE_unc->SetTitle("endcap uncorrected");

   mg_tauMu->Add(FR_OS_tauMu_EB);
   FR_OS_tauMu_EB->SetLineColor(kBlue);
   FR_OS_tauMu_EB->SetLineStyle(2);
   FR_OS_tauMu_EB->SetMarkerSize(0);
   FR_OS_tauMu_EB->SetTitle("barel corrected");
   mg_tauMu->Add(FR_OS_tauMu_EE);
   FR_OS_tauMu_EE->SetLineColor(kRed);
   FR_OS_tauMu_EE->SetLineStyle(2);
   FR_OS_tauMu_EE->SetMarkerSize(0);
   FR_OS_tauMu_EE->SetTitle("endcap corrected");
   mg_tauMu->Add(FR_OS_tauMu_EB_unc);
   FR_OS_tauMu_EB_unc->SetLineColor(kBlue);
   FR_OS_tauMu_EB_unc->SetLineStyle(1);
   FR_OS_tauMu_EB_unc->SetMarkerSize(0);
   FR_OS_tauMu_EB_unc->SetTitle("barel uncorrected");
   mg_tauMu->Add(FR_OS_tauMu_EE_unc);
   FR_OS_tauMu_EE_unc->SetLineColor(kRed);
   FR_OS_tauMu_EE_unc->SetLineStyle(1);
   FR_OS_tauMu_EE_unc->SetMarkerSize(0);
   FR_OS_tauMu_EE_unc->SetTitle("endcap uncorrected");

   mg_tauTau->Add(FR_OS_tauTau_EB);
   FR_OS_tauTau_EB->SetLineColor(kBlue);
   FR_OS_tauTau_EB->SetLineStyle(2);
   FR_OS_tauTau_EB->SetMarkerSize(0);
   FR_OS_tauTau_EB->SetTitle("barel corrected");
   mg_tauTau->Add(FR_OS_tauTau_EE);
   FR_OS_tauTau_EE->SetLineColor(kRed);
   FR_OS_tauTau_EE->SetLineStyle(2);
   FR_OS_tauTau_EE->SetMarkerSize(0);
   FR_OS_tauTau_EE->SetTitle("endcap corrected");
   mg_tauTau->Add(FR_OS_tauTau_EB_unc);
   FR_OS_tauTau_EB_unc->SetLineColor(kBlue);
   FR_OS_tauTau_EB_unc->SetLineStyle(1);
   FR_OS_tauTau_EB_unc->SetMarkerSize(0);
   FR_OS_tauTau_EB_unc->SetTitle("barel uncorrected");
   mg_tauTau->Add(FR_OS_tauTau_EE_unc);
   FR_OS_tauTau_EE_unc->SetLineColor(kRed);
   FR_OS_tauTau_EE_unc->SetLineStyle(1);
   FR_OS_tauTau_EE_unc->SetMarkerSize(0);
   FR_OS_tauTau_EE_unc->SetTitle("endcap uncorrected");  
 
   gStyle->SetEndErrorSize(0);
   
   TLegend *leg_ele,*leg_mu,*leg_tauE,*leg_tauMu,*leg_tauTau;
   CMS_lumi *lumi = new CMS_lumi;

   c_ele->cd();
   lumi->set_lumi(c_ele, _lumi, 0);
   mg_electrons->Draw("AP");
   mg_electrons->GetXaxis()->SetTitle("p_{T} [GeV]");
   mg_electrons->GetYaxis()->SetTitle("Fake Rate");
   mg_electrons->SetTitle("Electron fake rate");
   mg_electrons->SetMaximum(0.35);
   leg_ele = CreateLegend_FR("left",FR_OS_electron_EB_unc,FR_OS_electron_EB,FR_OS_electron_EE_unc,FR_OS_electron_EE);
   leg_ele->Draw();
   SavePlots(c_ele, "Plots/FR_OS_electrons");
   
   c_mu->cd();
   lumi->set_lumi(c_mu, _lumi, 0);
   mg_muons->Draw("AP");
   mg_muons->GetXaxis()->SetTitle("p_{T} [GeV]");
   mg_muons->GetYaxis()->SetTitle("Fake Rate");
   mg_muons->SetTitle("Muon fake rate");
   mg_muons->SetMaximum(0.35);
   leg_mu = CreateLegend_FR("left",FR_OS_muon_EB_unc,FR_OS_muon_EB,FR_OS_muon_EE_unc,FR_OS_muon_EE);
   leg_mu->Draw();
   SavePlots(c_mu, "Plots/FR_OS_muons");

   c_tauE->cd();
   lumi->set_lumi(c_tauE, _lumi, 0);
   mg_tauE->Draw("AP");
   mg_tauE->GetXaxis()->SetTitle("p_{T} [GeV]");
   mg_tauE->GetYaxis()->SetTitle("Fake Rate");
   mg_tauE->SetTitle("#tau fake rate, #tau_{e}#tau_{h} channel");
   mg_tauE->SetMaximum(0.05);
   leg_tauE = CreateLegend_FR("left",FR_OS_tauE_EB_unc,FR_OS_tauE_EB,FR_OS_tauE_EE_unc,FR_OS_tauE_EE);
   leg_tauE->Draw();
   SavePlots(c_tauE, "Plots/FR_OS_tauE");

   c_tauMu->cd();
   lumi->set_lumi(c_tauMu, _lumi, 0);
   mg_tauMu->Draw("AP");
   mg_tauMu->GetXaxis()->SetTitle("p_{T} [GeV]");
   mg_tauMu->GetYaxis()->SetTitle("Fake Rate");
   mg_tauMu->SetTitle("#tau fake rate, #tau_{e}#tau_{h} channel");
   mg_tauMu->SetMaximum(0.05);
   leg_tauMu = CreateLegend_FR("left",FR_OS_tauMu_EB_unc,FR_OS_tauMu_EB,FR_OS_tauMu_EE_unc,FR_OS_tauMu_EE);
   leg_tauMu->Draw();
   SavePlots(c_tauMu, "Plots/FR_OS_tauMu");

   c_tauTau->cd();
   lumi->set_lumi(c_tauTau, _lumi, 0);
   mg_tauTau->Draw("AP");
   mg_tauTau->GetXaxis()->SetTitle("p_{T} [GeV]");
   mg_tauTau->GetYaxis()->SetTitle("Fake Rate");
   mg_tauTau->SetTitle("#tau fake rate, #tau_{e}#tau_{h} channel");
   mg_tauTau->SetMaximum(0.05);
   leg_tauTau = CreateLegend_FR("left",FR_OS_tauTau_EB_unc,FR_OS_tauTau_EB,FR_OS_tauTau_EE_unc,FR_OS_tauTau_EE);
   leg_tauTau->Draw();
   SavePlots(c_tauTau, "Plots/FR_OS_tauTau");  

}
//===============================================================

//===============================================================
void OSmethod::RemoveNegativeBins1D(TH1F *h)
{
   for (int i_bin_x = 1; i_bin_x <= h->GetXaxis()->GetNbins(); i_bin_x++)
   {
      if( h->GetBinContent(i_bin_x) < 0.) h->SetBinContent(i_bin_x, 0);
   }
   
}
//===============================================================

//===============================================================
void OSmethod::RemoveNegativeBins2D(TH2F *h)
{
   for (int i_bin_x = 1; i_bin_x <= h->GetXaxis()->GetNbins(); i_bin_x++)
   {
      for (int i_bin_y = 1; i_bin_y <= h->GetYaxis()->GetNbins(); i_bin_y++)
      {
         if( h->GetBinContent(i_bin_x,i_bin_y) < 0.) h->SetBinContent(i_bin_x,i_bin_y,0);
      }
      
   }
   
}
//===============================================================

//===============================================================
void OSmethod::Set_pT_binning(int size, float *bins)
{
   _n_pT_bins = size;

   for (int i = 0; i < size; i++)
   {
      _pT_bins[i] = bins[i];
   }
}
//===============================================================

//===============================================================
void OSmethod::SetLumi(float lumi)
{
   _lumi = lumi;
}
//===============================================================


//==========================================================
int OSmethod::find_current_process( TString input_file_name )
{
   
   int current_process = -999;
   
   // Assign dataset to correct process
   if ( input_file_name.Contains("Data") )           current_process = Settings::Data;
   if ( input_file_name.Contains("WZ") )             current_process = Settings::WZ;
   if ( input_file_name.Contains("ZZTo4l") )         current_process = Settings::qqZZ;
   if ( input_file_name.Contains("DYJetsToLL") )     current_process = Settings::DY;
   if ( input_file_name.Contains("TTJets") )         current_process = Settings::ttbar;
   if ( input_file_name.Contains("TTTo2L2Nu") )      current_process = Settings::ttbar;
   
   return current_process;
}
//==========================================================


//=============================
int OSmethod::FindFinalState()
{
   int final_state = -999;

   if ( Z1Flav == -121 )
   {
      if ( Z2Flav == -169 ) final_state = Settings::fs2e2mu;
      else if (Z2Flav == -165) final_state = Settings::fs2eetau;
      else if (Z2Flav == -195) final_state = Settings::fs2emutau;
      else if (Z2Flav == -225) final_state = Settings::fs2etautau;
      else
         final_state = -1;//cerr << "[ERROR] in event " << RunNumber << ":" << LumiNumber << ":" << EventNumber << ", Z2Flav = " << Z2Flav << endl;
      }
   else if ( Z1Flav == -169 )
   {
      if ( Z2Flav == -169 ) final_state = Settings::fs4mu;
      else if (Z2Flav == -165) final_state = Settings::fs2muetau;
      else if (Z2Flav == -195) final_state = Settings::fs2mumutau;
      else if (Z2Flav == -225) final_state = Settings::fs2mutautau;
      else
         final_state = -1;//cerr << "[ERROR] in event " << RunNumber << ":" << LumiNumber << ":" << EventNumber << ", Z2Flav = " << Z2Flav << endl;
   }
   else
   {
      final_state = -1;//cerr << "[ERROR] in event " << RunNumber << ":" << LumiNumber << ":" << EventNumber << ", Z1Flav = " << Z1Flav << endl;
   }
   
   return final_state;
}
//=============================


//=================================
float OSmethod::calculate_K_factor(TString input_file_name)
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
int OSmethod::GENMatch(int ileg,TString input_file_name)
{
   if (abs(LepLepId->at(ileg))!=15) {
      cout<<"Only applicable to pdgId = 15, not "<<ileg<<endl;
      return 6;
   }
   if ( input_file_name.Contains("ZZ") || input_file_name.Contains("H") || input_file_name.Contains("ggTo") ) {
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

float OSmethod::calculate_TauIDSF_OneLeg(short DM, float pt, float eta, int genmatch, int flav)
{
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

float OSmethod::calculate_TauIDSF(TString input_file_name)
{
   float SF = 1;
/*   TString WP[3];
   if (Z2Flav==-165) {
      WP[0]="Medium";WP[1]="Medium";WP[2]="Tight";
   }
   else if (Z2Flav==-195 || Z2Flav==-225) {
      WP[0]="Tight";WP[1]="VLoose";WP[2]="Tight";
   }*/
   for (size_t ileg=2;ileg<=3;ileg++) {
      if (abs(LepLepId->at(ileg))!=15)
	 continue;
      int genmatch=GENMatch(ileg,input_file_name);
      SF=SF*calculate_TauIDSF_OneLeg(TauDecayMode->at(ileg),LepPt->at(ileg),LepEta->at(ileg),genmatch,abs(Z2Flav));
   }
   return SF;
}

//===================================================
bool OSmethod::GetVarLogX ( TString variable_name )
{
   //=============
   // M4l
   //=============
   if(variable_name == "M4l")                return bool(Plots::M4l().var_log_x);

   else
   {
      cout << "[ERROR] Wrong variable name choosen!" << endl;
      abort();
      return bool(Plots::M4l().var_log_x);
   }
}
//===================================================



//===================================================
bool OSmethod::GetVarLogY ( TString variable_name )
{
   //=============
   // M4l
   //=============
   if(variable_name == "M4l")                return bool(Plots::M4l().var_log_y);

   else
   {
      cout << "[ERROR] Wrong variable name choosen!" << endl;
      abort();
      return bool(Plots::M4l().var_log_y);
   }
}
//===================================================

//=========================================================================================================
TLegend* OSmethod::CreateLegend_FR( string position, TGraphErrors *EB_unc, TGraphErrors *EB_cor,TGraphErrors *EE_unc,TGraphErrors *EE_cor )
{
   TLegend *leg;
   leg = new TLegend( .64, .65, .97, .9 );
   if(position == "right") leg = new TLegend( .64, .65, .97, .9 );
   else if(position == "left") leg = new TLegend(.18,.65,.51,.9);
   
   leg->AddEntry( EB_unc, "barrel uncorrected", "l" );
   leg->AddEntry( EB_cor, "barrel corrected","l");
   leg->AddEntry( EE_unc, "endcap uncorrected", "l" );
   leg->AddEntry( EE_cor, "endcap corrected", "l" );
   
   return leg;
}
//=========================================================================================================

//=========================================================================================================
TLegend* OSmethod::CreateLegend_ZXcontr( string position, TH1F *h_2P2F_SR, TH1F *h_3P1F_SR,TH1F *h_3P1F_ZZ,TH1F *h_3P1F_SR_final,TH1F *total )
{
   TLegend *leg;
   leg = new TLegend( .64, .65, .97, .9 );
   if(position == "right") leg = new TLegend( .64, .65, .97, .9 );
   else if(position == "left") leg = new TLegend(.18,.65,.51,.9);
   
   leg->AddEntry( h_2P2F_SR, "2P2F", "l" );
   leg->AddEntry( h_3P1F_SR, "3P1F w/o removal","l");
   leg->AddEntry( h_3P1F_ZZ, "3P1F ZZ contr.", "l" );
   leg->AddEntry( h_3P1F_SR_final, "3P1F final", "l" );
   leg->AddEntry( total, "Z+X final", "l" );
   
   return leg;
}
//=========================================================================================================

//=========================================================================================================
TLegend* OSmethod::CreateLegend_2P2F( string position, TH1F *data, TH1F *WZ,TH1F *qqZZ,TH1F *DY,TH1F *ttbar )
{
   TLegend *leg;
   leg = new TLegend( .64, .65, .97, .9 );
   if(position == "right") leg = new TLegend( .64, .65, .97, .9 );
   else if(position == "left") leg = new TLegend(.18,.65,.51,.9);
   leg->SetFillColor(0);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   
   leg->AddEntry( data, "Data", "p" );
   leg->AddEntry( WZ,"WZ","f");
   leg->AddEntry( qqZZ, "Z#gamma*, ZZ", "f" );
   leg->AddEntry( DY, "Z + jets", "f" );
   leg->AddEntry( ttbar, "t#bar{t} + jets", "f" );
   
   return leg;
}
//=========================================================================================================

//=========================================================================================================
TLegend* OSmethod::CreateLegend_3P1F( string position, TH1F *data, TH1F *h_2P2F, TH1F *WZ,TH1F *qqZZ,TH1F *DY,TH1F *ttbar )
{
   TLegend *leg;
   leg = new TLegend( .64, .65, .97, .9 );
   if(position == "right") leg = new TLegend( .64, .65, .97, .9 );
   else if(position == "left") leg = new TLegend(.18,.65,.51,.9);
   leg->SetFillColor(0);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   
   leg->AddEntry( data, "Data", "p" );
   leg->AddEntry( h_2P2F, "2P2F extr.", "l" );
   leg->AddEntry( WZ,"WZ","f");
   leg->AddEntry( qqZZ, "Z#gamma*, ZZ", "f" );
   leg->AddEntry( DY, "Z + jets", "f" );
   leg->AddEntry( ttbar, "t#bar{t} + jets", "f" );
   
   return leg;
}
//=========================================================================================================



//=======================================
void OSmethod::SavePlots( TCanvas *c, TString name)
{
   c->SaveAs(name + ".pdf");
   //c->SaveAs(name + ".root");
   //c->SaveAs(name + ".eps");
   c->SaveAs(name + ".png");
   //gSystem->Exec("convert -density 300 -quality 100 " + name + ".eps " + name + ".png");
}
//=======================================





