// Include classes
#include <HTauTauHMuMu/AnalysisStep/test/ZpXEstimation/include/SSmethod.h>

// Constructor
//============================================================
SSmethod::SSmethod():Tree()
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
   
   _s_region.push_back("ZLL");
   
   // Z+X SS factors
   // Default: ORIGINAL OS/SS ratios evaluated on MC (never recomputed for Run II) 
   // OS/SS ratios evaluated on data and taken when computing FR in SS method
   _fs_ROS_SS.push_back(1.);//4mu
   _fs_ROS_SS.push_back(1.);//2muetau
   _fs_ROS_SS.push_back(1.);//2mumutau
   _fs_ROS_SS.push_back(1.);//2mutautau
   _fs_ROS_SS.push_back(1.);//2e2mu
   _fs_ROS_SS.push_back(1.);//2eetau
   _fs_ROS_SS.push_back(1.);//2emutau
   _fs_ROS_SS.push_back(1.);//2etautau


   vector<float> temp;
   for ( int i_fs = 0; i_fs < num_of_final_states; i_fs++ )
   {
      //for ( int i_cat = 0; i_cat < num_of_categories_stxs; i_cat++ )
      //{
         temp.push_back(0);
         _N_OS_events[i_fs] = 0.;
         _N_SS_events[i_fs] = 0.;
      //}
      _expected_yield_SR.push_back(0);
      _expected_yield_SR_up.push_back(0);
      _expected_yield_SR_dn.push_back(0);
      _number_of_events_CR.push_back(0);
   }
  
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
SSmethod::~SSmethod()
{
}
//====================


//================================================================================================
void SSmethod::Calculate_SSOS_Ratio( TString input_file_data_name, TString input_file_MC_name , bool subtractMC )
{
   input_file_data = TFile::Open( input_file_data_name);
   input_file_MC   = TFile::Open( input_file_MC_name);
   
   hCounters = (TH1F*)input_file_MC->Get("CRZLLTree/Counters");
   gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);
   
   //Loop over data CR to get the number of events in OS and SS
   input_tree_data = (TTree*)input_file_data->Get("CRZLLTree/candTree");
   Init( input_tree_data, input_file_data_name , true);
   _current_process = Settings::Data;
   
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
      
      _current_final_state = FindFinalState();
      if (_current_final_state<0) {continue;}
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
  */    
      if ((test_bit(CRflag, CRZLLss))) _N_SS_events[_current_final_state]+=1.0;
      if ((test_bit(CRflag, CRZLLos_2P2F)) || (test_bit(CRflag, CRZLLos_3P1F))) _N_OS_events[_current_final_state]+=1.0;

      _k_factor = calculate_K_factor(input_file_data_name);
      _TauIDSF = calculate_TauIDSF(input_file_data_name);
      _event_weight = (_lumi * 1000 * xsec * _k_factor * _TauIDSF * overallEventWeight * L1prefiringWeight) / gen_sum_weights;
   }
   
   //Loop over MC to estimate ZZTo4L events in OS
   if( subtractMC )
   {
      input_tree_MC = (TTree*)input_file_MC->Get("CRZLLTree/candTree");
      Init( input_tree_MC, input_file_MC_name , true);
      _current_process = Settings::qqZZ;
   
      if (fChain == 0) return;
   
      nentries = fChain->GetEntriesFast();
      cout<<"[INFO] Processing "<<input_file_MC_name<<", total event: "<<nentries<<endl;

      nbytes = 0, nb = 0;
   
      for (Long64_t jentry=0; jentry<nentries;jentry++)
      {
         Long64_t ientry = LoadTree(jentry);
	 if (ientry%50000==0) cout<<ientry<<endl;
         if (ientry < 0) break;
         nb = fChain->GetEntry(jentry);
         nbytes += nb;
         
         _current_final_state = FindFinalState();
         if (_current_final_state<0) {continue;}
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
  */       
         _k_factor = calculate_K_factor(input_file_data_name);
         _TauIDSF = calculate_TauIDSF(input_file_data_name);
	 _event_weight = (_lumi * 1000 * xsec * _k_factor * _TauIDSF * overallEventWeight * L1prefiringWeight) / gen_sum_weights;
         
         if ((test_bit(CRflag, CRZLLos_2P2F)) || (test_bit(CRflag, CRZLLos_3P1F))) _N_OS_events[_current_final_state]-=_event_weight;
      }
   }
   
   //Calculate inclusive numbers
   //for ( int i_cat = 0; i_cat < num_of_categories_stxs - 1; i_cat++ )
   //{
      for ( int i_fs = 0; i_fs < num_of_final_states - 1; i_fs++ )
      {
         _N_SS_events[Settings::fs4l]          += _N_SS_events[i_fs];   //calculate N events for inclusive 4l final state
         _N_OS_events[Settings::fs4l]          += _N_OS_events[i_fs];
         //_N_SS_events[i_fs] += _N_SS_events[i_fs][i_cat];   //calculate N events for inclusive category
         //_N_OS_events[i_fs][Settings::inclusive_stxs] += _N_OS_events[i_fs][i_cat];
       /*  
         if (false)//( MERGE_2E2MU )
         {
            _N_SS_events[Settings::fs2e2mu][i_cat]    += _N_SS_events[Settings::fs2mu2e][i_cat];   //merge 2e2mu and 2mu2e final states
            _N_OS_events[Settings::fs2e2mu][i_cat]    += _N_OS_events[Settings::fs2mu2e][i_cat];
            _N_SS_events[Settings::fs2mu2e][i_cat]    = 0.;
            _N_OS_events[Settings::fs2mu2e][i_cat]    = 0.;
         }*/
      }
   //}
   //for ( int i_fs = 0; i_fs < num_of_final_states - 1; i_fs++  )
   //{
   //   _N_SS_events[Settings::fs4l] += _N_SS_events[i_fs][Settings::inclusive_stxs];
   //   _N_OS_events[Settings::fs4l][Settings::inclusive_stxs] += _N_OS_events[i_fs][Settings::inclusive_stxs];
   //}
   
   // Print Z + X expected yields for inclusive category
   cout << endl;
   cout << "========================================================================================" << endl;
   cout << "[INFO] Control printout of OS/SS ratio calculation." << endl;
   cout << "========================================================================================" << endl;
   for ( int i_fs = 0; i_fs < num_of_final_states - 1; i_fs++ )
   {
      if (false) continue;//( MERGE_2E2MU && i_fs == Settings::fs2mu2e) continue;
      cout << "Final state: " << _s_final_state.at(i_fs) << endl;
      cout << _N_OS_events[i_fs]/_N_SS_events[i_fs]<< " +/- " << sqrt(1./_N_OS_events[i_fs] + 1./_N_SS_events[i_fs]) << endl;
   }
   
   cout << "[INFO] Total = " << _N_OS_events[Settings::fs4l]/_N_SS_events[Settings::fs4l] << " +/- " <<
   sqrt(1./_N_OS_events[Settings::fs4l] + 1./_N_SS_events[Settings::fs4l]) << endl;
   cout << "========================================================================================" << endl;
   cout << endl;
   
   if(true)
	{
	  for ( int i_fs = 0; i_fs < num_of_final_states - 1; i_fs++ ) _fs_ROS_SS[i_fs] = _N_OS_events[i_fs]/_N_SS_events[i_fs];
	  //_fs_ROS_SS[Settings::fs4mu]   = _N_OS_events[Settings::fs4mu]/_N_SS_events[Settings::fs4mu][Settings::inclusive_stxs];//4mu
	  //_fs_ROS_SS[Settings::fs4e]    = _N_OS_events[Settings::fs4e][Settings::inclusive_stxs]/_N_SS_events[Settings::fs4e][Settings::inclusive_stxs];//4e
	  //_fs_ROS_SS[Settings::fs2e2mu] = _N_OS_events[Settings::fs2e2mu][Settings::inclusive_stxs]/_N_SS_events[Settings::fs2e2mu][Settings::inclusive_stxs];//2e2mu
	  //_fs_ROS_SS[Settings::fs2mu2e] = _N_OS_events[Settings::fs2mu2e][Settings::inclusive_stxs]/_N_SS_events[Settings::fs2mu2e][Settings::inclusive_stxs];//2mu2e
	}

   cout << "[INFO] OS/SS ratios calculated." << endl;
      
}
//================================================================================================

//===============================================================================
void SSmethod::FillFRHistos( TString input_file_data_name )
{
   input_file_data = TFile::Open(input_file_data_name);
   
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
   Int_t _total_events[num_of_final_states];
   Int_t _afterSel_events[num_of_final_states];
   Int_t _failZ1MassCut[num_of_final_states];
   Int_t _failLepPtCut[num_of_final_states];
   Int_t _failEtaCut[num_of_final_states];
   Int_t _failSipVtxCut[num_of_final_states];
   Int_t _failMETCut[num_of_final_states];
   Int_t _passingSelection[num_of_final_states];
   Int_t _faillingSelection[num_of_final_states];
	
	for (int i_fs = 0; i_fs < num_of_final_states; i_fs++)
	{
		_total_events[i_fs] = 0.;
		_afterSel_events[i_fs] = 0.;
		_failZ1MassCut[i_fs] = 0.;
		_failLepPtCut[i_fs] = 0.;
		_failEtaCut[i_fs] = 0.;
		_failSipVtxCut[i_fs] = 0.;
		_failMETCut[i_fs] = 0.;
		_passingSelection[i_fs] = 0.;
		_faillingSelection[i_fs] = 0.;
	}

   for (Long64_t jentry=0; jentry<nentries;jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry%50000==0) cout<<ientry<<endl;
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
            
      if (abs(LepLepId->at(2)) == 11) _total_events[Settings::ele]++;
      else if (abs(LepLepId->at(2)) == 13) _total_events[Settings::mu]++;
      else {_total_events[Settings::tauE]++;_total_events[Settings::tauMu]++;_total_events[Settings::tauTau]++;}
      
      if ( Z1Mass < 40. || Z1Mass > 120.) {//(abs(LepLepId->at(2)) == 11) ? _failZ1MassCut[Settings::ele]++ : _failZ1MassCut[Settings::mu]++; continue;}
	if (abs(LepLepId->at(2)) == 11) _failZ1MassCut[Settings::ele]++;
	else if (abs(LepLepId->at(2)) == 13) _failZ1MassCut[Settings::mu]++;
	else if (abs(LepLepId->at(2)) == 15) {_failZ1MassCut[Settings::tauE]++;_failZ1MassCut[Settings::tauMu]++;_failZ1MassCut[Settings::tauTau]++;}
	continue;
      }
      //else if ( Z1Mass > 120. ) {(fabs(LepLepId->at(2)) == 11) ? _failZ1MassCut[Settings::ele]++ : _failZ1MassCut[Settings::mu]++; continue;}
      else if ( (LepPt->at(0) > LepPt->at(1)) && (LepPt->at(0) < 20. || LepPt->at(1) < 10.) ) {//(fabs(LepLepId->at(2)) == 11) ? _failLepPtCut[Settings::ele]++ : _failLepPtCut[Settings::mu]++; continue;}
        if (abs(LepLepId->at(2)) == 11) _failLepPtCut[Settings::ele]++;
        else if (abs(LepLepId->at(2)) == 13) _failLepPtCut[Settings::mu]++;
        else if (abs(LepLepId->at(2)) == 15) {_failLepPtCut[Settings::tauE]++;_failLepPtCut[Settings::tauMu]++;_failLepPtCut[Settings::tauTau]++;}
        continue;
      }
      else if ( (LepPt->at(1) > LepPt->at(0)) && (LepPt->at(1) < 20. || LepPt->at(0) < 10.) ) {//(fabs(LepLepId->at(2)) == 11) ? _failLepPtCut[Settings::ele]++ : _failLepPtCut[Settings::mu]++; continue;}
        if (abs(LepLepId->at(2)) == 11) _failLepPtCut[Settings::ele]++;
        else if (abs(LepLepId->at(2)) == 13) _failLepPtCut[Settings::mu]++;
        else if (abs(LepLepId->at(2)) == 15) {_failLepPtCut[Settings::tauE]++;_failLepPtCut[Settings::tauMu]++;_failLepPtCut[Settings::tauTau]++;}
        continue;
      }
      else if ( (LepPt->at(2) < 20 && abs(LepLepId->at(2))==15) ) {_failLepPtCut[Settings::tauE]++;_failLepPtCut[Settings::tauMu]++;_failLepPtCut[Settings::tauTau]++;}
      else if ( (fabs(LepEta->at(2)) > 2.5 )) {//(fabs(LepLepId->at(2)) == 11) ? _failEtaCut[Settings::ele]++ : _failEtaCut[Settings::mu]++; continue;}
        if (abs(LepLepId->at(2)) == 11) _failEtaCut[Settings::ele]++;
        else if (abs(LepLepId->at(2)) == 13) _failEtaCut[Settings::mu]++;
        else if (abs(LepLepId->at(2)) == 15) {_failEtaCut[Settings::tauE]++;_failEtaCut[Settings::tauMu]++;_failEtaCut[Settings::tauTau]++;}
        continue;
      }
      else if (abs(LepLepId->at(2))==15 && fabs(LepEta->at(2)) > 2.3) {_failEtaCut[Settings::tauE]++;_failEtaCut[Settings::tauMu]++;_failEtaCut[Settings::tauTau]++;continue;}
      else if ( (LepSIP->at(2) > 4. || Lepdxy->at(2) > 0.5 || Lepdz->at(2) > 1.0) && (abs(LepLepId->at(2)) == 11)) { _failSipVtxCut[Settings::ele]++; continue;} // Included dxy/dz cuts for ele
      else if ( (LepSIP->at(2) > 4. || Lepdxy->at(2) > 0.5 || Lepdz->at(2) > 1.0) && (abs(LepLepId->at(2)) == 13)) { _failSipVtxCut[Settings::mu]++; continue;}  // Included dxy/dz cuts for mu
      else if ( Lepdz->at(2) > 10. && abs(LepLepId->at(2)) == 15 ) {_failSipVtxCut[Settings::tauE]++;_failSipVtxCut[Settings::tauMu]++;_failSipVtxCut[Settings::tauTau]++;continue;}
      // NB: Included SIP cut on muons that was removed when it was included in the muon BDT                                                        

      else if ( PFMET > 25. ) {//(fabs(LepLepId->at(2)) == 11) ? _failMETCut[Settings::ele]++ : _failMETCut[Settings::mu]++; continue;}
        if (abs(LepLepId->at(2)) == 11) _failMETCut[Settings::ele]++;
        else if (abs(LepLepId->at(2)) == 13) _failMETCut[Settings::mu]++;
        else if (abs(LepLepId->at(2)) == 15) {_failMETCut[Settings::tauE]++;_failMETCut[Settings::tauMu]++;_failMETCut[Settings::tauTau]++;}
        continue;
      }

      else if (abs(LepLepId->at(2))==11 || abs(LepLepId->at(2))==13)
	{
	  (abs(LepLepId->at(2)) == 11) ? _afterSel_events[Settings::ele]++ : _afterSel_events[Settings::mu]++;
	  // Final event weight
	  _k_factor = calculate_K_factor(input_file_data_name);
	  _event_weight = (_lumi * 1000 * xsec * _k_factor * overallEventWeight * L1prefiringWeight) / gen_sum_weights;
	  
	  //if( LepisID->at(2) ) // Changed because we are not using BDT-based muon ID but PF+ISO            
          if(LepisID->at(2) && ((abs(LepLepId->at(2)) == 11) ? LepCombRelIsoPF->at(2) < 999999. : LepCombRelIsoPF->at(2) < 0.35))
            {
              (abs(LepLepId->at(2)) == 11) ? _passingSelection[Settings::ele]++ : _passingSelection[Settings::mu]++;
              if(abs(LepLepId->at(2)) == 11 ) passing[_current_process][Settings::ele]->Fill(LepPt->at(2), (fabs(LepEta->at(2)) < 1.479) ? 0.5 : 1.5 , (_current_process == Settings::Data) ? 1 :  _event_weight);
              else if(abs(LepLepId->at(2)) == 13 ) passing[_current_process][Settings::mu]->Fill(LepPt->at(2), (fabs(LepEta->at(2)) < 1.2) ? 0.5 : 1.5 , (_current_process == Settings::Data) ? 1 :  _event_weight);
            }
          else
            {
              (abs(LepLepId->at(2)) == 11) ? _faillingSelection[Settings::ele]++ : _faillingSelection[Settings::mu]++;
              if(abs(LepLepId->at(2)) == 11 ) failing[_current_process][Settings::ele]->Fill(LepPt->at(2), (fabs(LepEta->at(2)) < 1.479) ? 0.5 : 1.5 , (_current_process == Settings::Data) ? 1 :  _event_weight);
              else if(abs(LepLepId->at(2)) == 13 ) failing[_current_process][Settings::mu]->Fill(LepPt->at(2), (fabs(LepEta->at(2)) < 1.2) ? 0.5 : 1.5 , (_current_process == Settings::Data) ? 1 :  _event_weight);
            }
        }
      else// if (abs(LepLepId->at(2))==15)
	{
	  //cout<<"test"<<endl;
	  _afterSel_events[Settings::tauE]++;_afterSel_events[Settings::tauMu]++;_afterSel_events[Settings::tauTau]++;
	  _k_factor = calculate_K_factor(input_file_data_name);
	  _event_weight = (_lumi * 1000 * xsec * _k_factor * overallEventWeight * L1prefiringWeight) / gen_sum_weights;

	  if (LepisID->at(2) && TauVSmu->at(2)>=4 && TauVSe->at(2)>=3 && TauVSjet->at(2)>=6) {
	    _passingSelection[Settings::tauMu]++;_passingSelection[Settings::tauTau]++;
	    _TauIDSF=calculate_TauIDSF_OneLeg(TauDecayMode->at(2), LepPt->at(2), LepEta->at(2), GENMatch(2,input_file_data_name), 225);
	    passing[_current_process][Settings::tauMu]->Fill(LepPt->at(2), (fabs(LepEta->at(2)) < 1.479) ? 0.5 : 1.5 , (_current_process == Settings::Data) ? 1 :  _event_weight*_TauIDSF);
	    passing[_current_process][Settings::tauTau]->Fill(LepPt->at(2), (fabs(LepEta->at(2)) < 1.479) ? 0.5 : 1.5 , (_current_process == Settings::Data) ? 1 :  _event_weight*_TauIDSF);
	  }
	  else {
	    _faillingSelection[Settings::tauMu]++;_faillingSelection[Settings::tauTau]++;
            _TauIDSF=calculate_TauIDSF_OneLeg(TauDecayMode->at(2), LepPt->at(2), LepEta->at(2), GENMatch(2,input_file_data_name), 225);
            failing[_current_process][Settings::tauMu]->Fill(LepPt->at(2), (fabs(LepEta->at(2)) < 1.479) ? 0.5 : 1.5 , (_current_process == Settings::Data) ? 1 :  _event_weight*_TauIDSF);
            failing[_current_process][Settings::tauTau]->Fill(LepPt->at(2), (fabs(LepEta->at(2)) < 1.479) ? 0.5 : 1.5 , (_current_process == Settings::Data) ? 1 :  _event_weight*_TauIDSF);
          }
	  if (LepisID->at(2) && TauVSmu->at(2)>=4 && TauVSe->at(2)>=5 && TauVSjet->at(2)>=5) {
	    _passingSelection[Settings::tauE]++;
	    _TauIDSF=calculate_TauIDSF_OneLeg(TauDecayMode->at(2), LepPt->at(2), LepEta->at(2), GENMatch(2,input_file_data_name), 165);
	    passing[_current_process][Settings::tauE]->Fill(LepPt->at(2), (fabs(LepEta->at(2)) < 1.479) ? 0.5 : 1.5 , (_current_process == Settings::Data) ? 1 :  _event_weight*_TauIDSF);
	  }
	  else {
            _faillingSelection[Settings::tauE]++;
            _TauIDSF=calculate_TauIDSF_OneLeg(TauDecayMode->at(2), LepPt->at(2), LepEta->at(2), GENMatch(2,input_file_data_name), 165);
            failing[_current_process][Settings::tauE]->Fill(LepPt->at(2), (fabs(LepEta->at(2)) < 1.479) ? 0.5 : 1.5 , (_current_process == Settings::Data) ? 1 :  _event_weight*_TauIDSF);
          }

	}
   } // END events loop          
	
	// SS method: control printout for ele/mu in Z+L CR
	if( _current_process == Settings::Data)
	{
		cout << endl;
		cout << "========================================================================================" << endl;
		cout << "[INFO] Control printout for electrons in Z+L control region." << endl;
		cout << "========================================================================================" << endl;
		cout << "[INFO] Total number of events in Z+L control region = " << _total_events[Settings::ele] << endl;
		cout << "[INFO] Events after 40 < Z1 < 120 GeV cut  = " << _total_events[Settings::ele] - _failZ1MassCut[Settings::ele] << endl;
		cout << "[INFO] Events after LepPt > 20,10 GeV cut  = " << _total_events[Settings::ele] - _failZ1MassCut[Settings::ele] - _failLepPtCut[Settings::ele] << endl;
		cout << "[INFO] Events after LepEta cut  = " << _total_events[Settings::ele] - _failZ1MassCut[Settings::ele] - _failLepPtCut[Settings::ele] - _failEtaCut[Settings::ele] << endl;
		cout << "[INFO] Events after SIP < 4 cut + dxy/dz cuts  = " << _total_events[Settings::ele] - _failZ1MassCut[Settings::ele] - _failLepPtCut[Settings::ele] - _failEtaCut[Settings::ele] - _failSipVtxCut[Settings::ele] << endl;
		cout << "[INFO] Events after MET < 25 cut  = " << _total_events[Settings::ele] - _failZ1MassCut[Settings::ele] - _failLepPtCut[Settings::ele] - _failEtaCut[Settings::ele] - _failSipVtxCut[Settings::ele] - _failMETCut[Settings::ele] << endl;
		cout << "[INFO] Total events left = " << _passingSelection[Settings::ele] + _faillingSelection[Settings::ele] << endl;
		cout << "[INFO] Passing selection = " << _passingSelection[Settings::ele]  << endl;
		cout << "[INFO] Failling selection = " << _faillingSelection[Settings::ele] << endl;
		cout << "========================================================================================" << endl;
		cout << endl;
		
		cout << endl;
		cout << "========================================================================================" << endl;
		cout << "[INFO] Control printout for muons in Z+L control region." << endl;
		cout << "========================================================================================" << endl;
		cout << "[INFO] Total number of events in Z+L control region = " << _total_events[Settings::mu] << endl;
		cout << "[INFO] Events after 40 < Z1 < 120 GeV cut  = " << _total_events[Settings::mu] - _failZ1MassCut[Settings::mu] << endl;
		cout << "[INFO] Events after LepPt > 20,10 GeV cut  = " << _total_events[Settings::mu] - _failZ1MassCut[Settings::mu] - _failLepPtCut[Settings::mu] << endl;
		cout << "[INFO] Events after LepEta cut  = " << _total_events[Settings::mu] - _failZ1MassCut[Settings::mu] - _failLepPtCut[Settings::mu] - _failEtaCut[Settings::mu] << endl;
		cout << "[INFO] Events after SIP < 4 cut + dxy/dz cuts  = " << _total_events[Settings::mu] - _failZ1MassCut[Settings::mu] - _failLepPtCut[Settings::mu] - _failEtaCut[Settings::mu] - _failSipVtxCut[Settings::mu] << endl;
                cout << "[INFO] Events after MET < 25 cut  = " << _total_events[Settings::mu] - _failZ1MassCut[Settings::mu] - _failLepPtCut[Settings::mu] - _failEtaCut[Settings::mu] - _failSipVtxCut[Settings::mu] - _failMETCut[Settings::mu] << endl;
		cout << "[INFO] Total events left = " << _passingSelection[Settings::mu] + _faillingSelection[Settings::mu] << endl;
		cout << "[INFO] Passing selection = " << _passingSelection[Settings::mu]  << endl;
		cout << "[INFO] Failling selection = " << _faillingSelection[Settings::mu] << endl;
		cout << "========================================================================================" << endl;
		cout << endl;

                cout << endl;
                cout << "========================================================================================" << endl;
                cout << "[INFO] Control printout for tauEs in Z+L control region." << endl;
                cout << "========================================================================================" << endl;
                cout << "[INFO] Total number of events in Z+L control region = " << _total_events[Settings::tauE] << endl;
                cout << "[INFO] Events after 40 < Z1 < 120 GeV cut  = " << _total_events[Settings::tauE] - _failZ1MassCut[Settings::tauE] << endl;
                cout << "[INFO] Events after LepPt > 20,10 GeV cut  = " << _total_events[Settings::tauE] - _failZ1MassCut[Settings::tauE] - _failLepPtCut[Settings::tauE] << endl;
                cout << "[INFO] Events after LepEta cut  = " << _total_events[Settings::tauE] - _failZ1MassCut[Settings::tauE] - _failLepPtCut[Settings::tauE] - _failEtaCut[Settings::tauE] << endl;
                cout << "[INFO] Events after SIP < 4 cut + dxy/dz cuts  = " << _total_events[Settings::tauE] - _failZ1MassCut[Settings::tauE] - _failLepPtCut[Settings::tauE] - _failEtaCut[Settings::tauE] - _failSipVtxCut[Settings::tauE] << endl;
                cout << "[INFO] Events after MET < 25 cut  = " << _total_events[Settings::tauE] - _failZ1MassCut[Settings::tauE] - _failLepPtCut[Settings::tauE] - _failEtaCut[Settings::tauE] - _failSipVtxCut[Settings::tauE] - _failMETCut[Settings::tauE] << endl;
                cout << "[INFO] Total events left = " << _passingSelection[Settings::tauE] + _faillingSelection[Settings::tauE] << endl;
                cout << "[INFO] Passing selection = " << _passingSelection[Settings::tauE]  << endl;
                cout << "[INFO] Failling selection = " << _faillingSelection[Settings::tauE] << endl;
                cout << "========================================================================================" << endl;
                cout << endl;

                cout << endl;
                cout << "========================================================================================" << endl;
                cout << "[INFO] Control printout for tauMus in Z+L control region." << endl;
                cout << "========================================================================================" << endl;
                cout << "[INFO] Total number of events in Z+L control region = " << _total_events[Settings::tauMu] << endl;
                cout << "[INFO] Events after 40 < Z1 < 120 GeV cut  = " << _total_events[Settings::tauMu] - _failZ1MassCut[Settings::tauMu] << endl;
                cout << "[INFO] Events after LepPt > 20,10 GeV cut  = " << _total_events[Settings::tauMu] - _failZ1MassCut[Settings::tauMu] - _failLepPtCut[Settings::tauMu] << endl;
                cout << "[INFO] Events after LepEta cut  = " << _total_events[Settings::tauMu] - _failZ1MassCut[Settings::tauMu] - _failLepPtCut[Settings::tauMu] - _failEtaCut[Settings::tauMu] << endl;
                cout << "[INFO] Events after SIP < 4 cut + dxy/dz cuts  = " << _total_events[Settings::tauMu] - _failZ1MassCut[Settings::tauMu] - _failLepPtCut[Settings::tauMu] - _failEtaCut[Settings::tauMu] - _failSipVtxCut[Settings::tauMu] << endl;
                cout << "[INFO] Events after MET < 25 cut  = " << _total_events[Settings::tauMu] - _failZ1MassCut[Settings::tauMu] - _failLepPtCut[Settings::tauMu] - _failEtaCut[Settings::tauMu] - _failSipVtxCut[Settings::tauMu] - _failMETCut[Settings::tauMu] << endl;
                cout << "[INFO] Total events left = " << _passingSelection[Settings::tauMu] + _faillingSelection[Settings::tauMu] << endl;
                cout << "[INFO] Passing selection = " << _passingSelection[Settings::tauMu]  << endl;
                cout << "[INFO] Failling selection = " << _faillingSelection[Settings::tauMu] << endl;
                cout << "========================================================================================" << endl;
                cout << endl;

                cout << endl;
                cout << "========================================================================================" << endl;
                cout << "[INFO] Control printout for tauTaus in Z+L control region." << endl;
                cout << "========================================================================================" << endl;
                cout << "[INFO] Total number of events in Z+L control region = " << _total_events[Settings::tauTau] << endl;
                cout << "[INFO] Events after 40 < Z1 < 120 GeV cut  = " << _total_events[Settings::tauTau] - _failZ1MassCut[Settings::tauTau] << endl;
                cout << "[INFO] Events after LepPt > 20,10 GeV cut  = " << _total_events[Settings::tauTau] - _failZ1MassCut[Settings::tauTau] - _failLepPtCut[Settings::tauTau] << endl;
                cout << "[INFO] Events after LepEta cut  = " << _total_events[Settings::tauTau] - _failZ1MassCut[Settings::tauTau] - _failLepPtCut[Settings::tauTau] - _failEtaCut[Settings::tauTau] << endl;
                cout << "[INFO] Events after SIP < 4 cut + dxy/dz cuts  = " << _total_events[Settings::tauTau] - _failZ1MassCut[Settings::tauTau] - _failLepPtCut[Settings::tauTau] - _failEtaCut[Settings::tauTau] - _failSipVtxCut[Settings::tauTau] << endl;
                cout << "[INFO] Events after MET < 25 cut  = " << _total_events[Settings::tauTau] - _failZ1MassCut[Settings::tauTau] - _failLepPtCut[Settings::tauTau] - _failEtaCut[Settings::tauTau] - _failSipVtxCut[Settings::tauTau] - _failMETCut[Settings::tauTau] << endl;
                cout << "[INFO] Total events left = " << _passingSelection[Settings::tauTau] + _faillingSelection[Settings::tauTau] << endl;
                cout << "[INFO] Passing selection = " << _passingSelection[Settings::tauTau]  << endl;
                cout << "[INFO] Failling selection = " << _faillingSelection[Settings::tauTau] << endl;
                cout << "========================================================================================" << endl;
                cout << endl;
	}
	
   cout << "[INFO] Processing of " << input_file_data_name << " done." << endl;
}
//===============================================================================



//===============================================================================
void SSmethod::FillDataMCPlots( TString input_file_data_name )
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
      
      if (!(test_bit(CRflag, CRZLLss))) continue;
      
      //cout << "overallEventWeight=" << overallEventWeight << endl;
      //cout << "dataMCWeight=" << dataMCWeight << endl;
      _current_final_state = FindFinalState();
      if (_current_final_state<0) {continue;}
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
*/
      
      _k_factor = calculate_K_factor(input_file_data_name);
      _TauIDSF = calculate_TauIDSF(input_file_data_name);
      _event_weight = (_lumi * 1000 * xsec * _k_factor * _TauIDSF * overallEventWeight * L1prefiringWeight) / gen_sum_weights;
      //cout << "lumi = " << _lumi << " xsec = " << xsec << " k_factor = " << _k_factor << " SF+PU+GenWeight = " << overallEventWeight << " Sum_Weight = " << gen_sum_weights << endl;
   
      histos_1D[Settings::regZLL][_current_process][_current_final_state]->Fill(ZZMass,(_current_process == Settings::Data) ? 1 :  _event_weight);

   } // END events loop
   
   cout << "[INFO] Processing of " << input_file_data_name << " done." << endl;
}
//===============================================================================



//===============================================================================
void SSmethod::MakeHistogramsZX( TString input_file_data_name, TString  input_file_FR_name )
{
   
   FakeRates *FR = new FakeRates( input_file_FR_name );
   
   input_file_data = TFile::Open( input_file_data_name);
   input_tree_data = (TTree*)input_file_data->Get("CRZLLTree/candTree");
   Init( input_tree_data, input_file_data_name , true);
   
   _current_process = find_current_process(input_file_data_name);
   
   
   if (fChain == 0) return;
   
   Long64_t nentries = fChain->GetEntriesFast();
   cout<<"[INFO] Processing "<<input_file_data_name<<", total event: "<<nentries<<endl;
 
   Long64_t nbytes = 0, nb = 0;   
   //Int_t nevents_CRLLSS = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++)
   {
      
      Long64_t ientry = LoadTree(jentry);
      if (ientry%50000==0) cout<<ientry<<endl;
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      
      if ( !CRflag ) continue;
      if ( !test_bit(CRflag, CRZLLss) ) continue;
      //nevents_CRLLSS += 1;

      // Included SIP and dxy/dz cuts for 3rd and 4th lepton                
      if ( fabs(LepEta->at(2)) > 2.5 || fabs(LepEta->at(3)) > 2.5) {continue;}
      if ( (abs(LepLepId->at(2))==15 && fabs(LepEta->at(2)) > 2.3) || (abs(LepLepId->at(3))==15 && fabs(LepEta->at(3)) > 2.3) ) {continue;}
      if ( (LepSIP->at(2) > 4. || Lepdxy->at(2) > 0.5 || Lepdz->at(2) > 1.0) && (abs(LepLepId->at(2)) == 11 || abs(LepLepId->at(2)) == 13) ) {continue;} // Included dxy/dz cuts for 3rd LEPTON (ele and mu)             
      if ( (LepSIP->at(3) > 4. || Lepdxy->at(3) > 0.5 || Lepdz->at(3) > 1.0) && (abs(LepLepId->at(2)) == 11 || abs(LepLepId->at(2)) == 13) ) {continue;} // Included dxy/dz cuts for 4th LEPTON (ele and mu)
      if ( ( Lepdz->at(2) > 10 || LepPt->at(2) < 20) && abs(LepLepId->at(2))==15 ) {continue;}
      if ( ( Lepdz->at(3) > 10 || LepPt->at(3) < 20) && abs(LepLepId->at(3))==15 ) {continue;}

      if ( ZZMass < 70. ) continue;

      _current_final_state = FindFinalState();
      if (_current_final_state<0) {continue;}
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
      
  */    
      _k_factor = calculate_K_factor(input_file_data_name);
      _TauIDSF = calculate_TauIDSF(input_file_data_name);
      _event_weight = (_lumi * 1000 * xsec * _k_factor * _TauIDSF * overallEventWeight * L1prefiringWeight) / gen_sum_weights;

      int tauChannel=-1;
      if (abs(Z2Flav)==165) tauChannel=0;
      else if (abs(Z2Flav)==195) tauChannel=1;
      else if (abs(Z2Flav)==225) tauChannel=2;

      // Calculate yield
      _yield_SR    = _fs_ROS_SS.at(_current_final_state)*FR->GetFakeRate(LepPt->at(2),LepEta->at(2),LepLepId->at(2),tauChannel)*FR->GetFakeRate(LepPt->at(3),LepEta->at(3),LepLepId->at(3),tauChannel);
      _yield_SR_up = _fs_ROS_SS.at(_current_final_state)*FR->GetFakeRate_Up(LepPt->at(2),LepEta->at(2),LepLepId->at(2),tauChannel)*FR->GetFakeRate_Up(LepPt->at(3),LepEta->at(3),LepLepId->at(3),tauChannel);
      _yield_SR_dn = _fs_ROS_SS.at(_current_final_state)*FR->GetFakeRate_Dn(LepPt->at(2),LepEta->at(2),LepLepId->at(2),tauChannel)*FR->GetFakeRate_Dn(LepPt->at(3),LepEta->at(3),LepLepId->at(3),tauChannel);
      
      //if (_yield_SR<0) cout<<_yield_SR<<","<<tauChannel<<","<<_fs_ROS_SS.at(_current_final_state)<<","<<FR->GetFakeRate(LepPt->at(2),LepEta->at(2),LepLepId->at(2),tauChannel)<<","<<FR->GetFakeRate(LepPt->at(3),LepEta->at(3),LepLepId->at(3),tauChannel)<<endl;
      _expected_yield_SR[_current_final_state]    += _yield_SR;
      _expected_yield_SR_up[_current_final_state] += _yield_SR_up;
      _expected_yield_SR_dn[_current_final_state] += _yield_SR_dn;
      _number_of_events_CR[_current_final_state]++;
		
      // Fill m4l Z+X histograms
      histos_ZX[Settings::regZLL][_current_process][_current_final_state]->Fill(ZZMass,(_current_process == Settings::Data) ? _yield_SR :  _yield_SR*_event_weight);
      

   } // End events loop

   //std::cout << "###################################################\n";
   //std::cout << "# events CRLLSS = " << nevents_CRLLSS << '\n';   
   //std::cout << "###################################################\n";


   //for (  int i_cat = 0; i_cat < num_of_categories_stxs - 1; i_cat++  )
   //{
      for ( int i_fs = 0; i_fs < num_of_final_states - 1; i_fs++  )
      {
         _expected_yield_SR[Settings::fs4l]             += _expected_yield_SR[i_fs];   //calculate expected yield for inclusive 4l final state
         _number_of_events_CR[Settings::fs4l]           += _number_of_events_CR[i_fs];
         //_expected_yield_SR[i_fs][Settings::inclusive_stxs]    += _expected_yield_SR[i_fs][i_cat];   //calculate expected yield for inclusive category
         //_expected_yield_SR_up[i_fs][Settings::inclusive_stxs] += _expected_yield_SR_up[i_fs][i_cat];
         //_expected_yield_SR_dn[i_fs][Settings::inclusive_stxs] += _expected_yield_SR_dn[i_fs][i_cat];
         //_number_of_events_CR[i_fs][Settings::inclusive_stxs]  += _number_of_events_CR[i_fs][i_cat];
         
      /*   if (false)//( MERGE_2E2MU )
         {
            _expected_yield_SR[Settings::fs2e2mu][i_cat]       += _expected_yield_SR[Settings::fs2mu2e][i_cat];   //merge 2e2mu and 2mu2e final states
            _expected_yield_SR_up[Settings::fs2e2mu][i_cat]    += _expected_yield_SR_up[Settings::fs2mu2e][i_cat];
            _expected_yield_SR_dn[Settings::fs2e2mu][i_cat]    += _expected_yield_SR_dn[Settings::fs2mu2e][i_cat];
            _number_of_events_CR[Settings::fs2e2mu][i_cat]     += _number_of_events_CR[Settings::fs2mu2e][i_cat];
            _expected_yield_SR[Settings::fs2mu2e][i_cat]       = 0.;
            _expected_yield_SR_up[Settings::fs2mu2e][i_cat]    = 0.;
            _expected_yield_SR_dn[Settings::fs2mu2e][i_cat]    = 0.;
            _number_of_events_CR[Settings::fs2mu2e][i_cat]     = 0.;
         }*/
      }
   //}
   for ( int i_fs = 0; i_fs < num_of_final_states - 1; i_fs++  )
   {
      _expected_yield_SR[Settings::fs4l]    += _expected_yield_SR[i_fs];
      _expected_yield_SR_up[Settings::fs4l] += _expected_yield_SR_up[i_fs];
      _expected_yield_SR_dn[Settings::fs4l] += _expected_yield_SR_dn[i_fs];
   }
   
   // Print Z + X expected yields and uncertainties
   cout << endl;
   cout << "===================================================================================================================================" << endl;
   cout << "[INFO] Control printout for Z+X yields in final states derived with this fake rate file: " << input_file_FR_name << endl;
   cout << "===================================================================================================================================" << endl;
   for ( int i_fs = 0; i_fs < num_of_final_states - 1; i_fs++ )
   {
      if (false) continue;//( MERGE_2E2MU && i_fs == Settings::fs2mu2e) continue;
      //for ( int i_cat = 0; i_cat < num_of_categories_stxs; i_cat++)
      //{
        float stat    = _expected_yield_SR[i_fs]/sqrt(_number_of_events_CR[i_fs]); //statistical uncertainty
        float syst_up = _expected_yield_SR[i_fs]*((_expected_yield_SR_up[i_fs]/_expected_yield_SR[i_fs]) - 1.); //systematic uncertainty due to fake rate variation
        float syst_dn = _expected_yield_SR[i_fs]*(1. - (_expected_yield_SR_dn[i_fs]/_expected_yield_SR[i_fs]));
        float syst_comp = _expected_yield_SR[i_fs]*0.3; //background composition uncertainty of 30% measured in Run I

        float comb_up = sqrt(stat*stat + syst_up*syst_up + syst_comp*syst_comp);
        float comb_dn = sqrt(stat*stat + syst_dn*syst_dn + syst_comp*syst_comp);
			
			
	cout << "Final state: " << _s_final_state.at(i_fs) << endl;
	cout << _expected_yield_SR[i_fs] << " +/- " << comb_dn << "/" << comb_up << "(total):" << "  - " << stat << " (stat., evt: " <<_number_of_events_CR[i_fs] << ")" 
             << "   - " << syst_dn << "/" << syst_up << " (syst.)" <<  "   Relative uncertainty = " << (1. - comb_dn/_expected_yield_SR[i_fs]) << "/" << (1. + comb_up/_expected_yield_SR[i_fs]) << endl;
      //}
      cout << "==================================================================================================================================" << endl;
   }
   
   cout << "[INFO] Total = " << _expected_yield_SR[Settings::fs4l] << endl;
   cout << "==================================================================================================================================" << endl;
   cout << endl;	
   
   cout << "[INFO] Z+X histograms filled." << endl;
}
//===============================================================================




//===============================================================
void SSmethod::DeclareFRHistos()
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
void SSmethod::DeclareDataMCHistos()
{
   for (int i_reg = 0; i_reg < num_of_regions_ss; i_reg ++)
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
void SSmethod::DeclareZXHistos()
{
   for (int i_reg = 0; i_reg < num_of_regions_ss; i_reg ++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         for (int i_fs = 0; i_fs < num_of_final_states; i_fs++)
         {
            //for (int i_cat = 0; i_cat < num_of_categories_stxs; i_cat++)
            //{
               _histo_name = "ZX_M4l_" + _s_region.at(i_reg) + "_" + _s_process.at(i_proc) + "_" + _s_final_state.at(i_fs);// + "_" + _s_category_stxs.at(i_cat);
               _histo_labels = ";" + Plots::M4l().var_X_label + ";" + Plots::M4l().var_Y_label;
               histos_ZX[i_reg][i_proc][i_fs] = new TH1F(_histo_name, _histo_labels, Plots::M4l().var_N_bin, Plots::M4l().var_min, Plots::M4l().var_max);
            //}
         }
      }
   }
}
//===============================================================

//===============================================================
void SSmethod::SaveFRHistos( TString file_name, bool subtractWZ, bool remove_negative_bins)
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
void SSmethod::SaveDataMCHistos( TString file_name )
{
   FillDataMCInclusive();
   
   TFile* fOutHistos = TFile::Open(file_name, "recreate");
   fOutHistos->cd();
   
   for (int i_reg = 0; i_reg < num_of_regions_ss; i_reg ++)
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
void SSmethod::FillDataMCInclusive( )
{
   for (int i_reg = 0; i_reg < num_of_regions_ss; i_reg ++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         for (int i_fs = 0; i_fs < Settings::fs4l; i_fs++)
         {
            //for (int i_cat = 0; i_cat < Settings::inclusive_stxs; i_cat++)
            //{
               //histos_1D[i_reg][i_proc][i_fs][Settings::inclusive_stxs]->Add(histos_1D[i_reg][i_proc][i_fs][i_cat]);
               histos_1D[i_reg][i_proc][Settings::fs4l]    ->Add(histos_1D[i_reg][i_proc][i_fs]);
            //}
         }
      }
   }
/*   
   for (int i_reg = 0; i_reg < num_of_regions_ss; i_reg ++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         for (int i_fs = 0; i_fs < Settings::fs4l; i_fs++)
         {
            histos_1D[i_reg][i_proc][Settings::fs4l][Settings::inclusive_stxs]->Add(histos_1D[i_reg][i_proc][i_fs][Settings::inclusive_stxs]);
         }
      }
   }
*/   
   cout << "[INFO] All Data/MC histograms summed." << endl;
}
//===============================================================

//===============================================================
void SSmethod::SaveZXHistos( TString file_name )
{
   FillZXInclusive();
   
   TFile* fOutHistos = TFile::Open(file_name, "recreate");
   fOutHistos->cd();
   
   for (int i_reg = 0; i_reg < num_of_regions_ss; i_reg ++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         for (int i_fs = 0; i_fs < num_of_final_states; i_fs++)
         {
            //for (int i_cat = 0; i_cat < num_of_categories_stxs; i_cat++)
            //{
               histos_ZX[i_reg][i_proc][i_fs]->Write();
            //}
         }
      }
   }
   
   fOutHistos->Close();
   delete fOutHistos;
   
   cout << "[INFO] All Z+X histograms saved." << endl;
}
//===============================================================

//===============================================================
void SSmethod::FillZXInclusive( )
{
   for (int i_reg = 0; i_reg < num_of_regions_ss; i_reg ++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         for (int i_fs = 0; i_fs < Settings::fs4l; i_fs++)
         {
            //for (int i_cat = 0; i_cat < Settings::inclusive_stxs; i_cat++)
            //{
               //histos_ZX[i_reg][i_proc][i_fs][Settings::inclusive_stxs]->Add(histos_ZX[i_reg][i_proc][i_fs][i_cat]);
               histos_ZX[i_reg][i_proc][Settings::fs4l]    ->Add(histos_ZX[i_reg][i_proc][i_fs]);
            //}
         }
      }
   }
   /*
   for (int i_reg = 0; i_reg < num_of_regions_ss; i_reg ++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         for (int i_fs = 0; i_fs < Settings::fs4l; i_fs++)
         {
            histos_ZX[i_reg][i_proc][Settings::fs4l][Settings::inclusive_stxs]->Add(histos_ZX[i_reg][i_proc][i_fs][Settings::inclusive_stxs]);
         }
      }
   }
 */  
   cout << "[INFO] All Z+X histograms summed." << endl;
}
//===============================================================

//===============================================================
void SSmethod::GetFRHistos( TString file_name)
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
void SSmethod::GetDataMCHistos( TString file_name)
{
   TFile* histo_file = TFile::Open(file_name);
   
   for (int i_reg = 0; i_reg < num_of_regions_ss; i_reg ++)
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
void SSmethod::GetZXHistos( TString file_name)
{
   TFile* histo_file = TFile::Open(file_name);
   
   for (int i_reg = 0; i_reg < num_of_regions_ss; i_reg ++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         for (int i_fs = 0; i_fs < num_of_final_states; i_fs++)
         {
            //for (int i_cat = 0; i_cat < num_of_categories_stxs; i_cat++)
            //{
               _histo_name = "ZX_M4l_" + _s_region.at(i_reg) + "_" + _s_process.at(i_proc) + "_" + _s_final_state.at(i_fs);// + "_" + _s_category_stxs.at(i_cat);
               histos_ZX[i_reg][i_proc][i_fs] = (TH1F*)histo_file->Get(_histo_name);
            //}
         }
      }
   }
   
   cout << "[INFO] All Z+X histograms retrieved from file." << endl;
}

//===============================================================


//===============================================================
void SSmethod::ProduceFakeRates( TString file_name , TString input_file_data_name /*= "DONT_CORRECT"*/)
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
	 if ( (i_flav == Settings::tauE || i_flav == Settings::tauMu || i_flav == Settings::tauTau) && i_pT_bin <= 2) continue;// hadronic taus do not have 5 - 20 GeV bin
         temp_NP = passing[Settings::Total][i_flav]->IntegralAndError(passing[Settings::Total][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin]),passing[Settings::Total][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin+1]) - 1, 1, 1, temp_error_NP);
         temp_NF = failing[Settings::Total][i_flav]->IntegralAndError(failing[Settings::Total][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin]),failing[Settings::Total][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin+1]) - 1, 1, 1, temp_error_NF);
        
	 //cout<<i_flav<<": "<<temp_NP<<", "<<temp_NF<<endl; 
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
         
         //cout<<i_flav<<": "<<temp_NP<<", "<<temp_NF<<endl;

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
   FR_SS_muon_EB = new TGraphErrors (vector_X[Settings::corrected][Settings::EB][Settings::mu].size(),
                                     &(vector_X[Settings::corrected][Settings::EB][Settings::mu][0]),
                                     &(vector_Y[Settings::corrected][Settings::EB][Settings::mu][0]),
                                     &(vector_EX[Settings::corrected][Settings::EB][Settings::mu][0]),
                                     &(vector_EY[Settings::corrected][Settings::EB][Settings::mu][0]));
   FR_SS_muon_EB->SetName("FR_SS_muon_EB");
   
   FR_SS_muon_EE = new TGraphErrors (vector_X[Settings::corrected][Settings::EE][Settings::mu].size(),
                                     &(vector_X[Settings::corrected][Settings::EE][Settings::mu][0]),
                                     &(vector_Y[Settings::corrected][Settings::EE][Settings::mu][0]),
                                     &(vector_EX[Settings::corrected][Settings::EE][Settings::mu][0]),
                                     &(vector_EY[Settings::corrected][Settings::EE][Settings::mu][0]));
   FR_SS_muon_EE->SetName("FR_SS_muon_EE");

   FR_SS_tauE_EB = new TGraphErrors (vector_X[Settings::corrected][Settings::EB][Settings::tauE].size(),
                                     &(vector_X[Settings::corrected][Settings::EB][Settings::tauE][0]),
                                     &(vector_Y[Settings::corrected][Settings::EB][Settings::tauE][0]),
                                     &(vector_EX[Settings::corrected][Settings::EB][Settings::tauE][0]),
                                     &(vector_EY[Settings::corrected][Settings::EB][Settings::tauE][0]));
   FR_SS_tauE_EB->SetName("FR_SS_tauE_EB");

   FR_SS_tauE_EE = new TGraphErrors (vector_X[Settings::corrected][Settings::EE][Settings::tauE].size(),
                                     &(vector_X[Settings::corrected][Settings::EE][Settings::tauE][0]),
                                     &(vector_Y[Settings::corrected][Settings::EE][Settings::tauE][0]),
                                     &(vector_EX[Settings::corrected][Settings::EE][Settings::tauE][0]),
                                     &(vector_EY[Settings::corrected][Settings::EE][Settings::tauE][0]));
   FR_SS_tauE_EE->SetName("FR_SS_tauE_EE");

   FR_SS_tauMu_EB = new TGraphErrors (vector_X[Settings::corrected][Settings::EB][Settings::tauMu].size(),
                                     &(vector_X[Settings::corrected][Settings::EB][Settings::tauMu][0]),
                                     &(vector_Y[Settings::corrected][Settings::EB][Settings::tauMu][0]),
                                     &(vector_EX[Settings::corrected][Settings::EB][Settings::tauMu][0]),
                                     &(vector_EY[Settings::corrected][Settings::EB][Settings::tauMu][0]));
   FR_SS_tauMu_EB->SetName("FR_SS_tauMu_EB");

   FR_SS_tauMu_EE = new TGraphErrors (vector_X[Settings::corrected][Settings::EE][Settings::tauMu].size(),
                                     &(vector_X[Settings::corrected][Settings::EE][Settings::tauMu][0]),
                                     &(vector_Y[Settings::corrected][Settings::EE][Settings::tauMu][0]),
                                     &(vector_EX[Settings::corrected][Settings::EE][Settings::tauMu][0]),
                                     &(vector_EY[Settings::corrected][Settings::EE][Settings::tauMu][0]));
   FR_SS_tauMu_EE->SetName("FR_SS_tauMu_EE");

   FR_SS_tauTau_EB = new TGraphErrors (vector_X[Settings::corrected][Settings::EB][Settings::tauTau].size(),
                                     &(vector_X[Settings::corrected][Settings::EB][Settings::tauTau][0]),
                                     &(vector_Y[Settings::corrected][Settings::EB][Settings::tauTau][0]),
                                     &(vector_EX[Settings::corrected][Settings::EB][Settings::tauTau][0]),
                                     &(vector_EY[Settings::corrected][Settings::EB][Settings::tauTau][0]));
   FR_SS_tauTau_EB->SetName("FR_SS_tauTau_EB");

   FR_SS_tauTau_EE = new TGraphErrors (vector_X[Settings::corrected][Settings::EE][Settings::tauTau].size(),
                                     &(vector_X[Settings::corrected][Settings::EE][Settings::tauTau][0]),
                                     &(vector_Y[Settings::corrected][Settings::EE][Settings::tauTau][0]),
                                     &(vector_EX[Settings::corrected][Settings::EE][Settings::tauTau][0]),
                                     &(vector_EY[Settings::corrected][Settings::EE][Settings::tauTau][0]));
   FR_SS_tauTau_EE->SetName("FR_SS_tauTau_EE");   

//uncorrected
   FR_SS_electron_EB_unc = new TGraphErrors (vector_X[Settings::uncorrected][Settings::EB][Settings::ele].size(),
                                         &(vector_X[Settings::uncorrected][Settings::EB][Settings::ele][0]),
                                         &(vector_Y[Settings::uncorrected][Settings::EB][Settings::ele][0]),
                                         &(vector_EX[Settings::uncorrected][Settings::EB][Settings::ele][0]),
                                         &(vector_EY[Settings::uncorrected][Settings::EB][Settings::ele][0]));
   FR_SS_electron_EB_unc->SetName("FR_SS_electron_EB_unc");
   
   FR_SS_electron_EE_unc = new TGraphErrors (vector_X[Settings::uncorrected][Settings::EE][Settings::ele].size(),
                                         &(vector_X[Settings::uncorrected][Settings::EE][Settings::ele][0]),
                                         &(vector_Y[Settings::uncorrected][Settings::EE][Settings::ele][0]),
                                         &(vector_EX[Settings::uncorrected][Settings::EE][Settings::ele][0]),
                                         &(vector_EY[Settings::uncorrected][Settings::EE][Settings::ele][0]));
   FR_SS_electron_EE_unc->SetName("FR_SS_electron_EE_unc");
   
   FR_SS_muon_EB_unc = new TGraphErrors (vector_X[Settings::uncorrected][Settings::EB][Settings::mu].size(),
                                     &(vector_X[Settings::uncorrected][Settings::EB][Settings::mu][0]),
                                     &(vector_Y[Settings::uncorrected][Settings::EB][Settings::mu][0]),
                                     &(vector_EX[Settings::uncorrected][Settings::EB][Settings::mu][0]),
                                     &(vector_EY[Settings::uncorrected][Settings::EB][Settings::mu][0]));
   FR_SS_muon_EB_unc->SetName("FR_SS_muon_EB_unc");
   
   FR_SS_muon_EE_unc = new TGraphErrors (vector_X[Settings::uncorrected][Settings::EE][Settings::mu].size(),
                                     &(vector_X[Settings::uncorrected][Settings::EE][Settings::mu][0]),
                                     &(vector_Y[Settings::uncorrected][Settings::EE][Settings::mu][0]),
                                     &(vector_EX[Settings::uncorrected][Settings::EE][Settings::mu][0]),
                                     &(vector_EY[Settings::uncorrected][Settings::EE][Settings::mu][0]));
   FR_SS_muon_EE_unc->SetName("FR_SS_muon_EE_unc");

   FR_SS_tauE_EB_unc = new TGraphErrors (vector_X[Settings::uncorrected][Settings::EB][Settings::tauE].size(),
                                     &(vector_X[Settings::uncorrected][Settings::EB][Settings::tauE][0]),
                                     &(vector_Y[Settings::uncorrected][Settings::EB][Settings::tauE][0]),
                                     &(vector_EX[Settings::uncorrected][Settings::EB][Settings::tauE][0]),
                                     &(vector_EY[Settings::uncorrected][Settings::EB][Settings::tauE][0]));
   FR_SS_tauE_EB_unc->SetName("FR_SS_tauE_EB_unc");

   FR_SS_tauE_EE_unc = new TGraphErrors (vector_X[Settings::uncorrected][Settings::EE][Settings::tauE].size(),
                                     &(vector_X[Settings::uncorrected][Settings::EE][Settings::tauE][0]),
                                     &(vector_Y[Settings::uncorrected][Settings::EE][Settings::tauE][0]),
                                     &(vector_EX[Settings::uncorrected][Settings::EE][Settings::tauE][0]),
                                     &(vector_EY[Settings::uncorrected][Settings::EE][Settings::tauE][0]));
   FR_SS_tauE_EE_unc->SetName("FR_SS_tauE_EE_unc");

   FR_SS_tauMu_EB_unc = new TGraphErrors (vector_X[Settings::uncorrected][Settings::EB][Settings::tauMu].size(),
                                     &(vector_X[Settings::uncorrected][Settings::EB][Settings::tauMu][0]),
                                     &(vector_Y[Settings::uncorrected][Settings::EB][Settings::tauMu][0]),
                                     &(vector_EX[Settings::uncorrected][Settings::EB][Settings::tauMu][0]),
                                     &(vector_EY[Settings::uncorrected][Settings::EB][Settings::tauMu][0]));
   FR_SS_tauMu_EB_unc->SetName("FR_SS_tauMu_EB_unc");

   FR_SS_tauMu_EE_unc = new TGraphErrors (vector_X[Settings::uncorrected][Settings::EE][Settings::tauMu].size(),
                                     &(vector_X[Settings::uncorrected][Settings::EE][Settings::tauMu][0]),
                                     &(vector_Y[Settings::uncorrected][Settings::EE][Settings::tauMu][0]),
                                     &(vector_EX[Settings::uncorrected][Settings::EE][Settings::tauMu][0]),
                                     &(vector_EY[Settings::uncorrected][Settings::EE][Settings::tauMu][0]));
   FR_SS_tauMu_EE_unc->SetName("FR_SS_tauMu_EE_unc");

   FR_SS_tauTau_EB_unc = new TGraphErrors (vector_X[Settings::uncorrected][Settings::EB][Settings::tauTau].size(),
                                     &(vector_X[Settings::uncorrected][Settings::EB][Settings::tauTau][0]),
                                     &(vector_Y[Settings::uncorrected][Settings::EB][Settings::tauTau][0]),
                                     &(vector_EX[Settings::uncorrected][Settings::EB][Settings::tauTau][0]),
                                     &(vector_EY[Settings::uncorrected][Settings::EB][Settings::tauTau][0]));
   FR_SS_tauTau_EB_unc->SetName("FR_SS_tauTau_EB_unc");

   FR_SS_tauTau_EE_unc = new TGraphErrors (vector_X[Settings::uncorrected][Settings::EE][Settings::tauTau].size(),
                                     &(vector_X[Settings::uncorrected][Settings::EE][Settings::tauTau][0]),
                                     &(vector_Y[Settings::uncorrected][Settings::EE][Settings::tauTau][0]),
                                     &(vector_EX[Settings::uncorrected][Settings::EE][Settings::tauTau][0]),
                                     &(vector_EY[Settings::uncorrected][Settings::EE][Settings::tauTau][0]));
   FR_SS_tauTau_EE_unc->SetName("FR_SS_tauTau_EE_unc");	

   // Electron fake rates must be corrected using average number of missing hits
   //if ( input_file_data_name != "DONT_CORRECT" ) CorrectElectronFakeRate(input_file_data_name);
	
   FR_SS_electron_EB = new TGraphErrors (vector_X[Settings::corrected][Settings::EB][Settings::ele].size(),
				     &(vector_X[Settings::corrected][Settings::EB][Settings::ele][0]),
				     &(vector_Y[Settings::corrected][Settings::EB][Settings::ele][0]),
				     &(vector_EX[Settings::corrected][Settings::EB][Settings::ele][0]),
				     &(vector_EY[Settings::corrected][Settings::EB][Settings::ele][0]));
   FR_SS_electron_EB->SetName("FR_SS_electron_EB");
	
   FR_SS_electron_EE = new TGraphErrors (vector_X[Settings::corrected][Settings::EE][Settings::ele].size(),
				     &(vector_X[Settings::corrected][Settings::EE][Settings::ele][0]),
				     &(vector_Y[Settings::corrected][Settings::EE][Settings::ele][0]),
				     &(vector_EX[Settings::corrected][Settings::EE][Settings::ele][0]),
				     &(vector_EY[Settings::corrected][Settings::EE][Settings::ele][0]));
   FR_SS_electron_EE->SetName("FR_SS_electron_EE");


   PlotFR();
   
   TFile* fOutHistos = TFile::Open(file_name, "recreate");
   fOutHistos->cd();
   
   FR_SS_electron_EB->Write();
   FR_SS_electron_EE->Write();
   FR_SS_muon_EB->Write();
   FR_SS_muon_EE->Write();
   FR_SS_tauE_EB->Write();
   FR_SS_tauE_EE->Write();
   FR_SS_tauMu_EB->Write();
   FR_SS_tauMu_EE->Write();
   FR_SS_tauTau_EB->Write();
   FR_SS_tauTau_EE->Write();
   
   FR_SS_electron_EB_unc->Write();
   FR_SS_electron_EE_unc->Write();
   FR_SS_muon_EB_unc->Write();
   FR_SS_muon_EE_unc->Write();
   FR_SS_tauE_EB_unc->Write();
   FR_SS_tauE_EE_unc->Write();
   FR_SS_tauMu_EB_unc->Write();
   FR_SS_tauMu_EE_unc->Write();
   FR_SS_tauTau_EB_unc->Write();
   FR_SS_tauTau_EE_unc->Write();
   
   fOutHistos->Close();
   delete fOutHistos;
   
   cout << "[INFO] Fake rates produced and stored in a file." << endl;
}
//===============================================================


//========================================================================
void SSmethod::CorrectElectronFakeRate( TString input_file_data_name )
{
	TGraphErrors *FR_MissingHits_graph[num_of_eta_bins][99];
	
	Calculate_FR_nMissingHits(input_file_data_name, FR_MissingHits_graph);
	Fit_FRnMH_graphs(FR_MissingHits_graph);
	cout << "[INFO] All graphs fitted." << endl;
	Correct_Final_FR( input_file_data_name );
	cout << "[INFO] Electron fake rates have been corrected." << endl;
}
//========================================================================


//========================================================================
void SSmethod::Calculate_FR_nMissingHits( TString input_file_data_name, TGraphErrors *FR_MissingHits_graph[99][99] )
{
	input_file_data = TFile::Open( input_file_data_name);
	
	hCounters = (TH1F*)input_file_data->Get("CRZLTree/Counters");
	gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);
	
	input_tree_data = (TTree*)input_file_data->Get("CRZLTree/candTree");
	Init( input_tree_data, input_file_data_name , false);
	
	_current_process = find_current_process(input_file_data_name);
	
	for ( int i_pt = 0; i_pt < _n_pT_bins-2; i_pt++)
	{
		for ( int i_eta = 0; i_eta < num_of_eta_bins; i_eta++)
		{
			for ( int i_ZMass = 0; i_ZMass < num_of_z_mass_windows; i_ZMass++ )
			{
				_N_MissingHits[i_ZMass][i_eta][i_pt] = 0.;
				_N_Passing[i_ZMass][i_eta][i_pt] = 0.;
				_N_Failling[i_ZMass][i_eta][i_pt] = 0.;
			}
		}
	}
	
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
		
		if ( abs(LepLepId->at(2)) != 11 ) continue; // only electrons
		if ( (LepPt->at(0) > LepPt->at(1)) && (LepPt->at(0) < 20. || LepPt->at(1) < 10.) ) continue;
		if ( (LepPt->at(1) > LepPt->at(0)) && (LepPt->at(1) < 20. || LepPt->at(0) < 10.) ) continue;
		if ( LepSIP->at(2) > 4. || Lepdxy->at(2) > 0.5 || Lepdz->at(2) > 1.0) continue;
		if ( PFMET > 25. ) continue;
		else
		{
			_current_pT_bin = Find_Ele_pT_bin ( LepPt->at(2) );
			_current_eta_bin = Find_Ele_eta_bin ( LepEta->at(2));
			
			if ( fabs( Z1Mass - 91.1876 ) > 5. )//( (Z1Mass > 40.) && (Z1Mass < 120.) )
			{
				_N_MissingHits[Settings::_Window5_10][_current_eta_bin][_current_pT_bin] += LepMissingHit->at(2);
				if(LepisID->at(2) ) _N_Passing[Settings::_Window5_10][_current_eta_bin][_current_pT_bin] += 1.;
				else _N_Failling[Settings::_Window5_10][_current_eta_bin][_current_pT_bin] += 1.;
			}
			
			else if ( fabs( Z1Mass - 91.1876 ) > 2. )
			{
				_N_MissingHits[Settings::_Window2_5][_current_eta_bin][_current_pT_bin] += LepMissingHit->at(2);
				if(LepisID->at(2) ) _N_Passing[Settings::_Window2_5][_current_eta_bin][_current_pT_bin] += 1.;
				else _N_Failling[Settings::_Window2_5][_current_eta_bin][_current_pT_bin] += 1.;
			}
			
			else //if ( fabs( Z1Mass - 91.1876 ) > 2. )//( (Z1Mass > 60.) && (Z1Mass < 120.) )
			{
				_N_MissingHits[Settings::_Window0_2][_current_eta_bin][_current_pT_bin] += LepMissingHit->at(2);
				if(LepisID->at(2) ) _N_Passing[Settings::_Window0_2][_current_eta_bin][_current_pT_bin] += 1.;
				else _N_Failling[Settings::_Window0_2][_current_eta_bin][_current_pT_bin] += 1.;
			}
			/*
			TLorentzVector p1,p2,p3;
			p1.SetPtEtaPhiM(LepPt->at(0), LepEta->at(0), LepPhi->at(0), 0.);
			p2.SetPtEtaPhiM(LepPt->at(1), LepEta->at(1), LepPhi->at(1), 0.);
			p3.SetPtEtaPhiM(LepPt->at(2), LepEta->at(2), LepPhi->at(2), 0.);
			
			else //if ( fabs( ((p1+p2)+p3).M() - 91.2 ) < 5. )//3 lepton mass
			{
				_N_MissingHits[Settings::_Window0_2][_current_eta_bin][_current_pT_bin] += LepMissingHit->at(2);
				if(LepisID->at(2) ) _N_Passing[Settings::_Window0_2][_current_eta_bin][_current_pT_bin] += 1.;
				else _N_Failling[Settings::_Window0_2][_current_eta_bin][_current_pT_bin] += 1.;
			}*/
		}
	} // END events loop
	
	//Fill vectors to produce TGraphs
	
	vector<Float_t> vector_x[num_of_eta_bins][_n_pT_bins-2];
	vector<Float_t> vector_y[num_of_eta_bins][_n_pT_bins-2];
	vector<Float_t> vector_ex[num_of_eta_bins][_n_pT_bins-2];
	vector<Float_t> vector_ey[num_of_eta_bins][_n_pT_bins-2];
	
	for ( int i_pt = 0; i_pt < _n_pT_bins-2; i_pt++)
	{
		for ( int i_eta = 0; i_eta < num_of_eta_bins; i_eta++)
		{
			for ( int i_ZMass = 0; i_ZMass < num_of_z_mass_windows; i_ZMass++ )
			{
				cout<<i_eta<<"\t"<<i_pt<<"\t"<<i_ZMass<<"\t"<<_N_MissingHits[i_ZMass][i_eta][i_pt]/(_N_Passing[i_ZMass][i_eta][i_pt] + _N_Failling[i_ZMass][i_eta][i_pt])<<"\t"<<_N_Passing[i_ZMass][i_eta][i_pt]/(_N_Passing[i_ZMass][i_eta][i_pt] + _N_Failling[i_ZMass][i_eta][i_pt])<<endl;
				vector_x[i_eta][i_pt].push_back(_N_MissingHits[i_ZMass][i_eta][i_pt]/(_N_Passing[i_ZMass][i_eta][i_pt] + _N_Failling[i_ZMass][i_eta][i_pt]));
				vector_ex[i_eta][i_pt].push_back(sqrt(pow((1./pow(_N_Failling[i_ZMass][i_eta][i_pt]+_N_Passing[i_ZMass][i_eta][i_pt],1)),2)*_N_MissingHits[i_ZMass][i_eta][i_pt] + pow((_N_MissingHits[i_ZMass][i_eta][i_pt]/pow(_N_Failling[i_ZMass][i_eta][i_pt]+_N_Passing[i_ZMass][i_eta][i_pt],2)),2)*(_N_Failling[i_ZMass][i_eta][i_pt]+_N_Passing[i_ZMass][i_eta][i_pt])));
				
				vector_y[i_eta][i_pt].push_back(_N_Passing[i_ZMass][i_eta][i_pt]/(_N_Passing[i_ZMass][i_eta][i_pt] + _N_Failling[i_ZMass][i_eta][i_pt]));
				vector_ey[i_eta][i_pt].push_back(sqrt(pow((_N_Failling[i_ZMass][i_eta][i_pt]/pow(_N_Failling[i_ZMass][i_eta][i_pt]+_N_Passing[i_ZMass][i_eta][i_pt],2)),2)*_N_Passing[i_ZMass][i_eta][i_pt] + pow((_N_Passing[i_ZMass][i_eta][i_pt]/pow(_N_Failling[i_ZMass][i_eta][i_pt]+_N_Passing[i_ZMass][i_eta][i_pt],2)),2)*_N_Failling[i_ZMass][i_eta][i_pt]));
				
//				cout << "========================================" << endl;
//				cout << "[INFO] Z+L printout." << endl;
//				cout << "========================================" << endl;
//				cout << "Control region number: " << i_ZMass << endl;
//				cout << "pT bin = " << _pT_bins[i_pt+1] << " - " <<  _pT_bins[i_pt + 2] << endl;
//				cout << "eta bin = " << i_eta << endl;
//				cout << "NP = " << _N_Passing[i_ZMass][i_eta][i_pt] << endl;
//				cout << "NF = " << _N_Failling[i_ZMass][i_eta][i_pt] << endl;
//				cout << "MH = " << _N_MissingHits[i_ZMass][i_eta][i_pt] << endl;
//				cout << "avg_MH = " << vector_x[i_eta][i_pt][i_ZMass] << endl;
//				cout << "avg_MH error = " << vector_ex[i_eta][i_pt][i_ZMass] << endl;
//				cout << "FR = " << vector_y[i_eta][i_pt][i_ZMass] << endl;
//				cout << "FR error = " << vector_ey[i_eta][i_pt][i_ZMass] << endl;
			}
		}
		
	}
	
	for ( int i_pt = 0; i_pt < _n_pT_bins-2; i_pt++)
	{
		for ( int i_eta = 0; i_eta < num_of_eta_bins; i_eta++)
		{
			FR_MissingHits_graph[i_eta][i_pt] = new TGraphErrors(vector_x[i_eta][i_pt].size(),
																					&(vector_x[i_eta][i_pt][0]),
																					&(vector_y[i_eta][i_pt][0]),
																					&(vector_ex[i_eta][i_pt][0]),
																					&(vector_ey[i_eta][i_pt][0]));
		}
		
	}
	
}
//========================================================================


//============================================================================
void SSmethod::Fit_FRnMH_graphs(TGraphErrors *FR_MissingHits_graph[99][99])
{
	for ( int i_pt = 0; i_pt < _n_pT_bins-2; i_pt++)
	{
		for ( int i_eta = 0; i_eta < num_of_eta_bins; i_eta++)
		{
			TString func_name;
			func_name.Form("FR_MissingHits_func_eta_%d_pT_%d",i_eta,i_pt);
			Ele_FR_correction_function[i_eta][i_pt] = new TF1(func_name,"[0]*x+[1]",0,3);
			Ele_FR_correction_function[i_eta][i_pt]->SetParameter(0,1.);
			Ele_FR_correction_function[i_eta][i_pt]->SetParameter(1,0.);
			
			FR_MissingHits_graph[i_eta][i_pt]->Fit(Ele_FR_correction_function[i_eta][i_pt], "Q");
			
			TString graph_name;
			graph_name.Form("FR_MissingHits_graph_eta_%d_pT_%d",i_eta,i_pt);
			FR_MissingHits_graph[i_eta][i_pt]->SetName(graph_name);
			FR_MissingHits_graph[i_eta][i_pt]->GetXaxis()->SetTitle("<# Missing Hits>");
			FR_MissingHits_graph[i_eta][i_pt]->GetYaxis()->SetTitle("Fake Rate");
			TCanvas *c1 = new TCanvas(graph_name,graph_name,900,900);
			c1->cd();
			FR_MissingHits_graph[i_eta][i_pt]->Draw("AP");
			system("mkdir -p Fits");
			SavePlots(c1, "Fits/" + graph_name);
		}
	}
}
//============================================================================


//=============================================================
void SSmethod::Correct_Final_FR( TString input_file_data_name)
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
	
	float _N_MissingHits_ZLL[num_of_eta_bins][_n_pT_bins];
	float _N_Passing_ZLL[num_of_eta_bins][_n_pT_bins];
	float _N_Failling_ZLL[num_of_eta_bins][_n_pT_bins];
	
	float _avg_MissingHits_ZLL[num_of_eta_bins][_n_pT_bins];
	
	for ( int i_pt = 0; i_pt <= _n_pT_bins-2; i_pt++)
	{
		for ( int i_eta = 0; i_eta < num_of_eta_bins; i_eta++)
		{
			_N_MissingHits_ZLL[i_eta][i_pt] = 0.;
			_N_Passing_ZLL[i_eta][i_pt] = 0.;
			_N_Failling_ZLL[i_eta][i_pt] = 0.;
		}
	}
	
	for (Long64_t jentry=0; jentry<nentries;jentry++)
	{
		Long64_t ientry = LoadTree(jentry);
		if (ientry%50000==0) cout<<ientry<<endl;
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);
		nbytes += nb;
		
		if (!(test_bit(CRflag, CRZLLss))) continue;
		
		if ( abs(Z2Flav) != 121) continue; // only electrons
		if ( abs(Z1Flav) != 121) continue; // only 4e

		else
		{
			_current_pT_bin = Find_Ele_pT_bin ( LepPt->at(2) );
			_current_eta_bin = Find_Ele_eta_bin ( LepEta->at(2));

			_N_MissingHits_ZLL[_current_eta_bin][_current_pT_bin] += LepMissingHit->at(2);
			if(LepisID->at(2) ) _N_Passing_ZLL[_current_eta_bin][_current_pT_bin] += 1.;
			else _N_Failling_ZLL[_current_eta_bin][_current_pT_bin] += 1.;
			
			_current_pT_bin = Find_Ele_pT_bin ( LepPt->at(3) );
			_current_eta_bin = Find_Ele_eta_bin ( LepEta->at(3));
			
			_N_MissingHits_ZLL[_current_eta_bin][_current_pT_bin] += LepMissingHit->at(3);
			if(LepisID->at(3) ) _N_Passing_ZLL[_current_eta_bin][_current_pT_bin] += 1.;
			else _N_Failling_ZLL[_current_eta_bin][_current_pT_bin] += 1.;
		}
		
	} // END events loop
	
	for ( int i_pt = 0; i_pt < _n_pT_bins-2; i_pt++)
	{
		Float_t sigma_avgMH = 0;
		_avg_MissingHits_ZLL[Settings::EB][i_pt] = _N_MissingHits_ZLL[Settings::EB][i_pt]/(_N_Passing_ZLL[Settings::EB][i_pt] + _N_Failling_ZLL[Settings::EB][i_pt]);
		sigma_avgMH = sqrt(pow((1./pow(_N_Failling_ZLL[Settings::EB][i_pt]+_N_Passing_ZLL[Settings::EB][i_pt],1)),2)*_N_MissingHits_ZLL[Settings::EB][i_pt] + pow((_N_MissingHits_ZLL[Settings::EB][i_pt]/pow(_N_Failling_ZLL[Settings::EB][i_pt]+_N_Passing_ZLL[Settings::EB][i_pt],2)),2)*(_N_Failling_ZLL[Settings::EB][i_pt]+_N_Passing_ZLL[Settings::EB][i_pt]));
		
		vector_X[Settings::corrected][Settings::EB][Settings::ele][i_pt] = ((_pT_bins[i_pt + 1] + _pT_bins[i_pt + 2])/2);
		vector_Y[Settings::corrected][Settings::EB][Settings::ele][i_pt] = (Ele_FR_correction_function[Settings::EB][i_pt]->Eval(_avg_MissingHits_ZLL[Settings::EB][i_pt]));
		
		vector_EX[Settings::corrected][Settings::EB][Settings::ele][i_pt] = ((_pT_bins[i_pt + 2] - _pT_bins[i_pt + 1])/2);
		vector_EY[Settings::corrected][Settings::EB][Settings::ele][i_pt] = (Ele_FR_correction_function[Settings::EB][i_pt]->Eval(_avg_MissingHits_ZLL[Settings::EB][i_pt]) - Ele_FR_correction_function[Settings::EB][i_pt]->Eval(_avg_MissingHits_ZLL[Settings::EB][i_pt] - sigma_avgMH));
		
		
		_avg_MissingHits_ZLL[Settings::EE][i_pt] = _N_MissingHits_ZLL[Settings::EE][i_pt]/(_N_Passing_ZLL[Settings::EE][i_pt] + _N_Failling_ZLL[Settings::EE][i_pt]);
		sigma_avgMH = sqrt(pow((1./pow(_N_Failling_ZLL[Settings::EE][i_pt]+_N_Passing_ZLL[Settings::EE][i_pt],1)),2)*_N_MissingHits_ZLL[Settings::EE][i_pt] + pow((_N_MissingHits_ZLL[Settings::EE][i_pt]/pow(_N_Failling_ZLL[Settings::EE][i_pt]+_N_Passing_ZLL[Settings::EE][i_pt],2)),2)*(_N_Failling_ZLL[Settings::EE][i_pt]+_N_Passing_ZLL[Settings::EE][i_pt]));
		
		vector_X[Settings::corrected][Settings::EE][Settings::ele][i_pt] = ((_pT_bins[i_pt + 1] + _pT_bins[i_pt + 2])/2);
		vector_Y[Settings::corrected][Settings::EE][Settings::ele][i_pt] = (Ele_FR_correction_function[Settings::EE][i_pt]->Eval(_avg_MissingHits_ZLL[Settings::EE][i_pt]));
		
		vector_EX[Settings::corrected][Settings::EE][Settings::ele][i_pt] = ((_pT_bins[i_pt + 2] - _pT_bins[i_pt + 1])/2);
		vector_EY[Settings::corrected][Settings::EE][Settings::ele][i_pt] = (Ele_FR_correction_function[Settings::EE][i_pt]->Eval(_avg_MissingHits_ZLL[Settings::EE][i_pt]) - Ele_FR_correction_function[Settings::EE][i_pt]->Eval(_avg_MissingHits_ZLL[Settings::EE][i_pt] - sigma_avgMH));
		
//		cout << "========================================" << endl;
//		cout << "[INFO] Z+LL printout." << endl;
//		cout << "========================================" << endl;
//		cout << "pT bin = " << _pT_bins[i_pt + 1] << " - " <<  _pT_bins[i_pt + 2] << endl;
//		cout << "eta bin = " << Settings::EB << endl;
//		cout << "NP = " << _N_Passing_ZLL[Settings::EB][i_pt] << endl;
//		cout << "NF = " << _N_Failling_ZLL[Settings::EB][i_pt] << endl;
//		cout << "avg_MH = " << _avg_MissingHits_ZLL[Settings::EB][i_pt] << endl;
//		cout << "FR = " << _N_Passing_ZLL[Settings::EB][i_pt]/(_N_Passing_ZLL[Settings::EB][i_pt]+_N_Failling_ZLL[Settings::EB][i_pt]) << endl;
//		cout << "corr FR = " << vector_Y[Settings::corrected][Settings::EB][Settings::ele][i_pt] << endl;
//		cout << "========================================" << endl;
//		cout << "pT bin = " << _pT_bins[i_pt + 1] << " - " <<  _pT_bins[i_pt + 2] << endl;
//		cout << "eta bin = " << Settings::EE << endl;
//		cout << "NP = " << _N_Passing_ZLL[Settings::EE][i_pt] << endl;
//		cout << "NF = " << _N_Failling_ZLL[Settings::EE][i_pt] << endl;
//		cout << "avg_MH = " << _avg_MissingHits_ZLL[Settings::EE][i_pt] << endl;
//		cout << "FR = " << _N_Passing_ZLL[Settings::EE][i_pt]/(_N_Passing_ZLL[Settings::EE][i_pt]+_N_Failling_ZLL[Settings::EE][i_pt]) << endl;
//		cout << "corr FR = " << vector_Y[Settings::corrected][Settings::EE][Settings::ele][i_pt] << endl;
	}
	
}
//=============================================================


//===============================================================
void SSmethod::SubtractWZ()
{
   passing[Settings::Total][Settings::mu]->Add(passing[Settings::WZ][Settings::mu], -1.);
   failing[Settings::Total][Settings::mu]->Add(failing[Settings::WZ][Settings::mu], -1.);
   
   cout << "[INFO] WZ contribution subtracted." << endl;
   
}
//===============================================================

//===============================================================
void SSmethod::PlotFR()
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
   
   mg_electrons->Add(FR_SS_electron_EB);
   FR_SS_electron_EB->SetLineColor(kBlue);
   FR_SS_electron_EB->SetLineStyle(2);
   FR_SS_electron_EB->SetMarkerSize(0);
   FR_SS_electron_EB->SetTitle("barrel corrected");
   mg_electrons->Add(FR_SS_electron_EE);
   FR_SS_electron_EE->SetLineColor(kRed);
   FR_SS_electron_EE->SetLineStyle(2);
   FR_SS_electron_EE->SetMarkerSize(0);
   FR_SS_electron_EE->SetTitle("endcap corrected");
   mg_electrons->Add(FR_SS_electron_EB_unc);
   FR_SS_electron_EB_unc->SetLineColor(kBlue);
   FR_SS_electron_EB_unc->SetLineStyle(1);
   FR_SS_electron_EB_unc->SetMarkerSize(0);
   FR_SS_electron_EB_unc->SetTitle("barrel uncorrected");
   mg_electrons->Add(FR_SS_electron_EE_unc);
   FR_SS_electron_EE_unc->SetLineColor(kRed);
   FR_SS_electron_EE_unc->SetLineStyle(1);
   FR_SS_electron_EE_unc->SetMarkerSize(0);
   FR_SS_electron_EE_unc->SetTitle("endcap uncorrected");
   
   mg_muons->Add(FR_SS_muon_EB);
   FR_SS_muon_EB->SetLineColor(kBlue);
   FR_SS_muon_EB->SetLineStyle(2);
   FR_SS_muon_EB->SetMarkerSize(0);
   FR_SS_muon_EB->SetTitle("barrel corrected");
   mg_muons->Add(FR_SS_muon_EE);
   FR_SS_muon_EE->SetLineColor(kRed);
   FR_SS_muon_EE->SetLineStyle(2);
   FR_SS_muon_EE->SetMarkerSize(0);
   FR_SS_muon_EE->SetTitle("endcap corrected");
   mg_muons->Add(FR_SS_muon_EB_unc);
   FR_SS_muon_EB_unc->SetLineColor(kBlue);
   FR_SS_muon_EB_unc->SetLineStyle(1);
   FR_SS_muon_EB_unc->SetMarkerSize(0);
   FR_SS_muon_EB_unc->SetTitle("barrel uncorrected");
   mg_muons->Add(FR_SS_muon_EE_unc);
   FR_SS_muon_EE_unc->SetLineColor(kRed);
   FR_SS_muon_EE_unc->SetLineStyle(1);
   FR_SS_muon_EE_unc->SetMarkerSize(0);
   FR_SS_muon_EE_unc->SetTitle("endcap uncorrected");

   mg_tauE->Add(FR_SS_tauE_EB);
   FR_SS_tauE_EB->SetLineColor(kBlue);
   FR_SS_tauE_EB->SetLineStyle(2);
   FR_SS_tauE_EB->SetMarkerSize(0);
   FR_SS_tauE_EB->SetTitle("barel corrected");
   mg_tauE->Add(FR_SS_tauE_EE);
   FR_SS_tauE_EE->SetLineColor(kRed);
   FR_SS_tauE_EE->SetLineStyle(2);
   FR_SS_tauE_EE->SetMarkerSize(0);
   FR_SS_tauE_EE->SetTitle("endcap corrected");
   mg_tauE->Add(FR_SS_tauE_EB_unc);
   FR_SS_tauE_EB_unc->SetLineColor(kBlue);
   FR_SS_tauE_EB_unc->SetLineStyle(1);
   FR_SS_tauE_EB_unc->SetMarkerSize(0);
   FR_SS_tauE_EB_unc->SetTitle("barel uncorrected");
   mg_tauE->Add(FR_SS_tauE_EE_unc);
   FR_SS_tauE_EE_unc->SetLineColor(kRed);
   FR_SS_tauE_EE_unc->SetLineStyle(1);
   FR_SS_tauE_EE_unc->SetMarkerSize(0);
   FR_SS_tauE_EE_unc->SetTitle("endcap uncorrected");

   mg_tauMu->Add(FR_SS_tauMu_EB);
   FR_SS_tauMu_EB->SetLineColor(kBlue);
   FR_SS_tauMu_EB->SetLineStyle(2);
   FR_SS_tauMu_EB->SetMarkerSize(0);
   FR_SS_tauMu_EB->SetTitle("barel corrected");
   mg_tauMu->Add(FR_SS_tauMu_EE);
   FR_SS_tauMu_EE->SetLineColor(kRed);
   FR_SS_tauMu_EE->SetLineStyle(2);
   FR_SS_tauMu_EE->SetMarkerSize(0);
   FR_SS_tauMu_EE->SetTitle("endcap corrected");
   mg_tauMu->Add(FR_SS_tauMu_EB_unc);
   FR_SS_tauMu_EB_unc->SetLineColor(kBlue);
   FR_SS_tauMu_EB_unc->SetLineStyle(1);
   FR_SS_tauMu_EB_unc->SetMarkerSize(0);
   FR_SS_tauMu_EB_unc->SetTitle("barel uncorrected");
   mg_tauMu->Add(FR_SS_tauMu_EE_unc);
   FR_SS_tauMu_EE_unc->SetLineColor(kRed);
   FR_SS_tauMu_EE_unc->SetLineStyle(1);
   FR_SS_tauMu_EE_unc->SetMarkerSize(0);
   FR_SS_tauMu_EE_unc->SetTitle("endcap uncorrected");

   mg_tauTau->Add(FR_SS_tauTau_EB);
   FR_SS_tauTau_EB->SetLineColor(kBlue);
   FR_SS_tauTau_EB->SetLineStyle(2);
   FR_SS_tauTau_EB->SetMarkerSize(0);
   FR_SS_tauTau_EB->SetTitle("barel corrected");
   mg_tauTau->Add(FR_SS_tauTau_EE);
   FR_SS_tauTau_EE->SetLineColor(kRed);
   FR_SS_tauTau_EE->SetLineStyle(2);
   FR_SS_tauTau_EE->SetMarkerSize(0);
   FR_SS_tauTau_EE->SetTitle("endcap corrected");
   mg_tauTau->Add(FR_SS_tauTau_EB_unc);
   FR_SS_tauTau_EB_unc->SetLineColor(kBlue);
   FR_SS_tauTau_EB_unc->SetLineStyle(1);
   FR_SS_tauTau_EB_unc->SetMarkerSize(0);
   FR_SS_tauTau_EB_unc->SetTitle("barel uncorrected");
   mg_tauTau->Add(FR_SS_tauTau_EE_unc);
   FR_SS_tauTau_EE_unc->SetLineColor(kRed);
   FR_SS_tauTau_EE_unc->SetLineStyle(1);
   FR_SS_tauTau_EE_unc->SetMarkerSize(0);
   FR_SS_tauTau_EE_unc->SetTitle("endcap uncorrected");  

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
   leg_ele = CreateLegend_FR("left",FR_SS_electron_EB_unc,FR_SS_electron_EB,FR_SS_electron_EE_unc,FR_SS_electron_EE);
   leg_ele->Draw();
   system("mkdir -p Plots");
   SavePlots(c_ele, "Plots/FR_SS_electrons");
   
   c_mu->cd();
   lumi->set_lumi(c_mu, _lumi, 0);
   mg_muons->Draw("AP");
	mg_muons->GetXaxis()->SetTitle("p_{T} [GeV]");
	mg_muons->GetYaxis()->SetTitle("Fake Rate");
	mg_muons->SetTitle("Muon fake rate");
   mg_muons->SetMaximum(0.35);
   leg_mu = CreateLegend_FR("left",FR_SS_muon_EB_unc,FR_SS_muon_EB,FR_SS_muon_EE_unc,FR_SS_muon_EE);
   leg_mu->Draw();
   SavePlots(c_mu, "Plots/FR_SS_muons");
   
   c_tauE->cd();
   lumi->set_lumi(c_tauE, _lumi, 0);
   mg_tauE->Draw("AP");
   mg_tauE->GetXaxis()->SetTitle("p_{T} [GeV]");
   mg_tauE->GetYaxis()->SetTitle("Fake Rate");
   mg_tauE->SetTitle("#tau fake rate, #tau_{e}#tau_{h} channel");
   mg_tauE->SetMaximum(0.05);
   leg_tauE = CreateLegend_FR("left",FR_SS_tauE_EB_unc,FR_SS_tauE_EB,FR_SS_tauE_EE_unc,FR_SS_tauE_EE);
   leg_tauE->Draw();
   SavePlots(c_tauE, "Plots/FR_SS_tauE");

   c_tauMu->cd();
   lumi->set_lumi(c_tauMu, _lumi, 0);
   mg_tauMu->Draw("AP");
   mg_tauMu->GetXaxis()->SetTitle("p_{T} [GeV]");
   mg_tauMu->GetYaxis()->SetTitle("Fake Rate");
   mg_tauMu->SetTitle("#tau fake rate, #tau_{e}#tau_{h} channel");
   mg_tauMu->SetMaximum(0.05);
   leg_tauMu = CreateLegend_FR("left",FR_SS_tauMu_EB_unc,FR_SS_tauMu_EB,FR_SS_tauMu_EE_unc,FR_SS_tauMu_EE);
   leg_tauMu->Draw();
   SavePlots(c_tauMu, "Plots/FR_SS_tauMu");

   c_tauTau->cd();
   lumi->set_lumi(c_tauTau, _lumi, 0);
   mg_tauTau->Draw("AP");
   mg_tauTau->GetXaxis()->SetTitle("p_{T} [GeV]");
   mg_tauTau->GetYaxis()->SetTitle("Fake Rate");
   mg_tauTau->SetTitle("#tau fake rate, #tau_{e}#tau_{h} channel");
   mg_tauTau->SetMaximum(0.05);
   leg_tauTau = CreateLegend_FR("left",FR_SS_tauTau_EB_unc,FR_SS_tauTau_EB,FR_SS_tauTau_EE_unc,FR_SS_tauTau_EE);
   leg_tauTau->Draw();
   SavePlots(c_tauTau, "Plots/FR_SS_tauTau");  

}
//===============================================================

//========================================================================================================
void SSmethod::PlotDataMC( TString variable_name, TString folder )
{
   TCanvas *c;
   c = new TCanvas("ZLLss", variable_name, 600, 600);
   
   if ( GetVarLogX( variable_name) ) c->SetLogx();
   if ( GetVarLogY( variable_name) ) c->SetLogy();
   
   for( int i_fs = 0; i_fs <= Settings::fs4l ; i_fs++ )
   {
      //for ( int i_cat = 0; i_cat <= Settings::inclusive_stxs; i_cat++ )
      //{
         histos_1D[Settings::regZLL][Settings::WZ][i_fs]   ->SetFillColor(kMagenta-7);
         histos_1D[Settings::regZLL][Settings::qqZZ][i_fs] ->SetFillColor(kCyan+1);
         histos_1D[Settings::regZLL][Settings::DY][i_fs]   ->SetFillColor(kGreen+2);
         histos_1D[Settings::regZLL][Settings::ttbar][i_fs]->SetFillColor(kBlue-4);
         
         histos_1D[Settings::regZLL][Settings::WZ][i_fs]   ->SetLineColor(kMagenta-7);
         histos_1D[Settings::regZLL][Settings::qqZZ][i_fs] ->SetLineColor(kCyan+1);
         histos_1D[Settings::regZLL][Settings::DY][i_fs]   ->SetLineColor(kGreen+2);
         histos_1D[Settings::regZLL][Settings::ttbar][i_fs]->SetLineColor(kBlue-4);
         
         histos_1D[Settings::regZLL][Settings::Data][i_fs]->SetMarkerSize(0.8);
         histos_1D[Settings::regZLL][Settings::Data][i_fs]->SetMarkerStyle(20);
         histos_1D[Settings::regZLL][Settings::Data][i_fs]->SetBinErrorOption(TH1::kPoisson);
         histos_1D[Settings::regZLL][Settings::Data][i_fs]->SetLineColor(kBlack);
         
         THStack *stack = new THStack( "stack", "stack" );
			stack->Add(histos_1D[Settings::regZLL][Settings::qqZZ][i_fs]);
			stack->Add(histos_1D[Settings::regZLL][Settings::WZ][i_fs]);
			stack->Add(histos_1D[Settings::regZLL][Settings::ttbar][i_fs]);
         stack->Add(histos_1D[Settings::regZLL][Settings::DY][i_fs]);
			
         
         stack->Draw("HIST");
         
         float data_max = histos_1D[Settings::regZLL][Settings::Data][i_fs]->GetBinContent(histos_1D[Settings::regZLL][Settings::Data][i_fs]->GetMaximumBin());
         float data_max_error = histos_1D[Settings::regZLL][Settings::Data][i_fs]->GetBinErrorUp(histos_1D[Settings::regZLL][Settings::Data][i_fs]->GetMaximumBin());
         
         stack->SetMinimum(1e-5);
         stack->SetMaximum((data_max + data_max_error)*1.1);
         
         TString _fs_label;
	 if ( i_fs == Settings::fs4l) _fs_label = "m_{4#font[12]{l}} (GeV)";
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
         stack->GetYaxis()->SetTitle(histos_1D[Settings::regZLL][Settings::Data][i_fs]->GetYaxis()->GetTitle());
         stack->GetYaxis()->SetTitleSize(0.04);
         stack->GetYaxis()->SetLabelSize(0.04);
         
         stack->GetXaxis()->SetTitleOffset(1.2);
         stack->GetYaxis()->SetTitleOffset(1.25);
         
         histos_1D[Settings::regZLL][Settings::Data][i_fs]->Draw("SAME p E1 X0");
         
         TLegend *legend;
         legend  = CreateLegend_ZLL("right",histos_1D[Settings::regZLL][Settings::Data][i_fs],histos_1D[Settings::regZLL][Settings::WZ][i_fs],histos_1D[Settings::regZLL][Settings::qqZZ][i_fs],histos_1D[Settings::regZLL][Settings::DY][i_fs],histos_1D[Settings::regZLL][Settings::ttbar][i_fs]);
         legend->Draw();
         
         // Draw lumi
         CMS_lumi *lumi = new CMS_lumi;
         lumi->set_lumi(c, _lumi, 0);
         
         TString _out_file_name;
         _out_file_name = folder + "/" + variable_name + "_SS_" + _s_final_state.at(i_fs);// + "_" + _s_category_stxs.at(i_cat);
         SavePlots(c, _out_file_name);
         
      //}
      
   }
}
//========================================================================================================


//========================================================================================================
void SSmethod::PlotZX( TString variable_name, TString folder )
{
   TCanvas *c;
   c = new TCanvas("ZLLss", variable_name, 600, 600);
   
   if ( GetVarLogX( variable_name) ) c->SetLogx();
   if ( GetVarLogY( variable_name) ) c->SetLogy();
   
   for( int i_fs = 0; i_fs <= Settings::fs4l ; i_fs++ )
   {
      //for ( int i_cat = 0; i_cat <= Settings::inclusive_stxs; i_cat++ )
      //{
         histos_ZX[Settings::regZLL][Settings::Data][i_fs]->SetFillColor(kGreen+2);
         histos_ZX[Settings::regZLL][Settings::Data][i_fs]->SetLineColor(kGreen+2);
         
         histos_ZX[Settings::regZLL][Settings::Data][i_fs]->Draw("HIST");
         
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
         histos_ZX[Settings::regZLL][Settings::Data][i_fs]->GetXaxis()->SetTitle(_fs_label);
         histos_ZX[Settings::regZLL][Settings::Data][i_fs]->GetXaxis()->SetTitleSize(0.04);
         histos_ZX[Settings::regZLL][Settings::Data][i_fs]->GetXaxis()->SetLabelSize(0.04);
         histos_ZX[Settings::regZLL][Settings::Data][i_fs]->GetYaxis()->SetTitle(histos_ZX[Settings::regZLL][Settings::Data][i_fs]->GetYaxis()->GetTitle());
         histos_ZX[Settings::regZLL][Settings::Data][i_fs]->GetYaxis()->SetTitleSize(0.04);
         histos_ZX[Settings::regZLL][Settings::Data][i_fs]->GetYaxis()->SetLabelSize(0.04);
         
         histos_ZX[Settings::regZLL][Settings::Data][i_fs]->GetXaxis()->SetTitleOffset(1.2);
         histos_ZX[Settings::regZLL][Settings::Data][i_fs]->GetYaxis()->SetTitleOffset(1.25);
         
         // Draw lumi
         CMS_lumi *lumi = new CMS_lumi;
         lumi->set_lumi(c, _lumi, 0);
         
         TString _out_file_name;
         _out_file_name = folder + "/" + variable_name + "_ZX_SS_" + _s_final_state.at(i_fs);// + "_" + _s_category_stxs.at(i_cat);
         SavePlots(c, _out_file_name);
         
      //}
      
   }
}
//========================================================================================================


//========================================================================================================
void SSmethod::FitZX( TString variable_name, TString folder )
{
   TCanvas *c;
   TF1  *fit_function;
   CMS_lumi *lumi = new CMS_lumi;
	
   c = new TCanvas("Fits_ZLLss", variable_name, 600, 600);
	
   if ( GetVarLogX( variable_name) ) c->SetLogx();
   if ( GetVarLogY( variable_name) ) c->SetLogy();
	
   for( int i_fs = 0; i_fs <= Settings::fs4l ; i_fs++ )
   {
      //for ( int i_cat = 0; i_cat <= Settings::inclusive_stxs; i_cat++ )
      //{
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

         histos_ZX[Settings::regZLL][Settings::Data][i_fs]->GetXaxis()->SetTitle(_fs_label);
         histos_ZX[Settings::regZLL][Settings::Data][i_fs]->GetXaxis()->SetTitleSize(0.04);
         histos_ZX[Settings::regZLL][Settings::Data][i_fs]->GetXaxis()->SetLabelSize(0.04);
         histos_ZX[Settings::regZLL][Settings::Data][i_fs]->GetYaxis()->SetTitle(histos_ZX[Settings::regZLL][Settings::Data][i_fs]->GetYaxis()->GetTitle());
         histos_ZX[Settings::regZLL][Settings::Data][i_fs]->GetYaxis()->SetTitleSize(0.04);
         histos_ZX[Settings::regZLL][Settings::Data][i_fs]->GetYaxis()->SetLabelSize(0.04);
			
         histos_ZX[Settings::regZLL][Settings::Data][i_fs]->GetXaxis()->SetTitleOffset(1.2);
         histos_ZX[Settings::regZLL][Settings::Data][i_fs]->GetYaxis()->SetTitleOffset(1.25);
			
			gStyle->SetOptFit();
			gStyle->SetStatY(0.9);
			gStyle->SetStatX(0.95);
			gStyle->SetStatW(0.2);
			gStyle->SetStatH(0.1);
			
			fit_function = new TF1("fit_function","[0]*TMath::Landau(x, [1], [2])",70,800);
			fit_function->SetParNames("Constant","MPV","#sigma");
			fit_function->SetParameter(0,1.);
			fit_function->SetParameter(1,100.);
			fit_function->SetParameter(2,10.);

			histos_ZX[Settings::regZLL][Settings::Data][i_fs]->Fit("fit_function");
			histos_ZX[Settings::regZLL][Settings::Data][i_fs]->Draw("");
			
         // Draw lumi
         lumi->set_lumi(c, _lumi, 0);
			
         TString _out_file_name;
         _out_file_name = folder + "/" + variable_name + "_ZX_SS_fit_" + _s_final_state.at(i_fs);// + "_" + _s_category_stxs.at(i_cat);
         SavePlots(c, _out_file_name);
      //}
	}
	gStyle->SetOptFit(0);
}
//========================================================================================================



//===============================================================
void SSmethod::RemoveNegativeBins1D(TH1F *h)
{
   for (int i_bin_x = 1; i_bin_x <= h->GetXaxis()->GetNbins(); i_bin_x++)
   {
      if( h->GetBinContent(i_bin_x) < 0.) h->SetBinContent(i_bin_x, 0);
   }
   
}
//===============================================================

//===============================================================
void SSmethod::RemoveNegativeBins2D(TH2F *h)
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
void SSmethod::Set_pT_binning(int size, float *bins)
{
   _n_pT_bins = size;

   for (int i = 0; i < size; i++)
   {
      _pT_bins[i] = bins[i];
   }
}
//===============================================================

//===============================================================
void SSmethod::SetLumi(float lumi)
{
   _lumi = lumi;
}
//===============================================================


//==========================================================
int SSmethod::find_current_process( TString input_file_name )
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
int SSmethod::FindFinalState()
{
   int final_state = -999;

   if ( Z1Flav == -121 )
   {
      if ( abs(Z2Flav) == 169 ) final_state = Settings::fs2e2mu;
      else if (abs(Z2Flav) == 165) final_state = Settings::fs2eetau;
      else if (abs(Z2Flav) == 195) final_state = Settings::fs2emutau;
      else if (abs(Z2Flav) == 225) final_state = Settings::fs2etautau;
      else
         final_state = -2;//cerr << "[ERROR] in event " << RunNumber << ":" << LumiNumber << ":" << EventNumber << ", Z2Flav = " << Z2Flav << endl;
      }
   else if ( Z1Flav == -169 )
   {
      if ( abs(Z2Flav) == 169 ) final_state = Settings::fs4mu;
      else if (abs(Z2Flav) == 165) final_state = Settings::fs2muetau;
      else if (abs(Z2Flav) == 195) final_state = Settings::fs2mumutau;
      else if (abs(Z2Flav) == 225) final_state = Settings::fs2mutautau;
      else
         final_state = -2;//cerr << "[ERROR] in event " << RunNumber << ":" << LumiNumber << ":" << EventNumber << ", Z2Flav = " << Z2Flav << endl;
   }
   else
   {
      final_state = -1;//cerr << "[ERROR] in event " << RunNumber << ":" << LumiNumber << ":" << EventNumber << ", Z1Flav = " << Z1Flav << endl;
   }
   
   return final_state;
}
//=============================

//=========================================
int SSmethod::Find_Ele_pT_bin( Float_t pT )
{
	int bin = 0;
	
	for ( int i = 2; i <= _n_pT_bins ; i++)
	{
		if ( (pT > _pT_bins[i-1]) && (pT < _pT_bins[i]) ) bin = i - 2;
	}
	if ( pT > 80. ) bin = _n_pT_bins - 2;

	//cout << "PT = " << pT << " bin = " << bin << endl;
	return bin;
}
//=========================================

//=========================================
int SSmethod::Find_Ele_eta_bin( Float_t eta )
{
	int bin = 0;
	
	if (abs(eta) < 1.479) bin = Settings::EB;
	else bin = Settings::EE;
	
	//cout << "eta = " << eta << " bin = " << bin << endl;
	return bin;
}
//=========================================


//=================================
float SSmethod::calculate_K_factor(TString input_file_name)
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
int SSmethod::GENMatch(int ileg,TString input_file_name)
{
   if (abs(LepLepId->at(ileg))!=15) {
      cout<<"Only applicable to pdgId = 15, not "<<ileg<<endl;
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

float SSmethod::calculate_TauIDSF_OneLeg(short DM, float pt, float eta, int genmatch, int flav)
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

float SSmethod::calculate_TauIDSF(TString input_file_name)
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
bool SSmethod::GetVarLogX ( TString variable_name )
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
bool SSmethod::GetVarLogY ( TString variable_name )
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
TLegend* SSmethod::CreateLegend_FR( string position, TGraphErrors *EB_unc, TGraphErrors *EB_cor,TGraphErrors *EE_unc,TGraphErrors *EE_cor )
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
TLegend* SSmethod::CreateLegend_ZLL( string position, TH1F *data, TH1F *WZ,TH1F *qqZZ,TH1F *DY,TH1F *ttbar )
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

//=======================================
void SSmethod::SavePlots( TCanvas *c, TString name)
{
   c->SaveAs(name + ".pdf");
   //c->SaveAs(name + ".root");
   //c->SaveAs(name + ".eps");
   c->SaveAs(name + ".png");
   //gSystem->Exec("convert -density 300 -quality 100 " + name + ".eps " + name + ".png");
}
//=======================================
