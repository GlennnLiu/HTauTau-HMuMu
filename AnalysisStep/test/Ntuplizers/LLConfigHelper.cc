#include "LLConfigHelper.h"
#include "LLConfigHelper.h"
#include <HTauTauHMuMu/AnalysisStep/interface/bitops.h>

using namespace std;
using namespace edm;

LLConfigHelper::LLConfigHelper(const ParameterSet& pset) :
  PD(pset.getParameter<std::string>("PD")),
  isMC_(pset.getUntrackedParameter<bool>("isMC")),
  theSetup(pset.getParameter<int>("setup")),
  theSampleType(pset.getParameter<int>("sampleType")),
  skimPaths(pset.getParameter<std::vector<std::string> >("skimPaths")),
  MCFilter(pset.getParameter<std::string>("MCFilterPath")),
  anyTrigger(true) // take the OR of all trigger paths
  
{
  string channel = pset.getUntrackedParameter<string>("channel");
  theChannel = finalState(channel);
  
  // Check for inconsistent configurations
  if ( ( theSampleType!=2011 && theSampleType!=2012 && theSampleType!=2015 && theSampleType!=2016 && theSampleType!=2017 && theSampleType!=2018) ||
       ( theSetup!=2011 && theSetup!=2012 && theSetup!=2015 && theSetup!=2016 && theSetup!=2017 && theSetup!=2018) ||
       ( theSampleType!=theSetup ) // No sample rescaling supported as of now.
       // We may add exception for MC only when needed.
       ) {
    cout << "ERROR: LLConfigHelper: inconsistent setup: sampleType=" << theSampleType << ", setup=" << theSetup << ", isMC=" <<isMC_ << endl;
    abort();
  }
  
  
  if ((isMC_&&PD!="") || (!isMC_ && (PD!="DoubleMuon" && PD!="SingleMuon" && PD!="SingleElectron" && PD!="MuonEG" && PD!="Tau"))) {
    cout << "ERROR: LLConfigHelper: isMC: " << isMC_ << " PD: " << PD << endl;
    abort();
  }    

  if (!isMC_&&MCFilter!="") {
    cout << "ERROR: LLConfigHelper: MCFilter= " << MCFilter << " when isMC=0" << endl;
    abort();
  }
  
}

bool LLConfigHelper::passMCFilter(const edm::Event & event, edm::Handle<edm::TriggerResults> & trigRes){
  if (MCFilter=="") return true;
  return passFilter(event, trigRes, MCFilter);
}

bool 
LLConfigHelper::passSkim(const edm::Event & event, edm::Handle<edm::TriggerResults> & trigRes, short& trigworld){
  bool evtPassSkim = false;
  if (skimPaths.size()==0) {
    evtPassSkim=true;
  } else {
    for (vector<string>::const_iterator name = skimPaths.begin(); name!= skimPaths.end(); ++name) {
      if (passFilter(event, trigRes, *name)) {
	evtPassSkim = true; 
	break;
      }
    }
  }
  if (evtPassSkim) set_bit_16(trigworld,15);
  return evtPassSkim;
}

bool 
LLConfigHelper::passTrigger(const edm::Event & event, edm::Handle<edm::TriggerResults> & trigRes, short& trigworld){

  bool passSingleMu  = passFilter(event, trigRes, "triggerSingleMu");
  bool passSingleEle  = passFilter(event, trigRes, "triggerSingleEle");
  bool passDiMu  = passFilter(event, trigRes, "triggerDiMu");
  bool passDiTau = passFilter(event, trigRes, "triggerDiTau");
  bool passMuTau = passFilter(event, trigRes, "triggerMuTau");
  bool passEleTau = passFilter(event, trigRes, "triggerEleTau");
  bool passMuEle = passFilter(event, trigRes, "triggerMuEle");
  
  bool evtPassTrigger = false;

  // Check all triggers together if anyTrigger is specified (or for CRs)
  if (anyTrigger || theChannel=='SR') {
    if ((PD=="" && (passSingleMu || passSingleEle || passDiMu || passDiTau || passMuTau || passEleTau || passMuEle )) ||
	  ((PD=="DoubleMuon") && passDiMu) ||
    ((PD=="MuonEG") && passMuEle && !passDiMu) ||
    ((PD=="SingleMuon") && passSingleMu && !passDiMu && !passMuEle) ||
	  ((PD=="SingleElectron") && passSingleEle && !passDiMu && !passMuEle && !passSingleMu) ||
    ((PD=="Tau") && (passDiTau || passMuTau || passEleTau) && !passDiMu && !passMuEle && !passSingleMu && !passSingleEle)
	  ) {
      evtPassTrigger = true;
      } 
  }
  
  else {
    cout << "[ERROR]: LLConfigHelper: unexpected config " << theChannel << endl;
    abort();
  }
  
  if (evtPassTrigger) set_bit_16(trigworld,0);
  if (passDiMu) set_bit_16(trigworld,1);
  if (passMuEle) set_bit_16(trigworld,2);
  if (passSingleMu) set_bit_16(trigworld,3);
  if (passSingleEle) set_bit_16(trigworld,4);
  if (passDiTau) set_bit_16(trigworld,5);
  if (passMuTau) set_bit_16(trigworld,6);
  if (passEleTau) set_bit_16(trigworld,7);
  

  return evtPassTrigger;
}

// bool
// LLConfigHelper::passMETTrigger(const edm::Event & event, edm::Handle<edm::TriggerResults> & trigRes){

//   bool passMETFilter  = true;
	
// //  if (theSetup >= 2017) {
// //    passMETFilter = passFilter(event, trigRes, "triggerMETFilters");
// //  }
//   return passMETFilter;
// }


bool
LLConfigHelper::passFilter(const edm::Event & event, edm::Handle<edm::TriggerResults> & trigRes, const string& filterPath) {

  //if (event.id()==cachedEvtId) return;
  //cachedEvtId = event.id();

  // Initialize trigger results table
  if (trigRes.isValid()) {
    triggerResults = trigRes;
    triggerNames = &(event.triggerNames(*triggerResults));
  } else {
    cout << "ERROR: failed to get TriggerResults" << endl;
  }

  //  for (unsigned i=0; i<triggerNames->size(); i++) cout << triggerNames->triggerName(i) << endl;
  unsigned i =  triggerNames->triggerIndex(filterPath);
  
  if (i== triggerNames->size()){
    cout << "ERROR: LLConfigHelper::isTriggerBit: path does not exist! " << filterPath << endl;
    abort();
  }
  //  cout << " Trigger result for " << filterPath << " : accept=" << triggerResults->accept(i) << endl;
  return triggerResults->accept(i);

}




