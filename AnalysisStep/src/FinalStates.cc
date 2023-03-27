#include <HTauTauHMuMu/AnalysisStep/interface/FinalStates.h>

namespace {
  const unsigned nFS = 29;
}


std::string finalState(int iFS) {
  if (iFS<0 || iFS>=(int)nFS) return "None";
  const std::string finalStates[nFS] = {
      "mumu","etau","mutau","tautau","emu","ee"
  };
  return finalStates[iFS];			     	
}

Channel finalState(std::string sFS) {
  for (unsigned i=0; i<nFS; ++i) {
    if (sFS==finalState(i)) return (Channel) i;
  }
  return NONE;
}
