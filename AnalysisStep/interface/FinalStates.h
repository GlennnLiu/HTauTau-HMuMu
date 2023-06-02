#ifndef FinalStates_h
#define FinalStates_h
// Some standard labels for final states.

#include <string>

enum Channel {mumu=0, 
			etau=1, 
			mutau=2, 
			tautau=3,
			emu=4,
			ee=5,
			SR=6,
			CRTT=7,
			CRQCD=8,
			CRWJ=9,
			CRQCDvSR=10,
			CRQCDvSROS=11,
			CRQCDvSRSS=12,
			CRWJvSR=13,
			CRAPPOS=14,
			CRAPPSS=15,
			NONE = 99, BUGGY=666};

//Return string corresponding to integer code 
std::string finalState(int iFS);

//Return integer code corresponding to name
Channel finalState(std::string sFS);

#endif
