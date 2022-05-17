#ifndef utils_h
#define utils_h

#include <TLorentzVector.h>
#include <DataFormats/Math/interface/LorentzVector.h>

namespace htautauhmumu{
  TLorentzVector tlv(const math::XYZTLorentzVector& v);
}
int SetupToSqrts(int setup);



#endif

