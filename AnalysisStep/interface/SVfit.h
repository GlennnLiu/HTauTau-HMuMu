#ifndef SVfit_h
#define SVfit_h

// Standard libraries
#include <vector>
#include <string>
#include <cmath>

// ROOT libraries
#include <TLorentzVector.h>
#include "TMatrixD.h"

// ClassicSVfit libraries
#include <TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h>
#include <TauAnalysis/ClassicSVfit/interface/ClassicSVfitIntegrand.h>
#include <TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h>
#include <TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h>

using namespace classic_svFit;

// SVfit class
class SVfit {

  public:
    SVfit (int verbosity, TLorentzVector tau1, TLorentzVector tau2, TLorentzVector met, TMatrixD met_cov, int pairType, int DM1, int DM2);
    ~SVfit ();

    std::vector<double> FitAndGetResult();

  private:
    int verbosity_;
    std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons_;
    double METx_;
    double METy_;
    TMatrixD covMET_;
    double kappa_;

};

#endif // SVfit_h
