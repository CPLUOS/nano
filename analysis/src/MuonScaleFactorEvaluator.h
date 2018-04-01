#include <vector>
#include <utility>
#include <algorithm>
#include <cassert>
#include "TLorentzVector.h"

class MuonScaleFactorEvaluator
{
public:
  MuonScaleFactorEvaluator(){
    const unsigned int n = (pt_bins.size()-1)*(eta_bins.size()-1);
    // FIXME : check that these bins are monolothically increasing
    assert(values.size() == n);
    assert(errors.size() == n);

    // For cache
    width = pt_bins.size()-1;
  };
  double operator()(const double x, const double y, const double shift = 0){
    auto xbin = std::lower_bound(pt_bins.begin(), pt_bins.end(), x);
    if ( xbin == pt_bins.end() || xbin+1 == pt_bins.end() ) return 1;
    auto ybin = std::lower_bound(eta_bins.begin(), eta_bins.end(), y);
    if ( ybin == eta_bins.end() || ybin+1 == eta_bins.end() ) return 1;

    const int column = xbin-pt_bins.begin();
    const int row = ybin-eta_bins.begin();

    const int bin = row*width+column;
    const double value = values.at(bin);
    const double error = errors.at(bin);

    return std::max(0.0, value+shift*error);
  };
  double getScaleFactor(const TParticle& p, const int pid, const double shift = 0){
    const int aid = abs(p.GetPdgCode());
    if ( aid == pid ) {
      const double x = p.Pt(), y = p.Eta();
      
      auto xbin = std::lower_bound(pt_bins.begin(), pt_bins.end(), x);
      if ( xbin == pt_bins.end() || xbin+1 == pt_bins.end() ) return 1;
      auto ybin = std::lower_bound(eta_bins.begin(), eta_bins.end(), y);
      if ( ybin == eta_bins.end() || ybin+1 == eta_bins.end() ) return 1;

      const int column = xbin-pt_bins.begin();
      const int row = ybin-eta_bins.begin();

      const int bin = row*width+column;
      const double value = values.at(bin);
      const double error = errors.at(bin);

      return std::max(0.0, value+shift*error);
    }
    return 1;
  };

private:
  static const std::vector<double> pt_bins;
  static const std::vector<double> eta_bins;
  static const std::vector<double> values;
  static const std::vector<double> errors;

  int width;
};

//For Muons
const std::vector<double> MuonScaleFactorEvaluator::pt_bins = {
  10.000000,  20.000000,  25.000000,  30.000000,  40.000000,  50.000000,  60.000000,  120.000000,  10000.000000,
};

const std::vector<double> MuonScaleFactorEvaluator::eta_bins = {
  -2.400000,  -2.100000,  -1.600000,  -1.200000,  -0.900000,  -0.600000,  -0.300000,  -0.200000,
  0.000000,  0.200000,  0.300000,  0.600000,  0.900000,  1.200000,  1.600000,  2.100000,  2.400000,
};

const std::vector<double> MuonScaleFactorEvaluator::values = {
  0.992278,  0.992278,  0.992278,  0.992278,  0.992278,  0.992278,  0.992278,  0.992278,  0.995300,
  0.955047,  0.961484,  0.961262,  0.966372,  0.961815,  0.960081,  0.995300,  0.996584,  0.974587,
  0.978759,  0.982496,  0.985181,  0.980446,  0.984624,  0.996584,  0.997201,  0.975192,  0.979367,
  0.983106,  0.985792,  0.981053,  0.985234,  0.997201,  0.997698,  0.965529,  0.968893,  0.970796,
  0.970370,  0.971478,  0.971802,  0.997698,  0.997773,  0.969139,  0.974249,  0.976080,  0.979071,
  0.976589,  0.988743,  0.997773,  0.997067,  0.968451,  0.973557,  0.975387,  0.978375,  0.975896,
  0.988043,  0.997067,  0.997287,  0.968666,  0.973772,  0.975603,  0.978592,  0.976112,  0.988261,
  0.997287,  0.997287,  0.968666,  0.973772,  0.975603,  0.978592,  0.976112,  0.988261,  0.997287,
  0.997948,  0.969309,  0.974419,  0.976251,  0.979242,  0.976760,  0.988916,  0.997948,  0.998722,
  0.970063,  0.975178,  0.977010,  0.980005,  0.977520,  0.989683,  0.998722,  0.998437,  0.969787,
  0.974899,  0.976732,  0.979725,  0.977241,  0.989401,  0.998437,  0.997708,  0.965550,  0.968913,
  0.970813,  0.970386,  0.971492,  0.971818,  0.997708,  0.995418,  0.973447,  0.977615,  0.981347,
  0.984029,  0.979298,  0.983472,  0.995418,  0.995424,  0.973443,  0.977610,  0.981345,  0.984028,
  0.979303,  0.983477,  0.995424,  0.988592,  0.948495,  0.954918,  0.954743,  0.959830,  0.955324,
  0.953610,  0.988592,
};

const std::vector<double> MuonScaleFactorEvaluator::errors = {
  0.000661,  0.000661,  0.000661,  0.000661,  0.000661,  0.000661,  0.000661,  0.000661,  0.000177,
  0.010855,  0.026217,  0.010525,  0.011035,  0.010665,  0.011252,  0.000177,  0.000153,  0.010676,
  0.010598,  0.014068,  0.010595,  0.010591,  0.010707,  0.000153,  0.000173,  0.010676,  0.010598,
  0.014068,  0.010595,  0.010591,  0.010707,  0.000173,  0.000082,  0.012978,  0.010674,  0.010576,
  0.010503,  0.010558,  0.010714,  0.000082,  0.000079,  0.010801,  0.023378,  0.010543,  0.010555,
  0.010557,  0.010697,  0.000079,  0.000145,  0.010802,  0.023378,  0.010544,  0.010555,  0.010558,
  0.010697,  0.000145,  0.000063,  0.010801,  0.023378,  0.010543,  0.010555,  0.010557,  0.010697,
  0.000063,  0.000063,  0.010801,  0.023378,  0.010543,  0.010555,  0.010557,  0.010697,  0.000063,
  0.000145,  0.010802,  0.023378,  0.010544,  0.010555,  0.010558,  0.010697,  0.000145,  0.000077,
  0.010801,  0.023378,  0.010543,  0.010555,  0.010557,  0.010697,  0.000077,  0.000084,  0.010801,
  0.023378,  0.010543,  0.010555,  0.010557,  0.010697,  0.000084,  0.000170,  0.012979,  0.010675,
  0.010577,  0.010504,  0.010559,  0.010715,  0.000170,  0.000160,  0.010676,  0.010598,  0.014068,
  0.010595,  0.010591,  0.010707,  0.000160,  0.000165,  0.010676,  0.010598,  0.014068,  0.010595,
  0.010591,  0.010707,  0.000165,  0.000659,  0.010874,  0.026225,  0.010544,  0.011053,  0.010684,
  0.011270,  0.000659,
};