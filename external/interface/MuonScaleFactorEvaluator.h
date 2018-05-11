#ifndef MUONSCALE_H
#define MUONSCALE_H

#include <vector>
#include <utility>
#include <algorithm>
#include <cassert>
#include "TLorentzVector.h"
#include "TParticle.h"

class MuonScaleFactorEvaluator
{
public:
  MuonScaleFactorEvaluator();
  double operator()(const double x, const double y, const double shift = 0);
  double getScaleFactor(const TParticle& p, const int pid, const double shift = 0);

private:
  static const std::vector<double> pt_bins;
  static const std::vector<double> eta_bins;
  static const std::vector<double> values;
  static const std::vector<double> errors;

  int width;
};

#endif
