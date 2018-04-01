#include "BTagCalibrationStandalone.cc"

#include <string>
#include <vector>
#include <map>

#include "TFile.h"
#include "TH1D.h"
#include "TParticle.h"

using namespace std;

class BTagWeightEvaluator
{
public:
  BTagWeightEvaluator() {};
  void initCSVWeight();

  // For per-jet SF evaluation - useful for CSV weight
  double getSF(const TParticle& jet, const float CSVV2, const int hadronFlavour, const int unc);

  enum CSVUNC {
    CENTRAL, JES_UP, JES_DN,
    LF_UP, LF_DN, HF_UP, HF_DN,
    HFSTAT1_UP, HFSTAT1_DN, HFSTAT2_UP, HFSTAT2_DN,
    LFSTAT1_UP, LFSTAT1_DN, LFSTAT2_UP, LFSTAT2_DN,
    CFERR1_UP, CFERR1_DN, CFERR2_UP, CFERR2_DN
  };

private:

  std::vector<std::string> uncNames_;
  std::map<int, BTagCalibrationReader> readers_;
};

void BTagWeightEvaluator::initCSVWeight()
{
  string csvFileName = "btagSF_CSVv2_ichep2016.csv";
  std::string env = std::getenv("CMSSW_BASE");
  std::string csvFile = env+"/src/nano/analysis/data/btagSF/"+csvFileName;
  std::cout << csvFile << std::endl;
  BTagCalibration calib("csvv2", csvFile);
  uncNames_ = {
    "central", "up_jes", "down_jes",
    "up_lf", "down_lf", "up_hf", "down_hf",
    "up_hfstats1", "down_hfstats1", "up_hfstats2", "down_hfstats2",
    "up_lfstats1", "down_lfstats1", "up_lfstats2", "down_lfstats2",
    "up_cferr1", "down_cferr1", "up_cferr2", "down_cferr2"
  };
  for ( unsigned int i=0; i<uncNames_.size(); ++i ) {
    std::cout << "test : " << i << std::endl;
    readers_[i] = BTagCalibrationReader(&calib, BTagEntry::OP_RESHAPING, "iterativefit", uncNames_[i]);
  }
}

double BTagWeightEvaluator::getSF(const TParticle& jet, const float CSVV2, const int hadronFlavour, const int unc)
{
  const double pt = std::min(jet.Pt(), 999.);
  const double eta = jet.Eta();
  const double aeta = std::abs(eta);
  if ( pt <= 20 or aeta >= 2.4 ) return 1.0;

  double discr = CSVV2;
  const int flav = std::abs(hadronFlavour);

  if      ( discr < -1.0 ) discr = -0.05;
  else if ( discr >  1.0 ) discr = 1.0;

  // Special care for the flavour dependent SFs
  int uncKey = unc;
  if ( flav == 5 ) {
    if ( unc != LF_UP and unc != LF_DN and
        unc != HFSTAT1_UP and unc != HFSTAT1_DN and
        unc != HFSTAT2_UP and unc != HFSTAT2_DN ) uncKey = CENTRAL;
  }
  else if ( flav == 4 ) {
    if ( unc != CFERR1_UP and unc != CFERR1_DN and
        unc != CFERR2_UP and unc != CFERR2_DN ) uncKey = CENTRAL;
  }
  else {
    if ( unc != HF_UP and unc != HF_DN and
        unc != LFSTAT1_UP and unc != LFSTAT1_DN and
        unc != LFSTAT2_UP and unc != LFSTAT2_DN ) uncKey = CENTRAL;
  }

  auto readerItr = readers_.find(uncKey);
  if ( readerItr == readers_.end() ) return 1.0;

  BTagEntry::JetFlavor jf = BTagEntry::FLAV_UDSG;
  if      ( flav == 5 ) jf = BTagEntry::FLAV_B;
  else if ( flav == 4 ) jf = BTagEntry::FLAV_C;

  return readerItr->second.eval(jf, aeta, pt, discr);
  
  return 1.;
};
