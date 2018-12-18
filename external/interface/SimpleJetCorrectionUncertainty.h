// 
// This code is from copying 
//   https://github.com/cms-sw/cmssw/blob/CMSSW_9_4_X/CondFormats/JetMETObjects/interface/SimpleJetCorrectionUncertainty.h
//   in 4d7b6ef commit
// Copied by Byeonghak Ko
//


#ifndef SimpleJetCorrectionUncertainty_h
#define SimpleJetCorrectionUncertainty_h

//#include "CondFormats/Serialization/interface/Serializable.h" // MODIFIED


// MODIFIED
#define COND_SERIALIZABLE
#define COND_TRANSIENT


#include <string>
#include <vector>
class JetCorrectorParameters;

class SimpleJetCorrectionUncertainty 
{
 public:
  SimpleJetCorrectionUncertainty();
  SimpleJetCorrectionUncertainty(const std::string& fDataFile);
  SimpleJetCorrectionUncertainty(const JetCorrectorParameters& fParameters);
  ~SimpleJetCorrectionUncertainty();
  const JetCorrectorParameters& parameters() const {return *mParameters;}
  float uncertainty(const std::vector<float>& fX, float fY, bool fDirection) const;

 private:
  SimpleJetCorrectionUncertainty(const SimpleJetCorrectionUncertainty&) = delete;
  SimpleJetCorrectionUncertainty& operator= (const SimpleJetCorrectionUncertainty&) = delete;
  int findBin(const std::vector<float>& v, float x) const;
  float uncertaintyBin(unsigned fBin, float fY, bool fDirection) const;
  float linearInterpolation (float fZ, const float fX[2], const float fY[2]) const;
  JetCorrectorParameters* mParameters;
};

#endif

