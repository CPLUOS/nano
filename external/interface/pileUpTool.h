#ifndef PILEUPTOOL_H
#define PILEUPTOOL_H

#include <vector>
#include <numeric>

class pileUpTool {

  std::vector<float> m_weights;
  std::vector<float> m_weights_up;
  std::vector<float> m_weights_dn;

  static const std::vector<float> Moriond17MC;
  static const std::vector<float> Moriond17RD;
  static const std::vector<float> Moriond17RD_up;
  static const std::vector<float> Moriond17RD_dn;

 public:  
  pileUpTool();
  
  float getWeight(int nTrueInt, int sys = 0);
  
};

#endif
