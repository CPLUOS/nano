#ifndef vtsAnalysis_H
#define vtsAnalysis_H

#include "hadAnalysis.h"

class vtsAnalysis : public hadAnalysis 
{
private:
  std::map<unsigned int, int> qjMapForMC_;
  std::vector<int> qMC_;
  std::vector<int> genJet_;           
  std::vector<struct JetStat> recoJet_;
  TTree *hadt = 0;

  //functions
  void ResetBranch();
  void MakeBranch(TTree* t);
  void MatchingForMC();
  void HadronAnalysis();
  void hadronAnalysisWithHadTruth();

public:
  void setOutput(std::string outputName);

  vtsAnalysis(TTree *tree = 0, Bool_t isMC = false, Bool_t dl = false, Bool_t sle = false, Bool_t slm = false) : hadAnalysis(tree, isMC, dl, sle, slm) {}
  vtsAnalysis(TTree *nano = 0, Bool_t isMC = false, TTree *had = 0, TTree *hadTruth = 0, Bool_t dl = false, Bool_t sle = false, Bool_t slm = false) : hadAnalysis(nano, isMC, had, hadTruth, dl, sle, slm) {}
  ~vtsAnalysis() {}
  virtual void Loop();
};

#endif
