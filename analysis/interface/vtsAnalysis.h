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

  //functions
  void ResetBranch();
  void MatchingForMC();
  void HadronAnalysis();
  void MakeBranch(TTree* t);

public:
  void setOutput(std::string outputName);

  vtsAnalysis(TTree *tree=0, Bool_t isMC = false, Bool_t dl = false, Bool_t sle = false, Bool_t slm = false);
  ~vtsAnalysis();
  virtual void     Loop();
};

vtsAnalysis::vtsAnalysis(TTree *tree, Bool_t isMC, Bool_t dl, Bool_t sle, Bool_t slm) : hadAnalysis(tree, isMC, dl, sle, slm)
{ }

vtsAnalysis::~vtsAnalysis()
{ }

#endif
