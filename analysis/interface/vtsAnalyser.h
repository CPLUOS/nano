#ifndef vtsAnalyser_H
#define vtsAnalyser_H

#include "hadAnalyser.h"

class vtsAnalyser : public hadAnalyser 
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

  vtsAnalyser(TTree *tree=0, TTree *had=0, TTree *hadTruth=0, Bool_t isMC = false, Bool_t dl = false, Bool_t sle = false, Bool_t slm = false);
  vtsAnalyser(TTree *tree=0, Bool_t isMC=false, Bool_t dl=false, Bool_t sle=false, Bool_t slm=false) : vtsAnalyser(tree, 0, 0, isMC, dl, sle, slm) {}
  ~vtsAnalyser();
  virtual void     Loop();
};

vtsAnalyser::vtsAnalyser(TTree *tree, TTree *had, TTree *hadTruth, Bool_t isMC, Bool_t dl, Bool_t sle, Bool_t slm) : hadAnalyser(tree, had, hadTruth, isMC, dl, sle, slm)
{ }

vtsAnalyser::~vtsAnalyser()
{ }

#endif
