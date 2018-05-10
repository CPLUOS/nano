#ifndef topAnalysis_H
#define topAnalysis_H

#include "nanoAnalysis.h"
#include "TopTriggerSF.h"
//#include "TTbarModeDefs.h"

class topAnalysis : public nanoAnalysis
{
private:
  std::vector<Float_t> b_csvweights;
  float b_btagweight;

  TParticle GetTParticle(int pdgId, int idx);

public:
  enum TTLLChannel { CH_NOLL = 0, CH_MUEL, CH_ELEL, CH_MUMU };

  Bool_t lumiCheck();
  std::vector<TParticle> muonSelection();
  std::vector<TParticle> elecSelection();
  std::vector<TLorentzVector> recoleps;
  std::vector<TParticle> jetSelection();
  virtual std::vector<TParticle> bjetSelection() = 0;

public:
  Bool_t m_isDL, m_isSL_e, m_isSL_m;

//  void setOutput(std::string outputName);
//  void LoadModules(pileUpTool* pileUp, lumiTool* lumi);
//  void collectTMVAvalues();
  topAnalysis(TTree *tree=0, Bool_t isMC = false, Bool_t dl = false, Bool_t sle = false, Bool_t slm = false);
  ~topAnalysis();  
  virtual void Loop() = 0;
};

topAnalysis::topAnalysis(TTree *tree, Bool_t isMC, Bool_t dl, Bool_t sle, Bool_t slm) : nanoAnalysis(tree, isMC), m_isDL(dl), m_isSL_e(sle), m_isSL_m(slm)// Events(tree, isMC), m_isDL(dl), m_isSL_e(sle), m_isSL_m(slm)
{
}

topAnalysis::~topAnalysis()
{
}
#endif
