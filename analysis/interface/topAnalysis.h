#ifndef topAnalysis_H
#define topAnalysis_H

#include "nanoAnalysis.h"
#include "nano/external/interface/TopTriggerSF.h"
//#include "nano/external/interface/TTbarModeDefs.h"

class topAnalysis : public nanoAnalysis
{
protected:
  std::vector<Float_t> b_csvweights;
  float b_btagweight;
  //TParticle GetTParticle(int pdgId, int idx);

public:
  enum TTLLChannel { CH_NOLL = 0, CH_MUEL, CH_ELEL, CH_MUMU };

  //Bool_t lumiCheck();
  std::vector<TParticle> muonSelection();
  std::vector<TParticle> elecSelection();
  std::vector<TLorentzVector> recoleps;
  std::vector<TParticle> jetSelection();
  std::vector<TParticle> bjetSelection();

public:
  Bool_t m_isDL, m_isSL_e, m_isSL_m;

//  void setOutput(std::string outputName);
//  void LoadModules(pileUpTool* pileUp, lumiTool* lumi);
//  void collectTMVAvalues();
  topAnalysis(TTree *tree=0, Bool_t isMC = false, Bool_t dl = false, Bool_t sle = false, Bool_t slm = false);
  ~topAnalysis();  
  virtual void Loop() = 0;
};

#endif
