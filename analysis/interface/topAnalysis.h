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

  bool isDilep, isSemiLep;
  
public:
  std::vector<TParticle> muonSelection();
  std::vector<TParticle> elecSelection();
  std::vector<TParticle> vetoMuonSelection();
  std::vector<TParticle> vetoElecSelection();
  std::vector<TLorentzVector> recoleps;
  std::vector<TParticle> jetSelection();
  std::vector<TParticle> bjetSelection();

  topAnalysis(TTree *tree=0, Bool_t isMC = false, Bool_t isDilep = true, Bool_t isSemiLep = false);
  ~topAnalysis() {}
};

#endif
