#ifndef topObjectSelection_H
#define topObjectSelection_H

#include "nanoAnalyser.h"
#include "nano/external/interface/TopTriggerSF.h"
//#include "nano/external/interface/TTbarModeDefs.h"

class topObjectSelection : public nanoAnalyser
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

  topObjectSelection(TTree *tree=0, Bool_t isMC = false, Bool_t isDilep = true, Bool_t isSemiLep = false);
  ~topObjectSelection() {}
};

#endif
