#ifndef topObjectSelection_H
#define topObjectSelection_H

#include "nanoBase.h"
#include "nano/external/interface/TopTriggerSF.h"
//#include "nano/external/interface/TTbarModeDefs.h"

class topObjectSelection : public nanoBase
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

  topObjectSelection(TTree *tree=0, TTree *had=0, TTree *hadTruth=0, Bool_t isMC = false, Bool_t isDilep = true, Bool_t isSemiLep = false);
  topObjectSelection(TTree *tree=0, Bool_t isMC=false, Bool_t isDilep=true, Bool_t isSemiLep=false) : topObjectSelection(tree, 0, 0, isMC, isDilep, isSemiLep) {}
  ~topObjectSelection() {}
};

#endif
