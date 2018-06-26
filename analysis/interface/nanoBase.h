#ifndef nanoBase_H
#define nanoBase_H

#include "Events.h"

#include <TH1D.h>
#include <TLorentzVector.h>
#include <TParticle.h>

#include <TString.h>

#include "nano/external/interface/pileUpTool.h"
#include "nano/external/interface/RoccoR.h"
#include "nano/external/interface/lumiTool.h"

#include "nano/external/interface/MuonScaleFactorEvaluator.h"
#include "nano/external/interface/ElecScaleFactorEvaluator.h"

#include "nano/external/interface/BTagCalibrationStandalone.h"
//#include "nano/external/interface/TopTriggerSF.h"
//#include "nano/external/interface/TTbarModeDefs.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/DataLoader.h"
#include "TMVA/MethodCuts.h"

#define Branch_(type, name, suffix) t->Branch(#name, &(b_##name), #name "/" #suffix);
#define BranchI(name) Branch_(Int_t, name, I)
#define BranchF(name) Branch_(Float_t, name, F)
#define BranchO(name) Branch_(Bool_t, name, O)
#define BranchA_(type, name, size, suffix) t->Branch(#name, &(b_##name), #name"["#size"]/"#suffix);
#define BranchAI(name, size) BranchA_(Int_t, name, size, I);
#define BranchAF(name, size) BranchA_(Float_t, name, size, F);
#define BranchAO(name, size) BranchA_(Bool_t, name, size, O);
#define BranchVI(name) t->Branch(#name, "vector<int>", &(b_##name));
#define BranchVF(name) t->Branch(#name, "vector<float>", &(b_##name));
#define BranchVO(name) t->Branch(#name, "vector<bool>", &(b_##name));
#define BranchTLV(name) t->Branch(#name, "TLorentzVector", &(b_##name));

class nanoBase : public Events
{
public:
  nanoBase(TTree *tree=0, TTree *had=0, TTree *hadTruth=0, Bool_t isMC=false);
  nanoBase(TTree *tree, Bool_t isMC=false) : nanoBase(tree,0,0,isMC) {}
  virtual ~nanoBase();
  virtual void Loop() = 0;

  // Output Variables
  TFile* m_output;
  // Tree
  TTree* m_tree;

  pileUpTool* m_pileUp;
  lumiTool* m_lumi;
  RoccoR* m_rocCor;
  MuonScaleFactorEvaluator m_muonSF;
  ElecScaleFactorEvaluator m_elecSF;
  BTagCalibrationReader m_btagSF;

  Bool_t m_isMC;
};

#endif
