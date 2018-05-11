#ifndef nanoAnalysis_H
#define nanoAnalysis_H

#include "Events.h"

#include <TH1D.h>
#include <TLorentzVector.h>
#include <TParticle.h>

#include <TString.h>

#include "nano/external/src/pileUpTool.h"
#include "nano/external/src/RoccoR.h"
#include "nano/external/src/lumiTool.h"

#include "nano/external/src/MuonScaleFactorEvaluator.h"
#include "nano/external/src/ElecScaleFactorEvaluator.h"

#include "nano/external/src/BTagCalibrationStandalone.h"
//#include "nano/external/src/TopTriggerSF.h"
//#include "nano/external/src/TTbarModeDefs.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/DataLoader.h"
#include "TMVA/MethodCuts.h"

class nanoAnalysis : public Events
{
public:
  nanoAnalysis(TTree *tree=0, Bool_t isMC=false);
  virtual ~nanoAnalysis();
  virtual void Loop() = 0;

  //Output Variables
  TFile* m_output;
  //Tree
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
