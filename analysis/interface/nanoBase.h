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

#include "nano/external/interface/JetCorrectionUncertainty.h"
#include "nano/external/interface/JetCorrectorParameters.h"
#include "nano/external/interface/JetResolution.h"
#include <TRandom3.h>

#include "nano/external/interface/BTagCalibrationStandalone.h"
//#include "nano/external/interface/TopTriggerSF.h"
//#include "nano/external/interface/TTbarModeDefs.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/DataLoader.h"
#include "TMVA/MethodCuts.h"

class nanoBase : public Events
{
public:
  nanoBase(TTree *tree=0, TTree *had=0, TTree *hadTruth=0, Bool_t isMC=false);
  nanoBase(TTree *tree, Bool_t isMC=false) : nanoBase(tree,0,0,isMC) {}
  virtual ~nanoBase();
  virtual void Loop() = 0;
  
  // In uncertainty study we need to switch the kinematic variables of jets
  // The following variables are for this switch
  // In the topObjectSelection.cc these variables are used instead of Jet_pt, Jet_mass, and so on.
  // In default, these are same as the original ones, but when a user wants to study systematic uncertainty 
  // so that he/she needs to switch them to the evaluated ones, 
  // just touching them in anlalyser class will be okay, and this is for it.
  virtual void GetJetMassPt(UInt_t nIdx, 
    Float_t &fJetMass, Float_t &fJetPt, Float_t &fJetEta, Float_t &fJetPhi, UInt_t unFlag);
  
  Int_t GetMatchGenJet(UInt_t nIdxJet, Float_t fResolution);

  // Output Variables
  TFile* m_output;
  // Tree
  TTree* m_tree;

  pileUpTool* m_pileUp;
  lumiTool* m_lumi;
  RoccoR* m_rocCor;
  MuonScaleFactorEvaluator m_muonSF;
  ElecScaleFactorEvaluator m_elecSF;
  BTagCalibrationReader m_btagSF, m_btagSF_up, m_btagSF_dn;
  JetCorrectionUncertainty *m_jecUnc;
  JMENano::JetResolution m_jetResObj;
  JMENano::JetResolutionScaleFactor m_jetResSFObj;
  TRandom3 *m_rndEngine;
  
  Float_t m_fDRcone_JER;
  Float_t m_fResFactorMathcer;

  Bool_t m_isMC;
  
  enum {
    OptBit_JER_Up = 0, 
    OptBit_JER_Dn, 
    OptBit_JES_Up, 
    OptBit_JES_Dn
  };
  
  enum {
    OptFlag_JER_Up = ( 1 << OptBit_JER_Up ), 
    OptFlag_JER_Dn = ( 1 << OptBit_JER_Dn ), 
    OptFlag_JES_Up = ( 1 << OptBit_JES_Up ), 
    OptFlag_JES_Dn = ( 1 << OptBit_JES_Dn )
  };
};

#endif
