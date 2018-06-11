#ifndef topEventSelectionSL_H
#define topEventSelectionSL_H

#include "topObjectSelection.h"
#include "nano/external/interface/TopTriggerSF.h"
//#include "nano/external/interface/TTbarModeDefs.h"

class topEventSelectionSL : public topObjectSelection
{
protected:
  //Histogram
  TH1D* h_nevents;
  TH1D* h_genweights;
  TH1D* h_weights;
  TH1D* h_cutFlow;

  TH1D *h_cutFlowEl, *h_cutFlowMu;

  //Variables
  TLorentzVector b_lep, b_jet1, b_jet2, b_jet3, b_jet4;
  TParticle recolep;
  int b_lep_pid;
  float b_jet1_CSVInclV2, b_jet2_CSVInclV2, b_jet3_CSVInclV2, b_jet4_CSVInclV2;

  int b_nvertex, b_step, b_channel, b_njet, b_nbjet;
  float b_met, b_weight, b_genweight, b_puweight;
  float b_mueffweight, b_mueffweight_up, b_mueffweight_dn,
        b_eleffweight, b_eleffweight_up, b_eleffweight_dn;
  float b_tri, b_tri_up, b_tri_dn;

  //Triggers
  Bool_t b_trig_m, b_trig_e;

  //Tools
    // TH1D* hist_mc;
  MuonScaleFactorEvaluator muonSF_;
  ElecScaleFactorEvaluator elecSF_;
  
  Bool_t m_isSL_e, m_isSL_m;

public:

  enum TTSLChannel { CH_NOLL = 0, CH_EL, CH_MU };

  int EventSelection();

  void Reset();

  topEventSelectionSL(TTree *tree=0, TTree *had=0, TTree *hadTruth=0, Bool_t isMC = false, Bool_t sle = false, Bool_t slm = false);
  topEventSelectionSL(TTree *tree=0, Bool_t isMC=false, Bool_t sle=false, Bool_t slm=false) : topEventSelectionSL(tree, 0, 0, isMC, sle, slm) {}
  ~topEventSelectionSL();  
  virtual void Loop() = 0;
};

#endif
