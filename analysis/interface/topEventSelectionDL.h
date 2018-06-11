#ifndef topEventSelectionDL_H
#define topEventSelectionDL_H

#include "topObjectSelection.h"
#include "nano/external/interface/TopTriggerSF.h"
//#include "nano/external/interface/TTbarModeDefs.h"

class topEventSelectionDL : public topObjectSelection
{
protected:
  //Histogram
  TH1D* h_nevents;
  TH1D* h_genweights;
  TH1D* h_weights;
  TH1D* h_cutFlow;

  //Variables
  TLorentzVector b_lep1, b_lep2, b_dilep, b_jet1, b_jet2;
  TParticle recolep1, recolep2;
  int b_lep1_pid, b_lep2_pid;
  float b_jet1_CSVInclV2, b_jet2_CSVInclV2;

  int b_nvertex, b_step, b_channel, b_njet, b_nbjet;
  bool b_step1, b_step2, b_step3, b_step4, b_step5, b_step6, b_step7;  
  float b_met, b_weight, b_genweight, b_puweight;
  float b_mueffweight, b_mueffweight_up, b_mueffweight_dn,
        b_eleffweight, b_eleffweight_up, b_eleffweight_dn;
  float b_tri, b_tri_up, b_tri_dn;

  //Triggers
  Bool_t b_trig_m, b_trig_m2,  b_trig_e, b_trig_mm, b_trig_em, b_trig_ee;

  //Tools
//  TH1D* hist_mc;
  MuonScaleFactorEvaluator muonSF_;
  ElecScaleFactorEvaluator elecSF_;

  //For C meson
  std::vector<TLorentzVector> d0s;
  TLorentzVector b_d0;
  std::vector<float> b_d0_lepSV_lowM;
  std::vector<float> b_d0_lepSV_correctM;
  std::vector<float> b_d0_lepSV_dRM;
  
  Bool_t m_isDL, m_isSL_e, m_isSL_m;

public:

  enum TTLLChannel { CH_NOLL = 0, CH_MUEL, CH_ELEL, CH_MUMU };

  virtual int EventSelection();

  topEventSelectionDL(TTree *tree=0, TTree *had=0, TTree *hadTruth=0, Bool_t isMC = false, Bool_t dl = true, Bool_t sle = false, Bool_t slm = false);
  topEventSelectionDL(TTree *tree=0, Bool_t isMC=false, Bool_t dl=true, Bool_t sle=false, Bool_t slm=false) : topEventSelectionDL(tree, 0, 0, isMC, dl, sle, slm) {}
  ~topEventSelectionDL();  
  virtual void Loop() = 0;

  void Reset();
};

#endif
