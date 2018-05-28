#ifndef topAnalysis_H
#define topAnalysis_H

#include "nanoAnalysis.h"
#include "nano/external/interface/TopTriggerSF.h"
//#include "nano/external/interface/TTbarModeDefs.h"

class topAnalysis : public nanoAnalysis
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

  std::vector<Float_t> b_csvweights;

  int b_nvertex, b_step, b_channel, b_njet, b_nbjet;
  float b_met, b_weight, b_genweight, b_puweight, b_btagweight;
  float b_mueffweight, b_mueffweight_up, b_mueffweight_dn,
        b_eleffweight, b_eleffweight_up, b_eleffweight_dn;
  float b_tri, b_tri_up, b_tri_dn;

  //Triggers
  Bool_t b_trig_m, b_trig_m2,  b_trig_e, b_trig_mm, b_trig_em, b_trig_ee;

  //Tools
  MuonScaleFactorEvaluator muonSF_;
  ElecScaleFactorEvaluator elecSF_;
  
  Bool_t m_isDL, m_isSL_e, m_isSL_m;

public:

  enum TTLLChannel { CH_NOLL = 0, CH_MUEL, CH_ELEL, CH_MUMU };

  int EventSelection();

  std::vector<TParticle> muonSelection();
  std::vector<TParticle> elecSelection();
  std::vector<TLorentzVector> recoleps;
  std::vector<TParticle> jetSelection();
  std::vector<TParticle> bjetSelection();

  topAnalysis(TTree *tree=0, Bool_t isMC = false, Bool_t dl = false, Bool_t sle = false, Bool_t slm = false);
  ~topAnalysis();  
  virtual void Loop() = 0;
  void Reset();
};

#endif
