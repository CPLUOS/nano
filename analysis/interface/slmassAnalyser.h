#ifndef slmassAnalyser_H
#define slmassAnalyser_H

#include "topEventSelectionSL.h"

class slmassAnalyser : public topEventSelectionSL {
private: 
  //Varialbes
  float b_cme_dca, b_cme_angleXY, b_cme_angleXYZ, b_cme_jetDR, b_cme_legDR;
  float b_cme_lxy, b_cme_lxyE, b_cme_l3D, b_cme_l3DE;
  float b_cme_x, b_cme_y, b_cme_z, b_cme_pt, b_cme_chi2, b_cme_eta, b_cme_phi;
  float b_cme_dau1_chi2, b_cme_dau1_nHits, b_cme_dau1_pt, b_cme_dau1_ipsigXY, b_cme_dau1_ipsigZ;
  float b_cme_dau2_chi2, b_cme_dau2_nHits, b_cme_dau2_pt, b_cme_dau2_ipsigXY, b_cme_dau2_ipsigZ;
  float b_cme_jet_btagCMVA, b_cme_jet_btagCSVV2, b_cme_jet_btagDeepB, b_cme_jet_btagDeepC;
  float b_cme_mass;
  float b_cme_tmva_bdtg;
  int b_cme_pdgId;
  int b_cme_nMatched;
  float b_bdtg;
  int b_maxbIdx;
  //For C meson
  std::vector<TLorentzVector> hads;
  TLorentzVector b_had, b_had_vecSumDMLep;
  float b_had_sumM;
  //std::vector<float> b_d0_lepSV_lowM;
  //std::vector<float> b_d0_lepSV_correctM;
  //std::vector<float> b_d0_lepSV_dRM;
  //Making output branch
  void MakeBranch(TTree* t);
  void resetBranch();
  bool analysis();

  //TMVA
  TMVA::Reader* bdtg;

public :
  //set output file
  void setOutput(std::string outputName);
  void collectTMVAvalues();
  void cmesonSelection();

  slmassAnalyser(TTree *tree=0, TTree *had=0, TTree *hadTruth=0, Bool_t isMC = false, Bool_t sle = false, Bool_t slm = false);
  slmassAnalyser(TTree *tree, Bool_t isMC,  Bool_t sle=false, Bool_t slm=false) : slmassAnalyser(tree, 0, 0, isMC, sle, slm) {}
  ~slmassAnalyser();
  virtual void     Loop();

};

#endif
