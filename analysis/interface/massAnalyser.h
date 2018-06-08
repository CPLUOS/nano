#ifndef massAnalyser_H
#define massAnalyser_H

#include "topEventSelectionDL.h"

class massAnalyser : public topEventSelectionDL {
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
  std::vector<TLorentzVector> d0s;
  TLorentzVector b_d0;
  std::vector<float> b_d0_lepSV_lowM;
  std::vector<float> b_d0_lepSV_correctM;
  std::vector<float> b_d0_lepSV_dRM;
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

  massAnalyser(TTree *tree=0, TTree *had=0, TTree *hadTruth=0, Bool_t isMC = false, Bool_t dl = false, Bool_t sle = false, Bool_t slm = false);
  massAnalyser(TTree *tree, Bool_t isMC, Bool_t dl=false, Bool_t sle=false, Bool_t slm=false) : massAnalyser(tree, 0, 0, isMC, dl, sle, slm) {}
  ~massAnalyser();
  virtual void     Loop();

};

#endif
