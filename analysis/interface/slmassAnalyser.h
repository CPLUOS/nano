#ifndef slmassAnalyser_H
#define slmassAnalyser_H

#include "topEventSelectionSL.h"
#include <TClonesArray.h>

class slmassAnalyser : public topEventSelectionSL {
private: 
  // Varialbes


  float b_cme_dau1_chi2, b_cme_dau1_nHits, b_cme_dau1_pt; 
  float b_cme_dau2_chi2, b_cme_dau2_nHits, b_cme_dau2_pt; 
  float b_cme_l3DE, b_cme_dca, b_cme_jetDR, b_cme_legDR, b_cme_chi2, b_cme_jet_btagCSVV2;
  float b_cme_tmva, b_cme_diffMass, b_cme_mass, b_cme_sum;
  int b_cme_pdgId,b_cme_cut,b_cme_charge;
  float b_cme_x, b_cme_y, b_cme_z, b_cme_angleXY;
  float b_bdtg;
  int hadnum;

  float tmp_l3DE;
  float tmp_jetDR; 
  float tmp_legDR; 
  float tmp_dca; 
  float tmp_chi2; 
  float tmp_jet_btagCSVV2; 
  float tmp_mass;
  float tmp_tmva;

  //For C meson
  std::vector<TLorentzVector> hads;
  TClonesArray *b_hads;
  TLorentzVector b_had;
  //Making output branch
  void MakeBranch(TTree* t);
  void resetBranch();
  bool analysis();

  //TMVA
  TMVA::Reader* bdtg;
public :
  //set output file
  void setOutput(std::string outputName);
  void cmesonSelection();

  slmassAnalyser(TTree *tree=0, TTree *had=0, TTree *hadTruth=0, Bool_t isMC = false, Bool_t sle = false, Bool_t slm = false);
  slmassAnalyser(TTree *tree, Bool_t isMC,  Bool_t sle=false, Bool_t slm=false) : slmassAnalyser(tree, 0, 0, isMC, sle, slm) {}
  ~slmassAnalyser();
  virtual void     Loop();

  int SetCutValues();   


};

#endif
