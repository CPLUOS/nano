#ifndef vtsAnalyser_H
#define vtsAnalyser_H

#include "hadAnalyser.h"

class vtsAnalyser : public hadAnalyser 
{
private:
  // for Test()
  int nTotHadKS = 0; int nMatchedKS = 0; int nHadKSDauCut = 0; int nRealKSFromTop = 0; int nRealKSWithJet = 0; int nRealKSWithJetAndCut = 0; int nKSWithCut =0; int nKSWithDauCut =0;

  std::vector<float> b_hadTruth_pt_vec; std::vector<float> b_hadTruth_eta_vec; std::vector<float> b_hadTruth_phi_vec; std::vector<float> b_hadTruth_mass_vec;
  std::vector<float> b_hadTruth_lxy_vec; std::vector<float> b_hadTruth_lxySig_vec; std::vector<float> b_hadTruth_angleXY_vec; std::vector<float> b_hadTruth_angleXYZ_vec; std::vector<float> b_hadTruth_chi2_vec; std::vector<float> b_hadTruth_dca_vec;
  std::vector<float> b_hadTruth_l3D_vec; std::vector<float> b_hadTruth_l3DSig_vec; std::vector<float> b_hadTruth_legDR_vec; std::vector<int> b_hadTruth_pdgId_vec;
  std::vector<float> b_hadTruth_dau1_chi2_vec; std::vector<float> b_hadTruth_dau1_ipsigXY_vec; std::vector<float> b_hadTruth_dau1_ipsigZ_vec; std::vector<float> b_hadTruth_dau1_pt_vec;
  std::vector<float> b_hadTruth_dau2_chi2_vec; std::vector<float> b_hadTruth_dau2_ipsigXY_vec; std::vector<float> b_hadTruth_dau2_ipsigZ_vec; std::vector<float> b_hadTruth_dau2_pt_vec;
  std::vector<int> b_hadTruth_isFrom_vec;

  std::vector<int> b_hadTruth_isFrom_cut_vec; std::vector<float> b_hadTruth_x_cut_vec;
  std::vector<int> b_hadTruth_isFrom_nc_vec; std::vector<float> b_hadTruth_x_nc_vec;
  int b_hadTruth_isFrom;
  float b_hadTruth_x;
  std::vector<int> b_comb_isFrom_vec; std::vector<float> b_comb_x_vec;  
  int b_comb_isFrom;
  float b_comb_x;
  std::vector<int> b_had_isFrom_vec; std::vector<int> b_had_isFrom_vec_2; std::vector<float> b_had_x_vec;
  std::vector<int> b_had_isFrom_dc_vec; std::vector<int> b_had_isFrom_dc_vec_2; std::vector<float> b_had_x_dc_vec;


  std::map<unsigned int, int> qjMapForMC_;
  std::vector<int> tqMC_;
  std::vector<int> wqMC_;
  std::vector<int> genJet_;           
  std::vector<struct JetStat> recoJet_;

  //functions
  void ResetBranch();
  void MakeBranch(TTree* t);
  void MatchingForMC();
  void HadronAnalysis();

  int Test();

public:
  void setOutput(std::string outputName);

  vtsAnalyser(TTree *tree=0, TTree *had=0, TTree *hadTruth=0, Bool_t isMC = false, Bool_t dl = false, Bool_t sle = false, Bool_t slm = false);
  vtsAnalyser(TTree *tree=0, Bool_t isMC=false, Bool_t dl=false, Bool_t sle=false, Bool_t slm=false) : vtsAnalyser(tree, 0, 0, isMC, dl, sle, slm) {}
  ~vtsAnalyser();
  virtual void     Loop();
};

vtsAnalyser::vtsAnalyser(TTree *tree, TTree *had, TTree *hadTruth, Bool_t isMC, Bool_t dl, Bool_t sle, Bool_t slm) : hadAnalyser(tree, had, hadTruth, isMC, dl, sle, slm)
{ }

vtsAnalyser::~vtsAnalyser()
{ }

#endif
