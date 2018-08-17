#ifndef vtsAnalyser_H
#define vtsAnalyser_H

#include "hadAnalyser.h"

class vtsAnalyser : public hadAnalyser {
public:
  float m_jetConeSize = 0.4; float m_xCut = 0.2;

  vtsAnalyser(TTree *tree=0, TTree *had=0, TTree *hadTruth=0, Bool_t isMC = false, Bool_t dl = false, Bool_t sle = false, Bool_t slm = false);
  vtsAnalyser(TTree *tree=0, Bool_t isMC=false, Bool_t dl=false, Bool_t sle=false, Bool_t slm=false) : vtsAnalyser(tree, 0, 0, isMC, dl, sle, slm) {}

  ~vtsAnalyser() {}

  void setOutput(std::string outputName);
  virtual void Loop();

private:
  TTree *m_hadtrForTMVA;
  TTree *m_jettrForTMVA;
  TMVA::Reader *m_hadReader;

  bool b_passedEvent;
  int b_nJet, b_nSelJet, b_nSelJetEv;

  /* for MatchingForMC() */
  std::map<unsigned int, int> m_qjMapForMC; std::map<unsigned int, int> m_qgjMapForMC;
  std::map<unsigned int, int> m_closestRecJetForLep1; std::map<unsigned int, int> m_closestRecJetForLep2; 
  std::map<unsigned int, int> m_closestGenJetForLep1; std::map<unsigned int, int> m_closestGenJetForLep2;
  std::vector<int> m_tqMC; std::vector<int> m_wqMC;

  /// Additional information about global properties of the jet,
  /// relevant to our Vts analysis
  struct jetInfo {
    int idx; /// jet idx
    double pt; /// jet pt
    double drsj; /// DeltaR(s,jet)
    double drbj; /// DeltaR(b,jet)
    double drl1j; /// DeltaR(lep1,jet)
    double drl2j; /// DeltaR(lep2,jet)
  };

  std::vector<jetInfo> m_genJet;
  std::vector<jetInfo> m_recJet;

  float b_Jet_dr_closest_s, b_Jet_dr_closest_b;
  float b_SelJet_dr_closest_s, b_SelJet_dr_closest_b;
  float b_GenJet_dr_closest_s, b_GenJet_dr_closest_b;

  bool b_GenSJet,             b_GenBJet,             b_GenBothJet,             b_RecSJet,             b_RecBJet,             b_RecBothJet;
  bool b_GenSJetClosestToLep, b_GenBJetClosestToLep, b_GenBothJetClosestToLep, b_RecSJetClosestToLep, b_RecBJetClosestToLep, b_RecBothJetClosestToLep;
  bool b_GenSJetIsHighest,    b_GenBJetIsHighest,    b_GenBothJetIsHighest,    b_RecSJetIsHighest,    b_RecBJetIsHighest,    b_RecBothJetIsHighest;

  /* for GenHaronAnalysis() */
  std::vector<int>   b_genHadron_isGenFrom_vec; std::vector<bool> b_genHadron_isGenFromTop_vec; std::vector<bool> b_genHadron_inVol_vec;
  std::vector<float> b_genHadron_d_vec; std::vector<float> b_genHadron_pt_vec; std::vector<float> b_genHadron_eta_vec; std::vector<float> b_genHadron_phi_vec; std::vector<float> b_genHadron_mass_vec;
  std::vector<float> b_genHadron_dau1_pt_vec; std::vector<float> b_genHadron_dau1_eta_vec; std::vector<float> b_genHadron_dau1_phi_vec; std::vector<int> b_genHadron_dau1_pdgId_vec;
  std::vector<float> b_genHadron_dau2_pt_vec; std::vector<float> b_genHadron_dau2_eta_vec; std::vector<float> b_genHadron_dau2_phi_vec; std::vector<int> b_genHadron_dau2_pdgId_vec;
  std::vector<float> b_genHadron_x_closest_j_vec;  std::vector<float> b_genHadron_dr_closest_j_vec;  std::vector<float> b_genHadron_x_highest_j_vec;  std::vector<float> b_genHadron_dr_highest_j_vec; 
  std::vector<float> b_genHadron_x_closest_gj_vec; std::vector<float> b_genHadron_dr_closest_gj_vec; std::vector<float> b_genHadron_x_highest_gj_vec; std::vector<float> b_genHadron_dr_highest_gj_vec; 
  std::vector<int>   b_genHadron_isClosestPair_xOrder_j_vec;  std::vector<int> b_genHadron_isHighestPair_xOrder_j_vec;
  std::vector<int>   b_genHadron_isClosestPair_xOrder_gj_vec; std::vector<int> b_genHadron_isHighestPair_xOrder_gj_vec; 

  std::vector<int>   b_nSJet_vec; std::vector<int>  b_nBJet_vec; std::vector<int>   b_nGenSJet_vec; std::vector<int>  b_nGenBJet_vec;

  /* for GenAnalysis() */
  std::vector<int>   b_GenPart_isGenFrom_vec; std::vector<bool> b_GenPart_isGenFromTop_vec; std::vector<bool> b_GenPart_isGenFromW_vec; std::vector<bool> b_GenPart_isFromKstar_vec;
  std::vector<float> b_GenPart_d_vec; std::vector<float> b_GenPart_pt_vec; std::vector<float> b_GenPart_eta_vec; std::vector<float> b_GenPart_phi_vec; std::vector<float> b_GenPart_mass_vec;
  std::vector<float> b_GenPart_x_closest_j_vec;  std::vector<float> b_GenPart_dr_closest_j_vec;  std::vector<float> b_GenPart_x_highest_j_vec;  std::vector<float> b_GenPart_dr_highest_j_vec;
  std::vector<float> b_GenPart_x_closest_gj_vec; std::vector<float> b_GenPart_dr_closest_gj_vec; std::vector<float> b_GenPart_x_highest_gj_vec; std::vector<float> b_GenPart_dr_highest_gj_vec; 
  std::vector<int>   b_GenPart_isClosestPair_xOrder_j_vec;  std::vector<int> b_GenPart_isHighestPair_xOrder_j_vec;
  std::vector<int>   b_GenPart_isClosestPair_xOrder_gj_vec; std::vector<int> b_GenPart_isHighestPair_xOrder_gj_vec;

  /* for RecAnalysis() */
  std::vector<int>   b_hadTruth_pdgId_vec; std::vector<int> b_hadTruth_nMatched_vec; std::vector<int> b_hadTruth_isFrom_vec;
  std::vector<bool>  b_hadTruth_isHadFromTop_vec; std::vector<bool> b_hadTruth_isHadFromW_vec; std::vector<bool> b_hadTruth_isHadFromS_vec; std::vector<bool> b_hadTruth_isHadFromC_vec; std::vector<bool> b_hadTruth_isHadFromB_vec;
  std::vector<float> b_hadTruth_d_vec; std::vector<float> b_hadTruth_pt_vec; std::vector<float> b_hadTruth_eta_vec; std::vector<float> b_hadTruth_phi_vec; std::vector<float> b_hadTruth_mass_vec;
  std::vector<float> b_hadTruth_lxy_vec; std::vector<float> b_hadTruth_lxySig_vec; std::vector<float> b_hadTruth_l3D_vec; std::vector<float> b_hadTruth_l3DSig_vec; std::vector<float> b_hadTruth_legDR_vec;
  std::vector<float> b_hadTruth_angleXY_vec; std::vector<float> b_hadTruth_angleXYZ_vec; std::vector<float> b_hadTruth_chi2_vec; std::vector<float> b_hadTruth_dca_vec;
  std::vector<float> b_hadTruth_dau1_chi2_vec; std::vector<float> b_hadTruth_dau1_ipsigXY_vec; std::vector<float> b_hadTruth_dau1_ipsigZ_vec; std::vector<float> b_hadTruth_dau1_pt_vec;
  std::vector<float> b_hadTruth_dau2_chi2_vec; std::vector<float> b_hadTruth_dau2_ipsigXY_vec; std::vector<float> b_hadTruth_dau2_ipsigZ_vec; std::vector<float> b_hadTruth_dau2_pt_vec;
  std::vector<float> b_hadTruth_x_closest_j_vec;  std::vector<float> b_hadTruth_dr_closest_j_vec;  std::vector<float> b_hadTruth_x_highest_j_vec;  std::vector<float> b_hadTruth_dr_highest_j_vec;
  std::vector<float> b_hadTruth_x_closest_gj_vec; std::vector<float> b_hadTruth_dr_closest_gj_vec; std::vector<float> b_hadTruth_x_highest_gj_vec; std::vector<float> b_hadTruth_dr_highest_gj_vec;
  std::vector<int>   b_hadTruth_isClosestPair_xOrder_j_vec;  std::vector<int> b_hadTruth_isHighestPair_xOrder_j_vec;
  std::vector<int>   b_hadTruth_isClosestPair_xOrder_gj_vec; std::vector<int> b_hadTruth_isHighestPair_xOrder_gj_vec;

  /* for JetAnalysis() */
  float b_Jet_axis1, b_Jet_axis2, b_Jet_cpt1, b_Jet_cpt2, b_Jet_cpt3, b_Jet_npt1, b_Jet_npt2, b_Jet_npt3, b_Jet_ptD, b_Jet_delta;
  int   b_Jet_nmult, b_Jet_cmult;
  std::vector<int> b_Jet_isCorrectMat;

  /* for CollectVar() */
  float b_MET_pt, b_MET_phi, b_MET_sumEt;

  /* for FillHadTreeForTMVA() */
  int   b_Rec_pdgId,            b_Rec_nMatched,         b_Rec_isFrom;
  bool  b_Rec_isHadFromTop,     b_Rec_isHadFromW,       b_Rec_isHadFromS,  b_Rec_isHadFromC, b_Rec_isHadFromB;
  float b_Rec_d,                b_Rec_pt,               b_Rec_eta,         b_Rec_phi,        b_Rec_mass;
  float b_Rec_lxy,              b_Rec_lxySig,           b_Rec_l3D,         b_Rec_l3DSig,     b_Rec_legDR;
  float b_Rec_angleXY,          b_Rec_angleXYZ,         b_Rec_chi2,        b_Rec_dca;
  float b_Rec_dau1_chi2,        b_Rec_dau1_ipsigXY,     b_Rec_dau1_ipsigZ, b_Rec_dau1_pt;
  float b_Rec_dau2_chi2,        b_Rec_dau2_ipsigXY,     b_Rec_dau2_ipsigZ, b_Rec_dau2_pt;
  float b_Rec_bdt_score_pythia, b_Rec_bdt_score_herwig;

  /* for FillJetTreeForTMVA() */
  int   b_isSJet, b_isBJet, b_isHighest, b_isClosestToLep;
  int   b_cmult,  b_nmult;
  float b_pt,     b_eta,    b_phi,       b_mass;
  float b_c_x1,   b_c_x2,   b_c_x3;
  float b_n_x1,   b_n_x2,   b_n_x3;
  float b_axis1,  b_axis2,  b_ptD,       b_area;
  float b_CSVV2;

    /* for TMVA with KS info. */
  int   b_KS_nMatched_pythia,         b_KS_isFrom_pythia;
  bool  b_KS_isHadFromTop_pythia,     b_KS_isHadFromW_pythia,       b_KS_isHadFromS_pythia,  b_KS_isHadFromC_pythia, b_KS_isHadFromB_pythia;
  float b_KS_d_pythia,                b_KS_pt_pythia,               b_KS_eta_pythia,         b_KS_phi_pythia,        b_KS_mass_pythia;
  float b_KS_lxy_pythia,              b_KS_lxySig_pythia,           b_KS_l3D_pythia,         b_KS_l3DSig_pythia,     b_KS_legDR_pythia;
  float b_KS_angleXY_pythia,          b_KS_angleXYZ_pythia,         b_KS_chi2_pythia,        b_KS_dca_pythia;
  float b_KS_dau1_chi2_pythia,        b_KS_dau1_ipsigXY_pythia,     b_KS_dau1_ipsigZ_pythia, b_KS_dau1_pt_pythia;
  float b_KS_dau2_chi2_pythia,        b_KS_dau2_ipsigXY_pythia,     b_KS_dau2_ipsigZ_pythia, b_KS_dau2_pt_pythia;
  float b_KS_best_bdt_pythia;
 
  int   b_KS_nMatched_herwig,         b_KS_isFrom_herwig;
  bool  b_KS_isHadFromTop_herwig,     b_KS_isHadFromW_herwig,       b_KS_isHadFromS_herwig,  b_KS_isHadFromC_herwig, b_KS_isHadFromB_herwig;
  float b_KS_d_herwig,                b_KS_pt_herwig,               b_KS_eta_herwig,         b_KS_phi_herwig,        b_KS_mass_herwig;
  float b_KS_lxy_herwig,              b_KS_lxySig_herwig,           b_KS_l3D_herwig,         b_KS_l3DSig_herwig,     b_KS_legDR_herwig;
  float b_KS_angleXY_herwig,          b_KS_angleXYZ_herwig,         b_KS_chi2_herwig,        b_KS_dca_herwig;
  float b_KS_dau1_chi2_herwig,        b_KS_dau1_ipsigXY_herwig,     b_KS_dau1_ipsigZ_herwig, b_KS_dau1_pt_herwig;
  float b_KS_dau2_chi2_herwig,        b_KS_dau2_ipsigXY_herwig,     b_KS_dau2_ipsigZ_herwig, b_KS_dau2_pt_herwig;
  float b_KS_best_bdt_herwig;
 
 
  int b_jet_start, b_jet_end; 
  int b_had_start, b_had_end;


  /* functions */
  void ResetBranch();
  void MakeBranch();
  void MatchingForMC();
  void HadronAnalysis();
  void GenHadronAnalysis();
  void GenAnalysis();
  void RecAnalysis();
  void JetAnalysis();

  bool isGenFrom(int count, int idx, int & isFrom, bool & isFromTop, bool & isFromW, bool & isFromKstar);
  void CollectVar();

  void ResetForTMVA();
  void FillJetTreeForTMVA();
  void FillHadTreeForTMVA();
  void SetMVAReader();
};

vtsAnalyser::vtsAnalyser(TTree *tree, TTree *had, TTree *hadTruth, Bool_t isMC, Bool_t dl, Bool_t sle, Bool_t slm) :
  hadAnalyser(tree, had, hadTruth, isMC, dl, sle, slm)
{
}

#endif
