#ifndef vtsAnalyser_H
#define vtsAnalyser_H

#include "hadAnalyser.h"

class vtsAnalyser : public hadAnalyser {
public:
  float m_jetConeSize = 0.4; float m_xCut = 0.2; unsigned int m_jetDauArrSize = 350;
  std::vector<TParticle> m_selectedJet;

  vtsAnalyser(TTree *tree=0, TTree *had=0, TTree *hadTruth=0, Bool_t isMC = false, Bool_t dl = false, Bool_t sle = false, Bool_t slm = false);
  vtsAnalyser(TTree *tree=0, Bool_t isMC=false, Bool_t dl=false, Bool_t sle=false, Bool_t slm=false) : vtsAnalyser(tree, 0, 0, isMC, dl, sle, slm) {}
  vtsAnalyser(TTree *tree=0, Bool_t isMC=false, Bool_t dl=false, Bool_t sle=false, Bool_t slm=false, Bool_t isGenericMC=false) : hadAnalyser(tree, 0, 0, isMC, dl, sle, slm), m_isGenericMC(isGenericMC) {}
  ~vtsAnalyser() {}

  void setOutput(std::string outputName);
  virtual void Loop();

private:
  Bool_t m_isGenericMC = false;

  TTree *m_hadtrForTMVA;
  TTree *m_jettrForTMVA;
  TMVA::Reader *m_hadReader;
  TMVA::Reader *m_jetReader;
  TMVA::Reader *m_jksReader;

  bool b_passedEvent;
  int b_nJet, b_nSelJet, b_nSelJetEv;

  /* for MatchingForMC() */
  std::map<unsigned int, int> m_qjMapForMC;  // +-3 : jet matched to s-quark, +-5 : jet matched to b-quark, -1 : jet matched to two quarks (overlap)
  std::map<unsigned int, int> m_qgjMapForMC;
  std::map<unsigned int, int> m_closestRecJetForLep1; std::map<unsigned int, int> m_closestRecJetForLep2; 
  std::map<unsigned int, int> m_closestGenJetForLep1; std::map<unsigned int, int> m_closestGenJetForLep2;

  /// Additional information about global properties of the jet,
  /// relevant to our Vts analysis

  struct jetInfo {
    int idx; /// jet idx
    double pt; /// jet pt
    double dr1j; /// DeltaR(s,jet) if bbbar sample, then also DeltaR(b, jet)
    double dr2j; /// DeltaR(b,jet)
    double drl1j; /// DeltaR(lep1,jet)
    double drl2j; /// DeltaR(lep2,jet)
  };

  /// information about global properties of the gen quark
  struct genInfo {
    int idx; // gen idx
    int pdgId; // gen pdgId
    int status; // gen status
    int mom_idx; // mom idx
    int mom_pdgId; // mom pdgId
    int mom_status; // mom status
    TLorentzVector tlv; // gen 4-vector
  };

  std::vector<jetInfo> m_genJet; 
  std::vector<jetInfo> m_recJet;
  std::vector<genInfo> m_tqMC;   
  std::vector<genInfo> m_wqMC;

  int   b_matched1_jidx,    b_matched1_isOverlap;
  int   b_matched1_idx,     b_matched1_pdgId,     b_matched1_status;
  int   b_matched1_mom_idx, b_matched1_mom_pdgId, b_matched1_mom_status;
  float b_matched1_dr,      b_matched1_x;
  TLorentzVector b_matched1_tlv;

  int   b_matched2_jidx,    b_matched2_isOverlap;
  int   b_matched2_idx,     b_matched2_pdgId,     b_matched2_status;
  int   b_matched2_mom_idx, b_matched2_mom_pdgId, b_matched2_mom_status;
  float b_matched2_dr,      b_matched2_x;
  TLorentzVector b_matched2_tlv;

  int b_RecSJet,             b_RecBJet;            
  int b_RecSJetClosestToLep, b_RecBJetClosestToLep;
  int b_RecSJetIsHighest,    b_RecBJetIsHighest;

  /* for RecAnalysis() */
  std::vector<int>   b_hadTruth_pdgId_vec;        std::vector<int>   b_hadTruth_nMatched_vec;      std::vector<int>   b_hadTruth_isFrom_vec;
  std::vector<bool>  b_hadTruth_isHadFromTop_vec; std::vector<bool>  b_hadTruth_isHadFromW_vec;    std::vector<bool>  b_hadTruth_isHadFromS_vec;   std::vector<bool>  b_hadTruth_isHadFromC_vec;   std::vector<bool>  b_hadTruth_isHadFromB_vec;
  std::vector<float> b_hadTruth_d_vec;            std::vector<float> b_hadTruth_pt_vec;            std::vector<float> b_hadTruth_eta_vec;          std::vector<float> b_hadTruth_phi_vec;          std::vector<float> b_hadTruth_mass_vec;
  std::vector<float> b_hadTruth_lxy_vec;          std::vector<float> b_hadTruth_lxySig_vec;        std::vector<float> b_hadTruth_l3D_vec;          std::vector<float> b_hadTruth_l3DSig_vec;       std::vector<float> b_hadTruth_legDR_vec;
  std::vector<float> b_hadTruth_angleXY_vec;      std::vector<float> b_hadTruth_angleXYZ_vec;      std::vector<float> b_hadTruth_chi2_vec;         std::vector<float> b_hadTruth_dca_vec;
  std::vector<float> b_hadTruth_dau1_chi2_vec;    std::vector<float> b_hadTruth_dau1_ipsigXY_vec;  std::vector<float> b_hadTruth_dau1_ipsigZ_vec;  std::vector<float> b_hadTruth_dau1_pt_vec;
  std::vector<float> b_hadTruth_dau2_chi2_vec;    std::vector<float> b_hadTruth_dau2_ipsigXY_vec;  std::vector<float> b_hadTruth_dau2_ipsigZ_vec;  std::vector<float> b_hadTruth_dau2_pt_vec;
  std::vector<float> b_hadTruth_x_closest_j_vec;  std::vector<float> b_hadTruth_dr_closest_j_vec;  std::vector<float> b_hadTruth_x_highest_j_vec;  std::vector<float> b_hadTruth_dr_highest_j_vec;
  std::vector<float> b_hadTruth_x_closest_gj_vec; std::vector<float> b_hadTruth_dr_closest_gj_vec; std::vector<float> b_hadTruth_x_highest_gj_vec; std::vector<float> b_hadTruth_dr_highest_gj_vec;
  std::vector<int>   b_hadTruth_isClosestPair_xOrder_j_vec;  std::vector<int> b_hadTruth_isHighestPair_xOrder_j_vec;
  std::vector<int>   b_hadTruth_isClosestPair_xOrder_gj_vec; std::vector<int> b_hadTruth_isHighestPair_xOrder_gj_vec;

  /* for CollectVar() */
  float b_MET_pt, b_MET_phi, b_MET_sumEt;

  /* for FillHadTreeForTMVA() */
  int   b_Rec_pdgId,        b_Rec_nMatched,     b_Rec_isFrom;
  bool  b_Rec_isHadFromTop, b_Rec_isHadFromW,   b_Rec_isHadFromS,   b_Rec_isHadFromC,   b_Rec_isHadFromB;
  float b_Rec_d,            b_Rec_pt,           b_Rec_eta,          b_Rec_phi,          b_Rec_mass;
  float b_Rec_lxy,          b_Rec_lxySig,       b_Rec_l3D,          b_Rec_l3DSig,       b_Rec_legDR;
  float b_Rec_angleXY,      b_Rec_angleXYZ,     b_Rec_chi2,         b_Rec_dca;
  float b_Rec_dau1_chi2,    b_Rec_dau1_ipsigXY, b_Rec_dau1_ipsigZ,  b_Rec_dau1_pt;     
  float b_Rec_dau2_chi2,    b_Rec_dau2_ipsigXY, b_Rec_dau2_ipsigZ,  b_Rec_dau2_pt;     
  float b_Rec_bdt_score;

  /* for FillJetTreeForTMVA() */
  int   b_isSJet, b_isBJet, b_isOverlap, b_isHighest, b_isClosestToLep;
  float b_cmult,  b_nmult; // Data type was changed since JKS (and Jet) BDT require to provide variables for TMVA::reader as float
  float b_pt,     b_eta,    b_phi,       b_mass;
  float b_c_x1,   b_c_x2,   b_c_x3;
  float b_n_x1,   b_n_x2,   b_n_x3;
  float b_axis1,  b_axis2,  b_ptD,       b_area;
  float b_CSVV2;
  float b_dr1,    b_dr2;

  float b_dau_pt[350], b_dau_eta[350], b_dau_phi[350];  
  int   b_dau_charge[350];

    /* for TMVA with KS info. */
  float b_Jet_bdt_score,       b_JKS_bdt_score;

  int   b_KS_idx,              b_KS_nMatched,         b_KS_isFrom;
  bool  b_KS_isHadFromTop,     b_KS_isHadFromW,       b_KS_isHadFromS,  b_KS_isHadFromC, b_KS_isHadFromB;
  float b_KS_d,                b_KS_pt,               b_KS_eta,         b_KS_phi,        b_KS_mass;
  float b_KS_lxy,              b_KS_lxySig,           b_KS_l3D,         b_KS_l3DSig,     b_KS_legDR;
  float b_KS_angleXY,          b_KS_angleXYZ,         b_KS_chi2,        b_KS_dca;
  float b_KS_dau1_chi2,        b_KS_dau1_ipsigXY,     b_KS_dau1_ipsigZ, b_KS_dau1_pt;
  float b_KS_dau2_chi2,        b_KS_dau2_ipsigXY,     b_KS_dau2_ipsigZ, b_KS_dau2_pt;
  float b_KS_dr,               b_KS_x,                b_KS_best_bdt;
 
  int b_jet_start, b_jet_end; 
  int b_had_start, b_had_end;

  /* functions */
  void jetOrdering(TString param, std::vector<vtsAnalyser::jetInfo>& jets);

  void ResetBranch();
  void MakeBranch();
  void MatchingForMC();
  void RecAnalysis();

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
