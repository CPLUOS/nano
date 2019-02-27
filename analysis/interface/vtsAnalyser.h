#ifndef vtsAnalyser_H
#define vtsAnalyser_H

#include "hadAnalyser.h"

class vtsAnalyser : public hadAnalyser {
public:
  float CONESIZE = 0.4; float m_xCut = 0.2; unsigned int m_jetDauArrSize = 350;
  std::vector<TParticle> m_selectedJet;
  int m_nj = 0; int m_nj2 = 0;
  int m_ej = -1; int m_jj = -2;
  int m_nej = 0; int m_njj = 0;

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
  int b_ntotjet;

  /* for MatchingForMC() */
  std::map<unsigned int, int> m_qjMapForMC;  // +-3 : jet matched to s-quark, +-5 : jet matched to b-quark, -1 : jet matched to two quarks (overlap)

  /// Additional information about global properties of the jet,
  /// relevant to our Vts analysis
  struct jetInfo {
    unsigned int idx; /// jet idx
    double pt; /// jet pt
    double drsj; /// DeltaR(s,jet)
    double drbj; /// DeltaR(b,jet)
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

  struct ksInfo {
    unsigned int idx; // ks idx
    double dr; // DeltaR(ks, jet)
    double x; // ks_pt/jet_pt
  };

  std::vector<jetInfo> m_jetDeltaRs;
  std::vector<genInfo> m_tqMC;   
  std::vector<genInfo> m_wqMC;

  int   b_tq1_idx, b_tq1_pdgId;
  int   b_tq2_idx, b_tq2_pdgId;
  TLorentzVector b_tq1_tlv;
  TLorentzVector b_tq2_tlv;
  int   b_tq1_matched_jidx, b_tq1_matched_isOverlap, b_tq1_matched_dr, b_tq1_matched_x;
  int   b_tq2_matched_jidx, b_tq2_matched_isOverlap, b_tq2_matched_dr, b_tq2_matched_x;

  /* for JetAnalysis() */
  float b_Jet_axis1, b_Jet_axis2, b_Jet_cpt1, b_Jet_cpt2, b_Jet_cpt3, b_Jet_npt1, b_Jet_npt2, b_Jet_npt3, b_Jet_ptD, b_Jet_delta;
  int   b_Jet_nmult, b_Jet_cmult;
  std::vector<int> b_Jet_isCorrectMat;

  /* for FillJetTreeForTMVA() */
  int   b_isSJet, b_isBJet, b_isOverlap, b_isHighest, b_isClosestToLep;
  float b_cmult,  b_nmult; // Data type was changed since JKS (and Jet) BDT require to provide variables for TMVA::reader as float
  float b_pt,     b_eta,    b_phi,       b_mass;
  float b_c_x1,   b_c_x2,   b_c_x3;
  float b_n_x1,   b_n_x2,   b_n_x3;
  float b_axis1,  b_axis2,  b_ptD,       b_area;
  float b_CSVV2,  b_hadronFlavour;
  float b_dr1,    b_dr2;

  float b_dau_pt[350], b_dau_eta[350], b_dau_phi[350];  
  int   b_dau_charge[350];

    /* for TMVA with KS info. */
  float b_Jet_bdt_score,       b_JKS_bdt_score;

  int   b_Ks_pdgId;
  int   b_Ks_idx,              b_Ks_nMatched,         b_Ks_isFrom;
  bool  b_Ks_isHadFromTop,     b_Ks_isHadFromW,       b_Ks_isHadFromS,  b_Ks_isHadFromC, b_Ks_isHadFromB;
  float b_Ks_d,                b_Ks_pt,               b_Ks_eta,         b_Ks_phi,        b_Ks_mass;
  float b_Ks_lxy,              b_Ks_lxySig,           b_Ks_l3D,         b_Ks_l3DSig,     b_Ks_legDR;
  float b_Ks_angleXY,          b_Ks_angleXYZ,         b_Ks_chi2,        b_Ks_dca;
  float b_Ks_dau1_chi2,        b_Ks_dau1_ipsigXY,     b_Ks_dau1_ipsigZ, b_Ks_dau1_pt;
  float b_Ks_dau2_chi2,        b_Ks_dau2_ipsigXY,     b_Ks_dau2_ipsigZ, b_Ks_dau2_pt;
  float b_Ks_bdt_score;
  float b_Ks_dr,               b_Ks_x;
 
  int b_jet_start, b_jet_end; 
  int b_had_start, b_had_end;

  /* functions */
  void ResetBranch();
  void MakeBranch();
  void MakeBranchOfHadron(TTree* tr);
  void MatchingForMC();
  void GenAnalysis();
  void RecAnalysis();
  void JetAnalysis();

  bool isGenFrom(int count, int idx, int & isFrom, bool & isFromTop, bool & isFromW, bool & isFromKstar);

  void FillTMVATrees();
  void IdentifyJet(unsigned int jetIdx, unsigned int sIdx, unsigned int bIdx);
  int  FindMatchedHadron(TLorentzVector jet_tlv);
  void SetJetValues(int i);
  void SetHadronValues(int i);

  void SetMVAReader();
  void ResetHadTree();
  void ResetJetTree();

};

vtsAnalyser::vtsAnalyser(TTree *tree, TTree *had, TTree *hadTruth, Bool_t isMC, Bool_t dl, Bool_t sle, Bool_t slm) :
  hadAnalyser(tree, had, hadTruth, isMC, dl, sle, slm)
{
}

#endif
