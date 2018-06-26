#ifndef hadAnalyser_H
#define hadAnalyser_H

#include "topEventSelectionDL.h"

class hadAnalyser : public topEventSelectionDL 
{
protected:
  //Vriable
  int b_chk = 0;
  TLorentzVector b_had_tlv;

  int b_had_isFrom;
  bool b_had_isHadJetMatched;
  float b_had_d, b_had_x, b_had_dr;
  float b_had_lxy, b_had_lxySig, b_had_angleXY, b_had_angleXYZ, b_had_chi2, b_had_dca;
  float b_had_pt, b_had_eta, b_had_l3D, b_had_l3DSig, b_had_legDR, b_had_mass; 
  int b_had_pdgId; 

  float b_had_dau1_chi2, b_had_dau1_ipsigXY, b_had_dau1_ipsigZ, b_had_dau1_pt;
  float b_had_dau2_chi2, b_had_dau2_ipsigXY, b_had_dau2_ipsigZ, b_had_dau2_pt;

  float  b_Jet_btagCSVV2, b_Jet_btagCMVA, b_Jet_btagDeepB, b_Jet_btagDeepC;
  float  b_Jet_area, b_Jet_pt;
  int b_Jet_nConstituents, b_Jet_nElectrons, b_Jet_nMuons;

  int b_hadTruth_nMatched, b_hadTruth_nTrueDau;
  int b_hadTruth_isHadFromTsb;
  bool b_hadTruth_isHadFromTop, b_hadTruth_isHadFromW, b_hadTruth_isHadFromS, b_hadTruth_isHadFromC, b_hadTruth_isHadFromB;

public:
  //Struct
  struct HadStat {
    int idx = -1;
    int pdgId = -99;
    int label = -99;
    int jetIdx = -99;
    float x = -1;
    float dr = -1;
    bool isHadJetMatched = false;
    HadStat(int idx=-1, int pdgId=-99, int label=-99, int jetIdx=-99, float x=-1, float dr=-1, bool isHadJetMatched=false): idx(idx), pdgId(pdgId), label(label), jetIdx(jetIdx), x(x), dr(dr), isHadJetMatched(isHadJetMatched) {}
  };
  struct JetStat {
    int idx = -1;
    float dr = -1;
    int matchedQuark = -99;
    JetStat(int idx, float dr, int matchedQuark): idx(idx), dr(dr), matchedQuark(matchedQuark) {}
  };
  //Function
  Double_t GetD(float pt, float eta, float phi, float m, float vx, float vy, float vz);

  hadAnalyser(TTree *tree=0, TTree *had=0, TTree *hadTruth=0, Bool_t isMC = false, Bool_t dl = false, Bool_t sle = false, Bool_t slm = false);
  hadAnalyser(TTree *tree, Bool_t isMC, Bool_t dl, Bool_t sle, Bool_t slm) : hadAnalyser(tree, 0, 0, isMC, dl, sle, slm) {}
  ~hadAnalyser();
  virtual void     Loop();
};
#endif
