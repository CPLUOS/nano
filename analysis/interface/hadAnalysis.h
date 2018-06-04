#ifndef hadAnalysis_H
#define hadAnalysis_H

#include "dilepTopAnalysis.h"

class hadAnalysis : public dilepTopAnalysis 
{
protected:
  //Vriable
  int b_chk = 0;
  TLorentzVector b_had_tlv;

  int b_isFrom_had;
  bool b_isHadJetMatched_had;
  float b_d_had, b_x_had, b_dr_had;
  float b_lxy_had, b_lxySig_had, b_angleXY_had, b_angleXYZ_had, b_chi2_had, b_dca_had;
  float b_pt_had, b_eta_had, b_l3D_had, b_l3DSig_had, b_legDR_had, b_mass_had; 
  int b_pdgId_had; 

  float b_dau1_chi2_had, b_dau1_ipsigXY_had, b_dau1_ipsigZ_had, b_dau1_pt_had;
  float b_dau2_chi2_had, b_dau2_ipsigXY_had, b_dau2_ipsigZ_had, b_dau2_pt_had;

  float  b_btagCSVV2_Jet, b_btagCMVA_Jet, b_btagDeepB_Jet, b_btagDeepC_Jet;
  float  b_area_Jet, b_pt_Jet;
  int b_nConstituents_Jet, b_nElectrons_Jet, b_nMuons_Jet;

public:
  //Struct
  struct HadStat {
    int idx = -1;
    int pdgId = -99;
    float x = -1;
    float dr = -1;
    int label = -99;
    int jetIdx = -99;
    bool isHadJetMatched = false;
  };

  struct JetStat {
    int idx = -1;
    float dr = -1;
    int matchedQuark = -99;
  };
  //Function
  Double_t GetD(float pt, float eta, float phi, float m, float vx, float vy, float vz);

  hadAnalysis(TTree *tree=0, Bool_t isMC = false, Bool_t dl = false, Bool_t sle = false, Bool_t slm = false);
  ~hadAnalysis();
  virtual void     Loop();
};
#endif
