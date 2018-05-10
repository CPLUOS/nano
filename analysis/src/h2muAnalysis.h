#ifndef h2muAnalysis_H
#define h2muAnalysis_H
#include "nanoAnalysis.h"

class h2muAnalysis : public nanoAnalysis {
private: 
  //histogram
  TH1D* h_Event_Tot;
  TH1D* h_genweights;
  TH1D* h_weight;
  TH1D* h_XL;
  TH1D* h_nFH4;
  TH1D* h_Out;
  TH1D* h_Non;
 
  //Variables
  TLorentzVector b_Dilep, b_Mu1, b_Mu2;
  TLorentzVector b_lep;
  TLorentzVector b_nonbJet1, b_nonbJet2 ,b_nonbJet3 ,b_nonbJet4;

  std::vector<TLorentzVector> b_Mu_tlv, b_El_tlv, b_Jet_tlv, b_bJet_tlv, b_nonbJet_tlv;
  
  std::vector<Float_t> b_CSVv2;
  std::vector<Float_t> b_csvweights;
  
  Float_t b_genweight, b_weight;
  Float_t b_puweight, b_puweight_up, b_puweight_dn; 
  Float_t b_mueffweight, b_mueffweight_up, b_mueffweight_dn;
  Float_t b_btagweight, b_btagweight_up, b_btagweight_dn;

  Int_t b_Event_No, b_Event_Total;

  // BDT Variables
  Float_t b_Central_Jets, b_Forward_Jets; 
  // all
  Float_t b_Met_phi;
  Float_t b_Met;
  Float_t b_npvs;
  Float_t b_all_muEtaDiff, b_all_muPtDiff, b_all_muPhiDiff, b_all_muDR;
  Float_t b_all_Dilep_Pt, b_all_Dilep_Eta ,b_all_Dilep_Phi;
  Float_t b_DijetM1, b_DijetM2, b_DijetEta1, b_DijetEta2, DijetM_hold, DijetEta_hold;

  Float_t b_DiJetM12, b_DiJetM13, b_DiJetM14 ,b_DiJetM23 , b_DiJetM24, b_DiJetM34;
  Float_t b_minDR1, b_minDR2, b_minDR, b_XlepPt, b_mT2, b_mT;
  Float_t b_etaJ1, b_etaJ2;
  Float_t b_MVA_BDTXL, b_MVA_BDTFH, b_MVA_BDTnoB, b_MVA_BDTOut;
  Float_t b_CSV;
  Float_t DR_Hold, DR_Hold2, MuPT_Hold ,ElPT_Hold, mT_Hold;
  //Step and Cutflow
  TH1D* h_cutFlow;
  Int_t b_Step;
  
  //Channel
  Int_t b_channel;
  Float_t b_nlep, b_nmuon, b_nelec, b_njet, b_nbjet, b_charge, b_nnonbjet, b_nexLep;
  Bool_t keep;
  //triggers
  Bool_t b_trig_m, b_trig_m2,  b_trig_e, b_trig_mm, b_trig_em, b_trig_ee;
  Int_t b_FL, b_FH2, b_FH3, b_FH4, b_SL;
  Int_t b_XL, b_nFH4, b_Out, b_nonB;
  //weight files //
  std::string weightXL, weightFH, weightnoB, weightOut;
  //Tools
  std::vector<UInt_t> idxs;

  TMVA::Reader* bdt_XL;
  TMVA::Reader* bdt_FH;
  TMVA::Reader* bdt_noB;
  TMVA::Reader* bdt_Out;

  //Making output branch
  void MakeBranch(TTree* t);
  void ResetBranch();

  bool Analysis();
  //For Selection
  enum TTLLChannel { CH_NOLL = 0, CH_MUEL, CH_ELEL, CH_MUMU };

  bool hasOverLap(TLorentzVector cand, vector<TParticle> objects, Float_t rad);
  
  std::vector<TParticle> MuonSelection();
  std::vector<TParticle> ElectronSelection(std::vector<TParticle>);
  std::vector<TParticle> JetSelection(std::vector<TParticle>, std::vector<TParticle>);
  std::vector<TParticle> BJetSelection(std::vector<TParticle>, std::vector<TParticle>);
  std::vector<TParticle> nonbJetSelection(std::vector<TParticle>, std::vector<TParticle>);
  
  Double_t roccoR(TLorentzVector m, Int_t &q, Int_t &nGen, Int_t &nTrackerLayers);

public:
  //set output file
  h2muAnalysis(TTree *tree=0, Bool_t isMc = false);
  ~h2muAnalysis();
  void SetOutput(std::string outputName);
  virtual void Loop();
};

h2muAnalysis::h2muAnalysis(TTree *tree, Bool_t isMC) : nanoAnalysis(tree, isMC)
{
  string env = getenv("CMSSW_BASE");
  m_rocCor = new RoccoR(env+"/src/nano/analysis/data/rcdata.2016.v3/");
}

h2muAnalysis::~h2muAnalysis()
{
}

#endif
