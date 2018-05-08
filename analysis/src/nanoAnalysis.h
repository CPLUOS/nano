#define Events_cxx
#include "nano/analysis/src/Events_2016v3.h"

#include <TH1D.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include <TString.h>

#include "pileUpTool.h"
#include "RoccoR.h"
#include "lumiTool.h"
#include "MuonScaleFactorEvaluator.h"
#include "ElecScaleFactorEvaluator.h"
#include "BTagWeightEvaluator.h"
#include "TopTriggerSF.h"
#include "TTbarModeDefs.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/DataLoader.h"
#include "TMVA/MethodCuts.h"


//void Events::Loop(){}

void makeEventsClass(const char* filedir){
  TFile *f = TFile::Open(filedir);
  TTree *t = (TTree*) f->Get("Events");
  t->MakeClass(); // this will generate Events.h file.
}

#ifdef nanoAnalysis_cxx
class nanoAnalysis : public Events {
private: 
      //Output Variables
      TFile* m_output;
      
      //Tree
      TTree* m_tree;
      
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
      Float_t b_jes_up, b_jes_dn;
      Float_t b_cferr1_up, b_cferr1_dn, b_cferr2_up, b_cferr2_dn;
      Float_t b_hf_up, b_hf_dn, b_hfstats1_up, b_hfstats1_dn, b_hfstats2_up, b_hfstats2_dn;
      Float_t b_lf_up, b_lf_dn, b_lfstats1_up, b_lfstats1_dn, b_lfstats2_up, b_lfstats2_dn;

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
      pileUpTool* m_pileUp;
      lumiTool* m_lumi;
      RoccoR* m_rocCor;
      MuonScaleFactorEvaluator m_muonSF;
      ElecScaleFactorEvaluator m_elecSF;
      BTagWeightEvaluator m_btagSF;
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
      Bool_t LumiCheck();
      enum TTLLChannel { CH_NOLL = 0, CH_MUEL, CH_ELEL, CH_MUMU };

      bool hasOverLap(TLorentzVector cand, vector<TParticle> objects, Float_t rad);
      
      std::vector<TParticle> MuonSelection();
      std::vector<TParticle> ElectronSelection(std::vector<TParticle>);
      std::vector<TParticle> JetSelection(std::vector<TParticle>, std::vector<TParticle>);
      std::vector<TParticle> BJetSelection(std::vector<TParticle>, std::vector<TParticle>);
      std::vector<TParticle> nonbJetSelection(std::vector<TParticle>, std::vector<TParticle>);
      
      Double_t roccoR(TLorentzVector m, Int_t &q, Int_t &nGen, Int_t &nTrackerLayers);

public:
  Bool_t m_isMC;

  //set output file
  void SetOutput(std::string outputName);
  void LoadModules(pileUpTool* pileUp, lumiTool* lumi, RoccoR* rocCor);
  nanoAnalysis(TTree *tree=0, Bool_t isMc = false);
  ~nanoAnalysis();
  virtual void     Loop();
};

nanoAnalysis::nanoAnalysis(TTree *tree, Bool_t isMC) : Events(tree), m_isMC(isMC)
{
}

nanoAnalysis::~nanoAnalysis()
{
  m_output->Write();
  m_output->Close();
}
#endif

#ifdef vtsAnalysis_cxx
class vtsAnalysis : public Events {
private:
  TFile* outFile; TTree* outTree;

  TH1D* h_cutFlow;
  TParticle recolep1, recolep2;
  TLorentzVector recolep1_tlv, recolep2_tlv;
  std::vector<TLorentzVector> recoleps;
  //std::vector<TLorentzVector> recoleps[recolep1_tlv,recolep2_tlv];
  //std::vector<TLorentzVector> recoleps(2);

  int b_channel;
  int b_njet;
  float b_met;
  TLorentzVector b_dilep_tlv;
  enum TTLLChannel { CH_NOLL = 0, CH_MUEL, CH_ELEL, CH_MUMU };

  TLorentzVector b_had_tlv;
  int b_isFrom_had;
  float b_d_had, b_x_had;
  float b_lxy_had, b_lxySig_had, b_angleXY_had, b_angleXYZ_had, b_chi2_had, b_dca_had;
  float b_pt_had, b_eta_had, b_l3D_had, b_l3DSig_had, b_legDR_had, b_mass_had; 
  int b_pdgId_had; 
  float b_dau1_chi2_had, b_dau1_ipsigXY_had, b_dau1_ipsigZ_had, b_dau1_pt_had;
  float b_dau2_chi2_had, b_dau2_ipsigXY_had, b_dau2_ipsigZ_had, b_dau2_pt_had;

  float  b_btagCSVV2_Jet, b_btagCMVA_Jet, b_btagDeepB_Jet, b_btagDeepC_Jet;
  float  b_area_Jet, b_pt_Jet;
  int b_nConstituents_Jet, b_nElectrons_Jet, b_nMuons_Jet;

  std::map<unsigned int, int> qjMapForMC_;
  std::vector<int> qMC_;
  std::vector<int> genJet_;           
  std::vector<int> recoJet_;

  //functions
  void EventSelection();
  void ResetBranch();

  void MatchingForMC();
  void HadronAnalysis();

  TParticle GetTParticle(int pdgId, int idx);

  Double_t DeltaR(Double_t deta, Double_t dphi); 
  Double_t DeltaPhi(Double_t phi1, Double_t phi2);
  Double_t GetD(float pt, float eta, float phi, float m, float vx, float vy, float vz);

  std::vector<TParticle> muonSelection();
  std::vector<TParticle> elecSelection();
  std::vector<TParticle> jetSelection();

  struct HadStat {
    int idx = -1;
    int pdgId = -99;
    float x = -1;
    int label = -9;
    int jetIdx = -99;
  };

public:
  Bool_t isMC_ = false;
  void MakeTree(std::string outputName);
  vtsAnalysis(TTree *tree=0);
  ~vtsAnalysis();
  virtual void     Loop();
};
vtsAnalysis::vtsAnalysis(TTree *tree) : Events(tree)
{ }
vtsAnalysis::~vtsAnalysis(){ }
#endif

#ifdef topAnalysis_cxx
class topAnalysis : public Events {
private: 
  TFile* m_output;    
  //Tree
  TTree* m_tree;
    
  //histogram
  TH1D* h_nevents;
  TH1D* h_genweights;
  TH1D* h_weights;
  TH1D* h_cutFlow;
    
  //Variables
  TLorentzVector b_lep1, b_lep2, b_dilep, b_jet1, b_jet2;
  TParticle recolep1, recolep2;
  int b_lep1_pid, b_lep2_pid;
  float b_jet1_CSVInclV2, b_jet2_CSVInclV2;

  std::vector<Float_t> b_csvweights;

  int b_nvertex, b_step, b_channel, b_njet, b_nbjet;
  bool b_step1, b_step2, b_step3, b_step4, b_step5, b_step6, b_step7;  
  float b_met, b_weight, b_genweight, b_puweight, b_btagweight;
  float b_mueffweight, b_mueffweight_up, b_mueffweight_dn,
        b_eleffweight, b_eleffweight_up, b_eleffweight_dn;
  float b_tri, b_tri_up, b_tri_dn;
  float b_cme_dca, b_cme_angleXY, b_cme_angleXYZ, b_cme_jetDR, b_cme_legDR;
  float b_cme_lxy, b_cme_lxyE, b_cme_l3D, b_cme_l3DE;
  float b_cme_x, b_cme_y, b_cme_z, b_cme_pt, b_cme_chi2, b_cme_eta, b_cme_phi;
  float b_cme_dau1_chi2, b_cme_dau1_nHits, b_cme_dau1_pt, b_cme_dau1_ipsigXY, b_cme_dau1_ipsigZ;
  float b_cme_dau2_chi2, b_cme_dau2_nHits, b_cme_dau2_pt, b_cme_dau2_ipsigXY, b_cme_dau2_ipsigZ;
  float b_cme_jet_btagCMVA, b_cme_jet_btagCSVV2, b_cme_jet_btagDeepB, b_cme_jet_btagDeepC;
  float b_cme_mass;
  float b_cme_tmva_bdtg;
  float b_cme_pdgId;
  int b_cme_nMatched;
  float b_bdtg;
  int b_maxbIdx;
  

  //Triggers
  Bool_t b_trig_m, b_trig_m2,  b_trig_e, b_trig_mm, b_trig_em, b_trig_ee;
  //Tools
  TH1D* hist_mc;
  MuonScaleFactorEvaluator muonSF_;
  ElecScaleFactorEvaluator elecSF_;
  BTagWeightEvaluator m_btagSF;
  pileUpTool *m_pileUp;
  lumiTool* m_lumi;

  //Making output branch
  void MakeBranch(TTree* t);
  void resetBranch();
  bool analysis();

  //For Selection
  Bool_t lumiCheck();
  std::vector<TParticle> muonSelection();
  std::vector<TParticle> elecSelection();
  std::vector<TLorentzVector> recoleps;
  std::vector<TParticle> jetSelection();
  std::vector<TParticle> bjetSelection();

  //TMVA
  TMVA::Reader* bdtg;

public :
  Bool_t m_isMC, m_isDL, m_isSL_e, m_isSL_m;

  //set output file
  void setOutput(std::string outputName);
  void LoadModules(pileUpTool* pileUp, lumiTool* lumi);
  void collectTMVAvalues();
  topAnalysis(TTree *tree=0, Bool_t flag = false, Bool_t dl = false, Bool_t sle = false, Bool_t slm = false);
  ~topAnalysis();
  virtual void     Loop();

};

topAnalysis::topAnalysis(TTree *tree, Bool_t flag, Bool_t dl, Bool_t sle, Bool_t slm) : Events(tree),  m_isMC(flag), m_isDL(dl), m_isSL_e(sle), m_isSL_m(slm)
{
}


topAnalysis::~topAnalysis()
{
  m_output->Write();
  m_output->Close();
}
#endif

