#include "TopTriggerSF.h"
#include "TTbarModeDefs.h"

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
  vtsAnAlysis(TTree *tree=0);
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
  Bool_t m_isDL, m_isSL_e, m_isSL_m;

  //set output file
  void setOutput(std::string outputName);
  void LoadModules(pileUpTool* pileUp, lumiTool* lumi);
  void collectTMVAvalues();
  topAnalysis(TTree *tree=0, Bool_t isMC = false, Bool_t dl = false, Bool_t sle = false, Bool_t slm = false);
  ~topAnalysis();
  virtual void     Loop();

};

topAnalysis::topAnalysis(TTree *tree, Bool_t isMC, Bool_t dl, Bool_t sle, Bool_t slm) : Events(tree, isMC), m_isDL(dl), m_isSL_e(sle), m_isSL_m(slm)
{
}


topAnalysis::~topAnalysis()
{
  m_output->Write();
  m_output->Close();
}
#endif
