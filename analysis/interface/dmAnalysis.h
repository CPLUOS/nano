#ifndef dmAnalysis_H
#define dmAnalysis_H
#include "nanoBase.h"

class dmAnalysis : public nanoBase { 
private:

  TH1D* h_nevents;  
  TH1D* h_genweights;  
  TH1D* h_weights;  
 
  TLorentzVector b_mu1, b_mu2, b_dimu, b_met1, b_mume, b_elec, b_jet1, b_jet2, b_dijet;
  TParticle recolep1, recolep2, recolep, recomet, recojet1, recojet2;

  std::vector<Float_t> b_csvweights; 

  int b_nvertex, b_step, b_channel, b_nelec, b_ntau, b_njet, b_nbjet;
  bool b_step1, b_step2, b_step3, b_step4, b_step5, b_step6, b_step7;
  float b_met, b_weight, b_genweight, b_puweight,
        b_mpt, b_eta2, b_dEta;


  Bool_t b_trigger_dm;


  void MakeBranch(TTree* t);
  void resetBranch();
  bool analysis();

  enum DMChannel { CH_NOLL = 0, CH_ZJets, CH_WJets };

  Bool_t LumiCheck();
  std::vector<TParticle> muonSelection();
  std::vector<TParticle> elecvetoSelection();
  std::vector<TParticle> tauvetoSelection();
  std::vector<TLorentzVector> recoleps;
  std::vector<TLorentzVector> recojets;
  std::vector<TParticle> jetSelection();
  std::vector<TParticle> bjetvetoSelection();
  //std::vector<TParticle> MET();

public:

  dmAnalysis(TTree *tree=0, Bool_t isMC = false);
  ~dmAnalysis();
  void setOutput(std::string outputName);
  virtual void Loop();

};

dmAnalysis::dmAnalysis(TTree *tree, Bool_t isMC) : nanoBase(tree, 0, 0, isMC)
{
  std::string env = getenv("CMSSW_BASE");
}

dmAnalysis::~dmAnalysis()
{
}
#endif
