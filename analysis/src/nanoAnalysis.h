#define Events_cxx
#include "nano/analysis/src/Events.h"

#include <TH1D.h>
#include <TLorentzVector.h>
#include <TParticle.h>

#include "pileUpTool.h"
#include "RoccoR.h"
#include "lumiTool.h"
#include "MuonScaleFactorEvaluator.h"
#include "ElecScaleFactorEvaluator.h"
#include "BTagWeightEvaluator.h"

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
      TH1D* h_FL;
      TH1D* h_SL;
      TH1D* h_FH2;
      TH1D* h_FH3;
      TH1D* h_FH4;
      TH1D* h_Out;
 
      //Variables
      TLorentzVector b_Dilep, b_Mu1, b_Mu2;
      TLorentzVector b_lep1, b_lep2;

      std::vector<TLorentzVector> b_Mu_tlv, b_El_tlv, b_Jet_tlv, b_bJet_tlv;
      
      std::vector<Float_t> b_CSVv2;
      std::vector<Float_t> b_csvweights;
      
      Float_t b_genweight, b_weight;
      Float_t b_puweight, b_puweight_up, b_puweight_dn; 
      Float_t b_mueffweight, b_mueffweight_up, b_mueffweight_dn;
      Float_t b_btagweight, b_btagweight_up, b_btagweight_dn;
      Float_t b_lf_up, b_lf_dn, b_hfstats1_up, b_hfstats1_dn, b_hfstats2_up, b_hfstats2_dn;

      Int_t b_Event_No, b_Event_Total;
      
      Float_t b_Met_phi;
      Float_t b_Met;
      Int_t b_npvs;

      //Step and Cutflow
      TH1D* h_cutFlow;
      Int_t b_Step;
     
      //Channel
      Int_t b_channel, b_nlep, b_nmuon, b_nelec, b_njet, b_nbjet, b_charge;
      Bool_t keep;
      //triggers
      Bool_t b_trig_m, b_trig_m2,  b_trig_e, b_trig_mm, b_trig_em, b_trig_ee;
      Int_t b_FL, b_FH2, b_FH3, b_FH4, b_SL;

      //Tools
      pileUpTool* m_pileUp;
      lumiTool* m_lumi;
      RoccoR* m_rocCor;
      MuonScaleFactorEvaluator m_muonSF;
      ElecScaleFactorEvaluator m_elecSF;
      BTagWeightEvaluator m_btagSF;
      std::vector<UInt_t> idxs;
      
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
      std::vector<TParticle> BJetSelection(std::vector<TParticle>);
      
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
      TH1D* h_nEvents;

      Int_t b_nMuon;
      TLorentzVector b_Muon_tlv;

      //functions
      void Analysis();
      void ResetBranch();
      TTree* MakeTree();

public:
      void out(std::string outputName);
      vtsAnalysis(TTree *tree=0);
      ~vtsAnalysis();
      virtual void     Loop();
};
vtsAnalysis::vtsAnalysis(TTree *tree) : Events(tree)
{ }
vtsAnalysis::~vtsAnalysis(){ }
#endif
