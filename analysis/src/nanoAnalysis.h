#define Events_cxx
#include "nano/analysis/src/Events.h"

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
#include "TMVA/MethodCuts.h"

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
      Float_t b_lf_up, b_lf_dn, b_hfstats1_up, b_hfstats1_dn, b_hfstats2_up, b_hfstats2_dn;

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
  // TMVA Booking //
  weightXL = std::string(std::getenv("CMSSW_BASE"))+"/src/nano/analysis/test/h2mu/TMVA/newEvent/split/XL_50/weights/TMVAClassification_BDT.weights.xml";
  weightFH = std::string(std::getenv("CMSSW_BASE"))+"/src/nano/analysis/test/h2mu/TMVA/newEvent/split/nFH4_50/weights/TMVAClassification_BDT.weights.xml";
  weightOut = std::string(std::getenv("CMSSW_BASE"))+"/src/nano/analysis/test/h2mu/TMVA/newEvent/split/Out_100/weights/TMVAClassification_BDT.weights.xml";
  weightnoB = std::string(std::getenv("CMSSW_BASE"))+"/src/nano/analysis/test/h2mu/TMVA/newEvent/nonB_120/weights/TMVAClassification_BDT.weights.xml";

  bdt_XL = new TMVA::Reader();
  bdt_XL->AddVariable( "Met", &b_Met );
  bdt_XL->AddVariable( "all_muEtaDiff", &b_all_muEtaDiff );
  bdt_XL->AddVariable( "all_muPtDiff", &b_all_muPhiDiff );
  bdt_XL->AddVariable( "all_muPhiDiff", &b_all_muPhiDiff );
  bdt_XL->AddVariable( "all_muDR", &b_all_muDR );
  bdt_XL->AddVariable( "all_Dilep_Pt", &b_all_Dilep_Pt );
  bdt_XL->AddVariable( "all_Dilep_Eta", &b_all_Dilep_Eta );
  bdt_XL->AddVariable( "nelec", &b_nelec );
  bdt_XL->AddVariable( "nmuon", &b_nmuon );
  bdt_XL->AddVariable( "nnonbjet", &b_nnonbjet );
  bdt_XL->AddVariable( "nbjet", &b_nbjet );
  bdt_XL->AddVariable( "Central_Jets", &b_Central_Jets );
  bdt_XL->AddVariable( "Forward_Jets", &b_Forward_Jets );
  bdt_XL->AddVariable( "minDR", &b_minDR );
  bdt_XL->AddVariable( "XlepPt", &b_XlepPt );
  bdt_XL->AddVariable( "mT2", &b_mT2 );
  bdt_XL->AddVariable( "mT", &b_mT );
  bdt_XL->AddVariable( "DiJetM12", &b_DiJetM12 );
  bdt_XL->BookMVA("BDT", weightXL); 
  
  bdt_FH = new TMVA::Reader();
  bdt_FH->AddVariable( "Met", &b_Met );
  bdt_FH->AddVariable( "all_muEtaDiff", &b_all_muEtaDiff );
  bdt_FH->AddVariable( "all_muPtDiff", &b_all_muPhiDiff );
  bdt_FH->AddVariable( "all_muPhiDiff", &b_all_muPhiDiff );
  bdt_FH->AddVariable( "all_muDR", &b_all_muDR );
  bdt_FH->AddVariable( "all_Dilep_Pt", &b_all_Dilep_Pt );
  bdt_FH->AddVariable( "all_Dilep_Eta", &b_all_Dilep_Eta );
  bdt_FH->AddVariable( "nnonbjet", &b_nnonbjet );
  bdt_FH->AddVariable( "nbjet", &b_nbjet );
  bdt_FH->AddVariable( "Central_Jets", &b_Central_Jets );
  bdt_FH->AddVariable( "Forward_Jets", &b_Forward_Jets );
  bdt_FH->AddVariable( "minDR1", &b_minDR1 );
  bdt_FH->AddVariable( "minDR2", &b_minDR2 );
  bdt_FH->AddVariable( "mT2", &b_mT2 );
  bdt_FH->AddVariable( "mT", &b_mT );
  bdt_FH->AddVariable( "DiJetM12", &b_DiJetM12 );
  bdt_FH->AddVariable( "DiJetM13", &b_DiJetM13 );
  bdt_FH->AddVariable( "DiJetM14", &b_DiJetM14 );
  bdt_FH->AddVariable( "DiJetM23", &b_DiJetM23 );
  bdt_FH->AddVariable( "DiJetM24", &b_DiJetM24 );
  bdt_FH->AddVariable( "DiJetM34", &b_DiJetM34 );
  bdt_FH->BookMVA("BDT", weightFH); 
  
  bdt_Out = new TMVA::Reader();
  bdt_Out->AddVariable( "Met", &b_Met );
  bdt_Out->AddVariable( "all_muEtaDiff", &b_all_muEtaDiff );
  bdt_Out->AddVariable( "all_muPtDiff", &b_all_muPhiDiff );
  bdt_Out->AddVariable( "all_muPhiDiff", &b_all_muPhiDiff );
  bdt_Out->AddVariable( "all_muDR", &b_all_muDR );
  bdt_Out->AddVariable( "all_Dilep_Pt", &b_all_Dilep_Pt );
  bdt_Out->AddVariable( "all_Dilep_Eta", &b_all_Dilep_Eta );
  bdt_Out->AddVariable( "nnonbjet", &b_nnonbjet );
  bdt_Out->AddVariable( "nbjet", &b_nbjet );
  bdt_Out->AddVariable( "Central_Jets", &b_Central_Jets );
  bdt_Out->AddVariable( "Forward_Jets", &b_Forward_Jets );
  bdt_Out->AddVariable( "mT2", &b_mT2 );
  bdt_Out->AddVariable( "mT", &b_mT );
  bdt_Out->AddVariable( "DiJetM12", &b_DiJetM12 );
  bdt_Out->BookMVA("BDT", weightOut);

  bdt_noB = new TMVA::Reader();
  bdt_noB->AddVariable( "Met", &b_Met );
//  bdt_noB->AddVariable( "all_muEtaDiff", &b_all_muEtaDiff );
//  bdt_noB->AddVariable( "all_muPtDiff", &b_all_muPhiDiff );
//  bdt_noB->AddVariable( "all_muPhiDiff", &b_all_muPhiDiff );
//  bdt_noB->AddVariable( "all_muDR", &b_all_muDR );
  bdt_noB->AddVariable( "all_Dilep_Pt", &b_all_Dilep_Pt );
  bdt_noB->AddVariable( "all_Dilep_Eta", &b_all_Dilep_Eta );
 // bdt_noB->AddVariable( "nelec", &b_nelec );
 // bdt_noB->AddVariable( "nexLep", &b_nexLep );
 // bdt_noB->AddVariable( "nmuon", &b_nmuon );
  bdt_noB->AddVariable( "Central_Jets", &b_Central_Jets );
  bdt_noB->AddVariable( "Forward_Jets", &b_Forward_Jets );
  bdt_noB->AddVariable( "etaJ1", &b_etaJ1 );
  bdt_noB->AddVariable( "etaJ2", &b_etaJ2 );
  bdt_noB->AddVariable( "DijetM1", &b_DijetM1 );
  bdt_noB->AddVariable( "DijetM2", &b_DijetM2 );
  bdt_noB->AddVariable( "DijetEta1", &b_DijetEta1 );
  bdt_noB->AddVariable( "DijetEta2", &b_DijetEta2 );
  bdt_noB->BookMVA("BDT", weightnoB);

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

public :
  Bool_t m_isMC, m_isDL, m_isSL_e, m_isSL_m;

  //set output file
  void setOutput(std::string outputName);
  void LoadModules(pileUpTool* pileUp, lumiTool* lumi);
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

