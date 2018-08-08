//#define nanoAnalyser_cxx
#define Events_cxx
#include "DataFormats/HepMCCandidate/interface/GenStatusFlags.h"
//#include "nano/analysis/interface/nanoBase.h"
#include "nano/analysis/interface/topEventSelectionSL.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMatrixDSymfwd.h>
#include <TMatrixDSymEigen.h>
#include <iostream>
#include <fstream>
#include <cstdlib>


using namespace std;
/*
To compile:
cd $CMSSW_BASE ; scram b -j 8 ; cd -

To run:
singletopAnalyser [LIST_FILE_OF_FILES] ["MC" or "RD"] [idx of the first root file] [(idx+1) of the last root file]

To run (for the author):
singletopAnalyser -q [LIST_FILE_OF_FILES] ["MC" or "RD"] [idx of the first root file] [(idx+1) of the last root file]

To throw jobs to condor:
python makejobs.py [DEST DIR] j`python batch_nanoAOD.py n [# OF FILES PER JOB]`
*/


class singletopAnalyser : public topEventSelectionSL {
private: 
  Long64_t m_nEventIdx;
  
  TH1D *m_h1CosS;
  TH1D *m_h1CosSReco;
  TH1D *m_h1IsHighPtMuon, *m_h1DRMuon;
  
  Bool_t m_isSig;
  int m_nBkgType;
  Bool_t m_isFullGen;
  
  Int_t m_nIdxGenTop, m_nIdxGenAssoQ, m_nIdxGenLep, m_nIdxGenBFromT, m_nIdxGenW;
  Int_t m_nIDGenLep;
  
  Int_t m_nIdxGenNeu, m_nIDGenNeu;
  
  Int_t b_onlyGen;
  
  Float_t b_met_sumEt;
  
  Int_t b_filter_met;
  
  TLorentzVector b_gentop1;
  TLorentzVector b_genW;
  
  Int_t m_nIdxGenTop2nd;
  Int_t m_nIdxGenW2nd, m_nIdxGenBFromT2nd, m_nIdxGenLep2nd;
  Int_t m_nIDGenLep2nd;
  
  TLorentzVector b_ttbar_gentop2;
  TLorentzVector b_ttbar_genW2;
  
  TLorentzVector b_ttbar_genLep1;
  TLorentzVector b_ttbar_genLep2;
  Int_t b_ttbar_lep1_pdgId, b_ttbar_lep2_pdgId;
  
  TLorentzVector b_ttbar_genB1;
  TLorentzVector b_ttbar_genB2;
  
  Float_t b_ttbar_deltaR_B1_B2;
  
  Int_t b_ttbar_channel;
  Int_t b_ttbar_missedlep;
  
  TLorentzVector b_ttbar_lep1;
  TLorentzVector b_ttbar_lep2;
  Float_t b_ttbar_lep1_deltaRGen, b_ttbar_lep2_deltaRGen;
  Float_t b_ttbar_lep1_reliso, b_ttbar_lep2_reliso;
  
  TLorentzVector b_lep2NonIso;
  Int_t b_lep2NonIso_pid;
  Float_t b_lep2NonIso_reliso;
  
  TLorentzVector b_jet1, b_bjet1, b_bjet2;
  Float_t b_CSVjet1, b_CSVbjet1, b_CSVbjet2;
  
  Float_t b_jetQ1, b_jetQ2, b_jetQ3;
  Float_t b_jetC;
  
  TLorentzVector b_DiffLepMom11;
  TLorentzVector b_DiffLepMom12;
  TLorentzVector b_DiffLepMom21;
  TLorentzVector b_DiffLepMom22;
  
  TLorentzVector b_W1;
  TLorentzVector b_top1;
  TLorentzVector b_top1_lower;
  TLorentzVector b_top1_higher;
  TLorentzVector b_top1_imaginary;
  TLorentzVector b_top1_genneu;
  TLorentzVector b_top1_genWMass;
  TLorentzVector b_top1_genW;
  TLorentzVector b_top1_genWgenB;
  
  TLorentzVector b_recoW_lower;
  TLorentzVector b_recoW_higher;
  TLorentzVector b_recoNeu_lower;
  TLorentzVector b_recoNeu_higher;
  
  Int_t b_truthtop1_lowhighOn;
  Int_t b_truthtop1_imaginaryOn;
  Float_t b_truthtop1_lowhighNeuPT;
  
  TLorentzVector b_fatjet;
  Int_t b_nfatjet;
  Int_t b_nsubjet;
  
  Int_t b_nhadcand;
  
  Float_t b_cos_star_gen;
  Float_t b_cos_star_bjet_gen;
  Float_t b_cos_labf_gen;
  Float_t b_cos_labf_bjet_gen;
  
  Float_t b_cos_star;
  
public:
  enum {
    MY_FLAG_RECOW_WRESTRICTION, 
    MY_FLAG_RECOW_WRESTR_HIGH, 
    MY_FLAG_RECOW_WRESTR_IMG, 
    MY_FLAG_RECOW_WMASS_IN_GEN, 
    MY_FLAG_RECOW_NEU_IN_GEN, 
    MY_FLAG_RECOW_W_IN_GEN
  };
  
  enum {
    MY_FLAG_BKGTYPE_TTBAR = 1, 
    MY_FLAG_BKGTYPE_QCD, 
    MY_FLAG_BKGTYPE_WJETS, 
    MY_FLAG_BKGTYPE_ST_OTHERS
  };
  
  enum {
    // # of final lepton = 0
    MY_TYPE_TTBAR_FULLHADRONIC = 1, // Full-hadronic
    MY_TYPE_TTBAR_MONOTAU_HADRONIC, // Semileptonic
    MY_TYPE_TTBAR_DITAU_HADRONIC, // Dileptonic
    // # of final lepton = 1
    MY_TYPE_TTBAR_SEMILEPTONIC, // Semileptonic
    MY_TYPE_TTBAR_ONETAU_LEPTONIC, 
    MY_TYPE_TTBAR_ONELEPTON_ONETAU_HADRONIC, // Dileptonic
    MY_TYPE_TTBAR_DITAU_LEPTONIC_HADRONIC, 
    // # of final lepton = 2 (Only from Dileptonic)
    MY_TYPE_TTBAR_DILEPTONIC, 
    MY_TYPE_TTBAR_ONELEPTON_ONETAU_LEPTONIC, 
    MY_TYPE_TTBAR_DITAU_DILEPTONIC, 
  };
  
public: 
  //set output file
  void setOutput(std::string outputName);
  void MakeBranch(TTree* t);
  void resetBranch();
  void LoadModules(pileUpTool* pileUp, lumiTool* lumi);
  //void collectTMVAvalues();
  singletopAnalyser(TTree *tree=0, TTree *had=0, TTree *hadtruth=0, 
      Bool_t isMC = false, Bool_t isSig = false, int nBkgType = 0, Bool_t isFullGen = false);
  ~singletopAnalyser();
  virtual void     Loop();
  
  int DoMoreGenLvl(TLorentzVector &vec4Top, TLorentzVector &vec4Lep, TLorentzVector &vec4AssoQ);
  int RunEvt();
  
  TLorentzVector Get4VecGen(int nIdx);
  
  int GetIdxGenTop();
  int GetIdxGenAssoQuark();
  int GetIdxGenLepton(int nIsFor2nd = 0);
  int GetIdxGenBFromTop(int nIsFor2nd = 0);
  int GetIdxGenWFromTop(int nIsFor2nd = 0);
  int GetIdxGenNeutrino();
  
  int GetGenInfoSignal();
  
  int GetIdxGenTop2nd();
  int GetLeptonInfoMatching();
  
  int GetGenInfoTTbar();
  int GetGenInfoQCD();
  int GetGenInfoWJets();
  int GetGenInfoSTOthers();
  
  int GetJets();
  int CalcSphericity();
  
  TLorentzVector RecoWFromTop(double *pdDiffMET, int nFlag = MY_FLAG_RECOW_WRESTRICTION);
  int RecoTop();
  int CalcRecoCosStar();
};


singletopAnalyser::singletopAnalyser(TTree *tree, TTree *had, TTree *hadtruth, 
  Bool_t isMC, Bool_t isSig, int nBkgType, Bool_t isFullGen) : 
  topEventSelectionSL(tree, had, hadtruth, isMC), 
  m_isSig(isSig), m_nBkgType(nBkgType), m_isFullGen(isFullGen)
{
  m_output = NULL;
  //m_btagSF.initCSVWeight();
  std::string env = getenv("CMSSW_BASE");
  m_rocCor = new RoccoR(env + "/src/nano/analysis/data/rcdata.2016.v3/");
}


singletopAnalyser::~singletopAnalyser() {
  if ( m_output != NULL ) {
    m_output->Write();
    m_output->Close();
  }
}


void singletopAnalyser::setOutput(std::string outputName) {
  m_output = TFile::Open(outputName.c_str(), "recreate");
  
  m_tree = new TTree("event", "event");
  MakeBranch(m_tree);
  
  h_nevents = new TH1D("nevents", "nevents", 2, -1, 1);
  h_weights = new TH1D("weight", "weight", 2, -1, 1);
  h_genweights = new TH1D("genweight", "genweight", 2, -1, 1);
  
  h_cutFlow = new TH1D("cutFlow", "cutFlow", 7, 0, 7);
  h_cutFlowEl = new TH1D("cutFlowEl", "cutFlowEl", 7, 0, 7);
  h_cutFlowMu = new TH1D("cutFlowMu", "cutFlowMu", 7, 0, 7);
  
  m_h1CosS = new TH1D("cos_star", "cos_star", 10, -1, 1);
  m_h1CosSReco = new TH1D("cos_star_reco", "cos_star_reco", 10, -1, 1);
  m_h1IsHighPtMuon = new TH1D("isClosest", "isClosest", 2, 0, 2);
  m_h1DRMuon = new TH1D("deltaR_muon", "deltaR_muon", 50, 0, 5);
}


void singletopAnalyser::MakeBranch(TTree *t) {
  t->Branch("onlyGen", &b_onlyGen, "onlyGen/I");
  
  t->Branch("step", &b_step, "step/I");
  t->Branch("channel", &b_channel, "channel/I");
  t->Branch("njet", &b_njet, "njet/I");
  t->Branch("nbjet", &b_nbjet, "nbjet/I");
  
  t->Branch("met_sumEt", &b_met_sumEt, "met_sumEt/I");
  
  t->Branch("filter_met", &b_filter_met, "filter_met/I");
  
  t->Branch("gentop1", "TLorentzVector", &b_gentop1);
  t->Branch("genW", "TLorentzVector", &b_genW);
  
  t->Branch("ttbar_gentop2", "TLorentzVector", &b_ttbar_gentop2);
  t->Branch("ttbar_genW2", "TLorentzVector", &b_ttbar_genW2);
  
  t->Branch("ttbar_genLep1", "TLorentzVector", &b_ttbar_genLep1);
  t->Branch("ttbar_genLep2", "TLorentzVector", &b_ttbar_genLep2);
  t->Branch("ttbar_lep1_pdgId", &b_ttbar_lep1_pdgId, "ttbar_lep1_pdgId/I");
  t->Branch("ttbar_lep2_pdgId", &b_ttbar_lep2_pdgId, "ttbar_lep2_pdgId/I");
  
  t->Branch("ttbar_genB1", "TLorentzVector", &b_ttbar_genB1);
  t->Branch("ttbar_genB2", "TLorentzVector", &b_ttbar_genB2);
  
  t->Branch("ttbar_deltaR_B1_B2", &b_ttbar_deltaR_B1_B2, "ttbar_deltaR_B1_B2/F");
  
  t->Branch("ttbar_channel", &b_ttbar_channel, "ttbar_channel/I");
  t->Branch("ttbar_missedlep", &b_ttbar_missedlep, "ttbar_missedlep/I");
  
  t->Branch("ttbar_lep1", "TLorentzVector", &b_ttbar_lep1);
  t->Branch("ttbar_lep2", "TLorentzVector", &b_ttbar_lep2);
  t->Branch("ttbar_lep1_deltaRGen", &b_ttbar_lep1_deltaRGen, "ttbar_lep1_deltaRGen/F");
  t->Branch("ttbar_lep2_deltaRGen", &b_ttbar_lep2_deltaRGen, "ttbar_lep2_deltaRGen/F");
  t->Branch("ttbar_lep1_reliso", &b_ttbar_lep1_reliso, "ttbar_lep1_reliso/F");
  t->Branch("ttbar_lep2_reliso", &b_ttbar_lep2_reliso, "ttbar_lep2_reliso/F");
  
  t->Branch("lep", "TLorentzVector", &b_lep);
  t->Branch("jet1", "TLorentzVector", &b_jet1);
  t->Branch("bjet1", "TLorentzVector", &b_bjet1);
  t->Branch("bjet2", "TLorentzVector", &b_bjet2);
  
  t->Branch("isolep", &b_isolep, "isolep/F");
  
  t->Branch("CSVv2_jet1", &b_CSVjet1, "CSVv2_jet1/F");
  t->Branch("CSVv2_bjet1", &b_CSVbjet1, "CSVv2_bjet1/F");
  t->Branch("CSVv2_bjet2", &b_CSVbjet2, "CSVv2_bjet2/F");
  
  t->Branch("DiffLepMom11", "TLorentzVector", &b_DiffLepMom11);
  t->Branch("DiffLepMom12", "TLorentzVector", &b_DiffLepMom12);
  t->Branch("DiffLepMom21", "TLorentzVector", &b_DiffLepMom21);
  t->Branch("DiffLepMom22", "TLorentzVector", &b_DiffLepMom22);
  
  t->Branch("lep_pid", &b_lep_pid, "lep_pid/I");
  
  t->Branch("met", &b_met, "met/F");
  t->Branch("met_phi", &b_met_phi, "met_phi/F");
  
  t->Branch("lep2NonIso", "TLorentzVector", &b_lep2NonIso);
  t->Branch("lep2NonIso_pid", &b_lep2NonIso_pid, "lep2NonIso_pid/I");
  t->Branch("lep2NonIso_reliso", &b_lep2NonIso_reliso, "lep2NonIso_reliso/F");
  
  t->Branch("W1", "TLorentzVector", &b_W1);
  t->Branch("top1", "TLorentzVector", &b_top1);
  t->Branch("top1_lower", "TLorentzVector", &b_top1_lower);
  t->Branch("top1_higher", "TLorentzVector", &b_top1_higher);
  t->Branch("top1_imaginary", "TLorentzVector", &b_top1_imaginary);
  t->Branch("top1_genneu", "TLorentzVector", &b_top1_genneu);
  t->Branch("top1_genWMass", "TLorentzVector", &b_top1_genWMass);
  t->Branch("top1_genW", "TLorentzVector", &b_top1_genW);
  t->Branch("top1_genWgenB", "TLorentzVector", &b_top1_genWgenB);
  
  t->Branch("recoW_lower",  "TLorentzVector", &b_recoW_lower);
  t->Branch("recoW_higher", "TLorentzVector", &b_recoW_higher);
  t->Branch("recoNeu_lower",  "TLorentzVector", &b_recoNeu_lower);
  t->Branch("recoNeu_higher", "TLorentzVector", &b_recoNeu_higher);
  
  t->Branch("truthtop1_lowhighOn", &b_truthtop1_lowhighOn, "truthtop1_lowhighOn/I");
  t->Branch("truthtop1_imaginaryOn", &b_truthtop1_imaginaryOn, "truthtop1_imaginaryOn/I");
  t->Branch("truthtop1_lowhighNeuPT", &b_truthtop1_lowhighNeuPT, "truthtop1_lowhighNeuPT/F");
  
  t->Branch("jetQ1", &b_jetQ1, "jetQ1/F");
  t->Branch("jetQ2", &b_jetQ2, "jetQ2/F");
  t->Branch("jetQ3", &b_jetQ3, "jetQ3/F");
  t->Branch("jetC", &b_jetC, "jetC/F");
  
  t->Branch("fatjet", "TLorentzVector", &b_fatjet);
  t->Branch("nfatjet", &b_nfatjet, "nfatjet/I");
  t->Branch("nsubjet", &b_nsubjet, "nsubjet/I");
  
  t->Branch("nhadcand", &b_nhadcand, "nhadcand/I");
  
  t->Branch("trig_e", &b_trig_e, "trig_e/O");
  t->Branch("trig_m", &b_trig_m, "trig_m/O");
  
  t->Branch("tri", &b_tri, "tri/F");
  t->Branch("tri_up", &b_tri_up, "tri_up/F");
  t->Branch("tri_dn", &b_tri_dn, "tri_dn/F");
  
  t->Branch("npvs", &b_nvertex, "npvs/I");
  
  t->Branch("weight", &b_weight, "weight/F");
  t->Branch("puweight", &b_puweight, "puweight/F");
  t->Branch("genweight", &b_genweight, "genweight/F");
  t->Branch("btagweight", &b_btagweight, "btagweight/F");
  
  t->Branch("mueffweight", &b_mueffweight, "mueffweight/F");
  t->Branch("mueffweight_up", &b_mueffweight_up, "mueffweight_up/F");
  t->Branch("mueffweight_dn", &b_mueffweight_dn, "mueffweight_dn/F");
  
  t->Branch("eleffweight", &b_eleffweight, "eleffweight/F");
  t->Branch("eleffweight_up", &b_eleffweight_up, "eleffweight_up/F");
  t->Branch("eleffweight_dn", &b_eleffweight_dn, "eleffweight_dn/F");
  
  t->Branch("cos_star_gen",      &b_cos_star_gen,      "cos_star_gen/F");
  t->Branch("cos_star_bjet_gen", &b_cos_star_bjet_gen, "cos_star_bjet_gen/F");
  t->Branch("cos_labf_gen",      &b_cos_labf_gen,      "cos_labf_gen/F");
  t->Branch("cos_labf_bjet_gen", &b_cos_labf_bjet_gen, "cos_labf_bjet_gen/F");
  
  t->Branch("cos_star", &b_cos_star, "cos_star/F");
}


void singletopAnalyser::resetBranch() {
  topEventSelectionSL::Reset();
  
  b_onlyGen = 0;
  
  b_met_sumEt = -9;
  
  b_filter_met = 0;
  
  b_gentop1.SetPtEtaPhiM(0, 0, 0, 0);
  b_genW.SetPtEtaPhiM(0, 0, 0, 0);
  
  b_ttbar_gentop2.SetPtEtaPhiM(0, 0, 0, 0);
  b_ttbar_genW2.SetPtEtaPhiM(0, 0, 0, 0);
  
  b_ttbar_genLep1.SetPtEtaPhiM(0, 0, 0, 0);
  b_ttbar_genLep2.SetPtEtaPhiM(0, 0, 0, 0);
  b_ttbar_lep1_pdgId = b_ttbar_lep2_pdgId = 0;
  
  b_ttbar_genB1.SetPtEtaPhiM(0, 0, 0, 0);
  b_ttbar_genB2.SetPtEtaPhiM(0, 0, 0, 0);
  
  b_ttbar_deltaR_B1_B2 = 0;
  
  b_ttbar_channel = 0;
  b_ttbar_missedlep = 0;
  
  b_ttbar_lep1.SetPtEtaPhiM(0, 0, 0, 0);
  b_ttbar_lep2.SetPtEtaPhiM(0, 0, 0, 0);
  b_ttbar_lep1_deltaRGen = b_ttbar_lep2_deltaRGen = 0;
  b_ttbar_lep1_reliso = b_ttbar_lep2_reliso = 0;
  
  b_jet1.SetPtEtaPhiM(0, 0, 0, 0);
  b_bjet1.SetPtEtaPhiM(0, 0, 0, 0);
  b_bjet2.SetPtEtaPhiM(0, 0, 0, 0);
  
  b_CSVjet1 = -1.0;
  b_CSVbjet1 = -1.0;
  b_CSVbjet2 = -1.0;
  
  b_jetQ1 = b_jetQ2 = b_jetQ3 = 0;
  b_jetC = 0;
  
  b_lep2NonIso.SetPtEtaPhiM(0, 0, 0, 0);
  b_lep2NonIso_pid = 0;
  b_lep2NonIso_reliso = 0;
  
  b_DiffLepMom11.SetPtEtaPhiM(0, 0, 0, 0);
  b_DiffLepMom12.SetPtEtaPhiM(0, 0, 0, 0);
  b_DiffLepMom21.SetPtEtaPhiM(0, 0, 0, 0);
  b_DiffLepMom22.SetPtEtaPhiM(0, 0, 0, 0);
  
  b_W1.SetPtEtaPhiM(0, 0, 0, 0);
  b_top1.SetPtEtaPhiM(0, 0, 0, 0);
  b_top1_lower.SetPtEtaPhiM(0, 0, 0, 0);
  b_top1_higher.SetPtEtaPhiM(0, 0, 0, 0);
  b_top1_imaginary.SetPtEtaPhiM(0, 0, 0, 0);
  b_top1_genneu.SetPtEtaPhiM(0, 0, 0, 0);
  b_top1_genWMass.SetPtEtaPhiM(0, 0, 0, 0);
  b_top1_genW.SetPtEtaPhiM(0, 0, 0, 0);
  b_top1_genWgenB.SetPtEtaPhiM(0, 0, 0, 0);
  
  b_recoW_lower.SetPtEtaPhiM(0, 0, 0, 0);
  b_recoW_higher.SetPtEtaPhiM(0, 0, 0, 0);
  b_recoNeu_lower.SetPtEtaPhiM(0, 0, 0, 0);
  b_recoNeu_higher.SetPtEtaPhiM(0, 0, 0, 0);
  
  b_truthtop1_lowhighOn = 0;
  b_truthtop1_imaginaryOn = 0;
  b_truthtop1_lowhighNeuPT = 0.0;
  
  b_fatjet.SetPtEtaPhiM(0, 0, 0, 0);
  b_nfatjet = -1;
  b_nsubjet = -1;
  
  b_nhadcand = 0;
  
  b_cos_star_gen = b_cos_star_bjet_gen = -2;
  b_cos_labf_gen = b_cos_labf_bjet_gen = -2;
  b_cos_star = -2;
}


////////////////////////////////////////////////////////////////////////////////
// For Gen level test : Start
////////////////////////////////////////////////////////////////////////////////


int singletopAnalyser::GetIdxGenTop() {
  UInt_t i;
  
  // Finding the top quark
  for ( i = 0 ; i < nGenPart ; i++ ) {
    if ( GenPart_pdgId[ i ] == 6 || GenPart_pdgId[ i ] == -6 ) {
      if ( ( GenPart_statusFlags[ i ] & ( 1 << reco::GenStatusFlags::kIsFirstCopy ) ) == 0 ) continue;
      
      m_nIdxGenTop = i;
      break;
    }
  }
  
  return ( i < nGenPart ? 0 : -1 );
}


int singletopAnalyser::GetIdxGenAssoQuark() {
  UInt_t i;
  
  for ( i = 0 ; i < nGenPart ; i++ ) {
    int nPID = abs(GenPart_pdgId[ i ]);
    
    if ( 0 < nPID && nPID < 5 && GenPart_genPartIdxMother[ i ] == 0 ) {
      m_nIdxGenAssoQ = i;
      break;
    }
  }
  
  return ( i < nGenPart ? 0 : -1 );
}


int singletopAnalyser::GetIdxGenBFromTop(int nIsFor2nd) {
  UInt_t i;
  Int_t j;
  
  Int_t nIdxTop = ( nIsFor2nd == 0 ? m_nIdxGenTop : m_nIdxGenTop2nd );
  
  // Finding the b quark from the top quark
  for ( i = 0 ; i < nGenPart ; i++ ) {
    int nPID = GenPart_pdgId[ i ];
    if ( abs(nPID) != 5 ) continue;
    
    int nIdxMother = GenPart_genPartIdxMother[ i ];
    if ( nIdxMother < 0 ) continue;
    
    if ( GenPart_pdgId[ nIdxMother ] * nPID != 6 * 5 ) continue;
    
    // Checking whether the mother top quark is same as what we are looking at (nIdxTop)
    for ( j = nIdxMother ; j >= 0 ; j = GenPart_genPartIdxMother[ j ] ) {
      if ( j == nIdxTop ) break;
    }
    
    if ( j < 0 ) continue;
    
    m_nIdxGenBFromT = i;
    break;
  }
  
  return ( i < nGenPart ? 0 : -1 );
}


int singletopAnalyser::GetIdxGenWFromTop(int nIsFor2nd) {
  UInt_t i;
  Int_t j;
  
  Int_t nIdxTop = ( nIsFor2nd == 0 ? m_nIdxGenTop : m_nIdxGenTop2nd );
  
  // Finding the W from the top quark
  for ( i = 0 ; i < nGenPart ; i++ ) {
    int nPID = GenPart_pdgId[ i ];
    if ( abs(nPID) != 24 ) continue;
    
    int nIdxMother = GenPart_genPartIdxMother[ i ];
    if ( nIdxMother < 0 ) continue;
    
    if ( GenPart_pdgId[ nIdxMother ] * nPID != 6 * 24 ) continue;
    
    // Checking whether the mother top quark is same as what we are looking at (nIdxTop)
    for ( j = nIdxMother ; j >= 0 ; j = GenPart_genPartIdxMother[ j ] ) {
      if ( j == nIdxTop ) break;
    }
    
    if ( j < 0 ) continue;
    
    m_nIdxGenW = i;
    break;
  }
  
  return ( i < nGenPart ? 0 : -1 );
}


int singletopAnalyser::GetIdxGenLepton(int nIsFor2nd) {
  UInt_t i;
  Int_t j;
  
  Int_t nIdxW = ( nIsFor2nd == 0 ? m_nIdxGenW : m_nIdxGenW2nd );
  
  // Finding the lepton from the top quark
  for ( i = 0 ; i < nGenPart ; i++ ) {
    int nPID = GenPart_pdgId[ i ];
    if ( abs(nPID) != 11 && abs(nPID) != 13 ) continue;
    
    int nIdxMother = GenPart_genPartIdxMother[ i ];
    if ( nIdxMother < 0 ) continue;
    
    if ( abs(nPID) == 11 && GenPart_pdgId[ nIdxMother ] * nPID != 24 * -11 ) continue;
    if ( abs(nPID) == 13 && GenPart_pdgId[ nIdxMother ] * nPID != 24 * -13 ) continue;
    
    // Checking whether the mother top quark is same as what we are looking at (nIdxTop)
    for ( j = nIdxMother ; j >= 0 ; j = GenPart_genPartIdxMother[ j ] ) {
      if ( j == nIdxW ) break;
    }
    
    if ( j < 0 ) continue;
    
    m_nIdxGenLep = i;
    m_nIDGenLep = GenPart_pdgId[ i ];
    break;
  }
  
  return ( i < nGenPart ? 0 : -1 );
}


int singletopAnalyser::GetIdxGenNeutrino() {
  UInt_t i;
  Int_t j;
  
  // Finding the lepton from the top quark
  for ( i = 0 ; i < nGenPart ; i++ ) {
    int nPID = abs(GenPart_pdgId[ i ]);
    int nPassW = 0;
    
    if ( ( nPID == 12 || nPID == 14 ) ) {
      for ( j = GenPart_genPartIdxMother[ i ] ;; ) {
        int nParentID = abs(GenPart_pdgId[ j ]);
        
        if ( nParentID == 24 ) {
          int nGrandParentID = abs(GenPart_pdgId[ GenPart_genPartIdxMother[ j ] ]);
          if ( nGrandParentID == 6 ) nPassW = 1;
        }
        
        if ( nPassW == 1 && j == m_nIdxGenTop ) break;
        
        j = GenPart_genPartIdxMother[ j ];
        if ( j < 0 ) break;
      }
      
      if ( j < 0 ) continue;
      
      m_nIdxGenNeu = i;
      m_nIDGenNeu = GenPart_pdgId[ i ];
      break;
    }
  }
  
  return ( i < nGenPart ? 0 : -1 );
}


int singletopAnalyser::DoMoreGenLvl(TLorentzVector &vec4Top, TLorentzVector &vec4Lep, TLorentzVector &vec4AssoQ) {
  if ( GetIdxGenNeutrino() != 0 ) return 1;
  
  
  
  return 0;
}


TLorentzVector singletopAnalyser::Get4VecGen(int nIdx) {
  TLorentzVector vec4Res;
  
  vec4Res.SetPtEtaPhiM(GenPart_pt[ nIdx ], GenPart_eta[ nIdx ], 
    GenPart_phi[ nIdx ], GenPart_mass[ nIdx ]);
  
  return vec4Res;
}


int singletopAnalyser::GetGenInfoSignal() {
  // If any one of them holds, then it is reaaaaally strange!
  if ( GetIdxGenTop() != 0 ) return -1;
  if ( GetIdxGenAssoQuark() != 0 ) return -1;
  if ( GetIdxGenWFromTop() != 0 ) return -1;
  if ( GetIdxGenLepton() != 0 ) return -1;
  if ( GetIdxGenBFromTop() != 0 ) return -1;
  if ( GetIdxGenNeutrino() != 0 ) return -1;
  
  // Getting the 4-momentum of top quark IN GEN LEVEL; it will be used to back-boost
  // vec4TopRest.Boost(-vec4Top.BoostVector()); // vec4TopRest: Copy of vec4Top, it makes back-boost
  TLorentzVector vec4Top, vec4AssoQ, vec4B, vec4W, vec4Lep;
  TLorentzVector vec4AssoQRest, vec4LepRest, vec4BRest;
  
  vec4Top = Get4VecGen(m_nIdxGenTop);
  vec4AssoQ = Get4VecGen(m_nIdxGenAssoQ);
  vec4B   = Get4VecGen(m_nIdxGenBFromT);
  vec4W   = Get4VecGen(m_nIdxGenW);
  vec4Lep = Get4VecGen(m_nIdxGenLep);
  
  b_gentop1 = vec4Top;
  b_genW = vec4W;
  
  // Calculating theta*s
  
  vec4AssoQRest = vec4AssoQ;
  vec4LepRest = vec4Lep;
  vec4BRest   = vec4B;
  
  vec4AssoQRest.Boost(-vec4Top.BoostVector());
  vec4LepRest.Boost(-vec4Top.BoostVector());
  vec4BRest.Boost(-vec4Top.BoostVector());
  
  TVector3 vec3AssoQ;
  vec3AssoQ = vec4AssoQRest.Vect().Unit();
  
  b_cos_star_gen      = vec4LepRest.Vect().Unit().Dot(vec3AssoQ);
  b_cos_star_bjet_gen = vec4BRest.Vect().Unit().Dot(vec3AssoQ);
  m_h1CosS->Fill(b_cos_star_gen);
  
  // Calculating theta in lab frame
  
  vec3AssoQ = vec4AssoQ.Vect().Unit();
  b_cos_labf_gen      = vec4Lep.Vect().Unit().Dot(vec3AssoQ);
  b_cos_labf_bjet_gen = vec4B.Vect().Unit().Dot(vec3AssoQ);
  
  return 0;
}


// Background study: TTbar


int singletopAnalyser::GetIdxGenTop2nd() {
  UInt_t i;
  
  // Finding the top quark
  for ( i = m_nIdxGenTop + 1 ; i < nGenPart ; i++ ) {
    if ( GenPart_pdgId[ i ] == 6 || GenPart_pdgId[ i ] == -6 ) {
      if ( ( GenPart_statusFlags[ i ] & ( 1 << reco::GenStatusFlags::kIsFirstCopy ) ) == 0 ) continue;
      
      m_nIdxGenTop2nd = i;
      break;
    }
  }
  
  return ( i < nGenPart ? 0 : -1 );
}

 
int singletopAnalyser::GetGenInfoTTbar() {
  Int_t nIdxBuf;
  
  int nRes;
  
  if ( GetIdxGenTop() != 0 ) return 0;
  if ( GetIdxGenTop2nd() != 0 ) return 0;
  
  // Ordering the top quarks by ID
  if ( GenPart_pdgId[ m_nIdxGenTop ] < 0 ) {
    nIdxBuf = m_nIdxGenTop;
    m_nIdxGenTop = m_nIdxGenTop2nd;
    m_nIdxGenTop2nd = nIdxBuf;
  }
  
  // Getting the indices of b quarks from top quarks
  if ( GetIdxGenBFromTop() != 0 ) return 0;
  
  // Getting the 2nd b quark (NOTE: the idx of the found b quark is kept in m_nIdxGenBFromT)
  nIdxBuf = m_nIdxGenBFromT;
  if ( GetIdxGenBFromTop(1) != 0 ) return 0;
  m_nIdxGenBFromT2nd = m_nIdxGenBFromT;
  m_nIdxGenBFromT = nIdxBuf;
  
  if ( GetIdxGenWFromTop() != 0 ) return 0;
  
  nIdxBuf = m_nIdxGenW;
  if ( GetIdxGenWFromTop(1) != 0 ) return 0;
  m_nIdxGenW2nd = m_nIdxGenW;
  m_nIdxGenW = nIdxBuf;
  
  b_ttbar_channel = 1;
  
  // Getting the indices of leptops from (W from) top quarks
  nRes = GetIdxGenLepton();
  
  if ( nRes == 0 ) {
    b_ttbar_channel++;
    
    b_ttbar_lep1_pdgId = m_nIDGenLep;
    nIdxBuf = m_nIdxGenLep;
    
    // Getting the 2nd lepton
    nRes = GetIdxGenLepton(1);
    
    // If existing (NOTE: the idx of the found lepton is kept in m_nIdxGenLep)
    if ( nRes == 0 ) {
      m_nIdxGenLep2nd = m_nIdxGenLep;
      m_nIdxGenLep = nIdxBuf;
      
      m_nIDGenLep2nd = m_nIDGenLep;
      
      b_ttbar_channel++;
    }
  } else {
    nRes = GetIdxGenLepton(1);
    if ( nRes == 0 ) b_ttbar_channel++;
  }
  
  b_gentop1 = Get4VecGen(m_nIdxGenTop);
  b_ttbar_gentop2 = Get4VecGen(m_nIdxGenTop2nd);
  b_ttbar_genB1 = Get4VecGen(m_nIdxGenBFromT);
  b_ttbar_genB2 = Get4VecGen(m_nIdxGenBFromT2nd);
  
  if ( b_ttbar_channel >= 2 ) b_ttbar_genLep1 = Get4VecGen(m_nIdxGenLep);
  if ( b_ttbar_channel >= 3 ) b_ttbar_genLep2 = Get4VecGen(m_nIdxGenLep2nd);
  
  // Calculating further variables
  
  b_ttbar_deltaR_B1_B2 = b_ttbar_genB2.DeltaR(b_ttbar_genB1);
  
  GetLeptonInfoMatching();
  
  return 0;
}


int singletopAnalyser::GetLeptonInfoMatching() {
  Int_t n;
  UInt_t i;
  //Float_t fNearest;
  
  TLorentzVector vec4Res;
  Int_t nIdx, nIdxFound = 0;
  //Float_t fRelIso = 0;
  
  b_ttbar_missedlep = 0;
  
  for ( n = 0 ; n < b_ttbar_channel - 1 ; n++ ) {
    Int_t nIdxGenLep;
    TLorentzVector vec4Gen;
    //Int_t nID;
    
    if ( n == 0 ) {
      nIdxGenLep = m_nIdxGenLep;
      vec4Gen = b_ttbar_genLep1;
      //nID = b_ttbar_lep1_pdgId;
    } else {
      nIdxGenLep = m_nIdxGenLep2nd;
      vec4Gen = b_ttbar_genLep2;
      //nID = b_ttbar_lep2_pdgId;
    }
    
    nIdxFound = -1;
    
    for ( i = 0 ; i < nElectron ; i++ ) {
      for ( nIdx = Electron_genPartIdx[ i ] ; nIdx >= 0 ; nIdx = GenPart_genPartIdxMother[ nIdx ] ) {
        Int_t nIDTmp = abs(GenPart_pdgId[ nIdx ]);
        
        if ( nIDTmp != 11 ) break;
        
        if ( nIdx == nIdxGenLep ) {
          nIdxFound = i;
          break;
        }
      }
      
      if ( nIdxFound >= 0 ) break;
    }
    
    if ( nIdxFound >= 0 ) { // When electron is Found
      //printf("%i\n", nIdxFound);
    } else {
      for ( i = 0 ; i < nMuon ; i++ ) {
        for ( nIdx = Muon_genPartIdx[ i ] ; nIdx >= 0 ; nIdx = GenPart_genPartIdxMother[ nIdx ] ) {
          Int_t nIDTmp = abs(GenPart_pdgId[ nIdx ]);
          
          if ( nIDTmp != 11 && nIDTmp != 13 ) break;
          
          if ( nIdx == nIdxGenLep ) {
            nIdxFound = i;
            break;
          }
        }
          
        if ( nIdxFound >= 0 ) break;
      }
    }
    
    if ( nIdxFound >= 0 ) { // When muon is Found
      //printf("%i\n", nIdxFound);
    } else { // No matching reco lepton
      b_ttbar_missedlep++;
      /*printf("################# %i-th of %i LEPTON (%i) CANNOT BE FOUND! ###################\n", 
        n + 1, b_ttbar_channel - 1, nIdxGenLep);
      
      for ( Int_t i = 0 ; i < (int)nElectron ; i++ ) {
        printf("%i - %i\n", Electron_pdgId[ i ], Electron_genPartIdx[ i ]);
      }
      
      for ( Int_t i = 0 ; i < (int)nMuon ; i++ ) {
        printf("%i - %i\n", Muon_pdgId[ i ], Muon_genPartIdx[ i ]);
      }
      
      printf("================ GEN PARTICLES %10i ==========\n", (int)m_nEventIdx);
      printf("Top idx : %i, %i\n", m_nIdxGenTop, m_nIdxGenTop2nd);
      printf("W idx : %i, %i\n", m_nIdxGenW, m_nIdxGenW2nd);
      if ( b_ttbar_channel >= 2 ) printf("1st lepton idx : %i\n", m_nIdxGenLep);
      if ( b_ttbar_channel >= 3 ) printf("2nd lepton idx : %i\n", m_nIdxGenLep2nd);
      
      for ( Int_t i = 0 ; i < (int)nGenPart ; i++ ) {
        TLorentzVector vec4P4;
        vec4P4.SetPtEtaPhiM(GenPart_pt[ i ], GenPart_eta[ i ], GenPart_phi[ i ], GenPart_mass[ i ]);
        printf("%5i - %7i | %5i %4X %4X (%i %i) | %8.3f %8.3f %8.3f | %8.3f %8.3f %8.3f %8.3f\n", 
          i, GenPart_pdgId[ i ], 
          GenPart_genPartIdxMother[ i ], GenPart_status[ i ], GenPart_statusFlags[ i ], 
          ( GenPart_statusFlags[ i ] >> 12 ) & 0x3, ( GenPart_statusFlags[ i ] >> 7 ) & 0x3,
          GenPart_pt[ i ], GenPart_eta[ i ], GenPart_phi[ i ], 
          vec4P4.X(), vec4P4.Y(), vec4P4.Z(), GenPart_mass[ i ]);
      }
      
      printf("==================================================\n");*/
    }
  }
  
  return 0;
}


// Background study: QCD


int singletopAnalyser::GetGenInfoQCD() {
  return 0;
}


// Background study: WJets


int singletopAnalyser::GetGenInfoWJets() {
  return 0;
}


// Background study: Single top others (s-channel, tW)


int singletopAnalyser::GetGenInfoSTOthers() {
  return 0;
}


////////////////////////////////////////////////////////////////////////////////
// For Gen level test : End
////////////////////////////////////////////////////////////////////////////////


int singletopAnalyser::GetJets() {
  Int_t i;
  
  // Getting the first non-b-tagged jet
  
  for ( i = 0 ; i < b_njet ; i++ ) {
    if ( m_jetsCMVA[ i ] <= 0.9432 ) {
      b_jet1 = m_jets[ i ];
      b_CSVjet1 = m_jetsCMVA[ i ];
      break;
    }
  }
  
  // Finding associated b-tagged jet with the lepton
  
  Float_t fDRMin = 1048576;
  Int_t nIdx = -1;
  
  for ( Int_t i = 0 ; i < b_njet ; i++ ) {
    if ( m_jetsCMVA[ i ] <= 0.9432 ) continue;
    nIdx = i;
    break;
    
    Float_t fDR = b_lep.DeltaR(m_jets[ i ]);
    
    if ( fDRMin > fDR ) {
      fDRMin = fDR;
      nIdx = i;
    }
  }
  
  if ( nIdx >= 0 ) {
    b_bjet1 = m_jets[ nIdx ];
    b_CSVbjet1 = m_jetsCMVA[ i ];
  }
  
  // Getting the second b-tagged jet
  
  for ( i = 0 ; i < b_njet ; i++ ) {
    if ( i != nIdx && m_jetsCMVA[ i ] > 0.9432 ) {
      b_bjet2 = m_jets[ i ];
      b_CSVbjet2 = m_jetsCMVA[ i ];
      break;
    }
  }
  
  return 0;
}


int singletopAnalyser::CalcSphericity() {
  Int_t i, j;
  
  TVector3 arrvec3Jets[ 3 ];
  
  arrvec3Jets[ 0 ] = m_jets[ 0 ].Vect();
  arrvec3Jets[ 1 ] = m_jets[ 1 ].Vect();
  arrvec3Jets[ 2 ] = m_jets[ 2 ].Vect();
  
  Double_t arrSphericity[ 9 ];
  Double_t fNorm = 1.0 / 
    ( arrvec3Jets[ 0 ].Mag2() + arrvec3Jets[ 1 ].Mag2() + arrvec3Jets[ 2 ].Mag2() );
  
  // Setting entries of sphericity matrix
  for ( i = 0 ; i < 3 ; i++ )
  for ( j = 0 ; j < 3 ; j++ ) {
    arrSphericity[ 3 * i + j ]  = arrvec3Jets[ i ].Px() * arrvec3Jets[ j ].Px();
    arrSphericity[ 3 * i + j ] += arrvec3Jets[ i ].Py() * arrvec3Jets[ j ].Py();
    arrSphericity[ 3 * i + j ] += arrvec3Jets[ i ].Pz() * arrvec3Jets[ j ].Pz();
    arrSphericity[ 3 * i + j ] *= fNorm;
  }
  
  // Calculating eigenvalues
  TMatrixDSym matSphericity(3, arrSphericity);
  TMatrixDSymEigen eigenValues(matSphericity);
  
  Double_t dE1, dE2, dE3;
  
  dE1 = eigenValues.GetEigenValues()[ 0 ];
  dE2 = eigenValues.GetEigenValues()[ 1 ];
  dE3 = eigenValues.GetEigenValues()[ 2 ];
  
  // Align dEx and then putting them into b_jetQx
  if ( dE1 > dE2 ) {
    b_jetQ1 = dE1;
    b_jetQ3 = dE2;
  } else {
    b_jetQ1 = dE2;
    b_jetQ3 = dE1;
  }
  
  if ( b_jetQ1 < dE3 ) b_jetQ1 = dE3;
  if ( b_jetQ3 > dE3 ) b_jetQ3 = dE3;
  b_jetQ2 = dE1 + dE2 + dE3 - b_jetQ1 - b_jetQ3;
  
  b_jetC = 3.0 * ( b_jetQ1 * b_jetQ2 + b_jetQ2 * b_jetQ3 + b_jetQ3 * b_jetQ1 );
  
  return 0;
}


TLorentzVector singletopAnalyser::RecoWFromTop(double *pdDiffMET, int nFlag) {
  TLorentzVector vec4W;
  
  double dMW = 80.4;
  
  switch ( nFlag ) {
  case MY_FLAG_RECOW_WRESTRICTION: 
  case MY_FLAG_RECOW_WRESTR_HIGH: 
  case MY_FLAG_RECOW_WRESTR_IMG: 
  case MY_FLAG_RECOW_WMASS_IN_GEN: 
  {
    // Calculating the 4-momentum of neutrino
    // The following method comes from AN_2012_448, p. 11, in 
    // http://cms.cern.ch/iCMS/analysisadmin/cadilines?line=TOP-13-001
    // or https://indico.cern.ch/event/472719/contributions/2166641/attachments/
    // 1274021/1889392/2016_05_17_toplhc.pdf
    // But the choosing in two solutions is different; 
    // the one with smallest deltaR between reco W and reco B.
    // Prof. Kim pointed out that this method can have a bias on low pT of top. Check it out.
    if ( nFlag == MY_FLAG_RECOW_WMASS_IN_GEN ) {
      dMW = GenPart_mass[ m_nIdxGenW ];
    }
    
    double dELep, dPzLep, dMTSqr;
    double dPtMMul, dDetSqr;
    double dPzN1, dPzN2, dPzN = 0;
    double dMET, dMETPhi;
    
    dELep = b_lep.E();
    dPzLep = b_lep.Pz();
    dMTSqr = dELep * dELep - dPzLep * dPzLep;
    
    dMET = b_met;
    dMETPhi = b_met_phi;
    
    dPtMMul = 0.5 * dMW * dMW + dMET * b_lep.Pt() * cos(b_lep.Phi() - dMETPhi);
    dDetSqr = dPtMMul * dPtMMul - dMET * dMET * dMTSqr;
    
    if ( dDetSqr >= 0 ) {
      if ( dDetSqr > 0 ) b_truthtop1_lowhighOn = 1;
      
      dDetSqr = sqrt(dDetSqr);
      dPzN1 = ( dPzLep * dPtMMul + dELep * dDetSqr ) / dMTSqr;
      dPzN2 = ( dPzLep * dPtMMul - dELep * dDetSqr ) / dMTSqr;
      
      if ( nFlag != MY_FLAG_RECOW_WRESTR_HIGH ) {
        dPzN = ( abs(dPzN1) < abs(dPzN2) ? dPzN1 : dPzN2 );
      } else {
        dPzN = ( abs(dPzN1) > abs(dPzN2) ? dPzN1 : dPzN2 );
      }
    } else { // In this case, dDetSqr = 0, but dMET is changed to be adjusted with this
      dPzN = dPzLep * dPtMMul / dMTSqr; // Just setting dDetSqr = 0
      dMET = dPtMMul / sqrt(dMTSqr);
      
      b_truthtop1_imaginaryOn = 1;
    }
    
    // Keeping the modified MET to use on correction of the b-jet pT
    if ( pdDiffMET != NULL ) *pdDiffMET = dMET - b_met;
    
    TLorentzVector vec4Neu;
    vec4Neu.SetPxPyPzE(dMET * cos(dMETPhi), dMET * sin(dMETPhi), 
      dPzN, sqrt(dMET * dMET + dPzN * dPzN)); // Because m_neutrino ~ 0
    
    // End of calculating of 4-momentum of neutrino
    
    vec4W = b_lep + vec4Neu;
    
    break;
  }
  case MY_FLAG_RECOW_NEU_IN_GEN:
    vec4W = b_lep + Get4VecGen(m_nIdxGenNeu);
    break;
  case MY_FLAG_RECOW_W_IN_GEN:
    vec4W = Get4VecGen(m_nIdxGenW);
    break;
  }
  
  return vec4W;
}


int singletopAnalyser::RecoTop() {
  TLorentzVector vec4BJet = b_bjet1;
  TLorentzVector vec4WLow, vec4WHigh;
  double dMETMod;
  
  vec4WLow =  RecoWFromTop(&dMETMod);
  vec4WHigh = RecoWFromTop(&dMETMod, MY_FLAG_RECOW_WRESTR_HIGH);
  
  b_top1_lower  = vec4BJet + vec4WLow;
  b_top1_higher = vec4BJet + vec4WHigh;
  
  TLorentzVector vec4NeuLow  = vec4WLow  - b_lep;
  TLorentzVector vec4NeuHigh = vec4WHigh - b_lep;
  
  // 1 : Low, 2 : High
  b_DiffLepMom11 = b_lep + vec4NeuLow;
  b_DiffLepMom12 = b_lep + vec4NeuLow;
  b_DiffLepMom21 = b_lep + vec4NeuHigh;
  b_DiffLepMom22 = b_lep + vec4NeuHigh;
  
  b_DiffLepMom11.Boost(-vec4WLow.BoostVector());
  b_DiffLepMom12.Boost(-vec4WHigh.BoostVector());
  b_DiffLepMom21.Boost(-vec4WLow.BoostVector());
  b_DiffLepMom22.Boost(-vec4WHigh.BoostVector());
  
  if ( m_isMC && m_isSig ) {
    b_top1_imaginary = vec4BJet + RecoWFromTop(NULL, MY_FLAG_RECOW_WRESTR_IMG);
    b_top1_genWMass  = vec4BJet + RecoWFromTop(NULL, MY_FLAG_RECOW_WMASS_IN_GEN);
    b_top1_genneu = vec4BJet + RecoWFromTop(NULL, MY_FLAG_RECOW_NEU_IN_GEN);
    b_top1_genW = vec4BJet + RecoWFromTop(NULL, MY_FLAG_RECOW_W_IN_GEN);
    
    b_top1_genWgenB = Get4VecGen(m_nIdxGenBFromT) + RecoWFromTop(NULL, MY_FLAG_RECOW_W_IN_GEN);
    
    b_truthtop1_lowhighNeuPT = log(b_gentop1.DeltaR(b_top1_lower) / b_gentop1.DeltaR(b_top1_higher));
  }
  
  b_W1 = ( vec4BJet.DeltaR(vec4WLow) < vec4BJet.DeltaR(vec4WHigh) ? vec4WLow : vec4WHigh );
  b_top1 = vec4BJet + b_W1;
  
  b_recoW_lower  = vec4WLow;
  b_recoW_higher = vec4WHigh;
  b_recoNeu_lower  = vec4NeuLow;
  b_recoNeu_higher = vec4NeuHigh;
  
  return 0;
}


int singletopAnalyser::CalcRecoCosStar() {
  TLorentzVector vec4LepRest = b_lep, vec4JetRest = b_jet1;
  
  vec4LepRest.Boost(-b_top1.BoostVector());
  vec4JetRest.Boost(-b_top1.BoostVector());
  
  TVector3 vec3Lep = vec4LepRest.Vect().Unit();
  TVector3 vec3Jet = vec4JetRest.Vect().Unit();
  double dCos = vec3Lep.Dot(vec3Jet);
  
  m_h1CosSReco->Fill(dCos);
  
  b_cos_star = dCos;
  
  return 0;
}


int singletopAnalyser::RunEvt() {
  int nRes = 0;
  
  // Begin of jobs IN GEN LEVEL
  
  if ( m_isMC && m_isSig ) {
    nRes = GetGenInfoSignal();
  } else {
    switch ( m_nBkgType ) {
      case MY_FLAG_BKGTYPE_TTBAR:     nRes = GetGenInfoTTbar(); break;
      case MY_FLAG_BKGTYPE_QCD:       nRes = GetGenInfoQCD(); break;
      case MY_FLAG_BKGTYPE_WJETS:     nRes = GetGenInfoWJets(); break;
      case MY_FLAG_BKGTYPE_ST_OTHERS: nRes = GetGenInfoSTOthers(); break;
    }
  }
  
  //if ( nRes != 0 ) return nRes;
  
  // End of jobs IN GEN LEVEL
  
  //nRes = ObjectSelection();
  nRes = EventSelection();
  if ( nRes < 4 ) return 0;
  
  b_met_sumEt = MET_sumEt;
  b_filter_met = ( Flag_METFilters ? 1 : 0 );
  
  GetJets();
  if ( b_njet >= 3 ) CalcSphericity();
  
  b_fatjet.SetPtEtaPhiM(FatJet_pt[ 0 ], FatJet_eta[ 0 ], FatJet_phi[ 0 ], FatJet_mass[ 0 ]);
  b_nfatjet = nFatJet;
  b_nsubjet = nSubJet;
  
  b_nhadcand = nhad;
  
  /*if ( m_isMC && !m_isSig && m_nBkgType == MY_FLAG_BKGTYPE_TTBAR ) {
    if ( b_lep_pid != b_ttbar_lep1_pdgId || b_lep.DeltaR(b_ttbar_lep1) > b_lep.DeltaR(b_ttbar_lep2) ) {
      TLorentzVector vec4Swap;
      Int_t nIDSwap;
      Float_t fRelIsoSwap;
      
      vec4Swap = b_ttbar_lep1; b_ttbar_lep1 = b_ttbar_lep2; b_ttbar_lep2 = vec4Swap;
      nIDSwap = b_ttbar_lep1_pdgId; b_ttbar_lep1_pdgId = b_ttbar_lep2_pdgId; b_ttbar_lep2_pdgId = nIDSwap;
      fRelIsoSwap = b_ttbar_lep1_reliso; b_ttbar_lep1_reliso = b_ttbar_lep2_reliso; b_ttbar_lep2_reliso = fRelIsoSwap;
    }
  }*/
  if ( m_isMC && !m_isSig && m_nBkgType == MY_FLAG_BKGTYPE_TTBAR ) {
    if ( b_ttbar_channel > 1 ) {
      if ( abs(b_lep_pid) == 13 ) {
        Int_t nIdx = -1;
        for ( Int_t i = 0 ; i < (int)nMuon ; i++ ) {
          TLorentzVector mom;
          mom.SetPtEtaPhiM(Muon_pt[ i ], Muon_eta[ i ], Muon_phi[ i ], Muon_mass[ i ]);
          if ( mom.DeltaR(b_lep) < 0.001 ) {
            nIdx = i;
            break;
          }
        }
        
        b_ttbar_lep1_deltaRGen = ( !( ( abs(b_ttbar_lep1_pdgId) == 13 && nIdx == m_nIdxGenLep ) || 
          ( abs(b_ttbar_lep2_pdgId) == 13 && nIdx == m_nIdxGenLep2nd ) ) ? 1 : -1 );
      } else {
        
      }
      
      if ( m_nIdxGenLep ) {
      }
    }
  }
  
  // Reconstruction of top quark
  
  RecoTop();
  if ( ( b_njet == 2 && b_nbjet == 1 ) || ( b_njet == 3 && b_nbjet == 2 ) )
    CalcRecoCosStar();
  
  return 0;
}


void singletopAnalyser::Loop() {
  Long64_t i;
  
  int nRes;
  
  if ( fChain == 0 ) return;

  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  
  for ( i = 0 ; i < nentries ; i++ ) {
    //Prepare for new loop
    resetBranch();
    Long64_t ientry = LoadTree(i);
    if ( ientry < 0 ) break;
    
    nb = fChain->GetEntry(i);
    nbytes += nb;
    
    m_nEventIdx = i;
    nRes = RunEvt();
    
    if ( nRes == 0 || nRes > 0 ) { // ==0 : passed on reco, >0 : failed on reco but passed on gen
      if ( nRes != 0 ) b_onlyGen = 1;
      m_tree->Fill();
    }
  }
}


int main(int argc, char **argv) {
  int i;
  
  string env = getenv("CMSSW_BASE");
  string username = getenv("USER");
  
  bool bIsMC;
  int nIdxStart = 0, nIdxEnd = 1048575 * 1024;
  
  TFile *fRoot;
  int nTrial;
  
  int nBkgType;
  
  // Copied from h2muAnaylser.cc to make sure the compatibility
  if (argc > 1 && strcmp(argv[ 1 ], "-q") != 0 ) {
    string dirName = "root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/"+username+"/nanoAOD/"+std::string(argv[1])+"/"+std::string(argv[2]);
   // string dirName = env+("/src/nano/analysis/test/h2mu/Results/")+argv[1]+"/"+argv[2];
    string temp = argv[2];
    Bool_t isMC = false;
    Size_t found = temp.find("Run");
    if(found == string::npos) isMC = true;
    for (Int_t i = 3; i < argc; i++) {
      std::string inputName = argv[ i ];
      TFile *f = TFile::Open(argv[i], "read");

      TTree *tree;
      f->GetObject("Events", tree);
      
      nBkgType = 0;
      if ( inputName.find("TT") != std::string::npos ) nBkgType = singletopAnalyser::MY_FLAG_BKGTYPE_TTBAR;
      if ( inputName.find("QCD") != std::string::npos ) nBkgType = singletopAnalyser::MY_FLAG_BKGTYPE_QCD;
      if ( inputName.find("WJets") != std::string::npos ) nBkgType = singletopAnalyser::MY_FLAG_BKGTYPE_WJETS;
      if ( inputName.find("ST_") != std::string::npos ) nBkgType = singletopAnalyser::MY_FLAG_BKGTYPE_ST_OTHERS;

      temp = argv[i];
      found = temp.find_last_of('/');
      string outPutName = dirName+temp.substr(found);
      singletopAnalyser t(tree, NULL, NULL, isMC, 
        std::string(argv[i]).find("ST_t-channel") != std::string::npos, nBkgType, true);
      t.setOutput(outPutName);
      t.Loop();
    }
    
    return 0;
  }
  
  if ( argc < 3 ) {
    printf("Usage: \n"
      "(on single core)\n"
      " $ singletopAnalyser [DATASET NAME] [MC or RD] "
      "[idx of the first root file] [(idx+1) of the last root file]\n"
      "(on single core; for the author)\n"
      " $ singletopAnalyser -q [LIST_FILE_OF_FILES] [MC or RD] "
      "[idx of the first root file] [(idx+1) of the last root file]\n"
      "(for grid job (condor; in KISTI) \n"
      " $ python makejobs.py [DEST DIR] j`python batch_nanoAOD.py n [# OF FILES PER JOB]`\n");
    
    return 1;
  }
  
  if ( argc >= 4 ) nIdxStart = atoi(argv[ 4 ]);
  if ( argc >= 5 ) nIdxEnd = atoi(argv[ 5 ]);
  
  bIsMC = ( strcmp(argv[ 3 ], "RD") != 0 );
  
  std::ifstream fileList(argv[ 2 ]);
  std::string strLine;
  
  if ( !fileList.is_open() ) {
    perror("Error: Cannot open the list file\n");
    return 1;
  }
  
  // Setup for reading on cache
  TFile::SetCacheFileDir(".");
  
  for ( i = 0 ; getline(fileList, strLine) ; i++ ) {
    i--;
    if ( strLine.substr(0, 1) == "#" ) continue;
    if ( strLine.find_first_not_of(" \t") == string::npos ) continue;
    i++;
    
    if ( i <  nIdxStart ) continue;
    if ( i >= nIdxEnd )   break;
    
    cout << "Start : " << strLine << endl;
    
    fRoot = NULL;
    nTrial = 0;
    
    while ( fRoot == NULL ) {
      fRoot = TFile::Open(strLine.c_str(), "CACHEREAD"); // reading on cache
      //fRoot = TFile::Open(strLine.c_str());
      
      if ( fRoot == NULL ) sleep(15 * 1000);
      if ( nTrial++ > 5 ) break;
    }
    
    if ( fRoot == NULL ) {
      perror("Error: Failed to open the root file\n");
      return 1;
    }
    
    nBkgType = 0;
    if ( strLine.find("TT") != std::string::npos )    nBkgType = singletopAnalyser::MY_FLAG_BKGTYPE_TTBAR;
    if ( strLine.find("QCD") != std::string::npos )   nBkgType = singletopAnalyser::MY_FLAG_BKGTYPE_QCD;
    if ( strLine.find("WJets") != std::string::npos ) nBkgType = singletopAnalyser::MY_FLAG_BKGTYPE_WJETS;
    if ( strLine.find("ST_") != std::string::npos )   nBkgType = singletopAnalyser::MY_FLAG_BKGTYPE_ST_OTHERS;
    
    TTree *treeMain = (TTree *)fRoot->Get("Events");
    
    singletopAnalyser t(treeMain, treeMain, treeMain, bIsMC, 
      strLine.find("ST_t-channel") != std::string::npos, nBkgType, true);
    
    std::string strOutput = "out";
    strOutput = strOutput + std::to_string(i) + ".root";
    t.setOutput(strOutput);
    
    t.Loop();
    
    //fRoot->Close(); // Don't do that; it's done in ~Events()
    cout << "End   : " << strLine << endl;
  }

  return 0;
}


