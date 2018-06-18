//#define nanoAnalyser_cxx
#define Events_cxx
#include "DataFormats/HepMCCandidate/interface/GenStatusFlags.h"
#include "nano/analysis/interface/nanoBase.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <cstdlib>


using namespace std;
/*
To compile:
cd $CMSSW_BASE ; scram b -j 8 ; cd -

To run:
singletopAnalyser [LIST_FILE_OF_FILES] ["MC" or "RD"] [idx of the first root file] [(idx+1) of the last root file]

To throw jobs to condor:
python makejobs.py nanoAOD j`python batch_nanoAOD.py n 2`
*/


class singletopAnalyser : public nanoBase {
private: 
  TH1F *m_h1Weights, *m_h1GenWeights;
  
  TH1F *m_h1CosS;
  TH1F *m_h1CosSReco;
  TH1F *m_h1IsHighPtMuon, *m_h1DRMuon;
  
  //Bool_t m_isMC;
  Bool_t m_isSig;
  Bool_t m_isFullGen;
  
  Int_t m_nIdxGenTop, m_nIdxGenAssoQ, m_nIdxGenLep, m_nIdxGenBFromT, m_nIdxGenW;
  Int_t m_nIDGenLep;
  
  Int_t m_nIdxGenNeu, m_nIDGenNeu;
  
  Int_t b_onlyGen;
  
  Int_t b_objstep, b_step, b_channel, b_njets, b_nbjets;
  
  TLorentzVector b_gentop1;
  TLorentzVector b_genW;
  
  TParticle recolep1;
  
  TLorentzVector b_lep1, b_jet1, b_bjet1, b_bjet2;
  Int_t b_lep1_pid;
  Float_t b_met, b_met_phi;
  
  TLorentzVector b_DiffLepMom11;
  TLorentzVector b_DiffLepMom12;
  TLorentzVector b_DiffLepMom21;
  TLorentzVector b_DiffLepMom22;
  
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
  
  Bool_t b_tri;
  Float_t b_triw, b_triw_up, b_triw_dn;
  
  Int_t b_npvs;
  Float_t b_pvz, b_pv_ndof;
  
  Float_t b_weight, b_genweight;
  Float_t b_puweight, b_puweight_up, b_puweight_dn;
  Float_t b_btagweight, b_btagweight_up, b_btagweight_dn;
  Float_t b_mueffweight, b_mueffweight_up, b_mueffweight_dn;
  Float_t b_eleffweight, b_eleffweight_up, b_eleffweight_dn;
  std::vector<Float_t> b_csvweights;
  
  Float_t b_cos_star_gen;
  Float_t b_cos_star_bjet_gen;
  Float_t b_cos_labf_gen;
  Float_t b_cos_labf_bjet_gen;
  
  Float_t b_cos_star;
  
  
public: 
  //set output file
  void setOutput(std::string outputName);
  void MakeBranch(TTree* t);
  void resetBranch();
  void LoadModules(pileUpTool* pileUp, lumiTool* lumi);
  //void collectTMVAvalues();
  singletopAnalyser(TTree *tree=0, Bool_t flag = false, Bool_t isSig = false, Bool_t isFullGen = false);
  ~singletopAnalyser();
  virtual void     Loop();
  
  int DoMoreGenLvl(TLorentzVector &vec4Top, TLorentzVector &vec4Lep, TLorentzVector &vec4AssoQ);
  int RunEvt();
  
  int GetIdxGenTop();
  int GetIdxGenAssoQuark();
  int GetIdxGenLepton();
  int GetIdxGenBFromTop();
  int GetIdxGenWFromTop();
  int GetIdxGenNeutrino();
  
  bool hasOverLap(TLorentzVector cand, vector<TParticle> objects, Float_t rad);
  Double_t roccoR(TLorentzVector m, int &q, int &nGen, int &nTrackerLayers);
  
  vector<TParticle> muonSelection();
  vector<TParticle> elecSelection();
  vector<TParticle> jetSelection(vector<int> &vecIdx);
  vector<TParticle> bjetSelection(vector<int> &vecIdx);
  
  TLorentzVector Get4VecGen(int nIdx);
  TLorentzVector Get4VecElec(int nIdx);
  TLorentzVector Get4VecMuon(int nIdx);
  TLorentzVector Get4VecAJet(int nIdx);
  TLorentzVector Get4VecBJet(int nIdx);
  
#define MY_FLAG_RECOW_WRESTRICTION  0
#define MY_FLAG_RECOW_WRESTR_HIGH   1
#define MY_FLAG_RECOW_WRESTR_IMG    2
#define MY_FLAG_RECOW_WMASS_IN_GEN  3
#define MY_FLAG_RECOW_NEU_IN_GEN    4
#define MY_FLAG_RECOW_W_IN_GEN      5
  
  TLorentzVector RecoWFromTop(double *pdDiffMET, int nFlag = MY_FLAG_RECOW_WRESTRICTION);
  int RecoTop();
  int CalcRecoCosStar();
};


singletopAnalyser::singletopAnalyser(TTree *tree, Bool_t flag, Bool_t isSig, Bool_t isFullGen) : nanoBase(tree, flag), 
  m_isSig(isSig), m_isFullGen(isFullGen)
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
  
  m_h1Weights = new TH1F("weight", "weight", 2, -1, 1);
  m_h1GenWeights = new TH1F("genweight", "genweight", 2, -1, 1);
  
  m_h1CosS = new TH1F("cos_star", "cos_star", 10, -1, 1);
  m_h1CosSReco = new TH1F("cos_star_reco", "cos_star_reco", 10, -1, 1);
  m_h1IsHighPtMuon = new TH1F("isClosest", "isClosest", 2, 0, 2);
  m_h1DRMuon = new TH1F("deltaR_muon", "deltaR_muon", 50, 0, 5);
}


void singletopAnalyser::MakeBranch(TTree *t) {
  t->Branch("onlyGen", &b_onlyGen, "onlyGen/I");
  
  t->Branch("objstep", &b_objstep, "objstep/I");
  t->Branch("step", &b_step, "step/I");
  t->Branch("channel", &b_channel, "channel/I");
  t->Branch("njet", &b_njets, "njet/I");
  t->Branch("nbjet", &b_nbjets, "nbjet/I");
  
  t->Branch("gentop1", "TLorentzVector", &b_gentop1);
  t->Branch("genW", "TLorentzVector", &b_genW);
  
  t->Branch("lep1", "TLorentzVector", &b_lep1);
  t->Branch("jet1", "TLorentzVector", &b_jet1);
  t->Branch("bjet1", "TLorentzVector", &b_bjet1);
  t->Branch("bjet2", "TLorentzVector", &b_bjet2);
  
  t->Branch("DiffLepMom11", "TLorentzVector", &b_DiffLepMom11);
  t->Branch("DiffLepMom12", "TLorentzVector", &b_DiffLepMom12);
  t->Branch("DiffLepMom21", "TLorentzVector", &b_DiffLepMom21);
  t->Branch("DiffLepMom22", "TLorentzVector", &b_DiffLepMom22);
  
  t->Branch("lep1_pid", &b_lep1_pid, "lep1_pid/I");
  
  t->Branch("met", &b_met, "met/F");
  t->Branch("met_phi", &b_met_phi, "met_phi/F");
  
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
  
  t->Branch("trigger", &b_tri, "trigger/O");
  t->Branch("tri", &b_triw, "tri/F");
  t->Branch("tri_up", &b_triw_up, "tri_up/F");
  t->Branch("tri_dn", &b_triw_dn, "tri_dn/F");
  
  t->Branch("npvs", &b_npvs, "npvs/I");
  t->Branch("PVz", &b_pvz, "PVz/F");
  t->Branch("PV_nDoF", &b_pv_ndof, "PV_nDoF/F");
  
  t->Branch("weight", &b_weight, "weight/F");
  
  t->Branch("puweight", &b_puweight, "puweight/F");
  t->Branch("puweight_up", &b_puweight_up, "puweight_up/F");
  t->Branch("puweight_dn", &b_puweight_dn, "puweight_dn/F");
  
  t->Branch("genweight", &b_genweight, "genweight/F");
  t->Branch("csvweight", "std::vector<float>", &b_csvweights);
  
  t->Branch("btagweight", &b_btagweight, "btagweight/F");
  t->Branch("btagweight_up", &b_btagweight_up, "btagweight_up/F");
  t->Branch("btagweight_dn", &b_btagweight_dn, "btagweight_dn/F");
  
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
  b_onlyGen = 0;
  
  b_objstep = -1;
  b_step = -1;
  b_channel = 0;
  b_njets = 0;
  b_nbjets = 0;
  
  b_gentop1.SetPtEtaPhiM(0, 0, 0, 0);
  b_genW.SetPtEtaPhiM(0, 0, 0, 0);
  
  b_lep1.SetPtEtaPhiM(0, 0, 0, 0);
  b_jet1.SetPtEtaPhiM(0, 0, 0, 0);
  b_bjet1.SetPtEtaPhiM(0, 0, 0, 0);
  b_bjet2.SetPtEtaPhiM(0, 0, 0, 0);
  b_met = 0;
  b_met_phi = 0;
  b_lep1_pid = 0;
  
  b_DiffLepMom11.SetPtEtaPhiM(0, 0, 0, 0);
  b_DiffLepMom12.SetPtEtaPhiM(0, 0, 0, 0);
  b_DiffLepMom21.SetPtEtaPhiM(0, 0, 0, 0);
  b_DiffLepMom22.SetPtEtaPhiM(0, 0, 0, 0);
  
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
  
  b_tri = false;
  b_triw = b_triw_up = b_triw_dn = 1.0;
  
  b_npvs = 0;
  b_pvz = 0.0;
  b_pv_ndof = 0;
  
  b_weight = b_genweight = 1.0;
  b_puweight = b_puweight_up = b_puweight_dn = 1.0;
  b_btagweight = b_btagweight_up = b_btagweight_dn = 1.0;
  b_mueffweight = b_mueffweight_up = b_mueffweight_dn = 1.0;
  b_eleffweight = b_eleffweight_up = b_eleffweight_dn = 1.0;
  b_csvweights.clear();
  
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


int singletopAnalyser::GetIdxGenLepton() {
  UInt_t i;
  Int_t j;
  
  // Finding the lepton from the top quark
  for ( i = 0 ; i < nGenPart ; i++ ) {
    int nPID = abs(GenPart_pdgId[ i ]);
    int nPassW = 0;
    
    if ( ( nPID == 11 || nPID == 13 ) ) {
      for ( j = GenPart_genPartIdxMother[ i ] ; j >= 0 ; j = GenPart_genPartIdxMother[ j ] ) {
        int nParentID = abs(GenPart_pdgId[ j ]);
        
        if ( nParentID == 24 ) {
          int nGrandParentID = abs(GenPart_pdgId[ GenPart_genPartIdxMother[ j ] ]);
          if ( nGrandParentID == 6 ) {
            nPassW = 1;
            m_nIdxGenW = j;
          }
        }
        
        if ( nPassW == 1 && j == m_nIdxGenTop ) break;
      }
      
      if ( j < 0 ) continue;
      
      m_nIdxGenLep = i;
      m_nIDGenLep = GenPart_pdgId[ i ];
      break;
    }
  }
  
  return ( i < nGenPart ? 0 : -1 );
}


int singletopAnalyser::GetIdxGenBFromTop() {
  UInt_t i;
  
  // Finding the b quark from the top quark
  for ( i = 0 ; i < nGenPart ; i++ ) {
    int nPID = GenPart_pdgId[ i ];
    if ( abs(nPID) != 5 ) continue;
    
    int nIdxMother = GenPart_genPartIdxMother[ i ];
    if ( nIdxMother < 0 ) continue;
    
    if ( GenPart_pdgId[ nIdxMother ] * nPID != 6 * 5 ) continue;
    
    m_nIdxGenBFromT = i;
    break;
  }
  
  return ( i < nGenPart ? 0 : -1 );
}


int singletopAnalyser::GetIdxGenWFromTop() {
  UInt_t i;
  
  // Finding the W from the top quark
  for ( i = 0 ; i < nGenPart ; i++ ) {
    int nPID = GenPart_pdgId[ i ];
    if ( abs(nPID) != 24 ) continue;
    
    int nIdxMother = GenPart_genPartIdxMother[ i ];
    if ( nIdxMother < 0 ) continue;
    
    if ( GenPart_pdgId[ nIdxMother ] * nPID != 6 * 24 ) continue;
    
    m_nIdxGenW = i;
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


////////////////////////////////////////////////////////////////////////////////
// For Gen level test : End
////////////////////////////////////////////////////////////////////////////////


bool singletopAnalyser::hasOverLap(TLorentzVector cand, vector<TParticle> objects, Float_t rad) {
  for ( auto obj: objects ) {
    TLorentzVector mom;
    obj.Momentum(mom);
    if ( cand.DeltaR(mom) < rad ) return true;
  }
  return false;
}


Double_t singletopAnalyser::roccoR(TLorentzVector m, int &q, int &nGen, int &nTrackerLayers) {
  Float_t u1 = gRandom->Rndm();
  Float_t u2 = gRandom->Rndm();
  
  if ( !m_isMC ) {
    return m_rocCor->kScaleDT(q, m.Pt(), m.Eta(), m.Phi(), 0, 0);
  } else {
    if (nGen > -1) {
      return m_rocCor->kScaleFromGenMC(q, m.Pt(), m.Eta(), m.Phi(),
        nTrackerLayers, GenPart_pt[ nGen ], u1, 0, 0);
    } else {
      return m_rocCor->kScaleAndSmearMC(q, m.Pt(), m.Eta(), m.Phi(),
        nTrackerLayers, u1, u2, 0, 0);
    }
  }
}


vector<TParticle> singletopAnalyser::muonSelection() {
  UInt_t i;
  
  vector<TParticle> muons;
  
  for ( i = 0 ; i < nMuon ; i++ ) {
    if ( !Muon_trackerMu[ i ] ) continue;
    if ( !Muon_globalMu[ i ] ) continue;
    if ( !Muon_tightId[ i ] ) continue;
    if ( Muon_pfRelIso04_all[ i ] > 0.25 ) continue;
    
    if ( Muon_pt[ i ] < 10 ) continue;
    if ( std::abs(Muon_eta[ i ]) > 2.4 ) continue;
    
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Muon_pt[ i ], Muon_eta[ i ], Muon_phi[ i ], Muon_mass[ i ]);
    //mom = mom * roccoR(mom, Muon_charge[i], Muon_genPartIdx[i], Muon_nTrackerLayers[i]);
    
    auto muon = TParticle();
    muon.SetPdgCode(13 * Muon_charge[ i ]* -1);
    muon.SetMomentum(mom);
    
    muons.push_back(muon);
  }
  
  return muons;
}


vector<TParticle> singletopAnalyser::elecSelection() {
  UInt_t i;
  
  vector<TParticle> elecs;
  
  for ( i = 0 ; i < nElectron ; i++ ){
    if ( Electron_pt[ i ] < 10 ) continue;
    if ( std::abs(Electron_eta[ i ]) > 2.5 ) continue;
    if ( Electron_cutBased[ i ] < 3 ) continue;
    
    float el_scEta = Electron_deltaEtaSC[ i ] + Electron_eta[ i ];
    if ( std::abs(el_scEta) > 1.4442 && std::abs(el_scEta) < 1.566 ) continue;
    
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Electron_pt[ i ], Electron_eta[ i ], Electron_phi[ i ], Electron_mass[ i ]);
    
    if ( hasOverLap(mom, elecs, 0.3) ) continue;
   
    auto elec = TParticle();
    elec.SetPdgCode(11 * Electron_charge[ i ] * -1);
    elec.SetMomentum(mom);
    elec.SetWeight(el_scEta);
    
    elecs.push_back(elec);
  }
  
  return elecs;
}


vector<TParticle> singletopAnalyser::jetSelection(vector<int> &vecIdx) {
  UInt_t i, j;
  
  vector<TParticle> jets; 
  float Jet_SF_CSV[ 19 ] = {1.0,};
  
  for ( i = 0; i < nJet ; i++ ){
    if ( std::abs(Jet_eta[ i ]) > 4.7 ) continue;
    if ( std::abs(Jet_eta[ i ]) >= 2.4 && Jet_pt[ i ] < 30 ) continue;
    if ( std::abs(Jet_eta[ i ]) < 2.4 && Jet_pt[ i ] < 20 ) continue;
    if ( Jet_jetId[ i ] < 1 ) continue;
    
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Jet_pt[ i ], Jet_eta[ i ], Jet_phi[ i ], Jet_mass[ i ]);
    
    if ( mom.TLorentzVector::DeltaR(b_lep1) < 0.4 ) continue;
    
    auto jet = TParticle();
    jet.SetMomentum(mom);
    jets.push_back(jet);
    
    for ( j = 0 ; j < 19 ; j++ ) {   
      //Jet_SF_CSV[ j ] *= m_btagSF.getSF(jet, Jet_btagCSVV2[ i ], Jet_hadronFlavour[ i ], j);
    }
    
    vecIdx.push_back(i);
  }
  
  for ( i = 0 ; i < 19 ; i++ ) b_csvweights.push_back(Jet_SF_CSV[ i ]);
  b_btagweight = Jet_SF_CSV[ 0 ];
  
  return jets;
}


vector<TParticle> singletopAnalyser::bjetSelection(vector<int> &vecIdx) {
  UInt_t i, j;
  
  vector<TParticle> bjets; 
  
  for ( i = 0; i < nJet ; i++ ){
    if ( Jet_pt[ i ] < 20 ) continue;
    if ( std::abs(Jet_eta[ i ]) > 2.4 ) continue;
    if ( Jet_jetId[ i ] < 1 ) continue;
    if ( Jet_btagCSVV2[ i ] < 0.8484 ) continue;
    
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Jet_pt[ i ], Jet_eta[ i ], Jet_phi[ i ], Jet_mass[ i ]);
    if ( mom.TLorentzVector::DeltaR(b_lep1) < 0.4 ) continue;
    
    auto bjet = TParticle();
    bjet.SetMomentum(mom);
    bjets.push_back(bjet);
    
    for ( j = 0 ; j < vecIdx.size() ; j++ ) if ( vecIdx[ j ] >= 0 && (UInt_t)vecIdx[ j ] == i ) break;
    if ( j < vecIdx.size() ) vecIdx[ j ] = -1 * vecIdx[ j ] - 1;
  }
  
  return bjets;
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
    // The following method comes from 
    // https://indico.cern.ch/event/472719/contributions/2166641/attachments/
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
    
    dELep = b_lep1.E();
    dPzLep = b_lep1.Pz();
    dMTSqr = dELep * dELep - dPzLep * dPzLep;
    
    dMET = b_met;
    dMETPhi = b_met_phi;
    
    dPtMMul = 0.5 * dMW * dMW + dMET * b_lep1.Pt() * cos(b_lep1.Phi() - dMETPhi);
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
    
    vec4W = b_lep1 + vec4Neu;
    
    break;
  }
  case MY_FLAG_RECOW_NEU_IN_GEN:
    vec4W = b_lep1 + Get4VecGen(m_nIdxGenNeu);
    break;
  case MY_FLAG_RECOW_W_IN_GEN:
    vec4W = Get4VecGen(m_nIdxGenW);
    break;
  }
  
  return vec4W;
}


int singletopAnalyser::RecoTop() {
  TLorentzVector vec4BJet = ( b_nbjets == 1 || 
    b_bjet1.DeltaR(b_lep1) < b_bjet2.DeltaR(b_lep1) ? b_bjet1 : b_bjet2 );
  TLorentzVector vec4WLow, vec4WHigh;
  double dMETMod;
  
  vec4WLow =  RecoWFromTop(&dMETMod);
  vec4WHigh = RecoWFromTop(&dMETMod, MY_FLAG_RECOW_WRESTR_HIGH);
  
  b_top1_lower  = vec4BJet + vec4WLow;
  b_top1_higher = vec4BJet + vec4WHigh;
  
  TLorentzVector vec4NeuLow  = vec4WLow  - b_lep1;
  TLorentzVector vec4NeuHigh = vec4WHigh - b_lep1;
  
  // 1 : Low, 2 : High
  b_DiffLepMom11 = b_lep1 + vec4NeuLow;
  b_DiffLepMom12 = b_lep1 + vec4NeuLow;
  b_DiffLepMom21 = b_lep1 + vec4NeuHigh;
  b_DiffLepMom22 = b_lep1 + vec4NeuHigh;
  
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
  
  b_top1 = vec4BJet + ( vec4BJet.DeltaR(vec4WLow) < vec4BJet.DeltaR(vec4WHigh) ? vec4WLow : vec4WHigh );
  
  b_recoW_lower  = vec4WLow;
  b_recoW_higher = vec4WHigh;
  b_recoNeu_lower  = vec4NeuLow;
  b_recoNeu_higher = vec4NeuHigh;
  
  return 0;
}


int singletopAnalyser::CalcRecoCosStar() {
  TLorentzVector vec4LepRest = b_lep1, vec4JetRest = b_jet1;
  
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
  UInt_t i;
  
  // Begin of jobs IN GEN LEVEL
  
  if ( m_isMC && m_isSig ) {
    // If any one of them holds, then it is reaaaaally strange!
    if ( GetIdxGenTop() != 0 ) return -1;
    if ( GetIdxGenAssoQuark() != 0 ) return -1;
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
  }
  
  // End of jobs IN GEN LEVEL
  
  // Getting weights
  
  if ( m_isMC ) {
    Int_t nVtx = Pileup_nTrueInt;
    
    b_puweight = m_pileUp->getWeight(nVtx);
    b_puweight_up = m_pileUp->getWeight(nVtx, 1);
    b_puweight_dn = m_pileUp->getWeight(nVtx, -1);
    
    b_genweight = genWeight;
    b_weight = b_genweight * b_puweight;
    
    m_h1GenWeights->Fill(0.5, b_genweight);
    m_h1Weights->Fill(0.5, b_weight);
  } else {
    b_puweight = 1;
    b_genweight = 1;
    
    if ( !m_lumi->LumiCheck(run, luminosityBlock) ) return 0;
  }
  
  b_objstep = 1;
  
  // Checking pile-up configuration
  
  b_pvz = PV_z;
  b_npvs = PV_npvs;
  b_pv_ndof = PV_ndof;
  
  if ( fabs(PV_z) >= 24. ) return 0;
  if ( PV_npvs == 0 ) return 0;
  if ( PV_ndof < 4 ) return 0;
  
  b_objstep = 2;
  
  // Trigger
  
  b_tri = HLT_IsoTkMu24 || HLT_IsoMu24;
  
  Bool_t IsoMu24 = false;
  Bool_t IsoTkMu24 = false;
  
  for ( i = 0 ; i < nTrigObj ; i++ ) {
    if ( TrigObj_id[ i ] != 13 ) continue;
    if ( TrigObj_pt[ i ] < 24 ) continue;
    Int_t bits = TrigObj_filterBits[ i ];
    if ( bits & 0x2 ) IsoMu24 = true;
    if ( bits & 0x8 ) IsoTkMu24 = true;  
  }
  
  if ( !( IsoMu24 || IsoTkMu24 ) ) return 0;
  
  b_objstep = 3;
  
  // Gathering leptons
  
  auto muons = muonSelection();
  auto elecs = elecSelection();
  
  if ( muons.size() + elecs.size() < 1 ) return 0;
  
  double dPtMuon = ( muons.size() > 0 ? muons[ 0 ].Pt() : -1.0 );
  double dPtElec = ( elecs.size() > 0 ? elecs[ 0 ].Pt() : -1.0 );
  
  if ( muons.size() > 0 && dPtMuon > dPtElec ) {
    muons[ 0 ].Momentum(b_lep1);
    b_lep1_pid = muons[ 0 ].GetPdgCode();
    
    b_mueffweight = m_muonSF.getScaleFactor(muons[ 0 ], 13, 0);
    b_mueffweight_up = m_muonSF.getScaleFactor(muons[ 0 ], 13,  1);
    b_mueffweight_dn = m_muonSF.getScaleFactor(muons[ 0 ], 13, -1);
  } else {
    elecs[ 0 ].Momentum(b_lep1);
    b_lep1_pid = elecs[ 0 ].GetPdgCode();
  }
  
  // Gathering jets and b-jets
  
  vector<int> vecIdxJet;
  auto jets  = jetSelection(vecIdxJet);
  auto bjets = bjetSelection(vecIdxJet);
  
  if ( jets.size() <= 0 || bjets.size() <= 0 ) return 0;
  
  for ( i = 0 ; i < jets.size() ; i++ ) {
    if ( vecIdxJet[ i ] >= 0 ) break;
  }
  
  if ( jets.size() > 0 && i < jets.size() ) jets[ i ].Momentum(b_jet1);
  
  if ( bjets.size() > 0 ) bjets[ 0 ].Momentum(b_bjet1);
  if ( bjets.size() > 1 ) bjets[ 1 ].Momentum(b_bjet2);
  
  b_njets  = jets.size();
  b_nbjets = bjets.size();
  
  // Saving MET
  
  b_met     = MET_pt;
  b_met_phi = MET_phi;
  
  // Reconstruction of top quark
  
  RecoTop();
  if ( ( b_njets == 2 && b_nbjets == 1 ) || ( b_njets == 3 && b_nbjets == 2 ) )
    CalcRecoCosStar();
  
  b_objstep = 4;
  
  return 0;
}


void singletopAnalyser::Loop() {
  int nRes;
  
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //Prepare for new loop
    resetBranch();
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
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
  
  bool bIsMC;
  int nIdxStart = 0, nIdxEnd = 1048575 * 1024;
  
  if ( argc < 3 ) {
    printf("Usage: \n"
      "(on single core) singletopAnalyser [LIST_FILE_OF_FILES] [MC or RD] "
      "[idx of the first root file] [(idx+1) of the last root file]\n"
      "(for grid job (condor; in KISTI) python makejobs.py nanoAOD "
      "j`python batch_nanoAOD.py n [# OF FILES PER JOB]`\n");
    
    return 1;
  }
  
  if ( argc >= 4 ) nIdxStart = atoi(argv[ 3 ]);
  if ( argc >= 5 ) nIdxEnd = atoi(argv[ 4 ]);
  
  bIsMC = ( strcmp(argv[ 2 ], "RD") != 0 );
  
  std::ifstream fileList(argv[ 1 ]);
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
    
    TFile *fRoot = NULL;
    int nTrial = 0;
    
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
    
    singletopAnalyser t((TTree *)fRoot->Get("Events"), bIsMC, 
      strLine.find("ST_t-channel") != std::string::npos, true);
    
    std::string strOutput = "out";
    strOutput = strOutput + std::to_string(i) + ".root";
    t.setOutput(strOutput);
    
    t.Loop();
    
    //fRoot->Close(); // Don't do that; it's done in ~Events()
    cout << "End   : " << strLine << endl;
  }

  return 0;
}


