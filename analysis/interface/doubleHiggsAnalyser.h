#ifndef doubleHiggsAnalyser_H
#define doubleHiggsAnalyser_H
#include "nanoBase.h"
#include "nano/oxbridgekinetics/src/Mt2/Basic_Mt2_332_Calculator.h"
#include "nano/oxbridgekinetics/src/Mt2/ChengHanBisect_Mt2_332_Calculator.h"

class doubleHiggsAnalyser : public nanoBase {
private :
  Bool_t delphes_flag; 
  TTree *del_tree;

  TFile *out_file;
  TTree *out_tree;
  
  // Output Variables
  Float_t lester_MT2 = -99;
  Float_t basic_MT2_332 = -99;
  Float_t ch_bisect_MT2_332 = -99;
  Float_t lepton1_pt = -99;
  Float_t lepton2_pt = -99;
  Float_t missing_et = -99;
  Float_t mt = -99;

  Float_t muon1_mass = -99;
  Float_t muon2_mass = -99;
  Float_t muon_mass = -99;
  Float_t muon_px = -99;
  Float_t muon_py = -99;

  Int_t fromHiggs = 0;
  Int_t fromTop = 0;
  Int_t fromZ = 0;
  Int_t from1 = 0;
  Int_t from2 = 0;

  // TLorentzVectors
  TLorentzVector muon1;
  TLorentzVector muon2;
  TLorentzVector missing;

  // Delphes Variables
  TClonesArray *particles = 0;
  TClonesArray *missings = 0;
  TClonesArray *jets = 0;
  std::map<Float_t, int, std::greater<Float_t>> muons;

  // MT2 Calculators
  Mt2::Basic_Mt2_332_Calculator basic_mt2_332Calculator;
  Mt2::ChengHanBisect_Mt2_332_Calculator ch_bisect_mt2_332Calculator;

  double invis_mass     =  100; // GeV
  
public :
  // before loop settings
  void MakeOutputBranch(TTree *tree);
  void SetOutput(TString output_file_name);
  void SetBranchAddress();
  void SetNanoBranchAddress();
  void SetDelphesBranchAddress();
  void Initiate(TString output_file_name);
  // during loop
  void ResetVariables();
  void ResetNanoVariables();
  void ResetDelphesVariables();
  bool Analysis();
  bool NanoAnalysis();
  bool DelphesAnalysis();
    // calculating MT2
    double get_MT2();
    double get_Basic_332_MT2();
  // after loop
  void Finalize();
  virtual void Loop();

  doubleHiggsAnalyser(TTree *tree=0, Bool_t isMC = false);
  doubleHiggsAnalyser(TTree *tree=0, Bool_t isMC = false, Bool_t isDelphes = false);
  doubleHiggsAnalyser(TChain *tchain=0, Bool_t isMC = false, Bool_t isDelphes = false);
  ~doubleHiggsAnalyser();
};

// Constructors //////
doubleHiggsAnalyser::doubleHiggsAnalyser(TTree *tree, Bool_t isMC, Bool_t isDelphes) {
  if (!isDelphes) {
    throw std::invalid_argument("This constructor (tree, isMC, isDelphes) is for Delphes Object Analysis");
  }
  delphes_flag = isDelphes;
  del_tree = tree;
}

doubleHiggsAnalyser::doubleHiggsAnalyser(TTree *tree, Bool_t isMC) : nanoBase(tree, isMC) {
  std::string env = getenv("CMSSW_BASE");
  m_rocCor = new RoccoR(env+"/src/nano/analysis/data/rcdata.2016.v3/");
}

doubleHiggsAnalyser::doubleHiggsAnalyser(TChain *tchain, Bool_t isMC, Bool_t isDelphes) {
  if (!isDelphes) {
    throw std::invalid_argument("This constructor (tree, isMC, isDelphes) is for Delphes Object Analysis");
  }
  delphes_flag = isDelphes;
  TTree *tree = dynamic_cast<TTree*>(tchain);
  del_tree = tree;
}

// deconstructors
doubleHiggsAnalyser::~doubleHiggsAnalyser() { }

#endif
