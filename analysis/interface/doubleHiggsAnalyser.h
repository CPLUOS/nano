#ifndef doubleHiggsAnalyser_H
#define doubleHiggsAnalyser_H
#include "nanoBase.h"
#include "nano/oxbridgekinetics/src/Mt2/Basic_Mt2_332_Calculator.h"
#include "nano/oxbridgekinetics/src/Mt2/ChengHanBisect_Mt2_332_Calculator.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/DataLoader.h"
#include "TMVA/MethodCuts.h"

class doubleHiggsAnalyser : public nanoBase {
private :
  // Output Variables
  TFile *out_file;
  TTree *out_tree;
  //// MT2 variables
  Float_t basic_MT2_332_bbll = -99;
  Float_t basic_MT2_332_blbl = -99;
  Float_t basic_MT2_332_b = -99;
  Float_t basic_MT2_332_l = -99;
  Float_t mT = -99;
  //// lepton variables
  TLorentzVector lepton1;
  TLorentzVector lepton2;
  TLorentzVector leptonlepton;
  Float_t ll_M = -99;
  Float_t ll_Pt = -99;
  Float_t ll_deltaR = -99;
  Float_t ll_deltaPhi = -99;
  std::map<Float_t, std::pair<int,int>, std::greater<Float_t>> leptons;
  std::map<Float_t, std::pair<int,int>, std::greater<Float_t>>::iterator lepton_iter;
  //// bottom variables
  TLorentzVector bottom1;
  TLorentzVector bottom2;
  TLorentzVector bottombottom;
  Float_t bb_M = -99;
  Float_t bb_Pt = -99;
  Float_t bb_deltaR = -99;
  Float_t bb_deltaPhi = -99;
  std::map<Float_t, int, std::greater<Float_t>> bottoms;
  std::map<Float_t, int, std::greater<Float_t>>::iterator bottom_iter;
  ////lepton and bottom variables
  TLorentzVector bottomlepton11;
  TLorentzVector bottomlepton12;
  TLorentzVector bottomlepton21;
  TLorentzVector bottomlepton22;
  TLorentzVector bbll;
  std::vector<Float_t> bl_deltaR;
  Float_t bl_min_deltaR = -99;
  Float_t bbll_deltaR = -99;
  Float_t bbll_deltaPhi = -99;
  //// MET variables
  TLorentzVector missing;
  //// cut variables
  Int_t step = 0;
  //// higgsness and topness
  Float_t higgsness = -99;
  Float_t topness = -99;
  //// tmva variables
  Float_t tmva_bdtg_output = -99;
  bool tmva_flag = false;
  
  // Instances
  //// TMVA reader
  TMVA::Reader* bdtg_reader;
  /// MT2 Calculators
  Mt2::Basic_Mt2_332_Calculator basic_mt2_332Calculator;
  Mt2::ChengHanBisect_Mt2_332_Calculator ch_bisect_mt2_332Calculator;

public :
  // before loop
  void MakeOutputBranch(TTree *tree);
  void SetOutput(TString output_file_name);
  void SetBranchAddress();
  void SetTMVA(TString weight_file_path);
  void Initiate(TString output_file_name);
  // during loop
  void ResetVariables();
  bool Analysis();
  // loop
  virtual void Loop();
  // after loop
  void Finalize();

  doubleHiggsAnalyser(TTree *tree=0, Bool_t isMC = false);
  ~doubleHiggsAnalyser();
};

// Constructors //////
doubleHiggsAnalyser::doubleHiggsAnalyser(TTree *tree, Bool_t isMC) : nanoBase(tree, isMC) {
  std::string env = getenv("CMSSW_BASE");
  //m_rocCor = new RoccoR(env+"/src/nano/analysis/data/rcdata.2016.v3/");
}

// deconstructors
doubleHiggsAnalyser::~doubleHiggsAnalyser() { }

#endif
