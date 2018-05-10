#define Events_cxx
<<<<<<< HEAD
#ifndef nanoAnalysis_H
#define nanoAnalysis_H

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
#include "BTagCalibrationStandalone.cc"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/DataLoader.h"
#include "TMVA/MethodCuts.h"


class nanoAnalysis : public Events
{
public:
  nanoAnalysis(TTree *tree=0, Bool_t isMC=false);
  virtual ~nanoAnalysis();
  virtual void Loop() = 0;

  //Output Variables
  TFile* m_output;
  //Tree
  TTree* m_tree;

  pileUpTool* m_pileUp;
  lumiTool* m_lumi;
  RoccoR* m_rocCor;
  MuonScaleFactorEvaluator m_muonSF;
  ElecScaleFactorEvaluator m_elecSF;
  BTagCalibrationReader m_btagSF;

  Bool_t m_isMC;
};

nanoAnalysis::nanoAnalysis(TTree *tree, Bool_t isMC) : Events(tree), m_isMC(isMC)
{
  m_pileUp = new pileUpTool();
  string env = getenv("CMSSW_BASE");
  string username = getenv("USER");
  m_lumi = new lumiTool(env+"/src/nano/analysis/data/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt");
  string csvFileName = "CSVv2_Moriond17_B_H.csv";
  std::string csvFile = env+"/src/nano/analysis/data/btagSF/"+csvFileName;
  BTagCalibration calib("csvv2", csvFile);
  m_btagSF = BTagCalibrationReader(BTagEntry::OP_MEDIUM,"central",{"up","down"});
  m_btagSF.load(calib, BTagEntry::FLAV_B, "mujets");
}

nanoAnalysis::~nanoAnalysis()
{
 m_output->Write();
 m_output->Close();
}

void nanoAnalysis::Loop(){};

#endif
