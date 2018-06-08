#include "nano/analysis/interface/nanoBase.h"

using std::string;

nanoBase::nanoBase(TTree *tree, TTree *had, TTree *hadTruth, Bool_t isMC) : Events(tree, had, hadTruth), m_isMC(isMC) {
  m_pileUp = new pileUpTool();
  string env = getenv("CMSSW_BASE");
  m_lumi = new lumiTool(env+"/src/nano/analysis/data/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt");
  string csvFileName = "CSVv2_Moriond17_B_H.csv";
  std::string csvFile = env+"/src/nano/analysis/data/btagSF/"+csvFileName;
  BTagCalibration calib("csvv2", csvFile);
  m_btagSF = BTagCalibrationReader(BTagEntry::OP_MEDIUM,"central",{"up","down"});
  m_btagSF.load(calib, BTagEntry::FLAV_B, "mujets");
}

nanoBase::~nanoBase() {
 m_output->Write();
 m_output->Close();
}

void nanoBase::Loop(){};
