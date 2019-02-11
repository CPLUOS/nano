#include "nano/analysis/interface/nanoBase.h"
#include <string>
#include <fstream>

using std::string;

inline bool exists_test (string& name) {
  std::ifstream f(name.c_str());
  return f.good();
}

nanoBase::nanoBase(TTree *tree, TTree *had, TTree *hadTruth, Bool_t isMC) :
  Events(tree, had, hadTruth),
  m_output(0),
  m_tree(0),
  m_isMC(isMC)
{
  m_pileUp = new pileUpTool();
  string env = getenv("CMSSW_BASE");
  string lumi = env+"/src/nano/analysis/data/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt";

  if (!exists_test(lumi)) {
     std::cout << "Missing data file, run getFiles and try again" << std::endl;
     exit(50);
  }
  m_lumi = new lumiTool(lumi);
  string csvFileName = "CSVv2_Moriond17_B_H.csv";
  std::string csvFile = env+"/src/nano/analysis/data/btagSF/"+csvFileName;
  BTagCalibration calib("csvv2", csvFile);
  m_btagSF = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central");
  m_btagSF_up = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "up");
  m_btagSF_dn = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "down");
  m_btagSF.load(calib, BTagEntry::FLAV_B, "comb");
  m_btagSF_up.load(calib, BTagEntry::FLAV_B, "comb");
  m_btagSF_dn.load(calib, BTagEntry::FLAV_B, "comb");
  
  string strBaseData = env+"/src/nano/analysis/data/scaleFactor/2016/";
  m_strTrigSFEl = strBaseData+"HLT_Ele32_eta2p1_WPTight_Gsf_FullRunRange.root";
  m_strTrigSFMu = strBaseData+"EfficienciesAndSF_Run2016.root";
  m_strLeptonSFEl = strBaseData+"2016LegacyReReco_ElectronTight.root";
  m_strLeptonSFMu = strBaseData+"Run2016_SF_ID.root";
  
  if (!exists_test(m_strTrigSFEl)   || !exists_test(m_strTrigSFMu) || 
      !exists_test(m_strLeptonSFEl) || !exists_test(m_strLeptonSFMu))
  {
    std::cout << "Missing data file, run getFiles and try again" << std::endl;
    exit(50);
  }
  
  m_tableTrigSFEl.LoadData(m_strTrigSFEl, "SF");
  m_tableTrigSFMu.LoadData(m_strTrigSFMu, "IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio");
  
  m_tableLeptonSFEl.LoadData(m_strLeptonSFEl, "EGamma_SF2D");
  m_tableLeptonSFMu.LoadData(m_strLeptonSFMu, "NUM_TightID_DEN_genTracks_eta_pt");
}

nanoBase::~nanoBase()
{
  if (m_output) {
    m_output->Write();
    m_output->Close();
  }
}

void nanoBase::Loop()
{
}
