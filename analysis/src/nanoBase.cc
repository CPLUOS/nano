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
  
  m_jecUnc = NULL;
  m_rndEngine = NULL;
  
  if (isMC) {
    std::string env = getenv("CMSSW_BASE");
    
    std::string strPathJetResSFObj = env+"/src/nano/analysis/data/jer/"
      "Summer16_25nsV1_MC_SF_AK4PFchs.txt";
    std::string strPathJetResObj = env+"/src/nano/analysis/data/jer/"
      "Summer16_25nsV1_MC_PtResolution_AK4PFchs.txt";
    
    if (!exists_test(strPathJetResSFObj) || !exists_test(strPathJetResObj)) {
      std::cout << "Missing data file, run getFiles and try again" << std::endl;
      exit(50);
    }
    
    m_jetResSFObj = JMENano::JetResolutionScaleFactor(strPathJetResSFObj.c_str());
    m_jetResObj = JMENano::JetResolution(strPathJetResObj.c_str());
    
    m_rndEngine = new TRandom3(12345);
    
    std::string strPathJecUnc = env + "/src/nano/analysis/data/jer/"
      //"Summer16_23Sep2016V4_MC_Uncertainty_AK4PFchs.txt";
      "Summer16_23Sep2016V4_MC_UncertaintySources_AK4PFchs.txt";
    
    JetCorrectorParameters JetCorPar(strPathJecUnc, "Total");
    m_jecUnc = new JetCorrectionUncertainty(JetCorPar);
  }
  
  m_fDRcone_JER = 0.4; // For AK4 jets
  m_fResFactorMathcer = 3; // According to https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution
}

nanoBase::~nanoBase()
{
  if (m_output) {
    m_output->Write();
    m_output->Close();
  }
  
  if (m_jecUnc != NULL) delete m_jecUnc;
  if (m_rndEngine != NULL) delete m_rndEngine;
}

// In uncertainty study we need to switch the kinematic variables of jets
// In default, these are same as the original ones, but when a user wants to study systematic uncertainty 
// so that he/she needs to switch them to the evaluated ones, 
// just touching them in anlalyser class will be okay, and this is for it.
void nanoBase::GetJetMassPt(UInt_t nIdx, 
  Float_t &fJetMass, Float_t &fJetPt, Float_t &fJetEta, Float_t &fJetPhi, UInt_t unFlag) 
{
  Float_t fCorrFactor = 1.0;
  
  fJetMass = Jet_mass[nIdx];
  fJetPt   = Jet_pt[nIdx];
  fJetEta  = Jet_eta[nIdx];
  fJetPhi  = Jet_phi[nIdx];
  
  //if ( m_isMC && ( unFlag & ( OptFlag_JER_Up | OptFlag_JER_Dn | OptFlag_JES_Up | OptFlag_JES_Dn ) ) != 0 )
  if (m_isMC) {
    // Evaluating the central part cJER of the factor
    JMENano::JetParameters jetPars = {{JMENano::Binning::JetPt, fJetPt},
                                  {JMENano::Binning::JetEta, fJetEta},
                                  {JMENano::Binning::Rho, fixedGridRhoFastjetAll}};
    
    const double jetRes = m_jetResObj.getResolution(jetPars); // Note: this is relative resolution.
    const double cJER = m_jetResSFObj.getScaleFactor(jetPars, 
      ((unFlag & (OptFlag_JER_Up | OptFlag_JER_Dn)) == 0 ? Variation::NOMINAL : 
        ((unFlag & OptFlag_JER_Up) != 0 ? Variation::UP : Variation::DOWN)));
    
    // We need corresponding genJet
    Int_t nIdxGen = GetMatchGenJet(nIdx, fJetPt * jetRes);
    const double genJetPt = GenJet_pt[nIdxGen];
    
    if ( nIdxGen >= 0 ) {
      fCorrFactor = (genJetPt+(fJetPt-genJetPt)*cJER)/fJetPt;
    } else {
      const double smear = m_rndEngine->Gaus(0, 1);
      fCorrFactor = (cJER <= 1 ? 1 : 1+smear*jetRes*sqrt(cJER*cJER-1));
    }
    
    if ((unFlag & (OptFlag_JES_Up | OptFlag_JES_Dn)) != 0) { // JES
      // The evaluator needs pT and eta of the current jet
      m_jecUnc->setJetPt(fCorrFactor*fJetPt);
      m_jecUnc->setJetEta(fJetEta);
      
      if ((unFlag & OptFlag_JES_Up) != 0) {
        fCorrFactor *= 1+m_jecUnc->getUncertainty(true);
      } else {
        fCorrFactor *= 1-m_jecUnc->getUncertainty(true);
      }
    }
  }
  
  fJetMass *= fCorrFactor;
  fJetPt   *= fCorrFactor;
}


Int_t nanoBase::GetMatchGenJet(UInt_t nIdxJet, Float_t fResolution) {
  UInt_t i;
  
  double dEta, dPhi, dR;
  double dRFound = m_fDRcone_JER;
  UInt_t nIdxFound = 999;
  
  for (i = 0; i < nGenJet; i++) {
    dEta = Jet_eta[nIdxJet]-GenJet_eta[i];
    dPhi = std::abs(Jet_phi[nIdxJet]-GenJet_phi[i]);
    if (dPhi > (double)M_PI) dPhi -= (double)(2*M_PI);
    
    dR = std::sqrt(dEta*dEta+dPhi*dPhi);
    
    if (dR >= (double)m_fDRcone_JER*0.5) continue;
    if (dRFound > dR) {
      if (std::abs(Jet_pt[nIdxJet] - GenJet_pt[i]) >= m_fResFactorMathcer*fResolution) continue;
      
      dRFound = dR;
      nIdxFound = i;
    }
  }
 
  return (dRFound < 0.75*(double)m_fDRcone_JER ? (Int_t)nIdxFound : -1);
}

void nanoBase::Loop()
{
}
