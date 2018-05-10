#ifndef topAnalysis_H
#define topAnalysis_H

#include "nano/analysis/src/nanoAnalysis.h"
#include "TopTriggerSF.h"
//#include "TTbarModeDefs.h"

class topAnalysis : public nanoAnalysis
{
protected:
  std::vector<Float_t> b_csvweights;
  float b_btagweight;
  //TParticle GetTParticle(int pdgId, int idx);

public:
  enum TTLLChannel { CH_NOLL = 0, CH_MUEL, CH_ELEL, CH_MUMU };

  //Bool_t lumiCheck();
  std::vector<TParticle> muonSelection();
  std::vector<TParticle> elecSelection();
  std::vector<TLorentzVector> recoleps;
  std::vector<TParticle> jetSelection();
  std::vector<TParticle> bjetSelection();

public:
  Bool_t m_isDL, m_isSL_e, m_isSL_m;

//  void setOutput(std::string outputName);
//  void LoadModules(pileUpTool* pileUp, lumiTool* lumi);
//  void collectTMVAvalues();
  topAnalysis(TTree *tree=0, Bool_t isMC = false, Bool_t dl = false, Bool_t sle = false, Bool_t slm = false);
  ~topAnalysis();  
  virtual void Loop() = 0;
};

topAnalysis::topAnalysis(TTree *tree, Bool_t isMC, Bool_t dl, Bool_t sle, Bool_t slm) : nanoAnalysis(tree, isMC), m_isDL(dl), m_isSL_e(sle), m_isSL_m(slm)// Events(tree, isMC), m_isDL(dl), m_isSL_e(sle), m_isSL_m(slm)
{
}

topAnalysis::~topAnalysis()
{
}

vector<TParticle> topAnalysis::muonSelection()
{
  vector<TParticle> muons; 
  for (UInt_t i = 0; i < nMuon; ++i){
    if (!Muon_tightId[i]) continue;
    if (Muon_pt[i] < 20) continue;
    if (std::abs(Muon_eta[i]) > 2.4) continue;
    if (Muon_pfRelIso04_all[i] > 0.15) continue;
    if (!Muon_globalMu[i]) continue;
    if (!Muon_isPFcand[i]) continue;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]);
    auto muon = TParticle();
    muon.SetPdgCode(13*Muon_charge[i]*-1);
    muon.SetMomentum(mom);

    muons.push_back(muon);
  }
  return muons;
}


vector<TParticle> topAnalysis::elecSelection()
{
  vector<TParticle> elecs; 
  for (UInt_t i = 0; i < nElectron; ++i){
    if (Electron_pt[i] < 20) continue;
    if (std::abs(Electron_eta[i]) > 2.4) continue;
    if (Electron_cutBased[i] < 3) continue; 
    float el_scEta = Electron_deltaEtaSC[i] + Electron_eta[i];
    if ( std::abs(el_scEta) > 1.4442 &&  std::abs(el_scEta) < 1.566 ) continue;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Electron_pt[i], Electron_eta[i], Electron_phi[i], Electron_mass[i]);

    auto elec = TParticle();
    elec.SetPdgCode(11*Electron_charge[i]*-1);
    elec.SetMomentum(mom);
    elec.SetWeight(el_scEta);
    elecs.push_back(elec);
  }
  return elecs;
}

vector<TParticle> topAnalysis::jetSelection()
{
  vector<TParticle> jets;
  float Jet_SF_CSV[19] = {1.0,};
  for (UInt_t i = 0; i < nJet; ++i){
    if (Jet_pt[i] < 30) continue;
    if (std::abs(Jet_eta[i]) > 2.4) continue;
    if (Jet_jetId[i] < 1) continue;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
    bool hasOverLap = false;
    for (auto lep : recoleps){
        if (mom.TLorentzVector::DeltaR(lep) < 0.4) hasOverLap = true;
    }
    if (hasOverLap) continue;
    auto jet = TParticle();
    jet.SetMomentum(mom);
    jets.push_back(jet);
    for (UInt_t iu = 0; iu < 19; iu++) {
     // Jet_SF_CSV[iu] *= m_btagSF.getSF(jet, Jet_btagCSVV2[i], Jet_hadronFlavour[i], iu);
    }
  }
  for (UInt_t i =0; i<19; i++) b_csvweights.push_back(Jet_SF_CSV[i]);
  b_btagweight = Jet_SF_CSV[0];
  
  return jets;
}

vector<TParticle> topAnalysis::bjetSelection()
{
  vector<TParticle> bjets;
  for (UInt_t i = 0; i < nJet; ++i ) {
    if (Jet_pt[i] < 30) continue;
    if (std::abs(Jet_eta[i]) > 2.4) continue;
    if (Jet_jetId[i] < 1) continue;
    if (Jet_btagCSVV2[i] < 0.8484) continue;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
    bool hasOverLap = false;
    for (auto lep : recoleps) {
      if (mom.TLorentzVector::DeltaR(lep) < 0.4) hasOverLap = true;
    }
    if (hasOverLap) continue;
    auto bjet = TParticle();
    bjet.SetMomentum(mom);
    bjets.push_back(bjet);
  }
  return bjets;
}

#endif
