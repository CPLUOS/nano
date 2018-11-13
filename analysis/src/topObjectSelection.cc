#include "nano/analysis/interface/topObjectSelection.h"

using std::vector;

topObjectSelection::topObjectSelection(TTree *tree, TTree *had, TTree *hadTruth, Bool_t isMC, Bool_t _isDilep, Bool_t _isSemiLep) :
  nanoBase(tree, had, hadTruth, isMC),
  isDilep(_isDilep),
  isSemiLep(_isSemiLep)
{}

vector<TParticle> topObjectSelection::elecSelection() {
  vector<TParticle> elecs; 
  for (UInt_t i = 0; i < nElectron; ++i){
    if (Electron_pt[i] < 20) continue;
    if (isSemiLep) { if (Electron_pt[i] < 30) continue; }
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
    elec.SetFirstMother(i);
    elecs.push_back(elec);
  }
  return elecs;
}

vector<TParticle> topObjectSelection::muonSelection() {
  vector<TParticle> muons; 
  for (UInt_t i = 0; i < nMuon; ++i){
    if (!Muon_tightId[i]) continue;
    if (Muon_pt[i] < 20) continue;
    if (isSemiLep) { if (Muon_pt[i] < 26) continue; }
    if (std::abs(Muon_eta[i]) > 2.4) continue;
    if (isSemiLep) { if (std::abs(Muon_eta[i]) > 2.1) continue; }
    if (Muon_pfRelIso04_all[i] > 0.15) continue;
    if (!Muon_globalMu[i]) continue;
    if (!Muon_isPFcand[i]) continue;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]);
    auto muon = TParticle();
    muon.SetPdgCode(13*Muon_charge[i]*-1);
    muon.SetMomentum(mom);
    muon.SetFirstMother(i);
    muons.push_back(muon);
  }
  return muons;
}

/// TODO: Implement veto electron selection (currently empty on TTbarXSecSynchronization)
vector<TParticle> topObjectSelection::vetoElecSelection() {
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
    elec.SetFirstMother(i);
    elecs.push_back(elec);
  }
  return elecs;
}

vector<TParticle> topObjectSelection::vetoMuonSelection() {
  vector<TParticle> muons; 
  for (UInt_t i = 0; i < nMuon; ++i){
    // if (!Muon_looseId[i]) continue;

    // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Loose_Muon
    if (!Muon_isPFcand[i]) continue;
    if (!(Muon_globalMu[i] || Muon_trackerMu[i])) continue;

    
    if (Muon_pt[i] < 10) continue;
    if (std::abs(Muon_eta[i]) > 2.4) continue;
    if (Muon_pfRelIso04_all[i] > 0.25) continue;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]);
    auto muon = TParticle();
    muon.SetPdgCode(13*Muon_charge[i]*-1);
    muon.SetMomentum(mom);
    muon.SetFirstMother(i);
    muons.push_back(muon);
  }
  return muons;
}

vector<TParticle> topObjectSelection::genJetSelection() {
  vector<TParticle> jets;
  for (UInt_t i = 0; i < nJet; ++i){
    if (GenJet_pt[i] < 30) continue;
    if (std::abs(GenJet_eta[i]) > 2.4) continue;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(GenJet_pt[i], GenJet_eta[i], GenJet_phi[i], GenJet_mass[i]);
    bool hasOverLap = false;
    for (auto lep : recoleps){
        if (mom.TLorentzVector::DeltaR(lep) < 0.4) hasOverLap = true;
    }
    if (hasOverLap) continue;
    auto jet = TParticle();
    jet.SetMomentum(mom);
    jet.SetFirstMother(i);
    jets.push_back(jet);
  }
  return jets;
}

vector<TParticle> topObjectSelection::jetSelection() {
  double csvWgtHF = 1.0, csvWgtLF = 1.0, csvWgtC = 1.0, csvWgtTotal = 1.0;
  double csvWgtHF_up = 1.0, csvWgtLF_up = 1.0, csvWgtC_up = 1.0, csvWgtTotal_up = 1.0;
  double csvWgtHF_dn = 1.0, csvWgtLF_dn = 1.0, csvWgtC_dn = 1.0, csvWgtTotal_dn = 1.0;
  vector<TParticle> jets;
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
    jet.SetFirstMother(i);
    jets.push_back(jet);
    if (m_isMC) {
      b_btagCSVV2 = Jet_btagCSVV2[i];
      BTagEntry::JetFlavor JF = BTagEntry::FLAV_UDSG;
      if (abs(Jet_hadronFlavour[i]) == 5) {
        JF = BTagEntry::FLAV_B;
        double iCSVWgtHF = m_btagSF.eval_auto_bounds("central", JF, abs(Jet_eta[i]), Jet_pt[i], Jet_btagCSVV2[i]);
        double iCSVWgtHF_up = m_btagSF_up.eval_auto_bounds("up", JF, abs(Jet_eta[i]), Jet_pt[i], Jet_btagCSVV2[i]);
        double iCSVWgtHF_dn = m_btagSF_dn.eval_auto_bounds("down", JF, abs(Jet_eta[i]), Jet_pt[i], Jet_btagCSVV2[i]);
        if (iCSVWgtHF != 0) csvWgtHF *= iCSVWgtHF;
        if (iCSVWgtHF_up != 0) csvWgtHF_up *= iCSVWgtHF_up;
        if (iCSVWgtHF_dn != 0) csvWgtHF_dn *= iCSVWgtHF_dn;
      }
      else if (abs(Jet_hadronFlavour[i]) == 4) {
        JF = BTagEntry::FLAV_C;
        double iCSVWgtC = m_btagSF.eval_auto_bounds("central", JF, abs(Jet_eta[i]), Jet_pt[i], Jet_btagCSVV2[i]);
        double iCSVWgtC_up = m_btagSF_up.eval_auto_bounds("up", JF, abs(Jet_eta[i]), Jet_pt[i], Jet_btagCSVV2[i]);
        double iCSVWgtC_dn = m_btagSF_dn.eval_auto_bounds("down", JF, abs(Jet_eta[i]), Jet_pt[i], Jet_btagCSVV2[i]);
        if (iCSVWgtC != 0) csvWgtC *= iCSVWgtC;
        if (iCSVWgtC_up != 0) csvWgtC_up *= iCSVWgtC_up;
        if (iCSVWgtC_dn != 0) csvWgtC_dn *= iCSVWgtC_dn;
      }
      else {
        JF = BTagEntry::FLAV_UDSG;
        double iCSVWgtLF = m_btagSF.eval_auto_bounds("central", JF, abs(Jet_eta[i]), Jet_pt[i], Jet_btagCSVV2[i]);
        double iCSVWgtLF_up = m_btagSF_up.eval_auto_bounds("up", JF, abs(Jet_eta[i]), Jet_pt[i], Jet_btagCSVV2[i]);
        double iCSVWgtLF_dn = m_btagSF_dn.eval_auto_bounds("down", JF, abs(Jet_eta[i]), Jet_pt[i], Jet_btagCSVV2[i]);
        if (iCSVWgtLF != 0) csvWgtLF *= iCSVWgtLF;
        if (iCSVWgtLF_up != 0) csvWgtLF_up *= iCSVWgtLF_up;
        if (iCSVWgtLF_dn != 0) csvWgtLF_dn *= iCSVWgtLF_dn;
      } 
    }
  }
  csvWgtTotal = csvWgtHF * csvWgtC * csvWgtLF;  
  csvWgtTotal_up = csvWgtHF_up * csvWgtC_up * csvWgtLF_up;  
  csvWgtTotal_dn = csvWgtHF_dn * csvWgtC_dn * csvWgtLF_dn;  

  b_btagweight = csvWgtTotal;
  b_btagweight_up = csvWgtTotal_up;
  b_btagweight_dn = csvWgtTotal_dn;
  
  return jets;
}

vector<TParticle> topObjectSelection::bjetSelection() {
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
    bjet.SetFirstMother(i);
    bjets.push_back(bjet);
  }
  return bjets;
}
