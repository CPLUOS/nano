#include "nano/analysis/interface/topObjectSelection.h"

using std::vector;

topObjectSelection::topObjectSelection(TTree *tree, TTree *had, TTree *hadTruth, Bool_t isMC) :
  nanoBase(tree, had, hadTruth, isMC)
{}

vector<TParticle> topObjectSelection::elecSelection() {
  vector<TParticle> elecs; 
  for (UInt_t i = 0; i < nElectron; ++i){
    if (Electron_pt[i] < cut_ElectronPt) continue;
    if (std::abs(Electron_eta[i]) > cut_ElectronEta) continue;
    if (cut_ElectronIDType != NULL && cut_ElectronIDType[i] < cut_ElectronIDCut) continue; 
    float el_scEta = Electron_deltaEtaSC[i] + Electron_eta[i];
    if ( cut_ElectronSCEtaLower < std::abs(el_scEta) && std::abs(el_scEta) < cut_ElectronSCEtaUpper ) continue;
    if ( !additionalConditionForElectron(i) ) continue;
    
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
    if (cut_MuonIDType != NULL && !cut_MuonIDType[i]) continue;
    if (Muon_pt[i] < cut_MuonPt) continue;
    if (std::abs(Muon_eta[i]) > cut_MuonEta) continue;
    if (Muon_pfRelIso04_all[i] > cut_MuonRelIso04All) continue;
    if ( !additionalConditionForMuon(i) ) continue;
    
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
    if (Electron_pt[i] < cut_VetoElectronPt) continue;
    if (std::abs(Electron_eta[i]) > cut_VetoElectronEta) continue;
    if (cut_VetoElectronIDType != NULL && cut_VetoElectronIDType[i] < cut_VetoElectronIDCut) continue; 
    float el_scEta = Electron_deltaEtaSC[i] + Electron_eta[i];
    if ( cut_VetoElectronSCEtaLower < std::abs(el_scEta) && 
                                      std::abs(el_scEta) < cut_VetoElectronSCEtaUpper ) continue;
    if ( !additionalConditionForVetoElectron(i) ) continue;
    
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
    if (cut_VetoMuonIDType != NULL && !cut_VetoMuonIDType[i]) continue;
    if (Muon_pt[i] < cut_VetoMuonPt) continue;
    if (std::abs(Muon_eta[i]) > cut_VetoMuonEta) continue;
    if (Muon_pfRelIso04_all[i] > cut_VetoMuonRelIso04All) continue;
    if ( !additionalConditionForVetoMuon(i) ) continue;
    
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
    if (GenJet_pt[i] < cut_GenJetPt) continue;
    if (std::abs(GenJet_eta[i]) > cut_GenJetEta) continue;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(GenJet_pt[i], GenJet_eta[i], GenJet_phi[i], GenJet_mass[i]);
    bool hasOverLap = false;
    for (auto lep : recoleps){
        if (mom.TLorentzVector::DeltaR(lep) < cut_GenJetConeSizeOverlap) hasOverLap = true;
    }
    if (hasOverLap) continue;
    if ( !additionalConditionForGenJet(i) ) continue;
    
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
    if (Jet_pt[i] < cut_JetPt) continue;
    if (std::abs(Jet_eta[i]) > cut_JetEta) continue;
    if (Jet_jetId[i] < cut_JetID) continue;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
    bool hasOverLap = false;
    for (auto lep : recoleps){
        if (mom.TLorentzVector::DeltaR(lep) < cut_JetConeSizeOverlap) hasOverLap = true;
    }
    if (hasOverLap) continue;
    if ( !additionalConditionForJet(i) ) continue;
    
    auto jet = TParticle();
    jet.SetMomentum(mom);
    jet.SetFirstMother(i);
    jets.push_back(jet);
    if (m_isMC) {
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
    if (Jet_pt[i] < cut_BJetPt) continue;
    if (std::abs(Jet_eta[i]) > cut_BJetEta) continue;
    if (Jet_jetId[i] < cut_BJetID) continue;
    if (cut_BJetTypeBTag[i] < cut_BJetBTagCut) continue;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
    bool hasOverLap = false;
    for (auto lep : recoleps) {
      if (mom.TLorentzVector::DeltaR(lep) < cut_BJetConeSizeOverlap) hasOverLap = true;
    }
    if (hasOverLap) continue;
    if ( !additionalConditionForBJet(i) ) continue;
    
    auto bjet = TParticle();
    bjet.SetMomentum(mom);
    bjet.SetFirstMother(i);
    bjets.push_back(bjet);
  }
  return bjets;
}
