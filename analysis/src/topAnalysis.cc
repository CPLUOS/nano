#include "nano/analysis/interface/topAnalysis.h"

using std::vector;

topAnalysis::topAnalysis(TTree *tree, Bool_t isMC) : nanoAnalysis(tree, isMC) {}

vector<TParticle> topAnalysis::elecSelection() {
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

vector<TParticle> topAnalysis::muonSelection() {
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

/// TODO: Implement veto electron selection (currently empty on TTbarXSecSynchronization)
vector<TParticle> topAnalysis::vetoElecSelection() {
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

vector<TParticle> topAnalysis::vetoMuonSelection() {
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

    muons.push_back(muon);
  }
  return muons;
}

vector<TParticle> topAnalysis::jetSelection() {
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

vector<TParticle> topAnalysis::bjetSelection() {
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

void topAnalysis::Reset() {

  recolep1.Clear(); recolep2.Clear();
  b_lep1.SetPtEtaPhiM(0,0,0,0); b_lep2.SetPtEtaPhiM(0,0,0,0); b_dilep.SetPtEtaPhiM(0,0,0,0); b_jet1.SetPtEtaPhiM(0,0,0,0); b_jet2.SetPtEtaPhiM(0,0,0,0);

  b_lep1_pid = 0; b_lep2_pid = 0;
  b_jet1_CSVInclV2 = -1; b_jet2_CSVInclV2 = -1;
  b_csvweights.clear();

  b_nvertex = -1; b_step = -1; b_channel = 0; b_njet = -1; b_nbjet = -1;
  b_met = -9; b_weight = 1; b_genweight = 1; b_puweight = 1; b_btagweight = 1;
  b_mueffweight = 1;b_mueffweight_up = 1;b_mueffweight_dn = 1;
  b_eleffweight = 1;b_eleffweight_up = 1;b_eleffweight_dn = 1;

  b_tri = 0; b_tri_up = 0; b_tri_dn = 0;
  b_trig_m = false; b_trig_m2 = false; b_trig_e = false; b_trig_mm = false; b_trig_em = false; b_trig_ee = false;

  recoleps.clear();
}
