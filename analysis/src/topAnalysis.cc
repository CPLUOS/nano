#include "nano/analysis/interface/topAnalysis.h"

using std::vector;

topAnalysis::topAnalysis(TTree *tree, Bool_t isMC, Bool_t dl, Bool_t sle, Bool_t slm) : nanoAnalysis(tree, isMC), m_isDL(dl), m_isSL_e(sle), m_isSL_m(slm) {
}

topAnalysis::~topAnalysis() {
}

int topAnalysis::EventSelection() {
  h_cutFlow->Fill(0);

  //Run for MC
  if (m_isMC) {
    Int_t nvtx = Pileup_nTrueInt;
    b_puweight = m_pileUp->getWeight(nvtx);

    b_genweight = genWeight;
    h_genweights->Fill(0.5, b_genweight);
    b_weight = b_genweight * b_puweight;
  }
  else {
    b_puweight = 1;
    b_genweight = 0;
    if (!(m_lumi->LumiCheck(run, luminosityBlock))) return b_step;
  }

  h_nevents->Fill(0.5, b_genweight*b_puweight);

  h_cutFlow->Fill(1);

  if (std::abs(PV_z) >= 24.) return b_step;
  if (PV_npvs == 0) return b_step;
  if (PV_ndof < 4) return b_step;

  h_cutFlow->Fill(2);

  auto muons = muonSelection();
  auto elecs = elecSelection();

  if (muons.size() + elecs.size() != 2) return b_step;

  h_cutFlow->Fill(3);

  int mulpdg = 1;
  if (muons.size() == 2) {
      recolep1 = muons[0];
      recolep2 = muons[1];
      mulpdg = muons[0].GetPdgCode()*muons[1].GetPdgCode();
      b_channel = CH_MUMU;
  } else if (muons.size() == 1 && elecs.size() == 1) {
      recolep1 = muons[0];
      recolep2 = elecs[0];
      mulpdg = muons[0].GetPdgCode()*elecs[0].GetPdgCode();
      b_channel = CH_MUEL;
  } else if (elecs.size() == 2) {
      recolep1 = elecs[0];
      recolep2 = elecs[1];
      mulpdg = elecs[0].GetPdgCode()*elecs[1].GetPdgCode();
      b_channel = CH_ELEL;
  }

  recolep1.Momentum(b_lep1);
  recolep2.Momentum(b_lep2);

  recoleps.push_back(b_lep1);
  recoleps.push_back(b_lep2);

  b_dilep = b_lep1 + b_lep2;

  //Triggers
  b_trig_m = HLT_IsoTkMu24 || HLT_IsoMu24;
  b_trig_e = HLT_Ele27_WPTight_Gsf;
  b_trig_mm = HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL
    || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;
  b_trig_em = HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL
    || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
  b_trig_ee = HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;

  if (b_channel == CH_MUMU) {
    if (m_isMC) {
      if (!(b_trig_mm || b_trig_m)) return b_step;
    } else if (m_isDL) {
      if (!(b_trig_mm)) return b_step;
    } else if (m_isSL_m) {
      if (b_trig_mm||!b_trig_m) return b_step;
    }
  }

  if (b_channel == CH_MUEL) {
    if (m_isMC) {
      if (!(b_trig_em || b_trig_m || b_trig_e)) return b_step;
    } else if (m_isDL) {
      if (!(b_trig_em)) return b_step;
    } else if (m_isSL_e) {
      if (b_trig_em || !b_trig_e || b_trig_m) return b_step;
    } else if (m_isSL_m) {
      if (b_trig_em || b_trig_e || !b_trig_m) return b_step;
    }
  }

  if (b_channel == CH_ELEL) {
    if (m_isMC) {
      if (!(b_trig_ee || b_trig_e)) return b_step;
    } else if (m_isDL) {
      if (!b_trig_ee) return b_step;
    } else if (m_isSL_e) {
      if (b_trig_ee || !b_trig_e) return b_step;
    }
  }

  //leptonS
  b_mueffweight    = muonSF_.getScaleFactor(recolep1, 13, 0)*muonSF_.getScaleFactor(recolep2, 13,  0);
  b_mueffweight_up = muonSF_.getScaleFactor(recolep1, 13, 1)*muonSF_.getScaleFactor(recolep2, 13, 1);
  b_mueffweight_dn = muonSF_.getScaleFactor(recolep1, 13, -1)*muonSF_.getScaleFactor(recolep2, 13, -1);

  b_eleffweight    = elecSF_.getScaleFactor(recolep1, 11, 0)*elecSF_.getScaleFactor(recolep2, 11,  0);
  b_eleffweight_up = elecSF_.getScaleFactor(recolep1, 11, 1)*elecSF_.getScaleFactor(recolep2, 11, 1);
  b_eleffweight_dn = elecSF_.getScaleFactor(recolep1, 11, -1)*elecSF_.getScaleFactor(recolep2, 11, -1);

  b_tri = b_tri_up = b_tri_dn = 0;
  b_tri = computeTrigSF(recolep1, recolep2);
  b_tri_up = computeTrigSF(recolep1, recolep2, 1);
  b_tri_dn = computeTrigSF(recolep1, recolep2, -1);

  if (b_dilep.M() < 20. || mulpdg > 0) return b_step;
  b_step1 = true;
  b_step = 1;
  h_cutFlow->Fill(4);

  if (b_channel == CH_MUEL || b_dilep.M() < 76 || b_dilep.M() > 106) {
    b_step2 = true;
    b_step = 2;
    h_cutFlow->Fill(5);
  }

  b_met = MET_pt;

  if (b_channel == CH_MUEL || b_met > 40) {
    b_step3 = true;
    if (b_step == 2) {
      ++b_step;
      h_cutFlow->Fill(6);
    }
  }

  auto jets = jetSelection();
  b_njet = jets.size();

  if (b_njet >= 2) {
    b_step4 = true;
    if (b_step == 3) {
      ++b_step;
      h_cutFlow->Fill(7);
    }
  }

  auto bjets = bjetSelection();
  b_nbjet = bjets.size();

  if (b_nbjet > 0) {
    b_step5 = true;
    if (b_step == 4) {
      ++b_step;
      h_cutFlow->Fill(8);
    }
  }

  if (nhad < 1) return 0;

  TLorentzVector vecSumDMLep1, vecSumDMLep2;
  float fMDMLep1, fMDMLep2;
  float fDeltaEta, fDeltaPhi;
  float fSqrtdRMLep1, fSqrtdRMLep2;

  for (UInt_t i = 0; i < nhad; ++i) {
    TLorentzVector d0_tlv;
    d0_tlv.SetPtEtaPhiM(had_pt[i], had_eta[i], had_phi[i], had_mass[i]);
    d0s.push_back(d0_tlv);
    if (d0s.size() < 1) continue;
    sort(d0s.begin(), d0s.end(), [](const TLorentzVector& a, const TLorentzVector& b){return a.Pt() > b.Pt();});
    d0s.erase(d0s.begin()+1, d0s.end());
    b_d0 = d0s[0];

    vecSumDMLep1 = b_lep1 + b_d0;
    vecSumDMLep2 = b_lep2 + b_d0;
    fMDMLep1 = vecSumDMLep1.M();
    fMDMLep2 = vecSumDMLep2.M();
    fDeltaEta = b_lep1.Eta() - b_d0.Eta();
    fDeltaPhi = b_lep1.Phi() - b_d0.Phi();
    fSqrtdRMLep1 = fDeltaEta * fDeltaEta + fDeltaPhi * fDeltaPhi;

    fDeltaEta = b_lep2.Eta() - b_d0.Eta();
    fDeltaPhi = b_lep2.Phi() - b_d0.Phi();
    fSqrtdRMLep2 = fDeltaEta * fDeltaEta + fDeltaPhi * fDeltaPhi;

    b_d0_lepSV_lowM.push_back(( fMDMLep1 >= fMDMLep2 ? fMDMLep1 : fMDMLep2 ));
    b_d0_lepSV_dRM.push_back(( fSqrtdRMLep1 >= fSqrtdRMLep2 ? fMDMLep1 : fMDMLep2 ));
  }
  return b_step;
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

