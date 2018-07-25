#include "nano/analysis/interface/topEventSelectionDL.h"

using std::vector;

topEventSelectionDL::topEventSelectionDL(TTree *tree, TTree *had, TTree *hadTruth, Bool_t isMC, Bool_t dl, Bool_t sle, Bool_t slm) :
  topObjectSelection(tree, had, hadTruth, isMC, true, false),
  m_isDL(dl), m_isSL_e(sle), m_isSL_m(slm) {
}

topEventSelectionDL::~topEventSelectionDL() {
}

int topEventSelectionDL::EventSelection() {
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
  b_nvertex = PV_npvs;
  if (std::abs(PV_z) >= 24.) return b_step;
  if (PV_npvs == 0) return b_step;
  if (PV_ndof < 4) return b_step;

  h_cutFlow->Fill(2);

  auto muons = muonSelection();
  auto elecs = elecSelection();

  if (muons.size() + elecs.size() != 2) return b_step;
  b_step = 0;
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

  b_lep1_pid = recolep1.GetPdgCode();
  b_lep2_pid = recolep2.GetPdgCode();

  b_lep1_idx = recolep1.GetFirstMother();
  b_lep2_idx = recolep2.GetFirstMother();

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
  b_step = 1;
  h_cutFlow->Fill(4);

  if (b_channel == CH_MUEL || b_dilep.M() < 76 || b_dilep.M() > 106) {
    b_step = 2;
    h_cutFlow->Fill(5);
  }

  b_met = MET_pt;

  if (b_channel == CH_MUEL || b_met > 40) {
    if (b_step == 2) {
      ++b_step;
      h_cutFlow->Fill(6);
    }
  }

  auto jets = jetSelection();
  b_njet = jets.size();

  if (b_njet >= 2) {
    if (b_step == 3) {
      ++b_step;
      h_cutFlow->Fill(7);
    }
  }

  auto bjets = bjetSelection();
  b_nbjet = bjets.size();

  if (b_nbjet > 0) {
    if (b_step == 4) {
      ++b_step;
      h_cutFlow->Fill(8);
    }
  }

  if (nhad < 1) return 0;
  else return b_step;
}


void topEventSelectionDL::Reset() {

  recolep1.Clear(); recolep2.Clear();
  b_lep1.SetPtEtaPhiM(0,0,0,0); b_lep2.SetPtEtaPhiM(0,0,0,0); b_dilep.SetPtEtaPhiM(0,0,0,0); b_jet1.SetPtEtaPhiM(0,0,0,0); b_jet2.SetPtEtaPhiM(0,0,0,0);

  b_lep1_pid = 0; b_lep2_pid = 0; b_lep1_idx = -1; b_lep2_idx = -1;
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
