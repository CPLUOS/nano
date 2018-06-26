#include "nano/analysis/interface/topEventSelectionSL.h"

using std::vector;

topEventSelectionSL::topEventSelectionSL(TTree *tree, TTree *had, TTree *hadTruth, Bool_t isMC, Bool_t sle, Bool_t slm) :
  topObjectSelection(tree, had, hadTruth, isMC, false, true),
  h_nevents(0),
  h_genweights(0),
  h_cutFlow(0),
  h_cutFlowEl(0),
  h_cutFlowMu(0),
  m_isSL_e(sle),
  m_isSL_m(slm)
{
}

topEventSelectionSL::~topEventSelectionSL()
{
}

void topEventSelectionSL::Reset()
{
  b_step = 0;

  b_channel = -9; 
  b_njet = -9; b_nvertex = -9;

  b_lep.SetPtEtaPhiM(0,0,0,0); b_lep_pid = 0;

  b_met = -9; b_weight = 1; b_genweight = 1; b_puweight = 1; b_btagweight = 1;
  b_mueffweight = 1;b_mueffweight_up = 1;b_mueffweight_dn = 1;
  b_eleffweight = 1;b_eleffweight_up = 1;b_eleffweight_dn = 1;
  b_tri = 0;

  recoleps.clear();
  b_csvweights.clear();
}

int topEventSelectionSL::EventSelection()
{
  b_step = 0;
  if (h_cutFlow) h_cutFlow->Fill(0);
  if (h_cutFlowEl) h_cutFlowEl->Fill(0);
  if (h_cutFlowMu) h_cutFlowMu->Fill(0);

  //Run for MC
  if (m_isMC) {
    Int_t nvtx = Pileup_nTrueInt;
    b_puweight = m_pileUp->getWeight(nvtx);

    b_genweight = genWeight;
    if (h_genweights) h_genweights->Fill(0.5, b_genweight);
    b_weight = b_genweight * b_puweight;
  } else {
    b_puweight = 1;
    b_genweight = 0;
    if (!(m_lumi->LumiCheck(run, luminosityBlock))) return b_step;
  }

  if (h_nevents) h_nevents->Fill(0.5, b_genweight*b_puweight);

  if (h_cutFlow) h_cutFlow->Fill(1);

  if (std::abs(PV_z) >= 24.) return b_step;
  if (PV_npvs == 0) return b_step;
  if (PV_ndof < 4) return b_step;

  if (h_cutFlow) h_cutFlow->Fill(2);

  //Triggers
  b_trig_m = HLT_IsoTkMu24 || HLT_IsoMu24;
  b_trig_e = HLT_Ele27_WPTight_Gsf;

  // TODO Check trigger requirements (TTbarXSecSynchronization page doesn't have yet)
  
  // if (b_channel == CH_MU) {
  //   if (!b_trig_m) return b_step;
  // }

  // if (b_channel == CH_EL) {
  //   if (!b_trig_e) return b_step;
  // }

  //leptonS
  b_mueffweight    = muonSF_.getScaleFactor(recolep, 13, 0);
  b_mueffweight_up = muonSF_.getScaleFactor(recolep, 13, 1);
  b_mueffweight_dn = muonSF_.getScaleFactor(recolep, 13, -1);

  b_eleffweight    = elecSF_.getScaleFactor(recolep, 11, 0);
  b_eleffweight_up = elecSF_.getScaleFactor(recolep, 11, 1);
  b_eleffweight_dn = elecSF_.getScaleFactor(recolep, 11, -1);

  b_tri = b_tri_up = b_tri_dn = 0;
  b_tri = 1; //computeTrigSF(recolep1, recolep2);
  b_tri_up = 1; //computeTrigSF(recolep1, recolep2, 1);
  b_tri_dn = 1; //computeTrigSF(recolep1, recolep2, -1);

  b_met = MET_pt;

  auto muons = muonSelection();
  auto elecs = elecSelection();

  auto bjets = bjetSelection();
  b_nbjet = bjets.size();

  auto jets = jetSelection();
  b_njet = jets.size();

  if (muons.size() + elecs.size() != 1) return b_step;
  b_step = 1;
  if (h_cutFlow) h_cutFlow->Fill(3);

  TH1D * h_cutFlowLep = (elecs.size() == 1) ? h_cutFlowEl : h_cutFlowMu;
  if (h_cutFlowLep) h_cutFlowLep->Fill(1);
  
  if (muons.size() == 1) {
      recolep = muons[0];
      b_channel = CH_MU;
  } else if (elecs.size() == 1) {
      recolep = elecs[0];
      b_channel = CH_EL;
  }

  recolep.Momentum(b_lep);

  recoleps.push_back(b_lep);

  // Veto Leptons

  auto vetoMu = vetoMuonSelection();
  auto vetoEl = vetoElecSelection();

  if ((muons.size() == 0 && vetoMu.size() > 0) || (muons.size() == 1 && vetoMu.size() > 1))
    return b_step;
  if ((elecs.size() == 0 && vetoEl.size() > 0) || (elecs.size() == 1 && vetoEl.size() > 1))
    return b_step;
  
  b_step = 2;
  if (h_cutFlow) h_cutFlow->Fill(4);
  if (h_cutFlowLep) h_cutFlowLep->Fill(2);

  if (b_njet > 0) {
    b_step = 3;
    if (h_cutFlow) h_cutFlow->Fill(5);
    if (h_cutFlowLep) h_cutFlowLep->Fill(3);
  } else return b_step;
  
  if (b_nbjet > 0) {
    b_step = 4;
    if (h_cutFlow) h_cutFlow->Fill(6);
    if (h_cutFlowLep) h_cutFlowLep->Fill(4);
  } else return b_step;

  return b_step;
}
