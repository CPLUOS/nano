#include "nano/analysis/interface/topEventSelectionSL.h"

using std::vector;

topEventSelectionSL::topEventSelectionSL(TTree *tree, TTree *had, TTree *hadTruth, Bool_t isMC, Bool_t sle, Bool_t slm) :
  topObjectSelection(tree, had, hadTruth, isMC),
  h_nevents(0),
  h_genweights(0),
  h_cutFlow(0),
  h_cutFlowEl(0),
  h_cutFlowMu(0),
  m_isSL_e(sle),
  m_isSL_m(slm)
{
  SetCutValues();
}

topEventSelectionSL::~topEventSelectionSL()
{
}

int topEventSelectionSL::SetCutValues() {
  cut_ElectronPt = 30;
  cut_ElectronEta = 2.4;
  cut_ElectronIDType = Electron_cutBased;
  cut_ElectronIDCut = 3;
  cut_ElectronSCEtaLower = 1.4442;
  cut_ElectronSCEtaUpper = 1.566;
  cut_ElectronRelIso03All = 10000000000;
  
  cut_MuonIDType = Muon_tightId;
  cut_MuonPt = 26;
  cut_MuonEta = 2.1;
  cut_MuonRelIso04All = 0.15;
  
  cut_VetoElectronPt = 20;
  cut_VetoElectronEta = 2.4;
  cut_VetoElectronIDType = Electron_cutBased;
  cut_VetoElectronIDCut = 3;
  cut_VetoElectronSCEtaLower = 1.4442;
  cut_VetoElectronSCEtaUpper = 1.566;
  cut_VetoElectronRelIso03All = 10000000000;
  
  cut_VetoMuonIDType = NULL;
  cut_VetoMuonPt = 10;
  cut_VetoMuonEta = 2.4;
  cut_VetoMuonRelIso04All = 0.25;
  
  cut_GenJetPt = 30;
  cut_GenJetEta = 2.4;
  cut_GenJetConeSizeOverlap = 0.4;
  
  cut_JetID = 1;
  cut_JetPt = 30;
  cut_JetEta = 2.4;
  cut_JetConeSizeOverlap = 0.4;
  
  cut_BJetID = 1;
  cut_BJetPt = 30;
  cut_BJetEta = 2.4;
  cut_BJetConeSizeOverlap = 0.4;
  cut_BJetTypeBTag = Jet_btagCSVV2;
  cut_BJetBTagCut = 0.8484;
  
  return 0;
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
  b_trig_m = TrigForMu();
  b_trig_e = TrigForEl();
  if ( !( b_trig_m || b_trig_e ) ) return b_step;

  // TODO Check trigger requirements (TTbarXSecSynchronization page doesn't have yet)
  
  // if (b_channel == CH_MU) {
  //   if (!b_trig_m) return b_step;
  // }

  // if (b_channel == CH_EL) {
  //   if (!b_trig_e) return b_step;
  // }

  b_met = MET_pt;

  auto muons = muonSelection();
  auto elecs = elecSelection();

  if (muons.size() + elecs.size() != 1) return b_step;
  b_step = 1;
  if (h_cutFlow) h_cutFlow->Fill(3);

  TH1D * h_cutFlowLep = (elecs.size() == 1) ? h_cutFlowEl : h_cutFlowMu;
  if (h_cutFlowLep) h_cutFlowLep->Fill(1);
  
  if (muons.size() == 1) {
    recolep = muons[0];
    b_channel = CH_MU;
    
    b_mueffweight    = m_tableLeptonSFMu.getFactor(recolep.Pt(), recolep.Eta());
    b_mueffweight_up = m_tableLeptonSFMu.getFactor(recolep.Pt(), recolep.Eta(),  1);
    b_mueffweight_dn = m_tableLeptonSFMu.getFactor(recolep.Pt(), recolep.Eta(), -1);
  } else if (elecs.size() == 1) {
    recolep = elecs[0];
    b_channel = CH_EL;
    
    b_eleffweight    = m_tableLeptonSFEl.getFactor(recolep.Pt(), recolep.Eta());
    b_eleffweight_up = m_tableLeptonSFEl.getFactor(recolep.Pt(), recolep.Eta(),  1);
    b_eleffweight_dn = m_tableLeptonSFEl.getFactor(recolep.Pt(), recolep.Eta(), -1);
  }

  recolep.Momentum(b_lep);
  b_lep_pid = recolep.GetPdgCode();
  recoleps.push_back(b_lep);

  b_tri = b_tri_up = b_tri_dn = 0;
  computePtEtaTable &tableTrigSF = (std::abs(b_lep_pid) == 11 ? m_tableTrigSFEl : m_tableTrigSFMu);
  b_tri    = tableTrigSF.getFactor(recolep.Pt(), recolep.Eta());
  b_tri_up = tableTrigSF.getFactor(recolep.Pt(), recolep.Eta(),  1);
  b_tri_dn = tableTrigSF.getFactor(recolep.Pt(), recolep.Eta(), -1);

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

  auto bjets = bjetSelection();
  b_nbjet = bjets.size();

  auto jets = jetSelection();
  b_njet = jets.size();

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
