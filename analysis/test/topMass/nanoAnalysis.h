//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jan 18 21:36:33 2018 by ROOT version 6.10/09
// from TTree Events/Events
// found on file: /home/nanoAOD/run2_2016v3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/180117_180123/0000/nanoAOD_111.root
//////////////////////////////////////////////////////////

#ifndef nanoAnalysis_h
#define nanoAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include <TH1D.h>
#include <vector>

#include "../src/RoccoR.cc"
#include "../src/pileUpTool.h"
#include "../src/MuonScaleFactorEvaluator.h"
#include "../src/ElecScaleFactorEvaluator.h"
//#include "../plugins/jsoncpp.cpp"
#include "../src/lumiTool.h"
#include "../src/BTagWeightEvaluator.h"
// Header file for the classes stored in the TTree if any.

class nanoAnalysis {
  //Output Variables
private:
  TFile* m_output;      
  //Tree
  TTree* m_tree;
      
  //histogram
  TH1D* h_nevents;
  TH1D* h_genweights;
  TH1D* h_weights;
  TH1D* h_cutFlow;
      
  //Variables
  TLorentzVector b_lep1, b_lep2, b_dilep, b_jet1, b_jet2;
  TParticle recolep1, recolep2;
  int b_lep1_pid, b_lep2_pid;
  float b_jet1_CSVInclV2, b_jet2_CSVInclV2;
  std::vector<Float_t> b_csvweights;
  int b_nvertex, b_step, b_channel, b_njet, b_nbjet;
  bool b_step1, b_step2, b_step3, b_step4, b_step5, b_step6, b_step7;  
  float b_met, b_weight, b_genweight, b_puweight, b_btagweight;
  float b_mueffweight, b_mueffweight_up, b_mueffweight_dn,
        b_eleffweight, b_eleffweight_up, b_eleffweight_dn;

  // Tools
  RoccoR* m_rocCor;
  TH1D* hist_mc;
  MuonScaleFactorEvaluator muonSF_;
  ElecScaleFactorEvaluator elecSF_;
  BTagWeightEvaluator m_btagSF;
  Bool_t m_isMC;
  Bool_t m_isDL;
  Bool_t m_isSL_e;
  Bool_t m_isSL_m;

  pileUpTool *m_pileUp;
  lumiTool* m_lumi;
  //LumiMap
  std::map<UInt_t, std::vector<std::array<UInt_t, 2>>> lumiMap;

  void analysis();
  void mcAnalysis();
  
  //Making output branch
  void MakeBranch(TTree* t);
  void resetBranch();

  //For Selection
  Bool_t lumiCheck();
  std::vector<TParticle> muonSelection();
  Double_t roccoR(TLorentzVector m, int q, int nGen, int nTrackerLayers);
  enum TTLLChannel { CH_NOLL = 0, CH_MUEL, CH_ELEL, CH_MUMU }; 
  std::vector<TParticle> elecSelection();
  std::vector<TLorentzVector> recoleps;
  std::vector<TParticle> jetSelection();
  std::vector<TParticle> bjetSelection();

 public :
  //set output file
  void setOutput(std::string outputName);
  void LoadModules(pileUpTool* pileUp, lumiTool* lumi); 

  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // Fixed size dimensions of array or collections stored in the TTree if any.

  // Declaration of leaf types
  UInt_t          run;
  UInt_t          luminosityBlock;
  ULong64_t       event;
  Float_t         CaloMET_phi;
  Float_t         CaloMET_pt;
  Float_t         CaloMET_sumEt;
  UInt_t          ncmeson;
  Float_t         cmeson_dca[63];   //[ncmeson]
  Float_t         cmeson_angleXY[63];   //[ncmeson]
  Float_t         cmeson_angleXYZ[63];   //[ncmeson]
  Float_t         cmeson_trk_normalizedChi2[63];   //[ncmeson]
  Float_t         cmeson_trk_nHits[63];   //[ncmeson]
  Float_t         cmeson_trk_pt[63];   //[ncmeson]
  Float_t         cmeson_trk_ipsigXY[63];   //[ncmeson]
  Float_t         cmeson_trk_ipsigZ[63];   //[ncmeson]
  Float_t         cmeson_lxy[63];   //[ncmeson]
  Float_t         cmeson_lxySig[63];   //[ncmeson]
  Float_t         cmeson_l3D[63];   //[ncmeson]
  Float_t         cmeson_l3DSig[63];   //[ncmeson]
  Float_t         cmeson_jetDR[63];   //[ncmeson]
  Float_t         cmeson_legDR[63];   //[ncmeson]
  Float_t         cmeson_diffMass[63];   //[ncmeson]
  Int_t           cmeson_nJet[63];   //[ncmeson]
  Int_t           cmeson_mcMatch[63];   //[ncmeson]
  UInt_t          nElectron;
  Float_t         Electron_deltaEtaSC[6];   //[nElectron]
  Float_t         Electron_dr03EcalRecHitSumEt[6];   //[nElectron]
  Float_t         Electron_dr03HcalDepth1TowerSumEt[6];   //[nElectron]
  Float_t         Electron_dr03TkSumPt[6];   //[nElectron]
  Float_t         Electron_dxy[6];   //[nElectron]
  Float_t         Electron_dxyErr[6];   //[nElectron]
  Float_t         Electron_dz[6];   //[nElectron]
  Float_t         Electron_dzErr[6];   //[nElectron]
  Float_t         Electron_eCorr[6];   //[nElectron]
  Float_t         Electron_eInvMinusPInv[6];   //[nElectron]
  Float_t         Electron_energyErr[6];   //[nElectron]
  Float_t         Electron_eta[6];   //[nElectron]
  Float_t         Electron_hoe[6];   //[nElectron]
  Float_t         Electron_ip3d[6];   //[nElectron]
  Float_t         Electron_mass[6];   //[nElectron]
  Float_t         Electron_miniPFRelIso_all[6];   //[nElectron]
  Float_t         Electron_miniPFRelIso_chg[6];   //[nElectron]
  Float_t         Electron_mvaSpring16GP[6];   //[nElectron]
  Float_t         Electron_mvaSpring16HZZ[6];   //[nElectron]
  Float_t         Electron_pfRelIso03_all[6];   //[nElectron]
  Float_t         Electron_pfRelIso03_chg[6];   //[nElectron]
  Float_t         Electron_phi[6];   //[nElectron]
  Float_t         Electron_pt[6];   //[nElectron]
  Float_t         Electron_r9[6];   //[nElectron]
  Float_t         Electron_sieie[6];   //[nElectron]
  Float_t         Electron_sip3d[6];   //[nElectron]
  Float_t         Electron_mvaTTH[6];   //[nElectron]
  Int_t           Electron_charge[6];   //[nElectron]
  Int_t           Electron_cutBased[6];   //[nElectron]
  Int_t           Electron_cutBased_HLTPreSel[6];   //[nElectron]
  Int_t           Electron_jetIdx[6];   //[nElectron]
  Int_t           Electron_pdgId[6];   //[nElectron]
  Int_t           Electron_photonIdx[6];   //[nElectron]
  Int_t           Electron_tightCharge[6];   //[nElectron]
  Int_t           Electron_vidNestedWPBitmap[6];   //[nElectron]
  Bool_t          Electron_convVeto[6];   //[nElectron]
  Bool_t          Electron_cutBased_HEEP[6];   //[nElectron]
  Bool_t          Electron_isPFcand[6];   //[nElectron]
  UChar_t         Electron_lostHits[6];   //[nElectron]
  Bool_t          Electron_mvaSpring16GP_WP80[6];   //[nElectron]
  Bool_t          Electron_mvaSpring16GP_WP90[6];   //[nElectron]
  Bool_t          Electron_mvaSpring16HZZ_WPL[6];   //[nElectron]
  UInt_t          nFatJet;
  Float_t         FatJet_area[4];   //[nFatJet]
  Float_t         FatJet_btagCMVA[4];   //[nFatJet]
  Float_t         FatJet_btagCSVV2[4];   //[nFatJet]
  Float_t         FatJet_btagDeepB[4];   //[nFatJet]
  Float_t         FatJet_btagHbb[4];   //[nFatJet]
  Float_t         FatJet_eta[4];   //[nFatJet]
  Float_t         FatJet_mass[4];   //[nFatJet]
  Float_t         FatJet_msoftdrop[4];   //[nFatJet]
  Float_t         FatJet_n2b1[4];   //[nFatJet]
  Float_t         FatJet_n3b1[4];   //[nFatJet]
  Float_t         FatJet_phi[4];   //[nFatJet]
  Float_t         FatJet_pt[4];   //[nFatJet]
  Float_t         FatJet_tau1[4];   //[nFatJet]
  Float_t         FatJet_tau2[4];   //[nFatJet]
  Float_t         FatJet_tau3[4];   //[nFatJet]
  Float_t         FatJet_tau4[4];   //[nFatJet]
  Int_t           FatJet_subJetIdx1[4];   //[nFatJet]
  Int_t           FatJet_subJetIdx2[4];   //[nFatJet]
  UInt_t          nGenJetAK8;
  Float_t         GenJetAK8_eta[4];   //[nGenJetAK8]
  Float_t         GenJetAK8_mass[4];   //[nGenJetAK8]
  Float_t         GenJetAK8_phi[4];   //[nGenJetAK8]
  Float_t         GenJetAK8_pt[4];   //[nGenJetAK8]
  UInt_t          nGenJet;
  Float_t         GenJet_eta[20];   //[nGenJet]
  Float_t         GenJet_mass[20];   //[nGenJet]
  Float_t         GenJet_phi[20];   //[nGenJet]
  Float_t         GenJet_pt[20];   //[nGenJet]
  UInt_t          nGenPart;
  Float_t         GenPart_eta[102];   //[nGenPart]
  Float_t         GenPart_mass[102];   //[nGenPart]
  Float_t         GenPart_phi[102];   //[nGenPart]
  Float_t         GenPart_pt[102];   //[nGenPart]
  Int_t           GenPart_genPartIdxMother[102];   //[nGenPart]
  Int_t           GenPart_pdgId[102];   //[nGenPart]
  Int_t           GenPart_status[102];   //[nGenPart]
  Int_t           GenPart_statusFlags[102];   //[nGenPart]
  Float_t         Generator_scalePDF;
  Float_t         Generator_x1;
  Float_t         Generator_x2;
  Float_t         Generator_xpdf1;
  Float_t         Generator_xpdf2;
  Int_t           Generator_id1;
  Int_t           Generator_id2;
  UInt_t          nGenVisTau;
  Float_t         GenVisTau_eta[3];   //[nGenVisTau]
  Float_t         GenVisTau_mass[3];   //[nGenVisTau]
  Float_t         GenVisTau_phi[3];   //[nGenVisTau]
  Float_t         GenVisTau_pt[3];   //[nGenVisTau]
  Int_t           GenVisTau_charge[3];   //[nGenVisTau]
  Int_t           GenVisTau_genPartIdxMother[3];   //[nGenVisTau]
  Int_t           GenVisTau_status[3];   //[nGenVisTau]
  Float_t         genWeight;
  Float_t         LHEWeight_originalXWGTUP;
  UInt_t          nLHEPdfWeight;
  Float_t         LHEPdfWeight[1];   //[nLHEPdfWeight]
  UInt_t          nLHEScaleWeight;
  Float_t         LHEScaleWeight[9];   //[nLHEScaleWeight]
  UInt_t          nJet;
  Float_t         Jet_area[35];   //[nJet]
  Float_t         Jet_btagCMVA[35];   //[nJet]
  Float_t         Jet_btagCSVV2[35];   //[nJet]
  Float_t         Jet_btagDeepB[35];   //[nJet]
  Float_t         Jet_btagDeepC[35];   //[nJet]
  Float_t         Jet_chEmEF[35];   //[nJet]
  Float_t         Jet_chHEF[35];   //[nJet]
  Float_t         Jet_eta[35];   //[nJet]
  Float_t         Jet_mass[35];   //[nJet]
  Float_t         Jet_neEmEF[35];   //[nJet]
  Float_t         Jet_neHEF[35];   //[nJet]
  Float_t         Jet_phi[35];   //[nJet]
  Float_t         Jet_pt[35];   //[nJet]
  Float_t         Jet_qgl[35];   //[nJet]
  Float_t         Jet_rawFactor[35];   //[nJet]
  Float_t         Jet_bReg[35];   //[nJet]
  Int_t           Jet_electronIdx1[35];   //[nJet]
  Int_t           Jet_electronIdx2[35];   //[nJet]
  Int_t           Jet_jetId[35];   //[nJet]
  Int_t           Jet_muonIdx1[35];   //[nJet]
  Int_t           Jet_muonIdx2[35];   //[nJet]
  Int_t           Jet_nConstituents[35];   //[nJet]
  Int_t           Jet_nElectrons[35];   //[nJet]
  Int_t           Jet_nMuons[35];   //[nJet]
  Int_t           Jet_puId[35];   //[nJet]
  Float_t         LHE_HT;
  Float_t         LHE_HTIncoming;
  Float_t         LHE_Vpt;
  UChar_t         LHE_Njets;
  UChar_t         LHE_Nb;
  UChar_t         LHE_Nc;
  UChar_t         LHE_Nuds;
  UChar_t         LHE_Nglu;
  UChar_t         LHE_NpNLO;
  UChar_t         LHE_NpLO;
  Float_t         GenMET_phi;
  Float_t         GenMET_pt;
  Float_t         MET_MetUnclustEnUpDeltaX;
  Float_t         MET_MetUnclustEnUpDeltaY;
  Float_t         MET_covXX;
  Float_t         MET_covXY;
  Float_t         MET_covYY;
  Float_t         MET_phi;
  Float_t         MET_pt;
  Float_t         MET_significance;
  Float_t         MET_sumEt;
  UInt_t          nMuon;
  Float_t         Muon_dxy[10];   //[nMuon]
  Float_t         Muon_dxyErr[10];   //[nMuon]
  Float_t         Muon_dz[10];   //[nMuon]
  Float_t         Muon_dzErr[10];   //[nMuon]
  Float_t         Muon_eta[10];   //[nMuon]
  Float_t         Muon_ip3d[10];   //[nMuon]
  Float_t         Muon_mass[10];   //[nMuon]
  Float_t         Muon_miniPFRelIso_all[10];   //[nMuon]
  Float_t         Muon_miniPFRelIso_chg[10];   //[nMuon]
  Float_t         Muon_pfRelIso03_all[10];   //[nMuon]
  Float_t         Muon_pfRelIso03_chg[10];   //[nMuon]
  Float_t         Muon_pfRelIso04_all[10];   //[nMuon]
  Float_t         Muon_phi[10];   //[nMuon]
  Float_t         Muon_pt[10];   //[nMuon]
  Float_t         Muon_ptErr[10];   //[nMuon]
  Float_t         Muon_segmentComp[10];   //[nMuon]
  Float_t         Muon_sip3d[10];   //[nMuon]
  Float_t         Muon_mvaTTH[10];   //[nMuon]
  Int_t           Muon_charge[10];   //[nMuon]
  Int_t           Muon_jetIdx[10];   //[nMuon]
  Int_t           Muon_nStations[10];   //[nMuon]
  Int_t           Muon_nTrackerLayers[10];   //[nMuon]
  Int_t           Muon_pdgId[10];   //[nMuon]
  Int_t           Muon_tightCharge[10];   //[nMuon]
  Bool_t          Muon_globalMu[10];   //[nMuon]
  UChar_t         Muon_highPtId[10];   //[nMuon]
  Bool_t          Muon_isPFcand[10];   //[nMuon]
  Bool_t          Muon_mediumId[10];   //[nMuon]
  Bool_t          Muon_softId[10];   //[nMuon]
  Bool_t          Muon_tightId[10];   //[nMuon]
  Bool_t          Muon_trackerMu[10];   //[nMuon]
  UInt_t          nPhoton;
  Float_t         Photon_eCorr[7];   //[nPhoton]
  Float_t         Photon_energyErr[7];   //[nPhoton]
  Float_t         Photon_eta[7];   //[nPhoton]
  Float_t         Photon_hoe[7];   //[nPhoton]
  Float_t         Photon_mass[7];   //[nPhoton]
  Float_t         Photon_mvaID[7];   //[nPhoton]
  Float_t         Photon_pfRelIso03_all[7];   //[nPhoton]
  Float_t         Photon_pfRelIso03_chg[7];   //[nPhoton]
  Float_t         Photon_phi[7];   //[nPhoton]
  Float_t         Photon_pt[7];   //[nPhoton]
  Float_t         Photon_r9[7];   //[nPhoton]
  Float_t         Photon_sieie[7];   //[nPhoton]
  Int_t           Photon_charge[7];   //[nPhoton]
  Int_t           Photon_cutBased[7];   //[nPhoton]
  Int_t           Photon_electronIdx[7];   //[nPhoton]
  Int_t           Photon_jetIdx[7];   //[nPhoton]
  Int_t           Photon_pdgId[7];   //[nPhoton]
  Int_t           Photon_vidNestedWPBitmap[7];   //[nPhoton]
  Bool_t          Photon_electronVeto[7];   //[nPhoton]
  Bool_t          Photon_mvaID_WP80[7];   //[nPhoton]
  Bool_t          Photon_mvaID_WP90[7];   //[nPhoton]
  Bool_t          Photon_pixelSeed[7];   //[nPhoton]
  Int_t           Pileup_nPU;
  Float_t         Pileup_nTrueInt;
  Float_t         PuppiMET_phi;
  Float_t         PuppiMET_pt;
  Float_t         PuppiMET_sumEt;
  Float_t         RawMET_phi;
  Float_t         RawMET_pt;
  Float_t         RawMET_sumEt;
  Float_t         fixedGridRhoFastjetAll;
  Float_t         fixedGridRhoFastjetCentralCalo;
  Float_t         fixedGridRhoFastjetCentralNeutral;
  UInt_t          nGenDressedLepton;
  Float_t         GenDressedLepton_eta[3];   //[nGenDressedLepton]
  Float_t         GenDressedLepton_mass[3];   //[nGenDressedLepton]
  Float_t         GenDressedLepton_phi[3];   //[nGenDressedLepton]
  Float_t         GenDressedLepton_pt[3];   //[nGenDressedLepton]
  Int_t           GenDressedLepton_pdgId[3];   //[nGenDressedLepton]
  UInt_t          nSoftActivityJet;
  Float_t         SoftActivityJet_eta[6];   //[nSoftActivityJet]
  Float_t         SoftActivityJet_phi[6];   //[nSoftActivityJet]
  Float_t         SoftActivityJet_pt[6];   //[nSoftActivityJet]
  Float_t         SoftActivityJetHT;
  Float_t         SoftActivityJetHT10;
  Float_t         SoftActivityJetHT2;
  Float_t         SoftActivityJetHT5;
  Int_t           SoftActivityJetNjets10;
  Int_t           SoftActivityJetNjets2;
  Int_t           SoftActivityJetNjets5;
  UInt_t          nSubJet;
  Float_t         SubJet_btagCMVA[6];   //[nSubJet]
  Float_t         SubJet_btagCSVV2[6];   //[nSubJet]
  Float_t         SubJet_btagDeepB[6];   //[nSubJet]
  Float_t         SubJet_eta[6];   //[nSubJet]
  Float_t         SubJet_mass[6];   //[nSubJet]
  Float_t         SubJet_n2b1[6];   //[nSubJet]
  Float_t         SubJet_n3b1[6];   //[nSubJet]
  Float_t         SubJet_phi[6];   //[nSubJet]
  Float_t         SubJet_pt[6];   //[nSubJet]
  Float_t         SubJet_tau1[6];   //[nSubJet]
  Float_t         SubJet_tau2[6];   //[nSubJet]
  Float_t         SubJet_tau3[6];   //[nSubJet]
  Float_t         SubJet_tau4[6];   //[nSubJet]
  UInt_t          nTau;
  Float_t         Tau_chargedIso[5];   //[nTau]
  Float_t         Tau_dxy[5];   //[nTau]
  Float_t         Tau_dz[5];   //[nTau]
  Float_t         Tau_eta[5];   //[nTau]
  Float_t         Tau_footprintCorr[5];   //[nTau]
  Float_t         Tau_leadTkDeltaEta[5];   //[nTau]
  Float_t         Tau_leadTkDeltaPhi[5];   //[nTau]
  Float_t         Tau_leadTkPtOverTauPt[5];   //[nTau]
  Float_t         Tau_mass[5];   //[nTau]
  Float_t         Tau_neutralIso[5];   //[nTau]
  Float_t         Tau_phi[5];   //[nTau]
  Float_t         Tau_photonsOutsideSignalCone[5];   //[nTau]
  Float_t         Tau_pt[5];   //[nTau]
  Float_t         Tau_puCorr[5];   //[nTau]
  Float_t         Tau_rawAntiEle[5];   //[nTau]
  Float_t         Tau_rawIso[5];   //[nTau]
  Float_t         Tau_rawMVAnewDM[5];   //[nTau]
  Float_t         Tau_rawMVAoldDM[5];   //[nTau]
  Float_t         Tau_rawMVAoldDMdR03[5];   //[nTau]
  Int_t           Tau_charge[5];   //[nTau]
  Int_t           Tau_decayMode[5];   //[nTau]
  Int_t           Tau_jetIdx[5];   //[nTau]
  Int_t           Tau_rawAntiEleCat[5];   //[nTau]
  UChar_t         Tau_idAntiEle[5];   //[nTau]
  UChar_t         Tau_idAntiMu[5];   //[nTau]
  Bool_t          Tau_idDecayMode[5];   //[nTau]
  Bool_t          Tau_idDecayModeNewDMs[5];   //[nTau]
  UChar_t         Tau_idMVAnewDM[5];   //[nTau]
  UChar_t         Tau_idMVAoldDM[5];   //[nTau]
  UChar_t         Tau_idMVAoldDMdR03[5];   //[nTau]
  Float_t         TkMET_phi;
  Float_t         TkMET_pt;
  Float_t         TkMET_sumEt;
  UInt_t          nTrigObj;
  Float_t         TrigObj_pt[28];   //[nTrigObj]
  Float_t         TrigObj_eta[28];   //[nTrigObj]
  Float_t         TrigObj_phi[28];   //[nTrigObj]
  Float_t         TrigObj_l1pt[28];   //[nTrigObj]
  Float_t         TrigObj_l1pt_2[28];   //[nTrigObj]
  Float_t         TrigObj_l2pt[28];   //[nTrigObj]
  Int_t           TrigObj_id[28];   //[nTrigObj]
  Int_t           TrigObj_l1iso[28];   //[nTrigObj]
  Int_t           TrigObj_l1charge[28];   //[nTrigObj]
  Int_t           TrigObj_filterBits[28];   //[nTrigObj]
  Int_t           genTtbarId;
  UInt_t          nOtherPV;
  Float_t         OtherPV_z[3];   //[nOtherPV]
  Float_t         PV_ndof;
  Float_t         PV_x;
  Float_t         PV_y;
  Float_t         PV_z;
  Float_t         PV_chi2;
  Float_t         PV_score;
  Int_t           PV_npvs;
  UInt_t          nSV;
  Float_t         SV_dlen[9];   //[nSV]
  Float_t         SV_dlenSig[9];   //[nSV]
  Float_t         SV_pAngle[9];   //[nSV]
  Float_t         cmeson_chi2[63];   //[ncmeson]
  Float_t         cmeson_eta[63];   //[ncmeson]
  Float_t         cmeson_mass[63];   //[ncmeson]
  Float_t         cmeson_phi[63];   //[ncmeson]
  Float_t         cmeson_pt[63];   //[ncmeson]
  Float_t         cmeson_x[63];   //[ncmeson]
  Float_t         cmeson_y[63];   //[ncmeson]
  Float_t         cmeson_z[63];   //[ncmeson]
  Int_t           cmeson_ndof[63];   //[ncmeson]
  Int_t           cmeson_pdgId[63];   //[ncmeson]
  Int_t           Electron_genPartIdx[6];   //[nElectron]
  UChar_t         Electron_genPartFlav[6];   //[nElectron]
  Int_t           GenJetAK8_partonFlavour[4];   //[nGenJetAK8]
  UChar_t         GenJetAK8_hadronFlavour[4];   //[nGenJetAK8]
  Int_t           GenJet_partonFlavour[20];   //[nGenJet]
  UChar_t         GenJet_hadronFlavour[20];   //[nGenJet]
  Int_t           Jet_genJetIdx[35];   //[nJet]
  Int_t           Jet_hadronFlavour[35];   //[nJet]
  Int_t           Jet_partonFlavour[35];   //[nJet]
  Int_t           Muon_genPartIdx[10];   //[nMuon]
  UChar_t         Muon_genPartFlav[10];   //[nMuon]
  Int_t           Photon_genPartIdx[7];   //[nPhoton]
  UChar_t         Photon_genPartFlav[7];   //[nPhoton]
  Float_t         MET_fiducialGenPhi;
  Float_t         MET_fiducialGenPt;
  UChar_t         Electron_cleanmask[6];   //[nElectron]
  UChar_t         Jet_cleanmask[35];   //[nJet]
  UChar_t         Muon_cleanmask[10];   //[nMuon]
  UChar_t         Photon_cleanmask[7];   //[nPhoton]
  UChar_t         Tau_cleanmask[5];   //[nTau]
  Float_t         SV_chi2[9];   //[nSV]
  Float_t         SV_eta[9];   //[nSV]
  Float_t         SV_mass[9];   //[nSV]
  Float_t         SV_ndof[9];   //[nSV]
  Float_t         SV_phi[9];   //[nSV]
  Float_t         SV_pt[9];   //[nSV]
  Float_t         SV_x[9];   //[nSV]
  Float_t         SV_y[9];   //[nSV]
  Float_t         SV_z[9];   //[nSV]
  Int_t           Tau_genPartIdx[5];   //[nTau]
  UChar_t         Tau_genPartFlav[5];   //[nTau]
  Bool_t          L1simulation_step;
  Bool_t          HLTriggerFirstPath;
  Bool_t          HLT_AK8PFJet360_TrimMass30;
  Bool_t          HLT_AK8PFJet400_TrimMass30;
  Bool_t          HLT_AK8PFHT750_TrimMass50;
  Bool_t          HLT_AK8PFHT800_TrimMass50;
  Bool_t          HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20;
  Bool_t          HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087;
  Bool_t          HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087;
  Bool_t          HLT_AK8DiPFJet300_200_TrimMass30;
  Bool_t          HLT_AK8PFHT700_TrimR0p1PT0p03Mass50;
  Bool_t          HLT_AK8PFHT650_TrimR0p1PT0p03Mass50;
  Bool_t          HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20;
  Bool_t          HLT_AK8DiPFJet280_200_TrimMass30;
  Bool_t          HLT_AK8DiPFJet250_200_TrimMass30;
  Bool_t          HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20;
  Bool_t          HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20;
  Bool_t          HLT_CaloJet260;
  Bool_t          HLT_CaloJet500_NoJetID;
  Bool_t          HLT_Dimuon13_PsiPrime;
  Bool_t          HLT_Dimuon13_Upsilon;
  Bool_t          HLT_Dimuon20_Jpsi;
  Bool_t          HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf;
  Bool_t          HLT_DoubleEle25_CaloIdL_GsfTrkIdVL;
  Bool_t          HLT_DoubleEle33_CaloIdL;
  Bool_t          HLT_DoubleEle33_CaloIdL_MW;
  Bool_t          HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW;
  Bool_t          HLT_DoubleEle33_CaloIdL_GsfTrkIdVL;
  Bool_t          HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg;
  Bool_t          HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg;
  Bool_t          HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg;
  Bool_t          HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg;
  Bool_t          HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1;
  Bool_t          HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1;
  Bool_t          HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg;
  Bool_t          HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg;
  Bool_t          HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1;
  Bool_t          HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL;
  Bool_t          HLT_DoubleMu33NoFiltersNoVtx;
  Bool_t          HLT_DoubleMu38NoFiltersNoVtx;
  Bool_t          HLT_DoubleMu23NoFiltersNoVtxDisplaced;
  Bool_t          HLT_DoubleMu28NoFiltersNoVtxDisplaced;
  Bool_t          HLT_DoubleMu0;
  Bool_t          HLT_DoubleMu4_3_Bs;
  Bool_t          HLT_DoubleMu4_3_Jpsi_Displaced;
  Bool_t          HLT_DoubleMu4_JpsiTrk_Displaced;
  Bool_t          HLT_DoubleMu4_LowMassNonResonantTrk_Displaced;
  Bool_t          HLT_DoubleMu3_Trk_Tau3mu;
  Bool_t          HLT_DoubleMu4_PsiPrimeTrk_Displaced;
  Bool_t          HLT_Mu7p5_L2Mu2_Jpsi;
  Bool_t          HLT_Mu7p5_L2Mu2_Upsilon;
  Bool_t          HLT_Mu7p5_Track2_Jpsi;
  Bool_t          HLT_Mu7p5_Track3p5_Jpsi;
  Bool_t          HLT_Mu7p5_Track7_Jpsi;
  Bool_t          HLT_Mu7p5_Track2_Upsilon;
  Bool_t          HLT_Mu7p5_Track3p5_Upsilon;
  Bool_t          HLT_Mu7p5_Track7_Upsilon;
  Bool_t          HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing;
  Bool_t          HLT_Dimuon0er16_Jpsi_NoVertexing;
  Bool_t          HLT_Dimuon6_Jpsi_NoVertexing;
  Bool_t          HLT_Photon150;
  Bool_t          HLT_Photon90_CaloIdL_HT300;
  Bool_t          HLT_HT250_CaloMET70;
  Bool_t          HLT_DoublePhoton60;
  Bool_t          HLT_DoublePhoton85;
  Bool_t          HLT_Ele17_Ele8_Gsf;
  Bool_t          HLT_Ele20_eta2p1_WPLoose_Gsf_LooseIsoPFTau28;
  Bool_t          HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau29;
  Bool_t          HLT_Ele22_eta2p1_WPLoose_Gsf;
  Bool_t          HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1;
  Bool_t          HLT_Ele23_WPLoose_Gsf;
  Bool_t          HLT_Ele23_WPLoose_Gsf_WHbbBoost;
  Bool_t          HLT_Ele24_eta2p1_WPLoose_Gsf;
  Bool_t          HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20;
  Bool_t          HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1;
  Bool_t          HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau30;
  Bool_t          HLT_Ele25_WPTight_Gsf;
  Bool_t          HLT_Ele25_eta2p1_WPLoose_Gsf;
  Bool_t          HLT_Ele25_eta2p1_WPTight_Gsf;
  Bool_t          HLT_Ele27_WPLoose_Gsf;
  Bool_t          HLT_Ele27_WPLoose_Gsf_WHbbBoost;
  Bool_t          HLT_Ele27_WPTight_Gsf;
  Bool_t          HLT_Ele27_WPTight_Gsf_L1JetTauSeeded;
  Bool_t          HLT_Ele27_eta2p1_WPLoose_Gsf;
  Bool_t          HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1;
  Bool_t          HLT_Ele27_eta2p1_WPTight_Gsf;
  Bool_t          HLT_Ele30_WPTight_Gsf;
  Bool_t          HLT_Ele30_eta2p1_WPLoose_Gsf;
  Bool_t          HLT_Ele30_eta2p1_WPTight_Gsf;
  Bool_t          HLT_Ele32_WPTight_Gsf;
  Bool_t          HLT_Ele32_eta2p1_WPLoose_Gsf;
  Bool_t          HLT_Ele32_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1;
  Bool_t          HLT_Ele32_eta2p1_WPTight_Gsf;
  Bool_t          HLT_Ele35_WPLoose_Gsf;
  Bool_t          HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50;
  Bool_t          HLT_Ele36_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1;
  Bool_t          HLT_Ele45_WPLoose_Gsf;
  Bool_t          HLT_Ele45_WPLoose_Gsf_L1JetTauSeeded;
  Bool_t          HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50;
  Bool_t          HLT_Ele105_CaloIdVT_GsfTrkIdT;
  Bool_t          HLT_Ele30WP60_SC4_Mass55;
  Bool_t          HLT_Ele30WP60_Ele8_Mass55;
  Bool_t          HLT_HT200;
  Bool_t          HLT_HT275;
  Bool_t          HLT_HT325;
  Bool_t          HLT_HT425;
  Bool_t          HLT_HT575;
  Bool_t          HLT_HT410to430;
  Bool_t          HLT_HT430to450;
  Bool_t          HLT_HT450to470;
  Bool_t          HLT_HT470to500;
  Bool_t          HLT_HT500to550;
  Bool_t          HLT_HT550to650;
  Bool_t          HLT_HT650;
  Bool_t          HLT_Mu16_eta2p1_MET30;
  Bool_t          HLT_IsoMu16_eta2p1_MET30;
  Bool_t          HLT_IsoMu16_eta2p1_MET30_LooseIsoPFTau50_Trk30_eta2p1;
  Bool_t          HLT_IsoMu17_eta2p1;
  Bool_t          HLT_IsoMu17_eta2p1_LooseIsoPFTau20;
  Bool_t          HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1;
  Bool_t          HLT_DoubleIsoMu17_eta2p1;
  Bool_t          HLT_DoubleIsoMu17_eta2p1_noDzCut;
  Bool_t          HLT_IsoMu18;
  Bool_t          HLT_IsoMu19_eta2p1_LooseIsoPFTau20;
  Bool_t          HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1;
  Bool_t          HLT_IsoMu19_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg;
  Bool_t          HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20;
  Bool_t          HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg;
  Bool_t          HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg;
  Bool_t          HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg;
  Bool_t          HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg;
  Bool_t          HLT_IsoMu20;
  Bool_t          HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1;
  Bool_t          HLT_IsoMu21_eta2p1_LooseIsoPFTau50_Trk30_eta2p1_SingleL1;
  Bool_t          HLT_IsoMu21_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg;
  Bool_t          HLT_IsoMu22;
  Bool_t          HLT_IsoMu22_eta2p1;
  Bool_t          HLT_IsoMu24;
  Bool_t          HLT_IsoMu27;
  Bool_t          HLT_IsoTkMu18;
  Bool_t          HLT_IsoTkMu20;
  Bool_t          HLT_IsoTkMu22;
  Bool_t          HLT_IsoTkMu22_eta2p1;
  Bool_t          HLT_IsoTkMu24;
  Bool_t          HLT_IsoTkMu27;
  Bool_t          HLT_JetE30_NoBPTX3BX;
  Bool_t          HLT_JetE30_NoBPTX;
  Bool_t          HLT_JetE50_NoBPTX3BX;
  Bool_t          HLT_JetE70_NoBPTX3BX;
  Bool_t          HLT_L1SingleMu18;
  Bool_t          HLT_L2Mu10;
  Bool_t          HLT_L1SingleMuOpen;
  Bool_t          HLT_L1SingleMuOpen_DT;
  Bool_t          HLT_L2DoubleMu23_NoVertex;
  Bool_t          HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10;
  Bool_t          HLT_L2DoubleMu38_NoVertex_2Cha_Angle2p5_Mass10;
  Bool_t          HLT_L2Mu10_NoVertex_NoBPTX3BX;
  Bool_t          HLT_L2Mu10_NoVertex_NoBPTX;
  Bool_t          HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX;
  Bool_t          HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX;
  Bool_t          HLT_LooseIsoPFTau50_Trk30_eta2p1;
  Bool_t          HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80;
  Bool_t          HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90;
  Bool_t          HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110;
  Bool_t          HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120;
  Bool_t          HLT_PFTau120_eta2p1;
  Bool_t          HLT_PFTau140_eta2p1;
  Bool_t          HLT_VLooseIsoPFTau120_Trk50_eta2p1;
  Bool_t          HLT_VLooseIsoPFTau140_Trk50_eta2p1;
  Bool_t          HLT_Mu17_Mu8;
  Bool_t          HLT_Mu17_Mu8_DZ;
  Bool_t          HLT_Mu17_Mu8_SameSign;
  Bool_t          HLT_Mu17_Mu8_SameSign_DZ;
  Bool_t          HLT_Mu20_Mu10;
  Bool_t          HLT_Mu20_Mu10_DZ;
  Bool_t          HLT_Mu20_Mu10_SameSign;
  Bool_t          HLT_Mu20_Mu10_SameSign_DZ;
  Bool_t          HLT_Mu17_TkMu8_DZ;
  Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;
  Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;
  Bool_t          HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL;
  Bool_t          HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;
  Bool_t          HLT_Mu25_TkMu0_dEta18_Onia;
  Bool_t          HLT_Mu27_TkMu8;
  Bool_t          HLT_Mu30_TkMu11;
  Bool_t          HLT_Mu30_eta2p1_PFJet150_PFJet50;
  Bool_t          HLT_Mu40_TkMu11;
  Bool_t          HLT_Mu40_eta2p1_PFJet200_PFJet50;
  Bool_t          HLT_Mu20;
  Bool_t          HLT_TkMu17;
  Bool_t          HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL;
  Bool_t          HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;
  Bool_t          HLT_TkMu20;
  Bool_t          HLT_Mu24_eta2p1;
  Bool_t          HLT_TkMu24_eta2p1;
  Bool_t          HLT_Mu27;
  Bool_t          HLT_TkMu27;
  Bool_t          HLT_Mu45_eta2p1;
  Bool_t          HLT_Mu50;
  Bool_t          HLT_TkMu50;
  Bool_t          HLT_Mu38NoFiltersNoVtx_Photon38_CaloIdL;
  Bool_t          HLT_Mu42NoFiltersNoVtx_Photon42_CaloIdL;
  Bool_t          HLT_Mu28NoFiltersNoVtxDisplaced_Photon28_CaloIdL;
  Bool_t          HLT_Mu33NoFiltersNoVtxDisplaced_Photon33_CaloIdL;
  Bool_t          HLT_Mu23NoFiltersNoVtx_Photon23_CaloIdL;
  Bool_t          HLT_DoubleMu18NoFiltersNoVtx;
  Bool_t          HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Tight;
  Bool_t          HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Loose;
  Bool_t          HLT_Mu28NoFiltersNoVtx_DisplacedJet40_Loose;
  Bool_t          HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Tight;
  Bool_t          HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Loose;
  Bool_t          HLT_Mu38NoFiltersNoVtx_DisplacedJet60_Loose;
  Bool_t          HLT_Mu28NoFiltersNoVtx_CentralCaloJet40;
  Bool_t          HLT_PFHT300_PFMET100;
  Bool_t          HLT_PFHT300_PFMET110;
  Bool_t          HLT_PFHT550_4JetPt50;
  Bool_t          HLT_PFHT650_4JetPt50;
  Bool_t          HLT_PFHT750_4JetPt50;
  Bool_t          HLT_PFHT750_4JetPt70;
  Bool_t          HLT_PFHT750_4JetPt80;
  Bool_t          HLT_PFHT800_4JetPt50;
  Bool_t          HLT_PFHT850_4JetPt50;
  Bool_t          HLT_PFJet15_NoCaloMatched;
  Bool_t          HLT_PFJet25_NoCaloMatched;
  Bool_t          HLT_DiPFJet15_NoCaloMatched;
  Bool_t          HLT_DiPFJet25_NoCaloMatched;
  Bool_t          HLT_DiPFJet15_FBEta3_NoCaloMatched;
  Bool_t          HLT_DiPFJet25_FBEta3_NoCaloMatched;
  Bool_t          HLT_DiPFJetAve15_HFJEC;
  Bool_t          HLT_DiPFJetAve25_HFJEC;
  Bool_t          HLT_DiPFJetAve35_HFJEC;
  Bool_t          HLT_AK8PFJet40;
  Bool_t          HLT_AK8PFJet60;
  Bool_t          HLT_AK8PFJet80;
  Bool_t          HLT_AK8PFJet140;
  Bool_t          HLT_AK8PFJet200;
  Bool_t          HLT_AK8PFJet260;
  Bool_t          HLT_AK8PFJet320;
  Bool_t          HLT_AK8PFJet400;
  Bool_t          HLT_AK8PFJet450;
  Bool_t          HLT_AK8PFJet500;
  Bool_t          HLT_PFJet40;
  Bool_t          HLT_PFJet60;
  Bool_t          HLT_PFJet80;
  Bool_t          HLT_PFJet140;
  Bool_t          HLT_PFJet200;
  Bool_t          HLT_PFJet260;
  Bool_t          HLT_PFJet320;
  Bool_t          HLT_PFJet400;
  Bool_t          HLT_PFJet450;
  Bool_t          HLT_PFJet500;
  Bool_t          HLT_DiPFJetAve40;
  Bool_t          HLT_DiPFJetAve60;
  Bool_t          HLT_DiPFJetAve80;
  Bool_t          HLT_DiPFJetAve140;
  Bool_t          HLT_DiPFJetAve200;
  Bool_t          HLT_DiPFJetAve260;
  Bool_t          HLT_DiPFJetAve320;
  Bool_t          HLT_DiPFJetAve400;
  Bool_t          HLT_DiPFJetAve500;
  Bool_t          HLT_DiPFJetAve60_HFJEC;
  Bool_t          HLT_DiPFJetAve80_HFJEC;
  Bool_t          HLT_DiPFJetAve100_HFJEC;
  Bool_t          HLT_DiPFJetAve160_HFJEC;
  Bool_t          HLT_DiPFJetAve220_HFJEC;
  Bool_t          HLT_DiPFJetAve300_HFJEC;
  Bool_t          HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu140;
  Bool_t          HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu80;
  Bool_t          HLT_DiCentralPFJet170;
  Bool_t          HLT_SingleCentralPFJet170_CFMax0p1;
  Bool_t          HLT_DiCentralPFJet170_CFMax0p1;
  Bool_t          HLT_DiCentralPFJet220_CFMax0p3;
  Bool_t          HLT_DiCentralPFJet330_CFMax0p5;
  Bool_t          HLT_DiCentralPFJet430;
  Bool_t          HLT_PFHT125;
  Bool_t          HLT_PFHT200;
  Bool_t          HLT_PFHT250;
  Bool_t          HLT_PFHT300;
  Bool_t          HLT_PFHT350;
  Bool_t          HLT_PFHT400;
  Bool_t          HLT_PFHT475;
  Bool_t          HLT_PFHT600;
  Bool_t          HLT_PFHT650;
  Bool_t          HLT_PFHT800;
  Bool_t          HLT_PFHT900;
  Bool_t          HLT_PFHT200_PFAlphaT0p51;
  Bool_t          HLT_PFHT200_DiPFJetAve90_PFAlphaT0p57;
  Bool_t          HLT_PFHT200_DiPFJetAve90_PFAlphaT0p63;
  Bool_t          HLT_PFHT250_DiPFJetAve90_PFAlphaT0p55;
  Bool_t          HLT_PFHT250_DiPFJetAve90_PFAlphaT0p58;
  Bool_t          HLT_PFHT300_DiPFJetAve90_PFAlphaT0p53;
  Bool_t          HLT_PFHT300_DiPFJetAve90_PFAlphaT0p54;
  Bool_t          HLT_PFHT350_DiPFJetAve90_PFAlphaT0p52;
  Bool_t          HLT_PFHT350_DiPFJetAve90_PFAlphaT0p53;
  Bool_t          HLT_PFHT400_DiPFJetAve90_PFAlphaT0p51;
  Bool_t          HLT_PFHT400_DiPFJetAve90_PFAlphaT0p52;
  Bool_t          HLT_MET60_IsoTrk35_Loose;
  Bool_t          HLT_MET75_IsoTrk50;
  Bool_t          HLT_MET90_IsoTrk50;
  Bool_t          HLT_PFMET120_BTagCSV_p067;
  Bool_t          HLT_PFMET120_Mu5;
  Bool_t          HLT_PFMET170_NotCleaned;
  Bool_t          HLT_PFMET170_NoiseCleaned;
  Bool_t          HLT_PFMET170_HBHECleaned;
  Bool_t          HLT_PFMET170_JetIdCleaned;
  Bool_t          HLT_PFMET170_BeamHaloCleaned;
  Bool_t          HLT_PFMET170_HBHE_BeamHaloCleaned;
  Bool_t          HLT_PFMETTypeOne190_HBHE_BeamHaloCleaned;
  Bool_t          HLT_PFMET90_PFMHT90_IDTight;
  Bool_t          HLT_PFMET100_PFMHT100_IDTight;
  Bool_t          HLT_PFMET100_PFMHT100_IDTight_BeamHaloCleaned;
  Bool_t          HLT_PFMET110_PFMHT110_IDTight;
  Bool_t          HLT_PFMET120_PFMHT120_IDTight;
  Bool_t          HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight_BTagCSV_p067;
  Bool_t          HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight;
  Bool_t          HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq200;
  Bool_t          HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq460;
  Bool_t          HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240;
  Bool_t          HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq500;
  Bool_t          HLT_QuadPFJet_VBF;
  Bool_t          HLT_L1_TripleJet_VBF;
  Bool_t          HLT_QuadJet45_TripleBTagCSV_p087;
  Bool_t          HLT_QuadJet45_DoubleBTagCSV_p087;
  Bool_t          HLT_DoubleJet90_Double30_TripleBTagCSV_p087;
  Bool_t          HLT_DoubleJet90_Double30_DoubleBTagCSV_p087;
  Bool_t          HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160;
  Bool_t          HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6;
  Bool_t          HLT_DoubleJetsC112_DoubleBTagCSV_p026_DoublePFJetsC172;
  Bool_t          HLT_DoubleJetsC112_DoubleBTagCSV_p014_DoublePFJetsC112MaxDeta1p6;
  Bool_t          HLT_DoubleJetsC100_SingleBTagCSV_p026;
  Bool_t          HLT_DoubleJetsC100_SingleBTagCSV_p014;
  Bool_t          HLT_DoubleJetsC100_SingleBTagCSV_p026_SinglePFJetC350;
  Bool_t          HLT_DoubleJetsC100_SingleBTagCSV_p014_SinglePFJetC350;
  Bool_t          HLT_Photon135_PFMET100;
  Bool_t          HLT_Photon20_CaloIdVL_IsoL;
  Bool_t          HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_PFMET40;
  Bool_t          HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_VBF;
  Bool_t          HLT_Photon250_NoHE;
  Bool_t          HLT_Photon300_NoHE;
  Bool_t          HLT_Photon26_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon16_AND_HE10_R9Id65_Eta2_Mass60;
  Bool_t          HLT_Photon36_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon22_AND_HE10_R9Id65_Eta2_Mass15;
  Bool_t          HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_PFMET40;
  Bool_t          HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_VBF;
  Bool_t          HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_PFMET40;
  Bool_t          HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_VBF;
  Bool_t          HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_PFMET40;
  Bool_t          HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF;
  Bool_t          HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_PFMET40;
  Bool_t          HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_VBF;
  Bool_t          HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_PFMET40;
  Bool_t          HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_VBF;
  Bool_t          HLT_Mu8_TrkIsoVVL;
  Bool_t          HLT_Mu17_TrkIsoVVL;
  Bool_t          HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30;
  Bool_t          HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30;
  Bool_t          HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30;
  Bool_t          HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30;
  Bool_t          HLT_BTagMu_DiJet20_Mu5;
  Bool_t          HLT_BTagMu_DiJet40_Mu5;
  Bool_t          HLT_BTagMu_DiJet70_Mu5;
  Bool_t          HLT_BTagMu_DiJet110_Mu5;
  Bool_t          HLT_BTagMu_DiJet170_Mu5;
  Bool_t          HLT_BTagMu_Jet300_Mu5;
  Bool_t          HLT_BTagMu_AK8Jet300_Mu5;
  Bool_t          HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
  Bool_t          HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_L1JetTauSeeded;
  Bool_t          HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
  Bool_t          HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL;
  Bool_t          HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL;
  Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;
  Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
  Bool_t          HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;
  Bool_t          HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
  Bool_t          HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;
  Bool_t          HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL;
  Bool_t          HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ;
  Bool_t          HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;
  Bool_t          HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
  Bool_t          HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL;
  Bool_t          HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL;
  Bool_t          HLT_Mu37_Ele27_CaloIdL_GsfTrkIdVL;
  Bool_t          HLT_Mu27_Ele37_CaloIdL_GsfTrkIdVL;
  Bool_t          HLT_Mu8_DiEle12_CaloIdL_TrackIdL;
  Bool_t          HLT_Mu12_Photon25_CaloIdL;
  Bool_t          HLT_Mu12_Photon25_CaloIdL_L1ISO;
  Bool_t          HLT_Mu12_Photon25_CaloIdL_L1OR;
  Bool_t          HLT_Mu17_Photon22_CaloIdL_L1ISO;
  Bool_t          HLT_Mu17_Photon30_CaloIdL_L1ISO;
  Bool_t          HLT_Mu17_Photon35_CaloIdL_L1ISO;
  Bool_t          HLT_DiMu9_Ele9_CaloIdL_TrackIdL;
  Bool_t          HLT_TripleMu_5_3_3;
  Bool_t          HLT_TripleMu_12_10_5;
  Bool_t          HLT_Mu3er_PFHT140_PFMET125;
  Bool_t          HLT_Mu6_PFHT200_PFMET80_BTagCSV_p067;
  Bool_t          HLT_Mu6_PFHT200_PFMET100;
  Bool_t          HLT_Mu14er_PFMET100;
  Bool_t          HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL;
  Bool_t          HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;
  Bool_t          HLT_Ele12_CaloIdL_TrackIdL_IsoVL;
  Bool_t          HLT_Ele17_CaloIdL_GsfTrkIdVL;
  Bool_t          HLT_Ele17_CaloIdL_TrackIdL_IsoVL;
  Bool_t          HLT_Ele23_CaloIdL_TrackIdL_IsoVL;
  Bool_t          HLT_PFHT650_WideJetMJJ900DEtaJJ1p5;
  Bool_t          HLT_PFHT650_WideJetMJJ950DEtaJJ1p5;
  Bool_t          HLT_Photon22;
  Bool_t          HLT_Photon30;
  Bool_t          HLT_Photon36;
  Bool_t          HLT_Photon50;
  Bool_t          HLT_Photon75;
  Bool_t          HLT_Photon90;
  Bool_t          HLT_Photon120;
  Bool_t          HLT_Photon175;
  Bool_t          HLT_Photon165_HE10;
  Bool_t          HLT_Photon22_R9Id90_HE10_IsoM;
  Bool_t          HLT_Photon30_R9Id90_HE10_IsoM;
  Bool_t          HLT_Photon36_R9Id90_HE10_IsoM;
  Bool_t          HLT_Photon50_R9Id90_HE10_IsoM;
  Bool_t          HLT_Photon75_R9Id90_HE10_IsoM;
  Bool_t          HLT_Photon90_R9Id90_HE10_IsoM;
  Bool_t          HLT_Photon120_R9Id90_HE10_IsoM;
  Bool_t          HLT_Photon165_R9Id90_HE10_IsoM;
  Bool_t          HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90;
  Bool_t          HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70;
  Bool_t          HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55;
  Bool_t          HLT_Diphoton30_18_Solid_R9Id_AND_IsoCaloId_AND_HE_R9Id_Mass55;
  Bool_t          HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55;
  Bool_t          HLT_Dimuon0_Jpsi_Muon;
  Bool_t          HLT_Dimuon0_Upsilon_Muon;
  Bool_t          HLT_QuadMuon0_Dimuon0_Jpsi;
  Bool_t          HLT_QuadMuon0_Dimuon0_Upsilon;
  Bool_t          HLT_Rsq0p25_Calo;
  Bool_t          HLT_RsqMR240_Rsq0p09_MR200_4jet_Calo;
  Bool_t          HLT_RsqMR240_Rsq0p09_MR200_Calo;
  Bool_t          HLT_Rsq0p25;
  Bool_t          HLT_Rsq0p30;
  Bool_t          HLT_RsqMR240_Rsq0p09_MR200;
  Bool_t          HLT_RsqMR240_Rsq0p09_MR200_4jet;
  Bool_t          HLT_RsqMR270_Rsq0p09_MR200;
  Bool_t          HLT_RsqMR270_Rsq0p09_MR200_4jet;
  Bool_t          HLT_Rsq0p02_MR300_TriPFJet80_60_40_BTagCSV_p063_p20_Mbb60_200;
  Bool_t          HLT_Rsq0p02_MR400_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200;
  Bool_t          HLT_Rsq0p02_MR450_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200;
  Bool_t          HLT_Rsq0p02_MR500_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200;
  Bool_t          HLT_Rsq0p02_MR550_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200;
  Bool_t          HLT_HT200_DisplacedDijet40_DisplacedTrack;
  Bool_t          HLT_HT250_DisplacedDijet40_DisplacedTrack;
  Bool_t          HLT_HT350_DisplacedDijet40_DisplacedTrack;
  Bool_t          HLT_HT350_DisplacedDijet80_DisplacedTrack;
  Bool_t          HLT_HT350_DisplacedDijet80_Tight_DisplacedTrack;
  Bool_t          HLT_HT350_DisplacedDijet40_Inclusive;
  Bool_t          HLT_HT400_DisplacedDijet40_Inclusive;
  Bool_t          HLT_HT500_DisplacedDijet40_Inclusive;
  Bool_t          HLT_HT550_DisplacedDijet40_Inclusive;
  Bool_t          HLT_HT550_DisplacedDijet80_Inclusive;
  Bool_t          HLT_HT650_DisplacedDijet80_Inclusive;
  Bool_t          HLT_HT750_DisplacedDijet80_Inclusive;
  Bool_t          HLT_VBF_DisplacedJet40_DisplacedTrack;
  Bool_t          HLT_VBF_DisplacedJet40_DisplacedTrack_2TrackIP2DSig5;
  Bool_t          HLT_VBF_DisplacedJet40_TightID_DisplacedTrack;
  Bool_t          HLT_VBF_DisplacedJet40_Hadronic;
  Bool_t          HLT_VBF_DisplacedJet40_Hadronic_2PromptTrack;
  Bool_t          HLT_VBF_DisplacedJet40_TightID_Hadronic;
  Bool_t          HLT_VBF_DisplacedJet40_VTightID_Hadronic;
  Bool_t          HLT_VBF_DisplacedJet40_VVTightID_Hadronic;
  Bool_t          HLT_VBF_DisplacedJet40_VTightID_DisplacedTrack;
  Bool_t          HLT_VBF_DisplacedJet40_VVTightID_DisplacedTrack;
  Bool_t          HLT_PFMETNoMu90_PFMHTNoMu90_IDTight;
  Bool_t          HLT_PFMETNoMu100_PFMHTNoMu100_IDTight;
  Bool_t          HLT_PFMETNoMu110_PFMHTNoMu110_IDTight;
  Bool_t          HLT_PFMETNoMu120_PFMHTNoMu120_IDTight;
  Bool_t          HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_IDTight;
  Bool_t          HLT_MonoCentralPFJet80_PFMETNoMu100_PFMHTNoMu100_IDTight;
  Bool_t          HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight;
  Bool_t          HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight;
  Bool_t          HLT_Ele27_eta2p1_WPLoose_Gsf_HT200;
  Bool_t          HLT_Photon90_CaloIdL_PFHT500;
  Bool_t          HLT_DoubleMu8_Mass8_PFHT250;
  Bool_t          HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT250;
  Bool_t          HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT250;
  Bool_t          HLT_DoubleMu8_Mass8_PFHT300;
  Bool_t          HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300;
  Bool_t          HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300;
  Bool_t          HLT_Mu10_CentralPFJet30_BTagCSV_p13;
  Bool_t          HLT_DoubleMu3_PFMET50;
  Bool_t          HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13;
  Bool_t          HLT_Ele15_IsoVVVL_BTagCSV_p067_PFHT400;
  Bool_t          HLT_Ele15_IsoVVVL_PFHT350_PFMET50;
  Bool_t          HLT_Ele15_IsoVVVL_PFHT600;
  Bool_t          HLT_Ele15_IsoVVVL_PFHT350;
  Bool_t          HLT_Ele15_IsoVVVL_PFHT400_PFMET50;
  Bool_t          HLT_Ele15_IsoVVVL_PFHT400;
  Bool_t          HLT_Ele50_IsoVVVL_PFHT400;
  Bool_t          HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60;
  Bool_t          HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60;
  Bool_t          HLT_Mu15_IsoVVVL_BTagCSV_p067_PFHT400;
  Bool_t          HLT_Mu15_IsoVVVL_PFHT350_PFMET50;
  Bool_t          HLT_Mu15_IsoVVVL_PFHT600;
  Bool_t          HLT_Mu15_IsoVVVL_PFHT350;
  Bool_t          HLT_Mu15_IsoVVVL_PFHT400_PFMET50;
  Bool_t          HLT_Mu15_IsoVVVL_PFHT400;
  Bool_t          HLT_Mu50_IsoVVVL_PFHT400;
  Bool_t          HLT_Dimuon16_Jpsi;
  Bool_t          HLT_Dimuon10_Jpsi_Barrel;
  Bool_t          HLT_Dimuon8_PsiPrime_Barrel;
  Bool_t          HLT_Dimuon8_Upsilon_Barrel;
  Bool_t          HLT_Dimuon0_Phi_Barrel;
  Bool_t          HLT_Mu16_TkMu0_dEta18_Onia;
  Bool_t          HLT_Mu16_TkMu0_dEta18_Phi;
  Bool_t          HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx;
  Bool_t          HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx;
  Bool_t          HLT_Mu8;
  Bool_t          HLT_Mu17;
  Bool_t          HLT_Mu3_PFJet40;
  Bool_t          HLT_Ele8_CaloIdM_TrackIdM_PFJet30;
  Bool_t          HLT_Ele12_CaloIdM_TrackIdM_PFJet30;
  Bool_t          HLT_Ele17_CaloIdM_TrackIdM_PFJet30;
  Bool_t          HLT_Ele23_CaloIdM_TrackIdM_PFJet30;
  Bool_t          HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet140;
  Bool_t          HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165;
  Bool_t          HLT_PFHT400_SixJet30_DoubleBTagCSV_p056;
  Bool_t          HLT_PFHT450_SixJet40_BTagCSV_p056;
  Bool_t          HLT_PFHT400_SixJet30;
  Bool_t          HLT_PFHT450_SixJet40;
  Bool_t          HLT_Ele115_CaloIdVT_GsfTrkIdT;
  Bool_t          HLT_Mu55;
  Bool_t          HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15;
  Bool_t          HLT_Photon90_CaloIdL_PFHT600;
  Bool_t          HLT_PixelTracks_Multiplicity60ForEndOfFill;
  Bool_t          HLT_PixelTracks_Multiplicity85ForEndOfFill;
  Bool_t          HLT_PixelTracks_Multiplicity110ForEndOfFill;
  Bool_t          HLT_PixelTracks_Multiplicity135ForEndOfFill;
  Bool_t          HLT_PixelTracks_Multiplicity160ForEndOfFill;
  Bool_t          HLT_FullTracks_Multiplicity80;
  Bool_t          HLT_FullTracks_Multiplicity100;
  Bool_t          HLT_FullTracks_Multiplicity130;
  Bool_t          HLT_FullTracks_Multiplicity150;
  Bool_t          HLT_ECALHT800;
  Bool_t          HLT_DiSC30_18_EIso_AND_HE_Mass70;
  Bool_t          HLT_Photon125;
  Bool_t          HLT_MET100;
  Bool_t          HLT_MET150;
  Bool_t          HLT_MET200;
  Bool_t          HLT_Ele27_HighEta_Ele20_Mass55;
  Bool_t          HLT_L1FatEvents;
  Bool_t          HLT_Physics;
  Bool_t          HLT_L1FatEvents_part0;
  Bool_t          HLT_L1FatEvents_part1;
  Bool_t          HLT_L1FatEvents_part2;
  Bool_t          HLT_L1FatEvents_part3;
  Bool_t          HLT_Random;
  Bool_t          HLT_ZeroBias;
  Bool_t          HLT_AK4CaloJet30;
  Bool_t          HLT_AK4CaloJet40;
  Bool_t          HLT_AK4CaloJet50;
  Bool_t          HLT_AK4CaloJet80;
  Bool_t          HLT_AK4CaloJet100;
  Bool_t          HLT_AK4PFJet30;
  Bool_t          HLT_AK4PFJet50;
  Bool_t          HLT_AK4PFJet80;
  Bool_t          HLT_AK4PFJet100;
  Bool_t          HLT_HISinglePhoton10;
  Bool_t          HLT_HISinglePhoton15;
  Bool_t          HLT_HISinglePhoton20;
  Bool_t          HLT_HISinglePhoton40;
  Bool_t          HLT_HISinglePhoton60;
  Bool_t          HLT_EcalCalibration;
  Bool_t          HLT_HcalCalibration;
  Bool_t          HLT_GlobalRunHPDNoise;
  Bool_t          HLT_L1BptxMinus;
  Bool_t          HLT_L1BptxPlus;
  Bool_t          HLT_L1NotBptxOR;
  Bool_t          HLT_L1BeamGasMinus;
  Bool_t          HLT_L1BeamGasPlus;
  Bool_t          HLT_L1BptxXOR;
  Bool_t          HLT_L1MinimumBiasHF_OR;
  Bool_t          HLT_L1MinimumBiasHF_AND;
  Bool_t          HLT_HcalNZS;
  Bool_t          HLT_HcalPhiSym;
  Bool_t          HLT_HcalIsolatedbunch;
  Bool_t          HLT_ZeroBias_FirstCollisionAfterAbortGap;
  Bool_t          HLT_ZeroBias_FirstCollisionAfterAbortGap_copy;
  Bool_t          HLT_ZeroBias_FirstCollisionAfterAbortGap_TCDS;
  Bool_t          HLT_ZeroBias_IsolatedBunches;
  Bool_t          HLT_ZeroBias_FirstCollisionInTrain;
  Bool_t          HLT_ZeroBias_FirstBXAfterTrain;
  Bool_t          HLT_Photon500;
  Bool_t          HLT_Photon600;
  Bool_t          HLT_Mu300;
  Bool_t          HLT_Mu350;
  Bool_t          HLT_MET250;
  Bool_t          HLT_MET300;
  Bool_t          HLT_MET600;
  Bool_t          HLT_MET700;
  Bool_t          HLT_PFMET300;
  Bool_t          HLT_PFMET400;
  Bool_t          HLT_PFMET500;
  Bool_t          HLT_PFMET600;
  Bool_t          HLT_Ele250_CaloIdVT_GsfTrkIdT;
  Bool_t          HLT_Ele300_CaloIdVT_GsfTrkIdT;
  Bool_t          HLT_HT2000;
  Bool_t          HLT_HT2500;
  Bool_t          HLT_IsoTrackHE;
  Bool_t          HLT_IsoTrackHB;
  Bool_t          HLTriggerFinalPath;
  Bool_t          Flag_HBHENoiseFilter;
  Bool_t          Flag_HBHENoiseIsoFilter;
  Bool_t          Flag_CSCTightHaloFilter;
  Bool_t          Flag_CSCTightHaloTrkMuUnvetoFilter;
  Bool_t          Flag_CSCTightHalo2015Filter;
  Bool_t          Flag_globalTightHalo2016Filter;
  Bool_t          Flag_globalSuperTightHalo2016Filter;
  Bool_t          Flag_HcalStripHaloFilter;
  Bool_t          Flag_hcalLaserEventFilter;
  Bool_t          Flag_EcalDeadCellTriggerPrimitiveFilter;
  Bool_t          Flag_EcalDeadCellBoundaryEnergyFilter;
  Bool_t          Flag_goodVertices;
  Bool_t          Flag_eeBadScFilter;
  Bool_t          Flag_ecalLaserCorrFilter;
  Bool_t          Flag_trkPOGFilters;
  Bool_t          Flag_chargedHadronTrackResolutionFilter;
  Bool_t          Flag_muonBadTrackFilter;
  Bool_t          Flag_trkPOG_manystripclus53X;
  Bool_t          Flag_trkPOG_toomanystripclus53X;
  Bool_t          Flag_trkPOG_logErrorTooManyClusters;
  Bool_t          Flag_METFilters;

  // List of branches
  TBranch        *b_run;   //!
  TBranch        *b_luminosityBlock;   //!
  TBranch        *b_event;   //!
  TBranch        *b_CaloMET_phi;   //!
  TBranch        *b_CaloMET_pt;   //!
  TBranch        *b_CaloMET_sumEt;   //!
  TBranch        *b_ncmeson;   //!
  TBranch        *b_cmeson_dca;   //!
  TBranch        *b_cmeson_angleXY;   //!
  TBranch        *b_cmeson_angleXYZ;   //!
  TBranch        *b_cmeson_trk_normalizedChi2;   //!
  TBranch        *b_cmeson_trk_nHits;   //!
  TBranch        *b_cmeson_trk_pt;   //!
  TBranch        *b_cmeson_trk_ipsigXY;   //!
  TBranch        *b_cmeson_trk_ipsigZ;   //!
  TBranch        *b_cmeson_lxy;   //!
  TBranch        *b_cmeson_lxySig;   //!
  TBranch        *b_cmeson_l3D;   //!
  TBranch        *b_cmeson_l3DSig;   //!
  TBranch        *b_cmeson_jetDR;   //!
  TBranch        *b_cmeson_legDR;   //!
  TBranch        *b_cmeson_diffMass;   //!
  TBranch        *b_cmeson_nJet;   //!
  TBranch        *b_cmeson_mcMatch;   //!
  TBranch        *b_nElectron;   //!
  TBranch        *b_Electron_deltaEtaSC;   //!
  TBranch        *b_Electron_dr03EcalRecHitSumEt;   //!
  TBranch        *b_Electron_dr03HcalDepth1TowerSumEt;   //!
  TBranch        *b_Electron_dr03TkSumPt;   //!
  TBranch        *b_Electron_dxy;   //!
  TBranch        *b_Electron_dxyErr;   //!
  TBranch        *b_Electron_dz;   //!
  TBranch        *b_Electron_dzErr;   //!
  TBranch        *b_Electron_eCorr;   //!
  TBranch        *b_Electron_eInvMinusPInv;   //!
  TBranch        *b_Electron_energyErr;   //!
  TBranch        *b_Electron_eta;   //!
  TBranch        *b_Electron_hoe;   //!
  TBranch        *b_Electron_ip3d;   //!
  TBranch        *b_Electron_mass;   //!
  TBranch        *b_Electron_miniPFRelIso_all;   //!
  TBranch        *b_Electron_miniPFRelIso_chg;   //!
  TBranch        *b_Electron_mvaSpring16GP;   //!
  TBranch        *b_Electron_mvaSpring16HZZ;   //!
  TBranch        *b_Electron_pfRelIso03_all;   //!
  TBranch        *b_Electron_pfRelIso03_chg;   //!
  TBranch        *b_Electron_phi;   //!
  TBranch        *b_Electron_pt;   //!
  TBranch        *b_Electron_r9;   //!
  TBranch        *b_Electron_sieie;   //!
  TBranch        *b_Electron_sip3d;   //!
  TBranch        *b_Electron_mvaTTH;   //!
  TBranch        *b_Electron_charge;   //!
  TBranch        *b_Electron_cutBased;   //!
  TBranch        *b_Electron_cutBased_HLTPreSel;   //!
  TBranch        *b_Electron_jetIdx;   //!
  TBranch        *b_Electron_pdgId;   //!
  TBranch        *b_Electron_photonIdx;   //!
  TBranch        *b_Electron_tightCharge;   //!
  TBranch        *b_Electron_vidNestedWPBitmap;   //!
  TBranch        *b_Electron_convVeto;   //!
  TBranch        *b_Electron_cutBased_HEEP;   //!
  TBranch        *b_Electron_isPFcand;   //!
  TBranch        *b_Electron_lostHits;   //!
  TBranch        *b_Electron_mvaSpring16GP_WP80;   //!
  TBranch        *b_Electron_mvaSpring16GP_WP90;   //!
  TBranch        *b_Electron_mvaSpring16HZZ_WPL;   //!
  TBranch        *b_nFatJet;   //!
  TBranch        *b_FatJet_area;   //!
  TBranch        *b_FatJet_btagCMVA;   //!
  TBranch        *b_FatJet_btagCSVV2;   //!
  TBranch        *b_FatJet_btagDeepB;   //!
  TBranch        *b_FatJet_btagHbb;   //!
  TBranch        *b_FatJet_eta;   //!
  TBranch        *b_FatJet_mass;   //!
  TBranch        *b_FatJet_msoftdrop;   //!
  TBranch        *b_FatJet_n2b1;   //!
  TBranch        *b_FatJet_n3b1;   //!
  TBranch        *b_FatJet_phi;   //!
  TBranch        *b_FatJet_pt;   //!
  TBranch        *b_FatJet_tau1;   //!
  TBranch        *b_FatJet_tau2;   //!
  TBranch        *b_FatJet_tau3;   //!
  TBranch        *b_FatJet_tau4;   //!
  TBranch        *b_FatJet_subJetIdx1;   //!
  TBranch        *b_FatJet_subJetIdx2;   //!
  TBranch        *b_nGenJetAK8;   //!
  TBranch        *b_GenJetAK8_eta;   //!
  TBranch        *b_GenJetAK8_mass;   //!
  TBranch        *b_GenJetAK8_phi;   //!
  TBranch        *b_GenJetAK8_pt;   //!
  TBranch        *b_nGenJet;   //!
  TBranch        *b_GenJet_eta;   //!
  TBranch        *b_GenJet_mass;   //!
  TBranch        *b_GenJet_phi;   //!
  TBranch        *b_GenJet_pt;   //!
  TBranch        *b_nGenPart;   //!
  TBranch        *b_GenPart_eta;   //!
  TBranch        *b_GenPart_mass;   //!
  TBranch        *b_GenPart_phi;   //!
  TBranch        *b_GenPart_pt;   //!
  TBranch        *b_GenPart_genPartIdxMother;   //!
  TBranch        *b_GenPart_pdgId;   //!
  TBranch        *b_GenPart_status;   //!
  TBranch        *b_GenPart_statusFlags;   //!
  TBranch        *b_Generator_scalePDF;   //!
  TBranch        *b_Generator_x1;   //!
  TBranch        *b_Generator_x2;   //!
  TBranch        *b_Generator_xpdf1;   //!
  TBranch        *b_Generator_xpdf2;   //!
  TBranch        *b_Generator_id1;   //!
  TBranch        *b_Generator_id2;   //!
  TBranch        *b_nGenVisTau;   //!
  TBranch        *b_GenVisTau_eta;   //!
  TBranch        *b_GenVisTau_mass;   //!
  TBranch        *b_GenVisTau_phi;   //!
  TBranch        *b_GenVisTau_pt;   //!
  TBranch        *b_GenVisTau_charge;   //!
  TBranch        *b_GenVisTau_genPartIdxMother;   //!
  TBranch        *b_GenVisTau_status;   //!
  TBranch        *b_genWeight;   //!
  TBranch        *b_LHEWeight_originalXWGTUP;   //!
  TBranch        *b_nLHEPdfWeight;   //!
  TBranch        *b_LHEPdfWeight;   //!
  TBranch        *b_nLHEScaleWeight;   //!
  TBranch        *b_LHEScaleWeight;   //!
  TBranch        *b_nJet;   //!
  TBranch        *b_Jet_area;   //!
  TBranch        *b_Jet_btagCMVA;   //!
  TBranch        *b_Jet_btagCSVV2;   //!
  TBranch        *b_Jet_btagDeepB;   //!
  TBranch        *b_Jet_btagDeepC;   //!
  TBranch        *b_Jet_chEmEF;   //!
  TBranch        *b_Jet_chHEF;   //!
  TBranch        *b_Jet_eta;   //!
  TBranch        *b_Jet_mass;   //!
  TBranch        *b_Jet_neEmEF;   //!
  TBranch        *b_Jet_neHEF;   //!
  TBranch        *b_Jet_phi;   //!
  TBranch        *b_Jet_pt;   //!
  TBranch        *b_Jet_qgl;   //!
  TBranch        *b_Jet_rawFactor;   //!
  TBranch        *b_Jet_bReg;   //!
  TBranch        *b_Jet_electronIdx1;   //!
  TBranch        *b_Jet_electronIdx2;   //!
  TBranch        *b_Jet_jetId;   //!
  TBranch        *b_Jet_muonIdx1;   //!
  TBranch        *b_Jet_muonIdx2;   //!
  TBranch        *b_Jet_nConstituents;   //!
  TBranch        *b_Jet_nElectrons;   //!
  TBranch        *b_Jet_nMuons;   //!
  TBranch        *b_Jet_puId;   //!
  TBranch        *b_LHE_HT;   //!
  TBranch        *b_LHE_HTIncoming;   //!
  TBranch        *b_LHE_Vpt;   //!
  TBranch        *b_LHE_Njets;   //!
  TBranch        *b_LHE_Nb;   //!
  TBranch        *b_LHE_Nc;   //!
  TBranch        *b_LHE_Nuds;   //!
  TBranch        *b_LHE_Nglu;   //!
  TBranch        *b_LHE_NpNLO;   //!
  TBranch        *b_LHE_NpLO;   //!
  TBranch        *b_GenMET_phi;   //!
  TBranch        *b_GenMET_pt;   //!
  TBranch        *b_MET_MetUnclustEnUpDeltaX;   //!
  TBranch        *b_MET_MetUnclustEnUpDeltaY;   //!
  TBranch        *b_MET_covXX;   //!
  TBranch        *b_MET_covXY;   //!
  TBranch        *b_MET_covYY;   //!
  TBranch        *b_MET_phi;   //!
  TBranch        *b_MET_pt;   //!
  TBranch        *b_MET_significance;   //!
  TBranch        *b_MET_sumEt;   //!
  TBranch        *b_nMuon;   //!
  TBranch        *b_Muon_dxy;   //!
  TBranch        *b_Muon_dxyErr;   //!
  TBranch        *b_Muon_dz;   //!
  TBranch        *b_Muon_dzErr;   //!
  TBranch        *b_Muon_eta;   //!
  TBranch        *b_Muon_ip3d;   //!
  TBranch        *b_Muon_mass;   //!
  TBranch        *b_Muon_miniPFRelIso_all;   //!
  TBranch        *b_Muon_miniPFRelIso_chg;   //!
  TBranch        *b_Muon_pfRelIso03_all;   //!
  TBranch        *b_Muon_pfRelIso03_chg;   //!
  TBranch        *b_Muon_pfRelIso04_all;   //!
  TBranch        *b_Muon_phi;   //!
  TBranch        *b_Muon_pt;   //!
  TBranch        *b_Muon_ptErr;   //!
  TBranch        *b_Muon_segmentComp;   //!
  TBranch        *b_Muon_sip3d;   //!
  TBranch        *b_Muon_mvaTTH;   //!
  TBranch        *b_Muon_charge;   //!
  TBranch        *b_Muon_jetIdx;   //!
  TBranch        *b_Muon_nStations;   //!
  TBranch        *b_Muon_nTrackerLayers;   //!
  TBranch        *b_Muon_pdgId;   //!
  TBranch        *b_Muon_tightCharge;   //!
  TBranch        *b_Muon_globalMu;   //!
  TBranch        *b_Muon_highPtId;   //!
  TBranch        *b_Muon_isPFcand;   //!
  TBranch        *b_Muon_mediumId;   //!
  TBranch        *b_Muon_softId;   //!
  TBranch        *b_Muon_tightId;   //!
  TBranch        *b_Muon_trackerMu;   //!
  TBranch        *b_nPhoton;   //!
  TBranch        *b_Photon_eCorr;   //!
  TBranch        *b_Photon_energyErr;   //!
  TBranch        *b_Photon_eta;   //!
  TBranch        *b_Photon_hoe;   //!
  TBranch        *b_Photon_mass;   //!
  TBranch        *b_Photon_mvaID;   //!
  TBranch        *b_Photon_pfRelIso03_all;   //!
  TBranch        *b_Photon_pfRelIso03_chg;   //!
  TBranch        *b_Photon_phi;   //!
  TBranch        *b_Photon_pt;   //!
  TBranch        *b_Photon_r9;   //!
  TBranch        *b_Photon_sieie;   //!
  TBranch        *b_Photon_charge;   //!
  TBranch        *b_Photon_cutBased;   //!
  TBranch        *b_Photon_electronIdx;   //!
  TBranch        *b_Photon_jetIdx;   //!
  TBranch        *b_Photon_pdgId;   //!
  TBranch        *b_Photon_vidNestedWPBitmap;   //!
  TBranch        *b_Photon_electronVeto;   //!
  TBranch        *b_Photon_mvaID_WP80;   //!
  TBranch        *b_Photon_mvaID_WP90;   //!
  TBranch        *b_Photon_pixelSeed;   //!
  TBranch        *b_Pileup_nPU;   //!
  TBranch        *b_Pileup_nTrueInt;   //!
  TBranch        *b_PuppiMET_phi;   //!
  TBranch        *b_PuppiMET_pt;   //!
  TBranch        *b_PuppiMET_sumEt;   //!
  TBranch        *b_RawMET_phi;   //!
  TBranch        *b_RawMET_pt;   //!
  TBranch        *b_RawMET_sumEt;   //!
  TBranch        *b_fixedGridRhoFastjetAll;   //!
  TBranch        *b_fixedGridRhoFastjetCentralCalo;   //!
  TBranch        *b_fixedGridRhoFastjetCentralNeutral;   //!
  TBranch        *b_nGenDressedLepton;   //!
  TBranch        *b_GenDressedLepton_eta;   //!
  TBranch        *b_GenDressedLepton_mass;   //!
  TBranch        *b_GenDressedLepton_phi;   //!
  TBranch        *b_GenDressedLepton_pt;   //!
  TBranch        *b_GenDressedLepton_pdgId;   //!
  TBranch        *b_nSoftActivityJet;   //!
  TBranch        *b_SoftActivityJet_eta;   //!
  TBranch        *b_SoftActivityJet_phi;   //!
  TBranch        *b_SoftActivityJet_pt;   //!
  TBranch        *b_SoftActivityJetHT;   //!
  TBranch        *b_SoftActivityJetHT10;   //!
  TBranch        *b_SoftActivityJetHT2;   //!
  TBranch        *b_SoftActivityJetHT5;   //!
  TBranch        *b_SoftActivityJetNjets10;   //!
  TBranch        *b_SoftActivityJetNjets2;   //!
  TBranch        *b_SoftActivityJetNjets5;   //!
  TBranch        *b_nSubJet;   //!
  TBranch        *b_SubJet_btagCMVA;   //!
  TBranch        *b_SubJet_btagCSVV2;   //!
  TBranch        *b_SubJet_btagDeepB;   //!
  TBranch        *b_SubJet_eta;   //!
  TBranch        *b_SubJet_mass;   //!
  TBranch        *b_SubJet_n2b1;   //!
  TBranch        *b_SubJet_n3b1;   //!
  TBranch        *b_SubJet_phi;   //!
  TBranch        *b_SubJet_pt;   //!
  TBranch        *b_SubJet_tau1;   //!
  TBranch        *b_SubJet_tau2;   //!
  TBranch        *b_SubJet_tau3;   //!
  TBranch        *b_SubJet_tau4;   //!
  TBranch        *b_nTau;   //!
  TBranch        *b_Tau_chargedIso;   //!
  TBranch        *b_Tau_dxy;   //!
  TBranch        *b_Tau_dz;   //!
  TBranch        *b_Tau_eta;   //!
  TBranch        *b_Tau_footprintCorr;   //!
  TBranch        *b_Tau_leadTkDeltaEta;   //!
  TBranch        *b_Tau_leadTkDeltaPhi;   //!
  TBranch        *b_Tau_leadTkPtOverTauPt;   //!
  TBranch        *b_Tau_mass;   //!
  TBranch        *b_Tau_neutralIso;   //!
  TBranch        *b_Tau_phi;   //!
  TBranch        *b_Tau_photonsOutsideSignalCone;   //!
  TBranch        *b_Tau_pt;   //!
  TBranch        *b_Tau_puCorr;   //!
  TBranch        *b_Tau_rawAntiEle;   //!
  TBranch        *b_Tau_rawIso;   //!
  TBranch        *b_Tau_rawMVAnewDM;   //!
  TBranch        *b_Tau_rawMVAoldDM;   //!
  TBranch        *b_Tau_rawMVAoldDMdR03;   //!
  TBranch        *b_Tau_charge;   //!
  TBranch        *b_Tau_decayMode;   //!
  TBranch        *b_Tau_jetIdx;   //!
  TBranch        *b_Tau_rawAntiEleCat;   //!
  TBranch        *b_Tau_idAntiEle;   //!
  TBranch        *b_Tau_idAntiMu;   //!
  TBranch        *b_Tau_idDecayMode;   //!
  TBranch        *b_Tau_idDecayModeNewDMs;   //!
  TBranch        *b_Tau_idMVAnewDM;   //!
  TBranch        *b_Tau_idMVAoldDM;   //!
  TBranch        *b_Tau_idMVAoldDMdR03;   //!
  TBranch        *b_TkMET_phi;   //!
  TBranch        *b_TkMET_pt;   //!
  TBranch        *b_TkMET_sumEt;   //!
  TBranch        *b_nTrigObj;   //!
  TBranch        *b_TrigObj_pt;   //!
  TBranch        *b_TrigObj_eta;   //!
  TBranch        *b_TrigObj_phi;   //!
  TBranch        *b_TrigObj_l1pt;   //!
  TBranch        *b_TrigObj_l1pt_2;   //!
  TBranch        *b_TrigObj_l2pt;   //!
  TBranch        *b_TrigObj_id;   //!
  TBranch        *b_TrigObj_l1iso;   //!
  TBranch        *b_TrigObj_l1charge;   //!
  TBranch        *b_TrigObj_filterBits;   //!
  TBranch        *b_genTtbarId;   //!
  TBranch        *b_nOtherPV;   //!
  TBranch        *b_OtherPV_z;   //!
  TBranch        *b_PV_ndof;   //!
  TBranch        *b_PV_x;   //!
  TBranch        *b_PV_y;   //!
  TBranch        *b_PV_z;   //!
  TBranch        *b_PV_chi2;   //!
  TBranch        *b_PV_score;   //!
  TBranch        *b_PV_npvs;   //!
  TBranch        *b_nSV;   //!
  TBranch        *b_SV_dlen;   //!
  TBranch        *b_SV_dlenSig;   //!
  TBranch        *b_SV_pAngle;   //!
  TBranch        *b_cmeson_chi2;   //!
  TBranch        *b_cmeson_eta;   //!
  TBranch        *b_cmeson_mass;   //!
  TBranch        *b_cmeson_phi;   //!
  TBranch        *b_cmeson_pt;   //!
  TBranch        *b_cmeson_x;   //!
  TBranch        *b_cmeson_y;   //!
  TBranch        *b_cmeson_z;   //!
  TBranch        *b_cmeson_ndof;   //!
  TBranch        *b_cmeson_pdgId;   //!
  TBranch        *b_Electron_genPartIdx;   //!
  TBranch        *b_Electron_genPartFlav;   //!
  TBranch        *b_GenJetAK8_partonFlavour;   //!
  TBranch        *b_GenJetAK8_hadronFlavour;   //!
  TBranch        *b_GenJet_partonFlavour;   //!
  TBranch        *b_GenJet_hadronFlavour;   //!
  TBranch        *b_Jet_genJetIdx;   //!
  TBranch        *b_Jet_hadronFlavour;   //!
  TBranch        *b_Jet_partonFlavour;   //!
  TBranch        *b_Muon_genPartIdx;   //!
  TBranch        *b_Muon_genPartFlav;   //!
  TBranch        *b_Photon_genPartIdx;   //!
  TBranch        *b_Photon_genPartFlav;   //!
  TBranch        *b_MET_fiducialGenPhi;   //!
  TBranch        *b_MET_fiducialGenPt;   //!
  TBranch        *b_Electron_cleanmask;   //!
  TBranch        *b_Jet_cleanmask;   //!
  TBranch        *b_Muon_cleanmask;   //!
  TBranch        *b_Photon_cleanmask;   //!
  TBranch        *b_Tau_cleanmask;   //!
  TBranch        *b_SV_chi2;   //!
  TBranch        *b_SV_eta;   //!
  TBranch        *b_SV_mass;   //!
  TBranch        *b_SV_ndof;   //!
  TBranch        *b_SV_phi;   //!
  TBranch        *b_SV_pt;   //!
  TBranch        *b_SV_x;   //!
  TBranch        *b_SV_y;   //!
  TBranch        *b_SV_z;   //!
  TBranch        *b_Tau_genPartIdx;   //!
  TBranch        *b_Tau_genPartFlav;   //!
  TBranch        *b_L1simulation_step;   //!
  TBranch        *b_HLTriggerFirstPath;   //!
  TBranch        *b_HLT_AK8PFJet360_TrimMass30;   //!
  TBranch        *b_HLT_AK8PFJet400_TrimMass30;   //!
  TBranch        *b_HLT_AK8PFHT750_TrimMass50;   //!
  TBranch        *b_HLT_AK8PFHT800_TrimMass50;   //!
  TBranch        *b_HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20;   //!
  TBranch        *b_HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087;   //!
  TBranch        *b_HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087;   //!
  TBranch        *b_HLT_AK8DiPFJet300_200_TrimMass30;   //!
  TBranch        *b_HLT_AK8PFHT700_TrimR0p1PT0p03Mass50;   //!
  TBranch        *b_HLT_AK8PFHT650_TrimR0p1PT0p03Mass50;   //!
  TBranch        *b_HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20;   //!
  TBranch        *b_HLT_AK8DiPFJet280_200_TrimMass30;   //!
  TBranch        *b_HLT_AK8DiPFJet250_200_TrimMass30;   //!
  TBranch        *b_HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20;   //!
  TBranch        *b_HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20;   //!
  TBranch        *b_HLT_CaloJet260;   //!
  TBranch        *b_HLT_CaloJet500_NoJetID;   //!
  TBranch        *b_HLT_Dimuon13_PsiPrime;   //!
  TBranch        *b_HLT_Dimuon13_Upsilon;   //!
  TBranch        *b_HLT_Dimuon20_Jpsi;   //!
  TBranch        *b_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf;   //!
  TBranch        *b_HLT_DoubleEle25_CaloIdL_GsfTrkIdVL;   //!
  TBranch        *b_HLT_DoubleEle33_CaloIdL;   //!
  TBranch        *b_HLT_DoubleEle33_CaloIdL_MW;   //!
  TBranch        *b_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW;   //!
  TBranch        *b_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL;   //!
  TBranch        *b_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg;   //!
  TBranch        *b_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg;   //!
  TBranch        *b_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg;   //!
  TBranch        *b_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg;   //!
  TBranch        *b_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1;   //!
  TBranch        *b_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1;   //!
  TBranch        *b_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg;   //!
  TBranch        *b_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg;   //!
  TBranch        *b_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1;   //!
  TBranch        *b_HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL;   //!
  TBranch        *b_HLT_DoubleMu33NoFiltersNoVtx;   //!
  TBranch        *b_HLT_DoubleMu38NoFiltersNoVtx;   //!
  TBranch        *b_HLT_DoubleMu23NoFiltersNoVtxDisplaced;   //!
  TBranch        *b_HLT_DoubleMu28NoFiltersNoVtxDisplaced;   //!
  TBranch        *b_HLT_DoubleMu0;   //!
  TBranch        *b_HLT_DoubleMu4_3_Bs;   //!
  TBranch        *b_HLT_DoubleMu4_3_Jpsi_Displaced;   //!
  TBranch        *b_HLT_DoubleMu4_JpsiTrk_Displaced;   //!
  TBranch        *b_HLT_DoubleMu4_LowMassNonResonantTrk_Displaced;   //!
  TBranch        *b_HLT_DoubleMu3_Trk_Tau3mu;   //!
  TBranch        *b_HLT_DoubleMu4_PsiPrimeTrk_Displaced;   //!
  TBranch        *b_HLT_Mu7p5_L2Mu2_Jpsi;   //!
  TBranch        *b_HLT_Mu7p5_L2Mu2_Upsilon;   //!
  TBranch        *b_HLT_Mu7p5_Track2_Jpsi;   //!
  TBranch        *b_HLT_Mu7p5_Track3p5_Jpsi;   //!
  TBranch        *b_HLT_Mu7p5_Track7_Jpsi;   //!
  TBranch        *b_HLT_Mu7p5_Track2_Upsilon;   //!
  TBranch        *b_HLT_Mu7p5_Track3p5_Upsilon;   //!
  TBranch        *b_HLT_Mu7p5_Track7_Upsilon;   //!
  TBranch        *b_HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing;   //!
  TBranch        *b_HLT_Dimuon0er16_Jpsi_NoVertexing;   //!
  TBranch        *b_HLT_Dimuon6_Jpsi_NoVertexing;   //!
  TBranch        *b_HLT_Photon150;   //!
  TBranch        *b_HLT_Photon90_CaloIdL_HT300;   //!
  TBranch        *b_HLT_HT250_CaloMET70;   //!
  TBranch        *b_HLT_DoublePhoton60;   //!
  TBranch        *b_HLT_DoublePhoton85;   //!
  TBranch        *b_HLT_Ele17_Ele8_Gsf;   //!
  TBranch        *b_HLT_Ele20_eta2p1_WPLoose_Gsf_LooseIsoPFTau28;   //!
  TBranch        *b_HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau29;   //!
  TBranch        *b_HLT_Ele22_eta2p1_WPLoose_Gsf;   //!
  TBranch        *b_HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1;   //!
  TBranch        *b_HLT_Ele23_WPLoose_Gsf;   //!
  TBranch        *b_HLT_Ele23_WPLoose_Gsf_WHbbBoost;   //!
  TBranch        *b_HLT_Ele24_eta2p1_WPLoose_Gsf;   //!
  TBranch        *b_HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20;   //!
  TBranch        *b_HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1;   //!
  TBranch        *b_HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau30;   //!
  TBranch        *b_HLT_Ele25_WPTight_Gsf;   //!
  TBranch        *b_HLT_Ele25_eta2p1_WPLoose_Gsf;   //!
  TBranch        *b_HLT_Ele25_eta2p1_WPTight_Gsf;   //!
  TBranch        *b_HLT_Ele27_WPLoose_Gsf;   //!
  TBranch        *b_HLT_Ele27_WPLoose_Gsf_WHbbBoost;   //!
  TBranch        *b_HLT_Ele27_WPTight_Gsf;   //!
  TBranch        *b_HLT_Ele27_WPTight_Gsf_L1JetTauSeeded;   //!
  TBranch        *b_HLT_Ele27_eta2p1_WPLoose_Gsf;   //!
  TBranch        *b_HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1;   //!
  TBranch        *b_HLT_Ele27_eta2p1_WPTight_Gsf;   //!
  TBranch        *b_HLT_Ele30_WPTight_Gsf;   //!
  TBranch        *b_HLT_Ele30_eta2p1_WPLoose_Gsf;   //!
  TBranch        *b_HLT_Ele30_eta2p1_WPTight_Gsf;   //!
  TBranch        *b_HLT_Ele32_WPTight_Gsf;   //!
  TBranch        *b_HLT_Ele32_eta2p1_WPLoose_Gsf;   //!
  TBranch        *b_HLT_Ele32_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1;   //!
  TBranch        *b_HLT_Ele32_eta2p1_WPTight_Gsf;   //!
  TBranch        *b_HLT_Ele35_WPLoose_Gsf;   //!
  TBranch        *b_HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50;   //!
  TBranch        *b_HLT_Ele36_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1;   //!
  TBranch        *b_HLT_Ele45_WPLoose_Gsf;   //!
  TBranch        *b_HLT_Ele45_WPLoose_Gsf_L1JetTauSeeded;   //!
  TBranch        *b_HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50;   //!
  TBranch        *b_HLT_Ele105_CaloIdVT_GsfTrkIdT;   //!
  TBranch        *b_HLT_Ele30WP60_SC4_Mass55;   //!
  TBranch        *b_HLT_Ele30WP60_Ele8_Mass55;   //!
  TBranch        *b_HLT_HT200;   //!
  TBranch        *b_HLT_HT275;   //!
  TBranch        *b_HLT_HT325;   //!
  TBranch        *b_HLT_HT425;   //!
  TBranch        *b_HLT_HT575;   //!
  TBranch        *b_HLT_HT410to430;   //!
  TBranch        *b_HLT_HT430to450;   //!
  TBranch        *b_HLT_HT450to470;   //!
  TBranch        *b_HLT_HT470to500;   //!
  TBranch        *b_HLT_HT500to550;   //!
  TBranch        *b_HLT_HT550to650;   //!
  TBranch        *b_HLT_HT650;   //!
  TBranch        *b_HLT_Mu16_eta2p1_MET30;   //!
  TBranch        *b_HLT_IsoMu16_eta2p1_MET30;   //!
  TBranch        *b_HLT_IsoMu16_eta2p1_MET30_LooseIsoPFTau50_Trk30_eta2p1;   //!
  TBranch        *b_HLT_IsoMu17_eta2p1;   //!
  TBranch        *b_HLT_IsoMu17_eta2p1_LooseIsoPFTau20;   //!
  TBranch        *b_HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1;   //!
  TBranch        *b_HLT_DoubleIsoMu17_eta2p1;   //!
  TBranch        *b_HLT_DoubleIsoMu17_eta2p1_noDzCut;   //!
  TBranch        *b_HLT_IsoMu18;   //!
  TBranch        *b_HLT_IsoMu19_eta2p1_LooseIsoPFTau20;   //!
  TBranch        *b_HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1;   //!
  TBranch        *b_HLT_IsoMu19_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg;   //!
  TBranch        *b_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20;   //!
  TBranch        *b_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg;   //!
  TBranch        *b_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg;   //!
  TBranch        *b_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg;   //!
  TBranch        *b_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg;   //!
  TBranch        *b_HLT_IsoMu20;   //!
  TBranch        *b_HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1;   //!
  TBranch        *b_HLT_IsoMu21_eta2p1_LooseIsoPFTau50_Trk30_eta2p1_SingleL1;   //!
  TBranch        *b_HLT_IsoMu21_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg;   //!
  TBranch        *b_HLT_IsoMu22;   //!
  TBranch        *b_HLT_IsoMu22_eta2p1;   //!
  TBranch        *b_HLT_IsoMu24;   //!
  TBranch        *b_HLT_IsoMu27;   //!
  TBranch        *b_HLT_IsoTkMu18;   //!
  TBranch        *b_HLT_IsoTkMu20;   //!
  TBranch        *b_HLT_IsoTkMu22;   //!
  TBranch        *b_HLT_IsoTkMu22_eta2p1;   //!
  TBranch        *b_HLT_IsoTkMu24;   //!
  TBranch        *b_HLT_IsoTkMu27;   //!
  TBranch        *b_HLT_JetE30_NoBPTX3BX;   //!
  TBranch        *b_HLT_JetE30_NoBPTX;   //!
  TBranch        *b_HLT_JetE50_NoBPTX3BX;   //!
  TBranch        *b_HLT_JetE70_NoBPTX3BX;   //!
  TBranch        *b_HLT_L1SingleMu18;   //!
  TBranch        *b_HLT_L2Mu10;   //!
  TBranch        *b_HLT_L1SingleMuOpen;   //!
  TBranch        *b_HLT_L1SingleMuOpen_DT;   //!
  TBranch        *b_HLT_L2DoubleMu23_NoVertex;   //!
  TBranch        *b_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10;   //!
  TBranch        *b_HLT_L2DoubleMu38_NoVertex_2Cha_Angle2p5_Mass10;   //!
  TBranch        *b_HLT_L2Mu10_NoVertex_NoBPTX3BX;   //!
  TBranch        *b_HLT_L2Mu10_NoVertex_NoBPTX;   //!
  TBranch        *b_HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX;   //!
  TBranch        *b_HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX;   //!
  TBranch        *b_HLT_LooseIsoPFTau50_Trk30_eta2p1;   //!
  TBranch        *b_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80;   //!
  TBranch        *b_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90;   //!
  TBranch        *b_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110;   //!
  TBranch        *b_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120;   //!
  TBranch        *b_HLT_PFTau120_eta2p1;   //!
  TBranch        *b_HLT_PFTau140_eta2p1;   //!
  TBranch        *b_HLT_VLooseIsoPFTau120_Trk50_eta2p1;   //!
  TBranch        *b_HLT_VLooseIsoPFTau140_Trk50_eta2p1;   //!
  TBranch        *b_HLT_Mu17_Mu8;   //!
  TBranch        *b_HLT_Mu17_Mu8_DZ;   //!
  TBranch        *b_HLT_Mu17_Mu8_SameSign;   //!
  TBranch        *b_HLT_Mu17_Mu8_SameSign_DZ;   //!
  TBranch        *b_HLT_Mu20_Mu10;   //!
  TBranch        *b_HLT_Mu20_Mu10_DZ;   //!
  TBranch        *b_HLT_Mu20_Mu10_SameSign;   //!
  TBranch        *b_HLT_Mu20_Mu10_SameSign_DZ;   //!
  TBranch        *b_HLT_Mu17_TkMu8_DZ;   //!
  TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;   //!
  TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;   //!
  TBranch        *b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL;   //!
  TBranch        *b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;   //!
  TBranch        *b_HLT_Mu25_TkMu0_dEta18_Onia;   //!
  TBranch        *b_HLT_Mu27_TkMu8;   //!
  TBranch        *b_HLT_Mu30_TkMu11;   //!
  TBranch        *b_HLT_Mu30_eta2p1_PFJet150_PFJet50;   //!
  TBranch        *b_HLT_Mu40_TkMu11;   //!
  TBranch        *b_HLT_Mu40_eta2p1_PFJet200_PFJet50;   //!
  TBranch        *b_HLT_Mu20;   //!
  TBranch        *b_HLT_TkMu17;   //!
  TBranch        *b_HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL;   //!
  TBranch        *b_HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;   //!
  TBranch        *b_HLT_TkMu20;   //!
  TBranch        *b_HLT_Mu24_eta2p1;   //!
  TBranch        *b_HLT_TkMu24_eta2p1;   //!
  TBranch        *b_HLT_Mu27;   //!
  TBranch        *b_HLT_TkMu27;   //!
  TBranch        *b_HLT_Mu45_eta2p1;   //!
  TBranch        *b_HLT_Mu50;   //!
  TBranch        *b_HLT_TkMu50;   //!
  TBranch        *b_HLT_Mu38NoFiltersNoVtx_Photon38_CaloIdL;   //!
  TBranch        *b_HLT_Mu42NoFiltersNoVtx_Photon42_CaloIdL;   //!
  TBranch        *b_HLT_Mu28NoFiltersNoVtxDisplaced_Photon28_CaloIdL;   //!
  TBranch        *b_HLT_Mu33NoFiltersNoVtxDisplaced_Photon33_CaloIdL;   //!
  TBranch        *b_HLT_Mu23NoFiltersNoVtx_Photon23_CaloIdL;   //!
  TBranch        *b_HLT_DoubleMu18NoFiltersNoVtx;   //!
  TBranch        *b_HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Tight;   //!
  TBranch        *b_HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Loose;   //!
  TBranch        *b_HLT_Mu28NoFiltersNoVtx_DisplacedJet40_Loose;   //!
  TBranch        *b_HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Tight;   //!
  TBranch        *b_HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Loose;   //!
  TBranch        *b_HLT_Mu38NoFiltersNoVtx_DisplacedJet60_Loose;   //!
  TBranch        *b_HLT_Mu28NoFiltersNoVtx_CentralCaloJet40;   //!
  TBranch        *b_HLT_PFHT300_PFMET100;   //!
  TBranch        *b_HLT_PFHT300_PFMET110;   //!
  TBranch        *b_HLT_PFHT550_4JetPt50;   //!
  TBranch        *b_HLT_PFHT650_4JetPt50;   //!
  TBranch        *b_HLT_PFHT750_4JetPt50;   //!
  TBranch        *b_HLT_PFHT750_4JetPt70;   //!
  TBranch        *b_HLT_PFHT750_4JetPt80;   //!
  TBranch        *b_HLT_PFHT800_4JetPt50;   //!
  TBranch        *b_HLT_PFHT850_4JetPt50;   //!
  TBranch        *b_HLT_PFJet15_NoCaloMatched;   //!
  TBranch        *b_HLT_PFJet25_NoCaloMatched;   //!
  TBranch        *b_HLT_DiPFJet15_NoCaloMatched;   //!
  TBranch        *b_HLT_DiPFJet25_NoCaloMatched;   //!
  TBranch        *b_HLT_DiPFJet15_FBEta3_NoCaloMatched;   //!
  TBranch        *b_HLT_DiPFJet25_FBEta3_NoCaloMatched;   //!
  TBranch        *b_HLT_DiPFJetAve15_HFJEC;   //!
  TBranch        *b_HLT_DiPFJetAve25_HFJEC;   //!
  TBranch        *b_HLT_DiPFJetAve35_HFJEC;   //!
  TBranch        *b_HLT_AK8PFJet40;   //!
  TBranch        *b_HLT_AK8PFJet60;   //!
  TBranch        *b_HLT_AK8PFJet80;   //!
  TBranch        *b_HLT_AK8PFJet140;   //!
  TBranch        *b_HLT_AK8PFJet200;   //!
  TBranch        *b_HLT_AK8PFJet260;   //!
  TBranch        *b_HLT_AK8PFJet320;   //!
  TBranch        *b_HLT_AK8PFJet400;   //!
  TBranch        *b_HLT_AK8PFJet450;   //!
  TBranch        *b_HLT_AK8PFJet500;   //!
  TBranch        *b_HLT_PFJet40;   //!
  TBranch        *b_HLT_PFJet60;   //!
  TBranch        *b_HLT_PFJet80;   //!
  TBranch        *b_HLT_PFJet140;   //!
  TBranch        *b_HLT_PFJet200;   //!
  TBranch        *b_HLT_PFJet260;   //!
  TBranch        *b_HLT_PFJet320;   //!
  TBranch        *b_HLT_PFJet400;   //!
  TBranch        *b_HLT_PFJet450;   //!
  TBranch        *b_HLT_PFJet500;   //!
  TBranch        *b_HLT_DiPFJetAve40;   //!
  TBranch        *b_HLT_DiPFJetAve60;   //!
  TBranch        *b_HLT_DiPFJetAve80;   //!
  TBranch        *b_HLT_DiPFJetAve140;   //!
  TBranch        *b_HLT_DiPFJetAve200;   //!
  TBranch        *b_HLT_DiPFJetAve260;   //!
  TBranch        *b_HLT_DiPFJetAve320;   //!
  TBranch        *b_HLT_DiPFJetAve400;   //!
  TBranch        *b_HLT_DiPFJetAve500;   //!
  TBranch        *b_HLT_DiPFJetAve60_HFJEC;   //!
  TBranch        *b_HLT_DiPFJetAve80_HFJEC;   //!
  TBranch        *b_HLT_DiPFJetAve100_HFJEC;   //!
  TBranch        *b_HLT_DiPFJetAve160_HFJEC;   //!
  TBranch        *b_HLT_DiPFJetAve220_HFJEC;   //!
  TBranch        *b_HLT_DiPFJetAve300_HFJEC;   //!
  TBranch        *b_HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu140;   //!
  TBranch        *b_HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu80;   //!
  TBranch        *b_HLT_DiCentralPFJet170;   //!
  TBranch        *b_HLT_SingleCentralPFJet170_CFMax0p1;   //!
  TBranch        *b_HLT_DiCentralPFJet170_CFMax0p1;   //!
  TBranch        *b_HLT_DiCentralPFJet220_CFMax0p3;   //!
  TBranch        *b_HLT_DiCentralPFJet330_CFMax0p5;   //!
  TBranch        *b_HLT_DiCentralPFJet430;   //!
  TBranch        *b_HLT_PFHT125;   //!
  TBranch        *b_HLT_PFHT200;   //!
  TBranch        *b_HLT_PFHT250;   //!
  TBranch        *b_HLT_PFHT300;   //!
  TBranch        *b_HLT_PFHT350;   //!
  TBranch        *b_HLT_PFHT400;   //!
  TBranch        *b_HLT_PFHT475;   //!
  TBranch        *b_HLT_PFHT600;   //!
  TBranch        *b_HLT_PFHT650;   //!
  TBranch        *b_HLT_PFHT800;   //!
  TBranch        *b_HLT_PFHT900;   //!
  TBranch        *b_HLT_PFHT200_PFAlphaT0p51;   //!
  TBranch        *b_HLT_PFHT200_DiPFJetAve90_PFAlphaT0p57;   //!
  TBranch        *b_HLT_PFHT200_DiPFJetAve90_PFAlphaT0p63;   //!
  TBranch        *b_HLT_PFHT250_DiPFJetAve90_PFAlphaT0p55;   //!
  TBranch        *b_HLT_PFHT250_DiPFJetAve90_PFAlphaT0p58;   //!
  TBranch        *b_HLT_PFHT300_DiPFJetAve90_PFAlphaT0p53;   //!
  TBranch        *b_HLT_PFHT300_DiPFJetAve90_PFAlphaT0p54;   //!
  TBranch        *b_HLT_PFHT350_DiPFJetAve90_PFAlphaT0p52;   //!
  TBranch        *b_HLT_PFHT350_DiPFJetAve90_PFAlphaT0p53;   //!
  TBranch        *b_HLT_PFHT400_DiPFJetAve90_PFAlphaT0p51;   //!
  TBranch        *b_HLT_PFHT400_DiPFJetAve90_PFAlphaT0p52;   //!
  TBranch        *b_HLT_MET60_IsoTrk35_Loose;   //!
  TBranch        *b_HLT_MET75_IsoTrk50;   //!
  TBranch        *b_HLT_MET90_IsoTrk50;   //!
  TBranch        *b_HLT_PFMET120_BTagCSV_p067;   //!
  TBranch        *b_HLT_PFMET120_Mu5;   //!
  TBranch        *b_HLT_PFMET170_NotCleaned;   //!
  TBranch        *b_HLT_PFMET170_NoiseCleaned;   //!
  TBranch        *b_HLT_PFMET170_HBHECleaned;   //!
  TBranch        *b_HLT_PFMET170_JetIdCleaned;   //!
  TBranch        *b_HLT_PFMET170_BeamHaloCleaned;   //!
  TBranch        *b_HLT_PFMET170_HBHE_BeamHaloCleaned;   //!
  TBranch        *b_HLT_PFMETTypeOne190_HBHE_BeamHaloCleaned;   //!
  TBranch        *b_HLT_PFMET90_PFMHT90_IDTight;   //!
  TBranch        *b_HLT_PFMET100_PFMHT100_IDTight;   //!
  TBranch        *b_HLT_PFMET100_PFMHT100_IDTight_BeamHaloCleaned;   //!
  TBranch        *b_HLT_PFMET110_PFMHT110_IDTight;   //!
  TBranch        *b_HLT_PFMET120_PFMHT120_IDTight;   //!
  TBranch        *b_HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight_BTagCSV_p067;   //!
  TBranch        *b_HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight;   //!
  TBranch        *b_HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq200;   //!
  TBranch        *b_HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq460;   //!
  TBranch        *b_HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240;   //!
  TBranch        *b_HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq500;   //!
  TBranch        *b_HLT_QuadPFJet_VBF;   //!
  TBranch        *b_HLT_L1_TripleJet_VBF;   //!
  TBranch        *b_HLT_QuadJet45_TripleBTagCSV_p087;   //!
  TBranch        *b_HLT_QuadJet45_DoubleBTagCSV_p087;   //!
  TBranch        *b_HLT_DoubleJet90_Double30_TripleBTagCSV_p087;   //!
  TBranch        *b_HLT_DoubleJet90_Double30_DoubleBTagCSV_p087;   //!
  TBranch        *b_HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160;   //!
  TBranch        *b_HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6;   //!
  TBranch        *b_HLT_DoubleJetsC112_DoubleBTagCSV_p026_DoublePFJetsC172;   //!
  TBranch        *b_HLT_DoubleJetsC112_DoubleBTagCSV_p014_DoublePFJetsC112MaxDeta1p6;   //!
  TBranch        *b_HLT_DoubleJetsC100_SingleBTagCSV_p026;   //!
  TBranch        *b_HLT_DoubleJetsC100_SingleBTagCSV_p014;   //!
  TBranch        *b_HLT_DoubleJetsC100_SingleBTagCSV_p026_SinglePFJetC350;   //!
  TBranch        *b_HLT_DoubleJetsC100_SingleBTagCSV_p014_SinglePFJetC350;   //!
  TBranch        *b_HLT_Photon135_PFMET100;   //!
  TBranch        *b_HLT_Photon20_CaloIdVL_IsoL;   //!
  TBranch        *b_HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_PFMET40;   //!
  TBranch        *b_HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_VBF;   //!
  TBranch        *b_HLT_Photon250_NoHE;   //!
  TBranch        *b_HLT_Photon300_NoHE;   //!
  TBranch        *b_HLT_Photon26_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon16_AND_HE10_R9Id65_Eta2_Mass60;   //!
  TBranch        *b_HLT_Photon36_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon22_AND_HE10_R9Id65_Eta2_Mass15;   //!
  TBranch        *b_HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_PFMET40;   //!
  TBranch        *b_HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_VBF;   //!
  TBranch        *b_HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_PFMET40;   //!
  TBranch        *b_HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_VBF;   //!
  TBranch        *b_HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_PFMET40;   //!
  TBranch        *b_HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF;   //!
  TBranch        *b_HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_PFMET40;   //!
  TBranch        *b_HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_VBF;   //!
  TBranch        *b_HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_PFMET40;   //!
  TBranch        *b_HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_VBF;   //!
  TBranch        *b_HLT_Mu8_TrkIsoVVL;   //!
  TBranch        *b_HLT_Mu17_TrkIsoVVL;   //!
  TBranch        *b_HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30;   //!
  TBranch        *b_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30;   //!
  TBranch        *b_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30;   //!
  TBranch        *b_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30;   //!
  TBranch        *b_HLT_BTagMu_DiJet20_Mu5;   //!
  TBranch        *b_HLT_BTagMu_DiJet40_Mu5;   //!
  TBranch        *b_HLT_BTagMu_DiJet70_Mu5;   //!
  TBranch        *b_HLT_BTagMu_DiJet110_Mu5;   //!
  TBranch        *b_HLT_BTagMu_DiJet170_Mu5;   //!
  TBranch        *b_HLT_BTagMu_Jet300_Mu5;   //!
  TBranch        *b_HLT_BTagMu_AK8Jet300_Mu5;   //!
  TBranch        *b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;   //!
  TBranch        *b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_L1JetTauSeeded;   //!
  TBranch        *b_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;   //!
  TBranch        *b_HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL;   //!
  TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL;   //!
  TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;   //!
  TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;   //!
  TBranch        *b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;   //!
  TBranch        *b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;   //!
  TBranch        *b_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;   //!
  TBranch        *b_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL;   //!
  TBranch        *b_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ;   //!
  TBranch        *b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;   //!
  TBranch        *b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;   //!
  TBranch        *b_HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL;   //!
  TBranch        *b_HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL;   //!
  TBranch        *b_HLT_Mu37_Ele27_CaloIdL_GsfTrkIdVL;   //!
  TBranch        *b_HLT_Mu27_Ele37_CaloIdL_GsfTrkIdVL;   //!
  TBranch        *b_HLT_Mu8_DiEle12_CaloIdL_TrackIdL;   //!
  TBranch        *b_HLT_Mu12_Photon25_CaloIdL;   //!
  TBranch        *b_HLT_Mu12_Photon25_CaloIdL_L1ISO;   //!
  TBranch        *b_HLT_Mu12_Photon25_CaloIdL_L1OR;   //!
  TBranch        *b_HLT_Mu17_Photon22_CaloIdL_L1ISO;   //!
  TBranch        *b_HLT_Mu17_Photon30_CaloIdL_L1ISO;   //!
  TBranch        *b_HLT_Mu17_Photon35_CaloIdL_L1ISO;   //!
  TBranch        *b_HLT_DiMu9_Ele9_CaloIdL_TrackIdL;   //!
  TBranch        *b_HLT_TripleMu_5_3_3;   //!
  TBranch        *b_HLT_TripleMu_12_10_5;   //!
  TBranch        *b_HLT_Mu3er_PFHT140_PFMET125;   //!
  TBranch        *b_HLT_Mu6_PFHT200_PFMET80_BTagCSV_p067;   //!
  TBranch        *b_HLT_Mu6_PFHT200_PFMET100;   //!
  TBranch        *b_HLT_Mu14er_PFMET100;   //!
  TBranch        *b_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL;   //!
  TBranch        *b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;   //!
  TBranch        *b_HLT_Ele12_CaloIdL_TrackIdL_IsoVL;   //!
  TBranch        *b_HLT_Ele17_CaloIdL_GsfTrkIdVL;   //!
  TBranch        *b_HLT_Ele17_CaloIdL_TrackIdL_IsoVL;   //!
  TBranch        *b_HLT_Ele23_CaloIdL_TrackIdL_IsoVL;   //!
  TBranch        *b_HLT_PFHT650_WideJetMJJ900DEtaJJ1p5;   //!
  TBranch        *b_HLT_PFHT650_WideJetMJJ950DEtaJJ1p5;   //!
  TBranch        *b_HLT_Photon22;   //!
  TBranch        *b_HLT_Photon30;   //!
  TBranch        *b_HLT_Photon36;   //!
  TBranch        *b_HLT_Photon50;   //!
  TBranch        *b_HLT_Photon75;   //!
  TBranch        *b_HLT_Photon90;   //!
  TBranch        *b_HLT_Photon120;   //!
  TBranch        *b_HLT_Photon175;   //!
  TBranch        *b_HLT_Photon165_HE10;   //!
  TBranch        *b_HLT_Photon22_R9Id90_HE10_IsoM;   //!
  TBranch        *b_HLT_Photon30_R9Id90_HE10_IsoM;   //!
  TBranch        *b_HLT_Photon36_R9Id90_HE10_IsoM;   //!
  TBranch        *b_HLT_Photon50_R9Id90_HE10_IsoM;   //!
  TBranch        *b_HLT_Photon75_R9Id90_HE10_IsoM;   //!
  TBranch        *b_HLT_Photon90_R9Id90_HE10_IsoM;   //!
  TBranch        *b_HLT_Photon120_R9Id90_HE10_IsoM;   //!
  TBranch        *b_HLT_Photon165_R9Id90_HE10_IsoM;   //!
  TBranch        *b_HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90;   //!
  TBranch        *b_HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70;   //!
  TBranch        *b_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55;   //!
  TBranch        *b_HLT_Diphoton30_18_Solid_R9Id_AND_IsoCaloId_AND_HE_R9Id_Mass55;   //!
  TBranch        *b_HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55;   //!
  TBranch        *b_HLT_Dimuon0_Jpsi_Muon;   //!
  TBranch        *b_HLT_Dimuon0_Upsilon_Muon;   //!
  TBranch        *b_HLT_QuadMuon0_Dimuon0_Jpsi;   //!
  TBranch        *b_HLT_QuadMuon0_Dimuon0_Upsilon;   //!
  TBranch        *b_HLT_Rsq0p25_Calo;   //!
  TBranch        *b_HLT_RsqMR240_Rsq0p09_MR200_4jet_Calo;   //!
  TBranch        *b_HLT_RsqMR240_Rsq0p09_MR200_Calo;   //!
  TBranch        *b_HLT_Rsq0p25;   //!
  TBranch        *b_HLT_Rsq0p30;   //!
  TBranch        *b_HLT_RsqMR240_Rsq0p09_MR200;   //!
  TBranch        *b_HLT_RsqMR240_Rsq0p09_MR200_4jet;   //!
  TBranch        *b_HLT_RsqMR270_Rsq0p09_MR200;   //!
  TBranch        *b_HLT_RsqMR270_Rsq0p09_MR200_4jet;   //!
  TBranch        *b_HLT_Rsq0p02_MR300_TriPFJet80_60_40_BTagCSV_p063_p20_Mbb60_200;   //!
  TBranch        *b_HLT_Rsq0p02_MR400_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200;   //!
  TBranch        *b_HLT_Rsq0p02_MR450_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200;   //!
  TBranch        *b_HLT_Rsq0p02_MR500_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200;   //!
  TBranch        *b_HLT_Rsq0p02_MR550_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200;   //!
  TBranch        *b_HLT_HT200_DisplacedDijet40_DisplacedTrack;   //!
  TBranch        *b_HLT_HT250_DisplacedDijet40_DisplacedTrack;   //!
  TBranch        *b_HLT_HT350_DisplacedDijet40_DisplacedTrack;   //!
  TBranch        *b_HLT_HT350_DisplacedDijet80_DisplacedTrack;   //!
  TBranch        *b_HLT_HT350_DisplacedDijet80_Tight_DisplacedTrack;   //!
  TBranch        *b_HLT_HT350_DisplacedDijet40_Inclusive;   //!
  TBranch        *b_HLT_HT400_DisplacedDijet40_Inclusive;   //!
  TBranch        *b_HLT_HT500_DisplacedDijet40_Inclusive;   //!
  TBranch        *b_HLT_HT550_DisplacedDijet40_Inclusive;   //!
  TBranch        *b_HLT_HT550_DisplacedDijet80_Inclusive;   //!
  TBranch        *b_HLT_HT650_DisplacedDijet80_Inclusive;   //!
  TBranch        *b_HLT_HT750_DisplacedDijet80_Inclusive;   //!
  TBranch        *b_HLT_VBF_DisplacedJet40_DisplacedTrack;   //!
  TBranch        *b_HLT_VBF_DisplacedJet40_DisplacedTrack_2TrackIP2DSig5;   //!
  TBranch        *b_HLT_VBF_DisplacedJet40_TightID_DisplacedTrack;   //!
  TBranch        *b_HLT_VBF_DisplacedJet40_Hadronic;   //!
  TBranch        *b_HLT_VBF_DisplacedJet40_Hadronic_2PromptTrack;   //!
  TBranch        *b_HLT_VBF_DisplacedJet40_TightID_Hadronic;   //!
  TBranch        *b_HLT_VBF_DisplacedJet40_VTightID_Hadronic;   //!
  TBranch        *b_HLT_VBF_DisplacedJet40_VVTightID_Hadronic;   //!
  TBranch        *b_HLT_VBF_DisplacedJet40_VTightID_DisplacedTrack;   //!
  TBranch        *b_HLT_VBF_DisplacedJet40_VVTightID_DisplacedTrack;   //!
  TBranch        *b_HLT_PFMETNoMu90_PFMHTNoMu90_IDTight;   //!
  TBranch        *b_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight;   //!
  TBranch        *b_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight;   //!
  TBranch        *b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight;   //!
  TBranch        *b_HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_IDTight;   //!
  TBranch        *b_HLT_MonoCentralPFJet80_PFMETNoMu100_PFMHTNoMu100_IDTight;   //!
  TBranch        *b_HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight;   //!
  TBranch        *b_HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight;   //!
  TBranch        *b_HLT_Ele27_eta2p1_WPLoose_Gsf_HT200;   //!
  TBranch        *b_HLT_Photon90_CaloIdL_PFHT500;   //!
  TBranch        *b_HLT_DoubleMu8_Mass8_PFHT250;   //!
  TBranch        *b_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT250;   //!
  TBranch        *b_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT250;   //!
  TBranch        *b_HLT_DoubleMu8_Mass8_PFHT300;   //!
  TBranch        *b_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300;   //!
  TBranch        *b_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300;   //!
  TBranch        *b_HLT_Mu10_CentralPFJet30_BTagCSV_p13;   //!
  TBranch        *b_HLT_DoubleMu3_PFMET50;   //!
  TBranch        *b_HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13;   //!
  TBranch        *b_HLT_Ele15_IsoVVVL_BTagCSV_p067_PFHT400;   //!
  TBranch        *b_HLT_Ele15_IsoVVVL_PFHT350_PFMET50;   //!
  TBranch        *b_HLT_Ele15_IsoVVVL_PFHT600;   //!
  TBranch        *b_HLT_Ele15_IsoVVVL_PFHT350;   //!
  TBranch        *b_HLT_Ele15_IsoVVVL_PFHT400_PFMET50;   //!
  TBranch        *b_HLT_Ele15_IsoVVVL_PFHT400;   //!
  TBranch        *b_HLT_Ele50_IsoVVVL_PFHT400;   //!
  TBranch        *b_HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60;   //!
  TBranch        *b_HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60;   //!
  TBranch        *b_HLT_Mu15_IsoVVVL_BTagCSV_p067_PFHT400;   //!
  TBranch        *b_HLT_Mu15_IsoVVVL_PFHT350_PFMET50;   //!
  TBranch        *b_HLT_Mu15_IsoVVVL_PFHT600;   //!
  TBranch        *b_HLT_Mu15_IsoVVVL_PFHT350;   //!
  TBranch        *b_HLT_Mu15_IsoVVVL_PFHT400_PFMET50;   //!
  TBranch        *b_HLT_Mu15_IsoVVVL_PFHT400;   //!
  TBranch        *b_HLT_Mu50_IsoVVVL_PFHT400;   //!
  TBranch        *b_HLT_Dimuon16_Jpsi;   //!
  TBranch        *b_HLT_Dimuon10_Jpsi_Barrel;   //!
  TBranch        *b_HLT_Dimuon8_PsiPrime_Barrel;   //!
  TBranch        *b_HLT_Dimuon8_Upsilon_Barrel;   //!
  TBranch        *b_HLT_Dimuon0_Phi_Barrel;   //!
  TBranch        *b_HLT_Mu16_TkMu0_dEta18_Onia;   //!
  TBranch        *b_HLT_Mu16_TkMu0_dEta18_Phi;   //!
  TBranch        *b_HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx;   //!
  TBranch        *b_HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx;   //!
  TBranch        *b_HLT_Mu8;   //!
  TBranch        *b_HLT_Mu17;   //!
  TBranch        *b_HLT_Mu3_PFJet40;   //!
  TBranch        *b_HLT_Ele8_CaloIdM_TrackIdM_PFJet30;   //!
  TBranch        *b_HLT_Ele12_CaloIdM_TrackIdM_PFJet30;   //!
  TBranch        *b_HLT_Ele17_CaloIdM_TrackIdM_PFJet30;   //!
  TBranch        *b_HLT_Ele23_CaloIdM_TrackIdM_PFJet30;   //!
  TBranch        *b_HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet140;   //!
  TBranch        *b_HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165;   //!
  TBranch        *b_HLT_PFHT400_SixJet30_DoubleBTagCSV_p056;   //!
  TBranch        *b_HLT_PFHT450_SixJet40_BTagCSV_p056;   //!
  TBranch        *b_HLT_PFHT400_SixJet30;   //!
  TBranch        *b_HLT_PFHT450_SixJet40;   //!
  TBranch        *b_HLT_Ele115_CaloIdVT_GsfTrkIdT;   //!
  TBranch        *b_HLT_Mu55;   //!
  TBranch        *b_HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15;   //!
  TBranch        *b_HLT_Photon90_CaloIdL_PFHT600;   //!
  TBranch        *b_HLT_PixelTracks_Multiplicity60ForEndOfFill;   //!
  TBranch        *b_HLT_PixelTracks_Multiplicity85ForEndOfFill;   //!
  TBranch        *b_HLT_PixelTracks_Multiplicity110ForEndOfFill;   //!
  TBranch        *b_HLT_PixelTracks_Multiplicity135ForEndOfFill;   //!
  TBranch        *b_HLT_PixelTracks_Multiplicity160ForEndOfFill;   //!
  TBranch        *b_HLT_FullTracks_Multiplicity80;   //!
  TBranch        *b_HLT_FullTracks_Multiplicity100;   //!
  TBranch        *b_HLT_FullTracks_Multiplicity130;   //!
  TBranch        *b_HLT_FullTracks_Multiplicity150;   //!
  TBranch        *b_HLT_ECALHT800;   //!
  TBranch        *b_HLT_DiSC30_18_EIso_AND_HE_Mass70;   //!
  TBranch        *b_HLT_Photon125;   //!
  TBranch        *b_HLT_MET100;   //!
  TBranch        *b_HLT_MET150;   //!
  TBranch        *b_HLT_MET200;   //!
  TBranch        *b_HLT_Ele27_HighEta_Ele20_Mass55;   //!
  TBranch        *b_HLT_L1FatEvents;   //!
  TBranch        *b_HLT_Physics;   //!
  TBranch        *b_HLT_L1FatEvents_part0;   //!
  TBranch        *b_HLT_L1FatEvents_part1;   //!
  TBranch        *b_HLT_L1FatEvents_part2;   //!
  TBranch        *b_HLT_L1FatEvents_part3;   //!
  TBranch        *b_HLT_Random;   //!
  TBranch        *b_HLT_ZeroBias;   //!
  TBranch        *b_HLT_AK4CaloJet30;   //!
  TBranch        *b_HLT_AK4CaloJet40;   //!
  TBranch        *b_HLT_AK4CaloJet50;   //!
  TBranch        *b_HLT_AK4CaloJet80;   //!
  TBranch        *b_HLT_AK4CaloJet100;   //!
  TBranch        *b_HLT_AK4PFJet30;   //!
  TBranch        *b_HLT_AK4PFJet50;   //!
  TBranch        *b_HLT_AK4PFJet80;   //!
  TBranch        *b_HLT_AK4PFJet100;   //!
  TBranch        *b_HLT_HISinglePhoton10;   //!
  TBranch        *b_HLT_HISinglePhoton15;   //!
  TBranch        *b_HLT_HISinglePhoton20;   //!
  TBranch        *b_HLT_HISinglePhoton40;   //!
  TBranch        *b_HLT_HISinglePhoton60;   //!
  TBranch        *b_HLT_EcalCalibration;   //!
  TBranch        *b_HLT_HcalCalibration;   //!
  TBranch        *b_HLT_GlobalRunHPDNoise;   //!
  TBranch        *b_HLT_L1BptxMinus;   //!
  TBranch        *b_HLT_L1BptxPlus;   //!
  TBranch        *b_HLT_L1NotBptxOR;   //!
  TBranch        *b_HLT_L1BeamGasMinus;   //!
  TBranch        *b_HLT_L1BeamGasPlus;   //!
  TBranch        *b_HLT_L1BptxXOR;   //!
  TBranch        *b_HLT_L1MinimumBiasHF_OR;   //!
  TBranch        *b_HLT_L1MinimumBiasHF_AND;   //!
  TBranch        *b_HLT_HcalNZS;   //!
  TBranch        *b_HLT_HcalPhiSym;   //!
  TBranch        *b_HLT_HcalIsolatedbunch;   //!
  TBranch        *b_HLT_ZeroBias_FirstCollisionAfterAbortGap;   //!
  TBranch        *b_HLT_ZeroBias_FirstCollisionAfterAbortGap_copy;   //!
  TBranch        *b_HLT_ZeroBias_FirstCollisionAfterAbortGap_TCDS;   //!
  TBranch        *b_HLT_ZeroBias_IsolatedBunches;   //!
  TBranch        *b_HLT_ZeroBias_FirstCollisionInTrain;   //!
  TBranch        *b_HLT_ZeroBias_FirstBXAfterTrain;   //!
  TBranch        *b_HLT_Photon500;   //!
  TBranch        *b_HLT_Photon600;   //!
  TBranch        *b_HLT_Mu300;   //!
  TBranch        *b_HLT_Mu350;   //!
  TBranch        *b_HLT_MET250;   //!
  TBranch        *b_HLT_MET300;   //!
  TBranch        *b_HLT_MET600;   //!
  TBranch        *b_HLT_MET700;   //!
  TBranch        *b_HLT_PFMET300;   //!
  TBranch        *b_HLT_PFMET400;   //!
  TBranch        *b_HLT_PFMET500;   //!
  TBranch        *b_HLT_PFMET600;   //!
  TBranch        *b_HLT_Ele250_CaloIdVT_GsfTrkIdT;   //!
  TBranch        *b_HLT_Ele300_CaloIdVT_GsfTrkIdT;   //!
  TBranch        *b_HLT_HT2000;   //!
  TBranch        *b_HLT_HT2500;   //!
  TBranch        *b_HLT_IsoTrackHE;   //!
  TBranch        *b_HLT_IsoTrackHB;   //!
  TBranch        *b_HLTriggerFinalPath;   //!
  TBranch        *b_Flag_HBHENoiseFilter;   //!
  TBranch        *b_Flag_HBHENoiseIsoFilter;   //!
  TBranch        *b_Flag_CSCTightHaloFilter;   //!
  TBranch        *b_Flag_CSCTightHaloTrkMuUnvetoFilter;   //!
  TBranch        *b_Flag_CSCTightHalo2015Filter;   //!
  TBranch        *b_Flag_globalTightHalo2016Filter;   //!
  TBranch        *b_Flag_globalSuperTightHalo2016Filter;   //!
  TBranch        *b_Flag_HcalStripHaloFilter;   //!
  TBranch        *b_Flag_hcalLaserEventFilter;   //!
  TBranch        *b_Flag_EcalDeadCellTriggerPrimitiveFilter;   //!
  TBranch        *b_Flag_EcalDeadCellBoundaryEnergyFilter;   //!
  TBranch        *b_Flag_goodVertices;   //!
  TBranch        *b_Flag_eeBadScFilter;   //!
  TBranch        *b_Flag_ecalLaserCorrFilter;   //!
  TBranch        *b_Flag_trkPOGFilters;   //!
  TBranch        *b_Flag_chargedHadronTrackResolutionFilter;   //!
  TBranch        *b_Flag_muonBadTrackFilter;   //!
  TBranch        *b_Flag_trkPOG_manystripclus53X;   //!
  TBranch        *b_Flag_trkPOG_toomanystripclus53X;   //!
  TBranch        *b_Flag_trkPOG_logErrorTooManyClusters;   //!
  TBranch        *b_Flag_METFilters;   //!
      
      
  nanoAnalysis(TTree *tree=0, Bool_t flag = false, Bool_t dl = false, Bool_t sle = false, Bool_t slm = false);
  virtual ~nanoAnalysis();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef nanoAnalysis_cxx
nanoAnalysis::nanoAnalysis(TTree *tree, Bool_t flag, Bool_t dl, Bool_t sle, Bool_t slm) : fChain(0), m_isMC(flag), m_isDL(dl), m_isSL_e(sle), m_isSL_m(slm)
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/home/nanoAOD/run2_2016v3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/180117_180123/0000/nanoAOD_111.root");
    if (!f || !f->IsOpen()) {
      f = TFile::Open("/home/nanoAOD/run2_2016v3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/180117_180123/0000/nanoAOD_111.root");
    }
    f->GetObject("Events",tree);

  }
    

  Init(tree);
}

nanoAnalysis::~nanoAnalysis()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
  m_output->Write();
  m_output->Close();
}

Int_t nanoAnalysis::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t nanoAnalysis::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void nanoAnalysis::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("run", &run, &b_run);
  fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
  fChain->SetBranchAddress("event", &event, &b_event);
  fChain->SetBranchAddress("CaloMET_phi", &CaloMET_phi, &b_CaloMET_phi);
  fChain->SetBranchAddress("CaloMET_pt", &CaloMET_pt, &b_CaloMET_pt);
  fChain->SetBranchAddress("CaloMET_sumEt", &CaloMET_sumEt, &b_CaloMET_sumEt);
  fChain->SetBranchAddress("ncmeson", &ncmeson, &b_ncmeson);
  fChain->SetBranchAddress("cmeson_dca", cmeson_dca, &b_cmeson_dca);
  fChain->SetBranchAddress("cmeson_angleXY", cmeson_angleXY, &b_cmeson_angleXY);
  fChain->SetBranchAddress("cmeson_angleXYZ", cmeson_angleXYZ, &b_cmeson_angleXYZ);
  fChain->SetBranchAddress("cmeson_trk_normalizedChi2", cmeson_trk_normalizedChi2, &b_cmeson_trk_normalizedChi2);
  fChain->SetBranchAddress("cmeson_trk_nHits", cmeson_trk_nHits, &b_cmeson_trk_nHits);
  fChain->SetBranchAddress("cmeson_trk_pt", cmeson_trk_pt, &b_cmeson_trk_pt);
  fChain->SetBranchAddress("cmeson_trk_ipsigXY", cmeson_trk_ipsigXY, &b_cmeson_trk_ipsigXY);
  fChain->SetBranchAddress("cmeson_trk_ipsigZ", cmeson_trk_ipsigZ, &b_cmeson_trk_ipsigZ);
  fChain->SetBranchAddress("cmeson_lxy", cmeson_lxy, &b_cmeson_lxy);
  fChain->SetBranchAddress("cmeson_lxySig", cmeson_lxySig, &b_cmeson_lxySig);
  fChain->SetBranchAddress("cmeson_l3D", cmeson_l3D, &b_cmeson_l3D);
  fChain->SetBranchAddress("cmeson_l3DSig", cmeson_l3DSig, &b_cmeson_l3DSig);
  fChain->SetBranchAddress("cmeson_jetDR", cmeson_jetDR, &b_cmeson_jetDR);
  fChain->SetBranchAddress("cmeson_legDR", cmeson_legDR, &b_cmeson_legDR);
  fChain->SetBranchAddress("cmeson_diffMass", cmeson_diffMass, &b_cmeson_diffMass);
  fChain->SetBranchAddress("cmeson_nJet", cmeson_nJet, &b_cmeson_nJet);
  fChain->SetBranchAddress("cmeson_mcMatch", cmeson_mcMatch, &b_cmeson_mcMatch);
  fChain->SetBranchAddress("nElectron", &nElectron, &b_nElectron);
  fChain->SetBranchAddress("Electron_deltaEtaSC", Electron_deltaEtaSC, &b_Electron_deltaEtaSC);
  fChain->SetBranchAddress("Electron_dr03EcalRecHitSumEt", Electron_dr03EcalRecHitSumEt, &b_Electron_dr03EcalRecHitSumEt);
  fChain->SetBranchAddress("Electron_dr03HcalDepth1TowerSumEt", Electron_dr03HcalDepth1TowerSumEt, &b_Electron_dr03HcalDepth1TowerSumEt);
  fChain->SetBranchAddress("Electron_dr03TkSumPt", Electron_dr03TkSumPt, &b_Electron_dr03TkSumPt);
  fChain->SetBranchAddress("Electron_dxy", Electron_dxy, &b_Electron_dxy);
  fChain->SetBranchAddress("Electron_dxyErr", Electron_dxyErr, &b_Electron_dxyErr);
  fChain->SetBranchAddress("Electron_dz", Electron_dz, &b_Electron_dz);
  fChain->SetBranchAddress("Electron_dzErr", Electron_dzErr, &b_Electron_dzErr);
  fChain->SetBranchAddress("Electron_eCorr", Electron_eCorr, &b_Electron_eCorr);
  fChain->SetBranchAddress("Electron_eInvMinusPInv", Electron_eInvMinusPInv, &b_Electron_eInvMinusPInv);
  fChain->SetBranchAddress("Electron_energyErr", Electron_energyErr, &b_Electron_energyErr);
  fChain->SetBranchAddress("Electron_eta", Electron_eta, &b_Electron_eta);
  fChain->SetBranchAddress("Electron_hoe", Electron_hoe, &b_Electron_hoe);
  fChain->SetBranchAddress("Electron_ip3d", Electron_ip3d, &b_Electron_ip3d);
  fChain->SetBranchAddress("Electron_mass", Electron_mass, &b_Electron_mass);
  fChain->SetBranchAddress("Electron_miniPFRelIso_all", Electron_miniPFRelIso_all, &b_Electron_miniPFRelIso_all);
  fChain->SetBranchAddress("Electron_miniPFRelIso_chg", Electron_miniPFRelIso_chg, &b_Electron_miniPFRelIso_chg);
  fChain->SetBranchAddress("Electron_mvaSpring16GP", Electron_mvaSpring16GP, &b_Electron_mvaSpring16GP);
  fChain->SetBranchAddress("Electron_mvaSpring16HZZ", Electron_mvaSpring16HZZ, &b_Electron_mvaSpring16HZZ);
  fChain->SetBranchAddress("Electron_pfRelIso03_all", Electron_pfRelIso03_all, &b_Electron_pfRelIso03_all);
  fChain->SetBranchAddress("Electron_pfRelIso03_chg", Electron_pfRelIso03_chg, &b_Electron_pfRelIso03_chg);
  fChain->SetBranchAddress("Electron_phi", Electron_phi, &b_Electron_phi);
  fChain->SetBranchAddress("Electron_pt", Electron_pt, &b_Electron_pt);
  fChain->SetBranchAddress("Electron_r9", Electron_r9, &b_Electron_r9);
  fChain->SetBranchAddress("Electron_sieie", Electron_sieie, &b_Electron_sieie);
  fChain->SetBranchAddress("Electron_sip3d", Electron_sip3d, &b_Electron_sip3d);
  fChain->SetBranchAddress("Electron_mvaTTH", Electron_mvaTTH, &b_Electron_mvaTTH);
  fChain->SetBranchAddress("Electron_charge", Electron_charge, &b_Electron_charge);
  fChain->SetBranchAddress("Electron_cutBased", Electron_cutBased, &b_Electron_cutBased);
  fChain->SetBranchAddress("Electron_cutBased_HLTPreSel", Electron_cutBased_HLTPreSel, &b_Electron_cutBased_HLTPreSel);
  fChain->SetBranchAddress("Electron_jetIdx", Electron_jetIdx, &b_Electron_jetIdx);
  fChain->SetBranchAddress("Electron_pdgId", Electron_pdgId, &b_Electron_pdgId);
  fChain->SetBranchAddress("Electron_photonIdx", Electron_photonIdx, &b_Electron_photonIdx);
  fChain->SetBranchAddress("Electron_tightCharge", Electron_tightCharge, &b_Electron_tightCharge);
  fChain->SetBranchAddress("Electron_vidNestedWPBitmap", Electron_vidNestedWPBitmap, &b_Electron_vidNestedWPBitmap);
  fChain->SetBranchAddress("Electron_convVeto", Electron_convVeto, &b_Electron_convVeto);
  fChain->SetBranchAddress("Electron_cutBased_HEEP", Electron_cutBased_HEEP, &b_Electron_cutBased_HEEP);
  fChain->SetBranchAddress("Electron_isPFcand", Electron_isPFcand, &b_Electron_isPFcand);
  fChain->SetBranchAddress("Electron_lostHits", Electron_lostHits, &b_Electron_lostHits);
  fChain->SetBranchAddress("Electron_mvaSpring16GP_WP80", Electron_mvaSpring16GP_WP80, &b_Electron_mvaSpring16GP_WP80);
  fChain->SetBranchAddress("Electron_mvaSpring16GP_WP90", Electron_mvaSpring16GP_WP90, &b_Electron_mvaSpring16GP_WP90);
  fChain->SetBranchAddress("Electron_mvaSpring16HZZ_WPL", Electron_mvaSpring16HZZ_WPL, &b_Electron_mvaSpring16HZZ_WPL);
  fChain->SetBranchAddress("nFatJet", &nFatJet, &b_nFatJet);
  fChain->SetBranchAddress("FatJet_area", FatJet_area, &b_FatJet_area);
  fChain->SetBranchAddress("FatJet_btagCMVA", FatJet_btagCMVA, &b_FatJet_btagCMVA);
  fChain->SetBranchAddress("FatJet_btagCSVV2", FatJet_btagCSVV2, &b_FatJet_btagCSVV2);
  fChain->SetBranchAddress("FatJet_btagDeepB", FatJet_btagDeepB, &b_FatJet_btagDeepB);
  fChain->SetBranchAddress("FatJet_btagHbb", FatJet_btagHbb, &b_FatJet_btagHbb);
  fChain->SetBranchAddress("FatJet_eta", FatJet_eta, &b_FatJet_eta);
  fChain->SetBranchAddress("FatJet_mass", FatJet_mass, &b_FatJet_mass);
  fChain->SetBranchAddress("FatJet_msoftdrop", FatJet_msoftdrop, &b_FatJet_msoftdrop);
  fChain->SetBranchAddress("FatJet_n2b1", FatJet_n2b1, &b_FatJet_n2b1);
  fChain->SetBranchAddress("FatJet_n3b1", FatJet_n3b1, &b_FatJet_n3b1);
  fChain->SetBranchAddress("FatJet_phi", FatJet_phi, &b_FatJet_phi);
  fChain->SetBranchAddress("FatJet_pt", FatJet_pt, &b_FatJet_pt);
  fChain->SetBranchAddress("FatJet_tau1", FatJet_tau1, &b_FatJet_tau1);
  fChain->SetBranchAddress("FatJet_tau2", FatJet_tau2, &b_FatJet_tau2);
  fChain->SetBranchAddress("FatJet_tau3", FatJet_tau3, &b_FatJet_tau3);
  fChain->SetBranchAddress("FatJet_tau4", FatJet_tau4, &b_FatJet_tau4);
  fChain->SetBranchAddress("FatJet_subJetIdx1", FatJet_subJetIdx1, &b_FatJet_subJetIdx1);
  fChain->SetBranchAddress("FatJet_subJetIdx2", FatJet_subJetIdx2, &b_FatJet_subJetIdx2);
  fChain->SetBranchAddress("nGenJetAK8", &nGenJetAK8, &b_nGenJetAK8);
  fChain->SetBranchAddress("GenJetAK8_eta", GenJetAK8_eta, &b_GenJetAK8_eta);
  fChain->SetBranchAddress("GenJetAK8_mass", GenJetAK8_mass, &b_GenJetAK8_mass);
  fChain->SetBranchAddress("GenJetAK8_phi", GenJetAK8_phi, &b_GenJetAK8_phi);
  fChain->SetBranchAddress("GenJetAK8_pt", GenJetAK8_pt, &b_GenJetAK8_pt);
  fChain->SetBranchAddress("nGenJet", &nGenJet, &b_nGenJet);
  fChain->SetBranchAddress("GenJet_eta", GenJet_eta, &b_GenJet_eta);
  fChain->SetBranchAddress("GenJet_mass", GenJet_mass, &b_GenJet_mass);
  fChain->SetBranchAddress("GenJet_phi", GenJet_phi, &b_GenJet_phi);
  fChain->SetBranchAddress("GenJet_pt", GenJet_pt, &b_GenJet_pt);
  fChain->SetBranchAddress("nGenPart", &nGenPart, &b_nGenPart);
  fChain->SetBranchAddress("GenPart_eta", GenPart_eta, &b_GenPart_eta);
  fChain->SetBranchAddress("GenPart_mass", GenPart_mass, &b_GenPart_mass);
  fChain->SetBranchAddress("GenPart_phi", GenPart_phi, &b_GenPart_phi);
  fChain->SetBranchAddress("GenPart_pt", GenPart_pt, &b_GenPart_pt);
  fChain->SetBranchAddress("GenPart_genPartIdxMother", GenPart_genPartIdxMother, &b_GenPart_genPartIdxMother);
  fChain->SetBranchAddress("GenPart_pdgId", GenPart_pdgId, &b_GenPart_pdgId);
  fChain->SetBranchAddress("GenPart_status", GenPart_status, &b_GenPart_status);
  fChain->SetBranchAddress("GenPart_statusFlags", GenPart_statusFlags, &b_GenPart_statusFlags);
  fChain->SetBranchAddress("Generator_scalePDF", &Generator_scalePDF, &b_Generator_scalePDF);
  fChain->SetBranchAddress("Generator_x1", &Generator_x1, &b_Generator_x1);
  fChain->SetBranchAddress("Generator_x2", &Generator_x2, &b_Generator_x2);
  fChain->SetBranchAddress("Generator_xpdf1", &Generator_xpdf1, &b_Generator_xpdf1);
  fChain->SetBranchAddress("Generator_xpdf2", &Generator_xpdf2, &b_Generator_xpdf2);
  fChain->SetBranchAddress("Generator_id1", &Generator_id1, &b_Generator_id1);
  fChain->SetBranchAddress("Generator_id2", &Generator_id2, &b_Generator_id2);
  fChain->SetBranchAddress("nGenVisTau", &nGenVisTau, &b_nGenVisTau);
  fChain->SetBranchAddress("GenVisTau_eta", GenVisTau_eta, &b_GenVisTau_eta);
  fChain->SetBranchAddress("GenVisTau_mass", GenVisTau_mass, &b_GenVisTau_mass);
  fChain->SetBranchAddress("GenVisTau_phi", GenVisTau_phi, &b_GenVisTau_phi);
  fChain->SetBranchAddress("GenVisTau_pt", GenVisTau_pt, &b_GenVisTau_pt);
  fChain->SetBranchAddress("GenVisTau_charge", GenVisTau_charge, &b_GenVisTau_charge);
  fChain->SetBranchAddress("GenVisTau_genPartIdxMother", GenVisTau_genPartIdxMother, &b_GenVisTau_genPartIdxMother);
  fChain->SetBranchAddress("GenVisTau_status", GenVisTau_status, &b_GenVisTau_status);
  fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
  fChain->SetBranchAddress("LHEWeight_originalXWGTUP", &LHEWeight_originalXWGTUP, &b_LHEWeight_originalXWGTUP);
  fChain->SetBranchAddress("nLHEPdfWeight", &nLHEPdfWeight, &b_nLHEPdfWeight);
  fChain->SetBranchAddress("LHEPdfWeight", &LHEPdfWeight, &b_LHEPdfWeight);
  fChain->SetBranchAddress("nLHEScaleWeight", &nLHEScaleWeight, &b_nLHEScaleWeight);
  fChain->SetBranchAddress("LHEScaleWeight", LHEScaleWeight, &b_LHEScaleWeight);
  fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
  fChain->SetBranchAddress("Jet_area", Jet_area, &b_Jet_area);
  fChain->SetBranchAddress("Jet_btagCMVA", Jet_btagCMVA, &b_Jet_btagCMVA);
  fChain->SetBranchAddress("Jet_btagCSVV2", Jet_btagCSVV2, &b_Jet_btagCSVV2);
  fChain->SetBranchAddress("Jet_btagDeepB", Jet_btagDeepB, &b_Jet_btagDeepB);
  fChain->SetBranchAddress("Jet_btagDeepC", Jet_btagDeepC, &b_Jet_btagDeepC);
  fChain->SetBranchAddress("Jet_chEmEF", Jet_chEmEF, &b_Jet_chEmEF);
  fChain->SetBranchAddress("Jet_chHEF", Jet_chHEF, &b_Jet_chHEF);
  fChain->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
  fChain->SetBranchAddress("Jet_mass", Jet_mass, &b_Jet_mass);
  fChain->SetBranchAddress("Jet_neEmEF", Jet_neEmEF, &b_Jet_neEmEF);
  fChain->SetBranchAddress("Jet_neHEF", Jet_neHEF, &b_Jet_neHEF);
  fChain->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
  fChain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
  fChain->SetBranchAddress("Jet_qgl", Jet_qgl, &b_Jet_qgl);
  fChain->SetBranchAddress("Jet_rawFactor", Jet_rawFactor, &b_Jet_rawFactor);
  fChain->SetBranchAddress("Jet_bReg", Jet_bReg, &b_Jet_bReg);
  fChain->SetBranchAddress("Jet_electronIdx1", Jet_electronIdx1, &b_Jet_electronIdx1);
  fChain->SetBranchAddress("Jet_electronIdx2", Jet_electronIdx2, &b_Jet_electronIdx2);
  fChain->SetBranchAddress("Jet_jetId", Jet_jetId, &b_Jet_jetId);
  fChain->SetBranchAddress("Jet_muonIdx1", Jet_muonIdx1, &b_Jet_muonIdx1);
  fChain->SetBranchAddress("Jet_muonIdx2", Jet_muonIdx2, &b_Jet_muonIdx2);
  fChain->SetBranchAddress("Jet_nConstituents", Jet_nConstituents, &b_Jet_nConstituents);
  fChain->SetBranchAddress("Jet_nElectrons", Jet_nElectrons, &b_Jet_nElectrons);
  fChain->SetBranchAddress("Jet_nMuons", Jet_nMuons, &b_Jet_nMuons);
  fChain->SetBranchAddress("Jet_puId", Jet_puId, &b_Jet_puId);
  fChain->SetBranchAddress("LHE_HT", &LHE_HT, &b_LHE_HT);
  fChain->SetBranchAddress("LHE_HTIncoming", &LHE_HTIncoming, &b_LHE_HTIncoming);
  fChain->SetBranchAddress("LHE_Vpt", &LHE_Vpt, &b_LHE_Vpt);
  fChain->SetBranchAddress("LHE_Njets", &LHE_Njets, &b_LHE_Njets);
  fChain->SetBranchAddress("LHE_Nb", &LHE_Nb, &b_LHE_Nb);
  fChain->SetBranchAddress("LHE_Nc", &LHE_Nc, &b_LHE_Nc);
  fChain->SetBranchAddress("LHE_Nuds", &LHE_Nuds, &b_LHE_Nuds);
  fChain->SetBranchAddress("LHE_Nglu", &LHE_Nglu, &b_LHE_Nglu);
  fChain->SetBranchAddress("LHE_NpNLO", &LHE_NpNLO, &b_LHE_NpNLO);
  fChain->SetBranchAddress("LHE_NpLO", &LHE_NpLO, &b_LHE_NpLO);
  fChain->SetBranchAddress("GenMET_phi", &GenMET_phi, &b_GenMET_phi);
  fChain->SetBranchAddress("GenMET_pt", &GenMET_pt, &b_GenMET_pt);
  fChain->SetBranchAddress("MET_MetUnclustEnUpDeltaX", &MET_MetUnclustEnUpDeltaX, &b_MET_MetUnclustEnUpDeltaX);
  fChain->SetBranchAddress("MET_MetUnclustEnUpDeltaY", &MET_MetUnclustEnUpDeltaY, &b_MET_MetUnclustEnUpDeltaY);
  fChain->SetBranchAddress("MET_covXX", &MET_covXX, &b_MET_covXX);
  fChain->SetBranchAddress("MET_covXY", &MET_covXY, &b_MET_covXY);
  fChain->SetBranchAddress("MET_covYY", &MET_covYY, &b_MET_covYY);
  fChain->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
  fChain->SetBranchAddress("MET_pt", &MET_pt, &b_MET_pt);
  fChain->SetBranchAddress("MET_significance", &MET_significance, &b_MET_significance);
  fChain->SetBranchAddress("MET_sumEt", &MET_sumEt, &b_MET_sumEt);
  fChain->SetBranchAddress("nMuon", &nMuon, &b_nMuon);
  fChain->SetBranchAddress("Muon_dxy", Muon_dxy, &b_Muon_dxy);
  fChain->SetBranchAddress("Muon_dxyErr", Muon_dxyErr, &b_Muon_dxyErr);
  fChain->SetBranchAddress("Muon_dz", Muon_dz, &b_Muon_dz);
  fChain->SetBranchAddress("Muon_dzErr", Muon_dzErr, &b_Muon_dzErr);
  fChain->SetBranchAddress("Muon_eta", Muon_eta, &b_Muon_eta);
  fChain->SetBranchAddress("Muon_ip3d", Muon_ip3d, &b_Muon_ip3d);
  fChain->SetBranchAddress("Muon_mass", Muon_mass, &b_Muon_mass);
  fChain->SetBranchAddress("Muon_miniPFRelIso_all", Muon_miniPFRelIso_all, &b_Muon_miniPFRelIso_all);
  fChain->SetBranchAddress("Muon_miniPFRelIso_chg", Muon_miniPFRelIso_chg, &b_Muon_miniPFRelIso_chg);
  fChain->SetBranchAddress("Muon_pfRelIso03_all", Muon_pfRelIso03_all, &b_Muon_pfRelIso03_all);
  fChain->SetBranchAddress("Muon_pfRelIso03_chg", Muon_pfRelIso03_chg, &b_Muon_pfRelIso03_chg);
  fChain->SetBranchAddress("Muon_pfRelIso04_all", Muon_pfRelIso04_all, &b_Muon_pfRelIso04_all);
  fChain->SetBranchAddress("Muon_phi", Muon_phi, &b_Muon_phi);
  fChain->SetBranchAddress("Muon_pt", Muon_pt, &b_Muon_pt);
  fChain->SetBranchAddress("Muon_ptErr", Muon_ptErr, &b_Muon_ptErr);
  fChain->SetBranchAddress("Muon_segmentComp", Muon_segmentComp, &b_Muon_segmentComp);
  fChain->SetBranchAddress("Muon_sip3d", Muon_sip3d, &b_Muon_sip3d);
  fChain->SetBranchAddress("Muon_mvaTTH", Muon_mvaTTH, &b_Muon_mvaTTH);
  fChain->SetBranchAddress("Muon_charge", Muon_charge, &b_Muon_charge);
  fChain->SetBranchAddress("Muon_jetIdx", Muon_jetIdx, &b_Muon_jetIdx);
  fChain->SetBranchAddress("Muon_nStations", Muon_nStations, &b_Muon_nStations);
  fChain->SetBranchAddress("Muon_nTrackerLayers", Muon_nTrackerLayers, &b_Muon_nTrackerLayers);
  fChain->SetBranchAddress("Muon_pdgId", Muon_pdgId, &b_Muon_pdgId);
  fChain->SetBranchAddress("Muon_tightCharge", Muon_tightCharge, &b_Muon_tightCharge);
  fChain->SetBranchAddress("Muon_globalMu", Muon_globalMu, &b_Muon_globalMu);
  fChain->SetBranchAddress("Muon_highPtId", Muon_highPtId, &b_Muon_highPtId);
  fChain->SetBranchAddress("Muon_isPFcand", Muon_isPFcand, &b_Muon_isPFcand);
  fChain->SetBranchAddress("Muon_mediumId", Muon_mediumId, &b_Muon_mediumId);
  fChain->SetBranchAddress("Muon_softId", Muon_softId, &b_Muon_softId);
  fChain->SetBranchAddress("Muon_tightId", Muon_tightId, &b_Muon_tightId);
  fChain->SetBranchAddress("Muon_trackerMu", Muon_trackerMu, &b_Muon_trackerMu);
  fChain->SetBranchAddress("nPhoton", &nPhoton, &b_nPhoton);
  fChain->SetBranchAddress("Photon_eCorr", Photon_eCorr, &b_Photon_eCorr);
  fChain->SetBranchAddress("Photon_energyErr", Photon_energyErr, &b_Photon_energyErr);
  fChain->SetBranchAddress("Photon_eta", Photon_eta, &b_Photon_eta);
  fChain->SetBranchAddress("Photon_hoe", Photon_hoe, &b_Photon_hoe);
  fChain->SetBranchAddress("Photon_mass", Photon_mass, &b_Photon_mass);
  fChain->SetBranchAddress("Photon_mvaID", Photon_mvaID, &b_Photon_mvaID);
  fChain->SetBranchAddress("Photon_pfRelIso03_all", Photon_pfRelIso03_all, &b_Photon_pfRelIso03_all);
  fChain->SetBranchAddress("Photon_pfRelIso03_chg", Photon_pfRelIso03_chg, &b_Photon_pfRelIso03_chg);
  fChain->SetBranchAddress("Photon_phi", Photon_phi, &b_Photon_phi);
  fChain->SetBranchAddress("Photon_pt", Photon_pt, &b_Photon_pt);
  fChain->SetBranchAddress("Photon_r9", Photon_r9, &b_Photon_r9);
  fChain->SetBranchAddress("Photon_sieie", Photon_sieie, &b_Photon_sieie);
  fChain->SetBranchAddress("Photon_charge", Photon_charge, &b_Photon_charge);
  fChain->SetBranchAddress("Photon_cutBased", Photon_cutBased, &b_Photon_cutBased);
  fChain->SetBranchAddress("Photon_electronIdx", Photon_electronIdx, &b_Photon_electronIdx);
  fChain->SetBranchAddress("Photon_jetIdx", Photon_jetIdx, &b_Photon_jetIdx);
  fChain->SetBranchAddress("Photon_pdgId", Photon_pdgId, &b_Photon_pdgId);
  fChain->SetBranchAddress("Photon_vidNestedWPBitmap", Photon_vidNestedWPBitmap, &b_Photon_vidNestedWPBitmap);
  fChain->SetBranchAddress("Photon_electronVeto", Photon_electronVeto, &b_Photon_electronVeto);
  fChain->SetBranchAddress("Photon_mvaID_WP80", Photon_mvaID_WP80, &b_Photon_mvaID_WP80);
  fChain->SetBranchAddress("Photon_mvaID_WP90", Photon_mvaID_WP90, &b_Photon_mvaID_WP90);
  fChain->SetBranchAddress("Photon_pixelSeed", Photon_pixelSeed, &b_Photon_pixelSeed);
  fChain->SetBranchAddress("Pileup_nPU", &Pileup_nPU, &b_Pileup_nPU);
  fChain->SetBranchAddress("Pileup_nTrueInt", &Pileup_nTrueInt, &b_Pileup_nTrueInt);
  fChain->SetBranchAddress("PuppiMET_phi", &PuppiMET_phi, &b_PuppiMET_phi);
  fChain->SetBranchAddress("PuppiMET_pt", &PuppiMET_pt, &b_PuppiMET_pt);
  fChain->SetBranchAddress("PuppiMET_sumEt", &PuppiMET_sumEt, &b_PuppiMET_sumEt);
  fChain->SetBranchAddress("RawMET_phi", &RawMET_phi, &b_RawMET_phi);
  fChain->SetBranchAddress("RawMET_pt", &RawMET_pt, &b_RawMET_pt);
  fChain->SetBranchAddress("RawMET_sumEt", &RawMET_sumEt, &b_RawMET_sumEt);
  fChain->SetBranchAddress("fixedGridRhoFastjetAll", &fixedGridRhoFastjetAll, &b_fixedGridRhoFastjetAll);
  fChain->SetBranchAddress("fixedGridRhoFastjetCentralCalo", &fixedGridRhoFastjetCentralCalo, &b_fixedGridRhoFastjetCentralCalo);
  fChain->SetBranchAddress("fixedGridRhoFastjetCentralNeutral", &fixedGridRhoFastjetCentralNeutral, &b_fixedGridRhoFastjetCentralNeutral);
  fChain->SetBranchAddress("nGenDressedLepton", &nGenDressedLepton, &b_nGenDressedLepton);
  fChain->SetBranchAddress("GenDressedLepton_eta", GenDressedLepton_eta, &b_GenDressedLepton_eta);
  fChain->SetBranchAddress("GenDressedLepton_mass", GenDressedLepton_mass, &b_GenDressedLepton_mass);
  fChain->SetBranchAddress("GenDressedLepton_phi", GenDressedLepton_phi, &b_GenDressedLepton_phi);
  fChain->SetBranchAddress("GenDressedLepton_pt", GenDressedLepton_pt, &b_GenDressedLepton_pt);
  fChain->SetBranchAddress("GenDressedLepton_pdgId", GenDressedLepton_pdgId, &b_GenDressedLepton_pdgId);
  fChain->SetBranchAddress("nSoftActivityJet", &nSoftActivityJet, &b_nSoftActivityJet);
  fChain->SetBranchAddress("SoftActivityJet_eta", SoftActivityJet_eta, &b_SoftActivityJet_eta);
  fChain->SetBranchAddress("SoftActivityJet_phi", SoftActivityJet_phi, &b_SoftActivityJet_phi);
  fChain->SetBranchAddress("SoftActivityJet_pt", SoftActivityJet_pt, &b_SoftActivityJet_pt);
  fChain->SetBranchAddress("SoftActivityJetHT", &SoftActivityJetHT, &b_SoftActivityJetHT);
  fChain->SetBranchAddress("SoftActivityJetHT10", &SoftActivityJetHT10, &b_SoftActivityJetHT10);
  fChain->SetBranchAddress("SoftActivityJetHT2", &SoftActivityJetHT2, &b_SoftActivityJetHT2);
  fChain->SetBranchAddress("SoftActivityJetHT5", &SoftActivityJetHT5, &b_SoftActivityJetHT5);
  fChain->SetBranchAddress("SoftActivityJetNjets10", &SoftActivityJetNjets10, &b_SoftActivityJetNjets10);
  fChain->SetBranchAddress("SoftActivityJetNjets2", &SoftActivityJetNjets2, &b_SoftActivityJetNjets2);
  fChain->SetBranchAddress("SoftActivityJetNjets5", &SoftActivityJetNjets5, &b_SoftActivityJetNjets5);
  fChain->SetBranchAddress("nSubJet", &nSubJet, &b_nSubJet);
  fChain->SetBranchAddress("SubJet_btagCMVA", SubJet_btagCMVA, &b_SubJet_btagCMVA);
  fChain->SetBranchAddress("SubJet_btagCSVV2", SubJet_btagCSVV2, &b_SubJet_btagCSVV2);
  fChain->SetBranchAddress("SubJet_btagDeepB", SubJet_btagDeepB, &b_SubJet_btagDeepB);
  fChain->SetBranchAddress("SubJet_eta", SubJet_eta, &b_SubJet_eta);
  fChain->SetBranchAddress("SubJet_mass", SubJet_mass, &b_SubJet_mass);
  fChain->SetBranchAddress("SubJet_n2b1", SubJet_n2b1, &b_SubJet_n2b1);
  fChain->SetBranchAddress("SubJet_n3b1", SubJet_n3b1, &b_SubJet_n3b1);
  fChain->SetBranchAddress("SubJet_phi", SubJet_phi, &b_SubJet_phi);
  fChain->SetBranchAddress("SubJet_pt", SubJet_pt, &b_SubJet_pt);
  fChain->SetBranchAddress("SubJet_tau1", SubJet_tau1, &b_SubJet_tau1);
  fChain->SetBranchAddress("SubJet_tau2", SubJet_tau2, &b_SubJet_tau2);
  fChain->SetBranchAddress("SubJet_tau3", SubJet_tau3, &b_SubJet_tau3);
  fChain->SetBranchAddress("SubJet_tau4", SubJet_tau4, &b_SubJet_tau4);
  fChain->SetBranchAddress("nTau", &nTau, &b_nTau);
  fChain->SetBranchAddress("Tau_chargedIso", Tau_chargedIso, &b_Tau_chargedIso);
  fChain->SetBranchAddress("Tau_dxy", Tau_dxy, &b_Tau_dxy);
  fChain->SetBranchAddress("Tau_dz", Tau_dz, &b_Tau_dz);
  fChain->SetBranchAddress("Tau_eta", Tau_eta, &b_Tau_eta);
  fChain->SetBranchAddress("Tau_footprintCorr", Tau_footprintCorr, &b_Tau_footprintCorr);
  fChain->SetBranchAddress("Tau_leadTkDeltaEta", Tau_leadTkDeltaEta, &b_Tau_leadTkDeltaEta);
  fChain->SetBranchAddress("Tau_leadTkDeltaPhi", Tau_leadTkDeltaPhi, &b_Tau_leadTkDeltaPhi);
  fChain->SetBranchAddress("Tau_leadTkPtOverTauPt", Tau_leadTkPtOverTauPt, &b_Tau_leadTkPtOverTauPt);
  fChain->SetBranchAddress("Tau_mass", Tau_mass, &b_Tau_mass);
  fChain->SetBranchAddress("Tau_neutralIso", Tau_neutralIso, &b_Tau_neutralIso);
  fChain->SetBranchAddress("Tau_phi", Tau_phi, &b_Tau_phi);
  fChain->SetBranchAddress("Tau_photonsOutsideSignalCone", Tau_photonsOutsideSignalCone, &b_Tau_photonsOutsideSignalCone);
  fChain->SetBranchAddress("Tau_pt", Tau_pt, &b_Tau_pt);
  fChain->SetBranchAddress("Tau_puCorr", Tau_puCorr, &b_Tau_puCorr);
  fChain->SetBranchAddress("Tau_rawAntiEle", Tau_rawAntiEle, &b_Tau_rawAntiEle);
  fChain->SetBranchAddress("Tau_rawIso", Tau_rawIso, &b_Tau_rawIso);
  fChain->SetBranchAddress("Tau_rawMVAnewDM", Tau_rawMVAnewDM, &b_Tau_rawMVAnewDM);
  fChain->SetBranchAddress("Tau_rawMVAoldDM", Tau_rawMVAoldDM, &b_Tau_rawMVAoldDM);
  fChain->SetBranchAddress("Tau_rawMVAoldDMdR03", Tau_rawMVAoldDMdR03, &b_Tau_rawMVAoldDMdR03);
  fChain->SetBranchAddress("Tau_charge", Tau_charge, &b_Tau_charge);
  fChain->SetBranchAddress("Tau_decayMode", Tau_decayMode, &b_Tau_decayMode);
  fChain->SetBranchAddress("Tau_jetIdx", Tau_jetIdx, &b_Tau_jetIdx);
  fChain->SetBranchAddress("Tau_rawAntiEleCat", Tau_rawAntiEleCat, &b_Tau_rawAntiEleCat);
  fChain->SetBranchAddress("Tau_idAntiEle", Tau_idAntiEle, &b_Tau_idAntiEle);
  fChain->SetBranchAddress("Tau_idAntiMu", Tau_idAntiMu, &b_Tau_idAntiMu);
  fChain->SetBranchAddress("Tau_idDecayMode", Tau_idDecayMode, &b_Tau_idDecayMode);
  fChain->SetBranchAddress("Tau_idDecayModeNewDMs", Tau_idDecayModeNewDMs, &b_Tau_idDecayModeNewDMs);
  fChain->SetBranchAddress("Tau_idMVAnewDM", Tau_idMVAnewDM, &b_Tau_idMVAnewDM);
  fChain->SetBranchAddress("Tau_idMVAoldDM", Tau_idMVAoldDM, &b_Tau_idMVAoldDM);
  fChain->SetBranchAddress("Tau_idMVAoldDMdR03", Tau_idMVAoldDMdR03, &b_Tau_idMVAoldDMdR03);
  fChain->SetBranchAddress("TkMET_phi", &TkMET_phi, &b_TkMET_phi);
  fChain->SetBranchAddress("TkMET_pt", &TkMET_pt, &b_TkMET_pt);
  fChain->SetBranchAddress("TkMET_sumEt", &TkMET_sumEt, &b_TkMET_sumEt);
  fChain->SetBranchAddress("nTrigObj", &nTrigObj, &b_nTrigObj);
  fChain->SetBranchAddress("TrigObj_pt", TrigObj_pt, &b_TrigObj_pt);
  fChain->SetBranchAddress("TrigObj_eta", TrigObj_eta, &b_TrigObj_eta);
  fChain->SetBranchAddress("TrigObj_phi", TrigObj_phi, &b_TrigObj_phi);
  fChain->SetBranchAddress("TrigObj_l1pt", TrigObj_l1pt, &b_TrigObj_l1pt);
  fChain->SetBranchAddress("TrigObj_l1pt_2", TrigObj_l1pt_2, &b_TrigObj_l1pt_2);
  fChain->SetBranchAddress("TrigObj_l2pt", TrigObj_l2pt, &b_TrigObj_l2pt);
  fChain->SetBranchAddress("TrigObj_id", TrigObj_id, &b_TrigObj_id);
  fChain->SetBranchAddress("TrigObj_l1iso", TrigObj_l1iso, &b_TrigObj_l1iso);
  fChain->SetBranchAddress("TrigObj_l1charge", TrigObj_l1charge, &b_TrigObj_l1charge);
  fChain->SetBranchAddress("TrigObj_filterBits", TrigObj_filterBits, &b_TrigObj_filterBits);
  fChain->SetBranchAddress("genTtbarId", &genTtbarId, &b_genTtbarId);
  fChain->SetBranchAddress("nOtherPV", &nOtherPV, &b_nOtherPV);
  fChain->SetBranchAddress("OtherPV_z", OtherPV_z, &b_OtherPV_z);
  fChain->SetBranchAddress("PV_ndof", &PV_ndof, &b_PV_ndof);
  fChain->SetBranchAddress("PV_x", &PV_x, &b_PV_x);
  fChain->SetBranchAddress("PV_y", &PV_y, &b_PV_y);
  fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
  fChain->SetBranchAddress("PV_chi2", &PV_chi2, &b_PV_chi2);
  fChain->SetBranchAddress("PV_score", &PV_score, &b_PV_score);
  fChain->SetBranchAddress("PV_npvs", &PV_npvs, &b_PV_npvs);
  fChain->SetBranchAddress("nSV", &nSV, &b_nSV);
  fChain->SetBranchAddress("SV_dlen", SV_dlen, &b_SV_dlen);
  fChain->SetBranchAddress("SV_dlenSig", SV_dlenSig, &b_SV_dlenSig);
  fChain->SetBranchAddress("SV_pAngle", SV_pAngle, &b_SV_pAngle);
  fChain->SetBranchAddress("cmeson_chi2", cmeson_chi2, &b_cmeson_chi2);
  fChain->SetBranchAddress("cmeson_eta", cmeson_eta, &b_cmeson_eta);
  fChain->SetBranchAddress("cmeson_mass", cmeson_mass, &b_cmeson_mass);
  fChain->SetBranchAddress("cmeson_phi", cmeson_phi, &b_cmeson_phi);
  fChain->SetBranchAddress("cmeson_pt", cmeson_pt, &b_cmeson_pt);
  fChain->SetBranchAddress("cmeson_x", cmeson_x, &b_cmeson_x);
  fChain->SetBranchAddress("cmeson_y", cmeson_y, &b_cmeson_y);
  fChain->SetBranchAddress("cmeson_z", cmeson_z, &b_cmeson_z);
  fChain->SetBranchAddress("cmeson_ndof", cmeson_ndof, &b_cmeson_ndof);
  fChain->SetBranchAddress("cmeson_pdgId", cmeson_pdgId, &b_cmeson_pdgId);
  fChain->SetBranchAddress("Electron_genPartIdx", Electron_genPartIdx, &b_Electron_genPartIdx);
  fChain->SetBranchAddress("Electron_genPartFlav", Electron_genPartFlav, &b_Electron_genPartFlav);
  fChain->SetBranchAddress("GenJetAK8_partonFlavour", GenJetAK8_partonFlavour, &b_GenJetAK8_partonFlavour);
  fChain->SetBranchAddress("GenJetAK8_hadronFlavour", GenJetAK8_hadronFlavour, &b_GenJetAK8_hadronFlavour);
  fChain->SetBranchAddress("GenJet_partonFlavour", GenJet_partonFlavour, &b_GenJet_partonFlavour);
  fChain->SetBranchAddress("GenJet_hadronFlavour", GenJet_hadronFlavour, &b_GenJet_hadronFlavour);
  fChain->SetBranchAddress("Jet_genJetIdx", Jet_genJetIdx, &b_Jet_genJetIdx);
  fChain->SetBranchAddress("Jet_hadronFlavour", Jet_hadronFlavour, &b_Jet_hadronFlavour);
  fChain->SetBranchAddress("Jet_partonFlavour", Jet_partonFlavour, &b_Jet_partonFlavour);
  fChain->SetBranchAddress("Muon_genPartIdx", Muon_genPartIdx, &b_Muon_genPartIdx);
  fChain->SetBranchAddress("Muon_genPartFlav", Muon_genPartFlav, &b_Muon_genPartFlav);
  fChain->SetBranchAddress("Photon_genPartIdx", Photon_genPartIdx, &b_Photon_genPartIdx);
  fChain->SetBranchAddress("Photon_genPartFlav", Photon_genPartFlav, &b_Photon_genPartFlav);
  fChain->SetBranchAddress("MET_fiducialGenPhi", &MET_fiducialGenPhi, &b_MET_fiducialGenPhi);
  fChain->SetBranchAddress("MET_fiducialGenPt", &MET_fiducialGenPt, &b_MET_fiducialGenPt);
  fChain->SetBranchAddress("Electron_cleanmask", Electron_cleanmask, &b_Electron_cleanmask);
  fChain->SetBranchAddress("Jet_cleanmask", Jet_cleanmask, &b_Jet_cleanmask);
  fChain->SetBranchAddress("Muon_cleanmask", Muon_cleanmask, &b_Muon_cleanmask);
  fChain->SetBranchAddress("Photon_cleanmask", Photon_cleanmask, &b_Photon_cleanmask);
  fChain->SetBranchAddress("Tau_cleanmask", Tau_cleanmask, &b_Tau_cleanmask);
  fChain->SetBranchAddress("SV_chi2", SV_chi2, &b_SV_chi2);
  fChain->SetBranchAddress("SV_eta", SV_eta, &b_SV_eta);
  fChain->SetBranchAddress("SV_mass", SV_mass, &b_SV_mass);
  fChain->SetBranchAddress("SV_ndof", SV_ndof, &b_SV_ndof);
  fChain->SetBranchAddress("SV_phi", SV_phi, &b_SV_phi);
  fChain->SetBranchAddress("SV_pt", SV_pt, &b_SV_pt);
  fChain->SetBranchAddress("SV_x", SV_x, &b_SV_x);
  fChain->SetBranchAddress("SV_y", SV_y, &b_SV_y);
  fChain->SetBranchAddress("SV_z", SV_z, &b_SV_z);
  fChain->SetBranchAddress("Tau_genPartIdx", Tau_genPartIdx, &b_Tau_genPartIdx);
  fChain->SetBranchAddress("Tau_genPartFlav", Tau_genPartFlav, &b_Tau_genPartFlav);
  fChain->SetBranchAddress("L1simulation_step", &L1simulation_step, &b_L1simulation_step);
  fChain->SetBranchAddress("HLTriggerFirstPath", &HLTriggerFirstPath, &b_HLTriggerFirstPath);
  fChain->SetBranchAddress("HLT_AK8PFJet360_TrimMass30", &HLT_AK8PFJet360_TrimMass30, &b_HLT_AK8PFJet360_TrimMass30);
  fChain->SetBranchAddress("HLT_AK8PFJet400_TrimMass30", &HLT_AK8PFJet400_TrimMass30, &b_HLT_AK8PFJet400_TrimMass30);
  fChain->SetBranchAddress("HLT_AK8PFHT750_TrimMass50", &HLT_AK8PFHT750_TrimMass50, &b_HLT_AK8PFHT750_TrimMass50);
  fChain->SetBranchAddress("HLT_AK8PFHT800_TrimMass50", &HLT_AK8PFHT800_TrimMass50, &b_HLT_AK8PFHT800_TrimMass50);
  fChain->SetBranchAddress("HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20", &HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20, &b_HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20);
  fChain->SetBranchAddress("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087", &HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087, &b_HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087);
  fChain->SetBranchAddress("HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087", &HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087, &b_HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087);
  fChain->SetBranchAddress("HLT_AK8DiPFJet300_200_TrimMass30", &HLT_AK8DiPFJet300_200_TrimMass30, &b_HLT_AK8DiPFJet300_200_TrimMass30);
  fChain->SetBranchAddress("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50", &HLT_AK8PFHT700_TrimR0p1PT0p03Mass50, &b_HLT_AK8PFHT700_TrimR0p1PT0p03Mass50);
  fChain->SetBranchAddress("HLT_AK8PFHT650_TrimR0p1PT0p03Mass50", &HLT_AK8PFHT650_TrimR0p1PT0p03Mass50, &b_HLT_AK8PFHT650_TrimR0p1PT0p03Mass50);
  fChain->SetBranchAddress("HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20", &HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20, &b_HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20);
  fChain->SetBranchAddress("HLT_AK8DiPFJet280_200_TrimMass30", &HLT_AK8DiPFJet280_200_TrimMass30, &b_HLT_AK8DiPFJet280_200_TrimMass30);
  fChain->SetBranchAddress("HLT_AK8DiPFJet250_200_TrimMass30", &HLT_AK8DiPFJet250_200_TrimMass30, &b_HLT_AK8DiPFJet250_200_TrimMass30);
  fChain->SetBranchAddress("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20", &HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20, &b_HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20);
  fChain->SetBranchAddress("HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20", &HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20, &b_HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20);
  fChain->SetBranchAddress("HLT_CaloJet260", &HLT_CaloJet260, &b_HLT_CaloJet260);
  fChain->SetBranchAddress("HLT_CaloJet500_NoJetID", &HLT_CaloJet500_NoJetID, &b_HLT_CaloJet500_NoJetID);
  fChain->SetBranchAddress("HLT_Dimuon13_PsiPrime", &HLT_Dimuon13_PsiPrime, &b_HLT_Dimuon13_PsiPrime);
  fChain->SetBranchAddress("HLT_Dimuon13_Upsilon", &HLT_Dimuon13_Upsilon, &b_HLT_Dimuon13_Upsilon);
  fChain->SetBranchAddress("HLT_Dimuon20_Jpsi", &HLT_Dimuon20_Jpsi, &b_HLT_Dimuon20_Jpsi);
  fChain->SetBranchAddress("HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf", &HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf, &b_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf);
  fChain->SetBranchAddress("HLT_DoubleEle25_CaloIdL_GsfTrkIdVL", &HLT_DoubleEle25_CaloIdL_GsfTrkIdVL, &b_HLT_DoubleEle25_CaloIdL_GsfTrkIdVL);
  fChain->SetBranchAddress("HLT_DoubleEle33_CaloIdL", &HLT_DoubleEle33_CaloIdL, &b_HLT_DoubleEle33_CaloIdL);
  fChain->SetBranchAddress("HLT_DoubleEle33_CaloIdL_MW", &HLT_DoubleEle33_CaloIdL_MW, &b_HLT_DoubleEle33_CaloIdL_MW);
  fChain->SetBranchAddress("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW", &HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW, &b_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW);
  fChain->SetBranchAddress("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL", &HLT_DoubleEle33_CaloIdL_GsfTrkIdVL, &b_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL);
  fChain->SetBranchAddress("HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg", &HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg, &b_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg);
  fChain->SetBranchAddress("HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg", &HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg, &b_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg);
  fChain->SetBranchAddress("HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg", &HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg, &b_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg);
  fChain->SetBranchAddress("HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg", &HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg, &b_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg);
  fChain->SetBranchAddress("HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1", &HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1, &b_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1);
  fChain->SetBranchAddress("HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1", &HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1, &b_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1);
  fChain->SetBranchAddress("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg", &HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg, &b_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg);
  fChain->SetBranchAddress("HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg", &HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg, &b_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg);
  fChain->SetBranchAddress("HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1", &HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1, &b_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1);
  fChain->SetBranchAddress("HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL", &HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL, &b_HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL);
  fChain->SetBranchAddress("HLT_DoubleMu33NoFiltersNoVtx", &HLT_DoubleMu33NoFiltersNoVtx, &b_HLT_DoubleMu33NoFiltersNoVtx);
  fChain->SetBranchAddress("HLT_DoubleMu38NoFiltersNoVtx", &HLT_DoubleMu38NoFiltersNoVtx, &b_HLT_DoubleMu38NoFiltersNoVtx);
  fChain->SetBranchAddress("HLT_DoubleMu23NoFiltersNoVtxDisplaced", &HLT_DoubleMu23NoFiltersNoVtxDisplaced, &b_HLT_DoubleMu23NoFiltersNoVtxDisplaced);
  fChain->SetBranchAddress("HLT_DoubleMu28NoFiltersNoVtxDisplaced", &HLT_DoubleMu28NoFiltersNoVtxDisplaced, &b_HLT_DoubleMu28NoFiltersNoVtxDisplaced);
  fChain->SetBranchAddress("HLT_DoubleMu0", &HLT_DoubleMu0, &b_HLT_DoubleMu0);
  fChain->SetBranchAddress("HLT_DoubleMu4_3_Bs", &HLT_DoubleMu4_3_Bs, &b_HLT_DoubleMu4_3_Bs);
  fChain->SetBranchAddress("HLT_DoubleMu4_3_Jpsi_Displaced", &HLT_DoubleMu4_3_Jpsi_Displaced, &b_HLT_DoubleMu4_3_Jpsi_Displaced);
  fChain->SetBranchAddress("HLT_DoubleMu4_JpsiTrk_Displaced", &HLT_DoubleMu4_JpsiTrk_Displaced, &b_HLT_DoubleMu4_JpsiTrk_Displaced);
  fChain->SetBranchAddress("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced", &HLT_DoubleMu4_LowMassNonResonantTrk_Displaced, &b_HLT_DoubleMu4_LowMassNonResonantTrk_Displaced);
  fChain->SetBranchAddress("HLT_DoubleMu3_Trk_Tau3mu", &HLT_DoubleMu3_Trk_Tau3mu, &b_HLT_DoubleMu3_Trk_Tau3mu);
  fChain->SetBranchAddress("HLT_DoubleMu4_PsiPrimeTrk_Displaced", &HLT_DoubleMu4_PsiPrimeTrk_Displaced, &b_HLT_DoubleMu4_PsiPrimeTrk_Displaced);
  fChain->SetBranchAddress("HLT_Mu7p5_L2Mu2_Jpsi", &HLT_Mu7p5_L2Mu2_Jpsi, &b_HLT_Mu7p5_L2Mu2_Jpsi);
  fChain->SetBranchAddress("HLT_Mu7p5_L2Mu2_Upsilon", &HLT_Mu7p5_L2Mu2_Upsilon, &b_HLT_Mu7p5_L2Mu2_Upsilon);
  fChain->SetBranchAddress("HLT_Mu7p5_Track2_Jpsi", &HLT_Mu7p5_Track2_Jpsi, &b_HLT_Mu7p5_Track2_Jpsi);
  fChain->SetBranchAddress("HLT_Mu7p5_Track3p5_Jpsi", &HLT_Mu7p5_Track3p5_Jpsi, &b_HLT_Mu7p5_Track3p5_Jpsi);
  fChain->SetBranchAddress("HLT_Mu7p5_Track7_Jpsi", &HLT_Mu7p5_Track7_Jpsi, &b_HLT_Mu7p5_Track7_Jpsi);
  fChain->SetBranchAddress("HLT_Mu7p5_Track2_Upsilon", &HLT_Mu7p5_Track2_Upsilon, &b_HLT_Mu7p5_Track2_Upsilon);
  fChain->SetBranchAddress("HLT_Mu7p5_Track3p5_Upsilon", &HLT_Mu7p5_Track3p5_Upsilon, &b_HLT_Mu7p5_Track3p5_Upsilon);
  fChain->SetBranchAddress("HLT_Mu7p5_Track7_Upsilon", &HLT_Mu7p5_Track7_Upsilon, &b_HLT_Mu7p5_Track7_Upsilon);
  fChain->SetBranchAddress("HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing", &HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing, &b_HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing);
  fChain->SetBranchAddress("HLT_Dimuon0er16_Jpsi_NoVertexing", &HLT_Dimuon0er16_Jpsi_NoVertexing, &b_HLT_Dimuon0er16_Jpsi_NoVertexing);
  fChain->SetBranchAddress("HLT_Dimuon6_Jpsi_NoVertexing", &HLT_Dimuon6_Jpsi_NoVertexing, &b_HLT_Dimuon6_Jpsi_NoVertexing);
  fChain->SetBranchAddress("HLT_Photon150", &HLT_Photon150, &b_HLT_Photon150);
  fChain->SetBranchAddress("HLT_Photon90_CaloIdL_HT300", &HLT_Photon90_CaloIdL_HT300, &b_HLT_Photon90_CaloIdL_HT300);
  fChain->SetBranchAddress("HLT_HT250_CaloMET70", &HLT_HT250_CaloMET70, &b_HLT_HT250_CaloMET70);
  fChain->SetBranchAddress("HLT_DoublePhoton60", &HLT_DoublePhoton60, &b_HLT_DoublePhoton60);
  fChain->SetBranchAddress("HLT_DoublePhoton85", &HLT_DoublePhoton85, &b_HLT_DoublePhoton85);
  fChain->SetBranchAddress("HLT_Ele17_Ele8_Gsf", &HLT_Ele17_Ele8_Gsf, &b_HLT_Ele17_Ele8_Gsf);
  fChain->SetBranchAddress("HLT_Ele20_eta2p1_WPLoose_Gsf_LooseIsoPFTau28", &HLT_Ele20_eta2p1_WPLoose_Gsf_LooseIsoPFTau28, &b_HLT_Ele20_eta2p1_WPLoose_Gsf_LooseIsoPFTau28);
  fChain->SetBranchAddress("HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau29", &HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau29, &b_HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau29);
  fChain->SetBranchAddress("HLT_Ele22_eta2p1_WPLoose_Gsf", &HLT_Ele22_eta2p1_WPLoose_Gsf, &b_HLT_Ele22_eta2p1_WPLoose_Gsf);
  fChain->SetBranchAddress("HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1", &HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1, &b_HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1);
  fChain->SetBranchAddress("HLT_Ele23_WPLoose_Gsf", &HLT_Ele23_WPLoose_Gsf, &b_HLT_Ele23_WPLoose_Gsf);
  fChain->SetBranchAddress("HLT_Ele23_WPLoose_Gsf_WHbbBoost", &HLT_Ele23_WPLoose_Gsf_WHbbBoost, &b_HLT_Ele23_WPLoose_Gsf_WHbbBoost);
  fChain->SetBranchAddress("HLT_Ele24_eta2p1_WPLoose_Gsf", &HLT_Ele24_eta2p1_WPLoose_Gsf, &b_HLT_Ele24_eta2p1_WPLoose_Gsf);
  fChain->SetBranchAddress("HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20", &HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20, &b_HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20);
  fChain->SetBranchAddress("HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1", &HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1, &b_HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1);
  fChain->SetBranchAddress("HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau30", &HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau30, &b_HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau30);
  fChain->SetBranchAddress("HLT_Ele25_WPTight_Gsf", &HLT_Ele25_WPTight_Gsf, &b_HLT_Ele25_WPTight_Gsf);
  fChain->SetBranchAddress("HLT_Ele25_eta2p1_WPLoose_Gsf", &HLT_Ele25_eta2p1_WPLoose_Gsf, &b_HLT_Ele25_eta2p1_WPLoose_Gsf);
  fChain->SetBranchAddress("HLT_Ele25_eta2p1_WPTight_Gsf", &HLT_Ele25_eta2p1_WPTight_Gsf, &b_HLT_Ele25_eta2p1_WPTight_Gsf);
  fChain->SetBranchAddress("HLT_Ele27_WPLoose_Gsf", &HLT_Ele27_WPLoose_Gsf, &b_HLT_Ele27_WPLoose_Gsf);
  fChain->SetBranchAddress("HLT_Ele27_WPLoose_Gsf_WHbbBoost", &HLT_Ele27_WPLoose_Gsf_WHbbBoost, &b_HLT_Ele27_WPLoose_Gsf_WHbbBoost);
  fChain->SetBranchAddress("HLT_Ele27_WPTight_Gsf", &HLT_Ele27_WPTight_Gsf, &b_HLT_Ele27_WPTight_Gsf);
  fChain->SetBranchAddress("HLT_Ele27_WPTight_Gsf_L1JetTauSeeded", &HLT_Ele27_WPTight_Gsf_L1JetTauSeeded, &b_HLT_Ele27_WPTight_Gsf_L1JetTauSeeded);
  fChain->SetBranchAddress("HLT_Ele27_eta2p1_WPLoose_Gsf", &HLT_Ele27_eta2p1_WPLoose_Gsf, &b_HLT_Ele27_eta2p1_WPLoose_Gsf);
  fChain->SetBranchAddress("HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1", &HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1, &b_HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1);
  fChain->SetBranchAddress("HLT_Ele27_eta2p1_WPTight_Gsf", &HLT_Ele27_eta2p1_WPTight_Gsf, &b_HLT_Ele27_eta2p1_WPTight_Gsf);
  fChain->SetBranchAddress("HLT_Ele30_WPTight_Gsf", &HLT_Ele30_WPTight_Gsf, &b_HLT_Ele30_WPTight_Gsf);
  fChain->SetBranchAddress("HLT_Ele30_eta2p1_WPLoose_Gsf", &HLT_Ele30_eta2p1_WPLoose_Gsf, &b_HLT_Ele30_eta2p1_WPLoose_Gsf);
  fChain->SetBranchAddress("HLT_Ele30_eta2p1_WPTight_Gsf", &HLT_Ele30_eta2p1_WPTight_Gsf, &b_HLT_Ele30_eta2p1_WPTight_Gsf);
  fChain->SetBranchAddress("HLT_Ele32_WPTight_Gsf", &HLT_Ele32_WPTight_Gsf, &b_HLT_Ele32_WPTight_Gsf);
  fChain->SetBranchAddress("HLT_Ele32_eta2p1_WPLoose_Gsf", &HLT_Ele32_eta2p1_WPLoose_Gsf, &b_HLT_Ele32_eta2p1_WPLoose_Gsf);
  fChain->SetBranchAddress("HLT_Ele32_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1", &HLT_Ele32_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1, &b_HLT_Ele32_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1);
  fChain->SetBranchAddress("HLT_Ele32_eta2p1_WPTight_Gsf", &HLT_Ele32_eta2p1_WPTight_Gsf, &b_HLT_Ele32_eta2p1_WPTight_Gsf);
  fChain->SetBranchAddress("HLT_Ele35_WPLoose_Gsf", &HLT_Ele35_WPLoose_Gsf, &b_HLT_Ele35_WPLoose_Gsf);
  fChain->SetBranchAddress("HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50", &HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50, &b_HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50);
  fChain->SetBranchAddress("HLT_Ele36_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1", &HLT_Ele36_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1, &b_HLT_Ele36_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1);
  fChain->SetBranchAddress("HLT_Ele45_WPLoose_Gsf", &HLT_Ele45_WPLoose_Gsf, &b_HLT_Ele45_WPLoose_Gsf);
  fChain->SetBranchAddress("HLT_Ele45_WPLoose_Gsf_L1JetTauSeeded", &HLT_Ele45_WPLoose_Gsf_L1JetTauSeeded, &b_HLT_Ele45_WPLoose_Gsf_L1JetTauSeeded);
  fChain->SetBranchAddress("HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50", &HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50, &b_HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50);
  fChain->SetBranchAddress("HLT_Ele105_CaloIdVT_GsfTrkIdT", &HLT_Ele105_CaloIdVT_GsfTrkIdT, &b_HLT_Ele105_CaloIdVT_GsfTrkIdT);
  fChain->SetBranchAddress("HLT_Ele30WP60_SC4_Mass55", &HLT_Ele30WP60_SC4_Mass55, &b_HLT_Ele30WP60_SC4_Mass55);
  fChain->SetBranchAddress("HLT_Ele30WP60_Ele8_Mass55", &HLT_Ele30WP60_Ele8_Mass55, &b_HLT_Ele30WP60_Ele8_Mass55);
  fChain->SetBranchAddress("HLT_HT200", &HLT_HT200, &b_HLT_HT200);
  fChain->SetBranchAddress("HLT_HT275", &HLT_HT275, &b_HLT_HT275);
  fChain->SetBranchAddress("HLT_HT325", &HLT_HT325, &b_HLT_HT325);
  fChain->SetBranchAddress("HLT_HT425", &HLT_HT425, &b_HLT_HT425);
  fChain->SetBranchAddress("HLT_HT575", &HLT_HT575, &b_HLT_HT575);
  fChain->SetBranchAddress("HLT_HT410to430", &HLT_HT410to430, &b_HLT_HT410to430);
  fChain->SetBranchAddress("HLT_HT430to450", &HLT_HT430to450, &b_HLT_HT430to450);
  fChain->SetBranchAddress("HLT_HT450to470", &HLT_HT450to470, &b_HLT_HT450to470);
  fChain->SetBranchAddress("HLT_HT470to500", &HLT_HT470to500, &b_HLT_HT470to500);
  fChain->SetBranchAddress("HLT_HT500to550", &HLT_HT500to550, &b_HLT_HT500to550);
  fChain->SetBranchAddress("HLT_HT550to650", &HLT_HT550to650, &b_HLT_HT550to650);
  fChain->SetBranchAddress("HLT_HT650", &HLT_HT650, &b_HLT_HT650);
  fChain->SetBranchAddress("HLT_Mu16_eta2p1_MET30", &HLT_Mu16_eta2p1_MET30, &b_HLT_Mu16_eta2p1_MET30);
  fChain->SetBranchAddress("HLT_IsoMu16_eta2p1_MET30", &HLT_IsoMu16_eta2p1_MET30, &b_HLT_IsoMu16_eta2p1_MET30);
  fChain->SetBranchAddress("HLT_IsoMu16_eta2p1_MET30_LooseIsoPFTau50_Trk30_eta2p1", &HLT_IsoMu16_eta2p1_MET30_LooseIsoPFTau50_Trk30_eta2p1, &b_HLT_IsoMu16_eta2p1_MET30_LooseIsoPFTau50_Trk30_eta2p1);
  fChain->SetBranchAddress("HLT_IsoMu17_eta2p1", &HLT_IsoMu17_eta2p1, &b_HLT_IsoMu17_eta2p1);
  fChain->SetBranchAddress("HLT_IsoMu17_eta2p1_LooseIsoPFTau20", &HLT_IsoMu17_eta2p1_LooseIsoPFTau20, &b_HLT_IsoMu17_eta2p1_LooseIsoPFTau20);
  fChain->SetBranchAddress("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1", &HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1, &b_HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1);
  fChain->SetBranchAddress("HLT_DoubleIsoMu17_eta2p1", &HLT_DoubleIsoMu17_eta2p1, &b_HLT_DoubleIsoMu17_eta2p1);
  fChain->SetBranchAddress("HLT_DoubleIsoMu17_eta2p1_noDzCut", &HLT_DoubleIsoMu17_eta2p1_noDzCut, &b_HLT_DoubleIsoMu17_eta2p1_noDzCut);
  fChain->SetBranchAddress("HLT_IsoMu18", &HLT_IsoMu18, &b_HLT_IsoMu18);
  fChain->SetBranchAddress("HLT_IsoMu19_eta2p1_LooseIsoPFTau20", &HLT_IsoMu19_eta2p1_LooseIsoPFTau20, &b_HLT_IsoMu19_eta2p1_LooseIsoPFTau20);
  fChain->SetBranchAddress("HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1", &HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1, &b_HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1);
  fChain->SetBranchAddress("HLT_IsoMu19_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg", &HLT_IsoMu19_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg, &b_HLT_IsoMu19_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg);
  fChain->SetBranchAddress("HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20", &HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20, &b_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20);
  fChain->SetBranchAddress("HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg", &HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg, &b_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg);
  fChain->SetBranchAddress("HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg", &HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg, &b_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg);
  fChain->SetBranchAddress("HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg", &HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg, &b_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg);
  fChain->SetBranchAddress("HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg", &HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg, &b_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg);
  fChain->SetBranchAddress("HLT_IsoMu20", &HLT_IsoMu20, &b_HLT_IsoMu20);
  fChain->SetBranchAddress("HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1", &HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1, &b_HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1);
  fChain->SetBranchAddress("HLT_IsoMu21_eta2p1_LooseIsoPFTau50_Trk30_eta2p1_SingleL1", &HLT_IsoMu21_eta2p1_LooseIsoPFTau50_Trk30_eta2p1_SingleL1, &b_HLT_IsoMu21_eta2p1_LooseIsoPFTau50_Trk30_eta2p1_SingleL1);
  fChain->SetBranchAddress("HLT_IsoMu21_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg", &HLT_IsoMu21_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg, &b_HLT_IsoMu21_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg);
  fChain->SetBranchAddress("HLT_IsoMu22", &HLT_IsoMu22, &b_HLT_IsoMu22);
  fChain->SetBranchAddress("HLT_IsoMu22_eta2p1", &HLT_IsoMu22_eta2p1, &b_HLT_IsoMu22_eta2p1);
  fChain->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24, &b_HLT_IsoMu24);
  fChain->SetBranchAddress("HLT_IsoMu27", &HLT_IsoMu27, &b_HLT_IsoMu27);
  fChain->SetBranchAddress("HLT_IsoTkMu18", &HLT_IsoTkMu18, &b_HLT_IsoTkMu18);
  fChain->SetBranchAddress("HLT_IsoTkMu20", &HLT_IsoTkMu20, &b_HLT_IsoTkMu20);
  fChain->SetBranchAddress("HLT_IsoTkMu22", &HLT_IsoTkMu22, &b_HLT_IsoTkMu22);
  fChain->SetBranchAddress("HLT_IsoTkMu22_eta2p1", &HLT_IsoTkMu22_eta2p1, &b_HLT_IsoTkMu22_eta2p1);
  fChain->SetBranchAddress("HLT_IsoTkMu24", &HLT_IsoTkMu24, &b_HLT_IsoTkMu24);
  fChain->SetBranchAddress("HLT_IsoTkMu27", &HLT_IsoTkMu27, &b_HLT_IsoTkMu27);
  fChain->SetBranchAddress("HLT_JetE30_NoBPTX3BX", &HLT_JetE30_NoBPTX3BX, &b_HLT_JetE30_NoBPTX3BX);
  fChain->SetBranchAddress("HLT_JetE30_NoBPTX", &HLT_JetE30_NoBPTX, &b_HLT_JetE30_NoBPTX);
  fChain->SetBranchAddress("HLT_JetE50_NoBPTX3BX", &HLT_JetE50_NoBPTX3BX, &b_HLT_JetE50_NoBPTX3BX);
  fChain->SetBranchAddress("HLT_JetE70_NoBPTX3BX", &HLT_JetE70_NoBPTX3BX, &b_HLT_JetE70_NoBPTX3BX);
  fChain->SetBranchAddress("HLT_L1SingleMu18", &HLT_L1SingleMu18, &b_HLT_L1SingleMu18);
  fChain->SetBranchAddress("HLT_L2Mu10", &HLT_L2Mu10, &b_HLT_L2Mu10);
  fChain->SetBranchAddress("HLT_L1SingleMuOpen", &HLT_L1SingleMuOpen, &b_HLT_L1SingleMuOpen);
  fChain->SetBranchAddress("HLT_L1SingleMuOpen_DT", &HLT_L1SingleMuOpen_DT, &b_HLT_L1SingleMuOpen_DT);
  fChain->SetBranchAddress("HLT_L2DoubleMu23_NoVertex", &HLT_L2DoubleMu23_NoVertex, &b_HLT_L2DoubleMu23_NoVertex);
  fChain->SetBranchAddress("HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10", &HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10, &b_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10);
  fChain->SetBranchAddress("HLT_L2DoubleMu38_NoVertex_2Cha_Angle2p5_Mass10", &HLT_L2DoubleMu38_NoVertex_2Cha_Angle2p5_Mass10, &b_HLT_L2DoubleMu38_NoVertex_2Cha_Angle2p5_Mass10);
  fChain->SetBranchAddress("HLT_L2Mu10_NoVertex_NoBPTX3BX", &HLT_L2Mu10_NoVertex_NoBPTX3BX, &b_HLT_L2Mu10_NoVertex_NoBPTX3BX);
  fChain->SetBranchAddress("HLT_L2Mu10_NoVertex_NoBPTX", &HLT_L2Mu10_NoVertex_NoBPTX, &b_HLT_L2Mu10_NoVertex_NoBPTX);
  fChain->SetBranchAddress("HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX", &HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX, &b_HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX);
  fChain->SetBranchAddress("HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX", &HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX, &b_HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX);
  fChain->SetBranchAddress("HLT_LooseIsoPFTau50_Trk30_eta2p1", &HLT_LooseIsoPFTau50_Trk30_eta2p1, &b_HLT_LooseIsoPFTau50_Trk30_eta2p1);
  fChain->SetBranchAddress("HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80", &HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80, &b_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80);
  fChain->SetBranchAddress("HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90", &HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90, &b_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90);
  fChain->SetBranchAddress("HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110", &HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110, &b_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110);
  fChain->SetBranchAddress("HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120", &HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120, &b_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120);
  fChain->SetBranchAddress("HLT_PFTau120_eta2p1", &HLT_PFTau120_eta2p1, &b_HLT_PFTau120_eta2p1);
  fChain->SetBranchAddress("HLT_PFTau140_eta2p1", &HLT_PFTau140_eta2p1, &b_HLT_PFTau140_eta2p1);
  fChain->SetBranchAddress("HLT_VLooseIsoPFTau120_Trk50_eta2p1", &HLT_VLooseIsoPFTau120_Trk50_eta2p1, &b_HLT_VLooseIsoPFTau120_Trk50_eta2p1);
  fChain->SetBranchAddress("HLT_VLooseIsoPFTau140_Trk50_eta2p1", &HLT_VLooseIsoPFTau140_Trk50_eta2p1, &b_HLT_VLooseIsoPFTau140_Trk50_eta2p1);
  fChain->SetBranchAddress("HLT_Mu17_Mu8", &HLT_Mu17_Mu8, &b_HLT_Mu17_Mu8);
  fChain->SetBranchAddress("HLT_Mu17_Mu8_DZ", &HLT_Mu17_Mu8_DZ, &b_HLT_Mu17_Mu8_DZ);
  fChain->SetBranchAddress("HLT_Mu17_Mu8_SameSign", &HLT_Mu17_Mu8_SameSign, &b_HLT_Mu17_Mu8_SameSign);
  fChain->SetBranchAddress("HLT_Mu17_Mu8_SameSign_DZ", &HLT_Mu17_Mu8_SameSign_DZ, &b_HLT_Mu17_Mu8_SameSign_DZ);
  fChain->SetBranchAddress("HLT_Mu20_Mu10", &HLT_Mu20_Mu10, &b_HLT_Mu20_Mu10);
  fChain->SetBranchAddress("HLT_Mu20_Mu10_DZ", &HLT_Mu20_Mu10_DZ, &b_HLT_Mu20_Mu10_DZ);
  fChain->SetBranchAddress("HLT_Mu20_Mu10_SameSign", &HLT_Mu20_Mu10_SameSign, &b_HLT_Mu20_Mu10_SameSign);
  fChain->SetBranchAddress("HLT_Mu20_Mu10_SameSign_DZ", &HLT_Mu20_Mu10_SameSign_DZ, &b_HLT_Mu20_Mu10_SameSign_DZ);
  fChain->SetBranchAddress("HLT_Mu17_TkMu8_DZ", &HLT_Mu17_TkMu8_DZ, &b_HLT_Mu17_TkMu8_DZ);
  fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL);
  fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ);
  fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL, &b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL);
  fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ, &b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ);
  fChain->SetBranchAddress("HLT_Mu25_TkMu0_dEta18_Onia", &HLT_Mu25_TkMu0_dEta18_Onia, &b_HLT_Mu25_TkMu0_dEta18_Onia);
  fChain->SetBranchAddress("HLT_Mu27_TkMu8", &HLT_Mu27_TkMu8, &b_HLT_Mu27_TkMu8);
  fChain->SetBranchAddress("HLT_Mu30_TkMu11", &HLT_Mu30_TkMu11, &b_HLT_Mu30_TkMu11);
  fChain->SetBranchAddress("HLT_Mu30_eta2p1_PFJet150_PFJet50", &HLT_Mu30_eta2p1_PFJet150_PFJet50, &b_HLT_Mu30_eta2p1_PFJet150_PFJet50);
  fChain->SetBranchAddress("HLT_Mu40_TkMu11", &HLT_Mu40_TkMu11, &b_HLT_Mu40_TkMu11);
  fChain->SetBranchAddress("HLT_Mu40_eta2p1_PFJet200_PFJet50", &HLT_Mu40_eta2p1_PFJet200_PFJet50, &b_HLT_Mu40_eta2p1_PFJet200_PFJet50);
  fChain->SetBranchAddress("HLT_Mu20", &HLT_Mu20, &b_HLT_Mu20);
  fChain->SetBranchAddress("HLT_TkMu17", &HLT_TkMu17, &b_HLT_TkMu17);
  fChain->SetBranchAddress("HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL", &HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL, &b_HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL);
  fChain->SetBranchAddress("HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ", &HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ, &b_HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ);
  fChain->SetBranchAddress("HLT_TkMu20", &HLT_TkMu20, &b_HLT_TkMu20);
  fChain->SetBranchAddress("HLT_Mu24_eta2p1", &HLT_Mu24_eta2p1, &b_HLT_Mu24_eta2p1);
  fChain->SetBranchAddress("HLT_TkMu24_eta2p1", &HLT_TkMu24_eta2p1, &b_HLT_TkMu24_eta2p1);
  fChain->SetBranchAddress("HLT_Mu27", &HLT_Mu27, &b_HLT_Mu27);
  fChain->SetBranchAddress("HLT_TkMu27", &HLT_TkMu27, &b_HLT_TkMu27);
  fChain->SetBranchAddress("HLT_Mu45_eta2p1", &HLT_Mu45_eta2p1, &b_HLT_Mu45_eta2p1);
  fChain->SetBranchAddress("HLT_Mu50", &HLT_Mu50, &b_HLT_Mu50);
  fChain->SetBranchAddress("HLT_TkMu50", &HLT_TkMu50, &b_HLT_TkMu50);
  fChain->SetBranchAddress("HLT_Mu38NoFiltersNoVtx_Photon38_CaloIdL", &HLT_Mu38NoFiltersNoVtx_Photon38_CaloIdL, &b_HLT_Mu38NoFiltersNoVtx_Photon38_CaloIdL);
  fChain->SetBranchAddress("HLT_Mu42NoFiltersNoVtx_Photon42_CaloIdL", &HLT_Mu42NoFiltersNoVtx_Photon42_CaloIdL, &b_HLT_Mu42NoFiltersNoVtx_Photon42_CaloIdL);
  fChain->SetBranchAddress("HLT_Mu28NoFiltersNoVtxDisplaced_Photon28_CaloIdL", &HLT_Mu28NoFiltersNoVtxDisplaced_Photon28_CaloIdL, &b_HLT_Mu28NoFiltersNoVtxDisplaced_Photon28_CaloIdL);
  fChain->SetBranchAddress("HLT_Mu33NoFiltersNoVtxDisplaced_Photon33_CaloIdL", &HLT_Mu33NoFiltersNoVtxDisplaced_Photon33_CaloIdL, &b_HLT_Mu33NoFiltersNoVtxDisplaced_Photon33_CaloIdL);
  fChain->SetBranchAddress("HLT_Mu23NoFiltersNoVtx_Photon23_CaloIdL", &HLT_Mu23NoFiltersNoVtx_Photon23_CaloIdL, &b_HLT_Mu23NoFiltersNoVtx_Photon23_CaloIdL);
  fChain->SetBranchAddress("HLT_DoubleMu18NoFiltersNoVtx", &HLT_DoubleMu18NoFiltersNoVtx, &b_HLT_DoubleMu18NoFiltersNoVtx);
  fChain->SetBranchAddress("HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Tight", &HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Tight, &b_HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Tight);
  fChain->SetBranchAddress("HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Loose", &HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Loose, &b_HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Loose);
  fChain->SetBranchAddress("HLT_Mu28NoFiltersNoVtx_DisplacedJet40_Loose", &HLT_Mu28NoFiltersNoVtx_DisplacedJet40_Loose, &b_HLT_Mu28NoFiltersNoVtx_DisplacedJet40_Loose);
  fChain->SetBranchAddress("HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Tight", &HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Tight, &b_HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Tight);
  fChain->SetBranchAddress("HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Loose", &HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Loose, &b_HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Loose);
  fChain->SetBranchAddress("HLT_Mu38NoFiltersNoVtx_DisplacedJet60_Loose", &HLT_Mu38NoFiltersNoVtx_DisplacedJet60_Loose, &b_HLT_Mu38NoFiltersNoVtx_DisplacedJet60_Loose);
  fChain->SetBranchAddress("HLT_Mu28NoFiltersNoVtx_CentralCaloJet40", &HLT_Mu28NoFiltersNoVtx_CentralCaloJet40, &b_HLT_Mu28NoFiltersNoVtx_CentralCaloJet40);
  fChain->SetBranchAddress("HLT_PFHT300_PFMET100", &HLT_PFHT300_PFMET100, &b_HLT_PFHT300_PFMET100);
  fChain->SetBranchAddress("HLT_PFHT300_PFMET110", &HLT_PFHT300_PFMET110, &b_HLT_PFHT300_PFMET110);
  fChain->SetBranchAddress("HLT_PFHT550_4JetPt50", &HLT_PFHT550_4JetPt50, &b_HLT_PFHT550_4JetPt50);
  fChain->SetBranchAddress("HLT_PFHT650_4JetPt50", &HLT_PFHT650_4JetPt50, &b_HLT_PFHT650_4JetPt50);
  fChain->SetBranchAddress("HLT_PFHT750_4JetPt50", &HLT_PFHT750_4JetPt50, &b_HLT_PFHT750_4JetPt50);
  fChain->SetBranchAddress("HLT_PFHT750_4JetPt70", &HLT_PFHT750_4JetPt70, &b_HLT_PFHT750_4JetPt70);
  fChain->SetBranchAddress("HLT_PFHT750_4JetPt80", &HLT_PFHT750_4JetPt80, &b_HLT_PFHT750_4JetPt80);
  fChain->SetBranchAddress("HLT_PFHT800_4JetPt50", &HLT_PFHT800_4JetPt50, &b_HLT_PFHT800_4JetPt50);
  fChain->SetBranchAddress("HLT_PFHT850_4JetPt50", &HLT_PFHT850_4JetPt50, &b_HLT_PFHT850_4JetPt50);
  fChain->SetBranchAddress("HLT_PFJet15_NoCaloMatched", &HLT_PFJet15_NoCaloMatched, &b_HLT_PFJet15_NoCaloMatched);
  fChain->SetBranchAddress("HLT_PFJet25_NoCaloMatched", &HLT_PFJet25_NoCaloMatched, &b_HLT_PFJet25_NoCaloMatched);
  fChain->SetBranchAddress("HLT_DiPFJet15_NoCaloMatched", &HLT_DiPFJet15_NoCaloMatched, &b_HLT_DiPFJet15_NoCaloMatched);
  fChain->SetBranchAddress("HLT_DiPFJet25_NoCaloMatched", &HLT_DiPFJet25_NoCaloMatched, &b_HLT_DiPFJet25_NoCaloMatched);
  fChain->SetBranchAddress("HLT_DiPFJet15_FBEta3_NoCaloMatched", &HLT_DiPFJet15_FBEta3_NoCaloMatched, &b_HLT_DiPFJet15_FBEta3_NoCaloMatched);
  fChain->SetBranchAddress("HLT_DiPFJet25_FBEta3_NoCaloMatched", &HLT_DiPFJet25_FBEta3_NoCaloMatched, &b_HLT_DiPFJet25_FBEta3_NoCaloMatched);
  fChain->SetBranchAddress("HLT_DiPFJetAve15_HFJEC", &HLT_DiPFJetAve15_HFJEC, &b_HLT_DiPFJetAve15_HFJEC);
  fChain->SetBranchAddress("HLT_DiPFJetAve25_HFJEC", &HLT_DiPFJetAve25_HFJEC, &b_HLT_DiPFJetAve25_HFJEC);
  fChain->SetBranchAddress("HLT_DiPFJetAve35_HFJEC", &HLT_DiPFJetAve35_HFJEC, &b_HLT_DiPFJetAve35_HFJEC);
  fChain->SetBranchAddress("HLT_AK8PFJet40", &HLT_AK8PFJet40, &b_HLT_AK8PFJet40);
  fChain->SetBranchAddress("HLT_AK8PFJet60", &HLT_AK8PFJet60, &b_HLT_AK8PFJet60);
  fChain->SetBranchAddress("HLT_AK8PFJet80", &HLT_AK8PFJet80, &b_HLT_AK8PFJet80);
  fChain->SetBranchAddress("HLT_AK8PFJet140", &HLT_AK8PFJet140, &b_HLT_AK8PFJet140);
  fChain->SetBranchAddress("HLT_AK8PFJet200", &HLT_AK8PFJet200, &b_HLT_AK8PFJet200);
  fChain->SetBranchAddress("HLT_AK8PFJet260", &HLT_AK8PFJet260, &b_HLT_AK8PFJet260);
  fChain->SetBranchAddress("HLT_AK8PFJet320", &HLT_AK8PFJet320, &b_HLT_AK8PFJet320);
  fChain->SetBranchAddress("HLT_AK8PFJet400", &HLT_AK8PFJet400, &b_HLT_AK8PFJet400);
  fChain->SetBranchAddress("HLT_AK8PFJet450", &HLT_AK8PFJet450, &b_HLT_AK8PFJet450);
  fChain->SetBranchAddress("HLT_AK8PFJet500", &HLT_AK8PFJet500, &b_HLT_AK8PFJet500);
  fChain->SetBranchAddress("HLT_PFJet40", &HLT_PFJet40, &b_HLT_PFJet40);
  fChain->SetBranchAddress("HLT_PFJet60", &HLT_PFJet60, &b_HLT_PFJet60);
  fChain->SetBranchAddress("HLT_PFJet80", &HLT_PFJet80, &b_HLT_PFJet80);
  fChain->SetBranchAddress("HLT_PFJet140", &HLT_PFJet140, &b_HLT_PFJet140);
  fChain->SetBranchAddress("HLT_PFJet200", &HLT_PFJet200, &b_HLT_PFJet200);
  fChain->SetBranchAddress("HLT_PFJet260", &HLT_PFJet260, &b_HLT_PFJet260);
  fChain->SetBranchAddress("HLT_PFJet320", &HLT_PFJet320, &b_HLT_PFJet320);
  fChain->SetBranchAddress("HLT_PFJet400", &HLT_PFJet400, &b_HLT_PFJet400);
  fChain->SetBranchAddress("HLT_PFJet450", &HLT_PFJet450, &b_HLT_PFJet450);
  fChain->SetBranchAddress("HLT_PFJet500", &HLT_PFJet500, &b_HLT_PFJet500);
  fChain->SetBranchAddress("HLT_DiPFJetAve40", &HLT_DiPFJetAve40, &b_HLT_DiPFJetAve40);
  fChain->SetBranchAddress("HLT_DiPFJetAve60", &HLT_DiPFJetAve60, &b_HLT_DiPFJetAve60);
  fChain->SetBranchAddress("HLT_DiPFJetAve80", &HLT_DiPFJetAve80, &b_HLT_DiPFJetAve80);
  fChain->SetBranchAddress("HLT_DiPFJetAve140", &HLT_DiPFJetAve140, &b_HLT_DiPFJetAve140);
  fChain->SetBranchAddress("HLT_DiPFJetAve200", &HLT_DiPFJetAve200, &b_HLT_DiPFJetAve200);
  fChain->SetBranchAddress("HLT_DiPFJetAve260", &HLT_DiPFJetAve260, &b_HLT_DiPFJetAve260);
  fChain->SetBranchAddress("HLT_DiPFJetAve320", &HLT_DiPFJetAve320, &b_HLT_DiPFJetAve320);
  fChain->SetBranchAddress("HLT_DiPFJetAve400", &HLT_DiPFJetAve400, &b_HLT_DiPFJetAve400);
  fChain->SetBranchAddress("HLT_DiPFJetAve500", &HLT_DiPFJetAve500, &b_HLT_DiPFJetAve500);
  fChain->SetBranchAddress("HLT_DiPFJetAve60_HFJEC", &HLT_DiPFJetAve60_HFJEC, &b_HLT_DiPFJetAve60_HFJEC);
  fChain->SetBranchAddress("HLT_DiPFJetAve80_HFJEC", &HLT_DiPFJetAve80_HFJEC, &b_HLT_DiPFJetAve80_HFJEC);
  fChain->SetBranchAddress("HLT_DiPFJetAve100_HFJEC", &HLT_DiPFJetAve100_HFJEC, &b_HLT_DiPFJetAve100_HFJEC);
  fChain->SetBranchAddress("HLT_DiPFJetAve160_HFJEC", &HLT_DiPFJetAve160_HFJEC, &b_HLT_DiPFJetAve160_HFJEC);
  fChain->SetBranchAddress("HLT_DiPFJetAve220_HFJEC", &HLT_DiPFJetAve220_HFJEC, &b_HLT_DiPFJetAve220_HFJEC);
  fChain->SetBranchAddress("HLT_DiPFJetAve300_HFJEC", &HLT_DiPFJetAve300_HFJEC, &b_HLT_DiPFJetAve300_HFJEC);
  fChain->SetBranchAddress("HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu140", &HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu140, &b_HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu140);
  fChain->SetBranchAddress("HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu80", &HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu80, &b_HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu80);
  fChain->SetBranchAddress("HLT_DiCentralPFJet170", &HLT_DiCentralPFJet170, &b_HLT_DiCentralPFJet170);
  fChain->SetBranchAddress("HLT_SingleCentralPFJet170_CFMax0p1", &HLT_SingleCentralPFJet170_CFMax0p1, &b_HLT_SingleCentralPFJet170_CFMax0p1);
  fChain->SetBranchAddress("HLT_DiCentralPFJet170_CFMax0p1", &HLT_DiCentralPFJet170_CFMax0p1, &b_HLT_DiCentralPFJet170_CFMax0p1);
  fChain->SetBranchAddress("HLT_DiCentralPFJet220_CFMax0p3", &HLT_DiCentralPFJet220_CFMax0p3, &b_HLT_DiCentralPFJet220_CFMax0p3);
  fChain->SetBranchAddress("HLT_DiCentralPFJet330_CFMax0p5", &HLT_DiCentralPFJet330_CFMax0p5, &b_HLT_DiCentralPFJet330_CFMax0p5);
  fChain->SetBranchAddress("HLT_DiCentralPFJet430", &HLT_DiCentralPFJet430, &b_HLT_DiCentralPFJet430);
  fChain->SetBranchAddress("HLT_PFHT125", &HLT_PFHT125, &b_HLT_PFHT125);
  fChain->SetBranchAddress("HLT_PFHT200", &HLT_PFHT200, &b_HLT_PFHT200);
  fChain->SetBranchAddress("HLT_PFHT250", &HLT_PFHT250, &b_HLT_PFHT250);
  fChain->SetBranchAddress("HLT_PFHT300", &HLT_PFHT300, &b_HLT_PFHT300);
  fChain->SetBranchAddress("HLT_PFHT350", &HLT_PFHT350, &b_HLT_PFHT350);
  fChain->SetBranchAddress("HLT_PFHT400", &HLT_PFHT400, &b_HLT_PFHT400);
  fChain->SetBranchAddress("HLT_PFHT475", &HLT_PFHT475, &b_HLT_PFHT475);
  fChain->SetBranchAddress("HLT_PFHT600", &HLT_PFHT600, &b_HLT_PFHT600);
  fChain->SetBranchAddress("HLT_PFHT650", &HLT_PFHT650, &b_HLT_PFHT650);
  fChain->SetBranchAddress("HLT_PFHT800", &HLT_PFHT800, &b_HLT_PFHT800);
  fChain->SetBranchAddress("HLT_PFHT900", &HLT_PFHT900, &b_HLT_PFHT900);
  fChain->SetBranchAddress("HLT_PFHT200_PFAlphaT0p51", &HLT_PFHT200_PFAlphaT0p51, &b_HLT_PFHT200_PFAlphaT0p51);
  fChain->SetBranchAddress("HLT_PFHT200_DiPFJetAve90_PFAlphaT0p57", &HLT_PFHT200_DiPFJetAve90_PFAlphaT0p57, &b_HLT_PFHT200_DiPFJetAve90_PFAlphaT0p57);
  fChain->SetBranchAddress("HLT_PFHT200_DiPFJetAve90_PFAlphaT0p63", &HLT_PFHT200_DiPFJetAve90_PFAlphaT0p63, &b_HLT_PFHT200_DiPFJetAve90_PFAlphaT0p63);
  fChain->SetBranchAddress("HLT_PFHT250_DiPFJetAve90_PFAlphaT0p55", &HLT_PFHT250_DiPFJetAve90_PFAlphaT0p55, &b_HLT_PFHT250_DiPFJetAve90_PFAlphaT0p55);
  fChain->SetBranchAddress("HLT_PFHT250_DiPFJetAve90_PFAlphaT0p58", &HLT_PFHT250_DiPFJetAve90_PFAlphaT0p58, &b_HLT_PFHT250_DiPFJetAve90_PFAlphaT0p58);
  fChain->SetBranchAddress("HLT_PFHT300_DiPFJetAve90_PFAlphaT0p53", &HLT_PFHT300_DiPFJetAve90_PFAlphaT0p53, &b_HLT_PFHT300_DiPFJetAve90_PFAlphaT0p53);
  fChain->SetBranchAddress("HLT_PFHT300_DiPFJetAve90_PFAlphaT0p54", &HLT_PFHT300_DiPFJetAve90_PFAlphaT0p54, &b_HLT_PFHT300_DiPFJetAve90_PFAlphaT0p54);
  fChain->SetBranchAddress("HLT_PFHT350_DiPFJetAve90_PFAlphaT0p52", &HLT_PFHT350_DiPFJetAve90_PFAlphaT0p52, &b_HLT_PFHT350_DiPFJetAve90_PFAlphaT0p52);
  fChain->SetBranchAddress("HLT_PFHT350_DiPFJetAve90_PFAlphaT0p53", &HLT_PFHT350_DiPFJetAve90_PFAlphaT0p53, &b_HLT_PFHT350_DiPFJetAve90_PFAlphaT0p53);
  fChain->SetBranchAddress("HLT_PFHT400_DiPFJetAve90_PFAlphaT0p51", &HLT_PFHT400_DiPFJetAve90_PFAlphaT0p51, &b_HLT_PFHT400_DiPFJetAve90_PFAlphaT0p51);
  fChain->SetBranchAddress("HLT_PFHT400_DiPFJetAve90_PFAlphaT0p52", &HLT_PFHT400_DiPFJetAve90_PFAlphaT0p52, &b_HLT_PFHT400_DiPFJetAve90_PFAlphaT0p52);
  fChain->SetBranchAddress("HLT_MET60_IsoTrk35_Loose", &HLT_MET60_IsoTrk35_Loose, &b_HLT_MET60_IsoTrk35_Loose);
  fChain->SetBranchAddress("HLT_MET75_IsoTrk50", &HLT_MET75_IsoTrk50, &b_HLT_MET75_IsoTrk50);
  fChain->SetBranchAddress("HLT_MET90_IsoTrk50", &HLT_MET90_IsoTrk50, &b_HLT_MET90_IsoTrk50);
  fChain->SetBranchAddress("HLT_PFMET120_BTagCSV_p067", &HLT_PFMET120_BTagCSV_p067, &b_HLT_PFMET120_BTagCSV_p067);
  fChain->SetBranchAddress("HLT_PFMET120_Mu5", &HLT_PFMET120_Mu5, &b_HLT_PFMET120_Mu5);
  fChain->SetBranchAddress("HLT_PFMET170_NotCleaned", &HLT_PFMET170_NotCleaned, &b_HLT_PFMET170_NotCleaned);
  fChain->SetBranchAddress("HLT_PFMET170_NoiseCleaned", &HLT_PFMET170_NoiseCleaned, &b_HLT_PFMET170_NoiseCleaned);
  fChain->SetBranchAddress("HLT_PFMET170_HBHECleaned", &HLT_PFMET170_HBHECleaned, &b_HLT_PFMET170_HBHECleaned);
  fChain->SetBranchAddress("HLT_PFMET170_JetIdCleaned", &HLT_PFMET170_JetIdCleaned, &b_HLT_PFMET170_JetIdCleaned);
  fChain->SetBranchAddress("HLT_PFMET170_BeamHaloCleaned", &HLT_PFMET170_BeamHaloCleaned, &b_HLT_PFMET170_BeamHaloCleaned);
  fChain->SetBranchAddress("HLT_PFMET170_HBHE_BeamHaloCleaned", &HLT_PFMET170_HBHE_BeamHaloCleaned, &b_HLT_PFMET170_HBHE_BeamHaloCleaned);
  fChain->SetBranchAddress("HLT_PFMETTypeOne190_HBHE_BeamHaloCleaned", &HLT_PFMETTypeOne190_HBHE_BeamHaloCleaned, &b_HLT_PFMETTypeOne190_HBHE_BeamHaloCleaned);
  fChain->SetBranchAddress("HLT_PFMET90_PFMHT90_IDTight", &HLT_PFMET90_PFMHT90_IDTight, &b_HLT_PFMET90_PFMHT90_IDTight);
  fChain->SetBranchAddress("HLT_PFMET100_PFMHT100_IDTight", &HLT_PFMET100_PFMHT100_IDTight, &b_HLT_PFMET100_PFMHT100_IDTight);
  fChain->SetBranchAddress("HLT_PFMET100_PFMHT100_IDTight_BeamHaloCleaned", &HLT_PFMET100_PFMHT100_IDTight_BeamHaloCleaned, &b_HLT_PFMET100_PFMHT100_IDTight_BeamHaloCleaned);
  fChain->SetBranchAddress("HLT_PFMET110_PFMHT110_IDTight", &HLT_PFMET110_PFMHT110_IDTight, &b_HLT_PFMET110_PFMHT110_IDTight);
  fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight", &HLT_PFMET120_PFMHT120_IDTight, &b_HLT_PFMET120_PFMHT120_IDTight);
  fChain->SetBranchAddress("HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight_BTagCSV_p067", &HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight_BTagCSV_p067, &b_HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight_BTagCSV_p067);
  fChain->SetBranchAddress("HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight", &HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight, &b_HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight);
  fChain->SetBranchAddress("HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq200", &HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq200, &b_HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq200);
  fChain->SetBranchAddress("HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq460", &HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq460, &b_HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq460);
  fChain->SetBranchAddress("HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240", &HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240, &b_HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240);
  fChain->SetBranchAddress("HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq500", &HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq500, &b_HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq500);
  fChain->SetBranchAddress("HLT_QuadPFJet_VBF", &HLT_QuadPFJet_VBF, &b_HLT_QuadPFJet_VBF);
  fChain->SetBranchAddress("HLT_L1_TripleJet_VBF", &HLT_L1_TripleJet_VBF, &b_HLT_L1_TripleJet_VBF);
  fChain->SetBranchAddress("HLT_QuadJet45_TripleBTagCSV_p087", &HLT_QuadJet45_TripleBTagCSV_p087, &b_HLT_QuadJet45_TripleBTagCSV_p087);
  fChain->SetBranchAddress("HLT_QuadJet45_DoubleBTagCSV_p087", &HLT_QuadJet45_DoubleBTagCSV_p087, &b_HLT_QuadJet45_DoubleBTagCSV_p087);
  fChain->SetBranchAddress("HLT_DoubleJet90_Double30_TripleBTagCSV_p087", &HLT_DoubleJet90_Double30_TripleBTagCSV_p087, &b_HLT_DoubleJet90_Double30_TripleBTagCSV_p087);
  fChain->SetBranchAddress("HLT_DoubleJet90_Double30_DoubleBTagCSV_p087", &HLT_DoubleJet90_Double30_DoubleBTagCSV_p087, &b_HLT_DoubleJet90_Double30_DoubleBTagCSV_p087);
  fChain->SetBranchAddress("HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160", &HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160, &b_HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160);
  fChain->SetBranchAddress("HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6", &HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6, &b_HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6);
  fChain->SetBranchAddress("HLT_DoubleJetsC112_DoubleBTagCSV_p026_DoublePFJetsC172", &HLT_DoubleJetsC112_DoubleBTagCSV_p026_DoublePFJetsC172, &b_HLT_DoubleJetsC112_DoubleBTagCSV_p026_DoublePFJetsC172);
  fChain->SetBranchAddress("HLT_DoubleJetsC112_DoubleBTagCSV_p014_DoublePFJetsC112MaxDeta1p6", &HLT_DoubleJetsC112_DoubleBTagCSV_p014_DoublePFJetsC112MaxDeta1p6, &b_HLT_DoubleJetsC112_DoubleBTagCSV_p014_DoublePFJetsC112MaxDeta1p6);
  fChain->SetBranchAddress("HLT_DoubleJetsC100_SingleBTagCSV_p026", &HLT_DoubleJetsC100_SingleBTagCSV_p026, &b_HLT_DoubleJetsC100_SingleBTagCSV_p026);
  fChain->SetBranchAddress("HLT_DoubleJetsC100_SingleBTagCSV_p014", &HLT_DoubleJetsC100_SingleBTagCSV_p014, &b_HLT_DoubleJetsC100_SingleBTagCSV_p014);
  fChain->SetBranchAddress("HLT_DoubleJetsC100_SingleBTagCSV_p026_SinglePFJetC350", &HLT_DoubleJetsC100_SingleBTagCSV_p026_SinglePFJetC350, &b_HLT_DoubleJetsC100_SingleBTagCSV_p026_SinglePFJetC350);
  fChain->SetBranchAddress("HLT_DoubleJetsC100_SingleBTagCSV_p014_SinglePFJetC350", &HLT_DoubleJetsC100_SingleBTagCSV_p014_SinglePFJetC350, &b_HLT_DoubleJetsC100_SingleBTagCSV_p014_SinglePFJetC350);
  fChain->SetBranchAddress("HLT_Photon135_PFMET100", &HLT_Photon135_PFMET100, &b_HLT_Photon135_PFMET100);
  fChain->SetBranchAddress("HLT_Photon20_CaloIdVL_IsoL", &HLT_Photon20_CaloIdVL_IsoL, &b_HLT_Photon20_CaloIdVL_IsoL);
  fChain->SetBranchAddress("HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_PFMET40", &HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_PFMET40, &b_HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_PFMET40);
  fChain->SetBranchAddress("HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_VBF", &HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_VBF, &b_HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_VBF);
  fChain->SetBranchAddress("HLT_Photon250_NoHE", &HLT_Photon250_NoHE, &b_HLT_Photon250_NoHE);
  fChain->SetBranchAddress("HLT_Photon300_NoHE", &HLT_Photon300_NoHE, &b_HLT_Photon300_NoHE);
  fChain->SetBranchAddress("HLT_Photon26_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon16_AND_HE10_R9Id65_Eta2_Mass60", &HLT_Photon26_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon16_AND_HE10_R9Id65_Eta2_Mass60, &b_HLT_Photon26_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon16_AND_HE10_R9Id65_Eta2_Mass60);
  fChain->SetBranchAddress("HLT_Photon36_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon22_AND_HE10_R9Id65_Eta2_Mass15", &HLT_Photon36_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon22_AND_HE10_R9Id65_Eta2_Mass15, &b_HLT_Photon36_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon22_AND_HE10_R9Id65_Eta2_Mass15);
  fChain->SetBranchAddress("HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_PFMET40", &HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_PFMET40, &b_HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_PFMET40);
  fChain->SetBranchAddress("HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_VBF", &HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_VBF, &b_HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_VBF);
  fChain->SetBranchAddress("HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_PFMET40", &HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_PFMET40, &b_HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_PFMET40);
  fChain->SetBranchAddress("HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_VBF", &HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_VBF, &b_HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_VBF);
  fChain->SetBranchAddress("HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_PFMET40", &HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_PFMET40, &b_HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_PFMET40);
  fChain->SetBranchAddress("HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF", &HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF, &b_HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF);
  fChain->SetBranchAddress("HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_PFMET40", &HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_PFMET40, &b_HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_PFMET40);
  fChain->SetBranchAddress("HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_VBF", &HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_VBF, &b_HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_VBF);
  fChain->SetBranchAddress("HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_PFMET40", &HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_PFMET40, &b_HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_PFMET40);
  fChain->SetBranchAddress("HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_VBF", &HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_VBF, &b_HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_VBF);
  fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL", &HLT_Mu8_TrkIsoVVL, &b_HLT_Mu8_TrkIsoVVL);
  fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL, &b_HLT_Mu17_TrkIsoVVL);
  fChain->SetBranchAddress("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30);
  fChain->SetBranchAddress("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30);
  fChain->SetBranchAddress("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30);
  fChain->SetBranchAddress("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30);
  fChain->SetBranchAddress("HLT_BTagMu_DiJet20_Mu5", &HLT_BTagMu_DiJet20_Mu5, &b_HLT_BTagMu_DiJet20_Mu5);
  fChain->SetBranchAddress("HLT_BTagMu_DiJet40_Mu5", &HLT_BTagMu_DiJet40_Mu5, &b_HLT_BTagMu_DiJet40_Mu5);
  fChain->SetBranchAddress("HLT_BTagMu_DiJet70_Mu5", &HLT_BTagMu_DiJet70_Mu5, &b_HLT_BTagMu_DiJet70_Mu5);
  fChain->SetBranchAddress("HLT_BTagMu_DiJet110_Mu5", &HLT_BTagMu_DiJet110_Mu5, &b_HLT_BTagMu_DiJet110_Mu5);
  fChain->SetBranchAddress("HLT_BTagMu_DiJet170_Mu5", &HLT_BTagMu_DiJet170_Mu5, &b_HLT_BTagMu_DiJet170_Mu5);
  fChain->SetBranchAddress("HLT_BTagMu_Jet300_Mu5", &HLT_BTagMu_Jet300_Mu5, &b_HLT_BTagMu_Jet300_Mu5);
  fChain->SetBranchAddress("HLT_BTagMu_AK8Jet300_Mu5", &HLT_BTagMu_AK8Jet300_Mu5, &b_HLT_BTagMu_AK8Jet300_Mu5);
  fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
  fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_L1JetTauSeeded", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_L1JetTauSeeded, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_L1JetTauSeeded);
  fChain->SetBranchAddress("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
  fChain->SetBranchAddress("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL", &HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL, &b_HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL);
  fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL", &HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL, &b_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL);
  fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL);
  fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ);
  fChain->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL);
  fChain->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ);
  fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL, &b_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL);
  fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL", &HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL, &b_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL);
  fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ);
  fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL);
  fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
  fChain->SetBranchAddress("HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL", &HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL, &b_HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL);
  fChain->SetBranchAddress("HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL", &HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL, &b_HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL);
  fChain->SetBranchAddress("HLT_Mu37_Ele27_CaloIdL_GsfTrkIdVL", &HLT_Mu37_Ele27_CaloIdL_GsfTrkIdVL, &b_HLT_Mu37_Ele27_CaloIdL_GsfTrkIdVL);
  fChain->SetBranchAddress("HLT_Mu27_Ele37_CaloIdL_GsfTrkIdVL", &HLT_Mu27_Ele37_CaloIdL_GsfTrkIdVL, &b_HLT_Mu27_Ele37_CaloIdL_GsfTrkIdVL);
  fChain->SetBranchAddress("HLT_Mu8_DiEle12_CaloIdL_TrackIdL", &HLT_Mu8_DiEle12_CaloIdL_TrackIdL, &b_HLT_Mu8_DiEle12_CaloIdL_TrackIdL);
  fChain->SetBranchAddress("HLT_Mu12_Photon25_CaloIdL", &HLT_Mu12_Photon25_CaloIdL, &b_HLT_Mu12_Photon25_CaloIdL);
  fChain->SetBranchAddress("HLT_Mu12_Photon25_CaloIdL_L1ISO", &HLT_Mu12_Photon25_CaloIdL_L1ISO, &b_HLT_Mu12_Photon25_CaloIdL_L1ISO);
  fChain->SetBranchAddress("HLT_Mu12_Photon25_CaloIdL_L1OR", &HLT_Mu12_Photon25_CaloIdL_L1OR, &b_HLT_Mu12_Photon25_CaloIdL_L1OR);
  fChain->SetBranchAddress("HLT_Mu17_Photon22_CaloIdL_L1ISO", &HLT_Mu17_Photon22_CaloIdL_L1ISO, &b_HLT_Mu17_Photon22_CaloIdL_L1ISO);
  fChain->SetBranchAddress("HLT_Mu17_Photon30_CaloIdL_L1ISO", &HLT_Mu17_Photon30_CaloIdL_L1ISO, &b_HLT_Mu17_Photon30_CaloIdL_L1ISO);
  fChain->SetBranchAddress("HLT_Mu17_Photon35_CaloIdL_L1ISO", &HLT_Mu17_Photon35_CaloIdL_L1ISO, &b_HLT_Mu17_Photon35_CaloIdL_L1ISO);
  fChain->SetBranchAddress("HLT_DiMu9_Ele9_CaloIdL_TrackIdL", &HLT_DiMu9_Ele9_CaloIdL_TrackIdL, &b_HLT_DiMu9_Ele9_CaloIdL_TrackIdL);
  fChain->SetBranchAddress("HLT_TripleMu_5_3_3", &HLT_TripleMu_5_3_3, &b_HLT_TripleMu_5_3_3);
  fChain->SetBranchAddress("HLT_TripleMu_12_10_5", &HLT_TripleMu_12_10_5, &b_HLT_TripleMu_12_10_5);
  fChain->SetBranchAddress("HLT_Mu3er_PFHT140_PFMET125", &HLT_Mu3er_PFHT140_PFMET125, &b_HLT_Mu3er_PFHT140_PFMET125);
  fChain->SetBranchAddress("HLT_Mu6_PFHT200_PFMET80_BTagCSV_p067", &HLT_Mu6_PFHT200_PFMET80_BTagCSV_p067, &b_HLT_Mu6_PFHT200_PFMET80_BTagCSV_p067);
  fChain->SetBranchAddress("HLT_Mu6_PFHT200_PFMET100", &HLT_Mu6_PFHT200_PFMET100, &b_HLT_Mu6_PFHT200_PFMET100);
  fChain->SetBranchAddress("HLT_Mu14er_PFMET100", &HLT_Mu14er_PFMET100, &b_HLT_Mu14er_PFMET100);
  fChain->SetBranchAddress("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL, &b_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL);
  fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL);
  fChain->SetBranchAddress("HLT_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Ele12_CaloIdL_TrackIdL_IsoVL, &b_HLT_Ele12_CaloIdL_TrackIdL_IsoVL);
  fChain->SetBranchAddress("HLT_Ele17_CaloIdL_GsfTrkIdVL", &HLT_Ele17_CaloIdL_GsfTrkIdVL, &b_HLT_Ele17_CaloIdL_GsfTrkIdVL);
  fChain->SetBranchAddress("HLT_Ele17_CaloIdL_TrackIdL_IsoVL", &HLT_Ele17_CaloIdL_TrackIdL_IsoVL, &b_HLT_Ele17_CaloIdL_TrackIdL_IsoVL);
  fChain->SetBranchAddress("HLT_Ele23_CaloIdL_TrackIdL_IsoVL", &HLT_Ele23_CaloIdL_TrackIdL_IsoVL, &b_HLT_Ele23_CaloIdL_TrackIdL_IsoVL);
  fChain->SetBranchAddress("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5", &HLT_PFHT650_WideJetMJJ900DEtaJJ1p5, &b_HLT_PFHT650_WideJetMJJ900DEtaJJ1p5);
  fChain->SetBranchAddress("HLT_PFHT650_WideJetMJJ950DEtaJJ1p5", &HLT_PFHT650_WideJetMJJ950DEtaJJ1p5, &b_HLT_PFHT650_WideJetMJJ950DEtaJJ1p5);
  fChain->SetBranchAddress("HLT_Photon22", &HLT_Photon22, &b_HLT_Photon22);
  fChain->SetBranchAddress("HLT_Photon30", &HLT_Photon30, &b_HLT_Photon30);
  fChain->SetBranchAddress("HLT_Photon36", &HLT_Photon36, &b_HLT_Photon36);
  fChain->SetBranchAddress("HLT_Photon50", &HLT_Photon50, &b_HLT_Photon50);
  fChain->SetBranchAddress("HLT_Photon75", &HLT_Photon75, &b_HLT_Photon75);
  fChain->SetBranchAddress("HLT_Photon90", &HLT_Photon90, &b_HLT_Photon90);
  fChain->SetBranchAddress("HLT_Photon120", &HLT_Photon120, &b_HLT_Photon120);
  fChain->SetBranchAddress("HLT_Photon175", &HLT_Photon175, &b_HLT_Photon175);
  fChain->SetBranchAddress("HLT_Photon165_HE10", &HLT_Photon165_HE10, &b_HLT_Photon165_HE10);
  fChain->SetBranchAddress("HLT_Photon22_R9Id90_HE10_IsoM", &HLT_Photon22_R9Id90_HE10_IsoM, &b_HLT_Photon22_R9Id90_HE10_IsoM);
  fChain->SetBranchAddress("HLT_Photon30_R9Id90_HE10_IsoM", &HLT_Photon30_R9Id90_HE10_IsoM, &b_HLT_Photon30_R9Id90_HE10_IsoM);
  fChain->SetBranchAddress("HLT_Photon36_R9Id90_HE10_IsoM", &HLT_Photon36_R9Id90_HE10_IsoM, &b_HLT_Photon36_R9Id90_HE10_IsoM);
  fChain->SetBranchAddress("HLT_Photon50_R9Id90_HE10_IsoM", &HLT_Photon50_R9Id90_HE10_IsoM, &b_HLT_Photon50_R9Id90_HE10_IsoM);
  fChain->SetBranchAddress("HLT_Photon75_R9Id90_HE10_IsoM", &HLT_Photon75_R9Id90_HE10_IsoM, &b_HLT_Photon75_R9Id90_HE10_IsoM);
  fChain->SetBranchAddress("HLT_Photon90_R9Id90_HE10_IsoM", &HLT_Photon90_R9Id90_HE10_IsoM, &b_HLT_Photon90_R9Id90_HE10_IsoM);
  fChain->SetBranchAddress("HLT_Photon120_R9Id90_HE10_IsoM", &HLT_Photon120_R9Id90_HE10_IsoM, &b_HLT_Photon120_R9Id90_HE10_IsoM);
  fChain->SetBranchAddress("HLT_Photon165_R9Id90_HE10_IsoM", &HLT_Photon165_R9Id90_HE10_IsoM, &b_HLT_Photon165_R9Id90_HE10_IsoM);
  fChain->SetBranchAddress("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90", &HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90, &b_HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90);
  fChain->SetBranchAddress("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70", &HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70, &b_HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70);
  fChain->SetBranchAddress("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55", &HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55, &b_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55);
  fChain->SetBranchAddress("HLT_Diphoton30_18_Solid_R9Id_AND_IsoCaloId_AND_HE_R9Id_Mass55", &HLT_Diphoton30_18_Solid_R9Id_AND_IsoCaloId_AND_HE_R9Id_Mass55, &b_HLT_Diphoton30_18_Solid_R9Id_AND_IsoCaloId_AND_HE_R9Id_Mass55);
  fChain->SetBranchAddress("HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55", &HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55, &b_HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55);
  fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_Muon", &HLT_Dimuon0_Jpsi_Muon, &b_HLT_Dimuon0_Jpsi_Muon);
  fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon", &HLT_Dimuon0_Upsilon_Muon, &b_HLT_Dimuon0_Upsilon_Muon);
  fChain->SetBranchAddress("HLT_QuadMuon0_Dimuon0_Jpsi", &HLT_QuadMuon0_Dimuon0_Jpsi, &b_HLT_QuadMuon0_Dimuon0_Jpsi);
  fChain->SetBranchAddress("HLT_QuadMuon0_Dimuon0_Upsilon", &HLT_QuadMuon0_Dimuon0_Upsilon, &b_HLT_QuadMuon0_Dimuon0_Upsilon);
  fChain->SetBranchAddress("HLT_Rsq0p25_Calo", &HLT_Rsq0p25_Calo, &b_HLT_Rsq0p25_Calo);
  fChain->SetBranchAddress("HLT_RsqMR240_Rsq0p09_MR200_4jet_Calo", &HLT_RsqMR240_Rsq0p09_MR200_4jet_Calo, &b_HLT_RsqMR240_Rsq0p09_MR200_4jet_Calo);
  fChain->SetBranchAddress("HLT_RsqMR240_Rsq0p09_MR200_Calo", &HLT_RsqMR240_Rsq0p09_MR200_Calo, &b_HLT_RsqMR240_Rsq0p09_MR200_Calo);
  fChain->SetBranchAddress("HLT_Rsq0p25", &HLT_Rsq0p25, &b_HLT_Rsq0p25);
  fChain->SetBranchAddress("HLT_Rsq0p30", &HLT_Rsq0p30, &b_HLT_Rsq0p30);
  fChain->SetBranchAddress("HLT_RsqMR240_Rsq0p09_MR200", &HLT_RsqMR240_Rsq0p09_MR200, &b_HLT_RsqMR240_Rsq0p09_MR200);
  fChain->SetBranchAddress("HLT_RsqMR240_Rsq0p09_MR200_4jet", &HLT_RsqMR240_Rsq0p09_MR200_4jet, &b_HLT_RsqMR240_Rsq0p09_MR200_4jet);
  fChain->SetBranchAddress("HLT_RsqMR270_Rsq0p09_MR200", &HLT_RsqMR270_Rsq0p09_MR200, &b_HLT_RsqMR270_Rsq0p09_MR200);
  fChain->SetBranchAddress("HLT_RsqMR270_Rsq0p09_MR200_4jet", &HLT_RsqMR270_Rsq0p09_MR200_4jet, &b_HLT_RsqMR270_Rsq0p09_MR200_4jet);
  fChain->SetBranchAddress("HLT_Rsq0p02_MR300_TriPFJet80_60_40_BTagCSV_p063_p20_Mbb60_200", &HLT_Rsq0p02_MR300_TriPFJet80_60_40_BTagCSV_p063_p20_Mbb60_200, &b_HLT_Rsq0p02_MR300_TriPFJet80_60_40_BTagCSV_p063_p20_Mbb60_200);
  fChain->SetBranchAddress("HLT_Rsq0p02_MR400_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200", &HLT_Rsq0p02_MR400_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200, &b_HLT_Rsq0p02_MR400_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200);
  fChain->SetBranchAddress("HLT_Rsq0p02_MR450_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200", &HLT_Rsq0p02_MR450_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200, &b_HLT_Rsq0p02_MR450_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200);
  fChain->SetBranchAddress("HLT_Rsq0p02_MR500_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200", &HLT_Rsq0p02_MR500_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200, &b_HLT_Rsq0p02_MR500_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200);
  fChain->SetBranchAddress("HLT_Rsq0p02_MR550_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200", &HLT_Rsq0p02_MR550_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200, &b_HLT_Rsq0p02_MR550_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200);
  fChain->SetBranchAddress("HLT_HT200_DisplacedDijet40_DisplacedTrack", &HLT_HT200_DisplacedDijet40_DisplacedTrack, &b_HLT_HT200_DisplacedDijet40_DisplacedTrack);
  fChain->SetBranchAddress("HLT_HT250_DisplacedDijet40_DisplacedTrack", &HLT_HT250_DisplacedDijet40_DisplacedTrack, &b_HLT_HT250_DisplacedDijet40_DisplacedTrack);
  fChain->SetBranchAddress("HLT_HT350_DisplacedDijet40_DisplacedTrack", &HLT_HT350_DisplacedDijet40_DisplacedTrack, &b_HLT_HT350_DisplacedDijet40_DisplacedTrack);
  fChain->SetBranchAddress("HLT_HT350_DisplacedDijet80_DisplacedTrack", &HLT_HT350_DisplacedDijet80_DisplacedTrack, &b_HLT_HT350_DisplacedDijet80_DisplacedTrack);
  fChain->SetBranchAddress("HLT_HT350_DisplacedDijet80_Tight_DisplacedTrack", &HLT_HT350_DisplacedDijet80_Tight_DisplacedTrack, &b_HLT_HT350_DisplacedDijet80_Tight_DisplacedTrack);
  fChain->SetBranchAddress("HLT_HT350_DisplacedDijet40_Inclusive", &HLT_HT350_DisplacedDijet40_Inclusive, &b_HLT_HT350_DisplacedDijet40_Inclusive);
  fChain->SetBranchAddress("HLT_HT400_DisplacedDijet40_Inclusive", &HLT_HT400_DisplacedDijet40_Inclusive, &b_HLT_HT400_DisplacedDijet40_Inclusive);
  fChain->SetBranchAddress("HLT_HT500_DisplacedDijet40_Inclusive", &HLT_HT500_DisplacedDijet40_Inclusive, &b_HLT_HT500_DisplacedDijet40_Inclusive);
  fChain->SetBranchAddress("HLT_HT550_DisplacedDijet40_Inclusive", &HLT_HT550_DisplacedDijet40_Inclusive, &b_HLT_HT550_DisplacedDijet40_Inclusive);
  fChain->SetBranchAddress("HLT_HT550_DisplacedDijet80_Inclusive", &HLT_HT550_DisplacedDijet80_Inclusive, &b_HLT_HT550_DisplacedDijet80_Inclusive);
  fChain->SetBranchAddress("HLT_HT650_DisplacedDijet80_Inclusive", &HLT_HT650_DisplacedDijet80_Inclusive, &b_HLT_HT650_DisplacedDijet80_Inclusive);
  fChain->SetBranchAddress("HLT_HT750_DisplacedDijet80_Inclusive", &HLT_HT750_DisplacedDijet80_Inclusive, &b_HLT_HT750_DisplacedDijet80_Inclusive);
  fChain->SetBranchAddress("HLT_VBF_DisplacedJet40_DisplacedTrack", &HLT_VBF_DisplacedJet40_DisplacedTrack, &b_HLT_VBF_DisplacedJet40_DisplacedTrack);
  fChain->SetBranchAddress("HLT_VBF_DisplacedJet40_DisplacedTrack_2TrackIP2DSig5", &HLT_VBF_DisplacedJet40_DisplacedTrack_2TrackIP2DSig5, &b_HLT_VBF_DisplacedJet40_DisplacedTrack_2TrackIP2DSig5);
  fChain->SetBranchAddress("HLT_VBF_DisplacedJet40_TightID_DisplacedTrack", &HLT_VBF_DisplacedJet40_TightID_DisplacedTrack, &b_HLT_VBF_DisplacedJet40_TightID_DisplacedTrack);
  fChain->SetBranchAddress("HLT_VBF_DisplacedJet40_Hadronic", &HLT_VBF_DisplacedJet40_Hadronic, &b_HLT_VBF_DisplacedJet40_Hadronic);
  fChain->SetBranchAddress("HLT_VBF_DisplacedJet40_Hadronic_2PromptTrack", &HLT_VBF_DisplacedJet40_Hadronic_2PromptTrack, &b_HLT_VBF_DisplacedJet40_Hadronic_2PromptTrack);
  fChain->SetBranchAddress("HLT_VBF_DisplacedJet40_TightID_Hadronic", &HLT_VBF_DisplacedJet40_TightID_Hadronic, &b_HLT_VBF_DisplacedJet40_TightID_Hadronic);
  fChain->SetBranchAddress("HLT_VBF_DisplacedJet40_VTightID_Hadronic", &HLT_VBF_DisplacedJet40_VTightID_Hadronic, &b_HLT_VBF_DisplacedJet40_VTightID_Hadronic);
  fChain->SetBranchAddress("HLT_VBF_DisplacedJet40_VVTightID_Hadronic", &HLT_VBF_DisplacedJet40_VVTightID_Hadronic, &b_HLT_VBF_DisplacedJet40_VVTightID_Hadronic);
  fChain->SetBranchAddress("HLT_VBF_DisplacedJet40_VTightID_DisplacedTrack", &HLT_VBF_DisplacedJet40_VTightID_DisplacedTrack, &b_HLT_VBF_DisplacedJet40_VTightID_DisplacedTrack);
  fChain->SetBranchAddress("HLT_VBF_DisplacedJet40_VVTightID_DisplacedTrack", &HLT_VBF_DisplacedJet40_VVTightID_DisplacedTrack, &b_HLT_VBF_DisplacedJet40_VVTightID_DisplacedTrack);
  fChain->SetBranchAddress("HLT_PFMETNoMu90_PFMHTNoMu90_IDTight", &HLT_PFMETNoMu90_PFMHTNoMu90_IDTight, &b_HLT_PFMETNoMu90_PFMHTNoMu90_IDTight);
  fChain->SetBranchAddress("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight", &HLT_PFMETNoMu100_PFMHTNoMu100_IDTight, &b_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight);
  fChain->SetBranchAddress("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight", &HLT_PFMETNoMu110_PFMHTNoMu110_IDTight, &b_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight);
  fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight);
  fChain->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_IDTight, &b_HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_IDTight);
  fChain->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu100_PFMHTNoMu100_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu100_PFMHTNoMu100_IDTight, &b_HLT_MonoCentralPFJet80_PFMETNoMu100_PFMHTNoMu100_IDTight);
  fChain->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight, &b_HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight);
  fChain->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight, &b_HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight);
  fChain->SetBranchAddress("HLT_Ele27_eta2p1_WPLoose_Gsf_HT200", &HLT_Ele27_eta2p1_WPLoose_Gsf_HT200, &b_HLT_Ele27_eta2p1_WPLoose_Gsf_HT200);
  fChain->SetBranchAddress("HLT_Photon90_CaloIdL_PFHT500", &HLT_Photon90_CaloIdL_PFHT500, &b_HLT_Photon90_CaloIdL_PFHT500);
  fChain->SetBranchAddress("HLT_DoubleMu8_Mass8_PFHT250", &HLT_DoubleMu8_Mass8_PFHT250, &b_HLT_DoubleMu8_Mass8_PFHT250);
  fChain->SetBranchAddress("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT250", &HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT250, &b_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT250);
  fChain->SetBranchAddress("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT250", &HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT250, &b_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT250);
  fChain->SetBranchAddress("HLT_DoubleMu8_Mass8_PFHT300", &HLT_DoubleMu8_Mass8_PFHT300, &b_HLT_DoubleMu8_Mass8_PFHT300);
  fChain->SetBranchAddress("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300", &HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300, &b_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300);
  fChain->SetBranchAddress("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300", &HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300, &b_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300);
  fChain->SetBranchAddress("HLT_Mu10_CentralPFJet30_BTagCSV_p13", &HLT_Mu10_CentralPFJet30_BTagCSV_p13, &b_HLT_Mu10_CentralPFJet30_BTagCSV_p13);
  fChain->SetBranchAddress("HLT_DoubleMu3_PFMET50", &HLT_DoubleMu3_PFMET50, &b_HLT_DoubleMu3_PFMET50);
  fChain->SetBranchAddress("HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13", &HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13, &b_HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13);
  fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_BTagCSV_p067_PFHT400", &HLT_Ele15_IsoVVVL_BTagCSV_p067_PFHT400, &b_HLT_Ele15_IsoVVVL_BTagCSV_p067_PFHT400);
  fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT350_PFMET50", &HLT_Ele15_IsoVVVL_PFHT350_PFMET50, &b_HLT_Ele15_IsoVVVL_PFHT350_PFMET50);
  fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT600", &HLT_Ele15_IsoVVVL_PFHT600, &b_HLT_Ele15_IsoVVVL_PFHT600);
  fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT350", &HLT_Ele15_IsoVVVL_PFHT350, &b_HLT_Ele15_IsoVVVL_PFHT350);
  fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT400_PFMET50", &HLT_Ele15_IsoVVVL_PFHT400_PFMET50, &b_HLT_Ele15_IsoVVVL_PFHT400_PFMET50);
  fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT400", &HLT_Ele15_IsoVVVL_PFHT400, &b_HLT_Ele15_IsoVVVL_PFHT400);
  fChain->SetBranchAddress("HLT_Ele50_IsoVVVL_PFHT400", &HLT_Ele50_IsoVVVL_PFHT400, &b_HLT_Ele50_IsoVVVL_PFHT400);
  fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60", &HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60, &b_HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60);
  fChain->SetBranchAddress("HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60", &HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60, &b_HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60);
  fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_BTagCSV_p067_PFHT400", &HLT_Mu15_IsoVVVL_BTagCSV_p067_PFHT400, &b_HLT_Mu15_IsoVVVL_BTagCSV_p067_PFHT400);
  fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT350_PFMET50", &HLT_Mu15_IsoVVVL_PFHT350_PFMET50, &b_HLT_Mu15_IsoVVVL_PFHT350_PFMET50);
  fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT600", &HLT_Mu15_IsoVVVL_PFHT600, &b_HLT_Mu15_IsoVVVL_PFHT600);
  fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT350", &HLT_Mu15_IsoVVVL_PFHT350, &b_HLT_Mu15_IsoVVVL_PFHT350);
  fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT400_PFMET50", &HLT_Mu15_IsoVVVL_PFHT400_PFMET50, &b_HLT_Mu15_IsoVVVL_PFHT400_PFMET50);
  fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT400", &HLT_Mu15_IsoVVVL_PFHT400, &b_HLT_Mu15_IsoVVVL_PFHT400);
  fChain->SetBranchAddress("HLT_Mu50_IsoVVVL_PFHT400", &HLT_Mu50_IsoVVVL_PFHT400, &b_HLT_Mu50_IsoVVVL_PFHT400);
  fChain->SetBranchAddress("HLT_Dimuon16_Jpsi", &HLT_Dimuon16_Jpsi, &b_HLT_Dimuon16_Jpsi);
  fChain->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel", &HLT_Dimuon10_Jpsi_Barrel, &b_HLT_Dimuon10_Jpsi_Barrel);
  fChain->SetBranchAddress("HLT_Dimuon8_PsiPrime_Barrel", &HLT_Dimuon8_PsiPrime_Barrel, &b_HLT_Dimuon8_PsiPrime_Barrel);
  fChain->SetBranchAddress("HLT_Dimuon8_Upsilon_Barrel", &HLT_Dimuon8_Upsilon_Barrel, &b_HLT_Dimuon8_Upsilon_Barrel);
  fChain->SetBranchAddress("HLT_Dimuon0_Phi_Barrel", &HLT_Dimuon0_Phi_Barrel, &b_HLT_Dimuon0_Phi_Barrel);
  fChain->SetBranchAddress("HLT_Mu16_TkMu0_dEta18_Onia", &HLT_Mu16_TkMu0_dEta18_Onia, &b_HLT_Mu16_TkMu0_dEta18_Onia);
  fChain->SetBranchAddress("HLT_Mu16_TkMu0_dEta18_Phi", &HLT_Mu16_TkMu0_dEta18_Phi, &b_HLT_Mu16_TkMu0_dEta18_Phi);
  fChain->SetBranchAddress("HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx", &HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx, &b_HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx);
  fChain->SetBranchAddress("HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx", &HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx, &b_HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx);
  fChain->SetBranchAddress("HLT_Mu8", &HLT_Mu8, &b_HLT_Mu8);
  fChain->SetBranchAddress("HLT_Mu17", &HLT_Mu17, &b_HLT_Mu17);
  fChain->SetBranchAddress("HLT_Mu3_PFJet40", &HLT_Mu3_PFJet40, &b_HLT_Mu3_PFJet40);
  fChain->SetBranchAddress("HLT_Ele8_CaloIdM_TrackIdM_PFJet30", &HLT_Ele8_CaloIdM_TrackIdM_PFJet30, &b_HLT_Ele8_CaloIdM_TrackIdM_PFJet30);
  fChain->SetBranchAddress("HLT_Ele12_CaloIdM_TrackIdM_PFJet30", &HLT_Ele12_CaloIdM_TrackIdM_PFJet30, &b_HLT_Ele12_CaloIdM_TrackIdM_PFJet30);
  fChain->SetBranchAddress("HLT_Ele17_CaloIdM_TrackIdM_PFJet30", &HLT_Ele17_CaloIdM_TrackIdM_PFJet30, &b_HLT_Ele17_CaloIdM_TrackIdM_PFJet30);
  fChain->SetBranchAddress("HLT_Ele23_CaloIdM_TrackIdM_PFJet30", &HLT_Ele23_CaloIdM_TrackIdM_PFJet30, &b_HLT_Ele23_CaloIdM_TrackIdM_PFJet30);
  fChain->SetBranchAddress("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet140", &HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet140, &b_HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet140);
  fChain->SetBranchAddress("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165", &HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165, &b_HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165);
  fChain->SetBranchAddress("HLT_PFHT400_SixJet30_DoubleBTagCSV_p056", &HLT_PFHT400_SixJet30_DoubleBTagCSV_p056, &b_HLT_PFHT400_SixJet30_DoubleBTagCSV_p056);
  fChain->SetBranchAddress("HLT_PFHT450_SixJet40_BTagCSV_p056", &HLT_PFHT450_SixJet40_BTagCSV_p056, &b_HLT_PFHT450_SixJet40_BTagCSV_p056);
  fChain->SetBranchAddress("HLT_PFHT400_SixJet30", &HLT_PFHT400_SixJet30, &b_HLT_PFHT400_SixJet30);
  fChain->SetBranchAddress("HLT_PFHT450_SixJet40", &HLT_PFHT450_SixJet40, &b_HLT_PFHT450_SixJet40);
  fChain->SetBranchAddress("HLT_Ele115_CaloIdVT_GsfTrkIdT", &HLT_Ele115_CaloIdVT_GsfTrkIdT, &b_HLT_Ele115_CaloIdVT_GsfTrkIdT);
  fChain->SetBranchAddress("HLT_Mu55", &HLT_Mu55, &b_HLT_Mu55);
  fChain->SetBranchAddress("HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15", &HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15, &b_HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15);
  fChain->SetBranchAddress("HLT_Photon90_CaloIdL_PFHT600", &HLT_Photon90_CaloIdL_PFHT600, &b_HLT_Photon90_CaloIdL_PFHT600);
  fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity60ForEndOfFill", &HLT_PixelTracks_Multiplicity60ForEndOfFill, &b_HLT_PixelTracks_Multiplicity60ForEndOfFill);
  fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity85ForEndOfFill", &HLT_PixelTracks_Multiplicity85ForEndOfFill, &b_HLT_PixelTracks_Multiplicity85ForEndOfFill);
  fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity110ForEndOfFill", &HLT_PixelTracks_Multiplicity110ForEndOfFill, &b_HLT_PixelTracks_Multiplicity110ForEndOfFill);
  fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity135ForEndOfFill", &HLT_PixelTracks_Multiplicity135ForEndOfFill, &b_HLT_PixelTracks_Multiplicity135ForEndOfFill);
  fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity160ForEndOfFill", &HLT_PixelTracks_Multiplicity160ForEndOfFill, &b_HLT_PixelTracks_Multiplicity160ForEndOfFill);
  fChain->SetBranchAddress("HLT_FullTracks_Multiplicity80", &HLT_FullTracks_Multiplicity80, &b_HLT_FullTracks_Multiplicity80);
  fChain->SetBranchAddress("HLT_FullTracks_Multiplicity100", &HLT_FullTracks_Multiplicity100, &b_HLT_FullTracks_Multiplicity100);
  fChain->SetBranchAddress("HLT_FullTracks_Multiplicity130", &HLT_FullTracks_Multiplicity130, &b_HLT_FullTracks_Multiplicity130);
  fChain->SetBranchAddress("HLT_FullTracks_Multiplicity150", &HLT_FullTracks_Multiplicity150, &b_HLT_FullTracks_Multiplicity150);
  fChain->SetBranchAddress("HLT_ECALHT800", &HLT_ECALHT800, &b_HLT_ECALHT800);
  fChain->SetBranchAddress("HLT_DiSC30_18_EIso_AND_HE_Mass70", &HLT_DiSC30_18_EIso_AND_HE_Mass70, &b_HLT_DiSC30_18_EIso_AND_HE_Mass70);
  fChain->SetBranchAddress("HLT_Photon125", &HLT_Photon125, &b_HLT_Photon125);
  fChain->SetBranchAddress("HLT_MET100", &HLT_MET100, &b_HLT_MET100);
  fChain->SetBranchAddress("HLT_MET150", &HLT_MET150, &b_HLT_MET150);
  fChain->SetBranchAddress("HLT_MET200", &HLT_MET200, &b_HLT_MET200);
  fChain->SetBranchAddress("HLT_Ele27_HighEta_Ele20_Mass55", &HLT_Ele27_HighEta_Ele20_Mass55, &b_HLT_Ele27_HighEta_Ele20_Mass55);
  fChain->SetBranchAddress("HLT_L1FatEvents", &HLT_L1FatEvents, &b_HLT_L1FatEvents);
  fChain->SetBranchAddress("HLT_Physics", &HLT_Physics, &b_HLT_Physics);
  fChain->SetBranchAddress("HLT_L1FatEvents_part0", &HLT_L1FatEvents_part0, &b_HLT_L1FatEvents_part0);
  fChain->SetBranchAddress("HLT_L1FatEvents_part1", &HLT_L1FatEvents_part1, &b_HLT_L1FatEvents_part1);
  fChain->SetBranchAddress("HLT_L1FatEvents_part2", &HLT_L1FatEvents_part2, &b_HLT_L1FatEvents_part2);
  fChain->SetBranchAddress("HLT_L1FatEvents_part3", &HLT_L1FatEvents_part3, &b_HLT_L1FatEvents_part3);
  fChain->SetBranchAddress("HLT_Random", &HLT_Random, &b_HLT_Random);
  fChain->SetBranchAddress("HLT_ZeroBias", &HLT_ZeroBias, &b_HLT_ZeroBias);
  fChain->SetBranchAddress("HLT_AK4CaloJet30", &HLT_AK4CaloJet30, &b_HLT_AK4CaloJet30);
  fChain->SetBranchAddress("HLT_AK4CaloJet40", &HLT_AK4CaloJet40, &b_HLT_AK4CaloJet40);
  fChain->SetBranchAddress("HLT_AK4CaloJet50", &HLT_AK4CaloJet50, &b_HLT_AK4CaloJet50);
  fChain->SetBranchAddress("HLT_AK4CaloJet80", &HLT_AK4CaloJet80, &b_HLT_AK4CaloJet80);
  fChain->SetBranchAddress("HLT_AK4CaloJet100", &HLT_AK4CaloJet100, &b_HLT_AK4CaloJet100);
  fChain->SetBranchAddress("HLT_AK4PFJet30", &HLT_AK4PFJet30, &b_HLT_AK4PFJet30);
  fChain->SetBranchAddress("HLT_AK4PFJet50", &HLT_AK4PFJet50, &b_HLT_AK4PFJet50);
  fChain->SetBranchAddress("HLT_AK4PFJet80", &HLT_AK4PFJet80, &b_HLT_AK4PFJet80);
  fChain->SetBranchAddress("HLT_AK4PFJet100", &HLT_AK4PFJet100, &b_HLT_AK4PFJet100);
  fChain->SetBranchAddress("HLT_HISinglePhoton10", &HLT_HISinglePhoton10, &b_HLT_HISinglePhoton10);
  fChain->SetBranchAddress("HLT_HISinglePhoton15", &HLT_HISinglePhoton15, &b_HLT_HISinglePhoton15);
  fChain->SetBranchAddress("HLT_HISinglePhoton20", &HLT_HISinglePhoton20, &b_HLT_HISinglePhoton20);
  fChain->SetBranchAddress("HLT_HISinglePhoton40", &HLT_HISinglePhoton40, &b_HLT_HISinglePhoton40);
  fChain->SetBranchAddress("HLT_HISinglePhoton60", &HLT_HISinglePhoton60, &b_HLT_HISinglePhoton60);
  fChain->SetBranchAddress("HLT_EcalCalibration", &HLT_EcalCalibration, &b_HLT_EcalCalibration);
  fChain->SetBranchAddress("HLT_HcalCalibration", &HLT_HcalCalibration, &b_HLT_HcalCalibration);
  fChain->SetBranchAddress("HLT_GlobalRunHPDNoise", &HLT_GlobalRunHPDNoise, &b_HLT_GlobalRunHPDNoise);
  fChain->SetBranchAddress("HLT_L1BptxMinus", &HLT_L1BptxMinus, &b_HLT_L1BptxMinus);
  fChain->SetBranchAddress("HLT_L1BptxPlus", &HLT_L1BptxPlus, &b_HLT_L1BptxPlus);
  fChain->SetBranchAddress("HLT_L1NotBptxOR", &HLT_L1NotBptxOR, &b_HLT_L1NotBptxOR);
  fChain->SetBranchAddress("HLT_L1BeamGasMinus", &HLT_L1BeamGasMinus, &b_HLT_L1BeamGasMinus);
  fChain->SetBranchAddress("HLT_L1BeamGasPlus", &HLT_L1BeamGasPlus, &b_HLT_L1BeamGasPlus);
  fChain->SetBranchAddress("HLT_L1BptxXOR", &HLT_L1BptxXOR, &b_HLT_L1BptxXOR);
  fChain->SetBranchAddress("HLT_L1MinimumBiasHF_OR", &HLT_L1MinimumBiasHF_OR, &b_HLT_L1MinimumBiasHF_OR);
  fChain->SetBranchAddress("HLT_L1MinimumBiasHF_AND", &HLT_L1MinimumBiasHF_AND, &b_HLT_L1MinimumBiasHF_AND);
  fChain->SetBranchAddress("HLT_HcalNZS", &HLT_HcalNZS, &b_HLT_HcalNZS);
  fChain->SetBranchAddress("HLT_HcalPhiSym", &HLT_HcalPhiSym, &b_HLT_HcalPhiSym);
  fChain->SetBranchAddress("HLT_HcalIsolatedbunch", &HLT_HcalIsolatedbunch, &b_HLT_HcalIsolatedbunch);
  fChain->SetBranchAddress("HLT_ZeroBias_FirstCollisionAfterAbortGap", &HLT_ZeroBias_FirstCollisionAfterAbortGap, &b_HLT_ZeroBias_FirstCollisionAfterAbortGap);
  fChain->SetBranchAddress("HLT_ZeroBias_FirstCollisionAfterAbortGap_copy", &HLT_ZeroBias_FirstCollisionAfterAbortGap_copy, &b_HLT_ZeroBias_FirstCollisionAfterAbortGap_copy);
  fChain->SetBranchAddress("HLT_ZeroBias_FirstCollisionAfterAbortGap_TCDS", &HLT_ZeroBias_FirstCollisionAfterAbortGap_TCDS, &b_HLT_ZeroBias_FirstCollisionAfterAbortGap_TCDS);
  fChain->SetBranchAddress("HLT_ZeroBias_IsolatedBunches", &HLT_ZeroBias_IsolatedBunches, &b_HLT_ZeroBias_IsolatedBunches);
  fChain->SetBranchAddress("HLT_ZeroBias_FirstCollisionInTrain", &HLT_ZeroBias_FirstCollisionInTrain, &b_HLT_ZeroBias_FirstCollisionInTrain);
  fChain->SetBranchAddress("HLT_ZeroBias_FirstBXAfterTrain", &HLT_ZeroBias_FirstBXAfterTrain, &b_HLT_ZeroBias_FirstBXAfterTrain);
  fChain->SetBranchAddress("HLT_Photon500", &HLT_Photon500, &b_HLT_Photon500);
  fChain->SetBranchAddress("HLT_Photon600", &HLT_Photon600, &b_HLT_Photon600);
  fChain->SetBranchAddress("HLT_Mu300", &HLT_Mu300, &b_HLT_Mu300);
  fChain->SetBranchAddress("HLT_Mu350", &HLT_Mu350, &b_HLT_Mu350);
  fChain->SetBranchAddress("HLT_MET250", &HLT_MET250, &b_HLT_MET250);
  fChain->SetBranchAddress("HLT_MET300", &HLT_MET300, &b_HLT_MET300);
  fChain->SetBranchAddress("HLT_MET600", &HLT_MET600, &b_HLT_MET600);
  fChain->SetBranchAddress("HLT_MET700", &HLT_MET700, &b_HLT_MET700);
  fChain->SetBranchAddress("HLT_PFMET300", &HLT_PFMET300, &b_HLT_PFMET300);
  fChain->SetBranchAddress("HLT_PFMET400", &HLT_PFMET400, &b_HLT_PFMET400);
  fChain->SetBranchAddress("HLT_PFMET500", &HLT_PFMET500, &b_HLT_PFMET500);
  fChain->SetBranchAddress("HLT_PFMET600", &HLT_PFMET600, &b_HLT_PFMET600);
  fChain->SetBranchAddress("HLT_Ele250_CaloIdVT_GsfTrkIdT", &HLT_Ele250_CaloIdVT_GsfTrkIdT, &b_HLT_Ele250_CaloIdVT_GsfTrkIdT);
  fChain->SetBranchAddress("HLT_Ele300_CaloIdVT_GsfTrkIdT", &HLT_Ele300_CaloIdVT_GsfTrkIdT, &b_HLT_Ele300_CaloIdVT_GsfTrkIdT);
  fChain->SetBranchAddress("HLT_HT2000", &HLT_HT2000, &b_HLT_HT2000);
  fChain->SetBranchAddress("HLT_HT2500", &HLT_HT2500, &b_HLT_HT2500);
  fChain->SetBranchAddress("HLT_IsoTrackHE", &HLT_IsoTrackHE, &b_HLT_IsoTrackHE);
  fChain->SetBranchAddress("HLT_IsoTrackHB", &HLT_IsoTrackHB, &b_HLT_IsoTrackHB);
  fChain->SetBranchAddress("HLTriggerFinalPath", &HLTriggerFinalPath, &b_HLTriggerFinalPath);
  fChain->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, &b_Flag_HBHENoiseFilter);
  fChain->SetBranchAddress("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter, &b_Flag_HBHENoiseIsoFilter);
  fChain->SetBranchAddress("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter, &b_Flag_CSCTightHaloFilter);
  fChain->SetBranchAddress("Flag_CSCTightHaloTrkMuUnvetoFilter", &Flag_CSCTightHaloTrkMuUnvetoFilter, &b_Flag_CSCTightHaloTrkMuUnvetoFilter);
  fChain->SetBranchAddress("Flag_CSCTightHalo2015Filter", &Flag_CSCTightHalo2015Filter, &b_Flag_CSCTightHalo2015Filter);
  fChain->SetBranchAddress("Flag_globalTightHalo2016Filter", &Flag_globalTightHalo2016Filter, &b_Flag_globalTightHalo2016Filter);
  fChain->SetBranchAddress("Flag_globalSuperTightHalo2016Filter", &Flag_globalSuperTightHalo2016Filter, &b_Flag_globalSuperTightHalo2016Filter);
  fChain->SetBranchAddress("Flag_HcalStripHaloFilter", &Flag_HcalStripHaloFilter, &b_Flag_HcalStripHaloFilter);
  fChain->SetBranchAddress("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter, &b_Flag_hcalLaserEventFilter);
  fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
  fChain->SetBranchAddress("Flag_EcalDeadCellBoundaryEnergyFilter", &Flag_EcalDeadCellBoundaryEnergyFilter, &b_Flag_EcalDeadCellBoundaryEnergyFilter);
  fChain->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices, &b_Flag_goodVertices);
  fChain->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
  fChain->SetBranchAddress("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter, &b_Flag_ecalLaserCorrFilter);
  fChain->SetBranchAddress("Flag_trkPOGFilters", &Flag_trkPOGFilters, &b_Flag_trkPOGFilters);
  fChain->SetBranchAddress("Flag_chargedHadronTrackResolutionFilter", &Flag_chargedHadronTrackResolutionFilter, &b_Flag_chargedHadronTrackResolutionFilter);
  fChain->SetBranchAddress("Flag_muonBadTrackFilter", &Flag_muonBadTrackFilter, &b_Flag_muonBadTrackFilter);
  fChain->SetBranchAddress("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X, &b_Flag_trkPOG_manystripclus53X);
  fChain->SetBranchAddress("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X, &b_Flag_trkPOG_toomanystripclus53X);
  fChain->SetBranchAddress("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters, &b_Flag_trkPOG_logErrorTooManyClusters);
  fChain->SetBranchAddress("Flag_METFilters", &Flag_METFilters, &b_Flag_METFilters);
  Notify();
}

Bool_t nanoAnalysis::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void nanoAnalysis::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t nanoAnalysis::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}

#endif // #ifdef nanoAnalysis_cxx
