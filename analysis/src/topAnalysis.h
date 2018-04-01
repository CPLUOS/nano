//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jan 18 21:36:33 2018 by ROOT version 6.10/09
// from TTree Events/Events
// found on file: /home/nanoAOD/run2_2016v3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/180117_180123/0000/nanoAOD_111.root
//////////////////////////////////////////////////////////

#ifndef topAnalysis_h
#define topAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include <TH1D.h>
#include <vector>

#include "RoccoR.cc"
#include "pileUpTool.h"
#include "MuonScaleFactorEvaluator.h"
#include "ElecScaleFactorEvaluator.h"
#include "json.h"
#include "lumiTool.h"
#include "BTagWeightEvaluator.h"
//#include "TopTriggerSF.h"
//#include "TTbarModeDefs.h"
// Header file for the classes stored in the TTree if any.

class topAnalysis {
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
  //float b_tri, b_tri_up, b_tri_dn;
  //Triggers
  Bool_t b_trig_m, b_trig_m2,  b_trig_e, b_trig_mm, b_trig_em, b_trig_ee;


  // Tools
  RoccoR* m_rocCor;
  TH1D* hist_mc;
  MuonScaleFactorEvaluator muonSF_;
  ElecScaleFactorEvaluator elecSF_;
  BTagWeightEvaluator m_btagSF;

  pileUpTool *m_pileUp;
  lumiTool* m_lumi;
  //LumiMap
  std::map<UInt_t, std::vector<std::array<UInt_t, 2>>> lumiMap;

  
  //Making output branch
  void MakeBranch(TTree* t);
  void resetBranch();

  void analysis();
  //For Selection
  Bool_t lumiCheck();
  enum TTLLChannel { CH_NOLL = 0, CH_MUEL, CH_ELEL, CH_MUMU }; 
  std::vector<TParticle> muonSelection();
  std::vector<TParticle> elecSelection();
  std::vector<TLorentzVector> recoleps;
  std::vector<TParticle> jetSelection();
  std::vector<TParticle> bjetSelection();
  Double_t roccoR(TLorentzVector m, int q, int nGen, int nTrackerLayers);

 public :
  //set output file
  void setOutput(std::string outputName);
  void LoadModules(pileUpTool* pileUp, lumiTool* lumi); 

  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  Bool_t m_isMC;
  Bool_t m_isDL;
  Bool_t m_isSL_e;
  Bool_t m_isSL_m;

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
  Float_t         genWeight;
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
  Int_t           Pileup_nPU;
  Float_t         Pileup_nTrueInt;
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
  Int_t           Jet_genJetIdx[35];   //[nJet]
  Int_t           Jet_hadronFlavour[35];   //[nJet]
  Int_t           Jet_partonFlavour[35];   //[nJet]
  Int_t           Muon_genPartIdx[10];   //[nMuon]
  UChar_t         Muon_genPartFlav[10];   //[nMuon]
  UChar_t         Electron_cleanmask[6];   //[nElectron]
  UChar_t         Jet_cleanmask[35];   //[nJet]
  UChar_t         Muon_cleanmask[10];   //[nMuon]
  Bool_t          HLT_Ele27_WPTight_Gsf;
  Bool_t          HLT_IsoMu24;
  Bool_t          HLT_IsoTkMu24;
  Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;
  Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;
  Bool_t          HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL;
  Bool_t          HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;
  Bool_t          HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
  Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;
  Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
  Bool_t          HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;
  Bool_t          HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;

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
  TBranch        *b_genWeight;   //!
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
  TBranch        *b_Pileup_nPU;   //!
  TBranch        *b_Pileup_nTrueInt;   //!
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
  TBranch        *b_Jet_genJetIdx;   //!
  TBranch        *b_Jet_hadronFlavour;   //!
  TBranch        *b_Jet_partonFlavour;   //!
  TBranch        *b_Muon_genPartIdx;   //!
  TBranch        *b_Muon_genPartFlav;   //!
  TBranch        *b_Electron_cleanmask;   //!
  TBranch        *b_Jet_cleanmask;   //!
  TBranch        *b_Muon_cleanmask;   //!
  TBranch        *b_HLT_Ele27_WPTight_Gsf;   //!
  TBranch        *b_HLT_IsoMu24;   //!
  TBranch        *b_HLT_IsoTkMu24;   //!
  TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;   //!
  TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;   //!
  TBranch        *b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL;   //!
  TBranch        *b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;   //!
  TBranch        *b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
  TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;   //!
  TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;   //!
  TBranch        *b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;   //!
  TBranch        *b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;   //!
      
      
  topAnalysis(TTree *tree=0, Bool_t flag = false, Bool_t dl = false, Bool_t sle = false, Bool_t slm = false);
  virtual ~topAnalysis();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef topAnalysis_cxx
topAnalysis::topAnalysis(TTree *tree, Bool_t flag, Bool_t dl, Bool_t sle, Bool_t slm) : fChain(0), m_isMC(flag), m_isDL(dl), m_isSL_e(sle), m_isSL_m(slm)
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

topAnalysis::~topAnalysis()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
  m_output->Write();
  m_output->Close();
}

Int_t topAnalysis::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t topAnalysis::LoadTree(Long64_t entry)
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

void topAnalysis::Init(TTree *tree)
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
  fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
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
  fChain->SetBranchAddress("Pileup_nPU", &Pileup_nPU, &b_Pileup_nPU);
  fChain->SetBranchAddress("Pileup_nTrueInt", &Pileup_nTrueInt, &b_Pileup_nTrueInt);
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
  fChain->SetBranchAddress("Jet_genJetIdx", Jet_genJetIdx, &b_Jet_genJetIdx);
  fChain->SetBranchAddress("Jet_hadronFlavour", Jet_hadronFlavour, &b_Jet_hadronFlavour);
  fChain->SetBranchAddress("Jet_partonFlavour", Jet_partonFlavour, &b_Jet_partonFlavour);
  fChain->SetBranchAddress("Muon_genPartIdx", Muon_genPartIdx, &b_Muon_genPartIdx);
  fChain->SetBranchAddress("Muon_genPartFlav", Muon_genPartFlav, &b_Muon_genPartFlav);
  fChain->SetBranchAddress("Electron_cleanmask", Electron_cleanmask, &b_Electron_cleanmask);
  fChain->SetBranchAddress("Jet_cleanmask", Jet_cleanmask, &b_Jet_cleanmask);
  fChain->SetBranchAddress("Muon_cleanmask", Muon_cleanmask, &b_Muon_cleanmask);
  fChain->SetBranchAddress("HLT_Ele27_WPTight_Gsf", &HLT_Ele27_WPTight_Gsf, &b_HLT_Ele27_WPTight_Gsf);
  fChain->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24, &b_HLT_IsoMu24);
  fChain->SetBranchAddress("HLT_IsoTkMu24", &HLT_IsoTkMu24, &b_HLT_IsoTkMu24);
  fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL);
  fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ);
  fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL, &b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL);
  fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ, &b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ);
  fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
  fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL);
  fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ);
  fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL);
  fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
  Notify();
}

Bool_t topAnalysis::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void topAnalysis::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t topAnalysis::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef topAnalysis_cxx
