//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu May 31 16:46:05 2018 by ROOT version 6.10/09
// from TTree HadTruthEvents/HadTruthEvents
// found on file: /home/iwatson/tsW/mcATnlo/data/tt01j_bbar_1lp_FxFx_AOD/HADTRUTHAOD/1.root
//////////////////////////////////////////////////////////

#ifndef HadTruthEvents_h
#define HadTruthEvents_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class HadTruthEvents {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          luminosityBlock;
   ULong64_t       event;
   UInt_t          nGenPart;
   Float_t         GenPart_eta[5000];   //[nGenPart]
   Float_t         GenPart_mass[5000];   //[nGenPart]
   Float_t         GenPart_phi[5000];   //[nGenPart]
   Float_t         GenPart_pt[5000];   //[nGenPart]
   Int_t           GenPart_genPartIdxMother[5000];   //[nGenPart]
   Int_t           GenPart_pdgId[5000];   //[nGenPart]
   Int_t           GenPart_status[5000];   //[nGenPart]
   Int_t           GenPart_statusFlags[5000];   //[nGenPart]
   UInt_t          nhad_jet;
   Float_t         had_jet_btagCMVA[250];   //[nhad_jet]
   Float_t         had_jet_btagCSVV2[250];   //[nhad_jet]
   Float_t         had_jet_btagDeepB[250];   //[nhad_jet]
   Float_t         had_jet_btagDeepC[250];   //[nhad_jet]
   Float_t         had_jet_eta[250];   //[nhad_jet]
   Float_t         had_jet_mass[250];   //[nhad_jet]
   Float_t         had_jet_phi[250];   //[nhad_jet]
   Float_t         had_jet_pt[250];   //[nhad_jet]
   UInt_t          nhad;
   Float_t         had_jetDR[250];   //[nhad]
   Float_t         had_legDR[250];   //[nhad]
   Float_t         had_diffMass[250];   //[nhad]
   Float_t         had_lxy[250];   //[nhad]
   Float_t         had_lxyErr[250];   //[nhad]
   Float_t         had_l3D[250];   //[nhad]
   Float_t         had_l3DErr[250];   //[nhad]
   Float_t         had_dca[250];   //[nhad]
   Float_t         had_angleXY[250];   //[nhad]
   Float_t         had_angleXYZ[250];   //[nhad]
   Float_t         had_dau1_chi2[250];   //[nhad]
   Float_t         had_dau1_nHits[250];   //[nhad]
   Float_t         had_dau1_pt[250];   //[nhad]
   Float_t         had_dau1_ipsigXY[250];   //[nhad]
   Float_t         had_dau1_ipsigZ[250];   //[nhad]
   Float_t         had_dau2_chi2[250];   //[nhad]
   Float_t         had_dau2_nHits[250];   //[nhad]
   Float_t         had_dau2_pt[250];   //[nhad]
   Float_t         had_dau2_ipsigXY[250];   //[nhad]
   Float_t         had_dau2_ipsigZ[250];   //[nhad]
   Int_t           had_nJet[250];   //[nhad]
   Int_t           had_nDau[250];   //[nhad]
   Int_t           had_dau1_charge[250];   //[nhad]
   Int_t           had_dau2_charge[250];   //[nhad]
   Int_t           had_idx[250];   //[nhad]
   Int_t           had_dau1_idx[250];   //[nhad]
   Int_t           had_dau2_idx[250];   //[nhad]
   UInt_t          ngenHadron;
   Float_t         genHadron_dau1_pt[50];   //[ngenHadron]
   Float_t         genHadron_dau2_pt[50];   //[ngenHadron]
   Float_t         genHadron_dau1_eta[50];   //[ngenHadron]
   Float_t         genHadron_dau2_eta[50];   //[ngenHadron]
   Float_t         genHadron_dau1_phi[50];   //[ngenHadron]
   Float_t         genHadron_dau2_phi[50];   //[ngenHadron]
   Float_t         genHadron_vx[50];   //[ngenHadron]
   Float_t         genHadron_vy[50];   //[ngenHadron]
   Float_t         genHadron_vz[50];   //[ngenHadron]
   Int_t           genHadron_isGenHadFromTsb[50];   //[ngenHadron]
   Int_t           genHadron_dau1_pdgId[50];   //[ngenHadron]
   Int_t           genHadron_dau2_pdgId[50];   //[ngenHadron]
   Int_t           genHadron_isGenParticle[50];   //[ngenHadron]
   UChar_t         genHadron_isGenHadFromTop[50];   //[ngenHadron]
   UChar_t         genHadron_inVol[50];   //[ngenHadron]
   UChar_t         genHadron_isMatching[50];   //[ngenHadron]
   UChar_t         genHadron_isMatched[50];   //[ngenHadron]
   UInt_t          nhadTruth;
   Int_t           hadTruth_nMatched[250];   //[nhadTruth]
   Int_t           hadTruth_nTrueDau[250];   //[nhadTruth]
   Int_t           hadTruth_isHadFromTsb[250];   //[nhadTruth]
   UChar_t         hadTruth_isHadFromTop[250];   //[nhadTruth]
   UChar_t         hadTruth_isHadFromW[250];   //[nhadTruth]
   UChar_t         hadTruth_isHadFromS[250];   //[nhadTruth]
   UChar_t         hadTruth_isHadFromC[250];   //[nhadTruth]
   UChar_t         hadTruth_isHadFromB[250];   //[nhadTruth]
   UInt_t          nOtherPV;
   Float_t         OtherPV_z[5];   //[nOtherPV]
   Float_t         PV_ndof;
   Float_t         PV_x;
   Float_t         PV_y;
   Float_t         PV_z;
   Float_t         PV_chi2;
   Float_t         PV_score;
   Int_t           PV_npvs;
   Int_t           PV_npvsGood;
   UInt_t          nSV;
   Float_t         SV_dlen[10];   //[nSV]
   Float_t         SV_dlenSig[10];   //[nSV]
   Float_t         SV_pAngle[10];   //[nSV]
   Float_t         had_chi2[250];   //[nhad]
   Float_t         had_eta[250];   //[nhad]
   Float_t         had_mass[250];   //[nhad]
   Float_t         had_phi[250];   //[nhad]
   Float_t         had_pt[250];   //[nhad]
   Float_t         had_x[250];   //[nhad]
   Float_t         had_y[250];   //[nhad]
   Float_t         had_z[250];   //[nhad]
   Int_t           had_ndof[250];   //[nhad]
   Int_t           had_pdgId[250];   //[nhad]
   Float_t         genHadron_eta[50];   //[ngenHadron]
   Float_t         genHadron_mass[50];   //[ngenHadron]
   Float_t         genHadron_phi[38];   //[ngenHadron]
   Float_t         genHadron_pt[38];   //[ngenHadron]
   Float_t         genHadron_x[38];   //[ngenHadron]
   Float_t         genHadron_y[38];   //[ngenHadron]
   Float_t         genHadron_z[38];   //[ngenHadron]
   Int_t           genHadron_pdgId[38];   //[ngenHadron]
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
   Bool_t          HLT_Photon250_PFMET100;
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
   Bool_t          HLT_PixelTracks_Multiplicity250ForEndOfFill;
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
   Bool_t          HLT_L1FatHadTruthEvents;
   Bool_t          HLT_Physics;
   Bool_t          HLT_L1FatHadTruthEvents_part0;
   Bool_t          HLT_L1FatHadTruthEvents_part1;
   Bool_t          HLT_L1FatHadTruthEvents_part2;
   Bool_t          HLT_L1FatHadTruthEvents_part3;
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
   Bool_t          L1Reco_step;
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
   Bool_t          Flag_BadChargedCandidateFilter;
   Bool_t          Flag_BadPFMuonFilter;
   Bool_t          Flag_BadChargedCandidateSummer16Filter;
   Bool_t          Flag_BadPFMuonSummer16Filter;
   Bool_t          Flag_trkPOG_manystripclus53X;
   Bool_t          Flag_trkPOG_toomanystripclus53X;
   Bool_t          Flag_trkPOG_logErrorTooManyClusters;
   Bool_t          Flag_METFilters;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_event;   //!
   TBranch        *b_nGenPart;   //!
   TBranch        *b_GenPart_eta;   //!
   TBranch        *b_GenPart_mass;   //!
   TBranch        *b_GenPart_phi;   //!
   TBranch        *b_GenPart_pt;   //!
   TBranch        *b_GenPart_genPartIdxMother;   //!
   TBranch        *b_GenPart_pdgId;   //!
   TBranch        *b_GenPart_status;   //!
   TBranch        *b_GenPart_statusFlags;   //!
   TBranch        *b_nhad_jet;   //!
   TBranch        *b_had_jet_btagCMVA;   //!
   TBranch        *b_had_jet_btagCSVV2;   //!
   TBranch        *b_had_jet_btagDeepB;   //!
   TBranch        *b_had_jet_btagDeepC;   //!
   TBranch        *b_had_jet_eta;   //!
   TBranch        *b_had_jet_mass;   //!
   TBranch        *b_had_jet_phi;   //!
   TBranch        *b_had_jet_pt;   //!
   TBranch        *b_nhad;   //!
   TBranch        *b_had_jetDR;   //!
   TBranch        *b_had_legDR;   //!
   TBranch        *b_had_diffMass;   //!
   TBranch        *b_had_lxy;   //!
   TBranch        *b_had_lxyErr;   //!
   TBranch        *b_had_l3D;   //!
   TBranch        *b_had_l3DErr;   //!
   TBranch        *b_had_dca;   //!
   TBranch        *b_had_angleXY;   //!
   TBranch        *b_had_angleXYZ;   //!
   TBranch        *b_had_dau1_chi2;   //!
   TBranch        *b_had_dau1_nHits;   //!
   TBranch        *b_had_dau1_pt;   //!
   TBranch        *b_had_dau1_ipsigXY;   //!
   TBranch        *b_had_dau1_ipsigZ;   //!
   TBranch        *b_had_dau2_chi2;   //!
   TBranch        *b_had_dau2_nHits;   //!
   TBranch        *b_had_dau2_pt;   //!
   TBranch        *b_had_dau2_ipsigXY;   //!
   TBranch        *b_had_dau2_ipsigZ;   //!
   TBranch        *b_had_nJet;   //!
   TBranch        *b_had_nDau;   //!
   TBranch        *b_had_dau1_charge;   //!
   TBranch        *b_had_dau2_charge;   //!
   TBranch        *b_had_idx;   //!
   TBranch        *b_had_dau1_idx;   //!
   TBranch        *b_had_dau2_idx;   //!
   TBranch        *b_ngenHadron;   //!
   TBranch        *b_genHadron_dau1_pt;   //!
   TBranch        *b_genHadron_dau2_pt;   //!
   TBranch        *b_genHadron_dau1_eta;   //!
   TBranch        *b_genHadron_dau2_eta;   //!
   TBranch        *b_genHadron_dau1_phi;   //!
   TBranch        *b_genHadron_dau2_phi;   //!
   TBranch        *b_genHadron_vx;   //!
   TBranch        *b_genHadron_vy;   //!
   TBranch        *b_genHadron_vz;   //!
   TBranch        *b_genHadron_isGenHadFromTsb;   //!
   TBranch        *b_genHadron_dau1_pdgId;   //!
   TBranch        *b_genHadron_dau2_pdgId;   //!
   TBranch        *b_genHadron_isGenParticle;   //!
   TBranch        *b_genHadron_isGenHadFromTop;   //!
   TBranch        *b_genHadron_inVol;   //!
   TBranch        *b_genHadron_isMatching;   //!
   TBranch        *b_genHadron_isMatched;   //!
   TBranch        *b_nhadTruth;   //!
   TBranch        *b_hadTruth_nMatched;   //!
   TBranch        *b_hadTruth_nTrueDau;   //!
   TBranch        *b_hadTruth_isHadFromTsb;   //!
   TBranch        *b_hadTruth_isHadFromTop;   //!
   TBranch        *b_hadTruth_isHadFromW;   //!
   TBranch        *b_hadTruth_isHadFromS;   //!
   TBranch        *b_hadTruth_isHadFromC;   //!
   TBranch        *b_hadTruth_isHadFromB;   //!
   TBranch        *b_nOtherPV;   //!
   TBranch        *b_OtherPV_z;   //!
   TBranch        *b_PV_ndof;   //!
   TBranch        *b_PV_x;   //!
   TBranch        *b_PV_y;   //!
   TBranch        *b_PV_z;   //!
   TBranch        *b_PV_chi2;   //!
   TBranch        *b_PV_score;   //!
   TBranch        *b_PV_npvs;   //!
   TBranch        *b_PV_npvsGood;   //!
   TBranch        *b_nSV;   //!
   TBranch        *b_SV_dlen;   //!
   TBranch        *b_SV_dlenSig;   //!
   TBranch        *b_SV_pAngle;   //!
   TBranch        *b_had_chi2;   //!
   TBranch        *b_had_eta;   //!
   TBranch        *b_had_mass;   //!
   TBranch        *b_had_phi;   //!
   TBranch        *b_had_pt;   //!
   TBranch        *b_had_x;   //!
   TBranch        *b_had_y;   //!
   TBranch        *b_had_z;   //!
   TBranch        *b_had_ndof;   //!
   TBranch        *b_had_pdgId;   //!
   TBranch        *b_genHadron_eta;   //!
   TBranch        *b_genHadron_mass;   //!
   TBranch        *b_genHadron_phi;   //!
   TBranch        *b_genHadron_pt;   //!
   TBranch        *b_genHadron_x;   //!
   TBranch        *b_genHadron_y;   //!
   TBranch        *b_genHadron_z;   //!
   TBranch        *b_genHadron_pdgId;   //!
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
   TBranch        *b_HLT_Photon250_PFMET100;   //!
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
   TBranch        *b_HLT_PixelTracks_Multiplicity250ForEndOfFill;   //!
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
   TBranch        *b_HLT_L1FatHadTruthEvents;   //!
   TBranch        *b_HLT_Physics;   //!
   TBranch        *b_HLT_L1FatHadTruthEvents_part0;   //!
   TBranch        *b_HLT_L1FatHadTruthEvents_part1;   //!
   TBranch        *b_HLT_L1FatHadTruthEvents_part2;   //!
   TBranch        *b_HLT_L1FatHadTruthEvents_part3;   //!
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
   TBranch        *b_L1Reco_step;   //!
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
   TBranch        *b_Flag_BadChargedCandidateFilter;   //!
   TBranch        *b_Flag_BadPFMuonFilter;   //!
   TBranch        *b_Flag_BadChargedCandidateSummer16Filter;   //!
   TBranch        *b_Flag_BadPFMuonSummer16Filter;   //!
   TBranch        *b_Flag_trkPOG_manystripclus53X;   //!
   TBranch        *b_Flag_trkPOG_toomanystripclus53X;   //!
   TBranch        *b_Flag_trkPOG_logErrorTooManyClusters;   //!
   TBranch        *b_Flag_METFilters;   //!

   HadTruthEvents(TTree *tree=0);
   virtual ~HadTruthEvents();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef HadTruthEvents_cxx
HadTruthEvents::HadTruthEvents(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/home/iwatson/tsW/mcATnlo/data/tt01j_bbar_1lp_FxFx_AOD/HADTRUTHAOD/1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/home/iwatson/tsW/mcATnlo/data/tt01j_bbar_1lp_FxFx_AOD/HADTRUTHAOD/1.root");
      }
      f->GetObject("HadTruthEvents",tree);

   }
   Init(tree);
}

HadTruthEvents::~HadTruthEvents()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t HadTruthEvents::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t HadTruthEvents::LoadTree(Long64_t entry)
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

void HadTruthEvents::Init(TTree *tree)
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
   fChain->SetBranchAddress("nGenPart", &nGenPart, &b_nGenPart);
   fChain->SetBranchAddress("GenPart_eta", GenPart_eta, &b_GenPart_eta);
   fChain->SetBranchAddress("GenPart_mass", GenPart_mass, &b_GenPart_mass);
   fChain->SetBranchAddress("GenPart_phi", GenPart_phi, &b_GenPart_phi);
   fChain->SetBranchAddress("GenPart_pt", GenPart_pt, &b_GenPart_pt);
   fChain->SetBranchAddress("GenPart_genPartIdxMother", GenPart_genPartIdxMother, &b_GenPart_genPartIdxMother);
   fChain->SetBranchAddress("GenPart_pdgId", GenPart_pdgId, &b_GenPart_pdgId);
   fChain->SetBranchAddress("GenPart_status", GenPart_status, &b_GenPart_status);
   fChain->SetBranchAddress("GenPart_statusFlags", GenPart_statusFlags, &b_GenPart_statusFlags);
   fChain->SetBranchAddress("nhad_jet", &nhad_jet, &b_nhad_jet);
   fChain->SetBranchAddress("had_jet_btagCMVA", had_jet_btagCMVA, &b_had_jet_btagCMVA);
   fChain->SetBranchAddress("had_jet_btagCSVV2", had_jet_btagCSVV2, &b_had_jet_btagCSVV2);
   fChain->SetBranchAddress("had_jet_btagDeepB", had_jet_btagDeepB, &b_had_jet_btagDeepB);
   fChain->SetBranchAddress("had_jet_btagDeepC", had_jet_btagDeepC, &b_had_jet_btagDeepC);
   fChain->SetBranchAddress("had_jet_eta", had_jet_eta, &b_had_jet_eta);
   fChain->SetBranchAddress("had_jet_mass", had_jet_mass, &b_had_jet_mass);
   fChain->SetBranchAddress("had_jet_phi", had_jet_phi, &b_had_jet_phi);
   fChain->SetBranchAddress("had_jet_pt", had_jet_pt, &b_had_jet_pt);
   fChain->SetBranchAddress("nhad", &nhad, &b_nhad);
   fChain->SetBranchAddress("had_jetDR", had_jetDR, &b_had_jetDR);
   fChain->SetBranchAddress("had_legDR", had_legDR, &b_had_legDR);
   fChain->SetBranchAddress("had_diffMass", had_diffMass, &b_had_diffMass);
   fChain->SetBranchAddress("had_lxy", had_lxy, &b_had_lxy);
   fChain->SetBranchAddress("had_lxyErr", had_lxyErr, &b_had_lxyErr);
   fChain->SetBranchAddress("had_l3D", had_l3D, &b_had_l3D);
   fChain->SetBranchAddress("had_l3DErr", had_l3DErr, &b_had_l3DErr);
   fChain->SetBranchAddress("had_dca", had_dca, &b_had_dca);
   fChain->SetBranchAddress("had_angleXY", had_angleXY, &b_had_angleXY);
   fChain->SetBranchAddress("had_angleXYZ", had_angleXYZ, &b_had_angleXYZ);
   fChain->SetBranchAddress("had_dau1_chi2", had_dau1_chi2, &b_had_dau1_chi2);
   fChain->SetBranchAddress("had_dau1_nHits", had_dau1_nHits, &b_had_dau1_nHits);
   fChain->SetBranchAddress("had_dau1_pt", had_dau1_pt, &b_had_dau1_pt);
   fChain->SetBranchAddress("had_dau1_ipsigXY", had_dau1_ipsigXY, &b_had_dau1_ipsigXY);
   fChain->SetBranchAddress("had_dau1_ipsigZ", had_dau1_ipsigZ, &b_had_dau1_ipsigZ);
   fChain->SetBranchAddress("had_dau2_chi2", had_dau2_chi2, &b_had_dau2_chi2);
   fChain->SetBranchAddress("had_dau2_nHits", had_dau2_nHits, &b_had_dau2_nHits);
   fChain->SetBranchAddress("had_dau2_pt", had_dau2_pt, &b_had_dau2_pt);
   fChain->SetBranchAddress("had_dau2_ipsigXY", had_dau2_ipsigXY, &b_had_dau2_ipsigXY);
   fChain->SetBranchAddress("had_dau2_ipsigZ", had_dau2_ipsigZ, &b_had_dau2_ipsigZ);
   fChain->SetBranchAddress("had_nJet", had_nJet, &b_had_nJet);
   fChain->SetBranchAddress("had_nDau", had_nDau, &b_had_nDau);
   fChain->SetBranchAddress("had_dau1_charge", had_dau1_charge, &b_had_dau1_charge);
   fChain->SetBranchAddress("had_dau2_charge", had_dau2_charge, &b_had_dau2_charge);
   fChain->SetBranchAddress("had_idx", had_idx, &b_had_idx);
   fChain->SetBranchAddress("had_dau1_idx", had_dau1_idx, &b_had_dau1_idx);
   fChain->SetBranchAddress("had_dau2_idx", had_dau2_idx, &b_had_dau2_idx);
   fChain->SetBranchAddress("ngenHadron", &ngenHadron, &b_ngenHadron);
   fChain->SetBranchAddress("genHadron_dau1_pt", genHadron_dau1_pt, &b_genHadron_dau1_pt);
   fChain->SetBranchAddress("genHadron_dau2_pt", genHadron_dau2_pt, &b_genHadron_dau2_pt);
   fChain->SetBranchAddress("genHadron_dau1_eta", genHadron_dau1_eta, &b_genHadron_dau1_eta);
   fChain->SetBranchAddress("genHadron_dau2_eta", genHadron_dau2_eta, &b_genHadron_dau2_eta);
   fChain->SetBranchAddress("genHadron_dau1_phi", genHadron_dau1_phi, &b_genHadron_dau1_phi);
   fChain->SetBranchAddress("genHadron_dau2_phi", genHadron_dau2_phi, &b_genHadron_dau2_phi);
   fChain->SetBranchAddress("genHadron_vx", genHadron_vx, &b_genHadron_vx);
   fChain->SetBranchAddress("genHadron_vy", genHadron_vy, &b_genHadron_vy);
   fChain->SetBranchAddress("genHadron_vz", genHadron_vz, &b_genHadron_vz);
   fChain->SetBranchAddress("genHadron_isGenHadFromTsb", genHadron_isGenHadFromTsb, &b_genHadron_isGenHadFromTsb);
   fChain->SetBranchAddress("genHadron_dau1_pdgId", genHadron_dau1_pdgId, &b_genHadron_dau1_pdgId);
   fChain->SetBranchAddress("genHadron_dau2_pdgId", genHadron_dau2_pdgId, &b_genHadron_dau2_pdgId);
   fChain->SetBranchAddress("genHadron_isGenParticle", genHadron_isGenParticle, &b_genHadron_isGenParticle);
   fChain->SetBranchAddress("genHadron_isGenHadFromTop", genHadron_isGenHadFromTop, &b_genHadron_isGenHadFromTop);
   fChain->SetBranchAddress("genHadron_inVol", genHadron_inVol, &b_genHadron_inVol);
   fChain->SetBranchAddress("genHadron_isMatching", genHadron_isMatching, &b_genHadron_isMatching);
   fChain->SetBranchAddress("genHadron_isMatched", genHadron_isMatched, &b_genHadron_isMatched);
   fChain->SetBranchAddress("nhadTruth", &nhadTruth, &b_nhadTruth);
   fChain->SetBranchAddress("hadTruth_nMatched", hadTruth_nMatched, &b_hadTruth_nMatched);
   fChain->SetBranchAddress("hadTruth_nTrueDau", hadTruth_nTrueDau, &b_hadTruth_nTrueDau);
   fChain->SetBranchAddress("hadTruth_isHadFromTsb", hadTruth_isHadFromTsb, &b_hadTruth_isHadFromTsb);
   fChain->SetBranchAddress("hadTruth_isHadFromTop", hadTruth_isHadFromTop, &b_hadTruth_isHadFromTop);
   fChain->SetBranchAddress("hadTruth_isHadFromW", hadTruth_isHadFromW, &b_hadTruth_isHadFromW);
   fChain->SetBranchAddress("hadTruth_isHadFromS", hadTruth_isHadFromS, &b_hadTruth_isHadFromS);
   fChain->SetBranchAddress("hadTruth_isHadFromC", hadTruth_isHadFromC, &b_hadTruth_isHadFromC);
   fChain->SetBranchAddress("hadTruth_isHadFromB", hadTruth_isHadFromB, &b_hadTruth_isHadFromB);
   fChain->SetBranchAddress("nOtherPV", &nOtherPV, &b_nOtherPV);
   fChain->SetBranchAddress("OtherPV_z", OtherPV_z, &b_OtherPV_z);
   fChain->SetBranchAddress("PV_ndof", &PV_ndof, &b_PV_ndof);
   fChain->SetBranchAddress("PV_x", &PV_x, &b_PV_x);
   fChain->SetBranchAddress("PV_y", &PV_y, &b_PV_y);
   fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
   fChain->SetBranchAddress("PV_chi2", &PV_chi2, &b_PV_chi2);
   fChain->SetBranchAddress("PV_score", &PV_score, &b_PV_score);
   fChain->SetBranchAddress("PV_npvs", &PV_npvs, &b_PV_npvs);
   fChain->SetBranchAddress("PV_npvsGood", &PV_npvsGood, &b_PV_npvsGood);
   fChain->SetBranchAddress("nSV", &nSV, &b_nSV);
   fChain->SetBranchAddress("SV_dlen", SV_dlen, &b_SV_dlen);
   fChain->SetBranchAddress("SV_dlenSig", SV_dlenSig, &b_SV_dlenSig);
   fChain->SetBranchAddress("SV_pAngle", SV_pAngle, &b_SV_pAngle);
   fChain->SetBranchAddress("had_chi2", had_chi2, &b_had_chi2);
   fChain->SetBranchAddress("had_eta", had_eta, &b_had_eta);
   fChain->SetBranchAddress("had_mass", had_mass, &b_had_mass);
   fChain->SetBranchAddress("had_phi", had_phi, &b_had_phi);
   fChain->SetBranchAddress("had_pt", had_pt, &b_had_pt);
   fChain->SetBranchAddress("had_x", had_x, &b_had_x);
   fChain->SetBranchAddress("had_y", had_y, &b_had_y);
   fChain->SetBranchAddress("had_z", had_z, &b_had_z);
   fChain->SetBranchAddress("had_ndof", had_ndof, &b_had_ndof);
   fChain->SetBranchAddress("had_pdgId", had_pdgId, &b_had_pdgId);
   fChain->SetBranchAddress("genHadron_eta", genHadron_eta, &b_genHadron_eta);
   fChain->SetBranchAddress("genHadron_mass", genHadron_mass, &b_genHadron_mass);
   fChain->SetBranchAddress("genHadron_phi", genHadron_phi, &b_genHadron_phi);
   fChain->SetBranchAddress("genHadron_pt", genHadron_pt, &b_genHadron_pt);
   fChain->SetBranchAddress("genHadron_x", genHadron_x, &b_genHadron_x);
   fChain->SetBranchAddress("genHadron_y", genHadron_y, &b_genHadron_y);
   fChain->SetBranchAddress("genHadron_z", genHadron_z, &b_genHadron_z);
   fChain->SetBranchAddress("genHadron_pdgId", genHadron_pdgId, &b_genHadron_pdgId);
   fChain->SetBranchAddress("L1simulation_step", &L1simulation_step, &b_L1simulation_step);
   Notify();
}

Bool_t HadTruthEvents::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void HadTruthEvents::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t HadTruthEvents::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef HadTruthEvents_cxx
