#include "nano/analysis/interface/h2muAnalyser.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <cstdlib>

using namespace std;
/*
To compile:
g++ `root-config --cflags --glibs` -lEG h2muAnalyser.cc -o h2muAnalyser
*/

void h2muAnalyser::SetOutput(string outputName)
{
  m_output = TFile::Open(outputName.c_str(), "recreate");
  // TMVA Booking //
/*  weightXL = std::string(std::getenv("CMSSW_BASE"))+"/src/nano/analysis/test/h2mu/TMVA/newEvent/split/XL_50/weights/TMVAClassification_BDT.weights.xml";
  weightFH = std::string(std::getenv("CMSSW_BASE"))+"/src/nano/analysis/test/h2mu/TMVA/newEvent/split/nFH4_50/weights/TMVAClassification_BDT.weights.xml";
  weightOut = std::string(std::getenv("CMSSW_BASE"))+"/src/nano/analysis/test/h2mu/TMVA/newEvent/split/Out_100/weights/TMVAClassification_BDT.weights.xml";
  weightnoB = std::string(std::getenv("CMSSW_BASE"))+"/src/nano/analysis/test/h2mu/TMVA/newEvent/nonB_120/weights/TMVAClassification_BDT.weights.xml";

  bdt_XL = new TMVA::Reader();
  bdt_XL->AddVariable( "Met", &b_Met );
  bdt_XL->AddVariable( "all_muEtaDiff", &b_all_muEtaDiff );
  bdt_XL->AddVariable( "all_muPtDiff", &b_all_muPhiDiff );
  bdt_XL->AddVariable( "all_muPhiDiff", &b_all_muPhiDiff );
  bdt_XL->AddVariable( "all_muDR", &b_all_muDR );
  bdt_XL->AddVariable( "all_Dilep_Pt", &b_all_Dilep_Pt );
  bdt_XL->AddVariable( "all_Dilep_Eta", &b_all_Dilep_Eta );
  bdt_XL->AddVariable( "nelec", &b_nelec );
  bdt_XL->AddVariable( "nmuon", &b_nmuon );
  bdt_XL->AddVariable( "nnonbjet", &b_nnonbjet );
  bdt_XL->AddVariable( "nbjet", &b_nbjet );
  bdt_XL->AddVariable( "Central_Jets", &b_Central_Jets );
  bdt_XL->AddVariable( "Forward_Jets", &b_Forward_Jets );
  bdt_XL->AddVariable( "minDR", &b_minDR );
  bdt_XL->AddVariable( "XlepPt", &b_XlepPt );
  bdt_XL->AddVariable( "mT2", &b_mT2 );
  bdt_XL->AddVariable( "mT", &b_mT );
  bdt_XL->AddVariable( "DiJetM12", &b_DiJetM12 );
  bdt_XL->BookMVA("BDT", weightXL); 
  
  bdt_FH = new TMVA::Reader();
  bdt_FH->AddVariable( "Met", &b_Met );
  bdt_FH->AddVariable( "all_muEtaDiff", &b_all_muEtaDiff );
  bdt_FH->AddVariable( "all_muPtDiff", &b_all_muPhiDiff );
  bdt_FH->AddVariable( "all_muPhiDiff", &b_all_muPhiDiff );
  bdt_FH->AddVariable( "all_muDR", &b_all_muDR );
  bdt_FH->AddVariable( "all_Dilep_Pt", &b_all_Dilep_Pt );
  bdt_FH->AddVariable( "all_Dilep_Eta", &b_all_Dilep_Eta );
  bdt_FH->AddVariable( "nnonbjet", &b_nnonbjet );
  bdt_FH->AddVariable( "nbjet", &b_nbjet );
  bdt_FH->AddVariable( "Central_Jets", &b_Central_Jets );
  bdt_FH->AddVariable( "Forward_Jets", &b_Forward_Jets );
  bdt_FH->AddVariable( "minDR1", &b_minDR1 );
  bdt_FH->AddVariable( "minDR2", &b_minDR2 );
  bdt_FH->AddVariable( "mT2", &b_mT2 );
  bdt_FH->AddVariable( "mT", &b_mT );
  bdt_FH->AddVariable( "DiJetM12", &b_DiJetM12 );
  bdt_FH->AddVariable( "DiJetM13", &b_DiJetM13 );
  bdt_FH->AddVariable( "DiJetM14", &b_DiJetM14 );
  bdt_FH->AddVariable( "DiJetM23", &b_DiJetM23 );
  bdt_FH->AddVariable( "DiJetM24", &b_DiJetM24 );
  bdt_FH->AddVariable( "DiJetM34", &b_DiJetM34 );
  bdt_FH->BookMVA("BDT", weightFH); 
  
  bdt_Out = new TMVA::Reader();
  bdt_Out->AddVariable( "Met", &b_Met );
  bdt_Out->AddVariable( "all_muEtaDiff", &b_all_muEtaDiff );
  bdt_Out->AddVariable( "all_muPtDiff", &b_all_muPhiDiff );
  bdt_Out->AddVariable( "all_muPhiDiff", &b_all_muPhiDiff );
  bdt_Out->AddVariable( "all_muDR", &b_all_muDR );
  bdt_Out->AddVariable( "all_Dilep_Pt", &b_all_Dilep_Pt );
  bdt_Out->AddVariable( "all_Dilep_Eta", &b_all_Dilep_Eta );
  bdt_Out->AddVariable( "nnonbjet", &b_nnonbjet );
  bdt_Out->AddVariable( "nbjet", &b_nbjet );
  bdt_Out->AddVariable( "Central_Jets", &b_Central_Jets );
  bdt_Out->AddVariable( "Forward_Jets", &b_Forward_Jets );
  bdt_Out->AddVariable( "mT2", &b_mT2 );
  bdt_Out->AddVariable( "mT", &b_mT );
  bdt_Out->AddVariable( "DiJetM12", &b_DiJetM12 );
  bdt_Out->BookMVA("BDT", weightOut);

  bdt_noB = new TMVA::Reader();
  bdt_noB->AddVariable( "Met", &b_Met );
//  bdt_noB->AddVariable( "all_muEtaDiff", &b_all_muEtaDiff );
//  bdt_noB->AddVariable( "all_muPtDiff", &b_all_muPhiDiff );
//  bdt_noB->AddVariable( "all_muPhiDiff", &b_all_muPhiDiff );
//  bdt_noB->AddVariable( "all_muDR", &b_all_muDR );
  bdt_noB->AddVariable( "all_Dilep_Pt", &b_all_Dilep_Pt );
  bdt_noB->AddVariable( "all_Dilep_Eta", &b_all_Dilep_Eta );
 // bdt_noB->AddVariable( "nelec", &b_nelec );
 // bdt_noB->AddVariable( "nexLep", &b_nexLep );
 // bdt_noB->AddVariable( "nmuon", &b_nmuon );
  bdt_noB->AddVariable( "Central_Jets", &b_Central_Jets );
  bdt_noB->AddVariable( "Forward_Jets", &b_Forward_Jets );
  bdt_noB->AddVariable( "etaJ1", &b_etaJ1 );
  bdt_noB->AddVariable( "etaJ2", &b_etaJ2 );
  bdt_noB->AddVariable( "DijetM1", &b_DijetM1 );
  bdt_noB->AddVariable( "DijetM2", &b_DijetM2 );
  bdt_noB->AddVariable( "DijetEta1", &b_DijetEta1 );
  bdt_noB->AddVariable( "DijetEta2", &b_DijetEta2 );
  bdt_noB->BookMVA("BDT", weightnoB);
*/
  m_tree = new TTree("events", "events");
  MakeBranch(m_tree);
  
  h_Event_Tot = new TH1D("Event_total", "Event_total" ,1,0,1);
  h_genweights = new TH1D("genweight", "genweight" , 1,0,1);
  h_weight = new TH1D("weight", "weight", 1,0,1);
  h_cutFlow = new TH1D("cutflow", "cutflow", 11, -0.5, 10.5);
  h_XL = new TH1D("ExtraLep", "ExtraLep", 1, 0, 1);
  h_nFH4 = new TH1D("FullyHad4", "FullyHad4", 1, 0, 1);
  h_Out = new TH1D("Out", "Out", 1, 0, 1);
  h_Non = new TH1D("Non", "Non", 1, 0, 1);

}
void h2muAnalyser::MakeBranch(TTree* t)
{
  t->Branch("Event_No", &b_Event_No, "Event_No/I");
  t->Branch("Step", &b_Step, "Step/I");
  t->Branch("Dilep", "TLorentzVector", &b_Dilep);
  t->Branch("Mu1", "TLorentzVector", &b_Mu1);
  
  t->Branch("Mu2", "TLorentzVector", &b_Mu2);
  t->Branch("lep", "TLorentzVector", &b_lep);
  
  t->Branch("Jet_pu", "std::vector<float>", &b_Jet_pu);
  t->Branch("Jet_pt", "std::vector<float>", &b_Jet_pT);
  t->Branch("Jet_eta", "std::vector<float>", &b_Jet_Eta);
  t->Branch("Jet_phi", "std::vector<float>", &b_Jet_Phi);
  t->Branch("nonbJet_pt", "std::vector<float>", &b_nonbJet_pT);
 
  t->Branch("nonbJet_eta", "std::vector<float>", &b_nonbJet_Eta);
  t->Branch("nonbJet_phi", "std::vector<float>", &b_nonbJet_Phi);
  t->Branch("genweight", &b_genweight, "genweight/F");
  t->Branch("puweight", &b_puweight, "puweight/F");
  
  t->Branch("puweight_up", &b_puweight_up, "puweight_up/F");
  t->Branch("puweight_dn", &b_puweight_dn, "puweight_dn/F");
  t->Branch("weight", &b_weight, "weight/F");
  t->Branch("npvs", &b_npvs, "npvs/F");
  
  t->Branch("channel", &b_channel, "channel/I");
  t->Branch("charge", &b_charge, "charge/I");
  t->Branch("nlep", &b_nlep, "nlep/F");
  t->Branch("nmuon", &b_nmuon, "nmuon/F");
  
  t->Branch("nelec", &b_nelec, "nelec/F");
  t->Branch("njet", &b_njet, "njet/F");
  t->Branch("nnonbjet", &b_nnonbjet, "nnonbjet/F");
  t->Branch("nbjet", &b_nbjet, "nbjet/F");
  
  t->Branch("trig_m", &b_trig_m, "trig_m/O");
  t->Branch("Met", &b_Met, "Met/F");
  t->Branch("Met_phi", &b_Met_phi, "Met_phi/F");
  t->Branch("CSVv2", &b_CSVv2);
  
  t->Branch("nFH4", &b_nFH4, "nFH4/I");
  t->Branch("XL", &b_XL, "XL/I");
  t->Branch("Out", &b_Out, "Out/I");
  t->Branch("nonB", &b_nonB, "nonB/I");
  
  t->Branch("csvweight", "std::vector<float>", &b_csvweights);
//  t->Branch("btagweight", &b_btagweight, "btagweight/F");
//  t->Branch("btagweight_up", &b_btagweight_up, "btagweight_up/F");
//  t->Branch("btagweight_dn", &b_btagweight_dn, "btagweight_dn/F");
  
  /*t->Branch("jes_up", &b_jes_up, "jes_up/F");
  t->Branch("jes_dn", &b_jes_dn, "jes_dn/F");
  t->Branch("cferr1_up", &b_cferr1_up, "cferr1_up/F");
  t->Branch("cferr1_dn", &b_cferr1_dn, "cferr1_dn/F");
  
  t->Branch("cferr2_up", &b_cferr2_up, "cferr2_up/F");
  t->Branch("cferr2_dn", &b_cferr2_dn, "cferr2_dn/F");
  t->Branch("hf_up", &b_hf_up, "hf_up/F");
  t->Branch("hf_dn", &b_hf_dn, "hf_dn/F");
  
  t->Branch("hfstats1_up", &b_hfstats1_up, "hfstats1_up/F");
  t->Branch("hfstats1_dn", &b_hfstats1_dn, "hfstats1_dn/F");
  t->Branch("hfstats2_up", &b_hfstats2_up, "hfstats2_up/F");
  t->Branch("hfstats2_dn", &b_hfstats2_dn, "hfstats2_dn/F");
  
  t->Branch("lf_up", &b_lf_up, "lf_up/F");
  t->Branch("lf_dn", &b_lf_dn, "lf_dn/F");
  t->Branch("lfstats1_up", &b_lfstats1_up, "lfstats1_up/F");
  t->Branch("lfstats1_dn", &b_lfstats1_dn, "lfstats1_dn/F");
  
  t->Branch("lfstats2_up", &b_lfstats2_up, "lfstats2_up/F");
  t->Branch("lfstats2_dn", &b_lfstats2_dn, "lfstats2_dn/F");*/
  t->Branch("mueffweight", &b_mueffweight, "mueffweight/F");
  t->Branch("mueffweight_up", &b_mueffweight_up, "mueffweight_up/F");
  t->Branch("mueffweight_dn", &b_mueffweight_dn, "mueffweight_dn/F");
 
  ////////////////////BDT////////////////////////////
  t->Branch("all_muEtaDiff", &b_all_muEtaDiff, "all_muEtaDiff/F");
  t->Branch("all_muPtDiff", &b_all_muPtDiff, "all_muPtDiff/F");
  t->Branch("all_muPhiDiff", &b_all_muPhiDiff, "all_muPhiDiff/F");
  t->Branch("all_muDR", &b_all_muDR, "all_muDR/F");
  
  t->Branch("all_Dilep_Pt", &b_all_Dilep_Pt, "all_Dilep_Pt/F");
  t->Branch("all_Dilep_Eta", &b_all_Dilep_Eta, "all_Dilep_Eta/F");
  t->Branch("all_Dilep_Phi", &b_all_Dilep_Phi, "all_Dilep_Phi/F");
  
  t->Branch("DijetM1", &b_DijetM1, "DijetM1/F");
  t->Branch("DijetM2", &b_DijetM2, "DijetM2/F");
  t->Branch("DijetEta1",&b_DijetEta1, "DijetEta1/F");
  t->Branch("DijetEta2",&b_DijetEta2, "DijetEta2/F");
  
  t->Branch("DiJetM12", &b_DiJetM12, "DiJetM12/F");
  t->Branch("DiJetM13", &b_DiJetM13, "DiJetM13/F");
  t->Branch("DiJetM14", &b_DiJetM14, "DiJetM14/F");
  t->Branch("DiJetM23", &b_DiJetM23, "DiJetM23/F");
  
  t->Branch("DiJetM24", &b_DiJetM24, "DiJetM24/F");
  t->Branch("DiJetM34", &b_DiJetM34, "DiJetM34/F");
  t->Branch("minDR1", &b_minDR1, "minDR1/F");
  t->Branch("minDR2", &b_minDR2, "minDR2/F");
  
  t->Branch("minDR", &b_minDR, "minDR/F");
  t->Branch("XlepPt", &b_XlepPt, "XlepPt/F");
  t->Branch("mT2", &b_mT2, "mT2/F");
  t->Branch("mT", &b_mT, "mT/F");
  
  t->Branch("nexLep", &b_nexLep, "nexLep/F");
  t->Branch("etaJ1", &b_etaJ1, "etaJ1/F");
  t->Branch("etaJ2", &b_etaJ2, "etaJ2/F");
  t->Branch("Central_Jets", &b_Central_Jets, "Central_Jets/F");
  
  t->Branch("Forward_Jets", &b_Forward_Jets, "Forward_Jets/F");
  t->Branch("MVA_BDTXL", &b_MVA_BDTXL, "MVA_BDTXL/F");
  t->Branch("MVA_BDTFH", &b_MVA_BDTFH, "MVA_BDTFH/F");
  t->Branch("MVA_BDTnoB", &b_MVA_BDTnoB, "MVA_BDTnoB/F");
  
  t->Branch("MVA_BDTOut", &b_MVA_BDTOut, "MVA_BDTOut/F");
  t->Branch("CSV", &b_CSV, "CSV/F");
}

void h2muAnalyser::ResetBranch()
{  
  b_Event_No = 0;
  b_Step = 0;
  b_Dilep.SetPtEtaPhiM(0,0,0,0);
  b_Mu1.SetPtEtaPhiM(0,0,0,0);
  b_Mu2.SetPtEtaPhiM(0,0,0,0);
  b_lep.SetPtEtaPhiM(0,0,0,0);
  b_nonbJet1.SetPtEtaPhiM(0,0,0,0);
  b_nonbJet2.SetPtEtaPhiM(0,0,0,0);
  b_nonbJet3.SetPtEtaPhiM(0,0,0,0);
  b_nonbJet4.SetPtEtaPhiM(0,0,0,0);

  b_Mu_tlv.clear();
  b_El_tlv.clear();
  b_Jet_tlv.clear();
  b_bJet_tlv.clear();
  b_nonbJet_tlv.clear(); b_Jet_pT.clear(); b_Jet_Eta.clear(); b_Jet_pu.clear();
  b_Jet_pT.clear(); b_Jet_Eta.clear(); b_Jet_Phi.clear();
  b_nonbJet_pT.clear(); b_nonbJet_Eta.clear(); b_nonbJet_Phi.clear();
  b_csvweights.clear(); b_CSVv2.clear();
  b_Event_Total = 1;
  b_channel = -1;
  b_nlep = 0; b_nmuon = 0; b_nelec = 0; b_nnonbjet = 0; b_njet = 0; b_nbjet = 0;
  b_XL = 0; b_nFH4 = 0; b_Out = 0; b_nonB = 0; 
  // BDT //
  b_all_muEtaDiff = 0; b_all_muPtDiff = 0; b_all_muPhiDiff = 0; b_all_muDR = 0;
  b_all_Dilep_Pt = 0; b_all_Dilep_Eta = 0; b_all_Dilep_Phi = 0;
  b_Central_Jets = 0; b_Forward_Jets = 0;
  b_DiJetM12 = 0; b_DiJetM13 = 0; b_DiJetM14 = 0; b_DiJetM23 = 0; b_DiJetM24 = 0; b_DiJetM34 = 0;
  b_minDR1 = -2; b_minDR2 = -2;
  b_minDR = -2; b_XlepPt = 0; b_mT2 = -2; b_mT = -2; b_nexLep = 0; b_etaJ1 = -6; b_etaJ2 = -6;
  b_MVA_BDTXL = -999; b_MVA_BDTFH = -999; b_MVA_BDTnoB = -999; b_MVA_BDTOut = -999;
  DR_Hold = 0, DR_Hold2 = 0 ,mT_Hold = 0, MuPT_Hold = 0, ElPT_Hold = 0;
  b_DijetEta1 = 0, b_DijetEta2 = 0, DijetEta_hold = 0, DijetM_hold = 0, b_DijetM1 = 0, b_DijetM2 = 0;
  b_CSV = 0;
  b_npvs = 0;
}

bool h2muAnalyser::Analysis()
{
  h_Event_Tot->Fill(0.5, b_Event_Total);
  h_cutFlow->Fill(0);
  b_Step = 0;
  
  //Run for MC
  if(m_isMC){
    Int_t nvtx = Pileup_nTrueInt;
    b_puweight = m_pileUp->getWeight(nvtx);
    b_puweight_up = m_pileUp->getWeight(nvtx, 1);
    b_puweight_dn = m_pileUp->getWeight(nvtx, -1);
    b_genweight = genWeight;
    h_genweights->Fill(0.5, b_genweight);
    b_weight = b_genweight * b_puweight;
    h_weight->Fill(0.5, b_weight);
  } else {
    b_puweight = 1;
    b_genweight = 0;
    if (!(m_lumi->LumiCheck(run, luminosityBlock))) return false;
  }
  b_Step = 1;
  h_cutFlow->Fill(1);

  if (fabs(PV_z) >= 24.) return false;
  if (PV_npvs == 0) return false;
  if (PV_ndof < 4) return false;
  b_Step = 2;
  h_cutFlow->Fill(2);
  
  b_trig_m = HLT_IsoTkMu24 || HLT_IsoMu24;
//  if (!(HLT_IsoTkMu24 || HLT_IsoMu24)) return false;
  Bool_t IsoMu24 = false;
  Bool_t IsoTkMu24 = false;

  for (UInt_t i = 0; i < nTrigObj; ++i) {
    if (TrigObj_id[i] != 13) continue; 
    if (TrigObj_pt[i] < 24) continue;
    Int_t bits = TrigObj_filterBits[i];
    if (bits & 0x2) IsoMu24 = true;
    if (bits & 0x8) IsoTkMu24 = true;  
  }
  if (!(IsoMu24 || IsoTkMu24)) return false; 
  
  b_Step = 3;
  h_cutFlow->Fill(3);

  auto Muons = MuonSelection();
  auto Elecs = ElectronSelection(Muons);
 
  if (Muons.size() < 2) return false;
//if (Muons.size() != 2) return false;
  
  TParticle mu1;
  TParticle mu2; 
  b_Met = PuppiMET_pt;
  b_Met_phi = PuppiMET_phi;
  b_Step = 4;
  h_cutFlow->Fill(4);

  for (UInt_t i = 0; i < Muons.size(); i++) { 
    if ((b_Mu_tlv[0].Pt() > 26) || (b_Mu_tlv[i].Pt() > 26)) { 
      if ((Muons[0].GetPdgCode() * Muons[i].GetPdgCode()) < 0 ) { 
        b_Mu1 = b_Mu_tlv[0];
        b_Mu2 = b_Mu_tlv[i];

        b_charge = 1;
        mu1 = Muons[0];
        mu2 = Muons[i];
        break;
      }
    }
  }

  if (b_charge == 0) return false; 

  b_nlep = Muons.size() + Elecs.size();
  b_nmuon = Muons.size();
  b_nelec = Elecs.size();
  b_Step = 5;
  h_cutFlow->Fill(5);

  auto Jets = JetSelection(Muons, Elecs);
  auto BJets = BJetSelection(Muons, Elecs);
  auto nonbJets = nonbJetSelection(Muons, Elecs);
  
  b_njet = Jets.size();
  b_nnonbjet = nonbJets.size();
  b_nbjet = BJets.size();
  b_npvs = PV_npvs;

  for (UInt_t i = 0; i < Jets.size(); i++) {
     b_Jet_pT.push_back(b_Jet_tlv[i].Pt());
     b_Jet_Eta.push_back(b_Jet_tlv[i].Eta());
     b_Jet_Phi.push_back(b_Jet_tlv[i].Phi());
  }
  for (UInt_t i = 0; i < nonbJets.size(); i++) {
     b_nonbJet_pT.push_back(b_nonbJet_tlv[i].Pt());
     b_nonbJet_Eta.push_back(b_nonbJet_tlv[i].Eta());
     b_nonbJet_Phi.push_back(b_nonbJet_tlv[i].Phi());
  } 

  b_Dilep = b_Mu1 + b_Mu2;
  // Extra lep //
  if (Muons.size() + Elecs.size()  >= 3 && BJets.size() >= 1) {
     b_XL = 1;
     h_XL->Fill(0.5, b_XL); 

     b_CSV = b_CSVv2[0];
     if (Muons.size() > 2) {
       for (UInt_t k=1; k < Muons.size(); k++) {
         if (b_Mu_tlv[k].Pt() != b_Mu2.Pt()) {
           if (b_Mu_tlv[k].Pt() > MuPT_Hold) {
             MuPT_Hold = b_Mu_tlv[k].Pt();
           }
           mT_Hold = sqrt(2.0*b_Mu_tlv[k].Pt()*b_Met*(1-cos(fabs(b_Mu_tlv[k].Eta()-b_Met_phi))));
           if (mT_Hold > b_mT2) {
              b_mT2 = mT_Hold; 
           }
         }
       } 
     }
     if (Elecs.size() > 0) {
       ElPT_Hold = b_El_tlv[0].Pt();
       for (UInt_t k=0; k < Elecs.size(); k++) {
         mT_Hold = sqrt(2.0*b_El_tlv[k].Pt()*b_Met*(1-cos(fabs(b_El_tlv[k].Eta()-b_Met_phi))));
         if (mT_Hold > b_mT2) {
           b_mT2 = mT_Hold; 
         }
       }
     }
     if (ElPT_Hold > MuPT_Hold) {
       b_XlepPt = ElPT_Hold;
     } else {
       b_XlepPt = MuPT_Hold; 
     }

     if (nonbJets.size() >= 2) {
       nonbJets[0].Momentum(b_nonbJet1); 
       nonbJets[1].Momentum(b_nonbJet2); 
       b_DiJetM12 = b_nonbJet1.M() + b_nonbJet2.M();
     }

     for (UInt_t i = 0; i < BJets.size(); i++) {
       if (Elecs.size() > 0) {
         for (UInt_t m = 0; m < Elecs.size(); m++) { 
           DR_Hold = b_bJet_tlv[i].DeltaR(b_El_tlv[m]);
           if (DR_Hold < b_minDR || b_minDR < 0) {
             b_minDR = DR_Hold;
           }
         }
       }
       if (Muons.size() > 2) {
         for (UInt_t n = 1; n < Muons.size(); n++) {
           if (Muons[n].Pt() != b_Mu2.Pt()) {
             DR_Hold = b_bJet_tlv[i].DeltaR(b_Mu_tlv[n]);
             if (DR_Hold < b_minDR || b_minDR < 0) {
               b_minDR = DR_Hold;
             } 
           }
         }
       }
     }
     
     b_mT = sqrt(2.0*b_Mu1.Pt()*b_Met*(1-cos(fabs(b_Mu1.Eta()-b_Met_phi))));
  }
  // Hadronic //
  if (Elecs.size() == 0 && Muons.size() == 2 && BJets.size() >= 1 && nonbJets.size() >= 4) {
     b_nFH4 = 1;
     b_CSV = b_CSVv2[0];
     nonbJets[0].Momentum(b_nonbJet1); 
     nonbJets[1].Momentum(b_nonbJet2); 
     nonbJets[2].Momentum(b_nonbJet3); 
     nonbJets[3].Momentum(b_nonbJet4); 
     
     b_DiJetM12 = b_nonbJet1.M() + b_nonbJet2.M();
     b_DiJetM13 = b_nonbJet1.M() + b_nonbJet3.M();
     b_DiJetM14 = b_nonbJet1.M() + b_nonbJet4.M();
     b_DiJetM23 = b_nonbJet2.M() + b_nonbJet3.M();
     b_DiJetM24 = b_nonbJet2.M() + b_nonbJet4.M();
     b_DiJetM34 = b_nonbJet3.M() + b_nonbJet4.M();
     
     for (UInt_t i = 0; i < BJets.size(); i++) {
       DR_Hold = b_bJet_tlv[i].DeltaR(b_Mu1);
       DR_Hold2 = b_bJet_tlv[i].DeltaR(b_Mu2);
       if (DR_Hold < b_minDR1 || b_minDR1 < 0) {
         b_minDR1 = DR_Hold;
       }
       if (DR_Hold2 > b_minDR2 || b_minDR2 < 0) {
         b_minDR2 = DR_Hold2;
       }
     }
     
     b_mT2 = sqrt(2.0*b_Mu1.Pt()*b_Met*(1-cos(fabs(b_Mu1.Eta()-b_Met_phi))));
     b_mT = sqrt(2.0*b_Mu2.Pt()*b_Met*(1-cos(fabs(b_Mu2.Eta()-b_Met_phi))));
     h_nFH4->Fill(0.5, b_nFH4); 
  }
  //Outsider //
  if (BJets.size() >= 1 && b_nFH4 == 0 && b_XL == 0) {
     b_Out = 1;
     b_CSV = b_CSVv2[0];
     if (nonbJets.size() >= 2) {
     nonbJets[0].Momentum(b_nonbJet1);
     nonbJets[1].Momentum(b_nonbJet2);
     
     b_DiJetM12 = b_nonbJet1.M() + b_nonbJet2.M();
     }     

     b_mT2 = sqrt(2.0*b_Mu1.Pt()*b_Met*(1-cos(fabs(b_Mu1.Eta()-b_Met_phi))));
     b_mT = sqrt(2.0*b_Mu2.Pt()*b_Met*(1-cos(fabs(b_Mu2.Eta()-b_Met_phi))));
     
     h_Out->Fill(0.5,b_Out);
  }
  // nob // 
  if (BJets.size() == 0) {
     b_nonB = 1;
     b_nexLep = Muons.size() + Elecs.size() - 2;
     if (Jets.size() >= 2) {
       TLorentzVector j1; 
       TLorentzVector j2; 
       b_etaJ1 = b_Jet_tlv[0].Eta();
       b_etaJ2 = b_Jet_tlv[1].Eta();
       for (auto& J1 : b_Jet_tlv) {
         for (auto& J2 : b_Jet_tlv) {
           if (J1.M() != J2.M()) { 
             TLorentzVector Dijet_hold = J1 + J2;
             if (Dijet_hold.M() > b_DijetM2) {
               b_DijetM2 = Dijet_hold.M();
               b_DijetEta2 = Dijet_hold.Eta();
             }    
             if (b_DijetM1 == b_DijetM2) {
               b_DijetM2 = 0; 
             }    
             if (b_DijetM2 > b_DijetM1) {
               DijetM_hold = b_DijetM1;
               DijetEta_hold = b_DijetEta1;

               b_DijetM1 = b_DijetM2;
               b_DijetEta1 = b_DijetEta2;

               b_DijetM2 = DijetM_hold;
               b_DijetEta2 = DijetEta_hold; 
            }    
          }    
        }    
      }    
    }
    h_Non->Fill(0.5, b_nonB);
  }

  ////////// BDT //////////////////////////////
  /// BDT Varialbe Set ALL ///
  b_all_muEtaDiff = fabs(b_Mu1.Eta() - b_Mu2.Eta());
  b_all_muPtDiff = fabs(b_Mu1.Pt() - b_Mu2.Pt());
  //b_all_muPhiDiff = fabs(b_Mu1.Phi() - b_Mu2.Phi());
  b_all_muPhiDiff = b_Mu1.DeltaPhi(b_Mu2);
  b_all_muDR = b_Mu1.DeltaR(b_Mu2);
  b_all_Dilep_Pt = b_Dilep.Pt();   
  b_all_Dilep_Eta = b_Dilep.Eta();   
  b_all_Dilep_Phi = b_Dilep.Phi();   

  for (auto& J : Jets) {
    if (fabs(J.Eta()) < 2.4) {
       b_Central_Jets += 1;   
    } else {
       b_Forward_Jets += 1; 
    }
  }
  b_mueffweight = m_muonSF.getScaleFactor(mu1, 13, 0)*m_muonSF.getScaleFactor(mu2, 13, 0);
  b_mueffweight_up = m_muonSF.getScaleFactor(mu1, 13, +1)*m_muonSF.getScaleFactor(mu2, 13, +1);
  b_mueffweight_dn = m_muonSF.getScaleFactor(mu1, 13, -1)*m_muonSF.getScaleFactor(mu2, 13, -1);

  b_Event_No = 1;
 /* b_MVA_BDTXL = bdt_XL->EvaluateMVA("BDT");
  b_MVA_BDTFH = bdt_FH->EvaluateMVA("BDT");
  b_MVA_BDTnoB = bdt_noB->EvaluateMVA("BDT");
  b_MVA_BDTOut = bdt_Out->EvaluateMVA("BDT");
*/
  return true;
}


void h2muAnalyser::Loop()
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();
  
  // Events loop
  for (Long64_t iev=0; iev<nentries; iev++) {
    ResetBranch();
    fChain->GetEntry(iev);
    bool keep = Analysis();
    if (keep) { m_tree->Fill(); }
  }

}

int main(Int_t argc, Char_t** argv)
{
  string env = getenv("CMSSW_BASE");
  string username = getenv("USER");

  if (argc != 1) {
    string dirName = "root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/"+username+"/nanoAOD/"+std::string(argv[1])+"/"+std::string(argv[2]);
   // string dirName = env+("/src/nano/analysis/test/h2mu/Results/")+argv[1]+"/"+argv[2];
    string temp = argv[2];
    Bool_t isMC = false;
    Size_t found = temp.find("Run");
    if(found == string::npos) isMC = true;
    for (Int_t i = 3; i < argc; i++) {
      TFile *f = TFile::Open(argv[i], "read");

      TTree *tree;
      f->GetObject("Events", tree);

      temp = argv[i];
      found = temp.find_last_of('/');
      string outPutName = dirName+temp.substr(found);
      h2muAnalyser t(tree, isMC);
      t.SetOutput(outPutName);
      t.Loop();
    }
  } else {
  //  TFile *f = TFile::Open("/xrootd/store/group/nanoAOD/run2_2016v4/ttHToMuMu_M125_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6/180509_131219/0000/nanoAOD_112.root", "read");
  //  TFile *f = TFile::Open("/xrootd/store/group/nanoAOD/run2_2016v3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/180125_131129/0000/nanoAOD_100.root", "read");
    TFile *f = TFile::Open("/xrootd/store/group/nanoAOD/run2_2016v4/SingleMuon/Run2016C-07Aug17-v1/180504_150105/0000/nanoAOD_140.root", "read");
    TTree *tree;
    Bool_t isMC = false;
    string temp = "Run";
    Size_t found = temp.find("Run");
    if(found == string::npos) isMC = true;
    f->GetObject("Events", tree);

      h2muAnalyser t(tree, isMC);    
    t.SetOutput("test.root");
    t.Loop();
  }

  return 0;
}

//Object Selections
vector<TParticle> h2muAnalyser::MuonSelection()
{
  vector<TParticle> muons;
  for (UInt_t i = 0; i < nMuon; i++) {
    if (!Muon_trackerMu[i]) continue;
    if (!Muon_globalMu[i]) continue;
    if (!Muon_tightId[i]) continue;
   // if (!Muon_mediumId[i]) continue;
    if (Muon_pfRelIso04_all[i] > 0.25) continue;
    
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]);
    mom = mom * roccoR(mom, Muon_charge[i], Muon_genPartIdx[i], Muon_nTrackerLayers[i]);
    
    if (mom.Pt() < 10) continue;
    if (fabs(mom.Eta() > 2.4)) continue;
   
    auto muon = TParticle();
    muon.SetPdgCode(13*Muon_charge[i]*-1);
    muon.SetMomentum(mom);

    b_Mu_tlv.push_back(mom);
    muons.push_back(muon);
  }
  return muons;
}

vector<TParticle> h2muAnalyser::ElectronSelection(vector<TParticle> leptons)
{
  vector<TParticle> electrons;
  for(UInt_t i = 0; i < nElectron; i++) {
    if ( Electron_pt[i] < 10) continue;
    if (fabs(Electron_eta[i]) > 2.5 ) continue; //<~~~~~~~~~~~~~~ Higgs Electron pt == 10; Higgs Electron eta > 2.5  
    if( Electron_pfRelIso03_all[i] > 0.15) continue;
    if (Electron_cutBased[i] < 3) continue;
    float el_scEta = Electron_deltaEtaSC[i] + Electron_eta[i];
    if ( std::fabs(el_scEta) > 1.4442 &&  std::fabs(el_scEta) < 1.566 ) continue;

    TLorentzVector mom;
    mom.SetPtEtaPhiM(Electron_pt[i], Electron_eta[i], Electron_phi[i], Electron_mass[i]);

    if (hasOverLap(mom, leptons, 0.4)) continue;
    
    auto elec = TParticle();
    elec.SetPdgCode(11*Electron_charge[i]*-1);
    elec.SetMomentum(mom);
    
    b_El_tlv.push_back(mom);
    electrons.push_back(elec);
  }
  return electrons;
}
vector<TParticle> h2muAnalyser::JetSelection(vector<TParticle> Muons, vector<TParticle> Elecs)
{
  vector<TParticle> jets;
  for (UInt_t i = 0; i < nJet; i++) {
//    if (fabs(Jet_eta[i]) >= 2.4 &&  Jet_pt[i] < 30) continue;
//    if (fabs(Jet_eta[i]) < 2.4 &&  Jet_pt[i] < 20) continue;
    if (Jet_pt[i] < 30) continue;
//    if (fabs(Jet_eta[i]) > 2.4) continue; 
    if (fabs(Jet_eta[i]) > 4.7) continue; 
    if (Jet_jetId[i] < 1) continue;

    TLorentzVector mom;
    mom.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
    
    if (hasOverLap(mom, Muons, 0.4)) continue;
    if (hasOverLap(mom, Elecs, 0.4)) continue;
    
    auto jet = TParticle();
    jet.SetMomentum(mom);
    b_Jet_tlv.push_back(mom);
    jets.push_back(jet);

    // Jet Pu_ID 
    b_Jet_pu.push_back(Jet_puId[i]);
  }
  return jets;
}

vector<TParticle> h2muAnalyser::BJetSelection(vector<TParticle> Muons, vector<TParticle> Elecs)
{
  vector<TParticle> bJets;
//  b_btagweight = 1.0;
//  b_btagweight_up = 1.0;
//  b_btagweight_dn = 1.0;
  for (UInt_t i = 0; i < nJet; i++) {
    if (Jet_btagCSVV2[i] < 0.8484) continue;
    if (Jet_pt[i] < 20) continue;
    if (fabs(Jet_eta[i]) > 2.4) continue;
    
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
    
    if (hasOverLap(mom, Muons, 0.4)) continue;
    if (hasOverLap(mom, Elecs, 0.4)) continue;
    
    auto bjet = TParticle();
    bjet.SetMomentum(mom);
    b_bJet_tlv.push_back(mom);
    bJets.push_back(bjet);
    b_CSVv2.push_back(Jet_btagCSVV2[i]);
//    b_btagweight *= m_btagSF.eval_auto_bounds("central", BTagEntry::FLAV_B, Jet_eta[i], Jet_pt[i]);
//    b_btagweight_up *= m_btagSF.eval_auto_bounds("up", BTagEntry::FLAV_B, Jet_eta[i], Jet_pt[i]);
//    b_btagweight_dn *= m_btagSF.eval_auto_bounds("down", BTagEntry::FLAV_B, Jet_eta[i], Jet_pt[i]);
  }
  return bJets;
}

vector<TParticle> h2muAnalyser::nonbJetSelection(vector<TParticle> Muons, vector<TParticle> Elecs)
{
  vector<TParticle> nonbjets;
  for (UInt_t i = 0; i < nJet; i++) {
//    if (fabs(Jet_eta[i]) >= 2.4 &&  Jet_pt[i] < 30) continue;
//    if (fabs(Jet_eta[i]) < 2.4 &&  Jet_pt[i] < 20) continue;
    if (Jet_pt[i] < 30) continue; 
    if (fabs(Jet_eta[i]) > 4.7) continue; 
    if (Jet_jetId[i] < 1) continue;

    if (Jet_btagCSVV2[i] > 0.8484 && Jet_pt[i] > 20 && fabs(Jet_eta[i] < 2.4)) continue;

    TLorentzVector mom;
    mom.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
    
    if (hasOverLap(mom, Muons, 0.4)) continue;
    if (hasOverLap(mom, Elecs, 0.4)) continue;
    
    auto jet = TParticle();
    jet.SetMomentum(mom);
    b_nonbJet_tlv.push_back(mom);  
    nonbjets.push_back(jet);
  }
  return nonbjets;
}
bool h2muAnalyser::hasOverLap(TLorentzVector cand, vector<TParticle> objects, Float_t rad)
{
  for (auto obj: objects) {
    TLorentzVector mom;
    obj.Momentum(mom);
    if (cand.DeltaR(mom) < rad) {
      return true;
    }
  }
  return false;
}

Double_t h2muAnalyser::roccoR(TLorentzVector m, int &q, int &nGen, int &nTrackerLayers)
{
  Float_t u1 = gRandom->Rndm();
  Float_t u2 = gRandom->Rndm();
  if (!m_isMC) {
    return m_rocCor->kScaleDT(q, m.Pt(), m.Eta(), m.Phi(), 0, 0);
  } else {
    if (nGen > -1) {
      return m_rocCor->kScaleFromGenMC(q, m.Pt(), m.Eta(), m.Phi(),
				       nTrackerLayers, GenPart_pt[nGen],
				       u1, 0, 0);
    } else return m_rocCor->kScaleAndSmearMC(q, m.Pt(), m.Eta(), m.Phi(),
	        			     nTrackerLayers, u1, u2, 0, 0);
  }
}
