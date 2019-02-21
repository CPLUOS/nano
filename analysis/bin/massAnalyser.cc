#include "nano/analysis/interface/massAnalyser.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

massAnalyser::massAnalyser(TTree *tree, TTree *had, TTree *hadTruth, Bool_t isMC, Bool_t dl, Bool_t sle, Bool_t slm) : topEventSelectionDL(tree, had, hadTruth, isMC, dl, sle, slm)
{
}

massAnalyser::~massAnalyser() {
}

void massAnalyser::Loop() {
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntries();
  
  for (Long64_t iev=0; iev<nentries; iev++) {
    resetBranch();
    fChain->GetEntry(iev);
    int keep = EventSelection();
    if (keep >= 5) {
      collectTMVAvalues();
    }
    m_tree->Fill();
  }
}


int main(int argc, char* argv[]) {
  string env = getenv("CMSSW_BASE");
  string username = getenv("USER");

  if (argc <= 1) {
    TFile *f = TFile::Open("root://cms-xrdr.sdfarm.kr:1094///xrd/store/group/nanoAOD/run2_2016v4/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180430_152541/0000/nanoAOD_256.root", "read");
    TTree *tree;
    f->GetObject("Events", tree);

    massAnalyser t(tree, tree, 0, true, false, false, false);
    t.setOutput("test.root");
    t.Loop();
  }
  else {
    std::string strLine;
    
    std::string dirName = "root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/" + username + "/nanoAOD/" + std::string(argv[1]) + "/" + std::string(argv[2]);
    std::string temp = argv[2];
    
    Bool_t isDL = false;
    Size_t found_DL_1 = temp.find("Double");
    Size_t found_DL_2 = temp.find("MuonEG");
    if ((found_DL_1 != std::string::npos) || (found_DL_2 != std::string::npos) ) isDL = true;
    

    Bool_t isSL_e = false;
    Size_t found_SL_e = temp.find("SingleElectron");
    if (found_SL_e != std::string::npos) isSL_e = true;
    
    Bool_t isSL_m = false;
    Size_t found_SL_m = temp.find("SingleMuon");
    if (found_SL_m != std::string::npos) isSL_m = true;

    Bool_t isMC = false;
    Size_t found = temp.find("Run");
    if (found == std::string::npos) isMC = true;

    for (Int_t i = 3; i < argc; i++) {
      cerr << argv[i] << endl;
      TFile *f = TFile::Open(argv[i], "read");
      TTree *tree;                
      f->GetObject("Events", tree);
      
      temp = argv[i];   
      found = temp.find_last_of('/');
      std::string outPutName = dirName + temp.substr(found);
      massAnalyser t(tree, tree, 0, isMC, isDL, isSL_e, isSL_m);
      
      t.setOutput(outPutName);
      t.Loop();
    }
  }
  
  return 0;
}


void massAnalyser::setOutput(std::string outputName) {
  m_output = TFile::Open(outputName.c_str(), "recreate");
  m_tree = new TTree("event", "event");
  MakeBranch(m_tree);

  bdtg = new TMVA::Reader();
  //bdtg->AddVariable("cme_lxy", &b_cme_lxy);
  //bdtg->AddVariable("cme_lxyE", &b_cme_lxyE);
  //bdtg->AddVariable("cme_l3D", &b_cme_l3D);
  bdtg->AddVariable("cme_l3DE", &b_cme_l3DE);
  bdtg->AddVariable("cme_jetDR", &b_cme_jetDR);
  bdtg->AddVariable("cme_legDR", &b_cme_legDR);
  //bdtg->AddVariable("cme_nJet", &b_cme_nJet);
  bdtg->AddVariable("cme_dca", &b_cme_dca);
  //bdtg->AddVariable("cme_angleXY", &b_cme_angleXY);
  //bdtg->AddVariable("cme_angleXYZ", &b_cme_angleXYZ);
  //bdtg->AddVariable("cme_x", &b_cme_x);
  //bdtg->AddVariable("cme_y", &b_cme_y);
  //bdtg->AddVariable("cme_z", &b_cme_z);
  //bdtg->AddVariable("cme_pt", &b_cme_pt);
  bdtg->AddVariable("cme_chi2", &b_cme_chi2);
  //bdtg->AddVariable("cme_eta", &b_cme_eta);
  //bdtg->AddVariable("cme_phi", &b_cme_phi);
  //bdtg->AddVariable("cme_jet_btagCMVA", &b_cme_jet_btagCMVA);
  bdtg->AddVariable("cme_jet_btagCSVV2", &b_cme_jet_btagCSVV2);
  //bdtg->AddVariable("cme_jet_btagDeepB", &b_cme_jet_btagDeepB);
  //bdtg->AddVariable("cme_jet_btagDeepC", &b_cme_jet_btagDeepC);
  //bdtg->AddVariable("cme_dau1_chi2", &b_cme_dau1_chi2);
  //bdtg->AddVariable("cme_dau1_ipsigXY", &b_cme_dau1_ipsigXY);
  //bdtg->AddVariable("cme_dau1_ipsigZ", &b_cme_dau1_ipsigZ);
  //bdtg->AddVariable("cme_dau1_nHits", &b_cme_dau1_nHits);
  //bdtg->AddVariable("cme_dau1_pt", &b_cme_dau1_pt);
  //bdtg->AddVariable("cme_dau2_chi2", &b_cme_dau2_chi2);
  //bdtg->AddVariable("cme_dau2_ipsigXY", &b_cme_dau2_ipsigXY);
  //bdtg->AddVariable("cme_dau2_ipsigZ", &b_cme_dau2_ipsigZ);
  //bdtg->AddVariable("cme_dau2_nHits", &b_cme_dau2_nHits);
  //bdtg->AddVariable("cme_dau2_pt", &b_cme_dau2_pt);
  bdtg->AddSpectator("cme_mass", &b_cme_mass);
  bdtg->BookMVA("BDTG", "/home/jdj0715/nanoAOD/src/nano/analysis/test/topMass/cut/tmva/dataset/weights/TMVAClassification_BDTG13.weights.xml");
    
  h_nevents = new TH1D("nevents", "nevents", 1, 0, 1);
  h_genweights = new TH1D("genweight", "genweight", 1, 0, 1);
  h_weights = new TH1D("weight", "weight", 1, 0, 1);
  h_cutFlow = new TH1D("cutflow", "cutflow", 11, -0.5, 10.5);
}

void massAnalyser::MakeBranch(TTree* t) {
  t->Branch("nvertex", &b_nvertex, "nvertex/I");
  t->Branch("step", &b_step, "step/I");
  t->Branch("step1", &b_step1, "step1/O");
  t->Branch("step2", &b_step2, "step2/O");
  t->Branch("step3", &b_step3, "step3/O");
  t->Branch("step4", &b_step4, "step4/O");
  t->Branch("step5", &b_step5, "step5/O");
  t->Branch("channel", &b_channel, "channel/I");
  t->Branch("njet", &b_njet, "njet/I");
  t->Branch("nbjet", &b_nbjet, "nbjet/I");

  t->Branch("jet1", "TLorentzVector", &b_jet1);
  t->Branch("jet2", "TLorentzVector", &b_jet2);
  //t->Branch("bjet", "TLorentzVector", &b_bjet);
  t->Branch("lep1", "TLorentzVector", &b_lep1);
  //t->Branch("lep1_pid", &b_lep1_pid, "lep1_pid/I");    
  t->Branch("lep2", "TLorentzVector", &b_lep2);
  //t->Branch("lep2_pid", &b_lep2_pid, "lep2_pid/I");    
  t->Branch("dilep", "TLorentzVector", &b_dilep);
  t->Branch("tri", &b_tri, "tri/F");
  t->Branch("tri_up", &b_tri_up, "tri_up/F");
  t->Branch("tri_dn", &b_tri_dn, "tri_dn/F");
  t->Branch("met", &b_met, "met/F");
  //t->Branch("nhad", &nhad, "nhad/I");
  t->Branch("weight", &b_weight, "weight/F");
  t->Branch("puweight", &b_puweight, "puweight/F");
  t->Branch("genweight", &b_genweight, "genweight/F");
  t->Branch("csvweight", "std::vector<float>", &b_csvweights);
  t->Branch("btagweight", &b_btagweight, "btagweight/F");
  t->Branch("btagweight_up", &b_btagweight_up, "btagweight_up/F");
  t->Branch("btagweight_dn", &b_btagweight_dn, "btagweight_dn/F");
  t->Branch("mueffweight", &b_mueffweight, "mueffweight/F");
  t->Branch("eleffweight", &b_eleffweight, "eleffweight/F");
  t->Branch("PV_npvs", &PV_npvs, "PV_npvs/I");
  t->Branch("trig_m", &b_trig_m, "trig_m/O");
  t->Branch("trig_e", &b_trig_e, "trig_e/O");
  t->Branch("trig_mm", &b_trig_mm, "trig_mm/O");
  t->Branch("trig_em", &b_trig_em, "trig_em/O");
  t->Branch("trig_ee", &b_trig_ee, "trig_ee/O");
  
  t->Branch("cme_tmva_bdtg", &b_cme_tmva_bdtg, "cme_tmva_bdtg/F");
  t->Branch("cme_diffMass", &b_cme_diffMass, "cme_diffMass/F");
  t->Branch("cme_pdgId", &b_cme_pdgId, "cme_pdgId/I");
  
  t->Branch("d0","TLorentzVector",&b_d0);
  t->Branch("vecSumLep1AndD0","std::vector<float>",&b_vecSumLep1AndD0);
  t->Branch("vecSumLep2AndD0","std::vector<float>",&b_vecSumLep2AndD0);
  t->Branch("d0_lepSV_dRM","std::vector<float>",&b_d0_lepSV_dRM);
  t->Branch("d0_lepSV_dRM_wrong","std::vector<float>",&b_d0_lepSV_dRM_wrong);
  t->Branch("d0_lepSV_correctM","std::vector<float>",&b_d0_lepSV_correctM);
  
  t->Branch("cme_x", &b_cme_x, "cme_x/F");
  t->Branch("cme_y", &b_cme_y, "cme_y/F");
  t->Branch("cme_z", &b_cme_z, "cme_z/F");
}

void massAnalyser::resetBranch() {
  Reset();
  b_d0.Clear();
  d0s.clear();

  b_cme_mass = -999;
  b_cme_diffMass = -999;
  b_cme_pdgId = 0;
  b_cme_tmva_bdtg = -999;
  b_bdtg = -1; b_maxbIdx = -1;
  b_cme_jetDR = -999;
  b_cme_legDR = -999;
  b_cme_dca = -999;
  b_cme_chi2 = -999;
  b_cme_l3DE = -999;
  b_cme_jet_btagCSVV2 = -999;


  b_cme_x = -999;
  b_cme_y = -999;
  b_cme_z = -999;
  
  b_d0_lepSV_dRM.clear();
  b_d0_lepSV_dRM_wrong.clear();
  b_d0_lepSV_correctM.clear();

  b_vecSumLep1AndD0.clear();
  b_vecSumLep2AndD0.clear();
}

void massAnalyser::collectTMVAvalues() {
  for (UInt_t i=0; i < nhad; i++) {
    if (fabs(had_pdgId[i]) != 421) continue;
    if (had_pt[i] < 5) continue;
    if ((had_l3D[i]/had_l3DErr[i]) > 200) continue;
    if (fabs(had_jet_btagCSVV2[i]) < 0.05) continue;
    if (fabs(had_lxy[i]) < 0.1) continue;
    if (fabs(had_l3D[i]) < 0.2) continue;
    if (fabs(had_angleXY[i]) < 0.95) continue;
    if (fabs(had_dca[i]) > 1) continue;
    if (fabs(had_x[i]) > 6) continue;
    if (fabs(had_y[i]) > 6) continue;
    if (fabs(had_z[i]) > 18) continue;
    if (fabs(had_dau1_chi2[i]) > 5) continue;
    if (fabs(had_dau2_chi2[i]) > 5) continue;
    if (had_l3D[i] > 0 && had_l3DErr[i] > 0) {
      b_cme_l3DE = had_l3D[i] / had_l3DErr[i];
    }
    if (fabs(had_jetDR[i]) > 0.3) continue;
    if (fabs(had_legDR[i]) > 0.6) continue;
    if (b_cme_x == -999) continue;
    if (b_cme_y == -999) continue;
    if (b_cme_z == -999) continue;
    
    b_cme_jetDR = had_jetDR[i];
    b_cme_legDR = had_legDR[i];
    b_cme_dca = had_dca[i];
    b_cme_x = had_x[i];
    b_cme_y = had_y[i];
    b_cme_z = had_z[i];
    b_cme_chi2 = had_chi2[i];
    b_cme_jet_btagCSVV2 = had_jet_btagCSVV2[i];
    b_cme_mass = had_mass[i];
    b_cme_pdgId = had_pdgId[i];
    b_cme_tmva_bdtg = bdtg->EvaluateMVA("BDTG");
    if (b_cme_tmva_bdtg > b_bdtg) {
      b_maxbIdx = i;
      b_bdtg = b_cme_tmva_bdtg;
      b_cme_pdgId = had_pdgId[b_maxbIdx];
      b_cme_tmva_bdtg = b_bdtg;
      b_cme_diffMass = had_diffMass[b_maxbIdx];
    }
  }
  if (had_mass[b_maxbIdx] == 0) return;
  TLorentzVector d0_tlv;
  d0_tlv.SetPtEtaPhiM(had_pt[b_maxbIdx], had_eta[b_maxbIdx], had_phi[b_maxbIdx], had_mass[b_maxbIdx]);
  d0s.push_back(d0_tlv);
  b_d0 = d0s[0];
 
  TLorentzVector vecSumLep1AndD0, vecSumLep2AndD0;
  vecSumLep1AndD0 = b_lep1 + b_d0;
  vecSumLep2AndD0 = b_lep2 + b_d0;
  
  float deltaEta_1 = b_lep1.Eta()-b_d0.Eta();
  float deltaPhi_1 = b_lep1.Phi()-b_d0.Phi();
  float deltaEta_2 = b_lep2.Eta()-b_d0.Eta();
  float deltaPhi_2 = b_lep2.Phi()-b_d0.Phi();

  float distance1 = sqrt((deltaEta_1*deltaEta_1) + (deltaPhi_1*deltaPhi_1));
  float distance2 = sqrt((deltaEta_2*deltaEta_2) + (deltaPhi_2*deltaPhi_2));
  
  if (distance1 < distance2) {
      b_vecSumLep1AndD0.push_back(vecSumLep1AndD0.M());
  }
  else {
      b_vecSumLep2AndD0.push_back(vecSumLep2AndD0.M());
  }

  b_d0_lepSV_dRM.push_back(( distance1 < distance2 ? vecSumLep1AndD0.M() : vecSumLep2AndD0.M() ));
  b_d0_lepSV_dRM_wrong.push_back(( distance1 > distance2 ? vecSumLep1AndD0.M() : vecSumLep2AndD0.M() ));
  b_d0_lepSV_correctM.push_back(( b_lep1_pid * b_cme_pdgId < 0 ? vecSumLep1AndD0.M() : vecSumLep2AndD0.M() ));
}

