#include "nano/analysis/interface/slmassAnalyser.h"
#include "nano/analysis/interface/hadAnalyser.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <cstdlib>
using namespace std;

slmassAnalyser::slmassAnalyser(TTree *tree, TTree *had, TTree *hadTruth, Bool_t isMC,  Bool_t sle, Bool_t slm) : topEventSelectionSL(tree, had, hadTruth, isMC, sle, slm)
{
}


slmassAnalyser::~slmassAnalyser()
{
}


void slmassAnalyser::Loop() {
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntries();

  // Events loop
  for (Long64_t iev=0; iev<nentries; iev++) {
    Reset();
    resetBranch();
    fChain->GetEntry(iev);
    int keep = EventSelection();
    cmesonSelection();
    if (keep != 0) {
      collectTMVAvalues();
      m_tree->Fill();
    }
  }
}
int main(int argc, char* argv[]) {
  string env = getenv("CMSSW_BASE");
  string username = getenv("USER");

  if (argc != 1) {
    std::string dirName = "root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/" + username + "/nanoAOD/" + std::string(argv[1]) + "/" + std::string(argv[2]);
    std::string temp = argv[2];
    
    Bool_t isSL_e = false;
    Size_t found_SL_e = temp.find("SingleElectron");
    if (found_SL_e != std::string::npos) isSL_e = true;
    
    Bool_t isSL_m = false;
    Size_t found_SL_m = temp.find("SingleMuon");
    if (found_SL_m != std::string::npos) isSL_m = true;

    Bool_t isMC = false;
    Size_t found = temp.find("Run");
    if (found == std::string::npos) isMC = true;

    for(Int_t i = 3; i < argc; i++) {
      TFile *f = TFile::Open(argv[i], "read");
      TTree *tree;                  
      f->GetObject("Events", tree);
      
      temp = argv[i];   
      found = temp.find_last_of('/');
      std::string outPutName = dirName + temp.substr(found);
      slmassAnalyser t(tree, tree, 0, isMC, isSL_e, isSL_m);
      
      t.setOutput(outPutName);
      t.Loop();
    }
  }
  else {
    TFile *f = TFile::Open("/xrootd/store/group/nanoAOD/run2_2016v4/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180430_152541/0000/nanoAOD_256.root", "read");
    //TFile *f = TFile::Open("/cms/scratch/jdj0715/nanoAOD/src/nano/nanoAOD/prod/nanoAOD.root", "read");
    TTree *tree;
    f->GetObject("Events", tree);

    slmassAnalyser t(tree, tree, 0, true, false, false);
    t.setOutput("sltest.root");
    t.Loop();
  }
  return 0;
}

void slmassAnalyser::setOutput(std::string outputName) {
  m_output = TFile::Open(outputName.c_str(), "recreate");
  m_tree = new TTree("event", "event");
  MakeBranch(m_tree);
    
  h_nevents = new TH1D("nevents", "nevents", 1, 0, 1);
  h_genweights = new TH1D("genweight", "genweight", 1, 0, 1);
  h_weights = new TH1D("weight", "weight", 1, 0, 1);
  h_cutFlow = new TH1D("cutflow", "cutflow", 11, -0.5, 10.5);
  h_cutFlowEl = new TH1D("cutflowEl", "cutflow", 11, -0.5, 10.5);
  h_cutFlowMu = new TH1D("cutflowMu", "cutflow", 11, -0.5, 10.5);
}


void slmassAnalyser::MakeBranch(TTree* t) {
  t->Branch("nvertex", &b_nvertex, "nvertex/I");
  t->Branch("step", &b_step, "step/I");
  t->Branch("channel", &b_channel, "channel/I");
  t->Branch("njet", &b_njet, "njet/I");
  t->Branch("nbjet", &b_nbjet, "nbjet/I");
  m_tree->Branch("lep", "TLorentzVector", &b_lep); 
  t->Branch("had", "TLorentzVector", &b_had);
  t->Branch("tri", &b_tri, "tri/F");
  t->Branch("tri_up", &b_tri_up, "tri_up/F");
  t->Branch("tri_dn", &b_tri_dn, "tri_dn/F");
  t->Branch("met", &b_met, "met/F");
  t->Branch("nhad", &nhad, "nhad/i");
  t->Branch("weight", &b_weight, "weight/F");
  t->Branch("puweight", &b_puweight, "puweight/F");
  t->Branch("genweight", &b_genweight, "genweight/F");
  t->Branch("csvweight", "std::vector<float>", &b_csvweights);
  t->Branch("btagweight", &b_btagweight, "btagweight/F");
  t->Branch("mueffweight", &b_mueffweight, "mueffweight/F");
  t->Branch("eleffweight", &b_eleffweight, "eleffweight/F");
  t->Branch("PV_npvs", &PV_npvs, "PV_npvs/I");
  t->Branch("trig_m", &b_trig_m, "trig_m/O");
  t->Branch("trig_e", &b_trig_e, "trig_e/O");
  
  t->Branch("had_vecSumDMLep", "TLoremtzVector", &b_had_vecSumDMLep);
  t->Branch("had_sumM", &b_had_sumM, "had_sumM/F");

  t->Branch("cme_mass", &b_cme_mass, "cme_mass/F");
  t->Branch("cme_lxy", &b_cme_lxy, "cme_lxy/F");
  t->Branch("cme_lxyE", &b_cme_lxyE, "cme_lxyE/F");
  t->Branch("cme_x", &b_cme_x, "cme_x/F");
  t->Branch("cme_y", &b_cme_y, "cme_y/F");
  t->Branch("cme_z", &b_cme_z, "cme_z/F");
  t->Branch("cme_pt", &b_cme_pt, "cme_pt/F");
  t->Branch("cme_pdgId", &b_cme_pdgId, "cme_pdgId/I");
}


void slmassAnalyser::resetBranch() {
    hads.clear();
    b_cme_pdgId = 0;
    b_cme_mass = -999;
    b_cme_lxy = -999;
    b_cme_lxyE = -999;
    b_cme_x = -999;
    b_cme_y = -999;
    b_cme_z = -999;
    b_cme_pt = -999;
}


void slmassAnalyser::collectTMVAvalues() {
    for (UInt_t i = 0; i <nhad; ++i) {

        b_cme_lxy = had_lxy[i];
        b_cme_lxyE = had_lxy[i] / had_lxyErr[i];
        b_cme_l3D = had_l3D[i];
        b_cme_l3DE = had_l3D[i] / had_l3DErr[i];
        b_cme_jetDR = had_jetDR[i];
        b_cme_legDR = had_legDR[i];
        b_cme_dca = had_dca[i];
        b_cme_angleXY = had_angleXY[i];
        b_cme_angleXYZ = had_angleXYZ[i];
        b_cme_x = had_x[i];
        b_cme_y = had_y[i];
        b_cme_z = had_z[i];
        b_cme_pt = had_pt[i];
        b_cme_chi2 = had_chi2[i];
        b_cme_eta = had_eta[i];
        b_cme_phi = had_phi[i];
        b_cme_jet_btagCMVA = had_jet_btagCMVA[i];
        b_cme_jet_btagCSVV2 = had_jet_btagCSVV2[i];
        b_cme_jet_btagDeepB = had_jet_btagDeepB[i];
        b_cme_jet_btagDeepC = had_jet_btagDeepC[i];
        b_cme_dau1_chi2 = had_dau1_chi2[i];
        b_cme_dau1_ipsigXY = had_dau1_ipsigXY[i];
        b_cme_dau1_ipsigZ = had_dau1_ipsigZ[i];
        b_cme_dau1_nHits = had_dau1_nHits[i];
        b_cme_dau1_pt = had_dau1_pt[i];
        b_cme_dau2_chi2 = had_dau2_chi2[i];
        b_cme_dau2_ipsigXY = had_dau2_ipsigXY[i];
        b_cme_dau2_ipsigZ = had_dau2_ipsigZ[i];
        b_cme_dau2_nHits = had_dau2_nHits[i];
        b_cme_dau2_pt = had_dau2_pt[i];
        b_cme_mass = had_mass[i];
        b_cme_pdgId = had_pdgId[i]; 

    }
}

void slmassAnalyser::cmesonSelection() {
    if (nhad < 1) return;

    for (UInt_t i = 0; i <nhad; ++i) {
        TLorentzVector had_tlv;
        had_tlv.SetPtEtaPhiM(had_pt[i], had_eta[i], had_phi[i], had_mass[i]); 
        hads.push_back(had_tlv);
        if (hads.size() < 1) continue;
        sort(hads.begin(), hads.end(), [](const TLorentzVector& a, const TLorentzVector& b){return a.Pt() > b.Pt();});
        hads.erase(hads.begin()+1, hads.end());
        b_had = hads[0];
        
        b_had_vecSumDMLep = b_lep + b_had;
        b_had_sumM = b_had_vecSumDMLep.M();
        

    }
}


