#define vtsAnalysis_cxx
#include <iostream>
#include <vector>
#include <TLorentzVector.h>

#include "nano/analysis/src/nanoAnalysis.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;

TParticle vtsAnalysis::GetTParticle(int pdgId, int idx) {
  TLorentzVector tlv;
  if      (abs(pdgId) == 11) tlv.SetPtEtaPhiM(Electron_pt[idx], Electron_eta[idx], Electron_phi[idx], Electron_mass[idx]);
  else if (abs(pdgId) == 13) tlv.SetPtEtaPhiM(Muon_pt[idx], Muon_eta[idx], Muon_phi[idx], Muon_mass[idx]);
  else if (abs(pdgId) > 100) tlv.SetPtEtaPhiM(had_pt[idx], had_eta[idx], had_phi[idx], had_mass[idx]);
  else tlv.SetPtEtaPhiM(Jet_pt[idx], Jet_eta[idx], Jet_phi[idx], Jet_mass[idx]);
  return TParticle(pdgId, 0, 0, 0, 0, 0, tlv, tlv);
}

//Double_t vtsAnalysis::DeltaR(Double_t deta, Double_t dphi) {
//  return TMath::Sqrt(deta*deta + dphi*dphi);
//}

Double_t vtsAnalysis::DeltaPhi(Double_t phi1, Double_t phi2) {
  static const Double_t kPI = TMath::Pi();
  static const Double_t kTWOPI = 2*TMath::Pi();
  Double_t x = phi1 - phi2;
  if (TMath::IsNaN(x)) {
    std::cerr << "DeltaPhi function called with NaN" << std::endl;
    return x;
  }
  while (x >= kPI) x -= kTWOPI;
  while (x < -kPI) x += kTWOPI;
  return x;
}

Double_t vtsAnalysis::GetD(float pt, float eta, float phi, float m, float vx, float vy, float vz) {
  TLorentzVector tlv;
  tlv.SetPtEtaPhiM(pt,eta,phi,m);
  std::vector<Double_t> b = { tlv.Px()/(sqrt(tlv.Px()*tlv.Px()+tlv.Py()*tlv.Py()+tlv.Pz()*tlv.Pz())), tlv.Py()/(sqrt(tlv.Px()*tlv.Px()+tlv.Py()*tlv.Py()+tlv.Pz()*tlv.Pz())), tlv.Pz()/(sqrt(tlv.Px()*tlv.Px()+tlv.Py()*tlv.Py()+tlv.Pz()*tlv.Pz())) };
  std::vector<Double_t> a = { vx - PV_x, vy - PV_y, vz - PV_z };
  std::vector<Double_t> cross = { a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0] };
  return sqrt( cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2] );    
}

int main(Int_t argc, Char_t** argv) {

  string hostDir = "root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/";

  if (argc <= 1) {
    cout << "no input file is specified. running with default file." << endl;
    auto inFile = TFile::Open("/xrootd/store/group/nanoAOD/run2_2016v4/tsw/nanoAOD_1.root", "READ");
    auto inTree = (TTree*) inFile->Get("Events");
    cout << "c " << inTree << endl;
    vtsAnalysis ana(inTree);
    ana.MakeTree("nanotree.root");
    ana.Loop();
  }
  else{
    string jobName    = string(argv[1]);
    string sampleName = string(argv[2]);

    string outFileDir = hostDir+getenv("USER")+"/"+jobName+"/"+sampleName;
    for (Int_t i = 3; i < argc; i++) {
      auto inFileName = argv[i];
      TFile *inFile = TFile::Open(inFileName, "READ");
      TTree *inTree = (TTree*) inFile->Get("Events");
      vtsAnalysis ana(inTree);

      string outFileName = outFileDir+"/nanotree_"+to_string(i-3)+".root";
      ana.MakeTree(outFileName);
      ana.Loop();
    }
  }
}

void vtsAnalysis::Loop() {
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntries();

  // Events loop
  for (Long64_t iev=0; iev<nentries; iev++) {
    Long64_t entry = LoadTree(iev);
    fChain->GetEntry(entry);
    if (iev%10000 == 0) cout << iev << "/" << nentries << endl;

    ResetBranch();
    EventSelection();
    MatchingForMC();
    HadronAnalysis();
    outTree->Fill();
  }
//  for (int i =1; i<10; i++) { cout << h_cutFlow->GetBinContent(i) << "  "; }
  outFile->Write();
  outFile->Close();
  cout << endl;

}

void vtsAnalysis::MakeTree(std::string outFileName) {
  outFile = TFile::Open(outFileName.c_str(), "RECREATE");
  outTree = new TTree("events", "events");
  outTree->Branch("channel", &b_channel, "channel/I");
  outTree->Branch("njet", &b_njet, "njet/I");
  outTree->Branch("met", &b_met, "met/F");
  outTree->Branch("dilep_tlv", "TLorentzVector", &b_dilep_tlv);

  outTree->Branch("had_tlv", "TLorentzVector", &b_had_tlv);
  outTree->Branch("isFrom_had", &b_isFrom_had, "isFrom_had/I");
  outTree->Branch("d_had", &b_d_had , "d_had/F" );
  outTree->Branch("x_had", &b_x_had, "x_had/F");
  outTree->Branch("lxy_had", &b_lxy_had, "lxy_had/F");
  outTree->Branch("lxySig_had", &b_lxySig_had, "lxySig_had/F");
  outTree->Branch("angleXY_had", &b_angleXY_had, "angleXY_had/F");
  outTree->Branch("angleXYZ_had", &b_angleXYZ_had, "angleXYZ_had/F");
  outTree->Branch("chi2_had", &b_chi2_had, "chi2_had/F");
  outTree->Branch("dca_had", &b_dca_had, "dca_had/F");

  outTree->Branch("pt_had", &b_pt_had, "pt_had_had/F");
  outTree->Branch("eta_had", &b_eta_had, "eta_had/F");
  outTree->Branch("l3D_had", &b_l3D_had, "l3D_had/F");
  outTree->Branch("l3DSig_had", &b_l3DSig_had, "l3DSig_had/F");
  outTree->Branch("legDR_had", &b_legDR_had, "legDR_had/F");
  outTree->Branch("mass_had", &b_mass_had, "mass_had/F");
  outTree->Branch("pdgId_had", &b_pdgId_had, "pdgId_had/I");

  outTree->Branch("dau1_chi2_had", &b_dau1_chi2_had, "dau1_chi2_had/F");
  outTree->Branch("dau1_ipsigXY_had", &b_dau1_ipsigXY_had, "dau1_ipsigXY_had/F");
  outTree->Branch("dau1_ipsigZ_had", &b_dau1_ipsigZ_had, "dau1_ipsigZ_had/F");
  outTree->Branch("dau1_pt_had", &b_dau1_pt_had, "dau1_pt_had/F");

  outTree->Branch("dau2_chi2_had", &b_dau2_chi2_had, "dau2_chi2_had/F");
  outTree->Branch("dau2_ipsigXY_had", &b_dau2_ipsigXY_had, "dau2_ipsigXY_had/F");
  outTree->Branch("dau2_ipsigZ_had", &b_dau2_ipsigZ_had, "dau2_ipsigZ_had/F");
  outTree->Branch("dau2_pt_had", &b_dau2_pt_had, "dau2_pt_had/F");
}

void vtsAnalysis::ResetBranch() {
  b_channel = -9; b_njet = -9;
  b_met = -99;
  b_dilep_tlv.SetPtEtaPhiM(0,0,0,0);
  recoleps.clear();

  b_had_tlv.SetPtEtaPhiM(0,0,0,0);
  b_isFrom_had = -9;
  b_d_had = -1; b_x_had = -1;
  b_lxy_had = -1; b_lxySig_had = -1; b_angleXY_had = -1; b_angleXYZ_had = -1; b_chi2_had = -1; b_dca_had = -1;
  b_pt_had = -1; b_eta_had = -99; b_l3D_had = -1; b_l3DSig_had = -1; b_legDR_had = -1; b_mass_had = -99; b_pdgId_had = -99;
  b_dau1_chi2_had = -1; b_dau1_ipsigXY_had = -1; b_dau1_ipsigZ_had = -1; b_dau1_pt_had = -1;
  b_dau2_chi2_had = -1; b_dau2_ipsigXY_had = -1; b_dau2_ipsigZ_had = -1; b_dau2_pt_had = -1;


  qjMapForMC_.clear(); qMC_.clear(); genJet_.clear(); recoJet_.clear();
}


void vtsAnalysis::EventSelection() {
  //h_cutFlow->Fill(0);
  recolep1_tlv.SetPtEtaPhiM(0,0,0,0);
  recolep2_tlv.SetPtEtaPhiM(0,0,0,0);

  auto muons = muonSelection();
  auto elecs = elecSelection();
  if (muons.size()+ elecs.size() != 2) return;
  //h_cutFlow->Fill(1);

  if      (muons.size() == 2) { recolep1 = muons[0]; recolep2 = muons[1]; b_channel = CH_MUMU; }
  else if (elecs.size() == 2) { recolep1 = elecs[0]; recolep2 = elecs[1]; b_channel = CH_ELEL; }
  else { recolep1 = muons[0]; recolep2 = elecs[0]; b_channel = CH_MUEL; }

  recolep1.Momentum(recolep1_tlv);
  recolep2.Momentum(recolep2_tlv);
  recoleps.push_back(recolep1_tlv);
  recoleps.push_back(recolep2_tlv);
  b_dilep_tlv = recolep1_tlv + recolep2_tlv;

  auto jets = jetSelection();

  if (b_dilep_tlv.M() < 20.) return;
  //h_cutFlow->Fill(2);

  if (b_channel != CH_MUEL && 76 < b_dilep_tlv.M() && b_dilep_tlv.M() < 106) return;
  //h_cutFlow->Fill(3);

  b_met = MET_pt;
  b_njet = jets.size();

  if (b_channel != CH_MUEL && b_met < 40) return;
  //h_cutFlow->Fill(4);

  if (b_njet < 2) return;
  //h_cutFlow->Fill(5);

}

vector<TParticle> vtsAnalysis::muonSelection() {
  vector<TParticle> muons;
  for (UInt_t i = 0; i < nMuon; ++i) {
    if (!Muon_tightId[i]) continue;
    if (Muon_pt[i] < 20) continue;
    if (std::fabs(Muon_eta[i]) > 2.4) continue;
    if (Muon_pfRelIso04_all[i] > 0.15) continue;

    auto muon = GetTParticle(-13*Muon_charge[i], i);
    muons.push_back(muon);
  }
  return muons;
}

vector<TParticle> vtsAnalysis::elecSelection() {
  vector<TParticle> elecs;
  for (UInt_t i = 0; i < nElectron; ++i) {
    if (Electron_pt[i] < 20) continue;
    if (std::fabs(Electron_eta[i]) > 2.4) continue;
    if (Electron_cutBased[i] < 4) continue;
    float el_scEta = Electron_deltaEtaSC[i] + Electron_eta[i];
    if ( std::fabs(el_scEta) > 1.4442 &&  std::fabs(el_scEta) < 1.566 ) continue;

    auto elec = GetTParticle(-11*Electron_charge[i], i);
    elecs.push_back(elec);
  }
  return elecs;
}

vector<TParticle> vtsAnalysis::jetSelection() {
  vector<TParticle> jets;
  for (UInt_t i = 0; i < nJet; ++i) {
    if (Jet_pt[i] < 30) continue;
    if (std::fabs(Jet_eta[i]) > 2.4) continue;
    if (Jet_jetId[i] < 1) continue;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
    bool hasOverLap = false;
    for (auto lep : recoleps) { if (mom.DeltaR(lep) < 0.4) hasOverLap = true; }
    if (hasOverLap) continue;
    auto jet = TParticle();
    jet.SetMomentum(mom);
    jets.push_back(jet);
  }
  return jets;
}

void vtsAnalysis::MatchingForMC() {
  //Find s/b quark from Gen Info.  
  for (unsigned int i=0; i<nGenPart; ++i) {
    if (std::abs(GenPart_status[i] - 25) < 5 && abs(GenPart_pdgId[i]) == 3) {
      qMC_.push_back(i);
    }
    if (std::abs(GenPart_status[i] - 25) < 5 && abs(GenPart_pdgId[i]) == 5) {
      qMC_.push_back(i);
    }
  }
  if (qMC_.size() != 2) return;//(std::find(qMC_.begin(), qMC_.end(), -1) != qMC_.end()) return;

  auto q1 = qMC_[0];
  auto q2 = qMC_[1];
  TLorentzVector q1_tlv; 
  q1_tlv.SetPtEtaPhiM(GenPart_pt[q1], GenPart_eta[q1], GenPart_phi[q1], GenPart_mass[q1]);
  TLorentzVector q2_tlv; 
  q2_tlv.SetPtEtaPhiM(GenPart_pt[q2], GenPart_eta[q2], GenPart_phi[q2], GenPart_mass[q2]);

/*
  //Gen Particle & Gen Jet matching
  for (unsigned int j=0; j<nGenJet; ++j) {
    TLorentzVector gjet_tlv;
    gjet_tlv.SetPtEtaPhiM(GenJet_pt[j], GenJet_eta[j], GenJet_phi[j], GenJet_mass[j]);
    if (q1_tlv.DeltaR(gjet_tlv) < 0.3) { 
      if (abs(GenJet_partonFlavour[j] == 3) genJet_.push_back(j);
      if (abs(GenJet_partonFlavour[j] == 5) genJet_.push_back(j);
    }
    else if (q2_tlv.DeltaR(gjet_tlv) < 0.3) {
      if (abs(GenJet_partonFlavour[j] == 3) genJet_.push_back(j);
      if (abs(GenJet_partonFlavour[j] == 5) genJet_.push_back(j);
    }
  }
  if (genJet_.size() == 0) return;
*/

  //Gen Particle & Reco Jet matching
  for (unsigned int j=0; j<nJet;++j) {
    TLorentzVector jet_tlv;
    jet_tlv.SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
    if (q1_tlv.DeltaR(jet_tlv) < 0.3) {
      recoJet_.push_back(j);
      qjMapForMC_.insert({j, GenPart_pdgId[q1]});
    }
    else if (q2_tlv.DeltaR(jet_tlv) < 0.3) {
      recoJet_.push_back(j);
      qjMapForMC_.insert({j, GenPart_pdgId[q2]});
    }
  }
  if (recoJet_.size() == 0) return;

//  if (qMC_.size() != 2) std::cout << "qMC_ : " << qMC_.size() << " , idx : [ " << qMC_[0] << " , " << qMC_[1] << " ] " << std::endl;
//  if (recoJet_.size() != 2)std::cout << "reoJet_ : " << recoJet_.size() << " , idx : [ " << recoJet_[0] << " , " << recoJet_[1] << " ] " << std::endl;

  isMC_ = true;
}

void vtsAnalysis::HadronAnalysis() {
  std::vector<std::vector<struct HadStat>> JetCollection;
  std::vector<struct HadStat> HadInJet;
  struct HadStat Stat;
  unsigned int njets;
  if (isMC_) njets = recoJet_.size();
  else njets = nJet;
  //Reco Jet & Reco KS matching 
  for (unsigned int j=0; j<njets; ++j){
    HadInJet.clear();
    if(isMC_) j = recoJet_[j];
    TLorentzVector jet_tlv;
    jet_tlv.SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]); 
    for (unsigned int k=0; k<nhad; ++k) {
      if (abs(had_pdgId[k]) != 310 && abs(had_pdgId[k]) != 3122) continue;
      if ( (had_lxy[k]/had_lxyErr[k]) < 3 ) continue;
      if ( had_angleXY[k] < 0.95 ) continue;
      if ( had_chi2[k] > 5 ) continue;
      if ( had_dca[k] > 2 ) continue;
      //if ( !((had_lxy[k]/had_lxyErr[k]) > 15 && had_angleXY[k] > 0.98 && had_chi2[k] < 3 && had_dca[k] < 1) ) continue; // tight cut
      TLorentzVector had_tlv;
      had_tlv.SetPtEtaPhiM(had_pt[k], had_eta[k], had_phi[k], had_mass[k]);
      if (jet_tlv.DeltaR(had_tlv) < 0.3 && (had_pt[k]/Jet_pt[j]) > 0.15) {
        Stat.idx = k;
        Stat.pdgId = had_pdgId[k];
        Stat.x = had_pt[k]/Jet_pt[j];
        if (isMC_) Stat.label = qjMapForMC_[j];
        HadInJet.push_back(Stat);
      }
    }
    if (HadInJet.size() != 0 ) {
      if (HadInJet.size() > 1) sort(HadInJet.begin(), HadInJet.end(), [](struct HadStat a, struct HadStat b) {return a.x > b.x;}); // pick hadron with highest x in the jet
      JetCollection.push_back(HadInJet);
    }
  }
  if (JetCollection.size() != 0 ) {
    if (JetCollection.size() > 1) sort(JetCollection.begin(), JetCollection.end(), [](std::vector<struct HadStat> a, std::vector<struct HadStat> b) {return a[0].x > b[0].x;}); // pick jet-hadron pair with highest x
    auto idx = JetCollection[0][0].idx;
    b_had_tlv.SetPtEtaPhiM(had_pt[idx], had_eta[idx], had_phi[idx], had_mass[idx]);

    b_x_had = JetCollection[0][0].x;
    b_isFrom_had = JetCollection[0][0].label;
    b_d_had = GetD(had_pt[idx], had_eta[idx], had_phi[idx], had_mass[idx], had_x[idx], had_y[idx], had_z[idx]);
    b_lxy_had = had_lxy[idx];
    b_lxySig_had = had_lxy[idx]/had_lxyErr[idx];
    b_angleXY_had = had_angleXY[idx];
    b_angleXYZ_had = had_angleXY[idx];
    b_chi2_had = had_chi2[idx];
    b_dca_had = had_dca[idx];
    b_pt_had = had_pt[idx];
    b_eta_had = had_eta[idx];
    b_l3D_had = had_l3D[idx];
    b_l3DSig_had = had_l3D[idx]/had_l3DErr[idx];
    b_legDR_had = had_legDR[idx];
    b_mass_had = had_mass[idx];
    b_pdgId_had = had_pdgId[idx];
    b_dau1_chi2_had = had_dau1_chi2[idx]; 
    b_dau1_ipsigXY_had = had_dau1_ipsigXY[idx]; 
    b_dau1_ipsigZ_had = had_dau1_ipsigZ[idx]; 
    b_dau1_pt_had = had_dau1_pt[idx];
    b_dau2_chi2_had = had_dau2_chi2[idx]; 
    b_dau2_ipsigXY_had = had_dau1_ipsigXY[idx]; 
    b_dau2_ipsigZ_had = had_dau1_ipsigZ[idx]; 
    b_dau2_pt_had = had_dau1_pt[idx];
  }
}

