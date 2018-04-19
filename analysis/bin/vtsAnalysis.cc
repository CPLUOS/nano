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
  if(TMath::IsNaN(x)){
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


int main(Int_t argc, Char_t** argv){

  string hostDir = "root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/";

  if(argc <= 1){
    cout << "no input file is specified. running with default file." << endl;
    auto inFile = TFile::Open("/xrootd/store/group/nanoAOD/run2_2016v4/tsw/nanoAOD_1.root", "READ");
    auto inTree = (TTree*) inFile->Get("Events");
    vtsAnalysis ana(inTree);
    ana.MakeTree("nanotree.root");
    ana.Loop();
  }
  else{
    string jobName    = string(argv[1]);
    string sampleName = string(argv[2]);

    string outFileDir = hostDir+getenv("USER")+"/"+jobName+"/"+sampleName;
    for(Int_t i = 3; i < argc; i++){
      auto inFileName = argv[i];
      cout << inFileName << endl;
      TFile *inFile = TFile::Open(inFileName, "READ");
      TTree *inTree = (TTree*) inFile->Get("Events");
      vtsAnalysis ana(inTree);

      string outFileName = outFileDir+"/nanotree_"+to_string(i-3)+".root";
      ana.MakeTree(outFileName);
      ana.Loop();
    }
  }

}

void vtsAnalysis::Loop(){
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntries();

  // Events loop
  for (Long64_t iev=0; iev<nentries; iev++) {
    Long64_t entry = LoadTree(iev);
    fChain->GetEntry(entry);
    if (iev%10000 == 0) cout << iev << "/" << nentries << endl;

    ResetBranch();
    FindQuark();
    KshortAnalysis();
    EventSelection();
    outTree->Fill();
  }

  outFile->Write();
  outFile->Close();
  cout << endl;

}

void vtsAnalysis::MakeTree(std::string outFileName){
  outFile = TFile::Open(outFileName.c_str(), "RECREATE");
  outTree = new TTree("events", "events");
  outTree->Branch("channel", &b_channel, "channel/I");
  outTree->Branch("njet", &b_njet, "njet/I");
  outTree->Branch("met", &b_met, "met/F");
  outTree->Branch("dilep_tlv", "TLorentzVector", &b_dilep_tlv);

  outTree->Branch("KS_tlv", "TLorentzVector", &b_KS_tlv);
  outTree->Branch("isFrom_KS", &b_isFrom_KS, "isFrom_KS/I");
  outTree->Branch("d_KS", &b_d_KS , "d_KS/F" );
  outTree->Branch("x_KS", &b_x_KS, "x_KS/F");
  outTree->Branch("angleXY_KS", &b_angleXY_KS, "angleXY_KS/F");
  outTree->Branch("lxySig_KS", &b_lxySig_KS, "lxySig_KS/F");
  outTree->Branch("chi2_KS", &b_chi2_KS, "chi2_KS/F");
  outTree->Branch("dca_KS", &b_dca_KS, "dca_KS/F");
}

void vtsAnalysis::ResetBranch(){
  b_channel = -9; b_njet = -9;
  b_met = -99;
  b_dilep_tlv.SetPtEtaPhiM(0,0,0,0);
  recoleps.clear();

  b_KS_tlv.SetPtEtaPhiM(0,0,0,0);
  b_isFrom_KS = -9;
  b_d_KS = -1; b_x_KS = -1;
  b_angleXY_KS = -1; b_lxySig_KS = -1; b_chi2_KS = -1; b_dca_KS = -1;

  isS_ = false;
  q_ = -1; qb_ = -1;
  qj_ = -1; qbj_ = -1;
}

void vtsAnalysis::KshortAnalysis(){
  if (q_ == -1) return;

  TLorentzVector q_tlv;
  q_tlv.SetPtEtaPhiM(GenPart_pt[q_], GenPart_eta[q_], GenPart_phi[q_], GenPart_mass[q_]);
  TLorentzVector qb_tlv;
  qb_tlv.SetPtEtaPhiM(GenPart_pt[qb_], GenPart_eta[qb_], GenPart_phi[qb_], GenPart_mass[qb_]);
  for(unsigned int j=0; j<nJet;++j){
    TLorentzVector jet_tlv;
    jet_tlv.SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
    if ( q_tlv.DeltaR(jet_tlv) < 0.3) qj_ = j;
    if ( qb_tlv.DeltaR(jet_tlv) < 0.3) qbj_ = j;
  }
  if (qj_ == -1 || qbj_ == -1) return;
  //if (qj_ == -1 || qbj_ == -1) {cout << "No reco jet matched with gen s/b" << endl; return;}

  TLorentzVector qj_tlv;
  qj_tlv.SetPtEtaPhiM(Jet_pt[qj_], Jet_eta[qj_], Jet_phi[qj_], GenPart_mass[qj_]);
  TLorentzVector qbj_tlv;
  qbj_tlv.SetPtEtaPhiM(Jet_pt[qbj_], Jet_eta[qbj_], Jet_phi[qbj_], GenPart_mass[qbj_]);

  std::vector<struct HadStat> HadInJet;
  struct HadStat Stat;
  //Reconstructed hadrons
  for(unsigned int k=0; k<nhad; ++k){
    if (had_lxy[k]/had_lxyErr[k]>3 && had_angleXY[k]>0.95 && had_chi2[k]<5 &&had_dca[k]<2){ 
      TLorentzVector had_tlv;
      had_tlv.SetPtEtaPhiM(had_pt[k], had_eta[k], had_phi[k], had_mass[k]);
      if (abs(had_pdgId[k]) != 310) continue;
      //if (abs(had_pdgId[k]) != 310 && abs(had_pdgId[k]) != 3122) continue;
      auto dr = had_tlv.DeltaR(qj_tlv);
      auto drb = had_tlv.DeltaR(qbj_tlv);
      if(dr < 0.3 && (had_pt[k]/Jet_pt[qj_]) > 0.15){
        Stat.idx = k;
        Stat.pdgId = had_pdgId[k];
        Stat.x = had_pt[k]/Jet_pt[qj_];
        Stat.label = abs(GenPart_pdgId[q_]);
        HadInJet.push_back(Stat);
      }
      else if(drb < 0.3 && (had_pt[k]/Jet_pt[qbj_]) > 0.15){
        Stat.idx = k;
        Stat.pdgId = had_pdgId[k];
        Stat.x = had_pt[k]/Jet_pt[qbj_];
        Stat.label = abs(GenPart_pdgId[qb_]);
        HadInJet.push_back(Stat);
      }
    }
  }
  if (HadInJet.size() == 0) return;
  sort(HadInJet.begin(), HadInJet.end(), [](struct HadStat a, struct HadStat b){return a.x > b.x;}); // pick highest one
  auto idx = HadInJet[0].idx;
  auto pdgid = HadInJet[0].pdgId;
  b_KS_tlv.SetPtEtaPhiM(had_pt[idx], had_eta[idx], had_phi[idx], had_mass[idx]);

  if (abs(pdgid) == 310) {
    b_x_KS = HadInJet[0].x;
    b_isFrom_KS = HadInJet[0].label;
    b_d_KS = GetD(had_pt[idx], had_eta[idx], had_phi[idx], had_mass[idx], had_x[idx], had_y[idx], had_z[idx]);
    b_angleXY_KS = had_angleXY[idx];
    b_lxySig_KS = had_lxy[idx]/had_lxyErr[idx];
    b_chi2_KS = had_chi2[idx];
    b_dca_KS = had_dca[idx];
  }
}

void vtsAnalysis::FindQuark(){
  //Find s/b quark from Gen Info.
  for(unsigned int i=0; i<nGenPart; ++i){
    if (std::abs(GenPart_status[i] - 25) < 5 && GenPart_pdgId[i] == 3) { q_ = i; isS_ = true; }
    if (std::abs(GenPart_status[i] - 25) < 5 && GenPart_pdgId[i] == -3) { qb_ = i; }
    if (std::abs(GenPart_status[i] - 25) < 5 && GenPart_pdgId[i] == 5) { q_ = i; }
    if (std::abs(GenPart_status[i] - 25) < 5 && GenPart_pdgId[i] == -5) { qb_ = i; }
  }
}

void vtsAnalysis::EventSelection(){
  recolep1_tlv.SetPtEtaPhiM(0,0,0,0);
  recolep2_tlv.SetPtEtaPhiM(0,0,0,0);

  auto muons = muonSelection();
  auto elecs = elecSelection();
  if(muons.size()+ elecs.size() != 2) return;

  if      (muons.size() == 2){ recolep1 = muons[0]; recolep2 = muons[1]; b_channel = CH_MUMU; }
  else if (elecs.size() == 2){ recolep1 = elecs[0]; recolep2 = elecs[1]; b_channel = CH_ELEL; }
  else { recolep1 = muons[0]; recolep2 = elecs[0]; b_channel = CH_MUEL; }

  recolep1.Momentum(recolep1_tlv);
  recolep2.Momentum(recolep2_tlv);
  recoleps.push_back(recolep1_tlv);
  recoleps.push_back(recolep2_tlv);
  b_dilep_tlv = recolep1_tlv + recolep2_tlv;

  auto jets = jetSelection();

  if (b_dilep_tlv.M() < 20.) return;

  if (b_channel != CH_MUEL && 76 < b_dilep_tlv.M() && b_dilep_tlv.M() < 106) return;

  b_met = MET_pt;
  b_njet = jets.size();

  if (b_channel != CH_MUEL && b_met < 40) return;

  if (b_njet < 2) return;

}

vector<TParticle> vtsAnalysis::muonSelection(){
  vector<TParticle> muons;
  for (UInt_t i = 0; i < nMuon; ++i){
    if (!Muon_tightId[i]) continue;
    if (Muon_pt[i] < 20) continue;
    if (std::fabs(Muon_eta[i]) > 2.4) continue;
    if (Muon_pfRelIso04_all[i] > 0.15) continue;

    auto muon = GetTParticle(-13*Muon_charge[i], i);
    muons.push_back(muon);
  }
  return muons;
}

vector<TParticle> vtsAnalysis::elecSelection(){
  vector<TParticle> elecs;
  for (UInt_t i = 0; i < nElectron; ++i){
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

vector<TParticle> vtsAnalysis::jetSelection(){
  vector<TParticle> jets;
  for (UInt_t i = 0; i < nJet; ++i){
    if (Jet_pt[i] < 30) continue;
    if (std::fabs(Jet_eta[i]) > 2.4) continue;
    if (Jet_jetId[i] < 1) continue;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
    bool hasOverLap = false;
    for (auto lep : recoleps){ if (mom.DeltaR(lep) < 0.4) hasOverLap = true; }
    if (hasOverLap) continue;
    auto jet = TParticle();
    jet.SetMomentum(mom);
    jets.push_back(jet);
  }
  return jets;
}
