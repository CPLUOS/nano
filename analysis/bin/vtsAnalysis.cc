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
    else tlv.SetPtEtaPhiM(Jet_pt[idx], Jet_eta[idx], Jet_phi[idx], Jet_mass[idx]);

    return TParticle(pdgId, 0, 0, 0, 0, 0, tlv, tlv);
}

int main(Int_t argc, Char_t** argv){

  string hostDir = "root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/";

  if(argc <= 1){
    cout << "no input file is specified. running with default file." << endl;
    auto inFile = TFile::Open("nanoAOD_forsync.root", "READ");
    //auto inFile = TFile::Open("/xrootd/store/user/tt8888tt/nanoAOD_forsync.root", "READ");
    auto inTree = (TTree*) inFile->Get("Events");
    cout << "c " << inTree << endl;
    vtsAnalysis ana(inTree);
    //vtsAnalysis ana(0);
    ana.MakeTree("nanotree.root");
    ana.Loop();
  }
  else{
    string jobName    = string(argv[1]);
    string sampleName = string(argv[2]);

    string outFileDir = hostDir+getenv("USER")+"/"+jobName+"/"+sampleName;
    for(Int_t i = 3; i < argc; i++){
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

void vtsAnalysis::Loop(){
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntries();

  // Events loop
  for (Long64_t iev=0; iev<nentries; iev++) {
    Long64_t entry = LoadTree(iev);
    fChain->GetEntry(entry);
    if (iev%10000 == 0) cout << iev << "/" << nentries << endl;

    //ResetBranch();
    //EventSelection();
    //outTree->Fill();
  }

  for (int i =1; i<10; i++){ cout << h_cutFlow->GetBinContent(i) << "  "; }
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
}

void vtsAnalysis::ResetBranch(){
  b_channel = -9; b_njet = -9;
  b_met = -99;
  b_dilep_tlv.SetPtEtaPhiM(0,0,0,0);
  recoleps.clear();

}

void vtsAnalysis::EventSelection(){
  h_cutFlow->Fill(0);
  recolep1_tlv.SetPtEtaPhiM(0,0,0,0);
  recolep2_tlv.SetPtEtaPhiM(0,0,0,0);

  auto muons = muonSelection();
  auto elecs = elecSelection();
  if(muons.size()+ elecs.size() != 2) return;
  h_cutFlow->Fill(1);

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
  h_cutFlow->Fill(2);

  if (b_channel != CH_MUEL && 76 < b_dilep_tlv.M() && b_dilep_tlv.M() < 106) return;
  h_cutFlow->Fill(3);

  b_met = MET_pt;
  b_njet = jets.size();

  if (b_channel != CH_MUEL && b_met < 40) return;
  h_cutFlow->Fill(4);

  if (b_njet < 2) return;
  h_cutFlow->Fill(5);

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
