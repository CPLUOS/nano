#define vtsAnalysis_cxx
#include <iostream>
#include <vector>
#include <TLorentzVector.h>

#include "nano/analysis/src/nanoAnalysis.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;

int main(Int_t argc, Char_t** argv){

  if(argc <= 1){
    cout << "no input file is specified. running with default file." << endl;
    //auto inFile = TFile::Open("/cms/scratch/iwatson/tsW/cms/cmssw_ttswbw.root", "READ");
    //auto inTree = (TTree*) inFile->Get("Events");
    //vtsAnalysis ana(inTree);
    vtsAnalysis ana(0);
    ana.Loop();
  }
  else{
    for(Int_t i = 3; i < argc; i++){
      auto fileName = argv[i];
      TFile *inFile = TFile::Open(fileName, "READ");
      TTree *inTree = (TTree*) inFile->Get("Events");
      vtsAnalysis ana(inTree);
      ana.Loop();
    }
  }

}

void vtsAnalysis::Loop(){
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();

  outFile = TFile::Open("test.root", "RECREATE");
  h_nEvents = new TH1D("nEvents", "nEvents" ,1,0,1);
  outTree = MakeTree();

  // Events loop
  for (Long64_t iev=0; iev<nentries; iev++) {
    h_nEvents->Fill(0.5);
    Long64_t entry = LoadTree(iev);
    fChain->GetEntry(entry);
    if (iev%10000 == 0) cout << iev << "/" << nentries << endl;

    Analysis();
    outTree->Fill();
  }
  outFile->Write();
  outFile->Close();

}

TTree*  vtsAnalysis::MakeTree(){
  outTree = new TTree("events", "events");
  outTree->Branch("nMuon", &b_nMuon, "nMuon/I");
  outTree->Branch("Muon", "TLorentzVector", &b_Muon_tlv);
  return outTree;
}

void vtsAnalysis::ResetBranch(){
  b_nMuon = -9;
  b_Muon_tlv.SetPtEtaPhiM(0,0,0,0);
}

void vtsAnalysis::Analysis(){
  ResetBranch();

  b_nMuon = nMuon;
  if (nMuon <= 0) return;
  b_Muon_tlv.SetPtEtaPhiM(Muon_pt[0], Muon_eta[0], Muon_phi[0], Muon_mass[0]);
}
