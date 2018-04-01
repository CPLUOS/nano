#include <iostream>
#include <vector>
#include <TLorentzVector.h>

#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
//#include "classes/DelphesClasses.h"

using namespace std;

int main(){

  // Reading input
  auto inFile = TFile::Open("/cms/scratch/iwatson/tsW/cms/cmssw_ttswbw.root", "READ");
  auto inTree = (TTree*) inFile->Get("Events");

  map<TString,Int_t[500]>   mi;
  map<TString,Float_t[500]> mf;
  map<TString,Bool_t[500]>  mb;
  map<TString,UChar_t[500]> mc;
  map<TString,UInt_t>  mu; map<TString,ULong64_t>  ml;

  auto branches = inTree->GetListOfBranches();
  for (int i = 0; i < branches->GetEntries(); i++) {
    TBranch *br = (TBranch*) branches->At(i);
    TString name = br->GetName();
    if (name.Contains("HLT")) continue;
    
    auto typeName = string(br->GetLeaf(name)->GetTypeName());
    if      (typeName == string("Int_t"))     { inTree->SetBranchAddress(name, mi["b_"+name]); }
    else if (typeName == string("Float_t"))   { inTree->SetBranchAddress(name, mf["b_"+name]); }
    else if (typeName == string("Bool_t"))    { inTree->SetBranchAddress(name, mb["b_"+name]); }
    else if (typeName == string("UChar_t"))   { inTree->SetBranchAddress(name, mc["b_"+name]); }
    else if (typeName == string("UInt_t"))    { inTree->SetBranchAddress(name, &mu["b_"+name]); }
    else if (typeName == string("ULong64_t")) { inTree->SetBranchAddress(name, &ml["b_"+name]); }
    else std::cout << typeName << " is not avalible" << std::endl;
  }
  
  // Decalre out tree
  auto outFile = TFile::Open("out.root", "RECREATE");
  auto outTree = new TTree("tsW","tsW");
  float gen_Ks_pt, Ks_pt;
  int Ks_matchJet;
  outTree->Branch("gen_Ks_pt", &gen_Ks_pt, "gen_Ks_pt/F");
  outTree->Branch("Ks_pt", &Ks_pt, "Ks_pt/F");
  outTree->Branch("Ks_matchJet", &Ks_matchJet, "Ks_matchJet/I");

  // Events loop
  for (auto iev = 0; iev < inTree->GetEntries(); ++iev){
     inTree->GetEntry(iev);
     if (iev%1000 == 0) cout << iev << "/" << inTree->GetEntries() << endl;
     if (iev==10000) break;

     //gen loop
     for (unsigned int ig=0; ig<mu["b_nGenPart"]; ++ig) {
        if ( abs(mi["b_GenPart_pdgId"][ig]) != 310) continue;
        cout << mi["b_GenPart_pdgId"][ig] << endl;
        gen_Ks_pt = mf["b_GenPart_pt"][ig];
     }

     //reco loop
     for (unsigned int im=0; im<mu["b_nmeson"]; ++im) {
        if ( abs(mi["b_meson_pdgId"][im]) != 310) continue;
        Ks_pt = mf["b_meson_pt"][im];
     }
     outTree->Fill();
  }

  inFile->Close();

  outFile->Write();
  outFile->Close();
}
