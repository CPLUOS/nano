#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <TString.h>
#include <TSystem.h>
#include <TRefArray.h>
#include <TMatrixDfwd.h>
#include <TVectorD.h>
#include <TChain.h>

#include "nano/analysis/interface/doubleHiggsAnalyser.h"
#include "classes/DelphesClasses.h"
#include "nano/analysis/interface/lester_mt2_bisect.h"

#include "nano/oxbridgekinetics/src/Mt2/Basic_Mt2_332_Calculator.h"
#include "nano/oxbridgekinetics/src/Mt2/ChengHanBisect_Mt2_332_Calculator.h"

using namespace std;

// return mother particle(GenParticle) of particle 'p' among 'particles'.
GenParticle* getMother(TClonesArray* particles, GenParticle* p){
  if (p->M1==-1) return nullptr;
  auto mom = static_cast<GenParticle *>(particles->At(p->M1));
  while (mom->PID == p->PID){
    mom = static_cast<GenParticle *>(particles->At(mom->M1));
    if (mom->M1==-1) return nullptr;
  }
  return mom;
}

// return the PID of the mother particle of the particle at 'ip' among 'particles'
int isFrom(TClonesArray* particles, int ip){
   auto p = static_cast<GenParticle *>(particles->At(ip));
   // check if it's from Higgs
   auto mom = getMother(particles, p); 
   if (mom==nullptr) return 0;
   auto grmom = getMother(particles, mom);
   if (grmom==nullptr) return 0;
   return grmom->PID;
}

void doubleHiggsAnalyser::MakeOutputBranch(TTree *tree) {
  // MT2 variables
  tree->Branch("lester_MT2",&lester_MT2,"lester_MT2/F");
  tree->Branch("basic_MT2_332",&basic_MT2_332,"basic_MT2_332/F");
  tree->Branch("ch_bisect_MT2_332",&ch_bisect_MT2_332,"ch_bisect_MT2_332/F");

  tree->Branch("lepton1_pt",&lepton1_pt,"lepton1_pt/F");
  tree->Branch("lepton2_pt",&lepton2_pt,"lepton1_pt/F");
  tree->Branch("missing_et",&missing_et,"missing_et/F");
  tree->Branch("mt",&mt,"mt/F");

  // truth matching variables
  tree->Branch("fromHiggs",&fromHiggs,"fromHiggs/I");
  tree->Branch("fromTop",&fromTop,"fromTop/I");
  tree->Branch("fromZ",&fromZ,"fromZ/I");

  tree->Branch("lepton1MotherPID",&from1,"from1/I");
  tree->Branch("lepton2MotherPID",&from2,"from2/I");
}

void doubleHiggsAnalyser::SetOutput(TString output_file_name) {
  out_file = TFile::Open(output_file_name,"RECREATE");
  out_tree = new TTree("events","events"); 
}

void doubleHiggsAnalyser::SetBranchAddress() {
  if (delphes_flag) doubleHiggsAnalyser::SetDelphesBranchAddress();
  else doubleHiggsAnalyser::SetNanoBranchAddress();
}

void doubleHiggsAnalyser::SetNanoBranchAddress() {

}

void doubleHiggsAnalyser::SetDelphesBranchAddress() {
  del_tree->SetBranchAddress("Particle",&particles);
  del_tree->SetBranchAddress("MissingET",&missings);
  del_tree->SetBranchAddress("Jet",&jets);
}

void doubleHiggsAnalyser::Initiate(TString output_file_name) {
  // set output file
  doubleHiggsAnalyser::SetOutput(output_file_name);
  // make output branch
  doubleHiggsAnalyser::MakeOutputBranch(out_tree);
  // set branch address
  doubleHiggsAnalyser::SetBranchAddress();
}

void doubleHiggsAnalyser::ResetVariables() {
  if (delphes_flag) ResetDelphesVariables();
  else ResetNanoVariables();

  lester_MT2 = -99;
  basic_MT2_332 = -99;
  ch_bisect_MT2_332 = -99;
  lepton1_pt = -99;
  lepton2_pt = -99;
  missing_et = -99;
  mt = -99;
  fromHiggs = 0;
  fromTop = 0;
  fromZ = 0;
  from1 = 0;
  from2 = 0;
}

void doubleHiggsAnalyser::ResetNanoVariables() {

}

void doubleHiggsAnalyser::ResetDelphesVariables() {
  muons.clear();
}

bool doubleHiggsAnalyser::Analysis() {
  if (delphes_flag) return DelphesAnalysis();
  else return NanoAnalysis();
}

bool doubleHiggsAnalyser::NanoAnalysis() {
  return false;
}

bool doubleHiggsAnalyser::DelphesAnalysis() {
  //map<Float_t, int, greater<Float_t>> muons : map of <pt,index>:<K,V> of muon sorted by pt.
    
    auto m = static_cast<const MissingET *>(missings->At(0)); // There is always one MET object.
    if (m->MET<20) return false;
    missing.SetPtEtaPhiM(m->MET,m->Eta,m->Phi,0.19);
 
    bool jetflag;
    for (int ij = 0; ij < jets->GetEntries(); ij++){
      auto jet = static_cast<const Jet *>(jets->At(ij));
      if (abs(jet->Flavor) != 5) continue;
      else {jetflag = true; break;}
    }
    if (!jetflag) return false;
    
    
    for (int ip = 0; ip < particles->GetEntries(); ip++){
      auto p = static_cast<const GenParticle *>(particles->At(ip));
      if (abs(p->PID)!=13 || p->PT<20 || fabs(p->Eta)>2.5 || p->M1==-1) continue;
      muons.insert(make_pair(p->PT,ip));
    }
  
    if (muons.size()<2) {
      return false;
    }
    
    map<Float_t,int,greater<Float_t>>::iterator iter = muons.begin();
    auto mu1 = static_cast<const GenParticle *>(particles->At(iter->second));
    int from1 = isFrom(particles, iter->second);
    if (from1==25) fromHiggs++;
    if (from1==6) fromTop++;
    if (from1==32) fromZ++;
    muon1.SetPtEtaPhiM(mu1->PT,mu1->Eta,mu1->Phi,mu1->Mass);
    ++iter;
    auto mu2 = static_cast<const GenParticle *>(particles->At(iter->second));
    int from2 = isFrom(particles, iter->second);
    if (from2==25) fromHiggs++;
    if (from2==6) fromTop++;
    if (from2==32) fromZ++;
    muon2.SetPtEtaPhiM(mu2->PT,mu2->Eta,mu2->Phi,mu2->Mass);
    
    lepton1_pt = mu1->PT;
    lepton2_pt = mu2->PT;
    missing_et = m->MET;
    auto muons = muon1+muon2;
    mt = muons.M()*muons.M() + muons.Px()*muons.Px() + muons.Py()*muons.Py();

    // get MT2
    Mt2::LorentzTransverseVector vis_A(Mt2::TwoVector(muon1.Px(), muon1.Py()), muon1.M());
    Mt2::LorentzTransverseVector vis_B(Mt2::TwoVector(muon2.Px(),muon2.Py()), muon2.M());
    Mt2::TwoVector pT_Miss(missing.Px(), missing.Py());
    
    lester_MT2 = asymm_mt2_lester_bisect::get_mT2( // the simpliest way to calculate MT2
              muon1.M(),muon1.Px(),muon1.Py(),
              muon2.M(),muon2.Px(),muon2.Px(),
              missing.Px(),missing.Py(),
              missing.M(),missing.M());
    basic_MT2_332 = basic_mt2_332Calculator.mt2_332(vis_A, vis_B, pT_Miss, missing.M());
    if (basic_MT2_332 < 0.4) cout << "basic_MT2_332 is 0" << endl;
    ch_bisect_MT2_332 = ch_bisect_mt2_332Calculator.mt2_332(vis_A, vis_B, pT_Miss, missing.M()); 
  
  return true;
}

void doubleHiggsAnalyser::Loop() {
  TTree *t = 0;
  bool keep = false;
  if (delphes_flag) t = del_tree;
  else t = fChain;

  for (int iev = 0; iev < t->GetEntries(); iev++) {
    t->GetEntry(iev);
    keep = doubleHiggsAnalyser::Analysis();
    if (keep) out_tree->Fill();
    doubleHiggsAnalyser::ResetVariables();
  }
}

void doubleHiggsAnalyser::Finalize() {
  out_tree->Write();
  out_file->Close();
}

int main(Int_t argc, Char_t** argv)
{
    //TFile *f = TFile::Open("/cms/scratch/jlee/hh/hh_1.root");
    //TTree *tree;
    //f->GetObject("Delphes", tree);
    TChain *tree = new TChain("Delphes");

    tree->Add("/cms/scratch/jlee/hh/*.root");

    tree->SetBranchStatus("*",true);
    
    // for delphes format analysis
    doubleHiggsAnalyser ana(tree, true, true);
    // for nanoAOD format analysis
    // double HiggsAnalyser ana(tree, true);
    ana.Initiate("test.root"); // usage : Initiate("outputfilename.root")
    ana.Loop(); // Loop through the events and do doubleHIggsAnalyser::Analysis() per event.
    ana.Finalize(); // Write the tree and Close the file.

  return 0;
}

