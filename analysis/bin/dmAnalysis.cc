#define dmAnalysis_cxx
#include "nano/analysis/interface/dmAnalysis.h"
#include <TH2.h>
#include <TCanvas.h>
#include <iostream>
#include <cstdlib>


using namespace std;

vector<TParticle> dmAnalysis::muonSelection()
{
  vector<TParticle> muons;
  for (UInt_t i=0; i < nMuon; ++i){
    if (!Muon_tightId[i]) continue;
    if (Muon_pt[i] < 30) continue;
    if (std::abs(Muon_eta[i] > 2.1)) continue;
    if (Muon_pfRelIso04_all[i] > 0.15) continue;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]);
    auto muon = TParticle();
    muon.SetPdgCode(13*Muon_charge[i]*-1);
    muon.SetMomentum(mom);

    muons.push_back(muon);
  }
  return muons;
}


vector<TParticle> dmAnalysis::elecvetoSelection()
{
  vector<TParticle> elecs;
  for (UInt_t i=0; i < nElectron; ++i){
    if (Electron_cutBased[i] > 10) break;
    if (Electron_pt[i] > 10) break;
    if (std::abs(Electron_eta[i]) < 2.1) break;
    if (Electron_pfRelIso03_all[i] < 0.25) break;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Electron_pt[i], Electron_eta[i], Electron_phi[i], Electron_mass[i]);
    recoleps.push_back(mom);
  
    auto elec = TParticle();
    elec.SetPdgCode(11*Electron_charge[i]*-1);
    elec.SetMomentum(mom);

    elecs.push_back(elec);
  }
  return elecs;
}


vector<TParticle> dmAnalysis::tauvetoSelection()
{
  vector<TParticle> taus;
  for (UInt_t i=0; i < nTau; ++i){
    if (Tau_pt[i] > 20) break;
    if (std::abs(Tau_eta[i]) < 2.1) break;
    if (Tau_decayMode[i] == 1) break;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Tau_pt[i], Tau_eta[i], Tau_phi[i], Tau_mass[i]);
    bool hasOverLap = false;
    for (auto lep : recoleps){
        if (mom.TLorentzVector::DeltaR(lep) < 0.3) hasOverLap = true;
    }
    if (!hasOverLap) break;
    recoleps.push_back(mom); 
    auto tau = TParticle();
    tau.SetPdgCode(15*Tau_charge[i]*-1);
    tau.SetMomentum(mom);

    taus.push_back(tau);
    //recoleps.push_back(tau);
  }
  return taus;
}


vector<TParticle> dmAnalysis::jetSelection()
{
  vector<TParticle> jets;
  for (UInt_t i=0; i < nJet; ++i){
    if (Jet_jetId[i] < 1) continue;
    if (Jet_pt[i] < 60) continue;
    if (std::abs(Jet_eta[i]) > 5) continue;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
    auto jet = TParticle();
    jet.SetMomentum(mom);
    jets.push_back(jet);
  }
  return jets;
}


vector<TParticle> dmAnalysis::bjetvetoSelection()
{
  vector<TParticle> bjets;
  for (UInt_t i=0; i < nJet; ++i){
    if (Jet_btagCSVV2[i] > 0.8484) break;
    if (Jet_pt[i] > 30) break;
    if (std::abs(Jet_eta[i]) < 2.4) break;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
    auto bjet = TParticle();
    bjet.SetMomentum(mom);
    bjets.push_back(bjet);
  }
  return bjets;
}

//vector<TParticle> dmAnalysis::MET()
//{
//  //vector<TParticle> met;
//  TLorentzVector mom;
//  mom.SetPtEtaPhiM(MET_pt, 0, MET_phi, 0);
//  auto met = TParticle();
//  met.SetMomentum(mom);
//  return met;
//}


bool dmAnalysis::analysis()
{
  //bool m_isMC;
  //cout << m_isMC << endl;
  if (m_isMC) {
    cout << genWeight << endl;
    Int_t nvtx = Pileup_nTrueInt;
    b_puweight = m_pileUp->getWeight(nvtx);

    b_genweight = genWeight;
    h_genweights->Fill(0.5, b_genweight);
    b_weight = b_genweight * b_puweight;
    h_weights->Fill(0.5, b_weight);
    
  }
  else
  {
    b_puweight = 1;
    b_genweight = 0;
    //if (!(m_lumi->LumiCheck(run, luminosityBlock))) return false;
    //cout << "nono: " << b_puweight << ", " <<  b_genweight << endl;
  }
  //cout << "np" << endl;
  h_nevents->Fill(0.5,b_genweight*b_puweight);

  // if (std::abs(PV_z) >= 24.) return false;
  // if (PV_npvs == 0) return false;
  // if (PV_ndof < 4) return false;

  auto muons = muonSelection();
  //auto met = MET(); 
  //b_mpt = met.PT(); 
  
  TLorentzVector mom;
  mom.SetPtEtaPhiM(MET_pt, 0, MET_phi, 0);
  auto met = TParticle();
  met.SetMomentum(mom);
  
  b_mpt = met.Pt(); 

  if (muons.size() == 0) return false;

  //int mulpdg = -1;
  if (muons.size() > 1) {
    recolep1 = muons[0];
    recolep2 = muons[1];
    //mulpdg = muons[0].GetPdgCode()*muons[1].GetPdgCode();
    recolep1.Momentum(b_mu1);
    recolep2.Momentum(b_mu2);
  
    recoleps.push_back(b_mu1);
    recoleps.push_back(b_mu2);
   
    b_dimu = b_mu1 + b_mu2;
    b_channel = CH_ZJets;
  }
  
  if (muons.size() == 1) {
    recolep1 = muons[0];
    recomet = met;
    recolep1.Momentum(b_mu1);
    recomet.Momentum(b_met1);

    recoleps.push_back(b_mu1);
    b_mume = b_mu1 + b_met1;
    b_channel = CH_WJets;
  }


  //cout << b_channel << endl;
  // Triggers
  b_trigger_dm = HLT_IsoMu24;
  if (!(b_trigger_dm)) return false;

  if (b_channel == CH_ZJets) {
    if (abs(b_dimu.M()-90) > 30) return false; 
  }
  if (b_channel == CH_WJets) {
    if (abs(b_mume.Mt()-75) > 25) return false;
  }
  b_step1 = true;
  b_step = 1;
 
  b_nelec = nElectron; 

  auto elecs = elecvetoSelection();
  if (elecs.size() != 0) return false;
  //for (UInt_t i=0; i < nElectron; ++i){
  //  recolep = elecs[i];
  //  recolep.Momentum(b_elec);
  //  recoleps.push_back(b_elec);
  //}
  b_step2 = true;
  if (b_step == 1){
      ++b_step;
  }
 
  b_ntau = nTau; 

  auto taus = tauvetoSelection();
  if (taus.size() != 0) return false;
  b_step3 = true;
  if (b_step == 2){
    ++b_step;
  }
  
  auto jets = jetSelection();
  b_njet = jets.size(); 
  if (jets.size() < 2) return false;
  b_step4 = true;
  if (b_step == 3){
    ++b_step;
  }

  auto bjets = bjetvetoSelection();
  b_nbjet = bjets.size(); 
  if (b_nbjet != 0) return false;
  b_step5 = true;
  if (b_step == 4){
    ++b_step;
  }

  // CR2
  //if (b_mpt < 250) return false;
  //b_step6 = true;
  //if (b_step == 5){
  //  ++b_step;
  //}

  //// CR3
  //recojet1 = jets[0]; 
  //recojet2 = jets[1];

  //recojet1.Momentum(b_jet1);
  //recojet2.Momentum(b_jet2);


  //recojets.push_back(b_jet1);
  //recojets.push_back(b_jet2);

  //b_dijet = b_jet1 + b_jet2;
 
  //b_eta2 = recojet1.Eta()*recojet2.Eta();
  //b_dEta = recojet1.Eta()-recojet2.Eta();

  //if (b_eta2 >= 0) return false;
  //if (abs(b_dEta) < 3.8) return false;
  //if (b_dijet.M() < 1000) return false; 
  //b_step7 = true;
  //if (b_step == 6){
  //  ++b_step;
  //}

  return true;
}


void dmAnalysis::Loop()
{
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntries();
  cout << "# of the total events: " << nentries << endl;

  for (Long64_t iev=0; iev<nentries; iev++){
    resetBranch();
    fChain->GetEntry(iev);
    bool keep = analysis();
    if (keep ==1 ) cout << keep << endl;
    if (keep){
    m_tree->Fill();
    }
  }

  //m_output->Write();
}


int main(Int_t argc, Char_t** argv)
//int main(int argc, char* argv[])
{
  string env = getenv("CMSSW_BASE");
  string username = getenv("USER");
  
  if(argc != 1){
    string dirName = "root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/"+username+"/nanoAOD/"+std::string(argv[1])+"/"+std::string(argv[2]);
    string temp = argv[2];
    Bool_t isMC = false;
    Size_t found = temp.find("Run");
    if(found == string::npos) isMC = true;
    for(Int_t i = 3; i < argc; i++){
      TFile *f = TFile::Open(argv[i], "read");
      TTree *tree;
      f->GetObject("Events", tree);

      temp = argv[i];
      found = temp.find_last_of('/');
      std::string outPutName = dirName+temp.substr(found);
      dmAnalysis t(tree, isMC);

      //t.LoadModules(pileUp, lumi);
      t.setOutput(outPutName);
      t.Loop();
    }
  }

  else
  {
   TFile *f = TFile::Open("/xrootd/store/group/nanoAOD/run2_2016v5/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180610_152543/0000/nanoAOD_16.root", "read");
   //TFile *f = TFile::Open("/xrootd/store/group/nanoAOD/run2_2016v4/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180430_152541/0000/nanoAOD_256.root", "read");
    //TFile *f = TFile::Open("/xrootd/store/group/nanoAOD/run2_2016v5/VBF-C1N2_WZ_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUSummer16Fast_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180610_151009/0000/nanoAOD_167.root", "read");
    TTree *tree;
    Bool_t isMC = false;
    string temp = "Run";
    Size_t found = temp.find("Run");
    if(found == string::npos) isMC = true;
    f->GetObject("Events", tree);

    dmAnalysis t(tree, isMC);
    //t.LoadModules(pileUp, lumi);
    t.setOutput("test.root");
    t.Loop();
  }
  return 0;
}


void dmAnalysis::setOutput(std::string outputName)
{
  m_output = TFile::Open(outputName.c_str(), "recreate");
  m_tree = new TTree("event", "event");
  MakeBranch(m_tree);

  h_nevents = new TH1D("nevents", "nevents", 1, 0, 1);
  h_genweights = new TH1D("genweight", "genweight", 1, 0, 1);
  h_weights = new TH1D("weight", "weight", 1, 0, 1);

}


void dmAnalysis::MakeBranch(TTree* t)
{

  t->Branch("nvertex", &b_nvertex, "nvertex/I");
  t->Branch("step", &b_step, "step/I");
  t->Branch("channel", &b_channel, "channel/I");
  t->Branch("nelec", &b_nelec, "nelec/I");
  t->Branch("ntau", &b_ntau, "ntau/I");
  t->Branch("njet", &b_njet, "njet/I");
  t->Branch("nbjet", &b_nbjet, "nbjet/I");
  t->Branch("step1", &b_step1, "step1/O");
  t->Branch("step2", &b_step2, "step2/O");
  t->Branch("step3", &b_step3, "step3/O");
  t->Branch("step4", &b_step4, "step4/O");
  t->Branch("step5", &b_step5, "step5/O");
  t->Branch("step6", &b_step6, "step6/O");
  t->Branch("step7", &b_step7, "step7/O");
  t->Branch("tri", &b_trigger_dm, "tri/O");
  t->Branch("mpt", &b_mpt, "mpt/F");
  t->Branch("eta2", &b_eta2, "eta2/F");
  t->Branch("dEta", &b_dEta, "dEta/F");

  t->Branch("mu1", "TLorentzVector", &b_mu1);
  t->Branch("mu2", "TLorentzVector", &b_mu2);
  t->Branch("met1", "TLorentzVector", &b_met1);
  t->Branch("elec", "TLorentzVector", &b_elec);
  t->Branch("dimu", "TLorentzVector", &b_dimu);
  t->Branch("mume", "TLorentzVector", &b_mume);
  
  t->Branch("jet1", "TLorentzVector", &b_jet1);
  t->Branch("jet2", "TLorentzVector", &b_jet2);
  t->Branch("dijet", "TLorentzVector", &b_dijet);

  t->Branch("weight", &b_weight, "weight/F");
  t->Branch("puweight", &b_puweight, "puweight/F");
  t->Branch("genweight", &b_genweight, "genweight/F");
  t->Branch("csvweight", &b_csvweights, "csvweight/F");
  //dm->Branch("PV_npvs", &PV_npvs, "PV_npvs/F");
}


void dmAnalysis::resetBranch()
{
  b_mu1.SetPtEtaPhiM(0,0,0,0);
  b_mu2.SetPtEtaPhiM(0,0,0,0);
  b_dimu.SetPtEtaPhiM(0,0,0,0);
  b_met1.SetPtEtaPhiM(0,0,0,0);
  b_elec.SetPtEtaPhiM(0,0,0,0);
  b_mume.SetPtEtaPhiM(0,0,0,0);
  b_jet1.SetPtEtaPhiM(0,0,0,0);
  b_jet2.SetPtEtaPhiM(0,0,0,0);
  b_dijet.SetPtEtaPhiM(0,0,0,0);

  recoleps.clear();
  recojets.clear();
  b_csvweights.clear();
  b_trigger_dm = 0;
  //b_mpt = -9; 
  b_weight = 1; b_genweight = 1; b_puweight = 1;
  b_mpt = -9; b_eta2 = -9; b_dEta = -1;

  b_nvertex = 0; b_step = -1; b_channel = 0;
  b_nelec = 0; b_ntau = 0; b_njet = 0; b_nbjet = 0;
  b_step1 = 0; b_step2 = 0; b_step3 = 0; b_step4 = 0; b_step5 = 0; b_step6 = 0; b_step7 = 0;

}
