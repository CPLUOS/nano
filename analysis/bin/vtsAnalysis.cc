#include <iostream>
#include <vector>
#include <TLorentzVector.h>

#include "nano/analysis/src/hadAnalysis.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;
/*
TParticle topAnalysis::GetTParticle(int pdgId, int idx) {
  TLorentzVector tlv;
  if      (abs(pdgId) == 11) tlv.SetPtEtaPhiM(Electron_pt[idx], Electron_eta[idx], Electron_phi[idx], Electron_mass[idx]);
  else if (abs(pdgId) == 13) tlv.SetPtEtaPhiM(Muon_pt[idx], Muon_eta[idx], Muon_phi[idx], Muon_mass[idx]);
  else if (abs(pdgId) > 100) tlv.SetPtEtaPhiM(had_pt[idx], had_eta[idx], had_phi[idx], had_mass[idx]);
  else tlv.SetPtEtaPhiM(Jet_pt[idx], Jet_eta[idx], Jet_phi[idx], Jet_mass[idx]);
  return TParticle(pdgId, 0, 0, 0, 0, 0, tlv, tlv);
}
*/
//Double_t hadAnalysis::DeltaR(Double_t deta, Double_t dphi) {
//  return TMath::Sqrt(deta*deta + dphi*dphi);
//}

Double_t hadAnalysis::DeltaPhi(Double_t phi1, Double_t phi2) {
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

Double_t hadAnalysis::GetD(float pt, float eta, float phi, float m, float vx, float vy, float vz) {
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
    hadAnalysis ana(inTree,false,false,false,false);
    ana.setOutput("nanotree.root");
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
      hadAnalysis ana(inTree,false,false,false,false);

      string outFileName = outFileDir+"/nanotree_"+to_string(i-3)+".root";
      ana.setOutput(outFileName);
      ana.Loop();
    }
  }
}

void hadAnalysis::Loop() {
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntries();

  // Events loop
  for (Long64_t iev=0; iev<nentries; iev++) {
    Long64_t entry = LoadTree(iev);
    fChain->GetEntry(entry);
    if (iev%10000 == 0) cout << iev << "/" << nentries << endl;

//    std::cout << iev << " th event" << std::endl;

    ResetBranch();
    EventSelection();
    MatchingForMC();
    HadronAnalysis();
    m_tree->Fill();
  }
//  for (int i =1; i<10; i++) { cout << h_cutFlow->GetBinContent(i) << "  "; }
//  m_output->Write();
//  m_output->Close();
  cout << endl;

}

void hadAnalysis::setOutput(std::string outFileName)
{
  m_output = TFile::Open(outFileName.c_str(), "recreate");

  m_tree = new TTree("event", "event");
  MakeBranch(m_tree);
}

void hadAnalysis::MakeBranch(TTree* t)
{
  m_tree->Branch("channel", &b_channel, "channel/I");
  m_tree->Branch("njet", &b_njet, "njet/I");
  m_tree->Branch("met", &b_met, "met/F");
  m_tree->Branch("dilep_tlv", "TLorentzVector", &b_dilep_tlv);

  m_tree->Branch("had_tlv", "TLorentzVector", &b_had_tlv);
  m_tree->Branch("isFrom_had", &b_isFrom_had, "isFrom_had/I");
  m_tree->Branch("isHadJetMatched_had", &b_isHadJetMatched_had, "isHadJetMatched_had/O");
  m_tree->Branch("d_had", &b_d_had , "d_had/F" );
  m_tree->Branch("x_had", &b_x_had, "x_had/F");
  m_tree->Branch("dr_had", &b_dr_had, "dr_had/F"); // distance between hadron and jet-center
  m_tree->Branch("lxy_had", &b_lxy_had, "lxy_had/F");
  m_tree->Branch("lxySig_had", &b_lxySig_had, "lxySig_had/F");
  m_tree->Branch("angleXY_had", &b_angleXY_had, "angleXY_had/F");
  m_tree->Branch("angleXYZ_had", &b_angleXYZ_had, "angleXYZ_had/F");
  m_tree->Branch("chi2_had", &b_chi2_had, "chi2_had/F");
  m_tree->Branch("dca_had", &b_dca_had, "dca_had/F");

  m_tree->Branch("pt_had", &b_pt_had, "pt_had_had/F");
  m_tree->Branch("eta_had", &b_eta_had, "eta_had/F");
  m_tree->Branch("l3D_had", &b_l3D_had, "l3D_had/F");
  m_tree->Branch("l3DSig_had", &b_l3DSig_had, "l3DSig_had/F");
  m_tree->Branch("legDR_had", &b_legDR_had, "legDR_had/F");
  m_tree->Branch("mass_had", &b_mass_had, "mass_had/F");
  m_tree->Branch("pdgId_had", &b_pdgId_had, "pdgId_had/I");

  m_tree->Branch("dau1_chi2_had", &b_dau1_chi2_had, "dau1_chi2_had/F");
  m_tree->Branch("dau1_ipsigXY_had", &b_dau1_ipsigXY_had, "dau1_ipsigXY_had/F");
  m_tree->Branch("dau1_ipsigZ_had", &b_dau1_ipsigZ_had, "dau1_ipsigZ_had/F");
  m_tree->Branch("dau1_pt_had", &b_dau1_pt_had, "dau1_pt_had/F");

  m_tree->Branch("dau2_chi2_had", &b_dau2_chi2_had, "dau2_chi2_had/F");
  m_tree->Branch("dau2_ipsigXY_had", &b_dau2_ipsigXY_had, "dau2_ipsigXY_had/F");
  m_tree->Branch("dau2_ipsigZ_had", &b_dau2_ipsigZ_had, "dau2_ipsigZ_had/F");
  m_tree->Branch("dau2_pt_had", &b_dau2_pt_had, "dau2_pt_had/F");

  m_tree->Branch("btagCSVV2_Jet", &b_btagCSVV2_Jet, "btagCSVV2_Jet/F");
  m_tree->Branch("btagDeepB_Jet", &b_btagDeepB_Jet, "btagDeepB_Jet/F");
  m_tree->Branch("btagDeepC_Jet", &b_btagDeepC_Jet, "btagDeepC_Jet/F");
  m_tree->Branch("btagCMVA_Jet", &b_btagCMVA_Jet, "btagCMVA_Jet/F");

  m_tree->Branch("area_Jet", &b_area_Jet, "are_Jet/F");
  m_tree->Branch("pt_Jet", &b_pt_Jet, "pt_Jet/F");
  m_tree->Branch("nConstituents_Jet", &b_nConstituents_Jet, "nConstituents_Jet/I");
  m_tree->Branch("nElectrons_Jet", &b_nElectrons_Jet, "nElectrons_Jet/I");
  m_tree->Branch("nMuons_Jet", &b_nMuons_Jet, "nMuons_Jet/I");

}

void hadAnalysis::ResetBranch() {
  m_isMC = false;

  b_channel = -9; b_njet = -9;
  b_met = -99;
  b_dilep_tlv.SetPtEtaPhiM(0,0,0,0);
  recoleps.clear();

  b_had_tlv.SetPtEtaPhiM(0,0,0,0);
  b_isFrom_had = -99;
  b_isHadJetMatched_had = false;
  b_d_had = -1; b_x_had = -1; b_dr_had = -1;
  b_lxy_had = -1; b_lxySig_had = -1; b_angleXY_had = -1; b_angleXYZ_had = -1; b_chi2_had = -1; b_dca_had = -1;
  b_pt_had = -1; b_eta_had = -99; b_l3D_had = -1; b_l3DSig_had = -1; b_legDR_had = -1; b_mass_had = -99; b_pdgId_had = -99;
  b_dau1_chi2_had = -1; b_dau1_ipsigXY_had = -1; b_dau1_ipsigZ_had = -1; b_dau1_pt_had = -1;
  b_dau2_chi2_had = -1; b_dau2_ipsigXY_had = -1; b_dau2_ipsigZ_had = -1; b_dau2_pt_had = -1;

  b_btagCSVV2_Jet = -99; b_btagCMVA_Jet = -99; b_btagDeepB_Jet = -99; b_btagDeepC_Jet = -99;
  b_area_Jet = -1; b_pt_Jet = -1; b_nConstituents_Jet = -1; b_nElectrons_Jet = -1; b_nMuons_Jet = -1;

  qjMapForMC_.clear(); qMC_.clear(); genJet_.clear(); recoJet_.clear();
}


void hadAnalysis::EventSelection() {
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
/*
vector<TParticle> topAnalysis::muonSelection() {
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

vector<TParticle> topAnalysis::elecSelection() {
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

vector<TParticle> topAnalysis::jetSelection() {
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
*/
void hadAnalysis::MatchingForMC() {
  m_isMC = true;

  //Find s quark from Gen Info.  
  for (unsigned int i=0; i<nGenPart; ++i) {
    if (std::abs(GenPart_status[i] - 25) < 5 && abs(GenPart_pdgId[i]) == 3) {
      qMC_.push_back(i);
    }
  }
  if (qMC_.size() == 0) {
    ++b_chk;
    return;
  }

  auto q1 = qMC_[0];
  TLorentzVector q1_tlv; 
  q1_tlv.SetPtEtaPhiM(GenPart_pt[q1], GenPart_eta[q1], GenPart_phi[q1], GenPart_mass[q1]);

  unsigned int q2;
  TLorentzVector q2_tlv; 
  if (qMC_.size() == 2) {
    q2 = qMC_[1];
    q2_tlv.SetPtEtaPhiM(GenPart_pt[q2], GenPart_eta[q2], GenPart_phi[q2], GenPart_mass[q2]);
  }

  //Gen Particle & Reco Jet matching
  for (unsigned int j=0; j<nJet;++j) {
    struct JetStat Stat;
    TLorentzVector jet_tlv;
    jet_tlv.SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);

    if (q1_tlv.DeltaR(jet_tlv) < 0.3) {
      qjMapForMC_.insert({j, GenPart_pdgId[q1]});
      Stat.idx = j;
      Stat.dr = q1_tlv.DeltaR(jet_tlv);
      Stat.matchedQuark = qjMapForMC_[j];
      recoJet_.push_back(Stat);
    }
    else {
      if (qMC_.size() == 1) qjMapForMC_.insert({j, -9});
      else if (qMC_.size() == 2) {
        if (q2_tlv.DeltaR(jet_tlv) < 0.3) {
          qjMapForMC_.insert({j, GenPart_pdgId[q2]});
          Stat.idx = j;
          Stat.dr = q2_tlv.DeltaR(jet_tlv);
          Stat.matchedQuark = qjMapForMC_[j];
          recoJet_.push_back(Stat);
        }
        else qjMapForMC_.insert({j, -9});
      }
    }
  }

  if (recoJet_.size() == 0) return;
  else if (recoJet_.size() > 1) {
    if ( recoJet_[0].matchedQuark == recoJet_[1].matchedQuark ) {
      sort(recoJet_.begin(), recoJet_.end(), [](struct JetStat a, struct JetStat b) {return a.dr < b.dr;}); // pick the closest matched recojet
      qjMapForMC_[recoJet_[1].idx] = -9;
    }
  }
}

void hadAnalysis::HadronAnalysis() {
  std::vector<std::vector<struct HadStat>> JetCollection;
  std::vector<struct HadStat> Had;
  struct HadStat Stat;
  //Reco Jet & Reco KS matching 
  for (unsigned int j=0; j<nJet; ++j){
    Had.clear();
    TLorentzVector jet_tlv;
    jet_tlv.SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]); 
    for (unsigned int k=0; k<nhad; ++k) {
      if (abs(had_pdgId[k]) == 310) {
        if ( (had_lxy[k]/had_lxyErr[k]) < 3 ) continue;
        if ( had_angleXY[k] < 0.98 ) continue;
        if ( had_chi2[k] > 3 ) continue;
        if ( had_dca[k] > 1 ) continue;
      }
      else if (abs(had_pdgId[k]) == 3122) {
        if ( (had_lxy[k]/had_lxyErr[k]) < 3 ) continue;
        if ( had_angleXY[k] < 0.98 ) continue;
        if ( had_chi2[k] > 3 ) continue;
        if ( had_dca[k] > 1 ) continue;
      }
      else continue;
      TLorentzVector had_tlv;
      had_tlv.SetPtEtaPhiM(had_pt[k], had_eta[k], had_phi[k], had_mass[k]);
      if (jet_tlv.DeltaR(had_tlv) < 0.3 && (had_pt[k]/Jet_pt[j]) > 0.15) {
        Stat.idx = k;
        Stat.pdgId = had_pdgId[k];
        Stat.x = had_pt[k]/Jet_pt[j];
        if (m_isMC) Stat.label = qjMapForMC_[j];
        Stat.isHadJetMatched = true;
        Stat.dr = jet_tlv.DeltaR(had_tlv);
        Stat.jetIdx = j;
        Had.push_back(Stat);
      }
    }
    if (Had.size() != 0 ) {
      if (Had.size() > 1) sort(Had.begin(), Had.end(), [](struct HadStat a, struct HadStat b) {return a.x > b.x;}); // pick hadron with highest x in the jet
      JetCollection.push_back(Had);
    }
  }
  if (JetCollection.size() != 0 ) {
    if (JetCollection.size() > 1) sort(JetCollection.begin(), JetCollection.end(), [](std::vector<struct HadStat> a, std::vector<struct HadStat> b) {return a[0].x > b[0].x;}); // pick jet-hadron pair with highest x
    auto idx = JetCollection[0][0].idx;
    auto jidx = JetCollection[0][0].jetIdx;
    b_had_tlv.SetPtEtaPhiM(had_pt[idx], had_eta[idx], had_phi[idx], had_mass[idx]);

    b_x_had = JetCollection[0][0].x;
    b_isFrom_had = JetCollection[0][0].label;  // -99 : there is no matching between had and jet, -9 : there is t->s in the event,but not matched to jet, 0 : there is no t->s in the event (if no t->s and no matching between had-jet, then the event would be -99), +-3 : hadron is from t->s
    b_isHadJetMatched_had = JetCollection[0][0].isHadJetMatched;
    b_dr_had = JetCollection[0][0].dr;
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

    b_btagCSVV2_Jet = Jet_btagCSVV2[jidx];
    b_btagCMVA_Jet = Jet_btagCMVA[jidx]; 
    b_btagDeepB_Jet = Jet_btagDeepB[jidx];
    b_btagDeepC_Jet = Jet_btagDeepC[jidx];
    b_area_Jet = Jet_area[jidx]; 
    b_pt_Jet = Jet_pt[jidx];
    b_nConstituents_Jet = Jet_nConstituents[jidx]; 
    b_nElectrons_Jet = Jet_nElectrons[jidx];
    b_nMuons_Jet = Jet_nMuons[jidx];
  }
/*
  std::cout <<" b_chk : " << b_chk << std::endl;
  std::cout <<"JetCol : " << JetCollection.size() << std::endl;
  std::cout <<"isFrom : " << b_isFrom_had << std::endl;
  std::cout <<"isMatc : " << b_isHadJetMatched_had << std::endl;
*/
}
