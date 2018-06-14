#include <iostream>
#include <string>
#include <vector>
#include <TLorentzVector.h>

#include "nano/analysis/interface/vtsAnalyser.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;

string getFileName(const string& s) {
 char sep = '/';
 size_t i = s.rfind(sep, s.length());
 if (i != string::npos) {
    return(s.substr(i+1, s.length() - i));
 }
 return("");
}
string getDir (const string& path){
  size_t found = path.find_last_of("/\\");
  return(path.substr(0, found));
}
string getType (const string& path){
  size_t found = path.find("tt");
  return(path.substr(found));
}


int main(Int_t argc, Char_t** argv) {

  string hostDir = "root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/";

  if (argc <= 1) {
    cout << "no input file is specified. running with default file." << endl;
    auto inFile = TFile::Open("/xrootd/store/group/nanoAOD/run2_2016v4/tsw/nanoAOD_1.root", "READ");
    auto inTree = (TTree*) inFile->Get("Events");
    vtsAnalyser ana(inTree,inTree,0,true,false,false,false);
    ana.setOutput("nanotree.root");
    ana.Loop();
  } else if (argc == 5) {
    auto trn1 = argv[1];
    auto trn2 = argv[2];
    auto trn3 = argv[3];
    auto outn = argv[4]; 
    
    TFile *f1 = TFile::Open(trn1, "READ");
    TFile *f2 = TFile::Open(trn2, "READ");
    TFile *f3 = TFile::Open(trn3, "READ");
    TTree *t1 = (TTree*)f1->Get("Events");
    TTree *t2 = (TTree*)f2->Get("Events");
    TTree *t3 = (TTree*)f3->Get("Events");
   
    cout << t1 << " , " << t2 << " , " << t3 << endl; 
    vtsAnalyser ana(t1,t2,t3,true,false,false,false);
    ana.setOutput(outn);
    ana.Loop();
  } else {
    string jobName    = string(argv[1]);
    string sampleName = string(argv[2]);

    // temp
    Bool_t isMC = false;
    std::string temp = argv[1];
    Size_t found = temp.find("run");
    Size_t found2 = temp.find("NANOAOD");
    if (found == std::string::npos || found2 == std::string::npos) isMC = true;

    string outFileDir = hostDir+getenv("USER")+"/"+jobName+"/"+sampleName;
    for (Int_t i = 3; i < argc; i++) {
      auto inFileName = argv[i];
    // temp
      auto fileName = getFileName(argv[i]);
      auto dirName = getDir(getDir(argv[i]));
      auto sampleType = getType(dirName);
      //NANO
      TFile *inFile = TFile::Open(inFileName, "READ");
      TTree *inTree = (TTree*) inFile->Get("Events");
      TFile *hadFile;
      TTree *hadTree=0;
      TFile *hadTruthFile;
      TTree *hadTruthTree=0;
      TString hadFileName = dirName + "/HADAOD/" + fileName;
      TString hadTruthFileName = dirName + "/HADTRUTHAOD/" + fileName;
      if (found2 == std::string::npos) {
        hadFile = TFile::Open(hadFileName, "READ");
        hadTree = (TTree*) hadFile->Get("Events");
        hadTruthFile = TFile::Open(hadTruthFileName, "READ");
        hadTruthTree = (TTree*) hadTruthFile->Get("Events");
      }
      cout << "dirName : " << dirName << " fileName : " << fileName << endl;
      cout << "tree chk : had : " << hadTree << " , hadTruth : " << hadTruthTree << endl;
      vtsAnalyser ana(inTree,hadTree,hadTruthTree,isMC,false,false,false);
      string outFileName = outFileDir+"/nanotree_"+to_string(i-3)+".root";
      ana.setOutput(outFileName);
      ana.Loop();
    }
  }
}

void vtsAnalyser::Loop() {
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntries();

  // Events loop
  for (Long64_t iev=0; iev<nentries; iev++) {
    fChain->GetEntry(iev);
    if (h_fChain) h_fChain->GetEntry(iev);
    if (ht_fChain) ht_fChain->GetEntry(iev);
    cout << "event : " << iev << endl;
    cout << "nhadTruth : " << nhadTruth << endl;
    cout << "nhad      : " << nhad << endl;
    if (iev%10000 == 0) cout << iev << "/" << nentries << endl;
    ResetBranch();
    int PassedStep = EventSelection();
    if (PassedStep >= 0) {
      MatchingForMC();
      HadronAnalysis();
    }
    m_tree->Fill();
  }
}

void vtsAnalyser::setOutput(std::string outFileName) {
  m_output = TFile::Open(outFileName.c_str(), "recreate");
  m_tree = new TTree("event", "event");

  MakeBranch(m_tree);

  h_nevents = new TH1D("nevents", "nevents", 1, 0, 1);
  h_genweights = new TH1D("genweight", "genweight", 1, 0, 1);
  h_weights = new TH1D("weight", "weight", 1, 0, 1);
  h_cutFlow = new TH1D("cutflow", "cutflow", 11, -0.5, 10.5);
}

void vtsAnalyser::MakeBranch(TTree* t) {
  t->Branch("channel", &b_channel, "channel/I");
  t->Branch("njet", &b_njet, "njet/I");
  t->Branch("met", &b_met, "met/F");
  t->Branch("dilep_tlv", "TLorentzVector", &b_dilep);

  t->Branch("had_tlv", "TLorentzVector", &b_had_tlv);
  t->Branch("isFrom_had", &b_isFrom_had, "isFrom_had/I");
  t->Branch("isHadJetMatched_had", &b_isHadJetMatched_had, "isHadJetMatched_had/O");
  t->Branch("d_had", &b_d_had , "d_had/F" );
  t->Branch("x_had", &b_x_had, "x_had/F");
  t->Branch("dr_had", &b_dr_had, "dr_had/F"); // distance between hadron and jet-center
  t->Branch("lxy_had", &b_lxy_had, "lxy_had/F");
  t->Branch("lxySig_had", &b_lxySig_had, "lxySig_had/F");
  t->Branch("angleXY_had", &b_angleXY_had, "angleXY_had/F");
  t->Branch("angleXYZ_had", &b_angleXYZ_had, "angleXYZ_had/F");
  t->Branch("chi2_had", &b_chi2_had, "chi2_had/F");
  t->Branch("dca_had", &b_dca_had, "dca_had/F");

  t->Branch("pt_had", &b_pt_had, "pt_had_had/F");
  t->Branch("eta_had", &b_eta_had, "eta_had/F");
  t->Branch("l3D_had", &b_l3D_had, "l3D_had/F");
  t->Branch("l3DSig_had", &b_l3DSig_had, "l3DSig_had/F");
  t->Branch("legDR_had", &b_legDR_had, "legDR_had/F");
  t->Branch("mass_had", &b_mass_had, "mass_had/F");
  t->Branch("pdgId_had", &b_pdgId_had, "pdgId_had/I");

  t->Branch("dau1_chi2_had", &b_dau1_chi2_had, "dau1_chi2_had/F");
  t->Branch("dau1_ipsigXY_had", &b_dau1_ipsigXY_had, "dau1_ipsigXY_had/F");
  t->Branch("dau1_ipsigZ_had", &b_dau1_ipsigZ_had, "dau1_ipsigZ_had/F");
  t->Branch("dau1_pt_had", &b_dau1_pt_had, "dau1_pt_had/F");

  t->Branch("dau2_chi2_had", &b_dau2_chi2_had, "dau2_chi2_had/F");
  t->Branch("dau2_ipsigXY_had", &b_dau2_ipsigXY_had, "dau2_ipsigXY_had/F");
  t->Branch("dau2_ipsigZ_had", &b_dau2_ipsigZ_had, "dau2_ipsigZ_had/F");
  t->Branch("dau2_pt_had", &b_dau2_pt_had, "dau2_pt_had/F");

  t->Branch("btagCSVV2_Jet", &b_btagCSVV2_Jet, "btagCSVV2_Jet/F");
  t->Branch("btagDeepB_Jet", &b_btagDeepB_Jet, "btagDeepB_Jet/F");
  t->Branch("btagDeepC_Jet", &b_btagDeepC_Jet, "btagDeepC_Jet/F");
  t->Branch("btagCMVA_Jet", &b_btagCMVA_Jet, "btagCMVA_Jet/F");

  t->Branch("area_Jet", &b_area_Jet, "are_Jet/F");
  t->Branch("pt_Jet", &b_pt_Jet, "pt_Jet/F");
  t->Branch("nConstituents_Jet", &b_nConstituents_Jet, "nConstituents_Jet/I");
  t->Branch("nElectrons_Jet", &b_nElectrons_Jet, "nElectrons_Jet/I");
  t->Branch("nMuons_Jet", &b_nMuons_Jet, "nMuons_Jet/I");

}

void vtsAnalyser::ResetBranch() {
  Reset(); 
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

void vtsAnalyser::MatchingForMC() {
  //Find s quark from Gen Info.  
  for (unsigned int i=0; i<nGenPart; ++i) {
    if (std::abs(GenPart_status[i] - 25) >= 5 ) continue;
    if (abs(GenPart_pdgId[i]) == 3 || abs(GenPart_pdgId[i]) == 5) qMC_.push_back(i);
  }

  if (qMC_.size() < 2) { cout << " >>>>> it's not a tsWbW event. save nothing." << endl; return; }

  auto q1 = qMC_[0];
  TLorentzVector q1_tlv; 
  q1_tlv.SetPtEtaPhiM(GenPart_pt[q1], GenPart_eta[q1], GenPart_phi[q1], GenPart_mass[q1]);
  auto q2 = qMC_[1];
  TLorentzVector q2_tlv; 
  q2_tlv.SetPtEtaPhiM(GenPart_pt[q2], GenPart_eta[q2], GenPart_phi[q2], GenPart_mass[q2]);

  //Gen Particle & Reco Jet matching
  for (unsigned int j=0; j<nJet;++j) {
    TLorentzVector jet_tlv;
    jet_tlv.SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
    auto qIdx = q1;
    auto dr = q1_tlv.DeltaR(jet_tlv); 
    auto dr_ = q2_tlv.DeltaR(jet_tlv); 
    if ( dr_ < dr) {
      dr = dr_;
      qIdx = q2;
    }

    if (dr > 0.3) qjMapForMC_.insert({j, -9});
    else {
      qjMapForMC_.insert({j, GenPart_pdgId[qIdx]});
      JetStat stat(j, dr, qjMapForMC_[j]);
      recoJet_.push_back(stat);
    }
  }

  if (recoJet_.size() == 0) return;

  else if (recoJet_.size() > 1) {
    sort(recoJet_.begin(), recoJet_.end(), [](struct JetStat a, struct JetStat b) { return (a.dr < b.dr); } );

    bool idChanged = false;
    for (unsigned int i=1;i<recoJet_.size(); ++i) {
        if (idChanged) qjMapForMC_[recoJet_[i].idx] = -9;
        else { 
            if (recoJet_[0].matchedQuark != recoJet_[i].matchedQuark) idChanged = true;
            else qjMapForMC_[recoJet_[i].idx] = -9;
        }
    }
  }
}

void vtsAnalyser::HadronAnalysis() {
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
    b_isFrom_had = JetCollection[0][0].label;  // -99 : event that can't pass till step4(jet selection) or there is no matching between had and jet, -9 : there is t->qW in the event,but not matched to recoJet, 0 : there is no t->qW in the event (if no t->s and no matching between had-jet, then the event would be -99), +-3 : hadron is from t->sW, +-5 : hadron is from t->bW
    b_isHadJetMatched_had = JetCollection[0][0].isHadJetMatched;
    b_dr_had = JetCollection[0][0].dr;

    b_d_had = GetD(had_pt[idx], had_eta[idx], had_phi[idx], had_mass[idx], had_x[idx], had_y[idx], had_z[idx]);
    b_lxy_had = had_lxy[idx];
    b_lxySig_had = had_lxy[idx]/had_lxyErr[idx];
    b_angleXY_had = had_angleXY[idx];
    b_angleXYZ_had = had_angleXYZ[idx];
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
}

