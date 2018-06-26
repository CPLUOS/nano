#include <iostream>
#include <string>
#include <vector>
#include <TLorentzVector.h>
#include <sys/stat.h>
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
long getFileSize(const char* filename) {
  struct stat st;
  if ( stat(filename, &st) == 0) {
    return st.st_size;
  } else return 0;
}

int main(int argc, char* argv[]) {
  string hostDir = "root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/";

  if (argc <= 1) {
    cout << "no input file is specified. running with default file." << endl;
    auto inFile = TFile::Open("/xrootd/store/group/nanoAOD/run2_2016v4/tsw/nanoAOD_1.root", "READ");
    auto inTree = (TTree*) inFile->Get("Events");
    vtsAnalyser ana(inTree,inTree,0,true,false,false,false);
    ana.setOutput("nanotree.root");
    ana.Loop();
  } else {
    string jobName    = string(argv[1]);
    string sampleName = string(argv[2]);

    Bool_t isMC = (sampleName.find("Run") == std::string::npos);
    string outFileDir = hostDir+getenv("USER")+"/"+jobName+"/"+sampleName;
    for (Int_t i = 3; i < argc; i++) {
      auto inFileName = argv[i];
      if (string(inFileName).find("NANOAOD") == std::string::npos) {cout << " input file has to be nanoAOD " << endl; continue;}

      auto fileName = getFileName(argv[i]);
      auto dirName = getDir(getDir(argv[i]));
//      auto sampleType = getType(dirName);

      TFile *inFile = TFile::Open(inFileName, "READ");
      TTree *inTree = (TTree*) inFile->Get("Events");
      TString hadFileName = dirName + "/HADAOD/" + fileName;
      TString hadTruthFileName = dirName + "/HADTRUTHAOD/" + fileName;

//      TFile *hadFile = TFile::Open(hadFileName, "READ");
//      TTree *hadTree = (TTree*) hadFile->Get("Events");
      TFile *hadTruthFile = TFile::Open(hadTruthFileName, "READ");
      TTree *hadTruthTree = (TTree*) hadTruthFile->Get("Events");

      cout << "dirName : " << dirName << " fileName : " << fileName << endl;
//      cout << "tree chk : had : " << hadTree << " , hadTruth : " << hadTruthTree << endl;
      vtsAnalyser ana(inTree,hadTruthTree,hadTruthTree,isMC,false,false,false); // you don't need to use hadTree
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
//    cout << "nhadTruth : " << nhadTruth << endl;
//    cout << "nhad      : " << nhad << endl;

    if (iev%10000 == 0) cout << iev << "/" << nentries << endl;
    ResetBranch();
    int passedStep = EventSelection();
    if (passedStep >= -1) { // -1 : No event selection, 0 : PV cut, reco lepton cut and so on, 1~4 : step 1 ~ 4
      MatchingForMC();
      HadronAnalysis();
      Test();
    }
    m_tree->Fill();
  }
  if (Test()) {
    cout << "total KS : " << nTotHadKS << " matched (dau cut + real)  : " << nMatchedKS << " total dau_pt cut had KS : " << nHadKSDauCut << " total real KS : " << nRealKSFromTop << endl;
    cout << "total real KS with Jet : " << nRealKSWithJet << " total real KS with Jet + dR and x cut : " << nRealKSWithJetAndCut << "  total KS with cut : " << nKSWithCut << " total KS with cut + dau pt cut : " << nKSWithDauCut << endl;
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
  t->Branch("step", &b_step, "step/I");

  t->Branch("hadTruth_nMatched", &b_hadTruth_nMatched, "hadTruth_nMatched/I");
  t->Branch("hadTruth_nTrueDau", &b_hadTruth_nTrueDau, "hadTruth_nTrueDau/I");
  t->Branch("hadTruth_isHadFromTop", &b_hadTruth_isHadFromTop, "hadTruth_isHadFromTop/O");
  t->Branch("hadTruth_isHadFromTsb", &b_hadTruth_isHadFromTsb, "hadTruth_isHadFromTsb/I");
  t->Branch("hadTruth_isHadFromW", &b_hadTruth_isHadFromW, "hadTruth_isHadFromW/O");
  t->Branch("hadTruth_isHadFromS", &b_hadTruth_isHadFromS, "hadTruth_isHadFromS/O");
  t->Branch("hadTruth_isHadFromC", &b_hadTruth_isHadFromC, "hadTruth_isHadFromC/O");
  t->Branch("hadTruth_isHadFromB", &b_hadTruth_isHadFromB, "hadTruth_isHadFromB/O");

  // For Test
  t->Branch("hadTruth_pt_vec", "vector<float>", &b_hadTruth_pt_vec);
  t->Branch("hadTruth_eta_vec", "vector<float>", &b_hadTruth_eta_vec);
  t->Branch("hadTruth_phi_vec", "vector<float>", &b_hadTruth_phi_vec);
  t->Branch("hadTruth_mass_vec", "vector<float>", &b_hadTruth_mass_vec);
  t->Branch("hadTruth_lxy_vec", "vector<float>", &b_hadTruth_lxy_vec);
  t->Branch("hadTruth_lxySig_vec", "vector<float>", &b_hadTruth_lxySig_vec);
  t->Branch("hadTruth_angleXY_vec", "vector<float>", &b_hadTruth_angleXY_vec);
  t->Branch("hadTruth_angleXYZ_vec", "vector<float>", &b_hadTruth_angleXYZ_vec);
  t->Branch("hadTruth_chi2_vec", "vector<float>", &b_hadTruth_chi2_vec);
  t->Branch("hadTruth_dca_vec", "vector<float>", &b_hadTruth_dca_vec);
  t->Branch("hadTruth_l3D_vec", "vector<float>", &b_hadTruth_l3D_vec);
  t->Branch("hadTruth_l3DSig_vec", "vector<float>", &b_hadTruth_l3DSig_vec);
  t->Branch("hadTruth_legDR_vec", "vector<float>", &b_hadTruth_legDR_vec);
  t->Branch("hadTruth_pdgId_vec", "vector<int>", &b_hadTruth_pdgId_vec);
  t->Branch("hadTruth_dau1_chi2_vec", "vector<float>", &b_hadTruth_dau1_chi2_vec);
  t->Branch("hadTruth_dau1_ipsigXY_vec", "vector<float>", &b_hadTruth_dau1_ipsigXY_vec);
  t->Branch("hadTruth_dau1_ipsigZ_vec", "vector<float>", &b_hadTruth_dau1_ipsigZ_vec);
  t->Branch("hadTruth_dau1_pt_vec", "vector<float>", &b_hadTruth_dau1_pt_vec);
  t->Branch("hadTruth_dau2_chi2_vec", "vector<float>", &b_hadTruth_dau2_chi2_vec);
  t->Branch("hadTruth_dau2_ipsigXY_vec", "vector<float>", &b_hadTruth_dau2_ipsigXY_vec);
  t->Branch("hadTruth_dau2_ipsigZ_vec", "vector<float>", &b_hadTruth_dau2_ipsigZ_vec);
  t->Branch("hadTruth_dau2_pt_vec", "vector<float>", &b_hadTruth_dau2_pt_vec);
  t->Branch("hadTruth_isFrom_vec", "vector<int>", &b_hadTruth_isFrom_vec);

  t->Branch("hadTruth_isFrom_cut_vec", "vector<int>", &b_hadTruth_isFrom_cut_vec);
  t->Branch("hadTruth_isFrom_nc_vec", "vector<int>", &b_hadTruth_isFrom_nc_vec);
  t->Branch("hadTruth_isFrom", &b_hadTruth_isFrom, "hadTruth_isFrom/I");
  t->Branch("hadTruth_x_cut_vec", "vector<float>", &b_hadTruth_x_cut_vec);
  t->Branch("hadTruth_x_nc_vec", "vector<float>", &b_hadTruth_x_nc_vec);
  t->Branch("hadTruth_x", &b_hadTruth_x, "hadTruth_x/F");
  t->Branch("comb_isFrom_vec", "vector<int>", &b_comb_isFrom_vec);
  t->Branch("comb_isFrom", &b_comb_isFrom, "comb_isFrom/I");
  t->Branch("comb_x_vec", "vector<float>", &b_comb_x_vec);
  t->Branch("comb_x", &b_comb_x, "comb_x/F");
  t->Branch("had_isFrom_vec", "vector<int>", &b_had_isFrom_vec);
  t->Branch("had_isFrom_vec_2", "vector<int>", &b_had_isFrom_vec_2);
  t->Branch("had_x_vec", "vector<float>", &b_had_x_vec);
  t->Branch("had_isFrom_dc_vec", "vector<int>", &b_had_isFrom_dc_vec);
  t->Branch("had_isFrom_dc_vec_2", "vector<int>", &b_had_isFrom_dc_vec_2);
  t->Branch("had_x_dc_vec", "vector<float>", &b_had_x_dc_vec);

  t->Branch("had_tlv", "TLorentzVector", &b_had_tlv);
  t->Branch("had_isFrom", &b_had_isFrom, "had_isFrom/I");
  t->Branch("had_isHadJetMatched", &b_had_isHadJetMatched, "had_isHadJetMatched/O");
  t->Branch("had_d", &b_had_d , "had_d/F" );
  t->Branch("had_x", &b_had_x, "had_x/F");
  t->Branch("had_dr", &b_had_dr, "had_dr/F"); // distance between hadron and jet-center
  t->Branch("had_lxy", &b_had_lxy, "had_lxy/F");
  t->Branch("had_lxySig", &b_had_lxySig, "had_lxySig/F");
  t->Branch("had_angleXY", &b_had_angleXY, "had_angleXY/F");
  t->Branch("had_angleXYZ", &b_had_angleXYZ, "had_angleXYZ/F");
  t->Branch("had_chi2", &b_had_chi2, "had_chi2/F");
  t->Branch("had_dca", &b_had_dca, "had_dca/F");

  t->Branch("had_pt", &b_had_pt, "had_pt_had/F");
  t->Branch("had_eta", &b_had_eta, "had_eta/F");
  t->Branch("had_l3D", &b_had_l3D, "had_l3D/F");
  t->Branch("had_l3DSig", &b_had_l3DSig, "had_l3DSig/F");
  t->Branch("had_legDR", &b_had_legDR, "had_legDR/F");
  t->Branch("had_mass", &b_had_mass, "had_mass/F");
  t->Branch("had_pdgId", &b_had_pdgId, "had_pdgId/I");

  t->Branch("had_dau1_chi2", &b_had_dau1_chi2, "had_dau1_chi2/F");
  t->Branch("had_dau1_ipsigXY", &b_had_dau1_ipsigXY, "had_dau1_ipsigXY/F");
  t->Branch("had_dau1_ipsigZ", &b_had_dau1_ipsigZ, "had_dau1_ipsigZ/F");
  t->Branch("had_dau1_pt", &b_had_dau1_pt, "had_dau1_pt/F");

  t->Branch("had_dau2_chi2", &b_had_dau2_chi2, "had_dau2_chi2/F");
  t->Branch("had_dau2_ipsigXY", &b_had_dau2_ipsigXY, "had_dau2_ipsigXY/F");
  t->Branch("had_dau2_ipsigZ", &b_had_dau2_ipsigZ, "had_dau2_ipsigZ/F");
  t->Branch("had_dau2_pt", &b_had_dau2_pt, "had_dau2_pt/F");

  t->Branch("Jet_btagCSVV2", &b_Jet_btagCSVV2, "Jet_btagCSVV2/F");
  t->Branch("Jet_btagDeepB", &b_Jet_btagDeepB, "Jet_btagDeepB/F");
  t->Branch("Jet_btagDeepC", &b_Jet_btagDeepC, "Jet_btagDeepC/F");
  t->Branch("Jet_btagCMVA", &b_Jet_btagCMVA, "Jet_btagCMVA/F");

  t->Branch("Jet_area", &b_Jet_area, "Jet_area/F");
  t->Branch("Jet_pt", &b_Jet_pt, "Jet_pt/F");
  t->Branch("Jet_nConstituents", &b_Jet_nConstituents, "Jet_nConstituents/I");
  t->Branch("Jet_nElectrons", &b_Jet_nElectrons, "Jet_nElectrons/I");
  t->Branch("Jet_nMuons", &b_Jet_nMuons, "Jet_nMuons/I");
}

void vtsAnalyser::ResetBranch() {
  Reset();

  b_hadTruth_nMatched = -1; b_hadTruth_nTrueDau = -1;
  b_hadTruth_isHadFromTsb = -1;
  b_hadTruth_isHadFromTop = false; b_hadTruth_isHadFromW = false; b_hadTruth_isHadFromS = false; b_hadTruth_isHadFromC = false; b_hadTruth_isHadFromB = false;

  // For Test()
  b_hadTruth_isFrom_cut_vec.clear(); b_hadTruth_x_cut_vec.clear();
  b_hadTruth_isFrom_nc_vec.clear(); b_hadTruth_x_nc_vec.clear();
  b_hadTruth_isFrom = -99; b_hadTruth_x = -1;
  b_comb_isFrom_vec.clear(); b_comb_x_vec.clear();
  b_comb_isFrom = -99; b_comb_x = -1;
  b_had_isFrom_vec.clear(); b_had_isFrom_vec_2.clear(); b_had_x_vec.clear();
  b_had_isFrom_dc_vec.clear(); b_had_isFrom_dc_vec_2.clear(); b_had_x_dc_vec.clear();

  b_hadTruth_isFrom_vec.clear(); 
  b_hadTruth_pt_vec.clear(); b_hadTruth_eta_vec.clear(); b_hadTruth_phi_vec.clear(); b_hadTruth_mass_vec.clear();
  b_hadTruth_lxy_vec.clear(); b_hadTruth_lxySig_vec.clear(); b_hadTruth_angleXY_vec.clear(); b_hadTruth_angleXYZ_vec.clear(); b_hadTruth_chi2_vec.clear(); b_hadTruth_dca_vec.clear();
  b_hadTruth_l3D_vec.clear(); b_hadTruth_l3DSig_vec.clear(); b_hadTruth_legDR_vec.clear(); b_hadTruth_pdgId_vec.clear();
  b_hadTruth_dau1_chi2_vec.clear(); b_hadTruth_dau1_ipsigXY_vec.clear(); b_hadTruth_dau1_ipsigZ_vec.clear(); b_hadTruth_dau1_pt_vec.clear();
  b_hadTruth_dau2_chi2_vec.clear(); b_hadTruth_dau2_ipsigXY_vec.clear(); b_hadTruth_dau2_ipsigZ_vec.clear(); b_hadTruth_dau2_pt_vec.clear();

  b_had_tlv.SetPtEtaPhiM(0,0,0,0);
  b_had_isFrom = -99;
  b_had_isHadJetMatched = false;
  b_had_d = -1; b_had_x = -1; b_had_dr = -1;
  b_had_lxy = -1; b_had_lxySig = -1; b_had_angleXY = -1; b_had_angleXYZ = -1; b_had_chi2 = -1; b_had_dca = -1;
  b_had_pt = -1; b_had_eta = -99; b_had_l3D = -1; b_had_l3DSig = -1; b_had_legDR = -1; b_had_mass = -99; b_had_pdgId = -99;
  b_had_dau1_chi2 = -1; b_had_dau1_ipsigXY = -1; b_had_dau1_ipsigZ = -1; b_had_dau1_pt = -1;
  b_had_dau2_chi2 = -1; b_had_dau2_ipsigXY = -1; b_had_dau2_ipsigZ = -1; b_had_dau2_pt = -1;

  b_Jet_btagCSVV2 = -99; b_Jet_btagCMVA = -99; b_Jet_btagDeepB = -99; b_Jet_btagDeepC = -99;
  b_Jet_area = -1; b_Jet_pt = -1; b_Jet_nConstituents = -1; b_Jet_nElectrons = -1; b_Jet_nMuons = -1;

  qjMapForMC_.clear(); tqMC_.clear(); wqMC_.clear(); genJet_.clear(); recoJet_.clear();
}

void vtsAnalyser::MatchingForMC() {
  //Find s quark from Gen Info.  
  for (unsigned int i=0; i<nGenPart; ++i) {
    if (std::abs(GenPart_status[i] - 25) >= 5 ) continue;
    if (abs(GenPart_pdgId[i]) == 3 || abs(GenPart_pdgId[i]) == 5 || abs(GenPart_pdgId[i]) == 4) {
      if ( (abs(GenPart_pdgId[GenPart_genPartIdxMother[i]]) == 6 && GenPart_status[GenPart_genPartIdxMother[i]] == 62 ) ) { 
        tqMC_.push_back(i); 
      }
      if ( (abs(GenPart_pdgId[GenPart_genPartIdxMother[i]]) == 21 && GenPart_status[GenPart_genPartIdxMother[i]] == 21 ) ) {
        wqMC_.push_back(i);
      }
    }
  }
  if (tqMC_.size() < 2) { cout << " >>>>> it's not a tsWbW event. save nothing." << endl; return; }
  //Select quarks originated from t->qW
  auto tq1 = tqMC_[0];
  TLorentzVector tq1_tlv; 
  tq1_tlv.SetPtEtaPhiM(GenPart_pt[tq1], GenPart_eta[tq1], GenPart_phi[tq1], GenPart_mass[tq1]);
  auto tq2 = tqMC_[1];
  TLorentzVector tq2_tlv; 
  tq2_tlv.SetPtEtaPhiM(GenPart_pt[tq2], GenPart_eta[tq2], GenPart_phi[tq2], GenPart_mass[tq2]);
  //Select quarks from W->qq (not yet used)
  int wq1;
  TLorentzVector wq1_tlv;
  int wq2;
  TLorentzVector wq2_tlv;
  if (wqMC_.size() != 0) {
    wq1 = wqMC_[0];
    wq1_tlv.SetPtEtaPhiM(GenPart_pt[wq1], GenPart_eta[wq1], GenPart_phi[wq1], GenPart_mass[wq1]);
    if (wqMC_.size() == 2) {
      wq2 = wqMC_[1];
      wq2_tlv.SetPtEtaPhiM(GenPart_pt[wq2], GenPart_eta[wq2], GenPart_phi[wq2], GenPart_mass[wq2]);
    }
  }
  //Gen Particle & Reco Jet matching
  auto selectedJet = jetSelection();
  if (selectedJet.size() == 0) {
    cout << "Size of selectedJets is zero" << endl;
    return;
  }
  for (unsigned int ij=0; ij<selectedJet.size();++ij) {
    auto j = selectedJet[ij].GetFirstMother();
    TLorentzVector jet_tlv;
    jet_tlv.SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
    auto qIdx = tq1;
    auto dr = tq1_tlv.DeltaR(jet_tlv); 
    auto dr_ = tq2_tlv.DeltaR(jet_tlv); 
    if ( dr_ < dr) {
      dr = dr_;
      qIdx = tq2;
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
  //Reco Jet & Reco KS matching 
  auto selectedJet = jetSelection();
  if (selectedJet.size() == 0) {
    cout << "Size of selectedJets is zero" << endl;
    return;
  }
  for (unsigned int ij=0; ij<selectedJet.size();++ij) {
    auto j = selectedJet[ij].GetFirstMother();
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
/*
      else if (abs(had_pdgId[k]) == 3122) {
        if ( (had_lxy[k]/had_lxyErr[k]) < 3 ) continue;
        if ( had_angleXY[k] < 0.98 ) continue;
        if ( had_chi2[k] > 3 ) continue;
        if ( had_dca[k] > 1 ) continue;
      }
*/
      else continue;
      TLorentzVector had_tlv;
      had_tlv.SetPtEtaPhiM(had_pt[k], had_eta[k], had_phi[k], had_mass[k]);
      if (jet_tlv.DeltaR(had_tlv) < 0.3 && (had_pt[k]/Jet_pt[j]) > 0.15) {
        if (m_isMC) {
          Had.push_back(HadStat(k, had_pdgId[k], qjMapForMC_[j], j, had_pt[k]/Jet_pt[j], jet_tlv.DeltaR(had_tlv), true));
        } else {
          Had.push_back(move(HadStat(k, had_pdgId[k], 0, j, had_pt[k]/Jet_pt[j], jet_tlv.DeltaR(had_tlv), true)));
        }
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

    b_had_x = JetCollection[0][0].x;
    b_had_isFrom = JetCollection[0][0].label;  // -99 : event that can't pass till step4(jet selection) or there is no matching between had and jet, -9 : there is t->qW in the event,but not matched to recoJet, 0 : there is no t->qW in the event (if no t->s and no matching between had-jet, then the event would be -99), +-3 : hadron is from t->sW, +-5 : hadron is from t->bW
    b_had_isHadJetMatched = JetCollection[0][0].isHadJetMatched;
    b_had_dr = JetCollection[0][0].dr;

    b_had_d = GetD(had_pt[idx], had_eta[idx], had_phi[idx], had_mass[idx], had_x[idx], had_y[idx], had_z[idx]);
    b_had_lxy = had_lxy[idx];
    b_had_lxySig = had_lxy[idx]/had_lxyErr[idx];
    b_had_angleXY = had_angleXY[idx];
    b_had_angleXYZ = had_angleXYZ[idx];
    b_had_chi2 = had_chi2[idx];
    b_had_dca = had_dca[idx];
    b_had_pt = had_pt[idx];
    b_had_eta = had_eta[idx];
    b_had_l3D = had_l3D[idx];
    b_had_l3DSig = had_l3D[idx]/had_l3DErr[idx];
    b_had_legDR = had_legDR[idx];
    b_had_mass = had_mass[idx];
    b_had_pdgId = had_pdgId[idx];
    b_had_dau1_chi2 = had_dau1_chi2[idx]; 
    b_had_dau1_ipsigXY = had_dau1_ipsigXY[idx]; 
    b_had_dau1_ipsigZ = had_dau1_ipsigZ[idx]; 
    b_had_dau1_pt = had_dau1_pt[idx];
    b_had_dau2_chi2 = had_dau2_chi2[idx]; 
    b_had_dau2_ipsigXY = had_dau1_ipsigXY[idx]; 
    b_had_dau2_ipsigZ = had_dau1_ipsigZ[idx]; 
    b_had_dau2_pt = had_dau1_pt[idx];

    b_Jet_btagCSVV2 = Jet_btagCSVV2[jidx];
    b_Jet_btagCMVA = Jet_btagCMVA[jidx]; 
    b_Jet_btagDeepB = Jet_btagDeepB[jidx];
    b_Jet_btagDeepC = Jet_btagDeepC[jidx];
    b_Jet_area = Jet_area[jidx]; 
    b_Jet_pt = Jet_pt[jidx];
    b_Jet_nConstituents = Jet_nConstituents[jidx]; 
    b_Jet_nElectrons = Jet_nElectrons[jidx];
    b_Jet_nMuons = Jet_nMuons[jidx];

    b_hadTruth_nMatched = hadTruth_nMatched[idx];
    b_hadTruth_nTrueDau = hadTruth_nTrueDau[idx];
    b_hadTruth_isHadFromTop = hadTruth_isHadFromTop[idx];
    b_hadTruth_isHadFromTsb = hadTruth_isHadFromTsb[idx];
    b_hadTruth_isHadFromW= hadTruth_isHadFromW[idx];
    b_hadTruth_isHadFromS= hadTruth_isHadFromS[idx];
    b_hadTruth_isHadFromC= hadTruth_isHadFromC[idx];
    b_hadTruth_isHadFromB = hadTruth_isHadFromB[idx];
  }
}

//Not yet completed
int vtsAnalyser::Test() {
  for (unsigned int i=0; i<nhad; ++i) {
    if (had_pdgId[i] != 310) continue;
    ++nTotHadKS;
    if (hadTruth_nMatched[i] == 2) {
      b_hadTruth_pt_vec.push_back(had_pt[i]);
      b_hadTruth_eta_vec.push_back(had_eta[i]);
      b_hadTruth_phi_vec.push_back(had_phi[i]);
      b_hadTruth_mass_vec.push_back(had_mass[i]);
      b_hadTruth_lxy_vec.push_back(had_lxy[i]);
      b_hadTruth_lxySig_vec.push_back(had_lxy[i]/had_lxyErr[i]);
      b_hadTruth_angleXY_vec.push_back(had_angleXY[i]);
      b_hadTruth_angleXYZ_vec.push_back(had_angleXYZ[i]);
      b_hadTruth_chi2_vec.push_back(had_chi2[i]);
      b_hadTruth_dca_vec.push_back(had_dca[i]);
      b_hadTruth_l3D_vec.push_back(had_l3D[i]);
      b_hadTruth_l3DSig_vec.push_back(had_l3D[i]/had_l3DErr[i]);
      b_hadTruth_legDR_vec.push_back(had_legDR[i]);
      b_hadTruth_pdgId_vec.push_back(had_pdgId[i]);
      b_hadTruth_dau1_chi2_vec.push_back(had_dau1_chi2[i]);
      b_hadTruth_dau1_ipsigXY_vec.push_back(had_dau1_ipsigXY[i]);
      b_hadTruth_dau1_ipsigZ_vec.push_back(had_dau1_ipsigZ[i]);
      b_hadTruth_dau1_pt_vec.push_back(had_dau1_pt[i]);
      b_hadTruth_dau2_chi2_vec.push_back(had_dau2_chi2[i]);
      b_hadTruth_dau2_ipsigXY_vec.push_back(had_dau1_ipsigXY[i]);
      b_hadTruth_dau2_ipsigZ_vec.push_back(had_dau1_ipsigZ[i]);
      b_hadTruth_dau2_pt_vec.push_back(had_dau1_pt[i]);
      b_hadTruth_isFrom_vec.push_back(hadTruth_isHadFromTsb[i]);
      if (hadTruth_isHadFromTop[i]) ++nRealKSFromTop;
    }
    if (had_dau1_pt[i] >= 0.95 && had_dau2_pt[i] >= 0.95) {
      ++nHadKSDauCut;
    }
    if (hadTruth_nMatched[i] == 2 && (had_dau1_pt[i] >= 0.95 && had_dau2_pt[i] >= 0.95) && hadTruth_isHadFromTop[i]) {
      cout << " [ " << i << " ]  th had ==> dau1 pt : " << had_dau1_pt[i] << " , dau2 pt : " << had_dau2_pt[i] << endl;
      ++nMatchedKS;
    }
    TLorentzVector had_tlv;
    had_tlv.SetPtEtaPhiM(had_pt[i], had_eta[i], had_phi[i], had_mass[i]);
/*
    auto selectedJet = jetSelection();
    if (selectedJet.size() == 0) {
      cout << "Size of selectedJets is zero" << endl;
      return;
    }
    for (unsigned int ij=0; ij<selectedJet.size(); ++ij) { // Loop for jet which pass jetSelection()
      auto j = selectedJet[ij].GetFirstMother();
*/
    for (unsigned int j=0; j<nJet; ++j) { // Loop for all of recoJet
      TLorentzVector jet_tlv;
      jet_tlv.SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
      // For Real KS
      if (hadTruth_nMatched[i] == 2) {
        if (hadTruth_isHadFromTsb[i] == qjMapForMC_[j]) {
          cout << " =====> KS from : " << hadTruth_isHadFromTsb[i] << " , Jet : " << qjMapForMC_[j] << endl;
          ++nRealKSWithJet;
          b_hadTruth_isFrom_nc_vec.push_back(hadTruth_isHadFromTsb[i]);
          b_hadTruth_x_nc_vec.push_back(had_pt[i]/Jet_pt[j]);
          if ((jet_tlv.DeltaR(had_tlv) < 0.3) && (had_pt[i]/Jet_pt[j] > 0.15)) {
            ++nRealKSWithJetAndCut;
            b_hadTruth_isFrom_cut_vec.push_back(hadTruth_isHadFromTsb[i]);
            b_hadTruth_x_cut_vec.push_back(had_pt[i]/Jet_pt[j]);
          }
        } 
      }
      // For KS from Cuts
      if (had_pdgId[i] == 310) {
        if ( (had_lxy[i]/had_lxyErr[i]) < 3 ) continue;
        if ( had_angleXY[i] < 0.98 ) continue;
        if ( had_chi2[i] > 3 ) continue;
        if ( had_dca[i] > 1 ) continue;
        if ( (jet_tlv.DeltaR(had_tlv) < 0.3) && (had_pt[i]/Jet_pt[j] > 0.15)) {
          ++nKSWithCut;
          b_had_isFrom_vec.push_back(hadTruth_isHadFromTsb[i]);
          b_had_isFrom_vec_2.push_back(qjMapForMC_[j]);
          b_had_x_vec.push_back(had_pt[i]/Jet_pt[j]);
          if ( had_dau1_pt[i] >= 0.95 && had_dau2_pt[i] >= 0.95) {
            ++nKSWithDauCut;
            b_had_isFrom_dc_vec.push_back(hadTruth_isHadFromTsb[i]);
            b_had_isFrom_dc_vec_2.push_back(qjMapForMC_[j]);
            b_had_x_dc_vec.push_back(had_pt[i]/Jet_pt[j]);
          } 
        }
      } 
    }
  }
  return 1;
}
