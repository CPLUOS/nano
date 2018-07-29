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
    auto inFile = TFile::Open("/xrootd/store/group/nanoAOD/run2_2016v5/tsw/nanoAOD_111.root", "READ");
//    auto inFile = TFile::Open("/xrootd/store/group/nanoAOD/run2_2016v5/tsW_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6/180611_131219/0000/nanoAOD_000.root");
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
      if (!isMC) { 
        Bool_t isDL = false;
        if (string(inFileName).find("DoubleElectron") != std::string::npos || string(inFileName).find("DoubleMuon") != std::string::npos) isDL = true;
        Bool_t isSLE = (string(inFileName).find("SingleElectron") != std::string::npos);
        Bool_t isSLM = (string(inFileName).find("SingleMuon") != std::string::npos);
        TFile *inFile = TFile::Open(inFileName, "READ");
        TTree *inTree = (TTree*) inFile->Get("Events");
        vtsAnalyser ana(inTree,inTree,inTree,isMC,isDL,isSLE, isSLM);
        string outFileName = outFileDir+"/nanotree_"+to_string(i-3)+".root";
        ana.setOutput(outFileName);
        ana.Loop();
      }
      if (string(inFileName).find("NANOAOD") == std::string::npos) {
        cout << " input file is not tt###j_* sample" << endl;
        if (string(inFileName).find("run2") != std::string::npos) { 
          TFile *inFile = TFile::Open(inFileName, "READ");
          TTree *inTree = (TTree*) inFile->Get("Events");
          vtsAnalyser ana(inTree,inTree,inTree,isMC,false,false,false);
          string outFileName = outFileDir+"/nanotree_"+to_string(i-3)+".root";
          ana.setOutput(outFileName);
          ana.Loop();
          continue;
        } else {
          cout << " input file has to be nanoAOD " << endl;
          continue;
        }
      }
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
      string outFileName = outFileDir+"/nanotree_"+fileName;
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
    if (iev%10000 == 0) cout << iev << "/" << nentries << endl;
    ResetBranch();
    int passedStep = EventSelection();
    if (passedStep >= -1) { // -1 : No event selection, 0 : PV cut, reco lepton cut and so on, 1~4 : step 1 ~ 4
      MatchingForMC();
      HadronAnalysis();
      Test();
      GenHadronAnalysis();
      GenAnalysis();
      RecAnalysis();
    } else if (passedStep >= 0) {
      CollectVar();
    }
    m_tree->Fill();
  }
  cout << "no cut : " << ((float)m_k1/(float)m_ngj)*100 << endl;
  cout << "volcut : " << ((float)m_k2/(float)m_ngj)*100 << endl;
  cout << "daucut : " << ((float)m_k3/(float)m_ngj)*100 << endl;
  cout << "dr cut : " << ((float)m_k4/(float)m_ngj)*100 << endl;
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

  #define Branch_(type, name, suffix) t->Branch(#name, &(b_##name), #name "/" #suffix);
  #define BranchI(name) Branch_(Int_t, name, I)
  #define BranchF(name) Branch_(Float_t, name, F)
  #define BranchO(name) Branch_(Bool_t, name, O)
  #define BranchA_(type, name, size, suffix) t->Branch(#name, &(b_##name), #name"["#size"]/"#suffix);
  #define BranchAI(name, size) BranchA_(Int_t, name, size, I);
  #define BranchAF(name, size) BranchA_(Float_t, name, size, F);
  #define BranchAO(name, size) BranchA_(Bool_t, name, size, O);
  #define BranchV_(type, name) t->Branch(#name, "vector<"#type">", &(b_##name));
  #define BranchVI(name) BranchV_(Int_t, name); 
  #define BranchVF(name) BranchV_(Float_t, name);
  #define BranchVO(name) BranchV_(Bool_t, name);
  #define BranchTLV(name) t->Branch(#name, "TLorentzVector", &(b_##name));

  BranchI(nvertex); BranchI(channel); BranchI(njet); BranchF(met); BranchI(step);
  BranchI(hadTruth_nMatched); BranchI(hadTruth_nTrueDau); 
  BranchO(hadTruth_isHadFromTop); BranchI(hadTruth_isHadFromTsb); BranchO(hadTruth_isHadFromW); BranchO(hadTruth_isHadFromS); BranchO(hadTruth_isHadFromC); BranchO(hadTruth_isHadFromB);

  BranchTLV(had_tlv);
  BranchI(had_pdgId); BranchI(had_isFrom); BranchO(had_isHadJetMatched);
  BranchF(had_d); BranchF(had_pt); BranchF(had_eta); BranchF(had_phi); BranchF(had_mass);
  BranchF(had_lxy); BranchF(had_lxySig); BranchF(had_l3D); BranchF(had_l3DSig); BranchF(had_legDR);
  BranchF(had_angleXY); BranchF(had_angleXYZ); BranchF(had_chi2); BranchF(had_dca);
  BranchF(had_dau1_chi2); BranchF(had_dau1_ipsigXY); BranchF(had_dau1_ipsigZ); BranchF(had_dau1_pt);
  BranchF(had_dau2_chi2); BranchF(had_dau2_ipsigXY); BranchF(had_dau2_ipsigZ); BranchF(had_dau2_pt);
  BranchF(had_x); BranchF(had_dr);
  BranchF(Jet_btagCSVV2); BranchF(Jet_btagDeepB); BranchF(Jet_btagDeepC); BranchF(Jet_btagCMVA);
  BranchF(Jet_area); BranchF(Jet_pt); BranchI(Jet_nConstituents); BranchI(Jet_nElectrons); BranchI(Jet_nMuons);

  // For GenHadronAnalysis()
  BranchVI(genHadron_isGenFrom_vec); BranchVO(genHadron_isGenFromTop_vec); BranchVO(genHadron_inVol_vec);
  BranchVF(genHadron_d_vec); BranchVF(genHadron_pt_vec); BranchVF(genHadron_eta_vec); BranchVF(genHadron_phi_vec); BranchVF(genHadron_mass_vec);
  BranchVF(genHadron_dau1_pt_vec); BranchVF(genHadron_dau1_eta_vec); BranchVF(genHadron_dau1_phi_vec); BranchVI(genHadron_dau1_pdgId_vec);
  BranchVF(genHadron_dau2_pt_vec); BranchVF(genHadron_dau2_eta_vec); BranchVF(genHadron_dau2_phi_vec); BranchVI(genHadron_dau2_pdgId_vec);
  BranchVF(genHadron_x_closest_j_vec); BranchVF(genHadron_dr_closest_j_vec);
  BranchVF(genHadron_x_closest_gj_vec); BranchVF(genHadron_dr_closest_gj_vec); 
  BranchVF(genHadron_x_highest_j_vec); BranchVF(genHadron_dr_highest_j_vec);
  BranchVF(genHadron_x_highest_gj_vec); BranchVF(genHadron_dr_highest_gj_vec); 
  BranchVI(genHadron_isClosestPair_j_vec); BranchVI(genHadron_isClosestPair_gj_vec); BranchVI(genHadron_isHighestPair_j_vec); BranchVI(genHadron_isHighestPair_gj_vec);

  BranchVI(nSJet_vec); BranchVI(nBJet_vec); BranchVI(nGenSJet_vec); BranchVI(nGenBJet_vec);

  // For GenAnalysis()
  BranchVI(GenPart_isGenFrom_vec); BranchVO(GenPart_isGenFromTop_vec); BranchVO(GenPart_isGenFromW_vec); BranchVO(GenPart_isFromKstar_vec);
  BranchVF(GenPart_d_vec); BranchVF(GenPart_pt_vec); BranchVF(GenPart_eta_vec); BranchVF(GenPart_phi_vec); BranchVF(GenPart_mass_vec);
  BranchVF(GenPart_x_closest_j_vec); BranchVF(GenPart_dr_closest_j_vec);
  BranchVF(GenPart_x_closest_gj_vec); BranchVF(GenPart_dr_closest_gj_vec);
  BranchVF(GenPart_x_highest_j_vec); BranchVF(GenPart_dr_highest_j_vec);
  BranchVF(GenPart_x_highest_gj_vec); BranchVF(GenPart_dr_highest_gj_vec);
  BranchVI(GenPart_isClosestPair_j_vec); BranchVI(GenPart_isClosestPair_gj_vec); BranchVI(GenPart_isHighestPair_j_vec); BranchVI(GenPart_isHighestPair_gj_vec);

  // For RecAnalysis()
  BranchVI(hadTruth_pdgId_vec); BranchVI(hadTruth_nMatched_vec); BranchVI(hadTruth_isFrom_vec);
  BranchVO(hadTruth_isHadFromTop_vec); BranchVO(hadTruth_isHadFromW_vec); BranchVO(hadTruth_isHadFromS_vec); BranchVO(hadTruth_isHadFromC_vec); BranchVO(hadTruth_isHadFromB_vec);
  BranchVF(hadTruth_d_vec); BranchVF(hadTruth_pt_vec); BranchVF(hadTruth_eta_vec); BranchVF(hadTruth_phi_vec); BranchVF(hadTruth_mass_vec);
  BranchVF(hadTruth_lxy_vec); BranchVF(hadTruth_lxySig_vec); BranchVF(hadTruth_l3D_vec); BranchVF(hadTruth_l3DSig_vec); BranchVF(hadTruth_legDR_vec);
  BranchVF(hadTruth_angleXY_vec); BranchVF(hadTruth_angleXYZ_vec); BranchVF(hadTruth_chi2_vec); BranchVF(hadTruth_dca_vec);
  BranchVF(hadTruth_dau1_chi2_vec); BranchVF(hadTruth_dau1_ipsigXY_vec); BranchVF(hadTruth_dau1_ipsigZ_vec); BranchVF(hadTruth_dau1_pt_vec);
  BranchVF(hadTruth_dau2_chi2_vec); BranchVF(hadTruth_dau2_ipsigXY_vec); BranchVF(hadTruth_dau2_ipsigZ_vec); BranchVF(hadTruth_dau2_pt_vec);
  BranchVF(hadTruth_x_closest_j_vec); BranchVF(hadTruth_dr_closest_j_vec);
  BranchVF(hadTruth_x_closest_gj_vec); BranchVF(hadTruth_dr_closest_gj_vec);
  BranchVF(hadTruth_x_highest_j_vec); BranchVF(hadTruth_dr_highest_j_vec);
  BranchVF(hadTruth_x_highest_gj_vec); BranchVF(hadTruth_dr_highest_gj_vec);
  BranchVI(hadTruth_isClosestPair_j_vec); BranchVI(hadTruth_isClosestPair_gj_vec); BranchVI(hadTruth_isHighestPair_j_vec); BranchVI(hadTruth_isHighestPair_gj_vec);

  // For CollectVar()
  BranchF(MET_pt); BranchF(MET_phi); BranchF(MET_sumEt);
  BranchI(lep1_pid); BranchI(lep2_pid);
  BranchTLV(lep1); BranchTLV(lep2);
  BranchTLV(dilep);
}

void vtsAnalyser::ResetBranch() {
  Reset();

  b_hadTruth_nMatched = -1; b_hadTruth_nTrueDau = -1;
  b_hadTruth_isHadFromTsb = -1;
  b_hadTruth_isHadFromTop = false; b_hadTruth_isHadFromW = false; b_hadTruth_isHadFromS = false; b_hadTruth_isHadFromC = false; b_hadTruth_isHadFromB = false;

  // For GenHadronAnalysis()
  b_genHadron_isGenFrom_vec.clear(); b_genHadron_isGenFromTop_vec.clear(); b_genHadron_inVol_vec.clear();
  b_genHadron_d_vec.clear(); b_genHadron_pt_vec.clear(); b_genHadron_eta_vec.clear(); b_genHadron_phi_vec.clear(); b_genHadron_mass_vec.clear();
  b_genHadron_dau1_pt_vec.clear(); b_genHadron_dau1_eta_vec.clear(); b_genHadron_dau1_phi_vec.clear(); b_genHadron_dau1_pdgId_vec.clear();
  b_genHadron_dau2_pt_vec.clear(); b_genHadron_dau2_eta_vec.clear(); b_genHadron_dau2_phi_vec.clear(); b_genHadron_dau2_pdgId_vec.clear();
  b_genHadron_x_closest_j_vec.clear(); b_genHadron_dr_closest_j_vec.clear();
  b_genHadron_x_closest_gj_vec.clear(); b_genHadron_dr_closest_gj_vec.clear();
  b_genHadron_x_highest_j_vec.clear(); b_genHadron_dr_highest_j_vec.clear();
  b_genHadron_x_highest_gj_vec.clear(); b_genHadron_dr_highest_gj_vec.clear();
  b_genHadron_isClosestPair_j_vec.clear(); b_genHadron_isClosestPair_gj_vec.clear(); b_genHadron_isHighestPair_j_vec.clear(); b_genHadron_isHighestPair_gj_vec.clear();

  b_nSJet_vec.clear(); b_nBJet_vec.clear(); b_nGenSJet_vec.clear(); b_nGenBJet_vec.clear();

  // For GenAnalysis()
  b_GenPart_isGenFrom_vec.clear(); b_GenPart_isGenFromTop_vec.clear(); b_GenPart_isGenFromW_vec.clear(); b_GenPart_isFromKstar_vec.clear();
  b_GenPart_d_vec.clear(); b_GenPart_pt_vec.clear(); b_GenPart_eta_vec.clear(); b_GenPart_phi_vec.clear(); b_GenPart_mass_vec.clear();
  b_GenPart_x_closest_j_vec.clear(); b_GenPart_dr_closest_j_vec.clear();
  b_GenPart_x_closest_gj_vec.clear(); b_GenPart_dr_closest_gj_vec.clear();
  b_GenPart_x_highest_j_vec.clear(); b_GenPart_dr_highest_j_vec.clear();
  b_GenPart_x_highest_gj_vec.clear(); b_GenPart_dr_highest_gj_vec.clear();
  b_GenPart_isClosestPair_j_vec.clear(); b_GenPart_isClosestPair_gj_vec.clear(); b_GenPart_isHighestPair_j_vec.clear(); b_GenPart_isHighestPair_gj_vec.clear();

  // For RecAnalysis()
  b_hadTruth_pdgId_vec.clear(); b_hadTruth_nMatched_vec.clear(); b_hadTruth_isFrom_vec.clear(); 
  b_hadTruth_isHadFromTop_vec.clear(); b_hadTruth_isHadFromW_vec.clear(); b_hadTruth_isHadFromS_vec.clear(); b_hadTruth_isHadFromC_vec.clear(); b_hadTruth_isHadFromB_vec.clear(); 
  b_hadTruth_d_vec.clear(); b_hadTruth_pt_vec.clear(); b_hadTruth_eta_vec.clear(); b_hadTruth_phi_vec.clear(); b_hadTruth_mass_vec.clear();
  b_hadTruth_lxy_vec.clear(); b_hadTruth_lxySig_vec.clear(); b_hadTruth_l3D_vec.clear(); b_hadTruth_l3DSig_vec.clear(); b_hadTruth_legDR_vec.clear();
  b_hadTruth_angleXY_vec.clear(); b_hadTruth_angleXYZ_vec.clear(); b_hadTruth_chi2_vec.clear(); b_hadTruth_dca_vec.clear();
  b_hadTruth_dau1_chi2_vec.clear(); b_hadTruth_dau1_ipsigXY_vec.clear(); b_hadTruth_dau1_ipsigZ_vec.clear(); b_hadTruth_dau1_pt_vec.clear();
  b_hadTruth_dau2_chi2_vec.clear(); b_hadTruth_dau2_ipsigXY_vec.clear(); b_hadTruth_dau2_ipsigZ_vec.clear(); b_hadTruth_dau2_pt_vec.clear();
  b_hadTruth_x_closest_j_vec.clear(); b_hadTruth_dr_closest_j_vec.clear();
  b_hadTruth_x_closest_gj_vec.clear(); b_hadTruth_dr_closest_gj_vec.clear();
  b_hadTruth_x_highest_j_vec.clear(); b_hadTruth_dr_highest_j_vec.clear();
  b_hadTruth_x_highest_gj_vec.clear(); b_hadTruth_dr_highest_gj_vec.clear();
  b_hadTruth_isClosestPair_j_vec.clear(); b_hadTruth_isClosestPair_gj_vec.clear(); b_hadTruth_isHighestPair_j_vec.clear(); b_hadTruth_isHighestPair_gj_vec.clear();

  // For CollectVar()
  b_MET_pt = -1; b_MET_phi = -99; b_MET_sumEt = -1;

  b_had_tlv.SetPtEtaPhiM(0,0,0,0);
  b_had_pdgId = -99; b_had_isFrom = -99; b_had_isHadJetMatched = false;
  b_had_d = -1; b_had_x = -1; b_had_dr = -1;
  b_had_pt = -1; b_had_eta = -99; b_had_phi = -99; b_had_mass = -99;
  b_had_lxy = -1; b_had_lxySig = -1; b_had_l3D = -1; b_had_l3DSig = -1; b_had_legDR = -1;
  b_had_angleXY = -1; b_had_angleXYZ = -1; b_had_chi2 = -1; b_had_dca = -1;
  b_had_dau1_chi2 = -1; b_had_dau1_ipsigXY = -1; b_had_dau1_ipsigZ = -1; b_had_dau1_pt = -1;
  b_had_dau2_chi2 = -1; b_had_dau2_ipsigXY = -1; b_had_dau2_ipsigZ = -1; b_had_dau2_pt = -1;
  b_had_x = -1; b_had_dr = -1;
  b_Jet_btagCSVV2 = -99; b_Jet_btagCMVA = -99; b_Jet_btagDeepB = -99; b_Jet_btagDeepC = -99;
  b_Jet_area = -1; b_Jet_pt = -1; b_Jet_nConstituents = -1; b_Jet_nElectrons = -1; b_Jet_nMuons = -1;

  m_tqMC.clear(); m_wqMC.clear(); m_qjMapForMC.clear(); m_qgjMapForMC.clear(); m_recoJet.clear(); m_genJet.clear();
}

void vtsAnalyser::MatchingForMC() {
  //Find s quark from Gen Info. 

  // status is case of pythia 
  // top quark statusFlag(status) : 10497(62) 
  // s/b quark from top  statusFlag(status) : 22913(23)
  // W boson (t->W->q) statusFlag(status) : 14721(22)
  // W boson 1 (t->W1->W2->q) statusFlag(status) : 4481(22)  W boson 2 (t->W1->W2->q)  statusFlag(status) : 10497(52)
  // quark from W boson statusFlag(status) : 22913(23) or 4481(23)
  for (unsigned int i=0; i<nGenPart; ++i) {
    if (std::abs(GenPart_status[i]) != 23 ) continue;
    if (abs(GenPart_pdgId[i]) == 3 || abs(GenPart_pdgId[i]) == 5 || abs(GenPart_pdgId[i]) == 4) {
      if ( (abs(GenPart_pdgId[GenPart_genPartIdxMother[i]]) == 6 && GenPart_status[GenPart_genPartIdxMother[i]] == 62 ) ) { 
        m_tqMC.push_back(i);
      }
      if ( (abs(GenPart_pdgId[GenPart_genPartIdxMother[i]]) == 24 && (GenPart_status[GenPart_genPartIdxMother[i]] == 22 || GenPart_status[GenPart_genPartIdxMother[i]] == 52)) ) {
        m_wqMC.push_back(i);
      }
    }
  }
  if (m_tqMC.size() < 2) { cout << " >>>>> it's not a tsWbW event. save nothing." << endl; return; }
  //Select quarks originated from t->qW
  auto tq1 = m_tqMC[0];
  TLorentzVector tq1_tlv; 
  tq1_tlv.SetPtEtaPhiM(GenPart_pt[tq1], GenPart_eta[tq1], GenPart_phi[tq1], GenPart_mass[tq1]);
  auto tq2 = m_tqMC[1];
  TLorentzVector tq2_tlv; 
  tq2_tlv.SetPtEtaPhiM(GenPart_pt[tq2], GenPart_eta[tq2], GenPart_phi[tq2], GenPart_mass[tq2]);
  //Select quarks from W->qq (not yet used)
  int wq1;
  TLorentzVector wq1_tlv;
  int wq2;
  TLorentzVector wq2_tlv;
  if (m_wqMC.size() != 0) {
    wq1 = m_wqMC[0];
    wq1_tlv.SetPtEtaPhiM(GenPart_pt[wq1], GenPart_eta[wq1], GenPart_phi[wq1], GenPart_mass[wq1]);
    if (m_wqMC.size() == 2) {
      wq2 = m_wqMC[1];
      wq2_tlv.SetPtEtaPhiM(GenPart_pt[wq2], GenPart_eta[wq2], GenPart_phi[wq2], GenPart_mass[wq2]);
    }
  }
  // Gen Partcle & Reco Jet Matching
  std::vector<std::pair<unsigned int, float>> drReco1; std::vector<std::pair<unsigned int, float>> drReco2;
  auto selectedJet = jetSelection();
  if (selectedJet.size() == 0) {
    cout << "Size of selectedJets is zero" << endl;
    return;
  }
  for (unsigned int ij=0; ij<selectedJet.size();++ij) {
    auto j = selectedJet[ij].GetFirstMother();
    TLorentzVector jet_tlv;
    jet_tlv.SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
    auto dr1 = tq1_tlv.DeltaR(jet_tlv); auto dr2 = tq2_tlv.DeltaR(jet_tlv);
    drReco1.push_back({j, dr1}); drReco2.push_back({j, dr2});
  }
  sort(drReco1.begin(), drReco1.end(), [] (std::pair<unsigned int, float> a, std::pair<unsigned int, float> b) { return (a.second < b.second); } );
  sort(drReco2.begin(), drReco2.end(), [] (std::pair<unsigned int, float> a, std::pair<unsigned int, float> b) { return (a.second < b.second); } );
  m_qjMapForMC[drReco1[0].first] = GenPart_pdgId[tq1];
  m_qjMapForMC[drReco2[0].first] = GenPart_pdgId[tq2];

  // Gen Particle & Gen Jet matching
  std::vector<std::pair<unsigned int, float>> drGen1; std::vector<std::pair<unsigned int, float>> drGen2;
  for (unsigned int j=0; j<nGenJet; ++j) {
    TLorentzVector jet_tlv;
    jet_tlv.SetPtEtaPhiM(GenJet_pt[j], GenJet_eta[j], GenJet_phi[j], GenJet_mass[j]);
    auto dr1 = tq1_tlv.DeltaR(jet_tlv); auto dr2 = tq2_tlv.DeltaR(jet_tlv);
    drGen1.push_back({j, dr1}); drGen2.push_back({j, dr2});
  }
  sort(drGen1.begin(), drGen1.end(), [] (std::pair<unsigned int, float> a, std::pair<unsigned int, float> b) { return (a.second < b.second); } );
  sort(drGen2.begin(), drGen2.end(), [] (std::pair<unsigned int, float> a, std::pair<unsigned int, float> b) { return (a.second < b.second); } );
  if ( GenPart_pdgId[tq1] == GenJet_partonFlavour[drGen1[0].first] ) m_qgjMapForMC[drGen1[0].first] = GenPart_pdgId[tq1];
  if ( GenPart_pdgId[tq2] == GenJet_partonFlavour[drGen2[0].first] ) m_qgjMapForMC[drGen2[0].first] = GenPart_pdgId[tq2];
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
      if (jet_tlv.DeltaR(had_tlv) < m_jetConeSize && (had_pt[k]/Jet_pt[j]) > 0.15) {
        if (m_isMC) {
          Had.push_back(HadStat(k, had_pdgId[k], m_qjMapForMC[j], j, had_pt[k]/Jet_pt[j], jet_tlv.DeltaR(had_tlv), true));
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
    b_had_isFrom = JetCollection[0][0].label;  // -99 : event that can't pass till step4(jet selection) or there is no matching between had and jet, 0 : there is t->qW in the event,but not matched to recoJet (if no t->s and no matching between had-jet, then the event would be -99), +-3 : hadron is from t->sW, +-5 : hadron is from t->bW
    b_had_isHadJetMatched = JetCollection[0][0].isHadJetMatched;
    b_had_d = GetD(had_pt[idx], had_eta[idx], had_phi[idx], had_mass[idx], had_x[idx], had_y[idx], had_z[idx]);
    b_had_x = JetCollection[0][0].x;
    b_had_dr = JetCollection[0][0].dr;
    b_had_pt = had_pt[idx];
    b_had_eta = had_eta[idx];
    b_had_phi = had_phi[idx];
    b_had_mass = had_mass[idx];
    b_had_lxy = had_lxy[idx];
    b_had_lxySig = had_lxy[idx]/had_lxyErr[idx];
    b_had_angleXY = had_angleXY[idx];
    b_had_angleXYZ = had_angleXYZ[idx];
    b_had_chi2 = had_chi2[idx];
    b_had_dca = had_dca[idx];
    b_had_l3D = had_l3D[idx];
    b_had_l3DSig = had_l3D[idx]/had_l3DErr[idx];
    b_had_legDR = had_legDR[idx];
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

bool vtsAnalyser::isGenFrom(int count, int idx, int & isFrom, bool & isFromTop, bool & isFromW, bool & isFromKstar)
{
  isFrom = -99;
  isFromTop = false;
  isFromW = false;
  if ((abs(GenPart_pdgId[idx]) == 3 || abs(GenPart_pdgId[idx]) == 4 || abs(GenPart_pdgId[idx]) == 5) && GenPart_status[idx] == 23) {
    if (abs(GenPart_pdgId[GenPart_genPartIdxMother[idx]]) == 6 && GenPart_status[GenPart_genPartIdxMother[idx]] == 62) {
      isFrom = GenPart_pdgId[idx];
      isFromTop = true;
      isFromW = false;
      return true;
    } else if (abs(GenPart_pdgId[GenPart_genPartIdxMother[idx]]) == 24 && (GenPart_status[GenPart_genPartIdxMother[idx]] == 22 || GenPart_status[GenPart_genPartIdxMother[idx]] == 52)) {
      isFrom = GenPart_pdgId[idx];
      isFromW = true;
      int midx = 0;
      if (GenPart_status[GenPart_genPartIdxMother[idx]] == 22) midx = GenPart_genPartIdxMother[GenPart_genPartIdxMother[idx]];
      else if (GenPart_status[GenPart_genPartIdxMother[idx]] == 52) midx = GenPart_genPartIdxMother[GenPart_genPartIdxMother[GenPart_genPartIdxMother[idx]]];
      if (abs(GenPart_pdgId[midx]) == 6 && GenPart_status[midx] == 62) {
        isFromTop = true;
      } else {
        isFromTop = false;
      }
      return true;
    }
  }
  std::vector<int> kstar = {313, 323, 315, 325, 317, 327, 319, 329};
  if (abs(GenPart_pdgId[idx]) == 310) { 
    if (abs(GenPart_pdgId[GenPart_genPartIdxMother[idx]]) == 311) {
      if (std::find(kstar.begin(), kstar.end(), abs(GenPart_pdgId[GenPart_genPartIdxMother[GenPart_genPartIdxMother[idx]]])) != kstar.end()) isFromKstar = true;
      else isFromKstar = false;
    } else if (std::find(kstar.begin(), kstar.end(), abs(GenPart_pdgId[GenPart_genPartIdxMother[idx]])) != kstar.end()) {
      isFromKstar = true;
    } else isFromKstar = false;
  }
  auto mom = GenPart_genPartIdxMother[idx];
  ++count;
  if (mom != -1) return isGenFrom( count, mom, isFrom, isFromTop, isFromW, isFromKstar);
  else return false;
}

void vtsAnalyser::Test() {
  for (unsigned int j=0; j<nGenJet; ++j) {
    if (abs(GenJet_partonFlavour[j]) != 3) continue;
    m_ngj = m_ngj +1;
    TLorentzVector gjet_tlv;
    gjet_tlv.SetPtEtaPhiM(GenJet_pt[j], GenJet_eta[j], GenJet_phi[j], GenJet_mass[j]);

    auto maxPt = 0.; auto dr =0.; auto idx =0;
    for (unsigned int i=0; i<ngenHadron; ++i) {
      if (abs(genHadron_pdgId[i]) != 310) continue;
      if (genHadron_isGenHadFromTsb[i] != GenJet_partonFlavour[j]) continue;
      TLorentzVector ghad_tlv;
      ghad_tlv.SetPtEtaPhiM(genHadron_pt[i], genHadron_eta[i], genHadron_phi[i], genHadron_mass[i]);
      auto dRtmp = gjet_tlv.DeltaR(ghad_tlv);
      auto x = genHadron_pt[i]/GenJet_pt[j];

      if ( x < m_xcut) continue;
      if (genHadron_pt[i] > maxPt) {
        maxPt = genHadron_pt[i];
        dr = dRtmp;
        idx = i;
      }
      if (maxPt == 0) continue;
      m_k1 = m_k1+1;

      if (genHadron_inVol[idx] != 1) continue;
      m_k2 = m_k2+1;

      if (abs(genHadron_dau1_pdgId[idx]) != 211 || abs(genHadron_dau2_pdgId[idx]) != 211) continue;
      m_k3 = m_k3+1;

      if (dr > m_jetConeSize) continue;
      m_k4 = m_k4+1;
    }
  }
}

void vtsAnalyser::GenHadronAnalysis() {
  std::vector<std::pair<unsigned int, float>> drRec; std::vector<std::pair<unsigned int, float>> xRec;
  std::vector<std::pair<unsigned int, float>> drGen;  std::vector<std::pair<unsigned int, float>> xGen;
  int nGenHadKS = 0;
  b_genHadron_dr_closest_j_vec.resize(ngenHadron);
  b_genHadron_x_closest_j_vec.resize(ngenHadron);
  b_genHadron_dr_closest_gj_vec.resize(ngenHadron); 
  b_genHadron_x_closest_gj_vec.resize(ngenHadron);
  b_genHadron_dr_highest_j_vec.resize(ngenHadron); 
  b_genHadron_x_highest_j_vec.resize(ngenHadron);
  b_genHadron_dr_highest_gj_vec.resize(ngenHadron); 
  b_genHadron_x_highest_gj_vec.resize(ngenHadron);
  for (unsigned int i=0; i<ngenHadron; ++i) {
    if (abs(genHadron_pdgId[i]) != 310) continue;
    ++nGenHadKS;
    b_genHadron_isGenFrom_vec.push_back(genHadron_isGenHadFromTsb[i]);
    b_genHadron_isGenFromTop_vec.push_back(genHadron_isGenHadFromTop[i]);
    b_genHadron_inVol_vec.push_back(genHadron_inVol[i]);
    b_genHadron_d_vec.push_back(GetD(genHadron_pt[i], genHadron_eta[i], genHadron_phi[i], genHadron_mass[i], genHadron_x[i], genHadron_y[i], genHadron_z[i]));
    b_genHadron_pt_vec.push_back(genHadron_pt[i]);
    b_genHadron_eta_vec.push_back(genHadron_eta[i]);
    b_genHadron_phi_vec.push_back(genHadron_phi[i]);
    b_genHadron_mass_vec.push_back(genHadron_mass[i]);
    b_genHadron_dau1_pt_vec.push_back(genHadron_dau1_pt[i]);
    b_genHadron_dau1_eta_vec.push_back(genHadron_dau1_eta[i]);
    b_genHadron_dau1_phi_vec.push_back(genHadron_dau1_phi[i]);
    b_genHadron_dau1_pdgId_vec.push_back(genHadron_dau1_pdgId[i]);
    b_genHadron_dau2_pt_vec.push_back(genHadron_dau2_pt[i]);
    b_genHadron_dau2_eta_vec.push_back(genHadron_dau2_eta[i]);
    b_genHadron_dau2_phi_vec.push_back(genHadron_dau2_phi[i]);
    b_genHadron_dau2_pdgId_vec.push_back(genHadron_dau2_pdgId[i]);

    TLorentzVector gen_tlv;
    gen_tlv.SetPtEtaPhiM(genHadron_pt[i], genHadron_eta[i], genHadron_phi[i], genHadron_mass[i]);

    for (unsigned int j=0; j<nJet; ++j) { // Loop for all of recoJet
      TLorentzVector jet_tlv;
      jet_tlv.SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
      auto dr = jet_tlv.DeltaR(gen_tlv);
      auto x = genHadron_pt[i]/Jet_pt[j];
      if (nGenHadKS==1) {
        drRec.push_back({i, dr});
        xRec.push_back({i, x});
        if (abs(m_qjMapForMC[j]) == 3) b_nSJet_vec.push_back(j);
        else if (abs(m_qjMapForMC[j]) == 5) b_nBJet_vec.push_back(j);
      } else {
        if (genHadron_isGenHadFromTsb[i] == m_qjMapForMC[j]) {
          if (dr < drRec[j].second ) { // Find the closest gen KS for each Jets
            drRec[j] = {i, dr};
          } /*else if (dr == drRec[j].second) { // not completed
          }*/
          if (x > xRec[j].second ) { // Find the highest x gen KS for each Jets
            xRec[j] = {i, x};
          } /*else if (x == xRec[j].second) { // not completed
          }*/
        }
      }
    }
    for (unsigned int j=0; j<nGenJet; ++j) { // Loop for all of recoJet
      TLorentzVector jet_tlv;
      jet_tlv.SetPtEtaPhiM(GenJet_pt[j], GenJet_eta[j], GenJet_phi[j], GenJet_mass[j]);
      auto dr = jet_tlv.DeltaR(gen_tlv);
      auto x = genHadron_pt[i]/GenJet_pt[j];
      if (nGenHadKS==1) {
        drGen.push_back({i, dr});
        xGen.push_back({i, x});
        if (abs(m_qgjMapForMC[j]) == 3) b_nGenSJet_vec.push_back(j);
        else if (abs(m_qgjMapForMC[j]) == 5) b_nGenBJet_vec.push_back(j);
      } else {
        if (genHadron_isGenHadFromTsb[i] == m_qgjMapForMC[j]) {
          if (dr < drGen[j].second ) { // Find the closest gen KS for each GenJets
            drGen[j] = {i, dr};
          } /*else if (dr == drGen[j].second) { // not completed
          }*/
          if (x > xGen[j].second ) { // Find the highest x gen KS for each GenJets
            xGen[j] = {i, x};
          } /*else if (x == xGen[j].second) { // not completed
          }*/
        }
      }
    }
    b_genHadron_isClosestPair_j_vec.push_back(i); b_genHadron_isClosestPair_gj_vec.push_back(i);
    b_genHadron_isHighestPair_j_vec.push_back(i); b_genHadron_isHighestPair_gj_vec.push_back(i);
  }
  if (nGenHadKS !=0) {
    //  Save x and dr for the most closest(highest) KS for each jet
    for (unsigned int j=0; j<nJet; ++j) {
      int cRecIdx = drRec[j].first; int hRecIdx = xRec[j].first;
      TLorentzVector hgen_tlv; TLorentzVector hjet_tlv;
      hgen_tlv.SetPtEtaPhiM(genHadron_pt[hRecIdx], genHadron_eta[hRecIdx], genHadron_phi[hRecIdx], genHadron_mass[hRecIdx]);
      hjet_tlv.SetPtEtaPhiM(GenJet_pt[j], GenJet_eta[j], GenJet_phi[j], GenJet_mass[j]);
      b_genHadron_dr_closest_j_vec[cRecIdx] = drRec[j].second;
      b_genHadron_x_closest_j_vec[cRecIdx] = genHadron_pt[cRecIdx] / Jet_pt[j];
      b_genHadron_dr_highest_j_vec[hRecIdx] = hgen_tlv.DeltaR(hjet_tlv);
      b_genHadron_x_highest_j_vec[hRecIdx] = xRec[j].second;
//      cout << "highest x : ( " << xRec[j].first << " , " << xRec[j].second << " , " << hgen_tlv.DeltaR(hjet_tlv) << " ) " << endl;

    }
    for (unsigned int j=0; j<nGenJet; ++j) {
      int cGenIdx = drGen[j].first; int hGenIdx = xGen[j].first;
      TLorentzVector hgen_tlv; TLorentzVector hjet_tlv;
      hgen_tlv.SetPtEtaPhiM(genHadron_pt[hGenIdx], genHadron_eta[hGenIdx], genHadron_phi[hGenIdx], genHadron_mass[hGenIdx]);
      hjet_tlv.SetPtEtaPhiM(GenJet_pt[j], GenJet_eta[j], GenJet_phi[j], GenJet_mass[j]);
      b_genHadron_dr_closest_gj_vec[cGenIdx] = drGen[j].second;
      b_genHadron_x_closest_gj_vec[cGenIdx] = genHadron_pt[cGenIdx] / Jet_pt[j];
      b_genHadron_dr_highest_gj_vec[hGenIdx] = hgen_tlv.DeltaR(hjet_tlv);
      b_genHadron_x_highest_gj_vec[hGenIdx] = xGen[j].second;
//      cout << "highest x : ( " << xGen[j].first << " , " << xGen[j].second << " , " << hgen_tlv.DeltaR(hjet_tlv) << " ) " << endl;
    }
    // Give flag ( == index) for the most closest(highest) KS-jet pair per event
    std::sort(drRec.begin(), drRec.end(), [] (std::pair<unsigned int, float> a, std::pair<unsigned int, float> b) { return (a.second < b.second); } );
    std::sort(drGen.begin(), drGen.end(), [] (std::pair<unsigned int, float> a, std::pair<unsigned int, float> b) { return (a.second < b.second); } );
    std::sort(xRec.begin(), xRec.end(), [] (std::pair<unsigned int, float> a, std::pair<unsigned int, float> b) { return (a.second > b.second); } );
    std::sort(xGen.begin(), xGen.end(), [] (std::pair<unsigned int, float> a, std::pair<unsigned int, float> b) { return (a.second > b.second); } );
    int cRecIdx = drRec[0].first; int cGenIdx = drGen[0].first; int hRecIdx = xRec[0].first; int hGenIdx = xGen[0].first;
    std::replace_if(b_genHadron_isClosestPair_j_vec.begin(), b_genHadron_isClosestPair_j_vec.end(), [&] (int a) { return (a != cRecIdx); }, -1);
    std::replace_if(b_genHadron_isClosestPair_gj_vec.begin(), b_genHadron_isClosestPair_gj_vec.end(), [&] (int a) { return (a != cGenIdx); }, -1);
    std::replace_if(b_genHadron_isHighestPair_j_vec.begin(), b_genHadron_isHighestPair_j_vec.end(), [&] (int a) { return (a != hRecIdx); }, -1);
    std::replace_if(b_genHadron_isHighestPair_gj_vec.begin(), b_genHadron_isHighestPair_gj_vec.end(), [&] (int a) { return (a != hGenIdx); }, -1);
  }
}

void vtsAnalyser::GenAnalysis() {
  std::vector<std::pair<unsigned int, float>> drRec; std::vector<std::pair<unsigned int, float>> xRec;
  std::vector<std::pair<unsigned int, float>> drGen;  std::vector<std::pair<unsigned int, float>> xGen;
  int nGenKS = 0;
  b_GenPart_dr_closest_j_vec.resize(nGenPart);
  b_GenPart_x_closest_j_vec.resize(nGenPart);
  b_GenPart_dr_closest_gj_vec.resize(nGenPart);
  b_GenPart_x_closest_gj_vec.resize(nGenPart);
  b_GenPart_dr_highest_j_vec.resize(nGenPart);
  b_GenPart_x_highest_j_vec.resize(nGenPart);
  b_GenPart_dr_highest_gj_vec.resize(nGenPart);
  b_GenPart_x_highest_gj_vec.resize(nGenPart);
  for (unsigned int i=0; i<nGenPart; ++i) {
    if (abs(GenPart_pdgId[i]) != 310) continue;
    ++nGenKS;
    int count = 0; int isFrom = 0;
    bool isFromTop = false; bool isFromW = false; bool isFromKstar = false;
    TLorentzVector gen_tlv;
    gen_tlv.SetPtEtaPhiM(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i]);
    isGenFrom(count, i, isFrom, isFromTop, isFromW, isFromKstar); // MotherTracking
    b_GenPart_isFromKstar_vec.push_back(isFromKstar);
    b_GenPart_isGenFrom_vec.push_back(isFrom);
    b_GenPart_isGenFromTop_vec.push_back(isFromTop);
    b_GenPart_isGenFromW_vec.push_back(isFromW);
    b_GenPart_d_vec.push_back(GetD(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i], gen_tlv.X(), gen_tlv.Y(), gen_tlv.Z()));
    b_GenPart_pt_vec.push_back(GenPart_pt[i]);
    b_GenPart_eta_vec.push_back(GenPart_eta[i]);
    b_GenPart_phi_vec.push_back(GenPart_phi[i]);
    b_GenPart_mass_vec.push_back(GenPart_mass[i]);

    for (unsigned int j=0; j<nJet; ++j) { // Loop for all of recoJet
      TLorentzVector jet_tlv;
      jet_tlv.SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
      auto dr = jet_tlv.DeltaR(gen_tlv);
      auto x = GenPart_pt[i]/Jet_pt[j];
      if (nGenKS==1) {
        drRec.push_back({i, dr});
        xRec.push_back({i, x});
      } else { 
        if (isFrom == m_qjMapForMC[j]) {
          if (dr < drRec[j].second ) { // Find the closest gen KS for each Jets
            drRec[j] = {i, dr};
          } /*else if (dr == drRec[j].second) { // not completed
          }*/
          if (x > xRec[j].second ) { // Find the highest x gen KS for each Jets
            xRec[j] = {i, x};
          } /*else if (x == xRec[j].second) { // not completed
          }*/
        }
      }
    }
    for (unsigned int j=0; j<nGenJet; ++j) { // Loop for all of recoJet
      TLorentzVector jet_tlv;
      jet_tlv.SetPtEtaPhiM(GenJet_pt[j], GenJet_eta[j], GenJet_phi[j], GenJet_mass[j]);
      auto dr = jet_tlv.DeltaR(gen_tlv);
      auto x = GenPart_pt[i]/GenJet_pt[j];
      if (nGenKS==1) {
        drGen.push_back({i, dr});
        xGen.push_back({i, x});
      } else {
        if (isFrom == m_qgjMapForMC[j]) {
          if (dr < drGen[j].second ) { // Find the closest gen KS for each GenJets
            drGen[j] = {i, dr};
          } /*else if (dr == drGen[j].second) { // not completed
          }*/
          if (x > xGen[j].second ) { // Find the highest x gen KS for each GenJets
            xGen[j] = {i, x};
          } /*else if (x == xGen[j].second) { // not completed
          }*/
        }
      }
    }
    b_GenPart_isClosestPair_j_vec.push_back(i); b_GenPart_isClosestPair_gj_vec.push_back(i);
    b_GenPart_isHighestPair_j_vec.push_back(i); b_GenPart_isHighestPair_gj_vec.push_back(i);
  }
  if (nGenKS !=0) {
    //  Save x and dr for the most closest(highest) KS for each jet
    for (unsigned int j=0; j<nJet; ++j) {
      int cRecIdx = drRec[j].first; int hRecIdx = xRec[j].first;
      TLorentzVector hgen_tlv; TLorentzVector hjet_tlv;
      hgen_tlv.SetPtEtaPhiM(GenPart_pt[hRecIdx], GenPart_eta[hRecIdx], GenPart_phi[hRecIdx], GenPart_mass[hRecIdx]);
      hjet_tlv.SetPtEtaPhiM(GenJet_pt[j], GenJet_eta[j], GenJet_phi[j], GenJet_mass[j]);
      b_GenPart_dr_closest_j_vec[cRecIdx] = drRec[j].second;
      b_GenPart_x_closest_j_vec[cRecIdx] = GenPart_pt[cRecIdx] / Jet_pt[j];
      b_GenPart_dr_highest_j_vec[hRecIdx] = hgen_tlv.DeltaR(hjet_tlv);
      b_GenPart_x_highest_j_vec[hRecIdx] = xRec[j].second;
    }
    for (unsigned int j=0; j<nGenJet; ++j) {
      int cGenIdx = drGen[j].first; int hGenIdx = xGen[j].first;
      TLorentzVector hgen_tlv; TLorentzVector hjet_tlv;
      hgen_tlv.SetPtEtaPhiM(GenPart_pt[hGenIdx], GenPart_eta[hGenIdx], GenPart_phi[hGenIdx], GenPart_mass[hGenIdx]);
      hjet_tlv.SetPtEtaPhiM(GenJet_pt[j], GenJet_eta[j], GenJet_phi[j], GenJet_mass[j]);
      b_GenPart_dr_closest_gj_vec[cGenIdx] = drGen[j].second;
      b_GenPart_x_closest_gj_vec[cGenIdx] = GenPart_pt[cGenIdx] / Jet_pt[j];
      b_GenPart_dr_highest_gj_vec[hGenIdx] = hgen_tlv.DeltaR(hjet_tlv);
      b_GenPart_x_highest_gj_vec[hGenIdx] = xGen[j].second;
    }
    // Give flag ( == index) for the most closest(highest) KS-jet pair per event
    std::sort(drRec.begin(), drRec.end(), [] (std::pair<unsigned int, float> a, std::pair<unsigned int, float> b) { return (a.second < b.second); } );
    std::sort(drGen.begin(), drGen.end(), [] (std::pair<unsigned int, float> a, std::pair<unsigned int, float> b) { return (a.second < b.second); } );
    std::sort(xRec.begin(), xRec.end(), [] (std::pair<unsigned int, float> a, std::pair<unsigned int, float> b) { return (a.second > b.second); } );
    std::sort(xGen.begin(), xGen.end(), [] (std::pair<unsigned int, float> a, std::pair<unsigned int, float> b) { return (a.second > b.second); } );
    int cRecIdx = drRec[0].first; int cGenIdx = drGen[0].first; int hRecIdx = xRec[0].first; int hGenIdx = xGen[0].first;   
    std::replace_if(b_GenPart_isClosestPair_j_vec.begin(), b_GenPart_isClosestPair_j_vec.end(), [&] (int a) { return (a != cRecIdx); }, -1);
    std::replace_if(b_GenPart_isClosestPair_gj_vec.begin(), b_GenPart_isClosestPair_gj_vec.end(), [&] (int a) { return (a != cGenIdx); }, -1);
    std::replace_if(b_GenPart_isHighestPair_j_vec.begin(), b_GenPart_isHighestPair_j_vec.end(), [&] (int a) { return (a != hRecIdx); }, -1);
    std::replace_if(b_GenPart_isHighestPair_gj_vec.begin(), b_GenPart_isHighestPair_gj_vec.end(), [&] (int a) { return (a != hGenIdx); }, -1);
  }
}

void vtsAnalyser::RecAnalysis() {
  std::vector<std::pair<unsigned int, float>> drRec; std::vector<std::pair<unsigned int, float>> xRec;
  std::vector<std::pair<unsigned int, float>> drGen;  std::vector<std::pair<unsigned int, float>> xGen;
  int nRecKS = 0;
  b_hadTruth_dr_closest_j_vec.resize(nhad);
  b_hadTruth_x_closest_j_vec.resize(nhad);
  b_hadTruth_dr_closest_gj_vec.resize(nhad);
  b_hadTruth_x_closest_gj_vec.resize(nhad);
  b_hadTruth_dr_highest_j_vec.resize(nhad);
  b_hadTruth_x_highest_j_vec.resize(nhad);
  b_hadTruth_dr_highest_gj_vec.resize(nhad);
  b_hadTruth_x_highest_gj_vec.resize(nhad);
  for (unsigned int i=0; i<nhad; ++i) {
    if (had_pdgId[i] != 310) continue;
    ++nRecKS;
    b_hadTruth_isHadFromTop_vec.push_back(hadTruth_isHadFromTop[i]);
    b_hadTruth_isHadFromW_vec.push_back(hadTruth_isHadFromW[i]);
    b_hadTruth_isHadFromS_vec.push_back(hadTruth_isHadFromS[i]);
    b_hadTruth_isHadFromC_vec.push_back(hadTruth_isHadFromC[i]);
    b_hadTruth_isHadFromB_vec.push_back(hadTruth_isHadFromB[i]);
    b_hadTruth_isFrom_vec.push_back(hadTruth_isHadFromTsb[i]);
    b_hadTruth_nMatched_vec.push_back(hadTruth_nMatched[i]);
    b_hadTruth_d_vec.push_back(GetD(had_pt[i], had_eta[i], had_phi[i], had_mass[i], had_x[i], had_y[i], had_z[i]));
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

    TLorentzVector had_tlv;
    had_tlv.SetPtEtaPhiM(had_pt[i], had_eta[i], had_phi[i], had_mass[i]);
    for (unsigned int j=0; j<nJet; ++j) { // Loop for all of recoJet
      TLorentzVector jet_tlv;
      jet_tlv.SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
      auto dr = jet_tlv.DeltaR(had_tlv);
      auto x = had_pt[i]/Jet_pt[j];
      if (nRecKS==1) {
        drRec.push_back({i, dr});
        xRec.push_back({i, x});
      } else {
        if (hadTruth_isHadFromTsb[i] == m_qjMapForMC[j]) {
          if (dr < drRec[j].second ) { // Find the closest reco KS for each Jets
            drRec[j] = {i, dr};
          } /*else if (dr == drRec[j].second) { // not completed
          }*/
          if (x > xRec[j].second ) { // Find the highest x reco KS for each Jets
            xRec[j] = {i, x};
          } /*else if (x == xRec[j].second) { // not completed
          }*/
        }
      }
    }
    for (unsigned int j=0; j<nGenJet; ++j) { // Loop for all of recoJet
      TLorentzVector jet_tlv;
      jet_tlv.SetPtEtaPhiM(GenJet_pt[j], GenJet_eta[j], GenJet_phi[j], GenJet_mass[j]);
      auto dr = jet_tlv.DeltaR(had_tlv);
      auto x = had_pt[i]/GenJet_pt[j];
      if (nRecKS==1) {
        drGen.push_back({i, dr});
        xGen.push_back({i, x});
      } else {
        if (hadTruth_isHadFromTsb[i] == m_qgjMapForMC[j]) {
          if (dr < drGen[j].second ) { // Find the closest reco KS for each GenJets
            drGen[j] = {i, dr};
          } /*else if (dr == drGen[j].second) { // not completed
          }*/
          if (x > xGen[j].second ) { // Find the highest x reco KS for each GenJets
            xGen[j] = {i, x};
          } /*else if (x == xGen[j].second) { // not completed
          }*/
        }
      }
    }
    b_hadTruth_isClosestPair_j_vec.push_back(i); b_hadTruth_isClosestPair_gj_vec.push_back(i);
    b_hadTruth_isHighestPair_j_vec.push_back(i); b_hadTruth_isHighestPair_gj_vec.push_back(i);
  }
  if (nRecKS != 0) {
    //  Save x and dr for the most closest(highest) KS for each jet
    for (unsigned int j=0; j<nJet; ++j) {
      int cRecIdx = drRec[j].first; int hRecIdx = xRec[j].first;
      TLorentzVector hhad_tlv; TLorentzVector hjet_tlv;
      hhad_tlv.SetPtEtaPhiM(had_pt[hRecIdx], had_eta[hRecIdx], had_phi[hRecIdx], had_mass[hRecIdx]);
      hjet_tlv.SetPtEtaPhiM(GenJet_pt[j], GenJet_eta[j], GenJet_phi[j], GenJet_mass[j]);
      b_hadTruth_dr_closest_j_vec[cRecIdx] = drRec[j].second;
      b_hadTruth_x_closest_j_vec[cRecIdx] = had_pt[cRecIdx] / Jet_pt[j];
      b_hadTruth_dr_highest_j_vec[hRecIdx] = hhad_tlv.DeltaR(hjet_tlv);
      b_hadTruth_x_highest_j_vec[hRecIdx] = xRec[j].second;
    }
    for (unsigned int j=0; j<nGenJet; ++j) {
      int cGenIdx = drGen[j].first; int hGenIdx = xGen[j].first;
      TLorentzVector hhad_tlv; TLorentzVector hjet_tlv;
      hhad_tlv.SetPtEtaPhiM(had_pt[hGenIdx], had_eta[hGenIdx], had_phi[hGenIdx], had_mass[hGenIdx]);
      hjet_tlv.SetPtEtaPhiM(GenJet_pt[j], GenJet_eta[j], GenJet_phi[j], GenJet_mass[j]);
      b_hadTruth_dr_closest_gj_vec[cGenIdx] = drGen[j].second;
      b_hadTruth_x_closest_gj_vec[cGenIdx] = had_pt[cGenIdx] / Jet_pt[j];
      b_hadTruth_dr_highest_gj_vec[hGenIdx] = hhad_tlv.DeltaR(hjet_tlv);
      b_hadTruth_x_highest_gj_vec[hGenIdx] = xGen[j].second;
    }
    // Give flag ( == index) for the most closest(highest) KS-jet pair per event
    std::sort(drRec.begin(), drRec.end(), [] (std::pair<unsigned int, float> a, std::pair<unsigned int, float> b) { return (a.second < b.second); } );
    std::sort(drGen.begin(), drGen.end(), [] (std::pair<unsigned int, float> a, std::pair<unsigned int, float> b) { return (a.second < b.second); } );
    std::sort(xRec.begin(), xRec.end(), [] (std::pair<unsigned int, float> a, std::pair<unsigned int, float> b) { return (a.second > b.second); } );
    std::sort(xGen.begin(), xGen.end(), [] (std::pair<unsigned int, float> a, std::pair<unsigned int, float> b) { return (a.second > b.second); } );
    int cRecIdx = drRec[0].first; int cGenIdx = drGen[0].first; int hRecIdx = xRec[0].first; int hGenIdx = xGen[0].first; 
    std::replace_if(b_hadTruth_isClosestPair_j_vec.begin(), b_hadTruth_isClosestPair_j_vec.end(), [&] (int a) { return (a != cRecIdx); }, -1);
    std::replace_if(b_hadTruth_isClosestPair_gj_vec.begin(), b_hadTruth_isClosestPair_gj_vec.end(), [&] (int a) { return (a != cGenIdx); }, -1);
    std::replace_if(b_hadTruth_isHighestPair_j_vec.begin(), b_hadTruth_isHighestPair_j_vec.end(), [&] (int a) { return (a != hRecIdx); }, -1);
    std::replace_if(b_hadTruth_isHighestPair_gj_vec.begin(), b_hadTruth_isHighestPair_gj_vec.end(), [&] (int a) { return (a != hGenIdx); }, -1);
  }
}

void vtsAnalyser::CollectVar() {
  b_MET_pt = MET_pt;
  b_MET_phi = MET_phi;
  b_MET_sumEt = MET_sumEt;
}

