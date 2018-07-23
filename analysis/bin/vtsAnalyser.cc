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
      CollectVar();
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

  BranchI(channel); BranchI(njet); BranchF(met); BranchI(step);
  BranchI(hadTruth_nMatched); BranchI(hadTruth_nTrueDau); 
  BranchO(hadTruth_isHadFromTop); BranchI(hadTruth_isHadFromTsb); BranchO(hadTruth_isHadFromW); BranchO(hadTruth_isHadFromS); BranchO(hadTruth_isHadFromC); BranchO(hadTruth_isHadFromB);

  BranchTLV(had_tlv);
  BranchI(had_isFrom); BranchO(had_isHadJetMatched);
  BranchF(had_d); BranchF(had_x); BranchF(had_dr);
  BranchF(had_pt); BranchF(had_eta); BranchF(had_phi); BranchF(had_mass);
  BranchF(had_lxy); BranchF(had_lxySig); BranchF(had_angleXY); BranchF(had_angleXYZ); BranchF(had_chi2); BranchF(had_dca);
  BranchF(had_l3D); BranchF(had_l3DSig); BranchF(had_legDR); BranchI(had_pdgId);
  BranchF(had_dau1_chi2); BranchF(had_dau1_ipsigXY); BranchF(had_dau1_ipsigZ); BranchF(had_dau1_pt);
  BranchF(had_dau2_chi2); BranchF(had_dau2_ipsigXY); BranchF(had_dau2_ipsigZ); BranchF(had_dau2_pt);
  BranchF(Jet_btagCSVV2); BranchF(Jet_btagDeepB); BranchF(Jet_btagDeepC); BranchF(Jet_btagCMVA);
  BranchF(Jet_area); BranchF(Jet_pt); BranchI(Jet_nConstituents); BranchI(Jet_nElectrons); BranchI(Jet_nMuons);

  // For Test()
  BranchVI(GenPart_isGenFrom_vec); BranchVO(GenPart_isGenFromTop_vec); BranchVO(GenPart_isGenFromW_vec); BranchVO(GenPart_isFromKstar_vec);
  BranchVF(GenPart_d_vec);
  BranchVF(GenPart_pt_vec); BranchVF(GenPart_eta_vec); BranchVF(GenPart_phi_vec); BranchVF(GenPart_mass_vec);

  BranchVI(GenPart_isGenFrom_jetMat_vec); BranchVO(GenPart_isGenFromTop_jetMat_vec); BranchVO(GenPart_isGenFromW_jetMat_vec); BranchVO(GenPart_isFromKstar_jetMat_vec);
  BranchVF(GenPart_d_jetMat_vec); BranchVF(GenPart_x_jetMat_vec); BranchVF(GenPart_dr_jetMat_vec);
  BranchVF(GenPart_pt_jetMat_vec); BranchVF(GenPart_eta_jetMat_vec); BranchVF(GenPart_phi_jetMat_vec); BranchVF(GenPart_mass_jetMat_vec);

  BranchVI(GenPart_isGenFrom_gjetMat_vec); BranchVO(GenPart_isGenFromTop_gjetMat_vec); BranchVO(GenPart_isGenFromW_gjetMat_vec); BranchVO(GenPart_isFromKstar_gjetMat_vec);
  BranchVF(GenPart_d_gjetMat_vec); BranchVF(GenPart_x_gjetMat_vec); BranchVF(GenPart_dr_gjetMat_vec);
  BranchVF(GenPart_pt_gjetMat_vec); BranchVF(GenPart_eta_gjetMat_vec); BranchVF(GenPart_phi_gjetMat_vec); BranchVF(GenPart_mass_gjetMat_vec);

  BranchVI(GenPart_isFrom_j_vec); BranchVI(GenPart_nMatched_j_vec); BranchVO(GenPart_isInJet_j_vec); BranchVO(GenPart_isCorrectMat_j_vec); BranchVF(GenPart_d_j_vec); BranchVF(GenPart_x_j_vec); BranchVF(GenPart_dr_j_vec); BranchVF(GenPart_pt_j_vec);
  BranchVI(GenPart_isFrom_gj_vec); BranchVI(GenPart_nMatched_gj_vec); BranchVO(GenPart_isInJet_gj_vec); BranchVO(GenPart_isCorrectMat_gj_vec); BranchVF(GenPart_d_gj_vec); BranchVF(GenPart_x_gj_vec); BranchVF(GenPart_dr_gj_vec); BranchVF(GenPart_pt_gj_vec);

  BranchVI(hadTruth_nMatched_vec); BranchVI(hadTruth_isFrom_vec); BranchVO(hadTruth_isHadFromTop_vec); BranchVO(hadTruth_isHadFromW_vec); BranchVO(hadTruth_isHadFromS_vec); BranchVO(hadTruth_isHadFromC_vec); BranchVO(hadTruth_isHadFromB_vec);
  BranchVF(hadTruth_d_vec);
  BranchVF(hadTruth_pt_vec); BranchVF(hadTruth_eta_vec); BranchVF(hadTruth_phi_vec); BranchVF(hadTruth_mass_vec);
  BranchVF(hadTruth_lxy_vec); BranchVF(hadTruth_lxySig_vec); BranchVF(hadTruth_angleXY_vec); BranchVF(hadTruth_angleXYZ_vec); BranchVF(hadTruth_chi2_vec); BranchVF(hadTruth_dca_vec);
  BranchVF(hadTruth_l3D_vec); BranchVF(hadTruth_l3DSig_vec); BranchVF(hadTruth_legDR_vec); BranchVI(hadTruth_pdgId_vec);
  BranchVF(hadTruth_dau1_chi2_vec); BranchVF(hadTruth_dau1_ipsigXY_vec); BranchVF(hadTruth_dau1_ipsigZ_vec); BranchVF(hadTruth_dau1_pt_vec);
  BranchVF(hadTruth_dau2_chi2_vec); BranchVF(hadTruth_dau2_ipsigXY_vec); BranchVF(hadTruth_dau2_ipsigZ_vec); BranchVF(hadTruth_dau2_pt_vec);

  BranchVI(hadTruth_nMatched_jetMat_vec); BranchVI(hadTruth_isFrom_jetMat_vec); BranchVO(hadTruth_isHadFromTop_jetMat_vec); BranchVO(hadTruth_isHadFromW_jetMat_vec); BranchVO(hadTruth_isHadFromS_jetMat_vec); BranchVO(hadTruth_isHadFromC_jetMat_vec); BranchVO(hadTruth_isHadFromB_jetMat_vec); BranchVI(hadTruth_qjMapForMC_jetMat_vec);
  BranchVF(hadTruth_d_jetMat_vec); BranchVF(hadTruth_x_jetMat_vec); BranchVF(hadTruth_dr_jetMat_vec);
  BranchVF(hadTruth_pt_jetMat_vec); BranchVF(hadTruth_eta_jetMat_vec); BranchVF(hadTruth_phi_jetMat_vec); BranchVF(hadTruth_mass_jetMat_vec);
  BranchVF(hadTruth_lxy_jetMat_vec); BranchVF(hadTruth_lxySig_jetMat_vec); BranchVF(hadTruth_angleXY_jetMat_vec); BranchVF(hadTruth_angleXYZ_jetMat_vec); BranchVF(hadTruth_chi2_jetMat_vec); BranchVF(hadTruth_dca_jetMat_vec);
  BranchVF(hadTruth_l3D_jetMat_vec); BranchVF(hadTruth_l3DSig_jetMat_vec); BranchVF(hadTruth_legDR_jetMat_vec); BranchVI(hadTruth_pdgId_jetMat_vec);
  BranchVF(hadTruth_dau1_chi2_jetMat_vec); BranchVF(hadTruth_dau1_ipsigXY_jetMat_vec); BranchVF(hadTruth_dau1_ipsigZ_jetMat_vec); BranchVF(hadTruth_dau1_pt_jetMat_vec);
  BranchVF(hadTruth_dau2_chi2_jetMat_vec); BranchVF(hadTruth_dau2_ipsigXY_jetMat_vec); BranchVF(hadTruth_dau2_ipsigZ_jetMat_vec); BranchVF(hadTruth_dau2_pt_jetMat_vec);

  BranchVI(hadTruth_nMatched_gjetMat_vec); BranchVI(hadTruth_isFrom_gjetMat_vec); BranchVO(hadTruth_isHadFromTop_gjetMat_vec); BranchVO(hadTruth_isHadFromW_gjetMat_vec); BranchVO(hadTruth_isHadFromS_gjetMat_vec); BranchVO(hadTruth_isHadFromC_gjetMat_vec); BranchVO(hadTruth_isHadFromB_gjetMat_vec); BranchVI(hadTruth_qjMapForMC_gjetMat_vec);
  BranchVF(hadTruth_d_gjetMat_vec); BranchVF(hadTruth_x_gjetMat_vec); BranchVF(hadTruth_dr_gjetMat_vec);
  BranchVF(hadTruth_pt_gjetMat_vec); BranchVF(hadTruth_eta_gjetMat_vec); BranchVF(hadTruth_phi_gjetMat_vec); BranchVF(hadTruth_mass_gjetMat_vec);
  BranchVF(hadTruth_lxy_gjetMat_vec); BranchVF(hadTruth_lxySig_gjetMat_vec); BranchVF(hadTruth_angleXY_gjetMat_vec); BranchVF(hadTruth_angleXYZ_gjetMat_vec); BranchVF(hadTruth_chi2_gjetMat_vec); BranchVF(hadTruth_dca_gjetMat_vec);
  BranchVF(hadTruth_l3D_gjetMat_vec); BranchVF(hadTruth_l3DSig_gjetMat_vec); BranchVF(hadTruth_legDR_gjetMat_vec); BranchVI(hadTruth_pdgId_gjetMat_vec);
  BranchVF(hadTruth_dau1_chi2_gjetMat_vec); BranchVF(hadTruth_dau1_ipsigXY_gjetMat_vec); BranchVF(hadTruth_dau1_ipsigZ_gjetMat_vec); BranchVF(hadTruth_dau1_pt_gjetMat_vec);
  BranchVF(hadTruth_dau2_chi2_gjetMat_vec); BranchVF(hadTruth_dau2_ipsigXY_gjetMat_vec); BranchVF(hadTruth_dau2_ipsigZ_gjetMat_vec); BranchVF(hadTruth_dau2_pt_gjetMat_vec);

  BranchVI(hadTruth_isFrom_j_vec); BranchVI(hadTruth_nMatched_j_vec); BranchVO(hadTruth_isInJet_j_vec); BranchVO(hadTruth_isCorrectMat_j_vec); 
  BranchVF(hadTruth_d_j_vec); BranchVF(hadTruth_x_j_vec); BranchVF(hadTruth_dr_j_vec); 
  BranchVF(hadTruth_pt_j_vec);
  BranchVI(hadTruth_isFrom_gj_vec); BranchVI(hadTruth_nMatched_gj_vec); BranchVO(hadTruth_isInJet_gj_vec); BranchVO(hadTruth_isCorrectMat_gj_vec); 
  BranchVF(hadTruth_d_gj_vec); BranchVF(hadTruth_x_gj_vec); BranchVF(hadTruth_dr_gj_vec); 
  BranchVF(hadTruth_pt_gj_vec);

  BranchVI(hadTruth_idx_closest_j_vec); BranchVI(hadTruth_isFrom_closest_j_vec); BranchVO(hadTruth_isHadFromTop_closest_j_vec); BranchVI(hadTruth_nMatched_closest_j_vec); BranchVO(hadTruth_isInJet_closest_j_vec); BranchVO(hadTruth_isCorrectMat_closest_j_vec);
  BranchVF(hadTruth_d_closest_j_vec); BranchVF(hadTruth_x_closest_j_vec); BranchVF(hadTruth_dr_closest_j_vec);
  BranchVF(hadTruth_pt_closest_j_vec); BranchVF(hadTruth_eta_closest_j_vec); BranchVF(hadTruth_phi_closest_j_vec); BranchVF(hadTruth_mass_closest_j_vec);

  BranchVI(hadTruth_idx_highest_j_vec); BranchVI(hadTruth_isFrom_highest_j_vec); BranchVO(hadTruth_isHadFromTop_highest_j_vec); BranchVI(hadTruth_nMatched_highest_j_vec); BranchVO(hadTruth_isInJet_highest_j_vec); BranchVO(hadTruth_isCorrectMat_highest_j_vec);
  BranchVF(hadTruth_d_highest_j_vec); BranchVF(hadTruth_x_highest_j_vec); BranchVF(hadTruth_dr_highest_j_vec);
  BranchVF(hadTruth_pt_highest_j_vec); BranchVF(hadTruth_eta_highest_j_vec); BranchVF(hadTruth_phi_highest_j_vec); BranchVF(hadTruth_mass_highest_j_vec);

  BranchI(hadTruth_isFrom_closest_j); BranchI(hadTruth_isHadFromTop_closest_j); BranchI(hadTruth_nMatched_closest_j); BranchI(hadTruth_isInJet_closest_j); BranchI(hadTruth_isCorrectMat_closest_j);
  BranchF(hadTruth_d_closest_j); BranchF(hadTruth_x_closest_j); BranchF(hadTruth_dr_closest_j);
  BranchF(hadTruth_pt_closest_j); BranchF(hadTruth_eta_closest_j); BranchF(hadTruth_phi_closest_j); BranchF(hadTruth_mass_closest_j);

  BranchI(hadTruth_isFrom_highest_j); BranchI(hadTruth_isHadFromTop_highest_j);BranchI(hadTruth_nMatched_highest_j); BranchI(hadTruth_isInJet_highest_j); BranchI(hadTruth_isCorrectMat_highest_j);
  BranchF(hadTruth_d_highest_j); BranchF(hadTruth_x_highest_j); BranchF(hadTruth_dr_highest_j);
  BranchF(hadTruth_pt_highest_j); BranchF(hadTruth_eta_highest_j); BranchF(hadTruth_phi_highest_j); BranchF(hadTruth_mass_highest_j);

  // For CollectVar()
  BranchVF(lep_pt_vec); BranchVF(dilep_pt_vec); BranchVF(elec_pt_vec); BranchVF(mu_pt_vec);
  BranchVF(lep_eta_vec); BranchVF(dilep_eta_vec); BranchVF(elec_eta_vec); BranchVF(mu_eta_vec);
  BranchVF(lep_phi_vec); BranchVF(dilep_phi_vec); BranchVF(elec_phi_vec); BranchVF(mu_phi_vec);
  BranchVF(lep_mass_vec); BranchVF(dilep_mass_vec); BranchVF(elec_mass_vec); BranchVF(mu_mass_vec);
  BranchVF(lep_dxy_vec); BranchVF(elec_dxy_vec); BranchVF(mu_dxy_vec);
  BranchF(MET_pt); BranchF(MET_phi); BranchF(MET_sumEt);
}

void vtsAnalyser::ResetBranch() {
  Reset();

  b_hadTruth_nMatched = -1; b_hadTruth_nTrueDau = -1;
  b_hadTruth_isHadFromTsb = -1;
  b_hadTruth_isHadFromTop = false; b_hadTruth_isHadFromW = false; b_hadTruth_isHadFromS = false; b_hadTruth_isHadFromC = false; b_hadTruth_isHadFromB = false;

  // For Test()
  b_GenPart_isGenFrom_vec.clear(); b_GenPart_isGenFromTop_vec.clear(); b_GenPart_isGenFromW_vec.clear(); b_GenPart_isFromKstar_vec.clear();
  b_GenPart_d_vec.clear();
  b_GenPart_pt_vec.clear(); b_GenPart_eta_vec.clear(); b_GenPart_phi_vec.clear(); b_GenPart_mass_vec.clear();

  b_GenPart_isGenFrom_jetMat_vec.clear(); b_GenPart_isGenFromTop_jetMat_vec.clear(); b_GenPart_isGenFromW_jetMat_vec.clear(); b_GenPart_isFromKstar_jetMat_vec.clear();
  b_GenPart_d_jetMat_vec.clear(); b_GenPart_x_jetMat_vec.clear(); b_GenPart_dr_jetMat_vec.clear();
  b_GenPart_pt_jetMat_vec.clear(); b_GenPart_eta_jetMat_vec.clear(); b_GenPart_phi_jetMat_vec.clear(); b_GenPart_mass_jetMat_vec.clear();

  b_GenPart_isGenFrom_gjetMat_vec.clear(); b_GenPart_isGenFromTop_gjetMat_vec.clear(); b_GenPart_isGenFromW_gjetMat_vec.clear(); b_GenPart_isFromKstar_gjetMat_vec.clear();
  b_GenPart_d_gjetMat_vec.clear(); b_GenPart_x_gjetMat_vec.clear(); b_GenPart_dr_gjetMat_vec.clear();
  b_GenPart_pt_gjetMat_vec.clear(); b_GenPart_eta_gjetMat_vec.clear(); b_GenPart_phi_gjetMat_vec.clear(); b_GenPart_mass_gjetMat_vec.clear();

  b_GenPart_isFrom_j_vec.clear(); b_GenPart_nMatched_j_vec.clear(); b_GenPart_isInJet_j_vec.clear(); b_GenPart_isCorrectMat_j_vec.clear(); b_GenPart_d_j_vec.clear(); b_GenPart_x_j_vec.clear(); b_GenPart_dr_j_vec.clear(); b_GenPart_pt_j_vec.clear();
  b_GenPart_isFrom_gj_vec.clear(); b_GenPart_nMatched_gj_vec.clear(); b_GenPart_isInJet_gj_vec.clear(); b_GenPart_isCorrectMat_gj_vec.clear(); b_GenPart_d_gj_vec.clear(); b_GenPart_x_gj_vec.clear(); b_GenPart_dr_gj_vec.clear(); b_GenPart_pt_gj_vec.clear();

  b_hadTruth_nMatched_vec.clear(); b_hadTruth_isFrom_vec.clear(); b_hadTruth_isHadFromTop_vec.clear(); b_hadTruth_isHadFromW_vec.clear(); b_hadTruth_isHadFromS_vec.clear(); b_hadTruth_isHadFromC_vec.clear(); b_hadTruth_isHadFromB_vec.clear(); 
  b_hadTruth_d_vec.clear();
  b_hadTruth_pt_vec.clear(); b_hadTruth_eta_vec.clear(); b_hadTruth_phi_vec.clear(); b_hadTruth_mass_vec.clear();
  b_hadTruth_lxy_vec.clear(); b_hadTruth_lxySig_vec.clear(); b_hadTruth_angleXY_vec.clear(); b_hadTruth_angleXYZ_vec.clear(); b_hadTruth_chi2_vec.clear(); b_hadTruth_dca_vec.clear();
  b_hadTruth_l3D_vec.clear(); b_hadTruth_l3DSig_vec.clear(); b_hadTruth_legDR_vec.clear(); b_hadTruth_pdgId_vec.clear();
  b_hadTruth_dau1_chi2_vec.clear(); b_hadTruth_dau1_ipsigXY_vec.clear(); b_hadTruth_dau1_ipsigZ_vec.clear(); b_hadTruth_dau1_pt_vec.clear();
  b_hadTruth_dau2_chi2_vec.clear(); b_hadTruth_dau2_ipsigXY_vec.clear(); b_hadTruth_dau2_ipsigZ_vec.clear(); b_hadTruth_dau2_pt_vec.clear();

  b_hadTruth_nMatched_jetMat_vec.clear(); b_hadTruth_isFrom_jetMat_vec.clear(); b_hadTruth_isHadFromTop_jetMat_vec.clear(); b_hadTruth_isHadFromW_jetMat_vec.clear(); b_hadTruth_isHadFromS_jetMat_vec.clear(); b_hadTruth_isHadFromC_jetMat_vec.clear(); b_hadTruth_isHadFromB_jetMat_vec.clear(); b_hadTruth_qjMapForMC_jetMat_vec.clear();
  b_hadTruth_d_jetMat_vec.clear(); b_hadTruth_x_jetMat_vec.clear(); b_hadTruth_dr_jetMat_vec.clear();
  b_hadTruth_pt_jetMat_vec.clear(); b_hadTruth_eta_jetMat_vec.clear(); b_hadTruth_phi_jetMat_vec.clear(); b_hadTruth_mass_jetMat_vec.clear();
  b_hadTruth_lxy_jetMat_vec.clear(); b_hadTruth_lxySig_jetMat_vec.clear(); b_hadTruth_angleXY_jetMat_vec.clear(); b_hadTruth_angleXYZ_jetMat_vec.clear(); b_hadTruth_chi2_jetMat_vec.clear(); b_hadTruth_dca_jetMat_vec.clear();
  b_hadTruth_l3D_jetMat_vec.clear(); b_hadTruth_l3DSig_jetMat_vec.clear(); b_hadTruth_legDR_jetMat_vec.clear(); b_hadTruth_pdgId_jetMat_vec.clear();
  b_hadTruth_dau1_chi2_jetMat_vec.clear(); b_hadTruth_dau1_ipsigXY_jetMat_vec.clear(); b_hadTruth_dau1_ipsigZ_jetMat_vec.clear(); b_hadTruth_dau1_pt_jetMat_vec.clear();
  b_hadTruth_dau2_chi2_jetMat_vec.clear(); b_hadTruth_dau2_ipsigXY_jetMat_vec.clear(); b_hadTruth_dau2_ipsigZ_jetMat_vec.clear(); b_hadTruth_dau2_pt_jetMat_vec.clear();

  b_hadTruth_nMatched_gjetMat_vec.clear(); b_hadTruth_isFrom_gjetMat_vec.clear(); b_hadTruth_isHadFromTop_gjetMat_vec.clear(); b_hadTruth_isHadFromW_gjetMat_vec.clear(); b_hadTruth_isHadFromS_gjetMat_vec.clear(); b_hadTruth_isHadFromC_gjetMat_vec.clear(); b_hadTruth_isHadFromB_gjetMat_vec.clear(); b_hadTruth_qjMapForMC_gjetMat_vec.clear();
  b_hadTruth_d_gjetMat_vec.clear(); b_hadTruth_x_gjetMat_vec.clear(); b_hadTruth_dr_gjetMat_vec.clear();
  b_hadTruth_pt_gjetMat_vec.clear(); b_hadTruth_eta_gjetMat_vec.clear(); b_hadTruth_phi_gjetMat_vec.clear(); b_hadTruth_mass_gjetMat_vec.clear();
  b_hadTruth_lxy_gjetMat_vec.clear(); b_hadTruth_lxySig_gjetMat_vec.clear(); b_hadTruth_angleXY_gjetMat_vec.clear(); b_hadTruth_angleXYZ_gjetMat_vec.clear(); b_hadTruth_chi2_gjetMat_vec.clear(); b_hadTruth_dca_gjetMat_vec.clear();
  b_hadTruth_l3D_gjetMat_vec.clear(); b_hadTruth_l3DSig_gjetMat_vec.clear(); b_hadTruth_legDR_gjetMat_vec.clear(); b_hadTruth_pdgId_gjetMat_vec.clear();
  b_hadTruth_dau1_chi2_gjetMat_vec.clear(); b_hadTruth_dau1_ipsigXY_gjetMat_vec.clear(); b_hadTruth_dau1_ipsigZ_gjetMat_vec.clear(); b_hadTruth_dau1_pt_gjetMat_vec.clear();
  b_hadTruth_dau2_chi2_gjetMat_vec.clear(); b_hadTruth_dau2_ipsigXY_gjetMat_vec.clear(); b_hadTruth_dau2_ipsigZ_gjetMat_vec.clear(); b_hadTruth_dau2_pt_gjetMat_vec.clear();

  b_hadTruth_isFrom_j_vec.clear(); b_hadTruth_nMatched_j_vec.clear(); b_hadTruth_isInJet_j_vec.clear(); b_hadTruth_isCorrectMat_j_vec.clear(); 
  b_hadTruth_d_j_vec.clear(); b_hadTruth_x_j_vec.clear(); b_hadTruth_dr_j_vec.clear(); 
  b_hadTruth_pt_j_vec.clear();
  b_hadTruth_isFrom_gj_vec.clear(); b_hadTruth_nMatched_gj_vec.clear(); b_hadTruth_isInJet_gj_vec.clear(); b_hadTruth_isCorrectMat_gj_vec.clear(); 
  b_hadTruth_d_gj_vec.clear(); b_hadTruth_x_gj_vec.clear(); b_hadTruth_dr_gj_vec.clear(); 
  b_hadTruth_pt_gj_vec.clear();

  b_hadTruth_idx_closest_j_vec.clear();  b_hadTruth_isFrom_closest_j_vec.clear(); b_hadTruth_isHadFromTop_closest_j_vec.clear(); b_hadTruth_nMatched_closest_j_vec.clear(); b_hadTruth_isInJet_closest_j_vec.clear(); b_hadTruth_isCorrectMat_closest_j_vec.clear();
  b_hadTruth_d_closest_j_vec.clear(); b_hadTruth_x_closest_j_vec.clear(); b_hadTruth_dr_closest_j_vec.clear();
  b_hadTruth_pt_closest_j_vec.clear(); b_hadTruth_eta_closest_j_vec.clear(); b_hadTruth_phi_closest_j_vec.clear(); b_hadTruth_mass_closest_j_vec.clear();

  b_hadTruth_idx_highest_j_vec.clear(); b_hadTruth_isFrom_highest_j_vec.clear(); b_hadTruth_isHadFromTop_highest_j_vec.clear(); b_hadTruth_nMatched_highest_j_vec.clear(); b_hadTruth_isInJet_highest_j_vec.clear(); b_hadTruth_isCorrectMat_highest_j_vec.clear();
  b_hadTruth_d_highest_j_vec.clear(); b_hadTruth_x_highest_j_vec.clear(); b_hadTruth_dr_highest_j_vec.clear();
  b_hadTruth_pt_highest_j_vec.clear(); b_hadTruth_eta_highest_j_vec.clear(); b_hadTruth_phi_highest_j_vec.clear(); b_hadTruth_mass_highest_j_vec.clear();

  b_hadTruth_isFrom_closest_j = -99; b_hadTruth_isHadFromTop_closest_j = -1; b_hadTruth_nMatched_closest_j = -1; b_hadTruth_isInJet_closest_j = -1; b_hadTruth_isCorrectMat_closest_j = -1;
  b_hadTruth_d_closest_j = -1; b_hadTruth_x_closest_j = -1; b_hadTruth_dr_closest_j = -1;
  b_hadTruth_pt_closest_j = -1; b_hadTruth_eta_closest_j = -99; b_hadTruth_phi_closest_j = -99; b_hadTruth_mass_closest_j = -99;

  b_hadTruth_isFrom_highest_j = -99; b_hadTruth_isHadFromTop_highest_j = -1; b_hadTruth_nMatched_highest_j = -1; b_hadTruth_isInJet_highest_j = -1; b_hadTruth_isCorrectMat_highest_j = -1;
  b_hadTruth_d_highest_j = -1; b_hadTruth_x_highest_j = -1; b_hadTruth_dr_highest_j = -1;
  b_hadTruth_pt_highest_j = -1; b_hadTruth_eta_highest_j = -99; b_hadTruth_phi_highest_j = -99; b_hadTruth_mass_highest_j = -99;

  // For CollectVar()
  b_lep_pt_vec.clear(); b_dilep_pt_vec.clear(); b_elec_pt_vec.clear(); b_mu_pt_vec.clear();
  b_lep_eta_vec.clear(); b_dilep_eta_vec.clear(); b_elec_eta_vec.clear(); b_mu_eta_vec.clear();
  b_lep_phi_vec.clear(); b_dilep_phi_vec.clear(); b_elec_phi_vec.clear(); b_mu_phi_vec.clear();
  b_lep_mass_vec.clear(); b_dilep_mass_vec.clear(); b_elec_mass_vec.clear(); b_mu_mass_vec.clear();
  b_lep_dxy_vec.clear(); b_elec_dxy_vec.clear(); b_mu_dxy_vec.clear();
  b_MET_pt = -1; b_MET_phi = -99; b_MET_sumEt = -1;

  b_had_tlv.SetPtEtaPhiM(0,0,0,0);
  b_had_isFrom = -99;
  b_had_isHadJetMatched = false;
  b_had_d = -1; b_had_x = -1; b_had_dr = -1;
  b_had_pt = -1; b_had_eta = -99; b_had_phi = -99; b_had_mass = -99;
  b_had_lxy = -1; b_had_lxySig = -1; b_had_angleXY = -1; b_had_angleXYZ = -1; b_had_chi2 = -1; b_had_dca = -1;
  b_had_l3D = -1; b_had_l3DSig = -1; b_had_legDR = -1; b_had_pdgId = -99;
  b_had_dau1_chi2 = -1; b_had_dau1_ipsigXY = -1; b_had_dau1_ipsigZ = -1; b_had_dau1_pt = -1;
  b_had_dau2_chi2 = -1; b_had_dau2_ipsigXY = -1; b_had_dau2_ipsigZ = -1; b_had_dau2_pt = -1;

  b_Jet_btagCSVV2 = -99; b_Jet_btagCMVA = -99; b_Jet_btagDeepB = -99; b_Jet_btagDeepC = -99;
  b_Jet_area = -1; b_Jet_pt = -1; b_Jet_nConstituents = -1; b_Jet_nElectrons = -1; b_Jet_nMuons = -1;

  qjMapForMC_.clear(); qgjMapForMC_.clear(); tqMC_.clear(); wqMC_.clear(); genJet_.clear(); recoJet_.clear();
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
    if (std::abs(GenPart_status[i] - 25) >= 5 ) continue;
    if (abs(GenPart_pdgId[i]) == 3 || abs(GenPart_pdgId[i]) == 5 || abs(GenPart_pdgId[i]) == 4) {
      if ( (abs(GenPart_pdgId[GenPart_genPartIdxMother[i]]) == 6 && GenPart_status[GenPart_genPartIdxMother[i]] == 62 ) ) { 
        tqMC_.push_back(i);
//        cout << "TOP GENSTATUSFLAG(pdg, stat) : " << GenPart_statusFlags[GenPart_genPartIdxMother[i]] << " ( " << GenPart_pdgId[GenPart_genPartIdxMother[i]] << " , " << GenPart_status[GenPart_genPartIdxMother[i]] << ") , QUARK FLAG(pdg, stat) : " << GenPart_statusFlags[i] << " ( " << GenPart_pdgId[i] << " , " << GenPart_status[i] << " ) " <<  endl; 
      }
      if ( (abs(GenPart_pdgId[GenPart_genPartIdxMother[i]]) == 24 && (GenPart_status[GenPart_genPartIdxMother[i]] == 22 || GenPart_status[GenPart_genPartIdxMother[i]] == 52)) ) {
        wqMC_.push_back(i);
//        if (GenPart_status == 22) cout << "TOP GENSTATUSFLAG(pdg, stat) : " << GenPart_statusFlags[GenPart_genPartIdxMother[GenPart_genPartIdxMother[i]]] << " ( " << GenPart_pdgId[GenPart_genPartIdxMother[GenPart_genPartIdxMother[i]]] << " , " << GenPart_status[GenPart_genPartIdxMother[GenPart_genPartIdxMother[i]]] << " ) , W BOSON GENSTATUSFLAG(pdg, stat) : " << GenPart_statusFlags[GenPart_genPartIdxMother[i]] << " ( " << GenPart_pdgId[GenPart_genPartIdxMother[i]] << " , " << GenPart_status[GenPart_genPartIdxMother[i]] << " ) , QUARK FLAG(pdg, stat) : " << GenPart_statusFlags[i] << " ( " << GenPart_pdgId[i] << " , " << GenPart_status[i] << " ) " <<  endl;
//        if (GenPart_status == 52) cout << "TOP GENSTATUSFLAG(pdg, stat) : " << GenPart_statusFlags[GenPart_genPartIdxMother[GenPart_genPartIdxMother[GenPart_genPartIdxMother[i]]]] << " ( " << GenPart_pdgId[GenPart_genPartIdxMother[GenPart_genPartIdxMother[GenPart_genPartIdxMother[i]]]] << " , " << GenPart_status[GenPart_genPartIdxMother[GenPart_genPartIdxMother[GenPart_genPartIdxMother[i]]]] << " ) , W BOSON 1 GENSTATUSFLAG(pdg, stat) : " << GenPart_statusFlags[GenPart_genPartIdxMother[GenPart_genPartIdxMother[i]]] << " ( " << GenPart_pdgId[GenPart_genPartIdxMother[GenPart_genPartIdxMother[i]]] << " , " << GenPart_status[GenPart_genPartIdxMother[GenPart_genPartIdxMother[i]]] << " ) , W BOSON 2 GENSTATUSFLAG(pdg, stat) : " << GenPart_statusFlags[GenPart_genPartIdxMother[i]] << " ( " << GenPart_pdgId[GenPart_genPartIdxMother[i]] << " , " << GenPart_status[GenPart_genPartIdxMother[i]] << " ) , QUARK FLAG(pdg, stat) : " << GenPart_statusFlags[i] << " ( " << GenPart_pdgId[i] << " , " << GenPart_status[i] << " ) " <<  endl;
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
    if (dr > 0.5) qjMapForMC_.insert({j, -9});
    else {
      qjMapForMC_.insert({j, GenPart_pdgId[qIdx]});
      JetStat stat(j, dr, qjMapForMC_[j]);
      recoJet_.push_back(stat);
    }
  }

  if (recoJet_.size() > 1) {
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

  // Gen Particle & Gen Jet matching
  for (unsigned int j=0; j<nGenJet; ++j) {
    TLorentzVector jet_tlv;
    jet_tlv.SetPtEtaPhiM(GenJet_pt[j], GenJet_eta[j], GenJet_phi[j], GenJet_mass[j]);
    auto qIdx = tq1;
    auto dr = tq1_tlv.DeltaR(jet_tlv);
    auto dr_ = tq2_tlv.DeltaR(jet_tlv);
    if ( dr_ < dr) {
      dr = dr_;
      qIdx = tq2;
    }
    if (dr > 0.5) qgjMapForMC_.insert({j, -9});
    else {
      if ( GenPart_pdgId[qIdx] == GenJet_partonFlavour[j] ) qgjMapForMC_.insert({j, GenPart_pdgId[qIdx]});
      else qgjMapForMC_.insert({j,-9});
      JetStat stat(j, dr, qgjMapForMC_[j]);
      genJet_.push_back(stat);
//      cout << "[ " << j << " / " << nGenJet << " ] chk ==> qgjMapForMC_ : " << qgjMapForMC_[j] << " , partonFlavour : " << GenJet_partonFlavour[j] << endl;
    }
  }

  if (genJet_.size() > 1) {
    sort(genJet_.begin(), genJet_.end(), [](struct JetStat a, struct JetStat b) { return (a.dr < b.dr); } );

    bool idChanged = false;
    for (unsigned int i=1;i<genJet_.size(); ++i) {
      if (idChanged) qgjMapForMC_[genJet_[i].idx] = -9;
      else {
          if (genJet_[0].matchedQuark != genJet_[i].matchedQuark) idChanged = true;
          else qgjMapForMC_[genJet_[i].idx] = -9;
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
      if (jet_tlv.DeltaR(had_tlv) < 0.5 && (had_pt[k]/Jet_pt[j]) > 0.15) {
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
    b_had_isFrom = JetCollection[0][0].label;  // -99 : event that can't pass till step4(jet selection) or there is no matching between had and jet, -9 : there is t->qW in the event,but not matched to recoJet, 0 : there is no t->qW in the event (if no t->s and no matching between had-jet, then the event would be -99), +-3 : hadron is from t->sW, +-5 : hadron is from t->bW
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
//  cout << "[ " << count << " ] => gen pdgId(status) : " << GenPart_pdgId[idx] << " ( " << GenPart_status[idx] << " ) | mom pdgId(status) " <<  GenPart_pdgId[GenPart_genPartIdxMother[idx]] << " ( " << GenPart_status[GenPart_genPartIdxMother[idx]] << " ) " << endl;
  if ((abs(GenPart_pdgId[idx]) == 3 || abs(GenPart_pdgId[idx]) == 4 || abs(GenPart_pdgId[idx]) == 5) && GenPart_status[idx] == 23) {
    if (abs(GenPart_pdgId[GenPart_genPartIdxMother[idx]]) == 6 && GenPart_status[GenPart_genPartIdxMother[idx]] == 62) {
      isFrom = GenPart_pdgId[idx];
      isFromTop = true;
      isFromW = false;
//      cout << "gen KS is from top !" << endl;
      return true;
    } else if (abs(GenPart_pdgId[GenPart_genPartIdxMother[idx]]) == 24 && (GenPart_status[GenPart_genPartIdxMother[idx]] == 22 || GenPart_status[GenPart_genPartIdxMother[idx]] == 52)) {
      isFrom = GenPart_pdgId[idx];
      isFromW = true;
      int midx = 0;
      if (GenPart_status[GenPart_genPartIdxMother[idx]] == 22) midx = GenPart_genPartIdxMother[GenPart_genPartIdxMother[idx]];
      else if (GenPart_status[GenPart_genPartIdxMother[idx]] == 52) midx = GenPart_genPartIdxMother[GenPart_genPartIdxMother[GenPart_genPartIdxMother[idx]]];
      if (abs(GenPart_pdgId[midx]) == 6 && GenPart_status[midx] == 62) {
        isFromTop = true;
//        cout << "gen KS is from t->W !" << endl;
      } else {
        isFromTop = false;
//        cout << "gen KS is from W boson !" << endl;
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

//Not yet completed
int vtsAnalyser::Test() {
  for (unsigned int i=0; i<nGenPart; ++i) {
    if (abs(GenPart_pdgId[i]) != 310) continue;
    int count = 0;
    int isFrom = 0;
    bool isFromTop = false;
    bool isFromW = false;
    bool isFromKstar = false;
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
//    cout << "b_GenPart_isGenFrom_vec size : " << b_GenPart_isGenFrom_vec.size() << " , isFrom : " << isFrom << " , isFromTop : " << isFromTop << " , isFromW : " << isFromW << endl;
    for (unsigned int j=0; j<nJet; ++j) { // Loop for all of recoJet
      TLorentzVector jet_tlv;
      jet_tlv.SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
      // gen KS inside Jet 
      if ((jet_tlv.DeltaR(gen_tlv) < 0.5) && (GenPart_pt[i]/Jet_pt[j] > 0.15)) {
        if (isFrom == qjMapForMC_[j]) { // For jet with the same origin as gen KS
          b_GenPart_isFromKstar_jetMat_vec.push_back(isFromKstar);
          b_GenPart_isGenFrom_jetMat_vec.push_back(isFrom);
          b_GenPart_isGenFromTop_jetMat_vec.push_back(isFromTop);
          b_GenPart_isGenFromW_jetMat_vec.push_back(isFromW);
          b_GenPart_d_jetMat_vec.push_back(GetD(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i], gen_tlv.X(), gen_tlv.Y(), gen_tlv.Z()));
          b_GenPart_x_jetMat_vec.push_back(GenPart_pt[i]/Jet_pt[j]);
          b_GenPart_dr_jetMat_vec.push_back(jet_tlv.DeltaR(gen_tlv));
          b_GenPart_pt_jetMat_vec.push_back(GenPart_pt[i]);
          b_GenPart_eta_jetMat_vec.push_back(GenPart_eta[i]);
          b_GenPart_phi_jetMat_vec.push_back(GenPart_phi[i]);
          b_GenPart_mass_jetMat_vec.push_back(GenPart_mass[i]);

          b_GenPart_isFrom_j_vec.push_back(hadTruth_isHadFromTsb[i]);
          b_GenPart_nMatched_j_vec.push_back(hadTruth_nMatched[i]);
          b_GenPart_isInJet_j_vec.push_back(true);
          b_GenPart_isCorrectMat_j_vec.push_back(true);
          b_GenPart_d_j_vec.push_back(GetD(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i], gen_tlv.X(), gen_tlv.Y(), gen_tlv.Z()));
          b_GenPart_x_j_vec.push_back(GenPart_pt[i]/Jet_pt[j]);
          b_GenPart_dr_j_vec.push_back(jet_tlv.DeltaR(gen_tlv));
          b_GenPart_pt_j_vec.push_back(GenPart_pt[i]);
        } else{
          b_GenPart_isFrom_j_vec.push_back(hadTruth_isHadFromTsb[i]);
          b_GenPart_nMatched_j_vec.push_back(hadTruth_nMatched[i]);
          b_GenPart_isInJet_j_vec.push_back(true);
          b_GenPart_isCorrectMat_j_vec.push_back(false);
          b_GenPart_d_j_vec.push_back(GetD(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i], gen_tlv.X(), gen_tlv.Y(), gen_tlv.Z()));
          b_GenPart_x_j_vec.push_back(GenPart_pt[i]/Jet_pt[j]);
          b_GenPart_dr_j_vec.push_back(jet_tlv.DeltaR(gen_tlv));
          b_GenPart_pt_j_vec.push_back(GenPart_pt[i]);
        }
      } else {
        b_GenPart_isFrom_j_vec.push_back(hadTruth_isHadFromTsb[i]);
        b_GenPart_nMatched_j_vec.push_back(hadTruth_nMatched[i]);
        b_GenPart_isInJet_j_vec.push_back(false);
        b_GenPart_isCorrectMat_j_vec.push_back(false);
        b_GenPart_d_jetMat_vec.push_back(GetD(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i], gen_tlv.X(), gen_tlv.Y(), gen_tlv.Z()));
        b_GenPart_x_j_vec.push_back(GenPart_pt[i]/Jet_pt[j]);
        b_GenPart_dr_j_vec.push_back(jet_tlv.DeltaR(gen_tlv));
        b_GenPart_pt_j_vec.push_back(GenPart_pt[i]);
      }
    }
    for (unsigned int j=0; j<nGenJet; ++j) { // Loop for all of genJet
      TLorentzVector jet_tlv;
      jet_tlv.SetPtEtaPhiM(GenJet_pt[j], GenJet_eta[j], GenJet_phi[j], GenJet_mass[j]);
      // gen KS inside GenJet  
      if ((jet_tlv.DeltaR(gen_tlv) < 0.5) && (GenPart_pt[i]/Jet_pt[j] > 0.15)) {
        if (isFrom == qgjMapForMC_[j]) {//GenJet_partonFlavour_[j]) { // For jet with the same origin as gen KS
          b_GenPart_isFromKstar_gjetMat_vec.push_back(isFromKstar);
          b_GenPart_isGenFrom_gjetMat_vec.push_back(isFrom);
          b_GenPart_isGenFromTop_gjetMat_vec.push_back(isFromTop);
          b_GenPart_isGenFromW_gjetMat_vec.push_back(isFromW);
          b_GenPart_d_gjetMat_vec.push_back(GetD(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i], gen_tlv.X(), gen_tlv.Y(), gen_tlv.Z()));
          b_GenPart_x_gjetMat_vec.push_back(GenPart_pt[i]/GenJet_pt[j]);
          b_GenPart_dr_gjetMat_vec.push_back(jet_tlv.DeltaR(gen_tlv));
          b_GenPart_pt_gjetMat_vec.push_back(GenPart_pt[i]);
          b_GenPart_eta_gjetMat_vec.push_back(GenPart_eta[i]);
          b_GenPart_phi_gjetMat_vec.push_back(GenPart_phi[i]);
          b_GenPart_mass_gjetMat_vec.push_back(GenPart_mass[i]);
       
          b_GenPart_isFrom_gj_vec.push_back(hadTruth_isHadFromTsb[i]);
          b_GenPart_nMatched_gj_vec.push_back(hadTruth_nMatched[i]);
          b_GenPart_isInJet_gj_vec.push_back(true);
          b_GenPart_isCorrectMat_gj_vec.push_back(true);
          b_GenPart_d_gj_vec.push_back(GetD(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i], gen_tlv.X(), gen_tlv.Y(), gen_tlv.Z()));
          b_GenPart_x_gj_vec.push_back(GenPart_pt[i]/Jet_pt[j]);
          b_GenPart_dr_gj_vec.push_back(jet_tlv.DeltaR(gen_tlv));
          b_GenPart_pt_gj_vec.push_back(GenPart_pt[i]);
        } else{
          b_GenPart_isFrom_gj_vec.push_back(hadTruth_isHadFromTsb[i]);
          b_GenPart_nMatched_gj_vec.push_back(hadTruth_nMatched[i]);
          b_GenPart_isInJet_gj_vec.push_back(true);
          b_GenPart_isCorrectMat_gj_vec.push_back(false);
          b_GenPart_d_gj_vec.push_back(GetD(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i], gen_tlv.X(), gen_tlv.Y(), gen_tlv.Z()));
          b_GenPart_x_gj_vec.push_back(GenPart_pt[i]/Jet_pt[j]);
          b_GenPart_dr_gj_vec.push_back(jet_tlv.DeltaR(gen_tlv));
          b_GenPart_pt_gj_vec.push_back(GenPart_pt[i]);
        }
      } else {
        b_GenPart_isFrom_gj_vec.push_back(hadTruth_isHadFromTsb[i]);
        b_GenPart_nMatched_gj_vec.push_back(hadTruth_nMatched[i]);
        b_GenPart_isInJet_gj_vec.push_back(false);
        b_GenPart_isCorrectMat_gj_vec.push_back(false);
        b_GenPart_d_gj_vec.push_back(GetD(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i], gen_tlv.X(), gen_tlv.Y(), gen_tlv.Z()));
        b_GenPart_x_gj_vec.push_back(GenPart_pt[i]/Jet_pt[j]);
        b_GenPart_dr_gj_vec.push_back(jet_tlv.DeltaR(gen_tlv));
        b_GenPart_pt_gj_vec.push_back(GenPart_pt[i]);
      }     
    }  
  }

  std::vector<std::pair<unsigned int, float> > saveIdxDr; std::vector<std::pair<unsigned int, float> > saveIdxX;
  std::vector<std::pair<unsigned int, float> > x_j_sort; std::vector<std::pair<unsigned int, float> > dr_j_sort;
  std::map<unsigned int, bool> saveIsInJet; std::map<unsigned int, bool> saveIsCorrectMat;
  int nHadKS = 0;
  for (unsigned int i=0; i<nhad; ++i) {
    if (had_pdgId[i] != 310) continue;
    ++nTotHadKS;
    ++nHadKS;
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
    if (hadTruth_isHadFromTop[i] && hadTruth_nMatched[i] == 2) ++nRealKSFromTop;
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
      // KS inside Jet 
      if (hadTruth_isHadFromTsb[i] == qjMapForMC_[j] && hadTruth_nMatched[i] == 2) {cout << " =====> KS from : " << hadTruth_isHadFromTsb[i] << " , Jet : " << qjMapForMC_[j] << endl;++nRealKSWithJet;}
      if ((jet_tlv.DeltaR(had_tlv) < 0.5) && (had_pt[i]/Jet_pt[j] > 0.15)) {
        if (hadTruth_isHadFromTsb[i] == qjMapForMC_[j]) { // For jet with the same origin as KS
          if (hadTruth_nMatched[i] == 2) ++nRealKSWithJetAndCut;
          b_hadTruth_nMatched_jetMat_vec.push_back(hadTruth_nMatched[i]);
          b_hadTruth_isFrom_jetMat_vec.push_back(hadTruth_isHadFromTsb[i]);
          b_hadTruth_isHadFromTop_jetMat_vec.push_back(hadTruth_isHadFromTop[i]);
          b_hadTruth_isHadFromW_jetMat_vec.push_back(hadTruth_isHadFromW[i]);
          b_hadTruth_isHadFromS_jetMat_vec.push_back(hadTruth_isHadFromS[i]);
          b_hadTruth_isHadFromC_jetMat_vec.push_back(hadTruth_isHadFromC[i]);
          b_hadTruth_isHadFromB_jetMat_vec.push_back(hadTruth_isHadFromB[i]);
          b_hadTruth_qjMapForMC_jetMat_vec.push_back(qjMapForMC_[j]);
          b_hadTruth_d_jetMat_vec.push_back(GetD(had_pt[i], had_eta[i], had_phi[i], had_mass[i], had_x[i], had_y[i], had_z[i]));
          b_hadTruth_x_jetMat_vec.push_back(had_pt[i]/Jet_pt[j]);
          b_hadTruth_dr_jetMat_vec.push_back(jet_tlv.DeltaR(had_tlv));
          b_hadTruth_pt_jetMat_vec.push_back(had_pt[i]);
          b_hadTruth_eta_jetMat_vec.push_back(had_eta[i]);
          b_hadTruth_phi_jetMat_vec.push_back(had_phi[i]);
          b_hadTruth_mass_jetMat_vec.push_back(had_mass[i]);
          b_hadTruth_lxy_jetMat_vec.push_back(had_lxy[i]);
          b_hadTruth_lxySig_jetMat_vec.push_back(had_lxy[i]/had_lxyErr[i]);
          b_hadTruth_angleXY_jetMat_vec.push_back(had_angleXY[i]);
          b_hadTruth_angleXYZ_jetMat_vec.push_back(had_angleXYZ[i]);
          b_hadTruth_chi2_jetMat_vec.push_back(had_chi2[i]);
          b_hadTruth_dca_jetMat_vec.push_back(had_dca[i]);
          b_hadTruth_l3D_jetMat_vec.push_back(had_l3D[i]);
          b_hadTruth_l3DSig_jetMat_vec.push_back(had_l3D[i]/had_l3DErr[i]);
          b_hadTruth_legDR_jetMat_vec.push_back(had_legDR[i]);
          b_hadTruth_pdgId_jetMat_vec.push_back(had_pdgId[i]);
          b_hadTruth_dau1_chi2_jetMat_vec.push_back(had_dau1_chi2[i]);
          b_hadTruth_dau1_ipsigXY_jetMat_vec.push_back(had_dau1_ipsigXY[i]);
          b_hadTruth_dau1_ipsigZ_jetMat_vec.push_back(had_dau1_ipsigZ[i]);
          b_hadTruth_dau1_pt_jetMat_vec.push_back(had_dau1_pt[i]);
          b_hadTruth_dau2_chi2_jetMat_vec.push_back(had_dau2_chi2[i]);
          b_hadTruth_dau2_ipsigXY_jetMat_vec.push_back(had_dau1_ipsigXY[i]);
          b_hadTruth_dau2_ipsigZ_jetMat_vec.push_back(had_dau1_ipsigZ[i]);
          b_hadTruth_dau2_pt_jetMat_vec.push_back(had_dau1_pt[i]);

          saveIdxDr.push_back({i, jet_tlv.DeltaR(had_tlv)});
          saveIdxX.push_back({i, (had_pt[i]/Jet_pt[j])});
	  saveIsInJet.insert({i, true});
          saveIsCorrectMat.insert({i,true});

          b_hadTruth_isFrom_j_vec.push_back(hadTruth_isHadFromTsb[i]);
          b_hadTruth_nMatched_j_vec.push_back(hadTruth_nMatched[i]);
          b_hadTruth_isInJet_j_vec.push_back(true);
          b_hadTruth_isCorrectMat_j_vec.push_back(true);
          b_hadTruth_d_j_vec.push_back(GetD(had_pt[i], had_eta[i], had_phi[i], had_mass[i], had_x[i], had_y[i], had_z[i]));
          b_hadTruth_x_j_vec.push_back(had_pt[i]/Jet_pt[j]);
          b_hadTruth_dr_j_vec.push_back(jet_tlv.DeltaR(had_tlv));
          b_hadTruth_pt_j_vec.push_back(had_pt[i]);
        } else {
          saveIdxDr.push_back({i, jet_tlv.DeltaR(had_tlv)});
          saveIdxX.push_back({i, (had_pt[i]/Jet_pt[j])});
          saveIsInJet.insert({i, true});
          saveIsCorrectMat.insert({i,false});

          b_hadTruth_isFrom_j_vec.push_back(hadTruth_isHadFromTsb[i]);
          b_hadTruth_nMatched_j_vec.push_back(hadTruth_nMatched[i]);
          b_hadTruth_isInJet_j_vec.push_back(true);
          b_hadTruth_isCorrectMat_j_vec.push_back(false);
          b_hadTruth_d_j_vec.push_back(GetD(had_pt[i], had_eta[i], had_phi[i], had_mass[i], had_x[i], had_y[i], had_z[i]));
          b_hadTruth_x_j_vec.push_back(had_pt[i]/Jet_pt[j]);
          b_hadTruth_dr_j_vec.push_back(jet_tlv.DeltaR(had_tlv));
          b_hadTruth_pt_j_vec.push_back(had_pt[i]);
        }
      } else {
        saveIdxDr.push_back({i, jet_tlv.DeltaR(had_tlv)});
        saveIdxX.push_back({i, (had_pt[i]/Jet_pt[j])});
        saveIsInJet.insert({i, false});
        saveIsCorrectMat.insert({i,false});

        b_hadTruth_isFrom_j_vec.push_back(hadTruth_isHadFromTsb[i]);
        b_hadTruth_nMatched_j_vec.push_back(hadTruth_nMatched[i]);
        b_hadTruth_isInJet_j_vec.push_back(false);
        b_hadTruth_isCorrectMat_j_vec.push_back(false);
        b_hadTruth_d_j_vec.push_back(GetD(had_pt[i], had_eta[i], had_phi[i], had_mass[i], had_x[i], had_y[i], had_z[i]));
        b_hadTruth_x_j_vec.push_back(had_pt[i]/Jet_pt[j]);
        b_hadTruth_dr_j_vec.push_back(jet_tlv.DeltaR(had_tlv));
        b_hadTruth_pt_j_vec.push_back(had_pt[i]);
      }
    }
    for (unsigned int j=0; j<nGenJet; ++j) { // Loop for all of recoJet
      TLorentzVector jet_tlv;
      jet_tlv.SetPtEtaPhiM(GenJet_pt[j], GenJet_eta[j], GenJet_phi[j], GenJet_mass[j]);
      // KS inside GenJet  
      if ((jet_tlv.DeltaR(had_tlv) < 0.5) && (had_pt[i]/GenJet_pt[j] > 0.15)) {
        if (hadTruth_isHadFromTsb[i] == qgjMapForMC_[j]) { // For jet with the same origin as KS
          if (hadTruth_nMatched[i] == 2) ++nRealKSWithJetAndCut;
          b_hadTruth_nMatched_gjetMat_vec.push_back(hadTruth_nMatched[i]);
          b_hadTruth_isFrom_gjetMat_vec.push_back(hadTruth_isHadFromTsb[i]);
          b_hadTruth_isHadFromTop_gjetMat_vec.push_back(hadTruth_isHadFromTop[i]);
          b_hadTruth_isHadFromW_gjetMat_vec.push_back(hadTruth_isHadFromW[i]);
          b_hadTruth_isHadFromS_gjetMat_vec.push_back(hadTruth_isHadFromS[i]);
          b_hadTruth_isHadFromC_gjetMat_vec.push_back(hadTruth_isHadFromC[i]);
          b_hadTruth_isHadFromB_gjetMat_vec.push_back(hadTruth_isHadFromB[i]);
          b_hadTruth_qjMapForMC_gjetMat_vec.push_back(qjMapForMC_[j]);
          b_hadTruth_d_gjetMat_vec.push_back(GetD(had_pt[i], had_eta[i], had_phi[i], had_mass[i], had_x[i], had_y[i], had_z[i]));
          b_hadTruth_x_gjetMat_vec.push_back(had_pt[i]/GenJet_pt[j]);
          b_hadTruth_dr_gjetMat_vec.push_back(jet_tlv.DeltaR(had_tlv));
          b_hadTruth_pt_gjetMat_vec.push_back(had_pt[i]);
          b_hadTruth_eta_gjetMat_vec.push_back(had_eta[i]);
          b_hadTruth_phi_gjetMat_vec.push_back(had_phi[i]);
          b_hadTruth_mass_gjetMat_vec.push_back(had_mass[i]);
          b_hadTruth_lxy_gjetMat_vec.push_back(had_lxy[i]);
          b_hadTruth_lxySig_gjetMat_vec.push_back(had_lxy[i]/had_lxyErr[i]);
          b_hadTruth_angleXY_gjetMat_vec.push_back(had_angleXY[i]);
          b_hadTruth_angleXYZ_gjetMat_vec.push_back(had_angleXYZ[i]);
          b_hadTruth_chi2_gjetMat_vec.push_back(had_chi2[i]);
          b_hadTruth_dca_gjetMat_vec.push_back(had_dca[i]);
          b_hadTruth_l3D_gjetMat_vec.push_back(had_l3D[i]);
          b_hadTruth_l3DSig_gjetMat_vec.push_back(had_l3D[i]/had_l3DErr[i]);
          b_hadTruth_legDR_gjetMat_vec.push_back(had_legDR[i]);
          b_hadTruth_pdgId_gjetMat_vec.push_back(had_pdgId[i]);
          b_hadTruth_dau1_chi2_gjetMat_vec.push_back(had_dau1_chi2[i]);
          b_hadTruth_dau1_ipsigXY_gjetMat_vec.push_back(had_dau1_ipsigXY[i]);
          b_hadTruth_dau1_ipsigZ_gjetMat_vec.push_back(had_dau1_ipsigZ[i]);
          b_hadTruth_dau1_pt_gjetMat_vec.push_back(had_dau1_pt[i]);
          b_hadTruth_dau2_chi2_gjetMat_vec.push_back(had_dau2_chi2[i]);
          b_hadTruth_dau2_ipsigXY_gjetMat_vec.push_back(had_dau1_ipsigXY[i]);
          b_hadTruth_dau2_ipsigZ_gjetMat_vec.push_back(had_dau1_ipsigZ[i]);
          b_hadTruth_dau2_pt_gjetMat_vec.push_back(had_dau1_pt[i]);

          b_hadTruth_isFrom_j_vec.push_back(hadTruth_isHadFromTsb[i]);
          b_hadTruth_nMatched_j_vec.push_back(hadTruth_nMatched[i]);
          b_hadTruth_isInJet_gj_vec.push_back(true);
          b_hadTruth_isCorrectMat_gj_vec.push_back(true);
          b_hadTruth_d_gj_vec.push_back(GetD(had_pt[i], had_eta[i], had_phi[i], had_mass[i], had_x[i], had_y[i], had_z[i]));
          b_hadTruth_x_gj_vec.push_back(had_pt[i]/Jet_pt[j]);
          b_hadTruth_dr_gj_vec.push_back(jet_tlv.DeltaR(had_tlv));
          b_hadTruth_pt_gj_vec.push_back(had_pt[i]);
        } else {
          b_hadTruth_isFrom_j_vec.push_back(hadTruth_isHadFromTsb[i]);
          b_hadTruth_nMatched_j_vec.push_back(hadTruth_nMatched[i]);
          b_hadTruth_isInJet_gj_vec.push_back(true);
          b_hadTruth_isCorrectMat_gj_vec.push_back(false);
          b_hadTruth_d_gj_vec.push_back(GetD(had_pt[i], had_eta[i], had_phi[i], had_mass[i], had_x[i], had_y[i], had_z[i]));
          b_hadTruth_x_gj_vec.push_back(had_pt[i]/Jet_pt[j]);
          b_hadTruth_dr_gj_vec.push_back(jet_tlv.DeltaR(had_tlv));
          b_hadTruth_pt_gj_vec.push_back(had_pt[i]);
        }
      } else {
        b_hadTruth_isFrom_gj_vec.push_back(hadTruth_isHadFromTsb[i]);
        b_hadTruth_nMatched_gj_vec.push_back(hadTruth_nMatched[i]);
        b_hadTruth_isInJet_gj_vec.push_back(false);
        b_hadTruth_isCorrectMat_gj_vec.push_back(false);
        b_hadTruth_d_gj_vec.push_back(GetD(had_pt[i], had_eta[i], had_phi[i], had_mass[i], had_x[i], had_y[i], had_z[i]));
        b_hadTruth_x_gj_vec.push_back(had_pt[i]/Jet_pt[j]);
        b_hadTruth_dr_gj_vec.push_back(jet_tlv.DeltaR(had_tlv));
        b_hadTruth_pt_gj_vec.push_back(had_pt[i]);
      }
    }
  }
  if (nHadKS != 0) {
    for (unsigned int i=0; i<nJet; ++i) {
      for (int j=0; j<nHadKS; ++j) {
        dr_j_sort.push_back(saveIdxDr[(j*nJet)+i] );  
        x_j_sort.push_back(saveIdxX[(j*nJet)+i] );
      }
    }
    for (unsigned int i=0; i<nJet; ++i) {
      sort(dr_j_sort.begin()+nHadKS*i, dr_j_sort.begin()+nHadKS*(i+1), [](std::pair<unsigned int, float> a, std::pair<unsigned int, float> b) { return (a.second < b.second); } );
      sort(x_j_sort.begin()+nHadKS*i, x_j_sort.begin()+nHadKS*(i+1), [](std::pair<unsigned int, float> a, std::pair<unsigned int, float> b) { return (a.second > b.second); } );
    }
    for (unsigned int i=0; i<nJet; ++i) {
      auto idx_x = x_j_sort[nHadKS*i].first;
      auto idx_dr = dr_j_sort[nHadKS*i].first;
      TLorentzVector had_tlv; // For b_hadTruth_x_highest_j_vec
      had_tlv.SetPtEtaPhiM(had_pt[idx_x], had_eta[idx_x], had_phi[idx_x], had_mass[idx_x]);
      TLorentzVector jet_tlv;
      jet_tlv.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
  
      b_hadTruth_idx_closest_j_vec.push_back(idx_dr);
      b_hadTruth_nMatched_closest_j_vec.push_back(hadTruth_nMatched[idx_dr]);
      b_hadTruth_isFrom_closest_j_vec.push_back(hadTruth_isHadFromTsb[idx_dr]);
      b_hadTruth_isHadFromTop_closest_j_vec.push_back(hadTruth_isHadFromTop[idx_dr]);
      b_hadTruth_isInJet_closest_j_vec.push_back(saveIsInJet[idx_dr]);
      b_hadTruth_isCorrectMat_closest_j_vec.push_back(saveIsCorrectMat[idx_dr]);
      b_hadTruth_d_closest_j_vec.push_back(GetD(had_pt[idx_dr], had_eta[idx_dr], had_phi[idx_dr], had_mass[idx_dr], had_x[idx_dr], had_y[idx_dr], had_z[idx_dr]));
      b_hadTruth_x_closest_j_vec.push_back(had_pt[idx_dr]/Jet_pt[i]);
      b_hadTruth_dr_closest_j_vec.push_back(dr_j_sort[nHadKS*i].second);
      b_hadTruth_pt_closest_j_vec.push_back(had_pt[idx_dr]);
      b_hadTruth_eta_closest_j_vec.push_back(had_eta[idx_dr]);
      b_hadTruth_phi_closest_j_vec.push_back(had_phi[idx_dr]);
      b_hadTruth_mass_closest_j_vec.push_back(had_mass[idx_dr]);
  
      b_hadTruth_idx_highest_j_vec.push_back(idx_x);
      b_hadTruth_nMatched_highest_j_vec.push_back(hadTruth_nMatched[idx_x]);
      b_hadTruth_isFrom_highest_j_vec.push_back(hadTruth_isHadFromTsb[idx_x]);
      b_hadTruth_isHadFromTop_highest_j_vec.push_back(hadTruth_isHadFromTop[idx_x]);
      b_hadTruth_isInJet_highest_j_vec.push_back(saveIsInJet[idx_x]);
      b_hadTruth_isCorrectMat_highest_j_vec.push_back(saveIsCorrectMat[idx_x]);
      b_hadTruth_d_highest_j_vec.push_back(GetD(had_pt[idx_x], had_eta[idx_x], had_phi[idx_x], had_mass[idx_x], had_x[idx_x], had_y[idx_x], had_z[idx_x]));
      b_hadTruth_x_highest_j_vec.push_back(x_j_sort[nHadKS*i].second);
      b_hadTruth_dr_highest_j_vec.push_back(jet_tlv.DeltaR(had_tlv));
      b_hadTruth_pt_highest_j_vec.push_back(had_pt[idx_x]);
      b_hadTruth_eta_highest_j_vec.push_back(had_eta[idx_x]);
      b_hadTruth_phi_highest_j_vec.push_back(had_phi[idx_x]);
      b_hadTruth_mass_highest_j_vec.push_back(had_mass[idx_x]);
    }
    auto closest_dr = std::min_element(b_hadTruth_dr_closest_j_vec.begin(), b_hadTruth_dr_closest_j_vec.end());
    auto element_closest_dr = std::distance(b_hadTruth_dr_closest_j_vec.begin(), closest_dr);
    auto highest_x = std::max_element(b_hadTruth_x_highest_j_vec.begin(), b_hadTruth_x_highest_j_vec.end());
    auto element_highest_x = std::distance(b_hadTruth_x_highest_j_vec.begin(), highest_x);

    auto idx_dr = b_hadTruth_idx_closest_j_vec[element_closest_dr];
    auto idx_x = b_hadTruth_idx_highest_j_vec[element_highest_x];

    TLorentzVector had_tlv; // For b_hadTruth_x_highest_j_vec
    had_tlv.SetPtEtaPhiM(had_pt[idx_x], had_eta[idx_x], had_phi[idx_x], had_mass[idx_x]);
    TLorentzVector jet_tlv;
    jet_tlv.SetPtEtaPhiM(Jet_pt[element_highest_x], Jet_eta[element_highest_x], Jet_phi[element_highest_x], Jet_mass[element_highest_x]);

    b_hadTruth_nMatched_closest_j = hadTruth_nMatched[idx_dr];
    b_hadTruth_isFrom_closest_j = hadTruth_isHadFromTsb[idx_dr];
    b_hadTruth_isHadFromTop_closest_j = hadTruth_isHadFromTop[idx_dr];
    b_hadTruth_isInJet_closest_j = saveIsInJet[idx_dr];
    b_hadTruth_isCorrectMat_closest_j = saveIsCorrectMat[idx_dr];
    b_hadTruth_d_closest_j = GetD(had_pt[idx_dr], had_eta[idx_dr], had_phi[idx_dr], had_mass[idx_dr], had_x[idx_dr], had_y[idx_dr], had_z[idx_dr]);
    b_hadTruth_x_closest_j = had_pt[idx_dr]/Jet_pt[element_closest_dr];
    b_hadTruth_dr_closest_j = dr_j_sort[nHadKS*element_closest_dr].second;
    b_hadTruth_pt_closest_j = had_pt[idx_dr];
    b_hadTruth_eta_closest_j = had_eta[idx_dr];
    b_hadTruth_phi_closest_j = had_phi[idx_dr];
    b_hadTruth_mass_closest_j = had_mass[idx_dr];

    b_hadTruth_nMatched_highest_j = hadTruth_nMatched[idx_x];
    b_hadTruth_isFrom_highest_j = hadTruth_isHadFromTsb[idx_x];
    b_hadTruth_isHadFromTop_closest_j = hadTruth_isHadFromTop[idx_x];
    b_hadTruth_isInJet_highest_j = saveIsInJet[idx_x];
    b_hadTruth_isCorrectMat_highest_j = saveIsCorrectMat[idx_x];
    b_hadTruth_d_highest_j = GetD(had_pt[idx_x], had_eta[idx_x], had_phi[idx_x], had_mass[idx_x], had_x[idx_x], had_y[idx_x], had_z[idx_x]);
    b_hadTruth_x_highest_j = x_j_sort[nHadKS*element_highest_x].second;
    b_hadTruth_dr_highest_j = jet_tlv.DeltaR(had_tlv);
    b_hadTruth_pt_highest_j = had_pt[idx_x];
    b_hadTruth_eta_highest_j = had_eta[idx_x];
    b_hadTruth_phi_highest_j = had_phi[idx_x];
    b_hadTruth_mass_highest_j = had_mass[idx_x];
  } 
  return 1;
}

void vtsAnalyser::CollectVar() {
  auto muons = muonSelection();
  auto elecs = elecSelection();
  if (elecs.size() + muons.size() != 2) return; // this is for running correctly function regardless of return value of EventSelection()

  b_MET_pt = MET_pt;
  b_MET_phi = MET_phi;
  b_MET_sumEt = MET_sumEt;

  b_lep_pt_vec.push_back(b_lep1.Pt());
  b_lep_eta_vec.push_back(b_lep1.Eta());
  b_lep_phi_vec.push_back(b_lep1.Phi());
  b_lep_mass_vec.push_back(b_lep1.M());
  b_lep_pt_vec.push_back(b_lep2.Pt());
  b_lep_eta_vec.push_back(b_lep2.Eta());
  b_lep_phi_vec.push_back(b_lep2.Phi());
  b_lep_mass_vec.push_back(b_lep2.M());

  b_dilep_pt_vec.push_back(b_dilep.Pt());
  b_dilep_eta_vec.push_back(b_dilep.Eta());
  b_dilep_phi_vec.push_back(b_dilep.Phi());
  b_dilep_mass_vec.push_back(b_dilep.M());

  TLorentzVector elec;
  TLorentzVector mu;
  if (abs(b_lep1_pid) == 11) {
    elec = b_lep1;
    b_elec_pt_vec.push_back(elec.Pt());
    b_elec_eta_vec.push_back(elec.Eta());
    b_elec_phi_vec.push_back(elec.Phi());
    b_elec_mass_vec.push_back(elec.M());
    b_elec_dxy_vec.push_back(Electron_dxy[b_lep1_idx]);
    b_lep_dxy_vec.push_back(Electron_dxy[b_lep1_idx]);
  } else {
    mu = b_lep1;
    b_mu_pt_vec.push_back(mu.Pt());
    b_mu_eta_vec.push_back(mu.Eta());
    b_mu_phi_vec.push_back(mu.Phi());
    b_mu_mass_vec.push_back(mu.M());
    b_mu_dxy_vec.push_back(Muon_dxy[b_lep1_idx]);
    b_lep_dxy_vec.push_back(Muon_dxy[b_lep1_idx]);
  }
  if (abs(b_lep2_pid) == 11) {
    elec = b_lep2;
    b_elec_pt_vec.push_back(elec.Pt());
    b_elec_eta_vec.push_back(elec.Eta());
    b_elec_phi_vec.push_back(elec.Phi());
    b_elec_mass_vec.push_back(elec.M());
    b_elec_dxy_vec.push_back(Electron_dxy[b_lep2_idx]);
    b_lep_dxy_vec.push_back(Electron_dxy[b_lep2_idx]);
  } else {
    mu = b_lep2;
    b_mu_pt_vec.push_back(mu.Pt());
    b_mu_eta_vec.push_back(mu.Eta());
    b_mu_phi_vec.push_back(mu.Phi());
    b_mu_mass_vec.push_back(mu.M());
    b_mu_dxy_vec.push_back(Muon_dxy[b_lep2_idx]);
    b_lep_dxy_vec.push_back(Muon_dxy[b_lep2_idx]);
  }
}

