#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <TLorentzVector.h>
#include <sys/stat.h>
#include "nano/analysis/interface/vtsAnalyser.h"
#include "TFile.h"
#include "TTree.h"

#include "DataFormats/HepMCCandidate/interface/GenStatusFlags.h"

using namespace std;

string getFileName(const string& s)
{
 char sep = '/';
 size_t i = s.rfind(sep, s.length());
 if (i != string::npos) {
    return s.substr(i+1, s.length() - i);
 }
 return "";
}

string getDir (const string& path)
{
  size_t found = path.find_last_of("/\\");
  return path.substr(0, found);
}

string getType (const string& path)
{
  size_t found = path.find("tt");
  return path.substr(found);
}

long getFileSize(const char* filename)
{
  struct stat st;
  if (stat(filename, &st) == 0) {
    return st.st_size;
  }
  return 0;
}

int main(int argc, char* argv[])
{
  string hostDir = "root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/";

  if (argc <= 1) {
    cout << "no input file is specified. running with default file." << endl;
    auto inFile = TFile::Open("/xrootd/store/group/nanoAOD/run2_2016v5/tsw/nanoAOD_111.root", "READ");
    /* auto inFile = TFile::Open("/xrootd/store/group/nanoAOD/run2_2016v5/tsW_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6/180611_131219/0000/nanoAOD_000.root"); */
    auto inTree = (TTree*) inFile->Get("Events");
    vtsAnalyser ana(inTree,inTree,0,true,false,false,false);
    ana.setOutput("nanotree.root");
    ana.Loop();
  } else {
    string jobName    = string(argv[1]);
    string sampleName = string(argv[2]);
/*
    cout << "sampleName : " << sampleName << endl;
    cout << "isMC : " << isMC << endl;
    cout << "isGenericMC : " << isGenericMC << endl;
*/
    string outFileDir = hostDir + getenv("USER") + "/" + jobName + "/" + sampleName;
    for (Int_t i = 3; i < argc; i++) {
      auto inFileName = argv[i];
      auto fileName = getFileName(argv[i]);
      auto dirName = getDir(getDir(argv[i]));
      /* auto sampleType = getType(dirName); */

      Bool_t isMC = false;
      Bool_t isGenericMC = false;

      if (string(inFileName).find("Run2016") == std::string::npos && (string(inFileName).find("Run2017") == std::string::npos) && (string(inFileName).find("Run2018") == std::string::npos) ) {
        isMC = true;
      }
      if (string(inFileName).find("NANOAOD") != std::string::npos) isGenericMC = false;
      else isGenericMC = true;
      
      if (!isMC) { 
        Bool_t isDL = false;
        if (string(inFileName).find("DoubleMuon") != std::string::npos)     isDL = true;
        if (string(inFileName).find("DoubleEG") != std::string::npos)       isDL = true;
        if (string(inFileName).find("MuonEG") != std::string::npos)         isDL = true;
        cout << "inFileName : " << inFileName << endl;
        Bool_t isSLE = (string(inFileName).find("SingleElectron") != std::string::npos);
        Bool_t isSLM = (string(inFileName).find("SingleMuon") != std::string::npos);
        TFile *inFile = TFile::Open(inFileName, "READ");
        TTree *inTree = (TTree*) inFile->Get("Events");
        vtsAnalyser ana(inTree, inTree, inTree, isMC, isDL, isSLE, isSLM);
//        string outFileName = outFileDir+"/nanotree_"+to_string(i-3)+".root";
        string outFileName = outFileDir+"/rd_"+fileName;
        ana.setOutput(outFileName);
        ana.Loop();
        cout << inFileName << " : Analysis Complete " << endl;
      }
      if (isMC && isGenericMC) {
        cout << "inFileName : " << inFileName << endl;
        cout << " input file is not custumized ( e.g tt###j_* ) sample" << endl;
        if (string(inFileName).find("run2") != std::string::npos) { 
          TFile *inFile = TFile::Open(inFileName, "READ");
          TTree *inTree = (TTree*) inFile->Get("Events");
          vtsAnalyser ana(inTree,isMC,false,false,false,isGenericMC);
          string outFileName = outFileDir+"/mc_"+fileName;
          ana.setOutput(outFileName);
          ana.Loop();
          cout << inFileName << " : Analysis Complete " << endl;
        } else if (string(inFileName).find("nanoTest") != std::string::npos) { // Added temporarily
          TFile *inFile = TFile::Open(inFileName, "READ");
          TTree *inTree = (TTree*) inFile->Get("Events");
          vtsAnalyser ana(inTree,isMC,false,false,false,isGenericMC);
          string outFileName = outFileDir+"/mc_"+fileName;
          ana.setOutput(outFileName);
          ana.Loop();
          cout << inFileName << " : Analysis Complete " << endl;
        } else {
          cout << " input file has to be nanoAOD " << endl;
          continue;
        }
      } else if (isMC) {
        TFile *inFile = TFile::Open(inFileName, "READ");
        TTree *inTree = (TTree*) inFile->Get("Events");
        TString hadFileName      = dirName + "/HADAOD/" + fileName;
        TString hadTruthFileName = dirName + "/HADTRUTHAOD/" + fileName;
        TFile *hadFile      = TFile::Open(hadFileName, "READ");
        TFile *hadTruthFile = TFile::Open(hadTruthFileName, "READ");
        TTree *hadTree(0);      
        TTree *hadTruthTree(0);
        if (hadFile != NULL) hadTree = (TTree*) hadFile->Get("Events");
        if (hadTruthFile != NULL) {
          cout << ">>>>> Use NANOTREE and HADTRUTHTREE <<<<< " << endl;
          hadTruthTree = (TTree*) hadTruthFile->Get("Events");
        }
        cout << "dirName : " << dirName << " fileName : " << fileName << endl;
        vtsAnalyser ana(inTree,hadTree,hadTruthTree,isMC,false,false,false); // you don't need to use hadTree if hadTruthTree isn't NULL
        string outFileName;
        if (string(inFileName).find("herwig") == std::string::npos) outFileName = outFileDir+"/pythia_nanotree_"+fileName;
        else if (string(inFileName).find("herwig") != std::string::npos) outFileName = outFileDir+"/herwig_nanotree_"+fileName;
        else outFileName = outFileDir+"/nanotree_"+fileName;
        ana.setOutput(outFileName);
        ana.Loop();
        cout << inFileName << " : Analysis Complete " << endl;
        
      }
    }
  }
}

void vtsAnalyser::Loop() {
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntries();
  /* Events loop */
  for (Long64_t iev=0; iev<nentries; iev++) {
    b_had_start = -1; b_had_end = -1;
    b_jet_start = -1; b_jet_end = -1;
    fChain->GetEntry(iev);
    if (h_fChain) h_fChain->GetEntry(iev);
    if (ht_fChain) ht_fChain->GetEntry(iev);
//    cout << "event : " << iev << endl;
    if (iev%10000 == 0) cout << iev << "/" << nentries << endl;
    ResetBranch();
    b_nJet = nJet;
    b_nSelJet = jetSelection().size();
    EventSelection();
    if (b_step >= 4) { // -1 : No event selection, 0 : PV cut, reco lepton cut and so on, 1~4 : step 1 ~ 4
      b_passedEvent = true; 
      b_nSelJetEv = jetSelection().size();
      MatchingForMC();
      GenAnalysis();
      RecAnalysis();
      JetAnalysis();
      CollectVar();
      FillJetTreeForTMVA();
    } else b_passedEvent = false; 
    m_tree->Fill();
//    cout << " jet entry chk : " << b_jet_start << " " << b_jet_end << endl;
//    cout << " had entry chk : " << b_had_start << " " << b_had_end << endl;
  }
}

void vtsAnalyser::setOutput(std::string outFileName) {
  m_output = TFile::Open(outFileName.c_str(), "recreate");
  m_tree = new TTree("event", "event");
  m_hadtrForTMVA = new TTree("MVA_had", "MVA_had");
  m_jettrForTMVA = new TTree("MVA_jet", "MVA_jet");

  MakeBranch();
  SetMVAReader();

  h_nevents = new TH1D("nevents", "nevents", 1, 0, 1);
  h_genweights = new TH1D("genweight", "genweight", 1, 0, 1);
  h_weights = new TH1D("weight", "weight", 1, 0, 1);
  h_cutFlow = new TH1D("cutflow", "cutflow", 11, -0.5, 10.5);
}

void vtsAnalyser::MakeBranch() {
  m_hadtrForTMVA->Branch("pdgId",        &b_Rec_pdgId,        "pdgId/I");
  m_hadtrForTMVA->Branch("nMatched",     &b_Rec_nMatched,     "nMatched/I");
  m_hadtrForTMVA->Branch("isFrom",       &b_Rec_isFrom,       "isFrom/I");
  m_hadtrForTMVA->Branch("isHadFromTop", &b_Rec_isHadFromTop, "isHadFromTop/O");
  m_hadtrForTMVA->Branch("isHadFromW",   &b_Rec_isHadFromW,   "isHadFromW/O");
  m_hadtrForTMVA->Branch("isHadFromS",   &b_Rec_isHadFromS,   "isHadFromS/O");
  m_hadtrForTMVA->Branch("isHadFromC",   &b_Rec_isHadFromC,   "isHadFromC/O");
  m_hadtrForTMVA->Branch("isHadFromB",   &b_Rec_isHadFromB,   "isHadFromB/O");
  m_hadtrForTMVA->Branch("d",            &b_Rec_d,            "d/F");
  m_hadtrForTMVA->Branch("pt",           &b_Rec_pt,           "pt/F");
  m_hadtrForTMVA->Branch("eta",          &b_Rec_eta,          "eta/F");
  m_hadtrForTMVA->Branch("phi",          &b_Rec_phi,          "phi/F");
  m_hadtrForTMVA->Branch("mass",         &b_Rec_mass,         "mass/F");
  m_hadtrForTMVA->Branch("lxy",          &b_Rec_lxy,          "lxy/F");
  m_hadtrForTMVA->Branch("lxySig",       &b_Rec_lxySig,       "lxySig/F");
  m_hadtrForTMVA->Branch("l3D",          &b_Rec_l3D,          "l3D/F");
  m_hadtrForTMVA->Branch("l3DSig",       &b_Rec_l3DSig,       "l3DSig/F");
  m_hadtrForTMVA->Branch("legDR",        &b_Rec_legDR,        "legDR/F");
  m_hadtrForTMVA->Branch("angleXY",      &b_Rec_angleXY,      "angleXY/F");
  m_hadtrForTMVA->Branch("angleXYZ",     &b_Rec_angleXYZ,     "angleXYZ/F");
  m_hadtrForTMVA->Branch("chi2",         &b_Rec_chi2,         "chi2/F");
  m_hadtrForTMVA->Branch("dca",          &b_Rec_dca,          "dca/F");
  m_hadtrForTMVA->Branch("dau1_chi2",    &b_Rec_dau1_chi2,    "dau1_chi2/F");
  m_hadtrForTMVA->Branch("dau1_ipsigXY", &b_Rec_dau1_ipsigXY, "dau1_ipsigXY/F");
  m_hadtrForTMVA->Branch("dau1_ipsigZ",  &b_Rec_dau1_ipsigZ,  "dau1_ipsigZ/F");
  m_hadtrForTMVA->Branch("dau1_pt",      &b_Rec_dau2_pt,      "dau2_pt/F");
  m_hadtrForTMVA->Branch("dau2_chi2",    &b_Rec_dau2_chi2,    "dau2_chi2/F");
  m_hadtrForTMVA->Branch("dau2_ipsigXY", &b_Rec_dau2_ipsigXY, "dau2_ipsigXY/F");
  m_hadtrForTMVA->Branch("dau2_ipsigZ",  &b_Rec_dau2_ipsigZ,  "dau2_ipsigZ/F");
  m_hadtrForTMVA->Branch("dau2_pt",      &b_Rec_dau2_pt,      "dau2_pt/F");
  m_hadtrForTMVA->Branch("bdt_score",    &b_Rec_bdt_score,    "bdt_score/F");

  m_jettrForTMVA->Branch("isSJet",         &b_isSJet,         "isSJet/I");
  m_jettrForTMVA->Branch("isBJet",         &b_isBJet,         "isBJet/I");
  m_jettrForTMVA->Branch("isHighest",      &b_isHighest,      "isHighest/I");
  m_jettrForTMVA->Branch("isClosestToLep", &b_isClosestToLep, "isClosestToLep/I");
  m_jettrForTMVA->Branch("cmult",          &b_cmult,          "cmult/F");
  m_jettrForTMVA->Branch("nmult",          &b_nmult,          "nmult/F");
  m_jettrForTMVA->Branch("pt",             &b_pt,             "pt/F");
  m_jettrForTMVA->Branch("eta",            &b_eta,            "eta/F");
  m_jettrForTMVA->Branch("phi",            &b_phi,            "phi/F");
  m_jettrForTMVA->Branch("mass",           &b_mass,           "mass/F");
  m_jettrForTMVA->Branch("c_x1",           &b_c_x1,           "c_x1/F");
  m_jettrForTMVA->Branch("c_x2",           &b_c_x2,           "c_x2/F");
  m_jettrForTMVA->Branch("c_x3",           &b_c_x3,           "c_x3/F");
  m_jettrForTMVA->Branch("n_x1",           &b_n_x1,           "n_x1/F");
  m_jettrForTMVA->Branch("n_x2",           &b_n_x2,           "n_x2/F");
  m_jettrForTMVA->Branch("n_x3",           &b_n_x3,           "n_x3/F");
  m_jettrForTMVA->Branch("axis1",          &b_axis1,          "axis1/F");
  m_jettrForTMVA->Branch("axis2",          &b_axis2,          "axis2/F");
  m_jettrForTMVA->Branch("ptD",            &b_ptD,            "ptD/F");
  m_jettrForTMVA->Branch("area",           &b_area,           "area/F");
  m_jettrForTMVA->Branch("CSVV2",          &b_CSVV2,          "CSVV2/F");

  m_jettrForTMVA->Branch("dau_pt",      &b_dau_pt,      "dau_pt[350]/F");
  m_jettrForTMVA->Branch("dau_eta",     &b_dau_eta,     "dau_eta[350]/F");
  m_jettrForTMVA->Branch("dau_phi",     &b_dau_phi,     "dau_phi[350]/F");
  m_jettrForTMVA->Branch("dau_charge",  &b_dau_charge,  "dau_charge[350]/I");

  m_jettrForTMVA->Branch("Jet_bdt_score",   &b_Jet_bdt_score,   "Jet_bdt_score/F");
  m_jettrForTMVA->Branch("JKS_bdt_score",   &b_JKS_bdt_score,   "JKS_bdt_score/F");

  m_jettrForTMVA->Branch("KS_idx"     ,     &b_KS_idx,          "KS_idx/I");
  m_jettrForTMVA->Branch("KS_nMatched",     &b_KS_nMatched,     "KS_nMatched/I");
  m_jettrForTMVA->Branch("KS_isFrom",       &b_KS_isFrom,       "KS_isFrom/I");
  m_jettrForTMVA->Branch("KS_isHadFromTop", &b_KS_isHadFromTop, "KS_isHadFromTop/O");
  m_jettrForTMVA->Branch("KS_isHadFromW",   &b_KS_isHadFromW,   "KS_isHadFromW/O");
  m_jettrForTMVA->Branch("KS_isHadFromS",   &b_KS_isHadFromS,   "KS_isHadFromS/O");
  m_jettrForTMVA->Branch("KS_isHadFromC",   &b_KS_isHadFromC,   "KS_isHadFromC/O");
  m_jettrForTMVA->Branch("KS_isHadFromB",   &b_KS_isHadFromB,   "KS_isHadFromB/O");
  m_jettrForTMVA->Branch("KS_d",            &b_KS_d,            "KS_d/F");
  m_jettrForTMVA->Branch("KS_pt",           &b_KS_pt,           "KS_pt/F");
  m_jettrForTMVA->Branch("KS_eta",          &b_KS_eta,          "KS_eta/F");
  m_jettrForTMVA->Branch("KS_phi",          &b_KS_phi,          "KS_phi/F");
  m_jettrForTMVA->Branch("KS_mass",         &b_KS_mass,         "KS_mass/F");
  m_jettrForTMVA->Branch("KS_lxy",          &b_KS_lxy,          "KS_lxy/F");
  m_jettrForTMVA->Branch("KS_lxySig",       &b_KS_lxySig,       "KS_lxySig/F");
  m_jettrForTMVA->Branch("KS_l3D",          &b_KS_l3D,          "KS_l3D/F");
  m_jettrForTMVA->Branch("KS_l3DSig",       &b_KS_l3DSig,       "KS_l3DSig/F");
  m_jettrForTMVA->Branch("KS_legDR",        &b_KS_legDR,        "KS_legDR/F");
  m_jettrForTMVA->Branch("KS_angleXY",      &b_KS_angleXY,      "KS_angleXY/F");
  m_jettrForTMVA->Branch("KS_angleXYZ",     &b_KS_angleXYZ,     "KS_angleXYZ/F");
  m_jettrForTMVA->Branch("KS_chi2",         &b_KS_chi2,         "KS_chi2/F");
  m_jettrForTMVA->Branch("KS_dca",          &b_KS_dca,          "KS_dca/F");
  m_jettrForTMVA->Branch("KS_dau1_chi2",    &b_KS_dau1_chi2,    "KS_dau1_chi2/F");
  m_jettrForTMVA->Branch("KS_dau1_ipsigXY", &b_KS_dau1_ipsigXY, "KS_dau1_ipsigXY/F");
  m_jettrForTMVA->Branch("KS_dau1_ipsigZ",  &b_KS_dau1_ipsigZ,  "KS_dau1_ipsigZ/F");
  m_jettrForTMVA->Branch("KS_dau1_pt",      &b_KS_dau2_pt,      "KS_dau2_pt/F");
  m_jettrForTMVA->Branch("KS_dau2_chi2",    &b_KS_dau2_chi2,    "KS_dau2_chi2/F");
  m_jettrForTMVA->Branch("KS_dau2_ipsigXY", &b_KS_dau2_ipsigXY, "KS_dau2_ipsigXY/F");
  m_jettrForTMVA->Branch("KS_dau2_ipsigZ",  &b_KS_dau2_ipsigZ,  "KS_dau2_ipsigZ/F");
  m_jettrForTMVA->Branch("KS_dau2_pt",      &b_KS_dau2_pt,      "KS_dau2_pt/F");
  m_jettrForTMVA->Branch("KS_dr",           &b_KS_dr,           "KS_dr/F");
  m_jettrForTMVA->Branch("KS_x",            &b_KS_x,            "KS_x/F");
  m_jettrForTMVA->Branch("KS_best_bdt",     &b_KS_best_bdt,     "KS_best_bdt/F");
 
  #define Branch_(type, name, suffix) m_tree->Branch(#name, &(b_##name), #name "/" #suffix);
  #define BranchI(name) Branch_(Int_t, name, I)
  #define BranchF(name) Branch_(Float_t, name, F)
  #define BranchO(name) Branch_(Bool_t, name, O)
  #define BranchV_(type, name) m_tree->Branch(#name, "vector<"#type">", &(b_##name));
  #define BranchVI(name) BranchV_(Int_t, name); 
  #define BranchVF(name) BranchV_(Float_t, name);
  #define BranchVO(name) BranchV_(Bool_t, name);
  #define BranchTLV(name) m_tree->Branch(#name, "TLorentzVector", &(b_##name));

  BranchI(nvertex);   BranchI(channel); BranchI(njet)       BranchF(met);     BranchI(step); BranchO(passedEvent); 
  BranchI(nJet);      BranchI(nSelJet); BranchI(nSelJetEv); // njet is not nJet
  BranchI(jet_start); BranchI(jet_end); BranchI(had_start); BranchI(had_end);
  BranchI(hadTruth_nMatched);     BranchI(hadTruth_nTrueDau); 
  BranchO(hadTruth_isHadFromTop); BranchI(hadTruth_isHadFromTsb); BranchO(hadTruth_isHadFromW); BranchO(hadTruth_isHadFromS); BranchO(hadTruth_isHadFromC); BranchO(hadTruth_isHadFromB);

  /* weight */
  BranchF(genweight);   BranchF(puweight); 
  BranchF(eleffweight); BranchF(eleffweight_up); BranchF(eleffweight_dn); 
  BranchF(mueffweight); BranchF(mueffweight_up); BranchF(mueffweight_dn); 
  BranchF(btagweight);  BranchVF(csvweights);

  /* trigger */
  BranchF(tri);    BranchF(tri_up);  BranchF(tri_dn);
  BranchO(trig_m); BranchO(trig_m2); BranchO(trig_e); BranchO(trig_mm); BranchO(trig_em); BranchO(trig_ee);

  /* For MatchingForMC() */
  BranchI(matchedS_jidx);
  BranchI(matchedS_idx);     BranchI(matchedS_pdgId);     BranchI(matchedS_status);
  BranchI(matchedS_mom_idx); BranchI(matchedS_mom_pdgId); BranchI(matchedS_mom_status);
  BranchTLV(matchedS_tlv);

  BranchI(matchedB_jidx);
  BranchI(matchedB_idx);     BranchI(matchedB_pdgId);     BranchI(matchedB_status);
  BranchI(matchedB_mom_idx); BranchI(matchedB_mom_pdgId); BranchI(matchedB_mom_status);
  BranchTLV(matchedB_tlv);

  BranchF(Jet_dr_closest_s);    BranchF(Jet_dr_closest_b);
  BranchF(SelJet_dr_closest_s); BranchF(SelJet_dr_closest_b);
  BranchF(GenJet_dr_closest_s); BranchF(GenJet_dr_closest_b);

  BranchO(GenSJet);             BranchO(GenBJet);             BranchO(GenBothJet);             BranchO(RecSJet);             BranchO(RecBJet);             BranchO(RecBothJet); 
  BranchO(GenSJetClosestToLep); BranchO(GenBJetClosestToLep); BranchO(GenBothJetClosestToLep); BranchO(RecSJetClosestToLep); BranchO(RecBJetClosestToLep); BranchO(RecBothJetClosestToLep);  
  BranchO(GenSJetIsHighest);    BranchO(GenBJetIsHighest);    BranchO(GenBothJetIsHighest);    BranchO(RecSJetIsHighest);    BranchO(RecBJetIsHighest);    BranchO(RecBothJetIsHighest); 

  /* For GenAnalysis() */
  BranchVI(GenPart_isGenFrom_vec);    BranchVO(GenPart_isGenFromTop_vec);  BranchVO(GenPart_isGenFromW_vec);   BranchVO(GenPart_isFromKstar_vec);
  BranchVF(GenPart_d_vec);            BranchVF(GenPart_pt_vec);            BranchVF(GenPart_eta_vec);          BranchVF(GenPart_phi_vec);          BranchVF(GenPart_mass_vec);
  BranchVF(GenPart_x_closest_j_vec);  BranchVF(GenPart_dr_closest_j_vec);  BranchVF(GenPart_x_highest_j_vec);  BranchVF(GenPart_dr_highest_j_vec);
  BranchVF(GenPart_x_closest_gj_vec); BranchVF(GenPart_dr_closest_gj_vec); BranchVF(GenPart_x_highest_gj_vec); BranchVF(GenPart_dr_highest_gj_vec);
  BranchVI(GenPart_isClosestPair_xOrder_j_vec);  BranchVI(GenPart_isHighestPair_xOrder_j_vec); 
  BranchVI(GenPart_isClosestPair_xOrder_gj_vec); BranchVI(GenPart_isHighestPair_xOrder_gj_vec);

  /* For RecAnalysis() */
  BranchVI(hadTruth_pdgId_vec);        BranchVI(hadTruth_nMatched_vec);      BranchVI(hadTruth_isFrom_vec);
  BranchVO(hadTruth_isHadFromTop_vec); BranchVO(hadTruth_isHadFromW_vec);    BranchVO(hadTruth_isHadFromS_vec);   BranchVO(hadTruth_isHadFromC_vec);   BranchVO(hadTruth_isHadFromB_vec);
  BranchVF(hadTruth_d_vec);            BranchVF(hadTruth_pt_vec);            BranchVF(hadTruth_eta_vec);          BranchVF(hadTruth_phi_vec);          BranchVF(hadTruth_mass_vec);
  BranchVF(hadTruth_lxy_vec);          BranchVF(hadTruth_lxySig_vec);        BranchVF(hadTruth_l3D_vec);          BranchVF(hadTruth_l3DSig_vec);       BranchVF(hadTruth_legDR_vec);
  BranchVF(hadTruth_angleXY_vec);      BranchVF(hadTruth_angleXYZ_vec);      BranchVF(hadTruth_chi2_vec);         BranchVF(hadTruth_dca_vec);
  BranchVF(hadTruth_dau1_chi2_vec);    BranchVF(hadTruth_dau1_ipsigXY_vec);  BranchVF(hadTruth_dau1_ipsigZ_vec);  BranchVF(hadTruth_dau1_pt_vec);
  BranchVF(hadTruth_dau2_chi2_vec);    BranchVF(hadTruth_dau2_ipsigXY_vec);  BranchVF(hadTruth_dau2_ipsigZ_vec);  BranchVF(hadTruth_dau2_pt_vec);
  BranchVF(hadTruth_x_closest_j_vec);  BranchVF(hadTruth_dr_closest_j_vec);  BranchVF(hadTruth_x_highest_j_vec);  BranchVF(hadTruth_dr_highest_j_vec);
  BranchVF(hadTruth_x_closest_gj_vec); BranchVF(hadTruth_dr_closest_gj_vec); BranchVF(hadTruth_x_highest_gj_vec); BranchVF(hadTruth_dr_highest_gj_vec);
  BranchVI(hadTruth_isClosestPair_xOrder_j_vec);  BranchVI(hadTruth_isHighestPair_xOrder_j_vec);
  BranchVI(hadTruth_isClosestPair_xOrder_gj_vec); BranchVI(hadTruth_isHighestPair_xOrder_gj_vec);

  /* For JetAnalysis() */
  BranchVI(Jet_isCorrectMat);

  /* For CollectVar() */
  BranchF(MET_pt);   BranchF(MET_phi);  BranchF(MET_sumEt);
  BranchI(lep1_pid); BranchI(lep2_pid);
  BranchTLV(lep1);   BranchTLV(lep2);
  BranchTLV(dilep);
}

void vtsAnalyser::ResetBranch() {
  Reset();
  ResetForTMVA();  

  b_nJet = -1; b_nSelJet = -1; b_nSelJetEv = -1; 
  b_passedEvent = false;

  b_hadTruth_nMatched     = -1;    b_hadTruth_nTrueDau = -1;
  b_hadTruth_isHadFromTsb = -1;    b_hadTruth_isHadFromTop = false; b_hadTruth_isHadFromW = false; b_hadTruth_isHadFromS = false; b_hadTruth_isHadFromC = false; b_hadTruth_isHadFromB = false;

  /* For MatchingForMC() */
  m_tqMC.clear();       m_wqMC.clear(); 
  m_qjMapForMC.clear(); m_qgjMapForMC.clear();
  m_recJet.clear();     m_genJet.clear();
  m_closestRecJetForLep1.clear(); m_closestRecJetForLep2.clear(); 
  m_closestGenJetForLep1.clear(); m_closestGenJetForLep2.clear();

  b_matchedS_jidx    = -1;
  b_matchedS_idx     = -1; b_matchedS_pdgId     = -99; b_matchedS_status     = -1;
  b_matchedS_mom_idx = -1; b_matchedS_mom_pdgId = -99; b_matchedS_mom_status = -1;
  b_matchedS_tlv.SetPtEtaPhiM(0,0,0,0);

  b_matchedB_jidx    = -1;
  b_matchedB_idx     = -1; b_matchedB_pdgId     = -99; b_matchedB_status     = -1;
  b_matchedB_mom_idx = -1; b_matchedB_mom_pdgId = -99; b_matchedB_mom_status = -1;
  b_matchedB_tlv.SetPtEtaPhiM(0,0,0,0);

  b_Jet_dr_closest_s    = -1; b_Jet_dr_closest_b    = -1;
  b_SelJet_dr_closest_s = -1; b_SelJet_dr_closest_b = -1;
  b_GenJet_dr_closest_s = -1; b_GenJet_dr_closest_b = -1;
  b_GenSJet             = false; b_GenBJet             = false; b_GenBothJet             = false; b_RecSJet             = false; b_RecBJet             = false; b_RecBothJet             = false;
  b_GenSJetClosestToLep = false; b_GenBJetClosestToLep = false; b_GenBothJetClosestToLep = false; b_RecSJetClosestToLep = false; b_RecBJetClosestToLep = false; b_RecBothJetClosestToLep = false;
  b_GenSJetIsHighest    = false; b_GenBJetIsHighest    = false; b_GenBothJetIsHighest    = false; b_RecSJetIsHighest    = false; b_RecBJetIsHighest    = false; b_RecBothJetIsHighest    = false;

  /* For GenAnalysis() */
  b_GenPart_isGenFrom_vec.clear();    b_GenPart_isGenFromTop_vec.clear();  b_GenPart_isGenFromW_vec.clear();   b_GenPart_isFromKstar_vec.clear();
  b_GenPart_d_vec.clear();            b_GenPart_pt_vec.clear();            b_GenPart_eta_vec.clear();          b_GenPart_phi_vec.clear();          b_GenPart_mass_vec.clear();
  b_GenPart_x_closest_j_vec.clear();  b_GenPart_dr_closest_j_vec.clear();  b_GenPart_x_highest_j_vec.clear();  b_GenPart_dr_highest_j_vec.clear();
  b_GenPart_x_closest_gj_vec.clear(); b_GenPart_dr_closest_gj_vec.clear(); b_GenPart_x_highest_gj_vec.clear(); b_GenPart_dr_highest_gj_vec.clear();
  b_GenPart_isClosestPair_xOrder_j_vec.clear();  b_GenPart_isHighestPair_xOrder_j_vec.clear();
  b_GenPart_isClosestPair_xOrder_gj_vec.clear(); b_GenPart_isHighestPair_xOrder_gj_vec.clear();

  /* For RecAnalysis() */
  b_hadTruth_pdgId_vec.clear();        b_hadTruth_nMatched_vec.clear();      b_hadTruth_isFrom_vec.clear(); 
  b_hadTruth_isHadFromTop_vec.clear(); b_hadTruth_isHadFromW_vec.clear();    b_hadTruth_isHadFromS_vec.clear();   b_hadTruth_isHadFromC_vec.clear();   b_hadTruth_isHadFromB_vec.clear(); 
  b_hadTruth_d_vec.clear();            b_hadTruth_pt_vec.clear();            b_hadTruth_eta_vec.clear();          b_hadTruth_phi_vec.clear();          b_hadTruth_mass_vec.clear();
  b_hadTruth_lxy_vec.clear();          b_hadTruth_lxySig_vec.clear();        b_hadTruth_l3D_vec.clear();          b_hadTruth_l3DSig_vec.clear();       b_hadTruth_legDR_vec.clear();
  b_hadTruth_angleXY_vec.clear();      b_hadTruth_angleXYZ_vec.clear();      b_hadTruth_chi2_vec.clear();         b_hadTruth_dca_vec.clear();
  b_hadTruth_dau1_chi2_vec.clear();    b_hadTruth_dau1_ipsigXY_vec.clear();  b_hadTruth_dau1_ipsigZ_vec.clear();  b_hadTruth_dau1_pt_vec.clear();
  b_hadTruth_dau2_chi2_vec.clear();    b_hadTruth_dau2_ipsigXY_vec.clear();  b_hadTruth_dau2_ipsigZ_vec.clear();  b_hadTruth_dau2_pt_vec.clear();
  b_hadTruth_x_closest_j_vec.clear();  b_hadTruth_dr_closest_j_vec.clear();  b_hadTruth_x_highest_j_vec.clear();  b_hadTruth_dr_highest_j_vec.clear();
  b_hadTruth_x_closest_gj_vec.clear(); b_hadTruth_dr_closest_gj_vec.clear(); b_hadTruth_x_highest_gj_vec.clear(); b_hadTruth_dr_highest_gj_vec.clear();
  b_hadTruth_isClosestPair_xOrder_j_vec.clear();  b_hadTruth_isHighestPair_xOrder_j_vec.clear();
  b_hadTruth_isClosestPair_xOrder_gj_vec.clear(); b_hadTruth_isHighestPair_xOrder_gj_vec.clear();

  /* For JetAnalysis() */
  b_Jet_isCorrectMat.clear();

  /* For CollectVar() */
  b_MET_pt = -1; b_MET_phi = -99; b_MET_sumEt = -1;
}

void vtsAnalyser::MatchingForMC() {
  /* Find s quark from Gen Info. */
  /*
     status is case of pythia 
     top quark statusFlag(status) : 10497(62) 
     s/b quark from top  statusFlag(status) : 22913(23) or 4481(23)
     W boson (t->W->q) statusFlag(status) : 14721(22)
     W boson 1 (t->W1->W2->q) statusFlag(status) : 4481(22)  W boson 2 (t->W1->W2->q)  statusFlag(status) : 10497(52)
     quark from W boson statusFlag(status) : 22913(23) or 4481(23)
     For herwig sample
     top quark statusFlag(status) : 8193(11)
     s/b quark from top           : 4097(11)
  */

  // Finding s quark using GenStatusFlags + pdgId
  for (unsigned int i=0; i<nGenPart; ++i) {
    auto mom_pdgId  = GenPart_pdgId[GenPart_genPartIdxMother[i]];
    auto mom_status = GenPart_status[GenPart_genPartIdxMother[i]];
    if ((abs(GenPart_pdgId[i]) == 3 || abs(GenPart_pdgId[i]) == 5 || abs(GenPart_pdgId[i]) == 4) && ( (GenPart_statusFlags[i] & ( 1 << reco::GenStatusFlags::kIsFirstCopy)) != 0  ) ) {
      if ( (abs(mom_pdgId) == 6 && ( (GenPart_statusFlags[GenPart_genPartIdxMother[i]] & ( 1 << reco::GenStatusFlags::kIsLastCopy)) != 0 ) ) ) { 
        TLorentzVector tlv;
	tlv.SetPtEtaPhiM(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i]);
        m_tqMC.push_back({(int) i, GenPart_pdgId[i], GenPart_status[i], GenPart_genPartIdxMother[i], mom_pdgId, mom_status, tlv});
        cout << " qMC [ " << i << " / " << nGenPart << " ]  >>>>>>>>>>>> got it <<<<<<<<<<< " << endl;
        cout << " qMC [ " << i << " / " << nGenPart << " ]  >>>>> quark is : " << GenPart_pdgId[i] << " ( " << GenPart_status[i] << " ) ||| from : " << mom_pdgId << " ( " << mom_status << " ) " << endl;
      }
      if ( (abs(mom_pdgId) == 24 && (mom_status == 22 || mom_status == 52)) ) { // not completed yet since not used for now
        TLorentzVector tlv;
        tlv.SetPtEtaPhiM(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i]);
        m_wqMC.push_back({(int) i, GenPart_pdgId[i], GenPart_status[i], GenPart_genPartIdxMother[i], mom_pdgId, mom_status, tlv});
        cout << " wqMC [ " << i << " / " << nGenPart << " ]  >>>>>>>>>>>> got it <<<<<<<<<<< " << endl;
        cout << " wqMC [ " << i << " / " << nGenPart << " ]  >>>>> quark is : " << GenPart_pdgId[i] << " ( " << GenPart_status[i] << " ) ||| from : " << mom_pdgId << " ( " << mom_status << " ) " << endl;
      }
    }
  }

  if (m_tqMC.size() < 2) { 
    cout << " >>>>> it isn't a ttbar->qWqW event. save nothing." << endl; 
    if (!m_isMC) {
      cout << " >>>>> it isn't a MC sample ::: fill rec and gen jet info vector with (index, pt, dr(q1,j) = -1, dr(q2,j) = -1, dr(lep1,j) = -1, dr(lep2,j) = -1 <<<<<" << endl;
      for (unsigned int j=0; j<nGenJet; ++j) {
        m_genJet.push_back({(int) j, GenJet_pt[j], -1, -1, -1, -1});
      }
      auto selectedJet = jetSelection();
      if (selectedJet.size() != 0) {
        for (unsigned int ij=0; ij<selectedJet.size();++ij) {
          auto j = selectedJet[ij].GetFirstMother();
          m_recJet.push_back({j, Jet_pt[j], -1, -1, -1, -1});
        }
      }
    }
    return; 
  }
  /* Select quarks originated from t->qW */
  genInfo tq1; genInfo tq2; // tq1 : s-quark , tq2 : b-quark
  if      (abs(m_tqMC[0].pdgId) == 3) { tq1 = m_tqMC[0]; tq2 = m_tqMC[1]; } 
  else if (abs(m_tqMC[0].pdgId) == 5) { tq1 = m_tqMC[1]; tq2 = m_tqMC[0]; }
  TLorentzVector tq1_tlv = tq1.tlv; 
  TLorentzVector tq2_tlv = tq2.tlv;

  /* Select quarks from W->qq (not yet used) */
  int wq1;
  TLorentzVector wq1_tlv;
  int wq2;
  TLorentzVector wq2_tlv;
  if (m_wqMC.size() != 0) {
    cout << "number of quark from W-boson : " << m_wqMC.size() << endl;
    wq1 = m_wqMC[0].idx;
    wq1_tlv.SetPtEtaPhiM(GenPart_pt[wq1], GenPart_eta[wq1], GenPart_phi[wq1], GenPart_mass[wq1]);
    if (m_wqMC.size() == 2) {
      wq2 = m_wqMC[1].idx;
      wq2_tlv.SetPtEtaPhiM(GenPart_pt[wq2], GenPart_eta[wq2], GenPart_phi[wq2], GenPart_mass[wq2]);
    }
  }

  /* Gen Particle & Gen Jet matching */
  /* tydef jetInfo => {jet idx, jet pt, DeltaR(s,jet), DeltaR(b,jet), DeltaR(lep1,jet), DeltaR(lep2,jet) } */
  int noverlap = 0;
  for (unsigned int j=0; j<nGenJet; ++j) {
    TLorentzVector jet_tlv;
    jet_tlv.SetPtEtaPhiM(GenJet_pt[j], GenJet_eta[j], GenJet_phi[j], GenJet_mass[j]);
    auto drs = tq1_tlv.DeltaR(jet_tlv); auto drb = tq2_tlv.DeltaR(jet_tlv);
    auto drl1 = b_lep1.DeltaR(jet_tlv); auto drl2 = b_lep2.DeltaR(jet_tlv);
    if (drl1 < 0.4 || drl2 < 0.4) { ++noverlap; continue; }
    m_genJet.push_back({(int) j, GenJet_pt[j], drs, drb, drl1, drl2});
  }
  if (noverlap > 2) std::cout << " >>>> BAD GENJET OVERLAP REMOVAL FOR EVENT <<<< " << std::endl;
    
  /* Find closest gen jet to l+ or l- */
  sort(m_genJet.begin(), m_genJet.end(), [] (jetInfo a, jetInfo b) { return (a.drl1j < b.drl1j); } );
  m_closestGenJetForLep1[m_genJet[0].idx] = b_lep1_pid;
  sort(m_genJet.begin(), m_genJet.end(), [] (jetInfo a, jetInfo b) { return (a.drl2j < b.drl2j); } );
  m_closestGenJetForLep2[m_genJet[0].idx] = b_lep2_pid;
  /* Find closest gen jet to s or b */
  sort(m_genJet.begin(), m_genJet.end(), [] (jetInfo a, jetInfo b) { return (a.drsj < b.drsj); } );
  if (m_genJet[0].drsj <= m_jetConeSize) m_qgjMapForMC[m_genJet[0].idx] = tq1.pdgId; // if jet is inside dRCut == 0.4, then label the tag about s/b
  sort(m_genJet.begin(), m_genJet.end(), [] (jetInfo a, jetInfo b) { return (a.drbj < b.drbj); } );
  if (m_genJet[0].drbj <= m_jetConeSize) m_qgjMapForMC[m_genJet[0].idx] = tq2.pdgId;

  /* Gen Partcle & Selected Reco Jet Matching */
  auto selectedJet = jetSelection();
  if (selectedJet.size() != 0) {
    for (unsigned int ij=0; ij<selectedJet.size();++ij) {
      auto j = selectedJet[ij].GetFirstMother();
      TLorentzVector jet_tlv;
      jet_tlv.SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
      auto drs = tq1_tlv.DeltaR(jet_tlv); auto drb = tq2_tlv.DeltaR(jet_tlv);
      auto drl1 = b_lep1.DeltaR(jet_tlv); auto drl2 = b_lep2.DeltaR(jet_tlv);
      m_recJet.push_back({j, Jet_pt[j], drs, drb, drl1, drl2});
    }
    /* Find closest sel jet to l+ or l- */
    sort(m_recJet.begin(), m_recJet.end(), [] (jetInfo a, jetInfo b) { return (a.drl1j < b.drl1j); } ); // order by dR(jet, lepton1)
    m_closestRecJetForLep1[m_recJet[0].idx] = b_lep1_pid; // save std::map (jet idx as key, lepton1 pdgId as value)
    sort(m_recJet.begin(), m_recJet.end(), [] (jetInfo a, jetInfo b) { return (a.drl2j < b.drl2j); } ); // order by dR(jet, lepton2)
    m_closestRecJetForLep2[m_recJet[0].idx] = b_lep2_pid; // save std::map (jet idx as key, lepton2 pdgId as value)

    /* Find closest sel jet to s or b */ 
    sort(m_recJet.begin(), m_recJet.end(), [] (jetInfo a, jetInfo b) { return (a.drsj < b.drsj); } ); // order by dR(jet, gen s-quark)
    if (m_recJet[0].drsj <= m_jetConeSize) { 
      m_qjMapForMC[m_recJet[0].idx] = tq1.pdgId; // if jet is inside dRCut == 0.4, then label the tag about s
      b_matchedS_jidx               = m_recJet[0].idx;
      b_matchedS_idx                = tq1.idx;
      b_matchedS_pdgId              = tq1.pdgId;
      b_matchedS_status             = tq1.status;
      b_matchedS_mom_idx            = tq1.mom_idx;
      b_matchedS_mom_pdgId          = tq1.mom_pdgId;
      b_matchedS_mom_status         = tq1.mom_status;
      b_matchedS_tlv                = tq1.tlv;
    }
    sort(m_recJet.begin(), m_recJet.end(), [] (jetInfo a, jetInfo b) { return (a.drbj < b.drbj); } ); // order by dR(jet, gen b-quark)
    if (m_recJet[0].drbj <= m_jetConeSize) {
      m_qjMapForMC[m_recJet[0].idx] = tq2.pdgId; // if jet is inside dRCut == 0.4, then label the tag about b
      b_matchedB_jidx               = m_recJet[0].idx;
      b_matchedB_idx                = tq2.idx;
      b_matchedB_pdgId              = tq2.pdgId;
      b_matchedB_status             = tq2.status;
      b_matchedB_mom_idx            = tq2.mom_idx;
      b_matchedB_mom_pdgId          = tq2.mom_pdgId;
      b_matchedB_mom_status         = tq2.mom_status;
      b_matchedB_tlv                = tq2.tlv;
    }
  } else cout << "Size of selectedJets is zero" << endl;

  /* Save closest jet dr for s/b, Rec(Gen)S(B/Both)Jet, Rec(Gen)S(B/Both)Jet which is closest to l+ or l- and is in highest 2 pT jets */
  /* GenJet */
  int GenSJidx = -1, GenBJidx = -1;
  if (m_genJet.size() != 0) {
    sort(m_genJet.begin(), m_genJet.end(), [] (jetInfo a, jetInfo b) { return (a.drsj < b.drsj); } );
    b_GenJet_dr_closest_s = m_genJet[0].drsj;
    if (b_GenJet_dr_closest_s <= m_jetConeSize) {
      GenSJidx = m_genJet[0].idx;
      b_GenSJet = true;
      if (m_closestGenJetForLep1[m_genJet[0].idx] != 0 || m_closestGenJetForLep2[m_genJet[0].idx] != 0) b_GenSJetClosestToLep = true;
    }
    sort(m_genJet.begin(), m_genJet.end(), [] (jetInfo a, jetInfo b) { return (a.drbj < b.drbj); } );
    b_GenJet_dr_closest_b = m_genJet[0].drbj;
    if (b_GenJet_dr_closest_b <= m_jetConeSize) {
      GenBJidx = m_genJet[0].idx;
      b_GenBJet = true;
      if (m_closestGenJetForLep1[m_genJet[0].idx] != 0 || m_closestGenJetForLep2[m_genJet[0].idx] != 0) b_GenBJetClosestToLep = true;
    }
    sort(m_genJet.begin(), m_genJet.end(), [] (jetInfo a, jetInfo b) { return (a.pt > b.pt); } );
    if ((int)m_genJet[0].idx == GenSJidx || (int)m_genJet[1].idx == GenSJidx) b_GenSJetIsHighest = true;
    if ((int)m_genJet[0].idx == GenBJidx || (int)m_genJet[1].idx == GenBJidx) b_GenBJetIsHighest = true;
    if (b_GenSJet && b_GenBJet) {
      b_GenBothJet = true;
      if (b_GenSJetIsHighest && b_GenBJetIsHighest) b_GenBothJetIsHighest = true;
      if (b_GenSJetClosestToLep && b_GenBJetClosestToLep) b_GenBothJetClosestToLep = true;
    }
  }
  /* SelectedJet */
  int RecSJidx = -1, RecBJidx = -1;
  if (m_recJet.size() != 0) {
    sort(m_recJet.begin(), m_recJet.end(), [] (jetInfo a, jetInfo b) { return (a.drsj < b.drsj); } ); 
    b_SelJet_dr_closest_s = m_recJet[0].drsj;
    if (b_SelJet_dr_closest_s <= m_jetConeSize) {
      RecSJidx = m_recJet[0].idx;
      b_RecSJet = true;
      if (m_closestRecJetForLep1[m_recJet[0].idx] != 0 || m_closestRecJetForLep2[m_recJet[0].idx] != 0) b_RecSJetClosestToLep = true;
    }
    sort(m_recJet.begin(), m_recJet.end(), [] (jetInfo a, jetInfo b) { return (a.drbj < b.drbj); } );
    b_SelJet_dr_closest_b = m_recJet[0].drbj;
    if (b_SelJet_dr_closest_b <= m_jetConeSize) {
      RecBJidx = m_recJet[0].idx;
      b_RecBJet = true;
      if (m_closestRecJetForLep1[m_recJet[0].idx] != 0 || m_closestRecJetForLep2[m_recJet[0].idx] != 0) b_RecBJetClosestToLep = true;
    }
    sort(m_recJet.begin(), m_recJet.end(), [] (jetInfo a, jetInfo b) { return (a.pt > b.pt); } );
    if ((int)m_recJet[0].idx == RecSJidx || (int)m_recJet[1].idx == RecSJidx) b_RecSJetIsHighest = true;
    if ((int)m_recJet[0].idx == RecBJidx || (int)m_recJet[1].idx == RecBJidx) b_RecBJetIsHighest = true;
    if (b_RecSJet && b_RecBJet) {
      b_RecBothJet = true;
      if (b_RecSJetIsHighest && b_RecBJetIsHighest) b_RecBothJetIsHighest = true;
      if (b_RecSJetClosestToLep && b_RecBJetClosestToLep) b_RecBothJetClosestToLep = true;
    }
  }
}

bool vtsAnalyser::isGenFrom(int count, int idx, int & isFrom, bool & isFromTop, bool & isFromW, bool & isFromKstar) {
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

struct jetks { unsigned int idx; double dr; double x; };

void vtsAnalyser::GenAnalysis() {
  std::vector<jetks> recPair1; std::vector<jetks> genPair1;
  std::vector<jetks> recPair2; std::vector<jetks> genPair2;
  int nGenKS = 0;
  b_GenPart_dr_closest_j_vec.resize(nGenPart,-1);
  b_GenPart_x_closest_j_vec.resize(nGenPart,-1);
  b_GenPart_dr_closest_gj_vec.resize(nGenPart,-1);
  b_GenPart_x_closest_gj_vec.resize(nGenPart,-1);
  b_GenPart_dr_highest_j_vec.resize(nGenPart,-1);
  b_GenPart_x_highest_j_vec.resize(nGenPart,-1);
  b_GenPart_dr_highest_gj_vec.resize(nGenPart,-1);
  b_GenPart_x_highest_gj_vec.resize(nGenPart,-1);
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

    if (nJet !=0) {
      recPair1.resize(nJet, {0,99,-1}); recPair2.resize(nJet, {0,99,-1});
      for (unsigned int j=0; j<nJet; ++j) { // Loop for all of recoJet
        TLorentzVector jet_tlv;
        jet_tlv.SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
        auto dr = jet_tlv.DeltaR(gen_tlv);      auto x = GenPart_pt[i]/Jet_pt[j];
        auto drMin1 = recPair1[j].dr; auto xMax1 = recPair1[j].x;
        auto drMin2 = recPair2[j].dr; auto xMax2 = recPair2[j].x;
        if (isFrom != m_qjMapForMC[j]) continue;
        if (dr < drMin1) recPair1[j] = {i,dr,x};
        else if (dr == drMin1 && x > xMax1) recPair1[j] = {i,dr,x};
        if (x > xMax2) recPair2[j] = {i,dr,x};
        else if (x == xMax2 && dr < drMin2) recPair2[j] = {i,dr,x};
      }
    }
    if (nGenJet !=0) {
      genPair1.resize(nGenJet, {0,99,-1}); genPair2.resize(nGenJet, {0,99,-1});
      for (unsigned int j=0; j<nGenJet; ++j) { // Loop for all of recoJet
        TLorentzVector jet_tlv;
        jet_tlv.SetPtEtaPhiM(GenJet_pt[j], GenJet_eta[j], GenJet_phi[j], GenJet_mass[j]);
        auto dr = jet_tlv.DeltaR(gen_tlv); auto x = GenPart_pt[i]/GenJet_pt[j];
        auto drMin1 = genPair1[j].dr; auto xMax1 = genPair1[j].x;
        auto drMin2 = genPair2[j].dr; auto xMax2 = genPair2[j].x;
        if (isFrom != m_qgjMapForMC[j]) continue;
        if (dr < drMin1) genPair1[j] = {i,dr,x};
        else if (dr == drMin1 && x > xMax1) genPair1[j] = {i,dr,x};
        if (x > xMax2) genPair2[j] = {i,dr,x};
        else if (x == xMax2 && dr < drMin2) genPair2[j] = {i,dr,x};
      }
      b_GenPart_isClosestPair_xOrder_j_vec.push_back(i);  b_GenPart_isHighestPair_xOrder_j_vec.push_back(i);
      b_GenPart_isClosestPair_xOrder_gj_vec.push_back(i); b_GenPart_isHighestPair_xOrder_gj_vec.push_back(i);
    }
  }
  if (nGenKS !=0) {
    /*  Save x and dr for the most closest(highest) KS for each jet */
    if (nJet !=0) {
      for (unsigned int j=0; j<nJet; ++j) {
        int cRecIdx = recPair1[j].idx; int hRecIdx = recPair2[j].idx;
        b_GenPart_dr_closest_j_vec[cRecIdx] = recPair1[j].dr;
        b_GenPart_x_closest_j_vec[cRecIdx] = recPair1[j].x;
        b_GenPart_dr_highest_j_vec[hRecIdx] = recPair2[j].dr;
        b_GenPart_x_highest_j_vec[hRecIdx] = recPair2[j].x;
      }
      /* Give flag ( == index) for the most closest(highest) KS-jet pair per event by x ordering */
      std::sort(recPair1.begin(), recPair1.end(), [] (jetks  a, jetks b) { return (a.x > b.x); } );
      int hRecIdx1 = recPair1[0].idx;
      std::replace_if(b_GenPart_isClosestPair_xOrder_j_vec.begin(),  b_GenPart_isClosestPair_xOrder_j_vec.end(),  [&] (int a) { return (a != hRecIdx1); }, -1); 
      std::sort(recPair2.begin(), recPair2.end(), [] (jetks  a, jetks b) { return (a.x > b.x); } );
      int hRecIdx2 = recPair2[0].idx;
      std::replace_if(b_GenPart_isHighestPair_xOrder_j_vec.begin(),  b_GenPart_isHighestPair_xOrder_j_vec.end(),  [&] (int a) { return (a != hRecIdx2); }, -1);
    }
    if (nGenJet !=0) {
      for (unsigned int j=0; j<nGenJet; ++j) {
        int cGenIdx = genPair1[j].idx;int hGenIdx = genPair2[j].idx;
        b_GenPart_dr_closest_gj_vec[cGenIdx] = genPair1[j].dr;
        b_GenPart_x_closest_gj_vec[cGenIdx] = genPair1[j].x;
        b_GenPart_dr_highest_gj_vec[hGenIdx] = genPair2[j].dr;
        b_GenPart_x_highest_gj_vec[hGenIdx] = genPair2[j].x;
      }
      /* Give flag ( == index) for the most closest(highest) KS-jet pair per event by x ordering */
      std::sort(genPair1.begin(), genPair1.end(), [] (jetks  a, jetks b) { return (a.x > b.x); } );
      int hGenIdx1 = genPair1[0].idx;
      std::replace_if(b_GenPart_isClosestPair_xOrder_gj_vec.begin(), b_GenPart_isClosestPair_xOrder_gj_vec.end(), [&] (int a) { return (a != hGenIdx1); }, -1);
      std::sort(genPair2.begin(), genPair2.end(), [] (jetks  a, jetks b) { return (a.x > b.x); } );
      int hGenIdx2 = genPair2[0].idx;
      std::replace_if(b_GenPart_isHighestPair_xOrder_gj_vec.begin(), b_GenPart_isHighestPair_xOrder_gj_vec.end(), [&] (int a) { return (a != hGenIdx2); }, -1);

    }
  }
}

void vtsAnalyser::RecAnalysis() {
  b_had_start = m_hadtrForTMVA->GetEntries();

  std::vector<jetks> recPair1; std::vector<jetks> genPair1;
  std::vector<jetks> recPair2; std::vector<jetks> genPair2;
  int nRecKS = 0;
  b_hadTruth_dr_closest_j_vec.resize(nhad,-1);
  b_hadTruth_x_closest_j_vec.resize(nhad,-1);
  b_hadTruth_dr_closest_gj_vec.resize(nhad,-1);
  b_hadTruth_x_closest_gj_vec.resize(nhad,-1);
  b_hadTruth_dr_highest_j_vec.resize(nhad,-1);
  b_hadTruth_x_highest_j_vec.resize(nhad,-1);
  b_hadTruth_dr_highest_gj_vec.resize(nhad,-1);
  b_hadTruth_x_highest_gj_vec.resize(nhad,-1);

  for (unsigned int i=0; i<nhad; ++i) {
    if (had_pdgId[i] != 310) continue;
    ++nRecKS;
    if (!m_isGenericMC) {
      b_hadTruth_isHadFromTop_vec.push_back(hadTruth_isHadFromTop[i]);
      b_hadTruth_isHadFromW_vec.push_back(hadTruth_isHadFromW[i]);
      b_hadTruth_isHadFromS_vec.push_back(hadTruth_isHadFromS[i]);
      b_hadTruth_isHadFromC_vec.push_back(hadTruth_isHadFromC[i]);
      b_hadTruth_isHadFromB_vec.push_back(hadTruth_isHadFromB[i]);
      b_hadTruth_isFrom_vec.push_back(hadTruth_isHadFromTsb[i]);
      b_hadTruth_nMatched_vec.push_back(hadTruth_nMatched[i]);
    }
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
    b_hadTruth_l3DSig_vec.push_back(had_l3D[i] / had_l3DErr[i]);
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

    FillHadTreeForTMVA();

    TLorentzVector had_tlv;
    had_tlv.SetPtEtaPhiM(had_pt[i], had_eta[i], had_phi[i], had_mass[i]);
    recPair1.resize(nJet, {0,99,-1}); recPair2.resize(nJet, {0,99,-1});
    if (nJet !=0) {
      for (unsigned int j=0; j<nJet; ++j) { // Loop for all of recoJet
        TLorentzVector jet_tlv;
        jet_tlv.SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
        auto dr = jet_tlv.DeltaR(had_tlv); auto x = had_pt[i]/Jet_pt[j];
        auto drMin1 = recPair1[j].dr; auto xMax1 = recPair1[j].x;
        auto drMin2 = recPair2[j].dr; auto xMax2 = recPair2[j].x;
        if (hadTruth_isHadFromTsb[i] != m_qjMapForMC[j]) continue;
        if (dr < drMin1) recPair1[j] = {i,dr,x};
        else if (dr == drMin1 && x > xMax1) recPair1[j] = {i,dr,x};
        if (x > xMax2) recPair2[j] = {i,dr,x};
        else if (x == xMax2 && dr < drMin2) recPair2[j] = {i,dr,x};
      }
    }
    genPair1.resize(nGenJet, {0,99,-1}); genPair2.resize(nGenJet, {0,99,-1});
    if (nGenJet !=0) {
      for (unsigned int j=0; j<nGenJet; ++j) { // Loop for all of recoJet
        TLorentzVector jet_tlv;
        jet_tlv.SetPtEtaPhiM(GenJet_pt[j], GenJet_eta[j], GenJet_phi[j], GenJet_mass[j]);
        auto dr = jet_tlv.DeltaR(had_tlv); auto x = had_pt[i]/GenJet_pt[j];
        auto drMin1 = genPair1[j].dr; auto xMax1 = genPair1[j].x;
        auto drMin2 = genPair2[j].dr; auto xMax2 = genPair2[j].x;
        if (hadTruth_isHadFromTsb[i] != m_qgjMapForMC[j]) continue;
        if (dr < drMin1) genPair1[j] = {i,dr,x};
        else if (dr == drMin1 && x > xMax1) genPair1[j] = {i,dr,x};
        if (x > xMax2) genPair2[j] = {i,dr,x};
        else if (x == xMax2 && dr < drMin2) genPair2[j] = {i,dr,x};
      }
    }
    b_hadTruth_isClosestPair_xOrder_j_vec.push_back(i);  b_hadTruth_isHighestPair_xOrder_j_vec.push_back(i);
    b_hadTruth_isClosestPair_xOrder_gj_vec.push_back(i); b_hadTruth_isHighestPair_xOrder_gj_vec.push_back(i);
  }
  
  if (nRecKS !=0) {
    /*  Save x and dr for the most closest(highest) KS for each jet */
    if (nJet !=0) {
      for (unsigned int j=0; j<nJet; ++j) {
        int cRecIdx = recPair1[j].idx; int hRecIdx = recPair2[j].idx;
        b_hadTruth_dr_closest_j_vec[cRecIdx] = recPair1[j].dr;
        b_hadTruth_x_closest_j_vec[cRecIdx] = recPair1[j].x;
        b_hadTruth_dr_highest_j_vec[hRecIdx] = recPair2[j].dr;
        b_hadTruth_x_highest_j_vec[hRecIdx] = recPair2[j].x;
      }
      /* Give flag ( == index) for the most closest(highest) KS-jet pair per event by x ordering */
      auto xorder = [] (jetks  a, jetks b) { return (a.x > b.x); };
      int hRecIdx1 = std::max_element(recPair1.begin(), recPair1.end(), xorder)->idx;
      std::replace_if(b_hadTruth_isClosestPair_xOrder_j_vec.begin(),  b_hadTruth_isClosestPair_xOrder_j_vec.end(),  [&] (int a) { return (a != hRecIdx1); }, -1);
      int hRecIdx2 = std::max_element(recPair2.begin(), recPair2.end(), xorder)->idx;
      std::replace_if(b_hadTruth_isHighestPair_xOrder_j_vec.begin(),  b_hadTruth_isHighestPair_xOrder_j_vec.end(),  [&] (int a) { return (a != hRecIdx2); }, -1);
    }
    if (nGenJet !=0) {
      for (unsigned int j=0; j<nGenJet; ++j) {
        int cGenIdx = genPair1[j].idx;int hGenIdx = genPair2[j].idx;
        b_hadTruth_dr_closest_gj_vec[cGenIdx] = genPair1[j].dr;
        b_hadTruth_x_closest_gj_vec[cGenIdx] = genPair1[j].x;
        b_hadTruth_dr_highest_gj_vec[hGenIdx] = genPair2[j].dr;
        b_hadTruth_x_highest_gj_vec[hGenIdx] = genPair2[j].x;
      }
      /* Give flag ( == index) for the most closest(highest) KS-jet pair per event by x ordering */
      auto xorder = [] (jetks  a, jetks b) { return (a.x > b.x); };
      int hGenIdx1 = std::max_element(genPair1.begin(), genPair1.end(), xorder)->idx;
      std::replace_if(b_hadTruth_isClosestPair_xOrder_gj_vec.begin(), b_hadTruth_isClosestPair_xOrder_gj_vec.end(), [&] (int a) { return (a != hGenIdx1); }, -1);
      int hGenIdx2 = std::max_element(genPair2.begin(), genPair2.end(), xorder)->idx;
      std::replace_if(b_hadTruth_isHighestPair_xOrder_gj_vec.begin(), b_hadTruth_isHighestPair_xOrder_gj_vec.end(), [&] (int a) { return (a != hGenIdx2); }, -1);
    }
  }
  b_had_end = m_hadtrForTMVA->GetEntries();
}

void vtsAnalyser::JetAnalysis() {
  /* GenJet */
  

  /* SelectedJet */


  for (unsigned int i=0; i<nJet; ++i) {
    if (Jet_pt[i] < 30.) continue;
    if (abs(Jet_eta[i]) > 2.4) continue;
    if (Jet_jetId[i] < 1) continue;
    /* quark - jet matching check */
    if (abs(Jet_partonFlavour[i]) == 3 || abs(Jet_partonFlavour[i]) == 5) {
      if (Jet_partonFlavour[i] == m_qjMapForMC[i]) b_Jet_isCorrectMat.push_back(Jet_partonFlavour[i]);
      else b_Jet_isCorrectMat.push_back(0);
    }
  }
}

void vtsAnalyser::CollectVar() {
  b_MET_pt = MET_pt;
  b_MET_phi = MET_phi;
  b_MET_sumEt = MET_sumEt;
}

void vtsAnalyser::ResetForTMVA() {
  b_Jet_bdt_score   = -99;
  b_JKS_bdt_score   = -99;

  /* temp initialization of array */
  for (unsigned int i = 0; i < m_jetDauArrSize; ++i) {
    b_dau_pt[i]     = -9.;
    b_dau_eta[i]    = -9.;
    b_dau_phi[i]    = -9.;
    b_dau_charge[i] = -9;
  }

  b_KS_idx          = -99;   b_KS_nMatched     = -99;   b_KS_isFrom       = -99; 
  b_KS_isHadFromTop = false; b_KS_isHadFromW   = false; b_KS_isHadFromS   = false; b_KS_isHadFromC   = false; b_KS_isHadFromB   = false; 
  b_KS_d            = -99;   b_KS_pt           = -99;   b_KS_eta          = -99;   b_KS_phi          = -99;   b_KS_mass         = -99; 
  b_KS_lxy          = -99;   b_KS_lxySig       = -99;   b_KS_l3D          = -99;   b_KS_l3DSig       = -99;   b_KS_legDR        = -99; 
  b_KS_angleXY      = -99;   b_KS_angleXYZ     = -99;   b_KS_chi2         = -99;   b_KS_dca          = -99; 
  b_KS_dau1_chi2    = -99;   b_KS_dau1_ipsigXY = -99;   b_KS_dau1_ipsigZ  = -99;   b_KS_dau1_pt      = -99;   b_KS_dau2_chi2    = -99; 
  b_KS_dau2_ipsigXY = -99;   b_KS_dau2_ipsigZ  = -99;   b_KS_dau2_pt      = -99; 
  b_KS_dr           = -99;   b_KS_x            = -99;   b_KS_best_bdt     = -99; 
}

void vtsAnalyser::FillJetTreeForTMVA() {
  auto selectedJet = jetSelection();
  b_jet_start =  m_jettrForTMVA->GetEntries();
  if (selectedJet.size() != 0) {
    if (m_recJet.size() == 0) return;
    sort(m_recJet.begin(), m_recJet.end(), [] (jetInfo a, jetInfo b) { return (a.pt > b.pt); } ); // pT ordering
    auto highest_first_idx = m_recJet[0].idx; auto highest_second_idx = m_recJet[1].idx;
    sort(m_recJet.begin(), m_recJet.end(), [] (jetInfo a, jetInfo b) { return (a.drsj < b.drsj); } ); // dR(s,jet) ordering
    auto closest_s_idx = m_recJet[0].idx; auto closest_s_dr = m_recJet[0].drsj;
    sort(m_recJet.begin(), m_recJet.end(), [] (jetInfo a, jetInfo b) { return (a.drbj < b.drbj); } ); // dR(b,jet) ordering
    auto closest_b_idx = m_recJet[0].idx; auto closest_b_dr = m_recJet[0].drbj;
    sort(m_recJet.begin(), m_recJet.end(), [] (jetInfo a, jetInfo b) { return (a.drl1j < b.drl1j); } ); // dR(lep1,jet) ordering
    auto closest_lep1_idx = m_recJet[0].idx;
    sort(m_recJet.begin(), m_recJet.end(), [] (jetInfo a, jetInfo b) { return (a.drl2j < b.drl2j); } ); // dR(lep2,jet) ordering
    auto closest_lep2_idx = m_recJet[0].idx;

    for (unsigned int ij=0; ij<selectedJet.size();++ij) {
      ResetForTMVA();
      b_isSJet = 0; b_isBJet = 0; b_isHighest = 0; b_isClosestToLep = 0;

      auto j = selectedJet[ij].GetFirstMother();
      TLorentzVector jet_tlv; jet_tlv.SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);

      if ((j == (int) closest_s_idx) && (fabs(closest_s_dr) <= m_jetConeSize)) b_isSJet = 1;
      if ((j == (int) closest_b_idx) && (fabs(closest_b_dr) <= m_jetConeSize)) b_isBJet = 1;
      if ((j == (int) highest_first_idx) || (j == (int) highest_second_idx))   b_isHighest = 1;
      if ((j == (int) closest_lep1_idx)  || (j == (int) closest_lep2_idx))     b_isClosestToLep = 1;

      /* distinguish fake isClosestToLep */
      if (j == (int) closest_lep1_idx) {
        if ( m_recJet[0].drl1j == -1 ) b_isClosestToLep = -1;
      }
      if (j == (int) closest_lep2_idx) {
        if ( m_recJet[0].drl2j == -1 ) b_isClosestToLep = -1;
      }
      TLorentzVector tlv_j; 
      tlv_j.SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);      
    
      cout << " [ " << ij << " / " << selectedJet.size() << " ] >>>>>> [ " << setw(10) << tlv_j.Pt() << " " << setw(10) << tlv_j.Eta() << " " << setw(10) << tlv_j.Phi() << " " << setw(10) << tlv_j.M() << " ] >>>>>> ( " << setw(3) << j << " " << setw(3) << m_qjMapForMC[j] << " " << setw(3) << b_isSJet << " " << setw(3) << b_isBJet << " " << setw(3) << b_isHighest << " " << setw(3) << b_isClosestToLep << " ) " << endl;

      b_cmult = (float)jetID_cmult[j]; b_nmult = (float)jetID_nmult[j];
      b_pt = Jet_pt[j]; b_eta = Jet_eta[j]; b_phi = Jet_phi[j]; b_mass = Jet_mass[j];
      b_c_x1 = jetID_cpt1[j]/Jet_pt[j]; b_c_x2 = jetID_cpt2[j]/Jet_pt[j]; b_c_x3 = jetID_cpt3[j]/Jet_pt[j];
      b_n_x1 = jetID_npt1[j]/Jet_pt[j]; b_n_x2 = jetID_npt2[j]/Jet_pt[j]; b_n_x3 = jetID_npt3[j]/Jet_pt[j];
      b_axis1 = jetID_axis1[j]; b_axis2 = jetID_axis2[j]; b_ptD = jetID_ptD[j];
      b_area = Jet_area[j]; b_CSVV2 = Jet_btagCSVV2[j];

      /* Save jet daughter information */
      int ia = 0;
      for (auto didx = jetID_dauIdx1[j]; didx < jetID_dauIdx2[j]; ++didx) {
        b_dau_pt[ia]     = jetDau_pt[didx];
        b_dau_eta[ia]    = jetDau_eta[didx];
        b_dau_phi[ia]    = jetDau_phi[didx];
        b_dau_charge[ia] = jetDau_charge[didx];
        ia += 1;
      }

      std::pair<int, float> best_BDT = {-1, -99};
      
      if (b_had_start != -1 && b_had_end != -1) {
        for (auto ih=b_had_start; ih<b_had_end; ++ih) {
          m_hadtrForTMVA->GetEntry(ih);
          TLorentzVector had_tlv; had_tlv.SetPtEtaPhiM(b_Rec_pt, b_Rec_eta, b_Rec_phi, b_Rec_mass);
          if (jet_tlv.DeltaR(had_tlv) <= m_jetConeSize) {
            if (best_BDT.second < b_Rec_bdt_score) {best_BDT = {ih, b_Rec_bdt_score}; b_KS_dr = jet_tlv.DeltaR(had_tlv);} 
          }
        }
        int idx = best_BDT.first;
        if (idx != -1) {
          m_hadtrForTMVA->GetEntry(idx);
          b_KS_idx          = idx;
          b_KS_nMatched     = b_Rec_nMatched;
          b_KS_isFrom       = b_Rec_isFrom;
          b_KS_isHadFromTop = b_Rec_isHadFromTop;  
          b_KS_isHadFromW   = b_Rec_isHadFromW; 
          b_KS_isHadFromS   = b_Rec_isHadFromS;   
          b_KS_isHadFromC   = b_Rec_isHadFromC;   
          b_KS_isHadFromB   = b_Rec_isHadFromB;   
          b_KS_d            = b_Rec_d;      
          b_KS_pt           = b_Rec_pt;        
          b_KS_eta          = b_Rec_eta;         
          b_KS_phi          = b_Rec_phi;         
          b_KS_mass         = b_Rec_mass;        
          b_KS_lxy          = b_Rec_lxy;       
          b_KS_lxySig       = b_Rec_lxySig;      
          b_KS_l3D          = b_Rec_l3D;        
          b_KS_l3DSig       = b_Rec_l3DSig;       
          b_KS_legDR        = b_Rec_legDR;      
          b_KS_angleXY      = b_Rec_angleXY;            
          b_KS_angleXYZ     = b_Rec_angleXYZ;              
          b_KS_chi2         = b_Rec_chi2;        
          b_KS_dca          = b_Rec_dca;              
          b_KS_dau1_chi2    = b_Rec_dau1_chi2;            
          b_KS_dau1_ipsigXY = b_Rec_dau1_ipsigXY;              
          b_KS_dau1_ipsigZ  = b_Rec_dau1_ipsigZ;           
          b_KS_dau1_pt      = b_Rec_dau1_pt;            
          b_KS_dau2_chi2    = b_Rec_dau2_chi2;               
          b_KS_dau2_ipsigXY = b_Rec_dau2_ipsigXY;                
          b_KS_dau2_ipsigZ  = b_Rec_dau2_ipsigZ;                 
          b_KS_dau2_pt      = b_Rec_dau2_pt;
          b_KS_x            = b_Rec_pt/Jet_pt[j];
          b_KS_best_bdt     = b_Rec_bdt_score;          
        }
      }
      b_Jet_bdt_score = m_jetReader->EvaluateMVA("Jet_BDT_highest");
      b_JKS_bdt_score = m_jksReader->EvaluateMVA("JKS_BDT_highest");
      m_jettrForTMVA->Fill();
    }
  } else cout << ">>>> Size of selectedJets is zero <<<< " << endl;
  b_jet_end = m_jettrForTMVA->GetEntries();
}

void vtsAnalyser::FillHadTreeForTMVA() {
  if (!m_isGenericMC) {
    b_Rec_nMatched = b_hadTruth_nMatched_vec.back();
    b_Rec_isFrom = b_hadTruth_isFrom_vec.back();
    b_Rec_isHadFromTop = b_hadTruth_isHadFromTop_vec.back();
    b_Rec_isHadFromW = b_hadTruth_isHadFromW_vec.back();
    b_Rec_isHadFromS = b_hadTruth_isHadFromS_vec.back();
    b_Rec_isHadFromC = b_hadTruth_isHadFromC_vec.back();
    b_Rec_isHadFromB = b_hadTruth_isHadFromB_vec.back();
  }
  b_Rec_pdgId = b_hadTruth_pdgId_vec.back();
  b_Rec_d = b_hadTruth_d_vec.back();
  b_Rec_pt = b_hadTruth_pt_vec.back();
  b_Rec_eta = b_hadTruth_eta_vec.back();
  b_Rec_phi = b_hadTruth_phi_vec.back();
  b_Rec_mass = b_hadTruth_mass_vec.back();
  b_Rec_lxy = b_hadTruth_lxy_vec.back();
  b_Rec_lxySig = b_hadTruth_lxySig_vec.back();
  b_Rec_l3D = b_hadTruth_l3D_vec.back();
  b_Rec_l3DSig = b_hadTruth_l3DSig_vec.back();
  b_Rec_legDR = b_hadTruth_legDR_vec.back();
  b_Rec_angleXY = b_hadTruth_angleXY_vec.back();
  b_Rec_angleXYZ = b_hadTruth_angleXYZ_vec.back();
  b_Rec_chi2 = b_hadTruth_chi2_vec.back();
  b_Rec_dca = b_hadTruth_dca_vec.back();
  b_Rec_dau1_chi2 = b_hadTruth_dau1_chi2_vec.back();
  b_Rec_dau1_ipsigXY = b_hadTruth_dau1_ipsigXY_vec.back();
  b_Rec_dau1_ipsigZ = b_hadTruth_dau1_ipsigZ_vec.back();
  b_Rec_dau1_pt = b_hadTruth_dau1_pt_vec.back();
  b_Rec_dau2_chi2 = b_hadTruth_dau2_chi2_vec.back();
  b_Rec_dau2_ipsigXY = b_hadTruth_dau2_ipsigXY_vec.back();
  b_Rec_dau2_ipsigZ = b_hadTruth_dau2_ipsigZ_vec.back();
  b_Rec_dau2_pt = b_hadTruth_dau2_pt_vec.back();
  b_Rec_bdt_score = m_hadReader->EvaluateMVA("KS_BDT");

  m_hadtrForTMVA->Fill();
}

void vtsAnalyser::SetMVAReader() {

  m_hadReader = new TMVA::Reader();            

  m_hadReader->AddVariable("d",            &b_Rec_d);
  m_hadReader->AddVariable("pt",           &b_Rec_pt);
  m_hadReader->AddVariable("eta",          &b_Rec_eta);
  m_hadReader->AddVariable("phi",          &b_Rec_phi);
  m_hadReader->AddVariable("mass",         &b_Rec_mass);
  m_hadReader->AddVariable("lxy",          &b_Rec_lxy);
  m_hadReader->AddVariable("lxySig",       &b_Rec_lxySig);
  m_hadReader->AddVariable("l3D",          &b_Rec_l3D);
  m_hadReader->AddVariable("l3DSig",       &b_Rec_l3DSig);
  m_hadReader->AddVariable("legDR",        &b_Rec_legDR);
  m_hadReader->AddVariable("angleXY",      &b_Rec_angleXY);
  m_hadReader->AddVariable("angleXYZ",     &b_Rec_angleXYZ);
  m_hadReader->AddVariable("chi2",         &b_Rec_chi2);
  m_hadReader->AddVariable("dca",          &b_Rec_dca);
  m_hadReader->AddVariable("dau1_chi2",    &b_Rec_dau1_chi2);
  m_hadReader->AddVariable("dau1_ipsigXY", &b_Rec_dau1_ipsigXY);
  m_hadReader->AddVariable("dau1_ipsigZ",  &b_Rec_dau1_ipsigZ);
  m_hadReader->AddVariable("dau1_pt",      &b_Rec_dau1_pt);
  m_hadReader->AddVariable("dau2_chi2",    &b_Rec_dau2_chi2);
  m_hadReader->AddVariable("dau2_ipsigXY", &b_Rec_dau2_ipsigXY);
  m_hadReader->AddVariable("dau2_ipsigZ",  &b_Rec_dau2_ipsigZ);
  m_hadReader->AddVariable("dau2_pt",      &b_Rec_dau2_pt);

  m_hadReader->BookMVA("KS_BDT", "/cms/ldap_home/wjjang/wj_nanoAOD_CMSSW_9_4_4/src/nano/analysis/test/vts/tmva/dataset/Had/pp_real_vs_fake/weights/vts_dR_04_Had_BDT.weights.xml");

  m_jetReader = new TMVA::Reader();
  m_jetReader->AddVariable("pt",    &b_pt);
  m_jetReader->AddVariable("eta",   &b_eta);
  m_jetReader->AddVariable("phi",   &b_phi);
  m_jetReader->AddVariable("mass",  &b_mass);
  m_jetReader->AddVariable("c_x1",  &b_c_x1);
  m_jetReader->AddVariable("c_x2",  &b_c_x2);
  m_jetReader->AddVariable("c_x3",  &b_c_x3);
  m_jetReader->AddVariable("n_x1",  &b_n_x1);
  m_jetReader->AddVariable("n_x2",  &b_n_x2);
  m_jetReader->AddVariable("n_x3",  &b_n_x3);
  m_jetReader->AddVariable("cmult", &b_cmult);
  m_jetReader->AddVariable("nmult", &b_nmult);
  m_jetReader->AddVariable("axis1", &b_axis1);
  m_jetReader->AddVariable("axis2", &b_axis2);
  m_jetReader->AddVariable("ptD",   &b_ptD);
  m_jetReader->AddVariable("area",  &b_area);
  m_jetReader->AddVariable("CSVV2", &b_CSVV2);

  m_jetReader->BookMVA("Jet_BDT_highest", "/cms/ldap_home/wjjang/wj_nanoAOD_CMSSW_9_4_4/src/nano/analysis/test/vts/tmva/dataset/Jet/pp_combined_J_BDT_highest/weights/vts_dR_04_Jet_BDT.weights.xml");


  m_jksReader = new TMVA::Reader();
  m_jksReader->AddVariable("pt",    &b_pt);
  m_jksReader->AddVariable("eta",   &b_eta);
  m_jksReader->AddVariable("phi",   &b_phi);
  m_jksReader->AddVariable("mass",  &b_mass);
  m_jksReader->AddVariable("c_x1",  &b_c_x1);
  m_jksReader->AddVariable("c_x2",  &b_c_x2);
  m_jksReader->AddVariable("c_x3",  &b_c_x3);
  m_jksReader->AddVariable("n_x1",  &b_n_x1);
  m_jksReader->AddVariable("n_x2",  &b_n_x2);
  m_jksReader->AddVariable("n_x3",  &b_n_x3);
  m_jksReader->AddVariable("cmult", &b_cmult);
  m_jksReader->AddVariable("nmult", &b_nmult);
  m_jksReader->AddVariable("axis1", &b_axis1);
  m_jksReader->AddVariable("axis2", &b_axis2);
  m_jksReader->AddVariable("ptD",   &b_ptD);
  m_jksReader->AddVariable("area",  &b_area);
  m_jksReader->AddVariable("CSVV2", &b_CSVV2);

  m_jksReader->AddVariable("KS_d",                &b_KS_d);
  m_jksReader->AddVariable("KS_pt",               &b_KS_pt);
  m_jksReader->AddVariable("KS_eta",              &b_KS_eta);
  m_jksReader->AddVariable("KS_phi",              &b_KS_phi);
  m_jksReader->AddVariable("KS_mass",             &b_KS_mass);
  m_jksReader->AddVariable("KS_lxy",              &b_KS_lxy);
  m_jksReader->AddVariable("KS_lxySig",           &b_KS_lxySig);
  m_jksReader->AddVariable("KS_l3D",              &b_KS_l3D);
  m_jksReader->AddVariable("KS_l3DSig",           &b_KS_l3DSig);
  m_jksReader->AddVariable("KS_legDR",            &b_KS_legDR);
  m_jksReader->AddVariable("KS_angleXY",          &b_KS_angleXY);
  m_jksReader->AddVariable("KS_angleXYZ",         &b_KS_angleXYZ);
  m_jksReader->AddVariable("KS_chi2",             &b_KS_chi2);
  m_jksReader->AddVariable("KS_dca",              &b_KS_dca);
  m_jksReader->AddVariable("KS_dau1_chi2",        &b_KS_dau1_chi2);
  m_jksReader->AddVariable("KS_dau1_ipsigXY",     &b_KS_dau1_ipsigXY);
  m_jksReader->AddVariable("KS_dau1_ipsigZ",      &b_KS_dau1_ipsigZ);
  m_jksReader->AddVariable("KS_dau1_pt",          &b_KS_dau1_pt);
  m_jksReader->AddVariable("KS_dau2_chi2",        &b_KS_dau2_chi2);
  m_jksReader->AddVariable("KS_dau2_ipsigXY",     &b_KS_dau2_ipsigXY);
  m_jksReader->AddVariable("KS_dau2_ipsigZ",      &b_KS_dau2_ipsigZ);
  m_jksReader->AddVariable("KS_dau2_pt",          &b_KS_dau2_pt);
  m_jksReader->AddVariable("KS_best_bdt",         &b_KS_best_bdt);
  m_jksReader->AddVariable("KS_x",                &b_KS_x);

  m_jksReader->BookMVA("JKS_BDT_highest", "/cms/ldap_home/wjjang/wj_nanoAOD_CMSSW_9_4_4/src/nano/analysis/test/vts/tmva/dataset/JKS/pp_combined_JKS_BDT_highest/weights/vts_dR_04_Jet_BDT.weights.xml");

}
