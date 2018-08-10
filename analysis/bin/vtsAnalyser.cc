#include <iostream>
#include <string>
#include <vector>
#include <TLorentzVector.h>
#include <sys/stat.h>
#include "nano/analysis/interface/vtsAnalyser.h"
#include "TFile.h"
#include "TTree.h"

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
    Bool_t isMC = (sampleName.find("Run") == std::string::npos);
    string outFileDir = hostDir + getenv("USER") + "/" + jobName + "/" + sampleName;
    for (Int_t i = 3; i < argc; i++) {
      auto inFileName = argv[i];
      if (!isMC) { 
        Bool_t isDL = false;
        if (string(inFileName).find("DoubleElectron") != std::string::npos || string(inFileName).find("DoubleMuon") != std::string::npos) isDL = true;
        Bool_t isSLE = (string(inFileName).find("SingleElectron") != std::string::npos);
        Bool_t isSLM = (string(inFileName).find("SingleMuon") != std::string::npos);
        TFile *inFile = TFile::Open(inFileName, "READ");
        TTree *inTree = (TTree*) inFile->Get("Events");
        vtsAnalyser ana(inTree, inTree, inTree, isMC, isDL, isSLE, isSLM);
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
      /* auto sampleType = getType(dirName); */

      TFile *inFile = TFile::Open(inFileName, "READ");
      TTree *inTree = (TTree*) inFile->Get("Events");
      TString hadFileName = dirName + "/HADAOD/" + fileName;
      TString hadTruthFileName = dirName + "/HADTRUTHAOD/" + fileName;
/*
      TFile *hadFile = TFile::Open(hadFileName, "READ");
      TTree *hadTree = (TTree*) hadFile->Get("Events");
*/
      TFile *hadTruthFile = TFile::Open(hadTruthFileName, "READ");
      TTree *hadTruthTree = (TTree*) hadTruthFile->Get("Events");

      cout << "dirName : " << dirName << " fileName : " << fileName << endl;
      vtsAnalyser ana(inTree,hadTruthTree,hadTruthTree,isMC,false,false,false); // you don't need to use hadTree
      string outFileName = outFileDir+"/nanotree_"+fileName;
//      string outFileName = "/"+jobName+"/"+sampleName+"/nanotree_"+fileName;
            
      ana.setOutput(outFileName);
      ana.Loop();
    }
  }
}

void vtsAnalyser::Loop() {
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntries();

  /* Events loop */
  for (Long64_t iev=0; iev<nentries; iev++) {
    b_jet_start = -1; b_jet_end = -1;
    fChain->GetEntry(iev);
    if (h_fChain) h_fChain->GetEntry(iev);
    if (ht_fChain) ht_fChain->GetEntry(iev);
    cout << "event : " << iev << endl;
    if (iev%10000 == 0) cout << iev << "/" << nentries << endl;
    ResetBranch();
    b_nJet = nJet;
    b_nSelJet = jetSelection().size();
    int passedStep = EventSelection();
    if (passedStep >= 4) { // -1 : No event selection, 0 : PV cut, reco lepton cut and so on, 1~4 : step 1 ~ 4
      b_passedEvent = true; b_nSelJetEv = jetSelection().size();
      MatchingForMC();
      HadronAnalysis();
      GenHadronAnalysis();
      GenAnalysis();
      RecAnalysis();
      JetAnalysis();
      CollectVar();
      ScoreTMVA(outtrForTMVA);
    } else b_passedEvent = false; 
    m_tree->Fill();
    cout << " chk : " << b_jet_start << " " << b_jet_end << endl;
  }
}

void vtsAnalyser::setOutput(std::string outFileName) {
  m_output = TFile::Open(outFileName.c_str(), "recreate");
  m_tree = new TTree("event", "event");
  m_outtrForTMVA = new TTree("MVA_score", "MVA_score");
  MakeBranch(m_tree);

  h_nevents = new TH1D("nevents", "nevents", 1, 0, 1);
  h_genweights = new TH1D("genweight", "genweight", 1, 0, 1);
  h_weights = new TH1D("weight", "weight", 1, 0, 1);
  h_cutFlow = new TH1D("cutflow", "cutflow", 11, -0.5, 10.5);
}

void vtsAnalyser::MakeBranch(TTree* t) {
  outtrForTMVA->Branch("isSJet", &b_isSJet, "isSJet/I");
  outtrForTMVA->Branch("isBJet", &b_isBJet, "isBJet/I");
  outtrForTMVA->Branch("isHighest", &b_isHighest, "isHighest/I");
  outtrForTMVA->Branch("isClosestToLep", &b_isClosestToLep, "isClosestToLep/I");
  outtrForTMVA->Branch("cmult", &b_cmult, "cmult/I");
  outtrForTMVA->Branch("nmult", &b_nmult, "nmult/I");
  outtrForTMVA->Branch("pt", &b_pt, "pt/F");
  outtrForTMVA->Branch("eta", &b_eta, "eta/F");
  outtrForTMVA->Branch("phi", &b_phi, "phi/F");
  outtrForTMVA->Branch("mass", &b_mass, "mass/F");
  outtrForTMVA->Branch("c_x1", &b_c_x1, "c_x1/F");
  outtrForTMVA->Branch("c_x2", &b_c_x2, "c_x2/F");
  outtrForTMVA->Branch("c_x3", &b_c_x3, "c_x3/F");
  outtrForTMVA->Branch("n_x1", &b_n_x1, "n_x1/F");
  outtrForTMVA->Branch("n_x2", &b_n_x2, "n_x2/F");
  outtrForTMVA->Branch("n_x3", &b_n_x3, "n_x3/F");
  outtrForTMVA->Branch("axis1", &b_axis1, "axis1/F");
  outtrForTMVA->Branch("axis2", &b_axis2, "axis2/F");
  outtrForTMVA->Branch("ptD", &b_ptD, "ptD/F");
  outtrForTMVA->Branch("area", &b_area, "area/F");
  outtrForTMVA->Branch("CSVV2", &b_CSVV2, "CSVV2/F");

  
  #define Branch_(type, name, suffix) t->Branch(#name, &(b_##name), #name "/" #suffix);
  #define BranchI(name) Branch_(Int_t, name, I)
  #define BranchF(name) Branch_(Float_t, name, F)
  #define BranchO(name) Branch_(Bool_t, name, O)
  #define BranchV_(type, name) t->Branch(#name, "vector<"#type">", &(b_##name));
  #define BranchVI(name) BranchV_(Int_t, name); 
  #define BranchVF(name) BranchV_(Float_t, name);
  #define BranchVO(name) BranchV_(Bool_t, name);
  #define BranchTLV(name) t->Branch(#name, "TLorentzVector", &(b_##name));

  BranchI(nvertex); BranchI(channel); BranchI(njet) BranchF(met); BranchI(step); BranchO(passedEvent); BranchI(nJet); BranchI(nSelJet); BranchI(nSelJetEv); // njet is not nJet
  BranchI(jet_start); BranchI(jet_end);
  BranchI(hadTruth_nMatched); BranchI(hadTruth_nTrueDau); 
  BranchO(hadTruth_isHadFromTop); BranchI(hadTruth_isHadFromTsb); BranchO(hadTruth_isHadFromW); BranchO(hadTruth_isHadFromS); BranchO(hadTruth_isHadFromC); BranchO(hadTruth_isHadFromB);

  /* For MatchingForMC() */
  BranchF(Jet_dr_closest_s);    BranchF(Jet_dr_closest_b);
  BranchF(SelJet_dr_closest_s); BranchF(SelJet_dr_closest_b);
  BranchF(GenJet_dr_closest_s); BranchF(GenJet_dr_closest_b);

  BranchO(GenSJet);             BranchO(GenBJet);             BranchO(GenBothJet);             BranchO(RecSJet);             BranchO(RecBJet);             BranchO(RecBothJet); 
  BranchO(GenSJetClosestToLep); BranchO(GenBJetClosestToLep); BranchO(GenBothJetClosestToLep); BranchO(RecSJetClosestToLep); BranchO(RecBJetClosestToLep); BranchO(RecBothJetClosestToLep);  
  BranchO(GenSJetIsHighest);    BranchO(GenBJetIsHighest);    BranchO(GenBothJetIsHighest);    BranchO(RecSJetIsHighest);    BranchO(RecBJetIsHighest);    BranchO(RecBothJetIsHighest); 

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

  /* For GenHadronAnalysis() */
  BranchVI(genHadron_isGenFrom_vec); BranchVO(genHadron_isGenFromTop_vec); BranchVO(genHadron_inVol_vec);
  BranchVF(genHadron_d_vec); BranchVF(genHadron_pt_vec); BranchVF(genHadron_eta_vec); BranchVF(genHadron_phi_vec); BranchVF(genHadron_mass_vec);
  BranchVF(genHadron_dau1_pt_vec); BranchVF(genHadron_dau1_eta_vec); BranchVF(genHadron_dau1_phi_vec); BranchVI(genHadron_dau1_pdgId_vec);
  BranchVF(genHadron_dau2_pt_vec); BranchVF(genHadron_dau2_eta_vec); BranchVF(genHadron_dau2_phi_vec); BranchVI(genHadron_dau2_pdgId_vec);
  BranchVF(genHadron_x_closest_j_vec);  BranchVF(genHadron_dr_closest_j_vec);  BranchVF(genHadron_x_highest_j_vec);  BranchVF(genHadron_dr_highest_j_vec);
  BranchVF(genHadron_x_closest_gj_vec); BranchVF(genHadron_dr_closest_gj_vec); BranchVF(genHadron_x_highest_gj_vec); BranchVF(genHadron_dr_highest_gj_vec); 
  BranchVI(genHadron_isClosestPair_xOrder_j_vec);  BranchVI(genHadron_isHighestPair_xOrder_j_vec);  
  BranchVI(genHadron_isClosestPair_xOrder_gj_vec); BranchVI(genHadron_isHighestPair_xOrder_gj_vec);

  BranchVI(nSJet_vec); BranchVI(nBJet_vec); BranchVI(nGenSJet_vec); BranchVI(nGenBJet_vec);

  /* For GenAnalysis() */
  BranchVI(GenPart_isGenFrom_vec); BranchVO(GenPart_isGenFromTop_vec); BranchVO(GenPart_isGenFromW_vec); BranchVO(GenPart_isFromKstar_vec);
  BranchVF(GenPart_d_vec); BranchVF(GenPart_pt_vec); BranchVF(GenPart_eta_vec); BranchVF(GenPart_phi_vec); BranchVF(GenPart_mass_vec);
  BranchVF(GenPart_x_closest_j_vec);  BranchVF(GenPart_dr_closest_j_vec);  BranchVF(GenPart_x_highest_j_vec);  BranchVF(GenPart_dr_highest_j_vec);
  BranchVF(GenPart_x_closest_gj_vec); BranchVF(GenPart_dr_closest_gj_vec); BranchVF(GenPart_x_highest_gj_vec); BranchVF(GenPart_dr_highest_gj_vec);
  BranchVI(GenPart_isClosestPair_xOrder_j_vec);  BranchVI(GenPart_isHighestPair_xOrder_j_vec); 
  BranchVI(GenPart_isClosestPair_xOrder_gj_vec); BranchVI(GenPart_isHighestPair_xOrder_gj_vec);

  /* For RecAnalysis() */
  BranchVI(hadTruth_pdgId_vec); BranchVI(hadTruth_nMatched_vec); BranchVI(hadTruth_isFrom_vec);
  BranchVO(hadTruth_isHadFromTop_vec); BranchVO(hadTruth_isHadFromW_vec); BranchVO(hadTruth_isHadFromS_vec); BranchVO(hadTruth_isHadFromC_vec); BranchVO(hadTruth_isHadFromB_vec);
  BranchVF(hadTruth_d_vec); BranchVF(hadTruth_pt_vec); BranchVF(hadTruth_eta_vec); BranchVF(hadTruth_phi_vec); BranchVF(hadTruth_mass_vec);
  BranchVF(hadTruth_lxy_vec); BranchVF(hadTruth_lxySig_vec); BranchVF(hadTruth_l3D_vec); BranchVF(hadTruth_l3DSig_vec); BranchVF(hadTruth_legDR_vec);
  BranchVF(hadTruth_angleXY_vec); BranchVF(hadTruth_angleXYZ_vec); BranchVF(hadTruth_chi2_vec); BranchVF(hadTruth_dca_vec);
  BranchVF(hadTruth_dau1_chi2_vec); BranchVF(hadTruth_dau1_ipsigXY_vec); BranchVF(hadTruth_dau1_ipsigZ_vec); BranchVF(hadTruth_dau1_pt_vec);
  BranchVF(hadTruth_dau2_chi2_vec); BranchVF(hadTruth_dau2_ipsigXY_vec); BranchVF(hadTruth_dau2_ipsigZ_vec); BranchVF(hadTruth_dau2_pt_vec);
  BranchVF(hadTruth_x_closest_j_vec);  BranchVF(hadTruth_dr_closest_j_vec);  BranchVF(hadTruth_x_highest_j_vec);  BranchVF(hadTruth_dr_highest_j_vec);
  BranchVF(hadTruth_x_closest_gj_vec); BranchVF(hadTruth_dr_closest_gj_vec); BranchVF(hadTruth_x_highest_gj_vec); BranchVF(hadTruth_dr_highest_gj_vec);
  BranchVI(hadTruth_isClosestPair_xOrder_j_vec);  BranchVI(hadTruth_isHighestPair_xOrder_j_vec);
  BranchVI(hadTruth_isClosestPair_xOrder_gj_vec); BranchVI(hadTruth_isHighestPair_xOrder_gj_vec);

  /* For JetAnalysis() */
  BranchVI(Jet_isCorrectMat);

  /* For CollectVar() */
  BranchF(MET_pt); BranchF(MET_phi); BranchF(MET_sumEt);
  BranchI(lep1_pid); BranchI(lep2_pid);
  BranchTLV(lep1); BranchTLV(lep2);
  BranchTLV(dilep);
}

void vtsAnalyser::ResetBranch() {
  Reset();
  
  b_nJet = -1; b_nSelJet = -1; b_nSelJetEv = -1; 
  b_passedEvent = false;

  b_hadTruth_nMatched = -1; b_hadTruth_nTrueDau = -1;
  b_hadTruth_isHadFromTsb = -1;
  b_hadTruth_isHadFromTop = false; b_hadTruth_isHadFromW = false; b_hadTruth_isHadFromS = false; b_hadTruth_isHadFromC = false; b_hadTruth_isHadFromB = false;

  /* For MatchingForMC() */
  m_tqMC.clear(); m_wqMC.clear(); 
  m_qjMapForMC.clear(); m_qgjMapForMC.clear();
  m_recJet.clear(); m_genJet.clear();
  m_closestRecJetForLep1.clear(); m_closestRecJetForLep2.clear(); m_closestGenJetForLep1.clear(); m_closestGenJetForLep2.clear();

  b_Jet_dr_closest_s = -1;    b_Jet_dr_closest_b = -1;
  b_SelJet_dr_closest_s = -1; b_SelJet_dr_closest_b = -1;
  b_GenJet_dr_closest_s = -1; b_GenJet_dr_closest_b = -1;
  b_GenSJet = false;             b_GenBJet = false;             b_GenBothJet = false;             b_RecSJet = false;             b_RecBJet = false;             b_RecBothJet = false;
  b_GenSJetClosestToLep = false; b_GenBJetClosestToLep = false; b_GenBothJetClosestToLep = false; b_RecSJetClosestToLep = false; b_RecBJetClosestToLep = false; b_RecBothJetClosestToLep = false;
  b_GenSJetIsHighest = false;    b_GenBJetIsHighest = false;    b_GenBothJetIsHighest = false;    b_RecSJetIsHighest = false;    b_RecBJetIsHighest = false;    b_RecBothJetIsHighest = false;

  /* For GenHadronAnalysis() */
  b_genHadron_isGenFrom_vec.clear(); b_genHadron_isGenFromTop_vec.clear(); b_genHadron_inVol_vec.clear();
  b_genHadron_d_vec.clear(); b_genHadron_pt_vec.clear(); b_genHadron_eta_vec.clear(); b_genHadron_phi_vec.clear(); b_genHadron_mass_vec.clear();
  b_genHadron_dau1_pt_vec.clear(); b_genHadron_dau1_eta_vec.clear(); b_genHadron_dau1_phi_vec.clear(); b_genHadron_dau1_pdgId_vec.clear();
  b_genHadron_dau2_pt_vec.clear(); b_genHadron_dau2_eta_vec.clear(); b_genHadron_dau2_phi_vec.clear(); b_genHadron_dau2_pdgId_vec.clear();
  b_genHadron_x_closest_j_vec.clear();  b_genHadron_dr_closest_j_vec.clear();  b_genHadron_x_highest_j_vec.clear();  b_genHadron_dr_highest_j_vec.clear();
  b_genHadron_x_closest_gj_vec.clear(); b_genHadron_dr_closest_gj_vec.clear(); b_genHadron_x_highest_gj_vec.clear(); b_genHadron_dr_highest_gj_vec.clear();
  b_genHadron_isClosestPair_xOrder_j_vec.clear();  b_genHadron_isHighestPair_xOrder_j_vec.clear();
  b_genHadron_isClosestPair_xOrder_gj_vec.clear(); b_genHadron_isHighestPair_xOrder_gj_vec.clear();

  b_nSJet_vec.clear(); b_nBJet_vec.clear(); b_nGenSJet_vec.clear(); b_nGenBJet_vec.clear();

  /* For GenAnalysis() */
  b_GenPart_isGenFrom_vec.clear(); b_GenPart_isGenFromTop_vec.clear(); b_GenPart_isGenFromW_vec.clear(); b_GenPart_isFromKstar_vec.clear();
  b_GenPart_d_vec.clear(); b_GenPart_pt_vec.clear(); b_GenPart_eta_vec.clear(); b_GenPart_phi_vec.clear(); b_GenPart_mass_vec.clear();
  b_GenPart_x_closest_j_vec.clear();  b_GenPart_dr_closest_j_vec.clear();  b_GenPart_x_highest_j_vec.clear();  b_GenPart_dr_highest_j_vec.clear();
  b_GenPart_x_closest_gj_vec.clear(); b_GenPart_dr_closest_gj_vec.clear(); b_GenPart_x_highest_gj_vec.clear(); b_GenPart_dr_highest_gj_vec.clear();
  b_GenPart_isClosestPair_xOrder_j_vec.clear();  b_GenPart_isHighestPair_xOrder_j_vec.clear();
  b_GenPart_isClosestPair_xOrder_gj_vec.clear(); b_GenPart_isHighestPair_xOrder_gj_vec.clear();

  /* For RecAnalysis() */
  b_hadTruth_pdgId_vec.clear(); b_hadTruth_nMatched_vec.clear(); b_hadTruth_isFrom_vec.clear(); 
  b_hadTruth_isHadFromTop_vec.clear(); b_hadTruth_isHadFromW_vec.clear(); b_hadTruth_isHadFromS_vec.clear(); b_hadTruth_isHadFromC_vec.clear(); b_hadTruth_isHadFromB_vec.clear(); 
  b_hadTruth_d_vec.clear(); b_hadTruth_pt_vec.clear(); b_hadTruth_eta_vec.clear(); b_hadTruth_phi_vec.clear(); b_hadTruth_mass_vec.clear();
  b_hadTruth_lxy_vec.clear(); b_hadTruth_lxySig_vec.clear(); b_hadTruth_l3D_vec.clear(); b_hadTruth_l3DSig_vec.clear(); b_hadTruth_legDR_vec.clear();
  b_hadTruth_angleXY_vec.clear(); b_hadTruth_angleXYZ_vec.clear(); b_hadTruth_chi2_vec.clear(); b_hadTruth_dca_vec.clear();
  b_hadTruth_dau1_chi2_vec.clear(); b_hadTruth_dau1_ipsigXY_vec.clear(); b_hadTruth_dau1_ipsigZ_vec.clear(); b_hadTruth_dau1_pt_vec.clear();
  b_hadTruth_dau2_chi2_vec.clear(); b_hadTruth_dau2_ipsigXY_vec.clear(); b_hadTruth_dau2_ipsigZ_vec.clear(); b_hadTruth_dau2_pt_vec.clear();
  b_hadTruth_x_closest_j_vec.clear();  b_hadTruth_dr_closest_j_vec.clear();  b_hadTruth_x_highest_j_vec.clear();  b_hadTruth_dr_highest_j_vec.clear();
  b_hadTruth_x_closest_gj_vec.clear(); b_hadTruth_dr_closest_gj_vec.clear(); b_hadTruth_x_highest_gj_vec.clear(); b_hadTruth_dr_highest_gj_vec.clear();
  b_hadTruth_isClosestPair_xOrder_j_vec.clear();  b_hadTruth_isHighestPair_xOrder_j_vec.clear();
  b_hadTruth_isClosestPair_xOrder_gj_vec.clear(); b_hadTruth_isHighestPair_xOrder_gj_vec.clear();

  /* For JetAnalysis() */
  b_Jet_isCorrectMat.clear();

  /* For CollectVar() */
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
/*
 // Finding s quark using GenPart_status
  for (unsigned int i=0; i<nGenPart; ++i) {
    if (std::abs(GenPart_status[i]) != 23 ) continue;
    if (abs(GenPart_pdgId[i]) == 3 || abs(GenPart_pdgId[i]) == 5 || abs(GenPart_pdgId[i]) == 4) {
      if ( (abs(GenPart_pdgId[GenPart_genPartIdxMother[i]]) == 6 && GenPart_status[GenPart_genPartIdxMother[i]] == 62 ) ) { 
        m_tqMC.push_back(i); cout << " [ " << i << " / " << nGenPart << " ]  >>>>>>>>>>>> got it <<<<<<<<<<< " << endl;
      }
      if ( (abs(GenPart_pdgId[GenPart_genPartIdxMother[i]]) == 24 && (GenPart_status[GenPart_genPartIdxMother[i]] == 22 || GenPart_status[GenPart_genPartIdxMother[i]] == 52)) ) {
        m_wqMC.push_back(i);
      }
    }
  }
*/

  // Finding s quark using GenPart_statusFlag (the part of W boson isn't completed)  >>>>>>>> pythia
  for (unsigned int i=0; i<nGenPart; ++i) {
    if (std::abs(GenPart_statusFlags[i]) != 22913 && std::abs(GenPart_statusFlags[i]) != 4481) continue;
    if (abs(GenPart_pdgId[i]) == 3 || abs(GenPart_pdgId[i]) == 5 || abs(GenPart_pdgId[i]) == 4) {
      if ( (abs(GenPart_pdgId[GenPart_genPartIdxMother[i]]) == 6 && GenPart_statusFlags[GenPart_genPartIdxMother[i]] == 10497 ) ) { 
        m_tqMC.push_back(i);// cout << " [ " << i << " / " << nGenPart << " ]  >>>>>>>>>>>> got it <<<<<<<<<<< " << endl;
      }
      if ( (abs(GenPart_pdgId[GenPart_genPartIdxMother[i]]) == 24 && (GenPart_status[GenPart_genPartIdxMother[i]] == 22 || GenPart_status[GenPart_genPartIdxMother[i]] == 52)) ) {
        m_wqMC.push_back(i);
      }
    }
  }

/*
  // Finding s quark using GenPart_statusFlag (the part of W boson isn't completed)  >>>>>>>> herwig
  for (unsigned int i=0; i<nGenPart; ++i) {
    if (std::abs(GenPart_statusFlags[i]) != 4097) continue;
    if (abs(GenPart_pdgId[i]) == 3 || abs(GenPart_pdgId[i]) == 5 || abs(GenPart_pdgId[i]) == 4) {
//      cout << " [ " << i << " / " << nGenPart << " ] chk : " << GenPart_pdgId[i] << " " << GenPart_status[i] << " " << GenPart_statusFlags[i] << " || mom chk : idx : " << GenPart_genPartIdxMother[i] << " " <<  GenPart_pdgId[GenPart_genPartIdxMother[i]] << " " << GenPart_status[GenPart_genPartIdxMother[i]] << " " << GenPart_statusFlags[GenPart_genPartIdxMother[i]] << endl;
      if ( (abs(GenPart_pdgId[GenPart_genPartIdxMother[i]]) == 6 && GenPart_statusFlags[GenPart_genPartIdxMother[i]] == 8193 ) ) { 
        m_tqMC.push_back(i);// cout << " [ " << i << " / " << nGenPart << " ]  >>>>>>>>>>>> got it <<<<<<<<<<< " << endl;
      }
      if ( (abs(GenPart_pdgId[GenPart_genPartIdxMother[i]]) == 24 && (GenPart_status[GenPart_genPartIdxMother[i]] == 22 || GenPart_status[GenPart_genPartIdxMother[i]] == 52)) ) {
        m_wqMC.push_back(i);
      }
    }
  }
*/
  
  if (m_tqMC.size() < 2) { cout << " >>>>> it's not a tsWbW event. save nothing." << endl; return; }
  /* Select quarks originated from t->qW */
  int tq1 = -1; int tq2 = -1;
  if (abs(GenPart_pdgId[m_tqMC[0]]) == 3) { tq1 = m_tqMC[0]; tq2 = m_tqMC[1]; }
  else if (abs(GenPart_pdgId[m_tqMC[0]]) == 5) { tq1 = m_tqMC[1]; tq2 = m_tqMC[0]; }
  TLorentzVector tq1_tlv; 
  tq1_tlv.SetPtEtaPhiM(GenPart_pt[tq1], GenPart_eta[tq1], GenPart_phi[tq1], GenPart_mass[tq1]);
  TLorentzVector tq2_tlv; 
  tq2_tlv.SetPtEtaPhiM(GenPart_pt[tq2], GenPart_eta[tq2], GenPart_phi[tq2], GenPart_mass[tq2]);
  /* Select quarks from W->qq (not yet used) */
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
  if (m_genJet[0].drsj <= m_jetConeSize) m_qgjMapForMC[m_genJet[0].idx] = GenPart_pdgId[tq1]; // if jet is inside dRCut == 0.4, then label the tag about s/b
  sort(m_genJet.begin(), m_genJet.end(), [] (jetInfo a, jetInfo b) { return (a.drbj < b.drbj); } );
  if (m_genJet[0].drbj <= m_jetConeSize) m_qgjMapForMC[m_genJet[0].idx] = GenPart_pdgId[tq2];

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
    sort(m_recJet.begin(), m_recJet.end(), [] (jetInfo a, jetInfo b) { return (a.drl1j < b.drl1j); } );
    m_closestRecJetForLep1[m_recJet[0].idx] = b_lep1_pid;
    sort(m_recJet.begin(), m_recJet.end(), [] (jetInfo a, jetInfo b) { return (a.drl2j < b.drl2j); } );
    m_closestRecJetForLep2[m_recJet[0].idx] = b_lep2_pid;
    /* Find closest sel jet to s or b */ 
    sort(m_recJet.begin(), m_recJet.end(), [] (jetInfo a, jetInfo b) { return (a.drsj < b.drsj); } );
    if (m_recJet[0].drsj <= m_jetConeSize) m_qjMapForMC[m_recJet[0].idx] = GenPart_pdgId[tq1]; // if jet is inside dRCut == 0.4, then label the tag about s/b
    sort(m_recJet.begin(), m_recJet.end(), [] (jetInfo a, jetInfo b) { return (a.drbj < b.drbj); } );
    if (m_recJet[0].drbj <= m_jetConeSize) m_qjMapForMC[m_recJet[0].idx] = GenPart_pdgId[tq2];
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

void vtsAnalyser::HadronAnalysis() {
  std::vector<std::vector<struct HadStat>> JetCollection;
  std::vector<struct HadStat> Had;
  /* Reco Jet & Reco KS matching */
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

void vtsAnalyser::GenHadronAnalysis() {
  std::vector<jetks> recPair1; std::vector<jetks> genPair1;
  std::vector<jetks> recPair2; std::vector<jetks> genPair2;
  int nGenHadKS = 0;
  b_genHadron_dr_closest_j_vec.resize(ngenHadron,-1);
  b_genHadron_x_closest_j_vec.resize(ngenHadron,-1);
  b_genHadron_dr_closest_gj_vec.resize(ngenHadron,-1);
  b_genHadron_x_closest_gj_vec.resize(ngenHadron,-1);
  b_genHadron_dr_highest_j_vec.resize(ngenHadron,-1);
  b_genHadron_x_highest_j_vec.resize(ngenHadron,-1);
  b_genHadron_dr_highest_gj_vec.resize(ngenHadron,-1);
  b_genHadron_x_highest_gj_vec.resize(ngenHadron,-1);

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
    recPair1.resize(nJet, {0,99,-1}); recPair2.resize(nJet, {0,99,-1});
    for (unsigned int j=0; j<nJet; ++j) { // Loop for all of recoJet
      TLorentzVector jet_tlv;
      jet_tlv.SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
      auto dr = jet_tlv.DeltaR(gen_tlv);      auto x = genHadron_pt[i]/Jet_pt[j];
      auto drMin1 = recPair1[j].dr; auto xMax1 = recPair1[j].x;
      auto drMin2 = recPair2[j].dr; auto xMax2 = recPair2[j].x;
      if (nGenHadKS==1) {
        if (abs(m_qjMapForMC[j]) == 3) b_nSJet_vec.push_back(j);
        else if (abs(m_qjMapForMC[j]) == 5) b_nBJet_vec.push_back(j);
      }
      if (genHadron_isGenHadFromTsb[i] != m_qjMapForMC[j]) continue;
      //if (x > 1.) continue;
      if (dr < drMin1) recPair1[j] = {i,dr,x};
      else if (dr == drMin1 && x > xMax1) recPair1[j] = {i,dr,x};
      if (x > xMax2) recPair2[j] = {i,dr,x};
      else if (x == xMax2 && dr < drMin2) recPair2[j] = {i,dr,x};
    }
    genPair1.resize(nGenJet, {0,99,-1}); genPair2.resize(nGenJet, {0,99,-1});
    for (unsigned int j=0; j<nGenJet; ++j) { // Loop for all of recoJet
      TLorentzVector jet_tlv;
      jet_tlv.SetPtEtaPhiM(GenJet_pt[j], GenJet_eta[j], GenJet_phi[j], GenJet_mass[j]);
      auto dr = jet_tlv.DeltaR(gen_tlv); auto x = genHadron_pt[i]/GenJet_pt[j];
      auto drMin1 = genPair1[j].dr; auto xMax1 = genPair1[j].x;
      auto drMin2 = genPair2[j].dr; auto xMax2 = genPair2[j].x;
      if (nGenHadKS==1) {
        if (abs(m_qgjMapForMC[j]) == 3) b_nGenSJet_vec.push_back(j);
        else if (abs(m_qgjMapForMC[j]) == 5) b_nGenBJet_vec.push_back(j);
      }
      if (genHadron_isGenHadFromTsb[i] != m_qgjMapForMC[j]) continue;
      //if (x > 1.) continue;
      if (dr < drMin1) genPair1[j] = {i,dr,x};
      else if (dr == drMin1 && x > xMax1) genPair1[j] = {i,dr,x};
      if (x > xMax2) genPair2[j] = {i,dr,x};
      else if (x == xMax2 && dr < drMin2) genPair2[j] = {i,dr,x};
    }
    b_genHadron_isClosestPair_xOrder_j_vec.push_back(i);  b_genHadron_isHighestPair_xOrder_j_vec.push_back(i);
    b_genHadron_isClosestPair_xOrder_gj_vec.push_back(i); b_genHadron_isHighestPair_xOrder_gj_vec.push_back(i);
  }
  if (nGenHadKS !=0) {
    /*  Save x and dr for the most closest(highest) KS for each jet */
    for (unsigned int j=0; j<nJet; ++j) {
      int cRecIdx = recPair1[j].idx; int hRecIdx = recPair2[j].idx;
      b_genHadron_dr_closest_j_vec[cRecIdx] = recPair1[j].dr;
      b_genHadron_x_closest_j_vec[cRecIdx] = recPair1[j].x;
      b_genHadron_dr_highest_j_vec[hRecIdx] = recPair2[j].dr;
      b_genHadron_x_highest_j_vec[hRecIdx] = recPair2[j].x;
    }
    for (unsigned int j=0; j<nGenJet; ++j) {
      int cGenIdx = genPair1[j].idx;int hGenIdx = genPair2[j].idx;
      b_genHadron_dr_closest_gj_vec[cGenIdx] = genPair1[j].dr;
      b_genHadron_x_closest_gj_vec[cGenIdx] = genPair1[j].x;
      b_genHadron_dr_highest_gj_vec[hGenIdx] = genPair2[j].dr;
      b_genHadron_x_highest_gj_vec[hGenIdx] = genPair2[j].x;
    }

    /* Give flag ( == index) for the most closest(highest) KS-jet pair per event by x ordering */
    std::sort(recPair1.begin(), recPair1.end(), [] (jetks  a, jetks b) { return (a.x > b.x); } );
    std::sort(genPair1.begin(), genPair1.end(), [] (jetks  a, jetks b) { return (a.x > b.x); } );
    int hRecIdx1 = recPair1[0].idx; int hGenIdx1 = genPair1[0].idx;
    std::replace_if(b_genHadron_isClosestPair_xOrder_j_vec.begin(),  b_genHadron_isClosestPair_xOrder_j_vec.end(),  [&] (int a) { return (a != hRecIdx1); }, -1);
    std::replace_if(b_genHadron_isClosestPair_xOrder_gj_vec.begin(), b_genHadron_isClosestPair_xOrder_gj_vec.end(), [&] (int a) { return (a != hGenIdx1); }, -1);

    std::sort(recPair2.begin(), recPair2.end(), [] (jetks  a, jetks b) { return (a.x > b.x); } );
    std::sort(genPair2.begin(), genPair2.end(), [] (jetks  a, jetks b) { return (a.x > b.x); } );
    int hRecIdx2 = recPair2[0].idx; int hGenIdx2 = genPair2[0].idx;
    std::replace_if(b_genHadron_isHighestPair_xOrder_j_vec.begin(),  b_genHadron_isHighestPair_xOrder_j_vec.end(),  [&] (int a) { return (a != hRecIdx2); }, -1);
    std::replace_if(b_genHadron_isHighestPair_xOrder_gj_vec.begin(), b_genHadron_isHighestPair_xOrder_gj_vec.end(), [&] (int a) { return (a != hGenIdx2); }, -1);
  }
}

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

    recPair1.resize(nJet, {0,99,-1}); recPair2.resize(nJet, {0,99,-1});
    for (unsigned int j=0; j<nJet; ++j) { // Loop for all of recoJet
      TLorentzVector jet_tlv;
      jet_tlv.SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
      auto dr = jet_tlv.DeltaR(gen_tlv);      auto x = GenPart_pt[i]/Jet_pt[j];
      auto drMin1 = recPair1[j].dr; auto xMax1 = recPair1[j].x;
      auto drMin2 = recPair2[j].dr; auto xMax2 = recPair2[j].x;
      if (isFrom != m_qjMapForMC[j]) continue;
      if (x > 1.) continue;
      if (dr < drMin1) recPair1[j] = {i,dr,x};
      else if (dr == drMin1 && x > xMax1) recPair1[j] = {i,dr,x};
      if (x > xMax2) recPair2[j] = {i,dr,x};
      else if (x == xMax2 && dr < drMin2) recPair2[j] = {i,dr,x};
    }
    genPair1.resize(nGenJet, {0,99,-1}); genPair2.resize(nGenJet, {0,99,-1});
    for (unsigned int j=0; j<nGenJet; ++j) { // Loop for all of recoJet
      TLorentzVector jet_tlv;
      jet_tlv.SetPtEtaPhiM(GenJet_pt[j], GenJet_eta[j], GenJet_phi[j], GenJet_mass[j]);
      auto dr = jet_tlv.DeltaR(gen_tlv); auto x = GenPart_pt[i]/GenJet_pt[j];
      auto drMin1 = genPair1[j].dr; auto xMax1 = genPair1[j].x;
      auto drMin2 = genPair2[j].dr; auto xMax2 = genPair2[j].x;
      if (isFrom != m_qgjMapForMC[j]) continue;
      if (x > 1.) continue;
      if (dr < drMin1) genPair1[j] = {i,dr,x};
      else if (dr == drMin1 && x > xMax1) genPair1[j] = {i,dr,x};
      if (x > xMax2) genPair2[j] = {i,dr,x};
      else if (x == xMax2 && dr < drMin2) genPair2[j] = {i,dr,x};
    }
    b_GenPart_isClosestPair_xOrder_j_vec.push_back(i);  b_GenPart_isHighestPair_xOrder_j_vec.push_back(i);
    b_GenPart_isClosestPair_xOrder_gj_vec.push_back(i); b_GenPart_isHighestPair_xOrder_gj_vec.push_back(i);
  }
  if (nGenKS !=0) {
    /*  Save x and dr for the most closest(highest) KS for each jet */
    for (unsigned int j=0; j<nJet; ++j) {
      int cRecIdx = recPair1[j].idx; int hRecIdx = recPair2[j].idx;
      b_GenPart_dr_closest_j_vec[cRecIdx] = recPair1[j].dr;
      b_GenPart_x_closest_j_vec[cRecIdx] = recPair1[j].x;
      b_GenPart_dr_highest_j_vec[hRecIdx] = recPair2[j].dr;
      b_GenPart_x_highest_j_vec[hRecIdx] = recPair2[j].x;
    }
    for (unsigned int j=0; j<nGenJet; ++j) {
      int cGenIdx = genPair1[j].idx;int hGenIdx = genPair2[j].idx;
      b_GenPart_dr_closest_gj_vec[cGenIdx] = genPair1[j].dr;
      b_GenPart_x_closest_gj_vec[cGenIdx] = genPair1[j].x;
      b_GenPart_dr_highest_gj_vec[hGenIdx] = genPair2[j].dr;
      b_GenPart_x_highest_gj_vec[hGenIdx] = genPair2[j].x;
    }
    /* Give flag ( == index) for the most closest(highest) KS-jet pair per event by x ordering */
    std::sort(recPair1.begin(), recPair1.end(), [] (jetks  a, jetks b) { return (a.x > b.x); } );
    std::sort(genPair1.begin(), genPair1.end(), [] (jetks  a, jetks b) { return (a.x > b.x); } );
    int hRecIdx1 = recPair1[0].idx; int hGenIdx1 = genPair1[0].idx;
    std::replace_if(b_GenPart_isClosestPair_xOrder_j_vec.begin(),  b_GenPart_isClosestPair_xOrder_j_vec.end(),  [&] (int a) { return (a != hRecIdx1); }, -1);
    std::replace_if(b_GenPart_isClosestPair_xOrder_gj_vec.begin(), b_GenPart_isClosestPair_xOrder_gj_vec.end(), [&] (int a) { return (a != hGenIdx1); }, -1);

    std::sort(recPair2.begin(), recPair2.end(), [] (jetks  a, jetks b) { return (a.x > b.x); } );
    std::sort(genPair2.begin(), genPair2.end(), [] (jetks  a, jetks b) { return (a.x > b.x); } );
    int hRecIdx2 = recPair2[0].idx; int hGenIdx2 = genPair2[0].idx;
    std::replace_if(b_GenPart_isHighestPair_xOrder_j_vec.begin(),  b_GenPart_isHighestPair_xOrder_j_vec.end(),  [&] (int a) { return (a != hRecIdx2); }, -1);
    std::replace_if(b_GenPart_isHighestPair_xOrder_gj_vec.begin(), b_GenPart_isHighestPair_xOrder_gj_vec.end(), [&] (int a) { return (a != hGenIdx2); }, -1);
  }
}

void vtsAnalyser::RecAnalysis() {
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

    TLorentzVector had_tlv;
    had_tlv.SetPtEtaPhiM(had_pt[i], had_eta[i], had_phi[i], had_mass[i]);
    recPair1.resize(nJet, {0,99,-1}); recPair2.resize(nJet, {0,99,-1});
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
    genPair1.resize(nGenJet, {0,99,-1}); genPair2.resize(nGenJet, {0,99,-1});
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
    b_hadTruth_isClosestPair_xOrder_j_vec.push_back(i);  b_hadTruth_isHighestPair_xOrder_j_vec.push_back(i);
    b_hadTruth_isClosestPair_xOrder_gj_vec.push_back(i); b_hadTruth_isHighestPair_xOrder_gj_vec.push_back(i);
  }
  
  if (nRecKS !=0) {
    /*  Save x and dr for the most closest(highest) KS for each jet */
    for (unsigned int j=0; j<nJet; ++j) {
      int cRecIdx = recPair1[j].idx; int hRecIdx = recPair2[j].idx;
      b_hadTruth_dr_closest_j_vec[cRecIdx] = recPair1[j].dr;
      b_hadTruth_x_closest_j_vec[cRecIdx] = recPair1[j].x;
      b_hadTruth_dr_highest_j_vec[hRecIdx] = recPair2[j].dr;
      b_hadTruth_x_highest_j_vec[hRecIdx] = recPair2[j].x;
    }
    for (unsigned int j=0; j<nGenJet; ++j) {
      int cGenIdx = genPair1[j].idx;int hGenIdx = genPair2[j].idx;
      b_hadTruth_dr_closest_gj_vec[cGenIdx] = genPair1[j].dr;
      b_hadTruth_x_closest_gj_vec[cGenIdx] = genPair1[j].x;
      b_hadTruth_dr_highest_gj_vec[hGenIdx] = genPair2[j].dr;
      b_hadTruth_x_highest_gj_vec[hGenIdx] = genPair2[j].x;
    }

    /* Give flag ( == index) for the most closest(highest) KS-jet pair per event by x ordering */
    
    auto xorder = [] (jetks  a, jetks b) { return (a.x > b.x); };
    int hRecIdx1 = std::max_element(recPair1.begin(), recPair1.end(), xorder)->idx;
    int hGenIdx1 = std::max_element(genPair1.begin(), genPair1.end(), xorder)->idx;
    std::replace_if(b_hadTruth_isClosestPair_xOrder_j_vec.begin(),  b_hadTruth_isClosestPair_xOrder_j_vec.end(),  [&] (int a) { return (a != hRecIdx1); }, -1);
    std::replace_if(b_hadTruth_isClosestPair_xOrder_gj_vec.begin(), b_hadTruth_isClosestPair_xOrder_gj_vec.end(), [&] (int a) { return (a != hGenIdx1); }, -1);

    int hRecIdx2 = std::max_element(recPair2.begin(), recPair2.end(), xorder)->idx;
    int hGenIdx2 = std::max_element(genPair2.begin(), genPair2.end(), xorder)->idx;
    std::replace_if(b_hadTruth_isHighestPair_xOrder_j_vec.begin(),  b_hadTruth_isHighestPair_xOrder_j_vec.end(),  [&] (int a) { return (a != hRecIdx2); }, -1);
    std::replace_if(b_hadTruth_isHighestPair_xOrder_gj_vec.begin(), b_hadTruth_isHighestPair_xOrder_gj_vec.end(), [&] (int a) { return (a != hGenIdx2); }, -1);
  }
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

void vtsAnalyser::ScoreTMVA() {
  auto selectedJet = jetSelection();
  b_jet_start =  m_outtr->GetEntries();
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
      b_isSJet = 0; b_isBJet = 0; b_isHighest = 0; b_isClosestToLep = 0;

      auto j = selectedJet[ij].GetFirstMother();

      if ((j == (int) closest_s_idx) && (closest_s_dr <= m_jetConeSize)) b_isSJet = 1;
      if ((j == (int) closest_b_idx) && (closest_b_dr <= m_jetConeSize)) b_isBJet = 1;
      if ((j == (int) highest_first_idx) || (j == (int) highest_second_idx)) b_isHighest = 1;
      if (j == (int) closest_lep1_idx) b_isClosestToLep = 1;
      if (j == (int) closest_lep2_idx) b_isClosestToLep = 1;

      b_cmult = jetID_cmult[j]; b_nmult = jetID_nmult[j];
      b_pt = Jet_pt[j]; b_eta = Jet_eta[j]; b_phi = Jet_phi[j]; b_mass = Jet_mass[j];
      b_c_x1 = jetID_cpt1[j]/Jet_pt[j]; b_c_x2 = jetID_cpt2[j]/Jet_pt[j]; b_c_x3 = jetID_cpt3[j]/Jet_pt[j];
      b_n_x1 = jetID_npt1[j]/Jet_pt[j]; b_n_x2 = jetID_npt2[j]/Jet_pt[j]; b_n_x3 = jetID_npt3[j]/Jet_pt[j];
      b_axis1 = jetID_axis1[j]; b_axis2 = jetID_axis2[j]; b_ptD = jetID_ptD[j];
      b_area = Jet_area[j]; b_CSVV2 = Jet_btagCSVV2[j];

      m_outtr->Fill();
    }
  } else cout << ">>>> Size of selectedJets is zero <<<< " << endl;
  b_jet_end = m_outtr->GetEntries();
}

