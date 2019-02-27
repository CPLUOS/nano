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

string getDir (const string& path) { return path.substr(0, path.find_last_of("/\\")); }
string getType (const string& path) { return path.substr(path.find("tt")); }
string getFileName(const string& s) {
 char sep = '/';
 size_t i = s.rfind(sep, s.length());
 if (i != string::npos) {
    return s.substr(i+1, s.length() - i);
 }
 return "";
}

int main(int argc, char* argv[])
{
  string hostDir = "root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/";

  if (argc <= 1) {
    cout << "no input file is specified. running with default file." << endl;
    auto inFile = TFile::Open("/xrootd/store/group/nanoAOD/run2_2016v5/tsw/nanoAOD_111.root", "READ");
    auto inTree = (TTree*) inFile->Get("Events");
    vtsAnalyser ana(inTree,inTree,0,true,false,false,false);
    ana.setOutput("nanotree.root");
    ana.Loop();
  } else {
    string jobName    = string(argv[1]);
    string sampleName = string(argv[2]);
    string outFileDir = hostDir + getenv("USER") + "/" + jobName + "/" + sampleName;
    for (Int_t i = 3; i < argc; i++) {
      auto inFileName = argv[i];
      auto fileName = getFileName(argv[i]);
      auto dirName = getDir(getDir(argv[i]));

      /* auto sampleType = getType(dirName); */
      Bool_t isMC = false;
      Bool_t isGenericMC = false;
      if ( (string(inFileName).find("Run2016") == std::string::npos) &&
           (string(inFileName).find("Run2017") == std::string::npos) &&
           (string(inFileName).find("Run2018") == std::string::npos) ) isMC = true;
      if (string(inFileName).find("NANOAOD") != std::string::npos) isGenericMC = false;
      else isGenericMC = true;

      cout << "inFileName : " << inFileName << endl;
      TFile *inFile = TFile::Open(inFileName, "READ");
      TTree *inTree = (TTree*) inFile->Get("Events");
      if (!isMC) { 
        vtsAnalyser ana(inTree, inTree, inTree, isMC, true, false, false);
        string outFileName = outFileDir+"/rd_"+fileName;
        ana.setOutput(outFileName);
        ana.Loop();
      }
      else {
        if (isGenericMC) {
          vtsAnalyser ana(inTree,isMC,false,false,false,isGenericMC);
          string outFileName = outFileDir+"/mc_"+fileName;
          ana.setOutput(outFileName);
          ana.Loop();
        } else {
          TString hadFileName      = dirName + "/HADAOD/" + fileName;
          TString hadTruthFileName = dirName + "/HADTRUTHAOD/" + fileName;
          TFile *hadFile      = TFile::Open(hadFileName, "READ");
          TFile *hadTruthFile = TFile::Open(hadTruthFileName, "READ");
          TTree *hadTree(0);      
          TTree *hadTruthTree(0);
          if (hadFile != NULL) hadTree = (TTree*) hadFile->Get("Events");
          if (hadTruthFile != NULL) hadTruthTree = (TTree*) hadTruthFile->Get("Events");
          else cout << ">>>>> There dosen't exist hadTruth info. -----> Use NANOTREE and HADTREE <<<<< " << endl;

          vtsAnalyser ana(inTree,hadTree,hadTruthTree,isMC,false,false,false); // you don't need to use hadTree if hadTruthTree isn't NULL
          string outFileName = outFileDir+"/nanotree_"+fileName;
          if (string(inFileName).find("herwig") != std::string::npos) outFileName = outFileDir+"/herwig_nanotree_"+fileName;
          ana.setOutput(outFileName);
          ana.Loop();
        }
      }
      cout << inFileName << " : Analysis Complete " << endl;
    }
  }
}

void vtsAnalyser::Loop() {
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntries();
  /* Events loop */
  for (Long64_t iev=0; iev<nentries; iev++) {
    fChain->GetEntry(iev);
    if (h_fChain) h_fChain->GetEntry(iev);
    if (ht_fChain) ht_fChain->GetEntry(iev);
    if (iev%10000 == 0) cout << iev << "/" << nentries << endl;

    ResetBranch();
    EventSelection();
    m_selectedJet = jetSelection();
    b_ntotjet = nJet;
    b_passedEvent = (b_step >= 4); 

    b_had_start = -1; b_had_end = -1;
    b_jet_start = -1; b_jet_end = -1;
    if (m_isMC && b_step >= 4) { // -1 : No event selection, 0 : PV cut, reco lepton cut and so on, 1~4 : step 1 ~ 4
      MatchingForMC();
      FillTMVATrees();
    }

    m_tree->Fill();
  }
}

void vtsAnalyser::setOutput(std::string outFileName) {
  m_output       = TFile::Open(outFileName.c_str(), "recreate");
  m_tree         = new TTree("event", "event");
  m_hadtrForTMVA = new TTree("MVA_had", "MVA_had");
  m_jettrForTMVA = new TTree("MVA_jet", "MVA_jet");

  MakeBranch();
  SetMVAReader();

  h_nevents    = new TH1D("nevents", "nevents", 1, 0, 1);
  h_genweights = new TH1D("genweight", "genweight", 1, 0, 1);
  h_weights    = new TH1D("weight", "weight", 1, 0, 1);
  h_cutFlow    = new TH1D("cutflow", "cutflow", 11, -0.5, 10.5);
}

void vtsAnalyser::MakeBranchOfHadron(TTree* tr) {
  tr->Branch("Ks_pdgId",        &b_Ks_pdgId,        "Ks_pdgId/I");
  tr->Branch("Ks_nMatched",     &b_Ks_nMatched,     "Ks_nMatched/I");
  tr->Branch("Ks_isFrom",       &b_Ks_isFrom,       "Ks_isFrom/I");
  tr->Branch("Ks_isHadFromTop", &b_Ks_isHadFromTop, "Ks_isHadFromTop/O");
  tr->Branch("Ks_isHadFromW",   &b_Ks_isHadFromW,   "Ks_isHadFromW/O");
  tr->Branch("Ks_isHadFromS",   &b_Ks_isHadFromS,   "Ks_isHadFromS/O");
  tr->Branch("Ks_isHadFromC",   &b_Ks_isHadFromC,   "Ks_isHadFromC/O");
  tr->Branch("Ks_isHadFromB",   &b_Ks_isHadFromB,   "Ks_isHadFromB/O");
  tr->Branch("Ks_d",            &b_Ks_d,            "Ks_d/F");
  tr->Branch("Ks_pt",           &b_Ks_pt,           "Ks_pt/F");
  tr->Branch("Ks_eta",          &b_Ks_eta,          "Ks_eta/F");
  tr->Branch("Ks_phi",          &b_Ks_phi,          "Ks_phi/F");
  tr->Branch("Ks_mass",         &b_Ks_mass,         "Ks_mass/F");
  tr->Branch("Ks_lxy",          &b_Ks_lxy,          "Ks_lxy/F");
  tr->Branch("Ks_lxySig",       &b_Ks_lxySig,       "Ks_lxySig/F");
  tr->Branch("Ks_l3D",          &b_Ks_l3D,          "Ks_l3D/F");
  tr->Branch("Ks_l3DSig",       &b_Ks_l3DSig,       "Ks_l3DSig/F");
  tr->Branch("Ks_legDR",        &b_Ks_legDR,        "Ks_legDR/F");
  tr->Branch("Ks_angleXY",      &b_Ks_angleXY,      "Ks_angleXY/F");
  tr->Branch("Ks_angleXYZ",     &b_Ks_angleXYZ,     "Ks_angleXYZ/F");
  tr->Branch("Ks_chi2",         &b_Ks_chi2,         "Ks_chi2/F");
  tr->Branch("Ks_dca",          &b_Ks_dca,          "Ks_dca/F");
  tr->Branch("Ks_dau1_chi2",    &b_Ks_dau1_chi2,    "Ks_dau1_chi2/F");
  tr->Branch("Ks_dau1_ipsigXY", &b_Ks_dau1_ipsigXY, "Ks_dau1_ipsigXY/F");
  tr->Branch("Ks_dau1_ipsigZ",  &b_Ks_dau1_ipsigZ,  "Ks_dau1_ipsigZ/F");
  tr->Branch("Ks_dau1_pt",      &b_Ks_dau1_pt,      "Ks_dau1_pt/F");
  tr->Branch("Ks_dau2_chi2",    &b_Ks_dau2_chi2,    "Ks_dau2_chi2/F");
  tr->Branch("Ks_dau2_ipsigXY", &b_Ks_dau2_ipsigXY, "Ks_dau2_ipsigXY/F");
  tr->Branch("Ks_dau2_ipsigZ",  &b_Ks_dau2_ipsigZ,  "Ks_dau2_ipsigZ/F");
  tr->Branch("Ks_dau2_pt",      &b_Ks_dau2_pt,      "Ks_dau2_pt/F");
  tr->Branch("Ks_bdt_score",    &b_Ks_bdt_score,    "Ks_bdt_score/F");
}

void vtsAnalyser::MakeBranch() {
  m_jettrForTMVA->Branch("isSJet",         &b_isSJet,         "isSJet/I");
  m_jettrForTMVA->Branch("isBJet",         &b_isBJet,         "isBJet/I");
  m_jettrForTMVA->Branch("isOverlap",      &b_isOverlap,      "isOverlap/I");
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
  m_jettrForTMVA->Branch("hadronFlavour",  &b_hadronFlavour,  "hadronFlavour/F");

  m_jettrForTMVA->Branch("dau_pt",      &b_dau_pt,      "dau_pt[350]/F");
  m_jettrForTMVA->Branch("dau_eta",     &b_dau_eta,     "dau_eta[350]/F");
  m_jettrForTMVA->Branch("dau_phi",     &b_dau_phi,     "dau_phi[350]/F");
  m_jettrForTMVA->Branch("dau_charge",  &b_dau_charge,  "dau_charge[350]/I");

  m_jettrForTMVA->Branch("Jet_bdt_score",   &b_Jet_bdt_score,   "Jet_bdt_score/F");
  m_jettrForTMVA->Branch("JKS_bdt_score",   &b_JKS_bdt_score,   "JKS_bdt_score/F");
  m_jettrForTMVA->Branch("Ks_dr",           &b_Ks_dr,           "Ks_dr/F");
  m_jettrForTMVA->Branch("Ks_x",            &b_Ks_x,            "Ks_x/F");
 
  MakeBranchOfHadron(m_jettrForTMVA);
  MakeBranchOfHadron(m_hadtrForTMVA);

  #define Branch_(type, name, suffix) m_tree->Branch(#name, &(b_##name), #name "/" #suffix);
  #define BranchI(name) Branch_(Int_t, name, I)
  #define BranchF(name) Branch_(Float_t, name, F)
  #define BranchO(name) Branch_(Bool_t, name, O)
  #define BranchV_(type, name) m_tree->Branch(#name, "vector<"#type">", &(b_##name));
  #define BranchVI(name) BranchV_(Int_t, name); 
  #define BranchVF(name) BranchV_(Float_t, name);
  #define BranchVO(name) BranchV_(Bool_t, name);
  #define BranchTLV(name) m_tree->Branch(#name, "TLorentzVector", &(b_##name));

  BranchI(nvertex);     BranchI(step);     BranchI(channel);
  BranchI(njet);        BranchI(nbjet);    BranchF(met); 
  BranchI(jet_start);   BranchI(jet_end);  BranchI(had_start); BranchI(had_end);
  BranchO(passedEvent); BranchI(ntotjet);

  /* weight */
  BranchF(genweight);   BranchF(puweight); 
  BranchF(eleffweight); BranchF(eleffweight_up); BranchF(eleffweight_dn); 
  BranchF(mueffweight); BranchF(mueffweight_up); BranchF(mueffweight_dn); 
  BranchF(btagweight);  BranchVF(csvweights);

  /* trigger */
  BranchF(tri);    BranchF(tri_up);  BranchF(tri_dn);
  BranchO(trig_m); BranchO(trig_m2); BranchO(trig_e); BranchO(trig_mm); BranchO(trig_em); BranchO(trig_ee);

  /* For MatchingForMC() */
  BranchI(tq1_idx);     BranchI(tq1_pdgId);
  BranchI(tq2_idx);     BranchI(tq2_pdgId);
  BranchTLV(tq1_tlv);
  BranchTLV(tq2_tlv);
  BranchI(tq1_matched_jidx); BranchI(tq1_matched_isOverlap); BranchI(tq1_matched_dr); BranchI(tq1_matched_x);
  BranchI(tq2_matched_jidx); BranchI(tq2_matched_isOverlap); BranchI(tq2_matched_dr); BranchI(tq2_matched_x);

}

void vtsAnalyser::ResetBranch() {
  Reset();
  ResetHadTree();  ResetJetTree();

  b_ntotjet = -1;
  b_passedEvent = false;

  /* For MatchingForMC() */
  m_tqMC.clear();       m_wqMC.clear(); 
  m_qjMapForMC.clear();
  m_jetDeltaRs.clear();
  //m_closestRecJetForLep1.clear(); m_closestRecJetForLep2.clear(); 
  b_tq1_idx    = -99;   b_tq2_idx    = -99;
  b_tq1_pdgId  = -99;   b_tq2_pdgId  = -99;
  b_tq1_tlv.SetPtEtaPhiM(0,0,0,0);
  b_tq2_tlv.SetPtEtaPhiM(0,0,0,0);
  b_tq1_matched_jidx = -99; b_tq1_matched_isOverlap = -99; b_tq1_matched_dr = -99; b_tq1_matched_x = -99;
  b_tq2_matched_jidx = -99; b_tq2_matched_isOverlap = -99; b_tq2_matched_dr = -99; b_tq2_matched_x = -99;
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
    auto momIdx = GenPart_genPartIdxMother[i];
    if ((abs(GenPart_pdgId[i]) == 3 ||
         abs(GenPart_pdgId[i]) == 5 ||
         abs(GenPart_pdgId[i]) == 4)
         && ( (GenPart_statusFlags[i] & ( 1 << reco::GenStatusFlags::kIsFirstCopy)) != 0  ) ) {
      if ( (abs(GenPart_pdgId[momIdx]) == 6
           && ( (GenPart_statusFlags[momIdx] & ( 1 << reco::GenStatusFlags::kIsLastCopy)) != 0 ) ) ) { 
        TLorentzVector tlv;
	    tlv.SetPtEtaPhiM(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i]);
        m_tqMC.push_back({(int) i, GenPart_pdgId[i], GenPart_status[i], momIdx, GenPart_pdgId[momIdx], GenPart_status[momIdx], tlv});
      }
      if ( (abs(GenPart_pdgId[momIdx]) == 24
           && (GenPart_status[momIdx] == 22 || GenPart_status[momIdx] == 52)) ) { // not completed yet since not used for now
        TLorentzVector tlv;
        tlv.SetPtEtaPhiM(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i]);
        m_wqMC.push_back({(int) i, GenPart_pdgId[i], GenPart_status[i], GenPart_genPartIdxMother[i], GenPart_pdgId[momIdx], GenPart_status[momIdx], tlv});
      }
    }
  }

  if (m_tqMC.size() < 2) return;

  /* Select quarks originated from t->qW */
  genInfo tq1; genInfo tq2; // tq1 : s-quark , tq2 : b-quark
  if      (abs(m_tqMC[0].pdgId) == 3) { tq1 = m_tqMC[0]; tq2 = m_tqMC[1]; } 
  else if (abs(m_tqMC[0].pdgId) == 5) { tq1 = m_tqMC[1]; tq2 = m_tqMC[0]; }
  TLorentzVector tq1_tlv = tq1.tlv; 
  TLorentzVector tq2_tlv = tq2.tlv;

  b_tq1_idx    = tq1.idx;     b_tq2_idx    = tq2.idx;
  b_tq1_pdgId  = tq1.pdgId;   b_tq2_pdgId  = tq2.pdgId;
  b_tq1_tlv    = tq1.tlv;     b_tq2_tlv    = tq2.tlv;

  /* Select quarks from W->qq (not yet used) */
  genInfo wq1; genInfo wq2;
  TLorentzVector wq1_tlv;
  TLorentzVector wq2_tlv;
  if (m_wqMC.size() != 0) {
    cout << "number of quark from W-boson : " << m_wqMC.size() << endl;
    wq1 = m_wqMC[0];
    wq1_tlv = wq1.tlv;
    if (m_wqMC.size() == 2) {
      wq2 = m_wqMC[1];
      wq2_tlv = wq2.tlv;
    }
  }

  /* tydef jetInfo => {jet idx, jet pt, DeltaR(s,jet), DeltaR(b,jet), DeltaR(lep1,jet), DeltaR(lep2,jet) } */
  /* Gen Partcle & Selected Reco Jet Matching */
  for (unsigned int ij=0; ij<m_selectedJet.size();++ij) {
    unsigned int jidx = m_selectedJet[ij].GetFirstMother();
    TLorentzVector jet_tlv;
    jet_tlv.SetPtEtaPhiM(Jet_pt[jidx], Jet_eta[jidx], Jet_phi[jidx], Jet_mass[jidx]);
    auto dr1  = tq1_tlv.DeltaR(jet_tlv);  auto dr2  = tq2_tlv.DeltaR(jet_tlv);
    auto drl1 = b_lep1.DeltaR(jet_tlv);   auto drl2 = b_lep2.DeltaR(jet_tlv);
    m_jetDeltaRs.push_back({jidx, jet_tlv.Pt(), dr1, dr2, drl1, drl2});
  }

  /* Find closest sel jet to s or b */ 
  sort(m_jetDeltaRs.begin(), m_jetDeltaRs.end(), [] (jetInfo a, jetInfo b) { return (a.pt > b.pt); } ); 
  for (auto j : m_jetDeltaRs) {
      if (j.drsj < CONESIZE && j.drbj < CONESIZE) {
          m_qjMapForMC[j.idx]  = -9;
          b_tq1_matched_isOverlap = 1;
          b_tq2_matched_isOverlap = 1;
          break;
      }
      else if (j.drsj < CONESIZE && j.drbj > CONESIZE) {
          m_qjMapForMC[j.idx] = tq1.pdgId; 
          b_tq1_matched_dr    = m_jetDeltaRs[0].drsj;
          b_tq1_matched_x     = m_jetDeltaRs[0].pt/tq1.tlv.Pt();
          break;
      }
      else if (j.drbj < CONESIZE && j.drsj > CONESIZE) {
          m_qjMapForMC[j.idx] = tq2.pdgId; 
          b_tq2_matched_dr    = m_jetDeltaRs[0].drbj;
          b_tq2_matched_x     = m_jetDeltaRs[0].pt/tq2.tlv.Pt();
          break;
      }
      else cout << "No matched quark to this jet..." << endl;
  }
}

void vtsAnalyser::FillTMVATrees() {
  sort(m_jetDeltaRs.begin(), m_jetDeltaRs.end(), [] (jetInfo a, jetInfo b) { return (a.pt > b.pt); } ); // pT ordering
  auto highestPt = m_jetDeltaRs[0];  auto NhighestPt = m_jetDeltaRs[1];
  auto closestToTq1 = *max_element(m_jetDeltaRs.begin(), m_jetDeltaRs.end(), [] (jetInfo a, jetInfo b) { return (a.drsj < b.drsj); } ); // dR(tq1,jet) ordering
  auto closestToTq2 = *max_element(m_jetDeltaRs.begin(), m_jetDeltaRs.end(), [] (jetInfo a, jetInfo b) { return (a.drbj < b.drbj); } ); // dR(tq2,jet) ordering
  auto closestToLep1 = *max_element(m_jetDeltaRs.begin(), m_jetDeltaRs.end(), [] (jetInfo a, jetInfo b) { return (a.drl1j < b.drl1j); } ); // dR(lep1,jet) ordering
  auto closestToLep2 = *max_element(m_jetDeltaRs.begin(), m_jetDeltaRs.end(), [] (jetInfo a, jetInfo b) { return (a.drl2j < b.drl2j); } ); // dR(lep2,jet) ordering

  /* Fill jet tree */
  b_jet_start =  m_jettrForTMVA->GetEntries();
  for (unsigned int ij=0; ij<m_selectedJet.size();++ij) {
    ResetJetTree();
    SetJetValues(ij);
    if (m_isMC) { IdentifyJet(ij, closestToTq1.idx, closestToTq2.idx); }

    b_isHighest = 0; b_isClosestToLep = 0;
    if      (ij == highestPt.idx     || ij == NhighestPt.idx)    b_isHighest = 1;
    else if (ij == closestToLep1.idx || ij == closestToLep2.idx) b_isClosestToLep = 1;
    b_dr1   = closestToTq1.drsj;  b_dr2   = closestToTq2.drbj;

    TLorentzVector jet_tlv;
    jet_tlv.SetPtEtaPhiM(Jet_pt[ij], Jet_eta[ij], Jet_phi[ij], Jet_mass[ij]);
    auto ks_idx = FindMatchedHadron(jet_tlv);

    ResetHadTree();
    if (ks_idx != -1) {
      SetHadronValues(ks_idx);
      b_Ks_bdt_score = m_hadReader->EvaluateMVA("KS_BDT");
      b_Ks_x = had_pt[ks_idx]/jet_tlv.Pt();
      TLorentzVector had_tlv;
      had_tlv.SetPtEtaPhiM(had_pt[ks_idx], had_eta[ks_idx], had_phi[ks_idx], had_mass[ks_idx]);
      b_Ks_dr = jet_tlv.DeltaR(had_tlv);
    }

    m_jettrForTMVA->Fill();
  }
  b_jet_end = m_jettrForTMVA->GetEntries();

  /* Fill hadron tree */
  b_had_start = m_hadtrForTMVA->GetEntries();
  for (unsigned int j=0; j<nhad; ++j) {
    if (had_pdgId[j] != 310) continue;
    ResetHadTree();
    SetHadronValues(j);
    b_Ks_bdt_score = m_hadReader->EvaluateMVA("KS_BDT"); // calculate based on set hadron values
    m_hadtrForTMVA->Fill();
  }
  b_had_end = m_hadtrForTMVA->GetEntries();
}

void vtsAnalyser::IdentifyJet(unsigned int jetIdx, unsigned int sIdx, unsigned int bIdx) {
  auto jidx = m_jetDeltaRs[jetIdx];
  b_isSJet = 0; b_isBJet = 0;
  b_isOverlap = 0;
  if (jetIdx == sIdx && fabs(jidx.drsj) < CONESIZE) {
     b_tq1_matched_jidx = jetIdx;
     if (abs(m_qjMapForMC[jidx.idx]) == 3)      b_isSJet = 1;
     else if (abs(m_qjMapForMC[jidx.idx]) == 5) b_isBJet = 1;
     else b_isOverlap = 1; 
  }
  else if (jetIdx == bIdx && fabs(jidx.drbj) < CONESIZE) {
     b_tq2_matched_jidx = jetIdx;
     if (abs(m_qjMapForMC[jidx.idx]) == 3)      b_isSJet = 1; 
     else if (abs(m_qjMapForMC[jidx.idx]) == 5) b_isBJet = 1;
     else b_isOverlap = 1;
  }
}

int vtsAnalyser::FindMatchedHadron(TLorentzVector jet_tlv) {
  auto hidx = -1;
  auto maxBDTScore = -1;
  for (unsigned int j=0; j<nhad; ++j) {
    if (had_pdgId[j] != 310) continue;
    TLorentzVector had_tlv;
    had_tlv.SetPtEtaPhiM(had_pt[j], had_eta[j], had_phi[j], had_mass[j]);
    if (jet_tlv.DeltaR(had_tlv) > CONESIZE) continue;

    ResetHadTree();
    SetHadronValues(j);
    auto Ks_bdtScore = m_hadReader->EvaluateMVA("KS_BDT");
    if (maxBDTScore < Ks_bdtScore) {
      maxBDTScore = Ks_bdtScore;
      hidx = j;
    }
  }
  return hidx;
}

void vtsAnalyser::SetJetValues(int i) {
  /* Save jet information */
  b_cmult = (float) jetID_cmult[i];   b_nmult = (float) jetID_nmult[i];
  b_pt    = Jet_pt[i];  b_eta   = Jet_eta[i];  b_phi  = Jet_phi[i];  b_mass = Jet_mass[i];
  b_c_x1  = jetID_cpt1[i]/Jet_pt[i]; b_c_x2  = jetID_cpt2[i]/Jet_pt[i]; b_c_x3 = jetID_cpt3[i]/Jet_pt[i];
  b_n_x1  = jetID_npt1[i]/Jet_pt[i]; b_n_x2  = jetID_npt2[i]/Jet_pt[i]; b_n_x3 = jetID_npt3[i]/Jet_pt[i];
  b_axis1 = jetID_axis1[i];  b_axis2 = jetID_axis2[i];   b_ptD  = jetID_ptD[i];
  b_area  = Jet_area[i];     b_CSVV2 = Jet_btagCSVV2[i]; b_hadronFlavour = Jet_hadronFlavour[i];

  /* Save jet daughter information */
  int ia = 0;
  for (auto didx = jetID_dauIdx1[i]; didx < jetID_dauIdx2[i]; ++didx) {
    b_dau_pt[ia]     = jetDau_pt[didx];
    b_dau_eta[ia]    = jetDau_eta[didx];
    b_dau_phi[ia]    = jetDau_phi[didx];
    b_dau_charge[ia] = jetDau_charge[didx];
    ia += 1;
  }

  /* Save tmva value */
  //b_Jet_bdt_score = m_jetReader->EvaluateMVA("Jet_BDT_highest");
  //b_JKS_bdt_score = m_jksReader->EvaluateMVA("JKS_BDT_highest");
}

void vtsAnalyser::SetHadronValues(int i) {
  if (!m_isGenericMC) {
    b_Ks_isHadFromTop = hadTruth_isHadFromTop[i];
    b_Ks_isHadFromW = hadTruth_isHadFromW[i];
    b_Ks_isHadFromS = hadTruth_isHadFromS[i];
    b_Ks_isHadFromC = hadTruth_isHadFromC[i];
    b_Ks_isHadFromB = hadTruth_isHadFromB[i];
    b_Ks_isFrom = hadTruth_isHadFromTsb[i];
    b_Ks_nMatched = hadTruth_nMatched[i];
  }
  b_Ks_d = GetD(had_pt[i], had_eta[i], had_phi[i], had_mass[i], had_x[i], had_y[i], had_z[i]);
  b_Ks_pt = had_pt[i];
  b_Ks_eta = had_eta[i];
  b_Ks_phi = had_phi[i];
  b_Ks_mass = had_mass[i];
  b_Ks_lxy = had_lxy[i];
  b_Ks_lxySig = had_lxy[i]/had_lxyErr[i];
  b_Ks_angleXY = had_angleXY[i];
  b_Ks_angleXYZ = had_angleXYZ[i];
  b_Ks_chi2 = had_chi2[i];
  b_Ks_dca = had_dca[i];
  b_Ks_l3D = had_l3D[i];
  b_Ks_l3DSig = had_l3D[i] / had_l3DErr[i];
  b_Ks_legDR = had_legDR[i];
  b_Ks_pdgId = had_pdgId[i];
  b_Ks_dau1_chi2 = had_dau1_chi2[i];
  b_Ks_dau1_ipsigXY = had_dau1_ipsigXY[i];
  b_Ks_dau1_ipsigZ = had_dau1_ipsigZ[i];
  b_Ks_dau1_pt = had_dau1_pt[i];
  b_Ks_dau2_chi2 = had_dau2_chi2[i];
  b_Ks_dau2_ipsigXY = had_dau2_ipsigXY[i];
  b_Ks_dau2_ipsigZ = had_dau2_ipsigZ[i];
  b_Ks_dau2_pt = had_dau2_pt[i];
}

void vtsAnalyser::SetMVAReader() {
  #define TMVABranch_(reader, name) reader->AddVariable(#name, &(b_##name));
  #define hadTMVABranch(name) TMVABranch_(m_hadReader,name);
  #define jetTMVABranch(name) TMVABranch_(m_jetReader,name);
  #define jksTMVABranch(name) TMVABranch_(m_jksReader,name);

  m_hadReader = new TMVA::Reader();            
  hadTMVABranch(Ks_d); hadTMVABranch(Ks_pt); hadTMVABranch(Ks_eta); hadTMVABranch(Ks_phi);
  hadTMVABranch(Ks_lxy); hadTMVABranch(Ks_lxySig);
  hadTMVABranch(Ks_l3D); hadTMVABranch(Ks_l3DSig);
  hadTMVABranch(Ks_legDR); hadTMVABranch(Ks_angleXY); hadTMVABranch(Ks_angleXYZ);
  hadTMVABranch(Ks_chi2); hadTMVABranch(Ks_dca);
  hadTMVABranch(Ks_dau1_chi2); hadTMVABranch(Ks_dau1_ipsigXY); hadTMVABranch(Ks_dau1_ipsigZ); hadTMVABranch(Ks_dau1_pt);
  hadTMVABranch(Ks_dau2_chi2); hadTMVABranch(Ks_dau2_ipsigXY); hadTMVABranch(Ks_dau2_ipsigZ); hadTMVABranch(Ks_dau2_pt);
  m_hadReader->BookMVA("KS_BDT", "/cms/ldap_home/tt8888tt/nanoAOD_tmp/src/nano/analysis/test/vts/tmva/pp_real_vs_fake/weights/vts_dR_04_Had_BDT.weights.xml");

  m_jetReader = new TMVA::Reader();
  jetTMVABranch(pt);  jetTMVABranch(eta);  jetTMVABranch(phi);  jetTMVABranch(mass);
  jetTMVABranch(c_x1);  jetTMVABranch(c_x2);  jetTMVABranch(c_x3);
  jetTMVABranch(n_x1);  jetTMVABranch(n_x2);  jetTMVABranch(n_x3);
  jetTMVABranch(cmult); jetTMVABranch(nmult);
  jetTMVABranch(axis1); jetTMVABranch(axis2);
  jetTMVABranch(ptD); jetTMVABranch(area); jetTMVABranch(CSVV2);
  //m_jetReader->BookMVA("Jet_BDT_highest", "/cms/ldap_home/wjjang/wj_nanoAOD_CMSSW_9_4_4/src/nano/analysis/test/vts/tmva/dataset/Jet/pp_combined_J_BDT_highest/weights/vts_dR_04_Jet_BDT.weights.xml");

  m_jksReader = new TMVA::Reader();
  jksTMVABranch(pt); jksTMVABranch(eta); jksTMVABranch(phi); jksTMVABranch(mass);
  jksTMVABranch(c_x1); jksTMVABranch(c_x2); jksTMVABranch(c_x3);
  jksTMVABranch(n_x1); jksTMVABranch(n_x2); jksTMVABranch(n_x3);
  jksTMVABranch(cmult); jksTMVABranch(nmult);
  jksTMVABranch(axis1); jksTMVABranch(axis2);
  jksTMVABranch(ptD);  jksTMVABranch(area);  jksTMVABranch(CSVV2);

  jksTMVABranch(Ks_d); jksTMVABranch(Ks_pt); jksTMVABranch(Ks_eta); jksTMVABranch(Ks_phi);
  jksTMVABranch(Ks_lxy); jksTMVABranch(Ks_lxySig);
  jksTMVABranch(Ks_l3D); jksTMVABranch(Ks_l3DSig);
  jksTMVABranch(Ks_legDR); jksTMVABranch(Ks_angleXY); jksTMVABranch(Ks_angleXYZ);
  jksTMVABranch(Ks_chi2); jksTMVABranch(Ks_dca);
  jksTMVABranch(Ks_dau1_chi2); jksTMVABranch(Ks_dau1_ipsigXY); jksTMVABranch(Ks_dau1_ipsigZ); jksTMVABranch(Ks_dau1_pt);
  jksTMVABranch(Ks_dau2_chi2); jksTMVABranch(Ks_dau2_ipsigXY); jksTMVABranch(Ks_dau2_ipsigZ); jksTMVABranch(Ks_dau2_pt);

  jksTMVABranch(Ks_bdt_score);
  jksTMVABranch(Ks_x);
  //m_jksReader->BookMVA("JKs_BDT_highest", "/cms/ldap_home/wjjang/wj_nanoAOD_CMSSW_9_4_4/src/nano/analysis/test/vts/tmva/dataset/JKS/pp_combined_JKs_BDT_highest/weights/vts_dR_04_Jet_BDT.weights.xml");
}

void vtsAnalyser::ResetJetTree() {
  b_Jet_bdt_score   = -99;
  b_JKS_bdt_score   = -99;

  /* temp initialization of array */
  for (unsigned int i = 0; i < m_jetDauArrSize; ++i) {
    b_dau_pt[i]     = -9.;
    b_dau_eta[i]    = -9.;
    b_dau_phi[i]    = -9.;
    b_dau_charge[i] = -9;
  }

  b_Ks_dr = -99;
  b_Ks_x = -99;
}

void vtsAnalyser::ResetHadTree() {
  b_Ks_idx          = -99;   b_Ks_nMatched     = -99;   b_Ks_isFrom       = -99; 
  b_Ks_isHadFromTop = false; b_Ks_isHadFromW   = false; b_Ks_isHadFromS   = false; b_Ks_isHadFromC   = false; b_Ks_isHadFromB   = false; 
  b_Ks_d            = -99;   b_Ks_pt           = -99;   b_Ks_eta          = -99;   b_Ks_phi          = -99;   b_Ks_mass         = -99; 
  b_Ks_lxy          = -99;   b_Ks_lxySig       = -99;   b_Ks_l3D          = -99;   b_Ks_l3DSig       = -99;   b_Ks_legDR        = -99; 
  b_Ks_angleXY      = -99;   b_Ks_angleXYZ     = -99;   b_Ks_chi2         = -99;   b_Ks_dca          = -99; 
  b_Ks_dau1_chi2    = -99;   b_Ks_dau1_ipsigXY = -99;   b_Ks_dau1_ipsigZ  = -99;   b_Ks_dau1_pt      = -99;
  b_Ks_dau2_chi2    = -99;   b_Ks_dau2_ipsigXY = -99;   b_Ks_dau2_ipsigZ  = -99;   b_Ks_dau2_pt      = -99; 
  b_Ks_bdt_score    = -99;
}

