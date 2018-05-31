#include <iostream>
#include <vector>
#include <TLorentzVector.h>

#include "nano/analysis/interface/semiLepTopAnalysis.h"
#include "nano/analysis/interface/hadAnalysis.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;

class slCalibAnalysis : public semiLepTopAnalysis
{
private:
  //functions
  void MakeBranch(TTree* t);
  void HadronAnalysis();

  TLorentzVector b_had_tlv, b_jet_tlv;
  float b_x_had;
  float b_btagCSVV2_Jet;
  float b_lxy_had, b_lxySig_had, b_angleXY_had, b_angleXYZ_had, b_chi2_had, b_dca_had, b_dr_had;
  int b_pdgId_had, b_nConstituents_Jet;

public:
  void setOutput(std::string outputName);

  slCalibAnalysis(TTree *tree=0, Bool_t isMC = false, Bool_t sle = false, Bool_t slm = false) : semiLepTopAnalysis(tree, isMC, sle, slm) {}
  ~slCalibAnalysis() {}
  virtual void Loop();
};

void slCalibAnalysis::Loop() {
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntries();

  // Events loop
  for (Long64_t iev=0; iev<nentries; iev++) {
    Long64_t entry = LoadTree(iev);
    fChain->GetEntry(entry);
    if (iev%10000 == 0) cout << iev << "/" << nentries << endl;

    Reset();
    int PassedStep = EventSelection();
    if (PassedStep >= 0) {
      HadronAnalysis();
      m_tree->Fill();
    } else {
      m_tree->Fill();
    }
  }
}

void slCalibAnalysis::setOutput(std::string outFileName)
{
  m_output = TFile::Open(outFileName.c_str(), "recreate");
  m_tree = new TTree("event", "event");

  MakeBranch(m_tree);

  h_nevents = new TH1D("nevents", "nevents", 1, 0, 1);
  h_genweights = new TH1D("genweight", "genweight", 1, 0, 1);
  h_weights = new TH1D("weight", "weight", 1, 0, 1);
  h_cutFlow = new TH1D("cutflow", "cutflow", 11, -0.5, 10.5);
  h_cutFlowEl = new TH1D("cutflowEl", ";step;# events", 11, -0.5, 10.5);
  h_cutFlowMu = new TH1D("cutflowMu", ";step;# events", 11, -0.5, 10.5);
}

void slCalibAnalysis::MakeBranch(TTree* t)
{
  t->Branch("channel", &b_channel, "channel/I");
  t->Branch("step", &b_step, "step/I");
  t->Branch("njet", &b_njet, "njet/I");
  t->Branch("met", &b_met, "met/F");

  t->Branch("dr_had", &b_dr_had, "dr_had/F"); // distance between hadron and jet-center
  t->Branch("lxy_had", &b_lxy_had, "lxy_had/F");
  t->Branch("lxySig_had", &b_lxySig_had, "lxySig_had/F");
  t->Branch("angleXY_had", &b_angleXY_had, "angleXY_had/F");
  t->Branch("angleXYZ_had", &b_angleXYZ_had, "angleXYZ_had/F");
  t->Branch("chi2_had", &b_chi2_had, "chi2_had/F");
  t->Branch("dca_had", &b_dca_had, "dca_had/F");

  t->Branch("had_tlv", &b_had_tlv);
  t->Branch("jet_tlv", &b_jet_tlv);
  t->Branch("x_had", &b_had_tlv, "x_had/F");
  
  t->Branch("btagCSVV2_Jet", &b_btagCSVV2_Jet, "btagCSVV2_Jet/F");
  t->Branch("nConstituents_Jet", &b_nConstituents_Jet, "nConstituents_Jet/I");
}

void slCalibAnalysis::HadronAnalysis()
{
  
  std::vector<hadAnalysis::HadStat> JetCollection;
  std::vector<hadAnalysis::HadStat> Had;
  hadAnalysis::HadStat Stat;

  for (unsigned j=0; j < nJet; ++j) {
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
	//        if (m_isMC) Stat.label = qjMapForMC_[j];
        Stat.isHadJetMatched = true;
        Stat.dr = jet_tlv.DeltaR(had_tlv);
        Stat.jetIdx = j;
        Had.push_back(Stat);
      }
    }

    if (Had.size() != 0 ) {
      if (Had.size() > 1) sort(Had.begin(), Had.end(), [](hadAnalysis::HadStat a, hadAnalysis::HadStat b) {return a.x > b.x;}); // pick hadron with highest x in the jet
      JetCollection.push_back(Had[0]);
    }    
  }

  if (JetCollection.size() != 0 ) {
    if (JetCollection.size() > 1) sort(JetCollection.begin(), JetCollection.end(), [](hadAnalysis::HadStat a, hadAnalysis::HadStat b) {return a.x > b.x;}); // pick jet-hadron pair with highest x
    auto idx = JetCollection[0].idx;
    auto jidx = JetCollection[0].jetIdx;
    b_had_tlv.SetPtEtaPhiM(had_pt[idx], had_eta[idx], had_phi[idx], had_mass[idx]);
    b_jet_tlv.SetPtEtaPhiM(Jet_pt[jidx], Jet_eta[jidx], Jet_phi[jidx], Jet_mass[jidx]);

    b_x_had = JetCollection[0].x;
    // b_isFrom_had = JetCollection[0].label;  // -99 : event that can't pass till step4(jet selection) or there is no matching between had and jet, -9 : there is t->qW in the event,but not matched to recoJet, 0 : there is no t->qW in the event (if no t->s and no matching between had-jet, then the event would be -99), +-3 : hadron is from t->sW, +-5 : hadron is from t->bW
    // b_isHadJetMatched_had = JetCollection[0].isHadJetMatched;
    // b_dr_had = JetCollection[0].dr;

    // b_d_had = GetD(had_pt[idx], had_eta[idx], had_phi[idx], had_mass[idx], had_x[idx], had_y[idx], had_z[idx]);
    b_lxy_had = had_lxy[idx];
    b_lxySig_had = had_lxy[idx]/had_lxyErr[idx];
    b_angleXY_had = had_angleXY[idx];
    b_angleXYZ_had = had_angleXYZ[idx];
    b_chi2_had = had_chi2[idx];
    b_dca_had = had_dca[idx];
    // b_l3D_had = had_l3D[idx];
    // b_l3DSig_had = had_l3D[idx]/had_l3DErr[idx];
    // b_legDR_had = had_legDR[idx];
    b_pdgId_had = had_pdgId[idx];
    // b_dau1_chi2_had = had_dau1_chi2[idx]; 
    // b_dau1_ipsigXY_had = had_dau1_ipsigXY[idx]; 
    // b_dau1_ipsigZ_had = had_dau1_ipsigZ[idx]; 
    // b_dau1_pt_had = had_dau1_pt[idx];
    // b_dau2_chi2_had = had_dau2_chi2[idx]; 
    // b_dau2_ipsigXY_had = had_dau1_ipsigXY[idx]; 
    // b_dau2_ipsigZ_had = had_dau1_ipsigZ[idx]; 
    // b_dau2_pt_had = had_dau1_pt[idx];

    b_btagCSVV2_Jet = Jet_btagCSVV2[jidx];
    // b_btagCMVA_Jet = Jet_btagCMVA[jidx]; 
    // b_btagDeepB_Jet = Jet_btagDeepB[jidx];
    // b_btagDeepC_Jet = Jet_btagDeepC[jidx];
    // b_area_Jet = Jet_area[jidx]; 
    b_nConstituents_Jet = Jet_nConstituents[jidx]; 
    // b_nElectrons_Jet = Jet_nElectrons[jidx];
    // b_nMuons_Jet = Jet_nMuons[jidx];
  }
}

int main(int argc, char* argv[])
{
  string hostDir = "root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/";

  if (argc <= 1) {
    cout << "no input file is specified. running with default file." << endl;
    auto inFile = TFile::Open("/xrootd/store/group/nanoAOD/run2_2016v4/tsw/nanoAOD_1.root", "READ");
    auto inTree = (TTree*) inFile->Get("Events");
    slCalibAnalysis ana(inTree,true,false,false);
    ana.setOutput("nanotree.root");
    ana.Loop();
  } else {
    string inName = string(argv[1]);
    string outName = string(argv[2]);

    // string jobName    = string(argv[1]);
    // string sampleName = string(argv[2]);

    // temp
    Bool_t isMC = false;
    std::string temp = argv[1];
    Size_t found = temp.find("run");
    if (found == std::string::npos) isMC = true;

    // string outFileDir = hostDir+getenv("USER")+"/"+jobName+"/"+sampleName;
    // for (Int_t i = 3; i < argc; i++) {
    //   auto inFileName = argv[i];
    //   TFile *inFile = TFile::Open(inFileName, "READ");
    TFile *inFile = TFile::Open(inName.c_str(), "READ");
    TTree *inTree = (TTree*) inFile->Get("Events");
    slCalibAnalysis ana(inTree,isMC,false,false);
    // string outFileName = outFileDir+"/nanotree_"+to_string(i-3)+".root";
    ana.setOutput(outName.c_str());
    ana.Loop();
    // }
  }
}
