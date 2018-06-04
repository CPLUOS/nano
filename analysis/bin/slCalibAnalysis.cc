#include <iostream>
#include <vector>
#include <TLorentzVector.h>

#include "nano/analysis/interface/semiLepTopAnalysis.h"
#include "nano/analysis/interface/hadAnalysis.h"
#include "nano/analysis/interface/HadTruthEvents.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

using namespace std;

class slCalibAnalysis : public semiLepTopAnalysis
{
private:
  //functions
  void MakeBranch(TTree* t);
  void HadronAnalysis();

  TLorentzVector m_had_tlv, m_jet_tlv;
  float m_had_x;
  float m_jet_btagCSVV2;
  float m_had_lxy, m_had_lxySig, m_had_angleXY, m_had_angleXYZ, m_had_chi2, m_had_dca, m_had_dr;
  int m_had_pdgId, m_jet_nConstituents;
  TTree *hadt = 0;
  HadTruthEvents had;

public:
  void setOutput(std::string outputName);
  void setHadInput(TTree *t) { hadt = t; had.Init(hadt); }

  slCalibAnalysis(TTree *tree=0, Bool_t isMC = false, Bool_t sle = false, Bool_t slm = false) : semiLepTopAnalysis(tree, isMC, sle, slm) {}
  ~slCalibAnalysis() {}
  virtual void Loop();
};

void slCalibAnalysis::Loop() {
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntries();

  // Events loop
  for (Long64_t iev=0; iev<nentries; iev++) {
    fChain->GetEntry(iev);
    if (hadt) hadt->GetEntry(iev);
    if (iev%10000 == 0) cout << iev << "/" << nentries << " : " << fChain->GetCurrentFile()->GetName()
			     << (hadt ? hadt->GetCurrentFile()->GetName() : "")
			     << endl;
    if (hadt) {
      if (had.event != event || had.run != run || had.luminosityBlock != luminosityBlock) {
	std::cout << "Bad sync! " << event << " " << had.event << " " << run << " " << had.run
		  << " " << luminosityBlock << " " << had.luminosityBlock
		  << std::endl;
	exit(1);
      }
    }

    Reset();
    int PassedStep = EventSelection();
    if (PassedStep >= 3) {
      HadronAnalysis();
      m_tree->Fill();
    } else {
      // m_tree->Fill();
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

  t->Branch("had_dr", &m_had_dr, "had_dr/F"); // distance between hadron and jet-center
  t->Branch("had_lxy", &m_had_lxy, "had_lxy/F");
  t->Branch("had_lxySig", &m_had_lxySig, "had_lxySig/F");
  t->Branch("had_angleXY", &m_had_angleXY, "had_angleXY/F");
  t->Branch("had_angleXYZ", &m_had_angleXYZ, "had_angleXYZ/F");
  t->Branch("had_chi2", &m_had_chi2, "had_chi2/F");
  t->Branch("had_dca", &m_had_dca, "had_dca/F");

  t->Branch("had_tlv", &m_had_tlv);
  t->Branch("jet_tlv", &m_jet_tlv);
  t->Branch("had_x", &m_had_x, "had_x/F");
  
  t->Branch("jet_btagCSVV2", &m_jet_btagCSVV2, "jet_btagCSVV2/F");
  t->Branch("jet_nConstituents", &m_jet_nConstituents, "jet_nConstituents/I");
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
    m_had_tlv.SetPtEtaPhiM(had_pt[idx], had_eta[idx], had_phi[idx], had_mass[idx]);
    m_jet_tlv.SetPtEtaPhiM(Jet_pt[jidx], Jet_eta[jidx], Jet_phi[jidx], Jet_mass[jidx]);

    m_had_x = JetCollection[0].x;
    // m_had_isFrom = JetCollection[0].label;  // -99 : event that can't pass till step4(jet selection) or there is no matching between had and jet, -9 : there is t->qW in the event,but not matched to recoJet, 0 : there is no t->qW in the event (if no t->s and no matching between had-jet, then the event would be -99), +-3 : hadron is from t->sW, +-5 : hadron is from t->bW
    // m_had_isHadJetMatched = JetCollection[0].isHadJetMatched;
    // m_had_dr = JetCollection[0].dr;

    // m_had_d = GetD(had_pt[idx], had_eta[idx], had_phi[idx], had_mass[idx], had_x[idx], had_y[idx], had_z[idx]);
    m_had_lxy = had_lxy[idx];
    m_had_lxySig = had_lxy[idx]/had_lxyErr[idx];
    m_had_angleXY = had_angleXY[idx];
    m_had_angleXYZ = had_angleXYZ[idx];
    m_had_chi2 = had_chi2[idx];
    m_had_dca = had_dca[idx];
    // m_had_l3D = had_l3D[idx];
    // m_had_l3DSig = had_l3D[idx]/had_l3DErr[idx];
    // m_had_legDR = had_legDR[idx];
    m_had_pdgId = had_pdgId[idx];
    // m_dau1_had_chi2 = had_dau1_chi2[idx]; 
    // m_dau1_had_ipsigXY = had_dau1_ipsigXY[idx]; 
    // m_dau1_had_ipsigZ = had_dau1_ipsigZ[idx]; 
    // m_dau1_had_pt = had_dau1_pt[idx];
    // m_dau2_had_chi2 = had_dau2_chi2[idx]; 
    // m_dau2_had_ipsigXY = had_dau1_ipsigXY[idx]; 
    // m_dau2_had_ipsigZ = had_dau1_ipsigZ[idx]; 
    // m_dau2_had_pt = had_dau1_pt[idx];

    m_jet_btagCSVV2 = Jet_btagCSVV2[jidx];
    // m_jet_btagCMVA = Jet_btagCMVA[jidx]; 
    // m_jet_btagDeepB = Jet_btagDeepB[jidx];
    // m_jet_btagDeepC = Jet_btagDeepC[jidx];
    // m_jet_area = Jet_area[jidx]; 
    m_jet_nConstituents = Jet_nConstituents[jidx]; 
    // m_jet_nElectrons = Jet_nElectrons[jidx];
    // m_jet_nMuons = Jet_nMuons[jidx];
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
  } else if (argc > 4) {
    cout << "Usage: ./" << argv[0] << " <input file glob> <output filename> [<had file glob>]" << endl;
  } else {
    auto inName = argv[1];
    auto outName = argv[2];

    Bool_t isMC = false;
    auto temp = string(inName);
    auto found = temp.find("run");
    if (found == std::string::npos) isMC = true;

    TChain inTree{"Events"};
    inTree.Add(inName);
    TChain *hadTree = nullptr;
    if (argc > 3) {
      hadTree = new TChain{"Events"};
      hadTree->Add(argv[3]);
    }
    cout << "Running " << inTree.GetEntries() << " entries from " << inTree.GetListOfFiles()->GetEntries() << " files: " << inName << endl;
    slCalibAnalysis ana(&inTree,isMC,false,false);
    if (hadTree) ana.setHadInput(hadTree);
    ana.setOutput(outName);
    ana.Loop();
  }
}
