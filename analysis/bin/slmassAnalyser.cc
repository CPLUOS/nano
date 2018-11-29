#include "nano/analysis/interface/slmassAnalyser.h"
#include "nano/analysis/interface/hadAnalyser.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <cstdlib>
using namespace std;

slmassAnalyser::slmassAnalyser(TTree *tree, TTree *had, TTree *hadTruth, Bool_t isMC,  Bool_t sle, Bool_t slm) : topEventSelectionSL(tree, had, hadTruth, isMC, sle, slm)
{
SetCutValues();
}
//yebbi
slmassAnalyser::~slmassAnalyser()
{
}

int slmassAnalyser::SetCutValues(){
    topEventSelectionSL::SetCutValues();
    
    cut_ElectronPt = 34; 
    cut_ElectronEta = 2.1;
    cut_ElectronIDType = Electron_cutBased;
    cut_ElectronIDCut = 3;
    cut_ElectronSCEtaLower = 1.4442;
    cut_ElectronSCEtaUpper = 1.566;
    cut_ElectronRelIso03All = 0.1;
    
    cut_MuonIDType = Muon_tightId;
    cut_MuonPt = 26; 
    cut_MuonEta = 2.4;
    cut_MuonRelIso04All = 0.15;
    
    cut_VetoElectronPt = 15; 
    cut_VetoElectronEta = 2.4;
    cut_VetoElectronIDType = Electron_cutBased;
    cut_VetoElectronIDCut = 1;
    cut_VetoElectronSCEtaLower = 1.4442;
    cut_VetoElectronSCEtaUpper = 1.566;
    cut_VetoElectronRelIso03All = 0.26;
    
    cut_VetoMuonIDType = NULL;
    cut_VetoMuonPt = 15; 
    cut_VetoMuonEta = 2.4;
    cut_VetoMuonRelIso04All = 0.26;
    
    cut_JetID = 1;
    cut_JetPt = 30; 
    cut_JetEta = 2.4;
    cut_JetConeSizeOverlap = 0.4;
    
    cut_BJetID = 1;
    cut_BJetPt = 30; 
    cut_BJetEta = 2.4;
    cut_BJetConeSizeOverlap = 0.4;
    cut_BJetTypeBTag = Jet_btagCSVV2;
    cut_BJetBTagCut = 0.8484;
    
    return 0;

}



void slmassAnalyser::Loop() {
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "# of input entries: " << nentries << std::endl;
  // Events loop
  for (Long64_t iev=0; iev<nentries; iev++) {
    resetBranch();
    fChain->GetEntry(iev);
    int keep = EventSelection();
    if (keep != 0) {
        cmesonSelection();
        m_tree->Fill();
    }
  }

}

int main(int argc, char* argv[]) {
  string env = getenv("CMSSW_BASE");
  string username = getenv("USER");

  if (argc != 1) {
    std::string dirName = "root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/" + username + "/nanoAOD/" + std::string(argv[1]) + "/" + std::string(argv[2]);
    std::string temp = argv[2];
    
    Bool_t isSL_e = false;
    Size_t found_SL_e = temp.find("SingleElectron");
    if (found_SL_e != std::string::npos) isSL_e = true;
    
    Bool_t isSL_m = false;
    Size_t found_SL_m = temp.find("SingleMuon");
    if (found_SL_m != std::string::npos) isSL_m = true;

    Bool_t isMC = false;
    Size_t found = temp.find("Run");
    if (found == std::string::npos) isMC = true;

    for(Int_t i = 3; i < argc; i++) {
        TFile *f = TFile::Open(argv[i], "read");
        TTree *tree;                  
        f->GetObject("Events", tree);
    
        temp = argv[i];   
        found = temp.find_last_of('/');
        std::string outPutName = dirName + temp.substr(found);
        slmassAnalyser t(tree, tree, 0, isMC, isSL_e, isSL_m);
    
        t.setOutput(outPutName);
        t.Loop();
    }
  }
  else {
    //TFile *f = TFile::Open("/xrootd/store/group/nanoAOD/run2_2016v5/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180607_115926/0000/nanoAOD_256.root", "read");
    TFile *f = TFile::Open("/xrootd/store/group/nanoAOD/run2_2016v5/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180610_143635/0000/nanoAOD_256.root", "read");
    //TFile *f = TFile::Open("/xrootd/store/group/nanoAOD/run2_2016v5/SingleMuon/Run2016B-07Aug17_ver2-v1/180607_085033/0000/nanoAOD_430.root");
    //TFile *f = TFile::Open("/xrootd/store/group/nanoAOD/run2_2016v4/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180430_152541/0000/nanoAOD_256.root", "read");
    //TFile *f = TFile::Open("/xrootd/store/group/nanoAOD/run2_2016v5/WW_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180610_143956/0000/nanoAOD_2.ro 
    
    TTree *tree;
    f->GetObject("Events", tree);
    cout <<"out"<<std::endl;
    slmassAnalyser t(tree, tree, 0, true, false, false);
    t.setOutput("sltest.root");
    t.Loop();
  }
  return 0;
}

void slmassAnalyser::setOutput(std::string outputName) {
  int nTrial = 0;
  while (m_output == NULL ){
  m_output = TFile::Open(outputName.c_str(), "recreate");
    if(m_output == NULL) sleep(15*1000);
    if(nTrial++ >5) break;
  }

   
  m_tree = new TTree("event", "event");
  MakeBranch(m_tree);

  bdtg = new TMVA::Reader();
  bdtg->AddVariable("cme_l3DE", &tmp_l3DE);
  bdtg->AddVariable("cme_jetDR", &tmp_jetDR);
  bdtg->AddVariable("cme_legDR", &tmp_legDR);
  bdtg->AddVariable("cme_dca", &tmp_dca);
  bdtg->AddVariable("cme_chi2", &tmp_chi2);
  bdtg->AddVariable("cme_jet_btagCSVV2", &tmp_jet_btagCSVV2);
  bdtg->AddSpectator("cme_mass", &tmp_mass);
  bdtg->BookMVA("BDTG", "./Had_928.weights.xml");

  h_nevents = new TH1D("nevents", "nevents", 1, 0, 1);
  h_genweights = new TH1D("genweight", "genweight", 1, 0, 1);
  h_weights = new TH1D("weight", "weight", 1, 0, 1);
  h_cutFlow = new TH1D("cutflow", "cutflow", 11, -0.5, 10.5);
  h_cutFlowEl = new TH1D("cutflowEl", "cutflow", 11, -0.5, 10.5);
  h_cutFlowMu = new TH1D("cutflowMu", "cutflow", 11, -0.5, 10.5);
}


void slmassAnalyser::MakeBranch(TTree* t) {
  
  t->Branch("nvertex", &b_nvertex, "nvertex/I");
  t->Branch("step", &b_step, "step/I");
  t->Branch("channel", &b_channel, "channel/I");
  t->Branch("njet", &b_njet, "njet/I");
  t->Branch("nbjet", &b_nbjet, "nbjet/I");
  t->Branch("lep_pid", &b_lep_pid, "lep_pid/I");
  t->Branch("lep", "TLorentzVector", &b_lep); 
  t->Branch("had", "TLorentzVector", &b_had); 
  t->Branch("mass_sum", &b_sum, "mass_sum/F"); 
  t->Branch("tri", &b_tri, "tri/F");
  t->Branch("tri_up", &b_tri_up, "tri_up/F");
  t->Branch("tri_dn", &b_tri_dn, "tri_dn/F");
  t->Branch("met", &b_met, "met/F");
  t->Branch("nhad", &nhad, "nhad/I");
  t->Branch("weight", &b_weight, "weight/F");
  t->Branch("puweight", &b_puweight, "puweight/F");
  t->Branch("genweight", &b_genweight, "genweight/F");
  t->Branch("csvweight", "std::vector<float>", &b_csvweights);
  t->Branch("btagweight", &b_btagweight, "btagweight/F");
  t->Branch("mueffweight", &b_mueffweight, "mueffweight/F");
  t->Branch("eleffweight", &b_eleffweight, "eleffweight/F");
  t->Branch("PV_npvs", &PV_npvs, "PV_npvs/I");
  t->Branch("trig_m", &b_trig_m, "trig_m/I");
  t->Branch("trig_e", &b_trig_e, "trig_e/I");
  
//t->Branch("had_vecSumDMLep", "TLoremtzVector", &b_had_vecSumDMLep);
//t->Branch("had_sumM", &b_had_sumM, "had_sumM/F");
  //t->Branch("cme_nMatched", &b_cme_nMatched, "cme_nMatched/I");
  //t->Branch("cme_nTrueDau", &b_cme_nTrueDau, "cme_nTrueDau/I");
  t->Branch("cme_diffMass", &b_cme_diffMass, "cme_diffMass/F");
  t->Branch("cme_tmva", &b_cme_tmva, "cme_tmva/F");
  //t->Branch("cme_tmva_Jpsi", &b_cme_tmva_Jpsi, "cme_tmva_Jpsi/F");
  t->Branch("cme_mass", &b_cme_mass, "cme_mass/F");
  t->Branch("cme_pdgId", &b_cme_pdgId, "cme_pdgId/I");
  //t->Branch("cme_jet_btagDeepB", &b_cme_jet_btagDeepB,"cme_jet_btagDeepB/F");
  //t->Branch("cme_jet_btagDeepC", &b_cme_jet_btagDeepC,"cme_jet_btagDeepC/F");
  //t->Branch("cme_jet_btagCMVA", &b_cme_jet_btagCMVA,"cme_jet_btagCMVA/F");
  t->Branch("cme_jet_btagCSVV2", &b_cme_jet_btagCSVV2,"cme_jet_btagCSVV2/F");
  t->Branch("cme_l3DE", &b_cme_l3DE, "cme_l3DE/F");
  //t->Branch("cme_l3D", &b_cme_l3D,"cme_l3D/F");
  //t->Branch("cme_l3DErr", &b_cme_l3DErr, "cme_l3DErr/F");
  //t->Branch("cme_lxy", &b_cme_lxy, "cme_lxy/F");
  //t->Branch("cme_lxyE", &b_cme_lxyE,"cme_lxyE/F");
  t->Branch("cme_jetDR", &b_cme_jetDR, "cme_jetDR/F");
  t->Branch("cme_legDR", &b_cme_legDR, "cme_legDR/F");
  t->Branch("cme_dca", &b_cme_dca, "cme_dca/F");
  t->Branch("cme_chi2", &b_cme_chi2, "cme_chi2/F");
  t->Branch("cme_charge", &b_cme_charge, "cme_charge/I");
  t->Branch("cme_cut", &b_cme_cut, "cme_cut/I");
  //t->Branch("cme_nJet", &b_cme_nJet, "cme_nJet/I");
  
  t->Branch("cme_dau1_chi2", &b_cme_dau1_chi2,"cme_dau1_chi2/F");
  t->Branch("cme_dau1_nHits", &b_cme_dau1_nHits,"cme_dau1_nHits/F");
  t->Branch("cme_dau1_pt", &b_cme_dau1_pt,"cme_dau1_pt/F");
  t->Branch("cme_dau2_chi2", &b_cme_dau2_chi2,"cme_dau2_chi2/F");
  t->Branch("cme_dau2_nHits", &b_cme_dau2_nHits,"cme_dau2_nHits/F");
  t->Branch("cme_dau2_pt", &b_cme_dau2_pt,"cme_dau2_pt/F");
  ////t->Branch("cme_dau1_ipsigXY", &b_cme_dau1_ipsigXY,"cme_dau1_ipsigXY/F");
  ////t->Branch("cme_dau1_ipsigZ", &b_cme_dau1_ipsigZ,"cme_dau1_ipsigZ/F");
  //t->Branch("cme_dau2_ipsigXY", &b_cme_dau2_ipsigXY,"cme_dau2_ipsigXY/F");
  //t->Branch("cme_dau2_ipsigZ", &b_cme_dau2_ipsigZ,"cme_dau2_ipsigZ/F");
  //t->Branch("cme_jet_btagDeepC", &b_cme_jet_btagDeepC, "cme_jet_btagDeepC/I");
  //t->Branch("cme_pt", &b_cme_pt,"cme_pt/F");
  //t->Branch("cme_eta", &b_cme_eta,"cme_eta/F");
  //t->Branch("cme_phi", &b_cme_phi,"cme_phi/F");
  //t->Branch("cme_angleXYZ", &b_cme_angleXYZ, "cme_angleXYZ/F");
  //t->Branch("cme_angleXY", &b_cme_angleXY,"cme_angleXY/F");
  //t->Branch("cme_x", &b_cme_x, "cme_x/F");
  ///t->Branch("cme_y", &b_cme_y, "cme_y/F");
  //t->Branch("cme_z", &b_cme_z, "cme_z/F");
  
}


void slmassAnalyser::resetBranch() {
    Reset();
    b_cme_cut = 0; 
    b_cme_charge = 0; 
    b_had.SetPtEtaPhiM(0,0,0,0); 

    b_cme_diffMass = -88;
    b_cme_pdgId = 0;
    b_cme_tmva = -2;
    //b_cme_tmva_Jpsi = -2;
    b_cme_mass = -88;
    
    b_cme_lxy = -88;
    b_cme_l3DE = -88;
    b_cme_l3DErr = -88;
    b_cme_jetDR = -88;
    b_cme_legDR = -88;
    b_cme_nJet = -88;
    b_cme_dca = -88;
    b_cme_angleXYZ = -88;
    b_cme_x = -88;
    b_cme_y = -88;
    b_cme_z = -88;
    b_cme_chi2 = -88;
    b_cme_jet_btagDeepB = -88;
    b_cme_dau2_chi2 = -88;
    b_cme_dau1_chi2 = -88;
    b_cme_dau1_ipsigXY = -88;
    b_cme_dau2_ipsigXY = -88;
    b_cme_dau1_nHits = -88;
    b_cme_dau2_nHits = -88;
    b_cme_dau2_ipsigZ = -88;
  
    b_cme_lxyE = -88;
    b_cme_l3D = -88;
    b_cme_eta = -88;
    b_cme_phi = -88;
    b_cme_jet_btagCMVA = -88;
    b_cme_jet_btagCSVV2 = -88;
    b_cme_angleXY = -88;
    b_cme_pt =-88;
    b_cme_jet_btagDeepC = -88;
    b_cme_dau1_ipsigZ =-88;
    b_cme_dau1_pt =-88;
    b_cme_dau2_pt =-88;
}

void slmassAnalyser::cmesonSelection() {
    UInt_t hadnum = 100001;
    Float_t max_tmva = -2.0f;
    if (nhad < 1) return;
    for (UInt_t k = 0; k < nhad; ++k) {
        if (had_dau1_chi2[k] >4) continue;
        if (had_dau1_nHits[k] <3) continue;
        if (had_dau1_pt[k] <0.5) continue;
        if (had_dau2_chi2[k] >3) continue;
        if (had_dau2_nHits[k] <3) continue;
        if (had_dau2_pt[k] <0.5) continue;
        if (std::fabs(had_angleXY[k]) <0.95) continue;
        if (std::fabs(had_x[k]) >8) continue;
        if (std::fabs(had_y[k]) >8) continue;
        if (std::fabs(had_z[k]) >20) continue;

        if (had_l3DErr[k] <= 0) continue;
        tmp_l3DE = had_l3D[k] / had_l3DErr[k];
        if (TMath::IsNaN(tmp_l3DE)or tmp_l3DE>200) continue;
        
        tmp_jetDR = had_jetDR[k];
        if (had_jetDR[k] >0.3) continue;
        
        tmp_legDR = had_legDR[k];
        if (had_legDR[k] >0.6) continue;
        
        tmp_dca = had_dca[k];
        if (had_dca[k] >1) continue;
        
        tmp_chi2 = had_chi2[k];
        if (had_chi2[k] >10) continue;
        
        tmp_jet_btagCSVV2 = had_jet_btagCSVV2[k];
        if (had_jet_btagCSVV2[k] <0.05) continue;
        
        tmp_mass = had_mass[k];
        tmp_tmva = bdtg->EvaluateMVA("BDTG");
        
        if (max_tmva < tmp_tmva){
            max_tmva = tmp_tmva;
            hadnum = k;    
        }
        
    }     
    //cout <<hadnum<<std::endl;
    if (hadnum == 100001) return;
    b_cme_cut = 1;
    b_cme_tmva = max_tmva;
    b_cme_diffMass = had_diffMass[hadnum];
    b_cme_mass = had_mass[hadnum];
    b_cme_pdgId = had_pdgId[hadnum];
   
    b_cme_dau1_chi2 = had_dau1_chi2[hadnum];
    b_cme_dau1_nHits = had_dau1_nHits[hadnum];
    b_cme_dau1_pt = had_dau1_pt[hadnum]; 
    b_cme_dau2_chi2 = had_dau2_chi2[hadnum];
    b_cme_dau2_nHits = had_dau2_nHits[hadnum];
    b_cme_dau2_pt = had_dau2_pt[hadnum];
    b_cme_angleXY = had_angleXY[hadnum];
    b_cme_x = had_x[hadnum];
    b_cme_y = had_y[hadnum];
    b_cme_z = had_z[hadnum];

    b_cme_l3DE = had_l3D[hadnum] / had_l3DErr[hadnum];
    b_cme_jetDR = had_jetDR[hadnum];
    b_cme_legDR = had_legDR[hadnum];
    b_cme_dca = had_dca[hadnum];
    b_cme_chi2 = had_chi2[hadnum];
    b_cme_jet_btagCSVV2 = had_jet_btagCSVV2[hadnum];

    TLorentzVector had_tlv;
    had_tlv.SetPtEtaPhiM(had_pt[hadnum], had_eta[hadnum], had_phi[hadnum], had_mass[hadnum]);
    b_had = had_tlv;

    TLorentzVector sum_tlv;
    sum_tlv = b_lep + b_had;
    b_sum = sum_tlv.M();
    if( b_lep_pid * b_cme_pdgId < 0 ){
   // cout <<"-----------------"<<std::endl;
   // cout <<b_lep_pid<<std::endl;
   // cout <<b_cme_pdgId<<std::endl;
   // cout <<"-----------------"<<std::endl;
    b_cme_charge = -1; 
    }
    else {b_cme_charge = 1;}

}


