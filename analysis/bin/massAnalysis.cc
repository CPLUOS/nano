#include "nano/analysis/interface/massAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <cstdlib>
using namespace std;

bool massAnalysis::analysis()
{
  h_cutFlow->Fill(0);

  //Run for MC
  if (m_isMC) {
    Int_t nvtx = Pileup_nTrueInt;
    b_puweight = m_pileUp->getWeight(nvtx);
      
    b_genweight = genWeight;
    h_genweights->Fill(0.5, b_genweight);
    b_weight = b_genweight * b_puweight;
  }
  else
  {
    b_puweight = 1;
    b_genweight = 0;
    if (!(m_lumi->LumiCheck(run, luminosityBlock))) return false;
  }
  h_nevents->Fill(0.5,b_genweight*b_puweight); 
    
  h_cutFlow->Fill(1);
  if (std::abs(PV_z) >= 24.) return false;
  if (PV_npvs == 0) return false;
  if (PV_ndof < 4) return false;

  h_cutFlow->Fill(2);

  auto muons = muonSelection();
  auto elecs = elecSelection();

  if(muons.size()+ elecs.size() != 2) return false;

  h_cutFlow->Fill(3);

  int mulpdg = -1;
  if (muons.size() == 2){
      recolep1 = muons[0];
      recolep2 = muons[1];
      mulpdg = muons[0].GetPdgCode()*muons[1].GetPdgCode();
      b_channel = CH_MUMU;
  }
  
  if (muons.size() == 1 && elecs.size() == 1){
      recolep1 = muons[0];
      recolep2 = elecs[0];
      mulpdg = muons[0].GetPdgCode()*elecs[0].GetPdgCode();
      b_channel = CH_MUEL;
  }
  if (elecs.size() == 2){
      recolep1 = elecs[0];
      recolep2 = elecs[1];
      mulpdg = elecs[0].GetPdgCode()*elecs[1].GetPdgCode();
      b_channel = CH_ELEL;
  }

  recolep1.Momentum(b_lep1);
  recolep2.Momentum(b_lep2);

  recoleps.push_back(b_lep1);
  recoleps.push_back(b_lep2);

  b_dilep = b_lep1 + b_lep2;

  //Triggers
  b_trig_m = HLT_IsoTkMu24 || HLT_IsoMu24;
  b_trig_e = HLT_Ele27_WPTight_Gsf;  
  b_trig_mm = HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL
    || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;
  b_trig_em = HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL
    || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
  b_trig_ee = HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;


  if (b_channel == CH_MUMU){
    if (m_isMC){
      if (!(b_trig_mm || b_trig_m)) return false;
    }
    if (m_isDL){
      if (!(b_trig_mm)) return false;
    }
    if (m_isSL_m){
      if (b_trig_mm||!b_trig_m) return false;  
    }
  }

  if (b_channel == CH_MUEL){
    if (m_isMC){
      if (!(b_trig_em || b_trig_m || b_trig_e)) return false;
    }
    if (m_isDL){
      if (!(b_trig_em)) return false;
    }
    if (m_isSL_e) {
      if (b_trig_em || !b_trig_e || b_trig_m) return false;
    }
    if (m_isSL_m) {
      if (b_trig_em || b_trig_e || !b_trig_m) return false;
    } 
  }

  if (b_channel == CH_ELEL){
    if (m_isMC){
      if (!(b_trig_ee || b_trig_e)) return false;
    }
    if (m_isDL){
      if (!b_trig_ee) return false;
    }
    if (m_isSL_e){
      if (b_trig_ee || !b_trig_e) return false;
    }
  }

  //leptonSF
  b_mueffweight    = muonSF_.getScaleFactor(recolep1, 13, 0)*muonSF_.getScaleFactor(recolep2, 13,  0);
  b_mueffweight_up = muonSF_.getScaleFactor(recolep1, 13, +1)*muonSF_.getScaleFactor(recolep2, 13, +1);
  b_mueffweight_dn = muonSF_.getScaleFactor(recolep1, 13, -1)*muonSF_.getScaleFactor(recolep2, 13, -1);

  b_eleffweight    = elecSF_.getScaleFactor(recolep1, 11, 0)*elecSF_.getScaleFactor(recolep2, 11,  0);
  b_eleffweight_up = elecSF_.getScaleFactor(recolep1, 11, +1)*elecSF_.getScaleFactor(recolep2, 11, +1);
  b_eleffweight_dn = elecSF_.getScaleFactor(recolep1, 11, -1)*elecSF_.getScaleFactor(recolep2, 11, -1);

  
  b_tri = b_tri_up = b_tri_dn = 0;
  b_tri = computeTrigSF(recolep1, recolep2);
  b_tri_up = computeTrigSF(recolep1, recolep2,  1);
  b_tri_dn = computeTrigSF(recolep1, recolep2, -1);

  auto jets = jetSelection();
  auto bjets = bjetSelection();
  
  if (b_dilep.M() < 20. || mulpdg > 0) return false;
  b_step1 = true;
  b_step = 1;
  h_cutFlow->Fill(4);

  if (b_channel == CH_MUEL || b_dilep.M() < 76 || b_dilep.M() > 106){
    b_step2 = true;
    b_step = 2;
    h_cutFlow->Fill(5);
  }

  b_met = MET_pt;
  b_njet = jets.size();
  b_nbjet = bjets.size();

  if (b_channel == CH_MUEL || b_met > 40){ 
    b_step3 = true;
    if (b_step == 2){
      ++b_step;
      h_cutFlow->Fill(6);
    }
  }
  if (b_njet >= 2){
    b_step4 = true;
    if (b_step == 3){
      ++b_step;
      h_cutFlow->Fill(7);
    }
  }
  if (b_nbjet > 0){
    b_step5 = true;
    if (b_step == 4){
      ++b_step;
      h_cutFlow->Fill(8);
    }
  }
  return true;
}

void massAnalysis::Loop()
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //Prepare for new loop
    resetBranch();
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry); nbytes += nb;
    bool keep = analysis();
    //cout << keep << endl;
    if (keep){
      collectTMVAvalues();
      m_tree->Fill();
    }
  }
}

int main(int argc, char* argv[])
{
  string env = getenv("CMSSW_BASE");
  string username = getenv("USER");

  if(argc != 1)
  {
    std::string dirName = "root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/"+username+"/nanoAOD/"+std::string(argv[1])+"/"+std::string(argv[2]);
    std::string temp = argv[2];
    
    Bool_t isDL = false;
    Size_t found_DL = temp.find("Double");
    if(found_DL != std::string::npos) isDL = true;

    Bool_t isSL_e = false;
    Size_t found_SL_e = temp.find("SingleElectron");
    if(found_SL_e != std::string::npos) isSL_e = true;
    
    Bool_t isSL_m = false;
    Size_t found_SL_m = temp.find("SingleMuon");
    if(found_SL_m != std::string::npos) isSL_m = true;

    Bool_t isMC = false;
    Size_t found = temp.find("Run");
    if(found == std::string::npos) isMC = true;

    for(Int_t i = 3; i < argc; i++)
    {
      TFile *f = TFile::Open(argv[i], "read");
     
      TTree *tree;                  
      f->GetObject("Events", tree);
      
      temp = argv[i];   
      found = temp.find_last_of('/');
      std::string outPutName = dirName+temp.substr(found);
      massAnalysis t(tree, isMC, isDL, isSL_e, isSL_m);
      
      t.setOutput(outPutName);
      t.Loop();
    }
  }
  else
  {
    TFile *f = TFile::Open("/xrootd/store/group/nanoAOD/run2_2016v4/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180430_152541/0000/nanoAOD_256.root", "read");
    //TFile *f = TFile::Open("/xrootd/store/group/nanoAOD/run2_2016v2/DoubleMuon/Run2016B-18Apr2017_ver1-v1/171219_063744/0000/nanoAOD_23.root", "read");

    TTree *tree;
    f->GetObject("Events", tree);

    massAnalysis t(tree, true);
    t.setOutput("test.root");
    t.Loop();
  }
  return 0;
}

void massAnalysis::collectTMVAvalues()
{
  for(UInt_t i=0; i < nhad; ++i){
    //if(had_pdgId[i] != 421) continue;
    b_cme_lxy = had_lxy[i];
    b_cme_lxyE = had_lxy[i] / had_lxyErr[i];
    b_cme_l3D = had_l3D[i];
    b_cme_l3DE = had_l3D[i] / had_l3DErr[i];
    b_cme_jetDR = had_jetDR[i];
    b_cme_legDR = had_legDR[i];
    b_cme_dca = had_dca[i];
    b_cme_angleXY = had_angleXY[i];
    b_cme_angleXYZ = had_angleXYZ[i];
    b_cme_x = had_x[i];
    b_cme_y = had_y[i];
    b_cme_z = had_z[i];
    b_cme_pt = had_pt[i];
    b_cme_chi2 = had_chi2[i];
    b_cme_eta = had_eta[i];
    b_cme_phi = had_phi[i];
    b_cme_jet_btagCMVA = had_jet_btagCMVA[i];
    b_cme_jet_btagCSVV2 = had_jet_btagCSVV2[i];
    b_cme_jet_btagDeepB = had_jet_btagDeepB[i];
    b_cme_jet_btagDeepC = had_jet_btagDeepC[i];
    b_cme_dau1_chi2 = had_dau1_chi2[i];
    b_cme_dau1_ipsigXY = had_dau1_ipsigXY[i];
    b_cme_dau1_ipsigZ = had_dau1_ipsigZ[i];
    b_cme_dau1_nHits = had_dau1_nHits[i];
    b_cme_dau1_pt = had_dau1_pt[i];
    b_cme_dau2_chi2 = had_dau2_chi2[i];
    b_cme_dau2_ipsigXY = had_dau2_ipsigXY[i];
    b_cme_dau2_ipsigZ = had_dau2_ipsigZ[i];
    b_cme_dau2_nHits = had_dau2_nHits[i];
    b_cme_dau2_pt = had_dau2_pt[i];
    b_cme_mass = had_mass[i];
    b_cme_tmva_bdtg = bdtg->EvaluateMVA("BDTG");
    //b_cme_nMatched = hadTruth_nMatched[i];
    if (b_cme_tmva_bdtg > b_bdtg)
    {
        b_maxbIdx = i;
        b_bdtg = b_cme_tmva_bdtg;
    }
  }
      ////printf("i %d",b_maxbIdx);
  if (b_maxbIdx != -1 && b_bdtg != -1)
  {
    b_cme_mass = had_mass[b_maxbIdx];
    b_cme_tmva_bdtg = b_bdtg;
    b_cme_pdgId = had_pdgId[b_maxbIdx];
    m_tree->Fill();
  }

}


void massAnalysis::setOutput(std::string outputName)
{
  m_output = TFile::Open(outputName.c_str(), "recreate");

  m_tree = new TTree("event", "event");
  MakeBranch(m_tree);

  bdtg = new TMVA::Reader();
  bdtg->AddVariable( "cme_lxy", &b_cme_lxy );
  bdtg->AddVariable( "cme_lxyE", &b_cme_lxyE );
  bdtg->AddVariable( "cme_l3D", &b_cme_l3D );
  bdtg->AddVariable( "cme_l3DE", &b_cme_l3DE );
  bdtg->AddVariable( "cme_jetDR", &b_cme_jetDR );
  bdtg->AddVariable( "cme_legDR", &b_cme_legDR );
  bdtg->AddVariable( "cme_dca", &b_cme_dca );
  bdtg->AddVariable( "cme_angleXY", &b_cme_angleXY );
  bdtg->AddVariable( "cme_angleXYZ", &b_cme_angleXYZ );
  bdtg->AddVariable( "cme_x", &b_cme_x );
  bdtg->AddVariable( "cme_y", &b_cme_y );
  bdtg->AddVariable( "cme_z", &b_cme_z );
  bdtg->AddVariable( "cme_pt", &b_cme_pt );
  bdtg->AddVariable( "cme_chi2", &b_cme_chi2 );
  bdtg->AddVariable( "cme_eta", &b_cme_eta );
  bdtg->AddVariable( "cme_phi", &b_cme_phi );
  bdtg->AddVariable( "cme_jet_btagCMVA", &b_cme_jet_btagCMVA );
  bdtg->AddVariable( "cme_jet_btagCSVV2", &b_cme_jet_btagCSVV2 );
  bdtg->AddVariable( "cme_jet_btagDeepB", &b_cme_jet_btagDeepB );
  bdtg->AddVariable( "cme_jet_btagDeepC", &b_cme_jet_btagDeepC );
  bdtg->AddVariable( "cme_dau1_chi2", &b_cme_dau1_chi2 );
  bdtg->AddVariable( "cme_dau1_ipsigXY", &b_cme_dau1_ipsigXY );
  bdtg->AddVariable( "cme_dau1_ipsigZ", &b_cme_dau1_ipsigZ );
  bdtg->AddVariable( "cme_dau1_nHits", &b_cme_dau1_nHits );
  bdtg->AddVariable( "cme_dau1_pt", &b_cme_dau1_pt );
  bdtg->AddVariable( "cme_dau2_chi2", &b_cme_dau2_chi2 );
  bdtg->AddVariable( "cme_dau2_ipsigXY", &b_cme_dau2_ipsigXY );
  bdtg->AddVariable( "cme_dau2_ipsigZ", &b_cme_dau2_ipsigZ );
  bdtg->AddVariable( "cme_dau2_nHits", &b_cme_dau2_nHits );
  bdtg->AddVariable( "cme_dau2_pt", &b_cme_dau2_pt );
  bdtg->AddSpectator( "cme_mass", &b_cme_mass );
  bdtg->BookMVA("BDTG", "/cms/scratch/seulgi/nanoAOD/src/nano/analysis/test/topMass/cut/tmva/xml/TMVAClassification_BDTG.weights.xml");
  
  h_nevents = new TH1D("nevents", "nevents", 1, 0, 1);
  h_genweights = new TH1D("genweight", "genweight", 1, 0, 1);
  h_weights = new TH1D("weight", "weight", 1, 0, 1);
  h_cutFlow = new TH1D("cutflow", "cutflow", 11, -0.5, 10.5);
}

void massAnalysis::MakeBranch(TTree* t)
{
  t->Branch("nvertex", &b_nvertex, "nvertex/I");
  t->Branch("step", &b_step, "step/I");
  t->Branch("channel", &b_channel, "channel/I");
  t->Branch("njet", &b_njet, "njet/I");
  t->Branch("nbjet", &b_nbjet, "nbjet/I");
  t->Branch("step1", &b_step1, "step1/O");
  t->Branch("step2", &b_step2, "step2/O");
  t->Branch("step3", &b_step3, "step3/O");
  t->Branch("step4", &b_step4, "step4/O");
  t->Branch("step5", &b_step5, "step5/O");
  t->Branch("step6", &b_step6, "step6/O");
  t->Branch("step7", &b_step7, "step7/O");
  
  m_tree->Branch("lep1", "TLorentzVector", &b_lep1);
  m_tree->Branch("lep1_pid", &b_lep1_pid, "lep1_pid/I");    
  m_tree->Branch("lep2", "TLorentzVector", &b_lep2);
  m_tree->Branch("lep2_pid", &b_lep2_pid, "lep2_pid/I");    
  t->Branch("dilep", "TLorentzVector", &b_dilep);
  //m_tree->Branch("jet1", "TLorentzVector", &b_jet1);
  //m_tree->Branch("jet1_CSVInclV2", &b_jet1_CSVInclV2, "jet1_CSVInclV2/F");
  //m_tree->Branch("jet2", "TLorentzVector", &b_jet2);
  //m_tree->Branch("jet2_CSVInclV2", &b_jet2_CSVInclV2, "jet2_CSVInclV2/F");


  t->Branch("tri", &b_tri, "tri/F");
  t->Branch("tri_up", &b_tri_up, "tri_up/F");
  t->Branch("tri_dn", &b_tri_dn, "tri_dn/F");

  t->Branch("met", &b_met, "met/F");
  t->Branch("weight", &b_weight, "weight/F");
  t->Branch("puweight", &b_puweight, "puweight/F");
  t->Branch("genweight", &b_genweight, "genweight/F");
  t->Branch("csvweight", "std::vector<float>", &b_csvweights);
  t->Branch("btagweight", &b_btagweight, "btagweight/F");
  t->Branch("mueffweight", &b_mueffweight, "mueffweight/F");
  t->Branch("eleffweight", &b_eleffweight, "eleffweight/F");
  t->Branch("PV_npvs", &PV_npvs, "PV_npvs/I");

  t->Branch("trig_m", &b_trig_m, "trig_m/O");
  t->Branch("trig_e", &b_trig_e, "trig_e/O");
  t->Branch("trig_mm", &b_trig_mm, "trig_mm/O");
  t->Branch("trig_em", &b_trig_em, "trig_em/O");
  t->Branch("trig_ee", &b_trig_ee, "trig_ee/O");
  
  t->Branch("cme_tmva_bdtg", &b_cme_tmva_bdtg, "cme_tmva_bdtg/F");
  t->Branch("cme_mass", &b_cme_mass, "cme_mass/F");
  t->Branch("cme_pdgId", &b_cme_pdgId, "cme_pdgId/I");
  //t->Branch("nhad", &nhad, "nhad/I");
  //t->Branch("had_dca", &had_dca, "had_dca/F");
  //t->Branch("had_angleXY", &had_angleXY, "had_angleXY/F");
  //t->Branch("had_angleXYZ", &had_angleXYZ, "had_angleXYZ/F");
  //t->Branch("had_dau1_chi2", &had_dau1_chi2, "had_dau1_chi2/F");
  //t->Branch("had_dau1_nHits", &had_dau1_nHits, "had_dau1_nHits/F");
  //t->Branch("had_dau1_pt", &had_dau1_pt, "had_dau1_pt/F");
  //t->Branch("had_dau1_ipsigXY", &had_dau1_ipsigXY, "had_dau1_ipsigXY/F");
  //t->Branch("had_dau1_ipsigZ", &had_dau1_ipsigZ, "had_dau1_ipsigZ/F");
  //t->Branch("had_lxy", &had_lxy, "had_lxy/F");
  //t->Branch("had_lxyErr", &had_lxyErr, "had_lxyErr/F");
  //t->Branch("had_l3D", &had_l3D, "had_l3D/F");
  //t->Branch("had_l3DErr", &had_l3DErr, "had_l3DErr/F");
  //t->Branch("had_jetDR", &had_jetDR, "had_jetDR/F");
  //t->Branch("had_legDR", &had_legDR, "had_legDR/F");
  //t->Branch("had_diffMass", &had_diffMass, "had_diffMass/F");
  //t->Branch("had_nJet", &had_nJet, "had_nJet/I");
  //t->Branch("had_chi2", &had_chi2, "had_chi2/F");
  //t->Branch("had_eta", &had_eta, "had_eta/F");
  //t->Branch("had_mass", &had_mass, "had_mass/F");
  //t->Branch("had_phi", &had_phi, "had_phi/F");
  //t->Branch("had_pt", &had_pt, "had_pt/F");
  //t->Branch("had_x", &had_x, "had_x/F");
  //t->Branch("had_y", &had_y, "had_y/F");
  //t->Branch("had_z", &had_z, "had_z/F");
  //t->Branch("had_ndof", &had_ndof, "had_ndof/I");
  //t->Branch("had_pdgId", &had_pdgId, "had_pdgId/I");
}


void massAnalysis::resetBranch()
{
  b_lep1.SetPtEtaPhiM(0,0,0,0);
  b_lep2.SetPtEtaPhiM(0,0,0,0);
  b_dilep.SetPtEtaPhiM(0,0,0,0);

  //b_jet1.SetPtEtaPhiM(0,0,0,0);
  //b_jet2.SetPtEtaPhiM(0,0,0,0);
  recoleps.clear();
  b_csvweights.clear();
  b_lep1_pid = 0; b_lep2_pid = 0;
  b_jet1_CSVInclV2 = -1; b_jet2_CSVInclV2 = -1;

  b_nvertex = 0; b_step = -1; b_channel = 0; b_njet = 0; b_nbjet = 0;
  b_step1 = 0; b_step2 = 0; b_step3 = 0; b_step4 = 0; b_step5 = 0; b_step6 = 0; b_step7 = 0;
  b_met = -9; b_weight = 1; b_genweight = 1; b_puweight = 1; b_btagweight = 1;
  b_tri = 0;
  b_mueffweight = 1;b_mueffweight_up = 1;b_mueffweight_dn = 1;
  b_eleffweight = 1;b_eleffweight_up = 1;b_eleffweight_dn = 1;

  b_cme_mass = -999;
  b_cme_tmva_bdtg = -999;
  b_cme_pdgId = 0;
  b_bdtg = -1; b_maxbIdx = -1;

}
