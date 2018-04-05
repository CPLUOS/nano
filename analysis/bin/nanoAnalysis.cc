#define nanoAnalysis_cxx
#include "nano/analysis/src/nanoAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <cstdlib>

using namespace std;
/*
To compile:
g++ `root-config --cflags --glibs` -lEG nanoAnalysis.cc -o nanoAnalysis
*/

void nanoAnalysis::SetOutput(string outputName)
{
  m_output = TFile::Open(outputName.c_str(), "recreate");

  m_tree = new TTree("events", "events");
  MakeBranch(m_tree);
  
  h_Event_Tot = new TH1D("Event_total", "Event_total" ,1,0,1);
  h_genweights = new TH1D("genweight", "genweight" , 1,0,1);
  h_weight = new TH1D("weight", "weight", 1,0,1);
  h_cutFlow = new TH1D("cutflow", "cutflow", 11, -0.5, 10.5);
  h_FL = new TH1D("FullyLeptonic", "FullyLeptonic", 1, 0, 1);
  h_FH2 = new TH1D("FullyLeptonic2", "FullyLeptonic2", 1, 0, 1);
  h_FH3 = new TH1D("FullyLeptonic3", "FullyLeptonic3", 1, 0, 1);
  h_FH4 = new TH1D("FullyLeptonic4", "FullyLeptonic4", 1, 0, 1);
  h_SL = new TH1D("SemiLeptonic", "SemiLeptonic", 1, 0, 1);
  h_Out = new TH1D("Leftover", "Leftover", 1, 0, 1);
}

void nanoAnalysis::MakeBranch(TTree* t)
{
  t->Branch("Event_No", &b_Event_No, "Event_No/I");
  t->Branch("Step", &b_Step, "Step/I");
  t->Branch("Dilep", "TLorentzVector", &b_Dilep);
  t->Branch("Mu1", "TLorentzVector", &b_Mu1);
  t->Branch("Mu2", "TLorentzVector", &b_Mu2);
  t->Branch("lep1", "TLorentzVector", &b_lep1);
  t->Branch("lep2", "TLorentzVector", &b_lep2);
  t->Branch("genweight", &b_genweight, "genweight/F");
  t->Branch("puweight", &b_puweight, "puweight/F");
  t->Branch("puweight_up", &b_puweight_up, "puweight_up/F");
  t->Branch("puweight_dn", &b_puweight_dn, "puweight_dn/F");
  t->Branch("weight", &b_weight, "weight/F");
  t->Branch("npvs", &b_npvs, "npvs/F");
  t->Branch("channel", &b_channel, "channel/I");
  t->Branch("charge", &b_charge, "charge/I");
  t->Branch("nlep", &b_nlep, "nlep/F");
  t->Branch("nmuon", &b_nmuon, "nmuon/F");
  t->Branch("nelec", &b_nelec, "nelec/F");
  t->Branch("njet", &b_njet, "njet/F");
  t->Branch("nbjet", &b_nbjet, "nbjet/F");
  t->Branch("trig_m", &b_trig_m, "trig_m/O");
  t->Branch("trig_m2", &b_trig_m2, "trig_m2/O");
  t->Branch("Met", &b_Met, "Met/F");
  t->Branch("Met_phi", &b_Met_phi, "Met_phi/F");
  t->Branch("CSVv2", &b_CSVv2);
  t->Branch("FL", &b_FL, "FL/I");
  t->Branch("FH2", &b_FH2, "FH2/I");
  t->Branch("FH3", &b_FH3, "FH3/I");
  t->Branch("FH4", &b_FH4, "FH4/I");
  t->Branch("SL", &b_SL, "SL/I");
  t->Branch("csvweight", "std::vector<float>", &b_csvweights);
  t->Branch("csvweight", "std::vector<float>", &b_csvweights);
  t->Branch("btagweight", &b_btagweight, "btagweight/F");
  t->Branch("btagweight_up", &b_btagweight_up, "btagweight_up/F");
  t->Branch("btagweight_dn", &b_btagweight_dn, "btagweight_dn/F");
  t->Branch("lf__up", &b_lf_up, "lf_up/F");
  t->Branch("lf__dn", &b_lf_dn, "lf_dn/F");
  t->Branch("hfstats1_up", &b_hfstats1_up, "hfstats1_up/F");
  t->Branch("hfstats1_dn", &b_hfstats1_dn, "hfstats1_dn/F");
  t->Branch("hfstats2_up", &b_hfstats2_up, "hfstats2_up/F");
  t->Branch("hfstats2_dn", &b_hfstats2_dn, "hfstats2_dn/F");
  t->Branch("mueffweight", &b_mueffweight, "mueffweight/F");
  t->Branch("mueffweight_up", &b_mueffweight_up, "mueffweight_up/F");
  t->Branch("mueffweight_dn", &b_mueffweight_dn, "mueffweight_dn/F");
  /////////////////////BDT////////////////////////////
  t->Branch("all_muEtaDiff", &b_all_muEtaDiff, "all_muEtaDiff/F");
  t->Branch("all_muPtDiff", &b_all_muPtDiff, "all_muPtDiff/F");
  t->Branch("all_muPhiDiff", &b_all_muPhiDiff, "all_muPhiDiff/F");
  t->Branch("all_muDR", &b_all_muDR, "all_muDR/F");
  t->Branch("all_Dilep_Pt", &b_all_Dilep_Pt, "all_Dilep_Pt/F");
  t->Branch("all_Dilep_Eta", &b_all_Dilep_Eta, "all_Dilep_Eta/F");
  //FL//
  t->Branch("FL_lepEtaDiff", &b_FL_lepEtaDiff, "FL_lepEtaDiff/F");
  t->Branch("FL_lepPtDiff", &b_FL_lepPtDiff, "FL_lepPtDiff/F");
  t->Branch("FL_lepPhiDiff", &b_FL_lepPhiDiff, "FL_lepPhiDiff/F");
  t->Branch("FL_lep1Pt", &b_FL_lep1Pt, "FL_lep1Pt/F");
  t->Branch("FL_lep1Eta", &b_FL_lep1Eta, "FL_lep1Eta/F");
  t->Branch("FL_lep2Pt", &b_FL_lep2Pt, "FL_lep2Pt/F");
  t->Branch("FL_lep2Eta", &b_FL_lep2Eta, "FL_lep2Eta/F");
  t->Branch("FL_lepDR", &b_FL_lepDR, "FL_lepDR/F");
  t->Branch("FL_nJet", &b_FL_nJet, "FL_nJet/F");
  t->Branch("FL_nbJet", &b_FL_nbJet, "FL_nbJet/F");
  //FH//
  t->Branch("Central_Jets", &b_Central_Jets, "Central_Jets/F");
  t->Branch("Forward_Jets", &b_Forward_Jets, "Forward_Jets/F");
}

void nanoAnalysis::ResetBranch()
{
  b_Event_No = 0;
  b_Step = 0;
  b_Dilep.SetPtEtaPhiM(0,0,0,0);
  b_Mu1.SetPtEtaPhiM(0,0,0,0);
  b_Mu2.SetPtEtaPhiM(0,0,0,0);
  b_lep1.SetPtEtaPhiM(0,0,0,0);
  b_lep2.SetPtEtaPhiM(0,0,0,0);
  b_Mu_tlv.clear();
  b_El_tlv.clear();
  b_Jet_tlv.clear();
  b_bJet_tlv.clear();
  b_csvweights.clear();
  b_CSVv2.clear();
  b_Event_Total = 1;
  b_channel = -1;
  b_nlep = 0; b_nmuon = 0; b_nelec = 0;
  b_charge = 0;
  b_Met_phi = 0; b_Met = 0;
  b_FL = 0; b_SL = 0; b_FH2 = 0; b_FH3 = 0; b_FH4 = 0;
}

void nanoAnalysis::LoadModules(pileUpTool* pileUp, lumiTool* lumi, RoccoR* rocCor)
{
  //Get Modules
  m_rocCor = rocCor;
  m_lumi = lumi;
  m_pileUp = pileUp;
  m_btagSF.initCSVWeight();
}

bool nanoAnalysis::Analysis()
{
  h_Event_Tot->Fill(0.5, b_Event_Total);
  h_cutFlow->Fill(0);
  //Run for MC
  if(m_isMC){
    Int_t nvtx = Pileup_nTrueInt;
    b_puweight = m_pileUp->getWeight(nvtx);
    b_puweight_up = m_pileUp->getWeight(nvtx, 1);
    b_puweight_dn = m_pileUp->getWeight(nvtx, -1);
    b_genweight = genWeight;
    h_genweights->Fill(0.5, b_genweight);
    b_weight = b_genweight * b_puweight;
    h_weight->Fill(0.5, b_weight);
  }
  else {
    b_puweight = 1;
    b_genweight = 0;
    if(!(m_lumi->LumiCheck(run, luminosityBlock))) return false;
  }
  h_cutFlow->Fill(1);

  if (abs(PV_z) >= 24.) return false;
  if (PV_npvs == 0) return false;
  if (PV_ndof < 4) return false;
  h_cutFlow->Fill(2);
  auto Muons = MuonSelection();
  auto Elecs = ElectronSelection(Muons);
  auto Jets = JetSelection(Muons, Elecs);
  auto BJets = BJetSelection(Muons);
 
  if (Muons.size() < 2) return false;

  b_trig_m = HLT_IsoTkMu24 || HLT_IsoMu24;
  
  // make all the variables that you need to save here
  TParticle mu1;
  TParticle mu2; 
  b_Met = PuppiMET_pt;
  b_Met_phi = PuppiMET_phi;
  Bool_t IsoMu24 = false;
  Bool_t IsoTkMu24 = false;

  for (UInt_t i = 0; i < nTrigObj; ++i){
    if (TrigObj_id[i] != 13) continue;
    if (TrigObj_pt[i] < 24) continue;
    Int_t bits = TrigObj_filterBits[i];
    if (bits & 0x2) IsoMu24 = true;
    if (bits & 0x8) IsoTkMu24 = true;  
  }
  if (!(IsoMu24 || IsoTkMu24)) return false; 
  
  for(UInt_t i = 0; i < Muons.size(); i++)
  { 
    if( (b_Mu_tlv[0].Pt() > 26) || (b_Mu_tlv[0].Pt() > 26) )
    { 
      if( ( Muons[0].GetPdgCode() * Muons[i].GetPdgCode() ) < 0 )
      { 
        b_Mu1 = b_Mu_tlv[0];
        b_Mu2 = b_Mu_tlv[i];
        b_charge = 1;
        mu1 = Muons[0];
        mu2 = Muons[i];
        break;
      }
    }
  }

  if (b_charge == 0) return false; 

  b_Dilep = b_Mu1 + b_Mu2;
  b_nlep = Muons.size() + Elecs.size();
  b_nmuon = Muons.size();
  b_nelec = Elecs.size();
  b_njet = Jets.size();
  b_nbjet = BJets.size();
  b_npvs = PV_npvs;
  
  /////////////////////// Fully Leptonic //////////////////
  if (Muons.size() + Elecs.size() == 4)
  {
    Int_t mulpdg = 1;

    if (Muons.size() == 2)
    {
      Elecs[0].Momentum(b_lep1);
      Elecs[1].Momentum(b_lep2);
      mulpdg = Elecs[0].GetPdgCode()*Elecs[1].GetPdgCode();
      b_channel = CH_ELEL;
    }

    else if (Muons.size() == 3)
    {
      for (UInt_t i = 1; i < Muons.size(); i++)
      {
        if (Muons[i].Pt() != b_Mu2.Pt())
        {
          mulpdg = Muons[i].GetPdgCode()*Elecs[0].GetPdgCode();
          if (Muons[i].Pt() > Elecs[0].Pt())
          {
            Muons[i].Momentum(b_lep1);
            Elecs[0].Momentum(b_lep2);
          }
          else 
          {
            Muons[i].Momentum(b_lep2);
            Elecs[0].Momentum(b_lep1);
          }
          break;
        }
      }
      b_channel = CH_MUEL;
    }

    else if(Muons.size() == 4)
    {
      for (UInt_t i = 1; i < Muons.size(); i++)
      {
        if (Muons[i].Pt() != b_Mu2.Pt())
        {
          for (UInt_t j = 1; i < Muons.size(); j++)
          {
            if (Muons[j].Pt() != b_Mu2.Pt() && Muons[i].Pt() != Muons[j].Pt())
            {
              mulpdg = Muons[i].GetPdgCode()*Muons[j].GetPdgCode();
              if (Muons[i].Pt() > Muons[j].Pt())
              {
                Muons[i].Momentum(b_lep1);
                Muons[j].Momentum(b_lep2);
              }
              else 
              {
                Muons[i].Momentum(b_lep2);
                Muons[j].Momentum(b_lep1);
              }
              break;
            }
          }
        }
      }
      b_channel = CH_MUMU;
    }

    if (mulpdg < 0)
    {
      b_FL = 1;
      h_FL->Fill(0.5, b_FL);
    }
  }
  ////////////////////// Fully Hadronic //////////////////
  if (Elecs.size() == 0 && Muons.size() == 2 && Jets.size() >= 4)
  {
     if ((BJets.size() == 1)||(BJets.size() == 2))
     {
       b_FH4 = 1;
       h_FH4->Fill(0.5, b_FH4);
     }
  }
  if (Elecs.size() == 0 && Muons.size() == 2 && Jets.size() >= 3)
  {
     if ((BJets.size() == 1)||(BJets.size() == 2))
     {
       b_FH3 = 1;
       h_FH3->Fill(0.5, b_FH3);
     }
  }
  
  if (Elecs.size() == 0 && Muons.size() == 2 && Jets.size() >= 2)
  {
     if ((BJets.size() == 1)||(BJets.size() == 2))
     {
       b_FH2 = 1;
       h_FH2->Fill(0.5, b_FH2);
     }
  }

  /////////////////////// Semi Leptonic //////////////////
  if (Muons.size() + Elecs.size() == 3)
  {
    if (BJets.size() == 1)
    {
      if (Jets.size() >= 2)
      {
        b_SL = 1;
        h_SL->Fill(0.5, b_SL);
      }
    }
    else if (BJets.size() == 2)
    {
      if (Jets.size() >=1)
      {
        b_SL = 1;
        h_SL->Fill(0.5, b_SL);
      }
    }
  }
  
  
  ////////// BDT //////////////////////////////
  /// BDT Varialbe Set ALL ///
  b_all_muEtaDiff = abs(b_Mu1.Eta() - b_Mu2.Eta());
  b_all_muPtDiff = abs(b_Mu1.Pt() - b_Mu2.Pt());
  b_all_muPhiDiff = abs(b_Mu1.Phi() - b_Mu2.Phi());
  b_all_muDR = abs(b_Mu1.DeltaR(b_Mu2));
  b_all_Dilep_Pt = b_Dilep.Pt();   
  b_all_Dilep_Eta = b_Dilep.Eta();   

  //FL
  if (b_FL == 1)
  {
    /// lepton ///
    b_FL_lepEtaDiff = abs(b_lep1.Eta() - b_lep2.Eta());
    b_FL_lepPtDiff = abs(b_lep1.Pt() - b_lep2.Pt());
    b_FL_lepPhiDiff = abs(b_lep1.Phi() - b_lep2.Phi());
    b_FL_lepDR = abs(b_lep1.DeltaR(b_lep2));
    b_FL_lep1Pt = b_lep1.Pt();
    b_FL_lep1Eta = b_lep1.Eta();
    b_FL_lep2Pt = b_lep2.Pt();
    b_FL_lep2Eta = b_lep2.Eta();
  }
  //SL 
  //FH2
  //FH3
  //FH4
  if (b_FH4 == 1 || b_FH3 == 1 || b_FH2 ==1) 
  {
    for (UInt_t i; i > Jet.size(); i++)
    {
      if (abs(Jet[i].Eta()) < 2.4)
      {
         b_Central_Jets ++ 1;   
      }
      else 
      {
         b_Forward_Jets ++ 1; 
      }
    }
  }
  
  b_mueffweight = m_muonSF.getScaleFactor(mu1, 13, 0)*m_muonSF.getScaleFactor(mu2, 13, 0);
  b_mueffweight_up = m_muonSF.getScaleFactor(mu1, 13, +1)*m_muonSF.getScaleFactor(mu2, 13, +1);
  b_mueffweight_dn = m_muonSF.getScaleFactor(mu1, 13, -1)*m_muonSF.getScaleFactor(mu2, 13, -1);

  b_Event_No = 1;
  if (b_SL == 0 && b_FH2 == 0 && b_FL == 0) h_Out->Fill(0.5, 1);

  return true;
}


void nanoAnalysis::Loop()
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //Prepare for new loop
    ResetBranch();
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    bool keep = Analysis();
    if (keep){
      m_tree->Fill();
    }
  }
  
}

int main(Int_t argc, Char_t** argv)
{
  string env = getenv("CMSSW_BASE");
  string username = getenv("USER");
  RoccoR* rocCor = new RoccoR(env+"/src/nano/analysis/data/rcdata.2016.v3/");
  lumiTool* lumi = new lumiTool(env+"/src/nano/analysis/data/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt");
  pileUpTool* pileUp = new pileUpTool();

  if(argc != 1)
  {
    string dirName = "root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/"+username+"/nanoAOD/"+std::string(argv[1])+"/"+std::string(argv[2]);
    string temp = argv[2];
    Bool_t isMC = false;
    Size_t found = temp.find("Run");
    if(found == string::npos) isMC = true;
    for(Int_t i = 3; i < argc; i++)
    {
      TFile *f = TFile::Open(argv[i], "read");

      TTree *tree;
      f->GetObject("Events", tree);

      temp = argv[i];
      found = temp.find_last_of('/');
      string outPutName = dirName+temp.substr(found);
      nanoAnalysis t(tree, isMC);
      t.LoadModules(pileUp, lumi, rocCor);
      t.SetOutput(outPutName);
      t.Loop();
    }
  }
  else
  {
    //TFile *f = TFile::Open("root://cms-xrdr.sdfarm.kr:1094///xrd/store/group/nanoAOD/run2_2016v3/ttHToMuMu_M125_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6/180125_131219/0000/nanoAOD_004.root", "read");
    TFile *f = TFile::Open("root://cms-xrdr.sdfarm.kr:1094///xrd/store/group/nanoAOD/run2_2016v3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/180125_131129/0000/nanoAOD_100.root", "read");
    TTree *tree;
    f->GetObject("Events", tree);
    
    nanoAnalysis t(tree, true);
    t.LoadModules(pileUp, lumi, rocCor);
    t.SetOutput("test.root");
    t.Loop();
  }

  return 0;
}

//Object Selections
vector<TParticle> nanoAnalysis::MuonSelection()
{
  vector<TParticle> muons;
  for(UInt_t i = 0; i < nMuon; i++)
  {
    if (!Muon_trackerMu[i]) continue;
    if (!Muon_globalMu[i]) continue;
    if (!Muon_tightId[i]) continue;
    if (Muon_pfRelIso04_all[i] > 0.25) continue;
    
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]);
    mom = mom * roccoR(mom, Muon_charge[i], Muon_genPartIdx[i], Muon_nTrackerLayers[i]);
    
    if (mom.Pt() < 10) continue;
    if (abs(mom.Eta() > 2.4)) continue;
   
    auto muon = TParticle();
    muon.SetPdgCode(13*Muon_charge[i]*-1);
    muon.SetMomentum(mom);

    b_Mu_tlv.push_back(mom);
    muons.push_back(muon);
  }
  return muons;
}

vector<TParticle> nanoAnalysis::ElectronSelection(vector<TParticle> leptons)
{
  vector<TParticle> electrons;
  for(UInt_t i = 0; i < nElectron; i++)
  {
    if ( Electron_pt[i] < 10) continue;
    if (abs(Electron_eta[i]) > 2.5 ) continue; //<~~~~~~~~~~~~~~ Higgs Electron pt == 10; Higgs Electron eta > 2.5  
    //if( Electron_pfRelIso03_all[i] > 0.15 || Electron_cutBased[i] < 3 ) continue; //<~~~~~~~~~~~~~~ Top doesn't use isolation? 
    if (Electron_cutBased[i] < 3) continue;
    float el_scEta = Electron_deltaEtaSC[i] + Electron_eta[i];
    if ( std::abs(el_scEta) > 1.4442 &&  std::abs(el_scEta) < 1.566 ) continue;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Electron_pt[i], Electron_eta[i], Electron_phi[i], Electron_mass[i]);

    if (hasOverLap(mom, leptons, 0.3 )) continue;
    
    auto elec = TParticle();
    elec.SetPdgCode(11*Electron_charge[i]*-1);
    elec.SetMomentum(mom);

    b_El_tlv.push_back(mom);
    electrons.push_back(elec);
  }
  return electrons;
}
vector<TParticle> nanoAnalysis::JetSelection(vector<TParticle> Muons, vector<TParticle> Elecs)
{
  vector<TParticle> jets;
  float Jet_SF_CSV[19] = {1.0,};
  for(UInt_t i = 0; i < nJet; i++)
  {
    if (Jet_pt[i] < 30) continue;
    if (abs(Jet_eta[i]) > 4.7) continue; 
    if (Jet_jetId[i] < 1) continue;

    TLorentzVector mom;
    mom.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
    
    if (hasOverLap(mom, Muons, 0.4)) continue;
    if (hasOverLap(mom, Elecs, 0.4)) continue;
    
    auto jet = TParticle();
    jet.SetMomentum(mom);

    jets.push_back(jet);
    for (UInt_t iu = 0; iu < 19; iu++)
    {
      Jet_SF_CSV[iu] *= m_btagSF.getSF(jet, Jet_btagCSVV2[i], Jet_hadronFlavour[i], iu);
    }
  }
  for (UInt_t i =0; i<19; i++) b_csvweights.push_back(Jet_SF_CSV[i]);
  if(m_isMC){
    b_btagweight = Jet_SF_CSV[0];
    b_btagweight_up = Jet_SF_CSV[1];
    b_btagweight_dn = Jet_SF_CSV[2];

    b_lf_up = Jet_SF_CSV[3];
    b_lf_dn = Jet_SF_CSV[4];

    b_hfstats1_up = Jet_SF_CSV[7];
    b_hfstats1_dn = Jet_SF_CSV[8];

    b_hfstats2_up = Jet_SF_CSV[9];
    b_hfstats2_dn = Jet_SF_CSV[10];
  }
  return jets;
}

vector<TParticle> nanoAnalysis::BJetSelection(vector<TParticle> leptons)
{
  vector<TParticle> bJets;
  for (UInt_t i = 0; i < nJet; i++)
  {
    if (Jet_btagCSVV2[i] < 0.8484) continue;
    if (Jet_pt[i] < 20) continue;
    if (abs(Jet_eta[i]) > 2.4) continue;
    
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);

    auto bjet = TParticle();
    bjet.SetMomentum(mom);
    bJets.push_back(bjet);
    b_CSVv2.push_back(Jet_btagCSVV2[i]);
  }
  return bJets;
}

bool nanoAnalysis::hasOverLap(TLorentzVector cand, vector<TParticle> objects, Float_t rad)
{
  for (auto obj: objects){
    TLorentzVector mom;
    obj.Momentum(mom);
    if (cand.DeltaR(mom) < rad){
      return true;
    }
  }
  return false;
}

Double_t nanoAnalysis::roccoR(TLorentzVector m, int &q, int &nGen, int &nTrackerLayers)
{
  Float_t u1 = gRandom->Rndm();
  Float_t u2 = gRandom->Rndm();
  if (!m_isMC){
    return m_rocCor->kScaleDT(q, m.Pt(), m.Eta(), m.Phi(), 0, 0);
  }
  else {
    if (nGen > -1){
      return m_rocCor->kScaleFromGenMC(q, m.Pt(), m.Eta(), m.Phi(),
				       nTrackerLayers, GenPart_pt[nGen],
				       u1, 0, 0);
    }
    else
      return m_rocCor->kScaleAndSmearMC(q, m.Pt(), m.Eta(), m.Phi(),
					nTrackerLayers, u1, u2, 0, 0);
  }
  return 1.0;
}
