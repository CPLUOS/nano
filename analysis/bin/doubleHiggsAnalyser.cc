#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <TString.h>
#include <TSystem.h>
#include <TRefArray.h>
#include <TMatrixDfwd.h>
#include <TVectorD.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include "TMinuit.h"
#include "TError.h"

#include "nano/analysis/interface/doubleHiggsAnalyser.h"

#include "nano/external/interface/HiggsnessTopness.h"
#include "nano/oxbridgekinetics/src/Mt2/Basic_Mt2_332_Calculator.h"
#include "nano/oxbridgekinetics/src/Mt2/ChengHanBisect_Mt2_332_Calculator.h"

using namespace std;

void doubleHiggsAnalyser::MakeOutputBranch(TTree *tree) {
  // MT2 variables
  tree->Branch("basic_MT2_332_bbll",&basic_MT2_332_bbll,"basic_MT2_332_bbll/F");
  tree->Branch("basic_MT2_332_blbl",&basic_MT2_332_blbl,"basic_MT2_332_blbl/F");
  tree->Branch("basic_MT2_332_b",&basic_MT2_332_b,"basic_MT2_332_b/F");
  tree->Branch("basic_MT2_332_l",&basic_MT2_332_l,"basic_MT2_332_l/F");
  tree->Branch("mT",&mT,"mT/F");
  // lepton kinematic variables
  tree->Branch("lep1", "TLorentzVector", &lepton1);
  tree->Branch("lep2", "TLorentzVector", &lepton2);
  tree->Branch("ll", "TLorentzVector", &leptonlepton); 
  tree->Branch("ll_M", &ll_M, "ll_M/F"); 
  tree->Branch("ll_deltaR", &ll_deltaR, "ll_deltaR/F");
  tree->Branch("ll_deltaPhi", &ll_deltaPhi, "ll_deltaPhi/F");
  // bottom kinematic variables 
  tree->Branch("bot1", "TLorentzVector", &bottom1);
  tree->Branch("bot2", "TLorentzVector", &bottom1);
  tree->Branch("bb", "TLorentzVector", &bottombottom);
  tree->Branch("bb_deltaR", &bb_deltaR, "&bb_deltaR/F");
  tree->Branch("bb_deltaPhi", &bb_deltaPhi, "&bb_deltaPhi/F");
  // lepton + bottom
  tree->Branch("bl11", "TLorentzVector", &bottomlepton11); 
  tree->Branch("bl12", "TLorentzVector", &bottomlepton12); 
  tree->Branch("bl21", "TLorentzVector", &bottomlepton21); 
  tree->Branch("bl22", "TLorentzVector", &bottomlepton22); 
  // lepton bottom kinematic variables
  tree->Branch("bl_deltaR", "vector<Float_t>", &bl_deltaR);
  tree->Branch("bl_min_deltaR", &bl_min_deltaR, "bl_min_deltaR/F");
  tree->Branch("bbll_deltaR", &bbll_deltaR, "bbll_delatR/F");
  tree->Branch("bbll_deltaPhi", &bbll_deltaPhi, "bbll_deltaPhi/F");
  // missing et
  tree->Branch("MET","TLorentzVector",&missing);
  tree->Branch("bbll","TLorentzVector",&bbll);
  tree->Branch("topness",&topness,"topness/F");
  tree->Branch("higgsness",&higgsness,"higgsness/F");

  tree->Branch("step", &step, "step/I");
  tree->Branch("tmva_bdtg_output", &tmva_bdtg_output, "tmva_bdtg_output/F");
}

void doubleHiggsAnalyser::SetTMVA(TString weight_file_path) {
  bdtg_reader = new TMVA::Reader();
  bdtg_reader->AddVariable("ll_deltaR", &ll_deltaR);
  bdtg_reader->AddVariable("ll.Pt()", &ll_Pt);
  bdtg_reader->AddVariable("ll.M()", &ll_M);
  bdtg_reader->AddVariable("bb_deltaR", &bb_deltaR);
  bdtg_reader->AddVariable("bb.Pt()", &bb_Pt);
  bdtg_reader->AddVariable("bb.M()", &bb_M);
  bdtg_reader->AddVariable("bl_min_deltaR", &bl_min_deltaR);
  bdtg_reader->AddVariable("bbll_deltaR", &bbll_deltaR);
  bdtg_reader->AddVariable("bbll_deltaPhi", &bbll_deltaPhi);
  bdtg_reader->AddVariable("mT", &mT);
  bdtg_reader->AddVariable("basic_MT2_332_bbll", &basic_MT2_332_bbll);
  bdtg_reader->BookMVA("BDTG", weight_file_path);
  tmva_flag = true;
}

void doubleHiggsAnalyser::SetOutput(TString output_file_name) {
  out_file = TFile::Open(output_file_name,"RECREATE");
  out_tree = new TTree("events","events"); 
}

void doubleHiggsAnalyser::SetBranchAddress() {
}

void doubleHiggsAnalyser::Initiate(TString output_file_name) {
  // set output file
  doubleHiggsAnalyser::SetOutput(output_file_name);
  // make output branch
  doubleHiggsAnalyser::MakeOutputBranch(out_tree);
  // set branch address
  doubleHiggsAnalyser::SetBranchAddress();
}

void doubleHiggsAnalyser::ResetVariables() {
  //// MT2 variables
  basic_MT2_332_bbll = std::numeric_limits<float>::quiet_NaN();
  basic_MT2_332_blbl = std::numeric_limits<float>::quiet_NaN();
  basic_MT2_332_b = std::numeric_limits<float>::quiet_NaN();
  basic_MT2_332_l = std::numeric_limits<float>::quiet_NaN();
  mT = std::numeric_limits<float>::quiet_NaN();
  //// lepton variables
  lepton1.Clear();
  lepton2.Clear();
  leptonlepton.Clear();
  ll_M = std::numeric_limits<float>::quiet_NaN();
  ll_Pt = std::numeric_limits<float>::quiet_NaN();
  ll_deltaR = std::numeric_limits<float>::quiet_NaN();
  ll_deltaPhi = std::numeric_limits<float>::quiet_NaN();
  leptons.clear();
  //// bottom variables
  bottom1.Clear();
  bottom2.Clear();
  bottombottom.Clear();
  bb_M = std::numeric_limits<float>::quiet_NaN();
  bb_Pt = std::numeric_limits<float>::quiet_NaN();
  bb_deltaR = std::numeric_limits<float>::quiet_NaN();
  bb_deltaPhi = std::numeric_limits<float>::quiet_NaN();
  bottoms.clear();
  ////lepton and bottom variables
  bottomlepton11.Clear();
  bottomlepton12.Clear();
  bottomlepton21.Clear();
  bottomlepton22.Clear();
  bbll.Clear();
  bl_deltaR.clear();
  bl_min_deltaR = std::numeric_limits<float>::quiet_NaN();
  bbll_deltaR = std::numeric_limits<float>::quiet_NaN();
  bbll_deltaPhi = std::numeric_limits<float>::quiet_NaN();
  //// MET variables
  missing.Clear();
  //// cut variables
  step = 0;
  //// higgsness and topness
  higgsness = std::numeric_limits<float>::quiet_NaN();
  topness = std::numeric_limits<float>::quiet_NaN();
  //// tmva variables
  tmva_bdtg_output = std::numeric_limits<float>::quiet_NaN();
}

bool doubleHiggsAnalyser::Analysis() {
    doubleHiggsAnalyser::ResetVariables();
    
    // Missing ET //
    if (MET_pt < 20) return false;
    missing.SetPtEtaPhiM(MET_pt,0,MET_phi,0);

    // map<Float_t, std::pair<int,int>, greater<Float_t>> leptons : map of <pt,index>:<K,V> of leptons sorted by pt.
  
    // Muon Selection //
    for (UInt_t i = 0; i < nMuon; i++) {
      if (Muon_pt[i] < 10) continue;
      if (Muon_eta[i] > 2.4) continue;
      leptons.insert(make_pair(Muon_pt[i],make_pair(13*Muon_charge[i],i)));
    }
    // Electron Selection //
    
    if (leptons.size()<2) return false;
    //// find two leptons with highest pt with opposite charge
    lepton_iter = leptons.begin();
    auto l1_info = lepton_iter->second;
    int pid1 = l1_info.first;
    int index_l1 = l1_info.second;
    lepton_iter++;
    auto l2_info = lepton_iter->second;
    int pid2 = l2_info.first;
    int index_l2 = l2_info.second;
    lepton_iter++;
    while (pid1*pid2 > 0 && (lepton_iter!=leptons.end())) {
      l2_info = lepton_iter->second;
      pid2 = l2_info.first;
      index_l2 = l2_info.second;
      lepton_iter++;
    }
    if (pid1*pid2 > 0) return false;
    lepton1.SetPtEtaPhiM(Muon_pt[index_l1], Muon_eta[index_l1], Muon_phi[index_l1], Muon_mass[index_l1]);
    lepton2.SetPtEtaPhiM(Muon_pt[index_l2], Muon_eta[index_l2], Muon_phi[index_l2], Muon_mass[index_l2]);
    leptonlepton = lepton1 + lepton2;
    ll_deltaR = fabs(lepton1.DeltaR(lepton2));
    ll_deltaPhi = fabs(lepton1.DeltaPhi(lepton2));
   
    //// Cuts on leptons 
    if (ll_deltaR < 0.07 || ll_deltaR > 3.3) {return true;} step++;
    if (leptonlepton.M() < 5 || leptonlepton.M() > 100) {return true;} step++;
    
    // Bottom Selection //
    for (UInt_t i = 0; i < nJet; i++) {
      if (Jet_pt[i] < 30) continue;
      if (Jet_eta[i] > 2.4) continue;
      if (Jet_jetId[i] < 1) continue;
      if (Jet_btagCSVV2[i] < 0.8484) continue; // medium cut
      bottoms.insert(make_pair(Jet_pt[i],i));
    }
    
    if (bottoms.size()<2) return false;
    //// find two leptons with highest pt with opposite charge
    bottom_iter = bottoms.begin();
    int index_b1 = bottom_iter->second;
    bottom_iter++;
    int index_b2 = bottom_iter->second;
    bottom1.SetPtEtaPhiM(Jet_pt[index_b1], Jet_eta[index_b1], Jet_phi[index_b1], Jet_mass[index_b1]);
    bottom2.SetPtEtaPhiM(Jet_pt[index_b2], Jet_eta[index_b2], Jet_phi[index_b2], Jet_mass[index_b2]);
    bottombottom = bottom1 + bottom2;
    bb_deltaR = fabs(bottom1.DeltaR(bottom2));
    bb_deltaPhi = fabs(bottom1.DeltaPhi(bottom2));
   
    //// Cuts on bottoms 
    if (bb_deltaR > 5) {return true;} step++;
    if (bottombottom.M() < 22) {return true;} step++;

    // TLorentzVector combinations
    //// 1 bottom and 1 lepton
    bottomlepton11 = bottom1 + lepton1;
    bottomlepton12 = bottom1 + lepton2;
    bottomlepton21 = bottom2 + lepton1;
    bottomlepton22 = bottom2 + lepton2;
    //// 2 bottoms and 2 leptons and 1 missingET object combinations
    bbll = bottombottom + leptonlepton;
    bl_deltaR.push_back(lepton1.DeltaR(bottom1));
    bl_deltaR.push_back(lepton2.DeltaR(bottom1));
    bl_deltaR.push_back(lepton1.DeltaR(bottom2));
    bl_deltaR.push_back(lepton2.DeltaR(bottom2));
    bl_min_deltaR = *min_element(begin(bl_deltaR),end(bl_deltaR));
    bbll_deltaR = leptonlepton.DeltaR(bottombottom);
    bbll_deltaPhi = leptonlepton.DeltaPhi(bottombottom);
    
    // MT2 variables
    mT = sqrt(2*leptonlepton.Pt()*missing.E()*(1-cos(leptonlepton.Phi()-missing.Phi())));
    //// MT2(bb+ll)
    Mt2::LorentzTransverseVector vis_A(Mt2::TwoVector(leptonlepton.Px(), leptonlepton.Py()), leptonlepton.M());
    Mt2::LorentzTransverseVector vis_B(Mt2::TwoVector(bottombottom.Px(),bottombottom.Py()), bottombottom.M());
    Mt2::TwoVector pT_Miss(missing.Px(), missing.Py());
    //// MT2(bl+bl)
    Mt2::LorentzTransverseVector vis_A_blbl(Mt2::TwoVector(bottomlepton12.Px(), bottomlepton12.Py()), bottomlepton12.M());
    Mt2::LorentzTransverseVector vis_B_blbl(Mt2::TwoVector(bottomlepton21.Px(), bottomlepton21.Py()), bottomlepton21.M());
    //// MT2(b)
    Mt2::LorentzTransverseVector vis_A_b(Mt2::TwoVector(bottom1.Px(), bottom1.Py()), bottom1.M());
    Mt2::LorentzTransverseVector vis_B_b(Mt2::TwoVector(bottom2.Px(), bottom2.Py()), bottom2.M());
    Mt2::TwoVector pT_Miss_b(leptonlepton.Px(), leptonlepton.Py());
    //// MT2(l)
    Mt2::LorentzTransverseVector vis_A_l(Mt2::TwoVector(lepton1.Px(), lepton1.Py()), lepton1.M());
    Mt2::LorentzTransverseVector vis_B_l(Mt2::TwoVector(lepton2.Px(), lepton2.Py()), lepton2.M());
    Mt2::TwoVector pT_Miss_l(missing.Px(), missing.Py());

    basic_MT2_332_bbll = basic_mt2_332Calculator.mt2_332(vis_A, vis_B, pT_Miss, missing.M());
    basic_MT2_332_blbl = basic_mt2_332Calculator.mt2_332(vis_A_blbl, vis_B_blbl, pT_Miss, missing.M());
    basic_MT2_332_b = basic_mt2_332Calculator.mt2_332(vis_A_b, vis_B_b, pT_Miss_b, leptonlepton.M());
    basic_MT2_332_l = basic_mt2_332Calculator.mt2_332(vis_A_l, vis_B_l, pT_Miss_l, missing.M());
   
    // global TLorentzVector variables for higgsness, topness library 
    g_lepton1 = lepton1;
    g_lepton2 = lepton2;
    g_bottom1 = bottom1;
    g_bottom2 = bottom2;
    g_missing = missing;
    // higgsness and topness
    higgsness = GetHiggsness();
    topness = GetTopness();

    // tmva variables
    ll_Pt = leptonlepton.Pt();
    ll_M = leptonlepton.M();
    bb_Pt = bottombottom.Pt();
    bb_M = bottombottom.M(); 
    if(tmva_flag) tmva_bdtg_output = bdtg_reader->EvaluateMVA("BDTG");
    
  return true;
}

void doubleHiggsAnalyser::Loop() {
  TTree *t = 0;
  bool keep = false;
  t = fChain;
  int nevents = t->GetEntries();
  int proc = 0; int temp = 0;
  int nevent = 0;
  for (int iev = 0; iev < nevents; iev++) {
    t->GetEntry(iev);
    keep = doubleHiggsAnalyser::Analysis();
    if (keep) out_tree->Fill();
    nevent++; temp = nevent*100/nevents;
    if (temp != proc) {
      proc ++;
      cout << "#######################" << endl;
      cout << " proceeding : " << proc << " %" << endl;
      cout << "#######################" << endl;
    }
  }
  doubleHiggsAnalyser::Finalize(); // Write the tree and Close the file.         
}

void doubleHiggsAnalyser::Finalize() {
  out_tree->Write();
  out_file->Close();
}

int main(Int_t argc, Char_t** argv)
{
    string env = getenv("CMSSW_BASE");
    string username = getenv("USER");

    if (argc != 1) {
       string dirName = "root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/"+username+"/nanoAOD/"+std::string(argv[1])+"/"+std::string(argv[2]);
       string temp = argv[2];
       Bool_t isMC = false;
       Size_t found = temp.find("Run");
       if (found == string::npos) isMC = true;
       for (Int_t i = 3; i<argc; i++) {
           TFile *f = TFile::Open(argv[i], "read");
           TTree *tree;
           f->GetObject("Events",tree);
           temp = argv[i];
           found = temp.find_last_of('/');
           string outPutName = dirName+temp.substr(found);
           doubleHiggsAnalyser ana(tree, isMC);
           ana.Initiate(outPutName);
           TString weight_file_path = "/home/sunyoung/nanoAOD/src/nano/analysis/bin/weight/DoubleHiggs_BDTG.weights.xml";
           ana.SetTMVA(weight_file_path);
           ana.Loop(); // Loop through the events and do doubleHIggsAnalyser::Analysis() per event.
        }
    } else {
        TChain *tree = new TChain("Events");
        //TString data_path = "/xrootd_UOS/store/group/nanoAOD/run2_2016v5/GluGluToHHTo2B2VTo2L2Nu_node_SM_13TeV-madgraph-v2/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/181108_163530/0000/*.root";
        TString data_path = "/xrootd_UOS/store/group/nanoAOD/run2_2016v5/TTJets_Dilept_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180607_120252/0000/*.root";
        tree->Add(data_path);

        tree->SetBranchStatus("*",true);
        // for nanoAOD format analysis
        doubleHiggsAnalyser ana(tree, true);
        //ana.Initiate("HH_SM.root"); // usage : Initiate("outputfilename.root")
        ana.Initiate("TT.root"); // usage : Initiate("outputfilename.root")
        TString weight_file_path = "/cms/ldap_home/sunyoung/nanoAOD/src/nano/analysis/bin/weight/DoubleHiggs_BDTG.weights.xml";
        //ana.SetTMVA(weight_file_path);
        //ana.Loop(); // Loop through the events and do doubleHIggsAnalyser::Analysis() per event.
    }
  return 0;
}
