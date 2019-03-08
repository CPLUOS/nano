#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

void addJetVariable(TMVA::DataLoader *dataloader) {
  dataloader->AddVariable("pt",    'F');
  dataloader->AddVariable("eta",   'F');
  dataloader->AddVariable("phi",   'F');
  dataloader->AddVariable("mass",  'F');
  dataloader->AddVariable("c_x1",  'F');
  dataloader->AddVariable("c_x2",  'F');
  dataloader->AddVariable("c_x3",  'F');
  dataloader->AddVariable("n_x1",  'F');
  dataloader->AddVariable("n_x2",  'F');
  dataloader->AddVariable("n_x3",  'F');
  dataloader->AddVariable("cmult", 'I');
  dataloader->AddVariable("nmult", 'I');
  dataloader->AddVariable("axis1", 'F');
  dataloader->AddVariable("axis2", 'F');
  dataloader->AddVariable("ptD",   'F');
  dataloader->AddVariable("area",  'F');
  dataloader->AddVariable("CSVV2", 'F'); 
}

void addHadVariablePP(TMVA::DataLoader *dataloader) {
  dataloader->AddVariable("Ks_d",                'F');
  dataloader->AddVariable("Ks_pt",               'F');
  dataloader->AddVariable("Ks_eta",              'F');
  dataloader->AddVariable("Ks_phi",              'F');
  dataloader->AddVariable("Ks_mass",             'F');
  dataloader->AddVariable("Ks_lxy",              'F');
  dataloader->AddVariable("Ks_lxySig",           'F');
  dataloader->AddVariable("Ks_l3D",              'F');
  dataloader->AddVariable("Ks_l3DSig",           'F');
  dataloader->AddVariable("Ks_legDR",            'F');
  dataloader->AddVariable("Ks_angleXY",          'F');
  dataloader->AddVariable("Ks_angleXYZ",         'F');
  dataloader->AddVariable("Ks_chi2",             'F');
  dataloader->AddVariable("Ks_dca",              'F');
  dataloader->AddVariable("Ks_dau1_chi2",        'F');
  dataloader->AddVariable("Ks_dau1_ipsigXY",     'F');
  dataloader->AddVariable("Ks_dau1_ipsigZ",      'F');
  dataloader->AddVariable("Ks_dau1_pt",          'F');
  dataloader->AddVariable("Ks_dau2_chi2",        'F');
  dataloader->AddVariable("Ks_dau2_ipsigXY",     'F');
  dataloader->AddVariable("Ks_dau2_ipsigZ",      'F');
  dataloader->AddVariable("Ks_dau2_pt",          'F');
  dataloader->AddVariable("Ks_bdt_score",         'F');
  dataloader->AddVariable("Ks_x", 'F');
}

void addTree(TMVA::DataLoader *dataloader, TTree* sig, TTree* bkg, Double_t sigW = 1.0, Double_t bkgW = 1.0) {
  dataloader->AddSignalTree    ( sig, sigW );
  dataloader->AddBackgroundTree( bkg, bkgW );
}
void addTree(TMVA::DataLoader *dataloader, TTree* sigTrain, TTree* sigTest, TTree* bkgTrain, TTree* bkgTest, Double_t sigWTrain = 1.0, Double_t sigWTest = 1.0, Double_t bkgWTrain = 1.0, Double_t bkgWTest = 1.0) {
  dataloader->AddSignalTree    ( sigTrain, sigWTrain, "Training");
  dataloader->AddSignalTree    ( sigTest,  sigWTest,  "Test");
  dataloader->AddBackgroundTree( bkgTrain, bkgWTrain, "Training");
  dataloader->AddBackgroundTree( bkgTest,  bkgWTest,  "Test");
}

void addMethod(std::map<std::string,int> Use, TMVA::Factory *factory, TMVA::DataLoader *dataloader) {
  if (Use["MLPBFGS"])
    factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );
  if (Use["TMlpANN"])
    factory->BookMethod( dataloader, TMVA::Types::kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  );
  if (Use["BDTG"]) // Gradient Boost
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG",
                         "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );
  if (Use["BDT"])  // Adaptive Boost
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT",
                         "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );
}

void setJKS(std::map<std::string,int> Use,
            TMVA::Factory *factory,
            TMVA::DataLoader *dataloader,
            TTree* sigTrain,
            TTree* bkgTrain,
            float signalWeight = 1,
            float backgroundWeight = 1,
            TCut cut_sig = "isSJet",
            TCut cut_bkg = "!isSJet",
            TString opt = "SplitMode=Random:NormMode=NumEvents:!V",
            float sig_ratio=0.7,
            float bkg_ratio=0.7) {
  TString opt_t = "nTrain_Signal="+std::to_string(int(sigTrain->Draw("",cut_sig,"")*sig_ratio))+":nTrain_Background="+std::to_string(int(bkgTrain->Draw("",cut_bkg,"")*bkg_ratio))+":"+opt;
  addJetVariable(dataloader);
  if (string(dataloader->GetName()).find("pp_") != std::string::npos) addHadVariablePP(dataloader);
  else addHadVariablePP(dataloader);
  addTree(dataloader, sigTrain, bkgTrain, signalWeight, backgroundWeight);
  dataloader->PrepareTrainingAndTestTree(cut_sig, cut_bkg, opt_t );
  addMethod(Use, factory, dataloader);
}
void setJKS(std::map<std::string,int> Use, 
            TMVA::Factory *factory, 
            TMVA::DataLoader *dataloader, 
            TTree* sigTrain, 
            TTree* sigTest, 
            TTree* bkgTrain, 
            TTree* bkgTest, 
            float signalWeight = 1, 
            float backgroundWeight = 1, 
            TCut cut_sig = "isSJet", 
            TCut cut_bkg = "!isSJet", 
            TString opt = "SplitMode=Random:NormMode=NumEvents:!V",
            float sig_ratio=0.7,
            float bkg_ratio=0.7) {
  TString opt_t = "nTrain_Signal="+std::to_string(int(sigTrain->Draw("",cut_sig,"")*sig_ratio))+":nTrain_Background="+std::to_string(int(bkgTrain->Draw("",cut_bkg,"")*bkg_ratio))+":"+opt;
  addJetVariable(dataloader);
  if (string(dataloader->GetName()).find("pp_") != std::string::npos) addHadVariablePP(dataloader);
  else addHadVariablePP(dataloader);
  addTree(dataloader, sigTrain, sigTest, bkgTrain, bkgTest, signalWeight, backgroundWeight);
  dataloader->PrepareTrainingAndTestTree(cut_sig, cut_bkg, opt_t );
  addMethod(Use, factory, dataloader);
}

void setJ(std::map<std::string,int> Use, 
            TMVA::Factory *factory, 
            TMVA::DataLoader *dataloader, 
            TTree* sigTrain, 
            TTree* bkgTrain, 
            float signalWeight = 1, 
            float backgroundWeight = 1, 
            TCut cut_sig = "isSJet", 
            TCut cut_bkg = "!isSJet", 
            TString opt = "SplitMode=Random:NormMode=NumEvents:!V",
            float sig_ratio=0.7,
            float bkg_ratio=0.7) {
  TString opt_t = "nTrain_Signal="+std::to_string(int(sigTrain->Draw("",cut_sig,"")*sig_ratio))+":nTrain_Background="+std::to_string(int(bkgTrain->Draw("",cut_bkg,"")*bkg_ratio))+":"+opt;
  std::cout << opt_t << std::endl;
  addJetVariable(dataloader);
  addTree(dataloader, sigTrain, bkgTrain, signalWeight, backgroundWeight);
  dataloader->PrepareTrainingAndTestTree(cut_sig, cut_bkg, opt_t );
  addMethod(Use, factory, dataloader);
}
void setJ(std::map<std::string,int> Use, 
            TMVA::Factory *factory, 
            TMVA::DataLoader *dataloader, 
            TTree* sigTrain, 
            TTree* sigTest,
            TTree* bkgTrain,
            TTree* bkgTest,
            float signalWeight = 1,
            float backgroundWeight = 1,
            TCut cut_sig = "isSJet",
            TCut cut_bkg = "!isSJet",
            TString opt = "SplitMode=Random:NormMode=NumEvents:!V",
            float sig_ratio=0.7,
            float bkg_ratio=0.7) {
  TString opt_t = "nTrain_Signal="+std::to_string(int(sigTrain->Draw("",cut_sig,"")*sig_ratio))+":nTrain_Background="+std::to_string(int(bkgTrain->Draw("",cut_bkg,"")*bkg_ratio))+":"+opt;
  addJetVariable(dataloader);
  addTree(dataloader, sigTrain, sigTest, bkgTrain, bkgTest, signalWeight, backgroundWeight);
  dataloader->PrepareTrainingAndTestTree(cut_sig, cut_bkg, opt_t );
  addMethod(Use, factory, dataloader);
}

void setHad(std::map<std::string,int> Use,
            TMVA::Factory *factory,
            TMVA::DataLoader *dataloader,
            TTree* sigTrain,
            TTree* bkgTrain,
            float signalWeight = 1,
            float backgroundWeight = 1,
            TCut cut_sig = "isSJet",
            TCut cut_bkg = "!isSJet",
            TString opt = "SplitMode=Random:NormMode=NumEvents:!V",
            float sig_ratio=0.7,
            float bkg_ratio=0.7) {
  TString opt_t = "nTrain_Signal="+std::to_string(int(sigTrain->Draw("",cut_sig,"")*sig_ratio))+":nTrain_Background="+std::to_string(int(bkgTrain->Draw("",cut_bkg,"")*bkg_ratio))+":"+opt;
  addHadVariablePP(dataloader);
  addTree(dataloader, sigTrain, bkgTrain, signalWeight, backgroundWeight);
  dataloader->PrepareTrainingAndTestTree(cut_sig, cut_bkg, opt_t );
  addMethod(Use, factory, dataloader);
}
void setHad(std::map<std::string,int> Use,
            TMVA::Factory *factory,
            TMVA::DataLoader *dataloader,
            TTree* sigTrain,
            TTree* sigTest,
            TTree* bkgTrain,
            TTree* bkgTest,
            float signalWeight = 1,
            float backgroundWeight = 1,
            TCut cut_sig = "isSJet",
            TCut cut_bkg = "!isSJet",
            TString opt = "SplitMode=Random:NormMode=NumEvents:!V",
            float sig_ratio=0.7,
            float bkg_ratio=0.7) {
  TString opt_t = "nTrain_Signal="+std::to_string(int(sigTrain->Draw("",cut_sig,"")*sig_ratio))+":nTrain_Background="+std::to_string(int(bkgTrain->Draw("",cut_bkg,"")*bkg_ratio))+":"+opt;
  addHadVariablePP(dataloader);
  addTree(dataloader, sigTrain, sigTest, bkgTrain, bkgTest, signalWeight, backgroundWeight);
  dataloader->PrepareTrainingAndTestTree(cut_sig, cut_bkg, opt_t );
  addMethod(Use, factory, dataloader);
}

int vts_dR_04_Jet( TString myMethodList = "" )
{
  // This loads the library
  TMVA::Tools::Instance();

  // Default MVA methods to be trained + tested
  std::map<std::string,int> Use;

  Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
  Use["TMlpANN"]         = 0; // ROOT's own ANN
  // Boosted Decision Trees
  Use["BDT"]             = 1; // uses Adaptive Boost
  Use["BDTG"]            = 0; // uses Gradient Boost

  std::map<std::string, int> Opt;

  //
  // ---------------------------------------------------------------

  std::cout << std::endl;
  std::cout << "==> Start vts_dR_04_Jet" << std::endl;

  // Select methods (don't look at this code - not of interest)
  if (myMethodList != "") {
    for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

    std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
    for (UInt_t i=0; i<mlist.size(); i++) {
      std::string regMethod(mlist[i]);
      if (Use.find(regMethod) == Use.end()) {
        std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
        for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
        std::cout << std::endl;
        return 1;
      }
      Use[regMethod] = 1;
    }
  }

  // --------------------------------------------------------------------------------------------------

  // Here the preparation phase begins

  // Read training and test data
  // (it is also possible to use ASCII format as s_vs_b_highest -> see TMVA Users Guide)
  TFile *bbars_pythia(0);       TFile *bbars_herwig(0);
  TFile *bsbar_pythia(0);       TFile *bsbar_herwig(0);
  TFile *bbbar_pythia(0);
  TFile *bbars_bsbar_pythia(0); TFile *bbars_bsbar_herwig(0);
  TFile *bbars_bsbar_bbbar_pythia(0);

//  TString sample_bsbar_pythia = "/xrootd/store/user/wjjang/test/sum_tt012j/bsbar_pythia.root";
//  TString sample_bsbar_herwig = "/xrootd/store/user/wjjang/test/sum_tt012j/bsbar_herwig.root";

  TString sample_bbars_pythia             = "/xrootd/store/user/wjjang/test/sum_tt012j/20190305/tt012j_bbars_2l_FxFx_sum_630_combined.root";
  TString sample_bsbar_pythia             = "/xrootd/store/user/wjjang/test/sum_tt012j/20190305/tt012j_bsbar_2l_FxFx_sum_774_combined.root";
  TString sample_bbbar_pythia             = "/xrootd/store/user/wjjang/test/sum_tt012j/20190305/tt012j_bbbar_2l_FxFx_sum_201.root";
  TString sample_bbars_bsbar_pythia       = "/xrootd/store/user/wjjang/test/sum_tt012j/20190305/tt012j_bbars_bsbar_sum_630_774_combined.root";
  TString sample_bbars_bsbar_bbbar_pythia = "/xrootd/store/user/wjjang/test/sum_tt012j/20190305/tt012j_bbars_bsbar_bbbar_sum_630_774_201_combined.root";

//  TString sample_bbars_herwig       = "/xrootd/store/user/wjjang/test/sum_tt012j/20190105/tt012j_bbars_2l_FxFx_herwigpp_sum_49.root";
//  TString sample_bsbar_herwig       = "/xrootd/store/user/wjjang/test/sum_tt012j/20190105/tt012j_bsbar_2l_FxFx_herwigpp_sum_49.root";
//  TString sample_bbars_bsbar_herwig = "/xrootd/store/user/wjjang/test/sum_tt012j/20190105/tt012j_bbars_bsbar_herwigpp_sum_49_49.root";

  bbars_pythia = TFile::Open( sample_bbars_pythia );           //  bbars_herwig = TFile::Open( sample_bbars_herwig );
  bsbar_pythia = TFile::Open( sample_bsbar_pythia );           //  bsbar_herwig = TFile::Open( sample_bsbar_herwig );
  bbbar_pythia = TFile::Open( sample_bbbar_pythia );
  bbars_bsbar_pythia = TFile::Open( sample_bbars_bsbar_pythia ); //bbars_bsbar_herwig = TFile::Open( sample_bbars_bsbar_herwig );
  bbars_bsbar_bbbar_pythia = TFile::Open( sample_bbars_bsbar_bbbar_pythia );

  std::cout << "--- vts_dR_04_Jet       : Using input file 1 : " << bbars_pythia->GetName() << std::endl;
//  std::cout << "--- vts_dR_04_Jet       : Using input file 2 : " << bbars_herwig->GetName() << std::endl;
  std::cout << "--- vts_dR_04_Jet       : Using input file 3 : " << bsbar_pythia->GetName() << std::endl;
//  std::cout << "--- vts_dR_04_Jet       : Using input file 4 : " << bsbar_herwig->GetName() << std::endl;
  std::cout << "--- vts_dR_04_Jet       : Using input file 5 : " << bbbar_pythia->GetName() << std::endl;
  std::cout << "--- vts_dR_04_Jet       : Using input file 6 : " << bbars_bsbar_pythia->GetName() << std::endl;
//  std::cout << "--- vts_dR_04_Jet       : Using input file 7 : " << bbars_bsbar_herwig->GetName() << std::endl;
  std::cout << "--- vts_dR_04_Jet       : Using input file 8 : " << bbars_bsbar_bbbar_pythia->GetName() << std::endl;

  // Register the training and test trees
  TTree *bbars_pythia_tree           = (TTree*)bbars_pythia->Get("MVA_jet");
  TTree *bsbar_pythia_tree           = (TTree*)bsbar_pythia->Get("MVA_jet");
  TTree *bbbar_pythia_tree           = (TTree*)bbbar_pythia->Get("MVA_jet");
  TTree *bbars_bsbar_pythia_tree     = (TTree*)bbars_bsbar_pythia->Get("MVA_jet");
  TTree *bbars_bsbar_bbbar_pythia_tree     = (TTree*)bbars_bsbar_bbbar_pythia->Get("MVA_jet");
  
//  TTree *bbars_herwig_tree           = (TTree*)bbars_herwig->Get("MVA_jet");
//  TTree *bsbar_herwig_tree           = (TTree*)bsbar_herwig->Get("MVA_jet");
//  TTree *bbars_bsbar_herwig_tree     = (TTree*)bbars_bsbar_herwig->Get("MVA_jet");

  // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
  TString outfileName( "./output/vts_dR_04_Jet.root" );
//  TString outfileName( "vts_dR_04_Jet_20181022.root" );   
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

  TMVA::Factory *factory = new TMVA::Factory( "vts_dR_04_Jet", outputFile,
//                                              "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
                                              "!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification" );
  /* JKS Events */
  TMVA::DataLoader *pp_combined_JKS_BDT_highest = new TMVA::DataLoader("pp_combined_JKS_BDT_highest");
  TMVA::DataLoader *pp_combined_JKS_BDT_closest = new TMVA::DataLoader("pp_combined_JKS_BDT_closest");
  TMVA::DataLoader *pp_combined_JKS_BDT         = new TMVA::DataLoader("pp_combined_JKS_BDT");
  TMVA::DataLoader *pp_combined_JKS_BDT_highBDT = new TMVA::DataLoader("pp_combined_JKS_BDT_highBDT");

  TMVA::DataLoader *pp_s_vs_b_JKS_BDT_highest = new TMVA::DataLoader("pp_s_vs_b_JKS_BDT_highest");
  TMVA::DataLoader *pp_s_vs_b_JKS_BDT_closest = new TMVA::DataLoader("pp_s_vs_b_JKS_BDT_closest");
  TMVA::DataLoader *pp_s_vs_b_JKS_BDT         = new TMVA::DataLoader("pp_s_vs_b_JKS_BDT");
  TMVA::DataLoader *pp_s_vs_b_JKS_BDT_highBDT = new TMVA::DataLoader("pp_s_vs_b_JKS_BDT_highBDT");

  /* JKS Events wo KS */
  TMVA::DataLoader *JKS_wo_KS_highest        = new TMVA::DataLoader("JKS_wo_KS_highest");
  TMVA::DataLoader *JKS_wo_KS_highBDT        = new TMVA::DataLoader("JKS_wo_KS_highBDT");

  TMVA::DataLoader *JKS_wo_KS_s_vs_b_highest = new TMVA::DataLoader("JKS_wo_KS_s_vs_b_highest");
  TMVA::DataLoader *JKS_wo_KS_s_vs_b_highBDT = new TMVA::DataLoader("JKS_wo_KS_s_vs_b_highBDT");

  /* JKS Events wo Jet */
  TMVA::DataLoader *JKS_wo_Jet_highest        = new TMVA::DataLoader("JKS_wo_Jet_highest");
  TMVA::DataLoader *JKS_wo_Jet_highBDT        = new TMVA::DataLoader("JKS_wo_Jet_highBDT");

  TMVA::DataLoader *JKS_wo_Jet_s_vs_b_highest = new TMVA::DataLoader("JKS_wo_Jet_s_vs_b_highest");
  TMVA::DataLoader *JKS_wo_Jet_s_vs_b_highBDT = new TMVA::DataLoader("JKS_wo_Jet_s_vs_b_highBDT");

  /* Jet Events */
  TMVA::DataLoader *pp_combined_J_BDT_highest = new TMVA::DataLoader("pp_combined_J_BDT_highest");
  TMVA::DataLoader *pp_combined_J_BDT_closest = new TMVA::DataLoader("pp_combined_J_BDT_closest");
  TMVA::DataLoader *pp_combined_J_BDT         = new TMVA::DataLoader("pp_combined_J_BDT");
  TMVA::DataLoader *pp_combined_J_BDT_highBDT = new TMVA::DataLoader("pp_combined_J_BDT_highBDT");

  TMVA::DataLoader *pp_s_vs_b_J_BDT_highest = new TMVA::DataLoader("pp_s_vs_b_J_BDT_highest");
  TMVA::DataLoader *pp_s_vs_b_J_BDT_closest = new TMVA::DataLoader("pp_s_vs_b_J_BDT_closest");
  TMVA::DataLoader *pp_s_vs_b_J_BDT         = new TMVA::DataLoader("pp_s_vs_b_J_BDT");
  TMVA::DataLoader *pp_s_vs_b_J_BDT_highBDT = new TMVA::DataLoader("pp_s_vs_b_J_BDT_highBDT");

  /* All Events */
  TMVA::DataLoader *all_highest        = new TMVA::DataLoader("all_highest");
  TMVA::DataLoader *all_highBDT        = new TMVA::DataLoader("all_highBDT");

  TMVA::DataLoader *all_s_vs_b_highest = new TMVA::DataLoader("all_s_vs_b_highest");
  TMVA::DataLoader *all_s_vs_b_highBDT = new TMVA::DataLoader("all_s_vs_b_highBDT");

  /* All Events wo KS */
  TMVA::DataLoader *all_wo_KS_highest        = new TMVA::DataLoader("all_wo_KS_highest");
  TMVA::DataLoader *all_wo_KS_highBDT        = new TMVA::DataLoader("all_wo_KS_highBDT");

  TMVA::DataLoader *all_wo_KS_s_vs_b_highest = new TMVA::DataLoader("all_wo_KS_s_vs_b_highest");
  TMVA::DataLoader *all_wo_KS_s_vs_b_highBDT = new TMVA::DataLoader("all_wo_KS_s_vs_b_highBDT");

  //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
  //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";
  (TMVA::gConfig().GetVariablePlotting()).fMaxNumOfAllowedVariablesForScatterPlots = 99999;
  (TMVA::gConfig().GetVariablePlotting()).fNbinsXOfROCCurve = 1000;

  // global event weights per tree (see below for setting event-wise weights) 
     Double_t signalWeight     = 1.0; Double_t backgroundWeight = 1.0;

  /* Cuts for JKS >>>> Ks_bdt_score >= 0.0885 : jet includes Real KS || 20190304 bbars_bsbar_bbbar sample*/ 
  TCut cut_JKS_BDT_s_highest   = "isHighest && Ks_bdt_score >= 0.0885 && isSJet";
  TCut cut_JKS_BDT_b_highest   = "isHighest && Ks_bdt_score >= 0.0885 && isBJet";
  TCut cut_JKS_BDT_all_highest = "isHighest && Ks_bdt_score >= 0.0885 && !isSJet";

  TCut cut_JKS_BDT_s_closest   = "isClosestToLep && Ks_bdt_score >= 0.0885 && isSJet";
  TCut cut_JKS_BDT_b_closest   = "isClosestToLep && Ks_bdt_score >= 0.0885 && isBJet";
  TCut cut_JKS_BDT_all_closest = "isClosestToLep && Ks_bdt_score >= 0.0885 && !isSJet";

  TCut cut_JKS_BDT_s_highBDT   = "isHighBDT && Ks_bdt_score >= 0.0885 && isSJet";
  TCut cut_JKS_BDT_b_highBDT   = "isHighBDT && Ks_bdt_score >= 0.0885 && isBJet";
  TCut cut_JKS_BDT_all_highBDT = "isHighBDT && Ks_bdt_score >= 0.0885 && !isSJet";

  TCut cut_JKS_BDT_s           = "Ks_bdt_score >= 0.0885 && isSJet";
  TCut cut_JKS_BDT_b           = "Ks_bdt_score >= 0.0885 && isBJet";
  TCut cut_JKS_BDT_all         = "Ks_bdt_score >= 0.0885 && !isSJet";

  /* Cuts for J >>>> Ks_bdt_score < 0.0885 : jet includes Fake KS */ 
  TCut cut_J_BDT_s_highest   = "isHighest && Ks_bdt_score < 0.0885 && isSJet";
  TCut cut_J_BDT_b_highest   = "isHighest && Ks_bdt_score < 0.0885 && isBJet";
  TCut cut_J_BDT_all_highest = "isHighest && Ks_bdt_score < 0.0885 && !isSJet";

  TCut cut_J_BDT_s_closest   = "isClosestToLep && Ks_bdt_score < 0.0885 && isSJet";
  TCut cut_J_BDT_b_closest   = "isClosestToLep && Ks_bdt_score < 0.0885 && isBJet";
  TCut cut_J_BDT_all_closest = "isClosestToLep && Ks_bdt_score < 0.0885 && !isSJet";

  TCut cut_J_BDT_s_highBDT   = "isHighBDT && Ks_bdt_score < 0.0885 && isSJet";
  TCut cut_J_BDT_b_highBDT   = "isHighBDT && Ks_bdt_score < 0.0885 && isBJet";
  TCut cut_J_BDT_all_highBDT = "isHighBDT && Ks_bdt_score < 0.0885 && !isSJet";

  TCut cut_J_BDT_s           = "Ks_bdt_score < 0.0885 && isSJet";
  TCut cut_J_BDT_b           = "Ks_bdt_score < 0.0885 && isBJet";
  TCut cut_J_BDT_all         = "Ks_bdt_score < 0.0885 && !isSJet";

  TCut cut_s_highest   = "isHighest && isSJet";
  TCut cut_b_highest   = "isHighest && isBJet";
  TCut cut_non_highest = "isHighest && !isSJet && !isBJet";
  TCut cut_all_highest = "isHighest && !isSJet";

  TCut cut_s_closest   = "isClosestToLep && isSJet";
  TCut cut_b_closest   = "isClosestToLep && isBJet";
  TCut cut_non_closest = "isClosestToLep && !isSJet && !isBJet";
  TCut cut_all_closest = "isClosestToLep && !isSJet";

  TCut cut_s_highBDT   = "isHighBDT && isSJet";
  TCut cut_b_highBDT   = "isHighBDT && isBJet";
  TCut cut_non_highBDT = "isHighBDT && !isSJet && !isBJet";
  TCut cut_all_highBDT = "isHighBDT && !isSJet";

  TCut cut_s           = "isSJet";
  TCut cut_b           = "isBJet";
  TCut cut_all         = "!isSJet";

  /* JKS Events */
  setJKS(Use, factory, pp_combined_JKS_BDT_highest, bbars_bsbar_bbbar_pythia_tree, bbars_bsbar_bbbar_pythia_tree, signalWeight, backgroundWeight, cut_JKS_BDT_s_highest, cut_JKS_BDT_all_highest);
  setJKS(Use, factory, pp_combined_JKS_BDT_closest, bbars_bsbar_bbbar_pythia_tree, bbars_bsbar_bbbar_pythia_tree, signalWeight, backgroundWeight, cut_JKS_BDT_s_closest, cut_JKS_BDT_all_closest);
  setJKS(Use, factory, pp_s_vs_b_JKS_BDT_highest,   bbars_bsbar_bbbar_pythia_tree, bbars_bsbar_bbbar_pythia_tree, signalWeight, backgroundWeight, cut_JKS_BDT_s_highest, cut_JKS_BDT_b_highest);
  setJKS(Use, factory, pp_s_vs_b_JKS_BDT_closest,   bbars_bsbar_bbbar_pythia_tree, bbars_bsbar_bbbar_pythia_tree, signalWeight, backgroundWeight, cut_JKS_BDT_s_closest, cut_JKS_BDT_b_closest);
  /* JKS Events wo KS */
  setJ(  Use, factory, JKS_wo_KS_highest,           bbars_bsbar_bbbar_pythia_tree, bbars_bsbar_bbbar_pythia_tree, signalWeight, backgroundWeight, cut_JKS_BDT_s_highest, cut_JKS_BDT_all_highest);
  setJ(  Use, factory, JKS_wo_KS_s_vs_b_highest,    bbars_bsbar_bbbar_pythia_tree, bbars_bsbar_bbbar_pythia_tree, signalWeight, backgroundWeight, cut_JKS_BDT_s_highest, cut_JKS_BDT_b_highest);
  /* JKS Events wo Jet */
  setHad(Use, factory, JKS_wo_Jet_highest,          bbars_bsbar_bbbar_pythia_tree, bbars_bsbar_bbbar_pythia_tree, signalWeight, backgroundWeight, cut_JKS_BDT_s_highest, cut_JKS_BDT_all_highest);
  setHad(Use, factory, JKS_wo_Jet_s_vs_b_highest,   bbars_bsbar_bbbar_pythia_tree, bbars_bsbar_bbbar_pythia_tree, signalWeight, backgroundWeight, cut_JKS_BDT_s_highest, cut_JKS_BDT_b_highest);
  /* Jet Events */
  setJ(  Use, factory, pp_combined_J_BDT_highest,   bbars_bsbar_bbbar_pythia_tree, bbars_bsbar_bbbar_pythia_tree, signalWeight, backgroundWeight, cut_J_BDT_s_highest,   cut_J_BDT_all_highest);
  setJ(  Use, factory, pp_combined_J_BDT_closest,   bbars_bsbar_bbbar_pythia_tree, bbars_bsbar_bbbar_pythia_tree, signalWeight, backgroundWeight, cut_J_BDT_s_closest,   cut_J_BDT_all_closest);
  setJ(  Use, factory, pp_s_vs_b_J_BDT_highest,     bbars_bsbar_bbbar_pythia_tree, bbars_bsbar_bbbar_pythia_tree, signalWeight, backgroundWeight, cut_J_BDT_s_highest,   cut_J_BDT_b_highest);
  setJ(  Use, factory, pp_s_vs_b_J_BDT_closest,     bbars_bsbar_bbbar_pythia_tree, bbars_bsbar_bbbar_pythia_tree, signalWeight, backgroundWeight, cut_J_BDT_s_closest,   cut_J_BDT_b_closest);
  /* All Events */
  setJKS(Use, factory, all_highest,                 bbars_bsbar_bbbar_pythia_tree, bbars_bsbar_bbbar_pythia_tree, signalWeight, backgroundWeight, cut_s_highest,         cut_all_highest);
  setJKS(Use, factory, all_s_vs_b_highest,          bbars_bsbar_bbbar_pythia_tree, bbars_bsbar_bbbar_pythia_tree, signalWeight, backgroundWeight, cut_s_highest,         cut_b_highest);
  /* All Events wo KS */
  setJ(  Use, factory, all_wo_KS_highest,           bbars_bsbar_bbbar_pythia_tree, bbars_bsbar_bbbar_pythia_tree, signalWeight, backgroundWeight, cut_s_highest,         cut_all_highest);
  setJ(  Use, factory, all_wo_KS_s_vs_b_highest,    bbars_bsbar_bbbar_pythia_tree, bbars_bsbar_bbbar_pythia_tree, signalWeight, backgroundWeight, cut_s_highest,         cut_b_highest);

  // For an example of the category classifier usage, see: vts_dR_04_JetCategory
  //
  // --------------------------------------------------------------------------------------------------
  //  Now you can optimize the setting (configuration) of the MVAs using the set of training events
  // STILL EXPERIMENTAL and only implemented for BDT's !
  //
  //     factory->OptimizeAllMethods("SigEffAt001","Scan");
  //     factory->OptimizeAllMethods("ROCIntegral","FitGA");
  //
  // --------------------------------------------------------------------------------------------------

  // Now you can tell the factory to train, test, and evaluate the MVAs
  //
  // Train MVAs using the set of training events
  factory->TrainAllMethods();

  // Evaluate all MVAs using the set of test events
  factory->TestAllMethods();

  // Evaluate and compare performance of all configured MVAs
  factory->EvaluateAllMethods();

  // --------------------------------------------------------------

  // Save the output
  outputFile->Close();

  std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
  std::cout << "==> vts_dR_04_Jet is done!" << std::endl;

  delete factory;

  // Launch the GUI for the root macros
  if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );

  return 0;
}

int main( int argc, char** argv )
{
  // Select methods (don't look at this code - not of interest)
  TString methodList;
  for (int i=1; i<argc; i++) {
    TString regMethod(argv[i]);
    if(regMethod=="-b" || regMethod=="--batch") continue;
    if (!methodList.IsNull()) methodList += TString(",");
    methodList += regMethod;
  }
  return vts_dR_04_Jet(methodList);
}
