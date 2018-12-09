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
  dataloader->AddVariable("KS_d",                'F');
  dataloader->AddVariable("KS_pt",               'F');
  dataloader->AddVariable("KS_eta",              'F');
  dataloader->AddVariable("KS_phi",              'F');
  dataloader->AddVariable("KS_mass",             'F');
  dataloader->AddVariable("KS_lxy",              'F');
  dataloader->AddVariable("KS_lxySig",           'F');
  dataloader->AddVariable("KS_l3D",              'F');
  dataloader->AddVariable("KS_l3DSig",           'F');
  dataloader->AddVariable("KS_legDR",            'F');
  dataloader->AddVariable("KS_angleXY",          'F');
  dataloader->AddVariable("KS_angleXYZ",         'F');
  dataloader->AddVariable("KS_chi2",             'F');
  dataloader->AddVariable("KS_dca",              'F');
  dataloader->AddVariable("KS_dau1_chi2",        'F');
  dataloader->AddVariable("KS_dau1_ipsigXY",     'F');
  dataloader->AddVariable("KS_dau1_ipsigZ",      'F');
  dataloader->AddVariable("KS_dau1_pt",          'F');
  dataloader->AddVariable("KS_dau2_chi2",        'F');
  dataloader->AddVariable("KS_dau2_ipsigXY",     'F');
  dataloader->AddVariable("KS_dau2_ipsigZ",      'F');
  dataloader->AddVariable("KS_dau2_pt",          'F');
  dataloader->AddVariable("KS_best_bdt",         'F');
  dataloader->AddVariable("KS_x", 'F');
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
            TString opt = "SplitMode=Random:NormMode=NumEvents:!V") {
  addJetVariable(dataloader);
  if (string(dataloader->GetName()).find("pp_") != std::string::npos) addHadVariablePP(dataloader);
  addTree(dataloader, sigTrain, bkgTrain, signalWeight, backgroundWeight);
  dataloader->PrepareTrainingAndTestTree(cut_sig, cut_bkg, opt );
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
            TString opt = "SplitMode=Random:NormMode=NumEvents:!V") {
  addJetVariable(dataloader);
  if (string(dataloader->GetName()).find("pp_") != std::string::npos) addHadVariablePP(dataloader);
  addTree(dataloader, sigTrain, sigTest, bkgTrain, bkgTest, signalWeight, backgroundWeight);
  dataloader->PrepareTrainingAndTestTree(cut_sig, cut_bkg, opt );
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
            TString opt = "SplitMode=Random:NormMode=NumEvents:!V") {
  addJetVariable(dataloader);
  addTree(dataloader, sigTrain, bkgTrain, signalWeight, backgroundWeight);
  dataloader->PrepareTrainingAndTestTree(cut_sig, cut_bkg, opt );
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
            TString opt = "SplitMode=Random:NormMode=NumEvents:!V") {
  addJetVariable(dataloader);
  addTree(dataloader, sigTrain, sigTest, bkgTrain, bkgTest, signalWeight, backgroundWeight);
  dataloader->PrepareTrainingAndTestTree(cut_sig, cut_bkg, opt );
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
  Use["BDTG"]            = 1; // uses Gradient Boost

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

  TString sample_bbars_pythia             = "/xrootd/store/user/wjjang/test/sum_tt012j/20181204/tt012j_bbars_2l_FxFx_sum_630_combined.root";
  TString sample_bsbar_pythia             = "/xrootd/store/user/wjjang/test/sum_tt012j/20181204/tt012j_bsbar_2l_FxFx_sum_774_combined.root";
  TString sample_bbbar_pythia             = "/xrootd/store/user/wjjang/test/sum_tt012j/20181204/tt012j_bbbar_2l_FxFx_sum_201.root";
  TString sample_bbars_bsbar_pythia       = "/xrootd/store/user/wjjang/test/sum_tt012j/20181204/tt012j_bbars_bsbar_sum_630_774_combined.root";
  TString sample_bbars_bsbar_bbbar_pythia = "/xrootd/store/user/wjjang/test/sum_tt012j/20181204/tt012j_bbars_bsbar_bbbar_sum_630_774_201_combined.root";

  TString sample_bbars_herwig       = "/xrootd/store/user/wjjang/test/sum_tt012j/20181204/tt012j_bbars_2l_FxFx_herwigpp_sum_49.root";
  TString sample_bsbar_herwig       = "/xrootd/store/user/wjjang/test/sum_tt012j/20181204/tt012j_bsbar_2l_FxFx_herwigpp_sum_49.root";
  TString sample_bbars_bsbar_herwig = "/xrootd/store/user/wjjang/test/sum_tt012j/20181204/tt012j_bbars_bsbar_herwigpp_sum_49_49.root";

  bbars_pythia = TFile::Open( sample_bbars_pythia );             bbars_herwig = TFile::Open( sample_bbars_herwig );
  bsbar_pythia = TFile::Open( sample_bsbar_pythia );             bsbar_herwig = TFile::Open( sample_bsbar_herwig );
  bbbar_pythia = TFile::Open( sample_bbbar_pythia );
  bbars_bsbar_pythia = TFile::Open( sample_bbars_bsbar_pythia ); bbars_bsbar_herwig = TFile::Open( sample_bbars_bsbar_herwig );
  bbars_bsbar_bbbar_pythia = TFile::Open( sample_bbars_bsbar_bbbar_pythia );

  std::cout << "--- vts_dR_04_Jet       : Using input file 1 : " << bbars_pythia->GetName() << std::endl;
  std::cout << "--- vts_dR_04_Jet       : Using input file 2 : " << bbars_herwig->GetName() << std::endl;
  std::cout << "--- vts_dR_04_Jet       : Using input file 3 : " << bsbar_pythia->GetName() << std::endl;
  std::cout << "--- vts_dR_04_Jet       : Using input file 4 : " << bsbar_herwig->GetName() << std::endl;
  std::cout << "--- vts_dR_04_Jet       : Using input file 5 : " << bbbar_pythia->GetName() << std::endl;
  std::cout << "--- vts_dR_04_Jet       : Using input file 6 : " << bbars_bsbar_pythia->GetName() << std::endl;
  std::cout << "--- vts_dR_04_Jet       : Using input file 7 : " << bbars_bsbar_herwig->GetName() << std::endl;
  std::cout << "--- vts_dR_04_Jet       : Using input file 8 : " << bbars_bsbar_bbbar_pythia->GetName() << std::endl;

  // Register the training and test trees
  TTree *bbars_pythia_tree           = (TTree*)bbars_pythia->Get("MVA_jet");
  TTree *bsbar_pythia_tree           = (TTree*)bsbar_pythia->Get("MVA_jet");
  TTree *bbbar_pythia_tree           = (TTree*)bbbar_pythia->Get("MVA_jet");
  TTree *bbars_bsbar_pythia_tree     = (TTree*)bbars_bsbar_pythia->Get("MVA_jet");
  TTree *bbars_bsbar_bbbar_pythia_tree     = (TTree*)bbars_bsbar_bbbar_pythia->Get("MVA_jet");
  
  TTree *bbars_herwig_tree           = (TTree*)bbars_herwig->Get("MVA_jet");
  TTree *bsbar_herwig_tree           = (TTree*)bsbar_herwig->Get("MVA_jet");
  TTree *bbars_bsbar_herwig_tree     = (TTree*)bbars_bsbar_herwig->Get("MVA_jet");

  // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
  TString outfileName( "./output/vts_dR_04_Jet.root" );
//  TString outfileName( "vts_dR_04_Jet_20181022.root" );   
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

  TMVA::Factory *factory = new TMVA::Factory( "vts_dR_04_Jet", outputFile,
//                                              "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
                                              "!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification" );
  /* JKS BDT */
  TMVA::DataLoader *pp_combined_JKS_BDT_highest = new TMVA::DataLoader("pp_combined_JKS_BDT_highest");
  TMVA::DataLoader *pp_combined_JKS_BDT_closest = new TMVA::DataLoader("pp_combined_JKS_BDT_closest");
  TMVA::DataLoader *pp_combined_JKS_BDT         = new TMVA::DataLoader("pp_combined_JKS_BDT");

  TMVA::DataLoader *pp_s_vs_b_JKS_BDT_highest = new TMVA::DataLoader("pp_s_vs_b_JKS_BDT_highest");
  TMVA::DataLoader *pp_s_vs_b_JKS_BDT_closest = new TMVA::DataLoader("pp_s_vs_b_JKS_BDT_closest");
  TMVA::DataLoader *pp_s_vs_b_JKS_BDT         = new TMVA::DataLoader("pp_s_vs_b_JKS_BDT");

  TMVA::DataLoader *pp_s_vs_b_with_KS_highest   = new TMVA::DataLoader("pp_s_vs_bwith_KS_highest");
  TMVA::DataLoader *pp_s_vs_non_with_KS_highest = new TMVA::DataLoader("pp_s_vs_non_with_KS_highest");
  TMVA::DataLoader *pp_s_vs_all_with_KS_highest = new TMVA::DataLoader("pp_s_vs_all_with_KS_highest");
  TMVA::DataLoader *pp_s_vs_b_with_KS_closest   = new TMVA::DataLoader("pp_s_vs_b_with_KS_closest");
  TMVA::DataLoader *pp_s_vs_non_with_KS_closest = new TMVA::DataLoader("pp_s_vs_non_with_KS_closest");
  TMVA::DataLoader *pp_s_vs_all_with_KS_closest = new TMVA::DataLoader("pp_s_vs_all_with_KS_closest");
  TMVA::DataLoader *pp_s_vs_b_with_KS           = new TMVA::DataLoader("pp_s_vs_b_with_KS");
  TMVA::DataLoader *pp_s_vs_non_with_KS         = new TMVA::DataLoader("pp_s_vs_non_with_KS");
  TMVA::DataLoader *pp_s_vs_all_with_KS         = new TMVA::DataLoader("pp_s_vs_all_with_KS");

  /* J BDT */
  TMVA::DataLoader *pp_combined_J_BDT_highest = new TMVA::DataLoader("pp_combined_J_BDT_highest");
  TMVA::DataLoader *pp_combined_J_BDT_closest = new TMVA::DataLoader("pp_combined_J_BDT_closest");
  TMVA::DataLoader *pp_combined_J_BDT         = new TMVA::DataLoader("pp_combined_J_BDT");

  TMVA::DataLoader *pp_s_vs_b_J_BDT_highest = new TMVA::DataLoader("pp_s_vs_b_J_BDT_highest");
  TMVA::DataLoader *pp_s_vs_b_J_BDT_closest = new TMVA::DataLoader("pp_s_vs_b_J_BDT_closest");
  TMVA::DataLoader *pp_s_vs_b_J_BDT         = new TMVA::DataLoader("pp_s_vs_b_J_BDT");

  TMVA::DataLoader *pp_s_vs_b_highest         = new TMVA::DataLoader("pp_s_vs_b_highest");
  TMVA::DataLoader *pp_s_vs_non_highest       = new TMVA::DataLoader("pp_s_vs_non_highest");
  TMVA::DataLoader *pp_s_vs_all_highest       = new TMVA::DataLoader("pp_s_vs_all_highest");
  TMVA::DataLoader *pp_s_vs_b_closest         = new TMVA::DataLoader("pp_s_vs_b_closest");
  TMVA::DataLoader *pp_s_vs_non_closest       = new TMVA::DataLoader("pp_s_vs_non_closest");
  TMVA::DataLoader *pp_s_vs_all_closest       = new TMVA::DataLoader("pp_s_vs_all_closest");
  TMVA::DataLoader *pp_s_vs_b                 = new TMVA::DataLoader("pp_s_vs_b");
  TMVA::DataLoader *pp_s_vs_non               = new TMVA::DataLoader("pp_s_vs_non");
  TMVA::DataLoader *pp_s_vs_all               = new TMVA::DataLoader("pp_s_vs_all");

  //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
  //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";
  (TMVA::gConfig().GetVariablePlotting()).fMaxNumOfAllowedVariablesForScatterPlots = 99999;
  (TMVA::gConfig().GetVariablePlotting()).fNbinsXOfROCCurve = 1000;

  // global event weights per tree (see below for setting event-wise weights) 
     Double_t signalWeight     = 1.0; Double_t backgroundWeight = 1.0;

  /* Cuts for JKS >>>> KS_best_bdt >= 0.0783 : jet includes Real KS || 20181120 bbars_bsbar_bbbar sample*/ 
  TCut cut_JKS_BDT_s_highest   = "isHighest && KS_idx != -99 && KS_best_bdt >= 0.0783 && isSJet && !isBJet";
  TCut cut_JKS_BDT_b_highest   = "isHighest && KS_idx != -99 && KS_best_bdt >= 0.0783 && isBJet && !isSJet";
  TCut cut_JKS_BDT_all_highest = "isHighest && KS_idx != -99 && KS_best_bdt >= 0.0783 && !isSJet";

  TCut cut_JKS_BDT_s_closest   = "isClosestToLep && KS_idx != -99 && KS_best_bdt >= 0.0783 && isSJet && !isBJet";
  TCut cut_JKS_BDT_b_closest   = "isClosestToLep && KS_idx != -99 && KS_best_bdt >= 0.0783 && isBJet && !isSJet";
  TCut cut_JKS_BDT_all_closest = "isClosestToLep && KS_idx != -99 && KS_best_bdt >= 0.0783 && !isSJet";

  TCut cut_JKS_BDT_s           = "KS_idx != -99 && KS_best_bdt >= 0.0783 && isSJet && !isBJet";
  TCut cut_JKS_BDT_b           = "KS_idx != -99 && KS_best_bdt >= 0.0783 && isBJet && !isSJet";
  TCut cut_JKS_BDT_all         = "KS_idx != -99 && KS_best_bdt >= 0.0783 && !isSJet";

  /* Cuts for J >>>> KS_best_bdt < 0.0783 : jet includes Fake KS | KS_idx == -99 : KS isn't in the jet or there is no KS in a event */ 
  TCut cut_J_BDT_s_highest   = "isHighest && (KS_idx == -99 || ( KS_best_bdt >= -1. && KS_best_bdt < 0.0783 ) ) && isSJet && !isBJet";
  TCut cut_J_BDT_b_highest   = "isHighest && (KS_idx == -99 || ( KS_best_bdt >= -1. && KS_best_bdt < 0.0783 ) ) && isBJet && !isSJet";
  TCut cut_J_BDT_all_highest = "isHighest && (KS_idx == -99 || ( KS_best_bdt >= -1. && KS_best_bdt < 0.0783 ) ) && !isSJet";

  TCut cut_J_BDT_s_closest   = "isClosestToLep && (KS_idx == -99 || ( KS_best_bdt >= -1. && KS_best_bdt < 0.0783 ) ) && isSJet && !isBJet";
  TCut cut_J_BDT_b_closest   = "isClosestToLep && (KS_idx == -99 || ( KS_best_bdt >= -1. && KS_best_bdt < 0.0783 ) ) && isBJet && !isSJet";
  TCut cut_J_BDT_all_closest = "isClosestToLep && (KS_idx == -99 || ( KS_best_bdt >= -1. && KS_best_bdt < 0.0783 ) ) && !isSJet";

  TCut cut_J_BDT_s           = "(KS_idx == -99 || ( KS_best_bdt >= -1. && KS_best_bdt < 0.0783 ) ) && isSJet && !isBJet";
  TCut cut_J_BDT_b           = "(KS_idx == -99 || ( KS_best_bdt >= -1. && KS_best_bdt < 0.0783 ) ) && isBJet && !isSJet";
  TCut cut_J_BDT_all         = "(KS_idx == -99 || ( KS_best_bdt >= -1. && KS_best_bdt < 0.0783 ) ) && !isSJet";

  TCut cut_s_highest   = "isHighest && isSJet && !isBJet";
  TCut cut_b_highest   = "isHighest && isBJet && !isSJet";
  TCut cut_non_highest = "isHighest && !isSJet && !isBJet";
  TCut cut_all_highest = "isHighest && !isSJet";

  TCut cut_s_closest   = "isClosestToLep && isSJet && !isBJet";
  TCut cut_b_closest   = "isClosestToLep && isBJet && !isSJet";
  TCut cut_non_closest = "isClosestToLep && !isSJet && !isBJet";
  TCut cut_all_closest = "isClosestToLep && !isSJet";

  TCut cut_s           = "isSJet && !isBJet";
  TCut cut_b           = "isBJet && !isSJet";
  TCut cut_all         = "!isSJet";

  /* for 20181120 bbars_bsbar_bbbar samples ratio of train/test is about 7:3 */
  TString opt_JKS_BDT_nocondi = "nTrain_Signal=3100:nTrain_Background=78000:SplitMode=Random:NormMode=NumEvents:!V";
  TString opt_JKS_BDT_highest = "nTrain_Signal=2500:nTrain_Background=6000:SplitMode=Random:NormMode=NumEvents:!V";
  TString opt_JKS_BDT_closest = "nTrain_Signal=2200:nTrain_Background=5000:SplitMode=Random:NormMode=NumEvents:!V";

  TString opt_JKS_BDT_bnocondi = "nTrain_Signal=3100:nTrain_Background=37000:SplitMode=Random:NormMode=NumEvents:!V";
  TString opt_JKS_BDT_bhighest = "nTrain_Signal=2500:nTrain_Background=3500:SplitMode=Random:NormMode=NumEvents:!V";
  TString opt_JKS_BDT_bclosest = "nTrain_Signal=2200:nTrain_Background=3000:SplitMode=Random:NormMode=NumEvents:!V";

  TString opt_J_BDT_nocondi = "nTrain_Signal=70000:nTrain_Background=1700000:SplitMode=Random:NormMode=NumEvents:!V";
  TString opt_J_BDT_highest = "nTrain_Signal=58500:nTrain_Background=125000:SplitMode=Random:NormMode=NumEvents:!V";
  TString opt_J_BDT_closest = "nTrain_Signal=51000:nTrain_Background=100000:SplitMode=Random:NormMode=NumEvents:!V";

  TString opt_J_BDT_bnocondi = "nTrain_Signal=70000:nTrain_Background=680000:SplitMode=Random:NormMode=NumEvents:!V";
  TString opt_J_BDT_bhighest = "nTrain_Signal=58500:nTrain_Background=65000:SplitMode=Random:NormMode=NumEvents:!V";
  TString opt_J_BDT_bclosest = "nTrain_Signal=51000:nTrain_Background=58000:SplitMode=Random:NormMode=NumEvents:!V";

  setJKS(Use, factory, pp_combined_JKS_BDT_highest, bbars_bsbar_bbbar_pythia_tree, bbars_bsbar_bbbar_pythia_tree, signalWeight, backgroundWeight, cut_JKS_BDT_s_highest, cut_JKS_BDT_all_highest, opt_JKS_BDT_highest);
  setJKS(Use, factory, pp_s_vs_b_JKS_BDT_highest,   bbars_bsbar_bbbar_pythia_tree, bbars_bsbar_bbbar_pythia_tree, signalWeight, backgroundWeight, cut_JKS_BDT_s_highest, cut_JKS_BDT_b_highest,   opt_JKS_BDT_bhighest);
  setJ(  Use, factory, pp_combined_J_BDT_highest,   bbars_bsbar_bbbar_pythia_tree, bbars_bsbar_bbbar_pythia_tree, signalWeight, backgroundWeight, cut_J_BDT_s_highest,   cut_J_BDT_all_highest,   opt_J_BDT_highest);
  setJ(  Use, factory, pp_s_vs_b_J_BDT_highest,     bbars_bsbar_bbbar_pythia_tree, bbars_bsbar_bbbar_pythia_tree, signalWeight, backgroundWeight, cut_J_BDT_s_highest,   cut_J_BDT_b_highest,     opt_J_BDT_bhighest);

  setJKS(Use, factory, pp_combined_JKS_BDT_closest, bbars_bsbar_bbbar_pythia_tree, bbars_bsbar_bbbar_pythia_tree, signalWeight, backgroundWeight, cut_JKS_BDT_s_closest, cut_JKS_BDT_all_closest, opt_JKS_BDT_closest);
  setJKS(Use, factory, pp_s_vs_b_JKS_BDT_closest,   bbars_bsbar_bbbar_pythia_tree, bbars_bsbar_bbbar_pythia_tree, signalWeight, backgroundWeight, cut_JKS_BDT_s_closest, cut_JKS_BDT_b_closest,   opt_JKS_BDT_bclosest);
  setJ(  Use, factory, pp_combined_J_BDT_closest,   bbars_bsbar_bbbar_pythia_tree, bbars_bsbar_bbbar_pythia_tree, signalWeight, backgroundWeight, cut_J_BDT_s_closest,   cut_J_BDT_all_closest,   opt_J_BDT_closest);
  setJ(  Use, factory, pp_s_vs_b_J_BDT_closest,     bbars_bsbar_bbbar_pythia_tree, bbars_bsbar_bbbar_pythia_tree, signalWeight, backgroundWeight, cut_J_BDT_s_closest,   cut_J_BDT_b_closest,     opt_J_BDT_bclosest);

//  setJKS(Use, factory, pp_combined_JKS_BDT, bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight, cut_JKS_BDT_s, cut_JKS_BDT_all, opt_JKS_BDT_nocondi);
//  setJKS(Use, factory, pp_s_vs_b_JKS_BDT,   bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight, cut_JKS_BDT_s, cut_JKS_BDT_b,   opt_JKS_BDT_bnocondi);
//  setJ(  Use, factory, pp_combined_J_BDT,   bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight, cut_J_BDT_s,   cut_J_BDT_all,   opt_J_BDT_nocondi);
//  setJ(  Use, factory, pp_s_vs_b_J_BDT,     bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight, cut_J_BDT_s,   cut_J_BDT_b,     opt_J_BDT_bnocondi);

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
