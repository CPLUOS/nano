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
  dataloader->AddVariable("pt", 'F');
  dataloader->AddVariable("eta", 'F');
  dataloader->AddVariable("phi", 'F');
  dataloader->AddVariable("mass", 'F');
  dataloader->AddVariable("c_x1", 'F');
  dataloader->AddVariable("c_x2", 'F');
  dataloader->AddVariable("c_x3", 'F');
  dataloader->AddVariable("n_x1", 'F');
  dataloader->AddVariable("n_x2", 'F');
  dataloader->AddVariable("n_x3", 'F');
  dataloader->AddVariable("cmult", 'I');
  dataloader->AddVariable("nmult", 'I');
  dataloader->AddVariable("axis1", 'F');
  dataloader->AddVariable("axis2", 'F');
  dataloader->AddVariable("ptD", 'F');
  dataloader->AddVariable("area", 'F');
  dataloader->AddVariable("CSVV2", 'F'); 
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
  //

  std::map<std::string, int> Opt;

  Opt["diffBkg"]         = 0;
  Opt["diffGen"]         = 0;
  Opt["diffSam"]         = 0;
  Opt["pp"]              = 0;
  Opt["ph"]              = 0;
  Opt["hp"]              = 0;
  Opt["hp"]              = 0;


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

  TString sample_bbars_pythia = "/xrootd/store/user/wjjang/test/sum_tt012j/20180903/tt012j_bbars_2l_FxFx_sum_146.root";
  TString sample_bbars_herwig = "/xrootd/store/user/wjjang/test/sum_tt012j/20180903/tt012j_bbars_2l_FxFx_herwigpp_sum_49.root";

  TString sample_bsbar_pythia = "/xrootd/store/user/wjjang/test/sum_tt012j/20180903/tt012j_bsbar_2l_FxFx_sum_146.root";
  TString sample_bsbar_herwig = "/xrootd/store/user/wjjang/test/sum_tt012j/20180903/tt012j_bsbar_2l_FxFx_herwigpp_sum_49.root";

  TString sample_bbbar_pythia = "/xrootd/store/user/wjjang/test/sum_tt012j/20180903/tt012j_bbbar_2l_FxFx_sum_146.root";

  TString sample_bbars_bsbar_pythia = "/xrootd/store/user/wjjang/test/sum_tt012j/20180903/tt012j_bbars_bsbar_sum_146.root";
  TString sample_bbars_bsbar_herwig = "/xrootd/store/user/wjjang/test/sum_tt012j/20180903/tt012j_bbars_bsbar_herwigpp_sum_49.root";

  TString sample_bbars_bsbar_bbbar_pythia = "/xrootd/store/user/wjjang/test/sum_tt012j/20180903/tt012j_bbars_bsbar_bbbar_sum_146.root";

  sample_bbars_pythia = "/xrootd/store/user/wjjang/test/sum_tt012j/20180913/tt012j_bbars_2l_FxFx_sum_349.root";
  sample_bsbar_pythia = "/xrootd/store/user/wjjang/test/sum_tt012j/20180913/tt012j_bsbar_2l_FxFx_sum_350.root";
  sample_bbars_bsbar_pythia = "/xrootd/store/user/wjjang/test/sum_tt012j/20180913/tt012j_bbars_bsbar_sum_349_350.root";
  sample_bbars_bsbar_bbbar_pythia = "/xrootd/store/user/wjjang/test/sum_tt012j/20180913/tt012j_bbars_bsbar_bbbar_sum_349_350_201.root";

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
  TTree *bbars_herwig_tree           = (TTree*)bbars_herwig->Get("MVA_jet");

  TTree *bsbar_pythia_tree           = (TTree*)bsbar_pythia->Get("MVA_jet");
  TTree *bsbar_herwig_tree           = (TTree*)bsbar_herwig->Get("MVA_jet");

  TTree *bbbar_pythia_tree           = (TTree*)bbbar_pythia->Get("MVA_jet");

  TTree *bbars_bsbar_pythia_tree     = (TTree*)bbars_bsbar_pythia->Get("MVA_jet");
  TTree *bbars_bsbar_herwig_tree     = (TTree*)bbars_bsbar_herwig->Get("MVA_jet");

  TTree *bbars_bsbar_bbbar_pythia_tree     = (TTree*)bbars_bsbar_bbbar_pythia->Get("MVA_jet");

  // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
  TString outfileName( "vts_dR_04_Jet.root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

  TMVA::Factory *factory = new TMVA::Factory( "vts_dR_04_Jet", outputFile,
                                              "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

  TMVA::DataLoader *pp_combined_J_BDT_highest = new TMVA::DataLoader("pp_combined_J_BDT_highest");
  TMVA::DataLoader *pp_combined_J_BDT_closest = new TMVA::DataLoader("pp_combined_J_BDT_closest");
  TMVA::DataLoader *pp_combined_J_BDT         = new TMVA::DataLoader("pp_combined_J_BDT");

  TMVA::DataLoader *pp_bWsW_bWbW_highest = new TMVA::DataLoader("pp_bWsW_bWbW_highest");
  TMVA::DataLoader *pp_bWsW_bWbW_closest = new TMVA::DataLoader("pp_bWsW_bWbW_closest");
  TMVA::DataLoader *pp_bWsW_bWbW         = new TMVA::DataLoader("pp_bWsW_bWbW");

  TMVA::DataLoader *pp_J_BDT_highest = new TMVA::DataLoader("pp_J_BDT_highest");
  TMVA::DataLoader *pp_J_BDT_closest = new TMVA::DataLoader("pp_J_BDT_closest");
  TMVA::DataLoader *pp_J_BDT         = new TMVA::DataLoader("pp_J_BDT");
  TMVA::DataLoader *ph_J_BDT_highest = new TMVA::DataLoader("ph_J_BDT_highest");
  TMVA::DataLoader *ph_J_BDT_closest = new TMVA::DataLoader("ph_J_BDT_closest");
  TMVA::DataLoader *hp_J_BDT_highest = new TMVA::DataLoader("hp_J_BDT_highest");
  TMVA::DataLoader *hp_J_BDT_closest = new TMVA::DataLoader("hp_J_BDT_closest");
  TMVA::DataLoader *hh_J_BDT_highest = new TMVA::DataLoader("hh_J_BDT_highest");
  TMVA::DataLoader *hh_J_BDT_closest = new TMVA::DataLoader("hh_J_BDT_closest");

  TMVA::DataLoader *pp_s_vs_b_highest   = new TMVA::DataLoader("pp_s_vs_b_highest");
  TMVA::DataLoader *pp_s_vs_non_highest = new TMVA::DataLoader("pp_s_vs_non_highest");
  TMVA::DataLoader *pp_s_vs_all_highest = new TMVA::DataLoader("pp_s_vs_all_highest");
  TMVA::DataLoader *pp_s_vs_b_closest   = new TMVA::DataLoader("pp_s_vs_b_closest");
  TMVA::DataLoader *pp_s_vs_non_closest = new TMVA::DataLoader("pp_s_vs_non_closest");
  TMVA::DataLoader *pp_s_vs_all_closest = new TMVA::DataLoader("pp_s_vs_all_closest");
  TMVA::DataLoader *pp_s_vs_all         = new TMVA::DataLoader("pp_s_vs_all");

  TMVA::DataLoader *ph_s_vs_b_highest   = new TMVA::DataLoader("ph_s_vs_b_highest");
  TMVA::DataLoader *ph_s_vs_non_highest = new TMVA::DataLoader("ph_s_vs_non_highest");
  TMVA::DataLoader *ph_s_vs_all_highest = new TMVA::DataLoader("ph_s_vs_all_highest");
  TMVA::DataLoader *ph_s_vs_b_closest   = new TMVA::DataLoader("ph_s_vs_b_closest");
  TMVA::DataLoader *ph_s_vs_non_closest = new TMVA::DataLoader("ph_s_vs_non_closest");
  TMVA::DataLoader *ph_s_vs_all_closest = new TMVA::DataLoader("ph_s_vs_all_closest");

  TMVA::DataLoader *hp_s_vs_b_highest   = new TMVA::DataLoader("hp_s_vs_b_highest");
  TMVA::DataLoader *hp_s_vs_non_highest = new TMVA::DataLoader("hp_s_vs_non_highest");
  TMVA::DataLoader *hp_s_vs_all_highest = new TMVA::DataLoader("hp_s_vs_all_highest");
  TMVA::DataLoader *hp_s_vs_b_closest   = new TMVA::DataLoader("hp_s_vs_b_closest");
  TMVA::DataLoader *hp_s_vs_non_closest = new TMVA::DataLoader("hp_s_vs_non_closest");
  TMVA::DataLoader *hp_s_vs_all_closest = new TMVA::DataLoader("hp_s_vs_all_closest");

  TMVA::DataLoader *hh_s_vs_b_highest   = new TMVA::DataLoader("hh_s_vs_b_highest");
  TMVA::DataLoader *hh_s_vs_non_highest = new TMVA::DataLoader("hh_s_vs_non_highest");
  TMVA::DataLoader *hh_s_vs_all_highest = new TMVA::DataLoader("hh_s_vs_all_highest");
  TMVA::DataLoader *hh_s_vs_b_closest   = new TMVA::DataLoader("hh_s_vs_b_closest");
  TMVA::DataLoader *hh_s_vs_non_closest = new TMVA::DataLoader("hh_s_vs_non_closest");
  TMVA::DataLoader *hh_s_vs_all_closest = new TMVA::DataLoader("hh_s_vs_all_closest");

  //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
  //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";
  (TMVA::gConfig().GetVariablePlotting()).fMaxNumOfAllowedVariablesForScatterPlots = 99999;
  (TMVA::gConfig().GetVariablePlotting()).fNbinsXOfROCCurve = 1000;

  // global event weights per tree (see below for setting event-wise weights) 
  Double_t signalWeight     = 1.0; Double_t backgroundWeight = 1.0;

  //dataloader->SetBackgroundWeightExpression( "weight" );

  TString opt_7 = "nTrain_Signal=30000:nTrain_Background=65000:SplitMode=Random:NormMode=NumEvents:!V";

  TCut cut_J_BDT_s_highest   = "isHighest && KS_idx_pp != -99 && KS_best_bdt_pp < 0.0508 && isSJet";
  TCut cut_J_BDT_all_highest = "isHighest && KS_idx_pp != -99 && KS_best_bdt_pp < 0.0508 && !isSJet";
  TCut cut_J_BDT_s_closest   = "isClosestToLep && KS_idx_pp != -99 && KS_best_bdt_pp < 0.0508 && isSJet";
  TCut cut_J_BDT_all_closest = "isClosestToLep && KS_idx_pp != -99 && KS_best_bdt_pp < 0.0508 && !isSJet";
  TCut cut_J_BDT_s           = "KS_idx_pp != -99 && KS_best_bdt_pp < 0.0508 && isSJet";
  TCut cut_J_BDT_all         = "KS_idx_pp != -99 && KS_best_bdt_pp < 0.0508 && !isSJet";

  TCut cut_s_highest   = "isHighest && isSJet";
  TCut cut_b_highest   = "isHighest && isBJet";
  TCut cut_non_highest = "isHighest && !isSJet && !isBJet";
  TCut cut_all_highest = "isHighest && !isSJet";
  TCut cut_s_closest   = "isClosestToLep && isSJet";
  TCut cut_b_closest   = "isClosestToLep && isBJet";
  TCut cut_non_closest = "isClosestToLep && !isSJet && !isBJet";
  TCut cut_all_closest = "isClosestToLep && !isSJet";
  TCut cut_s           = "isSJet";
  TCut cut_all         = "!isSJet";

  addJetVariable(pp_combined_J_BDT);
  addTree(pp_combined_J_BDT,         bbars_bsbar_bbbar_pythia_tree, bbars_bsbar_bbbar_pythia_tree, signalWeight, backgroundWeight);
  pp_combined_J_BDT->PrepareTrainingAndTestTree(         cut_s,         cut_all,         opt_7 );
  addMethod(Use, factory, pp_combined_J_BDT);

  if (Opt["diffSam"]) {
    if (Opt["pp"]) {

      addJetVariable(pp_bWsW_bWbW_highest);
      addJetVariable(pp_bWsW_bWbW_closest);
      addJetVariable(pp_bWsW_bWbW);
      addTree(pp_bWsW_bWbW_highest, bbars_bsbar_pythia_tree, bbbar_pythia_tree, signalWeight, backgroundWeight);
      addTree(pp_bWsW_bWbW_closest, bbars_bsbar_pythia_tree, bbbar_pythia_tree, signalWeight, backgroundWeight);
      addTree(pp_bWsW_bWbW,         bbars_bsbar_pythia_tree, bbbar_pythia_tree, signalWeight, backgroundWeight);
      pp_bWsW_bWbW_highest->PrepareTrainingAndTestTree( cut_s_highest, cut_all_highest, opt_7 );
      pp_bWsW_bWbW_closest->PrepareTrainingAndTestTree( cut_s_closest, cut_all_closest, opt_7 );
      pp_bWsW_bWbW->PrepareTrainingAndTestTree(         cut_s,         cut_all,         opt_7 );
      addMethod(Use, factory, pp_bWsW_bWbW_highest);
      addMethod(Use, factory, pp_bWsW_bWbW_closest);
      addMethod(Use, factory, pp_bWsW_bWbW);
/*
      addJetVariable(pp_combined_J_BDT_highest);
      addJetVariable(pp_combined_J_BDT_closest);
      addJetVariable(pp_combined_J_BDT);
      addTree(pp_combined_J_BDT_highest, bbars_bsbar_bbbar_pythia_tree, bbars_bsbar_bbbar_pythia_tree, signalWeight, backgroundWeight);
      addTree(pp_combined_J_BDT_closest, bbars_bsbar_bbbar_pythia_tree, bbars_bsbar_bbbar_pythia_tree, signalWeight, backgroundWeight);
      addTree(pp_combined_J_BDT,         bbars_bsbar_bbbar_pythia_tree, bbars_bsbar_bbbar_pythia_tree, signalWeight, backgroundWeight);
      pp_combined_J_BDT_highest->PrepareTrainingAndTestTree( cut_s_highest, cut_all_highest, opt_7 );
      pp_combined_J_BDT_closest->PrepareTrainingAndTestTree( cut_s_closest, cut_all_closest, opt_7 );
      pp_combined_J_BDT->PrepareTrainingAndTestTree(         cut_s,         cut_all,         opt_7 );
      addMethod(Use, factory, pp_combined_J_BDT_highest);
      addMethod(Use, factory, pp_combined_J_BDT_closest);
      addMethod(Use, factory, pp_combined_J_BDT);
*/
    }
  }

  if (Opt["diffGen"]) {
    if (Opt["pp"]) {
      addJetVariable(pp_J_BDT_highest);
      addJetVariable(pp_J_BDT_closest);
      addJetVariable(pp_J_BDT);
      addTree(pp_J_BDT_highest, bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight);
      addTree(pp_J_BDT_closest, bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight);
      addTree(pp_J_BDT,         bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight);
      pp_J_BDT_highest->PrepareTrainingAndTestTree( cut_J_BDT_s_highest, cut_J_BDT_all_highest, opt_7 );
      pp_J_BDT_closest->PrepareTrainingAndTestTree( cut_J_BDT_s_closest, cut_J_BDT_all_closest, opt_7 );
      pp_J_BDT->PrepareTrainingAndTestTree(         cut_J_BDT_s,         cut_J_BDT_all,         opt_7 );
      addMethod(Use, factory, pp_J_BDT_highest);
      addMethod(Use, factory, pp_J_BDT_closest);
      addMethod(Use, factory, pp_J_BDT);
    }
    if (Opt["ph"]) {
      addJetVariable(ph_J_BDT_highest);
      addJetVariable(ph_J_BDT_closest);
      addTree(ph_J_BDT_highest, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, signalWeight, backgroundWeight);
      addTree(ph_J_BDT_closest, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, signalWeight, backgroundWeight);
      ph_J_BDT_highest->PrepareTrainingAndTestTree( cut_J_BDT_s_highest, cut_J_BDT_all_highest, opt_7 );
      ph_J_BDT_closest->PrepareTrainingAndTestTree( cut_J_BDT_s_closest, cut_J_BDT_all_closest, opt_7 );
      addMethod(Use, factory, ph_J_BDT_highest);
      addMethod(Use, factory, ph_J_BDT_closest);
    }
    if (Opt["hp"]) {
      addJetVariable(hp_J_BDT_highest);
      addJetVariable(hp_J_BDT_closest);
      addTree(hp_J_BDT_highest, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight);
      addTree(hp_J_BDT_closest, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight);
      hp_J_BDT_highest->PrepareTrainingAndTestTree( cut_J_BDT_s_highest, cut_J_BDT_all_highest, opt_7 );
      hp_J_BDT_closest->PrepareTrainingAndTestTree( cut_J_BDT_s_closest, cut_J_BDT_all_closest, opt_7 );
      addMethod(Use, factory, hp_J_BDT_highest);
      addMethod(Use, factory, hp_J_BDT_closest);
    }
    if (Opt["hh"]) {
      addJetVariable(hh_J_BDT_highest);
      addJetVariable(hh_J_BDT_closest);
      addTree(hh_J_BDT_highest, bbars_bsbar_herwig_tree, bbars_bsbar_herwig_tree, signalWeight, backgroundWeight);
      addTree(hh_J_BDT_closest, bbars_bsbar_herwig_tree, bbars_bsbar_herwig_tree, signalWeight, backgroundWeight);
      hh_J_BDT_highest->PrepareTrainingAndTestTree( cut_J_BDT_s_highest, cut_J_BDT_all_highest, opt_7 );
      hh_J_BDT_closest->PrepareTrainingAndTestTree( cut_J_BDT_s_closest, cut_J_BDT_all_closest, opt_7 );
      addMethod(Use, factory, hh_J_BDT_highest);
      addMethod(Use, factory, hh_J_BDT_closest);
    }
  }

  if (Opt["diffBkg"]) {
    if (Opt["pp"]) {
      addJetVariable(pp_s_vs_b_highest); 
      addJetVariable(pp_s_vs_non_highest); 
      addJetVariable(pp_s_vs_all_highest);
      addJetVariable(pp_s_vs_b_closest);
      addJetVariable(pp_s_vs_non_closest);
      addJetVariable(pp_s_vs_all_closest);
      addJetVariable(pp_s_vs_all);
      addTree(pp_s_vs_b_highest,   bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight);
      addTree(pp_s_vs_non_highest, bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight);
      addTree(pp_s_vs_all_highest, bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight);
      addTree(pp_s_vs_b_closest,   bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight);
      addTree(pp_s_vs_non_closest, bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight);
      addTree(pp_s_vs_all_closest, bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight);
      addTree(pp_s_vs_all,         bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight);
/*
      pp_s_vs_b_highest->PrepareTrainingAndTestTree(   cut_s_highest, cut_b_highest,   opt_7 );
      pp_s_vs_non_highest->PrepareTrainingAndTestTree( cut_s_highest, cut_non_highest, opt_7 );
      pp_s_vs_all_highest->PrepareTrainingAndTestTree( cut_s_highest, cut_all_highest, opt_7 );
      pp_s_vs_b_closest->PrepareTrainingAndTestTree(   cut_s_closest, cut_b_closest,   opt_7 );
      pp_s_vs_non_closest->PrepareTrainingAndTestTree( cut_s_closest, cut_non_closest, opt_7 );
      pp_s_vs_all_closest->PrepareTrainingAndTestTree( cut_s_closest, cut_all_closest, opt_7 );
*/
      pp_s_vs_all->PrepareTrainingAndTestTree(         cut_s,         cut_all,         opt_7 );
      addMethod(Use, factory, pp_s_vs_b_highest);
      addMethod(Use, factory, pp_s_vs_non_highest);
      addMethod(Use, factory, pp_s_vs_all_highest);
      addMethod(Use, factory, pp_s_vs_b_closest);
      addMethod(Use, factory, pp_s_vs_non_closest);
      addMethod(Use, factory, pp_s_vs_all_closest);
      addMethod(Use, factory, pp_s_vs_all);
    }
    if (Opt["ph"]) {
      addJetVariable(ph_s_vs_b_highest);
      addJetVariable(ph_s_vs_non_highest);
      addJetVariable(ph_s_vs_all_highest);
      addJetVariable(ph_s_vs_b_closest);
      addJetVariable(ph_s_vs_non_closest);
      addJetVariable(ph_s_vs_all_closest);
      addTree(ph_s_vs_b_highest,   bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, signalWeight, backgroundWeight);
      addTree(ph_s_vs_non_highest, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, signalWeight, backgroundWeight);
      addTree(ph_s_vs_all_highest, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, signalWeight, backgroundWeight);
      addTree(ph_s_vs_b_closest,   bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, signalWeight, backgroundWeight);
      addTree(ph_s_vs_non_closest, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, signalWeight, backgroundWeight);
      addTree(ph_s_vs_all_closest, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, signalWeight, backgroundWeight);
      ph_s_vs_b_highest->PrepareTrainingAndTestTree( cut_s_highest, cut_b_highest,     opt_7 );
      ph_s_vs_non_highest->PrepareTrainingAndTestTree( cut_s_highest, cut_non_highest, opt_7 );
      ph_s_vs_all_highest->PrepareTrainingAndTestTree( cut_s_highest, cut_all_highest, opt_7 );
      ph_s_vs_b_closest->PrepareTrainingAndTestTree( cut_s_closest, cut_b_closest,     opt_7 );
      ph_s_vs_non_closest->PrepareTrainingAndTestTree( cut_s_closest, cut_non_closest, opt_7 );
      ph_s_vs_all_closest->PrepareTrainingAndTestTree( cut_s_closest, cut_all_closest, opt_7 );
      addMethod(Use, factory, ph_s_vs_b_highest);
      addMethod(Use, factory, ph_s_vs_non_highest);
      addMethod(Use, factory, ph_s_vs_all_highest);
      addMethod(Use, factory, ph_s_vs_b_closest);
      addMethod(Use, factory, ph_s_vs_non_closest);
      addMethod(Use, factory, ph_s_vs_all_closest);
    }
    if (Opt["hp"]) {
      addJetVariable(hp_s_vs_b_highest);
      addJetVariable(hp_s_vs_non_highest);
      addJetVariable(hp_s_vs_all_highest);
      addJetVariable(hp_s_vs_b_closest);
      addJetVariable(hp_s_vs_non_closest);
      addJetVariable(hp_s_vs_all_closest);
      addTree(hp_s_vs_b_highest,   bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight);
      addTree(hp_s_vs_non_highest, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight);
      addTree(hp_s_vs_all_highest, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight);
      addTree(hp_s_vs_b_closest,   bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight);
      addTree(hp_s_vs_non_closest, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight);
      addTree(hp_s_vs_all_closest, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight);
      hp_s_vs_b_highest->PrepareTrainingAndTestTree( cut_s_highest, cut_b_highest,     opt_7 );
      hp_s_vs_non_highest->PrepareTrainingAndTestTree( cut_s_highest, cut_non_highest, opt_7 );
      hp_s_vs_all_highest->PrepareTrainingAndTestTree( cut_s_highest, cut_all_highest, opt_7 );
      hp_s_vs_b_closest->PrepareTrainingAndTestTree( cut_s_closest, cut_b_closest,     opt_7 );
      hp_s_vs_non_closest->PrepareTrainingAndTestTree( cut_s_closest, cut_non_closest, opt_7 );
      hp_s_vs_all_closest->PrepareTrainingAndTestTree( cut_s_closest, cut_all_closest, opt_7 );
      addMethod(Use, factory, hp_s_vs_b_highest);
      addMethod(Use, factory, hp_s_vs_non_highest);
      addMethod(Use, factory, hp_s_vs_all_highest);
      addMethod(Use, factory, hp_s_vs_b_closest);
      addMethod(Use, factory, hp_s_vs_non_closest);
      addMethod(Use, factory, hp_s_vs_all_closest);
    }
    if (Opt["hh"]) {
      addJetVariable(hh_s_vs_b_highest);
      addJetVariable(hh_s_vs_non_highest);
      addJetVariable(hh_s_vs_all_highest);
      addJetVariable(hh_s_vs_b_closest);
      addJetVariable(hh_s_vs_non_closest);
      addJetVariable(hh_s_vs_all_closest);
      addTree(hh_s_vs_b_highest,   bbars_bsbar_herwig_tree, bbars_bsbar_herwig_tree, signalWeight, backgroundWeight);
      addTree(hh_s_vs_non_highest, bbars_bsbar_herwig_tree, bbars_bsbar_herwig_tree, signalWeight, backgroundWeight);
      addTree(hh_s_vs_all_highest, bbars_bsbar_herwig_tree, bbars_bsbar_herwig_tree, signalWeight, backgroundWeight);
      addTree(hh_s_vs_b_closest,   bbars_bsbar_herwig_tree, bbars_bsbar_herwig_tree, signalWeight, backgroundWeight);
      addTree(hh_s_vs_non_closest, bbars_bsbar_herwig_tree, bbars_bsbar_herwig_tree, signalWeight, backgroundWeight);
      addTree(hh_s_vs_all_closest, bbars_bsbar_herwig_tree, bbars_bsbar_herwig_tree, signalWeight, backgroundWeight);
      hh_s_vs_b_highest->PrepareTrainingAndTestTree( cut_s_highest, cut_b_highest,     opt_7 ); 
      hh_s_vs_non_highest->PrepareTrainingAndTestTree( cut_s_highest, cut_non_highest, opt_7 ); 
      hh_s_vs_all_highest->PrepareTrainingAndTestTree( cut_s_highest, cut_all_highest, opt_7 );  
      hh_s_vs_b_closest->PrepareTrainingAndTestTree( cut_s_closest, cut_b_closest,     opt_7 );
      hh_s_vs_non_closest->PrepareTrainingAndTestTree( cut_s_closest, cut_non_closest, opt_7 );   
      hh_s_vs_all_closest->PrepareTrainingAndTestTree( cut_s_closest, cut_all_closest, opt_7 ); 
      addMethod(Use, factory, hh_s_vs_b_highest);
      addMethod(Use, factory, hh_s_vs_non_highest);
      addMethod(Use, factory, hh_s_vs_all_highest);
      addMethod(Use, factory, hh_s_vs_b_closest);
      addMethod(Use, factory, hh_s_vs_non_closest);
      addMethod(Use, factory, hh_s_vs_all_closest);
    }
  }

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
  delete pp_bWsW_bWbW_highest;
  delete pp_bWsW_bWbW_closest;
  delete pp_bWsW_bWbW;        
  delete pp_J_BDT_highest;
  delete pp_J_BDT_closest;
  delete pp_J_BDT;
  delete ph_J_BDT_highest;
  delete ph_J_BDT_closest;
  delete hp_J_BDT_highest;
  delete hp_J_BDT_closest;
  delete hh_J_BDT_highest;
  delete hh_J_BDT_closest;
  delete pp_s_vs_b_highest;
  delete pp_s_vs_non_highest;
  delete pp_s_vs_all_highest;
  delete pp_s_vs_b_closest;
  delete pp_s_vs_non_closest;
  delete pp_s_vs_all_closest;
  delete pp_s_vs_all;
  delete ph_s_vs_b_highest;
  delete ph_s_vs_non_highest;
  delete ph_s_vs_all_highest;
  delete ph_s_vs_b_closest;
  delete ph_s_vs_non_closest;
  delete ph_s_vs_all_closest;
  delete hp_s_vs_b_highest;
  delete hp_s_vs_non_highest;
  delete hp_s_vs_all_highest;
  delete hp_s_vs_b_closest;
  delete hp_s_vs_non_closest;
  delete hp_s_vs_all_closest;
  delete hh_s_vs_b_highest;
  delete hh_s_vs_non_highest;
  delete hh_s_vs_all_highest;
  delete hh_s_vs_b_closest;
  delete hh_s_vs_non_closest;
  delete hh_s_vs_all_closest;

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
