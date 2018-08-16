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

  Use["MLPBFGS"]         = 1; // Recommended ANN with optional training method
  Use["TMlpANN"]         = 1; // ROOT's own ANN
  // Boosted Decision Trees
  Use["BDT"]             = 1; // uses Adaptive Boost
  Use["BDTG"]            = 1; // uses Gradient Boost
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
  TFile *bbars_bsbar_pythia(0); TFile *bbars_bsbar_herwig(0);

//  TString sample_bsbar_pythia = "/xrootd/store/user/wjjang/test/sum_tt012j/bsbar_pythia.root";
//  TString sample_bsbar_herwig = "/xrootd/store/user/wjjang/test/sum_tt012j/bsbar_herwig.root";

  TString sample_bbars_pythia = "/xrootd/store/user/wjjang/test/sum_tt012j/20180813/tt012j_bbars_2l_FxFx_sum_146.root";
  TString sample_bbars_herwig = "/xrootd/store/user/wjjang/test/sum_tt012j/20180813/tt012j_bbars_2l_FxFx_herwigpp_sum_49.root";

  TString sample_bsbar_pythia = "/xrootd/store/user/wjjang/test/sum_tt012j/20180813/tt012j_bsbar_2l_FxFx_sum_146.root";
  TString sample_bsbar_herwig = "/xrootd/store/user/wjjang/test/sum_tt012j/20180813/tt012j_bsbar_2l_FxFx_herwigpp_sum_49.root";

  TString sample_bbars_bsbar_pythia = "/xrootd/store/user/wjjang/test/sum_tt012j/20180813/tt012j_bbars_bsbar_sum_146.root";
  TString sample_bbars_bsbar_herwig = "/xrootd/store/user/wjjang/test/sum_tt012j/20180813/tt012j_bbars_bsbar_herwigpp_sum_49.root";

  bbars_pythia = TFIle::Open( sample_bbars_pythia );             bbars_herwig = TFile::Open( sample_bbars_herwig );
  bsbar_pythia = TFile::Open( sample_bsbar_pythia );             bsbar_herwig = TFile::Open( sample_bsbar_herwig );
  bbars_bsbar_pythia = TFile::Open( sample_bbars_bsbar_pythia ); bbars_bsbar_herwig = TFile::Open( sample_bbars_bsbar_herwig );

  std::cout << "--- vts_dR_04_Jet       : Using input file 1 : " << bbars_pythia->GetName() << std::endl;
  std::cout << "--- vts_dR_04_Jet       : Using input file 2 : " << bbars_herwig->GetName() << std::endl;
  std::cout << "--- vts_dR_04_Jet       : Using input file 3 : " << bsbar_pythia->GetName() << std::endl;
  std::cout << "--- vts_dR_04_Jet       : Using input file 4 : " << bsbar_herwig->GetName() << std::endl;
  std::cout << "--- vts_dR_04_Jet       : Using input file 5 : " << bbars_bsbar_pythia->GetName() << std::endl;
  std::cout << "--- vts_dR_04_Jet       : Using input file 6 : " << bbars_bsbar_herwig->GetName() << std::endl;

  // Register the training and test trees
  TTree *bbars_pythia_tree           = (TTree*)bbars_pythia->Get("MVA_jet");
  TTree *bbars_herwig_tree           = (TTree*)bbars_herwig->Get("MVA_jet");

  TTree *bsbar_pythia_tree           = (TTree*)bsbar_pythia->Get("MVA_jet");
  TTree *bsbar_herwig_tree           = (TTree*)bsbar_herwig->Get("MVA_jet");

  TTree *bbars_bsbar_pythia_tree     = (TTree*)bbars_bsbar_pythia->Get("MVA_jet");
  TTree *bbars_bsbar_herwig_tree     = (TTree*)bbars_bsbar_herwig->Get("MVA_jet");

  // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
  TString outfileName( "vts_dR_04_Jet.root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

  TMVA::Factory *factory = new TMVA::Factory( "vts_dR_04_Jet", outputFile,
                                              "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

  TMVA::DataLoader *pp_s_vs_b_highest=   new TMVA::DataLoader("pp_s_vs_b_highest");
  TMVA::DataLoader *pp_s_vs_non_highest= new TMVA::DataLoader("pp_s_vs_non_highest");
  TMVA::DataLoader *pp_s_vs_all_highest= new TMVA::DataLoader("pp_s_vs_all_highest");
  TMVA::DataLoader *pp_s_vs_b_closest=   new TMVA::DataLoader("pp_s_vs_b_closest");
  TMVA::DataLoader *pp_s_vs_non_closest= new TMVA::DataLoader("pp_s_vs_non_closest");
  TMVA::DataLoader *pp_s_vs_all_closest= new TMVA::DataLoader("pp_s_vs_all_closest");

  TMVA::DataLoader *ph_s_vs_b_highest=   new TMVA::DataLoader("ph_s_vs_b_highest");
  TMVA::DataLoader *ph_s_vs_non_highest= new TMVA::DataLoader("ph_s_vs_non_highest");
  TMVA::DataLoader *ph_s_vs_all_highest= new TMVA::DataLoader("ph_s_vs_all_highest");
  TMVA::DataLoader *ph_s_vs_b_closest=   new TMVA::DataLoader("ph_s_vs_b_closest");
  TMVA::DataLoader *ph_s_vs_non_closest= new TMVA::DataLoader("ph_s_vs_non_closest");
  TMVA::DataLoader *ph_s_vs_all_closest= new TMVA::DataLoader("ph_s_vs_all_closest");

  TMVA::DataLoader *hp_s_vs_b_highest=   new TMVA::DataLoader("hp_s_vs_b_highest");
  TMVA::DataLoader *hp_s_vs_non_highest= new TMVA::DataLoader("hp_s_vs_non_highest");
  TMVA::DataLoader *hp_s_vs_all_highest= new TMVA::DataLoader("hp_s_vs_all_highest");
  TMVA::DataLoader *hp_s_vs_b_closest=   new TMVA::DataLoader("hp_s_vs_b_closest");
  TMVA::DataLoader *hp_s_vs_non_closest= new TMVA::DataLoader("hp_s_vs_non_closest");
  TMVA::DataLoader *hp_s_vs_all_closest= new TMVA::DataLoader("hp_s_vs_all_closest");

  TMVA::DataLoader *hh_s_vs_b_highest=   new TMVA::DataLoader("hh_s_vs_b_highest");
  TMVA::DataLoader *hh_s_vs_non_highest= new TMVA::DataLoader("hh_s_vs_non_highest");
  TMVA::DataLoader *hh_s_vs_all_highest= new TMVA::DataLoader("hh_s_vs_all_highest");
  TMVA::DataLoader *hh_s_vs_b_closest=   new TMVA::DataLoader("hh_s_vs_b_closest");
  TMVA::DataLoader *hh_s_vs_non_closest= new TMVA::DataLoader("hh_s_vs_non_closest");
  TMVA::DataLoader *hh_s_vs_all_closest= new TMVA::DataLoader("hh_s_vs_all_closest");

  //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
  //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";
  (TMVA::gConfig().GetVariablePlotting()).fMaxNumOfAllowedVariablesForScatterPlots = 99999;
  (TMVA::gConfig().GetVariablePlotting()).fNbinsXOfROCCurve = 1000;

  // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]

  addJetVariable(pp_s_vs_b_highest); 
  addJetVariable(pp_s_vs_non_highest); 
  addJetVariable(pp_s_vs_all_highest);
  addJetVariable(pp_s_vs_b_closest);
  addJetVariable(pp_s_vs_non_closest);
  addJetVariable(pp_s_vs_all_closest);

  addJetVariable(ph_s_vs_b_highest);
  addJetVariable(ph_s_vs_non_highest);
  addJetVariable(ph_s_vs_all_highest);
  addJetVariable(ph_s_vs_b_closest);
  addJetVariable(ph_s_vs_non_closest);
  addJetVariable(ph_s_vs_all_closest);

  addJetVariable(hp_s_vs_b_highest);
  addJetVariable(hp_s_vs_non_highest);
  addJetVariable(hp_s_vs_all_highest);
  addJetVariable(hp_s_vs_b_closest);
  addJetVariable(hp_s_vs_non_closest);
  addJetVariable(hp_s_vs_all_closest);

  addJetVariable(hh_s_vs_b_highest);
  addJetVariable(hh_s_vs_non_highest);
  addJetVariable(hh_s_vs_all_highest);
  addJetVariable(hh_s_vs_b_closest);
  addJetVariable(hh_s_vs_non_closest);
  addJetVariable(hh_s_vs_all_closest);

  // global event weights per tree (see below for setting event-wise weights)
  Double_t signalWeight     = 1.0; Double_t backgroundWeight = 1.0;

  // You can add an arbitrary number of signal or background trees
  addTree(pp_s_vs_b_highest,   bsbar_pythia_tree, bsbar_pythia_tree, signalWeight, backgroundWeight);
  addTree(pp_s_vs_non_highest, bsbar_pythia_tree, bsbar_pythia_tree, signalWeight, backgroundWeight);
  addTree(pp_s_vs_all_highest, bsbar_pythia_tree, bsbar_pythia_tree, signalWeight, backgroundWeight);
  addTree(pp_s_vs_b_closest,   bsbar_pythia_tree, bsbar_pythia_tree, signalWeight, backgroundWeight);
  addTree(pp_s_vs_non_closest, bsbar_pythia_tree, bsbar_pythia_tree, signalWeight, backgroundWeight);
  addTree(pp_s_vs_all_closest, bsbar_pythia_tree, bsbar_pythia_tree, signalWeight, backgroundWeight);

  addTree(ph_s_vs_b_highest,   bsbar_pythia_tree, bsbar_herwig_tree, bsbar_pythia_tree, bsbar_herwig_tree, signalWeight, backgroundWeight);
  addTree(ph_s_vs_non_highest, bsbar_pythia_tree, bsbar_herwig_tree, bsbar_pythia_tree, bsbar_herwig_tree, signalWeight, backgroundWeight);
  addTree(ph_s_vs_all_highest, bsbar_pythia_tree, bsbar_herwig_tree, bsbar_pythia_tree, bsbar_herwig_tree, signalWeight, backgroundWeight);
  addTree(ph_s_vs_b_closest,   bsbar_pythia_tree, bsbar_herwig_tree, bsbar_pythia_tree, bsbar_herwig_tree, signalWeight, backgroundWeight);
  addTree(ph_s_vs_non_closest, bsbar_pythia_tree, bsbar_herwig_tree, bsbar_pythia_tree, bsbar_herwig_tree, signalWeight, backgroundWeight);
  addTree(ph_s_vs_all_closest, bsbar_pythia_tree, bsbar_herwig_tree, bsbar_pythia_tree, bsbar_herwig_tree, signalWeight, backgroundWeight);

  addTree(hp_s_vs_b_highest,   bsbar_herwig_tree, bsbar_pythia_tree, bsbar_herwig_tree, bsbar_pythia_tree, signalWeight, backgroundWeight);
  addTree(hp_s_vs_non_highest, bsbar_herwig_tree, bsbar_pythia_tree, bsbar_herwig_tree, bsbar_pythia_tree, signalWeight, backgroundWeight);
  addTree(hp_s_vs_all_highest, bsbar_herwig_tree, bsbar_pythia_tree, bsbar_herwig_tree, bsbar_pythia_tree, signalWeight, backgroundWeight);
  addTree(hp_s_vs_b_closest,   bsbar_herwig_tree, bsbar_pythia_tree, bsbar_herwig_tree, bsbar_pythia_tree, signalWeight, backgroundWeight);
  addTree(hp_s_vs_non_closest, bsbar_herwig_tree, bsbar_pythia_tree, bsbar_herwig_tree, bsbar_pythia_tree, signalWeight, backgroundWeight);
  addTree(hp_s_vs_all_closest, bsbar_herwig_tree, bsbar_pythia_tree, bsbar_herwig_tree, bsbar_pythia_tree, signalWeight, backgroundWeight);

  addTree(hh_s_vs_b_highest,   bsbar_herwig_tree, bsbar_herwig_tree, signalWeight, backgroundWeight);
  addTree(hh_s_vs_non_highest, bsbar_herwig_tree, bsbar_herwig_tree, signalWeight, backgroundWeight);
  addTree(hh_s_vs_all_highest, bsbar_herwig_tree, bsbar_herwig_tree, signalWeight, backgroundWeight);
  addTree(hh_s_vs_b_closest,   bsbar_herwig_tree, bsbar_herwig_tree, signalWeight, backgroundWeight);
  addTree(hh_s_vs_non_closest, bsbar_herwig_tree, bsbar_herwig_tree, signalWeight, backgroundWeight);
  addTree(hh_s_vs_all_closest, bsbar_herwig_tree, bsbar_herwig_tree, signalWeight, backgroundWeight);

  //dataloader->SetBackgroundWeightExpression( "weight" );

  TCut cut_s_highest = "isSJet == 1 && isHighest == 1";
  TCut cut_b_highest = "isBJet == 1 && isHighest == 1";
  TCut cut_non_highest = "isSJet == 0 && isBJet == 0 && isHighest == 1";
  TCut cut_all_highest = "isSJet != 1 && isHighest == 1";

  TCut cut_s_closest = "isSJet == 1 && isClosestToLep == 1";
  TCut cut_b_closest = "isBJet == 1 && isClosestToLep == 1";
  TCut cut_non_closest = "isSJet == 0 && isBJet == 0 && isClosestToLep == 1";
  TCut cut_all_closest = "isSJet != 1 && isClosestToLep == 1";

  pp_s_vs_b_highest->PrepareTrainingAndTestTree( cut_s_highest, cut_b_highest,     "nTrain_Signal=3500:nTrain_Background=3500:SplitMode=Random:NormMode=NumEvents:!V" ); 
  pp_s_vs_non_highest->PrepareTrainingAndTestTree( cut_s_highest, cut_non_highest, "nTrain_Signal=3500:nTrain_Background=3500:SplitMode=Random:NormMode=NumEvents:!V" ); 
  pp_s_vs_all_highest->PrepareTrainingAndTestTree( cut_s_highest, cut_all_highest, "nTrain_Signal=3500:nTrain_Background=3500:SplitMode=Random:NormMode=NumEvents:!V" ); 
  pp_s_vs_b_closest->PrepareTrainingAndTestTree( cut_s_closest, cut_b_closest,     "nTrain_Signal=3500:nTrain_Background=3500:SplitMode=Random:NormMode=NumEvents:!V" );
  pp_s_vs_non_closest->PrepareTrainingAndTestTree( cut_s_closest, cut_non_closest, "nTrain_Signal=3500:nTrain_Background=3500:SplitMode=Random:NormMode=NumEvents:!V" );  
  pp_s_vs_all_closest->PrepareTrainingAndTestTree( cut_s_closest, cut_all_closest, "nTrain_Signal=3500:nTrain_Background=3500:SplitMode=Random:NormMode=NumEvents:!V" );  

  ph_s_vs_b_highest->PrepareTrainingAndTestTree( cut_s_highest, cut_b_highest,     "nTrain_Signal=3500:nTrain_Background=3500:SplitMode=Random:NormMode=NumEvents:!V" );
  ph_s_vs_non_highest->PrepareTrainingAndTestTree( cut_s_highest, cut_non_highest, "nTrain_Signal=3500:nTrain_Background=3500:SplitMode=Random:NormMode=NumEvents:!V" );
  ph_s_vs_all_highest->PrepareTrainingAndTestTree( cut_s_highest, cut_all_highest, "nTrain_Signal=3500:nTrain_Background=3500:SplitMode=Random:NormMode=NumEvents:!V" );
  ph_s_vs_b_closest->PrepareTrainingAndTestTree( cut_s_closest, cut_b_closest,     "nTrain_Signal=3500:nTrain_Background=3500:SplitMode=Random:NormMode=NumEvents:!V" );
  ph_s_vs_non_closest->PrepareTrainingAndTestTree( cut_s_closest, cut_non_closest, "nTrain_Signal=3500:nTrain_Background=3500:SplitMode=Random:NormMode=NumEvents:!V" );
  ph_s_vs_all_closest->PrepareTrainingAndTestTree( cut_s_closest, cut_all_closest, "nTrain_Signal=3500:nTrain_Background=3500:SplitMode=Random:NormMode=NumEvents:!V" );

  hp_s_vs_b_highest->PrepareTrainingAndTestTree( cut_s_highest, cut_b_highest,     "nTrain_Signal=3500:nTrain_Background=3500:SplitMode=Random:NormMode=NumEvents:!V" ); 
  hp_s_vs_non_highest->PrepareTrainingAndTestTree( cut_s_highest, cut_non_highest, "nTrain_Signal=3500:nTrain_Background=3500:SplitMode=Random:NormMode=NumEvents:!V" );  
  hp_s_vs_all_highest->PrepareTrainingAndTestTree( cut_s_highest, cut_all_highest, "nTrain_Signal=3500:nTrain_Background=3500:SplitMode=Random:NormMode=NumEvents:!V" ); 
  hp_s_vs_b_closest->PrepareTrainingAndTestTree( cut_s_closest, cut_b_closest,     "nTrain_Signal=3500:nTrain_Background=3500:SplitMode=Random:NormMode=NumEvents:!V" );
  hp_s_vs_non_closest->PrepareTrainingAndTestTree( cut_s_closest, cut_non_closest, "nTrain_Signal=3500:nTrain_Background=3500:SplitMode=Random:NormMode=NumEvents:!V" ); 
  hp_s_vs_all_closest->PrepareTrainingAndTestTree( cut_s_closest, cut_all_closest, "nTrain_Signal=3500:nTrain_Background=3500:SplitMode=Random:NormMode=NumEvents:!V" ); 

  hh_s_vs_b_highest->PrepareTrainingAndTestTree( cut_s_highest, cut_b_highest,     "nTrain_Signal=3500:nTrain_Background=3500:SplitMode=Random:NormMode=NumEvents:!V" ); 
  hh_s_vs_non_highest->PrepareTrainingAndTestTree( cut_s_highest, cut_non_highest, "nTrain_Signal=3500:nTrain_Background=3500:SplitMode=Random:NormMode=NumEvents:!V" ); 
  hh_s_vs_all_highest->PrepareTrainingAndTestTree( cut_s_highest, cut_all_highest, "nTrain_Signal=3500:nTrain_Background=3500:SplitMode=Random:NormMode=NumEvents:!V" );  
  hh_s_vs_b_closest->PrepareTrainingAndTestTree( cut_s_closest, cut_b_closest,     "nTrain_Signal=3500:nTrain_Background=3500:SplitMode=Random:NormMode=NumEvents:!V" );
  hh_s_vs_non_closest->PrepareTrainingAndTestTree( cut_s_closest, cut_non_closest, "nTrain_Signal=3500:nTrain_Background=3500:SplitMode=Random:NormMode=NumEvents:!V" );   
  hh_s_vs_all_closest->PrepareTrainingAndTestTree( cut_s_closest, cut_all_closest, "nTrain_Signal=3500:nTrain_Background=3500:SplitMode=Random:NormMode=NumEvents:!V" );  

  addMethod(Use, factory, pp_s_vs_b_highest);
  addMethod(Use, factory, pp_s_vs_non_highest);
  addMethod(Use, factory, pp_s_vs_all_highest);
  addMethod(Use, factory, pp_s_vs_b_closest);
  addMethod(Use, factory, pp_s_vs_non_closest);
  addMethod(Use, factory, pp_s_vs_all_closest);

  addMethod(Use, factory, ph_s_vs_b_highest);
  addMethod(Use, factory, ph_s_vs_non_highest);
  addMethod(Use, factory, ph_s_vs_all_highest);
  addMethod(Use, factory, ph_s_vs_b_closest);
  addMethod(Use, factory, ph_s_vs_non_closest);
  addMethod(Use, factory, ph_s_vs_all_closest);

  addMethod(Use, factory, hp_s_vs_b_highest);
  addMethod(Use, factory, hp_s_vs_non_highest);
  addMethod(Use, factory, hp_s_vs_all_highest);
  addMethod(Use, factory, hp_s_vs_b_closest);
  addMethod(Use, factory, hp_s_vs_non_closest);
  addMethod(Use, factory, hp_s_vs_all_closest);

  addMethod(Use, factory, hh_s_vs_b_highest);
  addMethod(Use, factory, hh_s_vs_non_highest);
  addMethod(Use, factory, hh_s_vs_all_highest);
  addMethod(Use, factory, hh_s_vs_b_closest);
  addMethod(Use, factory, hh_s_vs_non_closest);
  addMethod(Use, factory, hh_s_vs_all_closest);

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
  delete pp_s_vs_b_highest;
  delete pp_s_vs_non_highest;
  delete pp_s_vs_all_highest;
  delete pp_s_vs_b_closest;
  delete pp_s_vs_non_closest;
  delete pp_s_vs_all_closest;
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
