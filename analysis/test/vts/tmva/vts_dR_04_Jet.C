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
  dataloader->AddVariable("KS_d_pp",                'F');
  dataloader->AddVariable("KS_pt_pp",               'F');
  dataloader->AddVariable("KS_eta_pp",              'F');
  dataloader->AddVariable("KS_phi_pp",              'F');
  dataloader->AddVariable("KS_mass_pp",             'F');
  dataloader->AddVariable("KS_lxy_pp",              'F');
  dataloader->AddVariable("KS_lxySig_pp",           'F');
  dataloader->AddVariable("KS_l3D_pp",              'F');
  dataloader->AddVariable("KS_l3DSig_pp",           'F');
  dataloader->AddVariable("KS_legDR_pp",            'F');
  dataloader->AddVariable("KS_angleXY_pp",          'F');
  dataloader->AddVariable("KS_angleXYZ_pp",         'F');
  dataloader->AddVariable("KS_chi2_pp",             'F');
  dataloader->AddVariable("KS_dca_pp",              'F');
  dataloader->AddVariable("KS_dau1_chi2_pp",        'F');
  dataloader->AddVariable("KS_dau1_ipsigXY_pp",     'F');
  dataloader->AddVariable("KS_dau1_ipsigZ_pp",      'F');
  dataloader->AddVariable("KS_dau1_pt_pp",          'F');
  dataloader->AddVariable("KS_dau2_chi2_pp",        'F');
  dataloader->AddVariable("KS_dau2_ipsigXY_pp",     'F');
  dataloader->AddVariable("KS_dau2_ipsigZ_pp",      'F');
  dataloader->AddVariable("KS_dau2_pt_pp",          'F');
  dataloader->AddVariable("KS_best_bdt_pp",         'F');
  dataloader->AddVariable("KS_x_pp", 'F');
}

void addHadVariablePH(TMVA::DataLoader *dataloader) {
  dataloader->AddVariable("KS_d_ph",                'F');
  dataloader->AddVariable("KS_pt_ph",               'F');
  dataloader->AddVariable("KS_eta_ph",              'F');
  dataloader->AddVariable("KS_phi_ph",              'F');
  dataloader->AddVariable("KS_mass_ph",             'F');
  dataloader->AddVariable("KS_lxy_ph",              'F');
  dataloader->AddVariable("KS_lxySig_ph",           'F');
  dataloader->AddVariable("KS_l3D_ph",              'F');
  dataloader->AddVariable("KS_l3DSig_ph",           'F');
  dataloader->AddVariable("KS_legDR_ph",            'F');
  dataloader->AddVariable("KS_angleXY_ph",          'F');
  dataloader->AddVariable("KS_angleXYZ_ph",         'F');
  dataloader->AddVariable("KS_chi2_ph",             'F');
  dataloader->AddVariable("KS_dca_ph",              'F');
  dataloader->AddVariable("KS_dau1_chi2_ph",        'F');
  dataloader->AddVariable("KS_dau1_ipsigXY_ph",     'F');
  dataloader->AddVariable("KS_dau1_ipsigZ_ph",      'F');
  dataloader->AddVariable("KS_dau1_pt_ph",          'F');
  dataloader->AddVariable("KS_dau2_chi2_ph",        'F');
  dataloader->AddVariable("KS_dau2_ipsigXY_ph",     'F');
  dataloader->AddVariable("KS_dau2_ipsigZ_ph",      'F');
  dataloader->AddVariable("KS_dau2_pt_ph",          'F');
  dataloader->AddVariable("KS_best_bdt_ph",         'F');
  dataloader->AddVariable("KS_x_ph", 'F');
}

void addHadVariableHP(TMVA::DataLoader *dataloader) {
  dataloader->AddVariable("KS_d_hp",                'F');
  dataloader->AddVariable("KS_pt_hp",               'F');
  dataloader->AddVariable("KS_eta_hp",              'F');
  dataloader->AddVariable("KS_phi_hp",              'F');
  dataloader->AddVariable("KS_mass_hp",             'F');
  dataloader->AddVariable("KS_lxy_hp",              'F');
  dataloader->AddVariable("KS_lxySig_hp",           'F');
  dataloader->AddVariable("KS_l3D_hp",              'F');
  dataloader->AddVariable("KS_l3DSig_hp",           'F');
  dataloader->AddVariable("KS_legDR_hp",            'F');
  dataloader->AddVariable("KS_angleXY_hp",          'F');
  dataloader->AddVariable("KS_angleXYZ_hp",         'F');
  dataloader->AddVariable("KS_chi2_hp",             'F');
  dataloader->AddVariable("KS_dca_hp",              'F');
  dataloader->AddVariable("KS_dau1_chi2_hp",        'F');
  dataloader->AddVariable("KS_dau1_ipsigXY_hp",     'F');
  dataloader->AddVariable("KS_dau1_ipsigZ_hp",      'F');
  dataloader->AddVariable("KS_dau1_pt_hp",          'F');
  dataloader->AddVariable("KS_dau2_chi2_hp",        'F');
  dataloader->AddVariable("KS_dau2_ipsigXY_hp",     'F');
  dataloader->AddVariable("KS_dau2_ipsigZ_hp",      'F');
  dataloader->AddVariable("KS_dau2_pt_hp",          'F');
  dataloader->AddVariable("KS_best_bdt_hp",         'F');
  dataloader->AddVariable("KS_x_hp", 'F');
}

void addHadVariableHH(TMVA::DataLoader *dataloader) {
  dataloader->AddVariable("KS_d_hh",                'F');
  dataloader->AddVariable("KS_pt_hh",               'F');
  dataloader->AddVariable("KS_eta_hh",              'F');
  dataloader->AddVariable("KS_phi_hh",              'F');
  dataloader->AddVariable("KS_mass_hh",             'F');
  dataloader->AddVariable("KS_lxy_hh",              'F');
  dataloader->AddVariable("KS_lxySig_hh",           'F');
  dataloader->AddVariable("KS_l3D_hh",              'F');
  dataloader->AddVariable("KS_l3DSig_hh",           'F');
  dataloader->AddVariable("KS_legDR_hh",            'F');
  dataloader->AddVariable("KS_angleXY_hh",          'F');
  dataloader->AddVariable("KS_angleXYZ_hh",         'F');
  dataloader->AddVariable("KS_chi2_hh",             'F');
  dataloader->AddVariable("KS_dca_hh",              'F');
  dataloader->AddVariable("KS_dau1_chi2_hh",        'F');
  dataloader->AddVariable("KS_dau1_ipsigXY_hh",     'F');
  dataloader->AddVariable("KS_dau1_ipsigZ_hh",      'F');
  dataloader->AddVariable("KS_dau1_pt_hh",          'F');
  dataloader->AddVariable("KS_dau2_chi2_hh",        'F');
  dataloader->AddVariable("KS_dau2_ipsigXY_hh",     'F');
  dataloader->AddVariable("KS_dau2_ipsigZ_hh",      'F');
  dataloader->AddVariable("KS_dau2_pt_hh",          'F');
  dataloader->AddVariable("KS_best_bdt_hh",         'F');
  dataloader->AddVariable("KS_x_hh", 'F');
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
  if (string(dataloader->GetName()).find("ph_") != std::string::npos) addHadVariablePH(dataloader);
  if (string(dataloader->GetName()).find("hp_") != std::string::npos) addHadVariableHP(dataloader);
  if (string(dataloader->GetName()).find("hh_") != std::string::npos) addHadVariableHH(dataloader);
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
  if (string(dataloader->GetName()).find("ph_") != std::string::npos) addHadVariablePH(dataloader);
  if (string(dataloader->GetName()).find("hp_") != std::string::npos) addHadVariableHP(dataloader);
  if (string(dataloader->GetName()).find("hh_") != std::string::npos) addHadVariableHH(dataloader);
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
  Use["BDTG"]            = 0; // uses Gradient Boost

  std::map<std::string, int> Opt;

  Opt["diffBkg"]         = 0;
  Opt["diffGen"]         = 0;
  Opt["diffSam"]         = 1;
  Opt["pp"]              = 1;
  Opt["ph"]              = 0;
  Opt["hp"]              = 0;
  Opt["hp"]              = 0;

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

  TString sample_bbars_pythia = "/xrootd/store/user/wjjang/test/sum_tt012j/20181019/tt012j_bbars_2l_FxFx_sum_349.root";
  TString sample_bsbar_pythia = "/xrootd/store/user/wjjang/test/sum_tt012j/20181019/tt012j_bsbar_2l_FxFx_sum_774.root";
  TString sample_bbbar_pythia = "/xrootd/store/user/wjjang/test/sum_tt012j/20181019/tt012j_bbbar_2l_FxFx_sum_201.root";
  TString sample_bbars_bsbar_pythia = "/xrootd/store/user/wjjang/test/sum_tt012j/20181019/tt012j_bbars_bsbar_sum_349_774.root";
  TString sample_bbars_bsbar_bbbar_pythia = "/xrootd/store/user/wjjang/test/sum_tt012j/20181019/tt012j_bbars_bsbar_bbbar_sum_349_774_201.root";

  TString sample_bbars_herwig = "/xrootd/store/user/wjjang/test/sum_tt012j/20181019/tt012j_bbars_2l_FxFx_herwigpp_sum_49.root";
  TString sample_bsbar_herwig = "/xrootd/store/user/wjjang/test/sum_tt012j/20181019/tt012j_bsbar_2l_FxFx_herwigpp_sum_49.root";
  TString sample_bbars_bsbar_herwig = "/xrootd/store/user/wjjang/test/sum_tt012j/20181019/tt012j_bbars_bsbar_herwigpp_sum_49.root";

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
//  TString outfileName( "./output/vts_dR_04_Jet.root" );
  TString outfileName( "vts_dR_04_Jet_20181022.root" );   
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
   
  // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]

  //dataloader->SetBackgroundWeightExpression( "weight" );

  /* Cuts for JKS >>>> KS_best_bdt_pp >= 0.0726 : jet includes Real KS */
  TCut cut_JKS_BDT_s_highest   = "isHighest && KS_idx_pp != -99 && KS_best_bdt_pp >= 0.0726 && isSJet && !isBJet";
  TCut cut_JKS_BDT_b_highest = "isHighest && KS_idx_pp != -99 && KS_best_bdt_pp >= 0.0726 && isBJet && !isSJet";

  TCut cut_JKS_BDT_all_highest = "isHighest && KS_idx_pp != -99 && KS_best_bdt_pp >= 0.0726 && !isSJet";
  TCut cut_JKS_BDT_s_closest   = "isClosestToLep && KS_idx_pp != -99 && KS_best_bdt_pp >= 0.0726 && isSJet";
  TCut cut_JKS_BDT_b_closest = "isHighest && KS_idx_pp != -99 && KS_best_bdt_pp >= 0.0726 && isBJet";
  TCut cut_JKS_BDT_all_closest = "isClosestToLep && KS_idx_pp != -99 && KS_best_bdt_pp >= 0.0726 && !isSJet";
  TCut cut_JKS_BDT_s           = "KS_idx_pp != -99 && KS_best_bdt_pp >= 0.0726 && isSJet";
  TCut cut_JKS_BDT_b           = "KS_idx_pp != -99 && KS_best_bdt_pp >= 0.0726 && isBJet";
  TCut cut_JKS_BDT_all         = "KS_idx_pp != -99 && KS_best_bdt_pp >= 0.0726 && !isSJet";

  /* Cuts for J >>>> KS_best_bdt_pp < 0.0726 : jet includes Fake KS | KS_idx_pp == -99 : KS isn't in the jet or there is no KS in a event */ 
  TCut cut_J_BDT_s_highest   = "isHighest && (KS_idx_pp == -99 || ( KS_best_bdt_pp >= -1. && KS_best_bdt_pp < 0.0726 ) ) && isSJet";
  TCut cut_J_BDT_b_highest   = "isHighest && (KS_idx_pp == -99 || ( KS_best_bdt_pp >= -1. && KS_best_bdt_pp < 0.0726 ) ) && isBJet";
  TCut cut_J_BDT_all_highest = "isHighest && (KS_idx_pp == -99 || ( KS_best_bdt_pp >= -1. && KS_best_bdt_pp < 0.0726 ) ) && !isSJet";
  TCut cut_J_BDT_s_closest   = "isClosestToLep && (KS_idx_pp == -99 || ( KS_best_bdt_pp >= -1. && KS_best_bdt_pp < 0.0726 ) ) && isSJet";
  TCut cut_J_BDT_b_closest   = "isClosestToLep && (KS_idx_pp == -99 || ( KS_best_bdt_pp >= -1. && KS_best_bdt_pp < 0.0726 ) ) && isBJet";
  TCut cut_J_BDT_all_closest = "isClosestToLep && (KS_idx_pp == -99 || ( KS_best_bdt_pp >= -1. && KS_best_bdt_pp < 0.0726 ) ) && !isSJet";
  TCut cut_J_BDT_s           = "(KS_idx_pp == -99 || ( KS_best_bdt_pp >= -1. && KS_best_bdt_pp < 0.0726 ) ) && isSJet";
  TCut cut_J_BDT_b           = "(KS_idx_pp == -99 || ( KS_best_bdt_pp >= -1. && KS_best_bdt_pp < 0.0726 ) ) && isBJet";
  TCut cut_J_BDT_all         = "(KS_idx_pp == -99 || ( KS_best_bdt_pp >= -1. && KS_best_bdt_pp < 0.0726 ) ) && !isSJet";

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
/*
  TCut cut_s_highest   = "isHighest && KS_idx_pp != -99 && KS_nMatched_pp == 2 && isSJet";
  TCut cut_b_highest   = "isHighest && KS_idx_pp != -99 && KS_nMatched_pp == 2 && isBJet";
  TCut cut_non_highest = "isHighest && KS_idx_pp != -99 && KS_nMatched_pp == 2 && !isSJet && !isBJet";
  TCut cut_all_highest = "isHighest && KS_idx_pp != -99 && KS_nMatched_pp == 2 && !isSJet";
  TCut cut_s_closest   = "isClosestToLep && KS_idx_pp != -99 && KS_nMatched_pp == 2 && isSJet";
  TCut cut_b_closest   = "isClosestToLep && KS_idx_pp != -99 && KS_nMatched_pp == 2 && isBJet";
  TCut cut_non_closest = "isClosestToLep && KS_idx_pp != -99 && KS_nMatched_pp == 2 && !isSJet && !isBJet";
  TCut cut_all_closest = "isClosestToLep && KS_idx_pp != -99 && KS_nMatched_pp == 2 && !isSJet";
  TCut cut_s           = "KS_idx_pp != -99 && KS_nMatched_pp == 2 && isSJet";
  TCut cut_all         = "KS_idx_pp != -99 && KS_nMatched_pp == 2 && !isSJet";
*/

  /* for 20181019 samples ratio of train/test is about 7:3 */
  TString opt_JKS_BDT_nocondi = "nTrain_Signal=2400:nTrain_Background=5800:SplitMode=Random:NormMode=NumEvents:!V";
  TString opt_JKS_BDT_highest = "nTrain_Signal=2000:nTrain_Background=4000:SplitMode=Random:NormMode=NumEvents:!V";
  TString opt_JKS_BDT_closest = "nTrain_Signal=1700:nTrain_Background=3300:SplitMode=Random:NormMode=NumEvents:!V";

  TString opt_JKS_BDT_bnocondi = "nTrain_Signal=2400:nTrain_Background=580000:SplitMode=Random:NormMode=NumEvents:!V";
  TString opt_JKS_BDT_bhighest = "nTrain_Signal=2000:nTrain_Background=2100:SplitMode=Random:NormMode=NumEvents:!V";
  TString opt_JKS_BDT_bclosest = "nTrain_Signal=1700:nTrain_Background=2100:SplitMode=Random:NormMode=NumEvents:!V";

/*
  TString opt_JKS_nocondi     = "nTrain_Signal=5700000:nTrain_Background=6500000:SplitMode=Random:NormMode=NumEvents:!V";
  TString opt_JKS_highest     = "nTrain_Signal=2300000:nTrain_Background=4200000:SplitMode=Random:NormMode=NumEvents:!V";
  TString opt_JKS_closest     = "nTrain_Signal=2000000:nTrain_Background=3500000:SplitMode=Random:NormMode=NumEvents:!V";
*/
  TString opt_J_BDT_nocondi = "nTrain_Signal=57000:nTrain_Background=135000:SplitMode=Random:NormMode=NumEvents:!V";
  TString opt_J_BDT_highest = "nTrain_Signal=47000:nTrain_Background=87000:SplitMode=Random:NormMode=NumEvents:!V";
  TString opt_J_BDT_closest = "nTrain_Signal=41000:nTrain_Background=71000:SplitMode=Random:NormMode=NumEvents:!V";

  TString opt_J_BDT_bnocondi = "nTrain_Signal=57000:nTrain_Background=13500000:SplitMode=Random:NormMode=NumEvents:!V";
  TString opt_J_BDT_bhighest = "nTrain_Signal=47000:nTrain_Background=40000:SplitMode=Random:NormMode=NumEvents:!V";
  TString opt_J_BDT_bclosest = "nTrain_Signal=41000:nTrain_Background=35000:SplitMode=Random:NormMode=NumEvents:!V";
/*
  TString opt_J_nocondi     = "nTrain_Signal=3000000:nTrain_Background=6500000:SplitMode=Random:NormMode=NumEvents:!V";
  TString opt_J_highest     = "nTrain_Signal=2500000:nTrain_Background=4200000:SplitMode=Random:NormMode=NumEvents:!V";
  TString opt_J_closest     = "nTrain_Signal=2000000:nTrain_Background=3500000:SplitMode=Random:NormMode=NumEvents:!V";
*/
  if (Opt["diffSam"]) {
    if (Opt["pp"]) {
//      setJKS(Use, factory, pp_combined_JKS_BDT_highest, bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight, cut_JKS_BDT_s_highest, cut_JKS_BDT_all_highest, opt_JKS_BDT_highest);
//      setJKS(Use, factory, pp_combined_JKS_BDT_closest, bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight, cut_JKS_BDT_s_closest, cut_JKS_BDT_all_closest, opt_JKS_BDT_closest);

//      setJ(Use, factory, pp_combined_J_BDT_highest, bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight, cut_J_BDT_s_highest, cut_J_BDT_all_highest, opt_J_BDT_highest);
//      setJ(Use, factory, pp_combined_J_BDT_closest, bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight, cut_J_BDT_s_closest, cut_J_BDT_all_closest, opt_J_BDT_closest);

      setJKS(Use, factory, pp_s_vs_b_JKS_BDT_highest, bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight, cut_JKS_BDT_s_highest, cut_JKS_BDT_b_highest, opt_JKS_BDT_bhighest); 
//      setJKS(Use, factory, pp_s_vs_b_JKS_BDT_closest, bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight, cut_JKS_BDT_s_closest, cut_JKS_BDT_b_closest, opt_JKS_BDT_bclosest);

//      setJ(Use, factory, pp_s_vs_b_J_BDT_highest, bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight, cut_J_BDT_s_highest, cut_J_BDT_b_highest, opt_J_BDT_bhighest);
//      setJ(Use, factory, pp_s_vs_b_J_BDT_closest, bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight, cut_J_BDT_s_closest, cut_J_BDT_b_closest, opt_J_BDT_bclosest);

/*
      setJKS(Use, factory, pp_combined_JKS_BDT        , bbars_bsbar_bbbar_pythia_tree, bbars_bsbar_bbbar_pythia_tree, signalWeight, backgroundWeight, cut_JKS_BDT_s        , cut_JKS_BDT_all        , opt_JKS_BDT_nocondi);
      setJKS(Use, factory, pp_combined_JKS_BDT_highest, bbars_bsbar_bbbar_pythia_tree, bbars_bsbar_bbbar_pythia_tree, signalWeight, backgroundWeight, cut_JKS_BDT_s_highest, cut_JKS_BDT_all_highest, opt_JKS_BDT_highest);
      setJKS(Use, factory, pp_combined_JKS_BDT_closest, bbars_bsbar_bbbar_pythia_tree, bbars_bsbar_bbbar_pythia_tree, signalWeight, backgroundWeight, cut_JKS_BDT_s_closest, cut_JKS_BDT_all_closest, opt_JKS_BDT_closest);

      setJ(Use, factory, pp_combined_J_BDT        , bbars_bsbar_bbbar_pythia_tree, bbars_bsbar_bbbar_pythia_tree, signalWeight, backgroundWeight, cut_J_BDT_s        , cut_J_BDT_all        , opt_J_BDT_nocondi);
      setJ(Use, factory, pp_combined_J_BDT_highest, bbars_bsbar_bbbar_pythia_tree, bbars_bsbar_bbbar_pythia_tree, signalWeight, backgroundWeight, cut_J_BDT_s_highest, cut_J_BDT_all_highest, opt_J_BDT_highest);
      setJ(Use, factory, pp_combined_J_BDT_closest, bbars_bsbar_bbbar_pythia_tree, bbars_bsbar_bbbar_pythia_tree, signalWeight, backgroundWeight, cut_J_BDT_s_closest, cut_J_BDT_all_closest, opt_J_BDT_closest);
*/
    }
  }
  if (Opt["diffGen"]) {
    if (Opt["pp"]) {
    }
    if (Opt["ph"]) {
    }
    if (Opt["hp"]) {
    }
    if (Opt["hh"]) {
    }
  }
  if (Opt["diffBkg"]) {
    if (Opt["pp"]) {
/*
      setJKS(Use, factory, pp_s_vs_all_with_KS        , bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight, cut_s        , cut_all        , opt_JKS_nocondi);
      setJKS(Use, factory, pp_s_vs_all_with_KS_highest, bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight, cut_s_highest, cut_all_highest, opt_JKS_highest);
      setJKS(Use, factory, pp_s_vs_all_with_KS_closest, bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight, cut_s_closest, cut_all_closest, opt_JKS_closest);

      setJ(Use, factory, pp_s_vs_all        , bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight, cut_s        , cut_all        , opt_J_nocondi);
      setJ(Use, factory, pp_s_vs_all_highest, bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight, cut_s_highest, cut_all_highest, opt_J_highest);
      setJ(Use, factory, pp_s_vs_all_closest, bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight, cut_s_closest, cut_all_closest, opt_J_closest);
*/
    }
    if (Opt["ph"]) {
    }
    if (Opt["hp"]) {
    }
    if (Opt["hh"]) {
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
