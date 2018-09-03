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
  dataloader->AddVariable("KS_x_pp := KS_pt_pp/pt", 'F');
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
  dataloader->AddVariable("KS_x_ph := KS_pt_ph/pt", 'F');
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
  dataloader->AddVariable("KS_x_hp := KS_pt_hp/pt", 'F');
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
  dataloader->AddVariable("KS_x_hh := KS_pt_hh/pt", 'F');
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

int vts_dR_04_Jet_With_KS( TString myMethodList = "" )
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
  // ---------------------------------------------------------------

  std::cout << std::endl;
  std::cout << "==> Start vts_dR_04_Jet_With_KS" << std::endl;

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

  TString sample_bbars_pythia = "/xrootd/store/user/wjjang/test/sum_tt012j/20180903/tt012j_bbars_2l_FxFx_sum_146.root";
  TString sample_bbars_herwig = "/xrootd/store/user/wjjang/test/sum_tt012j/20180903/tt012j_bbars_2l_FxFx_herwigpp_sum_49.root";

  TString sample_bsbar_pythia = "/xrootd/store/user/wjjang/test/sum_tt012j/20180903/tt012j_bsbar_2l_FxFx_sum_146.root";
  TString sample_bsbar_herwig = "/xrootd/store/user/wjjang/test/sum_tt012j/20180903/tt012j_bsbar_2l_FxFx_herwigpp_sum_49.root";

  TString sample_bbars_bsbar_pythia = "/xrootd/store/user/wjjang/test/sum_tt012j/20180903/tt012j_bbars_bsbar_sum_146.root";
  TString sample_bbars_bsbar_herwig = "/xrootd/store/user/wjjang/test/sum_tt012j/20180903/tt012j_bbars_bsbar_herwigpp_sum_49.root";

  bbars_pythia = TFile::Open( sample_bbars_pythia );             bbars_herwig = TFile::Open( sample_bbars_herwig );
  bsbar_pythia = TFile::Open( sample_bsbar_pythia );             bsbar_herwig = TFile::Open( sample_bsbar_herwig );
  bbars_bsbar_pythia = TFile::Open( sample_bbars_bsbar_pythia ); bbars_bsbar_herwig = TFile::Open( sample_bbars_bsbar_herwig );

  std::cout << "--- vts_dR_04_Jet_With_KS       : Using input file 1 : " << bbars_pythia->GetName() << std::endl;
  std::cout << "--- vts_dR_04_Jet_With_KS       : Using input file 2 : " << bbars_herwig->GetName() << std::endl;
  std::cout << "--- vts_dR_04_Jet_With_KS       : Using input file 3 : " << bsbar_pythia->GetName() << std::endl;
  std::cout << "--- vts_dR_04_Jet_With_KS       : Using input file 4 : " << bsbar_herwig->GetName() << std::endl;
  std::cout << "--- vts_dR_04_Jet_With_KS       : Using input file 5 : " << bbars_bsbar_pythia->GetName() << std::endl;
  std::cout << "--- vts_dR_04_Jet_With_KS       : Using input file 6 : " << bbars_bsbar_herwig->GetName() << std::endl;

  // Register the training and test trees
  TTree *bbars_pythia_tree           = (TTree*)bbars_pythia->Get("MVA_jet");
  TTree *bbars_herwig_tree           = (TTree*)bbars_herwig->Get("MVA_jet");

  TTree *bsbar_pythia_tree           = (TTree*)bsbar_pythia->Get("MVA_jet");
  TTree *bsbar_herwig_tree           = (TTree*)bsbar_herwig->Get("MVA_jet");

  TTree *bbars_bsbar_pythia_tree     = (TTree*)bbars_bsbar_pythia->Get("MVA_jet");
  TTree *bbars_bsbar_herwig_tree     = (TTree*)bbars_bsbar_herwig->Get("MVA_jet");

  // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
  TString outfileName( "vts_dR_04_Jet_With_KS.root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

  TMVA::Factory *factory = new TMVA::Factory( "vts_dR_04_Jet_With_KS", outputFile,
//                                              "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
                                              "!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification" );

  TMVA::DataLoader *pp_JKS_BDT_highest = new TMVA::DataLoader("pp_JKS_BDT_highest");
  TMVA::DataLoader *pp_JKS_BDT_closest = new TMVA::DataLoader("pp_JKS_BDT_closest");
  TMVA::DataLoader *ph_JKS_BDT_highest = new TMVA::DataLoader("ph_JKS_BDT_highest");
  TMVA::DataLoader *ph_JKS_BDT_closest = new TMVA::DataLoader("ph_JKS_BDT_closest");
  TMVA::DataLoader *hp_JKS_BDT_highest = new TMVA::DataLoader("hp_JKS_BDT_highest");
  TMVA::DataLoader *hp_JKS_BDT_closest = new TMVA::DataLoader("hp_JKS_BDT_closest");
  TMVA::DataLoader *hh_JKS_BDT_highest = new TMVA::DataLoader("hh_JKS_BDT_highest");
  TMVA::DataLoader *hh_JKS_BDT_closest = new TMVA::DataLoader("hh_JKS_BDT_closest");

  TMVA::DataLoader *pp_s_vs_b_highest_with_KS=   new TMVA::DataLoader("pp_s_vs_b_highest_with_KS");
  TMVA::DataLoader *pp_s_vs_non_highest_with_KS= new TMVA::DataLoader("pp_s_vs_non_highest_with_KS");
  TMVA::DataLoader *pp_s_vs_all_highest_with_KS= new TMVA::DataLoader("pp_s_vs_all_highest_with_KS");
  TMVA::DataLoader *pp_s_vs_b_closest_with_KS=   new TMVA::DataLoader("pp_s_vs_b_closest_with_KS");
  TMVA::DataLoader *pp_s_vs_non_closest_with_KS= new TMVA::DataLoader("pp_s_vs_non_closest_with_KS");
  TMVA::DataLoader *pp_s_vs_all_closest_with_KS= new TMVA::DataLoader("pp_s_vs_all_closest_with_KS");

  TMVA::DataLoader *ph_s_vs_b_highest_with_KS=   new TMVA::DataLoader("ph_s_vs_b_highest_with_KS");
  TMVA::DataLoader *ph_s_vs_non_highest_with_KS= new TMVA::DataLoader("ph_s_vs_non_highest_with_KS");
  TMVA::DataLoader *ph_s_vs_all_highest_with_KS= new TMVA::DataLoader("ph_s_vs_all_highest_with_KS");
  TMVA::DataLoader *ph_s_vs_b_closest_with_KS=   new TMVA::DataLoader("ph_s_vs_b_closest_with_KS");
  TMVA::DataLoader *ph_s_vs_non_closest_with_KS= new TMVA::DataLoader("ph_s_vs_non_closest_with_KS");
  TMVA::DataLoader *ph_s_vs_all_closest_with_KS= new TMVA::DataLoader("ph_s_vs_all_closest_with_KS");

  TMVA::DataLoader *hp_s_vs_b_highest_with_KS=   new TMVA::DataLoader("hp_s_vs_b_highest_with_KS");
  TMVA::DataLoader *hp_s_vs_non_highest_with_KS= new TMVA::DataLoader("hp_s_vs_non_highest_with_KS");
  TMVA::DataLoader *hp_s_vs_all_highest_with_KS= new TMVA::DataLoader("hp_s_vs_all_highest_with_KS");
  TMVA::DataLoader *hp_s_vs_b_closest_with_KS=   new TMVA::DataLoader("hp_s_vs_b_closest_with_KS");
  TMVA::DataLoader *hp_s_vs_non_closest_with_KS= new TMVA::DataLoader("hp_s_vs_non_closest_with_KS");
  TMVA::DataLoader *hp_s_vs_all_closest_with_KS= new TMVA::DataLoader("hp_s_vs_all_closest_with_KS");

  TMVA::DataLoader *hh_s_vs_b_highest_with_KS=   new TMVA::DataLoader("hh_s_vs_b_highest_with_KS");
  TMVA::DataLoader *hh_s_vs_non_highest_with_KS= new TMVA::DataLoader("hh_s_vs_non_highest_with_KS");
  TMVA::DataLoader *hh_s_vs_all_highest_with_KS= new TMVA::DataLoader("hh_s_vs_all_highest_with_KS");
  TMVA::DataLoader *hh_s_vs_b_closest_with_KS=   new TMVA::DataLoader("hh_s_vs_b_closest_with_KS");
  TMVA::DataLoader *hh_s_vs_non_closest_with_KS= new TMVA::DataLoader("hh_s_vs_non_closest_with_KS");
  TMVA::DataLoader *hh_s_vs_all_closest_with_KS= new TMVA::DataLoader("hh_s_vs_all_closest_with_KS");

  //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
  //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";
  (TMVA::gConfig().GetVariablePlotting()).fMaxNumOfAllowedVariablesForScatterPlots = 99999;
  (TMVA::gConfig().GetVariablePlotting()).fNbinsXOfROCCurve = 1000;

  // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]

  addJetVariable(pp_JKS_BDT_highest);
  addJetVariable(pp_JKS_BDT_closest);
  addJetVariable(ph_JKS_BDT_highest);
  addJetVariable(ph_JKS_BDT_closest);
  addJetVariable(hp_JKS_BDT_highest);
  addJetVariable(hp_JKS_BDT_closest);
  addJetVariable(hh_JKS_BDT_highest);
  addJetVariable(hh_JKS_BDT_closest);

  addJetVariable(pp_s_vs_b_highest_with_KS); 
  addJetVariable(pp_s_vs_non_highest_with_KS); 
  addJetVariable(pp_s_vs_all_highest_with_KS);
  addJetVariable(pp_s_vs_b_closest_with_KS);
  addJetVariable(pp_s_vs_non_closest_with_KS);
  addJetVariable(pp_s_vs_all_closest_with_KS);

  addJetVariable(ph_s_vs_b_highest_with_KS);
  addJetVariable(ph_s_vs_non_highest_with_KS);
  addJetVariable(ph_s_vs_all_highest_with_KS);
  addJetVariable(ph_s_vs_b_closest_with_KS);
  addJetVariable(ph_s_vs_non_closest_with_KS);
  addJetVariable(ph_s_vs_all_closest_with_KS);

  addJetVariable(hp_s_vs_b_highest_with_KS);
  addJetVariable(hp_s_vs_non_highest_with_KS);
  addJetVariable(hp_s_vs_all_highest_with_KS);
  addJetVariable(hp_s_vs_b_closest_with_KS);
  addJetVariable(hp_s_vs_non_closest_with_KS);
  addJetVariable(hp_s_vs_all_closest_with_KS);

  addJetVariable(hh_s_vs_b_highest_with_KS);
  addJetVariable(hh_s_vs_non_highest_with_KS);
  addJetVariable(hh_s_vs_all_highest_with_KS);
  addJetVariable(hh_s_vs_b_closest_with_KS);
  addJetVariable(hh_s_vs_non_closest_with_KS);
  addJetVariable(hh_s_vs_all_closest_with_KS);

  addHadVariablePP(pp_JKS_BDT_highest);
  addHadVariablePP(pp_JKS_BDT_closest);
  addHadVariablePP(ph_JKS_BDT_highest);
  addHadVariablePP(ph_JKS_BDT_closest);
  addHadVariablePP(hp_JKS_BDT_highest);
  addHadVariablePP(hp_JKS_BDT_closest);
  addHadVariablePP(hh_JKS_BDT_highest);
  addHadVariablePP(hh_JKS_BDT_closest);

  addHadVariablePP(pp_s_vs_b_highest_with_KS);
  addHadVariablePP(pp_s_vs_non_highest_with_KS);
  addHadVariablePP(pp_s_vs_all_highest_with_KS);
  addHadVariablePP(pp_s_vs_b_closest_with_KS);
  addHadVariablePP(pp_s_vs_non_closest_with_KS);
  addHadVariablePP(pp_s_vs_all_closest_with_KS);

  addHadVariablePP(ph_s_vs_b_highest_with_KS);
  addHadVariablePP(ph_s_vs_non_highest_with_KS);
  addHadVariablePP(ph_s_vs_all_highest_with_KS);
  addHadVariablePP(ph_s_vs_b_closest_with_KS);
  addHadVariablePP(ph_s_vs_non_closest_with_KS);
  addHadVariablePP(ph_s_vs_all_closest_with_KS);

  addHadVariablePP(hp_s_vs_b_highest_with_KS);
  addHadVariablePP(hp_s_vs_non_highest_with_KS);
  addHadVariablePP(hp_s_vs_all_highest_with_KS);
  addHadVariablePP(hp_s_vs_b_closest_with_KS);
  addHadVariablePP(hp_s_vs_non_closest_with_KS);
  addHadVariablePP(hp_s_vs_all_closest_with_KS);

  addHadVariablePP(hh_s_vs_b_highest_with_KS);
  addHadVariablePP(hh_s_vs_non_highest_with_KS);
  addHadVariablePP(hh_s_vs_all_highest_with_KS);
  addHadVariablePP(hh_s_vs_b_closest_with_KS);
  addHadVariablePP(hh_s_vs_non_closest_with_KS);
  addHadVariablePP(hh_s_vs_all_closest_with_KS);
/*
  addHadVariableHH(pp_s_vs_b_highest_with_KS);
  addHadVariableHH(pp_s_vs_non_highest_with_KS);
  addHadVariableHH(pp_s_vs_all_highest_with_KS);
  addHadVariableHH(pp_s_vs_b_closest_with_KS);
  addHadVariableHH(pp_s_vs_non_closest_with_KS);
  addHadVariableHH(pp_s_vs_all_closest_with_KS);

  addHadVariableHH(ph_s_vs_b_highest_with_KS);
  addHadVariableHH(ph_s_vs_non_highest_with_KS);
  addHadVariableHH(ph_s_vs_all_highest_with_KS);
  addHadVariableHH(ph_s_vs_b_closest_with_KS);
  addHadVariableHH(ph_s_vs_non_closest_with_KS);
  addHadVariableHH(ph_s_vs_all_closest_with_KS);

  addHadVariableHH(hp_s_vs_b_highest_with_KS);
  addHadVariableHH(hp_s_vs_non_highest_with_KS);
  addHadVariableHH(hp_s_vs_all_highest_with_KS);
  addHadVariableHH(hp_s_vs_b_closest_with_KS);
  addHadVariableHH(hp_s_vs_non_closest_with_KS);
  addHadVariableHH(hp_s_vs_all_closest_with_KS);

  addHadVariableHH(hh_s_vs_b_highest_with_KS);
  addHadVariableHH(hh_s_vs_non_highest_with_KS);
  addHadVariableHH(hh_s_vs_all_highest_with_KS);
  addHadVariableHH(hh_s_vs_b_closest_with_KS);
  addHadVariableHH(hh_s_vs_non_closest_with_KS);
  addHadVariableHH(hh_s_vs_all_closest_with_KS);
*/
  // global event weights per tree (see below for setting event-wise weights)
  Double_t signalWeight     = 1.0; Double_t backgroundWeight = 1.0;

  // You can add an arbitrary number of signal or background trees
  addTree(pp_JKS_BDT_highest, bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight);
  addTree(pp_JKS_BDT_closest, bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight);
  addTree(ph_JKS_BDT_highest, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, signalWeight, backgroundWeight);
  addTree(ph_JKS_BDT_closest, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, signalWeight, backgroundWeight);
  addTree(hp_JKS_BDT_highest, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight);
  addTree(hp_JKS_BDT_closest, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight);
  addTree(hh_JKS_BDT_highest, bbars_bsbar_herwig_tree, bbars_bsbar_herwig_tree, signalWeight, backgroundWeight);
  addTree(hh_JKS_BDT_closest, bbars_bsbar_herwig_tree, bbars_bsbar_herwig_tree, signalWeight, backgroundWeight);

  addTree(pp_s_vs_b_highest_with_KS,   bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight);
  addTree(pp_s_vs_non_highest_with_KS, bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight);
  addTree(pp_s_vs_all_highest_with_KS, bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight);
  addTree(pp_s_vs_b_closest_with_KS,   bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight);
  addTree(pp_s_vs_non_closest_with_KS, bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight);
  addTree(pp_s_vs_all_closest_with_KS, bbars_bsbar_pythia_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight);

  addTree(ph_s_vs_b_highest_with_KS,   bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, signalWeight, backgroundWeight);
  addTree(ph_s_vs_non_highest_with_KS, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, signalWeight, backgroundWeight);
  addTree(ph_s_vs_all_highest_with_KS, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, signalWeight, backgroundWeight);
  addTree(ph_s_vs_b_closest_with_KS,   bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, signalWeight, backgroundWeight);
  addTree(ph_s_vs_non_closest_with_KS, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, signalWeight, backgroundWeight);
  addTree(ph_s_vs_all_closest_with_KS, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, signalWeight, backgroundWeight);

  addTree(hp_s_vs_b_highest_with_KS,   bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight);
  addTree(hp_s_vs_non_highest_with_KS, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight);
  addTree(hp_s_vs_all_highest_with_KS, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight);
  addTree(hp_s_vs_b_closest_with_KS,   bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight);
  addTree(hp_s_vs_non_closest_with_KS, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight);
  addTree(hp_s_vs_all_closest_with_KS, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, bbars_bsbar_herwig_tree, bbars_bsbar_pythia_tree, signalWeight, backgroundWeight);

  addTree(hh_s_vs_b_highest_with_KS,   bbars_bsbar_herwig_tree, bbars_bsbar_herwig_tree, signalWeight, backgroundWeight);
  addTree(hh_s_vs_non_highest_with_KS, bbars_bsbar_herwig_tree, bbars_bsbar_herwig_tree, signalWeight, backgroundWeight);
  addTree(hh_s_vs_all_highest_with_KS, bbars_bsbar_herwig_tree, bbars_bsbar_herwig_tree, signalWeight, backgroundWeight);
  addTree(hh_s_vs_b_closest_with_KS,   bbars_bsbar_herwig_tree, bbars_bsbar_herwig_tree, signalWeight, backgroundWeight);
  addTree(hh_s_vs_non_closest_with_KS, bbars_bsbar_herwig_tree, bbars_bsbar_herwig_tree, signalWeight, backgroundWeight);
  addTree(hh_s_vs_all_closest_with_KS, bbars_bsbar_herwig_tree, bbars_bsbar_herwig_tree, signalWeight, backgroundWeight);

  //dataloader->SetBackgroundWeightExpression( "weight" );

/*
  TCut cut_s_highest_with_KS   = "isSJet == 1 && isHighest == 1";
  TCut cut_b_highest_with_KS   = "isBJet == 1 && isHighest == 1";
  TCut cut_non_highest_with_KS = "abs(KS_isFrom_pythia) != 5 && isSJet == 0 && isBJet == 0 && isHighest == 1";
  TCut cut_all_highest_with_KS = "isSJet != 1 && isHighest == 1";

  TCut cut_s_closest_with_KS   = "isSJet == 1 && isClosestToLep == 1";
  TCut cut_b_closest_with_KS   = "isBJet == 1 && isClosestToLep == 1";
  TCut cut_non_closest_with_KS = "abs(KS_isFrom_pythia) != 5 && isSJet == 0 && isBJet == 0 && isClosestToLep == 1";
  TCut cut_all_closest_with_KS = "isSJet != 1 && isClosestToLep == 1";
*/

  TCut cut_s_highest_with_KS   = "isHighest && KS_idx_pp != -99 && KS_nMatched_pp == 2 && isSJet";
  TCut cut_b_highest_with_KS   = "isHighest && KS_idx_pp != -99 && KS_nMatched_pp != 2 && isBJet";
  TCut cut_non_highest_with_KS = "isHighest && KS_idx_pp != -99 && KS_nMatched_pp != 2 && !isSJet && !isBJet";
  TCut cut_all_highest_with_KS = "isHighest && KS_idx_pp != -99 && KS_nMatched_pp != 2 && !isSJet";

  TCut cut_s_closest_with_KS   = "isClosestToLep && KS_idx_pp != -99 && KS_nMatched_pp == 2 && isSJet";
  TCut cut_b_closest_with_KS   = "isClosestToLep && KS_idx_pp != -99 && KS_nMatched_pp != 2 && isBJet";
  TCut cut_non_closest_with_KS = "isClosestToLep && KS_idx_pp != -99 && KS_nMatched_pp != 2 && !isSJet && !isBJet";
  TCut cut_all_closest_with_KS = "isClosestToLep && KS_idx_pp != -99 && KS_nMatched_pp != 2 && !isSJet";

  TCut cut_JKS_BDT_s_highest   = "isHighest && KS_idx_pp != -99 && KS_best_bdt_pp > 0.0508 && isSJet";
  TCut cut_JKS_BDT_all_highest = "isHighest && KS_idx_pp != -99 && KS_best_bdt_pp > 0.0508 && !isSJet";
  TCut cut_JKS_BDT_s_closest   = "isClosestToLep && KS_idx_pp != -99 && KS_best_bdt_pp > 0.0508 && isSJet";
  TCut cut_JKS_BDT_all_closest = "isClosestToLep && KS_idx_pp != -99 && KS_best_bdt_pp > 0.0508 && !isSJet";

  TString opt_1 = "nTrain_Signal=3500:nTrain_Background=3500:SplitMode=Random:NormMode=NumEvents:!V";
  TString opt_2 = "nTrain_Signal=2000:nTrain_Background=2000:SplitMode=Random:NormMode=NumEvents:!V";
  TString opt_3 = "nTrain_Signal=600:nTrain_Background=7000:SplitMode=Random:NormMode=NumEvents:!V";
  TString opt_4 = "nTrain_Signal=400:nTrain_Background=2000:SplitMode=Random:NormMode=NumEvents:!V";
  TString opt_5 = "nTrain_Signal=350:nTrain_Background=700:SplitMode=Random:NormMode=NumEvents:!V";
 
  pp_JKS_BDT_highest->PrepareTrainingAndTestTree( cut_JKS_BDT_s_highest, cut_JKS_BDT_all_highest, opt_5 );
  pp_JKS_BDT_closest->PrepareTrainingAndTestTree( cut_JKS_BDT_s_closest, cut_JKS_BDT_all_closest, opt_5 );
  ph_JKS_BDT_highest->PrepareTrainingAndTestTree( cut_JKS_BDT_s_highest, cut_JKS_BDT_all_highest, opt_5 );
  ph_JKS_BDT_closest->PrepareTrainingAndTestTree( cut_JKS_BDT_s_closest, cut_JKS_BDT_all_closest, opt_5 );
  hp_JKS_BDT_highest->PrepareTrainingAndTestTree( cut_JKS_BDT_s_highest, cut_JKS_BDT_all_highest, opt_5 );
  hp_JKS_BDT_closest->PrepareTrainingAndTestTree( cut_JKS_BDT_s_closest, cut_JKS_BDT_all_closest, opt_5 );
  hh_JKS_BDT_highest->PrepareTrainingAndTestTree( cut_JKS_BDT_s_highest, cut_JKS_BDT_all_highest, opt_5 );
  hh_JKS_BDT_closest->PrepareTrainingAndTestTree( cut_JKS_BDT_s_closest, cut_JKS_BDT_all_closest, opt_5 );

  pp_s_vs_b_highest_with_KS->PrepareTrainingAndTestTree(   cut_s_highest_with_KS, cut_b_highest_with_KS,   opt_4 ); 
  pp_s_vs_non_highest_with_KS->PrepareTrainingAndTestTree( cut_s_highest_with_KS, cut_non_highest_with_KS, opt_4 ); 
  pp_s_vs_all_highest_with_KS->PrepareTrainingAndTestTree( cut_s_highest_with_KS, cut_all_highest_with_KS, opt_4 ); 
  pp_s_vs_b_closest_with_KS->PrepareTrainingAndTestTree(   cut_s_closest_with_KS, cut_b_closest_with_KS,   opt_4 );
  pp_s_vs_non_closest_with_KS->PrepareTrainingAndTestTree( cut_s_closest_with_KS, cut_non_closest_with_KS, opt_4 );  
  pp_s_vs_all_closest_with_KS->PrepareTrainingAndTestTree( cut_s_closest_with_KS, cut_all_closest_with_KS, opt_4 );  

  ph_s_vs_b_highest_with_KS->PrepareTrainingAndTestTree(   cut_s_highest_with_KS, cut_b_highest_with_KS,   opt_4 );
  ph_s_vs_non_highest_with_KS->PrepareTrainingAndTestTree( cut_s_highest_with_KS, cut_non_highest_with_KS, opt_4 );
  ph_s_vs_all_highest_with_KS->PrepareTrainingAndTestTree( cut_s_highest_with_KS, cut_all_highest_with_KS, opt_4 );
  ph_s_vs_b_closest_with_KS->PrepareTrainingAndTestTree(   cut_s_closest_with_KS, cut_b_closest_with_KS,   opt_4 );
  ph_s_vs_non_closest_with_KS->PrepareTrainingAndTestTree( cut_s_closest_with_KS, cut_non_closest_with_KS, opt_4 );
  ph_s_vs_all_closest_with_KS->PrepareTrainingAndTestTree( cut_s_closest_with_KS, cut_all_closest_with_KS, opt_4 );

  hp_s_vs_b_highest_with_KS->PrepareTrainingAndTestTree(   cut_s_highest_with_KS, cut_b_highest_with_KS,   opt_4 ); 
  hp_s_vs_non_highest_with_KS->PrepareTrainingAndTestTree( cut_s_highest_with_KS, cut_non_highest_with_KS, opt_4 );  
  hp_s_vs_all_highest_with_KS->PrepareTrainingAndTestTree( cut_s_highest_with_KS, cut_all_highest_with_KS, opt_4 ); 
  hp_s_vs_b_closest_with_KS->PrepareTrainingAndTestTree(   cut_s_closest_with_KS, cut_b_closest_with_KS,   opt_4 );
  hp_s_vs_non_closest_with_KS->PrepareTrainingAndTestTree( cut_s_closest_with_KS, cut_non_closest_with_KS, opt_4 ); 
  hp_s_vs_all_closest_with_KS->PrepareTrainingAndTestTree( cut_s_closest_with_KS, cut_all_closest_with_KS, opt_4 ); 

  hh_s_vs_b_highest_with_KS->PrepareTrainingAndTestTree(   cut_s_highest_with_KS, cut_b_highest_with_KS,   opt_4 ); 
  hh_s_vs_non_highest_with_KS->PrepareTrainingAndTestTree( cut_s_highest_with_KS, cut_non_highest_with_KS, opt_4 ); 
  hh_s_vs_all_highest_with_KS->PrepareTrainingAndTestTree( cut_s_highest_with_KS, cut_all_highest_with_KS, opt_4 );  
  hh_s_vs_b_closest_with_KS->PrepareTrainingAndTestTree(   cut_s_closest_with_KS, cut_b_closest_with_KS,   opt_4 );
  hh_s_vs_non_closest_with_KS->PrepareTrainingAndTestTree( cut_s_closest_with_KS, cut_non_closest_with_KS, opt_4 );   
  hh_s_vs_all_closest_with_KS->PrepareTrainingAndTestTree( cut_s_closest_with_KS, cut_all_closest_with_KS, opt_4 );  

  addMethod(Use, factory, pp_JKS_BDT_highest); 
  addMethod(Use, factory, pp_JKS_BDT_closest);
  addMethod(Use, factory, ph_JKS_BDT_highest);
  addMethod(Use, factory, ph_JKS_BDT_closest);
  addMethod(Use, factory, hp_JKS_BDT_highest);
  addMethod(Use, factory, hp_JKS_BDT_closest);
  addMethod(Use, factory, hh_JKS_BDT_highest);
  addMethod(Use, factory, hh_JKS_BDT_closest);

  addMethod(Use, factory, pp_s_vs_b_highest_with_KS);
  addMethod(Use, factory, pp_s_vs_non_highest_with_KS);
  addMethod(Use, factory, pp_s_vs_all_highest_with_KS);
  addMethod(Use, factory, pp_s_vs_b_closest_with_KS);
  addMethod(Use, factory, pp_s_vs_non_closest_with_KS);
  addMethod(Use, factory, pp_s_vs_all_closest_with_KS);

  addMethod(Use, factory, ph_s_vs_b_highest_with_KS);
  addMethod(Use, factory, ph_s_vs_non_highest_with_KS);
  addMethod(Use, factory, ph_s_vs_all_highest_with_KS);
  addMethod(Use, factory, ph_s_vs_b_closest_with_KS);
  addMethod(Use, factory, ph_s_vs_non_closest_with_KS);
  addMethod(Use, factory, ph_s_vs_all_closest_with_KS);

  addMethod(Use, factory, hp_s_vs_b_highest_with_KS);
  addMethod(Use, factory, hp_s_vs_non_highest_with_KS);
  addMethod(Use, factory, hp_s_vs_all_highest_with_KS);
  addMethod(Use, factory, hp_s_vs_b_closest_with_KS);
  addMethod(Use, factory, hp_s_vs_non_closest_with_KS);
  addMethod(Use, factory, hp_s_vs_all_closest_with_KS);

  addMethod(Use, factory, hh_s_vs_b_highest_with_KS);
  addMethod(Use, factory, hh_s_vs_non_highest_with_KS);
  addMethod(Use, factory, hh_s_vs_all_highest_with_KS);
  addMethod(Use, factory, hh_s_vs_b_closest_with_KS);
  addMethod(Use, factory, hh_s_vs_non_closest_with_KS);
  addMethod(Use, factory, hh_s_vs_all_closest_with_KS);

  // For an example of the category classifier usage, see: vts_dR_04_Jet_With_KSCategory
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
  std::cout << "==> vts_dR_04_Jet_With_KS is done!" << std::endl;

  delete factory;
  delete pp_JKS_BDT_highest;
  delete pp_JKS_BDT_closest;
  delete ph_JKS_BDT_highest;
  delete ph_JKS_BDT_closest;
  delete hp_JKS_BDT_highest;
  delete hp_JKS_BDT_closest;
  delete hh_JKS_BDT_highest;
  delete hh_JKS_BDT_closest;
  delete pp_s_vs_b_highest_with_KS;
  delete pp_s_vs_non_highest_with_KS;
  delete pp_s_vs_all_highest_with_KS;
  delete pp_s_vs_b_closest_with_KS;
  delete pp_s_vs_non_closest_with_KS;
  delete pp_s_vs_all_closest_with_KS;
  delete ph_s_vs_b_highest_with_KS;
  delete ph_s_vs_non_highest_with_KS;
  delete ph_s_vs_all_highest_with_KS;
  delete ph_s_vs_b_closest_with_KS;
  delete ph_s_vs_non_closest_with_KS;
  delete ph_s_vs_all_closest_with_KS;
  delete hp_s_vs_b_highest_with_KS;
  delete hp_s_vs_non_highest_with_KS;
  delete hp_s_vs_all_highest_with_KS;
  delete hp_s_vs_b_closest_with_KS;
  delete hp_s_vs_non_closest_with_KS;
  delete hp_s_vs_all_closest_with_KS;
  delete hh_s_vs_b_highest_with_KS;
  delete hh_s_vs_non_highest_with_KS;
  delete hh_s_vs_all_highest_with_KS;
  delete hh_s_vs_b_closest_with_KS;
  delete hh_s_vs_non_closest_with_KS;
  delete hh_s_vs_all_closest_with_KS;

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
  return vts_dR_04_Jet_With_KS(methodList);
}
