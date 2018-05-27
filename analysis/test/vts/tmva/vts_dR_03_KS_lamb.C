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

int vts_dR_03_KS_lamb( TString myMethodList = "" )
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
   std::cout << "==> Start vts_dR_03_KS_lamb" << std::endl;

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
   // (it is also possible to use ASCII format as input -> see TMVA Users Guide)
   TFile *input(0);
   TFile *input2(0);

   TString tsw = "/xrootd/store/user/wjjang/test/tsw/";
   TString tswf = "sum_235_dR_03_KS_lamb_no_ES.root";
   TString tswf2 = "old_235_dR_03_KS_lamb_no_ES.root";
   TString tswf3 = "sum_235_dR_03_KS_lamb_ES.root";
   TString tswf4 = "old_235_dR_03_KS_lamb_ES.root";

   TString tt01j_bbars = "/xrootd/store/user/wjjang/test/tt01j_bbars/";
   TString bbars = "sum_50_dR_03_KS_lamb_no_ES.root";
   TString bbars2 = "old_50_dR_03_KS_lamb_no_ES.root";
   TString bbars3 = "sum_50_dR_03_KS_lamb_ES.root";
   TString bbars4 = "old_50_dR_03_KS_lamb_ES.root";

   TString tt01j_bsbar = "/xrootd/store/user/wjjang/test/tt01j_bsbar/";
   TString bsbar = "sum_50_dR_03_KS_lamb_no_ES.root";
   TString bsbar2 = "old_50_dR_03_KS_lamb_no_ES.root";
   TString bsbar3 = "sum_50_dR_03_KS_lamb_ES.root";
   TString bsbar4 = "old_50_dR_03_KS_lamb_ES.root";

   input = TFile::Open( tsw+tswf );
   input2 = TFile::Open( tsw+tswf2 );

   std::cout << "--- vts_dR_03_KS_lamb       : Using input file: " << input->GetName() << std::endl;

   // Register the training and test trees
   TTree *signalTree     = (TTree*)input->Get("event");
   TTree *background     = (TTree*)input->Get("event");
   TTree *signalTree2     = (TTree*)input2->Get("event");
   TTree *background2     = (TTree*)input2->Get("event");

   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TString outfileName( "vts_dR_03_KS_lamb.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   TMVA::Factory *factory = new TMVA::Factory( "vts_dR_03_KS_lamb", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification" );
//                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

   TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset_"+tswf);
   TMVA::DataLoader *dataloader2=new TMVA::DataLoader("dataset_"+tswf2);

   //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
   //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";
   (TMVA::gConfig().GetVariablePlotting()).fMaxNumOfAllowedVariablesForScatterPlots = 99999;
   (TMVA::gConfig().GetVariablePlotting()).fNbinsXOfROCCurve = 1000;

   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]

//   dataloader->AddVariable( "mass_had", 'F' );

   dataloader->AddVariable( "log_d_had := log(d_had)", 'F' );
   dataloader->AddVariable( "log_x_had := log(x_had)", 'F' );
   dataloader->AddVariable( "log_dr_had := log(dr_had)", 'F' );
   dataloader->AddVariable( "log_lxySig_had := log(lxySig_had)", 'F' );
   dataloader->AddVariable( "log_lxy_had := log(lxy_had)", 'F' );
   dataloader->AddVariable( "log_chi2_had := log(chi2_had)", 'F' );
   dataloader->AddVariable( "log_angleXY_had := log(angleXY_had)", 'F' );
   dataloader->AddVariable( "log_dca_had := log(dca_had)", 'F' );

//   dataloader->AddVariable( "log_l3D_had := log(l3D_had)", 'F' );
//   dataloader->AddVariable( "log_l3DSig_had := log(l3DSig_had)", 'F' );
   dataloader->AddVariable( "log_legDR_had := log(legDR_had)", 'F' );
   dataloader->AddVariable( "log_dau1_ipsigXY := log(dau1_ipsigXY_had)", 'F' );
//   dataloader->AddVariable( "log_dau2_ipsigXY := log(dau2_ipsigXY_had)", 'F' );

   dataloader->AddVariable( "btagCSVV2_Jet", 'F' );
   dataloader->AddVariable( "btagCMVA_Jet", 'F' );
//   dataloader->AddVariable( "btagDeepB_Jet", 'F' );
//   dataloader->AddVariable( "btagDeepC_Jet", 'F' );

/*
   dataloader->AddVariable( "pt_Jet", 'F' );
   dataloader->AddVariable( "nConstituents_Jet", 'I' );
   dataloader->AddVariable( "nElectrons_Jet", 'I' );
   dataloader->AddVariable( "nMuons_Jet", 'I' );
*/

   dataloader->AddSpectator( "mass_had", 'F' );


//   dataloader2->AddVariable( "mass_had", 'F' );

   dataloader2->AddVariable( "log_d_had := log(d_had)", 'F' );
   dataloader2->AddVariable( "log_x_had := log(x_had)", 'F' );
//   dataloader2->AddVariable( "log_dr_had := log(dr_had)", 'F' );
   dataloader2->AddVariable( "log_lxySig_had := log(lxySig_had)", 'F' );
   dataloader2->AddVariable( "log_lxy_had := log(lxy_had)", 'F' );
   dataloader2->AddVariable( "log_chi2_had := log(chi2_had)", 'F' );
   dataloader2->AddVariable( "log_angleXY_had := log(angleXY_had)", 'F' );
   dataloader2->AddVariable( "log_dca_had := log(dca_had)", 'F' );

//   dataloader2->AddVariable( "log_l3D_had := log(l3D_had)", 'F' );
//   dataloader2->AddVariable( "log_l3DSig_had := log(l3DSig_had)", 'F' );
   dataloader2->AddVariable( "log_legDR_had := log(legDR_had)", 'F' );
   dataloader2->AddVariable( "log_dau1_ipsigXY := log(dau1_ipsigXY_had)", 'F' );
//   dataloader2->AddVariable( "log_dau2_ipsigXY := log(dau2_ipsigXY_had)", 'F' );

   dataloader2->AddVariable( "btagCSVV2_Jet", 'F' );
   dataloader2->AddVariable( "btagCMVA_Jet", 'F' );
//   dataloader2->AddVariable( "btagDeepB_Jet", 'F' );
//   dataloader2->AddVariable( "btagDeepC_Jet", 'F' );

/*
   dataloader2->AddVariable( "pt_Jet", 'F' );
   dataloader2->AddVariable( "nConstituents_Jet", 'I' );
   dataloader2->AddVariable( "nElectrons_Jet", 'I' );
   dataloader2->AddVariable( "nMuons_Jet", 'I' );
*/

   dataloader2->AddSpectator( "mass_had", 'F' );



   // global event weights per tree (see below for setting event-wise weights)
   Double_t signalWeight     = 1.0;
   Double_t backgroundWeight = 1.0;

   // You can add an arbitrary number of signal or background trees
   dataloader->AddSignalTree    ( signalTree,     signalWeight );
   dataloader->AddBackgroundTree( background, backgroundWeight );

   dataloader2->AddSignalTree    ( signalTree2,     signalWeight );
   dataloader2->AddBackgroundTree( background2, backgroundWeight );

   //dataloader->SetBackgroundWeightExpression( "weight" );

   TCut mycuts = "abs(isFrom_had) == 3 && fabs(mass_had - 0.5) < 0.05"; 
   TCut mycutb = "abs(isFrom_had) != 3 && abs(isFrom_had) != 99 && fabs(mass_had - 0.5) < 0.05"; 

   TCut mycuts2 = "abs(isFrom_had) == 3";
   TCut mycutb2 = "abs(isFrom_had) == 5";    


   dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,
//                                       "nTrain_Signal=150:nTrain_Background=700:SplitMode=Random:NormMode=NumEvents:!V" );
                                       "nTrain_Signal=2800:nTrain_Background=16000:SplitMode=Random:NormMode=NumEvents:!V" );

   dataloader2->PrepareTrainingAndTestTree( mycuts2, mycutb2,
//                                       "nTrain_Signal=150:nTrain_Background=700:SplitMode=Random:NormMode=NumEvents:!V" );
                                       "nTrain_Signal=2800:nTrain_Background=16000:SplitMode=Random:NormMode=NumEvents:!V" );
                                       

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

   if (Use["MLPBFGS"])
      factory->BookMethod( dataloader2, TMVA::Types::kMLP, "MLPBFGS2", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );
   if (Use["TMlpANN"])
      factory->BookMethod( dataloader2, TMVA::Types::kTMlpANN, "TMlpANN2", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  );
   if (Use["BDTG"]) // Gradient Boost
      factory->BookMethod( dataloader2, TMVA::Types::kBDT, "BDTG2",
                           "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );
   if (Use["BDT"])  // Adaptive Boost
      factory->BookMethod( dataloader2, TMVA::Types::kBDT, "BDT2",
                           "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );


   // For an example of the category classifier usage, see: vts_dR_03_KS_lambCategory
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
   std::cout << "==> vts_dR_03_KS_lamb is done!" << std::endl;

   delete factory;
   delete dataloader;
   delete dataloader2;
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
   return vts_dR_03_KS_lamb(methodList);
}
