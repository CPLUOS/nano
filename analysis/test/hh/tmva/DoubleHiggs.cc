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

int DoubleHiggs( TString myMethodList = "" )
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
   std::cout << "==> Start DoubleHiggs" << std::endl;

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

   TString directory = "/home/sunyoung/nanoAOD/src/nano/analysis/test/hh/result/";
   TString sig_file = "HH_SM.root";
   TString bg_file = "TT.root";

   input = TFile::Open( directory+sig_file );
   input2 = TFile::Open( directory+bg_file );

   std::cout << "double Higgs using file : " << input->GetName() << std::endl;

   // Register the training and test trees
   TTree *signalTree     = (TTree*)input->Get("events");
   TTree *backgroundTree     = (TTree*)input2->Get("events");

   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TString outfileName( "doubleHiggs_tmva" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   TMVA::Factory *factory = new TMVA::Factory( "DoubleHiggs", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification" );
//                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

   TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset_"+sig_file);
   TMVA::DataLoader *dataloader2=new TMVA::DataLoader("dataset_"+bg_file);

   //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
   //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";
   (TMVA::gConfig().GetVariablePlotting()).fMaxNumOfAllowedVariablesForScatterPlots = 99999;
   (TMVA::gConfig().GetVariablePlotting()).fNbinsXOfROCCurve = 1000;

   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]

   dataloader->AddVariable( "ll_deltaR", 'F' );
   dataloader->AddVariable( "ll.Pt()", 'F' );
   dataloader->AddVariable( "ll.M()", 'F' );
   dataloader->AddVariable( "bb_deltaR", 'F' );
   dataloader->AddVariable( "bb.Pt()", 'F' );
   dataloader->AddVariable( "bb.M()", 'F' );
   dataloader->AddVariable( "bl_min_deltaR", 'F' );
   dataloader->AddVariable( "bbll_deltaR", 'F' );
   dataloader->AddVariable( "bbll_deltaPhi", 'F' );
   dataloader->AddVariable( "mT", 'F');
   dataloader->AddVariable( "basic_MT2_332_bbll", 'F' );
   //dataloader->AddSpectator( "mass_had", 'F' );
   //dataloader2->AddVariable( "log_d_had := log(d_had)", 'F' );
   dataloader2->AddVariable( "ll_deltaR", 'F' );
   dataloader2->AddVariable( "ll.Pt()", 'F' );
   dataloader2->AddVariable( "ll.M()", 'F' );
   dataloader2->AddVariable( "bb_deltaR", 'F' );
   dataloader2->AddVariable( "bb.Pt()", 'F' );
   dataloader2->AddVariable( "bb.M()", 'F' );
   dataloader2->AddVariable( "bl_min_deltaR", 'F' );
   dataloader2->AddVariable( "bbll_deltaR", 'F' );
   dataloader2->AddVariable( "bbll_deltaPhi", 'F' );
   dataloader2->AddVariable( "mT", 'F');
   dataloader2->AddVariable( "basic_MT2_332_bbll", 'F' );
   //dataloader2->AddSpectator( "mass_had", 'F' );

   // global event weights per tree (see below for setting event-wise weights)
   Double_t signalWeight     = 1.0;
   Double_t backgroundWeight = 1.0;

   // You can add an arbitrary number of signal or background trees
   dataloader->AddSignalTree    ( signalTree,     signalWeight );
   dataloader->AddBackgroundTree( backgroundTree, backgroundWeight );

   //dataloader->SetBackgroundWeightExpression( "weight" );

   //TCut mycuts = "abs(isFrom_had) == 3 && fabs(mass_had - 0.5) < 0.05"; 
     TCut mycuts = "step==4";
     TCut mycutb = "step==4";
   //TCut mycutb = "abs(isFrom_had) != 3 && abs(isFrom_had) != 99 && fabs(mass_had - 0.5) < 0.05"; 

   //TCut mycuts2 = "abs(isFrom_had) == 3";
  // TCut mycutb2 = "abs(isFrom_had) == 5";    


   dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,
//                                       "nTrain_Signal=150:nTrain_Background=700:SplitMode=Random:NormMode=NumEvents:!V" );
                                       "nTrain_Signal=5000:nTrain_Background=30000:SplitMode=Random:NormMode=NumEvents:!V" );

   //if (Use["MLPBFGS"])
     // factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );
   //if (Use["TMlpANN"])
     // factory->BookMethod( dataloader, TMVA::Types::kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  );
   if (Use["BDTG"]) // Gradient Boost
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG",
                           "!H:!V:NTrees=700:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );
   if (Use["BDT"])  // Adaptive Boost
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=700:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );


   // For an example of the category classifier usage, see: DoubleHiggsCategory
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
   std::cout << "==> DoubleHiggs is done!" << std::endl;

   delete factory;
   delete dataloader;
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
   return DoubleHiggs(methodList);
}
