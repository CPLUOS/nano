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

int TMVAVts( TString myMethodList = "" )
{
   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

 //  Use["Likelihood"]      = 1;
/*
   Use["LikelihoodD"]     = 1; // the "D" extension indicates decorrelated input variables (see option strings)
   Use["LikelihoodPCA"]   = 1; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
   Use["LikelihoodKDE"]   = 1;
   Use["LikelihoodMIX"]   = 1;


   Use["PDERS"]           = 1;
   Use["PDERSD"]          = 1;
   Use["PDERSPCA"]        = 1;
   Use["PDEFoam"]         = 1;
   Use["PDEFoamBoost"]    = 1; // uses generalised MVA method boosting
*/
//   Use["KNN"]             = 1; // k-nearest neighbour method
//   Use["LD"]              = 1; // Linear Discriminant identical to Fisher
//   Use["Fisher"]          = 1;
/*
   Use["FisherG"]         = 1;
   Use["BoostedFisher"]   = 1; // uses generalised MVA method boosting
   Use["HMatrix"]         = 1;
*/
//   Use["MLP"]             = 1; // Recommended ANN
   Use["MLPBFGS"]         = 1; // Recommended ANN with optional training method
//   Use["MLPBNN"]          = 1; // Recommended ANN with BFGS training method and bayesian regulator
//   Use["CFMlpANN"]        = 1; // Depreciated ANN from ALEPH
   Use["TMlpANN"]         = 1; // ROOT's own ANN
//   Use["RuleFit"]         = 1;

   // Boosted Decision Trees
   Use["BDT"]             = 1; // uses Adaptive Boost
   Use["BDTG"]            = 1; // uses Gradient Boost
//   Use["BDTB"]            = 1; // uses Bagging
//   Use["BDTD"]            = 1; // decorrelation + Adaptive Boost
//   Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting

/*
   Use["MLPBFGS_sig"] = 1;
   Use["MLPBFGS_hidden_10"] = 1;
   Use["MLPBFGS_200"] = 1;
   Use["MLPBFGS_400"] = 1;
//   Use["MLPReLU"] = 1;

   Use["Grad_nT1000_sh0.05_BSF0.7_MD2"] = 1;
   Use["Grad_nT1000_sh0.05_BSF0.7_MD3"] = 1;
   Use["Grad_nT1000_sh0.05_BSF0.7_MD4"] = 1;
   Use["Grad_nT1000_sh0.05_BSF0.6_MD2"] = 1;
   Use["Grad_nT1000_sh0.05_BSF0.6_MD3"] = 1;
   Use["Grad_nT1000_sh0.05_BSF0.6_MD4"] = 1;
   Use["Grad_nT1000_sh0.05_BSF0.4_MD3"] = 1;  
   Use["Grad_nT1000_sh0.05_BSF0.4_MD4"] = 1;
   Use["Grad_nT1000_sh0.03_BSF0.7_MD3"] = 1;
   Use["Grad_nT1000_sh0.03_BSF0.7_MD4"] = 1; 
   Use["Grad_nT1000_sh0.03_BSF0.6_MD3"] = 1; 
   Use["Grad_nT1000_sh0.03_BSF0.6_MD4"] = 1; 
   Use["Grad_nT1000_sh0.03_BSF0.4_MD3"] = 1; 
   Use["Grad_nT1000_sh0.03_BSF0.4_MD4"] = 1;
*/

   //
   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVAVts" << std::endl;

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
   //TString fname = "/cms/ldap_home/wjjang/wj_nanoAOD_CMSSW_9_4_4/src/nano/analysis/test/Batch/nanotree.root";
   //TString fname = "/xrootd/store/user/tt8888tt/a/a/nanotree.root";
   TString fname = "/xrootd/store/user/wjjang/test/tsw/sum_235_dR_0.3_onlyKS_noCut.root";
   input = TFile::Open( fname );
   std::cout << "--- TMVAVts       : Using input file: " << input->GetName() << std::endl;

   // Register the training and test trees

   //TTree *signalTree     = (TTree*)input->Get("tsw");
   //TTree *background     = (TTree*)input->Get("tsw");

   TTree *signalTree     = (TTree*)input->Get("events");
   TTree *background     = (TTree*)input->Get("events");

   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TString outfileName( "TMVAVts.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   // Create the factory object. Later you can choose the methods
   // whose performance you'd like to investigate. The factory is
   // the only TMVA object you have to interact with
   //
   // The first argument is the base of the name of all the
   // weightfiles in the directory weight/
   //
   // The second argument is the output file for the training results
   // All TMVA output can be suppressed by removing the "!" (not) in
   // front of the "Silent" argument in the option string
   TMVA::Factory *factory = new TMVA::Factory( "TMVAVts", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification" );
//                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

   TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");
   //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
   //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";

   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]

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
   dataloader->AddVariable( "btagDeepB_Jet", 'F' );
   dataloader->AddVariable( "btagDeepC_Jet", 'F' );

/*
   dataloader->AddVariable( "pt_Jet", 'F' );
   dataloader->AddVariable( "nConstituents_Jet", 'I' );
   dataloader->AddVariable( "nElectrons_Jet", 'I' );
   dataloader->AddVariable( "nMuons_Jet", 'I' );
*/
   // global event weights per tree (see below for setting event-wise weights)
   Double_t signalWeight     = 1.0;
   Double_t backgroundWeight = 1.0;

   // You can add an arbitrary number of signal or background trees
   dataloader->AddSignalTree    ( signalTree,     signalWeight );
   dataloader->AddBackgroundTree( background, backgroundWeight );

   //dataloader->SetBackgroundWeightExpression( "weight" );

   TCut mycuts = "abs(isFrom_had) == 3 && (chi2_had < 3 && dca_had < 1 && angleXY_had > 0.98 && lxySig_had > 5)";// && l3DSig_had < 1000000000000. && lxySig_had < 1000000000000."; 
   TCut mycutb = "abs(isFrom_had) == 5 && (chi2_had < 3 && dca_had < 1 && angleXY_had > 0.98 && lxySig_had > 5)";// && l3DSig_had < 1000000000000. && lxySig_had < 1000000000000."; 

   dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,
                                       "nTrain_Signal=2800:nTrain_Background=16000:SplitMode=Random:NormMode=NumEvents:!V" );

   // Cut optimisation
   if (Use["Cuts"])
      factory->BookMethod( dataloader, TMVA::Types::kCuts, "Cuts",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );

   if (Use["CutsD"])
      factory->BookMethod( dataloader, TMVA::Types::kCuts, "CutsD",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" );

   if (Use["CutsPCA"])
      factory->BookMethod( dataloader, TMVA::Types::kCuts, "CutsPCA",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA" );

   if (Use["CutsGA"])
      factory->BookMethod( dataloader, TMVA::Types::kCuts, "CutsGA",
                           "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" );

   if (Use["CutsSA"])
      factory->BookMethod( dataloader, TMVA::Types::kCuts, "CutsSA",
                           "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );
   if (Use["Likelihood"])
      factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "Likelihood",
                           "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" );
   if (Use["LikelihoodD"])
      factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "LikelihoodD",
                           "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" );
   if (Use["LikelihoodPCA"])
      factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "LikelihoodPCA",
                           "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA" );
   if (Use["LikelihoodKDE"])
      factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "LikelihoodKDE",
                           "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50" );
   if (Use["LikelihoodMIX"])
      factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "LikelihoodMIX",
                           "!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50" );
   if (Use["PDERS"])
      factory->BookMethod( dataloader, TMVA::Types::kPDERS, "PDERS",
                           "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" );

   if (Use["PDERSD"])
      factory->BookMethod( dataloader, TMVA::Types::kPDERS, "PDERSD",
                           "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate" );

   if (Use["PDERSPCA"])
      factory->BookMethod( dataloader, TMVA::Types::kPDERS, "PDERSPCA",
                           "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA" );
   if (Use["PDEFoam"])
      factory->BookMethod( dataloader, TMVA::Types::kPDEFoam, "PDEFoam",
                           "!H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" );

   if (Use["PDEFoamBoost"])
      factory->BookMethod( dataloader, TMVA::Types::kPDEFoam, "PDEFoamBoost",
                           "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T" );

  if (Use["KNN"])
      factory->BookMethod( dataloader, TMVA::Types::kKNN, "KNN",
                           "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );

   if (Use["HMatrix"])
      factory->BookMethod( dataloader, TMVA::Types::kHMatrix, "HMatrix", "!H:!V:VarTransform=None" );

  if (Use["LD"])
      factory->BookMethod( dataloader, TMVA::Types::kLD, "LD", "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

   if (Use["Fisher"])
      factory->BookMethod( dataloader, TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

   if (Use["FisherG"])
      factory->BookMethod( dataloader, TMVA::Types::kFisher, "FisherG", "H:!V:VarTransform=Gauss" );

   if (Use["BoostedFisher"])
      factory->BookMethod( dataloader, TMVA::Types::kFisher, "BoostedFisher",
                           "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2:!Boost_DetailedMonitoring" );

   if (Use["FDA_MC"])
      factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_MC",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1" );

   if (Use["FDA_GA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
      factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_GA",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=100:Cycles=2:Steps=5:Trim=True:SaveBestGen=1" );

   if (Use["FDA_SA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
      factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_SA",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

   if (Use["FDA_MT"])
      factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_MT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );

   if (Use["FDA_GAMT"])
      factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_GAMT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );

   if (Use["FDA_MCMT"])
      factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_MCMT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20" );

   if (Use["MLP"])
      factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );

   if (Use["MLPBFGS"])
      factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );

   if (Use["MLPBNN"])
      factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=60:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" ); // BFGS training with bayesian regulators

  if (Use["CFMlpANN"])
      factory->BookMethod( dataloader, TMVA::Types::kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N"  ); // n_cycles:#nodes:#nodes:...

   if (Use["TMlpANN"])
      factory->BookMethod( dataloader, TMVA::Types::kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  ); // n_cycles:#nodes:#nodes:...

   if (Use["SVM"])
      factory->BookMethod( dataloader, TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );

   if (Use["RuleFit"])
      factory->BookMethod( dataloader, TMVA::Types::kRuleFit, "RuleFit",
                           "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );



   // Boosted Decision Trees
   if (Use["BDTG"]) // Gradient Boost
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG",
                           "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );

   if (Use["BDT"])  // Adaptive Boost
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

   if (Use["BDTB"]) // Bagging
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTB",
                           "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20" );

   if (Use["BDTD"]) // Decorrelation + Adaptive Boost
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTD",
                           "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate" );

   if (Use["BDTF"])  // Allow Using Fisher discriminant in node splitting for (strong) linearly correlated variables
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTF",
                           "!H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20" );

   if (Use["Grad_nT1000_sh0.05_BSF0.7_MD3"])
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "Grad_nT1000_sh0.05_BSF0.7_MD3",
                           "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.05:UseBaggedBoost:BaggedSampleFraction=0.7:nCuts=20:MaxDepth=3" );
   if (Use["Grad_nT1000_sh0.05_BSF0.7_MD4"])
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "Grad_nT1000_sh0.05_BSF0.7_MD4",
                           "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.05:UseBaggedBoost:BaggedSampleFraction=0.7:nCuts=20:MaxDepth=4" );
   if (Use["Grad_nT1000_sh0.05_BSF0.6_MD3"])
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "Grad_nT1000_sh0.05_BSF0.6_MD3",
                           "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.05:UseBaggedBoost:BaggedSampleFraction=0.6:nCuts=20:MaxDepth=3" );
   if (Use["Grad_nT1000_sh0.05_BSF0.6_MD4"])
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "Grad_nT1000_sh0.05_BSF0.6_MD4",
                           "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.05:UseBaggedBoost:BaggedSampleFraction=0.6:nCuts=20:MaxDepth=4" );
   if (Use["Grad_nT1000_sh0.05_BSF0.4_MD3"])
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "Grad_nT1000_sh0.05_BSF0.4_MD3",
                           "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.05:UseBaggedBoost:BaggedSampleFraction=0.4:nCuts=20:MaxDepth=3" );
   if (Use["Grad_nT1000_sh0.05_BSF0.4_MD4"])
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "Grad_nT1000_sh0.05_BSF0.4_MD4",
                           "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.05:UseBaggedBoost:BaggedSampleFraction=0.4:nCuts=20:MaxDepth=4" );

   if (Use["Grad_nT1000_sh0.03_BSF0.7_MD3"])
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "Grad_nT1000_sh0.03_BSF0.7_MD3",
                           "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.03:UseBaggedBoost:BaggedSampleFraction=0.7:nCuts=20:MaxDepth=3" );
   if (Use["Grad_nT1000_sh0.03_BSF0.7_MD4"])
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "Grad_nT1000_sh0.03_BSF0.7_MD4",
                           "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.03:UseBaggedBoost:BaggedSampleFraction=0.7:nCuts=20:MaxDepth=4" );
   if (Use["Grad_nT1000_sh0.03_BSF0.6_MD3"])
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "Grad_nT1000_sh0.03_BSF0.6_MD3",
                           "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.03:UseBaggedBoost:BaggedSampleFraction=0.6:nCuts=20:MaxDepth=3" );
   if (Use["Grad_nT1000_sh0.03_BSF0.6_MD4"])
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "Grad_nT1000_sh0.03_BSF0.6_MD4",
                           "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.03:UseBaggedBoost:BaggedSampleFraction=0.6:nCuts=20:MaxDepth=4" );
   if (Use["Grad_nT1000_sh0.03_BSF0.4_MD3"])
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "Grad_nT1000_sh0.03_BSF0.4_MD3",
                           "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.03:UseBaggedBoost:BaggedSampleFraction=0.4:nCuts=20:MaxDepth=3" );
   if (Use["Grad_nT1000_sh0.03_BSF0.4_MD4"])
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "Grad_nT1000_sh0.03_BSF0.4_MD4",
                           "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.03:UseBaggedBoost:BaggedSampleFraction=0.4:nCuts=20:MaxDepth=4" );


   if (Use["MLPBFGS_sig"])
      factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLPBFGS_sig", "H:!V:NeuronType=sigmoid:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );

   if (Use["MLPBFGS_hidden_10"])
      factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLPBFGS_hidden10", "H:!V:NeuronType=sigmoid:VarTransform=N:NCycles=600:HiddenLayers=N+10:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );

   if (Use["MLPBFGS_200"])
      factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLPBFGS_200", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=200:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );

   if (Use["MLPBFGS_400"])
      factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLPBFGS_400", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=400:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );

   // For an example of the category classifier usage, see: TMVAVtsCategory
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
   std::cout << "==> TMVAVts is done!" << std::endl;

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
   return TMVAVts(methodList);
}
