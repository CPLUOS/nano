// @(#)root/tmva $Id$
/**********************************************************************************
 * Project   : TMVA - a ROOT-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Root Macro: TMVAClassification                                                 *
 *                                                                                *
 * This macro provides examples for the training and testing of the               *
 * TMVA classifiers.                                                              *
 *                                                                                *
 * As input data is used a toy-MC sample consisting of four Gaussian-distributed  *
 * and linearly correlated input variables.                                       *
 *                                                                                *
 * The methods to be used can be switched on and off by means of booleans, or     *
 * via the prompt command, for example:                                           *
 *                                                                                *
 *    root -l ./TMVAClassification.C\(\"Fisher,Likelihood\"\)                     *
 *                                                                                *
 * (note that the backslashes are mandatory)                                      *
 * If no method given, a default set of classifiers is used.                      *
 *                                                                                *
 * The output file "TMVA.root" can be analysed with the use of dedicated          *
 * macros (simply say: root -l <macro.C>), which can be conveniently              *
 * invoked through a GUI that will appear at the end of the run of this macro.    *
 * Launch the GUI via the command:                                                *
 *                                                                                *
 *    root -l ./TMVAGui.C                                                         *
 *                                                                                *
 * You can also compile and run the example with the following commands           *
 *                                                                                *
 *    make                                                                        *
 *    ./TMVAClassification <Methods>                                              *
 *                                                                                *
 * where: <Methods> = "method1 method2"                                           *
 *        are the TMVA classifier names                                           *
 *                                                                                *
 * example:                                                                       *
 *    ./TMVAClassification Fisher LikelihoodPCA BDT                               *
 *                                                                                *
 * If no method given, a default set is of classifiers is used                    *
 **********************************************************************************/

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
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

int TMVAClassification( TString myMethodList = "" )
{
   // The explicit loading of the shared libTMVA is done in TMVAlogon.C, defined in .rootrc
   // if you use your private .rootrc, or run from a different directory, please copy the
   // corresponding lines from .rootrc

   // methods to be processed can be given as an argument; use format:
   //
   // mylinux~> root -l TMVAClassification.C\(\"myMethod1,myMethod2,myMethod3\"\)
   //
   // if you like to use a method via the plugin mechanism, we recommend using
   //
   // mylinux~> root -l TMVAClassification.C\(\"P_myMethod\"\)
   // (an example is given for using the BDT as plugin (see below),
   // but of course the real application is when you write your own
   // method based)

   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // --- Cut optimisation
   Use["Cuts"]            = 0;
   Use["CutsD"]           = 0;
   Use["CutsPCA"]         = 0;
   Use["CutsGA"]          = 0;
   Use["CutsSA"]          = 0;
   // 
   // --- 1-dimensional likelihood ("naive Bayes estimator")
   Use["Likelihood"]      = 0;
   Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
   Use["LikelihoodPCA"]   = 0; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
   Use["LikelihoodKDE"]   = 0;
   Use["LikelihoodMIX"]   = 0;
   //
   // --- Mutidimensional likelihood and Nearest-Neighbour methods
   Use["PDERS"]           = 0;
   Use["PDERSD"]          = 0;
   Use["PDERSPCA"]        = 0;
   Use["PDEFoam"]         = 0;
   Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
   Use["KNN"]             = 0; // k-nearest neighbour method
   //
   // --- Linear Discriminant Analysis
   Use["LD"]              = 0; // Linear Discriminant identical to Fisher
   Use["Fisher"]          = 0;
   Use["FisherG"]         = 0;
   Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
   Use["HMatrix"]         = 0;
   //
   // --- Function Discriminant analysis
   Use["FDA_GA"]          = 0; // minimisation of user-defined function using Genetics Algorithm
   Use["FDA_SA"]          = 0;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   Use["FDA_MCMT"]        = 0;
   //
   // --- Neural Networks (all are feed-forward Multilayer Perceptrons)
   Use["MLP"]             = 0; // Recommended ANN
   Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
   Use["MLPBNN"]          = 0; // Recommended ANN with BFGS training method and bayesian regulator
   Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
   Use["TMlpANN"]         = 0; // ROOT's own ANN
   //
   // --- Support Vector Machine 
   Use["SVM"]             = 0;
   // 
   // --- Boosted Decision Trees
   Use["BDT"]             = 1; // uses Adaptive Boost
   Use["BDTG"]            = 1; // uses Gradient Boost
   Use["BDTB"]            = 0; // uses Bagging
   Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
   Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting 
   // 
   // --- Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
   Use["RuleFit"]         = 0;
   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassification" << std::endl;

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

   // --- Here the preparation phase begins
   // Read training and test data
   // (it is also possible to use ASCII format as input -> see TMVA Users Guide)
   TString fname = "/xrootd/store/user/daniellee/nanoAOD/pre_Run2/results_merged/";
   ///Signal///
   TFile *ttH = TFile::Open( fname+"tth2mu_ttH.root" );
   
   ///Background///
   TFile *TTJets = TFile::Open( fname+"tth2mu_TTJets_DiLept.root" );
   TFile *TTZ = TFile::Open( fname+"tth2mu_TTZToLLNuNu.root" );
   TFile *DYJets = TFile::Open( fname+"tth2mu_DYJets.root" );
   TFile *TTW = TFile::Open( fname+"tth2mu_TTWJetsToLNu.root" );
   TFile *TT = TFile::Open( fname+"tth2mu_SingleTop_tW.root" );
   TFile *TTbar = TFile::Open( fname+"tth2mu_SingleTbar_tW.root" );

   TFile *WWW = TFile::Open( fname+"tth2mu_WWW.root" );
   TFile *WWZ = TFile::Open( fname+"tth2mu_WWZ.root" );
   TFile *WZZ = TFile::Open( fname+"tth2mu_WZZ.root" );
   TFile *ZZZ = TFile::Open( fname+"tth2mu_ZZZ.root" );

   TFile *WW = TFile::Open( fname+"tth2mu_WW.root" );
   TFile *WZTo3LNu = TFile::Open( fname+"tth2mu_WZTo3LNu_amcatnlo.root" );
   TFile *ZZTo2L2Q = TFile::Open( fname+"tth2mu_ZZTo2L2Q.root" );
   TFile *ZZTo4L_powheg = TFile::Open( fname+"tth2mu_ZZTo4L_powheg.root" );
   TFile *WJets = TFile::Open( fname+"tth2mu_WJets.root" );

   TFile *VBF = TFile::Open( fname+"tth2mu_VBF_HToMuMu.root" );
   TFile *GG = TFile::Open( fname+"tth2mu_GG_HToMuMu.root" );
   TFile *WM = TFile::Open( fname+"tth2mu_WMinusH_HToMuMu.root" );
   TFile *WP = TFile::Open( fname+"tth2mu_WPlusH_HToMuMu.root" );
   TFile *ZH = TFile::Open( fname+"tth2mu_ZH_HToMuMu.root" );
   
   // --- Register the training and test trees
   //signal//
   TTree *signal1 = (TTree*)ttH->Get("events");
   TTree *signal2 = (TTree*)VBF->Get("events");
   TTree *signal3 = (TTree*)GG->Get("events");
   TTree *signal4 = (TTree*)WM->Get("events");
   TTree *signal5 = (TTree*)WP->Get("events");
   TTree *signal6 = (TTree*)ZH->Get("events");
   //Background//
   TTree *background0 = (TTree*)TTJets->Get("events");
   TTree *background1 = (TTree*)TTZ->Get("events");
   TTree *background2 = (TTree*)DYJets->Get("events");
   TTree *background3 = (TTree*)TTW->Get("events");
   TTree *background4 = (TTree*)WWW->Get("events");
   TTree *background5 = (TTree*)WWZ->Get("events");
   TTree *background6 = (TTree*)WZZ->Get("events");
   TTree *background7 = (TTree*)ZZZ->Get("events");
   
   TTree *background8 = (TTree*)WW->Get("events");
   TTree *background9 = (TTree*)WZTo3LNu->Get("events");
   TTree *background10 = (TTree*)ZZTo2L2Q->Get("events");
   TTree *background11 = (TTree*)ZZTo4L_powheg->Get("events");
   TTree *background12 = (TTree*)WJets->Get("events");
   
   TTree *background13 = (TTree*)TT->Get("events");
   TTree *background14 = (TTree*)TTbar->Get("events");


   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TString Channel = "FH";
   TString outfileName( "FH.root" );
   if (Channel == "XL")
   {
      TString outfileName( "new_XL.root" );
   }
   else if (Channel == "FH")
   {
      TString outfileName( "new_FH.root" );
   }
   else if (Channel == "Out")
   {
      TString outfileName( "new_Out.root" );
   }
   else if (Channel == "Non")
   {
      TString outfileName( "new_non.root" );
   }
   cout << outfileName << endl;
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
   //TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
   //                                            "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification" );

   TMVA::DataLoader *dataloaderXL=new TMVA::DataLoader("XL");
   TMVA::DataLoader *dataloaderOut=new TMVA::DataLoader("Out");
   TMVA::DataLoader *dataloaderFH=new TMVA::DataLoader("nFH4");
   TMVA::DataLoader *dataloaderNon=new TMVA::DataLoader("nonB");
   // If you wish to modify default settings
   // (please check "src/Config.h" to see all available global options)
   //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
   //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";

   // Define the input variables that shall be used for the MVA training
   // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
   //factory->AddVariable( "myvar1 := var1+var2", 'F' );
   //factory->AddVariable( "myvar2 := var1-var2", "Expression 2", "", 'F' );
   //factory->AddVariable( "var3",                "Variable 3", "units", 'F' );
   //factory->AddVariable( "var4",                "Variable 4", "units", 'F' );
   //VARIABLES FOR Extra Lep CATEGORIES//
   dataloaderXL->AddVariable( "Met", 'F' );
   dataloaderXL->AddVariable( "all_muEtaDiff", 'F' );
   dataloaderXL->AddVariable( "all_muPtDiff", 'F' );
   dataloaderXL->AddVariable( "all_muPhiDiff", 'F' );
   dataloaderXL->AddVariable( "all_muDR", 'F' );
   dataloaderXL->AddVariable( "all_Dilep_Pt", 'F' );
   dataloaderXL->AddVariable( "all_Dilep_Eta", 'F' );
   dataloaderXL->AddVariable( "nelec", 'I' );
   dataloaderXL->AddVariable( "nmuon", 'I' );
   dataloaderXL->AddVariable( "nnonbjet", 'I' );
   dataloaderXL->AddVariable( "nbjet", 'I' );
   dataloaderXL->AddVariable( "Central_Jets", 'I' );
   dataloaderXL->AddVariable( "Forward_Jets", 'I' );
   dataloaderXL->AddVariable( "minDR", 'F' );
   dataloaderXL->AddVariable( "XlepPt", 'F' );
   dataloaderXL->AddVariable( "mT2", 'F' );
   dataloaderXL->AddVariable( "mT", 'F' );
   dataloaderXL->AddVariable( "DiJetM12", 'F' );
   //dataloaderXL->AddVariable( "CSV", 'F' );
   
   //VARIABLES FOR nFH4 CATEGORIES//
   dataloaderFH->AddVariable( "Met", 'F' );
   //dataloaderFH->AddVariable( "Met_phi", 'F' );
   dataloaderFH->AddVariable( "all_muEtaDiff", 'F' );
   dataloaderFH->AddVariable( "all_muPtDiff", 'F' );
   dataloaderFH->AddVariable( "all_muPhiDiff", 'F' );
   dataloaderFH->AddVariable( "all_muDR", 'F' );
   dataloaderFH->AddVariable( "all_Dilep_Pt", 'F' );
   dataloaderFH->AddVariable( "all_Dilep_Eta", 'F' );
   dataloaderFH->AddVariable( "nnonbjet", 'I' );
   dataloaderFH->AddVariable( "nbjet", 'I' );
   dataloaderFH->AddVariable( "Central_Jets", 'I' );
   dataloaderFH->AddVariable( "Forward_Jets", 'I' );
   dataloaderFH->AddVariable( "minDR1", 'F' );
   dataloaderFH->AddVariable( "minDR2", 'F' );
   dataloaderFH->AddVariable( "mT2", 'F' );
   dataloaderFH->AddVariable( "mT", 'F' );
   dataloaderFH->AddVariable( "DiJetM12", 'F' );
   dataloaderFH->AddVariable( "DiJetM13", 'F' );
   dataloaderFH->AddVariable( "DiJetM14", 'F' );
   dataloaderFH->AddVariable( "DiJetM23", 'F' );
   dataloaderFH->AddVariable( "DiJetM24", 'F' );
   dataloaderFH->AddVariable( "DiJetM34", 'F' );
   
   //VARIABLES FOR Out CATEGORIES//
   dataloaderOut->AddVariable( "Met", 'F' );
   //dataloaderOut->AddVariable( "Met_phi", 'F' );
   dataloaderOut->AddVariable( "all_muEtaDiff", 'F' );
   dataloaderOut->AddVariable( "all_muPtDiff", 'F' );
   dataloaderOut->AddVariable( "all_muPhiDiff", 'F' );
   dataloaderOut->AddVariable( "all_muDR", 'F' );
   dataloaderOut->AddVariable( "all_Dilep_Pt", 'F' );
   dataloaderOut->AddVariable( "all_Dilep_Eta", 'F' );
   dataloaderOut->AddVariable( "nnonbjet", 'I' );
   dataloaderOut->AddVariable( "nbjet", 'I' );
   dataloaderOut->AddVariable( "Central_Jets", 'I' );
   dataloaderOut->AddVariable( "Forward_Jets", 'I' );
   dataloaderOut->AddVariable( "mT2", 'F' );
   dataloaderOut->AddVariable( "mT", 'F' );
   dataloaderOut->AddVariable( "DiJetM12", 'F' );
   
   //VARIABLES FOR non CATEGORIES//
   dataloaderNon->AddVariable( "Met", 'F' );
   dataloaderNon->AddVariable( "all_muEtaDiff", 'F' );
   dataloaderNon->AddVariable( "all_muPtDiff", 'F' );
   dataloaderNon->AddVariable( "all_muPhiDiff", 'F' );
   dataloaderNon->AddVariable( "all_muDR", 'F' );
   dataloaderNon->AddVariable( "all_Dilep_Pt", 'F' );
   dataloaderNon->AddVariable( "all_Dilep_Eta", 'F' );
   dataloaderNon->AddVariable( "nelec", 'I' );
   dataloaderNon->AddVariable( "nexLep", 'I' );
   dataloaderNon->AddVariable( "nmuon", 'I' );
   dataloaderNon->AddVariable( "Central_Jets", 'I' );
   dataloaderNon->AddVariable( "Forward_Jets", 'I' );
   dataloaderNon->AddVariable( "etaJ1", 'F' );
   dataloaderNon->AddVariable( "etaJ2", 'F' );
   dataloaderNon->AddVariable( "DijetM1", 'F' );
   dataloaderNon->AddVariable( "DijetM2", 'F' );
   dataloaderNon->AddVariable( "DijetEta1", 'F' );
   dataloaderNon->AddVariable( "DijetEta2", 'F' );

   // You can add so-called "Spectator variables", which are not used in the MVA training,
   // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
   // input variables, the response values of all trained MVAs, and the spectator variables
   //factory->AddSpectator( "spec1 := var1*2",  "Spectator 1", "units", 'F' );
   //factory->AddSpectator( "spec2 := var1*3",  "Spectator 2", "units", 'F' );

   
   // global event weights per tree (see below for setting event-wise weights)
   
   Double_t Lumi = 35900;
   TH1D *Sig1 = (TH1D*)ttH->Get("weight");
   TH1D *Sig2 = (TH1D*)VBF->Get("weight");
   TH1D *Sig3 = (TH1D*)GG->Get("weight");
   TH1D *Sig4 = (TH1D*)WM->Get("weight");
   TH1D *Sig5 = (TH1D*)WP->Get("weight");
   TH1D *Sig6 = (TH1D*)ZH->Get("weight");
   TH1D *b0 = (TH1D*)TTJets->Get("weight");
   TH1D *b1 = (TH1D*)TTZ->Get("weight");
   TH1D *b2 = (TH1D*)DYJets->Get("weight");
   TH1D *b3 = (TH1D*)TTW->Get("weight");
   TH1D *b4 = (TH1D*)WWW->Get("weight");
   TH1D *b5 = (TH1D*)WWZ->Get("weight");
   TH1D *b6 = (TH1D*)WZZ->Get("weight");
   TH1D *b7 = (TH1D*)ZZZ->Get("weight");
   TH1D *b8 = (TH1D*)WW->Get("weight");
   TH1D *b9 = (TH1D*)WZTo3LNu->Get("weight");
   TH1D *b10 = (TH1D*)ZZTo2L2Q->Get("weight");
   TH1D *b11 = (TH1D*)ZZTo4L_powheg->Get("weight");
   TH1D *b12 = (TH1D*)WJets->Get("weight");
   TH1D *b13 = (TH1D*)TT->Get("weight");
   TH1D *b14 = (TH1D*)TTbar->Get("weight");

   Double_t signalWeight = (Lumi*0.0001113615)/Sig1->Integral(0,1);
   Double_t sig2 = (Lumi*0.0008208)/Sig2->Integral(0,1);
   Double_t sig3 = (Lumi*0.009618)/Sig3->Integral(0,1);
   Double_t sig4 = (Lumi*0.0001164)/Sig4->Integral(0,1);
   Double_t sig5 = (Lumi*0.0001858)/Sig5->Integral(0,1);
   Double_t sig6 = (Lumi*0.0002136)/Sig6->Integral(0,1);
   Double_t b0a = (Lumi*85.65)/b0->Integral(0,1);
   Double_t b1a = (Lumi*0.2529)/b1->Integral(0,1);
   Double_t b2a = (Lumi*5765.4)/b2->Integral(0,1);
   Double_t b3a = (Lumi*0.2043)/b3->Integral(0,1);
   Double_t b4a = (Lumi*0.2086)/b4->Integral(0,1);
   Double_t b5a = (Lumi*0.1651)/b5->Integral(0,1);
   Double_t b6a = (Lumi*0.05565)/b6->Integral(0,1);
   Double_t b7a = (Lumi*0.01398)/b7->Integral(0,1);
   Double_t b8a = (Lumi*16.523)/b8->Integral(0,1);
   Double_t b9a = (Lumi*5.26)/b9->Integral(0,1);
   Double_t b10a = (Lumi*3.22)/b10->Integral(0,1);
   Double_t b11a = (Lumi*1.212)/b11->Integral(0,1);
   Double_t b12a = (Lumi*0.2043)/b12->Integral(0,1);
   Double_t b13a = (Lumi*35.85)/b13->Integral(0,1);
   Double_t b14a = (Lumi*35.85)/b14->Integral(0,1);
   
   // You can add an arbitrary number of signal or background trees
   if (Channel == "XL")
   {  
      dataloaderXL->AddSignalTree( signal1,     signalWeight     );
      dataloaderXL->AddSignalTree( signal2 , sig2 );
      dataloaderXL->AddSignalTree( signal3 , sig3 );
      dataloaderXL->AddSignalTree( signal4 , sig4 );
      dataloaderXL->AddSignalTree( signal5 , sig5 );
      dataloaderXL->AddSignalTree( signal6 , sig6 );
      dataloaderXL->AddBackgroundTree( background0, b0a );
      dataloaderXL->AddBackgroundTree( background1, b1a);
      dataloaderXL->AddBackgroundTree( background2, b2a );
      dataloaderXL->AddBackgroundTree( background3, b3a );
      dataloaderXL->AddBackgroundTree( background4, b4a );
      dataloaderXL->AddBackgroundTree( background5, b5a );
      dataloaderXL->AddBackgroundTree( background6, b6a );
      dataloaderXL->AddBackgroundTree( background7, b7a );
      dataloaderXL->AddBackgroundTree( background8, b8a );
      dataloaderXL->AddBackgroundTree( background9, b9a );
      dataloaderXL->AddBackgroundTree( background10, b10a );
      dataloaderXL->AddBackgroundTree( background11, b11a );
      dataloaderXL->AddBackgroundTree( background12, b12a );
      dataloaderXL->AddBackgroundTree( background13, b13a );
      dataloaderXL->AddBackgroundTree( background14, b14a );
   }
   else if (Channel == "FH")
   {
      dataloaderFH->AddSignalTree( signal1,     signalWeight     );
      dataloaderFH->AddSignalTree( signal2 , sig2 );
      dataloaderFH->AddSignalTree( signal3 , sig3 );
      dataloaderFH->AddSignalTree( signal4 , sig4 );
      dataloaderFH->AddSignalTree( signal5 , sig5 );
      dataloaderFH->AddSignalTree( signal6 , sig6 );
      dataloaderFH->AddBackgroundTree( background0, b0a );
      dataloaderFH->AddBackgroundTree( background1, b1a);
      dataloaderFH->AddBackgroundTree( background2, b2a );
      dataloaderFH->AddBackgroundTree( background3, b3a );
      dataloaderFH->AddBackgroundTree( background4, b4a );
      dataloaderFH->AddBackgroundTree( background5, b5a );
      dataloaderFH->AddBackgroundTree( background6, b6a );
      dataloaderFH->AddBackgroundTree( background7, b7a );
      dataloaderFH->AddBackgroundTree( background8, b8a );
      dataloaderFH->AddBackgroundTree( background9, b9a );
      dataloaderFH->AddBackgroundTree( background10, b10a );
      dataloaderFH->AddBackgroundTree( background11, b11a );
      dataloaderFH->AddBackgroundTree( background12, b12a );
      dataloaderFH->AddBackgroundTree( background13, b13a );
      dataloaderFH->AddBackgroundTree( background14, b14a );
   }
   else if (Channel == "Out")
   {
      dataloaderOut->AddSignalTree( signal1,     signalWeight     );
      dataloaderOut->AddSignalTree( signal2 , sig2 );
      dataloaderOut->AddSignalTree( signal3 , sig3 );
      dataloaderOut->AddSignalTree( signal4 , sig4 );
      dataloaderOut->AddSignalTree( signal5 , sig5 );
      dataloaderOut->AddSignalTree( signal6 , sig6 );
      dataloaderOut->AddBackgroundTree( background0, b0a );
      dataloaderOut->AddBackgroundTree( background1, b1a);
      dataloaderOut->AddBackgroundTree( background2, b2a );
      dataloaderOut->AddBackgroundTree( background3, b3a );
      dataloaderOut->AddBackgroundTree( background4, b4a );
      dataloaderOut->AddBackgroundTree( background5, b5a );
      dataloaderOut->AddBackgroundTree( background6, b6a );
      dataloaderOut->AddBackgroundTree( background7, b7a );
      dataloaderOut->AddBackgroundTree( background8, b8a );
      dataloaderOut->AddBackgroundTree( background9, b9a );
      dataloaderOut->AddBackgroundTree( background10, b10a );
      dataloaderOut->AddBackgroundTree( background11, b11a );
      dataloaderOut->AddBackgroundTree( background12, b12a );
      dataloaderOut->AddBackgroundTree( background13, b13a );
      dataloaderOut->AddBackgroundTree( background14, b14a );
   }
   else if (Channel == "Non")
   {
      dataloaderNon->AddSignalTree( signal1,     signalWeight     );
      dataloaderNon->AddSignalTree( signal2 , sig2 );
      dataloaderNon->AddSignalTree( signal3 , sig3 );
      dataloaderNon->AddSignalTree( signal4 , sig4 );
      dataloaderNon->AddSignalTree( signal5 , sig5 );
      dataloaderNon->AddSignalTree( signal6 , sig6 );
      dataloaderNon->AddBackgroundTree( background0, b0a );
      dataloaderNon->AddBackgroundTree( background1, b1a);
      dataloaderNon->AddBackgroundTree( background2, b2a );
      dataloaderNon->AddBackgroundTree( background3, b3a );
      dataloaderNon->AddBackgroundTree( background4, b4a );
      dataloaderNon->AddBackgroundTree( background5, b5a );
      dataloaderNon->AddBackgroundTree( background6, b6a );
      dataloaderNon->AddBackgroundTree( background7, b7a );
      dataloaderNon->AddBackgroundTree( background8, b8a );
      dataloaderNon->AddBackgroundTree( background9, b9a );
      dataloaderNon->AddBackgroundTree( background10, b10a );
      dataloaderNon->AddBackgroundTree( background11, b11a );
      dataloaderNon->AddBackgroundTree( background12, b12a );
      dataloaderNon->AddBackgroundTree( background13, b13a );
      dataloaderNon->AddBackgroundTree( background14, b14a );
   } 
   //dataloaderNon->AddSignalTree    ( signal,     signalWeight     );
   //dataloaderNon->AddBackgroundTree( background0, b00 );
   //dataloaderNon->AddBackgroundTree( background1, b11);
   //dataloaderNon->AddBackgroundTree( background2, b22 );
   //dataloaderNon->AddBackgroundTree( background3, b33 );
   //dataloaderNon->AddSignalTree    ( signal,     signalWeight     );
   //dataloaderNon->AddBackgroundTree( background0, backgroundWeight );
   //dataloaderNon->AddBackgroundTree( background1, backgroundWeight );
   
   // To give different trees for training and testing, do as follows:
   //    factory->AddSignalTree( signalTrainingTree, signalTrainWeight, "Training" );
   //    factory->AddSignalTree( signalTestTree,     signalTestWeight,  "Test" );
   
   // Use the following code instead of the above two or four lines to add signal and background
   // training and test events "by hand"
   // NOTE that in this case one should not give expressions (such as "var1+var2") in the input
   //      variable definition, but simply compute the expression before adding the event
   //
   //     // --- begin ----------------------------------------------------------
   //     std::vector<Double_t> vars( 4 ); // vector has size of number of input variables
   //     Float_t  treevars[4], weight;
   //     
   //     // Signal
   //     for (UInt_t ivar=0; ivar<4; ivar++) signal->SetBranchAddress( Form( "var%i", ivar+1 ), &(treevars[ivar]) );
   //     for (UInt_t i=0; i<signal->GetEntries(); i++) {
   //        signal->GetEntry(i);
   //        for (UInt_t ivar=0; ivar<4; ivar++) vars[ivar] = treevars[ivar];
   //        // add training and test events; here: first half is training, second is testing
   //        // note that the weight can also be event-wise
   //        if (i < signal->GetEntries()/2.0) factory->AddSignalTrainingEvent( vars, signalWeight );
   //        else                              factory->AddSignalTestEvent    ( vars, signalWeight );
   //     }
   //   
   //     // Background (has event weights)
   //     background->SetBranchAddress( "weight", &weight );
   //     for (UInt_t ivar=0; ivar<4; ivar++) background->SetBranchAddress( Form( "var%i", ivar+1 ), &(treevars[ivar]) );
   //     for (UInt_t i=0; i<background->GetEntries(); i++) {
   //        background->GetEntry(i);
   //        for (UInt_t ivar=0; ivar<4; ivar++) vars[ivar] = treevars[ivar];
   //        // add training and test events; here: first half is training, second is testing
   //        // note that the weight can also be event-wise
   //        if (i < background->GetEntries()/2) factory->AddBackgroundTrainingEvent( vars, backgroundWeight*weight );
   //        else                                factory->AddBackgroundTestEvent    ( vars, backgroundWeight*weight );
   //     }
         // --- end ------------------------------------------------------------
   //
   // --- end of tree registration 

   // Set individual event weights (the variables must exist in the original TTree)
   //    for signal    : factory->SetSignalWeightExpression    ("weight1*weight2");
   //    for background: factory->SetBackgroundWeightExpression("weight1*weight2");
   //factory->SetBackgroundWeightExpression( "weight" );                 /////////////////////////////////<~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   //Apply additional cuts on the signal and background samples (can be different)
   //TCut mycuts = "FL == 1 && channel == 1 && Dilep.M() > 12 && FL_lep2Pt > 0 && FL_lep2Pt < 250 && Event_No % 2 != 0"; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
   //TCut mycutb = "FL == 1 && channel == 1 && Dilep.M() > 12 && FL_lep2Pt > 0 && FL_lep2Pt < 250"; // for example: TCut mycutb = "abs(var1)<0.5";
    TString cuts = "";
    TString cutb = "";
    
    if (Channel == "XL")
    {
       //cuts = "XL == 1 && XlepPt < 350 && Event_No % 2 != 0 && mT2 < 900"; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
       //cutb = "XL == 1 && XlepPt < 350 && mT2 < 900 "; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
       //cuts = "XL == 1 && Dilep.M() > 110 && Dilep.M() < 160 && Event_No % 2 != 0"; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
       cuts = "XL == 1 && minDR < 60000 && Dilep.M() > 110 && Dilep.M() < 160 && Event_No % 2 != 0"; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
      // cuts = "XL == 1 && minDR < 60000 && Dilep.M() > 110 && Dilep.M() < 160"; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
       cutb = "XL == 1 && minDR < 60000 && Dilep.M() > 110 && Dilep.M() < 160 "; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
    }
    else if (Channel == "FH")
    {
       cout << "Cut loaded" << endl;
       //cuts = "nFH4 == 1 && XlepPt < 350 && XmT2 < 350000 && HminDR2 > 0 && HminDR1 > 0 && HminDR1 < 4 && HminDR2 < 4 && Event_No % 2 != 0"; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
       //cutb = "nFH4 == 1 && XlepPt < 350 && XmT2 < 350000 && HminDR2 > 0 && HminDR1 > 0 && HminDR1 < 4 && HminDR2 < 4" ; // for example: TCut mycutb = "abs(var1)<0.5";
       //cuts = "nFH4 == 1 && Dilep.M() > 110 && Dilep.M() < 160 && Event_No % 2 != 0"; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
       cuts = "nFH4 == 1 && Dilep.M() > 110 && Dilep.M() < 160"; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
       cutb = "nFH4 == 1 && Dilep.M() > 110 && Dilep.M() < 160" ; // for example: TCut mycutb = "abs(var1)<0.5";
    }
    else if (Channel == "Out")
    {
       cuts = "Out == 1 && Dilep.M() > 110 && Dilep.M() < 160 && Event_No % 2 != 0"; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
       cutb = "Out == 1 && Dilep.M() > 110 && Dilep.M() < 160 "; // for example: TCut mycutb = "abs(var1)<0.5";
    }
    else if (Channel == "Non")
    {
       cuts = "nonB == 1 && Dilep.M() > 110 && Dilep.M() < 160 && Event_No % 2 != 0"; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
       cutb = "nonB == 1 && Dilep.M() > 110 && Dilep.M() < 160 "; // for example: TCut mycutb = "abs(var1)<0.5";
    }
    TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
    TCut mycutb = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
    mycuts = cuts; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
    mycutb = cutb; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
    cout << mycuts <<endl;
    cout << mycutb <<endl;

   // Tell the factory how to use the training and testing events
   //
   // If no numbers of events are given, half of the events in the tree are used 
   // for training, and the other half for testing:
   //    factory->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
   // To also specify the number of testing events, use:
   //    factory->PrepareTrainingAndTestTree( mycut,
   //                                         "NSigTrain=3000:NBkgTrain=3000:NSigTest=3000:NBkgTest=3000:SplitMode=Random:!V" );
   Int_t ntrain = 0;
   Int_t ntest = 0;
   //dataloaderXL->PrepareTrainingAndTestTree( mycuts, mycutb,
   //                                     "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
   //dataloaderFH->PrepareTrainingAndTestTree( mycuts, mycutb,
   ///                                     "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
   if (Channel == "XL")
   {
   dataloaderXL->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
   //                                     "nTrain_Signal=25456:nTrain_Background=6073:nTest_Signal=10910:nTest_Background=2603:SplitMode=Random:NormMode=NumEvents:!V" );
   }
   else if (Channel == "FH")
   {
   dataloaderFH->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
   //                                     "nTrain_Signal=42395:nTrain_Background=14315:nTest_Signal=18169:nTest_Background=6135:SplitMode=Random:NormMode=NumEvents:!V" );
   }
   else if (Channel == "Out")
   {
   dataloaderOut->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
     //                                   ("nTrain_Signal="+to_string(ntrain)+":nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V").c_str() );
   //                                     "nTrain_Signal=53771:nTrain_Background=67713:nTest_Signal=23045:nTest_Background=29020:SplitMode=Random:NormMode=NumEvents:!V" );
   }
   else if (Channel == "Non")
   {
   dataloaderNon->PrepareTrainingAndTestTree( mycuts, mycutb,
   //                                     "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
                                        "nTrain_Signal=399391:nTrain_Background=301148:nTest_Signal=171167:nTest_Background=129064:SplitMode=Random:NormMode=NumEvents:!V" );
   }

   // ---- Book MVA methods
   //
   // Please lookup the various method configuration options in the corresponding cxx files, eg:
   // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
   // it is possible to preset ranges in the option string in which the cut optimisation should be done:
   // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

   //// Cut optimisation
   //if (Use["Cuts"])
   //   factory->BookMethod( dataloaderXL,TMVA::Types::kCuts, "Cuts",
   //                        "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );

   //if (Use["CutsD"])
   //   factory->BookMethod( dataloader,TMVA::Types::kCuts, "CutsD",
   //                        "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" );

   //if (Use["CutsPCA"])
   //   factory->BookMethod( dataloader,TMVA::Types::kCuts, "CutsPCA",
   //                        "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA" );

   //if (Use["CutsGA"])
   //   factory->BookMethod( dataloader,TMVA::Types::kCuts, "CutsGA",
   //                        "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" );

   //if (Use["CutsSA"])
   //   factory->BookMethod( dataloader,TMVA::Types::kCuts, "CutsSA",
   //                        "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

   //// Likelihood ("naive Bayes estimator")
   //if (Use["Likelihood"])
   //   factory->BookMethod( dataloader,TMVA::Types::kLikelihood, "Likelihood",
   //                        "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" );

   //// Decorrelated likelihood
   //if (Use["LikelihoodD"])
   //   factory->BookMethod( dataloader,TMVA::Types::kLikelihood, "LikelihoodD",
   //                        "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" );

   //// PCA-transformed likelihood
   //if (Use["LikelihoodPCA"])
   //   factory->BookMethod( dataloader,TMVA::Types::kLikelihood, "LikelihoodPCA",
   //                        "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA" ); 

   //// Use a kernel density estimator to approximate the PDFs
   //if (Use["LikelihoodKDE"])
   //   factory->BookMethod( dataloader,TMVA::Types::kLikelihood, "LikelihoodKDE",
   //                        "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50" ); 

   //// Use a variable-dependent mix of splines and kernel density estimator
   //if (Use["LikelihoodMIX"])
   //   factory->BookMethod( dataloader,TMVA::Types::kLikelihood, "LikelihoodMIX",
   //                        "!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50" ); 

   //// Test the multi-dimensional probability density estimator
   //// here are the options strings for the MinMax and RMS methods, respectively:
   ////      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );
   ////      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );
   //if (Use["PDERS"])
   //   factory->BookMethod( dataloader,TMVA::Types::kPDERS, "PDERS",
   //                        "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" );

   //if (Use["PDERSD"])
   //   factory->BookMethod( dataloader,TMVA::Types::kPDERS, "PDERSD",
   //                        "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate" );

   //if (Use["PDERSPCA"])
   //   factory->BookMethod( dataloader,TMVA::Types::kPDERS, "PDERSPCA",
   //                        "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA" );

   //// Multi-dimensional likelihood estimator using self-adapting phase-space binning
   //if (Use["PDEFoam"])
   //   factory->BookMethod( dataloader,TMVA::Types::kPDEFoam, "PDEFoam",
   //                        "!H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" );

   //if (Use["PDEFoamBoost"])
   //   factory->BookMethod( dataloader,TMVA::Types::kPDEFoam, "PDEFoamBoost",
   //                        "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T" );

   //// K-Nearest Neighbour classifier (KNN)
   //if (Use["KNN"])
   //   factory->BookMethod( dataloader,TMVA::Types::kKNN, "KNN",
   //                        "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );

   //// H-Matrix (chi2-squared) method
   //if (Use["HMatrix"])
   //   factory->BookMethod( dataloader,TMVA::Types::kHMatrix, "HMatrix", "!H:!V:VarTransform=None" );

   //// Linear discriminant (same as Fisher discriminant)
   //if (Use["LD"])
   //   factory->BookMethod( dataloader,TMVA::Types::kLD, "LD", "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

   //// Fisher discriminant (same as LD)
   //if (Use["Fisher"])
   //   factory->BookMethod( dataloader,TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

   //// Fisher with Gauss-transformed input variables
   //if (Use["FisherG"])
   //   factory->BookMethod( dataloader,TMVA::Types::kFisher, "FisherG", "H:!V:VarTransform=Gauss" );

   //// Composite classifier: ensemble (tree) of boosted Fisher classifiers
   //if (Use["BoostedFisher"])
   //   factory->BookMethod( dataloader,TMVA::Types::kFisher, "BoostedFisher", 
   //                        "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2:!Boost_DetailedMonitoring" );

   //// Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
   //if (Use["FDA_MC"])
   //   factory->BookMethod( dataloader,TMVA::Types::kFDA, "FDA_MC",
   //                        "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1" );

   //if (Use["FDA_GA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
   //   factory->BookMethod( dataloader,TMVA::Types::kFDA, "FDA_GA",
   //                        "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1" );

   //if (Use["FDA_SA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
   //   factory->BookMethod( dataloader,TMVA::Types::kFDA, "FDA_SA",
   //                        "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

   //if (Use["FDA_MT"])
   //   factory->BookMethod( dataloader,TMVA::Types::kFDA, "FDA_MT",
   //                        "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );

   //if (Use["FDA_GAMT"])
   //   factory->BookMethod( dataloader,TMVA::Types::kFDA, "FDA_GAMT",
   //                        "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );

   //if (Use["FDA_MCMT"])
   //   factory->BookMethod( dataloader,TMVA::Types::kFDA, "FDA_MCMT",
   //                        "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20" );

   //// TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
   //if (Use["MLP"])
   //   factory->BookMethod( dataloader,TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );

   //if (Use["MLPBFGS"])
   //   factory->BookMethod( dataloader,TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );

   //if (Use["MLPBNN"])
   //   factory->BookMethod( dataloader,TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" ); // BFGS training with bayesian regulators

   //// CF(Clermont-Ferrand)ANN
   //if (Use["CFMlpANN"])
   //   factory->BookMethod( dataloader,TMVA::Types::kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=2000:HiddenLayers=N+1,N"  ); // n_cycles:#nodes:#nodes:...  

   //// Tmlp(Root)ANN
   //if (Use["TMlpANN"])
   //   factory->BookMethod( dataloader,TMVA::Types::kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  ); // n_cycles:#nodes:#nodes:...

   //// Support Vector Machine
   //if (Use["SVM"])
   //   factory->BookMethod( dataloader,TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );

   //// Boosted Decision Trees
   if (Use["BDTG"]) // Gradient Boost
   {
      if (Channel == "XL")
      {
         factory->BookMethod( dataloaderXL,TMVA::Types::kBDT, "BDTG",
                              "!H:!V:NTrees=300:BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=5:NegWeightTreatment=Pray" );
      }
      else if (Channel == "FH")
      {
         factory->BookMethod( dataloaderFH,TMVA::Types::kBDT, "BDTG",
                              "!H:!V:NTrees=300:BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=5:NegWeightTreatment=Pray" );
      }
      else if (Channel == "Out")
      {
         factory->BookMethod( dataloaderOut,TMVA::Types::kBDT, "BDTG",
                              "!H:!V:NTrees=500:BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=5:NegWeightTreatment=Pray" );
      }
      else if (Channel == "Non")
      {
         factory->BookMethod( dataloaderNon,TMVA::Types::kBDT, "BDTG",
                              "!H:!V:NTrees=500:BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=5:NegWeightTreatment=Pray" );
      }
   }
    //  factory->BookMethod( dataloaderNon,TMVA::Types::kBDT, "BDTG",
    //                       "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );
                          // "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.08:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=30:MaxDepth=4" );

   //if (Use["BDT"])  // Adaptive Boost
   //   factory->BookMethod( dataloaderNon,TMVA::Types::kBDT, "BDT",
   //                        "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );
   //if (Use["BDT"])  // Adaptive Boost
  // {
   //   if (Channel == "XL")
   //   {
   //   factory->BookMethod( dataloaderXL,TMVA::Types::kBDT, "BDT",
   //                        "!H:!V:NTrees=300:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );
   //   }
   //  else if (Channel == "FH")
   //   {
   //   factory->BookMethod( dataloaderFH,TMVA::Types::kBDT, "BDT",
   //                        "!H:!V:NTrees=300:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );
   //   }
   //   else if (Channel == "Out")
   //   {
   //   factory->BookMethod( dataloaderOut,TMVA::Types::kBDT, "BDT",
   //                        "!H:!V:NTrees=450:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );
   //   }
   //   else if (Channel == "Non")
   //   {
   //   factory->BookMethod( dataloaderNon,TMVA::Types::kBDT, "BDT",
   //                        "!H:!V:NTrees=450:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );
   //   }
   //}
   //if (Use["BDTB"]) // Bagging
   //   factory->BookMethod( dataloaderXL,TMVA::Types::kBDT, "BDTB",
   //                        "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20" );

   //if (Use["BDTD"]) // Decorrelation + Adaptive Boost
   //   factory->BookMethod( dataloader,TMVA::Types::kBDT, "BDTD",
   //                        "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate" );

   //if (Use["BDTF"])  // Allow Using Fisher discriminant in node splitting for (strong) linearly correlated variables
   //   factory->BookMethod( dataloader,TMVA::Types::kBDT, "BDTMitFisher",
   //                        "!H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20" );

   //// RuleFit -- TMVA implementation of Friedman's method
   //if (Use["RuleFit"])
   //   factory->BookMethod( dataloader,TMVA::Types::kRuleFit, "RuleFit",
   //                        "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );

   // For an example of the category classifier usage, see: TMVAClassificationCategory

   // --------------------------------------------------------------------------------------------------

   // ---- Now you can optimize the setting (configuration) of the MVAs using the set of training events

   // ---- STILL EXPERIMENTAL and only implemented for BDT's ! 
   // factory->OptimizeAllMethods("SigEffAt001","Scan");
   // factory->OptimizeAllMethods("ROCIntegral","FitGA");

   // --------------------------------------------------------------------------------------------------

   // ---- Now you can tell the factory to train, test, and evaluate the MVAs

   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

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
   return TMVAClassification(methodList); 
}
