#include <iostream>
#include <vector>
#include <TLorentzVector.h>

#include "nano/analysis/interface/topEventSelectionSL.h"
#include "nano/analysis/interface/hadAnalyser.h"
#include "nano/analysis/interface/HadTruthEvents.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

using namespace std;

class testAnalyser : public topEventSelectionSL
{
private:
  
public:
  int nEventsRan, nMuonsSeen, nElectronsSeen, nHadSeen;

  void setOutput(std::string outputName);

 testAnalyser(TTree *tree=0, TTree *had=0, TTree *hadTruth=0, Bool_t isMC = false, Bool_t sle = false, Bool_t slm = false) :
   topEventSelectionSL(tree, had, hadTruth, isMC, sle, slm),
   nEventsRan(0),
   nMuonsSeen(0),
   nElectronsSeen(0),
   nHadSeen(0)
  {
  }
  virtual void Loop();
  ~testAnalyser() {};
};

void testAnalyser::Loop()
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntries();

  // Events loop
  for (Long64_t iev=0; iev < nentries; iev++) {
    GetEntry(iev);

    Reset();
    // int PassedStep =
    EventSelection();

    nMuonsSeen += muonSelection().size();
    nElectronsSeen += elecSelection().size();
    nHadSeen += nhad;

    nEventsRan++;
  }
}

int main(int argc, char* argv[])
{
  std::cout << "Original data file: /xrootd/store/group/nanoAOD/run2_2016v5/SingleMuon/Run2016H-07Aug17-v1/180607_085729/0000/nanoAOD_197.root" << std::endl;
  
  auto inFile = TFile::Open("./data_test_v5.root");
  auto inTree = (TTree*) inFile->Get("Events");
  
  std::cout << " - Testing basic Nano over data file" << std::endl;
  testAnalyser ana(inTree,0,0,false);
  ana.Loop();
  
  std::cout << " - Testing had Nano over data file" << std::endl;
  testAnalyser anahad(inTree,inTree,0,false);
  anahad.Loop();

  std::cout << "Original MC file: /xrootd/store/group/nanoAOD/run2_2016v5/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/180607_115426/0000/nanoAOD_967.root" << std::endl;

  auto inFileMC = TFile::Open("./mc_test_v5.root");
  auto inTreeMC = (TTree*) inFileMC->Get("Events");

  std::cout << " - Testing basic Nano over mc file" << std::endl;
  testAnalyser anaMC(inTreeMC,0,0,false);
  anaMC.Loop();
  
  std::cout << " - Testing had Nano over mc file" << std::endl;
  testAnalyser anahadMC(inTreeMC,inTreeMC,0,false);
  anahadMC.Loop();
  
  std::cout << "Tests finished" << std::endl;


  std::cout << "\n\n:::SUMMARY:::\n" << std::endl;
  std::cout << "   DATA   " << std::endl;
  std::cout << "     - NEVENTS   : " << ana.nEventsRan <<  " ( had " << anahad.nEventsRan << " ) " << std::endl;
  std::cout << "     - NMUONS    : " << ana.nMuonsSeen <<  " ( had " << anahad.nMuonsSeen << " ) " << std::endl;
  std::cout << "     - NELECTRONS: " << ana.nElectronsSeen <<  " ( had " << anahad.nElectronsSeen << " ) " << std::endl;
  std::cout << "     - NHADS     : " << ana.nHadSeen <<  " ( had " << anahad.nHadSeen << " ) " << std::endl;
  std::cout << std::endl;
  std::cout << "   MC   " << std::endl;
  std::cout << "     - NEVENTS   : " << anaMC.nEventsRan <<  " ( had " << anahadMC.nEventsRan << " ) " << std::endl;
  std::cout << "     - NMUONS    : " << anaMC.nMuonsSeen <<  " ( had " << anahadMC.nMuonsSeen << " ) " << std::endl; 
  std::cout << "     - NELECTRONS: " << anaMC.nElectronsSeen <<  " ( had " << anahadMC.nElectronsSeen << " ) " << std::endl; 
  std::cout << "     - NHADS     : " << anaMC.nHadSeen <<  " ( had " << anahadMC.nHadSeen << " ) " << std::endl; 
  std::cout << std::endl;
}
