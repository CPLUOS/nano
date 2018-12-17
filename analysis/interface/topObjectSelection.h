#ifndef topObjectSelection_H
#define topObjectSelection_H

#include "nanoBase.h"
#include "nano/external/interface/TopTriggerSF.h"
//#include "nano/external/interface/TTbarModeDefs.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include <TRandom3.h>

class topObjectSelection : public nanoBase
{
protected:
  std::vector<Float_t> b_csvweights;
  float b_btagweight, b_btagweight_up, b_btagweight_dn;
  Float_t b_isolep;
  
  Float_t b_maxBDiscr_nonb;
  
  UInt_t m_unFlag;
  
  JetCorrectionUncertainty *jecUnc;
  
  JME::JetResolution jetResObj;
  JME::JetResolutionScaleFactor jetResSFObj;
  TRandom3 *rndEngine;
  
public:
  // YOU MUST SET UP ALL IN THE BELOW!!!
  // (SetCutValues() will force you to do it)
  Float_t cut_ElectronPt;
  Float_t cut_ElectronEta;
  Int_t  *cut_ElectronIDType; // For example, cut_ElectronIDType = Electron_cutBased;
  Int_t   cut_ElectronIDCut;
  Float_t cut_ElectronSCEtaLower;
  Float_t cut_ElectronSCEtaUpper;
  Float_t cut_ElectronRelIso03All;
  
  Bool_t *cut_MuonIDType; // For example, cut_MuonIDType = Muon_tightId;
  Float_t cut_MuonPt;
  Float_t cut_MuonEta;
  Float_t cut_MuonRelIso04All;
  
  Float_t cut_VetoElectronPt;
  Float_t cut_VetoElectronEta;
  Int_t  *cut_VetoElectronIDType; // For example, cut_VetoElectronIDType = Electron_cutBased;
  Int_t   cut_VetoElectronIDCut;
  Float_t cut_VetoElectronSCEtaLower;
  Float_t cut_VetoElectronSCEtaUpper;
  Float_t cut_VetoElectronRelIso03All;
  
  Bool_t *cut_VetoMuonIDType; // For example, cut_MuonIDType = NULL; or cut_MuonIDType = Muon_looseId;
  Float_t cut_VetoMuonPt;
  Float_t cut_VetoMuonEta;
  Float_t cut_VetoMuonRelIso04All;
  
  Float_t cut_GenJetPt;
  Float_t cut_GenJetEta;
  Float_t cut_GenJetConeSizeOverlap;
  
  Int_t   cut_JetID;
  Float_t cut_JetPt;
  Float_t cut_JetEta;
  Float_t cut_JetConeSizeOverlap;
  
  Int_t    cut_BJetID;
  Float_t  cut_BJetPt;
  Float_t  cut_BJetEta;
  Float_t  cut_BJetConeSizeOverlap;
  Float_t *cut_BJetTypeBTag; // For example, set it as cut_BJetTypeBTag = Jet_btagCSVV2;
  Float_t  cut_BJetBTagCut;
  
  Float_t m_fDRcone_JER;
  Float_t m_fResFactorMathcer;

public: 
  // Tip: If you want to use your own additional cut with the existing cut, 
  // instead of copying the existing code, use the following: 
  // bool [YOUR CLASS NAME]::additionalConditionFor[...]() {
  //   if ( !topObjectSelection::additionalConditionFor[...]() ) return false;
  //   [YOUR OWN CONDITIONS...]
  // }
  virtual bool additionalConditionForElectron(UInt_t nIdx) {return true;};
  virtual bool additionalConditionForMuon(UInt_t nIdx) {return Muon_isPFcand[ nIdx ] && Muon_globalMu[ nIdx ];};
  virtual bool additionalConditionForVetoElectron(UInt_t nIdx) {return true;};
  virtual bool additionalConditionForVetoMuon(UInt_t nIdx) {
    return Muon_isPFcand[ nIdx ] && ( Muon_globalMu[ nIdx ] || Muon_trackerMu[ nIdx ] );
  };
  virtual bool additionalConditionForGenJet(UInt_t nIdx) {return true;};
  
  virtual bool additionalConditionForJet(UInt_t nIdx, Float_t &fJetPt, Float_t &fJetEta, Float_t &fJetPhi, Float_t &fJetMass) 
    {return true;};
  virtual bool additionalConditionForBJet(UInt_t nIdx, Float_t &fJetPt, Float_t &fJetEta, Float_t &fJetPhi, Float_t &fJetMass) 
    {return true;};
  
  // In uncertainty study we need to switch the kinematic variables of jets
  // The following variables are for this switch
  // In the topObjectSelection.cc these variables are used instead of Jet_pt, Jet_mass, and so on.
  // In default, these are same as the original ones, but when a user wants to study systematic uncertainty 
  // so that he/she needs to switch them to the evaluated ones, 
  // just touching them in anlalyser class will be okay, and this is for it.
  virtual void GetJetMassPt(UInt_t nIdx, 
    Float_t &fJetMass, Float_t &fJetPt, Float_t &fJetEta, Float_t &fJetPhi);
  
  Int_t GetMatchGenJet(UInt_t nIdxJet, Float_t fResolution);
  
public:
  std::vector<TParticle> muonSelection();
  std::vector<TParticle> elecSelection();
  std::vector<TParticle> vetoMuonSelection();
  std::vector<TParticle> vetoElecSelection();
  std::vector<TLorentzVector> recoleps;
  std::vector<TParticle> jetSelection();
  std::vector<TParticle> bjetSelection();

  std::vector<TParticle> genJetSelection();

  topObjectSelection(TTree *tree=0, TTree *had=0, TTree *hadTruth=0, Bool_t isMC = false, UInt_t unFlag = 0);
  topObjectSelection(TTree *tree=0, Bool_t isMC=false, UInt_t unFlag = 0) : 
    topObjectSelection(tree, 0, 0, isMC, unFlag) {}
  ~topObjectSelection() {
    if ( jecUnc != NULL ) delete jecUnc;
    if ( rndEngine != NULL ) delete rndEngine;
  }
  
  // In this function you need to set all the cut conditions in the above
  // If you do not set this function up (so that you didn't set the cuts), the compiler will deny your code, 
  // so you can be noticed that you forgot the setting up.
  // And you don't need to run this function indivisually; it will be run in the creator of this class.
  virtual int SetCutValues() = 0;
  
  enum {
    OptBit_JER_Up = 0, 
    OptBit_JER_Dn, 
    OptBit_JES_Up, 
    OptBit_JES_Dn
  };
  
  enum {
    OptFlag_JER_Up = ( 1 << OptBit_JER_Up ), 
    OptFlag_JER_Dn = ( 1 << OptBit_JER_Dn ), 
    OptFlag_JES_Up = ( 1 << OptBit_JES_Up ), 
    OptFlag_JES_Dn = ( 1 << OptBit_JES_Dn )
  };
};

#endif
