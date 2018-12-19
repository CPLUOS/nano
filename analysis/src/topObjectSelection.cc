#include "nano/analysis/interface/topObjectSelection.h"

using std::vector;

topObjectSelection::topObjectSelection(TTree *tree, TTree *had, TTree *hadTruth, Bool_t isMC, UInt_t unFlag) : 
  nanoBase(tree, had, hadTruth, isMC), m_unFlag(unFlag), jecUnc(NULL), rndEngine(NULL)
{
  //if ( ( m_unFlag & ( OptFlag_JES_Up | OptFlag_JES_Dn | OptFlag_JER_Up | OptFlag_JER_Dn ) ) != 0 )
  if ( m_isMC ) {
    std::string env = getenv("CMSSW_BASE");
    
    std::string strPathJetResSFObj = env + "/src/nano/analysis/data/jer/"
      "Summer16_25nsV1_MC_SF_AK4PFchs.txt";
    jetResSFObj = JMENano::JetResolutionScaleFactor(strPathJetResSFObj.c_str());
    
    std::string strPathJetResObj = env + "/src/nano/analysis/data/jer/"
      "Summer16_25nsV1_MC_PtResolution_AK4PFchs.txt";
    jetResObj = JMENano::JetResolution(strPathJetResObj.c_str());
    
    rndEngine = new TRandom3(12345);
    
    if ( ( m_unFlag & ( OptFlag_JES_Up | OptFlag_JES_Dn ) ) != 0 ) {
      std::string strPathJecUnc = env + "/src/nano/analysis/data/jer/"
        //"Summer16_23Sep2016V4_MC_Uncertainty_AK4PFchs.txt";
        "Summer16_23Sep2016V4_MC_UncertaintySources_AK4PFchs.txt";
      
      JetCorrectorParameters JetCorPar(strPathJecUnc, "Total");
      jecUnc = new JetCorrectionUncertainty(JetCorPar);
    }
  }
  
  m_fDRcone_JER = 0.4; // For AK4 jets
  m_fResFactorMathcer = 3; // According to https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution
}

vector<TParticle> topObjectSelection::elecSelection() {
  vector<TParticle> elecs; 
  for (UInt_t i = 0; i < nElectron; ++i){
    if (Electron_pt[i] < cut_ElectronPt) continue;
    if (std::abs(Electron_eta[i]) > cut_ElectronEta) continue;
    if (cut_ElectronIDType != NULL && cut_ElectronIDType[i] < cut_ElectronIDCut) continue; 
    float el_scEta = Electron_deltaEtaSC[i] + Electron_eta[i];
    if ( cut_ElectronSCEtaLower < std::abs(el_scEta) && std::abs(el_scEta) < cut_ElectronSCEtaUpper ) continue;
    if ( Electron_pfRelIso03_all[ i ] > cut_ElectronRelIso03All ) continue;
    if ( !additionalConditionForElectron(i) ) continue;
    
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Electron_pt[i], Electron_eta[i], Electron_phi[i], Electron_mass[i]);

    auto elec = TParticle();
    elec.SetPdgCode(11*Electron_charge[i]*-1);
    elec.SetMomentum(mom);
    elec.SetWeight(el_scEta);
    elec.SetFirstMother(i);
    elecs.push_back(elec);
  }
  return elecs;
}

vector<TParticle> topObjectSelection::muonSelection() {
  vector<TParticle> muons; 
  for (UInt_t i = 0; i < nMuon; ++i){
    if (cut_MuonIDType != NULL && !cut_MuonIDType[i]) continue;
    if (Muon_pt[i] < cut_MuonPt) continue;
    if (std::abs(Muon_eta[i]) > cut_MuonEta) continue;
    if (Muon_pfRelIso04_all[i] > cut_MuonRelIso04All) continue;
    if ( !additionalConditionForMuon(i) ) continue;
    
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]);
    auto muon = TParticle();
    muon.SetPdgCode(13*Muon_charge[i]*-1);
    muon.SetMomentum(mom);
    muon.SetFirstMother(i);
    muons.push_back(muon);
  }
  return muons;
}

/// TODO: Implement veto electron selection (currently empty on TTbarXSecSynchronization)
vector<TParticle> topObjectSelection::vetoElecSelection() {
  vector<TParticle> elecs; 
  for (UInt_t i = 0; i < nElectron; ++i){
    if (Electron_pt[i] < cut_VetoElectronPt) continue;
    if (std::abs(Electron_eta[i]) > cut_VetoElectronEta) continue;
    if (cut_VetoElectronIDType != NULL && cut_VetoElectronIDType[i] < cut_VetoElectronIDCut) continue; 
    float el_scEta = Electron_deltaEtaSC[i] + Electron_eta[i];
    if ( cut_VetoElectronSCEtaLower < std::abs(el_scEta) && 
                                      std::abs(el_scEta) < cut_VetoElectronSCEtaUpper ) continue;
    if ( Electron_pfRelIso03_all[ i ] > cut_VetoElectronRelIso03All ) continue;
    if ( !additionalConditionForVetoElectron(i) ) continue;
    
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Electron_pt[i], Electron_eta[i], Electron_phi[i], Electron_mass[i]);
    auto elec = TParticle();
    elec.SetPdgCode(11*Electron_charge[i]*-1);
    elec.SetMomentum(mom);
    elec.SetWeight(el_scEta);
    elec.SetFirstMother(i);
    elecs.push_back(elec);
  }
  return elecs;
}

vector<TParticle> topObjectSelection::vetoMuonSelection() {
  vector<TParticle> muons; 
  for (UInt_t i = 0; i < nMuon; ++i){
    // if (!Muon_looseId[i]) continue;

    // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Loose_Muon
    if (cut_VetoMuonIDType != NULL && !cut_VetoMuonIDType[i]) continue;
    if (Muon_pt[i] < cut_VetoMuonPt) continue;
    if (std::abs(Muon_eta[i]) > cut_VetoMuonEta) continue;
    if (Muon_pfRelIso04_all[i] > cut_VetoMuonRelIso04All) continue;
    if ( !additionalConditionForVetoMuon(i) ) continue;
    
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]);
    auto muon = TParticle();
    muon.SetPdgCode(13*Muon_charge[i]*-1);
    muon.SetMomentum(mom);
    muon.SetFirstMother(i);
    muons.push_back(muon);
  }
  return muons;
}

vector<TParticle> topObjectSelection::genJetSelection() {
  vector<TParticle> jets;
  for (UInt_t i = 0; i < nJet; ++i){
    if (GenJet_pt[i] < cut_GenJetPt) continue;
    if (std::abs(GenJet_eta[i]) > cut_GenJetEta) continue;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(GenJet_pt[i], GenJet_eta[i], GenJet_phi[i], GenJet_mass[i]);
    bool hasOverLap = false;
    for (auto lep : recoleps){
        if (mom.TLorentzVector::DeltaR(lep) < cut_GenJetConeSizeOverlap) hasOverLap = true;
    }
    if (hasOverLap) continue;
    if ( !additionalConditionForGenJet(i) ) continue;
    
    auto jet = TParticle();
    jet.SetMomentum(mom);
    jet.SetFirstMother(i);
    jets.push_back(jet);
  }
  return jets;
}

vector<TParticle> topObjectSelection::jetSelection() {
  double csvWgtHF = 1.0, csvWgtLF = 1.0, csvWgtC = 1.0, csvWgtTotal = 1.0;
  double csvWgtHF_up = 1.0, csvWgtLF_up = 1.0, csvWgtC_up = 1.0, csvWgtTotal_up = 1.0;
  double csvWgtHF_dn = 1.0, csvWgtLF_dn = 1.0, csvWgtC_dn = 1.0, csvWgtTotal_dn = 1.0;
  vector<TParticle> jets;
  for (UInt_t i = 0; i < nJet; ++i){
    Float_t fJetMass, fJetPt, fJetEta, fJetPhi;
    GetJetMassPt(i, fJetMass, fJetPt, fJetEta, fJetPhi);
    
    if (fJetPt < cut_JetPt) continue;
    if (std::abs(fJetEta) > cut_JetEta) continue;
    if (Jet_jetId[i] < cut_JetID) continue;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(fJetPt, fJetEta, fJetPhi, fJetMass);
    bool hasOverLap = false;
    for (auto lep : recoleps){
        if (mom.TLorentzVector::DeltaR(lep) < cut_JetConeSizeOverlap) hasOverLap = true;
    }
    if (hasOverLap) continue;
    if ( !additionalConditionForJet(i, fJetPt, fJetEta, fJetPhi, fJetMass) ) continue;
    
    auto jet = TParticle();
    jet.SetMomentum(mom);
    jet.SetFirstMother(i);
    jets.push_back(jet);
    if (m_isMC) {
      BTagEntry::JetFlavor JF = BTagEntry::FLAV_UDSG;
      if (abs(Jet_hadronFlavour[i]) == 5) {
        JF = BTagEntry::FLAV_B;
        double iCSVWgtHF = m_btagSF.eval_auto_bounds("central", JF, abs(Jet_eta[i]), Jet_pt[i], Jet_btagCSVV2[i]);
        double iCSVWgtHF_up = m_btagSF_up.eval_auto_bounds("up", JF, abs(Jet_eta[i]), Jet_pt[i], Jet_btagCSVV2[i]);
        double iCSVWgtHF_dn = m_btagSF_dn.eval_auto_bounds("down", JF, abs(Jet_eta[i]), Jet_pt[i], Jet_btagCSVV2[i]);
        if (iCSVWgtHF != 0) csvWgtHF *= iCSVWgtHF;
        if (iCSVWgtHF_up != 0) csvWgtHF_up *= iCSVWgtHF_up;
        if (iCSVWgtHF_dn != 0) csvWgtHF_dn *= iCSVWgtHF_dn;
      }
      else if (abs(Jet_hadronFlavour[i]) == 4) {
        JF = BTagEntry::FLAV_C;
        double iCSVWgtC = m_btagSF.eval_auto_bounds("central", JF, abs(Jet_eta[i]), Jet_pt[i], Jet_btagCSVV2[i]);
        double iCSVWgtC_up = m_btagSF_up.eval_auto_bounds("up", JF, abs(Jet_eta[i]), Jet_pt[i], Jet_btagCSVV2[i]);
        double iCSVWgtC_dn = m_btagSF_dn.eval_auto_bounds("down", JF, abs(Jet_eta[i]), Jet_pt[i], Jet_btagCSVV2[i]);
        if (iCSVWgtC != 0) csvWgtC *= iCSVWgtC;
        if (iCSVWgtC_up != 0) csvWgtC_up *= iCSVWgtC_up;
        if (iCSVWgtC_dn != 0) csvWgtC_dn *= iCSVWgtC_dn;
      }
      else {
        JF = BTagEntry::FLAV_UDSG;
        double iCSVWgtLF = m_btagSF.eval_auto_bounds("central", JF, abs(Jet_eta[i]), Jet_pt[i], Jet_btagCSVV2[i]);
        double iCSVWgtLF_up = m_btagSF_up.eval_auto_bounds("up", JF, abs(Jet_eta[i]), Jet_pt[i], Jet_btagCSVV2[i]);
        double iCSVWgtLF_dn = m_btagSF_dn.eval_auto_bounds("down", JF, abs(Jet_eta[i]), Jet_pt[i], Jet_btagCSVV2[i]);
        if (iCSVWgtLF != 0) csvWgtLF *= iCSVWgtLF;
        if (iCSVWgtLF_up != 0) csvWgtLF_up *= iCSVWgtLF_up;
        if (iCSVWgtLF_dn != 0) csvWgtLF_dn *= iCSVWgtLF_dn;
      } 
    }
  }
  csvWgtTotal = csvWgtHF * csvWgtC * csvWgtLF;  
  csvWgtTotal_up = csvWgtHF_up * csvWgtC_up * csvWgtLF_up;  
  csvWgtTotal_dn = csvWgtHF_dn * csvWgtC_dn * csvWgtLF_dn;  

  b_btagweight = csvWgtTotal;
  b_btagweight_up = csvWgtTotal_up;
  b_btagweight_dn = csvWgtTotal_dn;
  
  return jets;
}

vector<TParticle> topObjectSelection::bjetSelection() {
  vector<TParticle> bjets;
  for (UInt_t i = 0; i < nJet; ++i ) {
    Float_t fJetMass, fJetPt, fJetEta, fJetPhi;
    GetJetMassPt(i, fJetMass, fJetPt, fJetEta, fJetPhi);
    
    if (fJetPt < cut_BJetPt) continue;
    if (std::abs(fJetEta) > cut_BJetEta) continue;
    if (Jet_jetId[i] < cut_BJetID) continue;
    if (cut_BJetTypeBTag[i] < cut_BJetBTagCut) continue;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(fJetPt, fJetEta, fJetPhi, fJetMass);
    bool hasOverLap = false;
    for (auto lep : recoleps) {
      if (mom.TLorentzVector::DeltaR(lep) < cut_BJetConeSizeOverlap) hasOverLap = true;
    }
    if (hasOverLap) continue;
    if ( !additionalConditionForBJet(i, fJetPt, fJetEta, fJetPhi, fJetMass) ) continue;
    
    auto bjet = TParticle();
    bjet.SetMomentum(mom);
    bjet.SetFirstMother(i);
    bjets.push_back(bjet);
  }
  return bjets;
}


// In uncertainty study we need to switch the kinematic variables of jets
// The following variables are for this switch
// In the topObjectSelection.cc these variables are used instead of Jet_pt, Jet_mass, and so on.
// In default, these are same as the original ones, but when a user wants to study systematic uncertainty 
// so that he/she needs to switch them to the evaluated ones, 
// just touching them in anlalyser class will be okay, and this is for it.
void topObjectSelection::GetJetMassPt(UInt_t nIdx, 
  Float_t &fJetMass, Float_t &fJetPt, Float_t &fJetEta, Float_t &fJetPhi) 
{
  Float_t fCorrFactor = 1.0;
  
  fJetMass = Jet_mass[ nIdx ];
  fJetPt   = Jet_pt[ nIdx ];
  fJetEta  = Jet_eta[ nIdx ];
  fJetPhi  = Jet_phi[ nIdx ];
  
  //if ( m_isMC && ( m_unFlag & ( OptFlag_JER_Up | OptFlag_JER_Dn | OptFlag_JES_Up | OptFlag_JES_Dn ) ) != 0 )
  if ( m_isMC ) {
    // Evaluating the central part cJER of the factor
    JMENano::JetParameters jetPars = {{JMENano::Binning::JetPt, fJetPt},
                                  {JMENano::Binning::JetEta, fJetEta},
                                  {JMENano::Binning::Rho, fixedGridRhoFastjetAll}};
    
    const double jetRes = jetResObj.getResolution(jetPars); // Note: this is relative resolution.
    const double cJER = jetResSFObj.getScaleFactor(jetPars, 
      ( ( m_unFlag & ( OptFlag_JER_Up | OptFlag_JER_Dn ) ) == 0 ? Variation::NOMINAL : 
        ( ( m_unFlag & OptFlag_JER_Up ) != 0 ? Variation::UP : Variation::DOWN ) ));
    
    // We need corresponding genJet
    //Int_t nIdxGen = Jet_genJetIdx[ nIdx ];
    Int_t nIdxGen = GetMatchGenJet(nIdx, fJetPt * jetRes);
    const double genJetPt = GenJet_pt[ nIdxGen ];
    
    // JER (nominal and up and down) - apply scaling method if matched genJet is found,
    //       apply gaussian smearing method if unmatched
    /*if ( nIdxGen >= 0 && //deltaR(genJet->p4(), jet.p4()) < 0.2 && // From CATTool
         std::abs(genJetPt - fJetPt) < jetRes * 3 * fJetPt ) 
      fCorrFactor = std::max(0., (genJetPt + ( fJetPt - genJetPt ) * cJER) / fJetPt);*/
    if ( nIdxGen >= 0 ) {
      fCorrFactor = ( genJetPt + ( fJetPt - genJetPt ) * cJER ) / fJetPt;
    } else {
      const double smear = rndEngine->Gaus(0, 1);
      fCorrFactor = ( cJER <= 1 ? 1 : 1 + smear * jetRes * sqrt(cJER * cJER - 1) );
    }
    
    if ( ( m_unFlag & ( OptFlag_JES_Up | OptFlag_JES_Dn ) ) != 0 ) { // JES
      // The evaluator needs pT and eta of the current jet
      jecUnc->setJetPt(fCorrFactor * fJetPt);
      jecUnc->setJetEta(fJetEta);
      
      if ( ( m_unFlag & OptFlag_JES_Up ) != 0 ) {
        fCorrFactor *= 1 + jecUnc->getUncertainty(true);
      } else {
        fCorrFactor *= 1 - jecUnc->getUncertainty(true);
      }
    }
  }
  
  fJetMass *= fCorrFactor;
  fJetPt   *= fCorrFactor;
}


Int_t topObjectSelection::GetMatchGenJet(UInt_t nIdxJet, Float_t fResolution) {
  UInt_t i;
  
  double dEta, dPhi, dR;
  double dRFound = m_fDRcone_JER;
  UInt_t nIdxFound = 999;
  
  for ( i = 0 ; i < nGenJet ; i++ ) {
    dEta = Jet_eta[ nIdxJet ] - GenJet_eta[ i ];
    dPhi = std::abs(Jet_phi[ nIdxJet ] - GenJet_phi[ i ]);
    if ( dPhi > (double)M_PI ) dPhi -= (double)( 2 * M_PI );
    
    dR = std::sqrt(dEta * dEta + dPhi * dPhi);
    
    if ( dR >= (double)m_fDRcone_JER * 0.5 ) continue;
    if ( dRFound > dR ) {
      if ( std::abs(Jet_pt[ nIdxJet ] - GenJet_pt[ i ]) >= m_fResFactorMathcer * fResolution ) continue;
      
      dRFound = dR;
      nIdxFound = i;
    }
  }
 
  return ( dRFound < 0.75 * (double)m_fDRcone_JER ? (Int_t)nIdxFound : -1 );
}


