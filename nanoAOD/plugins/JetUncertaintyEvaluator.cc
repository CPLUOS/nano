#include "JetMetUncertainty.h"


//#define debugMode
using namespace edm;
using namespace std;


JetMetUncertainty::JetMetUncertainty(const edm::ParameterSet &iConfig) : 
  src_(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("src"))),
  rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
  payloadName_(iConfig.getParameter<std::string>("payloadName")),
  jetResFilePath_(edm::FileInPath(iConfig.getParameter<std::string>("jetResFile")).fullPath()),
  jetResSFFilePath_(edm::FileInPath(iConfig.getParameter<std::string>("jetResSFFile")).fullPath())
{
  produces<nanoaod::FlatTable>();
}


void JetMetUncertainty::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


void JetMetUncertainty::beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup&) {
  edm::Service<edm::RandomNumberGenerator> rng;
  rng_ = &rng->getEngine(lumi.index());
}


void JetMetUncertainty::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  runOnMC_ = !iEvent.isRealData();
  
  Handle<pat::JetCollection> jets;
  iEvent.getByToken(src_, jets);
  
  // The following codes are from: 
  //  https://github.com/vallot/CATTools/blob/cat80x/CatProducer/plugins/CATJetProducer.cc
  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoToken_, rhoHandle);
  const double rho = *rhoHandle;
  
  if ( !payloadName_.empty() ) {
    // temp measure - payloadName should be AK4PFchs, but PHYS14_25_V2 does not have uncertainty
    edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
    iSetup.get<JetCorrectionsRecord>().get(payloadName_, JetCorParColl);
    JetCorrectorParameters const & JetCorPar = ( *JetCorParColl )[ "Uncertainty" ];
    
    jecUnc = new JetCorrectionUncertainty(JetCorPar);
  }
  
  JME::JetResolution jetResObj;
  JME::JetResolutionScaleFactor jetResSFObj;
  
  if ( runOnMC_ ) {
    jetResObj = JME::JetResolution(jetResFilePath_);
    jetResSFObj = JME::JetResolutionScaleFactor(jetResSFFilePath_);
  }
  
  std::vector<float> vecJetJESUp;
  std::vector<float> vecJetJESDn;
  std::vector<float> vecJetJERUp;
  std::vector<float> vecJetJERDn;
  
  for ( unsigned int i = 0 ; i < jets->size() ; i++ ) {
    const auto & jet = ( *jets )[ i ];
    
    double jetPt = jet.pt();
    double jetEta = jet.eta();
    
    // Computing JEC uncertainty
    float fJESUp = 1, fJESDn = 1;
    
    if ( !payloadName_.empty() ) {
      jecUnc->setJetEta(jetEta);
      jecUnc->setJetPt(jetPt); // here you must use the CORRECTED jet pt
      fJESUp = 1 + jecUnc->getUncertainty(true);
      
      jecUnc->setJetEta(jetEta);
      jecUnc->setJetPt(jetPt); // here you must use the CORRECTED jet pt
      fJESDn = 1 - jecUnc->getUncertainty(false);
    }
    
    // Computing JER uncertainty
    float fJERUp = 1, fJERDn = 1;
    
    if ( runOnMC_ ) {
      // adding genJet
      auto genJet = jet.genJetFwdRef();
      
      JME::JetParameters jetPars = {{JME::Binning::JetPt, jetPt},
                                    {JME::Binning::JetEta, jetEta},
                                    {JME::Binning::Rho, rho}};
      const double jetRes = jetResObj.getResolution(jetPars); // Note: this is relative resolution.
      const double cJERUp = jetResSFObj.getScaleFactor(jetPars, Variation::UP);
      const double cJERDn = jetResSFObj.getScaleFactor(jetPars, Variation::DOWN);
      
      // JER - apply scaling method if matched genJet is found,
      //       apply gaussian smearing method if unmatched
      if ( genJet.isNonnull() && deltaR(genJet->p4(), jet.p4()) < 0.2
           && std::abs(genJet->pt() - jetPt) < jetRes * 3 * jetPt ) 
      {
        const double genJetPt = genJet->pt();
        const double dPt = jetPt - genJetPt;
        
        fJERUp = std::max(0., (genJetPt + dPt * cJERUp) / jetPt);
        fJERDn = std::max(0., (genJetPt + dPt * cJERDn) / jetPt);
      } else {
        const double smear = CLHEP::RandGaussQ::shoot(rng_);
        
        fJERUp = ( cJERUp <= 1 ? 1 : 1 + smear * jetRes * sqrt(cJERUp * cJERUp - 1) );
        fJERDn = ( cJERDn <= 1 ? 1 : 1 + smear * jetRes * sqrt(cJERDn * cJERDn - 1) );
      }
    }
    
    vecJetJESUp.push_back(fJESUp);
    vecJetJESDn.push_back(fJESDn);
    vecJetJERUp.push_back(fJERUp);
    vecJetJERDn.push_back(fJERDn);
  }
  
  auto outJetCor = std::make_unique<nanoaod::FlatTable>(jets->size(), "Jet", false, true);
  //outJetCor->setDoc("Jet energy scale factor uncertainty");
  
  outJetCor->addColumn<float>("jes_up", vecJetJESUp, "Scale factor for JES up", 
    nanoaod::FlatTable::FloatColumn);
  outJetCor->addColumn<float>("jes_dn", vecJetJESDn, "Scale factor for JES down", 
    nanoaod::FlatTable::FloatColumn);
  
  outJetCor->addColumn<float>("jer_up", vecJetJERUp, "Scale factor for JER up", 
    nanoaod::FlatTable::FloatColumn);
  outJetCor->addColumn<float>("jer_dn", vecJetJERDn, "Scale factor for JER down", 
    nanoaod::FlatTable::FloatColumn);
  
  iEvent.put(move(outJetCor));
}

DEFINE_FWK_MODULE(JetMetUncertainty);


