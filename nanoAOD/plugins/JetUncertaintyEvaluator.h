#ifndef JetMetUncertainty_H
#define JetMetUncertainty_H

#include<memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

#include "JetMETCorrections/Modules/interface/JetResolution.h"

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGaussQ.h"

//#define debugMode

using namespace edm;
using namespace std;


class JetMetUncertainty: public edm::stream::EDProducer<> {
public:  
  explicit JetMetUncertainty(const edm::ParameterSet &iConfig);
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  void beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup&) override;
  void produce( edm::Event&, const edm::EventSetup& ) override;
  
  edm::EDGetTokenT<pat::JetCollection> src_;
  
  edm::EDGetTokenT<double> rhoToken_;
  std::string payloadName_;
  const std::string jetResFilePath_, jetResSFFilePath_;
  
  bool runOnMC_;
  
  JetCorrectionUncertainty *jecUnc;
  CLHEP::HepRandomEngine* rng_;
};

#endif // JetMetUncertainty_H


