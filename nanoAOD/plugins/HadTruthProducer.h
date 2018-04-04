#ifndef HadTruthProducer_H
#define HadTruthProducer_H
/*
  From cmssw RecoVertex/V0Producer/src/V0Fitter.cc
  matching based on Validation/RecoVertex/src/V0Validator.cc
*/

#include<memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "TLorentzVector.h"

using namespace edm;
using namespace std;

class HadTruthProducer : public edm::stream::EDProducer<>
{
public:
  explicit HadTruthProducer(const edm::ParameterSet & iConfig);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private: 
  void produce( edm::Event&, const edm::EventSetup& ) override;

  reco::LeafCandidate getCandidate(const TrackingParticle* tp) {
    return reco::LeafCandidate( tp->charge(), tp->p4(), tp->vertex(), tp->pdgId(), tp->status() );
  };
  bool isGenHadFrom(const reco::GenParticle* particle, int pdgId, int count,int & GenHadFromQuark, bool & GenHadFromTop);
  bool isHadFrom(const reco::GenParticleRef &particle, int pdgId, int count, int & hadFromQuark, bool & hadFromTop);
  void motherTracking(int PID, const TrackingVertex trackVertex, const TrackingParticle *decayTrk, int count, int & GenHadFromQuark, bool & GenHadFromTop);

  int trackingVertex_pdgId(const TrackingVertex* tv);
  const reco::GenParticleRef getMother(const TrackingParticleRef& tp);

  edm::EDGetTokenT<reco::RecoToSimCollection> recoRecoToSim_;
  edm::EDGetTokenT<reco::SimToRecoCollection> recoSimToReco_;
  edm::EDGetTokenT<reco::VertexCompositeCandidateCollection > hadronCands_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > genLabel_;
  edm::EDGetTokenT<TrackingVertexCollection> trackingVertexLabel_;
  edm::EDGetTokenT<TrackingParticleCollection> trackingParticleLabel_;
  
  const int pion_pdgId_ = 211, kaon_pdgId_ = 321, proton_pdgId_ = 2212;
  const float pion_m_ = 0.1396, kaon_m_ = 0.4937, proton_m_ = 0.938272;

  const int jpsi_pdgId_ = 443, d0_pdgId_ = 421, dstar_pdgId_ = 413;
  const float jpsi_m_ = 3.096, d0_m_ = 1.865, dstar_m_ = 2.010;

  const int kshort_pdgId_ = 310, lambda_pdgId_ = 3122;
  const float kshort_m_ = 0.4976, lambda_m_ = 1.11568;
  
  const int lambdab_pdgId_ = 5122;
  const float lambdab_m_ = 5.61958;
};

DEFINE_FWK_MODULE(HadTruthProducer);
#endif
