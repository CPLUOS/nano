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
  bool isHadFrom(const reco::GenParticleRef &particle, int pdgId, int count, int & hadFromQuark, bool & hadFromTop, bool & hadFromW, bool & hadFromS, bool & hadFromB, bool & hadFromC);
  void motherTracking(int PID, const TrackingVertex trackVertex, const TrackingParticle *decayTrk, int count, int & GenHadFromQuark, bool & GenHadFromTop);
  const reco::GenParticle* genChecking(int PID, const TrackingParticle *decayTrk);

  int trackingVertex_pdgId(const TrackingVertex* tv);
  const reco::GenParticleRef getMother(const TrackingParticleRef& tp);

  edm::EDGetTokenT<reco::RecoToSimCollection> recoRecoToSim_;
  edm::EDGetTokenT<reco::SimToRecoCollection> recoSimToReco_;
  edm::EDGetTokenT<reco::VertexCompositeCandidateCollection > hadronCands_;
  edm::EDGetTokenT<std::vector<std::vector<int>>> hadronIndices_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > genLabel_;
  edm::EDGetTokenT<TrackingVertexCollection> trackingVertexLabel_;
  edm::EDGetTokenT<TrackingParticleCollection> trackingParticleLabel_;
};

#endif
