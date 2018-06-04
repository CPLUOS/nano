#ifndef HadronProducer_H
#define HadronProducer_H
/*
  From cmssw RecoVertex/V0Producer/src/V0Fitter.cc
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

#include "DataFormats/Math/interface/deltaR.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "TLorentzVector.h"

//#define debugMode

using namespace edm;
using namespace std;

class HadronProducer : public edm::stream::EDProducer<> {
  typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
  typedef ROOT::Math::SVector<double, 3> SVector3;

public:
  struct hadronCandidate {
    reco::VertexCompositeCandidate vcc;
    pat::Jet jet;
    int nJet, nDau;
    float diffMass, lxy, lxyErr, l3D, l3DErr, dca, angleXY, angleXYZ;
    int idx;
    int dau1_idx = -1;
    int dau2_idx = -1;
  };

  class hadronCandidateCollection : public std::vector<hadronCandidate> {
    public:
      int current_index_ = 0;
      void push_back(const hadronCandidate& h) {
        hadronCandidate cand = h;
        cand.idx = this->current_index_;
        this->current_index_ += 1;
        this->std::vector<hadronCandidate>::push_back(cand);
      }
      void clear() {
        this->std::vector<hadronCandidate>::clear();
        this->current_index_ = 0;
      }
      iterator insert(hadronCandidateCollection::iterator position, hadronCandidateCollection::iterator first, hadronCandidateCollection::iterator last) {
        hadronCandidateCollection::iterator step = first;
        while(step < last) {
          step->idx = this->current_index_;
          this->current_index_ += 1;
          step++;
        }
        this->std::vector<hadronCandidate>::insert(position, first, last);
        return last;
      }
  };
  
  explicit HadronProducer(const edm::ParameterSet & iConfig);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  static const int W_pdgId_ = 24, S_pdgId_ = 3, B_pdgId_ = 5, C_pdgId_ = 4;
  
  static const int pion_pdgId_ = 211, kaon_pdgId_ = 321, proton_pdgId_ = 2212;
  static constexpr float pion_m_ = 0.1396, kaon_m_ = 0.4937, proton_m_ = 0.938272;

  static const int jpsi_pdgId_ = 443, d0_pdgId_ = 421, dstar_pdgId_ = 413;
  static constexpr float jpsi_m_ = 3.096, d0_m_ = 1.865, dstar_m_ = 2.010;

  static const int kshort_pdgId_ = 310, lambda_pdgId_ = 3122;
  static constexpr float kshort_m_ = 0.4976, lambda_m_ = 1.11568;
  
  static const int lambdab_pdgId_ = 5122;
  static constexpr float lambdab_m_ = 5.61958;
  
private:
  void produce( edm::Event&, const edm::EventSetup& ) override;

  reco::VertexCompositeCandidate fit(vector<reco::Candidate*>& cands,
				     reco::Vertex& pv, int pdgId,
				     float &dca, float &angleXY, float &angleXYZ);
    
  SVector3 getDistanceVector(int dim, reco::VertexCompositeCandidate& vertex,reco::Vertex& pv);
  pair<float, float> getDistance(int dim, reco::VertexCompositeCandidate& vertex,reco::Vertex& pv);
  
  vector<hadronCandidate> findJPsiCands(vector<reco::Candidate*> &leptons, reco::Vertex& pv, int nJet, const pat::Jet & aPatJet);
  vector<hadronCandidate> findD0Cands(vector<reco::Candidate*> &chargedHads, reco::Vertex& pv, int nJet, const pat::Jet & aPatJet);
  vector<hadronCandidate> findDStarCands(vector<HadronProducer::hadronCandidate>& d0cands, vector<reco::Candidate*> &chargedHads, reco::Vertex& pv, int nJet, const pat::Jet & aPatJet);

  vector<hadronCandidate> findKShortCands(vector<reco::Candidate*> &chargedHads, reco::Vertex& pv, int nJet);
  vector<hadronCandidate> findLambdaCands(vector<reco::Candidate*> &chargedHads, reco::Vertex& pv, int nJet);
  vector<hadronCandidate> findLambdaBCands(vector<hadronCandidate>& LambdaCands,vector<hadronCandidate>& jpsiCands, reco::Vertex& pv, int nJet, const pat::Jet & aPatJet);

  edm::EDGetTokenT<edm::View<pat::Jet> > jetLabel_;
  edm::EDGetTokenT<reco::VertexCollection> vertexLabel_;
  edm::EDGetTokenT<reco::CandidateView> pfCandidates_;
  edm::ESHandle<TransientTrackBuilder> trackBuilder_;

  // cuts on initial track selection
  float tkChi2Cut_;
  int tkNHitsCut_;
  float tkPtCut_;
  float tkIPSigXYCut_;
  float tkIPSigZCut_;
  // cuts on the vertex
  float vtxChi2Cut_;
  float vtxDecaySigXYCut_;
  float vtxDecaySigXYZCut_;
  // miscellaneous cuts
  float tkDCACut_;
  float cosThetaXYCut_;
  float cosThetaXYZCut_;  
};

#endif
