#include "HadTruthProducer.h"
//#define debugMode

HadTruthProducer::HadTruthProducer(const edm::ParameterSet & iConfig) :
  recoRecoToSim_(consumes<reco::RecoToSimCollection>(iConfig.getParameter<edm::InputTag>("recoRecoToSim"))),
  recoSimToReco_(consumes<reco::SimToRecoCollection>(iConfig.getParameter<edm::InputTag>("recoSimToReco"))),
  hadronCands_(consumes<reco::VertexCompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("hadronCands"))),
  genLabel_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genLabel"))),
  trackingVertexLabel_(consumes<TrackingVertexCollection>(iConfig.getParameter<edm::InputTag>("trackingVertexLabel"))),
  trackingParticleLabel_(consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("trackingParticleLabel")))
{
  produces<nanoaod::FlatTable>("hadTruth");
  produces<nanoaod::FlatTable>("genHadron");
  produces<std::vector<reco::LeafCandidate> >();
}

void
HadTruthProducer::produce( edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  Handle<reco::VertexCompositeCandidateCollection> hadronCands;
  iEvent.getByToken(hadronCands_,hadronCands);

  Handle<reco::RecoToSimCollection> recotosimH;
  iEvent.getByToken(recoRecoToSim_, recotosimH);
  auto recotosim = *recotosimH.product();
  Handle<reco::SimToRecoCollection> simtorecoH;
  iEvent.getByToken(recoSimToReco_, simtorecoH);
  auto simtoreco = *simtorecoH.product();
  
  edm::Handle<TrackingVertexCollection> trackingVertexs;
  iEvent.getByToken(trackingVertexLabel_, trackingVertexs);
  edm::Handle<TrackingParticleCollection> trackingParticles;
  iEvent.getByToken(trackingParticleLabel_, trackingParticles);

  vector<int> nmatchedv;

  vector<int> isHadFromTsb;
  vector<uint8_t> isHadFromTop;
  
  for (auto& cand : *hadronCands) {
    int count = 0;
    int hadFromQuark = -99;
    bool hadFromTop = false;

    // for dstar and lambdaB, need to match with grand mother
    reco::GenParticleRef trueHad;
    
    int numberOfDaughters = cand.numberOfDaughters();
    int nmatched = 0;
    for (int ndau =0; ndau < numberOfDaughters; ++ndau) {
      auto rcCand = dynamic_cast<const reco::RecoChargedCandidate*>(cand.daughter(ndau));
      RefToBase<reco::Track> track(rcCand->track());
      if (recotosim.find(track) != recotosim.end()) {
	
	TrackingParticleRef tpref = recotosim[track].begin()->first;
	if (rcCand->pdgId() == tpref->pdgId()) {
	  auto mother = getMother(tpref);
	  if (mother.isNull()) { continue; }
	  
	  if (trueHad.isNull()) { trueHad = mother; }
	  if (mother != trueHad) { continue; }
	  if (abs(mother->pdgId()) == cand.pdgId()) { nmatched++; }
	}
      }
    }
    nmatchedv.push_back(nmatched);

    if (nmatched == 2) {
      isHadFrom(trueHad, 6, count, hadFromQuark, hadFromTop);
    }

    isHadFromTsb.push_back(hadFromQuark);
    isHadFromTop.push_back(hadFromTop);
  }
  
  auto hadTruthTable = make_unique<nanoaod::FlatTable>(hadronCands->size(),"hadTruth",false);
  hadTruthTable->addColumn<int>("nMatched",nmatchedv,"no. of dau match",nanoaod::FlatTable::IntColumn);
  hadTruthTable->addColumn<int>("isHadFromTsb",isHadFromTsb,"Hadron from t->s/b",nanoaod::FlatTable::IntColumn);
  hadTruthTable->addColumn<uint8_t>("isHadFromTop",isHadFromTop,"Hadron from Top",nanoaod::FlatTable::UInt8Column);  
  iEvent.put(move(hadTruthTable),"hadTruth");
  
  auto candidates = make_unique<std::vector<reco::LeafCandidate>>();
  vector<int> isGenHadFromTsb;
  vector<uint8_t> isGenHadFromTop;
  vector<uint8_t> inVol;
  vector<int> dau1_pdgId, dau2_pdgId;
  vector<float> dau1_pt, dau1_eta, dau1_phi, dau2_pt, dau2_eta, dau2_phi;
  vector<float> vx, vy, vz;

  for (auto const& trackVertex : *trackingVertexs.product()) {
    if (trackVertex.eventId().bunchCrossing() != 0) continue;  // Consider only in-time events
    
    for (TrackingVertex::tp_iterator source = trackVertex.sourceTracks_begin(); source != trackVertex.sourceTracks_end(); ++source) {
      auto decayTrk = source->get();
      if (decayTrk->pdgId() != 310 && decayTrk->pdgId() != 3122) continue;

      int count = 0;
      int GenHadFromQuark = 0;
      bool GenHadFromTop = false;

      motherTracking(decayTrk->pdgId(), trackVertex, decayTrk, count, GenHadFromQuark, GenHadFromTop);

      inVol.push_back(trackVertex.inVolume());
      isGenHadFromTop.push_back(GenHadFromTop);
      isGenHadFromTsb.push_back(GenHadFromQuark);

      vx.push_back(trackVertex.position().x());
      vy.push_back(trackVertex.position().y());
      vz.push_back(trackVertex.position().z());

      candidates->push_back(getCandidate(decayTrk));
      if (trackVertex.nDaughterTracks() >= 2) { 
        auto dau1 = trackVertex.daughterTracks().at(0).get();
        auto dau2 = trackVertex.daughterTracks().at(1).get();
        dau1_pdgId.push_back(dau1->pdgId());
        dau2_pdgId.push_back(dau2->pdgId());
        dau1_pt.push_back(dau1->pt());
        dau2_pt.push_back(dau2->pt());
        dau1_eta.push_back(dau1->eta());
        dau2_eta.push_back(dau2->eta());
        dau1_phi.push_back(dau1->phi());
        dau2_phi.push_back(dau2->phi());
      }
      else {
        dau1_pdgId.push_back(0);
        dau2_pdgId.push_back(0);
        dau1_pt.push_back(-99);
        dau2_pt.push_back(-99);
        dau1_eta.push_back(-99);
        dau2_eta.push_back(-99);
        dau1_phi.push_back(-99);
        dau2_phi.push_back(-99);
      }
    }
  }

  auto genHadTable = make_unique<nanoaod::FlatTable>(candidates->size(),"genHadron",false);
  genHadTable->addColumn<int>("isGenHadFromTsb",isGenHadFromTsb,"KS/Lam from t->s/b",nanoaod::FlatTable::IntColumn);
  genHadTable->addColumn<uint8_t>("isGenHadFromTop",isGenHadFromTop,"KS/Lam from top",nanoaod::FlatTable::UInt8Column);
  genHadTable->addColumn<uint8_t>("inVol",inVol,"track in volume",nanoaod::FlatTable::UInt8Column);
  genHadTable->addColumn<int>("dau1_pdgId", dau1_pdgId,"first daughter PID",nanoaod::FlatTable::IntColumn); 
  genHadTable->addColumn<int>("dau2_pdgId", dau2_pdgId,"second daughter PID",nanoaod::FlatTable::IntColumn); 
  genHadTable->addColumn<float>("dau1_pt", dau1_pt,"first daughter PT",nanoaod::FlatTable::FloatColumn); 
  genHadTable->addColumn<float>("dau2_pt", dau2_pt,"second daughter PT",nanoaod::FlatTable::FloatColumn); 
  genHadTable->addColumn<float>("dau1_eta", dau1_eta,"first daughter eta",nanoaod::FlatTable::FloatColumn); 
  genHadTable->addColumn<float>("dau2_eta", dau2_eta,"second daughter eta",nanoaod::FlatTable::FloatColumn); 
  genHadTable->addColumn<float>("dau1_phi", dau1_phi,"first daughter phi",nanoaod::FlatTable::FloatColumn); 
  genHadTable->addColumn<float>("dau2_phi", dau2_phi,"second daughter phi",nanoaod::FlatTable::FloatColumn); 
  genHadTable->addColumn<float>("vx", vx,"vertex x postion",nanoaod::FlatTable::FloatColumn);
  genHadTable->addColumn<float>("vy", vy,"vertex y postion",nanoaod::FlatTable::FloatColumn);
  genHadTable->addColumn<float>("vz", vz,"vertex z postion",nanoaod::FlatTable::FloatColumn);

  iEvent.put(move(genHadTable),"genHadron");
  iEvent.put(move(candidates));
}

int HadTruthProducer::trackingVertex_pdgId(const TrackingVertex* tv)
{
  for (TrackingVertex::tp_iterator source = tv->sourceTracks_begin(); source != tv->sourceTracks_end(); ++source) {
    return source->get()->pdgId();
  }  
  return 0;
}

const reco::GenParticleRef HadTruthProducer::getMother(const TrackingParticleRef& tp)
{
  const TrackingVertexRef& tv = tp->parentVertex();
  if (tv->nSourceTracks()){
    for (TrackingVertex::tp_iterator source = tv->sourceTracks_begin(); source != tv->sourceTracks_end(); ++source) {
      auto mothers = source->get()->genParticles();
      if (!mothers.empty()){
	reco::GenParticleRefVector::const_iterator im = mothers.begin();
	return *im;
      }
    }
  }

  if (!tp->genParticles().empty()){
    auto genpart = tp->genParticles()[0];
    const reco::GenParticleRefVector& mothers = genpart->motherRefVector();
    if (!mothers.empty()){
      reco::GenParticleRefVector::const_iterator im = mothers.begin();
      return *im;
    }
  }

  return reco::GenParticleRef();
}


bool HadTruthProducer::isGenHadFrom(const reco::GenParticle* particle, int pdgId, int count,int & GenHadFromQuark, bool & GenHadFromTop)
{
  GenHadFromTop = false;
  if (abs(particle->pdgId()) == pdgId && particle->status() == 62) {
    auto dau2 = particle->daughter(1);
    GenHadFromQuark = dau2->pdgId();
    GenHadFromTop = true;
    return true;
  }
  
  const reco::GenParticleRefVector& mothers = particle->motherRefVector(); 
  count = count + 1;
  for (reco::GenParticleRefVector::const_iterator im = mothers.begin(); im != mothers.end(); ++im) {
    const reco::GenParticle& part = **im;
    return isGenHadFrom( &part, pdgId, count, GenHadFromQuark, GenHadFromTop);
  }
  return false;
}

bool HadTruthProducer::isHadFrom(const reco::GenParticleRef& particle, int pdgId, int count, int & hadFromQuark, bool & hadFromTop)
{
  if(abs(particle->pdgId()) == pdgId && particle->status() == 62) {
    auto dau2 = particle->daughter(1);
    hadFromQuark = dau2->pdgId();
    return hadFromTop = true;
  }
  count = count + 1;
  for(unsigned int im = 0; im < particle->numberOfMothers(); ++im) {
    const reco::GenParticleRef& mothers = particle->motherRef(im);
    if( isHadFrom( mothers, pdgId, count, hadFromQuark, hadFromTop) ) {
      return hadFromTop = true;
    }
  }
  return hadFromTop = false;
}

void HadTruthProducer::motherTracking(int PID, const TrackingVertex trackVertex, const TrackingParticle *decayTrk, int count, int & GenHadFromQuark, bool & GenHadFromTop)
{
  if (!decayTrk->genParticles().empty()) {
    for (TrackingParticle::genp_iterator igen = decayTrk->genParticle_begin(); igen != decayTrk->genParticle_end(); ++igen) {
      auto gen = igen->get();
      if (count != 0 || decayTrk->pdgId() == PID) {
        isGenHadFrom(gen, 6, count, GenHadFromQuark, GenHadFromTop);
      }
    }
  }
  else {
    count = count + 1;
    auto pv = decayTrk->parentVertex().get();
    for (TrackingVertex::tp_iterator pr = pv->sourceTracks_begin(); pr != pv->sourceTracks_end(); ++pr) {
      auto decayTrk2 = pr->get();
      motherTracking(PID, *pv, decayTrk2, count, GenHadFromQuark, GenHadFromTop);
    }
  }
}

void HadTruthProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
