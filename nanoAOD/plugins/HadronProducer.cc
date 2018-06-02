#include "HadronProducer.h"
//#define debugMode
using namespace edm;
using namespace std;

HadronProducer::HadronProducer(const edm::ParameterSet & iConfig) :
  jetLabel_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jetLabel"))),
  vertexLabel_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexLabel"))),
  pfCandidates_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("pfCandLabel")))
{
  // cuts on initial track selection
  tkChi2Cut_ = iConfig.getParameter<double>("tkChi2Cut");
  tkNHitsCut_ = iConfig.getParameter<int>("tkNHitsCut");
  tkPtCut_ = iConfig.getParameter<double>("tkPtCut");
  tkIPSigXYCut_ = iConfig.getParameter<double>("tkIPSigXYCut");
  tkIPSigZCut_ = iConfig.getParameter<double>("tkIPSigZCut");   
  // cuts on vertex
  vtxChi2Cut_ = iConfig.getParameter<double>("vtxChi2Cut");
  vtxDecaySigXYZCut_ = iConfig.getParameter<double>("vtxDecaySigXYZCut");
  vtxDecaySigXYCut_ = iConfig.getParameter<double>("vtxDecaySigXYCut");
  // miscellaneous cuts
  tkDCACut_ = iConfig.getParameter<double>("tkDCACut");
  cosThetaXYCut_ = iConfig.getParameter<double>("cosThetaXYCut");
  cosThetaXYZCut_ = iConfig.getParameter<double>("cosThetaXYZCut");
  
  produces<nanoaod::FlatTable>("had");
  produces<reco::VertexCompositeCandidateCollection>();
  produces<vector<pat::Jet>>("jet");
  produces<vector<vector<int>>>("index");
}

reco::VertexCompositeCandidate HadronProducer::fit(vector<reco::Candidate*>& cands,
						   reco::Vertex& pv, int pdgId,
						   float &dca, float &angleXY, float &angleXYZ)
{
  int charge = 0;
  vector<reco::TransientTrack> transientTracks;
  for (auto &dau : cands) {
    const reco::TransientTrack transientTrack = trackBuilder_->build(dau->bestTrack());
    transientTracks.emplace_back(transientTrack);
    //cout <<"no track ref "<<endl;
    charge += dau->charge();
  }

  if (transientTracks.size() < 2) {
    //cout <<"no tracks... something is wrong"<<endl;
    return reco::VertexCompositeCandidate();
  }

  // impactPointTSCP DCA
  // measure distance between tracks at their closest approach
  dca = -1;
  ClosestApproachInRPhi cApp;
  reco::TransientTrack tt1 = transientTracks[0], tt2 = transientTracks[1];
  cApp.calculate(tt1.impactPointTSCP().theState(), tt2.impactPointTSCP().theState());  
  if (!cApp.status()) return reco::VertexCompositeCandidate();
  
  dca = cApp.distance();

  if (dca > tkDCACut_) return reco::VertexCompositeCandidate();
  
  // the POCA should at least be in the sensitive volume
  GlobalPoint cxPt = cApp.crossingPoint();
  if (sqrt(cxPt.x()*cxPt.x() + cxPt.y()*cxPt.y()) > 120. || abs(cxPt.z()) > 300.)
    return reco::VertexCompositeCandidate();

  // the tracks should at least point in the same quadrant
  TrajectoryStateClosestToPoint const & posTSCP = tt1.trajectoryStateClosestToPoint(cxPt);
  TrajectoryStateClosestToPoint const & negTSCP = tt2.trajectoryStateClosestToPoint(cxPt);
  if (!posTSCP.isValid() || !negTSCP.isValid()) return reco::VertexCompositeCandidate();
  if (posTSCP.momentum().dot(negTSCP.momentum()) < 0) return reco::VertexCompositeCandidate();
  
  static KalmanVertexFitter m_kvf(true);
  
  TransientVertex tv = m_kvf.vertex(transientTracks);
  if (!tv.isValid()) return reco::VertexCompositeCandidate();
  
  reco::Vertex theVtx = tv;
  // loose cut on chi2
  if (theVtx.normalizedChi2() > vtxChi2Cut_) return reco::VertexCompositeCandidate();
  
  GlobalPoint vtxPos(theVtx.x(), theVtx.y(), theVtx.z());

  math::XYZTLorentzVector tlv;
  GlobalVector totalP;
  int i = 0;
  for (auto trk : tv.refittedTracks()) {
    TrajectoryStateClosestToPoint const & tscp = trk.trajectoryStateClosestToPoint(vtxPos);
    GlobalVector mom = tscp.momentum();
    double mass = cands[i]->mass();
    double energy = sqrt(mom.mag2() + mass*mass);    
    const math::XYZTLorentzVector lv(mom.x(), mom.y(), mom.z(), energy);
    totalP += mom;
    tlv += lv;
    ++i;
  }

  math::XYZPoint referencePos = pv.position();

  // 2D pointing angle
  double dx = theVtx.x()-referencePos.x();
  double dy = theVtx.y()-referencePos.y();
  double px = totalP.x();
  double py = totalP.y();
  angleXY = (dx*px+dy*py)/(sqrt(dx*dx+dy*dy)*sqrt(px*px+py*py));
  if (angleXY < cosThetaXYCut_) return reco::VertexCompositeCandidate();
  
  // 3D pointing angle
  double dz = theVtx.z()-referencePos.z();
  double pz = totalP.z();
  angleXYZ = (dx*px+dy*py+dz*pz)/(sqrt(dx*dx+dy*dy+dz*dz)*sqrt(px*px+py*py+pz*pz));
  if (angleXYZ < cosThetaXYZCut_) return reco::VertexCompositeCandidate();

  reco::Particle::Point vtx(theVtx.x(), theVtx.y(), theVtx.z());
  const reco::Vertex::CovarianceMatrix vtxCov(theVtx.covariance());

  reco::VertexCompositeCandidate secVert(charge, tlv, vtx, vtxCov, theVtx.chi2(), theVtx.ndof(), pdgId);
  for (auto dau : cands){
    const reco::PFCandidate *pfDau = dynamic_cast<const reco::PFCandidate*>(&(*dau));
    if (pfDau){
      auto recoDau = make_unique<reco::RecoChargedCandidate>(dau->charge(), dau->p4(), dau->vertex(), dau->pdgId());
      recoDau->setTrack(pfDau->trackRef());
      secVert.addDaughter(*recoDau);
    }
    else {
      secVert.addDaughter(*dau);
    }
  }
  auto sigXYcheck = getDistance(2,secVert,pv);
  if (sigXYcheck.first/sigXYcheck.second < vtxDecaySigXYCut_) return reco::VertexCompositeCandidate();
  auto sigXYZcheck = getDistance(3,secVert,pv);
  if (sigXYZcheck.first/sigXYZcheck.second < vtxDecaySigXYZCut_) return reco::VertexCompositeCandidate();

  return secVert;
}

HadronProducer::SVector3 HadronProducer::getDistanceVector(int dim, reco::VertexCompositeCandidate& vertex,reco::Vertex& pv)
{
  float z = 0.;
  if (dim == 3) z = vertex.vz() - pv.position().z();
  SVector3 distanceVector(vertex.vx() - pv.position().x(),
			  vertex.vy() - pv.position().y(),
			  z);
  return distanceVector;
}

pair<float, float> HadronProducer::getDistance(int dim, reco::VertexCompositeCandidate& vertex,reco::Vertex& pv)
{
  SMatrixSym3D totalCov = vertex.vertexCovariance() + pv.covariance();
  SVector3 distVecXYZ = getDistanceVector(dim, vertex, pv);
  float distMagXYZ = ROOT::Math::Mag(distVecXYZ);
  float sigmaDistMagXYZ = sqrt(ROOT::Math::Similarity(totalCov, distVecXYZ)) / distMagXYZ;

  return make_pair(distMagXYZ,sigmaDistMagXYZ);
}

void HadronProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void
HadronProducer::produce( edm::Event& iEvent, const edm::EventSetup& iSetup)
{  
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",trackBuilder_);
  
  Handle<reco::VertexCollection> recVtxs;
  iEvent.getByToken(vertexLabel_,recVtxs);
  reco::Vertex pv = recVtxs->at(0);
  math::XYZPoint primaryVertexPoint = pv.position();

  Handle<edm::View<pat::Jet> > jetHandle;
  iEvent.getByToken(jetLabel_, jetHandle);

  Handle<reco::CandidateView> pfCandidates;
  iEvent.getByToken(pfCandidates_, pfCandidates);

  hadronCandidateCollection hadronCandidates;

  vector<reco::Candidate*> chargedHadrons, leptons;  
  for (auto & pfcand : *pfCandidates){
    
    if (pfcand.charge() == 0) continue;
    //if ( pfcand.pt() < tkPtCut_ ) continue;
    const reco::Track * trk = pfcand.bestTrack();      
    if (trk == nullptr) continue;

    float ipsigXY = std::abs(trk->dxy(primaryVertexPoint)/trk->dxyError());
    if (ipsigXY < 2.0) continue;
    reco::Candidate* recoDau = pfcand.clone();
    
    if (abs(recoDau->pdgId()) == 11 || abs(recoDau->pdgId()) == 13)
      leptons.push_back(recoDau);
    else
      chargedHadrons.push_back(recoDau);        
  }

  // float old_XY = cosThetaXYCut_;
  // cosThetaXYCut_ = 0.9;
  
  // find KShort Cands
  auto KShortCands = findKShortCands(chargedHadrons, pv, -1);
  hadronCandidates.insert(hadronCandidates.end(), KShortCands.begin(), KShortCands.end());

  // find Lambda Cands
  auto LambdaCands = findLambdaCands(chargedHadrons, pv, -1);
  hadronCandidates.insert(hadronCandidates.end(), LambdaCands.begin(), LambdaCands.end());

  // cosThetaXYCut_ = old_XY;
  
  for (auto lep : leptons) delete lep;
  for (auto pion : chargedHadrons) delete pion;
  
  int njet = 0;
  for (const pat::Jet & aPatJet : *jetHandle) {
    if (aPatJet.pt() < 30 or abs(aPatJet.eta()) > 3 ) continue;

    chargedHadrons.clear(); leptons.clear();
    
    for( unsigned int idx = 0 ; idx < aPatJet.numberOfDaughters() ; ++idx) {
      auto dau = aPatJet.daughter(idx);

      if ( dau->charge() == 0 ) continue;      
      //if ( dau->pt() < tkPtCut_ ) continue;
      if (dau->bestTrack() == nullptr) continue;
      
      //const reco::Track * trk = dau->bestTrack();      
      // if (trk->normalizedChi2() > tkChi2Cut_) continue;
      // if (trk->numberOfValidHits() < tkNHitsCut_) continue;

      // float ipsigXY = std::abs(trk->dxy(primaryVertexPoint)/trk->dxyError());
      // //if (ipsigXY > tkIPSigXYCut_) continue;      
      // float ipsigZ = std::abs(trk->dz(primaryVertexPoint)/trk->dzError());
      //if (ipsigZ > tkIPSigZCut_) continue;
            
#ifdef debugMode
      cout <<"dau pt = "<< dau->pt() << ", eta = "<< dau->eta() << ", pid = "<< dau->pdgId()<<endl;
#endif

      reco::Candidate* recoDau = dau->clone();
      
      if ( abs(recoDau->pdgId()) == 11  || abs(recoDau->pdgId())==13)
	leptons.push_back(recoDau);
      else
	chargedHadrons.push_back(recoDau);

    }
    unsigned int dau_size = chargedHadrons.size() + leptons.size();
    
    if ( dau_size < 2 ) continue;
    
    // find JPsi Cands 
    auto jpsiCands = findJPsiCands(leptons, pv, njet, aPatJet);      
    hadronCandidates.insert(hadronCandidates.end(), jpsiCands.begin(), jpsiCands.end());
    
    // find D0 Cands
    auto d0Cands = findD0Cands(chargedHadrons, pv, njet, aPatJet);      
    hadronCandidates.insert(hadronCandidates.end(), d0Cands.begin(), d0Cands.end());

    // find dstar cands
    if (d0Cands.size()) {
      auto dStarCands = findDStarCands(d0Cands, chargedHadrons, pv, njet, aPatJet);
      hadronCandidates.insert(hadronCandidates.end(), dStarCands.begin(), dStarCands.end());	
    }
    
    // find LambdaB Cands
    auto LambdaBCands = findLambdaBCands(LambdaCands,jpsiCands, pv, njet, aPatJet);
    hadronCandidates.insert(hadronCandidates.end(), LambdaBCands.begin(), LambdaBCands.end());
    
    ++njet;
    for (auto lep : leptons) delete lep;
    for (auto pion : chargedHadrons) delete pion;
  }

  // saving all variables
  auto had_cands = make_unique<reco::VertexCompositeCandidateCollection>();
  auto had_jets = make_unique<vector<pat::Jet>>();
  auto had_indices = make_unique<vector<vector<int>>>();
  vector<int> had_nJet, had_nDau;
  vector<float> had_jetDR, had_legDR, had_diffMass;
  vector<float> had_lxy, had_lxyErr, had_l3D, had_l3DErr, had_dca, had_angleXY, had_angleXYZ;
  vector<float> had_dau1_chi2, had_dau1_nHits, had_dau1_pt, had_dau1_ipsigZ, had_dau1_ipsigXY;
  vector<float> had_dau2_chi2, had_dau2_nHits, had_dau2_pt, had_dau2_ipsigZ, had_dau2_ipsigXY;
  vector<int> had_idx, had_dau1_idx, had_dau2_idx, had_dau1_charge, had_dau2_charge;

  for (auto cand: hadronCandidates){
    had_cands->push_back(cand.vcc);
    had_jets->push_back(cand.jet);
    if (abs(cand.vcc.pdgId()) != lambdab_pdgId_) {
      const reco::Track* dau1 = cand.vcc.daughter(0)->bestTrack();
      if (dau1) {
	had_dau1_chi2.push_back(dau1->normalizedChi2());
	had_dau1_nHits.push_back(dau1->numberOfValidHits());
	had_dau1_pt.push_back(dau1->pt());
	had_dau1_ipsigZ.push_back(std::abs(dau1->dz(primaryVertexPoint)/dau1->dzError()));
	had_dau1_ipsigXY.push_back(std::abs(dau1->dxy(primaryVertexPoint)/dau1->dxyError()));
	had_dau1_charge.push_back(dau1->charge());
      } else {
	had_dau1_chi2.push_back(0);
	had_dau1_nHits.push_back(0);
	had_dau1_pt.push_back(0);
	had_dau1_ipsigZ.push_back(0);
	had_dau1_ipsigXY.push_back(0);
	had_dau1_charge.push_back(0);
      }
      
      const reco::Track* dau2 = cand.vcc.daughter(1)->bestTrack();
      if (dau2) {
	had_dau2_chi2.push_back(dau2->normalizedChi2());
	had_dau2_nHits.push_back(dau2->numberOfValidHits());
	had_dau2_pt.push_back(dau2->pt());
	had_dau2_ipsigZ.push_back(std::abs(dau2->dz(primaryVertexPoint)/dau2->dzError()));
	had_dau2_ipsigXY.push_back(std::abs(dau2->dxy(primaryVertexPoint)/dau2->dxyError()));
	had_dau2_charge.push_back(dau2->charge());
      } else {
	had_dau2_chi2.push_back(0);
	had_dau2_nHits.push_back(0);
	had_dau2_pt.push_back(0);
	had_dau2_ipsigZ.push_back(0);
	had_dau2_ipsigXY.push_back(0);
	had_dau2_charge.push_back(0);
      }
      if (dau1 && dau2) {
	had_legDR.push_back(reco::deltaR(*dau1, *dau2));
      } else {
	had_legDR.push_back(0.0);
      }
    } else {
	had_dau1_chi2.push_back(0);
	had_dau1_nHits.push_back(0);
	had_dau1_pt.push_back(0);
	had_dau1_ipsigZ.push_back(0);
	had_dau1_ipsigXY.push_back(0);
	had_dau1_charge.push_back(0);
	had_dau2_chi2.push_back(0);
	had_dau2_nHits.push_back(0);
	had_dau2_pt.push_back(0);
	had_dau2_ipsigZ.push_back(0);
	had_dau2_ipsigXY.push_back(0);
	had_dau2_charge.push_back(0);
	had_legDR.push_back(0.0);
    }
    
    had_nJet.push_back(cand.nJet);
    had_nDau.push_back(cand.nDau);      
    
    had_jetDR.push_back(reco::deltaR(cand.vcc, cand.jet));
    
    had_diffMass.push_back(cand.diffMass);
    had_lxy.push_back(cand.lxy);
    had_lxyErr.push_back(cand.lxyErr);
    had_l3D.push_back(cand.l3D);
    had_l3DErr.push_back(cand.l3DErr);
    had_dca.push_back(cand.dca);
    had_angleXY.push_back(cand.angleXY);
    had_angleXYZ.push_back(cand.angleXYZ);
    had_idx.push_back(cand.idx);
    had_dau1_idx.push_back(cand.dau1_idx);
    had_dau2_idx.push_back(cand.dau2_idx);
    
    vector<int> had_index;
    had_index.push_back(cand.idx);
    had_index.push_back(cand.dau1_idx);
    had_index.push_back(cand.dau2_idx);
    had_indices->push_back(had_index);
  }
  
  auto had_table = make_unique<nanoaod::FlatTable>(had_cands->size(),"had",false);
  had_table->addColumn<int>("nJet",had_nJet,"nJet of vertex cand",nanoaod::FlatTable::IntColumn);
  had_table->addColumn<int>("nDau",had_nDau,"nDau of vertex cand",nanoaod::FlatTable::IntColumn);

  had_table->addColumn<float>("jetDR",had_jetDR,"DR between jet",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("legDR",had_legDR,"DR between daugthers",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("diffMass",had_diffMass,"mass difference",nanoaod::FlatTable::FloatColumn);

  had_table->addColumn<float>("lxy",had_lxy,"2D decay length in cm",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("lxyErr",had_lxyErr,"2D decay length sigma in cm",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("l3D",had_l3D,"3D decay length in cm",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("l3DErr",had_l3DErr,"3D decay length sigma in cm",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("dca",had_dca,"distance of closest approach cm",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("angleXY",had_angleXY,"2D angle between vertex and tracks",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("angleXYZ",had_angleXYZ,"3D angle between vertex and tracks",nanoaod::FlatTable::FloatColumn);
  
  had_table->addColumn<float>("dau1_chi2",had_dau1_chi2,"dau1 chi2/ndof",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("dau1_nHits",had_dau1_nHits,"dau1 nHits",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("dau1_pt",had_dau1_pt,"dau1 Pt",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("dau1_ipsigXY",had_dau1_ipsigXY,"dau1 ipsigXY",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("dau1_ipsigZ",had_dau1_ipsigZ,"dau1 ipsigZ",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<int>("dau1_charge",had_dau1_charge,"dau1 charge",nanoaod::FlatTable::IntColumn);

  had_table->addColumn<float>("dau2_chi2",had_dau2_chi2,"dau2 chi2/ndof",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("dau2_nHits",had_dau2_nHits,"dau2 nHits",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("dau2_pt",had_dau2_pt,"dau2 Pt",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("dau2_ipsigXY",had_dau2_ipsigXY,"dau2 ipsigXY",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("dau2_ipsigZ",had_dau2_ipsigZ,"dau2 ipsigZ",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<int>("dau2_charge",had_dau2_charge,"dau2 charge",nanoaod::FlatTable::IntColumn);

  had_table->addColumn<int>("idx",had_idx,"index of itself",nanoaod::FlatTable::IntColumn);
  had_table->addColumn<int>("dau1_idx",had_dau1_idx,"index of dau1",nanoaod::FlatTable::IntColumn);
  had_table->addColumn<int>("dau2_idx",had_dau2_idx,"index of dau2",nanoaod::FlatTable::IntColumn);
  
  iEvent.put(move(had_table),"had");
  iEvent.put(move(had_cands));
  iEvent.put(move(had_jets),"jet");
  iEvent.put(move(had_indices),"index");
}

vector<HadronProducer::hadronCandidate> HadronProducer::findJPsiCands(vector<reco::Candidate*> &leptons, reco::Vertex& pv, int nJet, const pat::Jet & aPatJet)
{
  hadronCandidateCollection hadrons;
  // find jpsi to mumu or ee
  for (auto lep1Cand : leptons){
    if (lep1Cand->pdgId() > 0) continue; 
    for (auto lep2Cand : leptons){
      if (lep2Cand->pdgId() < 0) continue; 

      int pdgMul = lep1Cand->pdgId() * lep2Cand->pdgId();
      if ( pdgMul != -121 and pdgMul != -169 ) continue; 

      vector<reco::Candidate*> cands{lep1Cand, lep2Cand};

      hadronCandidate hc;

      reco::VertexCompositeCandidate cand = fit(cands, pv, jpsi_pdgId_,
						hc.dca, hc.angleXY, hc.angleXYZ);

      if (cand.numberOfDaughters() < 2) continue;
      if (fabs(cand.mass() - jpsi_m_) > 0.3) continue;
      
      hc.vcc = cand;
      hc.jet = aPatJet;
      
      auto d2 = getDistance(2,cand,pv);
      hc.lxy = d2.first;      
      hc.lxyErr = d2.second;
      auto d3 = getDistance(3,cand,pv);	
      hc.l3D = d3.first;
      hc.l3DErr = d3.second;

      hc.nJet = nJet;
      hc.nDau = 2;
      hc.diffMass = -9;

      hadrons.push_back(hc);
    }
  }
  return hadrons;
}

vector<HadronProducer::hadronCandidate> HadronProducer::findD0Cands(vector<reco::Candidate*> &chargedHads, reco::Vertex& pv, int nJet, const pat::Jet & aPatJet)
{
  hadronCandidateCollection hadrons;
  // find jpsi to mumu or ee
  for (auto &pion : chargedHads) {
    for (auto &kaon : chargedHads) {
      if (pion->charge() * kaon->charge() != -1) continue;

      pion->setMass(pion_m_);
      kaon->setMass(kaon_m_);

      vector<reco::Candidate*> cands{pion, kaon};

      hadronCandidate hc;

      // D0 -> K-pi+, D0bar -> K+pi-
      reco::VertexCompositeCandidate cand = fit(cands, pv, -kaon->charge()*d0_pdgId_,
						hc.dca, hc.angleXY, hc.angleXYZ);

      if (cand.numberOfDaughters() < 2) continue;
      if (fabs(cand.mass() - d0_m_) > 0.2) continue;
      
      hc.vcc = cand;
      hc.jet = aPatJet;
      
      auto d2 = getDistance(2,cand,pv);
      hc.lxy = d2.first;      
      hc.lxyErr = d2.second;
      auto d3 = getDistance(3,cand,pv);	
      hc.l3D = d3.first;
      hc.l3DErr = d3.second;

      hc.nJet = nJet;
      hc.nDau = 2;
      hc.diffMass = -9;

      hadrons.push_back(hc);
    }
  }
  return hadrons;
}

vector<HadronProducer::hadronCandidate> HadronProducer::findDStarCands(vector<HadronProducer::hadronCandidate>& d0cands, vector<reco::Candidate*> &chargedHads,
								       reco::Vertex& pv, int nJet, const pat::Jet & aPatJet)
{
  hadronCandidateCollection hadrons;
  // find jpsi to mumu or ee
  for (auto &d0 : d0cands) {
    for (auto &pion : chargedHads) {

      const reco::Track* dau1 = d0.vcc.daughter(0)->bestTrack();
      if (dau1 == pion->bestTrack()) continue;
      
      const reco::Track* dau2 = d0.vcc.daughter(1)->bestTrack();
      if (dau2 == pion->bestTrack()) continue;

      pion->setMass(pion_m_);
      
      // d0 first daughter should always be pion from findD0Cands
      if (fabs(d0.vcc.daughter(0)->mass() - pion_m_) > 0.0001) {
	cout <<"HadronProducer::findDStarCands first daughter is not pion "<< d0.vcc.daughter(0)->mass() << " " << pion_m_ << endl;
      }
      // D*+ -> [K- pi+]D0 pi+ (opposite signed kaon is suppressed by 2 OoM)
      // i.e. pions should be same charge
      if (d0.vcc.daughter(0)->pdgId() != pion->pdgId()) continue;
      

      hadronCandidate hc;

      vector<reco::Candidate*> cands{pion,
	  dynamic_cast<reco::Candidate*>(d0.vcc.daughter(0)),
	  dynamic_cast<reco::Candidate*>(d0.vcc.daughter(1))};
      //vector<reco::Candidate*> cands{pion, &d0.vcc};
      reco::VertexCompositeCandidate cand = fit(cands, pv, pion->charge()*dstar_pdgId_,
						hc.dca, hc.angleXY, hc.angleXYZ);
      
      if (cand.numberOfDaughters() < 2) continue;
      
      float diffMass_Dstar = cand.mass() - d0.vcc.mass();      
      if (fabs(diffMass_Dstar - (dstar_m_ - d0_m_)) > 0.2) continue;
      
      hc.vcc = cand;
      hc.jet = aPatJet;
      
      auto d2 = getDistance(2,cand,pv);
      hc.lxy = d2.first;      
      hc.lxyErr = d2.second;
      auto d3 = getDistance(3,cand,pv);	
      hc.l3D = d3.first;
      hc.l3DErr = d3.second;

      hc.nJet = nJet;
      hc.nDau = 3;
      hc.diffMass = diffMass_Dstar;

      hadrons.push_back(hc);
    }
  }
  
  return hadrons;
}

vector<HadronProducer::hadronCandidate> HadronProducer::findKShortCands(vector<reco::Candidate*> &chargedHads,
									reco::Vertex& pv, int nJet)
{
  hadronCandidateCollection hadrons;
  for (auto &pion1 : chargedHads) {
    // avoid double counting pions by explicit charge finding
    if (pion1->charge() != +1) continue;
    for (auto &pion2 : chargedHads) {
      if (pion2->charge() != -1) continue;

      pion1->setMass(pion_m_);
      pion2->setMass(pion_m_);

      if (fabs(kshort_m_ - (pion1->p4() + pion2->p4()).M()) > 1.0) continue;

      vector<reco::Candidate*> cands{pion1, pion2};

      hadronCandidate hc;

      reco::VertexCompositeCandidate cand = fit(cands, pv, kshort_pdgId_,
                                                hc.dca, hc.angleXY, hc.angleXYZ);

      if (cand.numberOfDaughters() < 2) continue;
      if (fabs(cand.mass() - kshort_m_) > 0.2) continue;

      hc.vcc = cand;
      hc.jet = pat::Jet();

      auto d2 = getDistance(2,cand,pv);
      hc.lxy = d2.first;
      hc.lxyErr = d2.second;
      auto d3 = getDistance(3,cand,pv);
      hc.l3D = d3.first;
      hc.l3DErr = d3.second;

      hc.nJet = nJet;
      hc.nDau = 2;
      hc.diffMass = -9;

      hadrons.push_back(hc);
    }
  }
  return hadrons;
}

vector<HadronProducer::hadronCandidate> HadronProducer::findLambdaCands(vector<reco::Candidate*> &chargedHads,
									reco::Vertex& pv, int nJet)
{
  hadronCandidateCollection hadrons;
  for (auto proton : chargedHads) {

    proton->setMass(proton_m_);
    if (!proton->bestTrack() || proton->bestTrack()->pt() < 0.4) continue;

    for (auto pion : chargedHads) {
      if ( proton->charge() * pion->charge() != -1 ) continue;

      pion->setMass(pion_m_);

      if (fabs(lambda_m_ - (proton->p4() + pion->p4()).M()) > 0.8) continue;
      
      vector<reco::Candidate*> cands{proton, pion};

      hadronCandidate hc;

      reco::VertexCompositeCandidate cand = fit(cands, pv,
						proton->charge()*lambda_pdgId_,
                                                hc.dca, hc.angleXY, hc.angleXYZ);

      if (cand.numberOfDaughters() < 2) continue;
      if (fabs(cand.mass() - lambda_m_) > 0.1) continue;

      hc.vcc = cand;
      hc.jet = pat::Jet();

      auto d2 = getDistance(2,cand,pv);
      hc.lxy = d2.first;
      hc.lxyErr = d2.second;
      auto d3 = getDistance(3,cand,pv);
      hc.l3D = d3.first;
      hc.l3DErr = d3.second;

      hc.nJet = nJet;
      hc.nDau = 2;
      hc.diffMass = -9;

      hadrons.push_back(hc);
    }
  }
  return hadrons;
}

vector<HadronProducer::hadronCandidate> HadronProducer::findLambdaBCands(vector<HadronProducer::hadronCandidate>& LambdaCands,vector<HadronProducer::hadronCandidate>& jpsiCands, reco::Vertex& pv, int nJet, const pat::Jet & aPatJet)
{
  hadronCandidateCollection hadrons;
  // find jpsi to mumu or ee
  for (auto lambda :LambdaCands) {
    for (auto jpsi : jpsiCands) {
      // if (lambda.vcc.pdgId()!=lambda_pdgId_ && lambda.vcc.pdgId()!=-lambda_pdgId_) continue;
      // if (lambda.vcc.pdgId() * jpsi.vcc.pdgId() != 0) continue;

      ///// TODO Can't vertex Lambda daughter with jpsi daughters due
      ///// to Lambda flight, need to write a fit function that can
      ///// vertex lambda (not a RecoChargedCandidate!) and jpsi
      ///// daughters together
      
      // vector<reco::Candidate*> cands{
      // 	  dynamic_cast<reco::Candidate*>(lambda.vcc.daughter(0)),
      // 	  dynamic_cast<reco::Candidate*>(lambda.vcc.daughter(1)),
      // 	  dynamic_cast<reco::Candidate*>(jpsi.vcc.daughter(0)),
      // 	  dynamic_cast<reco::Candidate*>(jpsi.vcc.daughter(1))
      // 	  };
      // reco::VertexCompositeCandidate cand = fit(cands, pv, lambdab_pdgId_,
      // 						hc.dca, hc.angleXY, hc.angleXYZ);

      math::XYZTLorentzVector tlv;
      tlv += lambda.vcc.p4();
      tlv += jpsi.vcc.p4();
      
      reco::VertexCompositeCandidate cand(0, tlv, jpsi.vcc.vertex(), jpsi.vcc.vertexCovariance(), jpsi.vcc.vertexChi2(), jpsi.vcc.vertexNdof(), (lambda.vcc.pdgId() > 0) ? lambdab_pdgId_ : -lambdab_pdgId_);

      // if (cand.numberOfDaughters() < 2) continue;
      if (fabs(cand.mass() - lambdab_m_) > 0.4) continue;
     
      cand.addDaughter(lambda.vcc);
      cand.addDaughter(jpsi.vcc);          
      
      hadronCandidate hc;
      
      hc.vcc = cand;
      hc.jet = aPatJet;
      
      auto d2 = getDistance(2,cand,pv);
      hc.lxy = d2.first;
      hc.lxyErr = d2.second;
      auto d3 = getDistance(3,cand,pv);	
      hc.l3D = d3.first;
      hc.l3DErr = d3.second;

      hc.nJet = nJet;
      hc.nDau = 2;
      hc.diffMass = -9;
     
      hc.dau1_idx = lambda.idx; 
      hc.dau2_idx = jpsi.idx;
      hadrons.push_back(hc);
    }
  }
  return hadrons;
}

DEFINE_FWK_MODULE(HadronProducer);
