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
}


reco::VertexCompositeCandidate HadronProducer::fit(vector<reco::Candidate*>& cands,
						   reco::Vertex& pv, int pdgId,
						   float &dca, float &angleXY, float &angleXYZ)
{
  int charge = 0;
  vector<reco::TransientTrack> transientTracks;
  for (auto dau : cands){
    const reco::TransientTrack transientTrack = trackBuilder_->build(dau->bestTrack());
    transientTracks.emplace_back(transientTrack);
    //cout <<"no track ref "<<endl;
    charge += dau->charge();
  }

  if (transientTracks.size() < 2){
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
  
  KalmanVertexFitter m_kvf(true);
  
  TransientVertex tv = m_kvf.vertex(transientTracks);
  if (!tv.isValid()) return reco::VertexCompositeCandidate();
  
  reco::Vertex theVtx = tv;
  // loose cut on chi2
  if (theVtx.normalizedChi2() > vtxChi2Cut_) return reco::VertexCompositeCandidate();

  // if (theVtx.normalizedChi2() < 0)
  //   std::cout << "-ve vtx: " << theVtx.normalizedChi2() << " chi2/ndof" << theVtx.chi2() << " / " << theVtx.ndof()
  // 	      << ". trknc2 " << tt1.normalizedChi2() << " " << tt2.normalizedChi2() << " dca:" << dca << " " << tv.trackWeight(tt1) << " " << tv.trackWeight(tt2)
  // 	      << ". (x,y,z) " << cxPt.x() << ", " << cxPt.y() << ", " << cxPt.z() << " " << theVtx.x() << ", " << theVtx.y() << ", " << theVtx.z()
  // 	      << std::endl;
  // else
  //   std::cout << "+ve vtx: " << theVtx.normalizedChi2() << " chi2/ndof" << theVtx.chi2() << " / " << theVtx.ndof()
  // 	      << ". trknc2 " << tt1.normalizedChi2() << " " << tt2.normalizedChi2() << " dca:" << dca << " " << tv.trackWeight(tt1) << " " << tv.trackWeight(tt2)
  // 	      << ". (x,y,z) " << cxPt.x() << ", " << cxPt.y() << ", " << cxPt.z() << " " << theVtx.x() << ", " << theVtx.y() << ", " << theVtx.z()
  // 	      << std::endl;    
  
  GlobalPoint vtxPos(theVtx.x(), theVtx.y(), theVtx.z());

  math::XYZTLorentzVector tlv;
  GlobalVector totalP;
  int i = 0;
  for (auto trk : tv.refittedTracks()){
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
  if (angleXY > cosThetaXYCut_) return reco::VertexCompositeCandidate();
  
  // 3D pointing angle
  double dz = theVtx.z()-referencePos.z();
  double pz = totalP.z();
  angleXYZ = (dx*px+dy*py+dz*pz)/(sqrt(dx*dx+dy*dy+dz*dz)*sqrt(px*px+py*py+pz*pz));
  if (angleXYZ > cosThetaXYZCut_) return reco::VertexCompositeCandidate();

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
  
  vector<hadronCandidate> hadronCandidates;

  vector<reco::Candidate*> chargedHadrons, leptons;  
  for (auto & pfcand : *pfCandidates){
    
    if ( pfcand.charge() == 0 ) continue;
    //if ( pfcand.pt() < tkPtCut_ ) continue;
    if (pfcand.bestTrack() == nullptr) continue;
    
    reco::Candidate* recoDau = pfcand.clone();
    
    if ( abs(recoDau->pdgId()) == 11  || abs(recoDau->pdgId())==13)
      leptons.push_back(recoDau);
    else
      chargedHadrons.push_back(recoDau);        
  }
  // find KShort Cands
  auto KShortCands = findKShortCands(chargedHadrons, pv, -1);
  hadronCandidates.insert(hadronCandidates.end(), KShortCands.begin(), KShortCands.end());

  // find Lambda Cands
  auto LambdaCands = findLambdaCands(chargedHadrons, pv, -1);
  hadronCandidates.insert(hadronCandidates.end(), LambdaCands.begin(), LambdaCands.end());

  
  int njet = 0;
  for (const pat::Jet & aPatJet : *jetHandle){
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
    if (d0Cands.size()){
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
  vector<int> had_nJet, had_nDau;
  vector<float> had_jetDR, had_legDR, had_diffMass;
  vector<float> had_lxy, had_lxySig, had_l3D, had_l3DSig, had_dca, had_angleXY, had_angleXYZ;
  vector<float> had_dau1_chi2, had_dau1_nHits, had_dau1_pt, had_dau1_ipsigZ, had_dau1_ipsigXY;
  vector<float> had_dau2_chi2, had_dau2_nHits, had_dau2_pt, had_dau2_ipsigZ, had_dau2_ipsigXY;
  
  for (auto cand: hadronCandidates){
    had_cands->push_back(cand.vcc);
    had_jets->push_back(cand.jet);
    const reco::Track* dau1 = cand.vcc.daughter(0)->bestTrack();
    had_dau1_chi2.push_back(dau1->normalizedChi2());
    had_dau1_nHits.push_back(dau1->numberOfValidHits());
    had_dau1_pt.push_back(dau1->pt());
    had_dau1_ipsigZ.push_back(std::abs(dau1->dz(primaryVertexPoint)/dau1->dzError()));
    had_dau1_ipsigXY.push_back(std::abs(dau1->dxy(primaryVertexPoint)/dau1->dxyError()));

    const reco::Track* dau2 = cand.vcc.daughter(1)->bestTrack();
    had_dau2_chi2.push_back(dau2->normalizedChi2());
    had_dau2_nHits.push_back(dau2->numberOfValidHits());
    had_dau2_pt.push_back(dau2->pt());
    had_dau2_ipsigZ.push_back(std::abs(dau2->dz(primaryVertexPoint)/dau2->dzError()));
    had_dau2_ipsigXY.push_back(std::abs(dau2->dxy(primaryVertexPoint)/dau2->dxyError()));
    
    had_nJet.push_back(cand.nJet);
    had_nDau.push_back(cand.nDau);      
    
    had_jetDR.push_back(reco::deltaR( cand.vcc, cand.jet));
    had_legDR.push_back(reco::deltaR( *dau1, *dau2));
    
    had_diffMass.push_back(cand.diffMass);
    had_lxy.push_back(cand.lxy);
    had_lxySig.push_back(cand.lxySig);
    had_l3D.push_back(cand.l3D);
    had_l3DSig.push_back(cand.l3DSig);
    had_dca.push_back(cand.dca);
    had_angleXY.push_back(cand.angleXY);
    had_angleXYZ.push_back(cand.angleXYZ);
  }
  
  auto had_table = make_unique<nanoaod::FlatTable>(had_cands->size(),"had",false);
  had_table->addColumn<int>("nJet",had_nJet,"nJet of vertex cand",nanoaod::FlatTable::IntColumn);
  had_table->addColumn<int>("nDau",had_nDau,"nDau of vertex cand",nanoaod::FlatTable::IntColumn);

  had_table->addColumn<float>("jetDR",had_jetDR,"DR between jet",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("legDR",had_legDR,"DR between daugthers",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("diffMass",had_diffMass,"mass difference",nanoaod::FlatTable::FloatColumn);

  had_table->addColumn<float>("lxy",had_lxy,"2D decay length in cm",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("lxySig",had_lxySig,"2D decay length sig in cm",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("l3D",had_l3D,"3D decay length in cm",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("l3DSig",had_l3DSig,"3D decay length sig in cm",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("dca",had_dca,"distance of closest approach cm",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("angleXY",had_angleXY,"2D angle between vertex and tracks",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("angleXYZ",had_angleXYZ,"3D angle between vertex and tracks",nanoaod::FlatTable::FloatColumn);
  
  had_table->addColumn<float>("dau1_chi2",had_dau1_chi2,"dau1 chi2/ndof",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("dau1_nHits",had_dau1_nHits,"dau1 nHits",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("dau1_pt",had_dau1_pt,"dau1 Pt",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("dau1_ipsigXY",had_dau1_ipsigXY,"dau1 ipsigXY",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("dau1_ipsigZ",had_dau1_ipsigZ,"dau1 ipsigZ",nanoaod::FlatTable::FloatColumn);

  had_table->addColumn<float>("dau2_chi2",had_dau2_chi2,"dau2 chi2/ndof",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("dau2_nHits",had_dau2_nHits,"dau2 nHits",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("dau2_pt",had_dau2_pt,"dau2 Pt",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("dau2_ipsigXY",had_dau2_ipsigXY,"dau2 ipsigXY",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("dau2_ipsigZ",had_dau2_ipsigZ,"dau2 ipsigZ",nanoaod::FlatTable::FloatColumn);

  iEvent.put(move(had_table),"had");
  iEvent.put(move(had_cands));
  iEvent.put(move(had_jets),"jet");
}

vector<HadronProducer::hadronCandidate> HadronProducer::findJPsiCands(vector<reco::Candidate*> &leptons, reco::Vertex& pv, int nJet, const pat::Jet & aPatJet)
{
  vector<hadronCandidate> hadrons;
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
      if (abs(cand.mass() - jpsi_m_) > 0.2) continue;
      
      hc.vcc = cand;
      hc.jet = aPatJet;
      
      auto d2 = getDistance(2,cand,pv);
      hc.lxy = d2.first;      
      hc.lxySig = d2.second;
      auto d3 = getDistance(3,cand,pv);	
      hc.l3D = d3.first;
      hc.l3DSig = d3.second;

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
  vector<hadronCandidate> hadrons;
  // find jpsi to mumu or ee
  for (auto pion : chargedHads){
    //cout <<"pion pt = "<< pion->pt() << ", eta = "<< pion->eta() << ", pid = "<< pion->pdgId()<<endl;
    for (auto kaon : chargedHads){
      //cout <<"kaon pt = "<< kaon->pt() << ", eta = "<< kaon->eta() << ", pid = "<< kaon->pdgId()<<endl;
      if ( pion->charge() * kaon->charge() != -1 ) continue;

      pion->setMass(pion_m_);
      pion->setPdgId(pion->charge()*pion_pdgId_);
      kaon->setMass(kaon_m_);
      kaon->setPdgId(kaon->charge()*kaon_pdgId_);

      vector<reco::Candidate*> cands{pion, kaon};

      hadronCandidate hc;

      reco::VertexCompositeCandidate cand = fit(cands, pv, d0_pdgId_,
						hc.dca, hc.angleXY, hc.angleXYZ);

      if (cand.numberOfDaughters() < 2) continue;
      if (abs(cand.mass() - d0_m_) > 0.2) continue;
      
      hc.vcc = cand;
      hc.jet = aPatJet;
      
      auto d2 = getDistance(2,cand,pv);
      hc.lxy = d2.first;      
      hc.lxySig = d2.second;
      auto d3 = getDistance(3,cand,pv);	
      hc.l3D = d3.first;
      hc.l3DSig = d3.second;

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
  vector<hadronCandidate> hadrons;
  // find jpsi to mumu or ee
  for (auto d0 : d0cands){
    
    for (auto pion : chargedHads){

      const reco::Track* dau1 = d0.vcc.daughter(0)->bestTrack();
      if (dau1 == pion->bestTrack()) continue;
      
      const reco::Track* dau2 = d0.vcc.daughter(1)->bestTrack();
      if (dau2 == pion->bestTrack()) continue;

      pion->setMass(pion_m_);
      pion->setPdgId(pion->charge()*pion_pdgId_);
      
      // d0 first daughter should always be pion from findD0Cands
      if (abs(d0.vcc.daughter(0)->pdgId()) != pion_pdgId_){
	cout <<"HadronProducer::findDStarCands first daughter is not pion "<< d0.vcc.daughter(0)->pdgId() <<endl;
      }
      // D*+ -> [K- pi+]D0 pi+ (opposite signed kaon is suppressed by 2 OoM)
      // i.e. pions should be opposite charge
      if (!(d0.vcc.daughter(0)->pdgId() + pion->pdgId() == 0)) continue;
      

      hadronCandidate hc;

      vector<reco::Candidate*> cands{pion,
	  dynamic_cast<reco::Candidate*>(d0.vcc.daughter(0)),
	  dynamic_cast<reco::Candidate*>(d0.vcc.daughter(1))};
      //vector<reco::Candidate*> cands{pion, &d0.vcc};
      reco::VertexCompositeCandidate cand = fit(cands, pv, dstar_pdgId_,
						hc.dca, hc.angleXY, hc.angleXYZ);
      
      if (cand.numberOfDaughters() < 2) continue;
      
      float diffMass_Dstar = cand.mass() - d0.vcc.mass();      
      if (abs(diffMass_Dstar - (dstar_m_ - d0_m_)) > 0.2) continue;
      
      hc.vcc = cand;
      hc.jet = aPatJet;
      
      auto d2 = getDistance(2,cand,pv);
      hc.lxy = d2.first;      
      hc.lxySig = d2.second;
      auto d3 = getDistance(3,cand,pv);	
      hc.l3D = d3.first;
      hc.l3DSig = d3.second;

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
  vector<hadronCandidate> hadrons;
  for (auto pion1 : chargedHads){
    // avoid double counting pions by explicit charge finding
    if (pion1->charge() != +1) continue;
    for (auto pion2 : chargedHads){
      if (pion2->charge() != -1) continue;

      pion1->setMass(pion_m_);
      pion1->setPdgId(pion1->charge()*pion_pdgId_);
      pion2->setMass(pion_m_);
      pion2->setPdgId(pion2->charge()*pion_pdgId_);

      vector<reco::Candidate*> cands{pion1, pion2};

      hadronCandidate hc;

      reco::VertexCompositeCandidate cand = fit(cands, pv, kshort_pdgId_,
                                                hc.dca, hc.angleXY, hc.angleXYZ);

      if (cand.numberOfDaughters() < 2) continue;
      if (abs(cand.mass() - kshort_m_) > 0.2) continue;

      hc.vcc = cand;
      hc.jet = pat::Jet();

      auto d2 = getDistance(2,cand,pv);
      hc.lxy = d2.first;
      hc.lxySig = d2.second;
      auto d3 = getDistance(3,cand,pv);
      hc.l3D = d3.first;
      hc.l3DSig = d3.second;

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
  vector<hadronCandidate> hadrons;
  for (auto proton : chargedHads){
    for (auto pion : chargedHads){
      if ( proton->charge() * pion->charge() != -1 ) continue;

      proton->setMass(proton_m_);
      proton->setPdgId(proton->charge()*proton_pdgId_);
      pion->setMass(pion_m_);
      pion->setPdgId(pion->charge()*pion_pdgId_);

      vector<reco::Candidate*> cands{proton, pion};

      hadronCandidate hc;

      reco::VertexCompositeCandidate cand = fit(cands, pv, lambda_pdgId_,
                                                hc.dca, hc.angleXY, hc.angleXYZ);

      if (cand.numberOfDaughters() < 2) continue;
      if (abs(cand.mass() - lambda_m_) > 0.2) continue;

      hc.vcc = cand;
      hc.jet = pat::Jet();

      auto d2 = getDistance(2,cand,pv);
      hc.lxy = d2.first;
      hc.lxySig = d2.second;
      auto d3 = getDistance(3,cand,pv);
      hc.l3D = d3.first;
      hc.l3DSig = d3.second;

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
  vector<hadronCandidate> hadrons;
  // find jpsi to mumu or ee
  for (auto lambda :LambdaCands){
    for (auto jpsi : jpsiCands){
      
      if (lambda.vcc.pdgId()!=lambda_pdgId_ && lambda.vcc.pdgId()!=-lambda_pdgId_) continue;
      if (lambda.vcc.pdgId() * jpsi.vcc.pdgId() != 0) continue;

      hadronCandidate hc;

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
      reco::VertexCompositeCandidate cand(0, tlv, jpsi.vcc.vertex(), jpsi.vcc.vertexCovariance(), jpsi.vcc.vertexChi2(), jpsi.vcc.vertexNdof(), lambdab_pdgId_);
      if (cand.numberOfDaughters() < 2) continue;
      if (abs(cand.mass() - lambdab_m_) > 0.2) continue;
      
      hc.vcc = cand;
      hc.jet = aPatJet;
      
      auto d2 = getDistance(2,cand,pv);
      hc.lxy = d2.first;
      hc.lxySig = d2.second;
      auto d3 = getDistance(3,cand,pv);	
      hc.l3D = d3.first;
      hc.l3DSig = d3.second;

      hc.nJet = nJet;
      hc.nDau = 2;
      hc.diffMass = -9;

      hadrons.push_back(hc);
    }
  }
  return hadrons;
}
