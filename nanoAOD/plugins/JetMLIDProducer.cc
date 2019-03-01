#include "JetMLIDProducer.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <cmath>

//#define debugMode
using namespace edm;
using namespace std;

JetMLIDProducer::JetMLIDProducer(const edm::ParameterSet & iConfig) :
  jetLabel_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jetLabel"))),
  vertexLabel_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexLabel"))),
  useQC( iConfig.getParameter<bool>("useQualityCuts"))
{
  produces<nanoaod::FlatTable>("jetID");
  produces<nanoaod::FlatTable>("jetDau");
}

void JetMLIDProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void
JetMLIDProducer::produce( edm::Event& iEvent, const edm::EventSetup& iSetup)
{    
  Handle<reco::VertexCollection> recVtxs;
  iEvent.getByToken(vertexLabel_,recVtxs);
  reco::Vertex pv = recVtxs->at(0);

  Handle<edm::View<pat::Jet> > jets;
  iEvent.getByToken(jetLabel_, jets);

  // saving all variables
  vector<float> jet_delta, jet_axis2, jet_axis1, jet_ptD,
    jet_cpt1, jet_cpt2, jet_cpt3,
    jet_npt1, jet_npt2, jet_npt3,
    jet_JetAngularity, jet_GeoMoment, jet_HalfPtMoment, jet_DRSquareMoment,
    jet_SmallDRPT, jet_MassMoment, jet_PTSquare, jet_MyMoment;

  vector<int> jet_cmult, jet_nmult;
  vector<int> jet_dauIdx1, jet_dauIdx2;
  vector<float> jetDau_pt, jetDau_eta, jetDau_phi, jetDau_E, jetDau_charge;

  int nDaughter = 0;
  for (auto jet = jets->begin();  jet != jets->end(); ++jet) {

    float sum_weight = 0., sum_deta = 0., sum_dphi = 0., sum_deta2 = 0., sum_dphi2 = 0., sum_detadphi = 0., sum_pt = 0.;
    int mult = 0, nmult = 0, cmult = 0;
    float pt_dr_log = 0;

    float JetAngularity = 0., GeoMoment = 0., HalfPtMoment = 0., DRSquareMoment = 0., SmallDRPT = 0., MassMoment = 0., PTSquare = 0., MyMoment = 0.;
    int ParticleCount = 0;

    std::vector<float> chargedPt;
    std::vector<float> neutralPt;

    jet_dauIdx1.push_back(nDaughter);
    //Loop over the jet constituents
    for(auto daughter : jet->getJetConstituentsQuick()){
      //packed candidate situation
      auto part = static_cast<const pat::PackedCandidate*>(daughter);

      if(part->charge()){
	chargedPt.push_back(part->pt());
	if(!(part->fromPV() > 1 && part->trackHighPurity()))
	  continue;
	if(useQC){
	  if((part->dz()*part->dz())/(part->dzError()*part->dzError()) > 25.)
	    continue;
	  if((part->dxy()*part->dxy())/(part->dxyError()*part->dxyError()) < 25.){
	    ++mult;
	    ++cmult;
	  }
	}
	else {
	  ++mult;
	  ++cmult;
	}
      }
      else {
	if(part->pt() < 1.0)
	  continue;
	neutralPt.push_back(part->pt());
	++mult;
	++nmult;
      }

      //Calculate pt_dr_log and some vars                                                                                                                                
      float dr = reco::deltaR(*jet, *part);
      float theta = (M_PI*dr)/(2*0.4);

      pt_dr_log += std::log(part->pt()/dr);

      JetAngularity += (daughter->pt()*daughter->pt())*(powf((std::sin(theta)),-2))*((1-powf((std::cos(theta)),-3)))/(jet->mass()); // a = -2 
      GeoMoment += (powf((part->pt()/jet->pt()),1))*(powf((dr/0.4),1));
      HalfPtMoment += (powf((part->pt()/jet->pt()),1.5))*(powf((dr/0.4),0));
      DRSquareMoment += (powf((part->pt()/jet->pt()),0))*(powf((dr/0.4),2));
      if(dr < 0.1)
	SmallDRPT += (powf((part->pt()/jet->pt()),1))*(powf((dr/0.4),0));
      MassMoment += (powf((part->pt()/jet->pt()),1))*(powf((dr/0.4),2));
      ParticleCount += (powf((part->pt()/jet->pt()),0))*(powf((dr/0.4),0));
      PTSquare += (powf((part->pt()/jet->pt()),2))*(powf((dr/0.4),0));
      MyMoment += (powf((part->pt()/jet->pt()),2))*(powf((dr/0.4),-1.5));
      
      float deta   = daughter->eta() - jet->eta();
      float dphi   = reco::deltaPhi(daughter->phi(), jet->phi());
      float partPt = daughter->pt();
      float weight = partPt*partPt;

      sum_weight   += weight;
      sum_pt       += partPt;
      sum_deta     += deta*weight;
      sum_dphi     += dphi*weight;
      sum_deta2    += deta*deta*weight;
      sum_detadphi += deta*dphi*weight;
      sum_dphi2    += dphi*dphi*weight;

      nDaughter += 1;
      jetDau_pt.push_back(daughter->pt());
      jetDau_eta.push_back(daughter->eta());
      jetDau_phi.push_back(daughter->phi());
      jetDau_E.push_back(daughter->energy());
      jetDau_charge.push_back(daughter->charge());
    }
    jet_dauIdx2.push_back(nDaughter);

    //Calculate axis2 and ptD
    float a = 0., b = 0., c = 0.;
    float ave_deta = 0., ave_dphi = 0., ave_deta2 = 0., ave_dphi2 = 0.;
    if(sum_weight > 0){
      ave_deta  = sum_deta/sum_weight;
      ave_dphi  = sum_dphi/sum_weight;
      ave_deta2 = sum_deta2/sum_weight;
      ave_dphi2 = sum_dphi2/sum_weight;
      a         = ave_deta2 - ave_deta*ave_deta;                          
      b         = ave_dphi2 - ave_dphi*ave_dphi;                          
      c         = -(sum_detadphi/sum_weight - ave_deta*ave_dphi);                
    }
    float delta = sqrt(fabs((a-b)*(a-b)+4*c*c));
    float axis2 = (a+b-delta > 0 ?  sqrt(0.5*(a+b-delta)) : 0);
    float axis1 = (a+b+delta > 0 ?  sqrt(0.5*(a+b+delta)) : 0);
    float ptD   = (sum_weight > 0 ? sqrt(sum_weight)/sum_pt : 0);

    jet_delta.push_back(delta);
    jet_axis2.push_back(axis2);
    jet_axis1.push_back(axis1);
    jet_ptD.push_back(ptD);

    jet_nmult.push_back(nmult);
    jet_cmult.push_back(cmult);

    sort(chargedPt.begin(), chargedPt.end(), [](float a, float b) {return (a > b);});
    sort(neutralPt.begin(), neutralPt.end(), [](float a, float b) {return (a > b);});
    
    while (chargedPt.size() < 3) chargedPt.push_back(0.);
    jet_cpt1.push_back(chargedPt[0]);
    jet_cpt2.push_back(chargedPt[0] + chargedPt[1]);
    jet_cpt3.push_back(chargedPt[0] + chargedPt[1] + chargedPt[2]);
    
    while (neutralPt.size() < 3) neutralPt.push_back(0.);
    jet_npt1.push_back(neutralPt[0]);
    jet_npt2.push_back(neutralPt[0] + neutralPt[1]);
    jet_npt3.push_back(neutralPt[0] + neutralPt[1] + neutralPt[2]);

    jet_JetAngularity.push_back(JetAngularity); 
    jet_GeoMoment.push_back(GeoMoment);
    jet_HalfPtMoment.push_back(HalfPtMoment);
    jet_DRSquareMoment.push_back(DRSquareMoment);
    jet_SmallDRPT.push_back(SmallDRPT);
    jet_MassMoment.push_back(MassMoment);
    jet_PTSquare.push_back(PTSquare);
    jet_MyMoment.push_back(MyMoment);
  }
  
  auto jetID_table = make_unique<nanoaod::FlatTable>(jets->size(),"jetID",false);
  jetID_table->addColumn<float>("delta",jet_delta,"delta",nanoaod::FlatTable::FloatColumn);
  jetID_table->addColumn<float>("axis2",jet_axis2,"axis2",nanoaod::FlatTable::FloatColumn);
  jetID_table->addColumn<float>("axis1",jet_axis1,"axis1",nanoaod::FlatTable::FloatColumn);
  jetID_table->addColumn<float>("ptD",jet_ptD,"ptD",nanoaod::FlatTable::FloatColumn);
  jetID_table->addColumn<int>("dauIdx1",jet_dauIdx1,"start index in daughter table",nanoaod::FlatTable::IntColumn);
  jetID_table->addColumn<int>("dauIdx2",jet_dauIdx2,"end index in daughter table",nanoaod::FlatTable::IntColumn);

  jetID_table->addColumn<float>("angularity",    jet_JetAngularity,  "angularity",     nanoaod::FlatTable::FloatColumn);
  jetID_table->addColumn<float>("geoMoment",     jet_GeoMoment,      "geoMoment",      nanoaod::FlatTable::FloatColumn);
  jetID_table->addColumn<float>("halfPtMoment",  jet_HalfPtMoment,   "halfPtMoment",   nanoaod::FlatTable::FloatColumn);
  jetID_table->addColumn<float>("dRSquareMoment",jet_DRSquareMoment, "dRSquareMoment", nanoaod::FlatTable::FloatColumn);
  jetID_table->addColumn<float>("smallDRPT",     jet_SmallDRPT,      "smallDRPT",      nanoaod::FlatTable::FloatColumn);
  jetID_table->addColumn<float>("massMoment",    jet_MassMoment,     "massMoment",     nanoaod::FlatTable::FloatColumn);
  jetID_table->addColumn<float>("pTSquare",      jet_PTSquare,       "pTSquare",       nanoaod::FlatTable::FloatColumn);
  jetID_table->addColumn<float>("myMoment",      jet_MyMoment,       "myMoment",       nanoaod::FlatTable::FloatColumn);

  jetID_table->addColumn<float>("cpt1",jet_cpt1,"cpt1",nanoaod::FlatTable::FloatColumn);
  jetID_table->addColumn<float>("cpt2",jet_cpt2,"cpt2",nanoaod::FlatTable::FloatColumn);
  jetID_table->addColumn<float>("cpt3",jet_cpt3,"cpt3",nanoaod::FlatTable::FloatColumn);
  jetID_table->addColumn<int>("cmult",jet_cmult,"cmult",nanoaod::FlatTable::IntColumn);

  jetID_table->addColumn<float>("npt1",jet_npt1,"npt1",nanoaod::FlatTable::FloatColumn);
  jetID_table->addColumn<float>("npt2",jet_npt2,"npt2",nanoaod::FlatTable::FloatColumn);
  jetID_table->addColumn<float>("npt3",jet_npt3,"npt3",nanoaod::FlatTable::FloatColumn);
  jetID_table->addColumn<int>("nmult",jet_nmult,"nmult",nanoaod::FlatTable::IntColumn);

  auto jetDau_table = make_unique<nanoaod::FlatTable>(nDaughter,"jetDau",false);
  jetDau_table->addColumn<float>("pt",jetDau_pt,"pt",nanoaod::FlatTable::FloatColumn);
  jetDau_table->addColumn<float>("eta",jetDau_eta,"eta",nanoaod::FlatTable::FloatColumn);
  jetDau_table->addColumn<float>("phi",jetDau_phi,"phi",nanoaod::FlatTable::FloatColumn);
  jetDau_table->addColumn<float>("E",jetDau_E,"E",nanoaod::FlatTable::FloatColumn);
  jetDau_table->addColumn<float>("charge",jetDau_charge,"charge",nanoaod::FlatTable::FloatColumn);

  iEvent.put(move(jetID_table),"jetID");
  iEvent.put(move(jetDau_table),"jetDau");
}

DEFINE_FWK_MODULE(JetMLIDProducer);
