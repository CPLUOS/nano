#include "nano/analysis/interface/hadAnalysis.h"

hadAnalysis::hadAnalysis(TTree *tree, Bool_t isMC, Bool_t dl, Bool_t sle, Bool_t slm) : dilepTopAnalysis(tree, isMC, dl, sle, slm)
{ }

hadAnalysis::~hadAnalysis()
{ }

Double_t hadAnalysis::GetD(float pt, float eta, float phi, float m, float vx, float vy, float vz) {
  TLorentzVector tlv;
  tlv.SetPtEtaPhiM(pt,eta,phi,m);
  TVector3 dv(vx - PV_x, vy - PV_y, vz - PV_z);
  return tlv.Vect().Unit().Cross(dv).Mag();
}

void hadAnalysis::Loop(){}
