#include "nano/analysis/interface/hadAnalyser.h"

hadAnalyser::hadAnalyser(TTree *tree, TTree *had, TTree *hadTruth, Bool_t isMC, Bool_t dl, Bool_t sle, Bool_t slm) : topEventSelectionDL(tree, had, hadTruth, isMC, dl, sle, slm)
{ }

hadAnalyser::~hadAnalyser()
{ }

Double_t hadAnalyser::GetD(float pt, float eta, float phi, float m, float vx, float vy, float vz) {
  TLorentzVector tlv;
  tlv.SetPtEtaPhiM(pt,eta,phi,m);
  TVector3 dv(vx - PV_x, vy - PV_y, vz - PV_z);
  return tlv.Vect().Unit().Cross(dv).Mag();
}

void hadAnalyser::Loop(){}
