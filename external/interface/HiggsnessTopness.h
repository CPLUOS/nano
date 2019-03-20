#include <TLorentzVector.h>
#include "TMinuit.h"
#include "TError.h"

using namespace std;

// global variables for topness higgsness
const double mtop = 173.;
const double mh = 125;
const double mw = 80.419;
const double mhw = 30;

const double sigmahlep = 2.;
const double sigmaWon = 5.;
const double sigmaWoff = 5.;
const double sigmat = 5.;

const int n_unknown_par = 4;

const double minuit_pmin = -5000;
const double minuit_pmax = 5000;

TLorentzVector g_missing, g_lepton1, g_lepton2, g_bottom1, g_bottom2;

TMinuit *ptMinuit;
TMinuit *ptMinuitT1;
TMinuit *ptMinuitT2;
Double_t minuit_arglist[10];
Int_t ierflag = 0;

double GetMass(double Px, double Py, double Pz, double E) {
  double f = E*E - Px*Px - Py*Py - Pz*Pz;
  double fout = 0;
  if ( f < 0. ) {
    fout = 0;
  } else {
    fout = sqrt(f);
  }
  return fout;
}

void fcnT1(  int& npar, double* deriv, double& f, double par[], int flag) {
  double chisq = 0;
  double Nu1Px = par[0];
  double Nu1Py = par[1];
  double Nu1Pz = par[2];
  double Nu1E = sqrt( Nu1Px*Nu1Px + Nu1Py*Nu1Py + Nu1Pz*Nu1Pz );

  double Nu2Px = g_missing.Px() - Nu1Px;
  double Nu2Py = g_missing.Py() - Nu1Py;
  double Nu2Pz = par[3];
  double Nu2E = sqrt( Nu2Px*Nu2Px + Nu2Py*Nu2Py + Nu2Pz*Nu2Pz );

  double Lep1E = g_lepton1.E();
  double Lep1Px = g_lepton1.Px();
  double Lep1Py = g_lepton1.Py();
  double Lep1Pz = g_lepton1.Pz();

  double Lep2E = g_lepton2.E();
  double Lep2Px = g_lepton2.Px();
  double Lep2Py = g_lepton2.Py();
  double Lep2Pz = g_lepton2.Pz();

  double Jet1E = g_bottom1.E();
  double Jet1Px = g_bottom1.Px();
  double Jet1Py = g_bottom1.Py();
  double Jet1Pz = g_bottom1.Pz();

  double Jet2E = g_bottom2.E();
  double Jet2Px = g_bottom2.Px();
  double Jet2Py = g_bottom2.Py();
  double Jet2Pz = g_bottom2.Pz();

  double jln111Px = Jet1Px + Lep1Px + Nu1Px;
  double jln111Py = Jet1Py + Lep1Py + Nu1Py;
  double jln111Pz = Jet1Pz + Lep1Pz + Nu1Pz;
  double jln111E = Jet1E + Lep1E + Nu1E;

  double mjln111 = GetMass( jln111Px, jln111Py, jln111Pz, jln111E );

  double deltaJLN111 = (mjln111*mjln111 - mtop*mtop) / (sigmat*sigmat);

  double jln222Px = Jet2Px + Lep2Px + Nu2Px;
  double jln222Py = Jet2Py + Lep2Py + Nu2Py;
  double jln222Pz = Jet2Pz + Lep2Pz + Nu2Pz;
  double jln222E = Jet2E + Lep2E + Nu2E;

  double mjln222 = GetMass( jln222Px, jln222Py, jln222Pz, jln222E );

  double deltaJLN222 = (mjln222*mjln222 - mtop*mtop) / (sigmat*sigmat);

  double ln11Px = Lep1Px + Nu1Px;
  double ln11Py = Lep1Py + Nu1Py;
  double ln11Pz = Lep1Pz + Nu1Pz;
  double ln11E = Lep1E + Nu1E;

  double mln11 = GetMass( ln11Px, ln11Py, ln11Pz, ln11E );
  double deltaLN11on = (mln11*mln11 - mw*mw) / (sigmaWon*sigmaWon);

  double ln22Px = Lep2Px + Nu2Px;
  double ln22Py = Lep2Py + Nu2Py;
  double ln22Pz = Lep2Pz + Nu2Pz;
  double ln22E = Lep2E + Nu2E;

  double mln22 = GetMass( ln22Px, ln22Py, ln22Pz, ln22E );

  double deltaLN22on = (mln22*mln22 - mw*mw) / (sigmaWon*sigmaWon);

  chisq = deltaJLN111*deltaJLN111 + deltaJLN222*deltaJLN222 + deltaLN11on*deltaLN11on + deltaLN22on*deltaLN22on;

  f = chisq;
}

void fcnT2(  int& npar, double* deriv, double& f, double par[], int flag) {
  double chisq = 0;
  double Nu1Px = par[0];
  double Nu1Py = par[1];
  double Nu1Pz = par[2];
  double Nu1E = sqrt( Nu1Px*Nu1Px + Nu1Py*Nu1Py + Nu1Pz*Nu1Pz );

  double Nu2Px = g_missing.Px() - Nu1Px;
  double Nu2Py = g_missing.Py() - Nu1Py;
  double Nu2Pz = par[3];
  double Nu2E = sqrt( Nu2Px*Nu2Px + Nu2Py*Nu2Py + Nu2Pz*Nu2Pz );

  double Lep1E = g_lepton1.E();
  double Lep1Px = g_lepton1.Px();
  double Lep1Py = g_lepton1.Py();
  double Lep1Pz = g_lepton1.Pz();

  double Lep2E = g_lepton2.E();
  double Lep2Px = g_lepton2.Px();
  double Lep2Py = g_lepton2.Py();
  double Lep2Pz = g_lepton2.Pz();

  double Jet1E = g_bottom1.E();
  double Jet1Px = g_bottom1.Px();
  double Jet1Py = g_bottom1.Py();
  double Jet1Pz = g_bottom1.Pz();

  double Jet2E = g_bottom2.E();
  double Jet2Px = g_bottom2.Px();
  double Jet2Py = g_bottom2.Py();
  double Jet2Pz = g_bottom2.Pz();

  double jln121Px = Jet1Px + Lep2Px + Nu1Px;
  double jln121Py = Jet1Py + Lep2Py + Nu1Py;
  double jln121Pz = Jet1Pz + Lep2Pz + Nu1Pz;
  double jln121E = Jet1E + Lep2E + Nu1E;

  double mjln121 = GetMass( jln121Px, jln121Py, jln121Pz, jln121E );

  double deltaJLN121 = (mjln121*mjln121 - mtop*mtop) / (sigmat*sigmat);

  double jln212Px = Jet2Px + Lep1Px + Nu2Px;
  double jln212Py = Jet2Py + Lep1Py + Nu2Py;
  double jln212Pz = Jet2Pz + Lep1Pz + Nu2Pz;
  double jln212E = Jet2E + Lep1E + Nu2E;

  double mjln212 = GetMass( jln212Px, jln212Py, jln212Pz, jln212E );

  double deltaJLN212 = (mjln212*mjln212 - mtop*mtop) / (sigmat*sigmat);

  double ln21Px = Lep2Px + Nu1Px;
  double ln21Py = Lep2Py + Nu1Py;
  double ln21Pz = Lep2Pz + Nu1Pz;
  double ln21E = Lep2E + Nu1E;

  double mln21 = GetMass( ln21Px, ln21Py, ln21Pz, ln21E );
  double deltaLN21on = (mln21*mln21 - mw*mw) / (sigmaWon*sigmaWon);

  double ln12Px = Lep1Px + Nu2Px;
  double ln12Py = Lep1Py + Nu2Py;
  double ln12Pz = Lep1Pz + Nu2Pz;
  double ln12E = Lep1E + Nu2E;

  double mln12 = GetMass( ln12Px, ln12Py, ln12Pz, ln12E );

  double deltaLN12on = (mln12*mln12 - mw*mw) / (sigmaWon*sigmaWon);

  chisq = deltaJLN121*deltaJLN121 + deltaJLN212*deltaJLN212 + deltaLN21on*deltaLN21on + deltaLN12on*deltaLN12on;

  f = chisq;
}

void fcnH(  int& npar, double* deriv, double& f, double par[], int flag) {
  double chisq = 0;
  double Nu1Px = par[0];
  double Nu1Py = par[1];
  double Nu1Pz = par[2];
  double Nu1E = sqrt( Nu1Px*Nu1Px + Nu1Py*Nu1Py + Nu1Pz*Nu1Pz );

  double Nu2Px = g_missing.Px() - Nu1Px;
  double Nu2Py = g_missing.Py() - Nu1Py;
  double Nu2Pz = par[3];
  double Nu2E = sqrt( Nu2Px*Nu2Px + Nu2Py*Nu2Py + Nu2Pz*Nu2Pz );

  double Lep1E = g_lepton1.E();
  double Lep1Px = g_lepton1.Px();
  double Lep1Py = g_lepton1.Py();
  double Lep1Pz = g_lepton1.Pz();

  double Lep2E = g_lepton2.E();
  double Lep2Px = g_lepton2.Px();
  double Lep2Py = g_lepton2.Py();
  double Lep2Pz = g_lepton2.Pz();

  double llnnPx = Lep1Px + Lep2Px + Nu1Px + Nu2Px;
  double llnnPy = Lep1Py + Lep2Py + Nu1Py + Nu2Py;
  double llnnPz = Lep1Pz + Lep2Pz + Nu1Pz + Nu2Pz;
  double llnnE = Lep1E + Lep2E + Nu1E + Nu2E;

  double mllnn = GetMass( llnnPx, llnnPy, llnnPz, llnnE );

  double deltaLLNN = (mllnn*mllnn - mh*mh) / (sigmahlep*sigmahlep);

  double ln11Px = Lep1Px + Nu1Px;
  double ln11Py = Lep1Py + Nu1Py;
  double ln11Pz = Lep1Pz + Nu1Pz;
  double ln11E = Lep1E + Nu1E;

  double mln11 = GetMass( ln11Px, ln11Py, ln11Pz, ln11E );

  double deltaLN11on = (mln11*mln11 - mw*mw) / (sigmaWon*sigmaWon);

  double deltaLN11off = (mln11*mln11 - mhw*mhw) / (sigmaWoff*sigmaWoff);

  double ln22Px = Lep2Px + Nu2Px;
  double ln22Py = Lep2Py + Nu2Py;
  double ln22Pz = Lep2Pz + Nu2Pz;
  double ln22E = Lep2E + Nu2E;

  double mln22 = GetMass( ln22Px, ln22Py, ln22Pz, ln22E );

  double deltaLN22on = (mln22*mln22 - mw*mw) / (sigmaWon*sigmaWon);

  double deltaLN22off = (mln22*mln22 - mhw*mhw) / (sigmaWoff*sigmaWoff);

  double ln12Px = Lep1Px + Nu2Px;
  double ln12Py = Lep1Py + Nu2Py;
  double ln12Pz = Lep1Pz + Nu2Pz;
  double ln12E = Lep1E + Nu2E;

  double mln12 = GetMass( ln12Px, ln12Py, ln12Pz, ln12E );

  double deltaLN12on = (mln12*mln12 - mw*mw) / (sigmaWoff*sigmaWoff);
  double deltaLN12off = (mln12*mln12 - mhw*mhw) / (sigmaWoff*sigmaWoff);

  double ln21Px = Lep2Px + Nu1Px;
  double ln21Py = Lep2Py + Nu1Py;
  double ln21Pz = Lep2Pz + Nu2Pz;
  double ln21E = Lep2E + Nu1E;

  double mln21 = GetMass( ln21Px, ln21Py, ln21Pz, ln21E );

  double deltaLN21on = (mln21*mln21 - mw*mw) / (sigmaWon*sigmaWon);
  double deltaLN21off = (mln21*mln21 - mhw*mhw) / (sigmaWoff*sigmaWoff);

 double  deltaCase1 = deltaLN11on*deltaLN11on + deltaLN22off*deltaLN22off;
 double  deltaCase2 = deltaLN12on*deltaLN12on + deltaLN21off*deltaLN21off;
 double  deltaCase3 = deltaLN22on*deltaLN22on + deltaLN11off*deltaLN11off;
 double  deltaCase4 = deltaLN21on*deltaLN21on + deltaLN12off*deltaLN12off;

  double deltaCaseTotal[] = {deltaCase1, deltaCase2, deltaCase3, deltaCase4};
  double deltaCaseMin = *min_element(deltaCaseTotal, deltaCaseTotal+4);

  //double nnPx = Nu1Px + Nu2Px;
  //double nnPy = Nu1Py + Nu2Py;
  //double nnPz = Nu1Pz + Nu2Pz;
  //double nnE = Nu1E + Nu2E;

  //double mnn = GetMass( nnPx, nnPy, nnPz, nnE );

  chisq = deltaLLNN*deltaLLNN + deltaCaseMin;

  f = chisq;
}

float GetHiggsness() {
  float higgsness = 0;
  ptMinuit = new TMinuit(n_unknown_par);
  ptMinuit->SetPrintLevel(-1);
  ptMinuit->SetFCN(fcnH);
  minuit_arglist[0] = 1;
  ptMinuit->mnexcm("SET ERR", minuit_arglist, 1, ierflag);
  Double_t vstart[4] = {0.0, 0.0, 0.0, 0.0};
  Double_t step[4] = {0.1, 0.1, 0.1, 0.1};
  TString par_name[4] = {"Nu1Px", "Nu1Py", "Nu1Pz", "Nu2Pz"};
  for (int i=0; i<4; i++) {
    ptMinuit->mnparm(i, par_name[i], vstart[i], step[i], minuit_pmin, minuit_pmax, ierflag);
  }

  minuit_arglist[0] = 500;
  minuit_arglist[1] = 1.;

  ptMinuit->mnexcm("MIGRAD", minuit_arglist, 2, ierflag);

  double outparH[n_unknown_par], errH[n_unknown_par];
  for (int i = 0; i<4; i++) {
    ptMinuit->GetParameter(i, outparH[i], errH[i]);
  }

  Double_t aminH, edmH, errdefH;
  Int_t nvparH, nparxH, icstatH;

  ptMinuit->mnstat(aminH, edmH, errdefH, nvparH, nparxH, icstatH);
  if ( aminH > 0 ) {
    higgsness = log10(aminH);
  }
  delete ptMinuit;
  return higgsness;
}

float GetTopness() {
  float topness = 0;
  ptMinuitT1 = new TMinuit(n_unknown_par);
  ptMinuitT1->SetPrintLevel(-1);
  ptMinuitT1->SetFCN(fcnT1);
  minuit_arglist[0] = 1;
  ptMinuitT1->mnexcm("SET ERR", minuit_arglist, 1, ierflag);
  Double_t vstart[4] = {0.0, 0.0, 0.0, 0.0};
  Double_t step[4] = {0.1, 0.1, 0.1, 0.1};
  TString par_name[4] = {"Nu1Px", "Nu1Py", "Nu1Pz", "Nu2Pz"};
  for (int i=0; i<4; i++) {
    ptMinuitT1->mnparm(i, par_name[i], vstart[i], step[i], minuit_pmin, minuit_pmax, ierflag);
  }

  minuit_arglist[0] = 500;
  minuit_arglist[1] = 1.;

  ptMinuitT1->mnexcm("MIGRAD", minuit_arglist, 2, ierflag);

  double outparT1[n_unknown_par], errT1[n_unknown_par];
  for (int i = 0; i<4; i++) {
    ptMinuitT1->GetParameter(i, outparT1[i], errT1[i]);
  }

  Double_t aminT1, edmT1, errdefT1;
  Int_t nvparT1, nparxT1, icstatT1;

  ptMinuitT1->mnstat(aminT1, edmT1, errdefT1, nvparT1, nparxT1, icstatT1);

  delete ptMinuitT1;

  ptMinuitT2 = new TMinuit(n_unknown_par);
  ptMinuitT2->SetPrintLevel(-1);
  ptMinuitT2->SetFCN(fcnT2);
  minuit_arglist[0] = 1;
  ptMinuitT2->mnexcm("SET ERR", minuit_arglist, 1, ierflag);
  //vstart[4] = {0.0, 0.0, 0.0, 0.0};
  //step[4] = {0.1, 0.1, 0.1, 0.1};
  //par_name[4] = {"Nu1Px", "Nu1Py", "Nu1Pz", "Nu2Pz"};
  for (int i=0; i<4; i++) {
    ptMinuitT2->mnparm(i, par_name[i], vstart[i], step[i], minuit_pmin, minuit_pmax, ierflag);
  }

  minuit_arglist[0] = 500;
  minuit_arglist[1] = 1.;

  ptMinuitT2->mnexcm("MIGRAD", minuit_arglist, 2, ierflag);

  double outparT2[n_unknown_par], errT2[n_unknown_par];
  for (int i = 0; i<4; i++) {
    ptMinuitT2->GetParameter(i, outparT2[i], errT2[i]);
  }

  Double_t aminT2, edmT2, errdefT2;
  Int_t nvparT2, nparxT2, icstatT2;

  ptMinuitT2->mnstat(aminT2, edmT2, errdefT2, nvparT2, nparxT2, icstatT2);

  delete ptMinuitT2;

  Double_t aminT;
  aminT = min(aminT1, aminT2);

  if( aminT > 0 )
  {
    topness = log10(aminT);
  }
  return topness;
}
