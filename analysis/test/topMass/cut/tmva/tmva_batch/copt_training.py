import ROOT, os, getopt, sys, array, math, glob, json
from ROOT import *
from array import array


## Make TTree ##
FileArg = sys.argv
tempdir = FileArg[1]
#Dirname = "/cms/scratch/seulgi/nanoAOD/src/nano/analysis/test/topMass/cut/batch/Results/%s/" %tempdir
Dirname = "/cms/scratch/yyoun/nanoAOD/src/nano/analysis/test/topMass/cut/tmva/tmva_batch/Results/copt_927/%s/" %tempdir
if not os.path.isdir(Dirname):
    os.makedirs(Dirname)

temp = FileArg[2].split('/').pop()
cattree = Dirname+temp

f = ROOT.TFile(cattree, "recreate")
Hadsig = ROOT.TTree("Hadsig", "Hadsig")
Hadbkg = ROOT.TTree("Hadbkg", "Hadbkg")
D0sig = ROOT.TTree("D0sig", "D0sig")
D0bkg = ROOT.TTree("D0bkg", "D0bkg")

## Variables ##
cme_dca = array("f",[0])
cme_angleXY = array("f",[0])
cme_angleXYZ = array("f",[0])
cme_lxy = array("f",[0])
cme_lxyErr = array("f",[0])
cme_lxyE = array("f",[0])
cme_l3D = array("f",[0])
cme_l3DErr = array("f",[0])
cme_l3DE = array("f",[0])
cme_jetDR = array("f",[0])
cme_legDR = array("f",[0])
cme_diffMass = array("f",[0])
cme_nJet = array("i",[0])
cmeTruth_nMatched = array("i",[0])
cmeTruth_nTrueDau = array("i",[0])
cme_chi2 = array("f",[0])
cme_eta = array("f",[0])
cme_phi = array("f",[0])
cme_mass = array("f",[0])
cme_pt = array("f",[0])
cme_x = array("f",[0])
cme_y = array("f",[0])
cme_z = array("f",[0])
cme_ndof = array("i",[0])
cme_pdgId = array("i",[0])

cme_jet_btagCMVA = array("f", [0])
cme_jet_btagCSVV2 = array("f", [0])
cme_jet_btagDeepB = array("f", [0])
cme_jet_btagDeepC = array("f", [0])
cme_nDau = array("i",[0])
cme_dau1_chi2 = array("f", [0])
cme_dau1_ipsigXY = array("f", [0])
cme_dau1_ipsigZ = array("f", [0])
cme_dau1_nHits = array("f", [0])
cme_dau1_pt = array("f", [0])
cme_dau2_chi2 = array("f", [0])
cme_dau2_ipsigXY = array("f", [0])
cme_dau2_ipsigZ = array("f", [0])
cme_dau2_nHits = array("f", [0])
cme_dau2_pt = array("f", [0])

# Branch ##
Hadsig.Branch("cme_dca", cme_dca, "cme_dca/F")
Hadsig.Branch("cme_diffMass", cme_diffMass, "cme_diffMass/F")
Hadsig.Branch("cme_angleXY", cme_angleXY, "cme_angleXY/F")
Hadsig.Branch("cme_angleXYZ", cme_angleXYZ, "cme_angleXYZ/F")
Hadsig.Branch("cme_lxy", cme_lxy, "cme_lxy/F")
Hadsig.Branch("cme_lxyErr", cme_lxyErr, "cme_lxyErr/F")
Hadsig.Branch("cme_lxyE", cme_lxyE, "cme_lxyE/F")
Hadsig.Branch("cme_l3D", cme_l3D, "cme_l3D/F")
Hadsig.Branch("cme_l3DErr", cme_l3DErr, "cme_l3DErr/F")
Hadsig.Branch("cme_l3DE", cme_l3DE, "cme_l3DE/F")
Hadsig.Branch("cme_jetDR", cme_jetDR, "cme_jetDR/F")
Hadsig.Branch("cme_legDR", cme_legDR, "cme_legDR/F")
Hadsig.Branch("cme_nJet", cme_nJet, "cme_nJet/I")
Hadsig.Branch("cme_chi2", cme_chi2, "cme_chi2/F")
Hadsig.Branch("cme_eta", cme_eta, "cme_eta/F")
Hadsig.Branch("cme_mass", cme_mass, "cme_mass/F")
Hadsig.Branch("cme_phi", cme_phi, "cme_phi/F")
Hadsig.Branch("cme_pt", cme_pt, "cme_pt/F")
Hadsig.Branch("cme_x", cme_x, "cme_x/F")
Hadsig.Branch("cme_y", cme_y, "cme_y/F")
Hadsig.Branch("cme_z", cme_z, "cme_z/F")
Hadsig.Branch("cme_ndof", cme_ndof, "cme_ndof/I")
Hadsig.Branch("cme_pdgId", cme_pdgId, "cme_pdgId/I")
Hadsig.Branch("cme_jet_btagCMVA", cme_jet_btagCMVA, "cme_jet_btagCMVA/F")
Hadsig.Branch("cme_jet_btagCSVV2", cme_jet_btagCSVV2, "cme_jet_btagCSVV2/F")
Hadsig.Branch("cme_jet_btagDeepB", cme_jet_btagDeepB, "cme_jet_btagDeepB/F")
Hadsig.Branch("cme_jet_btagDeepC", cme_jet_btagDeepC, "cme_jet_btagDeepC/F")
Hadsig.Branch("cme_nDau", cme_nDau, "cme_nDau/I")
Hadsig.Branch("cme_dau1_chi2", cme_dau1_chi2, "cme_dau1_chi2/F")
Hadsig.Branch("cme_dau1_pt", cme_dau1_pt, "cme_dau1_pt/F")
Hadsig.Branch("cme_dau1_ipsigXY", cme_dau1_ipsigXY, "cme_dau1_ipsigXY/F")
Hadsig.Branch("cme_dau1_ipsigZ", cme_dau1_ipsigZ, "cme_dau1_ipsigZ/F")
Hadsig.Branch("cme_dau1_nHits", cme_dau1_nHits, "cme_dau1_nHits/F")
Hadsig.Branch("cme_dau2_chi2", cme_dau2_chi2, "cme_dau2_chi2/F")
Hadsig.Branch("cme_dau2_pt", cme_dau2_pt, "cme_dau2_pt/F")
Hadsig.Branch("cme_dau2_ipsigXY", cme_dau2_ipsigXY, "cme_dau2_ipsigXY/F")
Hadsig.Branch("cme_dau2_ipsigZ", cme_dau2_ipsigZ, "cme_dau2_ipsigZ/F")
Hadsig.Branch("cme_dau2_nHits", cme_dau2_nHits, "cme_dau2_nHits/F")
Hadsig.Branch("cmeTruth_nMatched", cmeTruth_nMatched, "cmeTruth_nMatched/I")
Hadsig.Branch("cmeTruth_nTrueDau", cmeTruth_nTrueDau, "cmeTruth_nTrueDau/I")

Hadbkg.Branch("cme_dca", cme_dca, "cme_dca/F")
Hadbkg.Branch("cme_diffMass", cme_diffMass, "cme_diffMass/F")
Hadbkg.Branch("cme_angleXY", cme_angleXY, "cme_angleXY/F")
Hadbkg.Branch("cme_angleXYZ", cme_angleXYZ, "cme_angleXYZ/F")
Hadbkg.Branch("cme_lxy", cme_lxy, "cme_lxy/F")
Hadbkg.Branch("cme_lxyErr", cme_lxyErr, "cme_lxyErr/F")
Hadbkg.Branch("cme_lxyE", cme_lxyE, "cme_lxyE/F")
Hadbkg.Branch("cme_l3D", cme_l3D, "cme_l3D/F")
Hadbkg.Branch("cme_l3DErr", cme_l3DErr, "cme_l3DErr/F")
Hadbkg.Branch("cme_l3DE", cme_l3DE, "cme_l3DE/F")
Hadbkg.Branch("cme_jetDR", cme_jetDR, "cme_jetDR/F")
Hadbkg.Branch("cme_legDR", cme_legDR, "cme_legDR/F")
Hadbkg.Branch("cme_nJet", cme_nJet, "cme_nJet/I")
Hadbkg.Branch("cme_chi2", cme_chi2, "cme_chi2/F")
Hadbkg.Branch("cme_eta", cme_eta, "cme_eta/F")
Hadbkg.Branch("cme_mass", cme_mass, "cme_mass/F")
Hadbkg.Branch("cme_phi", cme_phi, "cme_phi/F")
Hadbkg.Branch("cme_pt", cme_pt, "cme_pt/F")
Hadbkg.Branch("cme_x", cme_x, "cme_x/F")
Hadbkg.Branch("cme_y", cme_y, "cme_y/F")
Hadbkg.Branch("cme_z", cme_z, "cme_z/F")
Hadbkg.Branch("cme_ndof", cme_ndof, "cme_ndof/I")
Hadbkg.Branch("cme_pdgId", cme_pdgId, "cme_pdgId/I")
Hadbkg.Branch("cme_jet_btagCMVA", cme_jet_btagCMVA, "cme_jet_btagCMVA/F")
Hadbkg.Branch("cme_jet_btagCSVV2", cme_jet_btagCSVV2, "cme_jet_btagCSVV2/F")
Hadbkg.Branch("cme_jet_btagDeepB", cme_jet_btagDeepB, "cme_jet_btagDeepB/F")
Hadbkg.Branch("cme_jet_btagDeepC", cme_jet_btagDeepC, "cme_jet_btagDeepC/F")
Hadbkg.Branch("cme_nDau", cme_nDau, "cme_nDau/I")
Hadbkg.Branch("cme_dau1_chi2", cme_dau1_chi2, "cme_dau1_chi2/F")
Hadbkg.Branch("cme_dau1_pt", cme_dau1_pt, "cme_dau1_pt/F")
Hadbkg.Branch("cme_dau1_ipsigXY", cme_dau1_ipsigXY, "cme_dau1_ipsigXY/F")
Hadbkg.Branch("cme_dau1_ipsigZ", cme_dau1_ipsigZ, "cme_dau1_ipsigZ/F")
Hadbkg.Branch("cme_dau1_nHits", cme_dau1_nHits, "cme_dau1_nHits/F")
Hadbkg.Branch("cme_dau2_chi2", cme_dau2_chi2, "cme_dau2_chi2/F")
Hadbkg.Branch("cme_dau2_pt", cme_dau2_pt, "cme_dau2_pt/F")
Hadbkg.Branch("cme_dau2_ipsigXY", cme_dau2_ipsigXY, "cme_dau2_ipsigXY/F")
Hadbkg.Branch("cme_dau2_ipsigZ", cme_dau2_ipsigZ, "cme_dau2_ipsigZ/F")
Hadbkg.Branch("cme_dau2_nHits", cme_dau2_nHits, "cme_dau2_nHits/F")
Hadbkg.Branch("cmeTruth_nMatched", cmeTruth_nMatched, "cmeTruth_nMatched/I")
Hadbkg.Branch("cmeTruth_nTrueDau", cmeTruth_nTrueDau, "cmeTruth_nTrueDau/I")

D0sig.Branch("cme_dca", cme_dca, "cme_dca/F")
D0sig.Branch("cme_diffMass", cme_diffMass, "cme_diffMass/F")
D0sig.Branch("cme_angleXY", cme_angleXY, "cme_angleXY/F")
D0sig.Branch("cme_angleXYZ", cme_angleXYZ, "cme_angleXYZ/F")
D0sig.Branch("cme_lxy", cme_lxy, "cme_lxy/F")
D0sig.Branch("cme_lxyErr", cme_lxyErr, "cme_lxyErr/F")
D0sig.Branch("cme_lxyE", cme_lxyE, "cme_lxyE/F")
D0sig.Branch("cme_l3D", cme_l3D, "cme_l3D/F")
D0sig.Branch("cme_l3DErr", cme_l3DErr, "cme_l3DErr/F")
D0sig.Branch("cme_l3DE", cme_l3DE, "cme_l3DE/F")
D0sig.Branch("cme_jetDR", cme_jetDR, "cme_jetDR/F")
D0sig.Branch("cme_legDR", cme_legDR, "cme_legDR/F")
D0sig.Branch("cme_nJet", cme_nJet, "cme_nJet/I")
D0sig.Branch("cme_chi2", cme_chi2, "cme_chi2/F")
D0sig.Branch("cme_eta", cme_eta, "cme_eta/F")
D0sig.Branch("cme_mass", cme_mass, "cme_mass/F")
D0sig.Branch("cme_phi", cme_phi, "cme_phi/F")
D0sig.Branch("cme_pt", cme_pt, "cme_pt/F")
D0sig.Branch("cme_x", cme_x, "cme_x/F")
D0sig.Branch("cme_y", cme_y, "cme_y/F")
D0sig.Branch("cme_z", cme_z, "cme_z/F")
D0sig.Branch("cme_ndof", cme_ndof, "cme_ndof/I")
D0sig.Branch("cme_pdgId", cme_pdgId, "cme_pdgId/I")
D0sig.Branch("cme_jet_btagCMVA", cme_jet_btagCMVA, "cme_jet_btagCMVA/F")
D0sig.Branch("cme_jet_btagCSVV2", cme_jet_btagCSVV2, "cme_jet_btagCSVV2/F")
D0sig.Branch("cme_jet_btagDeepB", cme_jet_btagDeepB, "cme_jet_btagDeepB/F")
D0sig.Branch("cme_jet_btagDeepC", cme_jet_btagDeepC, "cme_jet_btagDeepC/F")
D0sig.Branch("cme_nDau", cme_nDau, "cme_nDau/I")
D0sig.Branch("cme_dau1_chi2", cme_dau1_chi2, "cme_dau1_chi2/F")
D0sig.Branch("cme_dau1_pt", cme_dau1_pt, "cme_dau1_pt/F")
D0sig.Branch("cme_dau1_ipsigXY", cme_dau1_ipsigXY, "cme_dau1_ipsigXY/F")
D0sig.Branch("cme_dau1_ipsigZ", cme_dau1_ipsigZ, "cme_dau1_ipsigZ/F")
D0sig.Branch("cme_dau1_nHits", cme_dau1_nHits, "cme_dau1_nHits/F")
D0sig.Branch("cme_dau2_chi2", cme_dau2_chi2, "cme_dau2_chi2/F")
D0sig.Branch("cme_dau2_pt", cme_dau2_pt, "cme_dau2_pt/F")
D0sig.Branch("cme_dau2_ipsigXY", cme_dau2_ipsigXY, "cme_dau2_ipsigXY/F")
D0sig.Branch("cme_dau2_ipsigZ", cme_dau2_ipsigZ, "cme_dau2_ipsigZ/F")
D0sig.Branch("cme_dau2_nHits", cme_dau2_nHits, "cme_dau2_nHits/F")
D0sig.Branch("cmeTruth_nMatched", cmeTruth_nMatched, "cmeTruth_nMatched/I")
D0sig.Branch("cmeTruth_nTrueDau", cmeTruth_nTrueDau, "cmeTruth_nTrueDau/I")

D0bkg.Branch("cme_dca", cme_dca, "cme_dca/F")
D0bkg.Branch("cme_diffMass", cme_diffMass, "cme_diffMass/F")
D0bkg.Branch("cme_angleXY", cme_angleXY, "cme_angleXY/F")
D0bkg.Branch("cme_angleXYZ", cme_angleXYZ, "cme_angleXYZ/F")
D0bkg.Branch("cme_lxy", cme_lxy, "cme_lxy/F")
D0bkg.Branch("cme_lxyErr", cme_lxyErr, "cme_lxyErr/F")
D0bkg.Branch("cme_lxyE", cme_lxyE, "cme_lxyE/F")
D0bkg.Branch("cme_l3D", cme_l3D, "cme_l3D/F")
D0bkg.Branch("cme_l3DErr", cme_l3DErr, "cme_l3DErr/F")
D0bkg.Branch("cme_l3DE", cme_l3DE, "cme_l3DE/F")
D0bkg.Branch("cme_jetDR", cme_jetDR, "cme_jetDR/F")
D0bkg.Branch("cme_legDR", cme_legDR, "cme_legDR/F")
D0bkg.Branch("cme_nJet", cme_nJet, "cme_nJet/I")
D0bkg.Branch("cme_chi2", cme_chi2, "cme_chi2/F")
D0bkg.Branch("cme_eta", cme_eta, "cme_eta/F")
D0bkg.Branch("cme_mass", cme_mass, "cme_mass/F")
D0bkg.Branch("cme_phi", cme_phi, "cme_phi/F")
D0bkg.Branch("cme_pt", cme_pt, "cme_pt/F")
D0bkg.Branch("cme_x", cme_x, "cme_x/F")
D0bkg.Branch("cme_y", cme_y, "cme_y/F")
D0bkg.Branch("cme_z", cme_z, "cme_z/F")
D0bkg.Branch("cme_ndof", cme_ndof, "cme_ndof/I")
D0bkg.Branch("cme_pdgId", cme_pdgId, "cme_pdgId/I")
D0bkg.Branch("cme_jet_btagCMVA", cme_jet_btagCMVA, "cme_jet_btagCMVA/F")
D0bkg.Branch("cme_jet_btagCSVV2", cme_jet_btagCSVV2, "cme_jet_btagCSVV2/F")
D0bkg.Branch("cme_jet_btagDeepB", cme_jet_btagDeepB, "cme_jet_btagDeepB/F")
D0bkg.Branch("cme_jet_btagDeepC", cme_jet_btagDeepC, "cme_jet_btagDeepC/F")
D0bkg.Branch("cme_nDau", cme_nDau, "cme_nDau/I")
D0bkg.Branch("cme_dau1_chi2", cme_dau1_chi2, "cme_dau1_chi2/F")
D0bkg.Branch("cme_dau1_pt", cme_dau1_pt, "cme_dau1_pt/F")
D0bkg.Branch("cme_dau1_ipsigXY", cme_dau1_ipsigXY, "cme_dau1_ipsigXY/F")
D0bkg.Branch("cme_dau1_ipsigZ", cme_dau1_ipsigZ, "cme_dau1_ipsigZ/F")
D0bkg.Branch("cme_dau1_nHits", cme_dau1_nHits, "cme_dau1_nHits/F")
D0bkg.Branch("cme_dau2_chi2", cme_dau2_chi2, "cme_dau2_chi2/F")
D0bkg.Branch("cme_dau2_pt", cme_dau2_pt, "cme_dau2_pt/F")
D0bkg.Branch("cme_dau2_ipsigXY", cme_dau2_ipsigXY, "cme_dau2_ipsigXY/F")
D0bkg.Branch("cme_dau2_ipsigZ", cme_dau2_ipsigZ, "cme_dau2_ipsigZ/F")
D0bkg.Branch("cme_dau2_nHits", cme_dau2_nHits, "cme_dau2_nHits/F")
D0bkg.Branch("cmeTruth_nMatched", cmeTruth_nMatched, "cmeTruth_nMatched/I")
D0bkg.Branch("cmeTruth_nTrueDau", cmeTruth_nTrueDau, "cmeTruth_nTrueDau/I")

# Make histo for njsi
hnD0 = TH1F("nD0", "nD0", 20, -0.5, 19.5)
hnHad = TH1F("nHad", "nHad", 20, -0.5, 19.5)
hnHadAll = TH1F("nHadAll", "nHadAll", 20, -0.5, 19.5)

for i, hadFile in enumerate(FileArg[2:]):
    InFile = TNetXNGFile(hadFile)
    events = InFile.Get("Events")

    for iev, event in enumerate(events):
        nD0 = 0
        nHad = 0
        nHadAll = 0
        for k in range(event.nhad):
            nHadAll+=1 
            #if abs(event.had_pdgId[k]) != 421 and abs(event.had_pdgId[k]) != 443: continue
            cme_dca[0] = event.had_dca[k]
            cme_angleXY[0] = event.had_angleXY[k]
            cme_angleXYZ[0] = event.had_angleXYZ[k]
            cme_lxy[0] = event.had_lxy[k]
            cme_lxyErr[0] = event.had_lxyErr[k]
            cme_lxyE[0] = event.had_lxy[k] / event.had_lxyErr[k]
            cme_l3D[0] = event.had_l3D[k]
            cme_l3DErr[0] = event.had_l3DErr[k]
            cme_l3DE[0] = event.had_l3D[k] / event.had_l3DErr[k]
            cme_jetDR[0] = event.had_jetDR[k]
            cme_legDR[0] = event.had_legDR[k]
            cme_diffMass[0] = event.had_diffMass[k]
            cme_nJet[0] = event.had_nJet[k]
            cmeTruth_nMatched[0] = event.hadTruth_nMatched[k]
            cmeTruth_nTrueDau[0] = event.hadTruth_nTrueDau[k]
            cme_chi2[0] = event.had_chi2[k]
            cme_mass[0] = event.had_mass[k]
            cme_eta[0] = event.had_eta[k]
            cme_phi[0] = event.had_phi[k]
            cme_pt[0] = event.had_pt[k]
            cme_x[0] = event.had_x[k]
            cme_y[0] = event.had_y[k]
            cme_z[0] = event.had_z[k]
            cme_ndof[0] = event.had_ndof[k]
            cme_pdgId[0] = event.had_pdgId[k]
            cme_jet_btagCMVA[0] = event.had_jet_btagCMVA[k]
            cme_jet_btagCSVV2[0] = event.had_jet_btagCSVV2[k]
            cme_jet_btagDeepB[0] = event.had_jet_btagDeepB[k]
            cme_jet_btagDeepC[0] = event.had_jet_btagDeepC[k]
            cme_nDau[0] = event.had_nDau[k]
            cme_dau1_chi2[0] = event.had_dau1_chi2[k]
            cme_dau1_ipsigXY[0] = event.had_dau1_ipsigXY[k]
            cme_dau1_ipsigZ[0] = event.had_dau1_ipsigZ[k]
            cme_dau1_nHits[0] = event.had_dau1_nHits[k]
            cme_dau1_pt[0] = event.had_dau1_pt[k]
            cme_dau2_chi2[0] = event.had_dau2_chi2[k]
            cme_dau2_ipsigXY[0] = event.had_dau2_ipsigXY[k]
            cme_dau2_ipsigZ[0] = event.had_dau2_ipsigZ[k]
            cme_dau2_nHits[0] = event.had_dau2_nHits[k]
            cme_dau2_pt[0] = event.had_dau2_pt[k]
            
            if cme_l3DE[0] >200: continue
            if cme_jetDR[0] >0.3: continue
            if cme_legDR[0] >0.6 : continue
            if cme_dca[0] >1 : continue
            if cme_chi2[0] >10 : continue
            if cme_jet_btagCSVV2[0] < 0.05: continue 

            if cme_dau1_chi2[0] >4 : continue
            if cme_dau1_nHits[0] <3 : continue
            if cme_dau1_pt[0] <0.5 : continue
            if cme_dau2_chi2[0] >3 : continue
            if cme_dau2_nHits[0] <3 : continue
            if cme_dau2_pt[0] <0.5 : continue
            if cme_angleXY[0] <0.95 : continue
            if abs(cme_x[0]) >8 : continue
            if abs(cme_y[0]) >8 : continue
            if abs(cme_z[0]) >20: continue
            
            
            if cmeTruth_nMatched[0] == 2 and cmeTruth_nTrueDau[0] == 2: 
                Hadsig.Fill()
                nHad += 1
                if abs(cme_pdgId[0]) ==421:
                    D0sig.Fill()
                    nD0 += 1
            else: 
                Hadbkg.Fill() 
                nHad += 1
            #cmeTruth_nMatched[0] == 0:
                if abs(cme_pdgId[0]) == 421:
                    D0bkg.Fill()
                    nD0 += 1
        hnD0.Fill(nD0)
        hnHad.Fill(nHad)
        hnHadAll.Fill(nHadAll)


hnD0.Write()
hnHad.Write()
hnHadAll.Write()
f.Write()
f.Close()
