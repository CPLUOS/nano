import ROOT, os, getopt, sys, array, math, glob, json
from ROOT import *
from array import array


## Make TTree ##
FileArg = sys.argv
tempdir = FileArg[1]
#Dirname = "/cms/scratch/seulgi/nanoAOD/src/nano/analysis/test/topMass/cut/batch/Results/%s/" %tempdir
Dirname = "/cms/scratch/jdj0715/nanoAOD/src/nano/analysis/test/topMass/cut/tmva/tmva_batch/Results/%s/" %tempdir
if not os.path.isdir(Dirname):
    os.makedirs(Dirname)

temp = FileArg[2].split('/').pop()
cattree = Dirname+temp

f = ROOT.TFile(cattree, "recreate")
D0sig = ROOT.TTree("D0sig", "D0sig")
D0bkg = ROOT.TTree("D0bkg", "D0bkg")
Jpsisig = ROOT.TTree("Jpsisig", "Jpsisig")
Jpsibkg = ROOT.TTree("Jpsibkg", "Jpsibkg")

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

Jpsisig.Branch("cme_dca", cme_dca, "cme_dca/F")
Jpsisig.Branch("cme_diffMass", cme_diffMass, "cme_diffMass/F")
Jpsisig.Branch("cme_angleXY", cme_angleXY, "cme_angleXY/F")
Jpsisig.Branch("cme_angleXYZ", cme_angleXYZ, "cme_angleXYZ/F")
Jpsisig.Branch("cme_lxy", cme_lxy, "cme_lxy/F")
Jpsisig.Branch("cme_lxyErr", cme_lxyErr, "cme_lxyErr/F")
Jpsisig.Branch("cme_lxyE", cme_lxyE, "cme_lxyE/F")
Jpsisig.Branch("cme_l3D", cme_l3D, "cme_l3D/F")
Jpsisig.Branch("cme_l3DErr", cme_l3DErr, "cme_l3DErr/F")
Jpsisig.Branch("cme_l3DE", cme_l3DE, "cme_l3DE/F")
Jpsisig.Branch("cme_jetDR", cme_jetDR, "cme_jetDR/F")
Jpsisig.Branch("cme_legDR", cme_legDR, "cme_legDR/F")
Jpsisig.Branch("cme_nJet", cme_nJet, "cme_nJet/I")
Jpsisig.Branch("cme_chi2", cme_chi2, "cme_chi2/F")
Jpsisig.Branch("cme_eta", cme_eta, "cme_eta/F")
Jpsisig.Branch("cme_mass", cme_mass, "cme_mass/F")
Jpsisig.Branch("cme_phi", cme_phi, "cme_phi/F")
Jpsisig.Branch("cme_pt", cme_pt, "cme_pt/F")
Jpsisig.Branch("cme_x", cme_x, "cme_x/F")
Jpsisig.Branch("cme_y", cme_y, "cme_y/F")
Jpsisig.Branch("cme_z", cme_z, "cme_z/F")
Jpsisig.Branch("cme_ndof", cme_ndof, "cme_ndof/I")
Jpsisig.Branch("cme_pdgId", cme_pdgId, "cme_pdgId/I")
Jpsisig.Branch("cme_jet_btagCMVA", cme_jet_btagCMVA, "cme_jet_btagCMVA/F")
Jpsisig.Branch("cme_jet_btagCSVV2", cme_jet_btagCSVV2, "cme_jet_btagCSVV2/F")
Jpsisig.Branch("cme_jet_btagDeepB", cme_jet_btagDeepB, "cme_jet_btagDeepB/F")
Jpsisig.Branch("cme_jet_btagDeepC", cme_jet_btagDeepC, "cme_jet_btagDeepC/F")
Jpsisig.Branch("cme_nDau", cme_nDau, "cme_nDau/I")
Jpsisig.Branch("cme_dau1_chi2", cme_dau1_chi2, "cme_dau1_chi2/F")
Jpsisig.Branch("cme_dau1_pt", cme_dau1_pt, "cme_dau1_pt/F")
Jpsisig.Branch("cme_dau1_ipsigXY", cme_dau1_ipsigXY, "cme_dau1_ipsigXY/F")
Jpsisig.Branch("cme_dau1_ipsigZ", cme_dau1_ipsigZ, "cme_dau1_ipsigZ/F")
Jpsisig.Branch("cme_dau1_nHits", cme_dau1_nHits, "cme_dau1_nHits/F")
Jpsisig.Branch("cme_dau2_chi2", cme_dau2_chi2, "cme_dau2_chi2/F")
Jpsisig.Branch("cme_dau2_pt", cme_dau2_pt, "cme_dau2_pt/F")
Jpsisig.Branch("cme_dau2_ipsigXY", cme_dau2_ipsigXY, "cme_dau2_ipsigXY/F")
Jpsisig.Branch("cme_dau2_ipsigZ", cme_dau2_ipsigZ, "cme_dau2_ipsigZ/F")
Jpsisig.Branch("cme_dau2_nHits", cme_dau2_nHits, "cme_dau2_nHits/F")
Jpsisig.Branch("cmeTruth_nMatched", cmeTruth_nMatched, "cmeTruth_nMatched/I")
Jpsisig.Branch("cmeTruth_nTrueDau", cmeTruth_nTrueDau, "cmeTruth_nTrueDau/I")

Jpsibkg.Branch("cme_dca", cme_dca, "cme_dca/F")
Jpsibkg.Branch("cme_diffMass", cme_diffMass, "cme_diffMass/F")
Jpsibkg.Branch("cme_angleXY", cme_angleXY, "cme_angleXY/F")
Jpsibkg.Branch("cme_angleXYZ", cme_angleXYZ, "cme_angleXYZ/F")
Jpsibkg.Branch("cme_lxy", cme_lxy, "cme_lxy/F")
Jpsibkg.Branch("cme_lxyErr", cme_lxyErr, "cme_lxyErr/F")
Jpsibkg.Branch("cme_lxyE", cme_lxyE, "cme_lxyE/F")
Jpsibkg.Branch("cme_l3D", cme_l3D, "cme_l3D/F")
Jpsibkg.Branch("cme_l3DErr", cme_l3DErr, "cme_l3DErr/F")
Jpsibkg.Branch("cme_l3DE", cme_l3DE, "cme_l3DE/F")
Jpsibkg.Branch("cme_jetDR", cme_jetDR, "cme_jetDR/F")
Jpsibkg.Branch("cme_legDR", cme_legDR, "cme_legDR/F")
Jpsibkg.Branch("cme_nJet", cme_nJet, "cme_nJet/I")
Jpsibkg.Branch("cme_chi2", cme_chi2, "cme_chi2/F")
Jpsibkg.Branch("cme_eta", cme_eta, "cme_eta/F")
Jpsibkg.Branch("cme_mass", cme_mass, "cme_mass/F")
Jpsibkg.Branch("cme_phi", cme_phi, "cme_phi/F")
Jpsibkg.Branch("cme_pt", cme_pt, "cme_pt/F")
Jpsibkg.Branch("cme_x", cme_x, "cme_x/F")
Jpsibkg.Branch("cme_y", cme_y, "cme_y/F")
Jpsibkg.Branch("cme_z", cme_z, "cme_z/F")
Jpsibkg.Branch("cme_ndof", cme_ndof, "cme_ndof/I")
Jpsibkg.Branch("cme_pdgId", cme_pdgId, "cme_pdgId/I")
Jpsibkg.Branch("cme_jet_btagCMVA", cme_jet_btagCMVA, "cme_jet_btagCMVA/F")
Jpsibkg.Branch("cme_jet_btagCSVV2", cme_jet_btagCSVV2, "cme_jet_btagCSVV2/F")
Jpsibkg.Branch("cme_jet_btagDeepB", cme_jet_btagDeepB, "cme_jet_btagDeepB/F")
Jpsibkg.Branch("cme_jet_btagDeepC", cme_jet_btagDeepC, "cme_jet_btagDeepC/F")
Jpsibkg.Branch("cme_nDau", cme_nDau, "cme_nDau/I")
Jpsibkg.Branch("cme_dau1_chi2", cme_dau1_chi2, "cme_dau1_chi2/F")
Jpsibkg.Branch("cme_dau1_pt", cme_dau1_pt, "cme_dau1_pt/F")
Jpsibkg.Branch("cme_dau1_ipsigXY", cme_dau1_ipsigXY, "cme_dau1_ipsigXY/F")
Jpsibkg.Branch("cme_dau1_ipsigZ", cme_dau1_ipsigZ, "cme_dau1_ipsigZ/F")
Jpsibkg.Branch("cme_dau1_nHits", cme_dau1_nHits, "cme_dau1_nHits/F")
Jpsibkg.Branch("cme_dau2_chi2", cme_dau2_chi2, "cme_dau2_chi2/F")
Jpsibkg.Branch("cme_dau2_pt", cme_dau2_pt, "cme_dau2_pt/F")
Jpsibkg.Branch("cme_dau2_ipsigXY", cme_dau2_ipsigXY, "cme_dau2_ipsigXY/F")
Jpsibkg.Branch("cme_dau2_ipsigZ", cme_dau2_ipsigZ, "cme_dau2_ipsigZ/F")
Jpsibkg.Branch("cme_dau2_nHits", cme_dau2_nHits, "cme_dau2_nHits/F")
Jpsibkg.Branch("cmeTruth_nMatched", cmeTruth_nMatched, "cmeTruth_nMatched/I")
Jpsibkg.Branch("cmeTruth_nTrueDau", cmeTruth_nTrueDau, "cmeTruth_nTrueDau/I")



for i, hadFile in enumerate(FileArg[2:]):
    InFile = TNetXNGFile(hadFile)
    events = InFile.Get("Events")

    for iev, event in enumerate(events):
        for k in range(event.nhad):

            if abs(event.had_pdgId[k]) != 421 and abs(event.had_pdgId[k]) != 443: continue
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
            
            if cmeTruth_nMatched[0] == 2 and cmeTruth_nTrueDau[0] == 2:
                if abs(cme_pdgId[0]) == 421:
                    D0sig.Fill()
                else:
                    Jpsisig.Fill()
            elif cmeTruth_nMatched[0] == 0:
                if abs(cme_pdgId[0]) == 421:
                    D0bkg.Fill()
                else:
                    Jpsibkg.Fill()


f.Write()
f.Close()
