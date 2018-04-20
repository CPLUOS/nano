import ROOT, math, os, copy, getopt, sys, array, math, glob
from ROOT import *
from array import array


# True or Fake
def TORF(event):
    Tcmesons_var=[]
    Fcmesons_var=[]
    for i in range(event.nhad):
        if abs(event.had_pdgId[i]) == number: 
            if event.hadTruth_nMatched[i] == 2 && event.hadTruth_nTrueDau[i] == 2: 
                Tcmeson_var = event.had_var[i]
                Tcmesons_var.append(Tcmeson_var)
            elif event.hadTruth_nMatched[i] == 0: 
                Fcmeson_var = event.had_var[i]
                Fcmesons_var.append(Fcmeson_var)
        
    return Tcmesons_var, Fcmesons_var


# Making histograms
tstack_F = ROOT.THStack("TorF_var_whatcm", "TorF_var_whatcm")
tstack_T = ROOT.THStack("TorF_var_whatcm", "TorF_var_whatcm")

#datalumi = 38*1000
#ttbardir = "/xrootd/store/user/jlee/tsW_13TeV_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v4_hadAOD/"
#samples = [ttbardir]
##xsec = [831.76]

#hTlist = []
#hFlist = []
#for j, sampledir in enumerate(samples):
hF = copy.deepcopy(ROOT.TH1D("FAKE", "FAKE", 50, 0, .5))  
hT = copy.deepcopy(ROOT.TH1D("TRUE", "TRUE", 50, 0, .5))

legend = ROOT.TLegend(.85, .85, .95, .95)
legend.AddEntry(hT, "TRUE")
legend.AddEntry(hF, "FAKE")

nevents =0
#filelist = [l for l in os.listdir(sampledir) if "root" in l]
#for i, fileName in enumerate(filelist):
#    print fileName
inFile = ROOT.TFile("")
events = inFile.Get("Events")
nevents += events.GetEntries()

for iev, event in enumerate(events):

    selecT, selecF = TORF(event)

    for l in range(len(selecF)):
        hF.Fill(selecF[l])
    for l in range(len(selecT)):
        hT.Fill(selecT[l])

inFile.Close()

hF.SetTitle("Fake")
hT.SetTitle("True")

#scale = datalumi*xsec[j]/float(nevents)
scale1 = 1/hF.Integral(-100,10000)
scale2 = 1/hT.Integral(-100,10000)
#comp1 = hF.Integral()/hF.Integral(-100,10000)
#comp2 = hT.Integral()/hT.Integral(-100,10000)
#print "Fake: ", comp1
#print "True: ", comp2
hF.Scale(scale1)
hT.Scale(scale2)

outFile = ROOT.TFile("TF_var_whatcm.root", "RECREATE")
outFile.cd()

hF.SetLineColor(kRed)
hF.SetFillColor(kRed)
hF.SetFillStyle(3005)
hF.Write()

hT.SetLineColor(kBlue)
hT.SetFillColor(kBlue)
hT.SetFillStyle(3004)
hT.Write()

Fmax = hF.GetMaximum()
Tmax = hT.GetMaximum()

canv = ROOT.TCanvas()

if Fmax >= Tmax:
    hF.Draw("hist")
    hT.Draw("same")

if Fmax < Tmax:
    hT.Draw("hist")
    hF.Draw("same")

legend.Draw("same")
canv.Print("TF_var_whatcm.png")
    
outFile.Write()
outFile.Close()
