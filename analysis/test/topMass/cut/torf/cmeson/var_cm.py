import ROOT, math, os, copy, getopt, sys, array, math, glob
from ROOT import *
from array import array


# True or Fake
def TORF(event):
    Tcmesons_var=[]
    Fcmesons_var=[]
    for i in range(event.nhad):
        if event.had_pdgId[i] == 421:
            if event.hadTruth_nMatched[i] == 2 and event.hadTruth_nTrueDau[i] : 
                Tcmeson_var = event.had_var[i]
                Tcmesons_var.append(Tcmeson_var)
            elif event.hadTruth_nMatched[i] == 0: 
                Fcmeson_var = event.had_var[i]
                Fcmesons_var.append(Fcmeson_var)
        
    return Tcmesons_var, Fcmesons_var


    
ttbardir = "/xrootd/store/user/jlee/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/hadAOD/"
samples = [ttbardir]

for j, sampledir in enumerate(samples):
    hF = ROOT.TH1D("FAKE", "FAKE", 50, 0, .5)  
    hT = ROOT.TH1D("TRUE", "TRUE", 50, 0, .5)
    
    legend = ROOT.TLegend(.85, .85, .95, .95)
    legend.AddEntry(hT, "TRUE")
    legend.AddEntry(hF, "FAKE")
    
    nevents =0
    filelist = [l for l in os.listdir(sampledir) if "root" in l]
    for i, fileName in enumerate(filelist):
        print fileName
        #if i == 2: break
        inFile = ROOT.TFile(sampledir+fileName)
        events = inFile.Get("Events")
        nevents += events.GetEntries()
        
        #if i == 2 : break
        
        for iev, event in enumerate(events):
        
            selecT, selecF = TORF(event)
        
            for l in range(len(selecF)):
                hF.Fill(selecF[l])
            for l in range(len(selecT)):
                hT.Fill(selecT[l])
        
        inFile.Close()
    
    hF.SetTitle("TorF_var_d0")
    hT.SetTitle("TorF_var_d0")
    
    scale1 = 1/hF.Integral(-100,10000)
    scale2 = 1/hT.Integral(-100,10000)
    comp1 = hF.Integral()/hF.Integral(-100,10000)
    comp2 = hT.Integral()/hT.Integral(-100,10000)
    print "Fake: ", comp1
    print "True: ", comp2
    hF.Scale(scale1)
    hT.Scale(scale2)

hF.SetLineColor(kRed)
hF.SetFillColor(kRed)
hF.SetFillStyle(3005)

hT.SetLineColor(kBlue)
hT.SetFillColor(kBlue)
hT.SetFillStyle(3004)

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
hF.GetXaxis().SetTitle("var")
hF.GetYaxis().SetTitle("Normalized Entries")
canv.Print("TF_var_d0.png")
    
