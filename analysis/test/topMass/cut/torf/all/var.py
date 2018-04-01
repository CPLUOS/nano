import ROOT, math, os, copy, getopt, sys, array, math, glob
from ROOT import *
from array import array


# True or Fake
def TORF(event):
    Tcmesons_var=[]
    Fcmesons_var=[]
    for i in range(event.ncmeson):
        if event.cmeson_mcMatch[i] > 1: 
            Tcmeson_var = event.cmeson_var[i]
            #Tcmeson2_var = event.cmeson_varSig[i]
            #Tcmeson_var = Tcmeson1_var / Tcmeson2_var
            Tcmesons_var.append(Tcmeson_var)
        elif event.cmeson_mcMatch[i] < 1:
            Fcmeson_var = event.cmeson_var[i]
            #Fcmeson2_var = event.cmeson_varSig[i]
            #Fcmeson_var = Fcmeson1_var / Fcmeson2_var
            Fcmesons_var.append(Fcmeson_var)
    return Tcmesons_var, Fcmesons_var


# Making histograms
tstack_F = ROOT.THStack("TorF_var", "TorF_var")
tstack_T = ROOT.THStack("TorF_var", "TorF_var")

datalumi = 38*1000

#ttbardir1 = "/xrootd/store/user/jlee/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/171127_165639/0002/"
#ttbardir2 = "/xrootd/store/user/jlee/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/171127_165639/0006/"
#ttbardir = "/xrootd/store/user/seulgi/nano_171220/nanoall/"

ttbardir = "/xrootd/store/user/jlee/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180222_144032/0000/"
samples = [ttbardir]
#xsec = [831.76]

hFlist = []
hTlist = []
#leg = []
for j, sampledir in enumerate(samples):
    hF = copy.deepcopy(ROOT.TH1D("FAKE", "FAKE", 50, 0, .5))  
    hT = copy.deepcopy(ROOT.TH1D("TRUE", "TRUE", 50, 0, .5))
    
    legend = ROOT.TLegend(0.85, 0.85, 0.95, 0.95)
    legend.AddEntry(hT, "TRUE")
    legend.AddEntry(hF, "FAKE")

    nevents =0
    filelist = [l for l in os.listdir(sampledir) if "root" in l]
    for i, fileName in enumerate(filelist):
        print fileName
        inFile = ROOT.TFile(sampledir+fileName)
        events = inFile.Get("Events")
        nevents += events.GetEntries()

        if i == 100 : break

        for iev, event in enumerate(events):

            selecT, selecF = TORF(event)
        
            for l in range(len(selecF)):
                hF.Fill(selecF[l])
            for l in range(len(selecT)):
                hT.Fill(selecT[l])
       

        inFile.Close()

    hT.SetTitle(sampledir)
    hF.SetTitle(sampledir)
    scale1 = 1/hF.Integral()
    scale2 = 1/hT.Integral()
    hF.Scale(scale1)
    hT.Scale(scale2)
    hFlist.append(hF)
    hTlist.append(hT)


outFile = ROOT.TFile("TF_var.root", "RECREATE")
outFile.cd()

for i, hF in enumerate(hFlist):
    h.SetLineColor(kRed)
    h.SetFillColor(kRed)
    h.SetFillStyle(3005)
    tstack_F.Add(hF)
    hF.Write()

for i, hT in enumerate(hTlist):
    hT.SetLineColor(kBlue)
    hT.SetFillColor(kBlue)
    hT.SetFillStyle(3004)
    tstack_T.Add(hT)
    hT.Write()

Fmax = tstack_F.GetMaximum()
Tmax = tstack_T.GetMaximum()

tstack_F.Write()
tstack_T.Write()

canv = ROOT.TCanvas()

if Fmax >= Tmax:    
    tstack_F.Draw("hist")
    tstack_T.Draw("same")

if Fmax < Tmax:
    tstack_T.Draw("hist")
    tstack_F.Draw("same")
    
legend.Draw("same")
tstack_F.GetXaxis().SetTitle("var")
tstack_F.GetYaxis().SetTitle("Normalized Entries")
canv.Print("TF_var.png")

outFile.Write()
outFile.Close()
