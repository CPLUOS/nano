import ROOT, math


h = ROOT.THStack("h","d0 mass with tmva")

hs = ROOT.TH1D("sig", "sig", 60, 1.7, 2.)
hb = ROOT.TH1D("bkg", "bkg", 60, 1.7, 2.)
#hw = ROOT.TH1D("bkg", "bkg", 60, 1.7, 2.)


inFile = ROOT.TFile("/cms/scratch/seulgi/nanoAOD/src/nano/analysis/test/topMass/cut/batch/Results/results_merged/copt_2tmva.root")
d0s = inFile.Get("D0cand")



numsig = 0
numbkg = 0
for id0, d in enumerate(d0s):
    if d.cmeTruth_nTrueDau != 2: continue
    #if d.cme_mass < 1.86483-0.05 or d.cme_mass > 1.86483+0.05: continue
    if d.tmva_bdtg < 0.8855: continue
    if d.cmeTruth_nMatched == 2:
        numsig = numsig+1
        hs.Fill(d.cme_mass)
    elif d.cmeTruth_nMatched == 0:
        numbkg = numbkg+1
        hb.Fill(d.cme_mass)
    #elif d.cmeTruth_nMatched == 1:
    #    #numbkg = numbkg+1
    #    hw.Fill(d.cme_mass)


signi = float(numsig) /math.sqrt(float(numsig + numbkg))
print "numsig: ", numsig
print "numbkg: ", numbkg
print "signi: ", signi

#scale1 = 1/hs.Integral()
#print scale1
#scale2 = 1/hb.Integral()
#print scale2
##hs.Scale(scale1)
#hb.Scale(scale2)
#hs.Scale(scale1)

c = ROOT.TCanvas()

#hb.Draw("hist")
#hs.Draw("same") 
hb.SetFillColor(4)
hb.SetLineColor(4)
#hw.SetFillColor(5)
#hw.SetLineColor(5)
#hb.SetFillStyle(3004)
#hb.SetFillStyle(3004)
h.Add(hb)
#h.Add(hw)
hs.SetFillColor(2)
hs.SetLineColor(2)
#hs.SetFillStyle(3004)
h.Add(hs)
h.Draw()
#y = hs.GetMaximum() 
#print y
#hb.SetMinimum(0) 
#hb.Draw()
#h.Draw("same") 
c.Print("mass.png")
