import ROOT, math


h1 = ROOT.THStack("h1", "d0 mass before applying cut")
h2 = ROOT.THStack("h2", "d0 mass after applying cut")

hbs = ROOT.TH1D("before sig", "before sig", 60, 1.6, 2.2)
hbb = ROOT.TH1D("before bkg", "before bkg", 60, 1.6, 2.2)

hs = ROOT.TH1D("sig", "sig", 60, 1.6, 2.2)
hb = ROOT.TH1D("bkg", "bkg", 60, 1.6, 2.2)

inFile = ROOT.TFile("/cms/scratch/seulgi/nanorecent/src/nano/analysis/test/topMass/cut/batch/Results/results_merged/copt_cutbased.root")
d0s = inFile.Get("D0cand")

bfsig = 0
bfbkg = 0
afsig = 0
afbkg = 0

for id0, d in enumerate(d0s):
    if d.cme_pdgId != 421: continue
    cme_lxyS = d.cme_lxy / d.cme_lxySig
    cme_l3DS = d.cme_l3D / d.cme_l3DSig
    
    if d.cmeTruth_nMatched == 2:
       bfsig = bfsig+1
       hbs.Fill(d.cme_mass)
    if d.cmeTruth_nMatched == 0:
       bfbkg = bfbkg+1
       hbb.Fill(d.cme_mass)


    if cme_lxyS < 1: continue
    if cme_l3DS < 1: continue
    if d.cme_mass < 1.86483-0.05 or d.cme_mass > 1.86483+0.05: continue
    if d.cme_lxy < 0.05: continue
    if d.cme_l3D < 0.05: continue
    if d.cme_angleXY < 0.9: continue
    if d.cme_angleXYZ < 0.9: continue
    if d.cme_jet_btagCMVA < 0: continue
    if d.cme_dau1_ipsigXY < 0.5: continue
    if d.cme_dau1_ipsigZ < 0.5: continue
    if d.cme_dau2_ipsigXY < 0.5: continue
    if d.cme_dau2_ipsigZ < 0.5: continue
    if d.cme_dau2_pt < 1: continue

    if d.cmeTruth_nMatched == 2:
        afsig = afsig+1
        #print d.cme_mass
        hs.Fill(d.cme_mass)
    if d.cmeTruth_nMatched == 0:
        afbkg = afbkg+1
        hb.Fill(d.cme_mass)

print "bfsig: ", bfsig
print "bfbkg: ", bfbkg
print "afsig: ", afsig
print "afbkg: ", afbkg

bfratio = float(bfsig) /math.sqrt(float(bfsig + bfbkg))
afratio = float(afsig) /math.sqrt(float(afsig + afbkg))
print "bfratio: ", bfratio
print "afratio: ", afratio

c = ROOT.TCanvas()
c.SetCanvasSize(800,300)

hbs.SetFillColor(2)
hbb.SetFillColor(4)
hbs.SetLineColor(2)
hbb.SetLineColor(4)
h1.Add(hbb)
h1.Add(hbs)


hs.SetFillColor(2)
hb.SetFillColor(4)
hs.SetLineColor(2)
hb.SetLineColor(4)
h2.Add(hb)
h2.Add(hs)

c.Divide(2,1)
c.cd(1)
h1.Draw()
c.cd(2)
h2.Draw()

c.Print("bf.png")
