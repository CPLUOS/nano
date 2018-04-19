import ROOT, math, array
from array import array


#plotvar = ['cme_jetDR']
#plotvar = ['cme_lxy / cme_lxySig', 'cme_l3D', 'cme_angleXY', 'cme_dca', 'cme_chi2', 'cme_jetDR']
plotvar = ['cme_jet_btagCSVV2', 'cme_dau1_ipsigXY', 'cme_dau2_pt', 'cme_dau1_nHits', 'cme_dau2_chi2']

#binnings = ['(5000, 0, 100)', '(5000, 0, 100)', '(5000, -2.5, 2.5)', '(5000, 0, 5)', '(5000, 0, 20)', '(5000, 0, 0.5)'] 
binnings = ['(5000, 0, 2)', '(5000, 0, 20)', '(5000, 0, 50)', '(5000, 0, 50)', '(5000, 0, 10)']
#binnings = ['(5000, 0, 0.5)']

l = len(plotvar)

inFile = ROOT.TFile("/cms/scratch/seulgi/nanoAOD/src/nano/analysis/test/topMass/cut/batch/Results/results_merged/copt_cutbased.root")
d0s = inFile.Get("D0cand")

mg = ROOT.TMultiGraph()
legend = ROOT.TLegend(.65, .15, .85, .35)

for i in xrange(l):
    d0s.Draw(plotvar[i]+">>true"+binnings[i], "cmeTruth_nMatched == 2")
    d0s.Draw(plotvar[i]+">>fake"+binnings[i], "cmeTruth_nMatched == 0")
    h_true = ROOT.gROOT.FindObject("true")
    h_fake = ROOT.gROOT.FindObject("fake")
    t = h_true.Integral()
    f = h_fake.Integral()
    nbins = h_true.GetNbinsX()
     
    sigX = array("f", [])
    sigY = array("f", [])
    sigN = array("f", [])
    for j in range(nbins):
        sigN.append(h_true.GetXaxis().GetBinCenter(j+1))
        if plotvar[i] == 'cme_chi2' or plotvar[i] == 'cme_jetDR' or plotvar[i] == 'cme_dca':
        #if ("cme_dca" or "cme_chi2" or "cme_jetDR" in plotvar[i]):
            s = h_true.Integral(0, j)
            b = h_fake.Integral(0, j)
            #print s
        else:
            s = h_true.Integral(j, nbins)
            b = h_fake.Integral(j, nbins)
        #if s+b == 0: sigY.append(0); continue;
        sigX.append(float(b)/f)
        sigY.append(float(s)/t)

    significance = ROOT.TGraph(len(sigN), sigX, sigY)
    #significance.SetTitle(plotvar[i])
    significance.SetLineColor(2+i)
    legend.AddEntry(significance, plotvar[i],"l")
    mg.Add(significance)


c = ROOT.TCanvas()
c.SetGrid()
mg.Draw("AC")
legend.Draw("same")
mg.SetTitle("ROC curve; False positive rate (1-specificity); True positive(sensitivity)")
mg.GetXaxis().SetLimits(0.,1.)
#mg.SetMinimum(0.)
#mg.SetMaximum(1.)
c.Print("mg.png")

