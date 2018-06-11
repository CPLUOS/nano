from ROOT import *
import array, math

f = TFile("/xrootd/store/user/wjjang/test/vts_dR_03_KS_lamb.root")
t = f.Get("dataset_sum_235_dR_03_KS_lamb_no_ES.root/Method_BDT/BDT")
#f = TFile("/xrootd/store/user/wjjang/test/tmva_old/TMVAVts_10var_dR03_onlyKS_bTag.root")
#t = f.Get("dataset/Method_BDT/BDT")
h_sig = t.Get("MVA_BDT_S")
h_bkg = t.Get("MVA_BDT_B")

w = RooWorkspace()
x = RooRealVar("bdt","bdt",-1,1)
rooh_sig = RooDataHist("sigHist","sigHist",RooArgList(x), h_sig)
rooh_bkg = RooDataHist("bkgHist","bkgHist",RooArgList(x), h_bkg)

getattr(w, 'import')(x)
getattr(w, 'import')(rooh_sig)
getattr(w, 'import')(rooh_bkg)

fit_sig = w.factory("CBShape::sig(bdt, x0s[.005, -0.05, 0.2], sigmas[0.07, 0, 1], alphas[-2.5, -4, -1], ns[.2, 0, 2])")
#fit_sig = w.factory("Landau::sig(bdt, means[0.05,-0.1,0.15],sigmas[0.1,0.001,1])")
#fit_sig = w.factory("ArgusBG::sig(bdt,m0[0,0.2],cs[-0.5,-1,1])")
fit_bkg1 = w.factory("Gaussian::fit_bkg1(bdt, meanb1[-0.3,-0.4,-0.25],sigmab1[0.05,0.001,1])")
fit_bkg2 = w.factory("Gaussian::fit_bkg2(bdt, meanb2[-0.15,-0.2,0],sigmab2[0.05,0.001,1])")
fit_bkg = w.factory("SUM::bkg(bkg1[3,0,10]*fit_bkg1,bkg2[3,0,10]*fit_bkg2)")
#fit_sig = w.factory("CBShape::sig(bdt, x0s[.1, 0, 1], sigmas[.1, 0, 1], alphas[-2.5, -4, -1], ns[1, 0, 5])")
#fit_bkg = w.factory("CBShape::bkg(bdt, x0b[-0.25, -1, 1], sigmab[.03, 0, 1], alphab[-0.5, -2, 1], nb[10, 5, 15])")
fit_sig.fitTo(rooh_sig)
fit_bkg.fitTo(rooh_bkg)
tot = w.factory("SUM:tot(Nsig[100,0,1000000]*sig, Nbkg[10000,0,10000000]*bkg)")


fr = x.frame()
#tot.plotOn(fr)
rooh_sig.plotOn(fr, RooFit.MarkerColor(kBlue))
rooh_bkg.plotOn(fr, RooFit.MarkerColor(kRed))
fit_sig.plotOn(fr, RooFit.LineColor(kBlue))
fit_bkg.plotOn(fr, RooFit.LineColor(kRed))
c = TCanvas()
fr.Draw()
c.SaveAs("fit.png")

Ntoys = 20
Nexp = 69684853.
eps = 0.00048; epb = 0.00058;

hlist = []
Y = array.array("f",[])
for vts in xrange(0,14):
    vts = vts/100.
    Y.append(vts)
    h_vts = TH1D("","",200,-0.01,0.14)
    for i in range(Ntoys):
        w.var("Nsig").setVal(eps*Nexp*(2*vts**2*(1-vts**2)))
        w.var("Nbkg").setVal(epb*Nexp*(1.0-vts**2)**2)
        toy = tot.generate(RooArgSet(x))
        tot.fitTo(toy, RooFit.PrintLevel(-10))
        ns = w.var("Nsig").getVal()
        nb = w.var("Nbkg").getVal()
        ns_tru = ns / eps
        nb_tru = nb / epb
        vts2 = ns_tru/(ns_tru+nb_tru)
        h_vts.Fill(math.sqrt(vts2/(2-vts2)))
        #h_vts.Fill(math.sqrt(2*vts2/(vts2**2+1)))
    hlist.append(h_vts)

cnv = TCanvas()

hlist[0].SetStats(0)
hlist[0].Draw()
means = array.array("f",[])
stddevs = array.array("f",[])
hlist[0].SetMaximum(max(h.GetMaximum() for h in hlist)*1.05)
for i, h in enumerate(hlist):
    print h.GetMean()
    means.append(h.GetMean())
    stddevs.append(h.GetStdDev())
    h.SetLineColor(i+2)
    h.Draw("same")
cnv.SaveAs("toy_dist.png")

n = len(Y)
belt_1sig = TGraph(2*n)
belt_2sig = TGraph(2*n)
for i in range(n):
    belt_1sig.SetPoint(i,means[i]+stddevs[i],Y[i])
    belt_1sig.SetPoint(n+i,means[n-i-1]-stddevs[n-i-1],Y[n-i-1])
    belt_2sig.SetPoint(i,means[i]+2*stddevs[i],Y[i])
    belt_2sig.SetPoint(n+i,means[n-i-1]-2*stddevs[n-i-1],Y[n-i-1])

belt_1sig.SetFillColor(kRed-7)
belt_2sig.SetFillColor(kRed-10)
mean = TGraph(n, means, Y)
mean.SetLineWidth(2)
mean.SetLineColor(1)

leg = TLegend(0.65,0.15,0.8,0.3)
leg.AddEntry(belt_1sig, "#pm 1#sigma (68%)", "f")
leg.AddEntry(belt_2sig, "#pm 2#sigma (95%)", "f")
leg.AddEntry(mean, "Mean", "l")
h_ = TH1D("","",100,0,0.13)
h_.SetMaximum(0.13)
h_.SetMinimum(0)
h_.GetXaxis().SetTitle("|V_{ts}|")
h_.GetYaxis().SetTitle("|V_{ts}|'")
h_.SetTitle("Confidence Belt")
h_.SetStats(0)

cnv2 = TCanvas()
h_.Draw()
cnv2.SetGrid()
belt_2sig.Draw("fsame")
belt_1sig.Draw("fsame")
mean.Draw("same")
leg.Draw("same")
cnv2.SaveAs("confidence.png")

exit(0)
