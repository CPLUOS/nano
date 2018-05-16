from ROOT import *
import math, array

f = TFile("/xrootd/store/user/wjjang/test/TMVAVts_10var_dR03_onlyKS_bTag.root")
t = f.Get("dataset/Method_BDTG/BDTG")
h_s = t.Get("MVA_BDTG_S")
h_b = t.Get("MVA_BDTG_B")

w = RooWorkspace()
x = RooRealVar("bdt","bdt",-1,1)
sighist = RooDataHist("sigdatahist","sigdatahist",RooArgList(x), h_s);   
bkghist = RooDataHist("bkgdatahist","bkgdatahist",RooArgList(x), h_b);   
getattr(w, 'import')(x)
getattr(w, 'import')(sighist)
getattr(w, 'import')(bkghist)
w.factory("HistPdf::sig({bdt}, sigdatahist)")
w.factory("HistPdf::bkg({bdt}, bkgdatahist)")
tot = w.factory("SUM:tot(Nsig[100,0,1000000]*sig, Nbkg[10000,0,10000000]*bkg)")

Ntoys = 100
Nexp = 69684853.
val = 0.04

cnv = TCanvas()

hlist_ns=[]
hlist = []
hlist_ns = []
Y = array.array("f",[])
eps = 0.00048; epb = 0.00058;
for vts in xrange(0,14):
    vts = vts/100.
    Y.append(vts)
    h_vts = TH1D("","",100,-0.01,0.14)
    h_ns = TH1D("", "", 100, 0, 1000)
    h_nb = TH1D("", "", 100, 40000, 41000)
    tru_ns = eps*Nexp*(2*vts**2*(1-vts**2))
    tru_nb = (epb*Nexp*(1.0-vts**2)**2)
    for i in range(Ntoys):
        w.var("Nsig").setVal(eps*Nexp*(2*vts**2*(1-vts**2)))
        w.var("Nbkg").setVal(epb*Nexp*(1.0-vts**2)**2)
        toy = tot.generate(RooArgSet(x))
        w.var("Nsig").setVal(1+0.98*eps*Nexp*(2*vts**2*(1-vts**2)))
        w.var("Nbkg").setVal(1+1.02*epb*Nexp*(1.0-vts**2)**2)
        w.var("Nsig").setError(1+0.5*eps*Nexp*(2*vts**2*(1-vts**2)))
        w.var("Nbkg").setError(1+0.5*epb*Nexp*(1.0-vts**2)**2)
        tot.fitTo(toy, RooFit.PrintLevel(-10))
        ns = w.var("Nsig").getVal()
        nb = w.var("Nbkg").getVal()
        ns_tru = ns / eps
        nb_tru = nb / epb
        vts2 = ns_tru/(ns_tru+nb_tru)
        h_vts.Fill(math.sqrt(vts2/(2-vts2)))
        h_ns.Fill(ns)
        h_nb.Fill(nb)
    h_ns.Draw()
    l = TLine(tru_ns, 0, tru_ns, h_ns.GetMaximum())
    l.Draw()
    cnv.Print("h_ns_%f.png" % vts)
    h_nb.Draw()
    l = TLine(tru_nb, 0, tru_nb, h_nb.GetMaximum())
    l.Draw()
    cnv.Print("h_nb_%f.png" % vts)
    hlist_ns.append(h_ns)
    hlist.append(h_vts)

cnv = TCanvas()
hlist[0].SetStats(0)
hlist[0].Draw()
means = array.array("f",[])
stddevs = array.array("f",[])
hlist[0].SetMaximum(max(h.GetMaximum() for h in hlist)*1.05)
for i, h in enumerate(hlist):
    means.append(h.GetMean())
    #stddevs.append(1./math.sqrt(h.GetEntries()-1)*h.GetStdDev()**2)
    print 1./math.sqrt(Ntoys-1)*h.GetStdDev()**2
    print h.GetStdDev()

    stddevs.append(h.GetStdDev())
    h.SetLineColor(i+2)
    h.Draw("same")

cnv.SaveAs("a.png")

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

l = TLine(val,0,val,0.1)
l.SetLineWidth(2)
l.SetLineStyle(2)

leg = TLegend(0.65,0.15,0.8,0.3)
leg.AddEntry(belt_1sig, "#pm 68% CL", "f")
leg.AddEntry(belt_2sig, "#pm 95% CL", "f")
leg.AddEntry(mean, "Mean", "l")
h_ = TH1D("","",100,0,0.1)
h_.SetMaximum(0.1)
h_.SetMinimum(0)
h_.GetXaxis().SetTitle("measured |V_{ts}|")
h_.GetYaxis().SetTitle("generated |V_{ts}|")
h_.SetTitle("Confidence Belt")
h_.SetStats(0)

cnv2 = TCanvas()
h_.Draw()
belt_2sig.Draw("fsame")
belt_1sig.Draw("fsame")
mean.Draw("same")
l.Draw("same")
leg.Draw("same")
cnv2.SaveAs("b.png")

exit(0)
