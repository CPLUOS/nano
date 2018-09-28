from ROOT import *
import array, math

#f = TFile("/cms/ldap_home/wjjang/wj_nanoAOD_CMSSW_9_4_4/src/nano/analysis/test/vts/tmva/vts_dR_04_Jet_With_KS.root")
#t = f.Get("pp_combined_JKS_BDT/Method_BDT/BDT")
f = TFile("/cms/ldap_home/wjjang/wj_nanoAOD_CMSSW_9_4_4/src/nano/analysis/test/vts/tmva/vts_dR_04_Jet.root")
t = f.Get("pp_combined_J_BDT/Method_BDT/BDT")

h_sig = t.Get("MVA_BDT_S")
h_bkg = t.Get("MVA_BDT_B")

w = RooWorkspace()
x = RooRealVar("bdt","bdt", -1, 1)

rooh_sig = RooDataHist("sigHist", "sigHist", RooArgList(x), h_sig)
rooh_bkg = RooDataHist("bkgHist", "bkgHist", RooArgList(x), h_bkg)

getattr(w, 'import')(x)
getattr(w, 'import')(h_sig)
getattr(w, 'import')(h_bkg)
getattr(w, 'import')(rooh_sig)
getattr(w, 'import')(rooh_bkg)

fit_sig = w.factory("HistPdf::sig({bdt}, sigHist )")
fit_sig.fitTo(rooh_sig)

fit_bkg = w.factory("HistPdf::bkg({bdt}, bkgHist)")
fit_bkg.fitTo(rooh_bkg)

tot = w.factory("SUM:tot(Nsig[100,0,100000]*sig, Nbkg[10000,0,10000000]*bkg)")
print(type(fit_sig), type(fit_bkg))

fr = x.frame()
fit_leg = TLegend(0.7,0.7, 0.9, 0.9)
rooh_sig.plotOn(fr, RooFit.MarkerColor(kBlue))
rooh_bkg.plotOn(fr, RooFit.MarkerColor(kRed))
fit_sig.plotOn(fr, RooFit.LineColor(kBlue))
fit_bkg.plotOn(fr, RooFit.LineColor(kRed))
fleg_sig = fit_leg.AddEntry(rooh_sig, "Signal", "lp")
fleg_bkg = fit_leg.AddEntry(rooh_bkg, "Background","lp")
fleg_fit_sig = fit_leg.AddEntry(fit_sig, "Signal Fit", "l")
fleg_fit_bkg = fit_leg.AddEntry(fit_bkg, "Background Fit", "l")
fleg_sig.SetLineColor(kRed)
fleg_sig.SetMarkerStyle(8)
fleg_sig.SetMarkerColor(kRed)
fleg_bkg.SetLineColor(kBlue)
fleg_bkg.SetMarkerStyle(8)
fleg_bkg.SetMarkerColor(kBlue)
fleg_fit_sig.SetLineColor(kRed)
fleg_fit_bkg.SetLineColor(kBlue)

c = TCanvas()
fr.Draw()
fit_leg.Draw()
c.Draw()
c.SaveAs("fit.png")
c.SaveAs("fit.pdf")

Ntoys = 200    # number of toy MC generation
Ntot = 69684853. # total number of ttbar event (cross-section * luminosity)
BR_w_to_lep = 0.111 # Branching ratio of W->lepton + neutrino
BR_tau_to_elec = 0.178174 # Branching ratio of tau -> elec + neutrino
BR_tau_to_mu = 0.173936   # Branching ratio of tau -> mu + neutrino
ep_dilep = (BR_w_to_lep + BR_w_to_lep + BR_w_to_lep*(BR_tau_to_elec + BR_tau_to_mu))**2 # ratio of dilepton channel
ep_top = 0.1221 # ratio of event passing event selection for bbars+bsbar samples
ep_s_sb = 0.846761 # number of sWbW event with reco s jet / total number of sWbW event  + also passing event selection from bbars + bsbar samples
ep_ns_sb = 1 - ep_s_sb # number of sWbW event without reco s jet / total number of sWbW event from bbars + bsbar samples
N_ttbar_reco_dilep = Ntot*ep_dilep*ep_top
nSelJet_avg_bWsW = 2.868 # mean of nSelJet for bbars + bsbar samples
nSelJet_avg_bWbW = 2.84  # mean of nSelJet for bbbar sample

hlist = []
vtslist = []
Y = array.array("f",[])
for vts in xrange(0,14):
    vts = vts/100.
    vts_term_bWsW = 2*vts**2*(1-vts**2) # ratio of ttbar->bWsW
    vts_term_bWbW = (1-vts**2)**2
    Y.append(vts)
    h_vts = TH1D("","",100,-0.01,0.2)
    for i in range(Ntoys):
        w.var("Nsig").setVal(N_ttbar_reco_dilep*(vts_term_bWsW)*ep_s_sb*1) # N sig for bWsW : 1 means number of rec s jet per event (because used samples are bWsW)
        w.var("Nbkg").setVal(N_ttbar_reco_dilep*(vts_term_bWbW)*nSelJet_avg_bWbW  # N bkg for bWbW : because bkg means non-s, nSelJet avg for bWbW is multiflied and ep_bkg_bb is just 1 
                                + N_ttbar_reco_dilep*(vts_term_bWsW)*(ep_s_sb)*(nSelJet_avg_bWsW - 1) # N bkg for bWsW 1 : when rec s jet exist, avg of number of bkg is nSelJet_avg - 1 
				+ N_ttbar_reco_dilep*(vts_term_bWsW)*(ep_ns_sb)*(nSelJet_avg_bWsW) # N bkg for bWsW 2 : when rec s jet does not exist, avg of number of bkg is nSelJet_avg
                            )
        toy = tot.generate(RooArgSet(x))
        res = tot.fitTo(toy, RooFit.Save(), RooFit.PrintLevel(-10))
#	res.Print("r")
#        if res.status() is not 0 : 
#	  continue
        ns = w.var("Nsig").getVal()
        nb = w.var("Nbkg").getVal()
        h_vts.Fill(math.sqrt((ns/(2*ep_s_sb))/(N_ttbar_reco_dilep)))#*(1/math.sqrt(1-vts**2))) 
#        h_vts.Fill(math.sqrt((ns/(2*ep_s_sb))/(N_ttbar_reco_dilep))*(1/math.sqrt(1-vts**2)))
    hlist.append(h_vts)
    vtslist.append(vts)

cnv = TCanvas()
hlist[0].SetStats(0)
hlist[0].Draw()
means = array.array("f",[])
stddevs = array.array("f",[])
hlist[0].SetMaximum(max(h.GetMaximum() for h in hlist)*1.05)
toy_leg = TLegend(0.65,0.65,0.9,0.9)
for i, h in enumerate(hlist):
    print h.GetMean()
    means.append(h.GetMean())
    stddevs.append(h.GetStdDev())
    toy_leg.AddEntry(h, "vts = {0}".format(vtslist[i]))
    h.SetLineColor(i+2)
    h.Draw("same")
toy_leg.Draw()
cnv.SaveAs("toyMC.png")
cnv.SaveAs("toyMC.pdf")

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
cnv2.SaveAs("confidence.pdf")
