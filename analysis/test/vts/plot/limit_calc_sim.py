#!/usr/bin/env python3
# https://twiki.cern.ch/twiki/bin/view/RooStats/RooStatsExercisesMarch2015

import ROOT
from ROOT import RooStats, RooFit, kBlue, kRed, kBlack, kGreen, kYellow, RooWorkspace, TLegend, kDashed
#from ROOT import *
import datetime, os, sys

t_ymd = (str(datetime.datetime.now()).split()[0]).split("-")
now_time = t_ymd[0] + t_ymd[1] + t_ymd[2]

Ntot               = 113628397.1 # total number of ttbar event (cross-section * luminosity)
BR_w_to_e          = 0.1071 # Branching ratio of W->electon + neutrino
BR_w_to_mu         = 0.1063 # Branching ratio of W->muon + neutrino
BR_w_to_tau        = 0.1138 # Branching ratio of W->tau + neutrino
BR_tau_to_e        = 0.1782 # Branching ratio of tau -> elec + neutrino
BR_tau_to_mu       = 0.1739 # Branching ratio of tau -> mu + neutrino
ep_dilep           = (BR_w_to_e + BR_w_to_mu + BR_w_to_tau*(BR_tau_to_e + BR_tau_to_mu))**2 # ratio of dilepton channel

# from /xrootd/store/user/wjjang/test/sum_tt012j/20181204/tt012j_bbars_bsbar_bbbar_sum_630_774_201_combined.root
ep_top             = 0.2032 # ratio of event passing event selection(step4) 
ep_s_sb            = 0.8464 # number of event with reco s jet / ( total number of event passing event selection - number of bWbW event passing event selection) 
ep_ns_sb           = 1 - ep_s_sb # number of sWbW event without reco s jet / ( total number of event passing event selection - number of bWbW event passing event selection)
N_ttbar_reco_dilep = Ntot*ep_dilep*ep_top
#N_ttbar_reco_dilep = 5.32*N_ttbar_reco_dilep # Quick Check for semi-leptonic decay // semi-leptonic decay ratio is about 5.32 times higher than di-leptonic case
nSelJet_avg_bWsW   = 2.855 # mean of nSelJet from /xrootd/store/user/wjjang/test/sum_tt012j/20181204/tt012j_bbars_bsbar_sum_630_774_combined.root
nSelJet_avg_bWbW   = 2.84  # mean of nSelJet from /xrootd/store/user/wjjang/test/sum_tt012j/20181204/tt012j_bbbar_2l_FxFx_sum_201.root
vts                = 0.04133  # http://pdg.lbl.gov/2018/reviews/rpp2018-rev-ckm-matrix.pdf
vtb                = 0.999105 # http://pdg.lbl.gov/2018/reviews/rpp2018-rev-ckm-matrix.pdf
vts_term_bWsW      = 2*(vts**2)*(vtb**2) # + vts**4
vts_term_bWbW      = (vtb**2)**2


f = ROOT.TFile("/cms/ldap_home/wjjang/wj_nanoAOD_CMSSW_9_4_4/src/nano/analysis/test/vts/tmva/output/vts_dR_04_Jet.root")
t  = f.Get(str(sys.argv[1])+"/Method_BDT/BDT")
t2 = f.Get(str(sys.argv[2])+"/Method_BDT/BDT")
 
#t = f.Get("pp_combined_JKS_BDT_highest/Method_BDT/BDT")
#t = f.Get("pp_combined_J_BDT_highest/Method_BDT/BDT")

#t = f.Get("pp_combined_JKS_BDT_highBDT/Method_BDT/BDT")
#t = f.Get("pp_combined_J_BDT_highBDT/Method_BDT/BDT")

#t = f.Get("JKS_wo_KS_highest/Method_BDT/BDT")
#t = f.Get("JKS_wo_KS_highBDT/Method_BDT/BDT")

#t = f.Get("JKS_wo_KS_s_vs_b_highBDT/Method_BDT/BDT")
#t = f.Get("JKS_wo_KS_s_vs_b_highBDT/Method_BDT/BDT")

#t = f.Get("all_wo_KS_highest/Method_BDT/BDT")
#t = f.Get("all_wo_KS_highBDT/Method_BDT/BDT")

#t = f.Get("all_wo_KS_s_vs_b_highBDT/Method_BDT/BDT")
#t = f.Get("all_wo_KS_s_vs_b_highBDT/Method_BDT/BDT")

h_sig = t.Get("MVA_BDT_S")
h_bkg = t.Get("MVA_BDT_B")

h_sig2 = t2.Get("MVA_BDT_S")
h_bkg2 = t2.Get("MVA_BDT_B")

w = RooWorkspace()
c = ROOT.TCanvas("plots","plots",1920, 1080)
c.Divide(6,2)

w.factory("bdt[-1,1]")
rooh_sig  = ROOT.RooDataHist("sigHist", "sigHist", ROOT.RooArgList(w.var("bdt")), h_sig)
rooh_bkg  = ROOT.RooDataHist("bkgHist", "bkgHist", ROOT.RooArgList(w.var("bdt")), h_bkg)

#rooh_sig2 = ROOT.RooDataHist("sigHist2", "sigHist2", ROOT.RooArgList(w.var("bdt")), h_sig2)
#rooh_bkg2 = ROOT.RooDataHist("bkgHist2", "bkgHist2", ROOT.RooArgList(w.var("bdt")), h_bkg2)

#getattr(w, 'import')(h_sig)
#getattr(w, 'import')(h_bkg)
getattr(w, 'import')(rooh_sig)
getattr(w, 'import')(rooh_bkg)

#getattr(w, 'import')(h_sig2)
#getattr(w, 'import')(h_bkg2)
#getattr(w, 'import')(rooh_sig2)
#getattr(w, 'import')(rooh_bkg2)

if t.GetPath().find("_JKS_") is not -1 or t.GetPath().find("JKS_") is not -1:
    w.factory("Gaussian::sig_1(bdt, meanSig1[-0.0317964,-1, 1.], sigmaSig1[0.0654243,0, 1])")
    w.factory("Gaussian::sig_2(bdt, meanSig2[-0.1, -1, 1.], sigmaSig2[0.0656842, 0, 1] )")
    fit_sig = w.factory("SUM::sig(sig_1,f1[.5, 0, 1]*sig_2)")
    fit_sig.fitTo(rooh_sig)

    w.factory("CBShape::bkg_1(bdt,x0[-0.124641, -1, 1.], sigma[0.045158, 0, 1], alpha[-1.45058, -100, 100], n[4.33805,-1000 , 1000])")
    w.factory("Gaussian::bkg_2(bdt, meanBkg1[-0.10002, -1, 1.], sigmaBkg1[0.073841, 0, 1.] )")
    w.factory("Gaussian::bkg_3(bdt, meanBkg2[0.05, -1, 1.], sigmaBkg2[0.03841, 0, 1.] )")
    fit_bkg = w.factory("SUM::bkg(bkg_1, bf1[0.33, 0, 1]*bkg_2, bf2[0.33, 0, 1]*bkg_3)")
    fit_bkg.fitTo(rooh_bkg)
elif t.GetPath().find("_J_") is not -1 or t.GetPath().find("all_") is not -1:
    w.factory("Gaussian::sig_1(bdt, meanSig1[-0.0317964,-1, 1.], sigmaSig1[0.0654243,0, 1])")
    w.factory("Gaussian::sig_2(bdt, meanSig2[-0.1, -1, 1.], sigmaSig2[0.0656842, 0, 1] )")
    w.factory("Gaussian::sig_3(bdt, meanSig3[0, -1, 1], sigmaSig3[0, 0, 1])")
    fit_sig = w.factory("SUM::sig(sig_1, f1[.33, 0, 1]*sig_2, f2[.33, 0, 1]*sig_3)")
    fit_sig.fitTo(rooh_sig)

    w.factory("Gaussian::bkg_1(bdt, meanBkg0[-0.30002, -1, 1.], sigmaBkg0[0.073841, 0, 1.] )")
    w.factory("Gaussian::bkg_2(bdt, meanBkg1[-0.10002, -1, 1.], sigmaBkg1[0.073841, 0, 1.] )")
    w.factory("Gaussian::bkg_3(bdt, meanBkg2[0.05, -1, 1.], sigmaBkg2[0.03841, 0, 1.] )")
    fit_bkg = w.factory("SUM::bkg(bkg_1, bf1[0.33, 0, 1]*bkg_2, bf2[0.33, 0, 1]*bkg_3)")
    fit_bkg.fitTo(rooh_bkg)

w.var("bdt").setRange(-1,1)
#w.factory("bdt2[-1,1]")
rooh_sig2 = ROOT.RooDataHist("sigHist2", "sigHist2", ROOT.RooArgList(w.var("bdt")), h_sig2)
rooh_bkg2 = ROOT.RooDataHist("bkgHist2", "bkgHist2", ROOT.RooArgList(w.var("bdt")), h_bkg2)
getattr(w, 'import')(rooh_sig2)
getattr(w, 'import')(rooh_bkg2)

if t2.GetPath().find("_JKS_") is not -1 or t2.GetPath().find("JKS_") is not -1:
    w.factory("Gaussian::sig2_1(bdt, meanSig2_1[-0.0317964,-1, 1.], sigmaSig2_1[0.0654243,0, 1])")
    w.factory("Gaussian::sig2_2(bdt, meanSig2_2[-0.1, -1, 1.],      sigmaSig2_2[0.0656842, 0, 1] )")
    fit_sig2 = w.factory("SUM::sig2(sig2_1,f2_1[.5, 0, 1]*sig2_2)")
    fit_sig2.fitTo(rooh_sig2)

    w.factory("CBShape::bkg2_1(bdt,x2_0[-0.124641, -1, 1.], sigma2[0.045158, 0, 1], alpha2[-1.45058, -100, 100], n2[4.33805,-1000 , 1000])")
    w.factory("Gaussian::bkg2_2(bdt, meanBkg2_1[-0.10002, -1, 1.], sigmaBkg2_1[0.073841, 0, 1.] )")
    w.factory("Gaussian::bkg2_3(bdt, meanBkg2_2[0.05, -1, 1.], sigmaBkg2_2[0.03841, 0, 1.] )")
    fit_bkg2 = w.factory("SUM::bkg2(bkg2_1, bf2_1[0.33, 0, 1]*bkg2_2, bf2_2[0.33, 0, 1]*bkg2_3)")
    fit_bkg2.fitTo(rooh_bkg2)
elif t2.GetPath().find("_J_") is not -1 or t2.GetPath().find("all_") is not -1:
    w.factory("Gaussian::sig2_1(bdt, meanSig2_1[-0.0317964,-1, 1.], sigmaSig2_1[0.0654243,0, 1])")
    w.factory("Gaussian::sig2_2(bdt, meanSig2_2[-0.1, -1, 1.], sigmaSig2_2[0.0656842, 0, 1] )")
    w.factory("Gaussian::sig2_3(bdt, meanSig2_3[0, -1, 1], sigmaSig2_3[0, 0, 1])")
    fit_sig2 = w.factory("SUM::sig2(sig2_1, f2_1[.33, 0, 1]*sig2_2, f2_2[.33, 0, 1]*sig2_3)")
    fit_sig2.fitTo(rooh_sig2)

    w.factory("Gaussian::bkg2_1(bdt, meanBkg2_0[-0.30002, -1, 1.], sigmaBkg2_0[0.073841, 0, 1.] )")
    w.factory("Gaussian::bkg2_2(bdt, meanBkg2_1[-0.10002, -1, 1.], sigmaBkg2_1[0.073841, 0, 1.] )")
    w.factory("Gaussian::bkg2_3(bdt, meanBkg2_2[0.05, -1, 1.], sigmaBkg2_2[0.03841, 0, 1.] )")
    fit_bkg2 = w.factory("SUM::bkg2(bkg2_1, bf2_1[0.33, 0, 1]*bkg2_2, bf2_2[0.33, 0, 1]*bkg2_3)")
    fit_bkg2.fitTo(rooh_bkg2)


c.cd(1)
bfr = w.var("bdt").frame()
rooh_sig.plotOn(bfr, RooFit.MarkerColor(2))
rooh_bkg.plotOn(bfr, RooFit.MarkerColor(4))
fit_sig.plotOn(bfr, RooFit.LineColor(2))
fit_bkg.plotOn(bfr, RooFit.LineColor(4))
bfr.SetTitle("BDT distribution for JKS")
bfr.Draw()

c.cd(2)
bfr2 = w.var("bdt").frame()
rooh_sig2.plotOn(bfr2, RooFit.MarkerColor(2))
rooh_bkg2.plotOn(bfr2, RooFit.MarkerColor(4))
fit_sig2.plotOn(bfr2, RooFit.LineColor(2))
fit_bkg2.plotOn(bfr2, RooFit.LineColor(4))
bfr2.SetTitle("BDT distribution for Jet")
bfr2.Draw()

for i in range(ROOT.RooArgList(fit_sig.getVariables()).getSize()):
    if ROOT.RooArgList(fit_sig.getVariables()).at(i).GetName().find("bdt") == -1 and ROOT.RooArgList(fit_sig.getVariables()).at(i).GetName().find("nsig") == -1 and ROOT.RooArgList(fit_sig.getVariables()).at(i).GetName().find("nbkg") == -1 :
        ROOT.RooArgList(fit_sig.getVariables()).at(i).setConstant()
for i in range(ROOT.RooArgList(fit_bkg.getVariables()).getSize()):
    if ROOT.RooArgList(fit_bkg.getVariables()).at(i).GetName().find("bdt") == -1 and ROOT.RooArgList(fit_bkg.getVariables()).at(i).GetName().find("nsig") == -1 and ROOT.RooArgList(fit_bkg.getVariables()).at(i).GetName().find("nbkg") == -1 :
        ROOT.RooArgList(fit_bkg.getVariables()).at(i).setConstant()

for i in range(ROOT.RooArgList(fit_sig2.getVariables()).getSize()):
    if ROOT.RooArgList(fit_sig2.getVariables()).at(i).GetName().find("bdt") == -1 and ROOT.RooArgList(fit_sig2.getVariables()).at(i).GetName().find("nsig") == -1 and ROOT.RooArgList(fit_sig2.getVariables()).at(i).GetName().find("nbkg") == -1 :
        ROOT.RooArgList(fit_sig2.getVariables()).at(i).setConstant()
for i in range(ROOT.RooArgList(fit_bkg2.getVariables()).getSize()):
    if ROOT.RooArgList(fit_bkg2.getVariables()).at(i).GetName().find("bdt") == -1 and ROOT.RooArgList(fit_bkg2.getVariables()).at(i).GetName().find("nsig") == -1 and ROOT.RooArgList(fit_bkg2.getVariables()).at(i).GetName().find("nbkg") == -1 :
        ROOT.RooArgList(fit_bkg2.getVariables()).at(i).setConstant()

w.factory("expected_sig[10000]")
w.factory("mu[1, 0 , 2]")
w.var("expected_sig").setVal(N_ttbar_reco_dilep*vts_term_bWsW*ep_s_sb*1)
w.factory("expr::nsig('mu*expected_sig', {expected_sig, mu})")

w.factory("expected_bkg[100000]")
w.factory("nu[1, 0, 2]")
w.var("expected_bkg").setVal(N_ttbar_reco_dilep*(vts_term_bWbW)*nSelJet_avg_bWbW # N bkg for bWbW : because bkg means non-s, nSelJet avg for bWbW is multiflied and ep_bkg_bb is just 1 
                            + N_ttbar_reco_dilep*(vts_term_bWsW)*(ep_s_sb)*(nSelJet_avg_bWsW - 1) # N bkg for bWsW 1 : when rec s jet exist, avg of number of bkg is nSelJet_avg - 1 
                            + N_ttbar_reco_dilep*(vts_term_bWsW)*(ep_ns_sb)*nSelJet_avg_bWsW # N bkg for bWsW 2 : when rec s jet does not exist, avg of number of bkg is nSelJet_avg
                  )
w.factory("expr::nbkg('nu*expected_bkg', {expected_bkg, nu})")
w.factory("SUM::model(nsig*sig, nbkg*bkg)")
w.factory("SUM::model2(nsig*sig2, nbkg*bkg2)")

#w.factory("prod::sig_only(nsig,sig)")
#w.factory("prod::bkg_only(nbkg,bkg)")
#w.factory("expr::sig_only('nsig*sig', {nsig, sig})")
#w.factory("expr::bkg_only('nbkg*bkg', {nbkg, bkg})")

c.cd(3)
data = w.pdf("model").generate(ROOT.RooArgSet(w.var("bdt")))
data2 = w.pdf("model2").generate(ROOT.RooArgSet(w.var("bdt")))

sample = ROOT.RooCategory("sample","sample") ;
sample.defineType("JKS") ;
sample.defineType("Jet") ;
# Construct combined dataset in (bdt,sample)
combData = ROOT.RooDataSet("combData", "combined data", ROOT.RooArgSet(w.var("bdt")), RooFit.Index(sample), RooFit.Import("JKS",data), RooFit.Import("Jet",data2))

bfr3 = w.var("bdt").frame()#(RooFit.Bins(40))
data.plotOn(bfr3, RooFit.MarkerColor(2))
data2.plotOn(bfr3, RooFit.MarkerColor(4))
#w.pdf("model").plotOn(bfr3, RooFit.LineColor(1))
#w.pdf("model2").plotOn(bfr3, RooFit.LineColor(2))
#w.obj("sig_only").plotOn(f, RooFit.LineColor(2))
#w.obj("bkg_only").plotOn(f, RooFit.LineColor(4))
combData.plotOn(bfr3, RooFit.MarkerColor(1))
bfr3.SetTitle("Generated Dataset")
bfr3.Draw()

simPdf = ROOT.RooSimultaneous("simPdf", "simultaneous pdf", sample)
simPdf.addPdf(w.pdf("model"), "JKS")
simPdf.addPdf(w.pdf("model2"),"Jet")
simPdf.fitTo(combData)

getattr(w, 'import')(simPdf)

c.cd(4)
bfr4 = w.var("bdt").frame()
#combData.plotOn(bfr4, RooFit.Cut("sample==sample::JKS")) ;
# Plot "JKS" slice of simultaneous pdf. 
# NBL You _must_ project the sample index category with data using ProjWData 
# as a RooSimultaneous makes no prediction on the shape in the index category 
# and can thus not be integrated
simPdf.plotOn(bfr4, RooFit.Slice(sample,"JKS"), RooFit.ProjWData(ROOT.RooArgSet(sample), combData)) ;
simPdf.plotOn(bfr4, RooFit.Slice(sample,"JKS"), RooFit.Components("sig"), RooFit.ProjWData(ROOT.RooArgSet(sample),combData), RooFit.LineColor(2), RooFit.LineStyle(kDashed))
simPdf.plotOn(bfr4, RooFit.Slice(sample,"JKS"), RooFit.Components("bkg"), RooFit.ProjWData(ROOT.RooArgSet(sample),combData), RooFit.LineColor(4), RooFit.LineStyle(kDashed))
bfr4.SetTitle("JKS slice of sim pdf")
bfr4.Draw()

c.cd(5)
bfr5 = w.var("bdt").frame()
#combData.plotOn(bfr5, RooFit.Cut("sample==sample::Jet")) ;
# Plot "JKS" slice of simultaneous pdf. 
# NBL You _must_ project the sample index category with data using ProjWData 
# as a RooSimultaneous makes no prediction on the shape in the index category 
# and can thus not be integrated
simPdf.plotOn(bfr5, RooFit.Slice(sample,"Jet"), RooFit.ProjWData(ROOT.RooArgSet(sample), combData)) ;
simPdf.plotOn(bfr5, RooFit.Slice(sample,"Jet"), RooFit.Components("sig2"), RooFit.ProjWData(ROOT.RooArgSet(sample),combData), RooFit.LineColor(2), RooFit.LineStyle(kDashed))
simPdf.plotOn(bfr5, RooFit.Slice(sample,"Jet"), RooFit.Components("bkg2"), RooFit.ProjWData(ROOT.RooArgSet(sample),combData), RooFit.LineColor(4), RooFit.LineStyle(kDashed))
bfr5.SetTitle("Jet slice of sim pdf")
bfr5.Draw()

c.cd(6)
bfr6 = w.var("bdt").frame()
#combData.plotOn(bfr5, RooFit.Cut("sample==sample::Jet")) ;
# Plot "JKS" slice of simultaneous pdf. 
# NBL You _must_ project the sample index category with data using ProjWData 
# as a RooSimultaneous makes no prediction on the shape in the index category 
# and can thus not be integrated
simPdf.plotOn(bfr6, RooFit.ProjWData(ROOT.RooArgSet(sample), combData)) ;
simPdf.plotOn(bfr6, RooFit.Components("sig2"), RooFit.ProjWData(ROOT.RooArgSet(sample),combData), RooFit.LineColor(2), RooFit.LineStyle(kDashed))
simPdf.plotOn(bfr6, RooFit.Components("bkg2"), RooFit.ProjWData(ROOT.RooArgSet(sample),combData), RooFit.LineColor(4), RooFit.LineStyle(kDashed))
bfr6.SetTitle("sim pdf")
bfr6.Draw()

#exit()
"""
sct = RooSimWSTool(w)
model_sim = sct.build("model_sim","model",SplitParam("m","c"))
"""
mc = RooStats.ModelConfig("mc", w)
#mc.SetPdf(w.pdf("model"))
mc.SetPdf(simPdf)
mc.SetParametersOfInterest("mu")
mc.SetObservables("bdt")
mc.SetNuisanceParameters("nbkg")
#getattr(w, 'import')(mc)

c.cd(7)
bfr7 = w.var("bdt").frame()
#data = ROOT.RooStats.AsymptoticCalculator.MakeAsimovData(data, mc, ROOT.RooArgSet(w.var("bdt")), ROOT.RooArgSet())
data = ROOT.RooStats.AsymptoticCalculator.MakeAsimovData(combData, mc, ROOT.RooArgSet(w.var("bdt")), ROOT.RooArgSet())
dh = ROOT.RooDataHist("","",w.argSet("bdt"), data)#combData)
err_correction = [dh.set(dh.get(i), dh.weight(dh.get(i)), ROOT.TMath.Sqrt(dh.weight(dh.get(i)))) for i in range(0, dh.numEntries())]
dh.plotOn(bfr7)
bfr7.SetTitle("Asimov Dataset")
bfr7.Draw()

# exit(0)
# h = ROOT.TH1F("hobs", "hobs", 100, -1, 1); data.fillHistogram(h, ROOT.RooArgList(w.var("bdt")))
# h.Draw()
# dh = ROOT.RooDataHist("obs", "obs", ROOT.RooArgList(w.var("bdt")), h)

getattr(w, 'import')(dh)

#w.writeToFile('hgg.root')

# ## PLC
# pl = RooStats.ProfileLikelihoodCalculator(data, mc)
# pl.SetConfidenceLevel(.683)
# interval = pl.GetInterval()
# plot = RooStats.LikelihoodIntervalPlot(interval)
# plot.SetRange(0, 200)
# plot.Draw("")
#exit()
### HypoInv
bMod = mc.Clone()
bMod.SetName("bkg")
pod = w.var("mu").getVal()
w.var("mu").setVal(float(sys.argv[3]))
mu_val = str(sys.argv[3])
bMod.SetSnapshot(w.argSet("mu"))
w.var("mu").setVal(pod)
mc.SetSnapshot(w.argSet("mu"))

# Freq. Calcf
#fcalc = RooStats.FrequentistCalculator(dh, bMod, mc)
#fcalc.SetToys(300,300)

# Asymp. Calc
#fcalc = RooStats.AsymptoticCalculator(dh, bMod, mc)
print("dh", dh)
print("bMod", bMod)
print("mc", mc)
fcalc = RooStats.AsymptoticCalculator(dh, bMod, mc)
# fcalc.SetToys(300,300)

# Hypo Test Inverter
hi = RooStats.HypoTestInverter(fcalc)
hi.SetConfidenceLevel(0.95)
useCLs = False
hi.UseCLs(useCLs)

toymcs = hi.GetHypoTestCalculator().GetTestStatSampler()
# profile likelihood test statistics 
profll = RooStats.ProfileLikelihoodTestStat(mc.GetPdf())
# for CLs (bounded intervals) use one-sided profile likelihood
if (useCLs): profll.SetOneSided(True)

# ratio of profile likelihood - need to pass snapshot for the alt 
# ropl = RooStats.RatioOfProfiledLikelihoodsTestStat(mc.GetPdf(), bMod.GetPdf(), bMod.GetSnapshot())
# set the test statistic to use
toymcs.SetTestStatistic(profll);

c.cd(8)
c.SetLogy()
ROOT.gPad.SetLogy()
#gr.SetRangeUser(0.001, 1)
hi.SetFixedScan(40, 0, 2)
r = hi.GetInterval()
plot = RooStats.HypoTestInverterPlot("HTI_Result_Plot","Asymptotic Interval",r)

obs = plot.MakePlot()
obs.GetYaxis().SetRangeUser(0.00000001, 1.)
obs.Draw()
plot.Draw("same")

sig_list = {2:0.0455, 3:0.0027, 4:0.000063, 5:0.00000057}

for key, value in sig_list.items():
    sig = ROOT.TLine()
    sig.SetLineColor(2)
    sig.DrawLine(0, value, 2.2, value)
    sig_text = ROOT.TLatex()
    sig_text.SetTextColor(2)
    sig_text.DrawLatex(2.3, value, "{0} #sigma".format(key))

try:
    os.mkdir(now_time)
    print(now_time+" folder created")
except OSError:#FileExistsError:
    print(now_time+" folder already exists")
c.Draw()
c.SaveAs("./"+now_time+"/"+"limit_calc_"+str(sys.argv[1])+"_"+str(sys.argv[2])+"_mu_"+mu_val+"_"+now_time+".png")
