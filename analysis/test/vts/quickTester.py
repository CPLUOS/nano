import ROOT, os, array, math, random

filedir = "/cms/ldap_home/tt8888tt/CMSSW_10_0_0_pre2/src/nano/nanoAOD/prod/"

cut = "had_chi2<3&&had_dca<1&&had_angleXY>0.98&&had_lxy/had_lxyErr>5"
#cut += "&&had_dau1_chi2<10&&had_dau1_pt>0.35&&had_dau1_ipsigZ>-1&&had_dau1_ipsigXY>2"


f = ROOT.TFile(filedir+"nanoAOD_true.root")
tree = f.Get("Events")
tree.Draw("had_mass>>true(100,0.35,0.65)", "hadTruth_nMatched ==2&&had_pdgId==310&&"+cut)
tree.Draw("had_mass>>fake(100,0.35,0.65)", "hadTruth_nMatched !=2&&had_pdgId==310&&"+cut)

h_true = ROOT.gROOT.FindObject("true")
h_fake = ROOT.gROOT.FindObject("fake")
print h_true.Integral(0,101), h_fake.Integral(0,101)

cnv = ROOT.TCanvas()
h_true.SetLineColor(2)
h_true.SetFillColor(2)
h_fake.SetFillColor(4)
st = ROOT.THStack()
h_true.SetStats(0)
h_fake.SetStats(0)
st.Add(h_fake)
st.Add(h_true)
st.Draw()
st.GetHistogram().SetMaximum(max(h.GetMaximum() for h in [h_fake,h_true])*1.2)
st.GetHistogram().SetMinimum(0)
st.Draw()
cnv.SaveAs("mass.png")
