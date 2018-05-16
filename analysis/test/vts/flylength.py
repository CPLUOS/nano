import ROOT, os, array, math, random

filedir = "./"
#filedir = "/cms/ldap_home/tt8888tt/CMSSW_10_0_0_pre2/src/nano/nanoAOD/prod/"

#cut = "had_chi2<1&&had_dca<0.2&&had_lxy>0.25&&had_angleXY>0"
cut = "had_pt>0.5&&had_lxy>0&&had_chi2<5&&had_dca<2&&had_l3D>0.1&&had_angleXY>0.9&&had_lxy/had_lxySig>1"


f = ROOT.TFile(filedir+"nanoAOD.root")
tree = f.Get("Events")
"""
tree.Draw("had_mass>>true(100,0.35,0.65)", "hadTruth_nMatched ==2&&had_pdgId==310&&"+cut)
tree.Draw("had_mass>>fake(100,0.35,0.65)", "hadTruth_nMatched !=2&&had_pdgId==310&&"+cut)

h_true = ROOT.gROOT.FindObject("true")
h_fake = ROOT.gROOT.FindObject("fake")

cnv = ROOT.TCanvas()
h_true.SetLineColor(2)
h_true.SetFillColor(2)
h_fake.SetFillColor(4)
st = ROOT.THStack()
st.Add(h_fake)
st.Add(h_true)
st.Draw()
st.GetHistogram().SetMaximum(max(h.GetMaximum() for h in [h_fake,h_true])*1.2)
st.GetHistogram().SetMinimum(0)
st.Draw()
cnv.SaveAs("mass.png")
"""

KsLifetime = 8.954*(10**(-11)) #[s]
lightSpeed = 3*(10**8) #[m/s]

h_distance = ROOT.TH1D("distance", "distance", 50,0,2)
h_etaVsLength = ROOT.TH2D("Ks Fly length vs eta", "Ks Fly length vs eta", 10,0,2.4, 10,0,2)
h_flyLength_b = ROOT.TH1D("Ks Fly length (barrel)", "Ks Fly length (barrel)", 50,0,2)
h_flyLength_e = ROOT.TH1D("Ks Fly length (endcap)", "Ks Fly length (endcap)", 50,0,2)
outVolume = 0
tot = 0
for evt in tree:
    for i in range(evt.ngenHadron):
        if evt.genHadron_pdgId[i] != 310: continue
        tot +=1
        tlv = ROOT.TLorentzVector()
        tlv.SetPtEtaPhiM(evt.genHadron_pt[i],evt.genHadron_eta[i],evt.genHadron_phi[i],0.497)

        distance = math.sqrt((evt.genHadron_x[i] - evt.genHadron_vx[i])**2+(evt.genHadron_y[i] - evt.genHadron_vy[i])**2+(evt.genHadron_z[i] - evt.genHadron_vz[i])**2) #[cm]
        length = tlv.Beta()*tlv.Gamma()*ROOT.gRandom.Exp(KsLifetime)*lightSpeed
        #print "dis: %5.4f len: %10.4f [m]"%(distance*(10**(-2)), length)

        if abs(tlv.Eta())<1.5 and (abs(evt.genHadron_vx[i])>86.8 or abs(evt.genHadron_vy[i])>86.8 or abs(evt.genHadron_vz[i])>140): outVolume += 1
        if abs(tlv.Eta())>1.5 and (abs(evt.genHadron_vz[i])>220): outVolume += 1

        h_distance.Fill(distance*(10**(-2)))
        h_etaVsLength.Fill(abs(tlv.Eta()),length)
        if abs(tlv.Eta())<1.5 : h_flyLength_b.Fill(length)
        else : h_flyLength_e.Fill(length)
print outVolume/float(tot)
"""
cnv = ROOT.TCanvas()
h_distance.SetBinContent(50, h_distance.GetBinContent(50)+h_distance.GetBinContent(51))
h_distance.Draw()
h_distance.SetStats(0)
cnv.SaveAs("distance.png")

h_etaVsLength.Draw("colztext")
cnv.SaveAs("etaVsLength.png")

#cnv.SetLogy()
h_flyLength_b.SetBinContent(50, h_flyLength_b.GetBinContent(50)+h_flyLength_b.GetBinContent(51))
h_flyLength_e.SetBinContent(50, h_flyLength_e.GetBinContent(50)+h_flyLength_e.GetBinContent(51))
h_flyLength_b.SetLineWidth(2)
h_flyLength_e.SetLineWidth(2)
h_flyLength_b.SetLineColor(8)
h_flyLength_e.SetLineColor(4)
h_flyLength_b.SetMaximum(h_flyLength_b.GetMaximum()*1.2)
h_flyLength_b.Draw("hist")
h_flyLength_e.Draw("samehist")
h_flyLength_b.SetStats(0)
h_flyLength_b.GetYaxis().SetTitle("# of kshort")
h_flyLength_b.GetXaxis().SetTitle("distance [m]")
cnv.SaveAs("flyLength.png")

h_flyLength_b.Scale(1./h_flyLength_b.Integral(0,50))
h_flyLength_e.Scale(1./h_flyLength_e.Integral(0,50))
idxB = h_flyLength_b.GetXaxis().FindBin(0.868)
idxE = h_flyLength_e.GetXaxis().FindBin(2.200)
print "barrel: %.2f%% endcap: %.2f%%"%(h_flyLength_b.Integral(idxB,50)*100, h_flyLength_e.Integral(idxE,50)*100)
"""
