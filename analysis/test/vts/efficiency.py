import ROOT, os, copy, array

filedir = "/cms/ldap_home/wjjang/wj_nanoAOD_CMSSW_9_4_4/src/hep-tools/hadAOD/"

#plotvar = "eta"; binning = "12,-2.4,2.4"
plotvar = "pt"; binning = [40,0,40]
binning = [0,1,2,3,4,5,6,8,10,15,20,30,50]

extracut ="&&genHadron_dau1_pt>0.95&&genHadron_dau2_pt>0.95"
#extracut = ""

def getHist(filedir, plotvar, binning, cut):
    f = ROOT.TFile(filedir)
    if len(binning) == 3:
        h_tmp = ROOT.TH1D("tmp","tmp",binning[0],binning[1],binning[2])
    else:
        h_tmp = ROOT.TH1D("tmp","tmp",len(binning)-1, array.array('f',binning)) 
    f.Events.Project("tmp", "genHadron_%s"%plotvar, cut)
    return copy.deepcopy(h_tmp)

jetId = 3
cut = "abs(genHadron_pdgId)==310&&genHadron_isGenHadFromTop&&abs(genHadron_isGenHadFromTsb)==%d"%jetId

h_tot = getHist(filedir+"nanoAOD_0.root", plotvar, binning, cut)
h_tru = getHist(filedir+"nanoAOD_0.root", plotvar, binning, "genHadron_isMatched&&"+cut+extracut)
for i in range(1, 200):
    h_tot.Add(getHist(filedir+"nanoAOD_%d.root"%i, plotvar, binning, cut))
    h_tru.Add(getHist(filedir+"nanoAOD_%d.root"%i, plotvar, binning, "genHadron_isMatched&&"+cut+extracut))
teff_s = ROOT.TEfficiency(h_tru,h_tot)

jetId = 5
cut = "abs(genHadron_pdgId)==310&&genHadron_isGenHadFromTop&&abs(genHadron_isGenHadFromTsb)==%d"%jetId

h_tot = getHist(filedir+"nanoAOD_0.root", plotvar, binning, cut)
h_tru = getHist(filedir+"nanoAOD_0.root", plotvar, binning, "genHadron_isMatched&&"+cut+extracut)
for i in range(1, 200):
    h_tot.Add(getHist(filedir+"nanoAOD_%d.root"%i, plotvar, binning, cut))
    h_tru.Add(getHist(filedir+"nanoAOD_%d.root"%i, plotvar, binning, "genHadron_isMatched&&"+cut+extracut))
teff_b = ROOT.TEfficiency(h_tru,h_tot)

if len(binning) == 3:
    h_tmp = ROOT.TH1D("tmp","tmp",binning[0],binning[1],binning[2])
else:
    h_tmp = ROOT.TH1D("tmp","tmp",len(binning)-1, array.array('f',binning)) 

c = ROOT.TCanvas()
#c.SetGrid()
teff_s.SetLineColor(2)
teff_b.SetLineColor(4)
h_tmp.SetStats(0)
h_tmp.SetTitle("efficiency")
h_tmp.SetMinimum(0)
h_tmp.SetMaximum(0.22)
h_tmp.GetYaxis().SetTitle("Efficiency")
h_tmp.GetXaxis().SetTitle("genHadron %s"%plotvar)
h_tmp.Draw()
teff_s.Draw("e1same")
teff_b.Draw("e1same")
c.SaveAs("eff_%s.png"%plotvar)

