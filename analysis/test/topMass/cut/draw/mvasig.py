import ROOT, math, array


plotvar = 'tmva_bdtg'
binnings = '(2000, -1, 1)'

inFile = ROOT.TFile("/cms/scratch/seulgi/nanonow/src/nano/analysis/test/topMass/cut/batch/Results/results_merged/copt_tmva.root")
d0s = inFile.Get("D0cand")
print plotvar
d0s.Draw(plotvar+">>true"+binnings, "cmeTruth_nMatched == 2 && cme_mass > 1.81483 && cme_mass < 1.91483")
d0s.Draw(plotvar+">>fake"+binnings, "cmeTruth_nMatched == 0 && cme_mass > 1.81483 && cme_mass < 1.91483")


h_true = ROOT.gROOT.FindObject("true")
h_fake = ROOT.gROOT.FindObject("fake")
nbins = h_true.GetNbinsX()

sigX = array.array("f", [])
sigY = array.array("f", [])
for j in range(nbins):
    sigX.append(h_true.GetXaxis().GetBinCenter(j+1))
    #print "x-axis: ", h_true.GetXaxis().GetBinCenter(j+1)
    #if ("tmva_bdtg" in plotvar):
   #     print "plotvar: ", plotvar
   #     #print "mva: ", tmva_bdtg
   #     
    #s = h_true.Integral(0, j+1)
    #b = h_fake.Integral(0, j+1)
    #print "s: ", s
    #print "b: ", b
    #else:
    s = h_true.Integral(j, nbins)
    b = h_fake.Integral(j, nbins)
    #print "s: ", s
    #print "b: ", b
    #else:
    if s+b == 0: sigY.append(0); continue;
    sigY.append(float(s)/float(math.sqrt(s+b)))
    #print sigY
#print len(sigX)
#print sigY
significance = ROOT.TGraph(len(sigX), sigX, sigY)
n = significance.GetN()
y = significance.GetY()
x = significance.GetX()
locmax = ROOT.TMath.LocMax(n,y)
ymax = y[locmax]
xmax = x[locmax]
print ymax
print xmax
significance.SetTitle("Significance-tmva_bdtg; minimum tmva_bdtg; significance;")
c = ROOT.TCanvas()
significance.Draw()    
c.Print("cc.png")

