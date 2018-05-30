import ROOT, os, array, math

#filedir = "/cms/ldap_home/tt8888tt/CMSSW_10_0_0_pre2/src/nano/nanoAOD/prod/"
filedir = "/cms/ldap_home/wjjang/wj_nanoAOD_CMSSW_9_4_4/src/hep-tools/hadAOD/"

plotvars = ["mass", "pt",
            "chi2", "dca", "angleXY", "angleXYZ", "lxy/had_lxyErr","l3D/had_l3DErr", "legDR",
            "dau1_ipsigXY","dau1_ipsigZ","dau1_chi2", "dau1_pt",
            "dau2_ipsigXY","dau2_ipsigZ","dau2_chi2", "dau2_pt"]
plotvars = ["had_"+p for p in plotvars]
binnings = ["(1000,0.3,0.7)", "(1000,0,10)",
            "(1000,0,5)","(1000,0,2)","(1000,0,1.01)", "(1000,0,1.01)", "(1000,0,30)","(1000,0,30)","(1000,0,3.5)",
            "(1000,0,50)", "(1000,0,50)", "(1000,0,30)", "(1000,0,10)",
            "(1000,0,50)", "(1000,0,50)", "(1000,0,30)", "(1000,0,10)"]

pid = 310
rocs=[]
#cut = ""
cut = "&&had_dau1_pt>0.95&&had_dau2_pt>0.95"
cnv = ROOT.TCanvas()
for i, plotvar in enumerate(plotvars):
    f = ROOT.TFile(filedir+"nanoAOD_1.root")
    #f = ROOT.TFile(filedir+"nanoAOD.root")
    tree = f.Get("Events")
    tree.Draw(plotvar+">>true"+binnings[i], "hadTruth_nMatched ==2&&had_pdgId==%d"%pid+cut)
    tree.Draw(plotvar+">>fake"+binnings[i], "hadTruth_nMatched !=2&&had_pdgId==%d"%pid+cut)

    h_true = ROOT.gROOT.FindObject("true")
    h_fake = ROOT.gROOT.FindObject("fake")
    nbins = h_true.GetNbinsX()
    tot_true = h_true.Integral(0,nbins+1)
    tot_fake = h_fake.Integral(0,nbins+1)

    #significance
    """
    sigX = array.array("f",[])  
    sigY = array.array("f",[])  
    for j in range(nbins):
        sigX.append(h_true.GetXaxis().GetBinCenter(j+1))
        if ("chi2" in plotvar) or ("dca" in plotvar):
            if math.sqrt(h_true.Integral(0,j+1)+h_fake.Integral(0,j+1)) == 0: sigY.append(0)
            else : sigY.append(h_true.Integral(0,j+1)/math.sqrt(h_true.Integral(0,j+1)+h_fake.Integral(0,j+1)))
        else:
            if math.sqrt(h_true.Integral(j+1,nbins+1)+h_fake.Integral(j+1,nbins+1)) == 0: sigY.append(0)
            else : sigY.append(h_true.Integral(j+1,nbins+1)/math.sqrt(h_true.Integral(j+1,nbins+1)+h_fake.Integral(j+1,nbins+1)))
    significance = ROOT.TGraph(len(sigX), sigX, sigY)
    significance.SetTitle(plotvar)
    """

    h_true.Scale(1./tot_true)
    h_fake.Scale(1./tot_fake)
    print "true: "+str(h_true.Integral(0,nbins+1)), " fake: "+str(h_fake.Integral(0,nbins+1))

    #roc curve
    """
    rocX = array.array("f",[])  
    rocY = array.array("f",[])  
    for j in range(nbins):
        if ("chi2" in plotvar) or ("dca" in plotvar):
            rocX.append(h_true.Integral(0,j+1))
            rocY.append(1-h_fake.Integral(0,j+1))
        else:
            rocX.append(h_true.Integral(j+1,nbins+1))
            rocY.append(1-h_fake.Integral(j+1,nbins+1))
    gr = ROOT.TGraph(len(rocX), rocX, rocY)
    """

    h_true.SetLineWidth(2)
    h_fake.SetLineWidth(2)
    h_true.SetBinContent(1,h_true.GetBinContent(1)+h_true.GetBinContent(0))
    h_fake.SetBinContent(1,h_fake.GetBinContent(1)+h_fake.GetBinContent(0))
    h_true.SetBinContent(nbins,h_true.GetBinContent(nbins)+h_true.GetBinContent(nbins+1))
    h_fake.SetBinContent(nbins,h_fake.GetBinContent(nbins)+h_fake.GetBinContent(nbins+1))
    h_fake.Rebin(10)
    h_true.Rebin(10)
    h_fake.SetTitle(plotvar)
    h_fake.SetStats(0)
    h_fake.SetMaximum(max(h.GetMaximum() for h in [h_fake,h_true])*1.2)
    h_fake.SetMinimum(0)
    h_fake.Draw("hist")
    h_true.Draw("histsame")
    h_true.SetLineColor(2)
    cnv.SaveAs(plotvar+".png")
    if "l3D" in plotvar: cnv.SaveAs("l3DErr.png")
    if "lxy" in plotvar: cnv.SaveAs("lxyErr.png")
    #grs.append(gr)

    #cnv2 = ROOT.TCanvas()
    #significance.Draw()
    #cnv2.SaveAs(plotvar+"_sig.png")
    #if "Sig" in plotvar: cnv2.SaveAs("lxySig_sig.png")

"""
h_init = ROOT.TH1D("", "", 1000, 0, 1)
h_init.SetMinimum(0)
h_init.SetMaximum(1)
h_init.SetStats(0)
leg = ROOT.TLegend(0.2,0.2,0.4,0.4)

cnv2 = ROOT.TCanvas()
h_init.Draw()
for i, gr in enumerate(grs):
    gr.SetLineColor(i+1)
    gr.Draw("same")
    leg.AddEntry(gr, plotvars[i], "l")
leg.Draw()
cnv2.SaveAs("roc.png")
"""
