import math, array, ROOT, copy, CMS_lumi, tdrstyle
import PhysicsTools.PythonAnalysis.rootplot.core as rootplotcore
tdrstyle.setTDRStyle()
def defTH1(title, name, binning):
    if len(binning) == 3:
        hist = ROOT.TH1D(name, title, binning[0], binning[1], binning[2])
    else:
        hist = ROOT.TH1D(name, title, len(binning)-1, array.array('f', binning))
    return hist

def getTH1(title, binning, tree, plotvar, cut, scale = 0.):
    hist = defTH1(title, "name", binning)
    tree.Project("name", plotvar, cut)
    if hist.GetSumw2N() == 0:
        hist.Sumw2()
    if scale != 0:
        hist.Scale(scale)
    return copy.deepcopy(hist)

def makeTH1(filename, treename, title, binning, plotvar, cut, scale = 0.):
    tfile = ROOT.TFile(filename)
    tree  = tfile.Get(treename)
    
    hist = defTH1(title, "tmp", binning)
    for var in plotvar.split(','):
        hist.Add(getTH1(title, binning, tree, var, cut, scale))
        
    return copy.deepcopy(hist)

def getEntries(filename, treename):
    tfile = ROOT.TFile(filename)
    tree  = tfile.Get(treename)
    return tree.GetEntriesFast()

def getWeightedEntries(filename, treename, plotvar, weight):
    weighthist = makeTH1(filename, treename, '', [1, 0, 1], plotvar, weight)    
    return weighthist.Integral(-1,2)

def divide_canvas(canvas, ratio_fraction):
    margins = [ROOT.gStyle.GetPadTopMargin(), ROOT.gStyle.GetPadBottomMargin()]
    useable_height = 1 - (margins[0] + margins[1])
    canvas.Clear()
    pad = ROOT.TPad('mainPad', 'mainPad', 0., 0., 1., 1.)
    pad.SetFillStyle(4000)
    pad.Draw()
    pad.SetBottomMargin(margins[1] + ratio_fraction * useable_height)
    pad_ratio = ROOT.TPad('ratioPad', 'ratioPad', 0., 0., 1., 1.);
    pad_ratio.SetFillStyle(4000)
    pad_ratio.Draw()
    pad_ratio.SetTopMargin(margins[0] + (1 - ratio_fraction) * useable_height)
    return pad, pad_ratio

def makeCanvas(name, doRatio = False):
    H_ref = 600
    W_ref = 800
    if doRatio:
        H_ref = 800
    canv = ROOT.TCanvas(name,name,W_ref,H_ref)    
    return canv

def setMargins(canvas, doRatio = False):
    H_ref = 600
    W_ref = 800
    if doRatio:
        H_ref = 800
    W = W_ref
    H = H_ref
    T = 0.08*H_ref
    B = 0.12*H_ref 
    L = 0.12*W_ref
    R = 0.04*W_ref
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetFrameFillStyle(0)
    canvas.SetFrameBorderMode(0)
    canvas.SetLeftMargin( L/W )
    canvas.SetRightMargin( R/W )
    #canvas.SetTopMargin( T/H )
    #canvas.SetBottomMargin( B/H )
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    return canvas

def setDefAxis(axis, title, offset):
    axis.SetTitle(title)
    axis.SetTitleOffset(offset)
    axis.SetTitleColor(1)
    axis.SetTitleFont(42)
    axis.SetTitleSize(0.043)
    axis.SetLabelColor(1)
    axis.SetLabelFont(42)
    axis.SetLabelOffset(0.007)
    axis.SetLabelSize(0.03)
    axis.SetAxisColor(1)
    axis.SetTickLength(0.03)
    axis.SetNdivisions(510)
    #axis.SetStripDecimals(True)
    #axis.SetPadTickX(1)
    #axis.SetPadTickY(1)

def setDefTH1Style(th1, x_name, y_name):
    setDefAxis(th1.GetYaxis(),y_name, 1.2)
    setDefAxis(th1.GetXaxis(),x_name, 1)
    ROOT.gStyle.SetStripDecimals(True)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.cd()
    return th1
    
def drawTH1(name, cmsLumi, mclist, data, x_name, y_name, doLog=False, doRatio=True, ratioRange=0.45, legx=0.68, legfontsize=0.030):
    leg = ROOT.TLegend(legx,0.68,legx+0.2,0.91)
    leg.SetBorderSize(0)
    #leg.SetNColumns(2)
    leg.SetTextSize(legfontsize)
    leg.SetTextFont(42)
    leg.SetLineColor(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.AddEntry(data,"Data","lp")
    
    hs = ROOT.THStack("mcstack", "mcstack")    
    hratio = mclist[0].Clone("hratio")
    hratio.Reset()

    leghist = []
    for i, mc in enumerate(mclist):
        hnew = mc.Clone("hnew"+mc.GetName())
        hnew.Sumw2(False)
        hs.Add(hnew)
        hratio.Add(mc)
        inversed = mclist[len(mclist)-1-i]
        if not any(inversed.GetTitle() == s for s in leghist):
            leg.AddEntry(inversed, inversed.GetTitle(), "f")
            leghist.append(inversed.GetTitle())
                        
    #hratio = hs.GetStack().Last()
    hratio.Divide(data,hratio,1.,1.,"B")

    tdrstyle.setTDRStyle()

    setDefTH1Style(data, x_name, y_name)
    data.SetName('data')
    data.SetMaximum(data.GetMaximum()*1.8)
    if doLog:
        data.SetMaximum(data.GetMaximum()*100)
        #data.SetMinimum(10**-3)
    else:
        data.GetYaxis().SetTitleSize(0.04)
        data.GetYaxis().SetLabelSize(0.024)
        data.GetYaxis().SetTitleOffset(1.35)
        
    ratio_fraction = 0
    if doRatio:
        ratio_fraction = 0.3        
        data.GetXaxis().SetLabelSize(0)
        data.GetXaxis().SetTitleSize(0)
        setDefTH1Style(hratio, x_name, "Data/MC")
        hratio.GetYaxis().CenterTitle()
        hratio.GetYaxis().SetNdivisions(5)
            
    canv = makeCanvas(name, doRatio)
    pads=[canv]
    pads = rootplotcore.divide_canvas(canv, ratio_fraction)

    pads[0].cd()
    setMargins(pads[0], doRatio)
    if doLog:
        pads[0].SetLogy()

    data.Draw()
    hs.Draw("same")
    data.Draw("esamex0")
    leg.Draw("same")
    
    pads[0].Update()

    if doRatio:
        pads[1].cd()
        pads[1].SetGridy()
        setMargins(pads[1], doRatio)
        hratio.SetLineColor(1)
        hratio.Draw("e")
        hratio.SetMaximum(1.+ratioRange)
        hratio.SetMinimum(1.-ratioRange)

    for p in pads:
        p.RedrawAxis()
        p.Modified()
        p.Update()

    canv.cd()

    #iPos = 0 # in frame
    iPos = 11 # out frame
    if( iPos==0 ):
        cmsLumi.relPosX = 0.1
    cmsLumi.CMS_lumi(pads[0], 0, iPos)

    canv.Modified()
    canv.Update()
    return copy.deepcopy(canv)

def drellYanEstimation(mc_ee_in, mc_ee_out, mc_mm_in, mc_mm_out,
                       rd_ee_in, rd_mm_in, rd_em_in, kMM, kEE):
    #kMM = math.sqrt(rd_mm_in/rd_ee_in)/2.
    #kEE = math.sqrt(rd_ee_in/rd_mm_in)/2.

    rMC_mm = mc_mm_out/mc_mm_in
    rMC_ee = mc_ee_out/mc_ee_in
    print "rMC_mm  ", rMC_mm
    print "rMC_ee  ", rMC_ee
    
    nOutEst_mm = rMC_mm*(rd_mm_in - rd_em_in*kMM)
    nOutEst_ee = rMC_ee*(rd_ee_in - rd_em_in*kEE)
    return nOutEst_ee/mc_ee_out,nOutEst_mm/mc_mm_out

def findDataSet(name, datasets):
    for data in datasets:
        if data["name"] == name:
            return data
    return None

def adderrs(err1, err2, sign=1.):
    return math.sqrt(err1**2+sign*err2**2)

def table(mchistList, errList, signal_hist, signal_err):
    nums = {}
    errs = {}
    total = total_err = 0

    titles = list(set([mc.GetTitle() for mc in mchistList]))
    for t in titles:
        nums[t] = 0
        errs[t] = 0

    for i, mc in enumerate(mchistList):
        nbins = mc.GetSize()-2
        nums[mc.GetTitle()] += mc.Integral(0,nbins+1)
        errs[mc.GetTitle()] = adderrs(errs[mc.GetTitle()], errList[i])

        total += mc.Integral(0,nbins+1)
        total_err = adderrs(total_err, errList[i])
    
    nums['total'] = total
    errs['total'] = total_err

    bkg = total - signal_hist.Integral(0,signal_hist.GetSize()-1)
    bkg_err = adderrs(total_err, signal_err, -1)
    nums['bkg'] = bkg
    errs['bkg'] = bkg_err

    return nums, errs

def set_palette(name="", ncontours=999):
    """Set a color palette from a given RGB list
    stops, red, green and blue should all be lists of the same length
    see set_decent_colors for an example"""

    if name == "gray" or name == "grayscale":
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [1.00, 0.84, 0.61, 0.34, 0.00]
        green = [1.00, 0.84, 0.61, 0.34, 0.00]
        blue  = [1.00, 0.84, 0.61, 0.34, 0.00]
    # elif name == "whatever":
        # (define more palettes)
    else:
        # default palette, looks cool
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [0.00, 0.00, 0.87, 1.00, 0.51]
        green = [0.00, 0.81, 1.00, 0.20, 0.00]
        blue  = [0.51, 1.00, 0.12, 0.00, 0.00]

    s = array('d', stops)
    r = array('d', red)
    g = array('d', green)
    b = array('d', blue)

    npoints = len(s)
    TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
    gStyle.SetNumberContours(ncontours)

def overFlow(hist):
    nbins = hist.GetNbinsX()
    hist.SetBinContent(nbins, hist.GetBinContent(nbins)+hist.GetBinContent(nbins+1))
    hist.SetBinError(nbins, math.sqrt(hist.GetBinError(nbins)**2+hist.GetBinError(nbins+1)**2))
    hist.SetBinContent(1, hist.GetBinContent(1)+hist.GetBinContent(0))
    hist.SetBinError(1, math.sqrt(hist.GetBinError(1)**2+hist.GetBinError(0)**2))

def extraText(canv, position, content):
    canv.cd()
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextFont(42)
    tex.SetTextSize(0.04)
    tex.DrawLatex(position[0], position[1], content)
    canv.Update()

def sysUncertainty(filename, treename, binning, title, scale, cut, plotvar, h_nom, sysList):
    sysErrs_up = []
    sysErrs_dn = []
    for sys in sysList:
        if 'weight' not in sys:
            sysErrs_up.append(makeTH1(filename, "cattree/%s_u"%sys, title, binning, plotvar, cut, scale))
            sysErrs_dn.append(makeTH1(filename, "cattree/%s_d"%sys, title, binning, plotvar, cut, scale))
        else:
            sysErrs_up.append(makeTH1(filename, treename, title, binning, plotvar, cut.replace(sys,'%s_up'%sys), scale))
            sysErrs_dn.append(makeTH1(filename, treename, title, binning, plotvar, cut.replace(sys,'%s_dn'%sys), scale))

    for i in range(len(sysList)):
        sysErrs_up[i].Add(h_nom, -1)
        sysErrs_dn[i].Add(h_nom, -1)

    return sysErrs_up, sysErrs_dn


def drawRatioPlot(name, cmsLumi, mclist, data, x_name, y_name, doLog=False, doRatio=True, ratioRange=0.45, legx=0.68, legfontsize=0.030):
    leg = ROOT.TLegend(legx,0.68,legx+0.2,0.91)
    leg.SetBorderSize(0)
    #leg.SetNColumns(2)
    leg.SetTextSize(legfontsize)
    leg.SetTextFont(42)
    leg.SetLineColor(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.AddEntry(data,"Data","lp")
    
    hs = ROOT.THStack("mcstack", "mcstack")    
    hAllMC = mclist[0].Clone("hratio")
    hAllMC.Reset()

    leghist = []
    for i, mc in enumerate(mclist):
        hnew = mc.Clone("hnew"+mc.GetName())
        hnew.Sumw2(False)
        hs.Add(hnew)
        hAllMC.Add(mc)
        inversed = mclist[len(mclist)-1-i]
        if not any(inversed.GetTitle() == s for s in leghist):
            leg.AddEntry(inversed, inversed.GetTitle(), "f")
            leghist.append(inversed.GetTitle())
                        
    #hratio = hs.GetStack().Last()
    hratio = ROOT.TRatioPlot(data,hAllMC)
    hratio.SetSeparationMargin(0.0)

    tdrstyle.setTDRStyle()
    
    setDefTH1Style(hratio, x_name, y_name)

    canv = makeCanvas(name, False)
    hratio.Draw("")
    leg.Draw("same")

    cmsLumi.CMS_lumi(pads[0], 0, iPos)

    canv.Modified()
    canv.Update()
    return copy.deepcopy(canv)
