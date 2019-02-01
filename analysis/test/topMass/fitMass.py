#!/usr/bin/env python
import types, ROOT,nano.analysis.CMS_lumi, json, os, getopt, sys, copy
from nano.analysis.histoHelper import *
import DYestimation
ROOT.gROOT.SetBatch(True)

ROOT.gSystem.Load("libRooFit")




nFitType = 1


def myGetNameFit(): 
  if nFitType == 0 : # Landau
    return "landau"
  if nFitType == 1 : # Gaussian
    return "gaus"


def myGetFitFunc(name, title, xvar, centralvar, sigmavar):
  if nFitType == 0 : # Landau
    return ROOT.RooLandau(name, title, xvar, centralvar, sigmavar)
  if nFitType == 1 : # Gaussian
    return ROOT.RooGaussian(name, title, xvar, centralvar, sigmavar)
  

def myFitInvMass(hist, strName, x, binning, range, dicStyle):
  w = ROOT.RooWorkspace()
  tf1 = ROOT.TF1("f1_fit", myGetNameFit(), binning[ 1 ], binning[ 2 ])

  fMinRange = binning[ 1 ]
  fMaxRange = binning[ 2 ]
  

  if len(range) == 0 : 
    print("hehe")
    hist.Fit(tf1)
  else : 
    fMinRange = range[ 0 ]
    fMaxRange = range[ 1 ]
    
  
  xfitArg = ROOT.RooArgList(x)
  
  print("xfitarg ok")

  dh = ROOT.RooDataHist("dh", "data histogram", xfitArg, hist)

  #unbinned
  #tree = ROOT.TTree("tree", "tree")
  #dh = ROOT.RooDataSet("dh", "data histogram", tree, ROOT.RooArgSet(x))
  
  print("dh ok")

  dSigma = tf1.GetParameter(2)
  print("dSigma is", dSigma)
  #dSigmaError -> GetParError(2)*10.0 or just GetParError(2)?
  dSigmaError = tf1.GetParError(2)


  dMean = tf1.GetParameter(1)
 
  mean_1 = ROOT.RooRealVar("mean_1", "mean_1", dMean, fMinRange, fMaxRange)
  CB_central = ROOT.RooRealVar("mean" + dicStyle[ "suffix" ], "mean" + dicStyle[ "suffix" ], tf1.GetParameter(1), fMinRange, fMaxRange)
  print("fMinRange is", fMinRange)
  print("fMaxRange is", fMaxRange)
  print("mean is", tf1.GetParameter(1))
  #CB_sigma = ROOT.RooRealVar("sigma" + dicStyle[ "suffix" ], "sigma" + dicStyle[ "suffix" ], dSigma, dSigma - dSigmaError, dSigma + dSigmaError)
  CB_sigma = ROOT.RooRealVar("sigma" + dicStyle[ "suffix" ], "sigma" + dicStyle[ "suffix" ], tf1.GetParameter(2), fMinRange, fMaxRange)
  print("initial ok")
  

  M_t = dicHists[ strKey ][ "mass" ] #top quark mass

  #gauss = myGetFitFunc("gauss", "gauss", x, CB_central, CB_sigma)

  gammaParameter = 2.1 + 0.0026*(M_t) 
  gammaError = gammaParameter * 0.1


  betaParameter = 1.96 + 0.18*(M_t) 
  betaError = betaParameter * 0.1  

  muParameter = 8.0 + 0.01*(M_t) 
  muError = muParameter * 0.1  
  

  #gam  = ROOT.RooRealVar("gam", "gam", 1.3 + 0.0026*(M_t), fMinRange, fMaxRange)
  #gam = ROOT.RooRealVar("gam", "gam", gammaParameter, gammaParameter - gammaError , gammaParameter + gammaError)
  #gam  = ROOT.RooRealVar("gam", "gam", 2.1 + 0.0026*(M_t)) #1
  gam  = ROOT.RooRealVar("gam", "gam", 0.6 + 0.0112*(M_t))#2
  #gam  = ROOT.RooRealVar("gam", "gam", 0.6 + 0.0112*(M_t))#2
  #beta = ROOT.RooRealVar("beta", "beta", 25.5) #23
  #beta = ROOT.RooRealVar("beta", "beta", betaParameter, betaParameter-betaError, betaParameter+betaError)
  #beta = ROOT.RooRealVar("beta", "beta", 1.96 + 0.18*(M_t)) #1 1.96
  #beta = ROOT.RooRealVar("beta", "beta", 21.0 + 0.07*(M_t))#2
  beta = ROOT.RooRealVar("beta", "beta", 21.0 + 0.07*(M_t))#2
  #mu = ROOT.RooRealVar("mu", "mu", 20)
  #mu = ROOT.RooRealVar("mu", "mu", muParameter, muParameter-muError, muParameter+muError)
  #mu = ROOT.RooRealVar("mu", "mu", 8 + 0.01*(M_t)) #1
  #mu = ROOT.RooRealVar("mu", "mu", 17.0 - 0.047*(M_t))#2
  mu = ROOT.RooRealVar("mu", "mu", 17.0 - 0.047*(M_t))#2
  #alpha = 0.19 + 0.0016*(M_t)
  alpha = 0.20 + 0.0016*(M_t)
  minus_alpha = 1 - alpha
  mt = ROOT.RooRealVar("mt", "mt", M_t)
  #alpha_1 = ROOT.RooRealVar("alpha_1", "alpha_1", 0.00385, 0.000001, 0.00401)
  #alpha_1 = ROOT.RooRealVar("alpha_1", "alpha_1", 0.00401385264, 0, 0.005)
  #alpha_1 = ROOT.RooRealVar("alpha_1", "alpha_1", 0.385839, 0.375, 0.40) #0.6923895804
  #alpha_1 = ROOT.RooRealVar("alpha_1", "alpha_1", 0.358336, 0.348802, 0.37019) #0.6923895804
  alpha_1 = ROOT.RooRealVar("alpha_1", "alpha_1", 0.470005, 0.000001, 0.6) #0.6923895804
  #n_1 = ROOT.RooRealVar("n_1", "n_1", 0.015, 0.001, 0.03)
  #n_1 = ROOT.RooRealVar("n_1", "n_1", 48.9431, 29, 70)
  #n_1 = ROOT.RooRealVar("n_1", "n_1", 49.0584, 40, 65) #29 80,   40 70
  n_1 = ROOT.RooRealVar("n_1", "n_1", 12, 0, 20) #29 80,   40 70
  nsig = ROOT.RooRealVar("nsig","number of signal events", alpha) #0.46, 0.57
  nbkg = ROOT.RooRealVar("nbkg","number of background events", minus_alpha)
  #CB_alpha = ROOT.RooFormulaVar("CB_alpha", "CB_alpha", "pow(alpha_1*mt, 2)*(-1)", ROOT.RooArgList(mt, alpha_1))
  CB_alpha = ROOT.RooFormulaVar("CB_alpha", "CB_alpha", "(alpha_1)*(-1)", ROOT.RooArgList(alpha_1))
  #CB_alpha = ROOT.RooRealVar("CB_alpha", "CB_alpha", pow((0.0040135)*(M_t),2)*(-1)) #sigmarange
  #CB_n = ROOT.RooRealVar("CB_n", "CB_n", (M_t)*(0.020289)+0.5)
  #CB_n = ROOT.RooFormulaVar("CB_n", "CB_n", "(mt)*(n_1)", ROOT.RooArgList(mt, n_1)) #sigmarange
  CB_n = ROOT.RooFormulaVar("CB_n", "CB_n", "n_1", ROOT.RooArgList(n_1)) #sigmarange
  #CB_n = ROOT.RooRealVar("CB_n", "CB_n", (M_t)*0.073478+1.0)
  

  gamma = ROOT.RooGamma("gamma", "gamma", x, gam, beta, mu)
  crystalBall = ROOT.RooCBShape("crystalBall", "crystalBall shape", x, CB_central, CB_sigma, CB_alpha, CB_n)
  #model = ROOT.RooAddPdf("model","model", ROOT.RooArgList(gauss,gamma), ROOT.RooArgList(nsig,nbkg))
  model = crystalBall
  #model = ROOT.RooAddPdf("model","model", ROOT.RooArgList(gauss,gamma), ROOT.RooArgList(nsig))
  print("modeling ok")
  
  #fitResult = model.fitTo(dh, ROOT.RooFit.Extended(False), ROOT.RooFit.Save())
  fitResult = model.fitTo(dh, ROOT.RooFit.Extended(False), ROOT.RooFit.Minos(),  ROOT.RooFit.Save())
  
  top_mass_frame = x.frame()
  print("frame setting ok")
  strNameHisto = "roofithisto_" + strName
  dh.plotOn(top_mass_frame, ROOT.RooFit.Name(strNameHisto), ROOT.RooFit.MarkerColor(dicStyle[ "color" ]), ROOT.RooFit.MarkerStyle(dicStyle[ "marker" ]))
  print("dh ~ plot ok")


  #g.plotOn(top_mass_frame)
  strNameCurve = "fittingcurve_" + strName
  model.plotOn(top_mass_frame, ROOT.RooFit.Name(strNameCurve), ROOT.RooFit.LineColor(dicStyle[ "color" ]), ROOT.RooFit.LineStyle(dicStyle[ "line" ]))
  #gauss.paramOn(top_mass_frame)
  #model.plotOn(top_mass_frame, ROOT.RooFit.Components("gauss"), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineStyle(ROOT.kDashed))
  #model.plotOn(top_mass_frame, ROOT.RooFit.Components("gamma"), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.LineStyle(ROOT.kDashed))
  
  print "++++++++++++++++#########+++++++++++++++++"  
  print CB_central.getVal()
  #return {"frame":top_mass_frame, "histo":strNameHisto, "graph":strNameCurve, "peak_val":CB_central.getVal(), "peak_err":CB_central.getError(), "sigma_val":CB_sigma.getVal(), "sigma_err":CB_sigma.getError(), "chi2":top_mass_frame.chiSquare(), "CB_alpha_val":CB_alpha.getVal(), "CB_alpha_err":CB_alpha.getError()}
  return {"frame":top_mass_frame, "histo":strNameHisto, "graph":strNameCurve, "peak_val":CB_central.getVal(), "peak_err":CB_central.getError(), "sigma_val":CB_sigma.getVal(), "sigma_err":CB_sigma.getError(), "chi2":top_mass_frame.chiSquare(), "alpha_val":alpha_1.getVal(), "alpha_err":alpha_1.getError(), "n_val":n_1.getVal(), "n_err":n_1.getError()}

################################################################
## Getting inputs
################################################################
strFitType = "Landau"
strFitType = "Gaussian"
rangeFit = []

try:
    opts, args = getopt.getopt(sys.argv[2:], "hdnot:r:",["type","range"])
except getopt.GetoptError:          
    print 'Usage : ./fitmass.py <json> -t <type> -r <range>'
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print 'Usage : ./fitmass.py <json> -t <type> -r <range>'
        sys.exit()
    elif opt in ("-t", "--type"):
        strFitType = arg
    elif opt in ("-r", "--range"):
        rangeFit = eval(arg)

if strFitType == "Landau" :
  nFitType = 0
elif strFitType == "Gaussian" :
  nFitType = 1
  print "Gaussian go"
else :
  print "Error: wrong fitting type"
  sys.exit(2)

################################################################
## Reading the configuration file
################################################################
try:
  if len(sys.argv) is 1 : 
    print "Usage : %s [rootfilename]"
    sys.exit(2)

  fileDicHists = open(sys.argv[ 1 ])
  dicHists = json.load(fileDicHists)
  fileDicHists.close()
except IOError:
  print "Usage : %s [rootfilename]  # rootfilename must be correct"
  sys.exit(2)

listPartPath = sys.argv[ 1 ].split("/")
strPathPos = ""

for i in range(len(listPartPath) - 1) :
  strPathPos = strPathPos + listPartPath[ i ] + "/"
################################################################
## Initializing the result root file
################################################################
rootHists = ROOT.TFile.Open(strPathPos + dicHists["rootfilename"])
rootFits = ROOT.TFile.Open("fres/fittings_" + dicHists["rootfilename"], "RECREATE")

binning = dicHists["binning"]
#x_name = dicHists["x_name"]
x_name = "M_{lepton + D^{0}} [GeV]"
y_name = dicHists["y_name"]

canvasMain = makeCanvas("canvasMain")
strKeyData = ""
#strKeyMC = ""

dicHistoStyle = {
  "data":    {"label":"Data",   "marker":ROOT.kFullDotLarge, "color":ROOT.kBlack, "line":ROOT.kDashed}, 
  "1665":    {"label":"166.5",  "marker":ROOT.kCircle,   "color":ROOT.kRed,  "line":ROOT.kSolid}, 
  "1695":    {"label":"169.5",  "marker":ROOT.kCircle,   "color":ROOT.kPink,  "line":ROOT.kSolid}, 
  "1715": {"label":"171.5", "marker":ROOT.kCircle,       "color":ROOT.kGreen,  "line":ROOT.kSolid}, 
  "nominal": {"label":"172.5", "marker":ROOT.kCircle,    "color":ROOT.kBlue,  "line":ROOT.kSolid}, 
  "1735": {"label":"173.5", "marker":ROOT.kCircle,       "color":ROOT.kOrange,  "line":ROOT.kSolid}, 
  "1755":    {"label":"175.5",  "marker":ROOT.kCircle,   "color":ROOT.kOrange,  "line":ROOT.kSolid},
  "1785":    {"label":"178.5",  "marker":ROOT.kCircle,   "color":ROOT.kViolet,  "line":ROOT.kSolid}
}

################################################################
##  Getting peaks of TT samples
################################################################

for strKey in dicHists.keys():
  if type(dicHists[ strKey ]) != types.DictType or "type" not in dicHists[ strKey ] : continue
  if dicHists[ strKey ][ "type" ] not in ["TT_onlytt", "TT_withbkg"] : continue
  if dicHists[ strKey ][ "type" ] == "data" : strKeyData = strKey
  if dicHists[ strKey ][ "type" ] == "TT_onlytt" or "TT_withbkg" : strKeyData = strKey
  #print strKeyMC 
  
  # If str is missing, since strKey is in unicode type, it yields an error.
  #binned
  histCurr = ROOT.TH1D(rootHists.Get(str(strKey)))
  
  #unbinned will be added
    

  strStyleCurr = ""

  for strStyle in dicHistoStyle.keys() : 
    if strStyle in strKey : 
      strStyleCurr = strStyle
      break
  
  strSuffix = ""
  if dicHists[ strKey ][ "type" ] in ["TT_onlytt"] : 
    strSuffix = "_%0.1f"%(dicHists[ strKey ][ "mass" ])

  dicHistoStyle[ strStyleCurr ][ "suffix" ] = strSuffix
  dicHists[ strKeyData ][ "label" ] = dicHistoStyle[ strStyleCurr ][ "label" ]
  # -- Fitting by RooFit (no bkg)
 # print ""
 # print "#####################################################"
 # print "### Fitting (%s) by RooFit is beginning"%(strKey)
 # print "#####################################################"
 # 
  nFitTypeOrg = nFitType
  x_input = ROOT.RooRealVar("x_input" + strKeyData, histCurr.GetXaxis().GetTitle(), binning[ 1 ], binning[ 2 ])
  
  
  if "Gen" not in dicHists or dicHists[ "Gen" ] == 0: 
    print "##### Using the given distribution #####"
  else:
    print "##### Using Gaussian distribution for Gen #####"
    nFitType = 1

  nFitType = nFitTypeOrg
  
  dicHists[ strKey ][ "fitres" ] = myFitInvMass(histCurr, strKeyData, x_input, binning, rangeFit, dicHistoStyle[ strStyleCurr ])
 
  print("myFitInvMass done")

  # -- Dumping
  #print "------ Peak (%s) by RooFit : %f, Err : %f, chi : %f"%(strKey, 
  #  dicHists[ strKey ][ "fitres" ][ "peak_val" ], 
  #  dicHists[ strKey ][ "fitres" ][ "peak_err" ], 
  #  dicHists[ strKey ][ "fitres" ][ "chi2" ])
  #print "------      (sigma (%s) : %f, %f)"%(strKey, dicHists[ strKey ][ "fitres" ][ "sigma_val" ], 
  #  dicHists[ strKey ][ "fitres" ][ "sigma_err" ])
  
  # -- Setting the name and shapes of histogram
  frameCurr = dicHists[ strKey ][ "fitres" ][ "frame" ]
  #frameCurr = dicHists[ strKey ][ "fitres" ][ xframe ]
  #dicHists[ strKey ][ "fitres" ][ "x" ] = x
  #frameCurr = x.frame()
  
  #frameCurr = xframe
  print("frameCurr is", frameCurr)
  frameCurr.SetName(strKey)
  frameCurr.SetTitle(histCurr.GetTitle())
  print("strKey is", strKey)  


  setDefAxis(frameCurr.GetXaxis(), x_name, 1.0)
  setDefAxis(frameCurr.GetYaxis(), y_name, 1.2)
  frameCurr.SetMinimum(0)
  
  # -- Saving the histogram and the fitting curve
  rootFits.cd()
  frameCurr.Write()
  
  #top_mass_frame.Draw()
  #canvasMain.cd()
  #frameCurr.Draw()


#  strFilenameFits = "fres/fitting_plots_%s_%s.png"%(dicHists["rootfilename"], strType)
#  canvasMain.SaveAs(strFilenameFits)
################################################################
## Drawing all multiple plots
################################################################

for strType in [ "TT_onlytt", "TT_withbkg" ] : 

  # -- Making a legend
  #leg = ROOT.TLegend(0.76,0.73,0.76+0.20,0.91)
  leg = ROOT.TLegend(0.55,0.79,0.80,0.93)
  #leg.SetBorderSize(0)
  leg.SetBorderSize(1)
  leg.SetNColumns(2)
  leg.SetTextSize(0.05)
  leg.SetTextFont(42)
  leg.SetLineColor(0)
  #leg.SetFillColor(0)
  #leg.SetFillStyle(0)
  
  # -- Finding the maximum
  dMax = 0.0
  
  for strKey in dicHists.keys():
    if type(dicHists[ strKey ]) != types.DictType or "type" not in dicHists[ strKey ] : continue
    if dicHists[ strKey ][ "type" ] != strType : continue
    if dicHists[ strKey ][ "type" ] not in [ "TT_onlytt", "TT_withbkg"] : continue
    
    dMaxCurr = dicHists[ strKey ][ "fitres" ][ "frame" ].GetMaximum()
    #dMaxCurrx = dicHists[ strKey ][ "fitres" ][ "frame" ].GetMaximum()
    
    if dMax < dMaxCurr : dMax = dMaxCurr
    print dMaxCurr
    #print dMaxCurrx
    # -- Drawing the data plot first
  dicHists[ strKeyData ][ "fitres" ][ "frame" ].SetMaximum(dMax)
  print "+++++++++++++++++++++++++++++++++++++++++++++" 
  
  canvasMain.cd()
  dicHists[ strKeyData ][ "fitres" ][ "frame" ].Draw()
  
  #canvasMain.cd()
  #xframe.Draw()
   

  '''
  # -- Adding label into the legend
  strNameCurve = dicHists[ strKeyData ][ "fitres" ][ "graph" ]
  curveFit = dicHists[ strKeyData ][ "fitres" ][ "frame" ].findObject(strNameCurve)
  leg.AddEntry(curveFit, " ", "l")
  #
  strNameHisto = dicHists[ strKeyData ][ "fitres" ][ "histo" ]
  histRooFit = dicHists[ strKeyData ][ "fitres" ][ "frame" ].findObject(strNameHisto)
  leg.AddEntry(histRooFit, dicHists[ strKeyData ][ "label" ], "p")
  #
  '''
  #canvasMain.Clear()
  #canvasMain = makeCanvas("canvasMain")
  #canvasMain.Clone()

  for strKey in dicHists.keys():
    if type(dicHists[ strKey ]) != types.DictType or "type" not in dicHists[ strKey ] : continue
    if dicHists[ strKey ][ "type" ] != strType : continue
    if dicHists[ strKey ][ "type" ] not in ["TT_onlytt", "TT_withbkg"] : continue
    

    frameX = dicHists[ strKey ][ "fitres" ][ "frame" ]
    print("frameX is", frameX)    

    #frameX.SetMarkerStyle(ROOT.kStar)
    
    #line = ROOT.TLine(0, 2000, 250, 2000)
    #line.Draw()
    
    canvasMain.cd()
    frameX.Draw("same") #draw all TT_poweheg_mass
    #frameX.Draw() #draw only first one, [1665~1785]
    
    # -- Adding label into the legend
    strNameCurve = dicHists[ strKey ][ "fitres" ][ "graph" ]
    curveFit = dicHists[ strKey ][ "fitres" ][ "frame" ].findObject(strNameCurve)
    leg.AddEntry(curveFit, " ", "l")
    
    fittingValue = "  M_{top}  =  " + dicHists[ strKey ][ "label" ] + "  GeV"

    strNameHisto = dicHists[ strKey ][ "fitres" ][ "histo" ]
    histRooFit = dicHists[ strKey ][ "fitres" ][ "frame" ].findObject(strNameHisto)
    #leg.AddEntry(histRooFit, dicHists[ strKey ][ "label" ], "p")
    leg.AddEntry(histRooFit, fittingValue, "p")
    
    leg.Draw("same")
  
  strFilenameFits = "fres/fitting_plots_%s_%s.png"%(dicHists["rootfilename"], strType)
  canvasMain.SaveAs(strFilenameFits)
  print strFilenameFits

################################################################
##  Prepare to draw the linear plot
################################################################

dBinMinPeak = 165.0
dBinMaxPeak = 180.0
dSizeBin = 0.1

dDatMinPeak =  1048576.0
dDatMaxPeak = -1048576.0

for strKey in dicHists.keys() :
  if type(dicHists[ strKey ]) != types.DictType or "type" not in dicHists[ strKey ] : continue
  if dicHists[ strKey ][ "type" ] not in ["TT_onlytt", "TT_withbkg"] : continue
  
  dDatVal = dicHists[ strKey ][ "fitres" ][ "peak_val" ]
  dDatErr = dicHists[ strKey ][ "fitres" ][ "peak_err" ]
  dAlphaVal = dicHists[ strKey ][ "fitres" ][ "alpha_val" ]
  dAlphaErr = dicHists[ strKey ][ "fitres" ][ "alpha_err" ]
  dNVal = dicHists[ strKey ][ "fitres" ][ "n_val" ]
  dNErr = dicHists[ strKey ][ "fitres" ][ "n_err" ]

  
  if dDatMinPeak > dDatVal - dDatErr : dDatMinPeak = dDatVal - dDatErr
  if dDatMaxPeak < dDatVal + dDatErr : dDatMaxPeak = dDatVal + dDatErr

dDatMean = ( dDatMaxPeak + dDatMinPeak ) / 2
dDatMinPeak = dDatMinPeak - ( dDatMean - dDatMinPeak ) * 1.0
dDatMaxPeak = dDatMaxPeak + ( dDatMaxPeak - dDatMean ) * 1.0
#print "Range in linear plot : ", dDatMinPeak, dDatMean, dDatMaxPeak

dTan = 0.0
dYP  = 0.0

################################################################
##  Plotting the calibration curve by RooFit
################################################################

histPeak_nb = ROOT.TH1D("histPeak_onlyTT", "Peaks (no bkg)", 
  long(( dBinMaxPeak - dBinMinPeak ) / dSizeBin), dBinMinPeak, dBinMaxPeak)

histPeak_bk = ROOT.TH1D("histPeak_withbkg", "Peaks", 
  long(( dBinMaxPeak - dBinMinPeak ) / dSizeBin), dBinMinPeak, dBinMaxPeak)


listDicPeak = [
  #{"name":"no_bkg", "hist": histPeak_nb, "type":"TT_onlytt", 
  #  "label":"Without background", "marker":ROOT.kMultiply, "color":ROOT.kBlue}, 

{"name":"bkg", "hist": histPeak_bk, "type":"TT_withbkg", 
    "label":"With background", "marker":ROOT.kCircle, "color":ROOT.kRed}
]

for dicPeak in listDicPeak : 
 # print "##### Fitting by RooFit (%s)"%(dicPeak[ "name" ])
  histPeak = dicPeak[ "hist" ]
  
  histPeak.SetLineColor(1)

  for strKey in dicHists.keys() :
    if type(dicHists[ strKey ]) != types.DictType or "type" not in dicHists[ strKey ] : continue
    if dicHists[ strKey ][ "type" ] != dicPeak[ "type" ] : continue
    
    dMass = dicHists[ strKey ][ "mass" ]
    dDatVal = dicHists[ strKey ][ "fitres" ][ "peak_val" ]
    dDatErr = dicHists[ strKey ][ "fitres" ][ "peak_err" ]
    dAlphaVal = dicHists[ strKey ][ "fitres" ][ "alpha_val" ]
    dAlphaErr = dicHists[ strKey ][ "fitres" ][ "alpha_err" ]
    dNVal = dicHists[ strKey ][ "fitres" ][ "n_val" ]
    dNErr = dicHists[ strKey ][ "fitres" ][ "n_err" ]
    
    print dMass, "(" + strKey + ")", dDatVal, dDatErr
    print dAlphaVal, dAlphaErr, dNVal, dNErr

    nIdxBin = int(( dMass - dBinMinPeak ) / dSizeBin)
    histPeak.SetBinContent(nIdxBin, dDatVal)
    histPeak.SetBinError(nIdxBin, dDatErr)

  histPeak.SetMinimum(dDatMinPeak)
  histPeak.SetMaximum(dDatMaxPeak)

  tf1 = ROOT.TF1("f1_peaks_RooFit_" + dicPeak[ "name" ], "pol1", dBinMinPeak, dBinMaxPeak)
  histPeak.Fit(tf1)
  
  if dicPeak[ "name" ] == "bkg" : 
    dTan = tf1.GetParameter("p1")
    dYP  = tf1.GetParameter("p0")

  setDefAxis(histPeak.GetXaxis(), "M_{top} [GeV]", 1)
  setDefAxis(histPeak.GetYaxis(), "M_{lepton+D^{0}} [GeV]", 0.9)

  rootFits.cd()
  histPeak.Draw()
  histPeak.Write()

#top quark extract code
##
if dTan != 0 : 
    print "Top mass = %lf +/- %lf [GeV]"%(( dicHists[ strKeyData ][ "fitres" ][ "peak_val" ] - dYP ) / dTan, 
      dicHists[ strKeyData ][ "fitres" ][ "peak_err" ] / dTan)
else : 
    print "Top mass = (Unknown because of zeto-tan)"
##

################################################################
##  Drawing all calibration curves together
################################################################
canvasMain.cd()

# -- Making a legend
#leg = ROOT.TLegend(0.64,0.79,0.64+0.30,0.93)
leg = ROOT.TLegend(0.4,0.79,0.64+0.30,0.93)
leg.SetBorderSize(0)
#leg.SetNColumns(0)
leg.SetTextSize(0.040)
leg.SetTextFont(42)
leg.SetLineColor(0)
leg.SetFillColor(0)
leg.SetFillStyle(0)

strSign = ""

for dicPeak in listDicPeak : 
  histPeak = dicPeak[ "hist" ]
  
  histPeak.SetMarkerStyle(dicPeak[ "marker" ])
  histPeak.SetMarkerColor(dicPeak[ "color" ])
  try : 
    histPeak.GetFunction("f1_peaks_RooFit_" + dicPeak[ "name" ]).SetLineColor(dicPeak[ "color" ])
  except : 
    print dicPeak, histPeak.GetFunction("f1_peaks_RooFit_" + dicPeak[ "name" ])
  
  #ROOT.gStyle.SetOptFit(0)
  histPeak.Draw(strSign + "E1")
  strSign = "same"


  #leg.AddEntry(histPeak, dicPeak[ "label" ], "lp")
  leg.AddEntry(histPeak, "M_{lepton+D^{0}} = %.2f M_{top} + %.2f [GeV]"%(dTan, dYP), "lp")
  leg.Draw()
  #leg.Draw("same")
  #leg = ROOT.TLegend(0.76,0.73,0.76+0.20,0.91)

leg.Draw("same")

#real data band making code

###
polyData = ROOT.TH1D("histPeak_data", "Data", 
  int(( dBinMaxPeak - dBinMinPeak ) / dSizeBin), dBinMinPeak, dBinMaxPeak)

for i in range(int(( dBinMaxPeak - dBinMinPeak ) / dSizeBin) + 1) : 
  polyData.SetBinContent(i, dicHists[ strKeyData ][ "fitres" ][ "peak_val" ])
  polyData.SetBinError(i, dicHists[ strKeyData ][ "fitres" ][ "peak_err" ])

polyData.SetFillColorAlpha(ROOT.kBlack, 0.3)
polyData.SetMarkerStyle(ROOT.kDot)

polyData.Draw("sameE3")

leg.AddEntry(polyData, "Data", "f")

leg.Draw("same")
###


#leg.Draw()

canvasMain.SaveAs("fres/calibration_cuve_%s.png"%(dicHists["rootfilename"]))

################################################################
##  Everything is over; closing the file
################################################################

rootFits.Write()
rootFits.Close()
rootHists.Close()

print ''
print ''

