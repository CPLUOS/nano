#!/usr/bin/env python
import ROOT, nano.analysis.CMS_lumi, json, os, getopt, sys, copy
from nano.analysis.histoHelper import *
from ROOT import TLorentzVector
import DYestimation
ROOT.gROOT.SetBatch(True)

rootfileDir = "/xrootd/store/user/jdj0715/nanoAOD/test_1128_1/results_merged/topmass_"
datalumi = 35.9 # Run2016 B & C & D & E & F & G, v8-0-2 (H is not yet)
CMS_lumi.lumi_sqrtS = "%.1f fb^{-1}, #sqrt{s} = 13 TeV"%(datalumi)
datalumi = datalumi*1000 # due to fb

topMassList = [
  #  'TT_powheg_mtop1665',
  #  'TT_powheg_mtop1695',
  #  'TT_powheg_mtop1715',
  #  'TT_powheg_mtop1735',
  #  'TT_powheg_mtop1755',
  #  'TT_powheg_mtop1785',
    'TT_powheg'
    ] # it will be filled fully
#
mcfilelist = [
            "WJets",
            "SingleTbar_tW",
            "SingleTop_tW",
            "ZZ",
            "WW",
            "WZ",
            "DYJets",
            'DYJets_10to50'
             ]
channel_name = ['Combined', 'MuEl', 'ElEl','MuMu']
dMassNomial = 172.50
datasets = json.load(open("%s/src/nano/nanoAOD/data/dataset/dataset.json" % os.environ['CMSSW_BASE']))

## defalts
step = 5
channel = 0
cut = 'tri!=0&&cme_tmva_bdtg>0.6'
#cut = 'tri!=0&&abs(cme_pdgId==421)&&cme_tmva_bdtg>0.2'
weight = 'genweight*puweight*eleffweight*mueffweight*btagweight*tri'
weightTopPT = ''
binning = [25, 0, 250]
plotvar = 'd0_lepSV_dRM_false'
#binning = [50, 1.82, 1.9]
#plotvar = 'had.M()'
strType = ""
x_name = 'M_{lepton + D^{0}} [GeV]'
y_name = 'Events'
dolog = False
overflow = False
binNormalize = False
suffix = ''
strTypeSuffix = ""

################################################################
## get input
################################################################
try:
    opts, args = getopt.getopt(sys.argv[1:],"hdnoa:p:s:c:w:b:t:x:y:f:",["channel","plotvar","step","cut","weight","binning","type","binNormalize","overflow","x_name","y_name","dolog","suffix"])
except getopt.GetoptError:          
    print 'Usage : ./massPlot.py -c <cut> -w <weight> -b <binning> -p <plotvar> -x <x_name> -y <y_name> -d <dolog> -f <suffix>'
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print 'Usage : ./topDraw.py -c <cut> -w <weight> -b <binning> -p <plotvar> -x <x_name> -y <y_name> -d <dolog> -f <suffix>'
        sys.exit()
    elif opt in ("-a", "--channel"):
        channel = int(arg)
    elif opt in ("-p", "--plotvar"):
        plotvar = arg
    elif opt in ("-s", "--step"):
        step = int(arg)
    elif opt in ("-c", "--cut"):
        #cut = arg
        cut = "%s&&%s"%(cut,arg)
    elif opt in ("-w", "--weight"):
        weight = arg
    elif opt in ("-b", "--binning"):
        binning = eval(arg)
    elif opt in ("-t", "--type"):
        strType = arg
    elif opt in ("-x", "--x_name"):
        x_name = arg
    elif opt in ("-y", "--y_name"):
        y_name = arg
    elif opt in ("-d", "--dolog"):
        dolog = True
    elif opt in ("-o", "--overflow"):
        overflow = True
    elif opt in ("-n", "--binNormalize"):
        binNormalize = True
    elif opt in ("-f", "--suffix"):
        suffix = "_"+arg

#tname = "cattree/nom"
tnameData = "event"

################################################################
## Read type
################################################################
stepch_tcut = 'step>=%i%s'%(step, "&&channel==%i"%(channel) if channel != 0 else "")
tcutonly = '%s%s'%(stepch_tcut, "&&" + cut if cut != "" else "")
tcut = '(%s)*(%s)'%(tcutonly, weight + weightTopPT)
print "TCut =",tcut

## namming
if len(binning) <= 3:
  num = (binning[2]-binning[1])/float(binning[0])
  if num != 1:
    unit = "["+x_name.split('[')[1] if x_name.endswith(']') else ""
    y_name = y_name + "/%g%s"%(num,unit)

## DYestimation
#if not os.path.exists('./DYFactor.json'):
#  DYestimation.printDYFactor(rootfileDir, tname, datasets, datalumi, cut, weight, rdfilelist)# <------ This will create 'DYFactor.json' on same dir.
#DYestimation
if not os.path.exists('./DYFactor.json'):
    DYestimation.printDYFactor(rootfileDir, tname, datasets, datalumi, cut, weight, rdfilelist)#This will create 'DYFactor.json' on same dir.
dyratio=json.load(open('./DYFactor.json'))
#
## Initializing the result root file
strFilename = "invMass_%s%s%s"%(plotvar,suffix,strTypeSuffix)
outMassHist = ROOT.TFile.Open("invmass/" + strFilename + ".root","RECREATE")
dicListHist = {"rootfilename":strFilename + ".root", 
  "x_name":x_name, "y_name":y_name, "binning":binning, "Gen":0}

################################################################
## Saving MC histograms for backgrounds
################################################################
mchistList = []
for i, mcname in enumerate(mcfilelist):
  data = findDataSet(mcname, datasets)
  scale = datalumi*data["xsec"]
  colour = data["colour"]
  title = data["title"]
  if ('DYJets' in mcname) and channel!=0:
        scale = scale*dyratio[channel][step]

  rfname = rootfileDir + mcname +".root"
  tfile = ROOT.TFile(rfname)
  wentries = tfile.Get("nevents").Integral()
  scale = scale/wentries
  print "Bkg scales : ", mcname, scale
    
  mchist = makeTH1(rfname, tnameData, title, binning, plotvar, tcut, scale)
  mchist.SetLineColor(colour)
  mchist.SetFillColor(colour)
  mchistList.append(mchist)
  
  mchist.SetName("hist_bkg_" + mcname)
  outMassHist.cd()
  mchist.Write()

  dicListHist[mchist.GetName()] = {"type":"bkg_part"}
  
#overflow
if overflow:
  nbin = binning[0] if len(binning) == 3 else len(binning)-1
  for hist in mchistList:
    hist.SetBinContent(nbin, hist.GetBinContent(nbin+1))

#bin normalize
if binNormalize and len(binning)!=3:
  for hist in mchistList:
    for i in range(len(binning)):
      hist.SetBinContent(i, hist.GetBinContent(i)/hist.GetBinWidth(i))
      hist.SetBinError(i, hist.GetBinError(i)/hist.GetBinWidth(i))

hs_bkg = ROOT.THStack("bkg_hs","bkg_hs")
for hist in mchistList :
  hs_bkg.Add(hist)

hs_bkg.Draw()
bkgs = hs_bkg.GetStack().Last()

bkgs.SetName("hist_bkg")
bkgs.Draw()

outMassHist.cd()
bkgs.Write()

dicListHist[bkgs.GetName()] = {"type":"bkg"}

print "bkg entries: ",bkgs.GetEntries()

################################################################
##  Saving TT samples with side-band
################################################################
"""
arrInfoSideBandWin = [
  {"name":"center", "mv":0.0, }, 
]

for dMoveWin in [-1.0, 0.0, 1.0]:
  dicInfoWin = 
"""

for topMass in topMassList :
  massValue = topMass.split("mtop")[-1] if ( topMass.find("mtop") != -1 ) else "nominal"
  sum_hs =  hs_bkg.Clone()
  data = findDataSet(topMass, datasets)
  scale = datalumi*data["xsec"]
  colour = data["colour"]
  title = data["title"]

  dMassCurr = int(massValue) * 0.1 if massValue != "nominal" else dMassNomial

  rfname = rootfileDir + topMass + ".root"
  tfile = ROOT.TFile(rfname)
  wentries = tfile.Get("nevents").Integral()
  scale = scale/wentries
  print topMass, scale, wentries, colour, title
  
  #for dicInfoWin in arrInfoSideBandWin:
  mchist = makeTH1(rfname, tnameData, title, binning, plotvar, tcut, scale)
  mchist.SetLineColor(colour)
  mchist.SetFillColor(colour)
  print "topmass hsit : ",mchist.Integral()
  # -- Overflow
  if overflow:
    nbin = binning[0] if len(binning) == 3 else len(binning)-1
    mchist.SetBinContent(nbin, mchist.GetBinContent(nbin+1))

  # -- Bin normalize
  if binNormalize and len(binning) != 3:
    for i in range(len(binning)):
      mchist.SetBinContent(i, mchist.GetBinContent(i) / mchist.GetBinWidth(i))
      mchist.SetBinError(i, mchist.GetBinError(i) / mchist.GetBinWidth(i))

  # -- Getting the plot 
  sum_hs.Add( mchist )
  masshist = sum_hs.GetStack().Last()
 
  # -- Saving results into the root file
  masshist.Draw()
  print masshist.GetEntries()
  
  outMassHist.cd()
  mchist.SetName("hist_TT_onlytt_%s"%(massValue))
  mchist.Write()
  
  dicListHist[mchist.GetName()] = {"type":"TT_onlytt", "mass":dMassCurr}
  
  outMassHist.cd()
  masshist.SetName("hist_TT_withbkg_%s"%(massValue))
  masshist.Write()
  
  dicListHist[masshist.GetName()] = {"type":"TT_withbkg", "mass":dMassCurr}

if "correctM" in plotvar:
  dicListHist[ "Gen" ] = 1

if overflow:
  nbin = binning[0] if len(binning) == 3 else len(binning)-1
#    rdhist.SetBinError(i, rdhist.GetBinError(i)/rdhist.GetBinWidth(i))
#
outMassHist.cd()

outMassHist.Write()
outMassHist.Close()

################################################################
##  Saving informations about histograms
################################################################

fileDicHist = open("invmass/" + strFilename + ".json", "w")
fileDicHist.write(json.dumps(dicListHist))
fileDicHist.close()

print ''
print ''
