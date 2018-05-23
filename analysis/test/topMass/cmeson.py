#!/usr/bin/env python
import ROOT, nano.analysis.CMS_lumi, json, os, getopt, sys
from nano.analysis.histoHelper import *
from ROOT import TLorentzVector
import DYestimation
ROOT.gROOT.SetBatch(True)

json_used = 'Golden'
datalumi = 35.9#35.9fb-1
CMS_lumi.lumi_sqrtS = "%.1f fb^{-1}, #sqrt{s} = 13 TeV 25ns "%(datalumi)
datalumi = datalumi*1000
version = os.environ['CMSSW_VERSION']

mcfilelist = [
               'ZZ',
               'WW',
               'WZ',
               'WJets',
               "SingleTop_tW",
               "SingleTbar_tW",
               'DYJets',
               'DYJets_10to50',
               'TT_powheg',
             ]
#rdfilelist = []
rdfilelist = ['ch_me','ch_fullee','ch_mm']

rootfileDir = "/xrootd/store/user/jdj0715/nanoAOD/test13/results_merged/topmass_"
#rootfileDir = "/xrootd/store/user/dayoung/nano/test/results_merged/topmass_"

channel_name = ['MuEl', 'ElEl', 'MuMu']

datasets = json.load(open("%s/src/nano/nanoAOD/data/dataset/dataset.json" % os.environ['CMSSW_BASE']))

#default
step = 5
channel = 0
#cut = 'vecSumLepSV.M()>0&&ncmeson>0&&tri!=0&&cmeson_pdgId==421&&cmeson_lxy>0.1&&cmeson_l3D>0'
#cut = 'ncmeson>0&&tri!=0&&cmeson_pdgId==443&&cmeson_lxy>0.1&&cmeson_l3D>0.15'
#cut = 'tri !=0&&bjet.Pt()>30&&dilep.Pt()>20&&jet2.Pt()>30&&jet1.Pt()>30'
#cut = 'tri !=0&&abs(cme_pdgId)==443&&cme_tmva_bdtg>0.5'
#cut = 'nbjet<4&&njet<7&&tri !=0&&(lep1.Eta()<1.4442||1.566<lep1.Eta())&&(lep2.Eta()<1.4442||1.566<lep2.Eta())'
cut = 'tri !=0&&cme_pdgId==421&&cme_tmva_bdtg>0.3'
weight = 'genweight*puweight*eleffweight*mueffweight*tri'#*topPtWeight'
#plotvar = 'vecSumLepSV.M()'
#plotvar = 'njet'
plotvar = 'd0.M()'
#binning = [20,2.8,3.4]
#binning = [20,1.7,2.2]
#binning = [60,1,1000]
#binning = [20, 0, 200]
binning = [60, 1.8, 2.2]
#binning = [60, 20, 320]
x_name = 'Invariant mass(ll) [GeV]'
#x_name = 'Missing Energy [GeV]'
#x_name = 'Invariant mass of J/#psi [GeV]'
#x_name = 'mass_{J/#psi} [GeV]'
#x_name = 'Jet Multiplicity'
y_name = 'Events'
dolog = False
#dolog = True
tname = "event"
#get input
try:
    opts, args = getopt.getopt(sys.argv[1:],"hdc:w:b:p:x:y:j:a:s:",["cut","weight","binning","plotvar","x_name","y_name","json_used","dolog", "channel", "step"])
except getopt.GetoptError:          
    print 'Usage : ./topDraw.py -c <cut> -w <weight> -b <binning> -p <plotvar> -x <x_name> -y <y_name> -j <json_used> -d <dolog>'
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print 'Usage : ./topDraw.py -c <cut> -w <weight> -b <binning> -p <plotvar> -x <x_name> -y <y_name>  -j <json_used> -d <dolog>'
        sys.exit()
    elif opt in ("-c", "--cut"):
        cut = arg
    elif opt in ("-a", "--channel"):
        channel = int(arg)
    elif opt in ("-s", "--step"):
        step = int(arg)
    elif opt in ("-w", "--weight"):
        weight = arg
    elif opt in ("-b", "--binning"):
        binning = eval(arg)
    elif opt in ("-p", "--plotvar"):
        plotvar = arg
    elif opt in ("-x", "--x_name"):
        x_name = arg
    elif opt in ("-y", "--y_name"):
        y_name = arg
    elif opt in ("-j", "--json_used"):
        json_used = arg
    elif opt in ("-d", "--dolog"):
        dolog = True
print plotvar, x_name 


sig = 0
bg = 0

#cut define
stepch_tcut = 'step>=%i&&channel==%i'%(step,channel)
if channel == 0:
    stepch_tcut = 'step>=%i'%(step)

tcut = '(%s&&%s)*(%s)'%(stepch_tcut,cut,weight)
rd_tcut = '%s&&%s'%(stepch_tcut,cut)

#DYestimation
if not os.path.exists('./DYFactor.json'):
    DYestimation.printDYFactor(rootfileDir, tname, datasets, datalumi, cut, weight, rdfilelist)#This will create 'DYFactor.json' on same dir.
dyratio=json.load(open('./DYFactor.json'))

mchistList = []
for imc,mcname in enumerate(mcfilelist):
    data = findDataSet(mcname, datasets)
    scale = datalumi*data["xsec"]
    print"scale: %s " %(scale)
    colour = data["colour"]
    title = data["title"]
    if ('DYJets' in mcname) and channel!=0:
        scale = scale*dyratio[channel][step]

    rfname = rootfileDir + mcname +".root"
    print rfname
    tfile = ROOT.TFile(rfname)
    wentries = tfile.Get("nevents").Integral()
    scale = scale/wentries
    print "wentries:%s, scale:%s"%(wentries,scale)
    mchist = makeTH1(rfname, tname, title, binning, plotvar, tcut, scale)    
    mchist.SetFillColor(colour)
    mchist.SetLineColor(colour)
    mchistList.append(mchist)

    remchist = makeTH1(rfname, tname, title, binning, plotvar, tcut, scale)
    if "TT_powheg" in mcname:
      sig += remchist.Integral(remchist.FindBin(120), remchist.FindBin(130))
      TTH = remchist.Integral(remchist.FindBin(120), remchist.FindBin(130))
      print "TT_powheg"
    elif "DYJets" in mcname:
      sig += remchist.Integral(remchist.FindBin(120), remchist.FindBin(130))
      print "dy"
    elif "DYJets_10to50" in mcname:
      sig += remchist.Integral(remchist.FindBin(120), remchist.FindBin(130))
      print "dy50"
    elif "ZZ" in mcname:
      sig += remchist.Integral(remchist.FindBin(120), remchist.FindBin(130))
      print "ZZ"
    elif "WW" in mcname:
      sig += remchist.Integral(remchist.FindBin(120), remchist.FindBin(130))
      print "WW"
    elif "WZ" in mcname:
      sig += remchist.Integral(remchist.FindBin(120), remchist.FindBin(130))
      print "WZ"
    else:
      bg += remchist.Integral(remchist.FindBin(120), remchist.FindBin(130))
      print mcname
    if "WJets" in mcname:
      TTJets = remchist.Integral(remchist.FindBin(120), remchist.FindBin(130))
    if "SingleTop_tW" in mcname:
      SingleTop = remchist.Integral(remchist.FindBin(120), remchist.FindBin(130))
    if "SingleTbar_tW" in mcname:
      SingleTbar = remchist.Integral(remchist.FindBin(120), remchist.FindBin(130))


if dolog == True:
   minp = 0.05
   maxp = 1000000000






#data histo
if channel != 0:
    rfname = rootfileDir + rdfilelist[channel-1] +".root"
    rdhist = makeTH1(rfname, tname, 'data', binning, plotvar, rd_tcut)
else:
    rdhist = mchistList[0].Clone()
    rdhist.Reset()
    for i, rdfile in enumerate(rdfilelist):
        rfname = rootfileDir + rdfile +".root"
        rdtcut = 'channel==%d&&%s&&%s'%((i+1),stepch_tcut,cut)
        rdhist.Add(makeTH1(rfname, tname, 'data', binning, plotvar, rdtcut))
'''
rdhist.SetMinimum(minp)
rdhist.SetMaximum(maxp)
'''
#Drawing plots on canvas
var = plotvar.split(',')[0]
var = ''.join(i for i in var if not i.isdigit())
var = var.replace('.','_').lower()
var = var.replace('()','')

outfile = "%s_s%d_%s.png"%(channel_name[channel-1],step,var)
if channel == 0: outfile = "Dilepton_s%d_%s.png"%(step,var)
canv = drawTH1(outfile, CMS_lumi, mchistList, rdhist, x_name, y_name,dolog,True,0.5)
#canv = drawTH1(outfile, CMS_lumi, mchistList, rdhist, x_name, y_name,dolog)

mainPad = canv.GetPrimitive("mainPad")
mainPad.cd()
mainPad.GetPrimitive("data").Draw("esamex0")
extraText(canv, [0.4,0.85], outfile.split("_")[0])
canv.Update()

fIntMC = 0.0

for mchist in mchistList:
  fIntMC = fIntMC + mchist.Integral()

print "Integral (MC) : ", step, channel, fIntMC
print "Integral (data) : ", step, channel, rdhist.Integral()

canv.SaveAs(outfile)        
print outfile
