#!/usr/bin/env python
import ROOT, nano.analysis.CMS_lumi, json, os, getopt, sys
from nano.analysis.histoHelper import *
from ROOT import TLorentzVector
#import DYestimation
ROOT.gROOT.SetBatch(True)

json_used = 'Golden'
datalumi = 35.9#35.9fb-1
CMS_lumi.lumi_sqrtS = "%.1f fb^{-1}, #sqrt{s} = 13 TeV 25ns "%(datalumi)
datalumi = datalumi*1000
version = os.environ['CMSSW_VERSION']

mcfilelist = [
               'TT_powheg',
               'WJets',
               "SingleTop_tW",
               "SingleTbar_tW",
               'ZZ',
               'WW',
               'WZ',
               'DYJets',
               'DYJets_10to50',
             ]
rdfilelist = ['ch_me','ch_ee','ch_mm']

rootfileDir = "/xrootd/store/user/dayoung/nanoAOD/results_merged/topmass_"

channel_name = ['MuEl', 'ElEl', 'MuMu']

datasets = json.load(open("%s/src/nano/analysis/data/dataset/dataset.json" % os.environ['CMSSW_BASE']))

#defalts
step = 5
channel = 3
cut = 'dilep.M() > 20'
#weight = 'genweight*puweight'
#weight = 'genweight*puweight*mueffweight*eleffweight'
weight = 'genweight*puweight*mueffweight*eleffweight*btagweight'
#weight = 'genweight'
plotvar = 'dilep.M()'
binning = [60, 20, 320]
x_name = 'mass [GeV]'
y_name = 'Events'
dolog = False
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


#cut define
stepch_tcut = 'step>=%i&&channel==%i'%(step,channel)
if channel == 0:
    stepch_tcut = 'step>=%i'%(step)

tcut = '(%s&&%s)*(%s)'%(stepch_tcut,cut,weight)
rd_tcut = '%s&&%s'%(stepch_tcut,cut)
rdfname = rootfileDir + rdfilelist[0] +".root"

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


#Drawing plots on canvas
var = plotvar.split(',')[0]
var = ''.join(i for i in var if not i.isdigit())
var = var.replace('.','_').lower()
var = var.replace('()','')

outfile = "%s_s%d_%s.png"%(channel_name[channel-1],step,var)
if channel == 0: outfile = "Dilepton_s%d_%s.png"%(step,var)
canv = drawTH1(outfile, CMS_lumi, mchistList, rdhist, x_name, y_name,dolog,True,0.5)
canv.SaveAs(outfile)        
print outfile
