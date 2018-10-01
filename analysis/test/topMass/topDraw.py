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
            "ZZ",
            "WW",
            "WZ",
            "DYJets",
            'DYJets_10to50',
            "WToLNu_0J",
            "WToLNu_1J",
            "WToLNu_2J",
            #"WJetsToLNu_HT-100To200",
            #"WJetsToLNu_HT-2500ToInf",
            #"WJetsToLNu_HT-70To100",
            #"WJetsToLNu_HT-1200To2500",
            #"WJetsToLNu_HT-400To600",
            #"WJetsToLNu_HT-800To1200",
            #"WJetsToLNu_HT-200To400",
            #"WJetsToLNu_HT-600To800",    
            #"QCD_Pt-20to30_EMEnriched",
            #"QCD_Pt-30to50_EMEnriched",
            #"QCD_Pt-50to80_EMEnriched",
            #"QCD_Pt-80to120_EMEnriched",
            #"QCD_Pt-120to170_EMEnriched",
            #"QCD_Pt-170to300_EMEnriched",
            #"QCD_Pt-300toInf_EMEnriched",
        "QCD_Pt-15to20_MuEnriched",
        "QCD_Pt-20to30_MuEnriched",
        "QCD_Pt-30to50_MuEnriched",
        "QCD_Pt-50to80_MuEnriched",
        "QCD_Pt-80to120_MuEnriched",
        "QCD_Pt-120to170_MuEnriched",
        "QCD_Pt-170to300_MuEnriched",
        "QCD_Pt-300to470_MuEnriched",
        "QCD_Pt-470to600_MuEnriched",
        "QCD_Pt-600to800_MuEnriched",
        "QCD_Pt-800to1000_MuEnriched",
        "QCD_Pt-1000toInf_MuEnriched", 
        "SingleTop_tW",
        "SingleTbar_tW",
        "SingleTop_t-channel",
        "SingleTbar_t-channel",
        "SingleTop_s-channel",
         "TT_powheg",
            #'TT_powheg_mtop1665',
            #'TT_powheg_mtop1695',
            #'TT_powheg_mtop1715',
            #'TT_powheg_mtop1735',
            #'TT_powheg_mtop1755',
            #'TT_powheg_mtop1785',
             ]
#rdfilelist = ['ch_fullme','ch_fullee','ch_fullmm']
rdfilelist = ['ch_e','ch_m']

rootfileDir = "/xrootd/store/user/yyoun/nanoAOD/slmass_8142/results_merged/tth2mu_"
#rootfileDir = "/xrootd/store/user/dayoung/nano/test/results_merged/topmass_"

channel_name = ['El+Jets','Mu+Jets']

datasets = json.load(open("/cms/scratch/yyoun/nanoAOD/src/nano/nanoAOD/data/dataset/dataset.json"))
#% os.environ['CMSSW_BASE']))

#defalts
step =4
#cut = 'vecSumLepSV.M()>0&&ncmeson>0&&tri!=0&&cmeson_pdgId==421&&cmeson_lxy>0.1&&cmeson_l3D>0'
#cut = 'ncmeson>0&&tri!=0&&cmeson_pdgId==443&&cmeson_lxy>0.1&&cmeson_l3D>0.15'
#cut = 'tri !=0&&bjet.Pt()>30&&dilep.Pt()>20&&jet2.Pt()>30&&jet1.Pt()>30'


#=#MU
channel = 2
#PDGID : 421,413,443
cut = "trig_m!=0"#'njet >3 && nbjet >1'# && nbjet >0'# && cme_tmva_bdtg>0.7'# njet >=4 && nbjet == 2'# && had_sumM>0'#weight = 'genweight*puweight*mueffweight*tri'#*eleffweight*topPtWeight'
weight = 'btagweight*genweight*puweight*mueffweight*tri'#*eleffweight*topPtWeight'

#EL

#channel = 1
#cut = 'trig_e != 0 '#&& cme_pdgId == 421'
#weight = 'genweight*puweight*eleffweight*tri'#*eleffweight*topPtWeight'

#cut = 'tri != 0 '
#cut = 'nbjet<4&&njet<7&&tri !=0&&(lep1.Eta()<1.4442||1.566<lep1.Eta())&&(lep2.Eta()<1.4442||1.566<lep2.Eta())'
#cut = 'tri !=0'

#x_name = 'number of vertex'
#plotvar = 'PV_npvs'
#binning = [40, 0, 40]

#x_name = 'Invariant mass [GeV]'
#plotvar = 'lep.M()'
#binning = [10, -0.2, 0.2]

#x_name = 'Hadron mass [GeV]'
#plotvar = "had.M()"
#binning = [50, 1.6, 2.2]
#x_name = 'Hadron  mass [GeV]'
#plotvar = "had.M()"
#binning = [50, 1.9, 2.1]
#x_name = 'Hadron  mass [GeV]'
#plotvar = "had.M()"
#binning = [20, 2.9, 3.3]


################BASIC PLOT####################
x_name = 'Missing Energy [GeV]'
plotvar = 'met'
binning = [20, 0, 200]

#x_name = 'p_{T}^{Mu} [GeV]'
#x_name = 'p_{T}^{El} [GeV]'
#plotvar = 'lep.Pt()'
#binning = [20, 0, 200]

#x_name = '#eta'
#plotvar = 'lep.Eta()'
#binning = [25, -2.4, 2.4]

#x_name = '#phi'
#plotvar = 'lep.Phi()'
#binning = [30, -3, 3]
#x_name = 'Jet Multiplicity'
#x_name = 'Jet Multiplicity'
#plotvar = 'njet'
#binning = [10, 2, 12]
#
#x_name = 'b Jet Multiplicity'
#plotvar = 'nbjet'
#binning = [6, 0, 6]
#####################################

#plotvar = 'vecSumLepSV.M()'
#plotvar = 'jpsi.M()'
#binning = [60,2.8,3.4]
#binning = [8, 0, 8]
#x_name = 'Invariant mass_{J/#psi+L} [GeV]'
#x_name = 'mass_{J/#psi} [GeV]'
#x_name = 'b-jet Multiplicity'
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
#if not os.path.exists('./DYFactor.json'):
#    DYestimation.printDYFactor(rootfileDir, tname, datasets, datalumi, cut, weight, rdfilelist)#This will create 'DYFactor.json' on same dir.
#dyratio=json.load(open('./DYFactor.json'))

mchistList = []
for imc,mcname in enumerate(mcfilelist):
    data = findDataSet(mcname, datasets)
    if data is None:
        print mcname
    scale = datalumi*data["xsec"]
    print"scale: %s " %(scale)
    colour = data["colour"]
    title = data["title"]
 #   if ('DYJets' in mcname) and channel!=0:
 #       scale = scale*dyratio[channel][step]

    rfname = rootfileDir + mcname +".root"
    print rfname
    tfile = ROOT.TFile(rfname)
    wentries = tfile.Get("nevents").Integral()
    scale = scale/wentries
    if ('15to20_Mu' in mcname): scale= scale*0.003 
    if ('20to30_Mu' in mcname): scale= scale*0.0053 
    if ('30to50_Mu' in mcname): scale= scale*0.01182
    if ('50to80_Mu' in mcname): scale= scale*0.02276
    if ('80to120_Mu' in mcname): scale= scale*0.03844
    if ('120to170_Mu' in mcname): scale= scale*0.05362
    if ('170to300_Mu' in mcname): scale= scale*0.07335
    if ('300to470_Mu' in mcname): scale= scale*0.10196
    if ('470to600_Mu' in mcname): scale= scale*0.071
    if ('600to800_Mu' in mcname): scale= scale*0.13412 
    if ('800to1000_Mu' in mcname): scale= scale*0.144552
    if ('1000toInf_Mu' in mcname): scale= scale*0.15544
    #if ('20to30_EM' in mcname): scale= scale*0.0053 
    #if ('30to50_EM' in mcname): scale= scale*0.01182
    #if ('50to80_EM' in mcname): scale= scale*0.01182
    #if ('80to120_EM' in mcname): scale= scale*0.01182
    #if ('120to170_EM' in mcname): scale= scale*0.01182
    #if ('170to300_EM' in mcname): scale= scale*0.01182
    #if ('300toInf_EM' in mcname): scale= scale*0.15544
    print "wentries:%s, scale:%s"%(wentries,scale)
    mchist = makeTH1(rfname, tname, title, binning, plotvar, tcut, scale)    
    mchist.SetFillColor(colour)
    mchist.SetLineColor(colour)
    mchistList.append(mchist)

    remchist = makeTH1(rfname, tname, title, binning, plotvar, tcut, scale)

#if dolog == True:
#    minp = 0.05
#   maxp = 1000000000



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

rdhist.SetMinimum(0)
rdhist.SetMaximum()

#Drawing plots on canvas
var = plotvar.split(',')[0]
var = ''.join(i for i in var if not i.isdigit())
var = var.replace('.','_').lower()
var = var.replace('()','')

outfile = "%s_8142_%s_s%d.png"%(channel_name[channel-1],var,step)
if channel == 0: outfile = "SemiLepton_s%d_%s.png"%(step,var)
#canv = drawTH1(outfile, CMS_lumi, mchistList, rdhist, x_name, y_name,dolog,True,0.5)
canv = drawTH1(outfile, CMS_lumi, mchistList, rdhist, x_name, y_name,dolog)

mainPad = canv.GetPrimitive("mainPad")
mainPad.cd()
mainPad.GetPrimitive("data").Draw("esamex0")
extraText(canv, [0.4,0.85], outfile.split("_")[0])
canv.Update()

fIntMC = 0.0

for mchist in mchistList:
  fIntMC = fIntMC + mchist.Integral()

extraText(
    canv,
    [0.27,0.90],
    "Data: {} MC:{}".format(int(rdhist.Integral()), int(fIntMC)))
print "Integral (MC) : ", step, channel, fIntMC
print "Integral (data) : ", step, channel, rdhist.Integral()




canv.SaveAs(outfile)        
print outfile
