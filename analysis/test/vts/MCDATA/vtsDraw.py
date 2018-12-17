#!/usr/bin/env python
import ROOT, nano.analysis.CMS_lumi, json, os, getopt, sys
from nano.analysis.histoHelper import *
from ROOT import *
#import DYestimation
ROOT.gROOT.SetBatch(True)

json_used = 'Golden'
datalumi = 35.9#35.9fb-1
CMS_lumi.lumi_sqrtS = "%.1f fb^{-1}, #sqrt{s} = 13 TeV 25ns "%(datalumi)
datalumi = datalumi*1000
version = os.environ['CMSSW_VERSION']

mcFiles_topMass = [
               'ZZ',
               'WW',
               'WZ',
               'WJets',
               "SingleTop_tW",
               "SingleTbar_tW",
               'DYJets',
               'DYJets_10to50',
               'TT_powheg',
             ]  # MC used for topMass

mcFiles_vts_dyjet = [
                     'DYJets',
                     'DYJets_10to50',
                     'DYJets_MG',
#                     'DYJets_MG_10to50',
                    ]
mcFiles_vts_ttbar = [
                     'TT_powheg', #
                     'TTJets_aMC', #
                     'TTJets_DiLept_MG', #
                     'TTJets_DiLept',
#                     'TTJets_aMC_ext1', #
                     'TT_powheg_mtop1665',
                     'TT_powheg_mtop1695',
                     'TT_powheg_mtop1715',
                     'TT_powheg_mtop1735',
                     'TT_powheg_mtop1755',
                     'TT_powheg_mtop1785'
                    ]
mcFiles_vts_single = [
                      'SingleTop_tW',
                      'SingleTbar_tW',
                      'SingleTop_tW_noHadron',
                      'SingleTbar_tW_noHadron',
                      'SingleTop_t-channel', #
                      'SingleTbar_t-channel',
                      'SingleTop_s-channel'
                     ]
mcFiles_vts_diboson = [
                       'ZZ',
                       'WW',
                       'WZ',
                      ]
mcFiles_vts_diboson_2 = [
                         'ZZTo4L_powheg',
                         'WWTo2L2Nu',
                         'WZTo3LNu_amcatnlo',
                         'WZTo2LQQ',
                         'ZZTo2L2Nu',
                         'ZZTo2L2Q'
                        ]
mcFiles_vts_triboson = [
                        'WWW',
                        'WWZ',
                        'WZZ',
                        'ZZZ'
                       ]
mcFiles_vts_wjets = [
                     'WJets',
                     'WJetsToLNu_HT-100To200_ext1',
                     'WJetsToLNu_HT-200To400_ext1',
                     'WJetsToLNu_HT-400To600_ext1',
                     'WJetsToLNu_HT-600To800_ext1',
                     'WJetsToLNu_HT-800To1200_ext1',
                     'WJetsToLNu_HT-1200To2500_ext1',
                     'WJetsToLNu_HT-2500ToInf_ext1',
                     'WJetsToLNu_HT-70To100',
                     'WJetsToLNu_HT-100To200',
                     'WJetsToLNu_HT-200To400',
                     'WJetsToLNu_HT-400To600',
                     'WJetsToLNu_HT-600To800',
                     'WJetsToLNu_HT-800To1200',
                     'WJetsToLNu_HT-1200To2500',
                     'WJetsToLNu_HT-2500ToInf'
                    ]
mcFiles_vts_2 = [
                  'GG_HToMuMu',
                  'VBF_HToMuMu',
                  'WMinusH_HToMuMu',
                  'ZH_HToMuMu',

                  'TTZToLLNuNu',
                  'TTWJetsToLNu',

                  'tZq_ToLL',

                  'DYToLL_2J', #

                  'GluGlu_ZZTo4mu',
                  'GluGlu_ZZTo4e',
                  'GluGlu_ZZTo2e2mu',

                  'ttH_HToMuMu',

                  'tsW',

                  'VBF-C1N2_WZ',
                  'VBF-C1N2_leptonicDecays',
                  'VBF-C1N2_tauDecays',

                  'DYJetsToLL_M-50_HT-70to100',
                  'DYJetsToLL_M-50_HT-100to200',
                  'DYJetsToLL_M-50_HT-100to200_ext1',
                  'DYJetsToLL_M-50_HT-200to400',
                  'DYJetsToLL_M-50_HT-200to400_ext1',
                  'DYJetsToLL_M-50_HT-400to600',
                  'DYJetsToLL_M-50_HT-400to600_ext1',
                  'DYJetsToLL_M-50_HT-600to800',
                  'DYJetsToLL_M-50_HT-800to1200',
                  'DYJetsToLL_M-50_HT-1200to2500', #
                  'DYJetsToLL_M-50_HT-2500toInf',

                  'EWKZ2Jets_ZToLL',

                  'EWKWPlus2Jets_WToLNu',
                  'EWKWMinus2Jets_WToLNu',

                  'QCD_HT50to100',
                  'QCD_HT100to200', #
                  'QCD_HT200to300',
                  'QCD_HT300to500',
                  'QCD_HT500to700',
                  'QCD_HT700to1000',
                  'QCD_HT1000to1500',
                  'QCD_HT1500to2000',
                  'QCD_HT2000toInf',

                  'ZJetsTONuNu_HT-100To200',
                  'ZJetsTONuNu_HT-200To400',
                  'ZJetsTONuNu_HT-400To600',
                  'ZJetsTONuNu_HT-600To800',
                  'ZJetsTONuNu_HT-800To1200',
                  'ZJetsTONuNu_HT-1200To2500',
                  'ZJetsTONuNu_HT-2500ToInf',

                  'EWKZ2Jets_ZToNuNu',
                  'WWJJToLNuLNu_EWK',
                  'WZJJ_EWK',

                  'WToLNu_0J', #
                  'WToLNu_1J', #
                  'WToLNu_2J', #

                  'QCD_Pt-15to20_MuEnriched',
                  'QCD_Pt-20to30_MuEnriched',
                  'QCD_Pt-30to50_MuEnriched',
                  'QCD_Pt-50to80_MuEnriched',
                  'QCD_Pt-80to120_MuEnriched',
                  'QCD_Pt-120to170_MuEnriched',
                  'QCD_Pt-170to300_MuEnriched',
                  'QCD_Pt-300to470_MuEnriched',
                  'QCD_Pt-470to600_MuEnriched',
                  'QCD_Pt-600to800_MuEnriched',
                  'QCD_Pt-800to1000_MuEnriched',

                  'QCD_Pt-20to30_EMEnriched',
                  'QCD_Pt-30to50_EMEnriched',
                  'QCD_Pt-50to80_EMEnriched',
                  'QCD_Pt-80to120_EMEnriched', #
                  'QCD_Pt-120to170_EMEnriched', #
                  'QCD_Pt-170to300_EMEnriched',
                  'QCD_Pt-300toInf_EMEnriched'
                ]


#mcfilelist = mcFiles_vts_dyjet + mcFiles_vts_ttbar + mcFiles_vts_single+ mcFiles_vts_diboson + mcFiles_vts_diboson_2 + mcFiles_vts_triboson + mcFiles_vts_wjets + mcFiles_vts_2

mcfilelist = mcFiles_vts_diboson + mcFiles_vts_diboson_2 + mcFiles_vts_wjets + mcFiles_vts_single + mcFiles_vts_dyjet + mcFiles_vts_ttbar + mcFiles_vts_triboson + mcFiles_vts_2
#mcfilelist = mcFiles_topMass

#mcfilelist = ['tsW']
#mcfilelist =[]

dataFiles = [
             'SingleMuon_Run2016', 'SingleElectron_Run2016',
             'DoubleMuon_Run2016', 'DoubleEG_Run2016',
             'MuonEG_Run2016'
            ]
rdfilelist = [data+period for data in dataFiles for period in ["B","C","D","E","F","G","H"]]
rdfilelist = ['ch_fullme','ch_fullee','ch_fullmm']
#rdfilelist = []

rootfileDir = "/xrootd_user/wjjang/xrootd/sum/run2_2016v5/20181008/"
rdfileDir = "/xrootd_user/wjjang/xrootd/run2_2016v5/"

channel_name = ['MuEl', 'ElEl', 'MuMu']

mc_nev = 0
rd_nev = 0

datasets = json.load(open("%s/src/nano/nanoAOD/data/dataset/dataset.json" % os.environ['CMSSW_BASE']))

#default
step = 4
channel = 2 # 1: MuEl | 2: ElEl | 3: MuMu
#cut = 'vecSumLepSV.M()>0&&ncmeson>0&&tri!=0&&cmeson_pdgId==421&&cmeson_lxy>0.1&&cmeson_l3D>0'
#cut = 'ncmeson>0&&tri!=0&&cmeson_pdgId==443&&cmeson_lxy>0.1&&cmeson_l3D>0.15'
#cut = 'tri !=0&&bjet.Pt()>30&&dilep.Pt()>20&&jet2.Pt()>30&&jet1.Pt()>30'
#cut = 'tri !=0&&abs(cme_pdgId)==443&&cme_tmva_bdtg>0.5'
#cut = 'nbjet<4&&njet<7&&tri !=0&&(lep1.Eta()<1.4442||1.566<lep1.Eta())&&(lep2.Eta()<1.4442||1.566<lep2.Eta())'
cut = 'tri != 0' #'(passedEvent || !passedEvent)'
#weight = 'genweight*puweight*eleffweight*mueffweight*tri'#*topPtWeight'
weight = 'genweight*puweight*mueffweight*tri'#*topPtWeight'

#plotvar = 'njet'
plotvar = 'dilep.M()'
#plotvar = 'dilep.Pt()'
#binning = [20,2.8,3.4]
#binning = [20,1.7,2.2]
#binning = [60,1,1000]
#binning = [20, 0, 200]
#binning = [60, 1.7, 2.0]
binning = [60, 20, 320]
x_name = 'Invariant mass(ll) [GeV]'
#x_name = 'Transverse momentum(ll) [GeV]'
#x_name = 'Missing Energy [GeV]'
#x_name = 'Invariant mass of J/#psi [GeV]'
#x_name = 'mass_{J/#psi} [GeV]'
#x_name = 'Jet Multiplicity'
y_name = 'Events'
dolog = False
dolog = True
tname = "event"
#get input
try:
    opts, args = getopt.getopt(sys.argv[1:],"hdc:w:b:p:x:y:j:a:s:",["cut","weight","binning","plotvar","x_name","y_name","json_used","dolog", "channel", "step"])
except getopt.GetoptError:          
    print 'Usage : ./vtsDraw.py -c <cut> -w <weight> -b <binning> -p <plotvar> -x <x_name> -y <y_name> -j <json_used> -d <dolog>'
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print 'Usage : ./vtsDraw.py -c <cut> -w <weight> -b <binning> -p <plotvar> -x <x_name> -y <y_name>  -j <json_used> -d <dolog>'
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

"""
stepch_tcut = 'step>=%i'%(step)#'step>=%i && channel==%i'%(step,channel)
if channel == 0:
    stepch_tcut = 'step>=%i'%(step)
tcut = '%s'%(stepch_tcut)#'(%s&&%s)*(%s)'%(stepch_tcut,cut,weight)
rd_tcut = '%s'%(stepch_tcut)#'%s&&%s'%(stepch_tcut,cut)
"""
#DYestimation
#if not os.path.exists('./DYFactor.json'):
#    DYestimation.printDYFactor(rootfileDir, tname, datasets, datalumi, cut, weight, rdfilelist)#This will create 'DYFactor.json' on same dir.
#dyratio=json.load(open('./DYFactor.json'))

mchistList = []
for imc,mcname in enumerate(mcfilelist):
    data = findDataSet(mcname, datasets)
    scale = datalumi*data["xsec"]
    print"scale: %s " %(scale)
    colour = data["colour"]
    title = data["title"]
#    if ('DYJets' in mcname) and channel!=0:
#        scale = scale*dyratio[channel][step]

    rfname = rootfileDir + mcname +"_sum.root"
    print rfname
    tfile = ROOT.TFile(rfname)
    wentries = tfile.Get("nevents").Integral()
#    wentries = tfile.Get("event").GetEntries()
    mc_nev += wentries
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
rd_nev = 0
#data histo
rdhistlist = []
if channel != 0: # 1 : me , 2 : ee, 3 : mm
    path = os.listdir("/xrootd_user/wjjang/xrootd/run2_2016v5/")
    for sample in path:
        rd_dir = os.listdir("/xrootd_user/wjjang/xrootd/run2_2016v5/"+sample+"/")
        if channel == 1 :
            if sample.find("SingleElectron") == -1 and sample.find("SingleMuon") == -1 and sample.find("MuonEG") == -1 :
                continue
        if channel == 2 :
            if sample.find("SingleElectron") == -1 and sample.find("DoubleEG") == -1 :
                continue
        if channel == 3 :
            if sample.find("SingleMuon") == -1 and sample.find("DoubleMuon") == -1 :
                continue
        print(" >>>>>>>>>> rd plotting .... : " + sample + " <<<<<<<<<< ")
        for i, rd_file in enumerate(rd_dir) :
            rfname = "/xrootd_user/wjjang/xrootd/run2_2016v5/"+sample+"/"+rd_file
#            print(rfname)
#            rd_nev = rd_nev + TFile(rfname).Get("event").GetEntries()
            if i == 0 :
                rdh = makeTH1(rfname, tname, 'data', binning, plotvar, rd_tcut)
            else :
                rdh.Add(makeTH1(rfname, tname, 'data', binning, plotvar, rd_tcut))
        rdhistlist.append(rdh)
else:
    rdhist = mchistList[0].Clone()
    rdhist.Reset()
    for i, rdfile in enumerate(rdfilelist):
        rfname = rootfileDir + rdfile +"_sum.root"
#        rfname = "/xrootd_user/wjjang/xrootd/run2_2016v5/SingleElectron_Run2016B/rd_nanoAOD_100.root"
        rdtcut = 'channel==%d&&%s&&%s'%((i+1),stepch_tcut,cut)
        rdhist.Add(makeTH1(rfname, tname, 'data', binning, plotvar, rdtcut))

rdhist = rdhistlist[0]
for i, h in enumerate(rdhistlist):
    if i != 0 :
        rdhist.Add(rdhistlist[i])


rdhist.SetMinimum(minp)
rdhist.SetMaximum(maxp)

#Drawing plots on canvas
var = plotvar.split(',')[0]
var = ''.join(i for i in var if not i.isdigit())
var = var.replace('.','_').lower()
var = var.replace('()','')

outfile = "%s_s%d_%s.png"%(channel_name[channel-1],step,var)
if channel == 0: outfile = "Dilepton_s%d_%s.png"%(step,var)
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

print "Integral (MC) : ", step, channel, fIntMC
print "Integral (data) : ", step, channel, rdhist.Integral()
#print ("total entries for MC : ", mc_nev, " | total entries for RD : ", rd_nev)

canv.SaveAs(outfile)        
print outfile
