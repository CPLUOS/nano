#!/usr/bin/env python
import ROOT, nano.analysis.CMS_lumi, json, os, getopt, sys
from nano.analysis.histoHelper import *
from ROOT import TLorentzVector
#import DYestimation
ROOT.gROOT.SetBatch(True)
'''
h2muDraw.py -c 'll_m>50&&step>=5&&isTight==1&&filtered==1' -b [40,0,40] -p nvertex -x 'no. vertex' &
h2muDraw.py -c 'll_m>50&&step>=5&&isTight==1&&filtered==1' -b [200,0,200] -p ll_m -x 'mass [GeV]' &
h2muDraw.py -c 'll_m>50&&step>=5&&isTight==1&&filtered==1' -b [100,0,100] -p ll_pt -x 'diMuon p_{T}  [GeV]' &
h2muDraw.py -c 'll_m>50&&step>=5&&isTight==1&&filtered==1' -b [100,0,100] -p lep1_pt -x 'leading muon p_{T}  [GeV]' &
h2muDraw.py -c 'll_m>50&&step>=5&&isTight==1&&filtered==1' -b [100,0,100] -p met -x 'met  [GeV]' &
h2muDraw.py -c 'll_m>50&&step>=5&&isTight==1&&filtered==1' -b [100,-3,3] -p lep1_eta,lep2_eta -x '#eta' &
'''

json_used = 'Golden'
#datalumi = 36814 #35.9fb-1
datalumi =  35900 #35.9fb-1
#datalumi =  5746.010293 # B 5.74fb-1
#datalumi =  2572.903489 # C 5.74fb-1
#datalumi =  5222.873966 # G 5.74fb-1
#datalumi =  8360.481454 # H 5.74fb-1

#datalumi = 8360.481454 #35.9fb-1
version = os.environ['CMSSW_VERSION']
user = os.environ['USER']


rootfileDir = "/xrootd/store/user/{}/nanoAOD/v5_Pu/results_merged/tth2mu_".format(user)
#rootfileDir = "%s/src/nano/analysis/topMass/Results/results_merged/topmass_"% os.environ['CMSSW_BASE']
#rootfileDir = "/xrootd/store/user/pseudotop/ntuples/results_merged/v7-6-3/h2muAnalyzer_"
#rootfileDir = "%s/src/CATTools/CatAnalyzer/test/results_merged/h2muAnalyzer_" % os.environ['CMSSW_BASE']
#rootfileDir = "%s/cattuples/20160324_163101/results_merged/h2muAnalyzer_" % os.environ['HOME_SCRATCH']

CMS_lumi.lumi_sqrtS = "%.2f fb^{-1}, #sqrt{s} = 13 TeV 25ns "%(float(datalumi)/1000)
mcfilelist = [
             # 'ttH_HToMuMu',
             # 'WMinusH_HToMuMu',
             # 'WPlusH_HToMuMu',
             # 'ZH_HToMuMu',
             # 'VBF_HToMuMu',
             # 'GG_HToMuMu',
             # 'GluGlu_ZZTo4e',
             # 'GluGlu_ZZTo2e2mu',
             # 'GluGlu_ZZTo4mu',
             # 'ZZTo4L_powheg',
             # 'ZZTo2L2Nu_powheg',
             # 'TTZToLLNuNu',
             # "WWW",
             # "WWZ",
             # "WZZ",
             # "ZZZ",
             # 'ZZ',
             # 'WWTo2L2Nu',
               'WW',
             # 'WZTo2LQQ',
             # 'WZTo3LNu_powheg',
             # 'WZ',
             # "WWTo2L2Nu",
             # "WZTo3LNu_amcatnlo",
             # "WZTo2L2Q",
             # "ZZTo2L2Nu",
               "ZZTo2L2Q",
               "ZZTo4L_powheg",
               'WJets',
             # "tZq_ToLL",
               "ttWToLNu",
             # "SingleTop_tW_noHadron",
             # "SingleTbar_tW_noHadron",
               "TTJets_aMC",
               'TTWJetsToLNu', 
             #  "SingleTop_tW",
             #  "SingleTbar_tW",
               "TTJets_DiLept",
             # "TTJets_DiLept_Tune4",
               'TT_powheg',
             #  "DYToLL_2J",
               'DYJets',
             # 'DYJets_MG_10to50',
             # 'DYJets_MG2',
             # 'DYJets_2J',
             # 'DYJets_1J',
             # 'DYJets_0J',
             # 'DYJets_10to50',
             ]#ref : https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToMuMu
#mcfilelist = ['VBF_HToMuMu','WW','WZ','ZZ','TT_powheg','DYJets','DYJets_10to50']#,'WJets']
#mcfilelist = [ 'TTJets_aMC']
rdfilelist = [
             # 'SingleMuon_Run2016',
               'DoubleMuon_Run2016',#mumu
             # 'SingleMuon_Run2016C',#mumu
             # 'SingleMuon_Run2015C',
             #'SingleMuon_Run2015D'
             ]

datasets = json.load(open("%s/src/nano/nanoAOD/data/dataset/dataset.json" % os.environ['CMSSW_BASE']))
#cut_step = "(step>=5)"
#cut = 'Dilep.M()>=120&&Dilep.M()<=130&&nonB==1'
cut = 'Dilep.M()>=60'
#cut = 'filtered==1&&%s&&%s'%(cut_step,emu_pid)
#cut = 'channel==2'
print cut
#weight = 'genweight*puweight*mueffweight*eleffweight*tri'
#weight = 'weight*(mueffweight)'
#weight = 'genweight*puweight'
weight = 'genweight*puweight*mueffweight*btagweight'
#weight = 'weight'
#plotvar = 'njet'
plotvar = 'Dilep.M()'
binning = [50, -0.3, 0.3]#[150, 50, 200]
#binning = [50, 0, 1]#[150, 50, 200]
#x_name = 'Invariant Mass[GeV]'
x_name = 'mass'
y_name = 'Events'
dolog = True
f_name = 'Dlep'
try:
    opts, args = getopt.getopt(sys.argv[1:],"hdc:w:b:p:x:y:f:j:",["cut","weight","binning","plotvar","x_name","y_name","f_name","json_used","dolog"])
except getopt.GetoptError:          
    print 'Usage : ./cmesonDraw.py -c <cut> -w <weight> -b <binning> -p <plotvar> -x <x_name> -y <y_name> -f <f_name> -j <json_used> -d <dolog>'
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print 'Usage : ./cmesonDraw.py -c <cut> -w <weight> -b <binning> -p <plotvar> -x <x_name> -y <y_name> -f <f_name> -j <json_used> -d <dolog>'
        sys.exit()
    elif opt in ("-c", "--cut"):
        cut = arg
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
    elif opt in ("-f", "--f_name"):
        f_name = arg
    elif opt in ("-j", "--json_used"):
        json_used = arg
    elif opt in ("-d", "--dolog"):
        dolog = True
print plotvar,      x_name,      f_name

for a in ("B", "C", "D", "E", "F", "G", "H"):
   for b in rdfilelist:   
      if b.endswith(a):
         f_name = a+"_"+f_name
         print a
         print f_name
         break 
   
tname = "events"

mchistList = []

dolog = True
tcut = '(%s)*%s'%(cut,weight)
#tcut = '(%s)'%(cut)
rdfname = rootfileDir + rdfilelist[0] +".root"
"""
sig=[0,0,0,0,0,0]
bg=[0,0,0,0,0,0]
lumilist= [datalumi,300*1000,900*1000,3000*1000] 
"""

sysNameList = ["mueffweight", "puweight"] 
sysErr_up = []
sysErr_dn = []
staNameList = ["Lumi"]
staErr = []

for sysname in sysNameList:
   sysErr_up.append(defTH1(sysname+'_up', sysname+'_up', binning))
   sysErr_dn.append(defTH1(sysname+'_dn', sysname+'_dn', binning))
for staname in staNameList:
   staErr.append(defTH1(staname, staname, binning))

sig = 0
bg = 0
TTZ = 0
TTW = 0 
TTJets = 0
SingleTop = 0
SingleTbar = 0 
TTH = 0
cout = 0 
#if plotvar == 'Dilep.M()':
#    f2_txt = open("significance_%s.txt"%(f_name),"w")
for imc,mcname in enumerate(mcfilelist):
    data = findDataSet(mcname, datasets)
    scale = datalumi*data["xsec"]
    print"scale: %s " %(scale)
    colour = data["colour"]
    title = data["title"]
    #if 'DYJets' in mcname: 
        #scale = scale*dyratio[channel][step] 
    #    scale = scale*dyratio[channel][1] 
    rfname = rootfileDir + mcname +".root"
    print rfname
    tfile = ROOT.TFile(rfname)
    #wentries = tfile.Get("genweight").Integral()
    wentries = tfile.Get("weight").Integral()
    print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Wentires: %s"%(wentries)
    #wentries = tfile.Get("Events").Integral(0,1)
    #print wentries
    scale = scale/wentries
    print"new scale: %s" %(scale)
    mchist = makeTH1(rfname, tname, title, binning, plotvar, tcut, scale)    
   
    for i, sysname in enumerate(sysNameList):
      if "weight" in sysname:
         sysErr_up[i].Add(makeTH1(rfname, tname, title, binning, plotvar, tcut.replace(sysname,sysname+'_up'), scale))
         sysErr_dn[i].Add(makeTH1(rfname, tname, title, binning, plotvar, tcut.replace(sysname,sysname+'_dn'), scale))
    for i, staname in enumerate(staNameList):
         newscale = scale + (((datalumi*0.025)*data["xsec"])/wentries)
         staErr[i].Add(makeTH1(rfname, tname, title, binning, plotvar, tcut, newscale))

    mchist.SetFillColor(colour)
    mchist.SetLineColor(colour)

    remchist = makeTH1(rfname, tname, title, binning, plotvar, tcut, scale)
    if any(a in mcname for a in ("ttH", "HToMuMu")):
      sig += remchist.Integral(remchist.FindBin(120), remchist.FindBin(130))
    else: 
      bg += remchist.Integral(remchist.FindBin(120), remchist.FindBin(130))
      cout += 1
      print mcname
#    if "ttH" in mcname:  
#      TTH = remchist.Integral(remchist.FindBin(120), remchist.FindBin(130))
    if 'TTWJetsToLNu' in mcname:  
      TTW = remchist.Integral(remchist.FindBin(120), remchist.FindBin(130))
    if 'TTZToLLNuNu' in mcname:
      TTZ = remchist.Integral(remchist.FindBin(120), remchist.FindBin(130))
    if "TTJets" in mcname:
      TTJets = remchist.Integral(remchist.FindBin(120), remchist.FindBin(130))
    if "SingleTop_tW" in mcname:
      SingleTop = remchist.Integral(remchist.FindBin(120), remchist.FindBin(130))
    if "SingleTbar_tW" in mcname:
      SingleTbar = remchist.Integral(remchist.FindBin(120), remchist.FindBin(130))

    mchistList.append(mchist)

if "Dilep.M()" in plotvar:
   print "TTH          = {}".format(TTH)          
   print "TTW          = {}".format(TTW)          
   print "TTZ          = {}".format(TTZ)          
   print "TTJet        = {}".format(TTJets)          
   print "SingleTop    = {}".format(SingleTop)          
   print "SingleTbar   = {}".format(SingleTbar)          
   print "signal       = {}".format(sig) 
   print "background   = {}".format(bg)
   print "significance = {}".format(sig/math.sqrt(sig+bg))
   print cout

if plotvar == "Dilep.M()": cut = cut+'&&(Dilep.M()<120||Dilep.M()>130)'            
print "rdfname: %s\n tname: %s\n binning: %s\n plotvar: %s\n cut: %s\n"%(rdfname, tname, binning, plotvar, cut)
rdhist = makeTH1(rdfname, tname, 'data', binning, plotvar, cut)
nbins = rdhist.GetNbinsX()
if dolog == True: 
   minp = 0.005
   maxp = 1000000000
   if "FL" in f_name:
      maxp = 100
      minp = 0.0005
   if "FH" in f_name:
      maxp = 1000000000
      minp = 0.005
   if "XL" in f_name:
      maxp = 1000000
      minp = 0.005
   rdhist.SetMinimum(minp)  
   rdhist.SetMaximum(maxp)

### Error Band ###
errorBand = copy.deepcopy(rdhist)
errorBand.SetFillColor(14)
errorBand.SetFillStyle(3001)
errorBand.SetMarkerStyle(0)

h_nom = defTH1("nom","nom",binning)
for h in mchistList:
   h_nom.Add(h)

for i in  range(len(sysNameList)):
   sysErr_up[i].Add(h_nom, -1)
   sysErr_dn[i].Add(h_nom, -1)
   for j in range(1,nbins+1):
      maxErr = max(abs(sysErr_up[i].GetBinContent(j)), abs(sysErr_dn[i].GetBinContent(j)))
      sumErr = math.sqrt(errorBand.GetBinError(j)**2+maxErr**2)
      errorBand.SetBinError(j, sumErr)
for i in range(len(staNameList)):
   staErr[i].Add(h_nom, -1)
   for j in range (1,nbins+1):
      sumErr = math.sqrt(errorBand.GetBinError(j)**2+staErr[i].GetBinContent(j)**2)
      errorBand.SetBinError(j, sumErr)

canv = drawTH1(f_name, CMS_lumi, mchistList, rdhist, x_name, y_name,dolog)
mainPad = canv.GetPrimitive("mainPad")
mainPad.cd()
errorBand.Draw("e2same")
mainPad.GetPrimitive("data").Draw("esamex0")
canv.Update()
ratioPad = canv.GetPrimitive("ratioPad")
ratioPad.cd()
sysErrRatio = errorBand.Clone()
sysErrRatio.Divide(rdhist)
sysErrRatio.Draw("e2same")
ratioPad.GetPrimitive("hratio").Draw("esame")
canv.Update()
canv.SaveAs(f_name+".png")      

"""
    for l,lumi in enumerate(lumilist):
        rescale = lumi*data["xsec"]/wentries
        remchist = makeTH1(rfname, tname, title, binning, plotvar, tcut, rescale)    
        if "HToMuMu" in mcname:
            sig[l]+= remchist.Integral(remchist.FindBin(120),remchist.FindBin(130))
        else:
            bg[l]+= remchist.Integral(remchist.FindBin(120),remchist.FindBin(130))
    lumi3 = 0
    while ((sig[4]/math.sqrt(sig[4]+bg[4]))<3):
        rescale = lumi3*data["xsec"]/wentries
        remchist = makeTH1(rfname, tname, title, binning, plotvar, tcut, rescale)    
        if "HToMuMu" in mcname:
            sig[4]+= remchist.Integral(remchist.FindBin(120),remchist.FindBin(130))
        else:
            bg[4]+= remchist.Integral(remchist.FindBin(120),remchist.FindBin(130))
        if (imc < len(mcfilelist)-1):break
        lumi3 += 10
    lumilist.append(lumi3)
    lumi5 = 0
    while ((sig[j]/math.sqrt(sig[j]+bg[j]))<5):
        lumi5 += 10
        rescale = lumi5*data["xsec"]/wentries
        remchist = makeTH1(rfname, tname, title, binning, plotvar, tcut, rescale)    
        if "HToMuMu" in mcname:
            sig[5]+= remchist.Integral(remchist.FindBin(120),remchist.FindBin(130))
        else:
            bg[5]+= remchist.Integral(remchist.FindBin(120),remchist.FindBin(130))
 
print "="*50
print rfname
print "="*50
x_min = 110
f_txt_bw = open("bw_%s.txt"%(f_name),"w")
while (x_min<140):
  if plotvar == 'dilep.M()':# blind data around higgs mass
    f_txt = open("events_%s.txt"%(f_name),"w")
    print>>f_txt, "Run data : %s\n cut : %s\n # : \n %d\n"%(rdfilelist[0],f_name,rdhist.Integral(rdhist.FindBin(100),rdhist.FindBin(110)))
    value=[0,0,0,0]
    value[0],value[1],value[2],value[3] = drawBWFit("bw_"+f_name+".png",rdhist,88,94)
    print>>f_txt_bw, "==== %d ===="%(x_min)   
    print>>f_txt_bw, "f_name : %s\n mean : %f\n mean error : %f\n gamma : %f\n gamma error : %f\n"%(f_name,value[0],value[1],value[2],value[3])
    #f_txt2 = open("eventlist_%s_%s.txt"%(rdfilelist[0],f_name),"w")
    #print>>f_txt2, 
 #   print>>f2_txt, " cut : %s\n"%(f_name)
 #   for j in range(6):
 #       print>>f2_txt, "*"*(50)
 #       print>>f2_txt, " datalumi : %s\n sig : %s\n bg : %s\n significance : %s\n"%(lumilist[j],sig[j],bg[j],(sig[j]/math.sqrt(sig[j]+bg[j])))
   
    parameterization("fit_"+f_name+"_%d.png"%(x_min),CMS_lumi,rdhist, mchistList, x_min, binning[2], value[0], value[1], value[2], value[3])
    if 'SingleMuon' in rdfname:
        if len(binning) == 3:
            htmp = ROOT.TH1D("tmp", "tmp", binning[0], binning[1], binning[2])
        else:
            htmp = ROOT.TH1D("tmp", "tmp", len(binning)-1, array.array('f', binning))
        for i in range(binning[1],binning[2]):
            if (rdhist.FindBin(120)<=i<=rdhist.FindBin(130)):continue
            entries=rdhist.GetBinContent(i)
            htmp.SetBinContent(i,entries)
            
    #after blind the signal region.
    parameterization("fit_"+f_name+"_%d_signal_region_blinded.png"%(x_min), CMS_lumi,htmp, mchistList, x_min, binning[2], value[0], value[1], value[2], value[3],False,True)
    #f_txt2.close()
    f_txt.close()
    f2_txt.close()
    x_min +=10
    if x_min>110:break

f_txt_bw.close()
"""
