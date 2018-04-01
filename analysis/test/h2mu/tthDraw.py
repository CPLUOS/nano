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
h2muDraw.py -c 'll_m>50&&step>=5&&isTight==1&&filtered==1' -b [100,0,100] -p lep2_pt -x 'sub-leading muon p_{T}  [GeV]' &
h2muDraw.py -c 'll_m>50&&step>=5&&isTight==1&&filtered==1' -b [100,0,100] -p lep1_pt,lep2_pt -x 'muon p_{T}  [GeV]' &
h2muDraw.py -c 'll_m>50&&step>=5&&isTight==1&&filtered==1' -b [100,0,100] -p met -x 'met  [GeV]' &
h2muDraw.py -c 'll_m>50&&step>=5&&isTight==1&&filtered==1' -b [100,-3,3] -p lep1_eta,lep2_eta -x '#eta' &
'''

json_used = 'Golden'
datalumi = 36814 #35.9fb-1
#datalumi =  5975 #35.9fb-1

#datalumi = 8360.481454 #35.9fb-1
version = os.environ['CMSSW_VERSION']
user = os.environ['USER']


rootfileDir = "/xrootd/store/user/{}/nanoAOD/test/results_merged/tth2mu_".format(user)

#rootfileDir = "%s/src/nano/analysis/topMass/Results/results_merged/topmass_"% os.environ['CMSSW_BASE']
#rootfileDir = "/xrootd/store/user/pseudotop/ntuples/results_merged/v7-6-3/h2muAnalyzer_"
#rootfileDir = "%s/src/CATTools/CatAnalyzer/test/results_merged/h2muAnalyzer_" % os.environ['CMSSW_BASE']
#rootfileDir = "%s/cattuples/20160324_163101/results_merged/h2muAnalyzer_" % os.environ['HOME_SCRATCH']

CMS_lumi.lumi_sqrtS = "%.2f fb^{-1}, #sqrt{s} = 13 TeV 25ns "%(float(datalumi)/1000)
mcfilelist = [
#              'WMinusH_HToMuMu',
#              'WPlusH_HToMuMu',
#              'ZH_HToMuMu',
#              'VBF_HToMuMu',
#              'GG_HToMuMu',
               'ttH',
             # 'GluGluToZZTo2mu2tau',
             # 'GluGluToZZTo2e2mu',
             # 'GluGluToZZTo4mu',
             # 'ttZToLLNuNu',
             # 'ZZTo4L_powheg',
             # 'ZZTo2L2Nu_powheg',
     #          'ZZ',
             # 'WWTo2L2Nu',
     #          'WW',
#              'WZTo2LQQ',
             # 'WZTo3LNu_powheg',
     #          'WZ',
#              "WWTo2L2Nu",
#              "WZTo3LNu_amcatnlo",
             # "WZTo2L2Q",
#              "ZZTo2L2Nu",
#              "ZZTo2L2Q",
#              "ZZTo4L",
               "WWW",
               "WWZ",
               "WZZ",
               "ZZZ",
               "TTZToLLNuNu",
             # "ttWToLNu",
             # "SingleTop_tW_noHadron",
             # "SingleTbar_tW_noHadron",
#              "SingleTop_tW",
#              "SingleTbar_tW",
               "TTJets_aMC",
             # "TTJets_DiLept",
             # "TTJets_DiLept_Tune4",
             # 'TT_powheg',
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
             # 'DoubleMuon_Run2016',#mumu
              'SingleMuon_Run2016',#mumu
              #'SingleMuon_Run2015C',
             #'SingleMuon_Run2015D'
             ]

datasets = json.load(open("%s/src/nano/analysis/data/dataset/dataset.json" % os.environ['CMSSW_BASE']))
#cut_step = "(step>=5)"
#cut = 'lep1.Pt()>60&&lep2.Pt()>60&&dilep.M()>60&&step>=5'
#cut = 'dilep.M()>60&&step>4&&filtered&&MVA_BDT>-0.0246'
cut = 'Dilep.M()>60&&b_charge==1'
#cut = 'filtered==1&&%s&&%s'%(cut_step,emu_pid)
#cut = 'channel==2'
print cut
#weight = 'genweight*puweight*mueffweight*eleffweight*tri'
#weight = 'weight*(mueffweight)'
#weight = 'genweight*puweight'
weight = 'weight*mueffweight*btagweight'
#plotvar = 'met'
plotvar = 'Dilep.M()'
binning = [150, 50, 200]#[150, 50, 200]
#x_name = 'Invariant Mass[GeV]'
x_name = 'Invariant Mass'
y_name = 'Events'
dolog = True
f_name = 'Dilep_Test'
minp = 0.05
maxp = 1000000000
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
print plotvar, x_name, f_name
if "FL" in f_name:
   maxp = 100
   minp = 0.0005
if "FH" in f_name:
   maxp = 100000000
   minp = 0.0005
if "SL" in f_name:
   maxp = 10000000
   minp = 0.0005
   
tname = "events"

mchistList = []

dolog = True
tcut = '(%s)*%s'%(cut,weight)
#tcut = '(%s)'%(cut)
rdfname = rootfileDir + rdfilelist[0] +".root"

sig=[0,0,0,0,0,0]
bg=[0,0,0,0,0,0]
lumilist= [datalumi,300*1000,900*1000,3000*1000] 

#if plotvar == 'Dilep.M()':
#    f2_txt = open("significance_%s.txt"%(f_name),"w")
for imc,mcname in enumerate(mcfilelist):
    data = findDataSet(mcname, datasets)
   # if "DYJets_MG2" in mcname:
   #     scale = datalumi*6025.2
   # else:
    scale = datalumi*data["xsec"]
    print"scale: %s " %(scale)
    colour = data["colour"]
    title = data["title"]
    #if 'DYJets' in mcname: 
        #scale = scale*dyratio[channel][step] 
    #    scale = scale*dyratio[channel][1] 
#    if "HToMuMu" in mcname:
#        scale = scale*30.
#        title = title+" #times 30"
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
    #if "tt" in mcname: 
    #    scale = scale*1.07  
    mchist = makeTH1(rfname, tname, title, binning, plotvar, tcut, scale)    
    mchist.SetFillColor(colour)
    mchist.SetLineColor(colour)
    mchistList.append(mchist)

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
"""        
print "rdfname: %s\n tname: %s\n binning: %s\n plotvar: %s\n cut: %s\n"%(rdfname, tname, binning, plotvar, cut)
#rdhist = makeTH1(rdfname, tname, 'data', binning, plotvar, tcut+'&&(Dilep.M()<120||Dilep.M()>130)')
rdhist = makeTH1(rdfname, tname, 'data', binning, plotvar, cut+'&&(Dilep.M()<120||Dilep.M()>130)')
#rdhist.SetMinimum(minp)  
#rdhist.SetMaximum(maxp)
canv = drawTH1(f_name, CMS_lumi, mchistList, rdhist, x_name, y_name,dolog)
canv.SaveAs(f_name+".png")            
"""
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
