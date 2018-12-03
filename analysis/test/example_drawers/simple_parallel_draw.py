#!/usr/bin/env python


import ROOT, json, os, sys
from nano.analysis.parallel_draw_histos import parallel_draw_histos


listDatasets = [
  "SingleTop_t-channel", "SingleTbar_t-channel", 
  "SingleTop_s-channel", "SingleTop_tW", "SingleTbar_tW", 
  "TT_powheg", 
  "WToLNu_0J", "WToLNu_1J", "WToLNu_2J", 
  "DYJets_MG", "DYJets_MG_10to50", 
  "WW", "WZ", "ZZ", 
  "QCD_Pt-15to20_MuEnriched", "QCD_Pt-20to30_MuEnriched", "QCD_Pt-30to50_MuEnriched", 
  "QCD_Pt-50to80_MuEnriched", "QCD_Pt-80to120_MuEnriched", "QCD_Pt-120to170_MuEnriched", 
  "QCD_Pt-170to300_MuEnriched", "QCD_Pt-300to470_MuEnriched", "QCD_Pt-470to600_MuEnriched", 
  "QCD_Pt-600to800_MuEnriched", "QCD_Pt-800to1000_MuEnriched", "QCD_Pt-1000toInf_MuEnriched", 
  "SingleMuon_Run2016B", "SingleMuon_Run2016C", "SingleMuon_Run2016D", "SingleMuon_Run2016E", 
  "SingleMuon_Run2016F", "SingleMuon_Run2016G", "SingleMuon_Run2016H"
]

strPathDraw = "%s/src/nano/analysis/test/example_drawers"%os.environ[ "CMSSW_BASE" ]
strPathSrc = "%s/src/nano/analysis/test/example_drawers"%os.environ[ "CMSSW_BASE" ]

strNameSrc = "testlist"
bMulticore = False # Change this true if you want to run this on condor cluster

strSrcFormat  = "root://cms-xrdr.sdfarm.kr:1094///xrd"
strSrcFormat += "/store/user/quark2930/singletop/%(srcname)s/dir_%(dataset)s/res_%(idx)s.root"

strCut = "step >= 4 && trig_m > 0 && lep_pid == 13"

# Launching and running the parallel_draw_histos module
paraMain = parallel_draw_histos()

paraMain.SetSrcPath(os.path.join(strPathSrc, strNameSrc + ".txt"))
paraMain.SetSrcPathFormat(strSrcFormat)
paraMain.SetDirHistName()
paraMain.SetNameTree("event")
paraMain.SetCutMC(strCut)
paraMain.SetCutRD(strCut)
paraMain.SetWeightMC("genweight * puweight * tri * mueffweight * btagweight")
paraMain.SetWeightRD("1")
paraMain.SetDatasets(listDatasets)
paraMain.SetVars(
  {
    "lep.Pt()": {"bin": [20, 10, 210], "xaxis": "lep p_{T}", "yaxis": "Events"}, 
    "abs(lep.Eta())": {"bin": [20, 0, 2.5], "xaxis": "lep |#eta|", "yaxis": "Events"}, 
    "lep.Pt(): jet1.Eta()": {"bin": [20, -4.7, 4.7, 20, 0, 60], "xaxis": "jet' #eta", "yaxis": "lep p_{T}"}, 
  })
paraMain.SetPathDraw(strPathDraw)

paraMain.InitRun()

if bMulticore: 
  paraMain.RunOnCluster()
else: 
  paraMain.RunOnWorknode()


