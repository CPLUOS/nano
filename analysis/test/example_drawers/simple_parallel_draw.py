#!/usr/bin/env python


################################################################################
# 
# Author: Byeonghak Ko (quark2930@hotmail.com; bko@cern.ch)
# 
# simple_parallel_draw.py
# 
# This code is a simple demonstration of parallel_draw_histos class
# To run this, just run the following: 
#   ./simple_parallel_draw.py
# After running, a directory might have been created. In the directory, 
#  the root files hist_*.root contain the histograms.
# Users can stack them and draw together, with python/histoHelper.py.
# 
# If you want a more complicated one, please read 
#   nano/analysis/test/singletop/draw/parallel_draw_histos.py
# 
################################################################################


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

# As explained on SetSrcPathFormat() in python/parallel_draw_histos.py, 
# this is for extraction of the dataset name and the idx of a given ntuple file path.
# For example, if the path is 
#  root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/quark2930/singletop/testlist/dir_WW/res_7.root, 
# then "dataset" = "WW" and "idx" = "7", while "srcname" = "testlist" which is not used.
strSrcFormat  = "root://cms-xrdr.sdfarm.kr:1094///xrd"
strSrcFormat += "/store/user/quark2930/singletop/%(srcname)s/dir_%(dataset)s/res_%(idx)s.root"

strCut = "step >= 4 && trig_m > 0 && lep_pid == 13"

# Launching and running the parallel_draw_histos module
paraMain = parallel_draw_histos()

paraMain.SetSrcPath(os.path.join(strPathSrc, strNameSrc + ".txt"))
paraMain.SetSrcPathFormat(strSrcFormat)
paraMain.SetDatasets(listDatasets)

paraMain.SetNameTree("event")
paraMain.SetCutMC(strCut)
paraMain.SetCutRD(strCut)
paraMain.SetWeightMC("genweight * puweight * tri * mueffweight * btagweight")
paraMain.SetWeightRD("1")

# See the description on GetVars() in python/parallel_draw_histos.py for more informations
paraMain.SetVars(
  {
    "lep.Pt()": {"bin": [20, 10, 210], "xaxis": "lep p_{T}", "yaxis": "Events"}, 
    "abs(lep.Eta())": {"bin": [20, 0, 2.5], "xaxis": "lep |#eta|", "yaxis": "Events"}, 
    "lep.Pt(): jet1.Eta()": {"bin": [20, -4.7, 4.7, 20, 0, 60], "xaxis": "jet' #eta", "yaxis": "lep p_{T}"}, 
  })

paraMain.SetPathDraw(strPathDraw)
paraMain.SetDirHistName()

paraMain.InitRun()

if bMulticore: 
  paraMain.RunOnCluster()
else: 
  paraMain.RunOnWorknode()


