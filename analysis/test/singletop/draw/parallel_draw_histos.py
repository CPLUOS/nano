#!/usr/bin/env python


import ROOT, json, os, sys, copy
import getopt

from nano.analysis.parallel_draw_histos import parallel_draw_histos


strPathDraw = "%s/src/nano/analysis/test/singletop/draw"%os.environ[ "CMSSW_BASE" ]
strPathSrc = "%s/src/nano/analysis/test/singletop/rootlist/"%os.environ[ "CMSSW_BASE" ]


if len(sys.argv) <= 1: 
  print "Requiring arguments: the directory the ntuple is saved", 
  print "(not full of the path, just the end of it)"
  sys.exit(1)


strSpellNoDup  = "if [ `ps -ef | grep $UID | grep pytho[n] | grep draw_histos.py | wc -l` -gt 1 ]; "
strSpellNoDup += "then ./cannot_exists_one 2> /dev/null ; else ls > /dev/null ; fi"

if os.system(strSpellNoDup) != 0: 
  print "You are already running this code."
  sys.exit(0)

#strPathSrc = "/xrootd/store/user/quark2930/singletop/%s/"%sys.argv[ 1 ]
#strPathTMVA = "/xrootd/store/user/quark2930/singletop/tmva_%s/"
strPathListSet = os.path.join(strPathDraw, "listSet.json")

strNameSrc = sys.argv[ 1 ]

strNameTree = "event"
channel = 2 # 0: all, 1: el, 2: mu, 3: el + mu, (-): anti
bAllcharge = False
bNoChargeCut = False
strWeighthead = "gptmb"
nStep = 4

cut = ""

strCutAdd = ""

strNameTMVA = ""

bMulticore = True
bDefaultListSet = True

strHelpHowto = "Usage: ./[name of this py file] [dir name of roots] " + \
  "-a <channel> -c <cut> -w <weight> -n <local running> -l <listSet JSON file>"
try:
  opts, args = getopt.getopt(sys.argv[ 2: ], "hneAa:c:w:s:t:l:", 
    ["channel", "cut", "allcharge", "weight", "step", "tmva", "listset"])
except getopt.GetoptError:          
  sys.stderr.write(strHelpHowto + "\n")
  sys.exit(2)

for opt, arg in opts:
  if opt == '-h':
    print strHelpHowto
    sys.exit()
  elif opt in ("-a", "--channel"):
    channel = int(arg)
  elif opt in ("-A"):
    bNoChargeCut = True
  elif opt in ("-c", "--cut"):
    strCutAdd = arg
  elif opt in ("-e", "--allcharge"): 
    bAllcharge = True
  elif opt in ("-w", "--weight"):
    strWeighthead = arg
  elif opt in ("-s", "--step"):
    nStep = int(arg)
  elif opt in ("-t", "--tmva"):
    strNameTMVA = arg
  elif opt in ("-n"):
    bMulticore = False
  elif opt in ("-l", "--listset"):
    bDefaultListSet = False
    if not arg.startswith(strPathDraw): strPathListSet = os.path.join(strPathDraw, arg)
    else:                               strPathListSet = arg

if nStep > 0: cut += " && step >= %i"%nStep

if channel != 0: 
  strSign = "" if channel > 0 else "-"
  if abs(channel) == 1: 
    cut += " && trig_e > 0"
    if not bNoChargeCut: cut += " && lep_pid == %s11"%strSign if not bAllcharge else " && abs(lep_pid) == 11"
    
    strNameSrc += "_el"
    if strNameTMVA != "": strNameTMVA += "_el"
  if abs(channel) == 2: 
    cut += " && trig_m > 0"
    if not bNoChargeCut: cut += " && lep_pid == %s13"%strSign if not bAllcharge else " && abs(lep_pid) == 13"
    
    strNameSrc += "_mu"
    if strNameTMVA != "": strNameTMVA += "_mu"
  if abs(channel) == 3: 
    cut += " && ( lep_pid == %s11 || lep_pid == %s13)"%(strSign, strSign)

weight = ""
if "g" in strWeighthead: weight += " * genweight"
if "p" in strWeighthead: weight += " * puweight"
if "t" in strWeighthead: weight += " * tri"
if "m" in strWeighthead: weight += " * %seffweight"%("el" if abs(channel) == 1 else "mu")
if "b" in strWeighthead: weight += " * btagweight"
weight = weight[ 3: ] if weight != "" else "1.0"

cut += " && ( " + strCutAdd + " )" if len(strCutAdd) > 0 else ""
cut = cut[ 4: ]

#if not os.path.exists(strPathSrc): 
#  sys.stderr.write("Error: cannot find the directory of nano ntuples\n")
#  sys.exit(1)

dicSet = json.load(open(strPathListSet))

# Loading name of datasets
listSets = []
if channel == 0 or ( abs(channel) & 1 ) != 0: listSets += dicSet[ "listDatasets" ][ "RDEL" ]
if channel == 0 or ( abs(channel) & 2 ) != 0: listSets += dicSet[ "listDatasets" ][ "RDMU" ]

listSets += dicSet[ "listDatasets" ][ "SIG" ]
listSets += dicSet[ "listDatasets" ][ "BKG" ]
if channel == 0 or abs(channel) == 1: listSets += dicSet[ "listDatasets" ][ "BKG_EL" ]
if channel == 0 or abs(channel) == 2: listSets += dicSet[ "listDatasets" ][ "BKG_MU" ]

# Setting variables which will be drawn
# Before giving this, the tmva ones are needed to be removed if no TMVA ones
dicVars = dicSet[ "Vars" ]

if strNameTMVA == "": 
  dicVars = {strKey: dicVars[ strKey ] for strKey in dicVars.keys() if "tmva" not in strKey}

strPattern  = "root://cms-xrdr.sdfarm.kr:1094///xrd"
strPattern += "/store/user/quark2930/singletop/%(srcname)s/dir_%(dataset)s/res_%(idx)s.root"

# Launching and running the parallel_draw_histos module
paraMain = parallel_draw_histos()

paraMain.SetDirHistName()
#paraMain.SetSrcPath(os.path.join(strPathSrc, "../../example_drawers/testlist.txt"))
paraMain.SetSrcPath(os.path.join(strPathSrc, strNameSrc + ".txt"))
paraMain.SetSrcPathFormat(strPattern)
paraMain.SetNameTree(strNameTree)
paraMain.SetCutMC(cut)
paraMain.SetCutRD(cut)
paraMain.SetWeightMC(weight)
paraMain.SetWeightRD("1")
paraMain.SetDatasets(listSets)
paraMain.SetVars(dicVars)
paraMain.SetPathDraw(strPathDraw)

if strNameTMVA != "": 
  paraMain.SetFriendPath(os.path.join(strPathSrc, "tmva_" + strNameTMVA + ".txt"))
  paraMain.SetNameFriendTree("tmva")

paraMain.InitRun()

if bMulticore: 
  paraMain.RunOnCluster()
else: 
  paraMain.RunOnWorknode()


