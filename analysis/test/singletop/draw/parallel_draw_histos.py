#!/usr/bin/env python


import ROOT, json, os, sys, copy


def valToName(strVar): 
  return strVar.replace(".", "").replace("(", "").replace(")", "").replace(" ", "").replace("*", "")


def makeHisto(strName, strTitle, listBin):
  if listBin[ 0 ] >= 0: 
    if len(listBin) == 3: 
      return ROOT.TH1D(strName, strTitle, listBin[ 0 ], listBin[ 1 ], listBin[ 2 ])
    elif len(listBin) == 6: 
      return ROOT.TH2F(strName, strTitle, listBin[ 0 ], listBin[ 1 ], listBin[ 2 ], 
        listBin[ 3 ], listBin[ 4 ], listBin[ 5 ])
    elif len(listBin) == 9: 
      return ROOT.TH3F(strName, strTitle, listBin[ 0 ], listBin[ 1 ], listBin[ 2 ], 
        listBin[ 3 ], listBin[ 4 ], listBin[ 5 ], listBin[ 6 ], listBin[ 7 ], listBin[ 8 ])
  else: 
    if listBin[ 1 ] == "log1D": 
      import array
      
      nBin = listBin[ 2 ]
      fBinRatio = ( listBin[ 4 ] / listBin[ 3 ] ) ** ( 1.0 / nBin )
      
      listBinAct = []
      fBin = listBin[ 3 ]
      
      for i in range(nBin + 1): 
        listBinAct.append(fBin)
        fBin *= fBinRatio
      
      return ROOT.TH1D(strName, strTitle, nBin, array.array("d", listBinAct)) 

strPathDraw = "%s/src/nano/analysis/test/singletop/draw"%os.environ[ "CMSSW_BASE" ]
strPathThis = os.path.join(strPathDraw, "parallel_draw_histos.py")
strPathListSet = os.path.join(strPathDraw, "listSet.json")

if len(sys.argv) <= 1: 
  print "Requiring arguments: the directory the ntuple is saved", 
  print "(not full of the path, just the end of it)"
  sys.exit(1)

strTypeArgOneRoot = "--oneroot"
strTypeArgMerger = "--merger"

if sys.argv[ 1 ] == strTypeArgOneRoot: 
  strJSONQuery = " ".join(s for s in sys.argv[ 2: ])
  dicMain = json.loads(strJSONQuery.replace("(+(", "\"").replace("'", ""))
  
  strWeight = dicMain[ "weight" ].encode("ascii", "ignore")
  strCut    = dicMain[ "cut" ].encode("ascii", "ignore")
  
  strPathMain = dicMain[ "rootfile" ][ 0 ].encode("ascii", "ignore")
  strDirName = dicMain[ "rootfile" ][ 1 ].encode("ascii", "ignore")
  strNoFile = dicMain[ "rootfile" ][ 2 ].encode("ascii", "ignore")
  
  if "listset" in dicMain: strPathListSet = dicMain[ "listset" ].encode("ascii", "ignore")
  dicSet = json.load(open(strPathListSet))
  dicVar = dicSet[ "Vars" ]
  dicOutInfo = {}
  
  dicOutInfo[ "dataset" ] = dicMain[ "dataset" ].encode("ascii", "ignore")
  dicOutInfo[ "ntuple" ] = os.path.join(strPathMain, strDirName, strNoFile)
  dicOutInfo[ "no_file" ] = strNoFile.split("_")[ 1 ].split(".")[ 0 ]
  
  dicOutInfo[ "weight" ] = strWeight
  dicOutInfo[ "cut" ] = strCut
  
  strNameTree = dicMain[ "treename" ].encode("ascii", "ignore")
  fMain = ROOT.TFile.Open(dicOutInfo[ "ntuple" ])
  tree = fMain.Get(strNameTree)
  
  dicHist = {}
  strOutRoot = os.path.join(strPathDraw, dicMain[ "res_dir" ], 
    dicOutInfo[ "dataset" ], dicOutInfo[ "no_file" ] + ".root")
  
  # Getting # of whole entries
  dicOutInfo[ "entries" ] = "%lf"%( fMain.Get("nevents").Integral() )
  dicOutInfo[ "step0" ] = "%lf"%( fMain.Get("cutFlow").GetBinContent(1) )
  dicOutInfo[ "step1" ] = "%lf"%( fMain.Get("cutFlow").GetBinContent(2) )
  
  # Applying the cut
  # This way applies the cut only once for several histograms, so enhances the performances
  strNameAfterCut = strNameTree + "_aftercut"
  tree.Draw(">>" + strNameAfterCut, strCut, "entrylist")
  tree.SetEntryList(ROOT.gDirectory.Get(strNameAfterCut))
  
  # Drawing histograms
  for strVar in dicVar.keys(): 
    strNameVar = valToName(strVar)
    strNameHist = dicMain[ "histname" ].encode("ascii", "ignore") + "_" + strNameVar
    
    hTmp = makeHisto(strNameHist, strNameHist, dicVar[ strVar ][ "bin" ])
    tree.Project(strNameHist, strVar, strWeight)
    
    dicHist[ strVar ] = copy.deepcopy(hTmp)
    
    if dicHist[ strVar ].GetSumw2N() == 0:
      dicHist[ strVar ].Sumw2()
  
  #print "%s - %s : Drawing is done"%(datetime.datetime.now().strftime("%H:%M:%S"), dicMain[ "histname" ])
  
  fRootFile = ROOT.TFile(strOutRoot, "NEW")
  
  # Writing the string infos
  for strKey in dicOutInfo.keys(): ROOT.TNamed(strKey, dicOutInfo[ strKey ]).Write()
  for strVar in dicHist.keys(): dicHist[ strVar ].Write()
  
  fRootFile.Write()
  fRootFile.Close()
  
  os.system("touch %s"%(os.path.join(strPathDraw, dicMain[ "res_dir" ], 
    dicMain[ "dataset" ] + strNoFile.split("_")[ 1 ].split(".")[ 0 ] + ".done")))
elif sys.argv[ 1 ] == strTypeArgMerger:
  import shutil, time
  
  if len(sys.argv) < 5: 
    print "Wrong arguments"
    sys.exit(1)
  
  strDirHist = sys.argv[ 2 ]
  strDataset = sys.argv[ 3 ]
  nNumJobs = int(sys.argv[ 4 ])
  
  if nNumJobs <= 0: 
    print "Wrong arguments"
    sys.exit(1)
  
  strDirDataset = os.path.join(strPathDraw, strDirHist, strDataset)
  listRoot = []
  
  while True:
    listRoot = [ s for s in os.listdir(strDirDataset) if ".root" in s ]
    
    nDone = len(listRoot)
    if nDone >= nNumJobs: break
    time.sleep(1)
  
  fOrg = ROOT.TFile(os.path.join(strDirDataset, listRoot[ 0 ]))
  listRoot.pop(0)
  
  fEntries = float(fOrg.Get("entries").GetTitle())
  fStep0 = float(fOrg.Get("step0").GetTitle())
  fStep1 = float(fOrg.Get("step1").GetTitle())
  
  dicHisto = {}
  dicInfo  = {}
  i = 0
  
  while True:
    objkeyGet = fOrg.GetListOfKeys().At(i)
    if objkeyGet == None: break
    strNameObj = objkeyGet.GetName()
    
    if strNameObj.startswith("hist_"): 
      objGet = fOrg.Get(strNameObj)
      dicHisto[ strNameObj ] = copy.deepcopy(objGet)
      dicHisto[ strNameObj ].Sumw2(False)
    else: 
      dicInfo[ strNameObj ] = fOrg.Get(strNameObj).GetTitle()
    
    i += 1
  
  fOrg.Close()
  
  for strRootFile in listRoot: 
    fRead = ROOT.TFile(os.path.join(strDirDataset, strRootFile))
    
    fEntries += float(fRead.Get("entries").GetTitle())
    fStep0   += float(fRead.Get("step0").GetTitle())
    fStep1   += float(fRead.Get("step1").GetTitle())
    
    for strNameHist in dicHisto.keys(): 
      histRead = fRead.Get(strNameHist)
      histRead.Sumw2(False)
      dicHisto[ strNameHist ].Add(histRead)
    
    fRead.Close()
  
  strOut = os.path.join(strPathDraw, strDirHist, "hist_" + strDataset + ".root")
  fMerge = ROOT.TFile(strOut, "CREATE")
  
  dicInfo[ "entries" ] = "%lf"%fEntries
  dicInfo[ "step0" ] = "%lf"%fStep0
  dicInfo[ "step1" ] = "%lf"%fStep1
  
  for strKey in dicInfo.keys():  ROOT.TNamed(strKey, dicInfo[ strKey ]).Write()
  for strKey in dicHisto.keys(): dicHisto[ strKey ].Write()
  
  fMerge.Write()
  fMerge.Close()
  
  print "%s is done"%(strDataset)
else:
  import getopt, time, datetime
  
  strSpellNoDup  = "if [ `ps -ef | grep $UID | grep pytho[n] | grep draw_histos.py | wc -l` -gt 1 ]; "
  strSpellNoDup += "then ./cannot_exists_one 2> /dev/null ; else ls > /dev/null ; fi"
  
  if os.system(strSpellNoDup) != 0: 
    print "You are already running this code."
    sys.exit(0)
  
  strNameTree = "event"
  channel = 2 # 0 : all, 1 : el, 2 : mu, 3 : el + mu, (-) : anti
  allcharge = False
  
  #weight = "1"
  weight = "genweight * puweight * mueffweight * btagweight"
  #weight = "genweight * puweight * btagweight"
  #weight = "puweight"
  cut  = "step >= 4"
  #cut += " && ( njet - nbjet ) >= 1 && nbjet >= 1"
  #cut += " && njet == 2 && nbjet == 1"
  #cut += " && lep.Pt() >= 40"
  #cut += " && jet1.Pt() >= 40 && bjet1.Pt() >= 40"
  #cut += " && met > 50"
  
  strCutAdd = ""
  
  #bMulticore = False
  bMulticore = True
  bDefaultListSet = True
  
  strHelpHowto = "Usage : ./[name of this py file] [dir name of roots] " + \
    "-a <channel> -c <cut> -w <weight> -n <local running> -l <listSet JSON file>"
  try:
    opts, args = getopt.getopt(sys.argv[ 2: ], "hnea:c:d:w:l:", 
      ["channel", "cut", "cutadd", "allcharge", "weight", "listset"])
  except getopt.GetoptError:          
    sys.stderr.write(strHelpHowto + "\n")
    sys.exit(2)
  
  for opt, arg in opts:
    if opt == '-h':
      print strHelpHowto
      sys.exit()
    elif opt in ("-a", "--channel"):
      channel = int(arg)
    elif opt in ("-c", "--cut"):
      cut = arg
    elif opt in ("-d", "--cutadd"):
      strCutAdd = arg
    elif opt in ("-e", "--allcharge"): 
      allcharge = True
    elif opt in ("-w", "--weight"):
      weight = arg
    #elif opt in ("-m"):
    #  bMulticore = True
    elif opt in ("-n"):
      bMulticore = False
    elif opt in ("-l", "--listset"):
      bDefaultListSet = False
      if not arg.startswith(strPathDraw): strPathListSet = os.path.join(strPathDraw, arg)
      else:                               strPathListSet = arg
  
  if channel != 0: 
    strSign = "" if channel > 0 else "-"
    if abs(channel) == 1: 
      cut += " && trig_e > 0"
      cut += " && lep_pid == %s11"%strSign if not allcharge else " && abs(lep_pid) == 11"
    if abs(channel) == 2: 
      cut += " && trig_m > 0"
      cut += " && lep_pid == %s13"%strSign if not allcharge else " && abs(lep_pid) == 13"
    if abs(channel) == 3: 
      cut += " && ( lep_pid == %s11 || lep_pid == %s113)"%(strSign, strSign)
  
  cut += " && " if len(strCutAdd) > 0 else ""
  cut += strCutAdd
  
  strSrcPath = "/xrootd/store/user/quark2930/singletop/%s/"%sys.argv[ 1 ]
  
  if not os.path.exists(strSrcPath): 
    sys.stderr.write("Error: cannot find the directory of nano ntuples\n")
    sys.exit(1)
  
  dicSet = json.load(open(strPathListSet))
  
  # Loading name of datasets
  listSets = []
  if channel == 0 or ( abs(channel) & 1 ) != 0: listSets += dicSet[ "listDatasets" ][ "RDEL" ]
  if channel == 0 or ( abs(channel) & 2 ) != 0: listSets += dicSet[ "listDatasets" ][ "RDMU" ]
  
  listSets += dicSet[ "listDatasets" ][ "SIG" ]
  listSets += dicSet[ "listDatasets" ][ "BKG" ]
  if abs(channel) == 1: listSets += dicSet[ "listDatasets" ][ "BKG_EL" ]
  if abs(channel) == 2: listSets += dicSet[ "listDatasets" ][ "BKG_MU" ]
  
  dicOut = {}
  strDirHist = datetime.datetime.now().strftime("%y%m%d_%H%M%S")
  if not os.path.exists(os.path.join(strPathDraw, strDirHist)): 
    os.makedirs(os.path.join(strPathDraw, strDirHist))
    os.makedirs(os.path.join(strPathDraw, strDirHist, "logs"))
  
  nMaxProc = 20
  
  strJDSTemplate = """executable = parallel_draw_histos.py
universe   = vanilla
requirements = OpSysMajorVer == 6
arguments=%(flag)s %(json)s

log = %(path)s/condor.log

getenv     = True
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
output = %(path)s/logs/job_%(name)s.log
error = %(path)s/logs/job_%(name)s.err
queue 1
  """
  
  nNumJob = 0
  dicNumFiles = {}
  
  for strDataset in listSets: 
    listFiles = os.listdir(os.path.join(strSrcPath, "dir_" + strDataset))
    dicNumFiles[ strDataset ] = len(listFiles)
    
    strDirDataset = os.path.join(strPathDraw, strDirHist, strDataset)
    if not os.path.exists(strDirDataset): 
      os.makedirs(strDirDataset)
    
    for strFilename in listFiles: 
      dicOut[ "dataset" ] = strDataset
      dicOut[ "histname" ] = "hist_" + strDataset
      #dicOut[ "rootfile" ] = os.path.join(strSrcPath, strDataset + ".root")
      dicOut[ "rootfile" ] = [strSrcPath, "dir_" + strDataset, strFilename]
      dicOut[ "treename" ] = strNameTree
      dicOut[ "weight" ] = weight
      dicOut[ "cut" ] = cut
      dicOut[ "res_dir" ] = strDirHist
      
      if not bDefaultListSet: dicOut[ "listset" ] = strPathListSet
      
      if not bMulticore: 
        while True:
          nDone = len([ s for s in os.listdir(os.path.join(strPathDraw, strDirHist)) if ".done" in s ])
          if nDone + nMaxProc > nNumJob: break
          time.sleep(1)
        
        os.system("python %s %s '%s' &"%(strPathThis, strTypeArgOneRoot, 
          json.dumps(dicOut).replace("\"", "(+(")))
        nNumJob += 1
      else: 
        #dicOut[ "rootfile" ][ 0 ].replace("xrootd", "root://cms-xrdr.sdfarm.kr:1094///xrd")
        dicOut[ "rootfile" ][ 0 ] = "root://cms-xrdr.sdfarm.kr:1094///xrd" + dicOut[ "rootfile" ][ 0 ][ 7: ]
        
        strJDS = strJDSTemplate%{"flag": strTypeArgOneRoot, 
          "json": json.dumps(dicOut).replace("\"", "(+("), 
          "path": strDirHist, "name": strDataset + "_" + strFilename.replace(".root", "")}
        
        os.system("printf '%s' | condor_submit > /dev/null"%strJDS)
    
    if not bMulticore: 
      os.system("python %s %s %s %s %i &"%(strPathThis, strTypeArgMerger, 
        strDirHist, strDataset, dicNumFiles[ strDataset ]))
  
  if bMulticore: 
    print "All jobs have been thrown"
    
    for i, strDataset in enumerate(listSets): 
      os.system("python %s %s %s %s %i"%(strPathThis, strTypeArgMerger, 
        strDirHist, strDataset, dicNumFiles[ strDataset ]))
  
  while True:
    nDone = len([ s for s in os.listdir(os.path.join(strPathDraw, strDirHist)) if ".root" in s ])
    #if nDone >= nNumJob: break
    if nDone >= len(listSets): break
    time.sleep(1)
  
  print strDirHist
  os.system("date")


