import ROOT, json, os, getopt, sys, copy, time, datetime


def valToName(strVar): 
  return strVar.replace(".", "").replace("(", "").replace(")", "").replace(" ", "")


strPathDraw = "%s/src/nano/analysis/test/singletop/draw"%os.environ[ "CMSSW_BASE" ]
dicSet = json.load(open(os.path.join(strPathDraw, "listSet.json")))


if len(sys.argv) > 1 and sys.argv[ 1 ] == "part":   
  dicVar = dicSet[ "Vars" ]
  dicOutInfo = {}
  
  dicMain = json.loads(sys.argv[ 2 ])
  
  strWeight = dicMain[ "weight" ].encode("ascii", "ignore")
  strCut    = dicMain[ "cut" ].encode("ascii", "ignore")
  #strFullCut = "( %s ) * ( %s )"%(strWeight, strCut)
  
  dicOutInfo[ "ntuple" ] = dicMain[ "rootfile" ].encode("ascii", "ignore")
  dicOutInfo[ "weight" ] = strWeight
  dicOutInfo[ "cut" ] = strCut
  
  strNameTree = dicMain[ "treename" ].encode("ascii", "ignore")
  fMain = ROOT.TFile(dicMain[ "rootfile" ].encode("ascii", "ignore"))
  tree = fMain.Get(strNameTree)
  
  dicHist = {}
  strOutRoot = os.path.join(strPathDraw, dicMain[ "res_dir" ], 
    dicMain[ "histname" ].encode("ascii", "ignore") + ".root")
  
  # Getting # of whole entries
  dicOutInfo[ "entries" ] = "%i"%( fMain.Get("weight").Integral() )
  
  print "%s - %s : Starting cutting"%(datetime.datetime.now().strftime("%H:%M:%S"), dicMain[ "histname" ])
  
  # Applying the cut
  strNameAfterCut = strNameTree + "_aftercut"
  #tree.Draw(">>" + strNameAfterCut, strFullCut, "entrylist")
  tree.Draw(">>" + strNameAfterCut, strCut, "entrylist")
  tree.SetEntryList(ROOT.gDirectory.Get(strNameAfterCut))
  
  print "%s - %s : Cutting is done"%(datetime.datetime.now().strftime("%H:%M:%S"), dicMain[ "histname" ])
  
  # Drawing histograms
  for strVar in dicVar.keys(): 
    #if strVar != "lep1.Pt()" or strVar != "njet" or strVar != "met": continue
    strNameVar = valToName(strVar)
    strNameHist = dicMain[ "histname" ].encode("ascii", "ignore") + "_" + strNameVar
    listBin = dicVar[ strVar ][ "bin" ]
    
    hTmp = ROOT.TH1F(strNameHist, strNameHist, listBin[ 0 ], listBin[ 1 ], listBin[ 2 ])
    #tree.Project(strNameHist, strVar)
    tree.Project(strNameHist, strVar, strWeight)
    
    dicHist[ strVar ] = copy.deepcopy(hTmp)
    
    if dicHist[ strVar ].GetSumw2N() == 0:
      dicHist[ strVar ].Sumw2()
  
  print "%s - %s : Drawing is done"%(datetime.datetime.now().strftime("%H:%M:%S"), dicMain[ "histname" ])
  
  fRootFile = ROOT.TFile(strOutRoot, "NEW")
  
  # Writing the string infos
  for strKey in dicOutInfo.keys(): 
    strVal = dicOutInfo[ strKey ]
    dicOutInfo[ strKey ] = ROOT.TNamed(strKey, strVal)
    dicOutInfo[ strKey ].Write()
  
  for strVar in dicHist.keys(): dicHist[ strVar ].Write()
  
  fRootFile.Write()
  fRootFile.Close()
else:
  strSpellNoDup  = "if [ `ps -ef | grep $UID | grep pytho[n] | grep draw_histos.py | wc -l` -gt 1 ]; "
  strSpellNoDup += "then ./cannot_exists_one 2> /dev/null ; else ls > /dev/null ; fi"
  
  if os.system(strSpellNoDup) != 0: 
    print "You are already running this code."
    sys.exit(0)
  
  strNameTree = "event"
  channel = 2 # 0 : all, 1 : el, 2 : mu, 3 : el + mu, (-) : anti
  
  #weight = "1"
  weight = "genweight * puweight * mueffweight * btagweight"
  cut  = "step >= 4"
  #cut  = " && ( njet - nbjet ) >= 1 && nbjet >= 1"
  cut += " && lep.Pt() >= 26"
  cut += " && jet1.Pt() >= 40 && bjet.Pt() >= 40"
  cut += " && met > 45"
  #cut += " && njet == 2"
  #cut += " && jetC < 0.03 && ( 100 < top1.M() && top1.M() < 200 ) && njet <= 3 && abs(jet1.Eta()) > 1.2"
  
  strCutAdd = ""
  
  strHelpHowto = "Usage : ./[name of this py file] [dir name of roots] -a <channel> -c <cut> -w <weight>"
  try:
    opts, args = getopt.getopt(sys.argv[ 2: ], "hd:a:c:w:", ["channel", "cut", "cutadd", "weight"])
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
    elif opt in ("-w", "--weight"):
      weight = arg
  
  if channel != 0: 
    strSign = "" if channel > 0 else "-"
    if abs(channel) == 1: 
      cut += " && lep1_pid == %s11"%strSign
    if abs(channel) == 2: 
      cut += " && lep1_pid == %s13"%strSign
    if abs(channel) == 3: 
      cut += " && ( lep1_pid == %s11 || lep1_pid == %s113)"%(strSign, strSign)
  
  cut += " && " if len(strCutAdd) > 0 else ""
  cut += strCutAdd
  
  strSrcPath = "/xrootd/store/user/quark2930/singletop/%s/"%sys.argv[ 1 ]
  
  if not os.path.exists(strSrcPath): 
    sys.stderr.write("Error: cannot find the directory of nano ntuples\n")
    sys.exit(1)
  
  # Loading name of datasets
  listSetsPre = []
  if channel == 0 or ( abs(channel) & 1 ) != 0: listSetsPre += dicSet[ "listDatasets" ][ "RDEL" ]
  if channel == 0 or ( abs(channel) & 2 ) != 0: listSetsPre += dicSet[ "listDatasets" ][ "RDMU" ]
  listSetsPre += dicSet[ "listDatasets" ][ "SIG" ]
  listSetsPre += dicSet[ "listDatasets" ][ "BKG" ]
  
  # Expanding the names; for 'too many' datasets
  # 'Too many' dataset is divided in two parts.
  listSets = []
  
  for strSet in listSetsPre: 
    nTooMany = 0
    for strTooMany in dicSet[ "listDatasets" ][ "TooMany" ]: 
      if strTooMany in strSet: 
        nTooMany = 1
        break
    
    if nTooMany == 0: 
      listSets.append(strSet)
    else: 
      listSets.append(strSet + "1")
      listSets.append(strSet + "2")
  
  dicOut = {}
  strDirHist = datetime.datetime.now().strftime("%y%m%d_%H%M%S")
  if not os.path.exists(os.path.join(strPathDraw, strDirHist)): os.makedirs(os.path.join(strPathDraw, strDirHist))
  
  nMaxProc = 20
  
  for i, strDataset in enumerate(listSets): 
    dicOut[ "histname" ] = "hist_" + strDataset
    dicOut[ "rootfile" ] = os.path.join(strSrcPath, strDataset + ".root")
    dicOut[ "treename" ] = strNameTree
    dicOut[ "weight" ] = weight
    dicOut[ "cut" ] = cut
    dicOut[ "res_dir" ] = strDirHist
    
    while True:
      nDone = len([ s for s in os.listdir(os.path.join(strPathDraw, strDirHist)) if ".root" in s ])
      if nDone + nMaxProc > i: break
      time.sleep(1)
    
    os.system("python %s part '%s' &"%(os.path.join(strPathDraw, "draw_histos.py"), json.dumps(dicOut)))
  
  while True:
    nDone = len([ s for s in os.listdir(os.path.join(strPathDraw, strDirHist)) if ".root" in s ])
    if nDone >= len(listSets): break
    time.sleep(1)
  
  print strDirHist


