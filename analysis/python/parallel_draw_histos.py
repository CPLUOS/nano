#!/usr/bin/env python


import ROOT, json, os, sys, copy
import time, datetime


# It returns a list for my customized 'sscanf' for paths (why there is no built-in sscanf in python...)
# The return is a list of lists in the form of the following: 
#   ["the followed pattern", "name", "type", "the following pattern"]
# For example, if strFormat = "dir_%(dataset)s/res_%(idx)i.root", then the return will be
#   [["dir_", "dataset", "s", "/res_"], ["/res_", "idx", "i", ".root"], [".root", ""]]
# Note that "the following pattern" coincides with "the followed pattern" of the next item
# This function is totally for the pre-process of GetInfoFromFormat()
def initFormat(strFormat): 
  listFirst = strFormat.split("%") # Because the input is for a path, we don't need to care about %%
  listItems = []
  strHeadItem = listFirst[ 0 ] # The first of "the followed pattern"
  
  for strItem in listFirst[ 1: ]: 
    if strItem[ 0 ] != '(': return None # Format error!
    nEnd = strItem.find(')')
    if nEnd < 0: return None # Format error!
    
    strEndItem = strItem[ nEnd + 2: ] # Extracting "the following pattern"
    listItems.append([strHeadItem, strItem[ 1: nEnd ], strItem[ nEnd + 1: nEnd + 2 ], strEndItem])
    
    # "The following pattern" of this item will be "the followed pattern" in the next one
    strHeadItem = strEndItem
  
  listItems.append([strHeadItem, ""]) # The last item; no format matching
  
  return listItems


# From the result of initFormat we can extract informations from a string as like sscanf
# The following is an example: 
#   listA = initFormat("dir_%(dataset)s/res_%(idx)i.root")
#   dicB = GetInfoFromFormat("dir_SingleTop_t-channel/res_7.root", listA)
#   # The return is {"dataset": "SingleTop_t-channel", "idx": 7}
# These functions are used for extracting dataset name and the index of a given file path; 
#  see parallel_draw_histos.SetSrcPathFormat()
def GetInfoFromFormat(strIn, listFormat): 
  strCurr = strIn
  dicRes = {}
  
  for listItem in listFormat: 
    nIdxStartFind = strCurr.find(listItem[ 0 ]) # Seeking the position of "the followed pattern"
    if nIdxStartFind < 0: return None # Format error; failure to seek "the followed pattern"
    
    if listItem[ 1 ] == "": break # This is the end
    
    nIdxEndFind = strCurr.find(listItem[ 3 ]) # Seeking the position of "the following pattern"
    # Format error; failure to seek "the followed pattern"
    if nIdxEndFind < 0 or nIdxEndFind < nIdxStartFind: return None
    
    # Extracting the value; and changing the type if necessary
    valGet = strCurr[ nIdxStartFind + len(listItem[ 0 ]): nIdxEndFind ]
    if listItem[ 2 ] == "i": valGet = int(valGet)
    if listItem[ 2 ] == "f": valGet = float(valGet)
    
    # Inserting the value; and checking the consistency
    strKey = listItem[ 1 ]
    if strKey in dicRes and dicRes[ strKey ] != valGet: return None # Error: no consistency
    dicRes[ strKey ] = valGet
    
    # For next item
    strCurr = strCurr[ nIdxEndFind: ]
  
  return dicRes


def valToName(strVar): 
  return strVar.replace(".", "").replace("(", "").replace(")", "").replace(" ", "").replace("*", "").replace(":", "_").replace("+", "_add_").replace("/", "_div_")


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


# Usage of the following class:
#   from [PATH].parallel_draw_histos import parallel_draw_histos
#   paraMain = parallel_draw_histos()
#   
#   paraMain.SetSrcPath("[PATH OF LIST FILE]")
#   paraMain.SetSrcPathFormat("[PATH/...%(dastaset)s...%(idx)s...]")
#   paraMain.SetNameTree("[TREE NAME]")
#   paraMain.SetDirHistName("[NAME]")
#   paraMain.SetCut("[CUT]")
#   paraMain.SetWeight("[WEIGHT]")
#   paraMain.SetDatasets([LIST OF DATASETS])
#   paraMain.SetVars({DICT OF VARIABLES}) # See the description for this function
#   paraMain.SetPathDraw("[PATH]")
#   
#   paraMain.InitRun()
#   paraMain.RunOnCluster() or paraMain.RunOnWorknode()


class parallel_draw_histos:
  def __init__(self):
    self.strPathThis = os.path.join("%s/src/nano/analysis/python"%os.environ[ "CMSSW_BASE" ], 
      "parallel_draw_histos.py")
    
    self.strTypeArgOneRoot = "--oneroot"
    self.strTypeArgMerger = "--merger"
    
    self.strHeadDirXRD = "/xrootd"
    self.strAddrXRD = "root://cms-xrdr.sdfarm.kr:1094///xrd"
    
    self.strFriendPath = None
    self.strNameFriendTree = None
    
    self.nMaxProc = 20
   
  
  ### Variables for works
  ### Thus you MUST execute ALL the following Set[...] functions at least at once
  
  # strDirHist: The name of the directory in which the result will be
  # If you just execute this without any argument, it will set the name as the current time
  def GetDirHistName(self): return self.strDirHist
  def SetDirHistName(self, strName = ""): 
    self.strDirHist = strName if strName != "" else datetime.datetime.now().strftime("%y%m%d_%H%M%S")
  
  # strSrcPath: The path (MUTE BE in the xrd store) of a file containing the list of source ntuples
  # NOTE: All source root file are in directories of which the paths are in the following forms: 
  #   [Somewhere in XRD storage]...%(dataset)s.../[NAME OF THE ROOT FILE].root
  def GetSrcPath(self): return self.strSrcPath
  def SetSrcPath(self, strSrcPath): self.strSrcPath = strSrcPath
  
  # strSrcPathFormat: The pattern of paths of ntuple files
  # It is for extraction of the dataset name and the index of the file
  # For example, one can set this as 
  #   "root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/%(srcname)s/dir_%(dataset)s/res_%(idx)s.root", 
  # where the part of dataset name and the index number is replaced by %(dataset)s and %(idx)s, respectively.
  # These replacements are required. 
  # On the other hand, other formats, such as %(srcname)s (other name than 'srcname' is available), is not necessary.
  # However, to avoid cumbersome things (e.g., the name of main directory, which can be changed frequently)
  # the user can replace such parts by other patterns; %(srcname)s in this example is for that.
  # By the pattern, from the following path
  #   root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/quark2930/singletop/181105_1/dir_SingleTop_t-channel/res_7.root
  # "dataset" would be "SingleTop_t-channel" and "idx" be "7"
  # also "srcname" would be "quark2930/singletop/181105_1", which is not used.
  def GetSrcPathFormat(self): return self.strSrcPathFormat
  def SetSrcPathFormat(self, strSrcPathFormat): self.strSrcPathFormat = strSrcPathFormat
  
  # strNameTree: The name of the tree in the ntuple
  def GetNameTree(self): return self.strNameTree
  def SetNameTree(self, strNameTree): self.strNameTree = strNameTree
  
  # strCutMC: Cut condition for MC samples
  def GetCutMC(self): return self.strCutMC
  def SetCutMC(self, strCut): self.strCutMC = strCut
  
  # strCutRD: Cut condition for RD samples
  def GetCutRD(self): return self.strCutRD
  def SetCutRD(self, strCut): self.strCutRD = strCut
  
  # strWeight: Weight furmula for MC samples
  def GetWeightMC(self): return self.strWeightMC
  def SetWeightMC(self, strWeight): self.strWeightMC = strWeight
  
  # strWeight: Weight furmula for RD samples
  def GetWeightRD(self): return self.strWeightRD
  def SetWeightRD(self, strWeight): self.strWeightRD = strWeight
  
  # listDatasets: List of datasets
  def GetDatasets(self): return self.listDatasets
  def SetDatasets(self, listDatasets): self.listDatasets = listDatasets
  
  # dicVars: Dictionary containing infos for variables and histograms which will be drawn
  # 
  # The name of variable (formula is okay) and histogram in the source tree will be the key
  # For drawing a histogram of variables, the value contains a dictionary 
  #   with a list containing bin info with key "bin" (see makeHisto() above for how it will be given)
  # For drawing a histogram which is already drawn and stred in the source tree 
  #   just any simple string is okay for value; only name in the key is required
  # 
  # The following is an example: 
  # dicVars = {
  #   "var1": {"bin": [12, -0.5, 11.5]}, 
  #   "var2: var3": {"bin": [50, -1, 1, 50, 0, 200]}, 
  #   "var4.Pt()": {"bin": [50, 30, 230]}, 
  #   "var5 * var6.Eta()": {"bin": [50, -5, 5]}, 
  #   "histo_name": "histo"
  # }
  # 
  def GetVars(self): return self.dicVars
  def SetVars(self, dicVars): self.dicVars = dicVars
  
  # strPathDraw: The path for destination
  def GetPathDraw(self): return self.strPathDraw
  def SetPathDraw(self, strPathDraw): self.strPathDraw = strPathDraw
  
  
  ### The followings are for configuration of this parallelizer; no need to change
  
  # nMaxProc: The number of nodes you want to use; meaningless when you use cluster job
  def GetMaxProc(self): return self.nMaxProc
  def SetMaxProc(self, nMaxProc): self.nMaxProc = nMaxProc
  
  # strHeadDirXRD: The head path of xrd storage
  def GetHeadDirXRD(self): return self.strHeadDirXRD
  def SetHeadDirXRD(self, strHeadDirXRD): self.strHeadDirXRD = strHeadDirXRD
  
  # strAddrXRD: The address (in XRD protocol) of xrd storage
  def GetAddrXRD(self): return self.strAddrXRD
  def SetAddrXRD(self, strAddrXRD): self.strAddrXRD = strAddrXRD
  
  # strNameFriendTree: The name of the tree in the friend ntuple
  def GetNameFriendTree(self): return self.strNameFriendTree
  def SetNameFriendTree(self, strNameFriendTree): self.strNameFriendTree = strNameFriendTree
  
  # strFriendPath: The path (MUTE BE in the xrd store) in which the TMVA ntuples is located
  # NOTE: All source root file are in directories of which the paths are in the following forms: 
  #   [Somewhere in XRD storage]...%(dataset)s.../[NAME OF THE ROOT FILE].root
  def GetFriendPath(self): return self.strFriendPath
  def SetFriendPath(self, strFriendPath): self.strFriendPath = strFriendPath
  
  
  ### The followings are for internal variables
  
  # strFilenameDumpJSON: The path of dump JSON file
  def GetPathDumpJSON(self): return self.strFilenameDumpJSON
  def SetPathDumpJSON(self, strFilenameDumpJSON): self.strFilenameDumpJSON = strFilenameDumpJSON
  
  
  def InitRun(self): 
    # Preparing to throwing jobs: making basic directories
    if not os.path.exists(os.path.join(self.GetPathDraw(), self.GetDirHistName())): 
      os.makedirs(os.path.join(self.GetPathDraw(), self.GetDirHistName()))
      os.makedirs(os.path.join(self.GetPathDraw(), self.GetDirHistName(), "logs"))
    
    # Preparing to throwing jobs: making a card (in JSON format)
    dicOut = {}
    
    dicOut[ "dataset" ] = []
    dicOut[ "rootfile" ] = self.GetSrcPath()
    dicOut[ "filepattern" ] = self.GetSrcPathFormat()
    dicOut[ "treename" ] = self.GetNameTree()
    dicOut[ "cutMC" ] = self.GetCutMC()
    dicOut[ "cutRD" ] = self.GetCutRD()
    dicOut[ "weightMC" ] = self.GetWeightMC()
    dicOut[ "weightRD" ] = self.GetWeightRD()
    dicOut[ "Vars" ] = self.GetVars()
    dicOut[ "res_dir" ] = os.path.join(self.GetPathDraw(), self.GetDirHistName())
    
    strFriendPath = self.GetFriendPath()
    
    if strFriendPath is not None: 
      dicOut[ "friend" ] = self.GetFriendPath()
      dicOut[ "friendname" ] = self.GetNameFriendTree()
    
    self.listDatasetOfFile = []
    self.dicNumFiles = {}
    self.nTotalNumFiles = 0
    
    # Preparing to throwing jobs: inserting lists of source ntuples into the card
    fAllFiles = open(self.GetSrcPath())
    listAllFiles = fAllFiles.readlines()
    fAllFiles.close()
    
    # Initializing extraction of the pattern from paths
    listPattern = initFormat(self.GetSrcPathFormat())
    
    for strDataset in self.GetDatasets(): 
      listFiles = [ s for s in listAllFiles if GetInfoFromFormat(s, listPattern)[ "dataset" ] == strDataset ]
      
      self.listDatasetOfFile.append(strDataset)
      self.dicNumFiles[ strDataset ] = len(listFiles)
      self.nTotalNumFiles += len(listFiles)
      
      strDirDataset = os.path.join(self.GetPathDraw(), self.GetDirHistName(), strDataset)
      if not os.path.exists(strDirDataset): 
        os.makedirs(strDirDataset)
    
    self.SetPathDumpJSON(os.path.join(self.GetPathDraw(), self.GetDirHistName(), "config.json"))
    fWriteJSON = open(self.GetPathDumpJSON(), "w")
    json.dump(dicOut, fWriteJSON)
    fWriteJSON.close()
  
  
  def RunOnCluster(self): 
    # The card for condor submit
    strJDSTemplate  = "executable = " + self.strPathThis
    strJDSTemplate += """
universe   = vanilla
requirements = ( HasSingularity == true )
arguments=%(flag)s %(json)s $(Process)
accounting_group=group_cms
+SingularityImage = "/cvmfs/singularity.opensciencegrid.org/opensciencegrid/osgvo-el6:latest"
+SingularityBind = "/cvmfs, /cms, /share"

log = %(path)s/condor.log

getenv     = True
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
output = %(path)s/logs/job_$(Process).log
error = %(path)s/logs/job_$(Process).err
queue %(queue)i
    """
  
    # completing the submit card
    strJDS = strJDSTemplate%{"flag": self.strTypeArgOneRoot, 
      "json": self.GetPathDumpJSON(), "queue": self.nTotalNumFiles, "path": self.GetDirHistName()}
    
    # Submit!
    os.system("printf '%s' | condor_submit > /dev/null"%strJDS)
    
    print "All jobs have been thrown"
    
    dicMerge = self.dicNumFiles
    strDirDraw = os.path.join(self.GetPathDraw(), self.GetDirHistName())
    
    # Running codes for merging
    # Checking the progress (by counting the number of drawn root files) for each dataset, 
    # if it is done, then running the merger
    while len(dicMerge.keys()) > 0: 
      dicNext = {}
      
      for strDataset in dicMerge.keys(): 
        # Counting the result root files
        nNumDone = len([ s for s in os.listdir(strDirDraw) if strDataset in s and ".done" in s ])
        if nNumDone >= dicMerge[ strDataset ]: 
          # Merge!
          self.DoMerge(strDataset, self.dicNumFiles[ strDataset ])
          # If complete, it must be pulled out in the next dataset list
        else: 
          dicNext[ strDataset ] = dicMerge[ strDataset ]
      
      dicMerge = dicNext
      
    os.system("date")
    print self.GetDirHistName()
  
  
  def RunOnWorknode(self): 
    strDirHistFull = os.path.join(self.GetPathDraw(), self.GetDirHistName())
    nNumJob = 0
    
    for strDataset in self.listDatasetOfFile: 
      for iii in range(self.dicNumFiles[ strDataset ]): 
        while True:
          nDone = len([ s for s in os.listdir(strDirHistFull) if ".done" in s ])
          if nDone + self.nMaxProc > nNumJob: break
          time.sleep(1)
        
        os.system("python %s %s %s %i &"%(self.strPathThis, self.strTypeArgOneRoot, 
          self.GetPathDumpJSON(), nNumJob))
        nNumJob += 1
      
      # Merge!
      os.system("python %s %s %s %s %s %i &"%(self.strPathThis, self.strTypeArgMerger, 
        self.GetPathDraw(), self.strDirHist, strDataset, self.dicNumFiles[ strDataset ]))
    
    while True:
      nDone = len([ s for s in os.listdir(strDirHistFull) if ".root" in s ])
      #if nDone >= nNumJob: break
      if nDone >= len(self.dicNumFiles): break
      time.sleep(1)
    
    os.system("date")
    print self.GetDirHistName()
  
  
  def DoOneRoot(self, strPathJSON, nIdxJob): 
    fJSONQuery = open(strPathJSON, "r")
    dicMain = json.load(fJSONQuery)
    fJSONQuery.close()
    
    # Loading configurations
    strPathSrc = dicMain[ "rootfile" ].encode("ascii", "ignore")
    strPathFriend = None if "friend" not in dicMain else dicMain[ "friend" ].encode("ascii", "ignore")
    
    fListSrc = open(strPathSrc, "rb")
    #listInfoSrc = fListSrc.readlines()[ nIdxJob ].split()
    strSrc = fListSrc.readlines()[ nIdxJob ].splitlines(False)[ 0 ]
    fListSrc.close()
    
    strPattern = dicMain[ "filepattern" ].encode("ascii", "ignore")
    dicPattern = GetInfoFromFormat(strSrc, initFormat(strPattern))
    
    listInfoFriend = None
    if strPathFriend is not None: 
      fListFriend = open(strPathFriend, "rb")
      listInfoFriend = fListFriend.readlines()[ nIdxJob ].split()
      fListFriend.close()
    
    strDataset = dicPattern[ "dataset" ]
    
    if "Run" not in strDataset: 
      strWeight = dicMain[ "weightMC" ].encode("ascii", "ignore")
      strCut    = dicMain[ "cutMC" ].encode("ascii", "ignore")
    else: 
      strWeight = dicMain[ "weightRD" ].encode("ascii", "ignore")
      strCut    = dicMain[ "cutRD" ].encode("ascii", "ignore")
    
    strNoFile = dicPattern[ "idx" ]
    
    dicVar = dicMain[ "Vars" ]
    dicOutInfo = {}
    
    dicOutInfo[ "dataset" ] = strDataset
    dicOutInfo[ "ntuple" ] = strSrc
    dicOutInfo[ "no_file" ] = strNoFile
    
    dicOutInfo[ "weight" ] = strWeight
    dicOutInfo[ "cut" ] = strCut
    
    strNameFriend = None
    
    if strPathFriend is not None: 
      dicOutInfo[ "friend" ] = listInfoFriend[ 2 ]
      strNameFriend = dicMain[ "friendname" ].encode("ascii", "ignore")
    
    dicHist = {}
    
    # Loading the source ntuple root file
    strNameTree = dicMain[ "treename" ].encode("ascii", "ignore")
    
    fMain = ROOT.TFile.Open(dicOutInfo[ "ntuple" ])
    nTrial = 0
    
    while nTrial < 5: 
      if fMain: break
      
      time.sleep(10)
      fMain = ROOT.TFile.Open(dicOutInfo[ "ntuple" ])
      nTrial += 1
    
    if nTrial >= 5: 
      sys.stderr.write("#%i: Failed to open %s"%(nIdxJob, dicOutInfo[ "ntuple" ]))
      return
    
    tree = fMain.Get(strNameTree)
    
    # Adding the friend tree
    # If a file for friend tree is given, then it will be added
    # Especially, this is a powerful equipment for applying TMVA results
    if strPathFriend is not None: 
      tFriend = tree.AddFriend(strNameFriend, dicOutInfo[ "friend" ])
      nTrial = 0
      
      while nTrial < 5: 
        if not ( not tFriend.GetTree() ): break
        
        time.sleep(10)
        tFriend = tree.AddFriend(strNameTree, dicOutInfo[ "friend" ])
        nTrial += 1
      
      if nTrial >= 5: 
        sys.stderr.write("#%i: Failed to open %s"%(nIdxJob, dicOutInfo[ "friend" ]))
        return
    
    # Applying the cut
    # This way applies the cut only once for several histograms, so enhances the performances
    # Btw, if no entries in the existing tree, the following doesn't work, and also it is not needed
    if tree.GetEntries() > 0: 
      strNameAfterCut = strNameTree.replace("/", "") + "_aftercut"
      tree.Draw(">>" + strNameAfterCut, strCut, "entrylist")
      tree.SetEntryList(ROOT.gDirectory.Get(strNameAfterCut))
    
    # Drawing histograms
    for strVar in dicVar.keys(): 
      strNameVar = valToName(strVar)
      strNameHist  = "hist_%(dataset)s_%(var)s"%{"dataset": strDataset, "var": strNameVar}
      
      # Two cases: histogram drawn from the loaded tree, or histogram which exists already
      if type(dicVar[ strVar ]) is dict: 
        hTmp = makeHisto(strNameHist, strNameHist, dicVar[ strVar ][ "bin" ])
        tree.Project(strNameHist, strVar, strWeight)
        
        dicHist[ strVar ] = copy.deepcopy(hTmp)
        
        if dicHist[ strVar ].GetSumw2N() == 0:
          dicHist[ strVar ].Sumw2()
      else: 
        hTmp = copy.deepcopy(fMain.Get(strVar.encode("ascii", "ignore")))
        hTmp.SetName(strNameHist)
        
        dicHist[ strVar ] = hTmp
    
    # Now, ready to write
    strOutRoot = os.path.join(dicMain[ "res_dir" ], strDataset, strNoFile + ".root")
    fRootFile = ROOT.TFile(strOutRoot, "NEW")
    
    # Writing infos; the first is for string infos and the second is the histograms
    for strKey in dicOutInfo.keys(): ROOT.TNamed(strKey, dicOutInfo[ strKey ]).Write()
    for strVar in dicHist.keys(): dicHist[ strVar ].Write()
    
    # DO NOT FORGET!
    fRootFile.Write()
    fRootFile.Close()
    
    # Notifying that I'm done, in somewhat primitive way, but definite way
    os.system("touch %s"%(os.path.join(dicMain[ "res_dir" ], 
      strDataset + "_" + strNoFile + ".done")))
  
  
  def DoMerge(self, strDataset, nNumJobs): 
    import shutil, time
    
    if nNumJobs <= 0: 
      print "Wrong arguments"
      sys.exit(1)
    
    strDirDraw = os.path.join(self.GetPathDraw(), self.GetDirHistName())
    
    while True:
      if len([ s for s in os.listdir(strDirDraw) if strDataset in s and ".done" in s ]) >= nNumJobs: break
      time.sleep(1)
    
    strDirDataset = os.path.join(self.GetPathDraw(), self.GetDirHistName(), strDataset)
    listRoot = [ s for s in os.listdir(strDirDataset) if ".root" in s ]
    
    fOrg = ROOT.TFile(os.path.join(strDirDataset, listRoot[ 0 ]))
    listRoot.pop(0)
    
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
      
      for strNameHist in dicHisto.keys(): 
        histRead = fRead.Get(strNameHist)
        histRead.Sumw2(False)
        dicHisto[ strNameHist ].Add(histRead)
      
      fRead.Close()
    
    strOut = os.path.join(self.GetPathDraw(), self.GetDirHistName(), "hist_" + strDataset + ".root")
    fMerge = ROOT.TFile(strOut, "CREATE")
    
    for strKey in dicInfo.keys():  ROOT.TNamed(strKey, dicInfo[ strKey ]).Write()
    for strKey in dicHisto.keys(): dicHisto[ strKey ].Write()
    
    fMerge.Write()
    fMerge.Close()
    
    print "%s is done"%(strDataset)


if __name__ == "__main__": 
  paraMain = parallel_draw_histos()
  
  if len(sys.argv) < 2: 
    sys.stderr.write("Error: Too few arguments\n")
    sys.exit(1)
  
  if sys.argv[ 1 ] == paraMain.strTypeArgOneRoot: 
    if len(sys.argv) < 4: 
      sys.stderr.write("Error: Too few arguments\n")
      sys.exit(1)
    
    paraMain.DoOneRoot(sys.argv[ 2 ], int(sys.argv[ 3 ]))
  elif sys.argv[ 1 ] == paraMain.strTypeArgMerger: 
    if len(sys.argv) < 6: 
      sys.stderr.write("Error: Too few arguments\n")
      sys.exit(1)
    
    paraMain.SetPathDraw(sys.argv[ 2 ])
    paraMain.SetDirHistName(sys.argv[ 3 ])
    paraMain.DoMerge(sys.argv[ 4 ], int(sys.argv[ 5 ]))


