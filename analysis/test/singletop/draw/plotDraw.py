import ROOT, json, os, getopt, sys, time
import nano.analysis.CMS_lumi
from nano.analysis.histoHelper import *


def valToName(strVar): 
  return strVar.replace(".", "").replace("(", "").replace(")", "").replace(" ", "").replace("*", "")


listRDRun = [
  "2016B", "2016C", "2016D", "2016E", "2016F", "2016G", "2016H", 
]

#datalumi =  35900 #35.9fb-1
datalumi =  33274 #33.3fb-1

strSigTitle = "t-channel"
nSigColor = 2

strPathDraw = "%s/src/nano/analysis/test/singletop/draw"%os.environ[ "CMSSW_BASE" ]
dicSetDef = json.load(open(os.path.join(strPathDraw, "listSet.json")))
datasets = json.load(open("%s/src/nano/nanoAOD/data/dataset/dataset.json" % os.environ['CMSSW_BASE']))

channel = 2 # 0 : all, 1 : el, 2 : mu, 3 : el + mu, (-) : anti
plotvar = "lep.Pt()"
dolog = False
onlyPrintNum = False

try:
  opts, args = getopt.getopt(sys.argv[2:],"nlp:",["plotvar","dolog"])
except getopt.GetoptError:          
  print 'Usage : ./[name of the py file] [histogram directory] -p <plotvar> -d <dolog>'
  sys.exit(2)

for opt, arg in opts:
  if opt == '-h':
    print 'Usage : ./[name of the py file] [histogram directory] -p <plotvar> -d <dolog>'
    sys.exit()
  elif opt in ("-n", "--numberonly"):
    onlyPrintNum = True
    plotvar = "step*0+0.5"
  elif opt in ("-p", "--plotvar"):
    plotvar = arg
  elif opt in ("-l", "--dolog"):
    dolog = True

if plotvar not in dicSetDef[ "Vars" ].keys(): 
  print "This variable is not registered in listSet.json. Please add info about this variable."
  sys.exit(1)

binning = dicSetDef[ "Vars" ][ plotvar ][ "bin" ]
x_name  = dicSetDef[ "Vars" ][ plotvar ][ "xaxis" ]
y_name  = dicSetDef[ "Vars" ][ plotvar ][ "yaxis" ]

fMinPlot = None if "min" not in dicSetDef[ "Vars" ][ plotvar ] else dicSetDef[ "Vars" ][ plotvar ][ "min" ]
fMaxPlot = None if "max" not in dicSetDef[ "Vars" ][ plotvar ] else dicSetDef[ "Vars" ][ plotvar ][ "max" ]

strDirHist = sys.argv[ 1 ].replace("/", "")

# Determining what is the channel of the given histograms
listDirInFile = os.listdir(os.path.join(strPathDraw, strDirHist))
nIsChEl = 0
nIsChMu = 0

for strFile in listDirInFile: 
  if "Run" in strFile: 
    if "Electron" in strFile: nIsChEl = 1
    if "Muon"     in strFile: nIsChMu = 1

channel = nIsChEl + 2 * nIsChMu

strVarName = valToName(plotvar)

dicSet= {}

fStep0Entries = 0.0
fStep1Entries = 0.0

dicSet[ "SingleTbar_t-channel" ] = {"TYPE": "SIG"}
dicSet[ "SingleTop_t-channel" ] = {"TYPE": "SIG"}
for strKey in dicSetDef[ "listDatasets" ][ "BKG" ]: dicSet[ strKey ] = {"TYPE": "BKG"}

if channel == 0 or ( abs(channel) & 1 ) != 0: 
  for strKey in dicSetDef[ "listDatasets" ][ "BKG_EL" ]: dicSet[ strKey ] = {"TYPE": "BKG"}
  for strKey in dicSetDef[ "listDatasets" ][ "RDEL" ]: dicSet[ strKey ] = {"TYPE": "RD"}
if channel == 0 or ( abs(channel) & 2 ) != 0: 
  for strKey in dicSetDef[ "listDatasets" ][ "BKG_MU" ]: dicSet[ strKey ] = {"TYPE": "BKG"}
  for strKey in dicSetDef[ "listDatasets" ][ "RDMU" ]: dicSet[ strKey ] = {"TYPE": "RD"}

for strKey in dicSet.keys():
  strNameRoot = "hist_" + strKey.encode("ascii", "ignore")
  strNameHist = strNameRoot + "_" + strVarName
  
  dicSet[ strKey ][ "hist" ] = ROOT.TH1F(strNameHist + "_new", strNameHist + "_new", 
    binning[ 0 ], binning[ 1 ], binning[ 2 ])
  
  dicSet[ strKey ][ "entries" ] = 0
  dicSet[ strKey ][ "step0" ] = 0
  dicSet[ strKey ][ "step1" ] = 0
  
  fRoot = ROOT.TFile(os.path.join(strPathDraw, strDirHist, strNameRoot + ".root"))
  histLoad = fRoot.Get(strNameHist)
  histLoad.Sumw2(False)
  
  dicSet[ strKey ][ "hist" ].Add(histLoad)
  
  dicSet[ strKey ][ "entries" ] += float(fRoot.Get("entries").GetTitle())
  dicSet[ strKey ][ "step0" ] += float(fRoot.Get("step0").GetTitle())
  dicSet[ strKey ][ "step1" ] += float(fRoot.Get("step1").GetTitle())
  
  fRoot.Close()
  
  dicSet[ strKey ][ "finddataset" ] = findDataSet(strKey, datasets)
  
  # Scaling; the actual entry must be xsec * lumi ( * additional weight if there is)
  scale = 1.0
  if dicSet[ strKey ][ "TYPE" ] != "RD": 
    scale = datalumi * dicSet[ strKey ][ "finddataset" ][ "xsec" ] / dicSet[ strKey ][ "entries" ]
    if strKey in dicSetDef[ "Additional" ].keys(): 
      scale = scale * dicSetDef[ "Additional" ][ strKey ][ "AddWeight" ]
  else: 
    fStep0Entries += dicSet[ strKey ][ "step0" ]
    fStep1Entries += dicSet[ strKey ][ "step1" ]
  
  dicSet[ strKey ][ "scale" ] = scale
  dicSet[ strKey ][ "hist" ].Scale(scale)

# Calculating the luminosity
#for strRun in listRDRun: 

# To consider the lumi check
datalumi *= fStep1Entries / fStep0Entries

listMCSig = []
listMC = []
histRD  = ROOT.TH1F("histRD",  "histRD",  binning[ 0 ], binning[ 1 ], binning[ 2 ])
dicRDRun = {}

for strRDRun in listRDRun: 
  dicRDRun[ strRDRun ] = {}
  dicRDRun[ strRDRun ][ "hist" ] = ROOT.TH1F("histRD_" + strRDRun, "histRD", 
    binning[ 0 ], binning[ 1 ], binning[ 2 ])

# Inserting MC plots into stacks, merging RD plots, colouring, etc.
for strKey in dicSet.keys():
  dicSetInfo = dicSet[ strKey ][ "finddataset" ]
  
  if dicSet[ strKey ][ "TYPE" ] == "SIG": 
    # Decoration; titling, colouring
    dicSet[ strKey ][ "hist" ].SetTitle(strSigTitle)
    dicSet[ strKey ][ "hist" ].SetFillColor(nSigColor)
    dicSet[ strKey ][ "hist" ].SetLineColor(nSigColor)
    
    # Appending
    listMCSig.append(dicSet[ strKey ][ "hist" ])
  elif dicSet[ strKey ][ "TYPE" ] == "BKG": 
    # Loading cosmetric custom if there is
    strTitle = dicSetInfo[ "title" ]
    dicCustom = {} if strTitle not in dicSetDef[ "CosmetCustom" ] else dicSetDef[ "CosmetCustom" ][ strTitle ]
    
    # Decoration; more complicated than SIG because of custom
    strTitleR = dicSetInfo[ "title" ] if "title" not in dicCustom else dicCustom[ "title" ]
    nColor = dicSetInfo[ "colour" ] if "colour" not in dicCustom else dicCustom[ "colour" ]
    
    dicSet[ strKey ][ "hist" ].SetTitle(strTitleR)
    dicSet[ strKey ][ "hist" ].SetFillColor(nColor)
    dicSet[ strKey ][ "hist" ].SetLineColor(nColor)
    
    # Ordering (by title)!
    # Using a method giving an atttraction between samples with same title
    nIdx = 0
    for i in range(len(listMC)): 
      if listMC[ len(listMC) - 1 - i ].GetTitle() == dicSetInfo[ "title" ]: 
        nIdx = len(listMC) - 1 - i
        break
    
    listMC.insert(nIdx, dicSet[ strKey ][ "hist" ])
  elif dicSet[ strKey ][ "TYPE" ] == "RD": 
    # Well, no more than black dots for RD. Problem?
    histRD.Add(dicSet[ strKey ][ "hist" ])
    
    # No, I need to track the cutoff with respect to period
    for strRDRun in listRDRun: 
      if strRDRun in strKey: 
        dicRDRun[ strRDRun ][ "hist" ].Add(dicSet[ strKey ][ "hist" ])
        dicRDRun[ strRDRun ][ "step0" ] = dicSet[ strKey ][ "step0" ]
        dicRDRun[ strRDRun ][ "step1" ] = dicSet[ strKey ][ "step1" ]
        break

# Send the following samples under ground (in the stack)
listBack = ["QCD", "t#bar{t}+Jets"]
for strBack in listBack: 
  for i in range(len(listMC)): 
    if strBack == listMC[ i ].GetTitle(): listMC.insert(0, listMC.pop(i))

if dolog: listMC.reverse()

# Signal part must be on the top
listMC += listMCSig

# For only printing out the number of entries
if onlyPrintNum: 
  fNumEntries = 0
  strNameCurr = ""
  fNumTotal = 0
  
  for hCurr in listMC: 
    if hCurr.GetTitle() != strNameCurr: 
      if strNameCurr != "": 
        print "%s : %lf"%(strNameCurr, fNumEntries)
      
      fNumEntries = 0
      strNameCurr = hCurr.GetTitle()
    
    fNumEntries += hCurr.Integral(0, binning[ 0 ] + 1)
    fNumTotal += hCurr.Integral(0, binning[ 0 ] + 1)
  
  print "%s : %lf"%(strNameCurr, fNumEntries)
  print "Total MC : %lf"%(fNumTotal)
  
  listRDRun.sort()
  for strRDRun in listRDRun: 
    fNumEntries = dicRDRun[ strRDRun ][ "hist" ].Integral(0, binning[ 0 ] + 1)
    print "%s : %lf"%(strRDRun, fNumEntries)
  
  print "Total Data : %lf"%(histRD.Integral(0, binning[ 0 ]))
  
  sys.exit(0)

if fMinPlot is not None: histRD.SetMinimum(fMinPlot)
if fMaxPlot is not None: histRD.SetMaximum(fMaxPlot)

# Print it out!
strOutFile = "hist_%s_%s.png"%(strDirHist, strVarName)
print strOutFile, 

canvMain = drawTH1("asdf", CMS_lumi, listMC, histRD, x_name, y_name, dolog)
canvMain.SaveAs(strOutFile)


