#!/usr/bin/python3


import ROOT
import os, sys, ctypes, array


strOut = sys.argv[ 1 ] # The output file name
strNameHist = sys.argv[ 2 ] # The name of histogram
nOffsetArg = 3

listSrc = []
listLumi = []

for i in range(( len(sys.argv) - nOffsetArg) // 2): 
  listSrc.append(sys.argv[ 2 * i + nOffsetArg ])
  listLumi.append(float(sys.argv[ 2 * i + nOffsetArg + 1 ]))

fInvSumLumi = 1.0 / sum(listLumi)
for i in range(len(listLumi)): listLumi[ i ] *= fInvSumLumi # What we actually want is relative lumis

# Extracting binning info
listBinX = []
listBinY = []

fFirst = ROOT.TFile.Open(listSrc[ 0 ])
hFirst = fFirst.Get(strNameHist)
strTitle = hFirst.GetTitle()

hBinX = hFirst.ProjectionX()
hBinY = hFirst.ProjectionY()

nNumBinX = hBinX.GetNbinsX()
nNumBinY = hBinY.GetNbinsX()

for i in range(nNumBinX + 1): listBinX.append(hBinX.GetBinLowEdge(i + 1))
for i in range(nNumBinY + 1): listBinY.append(hBinY.GetBinLowEdge(i + 1))

fFirst.Close()

# Preparing the output file and directory
fNew = ROOT.TFile.Open(strOut, "RECREATE")

# Some of source contains histogram into a directory. It will be also duplicated
strDir = "/".join(strNameHist.split("/")[ : -1 ])
if strDir != "": 
  fNew.mkdir(strDir)
  fNew.cd(strDir)

# There will be data
hNew = ROOT.TH2F(strNameHist.split("/")[ -1 ], strTitle, 
  nNumBinX, array.array("d", listBinX), nNumBinY, array.array("d", listBinY))

# The main part: Reading all source and combining the info from them
# I used the following fomula: 
#   SF  = sum_i ( SF_i * (rel. lumi)_i )
#   err = sqrt( sum_i ( err_i * (rel. lumi)_i )^2 )
# First, SF and err^2 will be calculated in stored in hNew. After that, sqrt(err^2) will be stored.
fSrc = None
strPathCurr = ""
class NonConsistenceException(Exception): pass

try: 
  for nIdxSrc in range(len(listSrc)): 
    strPathCurr = listSrc[ nIdxSrc ]
    fSrc = ROOT.TFile.Open(strPathCurr)
    hSrc = fSrc.Get(strNameHist)
    
    if hSrc == None: raise NonConsistenceException
    
    # Btw, for preventing a mistake, I added a check whether histograms are consistence or not
    # The rule is easy: Does these histograms have same binning?
    if nIdxSrc != 0: # The reference is the first histo, where we already got the bin infos of this
      hBinXCurr = hSrc.ProjectionX()
      hBinYCurr = hSrc.ProjectionY()
      
      if nNumBinX != hBinXCurr.GetNbinsX() or nNumBinY != hBinYCurr.GetNbinsX(): 
        raise NonConsistenceException
      
      for i in range(nNumBinX + 1): 
        if listBinX[ i ] != hBinXCurr.GetBinLowEdge(i + 1): 
          raise NonConsistenceException
      
      for i in range(nNumBinY + 1): 
        if listBinY[ i ] != hBinYCurr.GetBinLowEdge(i + 1): 
          raise NonConsistenceException
    
    # Now, no problem. Time to calculate!
    for i in range(nNumBinX): 
      for j in range(nNumBinY): 
        fCX = 0.5 * ( listBinX[ i ] + listBinX[ i + 1 ] )
        fCY = 0.5 * ( listBinY[ j ] + listBinY[ j + 1 ] )
        
        nIdxXPre = ctypes.c_int()
        nIdxYPre = ctypes.c_int()
        nIdxZPre = ctypes.c_int()
        
        # There is a strange issue that the bin indices don't correspond to the displayed bins, 
        # which means, e.g., if I call hSrc.GetBinContent(1, 3) 
        # (see HLT_Ele32_eta2p1_WPTight_Gsf_FullRunRange.root)
        # then this method returns the content in (1, 1)-bin.
        # Even this phenomonon does not occur always; only some of histograms has this issue.
        # This complicated way to get the bin indices is for avoiding this problem.
        hSrc.GetBinXYZ(hSrc.FindBin(fCX, fCY), nIdxXPre, nIdxYPre, nIdxZPre)
        nIdxX = nIdxXPre.value
        nIdxY = nIdxYPre.value
        
        fSF = hNew.GetBinContent(i + 1, j + 1)
        fErr = hNew.GetBinError(i + 1, j + 1)
        
        fSF  += hSrc.GetBinContent(nIdxX, nIdxY) * listLumi[ nIdxSrc ]
        fErr += ( hSrc.GetBinError(nIdxX, nIdxY) * listLumi[ nIdxSrc ] ) ** 2
        
        hNew.SetBinContent(i + 1, j + 1, fSF)
        hNew.SetBinError(i + 1, j + 1, fErr)
    
    fSrc.Close()
  
except NonConsistenceException:
  sys.stderr.write("FATAL ERROR: Non-consistence of histograms; %s\n"%strPathCurr)
  
  fSrc.Close() # In here, one of source file had been opened but not closed
  fNew.Close()
  
  sys.exit(1)

# As mentioned, storing sqrt(err^2)
for nIdxX in range(nNumBinX): 
  for nIdxY in range(nNumBinY): 
    hNew.SetBinError(nIdxX + 1, nIdxY + 1, ( hNew.GetBinError(nIdxX + 1, nIdxY + 1) ) ** 0.5)

fNew.Write()
fNew.Close()


