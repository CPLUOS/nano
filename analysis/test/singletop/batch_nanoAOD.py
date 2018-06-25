import sys, os


def myGetNumList(strPath):
  listList = open(strPath).readlines()
  
  nNumRes = 0
  for strLine in listList: 
    if strLine[ 0 ] == "#": continue
    if strLine.strip() == "": continue
    nNumRes += 1
  
  return nNumRes

listMC = []
listRD = []
for s in open("rootlist/datasets_MC.txt").readlines(): listMC.append(s.replace("\n", ""))
for s in open("rootlist/datasets_RD.txt").readlines(): listRD.append(s.replace("\n", ""))

if len(sys.argv) < 3: sys.exit(-1)
if sys.argv[ 1 ] != "n" and len(sys.argv) < 4: sys.exit(-1)

nIdxJob = -1 if sys.argv[ 1 ] == "n" else int(sys.argv[ 1 ])
nNumPerJob = int(sys.argv[ 2 ])
strType = "" if sys.argv[ 1 ] == "n" else sys.argv[ 3 ]

strPathDataset = "../../../nanoAOD/data/dataset"
listDataset = os.listdir(strPathDataset)

listAll = listMC + listRD

i = 0
nNumPassJobs = 0

while True: 
  strTypeMCRD = "MC" if i < len(listMC) else "RD"
  
  strPathList = os.path.join(strPathDataset, "dataset_%s.txt"%listAll[ i ])
  nLenList = myGetNumList(strPathList)
  
  nNextNumPassJobs = nNumPassJobs + nLenList / nNumPerJob
  if nLenList % nNumPerJob != 0: nNextNumPassJobs += 1
  
  if nNumPassJobs <= nIdxJob and nIdxJob < nNextNumPassJobs:
    strOut = ""
    if strType == "sampname": strOut = listAll[ i ]
    if strType == "sampleno": strOut = "%i"%(nIdxJob - nNumPassJobs)
    if strType == "filename": strOut = strPathList
    if strType == "samptype": strOut = strTypeMCRD
    if strType == "idxstart": strOut = "%i"%(( nIdxJob - nNumPassJobs ) * nNumPerJob)
    if strType == "idxend":   strOut = "%i"%(( nIdxJob - nNumPassJobs + 1 ) * nNumPerJob)
    
    sys.stdout.write(strOut)
    break
  
  nNumPassJobs = nNextNumPassJobs;
  
  i += 1
  if i >= len(listAll): break

if nIdxJob == -1:
  strOut = "%i,%i"%(nNumPerJob, nNumPassJobs)
  sys.stdout.write(strOut)

sys.stdout.flush()


