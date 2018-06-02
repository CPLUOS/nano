#!/usr/bin/env python

import os
import pprint

pp = pprint.pprint
filename = 'dataset_hadron.txt'
f = open(filename , "w")
a = []
pathlist=  [#"/xrootd/store/group/nanoAOD/SingleMuon/run2_2016RD_NANO_Run2016B-18Apr2017_ver2-v1/171112_160845/0000/",
            # "/xrootd/store/group/nanoAOD/run2_2016v3/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180328_162849/0000/"
            "/xrootd/store/user/jlee/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/hadAOD/"
            ]   

for path in pathlist:
    temp = os.listdir(path)
    temp = map(lambda p: os.path.join(path, p), temp) 
    a +=  temp

for i in a:
    newstr = str.replace(i, '/xrootd','root://cms-xrdr.sdfarm.kr:1094///xrd')
    f.write(newstr + "\n")

pp(newstr)
f.close()    
