#!/usr/bin/env python 

import os, json, sys
import numpy as np
from math import ceil

mcFiles_h2mumu = [
                  'ttH',
                  "WWTo2L2Nu", "WZTo3LNu_amcatnlo", "WZTo2LQQ", "ZZTo2L2Nu", "ZZTo2L2Q",
                  "TTZToLLNuNu", "ZZTo4L_powheg", "ttWToLNu",
                  "WWW", "WWZ", "WZZ", "ZZZ",
                  "ZZ", "WZ", "WW"
                 ]

mcFiles_topmass = [
                   'TT_powheg',
                #   'SingleTop_tW', 'SingleTbar_tW',
                #   'DYJets','DYJets_10to50',
                #   'WJets', 'ZZ', 'WW', 'WZ',
                #   'TT_powheg_mtop1665', 'TT_powheg_mtop1695', 'TT_powheg_mtop1715',
                #   'TT_powheg_mtop1735', 'TT_powheg_mtop1755', 'TT_powheg_mtop1785',
                  ]

mcFiles_vts = ['WW']

dataFiles = [
             'SingleMuon_Run2016', 'SingleEG_Run2016',
             'DoubleMuon_Run2016', 'DoubleEG_Run2016',
             'MuonEG_Run2016', 'MuonEG_Run2016',
             'DoubleMuon_Run2016', 'DoubleEG_Run2016']
dataFiles = [data+period for period in ["B","C","D","E","F","G","H"] for data in dataFiles]

anaName = sys.argv[1]
analyser = anaName+'Analyser'
if   anaName == 'h2mu'    : RunFiles = mcFiles_h2mumu  + dataFiles
elif anaName == 'mass'    : RunFiles = mcFiles_topmass + dataFiles
elif anaName == 'vts'     : RunFiles = mcFiles_vts     + dataFiles
elif anaName == 'cutbased': RunFiles = mcFiles_topmass; analyser = "cutbased";
else: print "put right name of analysis (h2mu/mass/vts/cutbased)"
RunFiles = ['WW']

maxFiles = 10
setDir = "test"
cmsswBase = os.environ['CMSSW_BASE']
datadir = '%s/src/nano/nanoAOD/data/dataset/dataset_'%cmsswBase

for datasetName in RunFiles:
    fileList = datadir + datasetName + '.txt'
    jobName = anaName+'_'+datasetName 

    Dirname = "%s/src/nano/analysis/test/Batch/%s"%(cmsswBase,jobName)
    if os.path.isdir(Dirname):
        print "ERROR: output directory already existing."
        sys.exit()
    else: os.makedirs(Dirname)

    files = np.array([])
    for f in open(fileList).readlines():
        f = f.strip()
        f = f.strip('\',"')
        if len(f) < 5: continue
        if '#' == f[0] or '.root' != f[-5:]: continue
        files = np.append(files,[f])
    nFiles = len(files) 
    nSection = int(ceil(1.0*nFiles/maxFiles))
    for section in range(nSection):
        begin = section*maxFiles
        end = min(begin + maxFiles, nFiles)
        FileNames = files[begin:end]
        FileNamesStr = " ".join(str(i) for i in FileNames)

        print "@@ Writing run script..."
        jds = "%s/submit.jds" %Dirname 
        fout = open(jds, "w")
        print>>fout, "# Job description file for condor job"
        print>>fout, """executable = {0}/bin/slc6_amd64_gcc630/{1}
universe   = vanilla

log = condor.log

getenv     = True
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
output = {2}/job_{3}.log
error = {2}/job_{3}.err
queue""" .format(cmsswBase, analyser, Dirname, section)
        fout.close()

        subBatch = "condor_submit -batch-name {} -append 'arguments={} {} {}' {}".format(datasetName, setDir, datasetName, FileNamesStr, jds)
        os.system(subBatch)


