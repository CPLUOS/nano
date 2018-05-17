#!/usr/bin/env python 
# catGetDatasetInfo v7-4-4 # to make dataset lists
# sed -i 's/^\/store/root:\/\/cms-xrdr.sdfarm.kr:1094\/\/xrd\/store/g' *

import os, json, array, sys
import numpy as np
from math import ceil       
    
    
analysis = 'copt'

pythonCfg = 'copt_training.py'

RunFiles = [
             "hadron",
           ]

maxFiles = 300
#SetDir = "test"
#datadir = '{}/src/nano/nanoAOD/data/dataset/dataset_'.format(os.environ['CMSSW_BASE'])
datadir = '{}/src/nano/analysis/test/topMass/cut/tmva/tmva_batch/dataset_'.format(os.environ['CMSSW_BASE'])
#print datadir

for datasetName in RunFiles:
    fileList = datadir + datasetName + '.txt'
    print fileList
    jobName = analysis+'_'+datasetName 

    Dirname = "{}/src/nano/analysis/test/topMass/cut/tmva/tmva_batch/training_file/{}".format(os.environ['CMSSW_BASE'],jobName)
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
        print>>fout, """executable = {0}/src/nano/analysis/test/topMass/cut/tmva/tmva_batch/copt.sh
universe   = vanilla

log = condor.log

getenv     = True
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
output = {1}/job_{2}.log
error = {1}/job_{2}.err
queue""" .format(os.environ['CMSSW_BASE'], Dirname, section)
        fout.close()

        subBatch = "condor_submit -batch-name %s -append 'arguments=%s %s' %s" %(datasetName ,datasetName,FileNamesStr, jds)
        #subBatch = "condor_submit -batch-name {} -append 'arguments={} {} {} {}' {}".format(datasetName, analyser, SetDir, datasetName, FileNamesStr ,jds)    
        print subBatch
        os.system(subBatch)

