#!/usr/bin/env python 
# catGetDatasetInfo v7-4-4 # to make dataset lists
# sed -i 's/^\/store/root:\/\/cms-xrdr.sdfarm.kr:1094\/\/xrd\/store/g' *

import os, json, array, sys
import numpy as np
from math import ceil       
username = os.environ['USER']

analysis = 'top'

RunFiles = [
             'DoubleEG_Run2016B',
             'DoubleEG_Run2016Bv1',
             'DoubleEG_Run2016C',
             'DoubleEG_Run2016D',
             'DoubleEG_Run2016E',
             'DoubleEG_Run2016F',
             'DoubleEG_Run2016G',
             'DoubleEG_Run2016H',
             'MuonEG_Run2016B',
             'MuonEG_Run2016Bv1',
             'MuonEG_Run2016C',
             'MuonEG_Run2016D',
             'MuonEG_Run2016E',
             'MuonEG_Run2016F',
             'MuonEG_Run2016G',
             'MuonEG_Run2016H',
             'DoubleMuon_Run2016B',
             'DoubleMuon_Run2016Bv1',
             'DoubleMuon_Run2016C',
             'DoubleMuon_Run2016D',
             'DoubleMuon_Run2016E',
             'DoubleMuon_Run2016F',
             'DoubleMuon_Run2016G',
             'DoubleMuon_Run2016H',
             'SingleMuon_Run2016B', 
             'SingleMuon_Run2016Bv1',
             'SingleMuon_Run2016C',
             'SingleMuon_Run2016D',
             'SingleMuon_Run2016E',
             'SingleMuon_Run2016F',
             'SingleMuon_Run2016G',
             'SingleMuon_Run2016H',
             'SingleElectron_Run2016B',
             'SingleElectron_Run2016Bv1',
             'SingleElectron_Run2016C',
             'SingleElectron_Run2016D',
             'SingleElectron_Run2016E',
             'SingleElectron_Run2016F',
             'SingleElectron_Run2016G',
             'SingleElectron_Run2016H',
             'TT_powheg',
             'WJets',
             "SingleTop_tW",
             "SingleTbar_tW",
             'ZZ',
             'WW',
             'WZ',
             'DYJets',
             'DYJets_10to50',
            ]

datadir = '{}/src/nano/analysis/data/dataset/'.format(os.environ['CMSSW_BASE'])
#version = os.environ["CMSSW_VERSION"]


for i in RunFiles:
    datasetName = i
    fileList = datadir + 'dataset_' + datasetName + '.txt'
    jobName = analysis+'_'+datasetName 

    Dirname = "{}/src/nano/analysis/test/Batch/{}/".format(os.environ['CMSSW_BASE'],jobName)
    DirnameJDS = "{}/src/nano/analysis/test/Batch/".format(os.environ['CMSSW_BASE'])
    if os.path.isdir(Dirname):
        print "ERROR: output directory already existing."
        sys.exit()
    else: os.makedirs(Dirname)
    Dirname_ = "%s/src/nano/analysis/test/Results/%s/"%(os.environ['CMSSW_BASE'],datasetName)
    if not os.path.isdir(Dirname_):
        os.makedirs(Dirname_)

    Dirname_ = "%s/src/nano/analysis/test/Batch/NanoAOD/"%(os.environ['CMSSW_BASE'])
    if not os.path.isdir(Dirname_):
        os.makedirs(Dirname_)

    files = np.array([])
    for f in open(fileList).readlines():
        f = f.strip()
        f = f.strip('\',"')
        if len(f) < 5: continue
        if '#' == f[0] or '.root' != f[-5:]: continue
        files = np.append(files,[f])
    nFiles = len(files)     
    maxFiles = 10
    nSection = int(ceil(1.0*nFiles/maxFiles))
    count = 0
    for section in range(nSection):
        begin = section*maxFiles
        end = min(begin + maxFiles, nFiles)
        FileNames = files[begin:end]
        FileNamesStr = " ".join(str(i) for i in FileNames)

        print "@@ Writing run script..."
        jds = "%ssubmit.jds" %Dirname 
        fout = open(jds, "w")
        print>>fout, "# Job description file for condor job"
        print>>fout, """executable = {0}/src/nano/analysis/test/Batch/ttbar.sh
universe   = vanilla
log = condor.log
getenv     = True
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
output = {1}/job_{2}.log
error = {1}/job_{2}.err
transfer_input_files = NanoAOD
queue""" .format(os.environ['CMSSW_BASE'],Dirname, count)
        fout.close()
        count += 1 
        #jobName = analysis+'_'+datasetName
        subBatch = "condor_submit -batch-name {} -append 'arguments={} {}' {}".format(datasetName ,datasetName,FileNamesStr, jds)
        #print createbatch
        print subBatch 
            
        os.system(subBatch)
