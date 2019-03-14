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

mcFiles_dm = [
                "VBF-C1N2_WZ",
                "VBF-C1N2_leptonicDecays",

                "WW", "ZZ", "WZ",
                "WWJJToLNuLNu_EWK",
                "WZJJ_EWK",
                "ZZJJTo4L_EWK",

                "SingleTop_t-channel",
                "SingleTbar_t-channel",
                "SingleTop_tW",
                "SingleTbar_tW",

                "QCD_HT50to100",
                "QCD_HT100to200",
                "QCD_HT200to300",
                "QCD_HT300to500",
                "QCD_HT500to700",
                "QCD_HT700to1000",
                "QCD_HT1000to1500",
                "QCD_HT1500to2000",
                "QCD_HT2000toInf",

                "DYJets_MG",
                "DYJets_MG_10to50",

                "TT_powheg",

                #"WToLNu_0J",
                #"WToLNu_1J",
                #"WToLNu_2J",
                
                "WJets",
                "WJetsToLNu_HT-70To100",
                "WJetsToLNu_HT-100To200",
                "WJetsToLNu_HT-100To200_ext1",
                "WJetsToLNu_HT-200To400",
                "WJetsToLNu_HT-200To400_ext1",
                "WJetsToLNu_HT-400To600",
                "WJetsToLNu_HT-400To600_ext1",
                "WJetsToLNu_HT-600To800",
                "WJetsToLNu_HT-600To800_ext1",
                "WJetsToLNu_HT-800To1200",
                "WJetsToLNu_HT-800To1200_ext1",
                "WJetsToLNu_HT-1200To2500",
                "WJetsToLNu_HT-1200To2500_ext1",
                "WJetsToLNu_HT-2500ToInf",
                "WJetsToLNu_HT-2500ToInf_ext1",
                "EWKZ2Jets_ZToLL",
                "EWKWPlus2Jets_WToLNu",
                "EWKWMinus2Jets_WToLNu",
                "EWKZ2Jets_ZToNuNu",

                "ZJetsTONuNu_HT-100To200",
                "ZJetsTONuNu_HT-200To400",
                "ZJetsTONuNu_HT-400To600",
                "ZJetsTONuNu_HT-600To800",
                "ZJetsTONuNu_HT-800To1200",
                "ZJetsTONuNu_HT-1200To2500",
                "ZJetsTONuNu_HT-2500ToInf",

                #'DYJetsToLL_M-50_HT-70to100',
                #'DYJetsToLL_M-50_HT-100to200',
                #'DYJetsToLL_M-50_HT-100to200_ext1',
                #'DYJetsToLL_M-50_HT-200to400',
                #'DYJetsToLL_M-50_HT-200to400_ext1',
                #'DYJetsToLL_M-50_HT-400to600',
                #'DYJetsToLL_M-50_HT-400to600_ext1',
                #'DYJetsToLL_M-50_HT-600to800',
                #'DYJetsToLL_M-50_HT-800to1200',
                #'DYJetsToLL_M-50_HT-1200to2500',
                #'DYJetsToLL_M-50_HT-2500toInf',

              ]


mcFiles_vts = ['WW']

#dataFiles = [
#             'SingleMuon_Run2016', 'SingleEG_Run2016',
#             'DoubleMuon_Run2016', 'DoubleEG_Run2016',
#             'MuonEG_Run2016', 'MuonEG_Run2016',
#             'DoubleMuon_Run2016', 'DoubleEG_Run2016']
dataFiles = [
             'HTMHT_Run2016', 
             #'MET_Run2016',  
             #'SingleMuon_Run2016'
            ]
dataFiles = [data+period for period in ["B","C","D","E","F","G","H"] for data in dataFiles]

anaName = sys.argv[1]
analyser = anaName+'Analyser'
if   anaName == 'h2mu'    : RunFiles = mcFiles_h2mumu  + dataFiles;
elif anaName == 'mass'    : RunFiles = mcFiles_topmass + dataFiles;
elif anaName == 'vts'     : RunFiles = mcFiles_vts     + dataFiles;
#elif anaName == 'cutbased': RunFiles = mcFiles_topmass; analyser = "cutbased";
elif anaName == 'dm'      : RunFiles = dataFiles; #+ mcFiles_dm; # + dataFiles;
else: print "put right name of analysis (h2mu/mass/vts/cutbased/dm/dmBg)"
#RunFiles = ['WW']

maxFiles = 20
setDir = "test190129"
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

requirements = ( HasSingularity == true ) 
accounting_group=group_cms 
+SingularityImage = "/cvmfs/singularity.opensciencegrid.org/opensciencegrid/osgvo-el6:latest" 
+SingularityBind = "/cvmfs, /cms, /share"

getenv     = True
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
output = {2}/job_{3}.log
error = {2}/job_{3}.err
queue""" .format(cmsswBase, analyser, Dirname, section)
        fout.close()
        #print setDir, datasetName, FileNamesStr
        subBatch = "condor_submit -batch-name {} -append 'arguments={} {} {}' {}".format(datasetName, setDir, datasetName, FileNamesStr, jds)
        os.system(subBatch)


