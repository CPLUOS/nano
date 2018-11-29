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
                  'SingleTop_tW', 'SingleTbar_tW',
                'SingleTop_s-channel',
                'SingleTop_t-channel',
                'SingleTbar_t-channel',
                  'DYJets','DYJets_10to50',
                  'ZZ', 'WW', 'WZ', 
                 'WToLNu_0J',
                 'WToLNu_1J',
                 'WToLNu_2J',
                "QCD_Pt-15to20_MuEnriched",
                "QCD_Pt-20to30_MuEnriched",
                "QCD_Pt-30to50_MuEnriched",
                "QCD_Pt-50to80_MuEnriched",
                "QCD_Pt-80to120_MuEnriched",
                "QCD_Pt-120to170_MuEnriched",
                "QCD_Pt-170to300_MuEnriched",
                "QCD_Pt-300to470_MuEnriched",
                "QCD_Pt-470to600_MuEnriched",
                "QCD_Pt-600to800_MuEnriched",
                "QCD_Pt-800to1000_MuEnriched",
                "QCD_Pt-1000toInf_MuEnriched",
                "QCD_Pt-20to30_EMEnriched",
                "QCD_Pt-30to50_EMEnriched",
                "QCD_Pt-50to80_EMEnriched",
                "QCD_Pt-80to120_EMEnriched",
                "QCD_Pt-120to170_EMEnriched",
                "QCD_Pt-170to300_EMEnriched",
                "QCD_Pt-300toInf_EMEnriched",
                #  "DYJetsToLL_M-50_HT-70to100",
                #  "DYJetsToLL_M-50_HT-100to200",
                #  "DYJetsToLL_M-50_HT-200to400",
                #  "DYJetsToLL_M-50_HT-400to600",
                #  "DYJetsToLL_M-50_HT-600to800",
                #  "DYJetsToLL_M-50_HT-800to1200",
                #  "DYJetsToLL_M-50_HT-1200to2500",
                #  "DYJetsToLL_M-50_HT-2500toInf",
                #   "WJets",
                #  "WJetsToLNu_HT-70To100",
                #    "WJetsToLNu_HT-100To200",
                #    "WJetsToLNu_HT-200To400",
                #    "WJetsToLNu_HT-400To600",
                #    "WJetsToLNu_HT-600To800",
                #    "WJetsToLNu_HT-800To1200",
                #    "WJetsToLNu_HT-1200To2500",
                #    "WJetsToLNu_HT-2500ToInf",
                  ]

#mcFiles_vts = ['WW']

dataFiles = [
             'SingleMuon_Run2016', 'SingleElectron_Run2016',
             #'DoubleMuon_Run2016', 'DoubleEG_Run2016',
             #'MuonEG_Run2016', 'MuonEG_Run2016',
             #'DoubleMuon_Run2016', 'DoubleEG_Run2016'
]
dataFiles = [data+period for period in ["B","C","D","E","F","G","H"] for data in dataFiles]

anaName = sys.argv[1]
analyser = anaName+'Analyser'
if   anaName == 'h2mu'    : RunFiles = mcFiles_h2mumu  + dataFiles
elif anaName == 'slmass'    : RunFiles = mcFiles_topmass + dataFiles
elif anaName == 'vts'     : RunFiles = mcFiles_vts     + dataFiles
elif anaName == 'cutbased': RunFiles = mcFiles_topmass; analyser = "cutbased";
else: print "put right name of analysis (h2mu/mass/vts/cutbased)"
#RunFiles = ['WW']

maxFiles = 20
setDir = "slmass_928re"
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

        subBatch = "condor_submit -batch-name {} -append 'arguments={} {} {}' {}".format(datasetName, setDir, datasetName, FileNamesStr, jds)
        os.system(subBatch)


