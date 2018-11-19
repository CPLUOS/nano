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
                   'DYJets','DYJets_10to50',
                   'WJets', 'ZZ', 'WW', 'WZ',
                #   'TT_powheg_mtop1665', 'TT_powheg_mtop1695', 'TT_powheg_mtop1715',
                #   'TT_powheg_mtop1735', 'TT_powheg_mtop1755', 'TT_powheg_mtop1785',
                  ]

mcFiles_vts_dyjet = [
                     'DYJets',
                     'DYJets_10to50',
                     'DYJets_MG',
                     'DYJets_MG_10to50',
                    ]
mcFiles_vts_ttbar = [
                     'TT_powheg',
                     'TTJets_aMC',
                     'TTJets_DiLept_MG',
                     'TTJets_DiLept',
                     'TTJets_aMC_ext1', 
                     'TT_powheg_mtop1665',
                     'TT_powheg_mtop1695',
                     'TT_powheg_mtop1715',
                     'TT_powheg_mtop1735',
                     'TT_powheg_mtop1755',
                     'TT_powheg_mtop1785'
                    ]
mcFiles_vts_single = [
                      'SingleTop_tW',
                      'SingleTbar_tW',
                      'SingleTop_tW_noHadron',
                      'SingleTbar_tW_noHadron',
                      'SingleTop_t-channel',
                      'SingleTbar_t-channel',
                      'SingleTop_s-channel'
                     ]
mcFiles_vts_diboson = [
                       'ZZ',
                       'WW',     
                       'WZ', 
                      ]
mcFiles_vts_diboson_2 = [
                         'ZZTo4L_powheg',
                         'WWTo2L2Nu',     
                         'WZTo3LNu_amcatnlo', 
                         'WZTo2LQQ',
                         'ZZTo2L2Nu',
                         'ZZTo2L2Q'
                        ]
mcFiles_vts_triboson = [
                        'WWW',
                        'WWZ',     
                        'WZZ',
                        'ZZZ' 
                       ]
mcFiles_vts_wjets = [
                     'WJets',
                     'WJetsToLNu_HT-100To200_ext1',    
                     'WJetsToLNu_HT-200To400_ext1',    
                     'WJetsToLNu_HT-400To600_ext1',
                     'WJetsToLNu_HT-600To800_ext1',
                     'WJetsToLNu_HT-800To1200_ext1', 
                     'WJetsToLNu_HT-1200To2500_ext1',  
                     'WJetsToLNu_HT-2500ToInf_ext1',   
                     'WJetsToLNu_HT-70To100',
                     'WJetsToLNu_HT-100To200',         
                     'WJetsToLNu_HT-200To400',
                     'WJetsToLNu_HT-400To600',
                     'WJetsToLNu_HT-600To800',         
                     'WJetsToLNu_HT-800To1200',
                     'WJetsToLNu_HT-1200To2500',       
                     'WJetsToLNu_HT-2500ToInf'        
                    ]
mcFiles_vts = mcFiles_vts_dyjet + mcFiles_vts_ttbar + mcFiles_vts_single + mcFiles_vts_diboson + mcFiles_vts_diboson_2 + mcFiles_vts_triboson + mcFiles_vts_wjets

mcFiles_vts_2 = [
                  'GG_HToMuMu',
                  'VBF_HToMuMu',
                  'WMinusH_HToMuMu',
                  'ZH_HToMuMu',

                  'TTZToLLNuNu',
                  'TTWJetsToLNu',

                  'tZq_ToLL',

                  'DYToLL_2J',

                  'GluGlu_ZZTo4mu',
                  'GluGlu_ZZTo4e',
                  'GluGlu_ZZTo2e2mu',

                  'ttH_HToMuMu',

                  'tsW',

                  'VBF-C1N2_WZ',
                  'VBF-C1N2_leptonicDecays',
                  'VBF-C1N2_tauDecays',

                  'DYJetsToLL_M-50_HT-70to100',
                  'DYJetsToLL_M-50_HT-100to200',
                  'DYJetsToLL_M-50_HT-100to200_ext1',
                  'DYJetsToLL_M-50_HT-200to400',
                  'DYJetsToLL_M-50_HT-200to400_ext1',
                  'DYJetsToLL_M-50_HT-400to600',
                  'DYJetsToLL_M-50_HT-400to600_ext1',
                  'DYJetsToLL_M-50_HT-600to800',
                  'DYJetsToLL_M-50_HT-800to1200',
                  'DYJetsToLL_M-50_HT-1200to2500',
                  'DYJetsToLL_M-50_HT-2500toInf',

                  'EWKZ2Jets_ZToLL',

                  'EWKWPlus2Jets_WToLNu',
                  'EWKWMinus2Jets_WToLNu',

                  'QCD_HT50to100',
                  'QCD_HT100to200',
                  'QCD_HT200to300',
                  'QCD_HT300to500',
                  'QCD_HT500to700',
                  'QCD_HT700to1000',
                  'QCD_HT1000to1500',
                  'QCD_HT1500to2000',
                  'QCD_HT2000toInf',

                  'ZJetsTONuNu_HT-100To200',
                  'ZJetsTONuNu_HT-200To400',
                  'ZJetsTONuNu_HT-400To600',
                  'ZJetsTONuNu_HT-600To800',
                  'ZJetsTONuNu_HT-800To1200',
                  'ZJetsTONuNu_HT-1200To2500',
                  'ZJetsTONuNu_HT-2500ToInf',

                  'EWKZ2Jets_ZToNuNu',
                  'WWJJToLNuLNu_EWK',
                  'WZJJ_EWK',

                  'WToLNu_0J',
                  'WToLNu_1J',
                  'WToLNu_2J',

                  'QCD_Pt-15to20_MuEnriched',
                  'QCD_Pt-20to30_MuEnriched',
                  'QCD_Pt-30to50_MuEnriched',
                  'QCD_Pt-50to80_MuEnriched',
                  'QCD_Pt-80to120_MuEnriched',
                  'QCD_Pt-120to170_MuEnriched',
                  'QCD_Pt-170to300_MuEnriched',
                  'QCD_Pt-300to470_MuEnriched',
                  'QCD_Pt-470to600_MuEnriched',
                  'QCD_Pt-600to800_MuEnriched',
                  'QCD_Pt-800to1000_MuEnriched',

                  'QCD_Pt-20to30_EMEnriched',
                  'QCD_Pt-30to50_EMEnriched',
                  'QCD_Pt-50to80_EMEnriched',
                  'QCD_Pt-80to120_EMEnriched',
                  'QCD_Pt-120to170_EMEnriched',
                  'QCD_Pt-170to300_EMEnriched',
                  'QCD_Pt-300toInf_EMEnriched'
                ]
dataFiles = [
             'SingleMuon_Run2016', 'SingleElectron_Run2016',
             'DoubleMuon_Run2016', 'DoubleEG_Run2016',
             'MuonEG_Run2016'#, 'HTMHT_Run2016'
            ]
dataFiles = [data+period for period in ["B","C","D","E","F","G","H"] for data in dataFiles]

missing = [
           'DYJets_MG',
           'DYJets_MG_10to50',
           'TT_powheg',
           'TTJets_aMC',
           'TTJets_aMC_ext1',
           'SingleTop_t-channel',
           'DYToLL_2J',
           'QCD_HT100to200',
           'WToLNu_0J',
           'WToLNu_1J',
           'WToLNu_2J',
           'QCD_Pt-80to120_EMEnriched',
           'QCD_Pt-120to170_EMEnriched'
          ]

anaName = sys.argv[1]
analyser = anaName+'Analyser'
if   anaName == 'h2mu'    : RunFiles = mcFiles_h2mumu  + dataFiles
elif anaName == 'mass'    : RunFiles = mcFiles_topmass + dataFiles
elif anaName == 'vts'     : RunFiles = mcFiles_vts     + dataFiles
elif anaName == 'bTag'     : RunFiles = mcFiles_vts     + dataFiles
elif anaName == 'cutbased': RunFiles = mcFiles_topmass; analyser = "cutbased";
else: print "put right name of analysis (h2mu/mass/vts/cutbased)"
#RunFiles = ['WW']
#RunFiles = ['tt012j_bsbar_2l_FxFx_NANO', 'tt012j_bbbar_2l_FxFx_NANO', 'tt012j_bbars_2l_FxFx_NANO']
#RunFiles = [data+period for period in ["B","C","D","E","F","G","H"] for data in [ 'SingleMuon_Run2016', 'SingleElectron_Run2016', 'DoubleMuon_Run2016']]
#RunFiles = [data+period for period in ["B","C","D","E","F","G","H"] for data in [ 'SingleMuon_Run2016']]
#RunFiles = [data+period for period in ["B","C","D","E","F","G","H"] for data in [ 'SingleElectron_Run2016' ]]
#RunFiles = [data+period for period in ["B","C","D","E","F","G","H"] for data in [ 'DoubleMuon_Run2016' ]]
#RunFiles = [data+period for period in ["B","C","D","E","F","G","H"] for data in [ 'DoubleEG_Run2016' ]]
#RunFiles = [data+period for period in ["B","C","D","E","F","G","H"] for data in [ 'MuonEG_Run2016' ]]
#RunFiles = ['tt012j_bbars_2l_FxFx_NANO']
#RunFiles = ['tt012j_bsars_2l_FxFx_NANO']
#RunFiles = ['tt012j_bbbar_2l_FxFx_NANO']
#RunFiles = ['tt012j_bbars_2l_FxFx_NANO', 'tt012j_bsbar_2l_FxFx_NANO']
RunFiles = ['tt012j_bbars_2l_FxFx_NANO', 'tt012j_bsbar_2l_FxFx_NANO', 'tt012j_bbbar_2l_FxFx_NANO']
#RunFiles = ['tt012j_bbars_2l_FxFx_herwigpp_NANO', 'tt012j_bsbar_2l_FxFx_herwigpp_NANO']
RunFiles = ['tt012j_bbars_2l_FxFx_herwigpp_NANO', 'tt012j_bsbar_2l_FxFx_herwigpp_NANO', 'tt012j_bbars_2l_FxFx_NANO', 'tt012j_bsbar_2l_FxFx_NANO', 'tt012j_bbbar_2l_FxFx_NANO']
#RunFiles = mcFiles_vts_ttbar
#RunFiles = mcFiles_vts_wjets
RunFiles = mcFiles_vts
RunFiles = mcFiles_vts_2
RunFiles = dataFiles
#RunFiles = [ 'MuonEG_Run2016C' ]
#RunFiles = dataFiles + ["DYJets_MG_10to50"]
#RunFiles = mcFiles_vts+mcFiles_vts_2+dataFiles
#RunFiles = missing
#RunFiles = ['tsW']
RunFiles = ['tt012j_bbars_2l_FxFx_NANO_tt8888tt', 'tt012j_bsbar_2l_FxFx_NANO_tt8888tt', 'tt012j_bbars_2l_FxFx_NANO', 'tt012j_bsbar_2l_FxFx_NANO', 'tt012j_bbbar_2l_FxFx_NANO']
#RunFiles = ['tt012j_bbars_2l_FxFx_NANO', 'tt012j_bsbar_2l_FxFx_NANO', 'tt012j_bbbar_2l_FxFx_NANO']
#RunFiles = ['tt012j_bbars_2l_FxFx_NANO_tt8888tt', 'tt012j_bsbar_2l_FxFx_NANO_tt8888tt']

maxFiles = 2000
#maxFiles = 100
maxFiles = 50
#maxFiles = 30
maxFiles = 10
#maxFiles = 5
maxFiles = 1

setDir = "test"
#setDir = "run2_2016v5"
#setDir = "run2_2016v5_missing"
#setDir = "run2_2016v4"
#setDir = "tt012j_jetDau"
#setDir = "test2"

cmsswBase = os.environ['CMSSW_BASE']
datadir = '%s/src/nano/nanoAOD/data/dataset/dataset_'%cmsswBase

for datasetName in RunFiles:
    fileList = datadir + datasetName + '.txt'
    jobName = anaName+'_'+datasetName 

    Dirname = "%s/src/nano/analysis/test/Batch/%s"%(cmsswBase,jobName)
#    OutDir = "/xrootd/store/user/wjjang/Batch_log/%s"%(jobName)
    OutDir = "/cms/scratch/wjjang/Batch_log/%s"%(jobName)

    if os.path.isdir(OutDir):#os.path.isdir(Dirname):
        print "ERROR: output directory already existing."
        sys.exit()
    else: os.makedirs(OutDir)#os.makedirs(Dirname)

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
        jds = "%s/submit.jds" %OutDir#%Dirname 
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
output = /cms/scratch/wjjang/Batch_log/{2}/job_{3}.log
error = /cms/scratch/wjjang/Batch_log/{2}/job_{3}.err
queue""" .format(cmsswBase, analyser, jobName, section)#Dirname, section)
        fout.close()

        subBatch = "condor_submit -batch-name {} -append 'arguments={} {} {}' {}".format(datasetName, setDir, datasetName, FileNamesStr, jds)
        os.system(subBatch)


