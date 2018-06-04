#!/usr/bin/env python
# source /cvmfs/cms.cern.ch/crab3/crab.sh

import os,json,sys,shutil,time,getopt

requestName = ""
submit = False
psetName =""

try:
    opts, args = getopt.getopt(sys.argv[1:],"hsi:n:p:j:b:",["requestName","psetName"])
except getopt.GetoptError:          
    print 'Usage : ./submitCrab3.py -n <requestName> -p <psetName>'
    sys.exit(2)

for opt, arg in opts:
    if opt == '-h':
        print 'Usage : ./submitCrab3.py -n <requestName> -p <psetName>'
        sys.exit()
    elif opt in ("-n", "--requestName"):
        requestName = arg
    elif opt in ("-s"):
        submit = True
    elif opt in ("-p", "--psetName"):
        psetName = arg

if requestName == "" :
    print "requestName(-n) is mandantory"
    sys.exit(-1)
if psetName == "" :
    print "psetName(-n) is mandantory"
    sys.exit(-1)

from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.transferLogs    = False
config.General.transferOutputs = True

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = psetName
#config.JobType.maxMemoryMB = 2500

config.section_("Data")
config.Data.publication  = False
config.Data.splitting='FileBased'
config.Data.unitsPerJob=1

config.section_("Site")
#crab checkwrite --site=T3_KR_KISTI --lfn=/store/group/nanoAOD/
config.Site.storageSite = 'T3_KR_KISTI'
#config.Site.storageSite = 'T3_KR_UOS'
config.Data.outLFNDirBase = '/store/group/hadAOD/%s/'%(requestName)

from CRABAPI.RawCommand import crabCommand

from findParents import findParent
datasets = json.load(open("%s/src/nano/nanoAOD/data/dataset/dataset.json"%os.environ['CMSSW_BASE']))
for d in datasets:

    if d['name'] != "TT_powheg": continue

    dataset = d['DataSetName']
    if len( dataset ) == 0: continue

    if 'had_path' in d and d['had_path']: continue

    doHadron = d['doHadron']
    if int(doHadron) == 0: continue
        
    isMC = True
    if d['type'] == 'Data':
        isMC = False

    # skip wrong input file
    if 'MC' in psetName and not isMC:
        continue
    if 'RD' in psetName and isMC:
        continue

    dataset = findParent(d).strip()
    if isMC :
        label = dataset.split("/")[1]
    else :
        label = dataset.split("/")[1]+"_"+dataset.split("/")[2]

    dataRequestName = '%s_%s'%(requestName,label)
    outputDatasetTag = dataset.split("/")[2]
    
    config.Data.inputDataset = dataset
    config.General.requestName = dataRequestName
    config.Data.outputDatasetTag = outputDatasetTag
    config.JobType.pyCfgParams = ['doHadron=%s'%(doHadron)]
    
    if 'Fast' in dataset:
        config.JobType.pyCfgParams.append('fastSim=True')


    if os.path.exists('crab_'+dataRequestName):
        continue
        
    print config,
    if submit:
        print 'submitting!'
        try:
            crabCommand('submit', config = config)
        except:
            print 'already done'
        #from multiprocessing import Process
        #p = Process(target=submit, args=(config,))
        #p.start()
        #p.join()
        time.sleep(10)

    print '^'*80


if not submit:
    print "Dry run, not submitting job and only printing crab3 command"
    print "Add -s to submit job"
    print 'Usage : ./submitCrab3Had.py -n <requestName> -p <inputFile> -s'
