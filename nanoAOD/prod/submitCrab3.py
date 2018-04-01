#!/usr/bin/env python
import os,json,sys,shutil,time,getopt

def submitjob(requestName, psetName, dataset, submit,lumiMask=None):
    print 'v'*80
    print "creating job"
    print dataset

    isMC = False
    if 'MC' in psetName:
        isMC = True
    
    dataset = dataset.strip()
    if isMC :
        label = dataset.split("/")[1]
    else :
        label = dataset.split("/")[1]+"_"+dataset.split("/")[2]

    dataRequestName = '%s_%s'%(requestName,label)
    outputDatasetTag = dataset.split("/")[2]
    outLFNDirBase = '/store/group/nanoAOD/%s/'%(requestName)
    
    sendjob = "crab submit JobType.psetName='%s' General.requestName='%s' Data.outLFNDirBase='%s' Data.outputDatasetTag='%s' Data.inputDataset='%s'"%(psetName,dataRequestName,outLFNDirBase,outputDatasetTag,dataset)
    if not isMC:
       sendjob+=" Data.splitting='LumiBased' Data.unitsPerJob=20 Data.lumiMask='%s'"%(lumiMask)
       
    print sendjob
    
    if submit:
        print "submiting job"
        os.system(sendjob)
        time.sleep(5)
    print '^'*80

    
submitBlock = None
requestName = ""
datasets = []
inputFile =None
submit = False
psetName =""
lumiMask =""

try:
    opts, args = getopt.getopt(sys.argv[1:],"hsi:n:p:j:b:",["requestName","inputFile","psetName","lumiMask","submitBlock"])
except getopt.GetoptError:          
    print 'Usage : ./submitCrab3.py -n <requestName> -i <inputFile> -p <psetName> -j <lumiMask>'
    sys.exit(2)

for opt, arg in opts:
    if opt == '-h':
        print 'Usage : ./submitCrab3.py -n <requestName> -i <inputFile> -p <psetName> -j <lumiMask>'
        sys.exit()
    elif opt in ("-n", "--requestName"):
        requestName = arg
    elif opt in ("-s"):
        submit = True
    elif opt in ("-i", "--inputFile"):
        inputFile = arg
        if os.path.isfile(inputFile):
            lines = open(inputFile)
            datasets = lines.readlines()
        else:
            datasets.append(inputFile)
    elif opt in ("-p", "--psetName"):
        psetName = arg
    elif opt in ("-b", "--submitBlock"):
        submitBlock = arg
    elif opt in ("-j", "--lumiMask"):
        lumiMask = arg

if requestName == "" :
    print "requestName(-n) is mandantory"
    sys.exit(-1)

if psetName == "" :
    print "psetName(-n) is mandantory"
    sys.exit(-1)
    
if inputFile is None:
    datasets = json.load(open("%s/src/nano/nanoAOD/data/dataset/dataset.json"%os.environ['CMSSW_BASE']))
    for d in datasets:
        dataset = d['DataSetName']
        if len( dataset ) == 0: continue

        if 'MC' in psetName:
            isMC = True
            if d['type'] == 'Data':
                continue
        else :
            if d['type'] != 'Data':
                continue            

        #if os.path.exists('crab_%s_%s' % (requestName, dataset.split('/')[1])): continue
        #if os.path.exists('crab_%s_%s_%s' % (requestName, dataset.split('/')[1], dataset.split('/')[2])): continue
        #if len( d['path']) == 0:
            #print d['path'], len( d['path'])
        submitjob(requestName, psetName, dataset, submit, lumiMask)
        
        #if submitBlock == '1' and 'QCD' in dataset:
        #    continue
        #if submitBlock == '2' and 'QCD' not in dataset:
        #    continue

else:
    for dataset in datasets:
        if len(dataset) < 10:
            continue
        if dataset.startswith("#"):
            continue
        submitjob(requestName, psetName, dataset, submit)

if not submit:
    print "Dry run, not submitting job and only printing crab3 command"
    print "Add -s to submit job"
    print 'Usage : ./submitCrab3.py -n <requestName> -i <inputFile> -s'
