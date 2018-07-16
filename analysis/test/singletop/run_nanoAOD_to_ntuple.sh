#!/bin/bash
jobName=st_nanoAOD_to_ntuple

#if [ $# != 1 ]; then
#    echo "JOB SECTION NUMBER IS MISSING!!!"
#    exit 1
#fi
SECTION=`printf %03d $1`

idxsec=$1
numPerJob=$2

if [ _$CMS_PATH == _ ]; then
  export CMS_PATH=/cvmfs/cms.cern.ch
  source $CMS_PATH/cmsset_default.sh
fi

hostname
tar xzf job.tar.gz
cd singletop_nanoAOD/src/nano/analysis/test/singletop


scram build ProjectRename
eval `scram runtime -sh`

echo BEGIN `date` cmsRun job_${SECTION}_cfg.py

ls -al

filename=`python batch_nanoAOD.py $idxsec $numPerJob filename`
samptype=`python batch_nanoAOD.py $idxsec $numPerJob samptype`
idxstart=`python batch_nanoAOD.py $idxsec $numPerJob idxstart`
idxend=`python batch_nanoAOD.py $idxsec $numPerJob idxend`

echo "JOB : $numPerJob ($numPerJob; $idxsec) $filename $samptype $idxstart $idxend"
time singletopAnalyser -q $filename $samptype $idxstart $idxend
EXITCODE=$?

if [ $EXITCODE == 0 ]; then hadd res.root out*.root ; fi

ls -al
if [ $EXITCODE == 0 ]; then
  echo "FILENAME `python batch_nanoAOD.py $idxsec $numPerJob sampname` res_`python batch_nanoAOD.py $idxsec $numPerJob sampleno`.root"
  echo FINISHED `date` cmsRun job_${SECTION}_cfg.py
else
  rm -f core.*
  echo TERMINATED_$EXITCODE `date` cmsRun job_${SECTION}_cfg.py
  exit 1
fi
echo ENDED `date`
