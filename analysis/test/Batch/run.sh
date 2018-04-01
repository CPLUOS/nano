#!/bin/bash 

#RESULT_DIR="NanoAOD"
#DIR=`hostname`

#/bin/mkdir -p $RESULT_DIR/$DIR
#echo `hostname` >> $RESULT_DIR/$DIR/`hostname`.txt
#echo `date` >> $RESULT_DIR/$DIR/`hostname`.txt
#source /cvmfs/cms.cern.ch/cmsset_default.sh >> $RESULT_DIR/$DIR/`hostname`.txt
$SRT_CMSSW_BASE_SCRAMRTDEL/bin/slc6_amd64_gcc630/$1 ${@:2}
#echo "Done." >> $RESULT_DIR/$DIR/`hostname`.txt

