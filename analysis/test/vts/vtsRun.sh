#!/bin/bash 

#Temporary (not completed)

TT012J_BSBAR="/xrootd/store/user/iawatson/tt012j_bsbar_2l_FxFx"
TT012J_BBARS="/xrootd/store/user/iawatson/tt012j_bbars_2l_FxFx"


NANO=$TT012J_BBARS/NANOAOD
HAD=$TT012J_BBARS/HADAOD
HADTRUTH=$TT012J_BBARS/HADTRUTHAOD

RESULTPATH="root://cms-xrdr.sdfarm.kr:1094///xrd/store/user"

for f in $NANO/*.root; do
    if [[ $f -nt $HAD/`basename $f` ]]; then
        echo "There is no the corresponding file in HADAOD !" $f
    else
        echo "Processing " $f
        $SRT_CMSSW_BASE_SCRAMRTDEL/bin/slc6_amd64_gcc630/vtsAnalysis $f result/NANOHAD_`basename $f` $HAD/`basename $f`
#$RESULTPATH/${USER}/tt012j_bsbar_2l_FxFx/NANOHAD_`basename $f` $HAD/`basename $f`
    fi
done


