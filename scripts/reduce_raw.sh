#!/bin/bash
##
## reduc_script.sh
##  Oscar Morata 2018
##
##  Generic script to reduce the SCUBA2 observations in the Sh2-61 region
##
#

HOME_DIR=${HOME_DIR:-.}
RES_DIR=${RES_DIR:-results}
CFG_DIR=${CFG_DIR:-config}

#set -e


if [ $# -ne 1 ];then
    echo " ** ERROR:  name of configuration file is missing"
    exit 1
fi

inconf=$1
source $inconf


#echo;echo "  + file list:" $LIST_FILE
#echo "  + data directory:" $ORAC_DATA_IN
#echo "  + results directory:" $ORAC_DATA_OUT
#echo "  + reduction recipe:" $RECIPE
#echo "  + parameters:" $PARAMS;echo


case $FRQ in
    850) source ${ORAC_DIR}/etc/oracdr_scuba2_850.sh -cwd
	 ;;
    450) source ${ORAC_DIR}/etc/oracdr_scuba2_450.sh -cwd
	 ;;
    *) echo "ERROR: wrong value of FRQ"
       exit 1
       ;;
esac


source $inconf

ORAC_DATA_OUT=${RES_DIR}/${ORAC_DATA_OUT}
LIST_FILE=${HOME_DIR}/${LIST_FILE}


echo;echo "  + file list:" $LIST_FILE
echo "  + data directory:" $ORAC_DATA_IN
echo "  + results directory:" $ORAC_DATA_OUT
echo "  + reduction recipe:" $RECIPE
echo "  + parameters:" $PARAMS
echo


export ORAC_DATA_IN
export ORAC_DATA_OUT

if [[ ! -d ${ORAC_DATA_OUT} ]]
then
    mkdir -p ${ORAC_DATA_OUT}
fi

LOGFILE=${ORAC_DATA_OUT}/oracdr.log

oracdr  -files $LIST_FILE  -loop file  -nodisplay  \
	-log sf  -verbose $RECIPE  -recpars $PARAMS \
    |tee -a $LOGFILE



if [[ -n ${ORAC_DATA_OUT}/${REDUCED_FILE} && \
	  -n ${ORAC_DATA_OUT}/${OUTPUT_FILE} ]]
then
    echo "Copying to $OUTPUT_FILE"
    mv ${ORAC_DATA_OUT}/${REDUCED_FILE} ${ORAC_DATA_OUT}/${OUTPUT_FILE} |
       tee -a $LOGFILE
fi


reducbl=${ORAC_DATA_OUT}/reduced_blocks
if [[ ! -d  $reducbl ]]
then
    mkdir ${reducbl}
fi

echo;echo "   moving data to reduced_blocs..."
cd ${ORAC_DATA_OUT}
rm *256.png *64.png
mv s*.sdf s*.png  reduced_blocks/
mv log* reduced_blocks/
mv gs*.sdf gs*.png  gs*.FIT reduced_blocks/
echo "   ... done";echo

