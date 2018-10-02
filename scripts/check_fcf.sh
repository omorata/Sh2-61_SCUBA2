#!/bin/bash

picard () 
{ 
    ${ORAC_DIR}/etc/picard_start.sh ${1+"$@"}
}

set -e

LOGFILE=${RES_DIR}/checkcal.log


while getopts "d:i:" options
do
    case $options in
	d) RES_DIR=$OPTARG
	   ;;
	h) help
	   return 0
	   ;;
	i) INDATA=$OPTARG
	   ;;
	*) echo "ERROR: wrong option"
	   exit 1
	   ;;
    esac
done

touch $LOGFILE

infile=${RES_DIR}/$INDATA
echo $infile

export ORAC_DATA_OUT=${RES_DIR}

picard -log sf -nodisplay SCUBA2_CHECK_CAL $infile | tee -a $LOGFILE
