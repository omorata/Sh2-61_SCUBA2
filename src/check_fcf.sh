#!/bin/bash



picard () 
{ 
    ${ORAC_DIR}/etc/picard_start.sh ${1+"$@"}
}

set -e


while getopts "c:d:i:" options
do
    case $options in
	c) CONF_FILE=$OPTARG
	   ;;
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

LOGFILE=${RES_DIR}/checkcal.log
CHECK_DIR=${RES_DIR}/calcheck


source $CONF_FILE

touch $LOGFILE

if [[ -n $INDATA ]]
then
    infile=${RES_DIR}/$INDATA
    echo $infile
elif [[ -n ${CALIB_FILE} ]]
then
    infile=${RES_DIR}/${CALIB_FILE}
fi


if [[ ! -d ${CHECK_DIR} ]]
then
    mkdir ${CHECK_DIR}
fi

export ORAC_DATA_OUT=${CHECK_DIR}

picard -log sf -nodisplay SCUBA2_CHECK_CAL $infile | tee -a $LOGFILE
