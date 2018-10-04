#!/bin/bash
##
## reduc_script.sh
##
##  Oscar Morata 2018
##
##  Generic script to reduce the SCUBA2 observations in the Sh2-61 region
##
HOME_DIR=${HOME_DIR:-.}
RES_DIR=${RES_DIR:-results}
CFG_DIR=${CFG_DIR:-config}
BIN_DIR=${BIN_DIR:-scripts}

# flags for match filtering and reference reduction
#
MF=
REFERENCE=

set -e


if [ $# -ne 1 ];then
    echo " ** ERROR:  name of configuration file is missing"
    exit 1
fi

inconf=$1
source $inconf


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

export ORAC_DATA_IN
export ORAC_DATA_OUT

if [[ ! -d ${ORAC_DATA_OUT} ]]
then
    mkdir -p ${ORAC_DATA_OUT}
fi

LOGFILE=${ORAC_DATA_OUT}/oracdr.log

RESULT_FILE=${ORAC_DATA_OUT}/${OUTPUT_FILE}


echo
echo "  + file list:" $LIST_FILE |tee -a $LOGFILE
echo "  + data directory:" $ORAC_DATA_IN |tee -a $LOGFILE
echo "  + results directory:" $ORAC_DATA_OUT  |tee -a $LOGFILE
echo "  + reduction recipe:" $RECIPE  |tee -a $LOGFILE
echo "  + parameters:" $PARAMS  |tee -a $LOGFILE
echo


# if we have a reference reduced file, for future processing, we do not
# process the data again and copy the reference file as the reduced file
#
# this will be used to apply later a match filter with another beam
#
if [[ -z $REFERENCE ]]
then
    oracdr -files $LIST_FILE -loop file -nodisplay -log sf -verbose $RECIPE \
	   -recpars $PARAMS |tee -a $LOGFILE

    
    if [[ -n ${ORAC_DATA_OUT}/${REDUCED_FILE} && \
	      -n ${RESULT_FILE} ]]
    then
	echo "Copying to $OUTPUT_FILE"
	mv ${ORAC_DATA_OUT}/${REDUCED_FILE} ${RESULT_FILE} \
	    | tee -a $LOGFILE
    fi


    reducbl=${ORAC_DATA_OUT}/reduced_blocks
    if [[ ! -d  $reducbl ]]
    then
	mkdir ${reducbl}
    fi

    echo;echo "   moving data to reduced_blocks..."
    cd ${ORAC_DATA_OUT}
    rm *256.png *64.png
    mv s*.sdf s*.png  reduced_blocks/
    mv log* reduced_blocks/
    mv gs*.sdf gs*.png  gs*.FIT reduced_blocks/
    echo "   ... done";echo

else
    if [[ -f $REFERENCE ]]
    then
	echo;echo " Copying $REFERENCE to ${RESULT_FILE} ..."|tee -a $LOGFILE
	cp $REFERENCE ${RESULT_FILE}
    else
	echo;echo " ERROR: $REFERENCE could not be found"
	exit 2
    fi
fi


# To use match filtering, the run configuration file has to contain the
# variables:
#   MF=y
#   PSF_FILE=filename  (pointing to the file where the filter is)
#
#  and add the desired values to the params file
#
if [[ $MF == "y" ]]
then
    echo;echo;echo " ... applying match filter";echo |tee -a $LOGFILE
    
    OLDFILE=${RESULT_FILE/.sdf/-original.sdf}

    cp ${RESULT_FILE} $OLDFILE
    
    ${BIN_DIR}/convolve_beam.sh ${RESULT_FILE} ${PSF_FILE} ${RESULT_FILE} \
	| tee -a $LOGFILE

    if [[ -n ${CALIB_FILE} ]]
    then
	CALIB_REDUC=${ORAC_DATA_OUT}/${CALIB_FILE}
	${BIN_DIR}/convolve_beam.sh ${CALIB_REDUC} ${PSF_FILE} \
		  ${CALIB_REDUC} | tee -a $LOGFILE
    fi
    
fi
