#!/bin/bash
#

REDUCED_DIR=.
OUTFILE="out"
PARAMETER_FILE=""
RMS=n


set -e


##-- Functions ---------------------------------------------------------


help(){
    echo "Opcions:"
    echo "------------------------------------------------------"
    echo " -a action  action: add, cal, crop, mf, snr"
    echo " -d DIR     directory where reduced file is (default .)"
    echo " -h         help"
    echo " -i FILE    input file (without .sdf extension)"
    echo " -l FILE    list of files to process (for coadd)"
    echo " -o FILE    output file (without .sdf extension)"
    echo " -p FILE    parameter file"
    echo " -r         calculate rms of map (default:n)"
}



# redefine picard function
#
picard () 
{ 
    ${ORAC_DIR}/etc/picard_start.sh ${1+"$@"}
}




calibrate()
{
    DATA_FILE=${REDUCED_DIR}/$INFILE.sdf
    
    picard -log sf $OPT_PARAM CALIBRATE_SCUBA2_DATA $DATA_FILE |tee -a $LOGFILE

    mv ${REDUCED_DIR}/${INFILE}_uncal_cal.sdf ${REDUCED_DIR}/${INFILE}_cal.sdf 
}



# co-add maps
#
coadd()
{
    
    if [[ -z $LIST ]]
    then
	echo "  ERROR: missing listfile with names of files to coadd"
	return -1
    fi
    
    picard -log sf $OPT_PARAM MOSAIC_JCMT_IMAGES \
	   $(sed -e "s|^|${REDUCED_DIR}/|" $LIST) |tee -a $LOGFILE

    mv ${REDUCED_DIR}/*mos.sdf   ${REDUCED_DIR}/$OUTFILE.sdf
}



# crop map
#
crop(){
    DATA_FILE=${REDUCED_DIR}/$INFILE.sdf
    picard -log sf $OPT_PARAM CROP_SCUBA2_IMAGES $DATA_FILE |tee -a $LOGFILE


    OUTFILE=${REDUCED_DIR}/${INFILE}_crop.sdf
    # Calculate rms of the map.
    #
    if [[ $RMS == y ]]
    then
	${KAPPA_DIR}/stats $OUTFILE comp=err order |tee -a $LOGFILE
    fi
}



filter_match()
{
    echo "do something"
}



# make signal-to-noise map
#
make_snr_map()
{
    DATA_FILE=${REDUCED_DIR}/$INFILE.sdf
    OUTFILE=${REDUCED_DIR}/${INFILE}_snr.sdf

    ${KAPPA_DIR}/makesnr ${DATA_FILE} $OUTFILE |tee -a $LOGFILE
}


##-- End of functions --------------------------------------------------


while getopts "a:d:hi:l:o:p:r" options
do
    case $options in
	a) ACTION=$OPTARG
	   ;;
	d) REDUCED_DIR=$OPTARG
	   ;;
	h) help
	   exit 0
	   ;;
	i) INFILE=$OPTARG
	   ;;
	l) LIST=$OPTARG
	   ;;
	o) OUTFILE=$OPTARG
	   ;;
	p) PARAMETER_FILE=$OPTARG
	   ;;
	r) RMS=y
	   ;;
	*) echo "ERROR: wrong option"
	   exit 1
	   ;;
    esac
done




if [[ -z $PARAMETER_FILE && $ACTION != "snr" ]]
then
    echo "ERROR: parameter file  not defined"
    exit 2
fi

export ORAC_DATA_OUT=${REDUCED_DIR}

OPT_PARAM="-recpars $PARAMETER_FILE"
LOGFILE=${REDUCED_DIR}/process.log

touch $LOGFILE



if [[ -z $ACTION ]]
then
    echo "  ERROR: no action requested. Nothing to do"
    return 1
else
    case $ACTION in
	add) coadd
	     if [[ $? -eq -1 ]]
	     then
		 exit 1
	     fi
	     ;;
	cal) calibrate
	     ;;
	crop) crop
	      ;;
	mf) filter_match
	    ;;
	snr) make_snr_map
	     ;;
	*) echo "ERROR: unknown action"
	   exit 1
	
    esac
fi






#if [[ ${FILTER_MATCH} == "y" ]]
#then
#    echo do something
#fi



#convert to fits
#
#    ${CONVERT_DIR}/ndf2fits $OUTFILE.sdf $OUTFILE.fits |tee -a $LOGFILE

#fi

