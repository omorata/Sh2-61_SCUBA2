#!/bin/bash

## mkratios.sh
##
## O. Morata 2019
##
## make the map of the flux ratios between the JCMT SCUBA2 450 and 850
## micron maps
##

set -o errexit

while getopts "c:d:" opcions
do
    case $opcions in
	c) CFG_FILE=$OPTARG
	   ;;
	d) OUT_DIR=$OPTARG
	   ;;
	*) echo "WARNING: unknown option"
	   ;;
    esac
done

if [[ -f ${CFG_FILE} ]];then
    source ${CFG_FILE}
else
    echo;echo "ERROR: ${CFG_FILE} not found"
    exit 1
fi


if [[ -n $file450 ]];then
    in450=${RES_DIR}/$file450
else
    echo;echo "ERROR: file450 not defined"
    exit 1
fi

if [[ -n $file850 ]];then
    in850=${RES_DIR}/$file850
else
    echo;echo "ERROR: file450 not defined"
    exit 1
fi

snr450=${cut_450:-0.}
snr850=${cut_850:-0.}

echo "SNR: $snr450, $snr850"

loth=${lowthr:-10000000}
hith=${hitrh:-10000000}


if [[ ! -d $OUT_DIR ]];then
    mkdir -p $OUT_DIR
fi
output=${OUT_DIR}/$outfile


# clip below some value of SNR
#
err450=${OUT_DIR}/${in450##*/}
err850=${OUT_DIR}/${in850##*/}

tmperr450=${err450/.sdf/_err.sdf}
tmperr850=${err850/.sdf/_err.sdf}

echo
echo "  Clipping data for 450 micron: SNR < $snr450..."
${KAPPA_DIR}/errclip $in450 $tmperr450 $snr450 snr

echo
echo "  Clipping data for 850 micron: SNR < $snr850..."
${KAPPA_DIR}/errclip $in850 $tmperr850 $snr850 snr


# filter out values out of the threshold values
#
if [[ $loth != "None" ]];then
    tmpthresh450=${err450/.sdf/_thresh.sdf}
    tmpthresh850=${err850/.sdf/_thresh.sdf}

    echo
    echo "  Applying threshold to 450micron data..."
    ${KAPPA_DIR}/thresh $tmperr450 $tmpthresh450 $loth $hith bad bad

    echo
    echo "  Applying threshold to 850micron data..."
    ${KAPPA_DIR}/thresh $tmperr850 $tmpthresh850 $loth $hith bad bad
fi


# make ratio
#
echo
echo " Calculating ratio map..."

${KAPPA_DIR}/div $tmpthresh450 $tmpthresh850 $output


rm $tmperr450 $tmperr850
rm -f $tmpthresh450 $tmpthresh850
echo "...Done"
