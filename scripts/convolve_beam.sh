#!/bin/bash

## convolve_beam
##
## O. Morata 2018-19
##
## script to convolve JCMT SCUBA2 observations with the beam of the
## alternative frequency
##
## it multiplies by given scale (to get the right values in (m)Jy/beam)
##

set -o errexit

if [[ $# -eq 4 ]]
then
    infile=$1
    psf=$2
    scale=$3
    outfile=$4
else
    echo;echo " ERROR in convolve_beam: wrong number of input parameters"
fi


tmpfile=${infile/.sdf/_mf.sdf}
tmpfile2=${infile/.sdf/_mf_scale.sdf}


echo " Convolving $infile with $psf..."
${KAPPA_DIR}/convolve in=$infile psf=$psf out=$tmpfile norm=true \
	    xcentre=0 ycentre=0

${KAPPA_DIR}/cmult in=$tmpfile out=$tmpfile2 scalar=$scale

rm $tmpfile
mv $tmpfile2 $outfile
