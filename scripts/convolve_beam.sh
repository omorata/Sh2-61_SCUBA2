#!/bin/bash

## convolve_beam
##
## O. Morata 2018
##
## script to convolve JCMT SCUBA2 observations with the beam of the
## alternative frequency
##

if [[ $# -eq 3 ]]
then
    infile=$1
    psf=$2
    outfile=$3
else
    echo;echo " ERROR in convolve_beam: wrong number of input parameters"
fi


tmpfile=${infile/.sdf/_mf.sdf}


echo " Convolving $infile with $psf..."
${KAPPA_DIR}/convolve in=$infile psf=$psf out=$tmpfile norm=true \
	    xcentre=0 ycentre=0


mv $tmpfile $outfile
