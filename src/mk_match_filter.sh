#!/bin/bash

## mk_match_filter
##
## O. Morata 2018
##
## script to convolve JCMT SCUBA2 observations with the beam of the
## alternative frequency

#-- Functions -------------------------------------------------------------

picard () 
{ 
    ${ORAC_DIR}/etc/picard_start.sh ${1+"$@"}
}

#-- End of functions ---------------------------------------------------


if [[ $# -eq 3 ]]
then
    params=$1
    infile=$2
    outfile=$3
else
    echo;echo " ERROR in mk_match_filter: wrong number of input parameters"
fi

OPT_PARAM="-recpars $params"

picard -log sf -nodisplay ${OPT_PARAM} SCUBA2_MATCHED_FILTER $infile

mffile=${infile/.sdf/_mf.sdf}

mv $mffile $outfile

