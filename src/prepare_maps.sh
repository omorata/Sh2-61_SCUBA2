#!/bin/bash
##
##  prepare_maps.sh
##
##  script to get SCUBA2 850 and 450 micron maps ready to be processed
##  by other tools:
##     - we strip the extra axis
##     - we align the 450 micron map to the 850 micron map
##

CFG_DIR=${CFG_DIR:-.}
RES_DIR=${RES_DIR:-.}

output_dir=.

set -o errexit


##-- Functions ---------------------------------------------------------

help(){
    echo "prepare_maps.sh"
    echo " Options:"
    echo "--------------------------------------------------"
    echo " -f FILE   input file"
    echo " -h        help"
    echo " -o FILE   output file"
    echo " -r FILE   reference file for align"
    echo " -t task   task to execute (align, strip, tofits)"
}



convert_tofits(){
    # function to convert an sdf file into a fits file
    #
    
    sdffile=$1
    fitsfile=$2
    
    echo " >>> converting $sdffile to $fitsfile";echo
    str_fits=$fitsfile
    
    if [[ -f $fitsfile ]];then
	echo " *** $fitsfile already exists"
	echo "     Overwrite [x/n]"
	read overwrite

	if [[ $overwrite != "n" && $overwrite != "N" ]];then
	    str_fits="!$fitsfile"
	else
	    echo -e "\n  Nothing to do, then\n"
	    exit 0
	fi
    fi
    ${CONVERT_DIR}/ndf2fits $sdffile $str_fits
}



strip_data(){
    # strips the third axis from an sdf file
    #
    
    dfile=$1
    echo -e "\n Processing file $dfile"
    echo -e "  stripping third axis...\n"
	
    ofile=$2

    ${KAPPA_DIR}/ndfcopy in=\"$dfile\(:,:\)\" out=$ofile

}



align_maps(){
    # align maps to the reference, which is the first element of the
    # array
    #

    ifile=$1
    ref_file=$2
    ofile=$3
    echo;echo " >> Aligning $ifile to $ref"
    
    ${KAPPA_DIR}/wcsalign in=$ifile  out=$ofile  ref=$ref_file  \
		alignref=true accept
}

##-- End of functions --------------------------------------------------

# read cli options
#
while getopts "f:ho:r:t:" options
do
    case $options in
	f) infile=$OPTARG  ;;
	h) help
	   exit 0  ;;
	o) output_file=$OPTARG  ;;
	r) reference_file=$OPTARG  ;;
	t) task=$OPTARG  ;;
	*) echo " ERROR: wrong option"
	   help
	   exit 1  ;;
    esac
done


if [[ ! -f $infile ]];then
    echo -e "\n  ** ERROR: input file ${infile} not found\n"
    exit 1
fi

output_dir=$(dirname $infile)
if [[ ! -d ${output_dir} ]];then
    echo -e "\n  ** ERROR: output directory ${output_dir} not found\n"
    exit 1
fi

#if [[ -z $task ]]
#then
#    echo -e "\n  ** WARNING: task not defined; using tofits\n"
#    task="tofits"
#fi

case $task in
    "align") align_maps $infile $reference_file $output_file  ;;
    "strip") strip_data $infile $output_file  ;;
    "tofits") convert_tofits $infile  $output_file  ;;
    *) echo " wrong task option"
       exit 1
       ;;
esac

echo -e "\n ... done!"
