#!/bin/bash
##
##  prepare_maps.sh
##
##  script to get SCUBA2 850 and 450 micron maps ready to be processed
##  by other tools:
##     - we strip the extra axis\
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
    echo " -c file   configuration file"
    echo " -d DIR    configuration directory (default ${CFG_DIR})"
    echo " -h        help"
    echo " -i DIR    input directory (default ${RES_DIR})"
    echo " -o DIR    ouput directory"
}

check_file(){
# check if file $1 of type $2 exists
#
    if [[ ! -f $1 ]];then
	echo -e "\n  ** ERROR: $2 file $1 not found\n"
	exit 1
    fi
}



build_arrays(){
    #
    pfx=$1
    InArray=("${!2}")

    if [[ ${#InArray[@]} -gt 0 ]];then
	arr=()
	for e in $InArray
	do
	    arr+=($pfx/$e)
	done
    else
	echo "ERROR: no defined input data "
	exit 1
    fi
}



convert_tofits(){
    sdffile=$1
    fitsfile=${sdffile/\.sdf/\.fits}
    
    echo " >>> converting $sdffile to $fitsfile";echo
    ${CONVERT_DIR}/ndf2fits $sdffile $fitsfile
}



strip_data(){
    # strips the third axis from the sdf files
    #
    argIn=("${!1}")
    argOut=("${!2}")

    nelem=${#argIn[@]}
    idx=0

    while [ $idx -lt $nelem ]
    do
	inf=${argIn[$idx]}
	echo;echo " Processing file $inf"
	echo "  stripping third axis...";echo
	
	dfile=${RES_DIR}/$inf
	check_file $dfile "data"
	
	ofile=${output_dir}/${argOut[$idx]}

	if [[ $do_sec -eq "y"  &&  $idx -gt 0 ]];then
	    ofile=${ofile/\.sdf/_strp\.sdf}
	fi

	${KAPPA_DIR}/ndfcopy in=\"$dfile\(:,:\)\" out=$ofile

	if [[ $idx -eq 0 && $do_fits -eq "y" ]];then
	   convert_tofits $ofile
	fi
	   
	idx=$((idx+1))

    done
}



align_maps(){
    # align maps to the reference, which is the first element of the
    # array
    #
    argIn=("${!1}")

    nelem=${#argIn[@]}

    idx=1

    ref_file=${output_dir}/${argIn[0]}
    while [ $idx -lt $nelem ]
    do
	ofile=${output_dir}/${argIn[idx]}
	ofile_strip=${ofile/\.sdf/_strp\.sdf}

	check_file $ofile_strip "data"

	echo;echo " >> Aligning $ofile_strip to $ofile"
	
	${KAPPA_DIR}/wcsalign in=$ofile_strip \
		    out=$ofile \
		    ref=$ref_file alignref=true

	if [[ $do_fits -eq "y" ]];then
	   convert_tofits $ofile
	fi
	
	rm $ofile_strip
	
	idx=$((idx+1))
    done
}




##-- End of functions --------------------------------------------------

# read cli options
#
while getopts "c:d:hi:o:" options
do
    case $options in
	c) cfg_file=$OPTARG
	   ;;
	d) CFG_DIR=$OPTARG
	   ;;
	h) help
	   exit 0
	   ;;
	i) RES_DIR=$OPTARG
	   ;;
	o) output_dir=$OPTARG
	   ;;
	*) echo " ERROR: wrong option"
	   exit 1
	   ;;
    esac
done


# check and fill variables
#
if [[ -n $cfg_file ]];then
    cfg_file=${CFG_DIR}/${cfg_file}
    check_file $cfg_file "configuration"

    source $cfg_file
else
    echo "ERROR: configuration file ${cfg_file} is not defined"
    exit 1
fi

if [[ ! -d ${output_dir} ]];then
    echo -e "\n  ** ERROR: output directory ${output_dir} not found\n"
    exit 1
fi



arr_data=($in_data)
arr_out=($out_data)

if [[ $do_snr -eq "y" ]];then
    arr_snr=($in_snr)
    
    arr_out_snr=($out_snr)
fi


# strip third axis from data files
#
strip_data arr_data[@] arr_out[@]

if [[ do_snr -eq "y" ]];then
    strip_data arr_snr[@] arr_out_snr[@]
fi


# align maps
#
if [[ $do_sec -eq "y" ]];then
    align_maps arr_out[@]

    if [[ do_snr -eq "y" ]];then
	align_maps arr_out_snr[@]
    fi
fi
