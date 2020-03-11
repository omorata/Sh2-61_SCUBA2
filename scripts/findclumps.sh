#!/bin/bash
##
##  findclumps.sh
##
##  script to find clumps in a sdf file
##
##  Input parameters:
##
##   - configuration filename with the variables to run findclumps
##       the cfg_config variable should contain the name of another file
##       with the specific configuration parameter for the chosen clump
##       finding method
##   - data filename
##

# for makefile
#sh scripts/findclumps.sh -c fw-Sh2-61-j850r0_co_mb-snr6.sh -o z1 
#


## -c config file
## -f data_file
## -i input_directory
## -o output_directory
## -h help

CFG_DIR=${CFG_DIR:-.}
RES_DIR=${RES_DIR:-.}

set -o errexit


##-- Functions ---------------------------------------------------------

help(){
    echo "findclumps.sh"
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
	*) echo "ERROR: wrong option"
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
    echo "ERROR: configuration file is not defined"
    exit 1
fi


data_file=${RES_DIR}/${fclumps_infile}
check_file $data_file "data"


config_file=${CFG_DIR}/${fclumps_config}
check_file $config_file "configuration"

config_file="^"${config_file}


if [[ -d ${output_dir} ]];then
    out_file=${output_dir}/${fclumps_outfile}
    out_cat=${output_dir}/${fclumps_outcat}
    out_log=${output_dir}/${fclumps_logfile}
else
    echo -e "\n  ** ERROR: output directory ${output_dir} not found\n"
    exit 1
fi


# look for clumps
#
echo " >> Finding clumps...";echo
echo $test_config

${CUPID_DIR}/findclumps in=${data_file} \
	                out=${out_file} \
			outcat=${out_cat} \
			method=$fclumps_method \
			config=$config_file \
			deconv=$fclumps_deconv \
			logfile=${out_log} \
			wcspar=$fclumps_wcspar \
			shape=$fclumps_shape


# if we have used for example a SNR map to find clumps, extract the
# values from the original data 
#
if [[ ${fclumps_extract} -eq "y" ]];then
    echo;echo " >> Extracting clumps...";echo

    data_file=${RES_DIR}/${extract_data}

    check_file $data_file "data"
    
    out_extr_file=${output_dir}/${extract_out}
    out_extr_cat=${output_dir}/${extract_outcat}
    out_extr_log=${output_dir}/${extract_logfile}

    
    ${CUPID_DIR}/extractclumps data=$data_file \
		               deconv=$extract_deconv \
			       logfile=$out_extr_log \
			       mask=$out_file \
			       out=$out_extr_file \
			       outcat=$out_extr_cat \
			       shape=$extract_shape \
			       wcspar=$extract_wcspar
fi

echo " ... Done"
