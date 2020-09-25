#!/bin/bash
##
##  catalog_to_polygonfile.sh
##  O. Morata 2020
##
##  script to transform the output the shape of SCUBA2 clumps out of
##  print_catalog to the format usable by dbxmap.py
##
##

BIN_DIR=${BIN_DIR:-./src}

if [[ $# -ne 2 ]];then
    echo "  >> ERROR: wrong number of parameters"
    echo -e "\n    Usage: catalog_to_polygonfile.sh infile outfile\n"
    exit 1
fi

catalog=$1
outname=$2


# print the ID and Shape fields of the findclumps catalog
#
${BIN_DIR}/print_catalog.py  -i $catalog  -o $outname  -t 'findclumps' \
	  -f "['PIDENT', 'Shape']"


# modify output file to the format needed by dbxmap
#
sed -i 's/PIDENT  Shape/id type    coord       corners/' $outname
sed -i 's/J2000 TOPOCENTER /"world_deg" "/g' $outname
sed -i 's/$/"/g' $outname
sed -i 's/corners"/corners/' $outname

echo -e "... done!\n"
