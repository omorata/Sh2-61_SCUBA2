#!/bin/bash
##
##  polygons_to_plot.sh
##  O. Morata 2020
##
##  script to transform the output the shape of SCUBA2 clumps out of
##  print_catalog to the format usable by dbxmap.py
##
##

if [[ $# -ne 2 ]];then
    echo "  >> ERROR: wrong number of parameters"
    echo -e "\n    Usage: polygons_to_plot.sh infile outfile\n"
    exit 1
fi

catalog=$1
outname=$2

cp $catalog $outname

sed -i 's/PIDENT  Shape/id type    coord       corners/' $outname
sed -i 's/J2000 TOPOCENTER /"world_deg" "/g' $outname
sed -i 's/$/"/g' $outname
sed -i 's/corners"/corners/' $outname

echo -e "... done!\n"
