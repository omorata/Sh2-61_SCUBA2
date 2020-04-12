#!/usr/bin/env python3

import argparse as ap
import sys

import ClumpCatalog as cl


##-- Functions ---------------------------------------------------------

def read_command_line() :
    """ Read command line arguments."""
    
    parser = ap.ArgumentParser()

    parser.add_argument(
        '-f', dest='fields',  help='selected fields', 
        default=['PIDENT', 'Shape'], metavar='list')

    parser.add_argument(
        '-i', dest='infile',  help='input catalog', 
        default='cl1.fits', metavar='FILE')

    parser.add_argument(
        '-o', dest='outfile',  help='output file',
        default='shapes.out', metavar='FILE')

    return parser.parse_args()

##-- End of functions --------------------------------------------------

# default values
#
clfile = 'cl1.fits'
fout = 'shapes.out'
outfields = ['PIDENT', 'Shape']


args = read_command_line()

clfile = args.infile
fout = args.outfile
fields= args.fields


print("\n >> Reading catalog")
catlg = cl.ClumpCatalog.from_file([clfile])

print(" >> Printing shapes to file:", fout,"\n")
catlg.print_catalog(ctype='findclumps', fields=fields, filename=fout)

