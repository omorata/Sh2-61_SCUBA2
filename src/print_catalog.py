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
        default=None, metavar='list')

    parser.add_argument(
        '-i', dest='infile',  help='input catalog', 
        default=None, metavar='FILE')

    parser.add_argument(
        '-o', dest='outfile',  help='output file',
        default='catalog.out', metavar='FILE')

    parser.add_argument(
        '-t', dest='ctype', help='type of output: all, phys, findclumps',
        default='phys', metavar='var'
        )
    
    return parser.parse_args()

##-- End of functions --------------------------------------------------

#


args = read_command_line()


infile = args.infile
if not infile :
    print("  >> ERROR: missing input file")
    sys.exit(1)
    

print("\n >> Reading catalog")
catlg = cl.ClumpCatalog.from_file([infile])

# field definitions
#
phys_outfields = catlg.clumps[0].phys_names
findclumps_outfields = catlg.clumps[0].findclumps_names


fout = args.outfile

out_type = args.ctype

if args.fields != None:
    fields= args.fields

else:

    if out_type == 'phys':
        fields = phys_outfields
    
    elif out_type == 'findclumps':
        fields = findclumps_outfields

    elif out_type == 'all' :
        fields =  findclumps_outfieds
        fields.extend(phys_outfields)

    else:
        print("  >> ERROR: unknown output type", out_type)
        sys.exit(1)


    

print(" >> Printing shapes to file:", fout,"\n")
catlg.print_catalog(ctype=out_type, fields=fields, filename=fout)

