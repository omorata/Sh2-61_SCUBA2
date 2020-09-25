#!/usr/bin/env python3
#
""""stats_utils.py
    Oscar Morata 2020

    Useful functions to process the stats scripts
"""

import argparse
import os


## Functions -----------------------------------------------------------

def read_command_line() :
    """ Read command line arguments."""
    
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-c', dest='cfgfile',  help='configuration file', 
        default='config.yml', metavar='var')

    parser.add_argument(
        '-w', dest='wkdir',  help='working directory', 
        default='.', metavar='DIR')

    parser.add_argument(
        '-o', dest='odir',  help='output directory',
        default='.', metavar='DIR')

    
    return parser.parse_args()



def check_dir(tdir):
    """Check that directory tdir exists """
    
    if os.path.isdir(tdir) :
        return(tdir+'/')
    else:
        print("  +++ ERROR:", tdir, "is not found")
        sys.exit(1)

## End of functions ----------------------------------------------------
