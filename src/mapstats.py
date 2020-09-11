#!/usr/bin/env python3

""" mapstats.py
  Script to get statistics on maps of physical parameters

  O. Morata 2020

  Things to do:

   - output results:
     - x,y plots
     - KDEs

"""
import numpy as np
import numpy.ma as ma

from astropy.io import fits

import MapClass as maps

import matplotlib.pyplot as plt
from astropy.visualization import hist

import argparse as ap
import sys
import yaml

##-- Functions ---------------------------------------------------------

def read_command_line() :
    """ Read command line arguments."""
    
    parser = ap.ArgumentParser()

    parser.add_argument(
        '-c', dest='cfgfile',  help='configuration file', 
        default='config.yml', metavar='list')

    parser.add_argument(
        '-w', dest='wkdir',  help='working directory', 
        default='.', metavar='FILE')

    parser.add_argument(
        '-o', dest='odir',  help='output directory',
        default='.', metavar='FILE')

    
    return parser.parse_args()



def read_configfile(fname):
    """ Read the YAML configuration file."""
    
    print("  >> reading configuration file:", fname, "...")

    with open(fname, 'r') as ymlfile:
        try:
            cnfg = yaml.safe_load(ymlfile)
        except yaml.YAMLError as exc:
            print(exc)
        
    return cnfg



def read_clump_definition(infile):
    
    with fits.open(infile) as hdumask:
        clump_def = hdumask[0].data

    clump_idxs = clump_def.view(ma.MaskedArray)

    # mask out invalid values in clump mask and pixels not in any clump
    #
    clump_idxs_invalid = ma.masked_invalid(clump_idxs)
    incl = ma.masked_less(clump_idxs_invalid, 1)
    
    return incl



def plot_histo(data, cfg, outf):
    """Plot histogram of data

    Using configuration in cfg, and output to outf.
    It returns outbin for further processing.
    """

    meth = set_cnfgvalue(cfg, 'method', 'scott')
    httype = set_cnfgvalue(cfg, 'type', 'stepfilled')
    alpha = set_cnfgvalue(cfg, 'alpha', 0.8)
    dens = set_cnfgvalue(cfg, 'density', False)
    stack = set_cnfgvalue(cfg, 'stacked', True)
    cumul = set_cnfgvalue(cfg, 'cumulative', False)
    ylog = set_cnfgvalue(cfg, 'log', False)
    
    xlog = set_cnfgvalue(cfg, 'xlog', False)
    xlabel = set_cnfgvalue(cfg, 'xlabel', '')
    ylabel = set_cnfgvalue(cfg, 'ylabel', '')
    title = set_cnfgvalue(cfg, 'title', '')


    fig, ax = plt.subplots(1,1,figsize=(5,5))

    if xlog :
        plt.xscale('log')

    if xlabel:
        plt.xlabel(xlabel)
    if ylabel:
        plt.ylabel(ylabel)
    if title:
        plt.title(title)
        
    outbin = hist(data, bins=meth, ax=ax, histtype=httype,
             alpha=alpha, density=dens, range=(np.min(data),np.max(data)),
             stacked=stack, cumulative=cumul, log=ylog)
 
    fig.savefig(outf)

    return outbin


def set_cnfgvalue(cfg, tag, default):

    if tag in cfg:
        val = cfg[tag]
    else :
        val = default

    return val



def save_bins(binsfile, outb):
    """Save the histogram bins into binsfile"""
    
    sz = outb[0].size
    ab = np.zeros(sz, dtype=[('var1', float), ('var2', float)])
    ab['var1'] = outb[0]
    ab['var2'] = outb[1][:sz]
    
    np.savetxt(binsfile, ab, fmt='%7s   %s')



##-- End of functions --------------------------------------------------

args = read_command_line()


wdir = args.wkdir+'/'
outdir = args.odir+'/'


cnfg = read_configfile(args.cfgfile)



testmap = wdir+cnfg['infile']

varmap = maps.Map.from_fitsfile(testmap, name=testmap)

varmap = varmap.masked_invalid()
print("  >>> number of valid pixels in map:",varmap.count())


if 'clumpfile' in cnfg:
    cldefs = wdir+cnfg['clumpfile']

    print("  > reading clump definitions in", cldefs)
    inclumps = read_clump_definition(cldefs)

    if 'clumpcut' in cnfg:
        if len(cnfg['clumpcut']) % 2 != 0 :
            print("  ++ERROR: the filter in clumps does not have an even",
            "number of elements")
            sys.exit(1)
        
        clumpconds = cnfg['clumpcut']

    else:
        clumpconds = []
else:
    cldefs = None



if 'auxfiles' in cnfg:
    auxf = cnfg['auxfiles']
    numaux = len(auxf)

    if 'mask_aux' in cnfg :
        conds = cnfg['mask_aux']
        numcond = len(conds)

        if numcond != 3 * numaux:
            print("  ++ ERROR: wrong number of conditions for filters")
            sys.exit(1)
else:
    numaux  = 0


outfile = set_cnfgvalue(cnfg, 'outfile', 'out_histo.pdf')
outbins = set_cnfgvalue(cnfg, 'outbins', '')

# use data in clumps (if needed)
#
if cldefs:
    mapinclumps = varmap.masked_where(ma.getmask(inclumps))
    print("  >>> number of valid pixels in clumps:",mapinclumps.count())
    

    while len(clumpconds) > 0 :

        cond = clumpconds.pop(0)
        val = clumpconds.pop(0)

        if cond == ">" :
            mskcl = ma.masked_greater(inclumps, val)
        elif cond == "=" :
            mskcl = ma.masked_not_equal(inclumps, val)
        elif cond == "<" :
            mskcl = ma.masked_less(inclumps, val)
        else:
            print("  ++ wrong condition in filter clumps", cond)
            sys.exit(1)

        mapinclumps = mapinclumps.masked_where(ma.getmask(mskcl))
        print("  >>> number of valid pixels after clump filter:",
              mapinclumps.count())

else :
    mapinclumps = varmap.copy()


# filter using aux files
#
var_aux = mapinclumps.copy()

if numaux :

    for aux in auxf:

        dataset = conds.pop(0)
        cond = conds.pop(0)
        cut = conds.pop(0)
    
        auxfile = maps.Map.from_fitsfile(wdir+aux, name='aux file'+aux)

        if dataset == 'data':
            auxdata = auxfile.data[0]
        elif dataset == 'var' :
            auxdata = auxfile.data[1]
        else :
            print("  ++ ERROR: wrong definition of aux data set")
            sys.exit(1)

        if cond == ">" :
            mask_aux = ma.masked_greater(auxdata, cut)
        elif cond == "=" :
            mask_aux = ma.masked_values(auxdata, cut)
        elif cond == "<" :
            mask_aux = ma.masked_less(auxdata, cut)
        else:
            print("  ++ wrong condition in filter clumps", cond)
            sys.exit(1)

        var_aux = var_aux.masked_where(ma.getmask(mask_aux))
        print("  >>> number of valid pixels after filter aux:",
              var_aux.count())

        
if 'histo' in cnfg:
    hhh = var_aux.data[0]
    shdata = hhh[hhh.mask == False]

    outb = plot_histo(shdata, cnfg['histo'], outdir+outfile)

    if outbins :
        save_bins(outdir+outbins, outb)
 
else:
    print("\n  ++ histogram not defined")
    

print("....Done")

fig = plt.figure(figsize=(4,4))
plt.xscale('log')
plt.yscale('log')
plt.plot(shdata, shdata, '.', color='black')
fig.savefig("t.pdf")
