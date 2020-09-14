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
import astropy

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



def read_plot_configuration(cfg):
    """Read the plot configuration"""
    
    plot_cfg = {}

    plot_cfg['type'] = set_cnfgvalue(cfg,'type', 'xyplot')
    plot_cfg['size'] = set_cnfgvalue(cfg,'size', [7,7])
    plot_cfg['xlog'] = set_cnfgvalue(cfg, 'xlog', False)
    plot_cfg['ylog'] = set_cnfgvalue(cfg, 'ylog', False)
    plot_cfg['xlabel'] = set_cnfgvalue(cfg, 'xlabel', '')
    plot_cfg['ylabel'] = set_cnfgvalue(cfg, 'ylabel', '')
    plot_cfg['title'] = set_cnfgvalue(cfg, 'title', '')
    plot_cfg['colors'] = set_cnfgvalue(cfg, 'colors', ['skyblue', 'red',
                                                      'green'])

    return plot_cfg


def read_histo_configuration(cfg):

    hcfg = {}
    hcfg['meth'] = set_cnfgvalue(cfg, 'method', 'scott')
    hcfg['httype'] = set_cnfgvalue(cfg, 'type', 'stepfilled')
    hcfg['alpha'] = set_cnfgvalue(cfg, 'alpha', 0.8)
    hcfg['dens'] = set_cnfgvalue(cfg, 'density', False)
    hcfg['stack'] = set_cnfgvalue(cfg, 'stacked', True)
    hcfg['cumul'] = set_cnfgvalue(cfg, 'cumulative', False)
    hcfg['ylog'] = set_cnfgvalue(cfg, 'log', False)

    return hcfg



def plot_histo(data, cfg, plcfg, outf):
    """Plot histogram of data

    Using configuration in cfg, and output to outf.
    It returns outbin for further processing.
    """

    histocfg = read_histo_configuration(cfg)

    fig, ax = plt.subplots(1,1,figsize=plcfg['size'])

    if plcfg['xlog'] :
        plt.xscale('log')

    if plcfg['xlabel']:
        plt.xlabel(plcfg['xlabel'])
    if plcfg['ylabel'] :
        plt.ylabel(plcfg['ylabel'])
    if plcfg['title']:
        plt.title(plcfg['title'])


    plmin = np.min(data[0][0])
    plmax = np.max(data[0][0])

    if histocfg['meth'] == 'scott':
        w, bins = astropy.stats.scott_bin_width(data[0][0], return_bins=True)

    elif histocfg['meth'] == 'freedman' :
        w, bins = astropy.stats.freedman_bin_width(data[0][0],
                                                   return_bins=True)
    elif histocfg['meth']== 'knuth' :
        w, bins = astropy.stats.knuth_bin_width(data[0][0].ravel(),
                                                return_bins=True, quiet=False)

    for ds in range(np.shape(data)[0]):
        shdata = data[ds][0]
        outbin = hist(shdata, bins=bins, ax=ax, histtype=histocfg['httype'],
                      alpha=histocfg['alpha'], density=histocfg['dens'],
                      range=(plmin, plmax),
                      stacked=histocfg['stack'], cumulative=histocfg['cumul'],
                      log=histocfg['ylog'],
                      color=plcfg['colors'][ds])
 
    fig.savefig(outf)

    return outbin



def set_cnfgvalue(cfg, tag, default):
    """ Read configuration value of tag, or use default value"""

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



def read_infile(cfg, wkdir) :
    """ Read the input files"""

    ifiles = len(cfg['infile'])

    map = np.empty(ifiles, dtype=maps.Map)

    for i in range(ifiles):

        testmap = wkdir+cfg['infile'][i]

        map[i] = maps.Map.from_fitsfile(testmap, name=testmap)

        map[i] = map[i].masked_invalid()
        #print("  >>> number of valid pixels in map:",i,"-",varmap[i].count())

    return map, ifiles



def read_dataset(cfg, wkdir):
    """ Read and filter a dataset"""
    
    if not 'infile' in cfg:
        print("  ++ ERROR: No input file. Nothing to do\n     Bye!")
        sys.exit(1)

    varmap, ifiles = read_infile(cfg, wkdir) 

    cldefs, clumpconds, inclumps = read_clump_conds(cfg, wkdir)

    auxfiles, auxcond = read_aux_files(cfg, wkdir)

    mapinclumps = filter_clumps(varmap, ifiles, cldefs, clumpconds, inclumps)

    var_aux = filter_auxfiles(mapinclumps, ifiles, auxfiles, auxcond)

    vdata = unmask(var_aux)


    return vdata



def unmask(var) :

    var_unmask = []
                          
    for v in range(len(var)):
        v_tmp  = var[v].data[0]
        v_tmp = v_tmp[v_tmp.mask == False]
        var_unmask.append(v_tmp)

    return var_unmask
        


def read_clump_conds(cfg, wkdir):
    """Read the filters for clumps"""
    
    if 'clumpfile' in cfg:

        defs = wkdir+cfg['clumpfile']

        print("  > reading clump definitions in", defs)
        inclumps = read_clump_definition(defs)

        if 'clumpcut' in cfg:
            if len(cfg['clumpcut']) % 2 != 0 :
                print("  ++ERROR: the filter in clumps does not have an even",
                      "number of elements")
                sys.exit(1)
        
            conds = cfg['clumpcut']

        else:
            conds = []
    else:
        defs = None
        conds = None
        inclumps = None

    return defs, conds, inclumps



def read_aux_files(cfg, wkdir):
    """Read the filters in the aux files"""
    
    if 'auxfiles' in cfg:

        auxf = cfg['auxfiles']
        numaux = len(auxf)
        
        for i in range(numaux):
            auxf[i] = wkdir+auxf[i]

        if 'mask_aux' in cfg :
            conds = cfg['mask_aux']
            numcond = len(conds)

            if numcond != 3 * numaux:
                print("  ++ ERROR: wrong number of conditions for filters")
                sys.exit(1)
    else:
        #numaux  = 0
        auxf = None
        numcond = None

    return auxf, conds



def filter_clumps(varmap, ifiles, defs, conds, inclumps):
    """Filter data in clumps"""
    
    mapinclumps = np.empty(ifiles, dtype=maps.Map)

    print("inclumps")
    if defs:
        print("defs")
        
        for i in range(ifiles) :
            mapinclumps[i] = varmap[i].masked_where(ma.getmask(inclumps))
            print("  >>> number of valid pixels in clumps:", i, "=>",
                  mapinclumps[i].count())


        while len(conds) > 0 :
                
            cond = conds.pop(0)
            val = conds.pop(0)

            if cond == ">" :
                mskcl = ma.masked_greater(inclumps, val)
            elif cond == "=" :
                mskcl = ma.masked_not_equal(inclumps, val)
            elif cond == "<" :
                mskcl = ma.masked_less(inclumps, val)
            else:
                print("  ++ wrong condition in filter clumps", cond)
                sys.exit(1)


        
            for i in range(ifiles) :
                mapinclumps[i] = mapinclumps[i].masked_where(ma.getmask(mskcl))
                print("  >>> number of valid pixels after clump filter:",
                      i, "=>", mapinclumps[i].count())

    else :
        for i in range(ifiles) :
            mapinclumps[i] = varmap[i].copy()

    return mapinclumps



def filter_auxfiles(mapinclumps, ifiles, auxf, conds):
    """Filter using aux files"""
    
    var_aux = np.empty(ifiles, dtype=maps.Map)
    for i in range(ifiles) :
        var_aux[i] = mapinclumps[i].copy()

    if len(auxf) :

        for aux in auxf:

            dataset = conds.pop(0)
            cond = conds.pop(0)
            cut = conds.pop(0)
    
            auxfile = maps.Map.from_fitsfile(aux, name='aux file'+aux)

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
                
            for i in range(ifiles) :
                var_aux[i] = var_aux[i].masked_where(ma.getmask(mask_aux))
                print("  >>> number of valid pixels after filter aux:", i, "=>",
                      var_aux[i].count())

    return var_aux

##-- End of functions --------------------------------------------------

args = read_command_line()


wdir = args.wkdir+'/'
outdir = args.odir+'/'


cnfg = read_configfile(args.cfgfile)


ds_str = [ k for k in cnfg if 'dataset_' in k ]

if ds_str :
    nsets = len(ds_str)
    if nsets == 1 :
        decl_str = "dataset"
    else :
        decl_str = "datasets"
        
    print("  >> found", nsets, decl_str, "defined")

    vdata = []
    for ds in ds_str:
        vdata.append(read_dataset(cnfg[ds], wdir))
        
    
else :
    print("  ++ ERROR: no dataset definition found")
    sys.exit(1)



#if not 'infile' in cnfg:
#    print("  ++ ERROR: No input file. Nothing to do\n     Bye!")
#    sys.exit(1)

#varmap, ifiles = read_infile(cnfg, wdir) 



#if 'clumpfile' in cnfg:
#
#    cldefs = wdir+cnfg['clumpfile']
#
#    print("  > reading clump definitions in", cldefs)
#    inclumps = read_clump_definition(cldefs)
#
#    if 'clumpcut' in cnfg:
#        if len(cnfg['clumpcut']) % 2 != 0 :
#            print("  ++ERROR: the filter in clumps does not have an even",
#            "number of elements")
#            sys.exit(1)
#        
#        clumpconds = cnfg['clumpcut']
#
#    else:
#        clumpconds = []
#else:
#    cldefs = None



#if 'auxfiles' in cnfg:
#    auxf = cnfg['auxfiles']
#    numaux = len(auxf)
#
#    if 'mask_aux' in cnfg :
#        conds = cnfg['mask_aux']
#        numcond = len(conds)
#
#        if numcond != 3 * numaux:
#            print("  ++ ERROR: wrong number of conditions for filters")
#            sys.exit(1)
#else:
#    numaux  = 0


outfile = set_cnfgvalue(cnfg, 'outfile', 'out_histo.pdf')
outbins = set_cnfgvalue(cnfg, 'outbins', '')

# use data in clumps (if needed)
#
#if cldefs:
#
#    mapinclumps = np.empty(ifiles, dtype=maps.Map)
#    for i in range(ifiles) :
#        mapinclumps[i] = varmap[i].masked_where(ma.getmask(inclumps))
#        print("  >>> number of valid pixels in clumps:", i, "=>",
#              mapinclumps[i].count())
#
#
#    while len(clumpconds) > 0 :
#
#        cond = clumpconds.pop(0)
#        val = clumpconds.pop(0)
#
#        if cond == ">" :
#            mskcl = ma.masked_greater(inclumps, val)
#        elif cond == "=" :
#            mskcl = ma.masked_not_equal(inclumps, val)
#        elif cond == "<" :
#            mskcl = ma.masked_less(inclumps, val)
#        else:
#            print("  ++ wrong condition in filter clumps", cond)
#            sys.exit(1)
#
#
#        
#        for i in range(ifiles) :
#            mapinclumps[i] = mapinclumps[i].masked_where(ma.getmask(mskcl))
#            print("  >>> number of valid pixels after clump filter:", i, "=>",
#                  mapinclumps[i].count())
#
#else :
#    mapinclumps = np.empty(ifiles, dtype=maps.Map)
#    for i in range(ifiles) :
#        mapinclumps[i] = varmap[i].copy()


# filter using aux files
#
#var_aux = np.empty(ifiles, dtype=maps.Map)
#for i in range(ifiles) :
#    var_aux[i] = mapinclumps[i].copy()
#
#if numaux :
#
#    for aux in auxf:
#
#        dataset = conds.pop(0)
#        cond = conds.pop(0)
#        cut = conds.pop(0)
#    
#        auxfile = maps.Map.from_fitsfile(wdir+aux, name='aux file'+aux)
#
#        if dataset == 'data':
#            auxdata = auxfile.data[0]
#        elif dataset == 'var' :
#            auxdata = auxfile.data[1]
#        else :
#            print("  ++ ERROR: wrong definition of aux data set")
#            sys.exit(1)
#
#        if cond == ">" :
#            mask_aux = ma.masked_greater(auxdata, cut)
#        elif cond == "=" :
#            mask_aux = ma.masked_values(auxdata, cut)
#        elif cond == "<" :
#            mask_aux = ma.masked_less(auxdata, cut)
#        else:
#            print("  ++ wrong condition in filter clumps", cond)
#            sys.exit(1)
#
#        for i in range(ifiles) :
#            var_aux[i] = var_aux[i].masked_where(ma.getmask(mask_aux))
#            print("  >>> number of valid pixels after filter aux:", i, "=>",
#                  var_aux[i].count())


#test = cnfg['test']
#print("test", test)

#for n in test:
#    print("k", n)
            
#sys.exit(0)

#print(np.shape(vdata))
#print(outfile)
if 'plot' in cnfg:
    plcf = read_plot_configuration(cnfg['plot'])


    
if plcf['type'] == 'histo' and 'histo' in cnfg:
    
    outb = plot_histo(vdata, cnfg['histo'], plcf, outdir+outfile)

    if outbins :
        save_bins(outdir+outbins, outb)
 
else:
    print("\n  ++ histogram not defined")
    

print("....Done")


#vvv = var_aux[1].data[0]
#vvdata = vvv[vvv.mask == False]

fig = plt.figure(figsize=(7,7))
plt.xscale('log')
#plt.yscale('log')
plt.plot(vdata[0][1], vdata[0][0], '.', color='skyblue')
plt.plot(vdata[1][1], vdata[1][0], '.', color='red')
plt.plot(vdata[2][1], vdata[2][0], '.', color='green')
fig.savefig("t.pdf")
sys.exit(0)


maxcl = 55
ct = 0
vmsk = ma.masked_greater(inclumps, maxcl)
arr_size = vmsk.count()
print("size max:", arr_size)
dist = np.empty((2,arr_size))
for cl in [0]:

    vmsk = ma.masked_not_equal(inclumps, cl+1)
    #print(vmsk.count())
    
    v = var_aux[0].masked_where(ma.getmask(vmsk))
    #print(v.count())
    zz = v.data[0]

    #print(zz.count())
    ind = np.unravel_index(np.argmax(v.data[0], axis=None), v.data[0].shape)
    #print(ind)
    a = ma.where(zz > 2)



    nn = np.shape(a)[1]

    tt = {}
    for i in range(nn):
        x = a[0][i]
        y = a[1][i]
        d = 3 * np.sqrt((ind[0]-x)*(ind[0]-x)+(ind[1]-y)*(ind[1]-y))
        dist[0][ct] = d
        dist[1][ct] = zz[x][y]
        #print(x, y, dist[0][ct], dist[1][ct])

        ct += 1

        k_d = round(d,4)*10000
        #k_d = int(d)
        if k_d in tt :
            tt[k_d].append(zz[x][y])
        else:
            tt[k_d] = [zz[x][y]]

ks = sorted(tt.keys())
npt = len(ks)
avg = np.zeros((3,npt))
ct =0
#print(ks)
for key in ks:
    narr = np.asarray(tt[key])
    #print(narr)
    avg[0][ct] = key /10000
    avg[1][ct] = np.average(narr)
    avg[2][ct] = np.std(narr)
    #print(avg[0][ct], avg[1][ct], avg[2][ct])
    ct +=1
    
    
fig = plt.figure(figsize=(7,7))
plt.xscale('log')
#plt.yscale('log')

plt.xlim(0.02,100)
#plt.ylim(5e21,1e24)
plt.plot(dist[0], dist[1], '.', color='black')
#plt.plot(avg[0], avg[1], '.', color='red')
plt.errorbar(avg[0], avg[1], avg[2], color='blue')
#plt.plot(avg[0], avg[1]+avg[2], color='blue')
#plt.plot(avg[0], avg[1]-avg[2], color='blue')
fig.savefig("tx.pdf")
