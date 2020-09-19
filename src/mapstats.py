#!/usr/bin/env python3

""" mapstats.py
  Script to get statistics on maps of physical parameters

  O. Morata 2020

  Things to do:

   - output results:
     - KDEs

"""

import numpy as np
import numpy.ma as ma

from astropy.io import fits

import argparse
import os
import sys
import yaml

import MapClass as maps
import varplot 
import matplotlib.pyplot as plt

##-- Functions ---------------------------------------------------------

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

    #vdata = unmask(var_aux)

    #return vdata, inclumps

    return var_aux, inclumps


def unmask(var) :
    """Unmask map data"""
    
    new_vdata = []
    for ds in range(len(var)):

        new_fdata = []
        fmap = var[ds]

        for v in range(len(fmap)):
            v_tmp  = fmap[v].data[0]
            v_tmp = v_tmp[v_tmp.mask == False]
            new_fdata.append(v_tmp)

        new_vdata.append(new_fdata)

    return new_vdata
        


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
        conds = None

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

    if auxf :

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



def plot_vsdist(data, cfg, inclumps):


    maxcl = varplot.set_cnfgvalue(cfg, 'max_clumps', 100)
    pixsize = varplot.set_cnfgvalue(cfg, 'pixel_size', 1.)
    
    vmsk = ma.masked_greater(inclumps, maxcl)
    arr_size = vmsk.count()
    #print("size max:", arr_size)

    print(np.shape(data))
    fig = plt.figure(figsize=(7,7))
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.9,150)
    plt.ylim(2e21,1e24)

    
    for ds in range(np.shape(data)[0]) :
        varpl = data[ds][0]

        dist = np.empty((2,arr_size))
        ct = 0

        print("DS:", ds)
        for cl in range(maxcl):

            vmsk = ma.masked_not_equal(inclumps, cl+1)
        
            v_inclump = varpl.masked_where(ma.getmask(vmsk))
        
            data_v = v_inclump.data[0]
            #print(data_v)
            #print(zz.count())
            maxpos = np.unravel_index(np.argmax(data_v, axis=None),
                                   data_v.shape)

            valpos = ma.where(data_v)

            n_pix = np.shape(valpos)[1]
            print("npix", n_pix)
            if n_pix == 0 :
                continue

            dbins = {}
            for i in range(n_pix):
                x = valpos[0][i]
                y = valpos[1][i]
                dx = maxpos[0]-x
                dy = maxpos[1]-y
                dist[0][ct] = pixsize * np.sqrt(dx * dx + dy * dy)
                dist[1][ct] = data_v[x][y]


                k_d = int(dist[0][ct]/6) * 6
                #print(k_d)
                #k_d = int(dist[0][ct])
                if k_d in dbins :
                    dbins[k_d].append(data_v[x][y])
                else:
                    dbins[k_d] = [data_v[x][y]]

                ct += 1
                
            ks = sorted(dbins.keys())
            npt = len(ks)
            avg = np.zeros((3,npt))
            cpt =0
            #print(ks)
            for key in ks:
                narr = np.asarray(dbins[key])
                #print(narr)
                avg[0][cpt] = key 
                avg[1][cpt] = np.average(narr)
                avg[2][cpt] = np.std(narr)
                #print(avg[0][ct], avg[1][ct], avg[2][ct])
                cpt +=1
    
    
            plt.plot(dist[0], dist[1], '.', color='black')
            #plt.plot(avg[0], avg[1], '.', color='red')
            #plt.errorbar(avg[0], avg[1], avg[2], color='blue')
    #plt.plot(avg[0], avg[1]+avg[2], color='blue')
    #plt.plot(avg[0], avg[1]-avg[2], color='blue')
    outfile = varplot.set_cnfgvalue(cfg, 'outfile', 'plotvsdist.pf')
    fig.savefig(outfile)


##-- End of functions --------------------------------------------------

args = read_command_line()


wdir = check_dir(args.wkdir)
outdir = check_dir(args.odir)


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
        #vdata.append(read_dataset(cnfg[ds], wdir))
        vd, inclumps = read_dataset(cnfg[ds], wdir)
        vdata.append(vd)
    
else :
    print("  ++ ERROR: no dataset definition found")
    sys.exit(1)


#outbins = set_cnfgvalue(cnfg, 'outbins', '')
outbins = ''

if 'plot' in cnfg:
    plcf = varplot.read_plot_configuration(cnfg['plot'])

    
    if plcf['type'] == 'histo':
        if 'histo' in cnfg:
    
            outb = varplot.plot_histo(unmask(vdata), cnfg['histo'], plcf,
                                      outdir)

            if outbins :
                save_bins(outdir+outbins, outb)
 
        else:
            print("\n  ++ ERROR: histogram not defined")
            sys.exit(1)

            
    elif plcf['type'] == 'xyplot':

        if 'xyplot' in cnfg:
            varplot.plot_xy(unmask(vdata), cnfg['xyplot'], plcf, outdir)

        else :
            print("\n  ++ ERROR: xyplot not defined")
            sys.exit(1)

            
    elif plcf['type'] == 'vsdist' :
        if 'vsdist' in cnfg:

            plot_vsdist(vdata, cnfg['vsdist'], inclumps)

else:
    print(" ++  No plot defined")

    
print("....Done")


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
