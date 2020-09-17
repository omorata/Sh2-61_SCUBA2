#!/usr/bin/env python3

""" catgstats.py
  Script to get statistics from the catalog of clump physical parameters

  O. Morata 2020

  Things to do:
   - read command line arguments

   - process config file (if any)

   - use physpar and/or findclumps catalogs
   - select by clumps

   - other selection

   - output results:
     - x,y plots
     - histograms
     - KDEs

"""

import numpy as np
import numpy.ma as ma

import sys
import yaml

import ClumpCatalog as cl
import stats_utils as util
import varplot


##-- Functions ---------------------------------------------------------


def read_configfile(fname):
    """ Read the YAML configuration file."""
    
    print("  >> reading configuration file:", fname, "...")

    with open(fname, 'r') as ymlfile:
        try:
            cnfg = yaml.safe_load(ymlfile)
        except yaml.YAMLError as exc:
            print(exc)
        
    return cnfg



def read_dataset(cfg, wkdir):
    
    if not 'infile' in cfg:
        print("  ++ ERROR: No input file. Nothing to do\n     Bye!")
        sys.exit(1)

    varmap, ivar = read_catgvar(cfg, wkdir)

    auxfiles, auxpar, auxcond = read_aux(cfg, wkdir)

    varmap = filter_auxfiles(varmap, ivar, auxfiles, auxpar, auxcond)

    if hasattr(varmap, 'mask') :
        vdata = unmask(varmap)
        
        return vdata
    else :
        return varmap

    
    
def unmask(var) :

    var_unmask = []
                          
    for v in range(len(var)):
        v_tmp = var[v]
        #print(v_tmp, v_tmp.mask)
        v_tmp = v_tmp[v_tmp.mask == False]
        var_unmask.append(v_tmp)

    return var_unmask



def read_aux(cfg, wkdir):

    if 'auxfile' in cfg:
        auxf = cfg['auxfile']
        numaux = len(auxf)

        for i in range(numaux):
            auxf[i] = wkdir+auxf[i]
            #print("auxf", numaux)

        if 'auxpar' in cfg :
            auxpar = cfg['auxpar']
            if np.shape(auxpar)[0] != numaux:
                print("  ++ ERROR: the number of aux files is different",
                      "than the number of auxpar sets")
                sys.exit(1)
                
            #print("paraux", np.shape(auxpar))

        if 'conds' in cfg :
            auxconds = cfg['conds']
            #print("conds",  np.shape(auxconds))

            if np.shape(auxconds)[0] != numaux:
                print("  ++ ERROR: the number of aux files is different",
                      "from the number of condition sets")
                sys.exit(1)

            if np.shape(auxconds)[1] != 2 * np.shape(auxpar)[1]:
                print("  ++ ERROR: wrong number of conditions for filters")
                sys.exit(1)
    else:
        auxf = None
        auxpar = None
        auxconds = None

    return auxf, auxpar, auxconds


        
def read_catgvar(cfg, wkdir) :
    """Read the variables in par from the catalogue file infile"""

    ivar = len(cfg['params'])

    datparm = []

    for i in range(ivar):
        infile = wkdir+cfg['infile']
        
        catlg = cl.ClumpCatalog.from_file([infile])

        size = np.shape(catlg.clumps)[0]

        var_array = np.empty([size])

        for j in range(size):
            var_array[j] = catlg.clumps[j].record[cfg['params'][i]]

        datparm.append( var_array)

    return datparm, ivar



def filter_auxfiles(data, isets, auxf, auxpar, auxconds):

    
    if auxf:
        data = ma.asarray(data)

        for f in range(len(auxf)):
            catlg = cl.ClumpCatalog.from_file([auxf[f]])
            size = np.shape(catlg.clumps)[0]


            npars = len(auxpar[f])
            tmparr = np.zeros((npars, size))
            for i in range(size):
                for j in range(npars) :
                    tmparr[j,i] = catlg.clumps[i].record[auxpar[f][j]]


            tmparr = ma.asarray(tmparr)
            fcond = auxconds[f]
            col = 0
            while fcond:
                cond = fcond.pop(0)
                cut = fcond.pop(0)

                if cond == ">" :
                    mask_aux = ma.masked_greater(tmparr[col,:],cut )
                elif cond == "=":
                    mask_aux = ma.masked_equal(tmparr[col,:], cut)
                elif cond == "<" :
                    mask_aux = ma.masked_less(tmparr[col,:], cut)
                elif cond == "!=":
                    mask_aux = ma.masked_not_equal(tmparr[col,:], cut)
                col += 1
                
                tmparr.mask = mask_aux.mask

        if np.shape(tmparr.mask):
            data.mask = tmparr.mask
        else :
            data = data.data

    return data


##-- End of funtioncs --------------------------------------------------

args = util.read_command_line()


wdir = util.check_dir(args.wkdir)
outdir = util.check_dir(args.odir)


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
        #print("---", np.shape(vdata))
                         
else :
    print("  ++ ERROR: no dataset definition found")
    sys.exit(1)




if 'plot' in cnfg:
    plcf = varplot.read_plot_configuration(cnfg['plot'])

    
    if plcf['type'] == 'histo':
        if 'histo' in cnfg:
    
            outb = varplot.plot_histo(vdata, cnfg['histo'], plcf, outdir)

            #if outbins :
            #    save_bins(outdir+outbins, outb)
 
        else:
            print("\n  ++ ERROR: histogram not defined")
            sys.exit(1)
        
    elif plcf['type'] == 'xyplot':
        if 'xyplot' in cnfg:
            varplot.plot_xy(vdata, cnfg['xyplot'], plcf, outdir)

        else :
            print("\n  ++ ERROR: xyplot not defined")
            sys.exit(1)

else:
    print(" ++  No plot defined")

    
print("....Done")

#plot_histo(vars[:,show], 'catg.pdf')

#print("\n  ...done\n")

#a = 4
#b = 3
#x = vars[:,a]
#y = vars[:,b]

#c = slice(0,10)
#d = slice(30,100)

#fig = plt.figure(figsize=(4,4))
#plt.xscale('log')
#plt.yscale('log')
#plt.plot(x, y, '.', color='black')
#plt.plot(vars[c,a], vars[c,b], '.', color='red')
#plt.plot(vars[d,a], vars[d,b], '.', color='blue')
#fig.savefig("t.pdf")
