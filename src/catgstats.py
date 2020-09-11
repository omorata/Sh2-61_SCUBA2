#!/usr/bin/env python3

""" catgstats.py
  Script to get statistics from the catalog of clump physical parameters

  O. Morata 2020

  Things to do:
   - read command line arguments

   - process config file (if any)

   - select by clumps

   - other selection

   - output results:
     - x,y plots
     - histograms
     - KDEs

"""

import numpy as np

import matplotlib.pyplot as plt
from astropy.visualization import hist

import ClumpCatalog as cl


##-- Functions ---------------------------------------------------------

def read_catgvariables(fname, par) :
    """Read the variables in par from the catalogue file infile"""
    
    catlg = cl.ClumpCatalog.from_file([fname])

    size = np.shape(catlg.clumps)[0]
    parms = np.shape(par)[0]
    var_array = np.empty([size,parms])

    for i in range(size):
        for j in range(parms):
            var_array[i,j] = catlg.clumps[i].record[par[j]]

    return var_array



def plot_histo(show_array, outplot):
    
    max = np.max(show_array)
    min = np.min(show_array)
    fig, ax = plt.subplots(1,2,figsize=(10,4))

    for i, bins in enumerate(['scott', 'freedman']):
    #for i, bins in enumerate(['knuth', 'blocks']):
        a = hist(show_array, bins=bins, ax=ax[i], histtype='stepfilled',
                 alpha=0.9, density=False, range=(min,max), stacked=True,
                 cumulative=False, log=False)
        #print(a)

        fig.savefig(outplot)

##-- End of funtioncs --------------------------------------------------

infile = "results/analysis_maps/Sh2_61-j850r0_co_mb__j450r0_mb-fw_01-clump_table.fits"

params = ['weightT', 'npix', 'N', 'mass', 'deconv_radius', 'flux850', 'flux450',
          'flux450hi', 'var_N', 'var_mass']

show = 2

print("\n >> Reading catalog")
vars = read_catgvariables(infile, params)


plot_histo(vars[:,show], 'catg.pdf')

print("\n  ...done\n")

a = 4
b = 3
x = vars[:,a]
y = vars[:,b]

c = slice(0,10)
d = slice(30,100)

fig = plt.figure(figsize=(4,4))
plt.xscale('log')
plt.yscale('log')
plt.plot(x, y, '.', color='black')
plt.plot(vars[c,a], vars[c,b], '.', color='red')
plt.plot(vars[d,a], vars[d,b], '.', color='blue')
fig.savefig("t.pdf")
