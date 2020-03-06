#!/usr/bin/env python3

import ClumpCatalog as cl
import numpy as np
from astropy.io import fits

clfile = 'cl1.fits'
phyfile = 'cl2.fits'


catlg = cl.ClumpCatalog.fromfile([clfile,phyfile])
#catlg = cl.ClumpCatalog.fromfile([clfile])

#catlg.fillfrom_file(phyfile)

catlg.save_catalog("xred.fits", overwrite=True, ctype='phys')
catlg.save_catalog("xall.fits", overwrite=True, ctype='all')
    

catlg.print_catalog(ctype='phys')

