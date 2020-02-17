#!/usr/bin/env python3

import numpy as np
from astropy import units as u



class Param(object):

    def __init__(self, mu='2.8', d=100*u.pc, dtogas=160., beta=1.8,
                 beam=14.2*u.arcsec, pixsize=3.*u.arcsec, l450=450e-6*u.m,
                 l850=850e-6*u.m) :

        self.mu = mu
        self.mH = 1.6733e-27 * u.kg
        self.d = d
        self.pixsize = pixsize
        self.beam = beam
        
        
        self.pixarea = pixsize * pixsize

        pixelsbeam = np.pi / 4. / np.log(2.) * beam * beam / pixsize / pixsize

        self.flux_factor = 1e-26 / pixelsbeam / 1000.
