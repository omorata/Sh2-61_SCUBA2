#!/usr/bin/env python3
#
##  Param.py
##  O. Morata 2020
##
""" Definition of the class Param that contains several constants of the
    environment of the calculations.
"""


import numpy as np
from astropy import units as u
from astropy import constants as const


class Param(object):

    def __init__(self, mu='2.8', d=100*u.pc, dtogas=160., beta=1.8,
                 beam=14.2*u.arcsec, pixsize=3.*u.arcsec, l450=450e-6*u.m,
                 l850=850e-6*u.m) :

        self.mu = mu
        self.mH = 1.6733e-27 * u.kg
        self.dtogas = dtogas
        self.d = d.to(u.m)
        self.beta = beta
        self.l450 = l450
        self.l850 = l850
        self.beam = beam
        self.pixsize = pixsize
        
        self.pixarea = pixsize * pixsize

        self.pixelsbeam = np.pi / 4. / np.log(2.) * beam * beam / pixsize / pixsize

        hk = const.h * const.c / const.k_B
        self.hk850 = hk / l850 / u.K
        self.hk450 = hk / l450 / u.K

        self.nu = const.c / l850

        

    @classmethod
    def from_cfgfile(cls, cfg, section):
        cf = cfg[section]
        
        mu = float(cf['mu'])
        beta = float(cf['beta'])
        distance = float(cf['distance']) * u.pc
        dtogas = float(cf['dtogas'])
        beam = float(cf['beam']) * u.arcsec
        pixsize = float(cf['pixsize']) * u.arcsec

        new = cls(mu=mu, d=distance, dtogas=dtogas, beta=beta, beam=beam,
               pixsize=pixsize)

        return new


    
    def get_flux_factor(self, map) :
        """Calculate flux factor from the header of the main map."""
        
        map_units = map.header[0]['BUNIT']
        
        if map_units == 'mJy/beam' :
            self.flxunits = 1e-29
            self.flux_factor = self.flxunits / self.pixelsbeam 

        elif map_units == 'Jy/beam' :
            self.flxunits = 1e-26
            self.flux_factor = self.flxunits / self.pixelsbeam

        else :
            print("ERROR: undefined flux units in input map")
            sys.exit(1)
