#!/usr/bin/env python3
#

import sys

import numpy as np
#import numpy.ma as ma

#from astropy.io import fits
from astropy import units as u
from astropy import constants as const

#import math

#import argparse as ap
#import yaml

#import logging

# definition of the custom classes
#
import Param as par
#import ClumpCatalog as cl
import MapClass as maps


def calc_mapmass(flux, temp, par, calc_type='thin', solangle=0.):
    """Calculate the masses in the map.

       If calc_type = 'thin', follow the optically thin approximation. 
       If calc_type = 'tau', calculate opacity of the emission and then
                             calculate column densities and masses.
    """
    
    flux_mask = flux.masked_where(temp.getmask())

    bnu = planck_u(par.nu, temp.data[0])
    
    knu = absorption_coefficient("freq", par.nu, par.beta) / par.dtogas

    
    if calc_type == 'thin' :
        #mpmass = mapmass_h2_thin(flux_mask, temp, (par.d).to(u.m), knu, bnu,
        mpmass = mapmass_h2_thin(flux_mask, temp, par.d, knu, bnu,
                                 par.hk850)
        return mpmass
    
    elif calc_type == 'tau' :
        mass_tau, tau, cd_tau = mapmass_h2_tau(flux_mask, solangle, bnu, knu,
                                               par)
        
        return mass_tau, tau, cd_tau
        
    else :
        print(" >> ERROR: wrong mass calculation type", calc_type)
        sys.exit(1)
        


def mapmass_h2_thin(flux, temp, d, k, Bnu, hk) :
    """Calculate the pixel mass and mass variance (optically thin)."""

    mass = maps.Map.empty()

    ## calculate masses
    ##
    mass.data[0] = flux.data[0] * d * d / k / Bnu / const.M_sun

    ## calculate the variances
    ##
    hkt = hk / temp.data[0]
    fct = mass.data[0] / flux.data[0]

    term_varT= flux.data[0] * hkt / temp.data[0] / (1. - np.exp(-hkt))
    term_var = flux.data[1] + term_varT * term_varT * temp.data[1]

    mass.data[1] = fct * fct * term_var

    return mass



def mapmass_h2_tau(flux, solangle, bnu, knu, par) :
    """ Calculates the pixel opacity, column density and mass."""
    
    tau = maps.Map.empty()
    cd = maps.Map.empty()
    mass = maps.Map.empty()
        
    tau.data[0] = dust_opacity(flux.data[0], solangle, bnu)
    cd.data[0] = col_h2(tau.data[0], knu, par.mu, par.mH)
    mass.data[0] =  mass_h2(cd.data[0], par, solangle)
    
    return mass, tau, cd



def dust_opacity(flux, solangle, Bnu):
    """Calculate dust opacity."""

    trm1 = (flux / solangle.value / Bnu)
    tau = - np.log(1. - trm1)

    return tau



def planck_u(nu, temp) :
    """Planckian function."""
    
    f1 = 2. * const.h * nu * nu * nu / const.c / const.c

    ct = (const.h * nu / const.k_B).value
    vexp = np.divide(ct, temp)
    #vexp = (const.h * nu / const.k_B / t )

    f2 = np.exp(vexp) - 1.

    bnu = np.ma.masked_where(np.ma.getmask(temp), f1.value / f2)

    return bnu



def col_h2(tau, k, mu, mH):
    """Calculate the H2 column density, N(H2)."""
    
    colh2 = tau / mu / mH / k
    # convert to cm-2
    colh2 /= 10000.

    return colh2



def absorption_coefficient(valtype, val, beta) :
    """calculate k_lambda using Clarke et al. 2016 values

    IMPORTANT: Frequencies or wavelengths should have units
    """  

    if valtype == "freq" :
        lamb = const.c / val
        
    else :
        lamb = val

    lamb = lamb.to(u.micron)

    k_d = 0.051 * u.m * u.m / u.kg
    l_0 = 500 * u.micron
    k =  k_d * ( l_0 / lamb) ** (beta)

    return k



def mass_h2(N_h2, par, solangle):
    """Calculate H2 mass from the tau determined H_2 column density."""
    
    fct = get_col_factor(par, solangle)
    M_h2 = N_h2 / fct

    return M_h2



def get_col_factor(par, solangle):

    dist = (par.d).to(u.cm)
    col_factor = const.M_sun / par.mu / par.mH / dist / dist / solangle
    col_factor *= u.cm * u.cm * u.rad * u.rad
    return col_factor
