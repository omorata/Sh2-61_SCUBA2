#!/usr/bin/env python3
#

import numpy as np
import numpy.ma as ma

from scipy import optimize

# definition of the custom classes
#
import Param as par
import MapClass as maps


def get_maptemperature(mapratio, ini, pre, par):
    """Calculates the temperature of the map

    It calculates the temperature from the ratio of the two fluxes and
    also the variance in the temperature.
    """

    new = maps.Map.empty()
    
    new.data[0] = get_temperature(mapratio.data[0]/pre, ini,
                                  (par.hk850, par.hk450))


    new.data[1] = get_temp_variance(new.data[0], mapratio.data[1],
                                    par.hk850, par.hk450, pre)


    new.data[0] = ma.masked_where(ma.getmask(new.data[1]), new.data[0])
    
    return new
    


def get_temperature(ratio, ini_value, vargs):
    """ calculate the temperature from the flux ratio
    """

    g = lambda x, a, b, c: c * x**b - x**a + 1. - c
    
    root = optimize.newton(
        g, ini_value, args=(vargs[0], vargs[1], ratio), tol=1e-8, maxiter=250)

    mask_root = np.ma.masked_where(np.ma.getmask(ratio), root)

    t = np.ma.divide(1., np.ma.log(mask_root))
    
    return t



def get_temp_variance(temp, ratio_var, trmA, trmB, pre_fct):
    """ calculate the variance of the temperature from the variance in
        the flux ratio
    """

    r_var = ratio_var / pre_fct / pre_fct

    expA = np.ma.exp(np.ma.divide(trmA, temp))
    expB = np.ma.exp(np.ma.divide(trmB, temp))

    
    num = temp * temp * (expB - 1) * (expB - 1)
    den = expA * expB * (trmB-trmA) + trmA * expA - trmB * expB
    fct = (np.ma.divide(num,den))**2
    
    t_var = fct * r_var

    return t_var
