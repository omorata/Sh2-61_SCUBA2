#!/usr/bin/env python3
#

import numpy as np
import numpy.ma as ma

from scipy import optimize
import scipy.ndimage as ndimage

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



def fill_temperature(calc, manual, incl, Tdefault):
    """Fill the masked temperature pixels extrapolating the calculated
    ones
    """
    
    #Tdefault = 20.
    fp = np.array([[1,1,1], [1,0,1], [1,1,1]], np.uint8)
    dtype = [('x', int), ('y', int), ('dist', float)]

    w = np.array([0.707107,1,0.707107,1,1,0.707107,1,0.707107])
    #w = np.array([0.5,1,0.5,1,1,0.5,1,0.5])
    #w = np.array([1,1,1,1,1,1,1,1])

    print("\n  >> starting extrapolation of Tdust")

    #print(calc.data[0][323:326,300:325])
    #print(manual.data[0][323:326,300:325])
    #b = manual.getmask()
    #print(b[323:326,300:325])
    
    maxcl = np.int(np.nanmax(incl))

    for i in range(maxcl):
        cl = i + 1
        #print("Working on clump", cl)
        mskcl = ma.masked_not_equal(incl, cl)

        calccl = calc.masked_where(ma.getmask(mskcl))
        mancl = manual.masked_where(ma.getmask(mskcl))

        a = calccl.data[0]
        n = mancl.data[0]
        nm = ma.getmask(n)

        n[~n.mask] = -99
        a[~n.mask] = n[~n.mask]

        # get pixel positions of peak of the clump
        #
        maxpos = np.unravel_index(np.argmax(a, axis=None), a.shape)

        # get positions of pixels in the clump
        #
        valpos = ma.where(n)
        szarr = np.shape(valpos)[1]

        if a[maxpos[0],maxpos[1]] < 0 :
            a[maxpos[0],maxpos[1]] = Tdefault

        idx = 0
        sequ = []
        for i in range(szarr):
            x,y = (valpos[0][i],valpos[1][i])
            dx = maxpos[0] - x
            dy = maxpos[1] - y
            
            d = np.sqrt(dx*dx + dy*dy)
            if d < 1e-5 :
                continue
            
            sequ.append((x, y, d))
            idx += 1
            
        seqarr = np.array(sequ, dtype)
        sortedarr = np.sort(seqarr, order='dist')

        interpolated_values = []
        sizesorted = np.shape(sortedarr)[0]
        for j in range(sizesorted):
            x = sortedarr['x'][j]
            y = sortedarr['y'][j]

            nng = a[x-1:x+2,y-1:y+2]
            nng[nng.mask] = -100
        
            results = ndimage.generic_filter(nng, nn_func, footprint=fp,
                                             mode='constant',
                                             extra_arguments=(w,))
            
            if (not np.isfinite(results[1,1])):
                # if we get a nan, because of whatever, we use the
                # average of all the values so far (just because)
                #
                a[x,y] = np.average(interpolated_values)

            else :
                a[x,y] = results[1,1]
                interpolated_values.append(a[x,y])


        # update the manual array with the new values (important)
        #
        manual.data[0][~nm] = a[~nm]

    #print("+++++++++++++++++++++++++++")
    #print(calc.data[0][323:326,300:325])
    #print("++++++++++++++++++++++++++----+")
    print("MANUAL DATA")
    print(manual.data[0][323:326,300:325])
    #print("+++++++++++++++++++++++++++")
    #b = manual.getmask()
    #print(b[323:326,300:325])
    #print("+++++++++++++++++++++++++++")

    
    print("  .... extrapolation finished\n\n")
    return manual



def nn_func(a, weights):
    """Calculate weighted array"""

    a[a < 0] = 0

    b = a.copy()
    b[b > 0] = 1


    a = np.multiply(a,weights)
    c = np.multiply(b,weights)
    a = np.sum(a) / np.sum(c)

    return a
