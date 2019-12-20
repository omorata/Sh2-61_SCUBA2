#!/usr/bin/env python3
#

import sys

import numpy as np
import numpy.ma as ma

import astropy.wcs as wcs
from astropy.io import fits
from astropy import units as u
from astropy import constants as const

from scipy import optimize



##-- Functions ---------------------------------------------------------


def get_temperature(ratio, ini_value, vargs):
    """ calculate the temperature from the flux ratio
    """
    g = lambda x, a, b, c: c * x**b - x**a + 1. - c
    
    root = optimize.newton(g, ini_value,
                           args=(vargs[0], vargs[1], ratio), 
                           tol=1e-8, maxiter=250)

    mask_root = np.ma.masked_where(np.ma.getmask(ratio), root)

    t = 1. / np.log(mask_root)
    
    return t



def get_temp_variance_K(temp, ratio_var, trmA, trmB, pre_fct):
    """ calculate the variance of the temperature from the variance in
        the flux ratio
    """
    r_var = ratio_var / pre_fct / pre_fct
    expA = np.exp(trmA / temp)
    expB = np.exp(trmB / temp)
    
    num = temp * temp * (expB - 1) * (expB - 1)
    den = expA * expB * (trmB-trmA) + trmA * expA - trmB * expB
    fct = (num / den)**2
    
    t_var = fct * r_var
    return t_var



def calc_mass(flux, var_flux, temp_arr, var_temp, d, dtog, mH, mu, solangle,
              nu, beta, hk) :

    flux_mask = np.ma.masked_where(np.ma.getmask(var_temp), flux)

    print("    >> calculating Bnu...")
    bnu = planck_u(nu, temp_arr)
    print("       ... done")
    
    knu = absorption_coefficient("freq", nu, beta) / dtog
    
    print("    >> calculating dust opacities...")
    dust_op = dust_opacity(flux_mask, solangle, bnu)
    print("       ... done")


    print("    >> calculating column density...")
    dcol_h2 = col_h2(dust_op, (knu.si).value, mu, mH.value)
    print("       ... done")
    
    print("    >> calculating masses...")
    mass = mass_h2(dcol_h2, d.to(u.m), solangle, mu, mH.value)
    mass_array = np.ma.masked_where(np.ma.getmask(temp_arr), mass)
    print("       ... done")
    
    mass_th = mass_h2_thin(flux_mask, d.to(u.m), knu, bnu)
    mass_thin = np.ma.masked_where(np.ma.getmask(temp_arr), mass_th)


    var_mass_thin = get_mass_variance_thin(mass_thin, flux_mask,
                                           var_flux, temp_arr,
                                           var_temp, hk)
    
    return mass_array, mass_thin, var_mass_thin



def get_mass_variance_thin(mass, flux, var_flux, temp, var_temp, hk):

    hkt = hk / temp
    fct = mass / flux

    varT_term = flux * hkt / temp / (1. - np.exp(-hkt)) 
    varterm = var_flux + varT_term * varT_term * var_temp
    varM = fct * fct * varterm
    return varM
    


def dust_opacity(flux, solangle, Bnu):
    """ calculate dust opacity
    """
    trm1 = (flux / solangle.value / Bnu)
    tau = - np.log(1. - trm1)

    return tau



def planck_u(nu, temp) :
    """Planckian function
    """
    f1 = 2. * const.h * nu * nu * nu / const.c / const.c

    ct = (const.h * nu / const.k_B).value
    vexp = np.divide(ct, temp)
    #vexp = (const.h * nu / const.k_B / t )

    f2 = np.exp(vexp) - 1.

    bnu = np.ma.masked_where(np.ma.getmask(temp), f1.value / f2)

    return bnu



def col_h2(tau, k, mu, mH):
    """ column density of H2 N(H2)
    """
    colh2 = tau / mu / mH / k

    return colh2



def absorption_coefficient(valtype, val, beta) :
    """calculate k_lambda using Clarke et al. 2016 values
       frequencies or wavelenghts should have units
    """  
    if valtype == "freq" :
        lamb = const.c / val
        
    else :
        lamb = val

    lamb = lamb.to(u.micron)

    k_d = 0.051 * u.m * u.m / u.kg
    #k_d = 0.51 * u.cm * u.cm / u.g
    l_0 = 500 * u.micron
    k =  k_d * ( l_0 / lamb) ** (beta)

    return k



def mass_h2(N_h2, d, solangle, mu, mH):
    """calculate H2 mass
    """
    M_h2 = mu * mH * d * d * N_h2 *  solangle / const.M_sun

    return M_h2.value



def mass_h2_thin(flux, d, k, Bnu) :
    M_h2 = flux * d * d / k / Bnu / const.M_sun

    return M_h2



def filter_temperature(temp, vartemp, type_filter, cut):
    """ filter pixels using the temperature variance
    """
    if type_filter == "variance" :
        varfilter = np.ma.masked_where(vartemp > cut, vartemp)

    elif type_filter == "snr" :
        snrT = temp / np.sqrt(vartemp)

        varTsnr = np.ma.masked_where(snrT < cut, snrT)
        varfilter = np.ma.masked_where(np.ma.getmask(varTsnr), vartemp)

    elif type_filter == "None":
        varfilter = np.ma.asarray(vartemp)
        
    else :
        print(" ++ ERROR: unknown filter temperature option")
        sys.exit(1)
        
    return varfilter



def filter_mass(mass, var_mass, cut):
    """filter pixel dependeing on mass SNR
    """
    snrM = mass / np.sqrt(var_mass)

    filtermass = np.ma.masked_where(snrM < cut, snrM)

    return filtermass


##-- End of functions --------------------------------------------------


distance = 2500. * u.pc
dtogas = 160.

mH = 1.6733e-27 * u.kg
mu = 2.8

beta = 1.8

beam = 17.5461 * u.arcsec
pixsize = 3. * u.arcsec

pixelsbeam = np.pi / 4. / np.log(2.) * beam * beam / pixsize / pixsize

pixsize_rad = pixsize.to(u.radian)
pixsolangle = pixsize_rad * pixsize_rad

hk = const.h * const.c / const.k_B

l450 = 450e-6 * u.m
l850 = 850e-6 * u.m

f850 = const.c / l850

hk850 = hk / l850 / u.K
hk450 = hk / l450 / u.K


pre = (850. / 450.)**(3.+ beta)


print(" >> reading 850 micron data...")

hdu_data850 = fits.open('analysis_maps/Sh2_61-j850_r0_contamination_mf.fits')
data850 = hdu_data850[0].data
var850 = hdu_data850[1].data
#w1 = wcs.WCS(hdu_data850[0].header)

hdu_snr850 = fits.open('analysis_maps/Sh2_61-j850_r0_contamination_mf-snr.fits') 
snr850 = hdu_snr850[0].data



print(" >> reading clump mask...")

hdumask = fits.open('findclumps/Sh2-61-j850_r0_contamination_mf-fw-snr6-extr.FITS')
clump_mask = hdumask[0].data
n_clumps = np.int(np.nanmax(clump_mask))

clump_ma = clump_mask.view(ma.MaskedArray)



print(" >> reading 450 micron data...")
hdu_data450 = fits.open('analysis_maps/Sh2_61-j450_r0mf.fits')
data450 = hdu_data450[0].data
var450 = hdu_data450[1].data

hdu_snr450 = fits.open('analysis_maps/Sh2_61-j450_r0mf-snr.fits')
snr450 = hdu_snr450[0].data



# to avoid complains about NaNs
#
with np.errstate(invalid='ignore'):
    low450 = np.ma.masked_where(snr450 < 4., snr450)


tt = ma.masked_invalid(clump_ma)
inclumps = ma.masked_less(tt, 1)


low450_clumps = np.ma.masked_where(np.ma.getmask(low450), inclumps)

good450 = np.ma.masked_where(np.ma.getmask(low450_clumps), data450)

clumps850 = np.ma.masked_where(np.ma.getmask(inclumps), data850)

good850 = np.ma.masked_where(np.ma.getmask(low450_clumps), data850)


print(" >> calculating flux ratio...")
ratio = good450 / good850

var_ratio = (ratio * ratio) * (var850 / good850 / good850 +
                               var450 / good450 / good450)

print("  ...done")


print("good850", good850.count())

bad850 = np.ma.masked_where(~np.ma.getmask(ratio), inclumps)

#dbad = np.ma.masked_where(np.ma.getmask(bad850), data850)
#print("bad850", bad850.count())

#tot = ~np.ma.mask_or(dbad,good850)
#print(tot[325,300:325])

#t1 = np.ma.masked_where(np.ma.getmask(tot), data850)
#t = np.ma.masked_where(np.ma.getmask(inclumps),t1)
#print(t[325,300:325])
#print(t.count())

pre_ratio = ratio / pre
#print("pre", pre_ratio.count())
#print(np.shape(pre_ratio))
#print(pre_ratio[325,300:325])

t1 = np.ma.asarray(ratio)
ini_array = t1 / t1 * 3.

with np.errstate(invalid='ignore'):
    t_array = get_temperature(pre_ratio, ini_array, (hk850, hk450))

var_temp = get_temp_variance_K(t_array, var_ratio, hk850, hk450, pre)

varT_filter = filter_temperature(t_array, var_temp, "snr", 5.)


print("temp", t_array.count())
print("var_ratio", var_ratio.count())
print("varT", var_temp.count())
print("varTfilter", varT_filter.count())


#high_var = np.ma.masked_where(~np.ma.getmask(varTsnr), var_temp)
#print("varhigh", high_var.count())


#
# flux was in mJy, all still in SI
#
print(" convert flux to SI...")
flux_factor = 1e-26 / pixelsbeam / 1000.
S_850 = good850 *  flux_factor
varS_850 = var850 * flux_factor * flux_factor
print("... Done!")

print(" >> calculating masses...")
mass, thin_mass, var_thin_mass = calc_mass(S_850, varS_850,
                                           t_array, varT_filter,
                                           distance, dtogas, mH, mu,
                                           pixsolangle, f850, beta, hk850)
print(" ...Done!")

print("thin_mass", thin_mass.count())
print("var_thin_mass", var_thin_mass.count())

varMsnr = filter_mass(thin_mass, var_thin_mass, 2.)

#print("varMsnr", varMsnr.count())



print(" >> clump masses")
totpix = 0
for clump in range(1, n_clumps+1):
    mskcl = ma.masked_not_equal(clump_ma, clump)
    npix = mskcl.count()

    masscl = np.ma.masked_where(np.ma.getmask(mass), mskcl)

    masscl_thin = np.ma.masked_where(np.ma.getmask(varMsnr), mskcl)

    xxx = np.ma.masked_where(np.ma.getmask(masscl), mass)
    xxx_thin = np.ma.masked_where(np.ma.getmask(masscl_thin), thin_mass) 
    xxx_var_thin_mass = np.ma.masked_where(np.ma.getmask(masscl_thin), varMsnr)
    valid_pix = xxx_thin.count()
    totpix += valid_pix

    #print(clump, " pixels", valid_pix, " sum", np.nansum(xxx), " sumthin",
    #      np.nansum(xxx_thin), "std", np.sqrt(np.nansum(xxx_var_thin_mass)))

      
print("total pixels:", totpix)

