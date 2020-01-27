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

    t = np.ma.divide(1., np.ma.log(mask_root))
    
    return t



def get_temp_variance_K(temp, ratio_var, trmA, trmB, pre_fct):
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



def show_values(msk_arr, txt):
    """ auxiliar function to show some values of the msk_arr array
    """
    print("values for" , txt)
    print(msk_arr[324:326,300:325])
    print(txt, msk_arr.count())



def read_inputfiles(fn_data, fn_snr, txt):
    """ read input data and snr files
    """
    print(" >> reading",txt,"data...")
    
    hdu_data = fits.open(fn_data)
    data = hdu_data[0].data
    var = hdu_data[1].data
    #w1 = wcs.WCS(hdu_data850[0].header)

    hdu_snr = fits.open(fn_snr) 
    snr = hdu_snr[0].data

    return data, var, snr



def merge_masked_arrays(a, b) :
    """ merges two masked arrays a and b
    """
    merged = a.copy()

    merged[merged.mask] = b[merged.mask]
    return merged


##-- End of functions --------------------------------------------------


# parameters
#
mH = 1.6733e-27 * u.kg
mu = 2.8

l450 = 450e-6 * u.m
l850 = 850e-6 * u.m

distance = 2500. * u.pc
dtogas = 160.

beta = 1.8

beam = 17.5461 * u.arcsec
pixsize = 3. * u.arcsec

# default initial temperature estimate
#
ini_Testimate = 3.

sigma_cut450 = 3.
sigma_cutT = 3.
sigma_cutM = 2.

# calculate some constants
#
pixelsbeam = np.pi / 4. / np.log(2.) * beam * beam / pixsize / pixsize

pixsize_rad = pixsize.to(u.radian)
pixsolangle = pixsize_rad * pixsize_rad

hk = const.h * const.c / const.k_B

f850 = const.c / l850

hk850 = hk / l850 / u.K
hk450 = hk / l450 / u.K


pre = (850. / 450.)**(3.+ beta)


# file names
#
fname_850 = 'analysis_maps/Sh2_61-j850_r0_contamination_mf.fits'
fname_450 = 'analysis_maps/Sh2_61-j450_r0mf.fits'
fname_snr850 = 'analysis_maps/Sh2_61-j850_r0_contamination_mf-snr.fits'
fname_snr450 = 'analysis_maps/Sh2_61-j450_r0mf-snr.fits'

fname_clumps = 'findclumps/Sh2-61-j850_r0_contamination_mf-fw-snr6-extr.FITS'


data850, var850, snr850 = read_inputfiles(fname_850, fname_snr850, "850micron")
data450, var450, snr450 = read_inputfiles(fname_450, fname_snr450, "450micron")


print(" >> reading clump mask...")

hdumask = fits.open(fname_clumps)
clump_def = hdumask[0].data
n_clumps = np.int(np.nanmax(clump_def))

clump_idxs = clump_def.view(ma.MaskedArray)

# mask out invalid values in clump mask and pixels not in any clump
#
clump_idxs_invalid = ma.masked_invalid(clump_idxs)
inclumps = ma.masked_less(clump_idxs_invalid, 1)

#show_values(clump_idxs_invalid, "clump_idxs_invalid")
#show_values(inclumps,"inclumps")


# to avoid complains about NaNs
#
with np.errstate(invalid='ignore'):
    high450 = np.ma.masked_where(snr450 < sigma_cut450, snr450)



high450_idx = np.ma.masked_where(np.ma.getmask(high450), inclumps)

clumps_hi450 = np.ma.masked_where(np.ma.getmask(high450_idx), data450)

clumps_hi850 = np.ma.masked_where(np.ma.getmask(inclumps), data850)

doublef_cl850 = np.ma.masked_where(np.ma.getmask(high450_idx), data850)

singlef_cl850_idx = np.ma.masked_where(~np.ma.getmask(doublef_cl850), inclumps)

singlef_cl850 = np.ma.masked_where(np.ma.getmask(singlef_cl850_idx), data850) 



joint_array = merge_masked_arrays(doublef_cl850, singlef_cl850)



#show_values(high450,"low_450")
#show_values(high450_idx, "high450_idx")
#show_values(clumps_hi450, " clumps_hi450")
#show_values(clumps_hi850, " clumps_hi850")
#show_values(doublef_cl850, " doublef_cl850")
#show_values(singlef_cl850_idx, "singlef_cl850_idx")
#show_values(singlef_cl850, "singlef_cl850")
#show_values(joint_array, "joint")


print(" >> calculating flux ratio...")
ratio = np.ma.divide(clumps_hi450, doublef_cl850)

var_ratio = (ratio * ratio) * (var850 / doublef_cl850 / doublef_cl850 +
                               var450 / clumps_hi450 / clumps_hi450)

#show_values(ratio, "ratio")
#show_values(var_ratio, "var_ratio")
print("  ...done")





pre_ratio = ratio / pre

t1 = np.ma.asarray(ratio)
ini_array = t1 / t1 * ini_Testimate

with np.errstate(invalid='ignore'):
    t_array = get_temperature(pre_ratio, ini_array, (hk850, hk450))

var_temp = get_temp_variance_K(t_array, var_ratio, hk850, hk450, pre)


varT_filter = filter_temperature(t_array, var_temp, "snr", sigma_cutT)

show_values(t_array, "temp")
show_values(var_temp, "temp")
show_values(varT_filter, "temp")
#print("temp", t_array.count())
#print("var_ratio", var_ratio.count())
#print("varT", var_temp.count())
#print("varTfilter", varT_filter.count())


#high_var = np.ma.masked_where(~np.ma.getmask(varTsnr), var_temp)
#print("varhigh", high_var.count())


#
# flux was in mJy, all still in SI
#
print(" >> convert flux to SI...")
flux_factor = 1e-26 / pixelsbeam / 1000.
S_850 = doublef_cl850 *  flux_factor
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

varMsnr = filter_mass(thin_mass, var_thin_mass, sigma_cutM)

#print("varMsnr", varMsnr.count())



print(" >> clump masses")
totpix = 0
for clump in range(1, n_clumps+1):
    mskcl = ma.masked_not_equal(clump_idxs, clump)
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
