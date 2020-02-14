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

import math


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


    var_mass_thin = get_mass_variance_thin(mass_thin, flux_mask, var_flux,
                                           temp_arr, var_temp, hk)
    
    return mass_array, mass_thin, var_mass_thin



def get_mass_variance_thin(mass, flux, var_flux, temp, var_temp, hk):
    """Calculate the variance of the mass (optically thin emission).
    """
    
    hkt = hk / temp
    fct = mass / flux

    varT_term = flux * hkt / temp / (1. - np.exp(-hkt)) 
    varterm = var_flux + varT_term * varT_term * var_temp
    varM = fct * fct * varterm
    
    return varM
    


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

    return colh2



def absorption_coefficient(valtype, val, beta) :
    """calculate k_lambda using Clarke et al. 2016 values

    IMPORTANT: Frequencies or wavelenghts should have units
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



def mass_h2(N_h2, d, solangle, mu, mH):
    """Calculate H2 mass."""
    
    M_h2 = mu * mH * d * d * N_h2 *  solangle / const.M_sun

    return M_h2.value



def mass_h2_thin(flux, d, k, Bnu) :
    M_h2 = flux * d * d / k / Bnu / const.M_sun

    return M_h2



def filter_parameter(value, variance, type_filter, cut):
    """Filter pixels using the variance of the parameter

    It returns the array with the filtered variance
    """
    
    if type_filter == "variance" :
        varfilter = np.ma.masked_where(variance > cut, variance)

    elif type_filter == "snr" :
        snr = value / np.sqrt(variance)

        snr_cut = np.ma.masked_where(snr < cut, snr)
        varfilter = np.ma.masked_where(np.ma.getmask(snr_cut), variance)

    elif type_filter == "None":
        varfilter = np.ma.copy(variance)
        
    else :
        print(" ++ ERROR: unknown filter option")
        sys.exit(1)
        
    return varfilter



def show_values(msk_arr, txt):
    """Show some values of the msk_arr array."""
    
    print("values for" , txt)
    print(msk_arr[324:326,300:325])
    print(txt, msk_arr.count())
    print("  max:", np.nanmax(msk_arr))
    print("  min:", np.nanmin(msk_arr))



def read_fitsfile(fn_data, txt):
    """Read input data and header from a FITS file."""
    
    print(" >> reading",txt,"data...")
    
    with fits.open(fn_data) as hdu_data:
        #data = hdu_data[0].data
        #var = hdu_data[1].data

        data_info = [hdu_data[0].data, hdu_data[1].data]
        # reads header information
        #
        header_info = [hdu_data[0].header, hdu_data[1].header]
    
    return data_info, header_info



def merge_masked_arrays(a, b) :
    """Merge two masked arrays a and b."""
    
    merged = np.ma.copy(a)

    merged[merged.mask] = b[merged.mask]
    
    return merged



def flux450(f850, td, h850, h450, pre) :
    """Calculate the flux at 450micron from the flux at 850 micron and a
        given dust temperature
    """
    
    x = np.exp( 1. / td)
    f = pre * f850 * (x**h850 - 1) / (x**h450 - 1)
    return f



def get_variance_ratio(r, f850, f450, v850, v450):
    """Calculate the variance of the flux ratio."""
    
    var_r = (r * r) * (v850 / f850 / f850 + v450 / f450 / f450)
    return var_r



def calc_perc(array, p):
    """Calculate percentile for a masked array."""
    
    n_arr = np.ma.copy(array)
    size = n_arr.count()
    x = n_arr.flatten()
    x.sort()
    
    if size > 0 :
        return x[math.ceil((size * p / 100) - 1)]
    else:
        return None
    


def clump_weighted_avg(array, weights, clump_mask) :

    cl_a = np.ma.masked_where(np.ma.getmask(clump_mask), array)
    cl_w = np.ma.masked_where(np.ma.getmask(clump_mask), weights)
        
    return np.nansum(cl_a * cl_w) / np.nansum(cl_w)



def get_size(n, pixsize, beamsize) :
    """Calculate some geometrical parameters of the clump."""

    area = n* pixsize * pixsize
    ef_rad = np.sqrt(area / np.pi)
    dec_r = np.sqrt(4 * ef_rad * ef_rad -
                    (np.pi / 4. / np.log(2.) )* beamsize * beamsize) * 0.5
    
    return area, ef_rad, dec_r



def columndensity(mass, var_mass, area, d, mu, mH ) :
    """Calculate column density from a mass and an area

    mass in solar masses, area in arcsec^2, d in pc.
    Returns the column density in cm^-2 and the variance of the column
    density in cm^-4
    """

    arad = area.to(u.radian * u.radian)
    dcm = d.to(u.cm)

    solangle = 2. * np.pi * (1. - np.cos(2. * np.sqrt(arad/np.pi)))

    fact = const.M_sun / mu / mH / dcm / dcm / solangle
    
    col = mass * fact
    var_col = var_mass * fact * fact

    return col, var_col


##-- End of functions --------------------------------------------------


# constants
#
mH = 1.6733e-27 * u.kg
mu = 2.8

l450 = 450e-6 * u.m
l850 = 850e-6 * u.m

distance = 2500. * u.pc
dtogas = 160.


## Parameters
##

# file names
#
fname_850 = 'analysis_maps/Sh2_61-j850_r0_contamination_mf.fits'
fname_450 = 'analysis_maps/Sh2_61-j450_r0mf.fits'
fname_snr850 = 'analysis_maps/Sh2_61-j850_r0_contamination_mf-snr.fits'
fname_snr450 = 'analysis_maps/Sh2_61-j450_r0mf-snr.fits'

fname_clumps = 'findclumps/Sh2-61-j850_r0_contamination_mf-fw-snr6-extr.FITS'


beta = 1.8

beam = 17.5461 * u.arcsec
pixsize = 3. * u.arcsec

# default initial temperature estimate
#
ini_Testimate = 3.

sigma_cut450 = 4.

type_cutTd = "snr"
cut_Td = 3.

type_cutM = "snr"
sigma_cutM = 2.

manual_Tdust = 15.
##
## End of parameters


# calculate some variables
#
pixelsbeam = np.pi / 4. / np.log(2.) * beam * beam / pixsize / pixsize

pixarea = pixsize * pixsize

# input fluxes in mJy/beam
#
flux_factor = 1e-26 / pixelsbeam / 1000.

pixsize_rad = pixsize.to(u.radian)
pixsolangle = pixsize_rad * pixsize_rad

hk = const.h * const.c / const.k_B


f850 = const.c / l850

hk850 = hk / l850 / u.K
hk450 = hk / l450 / u.K

pre = (850. / 450.)**(3.+ beta)



print(" ++ Start")
print(" >> reading input files...")

data, header850 = read_fitsfile(fname_850, "850micron")
snr, header_snr = read_fitsfile(fname_snr850, "SNR 850micron")

data850 = data[0]  
var850 = data[1]
snr850 = snr[0]

data, header450 = read_fitsfile(fname_450, "450micron")
snr, header_snr = read_fitsfile(fname_snr450, "SNR 450micron")

data450 = data[0]  
var450 = data[1]
snr450 = snr[0]


print(" >> reading clump mask...")

hdumask = fits.open(fname_clumps)
clump_def = hdumask[0].data
n_clumps = np.int(np.nanmax(clump_def))

clump_idxs = clump_def.view(ma.MaskedArray)

# mask out invalid values in clump mask and pixels not in any clump
#
clump_idxs_invalid = ma.masked_invalid(clump_idxs)
inclumps = ma.masked_less(clump_idxs_invalid, 1)


# to avoid complains about NaNs
#
with np.errstate(invalid='ignore'):
    high450 = np.ma.masked_where(snr450 < sigma_cut450, snr450)



high450_idx = np.ma.masked_where(np.ma.getmask(high450), inclumps)
clumps_hi450 = np.ma.masked_where(np.ma.getmask(high450_idx), data450)

clumps_450 = np.ma.masked_where(np.ma.getmask(inclumps), data450)

clumps_hi850 = np.ma.masked_where(np.ma.getmask(inclumps), data850)

doublef_cl850 = np.ma.masked_where(np.ma.getmask(high450_idx), clumps_hi850)
singlef_cl850_idx = np.ma.masked_where(~np.ma.getmask(doublef_cl850), inclumps)
singlef_cl850 = np.ma.masked_where(np.ma.getmask(singlef_cl850_idx), data850) 


# definition of array to hold pixels where WE fix Tdust
#
manual_temp = np.ma.copy(singlef_cl850)


print(" >> calculating flux ratios...")

ratio = np.ma.divide(clumps_hi450, doublef_cl850)

var_ratio = get_variance_ratio(ratio, doublef_cl850, clumps_hi450, var850,
                               var450)

print("  ...done")


print(" >> calculating temperatures...")

pre_ratio = ratio / pre

t1 = np.ma.copy(ratio)
ini_array = np.ma.divide(t1, t1) * ini_Testimate

with np.errstate(invalid='ignore'):
    t_array = get_temperature(pre_ratio, ini_array, (hk850, hk450))


# array with the pixels where no temperature could be calculated
#
notemp = np.ma.masked_where(~np.ma.getmask(t_array), ratio)
manual_temp = merge_masked_arrays(manual_temp, notemp)


# The temperature array must have well defined temperatures and variance
# of the temperature. The pixels that do not match that go into the
# manual_temp array
#
var_temp = get_temp_variance(t_array, var_ratio, hk850, hk450, pre)

temp = np.ma.masked_where(np.ma.getmask(var_temp), t_array)


novartemp = np.ma.masked_where(~np.ma.getmask(var_temp), t_array)
manual_temp = merge_masked_arrays(manual_temp, novartemp)


varT_filter = filter_parameter(temp, var_temp, type_cutTd, cut_Td)
temp_filter = np.ma.masked_where(np.ma.getmask(varT_filter), temp)

novarfilter = np.ma.masked_where(~np.ma.getmask(varT_filter), temp)
manual_temp = merge_masked_arrays(manual_temp, novarfilter)

print("   ...done")


#
# flux was in mJy, all still in SI
#
print(" >> convert flux to SI...")
S_850 = doublef_cl850 *  flux_factor
varS_850 = var850 * flux_factor * flux_factor
print("   ...done")


print(" >> calculating masses...")
mass, thin_mass, var_thin_mass = calc_mass(S_850, varS_850,
                                           temp_filter, varT_filter,
                                           distance, dtogas, mH, mu,
                                           pixsolangle, f850, beta, hk850)
print("   ...done")


varMsnr = filter_parameter(thin_mass, var_thin_mass, type_cutM, sigma_cutM)

mass_850 = np.ma.masked_where(np.ma.getmask(varMsnr), thin_mass)

lowM = np.ma.masked_where(~np.ma.getmask(varMsnr), thin_mass)

manual_temp = merge_masked_arrays(manual_temp, lowM)

temp_filtermass = np.ma.masked_where(np.ma.getmask(varMsnr), temp_filter)


print(" >>\n >> processing fixed dust temperature pixels...")

f850_notemp = np.ma.masked_where(np.ma.getmask(manual_temp), clumps_hi850)
var850_notemp = np.ma.masked_where(np.ma.getmask(manual_temp), var850)
var450_notemp = np.ma.masked_where(np.ma.getmask(manual_temp), var450)

new_f450 = flux450(f850_notemp, manual_Tdust, hk850, hk450, pre)

#nt_ratio = np.ma.divide(new_f450, f850_notemp)

#var_ntratio = get_variance_ratio(nt_ratio, f850_notemp, new_f450,
#                                 var850_notemp, var450_notemp)

#varT_notemp = get_temp_variance(manual_temp, var_ntratio, hk850, hk450, pre)
x = np.ma.copy(new_f450)
varT_notemp = np.ma.divide(x, x) * 30.


#snrnotemp = manual_Tdust / np.ma.sqrt(var_ntratio)


S850_notemp = f850_notemp * flux_factor
varS_850notemp = var850_notemp * flux_factor * flux_factor

mass_notemp, thin_mass_notemp, var_thin_mass_notemp = calc_mass(S850_notemp,
                                                    varS_850notemp,
                                                    manual_Tdust,
                                                    varT_notemp,
                                                    distance, dtogas, mH, mu,
                                                    pixsolangle, f850,
                                                    beta, hk850)

print(np.nansum(thin_mass_notemp), "+-",
      np.ma.sqrt(np.nansum(var_thin_mass_notemp)))
print(np.nansum(mass_850), "+-", np.ma.sqrt(np.nansum(varMsnr)))


mass_th_total = np.ma.copy(mass_850)
mass_th_total = merge_masked_arrays(mass_th_total, thin_mass_notemp)

varM_th_total = np.ma.copy(varMsnr)
varM_th_total = merge_masked_arrays(varMsnr, var_thin_mass_notemp)

print("   ...done")


print(" >> clump masses")

totpix = [0, 0, 0, 0, 0]

print("clump  #  mass_th+-errM  mass_tot+errM")

for clump in range(1, n_clumps+1):
    mskcl = ma.masked_not_equal(clump_idxs, clump)
    npix = mskcl.count()

    str_out = '{0:2d} {1:4d}'.format(clump, npix)

    area, eff_radius, deconv_r = get_size(npix, pixsize, beam)

    str_out = ' {0} {1:7.1f} {2:8.2f} {3:8.2f}'.format(str_out,
                                                    area / u.arcsec / u.arcsec,
                                                    eff_radius / u.arcsec,
                                                        deconv_r / u.arcsec)


    cl_f850 = np.ma.masked_where(np.ma.getmask(mskcl), clumps_hi850)

    cl_S850 = cl_f850 * flux_factor
    cl_flux = np.nansum(cl_S850) / 1.e-26

    str_out = ' {0} {1:7.3f}'.format(str_out, cl_flux)

    cl_f450hi = np.ma.masked_where(np.ma.getmask(mskcl), clumps_hi450)
    cl_S450hi = cl_f450hi * flux_factor
    cl_flux450hi = np.nansum(cl_S450hi) / 1.e-26
    str_out = ' {0} {1:7.3f}'.format(str_out, cl_flux450hi)

    cl_f450 = np.ma.masked_where(np.ma.getmask(mskcl), clumps_450)
    cl_S450 = cl_f450 * flux_factor
    cl_flux450 = np.nansum(cl_S450) / 1.e-26
    str_out = ' {0} {1:7.3f}'.format(str_out, cl_flux450)

     
    #masscl = np.ma.masked_where(np.ma.getmask(mass), mskcl)

    cl_td_idx = np.ma.masked_where(np.ma.getmask(temp_filtermass), mskcl)

    cl_td = np.ma.masked_where(np.ma.getmask(cl_td_idx), temp_filtermass)
        
    cl_thin_mass_idx = np.ma.masked_where(np.ma.getmask(mass_850), mskcl)
    
    cl_thin_mass = np.ma.masked_where(np.ma.getmask(cl_thin_mass_idx),
                                      mass_850) 
    cl_var_thin_mass = np.ma.masked_where(np.ma.getmask(cl_thin_mass_idx),
                                          varMsnr)

    cl_N_th, cl_var_N_th = columndensity(cl_thin_mass, cl_var_thin_mass,
                                         pixarea, distance, mu, mH)

    #show_values(cl_N_th, "coldens")
    str_out = '{0} {1:4d}'.format(str_out, cl_thin_mass.count())
    str_out = '{0} {1:8.2f} ({2:6.2f})'.format(str_out,
                                        np.nansum(cl_thin_mass),
                                        np.ma.sqrt(np.nansum(cl_var_thin_mass)))

    
    cl_tot_mass_idx =  np.ma.masked_where(np.ma.getmask(mass_th_total), mskcl)
    cl_tot_mass = np.ma.masked_where(np.ma.getmask(cl_tot_mass_idx),
                                     mass_th_total) 
    cl_var_tot_mass = np.ma.masked_where(np.ma.getmask(cl_tot_mass_idx),
                                         varM_th_total)

    str_out = '{0} {1:8.2f} ({2:6.2f})'.format(str_out, np.nansum(cl_tot_mass),
                                    np.ma.sqrt(np.nansum(cl_var_tot_mass)))

    cl_N_tot, cl_var_N_tot = columndensity(cl_tot_mass, cl_var_tot_mass,
                                         pixarea, distance, mu, mH)
    
    clump_N, clump_varN = columndensity(np.nansum(cl_tot_mass),
                                        np.nansum(cl_var_tot_mass),
                                        area, distance, mu, mH)
    
    str_out = '{0} {1:9.3e} {2:9.3e}'.format(str_out, clump_N,
                                             np.sqrt(clump_varN))
    
    if cl_td.count() > 0 :

        maxT = np.nanmax(cl_td)
        minT = np.nanmin(cl_td)
        SavgT = clump_weighted_avg(temp_filtermass, clumps_hi450, cl_td_idx)

        SavgT = calc_perc(cl_td, 25)
        str_t = '///{0:6.1f} {1:6.1f} {2:6.1f}'.format(maxT, minT, SavgT)
        
    else :
        str_t = '///{0:20s}'.format(" ")

    print(str_out, str_t)
    

    totpix[0] += npix

    totpix[1] += cl_thin_mass.count()
    
    if cl_thin_mass.count() > 0 :
        totpix[2] += np.nansum(cl_thin_mass)

    totpix[3] += cl_tot_mass.count()

    if cl_tot_mass.count() > 0 :
        totpix[4] += np.nansum(cl_tot_mass)
    
    
      
#print("total pixels:", totpix)

print(" ++ ok")
