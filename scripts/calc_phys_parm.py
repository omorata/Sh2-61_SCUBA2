#!/usr/bin/env python3
#

import sys

import numpy as np
import numpy.ma as ma

from astropy.io import fits
from astropy import units as u
from astropy import constants as const

from scipy import optimize

import math

import argparse as ap
import yaml

# definition of the custom classes
#
import Param as par
import ClumpCatalog as cl
import MapClass as maps


##-- Functions ---------------------------------------------------------


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



def calc_mapmass(flux, temp, par):
    """Calculate the masses in the map."""
    
    flux_mask = flux.masked_where(temp.getmask())

    bnu = planck_u(par.nu, temp.data[0])
    
    knu = absorption_coefficient("freq", par.nu, par.beta) / par.dtogas
    
    mpmass = mapmass_h2_thin(flux_mask, temp, (par.d).to(u.m), knu, bnu,
                             par.hk850)

    return mpmass



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



def show_values(msk_arr, txt):
    """Show some values of the msk_arr array."""
    
    print("values for" , txt)
    print(msk_arr[324:326,300:325])
    print(txt, msk_arr.count())
    print("  max:", np.nanmax(msk_arr))
    print("  min:", np.nanmin(msk_arr))



def flux450(f850, td, h850, h450, pre) :
    """Calculate the flux at 450micron from the flux at 850 micron and a
        given dust temperature
    """
    
    x = np.exp( 1. / td)
    f = pre * f850 * (x**h850 - 1) / (x**h450 - 1)
    return f



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



def read_command_line() :
    """ Read command line arguments."""
    
    parser = ap.ArgumentParser()

    parser.add_argument(
        '-c', dest='cfgfile',  help='configuration file',
        default='config.cfg', metavar='FILE')

    parser.add_argument(
        '--clumps', action='store_true', default=False,
        help='process clumps')
    
    parser.add_argument(
        '--print', dest='print_table', action='store_true', default=False,
        help='print table')

    parser.add_argument('args', nargs=ap.REMAINDER)
    
    return parser.parse_args()



def read_configfile(fname):
    """ Read the YAML configuration file."""
    
    print("  >> reading configuration file:", fname, "...")

    with open(fname, 'r') as ymlfile:
        try:
            cnfg = yaml.safe_load(ymlfile)
        except yaml.YAMLError as exc:
            print(exc)
        
    return cnfg



def cfgval(cfg, section, key, status=''):

    if key in cfg[section] :
        return cfg[section][key]
    else :
        if status == 'required' :
            print(" >> ERROR: missing required key", key)
            sys.exit(1)
        else :
            return ''


        
def get_values(cfg, section, names=None, status='', type='str', altnames=None):
    """Read the values from the configuration file sections and keys."""

    vv = {}
    alt = False
    
    if altnames :
        alt = True
        
    for key in names :
        if alt:
            vkey = altnames.pop(0)
        else :
            vkey = key
            
        vv[vkey] = cfgval(cfg, section, key, status=status)
        if type == 'float' :
            vv[vkey] = float(vv[vkey])
            
    return vv


##-- End of functions --------------------------------------------------

l450 = 450e-6 * u.m
l850 = 850e-6 * u.m


print(" ++ Start")


args = read_command_line()



cnfg = read_configfile(args.cfgfile)

fname = get_values(cnfg, 'infiles', status='required',
                   names = ['data850', 'data450', 'snr850', 'snr450',
                            'clumpdef'])

fout = get_values(cnfg, 'outfiles', 
                   names = ['ratio', 'temperature', 'mass', 'clumptable'])


pr = par.Param.from_cfgfile(cnfg, 'phys_params')


defaults = get_values(cnfg, 'data_params', type='float',
                      names=['ini_temp','manual_Tdust','manual_varTdust'],
                      altnames=['iniTd', 'Td', 'varTd'])

typecut = get_values(cnfg, 'data_params', names=['type_cutTd', 'type_cutM'],
                     altnames=['Td', 'M'])

cuts = get_values(cnfg, 'data_params', names=['snr_450', 'cut_Td', 'cut_M'],
                  altnames=['snr450', 'Td', 'M'], type='float')



# calculate some variables
#
pixelsbeam = np.pi / 4. / np.log(2.) * pr.beam * pr.beam / pr.pixsize / pr.pixsize

pixarea = pr.pixsize * pr.pixsize

# input fluxes in mJy/beam
#
flux_factor = 1e-26 / pixelsbeam / 1000.

pixsize_rad = pr.pixsize.to(u.radian)
pixsolangle = pixsize_rad * pixsize_rad

hk = const.h * const.c / const.k_B


f850 = const.c / l850

hk850 = hk / l850 / u.K
hk450 = hk / l450 / u.K

pre = (850. / 450.)**(3.+ pr.beta)



print("  >> reading input files...")

map850 = maps.Map.from_fitsfile(fname['data850'], name="850micron")
mapsnr = maps.Map.from_fitsfile(fname['snr850'], name="SNR 850micron")
map450 = maps.Map.from_fitsfile(fname['data450'], name="450micron")
mapsnr450 = maps.Map.from_fitsfile(fname['snr450'], name="SNR 450micron")

header850 = map850.header


print("  >> reading clump mask...")

with fits.open(fname['clumpdef']) as hdumask:
    clump_def = hdumask[0].data

clump_idxs = clump_def.view(ma.MaskedArray)

# mask out invalid values in clump mask and pixels not in any clump
#
clump_idxs_invalid = ma.masked_invalid(clump_idxs)
inclumps = ma.masked_less(clump_idxs_invalid, 1)

maphi450 = mapsnr450.masked_less(cuts['snr450'])


test = ma.masked_where(maphi450.getmask(), inclumps)

mapclumpshi450 = map450.masked_where(ma.getmask(test))

mapclumps450 = map450.masked_where(ma.getmask(inclumps))


mapclumpshi850 = map850.masked_where(ma.getmask(inclumps))

mapdblf_cl850 = mapclumpshi850.masked_where(ma.getmask(test))
sf_cl850_idx = ma.masked_where(~mapdblf_cl850.getmask(), inclumps)
mapsf_cl850 = map850.masked_where(ma.getmask(sf_cl850_idx))


# definition of array to hold pixels where WE fix Tdust
#

mapmanual_temp = mapsf_cl850.copy()


print("  >> calculating flux ratios...")

map_ratio = maps.divide(mapclumpshi450, mapdblf_cl850)


if fout['ratio']:
    ok = map_ratio.save_fitsfile(
        fname=fout['ratio'], hdr_type='fluxratio', oldheader=map850.header,
        append=False, overwrite=True)

print("  ...done")


print("  >> calculating temperatures...")

ini_array = np.full_like(map_ratio.data[0], defaults['iniTd'])

with np.errstate(invalid='ignore'):
    maptemp = get_maptemperature(map_ratio, ini_array, pre, pr)

mapnotemp = map_ratio.masked_where(~maptemp.getmask())
mapmanual_temp = maps.merge_maps(mapmanual_temp, mapnotemp)


maptemp_filter = maps.filtermap(maptemp, typecut['Td'], cuts['Td'])
mapnotemp_filter = maptemp.masked_where(~maptemp_filter.getmask())
mapmanual_temp = maps.merge_maps(mapmanual_temp, mapnotemp_filter)

if fout['temperature'] :
    ok = maptemp_filter.save_fitsfile(
        fname=fout['temperature'], hdr_type='tdust', oldheader=map850.header,
        append=False, overwrite=True)

print("   ...done")


#
# flux was in mJy, all still in SI
#
print("  >> convert flux to SI...")

mapS_850 = mapdblf_cl850.cmult(flux_factor)
print("   ...done")


print("  >> calculating masses...")

mapmass = calc_mapmass(mapS_850, maptemp_filter, pr)



mapmass_filter = maps.filtermap(mapmass, typecut['M'], cuts['M'])

maptemp_filtermass = maptemp_filter.masked_where(mapmass_filter.getmask())

mapnot_filtermass = maptemp_filter.masked_where(~maptemp_filtermass.getmask())
mapmanual_temp = maps.merge_maps(mapmanual_temp, mapnot_filtermass)

print("   ...done")


print("  >>\n  >> processing fixed dust temperature pixels...")

mapf850_notemp = mapclumpshi850.masked_where(mapmanual_temp.getmask())

mapmanual_Tdust = maps.full_like(mapf850_notemp,
                                 (defaults['Td'], defaults['varTd']))

mapS850_notemp = mapf850_notemp.cmult(flux_factor)


mapmass_notemp = calc_mapmass(mapS850_notemp, mapmanual_Tdust, pr)


#print(np.nansum(mapmass_notemp.data[0]), "+-",
#      ma.sqrt(np.nansum(mapmass_notemp.data[1])) )

#print(np.nansum(mapmass_filter.data[0]), "+-",
#      ma.sqrt(np.nansum(mapmass_filter.data[1])) )

mapmass_total = maps.merge_maps(mapmass_filter, mapmass_notemp)

if fout['mass'] :
    ok = mapmass_total.save_fitsfile(
        fname=fout['mass'], hdr_type='mass', oldheader=map850.header,
        append=False, overwrite=True)

print("   ...done")


if args.clumps :
    print("  >> clump calculations")

    mapclumpcat = cl.ClumpCatalog.from_calcphys(
        idxs=clump_idxs, fluxes=[mapclumpshi850,mapclumpshi450,mapclumps450],
        temps=[maptemp_filtermass], mass=[mapmass_filter, mapmass_total],
        params=pr)

    if fout['clumptable']:
        mapclumpcat.save_catalog(fout['clumptable'], overwrite=True,
                                 ctype='phys')

        if args.print_table: 
            mapclumpcat.print_catalog(ctype='phys')

print(" ++ ok")
