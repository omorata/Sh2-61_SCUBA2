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

import logging

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
#        mpmass = mapmass_h2_thin(flux_mask, temp, (par.d).to(u.m), knu, bnu,
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



#def show_values(msk_arr, txt):
#    """Show some values of the msk_arr array."""
#    
#    print("values for" , txt)
#    print(msk_arr[324:326,300:325])
#    print(txt, msk_arr.count())
#    print("  max:", np.nanmax(msk_arr))
#    print("  min:", np.nanmin(msk_arr))



#def flux450(f850, td, h850, h450, pre) :
#    """Calculate the flux at 450micron from the flux at 850 micron and a
#        given dust temperature
#    """
#    
#    x = np.exp( 1. / td)
#    f = pre * f850 * (x**h850 - 1) / (x**h450 - 1)
#    return f



def get_col_factor(par, solangle):

    dist = (pr.d).to(u.cm)
    col_factor = const.M_sun / pr.mu / pr.mH / dist / dist / solangle
    col_factor *= u.cm * u.cm * u.rad * u.rad
    return col_factor



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
        '-l', dest='logfile',  help='log file',
        default='', metavar='FILE')

    parser.add_argument(
        '--print', dest='print_table', action='store_true', default=False,
        help='print table')
    
    parser.add_argument(
        '-w', dest='wkdir',  help='working directory',
        default='.', metavar='DIR')

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
    """ Read key, value pair from section in cfg file"""

    if not section in cfg:
        print("\n  WARNING: section \'", section, "\' not found\n")
        return None

    
    if key in cfg[section] :
        return cfg[section][key]
    else :
        if status == 'required' :
            print(" >> ERROR: missing required key", key)
            sys.exit(1)
        else :
          #  return ''
            return None

        
        
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



def set_outdefaults(opt):
    """ Sets the default values for opt"""
    
    if opt['tdust'] == None:
        opt['tdust'] = 'all'
    if opt['mass'] == None:
        opt['mass'] = 'all'
    if opt['N'] == None:
        opt['N'] = 'all'
    if opt['tau_opt'] == None:
        opt['tau_opt'] = 'thin'
    if opt['tau'] == None:
        opt['tau'] = 'all'


        
##-- End of functions --------------------------------------------------

range1 = (324,326,None,300,325,None)
range2 = (324,327,None,299,324,None)

print(" ++ Start")


args = read_command_line()

wdir = args.wkdir+'/'

if args.logfile:
    logg=True
else:
    logg=False

cnfg = read_configfile(args.cfgfile)

fname = get_values(cnfg, 'infiles', status='required',
                   names = ['data850', 'data450', 'snr850', 'snr450',
                            'clumpdef'])

fout = get_values(cnfg, 'outfiles', 
                  names = ['ratio', 'temperature', 'mass', 'N', 'clumptable',
                           'tempgood'])
#print(fout)

pr = par.Param.from_cfgfile(cnfg, 'phys_params')


defaults = get_values(cnfg, 'data_params', type='float',
                      names=['ini_temp','manual_Tdust','manual_varTdust'],
                      altnames=['iniTd', 'Td', 'varTd'])

typecut = get_values(cnfg, 'data_params', names=['type_cutTd', 'type_cutM'],
                     altnames=['Td', 'M'])

cuts = get_values(cnfg, 'data_params', names=['snr_450', 'cut_Td', 'cut_M'],
                  altnames=['snr450', 'Td', 'M'], type='float')

out_opts = get_values(cnfg, 'output_options',
                      names=['tdust', 'mass', 'N', 'tau_opt', 'tau'],
                      type='str')

set_outdefaults(out_opts)
# show settings for output

if logg: 
    logging.basicConfig(filename=args.logfile, filemode='w',
                        level=logging.INFO)
    logging.info('Started')

# calculate some variables
#

pixsize_rad = pr.pixsize.to(u.radian)
pixsolangle = pixsize_rad * pixsize_rad

pre = (850. / 450.)**(3.+ pr.beta)




print("  >> reading input files...")

if logg:
    logging.info('+ reading input files')

map850 = maps.Map.from_fitsfile(wdir+fname['data850'], name="850micron")
mapsnr = maps.Map.from_fitsfile(wdir+fname['snr850'], name="SNR 850micron")
map450 = maps.Map.from_fitsfile(wdir+fname['data450'], name="450micron")
mapsnr450 = maps.Map.from_fitsfile(wdir+fname['snr450'], name="SNR 450micron")

#header850 = map850.header


pr.get_flux_factor(map850)

    
print("  >> reading clump mask...")
if logg:
    logging.info('+ reading clump mask')

with fits.open(wdir+fname['clumpdef']) as hdumask:
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

#print(mapdblf_cl850.count())
#print(mapdblf_cl850.show(range1))
#print(mapsf_cl850.count())
#print(mapsf_cl850.show(range1))



# definition of array to hold pixels where WE fix Tdust
#

maskmanual = mapsf_cl850.copy()
#print("maskmanual")
#print(maskmanual.count())
#print(maskmanual.show(range1))


print("  >> calculating flux ratios...")
if logg:
    logging.info('+ calculating flux ratios')

map_ratio = maps.divide(mapclumpshi450, mapdblf_cl850)

#print("map_ratio")
#print(map_ratio.count())
#print(map_ratio.show(range1))


if fout['ratio']:
    ok = map_ratio.save_fitsfile(
        fname=wdir+fout['ratio'], hdr_type='fluxratio', oldheader=map850.header,
        append=False, overwrite=True)

print("  ...done")


print("  >> calculating temperatures...")

ini_array = np.full_like(map_ratio.data[0], defaults['iniTd'])

with np.errstate(invalid='ignore'):
    maptemp = get_maptemperature(map_ratio, ini_array, pre, pr)
    
mapnotemp = map_ratio.masked_where(~maptemp.getmask())

maskmanual = maps.merge_maps(maskmanual, mapnotemp)

#print("maskmanual")
#print(maskmanual.count())
#print(maskmanual.show(range1))


maptemp_filter = maps.filtermap(maptemp, typecut['Td'], cuts['Td'])
mapnotemp_filter = maptemp.masked_where(~maptemp_filter.getmask())
maskmanual = maps.merge_maps(maskmanual, mapnotemp_filter)




#print("  >>\n  >> processing fixed dust temperature pixels...")

#mapf850_notemp = mapclumpshi850.masked_where(maskmanual.getmask())

#mapmanual_Tdust = maps.full_like(mapf850_notemp,
#                                 (defaults['Td'], defaults['varTd']))

#temptotal = maps.merge_maps(maptemp_filter, mapmanual_Tdust)


#if fout['temperature'] :
#    ok =temptotal.save_fitsfile(
#        fname=wdir+fout['temperature'], hdr_type='tdust',
#        oldheader=map850.header, append=False, overwrite=True)

#print("   ...done")


#
# flux was in mJy, all still in SI
#
#print("  >> convert flux to SI...")

mapS_850 = mapdblf_cl850.cmult(pr.flux_factor)
#mapS850_notemp = mapf850_notemp.cmult(pr.flux_factor)

#print("   ...done")



if out_opts['tau_opt'] == 'thin' :

    print("  >> calculating masses, optically thin approximation...")

    mass_thin = calc_mapmass(mapS_850, maptemp_filter, pr)

    mass_calc = maps.filtermap(mass_thin, typecut['M'], cuts['M'])

    Tdust_calc = maptemp_filter.masked_where(mass_calc.getmask())

    mapnot_filtermass = maptemp_filter.masked_where(~Tdust_calc.getmask())
    maskmanual = maps.merge_maps(maskmanual, mapnot_filtermass)


    f850_manual = mapclumpshi850.masked_where(maskmanual.getmask())

    # convert fluxes to SI
    #
    S850_manual = f850_manual.cmult(pr.flux_factor)

    Tdust_manual = maps.full_like(f850_manual, (defaults['Td'],
                                                defaults['varTd']))

    mass_manual = calc_mapmass(S850_manual, Tdust_manual, pr)
    

    mass_total = maps.merge_maps(mass_calc, mass_manual)

    print("   ...done")



    # calculate array of column densities
    #
    column_factor = get_col_factor(pr, pixsolangle)

    coldens_total = mass_total.copy()
    coldens_total = coldens_total.cmult(column_factor)


        
elif out_opts['tau_opt'] == 'thick' :

    mass_tau, taus, coldns_tau = calc_mapmass(mapS_850, maptemp_filter, pr,
                                              calc_type='tau',
                                              solangle=pixsolangle)

    print("+++= mass_tau", mass_tau.count())
    print("+++= taus", taus.count())
    print("+++= coldns_tau", coldns_tau.count())
    
    #opac calc mass, tau, coldns manual_Tdust
    #
    # save all files at the end
    # todo: calculate variances of mass_tau, taus & coldns_tau

    
    #mass_total = maps.merge_maps(mass_tau, mapmass_notemp)
    #print("+++= mass_total", mass_total.count())

    if fout['mass'] :
        ok = mass_total.save_fitsfile(
            fname=wdir+fout['mass'], hdr_type='mass', oldheader=map850.header,
            append=False, overwrite=True)

    print("   ...done")

    
else:
    print(" >>> ERROR: the value of tau_opt", out_opts['tau_opt'],
          "is unknown")
    sys.exit(1)
    


Tdust_total = maps.merge_maps(Tdust_calc, Tdust_manual)


# save only calculated temperatures, if requested
#
if fout['tempgood'] :
    ok = maptemp_filter.save_fitsfile(
        fname=wdir+fout['tempgood'], hdr_type='tdust',
        oldheader=map850.header, append=False, overwrite=True)

# save temperature
#
if fout['temperature'] :
    ok =Tdust_total.save_fitsfile(
        fname=wdir+fout['temperature'], hdr_type='tdust',
        oldheader=map850.header, append=False, overwrite=True)

print("   ...done")



if fout['mass'] :
    ok = mass_total.save_fitsfile(
        fname=wdir+fout['mass'], hdr_type='mass', oldheader=map850.header,
        append=False, overwrite=True)

print("   ...done")



if fout['N']:
    ok = coldens_total.save_fitsfile(
        fname=wdir+fout['N'], hdr_type='column', oldheader=map850.header,
        append=False, overwrite=True)


    
if args.clumps :
    print("  >> clump calculations")

    mapclumpcat = cl.ClumpCatalog.from_calcphys(
        idxs=clump_idxs, fluxes=[mapclumpshi850,mapclumpshi450,mapclumps450],
        temps=[Tdust_calc], mass=[mass_calc, mass_total],
        params=pr)

    if fout['clumptable']:
        mapclumpcat.save_catalog(wdir+fout['clumptable'], overwrite=True,
                                 ctype='phys')

        if args.print_table: 
            mapclumpcat.print_catalog(ctype='phys')

print(" ++ ok")
if logg:
    logging.info('+ ok')
