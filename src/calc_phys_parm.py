#!/usr/bin/env python3
#

#import sys

import numpy as np
import numpy.ma as ma

from astropy.io import fits
from astropy import units as u

import argparse as ap
import yaml

import logging

# definition of the custom classes
#
import Param as par
import ClumpCatalog as cl
import MapClass as maps

import temperature as t
import mass as m

##-- Functions ---------------------------------------------------------


##def show_values(msk_arr, txt):
##    """Show some values of the msk_arr array."""
##    
##    print("values for" , txt)
##    print(msk_arr[324:326,300:325])
##    print(txt, msk_arr.count())
##    print("  max:", np.nanmax(msk_arr))
##    print("  min:", np.nanmin(msk_arr))



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


        
def print_outputsettings(opt):
    """ prints the output settings"""

    print("  ++ Output settings:")
    print("     -----------------")
    print("       Tdust: ", opt['tdust'])
    print("        mass: ", opt['mass'])
    print("           N: ", opt['N'])
    print("     tau_opt: ", opt['tau_opt'])
    if opt['tau_opt'] == 'thick' :
        print("         tau: ", opt['tau'])
    print("     -----------------")

        
##-- End of functions --------------------------------------------------


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

print_outputsettings(out_opts)


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
snr850 = maps.Map.from_fitsfile(wdir+fname['snr850'], name="SNR 850micron")
map450 = maps.Map.from_fitsfile(wdir+fname['data450'], name="450micron")
snr450 = maps.Map.from_fitsfile(wdir+fname['snr450'], name="SNR 450micron")

#header850 = map850.header


# set the pixel flux_factor for the main map
#
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

clumps_hi850 = map850.masked_where(ma.getmask(inclumps))


high_snr450 = snr450.masked_less(cuts['snr450'])
clumps_high_snr450 = ma.masked_where(high_snr450.getmask(), inclumps)
clumps_hi450 = map450.masked_where(ma.getmask(clumps_high_snr450))
clumps_450 = map450.masked_where(ma.getmask(inclumps))


# calculations are done on pixels with two detected frequencies
#
clumps_hihi850 = clumps_hi850.masked_where(ma.getmask(clumps_high_snr450))

# all the single detected frequency pixels are manual by default
#
clumps_lowhi850_idx = ma.masked_where(~clumps_hihi850.getmask(), inclumps)
maskmanual = map850.masked_where(ma.getmask(clumps_lowhi850_idx))



print("  >> calculating flux ratios...")
if logg:
    logging.info('+ calculating flux ratios')

flux_ratio = maps.divide(clumps_hi450, clumps_hihi850)



if fout['ratio']:
    ok = flux_ratio.save_fitsfile(
        fname=wdir+fout['ratio'], hdr_type='fluxratio', oldheader=map850.header,
        append=False, overwrite=True)

print("  ...done")


print("  >> calculating temperatures...")

Tdust_init = np.full_like(flux_ratio.data[0], defaults['iniTd'])

with np.errstate(invalid='ignore'):
    Td = t.get_maptemperature(flux_ratio, Tdust_init, pre, pr)
    
noTd = flux_ratio.masked_where(~Td.getmask())

maskmanual = maps.merge_maps(maskmanual, noTd)


# apply filter to Td
#
Tdust_filter = maps.filtermap(Td, typecut['Td'], cuts['Td'])
noTdust_filter = Td.masked_where(~Tdust_filter.getmask())
maskmanual = maps.merge_maps(maskmanual, noTdust_filter)




# flux was in mJy, all still in SI
#
mapS_850 = clumps_hihi850.cmult(pr.flux_factor)



if out_opts['tau_opt'] == 'thin' :

    print("  >> calculating masses, optically thin approximation...")

    mass_thin = m.calc_mapmass(mapS_850, Tdust_filter, pr)

    mass_calc = maps.filtermap(mass_thin, typecut['M'], cuts['M'])

    Tdust_calc = Tdust_filter.masked_where(mass_calc.getmask())

    nomass_calc = Tdust_filter.masked_where(~Tdust_calc.getmask())
    maskmanual = maps.merge_maps(maskmanual, nomass_calc)


    f850_manual = clumps_hi850.masked_where(maskmanual.getmask())

    # convert fluxes to SI
    #
    S850_manual = f850_manual.cmult(pr.flux_factor)

    Tdust_manual = maps.full_like(f850_manual, (defaults['Td'],
                                                defaults['varTd']))

    mass_manual = m.calc_mapmass(S850_manual, Tdust_manual, pr)
    

    mass_total = maps.merge_maps(mass_calc, mass_manual)

    print("   ...done")



    # calculate arrays of column densities
    #
    pixel_column_factor = m.get_col_factor(pr, pixsolangle)

    coldens_calc = mass_calc.copy()
    coldens_calc = coldens_calc.cmult(pixel_column_factor)

    coldens_manual = mass_manual.copy()
    coldens_manual = coldens_manual.cmult(pixel_column_factor)
    
    coldens_total = mass_total.copy()
    coldens_total = coldens_total.cmult(pixel_column_factor)


        
elif out_opts['tau_opt'] == 'thick' :

    mass_tau, taus, coldns_tau = m.calc_mapmass(mapS_850, Tdust_filter, pr,
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


# initialize output map variables
Tdust_out = None
mass_out = None
coldens_out = None

if fout['temperature'] :

    Tdust_out = maps.get_outmap([Tdust_total, Tdust_calc],
                                ['all', 'calc_only'], out_opts['tdust'])
        
        
    ok =Tdust_out.save_fitsfile(
        fname=wdir+fout['temperature'], hdr_type='tdust',
        oldheader=map850.header, append=False, overwrite=True)

print("   ...done")



if fout['mass'] :
    mass_out = maps.get_outmap([mass_total, mass_calc], ['all', 'calc_only'],
                               out_opts['mass'])
        
    ok = mass_out.save_fitsfile(
        fname=wdir+fout['mass'], hdr_type='mass', oldheader=map850.header,
        append=False, overwrite=True)

print("   ...done")



if fout['N']:
    coldens_out = maps.get_outmap([coldens_total, coldens_calc],
                                  ['all', 'calc_only'], out_opts['N'])
        
    ok = coldens_out.save_fitsfile(
        fname=wdir+fout['N'], hdr_type='column', oldheader=map850.header,
        append=False, overwrite=True)


    
if args.clumps :
    print("  >> clump calculations")

    if Tdust_out is None:
        Tdust_out = maps.get_outmap([Tdust_total, Tdust_calc],
                                    ['all', 'calc_only'], out_opts['tdust'])

    if mass_out is None:
        mass_out = maps.get_outmap([mass_total, mass_calc],
                                   ['all', 'calc_only'], out_opts['mass'])

    if coldens_out is None:
        coldens_out = maps.get_outmap([coldens_total, coldens_calc],
                                      ['all', 'calc_only'], out_opts['N'])

    
    mapclumpcat = cl.ClumpCatalog.from_calcphys(
        idxs=clump_idxs, fluxes=[clumps_hi850,clumps_hi450,clumps_450],
        temps=[Tdust_out], mass=[mass_out], coldens=[coldens_out],
        tpix=Tdust_calc, params=pr)
    
    if fout['clumptable']:
        mapclumpcat.save_catalog(wdir+fout['clumptable'], overwrite=True,
                                 ctype='phys')

        if args.print_table: 
            mapclumpcat.print_catalog(ctype='phys')

print(" ++ ok")
if logg:
    logging.info('+ ok')
