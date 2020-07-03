#!/usr/bin/env python3
##
## ClumpClass.py
## O. Morata 2020
##
""" ClumpClass.py

    Definition of the class Clump 

"""

import numpy as np
import numpy.ma as ma
import astropy.units as u
from astropy import constants as const



class Clump(object):

    phys_names = [
        'id', 'npix', 'nTpix', 'area', 'eff_radius', 'deconv_radius',
        'flux850', 'flux450hi', 'flux450',
        'mass', 'var_mass', 'masstot', 'var_masstot',
        'N', 'varN', 'maxT', 'minT', 'meanT', 'medianT', 'weightT']

    phys_formats = [
        'i4', 'i4' , 'i4', 'f8', 'f8', 'f8',
        'f8', 'f8', 'f8',
        'f8', 'f8', 'f8', 'f8',
        'f8', 'f8', 'f4', 'f4', 'f4', 'f4', 'f4']
    
    phys_units = [
        '', '', '', 'arcsec^2', 'arcsec', 'arcsec',
        'mJy', 'mJy', 'mJy',
        'Msol', 'Msol', 'Msol', 'Msol',
        'cm-2', 'cm-2', 'K', 'K', 'K', 'K', 'K' ]

    phys_prfmt = [
        '2d', '4d', '4d', '7.1f', '8.2f', '8.2f',
        '7.3f', '7.3f', '7.3f',
        '8.2f', '6.2f', '8.2f', '6.2f',
        '9.3e', '9.3e',
        '6.1f', '6.1f', '6.1f', '6.1f', '6.1f' ]

    
    findclumps_names = [
        'PIDENT', 'Peak1', 'Peak2', 'Cen1', 'Cen2', 'Size1', 'Size2', 'Sum',
        'Peak', 'Volume', 'Shape' ]

    findclumps_formats = [
        'i4', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8',
        'f8', 'f8', 'U500' ]

    findclumps_units = [
        '', 'deg', 'deg', 'deg', 'deg', 'arcsec', 'arcsec', 'mJy/beam',
        'mJy/beam', 'arcsec.arcsec', '' ]

    findclumps_prfmt = [
        '3d', '11.6f', '11.6f', '11.6f', '11.6f',  '11.6f', '11.6f', '14.6f',
        '12.6f', '14.6f', 's'
        ]


    list_names = []
    list_names.extend(findclumps_names)
    list_names.extend(phys_names)

    list_formats = []
    list_formats.extend(findclumps_formats)
    list_formats.extend(phys_formats)
 
    ddtype = { 'names' : list_names, 'formats' : list_formats}
    dtype_phys = { 'names' : phys_names, 'formats' : phys_formats}
   

    
    def __init__(self, id, record=None):

        self.id = id

        self.record = np.empty(1, dtype=self.ddtype)

        if record :
            self.fill_clump(record)



    @classmethod
    def from_calcphys(cls, id, **kwargs):

        new = cls(id)
        new.calc_phys(id, **kwargs)

        return new

    

    def fill_clump(self, rec):
        """Fills the fields of a clump from a FITS_rec object."""

        for f in self.list_names :
            try:
                self.record[f] = rec[f]
            except KeyError:
                pass


            
    def calc_phys(self, id, idxs='', fluxes='', temps='', mass='', params=''):
        """Fill a Clump Object from calcphys.

        fluxes, temps, and mass are suposed to be lists of type Maps()
        objects
        """
        
        fluxflds = ['flux850', 'flux450hi', 'flux450'] 
        massflds = ['mass', 'var_mass', 'masstot', 'var_masstot']

        cl = id + 1
        mskcl = ma.masked_not_equal(idxs, cl)
        npix = mskcl.count()

        self.record['id'] = cl
        self.record['npix'] = npix

        area, eff_radius, deconv_r = self.get_size(
            npix, params.pixsize, params.beam)
        
        self.record['area'] = area / u.arcsec / u.arcsec
        self.record['eff_radius'] = eff_radius / u.arcsec
        self.record['deconv_radius'] = deconv_r / u.arcsec

        numflux = np.shape(fluxes)[0]
        for flx in range(numflux):
            f = fluxes[flx].masked_where(ma.getmask(mskcl))

            S = f.cmult(params.flux_factor)
            S = S.filled(np.nan)

            field = fluxflds[flx]
            self.record[field] = np.nansum(S.data[0]) / 1e-26

                
        numasses = np.shape(mass)[0]
        for mss in range(numasses) :

            cl_m = mass[mss].masked_where(ma.getmask(mskcl))

            if mss == 0 :
                self.record['nTpix'] = cl_m.data[0].count()

            cl_m = cl_m.filled(np.nan)
                

            m = np.nansum(cl_m.data[0])
            var_m = np.nansum(cl_m.data[1])
            self.record[massflds[mss*2]] = m
            self.record[massflds[mss*2+1]] = var_m

            if mss == 1 :
                # column density ad-hoc
                #
                clump_N, clump_varN = self.columndensity(
                    m, var_m, area, params.d, params.mu, params.mH)

                self.record['N'] = clump_N
                self.record['varN'] = clump_varN

                    
        numtemps = np.shape(temps)[0]
        for t in range(numtemps) :
            cl_td = temps[t].masked_where(ma.getmask(mskcl))

            f = fluxes[0].masked_where(ma.getmask(mskcl))


            if cl_td.data[0].count() > 0 :
                self.record['maxT'] = np.nanmax(cl_td.data[0])
                self.record['minT'] = np.nanmin(cl_td.data[0])
                self.record['meanT'] = cl_td.data[0].mean()
                self.record['weightT'] = ma.average(cl_td.data[0],
                                                    weights=f.data[0])
                self.record['medianT'] = ma.median(cl_td.data[0])
            else:
                self.record['maxT'] = -99
                self.record['minT'] = -99



    def print_clump(self, ctype='phys', fields=None) :
        """Prints the information contained in a clump.

        ctype: ....
        """

        if ctype == 'phys' :
            header, out = self.extract_values(
                self.record, names=self.phys_names, fields=fields,
                print_formats=self.phys_prfmt)

        elif ctype == 'findclumps' :
            header, out = self.extract_values(
                self.record, names=self.findclumps_names, fields=fields,
                print_formats=self.findclumps_prfmt)

        elif ctype == 'all' :
            all_prfmt = self.findclumps_prfmt.copy()
            all_prfmt.extend(self.phys_prfmt)

            header, out = self.extract_values(
                self.record, names=self.list_names, fields=fields,
                print_formats=all_prfmt)

        return header, out



    @staticmethod
    def extract_values(record, names=None, fields=None, print_formats=None):

        header = ""
        out = ""
        ct = 0

        if fields == None :
            fields = names

        for ff in names :
            if ff not in fields :
                ct += 1
                continue

            prfmt = print_formats[ct]
            strfmt = '{0} {1:'+str(prfmt)+'}'
            out = strfmt.format(out, (record[ff])[0])
            header = '{0}  {1}'.format(header, ff)
            
            ct += 1

        return header, out



    @staticmethod
    def get_size(n, pixsize, beamsize) :
        """Calculate some geometrical parameters of the clump."""

        area = n * pixsize * pixsize
        ef_rad = np.sqrt(area / np.pi)
        dec_r = np.sqrt(4 * ef_rad * ef_rad -
                        (np.pi / 4. / np.log(2.) )* beamsize * beamsize) * 0.5

        dec_a = np.pi * dec_r * dec_r
        
        return area, ef_rad, dec_r


    
    @staticmethod
    def columndensity(mass, var_mass, area, d, mu, mH ) :
        """Calculate column density from a mass and an area

        mass in solar masses, area in arcsec^2, d in pc.
        Returns the column density in cm^-2 and the variance of the column
        density in cm^-4
        """

        arad = area.to(u.radian * u.radian)
        dcm = d.to(u.cm)

        solangle = 2. * np.pi * (1. - np.cos(np.sqrt(arad/np.pi)))
        
        fact = const.M_sun / mu / mH / dcm / dcm / solangle
        f = fact * u.cm * u.cm
        
        col = mass * f
        var_col = var_mass * f *f

        return col, var_col


    
    
