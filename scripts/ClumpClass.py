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
from astropy.io import fits



class Clump(object):

    phys_names = [
        'id', 'npix', 'nTpix', 'area', 'eff_radius', 'deconv_radius',
        'flux850', 'flux450hi', 'flux450',
        'mass', 'var_mass', 'masstot', 'var_masstot',
        'N', 'varN', 'maxT', 'minT']

    phys_formats = [
        'i4', 'i4' , 'i4', 'f8', 'f8', 'f8',
        'f8', 'f8', 'f8',
        'f8', 'f8', 'f8', 'f8',
        'f8', 'f8', 'f4', 'f4']
    
    phys_units = [
        '', '', '', 'arcsec^2', 'arcsec', 'arcsec',
        'mJy', 'mJy', 'mJy',
        'Msol', 'Msol', 'Msol', 'Msol',
        'cm-2', 'cm-2', 'K', 'K' ]

    phys_prfmt = [
        '2d', '4d', '4d', '7.1f', '8.2f', '8.2f',
        '7.3f', '7.3f', '7.3f',
        '8.2f', '6.2f', '8.2f', '6.2f',
        '9.3e', '9.3e',
        '6.1f', '6.1f' ]

    
    findclumps_names = [
        'PIDENT', 'Peak1', 'Peak2', 'Cen1', 'Cen2', 'Size1', 'Size2', 'Sum',
        'Peak', 'Volume', 'Shape' ]

    findclumps_formats = [
        'i4', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8',
        'f8', 'f8', 'U250' ]

    findclumps_units = [
        '', 'deg', 'deg', 'deg', 'deg', 'arcsec', 'arcsec', 'mJy/beam',
        'mJy/beam', 'arcsec.arcsec', '' ]


    list_names = []
    list_names.extend(findclumps_names)
    list_names.extend(phys_names)

    list_formats = []
    list_formats.extend(findclumps_formats)
    list_formats.extend(phys_formats)
 
    ddtype = { 'names' : list_names, 'formats' : list_formats}
    dtype_phys = { 'names' : phys_names, 'formats' : phys_formats}
   

    
    def __init__(self, id, record=None, calcphys=False, **kwargs):

        self.id = id

        self.record = np.empty(1, dtype=self.ddtype)

        if record :
            if not calcphys:
                self.fill_clump(record)
            else :
                print("todo")
                
        elif calcphys :
            self.calc_phys(id, **kwargs)



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
                    
            if cl_td.data[0].count() > 0 :
                self.record['maxT'] = np.nanmax(cl_td.data[0])
                self.record['minT'] = np.nanmin(cl_td.data[0])
            else:
                self.record['maxT'] = -99
                self.record['minT'] = -99



    def print_clump(self, ctype='phys') :
        """Prints the information contained in a clump.

        ctype: ....
        """

        out = ""
        ct = 0
        for ff in self.phys_names :
            
            prfmt = self.phys_prfmt[ct]
            strfmt = '{0} {1:'+str(prfmt)+'}'
            out = strfmt.format(out, (self.record[ff])[0])

            ct += 1

        return out

    

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

        solangle = 2. * np.pi * (1. - np.cos(2. * np.sqrt(arad/np.pi)))
        
        fact = const.M_sun / mu / mH / dcm / dcm / solangle
        f = fact * u.cm * u.cm
        
        col = mass * f
        var_col = var_mass * f *f

        return col, var_col


    
    
class ClumpCatalog (object):


    def __init__(self, ncl=0, clumps=None, findclump_file='', phys_file=''):
        """Initialization of the ClumpCatalog object."""

        if clumps :
            self.clumps = clumps
        else :
            self.clumps = []
            
        self.ncl = ncl
        self.fclfile = findclump_file
        self.phyfile = phys_file
        

        
    @classmethod
    def fromfile(cls, fname):
        """Creates a ClumpCatalog object from one or two files

        fname has to be a list
        """
        
        new = cls()
        
        for ff in fname :
            new.fillfrom_file(ff)

        return new



    @classmethod
    def from_calcphys(cls, **kwargs):
        """Generates a ClumpCatalog object from the calculation of
           physical parameters
        """

        new = cls()

        new.fillfrom_calcphys(**kwargs)
        
        return new

    
    
    def fillfrom_file(self, fname) :
        """Fill a ClassCatalog from the contents of file fname."""
            
        print("  >> reading", fname)
        
        with fits.open(fname) as hdul:
            units = np.shape(hdul)[0]
            if units == 2 :
                data = hdul[1].data
                
            elif units == 1 :
                print(" >> WARNING: something strange is going on.",
                      "There is only one HDU present; we'll try to open it")
                data = hdul[0].data
            elif units > 2:
                print(" >> More than 2 HDUs found. We open #1 only") 
                data = hdul[1].data

        if 'PIDENT' in data.dtype.names:
            self.fclfile = fname

        elif 'id' in data.dtype.names:
            self.phyfile = fname

        self.fillfrom_rec(data)

        #self.ncl = np.shape(self.clumps)[0]

        
        
    def id_isincatalog(self, i) :
        """Check if clump with id equal to i is in the catalog."""
        
        for cl in self.clumps :
            if i == cl.id :
                return cl
            
        return False

    
        
    def fillfrom_rec(self, rec):
        """Fill the a ClumpCatalog object from a FITS_rec object."""

        ncl = np.shape(rec)[0]
        newclumps = 0
        
        for cl in range(ncl) :
            if 'PIDENT' in rec.dtype.names:
                read_cl = rec['PIDENT'][cl] - 1
            elif 'id' in rec.dtype.names:
                read_cl = rec['id'][cl] - 1

            if not self.ncl :
                self.clumps.append(Clump(cl, record=rec[cl]))
                newclumps += 1

            else :
                clump = self.id_isincatalog(read_cl)
                if clump :
                    clump.fill_clump(rec[cl])
                else :
                    print(" >> WARNING: found new records")
                    self.clumps.append(Clump(self.ncl, record=rec[cl]))
                    newclumps += 1

        self.ncl += newclumps


        
    def fillfrom_calcphys(self, idxs=None, **kwargs):
        """Fill a ClumpCatalog Object from calcphys.

        fluxes, temps, and mass are suposed to be lists of type Maps()
        objects
        """
        
        n_clumps = np.int(np.nanmax(idxs))
        
        newclumps = 0
        
        for cl in range(1, n_clumps+1):
            id = cl - 1
            
            if not self.ncl :
                self.clumps.append(Clump(id, calcphys=True, idxs=idxs,
                                         **kwargs))
                newclumps += 1

            else :
                clump = self.id_isincatalog(cl)
                if clump :
                    clump.fill_clump(rec, calcphys=True, **kwargs)

        self.ncl += newclumps
        
        #print(np.shape(self.clumps))
        #for x in self.clumps:
        #    print(x.record)
                    

            
    def getview(self, ctype):
        """Returns a selected view of the catalog."""
        
        if ctype == 'all' :
            dt = Clump.ddtype
        elif ctype == 'phys' :
            dt = Clump.dtype_phys

        reduc = np.empty(self.ncl, dtype=dt)

        i = 0
        for clump in self.clumps:
            for fl in dt['names'] :
                reduc[fl][i] = (clump.record[fl])[0]
                
            i += 1
            
        return reduc

    
    
    def save_catalog(self, fname, overwrite=False, ctype='all') :
        """Save a fits table with the clump parameters.

        ctype can select a reduced view of the catalog
        """
        
        print(" >>> Saving fits table ", fname," ...")

        t = fits.BinTableHDU.from_columns(self.getview(ctype))
        t.writeto(fname, overwrite=overwrite)

        print(" >>> ...done")


        
    def print_catalog(self, ctype='phys', **kwargs) :
        """Prints the catalog

        ctype:
           - 'all' the fields
           - 'phys' calc_phys fields
           - 'custom' custom set of fields (to finish)
        """
        ##print header
        ##print units
        ## print clumps
        for clump in self.clumps:
            str_out = clump.print_clump(ctype, **kwargs)

            print(str_out)

