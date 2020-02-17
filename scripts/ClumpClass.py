#!/usr/bin/env python3
##
## ClumpClass.py
## O. Morata 2020
##
""" ClumpClass

    Definition of the class Clump 

"""

import numpy as np
import numpy.ma as ma
import astropy.units as u
from astropy import constants as const
from astropy.io import fits


class Clump (object):

    def __init__(self, idxs='', fluxes='', temps='', mass='',params='',
                 names='') :
        """Initialization of the Clump object."""

        self.n_clumps = np.int(np.nanmax(idxs))
        self.idxs = idxs
        self.fluxes = fluxes
        self.temps = temps
        self.mass = mass

        self.names= names
            
        self.params = params

        self.rec = Clump.createRec(self)
        self.fluxfields = ['flux850', 'flux450hi', 'flux450'] 
        self.massfields = ['th_mass', 'varth_mass', 'totmass', 'vartotmass']
        
        
        
    def createRec(self):
        """Creates an empty structured array."""
        
        list_names = ['id', 'npix', 'ngpix', 'area', 'eff_radius',
                      'deconv_radius', 'flux850', 'flux450hi', 'flux450',
                      'th_mass', 'varth_mass', 'totmass', 'vartotmass',
                      'N', 'varN', 'maxT', 'minT']

        list_formats = ['i4', 'i4' , 'i4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4',
                        'f4', 'f4', 'f4', 'f4', 'f8', 'f8', 'f4', 'f4']
        
        ddtype = { 'names' : list_names,
                   'formats' : list_formats}

        record = np.empty(self.n_clumps, dtype=ddtype)

        return record

    
        
    def make_table(self) :
        """Builds the catalogue table."""

        for clump in range(1, self.n_clumps+1):
            mskcl = ma.masked_not_equal(self.idxs, clump)
            npix = mskcl.count()


            self.rec[clump-1]['id'] = clump
            self.rec[clump-1]['npix'] = npix


            area, eff_radius, deconv_r = Clump.get_size(npix,
                                                        self.params.pixsize,
                                                        self.params.beam)
            self.rec[clump-1]['area'] = area / u.arcsec / u.arcsec
            self.rec[clump-1]['eff_radius'] = eff_radius / u.arcsec
            self.rec[clump-1]['deconv_radius'] = deconv_r / u.arcsec

            numflux = np.shape(self.fluxes)[0]

            for flx in range(numflux):
                f = np.ma.masked_where(np.ma.getmask(mskcl), self.fluxes[flx])
                S = f * self.params.flux_factor
                S = S.filled(np.nan)
                field = self.fluxfields[flx]
                self.rec[clump-1][field] = np.nansum(S) / 1e-26


            numasses = np.shape(self.mass)[0]

            for mss in range(numasses) :

                cl_m_idx = np.ma.masked_where(
                    np.ma.getmask(self.mass[mss].data[0]),mskcl)
                
                cl_m = np.ma.masked_where(np.ma.getmask(cl_m_idx),
                                          self.mass[mss].data[0])
                cl_m = cl_m.filled(np.nan)
                
                cl_varm = np.ma.masked_where(np.ma.getmask(cl_m_idx),
                                              self.mass[mss].data[1])
                cl_varm = cl_varm.filled(np.nan)
                
                self.rec[clump-1][self.massfields[mss*2]] = np.nansum(cl_m)
                self.rec[clump-1][self.massfields[mss*2+1]] = np.nansum(cl_varm)

                if mss == 1 :
                    # column density ad-hoc
                    #
            
                    clump_N, clump_varN = Clump.columndensity(
                        np.nansum(cl_m), np.nansum(cl_varm), area,
                        self.params.d, self.params.mu, self.params.mH)

                    self.rec[clump-1]['N'] = clump_N
                    self.rec[clump-1]['varN'] = clump_varN

            numtemps = np.shape(self.temps)[0]
            for t in range(numtemps) :
                cl_td_idx = ma.masked_where(ma.getmask(self.temps[t]),
                                            mskcl)

                cl_td = ma.masked_where(ma.getmask(cl_td_idx), self.temps[t])
                    
                if cl_td.count() > 0 :
                    self.rec[clump-1]['maxT'] = np.nanmax(cl_td)
                    self.rec[clump-1]['minT'] = np.nanmin(cl_td)
                else:
                    self.rec[clump-1]['maxT'] = -99
                    self.rec[clump-1]['minT'] = -99


        return 0

        
                    
    def print_table(self):
        """Prints catalogue table."""
        
        for i in range(np.shape(self.rec)[0]) :

            str_out = '{0:2d} {1:4d}'.format(self.rec['id'][i],
                                             self.rec['npix'][i])

            str_out = ' {0} {1:7.1f} {2:8.2f} {3:8.2f}'.format(
                str_out, self.rec['area'][i],
                self.rec['eff_radius'][i],
                self.rec['deconv_radius'][i])

            
            numflux = np.shape(self.fluxes)[0]
            for flx in range(numflux):
                ff = self.fluxfields[flx]
                str_out = ' {0} {1:7.3f}'.format(str_out, self.rec[ff][i])  

                
            numasses = np.shape(self.mass)[0]
            for mss in range(numasses) :
                m_field = self.massfields[mss*2]
                varm_field = self.massfields[mss*2+1]
                str_out = '{0} {1:8.2f} ({2:6.2f})'.format(
                    str_out, self.rec[i][m_field],
                    np.sqrt(self.rec[i][varm_field]))



            str_out = '{0} {1:9.3e} ({2:9.3e})'.format(
                str_out, self.rec[i]['N'], np.sqrt(self.rec[i]['varN']))

            if self.rec[i]['maxT'] >0 :
                str_out = '{0} ///{1:6.1f} {2:6.1f}'.format(str_out,
                                                            self.rec[i]['maxT'],
                                                            self.rec[i]['minT'])
            else:
                str_out = '{0} ///{1:13s}'.format(str_out, " ")
                
            print(str_out)




    def save_fitstable(self, fname, overwrite=False) :
        """Test of saving a fits table."""

        print(" >>> Saving fits table ", fname," ...")
        t = fits.BinTableHDU.from_columns(self.rec)
        t.writeto(fname, overwrite=overwrite)
        print(" >>> ...done")
    
        
            
            
    @staticmethod
    def get_size(n, pixsize, beamsize) :
        """Calculate some geometrical parameters of the clump."""

        area = n* pixsize * pixsize
        ef_rad = np.sqrt(area / np.pi)
        dec_r = np.sqrt(4 * ef_rad * ef_rad -
                        (np.pi / 4. / np.log(2.) )* beamsize * beamsize) * 0.5
    
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
