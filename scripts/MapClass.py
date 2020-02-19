#!/usr/bin/env python3
##
## MapClass.py
## O. Morata 2020
##
""" MapClass

    Definition of the class Map 

"""

import numpy as np
import numpy.ma as ma
from astropy.io import fits
from datetime import datetime



class Map (object):

    def __init__(self, name='', filename='', data=[None]*2, header=[None]*2):
        """Initialization of the Map object."""

        self.name = name
        self.fname = filename
        self.data= data
        self.header = header

        #if data == '' :
        #    if self.fname != '' :
        #        self.data, self.header = Map.read_fitsfile(self)
        #    else :
        #        self.data = [None] * 2
        #        self.header = [None] * 2
        #else :
        #    self.data = data


        
    @classmethod
    def empty(cls) :
        name = ''
        fname = ''
        data = [None]*2
        header = [None] * 2

        new = cls(name=name, filename=fname, data=data, header=header)

        return new


    
    @classmethod
    def fromfitsfile(cls, fname, name=''):
        """Read input data and header from a FITS file."""

        if name == '' :
            name = fname
            
        print(" >> reading",name,"data...")
    
        with fits.open(fname) as hdu_data:
            data_info = [hdu_data[0].data, hdu_data[1].data]
            header_info = [hdu_data[0].header, hdu_data[1].header]


        new = cls(name=name, filename=fname, data=data_info,
                  header=header_info)

        return new
    
        
#    def read_fitsfile(self):
#        """Read input data and header from a FITS file."""
#
#        if self.name == '' :
#            name = self.fname
#        else :
#            name = self.name
#            
#        print(" >> reading",name,"data...")
#    
#        with fits.open(self.fname) as hdu_data:
#            data_info = [hdu_data[0].data, hdu_data[1].data]
#            header_info = [hdu_data[0].header, hdu_data[1].header]
#
#            return data_info, header_info


        
    def save_fitsfile(self, oldheader='', append=False, overwrite=False,
                      hdr_type='', fname=''):

        if oldheader :
            self.header = Map.modify_header(oldheader, hdr_type)
        else :
            return 3

        if fname != '' :
            self.fname = fname
    
        print(" >> saving", self.fname, "...")

        data = self.data[0].filled(np.nan)
        var = self.data[1].filled(np.nan)

        hdu_data = fits.PrimaryHDU(data, header=self.header[0])
        hdu_variance = fits.ImageHDU(var, header=self.header[1])

        hdulist = fits.HDUList([hdu_data, hdu_variance])
        hdulist.writeto(self.fname, overwrite=overwrite)
        
        return 0


    
    def modify_header(old, htype) :
        """ modify a FITS header according to predefined types."""

        # get current time
        #
        timenow = str(datetime.utcnow()).split('.')[0]
        timenow = timenow.replace(' ','T')
    
        if htype == "fluxratio" :
            old[0]['LABEL'] = 'Flux ratio'
            old[0]['BUNIT'] = ''
            old[0]['history'] = 'flux ratio Sh2-61'
            
            old[1]['LABEL'] = 'Flux ratio variance'
            old[1]['BUNIT'] = ''

        elif htype == "tdust" :
            old[0]['LABEL'] = 'Tdust'
            old[0]['BUNIT'] = 'K'
            old[0]['history'] = 'dust temperature Sh2-61'
            
            old[1]['LABEL'] = 'Tdust variance'
            old[1]['BUNIT'] = 'K^2'

        elif htype == "mass" :
            old[0]['LABEL'] = 'Mass H_2'
            old[0]['BUNIT'] = 'M_sol'
            old[0]['history'] = 'H_2 mass Sh2-61'
        
            old[1]['LABEL'] = 'Mass H_2 variance'
            old[1]['BUNIT'] = '(M_sol)^2'

        elif htype == "column" :
            old[0]['LABEL'] = 'N(H_2)'
            old[0]['BUNIT'] = 'cm^-2'
            old[0]['history'] = 'H_2 column density Sh2-61'
        
            old[1]['LABEL'] = 'N(H_2) variance'
            old[1]['BUNIT'] = 'cm^-4'


            for h in range(2) :
                old[h]['DATE'] = timenow
                old[h]['ORIGIN'] = 'MapClass.py'
        
    
        return old



    def masked_less(self, val):
        new = Map.empty()

        new.data[0] = ma.masked_less(self.data[0], val)
        new.data[1] = ma.masked_less(self.data[1], val)

        return new



    def masked_greater(self, val):
        new = Map.empty()

        new.data[0] = ma.masked_greater(self.data[0], val)
        new.data[1] = ma.masked_greater(self.data[1], val)

        return new


    
    def getmask(self) :
        mask = ma.getmask(self.data[0])
        return mask


    
    def masked_where(self, mask) :
        new = Map.empty()
        
        new.data[0] = ma.masked_where(mask, self.data[0])
        new.data[1] = ma.masked_where(mask, self.data[1])

        return new


    
    def copy(self) :
        new = Map.empty()

        new.data[0] = ma.copy(self.data[0])
        new.data[1] = ma.copy(self.data[1])
        new.header = self.header.copy()

        return new


    def mult(self, factor):

        new = Map.empty()
        new.data[0] = self.data[0] * factor
        new.data[1] = self.data[1] * factor * factor
        new.header = self.header.copy()

        return new

    
def divide(a, b):
    """Divides one map by another.

    It also calculates the variance of the ratio
    """
    
    new = Map.empty()

    new.data[0] = ma.divide(a.data[0], b.data[0]) 
    new.data[1] = get_variance_ratio(new.data[0], b.data[0], a.data[0],
                                     b.data[1], a.data[1])

    return new


        
def get_variance_ratio(r, den, num, var_den, var_num):
    """Calculate the variance of a ratio."""
    
    var_r = (r * r) * (var_den / den / den + var_num / num / num)
    return var_r
