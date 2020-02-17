#!/usr/bin/env python3
##
## MapClass.py
## O. Morata 2020
##
""" MapClass

    Definition of the class Map 

"""

import numpy as np
from astropy.io import fits
from datetime import datetime



class Map (object):

    def __init__(self, name='', filename='', data='', header=''):
        """Initialization of the Map object."""

        self.name = name
        self.fname = filename

        if data == '' :
            self.data, self.header = Map.read_fitsfile(self)
        else :
            self.data = data
        
        
    def read_fitsfile(self):
        """Read input data and header from a FITS file."""
    
        print(" >> reading",self.name,"data...")
    
        with fits.open(self.fname) as hdu_data:
            data_info = [hdu_data[0].data, hdu_data[1].data]
            header_info = [hdu_data[0].header, hdu_data[1].header]

            return data_info, header_info


        
    def save_fitsfile(self, oldheader='', append=False, overwrite=False,
                      hdr_type=''):

        if oldheader :
            self.header = Map.modify_header(oldheader, hdr_type)
        else :
            return 3

    
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
    
