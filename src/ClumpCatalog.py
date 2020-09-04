#!/usr/bin/env python3
##
## ClumpCatalog.py
## O. Morata 2020
##
""" ClumpClass.py

    Definition of the class Clump 

"""

import numpy as np
from astropy.io import fits

import ClumpClass as cl


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
    def from_file(cls, fname):
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

        
        
    def id_isincatalog(self, i) :
        """Check if clump with id equal to i is in the catalog."""
        
        for clump in self.clumps :
            if i == clump.id :
                return clump
            
        return False

    
        
    def fillfrom_rec(self, rec):
        """Fill the a ClumpCatalog object from a FITS_rec object."""

        ncl = np.shape(rec)[0]
        newclumps = 0
        
        for icl in range(ncl) :
            if 'PIDENT' in rec.dtype.names:
                read_cl = rec['PIDENT'][icl] - 1
            elif 'id' in rec.dtype.names:
                read_cl = rec['id'][icl] - 1

            if not self.ncl :
                self.clumps.append(cl.Clump(icl, record=rec[icl]))
                newclumps += 1

            else :
                clump = self.id_isincatalog(read_cl)
                if clump :
                    clump.fill_clump(rec[icl])
                else :
                    print(" >> WARNING: found new records")
                    self.clumps.append(cl.Clump(self.ncl, record=rec[icl]))
                    newclumps += 1

        self.ncl += newclumps


        
    def fillfrom_calcphys(self, idxs=None, **kwargs):
        """Fill a ClumpCatalog Object from calcphys.

        fluxes, temps, mass, and coldens are suposed to be lists of type
        Maps() objects
        """
        
        n_clumps = np.int(np.nanmax(idxs))
        
        newclumps = 0
        
        for icl in range(1, n_clumps+1):
            id = icl - 1
            
            if not self.ncl :
                self.clumps.append(
                    cl.Clump.from_calcphys(id, idxs=idxs, **kwargs))
                newclumps += 1

            else :
                clump = self.id_isincatalog(icl)
                if clump :
                    # WARNING: this line has to be tested yet
                    clump.calc_phys(id, idxs=idxs, **kwargs)

        self.ncl += newclumps
                    

            
    def getview(self, ctype):
        """Returns a selected view of the catalog."""
        
        if ctype == 'all' :
            dt = cl.Clump.ddtype
        elif ctype == 'phys' :
            dt = cl.Clump.dtype_phys

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


        
    def print_catalog(self, ctype='phys', fields=None, filename='',
                      header=True) :
        """Prints the catalog

        ctype:
           - 'all' the fields
           - 'phys' calc_phys fields
           - 'custom' custom set of fields (to finish)
        """
        ##print units
        ## print clumpp

        if header :
            head_on = False
        else:
            head_on = True
        
        if filename:
            f = open(filename, 'w')

        
        for clump in self.clumps:
            str_header, str_out = clump.print_clump(ctype, fields)

            if filename :
                if not head_on:
                    f.write(str_header+'\n')
                    head_on = True
                f.write(str_out+'\n')
                
            else:
                if not head_on:
                    print(str_header)
                    head_on = True
                print(str_out)

