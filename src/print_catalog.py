#!/usr/bin/env python3

import ClumpCatalog as cl

clfile = 'cl1.fits'

catlg = cl.ClumpCatalog.from_file([clfile])

catlg.print_catalog(ctype='findclumps', fields=['PIDENT', 'Shape'],
                    filename='fw')

