#!/usr/bin/env python
# RTPhoS: create fake data stream output
# dates, target and calibration fluxes, errors
# thumbnails for target and comparison star
# Usage:
# RTPhoS_fakedata.py [tsleep] [tfilename] [cfilename]
#    1  1990-01-01|00:00:00  2447892.500058   5.00   33368.6875 647.047546   5.46   O     gauss01.fits 
import os
import sys
import argparse
import time
from datetime import datetime
import numpy as np
from astropy.io import fits
import random

tsleep   = float(sys.argv[1])
tfilename = sys.argv[2]
cfilename = sys.argv[3]
i = 1

try: 
    while True:
        tfile = open(tfilename, "a+")
        cfile = open(cfilename, "a+")
        now = datetime.utcnow()
        ## add two decimal places to seconds
        #micro = float(now.microsecond)
        #micro2 = str(int(round(micro,-4) / 10000.))
        now = now.strftime("%Y-%m-%d|%H:%M:%S")
        #now = now + '.' + micro2
        print now
        
        # Create random array image of 100x100 pixels
        t_image = np.random.random((100,100))
        c_image = np.random.random((100,100))
        t_image = t_image*500
        c_image = c_image*500
        hdu=fits.PrimaryHDU(t_image)
        hdu.writeto("target.fits", clobber='True')
        hdu=fits.PrimaryHDU(c_image)
        hdu.writeto("comp.fits", clobber='True')

        # example fake data output
        id            = str(i)
        obsid         = 'Vidojevica1.4m'
        bandpass      = 'R'
        UTCdatetime   = now
        BJD           = str(round( (2447892.500058+i*tsleep/(24*3600.)), 7))
        terr          = '5.0'
        targetflux    = str(round(10000. + np.random.random()*2000.,2))
        targetfluxerr = str(round(np.sqrt(float(targetflux)),2))
        compflux      = str(round(20000. + np.random.random()*1000.,2))
        compfluxerr   = str(round(np.sqrt(float(compflux)),2))
        seeing        = str(round(1.0 + np.random.random()*0.5,2))    
        flag          = '0'
        dummyname     = 'dummy.fits'
  
        t_table_data = id.ljust(7) + UTCdatetime.ljust(25) + BJD.ljust(16) + terr.ljust(4) +\
                     targetflux.ljust(12) +  targetfluxerr.ljust(8) + seeing.ljust(5) + \
                     flag.ljust(3) + dummyname.ljust(15) +'\n'
 
        c_table_data = id.ljust(7) + UTCdatetime.ljust(25) + BJD.ljust(16) + terr.ljust(4) +\
                     compflux.ljust(12) +  compfluxerr.ljust(8) + seeing.ljust(5) + \
                     flag.ljust(3) + dummyname.ljust(15) +'\n'

        tfile.write(t_table_data)
        cfile.write(c_table_data)
        i = i + 1
        tfile.close()
        cfile.close()
        time.sleep(tsleep)

except (KeyboardInterrupt, SystemExit):
    print "Process aborted by user."
    pass
