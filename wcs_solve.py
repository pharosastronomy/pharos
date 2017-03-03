#!/usr/bin/python
# M. Bogosavljevic, AOB, June 2015
#
# wcs_solve(path, filename)
# path - paths that contains files wcs.fits and axy.fits
# filename - the fits file itself

import os

def wcs_solve(path, filename):
 
       # NOTE: Change this later to some other directory
       indexexists = (os.path.isfile(path+'fieldindex.fits')) 
       wcsexists = (os.path.isfile(path+'wcs.fits') and os.path.isfile(path+'axy.fits'))
 
       if indexexists: 
           mycmd = 'solve-field --config fieldindex.cfg --downsample 2 --no-plots --no-tweak --overwrite ' + filename
           os.system(mycmd)
       else:
           if wcsexists:   
                  mycmd = 'wcs-xy2rd -w wcs.fits -i axy.fits -o rd.fits'
                  os.system(mycmd)

                  # NOTE: one caveat here is -P 4 is hardcoded, which is suitable for 
                  # images roughly 24 arcmin across
                  # see http://astrometry.net/doc/build-index.html#build-index
                  mycmd = 'build-astrometry-index -i rd.fits -o fieldindex.fits -P 4 -E'
                  os.system(mycmd)

                  mycmd = 'echo "add_path ./\nindex fieldindex.fits\ninparallel\n" > fieldindex.cfg'
                  os.system(mycmd)

                  mycmd = 'solve-field --config fieldindex.cfg --downsample 2 --no-plots --no-tweak --overwrite ' + filename
                  os.system(mycmd)
           else: 
               print (" Solve one frame for WCS first on nova.astrometry.net! ")
               print (" Proceeding without WCS solution ")
               return


if __name__ == "__main__":

   import sys
   path     = sys.argv[1]
   filename = sys.argv[2]
   
   wcs_solve(path, filename)
