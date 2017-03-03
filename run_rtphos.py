#!/usr/bin/python

"""
####################### Real-Time Photometry Server ############################

For information please contact:

Dr. Milan Bogosavljevic            |  Dr. Zach Ioannou
Astronomical Observatory Belgrade  |  Department of Physics
Belgrade, Serbia                   |  Sultan Qaboos University
                                   |  Muscat, Oman
milan@aob.rs                       |  zac@squ.edu.om

################################################################################

This python script monitors a directory for incoming images and performs
aperture and optimal photometry on targets chosen in a DS9 window....
<more description info and basic manual to go here...>


Included functions:
 
dict of floats - Gets a list of string values and makes a python dictionary
make_png       - Converts FITS images to .png images using f2n.py
gauss          - Computes a Gaussian profile for given parameters
fwhm_from_star - Calculates the FWHM of a stellar profile
get_comps_fwhm - Calculates a mean FWHM based on the number of comparison stars
stripdate      - Removes significant figures from a JD type date number
barytime       - Calculates Barycentric Dynamical Julian Date using barycor.f90
getoffsets     - Uses cross-corellation to calculate frame to frame offsets
run_photometry - Performs photometry by calling optimal.f90
positions      - Gets improved offsets by averaging offsets from each star
outputfiles    - Writes output to files
seekfits       - Performs all the necessary calibration and photometry steps
run_rtphos     - Initiates the enviroment variables and starts the sequence
main           - Main program start

Runs with: run_rtphos.py <xpa handle> rtphos.defaults

"""

# Initial imports
import pyregion
import astropy.io.fits as pyfits
from   astropy.time import Time
from   datetime import datetime
import pyds9
import numpy as np
from numpy import inf
import time, math
import ccdcalib
from   scipy.optimize import curve_fit
from   scipy import signal, ndimage
from   subprocess import call, Popen, PIPE
import f2n
import sys, select, os
import subprocess

#import matplotlib.pyplot as plt
#from matplotlib.colors import LogNorm
#from matplotlib import rcParams
#rcParams.update({'figure.autolayout': True})
#import matplotlib.patches as mpatches

# NOTE: where else can we put this?
#sys.path.append("~/pythoncode/f2n/f2n") # The directory that contains f2n.py and f2n_fonts !
#===============================================================================
##############################################################################
def dict_of_floats(list_of_strings, num_items):
    dict_of_floats={}
    xypos={}
    flags={}
    
    for i in range(num_items):
        for j in range(num_items):
            #dummy = [float(x) for x in list_of_strings[j].split()]
            # Keep values as strings since star name is also included
            dummy = [x for x in list_of_strings[j].split()]
            #print ("Dummy", dummy)
            dict_of_floats[j]=dummy[1:3]
            seeing = dummy[3]
            xypos[j] = dummy[4:6]
            flags[j] = dummy[6]
                        
    result = (dict_of_floats, seeing, xypos, flags)
    return (result)


##############################################################################
def make_png(png_image_name, frame_name, data, rbin):
# requires f2n installed

    png_image = f2n.f2nimage(numpyarray=data)  # give the image data to f2n class
    png_image.setzscale("flat","flat")  # works best to my liking
    png_image.rebin(rbin)
    png_image.makepilimage("lin")       # linear image scaling.
    # we can play with marking the star and comparisons in the image, like so
    # png_image.drawcircle(112, 101, r=15) 
    # TBD later, for now just label the frame
    png_image.writetitle(frame_name, colour=(200, 200, 0))
    png_image.tonet(png_image_name)     # write the png.


##############################################################################
# Define model function to be used for PSF fit to the stars:
# in this case its a Gauss with a center Xc, sigma, amplitude A and base level B
def gauss(x, *p):
    A, xc, sigma, B  = p
    return A*np.exp(-(x-xc)**2/(2.*sigma**2)) + B


##############################################################################
def fwhm_from_star(image):
# requires gauss
#   Calculate the azimuthal radial profile.
#   Input image should be a small stamp from the 2D image (already cropped)

    # Calculate the indices from the image
    y, x = np.indices(image.shape)
    center = np.array([(x.max()-x.min())/2.0, (x.max()-x.min())/2.0])
    r = np.hypot(x - center[0], y - center[1])
    
    # select only elements inclosed in a circle of radius r, not a box of side 2r
    mr = np.amax(r)/ 2.0**0.5
    take = np.where(r <= mr)
    cr  = r[take]
    im2 = image[take]
#    hdu=pyfits.PrimaryHDU(im2)
#    hdu.writeto("test.fits")

    # Get sorted radii
    ind = np.argsort(cr.flat)
    r_sorted = cr.flat[ind]
    i_sorted = im2.flat[ind]

    # mirror the curve 
    neg_r = -1 * r_sorted[::-1]
    neg_i = i_sorted[::-1]
    fully = np.concatenate([neg_i,i_sorted])
    fullx = np.concatenate([neg_r,r_sorted])

    # Fit a function 
    # p0 is the initial guess for the fitting coefficients
    # (A, xc, sigma, B in above in gauss function)
    # so make some educated guesses:
    n = len(fullx)
    sigma = 1.
    base = np.mean(i_sorted)
    peak = np.amax(fully) - base

    p0 = [peak, 0., sigma, base]
    coeff, var_matrix = curve_fit(gauss, fullx, fully, p0=p0)

#    For debugging: Uncomment to view the histogram and fit. Don't forget to
#    uncomment the import matplotlib statement at the top.
#    xfit=np.linspace(np.amin(fullx),np.amax(fullx),100)
#    fit = gauss(xfit,*coeff)
#    fwhm = int(2.3548 * coeff[2] * 100) / 100.0    
#    print "Peak position: ", coeff[1]
#    print "Best fit sigma and fwhm: ", coeff[2], fwhm
#    plt.plot(fullx, fully, 'g')
#    plt.plot(xfit,fit,'r')
#    plt.show()

    # convert sigma to fwhm:
    fwhm = int(2.3548 * coeff[2] * 100) / 100.0
    return fwhm


#############################################################################
def get_comps_fwhm(comparisons, xpapoint):
# requires fwhm_from_star, namesplit

    # create a ds9 object linked with an XPA point
    win = pyds9.DS9(xpapoint)
    # load the image from ds9 (not from disk!)
    hdu_link = win.get_pyfits() 
    image = hdu_link[0].data

    nc = len(comparisons)
    tot_fwhm = 0.
    # comparisons is alist of tuples, in each tuple first element has x,y,r
    print "* I will use", nc, "star(s) to get an estimate of the FWHM:"
    for i in range(0,nc):
        xt = comparisons[i][0][0]
        yt = comparisons[i][0][1]
        r  = comparisons[i][0][2]
        x1 = int(xt - r)
        x2 = int(xt + r)
        y1 = int(yt - r)
        y2 = int(yt + r)   
        # caution - I don't know why order of x and y is inverted here
        crop_image = image[y1:y2,x1:x2]
#        hdu=pyfits.PrimaryHDU(crop_image)
#        hdu.writeto("test.fits")
        fwhm = fwhm_from_star(crop_image)
        cname = (comparisons[i][1])
        print i+1, cname, fwhm
        tot_fwhm = tot_fwhm + abs(fwhm)
        
    mean_fwhm = tot_fwhm / nc
    print "* Calculated a mean FWHM of:",  ("%.2f" % abs(mean_fwhm)),"pixels"
    print

    return mean_fwhm

#############################################################################
def stripdate(longjd):

    # reduce a Julian Date type date format to 2 significant figures before the decimal.
    twosig = int(longjd/100.0)*100.0
    shortjd =  jd - twosig

    return (shortjd, twosig)


#############################################################################
def barytime(checklist, dirs):

    # Move to the reduced data directory
    prev_dir = os.path.abspath(os.curdir)
    os.chdir(dirs['reduced'])

    # Deconstruct the Date string from the DATE Keyword
    date  = checklist['DATE'].split('-')
    date = " ".join(date)
    # Deconstruct the Time string from the TIME Keyword
    time  = checklist['TIME'].split(':')
    time = " ".join(time)
    # Deconstruct the RA string from the RA Keyword
    if (checklist['RA'] != "Invalid"):
       checklist['RA']=checklist['RA'].lstrip() # Remove any leading spaces
       if ':' in checklist['RA']: ra  = checklist['RA'].split(':')
       if ' ' in checklist['RA']: ra  = checklist['RA'].split(' ')
    else:
       ra = "NaN"
    ra = " ".join(ra)
    # Deconstruct the Dec string from the DEC Keyword
    if (checklist['DEC'] != "Invalid"):
       checklist['DEC']=checklist['DEC'].lstrip() # Remove any leading spaces
       if ':' in checklist['DEC']: dec  = checklist['DEC'].split(':')
       if ' ' in checklist['DEC']: dec  = checklist['DEC'].split(' ')
    else:
       dec = "NaN"
    dec = " ".join(dec)
    # Get the Exposure time
    exposure = str(checklist['EXPOSURE'])
    if exposure=="Invalid": exposure = "0.0"

    equin  = '2000.0'
    utcorr = '0.0'

#   For Later:
#   Get Equinox form Header
#   Get Time Correction from Header

#   The input to barycor should be the following:
#   equin, rah, ram, ras, decd, decm, decs
#   year, month, day, hrs, mins, secs, utcorr, exposure
    datetimeinfo = []
    inputline1 = equin+" "+" "+ra+  " "+dec
    inputline2 = date +" "+" "+time+" "+" "+utcorr+" "+" "+exposure

    #print inputline1
    #print inputline2

    p = Popen(["barycor"], stdin=PIPE, stdout=PIPE)

    time_BDJD = p.communicate(inputline1+"\n"
                             +inputline2)[0]
                             
    os.chdir(prev_dir)
    return time_BDJD


#############################################################################
def getoffsets(dataref,data2red):
    
    # print "IN:", data2red[100,100]
    
    xshift = 0
    yshift = 0
    
    # Crop the image by 10 pixels on each side
    xsize1 = dataref.shape[0]
    ysize1 = dataref.shape[1]
    
    xstart = 9
    xend   = xsize1-10	
    ystart = 9
    yend   = ysize1-10
    
    croped1 = dataref[xstart:xend,ystart:yend]
    median1 = np.median(croped1)
    
    # Create image 1 mask and blur it using a Gausian filter of
    # FWHM of 1 pixel. Then any pixel with a value less than 100
    # is set to zero.
    
    # WARNING - This is completely arbitrary but usually
    # pixels that correspond to actual stars will have values
    # a lot greater than 100. (Will need to change this to the median sky value)
    
    mask1   = croped1
    mask1[mask1 < 1.5*median1] = 0.0
    blured1 = ndimage.gaussian_filter(croped1, sigma=1)
    mask1   = blured1
    mask1[mask1 < 100.0] = 0.0
    
    # Get the size of the croped masked arrays
    xsize = mask1.shape[0]
    ysize = mask1.shape[1]
    
    # Create collapsed image arrays for reference image dataref
    xvals1 = mask1.sum(axis=0)
    yvals1 = mask1.sum(axis=1)
    
    # Crop the data2 image by 10 pixes on each side
    xsize2 = data2red.shape[0]
    ysize2 = data2red.shape[1]
    xstart = 9
    xend   = xsize2-10
    ystart = 9
    yend   = ysize2-10
    croped2 = data2red[xstart:xend,ystart:yend]
    median2 = np.median(croped2)

    # Create image 2 mask and blur it
    mask2   = croped2
    mask2[mask2 < 1.5*median2] = 0.0
    blured2 = ndimage.gaussian_filter(croped2, sigma=1)
    mask2   = blured2
    mask2[mask2 < 100.0] = 0.0
    # print "OUT:", data2red[100,100]
    
    # Create collapsed image arrays for image 2
    xvals2=mask2.sum(axis=0)
    yvals2=mask2.sum(axis=1)
    
    # Calculate the x and y shift of the image in pixels using cross correlation
    xshift = (np.argmax(signal.correlate(xvals1,xvals2)))-(ysize-1)
    yshift = (np.argmax(signal.correlate(yvals1,yvals2)))-(xsize-1)
    
    return (xshift, yshift)


#############################################################################
def run_photometry(rtdefs, dirs, inputfile, origfile, psf_fwhm, thisoffset, \
                   runpass, alltargets, ntarg, ncomp):

    # Assign values
    filename   = inputfile+" "
    npsf       = str(ncomp)+" "    
    nstar      = str(ntarg)+" "   
    badskyskew = rtdefs['skyskew']+" "
    badskychi  = rtdefs['skyfit']+" "
    clip       = rtdefs['cradius']+" "
    aprad      = rtdefs['aradius']+" "
    iopt       = rtdefs['starnumber']+" "
    searchrad  = rtdefs['sradius']+" "
    adu        = rtdefs['gain']
    verbose    = rtdefs['verbose']
    fwhm       = str(psf_fwhm)+" "

    dummy = os.path.basename(origfile)
    dummy2 = os.path.splitext(dummy)
    flagfile = dummy2[0]+".flg "
    
    if (runpass==2):
       searchrad = "-1.0 "
       aprad     = "-1.0 "

    #verbose="Y"    # Uncomment for debuging

    # Apply the frame offset
    xoff = float(thisoffset[0])
    yoff = float(thisoffset[1])
    name=[]
    xc=[]
    yc=[]
    for i in range(0, ntarg+ncomp):
        xc.append(alltargets[i][0][0]-xoff)
        yc.append(alltargets[i][0][1]-yoff)
        name.append(alltargets[i][1])

    # Print X,Y coords for debuging
    #if (runpass==1):
    #   for i in range(0, ntarg+ncomp):
    #       print "X, Y 1st pass: ", i+1, xc[i], yc[i]
    #else:
    #   for i in range(0, ntarg+ncomp):
    #       print "X, Y 2nd pass: ", i+1, xc[i], yc[i]

    # Construct the optimal command line input
    input_txt=[]
    input_txt.append(filename+flagfile+npsf+nstar+verbose)
    input_txt.append(badskyskew+badskychi+fwhm+clip+aprad+iopt+searchrad+adu)
    for i in range(ntarg, ntarg+ncomp):
        input_txt.append(name[i]+" "+str(xc[i])+" "+str(yc[i]))
    for i in range(0, ntarg+ncomp):
        input_txt.append(name[i]+" "+str(xc[i])+" "+str(yc[i]))
        
    input_str = "\n".join(input_txt) + "\n"
    
    if (runpass==1):
       print "* Initiating photometry calculations..."
       print
       
    #print input_str
    #print
    
    os.chdir(dirs['reduced'])  # Move to the reduced image directory
    p = Popen(["optimal"], stdin=PIPE, stdout=PIPE)
    data_out = p.communicate(input_str)[0]
 
    if verbose=='Y':
        print " ### OPTIMAL OUTPUT START ###"
        print data_out
        print " ### OPTIMAL OUTPUT END ###"

    # Find the lines corresponding to the photometry output.
    # These are the lines after the line with 10 stars.
    results=data_out.split("\n")
    respos = results.index(' **********')
    del results[:respos+1]   # Trim the output list to contain only the results.

    total_records = len(results)-1
    total_stars=total_records/2

    optimal_data = results[0:total_stars]
    aperture_data = results[total_stars:total_records]

    # Debug
    #print results, len(results), respos
    #print ("Optimal data", optimal_data)
    #print ("Aperture data", aperture_data)

    # Gets a dictionary from a list of results (optimal or aperture)
    optimal_res=dict_of_floats(optimal_data, total_stars)
    optimal_stars=optimal_res[0]
    seeing = optimal_res[1]
    xypos = optimal_res[2]
    opflags = optimal_res[3]
    
    aperture_res=dict_of_floats(aperture_data, total_stars)
    if (runpass==1):
       aperture_stars=aperture_res[0]
       seeing = aperture_res[1]
       xypos = aperture_res[2]
       apflags = aperture_res[3]
    else:
       aperture_stars = optimal_stars
       #seeing = " - "
       #xypos = " - "
       apflags = opflags

    photometry_result = (optimal_stars, aperture_stars, seeing, xypos, opflags, apflags)

    os.chdir(dirs['data'])  # Move back to the raw data directory
    return photometry_result



#############################################################################
def positions(optimalist, xyposlist, initx, inity, ntarg):

    xoffs = []
    yoffs = []

    # Go around these loops for every entry and every target star
    for i in range(0,len(optimalist)):
   
        tmp_counts = float(optimalist[i][0])
        tmp_ecounts= float(optimalist[i][1])
        tmp_xpos   = float(xyposlist[i][0])
        tmp_ypos   = float(xyposlist[i][1])
        if (tmp_counts<=0.0 or tmp_ecounts<=0.0):
           tmp_sn = 0.0
        else:
           tmp_sn = tmp_counts/tmp_ecounts  
     
        # Find the x and y offsets from the initial values              
        if (tmp_ecounts >0.0 or tmp_sn>20.0):
           dummy = initx[i]-tmp_xpos
           xoffs.append(dummy)       
           dummy = inity[i]-tmp_ypos
           yoffs.append(dummy)
        
    # Find the median offset for each frame
    framepos = (np.median(xoffs), np.median(yoffs))      

    return framepos                     


#############################################################################
def outputfiles(dirs, alltargets, optimalist, aperatlist, opflaglist, apflaglist, xyposlist, seeing, \
                filterobs, frame_time, frame_timerr, pdatetime, filename, count, runpass):

#   Output text file format example.
#   No.    UTCdatetime       BJD            terr[s]  Flux       Flux_err    seeing flag  filename
#    1  1990-01-01|00:00:00  2447892.500058   5.00   33368.6875 647.047546   5.46   O     gauss01.fits 

    # Move to the reduced files directory
    prev_dir = os.path.abspath(os.curdir)
    os.chdir(dirs['reduced'])

    # Loop will write out files with the following outputformat:
    # seq. number, frame_time, frame_timerr, flux, flux error, seeing
    if (runpass==1):
       for i in range(0,len(alltargets)):
           # First write the optimal photometry data
           with open(alltargets[i][1]+".opt_tmp", "a") as outfile:
                outdata = (count, pdatetime, frame_time, frame_timerr*86400.0, optimalist[i][0], \
                           optimalist[i][1], float(seeing), filterobs, opflaglist[i], filename)
                fmtstring = '%5i %20s %15.6f %6.2f %12s %9s %6.2f %s %s %s \n'
                outfile.write(fmtstring % outdata)
#                outfile.write(str(count)+" "+pdatetime+" "+str(frame_time)+\
#                " "+ str(frame_timerr)+" "+str(optimalist[i][0])+" "+\
#                str(optimalist[i][1])+" "+str(seeing)+" "+filename+" \n")

           # Now write the aperture photometry data
           with open(alltargets[i][1]+".dat", "a") as outfile:
                #outdata = (count, aperatlist[i][0], aperatlist[i][1], xyposlist[i][0], xyposlist[i][1])
                outdata = (count, pdatetime, frame_time, frame_timerr*86400.0, aperatlist[i][0], \
                           aperatlist[i][1], float(seeing), filterobs, apflaglist[i], filename)
                fmtstring = '%5i %20s %15.6f %6.2f %12s %9s %6.2f %s %s %s \n'
                #fmtstring = '%5i %12s %9s %7s %7s \n'
                outfile.write(fmtstring % outdata)
#                outfile.write(str(count)+" "+pdatetime+" "+str(frame_time)+\
#                " "+ str(frame_timerr)+" "+str(aperatlist[i][0])+" "+\
#                str(aperatlist[i][1])+" "+str(seeing)+" "+filename+" \n")

    # Save the optimal data from the 2nd photometry pass.
    if (runpass==2):
       for i in range(0, len(alltargets)):
           with open(alltargets[i][1]+".opt", "a") as outfile:
                #outdata = (count, optimalist[i][0], optimalist[i][1], xyposlist[i][0], xyposlist[i][1] )
                outdata = (count, frame_time, frame_timerr*86400.0, optimalist[i][0], \
                           optimalist[i][1], float(seeing), opflaglist[i])
                fmtstring = '%5i %15.6f %6.2f %12s %9s %6.2f %s \n'
                #fmtstring = '%5i %12s %9s %7s %7s \n'
                outfile.write(fmtstring % outdata)

    ### Files are always appended. This might be a problem when running a
    ### new round of rtphos of the same data. Perhaps we should make rtphos
    ### to first check is *.opt and *.dat files exist in the output directory
    ### delete them and then proceed with the reduction(?)
    ### Should we place the output files in a new directory (e.g. output)?

    # Return to the previous directory
    os.chdir(prev_dir)
    return


#############################################################################
def seekfits(rtdefs, dataref, dirs, tsleep, comparisons, targets, psf_fwhm):
    
    # Combine targets and comparison star initial data (position, names, etc)
    alltargets = targets + comparisons
    ncomp = len(comparisons)  # number of comparison stars
    ntarg = len(targets)      # number of target stars
    allobj = ntarg+ncomp      # number of all objects

    frameslist=[]
    calib_frameslist = []
    frame_times = []
    frame_timerrs = []
    pdatetimes = []
    newoffsets= []    
    initx = []
    inity = []
    dist_dx = []
    dist_dy = []
    all_seeing = []    
    all_opdata = [[] for x in range(allobj)]
    all_apdata = [[] for x in range(allobj)]
    for i in range(0, allobj):
        all_opdata[i] = [[] for x in range(6)]
        all_apdata[i] = [[] for x in range(6)]
    
    # Initialize reduced data lists
    xdata=[]         # X-axis data (time)
    yseeing=[]       # Seeing data
    yrawtarget=[]    # Raw target counts
    yrawtargeterr=[] # Raw target error bars
    yrawcomp=[]      # Raw comparison counts
    yrawcomperr=[]   # Raw comparison error bars
    ydflux=[]        # Differential photometry counts
    ydfluxerr=[]     # Differential photometry error bars

    before = sorted(os.listdir(dirs['data']))
    # a switch to first check if files found on startup
    # need to be reduced or not (if the reduced output exists)
    reducebefore = True

    # counter needed just to print out progress
    count = 0

    while 1:
        print "*LISTENING*", dirs['data'], time.strftime('%X %x %Z')
        print "Hit <Enter> to exit initial photometry loop"
        if sys.stdin in select.select([sys.stdin], [], [], 0)[0]:
           line = raw_input()
           break
        after = sorted(os.listdir(dirs['data']))
        added = [f for f in after if not f in before]

        if reducebefore:
            added = added + before
            reducebefore = False

        if added:   
            for filein in added:                  
                # check if it is a fits file
                filename = dirs['data']+'/'+filein
                   
                # WARNING - hardcoded the '.fits' or '.fit' extensions
                if (filename.endswith('.fits') or filename.endswith('.fit')):
                    # Can load both data and header with this trick
                    data2, hdr = pyfits.getdata(filename, header=True) 
                    print "* Processing image: "+filename
                    junk, sfilename = os.path.split(filename)
                       
                    # Count the number of processed frames
                    count = count + 1 
                       
                    # Convert FITS into PNG to be used later to make
                    # a movie of the timeseries. 
                    # *** Maybe a good idea to do this after calibration but
                    # for now keeping it before. 
                    frame_name = os.path.splitext(os.path.basename(sfilename))[0]
                    png_image_name = dirs['png'] + frame_name + ".png"
                    #print sfilename
                    #print png_image_name
                    if not os.path.isfile(png_image_name):
                       print "PNG file of this frame does not exist (yet)"
                       # hardcoded rebin factor 2 here, TBD later
                       make_png(png_image_name, sfilename, dataref, 2)
                     
                    # First check that all the required header keywords are in
                    # the FITS file and then get the time stamp for this frame.
                    # If the RA and DEC of the image are in the headers then
                    # time will be in Barycentric Dynamical Julian Date. If not 
                    # then time will be in plain simple Julian Date.
                    checklist = ccdcalib.makechecklist(hdr)
                       
                    exposure = float(checklist['EXPOSURE'])   # in seconds
                    midexp = (exposure / 2.0) / 86400.        # in days
                    mdatetime = checklist['DATE']+" "+checklist['TIME']
                    pdatetime = checklist['DATE']+"|"+checklist['TIME']
                    # this is needed to get plot-able UTC time
                    sdatetime = datetime.strptime(mdatetime,  "%Y-%m-%d %H:%M:%S")
                    mdatetime = Time(mdatetime, format='iso', scale='utc')

                    if checklist['RA']=="Invalid" or checklist['DEC']=="Invalid":  
                        frame_time  = mdatetime.jd            # in JD
                    else:
                        time_BDJD = barytime(checklist, dirs) # in BDJD
                        time_BDJD = time_BDJD.split()
                        frame_time = float(time_BDJD[0])
                       
                    # Make frame time the middle of the exposure
                    frame_time = frame_time + midexp 
                    frame_times.append(frame_time)
                    frame_timerrs.append(midexp)
                    pdatetimes.append(pdatetime)

                    # Strip JD and reduced it to 2 significant figures.
                    stripjd = int(float(frame_time)/100.0)*100.0
                    twosig_time =  float(frame_time) - stripjd
                       
                    # Now initiate the calibration, offsets and photometry.
                    # First map all the saturated and non linear pixels of the
                    # current image. Then combine with the bad pixel mask to 
                    # make a final mask of unwanted pixels and flag them accordingly. 
                    # The flag file is a FITS file with an .flg extension. 
                    # See ccdcalib.py for more info.
                    ccdcalib.pixflag(rtdefs, dirs, filename, data2, hdr)
                                                         
                    # ccdcalib will either calibrate the image and place the
                    # calibrated image file in the '/reduced/' directory or
                    # if the image did not require calibration just copy the image
                    # to the '/reduced/' directory. In either case the image will
                    # have a prefix to indicate that ccdcalib has seen it.
                    calib_data = ccdcalib.calib(rtdefs, dirs, filename, data2, hdr)
                    (data2, hdr, calib_fname) = calib_data
                    
                    # Find offsets from reference frame (displayed on DS9)
                    thisoffset = (0,0)
                    if (count>1):
                       #print "Here1:", data2[100,100]
                       thisoffset = getoffsets(dataref,data2)
                       #print "Here2:", data2[100,100]   
                       
                    # For some reason the data2 array gets altered after the call
                    # to getoffsets. Need to find out why later. 
                    # For now just read it in again.    
                    data2, hdr = pyfits.getdata(filename, header=True) 
                    
                    print "* Frame Offsets: (x,y) ", thisoffset

                    # Call optimal and do the photometry.
                    frame_photometry = run_photometry(rtdefs, dirs, calib_fname, filename, \
                                                      psf_fwhm, thisoffset, 1, alltargets, ntarg, ncomp)
  
                    # Deconstruct photometry results from optimal.f90
                    (optimaldict, aperatdict, seeing, xypos, opflags, apflags) = frame_photometry
                    optimalist = optimaldict.values()
                    aperatlist = aperatdict.values()
                    xyposlist  = xypos.values()
                    opflaglist = opflags.values()
                    apflaglist = apflags.values()
                       
                    # Put the optimal photometry results into one big list of lists
                    # Each item in the large list corresponds to each target
                    # The lists in each target list correspond to:
                    # Counts, Error, X, Y, S/N, Flag
                    for j in range(0,len(optimalist)):
                        tmp_cts  = float(optimalist[j][0])
                        tmp_ects = float(optimalist[j][1])
                        if (tmp_ects<=0.0 or tmp_cts<=0.0):
                           tmp_sn = 0.0
                        else:
                           tmp_sn = tmp_cts/tmp_ects
                              
                        all_opdata[j][0].append(tmp_cts)
                        all_opdata[j][1].append(tmp_ects)
                        all_opdata[j][2].append(float(xyposlist[j][0]))
                        all_opdata[j][3].append(float(xyposlist[j][1]))
                        all_opdata[j][4].append(tmp_sn)
                        all_opdata[j][5].append(opflags[j])
                      
                    all_seeing.append(seeing)          

                    # Get better frame offsets based on optimal photometry positions
                    # and calculate the relative distances
                    ###### (FOR LATER): If optimal results not available
                    ######              keep the offsets from image cross-correlation routine.
                    if (count==1):
                       for i in range(0,len(optimalist)):
                           initx.append(float(xyposlist[i][0]))
                           inity.append(float(xyposlist[i][1]))
 
                    framepos = positions(optimalist, xyposlist, initx, inity, ntarg)
                    frameoffs = (framepos[0], framepos[1])                           
                    newoffsets.append(frameoffs)
                       
                    # Put the original and reduced filenames in a list
                    junk, sfilename = os.path.split(filename)
                    frameslist.append(filename)
                    calib_frameslist.append(calib_fname)
                       
                    # Get the Filter used for this image
                    filterobs = checklist['FILTER']   
                    if filterobs=="Invalid":
                       filterobs="INV"
                       
                    # Write the result to the output files
                    outputfiles(dirs, alltargets, optimalist, aperatlist, opflaglist, \
                    apflaglist, xyposlist, seeing, filterobs, frame_time, midexp, \
                    pdatetime, sfilename, count, 1)
                      
                    # Text output
                    print "================================================="
                    print "Filename: ", filename
                    print "Frame time           : ", frame_time, midexp
                    print "Refined frame offsets: ", "%.4f %.4f" % frameoffs
                    print "Optimal Photometry Results:"
                    for i in range(0,len(optimalist)):
                        print alltargets[i][1], optimalist[i][0], optimalist[i][1], \
                              seeing, xyposlist[i][0], xyposlist[i][1], opflaglist[i]
                    print "-------------------------------------------------"
                    print "Aperature Photometry Results:"
                    for i in range(0,len(aperatlist)):
                        print alltargets[i][1], aperatlist[i][0], aperatlist[i][1], \
                              seeing, xyposlist[i][0], xyposlist[i][1], apflaglist[i]
                    print
                    
                    # Live plotting
                    if rtdefs['liveplot']:
                       # Crop target and comparison star images and output to file
                       # for transmition.Images are overwritten with the same filename.
                       # image names are hardcoded to target.fits and comp.fits
                       # First remove inf values from the image array and set them to 0.
                       prev_dir = os.path.abspath(os.curdir)
                       os.chdir(dirs['reduced'])
                       data2[data2 == -inf] = 0.0             
                       # Crop a 100px square around the target and write to file
                       targetx = float(xyposlist[0][0])
                       targety = float(xyposlist[0][1])
                       target_crop = data2[int(targety-50):int(targety+50),int(targetx-50):int(targetx+50)]
                       hdu=pyfits.PrimaryHDU(target_crop)
                       hdu.writeto("target.fits", clobber='True')
                       # Crop a 100px square around the first comparison star and write to file
                       compx = float(xyposlist[1][0])
                       compy = float(xyposlist[1][1])
                       comp_crop = data2[int(compy-50):int(compy+50),int(compx-50):int(compx+50)]
                       hdu=pyfits.PrimaryHDU(comp_crop)
                       hdu.writeto("comp.fits", clobber='True')
                      
                       # Merge the target data with the data of the first comparison star
                       # and put them in a file to be read by the live plotting module.
                       # Since the file is of no other use the name is hardcoded here as
                       # liveplotdata.txt
                       port=5556      #HARDCODED HERE 
                       with open("liveplotdata.txt", "a") as outfile:
                            # output data written is:
                            # count, server port, filter, BJD, target name, target flux, target eflux, \
                            # comp name, comp flux, comp eflux, seeing
                            outdata = (count, port, filterobs, frame_time, alltargets[0][1], aperatlist[0][0], aperatlist[0][1],\
                                       alltargets[1][1], aperatlist[1][0], aperatlist[1][1], float(seeing))
                            fmtstring = '%5i %4i %3s %15.6f %8s %12s %9s %8s %12s %9s %6.2f \n'
                            outfile.write(fmtstring % outdata)
                    
                       # Constract the shell command to run the live plotting module.
                       # Need to find a better way but for now the only way for the code
                       # not to hang when launching a process is to create and run a
                       # shell script. Extremely ugly but it will do for now!
                       if count==1:
                          command1 = 'cp '+dirs['current']+'/liveplot.py '+'. \n'
                          command2 = 'xterm -hold -sb -sl 2000 -e tcsh -c "python liveplot.py '+str(port)+' '+filterobs+\
                                    ' -logfile "liveplotdata.txt" " &'
                          with open("rtphos_liveplot.csh", "w") as outfile:
                               outfile.write(command1)
                               outfile.write(command2)
                          
                          # Read that os.system is being derecated and replaced with subprocess...
                          # The same goes for all os calls. We should eventually replace all os calls.
                          subprocess.call("source ./rtphos_liveplot.csh", shell=True)
                   
                       os.chdir(prev_dir)    
                    
                    
                                        
                    # Fill the data lists
                    #xdata.append(twosig_time)
                    #yseeing.append(seeing)
                    #yrawtarget.append(float(optimalist[0][0]))
                    #yrawtargeterr.append(float(optimalist[0][1]))
                    #yrawcomp.append(float(optimalist[1][0]))
                    #yrawcomperr.append(float(optimalist[1][1]))
                    # Do the differential photometry calculations
                    #tcounts    = float(optimalist[0][0])
                    #terror     = float(optimalist[1][1])
                    #ccounts    = float(optimalist[1][0])
                    #cerror     = float(optimalist[1][1])
                    #ydfluxs    = (tcounts/ccounts)
                    #ydfluxerrs = ydfluxs*math.sqrt( ((terror/tcounts)**2.0) + \
                    #                              ((cerror/ccounts)**2.0) )
                    #ydflux.append(ydfluxs)                       
                    #ydfluxerr.append(ydfluxerrs)
                    
                    #time.sleep(10)
                    
                    
        before = after
        time.sleep(tsleep)   # Wait for tsleep seconds before repeating
           
    print
    print
    print "* Performing second pass for optimal photometry..." 
    print
    ###### WORK IN PROGRESS! THIS IS THE PROPER WAY OF RUNNING THE OPTIMAL
    ###### PHOTOMETRY CODE (READ THE PAPERS!)   
       
    # Setup lists for optimal's second pass
    avdx = []
    avdy = []
    comp_initx = []
    comp_inity = []
    targ_initx = []
    targ_inity = []
    targ_fixed_x = []
    targ_fixed_y = []
    deltaxs = [[] for x in range(ntarg)]
    deltays = [[] for x in range(ntarg)]
    avdx = [[] for x in range(ntarg)]
    avdy = [[] for x in range(ntarg)]
    for i in range(0, ntarg):
        deltaxs[i] = [[] for x in range(ncomp)]
        deltays[i] = [[] for x in range(ncomp)]

    # Calculate the average difference in position between each target and 
    # each comparison star. Since comparison stars are normally bright objects
    # taking the average distances from all frames provides a more accurate
    # position for the target stars. Especially useful for faint targets!
    for i in range(0, ntarg):
        targ_initx.append(all_opdata[i][2][0])
        targ_inity.append(all_opdata[i][3][0])
        dumtects = all_opdata[i][1]
        dumtargx = all_opdata[i][2]
        dumtargy = all_opdata[i][3]
        for j in range(ntarg,allobj):
            dumcects = all_opdata[j][1]
            dumcompx = all_opdata[j][2]
            dumcompy = all_opdata[j][3]
            dumsn    = all_opdata[j][4]
            for k in range(0, len(dumtargx)):
                # Only get the distance if star is visible and has a good S/N
                if (dumtects[k]>0.0 and dumcects[k]>0.0 and dumsn[k]>50.0):
                   if (i==0 and k==0):
                      comp_initx.append(dumcompx[k])
                      comp_inity.append(dumcompy[k])
                   dummy = dumtargx[k]-dumcompx[k]
                   deltaxs[i][j-ntarg].append(dummy)
                   dummy = dumtargy[k]-dumcompy[k]
                   deltays[i][j-ntarg].append(dummy)
                      
    # Put the calculated distances dx, dy in a list               
    for i in range(0,ntarg):
        for j in range(0, ncomp):
            dummy = deltaxs[i][j]
            avdx[i].append(np.mean(dummy))
            dummy = deltays[i][j]
            avdy[i].append(np.mean(dummy))
              
    # Fix the coordinates of the target stars to the more accurate ones
    # found above.
    for i in range(0,ntarg):
        dum_fixed_x =  np.asarray(comp_initx)+np.asarray(avdx[i])  
        dum_fixed_y =  np.asarray(comp_inity)+np.asarray(avdy[i])
        targ_fixed_x.append(np.mean(dum_fixed_x))
        targ_fixed_y.append(np.mean(dum_fixed_y))
        alltargets[i][0][0] = targ_fixed_x[i]
        alltargets[i][0][1] = targ_fixed_y[i]
           
    # Fix the coordinates of the comparison stars to those found by 
    # optimal in the first run. The more accurate offsets obtained earlier
    # assure a more accurate fix than the cross-correlation shift method.    
    for i in range(ntarg, allobj):
        alltargets[i][0][0] = all_opdata[i][2][0]
        alltargets[i][0][1] = all_opdata[i][3][0]

    # Print the new fixed positions
    print "* Starting frame is: " 
    print dirs['reduced']+calib_frameslist[0]
    print
    print "* Accurate star positions:"
    for i in range(0, allobj):
        print alltargets[i][1], alltargets[i][0][0], alltargets[i][0][1]
           
    print

    # Call the optimal photometry routine a second time to refine the photometry.
    for i in range(0, len(frameslist)):
        print "Processed file ",i+1," out of ", len(frameslist), \
              "- Frame Offset: ", "%.4f %.4f" % newoffsets[i]
           
        frame_photometry = run_photometry(rtdefs, dirs, calib_frameslist[i], frameslist[i], \
                               psf_fwhm, newoffsets[i], 2, alltargets, ntarg, ncomp)

        # Deconstruct photometry results from optimal.f90
        (optimaldict_2, aperatdict_2, seeing_2, xypos_2, opflags_2, apflags_2) = frame_photometry
        optimalist_2 = optimaldict_2.values()
        aperatlist_2 = aperatdict_2.values()
        xyposlist_2  = xypos_2.values()
        opflaglist_2 = opflags_2.values()
        apflaglist_2 = apflags_2.values()

        # Write the result to the output files
        count=i+1
        outputfiles(dirs, alltargets, optimalist_2, aperatlist_2, opflaglist_2, \
                    apflaglist_2, xyposlist_2, seeing_2, filterobs, frame_times[i], \
                    frame_timerrs[i], pdatetimes[i], "  ", count, 2)

    print
    print "==== RTPhoS END ==== " + time.strftime('%X %x %Z') 
    print
    print         

    return
    
#############################################################################
def run_rtphos(rtphosdir, xpapoint, pathdefs):
# requires get_comps_fwhm, seekfits
  
    print
    print
    print
    print "==== RTPhoS START ==== " + time.strftime('%X %x %Z')
    print

    ######## SETTING UP DIRECTORIES AND LINKS
    # create a ds9 object linked with an XPA point
    win = pyds9.DS9(xpapoint)
    # image name which is displayed in ds9 
    ref_filename = win.get("file")

    # Set up input and output directories
    path, filename = os.path.split(ref_filename) 

    # Read in default values from file
    with open (pathdefs, "r") as defsfile:
        defs = defsfile.readlines()
        print "* Read defaults from: "+ pathdefs
        print
        for l in range(0,len(defs)):
            print l, defs[l][:-1] # without the \n
        print "========================================================= "
        print

    ## If you would want to grep the defaults file and find match
    ## this would be one way to do it
    #Set data directory from defaults file
    #findthis = "# Raw data frames"
    #matching = [s for s in defs if findthis in s]
    #data_dir = matching[0].split()[0]

    # Will use positional assignement for defaults
    # These are all string types:    
    data_dir   = defs[1].split()[0]
    bias_dir   = defs[2].split()[0]
    dark_dir   = defs[3].split()[0]
    flat_dir   = defs[4].split()[0] 
    biaswc     = defs[5].split()[0]
    darkwc     = defs[6].split()[0]
    flatwc     = defs[7].split()[0]
    mbias      = defs[8].split()[0]
    mdark      = defs[9].split()[0]
    mflat      = defs[10].split()[0]
    psfs       = defs[11].split()[0]
    stars      = defs[12].split()[0]
    cprefix    = defs[13].split()[0]
    sradius    = defs[14].split()[0]
    aradius    = defs[15].split()[0]
    cradius    = defs[16].split()[0]
    starnumber = defs[17].split()[0]
    skyskew    = defs[18].split()[0]
    skyfit     = defs[19].split()[0]
    gain       = defs[20].split()[0]
    linlevel   = defs[21].split()[0]
    satslevel  = defs[22].split()[0]
    # These are numbers:
    verbose    =   int(defs[23].split()[0])
    tsleep     =   int(defs[24].split()[0])
    liveplot   =   int(defs[25].split()[0])

    # Current root working dir
    current_dir = os.path.abspath(os.path.join(data_dir, os.pardir))

    reduced_dir  = current_dir+"/reduced/"
    png_dir      = current_dir+"/png/"

    if not os.path.exists(reduced_dir): 
        os.makedirs(reduced_dir)
        print "* Created output reduced dir: " + reduced_dir
    if not os.path.exists(png_dir): 
        os.makedirs(png_dir)
        print "* Created output png dir: " + png_dir
    
    # Create the symbolic links required for running barycor.f90
    os.chdir(reduced_dir)
    print "* Creating symbolic links to JPL Ephemeris and leap data files..."
    call(['ln', '-s', rtphosdir+'/Timing/jpleph.dat', 'JPLEPH'])
    call(['ln', '-s', rtphosdir+'/Timing/leap.dat', 'leapdat'])
    os.chdir(data_dir) # Move back to the data directory
    print
    print 'ln', '-s', rtphosdir+'/Timing/jpleph.dat', 'JPLEPH'
    print

    # Convert verbose switch value to a string
    if verbose==1:
       verbose='Y'
    else:
       verbose='N'
    # Convert liveplot switch value to boolean
    if liveplot==1:
       liveplot=True
    else:
       liveplot=False

    liveplot=True

    # Make a dictionary with all the required directories.
    dirs = {'current':current_dir, 'bias':bias_dir, 'dark':dark_dir, \
            'flat':flat_dir, 'data':data_dir, 'reduced':reduced_dir, 'png':png_dir}

    # Make a dictionary with all the required parameters.
    rtdefs = {'biaswc':biaswc, 'darkwc':darkwc, 'flatwc':flatwc, 'mbias':mbias,\
               'mdark':mdark,   'mflat':mflat, 'psfs':psfs, 'stars':stars, 'cprefix':cprefix, \
              'sradius':sradius, 'aradius':aradius, 'cradius':cradius, 'starnumber':starnumber, \
              'skyskew':skyskew, 'skyfit':skyfit, 'gain':gain, 'linlevel':linlevel, \
              'satslevel':satslevel, 'verbose':verbose, 'tsleep':tsleep, 'liveplot':liveplot}

    ######## DS9 SOURCE IDENTIFICATION AND FWHM ESTIMATE ###########
    # Ceontrid regions in DS9
    win.set("regions select all")

    #### DISABLED DS9 CENTROIDING. IT IS TOO UNRELIABLE. BETTER TO RELY ON USER
    #### IN ANY CASE THE TARGETS WILL BE CENTROIDED BY OPTIMAL.F90 IN A MUCH
    #### MORE ROBUST WAY!
    # I do it twice with different radii on purpose
    #win.set("regions centroid radius 20")
    #win.set("regions centroid iteration 20")
    # not happy at all with DS9 centering so repeating it 20 times
    #for x in range(0, 19):
    #    win.set("regions centroid")
    #win.set("regions centroid radius 2")
    #win.set("regions centroid iteration 5")
    #win.set("regions centroid")
    # put back the default
    #win.set("regions centroid radius 20")
    #win.set("regions centroid iteration 20")

    # save regions file for later
    win.set("regions format ds9")
    win.set("regions save "+ ref_filename +".reg")
    
    # Get source (target, comparison) lists from regions selected
    sourcelist = win.get("regions selected")
    sources    = pyregion.parse(sourcelist)
    n = len(sources)
    print "* Reading target and comparison star lists..."
    print "* In DS9 I see "+ str(n) + " sources labeled:"

    for l in range(0,n):
       if (sources[l].comment is not None):
            sources[l].comment = sources[l].comment[sources[l].comment.find("{")+1:\
            sources[l].comment.find("}")]
            print l+1, sources[l].comment, sources[l].coord_list
       else:
            print "WARNING: You have some unlabeled sources!"
            print sources[l].coord_list
            print "====  Exiting. ==== "
            print
            return

    # find out which are the comparison stars (they must have "C-" in name)
    comparisons  = [(s.coord_list,s.comment) for s in sources if "C-" in s.comment]
    nc = len(comparisons)
    targets =  [(s.coord_list,s.comment) for s in sources if not("C-" in s.comment)]
    nt = len(targets)

    if nc == 0:
        print "WARNING: Must have one comparison star labeled as 'C-<name>' "
        raise Exception("Must have one comparison star labeled as 'C-<name>' ")
    
    if nt == 0:
        print "WARNING: No target selected!"
        raise Exception("No target selected!")
    
    print
    print "* Target stars:"
    print targets
    print
    print "* Comparison stars:"
    print comparisons
    print

    # get FWHM of stellar PSF using comparison stars
    psf_fwhm = get_comps_fwhm(comparisons, xpapoint)

    ##### START PIPELINE #############################################
    # This is where the pipeline looks at the data for the first time!
    print "* Calibrating the file displayed in DS9:"
    print ref_filename
    print
    dataref, hdr = pyfits.getdata(ref_filename, header=True)     

    # Check first DS9 image calibration and calibrate if required.   
    result = ccdcalib.calib(rtdefs, dirs, ref_filename, dataref, hdr)
    (dataref, hdr, calib_fname) = result

    print "####################################################################"
    print "Starting pipeline routine..."
    seekfits(rtdefs, dataref, dirs, tsleep, comparisons, targets, psf_fwhm)
    
    # Clean-up auxillary files
    if rtdefs['liveplot']:
       # Move to the reduced files directory
       prev_dir = os.path.abspath(os.curdir)
       os.chdir(dirs['reduced'])
       subprocess.call('rm -rf *_tmp *_diff liveplot.py comp.fits target.fits \
                       liveplotdata.txt rtphos_liveplot.csh', shell=True) 
       os.chdir(prev_dir)   
    

if  __name__ == "__main__":

    import sys
    rtphosdir, dummy = os.path.split(sys.argv[0])
    xpapoint       = sys.argv[1]
    pathdefs       = sys.argv[2]
    run_rtphos(rtphosdir, xpapoint, pathdefs)
    


    
    
    
