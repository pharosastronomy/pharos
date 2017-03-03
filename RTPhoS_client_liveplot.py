#!/usr/bin/env python
# RTPhoS live plotting client module which either
# a) receives data from server module live, in which case nothing is saved 
# or
# b) monitors and plots the updating of logged data created by RTPhoS_client_logger.py

# Usage:
# RTPhoS_client_liveplot.py {-port [port] | -logfile [filename]} -filter bandpass
# plus optional params and switches, do -h or --help to get options

# port should be a string like tcp://localhost:port or other address
# Examples:
# RTPhoS_client_liveplot.py -port tcp://localhost:5556 
# RTPhoS_client_liveplot.py -logfile mylogfile.dat 
# if both port and logfile are specified, the code will complain
# and continue with using the logfile, NOT port

import os
import sys
import argparse
import time, math
from datetime import datetime
import numpy as np
from astropy.io import fits
import zmq
import json
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
rcParams.update({'backend': 'TkAgg'})

####################
# minutes_before_now
def minutes_before_now(now, sometime):
# check how many minutes sincce time of data taken (assuming UTC)
#  expects sometime as 1990-01-01|00:00:00.00
    sometime2 = datetime.strptime(sometime, "%Y-%m-%d|%H:%M:%S")
    elapsedTime = sometime2 - now
    minutes = elapsedTime.total_seconds() / 60.
    return minutes

##################
# plot_initialize
def plot_initialize():
    # Initialize Plotting Environment
    plt.ion()
    fig = plt.figure(figsize=(12,12))
    ax1  = fig.add_subplot(421)
    ax2  = fig.add_subplot(422)
    ax3  = fig.add_subplot(412)
    ax4  = fig.add_subplot(413, sharex=ax3)
    ax5  = fig.add_subplot(414, sharex=ax3)

    # make plot settings once
    ax1.text(1, 5, 'Target', color='yellow', fontsize=10)
    ax1.set_xlabel("Physical X: ")
    ax1.set_ylabel("Physical Y: ") 
    ax1.set_xticklabels([])
    ax1.set_yticklabels([])
    
    ax2.text(1, 5, 'Comp', color='yellow', fontsize=10)
    ax2.set_xlabel("Physical X: ") 
    ax2.set_ylabel("Physical Y: ")
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])
    
    # Target & Comparison raw counts
    ax3.grid(True)
    ax3.margins(0.1,0.1)
    ax3.set_ylabel('Raw Counts')

    ax4.grid(True)
    ax4.margins(0.1,0.1)
    ax4.set_ylabel('Relative Flux')
    ax5.grid(True)
    ax5.margins(0.1,0.1)
    ax5.set_ylabel('Pixels')
    ax5.set_ylim([0,3])

    return fig, ax1, ax2, ax3, ax4, ax5
    
#############################
#
def diff_photom(targetflux, targetfluxerr, compflux, compfluxerr):
    # Do the differential photometry
    ydfluxval    = (targetflux/compflux)
    ydfluxerrval = ydfluxval*math.sqrt( ((targetfluxerr/targetflux)**2.0) + \
                                   ((compfluxerr/compflux)**2.0) )

    return ydfluxval, ydfluxerrval

#############################
#
def init_datapoins(tempfile, now, bjdswitch):
    # Initialize data lists
    xdata=[]         # X-axis data (time, either past minutes or BJD)
    yrawtarget = []  # Raw target counts
    yrawtargeterr=[] # Raw target error bars
    yrawcomp=[]      # Raw comparison counts
    yrawcomperr=[]   # Raw comparison error bars
    ydflux=[]        # Differential photometry counts
    ydfluxerr=[]     # Differential photometry error bars
    yseeing=[]       # Seeing data
    UTCtimes = []
    nline = 0
    with open(tempfile,'r') as fp:
        for line in fp:
            useline = line.strip()  
            columns = useline.split()
            # record current position in logfile
            
            # skip comment line
            if columns[0] == '#RTPhoS:':
                print "skipped:", columns
            else:
                nline = nline + 1
                obsid         = columns[0]
                serverport    = columns[1]
                UTCdatetime   = columns[2]
                bandpass      = columns[3]
                BJD           = np.float64(columns[4]) 
                targetflux    = np.float32(columns[5])
                targetfluxerr = np.float32(columns[6])
                compflux      = np.float32(columns[7])
                compfluxerr   = np.float32(columns[8])
                seeing        = np.float16(columns[9])

                # TBD:
                # use just data that matches -filter argument
                # something like matches = [x for x in a if x=='a']
                
                UTCtimes.append(UTCdatetime)

                if bjdswitch:
                    xdata.append(BJD)
                else:
                    # UTCtimes is a list, stores all actual times for
                    # recalculating minutes_before_now later
                    newx = minutes_before_now(now, UTCdatetime)
                    xdata.append(round(newx,3))

                yrawtarget.append(targetflux)
                yrawtargeterr.append(targetfluxerr)
                yrawcomp.append(compflux)
                yrawcomperr.append(compfluxerr)
                ydfluxval, ydfluxerrval = diff_photom(targetflux, targetfluxerr, \
                                                      compflux, compfluxerr)
                ydflux.append(ydfluxval)
                ydfluxerr.append(ydfluxerrval)
                yseeing.append(seeing)

    print  "RTPhoS live plot: initial read found %s rows of data in log" % str(nline)
    return xdata, yrawtarget, yrawtargeterr, yrawcomp, yrawcomperr, \
               ydflux, ydfluxerr, yseeing, UTCtimes

##############################
# plot_datalog
def init_plot_data(fig, ax1, ax2, ax3, ax4, ax5, \
                 xdata, yrawtarget, yrawtargeterr, yrawcomp, \
                      yrawcomperr, ydflux, ydfluxerr, yseeing, bjdswitch):


    xdata2 = []
    # convert to hours since min(BJD)
    lowestjd = np.floor(np.min(xdata))
    if bjdswitch:
        for k in range(0,len(xdata)):
            hours = (xdata[k]-lowestjd)*24.0
            xdata2.append(hours)
    else:
        xdata2 = xdata


    # Raw counts
    line3a, = ax3.plot(xdata2, yrawtarget, 'go')
    line3b, = ax3.plot(xdata2, yrawcomp, 'ro') 
    fig.legend([line3a, line3b,], ['Target', 'Comp'], loc='upper left', \
               bbox_to_anchor=[0.9,0.68], shadow=True,\
               numpoints=1, prop={'size':8})

    # Differential Photometry
    line4, = ax4.plot(xdata2, ydflux, 'bo')
    fig.legend([line4,], ['Diff'], loc='upper left', \
               bbox_to_anchor=[0.9,0.48], shadow=True, \
               numpoints=1, prop={'size':8})

    # Seeing 
    if bjdswitch:
        ax5.set_xlabel('Time [hours] + JD '+str(lowestjd))
    else:
        ax5.set_xlabel('Time [minutes, now = 0]')
 
    ax5.set_ylabel('Pixels')
    line5, = ax5.plot(xdata2, yseeing, 'bo') 
    fig.legend([line5,], ['Seeing'], loc='upper left', \
               bbox_to_anchor=[0.9,0.27], shadow=True, \
               numpoints=1, prop={'size':8})

    # plot the error bars and collect them in
    # a new variable, so that they can be erased later
    errorlines3a = []
    errorlines3b = []
    errorlines4 = []
    for k in range(0,len(xdata2)):
        l, = ax3.plot([xdata2[k],xdata2[k]],[yrawtarget[k]-yrawtargeterr[k],\
                                           yrawtarget[k]+yrawtargeterr[k]], 'g-')
        if bjdswitch is False: 
            errorlines3a.append(l)
 
        l, = ax3.plot([xdata2[k],xdata2[k]],[yrawcomp[k]-yrawcomperr[k],\
                                             yrawcomp[k]+yrawcomperr[k]], 'r-')
        if bjdswitch is False:
            errorlines3b.append(l)

        l, = ax4.plot([xdata2[k],xdata2[k]],[ydflux[k]-ydfluxerr[k],\
                                             ydflux[k]+ydfluxerr[k]], 'b-')
        if bjdswitch is False:
            errorlines4.append(l)

    return line3a, line3b, line4, line5, xdata2, errorlines3a, errorlines3b, errorlines4, lowestjd

#############################
# append_datapoints
def append_datapoints(columns, xdata, yrawtarget, yrawtargeterr, \
                      yrawcomp, yrawcomperr, ydflux, ydfluxerr,\
                      yseeing, UTCtimes, now, bjdswitch, lowestjd):
    # here xdata is already an array

    # TBD: this to check if it is the same as previous points
    obsid         = columns[0]
    serverport    = columns[1]
    UTCdatetime   = columns[2]
    bandpass      = columns[3]
    BJD           = np.float64(columns[4]) 
    targetflux    = np.float32(columns[5])
    targetfluxerr = np.float32(columns[6])
    compflux      = np.float32(columns[7])
    compfluxerr   = np.float32(columns[8])
    seeing        = np.float16(columns[9])

    # TBD:
    # use just data that matches -filter argument
    # something like matches = [x for x in a if x=='a']

    UTCtimes.append(UTCdatetime)

    if bjdswitch:
        hours = np.float64(BJD-lowestjd)*24.
        xdata.append(hours)
    else:
        # UTCtimes is a list, stores all actual times for
        # recalculating minutes_before_now later
        newx = minutes_before_now(now, UTCdatetime)
        xdata.append(round(newx,3))
 
    yrawtarget.append(targetflux)
    yrawtargeterr.append(targetfluxerr)
    yrawcomp.append(compflux)
    yrawcomperr.append(compfluxerr)
    ydfluxval, ydfluxerrval = diff_photom(targetflux, targetfluxerr, \
                                          compflux, compfluxerr)
    ydflux.append(ydfluxval)
    ydfluxerr.append(ydfluxerrval)
    yseeing.append(seeing)
    
    print "RTPhoS: datapoint at time %s added" % str(UTCdatetime)
    return  xdata, yrawtarget, yrawtargeterr, \
                      yrawcomp, yrawcomperr, ydflux, ydfluxerr,\
                      yseeing, UTCtimes

##################
# update_dataplots
def update_plots(ax1, ax2, ax3, ax4, ax5, line3a, line3b, line4, line5, \
                 errorlines3a, errorlines3b, errorlines4, \
                 xdata, yrawtarget, yrawtargeterr, yrawcomp,\
                 yrawcomperr, ydflux, ydfluxerr, yseeing, oldn, args):

    if os.path.exists("t1.fits") and os.path.exists("c1.fits"):
        timage = fits.getdata('t1.fits')
        cimage = fits.getdata('c1.fits')
        # Thumbnail image settings
        # Since these are random images need to remove all 0 and negative values.
        medianintens = np.median(timage)
        timage[timage==0] = medianintens # Remove zero values
        timage = np.log(timage)          # Use for log scale plotting
        timage = np.absolute(timage)     # Remove negative values 
        medianintens = np.median(timage)
        cimage[timage==0] = medianintens # Remove zero values
        medianintens = np.median(cimage)
        cimage[timage==0] = medianintens # Remove zero values
        cimage = np.log(cimage)          # Use for log scale plotting
        cimage = np.absolute(cimage)     # Remove negative values 
        medianintens = np.median(cimage)
        timage[timage==0] = medianintens # Remove zero values
        cropmin = np.amin(timage)        # Scale is set based on target image
        cropmax = np.amax(timage)
        # Attempt for a reasonable thumbnail intensity scale 
        maxintens = ((cropmax-cropmin)/2.0)+cropmin
        # Target thumbnail
        ax1.plot([50,50],[0,100],'r:')             # Plot cross-hairs
        ax1.plot([0,100],[50,50],'r:')             #      -""-
        ax1.imshow(timage, cmap='gray', norm=LogNorm(vmin=cropmin, vmax=maxintens))

        # Comparison thumbnail
        ax2.plot([50,50],[0,100],'r:')             # Plot cross-hairs
        ax2.plot([0,100],[50,50],'r:')             #      -""-
        ax2.imshow(timage, cmap='gray', norm=LogNorm(vmin=cropmin, vmax=maxintens))

    # Raw counts
    line3a.set_xdata(xdata)
    line3a.set_ydata(yrawtarget)
    line3b.set_xdata(xdata)
    line3b.set_ydata(yrawcomp)

    # Differential flux
    line4.set_xdata(xdata)
    line4.set_ydata(ydflux)

    # Seeing plot
    line5.set_xdata(xdata)
    line5.set_ydata(yseeing)
    
    safe_lines3a = errorlines3a
    safe_lines3b = errorlines3b
    safe_lines4 = errorlines4
    errorlines3a = []
    errorlines3b = []
    errorlines4 = []
    # remove previous error bars
    for k in range(0,oldn):
        ax3.lines.remove(safe_lines3a[k])
        ax3.lines.remove(safe_lines3b[k])
        ax4.lines.remove(safe_lines4[k])
    # replot all new ones
    for k in range(0,len(xdata)):
        l, = ax3.plot([xdata[k],xdata[k]], \
                      [yrawtarget[k]-yrawtargeterr[k],\
                       yrawtarget[k]+yrawtargeterr[k]], 'g-')
        errorlines3a.append(l)
        l, = ax3.plot([xdata[k],xdata[k]], \
                      [yrawcomp[k]-yrawcomperr[k],\
                       yrawcomp[k]+yrawcomperr[k]], 'r-')
        errorlines3b.append(l)
        l, = ax4.plot([xdata[k],xdata[k]],\
                      [ydflux[k]-ydfluxerr[k],\
                       ydflux[k]+ydfluxerr[k]], 'b-')
        errorlines4.append(l)

    return ax1, ax2, ax3, ax4, ax5, line3a, line3b, line4, line5, \
           errorlines3a, errorlines3b, errorlines4

###################
#
def update_plots_jd(ax1, ax2, ax3, ax4, ax5, line3a, line3b, line4, line5, \
                    xdata, yrawtarget, yrawtargeterr, yrawcomp,\
                    yrawcomperr, ydflux, ydfluxerr, yseeing, oldn):

    if os.path.exists("t1.fits") and os.path.exists("c1.fits"):
        timage = fits.getdata('t1.fits')
        cimage = fits.getdata('c1.fits')
        # Thumbnail image settings
        # Since these are random images need to remove all 0 and negative values.
        medianintens = np.median(timage)
        timage[timage==0] = medianintens # Remove zero values
        timage = np.log(timage)          # Use for log scale plotting
        timage = np.absolute(timage)     # Remove negative values 
        medianintens = np.median(timage)
        cimage[timage==0] = medianintens # Remove zero values
        medianintens = np.median(cimage)
        cimage[timage==0] = medianintens # Remove zero values
        cimage = np.log(cimage)          # Use for log scale plotting
        cimage = np.absolute(cimage)     # Remove negative values 
        medianintens = np.median(cimage)
        timage[timage==0] = medianintens # Remove zero values
        cropmin = np.amin(timage)        # Scale is set based on target image
        cropmax = np.amax(timage)
        # Attempt for a reasonable thumbnail intensity scale 
        maxintens = ((cropmax-cropmin)/2.0)+cropmin
        # Target thumbnail
        ax1.plot([50,50],[0,100],'r:')             # Plot cross-hairs
        ax1.plot([0,100],[50,50],'r:')             #      -""-
        ax1.imshow(timage, cmap='gray', norm=LogNorm(vmin=cropmin, vmax=maxintens))

        # Comparison thumbnail
        ax2.plot([50,50],[0,100],'r:')             # Plot cross-hairs
        ax2.plot([0,100],[50,50],'r:')             #      -""-
        ax2.imshow(timage, cmap='gray', norm=LogNorm(vmin=cropmin, vmax=maxintens))

    #print "Xdata:", len(xdata)
    #print "Yraw", len(yrawtarget)
    # Raw counts
    line3a.set_xdata(xdata)
    line3a.set_ydata(yrawtarget)
    line3b.set_xdata(xdata)
    line3b.set_ydata(yrawcomp)

    # Differential flux
    line4.set_xdata(xdata)
    line4.set_ydata(ydflux)

    # Seeing plot
    line5.set_xdata(xdata)
    line5.set_ydata(yseeing)
                                                                                       
    npoints_current = len(xdata)
    # plot just the extra few errorbars
    for k in range(oldn,npoints_current):
        ax3.plot([xdata[k],xdata[k]], \
                 [yrawtarget[k]-yrawtargeterr[k],\
                  yrawtarget[k]+yrawtargeterr[k]], 'g-')
        ax3.plot([xdata[k],xdata[k]], \
                 [yrawcomp[k]-yrawcomperr[k],\
                  yrawcomp[k]+yrawcomperr[k]], 'r-')
        ax4.plot([xdata[k],xdata[k]],\
                 [ydflux[k]-ydfluxerr[k],\
                  ydflux[k]+ydfluxerr[k]], 'b-')
 
    return ax1, ax2, ax3, ax4, ax5, line3a, line3b, line4, line5

##########################################################
### MAIN CODE ###
##########################################################
plt.ion()

# parse the arguments from the command line
parser = argparse.ArgumentParser(description='RTPhoS client live plotter')

# mandatory params:
# 
parser.add_argument('filter', type=str, \
                    help='Plot data mathing the bandpass chosen')
parser.add_argument('port', metavar='port', type=str, \
                    help='TCP/IP port to listen')
# optional
parser.add_argument('-logfile', type=str, nargs='?', \
                    help='Optional string, file name for the log to be saved')
parser.add_argument('-npoints', type=int, nargs='?', \
                    help='If reading logfile, plot just last npoints')
parser.add_argument('-minnow', type=int, nargs='?', \
                    help='Set this to 1 (True) to plot X axis as minutes from current time. Default 0 (False)')

# TBD:
#parser.add_argument('-pastmin', type=int, \
#                    help='If reading lofile, also plot last pastmin minutes of data')
args = parser.parse_args()

#if args.minnow is None:
#    bjdswitch = True
bjdswitch = True
port = args.port

if args.logfile is not None:
    uselog = True
    if port is not None:
        print "!RTPhoS: you have specified both port and logfile to monitor!"
        print "This is not possible. Continuing with using the logfile."
else:
    uselog = False
    print "RTPhoS: no log file set. Plotting live data only."
    if port is None: 
        print "Set -port or -logfile!"

# bandpass filter must be at the beginning of the messages sent
if args.filter is not None:
    print "RTPhoS: filtering messages, bandpass %s" % args.filter
    bandpass_filter = "{\"bandpass\": \"" + args.filter
else:
    print "-filter option not set!"

#################################
# if using live updating log file
################################
if uselog:
######
     try: 
        # open the log file
        # mae a current copy of the logfile as "._tmp"
        # read all data that is currently there
        # returns also the file position up to which it read
        tempfile = args.logfile+"_tmp"
        command = "cp -f "+args.logfile+" "+tempfile
        os.system(command)
        # get current time
        now = datetime.utcnow()
        xdata, yrawtarget, yrawtargeterr, yrawcomp, yrawcomperr, \
        ydflux, ydfluxerr, yseeing, UTCtimes = init_datapoins(tempfile, now, bjdswitch)

        # remember the time when file has been read last
        lastdatafiletime = os.stat(args.logfile).st_mtime
        
        # TBD:
        # now filter out the last -npoints
        # .....

        # initialize the plots
        fig, ax1, ax2, ax3, ax4, ax5 = plot_initialize()

        # if any data already present in the log, plot it
        line3a, line3b, line4, line5, xdata2, errorlines3a, errorlines3b, errorlines4, lowestjd = \
                     init_plot_data(fig, ax1, ax2, ax3, ax4, ax5, \
                                  xdata, yrawtarget, yrawtargeterr, yrawcomp, \
                                       yrawcomperr, ydflux, ydfluxerr, yseeing, bjdswitch)
        fig.canvas.show()
        npoints_current = len(xdata) 
######### while
        while True: 
            # check if input data file has been modified
            time.sleep(5)
            datafiletime = os.stat(args.logfile).st_mtime
            if datafiletime != lastdatafiletime:
                # make a new copy of the log and get a difference
                command = "rm "+args.logfile+"_diff & cp "+args.logfile+" "+args.logfile+"_tmp2 "
                os.system(command)
                command = "diff "+args.logfile+"_tmp "+args.logfile+"_tmp2"
                command = command + "| grep \">\" | sed 's/>//g' > "+args.logfile+"_diff"
                os.system(command)
                command = "mv -f "+args.logfile+"_tmp2 "+args.logfile+"_tmp"
                os.system(command)

                now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                print now + " RTPhoS logfile was updated."
                # update time of last file update
                lastdatafiletime = datafiletime
                # update delta-times in minutes for old points
                xdata = []
                now = datetime.utcnow()
                if not bjdswitch:
                   for t in UTCtimes:
                       newtime = minutes_before_now(now,t)
                       xdata.append(round(newtime,3)) 
                   
                # invoke shell commands to make a 
                fp = open(args.logfile+"_diff",'r')
                for line in fp:
                    useline = line.strip()  
                    columns = useline.split()
                    # skip comment line
                    if columns[0] == '#RTPhoS:':
                        print "skipped:", columns
                    else:
                        oldn = npoints_current
                        
                        xdata2, yrawtarget, yrawtargeterr, yrawcomp, yrawcomperr, ydflux, ydfluxerr,\
                            yseeing, UTCtimes = append_datapoints(columns, xdata2, yrawtarget, yrawtargeterr, \
                                                                  yrawcomp, yrawcomperr, ydflux, ydfluxerr,\
                                                                  yseeing, UTCtimes, now, bjdswitch, lowestjd)
                                 
                        if bjdswitch:
                            ax1, ax2, ax3, ax4, ax5, line3a, line3b, line4, line5, = update_plots_jd(ax1, ax2, ax3, ax4, ax5, line3a, line3b,\
                                                                                       line4, line5, xdata2, yrawtarget, \
                                                                                       yrawtargeterr, yrawcomp, yrawcomperr, 
                                                                                       ydflux, ydfluxerr, yseeing,\
                                                                                       oldn)
                        else:
                            npoints_current, ax1, ax2, ax3, ax4, ax5, line3a, line3b, line4, line5, \
                                errorlines3a, errorlines3b, errorlines4 = update_plots(ax1, ax2, ax3, ax4, ax5, line3a, line3b,\
                                                                                       line4, line5, errorlines3a, \
                                                                                       errorlines3b,\
                                                                                       errorlines4, xdata, yrawtarget, \
                                                                                       yrawtargeterr, yrawcomp, yrawcomperr, 
                                                                                       ydflux, ydfluxerr, yseeing,\
                                                                                       oldn)

                        npoints_current = len(xdata2)
                        fig.canvas.show()
                        print "RTPhoS: currently displayed %s points" % str(npoints_current)
############# end if datafiletime !=
            else: 
                # Sleep 1 second before attempting to parse input data file(s) again
                time.sleep(2)
######### end while
###### end try
     except (KeyboardInterrupt, SystemExit):
        print "Process aborted by user."

###############################
## if live plotting from server
###############################

else:
     # Socket to talk to server
     context = zmq.Context()
     socket = context.socket(zmq.SUB)
     socket.connect ("%s" % port)
     socket.setsockopt(zmq.SUBSCRIBE,bandpass_filter)
     print "RTPhoS: Collecting updates from broadcasting server ", port

     # initialize the plots
     fig, ax1, ax2, ax3, ax4, ax5 = plot_initialize()

     fig.canvas.show()

     # Initialize data lists
     xdata=[]         # X-axis data (time, either past minutes or BJD)
     xdata2 = []      # in hours since lowest BJD
     yrawtarget = []  # Raw target counts
     yrawtargeterr=[] # Raw target error bars
     yrawcomp=[]      # Raw comparison counts
     yrawcomperr=[]   # Raw comparison error bars
     ydflux=[]        # Differential photometry counts
     ydfluxerr=[]     # Differential photometry error bars
     yseeing=[]       # Seeing data
     UTCtimes = []

     try: 
        while True:
            # get the message packet, which is a dictionary
            messagedata =json.loads(socket.recv())
            now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            
            obsid         = str(messagedata['obsid'])
            serverport    = str(messagedata['port'])
            sendtimeUTC   = str(messagedata['sendtimeUTC'])
            bandpass      = str(messagedata['bandpass'])
            UTCdatetime   = str(messagedata['UTCdatetime'])
            BJD           = np.float64(messagedata['BJD'])
            targetflux    = np.float32(messagedata['targetflux'])
            targetfluxerr = np.float16(messagedata['targetfluxerr'])
            compflux      = np.float32(messagedata['compflux'])
            compfluxerr   = np.float16(messagedata['compfluxerr'])
            seeing        = np.float16(messagedata['seeing'])
            timage       = np.array(messagedata['thumbnail1'])
            cimage       = np.array(messagedata['thumbnail2'])

            oldn = len(xdata)
            UTCtimes.append(UTCdatetime)
                        
            if bjdswitch:
                xdata.append(BJD)
            else:
                # UTCtimes is a list, stores all actual times for
                # recalculating minutes_before_now later
                newx = minutes_before_now(now, UTCdatetime)
                xdata.append(round(newx,3))

            yrawtarget.append(targetflux)
            yrawtargeterr.append(targetfluxerr)
            yrawcomp.append(compflux)
            yrawcomperr.append(compfluxerr)
            ydfluxval, ydfluxerrval = diff_photom(targetflux, targetfluxerr, \
                                                  compflux, compfluxerr)
            ydflux.append(ydfluxval)
            ydfluxerr.append(ydfluxerrval)
            yseeing.append(seeing)
            
            print "RTPhoS: message received %s" % now
            print "Sender: %s" % obsid
            print "Sent time UTC : %s" % sendtimeUTC
            print "Observation time UTC: %s" % UTCdatetime
            now = datetime.utcnow()
            sometime = datetime.strptime(UTCdatetime, "%Y-%m-%d|%H:%M:%S")
            elapsedTime = sometime - now
            secondold = elapsedTime.total_seconds()
            print "Data is %s seconds old." % str(round(abs(secondold),1))

            # if this is the first data point
            print "Oldn:", oldn
            if oldn == 0:
                    line3a, line3b, line4, line5, xdata2, errorlines3a, errorlines3b, errorlines4, lowestjd = \
                     init_plot_data(fig, ax1, ax2, ax3, ax4, ax5, \
                                  xdata, yrawtarget, yrawtargeterr, yrawcomp, \
                                       yrawcomperr, ydflux, ydfluxerr, yseeing, bjdswitch)
            else:
                if bjdswitch:
                    hours = (BJD-lowestjd)*24.0
                    xdata2.append(hours)

                    ax1, ax2, ax3, ax4, ax5, line3a, line3b, line4, line5, = update_plots_jd(ax1, ax2, ax3, ax4, ax5, line3a, line3b,\
                                                                                         line4, line5, xdata2, yrawtarget, \
                                                                                         yrawtargeterr, yrawcomp, yrawcomperr, 
                                                                                         ydflux, ydfluxerr, yseeing,\
                                                                                         oldn)

     except (KeyboardInterrupt, SystemExit):
         print "Process aborted by user."

