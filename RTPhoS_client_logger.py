#!/usr/bin/env python
# RTPhoS client module which receives data from server module
# and logs into a flat ascii file.

# Usage:
# RTPhoS_client_logger.py [port] 
# {plus optional params and switches, do -h or --help to get options}

# port should be a string like tcp://localhost:port or other address
# Examples:
# RTPhoS_client_logger.py tcp://localhost:5556 
# RTPhoS_client_logger.py tcp://92.96.59.99:5556 -logfile mylog.dat

import os
import sys
import argparse
import zmq
import time
from datetime import datetime
import numpy as np
from astropy.io import fits
import json

##########################################################
### MAIN CODE ###
##########################################################

# parse the arguments from the command line
parser = argparse.ArgumentParser(description='RTPhoS client logger')

parser.add_argument('port', metavar='port', type=str, \
                    help='TCP/IP port to listen')
# optional
parser.add_argument('-logfile', type=str, nargs='?', \
                    help='Optional string, file name for the log to be saved')
parser.add_argument('-filter', type=str, nargs='?', \
                    help='If set, will only save lines mathing the bandpass given as string')
parser.add_argument('-thumbtarget', metavar='thumbtarget', type=str, nargs='?', \
                    help='If set to a string, will save thumbnails of target object [string.fits]')
parser.add_argument('-thumbcomp', metavar='thumpcomp', type=str, nargs='?', \
                    help='If set to a string, will save thumbnails of comparison star [string.fits])')

args = parser.parse_args()

# bandpass filter must be at the beginning of the messages sent
if args.filter is not None:
    print "RTPhoS: filtering messages, bandpass %s" % args.filter
    bandpass_filter = "{\"bandpass\": \"" + args.filter
else:
    bandpass_filter = ""

print "RTPhoS: Collecting updates from broadcasting server %s" % args.port
# Socket to listen to server
context = zmq.Context()
socket = context.socket(zmq.SUB)
socket.connect ("%s" % args.port)
socket.setsockopt(zmq.SUBSCRIBE,bandpass_filter)


now = datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")
if args.logfile is not None: 
    outfile = args.logfile
    print 'Logfile:', outfile
else:
    outfile = "RTPhoS_log_UTC"+ now.replace(" ", "_")+".txt"
    print 'Logfile:', outfile
print 'Filter: ', args.filter
print 'Thumbnail target name:', args.thumbtarget
print 'Thumbnail comparison name:', args.thumbcomp

myfile = open(outfile, "a")
myfile.write("#RTPhoS: Log start UTC time %s \n" % now)
myfile.close()

try: 
    while True:
        myfile = open(outfile, "a")
        # get the message packet, which is a dictionary
        messagedata =json.loads(socket.recv())
        now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        # as an example
        obsid         = str(messagedata['obsid'])
        serverport    = str(messagedata['port'])
        sendtimeUTC   = str(messagedata['sendtimeUTC'])
        sendtimeUTC   = sendtimeUTC.replace(" ","|")
        bandpass      = str(messagedata['bandpass'])
        UTCdatetime   = str(messagedata['UTCdatetime'])
        BJD           = str(messagedata['BJD'])
        targetflux    = str(messagedata['targetflux'])
        print messagedata['targetfluxerr']
        targetfluxerr = str(round(np.float(messagedata['targetfluxerr']),3))
        compflux      = str(messagedata['compflux'])
        compfluxerr   = str(messagedata['compfluxerr'])
        seeing        = str(messagedata['seeing'])
      
        print "RTPhoS: message received %s" % now
        print "Sender: %s" % obsid
        print "Sent time UTC : %s" % sendtimeUTC
        print "Observation time UTC: %s" % UTCdatetime

        table_data = obsid.ljust(20) + serverport.ljust(28) + UTCdatetime.ljust(25) + \
                     bandpass.ljust(5) + \
                     BJD.ljust(16) + targetflux.ljust(12) +  targetfluxerr.ljust(10) + \
                     compflux.ljust(12) + compfluxerr.ljust(10) + seeing.ljust(4) + '\n'

        # Write broadcast to text file
        myfile.write(table_data)
        myfile.close()
        
        # Write fits images. Only the latest transmitted images are saved.
        # Target star fits
        timage = messagedata['thumbnail1']
        if isinstance(timage, list):
                timage       = np.array(messagedata['thumbnail1'])                        
                hdu=fits.PrimaryHDU(timage)
                hdu.scale(type='float32', bzero=32768, bscale=1)
                hdu.writeto(args.thumbtarget, clobber='True')
        
        # Comparison star fits         
        cimage = messagedata['thumbnail2']        
        if isinstance(cimage, list):            
                cimage       = np.array(messagedata['thumbnail2'])                
                hdu=fits.PrimaryHDU(cimage)
                hdu.scale(type='float32', bzero=32768, bscale=1)
                hdu.writeto(args.thumbcomp, clobber='True')
        
        now = datetime.utcnow()
        sometime = datetime.strptime(UTCdatetime, "%Y-%m-%d|%H:%M:%S")
        elapsedTime = sometime - now
        secondold = elapsedTime.total_seconds()

        print "Data is %s seconds old." % str(round(abs(secondold),1))
        print "--- Message logged."

except (KeyboardInterrupt, SystemExit):
    print "Process aborted by user."
    pass
