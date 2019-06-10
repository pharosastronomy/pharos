#!/usr/bin/env python
# PHAROS Broadcast module module which watches output photometry files
# and broadcasts new data as it arrives over TCP/IP using the ZeroMQ library.

# Command line usage:
# PHAROS_Broadcast.py [port] [obsid] [band] [targetfile] [tsleep] 
# {plus optional params and switches, do -h or --help to get options}

import os
import sys
import argparse
import zmq
import time
from datetime import datetime
import numpy as np
from astropy.io import fits
import numpy as np
import json
from collections import OrderedDict

#########################################################################################
# file_len
# Returns the actual number of lines in a text file
# do not have lines with just whitespace!
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

#########################################################################################
# datatext_to_message
# Parses the chosen line of the text file and returns a dictionary
# the convention for nline is that lines are numbered from 1 to file_len(textfile)
# Used to read target data or comparison star data files with the format as:
#   No.    UTCdatetime       BJD            terr[s]  Flux       Flux_err    seeing flag  filename
#    1  1990-01-01|00:00:00  2447892.500058   5.00   33368.6875 647.047546   5.46   O     gauss01.fits 
def datatext_to_message(textfile, nline, public):

    with open(textfile,'r') as fp:
         i = 1
         for line in fp:
             if i == nline:
                 useline = line
                 break
             else:
                 i = i + 1

         useline = useline.strip()   
         columns = useline.split()
         thisentry    = columns[0]
         UTCdatetime  = columns[1]
         BJD          = columns[2] 
         flux         = columns[4]
         fluxerr      = columns[5]
         seeing       = columns[6]

         # Create a dictionary of the read values based on the privacy level
         if public==0:
            quality = 100.0*float(fluxerr)/float(flux)
            part_message = {"public":public, "linenumber":i, "UTCdatetime":UTCdatetime, \
                            "quality": quality, "seeing":seeing} 
         if public==1:
            part_message = {"public":public, "linenumber":i, "UTCdatetime":UTCdatetime, "BJD":BJD, \
                            "flux":flux, "fluxerr":fluxerr, "seeing":seeing} 

         print part_message
 
         return part_message

#########################################################################################
# format_and_broadcast
# Using the command line arguments in args, 
# check the public switch to determine how much data to broadcast. Check to see
# if the required data and thumbnails are available and constract the message to
# broadcast. If data is missing then NaN values will be sent. 

def format_for_broadcast(args, t_part_message, nline):
# nline here serves only to find the same line in comparison star data file
# if that is the preference set 

    compdone = False
    if args.thumbtarget is not None:
        print "Reading target fits file: %s" % args.thumbtarget
        timage_list = fits.getdata(args.thumbtarget)
        timage_list = timage_list.tolist()
    else:
        timage_list = 'NaN'
        
    if args.compfile is not None:
       c_part_message = datatext_to_message(args.compfile, nline)
       compdone = True
    if args.thumbcomp is not None:
        print "Reading comp fits file: %s" % args.thumbcomp
        cimage_list = fits.getdata(args.thumbcomp)
        cimage_list = cimage_list.tolist()
    else:
        cimage_list = 'NaN'

    UTCdatetime =  t_part_message['UTCdatetime']
    BJD = t_part_message['BJD']
    if args.public == 1:
        targetflux = t_part_message['flux']
        targetfluxerr = t_part_message['fluxerr']
    else:
        targetflux = 'NaN'
        # make target flux error a percentage in this case
        targetfluxerr = float(t_part_message['fluxerr']) / float(t_part_message['flux']) * 100.

    seeing = t_part_message['seeing']

    if compdone:
        compflux = c_part_message['flux']
        compfluxerr = c_part_message['fluxerr']
    else:
        compflux = 'NaN'
        compfluxerr = 'NaN'

    sendtimeUTC = datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")
    # create an ordered dictionary based on the privacy switch
    message = OrderedDict ( [("bandpass", args.band), ("sendtimeUTC", sendtimeUTC), ("obsid", args.obsid), \
                             ("port", args.port), \
                             ("UTCdatetime", UTCdatetime), ("BJD",BJD), ("targetflux", targetflux), \
                             ("targetfluxerr", targetfluxerr), ("compflux", compflux), \
                             ("compfluxerr", compfluxerr), ("seeing", seeing), \
                             ("thumbnail1", timage_list), ("thumbnail2", cimage_list)] )
 
    # convert to json message and broadcast 
    jsonmessage = json.dumps(message)
    return jsonmessage

##########################################################
### MAIN CODE ###
##########################################################

# parse the arguments from the command line
parser = argparse.ArgumentParser(description='PHAROS broadcasting module')

# required params
parser.add_argument('port', metavar='port', type=str, \
                    help='TCP/IP port to be used for the broadcast')
parser.add_argument('obsid', metavar='obsid', type=str, \
                    help='Observatory ID string to be used')
parser.add_argument('band', metavar='band', type=str, \
                    help='String stating which photometric band is being broadcast')
parser.add_argument('targetfile', metavar='targetfile', type=str, \
                    help='Target object data file to be monitored for updates and broadcast.')
parser.add_argument('tsleep', metavar='tsleep', type=float, \
                    help='Sleep time between checks of target data file')

# boolean switch: set it in order to broadcast the fluxes, otherwise they will be private
parser.add_argument('-public', type=int, default=0, \
                    help='Set this switch to 1 to broadcast the target flux. Default 0 (FALSE)')
# optional
parser.add_argument('--compfile', metavar='compfile', type=str, nargs='?', default='empty', \
                    help='Comparison star data to be monitored for updates and broadcast.')
parser.add_argument('--thumbtarget', metavar='thumbtarget', type=str, nargs='?', default='empty', \
                    help='fits thumbnail image of the target object')
parser.add_argument('--thumbcomp', metavar='thumbcomp', type=str, nargs = '?', default='empty', \
                    help='fits thumbnail image of the comparison star')
parser.add_argument('--path', metavar='path', type=str, nargs='?', default='./', \
                    help='path to data files to be read [default: current dir]')

args = parser.parse_args()

if args.compfile == 'empty':
    args.compfile = None

if args.thumbtarget == 'empty':
    args.thumbtarget = None

if args.thumbcomp == 'empty':
    args.thumbcomp = None

if args.path == None: 
    args.path = './'

#print args.compfile
#print args.thumbtarget
#print args.thumbcomp
#print args.path
print args.public

# connect to a publishing port
port = int(args.port)
context = zmq.Context()
socket = context.socket(zmq.PUB)
socket.bind("tcp://*:%s" % port)

lastdatafiletime = 0
nline = 1
firsttime = True 

## main loop ##
try:
  while True:
      # check if input data file has been modified
      datafiletime = os.stat(args.targetfile).st_mtime
      if datafiletime != lastdatafiletime:
          # update time of last file update
          lastdatafiletime = datafiletime
          # check the new length of the file 
          filelen = file_len(args.targetfile)
          #print "File has", filelen
          if firsttime:
              # read the just last line if it is the first time
              # otherwise, read all lines that are new since last visit
              print "PHAROS: Input file has", filelen, " existing lines"
              nline = filelen
              t_part_message = datatext_to_message(args.targetfile, nline, args.public)
              jsonmessage = format_for_broadcast(args,t_part_message,nline)
              socket.send(jsonmessage)
              print "PHAROS: Message sent ", datetime.now()
              #print jsonmessage
              firsttime = False
              # move counter to next line for next read
              nline = nline + 1
          else:
              if nline <= filelen:
                  for row in range(nline,filelen+1):
                      t_part_message = datatext_to_message(args.targetfile, row, args.public)
                      jsonmessage = format_for_broadcast(args,t_part_message,row)
                      socket.send(jsonmessage)
                      now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                      print "PHAROS: Message sent ", now
                      #print jsonmessage
                      # set nline to one line more for next read
                      nline = row + 1                  
      else: 
          # Sleep before attempting to parse input data file(s) again
          time.sleep(args.tsleep)
          now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
          print "PHAROS: Checking input files again...", now
except (KeyboardInterrupt, SystemExit):
  print "PHAROS: Data broadcast aborted by user."
  pass
