import sys
import zmq
import json
import numpy as np

# the mandatory arguments are port and bandpass 
# with bandpass argument you filter just the one bandpass you want to receive
# with this instance of the client

port =  sys.argv[1]

# bandpass filter must be at the beginning of the messages sent
if len(sys.argv) > 1:
 bandpass_filter = "{\"bandpass\": \"" + sys.argv[2]
else:
 bandpass_filter = "{\"bandpass\": \"R\""

print "Initiate with bandpass filter: ", sys.argv[2]

# should be either "tcp://localhost:port" or
# if not localhost, check inet address with ifconfig, 
# and then use something like
# tcp://192.168.1.103:5556 or tcp://92.96.59.99:port etc

# Socket to talk to server
context = zmq.Context()
socket = context.socket(zmq.SUB)
socket.connect ("%s" % port)
socket.setsockopt(zmq.SUBSCRIBE,bandpass_filter)
print "Collecting updates from broadcasting server ", port


try: 
    while True:
        # get the message packet, which is a dictionary
        messagedata =json.loads(socket.recv())
        print messagedata
        # as an example
        obsid         = messagedata['obsid']
        serverport    = messagedata['port']
        localUTtime   = messagedata['localUTtime']
        BJD           = messagedata['BJD']
        bandpass      = messagedata['bandpass']
        targetflux    = messagedata['targetflux']
        targetfluxerr = messagedata['targetfluxerr']
        compflux      = messagedata['compflux']
        compfluxerr   = messagedata['compfluxerr']
        seeing        = messagedata['seeing']      
        scidata       = np.array(messagedata['thumbnail1'])

        # for fun print to terminal
        print "obsid", obsid
        print "port", serverport
        print "localtime", localUTtime
        print "BJD", BJD
        print "bandpass", bandpass
        print "targetflux", targetflux
        print "targetfluxerr", targetfluxerr
        print "compflux", compflux
        print "compfluxerr", compfluxerr
        print "seeing", seeing
        print "scidata", scidata

except (KeyboardInterrupt, SystemExit):
    print "Process aborted by user."
    pass

