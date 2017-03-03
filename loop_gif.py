#!/usr/bin/python
# 
"""
M. Bogosavljevic, AOB, June 2015
this should listen to a directory (obtained from ref_filename)
subdirecrtory /png 
for images containing namestring wildcard
and find the last nframes there
and then create an animated gif from those last nframes
"""
import os, time
import glob

#===============================================================================
# Select last N frames of FITS files based on a filename wildcard.
#===============================================================================
def makefilelist(path,wildcard,nframes):

    # Move to the specified directory
    os.chdir(path)
    # Search for all files in the current directory that satisfy the wildcard provided
    filelist = sorted(glob.glob(wildcard))
    # print ("I see: ",filelist)
    lastn    = filelist[-nframes:]
    return lastn

#===============================================================================
# Main program
#===============================================================================
def loop_gif(ref_filename,wildcard,nframes,tsleep):

    # convert to int if string
    nframes = int(nframes)
    tsleep  = int(tsleep)

    # sanity check
    if nframes>30:
        print "Be reasonable and ask for less than 30 frames!"
        return

    path, filename = os.path.split(ref_filename)    
    png_dir = path + '/png'

    lastn    = makefilelist(png_dir,wildcard,nframes)
    mylist   = ' '.join(lastn)
    outfile = 'last'+str(nframes)+'loop.gif'

    if len(lastn)>=nframes:
        print "Creating animated gif from frames found:"
        mycmd = 'convert -delay 50 -loop 0 '+mylist+' '+outfile
        print mycmd
        os.system(mycmd)
        print "Done."
    else: 
        print "Not enough frames yet: "+str(len(lastn))+"<"+str(nframes)
    
    print "Checking for new files every "+str(tsleep)+" seconds"

    try:
        while 1:
            newlastn  = makefilelist(png_dir,wildcard,nframes)
            newlist   = ' '.join(newlastn)
            if (newlist != mylist) and len(newlastn)>=int(nframes):
                mycmd = 'convert -delay 50 -loop 0 '+newlist+' '+outfile
                print mycmd
                os.system(mycmd)
                mylist = newlist
                print "Waiting for new frames ... "
            time.sleep(tsleep)   # Wait for tsleep seconds before repeating
    except KeyboardInterrupt:
        pass



if __name__ == "__main__":

   import sys
   ref_filename    = sys.argv[1]
   wildcard        = sys.argv[2]
   nframes         = sys.argv[3]
   tsleep          = sys.argv[4]
   
   loop_gif(ref_filename,wildcard,nframes,tsleep)
