#!/usr/bin/python
# M. Bogosavljevic, AOB, June 2015

def write_rtphos_defaults( rtphos, pathdefs, pathlog, pathans, pathdata, pathbias, pathdark, pathflat, stringbias, stringdark, stringflat, namebias, namedark, nameflat, psfs, stars, cprefix, sradius, aradius, cradius, starnumber, skyskew, skyfit, gain, linlevel, satslevel, verbose, tsleep, fsleep, tfile, cfile, port, obsid, band, sslep, dpath, tdfile, public, cdfile, tfits, cfits, lport, lfile, lband, ltfits, lcfits ):

    defs_file = open(pathdefs, "w")

    defs_file.write(pathdefs   + "     # File containing new settings\n")
    defs_file.write(pathdata   + "     # Raw data frames\n")		
    defs_file.write(pathbias   + "     # Bias frames for calibration\n")	
    defs_file.write(pathdark   + "     # Dark frames for calibration\n")
    defs_file.write(pathflat   + "     # Flat frames for calibration\n")	
    defs_file.write(stringbias + "     # bias frames wildcard\n")		
    defs_file.write(stringdark + "     # dark frames wildcard\n")		
    defs_file.write(stringflat + "     # flat frames wildcard\n")		
    defs_file.write(namebias   + "     # Masterbias filename\n")		
    defs_file.write(namedark   + "     # Masterdark filename\n")		
    defs_file.write(nameflat   + "     # Masterflat filename\n")		
    defs_file.write(psfs       + "     # list of objects for psf \n")
    defs_file.write(stars      + "     # full list of objects \n")	
    defs_file.write(cprefix    + "     # Prefix for calibrated files\n")
    defs_file.write(sradius    + "     # Search radius for optimal centroiding\n")	
    defs_file.write(aradius    + "     # Aperature radius for photometry\n")	
    defs_file.write(cradius    + "     # PSF Clipping radius for optimal photom\n")	
    defs_file.write(starnumber + "     # Star number to optimize photometry for\n")	
    defs_file.write(skyskew    + "     # Sky profile skew for optimal photometry\n")	
    defs_file.write(skyfit     + "     # Sky fitting switch for optimal photometry\n")	
    defs_file.write(gain       + "     # Instrument Gain e-/ADU\n")	
    defs_file.write(linlevel   + "     # Instrument linearity level (counts)\n")	
    defs_file.write(satslevel  + "     # Instrument saturatoin level (counts) e-/ADU\n")	
    defs_file.write(verbose    + "     # Verbose output switch\n")	
    defs_file.write(tsleep     + "     # Sleep time before checking for new frames [s]\n")	

    defs_file.close()
    print "Wrote :"+pathdefs

    defs_file = open(pathans, "w")
    
    defs_file.write("param rtphosdef\n")
    defs_file.write(" rtphos  entry {RTPhoS location}       \""        +   rtphos + "\"\n" )
    defs_file.write(" pathdefs    entry {RTPhoS settings} \""              +   pathdefs   + "\"\n" )
    defs_file.write(" pathlog     entry {RTPhoS logfile}               \"" +   pathlog    + "\"\n" )
    defs_file.write(" pathans     entry {DS9 Analysis settings} \""        +   pathans    + "\"\n" )
    defs_file.write(" pathdata    entry {Raw data frames}	       \"" +   pathdata   + "\"\n" )   
    defs_file.write(" pathbias    entry {Bias frames for calibration}  \"" +   pathbias   + "\"\n" )    
    defs_file.write(" pathdark    entry {Dark frames for calibration}  \"" +   pathdark   + "\"\n" )    
    defs_file.write(" pathflat    entry {Flat frames for calibration}  \"" +   pathflat   + "\"\n" )   
    defs_file.write(" stringbias  entry {bias frames wildcard}	       \"" +   stringbias + "\"\n" )  
    defs_file.write(" stringdark  entry {dark frames wildcard}	       \"" +   stringdark + "\"\n" )  
    defs_file.write(" stringflat  entry {flat frames wildcard}	       \"" +   stringflat + "\"\n" )  
    defs_file.write(" namebias    entry {Masterbias filename}	       \"" +   namebias   + "\"\n" )  
    defs_file.write(" namedark    entry {Masterdark filename}	       \"" +   namedark   + "\"\n" )       
    defs_file.write(" nameflat    entry {Masterflat filename}	       \"" +   nameflat   + "\"\n" )
    defs_file.write(" psfs        entry {list of objects for psf}      \"" +   psfs       + "\"\n" )
    defs_file.write(" stars       entry {full list of objects}         \"" +   stars      + "\"\n" )
    defs_file.write(" cprefix     entry {Prefix for calibrated files}   \"" +  cprefix    + "\"\n" )
    defs_file.write(" sradius     entry {Search radius for optimal centroiding}     "   +  sradius + "\n" )
    defs_file.write(" aradius     entry {Aperature radius for photometry}           "   +  aradius + "\n" )
    defs_file.write(" cradius     entry {PSF Clipping radius for optimal photom}    "   +  cradius + "\n" )    
    defs_file.write(" starnumber  entry {Star number to optimize photometry for}    "   +  starnumber + "\n" )  
    defs_file.write(" skyskew     entry {Sky profile skew for optimal photometry}   "   +  skyskew   + "\n" )  
    defs_file.write(" skyfit      entry {Sky fitting switch for optimal photometry} "   +  skyfit    + "\n" )
    defs_file.write(" gain        entry {Instrument Gain e-/ADU}         " +  gain       + "\n" )
    defs_file.write(" linlevel    entry {Instrument linearity level (counts)}       " +  linlevel    + "\n" )
    defs_file.write(" gain        entry {Instrument saturation level (counts)}      " +  satslevel   + "\n" )
    defs_file.write(" verbose     checkbox {Verbose output switch}       " +  verbose    + "\n" )    
    defs_file.write(" tsleep   entry {Sleep time before checking for new frames [s]} " + tsleep + "\n" )
    defs_file.write("endparam\n")
    defs_file.write("\n")

    defs_file.write("param fakedata \n")
    defs_file.write(" fsleep entry {Fake data time period [sec]} " + fsleep + "\n")
    defs_file.write(" tfile  entry {Fake data target file} " + tfile + "\n")
    defs_file.write(" cfile  entry {Fake data comparison file} "+ cfile + "\n" )
    defs_file.write("endparam\n")
    defs_file.write("\n")

    defs_file.write("param server \n")
    defs_file.write(" port   entry {Broadcast to TCP/IP port}" + port + "\n")
    defs_file.write(" obsid  entry {Observatory identifier} " + obsid + "\n")
    defs_file.write(" band   entry {Filter passband identifier} " + band + "\n")
    defs_file.write(" ssleep entry {Sleep [sec]} " + ssleep + "\n") 
    defs_file.write(" dpath  entry {Data files path [default ./]} " + dpath + "\n")
    defs_file.write(" tdfile entry {Target object data file} " + tdfile + "\n")
    defs_file.write(" public checkbox {Check for broadcasting target flux} " + public + "\n")
    defs_file.write(" cdfile entry {Comparison star data file [optional]} " + cdfile + "\n")
    defs_file.write(" tfits  entry {Target stamp fits file [optional]} " + tfits + "\n")
    defs_file.write(" cfits  entry {Comparison stamp fits file [optional]} " + cfits + "\n")
    defs_file.write("endparam \n")

    defs_file.write("param logger \n")
    defs_file.write(" lport entry  {Receive from TCP/IP port} " + lport + "\n")
    defs_file.write(" lfile entry  {Log file name [optional]} " + lfile +  "\n")
    defs_file.write(" lband entry  {Filter passband to log [optional]} " + lband + "\n")
    defs_file.write(" ltfits entry {Target stamp fits file [optional]} " + ltfits + "\n")
    defs_file.write(" lcfits entry {Comparison stamp fits file [optional]} " + lcfits + "\n")
    defs_file.write("endparam \n")

    defs_file.write("RTPhoS Settings\n")
    defs_file.write("*\n")
    defs_file.write("menu\n")
    defs_file.write("$param(rtphosdef);  $param(gifparams); $param(sourceparams); xterm -hold -sb -sl 2000 -e bash -c \"$rtphos/write_rtphos_defaults.py $rtphos $pathdefs $pathlog $pathans $pathdata $pathbias $pathdark $pathflat $stringbias $stringdark $stringflat $namebias $namedark $nameflat $stars $psfs $cprefix $sradius $aradius $cradius $starnumber $skyskew $skyfit $gain $linlevel $satslevel $verbose $tsleep $fsleep $tfile $cfile $port $obsid $band $sslep $dpath $tdfile $public $cdfile $tfits $cfits $lport $lfile $lband $ltfits $lcfits \" \n")
    defs_file.write("\n")

    defs_file.write("RUN RTPhoS data reduction\n")
    defs_file.write("*\n")
    defs_file.write("menu\n")
    defs_file.write("$param(rtphosdef); xterm -hold -sb -sl 2000 -e \"tcsh -c \\\"python -u $rtphos/run_rtphos.py $xpa_method $pathdefs |& tee $pathlog\\\" \" \n")
    defs_file.write("\n")

    defs_file.write("RTPhoS Fake Data Module\n")
    defs_file.write("*\n")
    defs_file.write("menu\n")
    defs_file.write("$param(rtphosdef); $param(fakedata); xterm -hold -sb -sl 2000 -e tcsh -c \"$rtphos/RTPhoS_fakedata.py $fsleep $tfile $cfile \" \n")
    defs_file.write("\n")

    defs_file.write("RTPhoS Server Module\n")
    defs_file.write("*\n")
    defs_file.write("menu\n")
    defs_file.write("$param(rtphosdef); $param(server); xterm -hold -sb -sl 2000 -e tcsh -c \" $rtphos/RTPhoS_server_module.py $port $obsid $band $tdfile $tsleep -public $public --path $dpath --compfile $cdfile --thumbtarget $tfits --thumbcomp $cfits \" \n")
    defs_file.write("\n")

    defs_file.write("RTPhoS Client Logger Module\n")
    defs_file.write("*\n")
    defs_file.write("menu\n")
    defs_file.write(" $param(rtphosdef); $param(logger); xterm -hold -sb -sl 2000 -e tcsh -c \" $rtphos/RTPhoS_client_logger.py $lport -filter $lband -logfile $lfile -thumbtarget $ltfits -thumbcomp $lcfits \" \n")
    defs_file.write("\n")

    defs_file.close()

    print "Wrote :"+pathans


if  __name__ == "__main__":

    import sys
    rtphos        = sys.argv[1]
    pathdefs      = sys.argv[2]
    pathlog       = sys.argv[3]
    pathans        = sys.argv[4]
    pathdata       = sys.argv[5]
    pathbias       = sys.argv[6]
    pathdark       = sys.argv[7]
    pathflat       = sys.argv[8]
    stringbias     = sys.argv[9]
    stringdark     = sys.argv[10]
    stringflat     = sys.argv[11]
    namebias       = sys.argv[12]
    namedark       = sys.argv[13]
    nameflat       = sys.argv[14]
    psfs           = sys.argv[15]
    stars          = sys.argv[16]
    cprefix        = sys.argv[17]
    sradius        = sys.argv[18]
    aradius        = sys.argv[19]
    cradius        = sys.argv[20]
    starnumber     = sys.argv[21]
    skyskew        = sys.argv[22]
    skyfit         = sys.argv[23]
    gain           = sys.argv[24]
    linlevel       = sys.argv[25]
    satslevel      = sys.argv[26]
    verbose        = sys.argv[27]
    tsleep         = sys.argv[28]
    fsleep         = sys.argv[29]        
    tfile          = sys.argv[30]
    cfile          = sys.argv[31]
    port           = sys.argv[32]
    obsid          = sys.argv[33]
    band           = sys.argv[34]
    sslep          = sys.argv[35]
    dpath          = sys.argv[36]
    tdfile         = sys.argv[37]
    public         = sys.argv[38]
    cdfile         = sys.argv[39]
    tfits          = sys.argv[40]
    cfits          = sys.argv[41]
    lport          = sys.argv[42]
    lfile          = sys.argv[43]
    lband          = sys.argv[44]
    ltfits         = sys.argv[45]
    lcfits         = sys.argv[46]

    write_rtphos_defaults( rtphos, pathdefs, pathlog, pathans, pathdata, pathbias, pathdark, pathflat, stringbias, stringdark, stringflat, namebias, namedark, nameflat, psfs, stars, cprefix, sradius, aradius, cradius, starnumber, skyskew, skyfit, gain, linlevel, satslevel, verbose, tsleep, fsleep, tfile, cfile, port, obsid, band, sslep, dpath, tdfile, public, cdfile, tfits, cfits, lport, lfile, lband, ltfits, lcfits)

