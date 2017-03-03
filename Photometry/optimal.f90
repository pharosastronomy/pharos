program optimal_phot

!=============================================================================== 
! VERSION FOR RTPhoS - No User INPUT, all input controlled from run_RTPhos.py
!===============================================================================
! 
! Program to calculate the flux of stars in CCD images using Tim Naylor's optimal 
! extraction routines. For a thorough explanation of the algorithms read the 
! following two papers:
!
! "An Optimal Extraction Algorithm for Imaging Photometry"
! Naylor Tim, MNRAS, 296, 399, 1998
!
! "Optimal Photometry for Colour-Magnitude Diagrams and its Application on NGC 2547"
! Naylor, Tim et al. MNRAS, 335, 291, 2002
!
! Requirements:
! opt_extr.f90 - All the required routines for optimal and aperature photometry
! marq.f90     - Computes Marquardt routines needed for profile fitting
! subs.f90     - Various array and file sorting routines
! cfitsio      - NASA's FITS file manipulation library needs to be installed.
!
! Compiles with:
! edit optimal.csh and source it. (Edit needed to provide path for cfitsio library)
!
! Alternatively (I use g95 but any Fortran compiler will do):
! g95 -c opt_extr.f90 marq.f90
! g95 -o optimal optimal.f90 marq.o opt_extr.o -L/<Local Path>/cfitsio -lcfitsio
!
! Required improvements:
!
!

  use opt_extr

  implicit none

!----------------------------------------------------------------------------

! Variable declarations
! ------------------------------------------------------------------------------

! Various counters & test variables
  integer :: i, j, fnum, num
  integer :: progcnt, stepcnt

  real :: searchrad

  logical :: verbose

  character (len=1)  :: verbosein
  character (len=3)  :: suf  
  character (len=80) :: wformat
  character (len=80) :: dummy
  character (len=80) :: input(13)

! Variables for opening FITS files and reading the image data
  integer :: status, unit, readwrite, blocksize, naxis, group, nfound
  integer :: naxes(2),fpix(2), lpix(2), inc(2)
  integer :: low(2), high(2)

  real :: datamin, datamax, nullval
  real, allocatable :: numflag(:,:)  ! **** This is the image pixel flags array
  real, allocatable :: data(:,:)     ! **** This is the image data array

  logical :: anynul

  character (len=50) :: filename, flagfile
! This is the image pixel flag array in original opphot format
  character, allocatable, dimension(:,:) :: pix_flg 

! Variables regarding PSF stars
  integer :: npsf, ipsf, nfit

  real, dimension(3) :: shape_par !The 3 params defining the shape of the PSF 
  real, allocatable, dimension(:) :: dpsf, xpos0, ypos0  
  real, allocatable, dimension(:,:,:) :: newpsfpos

  double precision, allocatable, dimension(:,:) :: psfstars
  
  character (len=50) :: psfpos
  character (len=15), allocatable, dimension(:) :: psfnames

! Variables used for clipping the PSF star
  real :: clip_fwhm, cliprad, fwhm
  
  character(len=10) :: clipanswer, distanswer 

! Variables regarding Target stars
  integer :: istar, nstar, nframes, ntimes
     
  logical :: comp, iftimes
  real :: xcomp, ycomp
 
  real :: xpos, ypos, badpix, posfix
  real :: optflux, opterror, optnrm, peak, xfit, yfit, xerr, yerr
  real, allocatable, dimension(:) :: dpos, seeing
  real, allocatable, dimension(:,:) :: offsets
  real, allocatable, dimension(:,:,:) :: optres, apres, newxypos

  double precision, allocatable, dimension(:,:) :: stars, times
    
  character (len=50) :: starpos, offsetfile, timesfile

  character (len=15), allocatable, dimension(:) :: starnames
  character (len=50), allocatable, dimension(:) :: datafiles, timefiles
  character (len=1), allocatable, dimension(:,:,:) :: flagres

! Varialbes used for sky estimation
  real :: bad_sky_chi, bad_sky_skw, skynos, skycnt, skyerr
  
  character :: cflag
  character(len=10) :: bskyanswer 

  logical :: sky_stuffed

! Variables for aperature photometry
  integer :: ibox

  real :: aprad, adu
  real :: apx, apy
  real :: apflux, aperror

  character(len=10) :: apanswer

  logical :: aperature

! Variables for optimal photometry
  integer :: iopt  ! The ID# of the star whose values are to be optimized

  character(len=10) :: optanswer 

  logical :: optimal

! ***************************************************************
! Actual program starts here.
! ***************************************************************

! Initialize variables, set global constants
  npsf = 0         ! Number of PSF stars
  nstar = 0        ! Number of Target stars
  status = 0       ! Required by CFITSIO
  optimal=.true.   ! Switch for running the optimal photometry code
  aperature=.true. ! Switch for running the aperature photometry code

! Read input parameters
  read*, filename, flagfile, npsf, nstar, verbosein
  read*, bad_sky_skw, bad_sky_chi, fwhm, clip_fwhm, aprad, iopt, searchrad, adu

! Set the program output type
  verbose = .False.
  if (verbosein=='y'.or.verbosein=='Y') verbose=.True.

! ***************************************************************
! Read the co-ordinates for each PSF star.
! ***************************************************************
  if (npsf>0) then
     ! Allocate PSF star arrays
     allocate(psfnames(npsf))
     allocate(psfstars(3,npsf))           !star id, x, y pixel coordinates
     allocate(xpos0(npsf))
     allocate(ypos0(npsf))
     allocate(dpsf(npsf))

     do i=1, npsf
        psfstars(1,i)=i
        read*, psfnames(i), psfstars(2,i), psfstars(3,i)
     end do
     
     if (verbose) then
        print*, '@ Read ', npsf, 'PSF star location(s): '
        do i=1, npsf
           print*, "@ PSF Star ",psfnames(i)," at position (X,Y):", &
                   real(psfstars(2,i)), real(psfstars(3,i))    
        end do
        print*
     end if 
  
  else
     print*, '@ No PSF stars. Proceeding with Aperature Photometry only!'   
     print*
     optimal=.False.
     clip_fwhm = -1.0
  end if

! *******************************************************************
! Read the co-ordinates for each target star.
! *******************************************************************

  if (nstar>0) then
     ! Allocate target star arrays
     ! Since this is the pipeline version it runs every time for a single frame.
     ! Therefore ntimes must equal 1.
     ntimes = 1
     nstar=nstar+npsf
     allocate(starnames(nstar))
     allocate(dpos(nstar))
     allocate(stars(3,nstar))             
     allocate(optres(nstar,2,ntimes))   !star id, flux, eflux, no of points
     allocate(apres(nstar,2,ntimes))    !star id, flux, eflux, no of points
     allocate(flagres(nstar,2,ntimes))  !star id, optimal flag, aperature flag, no of points
     allocate(newxypos(nstar,2,ntimes)) !star id, xpos, ypos, no of points
     allocate(newpsfpos(npsf,2,ntimes)) !star id, xpos, ypos, no of points
     allocate(seeing(ntimes))           !seeing

     do i=1, nstar
        stars(1,i)=i
        read*, starnames(i), stars(2,i), stars(3,i)
     end do
     
     if (verbose) then
        print*, '@ Read ', nstar, 'star location(s): '
        do i=1, nstar
           print*, "@ Object ",starnames(i)," at position (X,Y):", &
                   real(stars(2,i)), real(stars(3,i))    
        end do
        print*
     end if 
     
  else
     print*, 'WARNING: Target Star Positions Not Found. Nothing to do.'
     print*, 'Exiting....!'   
     print*
     stop
  end if

! ***************************************************************
! Setup all the optimal and aperature photometry options
! ***************************************************************

  ! If not set, set the clipping radius for the PSF.
  ! Most of the signal-to-noise is obtained by setting the clipping
  ! radius to be 2*fwhm.  But there is little gain in speed by
  ! making it less than 3 pixels.
  if (clip_fwhm <= 0.0) then 
     optimal = .False.
     clip_fwhm=-1.0
  else
     if (clip_fwhm == 0.0 .or. clip_fwhm < fwhm) clip_fwhm=max(3.0/fwhm, 2.0)
  end if

  ! If not set, set the default aperture radius to the optimum (as shown in the paper),
  ! but don't let it be less than 2 pixels.
  if (aprad<0.0) then
     aperature = .False.
     aprad=-1.0
  else
     if (aprad == 0.0 .or. aprad < fwhm) aprad=max(2.0, 2.0*fwhm/3.0)
  end if

  ! Set the stage for the photometry to be optimized for a particular star or
  ! for a sky limited case.
  if (optimal) then
     optstar: do
              if (iopt > 0) then
                 do i=1, nstar
                    if (int(stars(1,i)) == iopt) then
                       iopt=i
                       exit optstar
                    end if
      	         end do
	         print*, "WARNING: The star to optimize for is not in the target list."
                 print*, "         (Optimal) Will proceed with sky limited optimization."
                 print*
                 iopt=-1
                 exit optstar
              else
                 iopt=-1
                 exit optstar
	      end if   
     end do optstar
  end if

! Set the search radius. If the radius is negative then the positions will be 
! considered as fixed. If not supplied or if the supplied one is less than 2*FWHM 
! then search radius is set to 2*FWHM
  dpsf=searchrad
  dpos=searchrad
  if (searchrad < 2.0*fwhm) then
     dpsf=2.0*fwhm
     dpos=2.0*fwhm
  end if
  if (searchrad <=0.0) then
     print*, "WARNING: Search radius is zero or negative. No centroiding will be performed"
     print*, "         and the object positions will be fixed to their input values."
     print*
     dpos = -1.0
  end if

! Summary of all the setup options.
  if (verbose) then
     print*, 'Photometry will be performed with the following options:'
     print*, '========================================================'
     print*, 'Number of PSF stars              :',npsf
     print*, 'Number of target stars           :',nstar
     print*, 'Sky skew flag limit              :',bad_sky_skw
     print*, 'Sky Chi^2 flat limit             :',bad_sky_chi
     print*, 'FWHM of image                    :',fwhm
     print*, 'PSF clipping radius              :',clip_fwhm
     print*, 'Aperature photometry radius      :',aprad
     print*, 'Star to be optimized (-ve if all):',iopt
     print*, 'Centroid search radius           :',dpsf(1)
     print*, 'Detector gain (e-/ADU)           :',adu
     print*, '======================================================='
     print*
  end if

! Check to see what kind of photometry will be performed.
  if (.not.optimal) then
     print*, '@ Optical Photometry will not be performed'
     flagres(:,1,:)='-'
  end if
  if (.not.aperature) then
     print*, '@ Aperature Photometry will not be performed'
     flagres(:,2,:)='-'
  end if

! If both the optimal and photometry options are negative exit the program.
  if ((.not.optimal).and.(.not.aperature)) then
     print*, '@ Nothing to do! Exiting...'
     stop
  end if

! ***************************************************************
! Begin the Photometry Processes
! ***************************************************************

  frame: do fnum=1, ntimes
  
         ! First read in the image data into a data array.
         ! Initialize CFITSIO required variables. 
         if (verbose) print*, "@ Filename: ", filename         
         nfound=0
         naxis=2
         group=1
         inc=1
         nullval=-999
         status=0
         readwrite=0
         call ftgiou(unit, status)
         call ftopen(unit, filename, readwrite, blocksize, status)
         call ftgknj(unit, 'NAXIS', 1, 2, naxes, nfound, status)
         if (nfound /= 2)then
            if (verbose) print *,'@ Could not read image. Moving to the next frame...'
            cycle frame
         end if
         ! Allocate the data array and determine the first and last pixels of the image.
         fpix=1
         lpix=naxes
	 low=fpix
         high=lpix
         allocate(data(1:naxes(1),1:naxes(2)))
         allocate(numflag(1:naxes(1),1:naxes(2)))
         allocate(pix_flg(1:naxes(1),1:naxes(2)))
         ! Read the data into the data array
         call ftgsve(unit, group, naxis, naxes, fpix, lpix, inc, nullval, data, anynul, status )
         print*, "Status 2: ", status
!         Print somethings from the FITS file to confirm that it read properly.
!         print*, naxis, naxes(1), naxes(2)
!         print*, fpix(1), fpix(2)
!         print*, lpix(1), lpix(2)
!         print*, "array value at x=100, y=100: ", data(100,100)
!         print*, "maximum value position: ", maxloc(data)
! 	      print*, "*********************"
         call ftclos(unit, status)
         call ftfiou(unit, status)
         ! Check for any error, and if so print out error messages.
         ! The PRINTERROR subroutine is listed at the end of this file.
         if (status .gt. 0 .and. verbose) call printerror(status)

!        Now do the same and read in the file containing the numeric pixel flags.
!        re-initilized cfitsio variables.
         nfound=0
         naxis=2
         group=1
         inc=1
         nullval=-999
         status=0
         readwrite=0

!        Open the flag file and read in the flags.
         call ftgiou(unit, status)
         call ftopen(unit, flagfile, readwrite, blocksize, status)
         call ftgknj(unit, 'NAXIS', 1, 2, naxes, nfound, status)
!        if there is a problem set all flags to OK.
         if (nfound /= 2)then
            if (verbose) print *,'@ Flag file cannot be read. All flags will be set to O'
            pix_flg = 'O'
         end if

!        Read in the numeric pixel flags
         call ftgsve(unit, group, naxis, naxes, fpix, lpix, inc, nullval, numflag, anynul, status )
         call ftclos(unit, status)
         call ftfiou(unit, status)
!        Check for any error, and if so print out error messages.
         if (status .gt. 0 .and. verbose) call printerror(status)

!        Print somethings from the pixel flag file to confirm that it read properly.
!         print*, "Pixel flag at x=100, y=100: ", numflag(100,100)
!         print*, "maximum value position: ", maxloc(numflag)
! 	     print*, "*********************"

!        Convert pixel flag numbers into opphot's scheme
         forall (i=1:naxes(1),j=1:naxes(2),numflag(i,j)==0) pix_flg(i,j)="O"
         forall (i=1:naxes(1),j=1:naxes(2),numflag(i,j)==1) pix_flg(i,j)="F"
         forall (i=1:naxes(1),j=1:naxes(2),numflag(i,j)==2) pix_flg(i,j)="S"
         forall (i=1:naxes(1),j=1:naxes(2),numflag(i,j)==3) pix_flg(i,j)="L"
!         forall (i=1:naxes(1),j=1:naxes(2),numflag(i,j)==4) pix_flg(i,j)="C"

!        Estimate the position of the PSF stars
         if (clip_fwhm > 0.0) then
            xpos0 = real(psfstars(2,:))
            ypos0 = real(psfstars(3,:))

            do i=1, npsf              
               if (verbose) print*, '@ Estimated position of PSF star ',i, 'is: ', &
                            xpos0(i), ypos0(i)
            end do

            if (verbose) print*, "@ Starting FWHM and Clipping radius is: ", fwhm, clip_fwhm

!           Fit the PSF stars
            call psf_calc(data, pix_flg, npsf, xpos0, ypos0, dpsf, &
                          adu, high, low, fwhm, shape_par, ipsf, nfit, verbose)

            if (ipsf == -1) then
               optimal = .False.
               if (verbose) then 
                  print*, 'WARNING: No good PSF star could be found within frame.'
                  print*, '         Setting optimal fluxes to zero for this frame.'
                  print*, '         ...continuing with aperature photometry only!'
               end if
               flagres(:,1,fnum)='J'
       	       optres(:,1,fnum)=0.0
	       optres(:,2,fnum)=0.0
            end if
	    
            if (optimal) then
               ! We now have a better estimate of fwhm, so use this instead.
               fwhm=1.665*sqrt(shape_par(1)*shape_par(2))
               seeing(fnum)=fwhm ! Save the estimate of the seeing to an array.
               cliprad = clip_fwhm*1.665*sqrt(shape_par(1)*shape_par(2)) 
               !cliprad = 3.0*fwhm            
               
    	       if (verbose) then
                  print*, "@ New FWHM and Clipping radius is: ", fwhm, cliprad
                  print*, '@ Fitted PSF star ', int(psfstars(1,ipsf)), &
                          ' which gave FWsHM ', 1.665*shape_par(1), 1.665*shape_par(2)
                  print*, '@ Rotated at an angle of ', 57.29*shape_par(3), &
                          ' degrees from the vertical.'
                  print*, '@ Chosen using ', nfit, ' stars.'
                  print*, '@ Will use a clipping radius of ', cliprad, ' pixels.'
               end if

               !Now get a handle on the flux in the star to be extracted optimally.
               if (iopt > 0) then
                  xpos = real(stars(2,iopt))
                  ypos = real(stars(3,iopt))
                  ! Call the extraction, with the normalisation explicitly set to one, and
                  ! return the peak flux in optnrm.
                  call extr(data, pix_flg, xpos, ypos, dpos(1), adu, high, low, fwhm, &
                            cliprad, shape_par, 0.0, .false., xcomp, ycomp, optflux, opterror, &
                            xfit, yfit, xerr, yerr, optnrm, cflag, skynos, verbose)

                  if (cflag /= 'O') then
                     optimal = .False.
                     if (verbose) then 
                        print*, 'WARNING: Star to be optimised has flag ', cflag
                        print*, '         will set its optimal flux to zero and continue'
                        print*, '         with aperature photometry!'
                     end if
                     flagres(:,1,fnum)=cflag
       	             optres(:,1,fnum)=0.0
	                 optres(:,2,fnum)=0.0
                     optnrm=0.0
                  end if

                  if (verbose) print*, '@ Extractions will be optimised for stars with ',&
                                     optnrm,' peak counts.'
                  optnrm=optnrm/(skynos*skynos)
               else                            
                  optnrm=0.0
               end if

            end if

            if (verbose) then
               print*
               print*
            end if

         end if

         ! Don't do optimal photometry, but we will need a FWHM for fitting
         ! the star to find its position.
         if (.not.optimal) then
            shape_par(1)=fwhm/1.665
            shape_par(2)=fwhm/1.665
            shape_par(3)=0.0
            ! And a cliprad.
            cliprad=2.0*fwhm
            ! And (believe it or not) a star flux to be optimised.
            optnrm = 0.0
         end if

         ! ***********************************************************
         ! Go around this loop for each star.
         ! ***********************************************************

         each_star: do istar=1, nstar

	        comp=.False.

            ! Calculate approximate position of star.
            xpos = real(stars(2,istar))
            ypos = real(stars(3,istar))

            if (optimal) then

               ! Call the optimal extraction routine.
               call extr(data, pix_flg, xpos, ypos, dpos(1), adu, & 
                         high, low, fwhm, cliprad, shape_par, optnrm, &
                         comp, xcomp, ycomp, optflux, opterror, &
                         xfit, yfit, xerr, yerr, peak, cflag, skynos,verbose)
		       if (verbose) then
		          print*, '@ New position of object ',trim(starnames(istar)),':', xfit, yfit
                          print*, '@ Opphot flux for object ',trim(starnames(istar)),&
                          ' is: ',optflux,'+/-',opterror, "(counts)"
		       end if
               ! Put the flux and its error in the output array and flag the result.	    
       	       optres(istar,1,fnum)=optflux
	       optres(istar,2,fnum)=opterror
               flagres(istar,1,fnum)= cflag
               ! Set the newly calculated positions.
               newxypos(istar,1,fnum)=xfit
               newxypos(istar,2,fnum)=yfit

               ! If the star was not in the frame go to the next one.
               if (opterror<0.0) then 
		  newxypos(istar,1,fnum)=-1.0
		  newxypos(istar,2,fnum)=-1.0
	          apres(istar,1,fnum)=0.0 
		  apres(istar,2,fnum)=-1.0
                  cycle each_star
               end if

	    end if

	    ! In addition to Optimal photometry do aperature photometry as well
	    ! using the centroided positions that we got from optimal photometry
	    ! xfit, yfit.
            if (aperature) then

               sky_stuffed = .false.

               ! Set the star position if not done by optimal fitting.
               if (.not.optimal) then
                  xfit = xpos
                  yfit = ypos
               end if

               !The box size is determined in the way described in the paper.
               ibox=int(sqrt(628.4*aprad*aprad+ 4.0*fwhm*fwhm))
               call skyfit(data, xfit, yfit, fwhm, ibox, low, high, &
                           pix_flg, skycnt, skyerr, skynos, cflag)

               ! Flag the data accordingly.
               ! Check to see if the sky background had problems
               if (cflag == 'I') sky_stuffed = .true.

               ! If the sky background failed set fluxes to zero. Otherwise carry
               ! out the aperture photometry.
               if (cflag == 'B') then
                  apres(istar,1,fnum)=0.0
                  apres(istar,2,fnum)=0.0
                  flagres(istar,2,fnum)=cflag
               else
                  apflux=sum_aper(data, aprad, skycnt, skyerr, skynos, &
                              xfit, yfit, adu, low, high, pix_flg, &
                              aperror, cflag)

                  if (sky_stuffed) cflag='I'
                  flagres(istar,2,fnum)=cflag
               end if

               ! Screen output for aperture photometry results.
               if (verbose) then

                   print*, '@ Apphot flux for object ',trim(starnames(istar)), &
                           ' is: ',apflux, '+/-',aperror, "(counts)"
                   print*, '@ Object is flagged as: ', cflag               
               end if

      	       apres(istar,1,fnum)=apflux
	       apres(istar,2,fnum)=aperror
               flagres(istar,2,fnum)=cflag

            else
            
      	       apres(istar,1,fnum)=-1.0
	       apres(istar,2,fnum)=-1.0
               flagres(istar,2,fnum)= "N/A"
               
            end if

         if (verbose) print*

         end do each_star

         deallocate(data)
         deallocate(pix_flg)

  end do frame

! The following print statement is used for run_rtphos.py to find the photometry
! results. It has 10 "*" so if this is changed it also needs to change the 
! appropriate line in run_rtphos.py   
  print*, "**********"

! ***************************************************************
! Write the output files.
! ***************************************************************
!  wformat='(I5," ",F13.4," ",F10.4," ",F7.4)'
! Write the optimal photometry output to screen.
! if (optimal) then
     do i=1, nstar
        do j=1, ntimes
           print*, starnames(i), optres(i,1,j), optres(i,2,j), seeing(j), newxypos(i,1,j), newxypos(i,2,j), flagres(i,1,j) 
        end do   
     end do
!  end if

! Write the aperature photometry output to screen.
!  if (aperature) then
     do i=1, nstar
        do j=1, ntimes
           print*, starnames(i), apres(i,1,j), apres(i,2,j), seeing(j), newxypos(i,1,j), newxypos(i,2,j), flagres(i,2,j) 
        end do   
     end do
!  end if

end program optimal_phot

!******************************************************************************
! printerror - Prints an error if the input FITS file is not properly read.
!******************************************************************************

subroutine printerror(status)

! This subroutine prints out the descriptive text corresponding to the
! error status value and prints out the contents of the internal
! error message stack generated by FITSIO whenever an error occurs.

! Variable declarations
! ---------------------
  integer :: status

  character (len=30) :: errtext
  character (len=80) :: errmessage*80
! ---------------------

! Check if status is OK (no error); if so, simply return
  if (status <= 0)return

! The FTGERR subroutine returns a descriptive 30-character text string that
! corresponds to the integer error status number.  A complete list of all
! the error numbers can be found in the back of the FITSIO User's Guide.
  call ftgerr(status, errtext)
  print *,'* FITSIO Error Status =',status,': ',errtext

! FITSIO usually generates an internal stack of error messages whenever
! an error occurs.  These messages provide much more information on the
! cause of the problem than can be provided by the single integer error
! status value.  The FTGMSG subroutine retrieves the oldest message from
! the stack and shifts any remaining messages on the stack down one
! position.  FTGMSG is called repeatedly until a blank message is
! returned, which indicates that the stack is empty.  Each error message
! may be up to 80 characters in length.  Another subroutine, called
! FTCMSG, is available to simply clear the whole error message stack in
! cases where one is not interested in the contents.
  call ftgmsg(errmessage)
  do while (errmessage /= ' ')
     print *,errmessage
     call ftgmsg(errmessage)
  end do

end subroutine printerror
