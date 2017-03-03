   module opt_extr
 
     ! This is the module that carries out the optimal extraction.

     ! Version 3.0

     ! There are two subroutines designed to be used by the programmer creating
     ! a programme to carry out optimal extraction.  psf_calc will calculate
     ! the model PSF, and opt_calc will carry out the optimal photometry.

     ! By Tim Naylor

     ! Modifications.
     ! Version 1.1 
     !   Has the subroutines skyhst and skymod.
     ! Version 1.2 
     !   Allows flux_aper to return the peak counts in the aperture.
     ! Version 2.0
     !   Allows the position angle of the Gaussian to run free.
     ! Version 2.1
     !   Bug cured whereby position angle was running free even after
     !   the shape parameters had been fixed.
     ! Version 2.2
     !   More comments to aid porting added.
     ! Version 2.3
     !   Returns an estimate of the error in the position.
     ! Version 2.4
     !   Negative glitches foul up the noise model as data-sky becomes
     !   negative, which can make skynos*skynos + (data-sky)/adu negative,
     !   and hence ask for the square root of a naegative number.  Obviously
     !   photometry around such a pixel is junk, but as we don't want the
     !   program to crash, it will set the noise to be the sky noise.
     ! Version 2.5
     !   Skyfit now spots if there are no good pixels available in the sky box.
     ! Version 2.6
     !   Minor Changes so that it compiles with the NAG Fortran 95
     !   complier, as well as the DEC Fortran 90 compiler.
     ! Version 2.7
     !   Bug cured whereby the rotation of the Gaussian fitted to the PSF
     !   star was limited to 0 -> 45 degrees, when it should be -45 -> +45.
     ! Version 2.8
     !   Minor change to opt_extr to speed up calculation of var for sky
     !     limited case.
     !   Cured a bug which crashed the program if the companion star was off
     !     the frame edge.
     !   Checks for high skewness in the sky histogram.
     !   Uses iterative sigma clipping to make initial sky and sigma estimates
     !     in skyfit.
     !   Lee Howells' fixes for a couple of bugs added.
     !   GFIT arrays normalised before going into curfit to over numerical
     !     overflows.
     ! Version 3.0
     !   General tidy-up of code, including changing many subroutine arguments.
     !   Introduction of pixel flags.
     !   Flags changed to characters rather than integers.
     ! Version 3.1
     !   Two bugs found that only realy show up when fitting a negative sky.
     !   The derivative for the position of the peak of the skewed Gaussian 
     !   was not calculated properly when the sky was negative, and the 
     !   clipped mean wasn't done properly either.
      
     implicit none

     ! The fitting routines are essentially 1D, but need to know the size of
     ! the second dimension of the image they are fitting.
     integer, private :: naxis2

     ! A local debugging variable, set true if you want lots of messages.
     logical, private :: debug=.True.

     ! The sky histogram, and the fit to it.  Placed here so they can be
     ! accessed by other programs through the routines skyhst and skymod.
     integer, private, parameter :: maxbin=100
     real, private, dimension(maxbin):: xhst, yhst
     ! And the parameters of the last skyfit.
     real :: ymax
     real, dimension(4) :: sky_par
     ! Values from the last sky determination.
     real :: sky_cnt, sky_err, sky_chi
     real, private, save :: bad_sky_chi=-1.0, bad_sky_skw=-1.0

     contains

       subroutine psf_calc(data, pix_flg, npsf, xpsf, ypsf, dpsf, &
       adu, high, low, fwhm, shape_par, ipsf, nfit, verbose)

       ! This subroutine calculates the point spread function.

       ! The idea is that this subroutine is given a list of PSF stars.  It
       ! will fit the first 49 it finds that are on the frame and not pixel
       ! flagged.  It will then find the star with the median FWHM, and
       ! return the parameters of the fit to that star.

       ! Inputs.

       ! The data array (can be adjusted if the fit is subtracted by gfit).
       real, dimension(:,:), intent(inout) :: data
       ! The pixel flag array.
       character, dimension(:,:), intent(in) :: pix_flg

       ! The number of available PSF stars.
       integer, intent(in) :: npsf
       ! The positions of the psf stars.
       real, dimension(:), intent(in) :: xpsf, ypsf
       ! The radius to be searched for the psf star.  If this is zero it
       ! will take the position as fixed.
       real, dimension(:), intent(in) :: dpsf

       ! The number of detected photons per count.
       real, intent(in) :: adu
       ! The highest and lowest pixel numbers to be used.
       integer, dimension(2), intent(in) :: high, low

       ! The approximate seeing, in pixels.
       real, intent(in) :: fwhm

       ! Option to print out more calculation details.
       logical :: verbose
               
       ! Outputs.

       ! A flag, set to the position in the list of star used, or -1 if no PSF 
       ! star could be found.
       integer, intent(out) :: ipsf

       ! The three parameters of the PSF.
       real, dimension(3), intent(out) :: shape_par

       ! The number of stars used when deciding which star to fit.
       integer, intent(out) :: nfit


       ! Locals.
       ! Two counters.
       integer :: i
       ! The size of the sky box.
       integer :: ibox
       ! The sky estimate, the error on the sky estimate, and the RMS of the sky.
       real :: skycnt, skyerr, skynos
       ! A flag returned by called subroutines.
       character :: cflag
       ! The parameters of the model.
       real, allocatable, dimension(:,:) :: a_par
       ! The (unused) error in position.
       real, dimension(2) :: e_pos
       ! The number of stars to be fitted.
       integer, parameter :: mfit=49
       ! An array of star ID numbers.
       integer, dimension(mfit) :: idstar

       ! Allocate the parameter array.
       allocate(a_par(6,mfit))

       ! Set it up so that if no psf star is found, this is signalled.
       ipsf = -1
       nfit=1
       psfstar: do i=1, npsf
         if (verbose) print*, '@ Trying star ', i, ' in the PSF list at ', xpsf(i), ypsf(i)
         ! Do an initial check to see if the star is likely to be in frame.
         if (nint(xpsf(i)-fwhm)<low(1) .or. nint(xpsf(i)+fwhm)>high(1) .or. &
             nint(ypsf(i)-fwhm)<low(2) .or. nint(ypsf(i)+fwhm)>high(2)) then
           ! The PSF star is too close to the frame edge.
	     if (verbose) then
               print*, 'WARNING: PSF star',i,'is too close or outside the CCD frame'
               print*, '         it will be ignored this time.'
	     end if
           cycle psfstar
         end if

         ! The sky box size is determined in the way described in the paper
         ibox=int(sqrt(453.2*fwhm*fwhm+ 4.0*fwhm*fwhm))
         call skyfit(data, xpsf(i), ypsf(i), fwhm, ibox, low, high, &
         pix_flg, skycnt, skyerr, skynos, cflag)
!         print*, '@ skyfit gave a sky of ', skycnt,' and a noise of ', skynos
         if (cflag /= 'O') then
           ! The the sky fit failed.  Can't think why this may be, but go to 
           ! next psf star anyway.
           ! print*, 'Skyfit failed with flag ', cflag
           cycle psfstar
         end if

         ! Set up for fitting the PSF star.
         ! Set position.
         a_par(5,nfit)=xpsf(i)
         a_par(6,nfit)=ypsf(i)
         ! Set Gaussian width.
         a_par(1,nfit)=fwhm/1.665
         a_par(2,nfit)=fwhm/1.665
         ! And the rotation (in radians).
         a_par(3,nfit)=0.0
         ! Set normalisation to 1, and hunt for best counts.
         a_par(4,nfit)=1.0
         if (debug) print*, '@ Approximate position ', a_par(5,nfit), &
         a_par(6,nfit)
         ! And set normalisation.          
         a_par(4,nfit)=max(data(nint(a_par(5,nfit)),&
         nint(a_par(6,nfit)))-skycnt, skynos) 
         call gfit(data, .false., dpsf(i), .false., skycnt, skynos, &
         adu, low, high, pix_flg, a_par(:,nfit), e_pos, cflag)
         if (cflag /= 'O') then
           if (debug) print*, '@ Fit failed because of flag ', cflag
           cycle psfstar
         end if
         ! Set the fourth parameter to be the peak signal-to-noise.
         a_par(4,nfit)=a_par(4,nfit)/skynos
         ! If you've got to here its a good star.
         idstar(nfit)=i
         nfit=nfit+1
         if (nfit-1 == mfit) exit psfstar
       end do psfstar
       nfit=nfit-1

       if (debug) print*, '@ Sucessfully fitted ', nfit, ' PSF stars.'

       if (nfit > 0) then
         ! Check that there are some stars with peak s:n > 10.
         if (maxval(a_par(4,1:nfit)) < 10.0 .or. nfit==1) then
           if (debug) print*, '@ No stars with S:N > 10.'
           ! Better take the first good star then.
           i=1
           nfit=1
         else
           ! Throw out all the low signal-to-noise stuff.
           i=1
           cleanup: do
             if (a_par(4,i) < 10.0) then
               if (i /= nfit) then
                 a_par(:,i:nfit-1)=a_par(:,i+1:nfit)
                 idstar(i:nfit-1)=idstar(i+1:nfit)
               end if
               nfit=nfit-1
             else
               i=i+1
             end if
             if (i-1 >= nfit) exit cleanup
           end do cleanup
           if (debug) print*,'@ After cleaning out poor S:N there are ', &
           nfit, 'stars.'
           i=imedian(a_par(1,1:nfit)*a_par(2,1:nfit), nfit)
         end if
         shape_par=a_par(1:3,i)
         ipsf=idstar(i)
         !write(49,*) nfit
         !do i=1, nfit
         !  write(49,*) i, 1.665*sqrt(a_par(1,i)*a_par(2,i))
         !end do
       end if

       ! deallocate the array...
       deallocate(a_par)

       if (debug) print*, '@ Done psf_calc.'

       end subroutine psf_calc


       integer function imedian(indata, ndata)

       real, dimension(:), intent(in) :: indata
       integer, intent(in) :: ndata
       
       integer :: k, l, m, locval
       real :: aval
       real, dimension(ndata) :: data
       integer, dimension(ndata) :: locator

        do k=1, ndata
          locator(k)=k
        end do
        data=indata(1:ndata)

        do 140 k=2, ndata
          aval=data(k)
          locval=locator(k)
          do l=1, k-1
            if (aval < data(l)) then
              do m=k, l+1, -1
                data(m)=data(m-1)
                locator(m)=locator(m-1)
              end do
              data(l)=aval
              locator(l)=locval
              goto 140
            end if
          end do
140     continue

        imedian=locator(ndata/2 + 1)

       end function imedian
       

       subroutine extr(data, pix_flg, xpos, ypos, dpos, &
       adu, high, low, fwhm, cliprad, shape_par, &
       optnrm, companion, xcomp, ycomp, flux, error, &
       xfit, yfit, xerr, yerr, peak, cflag, skynos_r, verbose)

       ! This subroutine calculates the flux optimally.     
                     
       ! Inputs.

       ! The data array, and the highest and lowest pixel numbers to be used.
       ! (data can be adjusted if the fit is subtracted by gfit.)
       real, dimension(:,:), intent(inout) :: data
       ! The pixel flag array.
       character, dimension(:,:), intent(in) :: pix_flg
       ! The position of the star.
       real, intent(in) :: xpos, ypos
       ! The radius to be searched for the star.  If this is zero it
       ! will take the position as fixed.
       real, intent(in) :: dpos
       ! The number of detected photons per count.
       real, intent(in) :: adu
       integer, dimension(2), intent(in) :: high, low
       ! The approximate seeing, in pixels.
       real, intent(in) :: fwhm
       ! The radius to clip the mask.
       real, intent(in) :: cliprad
       ! The three parameters of the PSF.
       real, dimension(3), intent(in) :: shape_par
       ! The peak flux the photometry is normalised for, divided by the sky
       ! noise squared.  Zero is sky limit.
       real, intent(in) :: optnrm
       ! Is there a companion star?
       logical, intent(in) :: companion
       ! What's its position?
       real, intent(in) :: xcomp, ycomp

       ! The RMS of the sky.
       real, optional, intent(out) :: skynos_r
       
       ! The returned values
       ! The flux and its error.
       real, intent(out) :: flux, error
       ! The fitted position of the star (set to xpos and ypos if position is
       ! fixed, or fit fails).
       real, intent(out) :: xfit, yfit
       ! And their error (set to zero if fit fails).
       real, intent(out) :: xerr, yerr
       ! The peak flux from the fitted Gaussian, or, if dpos is zero, the
       ! highest pixel.
       real, intent(out) :: peak
       ! A flag returned.  
       !     O if everything O.K.,  
       !     E if star too close to frame edge of the frame.
       !     B if the fit to the sky failed (inherited from skyfit).
       !     I if the sky is ill-determined (inherited from skyfit).
       !     Or the value of any pixel flags.
       character, intent(out) :: cflag

       ! Locals.

       ! The size of the sky box.
       integer :: ibox
       ! The sky estimate, the error on the sky estimate and the RMS of the sky.
       real :: skycnt, skyerr, skynos
       ! The fitted parameters.
       real, allocatable, dimension(:) :: a_par
       ! A counter.
       integer :: i
       ! The error in star position.
       real, dimension(2) :: e_pos
       ! The flag from the sky fitting
       character :: sky_flag

       ! Option to print out more calculation details.
       logical :: verbose

  !      debug=.True.

       if (debug) print*, 'Entering extr for star at ', xpos, ypos

       cflag='O'
       ! Set the flux and its error to zero, so if the fit fails the values
       ! are fairly meaningless.
       flux = 0.0
       error = 0.0
       ! The same for the fitted position errors.
       xerr = 0.0
       yerr = 0.0
       ! And preserve the position if no fit is made.
       xfit = xpos
       yfit = ypos

       ! Do an initial check to see if the star is likely to be in frame.
       if (nint(xpos)<low(1) .or. nint(xpos)>high(1) .or. &
           nint(ypos)<low(2) .or. nint(ypos)>high(2)) then
          if (verbose) then
             print*, '*** WARNING: Star is too close or outside the CCD frame'
             print*, '             it will be ignored this time.'
             print*, '             Setting fluxes to zero and putting negative errors'
	  end if
	   flux=0.0
           error=-1.0
           cflag = 'E'
           return
       else
         ! Measure the sky. The box size is determined in
         ! the way described in the paper
         ! Measure it first in a smaller box.
         ibox=int(sqrt(453.2*fwhm*fwhm+ 4.0*fwhm*fwhm))
         call skyfit(data, xpos, ypos, fwhm, ibox, low, high, pix_flg, &
         skycnt, skyerr, skynos, sky_flag)
         if (present(skynos_r)) skynos_r=skynos
         if (debug) print*, 'Returned from skyfit with flag ', sky_flag
                               
         ! Set up the parameters for the fit.
         if (companion) then
           allocate(a_par(9))
           ! Set the normalisation rather crudely, especially crude if the
           ! companion star is off the frame edge.
           a_par(7)=data( max(low(1),min(nint(xcomp),high(1))), &
                          max(low(2),min(nint(ycomp),high(2))) )-skycnt
           a_par(8)=xcomp
           a_par(9)=ycomp
         else
           allocate(a_par(6))
         end if
         do i=1, 3
           a_par(i)=shape_par(i)
         end do
         ! Set the position.
         a_par(5)=xpos
         a_par(6)=ypos
         if (dpos > 0.0) then
           ! If the position is free set normalisation to 1, and hunt for 
           ! best counts.
           a_par(4)=1.0
           ! Set the normalisation to the peak counts.
           a_par(4)=max(data(nint(a_par(5)),nint(a_par(6)))-skycnt, skynos)
           ! Reset the position.
           a_par(5)=xpos
           a_par(6)=ypos
         else
           a_par(4)=max(data(nint(a_par(5)),nint(a_par(6)))-skycnt, skynos)
         end if
         if (dpos > 0.0) then
           call gfit(data, .true., 0.0,  .false., skycnt, skynos, adu, &
           low, high, pix_flg, a_par, e_pos, cflag)
           call gfit(data, .true., dpos, .false., skycnt, skynos, adu, &
           low, high, pix_flg, a_par, e_pos, cflag)
           ! cflag is deliberately ignored.  We will check for saturation when
           ! we do the photometry.
           ! Update the positions.
           xfit=a_par(5)
           yfit=a_par(6)
           xerr=e_pos(1)
           yerr=e_pos(2)
         end if
         peak=a_par(4)

         ! Finally, do the optimal photometry.
         flux=sum_flux(data, pix_flg, skycnt, skyerr, skynos, &
         a_par, cliprad, optnrm, adu, low, high, error, cflag)
       
         sky_cnt=skycnt
         sky_err=skyerr

         !if (sqrt((xfit-xpos)**2.0 + (yfit-ypos)**2.0) > fwhm) then
         !  ! Probably something spurious on the wings of a star.
         !  cflag='D'
         !end if

         if (cflag == 'O') cflag=sky_flag

         deallocate (a_par)

       end if

       if (debug) print*, 'Exiting from s/r extr.'

       end subroutine extr
      
       subroutine lastsky(skycnt, skyerr)
      
         real, intent(out) :: skycnt, skyerr
       
         skycnt=sky_cnt
         skyerr=sky_err
       
       end subroutine lastsky


       subroutine opt_debug_on()

       ! Sets the debug flag, which will print lots of messages (and I
       ! mean lots...)

       debug=.False.

       end subroutine opt_debug_on



       ! ** Here endeth the "user" routines, and begineth the sky fitting 
       ! ** routines.


       subroutine skyfit(data, xpos, ypos, seeing, ibox, low, high, &
       pix_flg, skycnt, skyerr, skynos, cflag)

       ! This subroutine fits a skewed Gaussian to the pixel distribution
       ! in an area of sky.

       ! Inputs
       real, dimension(:,:), intent(in) :: data
       ! The first guess position of the star.
       real, intent(in) :: xpos, ypos
       ! The FWHM seeing.
       real, intent(in) :: seeing
       ! The length of the side of the sky box.
       integer, intent(in) :: ibox
       ! The lowest and highest usable pixel numbers.
       integer, intent(in) :: low(2), high(2)
       ! The pixel flag array.
       character, dimension(:,:), intent(in) :: pix_flg
       
       ! Outputs
       ! The estimated sky counts, error, and deviation.
       real, intent(out) :: skycnt, skyerr, skynos
       ! A return flag.  
       !      B means failed to fit sky.
       !      I If fit is unreliable because of skewness of sky.
       character, intent(inout) :: cflag
       
       ! Loop variables.
       integer :: icount, i, j
       ! The maximum size of the sky box allowing for edges.
       ! integer :: jbox
       ! The area to be used to estimate the sky.
       integer :: iyslo, iyshi, ixslo, ixshi
       ! The mean in the sky box, the deviation
       real :: skymen, skydev
       
       ! For the histogram.
       integer :: nbin
       real :: hstep, hmin
       integer, parameter :: intdat=1
       logical :: first
       integer :: iskymod
       real :: total
       !real :: total, b_chi
       integer :: n_reject

       integer :: jflag, iflag

       ! For finding the sky level.
       real, allocatable, dimension(:) :: clip_array

       ! Set default values for the output, to ensure they are defined for
       ! compilers which don't set things to zero.  I've chosen 1 for skyerr
       ! and skycnt as these are errors, and the code can divide by them
       ! later.
       cflag='O'
       skycnt=0.0
       skyerr=1.0
       skynos=1.0

       ! Form the sky box.  Note that the algorithm used here will not set the
       ! box symmetrically about the star if the star is near the array edge.
       ! If you want an symmetric box, you can play with the following few 
       ! lines of code, and the declaration of jbox.
       ! jbox=min(nint(xpos)-low(1), nint(ypos-low(2)))
       ! jbox=min(jbox, high(1)-nint(xpos))
       ! jbox=min(jbox, high(2)-nint(ypos))
       ! if (jbox < ibox) then 
       !   print*, 'Making box size ', jbox, 'was', ibox
       !   ibox=jbox
       ! end if
       ixslo=max(nint(xpos)-ibox/2,low(1))
       iyslo=max(nint(ypos)-ibox/2,low(2))
       ixshi=min(nint(xpos+ibox/2),high(1))
       iyshi=min(nint(ypos+ibox/2),high(2))

       ! Modal sky estimation.
       ! First we find the mean sky.
       icount=0
       allocate(clip_array((iyshi-iyslo+1)*(ixshi-ixslo+1)))
       skymen=0.0
       do j=iyslo, iyshi
         do i=ixslo, ixshi
           ! Are we in the object box?
           if ((real(i)-xpos)*(real(i)-xpos) + (real(j)-ypos) &
             *(real(j)-ypos) > 4.0*seeing*seeing) then
             if (pix_flg(i,j) == 'O') then
               icount=icount+1
               clip_array(icount)=data(i,j)
               !if (abs(data(i,j)) < 0.1) write(80,*) 1, i, j
             end if
           end if
         end do
       end do
       if (icount == 0) then
         cflag='I'
         deallocate(clip_array)
         return
       end if
       call clip_mean(clip_array(1:icount), skymen, skydev, iflag)
       if (iflag == 1) then
         cflag = 'B'
         deallocate(clip_array)
         return
       end if       
       if (debug) then
         print*, '@ clip_mean estimates sky is ', skymen, '+/-', skydev
         print*, '@ ', xpos, ypos
       end if
       deallocate(clip_array)

       ! Make a histogram of the sky pixel value distribution from
       ! mean-6sigma to mean+6sigma.  Thus the step size we require is;
       hstep = 12.0*skydev/real(maxbin)
       ! And the lowest value is;
       hmin=skymen-hstep*real(maxbin)/2.0
       if (intdat > 0) then
         ! Make sure the step is a multiple of an integer.
         hstep = real(intdat)*( int(hstep/real(intdat)) + 1 )
         ! And that the start value falls on an integer boundary.
         hmin = real(nint(hmin/hstep))*hstep
       end if
       first=.true.
470    if (debug) then 
         print*, '@ Histogram step size is ', hstep
         print*, '@ and the start value is ', hmin
       end if
       ! Now create the x array.
       do i=1, maxbin
         xhst(i)=hmin+real(i)*hstep
       end do
       ! And the y array.
       yhst=0.0
       n_reject=0
       do j=iyslo, iyshi
         do i=ixslo, ixshi
           ! Are we in the object box?
           if ((real(i)-xpos)*(real(i)-xpos) + (real(j)-ypos) &
             *(real(j)-ypos) > 4.0*seeing*seeing) then
             if (pix_flg(i,j) == 'O') then
               ! One way this can crash is if data is huge, then the value
               ! of nbin overuns an integer.
               nbin = nint(min(max(0.0,(data(i,j)-hmin)/hstep),real(maxbin+1)))
               if (nbin>=1 .and. nbin<=maxbin) &
               yhst(nbin)=yhst(nbin)+1.0
             else
               n_reject=n_reject+1
             end if
           end if
         end do
       end do
       ! Imagine if all the pixels are at values of n+0.25; they will all
       ! fall in the bin for value n, which skyfit will fit, giving a sky
       ! value that is in error by 0.25 of a pixel.  This cures the problem.
       ! First find the modal sky.
       iskymod = minval(maxloc(yhst))
       if (debug) print*, sum(yhst), ' pixels used ', n_reject, &
       ' dropped as flagged.'  
       if (sum(yhst) < n_reject) then
         ! Not enough data to fit.  For example they could all be on one side
         ! of the star.
         cflag='B'
         return
       end if
       ! Now total up all the X values in that bin.
       total = 0.0
       do j=iyslo, iyshi
         do i=ixslo, ixshi
           if ((real(i)-xpos)*(real(i)-xpos) + (real(j)-ypos) &
             *(real(j)-ypos) > 4.0*seeing*seeing) then
             if (pix_flg(i,j) == 'O') then
               nbin = nint(min(max(0.0,(data(i,j)-hmin)/hstep),real(maxbin+1)))
               if (nbin == iskymod) total = total + data(i,j)
             end if
           end if
         end do
       end do
       xhst = xhst + (total/yhst(iskymod)) - xhst(iskymod)
       ! A protection for the fitting routine.  If more than half of the
       ! counts have ended up in 1 bin, don't try and fit it (would you
       ! believe the answer anyway?).
       if (maxval(yhst) > 0.5*sum(yhst)) then
         cflag='B'
         skyerr=hstep
         skynos=hstep
         skycnt=yhst(minval(maxloc(xhst)))
         ! Don't even bother going around again.
         jflag=0
         sky_par(1)=0.0
       else
         ! Otherwise fit it.
         skynos=skydev
         skycnt=skwfit(xhst, yhst, maxbin, jflag, skyerr, skynos, sky_chi)
         ! Practically, the best the sky fitting routine can do is 0.1 of a 
         ! bin.
         skyerr = max(0.1*hstep,skyerr)
       end if
       if (jflag /= 0) then
         if (first) then
           ! It's worth changing the parameters a little to see if it'll
           ! fit then.
           hmin=hmin-3.0*hstep
           if (intdat > 0) then             
             hstep=hstep+real(intdat)
           else
             hstep=1.1*hstep
           end if                        
           first=.false.
           goto 470
         else
           cflag='B'
         end if
       end if
       
       ! Flag it if the chi-squared or the skewness of the sky fit was 
       ! too great to be reliable.
       if (sky_chi    > bad_sky_chi .and. bad_sky_chi > 0.0 ) cflag = 'I'
       if (sky_par(1) > bad_sky_skw .and. bad_sky_skw > 0.0 ) cflag = 'I'
       !print*, cflag, sky_chi, sky_par(1)
       !pause

       end subroutine skyfit

       subroutine set_bad_sky(bad_sky_chi_in, bad_sky_skw_in)
         real, intent(in) :: bad_sky_chi_in, bad_sky_skw_in

         bad_sky_chi=bad_sky_chi_in
         bad_sky_skw=bad_sky_skw_in

       end subroutine set_bad_sky
     
       subroutine skyhst(xdata, ydata, npts)

       ! Returns the histogram of the sky which was last fitted in skyfit,
       ! along with the number of data points.

       integer :: npts
       real, dimension(:) :: xdata, ydata

       integer :: i

       if (size(xdata) < maxbin) then
         print*, 'Programming error in calling S/R skyhst.'
         print*, 'xdata has less elements than the histogram.'
       else if (size(ydata) < maxbin) then
         print*, 'Programming error in calling S/R skyhst.'
         print*, 'ydata has less elements than the histogram.'
       else
         do i=1, maxbin
           xdata(i)=xhst(i)
           ydata(i)=ymax*yhst(i)
         end do
       end if
       npts=maxbin

       end subroutine skyhst

       
                                    
       real function skwfit(xdata, ydata, npts, icurf, skyerr, skynos, b_chi)

       ! The function fits the first npts points of the function ydata(xdata)
       ! with a skewed Gaussian, and then returns the fwhm and the peak
       ! position (skwfit).  icurf is the curfit status flag or set to -200
       ! if the fwhm became too small to fit with the binning given.  skynos is
       ! the sigma of the sky, on input it should be an estimate of this.

       ! The module marq is the marquardt routines.
       use marq

       ! Passed variables.

       integer, intent(in) :: npts
       integer, intent(out) :: icurf
       real, dimension(:), intent(in) :: xdata
       real, dimension(:), intent(inout) :: ydata
       real, intent(out) :: skyerr
       real, intent(inout) :: skynos
       real, intent(out) :: b_chi

       ! Locals.

       integer :: i, mode_w, npix, n_term, ixmax, iflag
       real :: d_sky_par(4), covar(4,4)
       real, allocatable, dimension(:) :: w
       real :: sum
       
       call reset()
       n_term=4
       ! Set skewness to zero.
       sky_par(1)=0.0
       ! And fix it
       iflag=fix_par(1, 4, n_term)
       ! And remove the limits which may have been set by earlier calls.
       iflag=limit_par(1,-1.0,1.0)
       ! Now limit the FWHM to be positive.
       iflag=limit_par( 4, 0.0, huge(sky_par(4)) )
       ! Set value of maximum.
       sum=0.0
       npix=0.0
       ymax=ydata(1)
       ixmax=1
       do i=1, npts
         if (ymax < ydata(i)) then
           ymax=ydata(i)
           ixmax=i
         end if
         sum=sum+ydata(i)*xdata(i)
         npix=npix+ydata(i)
       end do
       ! Set position.
       sky_par(2)=xdata(ixmax)
       ! And fix it.
       iflag=fix_par(2, 4, n_term)
       
       ! Set fwhm, to best guess based on sky deviation, but not less than
       ! the binning.
       sky_par(4)=max(skynos/0.416277,xdata(2)-xdata(1))
       
       ! Divide fit by ymax (to stop curfit overflows), and set weighting.
       if (debug) print*, '@ Dividing y values by ', ymax
       allocate(w(npts))
       !open(unit=23, file='hist.dat', status='unknown')
       do i=1, npts
         ! Make the errors on points with zero counts 1.
         w(i)=ymax*ymax
         ! Otherwise use the square root.
         if (ydata(i) > 0.0) w(i)=ymax*ymax/ydata(i)
         if (ydata(i) < 0.1*ymax) w(i)=0.0
         ydata(i)=ydata(i)/ymax
         !if (w(i) > tiny(w(i))) then
         !  write(23,*) xdata(i), ydata(i), sqrt(1.0/w(i))
         !else
         !  write(23,*) xdata(i), ydata(i), 0.0
         !end if
       end do
       !close(23)
       
       ! Set normalisation.
       sky_par(3)=1.0
       ! And fix it.
       iflag=fix_par(3, 4, n_term)
       
       ! print*, sky_par(1), sky_par(2), sky_par(3), sky_par(4), b_chi
       
       mode_w=-1
       icurf=1001
       if (debug) then
         icurf=100
         print*, '@ Parameters; skew, X(Y max), scaled Y max, fwhm.'
       end if
       call curfit(xdata, ydata, w, npts, mode_w, sky_par, d_sky_par, 4, &
       n_term, covar, 0.01, b_chi, -1.0, 20, icurf, skgauss, deriv_skg)
       
       ! If the sigma has gone smaller than the histogram step size
       ! it is dangerous to continue (the fit may crash).
       if (sky_par(4) < (xdata(npts)-xdata(2))/real(npts-1)) then
         icurf = -200
         if (debug) write(*,*) 'Binning was ', (xdata(npts)-xdata(2))/real(npts-1)
       else
         ! Free the normalisation.
         iflag=free_par(3, 4, n_term)
         ! And position.
         iflag=free_par(2, 4, n_term)       
         icurf=1001
         if (debug) then
           icurf=100
           print*, &
           '@ ** Now fitting with free position and normalisation.' 
           print*, '@ Parameters; skew, X(Y max), scaled Y max, fwhm.'
         end if
         call curfit(xdata, ydata, w, npts, mode_w, sky_par, d_sky_par, 4, &
         n_term, covar, 0.01, b_chi, -1.0, 20, icurf, skgauss, deriv_skg)
       
         if (icurf .ne. 0) then
           if (debug) print*, '@ Warning, CURFIT failed ', icurf
           ! Reset everything, and see if it can fit with non zero sky_par(1).
           sky_par(2)=xdata(ixmax)
           sky_par(3)=1.0
           sky_par(4)=max(skynos/0.416277,xdata(2)-xdata(1))
         end if

         ! If the sigma has gone smaller than the histogram step size
         ! it is dangerous to continue (the fit may crash).
         if (sky_par(4) < (xdata(npts)-xdata(2))/real(npts-1)) then
           icurf = -200
           if (debug) write(*,*) 'Binning was ', (xdata(npts)-xdata(2))/real(npts-1)
         else       
           ! Now add a little cockeyedness, make it similar to the setp used
           ! for calculating the derivatives in curfit.
           sky_par(1)=0.05
           ! And free it.
           iflag=free_par(1, 4, n_term)
           ! And force it positive.
           iflag=limit_par( 1, 0.0, huge(sky_par(1)) )
           ! Now limit the FWHM to be positive.
           iflag=limit_par( 4, 0.0, huge(sky_par(4)) )
       
           icurf=1001
           if (debug) then
             icurf=100
             print*, '@ ** Now fitting with free skew.'
           end if
           call curfit(xdata, ydata, w, npts, mode_w, sky_par, d_sky_par, 4, &
           n_term, covar, 0.0001, b_chi, -1.0, 20, icurf, skgauss, deriv_skg)
           if (icurf.ne.0 .and. debug) then
             print*, '@ Warning, CURFIT failed, flag was ', icurf
           end if

         end if

       end if
       
       if (debug) print*, '@ Errors were ', sqrt(abs(covar(1,1))), &
       sqrt(abs(covar(2,2))), sqrt(abs(covar(3,3))), &
       sqrt(abs(covar(4,4)))
       if (debug) print*,'@ final params', sky_par(1), sky_par(2), &
       sky_par(3), sky_par(4), b_chi
       
       skwfit=sky_par(2)
       skynos=0.416277*sky_par(4)
       skyerr=sqrt(abs(covar(2,2)))
       
       !open(unit=24, file='hist.fit', status='unknown')
       !do i=1, npts
       !  write(24,*) xdata(i), skgauss(xdata(i), sky_par, 4)
       !end do
       !close(24)

       deallocate(w)
              
       end function skwfit



       subroutine skymod(axfit, fit, npts)

       ! Returns the function fitted to the sky histogram.

       integer :: npts
       real, dimension(npts) :: fit, axfit

       integer :: i

       do i=1, npts
         axfit(i) = xhst(1) + real(i-1)*(xhst(maxbin)-xhst(1))/real(npts-1)
         fit(i)   = ymax*skgauss(axfit(i), sky_par, 4)
       end do

       end subroutine skymod


       
       real function skgauss(x, a, n_par)

       ! Returns a skewed Gaussian.

       integer :: n_par
       real :: x, a(n_par)

       real beta, x0, y0, fwhm
       real gamma, check, tmp
      
       beta=a(1)
       x0=a(2)
       y0=a(3)
       fwhm=a(4)
      
       if (abs(beta)<0.001 .or. (x-x0)<0.0) then
         ! Near enough a normal Gaussian, or below mid point.
         skgauss=y0*exp( -log(2.0) * ((2.0*(x-x0)/fwhm)**2.0) )
       else
         ! The function used here is a normal skewed Gaussian, except that
         ! we use the normal definition of fwhm and multiply it by
         ! sinh(beta)/(1.0 - exp(-1.0*beta)).  This means the left 
         ! hand side of the curve is the same as a Gaussian with that FWHM.
         gamma=fwhm*beta/(1.0 - exp(-1.0*beta))
         ! For the normal normalisation this would be
         ! gamma=fwhm*beta/sinh(beta)
         check=2.0*beta*(x-x0)/gamma
         if (check <= -1.0) then
           skgauss=0.0
         else
           tmp=log(1.0+check)
           skgauss=y0*exp( -log(2.0) * ((tmp/beta)**2.0) )
         end if
       end if
      
       end function skgauss
      
      
       subroutine deriv_skg(x, a, da, n_par, df_da)

       ! Returns the derivative of a skewed Gaussian.
      
       integer :: n_par
       real :: x, a(n_par), da(n_par), df_da(n_par)
      
       real, dimension(n_par) :: a1, a2
       real :: diff
       integer :: i, j
      
       param: do i=1, n_par
         do j=1, n_par
           a1(j)=a(j)
           a2(j)=a(j)
         end do
         if (i==1 .or. i==4) then
           ! Deal with the skew or FWHM.
           if (a(i) < 0.001) then
             ! The skew can't be negative.
             a1(i)=0.0
             a2(i)=0.005
             diff=-0.005
           else
             a1(i)=a1(i)-0.01*a(i)
             a2(i)=a2(i)+0.01*a(i)
             diff=-0.02*a(i)
           end if
         else if (i == 2) then
           ! Deal with x0.
           if (abs(a(i)-x) < 0.000001) then
             df_da(i)=0.0
             cycle param
           end if
           a1(i)=a1(i)-0.1*(a(i)-x)
           a2(i)=a2(i)+0.1*(a(i)-x)            
           diff=-0.2*(a(i)-x)
         else if (i == 3) then
           ! Normalisation,
           a1(i)=a1(i)-0.01*a(i)
           a2(i)=a2(i)+0.01*a(i)
           diff=-0.02*a(i)
         end if
         df_da(i)=skgauss(x,a1,n_par)-skgauss(x,a2,n_par)
         df_da(i)=df_da(i)/diff
       end do param

       end subroutine deriv_skg


       ! ** Here endeth the sky fitting routines.


       real function sum_flux(data, pix_flg, skycnt, skyerr, skynos, &
       a_par, cliprad, optnrm, adu, low, high, error, cflag)

       ! Inputs

       ! The data array.
       real, intent(in), dimension(:,:) :: data
       ! The pixel flag array.
       character, dimension(:,:), intent(in) :: pix_flg
       ! The estimated sky counts, the error in sky the determination, 
       ! and the sky noise.
       real, intent(in) :: skycnt, skyerr, skynos
       ! The profile shape specifiers.  (See function t_gauss for details.)
       !  1 and 2 are the FWsHM in orthogonal directions.
       !  3 specifies the angle of these directions wrt the pixel grid.
       !  4 is the peak flux.  Doesn't matter what is is (within reason),
       !    and will normally be set to the value resulting from the fit to
       !    the star.
       !  5 is the x position of the star.
       !  6 is the y position of the star.
       real, intent(in), dimension(:) :: a_par
       ! The radius at which the profile is clipped.
       real, intent(in) :: cliprad
       ! The peak flux of a star for which the extraction is to be normalised,
       ! divided by the RMS noise in the sky (skynos) for that star.  Should
       ! be set to zero for sky limited case normalisation.
       real, intent(in) :: optnrm
       ! The A-to-D conversion.
       real, intent(in) :: adu
       ! The active area of the detector.
       integer, intent(in) :: low(2), high(2)

       ! The outputs.
       ! The estimated error in sum_flux
       real, intent(out) :: error
       ! An error flag. Set to E if the star is off the edge of the frame, or
       ! the value of any pixel flags in the weight mask.  If the value is not
       ! 'O', then it is the value of the non 'O' flag nearest to the centre
       ! of the star.
       character, intent(out) :: cflag
                  
       ! Locals.
       ! The normalisation.
       real :: norm
       ! The distance of this pixel from the center of the profile.
       real :: dist, flg_dist
       ! Two counters.
       integer :: i, j
       ! The value of a Gaussian at a given point.
       real :: tgauss
       ! The weight mask and variance array.
       real, allocatable, dimension(:,:) :: weight, var_real
       ! The region being used.
       integer :: min1, max1, min2, max2
       ! The expected variance and profile.
       real :: var, profile
      

       ! The algorithm used is almost identical to that in Naylor (1998).
       ! The major change is that the variance array (V in the paper,
       ! var in the code here) is calculated as a fraction of the RMS of 
       ! the sky (skynos).  This is irrelevant for the weight mask, as 
       ! equation 10 of the paper shows it disappears in the normalisation.

       cflag='O'

       ! First calcuate the region of the image to be used.
       min1 = floor(a_par(5)-cliprad-0.5)
       max1=ceiling(a_par(5)+cliprad+0.5)
       min2 = floor(a_par(6)-cliprad-0.5)
       max2=ceiling(a_par(6)+cliprad+0.5)
       ! Uncomment these and comment the ones before to use all of the CCD image.
!       min1 = floor(a_par(5))
!       max1=ceiling(a_par(5))
!       min2 = floor(a_par(6))
!       max2=ceiling(a_par(6))

       if (min1<low(1) .or. min2<low(2) .or. &
         max1>high(1) .or. max2>high(2)) then

         ! The star is too close to the frame edge.  Put in values that
         ! won't give divide by zero errors.
         cflag='E'
         sum_flux=1.0
         error=-1.0

       else

         ! Calculate the real variances.
         allocate(var_real(min1:max1, min2:max2))
         var_real = skynos*skynos + (data(min1:max1,min2:max2)-skycnt)/adu
         ! If the data are very negative, this at least assures an answer.
         where (var_real <= 0.0) var_real = skynos*skynos

         ! Calculate the weight mask.
         norm=0.0
         flg_dist=huge(flg_dist)
         allocate(weight(min1:max1,min2:max2))
         do i=min1, max1
           do j=min2, max2
             dist=sqrt( (real(i)-a_par(5))*(real(i)-a_par(5)) + &
             (real(j)-a_par(6))*(real(j)-a_par(6)) ) 
             if (dist < cliprad+0.5) then
               ! This is a useful pixel.
               ! Find the estimated profile.
               tgauss = t_gauss(i, j, a_par, 2)
               ! What's the estimated variance here.
               var = skynos*skynos*(1.0 + optnrm*tgauss/(adu*a_par(4)))
               ! Calculate the value of the estimated profile at this point.
               ! Note that the value of the Gaussian at a point is used.
               ! Perfectionists may prefer to call the function i_gauss, that
               ! integrates over the pixel, but Tim has never found a case 
               ! where this helped.
               profile = tgauss/( 4.0*atan(1.0)*a_par(1)*a_par(2)*a_par(4) )
               if (dist > cliprad-0.5) profile = profile*(0.5-(dist-cliprad))
               weight(i,j) = profile / var           
               norm = norm + profile*profile / var
             else
               weight(i,j)=0.0
             end if
             ! Take the flag nearest to star centre.
             if (pix_flg(i,j)/='O' .and. weight(i,j)>0.0) then
               if (cflag=='O' .or. (cflag/='O' .and. dist<flg_dist)) then
                 cflag=pix_flg(i,j)
                 flg_dist=dist
               end if
             end if
           end do
         end do

         ! Normalise the weight mask.
         weight=weight/norm

         sum_flux = sum ( weight*(data(min1:max1,min2:max2)-skycnt) )
         ! Calculate the error, including that from the sky determination.
         error = sqrt( sum(var_real*weight*weight) + skyerr*skyerr*sum(weight) )

         deallocate(weight)
         deallocate(var_real)

       end if
          
       end function sum_flux
                  


       real function sum_aper(data, aprad, skycnt, skyerr, skynos, &
       xpos, ypos, adu, low, high, pix_flg, error, cflag)

       ! Performs aperture photometry.

       ! Inputs.
       ! The data array.
       real, intent(in), dimension(:,:) :: data
       ! The radius of the aperture, the sky counts, the
       ! estimated error in sky the determination, and the sky noise.
       real, intent(in) :: aprad, skycnt, skyerr, skynos
       ! The position of the center of the star.
       real, intent(in) :: xpos, ypos
       ! The A-to-D conversion.
       real, intent(in) :: adu
       ! The active area of the detector.
       integer, intent(in) :: low(2), high(2)
       ! The pixel flag array.
       character, dimension(:,:), intent(in) :: pix_flg

       ! The outputs.
       ! The estimated error in sum_aper
       real, intent(out) :: error
       ! An error flag set to 
       !    E if the star is off the edge of the frame.
       !    or the value of any pixel flags.
       character, intent(out) :: cflag
                  
       ! Locals.
       ! The distance of this pixel from the center of the profile.
       real :: dist
       ! Three counters.
       integer :: i, j
       ! The weight mask and variance array.
       real, allocatable, dimension(:,:) :: mask
       ! The region being used.
       integer :: min1, max1, min2, max2
       ! The maximum counts within the mask
       real :: ctsmax
       real :: flg_dist

       cflag='O'

       ! First calcuate the region of the image to be used.
       min1=floor(xpos-aprad-0.5)
       max1=ceiling(xpos+aprad+0.5)
       min2=floor(ypos-aprad-0.5)
       max2=ceiling(ypos+aprad+0.5)

       if (min1<low(1) .or. min2<low(2) .or. &
         max1>high(1) .or. max2>high(2)) then
         ! The star is too close to the frame edge.
         cflag = 'E'
         error = -1.0
         sum_aper = 0.0
         return
       end if

       ctsmax=data(nint(xpos),nint(ypos))
       ! Calculate the profile, putting it in mask.
       allocate(mask(min1:max1, min2:max2))
       flg_dist=huge(flg_dist)
       do i=min1, max1
         do j=min2, max2
           dist = sqrt( (real(i)-xpos)*(real(i)-xpos) + &
           (real(j)-ypos)*(real(j)-ypos) ) 
           if (dist < aprad+0.5) then
             ! This is a useful pixel.
             mask(i,j) = 1.0
             if (dist > aprad-0.5) then
               mask(i,j) = mask(i,j)*(0.5-(dist-aprad))
             end if
           else
             mask(i,j)=0.0
           end if
           if (mask(i,j) > 0.0) then
             ctsmax=max(data(i,j),ctsmax)
             ! Take the flag nearest to star centre.
             if (pix_flg(i,j)/='O' .and. mask(i,j)>0.0) then
               if (cflag=='O' .or. (cflag/='O' .and. dist<flg_dist)) then
                 cflag=pix_flg(i,j)
                 flg_dist=dist
               end if
             end if
           end if
         end do
       end do


       ! Find the total within the aperture.
       sum_aper = sum(mask*(data(min1:max1, min2:max2) - skycnt))
       ! Sum up the variance due to all the photons above the sky level.
       error = sum( (data(min1:max1,min2:max2)-skycnt)*mask*mask/adu)
       ! Now add the variance from calculated from the sky noise.
       error = error + sum(skynos*skynos*mask*mask)
       ! And finally the variance from the sky estimate.
       error = error + skyerr*skyerr*sum(mask)
       ! And change it into an error.  If a huge downward glitch makes the
       ! total in the mask negative, the abs makes sure the code doesn't 
       ! trip up.
       error=sqrt(abs(error))

       deallocate(mask)
          
       end function sum_aper


       subroutine gfit(data, fix_shape, dpos, fit_sub, skycnt, &
       skynos, adu, low, high, pix_flg, a_par, e_pos, cflag)
            
       ! Fits a profile to a star in the array data.
  
       ! The first guess of the parameters is in a_par.
       ! The subroutine returns a_par, the parameters of the fit.

       ! You can demand that all the shape parameters are fixed, and that
       ! only the position and normalisation are free using the logical
       ! fix_shape.  

       ! You can demand that the fit is subtracted from the data using the 
       ! logical fit_sub.

       ! You can fit as many profiles as you like, the number being set
       ! by the size of a_par, but the approximate error in position always
       ! applies to the first profile.

       ! The module marq is the marquardt routines.
       use marq   
              
       ! Inputs (data can be adjusted if the fit is subtracted).
       real, dimension(:,:), intent(inout) :: data
       logical, intent(in) :: fix_shape, fit_sub
       real, intent(in) :: dpos, skycnt, skynos, adu
       ! The active area of the array.
       integer, dimension(2) :: low, high
       ! The pixel flag array.
       character, dimension(:,:), intent(in) :: pix_flg
       
       ! Output.
       real, dimension(:), intent(inout) :: a_par
       ! The error in the X and Y position of the first profile.  Set to
       ! zero if the position is fixed.
       real, dimension(2), intent(out) :: e_pos
       ! A flag, normally O, 
       !         P if fit failed.
       !         or any pixel flag
       character, intent(out) :: cflag       
                  
       ! Locals.                              
       real, allocatable, dimension(:) :: x1, y1, w1
       real :: b_chi
       integer :: i, j, icount, icurf, jflag, ibad, ibox
       integer :: n_par, mode_w, n_term
       real, allocatable, dimension(:) :: da
       real, allocatable, dimension(:,:) :: covar
       integer :: ixcen, iycen
       real :: norm_fac
       
       ! The subset of the array used for fitting.
       integer :: ixbeg, ixend, iybeg, iyend, nsaxis(2)

       real :: dist, flg_dist

       ! The variables used when checking the rotation.
       ! real, allocatable, dimension(:) :: a_par_00
       ! real :: chi_00

       ! First we find out if the initial value of the peak flux 
       ! is less than zero.  Since we constrain it to be positive below
       ! a fit starting like this will just fail.
       if (a_par(4) <= 0.0) then
         e_pos=0.0
         cflag='P'
         return
       end if

       cflag='O'
       
       ! Set the position to the center of the sub-array.
       ixcen=nint(a_par(5))
       iycen=nint(a_par(6))
       
       ! Make a box whose length is twice the fwhm.
       ibox=int(2.0*sqrt(a_par(1)*a_par(2)*1.665*1.665)) + 1
       
       ! Put the data into the 1D array.
       ixbeg=max(low(1),ixcen-ibox)
       ixend=min(high(1),ixcen+ibox)
       iybeg=max(low(2),iycen-ibox)
       iyend=min(high(2),iycen+ibox)
       nsaxis(1)=ixend-ixbeg+1
       nsaxis(2)=iyend-iybeg+1
       naxis2=nsaxis(2)
       ! Allocate the 1D arrays.
       allocate(x1(nsaxis(1)*nsaxis(2)))
       allocate(y1(nsaxis(1)*nsaxis(2)))
       allocate(w1(nsaxis(1)*nsaxis(2)))
       
       icount=0
       ibad=0
       ! rewind(18)
       do i=ixbeg, ixend
         do j=iybeg, iyend
           icount=icount+1
           x1(icount) = real(icount)
           y1(icount) = data(i,j)-skycnt
           if (pix_flg(i,j) /= 'O') then
             w1(icount)=0.0
             ibad=ibad+1
             dist=sqrt( (real(i)-a_par(5))*(real(i)-a_par(5)) + &
             (real(j)-a_par(6))*(real(j)-a_par(6)) ) 
             if (cflag=='O' .or. (cflag/='O' .and. dist<flg_dist)) then
               cflag=pix_flg(i,j)
               flg_dist=dist
             end if
           else           
             ! Protect the next statement against huge values of skynos.
             if (10.0*skynos > sqrt(huge(skynos))) then
               w1(icount) = 0.0
             else
               w1(icount) = skynos*skynos + (data(i,j)-skycnt)/adu
               if (w1(icount) <= tiny(w1(icount))) then
                 ! This shouldn't happen, but if it does fudge the error so 
                 ! the user at least gets a fit.
                 w1(icount) = skynos*skynos 
               end if
               if (w1(icount) <= tiny(w1(icount))) then
                 ! Still dodgy, give up on this point.
                 w1(icount) = 0.0
               else 
                 w1(icount) = 1.0/w1(icount)
               end if
             end if
           end if
           ! write(*,*) i,j,icount,y1(icount),w1(icount)
         end do                                            
       end do

       if (sum(w1(1:icount)) < tiny(w1(1))) then
         if (cflag == 'O') cflag='P'
         e_pos=0.0
         return
         deallocate(x1, y1, w1)
       end if
       if (ibad > nint(3.142*a_par(1)*a_par(2)/(1.665*1.665))) then
         ! A large fraction of the seeing disc is bad pixels.
         ! cflag has already been set to what some of the cause was.
         e_pos=0.0
         return
         deallocate(x1, y1, w1)
       end if

       ! Normalise the arrays to a maximum value of 1.
       norm_fac=maxval(y1)
       if (norm_fac < tiny(norm_fac)) norm_fac=1.0
       y1=y1/norm_fac
       w1=w1*norm_fac*norm_fac
       ! And change any initial guesses of the flux we have.
       n_par=size(a_par)
       do i=4, n_par, 6
         a_par(i)=a_par(i)/norm_fac
       end do
       
       call reset()
       n_term=n_par
       ! if (3*(n_term/3)/=n_term .or. n_term<6) then
       !   print*, 'Programming error, n_term = ', n_term
       !   stop
       ! end if
       
       do i=4, n_par-2, 3
         a_par(i+1) = a_par(i+1) - real(ixbeg-1)
         a_par(i+2) = a_par(i+2) - real(iybeg-1)
       end do
       
       if (fix_shape) then
         jflag = fix_par(1, n_par, n_term)
         jflag = fix_par(2, n_par, n_term)
       else
         ! Limit the Gaussian width to be positive, more than 0.5 of a pixel,
         ! but not more than 10 times its current value, or twice the size
         ! of the box.
         jflag = limit_par( 1, 0.5, min(real(2*nsaxis(1)), 10.0*a_par(1)) )
         jflag = limit_par( 2, 0.5, min(real(2*nsaxis(2)), 10.0*a_par(2)) )
       end if
       ! print*, n_term
       ! Begin with the rotation fixed.
       jflag = fix_par(3, n_par, n_term)
       ! print*, n_term, dpos
       do i=4, n_par-2, 3
         ! Limit the flux to be positive.
         jflag=limit_par( i, 0.0, huge(a_par(i)) )
         ! And finally the limits on position.  Don't believe any fit
         ! which lies outside the box, or outside the constraint set 
         ! by the user.
         if (dpos > 0.0) then
           jflag = limit_par(i+1, max(a_par(i+1)-dpos, -0.5), &
                                  min(a_par(i+1)+dpos, real(nsaxis(1)+0.5)))
           jflag = limit_par(i+2, max(a_par(i+2)-dpos, -0.5), &
                                  min(a_par(i+2)+dpos, real(nsaxis(2)+0.5)))
         else            
           jflag = fix_par(i+1, n_par, n_term)
           jflag = fix_par(i+2, n_par, n_term)          
         end if
       end do
       
       mode_w=-1
       icurf=1001
       if (debug) then
         icurf=100
         print*, '@ ** Now fitting profile with rotation fixed.'
         print*, '@ initial parameters ', a_par
       end if
       
       allocate(da(n_par))
       allocate(covar(n_par,n_par))

       ! if (n_term == 1) print*, 'Warning n_term = 1.'
       call curfit(x1, y1, w1, nsaxis(1)*nsaxis(2), mode_w, a_par, da, n_par, &
       n_term, covar, 0.01, b_chi, -1.0, 20, icurf, gauss, d_gauss)

       if (debug) print*, '@ Returned from curfit.'

       if (.not. fix_shape) then
         ! Now let the rotation run free.
         jflag = free_par(3, n_par, n_term)
         jflag = limit_par(3, a_par(3)-0.78539816, a_par(3)+0.78539816)
         icurf=1001
         if (debug) then
           icurf=100
           print*, '@ ** Now fitting profile with rotation free.'
         end if
         call curfit(x1, y1, w1, nsaxis(1)*nsaxis(2), mode_w, a_par, da, n_par, &
         n_term, covar, 0.01, b_chi, -1.0, 20, icurf, gauss, d_gauss)
       end if

       ! print*, 'FWHM error', sqrt(abs(covar(3,3)))
       
       ! Undo the normalisation.
       do i=4, n_par, 6
         a_par(i)=a_par(i)*norm_fac
       end do

       !do i=1, nsaxis(1)*nsaxis(2)
       !  print*, i, x1(i), y1(i), w1(i), y1(i)-gauss(x1(i), a_par, n_par) + skycnt
       !end do

       if (fit_sub) then
         ! Make the output file the residuals.
         do i=1, nsaxis(1)*nsaxis(2)
           y1(i) = y1(i) - gauss(x1(i), a_par, n_par) + skycnt
         end do
       
         ! Recover the data from the 1D array.
         icount=0
         do i=ixbeg, ixend
           do j=iybeg, iyend
             icount=icount+1
             data(i,j) =  y1(icount)
           end do
         end do
       end if
       
       ! Return the correct values of a_par(5) and a_par(6)
       do i=4, n_par-2, 3
         a_par(i+1)=a_par(i+1)+real(ixbeg-1)
         a_par(i+2)=a_par(i+2)+real(iybeg-1)
       end do

       if (dpos > 0.0) then
         e_pos(1) = sqrt(abs(covar(5,5)))
         e_pos(2) = sqrt(abs(covar(6,6)))
       else
         e_pos=0.0
       end if
                     
       ! Deallocate the arrays.
       deallocate(da)
       deallocate(covar)
       
900    deallocate(x1)
       deallocate(y1)
       deallocate(w1)
       
       end subroutine gfit
       
       
       real function gauss(pix, a_par, n_par)

       ! This is the Gaussian used for fitting the profiles.

       ! The number of parameters.  There should be 3n+3 of these, where n is 
       ! the number of profiles.
       integer :: n_par
       ! The parameters, they are;
       ! a_par(1) - x width, when a_par(3) = 0.0
       ! a_par(2) - y width, when a_par(3) = 0.0
       ! a_par(3) - rotation
       ! a_par(4,7,10....) normalisation
       ! a_par(5,8,11....) x position
       ! a_par(6,9,12....) y position
       real :: pix, a_par(n_par)
       
       ! Locals
       real :: x, y, work, r, theta
       integer :: i, j, k

       ! Convert back from 1D pixels to 2D ones.
       i = nint(pix-1.0)/naxis2 + 1
       j = nint(pix) - (i-1)*naxis2
       
       gauss=0.0
       do k=4, n_par-2, 3
         x = real(i) - a_par(k+1)
         y = real(j) - a_par(k+2)
         ! Now we convert into r, theta co-ordinates.
         r = sqrt(x*x + y*y)
         ! If r is zero then the function will come out O.K. whatever theta is.
         if (abs(r) >= tiny(r)) then
           ! Otherwise find theta, getting the quadrant right.
           theta = sign(acos(x/r), y)
         else 
           ! Best set theta to something to avoid the NaN error.
           theta = 0.0
         end if
         ! And make the argument to the exponential in those terms.
         work = r*r * ( (cos(theta+a_par(3))/a_par(1))**2.0  &
                  +     (sin(theta+a_par(3))/a_par(2))**2.0 )
         gauss=gauss+a_par(k)*exp( -1.0 * (work) )
       end do
       
       end function gauss
       
       
       
       
       subroutine d_gauss(pix, a_par, da, n_par, df_da)
       
       integer :: n_par
       real :: pix
       real, dimension(n_par) :: a_par, da, df_da!, a_par2, a_par3
       
       real :: x, y, work, r, theta, work1, work2
       integer :: i, j, k

       ! Convert back from 1D pixels to 2D ones.
       i = nint(pix-1.0)/naxis2 + 1
       j = nint(pix) - (i-1)*naxis2
       
       do k=4, n_par-2, 3
       
       !   Find x and y.
         x = real(i) - a_par(k+1)
         y = real(j) - a_par(k+2)

         ! Convert into r and theta co-ordinates.
         r = sqrt(x*x + y*y)
         if (abs(r) >= tiny(r)) then
           theta = sign(acos(x/r), y)
         else
           ! Avoid the NaN error.
           theta=0.0
         end if
       
         work = r*r * ( (cos(theta+a_par(3))/a_par(1))**2.0  &
                  +     (sin(theta+a_par(3))/a_par(2))**2.0 )
         work = exp( -1.0 * (work) )
       
         if (k == 4) then
           ! Calculate the differentials w.r.t. the image shape parameters.
           df_da(1) = 2.0*a_par(4)*work*((r*cos(theta+a_par(3)))**2.0) / &
           (a_par(1)*a_par(1)*a_par(1))
           df_da(2) = 2.0*a_par(4)*work*((r*sin(theta+a_par(3)))**2.0) / &
           (a_par(2)*a_par(2)*a_par(2))
           df_da(3) = 2.0*a_par(4)*work*r*r*( &
           ( cos(theta+a_par(3))*sin(theta+a_par(3)) / (a_par(1)*a_par(1)) ) - &
           ( sin(theta+a_par(3))*cos(theta+a_par(3)) / (a_par(2)*a_par(2)) ) )
         end if
       
         ! Differential with respect to the normalisation.
         df_da(k) = work
         !a_par2=a_par
         !a_par3=a_par
         !a_par2(k+1)=a_par2(k+1)+0.01
         !a_par3(k+1)=a_par3(k+1)-0.01
         !df_da(k+1) = (gauss(pix, a_par2, n_par) - gauss(pix, a_par3, n_par))/0.02
         !a_par2=a_par
         !a_par3=a_par
         !a_par2(k+2)=a_par2(k+2)+0.01
         !a_par3(k+2)=a_par3(k+2)-0.01
         !df_da(k+2) = (gauss(pix, a_par2, n_par) - gauss(pix, a_par3, n_par))/0.02
         ! And w.r.t. the position.
         work1 = ( x*cos(a_par(3)) - y*sin(a_par(3)) ) / ( a_par(1)*a_par(1) )
         work2 = ( y*cos(a_par(3)) + x*sin(a_par(3)) ) / ( a_par(2)*a_par(2) )
         df_da(k+1) = 2.0*a_par(k) * work * &
         ( work2*sin(a_par(3)) + work1*cos(a_par(3)) )
         df_da(k+2) = 2.0*a_par(k) * work * & 
         ( work2*cos(a_par(3)) - work1*sin(a_par(3)) )
       
       end do
       
       end subroutine d_gauss
       
       
                        
       real function i_gauss(i, j, a_par, n_integ)
       
       ! Returns the integral under the pixel i, j
       
       ! Inputs.
       integer :: i, j
       real, dimension(:) :: a_par
       ! The number of points used for the integration is (2*n_integ +1) **2
       integer :: n_integ
       
       ! Locals.
       integer :: k, l
       real :: total, x, y, work, r, theta
       
       total=0.0
       do k=-n_integ, n_integ
         do l=-n_integ, n_integ 
           x = real(i) + real(k)/real(2*n_integ + 1) - a_par(5) 
           y = real(j) + real(k)/real(2*n_integ + 1) - a_par(6)
           ! Now we convert into r, theta co-ordinates.
           r = sqrt(x*x + y*y)
           ! If r is zero then the function will come out O.K. whatever theta is.
           ! Otherwise find theta, getting the quadrant right.
           if (abs(r) >= tiny(r)) then
             theta = sign(acos(x/r), y)
           else
             theta=0.0
           end if
           ! And make the argument to the exponential in those terms.
           work = r*r * ( (cos(theta+a_par(3))/a_par(1))**2.0  &
                    +     (sin(theta+a_par(3))/a_par(2))**2.0 )
           total = total + a_par(4)*exp( -1.0 * (work) )
         end do
       end do
       
       ! Normalise to points used.
       i_gauss = total / (real(2*n_integ + 1))**2.0
       
       
       end function i_gauss
       
       
       
       real function t_gauss(i, j, a_par, n_integ)
       
       ! Returns the value for the pixel i, j
       
       ! Inputs.
       integer :: i, j
       real, dimension(:) :: a_par
       ! The square root of the number of points used for the integration.
       ! (For compatability with i_gauss.)
       integer :: n_integ
       
       ! Locals.
       real :: x, y, work, r, theta
       
       x = real(i) - a_par(5) 
       y = real(j) - a_par(6)
       
       ! Now we convert into r, theta co-ordinates.
       r = sqrt(x*x + y*y)
       ! If r is zero then the function will come out O.K. whatever theta is.
       ! Otherwise find theta, getting the quadrant right.
       if (abs(r) >= tiny(r)) then
         theta = sign(acos(x/r), y)
       else
         theta=0.0
       end if
       ! And make the argument to the exponential in those terms.
       work = r*r * ( (cos(theta+a_par(3))/a_par(1))**2.0  &
                +     (sin(theta+a_par(3))/a_par(2))**2.0 )
       t_gauss = a_par(4)*exp( -1.0 * (work) )

       end function t_gauss

       

      subroutine clip_mean(data, mean, rms, iflag)

      real, dimension(:), intent(in)::data
      real, intent(out):: mean, rms
      integer, intent(out) :: iflag
      
      integer :: npoints, reject, i, iclip
      real :: mean_new, rms_new

      iflag=0

      npoints = size(data)
 
      if (npoints <= 1) then
        
        mean = 0.0
        ! Set the RMS to 1, just in case anyone divides by it.
        rms = 1.0
        if (debug) print*, '@ There were only ', npoints, ' points.'
        iflag=1
        
      else
        mean=sum(data)/real(npoints)
        rms=sqrt( sum((data-mean)**2.0)/real(npoints) )
        if (debug) print*, '@ First guess Mean and rms are ', mean, rms

        iclip=0
        clip: do
           reject=0
           mean_new=0.0
           rms_new = 0.0
           do i=1, npoints
             if (abs(data(i)-mean) < 2.0*rms) then
               mean_new = mean_new + data(i)
               rms_new  = rms_new  + (data(i)-mean)**2.0
             else
               reject=reject+1
             end if
           end do
           ! Have we rejected all the points?
           if (reject == npoints) then
             mean=0.0
             rms=1.0
             iflag=1
             if (debug) print*, '@ All data points rejected.'
             exit clip
           end if
           mean_new=mean_new/real(npoints-reject)
           if (abs(mean-mean_new) < tiny(mean)) exit clip
           ! (Remember the mean could be zero.)
           ! mean changed to abs(mean) to correct bug for negative skys.
           if (abs(mean-mean_new)/abs(mean) < 1.0e-06) exit clip
           iclip=iclip+1
           if (iclip > npoints) then
             print*, 'Stuck in endless loop in clip_mean, forcing exit.'
             exit clip
           end if
           ! Also jump out of the rms stars increasing again.
           rms_new = sqrt(rms_new/real(npoints-reject))
           if (rms_new>rms .and. iclip>1) exit clip
           rms=rms_new
           mean=mean_new
         end do clip

       end if

      end subroutine clip_mean

   end module opt_extr
