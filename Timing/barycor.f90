PROGRAM BARYCOR

!=============================================================================== 
! VERSION FOR RTPhoS - No User INPUT, all input controlled from run_RTPhos.py
!===============================================================================

! Program that takes a date and time and converts it into a Barycentric Dynamical
! Julian Date. The output values appear as:
!  BDJD, EXPOSURE, ERROR IN EXPOSURE
!
! The ERROR IN EXPOSURE is simply half the Exposure time of the corresponding frame. 
! The program assumes that the given time is the start of the exposure time
! it then shifts the time reported to the midle of the exposure.

! Compiles with:
!
! 	g95 -o barycor barycor.f90 libjpl.a /usr/lib/libf95.a
!
! Where:libjpl.a is the JPL ephemeris library (from STARLINK)
!       libf95.a is the gfortran (g95)library
!
! Any Fortran compiler will compile the code as long as the two libraries above
! are specified or are in the LD path of your system.
! If the compiler complaines for undefined references that probably means that
! your system does not have the libg2c library. Install (libg2c.so.0) and then
! either include in the LD path or link it directly using: -L/path-to-libg2c.so -lg2c
!
! Additional Requirements: 
! A file with a list of leap seconds (here leap.dat from
! the ARK Software Package ($ARKDIR/data/leap.dat). This code assumes tha the
! file is called leapdat. leapdat is a symbolic link to leap.dat created together
! with the symbolic link JPLEPH which points to the jpleph.dat file 
! (in my case: /opt/star-kapuahi/etc/jpleph.dat).
!
! This code runs through the run_RTPhoS Python script which creates the symbolic 
! links required, runs the code and then removes the links. 

  implicit none

! Variable declerations.
! ************************************************
  integer :: i
  integer :: num, year, month, day, hrs, mins, secs
  integer :: ihrs, imins, isecs, isig

  double precision :: ddays, utcorr, exposure

  integer :: jy,jm,jd,julian,idat,jdi
  real :: rah,ram,ras,jdr,jdrr
  real :: decd,decm,decs,sign
  double precision :: ra,dec,jdd,equin,BDT, hjd
  double precision :: V
! ************************************************

! Input
  read*, equin, rah, ram, ras, decd, decm, decs
  read*, year, month, day, hrs, mins, secs, utcorr, exposure

! Convert Right Ascension and Declination to decimal form.
  ra=dble(rah)+(dble(ram)+dble(ras)/60.0d0)/60.0d0
  if (decd .lt. 0) then
     dec=-( dble(abs(decd)) + ((dble(decm) + dble(decs)/60.0d0)/60.0d0) )
  else
     dec=dble(decd)+(dble(decm)+dble(decs)/60.0d0)/60.0d0
  end if

! Check that R.A. value is acceptable.
  if (ra<0.0.or.ra>360.0)then
     print*,'***R.A. NOT ALLOWED***'
!    Do something here to make sure the pipeline does not crash.
  end if

! Check that Declination value is acceptable.
  if (dec<-90.0.or.dec>90.0)then
     print*,'***DEC. NOT ALLOWED***'
!    Do something here to make sure the pipeline does not crash.
  end if

! Now do the conversion from UT time to BDJD.

! Calculate decimal day.
  ddays = (secs/86400.0d0) + (mins/1440.0d0) + (hrs/24.0d0)
! Apply Time Zone UT Correction.     
  ddays = ddays + (utcorr/24.0d0)

! Convert date and time to Julian Day.
  jdi = julian(year, month, day)
  jdrr = real(jdi)-0.5
  jdd = ddays + dble(jdrr)

! Perform the barycentric correction of the exposure start time.
  call dbary(jdd, ra, dec, equin, bdt, v)

! Ask how many decimal places to store and write the file header.
!  isig = (int(bdt/10.0d0**dble(isig)))*10.0d0**dble(isig)
!  bdt = bdt - dble(isig)

! Output
  print*, bdt+((dble(exposure)/86400.d0)/2.0d0), (dble(exposure)/86400.d0)/2.0d0

END PROGRAM BARYCOR

!---------------------------------------------------------------------

       Integer Function JULIAN( Year, Month, Day )

       Integer Year, Month, Day
       Integer YY, MM, DD

       If( Month .eq. 1 .or. Month .eq. 2 ) Then
         MM = Month + 9
         YY = Year - 1
       Else If( Month .ge. 3 .or. Month .le. 12 ) Then
         MM = Month - 3
         YY = Year
       Else
         Stop 'ERROR:: Non-existent month!'
       End If
       DD = Day + ( YY + 4712 ) * 1461 / 4 + ( MM * 306 + 5 ) / 10 + 59
       DD = DD - ( YY / 100 + 49 ) * 3 / 4 + 38

       JULIAN = DD

       End


!---------------------------------------------------------------------
 
        subroutine dbary(JD,ra,dec,equin,BDT,V)
       
!       Subroutine to account for leap seconds and calculate
!       Barycentric Correction.  MWS(19/7/95)

!       Need to input UTC in julian date, will convert to BDT
!       (NOT G.R. CORRECTIONS), V is radial velocity (km/s)
!       correction, not including rotation of the earth.
!       RA in decimal hours and Dec in decimal degrees.
!       Equinox in decimal years

        Double Precision ra,dec,JD,dt,dist,temp,equin
        Double Precision p_ra,p_dec

!       Components of direction vector.
        Double Precision b_(3)

!       Array which receives x,y,z,xdot,ydot,zdot (in AU and AU/day)
        Double Precision PV(6),TDT,BDT,V

        Double Precision Pi, Deg_Rad, Rad_Deg, U_L
        Common / math00 / Pi, Deg_Rad, Rad_Deg, U_L
        external barydis
!
        Pi = Atan( 1.0 ) * 4.0
        Deg_Rad = 1.8D+02 / Pi
        Rad_Deg = Pi / 1.8D+02
        U_L = 3.662422D+02 / 3.652422D+02


!       Precess coordinates to 2000.0 so PLEPH works
!       MJD(2000.0)=11545.0

!       print*,ra,dec,equin
        call prec(ra,dec,equin,11545.0d0,p_ra,p_dec)
!       print*,p_ra,p_dec,equin

!       To move from ra in hours to RA in degrees.
        p_ra=p_ra*15.0d0
        temp = cos(p_dec*rad_deg)
        b_( 1 ) = temp * cos( p_ra*rad_deg )
        b_( 2 ) = temp * sin( p_ra*rad_deg )
        b_( 3 ) = sin( p_dec*rad_deg )

        call barydis(JD,PV)

        dist=PV(1)*b_(1)+PV(2)*b_(2)+PV(3)*b_(3)
!       dist in AU, 499 etc is the light travel time for 1AU.
        dt=(dist*499.004782d0)

!       To calculate radial velocity shift
        V=PV(4)*b_(1)+PV(5)*b_(2)+PV(6)*b_(3)
!       To convert to km/s from AU/day
        V=V*(1.495978706d8/(24.0d0*3.6d3))

!       To account for leap seconds
        call leap(JD,TDT)
!       the next term corrects for the time the signal passed through
!       the barycenter.

        BDT=TDT+(dt/8.64d4)
!       print*,JD,dt,BDT
        end

!------------------------------------------------------------------------------

       subroutine barydis(JD,PV)

!       Programme to calculate the vector of the Barycenter-Earth vector in 
!       a cartesian right handed system with the x-axis aligned to the first
!       point in Aries and the z-axis to the North celestial pole.  Uses
!       Starlink/JPL package PLEPH (SUN/87.5)
!       M.Somers 17/7/1995

!       Array which receives x,y,z,xdot,ydot,zdot (in AU and AU/day)
       Double precision PV(6)

!       Barycentric Dynamical Time expressed as Modified Julian Date
!       (JD-2400000.5)
       Double Precision MBDT
       Double Precision JD

!       NP-The body whose coordinates are required.
!       NC-The origin of the coordinate system.
!       1=Mercury
!       2=Venus
!       etc..
!       9=Pluto
!       10=Moon
!       11=Sun
!       12=Solar System Barycenter
!       13=Earth-Moon Barycenter
       Integer NP,NC

!       true=success,false=TDB out of range or illegal NP,NC.
       Logical OK
       MBDT=JD-2400000.5
       NP=3
       NC=12
       Call PLEPH(MBDT,NP,NC,PV,OK)
!       print*,JD,MBDT,PV

       end

!------------------------------------------------------------------------------
       Subroutine PREC( Alpha, Delta, Epoch, MJD, p_Alpha, p_Delta )

       Implicit None

!       Precession programme using the "rigorous formula" in the
!       Astronomical Almanac (86).  Alpha in decimal hours, Delta
!       in decimal degrees, MJD = JD - 2440000, Epoch is in years

       Double Precision Alpha, Delta              ! initial coordinates
       Double Precision Epoch, MJD
       Double Precision p_Alpha, p_Delta       ! output coordinates

       Double Precision r_Alpha, r_Delta       ! in radian
       Double Precision Zeta, Zed, Theta, T
       Double Precision cos_Theta, sin_Theta
       Double Precision cos_Delta, sin_Delta
       Double Precision X, Y, Z
       Double Precision Work

       Double Precision Pi, Deg_Rad, Rad_Deg, U_L
       Common / math00 / Pi, Deg_Rad, Rad_Deg, U_L

!              Alpha, Delta to radians
       r_Alpha = ( Alpha * 1.5D+01 ) * Rad_Deg
       r_Delta = Delta * Rad_Deg
!              convert to J2000.0
       T = ( Epoch - 2.0D+03 ) * 1.0D-02
       Call PREC_ANGLES( T, Zeta, Zed, Theta )
       r_Alpha = r_Alpha - Zed
       cos_Theta = Cos( Theta )
       sin_Theta = Sin( Theta )
       cos_Delta = Cos( r_Delta )
       sin_Delta = Sin( r_Delta )
       Work = Cos( r_Alpha )
       X = Sin( r_Alpha ) * cos_Delta
       Y = sin_Theta * sin_Delta + Work * cos_Theta * cos_Delta
       Z = cos_Theta * sin_Delta - Work * sin_Theta * cos_Delta
       Call THREED_ANGLES( X, Y, Z, r_Alpha, r_Delta )
       r_Alpha = r_Alpha - Zeta
!              Next change to the desired epoch
       T = ( MJD - 1.1545D+04 ) / 3.6525D+04
       Call PREC_ANGLES( T, Zeta, Zed, Theta )
       r_Alpha = r_Alpha + Zeta
       cos_Theta = Cos( Theta )
       sin_Theta = Sin( Theta )
       cos_Delta = Cos( r_Delta )
       sin_Delta = Sin( r_Delta )
       Work = Cos( r_Alpha )
       X = Sin( r_Alpha ) * cos_Delta
       Y = Work * cos_Theta * cos_Delta - sin_Theta * sin_Delta
       Z = Work * sin_Theta * cos_Delta + cos_Theta * sin_Delta
       Call THREED_ANGLES( X, Y, Z, r_Alpha, r_Delta )
       r_Alpha = r_Alpha + Zed + Pi * 2.0
       r_Alpha = Mod( r_Alpha, Pi * 2.0 )
!              change into decimal hours/degrees
       p_Alpha = r_Alpha * Deg_Rad / 1.5D+01
       p_Delta = r_Delta * Deg_Rad

       End



       Subroutine PREC_ANGLES( T, Zeta, Zed, Theta )

       Implicit None

       Double Precision T       ! Time since J2000.0
       Double Precision Zeta, Zed, Theta

       Zeta = ( ( 8.72664626D-08 * T + 1.464331242D-06 ) * T+1.118086019D-02 ) * T
       Zed = ( ( 8.901179185D-08 * T + 5.307546255D-06 ) * T+1.118086019D-02 ) * T
       Theta = ( ( -2.024581932D-07 * T - 2.068215164D-06 ) * T+9.71717297D-03 ) * T

       End



       Subroutine THREED_ANGLES( X, Y, Z, r_Alpha, r_Delta )

       Implicit None

       Double Precision X, Y, Z       ! Orthogonal Vector
       Double Precision r_Alpha, r_Delta       ! Radians

       Double Precision cos_Delta

       Double Precision Pi, Deg_Rad, Rad_Deg, U_L
       Common / math00 / Pi, Deg_Rad, Rad_Deg, U_L

       r_Delta = Asin( Z )
       cos_Delta = Sqrt( 1.0 - Z * Z )
       If( X .ge. 0.0 ) Then
         r_Alpha = Acos( Y / cos_Delta )
       Else
         r_Alpha = Pi * 2.0D+00 - Acos( Y / cos_Delta )
       End If

       End


!------------------------------------------------------------------------------
       subroutine leap(JD,TDT)

       Double Precision JD,TDT
       real jdr,leapcount,jdrold,lcold,jd1
       character*50 last
       integer iunit

       iunit=23
       open(unit=iunit, file='leapdat', status='old')
       read(iunit,*)
       read(iunit,*)
       read(iunit,'(a)') last
       read(iunit,*,end=20) jdr,leapcount
       jd1=jdr
       if (JD.lt.dble(jd1) ) then
         print*,'WARNING:data is before the first leap second on'
        print*,' file: output is Barycentric Julian Date '
        print*,'as opposed to Barycentric Dynamical Time.'
        TDT=JD
        goto 30
        end if

10       jdrold=jdr
       lcold=leapcount
       read(iunit,*,end=20) jdr,leapcount
       if ((JD.gt.dble(jdrold)).and.(JD.le.dble(jdr)))then
!         If you are given a time in UT which is when a leap second
!         occurs, it os unclear what the ephemeris time is (since UT
!         (stops for 1 second).  Using .le. as opposed to .lt. above
!         makes an arbitary decision.
         TDT=JD+(dble(lcold+32.184d0)/8.64d4)
         goto 30
       end if
       goto 10

20       TDT=JD+(dble(lcold+32.184d0)/8.64d4)
!       if (JD.gt.dble(jdrold+730.0) ) then
!        print*,'WARNING: the data is over 2 years after the last'
!        print*,'leap second on file.'
!        print*,'leap second file last updated:'
!        print*,last
!       end if       
30       close(unit=iunit)
       end

!------------------------------------------------------------------------------

