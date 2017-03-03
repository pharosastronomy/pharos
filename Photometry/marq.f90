module marq
 
   ! These are Koji Mukai's versions of the Marquardt routines.  
   ! For a good reference as to how it works, see Numerical Recipes.
   ! The conversion to Fortran 90 was originally a bit basic, but it 
   ! was improved by Tim Naylor and Lee Howells.

   implicit none

   ! Define large numbers.
   Real, Parameter, private :: n_Large=-1.0*huge(n_Large), p_Large=huge(P_Large)

   integer, private, parameter :: m_Par=32

   ! "The taming factor".
   Real, private :: tame_Factor = 1.0

   Logical, private :: cyclic_Data = .False.
   Real, private :: cycle_Data = 0.0

   Integer, dimension(m_Par), private :: list_Cyclic=0
   Real, dimension(m_Par), private :: tab_Period=0.0, tab_Bottom=0.0

   ! Initializes the free_list common block so that no parameter is fixed.
   Integer, dimension(m_Par), private :: Index = (/ 1, 2, 3, 4, 5, 6, 7, 8, &
   9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, &
   27, 28, 29, 30, 31, 32 /)
   Integer, dimension(m_Par), private :: Table_Free=0

   ! Initializes the link_status so that no parameter is linked.
   Integer, dimension(m_Par), private :: l_Index = (/ 1, 2, 3, 4, 5, 6, 7, 8, &
   9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, &
   27, 28, 29, 30, 31, 32 /)
   Real, dimension(m_Par), private :: l_Offset=0.0
   
   ! Set up the default boundary within which parameters are confined.
   Real, dimension(m_par), private :: lo_A = n_Large, hi_A = p_Large
     
   Logical Debug
   Common / bug / Debug
       
   contains

       Subroutine COV_SORT( Covar, n_Par, n_Term )
 
       ! Repacks the covariance matrix into convenient form

       Real, intent(inout), dimension(:,:) ::  Covar
       Integer, intent(in) :: n_Par, n_Term

       Real Work
       Integer I, J, K, L

       ! Zero all elemsnts below diagonal
       Do J = 1, n_Par - 1
         Do I = J + 1, n_Par
           Covar( I, J ) = 0.0
         End Do
       End Do
                             ! Repack off-diagnoal to correct location
       Do I = 1, n_Term - 1
         Do J = I + 1, n_Term
           K = Index( I )
           L = Index( J )
           If( L .gt. K ) Then
             Covar( L, K ) = Covar( I, J )
           Else
             Covar( K, L ) = Covar( I, J )
           End If
         End Do
       End Do
                             ! Swap diagnoals using the top row
       Work = Covar( 1, 1 )
       Do J = 1, n_Par
         Covar( 1, J ) = Covar( J, J )
         Covar( J, J ) = 0.0
       End Do
       K = Index( 1 )
       Covar( K, K ) = Work
       Do J = 2, n_Term
         K = Index( J )
         Covar( K, K ) = Covar( 1, J )
       End Do
                             ! Finally, fill in above diagonals
       Do J = 2, n_Par
         Do I = 1, J - 1
           Covar( I, J ) = Covar( J, I )
         End Do
       End Do
 
       End Subroutine COV_SORT



       Subroutine CURFIT( X, Y, W, Points, mode_W, A, dA, n_Par, &
                                  n_Term, Covar, d_Chi, b_Chi, &
                                 sig_Clip, m_Rep, Flag, FUNCTN, FDERIV )
 
       ! A new version of CURFIT; see write-up

       ! Number of data points
       Integer, intent(in) :: Points
       ! X-data, Y-data
       Real, dimension(Points), intent(in) :: X, Y
       ! Weight
       real, dimension(Points), intent(inout) ::  W
       ! Weight mode.  1 Poisson, 0 equal weight, -1 explicit
       Integer, intent(in) :: mode_W
       ! Number of parameters
       Integer, intent(in) :: n_Par
       ! The parameters
       Real, intent(inout) :: A( n_Par )
       ! Step for numerical differenciation
       Real, intent(in) :: dA( n_Par )
       ! Number of free parameters
       Integer, intent(in) :: n_Term
       ! The covariance matrix
       Real, intent(out) :: Covar(:,:)
       ! Fractional Chi change for convergence
       Real, intent(in) :: d_Chi
       ! The final (best) reduced Chi squared
       Real, intent(out) :: b_Chi
       ! Iterative clipping flag.  If >0, do iterative clipping
       Real, intent(in) :: sig_Clip
       ! Maximum number of trials
       Integer, intent(in) :: m_Rep
       ! Controls output; error handling.
       Integer, intent(inout) :: Flag
       !   0: ok; 
       !  -1: not converged; 
       !  -2: lambda too big;
       ! -51: Too many iterations in svd_comp
       ! -91: initial parameter out of range; 
       ! -92: no D.O.F.
       ! -93: chi-squared is zero/
       ! User-supplied model function and its derivatives.       
       real, external :: FUNCTN
       external FDERIV

       Integer :: n_Free, J, K, Reject, flag_mrq
       Real Lambda, Temp
       Real o_Chi, Chi, Thres, imp_Chi
       Logical Interactive

!       Call SETBUG( )
       Debug = .False.

       ! Set all the outputs to real numbers, so alternate returns don't
       ! cause NaN errrors.
       covar=1.0
       b_chi=0.0
       
       If( Flag .gt. 1000 ) Then
         Interactive = .False.
       Else
         Interactive = .True.
       End If
       Flag = 0

       Do K = 1, n_Par
         If( A( K ) .le. lo_A( K ) .or. A( K ) .ge. hi_A( K ) ) Then
           Flag = -91
           Return
         End If
       End Do

       Do K = 1, n_Par
         If( list_Cyclic( K ) .eq. 1 ) Then       ! actively cyclic
           If( A( K ) .lt. tab_Bottom( K ) ) Then
             Write( *, 90 )
 90             Format( ' WARNING:: Parameter #', I2, ' out of ''cyclic'' range' / &
                      '           It will be modified before fitting...' )
             A( K ) = tab_Bottom( K ) + tab_Period( K ) &
                           + Mod( A( K ) - tab_Bottom( K ), tab_Period( K ) )
           Else If( A( K ) .ge. tab_Bottom( K ) + tab_Period( K ) ) Then
             Write( *, 90 )
             A( K ) = tab_Bottom( K ) &
                          + Mod( A( K ) - tab_Bottom( K ), tab_Period( K ) )
           End If
         End If
       End Do

       Reject = 1
       If( mode_W .gt. 0 ) Then       ! Poisson weighting
         Do J = 1, Points
           If( Y( J ) .gt. 0 ) Then
             W( J ) = 1.0 / Y( J )
           Else                            ! Ignore negative points
             W( J ) = 0.0
           End If
         End Do
       Else If( mode_W .eq. 0 ) Then       ! Equal weighting
         Do J = 1, Points
           W( J ) = 1.0
         End Do
       ! Else                            ! Predetermined (W, not sigma)
       End If

       Do While( Reject .gt. 0 )       ! Iterative Clipping Loop
         imp_Chi = 1000.0
                                   ! First, degree of freedom
         n_Free = 0
         Do J = 1, Points
           If( abs(W(J)) > tiny(W(J)) ) Then
             n_Free = n_Free + 1
           End If
         End Do
         n_Free = n_Free - n_Term
         If( n_Free .le. 0 ) Then
           Write( *, 120 ) Points, n_Term, n_Free
120           Format( ' ERROR:: No degrees of freedom in CURFIT' / &
                    '         # of points =', I5, ' # of free parameters =', &
                    I2, ' DOF =', I3 )
           Flag = -92
           Return
         End if
         ! First call to MARQUARDT
         Lambda = -0.001
         ! Make sure there's no junk in covar.
         covar = 0.0
         Call MARQUARDT( X, Y, W, Points, Lambda, A, dA, Covar, &
                            o_Chi, n_Par, n_Term, FUNCTN, FDERIV, flag )
         if (flag /= 0) return
         o_Chi = o_Chi / n_Free
         If( Interactive ) Then
           Write( *, 520 ) 0, o_Chi, Lambda, ( A( J ), J = 1, n_Par )
520           Format( ' Iteration # ', I2, ' Red. Chi-squ = ', 1PE10.3, &
                          ' Lambda = ', E8.1 / ( 5( ' ', E10.3 ) ) )
         End If

         ! Iteration
         Do K = 1, m_Rep
           Call MARQUARDT( X, Y, W, Points, Lambda, A, dA, &
                                 Covar, Chi, n_Par, n_Term, &
                                                      FUNCTN, FDERIV, flag )
           if (flag /= 0) return
           Chi = Chi / n_Free
           If( Interactive ) Then
             Write( *, 520 ) K, Chi, Lambda, ( A( J ), J = 1, n_Par )
           End If
           if (chi <= tiny(chi)) then
             print*, ' WARNING:: Chi-squared is zero.'
             flag = -93
             goto 300
           end if
           ! Small bug (which has a huge effect) spotted by Timn and Tariq.
           ! This line used to read .gt., which meant that if the fit had
           ! converged so well that o_chi was equal to chi you never tested
           ! against the convergence condition.
           ! If( o_Chi .ge. Chi ) Then
           ! But even this caused the g95 compiler problems, so I had to
           ! use this slightly laxer condition.
           if (o_chi - chi >= -1.0*epsilon(chi)) then
             imp_Chi = ( o_Chi - Chi ) / Chi
             If( imp_Chi .lt. d_Chi ) Goto 300
             o_Chi = Chi
           End If
           If( Lambda .ge. 1.0E+30 ) Then
             print*, o_chi, chi, o_chi-chi
             Write( *, 530 ) K, imp_Chi
530             Format( ' WARNING:: After ', I3, 'th iteration, lambda', &
                      ' is becoming too big' / &
                      '           The last fractional improvement in reduced', &
                      ' chi-squared is ', 1PE8.2 )
             If( imp_Chi .lt. 1.0E-03 ) Then
              Write( *, 531 )
531              Format( '           I would consider it as having converged!' )
             Else If( imp_Chi .lt. 0.1 ) Then
              Write( *, 532 )
532              Format( '           I might consider it as having converged' )
             Else
              Write( *, 533 )
533              Format( ' ERROR:: Curfit should not behave like this!' )
             End If
             Flag = -2
             Goto 300
           End If
         End Do
         ! .not. converged
         Write( *, 540 ) m_Rep
540         Format( ' WARNING:: Not converged after ', I2, ' iterations!' )
         Flag = -1

300      Continue
         if (debug) print*, '@ Converged, with flag ', flag
         b_Chi = Chi
         ! Check for clipping
         Reject = 0
         If( sig_Clip .gt. 0.0 ) Then
           Thres = Chi * sig_Clip * sig_Clip
           Do J = 1, Points
             Temp = Y( J ) - FUNCTN( X( J ), A, n_Par )
             Temp = Temp * Temp * W( J )
             If( Temp .gt. Thres ) Then
              Reject = Reject + 1
              W( J ) = 0.0
             End If
           End Do
         End If
         If( Reject .gt. 0 ) Then
           Write( *, 550 ) Reject
550           Format( ' * ', I3, ' points are rejected as outliers.' / &
                          '   Restarting fit again...' )
         End If
       End Do
       ! Final call...calculate covariances
       if (debug) print*, '@ Calculating co-variances.'
       Lambda = 0.0
       Call MARQUARDT( X, Y, W, Points, Lambda, A, dA, Covar,&
                             Chi, n_Par, n_Term, FUNCTN, FDERIV, flag_mrq )
       if (debug) print*, '@ done'
       if (flag == 0) flag = flag_mrq
       Do K = 1, n_Par
         If( list_Cyclic( K ) .ne. 0 ) Then       ! Cyclic
           If( A( K ) .lt. tab_Bottom( K ) ) Then
             A( K ) = tab_Bottom( K ) + tab_Period( K ) &
                           + Mod( A( K ) - tab_Bottom( K ), tab_Period( K ) )
           Else
             A( K ) = tab_Bottom( K ) &
                          + Mod( A( K ) - tab_Bottom( K ), tab_Period( K ) )
           End If
         End If
       End Do

       End Subroutine CURFIT



       Integer Function CYC_DAT( Cycle )

       Real Cycle

       If( Cycle .le. 0.0 ) Then
         Write( *, 100 )
100         Format( ' ERROR:: unacceptable periodicity' )
         cyclic_Data = .False.
         cycle_Data = 0.0
         CYC_DAT = -1
       Else
         cyclic_Data = .True.
         cycle_Data = Cycle
         CYC_DAT = 1
       End If

       End Function CYC_DAT



       Integer Function CYC_PAR( N, Lo, Hi )

       ! Select the relevant parameter
       Integer N
       ! Specify the range for cyclic behaviour
       Real Lo, Hi       

       Integer K

       If( Hi .le. Lo ) Then
         Write( *, 100 )
100         Format( ' ERROR:: Range of cyclic behaviour is not positve' )
         CYC_PAR = -1
       Else If( N .le. 0 .or. N .gt. m_Par ) Then
         Write( *, 110 )
110         Format( ' ERROR:: A non-existent parameter has bee specified' )
         CYC_PAR = -1
       Else If( list_Cyclic( N ) .ne. 0 ) Then
         Write( *, 120 ) N
120         Format( ' ERROR:: Parameter ', I2, ' is already declared cyclic' )
         CYC_PAR = -1
       Else If( l_Index( N ) .ne. N ) Then
         Write( *, 130 ) N
130         Format( ' ERROR:: Parameter ', I2, ' is linked' )
         CYC_PAR = -1
       Else If( lo_A( N ) .ne. n_Large .or. hi_A( N ) .ne. p_Large ) Then
         If( hi_A( N ) - lo_A( N ) .ge. Hi - Lo ) Then
           Write( *, 140 )
140           Format( ' ERROR:: Conflict between LIMIT_PAR and CYC_PAR' )
           CYC_PAR = -1
         Else
           Do K = 1, m_Par
             If( l_Index( K ) .eq. -N ) Then
              Write( *, 150 )
150               Format( ' ERROR:: Illegal attempt ', &
                                         'to use CYC_PAR on a LINKed parameter' )
              CYC_PAR = -1
              Return
             Else If( l_Index( K ) .eq. N .and. K .ne. N ) Then
              list_Cyclic( K ) = -1
              tab_Bottom( K ) = Lo
              tab_Period( K ) = Hi - Lo
             End If
           End Do
           Write( *, 160 )
160           Format( ' WARNING:: This parameter has already been LIMITed' )
           list_Cyclic( N ) = -1
           tab_Bottom( N ) = Lo
           tab_Period( N ) = Hi - Lo
           CYC_PAR = 0
         End If
       Else
         Do K = 1, m_Par
           If( l_Index( K ) .eq. -N ) Then
             Write( *, 160 )
             CYC_PAR = -1
             Return
           Else If( l_Index( K ) .eq. N .and. K .ne. N ) Then
             list_Cyclic( K ) = -1
             tab_Bottom( K ) = Lo
             tab_Period( K ) = Hi - Lo
             lo_A( K ) = Lo - tab_Period( K ) * 0.25 + l_Offset( K )
             hi_A( K ) = Lo + tab_Period( K ) * 1.25 + l_Offset( K )
           End If
         End Do
         list_Cyclic( N ) = 1
         tab_Bottom( N ) = Lo
         tab_Period( N ) = Hi - Lo
         lo_A( N ) = Lo - tab_Period( N ) * 0.25
         hi_A( N ) = Lo + tab_Period( N ) * 1.25
         CYC_PAR = 1
       End If

       End Function CYC_PAR



       Subroutine D_SCALE( Factor )

       Real Factor

       If( Factor .le. 0.0 ) Then
         tame_Factor = 1.0
       Else
         tame_Factor = Factor
       End If

       End Subroutine D_SCALE


       Integer Function FIX_PAR( N, n_Par, n_Term )

       ! Used to fix parameter #N throughout a fit
       !                          28 July 1987, KM

       ! Slave parameter
       Integer N
       ! Total number of parameters
       Integer n_Par
       ! Number of free parameters
       Integer n_Term              

       Integer I, K

       ! First, make sure it's not a linked parameter
       If( l_Index( N ) .ne. N ) Then
         Write( *, 100 ) N
100         Format( ' ERROR:: Parameter #', I2, ' is linked!' / &
                          '           This call to FIX_PAR is ignored' )
         FIX_PAR = -1              ! Error status
         Return
       End If
       Do K = 1, n_Par
         If( Abs( l_Index( K ) ) .eq. N .and. K .ne. N ) Then
           Write( *, 100 ) N
           FIX_PAR = -1       ! Error status
           Return
         End If
       End Do
!              Is it already fixed?
       If( Table_Free( N ) .eq. 1 ) Then
         Write( *, 200 ) N
200         Format( ' WARNING:: Parameter #', I2, ' is already fixed!' / &
                          '           This call to FIX_PAR is ignored' )
         FIX_PAR = 0              ! Warning status
         Return
       End If
       Table_Free( N ) = 1
       I = 1
       Do While( I .le. n_Term .and. Index( I ) .ne. N )
         I = I + 1
       End Do
       n_Term = n_Term - 1       ! One less free parameter
       Do K = I, n_Term
         Index( K ) = Index( K + 1 )
       End Do

       FIX_PAR = 1              ! OK status

       End Function FIX_PAR



       Integer Function FREE_PAR( N, n_Par, n_Term )

       ! Parameter N (>0) becomes free again!
       ! Or if N=0, all parameters will become free

       ! Fixed ---> Free
       Integer N
       ! Number of parameters
       Integer n_Par
       ! Number of free parameters
       Integer n_Term       

       Integer K

       If( N .ge. 1 .and. N .le. n_Par ) Then       ! Single mode
         If( l_Index( N ) .ne. N ) Then
           Write( *, 100 ) N
100           Format( ' ERROR:: Parameter #', I2, ' is linked!' /  &
                           '           Use SEPARATE_PAR to separate' )
           FREE_PAR = -1
           Return
         Else If( Table_Free( N ) .eq. 0 ) Then
           Write( *, 200 ) N
200           Format( ' WARNING:: Parameter #', I2, ' is free!' /  &
                   '           This call to FREE_PAR is ignored' )
           FREE_PAR = 0
           Return
         End If
         n_Term = n_Term + 1       ! One more free parameter
         Table_Free( N ) = 0
         Index( n_Term ) = N       ! N must now become free
         FREE_PAR = 1
       Else If( N .eq. 0 ) Then              ! All-clear mode
         Do K = 1, n_Par
           If( Table_Free( K ) .eq. 1 ) Then
             n_Term = n_Term + 1
             Table_Free( K ) = 0
             Index( n_Term ) = K
           End If
         End Do
         FREE_PAR = 1
       Else                                   ! Uninown mode
         Write( *, 300 )
300         Format( ' ERROR:: invalid argument!',  &
                           ' This call to FREE_PAR is ignored' )
         FREE_PAR = -1
       End If

       End Function FREE_PAR



       Subroutine GAUSS_JORDAN( Matrix, Order, Vector, Flag )

!       Solves the linear equation Matrix.Answer=Vector
!       and replaces Vector by Answer, and Matrix by its inverse.
!
!                                          7/July/1987, KM
!                     Modified a bit              August 1 1989, KM


       Integer, intent(in) :: Order                     ! Used size
       Real, dimension(:,:), intent(inout) :: matrix    ! Matrix of interest
       Real, dimension(:), intent(inout) :: vector        ! Vector of interest
       Integer, intent(out) :: Flag                     ! -1 if singular

       Integer Pivot( Order )       ! Bookkeeping
       Integer r_Index( Order )       ! variables for
       Integer c_Index( Order )       ! pivoting

       Real Big, inv_Pivot, Work
       Integer Row, Col, I, J, K, L, M

       Flag = 0
       
       Pivot=0

       ! Main loop over the columns to be reduced
       Do I = 1, Order
         Big = 0.0
         ! Outer loop to search for a pivot element
         Do J = 1, Order
           If( Pivot( J ) .ne. 1 ) Then
             Do K = 1, Order
              If( Pivot( K ) .eq. 0 ) Then
                If( Abs( Matrix( J, K ) ) .ge. Big ) Then
                  Big = Abs( Matrix( J, K ) )
                  Row = J
                  Col = K
                End If
              Else If( Pivot( K ) .gt. 1 ) Then
                Flag = -1
                Return
              End If
             End Do
           End If
         End Do
         Pivot( Col ) = Pivot( Col ) + 1

!         We now have the pivot element, so exchange rows to put the
!         pivot element on the diagonal.  Changes are labeled using
!         *_Index arrays (not actually interchanged).

         If( Row .ne. Col ) Then
           Do L = 1, Order
             Work = Matrix( Row, L )
             Matrix( Row, L ) = Matrix( Col, L )
             Matrix( Col, L ) = Work
           End Do
           Work = Vector( Row )
           Vector( Row ) = Vector( Col )
           Vector( Col ) = Work
         End If
         r_Index( I ) = Row       ! Now ready to divide the pivot row
         c_Index( I ) = Col       ! by the pivot element
         If( Matrix( Col, Col ) <= tiny(matrix(col,col))) Then       
           ! New - protect.  
           ! Used to have a number in there.
           Flag = -1
           Return
         End If
         inv_Pivot = 1.0 / Matrix( Col, Col )
         Matrix( Col, Col ) = 1.0
         Do L = 1, Order
           Matrix( Col, L ) = Matrix( Col, L ) * inv_Pivot
         End Do
         Vector( Col ) = Vector( Col ) * inv_Pivot
                             ! Next, reduce the rows
         Do M = 1, Order
           If( M .ne. Col ) Then       ! ...except for the pivot one
             Work = Matrix( M, Col )
             Matrix( M, Col ) = 0.0
             Do L = 1, Order
              Matrix( M, L ) &
                           = Matrix( M, L ) - Matrix( Col, L ) * Work
             End Do
             Vector( M ) = Vector( M ) - Vector( Col ) * Work
           End If
         End Do
       End Do              ! Endof the main loop over columns
                      ! Now unscramble the solution
       Do L = Order, 1, -1
         If( r_Index( L ) .ne. c_Index( L ) ) Then
           Do K = 1, Order
             Work = Matrix( K, r_Index( L ) )
             Matrix( K, r_Index( L ) ) = Matrix( K, c_Index( L ) )
             Matrix( K, c_Index( L ) ) = Work
           End Do
         End If
       End Do

       End Subroutine GAUSS_JORDAN




       Integer Function LIMIT_PAR( N, Lo, Hi )

       ! Set up the limit within which parameter A( N ) must change
       ! It does NOT check if the current value of A( N ) is within the bound.
       ! The default values are +/-1.0E+36.

       Integer N       ! Parameter of interest
       Real Lo              ! Lower boundary
       Real Hi              ! Upper boundary

       Integer K

       If( Lo .ge. Hi ) Then
         Write( *, 100 )
100         Format( ' ERROR:: Upper limit is smaller than lower limit' / &
                           '           This call to LIMIT_PAR is ignored' )
         LIMIT_PAR = -1
       Else If( N .ge. 1 .and. N .le. m_Par ) Then
         If( list_Cyclic( N ) .ne. 0 ) Then       ! Cyclic
           If( l_Index( N ) .ne. N ) Then
             Write( *, 105 )
105             Format( ' ERROR:: Sorry, you cannot do that' )
             LIMIT_PAR = -1
             Return
           Else If( Lo .ge. tab_Bottom( N ) &
                           .and. Hi .lt. tab_Bottom( N ) + tab_Period( N ) ) Then
             Write( *, 110 )
110             Format( ' WARNING:: LIMITing a cyclic parameter' )
             lo_A( N ) = Lo
             hi_A( N ) = Hi
             list_Cyclic( N ) = -1
             LIMIT_PAR = 0
           Else
             If( Hi - Lo .lt. tab_Period( N ) ) Then
              Write( *, 110 )
              Lo_A( N ) = Lo
              hi_A( N ) = Hi
              list_Cyclic( N ) = -1
               LIMIT_PAR = 0
             Else
              Write( *, 120 )
120              Format( ' ERROR:: Illegal attempt to LIMIT a cyclic parameter' )
              LIMIT_PAR = -1
             End If
           End If
           Do K = 1, N
             If( l_Index( K ) .eq. N .and. K .ne. N ) Then
              lo_A( K ) = lo_A( N ) + l_Offset( K )
              hi_A( K ) = hi_A( N ) + l_Offset( K )
             End If
           End Do
         Else
           Do K = 1, N
             If( l_Index( K ) .eq. -N ) Then
              lo_A( K ) = lo_A( N ) * l_Offset( K )
              hi_A( K ) = hi_A( N ) * l_Offset( K )
             Else If( l_Index( K ) .eq. N .and. K .ne. N ) Then
              lo_A( K ) = lo_A( N ) + l_Offset( K )
              hi_A( K ) = hi_A( N ) + l_Offset( K )
             End If
           End Do
           lo_A( N ) = Lo
           hi_A( N ) = Hi
           LIMIT_PAR = 1
         End If
       Else
         Write( *, 200 ) N
200         Format( ' ERROR:: No such parameter #', I2 / &
                           '           This call to LIMIT_PAR is ignored' )
         LIMIT_PAR = -1
       End If

       End Function LIMIT_PAR




       Integer Function NOLIMIT_PAR( N )

       Integer N       ! Parameter of interest

       Integer K

       If( N .ge. 1 .and. N .le. m_Par ) Then
         If( list_Cyclic( N ) .eq. 1 ) Then
           Write( *, 110 ) N
110           Format( ' ERROR:: No limit has been set for cyclic parameter', I3 )
           NOLIMIT_PAR = -1
         Else If( list_Cyclic( N ) .eq. 0 ) Then
           lo_A( N ) = n_Large
           hi_A( N ) = p_Large
           Do K = 1, m_Par
             If( l_Index( K ) .eq. N .and. K .ne. N ) Then
              lo_A( K ) = lo_A( N ) + l_Offset( K )
              hi_A( K ) = hi_A( N ) + l_Offset( K )
             Else If( l_Index( K ) .eq. -N ) Then
              lo_A( K ) = lo_A( N ) * l_Offset( K )
              hi_A( K ) = hi_A( N ) * l_Offset( K )
             End If
           End Do
           NOLIMIT_PAR = 1
         Else              ! Limited cyclic parameter
           If( l_Index( N ) .eq. N ) Then
             lo_A( N ) = tab_Bottom( N ) - tab_Period( N ) * 0.25
             hi_A( N ) = tab_Bottom( N ) + tab_Period( N ) * 1.25
             Do K = 1, m_Par
              If( l_Index( K ) .eq. N .and. K .ne. N ) Then
                lo_A( K ) = lo_A( N ) + l_Offset( K )
                hi_A( K ) = hi_A( N ) + l_Offset( K )
              End If
             End Do
           Else
             Write( *, 120 )
120             Format( ' ERROR:: Illegal attempt on a cyclic slave parameter' )
             NOLIMIT_PAR = -1
             Return
           End If
           list_Cyclic( N ) = 1
           NOLIMIT_PAR = 1
         End If
       Else If( N .eq. 0 ) Then
         NOLIMIT_PAR = 0
         Do K = 1, m_Par
           If( list_Cyclic( K ) .eq. 1 ) Then
             Write( *, 110 ) K
             NOLIMIT_PAR = -1
           Else If( list_Cyclic( K ) .eq. 0 ) Then
             If( l_Index( K ) .eq. K ) Then
              lo_A( K ) = n_Large
              hi_A( K ) = p_Large
             Else If( l_Index( K ) .gt. 0 ) Then
              lo_A( K ) = lo_A( l_Index( K ) ) + l_Offset( K )
              hi_A( K ) = hi_A( l_Index( K ) ) + l_Offset( K )
             Else
              lo_A( K ) = lo_A( -l_Index( K ) ) * l_Offset( K )
              hi_A( K ) = hi_A( -l_Index( K ) ) * l_Offset( K )
             End If
           Else
             If( l_Index( K ) .eq. K ) Then
              lo_A( K ) = tab_Bottom( K ) - tab_Period( K ) * 0.25
              hi_A( K ) = tab_Bottom( K ) + tab_Period( K ) * 1.25
             Else
              lo_A( K ) = tab_Bottom( l_Index( K ) ) &
                           - tab_Period( l_Index( K ) ) * 0.25 + l_Offset( K )
              hi_A( K ) = tab_Bottom( l_Index( K ) ) &
                           + tab_Period( l_Index( K ) ) * 1.25 + l_Offset( K )
             End If
             list_Cyclic( K ) = 1
           End If
         End Do
       Else
         Write( *, 200 ) N
200         Format( ' ERROR:: No such parameter #', I3 / &
                          '         This call to NOLIMIT_PAR is ignored' )
         NOLIMIT_PAR = -1
       End If

       End Function NOLIMIT_PAR



       Integer Function LINK_PAR( Master, Slave, Cnst, Mode, A, n_Par, n_Term )

       ! Used to `link' parameter #N to parameter #Master such that
       ! A( Slave ) = A( Master ) + Cnst (Mode = 0)
       ! or A( Slave ) = Cnst * A( Master ) (Mode = 1)
       ! at any time during the fit
       ! e.g. fiiting a resolved doublet whose separation is known.
       ! Both Master and Slave must be on the list of free parameters (Index)
       ! prior to the call.  Designed to be called as the final part
       ! of parameter initialization.
       !                                   12 August 1987, KM

       Integer Master              ! Master parameter
       Integer Slave              ! Slave parameter
       Real Cnst              ! Difference or Ratio
       Integer Mode              ! Selects difference or ratio
       Integer n_Par              ! Total number of parameters
       Real A( n_Par )              ! Parameters
       Integer n_Term              ! Number of free parameters

       Integer J, K

       ! First, look for Master and Slave in Index
       LINK_PAR = 1
       If( Master .le. 0 .or. Master .gt. n_Par ) Then
         Write( *, 110 ) Master
110         Format( ' ERROR::', I3, ' is out of range; ', &
                                  'this call to LINK_PAR is ignored' )
         LINK_PAR = -1
         Return
       Else If( Slave .le. 0 .or. Slave .gt. n_Par ) Then
         Write( *, 110 ) Slave
         LINK_PAR = -1
         Return
       Else If( Table_Free( Master ) .eq. 1 ) Then
         Write( *, 120 ) Master
120         Format( ' ERROR:: Parameter #', I3, ' is fixed' / &
                           '         This call to LINK_PAR is ignored' )
         LINK_PAR = -1
         Return
       Else If( Table_Free( Slave ) .eq. 1 ) Then
         Write( *, 120 ) Slave
         LINK_PAR = -1
         Return
       Else If( Table_Free( Master ) .eq. -1 ) Then
         If( Mode .eq. 0 ) Then
           Write( *, 130 ) Master, l_Index( Master )
130           Format( ' WARNING:: Parameter #', I3, ' is linked to', &
                           ' parameter #', I3 )
           Master = l_Index( Master )
           Cnst = Cnst + l_Offset( Master )
           LINK_PAR = 0
         Else
           Write( *, 140 )
140           Format( ' ERROR:: Link mode mismatch' / &
                           '         This call to LINK_PAR is ignored' )
           LINK_PAR = -1
           Return
         End If
       Else If( Table_Free( Master ) .eq. -2 ) Then
         If( Mode .eq. 1 ) Then
           Write( *, 130 ) Master, -l_Index( Master )
           Master = -l_Index( Master )
           Cnst = Cnst * l_Offset( Master )
           LINK_PAR = 0
         Else
           Write( *, 140 )
           LINK_PAR = -1
           Return
         End If
       Else If( Table_Free( Slave ) .lt. 0 ) Then
         Write( *, 150 ) Slave
150         Format( ' ERROR:: Parameter #', I3, ' is already linked' / &
                           '         This call to LINK_PAR is ignored' )
         LINK_PAR = -1
         Return
       Else If( list_Cyclic( Master ) .ne. 0 ) Then
         If( list_Cyclic( Slave ) .eq. 0 ) Then
           Write( *, 160 )
160           Format( ' ERROR:: Illegal attempt to link ', &
                                  'a cyclic parameter with a non-cyclic one' )
           LINK_PAR = -1
           Return
         Else If( list_Cyclic( Slave ) .eq. 1 ) Then
           If( Mode .eq. 1 ) Then
             Write( *, 170 )
170             Format( ' ERROR:: Illegal attempt ', &
                                  'to fix the ratio of cyclic parameters' )
             LINK_PAR = -1
             Return
           End If
           lo_A( Slave ) = lo_A( Master ) + Cnst
           hi_A( Slave ) = hi_A( MAster ) + Cnst
           list_Cyclic( Slave ) = -1
         Else
           Write( *, 180 )
180           Format( ' ERROR:: Sorry, you cannot do that' )
           LINK_PAR = -1
           Return
         End If
       Else If( list_Cyclic( Slave ) .ne. 0 ) Then
         Write( *, 160 )
         LINK_PAR = -1
         Return
       End If

       If( Mode .eq. 0 ) Then
         l_Index( Slave ) = Master
         l_Offset( Slave ) = Cnst
         A( Slave ) = A( Master ) + Cnst
         Table_Free( Slave ) = -1
       Else If( Mode .eq. 1 ) Then
         l_Index( Slave ) = -Master
         l_Offset( Slave ) = Cnst
         A( Slave ) = A( Master ) * Cnst
         Table_Free( Slave ) = -2
       End If
       n_Term = n_Term - 1       ! One less free parameter
       J = 1
       Do While( J .le. n_Term .and. Index( J ) .ne. Slave )
         J = J + 1
       End Do
       Do K = J, n_Term
         Index( K ) = Index( K + 1 )       ! Adjust Index
       End Do

       End Function LINK_PAR




       Integer Function SEPARATE_PAR( Slave, n_Par, n_Term )

       ! Parameter Slave (>0) becomes independent again!
       ! Or if Slave=0, all parameters will become independent

       Integer Slave       ! Slave ---> Free Man
       Integer n_Par       ! Number of parameters
       Integer n_Term       ! Number of free parameters

       Integer K

       SEPARATE_PAR = 1
       If( Slave .ge. 1 .and. Slave .le. n_Par ) Then       ! Single mode
         If( l_Index( Slave ) .eq. Slave ) Then
           Write( *, 100 ) Slave
100           Format( ' ERROR:: Parameter #', I2, ' is not linked!'/ &
                                  '         Call to SEPARATE_PAR ignored' )
           SEPARATE_PAR = -1
           Return
         End If
         l_Index( Slave ) = Slave       ! Now un-linked
         n_Term = n_Term + 1              ! One more free parameter
         Index( n_Term ) = Slave       ! Slave must now become free
         Table_Free( Slave ) = 0
         If( list_Cyclic( Slave ) .eq. -1 ) Then
           lo_A( Slave ) = tab_Bottom( Slave ) - tab_Period( Slave ) * 0.25
           hi_A( Slave ) = tab_Bottom( Slave ) + tab_Period( Slave ) * 1.25
           list_Cyclic( Slave ) = 1
         End If
       Else If( Slave .eq. 0 ) Then              ! All-clear mode
         Do K = 1, n_Par
           If( l_Index( K ) .ne. K ) Then       ! Found a linked par
             l_Index( K ) = K
             n_Term = n_Term + 1
             Index( n_Term ) = K
             Table_Free( K ) = 0
             If( list_Cyclic( K ) .eq. -1 ) Then
              lo_A( K ) = tab_Bottom( K ) - tab_Period( K ) * 0.25
              hi_A( K ) = tab_Bottom( K ) + tab_Period( K ) * 1.25
              list_Cyclic( K ) = 1
             End If
           End If
         End Do
       Else                                   ! Unknown mode
         Write( *, 200 )
200         Format( ' ERROR:: invalid argument!', &
                                  ' This call to SEPARATE_PAR ignored' )
         SEPARATE_PAR = -1
       End If

       End Function SEPARATE_PAR



       Subroutine MARQUARDT( X, Y, W, Points, Lambda, A, dA, &
         Covar, Chi, n_Par, n_Term, FUNCTN, FDERIV, flag_marq )


       ! Number of data points
       Integer, intent(in) :: Points
       ! Data and weights
       Real, dimension(Points), intent(in) :: X, Y, W
       ! Fit control variable
       Real, intent(inout) :: Lambda
       ! # of params of model
       Integer, intent(in) :: n_Par
       ! Parameters and increments
       Real, dimension(n_Par), intent(inout) :: A
       Real, dimension(n_Par), intent(in) :: dA
       ! Covariance matrix
       Real, dimension(:,:), intent(out) :: Covar
       ! Chi squared
       Real, intent(out) :: Chi
       ! # of free paras
       Integer, intent(in) :: n_Term
       ! Model Function and its derivative
       Real, external :: FUNCTN              
       External FDERIV
       ! A flag, normally zero.  Inherited from svd_comp.
       integer, intent(out) :: flag_marq                    

       Real, save :: Alpha( m_Par, m_Par )       ! Curvature matrix
       Real, save :: try_A( m_Par )              ! Trial paras
       Real, save :: Beta( m_Par )
       Real, save :: incr_A( m_Par )
       Real, save :: svd_U( m_Par, m_Par )
       Real, save :: svd_W( m_Par )
       Real, save :: svd_V( m_Par, m_Par )
       Integer, save :: I, J, K, Flag
       Real, save :: o_Chi

       flag_marq=0

       If( Lambda .lt. 0.0 ) Then
         ! First call: Initialize and return
         Lambda = 0.001       ! Initialize Lambda
         Call MRQ_COF( X, Y, W, Points, A, dA, n_Par, n_Term, &
                           Alpha, Beta, Chi, FUNCTN, FDERIV )
         o_Chi = Chi
         Do K = 1, n_Par
           try_A( K ) = A( K )
         End Do
         Return
       End If
       ! Change diagonal elements
       Do I = 1, n_Term
         Do J = 1, n_Term
           Covar( I, J ) = Alpha( I, J )
         End Do
         Covar( I, I ) = Alpha( I, I ) * ( 1.0 + Lambda )
         incr_A( I ) = Beta( I )
       End Do
       Call GAUSS_JORDAN( Covar, n_Term, incr_A, Flag )

       If( Lambda .eq. 0.0 ) Then       ! Final call...error calc only
         If( Flag .ne. 0 ) Then
           Do I = 1, n_Term              ! Recalculate COVAR
             Do J = 1, n_Term
              svd_U( I, J ) = Alpha( I, J )
             End Do
             svd_U( I, I ) = Alpha( I, I ) * ( 1.0 + Lambda )
             incr_A( I ) = Beta( I )
           End Do
           If( Debug ) &
           Print *, '@ Singular Matrix trying Singular Value Decomposition.'
           ! You may like to try svdcmp instead.
           Call SVD_COMP( svd_U, n_Term, n_Term, svd_W, &
           svd_V, flag_marq )
           if (flag_marq /= 0) return
           Call SVD_REC( svd_U, svd_W, svd_V, n_Term, n_Term )
           Do I = 1, n_Term
             Do J = 1, n_Term
              Covar( I, J ) = svd_U( I, J )
             End Do
           End Do
         End If
         Call COV_SORT( Covar, n_Par, n_Term )
         Return
       End If

       If( Flag .ne. 0 ) Then              ! Singular Matrix
         If( Debug ) Then
           Print *, '@ Singular Matrix...'
           Do I = 1, n_Term
             Print *, ( Alpha( I, J ), J = 1, n_Term )
           End Do
           Print *, ' ...trying Singular Value Decomposition'
         End If
         Do I = 1, n_Term               ! Recalculate COVAR
           Do J = 1, n_Term
             svd_U( I, J ) = Alpha( I, J )
           End Do
           svd_U( I, I ) = Alpha( I, I ) * ( 1.0 + Lambda )
           incr_A( I ) = Beta( I )
         End Do
         ! You may like to try svdcmp instead.
         Call SVD_COMP( svd_U, n_Term, n_Term, svd_W, &
         svd_V, flag_marq )
         if (flag_marq /= 0) return
         Call SV_BCSB ( svd_U, svd_W, svd_V, n_Term, n_Term, incr_A )
       End If

       Do I = 1, n_Term
         incr_A( I ) = incr_A( I ) * tame_Factor
       End Do

       Call UPDATE_PARAM( try_A, incr_A, n_Par, n_Term )

       Call MRQ_COF( X, Y, W, Points, try_A, dA, n_Par, n_Term, &
                          Covar, incr_A, Chi, FUNCTN, FDERIV )

       ! This line used to read...
       ! If( Chi .lt. o_Chi ) Then       
       ! But is clearly better done like this.
       if (o_chi - chi > -1.0*epsilon(chi)) then
         ! Success
         Lambda = Lambda * 0.1
         o_Chi = Chi
         Do I = 1, n_Term
           Do J = 1, n_Term
             Alpha( I, J ) = Covar( I, J )
           End Do
           Beta( I ) = incr_A( I )
         End Do
         Do J = 1, n_Par
           A( J ) = try_A( J )
         End Do
       Else                            ! Failure
         Lambda = Lambda * 100.0
         Chi = o_Chi
         Do I = 1, n_Term
           J = Index( I )
           try_A( J ) = A( J )
         End Do
       End If

       End Subroutine MARQUARDT



       Subroutine MRQ_COF( X, Y, W, Points, A, dA, n_Par, n_Term, &
                                 Alpha, Beta, Chi, FUNCTN, FDERIV )

       ! Calculates Alpha, Beta and Chi.
       !                 8 July 1987, KM

       ! Parameters
       Integer, intent(in) :: Points
       Real, dimension(Points), intent(in) :: X, Y, W
       Integer, intent(in) :: n_Par
       Real, dimension(n_par), intent(in) :: A, dA
       Integer, intent(in) :: n_Term
       Real, dimension(:,:), intent(out) :: Alpha 
       Real, dimension(n_Par), intent(out) :: Beta
       Real, intent(out) :: Chi
       Real, external :: FUNCTN
       External FDERIV

       ! Locals
       Real, allocatable, dimension(:) :: dF_dA
       Real :: dY, Work
       Integer :: I, J, K
       
       ! Initialize (symmetric) Alpha, Beta
       Do I = 1, n_Term
         Do J = 1, I
           Alpha( I, J ) = 0.0
         End Do
         Beta( I ) = 0.0
       End Do
       Chi = 0.0

       ! Create array of differentials.
       allocate (df_da(n_par))
       
       ! Summation loop over all data
       Do K = 1, Points
         If( W( K ) .gt. 0.0 ) Then
           dY = Y( K ) - FUNCTN( X( K ), A, n_Par )
           If( cyclic_Data ) Then
             dY = dY - Nint( dY / cycle_Data ) * cycle_Data
           End If
           Call FDERIV( X( K ), A, dA, n_Par, dF_dA )
           Do I = 1, n_Par
             If( Table_Free( I ) .eq. -1 ) Then
              J = l_Index( I )
              dF_dA( J ) = dF_dA( J ) + dF_dA( I )
             Else If( Table_Free( I ) .eq. -2 ) Then
              J = -l_Index( I )
              dF_dA( J ) = dF_dA( J ) + dF_dA( I ) * l_Offset( I )
             End If
           End Do
           Do I = 1, n_Term
             Work = dF_dA( Index( I ) ) * W( K )
             Do J = 1, I
              Alpha( I, J ) = Alpha( I, J ) + Work * dF_dA( Index( J ) )
             End Do
             Beta( I ) = Beta( I ) + dY * Work
           End Do
           Chi = Chi + dY * dY * W( K )
         End If
       End Do
       ! Finished with df_da
       deallocate(df_da)
       
       ! Fill in the symmetric side
       Do J = 2, n_Term
         Do I = 1, J - 1
           Alpha( I, J ) = Alpha( J, I )
         End Do
       End Do

       End Subroutine MRQ_COF



       Subroutine RESET( )

       Integer K

       Do K = 1, m_Par
         Index( K ) = K
         Table_Free( K ) = 0
         l_Index( K ) = K
         l_Offset( K ) = 0.0
         lo_A( K ) = n_Large
         hi_A( K ) = p_Large
         list_Cyclic( K ) = 0
         tab_Period( K ) = 0.0
         tab_Bottom( K ) = 0.0
       End Do
       cyclic_Data = .False.
       cycle_Data = 0.0
       tame_Factor = 1.0

     End subroutine reset



     Subroutine SVD_COMP( Matrix, size_X, size_Y, W, V, flag_svd )

       ! Singular Value Decomposition.

       ! Given an arbitrary matrix A, computes its singular value decomposition
       ! A = U W Vt; the matrix U replaces A on output.  The diagonal matrix
       ! W is output as a vector, and the matrix V (not its transpose Vt) is
       ! calculated.  size_X must be greater or equal to dim_Y; if it is smaller,
       ! then A should be filled up to square with zero rows.

       ! The Numerical Recipes used to look to see if certain numbers were
       ! close to zero, by adding them to a large number, and looking for
       ! equality with that large number.  Fortran 95 allows you to use the
       ! tiny function.

       ! Matrix to be decomposed / U
       real, dimension(:,:), intent(inout) :: Matrix
       ! Actual size of the matrix
       integer, intent(in) :: size_X, size_Y              
       ! Diagonal matrix
       real, dimension(:), intent(out) :: W
       ! V as in A = U.W.Vt
       real, dimension(:,:), intent(out) :: V
       ! A flag, set to -51 if 30 iterations are tried.
       integer, intent(out) :: flag_svd

       Real RV1( m_Par )
       Real Scale
       Real C, G, H, S, X, Y, Z
       double precision :: F
       Integer I, J, K, L, N, Loop
       real :: anorm

       flag_svd=0

       If( size_X .lt. size_Y ) Stop 'ERROR:: unacceptable size'
       If( size_Y .gt. m_Par ) Stop 'ERROR:: matrix too big'

       ! Householder reduction to bidiagonal form

       G = 0.0
       Scale = 0.0
       anorm=0.0
       Do I = 1, size_Y
         L = I + 1
         RV1( I ) = Scale * G
         G = 0.0
         S = 0.0
         Scale = 0.0
         Do K = I, size_X
           Scale = Scale + Abs( Matrix( K, I ) )
         End Do
         If( abs(Scale) > tiny(Scale) ) Then
           Do K = I, size_X
             Matrix( K, I ) = Matrix( K, I ) / Scale
             S = S + Matrix( K, I ) * Matrix( K, I )
           End Do
           F = Matrix( I, I )
           G = -Sign( Sqrt( S ), real(F) )
           H = F * G - S
           Matrix( I, I ) = F - G
           If( I .ne. size_Y ) Then
             Do J = L, size_Y
              S = 0.0
              Do K = I, size_X
                S = S + Matrix( K, I ) * Matrix( K, J )
              End Do
              F = S / H
              Do K = I, size_X
                Matrix( K, J ) = Matrix( K, J ) + F * Matrix( K, I )
              End Do
             End Do
           End If
           Do K = I, size_X
             Matrix( K, I ) = Scale * Matrix( K, I )
           End Do
         End If
         W( I ) = Scale * G
         G = 0.0
         S = 0.0
         Scale = 0.0
         If( I .ne. size_Y ) Then
           Do K = L, size_Y
             Scale = Scale + Abs( Matrix( I, K ) )
          End Do
           If( abs(Scale) > tiny(Scale) ) Then
             Do K = L, size_Y
               Matrix( I, K ) = Matrix( I, K ) / Scale
               S = S + Matrix( I, K ) * Matrix( I, K )
             End Do
             F = Matrix( I, L )
             G = -Sign( Sqrt( S ), real(F) )
             H = F * G - S
             Matrix( I, L ) = F - G
             Do K = L, size_Y
              RV1( K ) = Matrix( I, K ) / H
             End Do
             If( I .ne. size_X ) Then
              Do J = L, size_X
                S = 0.0
                Do K = L, size_Y
                  S = S + Matrix( J, K ) * Matrix( I, K )
                End Do
                Do K = L, size_Y
                  Matrix( J, K ) = Matrix( J, K ) + S * RV1( K )
                End Do
              End Do
             End If
             Do K = L, size_Y
               Matrix( I, K ) = Scale * Matrix( I, K )
             End Do
           End If
         End If
         anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
       End Do

       ! Accumulation of right-hand transformations
       Do I = size_Y, 1, -1
         If( I .lt. size_Y ) Then
           if (abs(G) > tiny(G))then
             Do J = L, size_Y       
               ! Double division to avoid possible underflow
               V( J, I ) = ( Matrix( I, J ) / Matrix( I, L ) ) / G
             End Do
             Do J = L, size_Y
               S = 0.0
               Do K = L, size_Y
                 S = S + Matrix( I, K ) * V( K, J )
               End Do
               Do K = L, size_Y
                 V( K, J ) = V( K, J ) + S * V( K, I )
               End Do
             End Do
           End If
           Do J = L, size_Y
             V( I, J ) = 0.0
             V( J, I ) = 0.0
           End Do
         End If
         V( I, I ) = 1.0
         G = RV1( I )
         L = I
       End Do

       ! Accumulation of left-hand transformations

       Do I = size_Y, 1, -1
         L = I + 1
         G = W( I )
         If( I .lt. size_Y ) Then
           Do J = L, size_Y
             Matrix( I, J ) = 0.0
           End Do
         End If
         if (abs(G) > tiny(G))then
           G = 1.0 / G
           If( I .ne. size_Y ) Then
              Do J = L, size_Y
                 S = 0.0
              Do K = L, size_X
                S = S + Matrix( K, I ) * Matrix( K, J )
              End Do
              F = ( S / Matrix( I, I ) ) * G
              Do K = I, size_X
                Matrix( K, J ) = Matrix( K, J ) + F * Matrix( K, I )
              End Do
             End Do
           End If
           Do J = I, size_X
             Matrix( J, I ) = Matrix( J, I ) * G
           End Do
         Else
           Do J = I, size_X
             Matrix( J, I ) = 0.0
           End Do
         End If
         Matrix( I, I ) = Matrix( I, I ) + 1.0
       End Do

       ! Diagonalization of the bidiagonal form
       Do K = size_Y, 1, -1
         Do Loop = 1, 30
           Do L = K, 1, -1
             N = L - 1
             if (abs(rv1(l)) < epsilon(anorm))  goto 200
             if (abs(w(n))   < epsilon(anorm))  goto 100
           End Do
100          Continue
           C = 0.0
           S = 1.0
           Do I = L, K
             F = S * RV1( I )
             if (abs(f) > epsilon(anorm)) then
               G = W( I )
               H = Sqrt( F * F + G * G )
               W( I ) = H
               H = 1.0 / H
               C = G * H
               S = -F * H
               Do J = 1, size_X
                 Y = Matrix( J, N )
                 Z = Matrix( J, I )
                 Matrix( J, N ) = Y * C + Z * S
                 Matrix( J, I ) = -Y * S + Z * C
               End Do
             End If
           End Do
200          Continue
           Z = W( K )
           If( L .eq. K ) Then              
             ! Convergence
             If( Z .lt. 0.0 ) Then       
               ! Singular value is made nonnegative
               W( K ) = -Z
               Do J = 1, size_Y
                 V( J, K ) = -V( J, K )
               End Do
             End If
             Goto 300
           End If
           If ( Loop == 30 ) then
             if (debug) print*, '@ ERROR:: No convergence in 30 iterations'
             !flag_svd = -51
             !return
           end if
           ! Shift from bottom 2-by-2 minor
           X = W( L )               
           N = K - 1
           Y = W( N )
           G = RV1( N )
           H = RV1( K )
           F = ( ( Y - Z ) * ( Y + Z ) + ( G - H ) * ( G + H ) ) &
                                                      / ( 2.0 * H * Y )
           if (abs(F) < Sqrt(huge(G)) ) then 
             G = Sqrt( F * F + 1.0 )
           else 
             G=F
           End if 
           F = ( ( X - Z ) * ( X + Z ) &
                          + H * ( Y / ( F + Sign( G, real(F) ) ) - H ) ) / X
           ! Next QR transformation
           C = 1.0
           S = 1.0
           Do J = L, N
             I = J + 1
             G = RV1( I )
             Y = W( I )
             H = S * G
             G = C * G
             Z = Sqrt( F * F + H * H )
             RV1( J ) = Z
             C = F / Z
             S = H / Z
             F = X * C + G * S
             G = -X * S + G * C
             H = Y * S
             Y = Y * C
             Do N = 1, size_Y
              X = V( N, J )
              Z = V( N, I )
              V( N, J ) = X * C + Z * S
              V( N, I ) = -X * S + Z * C
             End Do
             Z = Sqrt( F * F + H * H )
             W( J ) = Z
             ! Rotation can be arbitrary if Z=0
             if (abs(Z) > tiny(Z)) then
               Z = 1.0 / Z
               C = F * Z
               S = H * Z
             End If
             F = C * G + S * Y
             X = -S * G + C * Y
             Do N = 1, size_X
              Y = Matrix( N, J )
              Z = Matrix( N, I )
              Matrix( N, J ) = Y * C + Z * S
              Matrix( N, I ) = -Y * S + Z * C
             End Do
           End Do
           RV1( L ) = 0.0
           RV1( K ) = F
           W( K ) = X
         End Do
300      Continue
       End Do

       End Subroutine SVD_COMP



       Subroutine SV_BCSB( U, W, V, size_X, size_Y, B )

       ! U, W, V as returned by SVD_COMP, B is A.X on input, X on output.
       real, dimension(:,:) :: U, V
       real, dimension(:) :: W, B
       Integer size_X, size_Y       ! Actual sizes used

       Real Sum, Work( m_Par ), max_W, min_W
       Integer I, J, K
!                            ZERO the small W's
       max_W = W( 1 )
       Do J = 2, size_Y
         max_W = Max( max_W, W( J ) )
       End Do
       min_W = max_W * 1.0E-06
       Do J = 1, size_Y
         If( W( J ) .lt. min_W ) W( J ) = 0.0
       End Do

       Do J = 1, size_Y
         Sum = 0.0
         If( W( J ) .ne. 0.0 ) Then
           Do I = 1, size_X
             Sum = Sum + U( I, J ) * B( I )
           End Do
           Sum = Sum /  W( J )
         End If
         Work( J ) = Sum
       End Do

       Do J = 1, size_Y
         Sum = 0.0
         Do K = 1, size_Y
           Sum = Sum + V( J, K ) * Work( K )
         End Do
         B( J ) = Sum
       End Do

       End Subroutine SV_BCSB



      Subroutine SVD_REC( U, W, V, size_X, size_Y )

       ! U, W, V as returned by SVD_COMP
       ! U is replaced by the product Ut.W-1.Vt
       real, dimension(:,:) :: U, V       
       real, dimension(:) :: W
       Integer size_X, size_Y       ! Actual sizes used
       Real Work( m_Par, m_Par ), max_W, min_W
       Integer I, J, K
       
       max_W = maxval(W(1:size_Y))
       ! Handle the small W's
       min_W = max(max_W * 1.0E-06,tiny(max_W))
       Do J = 1, size_Y
         W( J ) = 1.0 / Max( W( J ), min_W )
       End Do

       Do I = 1, size_X
         Do J = 1, size_Y
           Work( I, J ) = 0.0
           Do K = 1, size_Y
             Work( I, J ) = Work( I, J ) + V( J, K ) * W( K ) * U( I, K )
           End Do
         End Do
       End Do

       Do I = 1, size_X
         Do J = 1, size_Y
           U( I, J ) = Work( I, J )
         End Do
       End Do

!       now answer(1...dimy)=sum(u(i,j)*question(i)),i=1...dimx

       End Subroutine SVD_REC



       Subroutine UPDATE_PARAM( try_A, incr_A, n_Par, n_Term )

       ! Update trial parameters within FIT2D. Checks for boundaries, and
       ! scales incr_A down if necessary to keep try_A within a sensible range.
       !                                   13 July 1987, KM
       !                            revised 23 Dec  1987

       ! ARGUMENTS
       Integer n_Par              ! Number of parameters
       Integer n_Term              ! Number of free parameters
       Real try_A( n_Par )       ! The current parameter values
       Real incr_A( n_Par )       ! To be changed by...

       Real Ratio, m_Ratio, Scale
       Real Labour
       Integer J, K, L
       Integer r_Index( m_Par )

       ! Is any parameter going off-boundary?
       m_Ratio = 0.0
       Do K = 1, n_Term
         ! Un-linked parameters first
         J = Index( K )
         r_Index( J ) = K
         If( incr_A( K ) .gt. 0.0 ) Then
           If( hi_A( J ) .le. try_A( J ) ) Then
             If( Debug ) Then
              Write( *, 99 ) J
 99              Format( ' @ Parameter #', I2, ' was about to go off-limit' / &
                           '    NOT UPDATING PARAMETERS DURING THIS ITERATION' )
             End If
             Return
           End If
           Ratio = incr_A( K ) / ( hi_A( J ) - try_A( J ) )
         Else If( incr_A( K ) .eq. 0.0 ) Then
           Ratio = 0.0
         Else
           If( lo_A( J ) .ge. try_A( J ) ) Then
             If( Debug ) Then
              Write( *, 99 ) J
             End If
             Return
           End If
           Ratio = incr_A( K ) / ( lo_A( J ) - try_A( J ) )
         End If
         If( Ratio .ge. 1.0 ) Then
           m_Ratio = Max( m_Ratio, Ratio )
           If( Debug ) Then
             Write( *, 100 ) J
100             Format( ' @ Parameter #', I2, ' was about to go off-limit' )
           End If
         End If
       End Do
       Do K = 1, n_Par                            ! Linked parameters next
         If( Table_Free( K ) .lt. 0 ) Then       ! Linked
           If( Table_Free( K ) .eq. -1 ) Then       ! Mode = 0
             Labour = incr_A( r_Index( l_Index( K ) ) )
           Else                            ! Mode = 1
             Labour = incr_A( r_Index( -l_Index( K ) ) ) * l_Offset( K )
           End If
           If( Labour .gt. 0.0 ) Then
             Ratio = Labour / ( hi_A( K ) - try_A( K ) )
           Else If( Labour .eq. 0.0 ) Then
             Ratio = 0.0
           Else
             Ratio = Labour / ( lo_A( K ) - try_A( K ) )
           End If
           If( Ratio .ge. 1.0 ) Then
             If( Debug ) Then
              Write( *, 100 ) K
             End If
             m_Ratio = Max( m_Ratio, Ratio )
           End If
         End If
       End Do
                             ! Rescale incr_A if necessary
       If( m_Ratio .ge. 1.0 ) Then
         Scale = 0.95 / m_Ratio
         If( Debug ) Then
           Write( *, 110 ) Scale
110           Format( ' @ Changes in trial parameters are scaled down by ', &
                                                              1PE10.3 )
         End If
       Else
         Scale = 1.0
       End If
                             ! Now calculate the new trial parameters
       Do K = 1, n_Term
         J = Index( K )
         try_A( J ) = try_A( J ) + incr_A( K ) * Scale
       End Do
                      ! Linked changes
       Do K = 1, n_Par
         If( Table_Free( K ) .eq. -1 ) Then
           L = l_Index( K )
           try_A( K ) = try_A( L ) + l_Offset( K )
         Else  If( Table_Free( K ) .eq. -2 ) Then
           L = -l_Index( K )
           try_A( K ) = try_A( L ) * l_Offset( K )
         End If
       End Do
                      ! Actively cyclic
       Do K = 1, n_Par
         If( list_Cyclic( K ) .eq. 1 ) Then
           If( try_A( K ) .lt. tab_Bottom( K ) ) Then
             try_A( K ) = tab_Bottom( K ) + tab_Period( K ) &
                           + Mod( try_A( K ) - tab_Bottom( K ), tab_Period( K ) )
           Else
             try_A( K ) = tab_Bottom( K ) &
                           + Mod( try_A( K ) - tab_Bottom( K ), tab_Period( K ) )
           End If
         End If
       End Do

       End Subroutine UPDATE_PARAM


       ! Here is a lightly modified version of the Numerical Recipes 
       ! svdcmp, so you can call it directly instead of the ARK one.
       ! A futher improvement would be to change the statements with
       ! anorm in them to use the epsilon function. 

      SUBROUTINE svdcmp(a, m, n, w, v, flag_svd)

      INTEGER m,n
      real, dimension(:,:) :: a, v
      real, dimension(:) :: w
      ! A flag, set to -51 if 30 iterations are tried.
      integer, intent(out) :: flag_svd

      INTEGER i,its,j,jj,k,l,nm
      REAL anorm,c,f,g,h,s,scale,x,y,z,rv1(m_par)

      flag_svd=0

      g=0.0
      scale=0.0
      anorm=0.0
      do 25 i=1,n
        l=i+1
        rv1(i)=scale*g
        g=0.0
        s=0.0
        scale=0.0
        if(i.le.m)then
          do 11 k=i,m
            scale=scale+abs(a(k,i))
11        continue
          if(abs(scale) > tiny(scale))then
            do 12 k=i,m
              a(k,i)=a(k,i)/scale
              s=s+a(k,i)*a(k,i)
12          continue
            f=a(i,i)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,i)=f-g
            do 15 j=l,n
              s=0.0
              do 13 k=i,m
                s=s+a(k,i)*a(k,j)
13            continue
              f=s/h
              do 14 k=i,m
                a(k,j)=a(k,j)+f*a(k,i)
14            continue
15          continue
            do 16 k=i,m
              a(k,i)=scale*a(k,i)
16          continue
          endif
        endif
        w(i)=scale *g
        g=0.0
        s=0.0
        scale=0.0
        if((i.le.m).and.(i.ne.n))then
          do 17 k=l,n
            scale=scale+abs(a(i,k))
17        continue
          if(abs(scale) > tiny(scale))then
            do 18 k=l,n
              a(i,k)=a(i,k)/scale
              s=s+a(i,k)*a(i,k)
18          continue
            f=a(i,l)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,l)=f-g
            do 19 k=l,n
              rv1(k)=a(i,k)/h
19          continue
            do 23 j=l,m
              s=0.0
              do 21 k=l,n
                s=s+a(j,k)*a(i,k)
21            continue
              do 22 k=l,n
                a(j,k)=a(j,k)+s*rv1(k)
22            continue
23          continue
            do 24 k=l,n
              a(i,k)=scale*a(i,k)
24          continue
          endif
        endif
        anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
25    continue
      do 32 i=n,1,-1
        if(i.lt.n)then
          if(abs(g) > tiny(g))then
            do 26 j=l,n
              v(j,i)=(a(i,j)/a(i,l))/g
26          continue
            do 29 j=l,n
              s=0.0
              do 27 k=l,n
                s=s+a(i,k)*v(k,j)
27            continue
              do 28 k=l,n
                v(k,j)=v(k,j)+s*v(k,i)
28            continue
29          continue
          endif
          do 31 j=l,n
            v(i,j)=0.0
            v(j,i)=0.0
31        continue
        endif
        v(i,i)=1.0
        g=rv1(i)
        l=i
32    continue
      do 39 i=min(m,n),1,-1
        l=i+1
        g=w(i)
        do 33 j=l,n
          a(i,j)=0.0
33      continue
        if(abs(g) > tiny(g))then
          g=1.0/g
          do 36 j=l,n
            s=0.0
            do 34 k=l,m
              s=s+a(k,i)*a(k,j)
34          continue
            f=(s/a(i,i))*g
            do 35 k=i,m
              a(k,j)=a(k,j)+f*a(k,i)
35          continue
36        continue
          do 37 j=i,m
            a(j,i)=a(j,i)*g
37        continue
        else
          do 38 j= i,m
            a(j,i)=0.0
38        continue
        endif
        a(i,i)=a(i,i)+1.0
39    continue
      do 49 k=n,1,-1
        do 48 its=1,30
          do 41 l=k,1,-1
            nm=l-1
            if((abs(rv1(l))+anorm).eq.anorm)  goto 2
            if((abs(w(nm))+anorm).eq.anorm)  goto 1
41        continue
1         c=0.0
          s=1.0
          do 43 i=l,k
            f=s*rv1(i)
            rv1(i)=c*rv1(i)
            if((abs(f)+anorm).eq.anorm) goto 2
            g=w(i)
            h=pythag(f,g)
            w(i)=h
            h=1.0/h
            c= (g*h)
            s=-(f*h)
            do 42 j=1,m
              y=a(j,nm)
              z=a(j,i)
              a(j,nm)=(y*c)+(z*s)
              a(j,i)=-(y*s)+(z*c)
42          continue
43        continue
2         z=w(k)
          if(l.eq.k)then
            if(z.lt.0.0)then
              w(k)=-z
              do 44 j=1,n
                v(j,k)=-v(j,k)
44            continue
            endif
            goto 3
          endif
          if (its == 30) then
            if (debug) print*, '@ ERROR:: No convergence in 30 iterations'
            flag_svd = -51
            return
          end if
          x=w(l)
          nm=k-1
          y=w(nm)
          g=rv1(nm)
          h=rv1(k)
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
          g=pythag(f,1.0)
          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
          c=1.0
          s=1.0
          do 47 j=l,nm
            i=j+1
            g=rv1(i)
            y=w(i)
            h=s*g
            g=c*g
            z=pythag(f,h)
            rv1(j)=z
            c=f/z
            s=h/z
            f= (x*c)+(g*s)
            g=-(x*s)+(g*c)
            h=y*s
            y=y*c
            do 45 jj=1,n
              x=v(jj,j)
              z=v(jj,i)
              v(jj,j)= (x*c)+(z*s)
              v(jj,i)=-(x*s)+(z*c)
45          continue
            z=pythag(f,h)
            w(j)=z
            if(abs(z) > tiny(z))then
              z=1.0/z
              c=f*z
              s=h*z
            endif
            f= (c*g)+(s*y)
            x=-(s*g)+(c*y)
            do 46 jj=1,m
              y=a(jj,j)
              z=a(jj,i)
              a(jj,j)= (y*c)+(z*s)
              a(jj,i)=-(y*s)+(z*c)
46          continue
47        continue
          rv1(l)=0.0
          rv1(k)=f
          w(k)=x
48      continue
3       continue
49    continue
      return
      END subroutine svdcmp

      real FUNCTION pythag(a,b)

      REAL a,b
      REAL absa,absb
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
        pythag=absa*sqrt(1.+(absb/absa)**2)
      else
        if(absb.eq.0.)then
          pythag=0.
        else
          pythag=absb*sqrt(1.+(absa/absb)**2)
        endif
      endif
      return
      END function pythag

  end module marq


