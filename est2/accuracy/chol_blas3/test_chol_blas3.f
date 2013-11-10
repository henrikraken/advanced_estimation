      program main

*     Test cholesky that uses level 3 BLAS on large problem

      implicit none


* Another concern is writing your programs so they can be carried from one compiler to the next without the need for
* reprogramming compiler dependent statements.  In choosing kind numbers there are ways to make your codes "portable"
* (meaning they can be executed on any Fortran 90 compiler without compiler dependent errors).  You can do this by
* declaring a kind number for single and double-precision in the declaration section of your program.
* 
*     INTEGER, PARAMETER  ::  double = SELECTED_REAL_KIND (12)
* 
*     REAL (KIND = double) :: A, B
* 
* In this series of commands the parameter double has been assigned the kind value for double-precision (p = 12 will
* require double precision).  The kind for A and B has been assigned the number for double-precision.  Code written with
* such declarations are portable without corrections. 

      integer, parameter  ::  double = selected_real_kind (12)



      integer n_dim, m_dim

      parameter (n_dim=5000)

      parameter (m_dim=16000)

      integer m, n, i, j, nprint_m, nprint_n, nprint, nblock, mblock
      integer obs_remaining, obs_this_time
      integer i_obs

      real (kind = double) :: h_y( m_dim, n_dim+1 ), xhat(n_dim)
*     real (kind = double) :: xhat(n_dim)

      real (kind = double) :: htwh( n_dim, n_dim ), htwy(n_dim)
      real (kind = double) :: diff, max_diff, sum, rms
      real (kind = double) :: mag_xhat, mag_deltaxhat, rel_err
      real (kind = double) :: start1, end1, time1
      real (kind = double) :: start2, end2, time2
      real (kind = double) :: end3, time3
      real (kind = double) :: Lambda_max, Lambda_min, C

      real (kind = double) :: secondr
      real (kind = double) :: start1b, end1b, time1b
      real (kind = double) :: start2b, end2b, time2b
      real (kind = double) :: ops, speed
      real (kind = double) :: y, h_row( n_dim)
      real (kind = double) :: dot
      real (kind = double) :: seed, n_real, m_real
      real (kind = double) :: epsilon
      real (kind = double) :: x_i, mag_x
      real (kind = double) :: max_abs_rel_err

      integer count_rate, count_max
      integer count_start_1, count_end_1
      integer count_start_2, count_end_2
      integer count_start_3, count_end_3
      integer max_abs_rel_err_i




      print *, ' '
      print *, 'Cholesky BLAS3'

*.....Set m, n

      print *, 'Enter value for m (no max )'
      read *, m
      print *, 'Enter value for n (max of ', n_dim, ')'
      read *, n
      print *, 'Enter nblock size'
      read *, nblock
      print *, 'Enter mblock size'
      read *, mblock
      print *, 'Enter epsilon'
      read *, epsilon
      print *, 'm= ', m, '  n= ', n, '  nblock= ', nblock,
     >         '  mblock= ', mblock, '  epsilon= ', epsilon

      if ( n .gt. n_dim ) then
        print *, 'ERROR:  n= ', n, ' > n_dim= ', n_dim,
     >           ' in test_chol_blas3'
        call abort2
      endif

      if ( mblock .gt. m_dim ) then
        print *, 'mblock= ', mblock, ' > m_dim= ', m_dim,
     >           ' in test_chol_blas3'
        call abort2
      endif

*.....Set values that tell how many rows and columns at the
*     top of matrices to print.

      nprint = 0
      nprint = 7
      nprint_m = min0( nprint, m )
      nprint_n = min0( nprint, n )
*     nprint_m = 10
*     nprint_n = 10

      call system_clock(count_start_1, count_rate, count_max)

*.....Zero out HTWH, HTWY.

      do j = 1, n
        htwy( j ) = 0.d0
        do i = 1, n
          htwh( i, j ) = 0.d0
        end do
      end do

*.....Loop to read in obs.
*     Read in mblock obs. at one time.

      i_obs = 0
      obs_remaining = m

      do while ( obs_remaining > 0 )

        obs_this_time = min0( mblock, obs_remaining )
        obs_remaining = obs_remaining - obs_this_time

*       print *, 'Obs. added to normal eqn. = ', obs_this_time
        do i = 1, obs_this_time
          i_obs = i_obs + 1
          call set_h_y_i( h_row, y, n, i_obs, seed, epsilon )
*         print 1003, i_obs, ( h_row( j ), j = 1, nprint_n ), y
          call dcopy( n, h_row, 1, h_y( i, 1 ), m_dim )
          h_y( i, n+1 ) = y
        end do

*       print *, ' '
*       print *, 'h_y'
*       do i = 1, nprint_m
*         print 1003, i, ( h_y( i, j ), j = 1, nprint_n ), h_y(i, n+1)
*       end do

*.......Add these obs into htwh.
*       Only need lower triangular part.

        call dsyrk( 'Lower', 'Transpose', n, obs_this_time,
     >               1.d0, h_y, m_dim, 1.d0, htwh, n_dim )

*.......Add these obs into htwy.

        call dgemv( 'Transpose', obs_this_time, n, 1.d0, h_y, m_dim,
     >               h_y( 1, n+1 ), 1, 1.d0, htwy, 1 )

*       do i=1,n
*         dot = 0.d0
*         print *, 'Taking the dot product of this row of h_y_T'
*         print 1004, i, (h_y(j,i), j=1,obs_this_time)
*         print *, 'with the column vector y, here printed as a row'
*         print 1005, (h_y(j,n+1), j=1,obs_this_time)
*         do j=1,obs_this_time
*           dot = dot + h_y(j,i) * h_y(j,n+1)
*         end do
*         htwy(i) = htwy(i) + dot
*       end do

      end do

      call system_clock(count_end_1, count_rate, count_max)
      time1 = dble(count_end_1 - count_start_1) / dble(count_rate)

*.....Calculate speed ( Mflops )

      n_real = dble(n)
      m_real = dble(m)
      ops = m_real*n_real*n_real
      speed = ( ops / time1 ) / 1.d+6

      print 1400, ops, time1, speed

*.....Print by rows.

      print *, ' '
      print *, 'htwh'
      do i = 1, nprint_n
        print 1002, i, ( htwh( i, j ), j = 1, i )
      end do

      print *, ' '
      print *, 'htwy'
      do i = 1, nprint_n
        print 1012, htwy( i )
      end do

*.....Cholesky

      call system_clock(count_start_2, count_rate, count_max)

      call chol_blas3( htwh, n_dim, n, nblock )

      call system_clock(count_end_2, count_rate, count_max)
      time2 = dble(count_end_2 - count_start_2) / dble(count_rate)

*.....Calculate speed ( Mflops )

      n_real = dble(n)
      ops = n_real*n_real*n_real / 3.d0
      speed = ( ops / time2 ) / 1.d+6

      print 1410, ops, time2, speed

*.....Print by rows.

      print *, ' '
      print *, 'L'
      do i = 1, nprint_n
        print 1002, i, ( htwh( i, j ), j = 1, i )
      end do

*.....Solve for xhat.

      call dcopy( n, htwy, 1, xhat, 1 )

      call dtrsv( 'Lower', 'No trans', 'Non-unit',
     >             n, htwh, n_dim, xhat, 1 )

      call dtrsv( 'Lower', 'Trans', 'Non-unit',
     >             n, htwh, n_dim, xhat, 1 )

      call system_clock(count_end_3, count_rate, count_max)
      time3 = dble(count_end_3 - count_start_1) / dble(count_rate)

      print 6000, time3, time2, speed

* *.....Check elements of xhat.  RMS of error in xhat.
* 
*       rms = 0.d0
*       do i = 1, n
*         x_i = dble(i)
*         diff = xhat( i ) - x_i
*         rms = rms + diff*diff
*       end do
*       rms = dsqrt( rms / dble( n ) )
* 
*       print 1200, epsilon, rms
* 
* *.....Check elements of xhat.  Ratio of mag of error in xhat to mag of xhat.
* 
*       mag_xhat = 0.d0
*       mag_deltaxhat = 0.d0
*       do i = 1, n
*         x_i = dble(i)
*         diff = xhat( i ) - x_i
*         mag_deltaxhat = mag_deltaxhat + diff*diff
*         mag_xhat = mag_xhat + xhat( i ) * xhat( i )
*       end do
*       mag_xhat= dsqrt(mag_xhat)
*       mag_deltaxhat= dsqrt(mag_deltaxhat)
*       rel_err = mag_deltaxhat / mag_xhat
* 
*       print 1201, epsilon, rel_err
* 
* *.....Check elements of xhat.  RMS of (error in xhat(i)) / xhat(i).
* 
*       rms = 0.d0
*       do i = 1, n
*         x_i = dble(i)
*         diff = ( xhat( i ) - x_i ) / xhat( i )
*         rms = rms + diff*diff
*       end do
*       rms = dsqrt( rms / dble( n ) )
* 
*       print 1202, epsilon, rms
* 
* 
* *.....Check elements of xhat.  Ratio of mag of error in xhat to mag of x.
* 
*       mag_x = 0.d0
*       mag_deltaxhat = 0.d0
*       do i = 1, n
*         x_i = dble(i)
*         diff = xhat( i ) - x_i
*         mag_deltaxhat = mag_deltaxhat + diff*diff
*         mag_x = mag_x + x_i * x_i
*       end do
*       mag_x= dsqrt(mag_x)
*       mag_deltaxhat= dsqrt(mag_deltaxhat)
*       rel_err = mag_deltaxhat / mag_x
* 
*       print 1203, epsilon, rel_err

* *.....Check elements of xhat.  RMS of (error in xhat(i)) / x(i).
* 
*       rms = 0.d0
*       do i = 1, n
*         x_i = dble(i)
*         diff = ( xhat( i ) - x_i ) / x_i
*         rms = rms + diff*diff
*       end do
*       rms = dsqrt( rms / dble( n ) )
* 
*       print 1204, epsilon, rms


*.....Check elements of xhat.  Max of abs((error in xhat(i)) / x(i))

      max_abs_rel_err = 0.d0
      do i = 1, n
        x_i = dble(i)
*       diff = ( xhat( i ) - x_i ) / x_i
*       print 1206, diff
        diff = dabs( ( xhat( i ) - x_i ) / x_i )
        if ( diff .gt. max_abs_rel_err ) then
           max_abs_rel_err=diff
           max_abs_rel_err_i=i
        endif
      end do

      print 1205, epsilon, max_abs_rel_err, max_abs_rel_err_i







*-* *.....Calculate covariance
*-* 
*-* *.....Set identity
*-* 
*-*       call set_to_ident ( h_y, n_dim, n )
*-* 
*-* *.....First Trsm
*-* 
*-*       print *, ' '
*-*       print *, 'First trsm call'
*-*       call dtrsm( 'Left', 'Lower', 'No trans', 'Non-unit',
*-*      >             n, n, 1.d0, htwh, n_dim, h_y, n_dim )
*-* 
*-* *.....Second Trsm
*-* 
*-*       print *, ' '
*-*       print *, 'Second trsm call'
*-*       call dtrsm( 'Left', 'Lower', 'Trans', 'Non-unit',
*-*      >             n, n, 1.d0, htwh, n_dim, h_y, n_dim )
*-* 
*-*       call cpu_time(end1)
*-*       end1b  = secondr( )
*-*       time1 = end1 - start1
*-*       time1b = end1b - start1b
*-*       print 1100, time2, time1
*-*       print 1101, time2b, time1b
*-* 
*-* *.....Calculate sum of elements of P
*-* 
*-*       call sum_matrix( h_y, n, n, n_dim, sum )
*-* 
*-*       print *, ' '
*-*       print *, 'checksum of P = ', sum
*-* 
*-*       if ( nprint .gt. 0 ) then
*-* 
*-*         print *, ' '
*-*         print *, 'P'
*-*         do i = 1, nprint_n
*-*           print 1002, i, ( h_y( i, j ), j = 1, nprint_n )
*-*         end do
*-* 
*-*       endif
*-* 
*-* *.....Calculate correlation matrix, max sigma, max corr. coef.
*-* 
*-*       call corr ( h_y, n, n_dim )
*-* 
*-*       if ( nprint .gt. 0 ) then
*-* 
*-*         print *, ' '
*-*         print *, 'Correlation Matrix'
*-*         do i = 1, nprint_n
*-*           print 1002, i, ( h_y( i, j ), j = 1, i )
*-*         end do
*-* 
*-*       endif


 1000 format ( 'Row ', i2, ':', 1p10e12.4, / (6x, 1p10e12.4)  )
 1001 format ( 'Row ', i2, ':', 10f8.4, / (6x,10e8.4)  )
 1002 format ( 'Row ', i2, ':', 1p8e9.1, / (7x, 1p9e9.1)  )
 1003 format ( 'HRow ', i2, ':', 1p11e9.1, / (7x, 1p11e9.1)  )
 1004 format ( 'Row ', i2, ':', 1p11e9.1, / (7x, 1p11e9.1)  )
 1005 format ( 1p11e9.1, / (7x, 1p11e9.1)  )
 1010 format ( 1pe12.4 )
 1011 format ( f8.4 )
 1012 format ( 1pe9.1 )
 1100 format ( / 1x, 'Cholesky BLAS3 time =', 1p1e11.3,
     >           1x, 'Total time =',  1p1e11.3 / )
 1101 format ( / 1x, 'Cholesky BLAS3 time =', 1p1e11.3,
     >           1x, 'Total time =',  1p1e11.3 / )
 1200 format (//' Epsilon= ', 1pe10.2, 3x,
     >          ' RMS of xhat(i) - x(i) =            ', 1pe10.2   )
 1201 format (  ' Epsilon= ', 1pe10.2, 3x,
     >          ' |xhat-x| / |xhat| =                ', 1pe10.2   )
 1202 format (  ' Epsilon= ', 1pe10.2, 3x,
     >          ' RMS of ( xhat(i)-x(i) / xhat(i) ) =', 1pe10.2   )
 1203 format (  ' Epsilon= ', 1pe10.2, 3x,
     >          ' |xhat-x| / |x| =                   ', 1pe10.2   )
 1204 format (  ' Epsilon= ', 1pe10.2, 3x,
     >          ' RMS of ( xhat(i)-x(i) / x(i) ) =   ', 1pe10.2 / )
 1205 format (  ' Epsilon= ', 1pe10.2, 3x,
     >          ' Max of abs(( xhat(i)-x(i) ) / x(i)) = ', 1pe10.2,
     >          ' at i = ', i5 / )
 1206 format(   ' DIFF ',  1pe10.2 )
 1400 format (// 1x, 'Accumulation        :',
     >           1x, 'Ops = ',   1pe9.3,
     >           2x, 'Time = ',  1pe9.3, ' sec',
     >           2x, 'Speed = ', 1pe9.3, ' Mf'    )
 1410 format (// 1x, 'Cholesky BLAS3      :',
     >           1x, 'Ops = ',   1pe9.3,
     >           2x, 'Time = ',  1pe9.3, ' sec',
     >           2x, 'Speed = ', 1pe9.3, ' Mf'    )
 2000 format ( / a80 )
 3000 format ( / 1x, 'Cholesky BLAS3: time=', 1p1e10.3,
     .           3x, 'Speed (Mflops)=', 1p1e10.3 //  )
 4000 format ( / 1x, 'Accumulation: time=', 1p1e10.3,
     .           3x, 'Speed (Mflops)=', 1p1e10.3 //  )
 6000 format (// 1x, 'Total time=',     1p1e10.3,
     .           3x, 'Cholesky BLAS3: time=', 1p1e10.3,
     .           2x, 'Speed (Mflops)=',     1p1e10.3 // )


      stop
      end
