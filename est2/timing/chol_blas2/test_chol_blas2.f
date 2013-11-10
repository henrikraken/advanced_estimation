      program main

*     Test cholesky that uses level 2 BLAS on large problem


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

*     parameter (m_dim=100)

      integer m, n, i, j, nprint_m, nprint_n, nprint, nblock
      integer counter, obs_remaining, obs_this_time
      integer i_obs

*.....Use h_y later for space for covariance.
*     Nope, I don't do that here, because I am just testing 
*     cholesky for now.

*     real h_y( n_dim, n_dim+1 ), xhat(n_dim)
*     real h_y( m_dim, n_dim+1 ), xhat(n_dim)
*     real xhat(n_dim)
      real (kind = double) :: xhat(n_dim)

*     Need a dimension statement for h_y, even though not using
*     it in this version, since it is after the stop statement.
      real (kind = double) :: h_y( 1 , 1 )

      real (kind = double) :: htwh( n_dim, n_dim ), htwy(n_dim)
      real (kind = double) :: diff, max_diff, sum, rms
      real (kind = double) :: start1, end1, time1
      real (kind = double) :: start2, end2, time2
      real (kind = double) :: end3, time3
      real (kind = double) :: Lambda_max, Lambda_min, C

      real (kind = double) :: secondr
      real (kind = double) :: start1b, end1b, time1b, start2b, end2b, time2b
      real (kind = double) :: ops, speed
      real (kind = double) :: y, h_row( n_dim)
      real (kind = double) :: seed, n_real, m_real

      integer count_rate, count_max
      integer count_start_1, count_end_1
      integer count_start_2, count_end_2
      integer count_start_3, count_end_3



      print *, ' '
      print *, 'Cholesky BLAS2'

*.....Set m, n

      print *, 'Enter value for m (no max )'
      read *, m
      print *, 'Enter value for n (max of ', n_dim, ')'
      read *, n
      print *, 'm= ', m, '  n= ', n

      if ( n .gt. n_dim ) then
        print *, 'ERROR:   n= ', n, ' > n_dim= ', n_dim,
     >           ' in test_chol_blas2'
        call abort2
      endif

*.....Set values that tell how many rows and columns at the
*     top of matrices to print.

      nprint = 0
      nprint = 7
      nprint_m = min0( nprint, m )
      nprint_n = min0( nprint, n )

      call system_clock(count_start_1, count_rate, count_max)

*.....Zero out HTWH, HTWY.

      do j = 1, n
        htwy( j ) = 0.0d0
        do i = 1, n
          htwh( i, j ) = 0.0d0
        end do
      end do

*.....Loop to read in obs.

      do i_obs = 1, m

        call set_h_y_i( h_row, y, n, i_obs, seed )

*.......Add this obs into htwh.
*       Only need lower triangular part.
*       htwh = htwh + h_row_trans * h_row

        call dsyr( 'Lower',  n, 1.d0, h_row, 1, htwh, n_dim )

*.......Add this obs into htwy.
*       htwy = htwy + h_row * y

        call daxpy ( n, y, h_row, 1, htwy, 1 )

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

      call chol_blas2( htwh, n_dim, n )

      call system_clock(count_end_2, count_rate, count_max)
      time2 = dble(count_end_2 - count_start_2) / dble(count_rate)

*.....Calculate speed ( Mflops )

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

*.....Check elements of xhat.

      rms = 0.0d0
      do i = 1, n
        diff = xhat( i ) - dble( i )
        rms = rms + diff*diff
      end do
      rms = dsqrt( rms / dble( n ) )

      print 1200, rms
      print 6000, time3, time2, speed


      stop




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
 1010 format ( 1pe12.4 )
 1011 format ( f8.4 )
 1012 format ( 1pe9.1 )
 1100 format ( / 1x, 'Cholesky BLAS2 time =', 1p1e11.3,
     >           1x, 'Total time =',  1p1e11.3 / )
 1101 format ( / 1x, 'Cholesky BLAS2 time =', 1p1e11.3,
     >           1x, 'Total time =',  1p1e11.3 / )
 1200 format (//' Magnitude of error in xhat= ', 1pe10.2 / )
 1400 format (// 1x, 'Accumulation        :',
     >           1x, 'Ops = ',   1pe9.3,
     >           2x, 'Time = ',  1pe9.3, ' sec',
     >           2x, 'Speed = ', 1pe9.3, ' Mf'    )
 1410 format (// 1x, 'Cholesky BLAS2      :',
     >           1x, 'Ops = ',   1pe9.3,
     >           2x, 'Time = ',  1pe9.3, ' sec',
     >           2x, 'Speed = ', 1pe9.3, ' Mf'    )
 2000 format ( / a80 )
 3000 format ( / 1x, 'Cholesky BLAS2: time=', 1p1e10.3,
     .           3x, 'Speed (Mflops)=', 1p1e10.3 //  )
 4000 format ( / 1x, 'Accumulation: time=', 1p1e10.3,
     .           3x, 'Speed (Mflops)=', 1p1e10.3 //  )
 6000 format (// 1x, 'Total time=',     1p1e10.3,
     .           3x, 'Cholesky BLAS2: time=', 1p1e10.3,
     .           2x, 'Speed (Mflops)=',     1p1e10.3 // )


      stop
      end
