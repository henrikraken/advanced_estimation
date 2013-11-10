      program main

*     Test Householder BLAS3.

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

      parameter (n_dim=4000)
      parameter (m_dim=16000)

      integer m, n, i, j, nprint_m, nprint_n, nprint
      integer mcount, m_batch, ldh, ldp, m_batch_max
      integer this_batch_size, obs_remaining, row, rows, obs, rows_upper
      integer nblock

      real (kind = double) :: h_y( m_dim, n_dim+1 )
      common /h/ h_y
      real (kind = double) :: h_row( n_dim ), y

      real (kind = double) :: xhat(n_dim)
      real (kind = double) :: diff, max_diff, sum, rms
      real (kind = double) :: mag_xhat, mag_deltaxhat, rel_err
      real (kind = double) :: start1, end1, time1, start2, end2, time2
      real (kind = double) :: ops, speed
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
      print *, 'Householder BLAS3'

*.....Set m, n

      print *, 'Enter value for m (max = ', m_dim, ')'
      read *, m
      print *, 'Enter value for n (max = ', n_dim, ')'
      read *, n
      print *, 'Enter value for nblock (max = ', n, ')'
      read *, nblock
      print *, 'Enter epsilon'
      read *, epsilon
      print *, 'm= ', m, '  n= ', n, '  nblock= ', nblock,
     >         '  epsilon= ', epsilon

      if ( m .gt. m_dim ) then
        print *, 'ERROR:  m= ', m, ' > m_dim= ', m_dim,
     >           ' in test_house_blas3'
        call abort2
      endif

      if ( n .gt. n_dim ) then
        print *, 'ERROR:  n= ', n, ' > n_dim= ', n_dim,
     >           ' in test_house_blas3'
        call abort2
      endif

      if ( nblock .gt. n ) then
        print *, 'ERROR:  nblock= ', nblock, ' > n= ', n,
     >           ' in test_house_blas3'
        call abort2
      endif

*.....Set values that tell how many rows and columns at the
*     top of matrices to print.

      nprint = 0
      nprint = 7
      nprint_m = min0( nprint, m )
      nprint_n = min0( nprint, n )

*.....Leading dimension of h

      ldh = m_dim

*.....Fill h_y.

      call system_clock(count_start_1, count_rate, count_max)

      do i = 1, m
        call set_h_y_i( h_row, y, n, i, seed, epsilon )
        do j = 1, n
          h_y( i, j ) = h_row( j ) 
        end do
        h_y( i, n+1 ) = y
      end do
      obs = n
 
*.....Print by rows.
 
      if ( nprint .gt. 0 ) then

        print *, ' '
        print *, 'h'
        do i = 1, nprint_m
          print 1002, i, ( h_y( i, j ), j = 1, nprint_n )
        end do

        print *, ' '
        print *, 'y'
        do i = 1, nprint_m
          print 1012, h_y( i, n+1 )
        end do

      endif

*.....Upper triangularize h.

      call system_clock(count_start_2, count_rate, count_max)

      call house_blas3 (h_y, ldh, m, n, nblock)
  
      call system_clock(count_end_2, count_rate, count_max)
      time2 = dble(count_end_2 - count_start_2) / dble(count_rate)

*.....Calculate speed ( Mflops )

      n_real = dble(n)
      m_real = dble(m)
      ops = ( 2.d0 * m_real * n_real * n_real )
     >    - (2.d0/3.d0) * ( n_real * n_real * n_real )
      speed = ( ops / time2 ) / 1.d+6
      print 1400, ops, time2, speed

      if ( nprint .gt. 0 ) then
        print *, ' '
        print *, 'R'
        do i = 1, nprint_m
          print 1002, i, ( h_y( i, j ), j = 1, nprint_n )
        end do
      endif

*.....Solve for xhat.

      call dcopy( n, h_y( 1, n+1 ), 1, xhat, 1 )
      call dtrsv( 'Upper', 'No trans', 'Non-unit',
     >             n, h_y, ldh, xhat, 1 )
  
      call system_clock(count_end_1, count_rate, count_max)
      time1 = dble(count_end_1 - count_start_1) / dble(count_rate)

      print 6000, time1, time2, speed

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


*.....Check elements of xhat.  Max of abs((error in xhat(i)) / x(i).)

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






 1002 format ( 'Row ', i2, ':', 1p8e9.1, (7x, 1p9e9.1)  )
 1012 format ( 1pe9.1 )
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
 1400 format (// 1x, 'Householder BLAS3   :',
     >           1x, 'Ops = ',   1pe9.3,
     >           2x, 'Time = ',  1pe9.3, ' sec',
     >           2x, 'Speed = ', 1pe9.3, ' Mf'    )
 6000 format (// 1x, 'Total time=',   1p1e10.3,
     .           3x, 'Householder BLAS3: time=', 1p1e10.3,
     .           2x, 'Speed (Mflops)=',   1p1e10.3 // )
  

      stop
      end
