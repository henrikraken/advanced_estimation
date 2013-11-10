      program main

*     Test givens.

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

      real (kind = double) :: h_y_T( n_dim+1, m_dim )
      common /h/ h_y_T
      real (kind = double) :: h_row( n_dim ), y

      real (kind = double) :: xhat(n_dim)
      real (kind = double) :: diff, max_diff, sum, rms
      real (kind = double) :: start1, end1, time1, start2, end2, time2
      real (kind = double) :: ops, speed
      real (kind = double) :: seed, n_real, m_real

      integer count_rate, count_max
      integer count_start_1, count_end_1
      integer count_start_2, count_end_2
      integer count_start_3, count_end_3


      print *, ' '
      print *, 'Givens H_trans'

*.....Set m, n

      print *, 'Enter value for m (max = ', m_dim, ')'
      read *, m
      print *, 'Enter value for n (max = ', n_dim, ')'
      read *, n

      if ( m .gt. m_dim ) then
        print *, 'ERROR:  m= ', m, ' > m_dim= ', m_dim,
     >           ' in test_givens_h_trans'
        call abort2
      endif

      if ( n .gt. n_dim ) then
        print *, 'ERROR:  n= ', n, ' > n_dim= ', n_dim,
     >           ' in test_givens_h_trans'
        call abort2
      endif

*.....Set values that tell how many rows and columns at the
*     top of matrices to print.

      nprint = 0
      nprint = 7
      nprint_m = min0( nprint, m )
      nprint_n = min0( nprint, n )

*.....Leading dimension of h_y_T

      ldh = n_dim+1

*.....Fill h_y_T.

      call system_clock(count_start_1, count_rate, count_max)

      do i = 1, m
        call set_h_y_i( h_row, y, n, i, seed )
        do j = 1, n
          h_y_T( j, i ) = h_row( j ) 
        end do
        h_y_T( n+1, i ) = y
      end do
      obs = n
 
*.....Print by rows.
 
      if ( nprint .gt. 0 ) then

        print *, ' '
        print *, 'h'
        do i = 1, nprint_m
          print 1002, i, ( h_y_T( j, i ), j = 1, nprint_n )
        end do

        print *, ' '
        print *, 'y'
        do i = 1, nprint_m
          print 1012, h_y_T( n+1, i )
        end do

      endif

*.....Upper triangularize h.

      call system_clock(count_start_2, count_rate, count_max)

      call givens_h_trans( h_y_T, ldh, m, n )
  
      call system_clock(count_end_2, count_rate, count_max)
      time2 = dble(count_end_2 - count_start_2) / dble(count_rate)

*.....Calculate speed ( Mflops )

      n_real = dble(n)
      m_real = dble(m)
      ops = ( ( 3.d0 * m_real ) - n_real ) * n_real * n_real
      speed = ( ops / time2 ) / 1.d+6
      print 1400, ops, time2, speed

      if ( nprint .gt. 0 ) then
        print *, ' '
        print *, 'R'
        do i = 1, nprint_m
          print 1002, i, ( h_y_T( j, i ), j = 1, nprint_n )
        end do
      endif

*.....Solve for xhat.

      call dcopy( n, h_y_T( n+1, 1 ), ldh, xhat, 1 )
      call dtrsv( 'Lower', 'Trans', 'Non-unit',
     >             n, h_y_T, ldh, xhat, 1 )
  
      call system_clock(count_end_1, count_rate, count_max)
      time1 = dble(count_end_1 - count_start_1) / dble(count_rate)

*.....Check elements of xhat.

      rms = 0.d0
      do i = 1, n
        diff = xhat( i ) - dble( i )
        rms = rms + diff*diff
      end do
      rms = dsqrt( rms / dble( n ) )

      print 1200, rms
      print 6000, time1, time2, speed


 1002 format ( 'Row ', i2, ':', 1p8e9.1, (7x, 1p9e9.1)  )
 1012 format ( 1pe9.1 )
 1200 format (//' Magnitude of error in xhat= ', 1pe10.2 / )
 1400 format (// 1x, 'Givens H_transpose  :',
     >           1x, 'Ops = ',   1pe9.3,
     >           2x, 'Time = ',  1pe9.3, ' sec',
     >           2x, 'Speed = ', 1pe9.3, ' Mf'    )
 6000 format (// 1x, 'Total time=',   1p1e10.3,
     .           3x, 'Givens H_trans: time=', 1p1e10.3,
     .           2x, 'Speed (Mflops)=',   1p1e10.3 // )
  

      stop
      end
