      subroutine set_h_y_i( h_row, y, n, i_obs, seed )

* Purpose: Set the row of h and the value of y
*          for obs number i_obs.
*
* Returns the row of h in the vector h_row.


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
* In this series of commands the parameter double has been assigned the kind value for double-precision,
* since SELECTED_REAL_KIND (12) means select the kind value to give 12 digits of accuracy, which will
* require double precision.  The kind for A and B has been assigned the number for double-precision.  Code written with
* such declarations are portable without corrections. 

      integer, parameter  ::  single = selected_real_kind (6)
      integer, parameter  ::  double = selected_real_kind (12)

*     seed needs to be dimension 2 for F90.
      integer, parameter  ::  seed_dim = 2



      integer n, i_obs
      real (kind = double) :: h_row( n ), y
      real (kind = single) :: random_j
      integer seed(seed_dim)

      integer i, j, k



*     print *, 'single, double, seed_dim print= ', 
*    >          single, double, seed_dim



*.....If this is first obs:
*     Initialize random number generator.  
*     Check dimension of seed.
*     Set seed.

      if ( i_obs .eq. 1 ) then
        call random_seed
        call random_seed(size=k)
        write(*,*) ' Number of integers for random starting value = ', k
        if ( k .ne. seed_dim ) then
            print *, 'ERROR:  Number of integers for random starting',
     >              ' value = ', k, ' not = seed_dim= ', seed_dim,
     >              ' in set_h_y_i'
            call abort2
        endif
        seed(1)=12345
        seed(2)=67890
        call random_seed(put=seed(1:seed_dim))
*       print 1400, i_obs, n, seed
 1400 format ( 'Set seed: i_obs, n, seed= ', 4i10 )
      endif

*.....Set row H with random numbers between -10 and 10.
*     Random number generator returns random numbers between 0 and 1.
*     Then mult. by 20 and subtract 10 to get
*     random numbers between -10 and 10.

      do j = 1, n
        call random_number(random_j)
        h_row( j ) = ( 20.d0*dble(random_j) ) - 10.d0
*       print *, j, random_j, h_row(j)
      enddo

*.....Calculate y using x(j) = j
*     That is, the true value of x(j) is j.

      y = 0.0d0

      do j = 1, n
        y = y + h_row( j ) * dble(j)
      enddo


      return
      end
