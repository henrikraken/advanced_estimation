      subroutine house_blas2( h_y, ldh, m, n, u, scr)

* Purpose: Zero out elements in h_y below the diagonal in the
*          first n columns using Householder transformations.
*          h_y contains the m by n H matrix, plus y in the last column,
*          so the dimensions of h_y are m by (n+1).

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



      integer m, n, ldh
      real (kind = double) :: h_y(ldh, n+1), u(m), scr(n+1)

      integer i, rows, columns
      real (kind = double) :: beta


*.....Loop through diagonal elements of H, which is in the first
*     n columns of h_y.
*     Zero out elements below diagonal using Householder
*     transformations.

      do i = 1, n

*.......rows = rows in lower right portion of h_y for this step.
*.......columns = columns in lower right portion of h_y for this step.

        rows = ldh - i +1
        columns = n +1 - i +1


*.......Send this matrix into house_0.

        call house_0( h_y( i, i), ldh, rows, columns, u, scr )

      end do


 1000 format ( 'Row ', i2, ':', 1p10e12.4, / (6x, 1p10e12.4)  )
 1001 format ( 'Row ', i2, ':', 10f8.4, / (6x,10e8.4)  )
 1010 format ( 1pe12.4 )


      return
      end


      subroutine house_0( A, ldA, rows_A, cols_A, u, gamma )

* Purpose: Calculate Householder to zero out the first column of A below
*          the 1,1 element.
*          Then apply this Householder to A.

      implicit none

      integer, parameter  ::  double = selected_real_kind (12)

      integer ldA, rows_A, cols_A
      real (kind = double) :: A(ldA, cols_A), u(rows_A), gamma(cols_A)

      real (kind = double) :: beta

*.....Calculate u and beta for Householder for first column.
*     Apply Householder to first column

      call house_1( A( 1, 1 ), rows_A, u, beta )

*.....Apply Householder to rest of A.

      call house_2( A( 1, 2 ), ldA, rows_A, cols_A - 1, u, beta, gamma )


      return
      end


      subroutine house_1( v, rows_v, u, beta )

* Purpose: Calculate u and beta for Householder transformation
*          for vector v.
*          Also apply Householder to v.

      implicit none

      integer, parameter  ::  double = selected_real_kind (12)

      integer rows_v
      real (kind = double) :: v( rows_v ), u( rows_v ), beta

      real (kind = double) :: mag_v, mag_u_sq, dot, sigma
      real (kind = double) :: dnrm2
      integer i


*.....Magnitude of v

      mag_v = dnrm2 (rows_v, v, 1 )

*.....sigma

      sigma = sign( mag_v, v( 1 ) )

*.....u

      call dcopy( rows_v, v(1),1,u(1),1)
      u( 1 ) = u( 1 ) + sigma

*.....beta

      beta = 1.d0 / ( sigma * u( 1 ) )

*.....Resulting v vector.

      v( 1 ) = -sigma
      do i = 2, rows_v
        v( i ) = 0.d0
      end do


      return
      end


      subroutine house_2( B, ldB, rows_B, cols_B, u, beta, gamma )

* Purpose: Multiply a matrix by a Householder.

      implicit none

      integer, parameter  ::  double = selected_real_kind (12)

      integer ldB, rows_B, cols_B
      real (kind = double) :: B( ldB, cols_B), u( rows_B )
      real (kind = double) :: beta, gamma(cols_B)

      real (kind = double) :: dot
      real (kind = double) :: sdot

*     real (kind = double) :: dot_product
*     external dot_product


*     Calculate T * B,
*     where T = Householder = I - beta * u * u'
*     where u' = u transpose
*
*     T * B = ( I - beta * u * u' ) * B 
*           = B - beta * u * u' * B 
*           = B - u * gamma

*     First, calculate gamma = beta * u' * B.
*     This is a matrix vector multiply.
*
*     NOTE: gemv has no capability to calculate beta * u' * B.
*     Instead, it can calcluate beta * B' * u.  Both beta * u' * B and
*     beta * B' * u are vectors.  Technically, beta * u' * B is a row
*     vector, while beta * B' * u is a column vector.  But in Fortran,
*     they are both just vectors.  You yourself have to keep track of
*     whether a Fortran vector represents a row vector or a column
*     vector.  So calculating beta * B' * u
*     and putting it into the gamma vector would give the same
*     Fortran gamma vector as calculating beta * B' * u and putting it
*     into the gamma vector.  So we use gemv to calculate beta * B' * u
*     and put it into the gamma vector, and we remember that represents
*     a row vector.  More on this in the dger step below.

      call dgemv('Transpose', rows_B, cols_B, beta, B, ldB, u, 1,
     >            0, gamma,1 )

*.....Now calculate
*
*     T * B = B - beta * u * u' * B = B - u * gamma
*
*     and overwrite B with this result.
*     This is a rank-1 update to B.
*
*     NOTE:  See note above at the gemv step.
*     gamma represents a row vector.  But Fortran has no way to specify
*     that.  So it is just a vector.  We want to calculate
*     B - u * gamma, where u is a column vector and gamma is a row
*     vector.  But there is no way to tell BLAS that gamma is a row
*     vector.  It thinks it is just a vector.  So in Fortran, we have
*     to treat gamma like it is a column vector, and ask dger to
*     calcluate B - u * gamma'.

      call dger( rows_B, cols_B, -1.d0, u, 1, gamma, 1, B, ldB )

 1001 format ( 'Row ', i2, ':', 10f8.4, / (6x,10e8.4)  )


      return
      end
