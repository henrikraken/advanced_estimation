      subroutine chol_blas2( a, leading_dim, n )

* Purpose: Cholesky decomposition of A into L * (L transpose).
*          L is lower triangular.
*          L overwrites lower triangular part of A.
*          Level-2 BLAS version.


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



      integer n, leading_dim, size
      real (kind = double) :: a(leading_dim, n), alpha

      integer i, j, k

      do j = 1, n-1

        a( j, j ) = dsqrt( a( j, j ) )

        size = n - j

        alpha = 1.d0/a( j, j )

        call dscal( size, alpha, a( j+1, j ), 1 )

        call dsyr( 'Lower', size, -1.d0, a( j+1, j ), 1, a( j+1, j+1 ),
     >              leading_dim )


      end do

      a( n, n ) = dsqrt( a( n, n ) )

      return
      end
