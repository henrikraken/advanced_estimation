      subroutine givens_h_trans( h_y_T, leading_dim, m, n )

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



      integer m, n, leading_dim
      real (kind = double) :: h_y_T(leading_dim, m)

      integer i, j, k, p
      real (kind = double) :: S, C, temp, tempinv

*.....Loop through rows.

      do 200 k = 2, m

*.......Zero out elements in this row.

*       Zero out elements up to element before the diag,
*       if this row is in top nxn part of matrix.
*       Otherwise, zero out all elements in this row.

        if ( k .le. n) then
          p = k - 1
        else
          p = n
        endif

        do 100 i = 1, p

*.........If this element is not already zero, do Givens rotation.

          if ( h_y_T( i, k ) .ne. 0.d0 ) then

            temp        = dsqrt ( h_y_T( i, i )**2 + h_y_T( i, k )**2 )
            tempinv     = 1.d0 / temp
            S           = h_y_T( i, k ) * tempinv
            C           = h_y_T( i, i ) * tempinv
            h_y_T( i, i ) = temp
            h_y_T( i, k ) = 0.d0

*...........Apply this Givens rotation to remaining elements in
*           rows i and k

            do j = i+1, n+1

              temp        =   C * h_y_T( j, i )  +  S * h_y_T( j, k )
              h_y_T( j, k ) = - S * h_y_T( j, i )  +  C * h_y_T( j, k )
              h_y_T( j, i ) = temp

            enddo 

          endif
 100    continue
 200  continue

      return
      end
