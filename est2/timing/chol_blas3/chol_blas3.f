      subroutine chol_blas3( A, lda, n, nblock )

* Purpose: Cholesky decomposition of A into L * (L').
*          L is lower triangular.
*          L' is L transpose.
*          L overwrites lower triangular part of A.
*          Level-3 BLAS version.
*
*  A - input positive definite matrix
*  lda - leading dimension of A from dimension
*                statement in calling program.
*  n - size of A
*  nblock - size of partition of A

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



      integer n, lda, size, nblock
      integer a11_start, a11_size,
     >        a22_start, a22_size,
     >        a21_start_row, a21_start_col, a21_rows, a21_cols,
     >        acur_size

      real (kind = double) :: a(lda, n), alpha

      integer i, j, k

*                       __        __
*                       |          |
*.....Partition A into  | A11  *   |
*                       |          |
*                       | A21  A22 |
*                       |          |
*                       --        --
*
* Factor A into LL', where L is partitioned like A.
*
*  __        __    __        __  __          __
*  |          |    |          |  |            |
*  | A11  *   |    | L11  0   |  | L11'  L21' |
*  |          | =  |          |  |            |
*  | A21  A22 |    | L21  L22 |  |  0    L22' |
*  |          |    |          |  |            |
*  --        --    --        --  --          --
*
* where * means this is a symmetric part that is not needed.
* and where ' means transpose.
*
* From the equation A = LL', you get the equations to solve:
*
*  A11 = L11 * L11'
*  A21 = L21 * L11'
*  A22 = L21 * L21' + L22 * L22'
*
* Factor A11 using level 2 Cholesky to give L11.
* L11 overwrites A11.
*
* Calculate L21 using triangular solve with multiple right
* hand sides (trsm).
* L21 overwrites A21.
*
* Calculate L22 * L22' = A22 - L21 * L21' using symmetric
* rank-k update (syrk).
* L22 * L22' overwrites A22.
*
* Then begin again with A22 as A. 
* Call this ACURRENT, or ACUR.
*
      a11_start = 1
      acur_size = n

      do while ( acur_size > 0 )

*.......Calculate size of a11.
*       Make it nblock x nblock if there is that much of A left.
*       Otherwise, make it as big as the remaining part of A.

        a11_size = min0( nblock, acur_size )

*.......Calculate A21 and A22 locations and sizes.

        a22_start = a11_start + a11_size
        a22_size  = n - a22_start + 1
        a21_rows = a22_size
        a21_cols = a11_size
        a21_start_row = a22_start
        a21_start_col = a11_start

*.......Do a Cholesky using level 2 BLAS on A11 to get L11.
*       This overwrites A11.

        call chol_blas2( a( a11_start, a11_start ), lda, a11_size )

*.......If there is more of A left, update it.

        if ( a22_size .gt. 0) then

*.........Calculate L21 using triangular solve with multiple
*         right hand sides.
*         This overwrites A21.

          call dtrsm( 'Right', 'Lower', 'Transpose', 'Non-unit',
     >                 a21_rows, a21_cols, 1.d0,
     >                 a(a11_start,a11_start), lda,
     >                 a(a21_start_row,a21_start_col), lda )

*.........Calculate L22 * L22' using symmetric rank-k update.
*         This overwrites A22.

          call dsyrk ('Lower', 'No Trans', a22_size, a21_cols,
     >                 -1.d0,a(a21_start_row,a21_start_col), lda,
     >                 1.d0 ,a(a22_start,a22_start), lda )

        end if

*.......Start again, factoring A22.

        a11_start = a22_start
        acur_size = a22_size

      end do


      return
      end
