      subroutine house_blas3 (h_y, leading_dim, m, n, r)


* Purpose: Zero out elements in h_y below the diagonal in the 
*          first n columns using householder blockqr transformations.
*          h_y contains the m by n matrix H, plus y in the last column,
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



      integer m, n, leading_dim, r
      real (kind = double) :: h_y (leading_dim, n+1), V (m), scr (n)
      real (kind = double) :: W(m, r), Y(m, r)
      real (kind = double) :: temp (r, n+1)

      integer i, j, c, rows, columns
      real (kind = double) :: sigma, beta (r)



*..... Loop through diagonal elements of H, which is in the first
*      n columns of h_y. Zero out elements below diagonal using 
*      blockqr householder transformations.  
  
      i = 1
      do
          
           if (i .gt. n) then
              exit
           endif


*..... The last block may have less than r columns.

           if (n-i+1 .lt. r) then
              r = n - i + 1
           endif


           call zero_out (W, m, r)
           call zero_out (Y, m, r)


*..... Do loop within the block.

           c = 1

           do j = i, i + r - 1
      
             rows = m - j + 1      
             columns = i + r - j


*..... Calculate householder vector Y, sigma, beta for this column.

*            call house (h_y (j:m, j), rows, Y (j:m, c), sigma, 
*    >                   beta (c))

             call house (h_y( j, j ), rows, Y( j, c ), sigma, 
     >                   beta( c ) )

             
*..... Apply householder transformation only to other columns in block 
*      (width r), not the entire remaining h_y matrix.

*            call householder_2 (h_y (j:m, j:(i+r-1)), rows, rows, 
*    >                  columns, Y (j:m, c), sigma, beta (c), scr)

             call householder_2 ( h_y( j, j ), leading_dim, rows, 
     >                  columns, Y( j, c ), sigma, beta( c ), scr)

             c = c + 1
     
           enddo


*..... Create the W and Y transform matrices using householders and
*      beta found above. 


*..... First column of W is -beta(1)*Y(1:m,1). W was initialized to be zero.

*          call daxpy (m, -beta(1), Y (1:m, 1), 1, W (1:m, 1), 1)

           call daxpy (m, -beta(1), Y( 1, 1 ), 1, W( 1, 1 ), 1)


*..... Set first column of temp as zero.

           call zero_out (temp, r, 1)


*..... Calculate W for the next blocks.
 
           do j = 2, r

*..... First, set V = Y(:, j)

*          call dcopy (m, Y(1:m, j), 1, V ,1)

           call dcopy (m, Y( 1, j ), 1, V ,1)


*..... Calcuate z in the W update; W = [W z], where, z = -beta(j)*
*      (V + W*Y'*V). store Y'*V in the first column of temp.

*          call dgemv ('T', m, j-1, 1.d0, Y (1:m, 1:(j-1)), m, V, 
*    >                   1, 0.d0, temp (1:(j-1), 1), 1)

           call dgemv ('T', m, j-1, 1.d0, Y( 1, 1 ), m, V, 
     >                   1, 0.d0, temp( 1, 1 ), 1 )


*..... Now, z = -beta(j)*(V + W*temp) = -beta(j)*V - beta(j)*W*temp.

*          call daxpy (m, -beta(j), V, 1, W(1:m, j), 1)  

           call daxpy (m, -beta(j), V, 1, W( 1, j ), 1)  

*          call dgemv ('N', m, j-1, -beta(j), W(1:m, 1:(j-1)), 
*    >                   m, temp (1:(j-1), 1), 1, 1.d0, W(1:m, j), 1)

           call dgemv ('N', m, j-1, -beta(j), W( 1, 1 ), 
     >                   m, temp( 1, 1 ), 1, 1.d0, W( 1, j ), 1)


           enddo


*..... Set temp as zero.

           call zero_out (temp, r, n+1)         
  

*..... Apply WY transform to the rest of h_y matrix (including last column y); 
*      h_y = (I + W*Y')'*h_y = h_y + Y*(W'*h_y).

           rows = m - i + 1
           columns = n + 1 - (i + r - 1)


*..... Store W'*h_y in the temp.

*          call dgemm ('T', 'N', r, columns, rows, 1.d0, W (i:m, 1:r), 
*    >                 rows, h_y (i:m, (i+r):(n+1)), rows, 0.d0, 
*    >                 temp (1:r, (i+r):(n+1)), r)

           call dgemm ('T', 'N', r, columns, rows, 1.d0, W( i, 1 ), 
     >                 m, h_y( i, (i+r) ), leading_dim, 0.d0, 
     >                 temp( 1, (i+r) ), r)


*..... h_y = h_y + Y*temp

*          call dgemm ('N', 'N', rows, columns, r, 1.d0, Y (i:m, 1:r), 
*    >                 rows, temp (1:r, (i+r):(n+1)), r, 1.d0, 
*    >                 h_y (i:m, (i+r):(n+1)), rows)

           call dgemm ('N', 'N', rows, columns, r, 1.d0, Y( i, 1 ), 
     >                 m, temp( 1, (i+r) ), r, 1.d0, 
     >                 h_y ( i, (i+r) ), leading_dim )


      i = i + r


      enddo                      


       
      return
      end



****************************************************************************

      subroutine zero_out (x, m, n)

*..... Zero out a (m by n) matrix.

      implicit none
      integer, parameter  ::  double = selected_real_kind (12)

      integer i, j, m, n
      real (kind = double) :: x(m,n)

      do i = 1, m
           do j = 1, n
               x(i,j) = 0.d0
           enddo
      enddo

      return
      end      



*****************************************************************************

      subroutine house (x, m, v, sigma, beta)

*..... Calculate householder vector v corresponding to the vector x. 

      implicit none
      integer, parameter  ::  double = selected_real_kind (12)

      integer m
      real (kind = double) :: x(m), v(m), sigma, ddot, beta, mu

      sigma = ddot (m, x, 1, x, 1) - x(1)*x(1)
      v(1) = 1.d0
      v(2:m) = x(2:m)
      
      if (sigma .eq. 0.d0) then
         beta = 0.d0
      else
         mu = dsqrt (x(1)*x(1) + sigma)
         if (x(1) .le. 0.d0) then
            v(1) = x(1) - mu
         else 
            v(1) = -sigma/(x(1) + mu)
         endif

         beta = 2.d0*v(1)**2/(sigma + v(1)**2)
         v(1:m) = v(1:m)/v(1)
      endif

      return
      end  



*******************************************************************************

      subroutine householder_2 ( h, leading_dim, m, n, u, sigma, beta,
     >                          gamma )

* Purpose: Multiply a matrix by a householder.
*          It is assumed that the householder has been
*          calculated to zero out the elements below the diagonal
*          element in column 1.

      implicit none
      integer, parameter  ::  double = selected_real_kind (12)

      integer leading_dim, m, n
      real (kind = double) :: h( leading_dim, n)
      real (kind = double) :: u( m ), sigma, beta, gamma(n)

      integer i__, j__
      real (kind = double) :: dot
      real (kind = double) :: sdot

*.....Multiply remaining columns by householder.
*     First, calculate beta * (u dot columnj)
*     for each of the columns 1 through n.
*     This can be done at one time using gemv.
*     Store results in gamma vector in elements 2 through n.

*     Calculate T * h,
*     where T = householder = I - beta * u * u'
*     where u' = u transpose
*
*     T * h = ( I - beta * u * u' ) * h 
*           = h - beta * u * u' * h 

*     First, calculate gamma = beta * u' * h.
*     This is a matrix vector multiply.

      call dgemv( 'Trans', m, n, beta, h, leading_dim, u, 1,
     >            0.d0, gamma, 1 )

*     do j__ = 1, n
*       gamma(j__) = 0.d0
*       do i__ = 1, m
*         gamma(j__) = gamma(j__) + u(i__) * h(i__,j__)
*       end do
*       gamma(j__) = beta * gamma(j__)
*     end do

*.....Now calculate
*
*     h = h - beta * u * u' * h = h - * u * gamma
*     
*     This is a rank-1 update to h.

      call dger( m, n, -1.d0, u, 1, gamma, 1, h, leading_dim )

      return
      end
