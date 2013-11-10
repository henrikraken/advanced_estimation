      subroutine set_to_ident( a, leading_dim, n )

* Purpose: Set matrix to identity

      implicit none

      integer n, leading_dim
      real a(leading_dim, n)

      integer i, j

      do j = 1, n
        do i = 1, n

          a( i, j ) = 0.0

        enddo
      enddo

      do j = 1, n

        a( j, j ) = 1.0

      enddo

      return
      end
