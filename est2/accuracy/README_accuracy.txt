This file is ~ase381a/est2/accuracy/README_accuracy.txt

* There are two directories in the accuracy directory: chol_blas3 and
house_blas3

* Each of the directories in the accuracy directory has a Makefile, a
run script for running the case in that directory, and the driver
program and all of the necessary subroutines, except for the subroutine
or subroutines that you are to code up.

You need to add your algorithm code into each directory.  Here is what
you need to add:

chol_blas3 - Add your Cholesky BLAS3 subroutine called chol_blas3 in a
file called chol_blas3.f .  In addition, put a copy of chol_blas2.f in
this directory, because your chol_blas3 subroutine will call chol_blas2.
Do this with the following commands:

cd
cd est2
cp timing/chol_blas3/chol_blas3.f accuracy/chol_blas3
cp timing/chol_blas2/chol_blas2.f accuracy/chol_blas3

house_blas3 - You don't need to add anything.

* I create a near-singular matrix by setting the values in col 1 of H as
follows:

col 1 =  (epsilon)*(col 1) + ( 1-epsilon)*(col 2)

For epsilon=1, col 1 does not change.

As epsilon nears zero, col 1 and col 2 become more nearly equal.  If
two cols are identical, then there are not n linearly independent cols
(observations), and the system is not full rank (rank n).  Thus, as
epsilon goes to zero, H becomes more nearly not full rank, and the
system gets closer to being singular.

You are to make H nearly rank deficient in some other way.
Some possibilities:

** Use two other cols besides col 1 and col 2 in the above equation.

** Replace col i with (epsilon)*(col i), where 0 < i < n.
So as epsilon gets smaller, col i becomes more nearly equal zero.

** Replace col i with
(epsilon)*(col i) + (epsilon-1)*( alpha*(col j) + beta*(col k))
So as epsilon gets smaller, col i becomes more nearly equal to a
linear combination of col j and col k.

To do this, you will need to change the code in the routine set_h_y_i.f.
The code to change looks like this:

*.....Make col 1 of H equal to (epsilon)*(col 1) + ( 1-epsilon)*(col 2)

        h_row( 1 ) = epsilon * h_row( 1 ) + ( 1.d0 - epsilon ) * h_row( 2 )

Your changed code needs to be in columns 7-72.

You can continue a code line by putting a ">" in col 6 of lines after
line 1.  For example:

        h_row( 1 ) = epsilon * h_row( 1 ) + ( 1.d0 - epsilon ) *
     >               (  h_row( 2 ) + h_row( 3 ) + h_row( 4 ) +
     >                  h_row( 5 ) + h_row( 6 ) + h_row( 7 )   )


The set_h_y_i.f routine is in both the accuracy/chol_blas3 and
accuracy/house_blas3 dirs.  They are identical.  Change the code in one
and copy that to the other.  For example, change the code in 
accuracy/chol_blas3/set_h_y_i.f and copy that file to
accuracy/house_blas3/set_h_y_i.f.

* I have created scripts and driver routines that use the Block Cholesky
and Block Householder routines to process the observations.  The
program takes epsilon as input.

The program will print out the error in xhat.  The actual quantity it prints is
the rms of ( xhat(i)-x(i) / x(i) ), where x(i) is the true value.

You should see that the errors in xhat increase as epsilon decreases.
You should also see that the errors are larger for the normal equation
method (Cholesky).

Note on your plot which method you use to make H nearly rank
deficient.  We can then see if the results vary depending on the method
for making H nearly rank deficient.

* My job.runall is set up for the small case (n=1800, m=3600).  It will
run this case for epsilon values starting at 1.0 and going down to 1.e-15.

* Use the best values of nblock and mblock you found for Cholesky BLAS3,
and for nblock you found for Householder BLAS3.  Put these into
job.runall before you submit it.

Also submit job.runall using two other values of nblock for Cholesky BLAS3
and for House BLAS3.  This will change the order that the arithmetic is
done, and we can see if the error changes.

So you will have 3 cases for Cholesky BLAS3 and 3 for House BLAS3.  Plot
the error in xhat vs. epsilon for the these 6 cases, all on one plot,

* Something to note about the accuracy assignment for the Cholesky is
that you may get a core dump (a system abort) for small values of
epsilon, or get a result of NaN (not a number) for the error in xhat.
It's not a problem.  Just plot the points for all of the epsilons for
which you get a result.

Here is what happens.  The Cholesky routine takes the square root of the
diagonal element at the top of each loop.  If the matrix is not singular,
that element should never be negative.  However, as epsilon gets very
small, the matrix gets very close to being singular, in which case roundoff
can result in that diagonal element becoming negative.  When the routine
tries to take the sqrt of a negative number, you get an abort.

Since each of you is making the matrix singular by a different method,
some of you may get this situation, and some may not.  If it happens to
you, it is OK.  Just plot the points for the values of epsilon that
succeeded.
