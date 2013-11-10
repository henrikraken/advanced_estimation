This file is ~ase381a/est2/README_assignment.txt


Copy the est2 directory, along with its subdirectories and files, to
your ase account by giving the following commands:

cd
cp -r ~ase381a/est2 .



(1) TIMING

This part of the assignment will compare the speed of various normal
equation methods and orthogonal decomposition methods.

See the directory ~ase381a/est2/timing and the
README file there for more detailed information.

You will use the following routines:

* Givens 
* Givens operating on H transpose (to take advantage of computer cache)
* Cholesky (BLAS2)
* Householder (BLAS2)
* Block Cholesky (BLAS3)
* Block Householder (BLAS3)

Run each routine on the problem I set up in
~ase381a/est2/timing.  We'll work with the following
problem sizes:

 n=1800 and m=3600
 n=4000 and m=16000

For each problem size, you will run each routine 3 times.  I have a
script and a driver routine for each of the routines above.  My driver
routine will print out the speed of each routine and the total time to
calculate the estimate.

After you have made 3 runs per routine, plot routine speed vs. the total time.
So total time on the x-axis and routine speed on the y-axis.

The block routines, which use BLAS3 (matrix-matrix operations) should
have the fastest speed and the lowest total time.

Make a plot for n=1800 and m=3600.  Make another plot for n=4000 and m=16000
(Givens is so slow that my job deck does not run Givens on this large case).
See if results change with problem size.



(2) ACCURACY

This part of the assignment will compare the accuracy of the normal
equation method (using Block Cholesky) to the orthogonal
transformation method (using Block Householder).

See the directory ~ase381a/est2/accuracy and the
README file there for more detailed information.

I create a near-singular matrix by setting the values in col 1 of H as
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


* I have created scripts and driver routines that use the Block Cholesky
and Block Householder routines to process the observations.  The
program takes epsilon as input.  I have things set up to run for
epsilon values starting at 1.0 and going down to 1.e-15.



The program will print out the error in xhat (the rms of the difference
between the elements in xhat and the correct xhat).

You should see that the errors in xhat increase as epsilon decreases.
You should also see that the errors are larger for the normal equation
method (Cholesky).

Do this using your optimum values of nblock for Cholesky BLAS3
and for House BLAS3.

But also do this for two other values of nblock for Cholesky BLAS3 and
for House BLAS3.  This will change the order that the arithmetic is
done, and we can see if the error changes.

So you will have 3 cases for Cholesky BLAS3 and 3 for House BLAS3.  On
one plot, plot the error in xhat vs. epsilon for the these 6 cases.

Note on your plot which method you use to make H nearly rank
deficient.  We can then see if the results vary depending on the method
for making H nearly rank deficient.




(3) Notes on getting the programs to compile.

You have to load the mkl module in order for the programs to compile.
To have this done automatically, give the following commands:

cd
echo "module load mkl" >> .profile_user
chmod u+x .profile_user
echo "module load mkl" >> .login_user
chmod u+x .login_user

This will create a .profile_user file and a .login_user file in your home dir
with the command to load the mkl module.  This will execute when you login.
Creating both of these files will make things work, regardless of which
interactive shell (bash or csh or tcsh) that you are using.

Another thing to note is that I had to take the make command out of the run
script.  So before you give a run command, you have to compile by giving the
make command, like this:

cd est2/timing/chol_blas2
make
run 1000 2000



(4) Turning in your assignments

To get the data for plotting, use grep.

For example, for the timing project, for the n=1800 m=3600 case,
you can grep for Total like this::

login1$ pwd
/home1/00567/ase381b/est2/timing/chol_blas2
login1$ grep -h Total out.chol_blas2.n_1800.m_3600.*
 Total time= 1.073E+01   Cholesky BLAS2: time= 1.052E+00  Speed (Mflops)= 1.848E+03
 Total time= 1.220E+01   Cholesky BLAS2: time= 8.267E-01  Speed (Mflops)= 2.352E+03
 Total time= 8.732E+00   Cholesky BLAS2: time= 9.370E-01  Speed (Mflops)= 2.075E+03

You can send that into an output file like this:

login1$ grep -h Total out.chol_blas2.n_1800.m_3600.* > data.chol_blas2.1800.3600
login1$ cat data.chol_blas2.1800.3600 
 Total time= 1.073E+01   Cholesky BLAS2: time= 1.052E+00  Speed (Mflops)= 1.848E+03
 Total time= 1.220E+01   Cholesky BLAS2: time= 8.267E-01  Speed (Mflops)= 2.352E+03
 Total time= 8.732E+00   Cholesky BLAS2: time= 9.370E-01  Speed (Mflops)= 2.075E+03


You can do the same grep for the output files for the larger size problem for
chol_blas2, and for the other cases like chol_blas3, etc.  That will put all the
values into files, then you can plot them.

For the accuracy project, you want to grep for Eps.  For example:

login2$ pwd
/home1/00567/ase381b/est2/accuracy/chol_blas3
login2$ grep -h Eps out.chol_blas3.n_1800.m_3600.nblock_500.mblock_500.eps_1.e* > data.chol_blas3.500.500
login2$ cat out.chol_blas3.500.500 
 Epsilon=   1.00E+00    Max of abs(( xhat(i)-x(i) ) / x(i)) =   4.68E-12 at i =     1
 Epsilon=   1.00E-01    Max of abs(( xhat(i)-x(i) ) / x(i)) =   2.66E-10 at i =     1
 Epsilon=   1.00E-02    Max of abs(( xhat(i)-x(i) ) / x(i)) =   1.52E-07 at i =     1
 Epsilon=   1.00E-03    Max of abs(( xhat(i)-x(i) ) / x(i)) =   2.61E-06 at i =     1
 Epsilon=   1.00E-04    Max of abs(( xhat(i)-x(i) ) / x(i)) =   4.29E-04 at i =     1
 Epsilon=   1.00E-05    Max of abs(( xhat(i)-x(i) ) / x(i)) =   1.22E-01 at i =     1
 Epsilon=   1.00E-06    Max of abs(( xhat(i)-x(i) ) / x(i)) =   4.21E+00 at i =     1
 Epsilon=   1.00E-07    Max of abs(( xhat(i)-x(i) ) / x(i)) =   5.84E+03 at i =     1
 Epsilon=   1.00E-08    Max of abs(( xhat(i)-x(i) ) / x(i)) =   7.85E+01 at i =     1
 Epsilon=   1.00E-09    Max of abs(( xhat(i)-x(i) ) / x(i)) =   0.00E+00 at i = 12986
 Epsilon=   1.00E-10    Max of abs(( xhat(i)-x(i) ) / x(i)) =   4.09E+02 at i =     1
 Epsilon=   1.00E-11    Max of abs(( xhat(i)-x(i) ) / x(i)) =   0.00E+00 at i = 16185
 Epsilon=   1.00E-12    Max of abs(( xhat(i)-x(i) ) / x(i)) =   4.75E+02 at i =     1
 Epsilon=   1.00E-13    Max of abs(( xhat(i)-x(i) ) / x(i)) =   2.24E+02 at i =     1
 Epsilon=   1.00E-14    Max of abs(( xhat(i)-x(i) ) / x(i)) =   0.00E+00 at i = 12928
 Epsilon=   1.00E-15    Max of abs(( xhat(i)-x(i) ) / x(i)) =   9.15E+02 at i =     2

Do the same for your other two block-size cases for chol_blas3, and for your 3 block-size cases
for house_blas3.  Then you will have your data in 6 output files, and you can plot it.

When you have made your plots, email them to me.
