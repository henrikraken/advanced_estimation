This file is ~ase381a/est2/timing/README_timing.txt

-----------
Timing runs
-----------

* The contents of the timing directory are

Directories to run timing tests on the different algorithms:

chol_blas2
chol_blas3
givens
givens_h_trans
house_blas2
house_blas3

The job deck to run all timing runs is in the timing directory and is:

  job.runall


* Each of the directories in the timing directory has a Makefile, a run
script for running the case in that directory, and the driver program
and all of the necessary subroutines, except for the following things
that you need to do, listed for each case.

In some cases, you will need to fill in missing parameters in some of
the subroutines, as we discussed in class.  Missing parameters in the
code are indicated by "---".


chol_blas2 - You don't need to do anything here.

chol_blas3 - Fill in the missing parameters in chol_blas3.f.

givens - You don't need to do anything here.

givens_h_trans - Add your Givens H_transpose subroutine called
givens_h_trans in a file called givens_h_trans.f

house_blas2 - Fill in the missing parameters in house_blas2.f

house_blas3 - You don't need to do anything here.


* Each directory has a run file for making the program in that
directory and running it.  Each run file is set up to take command line
arguments for n and m, plus nblock and mblock parameters where
needed.  So the command line arguments for each run file are:


chol_blas2:	run n m
chol_blas3:	run n m nblock mblock
givens:		run n m
givens_h_trans:	run n m
house_blas2:	run n m
house_blas3:	run n m nblock


For the project, we will work with a small case and a large case.
For the small case: n=1800, m=3600.
For the large case: n=4000, m=16000.
You must find good values for nblock and mblock.

NOTE, don't run givens for the large case.  It takes too long.  We will
run givens_h_trans on this large case.



* To debug, you will want to use smaller values for m and n.  For
example, for the small case, my cholesky BLAS2 takes about 8 seconds.
But if you use n=500, m=2000, that case takes less than 1 second.  So you
could debug with this very small case and get quick turnaround.
To debug, you would issue one of the above run commands in whichever
directory you are working in.  So if you are debugging your chol_blas3,
in that directory you would give the interactive command

run 500 2000 10 10

for example.  The values of 10 and 10 for nblock and mblock are not the
optimum values, but they are good enough for debugging a small case.

TACC does not like users to run interactive programs on lonestar.  However, it
is OK for you to run interactively on the grace node on lonestar.  So when you
are debugging and running a "run" script interactively, please do so on the
grace node.  That is, login to grace.lonestar.tacc.utexas.edu instead of
lonestar.tacc.utexas.edu.

* The run files for chol_blas3 and house_blas3 have an additional
parameter called nblock.  This is how many columns at a time the
algorithm operates on.  This affects the algorithm speed, which is
given in each output file on lines like the following:

 Cholesky BLAS3      : Ops = 2.133E+10  Time = 1.733E+00 sec  Speed = 1.231E+04 Mf

or

 Householder BLAS3   : Ops = 4.693E+11  Time = 4.424E+01 sec  Speed = 1.061E+04 Mf

Before running the timing case, you are to find the value for this
parameter which gives the fastest algorithm speed.  Cholesky and
Householder may have different optimal nblock values.  I have set
nblock to 10 in job.runall, which is not the optimal value.  To find
the optimal value, run these two cases with various values for nblock,
and see which value of nblock gives the highest algorithm speed.
You can run the small case interactively to do this.


* The run file for chol_blas3 has an additional parameter called
mblock.  This is how many obs at a time to accumulate into HTWH and
HTWY.  This affects the accumulation speed, which is given in each
output file on a line like the following:

 Accumulation        : Ops = 2.560E+11  Time = 2.148E+01 sec  Speed = 1.192E+04 Mf

Before running the timing cases, you are to find the value for this
parameter which gives the fastest accumulation speed.  I have set
mblock to 10 in job.runall, which is not the optimal value.
The values of mblock and nblock are independent.  That is, mblock affects
only the accumulation speed, and nblock affects only the Cholesky speed.
You can run the small case interactively to do this.



* To aid in debugging your routines, I print out the relative error in xhat.

relative error = RMS of ( xhat(i) - x(i) )
The relative error should be much less than 1.0.  Here is one value I got:

 Magnitude of error in xhat=   5.30E-12

* After you have all of your algorithms running and have determined some good
values for nblock and mblock for chol_blas3, and for nblock for house_blas3,
then set those values in job.runall.  Then submit job.runall by giving the
following command in the timing directory:

qsub job.runall

This command will run the job in the batch queue.

The job will go into each directory in the timing directory and run the
small case and the large case in that directory.  It will do this for
all directories in the timing directory.  It will repeat this process
3 times.


* You can monitor your job with these commands:

showq -u           Show your jobs.

If you need to kill your batch job, get the job ID using one of the
above commands.  For example, I have a job running under my byab305
account.  Here is what I get from the above commands:

lslogin1132> showq -u

SUMMARY OF JOBS FOR USER: <byab305>

ACTIVE JOBS--------------------
JOBID     JOBNAME    USERNAME      STATE   PROC  REMAINING       STARTTIME

901917           cp   byab305    Running      1    0:00:53  Thu Oct  7 14:53:31

You have  1 Active jobs utilizing:    1 of 5464 Compute     Processors ( 0.02%)

Total Jobs: 1     Active Jobs: 1     Idle Jobs: 0     Blocked Jobs: 0   


The job ID for my job is 901917.
Then you can kill the job using the qdel command.  For example, to
kill the above job, I would give the command

  qdel 901917


* Algorithm speed and total time are printed in the output files
on lines like the following:

 Total time= 2.352E+01   Cholesky BLAS3: time= 1.733E+00  Speed (Mflops)= 1.231E+04

Plot algorithm speed vs. total time for the small case on one
plot.  Do the same on another plot for the large case.  So total
time on the x-axis and algorithm speed on the y-axis.

Make a plot for n=1800 m=3600.  Make another plot for n=4000 m=16000.
See if results change with problem size.
