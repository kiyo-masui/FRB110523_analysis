README for FRB MCMC code.  
Author:  Jonathan Sievers
25 October 2015

This code estimates best-fit parameters and error distributions for
FRB110523 using Markov Chain Monte Carlo techniques.  The code is
written in C and GNU Octave, the open-source MATLAB approximation.
While most of the code will translate directly, if you wish to use
MATLAB instead of Octave, you'll have to translate the interpreter/C
interface routines yourself.

Preqrequisites:  Octave, C/C++ compilers (preferable OpenMP-aware),
and libcerf, a complex error function library.  Libcerf can be
obtained from http://apps.jcns.fz-juelich.de/src/libcerf/.  You will
also need to convert the raw data into a filtered/formatted .npy files
using analyse_burst.py included in the Masui directory in this
repository.  

Compilation:  First, the C code burst_fill.c must be compiled into a
shared library.  My compile line is:
gcc -I/home/sievers/local/include -O3 -fopenmp -fPIC --shared --std=c99  burst_fill.c -o libburst_fill.so
where the include path can be updated (or omitted) as-needed.  Note that cerf.h
must be either in a system location or in the -I search path.

Next the octave wrapper routines must be compiled.  My compile line
is:
mkoctfile burst_routines.cpp -I. -L. -L/home/sievers/local/lib -lburst_fill -lcerf -lgomp -v
where again the path to libcerf.so must be included in the library
search path defined by -L.  

Running: After generating the needed .npy files, the shared C library
and the octave interface wrapper burst_routines.oct, make sure that
the locations of libburst_fill.so and libcerf.so are in your
LD_LIBRARY_PATH.  The driver routines are named fit_burst_*.m.  When
run, they will call the other octave/matlab routines as-needed.  The
code will by default try to write into a ./chains directory, so either
create that or set the output file root in the driver routine, defined
by myopts.outroot, to an appropriate location.  While you are free to
come up with your own sampling steps, I include good starting
positions and parameter covariances in the "initial_conditions"
directory, which you should feel free to change as desired.
Similarly, the chain lengths are currently hardwired by myopts.nstep,
which can be changed as-needed.  The likelihood is threaded; using 5
Ivy Bridge cores at 2.8 GHz typical likelihood evaluation times are
0.02-0.07 seconds (depending on what model is being fit).  Using the
given starting positions/covariances, you should expect to have
reasonably converged chains in just a few minutes.

A note on MPI:  I run multiple chains to check convergence.  If you
desire, you can do the same, you'll need to make sure the c++ compiler
is set to be mpic++, and you'll also ned to wrap the MPI function
MPI_Init and MPI_Comm_rank.  The different chains don't communicate,
and the MPI calls are only used to set the processes to have
non-conflicting output names.  If you are happy running single chains,
create a dummy file mpi_comm_rank.m that returns 0.  If you wish to
run multiple chains without using MPI, a cheezy hack would be to have
mpi_comm_rank.m return a random integer from a range large enough that
collisions between chain names would be unlikely.  
