.. _usage-main:

*****
Usage
*****

Inputs
======

Running a PSO inverse design calculation requires several files:

* A Parameter File
* For each phase:
    
    * A model file
    * A template PSCF Parameter File

Each of these is described on the Parameter File page.

Command Line Usage
==================

PscfInverse must be run in an MPI environment to allow parallel execution
of SCFT calculations. This means that an MPI implementation must be
available on the system used to run the calculation. On many
supercomputing systems, the MPI implementation can be made available
with (assuming use of the OpenMPI) ``module load ompi``.
After the MPI implementation is loaded, the program will be executed
with the following command:

::

    mpiexec -np [MPI_NODES] python -m pscfinverse -f param

This command can be considered in three parts. The first three
entries (``mpiexec -np [MPI_NODES]``) tell the computer to
execute the program in an MPI environment (``mpiexec``),
and that the MPI environment should use a specific number (*[MPI_NODES*)
of cores/nodes/workers for the execution. 
(On typical high-performance computing systems, mpiexec will automatically
determine the specific hardware available to the calculation and will
spread the MPI "nodes" accross this hardware.)

The next three terms (``python -m pscfinverse``) tell the system to start a python
interpreter and execute the ``__main__`` method of the ``pscfinverse`` python
package, which is, of course, the main method of PscfInverse.

The final two terms (``-f param``) are arguments passed to the main method of
PscfInverse. The flag ``-f`` tells PscfInverse to expect the name of the 
parameter file for this run, and ``param`` is the name of the parameter file.

Output from the program can be redirected directly to a file with the following
command:

::

    mpiexec -np [MPI_NODES] python -m pscfinverse -f param > logfile

which is identical to the first command except that ``> logfile`` tells the
operating system to redirect the output to a file named *logfile*.

