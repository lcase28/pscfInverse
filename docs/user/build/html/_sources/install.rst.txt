
************
Installation
************

Software Requirements
=====================

Python Environment
------------------

This program is written for use on Unix-type operating systems.
Use of PscfInverse requires at least Python 3.5 or later, but Python 3.7 is the minimum
recommended version, as this is the Python version used in development.
The following libraries should be included in most standard Python installations,
but are required for operation of the software

 * abc
 * collections
 * contextlib
 * copy
 * enum
 * functools
 * io
 * itertools
 * multiprocessing
 * numbers
 * os
 * pathlib
 * psutil
 * re
 * signal
 * subprocess
 * time

In addition, the following non-standard libraries are also required:

 * dask
 * numpy
 * networkx

Users can verify that these libraries are available in the Python environment by
opening an interactive python terminal and attempting to import each.

Finally, installation of the field generation tool, pscfFieldGen, is
also required in the python environment. Information about these can be found
on the `PSCF website <https://pscf.cems.umn.edu/>`

Other Dependencies
------------------

The open-source SCFT software, PSCF, must be installed on the system.

PscfInverse must be run in an MPI environment using, for example, OpenMPI.
Such MPI implementations are typically already installed in most
supercomputing systems.

Hardware Requirements
=====================

Execution of PscfInverse is a very computationally intensive process, typically
requiring more computer power than is available on a typical desktop system
in order to run in a reasonable time. Users are advised to run these PSO searches
on high-performance computing clusters, or dedicated high-performance computer systems.
This is because PSO runs typically require the execution of hundreds or thousands
of SCFT calculations and PscfInverse searches will complete most quickly when these
SCFT calculations can be significantly parallelized using the high numbers of processing
cores available on supercomputing systems.


Obtaining the Source Code
=========================

Presently, the source code is available for download from the PSCF
website. After downloading the .zip file, extract the contents to
the directory you wish to use as the root directory for the install.
In descriptions that follow, the filler ``path/to/root`` will be used
to indicate the path to this directory.
After extracting the source, the repository should be contained in 
a new directory called ``path/to/root/pscfInverse``.

Modifying Python Paths
======================

To allow Python to find this library, it must be added to the python search path.

PythonPath
----------

Many python installations make use of the environment variable PYTHONPATH when
searching for libraries. To add PscfInverse to this search path, use the following
command:

::

    PYTHONPATH=$PYTHONPATH:/path/to/root/pscfInverse

Executing this on the command line only modifies the path until the end of
the terminal session. To make the change permanent,
add the above command to the file ~/.bashrc (on linux) or to ~/.profile (on Mac OS).

Anaconda Python
---------------

For Anaconda Python and other conda-managed environments, changes to PYTHONPATH are not
reliably reflected in the python interpreter's path.

To install to a conda-managed environment, first activate the desired environment using
``conda activate`` before proceeding. From here, the easiest way to install PscfInverse
is using the conda-develop command included in the conda-build package.
Install conda-build using

::

    conda install conda-build

After completing that installation, use the command

::

    conda-develop /path/to/root/pscfInverse

After this, PscfInverse should now be available in your python environment.

