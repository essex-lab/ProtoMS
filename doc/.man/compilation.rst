*****************************
Compilation and Installation
*****************************

The ProtoMS package supplies the following files and directories;

* **CMakeLists.txt** This file configures ProtoMS prior to compiling.

* **data** This directory contains a number of useful files, e.g. pre-equilibrated boxes and some template files

* **doc** This directory contains documentation

* **README** File that contains brief installation instructions for ProtoMS, and any last minute addendums or errata that arrived too late to make it into the manual!

* **parameter** This directory contains all of the standard parameter files that describe the standard forcefields implemented in ProtoMS.

* **protoms.py** A tool to setup common ProtoMS calculations

* **src** This directory contains all of the source code for the main program

* **tools** This directory contains numerous useful scripts to setup and analyse ProtoMS simualtions.

* **tutorial** This directory contains a number of examples that demonstrate applications of ProtoMS.


.. _fortran77:

====================
Programming Language
====================

ProtoMS is written in slightly extended Fortran 77. The extensions used are

* The maximum line length is up to 132 characters, rather than 72.

* Variable, subroutine and function names are greater than 6 characters.

* ``do/enddo`` loops are used rather than ``do/continue``.

* Fortran ``include`` is used to include the contents of other files.

* The ``flush``, ``getarg`` and ``getenv`` non-standard intrinsic functions are used.

* ProtoMS performs string manipulation using the ``len`` function. In addition, the string manipulation assumes the same string handling behaviour as the GNU Fortran compiler (g77), so there is the possibility of strange formatting bugs when using different compilers.

* The ``Date`` and ``Time`` Fortran 90 intrinsic subroutine is used to get the current time. This is used to provide a default seed to the random number generator. This can be removed by commenting out the relavant lines in ``getoptions.F``, though you will need to provide a random number seed manually.

The tools are written in python (https://www.python.org/) and should be compatible with the standard implementation version 2.7. It is known that is does not work with 2.6 and it has not been tested with version 3.0 or more. 

=================
Building ProtoMS
=================

ProtoMS has been written using the GNU Fortran compilers (https://gcc.gnu.org/), on the Linux operating system. ProtoMS is thus known to work well with this compiler and Linux. ProtoMS has also been compiled and tested using the Intel Fortran Compiler. ProtoMS has been compiled with other compilers but not extensively tested. It is therefore advised to use either GNU or Intel compilers with ProtoMS.

You also need an MPI package to perform simulations that require multiple processes, e.g. replica exchange. Such libraries should be available on most modern computers and clusters. The MPI compilers in the GNU package is called mpiff77. ProtoMS has been compiled with both OpenMPI (http://www.open-mpi.org/) and MPICH (https://www.mpich.org/). However, note that this is *not* a requirement any more. ProtoMS will compile without OpenMPI, but you won't be able to run for instance replica exchange.

Building ProtoMS is done with *cmake* (http://www.cmake.org/), thus you need this package installed on your machine. To build ProtoMS type the following in a terminal::

  mkdir build
  cd build
  cmake ..
  make install

and *cmake* will perform the necessary checks before it continues with the installation of ProtoMS. The executable will be placed in the top level of the folder hierarchy.

We recommend to set the environmental variable ``$PROTOMSHOME`` to the installation directory of ProtoMS. This variable is used as a shortcut in the tutorials and by the Python tools. ProtoMS is also able to substitue this variable when it is used in ProtoMS command files.

=================
Requirements
=================

The *MC program* has the following requirements:

* Fortan compiler, GNU (https://gcc.gnu.org/) or Intel recommended
* Python (https://www.python.org/), required to compile and run ProtoMS
* CMake (http://www.cmake.org/), required to compile ProtoMS

Optional:

* MPI, recommended OpenMPI (http://www.open-mpi.org/)  or MPICH (https://www.mpich.org/) - some functions will be unavailable without MPI


The *ProtoMS tools* have the following requirements:

* Python, version 2.7
* NumPy (http://www.numpy.org/)
* SciPy (http://www.scipy.org/)
* Matplotlib (http://www.matplotlib.org/)

Optional:

* AmberTools (http://www.ambermd.org/)          : Required to paramterise small molecules
* pymbar (https://github.com/choderalab/pymbar) : Required to perform MBAR calculations

