*****************************
Compilation and Installation
*****************************

The ProtoMS package supplies the following files and directories;

* **data** This directory contains a number of useful files, e.g. pre-equilibrated boxes and some template files

* **doc** This directory contains documentation

* **README** File that contains brief installation instructions for ProtoMS, and any last minute addendums or errata that arrived too late to make it into the manual!

* **parameter** This directory contains all of the standard parameter files that describe the standard forcefields implemented in ProtoMS.

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

The tools are written in python and should be compatible with the standard implementation version 2.7. It is known that is does not work with 2.6 and it has not been tested with version 3.0 or more.

=================
Building ProtoMS
=================

ProtoMS has been written using the GNU Fortran compilers, on the Linux operating system. ProtoMS is thus known to work well with this compiler and Linux. ProtoMS has also been compiled and tested using the Intel Fortran Compiler. ProtoMS has been compiled with other compilers but not extensively tested. It is therefore advised to use either GNU or Intel compilers with ProtoMS.

You also need the OpenMPI package as ProtoMS contains instructions that operates on multiple processes. Such libraries should be available on most modern computers and clusters. The MPI compilers in the GNU package is called mpiff77

Building ProtoMS should be straightforward if you have a Fortran compiler that supports the extensions described in section 2.1, and a version of ``make`` that supports the GNU Makefile format. Simply go into the ``src`` directory and edit the ``Makefile`` that you find there. This file contains a lot of comments to help you edit the file, and all you should need to do is edit the compilation flags to best optimise ProtoMS to your system. Once you have edited the ``Makefile`` you can then run ``make``. After about a minute, the compilation should hopefully finish, and the ProtoMS executable placed in the top directory. The executable will be called ``protoms3`` on UNIX/Linux. This executable should be run from the command line, or via a script.


