*****************************
Compilation and Installation
*****************************

The ProtoMS package supplies the following files and directories;

* **README** File that contains brief installation instructions for ProtoMS, and any last minute addendums or errata that arrived too late to make it into the manual!

* **tutorials** This directory contains a number of examples that demonstrate applications of ProtoMS.

* **interfaces** This directory contains interfaces for ProtoMS to common scripting languages. Currently interfaces to the Perl and Python languages are provided.

* **manual.pdf** This manual!

* **parameter** This directory contains all of the standard parameter files that describe the standard forcefields implemented in ProtoMS.

* **src** This directory contains all of the source code.

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

ProtoMS has been written using the GNU Fortran compiler, g77, version 3.3.4, on the Linux operating system. ProtoMS is thus known to work well with g77 and Linux. ProtoMS has also been compiled and tested using the Portland Group Fortran Compiler. ProtoMS has been compiled with other compilers but not extensively tested. It is therefore advised to use g77 of pgf77 with ProtoMS.

=================
Building ProtoMS
=================

Building ProtoMS should be straightforward if you have a Fortran compiler that supports the extensions described in section 2.1, and a version of ``make`` that supports the GNU Makefile format. Simply go into the ``src`` directory and edit the ``Makefile`` that you find there. This file contains a lot of comments to help you edit the file, and all you should need to do is edit the compilation flags to best optimise ProtoMS to your system. Once you have edited the ``Makefile`` you can then run ``make``. After about 5 minutes, the compilation should hopefully finish, and the ProtoMS executable placed in the top directory. The executable will be called ``protoms3`` on UNIX/Linux, and ``protoms2.exe`` on Windows. This executable should be run from the command line, or via a script. You should then change into the example directory and try out some of the examples. You should also try some of the tests as well to ensure that your version of ProtoMS is working correctly.


