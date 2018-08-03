*****************************
Compilation and Installation
*****************************


The ProtoMS package supplies the following files and directories;


* **CMakeLists.txt** This file configures ProtoMS prior to compiling.


* **data** This directory contains a number of useful files, e.g. pre-equilibrated boxes and some template files


* **doc** This directory contains documentation


* **README** File that contains brief installation instructions for ProtoMS, and any last minute addenda or errata that arrived too late to make it into the manual!


* **parameter** This directory contains all of the standard parameter files that describe the standard force fields implemented in ProtoMS.


* **protoms.py** A tool to setup common ProtoMS calculations


* **python** This directory contains the python library component of ProtoMS


* **src** This directory contains all of the source code for the main program


* **tools** This directory contains numerous useful scripts to setup and analyse ProtoMS simulations.


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


* The ``Date`` and ``Time`` Fortran 90 intrinsic subroutine is used to get the current time. This is used to provide a default seed to the random number generator. This can be removed by commenting out the relevant lines in ``getoptions.F``, though you will need to provide a random number seed manually.


=================
Requirements
=================


The *MC program* has the following requirements:


* Fortan compiler, GNU (https://gcc.gnu.org/) or Intel recommended

* Python (https://www.python.org/), required to compile and run ProtoMS

* CMake (http://www.cmake.org/), required to compile ProtoMS


Optional:


* MPI, recommended OpenMPI (http://www.open-mpi.org/) or MPICH (https://www.mpich.org/) - some functions will be unavailable without MPI



The *ProtoMS tools* have the following requirements:


* Python, version 2.7, version 3.5 or newer

* NumPy (http://www.numpy.org/)

* SciPy (http://www.scipy.org/)

* Matplotlib (http://www.matplotlib.org/)

* pymbar (https://github.com/choderalab/pymbar)


Optional:


* AmberTools (http://www.ambermd.org/)          : Required to parameterise small molecules

The tools are written in python (https://www.python.org/) and are compatible with version 2.7 as well as version 3.5, 3.6 and 3.7.


A docker image that contains all of the installation dependencies is available. See :ref:`using_docker` below.


==================
Installing ProtoMS
==================


ProtoMS has been written using the GNU Fortran compilers (https://gcc.gnu.org/), on the Linux operating system. ProtoMS is thus known to work well with this compiler and Linux. ProtoMS has also been compiled and tested using the Intel Fortran Compiler. ProtoMS has been compiled with other compilers but not extensively tested. It is therefore strongly advised to use either GNU or Intel compilers with ProtoMS.


You also need an MPI package to perform simulations that require multiple processes, e.g. replica exchange. Such libraries should be available on most modern computers and clusters. The MPI compilers in the GNU package is called mpiff77. ProtoMS has been compiled with both OpenMPI (http://www.open-mpi.org/) and MPICH (https://www.mpich.org/). However, note that this is *not* a requirement any more. ProtoMS will compile without OpenMPI, but you won't be able to run for instance replica exchange.


Building ProtoMS is done with *cmake* http://www.cmake.org/, thus you need this package installed on your machine. Before proceeding it is important that your environment is properly configured. In particular, since version 3.4, ProtoMS installs the Python package *protomslib* into your python environment. If you are using a virtual environment this must be activated so that cmake can locate the correct python interpreter. To prepare the build, type the following in a terminal from the root directory of the code::


  mkdir build

  cd build

  cmake ..


At this point you should check the output of cmake. Unless you're expecting it not to, cmake should have found an appropriate Fortran compiler,  MPI library and Python interpreter. Check that the paths and versions of these correspond to those you expect. If they do not, see :ref:`custom_build` for details on how to customise these. Also note that if cmake has found the system python interpreter (usually /usr/bin/python) it will attempt to install protomslib into a system location requiring root access. Again :ref:`custom_build` covers how to change the Python installation target. If you're happy with what cmake has found then type::


  make install


and *cmake* will perform the necessary checks before it continues with the installation of ProtoMS. The executable will be placed in the top level of the folder hierarchy.


In order for ProtoMS to find the relevant parameter files it is necessary to set the environmental variable ``$PROTOMSHOME`` to the installation directory of ProtoMS. This variable is used as a shortcut in the tutorials and by the Python tools. ProtoMS is also able to substitute this variable when it is used in ProtoMS command files.


Once building is complete it is highly recommended to run the test suite that comes with ProtoMS to test that the build was successful. From the build directory created above simply type::


  ctest -V


All tests should be expected to pass and the above command will provide detailed output. The most common reason for failures is the need to set the correct environment variables. Notably ``$AMBERHOME`` for the setup tests and ``$PROTOMSHOME``, as described above. Another reason for occasional failures is slight formatting and rounding differences between compilers, this can lead to values differing at the final decimal place in results files and such failures can be safely ignored.


.. _custom_build:


======================
Customising the Build
======================


**Tips on using cmake**


The job of cmake is to attempt to locate all of the necessary dependencies for the installation and create a Makefile that will compile ProtoMS. It searches your system for the required components and sets a number of internal variables that store their locations. After being run cmake stores its output in the build directory in a file called CMakeCache.txt. This can be useful after the fact to check which dependencies were found but equally if being run subsequently cmake will prefer to use cached values instead of updating dependencies. For this reason it can be a good idea to delete CMakeCache.txt if you find you need to run cmake more than once or cmake does not appear to be behaving as expected.


**Manually specify cmake variables**


The locations that cmake will search for dependencies are quite comprehensive, however they are also dependent on the system in use and the value of current environment variables. Thus cmake may not be able to find the required libraries even if they're present in your system or may find the wrong versions. To coerce cmake into finding the relevant dependencies you can try:


 1. Setting environment variables - The $PATH environment variable is checked by cmake for relevant executables e.g. gfortran, mpirun. Prepending to or rearranging entries in the PATH makes dependencies discoverable by cmake. The FC environment variable is a standard method for manually specifying the Fortran compiler.

 2. Manually setting cmake variables - Whilst cmake attempts to automatically discover correct values for dependencies you may find that setting them manually is easier. This can be performed interactively using the ``ccmake`` utility. If you execute ``ccmake ..`` from the build directory you will be presented with a interface showing the current value of cmake variables. Press ``t`` to see more values. You can edit values from this menu before pressing ``c`` to configure (any problems should be flagged by cmake here) and ``g`` to generate a new Makefile and exit.

 3. Manually setting cmake variables on the command line - If you prefer the value of any cmake variable can be specified directly from the command line. The ``-D`` flag to cmake can be used repeatedly for this purpose. For instance - ``cmake -DCMAKE_Fortran_COMPILER=gfortran ..`` - sets the value of the  variable CMAKE_Fortran_COMPILER to gfortran. You can use ``ccmake`` to determine the names of variables to set.

  

**Installation of protomslib**


You can customise the installation of the python library component by specify a value for the cmake variable PYTHON_INSTALL_OPTIONS (see above). The value of this variable will be appended like so to the command below which is executed by cmake::


  python setup.py install $PYTHON_INSTALL_OPTIONS


To see the available options you can run::


  python $PROTOMSHOME/python/setup.py install --help


The most frequently useful options are ``--user``, that requests an installation into ``$HOME/.local``, and --prefix that allows an installation root directory to be specified manually.


.. _using_docker:


============
Using Docker
============


Version 3.4 of ProtoMS is also available via docker. Downloading and running the image can be accomplished easily with the command::


  docker run -it jessexgroup/protoms:3.4


The image is based on the python:3.6-slim image with additional installation of the relevant python dependencies as well as amber tools 18. To construct your own docker images from scratch see ``Dockerfile_test`` and ``Dockerfile`` in the root ProtoMS directory and the instructions therein. This will allow you to use newer versions of the dependencies than are available via the public image.
