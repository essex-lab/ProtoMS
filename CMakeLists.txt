# CMake file for ProtoMS

# ---------------------------------------------- #
# Command line options:
# 
# MPI=YES                       ! defaults to NO
#
# ---------------------------------------------- #

cmake_minimum_required(VERSION 2.8.7)

# Project name and language - enables Fortran dependency checking
set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR} CACHE PATH "Setting install dir to the src dir")
project(ProtoMS Fortran)
set(prog_name "protoms3")
set(ProtoMS_VERSION_MAJOR 3)
set(ProtoMS_VERSION_MINOR 2)

# Check build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose build type: Debug, Release"
        FORCE)
endif()

if(NOT CMAKE_Fortran_FLAGS_Release)
    set(CMAKE_Fortran_FLAGS_Release "-O3")
endif()
if(NOT CMAKE_Fortran_FLAGS_Debug)
    set(CMAKE_Fortran_FLAGS_Debug "-g")
endif()

if(CMAKE_BUILD_TYPE MATCHES "Release")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_Release}")
endif()
if(CMAKE_BUILD_TYPE MATCHES "Debug")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_Debug}")
endif()

message("-- COMPILER: ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")

# Find MPI libraries

option(MPI "Use MPI for tools requiring parallel runs" NO)
option(FIND_MPI "Attempt to automatically find MPI information" YES)
option(MPI_USE_MPIEXEC "Use CMake's FindMPI to find mpiexec" NO)
option(OPENMP "Use OPENMP for parallel runs" NO)

if(MPI)
   add_definitions(-DUSE_MPI)
   if(FIND_MPI)
      find_package(MPI REQUIRED)
      include_directories(${MPI_INCLUDE_PATH})
      set(COMPILER_WRAPPER mpif77)
   endif()
   message(STATUS "MPI enabled")
endif()

# Find Fortran MPI libraries
option(AVOID_MPI "Don't use MPI even if available" OFF)
find_package(MPI)
if(MPI_Fortran_FOUND AND NOT AVOID_MPI)
    message("-- Using MPI")
    add_definitions(${MPI_Fortran_COMPILE_FLAGS})
    include_directories(${MPI_Fortran_INCLUDE_PATH})
    link_directories(${MPI_Fortran_LIBRARIES})

    # Define WITH_MPI_ for C preprocessor - conditional compilation
    add_definitions(-DWITH_MPI_)
else()
    message("-- Not using MPI")
endif()


# Find Python and create dynamic include files
find_package(PythonInterp 2.7 REQUIRED)
message("-- Generating compilation info...")
execute_process(COMMAND python
                ${CMAKE_CURRENT_SOURCE_DIR}/src/getid.py
                ${CMAKE_CURRENT_SOURCE_DIR}/src/idme.inc)
execute_process(COMMAND python
                ${CMAKE_CURRENT_SOURCE_DIR}/src/getsep.py
                ${CMAKE_CURRENT_SOURCE_DIR}/src/sep.inc)

# Set flags for GNU compiler
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    message("-- Using GNU compiler flags")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffixed-line-length-132 -ffixed-form -fimplicit-none")

# Set flags for Intel compiler
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
    message("-- Using Intel compiler flags")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -132")

# Set flags for IBM compiler
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "XL")
    message("-- Using IBM compiler flags")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qfixed=132")

# Compiler not recognised
else()
    message("-- Your compiler is not supported and has not been tested, please consider using
             GNU gfortran or the Intel Fortran compiler.")
    message("-- However, you may wish to attempt compilation by manually providing the correct
             flag to enable a line length of 132 in the CMake variable CMAKE_Fortran_FLAGS.")
endif()

message("    ${CMAKE_Fortran_FLAGS}")

# Collect all .F files in src directory and sort out .inc includes
# Could replace this GLOB with an actual list of files
# Dependencies are resolved automatically
file(GLOB SOURCE_FILES src/*.F)
include_directories(src)

# Exclude files that require MPI if it's not available
if(NOT MPI_Fortran_FOUND OR AVOID_MPI)
    message("-- Excluding files requiring MPI")
    message("-- Replica exchange will not be available")
    set(FILES_REQUIRING_MPI
        "setupreplicas.F"
        "lambdare.F"
        "assignlambdare.F"
        "assigntemperaturere.F"
        "temperaturere.F")
    foreach(MPI_FILE IN ITEMS ${FILES_REQUIRING_MPI})
        list(REMOVE_ITEM SOURCE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/${MPI_FILE}")
    endforeach(MPI_FILE)
endif()

# Create target from source files and link to MPI library if found
add_executable(${prog_name} ${SOURCE_FILES})
if(MPI_Fortran_FOUND AND NOT AVOID_MPI)
    target_link_libraries(${prog_name} ${MPI_Fortran_LIBRARIES})
endif()

# Setup install directory
# Optional 'make install' step
install(TARGETS protoms3
        RUNTIME DESTINATION "./")

# Add testing
# Needs actual unit/integration tests here
add_test(TestPrintVersion protoms3 -v)

###################################
# Python Unit tests
###################################
set(PYTHON_EXECUTABLE ${PYTHON_EXECUTABLE} ${PYTHON_DASH_B})

add_test(test_install_dependencies ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_install_dependencies.py)
add_test(test_ligand_setup ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_ligand_setup.py)
add_test(test_mpi ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_mpi_install.py)
add_test(test_parameters_ff ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_parameters_ff.py)
add_test(test_path ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_path.py)
add_test(test_setup_protoms ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_setup_protoms.py)
add_test(test_tools_protoms ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_tools_protoms.py)

add_test(test_prot_setup ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_setup/test_prot_setup.py)
add_test(test_equil_setup ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_equil/test_equil_prot.py)
add_test(test_sampling_prot ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_sampling/test_sampling_prot.py)
add_test(test_gcmc_sim ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_gcmc/test_gcmc_sim.py)
add_test(test_jaws1_sim ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_jaws1/test_jaws1_sim.py)
add_test(test_jaws2_sim ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_jaws2/test_jaws2_sim.py)

enable_testing()
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
                  DEPENDS protoms3)
