# CMake file for ProtoMS

cmake_minimum_required(VERSION 2.8.7)

# Project name and language - enables Fortran dependency checking
set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR} CACHE PATH "Setting install dir to the src dir")
project(ProtoMS Fortran)
set(PROG_NAME "protoms3")
set(ProtoMS_VERSION_MAJOR 3)
set(ProtoMS_VERSION_MINOR 3)

# Check build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose build type: Debug, Release"
        FORCE)
endif()

message("-- COMPILER: ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")

# CMake options
option(AVOID_MPI "Don't use MPI even if available" OFF)
option(COVERAGE "Prepare for coverage tests" OFF)

# Find Fortran MPI libraries
find_package(MPI)
if(NOT MPI_Fortran_FOUND AND NOT AVOID_MPI AND CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" AND EXISTS "/local/software/openmpi/1.6.4/intel/lib/libmpi_f90.so")
    # Attempt fix for Southampton RedHat configuration
    message("-- Attempting fix for Southampton RedHat configuration")
    set(MPIEXEC "/local/software/openmpi/1.6.4/intel/bin/mpiexec")
    set(MPI_Fortran_COMPILER "/local/software/openmpi/1.6.4/intel/bin/mpif90")
    set(MPI_Fortran_INCLUDE_PATH "/local/software/openmpi/1.6.4/intel/include")
    set(MPI_Fortran_LIBRARIES "/local/software/openmpi/1.6.4/intel/lib/libmpi_f90.so")
    find_package(MPI)
    if(MPI_Fortran_FOUND)
        message("-- Fix successful")
    else()
        message("-- Fix unsuccessful")
    endif()
endif()

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
execute_process(COMMAND ${PYTHON_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/src/getid.py
                ${CMAKE_CURRENT_SOURCE_DIR}/src/idme.inc)
execute_process(COMMAND ${PYTHON_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/src/getsep.py
                ${CMAKE_CURRENT_SOURCE_DIR}/src/sep.inc)

# Set flags for GNU compiler
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    message("-- Using GNU compiler flags")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffixed-line-length-132 -ffixed-form -fimplicit-none")
    if(COVERAGE)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fprofile-arcs -ftest-coverage")
    endif()
    if(CMAKE_BUILD_TYPE MATCHES "Release")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -fgcse-lm -fgcse-sm -fgcse-las -fsched-spec-load")
    elseif(CMAKE_BUILD_TYPE MATCHES "Debug")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g")
    endif()

# Set flags for Intel compiler
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
    message("-- Using Intel compiler flags")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -132")
    if(COVERAGE)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -prof-gen=srcpos")
    endif()
    if(CMAKE_BUILD_TYPE MATCHES "Release")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3")
    elseif(CMAKE_BUILD_TYPE MATCHES "Debug")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g")
    endif()

# Set flags for IBM compiler
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "XL")
    message("-- Using IBM compiler flags")
    message("-- This compiler is not supported but has been tested in the past - proceed with caution")
    if(COVERAGE)
        message("-- Coverage testing not supported using IBM compiler, please add coverage flags yourself")
    endif()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qfixed=132")
    if(CMAKE_BUILD_TYPE MATCHES "Release")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3")
    elseif(CMAKE_BUILD_TYPE MATCHES "Debug")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g")
    endif()

# Set flags for Portland compiler
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "PGI")
    message("-- Using PGI compiler flags")
    message("-- This compiler is not supported but has been tested in the past - proceed with caution")
    if(COVERAGE)
        message("-- Coverage testing not supported using Portland compiler, please add coverage flags yourself")
    endif()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mextend")
    if(CMAKE_BUILD_TYPE MATCHES "Release")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2")
    elseif(CMAKE_BUILD_TYPE MATCHES "Debug")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g")
    endif()
    
# Compiler not recognised
else()
    message("-- Your compiler is not supported and has not been tested, please consider using GNU gfortran or the Intel Fortran compiler.")
    message("-- However, you may wish to attempt compilation by manually providing the correct flag to enable a line length of 132 in the CMake variable CMAKE_Fortran_FLAGS.")
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
        "assignlambdare.F"
        "assigntemperaturere.F"
        "assignadamre.F"
        "lambdare.F"
        "gcmcre.F"
        "temperaturere.F")
    foreach(MPI_FILE IN ITEMS ${FILES_REQUIRING_MPI})
        list(REMOVE_ITEM SOURCE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/${MPI_FILE}")
    endforeach(MPI_FILE)
endif()

# Create target from source files and link to MPI library if found
add_executable(${PROG_NAME} ${SOURCE_FILES})
if(MPI_Fortran_FOUND AND NOT AVOID_MPI)
    target_link_libraries(${PROG_NAME} ${MPI_Fortran_LIBRARIES})
endif()

# Setup install directory
# Optional 'make install' step
install(TARGETS protoms3
        RUNTIME DESTINATION "./")
              
if(COVERAGE)
    #Enabling creation and submission of dashboard
    ##configure_file("${CMAKE_CURRENT_SOURCE_DIR}/CTestConfig.cmake" 
    #               @only
    #              )

    include(CTest)
    include(CTestConfig.cmake)
    message("-- Doing coverage testing")
endif()

###################################
# Python Unit tests
###################################
add_test(NAME TestPrintVersion COMMAND protoms3 -v)

set(PYTHON_EXECUTABLE ${PYTHON_EXECUTABLE} ${PYTHON_DASH_B})

add_test(NAME test_install_dependencies 
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_install_dependencies.py)

add_test(NAME test_ligand_setup 
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_ligand_setup.py)

add_test(NAME test_parameters_ff 
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_parameters_ff.py)

add_test(NAME test_path 
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_path.py)

add_test(NAME test_tools_protoms 
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_tools_protoms.py)

add_test(NAME ProtSetupTest
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_prot_setup.py)

add_test(NAME EquilSetupTest
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_equil.py EquilSetupTest)
add_test(NAME EquilSimulationTest
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_equil.py EquilSimulationTest)

add_test(NAME EnergiesSimulationTip3pTest
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_energies.py EnergiesSimulationTip3pTest)
add_test(NAME EnergiesSimulationTip4pTest
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_energies.py EnergiesSimulationTip4pTest)
add_test(NAME EnergiesSimulationSoftSolventTest
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_soft_energies.py EnergiesSimulationSoftSolventTest)
add_test(NAME EnergiesSimulationSoftProteinTest
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_soft_energies.py EnergiesSimulationSoftProteinTest)
add_test(NAME EnergiesSimulationSoftGcsoluteTest
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_soft_energies.py EnergiesSimulationSoftGcsoluteTest)
add_test(NAME EnergiesSimulationSoftSoluteTest
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_soft_energies.py EnergiesSimulationSoftSoluteTest)

add_test(NAME SamplingSetupTest
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_sampling.py SamplingSetupTest)
add_test(NAME SamplingSimulationTest
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_sampling.py SamplingSimulationTest)

add_test(NAME Jaws1SetupTest
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_jaws1.py Jaws1SetupTest)
add_test(NAME Jaws1SimulationTest
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_jaws1.py Jaws1SimulationTest)

  add_test(NAME FreeEnergyAnalysis
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_free_energy_analysis.py)
  
# Skip MPI tests is not using MPI
if(MPI_Fortran_FOUND AND NOT AVOID_MPI)
    add_test(NAME test_mpi
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_mpi_install.py)

    add_test(NAME GcmcSetupBoxTest
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_gcmc.py GcmcSetupBoxTest)
    add_test(NAME GcmcSetupTest
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_gcmc.py GcmcSetupTest)
    add_test(NAME GcmcSimulationTest
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_gcmc.py GcmcSimulationTest)

    add_test(NAME Jaws2SetupTest
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_jaws2.py Jaws2SetupTest)
    add_test(NAME Jaws2SimulationTest
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_jaws2.py Jaws2SimulationTest)

    add_test(NAME RetiSnglSetupTest
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_reti_sngl.py RetiSnglSetupTest)
    add_test(NAME RetiSnglSimulationFreeTest
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_reti_sngl.py RetiSnglSimulationFreeTest)
    add_test(NAME RetiSnglSimulationGasTest
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_reti_sngl.py RetiSnglSimulationGasTest)

    add_test(NAME RetiDblSetupTest
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_reti_dbl.py RetiDblSetupTest)
    add_test(NAME RetiDblSimulationTest
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_reti_dbl.py RetiDblSimulationTest)
endif()

enable_testing()
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
                  DEPENDS protoms3)
              add_custom_target(pytest COMMAND "py.test ${CMAKE_CURRENT_SOURCE_DIR}/tests/")

add_custom_target(clean-test COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/clean-test.cmake)

if(COVERAGE)
    if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
        add_custom_target(coverage COMMAND ${CMAKE_CTEST_COMMAND} -D Experimental)
    elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
        message("-- Coverage not yet automated with Intel compiler, see testing readme.")
    endif()
endif()
