# CMake file for ProtoMS
cmake_minimum_required(VERSION 2.8.10.2)

# Project name and language - enables Fortran dependency checking
project(ProtoMS Fortran)
set(prog_name "protoms3")
set(install_dir "protoms3")
set(ProtoMS_VERSION_MAJOR 3)
set(ProtoMS_VERSION_MINOR 0)

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
    set(CMAKE_Fortran_FLAGS_Debug "-Og -g")
endif()

if(CMAKE_BUILD_TYPE MATCHES "Release")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_Release}")
endif()
if(CMAKE_BUILD_TYPE MATCHES "Debug")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_Debug}")
endif()

message("-- COMPILER: ${CMAKE_Fortran_COMPILER_ID}")

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
find_package(PythonInterp REQUIRED)
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
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m64 -fomit-frame-pointer -ffixed-line-length-none")

# Set flags for Intel compiler
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
    message("-- Using Intel compiler flags")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -132")

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
if(MPI_Fortran_FOUND)
    target_link_libraries(${prog_name} ${MPI_Fortran_LIBRARIES})
endif()

# Setup install directory
# Optional 'make install' step
install(TARGETS protoms3
        RUNTIME DESTINATION ${install_dir})
install(PROGRAMS "protoms.py"
        DESTINATION ${install_dir})
install(DIRECTORY data
        DESTINATION ${install_dir})
install(DIRECTORY doc
        DESTINATION ${install_dir})
install(DIRECTORY parameter
        DESTINATION ${install_dir})
install(DIRECTORY tools
        DESTINATION ${install_dir})

# Add testing
# Needs actual unit/integration tests here
add_test(TestPrintVersion protoms3 -v)

enable_testing()
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
                  DEPENDS protoms3)
