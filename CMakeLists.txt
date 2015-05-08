# CMake file for ProtoMS
cmake_minimum_required(VERSION 2.8.12.2)

# Project name and language - enables Fortran dependency checking
project(ProtoMS Fortran)
set(prog_name "protoms3")

message("COMPILER: ${CMAKE_Fortran_COMPILER_ID}")

# Compile options
option(STATIC "Compile statically linked executable" OFF)
if(STATIC)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
endif()

# Find python and create dynamic header
find_package(PythonInterp)
if(PYTHONINTERP_FOUND)
    message("Generating compilation info...")
    execute_process(COMMAND python
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/getid.py
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/idme.inc)
else()
    # I don't think the original makefile does anything here
    # Don't see how it compiles without python
endif()

# Set flags for GNU compiler
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m64 -fomit-frame-pointer -ffixed-line-length-none")

# Other compilers - I don't have these installed so don't know the names to check for
# Set flags for g77 compiler
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "g77")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -malign-double -ffixed-line-length-none -Wall")

# Set flags for Intel compiler
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "IFORT")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -132")

# Set flags for Portland Compiler
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "pgf77")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fast -Mextend")
endif()

# Find Fortran MPI libraries - raise error if not found
find_package(MPI REQUIRED)
add_definitions(${MPI_Fortran_COMPILE_FLAGS})
include_directories(${MPI_Fortran_INCLUDE_PATH})
link_directories(${MPI_Fortran_LIBRARIES})

# Collect all .F files in src directory and sort out .inc includes
# Could (probably should) replace this GLOB with an actual list of files
# Dependencies are resolved automatically
file(GLOB SOURCE_FILES src/*.F)
include_directories(src)

# Create target from source files and link to MPI library
add_executable (${prog_name} ${SOURCE_FILES})
target_link_libraries(${prog_name} ${MPI_Fortran_LIBRARIES})
