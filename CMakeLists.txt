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
        "Choose build type: None, Debug, Release"
        FORCE)
endif()

if(NOT CMAKE_Fortran_FLAGS_RELEASE)
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
endif()
if(NOT CMAKE_Fortran_FLAGS_DEBUG)
    set(CMAKE_Fortran_FLAGS_RELEASE "-Og -g")
endif()

if(CMAKE_BUILD_TYPE MATCHES "Release")
    set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS_RELEASE} ${CMAKE_Fortran_FLAGS})
endif()
if(CMAKE_BUILD_TYPE MATCHES "Debug")
    set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS_DEBUG} ${CMAKE_Fortran_FLAGS})
endif()

message("COMPILER: ${CMAKE_Fortran_COMPILER_ID}")

# Compile options
option(STATIC "Compile statically linked executable" OFF)
if(STATIC)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
endif()

# Find Fortran MPI libraries - raise error if not found
find_package(MPI REQUIRED)
add_definitions(${MPI_Fortran_COMPILE_FLAGS})
include_directories(${MPI_Fortran_INCLUDE_PATH})
link_directories(${MPI_Fortran_LIBRARIES})


# Find python and create dynamic include file
find_package(PythonInterp REQUIRED)
message("Generating compilation info...")
execute_process(COMMAND python
                ${CMAKE_CURRENT_SOURCE_DIR}/src/getid.py
                ${CMAKE_CURRENT_SOURCE_DIR}/src/idme.inc)

# Set flags for GNU compiler
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    message("Using GNU compiler flags")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m64 -fomit-frame-pointer -ffixed-line-length-none")

# Set flags for Intel compiler
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
    message("Using Intel compiler flags")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -132")

else()
    message("Your compiler is not supported and has not been tested, please consider using
             GNU gfortran or the Intel Fortran compiler.")
    message("However, you may wish to attempt compilation by manually providing the correct
             flag to enable a line length of 132 in the CMake variable CMAKE_Fortran_FLAGS.")
endif()

message("Using flags:")
message("${CMAKE_Fortran_FLAGS}")

# Collect all .F files in src directory and sort out .inc includes
# Could (probably should) replace this GLOB with an actual list of files
# Dependencies are resolved automatically
file(GLOB SOURCE_FILES src/*.F)
include_directories(src)

# Create target from source files and link to MPI library
add_executable (${prog_name} ${SOURCE_FILES})
target_link_libraries(${prog_name} ${MPI_Fortran_LIBRARIES})

# Get documentation and tools
file(GLOB TOOL_FILES tools/*)
file(GLOB_RECURSE DOC_FILES doc/*)

# Setup install directory
install(TARGETS protoms3
        RUNTIME DESTINATION ${install_dir}/bin)
install(PROGRAMS "protoms.py"
        DESTINATION ${install_dir}/bin)
install(PROGRAMS ${TOOL_FILES}
        DESTINATION ${install_dir}/bin)
install(FILES ${DOC_FILES}
        DESTINATION ${install_dir}/doc)
