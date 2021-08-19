cmake_minimum_required(VERSION 3.12)

################################################################################
# Define our package information.
################################################################################

set(PACKAGE_VERSION_MAJOR "2")
set(PACKAGE_VERSION_MINOR "0")
set(PACKAGE_VERSION_PATCH "0")
set(PACKAGE_NAME "IntelQS")
set(PACKAGE_DESCRIPTION
    "Intel high performance quantum computing pure state simulator framework."
)
set(PACKAGE_AUTHOR "Intel Corporation")
set(PACKAGE_AUTHOR_EMAIL "AQUA, QAPG and IAGS Groups at Intel Corporation")

project(intel_qs CXX)

# Out-of-source builds only
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_DISABLE_SOURCE_CHANGES  ON)

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt.")
endif()

################################################################################
# Set the location of helper CMake modules in the source tree.
################################################################################

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

################################################################################
# Set some compiler options related to C++.
################################################################################

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

################################################################################
# Options for the build.
################################################################################

option(IqsMPI "Enable MPI?" OFF)
option(IqsMKL "Enable MKL?" OFF)
option(IqsPython "Enable Python wrap?" OFF)
option(IqsUtest "Enable unit test?" OFF)
option(BuildExamples "Build the Examples" OFF)
option(BuildInterface "Build the QASM Interface" OFF)
option(IqsNative "Enable the latest vector instructions?" OFF)

option(EigenSubmodule "Eigen as submodule?" ON)

################################################################################
# Enable MPI, yes or not?
################################################################################

if(IqsMPI)
    message(STATUS "MPI functionality: enabled.")

    if(IqsPython)
        message(STATUS "Python interface with MPI: enabled." )
    endif()

    if(IqsUtest)
        message(STATUS "UnitTest tested for compatibility with MPI.")
    endif()
    
    find_package(MPI REQUIRED) 
    # Locate an MPI compiler and set appropriate flags.
    message(STATUS "@@ before looking for MPI,    CXX = $ENV{CXX}.")
    message(STATUS "@@ before looking for MPI, MPICXX = ${MPI_CXX_COMPILER}.")
    
    message(STATUS "@@ after  looking for MPI,    CXX = $ENV{CXX}.")
    message(STATUS "@@ after  looking for MPI, MPICXX = ${MPI_CXX_COMPILER}.")

    # Check which MPI compiler is available.
    if(MPI_CXX_FOUND)
        message(STATUS "Found a suitable MPI compiler ${MPI_CXX_COMPILER}.")
        message(STATUS "Compiler vendor is [${CMAKE_CXX_COMPILER_ID}]")
        message(STATUS "To run MPI programs, use ${MPIEXEC_EXECUTABLE}")
        message(STATUS "Include path: ${MPI_CXX_INCLUDE_DIRS}")
        message(STATUS "Compile flags: ${MPI_CXX_COMPILE_FLAGS}")
        message(STATUS "Link flags: ${MPI_CXX_LINK_FLAGS}")
        message(STATUS "Libraries: ${MPI_CXX_LIBRARIES}")

        if(MPI_CXX_VERSION VERSION_LESS "3.0.0")
            message(AUTHOR_WARNING "Intel-QS requires that you use MPICH-3.")
            message(FATAL_ERROR "")
        endif()

        # GNU compilers cannot find the mpi.h header and link libraries because these
        # are not usually set in the environment of the shell. 
        if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
            if(MPI_CXX_VERSION VERSION_LESS "3")
                message(AUTHOR_WARNING "Intel-QS with MPI requires at least version 3.1")
                message(FATAL_ERROR "Unable to compile Intel-QS with MPICH - ${MPI_CXX_VERSION}")
            endif()

            message(STATUS "GCC MPI include path ${MPI_CXX_INCLUDE_DIRS}")
            message(STATUS "--- MPI include path ${MPI_CXX_INCLUDE_PATH}")
            message(STATUS "GCC MPI link libraries path ${MPI_CXX_LIBRARIES}")
            message(STATUS "GCC MPI flags ${MPI_CXX_COMPILE_OPTIONS}")
            message(STATUS "--- MPI flags ${MPI_CXX_COMPILE_DEFINITIONS}")

            add_compile_options(${MPI_CXX_COMPILE_OPTIONS})
            add_compile_definitions(${MPI_CXX_COMPILE_DEFINITIONS})
        
        # Intel compilers.
        elseif(CMAKE_CXX_COMPILER_ID MATCHES Intel)
            message(STATUS "Intel compiler has already configured includes and libs.")
        else()
            message(FATAL_ERROR "Currently, unsupported compiler environment.")
        endif()

        include_directories(${MPI_CXX_INCLUDE_DIRS})
        link_libraries("${MPI_CXX_LIBRARIES}")
    else()
            message(FATAL_ERROR "No MPI-compatible compiler has been found.")
    endif()

    # IntelQS pre-processor defintion indicating use of MPI.
    add_compile_definitions(INTELQS_HAS_MPI)

else()
    message(STATUS "MPI functionality: disabled.")
endif()


################################################################################
# Enable native compilation, yes or not?
# This option allows AVX 256 and 512 instructions to be used.
################################################################################
if (IqsNative)
    if(CMAKE_CXX_COMPILER_ID MATCHES Intel)
        message(STATUS "Setting xhost...")
        add_compile_options(-xhost)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES GNU)
        message(STATUS "Setting march=native...")
        add_compile_options(-march=native)
    endif()
endif()


################################################################################
# Enable compiler optimizations to fix GCC performance.
################################################################################
if (CMAKE_CXX_COMPILER_ID MATCHES GNU)
   add_compile_options(-O3)
endif()

################################################################################
# Locate MKL if it is already configured through Intel (mklvars.sh) scripts.
# - Follow the Intel compiler documentation to ensure that all of the
#   appropriate variables are set correctly in the environment.
################################################################################

# The default build with GCC   compiler is without MKL.
# The default build with Intel compiler is with    MKL.
if(CMAKE_CXX_COMPILER_ID MATCHES Intel)
    set(IqsMKL ON)
endif()

if (IqsMKL)
    if (NOT (DEFINED ENV{MKLROOT}))
        message(AUTHOR_WARNING "Intel-QS requires that the Intel MKL libraries are "
                               "present and configured correctly for your system."
                               "\nPlease invoke the following command to "
                               "configure your environment appropriately."
                               "\nbash> source $(YOUR_MKL_INSTALL)/bin/mklvars.sh intel64 ilp64"
               )
        message(FATAL_ERROR "Intel MKL was not found.")
    else()
        message(STATUS "MKLROOT variable detected...configuring MKL builds.")
        find_package(MKL REQUIRED)
#        set(MKL_USE_sdl true)
#        set(MKL_USE_interface "ilp64")
        message(STATUS "MKLROOT value = $ENV{MKLROOT}")
        message(STATUS "MKL definitions = ${MKL_DEFINITIONS}")
        message(STATUS "MKL include dir = ${MKL_INCLUDE_DIR}")
        message(STATUS "MKL libraries = ${MKL_LIBRARIES}")
        add_compile_definitions(USE_MKL)
        add_compile_definitions(MKL_ILP64)
        add_compile_options(${MKL_DEFINITIONS})
        include_directories(${MKL_INCLUDE_DIR})
    endif()
endif() # for the MKL configuration

############################################################################
# Ad-hoc modifications to include the libraries needed by VSL random number generators.
############################################################################

if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    # Ad-hoc modifications to compile and link with g++.
    add_compile_options(-O3)
    add_compile_options(-m64)
#    add_compile_options(-Wall)
    # MKL libraries: optional
    if (IqsMKL)
        message(STATUS "With GNU compilers, the use of MKL needs appropriate flags.")
        link_libraries(-Wl,--start-group $ENV{MKLROOT}/lib/intel64/libmkl_cdft_core.a
                                         $ENV{MKLROOT}/lib/intel64/libmkl_intel_ilp64.a
                                         $ENV{MKLROOT}/lib/intel64/libmkl_gnu_thread.a
                                         $ENV{MKLROOT}/lib/intel64/libmkl_core.a
                                         $ENV{MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_ilp64.a
                       -Wl,--end-group )
    endif()
    link_libraries( -lgomp -lpthread -lm -ldl)
    link_libraries( -lstdc++)

elseif(CMAKE_CXX_COMPILER_ID MATCHES Intel)
    message(STATUS "With Intel compilers, the use of MKL needs appropriate flags.")
    message(STATUS "the value of CXX is $ENV{CXX}")
    link_libraries("${MKL_LIBRARIES}")

    if(NOT IqsMKL)
        message(FATAL_ERROR "MKL is currently required when Intel compilers are used.")
        # This condition can be relaxed, but it requires ad-hoc determination of the
        # libraries to link in case MKL is not used.
    endif()

    if(IqsMPI)
        link_libraries( -Wl,--start-group
                        $ENV{MKLROOT}/lib/intel64/libmkl_cdft_core.a
                        $ENV{MKLROOT}/lib/intel64/libmkl_intel_ilp64.a
                        $ENV{MKLROOT}/lib/intel64/libmkl_intel_thread.a
                        $ENV{MKLROOT}/lib/intel64/libmkl_core.a
                        $ENV{MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_ilp64.a
                        -Wl,--end-group -liomp5 -lpthread -lm -ldl)
    else()
#dynamic link: not working
#        link_libraries( -liomp5 -lpthread -lm -ldl)
#        add_compile_options(-mkl=parallel)
#static link
        link_libraries( -Wl,--start-group
                        $ENV{MKLROOT}/lib/intel64/libmkl_intel_ilp64.a
                        $ENV{MKLROOT}/lib/intel64/libmkl_intel_thread.a
                        $ENV{MKLROOT}/lib/intel64/libmkl_core.a
                        -Wl,--end-group -liomp5 -lpthread -lm -ldl)
    endif()

elseif(CMAKE_CXX_COMPILER_ID MATCHES Clang)
    # Ad-hoc modifications to compile and link with clang++.
    add_compile_options(-O3)
    add_compile_options(-m64)
    if (IqsMKL)
        message(FATAL_ERROR "Currently, the code does not provide MKL support with clang.")
    endif()
    if (IqsMPI)
        message(FATAL_ERROR "Currently, the code does not provide MPI support with clang.")
    endif()
    link_libraries( -lgomp -lpthread -lm -ldl)
    link_libraries( -lstdc++)

else()
    message(STATUS "CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}")
    message(FATAL_ERROR "Currently, unsupported compiler environment.")
endif()

################################################################################
# Select how you would like to use OpenMP in conjunction with MPI.
################################################################################

find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    message(STATUS "OpenMP LIBS = ${OpenMP_CXX_LIBRARIES} FLAGS=${OpenMP_CXX_FLAGS}")
    add_compile_options("${OpenMP_CXX_FLAGS}")
    link_libraries("${OpenMP_CXX_LIBRARIES}")
elseif(CMAKE_CXX_COMPILER_ID MATCHES Clang)
    message(STATUS "TODO: untested OpenMP with clang")
else()
    message(FATAL_ERROR "Unable to locate OpenMP in any standard locations.")
endif()

################################################################################
# Configure standard IntelQS library preprocessor build flags.
################################################################################

add_compile_definitions(USE_MM_MALLOC)

################################################################################
# List the object files that will be contained in the statically linked
# library archive. Also, place the archive in the lib subdir. 
################################################################################

add_subdirectory(src)

if(BuildExamples)
    message(STATUS "Build examples : enabled.")
    add_subdirectory(examples)
else()
    message(STATUS "Build examples : disabled.")
endif()

add_subdirectory(tutorials)
add_subdirectory(benchmarks)

if (BuildInterface)
    add_subdirectory(interface)
endif()


# INFO: when below flag is defined, checks of the form "assert(condition)" are skipped.
#target_compile_definitions(intel_qs NDEBUG)

################################################################################
# Build the suite of unit-tests for IQS.
# Approach suggested in:
# https://github.com/google/googletest/blob/master/googletest/README.md
################################################################################

# Add the unit test contained in /unit_test
# If required, googletest framework is downloaded and configured.
if (IqsUtest)
    message(STATUS "Unit test : enabled.")
    add_subdirectory(unit_test)
else()
    message(STATUS "Unit test : disabled.")
endif()

################################################################################
# If without MPI, build the Python library for the IQS wrapper.
################################################################################

if (NOT IqsPython)
    message(STATUS "Python wrapper: disabled.")
else()
    message(STATUS "Python wrapper: enabled.")
    add_subdirectory(pybind11)
endif()

################################################################################
# Save user's options in 'config.status' file

set(CONFIG_STATUS_CONTENT "User-defined cmake options:\n\n")
STRING(CONCAT CONFIG_STATUS_CONTENT ${CONFIG_STATUS_CONTENT} "IqsMKL    = " ${IqsMKL} "\n")
STRING(CONCAT CONFIG_STATUS_CONTENT ${CONFIG_STATUS_CONTENT} "IqsMPI    = " ${IqsMPI} "\n")
STRING(CONCAT CONFIG_STATUS_CONTENT ${CONFIG_STATUS_CONTENT} "IqsUtest  = " ${IqsUtest} "\n")
STRING(CONCAT CONFIG_STATUS_CONTENT ${CONFIG_STATUS_CONTENT} "IqsPython = " ${IqsPython} "\n")
STRING(CONCAT CONFIG_STATUS_CONTENT ${CONFIG_STATUS_CONTENT} "IqsNative = " ${IqsNative} "\n\n")
STRING(CONCAT CONFIG_STATUS_CONTENT ${CONFIG_STATUS_CONTENT} "BuildExamples  = " ${BuildExamples} "\n")
STRING(CONCAT CONFIG_STATUS_CONTENT ${CONFIG_STATUS_CONTENT} "BuildInterface = " ${BuildInterface})
FILE(WRITE ${CMAKE_CURRENT_BINARY_DIR}/config.status ${CONFIG_STATUS_CONTENT} )

################################################################################
