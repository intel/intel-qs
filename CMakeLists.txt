cmake_minimum_required(VERSION 3.10)

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

project(IQS)

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
option(BuildExamples "Build the Examples" ON)

################################################################################
# Enable MPI, yes or not?
################################################################################

if(IqsMPI)
    message(STATUS "MPI functionality: enabled.")

    if(IqsPython)
        message(STATUS "Python not compatible with MPI: disabled.")
        set(IqsPython OFF)
    endif()

    if(IqsUtest)
        message(STATUS "UnitTest tested for compatibility with MPI.")
#        message(STATUS "UnitTest not compatible with MPI: disabled.")
#        set(IqsUtest OFF)
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
else()
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
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY lib)

add_subdirectory(src)

if(BuildExamples)
  add_subdirectory(examples)
endif()


#include_directories(${PROJECT_SOURCE_DIR}/include)

#add_library(intel_qs SHARED
#    src/interface_api_memory.cpp
#    src/interface_api_qasm.cpp
#    src/interface_api_qubitid.cpp
#    src/interface_api_version.cpp
#    src/highperfkernels.cpp
#    src/qureg_apply1qubitgate.cpp
#    src/qureg_apply2qubitgate.cpp
#    src/qureg_applyctrl1qubitgate.cpp
#    src/qureg_applydiag.cpp
#    src/qureg_applyswap.cpp
#    src/qureg_applytoffoli.cpp
#    src/qureg_expectval.cpp
#    src/qureg_fusion.cpp
#    src/qureg_init.cpp
#    src/qureg_measure.cpp
#    src/qureg_noisysimul.cpp
#    src/qureg_permute.cpp
#    src/qureg_utils.cpp
#    src/qureg_version.cpp
#    src/mpi_env.cpp
#    src/mpi_exception.cpp
#    src/utils.cpp
#    src/rng_utils.cpp
#)

#set_target_properties(intel_qs PROPERTIES POSITION_INDEPENDENT_CODE 1)


#target_include_directories(intel_qs INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

# INFO: when below flag is defined, checks of the form "assert(condition)" are skipped.
#target_compile_definitions(intel_qs NDEBUG)

################################################################################
# Build all of the executables (.exe) contained in the subfolder.
################################################################################

#set(EXECUTABLE_OUTPUT_PATH bin)

################################################################################
# Build the executables (.exe) contained in the unit_test/ subfolder.
# No googletest framework.
################################################################################

#add_executable(mpi_env_test.exe unit_test/mpi_env_test.cpp
#               src/mpi_env.cpp
#               src/mpi_exception.cpp)

#add_executable(mpi_speedup_test.exe unit_test/mpi_speedup_test.cpp
#               src/mpi_env.cpp
#               src/mpi_exception.cpp)
#target_link_libraries(mpi_speedup_test.exe intel_qs)

#add_executable(openmp_test.exe unit_test/openmp_test.cpp)

#target_link_libraries(mpi_test.exe intel_qs)

################################################################################
# Build the executables (.exe) contained in the examples/ subfolder.
################################################################################

#if(BuildExamples)
#    add_executable(get_started_with_IQS.exe tutorials/get_started_with_IQS.cpp)
#    target_link_libraries(get_started_with_IQS.exe iqs)

#    add_executable(get_started_with_noisy_IQS.exe tutorials/get_started_with_noisy_IQS.cpp)
#    target_link_libraries(get_started_with_noisy_IQS.exe intel_qs)

#    add_executable(benchgates.exe examples/benchgates.cpp)
#    target_link_libraries(benchgates.exe intel_qs)

#    add_executable(quantum_fourier_transform.exe  examples/quantum_fourier_transform.cpp)
#    target_link_libraries(quantum_fourier_transform.exe intel_qs)

#    add_executable(expect_value_test.exe examples/expect_value_test.cpp)
#    target_link_libraries(expect_value_test.exe intel_qs)

#    add_executable(heisenberg_dynamics.exe examples/heisenberg_dynamics.cpp)
#    target_link_libraries(heisenberg_dynamics.exe intel_qs)

#    add_executable(noisy_circuit_test.exe examples/noisy_circuit_test.cpp)
#    target_link_libraries(noisy_circuit_test.exe intel_qs)

#    add_executable(test_of_custom_gates.exe examples/test_of_custom_gates.cpp)
#    target_link_libraries(test_of_custom_gates.exe intel_qs)

#    add_executable(basic_code_for_scaling.exe benchmarks/basic_code_for_scaling.cpp)
#    target_link_libraries(basic_code_for_scaling.exe intel_qs)
#endif()

################################################################################
# Build the QASM interface executables (.exe) contained in the src/ folder.
################################################################################

#if(BuildInterface)
#    add_executable(intelqs_interface_main.exe src/intelqs_interface_main.cpp)
#    target_link_libraries(intelqs_interface_main.exe intel_qs)
#endif()

################################################################################
# Build the suite of unit-tests for IQS.
# Approach suggested in:
# https://github.com/google/googletest/blob/master/googletest/README.md
################################################################################

# Only download GoogleTest if asked to do so for unit testing purposes.
#if (IqsUtest)

    # Download and unpack googletest at configure time
#    configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
#    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
#      RESULT_VARIABLE result
#      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
#    if(result)
#      message(FATAL_ERROR "CMake step for googletest failed: ${result}")
#    endif()
#    execute_process(COMMAND ${CMAKE_COMMAND} --build .
#      RESULT_VARIABLE result
#      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
#    if(result)
#      message(FATAL_ERROR "Build step for googletest failed: ${result}")
#    endif()

    # Prevent overriding the parent project's compiler/linker
    # settings on Windows
#    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

    # Add googletest directly to our build. This defines
    # the gtest and gtest_main targets.
#    add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
#                     ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
#                     EXCLUDE_FROM_ALL)

    # The gtest/gtest_main targets carry header search path
    # dependencies automatically when using CMake 2.8.11 or
    # later. Otherwise we have to add them here ourselves.
#    if (CMAKE_VERSION VERSION_LESS 2.8.11)
#      include_directories("${gtest_SOURCE_DIR}/include")
#    endif()

#endif()

# Now simply link against gtest or gtest_main as needed. Eg
if (NOT IqsUtest)
    message(STATUS "Unit test: disabled.")
else()
    message(STATUS "Unit test: enabled.")
    add_executable(utest unit_test/suite_of_tests.cpp
                   src/qaoa_features.cpp)
    target_link_libraries(utest gtest_main intel_qs)
    #add_test(NAME utest_test COMMAND utest)
endif()


################################################################################
# If without MPI, build the Python library for the IQS wrapper.
################################################################################

if (NOT IqsPython)
    message(STATUS "Python wrapper: disabled.")
else()
    message(STATUS "Python wrapper: enabled.")
#    set(EXECUTABLE_OUTPUT_PATH lib)
#    add_executable(intelqs_py`python3-config --extension-suffix` pybind11/intelqs_py.cpp)
#    target_link_libraries(intelqs_py`python3-config --extension-suffix` intel_qs)
    find_package(pybind11 REQUIRED)
    pybind11_add_module(intelqs_py SHARED pybind11/intelqs_py.cpp
                        src/qaoa_features.cpp)
    target_link_libraries(intelqs_py PRIVATE intel_qs)
#    target_compile_definitions(intelqs_py NDEBUG)
#    add_link_options(-shared)
endif()
