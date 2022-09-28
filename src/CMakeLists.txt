################################################################################
# Construction of IQS as a shared library
################################################################################

if (IqsNoise)
    # Download and unpack eigen library at configure time.
    cmake_minimum_required(VERSION 3.12.0)
    include(ExternalProject)
    ExternalProject_Add(deps-eigen
      GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
      PREFIX            "${CMAKE_BINARY_DIR}/deps/eigen-prefix/"
      CONFIGURE_COMMAND ""
      BUILD_COMMAND     ""
      INSTALL_COMMAND   ""
      GIT_SHALLOW 1
    )
    message(STATUS "Eigen library is downloaded in: ${PROJECT_BINARY_DIR}/deps/eigen-prefix/src/eigen/")
endif()

################################################################################

set(IQS_FILES 
  chi_matrix.cpp
  highperfkernels.cpp
  spec_kernels.cpp
  qureg_apply1qubitgate.cpp
  qureg_apply2qubitgate.cpp
  qureg_applyctrl1qubitgate.cpp
  qureg_applydiag.cpp
  qureg_apply_channel.cpp
  qureg_applyswap.cpp
  qureg_applytoffoli.cpp
  qureg_expectval.cpp
  qureg_fusion.cpp
  qureg_init.cpp
  qureg_measure.cpp
  qureg_noisysimul.cpp
  qureg_permute.cpp
  qureg_utils.cpp
  qureg_version.cpp
  mpi_env.cpp
  mpi_exception.cpp
  utils.cpp
  rng_utils.cpp
  gate_spec.cpp
  CACHE INTERNAL "" FORCE)

if(IqsBuildAsStatic)
  add_library(iqs STATIC ${IQS_FILES})
else()
  add_library(iqs SHARED ${IQS_FILES})
endif()

target_include_directories(iqs PUBLIC ../include)
target_include_directories(iqs INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>)

if (IqsNoise)
    add_dependencies(iqs deps-eigen)
    target_include_directories(iqs PUBLIC ${CMAKE_BINARY_DIR}/deps/eigen-prefix/src/deps-eigen/Eigen)
    # When using ICC, there is a verbose warning from the external library 'eigen'. Disable the warning:
    if(CMAKE_CXX_COMPILER_ID MATCHES Intel)
        target_compile_options(iqs PUBLIC -wd2196)
    endif()
endif()

set_target_properties(iqs
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib"
    POSITION_INDEPENDENT_CODE 1
)

#add_subdirectory(interface)
