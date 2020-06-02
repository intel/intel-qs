################################################################################
# Construction of IQS as a shared library
################################################################################

set(IQS_FILES 
  highperfkernels.cpp
  qureg_apply1qubitgate.cpp
  qureg_apply2qubitgate.cpp
  qureg_applyctrl1qubitgate.cpp
  qureg_applydiag.cpp
  qureg_applyswap.cpp
  qureg_applytoffoli.cpp
  qureg_expectval.cpp
  qureg_fusion.cpp
  qureg_init.cpp
  qureg_ncu.cpp
  qureg_measure.cpp
  qureg_noisysimul.cpp
  qureg_permute.cpp
  qureg_utils.cpp
  qureg_version.cpp
  mpi_env.cpp
  mpi_exception.cpp
  utils.cpp
  rng_utils.cpp
  CACHE INTERNAL "" FORCE)

add_library(iqs SHARED ${IQS_FILES})

target_include_directories(iqs PUBLIC ${CMAKE_SOURCE_DIR}/include)

target_include_directories(iqs INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>)

set_target_properties(iqs
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib"
    POSITION_INDEPENDENT_CODE 1
)

#add_subdirectory(interface)
