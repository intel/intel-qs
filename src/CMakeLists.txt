################################################################################
# Set some compiler options related to C++.
################################################################################

set(IQS_FILES 
  interface_api_memory.cpp
  interface_api_qasm.cpp
  interface_api_qubitid.cpp
  interface_api_version.cpp
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

#add_library(iqs_shared SHARED ${IQS_FILES})
#target_include_directories(iqs_shared PUBLIC ../include)

add_library(intel_qs STATIC ${IQS_FILES})
target_include_directories(intel_qs PUBLIC ../include)

target_include_directories(intel_qs INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>)

set_target_properties(intel_qs
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib"
    POSITION_INDEPENDENT_CODE 1
)
