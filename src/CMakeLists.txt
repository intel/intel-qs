cmake_minimum_required(VERSION 3.12)
project(iqs)

################################################################################
# Set some compiler options related to C++.
################################################################################

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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
  CACHE INTERNAL "" FORCE)

add_library(iqs SHARED ${IQS_FILES})
add_library(iqstatic STATIC ${IQS_FILES})
set_target_properties(iqs PROPERTIES POSITION_INDEPENDENT_CODE 1)

target_include_directories(iqs INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)





