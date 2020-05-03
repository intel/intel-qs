set(IQS_INTERFACE_FILES
  intelqs_interface_main.cpp
  interface_api_memory.cpp
  interface_api_qasm.cpp
  interface_api_qubitid.cpp
  interface_api_version.cpp
  CACHE INTERNAL "" FORCE)

add_executable(iqs_interface.exe ${IQS_INTERFACE_FILES})
target_link_libraries(iqs_interface.exe iqs)

set_target_properties( iqs_interface.exe
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../bin"
)



