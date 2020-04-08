add_executable(basic_code_for_scaling.exe basic_code_for_scaling.cpp)
target_link_libraries(basic_code_for_scaling.exe intel_qs)

set_target_properties( basic_code_for_scaling.exe
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
)


