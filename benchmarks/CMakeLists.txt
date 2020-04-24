add_executable(basic_code_for_scaling.exe basic_code_for_scaling.cpp)
target_link_libraries(basic_code_for_scaling.exe PRIVATE iqs)

set_target_properties( basic_code_for_scaling.exe
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
)


