add_executable(get_started_with_IQS.exe get_started_with_IQS.cpp)
target_link_libraries(get_started_with_IQS.exe intel_qs)

add_executable(get_started_with_noisy_IQS.exe get_started_with_noisy_IQS.cpp)
target_link_libraries(get_started_with_noisy_IQS.exe intel_qs)

#    add_executable(basic_code_for_scaling.exe benchmarks/basic_code_for_scaling.cpp)
#    target_link_libraries(basic_code_for_scaling.exe intel_qs)

set_target_properties( get_started_with_IQS.exe
                       get_started_with_noisy_IQS.exe
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin"
)
