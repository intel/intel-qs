add_executable(get_started_with_IQS.exe get_started_with_IQS.cpp)
target_link_libraries(get_started_with_IQS.exe PRIVATE iqs)

add_executable(get_started_with_noisy_IQS.exe get_started_with_noisy_IQS.cpp)
target_link_libraries(get_started_with_noisy_IQS.exe iqs)

set_target_properties( get_started_with_IQS.exe
                       get_started_with_noisy_IQS.exe
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
)
