add_executable(benchgates.exe benchgates.cpp)
target_link_libraries(benchgates.exe PRIVATE intel_qs)

add_executable(quantum_fourier_transform.exe quantum_fourier_transform.cpp)
target_link_libraries(quantum_fourier_transform.exe intel_qs)

add_executable(expect_value_test.exe expect_value_test.cpp)
target_link_libraries(expect_value_test.exe intel_qs)

add_executable(heisenberg_dynamics.exe heisenberg_dynamics.cpp)
target_link_libraries(heisenberg_dynamics.exe intel_qs)

add_executable(noisy_circuit_test.exe noisy_circuit_test.cpp)
target_link_libraries(noisy_circuit_test.exe intel_qs)

add_executable(test_of_custom_gates.exe test_of_custom_gates.cpp)
target_link_libraries(test_of_custom_gates.exe intel_qs)

set_target_properties( benchgates.exe
                       quantum_fourier_transform.exe
                       expect_value_test.exe
                       heisenberg_dynamics.exe
                       noisy_circuit_test.exe
                       test_of_custom_gates.exe
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin"
)


