include_directories(../include)

# Build the suite of unit-tests for the project.
# Approach suggested in:
# https://github.com/google/googletest/blob/master/googletest/README.md

# Download and unpack googletest at configure time.
configure_file(../cmake/gtest.cmake.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
RESULT_VARIABLE result
WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker settings on Windows.
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)

################################################################################

# Link against gtest or gtest_main as needed.
add_executable(utest suite_of_tests.cpp
               ../src/qaoa_features.cpp)
target_include_directories(utest PRIVATE include)
target_link_libraries(utest PRIVATE intel_qs gtest_main)

add_executable(mpi_env_test.exe mpi_env_test.cpp
               ../src/mpi_env.cpp
               ../src/mpi_exception.cpp)
target_link_libraries(mpi_env_test.exe intel_qs)

add_executable(mpi_speedup_test.exe mpi_speedup_test.cpp
               ../src/mpi_env.cpp
               ../src/mpi_exception.cpp)
target_link_libraries(mpi_speedup_test.exe intel_qs)

add_executable(openmp_test.exe openmp_test.cpp)

################################################################################

set_target_properties( utest
                       mpi_env_test.exe
                       mpi_speedup_test.exe
                       openmp_test.exe
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin"
)

################################################################################
