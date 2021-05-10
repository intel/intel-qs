################################################################################
# Build the Python library for the IQS wrapper.

execute_process (
    COMMAND bash -c "git submodule add -b stable git@github.com:pybind/pybind11.git extern/pybind11"
)

execute_process (
    COMMAND bash -c "git submodule update --init"
)

#find_package(pybind11 REQUIRED)
add_subdirectory(extern/pybind11)

if (NOT IqsMPI)
    set(NAME_PY_LIB intelqs_py)
else()
    set(NAME_PY_LIB intelqs_py)
endif()

pybind11_add_module(${NAME_PY_LIB} SHARED
                    intelqs_py.cpp
                    ../src/qaoa_features.cpp)

#target_include_directories(${NAME_PY_LIB} PRIVATE ../deps/pybind11/include)
target_link_libraries(${NAME_PY_LIB} PRIVATE iqs)

set_target_properties(${NAME_PY_LIB}
    PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib"
)

################################################################################
