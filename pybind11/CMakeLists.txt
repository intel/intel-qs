################################################################################
# Build the Python library for the IQS wrapper.

execute_process (
	COMMAND bash -c "git submodule add -b stable git@github.com:pybind/pybind11.git ${PROJECT_SOURCE_DIR}/extern/pybind11"
    COMMAND bash -c "git submodule update --init"
)

add_subdirectory(${PROJECT_SOURCE_DIR}/extern/pybind11 ${PROJECT_BINARY_DIR}/extern/pybind11)

if (NOT IqsMPI)
    set(NAME_PY_LIB intelqs_py)
else()
    set(NAME_PY_LIB intelqs_py)
endif()

pybind11_add_module(${NAME_PY_LIB} SHARED
                    intelqs_py.cpp
                    ../src/chi_matrix.cpp
                    ../src/qaoa_features.cpp)

target_link_libraries(${NAME_PY_LIB} PRIVATE iqs)

set_target_properties(${NAME_PY_LIB}
    PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib"
)

################################################################################
