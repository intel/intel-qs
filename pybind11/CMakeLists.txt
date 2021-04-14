################################################################################
# Build the Python library for the IQS wrapper.

find_package(pybind11 REQUIRED)

if (NOT IqsMPI)
    set(NAME_PY_LIB intelqs_py)
else()
    set(NAME_PY_LIB intelqs_py)
endif()

pybind11_add_module(${NAME_PY_LIB} SHARED
                    intelqs_py.cpp
                    ../src/qaoa_features.cpp)

target_link_libraries(${NAME_PY_LIB} PRIVATE iqs)

set_target_properties(${NAME_PY_LIB}
    PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib"
)

################################################################################
