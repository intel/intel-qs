################################################################################
# Build the Python library for the IQS wrapper.

#find_package(pybind11 REQUIRED)
pybind11_add_module(intelqs_py SHARED
                    intelqs_py.cpp
                    ../src/qaoa_features.cpp)

target_link_libraries(intelqs_py PRIVATE iqs)

set_target_properties( intelqs_py
    PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib"
)

################################################################################
