# Directory structure of Intel-QS

```
.
|-- README.md
|-- CMakeLists.txt
|-- <files with auxiliary information>
|
|-- \docs
|   |-- <documentation files>
|
|-- \benchmarks
|   |-- CMakeLists.txt
|   |-- <benchmarking files and scripts>
|
|-- \interface
|   |-- CMakeLists.txt
|   |-- \include
|       |-- <header files for the interface application>
|   |-- \src
|       |-- CMakeLists.txt
|       |-- <source files for the interface application>
|
|-- \include
|   |-- qureg.hpp
|   |-- <header files for utility functions>
|
|-- \src
|   |-- <implementation of qureg methods>
|   |-- <implementation of utility functions>
|
|-- \pybind11
|   |-- intelqs_py.cpp
|
|-- \cmake
|   |-- FindMKL.cmake
|   |-- gtest.cmake.in
|
|-- \tutorials
|   |-- CMakeLists.txt
|   |-- get_started_with_IQS.cpp
|   |-- get_started_with_IQS.ipynb
|   |-- get_started_with_noisy_IQS.cpp
|
|-- \unit_test
|   |-- suite_of_tests.cpp
|   |-- \include
|   |   |-- <header files with unit tests>
|   |-- \data
|       |-- <files used in testing>
|
|-- \examples
|   |-- <C++ main files of small applications>
|
|-- \notebooks
    |-- <iPython notebooks of small applications>
```
