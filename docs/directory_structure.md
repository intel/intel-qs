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
|
|-- \tutorials
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
    |-- <C++ main files of small applications>
    |-- <iPython notebooks of small applications>
```
