# Change Log

All notable changes to this project will be documented in this file.

----

## [ next ] - [ TBD ]

### Added
- Python bindings with MPI
- new gate specialization for OpenMP performance improvement
- distributed implementation of SWAP-type gates
- class `Permutation` and `PermuteQubits` method to manipulate the representation
  of quantum states and possibly reduce the overhead from MPI communication

### Changed
- CMake build is now modular

### Removed

### Fixed
- when compiling with GCC, the option `-O3` is used for performance improvement
- several minor bugs

----

## [ 2.0.0 ] - [ 2020-12-12 ]

### Added
- pool of state functionality: parallel execution of multiple circuits
- Python bindings based on [PyBind11](https://github.com/pybind/pybind11) (no MPI)
- tutorials and a few new examples
- methods to simulate with noise and decoherence within QubitRegister class
- benchmarks to reproduce weak and strong scaling of 1-qubit gate runtime
- methods specialized for emulation of QAOA circuits

### Changed
- building process based on CMake instead of Make
- main class name changed from `QbitRegister` to `QubitRegister`
- MPI environment extended to support the "pool of state" functionality
- examples updated to new syntax and methods

### Removed

### Fixed

----

## [1.0.0 ] - [ 2017-11-06 ]

Hosted by the (now deprecated) repository <https://github.com/intel/Intel-QS>.

### Added
- implementation of customized MPI environment
- implementation of `QbitRegister` class for simulation of quantum circuits
- methods to simulate 1-qubit gates and controlled-1-qubit gates,
  state preparation and measurement

### Changed

### Removed

### Fixed

