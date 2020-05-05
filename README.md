![C++ build with CMake](https://github.com/iqusoft/intel-qs/workflows/C++%20build%20with%20CMake/badge.svg?branch=test%2Fgithub-action)
[![arXiv](https://img.shields.io/static/v1?label=arXiv&message=2001.10554&color=success)](https://arxiv.org/abs/2001.10554)
[![arXiv](https://img.shields.io/static/v1?label=arXiv&message=1601.07195&color=inactive)](https://arxiv.org/abs/1601.07195)

# Intel Quantum Simulator

Intel Quantum Simulator (Intel-QS), also known as qHiPSTER (The Quantum High Performance Software Testing Environment),
is a simulator of quantum circuits optimized to take maximum advantage of multi-core and multi-nodes architectures.
It is based on a complete representation of the qubit state, but avoids the explicit representation of gates and
other quantum operations in terms of matrices.
Intel-QS uses MPI (message-passing-interface) protocols to handle the communication between distributed
resources that are used to store and manipulate the quantum state.

----
## Build instructions

Intel-QS builds as a shared library which, once linked to the application program, allows to take advantage
of the high-performance implementation of circuit simulations.
The library can be built on a variety of different systems, from laptop to HPC server systems.

The directory structure of the repository can be found in
[intel-qs/docs/directory_structure.md](/docs/directory_structure.md).

The library object is: `/builb/lib/libiqs.so`


### Requirements

The following packages are required by the installation:

*  CMake tools version 3.12+
*  MPICH3 library for enabling the distributed communication
*  optional: MKL for distributed random number generation
*  optional: PyBind11 (installed via conda, not pip) required by the Python bunding of Intel-QS

The first step is cloning the repository:
```bash
  git clone https://github.com/iqusoft/intel-qs.git
  cd intel-qs
```

### Use Intel Parallel Studio compilers to build Intel-QS

If you wish to build Intel-QS using the latest Intel compiler technologies, then
you need to configure your environment properly according to that tool's documentation.
Assuming that you have installed Intel Parallel Studio in the standard location on your
system, you should invoke the following scripts through the source command on Linux.
```bash
  source /opt/intel/bin/compilervars.sh -arch intel64 -platform linux
  source /opt/intel/compiler_and_libraries/linux/mpi/intel64/bin/mpivars.sh
```

Now, use CMake to generate the appropriate makefiles to use the Intel Parallel Studio compilers.
The installation follows the out-of-source building and requires the creation of the directory `build`.
This directory is used to collect all the files generated during the installation process.
```bash
  mkdir build
  cd build
  CXX=mpiicpc cmake -DIqsMPI=ON -DIqsUtest=ON ..
  make
```
By default, MKL is required when Intel compilers are used.

To re-build Intel-QS with different settings or options, we recommend to delete all content of the
`build` directory and then restart from the CMake command.


### Use standard GNU tools to build Intel-QS

If you wish to build Intel-QS using only standard GNU compilers type:
```bash
  mkdir build
  cd build
  CXX=g++ cmake -DIqsMPI=OFF ..
  make
```
By default, MKL is not required when GNU compilers are used.
Optionally, MPI can be included by setting the option `-DIqsMPI=ON` instead. You must ensure
that you have at least version 3.1 of MPICH installed for the build to succeed.
https://www.mpich.org


### Enable MPI protocol for distributed memory use

The above installation enables MPI functionalities to deploy Intel-QS on High Performance
Computing and Cloud Computing infrastructures. There is the option of disabling MPI:
simply set the CMake option selection to `-DIqsMPI=OFF`
(or just omit the option selection since MPI is disabled by default in the CMake build).


### Enable Latest Vector Capability

To compile with the latest instruction set supported by your architecture, there is the option `-DIqsNative`. 
Compiled with `-DIqsNative=ON`, the latest vector instructions available on your machine, e.g. AVX2, AVX512, are used.
By default, `-DIqsNative=OFF`.

If the machine you compile and the machine you run have different vector capabilities, turning on `IqsNative=ON` might cause run-time problems.

Underneath, this option uses [`-xhost`](https://software.intel.com/en-us/cpp-compiler-developer-guide-and-reference-xhost-qxhost)
with Intel compilers and [`-march=native`](https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html) with GNU compilers.


### Enable Python binding (only available without MPI)

By default, whenever MPI is disabled, the building process includes the Python binding for
Intel-QS. The binding code uses the Pybind11 library which needs to be installed via 'conda'
(and not simply with pip) to include the relevant information in CMake.
See [this page](https://github.com/pybind/pybind11/issues/1628) for more info on this issue.

To disable the Python wrap, even without MPI, set the CMake option selection to
`-DIqsPython=OFF`.


### Unit test

By default, with MPI either enabled or disabled, the building process includes a suite
of unit tests written in the [googletest framework](https://github.com/google/googletest).
Following the recommended integration, the CMake building process automatically downloads
the up-to-date repository of gtest and installs it in the `build` path.

To disable the unit tests, set the CMake option selection to `-DIqsUtest=OFF`.

To run the unit tests, from `/build` launch the executable `./bin/utest`.


### Recommended build for HPC.

The recommended building process requires
[Intel Math Kernel Library](https://software.intel.com/en-us/mkl)
and the [MPI-ICPC compiler](https://software.intel.com/en-us/node/528770).

When the program is run in hybrid configuration (OpenMP+MPI), we recommend to manage
the OpenMP affinity directly. Affinity settings can be set using the syntax:
`KMP_AFFINITY=compact,1,0,granularity=fine`.
A quick look at the options can be found at
[this page](https://www.nas.nasa.gov/hecc/support/kb/using-intel-openmp-thread-affinity-for-pinning_285.html).


----
## Docker: build image and run/execute container

`Dockerfile` includes the instructions to build the docker image of an Ubuntu machine
with Intel-QS already installed. The image can be 'run' to create a container.
The container can be 'executed' to login into the machine.

```bash
  docker build -t qhipster .
  docker run -d -t qhipster
  docker ps
  docker exec -itd <container_id> /bin/bash
```

If Docker is used on a Windows host machine, the last line should be substituted by:
`winpty docker exec -itd <container_id> //bin/bash`.


----
## Getting started with Intel-QS

The simplest way of familiarize with the Intel Quantum Simulator is by exploring
the tutorials provided in the directory `tutorials/`.
In particular, the code `tutorials/get_started_with_IQS.cpp` provides step-by-step
description of the main commands to:
define a qubit register object, perform quantum gates, measure one or multiple qubits.

If the Python bindings were enabled, the same learning can be performed using the iPython
notebook `tutorials/get_started_with_IQS.ipynb`.


----
## How to contribute

Thanks for your interest in the project! We welcome pull requests from developers
of all skill levels. If you would like to contribute to Intel-QS, please take a
look to our [contributing policy](CONTRIBUTING.md) and also to the 
[code of conduct](CODE_OF_CONDUCT.md). 
For any bug, we use GitHub issues [GitHub issues](https://github.com/iqusoft/intel-qs/issues). Please submit your request there.

----
## How to contact us

If you have a question or want to discuss something, feel free to send an email to
[Justin Hogaboam](justin.w.hogaboam@intel.com),
[Gian Giacomo Guerreschi](gian.giacomo.guerreschi@intel.com), or to
[Fabio Baruffa](fabio.baruffa@intel.com).

----
## How to cite

When using Intel Quantum Simulator for research projects, please cite:

   Gian Giacomo Guerreschi, Justin Hogaboam, Fabio Baruffa, Nicolas P. D. Sawaya
   *Intel Quantum Simulator: A cloud-ready high-performance simulator of quantum circuits*
   [arXiv:2001.10554](https://arxiv.org/abs/2001.10554)

The original implementation is described here: 

   Mikhail Smelyanskiy, Nicolas P. D. Sawaya, Alán Aspuru-Guzik
   *qHiPSTER: The Quantum High Performance Software Testing Environment*
   [arXiv:1601.07195](https://arxiv.org/abs/1601.07195)
