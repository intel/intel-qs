![C++ build with CMake](https://github.com/iqusoft/intel-qs/workflows/C++%20build%20with%20CMake/badge.svg)
![Python build (no MPI)](https://github.com/iqusoft/intel-qs/workflows/Python%20build%20(no%20MPI)/badge.svg)
[![Published Dockerfile](https://img.shields.io/badge/docker%20build-passing-181717?style=flat-square&logo=github&labelColor=black&color=brightgreen)](https://github.com/iqusoft/intel-qs/blob/development/Dockerfile)
[![Quantum Science and Technology](https://img.shields.io/static/v1?label=QST&message=doi:10.1088/2058-9565/ab8505&color=success)](https://iopscience.iop.org/article/10.1088/2058-9565/ab8505)
[![arXiv](https://img.shields.io/static/v1?label=arXiv&message=1601.07195&color=success)](https://arxiv.org/abs/1601.07195)


# Intel Quantum Simulator

Intel Quantum Simulator (Intel-QS), also known as qHiPSTER (The Quantum High Performance Software Testing Environment),
is a simulator of quantum circuits optimized to take maximum advantage of multi-core and multi-nodes architectures.
It is based on a complete representation of the qubit state, but avoids the explicit representation of gates and
other quantum operations in terms of matrices.
Intel-QS uses the MPI (message-passing-interface) protocol to handle communication between the distributed
resources used to store and manipulate quantum states.


## Temporary notice: backward compatibility of April 2021 release

Intel-QS team is aware of the importance of backward compatibility. We do our best to assure it.
In the latest release we adopted good-coding practices and moved a few classes and methods under
the namespace `iqs`. This may cause disruption in older programs. The fix is simple, add `iqs::`
in front of declaration of objects like `QubitRegister`. Other namespaces like `qhipster` have
been susbtituted with namespace `iqs` too.


## Build instructions

Intel-QS builds as a shared library which, once linked to the application program, allows to take advantage
of the high-performance implementation of circuit simulations.
The library can be built on a variety of different systems, from laptop to HPC server systems.

The directory structure of the repository can be found in
[docs/directory_structure.md](/docs/directory_structure.md).

The complete guide to the installation can be found in
[docs/install_guide.md](/docs/install_guide.md).

At the end of the installation, the library object will be: `/builb/lib/libiqs.so`


### Requirements

The following packages are required by the installation:

*  CMake tools version 3.12+
*  MPICH3 library for enabling the distributed communication
*  optional: MKL for distributed random number generation
*  optional: PyBind11 (installed via conda, not pip) required by the Python binding of Intel-QS
*  optional: GoogleTest (automatically installed if needed during the build) required by the unit tests

The first step is cloning the repository:
```bash
  git clone https://github.com/iqusoft/intel-qs.git
  cd intel-qs
```


### Use standard GNU tools to build Intel-QS

Here we describe the basic build using the open-source GNU compiler.
For high-performance computing applications, we suggest adopting the
recommended build detailed in the [installation guide](/docs/install_guide.md).
The installation follows the out-of-source building and requires the creation of the directory `build`.
This directory is used to collect all the files generated during the installation process.

```bash
  mkdir build
  cd build
  CXX=g++ cmake -DIqsMPI=ON -DIqsUtest=ON -DIqsPython=ON -DBuildExample=ON ..
  make
```
The install is customizable and, above, we have chosen to use MPI, compile the
unit tests (based on [GoogleTest framework](https://github.com/google/googletest)),
create a Python library via [PyBind11](https://github.com/pybind/pybind11),
and compile a set of C++ examples.

To re-build Intel-QS with different settings or options, we recommend to delete all content of the
`build` directory and then restart from the CMake command.



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
`winpty docker exec -it <container_id> //bin/bash`.

More detailed instructions can be found in 
[intel-qs/docs/docker_guide.md](/docs/docker_guide.md),
together with instructions to launch a Jupyter notebook from within the container.



## Getting started with Intel-QS

The simplest way of familiarize with the Intel Quantum Simulator is by exploring
the tutorials provided in the directory `tutorials/`.
In particular, the code `tutorials/get_started_with_IQS.cpp` provides step-by-step
description of the main commands to:
define a qubit register object, perform quantum gates, measure one or multiple qubits.

If the Python bindings were enabled, the same learning can be performed using the iPython
notebook `tutorials/get_started_with_IQS.ipynb`.



## How to contribute or contact us

Thanks for your interest in the project! We welcome pull requests from developers
of all skill levels. If you would like to contribute to Intel-QS, please take a
look to our [contributing policy](CONTRIBUTING.md) and also to the 
[code of conduct](CODE_OF_CONDUCT.md). 
For any bug, we use GitHub issues [GitHub issues](https://github.com/iqusoft/intel-qs/issues). Please submit your request there.

If you have a question or want to discuss something, feel free to send an email to
[Justin Hogaboam](justin.w.hogaboam@intel.com),
[Gian Giacomo Guerreschi](gian.giacomo.guerreschi@intel.com), or to
[Fabio Baruffa](fabio.baruffa@intel.com).



## How to cite

When using Intel Quantum Simulator for research projects, please cite:

   Gian Giacomo Guerreschi, Justin Hogaboam, Fabio Baruffa, Nicolas P. D. Sawaya
   *Intel Quantum Simulator: A cloud-ready high-performance simulator of quantum circuits*
   [Quantum Sci. Technol. 5, 034007 (2020)](https://doi.org/10.1088/2058-9565/ab8505)

The original implementation is described here: 

   Mikhail Smelyanskiy, Nicolas P. D. Sawaya, Alán Aspuru-Guzik
   *qHiPSTER: The Quantum High Performance Software Testing Environment*
   [arXiv:1601.07195](https://arxiv.org/abs/1601.07195)
